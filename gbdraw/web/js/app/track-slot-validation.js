import { parseDepthTrackIndexIdentity } from './depth-track-state.js';

const normalizedString = (value) => String(value ?? '').trim().toLowerCase();

const normalizedOptionalString = (value) => {
  const text = normalizedString(value);
  return text || null;
};

const LINEAR_GENERIC_LAYOUT_KEYS = new Set([
  'id', 'renderer', 'type', 'side', 'h', 'height', 'spacing',
  'z', 'z_index', 'zindex', 'enabled', 'show', 'visible'
]);
const CIRCULAR_GENERIC_LAYOUT_KEYS = new Set([
  'side', 'r', 'radius', 'w', 'width', 'spacing', 'inner_gap_px',
  'outer_gap_px', 'z', 'z_index', 'zindex', 'strict', 'compress',
  'reserve', 'enabled', 'show', 'visible', 'ri', 'inner',
  'inner_radius', 'ro', 'outer', 'outer_radius', 'gap', 'gap_after',
  'gapafter', 'innerradius', 'outerradius'
]);

const parseSlotZIndex = (value, fieldName) => {
  if (value === null || value === undefined || value === '') return 0;
  if (typeof value === 'boolean') {
    throw new Error(`${fieldName} must be an integer.`);
  }
  const numeric = (
    typeof value === 'number' ||
    (typeof value === 'string' && /^[+-]?\d+$/.test(value.trim()))
  ) ? Number(value) : Number.NaN;
  if (!Number.isSafeInteger(numeric)) {
    throw new Error(`${fieldName} must be an integer.`);
  }
  return numeric;
};

const parseOptionalLinearPx = (value, fieldName, { allowZero }) => {
  if (value === null || value === undefined || value === '') return null;
  if (typeof value === 'boolean' || (typeof value !== 'number' && typeof value !== 'string')) {
    throw new Error(`${fieldName} must be a number of pixels.`);
  }
  const text = String(value).trim();
  if (!text || /%$/i.test(text)) {
    throw new Error(`${fieldName} only accepts px or unitless px values.`);
  }
  const numericText = /px$/i.test(text) ? text.slice(0, -2).trim() : text;
  const numeric = Number(numericText);
  if (!Number.isFinite(numeric)) {
    throw new Error(`${fieldName} must be a number of pixels.`);
  }
  if (numeric < 0 || (!allowZero && numeric === 0)) {
    throw new Error(`${fieldName} must be ${allowZero ? 'nonnegative' : 'positive'}.`);
  }
  return numeric;
};

const parseOptionalCircularScalar = (value, fieldName) => {
  if (value === null || value === undefined || value === '') return null;
  if (typeof value === 'boolean' || (typeof value !== 'number' && typeof value !== 'string')) {
    throw new Error(`${fieldName} must be a finite scalar.`);
  }
  const text = String(value).trim();
  let numericText = text;
  if (/px$/i.test(text)) numericText = text.slice(0, -2).trim();
  else if (/%$/.test(text)) numericText = text.slice(0, -1).trim();
  const numeric = Number(numericText);
  if (!text || !Number.isFinite(numeric)) {
    throw new Error(`${fieldName} must be a finite scalar.`);
  }
  return numeric;
};

const parseOptionalCircularGap = (value, fieldName) => {
  if (value === null || value === undefined || value === '') return null;
  if (typeof value === 'boolean' || (typeof value !== 'number' && typeof value !== 'string')) {
    throw new Error(`${fieldName} must be a nonnegative numeric pixel value.`);
  }
  const text = String(value).trim();
  if (!text || /(?:px|%)$/i.test(text)) {
    throw new Error(`${fieldName} must be a numeric pixel value without a unit.`);
  }
  const numeric = Number(text);
  if (!Number.isFinite(numeric) || numeric < 0) {
    throw new Error(`${fieldName} must be a nonnegative numeric pixel value.`);
  }
  return numeric;
};

const validateSlotGeometry = (slot, id, layoutKind) => {
  if (layoutKind === 'linear') {
    parseOptionalLinearPx(slot.height, `Linear track slot '${id}' height`, { allowZero: false });
    parseOptionalLinearPx(slot.spacing, `Linear track slot '${id}' spacing`, { allowZero: true });
    return;
  }
  if (layoutKind !== 'circular') {
    throw new Error(`Unsupported track-slot layout kind: ${layoutKind}.`);
  }
  const spacing = parseOptionalCircularScalar(
    slot.spacing,
    `Circular track slot '${id}' spacing`
  );
  parseOptionalCircularScalar(slot.radius, `Circular track slot '${id}' radius`);
  parseOptionalCircularScalar(slot.width, `Circular track slot '${id}' width`);
  const innerGap = parseOptionalCircularGap(
    slot.inner_gap_px,
    `Circular track slot '${id}' inner_gap_px`
  );
  const outerGap = parseOptionalCircularGap(
    slot.outer_gap_px,
    `Circular track slot '${id}' outer_gap_px`
  );
  if (spacing !== null && (innerGap !== null || outerGap !== null)) {
    throw new Error(
      `Circular track slot '${id}' cannot combine spacing with inner_gap_px or outer_gap_px.`
    );
  }
};

const validateGenericParams = (params, id, layoutKind) => {
  const genericKeys = layoutKind === 'linear'
    ? LINEAR_GENERIC_LAYOUT_KEYS
    : CIRCULAR_GENERIC_LAYOUT_KEYS;
  const invalidKey = Object.keys(params).find((key) => genericKeys.has(normalizedString(key)));
  if (invalidKey) {
    throw new Error(
      `${layoutKind === 'linear' ? 'Linear' : 'Circular'} track slot '${id}' stores generic layout field '${invalidKey}' in params.`
    );
  }
};

/**
 * Validate structural and binding invariants shared by canonical slots and
 * persisted UI slots.
 *
 * The returned objects are detached normalized views used only for validation;
 * callers remain responsible for committing their original state after every
 * preflight check succeeds.
 */
export const validateTrackSlotBindingInvariants = (
  slots,
  {
    modeLabel,
    layoutKind,
    supportedRenderers,
    supportedSides,
    anchorlessRenderers,
    depthTrackCount = null
  }
) => {
  if (!Array.isArray(slots)) {
    throw new Error(`${modeLabel} Custom Track Slots must be an array.`);
  }

  const rendererSet = new Set(supportedRenderers || []);
  const sideSet = new Set(supportedSides || []);
  const anchorlessSet = new Set(anchorlessRenderers || []);
  const validateDepthRange = depthTrackCount !== null && depthTrackCount !== undefined;
  if (
    validateDepthRange &&
    (!Number.isSafeInteger(depthTrackCount) || depthTrackCount < 0)
  ) {
    throw new Error(`${modeLabel} Depth track count must be a non-negative integer.`);
  }

  const ids = new Set();
  let enabledFeatureCount = 0;
  const normalized = slots.map((slot, index) => {
    if (!slot || typeof slot !== 'object' || Array.isArray(slot)) {
      throw new Error(`${modeLabel} track slot #${index + 1} must be an object.`);
    }
    if (typeof slot.id !== 'string' || !slot.id.trim()) {
      throw new Error(`${modeLabel} track slot #${index + 1} requires a non-empty id.`);
    }
    const id = slot.id.trim();
    if (ids.has(id)) {
      throw new Error(`Duplicate ${modeLabel} track slot id: ${id}.`);
    }
    ids.add(id);

    if (typeof slot.renderer !== 'string') {
      throw new Error(`${modeLabel} track slot '${id}' renderer must be a string.`);
    }
    const renderer = normalizedString(slot.renderer);
    if (!rendererSet.has(renderer)) {
      throw new Error(`${modeLabel} track slot '${id}' has unsupported renderer=${JSON.stringify(slot.renderer)}.`);
    }
    if (
      Object.prototype.hasOwnProperty.call(slot, 'enabled') &&
      typeof slot.enabled !== 'boolean'
    ) {
      throw new Error(`${modeLabel} track slot '${id}' enabled must be boolean.`);
    }
    const enabled = slot.enabled !== false;

    validateSlotGeometry(slot, id, layoutKind);

    let side = null;
    if (slot.side !== null && slot.side !== undefined && slot.side !== '') {
      if (typeof slot.side !== 'string') {
        throw new Error(`${modeLabel} track slot '${id}' side must be a string or null.`);
      }
      side = normalizedString(slot.side);
      if (!sideSet.has(side)) {
        throw new Error(`${modeLabel} track slot '${id}' has unsupported side=${JSON.stringify(slot.side)}.`);
      }
    }

    const hasParams = Object.prototype.hasOwnProperty.call(slot, 'params');
    const paramsValue = hasParams ? slot.params : {};
    if (!paramsValue || typeof paramsValue !== 'object' || Array.isArray(paramsValue)) {
      throw new Error(`${modeLabel} track slot '${id}' params must be an object.`);
    }
    const params = { ...paramsValue };
    validateGenericParams(params, id, layoutKind);
    const z = parseSlotZIndex(slot.z, `${modeLabel} track slot '${id}' z`);

    if (enabled && layoutKind === 'linear') {
      side = side ?? (renderer === 'annotations' ? 'above' : 'below');
      if (side === 'overlay' && !['features', 'annotations'].includes(renderer)) {
        throw new Error(
          `Linear track slot '${id}' uses side=overlay, which is only supported for features and annotations.`
        );
      }
      if (renderer === 'features') {
        enabledFeatureCount += 1;
        if (enabledFeatureCount > 1) {
          throw new Error('Linear Custom Track Slots support only one enabled features slot.');
        }
      }
    } else if (enabled && layoutKind === 'circular') {
      side = side ?? (renderer === 'annotations' ? 'outside' : 'inside');
      if (renderer === 'sequence_conservation' && side === 'overlay') {
        throw new Error(`Circular sequence_conservation slot '${id}' cannot use side=overlay.`);
      }
      if (renderer === 'features') {
        const lane = normalizedOptionalString(params.lane_direction ?? params.lanes);
        if (lane !== null && !['inside', 'outside', 'split'].includes(lane)) {
          throw new Error(
            `Circular features slot '${id}' has unsupported lane_direction=${JSON.stringify(lane)}.`
          );
        }
        if (lane !== null) {
          const expectedSide = lane === 'split' ? 'overlay' : lane;
          if (side !== expectedSide) {
            throw new Error(
              `Circular features slot '${id}' has conflicting side=${JSON.stringify(side)} and lane_direction=${JSON.stringify(lane)}.`
            );
          }
        }
      }
      if (renderer === 'ticks') {
        const obsoleteTickKey = ['axis', 'label_side', 'tick_side']
          .find((key) => Object.prototype.hasOwnProperty.call(params, key));
        if (obsoleteTickKey) {
          throw new Error(
            `Circular ticks slot '${id}' uses obsolete parameter '${obsoleteTickKey}'.`
          );
        }
      }
    }

    if (enabled && renderer === 'annotations') {
      const anchorId = String(params.anchor_slot || '').trim();
      if (side === 'overlay' && !anchorId) {
        throw new Error(`${modeLabel} annotation slot '${id}' with side=overlay requires anchor_slot.`);
      }
      if (side !== 'overlay' && anchorId) {
        throw new Error(`${modeLabel} annotation slot '${id}' uses anchor_slot without side=overlay.`);
      }
      const layer = normalizedOptionalString(params.layer) ?? 'foreground';
      if (!['foreground', 'underlay'].includes(layer)) {
        throw new Error(
          `${modeLabel} annotation slot '${id}' layer must be foreground or underlay.`
        );
      }
      params.layer = layer;
    }

    if (enabled && renderer === 'depth') {
      const trackIndex = parseDepthTrackIndexIdentity(
        Object.prototype.hasOwnProperty.call(params, 'track_index')
          ? params.track_index
          : 0,
        `${modeLabel} Depth slot '${id}' track_index`
      );
      params.track_index = trackIndex;
      if (validateDepthRange && trackIndex >= depthTrackCount) {
        const available = depthTrackCount === 0
          ? 'no global Depth series are available'
          : `available range is 0..${depthTrackCount - 1}`;
        throw new Error(
          `${modeLabel} Depth slot '${id}' references logical track index ${trackIndex}; ${available}.`
        );
      }
    }

    return { id, renderer, enabled, side, z, params };
  });

  const enabledById = new Map(
    normalized.filter((slot) => slot.enabled).map((slot) => [slot.id, slot])
  );
  normalized.forEach((slot) => {
    if (!slot.enabled || slot.renderer !== 'annotations' || slot.side !== 'overlay') return;
    const anchorId = String(slot.params.anchor_slot || '').trim();
    const anchor = enabledById.get(anchorId);
    if (!anchor) {
      throw new Error(
        `${modeLabel} annotation slot '${slot.id}' references unknown anchor_slot=${JSON.stringify(anchorId)}.`
      );
    }
    if (anchor.renderer === 'annotations' || anchorlessSet.has(anchor.renderer)) {
      throw new Error(
        `${modeLabel} annotation slot '${slot.id}' anchor '${anchor.id}' has no eligible drawable band.`
      );
    }
    const layer = normalizedString(slot.params.layer || 'foreground');
    if (layer === 'underlay' && slot.z >= anchor.z) {
      throw new Error(
        `${modeLabel} annotation underlay slot '${slot.id}' must have z less than anchor '${anchor.id}'.`
      );
    }
    if (layer === 'foreground' && slot.z <= anchor.z) {
      throw new Error(
        `${modeLabel} annotation foreground slot '${slot.id}' must have z greater than anchor '${anchor.id}'.`
      );
    }
  });

  return normalized;
};
