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
const OBSOLETE_CIRCULAR_TRACK_SLOT_KEYS = new Set([
  'spacing', 'strict', 'compress', 'reserve'
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
  if (value && typeof value === 'object' && !Array.isArray(value)) {
    const numeric = Number(value.value);
    const unit = normalizedString(value.unit);
    if (!Number.isFinite(numeric) || numeric <= 0 || !['px', 'factor'].includes(unit)) {
      throw new Error(`${fieldName} must be a positive finite px or factor scalar.`);
    }
    return numeric;
  }
  if (typeof value === 'boolean' || (typeof value !== 'number' && typeof value !== 'string')) {
    throw new Error(`${fieldName} must be a positive finite px or factor scalar.`);
  }
  const text = String(value).trim();
  let numericText = text;
  if (/px$/i.test(text)) numericText = text.slice(0, -2).trim();
  else if (/%$/.test(text)) numericText = text.slice(0, -1).trim();
  const numeric = Number(numericText);
  if (!text || !Number.isFinite(numeric) || numeric <= 0) {
    throw new Error(`${fieldName} must be a positive finite px or factor scalar.`);
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
  const obsoleteKey = Object.keys(slot).find((key) => (
    OBSOLETE_CIRCULAR_TRACK_SLOT_KEYS.has(normalizedString(key))
  ));
  if (obsoleteKey) {
    throw new Error(
      `Circular track slot '${id}' uses obsolete field '${obsoleteKey}'. ` +
      'Use inner_gap_px and outer_gap_px for physical gaps.'
    );
  }
  parseOptionalCircularScalar(slot.radius, `Circular track slot '${id}' radius`);
  parseOptionalCircularScalar(slot.width, `Circular track slot '${id}' width`);
  parseOptionalCircularGap(
    slot.inner_gap_px,
    `Circular track slot '${id}' inner_gap_px`
  );
  parseOptionalCircularGap(
    slot.outer_gap_px,
    `Circular track slot '${id}' outer_gap_px`
  );
};

const validateGenericParams = (params, id, layoutKind) => {
  const genericKeys = layoutKind === 'linear'
    ? LINEAR_GENERIC_LAYOUT_KEYS
    : CIRCULAR_GENERIC_LAYOUT_KEYS;
  const invalidKey = Object.keys(params).find((key) => genericKeys.has(normalizedString(key)));
  if (invalidKey) {
    if (
      layoutKind === 'circular' &&
      OBSOLETE_CIRCULAR_TRACK_SLOT_KEYS.has(normalizedString(invalidKey))
    ) {
      throw new Error(
        `Circular track slot '${id}' uses obsolete field 'params.${invalidKey}'. ` +
        'Use inner_gap_px and outer_gap_px for physical gaps.'
      );
    }
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

const CUSTOM_TRACK_RENDERERS = Object.freeze({
  circular: new Set([
    'features',
    'ticks',
    'dinucleotide_content',
    'dinucleotide_skew',
    'depth',
    'sequence_conservation',
    'annotations',
    'spacer'
  ]),
  linear: new Set([
    'features',
    'dinucleotide_content',
    'dinucleotide_skew',
    'depth',
    'annotations',
    'spacer'
  ])
});

const CUSTOM_TRACK_SIDES = Object.freeze({
  circular: new Set(['inside', 'outside', 'overlay']),
  linear: new Set(['above', 'below', 'overlay'])
});

const ANNOTATION_MARKS = new Set(['line', 'bracket', 'band', 'highlight']);
const INELIGIBLE_ANCHOR_RENDERERS = new Set(['annotations', 'ticks', 'spacer']);
const CONSERVATION_MANAGER = 'circular_conservation';
const PARAM_OWNER = new Map([
  ['lane_direction', 'features'],
  ['lanes', 'features'],
  ['tick_label_layout', 'ticks'],
  ['preset', 'ticks'],
  ['axis', 'ticks'],
  ['label_side', 'ticks'],
  ['tick_side', 'ticks'],
  ['nt', 'dinucleotide'],
  ['dinucleotide', 'dinucleotide'],
  ['positive_color', 'dinucleotide_skew'],
  ['negative_color', 'dinucleotide_skew'],
  ['high_color', 'dinucleotide_skew'],
  ['low_color', 'dinucleotide_skew'],
  ['track_index', 'depth-or-conservation'],
  ['source_index', 'sequence_conservation'],
  ['series_key', 'sequence_conservation'],
  ['managed', 'sequence_conservation'],
  ['set_id', 'annotations'],
  ['marks', 'annotations'],
  ['lane_gap_px', 'annotations'],
  ['padding_px', 'annotations'],
  ['overflow', 'annotations'],
  ['show_labels', 'annotations'],
  ['anchor_slot', 'annotations'],
  ['layer', 'annotations'],
  ['cover_anchor', 'annotations'],
  ['style_override', 'annotations']
]);

const cloneValidationSlot = (slot) => {
  if (!slot || typeof slot !== 'object' || Array.isArray(slot)) return slot;
  return {
    ...slot,
    params: (
      slot.params && typeof slot.params === 'object' && !Array.isArray(slot.params)
        ? { ...slot.params }
        : slot.params
    )
  };
};

const issueForRow = ({
  code,
  field,
  message,
  rowIndex,
  slotId = null
}) => ({
  code,
  field,
  message,
  rowIndex,
  slotId,
  severity: 'error'
});

const globalIssue = ({ code, field, message }) => ({
  code,
  field,
  message,
  severity: 'error'
});

const addRowIssue = (rowIssues, index, issue) => {
  const current = rowIssues.get(index) || [];
  current.push(issueForRow({ ...issue, rowIndex: index }));
  rowIssues.set(index, current);
};

const customTrackModeLabel = (mode) => (
  mode === 'circular' ? 'Circular' : 'Linear'
);

const defaultSideForSlot = (slot, mode, trackType) => {
  const renderer = normalizedString(slot?.renderer);
  if (mode === 'linear') {
    if (renderer === 'annotations') return 'above';
    if (renderer !== 'features') return 'below';
    const layout = normalizedString(trackType);
    if (layout === 'above') return 'above';
    if (layout === 'below') return 'below';
    return 'overlay';
  }
  if (renderer === 'annotations') return 'outside';
  if (renderer !== 'features') return 'inside';
  const lane = normalizedOptionalString(
    slot?.params?.lane_direction ?? slot?.params?.lanes
  );
  if (lane === 'split') return 'overlay';
  if (lane === 'inside' || lane === 'outside') return lane;
  const preset = normalizedString(trackType);
  if (preset === 'middle') return 'overlay';
  if (preset === 'spreadout') return 'outside';
  return 'inside';
};

const inferredCustomTrackAxisIndex = (slots, mode, trackType) => {
  const onAxisIndex = slots.findIndex((slot) => (
    slot?.enabled !== false &&
    normalizedString(slot?.side || defaultSideForSlot(slot, mode, trackType)) === 'overlay' &&
    normalizedString(slot?.renderer) !== 'annotations'
  ));
  if (onAxisIndex >= 0) return onAxisIndex;

  const featureIndex = slots.findIndex((slot) => (
    slot?.enabled !== false && normalizedString(slot?.renderer) === 'features'
  ));
  if (featureIndex >= 0) {
    const side = normalizedString(
      slots[featureIndex]?.side || defaultSideForSlot(slots[featureIndex], mode, trackType)
    );
    return side === (mode === 'linear' ? 'above' : 'outside')
      ? featureIndex + 1
      : featureIndex;
  }

  const outerSide = mode === 'linear' ? 'above' : 'outside';
  const firstInner = slots.findIndex((slot) => (
    normalizedString(slot?.side || defaultSideForSlot(slot, mode, trackType)) !== outerSide
  ));
  return firstInner < 0 ? slots.length : firstInner;
};

const featureUnderlaysVisible = (value) => {
  if (Array.isArray(value)) return value.length > 0;
  if (typeof value === 'number') return value > 0;
  return Boolean(value);
};

const nonnegativeFiniteNumber = (value) => {
  if (value === null || value === undefined || value === '') return null;
  if (typeof value === 'boolean') return Number.NaN;
  const numeric = Number(value);
  return Number.isFinite(numeric) && numeric >= 0 ? numeric : Number.NaN;
};

const paramsAllowedForRenderer = (renderer, key) => {
  const owner = PARAM_OWNER.get(key);
  if (!owner) return true;
  if (owner === renderer) return true;
  if (owner === 'dinucleotide') {
    return renderer === 'dinucleotide_content' || renderer === 'dinucleotide_skew';
  }
  if (owner === 'depth-or-conservation') {
    return renderer === 'depth' || renderer === 'sequence_conservation';
  }
  return false;
};

const normalizedInventoryIds = (values) => new Set(
  (Array.isArray(values) ? values : [])
    .map((value) => String(value?.id ?? value ?? '').trim())
    .filter(Boolean)
);

const conservationInventory = (series) => {
  const keys = new Set();
  const indexes = new Set();
  const indexesByKey = new Map();
  (Array.isArray(series) ? series : []).forEach((entry, index) => {
    const sourceKey = String(
      entry?.sourceKey ?? entry?.series_key ?? entry?.seriesKey ?? ''
    ).trim();
    if (sourceKey) {
      keys.add(sourceKey);
      if (!indexesByKey.has(sourceKey)) indexesByKey.set(sourceKey, new Set());
    }
    [
      entry?.sourceIndex,
      entry?.source_index,
      entry?.orderIndex,
      entry?.order_index,
      index
    ].forEach((value) => {
      const numeric = Number(value);
      if (!Number.isSafeInteger(numeric) || numeric < 0) return;
      indexes.add(numeric);
      if (sourceKey) indexesByKey.get(sourceKey).add(numeric);
    });
  });
  return { keys, indexes, indexesByKey };
};

const collectCustomTrackIssues = ({
  mode,
  slots,
  axisIndex,
  trackType,
  depthTrackCount,
  annotationSetIds,
  visibleFeatureUnderlays,
  conservationSeries
}) => {
  const rowIssues = new Map();
  const globalIssues = [];
  const draftSlots = Array.isArray(slots)
    ? slots.map(cloneValidationSlot)
    : [];
  const modeLabel = customTrackModeLabel(mode);

  if (!Array.isArray(slots)) {
    globalIssues.push(globalIssue({
      code: 'slots_not_array',
      field: 'slots',
      message: `${modeLabel} Custom Track Slots must be an array.`
    }));
    return {
      draftSlots,
      enabledSlots: [],
      emittedAxisIndex: null,
      rowIssues,
      globalIssues
    };
  }

  const supportedRenderers = CUSTOM_TRACK_RENDERERS[mode];
  const supportedSides = CUSTOM_TRACK_SIDES[mode];
  if (!supportedRenderers || !supportedSides) {
    globalIssues.push(globalIssue({
      code: 'mode_unsupported',
      field: 'mode',
      message: `Unsupported Custom Track mode: ${String(mode)}.`
    }));
    return {
      draftSlots,
      enabledSlots: [],
      emittedAxisIndex: null,
      rowIssues,
      globalIssues
    };
  }

  const annotationInventory = annotationSetIds === null || annotationSetIds === undefined
    ? null
    : normalizedInventoryIds(annotationSetIds);
  const conservation = conservationInventory(conservationSeries);
  const depthCountProvided = depthTrackCount !== null && depthTrackCount !== undefined;
  const validDepthCount = (
    !depthCountProvided ||
    (Number.isSafeInteger(depthTrackCount) && depthTrackCount >= 0)
  );
  if (!validDepthCount) {
    globalIssues.push(globalIssue({
      code: 'depth_count_invalid',
      field: 'depthTrackCount',
      message: `${modeLabel} Depth track count must be a non-negative integer.`
    }));
  }

  const enabledIndexes = [];
  draftSlots.forEach((slot, index) => {
    if (!slot || typeof slot !== 'object' || Array.isArray(slot)) {
      addRowIssue(rowIssues, index, {
        code: 'slot_not_object',
        field: 'slot',
        message: `${modeLabel} track row ${index + 1} must be an object.`
      });
      return;
    }
    if (
      Object.prototype.hasOwnProperty.call(slot, 'enabled') &&
      typeof slot.enabled !== 'boolean'
    ) {
      addRowIssue(rowIssues, index, {
        code: 'enabled_not_boolean',
        field: 'enabled',
        slotId: String(slot.id || '').trim() || null,
        message: `${modeLabel} track row ${index + 1} Enabled must be true or false.`
      });
    }
    if (slot.enabled !== false) enabledIndexes.push(index);
  });
  if (mode === 'linear' && enabledIndexes.length === 0) {
    globalIssues.push(globalIssue({
      code: 'linear_slots_empty',
      field: 'slots',
      message: 'Linear Custom Track Slots require at least one enabled row.'
    }));
  }

  let resolvedAxisIndex = null;
  if (axisIndex === null || axisIndex === undefined || axisIndex === '') {
    resolvedAxisIndex = inferredCustomTrackAxisIndex(draftSlots, mode, trackType);
  } else if (
    typeof axisIndex !== 'boolean' &&
    Number.isSafeInteger(Number(axisIndex)) &&
    Number(axisIndex) >= 0 &&
    Number(axisIndex) <= draftSlots.length
  ) {
    resolvedAxisIndex = Number(axisIndex);
  } else {
    globalIssues.push(globalIssue({
      code: 'axis_out_of_range',
      field: 'axisIndex',
      message: `${modeLabel} Axis boundary must be an integer from 0 to ${draftSlots.length}.`
    }));
  }
  const emittedAxisIndex = resolvedAxisIndex === null
    ? null
    : draftSlots
        .slice(0, resolvedAxisIndex)
        .filter((slot) => slot && typeof slot === 'object' && slot.enabled !== false)
        .length;

  const idsToIndexes = new Map();
  enabledIndexes.forEach((index) => {
    const slot = draftSlots[index];
    if (!slot || typeof slot !== 'object' || Array.isArray(slot)) return;
    const id = typeof slot.id === 'string' ? slot.id.trim() : '';
    if (!id) {
      addRowIssue(rowIssues, index, {
        code: 'id_required',
        field: 'id',
        message: `${modeLabel} track row ${index + 1} requires a non-empty ID.`
      });
      return;
    }
    const indexes = idsToIndexes.get(id) || [];
    indexes.push(index);
    idsToIndexes.set(id, indexes);
  });
  idsToIndexes.forEach((indexes, id) => {
    if (indexes.length < 2) return;
    indexes.forEach((index) => addRowIssue(rowIssues, index, {
      code: 'id_duplicate',
      field: 'id',
      slotId: id,
      message: `${modeLabel} track ID '${id}' is used by more than one enabled row.`
    }));
  });

  const normalizedEnabled = [];
  const rowIndexByNormalizedSlot = new WeakMap();
  enabledIndexes.forEach((index) => {
    const slot = draftSlots[index];
    if (!slot || typeof slot !== 'object' || Array.isArray(slot)) return;
    const id = typeof slot.id === 'string' ? slot.id.trim() : '';
    const slotLabel = id ? `'${id}'` : `row ${index + 1}`;
    const renderer = typeof slot.renderer === 'string'
      ? normalizedString(slot.renderer)
      : '';
    if (!renderer || !supportedRenderers.has(renderer)) {
      addRowIssue(rowIssues, index, {
        code: 'renderer_unsupported',
        field: 'renderer',
        slotId: id || null,
        message: `${modeLabel} track ${slotLabel} has an unsupported renderer.`
      });
    }

    let explicitSide = null;
    if (slot.side !== null && slot.side !== undefined && slot.side !== '') {
      explicitSide = typeof slot.side === 'string' ? normalizedString(slot.side) : null;
      if (!explicitSide || !supportedSides.has(explicitSide)) {
        addRowIssue(rowIssues, index, {
          code: 'side_unsupported',
          field: 'side',
          slotId: id || null,
          message: `${modeLabel} track ${slotLabel} has an unsupported side.`
        });
        explicitSide = null;
      }
    }

    const params = (
      slot.params && typeof slot.params === 'object' && !Array.isArray(slot.params)
        ? { ...slot.params }
        : null
    );
    if (!params) {
      addRowIssue(rowIssues, index, {
        code: 'params_not_object',
        field: 'params',
        slotId: id || null,
        message: `${modeLabel} track ${slotLabel} Params must be an object.`
      });
    }

    try {
      validateSlotGeometry(slot, id || `#${index + 1}`, mode);
    } catch (error) {
      addRowIssue(rowIssues, index, {
        code: 'geometry_invalid',
        field: 'geometry',
        slotId: id || null,
        message: error instanceof Error ? error.message : String(error)
      });
    }
    try {
      parseSlotZIndex(slot.z, `${modeLabel} track ${slotLabel} z`);
    } catch (error) {
      addRowIssue(rowIssues, index, {
        code: 'z_invalid',
        field: 'z',
        slotId: id || null,
        message: error instanceof Error ? error.message : String(error)
      });
    }
    if (params) {
      try {
        validateGenericParams(params, id || `#${index + 1}`, mode);
      } catch (error) {
        addRowIssue(rowIssues, index, {
          code: 'generic_param',
          field: 'params',
          slotId: id || null,
          message: error instanceof Error ? error.message : String(error)
        });
      }
      const mismatchedParam = Object.keys(params).find(
        (key) => !paramsAllowedForRenderer(renderer, normalizedString(key))
      );
      if (mismatchedParam) {
        addRowIssue(rowIssues, index, {
          code: 'renderer_param_mismatch',
          field: `params.${mismatchedParam}`,
          slotId: id || null,
          message: `${modeLabel} track ${slotLabel} retains '${mismatchedParam}' from another renderer.`
        });
      }
    }

    const enabledSlotIndex = normalizedEnabled.length;
    const axisDerivedSide = emittedAxisIndex === null
      ? defaultSideForSlot(slot, mode, trackType)
      : (
          mode === 'linear'
            ? (enabledSlotIndex < emittedAxisIndex ? 'above' : 'below')
            : (enabledSlotIndex < emittedAxisIndex ? 'outside' : 'inside')
        );
    let side = explicitSide;
    if (mode === 'circular' && renderer === 'features' && params) {
      const lane = normalizedOptionalString(params.lane_direction ?? params.lanes);
      if (lane !== null && !['inside', 'outside', 'split'].includes(lane)) {
        addRowIssue(rowIssues, index, {
          code: 'feature_lane',
          field: 'params.lane_direction',
          slotId: id || null,
          message: `Circular Features track ${slotLabel} has an invalid lane direction.`
        });
      } else if (lane !== null) {
        const laneSide = lane === 'split' ? 'overlay' : lane;
        if (explicitSide !== null && explicitSide !== laneSide) {
          addRowIssue(rowIssues, index, {
            code: 'feature_lane_side_conflict',
            field: 'side',
            slotId: id || null,
            message: `Circular Features track ${slotLabel} side conflicts with its lane direction.`
          });
        }
        side = laneSide;
      }
    }
    side = side || axisDerivedSide;

    if (mode === 'linear' && side === 'overlay' && !['features', 'annotations'].includes(renderer)) {
      addRowIssue(rowIssues, index, {
        code: 'overlay_renderer_unsupported',
        field: 'side',
        slotId: id || null,
        message: `Linear track ${slotLabel} can use On Axis only for Features or Annotations.`
      });
    }
    if (mode === 'circular' && renderer === 'sequence_conservation' && side === 'overlay') {
      addRowIssue(rowIssues, index, {
        code: 'conservation_overlay',
        field: 'side',
        slotId: id || null,
        message: `Circular comparison track ${slotLabel} cannot use Overlay.`
      });
    }

    if (mode === 'circular' && renderer === 'ticks' && params) {
      const obsoleteTickKey = ['axis', 'label_side', 'tick_side']
        .find((key) => Object.prototype.hasOwnProperty.call(params, key));
      if (obsoleteTickKey) {
        addRowIssue(rowIssues, index, {
          code: 'ticks_obsolete_param',
          field: `params.${obsoleteTickKey}`,
          slotId: id || null,
          message: `Circular Ticks track ${slotLabel} uses obsolete '${obsoleteTickKey}'.`
        });
      }
    }

    if (renderer === 'depth' && params) {
      let trackIndex = null;
      try {
        trackIndex = parseDepthTrackIndexIdentity(
          Object.prototype.hasOwnProperty.call(params, 'track_index')
            ? params.track_index
            : 0,
          `${modeLabel} Depth track ${slotLabel} track_index`
        );
      } catch (error) {
        addRowIssue(rowIssues, index, {
          code: 'depth_track_index',
          field: 'params.track_index',
          slotId: id || null,
          message: error instanceof Error ? error.message : String(error)
        });
      }
      if (trackIndex !== null && validDepthCount && depthCountProvided) {
        if (depthTrackCount === 0) {
          addRowIssue(rowIssues, index, {
            code: 'depth_source_missing',
            field: 'params.track_index',
            slotId: id || null,
            message: `${modeLabel} Depth track ${slotLabel} has no logical Depth source.`
          });
        } else if (trackIndex >= depthTrackCount) {
          addRowIssue(rowIssues, index, {
            code: 'depth_track_index_range',
            field: 'params.track_index',
            slotId: id || null,
            message: `${modeLabel} Depth track ${slotLabel} references index ${trackIndex}; available indexes are 0..${depthTrackCount - 1}.`
          });
        }
      }
    }

    if (renderer === 'annotations' && params) {
      const setId = String(params.set_id || '').trim();
      if (!setId) {
        addRowIssue(rowIssues, index, {
          code: 'annotation_set_required',
          field: 'params.set_id',
          slotId: id || null,
          message: `${modeLabel} Annotation track ${slotLabel} requires an annotation set.`
        });
      } else if (annotationInventory && !annotationInventory.has(setId)) {
        addRowIssue(rowIssues, index, {
          code: 'annotation_set_unknown',
          field: 'params.set_id',
          slotId: id || null,
          message: `${modeLabel} Annotation track ${slotLabel} references unknown set '${setId}'.`
        });
      }

      if (
        Object.prototype.hasOwnProperty.call(params, 'marks') &&
        params.marks !== null &&
        (
          !Array.isArray(params.marks) ||
          params.marks.some((mark) => !ANNOTATION_MARKS.has(normalizedString(mark)))
        )
      ) {
        addRowIssue(rowIssues, index, {
          code: 'annotation_marks',
          field: 'params.marks',
          slotId: id || null,
          message: `${modeLabel} Annotation track ${slotLabel} has invalid mark filters.`
        });
      }
      [
        ['lane_gap_px', 'annotation_lane_gap', 'lane gap'],
        ['padding_px', 'annotation_padding', 'padding']
      ].forEach(([field, code, label]) => {
        if (!Object.prototype.hasOwnProperty.call(params, field)) return;
        if (Number.isNaN(nonnegativeFiniteNumber(params[field]))) {
          addRowIssue(rowIssues, index, {
            code,
            field: `params.${field}`,
            slotId: id || null,
            message: `${modeLabel} Annotation track ${slotLabel} ${label} must be a nonnegative finite number.`
          });
        }
      });
      if (
        Object.prototype.hasOwnProperty.call(params, 'cover_anchor') &&
        typeof params.cover_anchor !== 'boolean'
      ) {
        addRowIssue(rowIssues, index, {
          code: 'annotation_cover_anchor',
          field: 'params.cover_anchor',
          slotId: id || null,
          message: `${modeLabel} Annotation track ${slotLabel} Cover anchor must be true or false.`
        });
      }
      if (
        Object.prototype.hasOwnProperty.call(params, 'overflow') &&
        !['error', 'compress', 'clip'].includes(normalizedString(params.overflow))
      ) {
        addRowIssue(rowIssues, index, {
          code: 'annotation_overflow',
          field: 'params.overflow',
          slotId: id || null,
          message: `${modeLabel} Annotation track ${slotLabel} has invalid overflow behavior.`
        });
      }
      const layer = normalizedOptionalString(params.layer) ?? 'foreground';
      if (!['foreground', 'underlay'].includes(layer)) {
        addRowIssue(rowIssues, index, {
          code: 'annotation_layer',
          field: 'params.layer',
          slotId: id || null,
          message: `${modeLabel} Annotation track ${slotLabel} has an invalid layer.`
        });
      }
      const anchorId = String(params.anchor_slot || '').trim();
      if (side === 'overlay' && !anchorId) {
        addRowIssue(rowIssues, index, {
          code: 'annotation_anchor_required',
          field: 'params.anchor_slot',
          slotId: id || null,
          message: `${modeLabel} Annotation track ${slotLabel} requires an Overlay anchor.`
        });
      } else if (side !== 'overlay' && anchorId) {
        addRowIssue(rowIssues, index, {
          code: 'annotation_anchor_without_overlay',
          field: 'params.anchor_slot',
          slotId: id || null,
          message: `${modeLabel} Annotation track ${slotLabel} has an anchor but is not Overlay.`
        });
      }
    }

    if (renderer === 'sequence_conservation' && params) {
      if (String(params.managed || '').trim() !== CONSERVATION_MANAGER) {
        addRowIssue(rowIssues, index, {
          code: 'conservation_unmanaged',
          field: 'params.managed',
          slotId: id || null,
          message: `Circular comparison track ${slotLabel} is not owned by Comparison Series.`
        });
      }
      const sourceKey = String(params.series_key || '').trim();
      const sourceIndex = Number(params.source_index);
      const keyValid = sourceKey && conservation.keys.has(sourceKey);
      const indexValid = (
        Number.isSafeInteger(sourceIndex) &&
        sourceIndex >= 0 &&
        conservation.indexes.has(sourceIndex) &&
        conservation.indexesByKey.get(sourceKey)?.has(sourceIndex)
      );
      if (!keyValid || !indexValid) {
        addRowIssue(rowIssues, index, {
          code: 'conservation_source_missing',
          field: sourceKey ? 'params.source_index' : 'params.series_key',
          slotId: id || null,
          message: `Circular comparison track ${slotLabel} is not bound to an active Comparison Series source.`
        });
      }
    }

    const normalizedSlot = {
      ...slot,
      id,
      renderer,
      enabled: true,
      side,
      z: Number.isSafeInteger(Number(slot.z)) ? Number(slot.z) : 0,
      params: params || {}
    };
    normalizedEnabled.push(normalizedSlot);
    rowIndexByNormalizedSlot.set(normalizedSlot, index);
  });

  if (mode === 'linear') {
    const featureIndexes = enabledIndexes.filter(
      (index) => normalizedString(draftSlots[index]?.renderer) === 'features'
    );
    featureIndexes.slice(1).forEach((index) => addRowIssue(rowIssues, index, {
      code: 'features_multiple',
      field: 'renderer',
      slotId: String(draftSlots[index]?.id || '').trim() || null,
      message: 'Linear Custom Track Slots support only one enabled Features row.'
    }));
  }
  if (featureUnderlaysVisible(visibleFeatureUnderlays)) {
    const featureCount = enabledIndexes.filter(
      (index) => normalizedString(draftSlots[index]?.renderer) === 'features'
    ).length;
    if (featureCount !== 1) {
      globalIssues.push(globalIssue({
        code: 'feature_underlay_features_count',
        field: 'slots',
        message: `Visible feature underlays require exactly one enabled ${modeLabel} Features row.`
      }));
    }
  }

  const enabledById = new Map(
    normalizedEnabled
      .filter((slot) => slot.id)
      .map((slot) => [slot.id, slot])
  );
  normalizedEnabled.forEach((slot) => {
    if (slot.renderer !== 'annotations' || slot.side !== 'overlay') return;
    const rowIndex = rowIndexByNormalizedSlot.get(slot) ?? 0;
    const anchorId = String(slot.params.anchor_slot || '').trim();
    if (!anchorId) return;
    const anchor = enabledById.get(anchorId);
    if (!anchor) {
      addRowIssue(rowIssues, rowIndex, {
        code: 'annotation_anchor_unknown',
        field: 'params.anchor_slot',
        slotId: slot.id || null,
        message: `${modeLabel} Annotation track '${slot.id}' references missing or disabled anchor '${anchorId}'.`
      });
      return;
    }
    if (INELIGIBLE_ANCHOR_RENDERERS.has(anchor.renderer)) {
      addRowIssue(rowIssues, rowIndex, {
        code: 'annotation_anchor_ineligible',
        field: 'params.anchor_slot',
        slotId: slot.id || null,
        message: `${modeLabel} Annotation track '${slot.id}' anchor '${anchorId}' has no eligible drawable band.`
      });
      return;
    }
    const layer = normalizedString(slot.params.layer || 'foreground');
    if (layer === 'underlay' && slot.z >= anchor.z) {
      addRowIssue(rowIssues, rowIndex, {
        code: 'annotation_underlay_z',
        field: 'z',
        slotId: slot.id || null,
        message: `${modeLabel} Annotation underlay '${slot.id}' must have z below anchor '${anchorId}'.`
      });
    } else if (layer === 'foreground' && slot.z <= anchor.z) {
      addRowIssue(rowIssues, rowIndex, {
        code: 'annotation_foreground_z',
        field: 'z',
        slotId: slot.id || null,
        message: `${modeLabel} Annotation foreground '${slot.id}' must have z above anchor '${anchorId}'.`
      });
    }
  });

  if (emittedAxisIndex !== null) {
    normalizedEnabled.forEach((slot, slotIndex) => {
      if (!supportedRenderers.has(slot.renderer)) return;
      const rowIndex = rowIndexByNormalizedSlot.get(slot);
      if (!Number.isInteger(rowIndex)) return;
      let conflictsWithAxis = false;
      if (mode === 'linear') {
        const derivedSide = slotIndex < emittedAxisIndex ? 'above' : 'below';
        if (slot.side === 'overlay') {
          conflictsWithAxis = (
            slot.renderer === 'features' && slotIndex !== emittedAxisIndex
          );
        } else {
          conflictsWithAxis = slot.side !== derivedSide;
        }
      } else {
        const derivedSide = slotIndex < emittedAxisIndex ? 'outside' : 'inside';
        const featureLane = slot.renderer === 'features'
          ? normalizedOptionalString(slot.params.lane_direction ?? slot.params.lanes)
          : null;
        const isOverlayException = (
          (slot.renderer === 'features' && featureLane === 'split' && slot.side === 'overlay') ||
          (['ticks', 'annotations'].includes(slot.renderer) && slot.side === 'overlay')
        );
        const featureLaneConflictsWithAxis = (
          slot.renderer === 'features' &&
          ['inside', 'outside'].includes(featureLane) &&
          featureLane !== derivedSide
        );
        conflictsWithAxis = !isOverlayException && (
          slot.side !== derivedSide || featureLaneConflictsWithAxis
        );
      }
      if (!conflictsWithAxis) return;
      addRowIssue(rowIssues, rowIndex, {
        code: 'axis_side_conflict',
        field: 'side',
        slotId: slot.id || null,
        message: `${modeLabel} track '${slot.id}' side conflicts with the enabled Axis boundary.`
      });
    });
  }

  return {
    draftSlots,
    enabledSlots: normalizedEnabled,
    emittedAxisIndex,
    rowIssues,
    globalIssues
  };
};

/**
 * Validate the Web Custom Track draft without mutating it or throwing.
 *
 * Python remains authoritative for the complete typed schema. This validator
 * intentionally covers only constraints that the Web editor can create and
 * explain locally.
 */
export const validateCustomTrackPlan = (options = {}) => collectCustomTrackIssues({
  mode: normalizedString(options.mode),
  slots: options.slots,
  axisIndex: options.axisIndex,
  trackType: options.trackType,
  depthTrackCount: options.depthTrackCount,
  annotationSetIds: options.annotationSetIds,
  visibleFeatureUnderlays: options.visibleFeatureUnderlays,
  conservationSeries: options.conservationSeries
});

export class CustomTrackPlanValidationError extends Error {
  constructor(issues) {
    const issueList = Array.isArray(issues) ? issues : [];
    const first = issueList[0];
    const remaining = Math.max(0, issueList.length - 1);
    super(
      first
        ? `${first.message}${remaining ? ` (${remaining} more track issue${remaining === 1 ? '' : 's'}.)` : ''}`
        : 'Custom Track Slots are invalid.'
    );
    this.name = 'CustomTrackPlanValidationError';
    this.code = 'CUSTOM_TRACK_PLAN_INVALID';
    this.issues = issueList;
  }
}

export const customTrackPlanIssues = (plan) => {
  if (!plan || typeof plan !== 'object') return [];
  const rowIssues = plan.rowIssues instanceof Map
    ? Array.from(plan.rowIssues.values()).flat()
    : [];
  return [...rowIssues, ...(Array.isArray(plan.globalIssues) ? plan.globalIssues : [])];
};

export const assertValidCustomTrackPlan = (plan) => {
  const issues = customTrackPlanIssues(plan);
  if (issues.length > 0) throw new CustomTrackPlanValidationError(issues);
  return plan;
};
