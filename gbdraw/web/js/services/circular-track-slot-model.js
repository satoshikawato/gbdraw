// @ts-check
// The Circular track-slot model: renderer vocabulary, default stacks, slot
// normalization and legacy migration, axis and feature placement, and the
// request payload, as the request, Session and Gallery services read them. The
// Circular track-slot editor (app/circular-track-slots.js) edits it.
import { parseDepthTrackIndexIdentity } from './depth-track-state.js';
import { resolveColorToHex } from '../utils/color-utils.js';
import { normalizeOptionalText } from './track-slot-display.js';
import { parseOptionalCircularScalar, parseOptionalPixel } from './track-slot-validation.js';
import { diagnosticError } from '../utils/error-normalization.js';

const SUPPORTED_RENDERERS = [
  'features',
  'ticks',
  'dinucleotide_content',
  'dinucleotide_skew',
  'depth',
  'sequence_conservation',
  'annotations',
  'spacer'
];

export const DEFAULT_SLOT_IDS = {
  features: 'features',
  ticks: 'ticks',
  dinucleotide_content: 'gc_content',
  dinucleotide_skew: 'gc_skew',
  depth: 'depth',
  sequence_conservation: 'conservation',
  annotations: 'annotations',
  spacer: 'spacer'
};

const TICK_LABEL_LAYOUTS = [
  'label_out_tick_in',
  'label_in_tick_out',
  'tick_only',
  'label_only'
];
export const DEFAULT_TICK_LABEL_LAYOUT = 'label_out_tick_in';
const DEFAULT_INNER_TICK_LABEL_LAYOUT = 'label_in_tick_out';
const OBSOLETE_CIRCULAR_TRACK_SLOT_KEYS = new Set([
  'spacing',
  'strict',
  'compress',
  'reserve'
]);
export const CIRCULAR_TRACK_PRESETS = ['tuckin', 'middle', 'spreadout'];

export const normalizeCircularTrackPreset = (value, fallback = 'tuckin') => {
  const text = String(value || fallback).trim().toLowerCase();
  return CIRCULAR_TRACK_PRESETS.includes(text) ? text : fallback;
};

export const laneDirectionForPreset = (preset) => {
  const normalized = normalizeCircularTrackPreset(preset);
  if (normalized === 'middle') return 'split';
  if (normalized === 'spreadout') return 'outside';
  return 'inside';
};

export const normalizeLaneDirection = (value, fallback = 'inside') => {
  const text = String(value || fallback).trim().toLowerCase();
  return ['inside', 'outside', 'split'].includes(text) ? text : fallback;
};

export const sideForLaneDirection = (value) => {
  const lane = normalizeLaneDirection(value);
  return lane === 'split' ? 'overlay' : lane;
};

export const laneDirectionForSide = (value) => {
  const side = normalizePlacement(value);
  if (side === 'outside') return 'outside';
  if (side === 'overlay') return 'split';
  return 'inside';
};

export const normalizeNt = (value, fallback = 'GC') => {
  const text = String(value || '').trim().toUpperCase();
  return text || fallback;
};

const cleanToken = (value, fallback) => {
  const text = String(value || '').trim();
  return text || fallback;
};

export const normalizeColorParam = (value) => {
  const text = normalizeOptionalText(value);
  if (text === null) return null;
  return resolveColorToHex(text);
};

const normalizeSkewColorParams = (params) => {
  if (normalizeOptionalText(params.positive_color) === null && normalizeOptionalText(params.high_color) !== null) {
    params.positive_color = params.high_color;
  }
  if (normalizeOptionalText(params.negative_color) === null && normalizeOptionalText(params.low_color) !== null) {
    params.negative_color = params.low_color;
  }
  delete params.high_color;
  delete params.low_color;
  const positiveColor = normalizeColorParam(params.positive_color);
  if (positiveColor === null) delete params.positive_color;
  else params.positive_color = positiveColor;
  const negativeColor = normalizeColorParam(params.negative_color);
  if (negativeColor === null) delete params.negative_color;
  else params.negative_color = negativeColor;
  return params;
};

const normalizeGapText = (value, field) => {
  try {
    const numeric = parseOptionalPixel(value, `Circular track slot ${field}`, { allowZero: true });
    return numeric === null ? null : String(numeric);
  } catch {
    return value; // Invalid drafts remain visible and cannot become auto/zero.
  }
};

export const normalizePlacement = (value, fallback = 'inside') => {
  const text = String(value || fallback).trim().toLowerCase();
  return ['inside', 'outside', 'overlay'].includes(text) ? text : fallback;
};

export const normalizeOptionalPlacement = (value) => {
  if (value === null || value === undefined || value === '') return null;
  return normalizePlacement(value);
};

export const normalizeTrackIndex = (value) => {
  if (value === null || value === undefined || value === '') return null;
  const numeric = Number(value);
  if (!Number.isInteger(numeric) || numeric < 0) return null;
  return numeric;
};

export const normalizeTickLabelLayout = (value, fallback = DEFAULT_TICK_LABEL_LAYOUT) => {
  const text = String(value || fallback).trim().toLowerCase();
  return TICK_LABEL_LAYOUTS.includes(text) ? text : fallback;
};

export const tickLabelLayoutFromSides = (labelSide, tickSide, fallback = DEFAULT_TICK_LABEL_LAYOUT) => {
  const label = String(labelSide || '').trim().toLowerCase();
  const tick = String(tickSide || '').trim().toLowerCase();
  if (label === 'outside' && tick === 'inside') return 'label_out_tick_in';
  if (label === 'inside' && tick === 'outside') return 'label_in_tick_out';
  if (label === 'none' && ['inside', 'outside', 'both'].includes(tick)) return 'tick_only';
  if (['inside', 'outside'].includes(label) && tick === 'none') return 'label_only';
  return fallback;
};

const defaultPresetTickLabelLayout = () => DEFAULT_INNER_TICK_LABEL_LAYOUT;

const isAutoOrientableTickLabelLayout = (value) => {
  const layout = normalizeTickLabelLayout(value);
  return [DEFAULT_TICK_LABEL_LAYOUT, DEFAULT_INNER_TICK_LABEL_LAYOUT].includes(layout);
};

const syncDefaultTickLayoutsForFeatureRelation = (slots) => {
  const featureIndex = slots.findIndex((slot) => slot?.enabled !== false && slot?.renderer === 'features');
  if (featureIndex < 0) return;
  slots.forEach((slot, index) => {
    if (!slot || slot.enabled === false || slot.renderer !== 'ticks') return;
    slot.params = cloneParams(slot.params);
    if (!isAutoOrientableTickLabelLayout(slot.params.tick_label_layout)) return;
    slot.params.tick_label_layout = index > featureIndex
      ? DEFAULT_INNER_TICK_LABEL_LAYOUT
      : DEFAULT_TICK_LABEL_LAYOUT;
  });
};

export const cloneParams = (params = {}) => {
  if (!params || typeof params !== 'object' || Array.isArray(params)) return {};
  return { ...params };
};

const obsoleteCircularTrackSlotKey = (source) => {
  for (const key of Object.keys(source || {})) {
    if (OBSOLETE_CIRCULAR_TRACK_SLOT_KEYS.has(String(key).toLowerCase())) {
      return key;
    }
  }
  for (const key of Object.keys(source?.params || {})) {
    if (OBSOLETE_CIRCULAR_TRACK_SLOT_KEYS.has(String(key).toLowerCase())) {
      return `params.${key}`;
    }
  }
  return null;
};

const assertCurrentCircularTrackSlotShape = (source) => {
  const obsoleteKey = obsoleteCircularTrackSlotKey(source);
  if (!obsoleteKey) return;
  throw new Error(
    `Circular track slot field '${obsoleteKey}' is obsolete. ` +
    'Use inner_gap_px and outer_gap_px for physical gaps.'
  );
};

export const migrateLegacyCircularTrackSlot = (slot) => {
  if (!slot || typeof slot !== 'object' || Array.isArray(slot)) return slot;
  const source = { ...slot };
  const params = cloneParams(source.params);
  const topLevelSpacing = Object.prototype.hasOwnProperty.call(source, 'spacing')
    ? source.spacing
    : undefined;
  const paramSpacing = Object.prototype.hasOwnProperty.call(params, 'spacing')
    ? params.spacing
    : undefined;
  const legacySpacing = normalizeOptionalText(topLevelSpacing) !== null
    ? topLevelSpacing
    : paramSpacing;

  Object.keys(source).forEach((key) => {
    if (OBSOLETE_CIRCULAR_TRACK_SLOT_KEYS.has(key.toLowerCase())) {
      delete source[key];
    }
  });
  Object.keys(params).forEach((key) => {
    if (OBSOLETE_CIRCULAR_TRACK_SLOT_KEYS.has(key.toLowerCase())) {
      delete params[key];
    }
  });
  source.params = params;
  if (normalizeOptionalText(legacySpacing) !== null) {
    if (normalizeOptionalText(source.inner_gap_px ?? source.innerGapPx) === null) {
      source.inner_gap_px = legacySpacing;
    }
    if (normalizeOptionalText(source.outer_gap_px ?? source.outerGapPx) === null) {
      source.outer_gap_px = legacySpacing;
    }
  }
  return source;
};

export const migrateLegacyCircularTrackSlotSpec = (spec) => {
  const text = String(spec || '').trim();
  const atIndex = text.indexOf('@');
  if (atIndex < 0) return text;

  const head = text.slice(0, atIndex).trim();
  const retained = [];
  /** @type {string | null} */
  let legacySpacing = null;
  let hasInnerGap = false;
  let hasOuterGap = false;
  text.slice(atIndex + 1).split(',').forEach((token) => {
    const equalsIndex = token.indexOf('=');
    if (equalsIndex < 0) return;
    const key = token.slice(0, equalsIndex).trim();
    const normalizedKey = key.toLowerCase();
    const rawValue = token.slice(equalsIndex + 1).trim();
    if (normalizedKey === 'spacing') {
      legacySpacing = rawValue;
      return;
    }
    if (OBSOLETE_CIRCULAR_TRACK_SLOT_KEYS.has(normalizedKey)) return;
    if (normalizedKey === 'inner_gap_px') hasInnerGap = true;
    if (normalizedKey === 'outer_gap_px') hasOuterGap = true;
    retained.push(`${key}=${rawValue}`);
  });
  if (legacySpacing !== null && legacySpacing !== '') {
    if (!hasInnerGap) retained.push(`inner_gap_px=${legacySpacing}`);
    if (!hasOuterGap) retained.push(`outer_gap_px=${legacySpacing}`);
  }
  return retained.length > 0 ? `${head}@${retained.join(',')}` : head;
};

export const normalizeSlotSide = (value) => normalizeOptionalPlacement(value);

export const resolveCircularTrackFeaturePlacement = (slot, preset = 'tuckin') => {
  const params = cloneParams(slot?.params);
  const rawLane = normalizeOptionalText(params.lane_direction ?? params.lanes);
  const explicitSide = normalizeOptionalPlacement(slot?.side);
  const laneDirection = rawLane === null
    ? (
        explicitSide === null
          ? laneDirectionForPreset(preset)
          : laneDirectionForSide(explicitSide)
      )
    : normalizeLaneDirection(rawLane);
  return {
    laneDirection,
    side: sideForLaneDirection(laneDirection)
  };
};

const featureLaneForSlot = (slot, preset = 'tuckin') => (
  resolveCircularTrackFeaturePlacement(slot, preset).laneDirection
);

const syncTickParamsForPlacement = (slot, side) => {
  if (!slot || slot.renderer !== 'ticks') return;
  void side;
  slot.params = cloneParams(slot.params);
  slot.params.tick_label_layout = normalizeTickLabelLayout(slot.params.tick_label_layout);
  delete slot.params.label_side;
  delete slot.params.tick_side;
};

export const syncSlotPlacementFromSide = (slot, side) => {
  if (!slot) return;
  const placement = normalizePlacement(side);
  slot.side = placement;
  if (slot.renderer === 'features') {
    slot.params = cloneParams(slot.params);
    slot.params.lane_direction = laneDirectionForSide(placement);
  } else if (slot.renderer === 'ticks') {
    syncTickParamsForPlacement(slot, placement);
  }
};

export const effectiveSlotPlacement = (slot, preset = 'tuckin') => {
  if (!slot) return 'inside';
  if (slot.renderer === 'features') {
    return sideForLaneDirection(featureLaneForSlot(slot, preset));
  }
  return normalizePlacement(slot.side, 'inside');
};

export const inferLegacyAxisIndexFromFeature = (slots, preset = 'tuckin') => {
  const onAxisIndex = slots.findIndex((slot) => slot?.renderer !== 'annotations' && effectiveSlotPlacement(slot, preset) === 'overlay');
  if (onAxisIndex >= 0) return onAxisIndex;
  const featureIndex = slots.findIndex((slot) => slot?.enabled !== false && slot?.renderer === 'features');
  if (featureIndex < 0) {
    const firstInside = slots.findIndex((slot) => normalizePlacement(slot?.side, 'inside') !== 'outside');
    return firstInside < 0 ? slots.length : firstInside;
  }
  const featureSide = sideForLaneDirection(featureLaneForSlot(slots[featureIndex], preset));
  if (featureSide === 'outside') return featureIndex + 1;
  return featureIndex;
};

export const clampCircularTrackAxisIndex = (value, slotCount) => {
  const length = Math.max(0, Number(slotCount) || 0);
  const numeric = Number(value);
  if (!Number.isInteger(numeric)) return null;
  return Math.min(Math.max(numeric, 0), length);
};

export const syncSlotsFromAxisIndex = (slots, axisIndex, preset = 'tuckin') => {
  const axis = clampCircularTrackAxisIndex(axisIndex, slots.length);
  const resolvedAxis = axis === null ? inferLegacyAxisIndexFromFeature(slots, preset) : axis;
  slots.forEach((slot, index) => {
    if (!slot) return;
    if (slot.renderer === 'annotations' && effectiveSlotPlacement(slot, preset) === 'overlay') return;
    if (effectiveSlotPlacement(slot, preset) === 'overlay') {
      syncSlotPlacementFromSide(slot, 'overlay');
      return;
    }
    if (slot.renderer === 'features' && featureLaneForSlot(slot, preset) === 'split') {
      syncSlotPlacementFromSide(slot, 'overlay');
      return;
    }
    const side = index < resolvedAxis ? 'outside' : 'inside';
    syncSlotPlacementFromSide(slot, side);
  });
  syncDefaultTickLayoutsForFeatureRelation(slots);
  return resolvedAxis;
};

export const enforceSingleOnAxisSlot = (slots, axisIndex, preset = 'tuckin') => {
  const onAxisIndices = slots
    .map((slot, index) => (
      slot?.renderer !== 'annotations' && effectiveSlotPlacement(slot, preset) === 'overlay' ? index : null
    ))
    .filter((index) => Number.isInteger(index));
  if (onAxisIndices.length === 0) {
    return clampCircularTrackAxisIndex(axisIndex, slots.length) ?? inferLegacyAxisIndexFromFeature(slots, preset);
  }

  const clampedAxis = clampCircularTrackAxisIndex(axisIndex, slots.length);
  const keepIndex = onAxisIndices.includes(clampedAxis) ? clampedAxis : onAxisIndices[0];
  onAxisIndices.forEach((index) => {
    if (index === keepIndex) return;
    syncSlotPlacementFromSide(slots[index], index < keepIndex ? 'outside' : 'inside');
  });
  return keepIndex;
};

/** @param {number | null} [axisIndex] */
export const applyCircularTrackOrderPlacements = (slots, defaultNt = 'GC', preset = 'tuckin', axisIndex = null) => {
  const normalized = normalizeCircularTrackSlots(slots, defaultNt, preset);
  const resolvedAxis = syncSlotsFromAxisIndex(normalized, axisIndex, preset);
  enforceSingleOnAxisSlot(normalized, resolvedAxis, preset);
  return normalized;
};

/** @param {{ id: string, renderer: string, enabled?: boolean, width?: any, radius?: any, inner_gap_px?: any, outer_gap_px?: any, side?: string | null, z?: number, params?: Record<string, any> }} slot */
export const makeSlot = ({
  id,
  renderer,
  enabled = true,
  width = null,
  radius = null,
  inner_gap_px = null,
  outer_gap_px = null,
  side = null,
  z = 0,
  params = {}
}) => ({
  id,
  renderer,
  enabled,
  width,
  radius,
  inner_gap_px,
  outer_gap_px,
  side: normalizeSlotSide(side),
  z,
  params: cloneParams(params)
});

const paramsMatchExactly = (params, expected = {}) => {
  const actualEntries = Object.entries(cloneParams(params))
    .filter(([, value]) => normalizeOptionalText(value) !== null)
    .map(([key, value]) => [String(key), normalizeOptionalText(value)]);
  /** @type {[string, string | null][]} */
  const expectedEntries = Object.entries(expected)
    .filter(([, value]) => normalizeOptionalText(value) !== null)
    .map(([key, value]) => [String(key), normalizeOptionalText(value)]);
  if (actualEntries.length !== expectedEntries.length) return false;
  const actualMap = Object.fromEntries(actualEntries);
  return expectedEntries.every(([key, value]) => actualMap[key] === value);
};

const hasBlankSlotGeometry = (source) =>
  normalizeOptionalText(source.width) === null &&
  normalizeOptionalText(source.radius) === null &&
  normalizeOptionalText(source.inner_gap_px) === null &&
  normalizeOptionalText(source.outer_gap_px) === null;

const hasDefaultSlotFlags = (source) =>
  source.enabled !== false &&
  Number(source.z || 0) === 0;

const isLegacyDefaultWebSlotShape = (source, renderer, defaultNt = 'GC', preset = 'tuckin') => {
  if (!source || typeof source !== 'object' || Array.isArray(source)) return false;
  if (!hasBlankSlotGeometry(source) || !hasDefaultSlotFlags(source)) return false;

  const normalizedPreset = normalizeCircularTrackPreset(preset);
  const normalizedId = String(source.id || '').trim();
  const side = normalizeOptionalPlacement(source.side);
  const params = cloneParams(source.params);

  if (renderer === 'features' && normalizedId === 'features') {
    const laneDirection = laneDirectionForPreset(normalizedPreset);
    return (
      side === sideForLaneDirection(laneDirection) &&
      paramsMatchExactly(params, { lane_direction: laneDirection })
    );
  }

  if (renderer === 'ticks' && normalizedId === 'ticks') {
    const tickLayout = normalizeTickLabelLayout(params.tick_label_layout);
    return (
      side === 'inside' &&
      [DEFAULT_TICK_LABEL_LAYOUT, defaultPresetTickLabelLayout()].includes(tickLayout) &&
      paramsMatchExactly(
        {
          ...params,
          tick_label_layout: DEFAULT_TICK_LABEL_LAYOUT
        },
        {
          tick_label_layout: DEFAULT_TICK_LABEL_LAYOUT,
          preset: normalizedPreset
        }
      )
    );
  }

  if (renderer === 'depth' && normalizedId === 'depth') {
    return side === 'inside' && paramsMatchExactly(params, {});
  }

  if (renderer === 'dinucleotide_content' && normalizedId === 'gc_content') {
    return (
      side === 'inside' &&
      paramsMatchExactly(params, { nt: normalizeNt(defaultNt) })
    );
  }

  if (renderer === 'dinucleotide_skew' && normalizedId === 'gc_skew') {
    return (
      side === 'inside' &&
      paramsMatchExactly(params, { nt: normalizeNt(defaultNt) })
    );
  }

  return false;
};

/**
 * @typedef {object} DefaultCircularTrackSlotsOptions
 * @property {string} [nt] dinucleotide of the GC tracks
 * @property {boolean} [showDepth]
 * @property {number} [depthTrackCount]
 * @property {boolean} [showGc]
 * @property {boolean} [showSkew]
 * @property {boolean} [showTicks]
 * @property {string} [preset] circular track preset (`tuckin`, `middle`, `spreadout`)
 */

/**
 * @param {DefaultCircularTrackSlotsOptions} [options]
 * @returns {Record<string, any>[]}
 */
export const createDefaultCircularTrackSlots = ({
  nt = 'GC',
  showDepth = false,
  depthTrackCount = 1,
  showGc = true,
  showSkew = true,
  showTicks = true,
  preset = 'tuckin'
} = {}) => {
  void nt;
  void preset;
  const slots = [
    makeSlot({
      id: 'features',
      renderer: 'features'
    })
  ];
  if (showTicks) {
    slots.push(makeSlot({
      id: 'ticks',
      renderer: 'ticks',
      params: {
        tick_label_layout: defaultPresetTickLabelLayout()
      }
    }));
  }
  if (showDepth) {
    const count = Math.max(1, Number(depthTrackCount) || 1);
    if (count === 1) {
      slots.push(makeSlot({ id: 'depth', renderer: 'depth' }));
    } else {
      for (let index = 0; index < count; index += 1) {
        slots.push(makeSlot({
          id: `depth_${index + 1}`,
          renderer: 'depth',
          params: {
            track_index: index
          }
        }));
      }
    }
  }
  if (showGc) {
    slots.push(makeSlot({
      id: 'gc_content',
      renderer: 'dinucleotide_content'
    }));
  }
  if (showSkew) {
    slots.push(makeSlot({
      id: 'gc_skew',
      renderer: 'dinucleotide_skew'
    }));
  }
  return slots;
};

// Shortcut -> diagnostic field (the Web control's adv key).
const CIRCULAR_GEOMETRY_SHORTCUT_FIELDS = Object.freeze({
  featureWidth: 'feature_width_circular',
  depthWidth: 'depth_width_circular',
  gcContentWidth: 'gc_content_width_circular',
  gcContentRadius: 'gc_content_radius_circular',
  gcSkewWidth: 'gc_skew_width_circular',
  gcSkewRadius: 'gc_skew_radius_circular'
});

export const normalizeCircularGeometryShortcuts = (values = {}) => {
  const normalized = {};
  for (const [field, diagnosticField] of Object.entries(CIRCULAR_GEOMETRY_SHORTCUT_FIELDS)) {
    const raw = values?.[field];
    if (raw === null || raw === undefined || raw === '') {
      normalized[field] = null;
      continue;
    }
    const numeric = typeof raw === 'boolean' ? NaN : Number(raw);
    if (!Number.isFinite(numeric) || numeric <= 0) {
      throw diagnosticError('INPUT_INVALID', { field: diagnosticField, reason: 'POSITIVE_OR_AUTO' });
    }
    normalized[field] = numeric;
  }
  return normalized;
};

export const hasCircularGeometryShortcuts = (values = {}) => (
  Object.values(normalizeCircularGeometryShortcuts(values))
    .some((value) => value !== null)
);

export const applyCircularGeometryShortcuts = (slots, values = {}) => {
  const geometry = normalizeCircularGeometryShortcuts(values);
  return (Array.isArray(slots) ? slots : []).map((slot) => {
    const next = {
      ...slot,
      params: cloneParams(slot?.params)
    };
    const renderer = String(next.renderer || '').trim();
    if (renderer === 'features' && geometry.featureWidth !== null) {
      next.width = `${geometry.featureWidth}px`;
    } else if (renderer === 'depth' && geometry.depthWidth !== null) {
      next.width = `${geometry.depthWidth}px`;
    } else if (renderer === 'dinucleotide_content') {
      if (geometry.gcContentWidth !== null) {
        next.width = `${geometry.gcContentWidth}px`;
      }
      if (geometry.gcContentRadius !== null) {
        next.radius = String(geometry.gcContentRadius);
      }
    } else if (renderer === 'dinucleotide_skew') {
      if (geometry.gcSkewWidth !== null) {
        next.width = `${geometry.gcSkewWidth}px`;
      }
      if (geometry.gcSkewRadius !== null) {
        next.radius = String(geometry.gcSkewRadius);
      }
    }
    return next;
  });
};

export const normalizeCircularTrackSlot = (slot, index = 0, defaultNt = 'GC', preset = 'tuckin') => {
  const source = slot && typeof slot === 'object' && !Array.isArray(slot) ? slot : {};
  assertCurrentCircularTrackSlotShape(source);
  const renderer = SUPPORTED_RENDERERS.includes(source.renderer) ? source.renderer : 'dinucleotide_skew';
  const fallbackId = DEFAULT_SLOT_IDS[renderer] || `slot_${index + 1}`;
  const inheritsPresetDefaults = isLegacyDefaultWebSlotShape(source, renderer, defaultNt, preset);
  const params = inheritsPresetDefaults ? {} : cloneParams(source.params);
  [
    'side',
    'r',
    'radius',
    'w',
    'width',
    'inner_gap_px',
    'outer_gap_px'
  ].forEach((key) => {
    delete params[key];
  });

  let side = inheritsPresetDefaults ? null : normalizeSlotSide(source.side);
  const radius = source.radius ?? null;
  const innerGapPx = normalizeGapText(source.inner_gap_px ?? source.innerGapPx, 'inner_gap_px');
  const outerGapPx = normalizeGapText(source.outer_gap_px ?? source.outerGapPx, 'outer_gap_px');

  if (renderer === 'dinucleotide_content' || renderer === 'dinucleotide_skew') {
    const nt = normalizeOptionalText(params.nt ?? params.dinucleotide);
    if (nt === null) {
      delete params.nt;
      delete params.dinucleotide;
    } else {
      params.nt = normalizeNt(nt);
      delete params.dinucleotide;
    }
    if (renderer === 'dinucleotide_skew') {
      normalizeSkewColorParams(params);
    }
  }
  if (renderer === 'ticks') {
    delete params['axis'];
    const legacyLayout = (
      normalizeOptionalText(params.tick_label_layout) === null
        ? tickLabelLayoutFromSides(params.label_side, params.tick_side)
        : params.tick_label_layout
    );
    params.tick_label_layout = normalizeTickLabelLayout(legacyLayout);
    delete params.label_side;
    delete params.tick_side;
    if (normalizeOptionalText(params.preset) === null) {
      delete params.preset;
    } else {
      params.preset = normalizeCircularTrackPreset(params.preset);
    }
  }
  if (renderer === 'features') {
    const rawLaneDirection = normalizeOptionalText(params.lane_direction ?? params.lanes);
    delete params.lanes;
    if (rawLaneDirection === null) {
      delete params.lane_direction;
      if (side !== null) {
        const placement = resolveCircularTrackFeaturePlacement({ side, params }, preset);
        params.lane_direction = placement.laneDirection;
        side = placement.side;
      }
    } else {
      const placement = resolveCircularTrackFeaturePlacement({
        side,
        params: { ...params, lane_direction: rawLaneDirection }
      }, preset);
      params.lane_direction = placement.laneDirection;
      side = placement.side;
    }
  }
  if (renderer === 'depth') {
    const parsedTrackIndex = normalizeTrackIndex(params.track_index);
    const idMatch = cleanToken(source.id, fallbackId).match(/^depth_(\d+)$/);
    if (parsedTrackIndex !== null) {
      params.track_index = parsedTrackIndex;
    } else if (idMatch) {
      params.track_index = Math.max(0, Number(idMatch[1]) - 1);
    } else {
      delete params.track_index;
    }
  }
  if (renderer === 'annotations') {
    params.set_id = String(params.set_id || '').trim();
    params.overflow = ['error', 'compress', 'clip'].includes(String(params.overflow || '').toLowerCase())
      ? String(params.overflow).toLowerCase()
      : 'error';
    params.show_labels = params.show_labels !== false && String(params.show_labels).toLowerCase() !== 'false';
    params.layer = String(params.layer || '').toLowerCase() === 'underlay' ? 'underlay' : 'foreground';
    if (Array.isArray(params.marks)) {
      params.marks = Array.from(new Set(
        params.marks
          .map((mark) => String(mark || '').trim().toLowerCase())
          .filter(Boolean)
      ));
      if (params.marks.length === 0) delete params.marks;
    }
    for (const [field, defaultValue] of [['lane_gap_px', 3], ['padding_px', 2]]) {
      if (params[field] === null || params[field] === undefined || params[field] === '') {
        delete params[field];
        continue;
      }
      const numeric = Number(params[field]);
      if (Number.isFinite(numeric) && numeric >= 0) {
        if (numeric === defaultValue) delete params[field];
        else params[field] = numeric;
      }
    }
    if (params.cover_anchor === true || String(params.cover_anchor).toLowerCase() === 'true') {
      params.cover_anchor = true;
    } else if (params.cover_anchor === false || String(params.cover_anchor).toLowerCase() === 'false') {
      delete params.cover_anchor;
    }
    if (side === null) side = 'outside';
  }

  return makeSlot({
    id: cleanToken(source.id, fallbackId),
    renderer,
    enabled: source.enabled !== false,
    width: source.width ?? null,
    radius,
    inner_gap_px: innerGapPx,
    outer_gap_px: outerGapPx,
    side,
    z: Number.isFinite(Number(source.z)) ? Number(source.z) : 0,
    params
  });
};

export const normalizeCircularTrackSlots = (slots, defaultNt = 'GC', preset = 'tuckin') => {
  const base = Array.isArray(slots)
    ? slots
    : createDefaultCircularTrackSlots({ nt: defaultNt, preset });
  return base.map((slot, index) => normalizeCircularTrackSlot(slot, index, defaultNt, preset));
};

export const parseCircularTrackSlotSpec = (spec, index = 0, defaultNt = 'GC', preset = 'tuckin') => {
  const text = String(spec || '').trim();
  const atIndex = text.indexOf('@');
  const head = (atIndex < 0 ? text : text.slice(0, atIndex)).trim();
  const separatorIndex = head.indexOf(':');
  const source = {
    id: separatorIndex < 0 ? '' : head.slice(0, separatorIndex).trim(),
    renderer: separatorIndex < 0 ? head : head.slice(separatorIndex + 1).trim(),
    enabled: true,
    params: {}
  };

  if (atIndex >= 0) {
    text.slice(atIndex + 1).split(',').forEach((token) => {
      const equalsIndex = token.indexOf('=');
      if (equalsIndex < 0) return;
      const key = token.slice(0, equalsIndex).trim();
      const rawValue = token.slice(equalsIndex + 1).trim();
      if (!key) return;
      if (OBSOLETE_CIRCULAR_TRACK_SLOT_KEYS.has(key.toLowerCase())) {
        throw new Error(
          `Circular track slot field '${key}' is obsolete. ` +
          'Use inner_gap_px and outer_gap_px for physical gaps.'
        );
      }
      const value = rawValue === 'true' ? true : (rawValue === 'false' ? false : rawValue);
      if (key === 'enabled') source.enabled = value !== false;
      else if (key === 'w' || key === 'width') source.width = rawValue;
      else if (key === 'r' || key === 'radius') source.radius = rawValue;
      else if (key === 'inner_gap_px') source.inner_gap_px = rawValue;
      else if (key === 'outer_gap_px') source.outer_gap_px = rawValue;
      else if (key === 'side') source.side = rawValue;
      else if (key === 'z') source.z = Number(rawValue);
      else source.params[key] = value;
    });
  }

  if (
    source.enabled !== false &&
    source.renderer === 'depth' &&
    Object.prototype.hasOwnProperty.call(source.params, 'track_index')
  ) {
    source.params.track_index = parseDepthTrackIndexIdentity(
      source.params.track_index,
      `Circular Depth slot '${source.id || `#${index + 1}`}' track_index`
    );
  }

  source.inner_gap_px = parseOptionalPixel(
    source.inner_gap_px, `Circular track '${source.id}' inner_gap_px`, { allowZero: true }
  );
  source.outer_gap_px = parseOptionalPixel(
    source.outer_gap_px, `Circular track '${source.id}' outer_gap_px`, { allowZero: true }
  );
  return normalizeCircularTrackSlot(source, index, defaultNt, preset);
};

export const parseCircularTrackSlotSpecs = (specs, defaultNt = 'GC', preset = 'tuckin') => (
  Array.isArray(specs)
    ? specs.map((spec, index) => parseCircularTrackSlotSpec(spec, index, defaultNt, preset))
    : []
);

const canonicalAnnotationParams = (params) => {
  const next = { ...params };
  if (!Array.isArray(next.marks) || next.marks.length === 0) {
    delete next.marks;
  } else {
    next.marks = Array.from(new Set(
      next.marks.map((mark) => String(mark || '').trim().toLowerCase()).filter(Boolean)
    ));
  }
  const laneGap = next.lane_gap_px === null || next.lane_gap_px === undefined || next.lane_gap_px === ''
    ? null
    : Number(next.lane_gap_px);
  if (laneGap === null || laneGap === 3) delete next.lane_gap_px;
  else next.lane_gap_px = laneGap;
  const padding = next.padding_px === null || next.padding_px === undefined || next.padding_px === ''
    ? null
    : Number(next.padding_px);
  if (padding === null || padding === 2) delete next.padding_px;
  else next.padding_px = padding;
  if (next.cover_anchor !== true) delete next.cover_anchor;
  if (String(next.overflow || '').trim().toLowerCase() === 'error') delete next.overflow;
  if (next.show_labels !== false) delete next.show_labels;
  if (String(next.layer || '').trim().toLowerCase() === 'foreground') delete next.layer;
  if (!String(next.anchor_slot || '').trim()) delete next.anchor_slot;
  return next;
};

const canonicalCircularParams = (slot) => {
  const params = cloneParams(slot?.params);
  Object.keys(params).forEach((key) => {
    if (params[key] === null || params[key] === undefined || key.startsWith('_')) {
      delete params[key];
    }
  });
  if (slot?.renderer === 'annotations') return canonicalAnnotationParams(params);
  return params;
};

/**
 * Encode one Web draft row as the canonical CircularTrackSlot object.
 *
 * This path is intentionally structured so nested annotation style overrides
 * and mark arrays survive request/session round trips.
 */
export const buildCircularTrackSlotPayload = (
  slot,
  defaultNt = 'GC',
  preset = 'tuckin'
) => {
  const normalized = normalizeCircularTrackSlot(slot, 0, defaultNt, preset);
  const params = canonicalCircularParams(normalized);
  let side = normalized.side;
  if (normalized.renderer === 'features') {
    const placement = resolveCircularTrackFeaturePlacement(normalized, preset);
    side = placement.side;
    params.lane_direction = placement.laneDirection;
  }
  return {
    kind: 'circularTrackSlot',
    id: normalized.id,
    renderer: normalized.renderer,
    enabled: normalized.enabled,
    side,
    radius: parseOptionalCircularScalar(normalized.radius, `Circular track '${normalized.id}' radius`),
    width: parseOptionalCircularScalar(normalized.width, `Circular track '${normalized.id}' width`),
    z: Number(normalized.z) || 0,
    params,
    innerGapPx: parseOptionalPixel(
      normalized.inner_gap_px,
      `Circular track '${normalized.id}' inner_gap_px`,
      { allowZero: true }
    ),
    outerGapPx: parseOptionalPixel(
      normalized.outer_gap_px,
      `Circular track '${normalized.id}' outer_gap_px`,
      { allowZero: true }
    )
  };
};

export { SUPPORTED_RENDERERS as CIRCULAR_TRACK_RENDERERS };
