// @ts-check
// The Circular track-slot model: renderer vocabulary, default stacks, slot
// normalization and legacy migration, axis and feature placement, and the
// request payload, as the request, Session and Gallery services read them. The
// Circular track-slot editor (app/circular-track-slots.js) edits it.
import { parseDepthTrackIndexIdentity } from './depth-track-state.js';
import { resolveColorToHex } from '../utils/color-utils.js';
import { normalizeOptionalText } from './track-slot-display.js';
import { canonicalAnnotationParams, normalizeAnnotationSlotParams, parseOptionalCircularScalar, parseOptionalPixel } from './track-slot-validation.js';
import { isUnsetValue, noteReplacedLeaf } from '../utils/field-reading.js';
import { diagnosticError } from '../utils/error-normalization.js';

/** @import { TrackSlotRow } from './track-slot-validation.js' */
/** @import { FieldLeaf } from '../utils/field-reading.js' */

/**
 * A Circular track-slot row as the editor and the request projection hold it.
 * @typedef {object} CircularTrackSlot
 * @property {string} id
 * @property {string} renderer
 * @property {boolean} enabled
 * @property {unknown} width A Width draft: text, a number, `{ value, unit }`, or null for Auto.
 * @property {unknown} radius A Radius draft, as the Width.
 * @property {unknown} inner_gap_px Text or a number of pixels; null for Auto.
 * @property {unknown} outer_gap_px As inner_gap_px.
 * @property {string | null} side
 * @property {number} z
 * @property {Record<string, unknown>} params
 */

/**
 * What the placement readers take of a row.
 * @typedef {object} SlotPlacementSource
 * @property {unknown} [enabled]
 * @property {unknown} [renderer]
 * @property {unknown} [side]
 * @property {unknown} [params]
 */

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

// The renderer names a slot may use beside the supported ones, as Python reads
// them on a spec string and on a slot object (gbdraw/tracks/circular.py
// `_RENDERER_ALIASES`, `_normalize_renderer`; OV-489).
/** @type {Readonly<Record<string, string>>} */
const RENDERER_ALIASES = Object.freeze({
  gc_content: 'dinucleotide_content',
  content: 'dinucleotide_content',
  gc_skew: 'dinucleotide_skew',
  skew: 'dinucleotide_skew',
  conservation: 'sequence_conservation',
  circular_comparison: 'sequence_conservation'
});
// The renderer a Circular slot names: the one change point of the alias lookup,
// for the slot normalizer and for Load's projection of a request's slot objects.
/** @param {unknown} name @returns {string} */
export const circularRendererName = (name) => {
  const named = String(name ?? '').trim().toLowerCase();
  return RENDERER_ALIASES[named] || named;
};

/** @type {Record<string, string>} */
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

/**
 * @param {unknown} value
 * @param {string} [fallback]
 * @returns {string}
 */
export const normalizeCircularTrackPreset = (value, fallback = 'tuckin') => {
  const text = String(value || fallback).trim().toLowerCase();
  return CIRCULAR_TRACK_PRESETS.includes(text) ? text : fallback;
};

/** @param {unknown} preset */
export const laneDirectionForPreset = (preset) => {
  const normalized = normalizeCircularTrackPreset(preset);
  if (normalized === 'middle') return 'split';
  if (normalized === 'spreadout') return 'outside';
  return 'inside';
};

/**
 * @param {unknown} value
 * @param {string} [fallback]
 * @returns {string}
 */
export const normalizeLaneDirection = (value, fallback = 'inside') => {
  const text = String(value || fallback).trim().toLowerCase();
  return ['inside', 'outside', 'split'].includes(text) ? text : fallback;
};

/** @param {unknown} value */
export const sideForLaneDirection = (value) => {
  const lane = normalizeLaneDirection(value);
  return lane === 'split' ? 'overlay' : lane;
};

/** @param {unknown} value */
export const laneDirectionForSide = (value) => {
  const side = normalizePlacement(value);
  if (side === 'outside') return 'outside';
  if (side === 'overlay') return 'split';
  return 'inside';
};

/**
 * @param {unknown} value
 * @param {string} [fallback]
 * @returns {string}
 */
export const normalizeNt = (value, fallback = 'GC') => {
  const text = String(value || '').trim().toUpperCase();
  return text || fallback;
};

/**
 * @param {unknown} value
 * @param {string} fallback
 * @returns {string}
 */
const cleanToken = (value, fallback) => {
  const text = String(value || '').trim();
  return text || fallback;
};

/** @param {unknown} value */
export const normalizeColorParam = (value) => {
  const text = normalizeOptionalText(value);
  if (text === null) return null;
  return resolveColorToHex(text);
};

/**
 * @param {Record<string, unknown>} params Edited in place.
 * @returns {Record<string, unknown>} The same object.
 */
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

/**
 * @param {unknown} value
 * @param {string} field
 */
const normalizeGapText = (value, field) => {
  try {
    const numeric = parseOptionalPixel(value, `Circular track slot ${field}`, { allowZero: true });
    return numeric === null ? null : String(numeric);
  } catch {
    return value; // Invalid drafts remain visible and cannot become auto/zero.
  }
};

/**
 * @param {unknown} value
 * @param {string} [fallback]
 * @returns {string}
 */
export const normalizePlacement = (value, fallback = 'inside') => {
  const text = String(value || fallback).trim().toLowerCase();
  return ['inside', 'outside', 'overlay'].includes(text) ? text : fallback;
};

/**
 * @param {unknown} value
 * @returns {string | null}
 */
export const normalizeOptionalPlacement = (value) => {
  if (value === null || value === undefined || value === '') return null;
  return normalizePlacement(value);
};

/**
 * @param {unknown} value
 * @returns {number | null}
 */
export const normalizeTrackIndex = (value) => {
  if (value === null || value === undefined || value === '') return null;
  const numeric = Number(value);
  if (!Number.isInteger(numeric) || numeric < 0) return null;
  return numeric;
};

/**
 * @param {unknown} value
 * @param {string} [fallback]
 * @returns {string}
 */
export const normalizeTickLabelLayout = (value, fallback = DEFAULT_TICK_LABEL_LAYOUT) => {
  const text = String(value || fallback).trim().toLowerCase();
  return TICK_LABEL_LAYOUTS.includes(text) ? text : fallback;
};

/**
 * @param {unknown} labelSide
 * @param {unknown} tickSide
 * @param {string} [fallback]
 * @returns {string}
 */
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

/** @param {unknown} value */
const isAutoOrientableTickLabelLayout = (value) => {
  const layout = normalizeTickLabelLayout(value);
  return [DEFAULT_TICK_LABEL_LAYOUT, DEFAULT_INNER_TICK_LABEL_LAYOUT].includes(layout);
};

/** @param {(CircularTrackSlot | null | undefined)[]} slots Edited in place. */
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

/**
 * @param {unknown} [params]
 * @returns {Record<string, unknown>} A copy; empty for a value that is no plain object.
 */
export const cloneParams = (params = {}) => {
  if (!params || typeof params !== 'object' || Array.isArray(params)) return {};
  return { ...params };
};

/**
 * @param {Record<string, unknown>} source
 * @returns {string | null}
 */
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

/** @param {Record<string, unknown>} source */
const assertCurrentCircularTrackSlotShape = (source) => {
  const obsoleteKey = obsoleteCircularTrackSlotKey(source);
  if (!obsoleteKey) return;
  throw new Error(
    `Circular track slot field '${obsoleteKey}' is obsolete. ` +
    'Use inner_gap_px and outer_gap_px for physical gaps.'
  );
};

/**
 * @param {TrackSlotRow} slot
 * @returns {TrackSlotRow} A copy without the obsolete keys; the argument itself when it is no object.
 */
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

/**
 * A legacy Circular slot spec under current tokens, and the legacy spacing
 * token it consumed: the last token of the first of `spacingKeys` the spec
 * holds, as stored, with its value and whether it applied to a gap the spec
 * leaves unset. `gapOf` gives the gap value a spacing fills, or null when it
 * fills none.
 * @param {unknown} spec
 * @param {readonly string[]} [spacingKeys] the keys read as the spacing, in lower case, the first one first
 * @param {(value: string) => string | null} [gapOf]
 * @returns {{ spec: string, spacing: { stored: string, value: string, filled: boolean } | null }}
 */
export const readLegacyCircularTrackSlotSpec = (spec, spacingKeys = ['spacing'], gapOf = (value) => value) => {
  const text = String(spec || '').trim();
  const atIndex = text.indexOf('@');
  if (atIndex < 0) return { spec: text, spacing: null };

  const head = text.slice(0, atIndex).trim();
  const tokens = text.slice(atIndex + 1).split(',').flatMap((token) => {
    const equalsIndex = token.indexOf('=');
    if (equalsIndex < 0) return [];
    const key = token.slice(0, equalsIndex).trim();
    return [{ stored: token.trim(), key, normalizedKey: key.toLowerCase(), value: token.slice(equalsIndex + 1).trim() }];
  });
  const spacingKey = spacingKeys.find((name) => tokens.some(({ normalizedKey }) => normalizedKey === name));
  const consumed = tokens.filter(({ normalizedKey }) => normalizedKey === spacingKey).at(-1);
  const retained = tokens.filter(({ normalizedKey }) => (
    normalizedKey !== spacingKey && !OBSOLETE_CIRCULAR_TRACK_SLOT_KEYS.has(normalizedKey)
  ));
  const filled = consumed && consumed.value !== ''
    ? ['inner_gap_px', 'outer_gap_px'].filter((name) => !retained.some(({ normalizedKey }) => normalizedKey === name))
    : [];
  const gap = consumed && filled.length > 0 ? gapOf(consumed.value) : null;
  const options = [
    ...retained.map(({ key, value }) => `${key}=${value}`),
    ...(gap === null ? [] : filled.map((name) => `${name}=${gap}`))
  ];
  return {
    spec: options.length > 0 ? `${head}@${options.join(',')}` : head,
    spacing: consumed ? { stored: consumed.stored, value: consumed.value, filled: filled.length > 0 } : null
  };
};

/** @param {unknown} spec */
export const migrateLegacyCircularTrackSlotSpec = (spec) => readLegacyCircularTrackSlotSpec(spec).spec;

/** @param {unknown} value */
export const normalizeSlotSide = (value) => normalizeOptionalPlacement(value);

/**
 * @param {SlotPlacementSource | null | undefined} slot
 * @param {unknown} [preset]
 */
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

/**
 * @param {SlotPlacementSource | null | undefined} slot
 * @param {unknown} [preset]
 */
const featureLaneForSlot = (slot, preset = 'tuckin') => (
  resolveCircularTrackFeaturePlacement(slot, preset).laneDirection
);

/**
 * @param {CircularTrackSlot | null | undefined} slot
 * @param {unknown} side
 */
const syncTickParamsForPlacement = (slot, side) => {
  if (!slot || slot.renderer !== 'ticks') return;
  void side;
  slot.params = cloneParams(slot.params);
  slot.params.tick_label_layout = normalizeTickLabelLayout(slot.params.tick_label_layout);
  delete slot.params.label_side;
  delete slot.params.tick_side;
};

/**
 * @param {CircularTrackSlot | null | undefined} slot Edited in place.
 * @param {unknown} side
 */
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

/**
 * @param {SlotPlacementSource | null | undefined} slot
 * @param {unknown} [preset]
 * @returns {string}
 */
export const effectiveSlotPlacement = (slot, preset = 'tuckin') => {
  if (!slot) return 'inside';
  if (slot.renderer === 'features') {
    return sideForLaneDirection(featureLaneForSlot(slot, preset));
  }
  return normalizePlacement(slot.side, 'inside');
};

/**
 * @param {(SlotPlacementSource | null | undefined)[]} slots
 * @param {unknown} [preset]
 * @returns {number}
 */
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

/**
 * @param {unknown} value
 * @param {unknown} slotCount
 * @param {FieldLeaf[] | null} [losses] Where a value it replaces is recorded (OV-758; `at` '').
 * @returns {number | null} null when value is no integer.
 */
export const clampCircularTrackAxisIndex = (value, slotCount, losses = null) => {
  const length = Math.max(0, Number(slotCount) || 0);
  const numeric = Number(value);
  const read = Number.isInteger(numeric) ? Math.min(Math.max(numeric, 0), length) : null;
  noteReplacedLeaf(losses, '', value, read);
  return read;
};

/**
 * @param {(CircularTrackSlot | null | undefined)[]} slots Edited in place.
 * @param {unknown} axisIndex
 * @param {unknown} [preset]
 * @returns {number} The Axis boundary used.
 */
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

/**
 * @param {(CircularTrackSlot | null | undefined)[]} slots Edited in place.
 * @param {unknown} axisIndex
 * @param {unknown} [preset]
 * @returns {number}
 */
export const enforceSingleOnAxisSlot = (slots, axisIndex, preset = 'tuckin') => {
  // The filter keeps exactly the integers of the map, which are the non-null indices.
  const onAxisIndices = /** @type {number[]} */ (slots
    .map((slot, index) => (
      slot?.renderer !== 'annotations' && effectiveSlotPlacement(slot, preset) === 'overlay' ? index : null
    ))
    .filter((index) => Number.isInteger(index)));
  if (onAxisIndices.length === 0) {
    return clampCircularTrackAxisIndex(axisIndex, slots.length) ?? inferLegacyAxisIndexFromFeature(slots, preset);
  }

  const clampedAxis = clampCircularTrackAxisIndex(axisIndex, slots.length);
  const keepIndex = clampedAxis !== null && onAxisIndices.includes(clampedAxis) ? clampedAxis : onAxisIndices[0];
  onAxisIndices.forEach((index) => {
    if (index === keepIndex) return;
    syncSlotPlacementFromSide(slots[index], index < keepIndex ? 'outside' : 'inside');
  });
  return keepIndex;
};

/**
 * @param {unknown} slots
 * @param {string} [defaultNt]
 * @param {string} [preset]
 * @param {number | null} [axisIndex]
 * @returns {CircularTrackSlot[]}
 */
export const applyCircularTrackOrderPlacements = (slots, defaultNt = 'GC', preset = 'tuckin', axisIndex = null) => {
  const normalized = normalizeCircularTrackSlots(slots, defaultNt, preset);
  const resolvedAxis = syncSlotsFromAxisIndex(normalized, axisIndex, preset);
  enforceSingleOnAxisSlot(normalized, resolvedAxis, preset);
  return normalized;
};

/**
 * @param {{ id: string, renderer: string, enabled?: boolean, width?: unknown, radius?: unknown, inner_gap_px?: unknown, outer_gap_px?: unknown, side?: string | null, z?: number, params?: unknown }} slot
 * @returns {CircularTrackSlot}
 */
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

/**
 * @param {Record<string, unknown>} params
 * @param {Record<string, unknown>} [expected]
 * @returns {boolean}
 */
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

/** @param {Record<string, unknown>} source */
const hasBlankSlotGeometry = (source) =>
  normalizeOptionalText(source.width) === null &&
  normalizeOptionalText(source.radius) === null &&
  normalizeOptionalText(source.inner_gap_px) === null &&
  normalizeOptionalText(source.outer_gap_px) === null;

/** @param {Record<string, unknown>} source */
const hasDefaultSlotFlags = (source) =>
  source.enabled !== false &&
  Number(source.z || 0) === 0;

/**
 * @param {TrackSlotRow} source
 * @param {string} renderer
 * @param {string} [defaultNt]
 * @param {string} [preset]
 * @returns {boolean}
 */
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
 * @returns {CircularTrackSlot[]}
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

/**
 * @param {Record<string, unknown>} [values] The geometry shortcuts of the form, by field.
 * @returns {Record<string, number | null>}
 */
export const normalizeCircularGeometryShortcuts = (values = {}) => {
  /** @type {Record<string, number | null>} */
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

/** @param {Record<string, unknown>} [values] */
export const hasCircularGeometryShortcuts = (values = {}) => (
  Object.values(normalizeCircularGeometryShortcuts(values))
    .some((value) => value !== null)
);

/**
 * @param {TrackSlotRow[]} slots
 * @param {Record<string, unknown>} [values]
 * @returns {Record<string, unknown>[]} Copies of the rows with the shortcut widths and radii.
 */
export const applyCircularGeometryShortcuts = (slots, values = {}) => {
  const geometry = normalizeCircularGeometryShortcuts(values);
  return (Array.isArray(slots) ? slots : []).map((slot) => {
    /** @type {Record<string, unknown>} */
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

/**
 * @param {TrackSlotRow} slot
 * @param {number} [index]
 * @param {string} [defaultNt]
 * @param {string} [preset]
 * @param {FieldLeaf[] | null} [losses] Where a value it replaces is recorded (OV-758;
 *   `at` under the stack, `[1].renderer`); a migration that keeps the meaning records nothing.
 * @returns {CircularTrackSlot}
 */
export const normalizeCircularTrackSlot = (slot, index = 0, defaultNt = 'GC', preset = 'tuckin', losses = null) => {
  /** @type {Record<string, unknown>} */
  const source = slot && typeof slot === 'object' && !Array.isArray(slot) ? slot : {};
  assertCurrentCircularTrackSlotShape(source);
  const resolved = circularRendererName(source.renderer);
  const renderer = SUPPORTED_RENDERERS.includes(resolved) ? resolved : 'dinucleotide_skew';
  if (renderer !== resolved) noteReplacedLeaf(losses, `[${index}].renderer`, source.renderer, renderer);
  const fallbackId = DEFAULT_SLOT_IDS[renderer] || `slot_${index + 1}`;
  const inheritsPresetDefaults = isLegacyDefaultWebSlotShape(source, renderer, defaultNt, preset);
  /** @type {Record<string, unknown>} */
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
    const storedLayout = params.tick_label_layout;
    params.tick_label_layout = normalizeTickLabelLayout(legacyLayout);
    noteReplacedLeaf(losses, `[${index}].params.tick_label_layout`, storedLayout, params.tick_label_layout);
    delete params.label_side;
    delete params.tick_side;
    if (normalizeOptionalText(params.preset) === null) {
      delete params.preset;
    } else {
      const storedPreset = params.preset;
      params.preset = normalizeCircularTrackPreset(params.preset);
      noteReplacedLeaf(losses, `[${index}].params.preset`, storedPreset, params.preset);
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
      noteReplacedLeaf(losses, `[${index}].params.${isUnsetValue(params.lane_direction) ? 'lanes' : 'lane_direction'}`,
        rawLaneDirection, placement.laneDirection);
      params.lane_direction = placement.laneDirection;
      side = placement.side;
    }
  }
  if (renderer === 'depth') {
    const storedTrackIndex = params.track_index;
    const parsedTrackIndex = normalizeTrackIndex(params.track_index);
    const idMatch = cleanToken(source.id, fallbackId).match(/^depth_(\d+)$/);
    if (parsedTrackIndex !== null) {
      params.track_index = parsedTrackIndex;
    } else if (idMatch) {
      params.track_index = Math.max(0, Number(idMatch[1]) - 1);
    } else {
      delete params.track_index;
    }
    if (parsedTrackIndex === null) noteReplacedLeaf(losses, `[${index}].params.track_index`, storedTrackIndex, params.track_index);
  }
  if (renderer === 'annotations') {
    normalizeAnnotationSlotParams(params, `[${index}].params`, losses);
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

/**
 * @param {unknown} slots
 * @param {string} [defaultNt]
 * @param {string} [preset]
 * @param {FieldLeaf[] | null} [losses] As `normalizeCircularTrackSlot`'s; a stack that is no list is recorded at ''.
 * @returns {CircularTrackSlot[]}
 */
export const normalizeCircularTrackSlots = (slots, defaultNt = 'GC', preset = 'tuckin', losses = null) => {
  const base = Array.isArray(slots)
    ? slots
    : createDefaultCircularTrackSlots({ nt: defaultNt, preset });
  if (base !== slots) noteReplacedLeaf(losses, '', slots, base);
  return base.map((slot, index) => normalizeCircularTrackSlot(slot, index, defaultNt, preset, losses));
};

/**
 * @param {unknown} spec `id:renderer@key=value,...`
 * @param {number} [index]
 * @param {string} [defaultNt]
 * @param {string} [preset]
 * @returns {CircularTrackSlot}
 */
export const parseCircularTrackSlotSpec = (spec, index = 0, defaultNt = 'GC', preset = 'tuckin') => {
  const text = String(spec || '').trim();
  const atIndex = text.indexOf('@');
  const head = (atIndex < 0 ? text : text.slice(0, atIndex)).trim();
  const separatorIndex = head.indexOf(':');
  /** @type {{ id: string, renderer: string, enabled: boolean, width?: string, radius?: string, inner_gap_px?: string | number | null, outer_gap_px?: string | number | null, side?: string, z?: number, params: Record<string, unknown> }} */
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

/**
 * @param {unknown} specs
 * @param {string} [defaultNt]
 * @param {string} [preset]
 * @returns {CircularTrackSlot[]}
 */
export const parseCircularTrackSlotSpecs = (specs, defaultNt = 'GC', preset = 'tuckin') => (
  Array.isArray(specs)
    ? specs.map((spec, index) => parseCircularTrackSlotSpec(spec, index, defaultNt, preset))
    : []
);

/**
 * @param {Pick<CircularTrackSlot, 'renderer' | 'params'>} slot
 * @returns {Record<string, unknown>}
 */
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
 * @param {TrackSlotRow} slot
 * @param {string} [defaultNt]
 * @param {string} [preset]
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
