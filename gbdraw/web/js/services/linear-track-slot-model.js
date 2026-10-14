// @ts-check
// The Linear track-slot model: renderer vocabulary, default stacks, slot
// normalization and schema migration, axis placement, and the request payload,
// as the request, Session and Gallery services read them. The Linear
// track-slot editor (app/linear-track-slots.js) edits it.
import { parseDepthTrackIndexIdentity } from './depth-track-state.js';
import { resolveColorToHex } from '../utils/color-utils.js';
import { normalizeOptionalText } from './track-slot-display.js';
import { requireCurrentLinearTrackLayout } from './current-option-values.js';
import { normalizeAnnotationSlotParams, parseOptionalPixel } from './track-slot-validation.js';
import { noteReplacedLeaf } from '../utils/field-reading.js';

/** @import { FieldLeaf } from '../utils/field-reading.js' */

/**
 * The parameters of a slot. The renderer decides which names are read; any other
 * name a Session slot carries is kept (Python owns the option set, R7).
 * @typedef {Record<string, unknown>} LinearTrackSlotParams
 */

/**
 * A Linear track slot row as `normalizeLinearTrackSlots` returns it (the Web
 * draft form: `height` and `spacing` are px text, or an invalid draft kept as typed).
 * @typedef {object} LinearTrackSlot
 * @property {string} id
 * @property {string} renderer
 * @property {boolean} enabled
 * @property {string} [depth_binding_error]
 * @property {string} side `above`, `below`, or `overlay`.
 * @property {unknown} height
 * @property {unknown} spacing
 * @property {number} z
 * @property {LinearTrackSlotParams} params
 */

/**
 * A slot row before normalization: a draft row, a Session slot, or a parsed spec.
 * @typedef {object} LinearTrackSlotInput
 * @property {unknown} [id]
 * @property {unknown} [renderer]
 * @property {unknown} [enabled]
 * @property {unknown} [depth_binding_error]
 * @property {unknown} [side]
 * @property {unknown} [height]
 * @property {unknown} [spacing]
 * @property {unknown} [z]
 * @property {unknown} [params]
 */

/** @typedef {{ renderer?: unknown, side?: unknown }} LinearSlotPlacementInput The fields that place a slot. */

/**
 * A canonical slot spec as `canonicalLinearTrackSlotSource` reads it.
 * @typedef {{
 *   id: string, renderer: string, enabled: boolean, side?: string | null, height?: string,
 *   spacing?: string, z?: unknown, params: LinearTrackSlotParams
 * }} LinearSlotSpec
 */

const SUPPORTED_RENDERERS = [
  'features',
  'dinucleotide_content',
  'dinucleotide_skew',
  'depth',
  'annotations',
  'spacer'
];

/** @type {Record<string, string>} */
export const RENDERER_ALIASES = {
  gc_content: 'dinucleotide_content',
  content: 'dinucleotide_content',
  gc_skew: 'dinucleotide_skew',
  skew: 'dinucleotide_skew'
};

/** @type {Record<string, string>} */
export const DEFAULT_SLOT_IDS = {
  features: 'features',
  dinucleotide_content: 'gc_content',
  dinucleotide_skew: 'gc_skew',
  depth: 'depth',
  annotations: 'annotations',
  spacer: 'spacer'
};

export const LINEAR_TRACK_SLOT_SCHEMA_VERSION = 2;
export const LEGACY_LINEAR_TRACK_SLOT_SCHEMA_VERSION = 1;

/**
 * @param {unknown} [params]
 * @returns {LinearTrackSlotParams}
 */
export const cloneParams = (params = {}) => {
  if (!params || typeof params !== 'object' || Array.isArray(params)) return {};
  return { ...params };
};

/** @param {unknown} value */
export const normalizeColorParam = (value) => {
  const text = normalizeOptionalText(value);
  if (text === null) return null;
  return resolveColorToHex(text);
};

/** @param {LinearTrackSlotParams} params */
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

/** @param {unknown} value @param {string} [fallback] */
export const normalizeRenderer = (value, fallback = 'features') => {
  const text = String(value || fallback).trim().toLowerCase();
  const renderer = RENDERER_ALIASES[text] || text;
  return SUPPORTED_RENDERERS.includes(renderer) ? renderer : fallback;
};

/** @param {unknown} value @param {string} [fallback] */
const normalizeSide = (value, fallback = 'below') => {
  const text = String(value || fallback).trim().toLowerCase();
  return ['above', 'below', 'overlay'].includes(text) ? text : fallback;
};

/** @param {unknown} value @param {string} [fallback] */
export const normalizePlacement = (value, fallback = 'below') => normalizeSide(value, fallback);

/** @param {unknown} [trackLayout] */
const sideForLinearTrackLayout = (trackLayout = 'middle') => {
  const normalized = requireCurrentLinearTrackLayout(trackLayout);
  if (normalized === 'above') return 'above';
  if (normalized === 'below') return 'below';
  return 'overlay';
};

/** @param {unknown} value @param {string} [fallback] */
export const normalizeNt = (value, fallback = 'GC') => {
  const text = String(value || '').trim().toUpperCase();
  return text || fallback;
};

/** @param {unknown} value */
export const normalizeTrackIndex = (value) => {
  if (value === null || value === undefined || value === '') return null;
  const numeric = Number(value);
  if (!Number.isInteger(numeric) || numeric < 0) return null;
  return numeric;
};

/** @param {unknown} value @param {string} field */
const normalizePxText = (value, field) => {
  try {
    const numeric = parseOptionalPixel(value, `Linear track slot ${field}`, { allowZero: field === 'spacing' });
    return numeric === null ? '' : `${numeric}px`;
  } catch {
    return value; // Keep an invalid draft intact for row feedback and submission failure.
  }
};

/**
 * @param {unknown} renderer
 * @param {LinearTrackSlotInput} [overrides]
 * @returns {LinearTrackSlot}
 */
export const defaultSlot = (renderer, overrides = {}) => {
  const normalizedRenderer = normalizeRenderer(renderer);
  const params = cloneParams(overrides.params);
  return {
    id: String(overrides.id || DEFAULT_SLOT_IDS[normalizedRenderer] || normalizedRenderer),
    renderer: normalizedRenderer,
    enabled: overrides.enabled !== false,
    side: normalizeSide(
      overrides.side,
      normalizedRenderer === 'features' ? 'overlay' : (normalizedRenderer === 'annotations' ? 'above' : 'below')
    ),
    height: normalizePxText(overrides.height, 'height'),
    spacing: normalizePxText(overrides.spacing, 'spacing'),
    z: Number.isInteger(Number(overrides.z)) ? Number(overrides.z) : 0,
    params
  };
};

/**
 * @typedef {object} DefaultLinearTrackSlotsOptions
 * @property {boolean} [showDepth]
 * @property {number} [depthTrackCount]
 * @property {boolean} [showGc]
 * @property {boolean} [showSkew]
 * @property {string} [nt] dinucleotide of the GC tracks
 * @property {unknown} [trackLayout] linear track layout (`above`, `middle`, `below`), checked when the side is read
 */

/**
 * @param {DefaultLinearTrackSlotsOptions} [options]
 * @returns {LinearTrackSlot[]}
 */
export const createDefaultLinearTrackSlots = ({
  showDepth = false,
  depthTrackCount = 1,
  showGc = false,
  showSkew = false,
  nt = 'GC',
  trackLayout = 'middle'
} = {}) => {
  const slots = [
    defaultSlot('features', {
      id: 'features',
      side: sideForLinearTrackLayout(trackLayout)
    })
  ];
  if (showDepth) {
    const count = Math.max(1, Number(depthTrackCount) || 1);
    for (let index = 0; index < count; index += 1) {
      slots.push(defaultSlot('depth', {
        id: count === 1 ? 'depth' : `depth_${index + 1}`,
        side: 'below',
        params: { track_index: index }
      }));
    }
  }
  if (showGc) {
    slots.push(defaultSlot('dinucleotide_content', {
      id: 'gc_content',
      side: 'below',
      params: { nt: normalizeNt(nt) }
    }));
  }
  if (showSkew) {
    slots.push(defaultSlot('dinucleotide_skew', {
      id: 'gc_skew',
      side: 'below',
      params: { nt: normalizeNt(nt) }
    }));
  }
  return slots;
};

/**
 * @param {unknown} value
 * @param {unknown} slotCount
 * @param {FieldLeaf[] | null} [losses] Where a value it replaces is recorded (OV-758; `at` '').
 */
export const clampLinearTrackAxisIndex = (value, slotCount, losses = null) => {
  if (value === null || value === undefined || value === '') return null;
  const numeric = Number(value);
  const read = Number.isInteger(numeric) ? Math.max(0, Math.min(Number(slotCount) || 0, numeric)) : null;
  noteReplacedLeaf(losses, '', value, read);
  return read;
};

/**
 * @param {(LinearTrackSlotInput | null | undefined)[] | null | undefined} slots
 * @param {string} [nt]
 * @param {unknown} [trackLayout]
 * @param {FieldLeaf[] | null} [losses] Where a value it replaces is recorded (OV-758; `at` under
 *   the stack, `[1].side`); a migration that keeps the meaning records nothing, and an empty stack holds nothing.
 * @returns {LinearTrackSlot[]}
 */
export const normalizeLinearTrackSlots = (slots, nt = 'GC', trackLayout = 'middle', losses = null) => {
  /** @type {(LinearTrackSlotInput | null | undefined)[]} */
  const source = Array.isArray(slots) && slots.length > 0
    ? slots
    : createDefaultLinearTrackSlots({ showGc: true, showSkew: true, nt, trackLayout });
  if (!Array.isArray(slots)) noteReplacedLeaf(losses, '', slots, source);
  /** @type {Set<string>} */
  const usedIds = new Set();
  let hasFeature = false;
  // The filter keeps the slot objects, so every row of the result is one.
  return /** @type {LinearTrackSlotInput[]} */ (source
    .filter((slot) => slot && typeof slot === 'object' && !Array.isArray(slot)))
    .map((slot, index) => {
      const renderer = normalizeRenderer(slot.renderer);
      if (!normalizeRenderer(slot.renderer, '')) noteReplacedLeaf(losses, `[${index}].renderer`, slot.renderer, renderer);
      const params = cloneParams(slot.params);
      if (renderer === 'features') {
        hasFeature = true;
      }
      if (renderer === 'depth') {
        const storedTrackIndex = params.track_index;
        const trackIndex = normalizeTrackIndex(storedTrackIndex);
        if (trackIndex === null) {
          if (slot.enabled === false) {
            delete params.track_index;
          } else {
            params.track_index = 0;
          }
          noteReplacedLeaf(losses, `[${index}].params.track_index`, storedTrackIndex, params.track_index);
        } else {
          params.track_index = trackIndex;
        }
      }
      if (renderer === 'dinucleotide_content' || renderer === 'dinucleotide_skew') {
        params.nt = normalizeNt(params.nt ?? params.dinucleotide, nt);
        delete params.dinucleotide;
        if (renderer === 'dinucleotide_skew') {
          normalizeSkewColorParams(params);
        }
      }
      if (renderer === 'annotations') normalizeAnnotationSlotParams(params, `[${index}].params`, losses);
      let id = String(slot.id || DEFAULT_SLOT_IDS[renderer] || `slot_${index + 1}`).trim();
      if (!id) id = `slot_${index + 1}`;
      if (usedIds.has(id)) id = `${id}_${index + 1}`;
      usedIds.add(id);
      const side = renderer === 'features'
        ? normalizeSide(slot.side, sideForLinearTrackLayout(trackLayout))
        : normalizeSide(slot.side, renderer === 'annotations' ? 'above' : 'below');
      const readSide = renderer !== 'features' && renderer !== 'annotations' && side === 'overlay' ? 'below' : side;
      noteReplacedLeaf(losses, `[${index}].side`, slot.side, readSide);
      return {
        id,
        renderer,
        enabled: slot.enabled !== false,
        ...(slot.depth_binding_error
          ? { depth_binding_error: String(slot.depth_binding_error) }
          : {}),
        side: readSide,
        height: normalizePxText(slot.height, 'height'),
        spacing: normalizePxText(slot.spacing, 'spacing'),
        z: Number.isInteger(Number(slot.z)) ? Number(slot.z) : 0,
        params
      };
    })
    .filter((slot) => slot.renderer !== 'features' || (hasFeature && slot.id));
};

/** @param {LinearSlotPlacementInput | null | undefined} slot */
export const effectiveLinearSlotPlacement = (slot) => {
  if (!slot) return 'below';
  const renderer = normalizeRenderer(slot.renderer);
  const side = normalizePlacement(slot.side, renderer === 'features' ? 'overlay' : 'below');
  if ((renderer === 'features' || renderer === 'annotations') && side === 'overlay') return 'overlay';
  return side === 'above' ? 'above' : 'below';
};

/** @param {(LinearSlotPlacementInput | null | undefined)[] | null | undefined} slots */
export const inferLinearTrackAxisIndexFromSlots = (slots) => {
  const normalizedSlots = Array.isArray(slots) ? slots : [];
  const featureIndex = normalizedSlots.findIndex((slot) => slot?.renderer === 'features');
  if (featureIndex >= 0) {
    const featurePlacement = effectiveLinearSlotPlacement(normalizedSlots[featureIndex]);
    if (featurePlacement === 'above') return featureIndex + 1;
    if (featurePlacement === 'overlay') return featureIndex;
    return featureIndex;
  }
  return normalizedSlots.filter((slot) => effectiveLinearSlotPlacement(slot) === 'above').length;
};

/**
 * @param {(LinearSlotPlacementInput | null | undefined)[] | null | undefined} slots
 * @param {number | null} [axisIndex]
 */
export const resolveLinearTrackAxisIndex = (slots, axisIndex = null) => {
  const normalizedSlots = Array.isArray(slots) ? slots : [];
  const clamped = clampLinearTrackAxisIndex(axisIndex, normalizedSlots.length);
  return clamped === null ? inferLinearTrackAxisIndexFromSlots(normalizedSlots) : clamped;
};

/** @param {LinearSlotPlacementInput | null | undefined} slot @param {unknown} placement */
export const syncLinearSlotPlacementFromSide = (slot, placement) => {
  if (!slot) return;
  const renderer = normalizeRenderer(slot.renderer);
  const normalizedPlacement = normalizePlacement(placement, renderer === 'features' ? 'overlay' : 'below');
  slot.side = renderer === 'features' || renderer === 'annotations' || normalizedPlacement !== 'overlay'
    ? normalizedPlacement
    : 'below';
};

/**
 * @param {(LinearSlotPlacementInput | null | undefined)[] | null | undefined} slots
 * @param {number | null | undefined} axisIndex
 */
export const syncLinearSlotsFromAxisIndex = (slots, axisIndex) => {
  const normalizedSlots = Array.isArray(slots) ? slots : [];
  const resolvedAxis = resolveLinearTrackAxisIndex(normalizedSlots, axisIndex);
  normalizedSlots.forEach((slot, index) => {
    if (!slot) return;
    if (normalizeRenderer(slot.renderer) === 'annotations' && effectiveLinearSlotPlacement(slot) === 'overlay') return;
    const isOnAxisFeature = (
      index === resolvedAxis &&
      normalizeRenderer(slot.renderer) === 'features' &&
      effectiveLinearSlotPlacement(slot) === 'overlay'
    );
    syncLinearSlotPlacementFromSide(
      slot,
      isOnAxisFeature ? 'overlay' : (index < resolvedAxis ? 'above' : 'below')
    );
  });
  return resolvedAxis;
};

/**
 * @param {(LinearSlotPlacementInput | null | undefined)[] | null | undefined} slots
 * @param {number | null | undefined} axisIndex
 */
export const enforceSingleLinearOnAxisSlot = (slots, axisIndex) => {
  const normalizedSlots = Array.isArray(slots) ? slots : [];
  const onAxisIndices = normalizedSlots
    .map((slot, index) => (
      normalizeRenderer(slot?.renderer) === 'features' &&
      effectiveLinearSlotPlacement(slot) === 'overlay'
        ? index
        : null
    ))
    .filter((index) => index !== null);
  if (onAxisIndices.length === 0) {
    return resolveLinearTrackAxisIndex(normalizedSlots, axisIndex);
  }

  const clampedAxis = clampLinearTrackAxisIndex(axisIndex, normalizedSlots.length);
  const keepIndex = clampedAxis !== null && onAxisIndices.includes(clampedAxis) ? clampedAxis : onAxisIndices[0];
  onAxisIndices.forEach((index) => {
    if (index === keepIndex) return;
    syncLinearSlotPlacementFromSide(normalizedSlots[index], index < keepIndex ? 'above' : 'below');
  });
  return keepIndex;
};

/** @param {unknown} value @param {string} fieldName @param {{ allowZero: boolean }} options */
const linearScalarPayload = (value, fieldName, { allowZero }) => {
  const numeric = parseOptionalPixel(value, fieldName, { allowZero });
  return numeric === null ? null : { value: numeric, unit: 'px' };
};

/** @param {LinearTrackSlotParams} params */
const canonicalLinearAnnotationParams = (params) => {
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

/** @param {{ params?: unknown, renderer?: unknown } | null | undefined} slot */
const canonicalLinearParams = (slot) => {
  const params = cloneParams(slot?.params);
  Object.keys(params).forEach((key) => {
    if (params[key] === null || params[key] === undefined || key.startsWith('_')) {
      delete params[key];
    }
  });
  if (slot?.renderer === 'annotations') return canonicalLinearAnnotationParams(params);
  return params;
};

/**
 * Encode one Web draft row as the canonical LinearTrackSlot object.
 * @param {LinearTrackSlotInput} slot
 */
export const buildLinearTrackSlotPayload = (slot) => {
  const normalized = normalizeLinearTrackSlots([slot])[0];
  if (!normalized) throw new Error('Linear track slot must be an object.');
  return {
    kind: 'linearTrackSlot',
    id: normalized.id,
    renderer: normalized.renderer,
    enabled: normalized.enabled,
    side: normalized.side,
    height: linearScalarPayload(
      normalized.height,
      `Linear track '${normalized.id}' height`,
      { allowZero: false }
    ),
    spacing: linearScalarPayload(
      normalized.spacing,
      `Linear track '${normalized.id}' spacing`,
      { allowZero: true }
    ),
    z: Number(normalized.z) || 0,
    params: canonicalLinearParams(normalized)
  };
};

/** @param {unknown} value */
const parseLinearTrackSlotRenderer = (value) => {
  const text = String(value ?? '').trim().toLowerCase();
  const renderer = RENDERER_ALIASES[text] || text;
  if (!SUPPORTED_RENDERERS.includes(renderer)) {
    throw new Error(`Unsupported linear track renderer: ${value}.`);
  }
  return renderer;
};

/** @param {unknown} value */
const parseLinearTrackSlotBoolean = (value) => {
  const text = String(value ?? '').trim().toLowerCase();
  if (['1', 'true', 'yes', 'on'].includes(text)) return true;
  if (['0', 'false', 'no', 'off'].includes(text)) return false;
  throw new Error(`Invalid linear track slot boolean: ${value}.`);
};

/** @param {string | { value?: unknown, unit?: unknown } | null | undefined} value @param {string} field */
const parseLinearTrackSlotPx = (value, field) => {
  if (value === null || value === undefined) return '';
  /** @type {unknown} */
  let rawValue = value;
  let rawUnit = 'px';
  if (value && typeof value === 'object' && !Array.isArray(value)) {
    const keys = Object.keys(value);
    if (keys.length !== 2 || !keys.includes('value') || !keys.includes('unit')) {
      throw new Error(`Canonical linear track slot ${field} must be a ScalarSpec.`);
    }
    rawValue = value.value;
    rawUnit = String(value.unit || '').trim().toLowerCase();
  }
  if (rawUnit !== 'px') {
    throw new Error(`Linear track slot ${field} only accepts px values.`);
  }
  const numeric = parseOptionalPixel(rawValue, `Linear track slot ${field}`, { allowZero: field === 'spacing' });
  return numeric === null ? '' : `${numeric}px`;
};

/** @param {unknown} value @param {string} field */
const parseStructuredLinearTrackSlotPx = (value, field) => {
  if (value === null) return '';
  if (!value || typeof value !== 'object' || Array.isArray(value)) {
    throw new Error(`Canonical linear track slot ${field} must be a ScalarSpec or null.`);
  }
  // `value` is an object that is not an array here.
  if (typeof /** @type {{ value?: unknown }} */ (value).value !== 'number') {
    throw new Error(`Canonical linear track slot ${field} ScalarSpec value must be numeric.`);
  }
  return parseLinearTrackSlotPx(/** @type {{ value?: unknown, unit?: unknown }} */ (value), field);
};

/**
 * @param {Record<string, unknown>} spec
 * @returns {LinearSlotSpec}
 */
const parseStructuredLinearTrackSlot = (spec) => {
  if (!spec || typeof spec !== 'object' || Array.isArray(spec) || spec.kind !== 'linearTrackSlot') {
    throw new Error('Canonical linear track slots must be strings or linearTrackSlot objects.');
  }
  const requiredKeys = ['kind', 'id', 'renderer', 'enabled', 'side', 'height', 'spacing', 'z', 'params'];
  const missingKey = requiredKeys.find((key) => !Object.prototype.hasOwnProperty.call(spec, key));
  if (missingKey) {
    throw new Error(`Canonical linear track slot is missing ${missingKey}.`);
  }
  const extraKey = Object.keys(spec).find((key) => !requiredKeys.includes(key));
  if (extraKey) {
    throw new Error(`Canonical linear track slot has unsupported field ${extraKey}.`);
  }
  if (typeof spec.id !== 'string' || typeof spec.renderer !== 'string') {
    throw new Error('Canonical linear track slot id and renderer must be strings.');
  }
  const id = spec.id.trim();
  if (!id) throw new Error('Canonical linear track slot id is required.');
  if (typeof spec.enabled !== 'boolean') {
    throw new Error('Canonical linear track slot enabled must be boolean.');
  }
  if (spec.side !== null && typeof spec.side !== 'string') {
    throw new Error('Canonical linear track slot side must be a string or null.');
  }
  const side = spec.side === null ? null : spec.side.trim().toLowerCase();
  if (side !== null && !['above', 'below', 'overlay'].includes(side)) {
    throw new Error(`Unsupported linear track slot side: ${spec.side}.`);
  }
  if (!Number.isInteger(spec.z)) {
    throw new Error('Canonical linear track slot z must be an integer.');
  }
  if (!spec.params || typeof spec.params !== 'object' || Array.isArray(spec.params)) {
    throw new Error('Canonical linear track slot params must be an object.');
  }
  const renderer = parseLinearTrackSlotRenderer(spec.renderer);
  const params = cloneParams(spec.params);
  if (
    spec.enabled &&
    renderer === 'depth' &&
    Object.prototype.hasOwnProperty.call(params, 'track_index')
  ) {
    params.track_index = parseDepthTrackIndexIdentity(
      params.track_index,
      `Depth slot '${id}' track_index`
    );
  }
  return {
    id,
    renderer,
    enabled: spec.enabled,
    side,
    height: parseStructuredLinearTrackSlotPx(spec.height, 'height'),
    spacing: parseStructuredLinearTrackSlotPx(spec.spacing, 'spacing'),
    z: spec.z,
    params
  };
};

/**
 * @param {unknown} spec
 * @returns {LinearSlotSpec}
 */
const parseStringLinearTrackSlot = (spec) => {
  if (typeof spec !== 'string') {
    throw new Error('Canonical linear track slots must be strings or linearTrackSlot objects.');
  }
  const text = spec.trim();
  if (!text) throw new Error('Canonical linear track slot cannot be empty.');
  const atIndex = text.indexOf('@');
  const head = (atIndex < 0 ? text : text.slice(0, atIndex)).trim();
  const separatorIndex = head.indexOf(':');
  if (separatorIndex < 0) {
    throw new Error(`Canonical linear track slot requires '<slot_id>:<renderer>': ${text}.`);
  }
  /** @type {LinearSlotSpec} */
  const source = {
    id: head.slice(0, separatorIndex).trim(),
    renderer: parseLinearTrackSlotRenderer(head.slice(separatorIndex + 1)),
    enabled: true,
    params: {}
  };
  if (!source.id) throw new Error('Canonical linear track slot id is required.');

  if (atIndex >= 0) {
    text.slice(atIndex + 1).split(',').forEach((token) => {
      if (!token.trim()) return;
      const equalsIndex = token.indexOf('=');
      if (equalsIndex < 0) throw new Error(`Invalid linear track slot option: ${token}.`);
      const key = token.slice(0, equalsIndex).trim().toLowerCase();
      const rawValue = token.slice(equalsIndex + 1).trim();
      if (!key) throw new Error(`Invalid linear track slot option: ${token}.`);
      if (key === 'id') source.id = rawValue;
      else if (key === 'renderer' || key === 'type') source.renderer = parseLinearTrackSlotRenderer(rawValue);
      else if (key === 'enabled' || key === 'show' || key === 'visible') {
        source.enabled = parseLinearTrackSlotBoolean(rawValue);
      } else if (key === 'h' || key === 'height') source.height = parseLinearTrackSlotPx(rawValue, 'height');
      else if (key === 'spacing') source.spacing = parseLinearTrackSlotPx(rawValue, 'spacing');
      else if (key === 'side') {
        const side = String(rawValue).trim().toLowerCase();
        if (!['above', 'below', 'overlay'].includes(side)) {
          throw new Error(`Unsupported linear track slot side: ${rawValue}.`);
        }
        source.side = side;
      } else if (key === 'z' || key === 'z_index' || key === 'zindex') {
        if (!rawValue) throw new Error(`Invalid linear track slot z: ${rawValue}.`);
        const z = Number(rawValue);
        if (!Number.isInteger(z)) throw new Error(`Invalid linear track slot z: ${rawValue}.`);
        source.z = z;
      } else if (key === 'nt' || key === 'dinucleotide') {
        source.params.nt = rawValue.toUpperCase();
      } else if (key === 'track_index') {
        source.params.track_index = parseDepthTrackIndexIdentity(
          rawValue,
          'Linear Depth slot track_index'
        );
      } else {
        source.params[key] = rawValue;
      }
    });
  }
  source.id = String(source.id || '').trim();
  if (!source.id) throw new Error('Canonical linear track slot id is required.');
  return source;
};

const GENERIC_LINEAR_TRACK_SLOT_FIELDS = new Set([
  'id', 'renderer', 'type', 'enabled', 'show', 'visible', 'side',
  'h', 'height', 'spacing', 'z', 'z_index', 'zindex'
]);

/** @param {LinearSlotSpec} source */
const validateCanonicalLinearTrackSlotSource = (source) => {
  if (
    source.side === 'overlay' &&
    source.renderer !== 'features' &&
    source.renderer !== 'annotations'
  ) {
    throw new Error('side=overlay is only supported for features and annotations slots.');
  }
  const genericParam = Object.keys(source.params || {}).find(
    (key) => GENERIC_LINEAR_TRACK_SLOT_FIELDS.has(String(key).trim().toLowerCase())
  );
  if (genericParam) {
    throw new Error(`Canonical linear track slot stores generic field ${genericParam} in params.`);
  }
  return source;
};

/** @param {unknown} spec */
const canonicalLinearTrackSlotSource = (spec) => validateCanonicalLinearTrackSlotSource(
  spec && typeof spec === 'object' && !Array.isArray(spec)
    ? parseStructuredLinearTrackSlot(/** @type {Record<string, unknown>} */ (spec)) // a non-array object
    : parseStringLinearTrackSlot(spec)
);

// `nt`: the dinucleotide of a content or skew slot that names none: the
// request's, as Python reads it (`slot_dinucleotide`; OV-490).
/** @param {unknown} spec @param {string} [nt] */
export const parseLinearTrackSlotSpec = (spec, nt = 'GC') => {
  const normalized = normalizeLinearTrackSlots([canonicalLinearTrackSlotSource(spec)], nt);
  return normalized[0] || null;
};

/** @param {unknown} specs @param {string} [nt] */
export const parseLinearTrackSlotSpecs = (specs, nt = 'GC') => {
  if (specs === null || specs === undefined) return [];
  if (!Array.isArray(specs)) throw new Error('Canonical linear track slots must be an array.');
  if (specs.length === 0) throw new Error('Canonical linear track slot list cannot be empty.');
  const parsed = specs.map((spec) => parseLinearTrackSlotSpec(spec, nt)).filter(Boolean);
  const ids = new Set();
  parsed.forEach((slot) => {
    if (ids.has(slot.id)) throw new Error(`Duplicate canonical linear track slot id: ${slot.id}.`);
    ids.add(slot.id);
  });
  if (parsed.filter((slot) => slot.enabled !== false && slot.renderer === 'features').length > 1) {
    throw new Error('Canonical linear track slots support only one enabled features slot.');
  }
  return parsed;
};

/**
 * @param {unknown} slots Anything but an array comes back as it is; an array of rows comes back as an array.
 * @param {number} [schemaVersion]
 * @returns {unknown}
 */
export const migrateLinearTrackSlotsToCurrentSchema = (
  slots,
  schemaVersion = LINEAR_TRACK_SLOT_SCHEMA_VERSION
) => {
  if (
    !Number.isInteger(schemaVersion) ||
    ![LEGACY_LINEAR_TRACK_SLOT_SCHEMA_VERSION, LINEAR_TRACK_SLOT_SCHEMA_VERSION].includes(schemaVersion)
  ) {
    throw new Error(`Unsupported linear track slot schema version: ${schemaVersion}.`);
  }
  const version = schemaVersion;
  if (!Array.isArray(slots)) return slots;
  return slots.map((slot) => {
    if (!slot || typeof slot !== 'object' || Array.isArray(slot)) return slot;
    const migrated = {
      ...slot,
      params: cloneParams(slot.params)
    };
    if (
      version === LEGACY_LINEAR_TRACK_SLOT_SCHEMA_VERSION &&
      normalizeRenderer(slot.renderer) === 'features'
    ) {
      delete migrated.height;
      delete migrated.spacing;
    }
    return migrated;
  });
};

/**
 * @param {(LinearTrackSlotInput | null | undefined)[] | null | undefined} slots
 * @param {number | null} [axisIndex]
 * @param {string} [nt]
 * @param {unknown} [trackLayout]
 */
export const applyLinearTrackOrderPlacements = (slots, axisIndex = null, nt = 'GC', trackLayout = 'middle') => {
  const normalized = normalizeLinearTrackSlots(slots, nt, trackLayout);
  const resolvedAxis = syncLinearSlotsFromAxisIndex(normalized, axisIndex);
  enforceSingleLinearOnAxisSlot(normalized, resolvedAxis);
  return normalized;
};

export { SUPPORTED_RENDERERS as LINEAR_TRACK_RENDERERS };
