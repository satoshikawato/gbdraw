// @ts-check
// The Linear track-slot model: renderer vocabulary, default stacks, slot
// normalization and schema migration, axis placement, and the request payload,
// as the request, Session and Gallery services read them. The Linear
// track-slot editor (app/linear-track-slots.js) edits it.
import { parseDepthTrackIndexIdentity } from './depth-track-state.js';
import { resolveColorToHex } from '../utils/color-utils.js';
import { normalizeOptionalText } from './track-slot-display.js';
import { requireCurrentLinearTrackLayout } from './current-option-values.js';
import { parseOptionalPixel } from './track-slot-validation.js';

const SUPPORTED_RENDERERS = [
  'features',
  'dinucleotide_content',
  'dinucleotide_skew',
  'depth',
  'annotations',
  'spacer'
];

export const RENDERER_ALIASES = {
  gc_content: 'dinucleotide_content',
  content: 'dinucleotide_content',
  gc_skew: 'dinucleotide_skew',
  skew: 'dinucleotide_skew'
};

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

export const cloneParams = (params = {}) => {
  if (!params || typeof params !== 'object' || Array.isArray(params)) return {};
  return { ...params };
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

export const normalizeRenderer = (value, fallback = 'features') => {
  const text = String(value || fallback).trim().toLowerCase();
  const renderer = RENDERER_ALIASES[text] || text;
  return SUPPORTED_RENDERERS.includes(renderer) ? renderer : fallback;
};

const normalizeSide = (value, fallback = 'below') => {
  const text = String(value || fallback).trim().toLowerCase();
  return ['above', 'below', 'overlay'].includes(text) ? text : fallback;
};

export const normalizePlacement = (value, fallback = 'below') => normalizeSide(value, fallback);

const sideForLinearTrackLayout = (trackLayout = 'middle') => {
  const normalized = requireCurrentLinearTrackLayout(trackLayout);
  if (normalized === 'above') return 'above';
  if (normalized === 'below') return 'below';
  return 'overlay';
};

export const normalizeNt = (value, fallback = 'GC') => {
  const text = String(value || '').trim().toUpperCase();
  return text || fallback;
};

export const normalizeTrackIndex = (value) => {
  if (value === null || value === undefined || value === '') return null;
  const numeric = Number(value);
  if (!Number.isInteger(numeric) || numeric < 0) return null;
  return numeric;
};

const normalizePxText = (value, field) => {
  try {
    const numeric = parseOptionalPixel(value, `Linear track slot ${field}`, { allowZero: field === 'spacing' });
    return numeric === null ? '' : `${numeric}px`;
  } catch {
    return value; // Keep an invalid draft intact for row feedback and submission failure.
  }
};

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
 * @property {string} [trackLayout] linear track layout (`above`, `middle`, `below`)
 */

/**
 * @param {DefaultLinearTrackSlotsOptions} [options]
 * @returns {Record<string, any>[]}
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

export const clampLinearTrackAxisIndex = (value, slotCount) => {
  if (value === null || value === undefined || value === '') return null;
  const numeric = Number(value);
  if (!Number.isInteger(numeric)) return null;
  return Math.max(0, Math.min(Number(slotCount) || 0, numeric));
};

export const normalizeLinearTrackSlots = (slots, nt = 'GC', trackLayout = 'middle') => {
  const source = Array.isArray(slots) && slots.length > 0
    ? slots
    : createDefaultLinearTrackSlots({ showGc: true, showSkew: true, nt, trackLayout });
  const usedIds = new Set();
  let hasFeature = false;
  return source
    .filter((slot) => slot && typeof slot === 'object' && !Array.isArray(slot))
    .map((slot, index) => {
      const renderer = normalizeRenderer(slot.renderer);
      const params = cloneParams(slot.params);
      if (renderer === 'features') {
        hasFeature = true;
      }
      if (renderer === 'depth') {
        const trackIndex = normalizeTrackIndex(params.track_index);
        if (trackIndex === null) {
          if (slot.enabled === false) {
            delete params.track_index;
          } else {
            params.track_index = 0;
          }
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
      }
      let id = String(slot.id || DEFAULT_SLOT_IDS[renderer] || `slot_${index + 1}`).trim();
      if (!id) id = `slot_${index + 1}`;
      if (usedIds.has(id)) id = `${id}_${index + 1}`;
      usedIds.add(id);
      const side = renderer === 'features'
        ? normalizeSide(slot.side, sideForLinearTrackLayout(trackLayout))
        : normalizeSide(slot.side, renderer === 'annotations' ? 'above' : 'below');
      return {
        id,
        renderer,
        enabled: slot.enabled !== false,
        ...(slot.depth_binding_error
          ? { depth_binding_error: String(slot.depth_binding_error) }
          : {}),
        side: renderer !== 'features' && renderer !== 'annotations' && side === 'overlay' ? 'below' : side,
        height: normalizePxText(slot.height, 'height'),
        spacing: normalizePxText(slot.spacing, 'spacing'),
        z: Number.isInteger(Number(slot.z)) ? Number(slot.z) : 0,
        params
      };
    })
    .filter((slot) => slot.renderer !== 'features' || (hasFeature && slot.id));
};

export const effectiveLinearSlotPlacement = (slot) => {
  if (!slot) return 'below';
  const renderer = normalizeRenderer(slot.renderer);
  const side = normalizePlacement(slot.side, renderer === 'features' ? 'overlay' : 'below');
  if ((renderer === 'features' || renderer === 'annotations') && side === 'overlay') return 'overlay';
  return side === 'above' ? 'above' : 'below';
};

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

/** @param {number | null} [axisIndex] */
export const resolveLinearTrackAxisIndex = (slots, axisIndex = null) => {
  const normalizedSlots = Array.isArray(slots) ? slots : [];
  const clamped = clampLinearTrackAxisIndex(axisIndex, normalizedSlots.length);
  return clamped === null ? inferLinearTrackAxisIndexFromSlots(normalizedSlots) : clamped;
};

export const syncLinearSlotPlacementFromSide = (slot, placement) => {
  if (!slot) return;
  const renderer = normalizeRenderer(slot.renderer);
  const normalizedPlacement = normalizePlacement(placement, renderer === 'features' ? 'overlay' : 'below');
  slot.side = renderer === 'features' || renderer === 'annotations' || normalizedPlacement !== 'overlay'
    ? normalizedPlacement
    : 'below';
};

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

const linearScalarPayload = (value, fieldName, { allowZero }) => {
  const numeric = parseOptionalPixel(value, fieldName, { allowZero });
  return numeric === null ? null : { value: numeric, unit: 'px' };
};

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

const parseLinearTrackSlotRenderer = (value) => {
  const text = String(value ?? '').trim().toLowerCase();
  const renderer = RENDERER_ALIASES[text] || text;
  if (!SUPPORTED_RENDERERS.includes(renderer)) {
    throw new Error(`Unsupported linear track renderer: ${value}.`);
  }
  return renderer;
};

const parseLinearTrackSlotBoolean = (value) => {
  const text = String(value ?? '').trim().toLowerCase();
  if (['1', 'true', 'yes', 'on'].includes(text)) return true;
  if (['0', 'false', 'no', 'off'].includes(text)) return false;
  throw new Error(`Invalid linear track slot boolean: ${value}.`);
};

const parseLinearTrackSlotPx = (value, field) => {
  if (value === null || value === undefined) return '';
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

const parseStructuredLinearTrackSlotPx = (value, field) => {
  if (value === null) return '';
  if (!value || typeof value !== 'object' || Array.isArray(value)) {
    throw new Error(`Canonical linear track slot ${field} must be a ScalarSpec or null.`);
  }
  if (typeof value.value !== 'number') {
    throw new Error(`Canonical linear track slot ${field} ScalarSpec value must be numeric.`);
  }
  return parseLinearTrackSlotPx(value, field);
};

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

const canonicalLinearTrackSlotSource = (spec) => validateCanonicalLinearTrackSlotSource(
  spec && typeof spec === 'object' && !Array.isArray(spec)
    ? parseStructuredLinearTrackSlot(spec)
    : parseStringLinearTrackSlot(spec)
);

export const parseLinearTrackSlotSpec = (spec) => {
  const normalized = normalizeLinearTrackSlots([canonicalLinearTrackSlotSource(spec)]);
  return normalized[0] || null;
};

export const parseLinearTrackSlotSpecs = (specs) => {
  if (specs === null || specs === undefined) return [];
  if (!Array.isArray(specs)) throw new Error('Canonical linear track slots must be an array.');
  if (specs.length === 0) throw new Error('Canonical linear track slot list cannot be empty.');
  const parsed = specs.map((spec) => parseLinearTrackSlotSpec(spec)).filter(Boolean);
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

/** @param {number | null} [axisIndex] */
export const applyLinearTrackOrderPlacements = (slots, axisIndex = null, nt = 'GC', trackLayout = 'middle') => {
  const normalized = normalizeLinearTrackSlots(slots, nt, trackLayout);
  const resolvedAxis = syncLinearSlotsFromAxisIndex(normalized, axisIndex);
  enforceSingleLinearOnAxisSlot(normalized, resolvedAxis);
  return normalized;
};

export { SUPPORTED_RENDERERS as LINEAR_TRACK_RENDERERS };
