// @ts-check
import { parseDepthTrackIndexIdentity } from './depth-track-state.js';
import { noteReplacedLeaf } from '../utils/field-reading.js';
import { diagnosticError } from '../utils/error-normalization.js';
import { msg, t } from '../utils/i18n.js';

/** @import { FieldLeaf } from '../utils/field-reading.js' */
/** @import { ErrorTextPart } from '../utils/error-normalization.js' */

/**
 * A `{ value, unit }` Width or Radius draft. Both fields are checked on read.
 * @typedef {object} CircularMeasureObject
 * @property {unknown} [value]
 * @property {unknown} [unit]
 */

/**
 * A row of a Custom Track Slots array as the plan validator reads it: keys of
 * unknown values. An element that is no object is reported as `slot_not_object`
 * and every key of it reads as undefined, as it does on a missing row.
 * @typedef {Record<string, unknown> | null | undefined} TrackSlotRow
 */

// An issue has two texts from one place (D-51): `message`, the English that
// Generate throws and the error normalizer reads, and `parts`, the sentence the
// slot editors show (`errorText`), one per issue code and worded for the row
// it sits in, without the mode or the row's name.
/**
 * @typedef {object} TrackSlotGlobalIssue A problem of the whole plan.
 * @property {string} code
 * @property {string} field
 * @property {string} message
 * @property {ErrorTextPart[]} parts
 * @property {'error'} severity
 */

/**
 * @typedef {object} TrackSlotRowIssue A problem of one row.
 * @property {string} code
 * @property {string} field
 * @property {string} message
 * @property {ErrorTextPart[]} parts
 * @property {number} rowIndex
 * @property {string | null} slotId
 * @property {'error'} severity
 */

/** @typedef {TrackSlotGlobalIssue | TrackSlotRowIssue} TrackSlotIssue */

/**
 * An enabled row the plan validator normalized: the row's own keys with these
 * overwritten.
 * @typedef {Record<string, unknown> & { id: string, renderer: string, enabled: boolean, side: string, z: number, params: Record<string, unknown> }} EnabledTrackSlot
 */

/**
 * @typedef {object} CustomTrackPlan The validated Custom Track draft.
 * @property {TrackSlotRow[]} draftSlots Detached copies of the rows.
 * @property {EnabledTrackSlot[]} enabledSlots
 * @property {number | null} emittedAxisIndex null when the Axis boundary is invalid.
 * @property {Map<number, TrackSlotRowIssue[]>} rowIssues By row index.
 * @property {TrackSlotGlobalIssue[]} globalIssues
 */

/**
 * @typedef {object} CustomTrackPlanInput What the plan validator reads of a draft.
 * @property {string} [mode] 'circular' or 'linear'; any other value is reported.
 * @property {unknown} [slots] The rows; a non-array is reported.
 * @property {unknown} [axisIndex] The Axis boundary; blank infers it.
 * @property {unknown} [trackType] The track preset or Linear layout.
 * @property {number | null} [depthTrackCount] null or undefined skips the Depth range check.
 * @property {unknown} [depthSourcedTrackIndexes] The Depth series that have a source; a non-array skips that check.
 * @property {unknown} [annotationSetIds] The annotation sets (ids or { id } entries); null or undefined skips the check.
 * @property {unknown} [conservationSeries] The Comparison Series entries.
 */

/**
 * @param {unknown} value
 * @returns {string}
 */
const normalizedString = (value) => String(value ?? '').trim().toLowerCase();

/**
 * @param {unknown} value
 * @returns {string | null}
 */
const normalizedOptionalString = (value) => {
  const text = normalizedString(value);
  return text || null;
};

/**
 * One shown sentence as the parts of an issue or fault.
 * @param {string} text A catalog key, marked with msg.
 * @param {Readonly<Record<string, string | number>>} [params]
 * @returns {ErrorTextPart[]}
 */
const sentence = (text, params) => [params ? { text, params } : { text }];
const CIRCULAR_GAPS_HINT = msg('Use inner_gap_px and outer_gap_px for physical gaps.');
/** The obsolete Circular field a row or its params uses. @param {string} field */
const obsoleteFieldParts = (field) => [
  ...sentence(msg("This row uses the obsolete field '{field}'."), { field }),
  ...sentence(CIRCULAR_GAPS_HINT)
];

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

/**
 * @param {unknown} value
 * @param {string} fieldName
 * @returns {number}
 */
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

// The decimal number text every track geometry field reads: Python's float
// grammar without its digit-group underscores. Number() also reads 0x10,
// 0b11 and 0o7, which the CLI rejects (TK-12).
const DECIMAL_NUMBER = '[+-]?(?:[0-9]+(?:\\.[0-9]*)?|\\.[0-9]+)(?:[eE][+-]?[0-9]+)?';
const DECIMAL_TEXT = new RegExp(`^${DECIMAL_NUMBER}$`);
// Pure physical-pixel text grammar; factor/% scalars keep their separate owner.
const PIXEL_TEXT = new RegExp(`^(${DECIMAL_NUMBER})(?:\\s*px)?$`, 'i');

/**
 * @param {unknown} value A number, or decimal number text.
 * @returns {number} NaN for any other value.
 */
const decimalValue = (value) => {
  if (typeof value === 'number') return value;
  const text = typeof value === 'string' ? value.trim() : '';
  return DECIMAL_TEXT.test(text) ? Number(text) : Number.NaN;
};

/**
 * @param {unknown} value A number, or text with an optional px suffix; null, undefined and blank text mean auto.
 * @param {string} fieldName The field the message names.
 * @param {{ allowZero: boolean }} options
 * @returns {number | null} null for auto.
 */
export const parseOptionalPixel = (value, fieldName, { allowZero }) => {
  if (value === null || value === undefined) return null;
  const invalid = () => new Error(
    `${fieldName} must be ${allowZero ? 'nonnegative' : 'positive'} finite number of pixels (px optional).`
  );
  let numeric;
  if (typeof value === 'number') {
    numeric = value;
  } else if (typeof value === 'string') {
    const text = value.trim();
    if (!text) return null;
    const match = PIXEL_TEXT.exec(text);
    if (!match) throw invalid();
    numeric = Number(match[1]);
  } else {
    throw invalid();
  }
  if (!Number.isFinite(numeric) || numeric < 0 || (!allowZero && numeric === 0)) {
    throw invalid();
  }
  return numeric === 0 ? 0 : numeric;
};

/**
 * @param {unknown} value A Width or Radius draft: text, a number, or { value, unit }.
 * @param {string} [fieldName] The field the message names.
 * @returns {{ value: number, unit: string } | null} null for auto.
 */
export const parseOptionalCircularScalar = (value, fieldName = 'Circular measure') => {
  if (value === null || value === undefined || value === '') return null;
  const invalid = () => new Error(`${fieldName} must be a number greater than 0, in px or ×R.`);
  if (value && typeof value === 'object' && !Array.isArray(value)) {
    // A non-array object: its value and unit are read as unknown and checked below.
    const measure = /** @type {CircularMeasureObject} */ (value);
    if (typeof measure.value !== 'number' && typeof measure.value !== 'string') throw invalid();
    const numeric = decimalValue(measure.value);
    const unit = normalizedString(measure.unit);
    if (!Number.isFinite(numeric) || numeric <= 0 || !['px', 'factor'].includes(unit)) {
      throw invalid();
    }
    return { value: numeric, unit };
  }
  if (typeof value === 'boolean' || (typeof value !== 'number' && typeof value !== 'string')) {
    throw invalid();
  }
  const text = String(value).trim();
  const isPx = /px$/i.test(text);
  const isPercent = /%$/.test(text);
  let numericText = text;
  if (isPx) numericText = text.slice(0, -2).trim();
  else if (isPercent) numericText = text.slice(0, -1).trim();
  const numeric = decimalValue(numericText);
  const resolved = isPercent ? numeric / 100 : numeric;
  if (!text || !Number.isFinite(resolved) || resolved <= 0) {
    throw invalid();
  }
  return { value: resolved, unit: isPx ? 'px' : 'factor' };
};

/** A Circular row's Width or Radius field. @typedef {'width' | 'radius'} CircularMeasureField */
const CIRCULAR_MEASURE_INVALID = Object.freeze({
  radius: msg('Radius must be a number greater than 0, in px or ×R.'),
  width: msg('Width must be a number greater than 0, in px or ×R.')
});
/**
 * The shown sentence of an invalid Circular Width or Radius: the row issue
 * and the row's own field show it (TK-12).
 * @param {CircularMeasureField} field
 * @returns {ErrorTextPart[]}
 */
export const circularMeasureInvalidParts = (field) => sentence(CIRCULAR_MEASURE_INVALID[field]);

const ANNOTATION_OVERFLOWS = ['error', 'compress', 'clip'];
/**
 * An annotation slot's params as both stacks read them, edited in place. What
 * the reading replaces is recorded in `losses` (OV-758).
 * @param {Record<string, unknown>} params
 * @param {string} at Where the params are under the stack (`[2].params`).
 * @param {FieldLeaf[] | null} [losses]
 */
export const normalizeAnnotationSlotParams = (params, at, losses = null) => {
  /** @type {[string, unknown][]} */
  const read = [
    ['set_id', String(params.set_id || '').trim()],
    ['overflow', ANNOTATION_OVERFLOWS.find((overflow) => overflow === String(params.overflow || '').toLowerCase()) ?? 'error'],
    ['show_labels', params.show_labels !== false && String(params.show_labels).toLowerCase() !== 'false'],
    ['layer', String(params.layer || '').toLowerCase() === 'underlay' ? 'underlay' : 'foreground']
  ];
  read.forEach(([field, value]) => {
    noteReplacedLeaf(losses, `${at}.${field}`, params[field], value);
    params[field] = value;
  });
  if (Array.isArray(params.marks)) {
    params.marks = Array.from(new Set(
      params.marks
        .map((mark) => String(mark || '').trim().toLowerCase())
        .filter(Boolean)
    ));
    // params.marks is the array assigned just above.
    if (/** @type {unknown[]} */ (params.marks).length === 0) delete params.marks;
  }
  for (const [field, defaultValue] of /** @type {[string, number][]} */ ([['lane_gap_px', 3], ['padding_px', 2]])) {
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
};

/**
 * A slot an Overlay annotation slot names by id (Python
 * `_OVERLAY_SLOT_REFERENCES`, gbdraw/tracks/parsing.py): `anchor_slot` is the
 * slot it is drawn on, `to_slot` (Extend to) the slot its span reaches.
 * @typedef {'anchor_slot' | 'to_slot'} OverlaySlotReference
 */
/**
 * The shown sentences of the row issues about one reference (D-51).
 * @typedef {object} OverlayReferenceText
 * @property {string} withoutOverlay
 * @property {string} unknown `{target}` is the slot id it names.
 * @property {string} ineligible `{target}` is the slot id it names.
 */
/**
 * Each reference with its name in the binding check's messages (Python's
 * wording), its issue codes (`annotation_<code>_unknown`, ...), its name in
 * the issue messages and the sentences the editors show for those issues.
 * @type {ReadonlyArray<{ key: OverlaySlotReference, role: string, code: string, label: string, text: OverlayReferenceText }>}
 */
const OVERLAY_SLOT_REFERENCES = Object.freeze([
  { key: 'anchor_slot', role: 'anchor', code: 'anchor', label: 'anchor', text: {
    withoutOverlay: msg('This row has an anchor but is not placed on it.'),
    unknown: msg("Anchor '{target}' is missing or disabled."),
    ineligible: msg("Anchor '{target}' has no drawable band.")
  } },
  { key: 'to_slot', role: 'to_slot', code: 'to_slot', label: 'Extend to slot', text: {
    withoutOverlay: msg('This row has an Extend to slot but is not placed on an anchor.'),
    unknown: msg("Extend to slot '{target}' is missing or disabled."),
    ineligible: msg("Extend to slot '{target}' has no drawable band.")
  } }
]);
/** The params only an Overlay annotation slot keeps: leaving Overlay drops them. */
export const OVERLAY_ANNOTATION_PARAMS = Object.freeze(['anchor_slot', 'to_slot', 'cover_anchor']);

/** A row's z as the plan check reads it: a draft that is no integer counts as 0. @param {unknown} value */
const rowZ = (value) => (Number.isSafeInteger(Number(value)) ? Number(value) : 0);

/**
 * Whether an Overlay annotation slot's z meets its layer against the z of a
 * slot it names (Python validate_overlay_annotation_anchors): an Underlay
 * draws below it, a Foreground above. Another layer is refused on its own.
 * @param {string} layer @param {number} z @param {number} targetZ
 */
const overlayLayerZHolds = (layer, z, targetZ) => (
  layer === 'underlay' ? z < targetZ : layer !== 'foreground' || z > targetZ
);

/**
 * The z an Overlay annotation row needs for its layer (OV-1297, LP-D91). The
 * GUI has no z control, so both track editors set it when Placement Overlay,
 * Layer, the anchor or Extend to changes: a z that meets the layer rule
 * against every slot the row names is kept (a loaded Session keeps its z),
 * else it is one below the lowest (Underlay) or one above the highest
 * (Foreground) of their z. A disabled slot counts too, so the row stays valid
 * when it is enabled again; the plan check refuses a disabled target anyway.
 * @param {{ z?: unknown, params?: Record<string, unknown> | null }} slot
 * @param {ReadonlyArray<{ id?: unknown, z?: unknown } | null | undefined>} slots The stack the row is in.
 * @returns {number}
 */
export const overlayAnnotationZ = (slot, slots) => {
  const params = slot.params || {};
  const layer = normalizedString(params.layer || 'foreground');
  const z = rowZ(slot.z);
  const targetZs = OVERLAY_SLOT_REFERENCES.flatMap(({ key }) => {
    const targetId = String(params[key] || '').trim();
    const target = targetId
      ? slots.find((candidate) => candidate && candidate !== slot && String(candidate.id || '').trim() === targetId)
      : null;
    return target ? [rowZ(target.z)] : [];
  });
  if (targetZs.every((targetZ) => overlayLayerZHolds(layer, z, targetZ))) return z;
  return layer === 'underlay' ? Math.min(...targetZs) - 1 : Math.max(...targetZs) + 1;
};

/**
 * An edit of an Overlay annotation row's Layer, anchor or Extend to, as both
 * track editors' selects make it: an empty value (Select anchor, None) unsets
 * the param, and the row takes the z its layer needs (`overlayAnnotationZ`).
 * @param {{ z: number, params: Record<string, unknown> }} slot
 * @param {Parameters<typeof overlayAnnotationZ>[1]} slots The stack the row is in.
 * @param {OverlaySlotReference | 'layer'} key
 * @param {unknown} value
 */
export const setOverlayAnnotationParam = (slot, slots, key, value) => {
  const text = String(value ?? '').trim();
  slot.params = { ...slot.params };
  if (text) slot.params[key] = text;
  else delete slot.params[key];
  slot.z = overlayAnnotationZ(slot, slots);
};

/**
 * The facts of a drawing that decide whether Add track can add a row.
 * @typedef {object} TrackAddFacts
 * @property {boolean} supported The mode draws this renderer.
 * @property {number} annotationSetCount Region annotation sets loaded in the drawing.
 * @property {number} depthSeriesCount Depth series with a loaded TSV.
 * @property {boolean} [featuresRowEnabled] The stack already has an enabled Features row (Linear draws one).
 */
/**
 * Why Add track cannot add a row of `renderer`, or '' when it can (OV-1254).
 * Both track editors read their availability and the button title here.
 * @param {string} renderer
 * @param {TrackAddFacts} facts
 * @returns {string}
 */
export const trackAddBlockReason = (renderer, { supported, annotationSetCount, depthSeriesCount, featuresRowEnabled = false }) => {
  if (!supported) return t('Choose a track type this diagram can draw');
  if (renderer === 'annotations' && annotationSetCount <= 0) return t('Import or add a Region annotation set in Region Annotations first');
  if (renderer === 'depth' && depthSeriesCount <= 0) return t('Load a Depth TSV first');
  if (renderer === 'features' && featuresRowEnabled) return t('The stack already has an enabled Features row; disable or delete it first');
  return '';
};

/**
 * An annotation slot's params as both modes write them to the request: a
 * default or an empty value is left out.
 * @param {Record<string, unknown>} params
 * @returns {Record<string, unknown>}
 */
export const canonicalAnnotationParams = (params) => {
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
  OVERLAY_SLOT_REFERENCES.forEach(({ key }) => {
    if (!String(next[key] || '').trim()) delete next[key];
  });
  return next;
};

/**
 * A leaf of a track slot row the binding check refuses, with the error it
 * raises: `field` is its place in the row (`width`, `params.side`).
 * @typedef {{ field: string, error: Error }} TrackSlotLeafFault
 */
/**
 * A fault the plan validator also reports as a row issue, with the sentence
 * the editors show for it (D-51).
 * @typedef {TrackSlotLeafFault & { parts: ErrorTextPart[] }} TrackSlotShownFault
 */

/**
 * @template {TrackSlotLeafFault} Fault
 * @param {Fault[]} faults Appended to.
 * @param {Omit<Fault, 'error'>} fault The fault the reader's error completes.
 * @param {() => unknown} read A reader that throws for a value it refuses.
 */
const checkLeaf = (faults, fault, read) => {
  try {
    read();
  } catch (error) {
    // `fault` and its error are a whole Fault.
    faults.push(/** @type {Fault} */ ({ ...fault, error: /** @type {Error} */ (error) }));
  }
};

/**
 * The geometry fields of a row the binding check refuses, in its order.
 * @param {Record<string, unknown>} slot The row.
 * @param {string} id
 * @param {string} layoutKind 'linear' or 'circular'; anything else is a fault.
 * @returns {TrackSlotShownFault[]}
 */
const slotGeometryFaults = (slot, id, layoutKind) => {
  /** @type {TrackSlotShownFault[]} */
  const faults = [];
  if (layoutKind === 'linear') {
    checkLeaf(faults, { field: 'height', parts: sentence(msg('Height must be a positive number of pixels (px optional).')) },
      () => parseOptionalPixel(slot.height, `Linear track slot '${id}' height`, { allowZero: false }));
    checkLeaf(faults, { field: 'spacing', parts: sentence(msg('Spacing must be a nonnegative number of pixels (px optional).')) },
      () => parseOptionalPixel(slot.spacing, `Linear track slot '${id}' spacing`, { allowZero: true }));
    return faults;
  }
  if (layoutKind !== 'circular') {
    return [{
      field: 'geometry',
      error: new Error(`Unsupported track-slot layout kind: ${layoutKind}.`),
      parts: sentence(msg('Unsupported track-slot layout kind: {kind}.'), { kind: layoutKind })
    }];
  }
  Object.keys(slot).filter((key) => OBSOLETE_CIRCULAR_TRACK_SLOT_KEYS.has(normalizedString(key))).forEach((key) => faults.push({
    field: key,
    error: new Error(`Circular track slot '${id}' uses obsolete field '${key}'. ${CIRCULAR_GAPS_HINT}`),
    parts: obsoleteFieldParts(key)
  }));
  for (const field of /** @type {const} */ (['radius', 'width'])) {
    checkLeaf(faults, { field, parts: circularMeasureInvalidParts(field) }, () => {
      try {
        parseOptionalCircularScalar(slot[field], `Circular track slot '${id}' ${field}`);
      } catch (error) {
        // The row's Width or Radius field shows this failure itself (TK-12).
        throw Object.assign(/** @type {Error} */ (error), { field });
      }
    });
  }
  for (const [field, text] of [
    ['inner_gap_px', msg('Inner gap must be a nonnegative number of pixels (px optional).')],
    ['outer_gap_px', msg('Outer gap must be a nonnegative number of pixels (px optional).')]
  ]) {
    checkLeaf(faults, { field, parts: sentence(text) },
      () => parseOptionalPixel(slot[field], `Circular track slot '${id}' ${field}`, { allowZero: true }));
  }
  return faults;
};

/**
 * The generic layout fields a row stores in its params, which the binding check refuses.
 * @param {Record<string, unknown>} params
 * @param {string} id
 * @param {string} layoutKind 'linear' selects the Linear keys; any other value the Circular keys.
 * @returns {TrackSlotShownFault[]}
 */
const genericParamFaults = (params, id, layoutKind) => {
  const genericKeys = layoutKind === 'linear' ? LINEAR_GENERIC_LAYOUT_KEYS : CIRCULAR_GENERIC_LAYOUT_KEYS;
  return Object.keys(params).filter((key) => genericKeys.has(normalizedString(key))).map((key) => (
    layoutKind === 'circular' && OBSOLETE_CIRCULAR_TRACK_SLOT_KEYS.has(normalizedString(key))
      ? {
          field: `params.${key}`,
          error: new Error(`Circular track slot '${id}' uses obsolete field 'params.${key}'. ${CIRCULAR_GAPS_HINT}`),
          parts: obsoleteFieldParts(`params.${key}`)
        }
      : {
          field: `params.${key}`,
          error: new Error(`${layoutKind === 'linear' ? 'Linear' : 'Circular'} track slot '${id}' stores generic layout field '${key}' in params.`),
          parts: sentence(msg("The layout field '{field}' belongs to the row, not to its params."), { field: key })
        }
  ));
};

/**
 * The leaves of one row (an object with an id) that the binding check
 * refuses, in the order it checks them. Generate and Save throw the first
 * (`validateTrackSlotBindingInvariants`); Load drops each and lists it
 * (`readImportedTrackSlots`).
 * @param {Record<string, unknown>} slot The row.
 * @param {string} id
 * @param {{ modeLabel: string, layoutKind: string, sideSet: ReadonlySet<string> }} options
 * @returns {TrackSlotLeafFault[]}
 */
export const trackSlotLeafFaults = (slot, id, { modeLabel, layoutKind, sideSet }) => {
  /** @type {TrackSlotLeafFault[]} */
  const faults = [];
  if (Object.hasOwn(slot, 'enabled') && typeof slot.enabled !== 'boolean') {
    faults.push({ field: 'enabled', error: new Error(`${modeLabel} track slot '${id}' enabled must be boolean.`) });
  }
  faults.push(...slotGeometryFaults(slot, id, layoutKind));
  if (slot.side !== null && slot.side !== undefined && slot.side !== '') {
    if (typeof slot.side !== 'string') {
      faults.push({ field: 'side', error: new Error(`${modeLabel} track slot '${id}' side must be a string or null.`) });
    } else if (!sideSet.has(normalizedString(slot.side))) {
      faults.push({ field: 'side', error: new Error(`${modeLabel} track slot '${id}' has unsupported side=${JSON.stringify(slot.side)}.`) });
    }
  }
  const params = Object.hasOwn(slot, 'params') ? slot.params : {};
  if (!params || typeof params !== 'object' || Array.isArray(params)) {
    faults.push({ field: 'params', error: new Error(`${modeLabel} track slot '${id}' params must be an object.`) });
  } else {
    faults.push(...genericParamFaults(/** @type {Record<string, unknown>} */ (params), id, layoutKind));
  }
  checkLeaf(faults, { field: 'z' }, () => parseSlotZIndex(slot.z, `${modeLabel} track slot '${id}' z`));
  return faults;
};

/**
 * @typedef {object} TrackSlotBindingOptions
 * @property {string} modeLabel 'Circular' or 'Linear', for the messages.
 * @property {string} layoutKind 'linear' or 'circular'.
 * @property {Iterable<string> | null | undefined} supportedRenderers
 * @property {Iterable<string> | null | undefined} supportedSides
 * @property {Iterable<string> | null | undefined} anchorlessRenderers Renderers an annotation cannot anchor to.
 * @property {number | null | undefined} [depthTrackCount] null or undefined skips the Depth range check.
 */

/**
 * A row the binding validator normalized.
 * @typedef {object} BoundTrackSlot
 * @property {string} id
 * @property {string} renderer
 * @property {boolean} enabled
 * @property {string | null} side
 * @property {number} z
 * @property {Record<string, unknown>} params
 */

/**
 * Validate structural and binding invariants shared by canonical slots and
 * persisted UI slots.
 *
 * The returned objects are detached normalized views used only for validation;
 * callers remain responsible for committing their original state after every
 * preflight check succeeds.
 * @param {unknown} slots Anything; a non-array throws.
 * @param {TrackSlotBindingOptions} options
 * @returns {BoundTrackSlot[]}
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
    const [leafFault] = trackSlotLeafFaults(slot, id, { modeLabel, layoutKind, sideSet });
    if (leafFault) throw leafFault.error;
    const enabled = slot.enabled !== false;
    /** @type {string | null} */
    let side = typeof slot.side === 'string' && slot.side !== '' ? normalizedString(slot.side) : null;
    // An object (`trackSlotLeafFaults`).
    const params = { .../** @type {Record<string, unknown>} */ (Object.hasOwn(slot, 'params') ? slot.params : {}) };
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
      if (renderer === 'features') {
        const lane = normalizedOptionalString(params.lane_direction ?? params.lanes);
        if (lane !== null && !['inside', 'outside', 'split'].includes(lane)) {
          throw new Error(
            `Circular features slot '${id}' has unsupported lane_direction=${JSON.stringify(lane)}.`
          );
        }
        side = side ?? (lane === 'split' ? 'overlay' : (lane ?? 'inside'));
        if (lane !== null) {
          const expectedSide = lane === 'split' ? 'overlay' : lane;
          if (side !== expectedSide) {
            throw new Error(
              `Circular features slot '${id}' has conflicting side=${JSON.stringify(side)} and lane_direction=${JSON.stringify(lane)}.`
            );
          }
        }
      } else {
        side = side ?? (renderer === 'annotations' ? 'outside' : 'inside');
      }
      if (renderer === 'sequence_conservation' && side === 'overlay') {
        throw new Error(`Circular sequence_conservation slot '${id}' cannot use side=overlay.`);
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
        throw diagnosticError('TRACK_INVALID', { reason: 'OVERLAY_ANCHOR', field: 'anchor_slot', slotIndex: index });
      }
      const offOverlay = side === 'overlay'
        ? undefined
        : OVERLAY_SLOT_REFERENCES.find(({ key }) => String(params[key] || '').trim());
      if (offOverlay) {
        throw diagnosticError('TRACK_INVALID', { reason: 'OVERLAY_SIDE', field: 'side', slotIndex: index });
      }
      const layer = normalizedOptionalString(params.layer) ?? 'foreground';
      if (!['foreground', 'underlay'].includes(layer)) {
        throw diagnosticError('TRACK_INVALID', { reason: 'LAYER', field: 'layer', slotIndex: index });
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
  normalized.forEach((slot, slotIndex) => {
    if (!slot.enabled || slot.renderer !== 'annotations' || slot.side !== 'overlay') return;
    const toSlot = String(slot.params.to_slot || '').trim();
    if (toSlot && toSlot === String(slot.params.anchor_slot || '').trim()) {
      throw diagnosticError('TRACK_INVALID', { reason: 'OVERLAY_SPAN', field: 'to_slot', slotIndex });
    }
    const layer = normalizedString(slot.params.layer || 'foreground');
    for (const { key } of OVERLAY_SLOT_REFERENCES) {
      const targetId = String(slot.params[key] || '').trim();
      if (key === 'to_slot' && !targetId) continue;
      const target = enabledById.get(targetId);
      if (!target) {
        throw diagnosticError('TRACK_INVALID', {
          reason: key === 'to_slot' ? 'OVERLAY_SPAN' : 'OVERLAY_ANCHOR', field: key, slotIndex
        });
      }
      if (target.renderer === 'annotations' || anchorlessSet.has(target.renderer)) {
        throw diagnosticError('TRACK_INVALID', {
          reason: key === 'to_slot' ? 'OVERLAY_SPAN' : 'OVERLAY_ANCHOR', field: key, slotIndex
        });
      }
      if (!overlayLayerZHolds(layer, slot.z, target.z)) {
        throw diagnosticError('TRACK_INVALID', { reason: 'DRAW_ORDER', field: 'z', slotIndex });
      }
    }
  });

  return normalized;
};

/** @type {Readonly<Record<string, ReadonlySet<string>>>} */
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

/** @type {Readonly<Record<string, ReadonlySet<string>>>} */
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
  ['to_slot', 'annotations'],
  ['layer', 'annotations'],
  ['cover_anchor', 'annotations'],
  ['style_override', 'annotations']
]);

/**
 * @param {TrackSlotRow} slot
 * @returns {TrackSlotRow}
 */
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

/**
 * What a check reports of a row issue.
 * @typedef {{ code: string, field: string, message: string, parts: ErrorTextPart[], slotId?: string | null }} RowIssueInput
 */

/**
 * @param {RowIssueInput & { rowIndex: number }} issue
 * @returns {TrackSlotRowIssue}
 */
const issueForRow = ({
  code,
  field,
  message,
  parts,
  rowIndex,
  slotId = null
}) => ({
  code,
  field,
  message,
  parts,
  rowIndex,
  slotId,
  severity: 'error'
});

/**
 * @param {{ code: string, field: string, message: string, parts: ErrorTextPart[] }} issue
 * @returns {TrackSlotGlobalIssue}
 */
const globalIssue = ({ code, field, message, parts }) => ({
  code,
  field,
  message,
  parts,
  severity: 'error'
});

/**
 * @param {Map<number, TrackSlotRowIssue[]>} rowIssues
 * @param {number} index
 * @param {RowIssueInput} issue
 * @returns {void}
 */
const addRowIssue = (rowIssues, index, issue) => {
  const current = rowIssues.get(index) || [];
  current.push(issueForRow({ ...issue, rowIndex: index }));
  rowIssues.set(index, current);
};

/** @param {string} mode */
const customTrackModeLabel = (mode) => (
  mode === 'circular' ? 'Circular' : 'Linear'
);

/**
 * @param {TrackSlotRow} slot
 * @param {string} mode
 * @param {unknown} trackType
 * @returns {string}
 */
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
  // A params value that is no object (a params_not_object issue) reads both keys as undefined.
  const lane = normalizedOptionalString(
    /** @type {Record<string, unknown> | undefined} */ (slot?.params)?.lane_direction ??
    /** @type {Record<string, unknown> | undefined} */ (slot?.params)?.lanes
  );
  if (lane === 'split') return 'overlay';
  if (lane === 'inside' || lane === 'outside') return lane;
  const preset = normalizedString(trackType);
  if (preset === 'middle') return 'overlay';
  if (preset === 'spreadout') return 'outside';
  return 'inside';
};

/**
 * @param {TrackSlotRow[]} slots
 * @param {string} mode
 * @param {unknown} trackType
 * @returns {number}
 */
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

/**
 * @param {unknown} value
 * @returns {number | null} NaN for an invalid value.
 */
const nonnegativeFiniteNumber = (value) => {
  if (value === null || value === undefined || value === '') return null;
  if (typeof value === 'boolean') return Number.NaN;
  const numeric = Number(value);
  return Number.isFinite(numeric) && numeric >= 0 ? numeric : Number.NaN;
};

/**
 * @param {string} renderer
 * @param {string} key
 * @returns {boolean}
 */
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

// The params that describe a row's track: `legend_label` (a Depth row's is its
// series title) and a comparison ring's `label`, `color` and `fileName`. Every
// renderer takes them (Python carries them, LP-D70), so they are in no
// PARAM_OWNER entry, and a renderer change drops them (OV-1186).
const TRACK_DESCRIPTION_PARAMS = new Set(['legend_label', 'label', 'color', 'fileName']);

/**
 * The params a row keeps when its renderer changes (TK-04): those the new
 * renderer accepts under PARAM_OWNER, without the params that describe the
 * previous track (`TRACK_DESCRIPTION_PARAMS`), so the new renderer starts
 * without them.
 * @param {unknown} params
 * @param {string} renderer
 * @returns {Record<string, unknown>}
 */
export const paramsKeptOnRendererChange = (params, renderer) => Object.fromEntries(
  Object.entries(params && typeof params === 'object' && !Array.isArray(params) ? params : {})
    .filter(([key]) => !TRACK_DESCRIPTION_PARAMS.has(key) && paramsAllowedForRenderer(renderer, normalizedString(key)))
);

/**
 * The ID of a track slot or annotation set that Add creates, in both modes: the
 * default ID of its kind (a slot's renderer, `annotations` for a set) while
 * nothing uses it, then that ID numbered `_2`, `_3`, and so on, taking the
 * first unused number. IDs compare after trimming.
 * @param {string} baseId
 * @param {ReadonlyArray<{ id?: unknown } | null | undefined>} items
 * @returns {string}
 */
export const nextNumberedId = (baseId, items) => {
  const usedIds = new Set(items.map((item) => String(item?.id || '').trim()).filter(Boolean));
  let id = baseId;
  for (let suffix = 2; usedIds.has(id); suffix += 1) id = `${baseId}_${suffix}`;
  return id;
};

/**
 * The ID of a track slot or annotation set that Duplicate creates, in both
 * modes: the source's own ID without one trailing `_<digits>`, numbered with
 * the first unused number from 2. A copy never takes the bare base, so copying
 * `my_skew` gives `my_skew_2`, and a copy of `my_skew_2` then gives `my_skew_3`,
 * not `my_skew_2_2`.
 * @param {string} sourceId
 * @param {ReadonlyArray<{ id?: unknown } | null | undefined>} items
 * @returns {string}
 */
export const nextCopyId = (sourceId, items) => {
  const trimmed = sourceId.trim();
  const baseId = trimmed.replace(/_\d+$/, '') || trimmed;
  return nextNumberedId(baseId, [...items, { id: baseId }]);
};

/**
 * @param {unknown} values
 * @returns {Set<string>}
 */
const normalizedInventoryIds = (values) => new Set(
  (Array.isArray(values) ? values : [])
    .map((value) => String(value?.id ?? value ?? '').trim())
    .filter(Boolean)
);

/**
 * @param {unknown} series
 * @returns {{ keys: Set<string>, indexes: Set<number>, indexesByKey: Map<string, Set<number>> }}
 */
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

/**
 * @param {CustomTrackPlanInput & { mode: string }} input
 * @returns {CustomTrackPlan}
 */
const collectCustomTrackIssues = ({
  mode,
  slots,
  axisIndex,
  trackType,
  depthTrackCount,
  depthSourcedTrackIndexes,
  annotationSetIds,
  conservationSeries
}) => {
  /** @type {Map<number, TrackSlotRowIssue[]>} */
  const rowIssues = new Map();
  /** @type {TrackSlotGlobalIssue[]} */
  const globalIssues = [];
  const draftSlots = Array.isArray(slots)
    ? slots.map(cloneValidationSlot)
    : [];
  const modeLabel = customTrackModeLabel(mode);

  if (!Array.isArray(slots)) {
    globalIssues.push(globalIssue({
      code: 'slots_not_array',
      field: 'slots',
      message: `${modeLabel} Custom Track Slots must be an array.`,
      parts: sentence(msg('Custom Track Slots must be a list of rows.'))
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
      message: `Unsupported Custom Track mode: ${String(mode)}.`,
      parts: sentence(msg('Unsupported Custom Track mode: {mode}.'), { mode: String(mode) })
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
  // A logical series inside the count may still have no source in any record.
  const sourcedDepthIndexes = Array.isArray(depthSourcedTrackIndexes)
    ? new Set(depthSourcedTrackIndexes)
    : null;
  if (!validDepthCount) {
    globalIssues.push(globalIssue({
      code: 'depth_count_invalid',
      field: 'depthTrackCount',
      message: `${modeLabel} Depth track count must be a non-negative integer.`,
      parts: sentence(msg('The Depth track count must be a non-negative integer.'))
    }));
  }

  /** @type {number[]} */
  const enabledIndexes = [];
  draftSlots.forEach((slot, index) => {
    if (!slot || typeof slot !== 'object' || Array.isArray(slot)) {
      addRowIssue(rowIssues, index, {
        code: 'slot_not_object',
        field: 'slot',
        message: `${modeLabel} track row ${index + 1} must be an object.`,
        parts: sentence(msg('This row must be an object.'))
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
        message: `${modeLabel} track row ${index + 1} Enabled must be true or false.`,
        parts: sentence(msg('Enabled must be true or false.'))
      });
    }
    if (slot.enabled !== false) enabledIndexes.push(index);
  });
  if (mode === 'linear' && enabledIndexes.length === 0) {
    globalIssues.push(globalIssue({
      code: 'linear_slots_empty',
      field: 'slots',
      message: 'Linear Custom Track Slots require at least one enabled row.',
      parts: sentence(msg('Custom Track Slots need at least one enabled row.'))
    }));
  }

  /** @type {number | null} */
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
      message: `${modeLabel} Axis boundary must be an integer from 0 to ${draftSlots.length}.`,
      parts: sentence(msg('The Axis boundary must be an integer from 0 to {max}.'), { max: draftSlots.length })
    }));
  }
  const emittedAxisIndex = resolvedAxisIndex === null
    ? null
    : draftSlots
        .slice(0, resolvedAxisIndex)
        .filter((slot) => slot && typeof slot === 'object' && slot.enabled !== false)
        .length;

  /** @type {Map<string, number[]>} */
  const idsToIndexes = new Map();
  enabledIndexes.forEach((index) => {
    const slot = draftSlots[index];
    if (!slot || typeof slot !== 'object' || Array.isArray(slot)) return;
    const id = typeof slot.id === 'string' ? slot.id.trim() : '';
    if (!id) {
      addRowIssue(rowIssues, index, {
        code: 'id_required',
        field: 'id',
        message: `${modeLabel} track row ${index + 1} requires a non-empty ID.`,
        parts: sentence(msg('This row needs a non-empty ID.'))
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
      message: `${modeLabel} track ID '${id}' is used by more than one enabled row.`,
      parts: sentence(msg("Track ID '{id}' is used by more than one enabled row."), { id })
    }));
  });

  /** @type {EnabledTrackSlot[]} */
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
        message: `${modeLabel} track ${slotLabel} has an unsupported renderer.`,
        parts: sentence(msg('This row has an unsupported renderer.'))
      });
    }

    /** @type {string | null} */
    let explicitSide = null;
    if (slot.side !== null && slot.side !== undefined && slot.side !== '') {
      explicitSide = typeof slot.side === 'string' ? normalizedString(slot.side) : null;
      if (!explicitSide || !supportedSides.has(explicitSide)) {
        addRowIssue(rowIssues, index, {
          code: 'side_unsupported',
          field: 'side',
          slotId: id || null,
          message: `${modeLabel} track ${slotLabel} has an unsupported side.`,
          parts: sentence(msg('This row has an unsupported side.'))
        });
        explicitSide = null;
      }
    }

    /** @type {Record<string, unknown> | null} */
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
        message: `${modeLabel} track ${slotLabel} Params must be an object.`,
        parts: sentence(msg('Params must be an object.'))
      });
    }

    const [geometryFault] = slotGeometryFaults(slot, id || `#${index + 1}`, mode);
    if (geometryFault) {
      // The cast: a Width or Radius error names its field (TK-12).
      const { field } = /** @type {Error & { field?: unknown }} */ (geometryFault.error);
      addRowIssue(rowIssues, index, {
        code: 'geometry_invalid',
        field: typeof field === 'string' ? field : 'geometry',
        slotId: id || null,
        message: geometryFault.error.message,
        parts: geometryFault.parts
      });
    }
    try {
      parseSlotZIndex(slot.z, `${modeLabel} track ${slotLabel} z`);
    } catch (error) {
      addRowIssue(rowIssues, index, {
        code: 'z_invalid',
        field: 'z',
        slotId: id || null,
        message: error instanceof Error ? error.message : String(error),
        parts: sentence(msg('The z value must be an integer.'))
      });
    }
    if (params) {
      const [genericFault] = genericParamFaults(params, id || `#${index + 1}`, mode);
      if (genericFault) {
        addRowIssue(rowIssues, index, {
          code: 'generic_param', field: 'params', slotId: id || null, message: genericFault.error.message, parts: genericFault.parts
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
          message: `${modeLabel} track ${slotLabel} retains '${mismatchedParam}' from another renderer.`,
          parts: sentence(msg("This row keeps '{param}' from another renderer."), { param: mismatchedParam })
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
          message: `Circular Features track ${slotLabel} has an invalid lane direction.`,
          parts: sentence(msg('This row has an invalid lane direction.'))
        });
      } else if (lane !== null) {
        const laneSide = lane === 'split' ? 'overlay' : lane;
        if (explicitSide !== null && explicitSide !== laneSide) {
          addRowIssue(rowIssues, index, {
            code: 'feature_lane_side_conflict',
            field: 'side',
            slotId: id || null,
            message: `Circular Features track ${slotLabel} side conflicts with its lane direction.`,
            parts: sentence(msg("This row's side conflicts with its lane direction."))
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
        message: `Linear track ${slotLabel} can use On Axis only for Features or Annotations.`,
        parts: sentence(msg('Only Features and Annotations rows can use On axis.'))
      });
    }
    if (mode === 'circular' && renderer === 'sequence_conservation' && side === 'overlay') {
      addRowIssue(rowIssues, index, {
        code: 'conservation_overlay',
        field: 'side',
        slotId: id || null,
        message: `Circular comparison track ${slotLabel} cannot use Overlay.`,
        parts: sentence(msg('A comparison row cannot use Overlay.'))
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
          message: `Circular Ticks track ${slotLabel} uses obsolete '${obsoleteTickKey}'.`,
          parts: sentence(msg("This row uses the obsolete parameter '{param}'."), { param: obsoleteTickKey })
        });
      }
    }

    if (renderer === 'depth' && params) {
      /** @type {number | null} */
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
          message: error instanceof Error ? error.message : String(error),
          parts: sentence(msg('The Depth track index must be a non-negative integer.'))
        });
      }
      if (trackIndex !== null && validDepthCount && depthCountProvided) {
        if (trackIndex >= depthTrackCount && depthTrackCount > 0) {
          addRowIssue(rowIssues, index, {
            code: 'depth_track_index_range',
            field: 'params.track_index',
            slotId: id || null,
            message: `${modeLabel} Depth track ${slotLabel} references index ${trackIndex}; available indexes are 0..${depthTrackCount - 1}.`,
            parts: sentence(msg('Depth track index {index} is not available; the available indexes are 0–{last}.'), { index: trackIndex, last: depthTrackCount - 1 })
          });
        } else if (depthTrackCount === 0 || (sourcedDepthIndexes && !sourcedDepthIndexes.has(trackIndex))) {
          addRowIssue(rowIssues, index, {
            code: 'depth_source_missing',
            field: 'params.track_index',
            slotId: id || null,
            message: `${modeLabel} Depth track ${slotLabel} has no logical Depth source.`,
            parts: sentence(msg('This row has no Depth source.'))
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
          message: `${modeLabel} Annotation track ${slotLabel} requires an annotation set.`,
          parts: sentence(msg('This row needs an annotation set.'))
        });
      } else if (annotationInventory && !annotationInventory.has(setId)) {
        addRowIssue(rowIssues, index, {
          code: 'annotation_set_unknown',
          field: 'params.set_id',
          slotId: id || null,
          message: `${modeLabel} Annotation track ${slotLabel} references unknown set '${setId}'.`,
          parts: sentence(msg("Annotation set '{set}' does not exist."), { set: setId })
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
          message: `${modeLabel} Annotation track ${slotLabel} has invalid mark filters.`,
          parts: sentence(msg('This row has invalid mark filters.'))
        });
      }
      [
        ['lane_gap_px', 'annotation_lane_gap', 'lane gap', msg('Lane gap must be a nonnegative finite number.')],
        ['padding_px', 'annotation_padding', 'padding', msg('Padding must be a nonnegative finite number.')]
      ].forEach(([field, code, label, text]) => {
        if (!Object.prototype.hasOwnProperty.call(params, field)) return;
        if (Number.isNaN(nonnegativeFiniteNumber(params[field]))) {
          addRowIssue(rowIssues, index, {
            code,
            field: `params.${field}`,
            slotId: id || null,
            message: `${modeLabel} Annotation track ${slotLabel} ${label} must be a nonnegative finite number.`,
            parts: sentence(text)
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
          message: `${modeLabel} Annotation track ${slotLabel} Cover anchor must be true or false.`,
          parts: sentence(msg('Cover anchor must be true or false.'))
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
          message: `${modeLabel} Annotation track ${slotLabel} has invalid overflow behavior.`,
          parts: sentence(msg('This row has an invalid overflow behavior.'))
        });
      }
      const layer = normalizedOptionalString(params.layer) ?? 'foreground';
      if (!['foreground', 'underlay'].includes(layer)) {
        addRowIssue(rowIssues, index, {
          code: 'annotation_layer',
          field: 'params.layer',
          slotId: id || null,
          message: `${modeLabel} Annotation track ${slotLabel} has an invalid layer.`,
          parts: sentence(msg('This row has an invalid layer.'))
        });
      }
      const anchorId = String(params.anchor_slot || '').trim();
      if (side === 'overlay' && !anchorId) {
        addRowIssue(rowIssues, index, {
          code: 'annotation_anchor_required',
          field: 'params.anchor_slot',
          slotId: id || null,
          message: `${modeLabel} Annotation track ${slotLabel} requires an Overlay anchor.`,
          parts: sentence(msg('This row needs an anchor.'))
        });
      } else if (side !== 'overlay') {
        OVERLAY_SLOT_REFERENCES.filter(({ key }) => String(params[key] || '').trim()).forEach(({ key, code, label, text }) => {
          addRowIssue(rowIssues, index, {
            code: `annotation_${code}_without_overlay`,
            field: `params.${key}`,
            slotId: id || null,
            message: `${modeLabel} Annotation track ${slotLabel} has an ${label} but is not Overlay.`,
            parts: sentence(text.withoutOverlay)
          });
        });
      }
    }

    if (renderer === 'sequence_conservation' && params) {
      if (String(params.managed || '').trim() !== CONSERVATION_MANAGER) {
        addRowIssue(rowIssues, index, {
          code: 'conservation_unmanaged',
          field: 'params.managed',
          slotId: id || null,
          message: `Circular comparison track ${slotLabel} is not owned by Comparison Series.`,
          parts: sentence(msg('This comparison row is not managed by Series.'))
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
          message: `Circular comparison track ${slotLabel} is not bound to an active Comparison Series source.`,
          parts: sentence(msg('This comparison row is not bound to an active Series source.'))
        });
      }
    }

    const normalizedSlot = {
      ...slot,
      id,
      renderer,
      enabled: true,
      side,
      z: rowZ(slot.z),
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
      message: 'Linear Custom Track Slots support only one enabled Features row.',
      parts: sentence(msg('Only one enabled Features row is allowed.'))
    }));
  }
  // Whether the plan needs a Features row for underlays depends on the
  // record's features, which only Python reads; it reports TRACK_INVALID
  // FEATURES_COUNT, as the CLI does (TK-13).

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
    // The checks of Python's validate_overlay_annotation_anchors, in its
    // order; like Python, a row reports its first fault.
    if (String(slot.params.to_slot || '').trim() === anchorId) {
      addRowIssue(rowIssues, rowIndex, {
        code: 'annotation_to_slot_same',
        field: 'params.to_slot',
        slotId: slot.id || null,
        message: `${modeLabel} Annotation track '${slot.id}' extends to its own anchor '${anchorId}'; choose another slot or None.`,
        parts: sentence(msg("This row extends to its own anchor '{anchor}'; choose another slot or None."), { anchor: anchorId })
      });
      return;
    }
    const layer = normalizedString(slot.params.layer || 'foreground');
    for (const { key, code, label, text } of OVERLAY_SLOT_REFERENCES) {
      const targetId = String(slot.params[key] || '').trim();
      if (!targetId) continue;
      const target = enabledById.get(targetId);
      /** @type {RowIssueInput | null} */
      let issue = null;
      if (!target) {
        issue = {
          code: `annotation_${code}_unknown`,
          field: `params.${key}`,
          slotId: slot.id || null,
          message: `${modeLabel} Annotation track '${slot.id}' references missing or disabled ${label} '${targetId}'.`,
          parts: sentence(text.unknown, { target: targetId })
        };
      } else if (INELIGIBLE_ANCHOR_RENDERERS.has(target.renderer)) {
        issue = {
          code: `annotation_${code}_ineligible`,
          field: `params.${key}`,
          slotId: slot.id || null,
          message: `${modeLabel} Annotation track '${slot.id}' ${label} '${targetId}' has no eligible drawable band.`,
          parts: sentence(text.ineligible, { target: targetId })
        };
      } else if (!overlayLayerZHolds(layer, slot.z, target.z)) {
        issue = {
          code: `annotation_${layer}_z`,
          field: 'z',
          slotId: slot.id || null,
          message: `${modeLabel} Annotation ${layer} '${slot.id}' must have z ${layer === 'underlay' ? 'below' : 'above'} ${label} '${targetId}'.`,
          parts: sentence(layer === 'underlay'
            ? msg("As an underlay, this row must be drawn below '{target}'.")
            : msg("As a foreground annotation, this row must be drawn above '{target}'."), { target: targetId })
        };
      }
      if (issue) {
        addRowIssue(rowIssues, rowIndex, issue);
        return;
      }
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
          featureLane !== null && ['inside', 'outside'].includes(featureLane) &&
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
        message: `${modeLabel} track '${slot.id}' side conflicts with the enabled Axis boundary.`,
        parts: sentence(msg("This row's side conflicts with the Axis boundary."))
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
 * @param {CustomTrackPlanInput} [options]
 * @returns {CustomTrackPlan}
 */
export const validateCustomTrackPlan = (options = {}) => collectCustomTrackIssues({
  mode: normalizedString(options.mode),
  slots: options.slots,
  axisIndex: options.axisIndex,
  trackType: options.trackType,
  depthTrackCount: options.depthTrackCount,
  depthSourcedTrackIndexes: options.depthSourcedTrackIndexes,
  annotationSetIds: options.annotationSetIds,
  conservationSeries: options.conservationSeries
});

export class CustomTrackPlanValidationError extends Error {
  /** @param {TrackSlotIssue[]} issues */
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

/**
 * @param {CustomTrackPlan | null | undefined} plan
 * @returns {TrackSlotIssue[]}
 */
export const customTrackPlanIssues = (plan) => {
  if (!plan || typeof plan !== 'object') return [];
  const rowIssues = plan.rowIssues instanceof Map
    ? Array.from(plan.rowIssues.values()).flat()
    : [];
  return [...rowIssues, ...(Array.isArray(plan.globalIssues) ? plan.globalIssues : [])];
};

/**
 * @param {CustomTrackPlan} plan
 * @returns {CustomTrackPlan} The plan itself; throws when it has an issue.
 */
export const assertValidCustomTrackPlan = (plan) => {
  const issues = customTrackPlanIssues(plan);
  if (issues.length > 0) throw new CustomTrackPlanValidationError(issues);
  return plan;
};
