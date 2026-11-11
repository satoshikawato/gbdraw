// @ts-check
// The record display model: source-bound record display rows, the requested
// display start and topology, and the record orientation intent, as the
// request and Session services read them. The Record Display controls
// (app/record-display-options.js) edit it.
import { regionWithinRecord, resolveDisambiguatedRecordSelection } from './record-options.js';
import { diagnosticError } from '../utils/error-normalization.js';
import { isPlainObject } from '../utils/contract-shape.js';
import { cloneJsonData } from './json-clone.js';
import { msg, t, tm } from '../utils/i18n.js';

/** @import { RecordEntryInput } from './record-options.js' */
/** @import { FieldReading } from '../utils/field-reading.js' */
/** @import { DrawingId } from './drawing-view.js' */
/** @typedef {'circular' | 'linear'} RecordDisplayScope */
// The topology a file declares, one whole key each (`detectedTopology` of services/record-discovery.js).
const DETECTED_TOPOLOGY_NOTES = Object.freeze({
  circular: msg('(file: circular)'), linear: msg('(file: linear)'), unknown: msg('(file: topology unknown)')
});
/** @param {unknown} topology a record row's `detectedTopology` */
export const detectedTopologyNote = (topology) => tm(topology === 'circular' || topology === 'linear'
  ? DETECTED_TOPOLOGY_NOTES[topology] : DETECTED_TOPOLOGY_NOTES.unknown);
/**
 * @typedef {{
 *   anchor: 'five-prime' | 'midpoint' | 'three-prime' | null,
 *   biologicalFeatureId: string,
 *   offsetBp: number,
 *   orientForward: boolean,
 *   placement: 'anchor' | 'feature-end',
 *   recordKey: string,
 *   schema: 1
 * }} AnchorIntent The feature a record's display is anchored to.
 */
/**
 * The display a record row asks for; a draft holds every field.
 * @typedef {{
 *   topologyOverride?: boolean | null,
 *   startCoordinate?: number | null,
 *   reverseComplementOverride?: boolean | null
 * }} RecordDisplayRequested
 */
/**
 * The crop a record row draws, as its Display start reads it: the span of the first window the row draws (1-based,
 * source coordinates; a Start after the End runs across the origin of a circular record) and that window's
 * orientation (`recordDisplayCrop`).
 * @typedef {{ start: number, end: number, reverseComplement: boolean }} RecordDisplayCrop
 */
/**
 * What the row's source says about the draft's record: the row's cropped and
 * reverse state, and the crop it draws (null while it draws none, or while the
 * crop cannot be read yet).
 * @typedef {{ cropped?: unknown, reverse?: unknown, crop?: RecordDisplayCrop | null }} RecordDisplayContext
 */
/**
 * A discovered record row as the surface reads it.
 * @typedef {RecordDisplayContext & {
 *   recordId?: unknown, recordLength?: number | null, detectedTopology?: unknown
 * }} RecordDisplayRowView
 */
/**
 * The identity of a record row or draft: its source input and record selector.
 * @typedef {object} RecordDisplayDraftIdentity
 * @property {string} sourceUid `circular`, or a Linear File card's UID
 * @property {string} selector
 */
/**
 * A draft of the display of one record (a drawing holds those of its own mode).
 * @typedef {RecordDisplayDraftIdentity & {
 *   recordId?: unknown, topologyOverride: boolean | null, startCoordinate: number | null,
 *   reverseComplementOverride: boolean | null, anchorIntent: unknown
 * }} RecordDisplayDraft
 */
/**
 * A record row of one source input, as `buildRecordDisplayRows` builds it with the input's crop and orientation.
 * @typedef {RecordDisplayDraftIdentity & RecordDisplayRowView & {
 *   key: string, scope: string, source: unknown, paired?: unknown
 * }} RecordDisplayRow
 */

/** @param {unknown} selector `#N`, as text (OV-639: a list is no selector). */
const isRecordSelectorText = (selector) => typeof selector === 'string' && /^#[1-9]\d*$/.test(selector);

// A record display row or draft names its source input and record. A drawing
// holds the drafts of its own inputs (PD-OI-090): `circular`, or a
// Linear File card's UID.
/** @param {{ sourceUid?: string, selector: string }} identity */
export const recordDisplayKey = ({ sourceUid, selector }) => {
  if (!sourceUid || !isRecordSelectorText(selector)) {
    throw new Error('Record display requires a source UID and exact record selector.');
  }
  return JSON.stringify([sourceUid, selector]);
};

/**
 * @template {RecordEntryInput} T
 * @param {{
 *   scope: RecordDisplayScope, sourceUid: string, source: unknown,
 *   records: readonly T[] | null | undefined, selector?: unknown
 * }} input
 */
export const buildRecordDisplayRows = ({ scope, sourceUid, source, records, selector = '' }) => {
  const selection = resolveDisambiguatedRecordSelection(records, selector);
  if (!['unspecified', 'resolved'].includes(selection.status)) {
    throw new Error(`Record selector is ${selection.status}: ${selector}.`);
  }
  return (selection.record ? [selection.record] : selection.entries).map((record) => {
    const row = { ...record, scope, sourceUid, source };
    return { ...row, key: recordDisplayKey(row) };
  });
};

/**
 * @template {{ sourceUid: string, selector: string }} D
 * @param {readonly D[]} drafts
 * @param {readonly { sourceUid?: string, selector: string }[]} discoveredRows
 * @param {readonly string[]} [replacedSourceUids]
 * @returns {D[]}
 */
export const reconcileRecordDisplayDrafts = (drafts, discoveredRows, replacedSourceUids = []) => {
  // Pass all discovered rows, including inactive selectors. A source replacement
  // purges its drafts even when the existing source card keeps its UID.
  const keys = new Set(discoveredRows.map(recordDisplayKey));
  return drafts.filter((draft) => !replacedSourceUids.includes(draft.sourceUid)
    && keys.has(recordDisplayKey(draft)));
};

// OV-278: a draft projected from a request that selects its record by ID (a
// CLI, Python or Gallery Session, which has no Web draft) has the selector
// '#1' whatever the record's place, because the projection does not read the
// file. Once the source's rows are known, it belongs to the row of its source
// with its record ID; a request selects only an ID that is unique in its file
// (`buildDisambiguatedRecordEntries`). The drafts are the drawing `drawingId`'s,
// so only that drawing's rows bind them: a drawing Find in other files creates
// has cards of its source drawing's uids. Returns `drafts` when no draft moves.
/**
 * @template {Array<RecordDisplayDraftIdentity & { recordId?: unknown }>} Drafts the record ID is only compared
 * @param {Drafts} drafts
 * @param {readonly { key: string, drawingId: DrawingId, sourceUid: string, selector: string, recordId: string }[]} rows
 * @param {DrawingId} drawingId The drawing whose drafts they are.
 * @returns {Drafts}
 */
export const bindRecordDisplayDrafts = (drafts, rows, drawingId) => {
  const taken = new Set(drafts.map(recordDisplayKey));
  let moved = false;
  const bound = drafts.map((draft) => {
    const sourceRows = rows.filter((row) => row.drawingId === drawingId && row.sourceUid === draft.sourceUid);
    if (!draft.recordId || sourceRows.some((row) => row.selector === draft.selector
      && row.recordId === draft.recordId)) return draft;
    const named = sourceRows.filter((row) => row.recordId === draft.recordId);
    if (named.length !== 1 || taken.has(named[0].key)) return draft;
    taken.add(named[0].key);
    moved = true;
    return { ...draft, selector: named[0].selector };
  });
  return moved ? /** @type {Drafts} */ (bound) : drafts;
};

/**
 * @param {unknown} value
 * @returns {number | null}
 */
export const parseRecordDisplayStart = (value) => {
  if (value === null || (typeof value === 'string' && !value.trim())) return null;
  if (!['string', 'number'].includes(typeof value)) throw new Error('Display start must be an integer.');
  const number = Number(value);
  if (!Number.isSafeInteger(number) || number < 1) {
    throw new Error('Display start must be a positive integer.');
  }
  return number;
};

const ANCHOR_INTENT_KEYS = [
  'anchor',
  'biologicalFeatureId',
  'offsetBp',
  'orientForward',
  'placement',
  'recordKey',
  'schema'
];

/** @param {unknown} text Text without a NUL. */
const isIdText = (text) => typeof text === 'string' && Boolean(text) && !text.includes('\0');
/** @param {unknown} value An anchor intent (`validateAnchorIntent`) or not. */
const isAnchorIntent = (value) => {
  if (!value || typeof value !== 'object' || Array.isArray(value)) return false;
  const intent = /** @type {Record<string, unknown>} */ (value);
  return Object.keys(intent).sort().join(',') === ANCHOR_INTENT_KEYS.join(',')
    && intent.schema === 1 && isIdText(intent.recordKey) && isIdText(intent.biologicalFeatureId)
    && (/** @type {readonly unknown[]} */ (['anchor', 'feature-end'])).includes(intent.placement)
    && (intent.placement === 'anchor'
      ? (/** @type {readonly unknown[]} */ (['five-prime', 'midpoint', 'three-prime'])).includes(intent.anchor)
      : intent.anchor === null)
    && Number.isSafeInteger(intent.offsetBp)
    && typeof intent.orientForward === 'boolean';
};
/**
 * The value is untrusted; every field read below is checked before it is trusted.
 * @param {Partial<Record<keyof AnchorIntent, unknown>> | null | undefined} intent
 * @returns {AnchorIntent | null}
 */
export const validateAnchorIntent = (intent) => {
  if (intent === null) return null;
  if (!isAnchorIntent(intent)) throw new Error('Invalid record display anchor intent.');
  return /** @type {AnchorIntent} */ (intent); // every field was checked above
};

/**
 * @param {unknown} drafts
 * @returns {unknown}
 */
export const migrateLegacyRecordDisplayDrafts = (drafts) => {
  if (!Array.isArray(drafts)) return drafts;
  return drafts.map((draft) => {
    if (!draft || Object.hasOwn(draft, 'reverseComplementOverride')
      || Object.hasOwn(draft, 'anchorIntent')) return draft;
    return { ...draft, reverseComplementOverride: null, anchorIntent: null };
  });
};

// A Session 41-44 draft also named each row's mode (`scoped`); the Session 46
// split moves each row into its mode's drawing. A row names its source input
// and a record by a `#N` selector, as text (OV-575, OV-639).
const DRAFT_FIELDS = 'anchorIntent,recordId,reverseComplementOverride,selector,sourceUid,startCoordinate,topologyOverride';
const SCOPED_DRAFT_FIELDS = 'anchorIntent,recordId,reverseComplementOverride,scope,selector,sourceUid,startCoordinate,topologyOverride';
/** @param {unknown} value @param {boolean} scoped */
const isRecordDisplayDraft = (value, scoped) => {
  if (!value || typeof value !== 'object' || Array.isArray(value)) return false;
  const draft = /** @type {Record<string, unknown>} */ (value);
  const { startCoordinate } = draft;
  return Object.keys(draft).sort().join(',') === (scoped ? SCOPED_DRAFT_FIELDS : DRAFT_FIELDS)
    && (!scoped || draft.scope === 'circular' || draft.scope === 'linear')
    && typeof draft.recordId === 'string' && !draft.recordId.includes('\0')
    && isIdText(draft.sourceUid) && isRecordSelectorText(draft.selector)
    && (draft.topologyOverride === null || typeof draft.topologyOverride === 'boolean')
    && (draft.reverseComplementOverride === null || typeof draft.reverseComplementOverride === 'boolean')
    && (startCoordinate === null || (Number.isSafeInteger(startCoordinate) && /** @type {number} */ (startCoordinate) >= 1))
    && (draft.anchorIntent === null || isAnchorIntent(draft.anchorIntent));
};

/**
 * The record display drafts Load reads, row by row (OV-639, OV-611; E-OD04): a
 * row that is no draft, or that repeats an earlier row's record, is dropped;
 * a value that is no list is read as none.
 * @param {unknown} drafts
 * @param {{ scoped?: boolean }} [options]
 * @returns {FieldReading<RecordDisplayDraft[]>}
 */
export const readRecordDisplayDrafts = (drafts, { scoped = false } = {}) => {
  if (!Array.isArray(drafts)) return { value: [], lost: drafts !== undefined && drafts !== null };
  /** @type {Set<string>} */
  const identities = new Set();
  const kept = /** @type {RecordDisplayDraft[]} */ (drafts.filter((draft) => {
    if (!isRecordDisplayDraft(draft, scoped)) return false;
    const key = JSON.stringify([scoped ? draft.scope : '', draft.sourceUid, draft.selector]);
    if (identities.has(key)) return false;
    identities.add(key);
    return true;
  }));
  const dropped = drafts.length - kept.length;
  return { value: kept, lost: dropped > 0, ...(dropped ? { dropped } : {}) };
};

/**
 * Save keeps only drafts Load reads without a loss.
 * @param {unknown} drafts
 * @param {{ scoped?: boolean }} [options]
 */
export const validateRecordDisplayDrafts = (drafts, { scoped = false } = {}) => {
  if (!Array.isArray(drafts) || readRecordDisplayDrafts(drafts, { scoped }).lost) {
    throw diagnosticError('INPUT_INVALID', { field: 'schema', reason: 'FIELDS' });
  }
  return drafts;
};

/**
 * The crop a window draws, as a record row's Display start reads it; null while either end is not stated.
 * @param {{ start: number | null, end: number | null, reverseComplement: boolean }} window
 * @returns {RecordDisplayCrop | null}
 */
export const recordDisplayCrop = ({ start, end, reverseComplement }) => (
  start === null || end === null ? null : { start, end, reverseComplement }
);

/**
 * The coordinate a crop's drawing starts at (Circular 12 o'clock, the Linear left end): its End when it is drawn
 * reverse-complemented, else its Start. An End past the end of a record of known length is drawn to the record end
 * (`regionWithinRecord`); a window across the origin is drawn as stated.
 * @param {RecordDisplayCrop} crop
 * @param {number | null | undefined} recordLength
 */
const cropDisplayStart = ({ start, end, reverseComplement }, recordLength) => {
  if (!reverseComplement) return start;
  return start <= end ? regionWithinRecord(start, end, recordLength)?.drawnEnd ?? end : end;
};

/**
 * @param {RecordDisplayRowView} row
 * @param {RecordDisplayRequested} [draft] the row's requested display
 * @param {RecordDisplayContext} [context]
 */
export const recordDisplaySurface = (row, draft = {}, { cropped = false, reverse = false, crop = null } = {}) => {
  const override = draft.topologyOverride ?? null;
  if (override !== null && typeof override !== 'boolean') throw new Error('Topology override must be boolean or null.');
  const effectiveCircular = override ?? row.detectedTopology === 'circular';
  const lengthKnown = Number.isSafeInteger(row.recordLength) && /** @type {number} */ (row.recordLength) > 0;
  const startEnabled = effectiveCircular && lengthKnown && !cropped;
  const cropStart = crop === null ? null : cropDisplayStart(crop, row.recordLength);
  const currentStart = cropStart ?? (reverse ? (lengthKnown ? row.recordLength : null) : 1);
  return {
    effectiveCircular,
    startEnabled,
    disabledReason: cropped ? t('Display start is unavailable for a cropped record.')
      : !lengthKnown ? t('Record length is unavailable.')
        : !effectiveCircular ? t('Display start requires a circular record.') : '',
    currentStart,
    // OV-374: the Display start field shows the requested start, kept while
    // the field is disabled, or else the current one; never an empty field.
    // A crop disables it and the drawing starts where the crop does (OV-3341).
    displayStart: cropStart ?? draft.startCoordinate ?? currentStart
  };
};

/**
 * @param {RecordDisplayRowView} row
 * @param {RecordDisplayRequested} [draft]
 * @param {RecordDisplayContext} [context]
 */
export const requestedRecordDisplay = (row, draft = {}, context = {}) => {
  const surface = recordDisplaySurface(row, draft, context);
  const start = parseRecordDisplayStart(draft.startCoordinate ?? null);
  if (surface.startEnabled && start !== null && start > /** @type {number} */ (row.recordLength)) { // startEnabled holds only for a known length
    throw new Error(`Display start for ${row.recordId} must be between 1 and ${row.recordLength}.`);
  }
  return { isCircular: draft.topologyOverride ?? null,
    startCoordinate: surface.startEnabled ? start : null };
};

/**
 * The windows a stored schema <= 9 record's `region` reads as at schema 10: one window, none for a value that is no
 * object (a missing field is null; only integers can be reversed). A region with the larger number first is the reversed span (v0.14.5); only a schema 10 window whose start
 * is after its end crosses the origin, so the window states the span in order (Python's `_decode_region`). The one
 * reader of a stored region's span and orientation, for the request promotion (`promoteRecordRegionsToWindows`,
 * session-request.js) and the protein-frame check that reads a stored request before it (`storedRecordReverses`,
 * main-session-comparison-frame.js).
 * @param {unknown} region
 * @returns {{ start: unknown, end: unknown, reverseComplement: unknown, column: null }[]}
 */
export const storedRegionWindows = (region) => {
  if (!isPlainObject(region)) return [];
  // A missing field reads as null, and only two integers can be the reversed span, as Python reads them.
  const start = region.start ?? null;
  const end = region.end ?? null;
  const reversed = Number.isInteger(start) && Number.isInteger(end) && /** @type {number} */ (start) > /** @type {number} */ (end); // both integers here
  return [reversed
    ? { start: end, end: start, reverseComplement: true, column: null }
    : { start, end, reverseComplement: region.reverseComplement ?? null, column: null }];
};

/**
 * The orientation and windows of a request record of schema 10, as the readers below read them.
 * @typedef {{ windows: readonly { reverseComplement: boolean }[], presentation: { reverseComplement: boolean } }} CanonicalRecordOrientation
 */

/**
 * Whether a request record is drawn reverse-complemented: by its presentation, or by its windows, all reversed
 * (Python's `RecordInput.is_reversed`).
 * @param {CanonicalRecordOrientation} record
 */
export const canonicalRecordReverseComplement = (record) => record.presentation.reverseComplement
  || (record.windows.length > 0 && record.windows.every((window) => window.reverseComplement));

/**
 * Draws a request record in the absolute direction `reverseComplement`: on its windows when it has them, else on
 * the record (Python's `RecordInput.oriented`).
 * @param {{ windows: { reverseComplement: boolean }[], presentation: { reverseComplement: boolean } }} record
 * @param {boolean} reverseComplement
 */
export const writeCanonicalRecordReverseComplement = (record, reverseComplement) => {
  if (typeof reverseComplement !== 'boolean') throw new Error('Record orientation must be boolean.');
  if (record.windows.length) {
    record.windows.forEach((window) => { window.reverseComplement = reverseComplement; });
    record.presentation.reverseComplement = false;
  } else record.presentation.reverseComplement = reverseComplement;
};

/**
 * The orientation fields of a crop form, a Linear File card (`windows`, `region_reverse`, `linearCardCropForm`) or the Circular Region block
 * (`circular_windows`, `circular_reverse`): the windows it lists (several, or one across the origin), each with its
 * own orientation, and its Reverse complement, which orients only a record it lists no window of, as the request
 * writers read them (Load folds a record-wide reverse complement into the listed windows).
 * @template {{ reverseComplement: boolean }} W
 * @typedef {{ windows: readonly W[], reverse: boolean }} CropFormOrientation
 */
/**
 * @template {{ reverseComplement: boolean }} W
 * @param {CropFormOrientation<W>} form
 * @returns {{ windows: W[], presentation: { reverseComplement: boolean } }} new window objects
 */
const cropFormRecord = ({ windows, reverse }) => ({
  windows: windows.map((window) => ({ ...window })), presentation: { reverseComplement: windows.length ? false : reverse }
});

/**
 * Whether a crop form draws its record reverse-complemented (`canonicalRecordReverseComplement`).
 * @template {{ reverseComplement: boolean }} W
 * @param {CropFormOrientation<W>} form
 */
export const cropFormReverseComplement = (form) => canonicalRecordReverseComplement(cropFormRecord(form));

/**
 * The orientation fields that draw a crop form's record in the absolute direction `reverseComplement`
 * (`writeCanonicalRecordReverseComplement`): its listed windows, as new objects, else its Reverse complement.
 * @template {{ reverseComplement: boolean }} W
 * @param {CropFormOrientation<W>} form
 * @param {boolean} reverseComplement
 * @returns {CropFormOrientation<W> & { windows: W[] }}
 */
export const orientedCropForm = (form, reverseComplement) => {
  const record = cropFormRecord(form);
  writeCanonicalRecordReverseComplement(record, reverseComplement);
  return { windows: record.windows, reverse: record.presentation.reverseComplement };
};

/**
 * A Linear File card's orientation fields: the windows it lists and its Reverse complement.
 * @template {{ reverseComplement: boolean }} W
 * @typedef {{ windows: W[], region_reverse: boolean }} LinearCardOrientation
 */
/**
 * The crop form of a Linear File card, as `cropFormReverseComplement` and `orientedCropForm` read it.
 * @template {{ reverseComplement: boolean }} W
 * @param {LinearCardOrientation<W>} card
 * @returns {CropFormOrientation<W>}
 */
export const linearCardCropForm = (card) => ({ windows: card.windows, reverse: card.region_reverse });

/**
 * A Linear File card's orientation fields as a checkpoint or a History step keeps them: the windows it lists, copied,
 * and its Reverse complement. Each window keeps its own orientation.
 * @template {{ reverseComplement: boolean }} W
 * @param {LinearCardOrientation<W>} card
 * @returns {LinearCardOrientation<W>}
 */
export const captureLinearCardOrientation = (card) => ({
  windows: cloneJsonData(card.windows), region_reverse: card.region_reverse
});

/**
 * @param {RecordDisplayRowView} row
 * @param {RecordDisplayRequested} [draft]
 * @param {RecordDisplayContext} [context]
 */
export const effectiveRecordReverseComplement = (row, draft = {}, context = {}) => {
  const inherited = Boolean(context.reverse ?? row.reverse);
  if (Boolean(context.cropped ?? row.cropped)) return inherited;
  const override = draft.reverseComplementOverride ?? null;
  if (override !== null && typeof override !== 'boolean') {
    throw new Error('Reverse-complement override must be boolean or null.');
  }
  return override ?? inherited;
};

/**
 * @param {RecordDisplayRowView} row
 * @param {RecordDisplayRequested} [draft]
 * @param {RecordDisplayContext} [context]
 */
export const requestedRecordTransform = (row, draft = {}, context = {}) => ({
  display: requestedRecordDisplay(row, draft, context),
  reverseComplement: effectiveRecordReverseComplement(row, draft, context)
});

/**
 * @param {unknown} value
 * @returns {boolean | null}
 */
export const requireReverseComplementOverride = (value) => {
  if (value !== null && typeof value !== 'boolean') {
    throw new Error('Reverse-complement override must be boolean or null.');
  }
  return value;
};
