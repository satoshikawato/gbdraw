// @ts-check
// The record display model: source-bound record display rows, the requested
// display start and topology, and the record orientation intent, as the
// request and Session services read them. The Record Display controls
// (app/record-display-options.js) edit it.
import { resolveDisambiguatedRecordSelection } from './record-options.js';

/** @import { RecordEntryInput } from './record-options.js' */
/** @typedef {'circular' | 'linear'} RecordDisplayScope */
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
 * What the row's source says about the draft's record: the row's cropped and
 * reverse state.
 * @typedef {{ cropped?: unknown, reverse?: unknown }} RecordDisplayContext
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

// A record display row or draft names its source input and record. A drawing
// holds the drafts of its own mode's inputs (PD-OI-086): `circular`, or a
// Linear File card's UID.
/** @param {{ sourceUid?: string, selector: string }} identity */
export const recordDisplayKey = ({ sourceUid, selector }) => {
  if (!sourceUid || !/^#[1-9]\d*$/.test(selector)) {
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
// (`buildDisambiguatedRecordEntries`). Returns `drafts` when no draft moves.
/**
 * @template {Array<RecordDisplayDraftIdentity & { recordId?: unknown }>} Drafts the record ID is only compared
 * @param {Drafts} drafts
 * @param {Array<{ key: string, sourceUid: string, selector: string, recordId: string }>} rows
 * @returns {Drafts}
 */
export const bindRecordDisplayDrafts = (drafts, rows) => {
  const taken = new Set(drafts.map(recordDisplayKey));
  let moved = false;
  const bound = drafts.map((draft) => {
    const sourceRows = rows.filter((row) => row.sourceUid === draft.sourceUid);
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

/**
 * The value is untrusted; every field read below is checked before it is trusted.
 * @param {Partial<Record<keyof AnchorIntent, unknown>> | null | undefined} intent
 * @returns {AnchorIntent | null}
 */
export const validateAnchorIntent = (intent) => {
  if (intent === null) return null;
  if (!intent || Object.keys(intent).sort().join(',') !== ANCHOR_INTENT_KEYS.join(',')
    || intent.schema !== 1
    || typeof intent.recordKey !== 'string' || !intent.recordKey || intent.recordKey.includes('\0')
    || typeof intent.biologicalFeatureId !== 'string' || !intent.biologicalFeatureId
    || intent.biologicalFeatureId.includes('\0')
    || !(/** @type {readonly unknown[]} */ (['anchor', 'feature-end'])).includes(intent.placement)
    || (intent.placement === 'anchor'
      ? !(/** @type {readonly unknown[]} */ (['five-prime', 'midpoint', 'three-prime'])).includes(intent.anchor)
      : intent.anchor !== null)
    || !Number.isSafeInteger(intent.offsetBp)
    || typeof intent.orientForward !== 'boolean') {
    throw new Error('Invalid record display anchor intent.');
  }
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
// split moves each row into its mode's drawing. A row's selector is text
// (OV-575); `textSelector: false` reads a row of a Session Python does not
// check (`validateCurrentWriterActiveConfig` owns which), whose selector is
// only matched against `#N` by `recordDisplayKey`.
/**
 * @param {unknown} drafts
 * @param {{ scoped?: boolean, textSelector?: boolean }} [options]
 */
export const validateRecordDisplayDrafts = (drafts, { scoped = false, textSelector = true } = {}) => {
  if (!Array.isArray(drafts)) throw new Error('Record display drafts must be an array.');
  const fields = scoped
    ? 'anchorIntent,recordId,reverseComplementOverride,scope,selector,sourceUid,startCoordinate,topologyOverride'
    : 'anchorIntent,recordId,reverseComplementOverride,selector,sourceUid,startCoordinate,topologyOverride';
  const identities = new Set();
  for (const draft of drafts) {
    if (!draft || Object.keys(draft).sort().join(',') !== fields
      || (scoped && !['circular', 'linear'].includes(draft.scope))
      || typeof draft.recordId !== 'string' || draft.recordId.includes('\0')
      || typeof draft.sourceUid !== 'string' || draft.sourceUid.includes('\0')
      || (textSelector && typeof draft.selector !== 'string')
      || (draft.topologyOverride !== null && typeof draft.topologyOverride !== 'boolean')
      || (draft.reverseComplementOverride !== null
        && typeof draft.reverseComplementOverride !== 'boolean')
      || (draft.startCoordinate !== null && (!Number.isSafeInteger(draft.startCoordinate) || draft.startCoordinate < 1))) {
      throw new Error('Invalid record display draft; only source-bound requested intent is supported.');
    }
    validateAnchorIntent(draft.anchorIntent);
    const key = scoped ? JSON.stringify([draft.scope, recordDisplayKey(draft)]) : recordDisplayKey(draft);
    if (identities.has(key)) throw new Error('Duplicate record display draft identity.');
    identities.add(key);
  }
  return drafts;
};

/**
 * @param {RecordDisplayRowView} row
 * @param {RecordDisplayRequested} [draft] the row's requested display
 * @param {RecordDisplayContext} [context]
 */
export const recordDisplaySurface = (row, draft = {}, { cropped = false, reverse = false } = {}) => {
  const override = draft.topologyOverride ?? null;
  if (override !== null && typeof override !== 'boolean') throw new Error('Topology override must be boolean or null.');
  const effectiveCircular = override ?? row.detectedTopology === 'circular';
  const lengthKnown = Number.isSafeInteger(row.recordLength) && /** @type {number} */ (row.recordLength) > 0;
  const startEnabled = effectiveCircular && lengthKnown && !cropped;
  const currentStart = reverse ? (lengthKnown ? row.recordLength : null) : 1;
  return {
    effectiveCircular,
    startEnabled,
    disabledReason: cropped ? 'Display start is unavailable for a cropped record.'
      : !lengthKnown ? 'Record length is unavailable.'
        : !effectiveCircular ? 'Display start requires a circular record.' : '',
    currentStart,
    // OV-374: the Display start field shows the requested start, kept while
    // the field is disabled, or else the current one; never an empty field.
    displayStart: draft.startCoordinate ?? currentStart
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
 * @param {{ region?: { reverseComplement?: unknown } | null, presentation?: { reverseComplement?: unknown } | null } | null | undefined} record
 */
export const canonicalRecordReverseComplement = (record) => Boolean(
  record?.region ? record.region.reverseComplement : record?.presentation?.reverseComplement
);

/**
 * @param {{ region?: { reverseComplement?: boolean } | null, presentation: { reverseComplement?: boolean } }} record
 * @param {boolean} reverseComplement
 */
export const writeCanonicalRecordReverseComplement = (record, reverseComplement) => {
  if (typeof reverseComplement !== 'boolean') throw new Error('Record orientation must be boolean.');
  if (record.region) {
    record.region.reverseComplement = reverseComplement;
    record.presentation.reverseComplement = false;
  } else record.presentation.reverseComplement = reverseComplement;
};

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
