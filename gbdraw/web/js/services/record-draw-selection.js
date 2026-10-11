// @ts-check
import { downloadSafeName } from '../utils/download-names.js';
// Which records of a drawing are drawn (D-01..D-10, record-selection). A
// drawing holds the keys of its OFF records in `recordsOff`: a Linear record
// card's `uid`, or a Circular record's source selector `#N` (1-based position
// in the input file). This module is the one place that answers "is this
// record drawn"; it holds no state.

/**
 * @typedef {string} RecordDrawKey A Linear card uid or a Circular source selector `#N`.
 * @typedef {'file' | 'id-asc' | 'id-desc' | 'length-desc' | 'length-asc'} RecordListSort
 * @typedef {{ key: RecordDrawKey, recordId: string, length: number | null, hasSettings?: boolean }} RecordListEntry
 *   One record of a source, in file order.
 * @typedef {{ key: RecordDrawKey, recordId: string, length: number | null, drawn: boolean,
 *   hasSettings: boolean, position: number }} RecordListRow One row of the record list (`position`: file order).
 */

export const RECORD_LIST_SORTS = Object.freeze(
  /** @type {RecordListSort[]} */ (['file', 'id-asc', 'id-desc', 'length-desc', 'length-asc'])
);

const CIRCULAR_RECORD_KEY = /^#[1-9][0-9]*$/;

/** @param {unknown} value */
const cleanKey = (value) => String(value ?? '').trim();

/** @param {Iterable<unknown> | null | undefined} recordsOff */
const offSet = (recordsOff) => new Set([...(recordsOff || [])].map(cleanKey).filter(Boolean));

/**
 * The Linear record cards a Generate draws, in file order.
 * @template {{ uid?: unknown }} T
 * @param {readonly T[] | null | undefined} sequences
 * @param {Iterable<unknown> | null | undefined} recordsOff
 * @returns {T[]}
 */
export const drawnLinearSequences = (sequences, recordsOff) => {
  const list = Array.isArray(sequences) ? sequences : [];
  const off = offSet(recordsOff);
  return off.size === 0 ? [...list] : list.filter((sequence) => !off.has(cleanKey(sequence?.uid)));
};

/**
 * The uids of the OFF cards that still exist.
 * @param {readonly { uid?: unknown }[] | null | undefined} sequences
 * @param {Iterable<unknown> | null | undefined} recordsOff
 * @returns {Set<string>}
 */
export const omittedLinearUids = (sequences, recordsOff) => {
  const off = offSet(recordsOff);
  return new Set((Array.isArray(sequences) ? sequences : [])
    .map((sequence) => cleanKey(sequence?.uid))
    .filter((uid) => off.has(uid)));
};

/**
 * The key of a discovered Circular record: its source selector.
 * @param {{ selector?: unknown, sourceIndex?: unknown }} entry
 * @returns {RecordDrawKey}
 */
export const circularRecordDrawKey = (entry) => {
  const selector = cleanKey(entry?.selector);
  if (CIRCULAR_RECORD_KEY.test(selector)) return selector;
  return Number.isInteger(entry?.sourceIndex) ? `#${Number(entry.sourceIndex) + 1}` : selector;
};

/**
 * The Circular records a Multi-Record Canvas or separate diagrams draw; each
 * entry keeps its `sourceIndex`.
 * @template {{ selector?: unknown, sourceIndex?: unknown }} T
 * @param {readonly T[] | null | undefined} entries
 * @param {Iterable<unknown> | null | undefined} recordsOff
 * @returns {T[]}
 */
export const circularRecordsToDraw = (entries, recordsOff) => {
  const list = Array.isArray(entries) ? entries : [];
  const off = offSet(recordsOff);
  return off.size === 0 ? [...list] : list.filter((entry) => !off.has(circularRecordDrawKey(entry)));
};

/**
 * The request key of a Circular record drawn with others (a canvas or a
 * batch): its file position, so it keeps its key while others go OFF.
 * @param {number} sourceIndex 0-based position in the input file.
 */
export const circularRecordRequestKey = (sourceIndex) => `record-${sourceIndex + 1}`;

/**
 * The request key of a Circular record drawn alone (a single presentation):
 * its preserved key, or one built from its record ID and selector.
 * @param {{ recordKey?: unknown, recordId?: unknown, selector?: unknown } | null | undefined} record
 */
export const circularSingleRecordKey = (record) => {
  const preserved = String(record?.recordKey || '').trim();
  if (preserved) return preserved;
  const recordId = downloadSafeName(record?.recordId, 'record');
  const selector = downloadSafeName(record?.selector, '1');
  return `circular-${recordId}-${selector}`;
};

/**
 * The request record keys of a drawing's OFF records: a Linear card's uid, or
 * a Circular record's `record-N`.
 * @param {string} mode 'circular', or a Linear mode.
 * @param {Iterable<unknown> | null | undefined} recordsOff
 * @returns {string[]}
 */
export const offRecordRequestKeys = (mode, recordsOff) => [...offSet(recordsOff)]
  .map((key) => (mode !== 'circular' ? key
    : CIRCULAR_RECORD_KEY.test(key) ? circularRecordRequestKey(Number(key.slice(1)) - 1) : ''))
  .filter(Boolean);

/**
 * Whether a key has the shape of its mode's record key.
 * @param {'circular' | 'linear'} mode
 * @param {unknown} key
 */
export const isRecordDrawKey = (mode, key) => (
  typeof key === 'string' && key.trim() === key && key !== ''
  && (mode !== 'circular' || CIRCULAR_RECORD_KEY.test(key))
);

/**
 * `recordsOff` after turning `keys` ON (`drawn`) or OFF; order kept, no duplicates.
 * @param {Iterable<unknown> | null | undefined} recordsOff
 * @param {Iterable<unknown>} keys
 * @param {boolean} drawn
 * @returns {RecordDrawKey[]}
 */
export const nextRecordsOff = (recordsOff, keys, drawn) => {
  const changed = offSet(keys);
  const current = offSet(recordsOff);
  return drawn
    ? [...current].filter((key) => !changed.has(key))
    : [...current, ...[...changed].filter((key) => !current.has(key))];
};

/**
 * `recordsOff` without the keys of records that no longer exist.
 * @param {Iterable<unknown> | null | undefined} recordsOff
 * @param {Iterable<unknown>} liveKeys
 * @returns {RecordDrawKey[]}
 */
export const pruneRecordsOff = (recordsOff, liveKeys) => {
  const live = offSet(liveKeys);
  return [...offSet(recordsOff)].filter((key) => live.has(key));
};

/**
 * Whether turning `keys` OFF would leave a source (a Linear File, or the
 * Circular input) with no drawn record (D-06). Turning records ON never does.
 * @param {{ sourceKeys: Iterable<unknown>, recordsOff: Iterable<unknown> | null | undefined,
 *   keys: Iterable<unknown>, drawn: boolean }} change
 * @returns {boolean}
 */
export const offChangeLeavesSourceEmpty = ({ sourceKeys, recordsOff, keys, drawn }) => {
  if (drawn) return false;
  const source = [...offSet(sourceKeys)];
  if (source.length === 0) return false;
  const off = new Set(nextRecordsOff(recordsOff, keys, false));
  return source.every((key) => off.has(key));
};

const chunkPattern = /(\d+)/;
/**
 * Natural order of record IDs: digit runs compare as numbers (contig_2 before
 * contig_10), the rest case-insensitively.
 * @param {string} left
 * @param {string} right
 */
export const compareRecordIdsNaturally = (left, right) => {
  const a = String(left).split(chunkPattern);
  const b = String(right).split(chunkPattern);
  for (let index = 0; index < Math.min(a.length, b.length); index += 1) {
    if (a[index] === b[index]) continue;
    if (index % 2 === 1) {
      const difference = Number(a[index]) - Number(b[index]);
      if (difference !== 0) return difference;
      // Equal numbers with different zero padding: the shorter run first.
      return a[index].length - b[index].length;
    }
    const folded = a[index].toLowerCase().localeCompare(b[index].toLowerCase(), 'en');
    if (folded !== 0) return folded;
    return a[index] < b[index] ? -1 : 1;
  }
  return a.length - b.length;
};

/** @param {number | null} value */
const knownLength = (value) => Number.isFinite(value) && Number(value) > 0;

/**
 * The rows of a source's record list: filtered by a case-insensitive
 * substring of the record ID and sorted for display (D-10). Sorting changes
 * the list only, never the drawing order.
 * @param {{ records: readonly RecordListEntry[] | null | undefined,
 *   recordsOff: Iterable<unknown> | null | undefined, query?: string, sort?: string }} input
 * @returns {RecordListRow[]}
 */
export const recordListRows = ({ records, recordsOff, query = '', sort = 'file' }) => {
  const off = offSet(recordsOff);
  const needle = String(query || '').trim().toLowerCase();
  const rows = (Array.isArray(records) ? records : []).map((record, position) => ({
    key: cleanKey(record?.key),
    recordId: String(record?.recordId ?? ''),
    length: knownLength(record?.length) ? Number(record.length) : null,
    drawn: !off.has(cleanKey(record?.key)),
    hasSettings: record?.hasSettings === true,
    position
  })).filter((row) => !needle || row.recordId.toLowerCase().includes(needle));
  /** @type {(left: RecordListRow, right: RecordListRow) => number} */
  let compare;
  switch (sort) {
    case 'id-asc': compare = (left, right) => compareRecordIdsNaturally(left.recordId, right.recordId); break;
    case 'id-desc': compare = (left, right) => compareRecordIdsNaturally(right.recordId, left.recordId); break;
    case 'length-desc':
    case 'length-asc': {
      const sign = sort === 'length-desc' ? -1 : 1;
      // A record without a known length sorts last either way.
      compare = (left, right) => (left.length === null || right.length === null
        ? Number(left.length === null) - Number(right.length === null)
        : sign * (left.length - right.length));
      break;
    }
    default: return rows;
  }
  return rows.sort((left, right) => compare(left, right) || left.position - right.position);
};

// D-08 Delete settings: what a record keeps while it is OFF and what Delete
// settings removes. The card fields return to these values (the record
// fields Reset Settings resets, plus LOSAT Gencode).
export const RECORD_SETTINGS_DEFAULTS = Object.freeze({
  definition: '', record_subtitle: '', region_start: null, region_end: null, region_reverse: false, losat_gencode: 1
});

/**
 * Whether a Linear record card holds a value that Delete settings resets.
 * @param {Record<string, unknown> | null | undefined} card
 */
export const linearCardHasSettings = (card) => Object.entries(RECORD_SETTINGS_DEFAULTS)
  .some(([field, fallback]) => (card?.[field] ?? fallback) !== fallback);

/**
 * One record whose edits Delete settings finds or removes.
 * @typedef {object} RecordSettingsOwner
 * @property {RecordDrawKey} key
 * @property {readonly string[]} requestKeys The request record keys its feature edits use: a Linear
 *   card's uid (which also owns its `<uid>:<n>` records), a Circular record's `record-N` and its
 *   single-record key.
 * @property {boolean} [ownsExpansions] Linear: the key also owns `<key>:<n>`.
 * @property {readonly string[]} bindingKeys The record-catalog keys an annotation binds it by.
 * @property {string} displaySource The `sourceUid` of its record display drafts.
 * @property {string | null} displaySelector Circular: its `#N`; Linear: null (every draft of the card).
 *
 * The edits of one drawing, mutated in place by `removeRecordEdits`.
 * @typedef {object} RecordEditData
 * @property {Record<string, { recordKey?: unknown }>} [featureOverrides]
 * @property {Record<string, { recordKey?: unknown }>} [featurePlacementOverrides]
 * @property {Record<string, unknown>} [featureStrokeOverrides] Keyed `<recordKey>\0…`.
 * @property {{ annotations: { target?: { kind?: string, recordKey?: unknown },
 *   metadata?: Record<string, unknown> | null }[] }[]} [annotationSets]
 * @property {{ sourceUid?: unknown, selector?: unknown }[]} [recordDisplayDrafts]
 * @property {{ queryUid?: unknown, subjectUid?: unknown }[]} [comparisonEdges] Linear edges, by card uid.
 * @property {string} annotationBindingField The annotation metadata field of a record binding.
 */

/** @param {readonly RecordSettingsOwner[]} owners */
const recordEditMatcher = (owners) => {
  /** @type {Map<string, RecordDrawKey>} */
  const byRequestKey = new Map();
  /** @type {Map<string, RecordDrawKey>} */
  const byBinding = new Map();
  /** @type {Map<string, RecordDrawKey>} */
  const byDisplay = new Map();
  /** @type {Set<string>} */
  const expanding = new Set();
  owners.forEach((owner) => {
    owner.requestKeys.forEach((requestKey) => {
      if (!requestKey) return;
      byRequestKey.set(requestKey, owner.key);
      if (owner.ownsExpansions) expanding.add(requestKey);
    });
    owner.bindingKeys.forEach((bindingKey) => { if (bindingKey) byBinding.set(bindingKey, owner.key); });
    byDisplay.set(JSON.stringify([owner.displaySource, owner.displaySelector]), owner.key);
  });
  /** @param {unknown} value @returns {RecordDrawKey | undefined} */
  const ofRequestKey = (value) => {
    const recordKey = cleanKey(value);
    const direct = byRequestKey.get(recordKey);
    if (direct !== undefined) return direct;
    const match = /^(.*):[1-9]\d*$/.exec(recordKey);
    return match && expanding.has(match[1]) ? byRequestKey.get(match[1]) : undefined;
  };
  return {
    ofRequestKey,
    /** @param {unknown} value */
    ofBinding: (value) => byBinding.get(cleanKey(value)),
    /** @param {{ sourceUid?: unknown, selector?: unknown }} draft */
    ofDisplayDraft: (draft) => byDisplay.get(JSON.stringify([cleanKey(draft?.sourceUid), null]))
      ?? byDisplay.get(JSON.stringify([cleanKey(draft?.sourceUid), cleanKey(draft?.selector)])),
    /** @param {{ queryUid?: unknown, subjectUid?: unknown }} edge */
    ofEdge: (edge) => byRequestKey.get(cleanKey(edge?.queryUid)) ?? byRequestKey.get(cleanKey(edge?.subjectUid))
  };
};

/**
 * Walks the edits of `owners` in a drawing's data: `visit(key, remove)` for each.
 * @param {readonly RecordSettingsOwner[]} owners
 * @param {RecordEditData} data
 * @param {(key: RecordDrawKey, remove: () => void) => void} visit
 */
const forEachRecordEdit = (owners, data, visit) => {
  const match = recordEditMatcher(owners);
  /** @param {Record<string, { recordKey?: unknown }> | undefined} rows */
  const rowsOf = (rows) => Object.entries(rows || {}).forEach(([draftKey, row]) => {
    const key = match.ofRequestKey(row?.recordKey);
    if (key !== undefined) visit(key, () => { delete /** @type {Record<string, unknown>} */ (rows)[draftKey]; });
  });
  rowsOf(data.featureOverrides);
  rowsOf(data.featurePlacementOverrides);
  Object.keys(data.featureStrokeOverrides || {}).forEach((strokeKey) => {
    const key = strokeKey.includes('\0') ? match.ofRequestKey(strokeKey.split('\0')[0]) : undefined;
    if (key !== undefined) visit(key, () => { delete /** @type {Record<string, unknown>} */ (data.featureStrokeOverrides)[strokeKey]; });
  });
  (data.annotationSets || []).forEach((set) => {
    [...set.annotations].forEach((item) => {
      const key = item?.target?.kind === 'featureIdentity'
        ? match.ofRequestKey(item.target.recordKey)
        : match.ofBinding(item?.metadata?.[data.annotationBindingField]);
      if (key !== undefined) visit(key, () => { set.annotations.splice(set.annotations.indexOf(item), 1); });
    });
  });
  /**
   * @template T
   * @param {T[] | undefined} list
   * @param {(entry: T) => RecordDrawKey | undefined} of
   */
  const listOf = (list, of) => [...(list || [])].forEach((entry) => {
    const key = of(entry);
    if (key !== undefined && list) visit(key, () => { list.splice(list.indexOf(entry), 1); });
  });
  listOf(data.recordDisplayDrafts, match.ofDisplayDraft);
  listOf(data.comparisonEdges, match.ofEdge);
};

/**
 * The records of `owners` that have a feature edit, Feature placement, stroke,
 * annotation, record display draft, or comparison pair in the drawing (D-08).
 * @param {readonly RecordSettingsOwner[]} owners
 * @param {RecordEditData} data
 * @returns {Set<RecordDrawKey>}
 */
export const recordsWithEdits = (owners, data) => {
  /** @type {Set<RecordDrawKey>} */
  const keys = new Set();
  forEachRecordEdit(owners, data, (key) => keys.add(key));
  return keys;
};

/**
 * Removes those edits in place (Delete settings, D-08); returns how many.
 * @param {readonly RecordSettingsOwner[]} owners
 * @param {RecordEditData} data
 */
export const removeRecordEdits = (owners, data) => {
  let removed = 0;
  forEachRecordEdit(owners, data, (_key, remove) => { remove(); removed += 1; });
  return removed;
};
