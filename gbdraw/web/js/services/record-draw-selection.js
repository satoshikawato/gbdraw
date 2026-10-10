// @ts-check
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
