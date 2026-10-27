// @ts-check
import { downloadSafeName } from '../utils/download-names.js';
import { CLI_RECORD_KEY_PREFIX, EXPANDED_RECORD_KEY_SEPARATOR } from '../contracts/schema.generated.js';
import { featureEditsOf, setFeatureOverrideRow, setFeaturePlacementRow } from './drawing/feature-edits.js';
/** @import { FeatureEditsSlice } from './drawing/feature-edits.js' */
import { paintOf, setFeatureStroke } from './drawing/paint.js';
/** @import { PaintSlice } from './drawing/paint.js' */
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
 * The CLI's key of its n-th record input, `record-<n>`
 * (gbdraw/api/record_planning.py; a twin, tools/twins.json
 * `concern:feature-instance-record-key`).
 * @param {number} position 1-based position of the record input.
 */
export const cliRecordKey = (position) => `${CLI_RECORD_KEY_PREFIX}${position}`;

/**
 * The request key of a Circular record drawn with others (a canvas or a
 * batch): its file position as the CLI keys it, so it keeps its key while
 * others go OFF.
 * @param {number} sourceIndex 0-based position in the input file.
 */
export const circularRecordRequestKey = (sourceIndex) => cliRecordKey(sourceIndex + 1);

/**
 * The record key of the `k`-th record (1-based) of a record input keyed
 * `baseKey` that expands to several records, as Python names it
 * (`resolve_record_inputs` in gbdraw/api/record_planning.py): `<baseKey>:<k>`;
 * `baseKey` itself when `k` is null (the input draws one record). The one
 * JavaScript composer of that form (a twin, tools/twins.json
 * `concern:feature-instance-record-key`).
 * @param {string} baseKey
 * @param {number | null} k
 */
export const expandedRecordKey = (baseKey, k) => (k === null ? baseKey : `${baseKey}${EXPANDED_RECORD_KEY_SEPARATOR}${k}`);

/**
 * The `k` of a record key `expandedRecordKey(baseKey, k)` composes, or null
 * when `recordKey` is not one of `baseKey`'s expansions.
 * @param {string} recordKey
 * @param {string} baseKey
 * @returns {number | null}
 */
export const recordKeyExpansion = (recordKey, baseKey) => {
  const prefix = `${baseKey}${EXPANDED_RECORD_KEY_SEPARATOR}`;
  const k = recordKey.startsWith(prefix) ? recordKey.slice(prefix.length) : '';
  return /^[1-9]\d*$/.test(k) ? Number(k) : null;
};

/**
 * A record's identity: the name its feature edits use in every view, from its
 * input and its position in the input's file (Owner D-02, D-04). Circular: the
 * canvas key `record-<k>`. Linear: the key the record has when its card draws
 * every record of its file, `<card uid>:<k>`, or the card uid for a file of one
 * record. A record's rows live under its request key while it is drawn and
 * under its identity while it is not (R2).
 * @param {{ mode: string, uid: string, sourceIndex: number, recordCount: number }} record
 *   `sourceIndex`: 0-based position in the file; `recordCount`: the file's records.
 */
export const recordIdentity = ({ mode, uid, sourceIndex, recordCount }) => (mode === 'circular'
  ? circularRecordRequestKey(sourceIndex)
  : expandedRecordKey(uid, recordCount > 1 ? sourceIndex + 1 : null));

/**
 * The request key of a Circular record drawn alone (a single presentation):
 * its preserved key, or one built from its record ID and selector. A preserved
 * canvas key that names another record (`record-<n>` of a CLI or Python
 * Session) is not kept: it is that record's identity (DESIGN A3).
 * @param {{ recordKey?: unknown, recordId?: unknown, selector?: unknown, sourceIndex?: unknown } | null | undefined} record
 */
export const circularSingleRecordKey = (record) => {
  const preserved = String(record?.recordKey || '').trim();
  const position = preserved.startsWith(CLI_RECORD_KEY_PREFIX) ? preserved.slice(CLI_RECORD_KEY_PREFIX.length) : '';
  const namesAnotherRecord = /^[1-9]\d*$/.test(position) && Number.isInteger(record?.sourceIndex)
    && Number(position) !== Number(record?.sourceIndex) + 1;
  if (preserved && !namesAnotherRecord) return preserved;
  const recordId = downloadSafeName(record?.recordId, 'record');
  const selector = downloadSafeName(record?.selector, '1');
  return `circular-${recordId}-${selector}`;
};

/**
 * The request keys that can name a Circular record: the key it has when drawn
 * with others and the key it has when drawn alone.
 * @param {{ sourceIndex: number, recordKey?: unknown, recordId?: unknown, selector?: unknown }} record
 * @returns {string[]}
 */
export const circularRecordRequestKeys = (record) => [
  circularRecordRequestKey(record.sourceIndex), circularSingleRecordKey(record)
];

/**
 * The records a drawing's OFF records keep their rows under, as request
 * records a row is matched against (`rowBelongsToRequest`): a Linear card's
 * uid with its `<uid>:<k>` identities, or a Circular record's `record-N`.
 * @param {string} mode 'circular', or a Linear mode.
 * @param {Iterable<unknown> | null | undefined} recordsOff
 * @returns {{ recordKey: string, cardinality: 'all' | 'exactly_one' }[]}
 */
export const offRecordRequestRecords = (mode, recordsOff) => [...offSet(recordsOff)]
  .map((key) => (mode !== 'circular' ? { recordKey: key, cardinality: /** @type {const} */ ('all') }
    : { recordKey: CIRCULAR_RECORD_KEY.test(key) ? circularRecordRequestKey(Number(key.slice(1)) - 1) : '',
      cardinality: /** @type {const} */ ('exactly_one') }))
  .filter(({ recordKey }) => recordKey);

/**
 * A record of a source file as a record list names it.
 * @typedef {{ sourceIndex?: unknown, selector?: unknown, recordId?: unknown, record_id?: unknown }} RecordSourceEntry
 * The record lists the identities of a request's records are read from: a
 * Linear card's records by card uid, the Circular input's records; each list
 * in file order, as `buildDisambiguatedRecordEntries` (services/record-options.js)
 * gives it, so `selector` and `recordId` are those a request key is built from.
 * @typedef {{ mode: string, linear?: ReadonlyMap<string, readonly RecordSourceEntry[]>,
 *   circular?: readonly RecordSourceEntry[] }} RecordIdentitySources
 * A record selector or region of a request record.
 * @typedef {{ kind?: string, index?: number, value?: string } | null | undefined} RequestRecordSelector
 * A record of a render request, as its identities are read from it.
 * @typedef {{ recordKey: string, cardinality?: string, selector?: RequestRecordSelector,
 *   region?: { selector?: RequestRecordSelector } | null }} IdentityRequestRecord
 * One drawn record of a request: the record key its features have in that
 * request, its identity, and (Circular) the other keys a 0.14.2 view switch
 * may have left its rows under.
 * @typedef {{ requestKey: string, identity: string, aliases: readonly string[] }} DrawnRecordIdentity
 */

/**
 * The 0-based positions of a source's records, by record ID (a duplicated ID names none).
 * @param {readonly RecordSourceEntry[]} entries
 */
const positionsOf = (entries) => {
  const byId = new Map();
  entries.forEach((entry, index) => {
    const recordId = cleanKey(entry?.recordId ?? entry?.record_id);
    const sourceIndex = Number.isInteger(entry?.sourceIndex) ? Number(entry.sourceIndex) : index;
    byId.set(recordId, byId.has(recordId) ? null : sourceIndex);
  });
  return byId;
};

/**
 * The 0-based position in its file of the one record a request record selects,
 * or null when the selection does not name one record of `entries`.
 * @param {IdentityRequestRecord} record
 * @param {readonly RecordSourceEntry[]} entries
 */
const selectedPosition = (record, entries) => {
  const selector = record.region?.selector || record.selector;
  if (selector?.kind === 'recordIndex' && Number.isInteger(selector.index)) {
    return Number(selector.index) < entries.length ? Number(selector.index) : null;
  }
  if (selector?.kind === 'recordId') return positionsOf(entries).get(cleanKey(selector.value)) ?? null;
  return entries.length === 1 ? 0 : null;
};

/**
 * The drawn records of a request with their identities (DESIGN 3, A1): each
 * record's position comes from the request's own selector, read against the
 * record lists of that request's sources, never from a card's current
 * selection. A record whose source list is unknown has no identity, so its
 * rows keep their key. A Linear card that draws every record of a file of
 * several records draws each one under its identity.
 * @param {readonly IdentityRequestRecord[]} records
 * @param {RecordIdentitySources} sources
 * @returns {DrawnRecordIdentity[]}
 */
export const drawnRecordIdentities = (records, sources) => records.flatMap((record) => {
  const recordKey = cleanKey(record?.recordKey);
  if (sources.mode === 'circular') {
    const entries = sources.circular || [];
    const position = selectedPosition(record, entries);
    if (position === null) return [];
    const entry = entries[position] || {};
    return [{
      requestKey: recordKey,
      identity: recordIdentity({ mode: 'circular', uid: '', sourceIndex: position, recordCount: entries.length }),
      aliases: circularRecordRequestKeys({
        sourceIndex: position, recordId: entry.recordId ?? entry.record_id, selector: entry.selector
      })
    }];
  }
  const linear = sources.linear || new Map();
  const base = linear.has(recordKey) ? recordKey
    : [...linear.keys()].find((uid) => recordKeyExpansion(recordKey, uid) !== null);
  const entries = base === undefined ? undefined : linear.get(base);
  if (base === undefined || !entries || entries.length === 0) return [];
  const identityOf = (/** @type {number} */ sourceIndex) => recordIdentity({
    mode: 'linear', uid: base, sourceIndex, recordCount: entries.length
  });
  if (record.cardinality === 'all' && !record.selector && !record.region) {
    return entries.map((/** @type {RecordSourceEntry} */ _entry, /** @type {number} */ sourceIndex) => {
      const identity = identityOf(sourceIndex);
      return { requestKey: identity, identity, aliases: [identity] };
    });
  }
  const expansion = recordKeyExpansion(recordKey, base);
  const position = expansion !== null ? expansion - 1 : selectedPosition(record, entries);
  if (position === null || position >= entries.length) return [];
  const identity = identityOf(position);
  return [{ requestKey: recordKey, identity, aliases: [identity] }];
});

/**
 * The record key changes of one Generate (DESIGN 2-3, A3): `from` key ->
 * `to` key, in priority order. A record's rows move from its key in the
 * previous request (its identity when it was not drawn) to its key in the new
 * request (its identity when it is not drawn). Then the aliases of each record
 * the new request draws bring the rows a 0.14.2 view switch left under its
 * other keys, after its own rows: a row moved from its previous key, or else
 * the row already under its key (listed as a move to itself), wins field by
 * field; a row left under the new key by an earlier view comes last. A key that names another record in either
 * request is never taken as an alias, and a record does not wait under an
 * identity another record uses as its request key. The moves apply together
 * from one snapshot, so a permutation is safe.
 * @param {readonly DrawnRecordIdentity[]} previousDrawn
 * @param {readonly DrawnRecordIdentity[]} nextDrawn
 * @returns {Map<string, string>}
 */
export const recordKeyMoves = (previousDrawn, nextDrawn) => {
  const previousKeyOf = new Map(previousDrawn.map(({ identity, requestKey }) => [identity, requestKey]));
  const nextKeyOf = new Map(nextDrawn.map(({ identity, requestKey }) => [identity, requestKey]));
  const previousKeys = new Set(previousKeyOf.values());
  const nextKeys = new Set(nextKeyOf.values());
  /** @type {Map<string, string>} */
  const moves = new Map();
  new Set([...previousKeyOf.keys(), ...nextKeyOf.keys()]).forEach((identity) => {
    const from = previousKeyOf.get(identity) ?? identity;
    const to = nextKeyOf.get(identity) ?? identity;
    if (from === to) return;
    // A row under another record's request key belongs to that record.
    if (!previousKeyOf.has(identity) && previousKeys.has(from)) return;
    if (!nextKeyOf.has(identity) && nextKeys.has(to)) return;
    moves.set(from, to);
  });
  const movedInto = new Set(moves.values());
  nextDrawn.forEach(({ requestKey, aliases }) => aliases.forEach((alias) => {
    if (alias === requestKey || moves.has(alias) || nextKeys.has(alias) || previousKeys.has(alias)) return;
    // The rows already under the record's key, which no move brings there,
    // are its own: they come first (`requestKey` -> `requestKey`).
    if (!movedInto.has(requestKey) && !moves.has(requestKey)) moves.set(requestKey, requestKey);
    moves.set(alias, requestKey);
  }));
  return moves;
};

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
 * The edits of one drawing, which `removeRecordEdits` removes: its per-feature
 * edits and Feature placements through their slice (`drawing.slices.featureEdits`).
 * @typedef {object} RecordEditData
 * @property {FeatureEditsSlice} [featureEdits]
 * @property {PaintSlice} [paint] Its paint (`drawing.slices.paint`): the feature strokes, keyed `<recordKey>\0…`.
 * @property {readonly { readonly annotations: readonly { readonly target?: { readonly kind?: string, readonly recordKey?: unknown },
 *   readonly metadata?: Readonly<Record<string, unknown>> | null }[] }[]} [annotationSets]
 * @property {(annotation: object) => void} [removeAnnotation] Removes an annotation of `annotationSets` (the writer of their slice).
 * @property {readonly { readonly recordKey: string, readonly remove: () => void }[]} [colorRuleEdits] The specific
 *   color rows that name a feature instance, by the instance's record key; `remove` takes the row out of the table.
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
    const base = [...expanding].find((key) => recordKeyExpansion(recordKey, key) !== null);
    return base === undefined ? undefined : byRequestKey.get(base);
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
  /**
   * @param {Readonly<Record<string, { readonly recordKey: string }>>} rows
   * @param {(draftKey: string) => void} remove
   */
  const rowsOf = (rows, remove) => Object.entries(rows).forEach(([draftKey, row]) => {
    const key = match.ofRequestKey(row.recordKey);
    if (key !== undefined) visit(key, () => remove(draftKey));
  });
  const { featureEdits } = data;
  if (featureEdits) {
    const { featureOverrides, featurePlacementOverrides } = featureEditsOf(featureEdits);
    rowsOf(featureOverrides, (draftKey) => setFeatureOverrideRow(featureEdits, draftKey, null));
    rowsOf(featurePlacementOverrides, (draftKey) => setFeaturePlacementRow(featureEdits, draftKey, null));
  }
  const { paint } = data;
  Object.keys(paint ? paintOf(paint).featureStrokeOverrides : {}).forEach((strokeKey) => {
    const key = strokeKey.includes('\0') ? match.ofRequestKey(strokeKey.split('\0')[0]) : undefined;
    if (paint && key !== undefined) visit(key, () => setFeatureStroke(paint, strokeKey, null));
  });
  (data.colorRuleEdits || []).forEach(({ recordKey, remove }) => {
    const key = match.ofRequestKey(recordKey);
    if (key !== undefined) visit(key, remove);
  });
  (data.annotationSets || []).forEach((set) => {
    [...set.annotations].forEach((item) => {
      const key = item?.target?.kind === 'featureIdentity'
        ? match.ofRequestKey(item.target.recordKey)
        : match.ofBinding(item?.metadata?.[data.annotationBindingField]);
      if (key !== undefined) visit(key, () => data.removeAnnotation?.(item));
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
