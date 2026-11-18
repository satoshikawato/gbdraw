// @ts-check
import { downloadSafeName } from '../utils/download-names.js';
/** @import { StoredRecordTranslation } from './linear-record-layout.js' */
import {
  CLI_RECORD_KEY_PREFIX, EXPANDED_RECORD_KEY_SEPARATOR, RECORD_INDEX_SELECTOR_PATTERN, RECORD_ROW_PATTERN
} from '../contracts/schema.generated.js';
import { featureEditsOf, setFeatureOverrideRow, setFeaturePlacementRow } from './drawing/feature-edits.js';
/** @import { FeatureEditsSlice } from './drawing/feature-edits.js' */
import { paintOf, setFeatureStroke } from './drawing/paint.js';
/** @import { PaintSlice } from './drawing/paint.js' */
/** @import { InputFileRef } from './history-files.js' */
/** @import { CircularCardMark } from './drawing/inputs.js' */
/** @import { CanonicalRecord } from '../contracts/request.generated.js' */
// Which records of a drawing are drawn (D-01..D-10, record-selection). A
// drawing holds the keys of its OFF records in `recordsOff`: a Linear record
// card's `uid`, or a Circular record's identity (`recordIdentity`: its source
// card's uid and its position in that card's file). This module is the one
// place that answers "is this record drawn"; it holds no state.

/**
 * @typedef {string} RecordDrawKey A Linear card uid or a Circular record identity.
 * @typedef {'file' | 'id-asc' | 'id-desc' | 'length-desc' | 'length-asc'} RecordListSort
 * @typedef {{ key: RecordDrawKey, recordId: string, length: number | null, hasSettings?: boolean }} RecordListEntry
 *   One record of a source, in file order.
 * @typedef {{ key: RecordDrawKey, recordId: string, length: number | null, drawn: boolean,
 *   hasSettings: boolean, position: number }} RecordListRow One row of the record list (`position`: file order).
 */

export const RECORD_LIST_SORTS = Object.freeze(
  /** @type {RecordListSort[]} */ (['file', 'id-asc', 'id-desc', 'length-desc', 'length-asc'])
);

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
 * The key of a discovered Circular record: its identity.
 * @param {{ identity?: unknown }} entry
 * @returns {RecordDrawKey}
 */
export const circularRecordDrawKey = (entry) => cleanKey(entry?.identity);

/**
 * The Circular records a Multi-Record Canvas or separate diagrams draw.
 * @template {{ identity?: unknown }} T
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
const CLI_RECORD_KEY = new RegExp(`^${CLI_RECORD_KEY_PREFIX}[1-9][0-9]*$`);
/** Whether a key is one `cliRecordKey` composes. @param {string} key */
export const isCliRecordKey = (key) => CLI_RECORD_KEY.test(key);

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
 * The record key Generate gives a Linear card (`buildRecords` in session-request.js); a card that
 * draws several records splits into its `expandedRecordKey` members.
 * @param {{ uid?: unknown }} seq @param {number} index
 */
export const linearSeqRecordKey = (seq, index) => String(seq.uid || `record-${index + 1}`);

/**
 * Whether a record key belongs to a Linear card: the card's key, or a split member `<card key>:<k>` of a card that
 * draws several records (`recordKeyExpansion`). The one matcher of a record key to its card
 * (`linearRecordCardIndexes` in session-request.js, `movesARecordSideways`).
 * @param {string} cardKey @param {string} recordKey
 */
export const cardHoldsRecordKey = (cardKey, recordKey) => recordKey === cardKey || recordKeyExpansion(recordKey, cardKey) !== null;

/**
 * A request record as `sameSelectedRecord` reads it: a selector names a record by its ID or its place in its file.
 * @typedef {{ cardinality?: string, selector?: CanonicalRecord['selector'] }} SelectingRecord
 */
/** @param {CanonicalRecord['selector']} a @param {CanonicalRecord['selector']} b */
export const sameRecordSelector = (a, b) => (a === null || b === null ? a === b
  : a.kind === 'recordIndex' ? b.kind === 'recordIndex' && a.index === b.index
    : b.kind === 'recordId' && a.value === b.value);
/**
 * Whether two request records of one source select the same record of it, whatever form each names it in: the same
 * named record, or the record of a file of one record, named in one and in the CLI's form in the other (an
 * `exactly_one` record without a selector, which Python draws only from a file of one record, so the other names
 * that record). The one rule of "which record a Linear record selects" for the card's live match
 * (`cardDrawsRecord`), the reuse of resolved protein results (`canReuseResolvedProteinArtifacts`) and the alignment
 * reset receipt Generate reads (`validateGeneratedSimilarityAlignmentResetReceipt`).
 * @param {SelectingRecord} a @param {SelectingRecord} b
 */
export const sameSelectedRecord = (a, b) => {
  const [selectorA, selectorB] = [a.selector ?? null, b.selector ?? null];
  return sameRecordSelector(selectorA, selectorB)
    || ((selectorA === null) !== (selectorB === null) && a.cardinality === 'exactly_one' && b.cardinality === 'exactly_one');
};

/**
 * Whether a drawn record's stored translation moves it sideways (x other than 0), whatever wrote it: Align at
 * feature, a Similarity alignment, or a loaded Session (an anchored CLI figure, the Python API). One rule for every
 * source (PD-OI-113 scenario revision 2); whether Wrap Rows cuts the record is Python's (`refuse_undrawable_wrap`,
 * whose `x != 0` this repeats: a twin with no generator, held by tests/web/linear-wrap-form.test.mjs). It reads the
 * entries of the drawn cards' records, a split member `<card key>:<k>` included (`cardHoldsRecordKey`), so an entry
 * of a record no longer drawn (its File removed, the record turned off) is not read.
 * @param {readonly StoredRecordTranslation[]} stored
 * @param {readonly { uid?: unknown }[]} drawnSeqs The drawn Linear cards (`drawnLinearSeqs`).
 * @returns {boolean}
 */
export const movesARecordSideways = (stored, drawnSeqs) => {
  const cards = drawnSeqs.map(linearSeqRecordKey);
  return stored.some((entry) => typeof entry?.x === 'number' && entry.x !== 0
    && cards.some((cardKey) => cardHoldsRecordKey(cardKey, String(entry.recordKey))));
};

/**
 * A record's identity: the name its feature edits use in every view, from its
 * input and its position in the input's file (Owner D-02, D-04), in both modes:
 * the key the record has when its card (a Linear record card, a Circular source
 * card) draws every record of its file, `<card uid>:<k>`, or the card uid for a
 * file of one record. A record's rows live under its request key while it is
 * drawn and under its identity while it is not (R2).
 * @param {{ uid: string, sourceIndex: number, recordCount: number }} record
 *   `sourceIndex`: 0-based position in the file; `recordCount`: the file's records.
 */
export const recordIdentity = ({ uid, sourceIndex, recordCount }) => (
  expandedRecordKey(uid, recordCount > 1 ? sourceIndex + 1 : null)
);

/**
 * Which record of a card's new file continues each record of its file then,
 * after a Replace in the same card (D-05, refining D-20): the one rule for
 * every row a record holds (feature edits, color rule instances, annotations,
 * strokes: `recordKeyMoves`; a record-major Depth matrix: `resliceCardRows`).
 * When every record ID of both files is unique, a record continues as the
 * record of the new file with its ID, so a reordered file keeps each record's
 * rows; a record whose ID the new file lacks continues as the record at its
 * place when no ID names that one (one record replacing one record, a new
 * version of the same genome, D-32). Otherwise (a duplicated or missing ID)
 * every record continues by place. -1: the record leaves.
 * @param {readonly unknown[]} previous the record IDs of the file then, in file order
 * @param {readonly unknown[]} next the record IDs of the new file, in file order
 * @returns {number[]} for each record then, the 0-based place of its record in the new file, or -1
 */
export const continuedRecordPlaces = (previous, next) => {
  const ids = (/** @type {readonly unknown[]} */ list) => list.map(cleanKey);
  const [before, after] = [ids(previous), ids(next)];
  const unique = (/** @type {string[]} */ list) => list.every(Boolean) && new Set(list).size === list.length;
  if (!unique(before) || !unique(after)) return before.map((_id, place) => (place < after.length ? place : -1));
  const kept = new Set(before);
  return before.map((id, place) => {
    const byId = after.indexOf(id);
    if (byId >= 0) return byId;
    return place < after.length && !kept.has(after[place]) ? place : -1;
  });
};

/**
 * The records a drawing's OFF records keep their rows under, as request
 * records a row is matched against (`rowBelongsToRequest`): a Linear card's
 * uid with its `<uid>:<k>` identities, or a Circular record's identity.
 * @param {Iterable<unknown> | null | undefined} recordsOff
 * @returns {{ recordKey: string, cardinality: 'all' }[]}
 */
export const offRecordRequestRecords = (recordsOff) => [...offSet(recordsOff)]
  .map((key) => ({ recordKey: key, cardinality: /** @type {const} */ ('all') }));

/**
 * A record of a source file as a record list names it.
 * @typedef {{ sourceIndex?: unknown, selector?: unknown, recordId?: unknown, record_id?: unknown,
 *   recordLength?: unknown }} RecordSourceEntry
 * The record lists the identities of a request's records are read from: each
 * card's records (a Linear record card's, a Circular source card's) by card
 * uid, in drawing order; each list in file order, as
 * `buildDisambiguatedRecordEntries` (services/record-options.js) gives it, so
 * `selector` and `recordId` are those a request key is built from.
 * `legacyOrder`: the uids of the cards whose records an older Session's flat
 * record positions count over, in their order then (`legacyRecordOrder`);
 * none when absent.
 * @typedef {{ mode: string, cards: ReadonlyMap<string, readonly RecordSourceEntry[]>,
 *   legacyOrder?: readonly string[] }} RecordIdentitySources
 * A record selector or region of a request record.
 * @typedef {{ kind?: string, index?: number, value?: string } | null | undefined} RequestRecordSelector
 * A record of a render request, as its identities are read from it.
 * @typedef {{ recordKey: string, cardinality?: string, selector?: RequestRecordSelector,
 *   windows: readonly unknown[] }} IdentityRequestRecord
 * One drawn record of a request: the record key its features have in that
 * request, its identity, the other keys an older Session may have left its rows
 * under (`legacyCircularRecordKeys`), the entry of its source list it is, and
 * the card it is a record of with the record IDs of that card's file in file
 * order (`cardRecordIds`), by which a Replace's records continue
 * (`continuedRecordPlaces`).
 * @typedef {{ requestKey: string, identity: string, aliases: readonly string[],
 *   entry: RecordSourceEntry, card: string, cardRecordIds: readonly string[] }} DrawnRecordIdentity
 */

/** The record ID of a source's record. @param {RecordSourceEntry} entry */
const entryRecordId = (entry) => cleanKey(entry?.recordId ?? entry?.record_id);

/**
 * The 0-based positions of a source's records, by record ID (a duplicated ID names none).
 * @param {readonly RecordSourceEntry[]} entries
 */
const positionsOf = (entries) => {
  const byId = new Map();
  entries.forEach((entry, index) => {
    const recordId = entryRecordId(entry);
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
export const selectedPosition = (record, entries) => {
  const selector = record.selector;
  if (selector?.kind === 'recordIndex' && Number.isInteger(selector.index)) {
    return Number(selector.index) < entries.length ? Number(selector.index) : null;
  }
  if (selector?.kind === 'recordId') return positionsOf(entries).get(cleanKey(selector.value)) ?? null;
  return entries.length === 1 ? 0 : null;
};

/**
 * The uids of the leading cards of `marked` that `inputs` (the drawing's
 * cards now, `circularSourceInputs`) still hold with the same files, in the
 * order of `marked`: the walk stops at the first card gone or whose file
 * changed, past which no position of then names a record of now.
 * @param {readonly CircularCardMark[]} marked
 * @param {readonly { uid: string, primaryFile: InputFileRef | null, pairedFile: InputFileRef | null }[]} inputs
 * @returns {string[]}
 */
export const heldCardPrefix = (marked, inputs) => {
  /** @type {string[]} */
  const held = [];
  for (const card of marked) {
    const now = inputs.find((input) => input.uid === card.uid);
    if (!now || now.primaryFile !== card.primaryFile || now.pairedFile !== card.pairedFile) break;
    held.push(card.uid);
  }
  return held;
};

/**
 * Whether a Session written before its Circular input was a list of File
 * cards (Session 50 and older) named its records by their flat position over
 * the cards it moved to (D-12): one card, or several whose uids are not the
 * CLI's per-file record keys (`record-<n>`), whose keys are already identities.
 * @param {readonly { uid: string }[]} cards
 */
export const cardsTakeLegacyRecordPositions = (cards) => cards.length === 1
  || (cards.length > 1 && !cards.every(({ uid }) => isCliRecordKey(uid)));

/**
 * The keys a Circular record's rows may wait under from a Session written
 * before its input was a list of File cards (Session 50 and older):
 * `record-<n>` and `circular-<record ID>-#<n>`, `n` its 1-based position in
 * the records of the cards that Session moved to, taken together in their
 * order at Load (`legacyOrder`), as the one input file of then held them.
 * The walk stops at a card whose records are not read. They are aliases
 * only: the first Generate moves their rows to the record's key
 * (`recordKeyMoves`), and no other reader takes them. A key that is the
 * identity of a record of the drawing is never an alias. Without a Load-time
 * order (a Session 51, a Session of the CLI's per-file keys) there are none.
 * @param {RecordIdentitySources} sources
 * @returns {Map<string, string[]>} each identity's legacy keys
 */
export const legacyCircularRecordKeys = (sources) => {
  const identities = new Set([...sources.cards].flatMap(([uid, entries]) => cardIdentities(uid, entries).map(({ identity }) => identity)));
  return new Map(legacyFlatRecords(sources).map(({ identity, entry }, position) => [identity, [
    cliRecordKey(position + 1),
    `circular-${downloadSafeName(entry.recordId ?? entry.record_id, 'record')}-${downloadSafeName(`#${position + 1}`, '1')}`
  ].filter((key) => !identities.has(key))]));
};

/** The identities of a card's records. @param {string} uid @param {readonly RecordSourceEntry[]} entries */
const cardIdentities = (uid, entries) => entries.map((entry, sourceIndex) => ({
  identity: recordIdentity({ uid, sourceIndex, recordCount: entries.length }), entry
}));

/**
 * The records an older Session's flat positions name, in that order: the
 * records of the cards as loaded (`legacyOrder`), up to the first card whose
 * records are not read. The one walk of those positions.
 * @param {RecordIdentitySources} sources
 * @returns {{ identity: string, entry: RecordSourceEntry }[]}
 */
const legacyFlatRecords = (sources) => {
  /** @type {{ identity: string, entry: RecordSourceEntry }[]} */
  const flat = [];
  for (const uid of sources.legacyOrder ?? []) {
    const entries = sources.cards.get(uid) ?? [];
    if (entries.length === 0) break;
    flat.push(...cardIdentities(uid, entries));
  }
  return flat;
};

const RECORD_INDEX_SELECTOR = new RegExp(`^(?:${RECORD_INDEX_SELECTOR_PATTERN})$`);
/**
 * The N of a `#N` record selector text as Python reads it (`record_index_number`,
 * gbdraw/io/record_select.py), 0 included, which a reader refuses; null for any
 * other value. A reader that names that refusal reads it; any other reads
 * `recordSelectorNumber`.
 * @param {unknown} key
 * @returns {number | null}
 */
export const recordIndexNumber = (key) => {
  const match = typeof key === 'string' ? RECORD_INDEX_SELECTOR.exec(key) : null;
  return match ? Number(match[1]) : null;
};
/**
 * The record, from 1, that a `#N` record selector names, as Python's renderer
 * reads it (`record_selector_number`, gbdraw/io/record_select.py, whose pattern
 * the generated contracts carry): spaces around it and after `#` and leading
 * zeros are allowed, so `#01` is record 1; null for a record ID or any other
 * value, and for `#0`, which the renderer refuses. The one Web reader of a
 * record selector's `#N` (OV-2960): the Session move, a kept `#N` (OV-2869),
 * a committed Record Order token, a Linear card's selector, an annotation
 * table's record column, and the sequence readers of LOSAT and of a match popup.
 * @param {unknown} key
 * @returns {number | null}
 */
export const recordSelectorNumber = (key) => {
  const number = recordIndexNumber(key);
  return number !== null && Number.isSafeInteger(number) && number >= 1 ? number : null;
};
const RECORD_ROW = new RegExp(`^(?:${RECORD_ROW_PATTERN})$`);
const MAX_SAFE_ROW = BigInt(Number.MAX_SAFE_INTEGER);
/**
 * The rows that row texts name, as Python reads a Record Order row (`record_row_number`, gbdraw/layout/record_placement.py,
 * whose pattern the generated contracts carry): spaces around it, an optional `+` and leading zeros are allowed, so
 * `+01` is row 1; at most 4300 digits, as `int()` converts; null for a row of zero (`0`, `+00`), any other text (`1.0`,
 * `0x1`, `1e0`) and a null text. The renderer reads only the order of the rows (`resolve_record_row_positions` ranks
 * them), so the rows of one list read as numbers in that exact order: each row's value when every row is a safe
 * integer, else each row's rank among the list's rows (OV-3313). The one Web reader of a row.
 * @param {readonly (string | null)[]} texts One list's row texts.
 * @returns {(number | null)[]}
 */
export const recordRowNumbers = (texts) => {
  const values = texts.map((text) => {
    const match = text === null ? null : RECORD_ROW.exec(text);
    const value = match ? BigInt(match[1]) : 0n;
    return value > 0n ? value : null;
  });
  /** @type {bigint[]} */
  const read = values.flatMap((value) => (value === null ? [] : [value]));
  if (read.every((value) => value <= MAX_SAFE_ROW)) return values.map((value) => (value === null ? null : Number(value)));
  const ranks = [...new Set(read)].sort((left, right) => (left < right ? -1 : left > right ? 1 : 0));
  return values.map((value) => (value === null ? null : ranks.indexOf(value) + 1));
};
/**
 * The rows of Record Order tokens `<selector>@<row>` (each row the text after the token's last `@`), as Python's
 * renderer reads them (`parse_record_row_position`, `recordRowNumbers`); null for a row the renderer refuses and for a
 * token that is no text or has no `@`.
 * @param {readonly unknown[]} tokens One list's tokens.
 * @returns {(number | null)[]}
 */
export const recordOrderRows = (tokens) => recordRowNumbers(tokens.map((token) => (
  typeof token === 'string' && token.includes('@') ? token.slice(token.lastIndexOf('@') + 1) : null
)));
/**
 * Whether a key is a record selector `#N` that a Session 46-50 Session move
 * kept as it was (OV-2869): the N-th record of all the drawing's records.
 * @param {unknown} key
 */
export const isFlatRecordSelector = (key) => recordSelectorNumber(key) !== null;

/**
 * The Circular record choices a Session names records in: the OFF records,
 * the Record select value, and the Record Order rows.
 * @template {{ selector: string }} P
 * @typedef {{ recordsOff: readonly string[], selector: string, positions: readonly P[] }} CircularRecordChoices
 */

/**
 * The choices with each flat `#N` (`isFlatRecordSelector`) placed at the
 * identity of the N-th record over the cards as loaded (`legacyFlatRecords`,
 * the walk the legacy keys take): a `#N` past them names no record now, so
 * an OFF record or a Record Order row is dropped and the Record select
 * cleared. Other keys stay.
 * @template {{ selector: string }} P
 * @param {CircularRecordChoices<P>} choices
 * @param {RecordIdentitySources} sources
 * @returns {CircularRecordChoices<P>}
 */
export const placeFlatRecordSelectors = ({ recordsOff, selector, positions }, sources) => {
  const flat = legacyFlatRecords(sources);
  /** @param {string} key @returns {string | null} */
  const place = (key) => {
    const number = recordSelectorNumber(key);
    return number === null ? key : flat[number - 1]?.identity ?? null;
  };
  return {
    recordsOff: recordsOff.flatMap((key) => place(key) ?? []),
    selector: place(selector) ?? '',
    positions: positions.flatMap((row) => {
      const placed = place(row.selector);
      return placed === null ? [] : [{ ...row, selector: placed }];
    })
  };
};

/**
 * The one record a request record names of its card's file, by its place or its record ID.
 * @typedef {{ kind: 'recordIndex', index: number } | { kind: 'recordId', value: string }} NamedRecordSelector
 * A record a Circular request names of its card's file: its request key, by which the draft names it (its Record
 * Order row, the request records before its card is read) until the card's read names it by its identity, and its
 * selector.
 * @typedef {{ recordKey: string, selector: NamedRecordSelector }} NamedCardRecord
 * The records a Circular File card draws of its file when a Load or a copy installed a request whose records of the
 * card each name one record (OV-3530: Python's copy of some records of a file; OV-3531: a CLI records table): the
 * drawing's settings hold them beside `recordsOff` (`recordsNamed`, never saved) until the card's read places them
 * (`placeNamedCardRecords`).
 * @typedef {{ uid: string, named: readonly NamedCardRecord[] }} NamedCardRecords
 */

/** @param {unknown} value @returns {value is NamedRecordSelector} */
const isNamedRecordSelector = (value) => {
  const selector = /** @type {{ kind?: unknown, index?: unknown, value?: unknown } | null} */ (
    typeof value === 'object' && value !== null ? value : null);
  return (selector?.kind === 'recordIndex' && Number.isSafeInteger(selector.index) && Number(selector.index) >= 0)
    || (selector?.kind === 'recordId' && typeof selector.value === 'string' && selector.value !== '');
};
/**
 * The record a request record's selector names of its card's file (`NamedRecordSelector`), or null when it names
 * none.
 * @param {unknown} selector
 * @returns {NamedRecordSelector | null}
 */
export const namedRecordSelectorOf = (selector) => {
  if (!isNamedRecordSelector(selector)) return null;
  return selector.kind === 'recordIndex' ? { kind: 'recordIndex', index: selector.index } : { kind: 'recordId', value: selector.value };
};
/**
 * A drawing's named card records as its settings or a History step hold them (`recordsNamed`); a value of another
 * shape names none.
 * @param {unknown} value
 * @returns {NamedCardRecords[]}
 */
export const readNamedCardRecords = (value) => (Array.isArray(value) ? value : []).flatMap((/** @type {unknown} */ card) => {
  const { uid, named } = /** @type {{ uid?: unknown, named?: unknown }} */ (typeof card === 'object' && card !== null ? card : {});
  if (typeof uid !== 'string' || uid === '' || !Array.isArray(named) || named.length === 0) return [];
  const records = named.flatMap((/** @type {unknown} */ record) => {
    const { recordKey, selector } = /** @type {{ recordKey?: unknown, selector?: unknown }} */ (
      typeof record === 'object' && record !== null ? record : {});
    const parsed = namedRecordSelectorOf(selector);
    return typeof recordKey === 'string' && recordKey !== '' && parsed ? [{ recordKey, selector: parsed }] : [];
  });
  return records.length === named.length ? [{ uid, named: records }] : [];
});

/**
 * What a drawing's named card records (`NamedCardRecords`) say against its read records: for each card whose
 * records are read and whose every named record is one of them, the identities of the file's other records (`off`)
 * and the identity of each named record by its request key (`identityOfKey`); `read` lists those cards. A card not
 * read, or one that names a record its file lacks, says nothing: it draws its named records alone through its
 * request records (`circularRecordPlan` in services/session-request.js). The one reading of the named records: the
 * records a Circular drawing draws (`circularRecordsOff`) and the placement (`placeNamedCardRecords`, which the
 * places a request names read too) read it.
 * @param {readonly NamedCardRecords[]} cards
 * @param {readonly (RecordSourceEntry & { uid?: unknown, value: string })[]} entries The drawing's read records,
 *   named (`circularRecordEntries`: `value` is a record's identity).
 * @returns {{ off: string[], identityOfKey: Map<string, string>, read: Set<string> }}
 */
const namedCardRecordReading = (cards, entries) => {
  /** @type {string[]} */
  const off = [];
  /** @type {Map<string, string>} */
  const identityOfKey = new Map();
  /** @type {Set<string>} */
  const read = new Set();
  cards.forEach(({ uid, named }) => {
    const cardEntries = entries.filter((entry) => cleanKey(entry.uid) === uid);
    const drawn = named.map(({ recordKey, selector }) => {
      const position = selectedPosition({ recordKey, selector, windows: [] }, cardEntries);
      return cardEntries.find((entry) => position !== null && Number(entry.sourceIndex) === position);
    });
    if (cardEntries.length === 0 || drawn.includes(undefined)) return;
    read.add(uid);
    named.forEach(({ recordKey }, index) => identityOfKey.set(recordKey, /** @type {{ value: string }} */ (drawn[index]).value));
    off.push(...cardEntries.filter((entry) => !drawn.includes(entry)).map((entry) => entry.value));
  });
  return { off, identityOfKey, read };
};

/**
 * The records a Circular drawing leaves OFF as its request draws them: its `recordsOff` and, of each read card that
 * names some records of its file (`recordsNamed`), the file's other records (`namedCardRecordReading`). The one
 * answer to "is this record drawn" for the request (`resolveCircularRequestRecordSet`), the placement, and every
 * control that shows it (the Records list, the Draw this record checkbox, the drawn count), so the screen shows what
 * the figure draws before and after the placement, also after an Undo across it.
 * @param {{ readonly recordsOff: readonly string[], readonly recordsNamed: readonly NamedCardRecords[] }} drawing
 * @param {readonly (RecordSourceEntry & { uid?: unknown, value: string })[]} entries The drawing's read records,
 *   named (`circularRecordEntries`).
 * @returns {string[]}
 */
export const circularRecordsOff = ({ recordsOff, recordsNamed }, entries) => {
  const { off } = namedCardRecordReading(recordsNamed, entries);
  return [...recordsOff, ...off.filter((key) => !recordsOff.includes(key))];
};

/**
 * The choices and named card records after a settled read of the drawing's cards places the named records
 * (`namedCardRecordReading`): the other records of each read card's file are OFF, a Record Order row that names a
 * named record by its request key names it by its identity (OV-3534), and the read cards' named records go
 * (`cards`: those left, a card not read or not placed). It changes no record the drawing draws, so a placement that
 * runs twice changes nothing more.
 * @template {{ selector: string }} P
 * @param {CircularRecordChoices<P>} choices
 * @param {readonly NamedCardRecords[]} cards
 * @param {readonly (RecordSourceEntry & { uid?: unknown, value: string })[]} entries The drawing's read records, named.
 * @returns {{ choices: CircularRecordChoices<P>, cards: NamedCardRecords[] }}
 */
export const placeNamedCardRecords = ({ recordsOff, selector, positions }, cards, entries) => {
  const { identityOfKey, read } = namedCardRecordReading(cards, entries);
  return {
    choices: {
      recordsOff: circularRecordsOff({ recordsOff, recordsNamed: cards }, entries), selector,
      positions: positions.map((row) => {
        const identity = identityOfKey.get(row.selector);
        return identity === undefined ? row : { ...row, selector: identity };
      })
    },
    cards: cards.filter(({ uid }) => !read.has(uid))
  };
};

/**
 * Places the records File cards name of their files (`named`) in a draft config, as a settled read of the cards
 * places them (`placeNamedCardRecords`): the other records of each card `entries` reads OFF, and the places that named
 * them by request key under their identities. An empty `recordsOff` is left out, as a Session writes it (omission =
 * every record ON). The one placement of a draft config: Save's, a copy's and a Gallery publication's.
 * @param {{ recordsOff?: unknown, adv?: unknown }} config Changed in place.
 * @param {readonly NamedCardRecords[]} named
 * @param {readonly (RecordSourceEntry & { uid?: unknown, value: string })[]} entries The cards' read records, named
 *   (`circularRecordEntries`).
 * @returns {NamedCardRecords[]} The cards it could not place: not read, or naming a record their file lacks.
 */
export const placeNamedRecordsInDraft = (config, named, entries) => {
  const adv = /** @type {{ multi_record_positions?: unknown }} */ (typeof config.adv === 'object' && config.adv !== null ? config.adv : {});
  /** @type {{ selector: string }[]} */
  const positions = Array.isArray(adv.multi_record_positions) ? adv.multi_record_positions : [];
  const placed = placeNamedCardRecords({
    recordsOff: Array.isArray(config.recordsOff) ? config.recordsOff.map(String) : [], selector: '', positions
  }, named, entries);
  if (placed.choices.recordsOff.length > 0) config.recordsOff = placed.choices.recordsOff;
  else delete config.recordsOff;
  if (Array.isArray(adv.multi_record_positions)) adv.multi_record_positions = placed.choices.positions;
  return placed.cards;
};

/**
 * The named records of the File cards a paste adds that its read could not place (`placeNamedCardRecords`), under
 * each card's uid in the target (`uids`, from the read request's): each named record keeps its request key under the
 * card, `<card uid>:<request key>`, which no key of the target holds (the card's uid is new), until the card's read
 * names it by its identity, as the request key of every named record does. `moves` takes each read request key to
 * its new key (the places of the read request). A card the paste does not add is left out.
 * @param {readonly NamedCardRecords[]} cards
 * @param {ReadonlyMap<string, string>} uids
 * @returns {{ cards: NamedCardRecords[], moves: Map<string, string> }}
 */
export const pastedNamedCardRecords = (cards, uids) => {
  /** @type {Map<string, string>} */
  const moves = new Map();
  const pasted = cards.flatMap(({ uid, named }) => {
    const card = uids.get(uid);
    if (card === undefined) return [];
    return [{ uid: card, named: named.map(({ recordKey, selector }) => {
      const key = `${card}:${recordKey}`;
      moves.set(recordKey, key);
      return { recordKey: key, selector };
    }) }];
  });
  return { cards: pasted, moves };
};

/**
 * The drawn records of a request with their identities (DESIGN 3, A1): each
 * record's position comes from the request's own selector, read against the
 * record lists of that request's sources, never from a card's current
 * selection. A record whose source list is unknown has no identity, so its
 * rows keep their key. A card that draws every record of a file of several
 * records draws each one under its identity. One rule for both modes; a
 * Circular record also brings the rows of its legacy keys.
 * @param {readonly IdentityRequestRecord[]} records
 * @param {RecordIdentitySources} sources
 * @returns {DrawnRecordIdentity[]}
 */
export const drawnRecordIdentities = (records, sources) => {
  const legacy = sources.mode === 'circular' ? legacyCircularRecordKeys(sources) : new Map();
  /** @param {string} identity */
  const aliasesOf = (identity) => [identity, ...(legacy.get(identity) || [])];
  return records.flatMap((record) => {
    const recordKey = cleanKey(record?.recordKey);
    const cards = sources.cards;
    const base = cards.has(recordKey) ? recordKey
      : [...cards.keys()].find((uid) => recordKeyExpansion(recordKey, uid) !== null);
    const entries = base === undefined ? undefined : cards.get(base);
    if (base === undefined || !entries || entries.length === 0) return [];
    const identityOf = (/** @type {number} */ sourceIndex) => recordIdentity({
      uid: base, sourceIndex, recordCount: entries.length
    });
    const card = { card: base, cardRecordIds: entries.map(entryRecordId) };
    if (record.cardinality === 'all' && !record.selector && record.windows.length === 0) {
      return entries.map((/** @type {RecordSourceEntry} */ entry, /** @type {number} */ sourceIndex) => {
        const identity = identityOf(sourceIndex);
        return { requestKey: identity, identity, aliases: aliasesOf(identity), entry, ...card };
      });
    }
    const expansion = recordKeyExpansion(recordKey, base);
    const position = expansion !== null ? expansion - 1 : selectedPosition(record, entries);
    if (position === null || position >= entries.length) return [];
    const identity = identityOf(position);
    return [{ requestKey: recordKey, identity, aliases: aliasesOf(identity), entry: entries[position], ...card }];
  });
};

/**
 * The records of the cards whose file's record IDs differ between two
 * requests (a Replace in the same card): `then`, each record's identity then
 * with the identity of the record it continues as now (null: it leaves), by
 * `continuedRecordPlaces`; `now`, the identities of those cards' records now.
 * @param {readonly DrawnRecordIdentity[]} previousDrawn
 * @param {readonly DrawnRecordIdentity[]} nextDrawn
 */
const continuedCardRecords = (previousDrawn, nextDrawn) => {
  const cardsOf = (/** @type {readonly DrawnRecordIdentity[]} */ drawn) => new Map(drawn.map(({ card, cardRecordIds }) => [card, cardRecordIds]));
  const nextCards = cardsOf(nextDrawn);
  /** @type {Map<string, string | null>} */
  const then = new Map();
  /** @type {Set<string>} */
  const now = new Set();
  cardsOf(previousDrawn).forEach((ids, uid) => {
    const after = nextCards.get(uid);
    if (!after || (after.length === ids.length && after.every((id, place) => id === ids[place]))) return;
    const places = continuedRecordPlaces(ids, after);
    ids.forEach((_id, place) => then.set(recordIdentity({ uid, sourceIndex: place, recordCount: ids.length }),
      places[place] < 0 ? null : recordIdentity({ uid, sourceIndex: places[place], recordCount: after.length })));
    after.forEach((_id, place) => now.add(recordIdentity({ uid, sourceIndex: place, recordCount: after.length })));
  });
  return { then, now };
};

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
 * from one snapshot, so a permutation is safe. A row never leaves or enters a
 * key whose source the Generate replaces (`unbound`): the replaced source no
 * longer binds that record (#494), and the prune decides its rows. The one
 * exception is a Replace in the same card (D-05, D-20): the records of a card
 * whose file's record IDs changed continue as `continuedRecordPlaces` says (by
 * record ID, else by place), and their rows move with them; a record that
 * leaves moves none.
 * @param {readonly DrawnRecordIdentity[]} previousDrawn
 * @param {readonly DrawnRecordIdentity[]} nextDrawn
 * @param {ReadonlySet<string>} [unbound]
 * @returns {Map<string, string>}
 */
export const recordKeyMoves = (previousDrawn, nextDrawn, unbound = new Set()) => {
  const previousKeyOf = new Map(previousDrawn.map(({ identity, requestKey }) => [identity, requestKey]));
  const nextKeyOf = new Map(nextDrawn.map(({ identity, requestKey }) => [identity, requestKey]));
  const previousKeys = new Set(previousKeyOf.values());
  const nextKeys = new Set(nextKeyOf.values());
  const continued = continuedCardRecords(previousDrawn, nextDrawn);
  // Each record with its identity then and now: a record of a changed card as it continues, any other as it is.
  /** @type {[string, string][]} */
  const records = [...continued.then].flatMap(([then, now]) => (now === null ? [] : [/** @type {[string, string]} */ ([then, now])]));
  new Set([...previousKeyOf.keys(), ...nextKeyOf.keys()]).forEach((identity) => {
    if (!continued.then.has(identity) && !continued.now.has(identity)) records.push([identity, identity]);
  });
  /** @type {Map<string, string>} */
  const moves = new Map();
  records.forEach(([then, now]) => {
    const from = previousKeyOf.get(then) ?? then;
    const to = nextKeyOf.get(now) ?? now;
    if (from === to || (!continued.then.has(then) && (unbound.has(from) || unbound.has(to)))) return;
    // A row under another record's request key belongs to that record.
    if (!previousKeyOf.has(then) && previousKeys.has(from)) return;
    if (!nextKeyOf.has(now) && nextKeys.has(to)) return;
    moves.set(from, to);
  });
  const movedInto = new Set(moves.values());
  nextDrawn.forEach(({ requestKey, aliases }) => aliases.forEach((alias) => {
    if (alias === requestKey || moves.has(alias) || nextKeys.has(alias) || previousKeys.has(alias)
      || unbound.has(alias) || unbound.has(requestKey)) return;
    // The rows already under the record's key, which no move brings there,
    // are its own: they come first (`requestKey` -> `requestKey`).
    if (!movedInto.has(requestKey) && !moves.has(requestKey)) moves.set(requestKey, requestKey);
    moves.set(alias, requestKey);
  }));
  return moves;
};

/**
 * The rows a list keeps one per record (a record-major Depth matrix) after the
 * records of a drawing change (a source card added, removed or moved): each
 * record of `next` takes the row its identity had in `previous`, and a record
 * new to the list an empty row. It reads records by identity, never by place,
 * from the same lists `recordKeyMoves` moves rows between. Rows past the end of
 * a shorter list are the records that have none yet.
 * @template R
 * @param {readonly R[]} rows one per identity of `previous`, in its order
 * @param {readonly (string | null)[]} previous null: a record whose row no record of `next` takes
 * @param {readonly string[]} next
 * @param {() => R} empty
 * @returns {R[]}
 */
export const resliceRecordRows = (rows, previous, next, empty) => next.map((identity) => {
  const index = previous.indexOf(identity);
  return index >= 0 && index < rows.length ? rows[index] : empty();
});

/**
 * A File card at one moment with the record IDs of its file, in file order.
 * @typedef {{ uid: string, recordIds: readonly string[] }} CardRecords
 */

/**
 * The rows of a list that keeps one row per record of a drawing's cards (a
 * record-major Depth matrix) after its cards change: `marked` are the cards
 * the rows stand for, in row order, `now` the drawing's cards. A record keeps
 * the row of the record of its card it continues (`continuedRecordPlaces`, the
 * rule a Generate moves the record's other rows by), so a card moved or
 * removed takes its rows along, and a Replace in the same card keeps each
 * record's row with the record of the new file with its record ID, or at its
 * place (D-05, D-20). A record new to the list gets an empty row (`resliceRecordRows`).
 * @template R
 * @param {readonly R[]} rows
 * @param {readonly CardRecords[]} marked
 * @param {readonly CardRecords[]} now
 * @param {() => R} empty
 * @returns {R[]}
 */
export const resliceCardRows = (rows, marked, now, empty) => resliceRecordRows(
  rows, continuedCardRows(marked, now).map(({ key }) => key), cardRecordKeys(now), empty
);

/**
 * Each record of `cards`, in row order: a record is its card's record at its
 * place in the file (`k`, its key `<uid>:<k>`).
 * @param {readonly CardRecords[]} cards
 * @returns {{ uid: string, k: number, key: string }[]}
 */
const cardRecords = (cards) => cards.flatMap(({ uid, recordIds }) => recordIds.map(
  (_id, index) => ({ uid, k: index + 1, key: expandedRecordKey(uid, index + 1) })
));
/** @param {readonly CardRecords[]} cards */
const cardRecordKeys = (cards) => cardRecords(cards).map(({ key }) => key);

/**
 * Each record of `marked`, in row order, with the key (`cardRecords`) of the
 * record of `now` it continues as, null when it leaves (its card removed, or
 * a Replace whose new file has no record for it). The one mapping of a card
 * edit's rows (`resliceCardRows`, `cardRowsLeaving`).
 * @param {readonly CardRecords[]} marked
 * @param {readonly CardRecords[]} now
 * @returns {{ uid: string, k: number, key: string | null }[]}
 */
const continuedCardRows = (marked, now) => marked.flatMap(({ uid, recordIds }) => {
  const after = now.find((card) => card.uid === uid);
  const places = after ? continuedRecordPlaces(recordIds, after.recordIds) : recordIds.map(() => -1);
  return recordIds.map((_id, index) => ({
    uid, k: index + 1, key: places[index] < 0 ? null : expandedRecordKey(uid, places[index] + 1)
  }));
});

/**
 * The rows `resliceCardRows` leaves out when the cards change from `marked`
 * to `now`: those of the records that continue as no record of `now`, with
 * the record each stands for.
 * @template R
 * @param {readonly R[]} rows one per record of `marked`, in its order
 * @param {readonly CardRecords[]} marked
 * @param {readonly CardRecords[]} now
 * @returns {{ uid: string, k: number, row: R }[]}
 */
export const cardRowsLeaving = (rows, marked, now) => continuedCardRows(marked, now)
  .flatMap(({ uid, k, key }, index) => (key !== null || index >= rows.length ? [] : [{ uid, k, row: rows[index] }]));

/**
 * Whether a key has the shape of a record draw key: a card uid or a record
 * identity, the same shape in both modes.
 * @param {'circular' | 'linear'} _mode
 * @param {unknown} key
 */
export const isRecordDrawKey = (_mode, key) => typeof key === 'string' && key.trim() === key && key !== '';

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

/**
 * A Linear card's window fields when it draws its whole record: no window, forward. The one owner of that value for
 * every write that clears a card's windows (Reset Settings, Delete settings, a replaced source file); the edits of a
 * card's windows are services/linear-card-windows.js. Each call returns new objects.
 * @returns {{ region_start: null, region_end: null, region_reverse: false, windows: never[] }}
 */
export const wholeRecordCardWindows = () => ({ region_start: null, region_end: null, region_reverse: false, windows: [] });

// D-08 Delete settings: what a record keeps while it is OFF and what Delete
// settings removes. The card fields return to these values (the record
// fields Reset Settings resets, plus LOSAT Gencode).
export const RECORD_SETTINGS_DEFAULTS = Object.freeze({
  definition: '', record_subtitle: '', ...wholeRecordCardWindows(), losat_gencode: 1
});

/**
 * Whether a Linear record card holds a value that Delete settings resets.
 * @param {Record<string, unknown> | null | undefined} card
 */
export const linearCardHasSettings = (card) => Object.entries(RECORD_SETTINGS_DEFAULTS)
  .some(([field, fallback]) => {
    const value = card?.[field] ?? fallback;
    // A window list holds a setting when it lists a window (OV-2993: `[] !== []` held for every card).
    return Array.isArray(fallback) ? Array.isArray(value) && value.length > 0 : value !== fallback;
  });

/**
 * One record whose edits Delete settings finds or removes.
 * @typedef {object} RecordSettingsOwner
 * @property {RecordDrawKey} key
 * @property {readonly string[]} requestKeys The request record keys its feature edits use: a Linear
 *   card's uid (which also owns its `<uid>:<n>` records), a Circular record's identity and its
 *   legacy keys (`legacyCircularRecordKeys`).
 * @property {boolean} [ownsExpansions] Linear: the key also owns `<key>:<n>`.
 * @property {readonly string[]} bindingKeys The record-catalog keys an annotation binds it by.
 * @property {string} displaySource The `sourceUid` of its record display drafts: its card's uid.
 * @property {string | null} displaySelector Circular: its `#k` in its card's file; Linear: null (every draft of the card).
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
 * @property {readonly { readonly recordKey: string, readonly remove: () => void }[]} [recordRangeEdits] The crops a
 *   Circular grid or batch keeps for its records, by the record's identity; `remove` drops the crop.
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
 * The record-catalog key an annotation binds the request record `recordKey` by (a Result's record key): the binding
 * key of the record that draws it. A Linear card that draws every record of its file owns `<uid>:<k>`, the k-th of
 * its records in the file; its binding keys follow the file's order. '' when no record draws it.
 * @param {readonly RecordSettingsOwner[]} owners
 * @param {string} recordKey
 * @returns {string}
 */
export const bindingKeyOfRequestRecord = (owners, recordKey) => {
  for (const owner of owners) {
    if (owner.requestKeys.includes(recordKey)) return owner.bindingKeys.length === 1 ? owner.bindingKeys[0] : '';
    if (!owner.ownsExpansions) continue;
    const k = owner.requestKeys.map((key) => recordKeyExpansion(recordKey, key)).find((value) => value !== null);
    if (k != null) return owner.bindingKeys[k - 1] ?? '';
  }
  return '';
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
  [...(data.colorRuleEdits || []), ...(data.recordRangeEdits || [])].forEach(({ recordKey, remove }) => {
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
 * annotation, record display draft, grid or batch crop, or comparison pair in the drawing (D-08).
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
