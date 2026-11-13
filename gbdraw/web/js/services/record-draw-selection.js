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
 * under (`legacyCircularRecordKeys`), and the entry of its source list it is.
 * @typedef {{ requestKey: string, identity: string, aliases: readonly string[],
 *   entry: RecordSourceEntry }} DrawnRecordIdentity
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
    if (record.cardinality === 'all' && !record.selector && record.windows.length === 0) {
      return entries.map((/** @type {RecordSourceEntry} */ entry, /** @type {number} */ sourceIndex) => {
        const identity = identityOf(sourceIndex);
        return { requestKey: identity, identity, aliases: aliasesOf(identity), entry };
      });
    }
    const expansion = recordKeyExpansion(recordKey, base);
    const position = expansion !== null ? expansion - 1 : selectedPosition(record, entries);
    if (position === null || position >= entries.length) return [];
    const identity = identityOf(position);
    return [{ requestKey: recordKey, identity, aliases: aliasesOf(identity), entry: entries[position] }];
  });
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
 * exception is the first record a Replace in the same card keeps (D-20).
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
  // D-20: a Replace in the same card keeps its first record, so the record of
  // a one-record file (`<uid>`) and the first record of a file of several
  // (`<uid>:1`) are one record when each request has only one of them. Its
  // rows move although its source was replaced; the prune keeps the edits of
  // the features the new file still has.
  /** @type {Set<string>} */
  const replacedFirstRecords = new Set();
  [...nextKeyOf.keys()].filter((identity) => !previousKeyOf.has(identity)).forEach((identity) => {
    const other = [...previousKeyOf.keys()].find((candidate) => !nextKeyOf.has(candidate)
      && (recordKeyExpansion(identity, candidate) === 1 || recordKeyExpansion(candidate, identity) === 1));
    if (other === undefined) return;
    previousKeyOf.set(identity, /** @type {string} */ (previousKeyOf.get(other)));
    previousKeyOf.delete(other);
    replacedFirstRecords.add(identity);
  });
  /** @type {Map<string, string>} */
  const moves = new Map();
  new Set([...previousKeyOf.keys(), ...nextKeyOf.keys()]).forEach((identity) => {
    const from = previousKeyOf.get(identity) ?? identity;
    const to = nextKeyOf.get(identity) ?? identity;
    if (from === to || (!replacedFirstRecords.has(identity) && (unbound.has(from) || unbound.has(to)))) return;
    // A row under another record's request key belongs to that record.
    if (!previousKeyOf.has(identity) && previousKeys.has(from)) return;
    if (!nextKeyOf.has(identity) && nextKeys.has(to)) return;
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
 * @param {readonly string[]} previous
 * @param {readonly string[]} next
 * @param {() => R} empty
 * @returns {R[]}
 */
export const resliceRecordRows = (rows, previous, next, empty) => next.map((identity) => {
  const index = previous.indexOf(identity);
  return index >= 0 && index < rows.length ? rows[index] : empty();
});

/**
 * The rows of a list that keeps one row per record of a drawing's cards (a
 * record-major Depth matrix) after its cards change: `marked` are the cards
 * the rows stand for, in row order, `now` the drawing's cards. A record keeps
 * the row of its card's record at the same place in the file (`<uid>:<k>`),
 * so a card moved or removed takes its rows along, and a Replace in the same
 * card keeps the rows of the records the new file has at those places (D-20;
 * the record of a one-record file is the first record of a file of several).
 * A record new to the list gets an empty row (`resliceRecordRows`).
 * @template R
 * @param {readonly R[]} rows
 * @param {readonly { uid: string, recordCount: number }[]} marked
 * @param {readonly { uid: string, recordCount: number }[]} now
 * @param {() => R} empty
 * @returns {R[]}
 */
export const resliceCardRows = (rows, marked, now, empty) => {
  /** @param {readonly { uid: string, recordCount: number }[]} cards */
  const keysOf = (cards) => cards.flatMap(({ uid, recordCount }) => Array.from({ length: recordCount }, (_unused, index) => expandedRecordKey(uid, index + 1)));
  return resliceRecordRows(rows, keysOf(marked), keysOf(now), empty);
};

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
