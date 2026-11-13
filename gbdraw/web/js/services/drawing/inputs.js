// @ts-check
// The input bindings: which stored file each input field names. A field holds
// an `InputFileRef` of the project's file store (`state.inputFiles`), never a
// File. The writers here take the object they write (the project's bindings
// now; a drawing's inputs slice later), so the store registers every bound
// file through one owner. The bound IDs keep their files stored: the store
// releases a file only when no History step and no binding names it.

/** @import { InputFileCandidate, InputFileRef, InputFileStore } from '../history-files.js' */
/** @import { RequestCanonicalComparison } from '../session-request.js' */
/** @import { LinearComparisonEdge, LinearComparisonPlan } from '../linear-comparisons.js' */
/** @import { DepthFileSlots } from '../depth-track-state.js' */
/** @import { CircularRecordDiscovery } from '../record-discovery.js' */
/** @import { FieldReader } from '../../utils/field-reading.js' */
/** @import { SliceHandle } from './slice-handle.js' */
import { CIRCULAR_SOURCE_FILE_SLOTS } from '../../contracts/schema.generated.js';
import { isInputFileRef } from '../history-files.js';
import { diagnosticError } from '../../utils/error-normalization.js';
import { reactive, ref } from '../../utils/reactivity.js';
import { isPositiveInteger, numberOr } from '../../utils/field-reading.js';
import { createSliceHandle, sliceMembers } from './slice-handle.js';

// The Linear record rows (`LinearSeq`): one owner builds and normalizes them.
let linearSeqUidCounter = 0;

const generateLinearSeqUid = () => {
  linearSeqUidCounter += 1;
  if (globalThis.crypto?.randomUUID) {
    return `linear-seq-${globalThis.crypto.randomUUID()}`;
  }
  return `linear-seq-${Date.now()}-${linearSeqUidCounter}`;
};

/** @param {unknown} value */
const normalizeLinearSeqNumber = (value) => {
  if (value === null || value === undefined || value === '') return null;
  const numeric = Number(value);
  return Number.isFinite(numeric) ? numeric : null;
};

/**
 * One linear record slot as `createLinearSeq` builds it. The file fields hold
 * `InputFileRef`s of `state.inputFiles`, or null.
 * @typedef {object} LinearSeq
 * @property {string} uid
 * @property {InputFileRef | null} gb
 * @property {InputFileRef | null} gff
 * @property {InputFileRef | null} fasta
 * @property {DepthFileSlots<InputFileRef>} depth one file, or a list indexed by series
 * @property {number | null} losat_gencode the card's LOSAT translation table; null when the field is cleared
 *   (`v-blank-null`), which the search passes as no table (`searchedLosatGencode`)
 * @property {string} definition
 * @property {string} record_subtitle
 * @property {string} file_definition
 * @property {string} file_subtitle
 * @property {string} inferred_definition
 * @property {string} region_record_id
 * @property {unknown} region_start a draft: `createLinearSeq` and the crop control write a number or null,
 *   a Session load stores the saved value as it is; the request projection validates it (R7)
 * @property {unknown} region_end a draft, as `region_start`
 * @property {boolean} region_reverse
 * @property {LinearSeqWindow[]} windows the windows the card's crop fields cannot hold (`cardCropWindow`): those of a
 *   record drawn in several (D-01), in the order the request gave them or the card's Windows list set, or one
 *   window across the origin of a circular record (D-34); empty for a record drawn whole or in one window that
 *   `region_*` hold. The card edits them through `services/linear-card-windows.js` (PD-OI-124).
 */

/**
 * One window of a record drawn in several: its span in source coordinates (1-based, both ends included; both
 * empty draws the whole record, PD-OI-124), its orientation, its column in the record's row (D-08), and
 * `acrossOrigin`, the card's "Across the origin" box: ticked, a Start after the End runs across the origin of a
 * circular record (D-45). The request validates it (R7).
 * @typedef {object} LinearSeqWindow
 * @property {number | null} start
 * @property {number | null} end
 * @property {boolean} reverseComplement
 * @property {number | null} column
 * @property {boolean} acrossOrigin
 */

/**
 * The windows a slot, the Circular form or a Session stores, read as `LinearSeqWindow`s; a value that is no list reads
 * as none. A window without `acrossOrigin` (a request's, or a slot saved before the box) crosses the origin when its
 * start is after its end: only a producer writes such a window (D-45).
 * @param {unknown} value
 * @returns {LinearSeqWindow[]}
 */
export const readLinearSeqWindows = (value) => (Array.isArray(value) ? value : []).map((item) => {
  const window = /** @type {Partial<Record<keyof LinearSeqWindow, unknown>>} */ (
    item && typeof item === 'object' && !Array.isArray(item) ? item : {}
  );
  const start = normalizeLinearSeqNumber(window.start);
  const end = normalizeLinearSeqNumber(window.end);
  return {
    start,
    end,
    reverseComplement: Boolean(window.reverseComplement),
    column: normalizeLinearSeqNumber(window.column),
    acrossOrigin: Object.hasOwn(window, 'acrossOrigin')
      ? Boolean(window.acrossOrigin) : start !== null && end !== null && start > end
  };
});

/**
 * The window a card's own crop fields (`region_*`) hold: the only window of `windows`, unless it crosses the origin
 * of a circular record (its "Across the origin" box, or, for a request's window, its start after its end, D-34),
 * which those fields would read as a reversed span (an End before the Start reverses the crop). Null when `windows`
 * stays the card's list. It reads the card's windows and a request record's, whose endpoints are unvalidated there.
 * @template {{ start: unknown, end: unknown, acrossOrigin?: boolean }} W
 * @param {readonly W[]} windows
 * @returns {W | null}
 */
export const cardCropWindow = (windows) => {
  const [only] = windows.length === 1 ? windows : [null];
  const crosses = only && (typeof only.acrossOrigin === 'boolean' ? only.acrossOrigin
    : typeof only.start === 'number' && typeof only.end === 'number' && only.start > only.end);
  return only && !crosses ? only : null;
};

/**
 * A Linear card's genetic code as stored (a draft, a History step, a Session):
 * null is a cleared field (`v-blank-null`), which stays cleared (OV-3463); a
 * whole number from 1 is that table; anything else, a card without one
 * included, is 1. A '' is a cleared field as gbdraw 0.14.5 saved it (Session
 * 50), which always loaded as 1 and still does.
 * @type {FieldReader<number | null>}
 */
export const readCardGencode = (value) => (value === null ? { value: null, lost: false } : numberOr(isPositiveInteger, 1)(value));

/**
 * @param {Partial<LinearSeq> | Record<string, unknown> | null} [overrides] a slot,
 *   or unvalidated data (a session, an import) that is normalized here
 * @returns {LinearSeq}
 */
export const createLinearSeq = (overrides = {}) => {
  const source = /** @type {Partial<LinearSeq>} */ (
    overrides && typeof overrides === 'object' && !Array.isArray(overrides) ? overrides : {}
  );
  const rawUid = String(source.uid ?? '').trim();
  // One window is the card's own crop, which `region_*` hold and the card edits; `windows` holds the rest.
  const stored = readLinearSeqWindows(source.windows);
  const onlyWindow = cardCropWindow(stored);
  const windows = onlyWindow ? [] : stored;
  return {
    uid: rawUid || generateLinearSeqUid(),
    gb: source.gb ?? null,
    gff: source.gff ?? null,
    fasta: source.fasta ?? null,
    depth: source.depth ?? null,
    losat_gencode: readCardGencode(source.losat_gencode).value,
    definition: String(source.definition ?? ''),
    record_subtitle: String(source.record_subtitle ?? ''),
    file_definition: String(source.file_definition ?? ''),
    file_subtitle: String(source.file_subtitle ?? ''),
    inferred_definition: String(source.inferred_definition ?? ''),
    region_record_id: String(source.region_record_id ?? ''),
    region_start: onlyWindow ? onlyWindow.start : windows.length ? null : normalizeLinearSeqNumber(source.region_start),
    region_end: onlyWindow ? onlyWindow.end : windows.length ? null : normalizeLinearSeqNumber(source.region_end),
    region_reverse: onlyWindow ? onlyWindow.reverseComplement : windows.length ? false : Boolean(source.region_reverse),
    windows
  };
};

/**
 * @param {readonly (Partial<LinearSeq> | Record<string, unknown> | null)[] | null | undefined} items
 *   slots, or unvalidated data (a session, an import), as `createLinearSeq` takes them
 * @returns {LinearSeq[]}
 */
export const normalizeLinearSeqList = (items) => {
  const baseItems = Array.isArray(items) && items.length > 0 ? items : [null];
  const seenUids = new Set();
  const normalized = baseItems.map((item) => {
    const next = createLinearSeq(item);
    if (!next.uid || seenUids.has(next.uid)) {
      next.uid = generateLinearSeqUid();
    }
    seenUids.add(next.uid);
    return next;
  });
  return normalized;
};

/** @param {LinearSeq | null | undefined} seq */
const hasLinearSeqPrimaryInput = (seq) => Boolean(seq?.gb || seq?.gff || seq?.fasta);

/**
 * @param {readonly (Partial<LinearSeq> | Record<string, unknown> | null)[] | null | undefined} items
 * @returns {LinearSeq[]}
 */
export const collapseEmptyLinearSeqList = (items) => {
  const previous = normalizeLinearSeqList(items);
  if (previous.length <= 1) return previous;
  const collapsed = previous.filter((seq) => hasLinearSeqPrimaryInput(seq));
  if (collapsed.length === previous.length) return previous;
  return normalizeLinearSeqList(collapsed);
};


// The Circular input sources (`CircularSource`): one card per file, in drawing
// order, as the CLI takes one record input per `--gbk` (or `--gff`/`--fasta`)
// value. The card's uid is its record key (a CLI Session's `record-<n>`, a
// random one for a card made in the Web), so a record's identity is `uid:k`,
// or `uid` for a file of one record, as a Linear card's (D-08).
let circularSourceUidCounter = 0;

/**
 * The uid of a Circular File card made in the Web, which keys its records:
 * `circular-source-<random UUID>`, or `circular-source-<time>-<n>` where no
 * random UUID is available.
 */
const CIRCULAR_SOURCE_UID_PREFIX = 'circular-source-';
/** Matches each uid `generateCircularSourceUid` composes. */
export const CIRCULAR_SOURCE_UID = Object.freeze(new RegExp(
  `^${CIRCULAR_SOURCE_UID_PREFIX}(?:[0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12}|[0-9]+-[0-9]+)$`));

const generateCircularSourceUid = () => {
  circularSourceUidCounter += 1;
  return `${CIRCULAR_SOURCE_UID_PREFIX}${globalThis.crypto?.randomUUID
    ? globalThis.crypto.randomUUID() : `${Date.now()}-${circularSourceUidCounter}`}`;
};

/**
 * A file slot of a Circular File card (`CIRCULAR_SOURCE_FILE_SLOTS`, Python's `CircularSourceCard`).
 * @typedef {(typeof CIRCULAR_SOURCE_FILE_SLOTS)[number]} CircularSourceFileSlot
 */
/**
 * One Circular input source: the drawing's input type picks the files it
 * reads (`gb`, or `gff` with `fasta`, or `fasta`); the others stay as they were.
 * @typedef {{ uid: string } & Record<CircularSourceFileSlot, InputFileRef | null>} CircularSource
 */

/**
 * @param {Partial<CircularSource>} [overrides]
 * @returns {CircularSource}
 */
export const createCircularSource = (overrides = {}) => ({
  uid: String(overrides.uid ?? '').trim() || generateCircularSourceUid(),
  gb: overrides.gb ?? null,
  gff: overrides.gff ?? null,
  fasta: overrides.fasta ?? null
});

/**
 * A drawing's bindings, as the retention walk reads them: its uploaded
 * comparison files and, from the inputs slice, its input fields and record
 * rows.
 * @typedef {object} BoundDrawingInputs
 * @property {{ readonly edges?: readonly { readonly id: string, readonly file?: unknown }[] } | null} [linearComparisonPlan]
 * @property {Readonly<Record<string, unknown>>} files
 * @property {readonly (Readonly<Record<string, unknown>> & { readonly uid: string })[]} linearSeqs
 * @property {readonly Readonly<CircularSource>[]} circularSources
 */
/**
 * The bindings that name stored files: every drawing's bindings, its five
 * card files among its input fields.
 * @typedef {object} InputBindings
 * @property {() => readonly BoundDrawingInputs[]} drawingList
 */


/**
 * A Linear comparison a Session committed, as the bindings hold it: the
 * request entry's fields, with a resource-backed entry's file by ref.
 * @typedef {Omit<RequestCanonicalComparison, 'file'> & { file?: InputFileRef | null }} BoundCanonicalComparison
 */
/**
 * A binding's value: one ref, or a list of them (a Depth matrix nests lists).
 * @typedef {InputFileRef | null | undefined | readonly BoundFileValue[]} BoundFileValue
 */

/**
 * The input file fields of a drawing (its inputs slice's `files`): each names
 * its stored files by `InputFileRef` of `state.inputFiles`.
 * @typedef {object} DrawingInputFiles
 * @property {InputFileRef | (InputFileRef | null)[] | (InputFileRef | null)[][] | null} c_depth one record's series, or a record-major matrix
 * @property {InputFileRef[]} c_conservation_blasts
 * @property {string | null} c_conservation_blasts_source a flag (`'losat-cache'`), not a file
 * @property {InputFileRef[][]} c_conservation_fastas the files of each LOSAT ring's comparison genome, one list per
 *   ring in source order (a LOSAT-cache replay keeps a ring without files as an empty list beside its table)
 * @property {InputFileRef[][]} c_conservation_sequence_sources the comparison sequence files of each ring whose
 *   table was uploaded or replayed (the span export reads them), one list per table place
 * @property {BoundCanonicalComparison[]} linearCanonicalComparisons
 * @property {InputFileRef | null} d_color the Default colors card's file
 * @property {InputFileRef | null} t_color the Specific rules card's file
 * @property {InputFileRef | null} blacklist the Blacklist card's file
 * @property {InputFileRef | null} whitelist the Whitelist card's file
 * @property {InputFileRef | null} qualifier_priority the Priority card's file
 * @property {readonly CircularCardMark[] | null} c_depth_rows the cards whose records the rows of a
 *   record-major `c_depth` stand for, in row order (not saved; History keeps it with the matrix);
 *   null: the rows stand for the records as the drawing's cards list them when they are first read
 * @property {readonly CircularCardMark[] | null} c_legacy_record_cards the cards a Session 50 or older
 *   moved its Circular input to, as Load installed them (`cardsTakeLegacyRecordPositions`; not saved);
 *   null: no record of the drawing has a flat position of then
 */
/**
 * A Circular File card as it was at one moment (a Load, the write of a Depth
 * matrix): its uid, the files its input type reads, and its record count when
 * that was known.
 * @typedef {{ uid: string, primaryFile: InputFileRef | null, pairedFile: InputFileRef | null,
 *   recordCount: number | null }} CircularCardMark
 */
/**
 * The five card files' fields (Default colors, Specific rules, Blacklist,
 * Whitelist, Priority): each drawing binds its own (R16-4).
 * @typedef {'d_color' | 't_color' | 'blacklist' | 'whitelist' | 'qualifier_priority'} CardFileSlot
 */
/** @type {readonly CardFileSlot[]} */
const CARD_FILE_SLOTS = Object.freeze(['d_color', 't_color', 'blacklist', 'whitelist', 'qualifier_priority']);

/** The `DrawingInputFiles` fields that hold a list of refs. @typedef {'c_conservation_blasts'} InputFileListSlot */

/**
 * The files of each ring's comparison genome, one list per place, an empty
 * place kept as an empty list: a drawing holds its ring genomes so
 * (`c_conservation_fastas`, `c_conservation_sequence_sources`). A place that
 * is one file or null is the form a Session wrote before a ring took several
 * files (Session 50 and older): the genome of that one file, or none.
 * Resource IDs read from JSON are read by `ringGenomeIdLists`.
 * @template T
 * @param {readonly RingGenomePlace<T>[] | null | undefined} places
 * @returns {T[][]}
 */
export const ringGenomeFileLists = (places) => (places || []).map((place) => (
  isRingGenomeFileList(place)
    ? place.filter(/** @returns {file is T} */ (file) => Boolean(file))
    : (place ? [place] : [])
));
/**
 * A place of a ring genome list: one file (Session 50 and older), a ring's files, or none.
 * @template T
 * @typedef {T | readonly (T | null | undefined)[] | null | undefined} RingGenomePlace
 */
/**
 * @template T
 * @param {RingGenomePlace<T>} place
 * @returns {place is readonly (T | null | undefined)[]}
 */
const isRingGenomeFileList = (place) => Array.isArray(place);
/**
 * The resource IDs of a ring genome list read from JSON (a Session's
 * `webMetadata`), one list per ring, in the forms `ringGenomeFileLists`
 * reads; what is not a resource ID names no file, as before.
 * @param {unknown} value
 * @returns {string[][]}
 */
export const ringGenomeIdLists = (value) => (Array.isArray(value) ? value : []).map((/** @type {unknown} */ place) => (
  (Array.isArray(place) ? place : [place]).filter(/** @returns {id is string} */ (id) => typeof id === 'string' && id !== '')
));

/**
 * The refs of a file or a list of files, each place kept (null for no file).
 * @param {Pick<InputFileStore, 'registerFile'>} store
 * @param {InputFileCandidate | readonly InputFileCandidate[]} list
 * @returns {(InputFileRef | null)[]}
 */
export const registerInputFiles = (store, list) => (
  (Array.isArray(list) ? list : (list ? [list] : [])).map((file) => store.registerFile(file))
);

// The fields of `files` that bind input files; the other members are flags
// (`c_conservation_blasts_source`) or generated data.
const BOUND_FILE_KEYS = Object.freeze([
  'c_depth', 'c_conservation_blasts', 'c_conservation_fastas',
  'c_conservation_sequence_sources', 'd_color', 't_color', 'blacklist', 'whitelist', 'qualifier_priority'
]);
const LINEAR_SEQ_FILE_KEYS = Object.freeze(['gb', 'gff', 'fasta', 'depth']);
// The fields of `files` that hold one list of genome files per ring (`boundRingGenomeLists`).
const RING_GENOME_KEYS = new Set(['c_conservation_fastas', 'c_conservation_sequence_sources']);
/**
 * A bound ring genome field as one list of refs per ring (`ringGenomeFileLists`):
 * a place nested deeper than a ring's files names no file.
 * @param {BoundFileValue} bound
 * @returns {InputFileRef[][]}
 */
const boundRingGenomeLists = (bound) => (isBoundFileList(bound) ? bound : []).map((place) => (
  isBoundFileList(place) ? place : [place]
).filter(/** @returns {file is InputFileRef} */ (file) => Boolean(file) && !isBoundFileList(file)));
/**
 * @param {BoundFileValue} value
 * @returns {value is readonly BoundFileValue[]}
 */
const isBoundFileList = (value) => Array.isArray(value);
// The fields of `files` that hold card marks (`CircularCardMark`), whose files bind no input.
const CARD_MARK_KEYS = Object.freeze(['c_depth_rows', 'c_legacy_record_cards']);

/**
 * The kinds of file an uploader takes; a binding slot's kind is the kind of
 * the files it binds. Each of the five cards is its own kind.
 * @typedef {'genbank' | 'gff' | 'fasta' | 'depth' | 'comparison' | CardFileSlot} InputFileKind
 */
/**
 * The slots an input file replacement writes (PD-OI-090 〔入力ファイル〕): every
 * uploader slot except the five cards'. `row.` names a Linear record row's field,
 * `source.` a Circular source card's, `edge` a Linear comparison's uploaded file.
 * @typedef {'c_depth' | 'c_conservation_blasts' | 'source.gb' | 'source.gff' | 'source.fasta'
 *   | 'row.gb' | 'row.gff' | 'row.fasta' | 'row.depth' | 'edge'} ReplaceableInputSlot
 */
/**
 * The kind of each slot that an uploader binds; a bound slot missing here
 * (the conservation sequences, a committed comparison) is offered to no
 * uploader.
 * @type {Readonly<Record<ReplaceableInputSlot | CardFileSlot, InputFileKind>>}
 */
const FILE_KIND_OF_SLOT = Object.freeze({
  'source.gb': 'genbank', 'row.gb': 'genbank', 'source.gff': 'gff', 'row.gff': 'gff',
  'source.fasta': 'fasta', 'row.fasta': 'fasta',
  c_depth: 'depth', 'row.depth': 'depth', c_conservation_blasts: 'comparison', edge: 'comparison',
  d_color: 'd_color', t_color: 't_color', blacklist: 'blacklist', whitelist: 'whitelist',
  qualifier_priority: 'qualifier_priority'
});

/**
 * Where a drawing binds a ref: the slot (`FILE_KIND_OF_SLOT`'s keys, `committed` for a committed
 * comparison), the Linear record row of a `row.` slot or the Circular card of a `source.` slot by
 * its uid, and the Linear comparison of
 * `edge` by its id (null for the other slots).
 * @typedef {object} BoundFilePlace
 * @property {string} slot
 * @property {string | null} rowUid
 * @property {string | null} edgeId
 */
/**
 * Every ref one drawing binds, with its place: the one walk of a drawing's bindings.
 * @param {BoundDrawingInputs} drawing
 * @param {(ref: InputFileRef, place: BoundFilePlace) => void} visit
 */
const visitDrawingFiles = (drawing, visit) => {
  /** @param {BoundFilePlace} place @param {unknown} value */
  const add = (place, value) => {
    if (Array.isArray(value)) {
      value.forEach((entry) => add(place, entry));
      return;
    }
    if (isInputFileRef(value)) visit(value, place);
  };
  /** @param {string} slot */
  const field = (slot) => ({ slot, rowUid: null, edgeId: null });
  BOUND_FILE_KEYS.forEach((key) => add(field(key), drawing.files[key]));
  (Array.isArray(drawing.files.linearCanonicalComparisons) ? drawing.files.linearCanonicalComparisons : [])
    .forEach((/** @type {{ file?: unknown } | null | undefined} */ comparison) => add(field('committed'), comparison?.file));
  drawing.linearSeqs.forEach((sequence) => LINEAR_SEQ_FILE_KEYS.forEach((key) => (
    add({ slot: `row.${key}`, rowUid: sequence.uid, edgeId: null }, sequence[key]))));
  drawing.circularSources.forEach((source) => CIRCULAR_SOURCE_FILE_SLOTS.forEach((key) => (
    add({ slot: `source.${key}`, rowUid: source.uid, edgeId: null }, source[key]))));
  (drawing.linearComparisonPlan?.edges || []).forEach((edge) => add({ slot: 'edge', rowUid: null, edgeId: edge.id }, edge.file));
};
/**
 * Every ref the bindings name, with the slot that names it, through the walk of each drawing.
 * @param {InputBindings} bindings
 * @param {(ref: InputFileRef, slot: string) => void} visit
 */
const visitBoundFiles = (bindings, visit) => {
  bindings.drawingList().forEach((drawing) => visitDrawingFiles(drawing, (ref, place) => visit(ref, place.slot)));
};

/**
 * A place an input file replacement writes: a drawing by its ID, and the place in it.
 * @template {string} Id the drawing ID type of the caller
 * @typedef {object} InputReplacementSite
 * @property {Id} drawingId
 * @property {ReplaceableInputSlot} slot
 * @property {string | null} rowUid
 * @property {string | null} edgeId
 */
/** @param {string} slot @returns {slot is ReplaceableInputSlot | CardFileSlot} */
const isUploaderSlot = (slot) => Object.hasOwn(FILE_KIND_OF_SLOT, slot);
/** @param {string} slot @returns {slot is ReplaceableInputSlot} */
const isReplaceableSlot = (slot) => isUploaderSlot(slot) && !CARD_FILE_SLOTS.some((card) => card === slot);
/**
 * The places where the drawings bind the file `fileId` that a replacement of it writes
 * (PD-OI-090 〔入力ファイル〕, R17-3): every uploader slot except the five cards' (each
 * drawing binds its own card files, so a card never asks), one site per place. A committed
 * comparison and the conservation sequences have no uploader and are not offered.
 * @template {string} Id
 * @param {readonly { drawingId: Id, drawing: BoundDrawingInputs }[]} drawings
 * @param {string} fileId
 * @returns {InputReplacementSite<Id>[]}
 */
export const projectInputReplacement = (drawings, fileId) => {
  /** @type {InputReplacementSite<Id>[]} */
  const sites = [];
  /** @type {Set<string>} */
  const seen = new Set();
  drawings.forEach(({ drawingId, drawing }) => visitDrawingFiles(drawing, (ref, place) => {
    const { slot, rowUid, edgeId } = place;
    if (ref.fileId !== fileId || !isReplaceableSlot(slot)) return;
    const key = JSON.stringify([drawingId, slot, rowUid, edgeId]);
    if (seen.has(key)) return;
    seen.add(key);
    sites.push({ drawingId, slot, rowUid, edgeId });
  }));
  return sites;
};

/**
 * The IDs of every file a binding names.
 * @param {InputBindings} bindings
 * @returns {Set<string>}
 */
export const boundInputFileIds = (bindings) => {
  /** @type {Set<string>} */
  const fileIds = new Set();
  visitBoundFiles(bindings, (ref) => fileIds.add(ref.fileId));
  return fileIds;
};

/**
 * The files in this Session that an uploader of `kind` offers (PD-OI-090
 * 〔入力ファイル〕, T-B4-6): every ref any drawing binds in a slot of that
 * kind, once per file, in binding order. Choosing one binds the ref itself,
 * so no file is stored twice. An uploader that takes several kinds (the
 * Circular Add files zone of GFF3 + FASTA) offers the files of each.
 * @param {InputBindings} bindings
 * @param {InputFileKind | readonly InputFileKind[]} kind
 * @returns {InputFileRef[]}
 */
export const sessionFilesOfKind = (bindings, kind) => {
  const kinds = typeof kind === 'string' ? [kind] : kind;
  /** @type {Map<string, InputFileRef>} */
  const refs = new Map();
  visitBoundFiles(bindings, (ref, slot) => {
    if (isUploaderSlot(slot) && kinds.includes(FILE_KIND_OF_SLOT[slot]) && !refs.has(ref.fileId)) refs.set(ref.fileId, ref);
  });
  return [...refs.values()];
};

/**
 * The kind of the files a Circular source card's slot binds.
 * @param {'gb' | 'gff' | 'fasta'} slot
 * @returns {InputFileKind}
 */
export const circularSlotFileKind = (slot) => FILE_KIND_OF_SLOT[`source.${slot}`];

/**
 * The store's release for History's budget: it keeps the files History names
 * and every bound file.
 * @param {InputBindings} bindings
 * @param {Pick<InputFileStore, 'retainOnly'>} store
 * @returns {(historyFileIds: Iterable<string>) => void}
 */
export const retainHistoryAndBoundFiles = (bindings, store) => (historyFileIds) => {
  const retained = boundInputFileIds(bindings);
  for (const fileId of historyFileIds) retained.add(fileId);
  store.retainOnly(retained);
};

/**
 * A loaded value of a binding: a File or Session resource view (a Load
 * candidate), a ref (a rolled-back snapshot), or a list of them.
 * @typedef {InputFileCandidate | InputFileRef | readonly LoadedFileValue[]} LoadedFileValue
 */

/**
 * The refs of a loaded value, each place and each nested list kept.
 * @param {Pick<InputFileStore, 'registerFile'>} store
 * @param {LoadedFileValue} value
 * @returns {BoundFileValue}
 */
export const registerInputFileTree = (store, value) => (
  Array.isArray(value)
    ? value.map((entry) => registerInputFileTree(store, entry))
    : store.registerFile(/** @type {InputFileCandidate} */ (value))
);

/**
 * The comparisons a Session committed, as the bindings hold them: a
 * resource-backed entry's loaded file is registered and named by its ref.
 * @param {Pick<InputFileStore, 'registerFile'>} store
 * @param {readonly (Omit<BoundCanonicalComparison, 'file'> & { file?: unknown })[]} comparisons
 * @returns {BoundCanonicalComparison[]}
 */
export const bindCanonicalComparisons = (store, comparisons) => comparisons.map((comparison) => (
  Object.hasOwn(comparison, 'file')
    ? { ...comparison, file: store.registerFile(/** @type {InputFileCandidate} */ (comparison.file)) }
    : /** @type {BoundCanonicalComparison} */ (comparison)
));

/**
 * The input fields of loaded or captured `files`, each file registered and
 * named by its ref (a released ref names nothing: null, and a place of
 * `c_conservation_blasts`, a list of files only, is dropped). The one
 * registration of a drawing's fields for Load and History.
 * @param {Pick<InputFileStore, 'registerFile'>} store
 * @param {Readonly<Record<string, unknown>>} loaded
 * @returns {Record<string, unknown>}
 */
const bindFileFields = (store, loaded) => {
  /** @type {Record<string, unknown>} */
  const files = { ...loaded };
  BOUND_FILE_KEYS.forEach((key) => {
    if (!Object.hasOwn(files, key)) return;
    const bound = registerInputFileTree(store, /** @type {LoadedFileValue} */ (files[key]));
    files[key] = RING_GENOME_KEYS.has(key) ? boundRingGenomeLists(bound) : bound;
  });
  if (Array.isArray(files.c_conservation_blasts)) files.c_conservation_blasts = files.c_conservation_blasts.filter((ref) => ref !== null);
  if (Array.isArray(files.linearCanonicalComparisons)) {
    files.linearCanonicalComparisons = bindCanonicalComparisons(store, files.linearCanonicalComparisons);
  }
  // A card mark names its files by their refs, as the cards do.
  CARD_MARK_KEYS.forEach((key) => {
    const marks = /** @type {readonly CircularCardMark[] | null | undefined} */ (files[key]);
    if (Array.isArray(marks)) {
      files[key] = marks.map((mark) => ({ ...mark, primaryFile: store.registerFile(mark.primaryFile), pairedFile: store.registerFile(mark.pairedFile) }));
    }
  });
  return files;
};

/**
 * A Linear record row with its files registered and named by their refs.
 * @param {Pick<InputFileStore, 'registerFile'>} store
 * @param {LoadedLinearSeq} seq
 * @returns {LinearSeq}
 */
const bindLinearSeqRow = (store, seq) => ({
  ...seq,
  gb: store.registerFile(/** @type {InputFileCandidate} */ (seq.gb)),
  gff: store.registerFile(/** @type {InputFileCandidate} */ (seq.gff)),
  fasta: store.registerFile(/** @type {InputFileCandidate} */ (seq.fasta)),
  depth: /** @type {DepthFileSlots<InputFileRef>} */ (registerInputFileTree(store, seq.depth))
});

/**
 * A Circular source card with its files registered and named by their refs.
 * @param {Pick<InputFileStore, 'registerFile'>} store
 * @param {LoadedCircularSource} source
 * @returns {CircularSource}
 */
const bindCircularSourceRow = (store, source) => ({
  uid: source.uid,
  gb: store.registerFile(/** @type {InputFileCandidate} */ (source.gb)),
  gff: store.registerFile(/** @type {InputFileCandidate} */ (source.gff)),
  fasta: store.registerFile(/** @type {InputFileCandidate} */ (source.fasta))
});

/**
 * A Circular source card as a Load candidate or a rolled-back snapshot holds it.
 * @typedef {{ uid: string } & Record<CircularSourceFileSlot, LoadedFileValue | null>} LoadedCircularSource
 */
/**
 * A Linear record row as a Load candidate or a rolled-back snapshot holds it.
 * @typedef {Omit<LinearSeq, 'gb' | 'gff' | 'fasta' | 'depth'> & {
 *   gb: LoadedFileValue, gff: LoadedFileValue, fasta: LoadedFileValue, depth: LoadedFileValue
 * }} LoadedLinearSeq
 */
/**
 * A comparison plan as a Load candidate or a rolled-back snapshot holds it.
 * @typedef {Omit<LinearComparisonPlan, 'edges'> & {
 *   edges: (Omit<LinearComparisonEdge, 'file'> & { file: LoadedFileValue })[]
 * }} LoadedComparisonPlan
 */
/**
 * The inputs a Load installs (its candidate, which holds the loaded Files) or
 * a rollback restores (refs): the drawing's fields and the card files.
 * @typedef {object} LoadedInputs
 * @property {Readonly<Record<string, unknown>>} files
 * @property {readonly LoadedLinearSeq[]} linearSeqs
 * @property {readonly LoadedCircularSource[]} circularSources
 * @property {LoadedComparisonPlan} linearComparisonPlan
 * @property {{ sources?: readonly { uid: string, primaryFile: LoadedFileValue, pairedFile: LoadedFileValue, error: CircularRecordDiscovery['error'] }[] }} circularRecordDiscovery
 */
/**
 * The bindings a Load installs: every file it holds is registered here, at
 * install, and each field holds its ref (a ref of a rollback stays bound).
 * The Load candidate's own reads before the commit read its Files directly.
 * @param {Pick<InputFileStore, 'registerFile'>} store
 * @param {LoadedInputs} loaded
 */
export const bindLoadedInputs = (store, loaded) => {
  /** @param {unknown} value */
  const one = (value) => store.registerFile(/** @type {InputFileCandidate} */ (value));
  const files = bindFileFields(store, loaded.files);
  const linearSeqs = loaded.linearSeqs.map((seq) => bindLinearSeqRow(store, seq));
  const circularSources = loaded.circularSources.map((source) => bindCircularSourceRow(store, source));
  /** @type {LinearComparisonPlan} */
  const linearComparisonPlan = {
    ...loaded.linearComparisonPlan,
    edges: loaded.linearComparisonPlan.edges.map((edge) => ({ ...edge, file: one(edge.file) }))
  };
  const discovery = loaded.circularRecordDiscovery;
  return {
    files,
    linearSeqs,
    circularSources,
    linearComparisonPlan,
    circularRecordDiscovery: {
      ...discovery,
      ...(discovery.sources ? { sources: discovery.sources.map((source) => ({
        ...source, primaryFile: one(source.primaryFile), pairedFile: one(source.pairedFile)
      })) } : {})
    }
  };
};

// The inputs slice of a drawing (DESIGN §2): its input fields (`files`, its
// five card files among them, R16-4), its Linear record rows
// (`linearSeqs`), its Circular source cards (`circularSources`), its input
// type, and its Circular record catalog. A drawing
// shows them under their member names; its file fields and rows are written
// only through the writers below, which take the slice's handle
// (`drawing.slices.inputs`); the input type and the record catalog stay
// written by their controls and discovery, as settings are.

/** @typedef {import('../input-types.js').InputType} InputType */
/** @typedef {import('../input-types.js').InputSlot} InputSlot */
/** One Circular record as discovery lists it. @typedef {import('../record-discovery.js').CircularListRecord} CircularRecordEntry */
/**
 * What `capture` returns and `install` writes.
 * @typedef {object} InputsData
 * @property {DrawingInputFiles} files
 * @property {LinearSeq[]} linearSeqs
 * @property {CircularSource[]} circularSources
 * @property {InputType} inputType
 * @property {CircularRecordEntry[]} circularRecordList
 * @property {CircularRecordDiscovery} circularRecordDiscovery
 */
/**
 * The members as a drawing holds them.
 * @typedef {object} InputsMembers
 * @property {DrawingInputFiles} files
 * @property {LinearSeq[]} linearSeqs
 * @property {CircularSource[]} circularSources
 * @property {import('../../utils/reactivity.js').Ref<InputType>} inputType
 * @property {import('../../utils/reactivity.js').Ref<CircularRecordEntry[]>} circularRecordList
 * @property {CircularRecordDiscovery} circularRecordDiscovery
 */
/** @typedef {Readonly<InputsMembers>} InputsView */
/** @typedef {SliceHandle<'inputs'>} InputsSlice The handle a writer takes. */

/** @returns {DrawingInputFiles} */
const emptyDrawingInputFiles = () => ({
  c_depth: null,
  c_conservation_blasts: [], c_conservation_blasts_source: null, c_conservation_fastas: [],
  c_conservation_sequence_sources: [], linearCanonicalComparisons: [],
  d_color: null, t_color: null, blacklist: null, whitelist: null, qualifier_priority: null,
  c_depth_rows: null, c_legacy_record_cards: null
});
const DRAWING_FILE_KEYS = /** @type {(keyof DrawingInputFiles)[]} */ (Object.keys(emptyDrawingInputFiles()));

/**
 * The record of a drawing with nothing to read: a new drawing, Load before it
 * binds files, and a discovery without its complete input all write this one.
 * @returns {CircularRecordDiscovery}
 */
export const idleCircularRecordDiscovery = () => ({
  status: 'idle', error: '', inputType: '', sources: []
});

// An inputs handle exists only through `create`, which gives it these members.
/** @param {InputsSlice} slice */
const membersOf = (slice) => /** @type {InputsMembers} */ (sliceMembers(slice, 'inputs'));

/**
 * A copy of a bound value: lists copied at every depth (a Depth matrix, the ring genomes), refs
 * by reference (they are frozen).
 * @template T
 * @param {T} value
 * @returns {T}
 */
const copyBound = (value) => /** @type {T} */ (Array.isArray(value) ? value.map(copyBound) : value);

/** @param {DrawingInputFiles} files @returns {DrawingInputFiles} */
const copyFiles = (files) => ({
  ...files,
  c_depth: copyBound(files.c_depth),
  c_conservation_blasts: [...files.c_conservation_blasts],
  c_conservation_fastas: copyBound(files.c_conservation_fastas),
  c_conservation_sequence_sources: copyBound(files.c_conservation_sequence_sources),
  linearCanonicalComparisons: files.linearCanonicalComparisons.map((comparison) => ({ ...comparison }))
});
/** @param {LinearSeq} seq @returns {LinearSeq} */
const copyRow = (seq) => ({ ...seq, depth: copyBound(seq.depth) });

// Each member's build, copy, and write.
/** @typedef {keyof InputsMembers} InputsMember */
const MEMBER_OPS = Object.freeze({
  files: {
    create: () => /** @type {DrawingInputFiles} */ (reactive(emptyDrawingInputFiles())),
    /** @param {InputsMembers} members */
    capture: (members) => copyFiles(members.files),
    /** @param {InputsMembers} members @param {InputsData} data */
    install: (members, data) => {
      const files = copyFiles({ ...emptyDrawingInputFiles(), ...data.files });
      DRAWING_FILE_KEYS.forEach((key) => { /** @type {Record<string, unknown>} */ (members.files)[key] = files[key]; });
    }
  },
  linearSeqs: {
    create: () => /** @type {LinearSeq[]} */ (reactive(normalizeLinearSeqList([]))),
    /** @param {InputsMembers} members */
    capture: (members) => members.linearSeqs.map(copyRow),
    /** @param {InputsMembers} members @param {InputsData} data */
    install: (members, data) => {
      members.linearSeqs.splice(0, members.linearSeqs.length, ...data.linearSeqs.map(copyRow));
    }
  },
  circularSources: {
    create: () => /** @type {CircularSource[]} */ (reactive([])),
    /** @param {InputsMembers} members */
    capture: (members) => members.circularSources.map((source) => ({ ...source })),
    /** @param {InputsMembers} members @param {InputsData} data */
    install: (members, data) => {
      members.circularSources.splice(0, members.circularSources.length, ...data.circularSources.map((source) => ({ ...source })));
    }
  },
  inputType: {
    create: () => /** @type {import('../../utils/reactivity.js').Ref<InputType>} */ (ref('gb')),
    /** @param {InputsMembers} members */
    capture: (members) => members.inputType.value,
    /** @param {InputsMembers} members @param {InputsData} data */
    install: (members, data) => { members.inputType.value = data.inputType; }
  },
  circularRecordList: {
    create: () => /** @type {import('../../utils/reactivity.js').Ref<CircularRecordEntry[]>} */ (ref([])),
    /** @param {InputsMembers} members */
    capture: (members) => members.circularRecordList.value.map((entry) => ({ ...entry })),
    /** @param {InputsMembers} members @param {InputsData} data */
    install: (members, data) => { members.circularRecordList.value = data.circularRecordList.map((entry) => ({ ...entry })); }
  },
  circularRecordDiscovery: {
    create: () => /** @type {CircularRecordDiscovery} */ (reactive(idleCircularRecordDiscovery())),
    /** @param {InputsMembers} members */
    capture: (members) => ({
      ...members.circularRecordDiscovery,
      sources: members.circularRecordDiscovery.sources.map((source) => ({ ...source }))
    }),
    /** @param {InputsMembers} members @param {InputsData} data */
    install: (members, data) => {
      Object.assign(members.circularRecordDiscovery, idleCircularRecordDiscovery(), data.circularRecordDiscovery, {
        sources: data.circularRecordDiscovery.sources.map((source) => ({ ...source }))
      });
    }
  }
});
/** @type {readonly InputsMember[]} */
const MEMBERS = Object.freeze(['files', 'linearSeqs', 'circularSources', 'inputType', 'circularRecordList', 'circularRecordDiscovery']);

/**
 * A new slice, or one over a drawing's existing members (a test fixture).
 * @param {Partial<InputsMembers>} [members]
 * @returns {{ handle: InputsSlice, view: InputsView }}
 */
const create = (members = {}) => {
  const data = /** @type {InputsMembers} */ (Object.fromEntries(MEMBERS.map((member) => (
    [member, members[member] ?? MEMBER_OPS[member].create()]
  ))));
  return { handle: createSliceHandle('inputs', data), view: Object.freeze({ ...data }) };
};

/**
 * @param {InputsSlice} slice
 * @returns {InputsData}
 */
const capture = (slice) => {
  const members = membersOf(slice);
  return /** @type {InputsData} */ (Object.fromEntries(MEMBERS.map((member) => (
    [member, MEMBER_OPS[member].capture(members)]
  ))));
};

// Writes each member `data` holds; a member it leaves out stays. The rows are
// installed as given: their builders (Load, a History step) normalize them.
/**
 * @param {InputsSlice} slice
 * @param {Partial<InputsData>} data
 */
const install = (slice, data) => {
  const members = membersOf(slice);
  MEMBERS.forEach((member) => {
    if (data[member]) MEMBER_OPS[member].install(members, /** @type {InputsData} */ (data));
  });
};

// Reset keeps every input binding, the chosen records, and the Depth
// assignment (PD-OI-070), but removes the drawing's five card files with the
// settings they loaded (D-15-78).
/** @param {InputsSlice} slice */
const reset = (slice) => {
  const { files } = membersOf(slice);
  CARD_FILE_SLOTS.forEach((slot) => { files[slot] = null; });
};

/**
 * What a History step holds of a drawing's inputs: a capture (`capture`) of
 * its file fields, rows, and input type, whose refs may have been released
 * since. A field `files` leaves out binds nothing.
 * @typedef {object} CapturedInputs
 * @property {Partial<DrawingInputFiles>} files
 * @property {readonly LinearSeq[]} linearSeqs
 * @property {readonly CircularSource[]} circularSources
 * @property {InputType} [inputType] left out, the input type stays
 */
/**
 * Install a History step's capture: each ref is registered with the store (a
 * released ref binds nothing), the rows are normalized as every row writer
 * leaves them, and the `files` fields `keep` names stay as they are (an
 * intent leaves out the comparisons a Generate wrote). This is the one writer
 * of a captured slice.
 * @param {Pick<InputFileStore, 'registerFile'>} store
 * @param {InputsSlice} slice
 * @param {CapturedInputs} data
 * @param {{ keep?: readonly (keyof DrawingInputFiles)[] }} [options]
 */
export const installBoundInputs = (store, slice, data, { keep = [] } = {}) => {
  const current = membersOf(slice).files;
  install(slice, {
    files: /** @type {DrawingInputFiles} */ ({
      ...bindFileFields(store, data.files),
      ...Object.fromEntries(keep.map((key) => [key, current[key]]))
    }),
    linearSeqs: normalizeLinearSeqList(data.linearSeqs.map((seq) => bindLinearSeqRow(store, seq))),
    circularSources: data.circularSources.map((source) => bindCircularSourceRow(store, source)),
    ...(data.inputType ? { inputType: data.inputType } : {})
  });
};

/**
 * Bind a file to one of the drawing's five card fields: the store registers
 * it and the field holds its ref (null unbinds the card). This is the one
 * writer of a card file field.
 * @param {Pick<InputFileStore, 'registerFile'>} store
 * @param {InputsSlice} slice
 * @param {CardFileSlot} slot
 * @param {InputFileCandidate} file
 * @returns {InputFileRef | null}
 */
export const bindInputFile = (store, slice, slot, file) => {
  const ref = store.registerFile(file);
  membersOf(slice).files[slot] = ref;
  return ref;
};

/**
 * Bind the files of one of the drawing's list fields; a place without a file
 * is dropped.
 * @param {Pick<InputFileStore, 'registerFile'>} store
 * @param {InputsSlice} slice
 * @param {InputFileListSlot} slot
 * @param {InputFileCandidate | readonly InputFileCandidate[]} list
 * @returns {InputFileRef[]}
 */
export const bindDrawingInputFileList = (store, slice, slot, list) => {
  const refs = registerInputFiles(store, list).filter((entry) => entry !== null);
  membersOf(slice).files[slot] = refs;
  return refs;
};

/**
 * Bind the drawing's Depth matrix, each place kept (null for none), or its
 * rings' comparison genome files or sequence sources, one list per ring
 * (`boundRingGenomeLists`).
 * @param {Pick<InputFileStore, 'registerFile'>} store
 * @param {InputsSlice} slice
 * @param {'c_depth' | 'c_conservation_fastas' | 'c_conservation_sequence_sources'} slot
 * @param {LoadedFileValue} value
 */
export const bindDrawingInputFileTree = (store, slice, slot, value) => {
  const bound = value === null || value === undefined ? (slot === 'c_depth' ? null : []) : registerInputFileTree(store, value);
  /** @type {Record<string, unknown>} */ (membersOf(slice).files)[slot] = slot === 'c_depth' ? bound : boundRingGenomeLists(bound);
};

/**
 * Mark the cards whose records the rows of the drawing's record-major Depth
 * matrix stand for (`c_depth_rows`, OV-2868): a card edit marks the cards as
 * last read, and a settled read marks the cards it read once the rows follow
 * them. The one writer of the field besides Load and a History restore.
 * @param {InputsSlice} slice
 * @param {readonly CircularCardMark[] | null} marks
 */
export const markCircularDepthRows = (slice, marks) => {
  membersOf(slice).files.c_depth_rows = marks;
};

/**
 * Whether the drawing's conservation BLAST files were derived from the LOSAT
 * cache (`'losat-cache'`) or uploaded (null).
 * @param {InputsSlice} slice
 * @param {'losat-cache' | null} source
 */
export const setConservationBlastsSource = (slice, source) => {
  membersOf(slice).files.c_conservation_blasts_source = source;
};

/**
 * The Linear comparisons the drawing committed, a resource-backed entry's
 * file by ref.
 * @param {Pick<InputFileStore, 'registerFile'>} store
 * @param {InputsSlice} slice
 * @param {readonly (Omit<BoundCanonicalComparison, 'file'> & { file?: unknown })[]} comparisons
 */
export const bindDrawingCanonicalComparisons = (store, slice, comparisons) => {
  membersOf(slice).files.linearCanonicalComparisons = bindCanonicalComparisons(store, comparisons);
};

// The card lists of a drawing: its Linear record rows and its Circular source
// cards. One owner adds, removes, moves, and replaces the rows of either list
// by uid; a Linear drawing always holds a row, a Circular one may hold none
// (its empty state).
/** @typedef {'linearSeqs' | 'circularSources'} InputRowList */
/**
 * @template {InputRowList} L
 * @typedef {L extends 'linearSeqs' ? LinearSeq : CircularSource} InputRowOf
 */
/**
 * @template {InputRowList} L
 * @param {InputsSlice} slice
 * @param {L} list
 * @returns {InputRowOf<L>[]}
 */
const rowsOf = (slice, list) => /** @type {InputRowOf<L>[]} */ (membersOf(slice)[list]);

/**
 * The index of the row of one uid; an unknown uid is refused.
 * @param {InputsSlice} slice
 * @param {InputRowList} list
 * @param {string} uid
 */
const rowIndexOf = (slice, list, uid) => {
  const index = rowsOf(slice, list).findIndex((row) => row.uid === uid);
  if (index < 0) throw diagnosticError('INPUT_INVALID', { field: list, reason: 'UNKNOWN_ROW', uid }, { stage: 'drawing-inputs' });
  return index;
};

/**
 * @template {InputRowList} L
 * @param {InputsSlice} slice
 * @param {L} list
 * @param {InputRowOf<L>} row
 * @param {number} [index]
 */
const addRow = (slice, list, row, index) => {
  const rows = rowsOf(slice, list);
  if (rows.some((entry) => entry.uid === row.uid)) throw diagnosticError('INPUT_INVALID', { field: list, reason: 'DUPLICATE_ROW', uid: row.uid }, { stage: 'drawing-inputs' });
  rows.splice(index === undefined ? rows.length : index, 0, row);
};

/** @param {InputsSlice} slice @param {InputRowList} list @param {string} uid */
const removeRow = (slice, list, uid) => {
  const rows = rowsOf(slice, list);
  rows.splice(rowIndexOf(slice, list, uid), 1);
  if (list === 'linearSeqs' && rows.length === 0) rows.push(createLinearSeq());
};

/** @param {InputsSlice} slice @param {InputRowList} list @param {string} uid @param {number} index */
const moveRow = (slice, list, uid, index) => {
  const rows = rowsOf(slice, list);
  const [row] = rows.splice(rowIndexOf(slice, list, uid), 1);
  rows.splice(Math.max(0, Math.min(index, rows.length)), 0, row);
};

/**
 * @template {InputRowList} L
 * @param {InputsSlice} slice
 * @param {L} list
 * @param {readonly InputRowOf<L>[]} rows
 */
const replaceRows = (slice, list, rows) => {
  const current = rowsOf(slice, list);
  current.splice(0, current.length, ...rows);
};

/**
 * Bind a file to one field of a Linear record row, found by its uid.
 * @param {Pick<InputFileStore, 'registerFile'>} store
 * @param {InputsSlice} slice
 * @param {string} uid
 * @param {'gb' | 'gff' | 'fasta' | 'depth'} key
 * @param {LoadedFileValue | null} file
 */
export const bindLinearSeqFile = (store, slice, uid, key, file) => {
  const row = rowsOf(slice, 'linearSeqs')[rowIndexOf(slice, 'linearSeqs', uid)];
  if (key === 'depth') row.depth = /** @type {DepthFileSlots<InputFileRef>} */ (file === null ? null : registerInputFileTree(store, file));
  else row[key] = store.registerFile(/** @type {InputFileCandidate} */ (file));
};

/**
 * Bind a file (or none) to one slot of a Circular source card, found by its uid.
 * @param {Pick<InputFileStore, 'registerFile'>} store
 * @param {InputsSlice} slice
 * @param {string} uid
 * @param {'gb' | 'gff' | 'fasta'} key
 * @param {InputFileCandidate | null} file
 * @returns {InputFileRef | null}
 */
export const bindCircularSourceFile = (store, slice, uid, key, file) => {
  const row = rowsOf(slice, 'circularSources')[rowIndexOf(slice, 'circularSources', uid)];
  row[key] = store.registerFile(/** @type {InputFileCandidate} */ (file));
  return row[key];
};

/**
 * Add a row (built by `createLinearSeq`) at `index`, or at the end.
 * @param {InputsSlice} slice
 * @param {LinearSeq} row
 * @param {number} [index]
 */
export const addLinearSeq = (slice, row, index) => addRow(slice, 'linearSeqs', row, index);
/**
 * Add a card (built by `createCircularSource`) at `index`, or at the end.
 * @param {InputsSlice} slice
 * @param {CircularSource} row
 * @param {number} [index]
 */
export const addCircularSource = (slice, row, index) => addRow(slice, 'circularSources', row, index);

/**
 * Remove the row of one uid; the last Linear row is replaced by an empty one,
 * so a Linear drawing always holds a row.
 * @param {InputsSlice} slice
 * @param {string} uid
 */
export const removeLinearSeq = (slice, uid) => removeRow(slice, 'linearSeqs', uid);
/** @param {InputsSlice} slice @param {string} uid */
export const removeCircularSource = (slice, uid) => removeRow(slice, 'circularSources', uid);

/**
 * Move the row of one uid to `index`.
 * @param {InputsSlice} slice
 * @param {string} uid
 * @param {number} index
 */
export const moveLinearSeq = (slice, uid, index) => moveRow(slice, 'linearSeqs', uid, index);
/** @param {InputsSlice} slice @param {string} uid @param {number} index */
export const moveCircularSource = (slice, uid, index) => moveRow(slice, 'circularSources', uid, index);

/**
 * Replace the drawing's rows by a list a transition built (normalized by
 * `normalizeLinearSeqList`).
 * @param {InputsSlice} slice
 * @param {readonly LinearSeq[]} rows
 */
export const replaceLinearSeqs = (slice, rows) => replaceRows(slice, 'linearSeqs', rows);
/** @param {InputsSlice} slice @param {readonly CircularSource[]} rows */
export const replaceCircularSources = (slice, rows) => replaceRows(slice, 'circularSources', rows);

export const INPUTS_SLICE = Object.freeze({
  members: MEMBERS,
  create, capture, install, reset
});
