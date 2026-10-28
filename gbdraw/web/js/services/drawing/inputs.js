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
/** @import { SliceHandle } from './slice-handle.js' */
import { isInputFileRef } from '../history-files.js';
import { diagnosticError } from '../../utils/error-normalization.js';
import { reactive, ref } from '../../utils/reactivity.js';
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
 * @property {number | ''} losat_gencode `v-model.number` writes an empty string when the field is cleared
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
 */

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
  const rawLosatGencode = Number(source.losat_gencode);
  return {
    uid: rawUid || generateLinearSeqUid(),
    gb: source.gb ?? null,
    gff: source.gff ?? null,
    fasta: source.fasta ?? null,
    depth: source.depth ?? null,
    losat_gencode: Number.isFinite(rawLosatGencode) && rawLosatGencode > 0 ? rawLosatGencode : 1,
    definition: String(source.definition ?? ''),
    record_subtitle: String(source.record_subtitle ?? ''),
    file_definition: String(source.file_definition ?? ''),
    file_subtitle: String(source.file_subtitle ?? ''),
    inferred_definition: String(source.inferred_definition ?? ''),
    region_record_id: String(source.region_record_id ?? ''),
    region_start: normalizeLinearSeqNumber(source.region_start),
    region_end: normalizeLinearSeqNumber(source.region_end),
    region_reverse: Boolean(source.region_reverse)
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


/**
 * A drawing's bindings, as the retention walk reads them: its uploaded
 * comparison files and, from the inputs slice, its input fields and record
 * rows.
 * @typedef {object} BoundDrawingInputs
 * @property {{ readonly edges?: readonly { readonly id: string, readonly file?: unknown }[] } | null} [linearComparisonPlan]
 * @property {Readonly<Record<string, unknown>>} files
 * @property {readonly (Readonly<Record<string, unknown>> & { readonly uid: string })[]} linearSeqs
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
 * @property {InputFileRef | null} c_gb
 * @property {InputFileRef | null} c_gff
 * @property {InputFileRef | null} c_fasta
 * @property {InputFileRef | (InputFileRef | null)[] | (InputFileRef | null)[][] | null} c_depth one record's series, or a record-major matrix
 * @property {InputFileRef[]} c_conservation_blasts
 * @property {string | null} c_conservation_blasts_source a flag (`'losat-cache'`), not a file
 * @property {(InputFileRef | null)[]} c_conservation_fastas
 * @property {(InputFileRef | null)[]} c_conservation_sequence_sources
 * @property {BoundCanonicalComparison[]} linearCanonicalComparisons
 * @property {InputFileRef | null} d_color the Default colors card's file
 * @property {InputFileRef | null} t_color the Specific rules card's file
 * @property {InputFileRef | null} blacklist the Blacklist card's file
 * @property {InputFileRef | null} whitelist the Whitelist card's file
 * @property {InputFileRef | null} qualifier_priority the Priority card's file
 */
/**
 * The five card files' fields (Default colors, Specific rules, Blacklist,
 * Whitelist, Priority): each drawing binds its own (R16-4).
 * @typedef {'d_color' | 't_color' | 'blacklist' | 'whitelist' | 'qualifier_priority'} CardFileSlot
 */
/** @type {readonly CardFileSlot[]} */
const CARD_FILE_SLOTS = Object.freeze(['d_color', 't_color', 'blacklist', 'whitelist', 'qualifier_priority']);

/** The `DrawingInputFiles` fields that hold a list of refs. @typedef {'c_conservation_blasts' | 'c_conservation_fastas'} InputFileListSlot */

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
  'c_gb', 'c_gff', 'c_fasta', 'c_depth', 'c_conservation_blasts', 'c_conservation_fastas',
  'c_conservation_sequence_sources', 'd_color', 't_color', 'blacklist', 'whitelist', 'qualifier_priority'
]);
const LINEAR_SEQ_FILE_KEYS = Object.freeze(['gb', 'gff', 'fasta', 'depth']);

/**
 * The kinds of file an uploader takes; a binding slot's kind is the kind of
 * the files it binds. Each of the five cards is its own kind.
 * @typedef {'genbank' | 'gff' | 'fasta' | 'depth' | 'comparison' | CardFileSlot} InputFileKind
 */
/**
 * The kind of each slot that an uploader binds; a bound slot missing here
 * (the conservation sequences, a committed comparison) is offered to no
 * uploader. `row.` names a Linear record row's field, `edge` a Linear
 * comparison's uploaded file.
 * @type {Readonly<Record<string, InputFileKind>>}
 */
const FILE_KIND_OF_SLOT = Object.freeze({
  c_gb: 'genbank', 'row.gb': 'genbank', c_gff: 'gff', 'row.gff': 'gff', c_fasta: 'fasta', 'row.fasta': 'fasta',
  c_depth: 'depth', 'row.depth': 'depth', c_conservation_blasts: 'comparison', edge: 'comparison',
  d_color: 'd_color', t_color: 't_color', blacklist: 'blacklist', whitelist: 'whitelist',
  qualifier_priority: 'qualifier_priority'
});

/**
 * Where a drawing binds a ref: the slot (`FILE_KIND_OF_SLOT`'s keys, `committed` for a committed
 * comparison), the Linear record row of a `row.` slot by its uid, and the Linear comparison of
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
 * @typedef {BoundFilePlace & { drawingId: Id }} InputReplacementSite
 */
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
    const kind = FILE_KIND_OF_SLOT[place.slot];
    if (ref.fileId !== fileId || !kind || CARD_FILE_SLOTS.some((slot) => slot === kind)) return;
    const key = JSON.stringify([drawingId, place.slot, place.rowUid, place.edgeId]);
    if (seen.has(key)) return;
    seen.add(key);
    sites.push({ drawingId, ...place });
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
 * so no file is stored twice.
 * @param {InputBindings} bindings
 * @param {InputFileKind} kind
 * @returns {InputFileRef[]}
 */
export const sessionFilesOfKind = (bindings, kind) => {
  /** @type {Map<string, InputFileRef>} */
  const refs = new Map();
  visitBoundFiles(bindings, (ref, slot) => {
    if (FILE_KIND_OF_SLOT[slot] === kind && !refs.has(ref.fileId)) refs.set(ref.fileId, ref);
  });
  return [...refs.values()];
};

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
    if (Object.hasOwn(files, key)) files[key] = registerInputFileTree(store, /** @type {LoadedFileValue} */ (files[key]));
  });
  if (Array.isArray(files.c_conservation_blasts)) files.c_conservation_blasts = files.c_conservation_blasts.filter((ref) => ref !== null);
  if (Array.isArray(files.linearCanonicalComparisons)) {
    files.linearCanonicalComparisons = bindCanonicalComparisons(store, files.linearCanonicalComparisons);
  }
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
 * @property {LoadedComparisonPlan} linearComparisonPlan
 * @property {{ primaryFile?: LoadedFileValue, pairedFile?: LoadedFileValue }} circularRecordDiscovery
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
  /** @type {LinearComparisonPlan} */
  const linearComparisonPlan = {
    ...loaded.linearComparisonPlan,
    edges: loaded.linearComparisonPlan.edges.map((edge) => ({ ...edge, file: one(edge.file) }))
  };
  const discovery = loaded.circularRecordDiscovery;
  return {
    files,
    linearSeqs,
    linearComparisonPlan,
    circularRecordDiscovery: {
      ...discovery,
      ...(Object.hasOwn(discovery, 'primaryFile') ? { primaryFile: one(discovery.primaryFile) } : {}),
      ...(Object.hasOwn(discovery, 'pairedFile') ? { pairedFile: one(discovery.pairedFile) } : {})
    }
  };
};

// The inputs slice of a drawing (DESIGN §2): its input fields (`files`, its
// five card files among them, R16-4), its Linear record rows
// (`linearSeqs`), its input type, and its Circular record catalog. A drawing
// shows them under their member names; its file fields and rows are written
// only through the writers below, which take the slice's handle
// (`drawing.slices.inputs`); the input type and the record catalog stay
// written by their controls and discovery, as settings are.

/** @typedef {'gb' | 'gff'} InputType */
/** One Circular record as discovery lists it. @typedef {Readonly<Record<string, unknown>>} CircularRecordEntry */
/**
 * What `capture` returns and `install` writes.
 * @typedef {object} InputsData
 * @property {DrawingInputFiles} files
 * @property {LinearSeq[]} linearSeqs
 * @property {InputType} inputType
 * @property {CircularRecordEntry[]} circularRecordList
 * @property {CircularRecordDiscovery} circularRecordDiscovery
 */
/**
 * The members as a drawing holds them.
 * @typedef {object} InputsMembers
 * @property {DrawingInputFiles} files
 * @property {LinearSeq[]} linearSeqs
 * @property {import('../../utils/reactivity.js').Ref<InputType>} inputType
 * @property {import('../../utils/reactivity.js').Ref<CircularRecordEntry[]>} circularRecordList
 * @property {CircularRecordDiscovery} circularRecordDiscovery
 */
/** @typedef {Readonly<InputsMembers>} InputsView */
/** @typedef {SliceHandle<'inputs'>} InputsSlice The handle a writer takes. */

/** @returns {DrawingInputFiles} */
const emptyDrawingInputFiles = () => ({
  c_gb: null, c_gff: null, c_fasta: null, c_depth: null,
  c_conservation_blasts: [], c_conservation_blasts_source: null, c_conservation_fastas: [],
  c_conservation_sequence_sources: [], linearCanonicalComparisons: [],
  d_color: null, t_color: null, blacklist: null, whitelist: null, qualifier_priority: null
});
const DRAWING_FILE_KEYS = /** @type {(keyof DrawingInputFiles)[]} */ (Object.keys(emptyDrawingInputFiles()));

/**
 * The record of a drawing with nothing to read: a new drawing, Load before it
 * binds files, and a discovery without its complete input all write this one.
 * @returns {CircularRecordDiscovery}
 */
export const idleCircularRecordDiscovery = () => ({
  status: 'idle', error: '', inputType: '', primaryFile: null, pairedFile: null, canonicalRecordIdentities: []
});

// An inputs handle exists only through `create`, which gives it these members.
/** @param {InputsSlice} slice */
const membersOf = (slice) => /** @type {InputsMembers} */ (sliceMembers(slice, 'inputs'));

/**
 * A copy of a bound value: lists copied at every depth (a Depth matrix), refs
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
  c_conservation_fastas: [...files.c_conservation_fastas],
  c_conservation_sequence_sources: [...files.c_conservation_sequence_sources],
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
      canonicalRecordIdentities: members.circularRecordDiscovery.canonicalRecordIdentities.map((entry) => ({ ...entry }))
    }),
    /** @param {InputsMembers} members @param {InputsData} data */
    install: (members, data) => {
      Object.assign(members.circularRecordDiscovery, idleCircularRecordDiscovery(), data.circularRecordDiscovery, {
        canonicalRecordIdentities: data.circularRecordDiscovery.canonicalRecordIdentities.map((entry) => ({ ...entry }))
      });
    }
  }
});
/** @type {readonly InputsMember[]} */
const MEMBERS = Object.freeze(['files', 'linearSeqs', 'inputType', 'circularRecordList', 'circularRecordDiscovery']);

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
 * Bind a file to one of the drawing's single-file fields.
 * @param {Pick<InputFileStore, 'registerFile'>} store
 * @param {InputsSlice} slice
 * @param {'c_gb' | 'c_gff' | 'c_fasta'} slot
 * @param {InputFileCandidate} file
 * @returns {InputFileRef | null}
 */
export const bindDrawingInputFile = (store, slice, slot, file) => {
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
 * Bind the drawing's Depth matrix, or its conservation comparison sequences
 * or sequence sources, each place kept (null for none).
 * @param {Pick<InputFileStore, 'registerFile'>} store
 * @param {InputsSlice} slice
 * @param {'c_depth' | 'c_conservation_fastas' | 'c_conservation_sequence_sources'} slot
 * @param {LoadedFileValue} value
 */
export const bindDrawingInputFileTree = (store, slice, slot, value) => {
  const bound = value === null || value === undefined ? (slot === 'c_depth' ? null : []) : registerInputFileTree(store, value);
  /** @type {Record<string, unknown>} */ (membersOf(slice).files)[slot] = bound;
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

/**
 * The drawing's row of one uid; an unknown uid is refused.
 * @param {InputsSlice} slice
 * @param {string} uid
 */
const rowIndexOf = (slice, uid) => {
  const index = membersOf(slice).linearSeqs.findIndex((seq) => seq.uid === uid);
  if (index < 0) throw diagnosticError('INPUT_INVALID', { field: 'linearSeqs', reason: 'UNKNOWN_ROW', uid }, { stage: 'drawing-inputs' });
  return index;
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
  const row = membersOf(slice).linearSeqs[rowIndexOf(slice, uid)];
  if (key === 'depth') row.depth = /** @type {DepthFileSlots<InputFileRef>} */ (file === null ? null : registerInputFileTree(store, file));
  else row[key] = store.registerFile(/** @type {InputFileCandidate} */ (file));
};

/**
 * Add a row (built by `createLinearSeq`) at `index`, or at the end.
 * @param {InputsSlice} slice
 * @param {LinearSeq} row
 * @param {number} [index]
 */
export const addLinearSeq = (slice, row, index) => {
  const rows = membersOf(slice).linearSeqs;
  if (rows.some((seq) => seq.uid === row.uid)) throw diagnosticError('INPUT_INVALID', { field: 'linearSeqs', reason: 'DUPLICATE_ROW', uid: row.uid }, { stage: 'drawing-inputs' });
  rows.splice(index === undefined ? rows.length : index, 0, row);
};

/**
 * Remove the row of one uid; the last row is replaced by an empty one, so a
 * drawing always holds a row.
 * @param {InputsSlice} slice
 * @param {string} uid
 */
export const removeLinearSeq = (slice, uid) => {
  const rows = membersOf(slice).linearSeqs;
  rows.splice(rowIndexOf(slice, uid), 1);
  if (rows.length === 0) rows.push(createLinearSeq());
};

/**
 * Move the row of one uid to `index`.
 * @param {InputsSlice} slice
 * @param {string} uid
 * @param {number} index
 */
export const moveLinearSeq = (slice, uid, index) => {
  const rows = membersOf(slice).linearSeqs;
  const [row] = rows.splice(rowIndexOf(slice, uid), 1);
  rows.splice(Math.max(0, Math.min(index, rows.length)), 0, row);
};

/**
 * Replace the drawing's rows by a list a transition built (normalized by
 * `normalizeLinearSeqList`).
 * @param {InputsSlice} slice
 * @param {readonly LinearSeq[]} rows
 */
export const replaceLinearSeqs = (slice, rows) => {
  const members = membersOf(slice);
  members.linearSeqs.splice(0, members.linearSeqs.length, ...rows);
};

export const INPUTS_SLICE = Object.freeze({
  members: MEMBERS,
  create, capture, install, reset
});
