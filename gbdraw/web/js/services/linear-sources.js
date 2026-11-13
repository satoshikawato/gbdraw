// @ts-check
/** @import { LinearSeq } from '../state.js' */
/** @import { InputFileRef } from './history-files.js' */
import { inputTypeSpec } from './input-types.js';
/** @import { InputSlot } from './input-types.js' */
import { depthFileSlotsFromValue } from './depth-track-state.js';
import { resolveDisambiguatedRecordSelection } from './record-options.js';
import { plainTextLinearRecordLabel } from './linear-comparisons.js';
import { LOSAT_GENETIC_CODES, LOSAT_GENETIC_CODE_TEXT } from '../contracts/schema.generated.js';
import { diagnosticError } from '../utils/error-normalization.js';

// Two Session resource views of one stored resource are one source (`sourceId`).
/** @param {InputFileRef | null | undefined} file */
const sourceFileIdentity = (file) => file?.sourceId || null;

/** The texts of a record card, as a draft or a request row holds them. @typedef {Partial<Pick<LinearSeq, 'definition' | 'record_subtitle' | 'file_definition' | 'file_subtitle' | 'inferred_definition'>>} LinearSeqTexts */
/**
 * An entry of a File's `records`: a record card and its position in the grouped list.
 * @template {Partial<LinearSeq>} [S=LinearSeq]
 * @typedef {{ sequence: S, index: number }} LinearSourceRecord
 */
/**
 * A File's record cards: the cards that read the same files.
 * @template {Partial<LinearSeq>} [S=LinearSeq]
 * @typedef {{ uid: S['uid'], sequence: S, index: number, files: unknown[], records: LinearSourceRecord<S>[] }} LinearSourceGroup
 */
/**
 * What the File accessors read of a group.
 * @typedef {{ records?: readonly LinearSourceRecord<Partial<LinearSeq>>[], sequence?: Partial<LinearSeq> | null }} LinearSourceView
 */

/**
 * @template {Partial<LinearSeq>} S
 * @param {readonly S[]} sequences
 * @returns {LinearSourceGroup<S>[]}
 */
export const groupLinearSourceRecords = (sequences) => {
  /** @type {LinearSourceGroup<S>[]} */
  const groups = [];
  /** @type {Map<unknown, LinearSourceGroup<S>[]>} */
  const bySource = new Map();
  sequences.forEach((sequence, index) => {
    const files = [sequence.gb, sequence.gff, sequence.fasta].map(sourceFileIdentity);
    const key = files.find(Boolean) || sequence.uid;
    const candidates = bySource.get(key) || [];
    let group = candidates.find((entry) => entry.files.every((file, position) => file === files[position]));
    if (!group) {
      group = { uid: sequence.uid, sequence, index, files, records: [] };
      candidates.push(group);
      bySource.set(key, candidates);
      groups.push(group);
    }
    group.records.push({ sequence, index });
  });
  return groups;
};

/** @param {LinearSourceView | null | undefined} source */
export const linearSourceHasPrimaryInput = (source) => Boolean(
  source?.records?.some(({ sequence }) => sequence?.gb || sequence?.gff || sequence?.fasta)
);

/** @param {LinearSourceView | null | undefined} source */
export const isPristineLinearSource = (source) => {
  const records = Array.isArray(source?.records) ? source.records : [];
  if (records.length !== 1) return false;
  const sequence = records[0]?.sequence;
  if (!sequence || linearSourceHasPrimaryInput(source)) return false;
  const depth = depthFileSlotsFromValue(sequence.depth);
  return (
    depth.every((file) => !file) &&
    String(sequence.definition || '') === '' &&
    String(sequence.record_subtitle || '') === '' &&
    String(sequence.file_definition || '') === '' &&
    String(sequence.file_subtitle || '') === '' &&
    String(sequence.region_record_id || '') === '' &&
    sequence.region_start == null &&
    sequence.region_end == null &&
    !sequence.region_reverse &&
    !(Array.isArray(sequence.windows) && sequence.windows.length)
  );
};

// The Linear File cards shown folded to their header line, by source uid, so
// the state follows a card through Move up and Move down. It is view state:
// no Session, History step, Undo, or Reset Settings holds or changes it.
/** @typedef {Set<string>} FoldedLinearSources */
/** @typedef {LinearSourceView & { uid: string }} FoldableLinearSource */

/**
 * What a card's disclosure reports: the user's press on the header line, or a
 * path that opened the card to show a setting.
 * @param {FoldedLinearSources} folded
 * @param {string} uid
 * @param {boolean} isFolded
 */
export const setLinearSourceFolded = (folded, uid, isFolded) => {
  if (isFolded) folded.add(uid);
  else folded.delete(uid);
};

/**
 * The fold state a Session load that was rolled back returns to.
 * @param {FoldedLinearSources} folded
 * @param {readonly string[]} uids
 */
export const restoreFoldedLinearSources = (folded, uids) => {
  folded.clear();
  uids.forEach((uid) => folded.add(uid));
};

/** @typedef {InputSlot} LinearFileField */

/**
 * The file slots a File card shows for an input type, in the order of its
 * uploaders.
 * @param {string} inputType
 * @returns {readonly LinearFileField[]}
 */
export const linearCardFileFields = (inputType) => {
  const { primary, paired } = inputTypeSpec(inputType);
  return paired ? [primary, paired] : [primary];
};

/**
 * A card holds its input when every slot it shows has a file. The folded line
 * and the load rule read this one answer.
 * @param {LinearSourceView} source
 * @param {string} inputType
 */
export const linearCardInputComplete = (source, inputType) => linearCardFileFields(inputType)
  .every((field) => Boolean(source.sequence?.[field]));

/**
 * The one rule of which cards a loaded Session shows folded: those that hold
 * their input. A card without a file, or with a slot still empty, stays open.
 * @param {FoldedLinearSources} folded
 * @param {readonly FoldableLinearSource[]} sources
 * @param {string} inputType
 */
export const foldLoadedLinearSources = (folded, sources, inputType) => restoreFoldedLinearSources(
  folded, sources.filter((source) => linearCardInputComplete(source, inputType)).map(({ uid }) => uid)
);

/**
 * @typedef {{ allowed: false, reason: string, sourceIndex?: number }
 *   | { allowed: true, intent: string, sourceIndex: number, insertionIndex: number,
 *       sourceUid: string, recordCount: number, removedUids: string[], retainedSequences: LinearSeq[] }
 * } LinearSourceRemovalPlan
 */

/**
 * @param {{ sequences: readonly LinearSeq[] | null | undefined, sourceUid: unknown, intent: string }} options
 * @returns {LinearSourceRemovalPlan}
 */
export const planLinearSourceRemoval = ({ sequences, sourceUid, intent }) => {
  const ordered = Array.from(sequences || []);
  const groups = groupLinearSourceRecords(ordered);
  const sourceIndex = groups.findIndex((source) => source.uid === sourceUid);
  if (sourceIndex < 0) return { allowed: false, reason: 'missing-source' };
  if (!['clear', 'delete'].includes(intent)) return { allowed: false, reason: 'invalid-intent' };
  if (intent === 'delete' && groups.length <= 1) {
    return { allowed: false, reason: 'sole-source', sourceIndex };
  }

  const source = groups[sourceIndex];
  const removedUids = new Set(source.records.map(({ sequence }) => sequence.uid));
  return {
    allowed: true,
    intent,
    sourceIndex,
    insertionIndex: source.index,
    sourceUid: source.uid,
    recordCount: source.records.length,
    removedUids: [...removedUids],
    retainedSequences: ordered.filter((sequence) => !removedUids.has(sequence.uid))
  };
};

/**
 * @param {LinearSourceView | null | undefined} source
 * @param {unknown} trackIndex
 */
export const linearSourceDepthStatus = (source, trackIndex) => {
  const index = Number(trackIndex);
  const records = Array.isArray(source?.records) ? source.records : [];
  const files = records.map(({ sequence }) => (
    Number.isInteger(index) && index >= 0
      ? depthFileSlotsFromValue(sequence?.depth)[index] || null
      : null
  ));
  const selected = files.filter(Boolean);
  const commonIdentity = selected.length === records.length && selected.length > 0
    ? sourceFileIdentity(selected[0])
    : null;
  const common = commonIdentity && selected.every((file) => sourceFileIdentity(file) === commonIdentity);
  return {
    state: selected.length === 0 ? 'empty' : common ? 'common' : 'mixed',
    file: common ? selected[0] : null,
    selectedCount: selected.length,
    recordCount: records.length
  };
};

/**
 * @template {Partial<LinearSeq>} S
 * @param {readonly S[] | null | undefined} sequences
 * @param {number} groupIndex
 * @param {number} direction
 * @returns {S[]}
 */
export const moveLinearSourceGroup = (sequences, groupIndex, direction) => {
  const ordered = Array.from(sequences || []);
  const index = groupIndex;
  const offset = direction;
  if (!Number.isInteger(index) || ![-1, 1].includes(offset)) return ordered;

  const groups = groupLinearSourceRecords(ordered);
  const target = index + offset;
  if (index < 0 || index >= groups.length || target < 0 || target >= groups.length) return ordered;

  [groups[index], groups[target]] = [groups[target], groups[index]];
  return groups.flatMap(({ records }) => records.map(({ sequence }) => sequence));
};

/** @param {LinearSourceView | null | undefined} source */
export const getLinearSourceDefaultDefinition = (source) => {
  return String(source?.records?.[0]?.sequence?.file_definition ?? source?.sequence?.file_definition ?? '');
};

/**
 * @param {LinearSourceView | null | undefined} source
 * @param {unknown} value
 */
export const setLinearSourceDefaultDefinition = (source, value) => {
  const text = String(value ?? '');
  (source?.records || []).forEach(({ sequence }) => {
    sequence.file_definition = text;
  });
  if (source?.sequence) {
    source.sequence.file_definition = text;
  }
};

/** @param {LinearSourceView | null | undefined} source */
export const getLinearSourceDefaultSubtitle = (source) => {
  return String(source?.records?.[0]?.sequence?.file_subtitle ?? source?.sequence?.file_subtitle ?? '');
};

/**
 * @param {LinearSourceView | null | undefined} source
 * @param {unknown} value
 */
export const setLinearSourceDefaultSubtitle = (source, value) => {
  const text = String(value ?? '');
  (source?.records || []).forEach(({ sequence }) => {
    sequence.file_subtitle = text;
  });
  if (source?.sequence) {
    source.sequence.file_subtitle = text;
  }
};

// The one decision of which name fields a File card shows: one block. A File
// with one row and no File-level value shows that record's own names on the
// card, and this returns the record. Any other File shows the File-level
// block, so no stored value is hidden, and each record's names in its Record
// options; this returns null.
// It reads the rows the File has, not whether its records were inspected.
/** @param {{ records: readonly LinearSourceRecord<Partial<LinearSeq>>[] }} source @returns {LinearSourceRecord<Partial<LinearSeq>> | null} */
export const linearFileCardRecord = (source) => {
  const fileValue = Boolean(getLinearSourceDefaultDefinition(source) || getLinearSourceDefaultSubtitle(source));
  return source.records.length === 1 && !fileValue ? source.records[0] : null;
};

// A record draws its own Definition, else the File default the user entered,
// else the definition inferred from that record (D-12).
/**
 * @param {LinearSeqTexts} sequence
 * @param {LinearSourceView | null} [source]
 */
export const resolveLinearRecordEffectiveDefinition = (sequence, source = null) => {
  const own = String(sequence?.definition ?? '').trim();
  if (own) return sequence.definition;
  const fileDef = String(sequence?.file_definition ?? (source ? getLinearSourceDefaultDefinition(source) : '')).trim();
  return fileDef || String(sequence?.inferred_definition ?? '').trim();
};

// How panels and lists name each record (OV-436): its drawn Definition as
// plain text, else the Definition another record of its File draws, else its
// accession or file name. A File that supplies several records adds each
// record's accession, so two chromosomes of one genome stay apart. The
// diagram's own record headings are not changed.
/**
 * @param {readonly Partial<LinearSeq>[]} sequences
 * @returns {Map<string, { definition: string, name: string, accession: string, label: string }>}
 */
export const linearRecordNames = (sequences) => {
  const names = new Map();
  const plain = (/** @type {unknown} */ value) => (String(value ?? '').trim() ? plainTextLinearRecordLabel(value) : '');
  groupLinearSourceRecords(sequences).forEach((source) => {
    const definitions = source.records.map((/** @type {LinearSourceRecord<Partial<LinearSeq>>} */ { sequence }) => plain(resolveLinearRecordEffectiveDefinition(sequence, source)));
    const fileDefinition = definitions.find(Boolean) || '';
    source.records.forEach((/** @type {LinearSourceRecord<Partial<LinearSeq>>} */ { sequence }, /** @type {number} */ position) => {
      const recordId = String(sequence?.region_record_id ?? '').trim();
      const definition = definitions[position] || fileDefinition;
      const name = definition || recordId || plain(sequence?.gb?.name || sequence?.gff?.name || sequence?.fasta?.name);
      const accession = source.records.length > 1 && recordId !== name ? recordId : '';
      names.set(String(sequence?.uid ?? ''), { definition, name, accession, label: [name, accession].filter(Boolean).join(' · ') || 'Record' });
    });
  });
  return names;
};

// The definition inferred from the record a row selects; an automatic selector
// names a record only when its File has exactly one.
/**
 * @param {readonly { inferredDefinition?: unknown, sourceIndex?: unknown, selector?: unknown, recordId?: unknown }[] | null | undefined} records
 * @param {unknown} selector
 */
export const inferredDefinitionForRecord = (records, selector) => {
  const list = Array.isArray(records) ? records : [];
  const record = String(selector ?? '').trim()
    ? resolveDisambiguatedRecordSelection(list, selector).record
    : (list.length === 1 ? list[0] : null);
  return String(record?.inferredDefinition || '');
};

/**
 * @param {LinearSeqTexts} sequence
 * @param {LinearSourceView | null} [source]
 */
export const resolveLinearRecordEffectiveSubtitle = (sequence, source = null) => {
  const own = String(sequence?.record_subtitle ?? '').trim();
  if (own) return sequence.record_subtitle;
  const fileSub = String(sequence?.file_subtitle ?? (source ? getLinearSourceDefaultSubtitle(source) : '')).trim();
  return fileSub || '';
};

/**
 * Whether LOSAT reads a typed TLOSATX genetic code: one of the NCBI tables in Python's `LOSAT_GENETIC_CODES`. The one
 * rule for a typed code, of Generate's refusal (`requireLosatGeneticCode`) and of the File card's mark
 * (`losatGencodeProblem`, OV-3326).
 * @param {number} code
 */
export const losatReadsGeneticCode = (code) => LOSAT_GENETIC_CODES.includes(code);

/**
 * Where a TLOSATX genetic code was typed: a Linear card's LOSAT Gencode (its number among the cards, and its identity,
 * which the refusal's action focuses), a comparison ring's Comparison gencode (its number among the rings), or the
 * Circular Reference gencode.
 * @typedef {{ field: 'losat_gencode', inputOrdinal: number, cardUid: string }
 *   | { field: 'conservation_losat_gencodes', inputOrdinal: number }
 *   | { field: 'reference_gencode' }} LosatGeneticCodeField
 */

/**
 * A TLOSATX genetic code a search passes. A code that is not one of the NCBI tables LOSAT reads (Python's
 * `LOSAT_GENETIC_CODES`) stops LOSAT with a trap, so Generate refuses it before any search, naming the field it was
 * typed in (OV-3326).
 * @param {number} code
 * @param {LosatGeneticCodeField} at
 * @returns {number}
 */
export const requireLosatGeneticCode = (code, at) => {
  if (losatReadsGeneticCode(code)) return code;
  throw diagnosticError('COMPARISON_INPUT', { ...at, reason: 'GENETIC_CODE', geneticCodes: LOSAT_GENETIC_CODE_TEXT });
};

/**
 * The explicit LOSAT translation table typed on one Linear card, or null for the default.
 * @param {{ losat_gencode?: unknown } | null | undefined} sequence
 * @returns {number | null}
 */
export const losatRecordGencode = (sequence) => {
  const raw = sequence?.losat_gencode;
  if (raw === null || raw === undefined || raw === '') return null;
  const value = Number(raw);
  return Number.isFinite(value) ? value : null;
};

/**
 * What a Linear card's LOSAT Gencode lacks under TLOSATX, as the File card marks it: `'empty'` when no table is typed
 * (the search passes none, so LOSAT uses its default table 1; the card asks for one), `'unread'` when the typed table
 * is not one LOSAT reads (`losatReadsGeneticCode`; Generate refuses it), else null.
 * @param {{ losat_gencode?: unknown } | null | undefined} sequence
 * @returns {'empty' | 'unread' | null}
 */
export const losatGencodeProblem = (sequence) => {
  const code = losatRecordGencode(sequence);
  if (code === null) return 'empty';
  return losatReadsGeneticCode(code) ? null : 'unread';
};

/**
 * The translation table a TLOSATX search passes for one Linear card, or null for the default; a table LOSAT does not
 * read is refused (`requireLosatGeneticCode`), naming the card.
 * @param {{ losat_gencode?: unknown } | null | undefined} sequence
 * @param {{ inputOrdinal: number, cardUid: string }} card The card's number among the drawing's cards, which the
 *   refusal names, and its identity, whose LOSAT Gencode the refusal's action focuses.
 * @returns {number | null}
 */
export const searchedLosatGencode = (sequence, card) => {
  const code = losatRecordGencode(sequence);
  return code === null ? null : requireLosatGeneticCode(code, { field: 'losat_gencode', ...card });
};

/**
 * A search's arguments as the raw cache spells them (`LOSAT_SEARCH_FLAGS`), in Python's order (`losat_cache_args`):
 * each flag with a value, none without one.
 * @param {readonly [string, string | number | null][]} options
 * @returns {string[]}
 */
const losatArgs = (options) => options.flatMap(([flag, value]) => (value === null || value === '' ? [] : [flag, String(value)]));

/**
 * The arguments of a LOSATN or TLOSATX search: LOSATN's task, TLOSATX's translation tables (none for a default
 * table), and the search space of the genome pair when Python plans one (`plan_genome_scales`, OV-3143). Generate and
 * the job estimate build them here. Only TLOSATX reads the tables (`gencodes`, query then database): LOSATN reads none,
 * so a Gencode left on a card (hidden under LOSATN) is never read or refused for it.
 * @param {{ program: 'blastn' | 'tblastx', task: string, gencodes: () => readonly [number | null, number | null],
 *   searchSpace: string | null }} search
 * @returns {string[]}
 */
export const nucleotideLosatArgs = ({ program, task, gencodes, searchSpace }) => {
  if (program === 'blastn') return losatArgs([['--task', task], ['--searchsp', searchSpace]]);
  const [queryGencode, dbGencode] = gencodes();
  return losatArgs([['--query-gencode', queryGencode], ['--db-gencode', dbGencode], ['--searchsp', searchSpace]]);
};

/**
 * The arguments of a LOSATP search: one HSP per subject protein for Pairwise, the candidate limit, and the proteins of
 * the subject genome as the database size when Python plans one (`plan_genome_scales`, PD-OI-018 revision 6).
 * Generate and the job estimate build them here.
 * @param {{ pairwise: boolean, candidateLimit: number | null, databaseSize: string | null }} search
 * @returns {string[]}
 */
export const proteinLosatArgs = ({ pairwise, candidateLimit, databaseSize }) => losatArgs([
  ['--max-hsps-per-subject', pairwise ? 1 : null], ['--max-target-seqs', candidateLimit], ['--dbsize', databaseSize]
]);

/** @typedef {'self' | 'within-source' | 'between-sources'} LosatSourceScope */
/**
 * One search of the plan: the record indexes of the query and subject sides and the record-pair specs it answers.
 * @template Spec
 * @typedef {{ queryIndexes: number[], subjectIndexes: number[], args: readonly unknown[], scope: LosatSourceScope, specs: Spec[] }} LosatSourceJob
 */

// The single LOSAT job plan, used by Generate and by the Settings job-count
// estimate (PD-OI-018 revision 4, N-10). Record-pair evidence and source-file
// execution have different cardinalities: between two source files the query
// source searches the whole subject source. A record never searches a database
// that holds itself unless that self comparison was requested: a comparison
// within one source searches that source without the query record, and a
// requested self comparison searches the record alone. Explicit translation
// tables may require compatible subsets within one source.
// With records turned OFF (D-07), `sequences` lists the drawn records first,
// whose indexes the specs use, then the OFF records (`losatSearchSequences`):
// an OFF record is never a query, but stays in its source's database, so
// turning a record ON or OFF changes no other record's E-value. An OFF record
// has no row and so no genome: it joins a database under the arguments of the
// drawn subject it is searched with (`buildArgs`'s third index), so the search
// space of that subject's genome pair (OV-3143) does not keep it out.
/**
 * @template {Partial<LinearSeq>} S
 * @template {{ queryIndex: number, subjectIndex: number }} Spec
 * @param {{
 *   sequences: readonly S[], specs: readonly Spec[],
 *   buildArgs: (queryIndex: number, subjectIndex: number, drawnSubjectIndex?: number) => readonly unknown[],
 *   drawnCount?: number
 * }} options `buildArgs(query, subject, drawnSubject)` gives a search's arguments: those `subject` takes as a member
 *   of a database searched with the drawn record `drawnSubject` (by default `subject` itself).
 * @returns {{ jobs: LosatSourceJob<Spec>[], bySpec: Map<Spec, LosatSourceJob<Spec>> }}
 */
export const planLosatSourceJobs = ({ sequences, specs, buildArgs, drawnCount = sequences.length }) => {
  /** @type {Map<number, LinearSourceGroup<S>>} */
  const byRecord = new Map();
  groupLinearSourceRecords(sequences).forEach((source) => {
    source.records.forEach(({ index }) => byRecord.set(index, source));
  });
  /** @type {Map<string, LosatSourceJob<Spec>>} */
  const jobs = new Map();
  /** @type {Map<Spec, LosatSourceJob<Spec>>} */
  const bySpec = new Map();
  for (const spec of specs) {
    const { queryIndex, subjectIndex } = spec;
    const args = buildArgs(queryIndex, subjectIndex);
    const argsText = JSON.stringify(args);
    // A spec addresses a listed record: `buildLosatJobSpecs` numbers them within the listed sequences.
    const querySource = /** @type {LinearSourceGroup<S>} */ (byRecord.get(queryIndex));
    const subjectSource = /** @type {LinearSourceGroup<S>} */ (byRecord.get(subjectIndex));
    /** @type {LosatSourceScope} */
    const scope = queryIndex === subjectIndex ? 'self'
      : querySource === subjectSource ? 'within-source' : 'between-sources';
    const key = JSON.stringify([querySource.uid, subjectSource.uid, args, scope,
      ...(scope === 'between-sources' ? [] : [queryIndex])]);
    let job = jobs.get(key);
    if (!job) {
      const queryIndexes = scope === 'between-sources'
        ? querySource.records.map(({ index }) => index).filter(
            (index) => index < drawnCount && JSON.stringify(buildArgs(index, subjectIndex)) === argsText
          )
        : [queryIndex];
      const subjectIndexes = scope === 'self'
        ? [subjectIndex]
        : subjectSource.records.map(({ index }) => index).filter((index) => (
            index !== queryIndex
            && JSON.stringify(index < drawnCount ? buildArgs(queryIndex, index) : buildArgs(queryIndex, index, subjectIndex))
              === argsText
          ));
      job = { queryIndexes, subjectIndexes, args, scope, specs: [] };
      jobs.set(key, job);
    }
    job.specs.push(spec);
    bySpec.set(spec, job);
  }
  return { jobs: [...jobs.values()], bySpec };
};

/**
 * The records a LOSAT search reads: the drawn records, whose positions the
 * request and the specs use, then the OFF records, which are database members
 * only (`planLosatSourceJobs`).
 * @template {{ uid?: unknown }} T
 * @param {readonly T[]} drawn The drawn records: the drawn units of a Generate (`linearUnitSearch`).
 * @param {readonly T[]} sequences Every record card, in file order.
 * @param {readonly T[]} [drawnCards] The cards of the drawn records, when `drawn` holds units made from them.
 * @returns {{ sequences: T[], drawnCount: number }}
 */
export const losatSearchSequences = (drawn, sequences, drawnCards = drawn) => {
  const drawnSet = new Set(drawnCards);
  return { sequences: [...drawn, ...sequences.filter((sequence) => !drawnSet.has(sequence))], drawnCount: drawn.length };
};

/** @typedef {{ indexes: number[], ids: Map<string, { index: number, originalId: string }>, fasta: string, sequenceKey: string, hash: string }} LosatSourceSide A side's records as one FASTA, and the identifiers it was renamed to. */
/**
 * A batch of the plan. `searchable` says both sides hold a sequence: a side of records without a CDS protein holds
 * none, so its pairs find no hit and LOSAT, which refuses an empty database, does not run (OV-3493; Python
 * `LosatBatch.searchable`, held equal by the shared job-plan vectors).
 * @template Spec
 * @typedef {{ query: LosatSourceSide, subject: LosatSourceSide, searchContext: string | null, args: readonly unknown[], scope: LosatSourceScope, specs: Spec[], searchable: boolean }} LosatSourceBatch
 */
/**
 * A side's renamed FASTA and digest, as `prepareLosatSourceBatches` built it from these record entries with this
 * digest function (OV-1565): Generate reuses the record entries it extracted before, so an unchanged side is not
 * rebuilt or rehashed. `ids` lists each renamed identifier with its record's position in the side and original ID.
 * @typedef {{ entries: readonly object[], hashText: (text: string) => Promise<string>, fasta: string, hash: string,
 *   ids: readonly (readonly [string, number, string])[] }} PreparedLosatSide
 */
/**
 * The sides the last call prepared, by kind and card uids (the other inputs of the renamed FASTA). Only the last
 * call's sides stay, so the hold is one Generate's sides whatever their count.
 * @type {Map<string, PreparedLosatSide>}
 */
let lastPreparedLosatSides = new Map();

/**
 * @template {Partial<LinearSeq>} S
 * @template {{ queryIndex: number, subjectIndex: number }} Spec
 * @param {{
 *   sequences: readonly S[], specs: readonly Spec[],
 *   getEntry: (index: number) => Promise<{ fasta: string }>,
 *   getDigest: (index: number) => Promise<string>,
 *   buildArgs: (queryIndex: number, subjectIndex: number, drawnSubjectIndex?: number) => readonly unknown[],
 *   hashText: (text: string) => Promise<string>, protein?: boolean, drawnCount?: number
 * }} options `getDigest(index)` is `hashText` of `getEntry(index).fasta`, read from its owner: a single-record side
 *   whose FASTA is its entry's text takes that digest instead of hashing the text again (OV-1598).
 * @returns {Promise<{ batches: LosatSourceBatch<Spec>[], bySpec: Map<Spec, LosatSourceBatch<Spec>> }>}
 */
export const prepareLosatSourceBatches = async ({
  sequences, specs, getEntry, getDigest, buildArgs, hashText, protein, drawnCount = sequences.length
}) => {
  const plan = planLosatSourceJobs({ sequences, specs, buildArgs, drawnCount });
  /** @type {Map<string, LosatSourceSide>} */
  const sides = new Map();
  /** @type {Map<string, PreparedLosatSide>} */
  const preparedSides = new Map();
  /** @param {number[]} indexes */
  const prepareSide = async (indexes) => {
    const ordered = [...indexes].sort((a, b) => String(sequences[a].uid).localeCompare(String(sequences[b].uid)));
    const key = ordered.join(',');
    // `sides` holds the sides prepared under a key.
    if (sides.has(key)) return /** @type {LosatSourceSide} */ (sides.get(key));
    /** @type {{ fasta: string }[]} */
    const entries = [];
    for (const index of ordered) entries.push(await getEntry(index));
    const uids = ordered.map((index) => String(sequences[index].uid));
    const preparedKey = JSON.stringify([Boolean(protein), uids]);
    let prepared = preparedSides.get(preparedKey) || lastPreparedLosatSides.get(preparedKey) || null;
    if (!prepared || prepared.hashText !== hashText
      || prepared.entries.some((entry, position) => entry !== entries[position])) {
      /** @type {Map<string, [string, number, string]>} */
      const renamed = new Map();
      /** @type {string[]} */
      const parts = [];
      for (const [position, entry] of entries.entries()) {
        const prefix = !protein && ordered.length > 1 ? `n_${await hashText(uids[position])}` : '';
        let ordinal = 0;
        const fasta = entry.fasta.replace(/^>(\S+)([^\r\n]*)/gm, (/** @type {string} */ _header, /** @type {string} */ originalId, /** @type {string} */ description) => {
          const id = prefix ? `${prefix}_${ordinal++}` : originalId;
          if (renamed.has(id)) throw new Error(`Duplicate LOSAT source identifier: ${id}`);
          renamed.set(id, [id, position, originalId]);
          return `>${id}${description}`;
        });
        parts.push(fasta.endsWith('\n') ? fasta : `${fasta}\n`);
      }
      const fasta = parts.join('');
      const hash = entries.length === 1 && fasta === entries[0].fasta ? await getDigest(ordered[0]) : await hashText(fasta);
      prepared = { entries, hashText, fasta, hash, ids: [...renamed.values()] };
    }
    preparedSides.set(preparedKey, prepared);
    const side = {
      indexes: ordered,
      ids: new Map(prepared.ids.map(([id, position, originalId]) => [id, { index: ordered[position], originalId }])),
      fasta: prepared.fasta,
      sequenceKey: `source:${prepared.hash}`,
      hash: prepared.hash
    };
    sides.set(key, side);
    return side;
  };
  /** @type {LosatSourceBatch<Spec>[]} */
  const batches = [];
  /** @type {Map<Spec, LosatSourceBatch<Spec>>} */
  const bySpec = new Map();
  for (const job of plan.jobs) {
    const query = await prepareSide(job.queryIndexes);
    const subject = await prepareSide(job.subjectIndexes);
    // The searched database is part of raw-cache identity.
    const searchContext = query.indexes.length > 1 || subject.indexes.length > 1
      ? await hashText(JSON.stringify([query.hash, subject.hash])) : null;
    const searchable = query.ids.size > 0 && subject.ids.size > 0;
    const batch = { query, subject, searchContext, args: job.args, scope: job.scope, specs: job.specs, searchable };
    batches.push(batch);
    job.specs.forEach((spec) => bySpec.set(spec, batch));
  }
  lastPreparedLosatSides = preparedSides;
  return { batches, bySpec };
};

/**
 * Runs the source jobs whose batch is searchable and answers each of the others with no rows. The runner gets each
 * searched job without its source-job bookkeeping (`batch`, `members`).
 * @template {{ cacheKey: string, batch: { searchable: boolean }, members?: unknown }} J
 * @template O
 * @template {{ cacheKey?: string, text: string }} R
 * @param {readonly J[]} jobs
 * @param {(searched: Omit<J, 'batch' | 'members'>[], options: O) => Promise<R[]>} run The search of the searchable jobs.
 * @param {O} options The runner's options.
 * @returns {Promise<Array<R | { cacheKey: string, text: string }>>}
 */
export const runSearchableLosatSourceJobs = async (jobs, run, options) => {
  const searched = jobs.filter((job) => job.batch.searchable).map(({ batch, members, ...job }) => job);
  const answered = jobs.filter((job) => !job.batch.searchable).map(({ cacheKey }) => ({ cacheKey, text: '' }));
  return [...answered, ...(searched.length > 0 ? await run(searched, options) : [])];
};

/**
 * @template K
 * @param {string} text
 * @param {{ query: LosatSourceSide, subject: LosatSourceSide }} batch
 * @param {readonly { queryIndex: number, subjectIndex: number, cacheKey: K }[]} jobs
 * @returns {{ cacheKey: K, text: string }[]}
 */
export const splitLosatSourceResult = (text, batch, jobs) => {
  /** @type {Map<string, string[]>} */
  const rows = new Map(jobs.map((job) => [`${job.queryIndex}:${job.subjectIndex}`, []]));
  for (const line of String(text).split(/\r?\n/)) {
    if (!line.trim() || line.trimStart().startsWith('#')) continue;
    const columns = line.split('\t');
    const query = batch.query.ids.get(columns[0]);
    const subject = batch.subject.ids.get(columns[1]);
    if (columns.length !== 12 || !query || !subject) {
      throw new Error('LOSAT source result contains malformed or unrecognized record endpoints.');
    }
    const target = rows.get(`${query.index}:${subject.index}`);
    if (target) {
      columns[0] = query.originalId;
      columns[1] = subject.originalId;
      target.push(columns.join('\t'));
    }
  }
  return jobs.map((job) => ({
    cacheKey: job.cacheKey,
    // `rows` holds an entry for every job.
    text: /** @type {string[]} */ (rows.get(`${job.queryIndex}:${job.subjectIndex}`)).map((line) => `${line}\n`).join('')
  }));
};
