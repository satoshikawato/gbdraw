// @ts-check
/** @import { LinearSeq } from '../state.js' */
/** @import { InputFileRef } from './history-files.js' */
import { depthFileSlotsFromValue } from './depth-track-state.js';
import { resolveDisambiguatedRecordSelection } from './record-options.js';
import { plainTextLinearRecordLabel } from './linear-comparisons.js';

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
    !sequence.region_reverse
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

/** @typedef {'gb' | 'gff' | 'fasta'} LinearFileField */

/**
 * The file slots a File card shows for an input type, in the order of its
 * uploaders.
 * @param {string} inputType
 * @returns {readonly LinearFileField[]}
 */
export const linearCardFileFields = (inputType) => (inputType === 'gff' ? ['gff', 'fasta'] : ['gb']);

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

// The explicit LOSAT translation table of one record, or null for the default.
/** @param {{ losat_gencode?: unknown } | null | undefined} sequence */
export const losatRecordGencode = (sequence) => {
  const raw = sequence?.losat_gencode;
  if (raw === null || raw === undefined || raw === '') return null;
  const value = Number(raw);
  return Number.isFinite(value) ? value : null;
};

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
// turning a record ON or OFF changes no other record's E-value.
/**
 * @template {Partial<LinearSeq>} S
 * @template {{ queryIndex: number, subjectIndex: number }} Spec
 * @param {{
 *   sequences: readonly S[], specs: readonly Spec[],
 *   buildArgs: (queryIndex: number, subjectIndex: number) => readonly unknown[], drawnCount?: number
 * }} options
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
            index !== queryIndex && JSON.stringify(buildArgs(queryIndex, index)) === argsText
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
 * @param {readonly T[]} drawn
 * @param {readonly T[]} sequences Every record card, in file order.
 * @returns {{ sequences: T[], drawnCount: number }}
 */
export const losatSearchSequences = (drawn, sequences) => {
  const drawnSet = new Set(drawn);
  return { sequences: [...drawn, ...sequences.filter((sequence) => !drawnSet.has(sequence))], drawnCount: drawn.length };
};

/** @typedef {{ indexes: number[], ids: Map<string, { index: number, originalId: string }>, fasta: string, sequenceKey: string, hash: string }} LosatSourceSide A side's records as one FASTA, and the identifiers it was renamed to. */
/**
 * @template Spec
 * @typedef {{ query: LosatSourceSide, subject: LosatSourceSide, searchContext: string | null, args: readonly unknown[], scope: LosatSourceScope, specs: Spec[] }} LosatSourceBatch
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
 *   buildArgs: (queryIndex: number, subjectIndex: number) => readonly unknown[],
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
    const batch = { query, subject, searchContext, args: job.args, scope: job.scope, specs: job.specs };
    batches.push(batch);
    job.specs.forEach((spec) => bySpec.set(spec, batch));
  }
  lastPreparedLosatSides = preparedSides;
  return { batches, bySpec };
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
