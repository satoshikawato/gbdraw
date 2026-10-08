// @ts-check
import { getSessionResourceSource } from './file-content-cache.js';
import { depthFileSlotsFromValue } from './depth-track-state.js';
import { resolveDisambiguatedRecordSelection } from './record-options.js';

// A composite backing carries `descriptors` and no single `descriptor`.
const sourceFileIdentity = (file) => (
  /** @type {{ descriptor?: any } | null} */ (getSessionResourceSource(file))?.descriptor || file
);

export const groupLinearSourceRecords = (sequences) => {
  const groups = [];
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

export const linearSourceHasPrimaryInput = (source) => Boolean(
  source?.records?.some(({ sequence }) => sequence?.gb || sequence?.gff || sequence?.fasta)
);

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

/**
 * @typedef {{ allowed: false, reason: string, sourceIndex?: number }
 *   | { allowed: true, intent: string, sourceIndex: number, insertionIndex: number,
 *       sourceUid: any, recordCount: number, removedUids: any[], retainedSequences: any[] }
 * } LinearSourceRemovalPlan
 */

/** @returns {LinearSourceRemovalPlan} */
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

export const getLinearSourceDefaultDefinition = (source) => {
  return String(source?.records?.[0]?.sequence?.file_definition ?? source?.sequence?.file_definition ?? '');
};

export const setLinearSourceDefaultDefinition = (source, value) => {
  const text = String(value ?? '');
  (source?.records || []).forEach(({ sequence }) => {
    sequence.file_definition = text;
  });
  if (source?.sequence) {
    source.sequence.file_definition = text;
  }
};

export const getLinearSourceDefaultSubtitle = (source) => {
  return String(source?.records?.[0]?.sequence?.file_subtitle ?? source?.sequence?.file_subtitle ?? '');
};

export const setLinearSourceDefaultSubtitle = (source, value) => {
  const text = String(value ?? '');
  (source?.records || []).forEach(({ sequence }) => {
    sequence.file_subtitle = text;
  });
  if (source?.sequence) {
    source.sequence.file_subtitle = text;
  }
};

// A record draws its own Definition, else the File default the user entered,
// else the definition inferred from that record (D-12).
export const resolveLinearRecordEffectiveDefinition = (sequence, source = null) => {
  const own = String(sequence?.definition ?? '').trim();
  if (own) return sequence.definition;
  const fileDef = String(sequence?.file_definition ?? (source ? getLinearSourceDefaultDefinition(source) : '')).trim();
  return fileDef || String(sequence?.inferred_definition ?? '').trim();
};

// The definition inferred from the record a row selects; an automatic selector
// names a record only when its File has exactly one.
export const inferredDefinitionForRecord = (records, selector) => {
  const list = Array.isArray(records) ? records : [];
  const record = String(selector ?? '').trim()
    ? resolveDisambiguatedRecordSelection(list, selector).record
    : (list.length === 1 ? list[0] : null);
  return String(record?.inferredDefinition || '');
};

export const resolveLinearRecordEffectiveSubtitle = (sequence, source = null) => {
  const own = String(sequence?.record_subtitle ?? '').trim();
  if (own) return sequence.record_subtitle;
  const fileSub = String(sequence?.file_subtitle ?? (source ? getLinearSourceDefaultSubtitle(source) : '')).trim();
  return fileSub || '';
};

// The explicit LOSAT translation table of one record, or null for the default.
export const losatRecordGencode = (sequence) => {
  const raw = sequence?.losat_gencode;
  if (raw === null || raw === undefined || raw === '') return null;
  const value = Number(raw);
  return Number.isFinite(value) ? value : null;
};

// The single LOSAT job plan, used by Generate and by the Settings job-count
// estimate (PD-OI-018 revision 4, N-10). Record-pair evidence and source-file
// execution have different cardinalities: between two source files the query
// source searches the whole subject source. A record never searches a database
// that holds itself unless that self comparison was requested: a comparison
// within one source searches that source without the query record, and a
// requested self comparison searches the record alone. Explicit translation
// tables may require compatible subsets within one source.
export const planLosatSourceJobs = ({ sequences, specs, buildArgs }) => {
  const byRecord = new Map();
  groupLinearSourceRecords(sequences).forEach((source) => {
    source.records.forEach(({ index }) => byRecord.set(index, source));
  });
  const jobs = new Map();
  const bySpec = new Map();
  for (const spec of specs) {
    const { queryIndex, subjectIndex } = spec;
    const args = buildArgs(queryIndex, subjectIndex);
    const argsText = JSON.stringify(args);
    const querySource = byRecord.get(queryIndex);
    const subjectSource = byRecord.get(subjectIndex);
    const scope = queryIndex === subjectIndex ? 'self'
      : querySource === subjectSource ? 'within-source' : 'between-sources';
    const key = JSON.stringify([querySource.uid, subjectSource.uid, args, scope,
      ...(scope === 'between-sources' ? [] : [queryIndex])]);
    let job = jobs.get(key);
    if (!job) {
      const queryIndexes = scope === 'between-sources'
        ? querySource.records.map(({ index }) => index).filter(
            (index) => JSON.stringify(buildArgs(index, subjectIndex)) === argsText
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

export const prepareLosatSourceBatches = async ({
  sequences, specs, getEntry, buildArgs, hashText, protein
}) => {
  const plan = planLosatSourceJobs({ sequences, specs, buildArgs });
  const sides = new Map();
  const prepareSide = async (indexes) => {
    const ordered = [...indexes].sort((a, b) => String(sequences[a].uid).localeCompare(String(sequences[b].uid)));
    const key = ordered.join(',');
    if (sides.has(key)) return sides.get(key);
    const ids = new Map();
    const parts = [];
    for (const index of ordered) {
      const entry = await getEntry(index);
      let fasta = entry.fasta;
      const prefix = !protein && ordered.length > 1
        ? `n_${await hashText(String(sequences[index].uid))}` : '';
      let ordinal = 0;
      fasta = fasta.replace(/^>(\S+)([^\r\n]*)/gm, (_header, originalId, description) => {
        const id = prefix ? `${prefix}_${ordinal++}` : originalId;
        if (ids.has(id)) throw new Error(`Duplicate LOSAT source identifier: ${id}`);
        ids.set(id, { index, originalId });
        return `>${id}${description}`;
      });
      parts.push(fasta.endsWith('\n') ? fasta : `${fasta}\n`);
    }
    const fasta = parts.join('');
    const hash = await hashText(fasta);
    const side = { indexes: ordered, ids, fasta, sequenceKey: `source:${hash}`, hash };
    sides.set(key, side);
    return side;
  };
  const batches = [];
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
  return { batches, bySpec };
};

export const splitLosatSourceResult = (text, batch, jobs) => {
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
    text: rows.get(`${job.queryIndex}:${job.subjectIndex}`).map((line) => `${line}\n`).join('')
  }));
};
