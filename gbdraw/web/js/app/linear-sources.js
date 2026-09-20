import { getSessionResourceSource } from '../services/file-content-cache.js';

export const groupLinearSourceRecords = (sequences) => {
  const groups = [];
  const bySource = new Map();
  const identity = (file) => getSessionResourceSource(file)?.descriptor || file;
  sequences.forEach((sequence, index) => {
    const files = [sequence.gb, sequence.gff, sequence.fasta].map(identity);
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

export const resolveLinearRecordEffectiveDefinition = (sequence, source = null) => {
  const own = String(sequence?.definition ?? '').trim();
  if (own) return sequence.definition;
  const fileDef = String(sequence?.file_definition ?? (source ? getLinearSourceDefaultDefinition(source) : '')).trim();
  return fileDef || '';
};

export const resolveLinearRecordEffectiveSubtitle = (sequence, source = null) => {
  const own = String(sequence?.record_subtitle ?? '').trim();
  if (own) return sequence.record_subtitle;
  const fileSub = String(sequence?.file_subtitle ?? (source ? getLinearSourceDefaultSubtitle(source) : '')).trim();
  return fileSub || '';
};

// Record-pair evidence and source-file execution have different cardinalities.
// Explicit translation tables may require compatible subsets within one source.
export const prepareLosatSourceBatches = async ({
  sequences, specs, getEntry, buildArgs, hashText, protein, excludeSelfComparisons = false
}) => {
  const byRecord = new Map();
  groupLinearSourceRecords(sequences).forEach((source) => {
    source.records.forEach(({ index }) => byRecord.set(index, source));
  });
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
  const batches = new Map();
  const bySpec = new Map();
  for (const spec of specs) {
    const args = buildArgs(spec.queryIndex, spec.subjectIndex);
    const argsText = JSON.stringify(args);
    const querySource = byRecord.get(spec.queryIndex);
    const subjectSource = byRecord.get(spec.subjectIndex);
    const separateRecords = excludeSelfComparisons && querySource === subjectSource;
    const key = JSON.stringify([querySource.uid, subjectSource.uid, args,
      ...(separateRecords ? [spec.queryIndex, spec.subjectIndex] : [])]);
    let batch = batches.get(key);
    if (!batch) {
      const query = await prepareSide((separateRecords ? [spec.queryIndex] : querySource.records.map(({ index }) => index)).filter(
        (index) => JSON.stringify(buildArgs(index, spec.subjectIndex)) === argsText
      ));
      const subject = await prepareSide((separateRecords ? [spec.subjectIndex] : subjectSource.records.map(({ index }) => index)).filter(
        (index) => JSON.stringify(buildArgs(spec.queryIndex, index)) === argsText
      ));
      const searchContext = query.indexes.length > 1 || subject.indexes.length > 1
        ? await hashText(JSON.stringify([query.hash, subject.hash])) : null;
      batch = { query, subject, searchContext, args, specs: [] };
      batches.set(key, batch);
    }
    batch.specs.push(spec);
    bySpec.set(spec, batch);
  }
  return { batches: [...batches.values()], bySpec };
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
