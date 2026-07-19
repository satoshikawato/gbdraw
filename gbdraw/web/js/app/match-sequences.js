import { formatFastaEntry } from './feature-sequence-fasta.js';

const IUPAC_COMPLEMENT = Object.freeze({
  A: 'T', C: 'G', G: 'C', T: 'A', U: 'A',
  R: 'Y', Y: 'R', S: 'S', W: 'W', K: 'M', M: 'K',
  B: 'V', D: 'H', H: 'D', V: 'B', N: 'N',
  '-': '-'
});

const text = (value) => String(value ?? '').trim();
const normalizedSequence = (value) => String(value ?? '').replace(/\s+/g, '').toUpperCase();
const optionalInteger = (value) => {
  if (value === null || value === undefined || value === '') return null;
  const numeric = Number(value);
  return Number.isInteger(numeric) ? numeric : null;
};
const safePart = (value, fallback = 'sequence') => {
  const cleaned = text(value).replace(/[^A-Za-z0-9_.-]+/g, '_').replace(/^_+|_+$/g, '');
  return cleaned || fallback;
};

const aliasesForSource = (source) => {
  const values = [source?.recordId, ...(Array.isArray(source?.aliases) ? source.aliases : [])];
  const aliases = new Set();
  values.forEach((value) => {
    const full = text(value);
    if (!full) return;
    aliases.add(full);
    const first = full.split(/\s+/)[0];
    if (first) aliases.add(first);
  });
  return aliases;
};

export const reverseComplementNucleotide = (sequence) => {
  const input = normalizedSequence(sequence);
  let output = '';
  for (let index = input.length - 1; index >= 0; index -= 1) {
    output += IUPAC_COMPLEMENT[input[index]] || 'N';
  }
  return output;
};

export const validateMatchCoordinates = (startRaw, endRaw, sequenceLength) => {
  const start = Number(startRaw);
  const end = Number(endRaw);
  const length = Number(sequenceLength);
  if (!Number.isInteger(start) || !Number.isInteger(end)) {
    return { valid: false, reason: 'Coordinates must be whole 1-based values.' };
  }
  if (start < 1 || end < 1) {
    return { valid: false, reason: 'Coordinates must be greater than or equal to 1.' };
  }
  if (!Number.isInteger(length) || length < 1) {
    return { valid: false, reason: 'The sequence source is empty.' };
  }
  if (start > length || end > length) {
    return { valid: false, reason: `Coordinates exceed the sequence length (${length.toLocaleString()} bp).` };
  }
  return { valid: true, start, end, orientation: start <= end ? '+' : '-' };
};

export const extractMatchedSpan = (sequence, startRaw, endRaw) => {
  const normalized = normalizedSequence(sequence);
  const validation = validateMatchCoordinates(startRaw, endRaw, normalized.length);
  if (!validation.valid) return { ...validation, sequence: '' };
  const low = Math.min(validation.start, validation.end);
  const high = Math.max(validation.start, validation.end);
  const sliced = normalized.slice(low - 1, high);
  return {
    ...validation,
    sequence: validation.orientation === '-' ? reverseComplementNucleotide(sliced) : sliced,
    sequenceLength: high - low + 1
  };
};

export const createSequenceSourceRegistry = (initialSources = []) => {
  const sources = new Map();

  const register = (source) => {
    const key = text(source?.key);
    const sequence = normalizedSequence(source?.sequence);
    if (!key || !sequence) return null;
    const entry = {
      key,
      recordId: text(source?.recordId) || key,
      aliases: Array.from(aliasesForSource(source)),
      sequence,
      origin: text(source?.origin),
      recordIndex: optionalInteger(source?.recordIndex),
      sourceIndex: optionalInteger(source?.sourceIndex)
    };
    sources.set(key, entry);
    return entry;
  };

  const reset = (nextSources = []) => {
    sources.clear();
    (Array.isArray(nextSources) ? nextSources : []).forEach(register);
  };

  const resolve = (sourceKey, recordId, context = {}) => {
    const explicitKey = text(sourceKey);
    if (explicitKey) {
      const direct = sources.get(explicitKey);
      return direct
        ? { source: direct, reason: '' }
        : { source: null, reason: 'The sequence source for this match is unavailable.' };
    }

    const expectedOrigin = text(context.origin);
    const expectedSourceIndex = optionalInteger(context.sourceIndex);
    const expectedRecordIndex = optionalInteger(context.recordIndex);
    let candidates = Array.from(sources.values()).filter((source) => {
      if (expectedOrigin && source.origin !== expectedOrigin) return false;
      if (expectedSourceIndex !== null && source.sourceIndex !== expectedSourceIndex) return false;
      if (expectedRecordIndex !== null && source.recordIndex !== expectedRecordIndex) return false;
      return true;
    });

    if (expectedRecordIndex !== null && candidates.length === 1) {
      return { source: candidates[0], reason: '' };
    }
    const wanted = text(recordId);
    if (!wanted) return { source: null, reason: 'The match does not identify a sequence record.' };
    const exact = candidates.filter((source) => source.recordId === wanted);
    if (exact.length === 1) return { source: exact[0], reason: '' };
    if (exact.length > 1) {
      return { source: null, reason: `Record ID “${wanted}” is ambiguous in the available sequence sources.` };
    }
    const alias = candidates.filter((source) => aliasesForSource(source).has(wanted));
    if (alias.length === 1) return { source: alias[0], reason: '' };
    if (alias.length > 1) {
      return { source: null, reason: `Record alias “${wanted}” is ambiguous in the available sequence sources.` };
    }
    return { source: null, reason: `No sequence source matched record “${wanted}”.` };
  };

  (Array.isArray(initialSources) ? initialSources : []).forEach(register);
  return { sources, register, reset, resolve, values: () => Array.from(sources.values()) };
};

const normalizeResolution = (resolved) => {
  if (!resolved) return { source: null, reason: '' };
  if (resolved.source !== undefined) return resolved;
  return { source: resolved, reason: '' };
};

export const buildMatchSequenceEntry = (
  span,
  {
    matchId = 'match',
    resolveSequenceSource,
    context = {},
    unavailableReason = ''
  } = {}
) => {
  const normalizedSpan = {
    role: text(span?.role) === 'subject' ? 'subject' : 'query',
    sourceKey: text(span?.sourceKey),
    recordId: text(span?.recordId),
    start: span?.start,
    end: span?.end,
    molecule: 'nucleotide',
    displayRole: text(span?.displayRole) || (text(span?.role) === 'subject' ? 'Subject' : 'Query')
  };
  const unavailable = (reason) => ({
    span: normalizedSpan,
    orientation: Number(normalizedSpan.start) <= Number(normalizedSpan.end) ? '+' : '-',
    sequenceLength: 0,
    fasta: '',
    filename: '',
    available: false,
    unavailableReason: reason || unavailableReason || 'Sequence is unavailable.'
  });
  if (unavailableReason) return unavailable(unavailableReason);
  if (typeof resolveSequenceSource !== 'function') return unavailable('Sequence source resolver is unavailable.');
  const resolution = normalizeResolution(resolveSequenceSource(
    normalizedSpan.sourceKey,
    normalizedSpan.recordId,
    context
  ));
  if (!resolution.source) return unavailable(resolution.reason);
  const extracted = extractMatchedSpan(
    resolution.source.sequence,
    normalizedSpan.start,
    normalizedSpan.end
  );
  if (!extracted.valid) return unavailable(extracted.reason);
  const recordId = normalizedSpan.recordId || resolution.source.recordId;
  const header = [
    `${safePart(matchId, 'match')}_${normalizedSpan.role}`,
    `record=${recordId}`,
    `coords=${normalizedSpan.start}..${normalizedSpan.end}`,
    `strand=${extracted.orientation}`
  ].join('|');
  const fasta = `${formatFastaEntry({ id: header, sequence: extracted.sequence })}\n`;
  return {
    span: { ...normalizedSpan, recordId },
    orientation: extracted.orientation,
    sequenceLength: extracted.sequenceLength,
    fasta,
    filename: `${safePart(matchId)}_${normalizedSpan.role}_${safePart(recordId)}_${normalizedSpan.start}-${normalizedSpan.end}.fna`,
    available: true,
    unavailableReason: ''
  };
};

export const buildMatchSequenceBundle = (
  spans,
  {
    matchId = 'match',
    resolveSequenceSource,
    contextForSpan = () => ({}),
    unavailableReasonForSpan = () => ''
  } = {}
) => {
  const entries = (Array.isArray(spans) ? spans : []).map((span) => buildMatchSequenceEntry(span, {
    matchId,
    resolveSequenceSource,
    context: contextForSpan(span),
    unavailableReason: unavailableReasonForSpan(span)
  }));
  const availableEntries = entries.filter((entry) => entry.available);
  return {
    entries,
    combinedFasta: availableEntries.length === entries.length && entries.length > 1
      ? availableEntries.map((entry) => entry.fasta).join('')
      : '',
    combinedFilename: `${safePart(matchId)}_both.fna`
  };
};

export const parseSequenceRecords = (value) => {
  const input = String(value ?? '');
  const records = [];
  let current = null;
  input.split(/\r?\n/).forEach((line) => {
    if (line.startsWith('>')) {
      if (current) records.push({ ...current, sequence: current.parts.join('').toUpperCase() });
      const header = line.slice(1).trim();
      current = { recordId: header.split(/\s+/)[0] || `record_${records.length + 1}`, header, parts: [] };
    } else if (current && line.trim()) {
      current.parts.push(line.replace(/\s+/g, ''));
    }
  });
  if (current) records.push({ ...current, sequence: current.parts.join('').toUpperCase() });
  return records.filter((record) => record.sequence);
};
