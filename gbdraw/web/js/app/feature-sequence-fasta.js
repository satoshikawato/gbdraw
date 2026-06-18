const FASTA_LINE_WIDTH = 60;

const normalizeStringArray = (value) => {
  if (Array.isArray(value)) {
    return value
      .filter((item) => item !== null && item !== undefined)
      .map((item) => String(item));
  }
  if (value === null || value === undefined || value === '') return [];
  return [String(value)];
};

const normalizeHeaderText = (value) =>
  String(value || '').replace(/\s+/g, ' ').trim();

const normalizeFastaId = (value) => {
  const text = normalizeHeaderText(value).replace(/^>+/, '').replace(/\s+/g, '_');
  return text || '';
};

const camelizeKey = (key) =>
  String(key || '').replace(/_([a-z])/g, (_, char) => char.toUpperCase());

const getQualifierValues = (feature, key) => {
  const normalizedKey = String(key || '').trim().toLowerCase();
  const qualifiers = feature?.qualifiers && typeof feature.qualifiers === 'object' && !Array.isArray(feature.qualifiers)
    ? feature.qualifiers
    : {};
  if (!normalizedKey) return [];
  if (Object.prototype.hasOwnProperty.call(qualifiers, normalizedKey)) {
    return normalizeStringArray(qualifiers[normalizedKey]);
  }
  const matchingKey = Object.keys(qualifiers).find((candidate) => candidate.toLowerCase() === normalizedKey);
  return matchingKey ? normalizeStringArray(qualifiers[matchingKey]) : [];
};

const getFeatureValues = (feature, key) => {
  if (!feature || typeof feature !== 'object') return [];
  const candidateKeys = [key, camelizeKey(key)];
  for (const candidateKey of candidateKeys) {
    if (Object.prototype.hasOwnProperty.call(feature, candidateKey)) {
      const values = normalizeStringArray(feature[candidateKey]);
      if (values.length > 0) return values;
    }
  }
  return getQualifierValues(feature, key);
};

const firstFeatureText = (feature, keys) => {
  for (const key of keys) {
    const value = getFeatureValues(feature, key)
      .map((entry) => normalizeHeaderText(entry))
      .find(Boolean);
    if (value) return value;
  }
  return '';
};

const sequenceText = (value) => String(value || '').replace(/\s+/g, '').trim();

export const wrapFastaSequence = (sequence, width = FASTA_LINE_WIDTH) => {
  const text = sequenceText(sequence);
  if (!text) return '';
  const lines = [];
  const lineWidth = Math.max(1, Number(width) || FASTA_LINE_WIDTH);
  for (let i = 0; i < text.length; i += lineWidth) {
    lines.push(text.slice(i, i + lineWidth));
  }
  return lines.join('\n');
};

export const formatFastaEntry = ({ id, description, sequence }) => {
  const normalizedId = normalizeFastaId(id) || 'sequence';
  const wrappedSequence = wrapFastaSequence(sequence);
  if (!wrappedSequence) return '';
  const headerDescription = normalizeHeaderText(description);
  const header = headerDescription ? `>${normalizedId} ${headerDescription}` : `>${normalizedId}`;
  return `${header}\n${wrappedSequence}`;
};

const getFeatureOrganism = (feature) =>
  firstFeatureText(feature, ['organism', 'source_organism', 'record_organism']);

const getFeatureDescription = (feature) => {
  const label = firstFeatureText(feature, [
    'product',
    'protein',
    'gene',
    'locus_tag',
    'name',
    'standard_name',
    'note',
    'type'
  ]);
  const organism = getFeatureOrganism(feature);
  if (!organism) return label;
  if (!label) return `[${organism}]`;
  return /\[[^\]]+\]\s*$/.test(label) ? label : `${label} [${organism}]`;
};

const strandIsMinus = (strand) => String(strand || '').trim() === '-';

const formatLocationToken = (startRaw, endRaw, strandRaw) => {
  const start = Number(startRaw);
  const end = Number(endRaw);
  if (!Number.isFinite(start) || !Number.isFinite(end)) return '';
  const startOneBased = Math.max(1, Math.trunc(start) + 1);
  const endOneBased = Math.max(startOneBased, Math.trunc(end));
  return strandIsMinus(strandRaw)
    ? `c${endOneBased}-${startOneBased}`
    : `${startOneBased}-${endOneBased}`;
};

const getFeatureLocationToken = (feature) => {
  const parts = Array.isArray(feature?.location_parts) ? feature.location_parts : [];
  const partTokens = parts
    .map((part) => formatLocationToken(part?.start, part?.end, part?.strand || feature?.strand))
    .filter(Boolean);
  if (partTokens.length === 1) return partTokens[0];
  if (partTokens.length > 1) return `join(${partTokens.join(',')})`;

  const directToken = formatLocationToken(feature?.start, feature?.end, feature?.strand);
  if (directToken) return directToken;

  const location = normalizeFastaId(feature?.location || '');
  return location || 'feature';
};

const getNucleotideFastaId = (feature) => {
  const recordId = firstFeatureText(feature, ['record_id', 'recordId']) || 'record';
  return `${recordId}:${getFeatureLocationToken(feature)}`;
};

const getProteinFastaId = (feature) =>
  firstFeatureText(feature, [
    'source_protein_id',
    'protein_id',
    'locus_tag',
    'gene',
    'gene_id',
    'old_locus_tag',
    'gff_id',
    'transcript_id',
    'name'
  ]) || getNucleotideFastaId(feature);

export const buildFeatureSequenceFastas = (feature, sequences = {}) => {
  const nucleotideSequence = sequenceText(
    sequences.nucleotideSequence ?? sequences.nucleotide_sequence ?? feature?.nucleotide_sequence ?? feature?.nucleotideSequence
  );
  const aminoAcidSequence = sequenceText(
    sequences.aminoAcidSequence ?? sequences.amino_acid_sequence ?? feature?.amino_acid_sequence ?? feature?.aminoAcidSequence
  );
  const description = getFeatureDescription(feature);
  return {
    nucleotideFasta: formatFastaEntry({
      id: getNucleotideFastaId(feature),
      description,
      sequence: nucleotideSequence
    }),
    aminoAcidFasta: formatFastaEntry({
      id: getProteinFastaId(feature),
      description,
      sequence: aminoAcidSequence
    })
  };
};
