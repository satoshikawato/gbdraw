// @ts-check
import {
  isInternalProteinDisplayId,
  locationPartBound,
  normalizeMemberStrand,
  normalizeStringArray,
  resolveDisplayProteinId
} from './feature-utils.js';
/** @import { FeatureLike } from './feature-utils.js' */

const FASTA_LINE_WIDTH = 60;

/**
 * The feature-like record a FASTA header reads: a catalog row, a popup payload,
 * or an orthogroup or pairwise member payload. Every field is read by name and
 * may be missing.
 * @typedef {{ [key: string]: unknown }} FastaFeature
 */

/**
 * A feature's sequences as Python read them (`FeatureSequencesOf` of `diagram-worker-protocol.js`); the catalog
 * stores none (D5).
 * @typedef {object} FastaSequences
 * @property {string} nucleotideSequence
 * @property {string} aminoAcidSequence
 */

/** @param {unknown} value */
const normalizeHeaderText = (value) =>
  String(value || '').replace(/\s+/g, ' ').trim();

/** @param {unknown} value */
const normalizeFastaId = (value) => {
  const text = normalizeHeaderText(value).replace(/^>+/, '').replace(/\s+/g, '_');
  return text || '';
};

/** @param {unknown} key */
const camelizeKey = (key) =>
  String(key || '').replace(/_([a-z])/g, (_, char) => char.toUpperCase());

/**
 * @param {FastaFeature | null | undefined} feature
 * @param {string} key
 * @returns {string[]}
 */
const getQualifierValues = (feature, key) => {
  const normalizedKey = String(key || '').trim().toLowerCase();
  // The checks leave an object that is not an array; its values are read as opaque.
  const qualifiers = /** @type {Record<string, unknown>} */ (
    feature?.qualifiers && typeof feature.qualifiers === 'object' && !Array.isArray(feature.qualifiers)
      ? feature.qualifiers
      : {}
  );
  if (!normalizedKey) return [];
  if (Object.prototype.hasOwnProperty.call(qualifiers, normalizedKey)) {
    return normalizeStringArray(qualifiers[normalizedKey]);
  }
  const matchingKey = Object.keys(qualifiers).find((candidate) => candidate.toLowerCase() === normalizedKey);
  return matchingKey ? normalizeStringArray(qualifiers[matchingKey]) : [];
};

/**
 * @param {FastaFeature | null | undefined} feature
 * @param {string} key
 * @returns {string[]}
 */
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

/**
 * @param {FastaFeature | null | undefined} feature
 * @param {string[]} keys
 */
const firstFeatureText = (feature, keys) => {
  for (const key of keys) {
    const value = getFeatureValues(feature, key)
      .map((entry) => normalizeHeaderText(entry))
      .find(Boolean);
    if (value) return value;
  }
  return '';
};

/**
 * @param {FastaFeature | null | undefined} feature
 * @param {string[]} keys
 */
const firstFeatureDisplayText = (feature, keys) => {
  for (const key of keys) {
    const value = getFeatureValues(feature, key)
      .map((entry) => normalizeHeaderText(entry))
      .find((entry) => entry && !isInternalProteinDisplayId(entry));
    if (value) return value;
  }
  return '';
};

/** @param {unknown} value */
const sequenceText = (value) => String(value || '').replace(/\s+/g, '').trim();

/**
 * @param {unknown} sequence
 * @param {unknown} [width]
 */
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

/**
 * @param {{ id: string, description?: string, sequence: string }} entry
 */
export const formatFastaEntry = ({ id, description, sequence }) => {
  const normalizedId = normalizeFastaId(id) || 'sequence';
  const wrappedSequence = wrapFastaSequence(sequence);
  if (!wrappedSequence) return '';
  const headerDescription = normalizeHeaderText(description);
  const header = headerDescription ? `>${normalizedId} ${headerDescription}` : `>${normalizedId}`;
  return `${header}\n${wrappedSequence}`;
};

/** @param {FastaFeature | null | undefined} feature */
const getFeatureOrganism = (feature) =>
  firstFeatureText(feature, ['organism', 'source_organism', 'record_organism']);

/** @param {FastaFeature | null | undefined} feature */
const getFeatureDescription = (feature) => {
  const label = firstFeatureDisplayText(feature, [
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

/** @param {unknown} strand */
const strandIsMinus = (strand) => String(strand || '').trim() === '-';

/**
 * @param {unknown} part A location part, or the feature read as its one part.
 * @param {unknown} strandRaw
 */
const formatLocationToken = (part, strandRaw) => {
  const start = locationPartBound(part, 'start');
  const end = locationPartBound(part, 'end');
  if (!Number.isFinite(start) || !Number.isFinite(end)) return '';
  const startOneBased = Math.max(1, Math.trunc(start) + 1);
  const endOneBased = Math.max(startOneBased, Math.trunc(end));
  return strandIsMinus(strandRaw)
    ? `c${endOneBased}-${startOneBased}`
    : `${startOneBased}-${endOneBased}`;
};

/** @param {FastaFeature | null | undefined} feature */
const getFeatureLocationToken = (feature) => {
  const parts = Array.isArray(feature?.location_parts) ? feature.location_parts : [];
  const partTokens = parts
    .map((part) => formatLocationToken(part, part?.strand || feature?.strand))
    .filter(Boolean);
  if (partTokens.length === 1) return partTokens[0];
  if (partTokens.length > 1) return `join(${partTokens.join(',')})`;

  const directToken = formatLocationToken(feature, feature?.strand);
  if (directToken) return directToken;

  const location = normalizeFastaId(feature?.location || '');
  return location || 'feature';
};

/** @param {FastaFeature | null | undefined} feature */
const getNucleotideFastaId = (feature) => {
  const recordId = firstFeatureText(feature, ['record_id', 'recordId']) || 'record';
  return `${recordId}:${getFeatureLocationToken(feature)}`;
};

/** @param {FastaFeature | null | undefined} feature */
const getProteinFastaId = (feature) =>
  firstFeatureDisplayText(feature, [
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

/**
 * The header fields of a Similarity group member's FASTA entries, for the orthogroup panel and the match popup alike:
 * the member's catalog row first, the member's own fields for what the row lacks, and the protein ID the popups show
 * (`resolveDisplayProteinId`); the header skips an internal source protein ID. The bases are read from Python when a
 * FASTA action is pressed.
 * @param {FeatureLike | null | undefined} member
 * @param {FeatureLike | null | undefined} row The member's catalog row, when it resolved to one.
 * @param {string} locationText The member's location as its popup writes it: the protein ID's last fallback.
 * @returns {FastaFeature}
 */
export const memberFastaFeature = (member, row, locationText) => {
  const sourceFeature = row && typeof row === 'object' ? row : {};
  const displayProteinId = resolveDisplayProteinId(sourceFeature, member, locationText);
  return {
    ...sourceFeature,
    record_id: sourceFeature.record_id || sourceFeature.recordId || member?.recordId || member?.record_id,
    start: sourceFeature.start ?? member?.start,
    end: sourceFeature.end ?? member?.end,
    strand: sourceFeature.strand || normalizeMemberStrand(member?.strand),
    source_protein_id: sourceFeature.source_protein_id || sourceFeature.sourceProteinId || member?.sourceProteinId || member?.source_protein_id,
    protein_id: displayProteinId,
    product: sourceFeature.product || member?.product,
    note: sourceFeature.note || member?.note,
    gene: sourceFeature.gene || member?.gene,
    organism: sourceFeature.organism || member?.organism
  };
};

/**
 * The FASTA texts of a feature: the header from the feature's fields, the entries from its read sequences.
 * @param {FastaFeature | null | undefined} feature
 * @param {FastaSequences} sequences
 */
export const buildFeatureSequenceFastas = (feature, sequences) => {
  const nucleotideSequence = sequenceText(sequences.nucleotideSequence);
  const aminoAcidSequence = sequenceText(sequences.aminoAcidSequence);
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
