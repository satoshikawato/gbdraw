import { resolveDisplayProteinId } from './feature-utils.js';
import { buildFeatureSequenceFastas } from './feature-sequence-fasta.js';
import {
  groupMetadataScopeLabel,
  normalizeGroupMetadataScope
} from './losat-normalization.js';

export const PAIRWISE_MATCH_SELECTOR = [
  'path[data-gbdraw-pairwise-match-id]',
  'path[data-match-kind]',
  'path[data-pairwise-match-style]'
].join(', ');

const MATCH_KIND_TITLES = {
  pairwise: 'Pairwise match',
  orthogroup: 'Orthogroup match',
  collinear: 'Collinearity block'
};

const MATCH_KIND_ALIASES = new Set(['pairwise', 'orthogroup', 'collinear']);
const normalizeText = (value) => String(value === null || value === undefined ? '' : value).trim();
const normalizeSequence = (value) => String(value ?? '').replace(/\s+/g, '').toUpperCase();

const attr = (element, name) => normalizeText(element?.getAttribute?.(name));

const firstText = (...values) => {
  for (const value of values) {
    const text = normalizeText(value);
    if (text) return text;
  }
  return '';
};

const splitMetadataValues = (value) => normalizeText(value)
  .split(';')
  .map((entry) => normalizeText(entry))
  .filter(Boolean);

const uniqueMetadataValues = (value) => {
  const seen = new Set();
  const values = [];
  splitMetadataValues(value).forEach((entry) => {
    if (seen.has(entry)) return;
    seen.add(entry);
    values.push(entry);
  });
  return values;
};

const generatedProteinIdPattern = /^(?:gbd_r\d+_cds\d+|p_.+_\d+_\d+_-?\d+_[0-9a-f]{12}(?:_\d+)?)$/i;
const generatedUnitIdPattern = /^gbd_r\d+_unit\d+$/i;

const isInternalDisplayId = (value) => {
  const text = normalizeText(value);
  return Boolean(text && (
    generatedProteinIdPattern.test(text) ||
    generatedUnitIdPattern.test(text)
  ));
};

const addUniqueDisplayText = (values, value) => {
  const text = normalizeText(value);
  if (!text || isInternalDisplayId(text) || values.includes(text)) return;
  values.push(text);
};

const firstSequenceText = (...values) => {
  for (const value of values) {
    const sequence = normalizeSequence(value);
    if (sequence) return sequence;
  }
  return '';
};

const makeSafeFilename = (value, fallback = 'orthogroup') => {
  const cleaned = normalizeText(value).replace(/[^\w.-]+/g, '_').replace(/^_+|_+$/g, '');
  return cleaned || fallback;
};

const sequenceKindLabel = (sequenceKind) => (sequenceKind === 'aa' ? 'aa' : 'nt');
const sequenceExtension = (sequenceKind) => (sequenceKindLabel(sequenceKind) === 'aa' ? 'faa' : 'fna');

const normalizeMatchKind = (value, element = null) => {
  const normalized = normalizeText(value).toLowerCase();
  if (MATCH_KIND_ALIASES.has(normalized)) return normalized;
  if (attr(element, 'data-collinearity-block-id')) return 'collinear';
  if (attr(element, 'data-orthogroup-id')) return 'orthogroup';
  return 'pairwise';
};

const addRow = (rows, label, value) => {
  const text = normalizeText(value);
  if (!text) return;
  rows.push({ label, value: text });
};

const intervalText = (start, end) => {
  const startText = normalizeText(start);
  const endText = normalizeText(end);
  if (!startText && !endText) return '';
  if (startText && endText) return `${startText}..${endText}`;
  return startText || endText;
};

const featureLocationText = (feature) => {
  if (!feature) return '';
  const direct = normalizeText(feature.location);
  if (direct) return direct;
  const start = Number(feature.start);
  const end = Number(feature.end);
  const startText = Number.isFinite(start) ? String(start + 1) : normalizeText(feature.start);
  const endText = Number.isFinite(end) ? String(end) : normalizeText(feature.end);
  const strand = normalizeText(feature.strand);
  const range = startText && endText ? `${startText}..${endText}` : startText || endText;
  return range && strand ? `${range} (${strand})` : range;
};

const qualifierFirstValue = (feature, key) => {
  const normalizedKey = normalizeText(key).toLowerCase();
  if (!feature || !normalizedKey) return '';
  const direct = feature[normalizedKey];
  if (direct !== null && direct !== undefined && direct !== '') return firstArrayValue(direct);
  const qualifiers = feature.qualifiers && typeof feature.qualifiers === 'object' && !Array.isArray(feature.qualifiers)
    ? feature.qualifiers
    : {};
  if (qualifiers[normalizedKey] !== null && qualifiers[normalizedKey] !== undefined) {
    return firstArrayValue(qualifiers[normalizedKey]);
  }
  const matchingKey = Object.keys(qualifiers).find((candidate) => candidate.toLowerCase() === normalizedKey);
  return matchingKey ? firstArrayValue(qualifiers[matchingKey]) : '';
};

const firstArrayValue = (value) => {
  if (Array.isArray(value)) {
    for (const entry of value) {
      const text = normalizeText(entry);
      if (text) return text;
    }
    return '';
  }
  return normalizeText(value);
};

const getOrthogroupById = (orthogroups, orthogroupId) => {
  const id = normalizeText(orthogroupId);
  if (!id) return null;
  return (Array.isArray(orthogroups) ? orthogroups : [])
    .find((entry) => normalizeText(entry?.id || entry?.orthogroupId || entry?.orthogroup_id) === id) || null;
};

const memberFeatureSvgId = (member) => firstText(member?.featureSvgId, member?.feature_svg_id);
const memberRecordIndex = (member) => {
  const recordIndex = Number(member?.recordIndex ?? member?.record_index);
  return Number.isInteger(recordIndex) ? recordIndex : null;
};

const featureRecordIndex = (feature) => {
  const recordIndex = Number(feature?.fileIdx ?? feature?.recordIndex ?? feature?.record_index);
  return Number.isInteger(recordIndex) ? recordIndex : null;
};

const getFeatureLookupValues = (featureLookup) => {
  if (!featureLookup) return [];
  if (typeof featureLookup.values === 'function') return Array.from(featureLookup.values());
  if (Array.isArray(featureLookup)) return featureLookup;
  if (typeof featureLookup === 'object') return Object.values(featureLookup);
  return [];
};

const getFeatureForMember = (member, featureLookup) => {
  const svgId = memberFeatureSvgId(member);
  if (!svgId || !featureLookup) return null;
  const recordIndex = memberRecordIndex(member);
  const renderedSvgId = recordIndex === null ? '' : `${svgId}_record_${recordIndex + 1}`;
  const direct = featureLookup.get?.(renderedSvgId) || featureLookup.get?.(svgId) || null;
  if (direct && (recordIndex === null || featureRecordIndex(direct) === recordIndex)) return direct;
  return getFeatureLookupValues(featureLookup).find((feature) => {
    const featureSvgIds = [
      feature?.svg_id,
      feature?.svgId,
      feature?.stable_svg_id,
      feature?.stableSvgId
    ].map(normalizeText).filter(Boolean);
    if (!featureSvgIds.includes(svgId)) return false;
    return recordIndex === null || featureRecordIndex(feature) === recordIndex;
  }) || direct || null;
};

const groupHasFeatureSvgId = (group, featureSvgId) => {
  const ids = splitMetadataValues(featureSvgId);
  if (ids.length === 0) return false;
  return (Array.isArray(group?.members) ? group.members : [])
    .some((member) => ids.includes(memberFeatureSvgId(member)));
};

const getOrthogroupForMatch = (
  orthogroups,
  {
    orthogroupId,
    queryFeatureSvgId,
    subjectFeatureSvgId
  } = {}
) => {
  const direct = getOrthogroupById(orthogroups, orthogroupId);
  if (direct) return direct;
  const groups = Array.isArray(orthogroups) ? orthogroups : [];
  return groups.find((group) => (
    groupHasFeatureSvgId(group, queryFeatureSvgId) ||
    groupHasFeatureSvgId(group, subjectFeatureSvgId)
  )) || null;
};

const featureOrthogroupId = (feature) => firstText(
  feature?.orthogroupId,
  feature?.orthogroup_id,
  feature?.orthogroupMember?.orthogroupId,
  feature?.orthogroupMember?.orthogroup_id,
  feature?.orthogroup_member?.orthogroupId,
  feature?.orthogroup_member?.orthogroup_id
);

const featureProduct = (feature) => firstText(
  feature?.product,
  qualifierFirstValue(feature, 'product'),
  feature?.orthogroupMember?.product,
  feature?.orthogroup_member?.product,
  feature?.note,
  qualifierFirstValue(feature, 'note')
);

const featureToOrthogroupMember = (feature) => {
  if (!feature) return null;
  const member = feature?.orthogroupMember || feature?.orthogroup_member || {};
  return {
    recordId: firstText(member?.recordId, member?.record_id, feature?.record_id),
    start: member?.start ?? feature?.start,
    end: member?.end ?? feature?.end,
    strand: firstText(member?.strand, feature?.strand),
    proteinId: firstText(member?.proteinId, member?.protein_id, feature?.proteinId, feature?.protein_id),
    sourceProteinId: firstText(
      member?.sourceProteinId,
      member?.source_protein_id,
      feature?.sourceProteinId,
      feature?.source_protein_id,
      qualifierFirstValue(feature, 'protein_id')
    ),
    product: firstText(member?.product, featureProduct(feature)),
    note: firstText(member?.note, feature?.note, qualifierFirstValue(feature, 'note')),
    featureSvgId: firstText(member?.featureSvgId, member?.feature_svg_id, feature?.svg_id)
  };
};

const buildFallbackOrthogroup = ({
  orthogroupId,
  queryFeature,
  subjectFeature,
  featureLookup
}) => {
  const id = normalizeText(orthogroupId);
  if (!id) return null;
  const features = getFeatureLookupValues(featureLookup);
  const matchingFeatures = features.filter((feature) => featureOrthogroupId(feature) === id);
  const fallbackFeatures = matchingFeatures.length
    ? matchingFeatures
    : [queryFeature, subjectFeature].filter(Boolean);
  if (!fallbackFeatures.length) return null;
  const members = fallbackFeatures
    .map(featureToOrthogroupMember)
    .filter(Boolean);
  const firstFeature = fallbackFeatures[0] || {};
  const recordCoverageFallback = new Set(
    fallbackFeatures
      .map((feature) => firstText(feature?.record_id, feature?.recordId, feature?.record_idx, feature?.recordIndex))
      .filter(Boolean)
  ).size;
  return {
    id,
    name: featureProduct(firstFeature),
    member_count: firstText(firstFeature?.orthogroupMemberCount, firstFeature?.orthogroup_member_count, members.length),
    record_coverage_count: firstText(firstFeature?.orthogroupRecordCoverage, firstFeature?.orthogroup_record_coverage, recordCoverageFallback),
    members
  };
};

const overrideValue = (overrides, key) => {
  const normalizedKey = normalizeText(key);
  if (!normalizedKey || !overrides) return '';
  if (overrides instanceof Map) return normalizeText(overrides.get(normalizedKey));
  return normalizeText(overrides[normalizedKey]);
};

const section = (title, rows, extras = {}) => ({
  title,
  rows: rows.filter((row) => normalizeText(row.value)),
  ...extras
});

const orthogroupTitle = (orthogroupId, displayName) => {
  const id = normalizeText(orthogroupId);
  const name = normalizeText(displayName);
  if (id && name) return `${id}:${name}`;
  return id || name || MATCH_KIND_TITLES.orthogroup;
};

const memberLocationText = (member) => {
  if (!member || typeof member !== 'object') return '';
  const start = Number(member.start);
  const end = Number(member.end);
  const startText = Number.isFinite(start) ? String(start + 1) : normalizeText(member.start);
  const endText = Number.isFinite(end) ? String(end) : normalizeText(member.end);
  const strand = normalizeText(member.strand);
  const range = startText && endText ? `${startText}..${endText}` : startText || endText;
  return range && strand ? `${range} (${strand})` : range;
};

const memberSequence = (member, feature, sequenceKind) => {
  if (sequenceKindLabel(sequenceKind) === 'aa') {
    return firstSequenceText(
      member?.aminoAcidSequence,
      member?.amino_acid_sequence,
      member?.proteinSequence,
      member?.sequence,
      feature?.aminoAcidSequence,
      feature?.amino_acid_sequence
    );
  }
  return firstSequenceText(
    member?.nucleotideSequence,
    member?.nucleotide_sequence,
    feature?.nucleotideSequence,
    feature?.nucleotide_sequence
  );
};

const normalizeMemberStrand = (strand) => {
  if (strand === -1 || String(strand).trim() === '-1') return '-';
  if (strand === 1 || String(strand).trim() === '1') return '+';
  return normalizeText(strand);
};

const buildMemberFeaturePayload = (member, feature, nucleotideSequence, aminoAcidSequence) => {
  const sourceFeature = feature && typeof feature === 'object' ? feature : {};
  return {
    ...sourceFeature,
    record_id: sourceFeature.record_id || sourceFeature.recordId || member?.recordId || member?.record_id,
    start: sourceFeature.start ?? member?.start,
    end: sourceFeature.end ?? member?.end,
    strand: sourceFeature.strand || normalizeMemberStrand(member?.strand),
    source_protein_id: sourceFeature.source_protein_id || sourceFeature.sourceProteinId || member?.sourceProteinId || member?.source_protein_id,
    protein_id: sourceFeature.protein_id || sourceFeature.proteinId || member?.proteinId || member?.protein_id,
    product: sourceFeature.product || member?.product,
    note: sourceFeature.note || member?.note,
    gene: sourceFeature.gene || member?.gene,
    organism: sourceFeature.organism || member?.organism,
    nucleotide_sequence: nucleotideSequence,
    amino_acid_sequence: aminoAcidSequence
  };
};

const memberFastaText = (member, feature, nucleotideSequence, aminoAcidSequence, sequenceKind) => {
  const fastaFeature = buildMemberFeaturePayload(member, feature, nucleotideSequence, aminoAcidSequence);
  const fastas = buildFeatureSequenceFastas(fastaFeature, {
    nucleotideSequence,
    aminoAcidSequence
  });
  const text = sequenceKindLabel(sequenceKind) === 'aa' ? fastas.aminoAcidFasta : fastas.nucleotideFasta;
  return text ? `${text}\n` : '';
};

const orthogroupSequenceFilename = (orthogroupId, displayName, sequenceKind) => {
  const id = normalizeText(orthogroupId) || 'orthogroup';
  const name = makeSafeFilename(displayName || id, id);
  return `${makeSafeFilename(`${id}_${name}_${sequenceKindLabel(sequenceKind)}`)}.${sequenceExtension(sequenceKind)}`;
};

const orthogroupMemberSequenceFilename = (member, orthogroupId, sequenceKind) => {
  const id = normalizeText(orthogroupId) || 'orthogroup';
  const memberId = firstText(
    member?.sourceProteinId,
    member?.source_protein_id,
    member?.proteinId,
    member?.protein_id,
    memberFeatureSvgId(member),
    'member'
  );
  return `${makeSafeFilename(`${id}_${memberId}_${sequenceKindLabel(sequenceKind)}`)}.${sequenceExtension(sequenceKind)}`;
};

const buildOrthogroupMemberRows = (group, featureLookup, orthogroupId) => {
  const members = Array.isArray(group?.members) ? group.members : [];
  return members
    .map((member) => {
      const feature = getFeatureForMember(member, featureLookup);
      const nucleotideSequence = memberSequence(member, feature, 'nt');
      const aminoAcidSequence = memberSequence(member, feature, 'aa');
      return {
        record: firstText(member?.recordId, member?.record_id),
        coordinates: memberLocationText(member),
        proteinId: resolveDisplayProteinId(null, member),
        role: firstText(member?.role, member?.memberRole, member?.member_role),
        confidence: firstText(member?.confidence, member?.memberConfidence, member?.member_confidence),
        assignmentReason: firstText(member?.assignmentReason, member?.assignment_reason),
        productOrNote: firstText(member?.product, member?.note),
        ntFasta: memberFastaText(member, feature, nucleotideSequence, aminoAcidSequence, 'nt'),
        aaFasta: memberFastaText(member, feature, nucleotideSequence, aminoAcidSequence, 'aa'),
        ntFilename: orthogroupMemberSequenceFilename(member, orthogroupId, 'nt'),
        aaFilename: orthogroupMemberSequenceFilename(member, orthogroupId, 'aa')
      };
    })
    .filter((row) => row.record || row.coordinates || row.proteinId || row.role || row.confidence || row.assignmentReason || row.productOrNote || row.ntFasta || row.aaFasta);
};

const getGroupMemberForFeatureSvgId = (group, featureSvgId) => {
  const id = normalizeText(featureSvgId);
  if (!id) return null;
  const members = Array.isArray(group?.members) ? group.members : [];
  return members.find((member) => memberFeatureSvgId(member) === id) || null;
};

const resolveFeatureSectionProteinIds = ({
  feature,
  featureSvgIds,
  featureLookup,
  group,
  fallbackProteinIds,
  locusId,
  displayName
}) => {
  const values = [];
  const addFeatureProteinId = (candidateFeature, member = null) => {
    const text = resolveDisplayProteinId(candidateFeature, member, '');
    addUniqueDisplayText(values, text);
  };

  const ids = splitMetadataValues(featureSvgIds);
  ids.forEach((featureSvgId) => {
    const candidateFeature = featureLookup?.get?.(featureSvgId) || null;
    const member = getGroupMemberForFeatureSvgId(group, featureSvgId);
    addFeatureProteinId(candidateFeature, member);
  });
  if (feature) {
    addFeatureProteinId(feature, ids.length === 1 ? getGroupMemberForFeatureSvgId(group, ids[0]) : null);
  }
  if (values.length === 0) {
    splitMetadataValues(locusId).forEach((value) => addUniqueDisplayText(values, value));
  }
  if (values.length === 0) {
    splitMetadataValues(displayName).forEach((value) => addUniqueDisplayText(values, value));
  }
  if (values.length === 0) {
    splitMetadataValues(fallbackProteinIds).forEach((value) => addUniqueDisplayText(values, value));
  }
  return values.join('; ');
};

const firstNonInternalDisplayText = (...values) => {
  for (const value of values) {
    const text = normalizeText(value);
    if (text && !isInternalDisplayId(text)) return text;
  }
  return '';
};

const featureRowDisplayName = (feature, fallbackDisplayName) => firstText(
  fallbackDisplayName,
  feature?.displayLabel,
  feature?.display_label,
  feature?.label,
  feature?.gene,
  qualifierFirstValue(feature, 'gene'),
  feature?.locus_tag,
  feature?.locusTag,
  qualifierFirstValue(feature, 'locus_tag'),
  featureProduct(feature)
);

const featureRowLocusId = (feature, fallbackLocusId) => firstText(
  feature?.locus_tag,
  feature?.locusTag,
  qualifierFirstValue(feature, 'locus_tag'),
  feature?.gene_id,
  feature?.geneId,
  qualifierFirstValue(feature, 'gene_id'),
  fallbackLocusId
);

const featureRowSubLabel = (primaryLabel, ...values) => {
  const seen = new Set(splitMetadataValues(primaryLabel));
  const labels = [];
  values.forEach((value) => {
    splitMetadataValues(value).forEach((entry) => {
      if (seen.has(entry)) return;
      seen.add(entry);
      labels.push(entry);
    });
  });
  return labels.join(' / ');
};

const buildFeatureListRows = ({
  featureSvgIds,
  featureLookup,
  group,
  recordId,
  interval,
  proteinId,
  locusId,
  displayName
}) => {
  const featureIds = uniqueMetadataValues(featureSvgIds);
  const proteinIds = splitMetadataValues(proteinId);
  const locusIds = splitMetadataValues(locusId);
  const displayNames = splitMetadataValues(displayName);
  const count = Math.max(featureIds.length, proteinIds.length, locusIds.length, displayNames.length);
  if (count === 0) return [];

  return Array.from({ length: count }, (_unused, index) => {
    const svgId = featureIds[index] || '';
    const feature = svgId ? featureLookup?.get?.(svgId) || null : null;
    const member = svgId ? getGroupMemberForFeatureSvgId(group, svgId) : null;
    const fallbackProteinId = firstNonInternalDisplayText(
      locusIds[index],
      displayNames[index],
      proteinIds[index]
    );
    const resolvedProteinId = resolveDisplayProteinId(feature, member, '');
    const displayProteinId = firstText(
      isInternalDisplayId(resolvedProteinId) ? '' : resolvedProteinId,
      fallbackProteinId,
      resolvedProteinId,
      proteinIds[index]
    );
    const rowRecord = firstText(feature?.record_id, feature?.recordId, member?.recordId, member?.record_id, recordId);
    const rowLocation = firstText(
      featureLocationText(feature),
      memberLocationText(member),
      count === 1 ? interval : ''
    );
    const rowLocusId = featureRowLocusId(feature, locusIds[index]);
    const rowDisplayName = featureRowDisplayName(feature, displayNames[index]);
    const product = featureProduct(feature);
    const label = firstText(displayProteinId, rowDisplayName, rowLocusId, product, svgId, `Feature ${index + 1}`);
    const subLabel = featureRowSubLabel(label, rowLocusId, rowDisplayName);
    const copyText = [
      rowRecord,
      rowLocation,
      displayProteinId,
      rowLocusId,
      rowDisplayName,
      product
    ].join('\t');
    return {
      key: `${svgId || 'feature'}-${index}`,
      svgId,
      feature,
      canOpen: Boolean(feature?.svg_id),
      label,
      record: rowRecord,
      location: rowLocation,
      proteinId: displayProteinId,
      locusId: rowLocusId,
      displayName: rowDisplayName,
      subLabel,
      product,
      type: firstText(feature?.type),
      copyText
    };
  }).filter((row) => (
    row.svgId ||
    row.record ||
    row.location ||
    row.proteinId ||
    row.locusId ||
    row.displayName ||
    row.product
  ));
};

const orthogroupMemberCopyText = (memberRows) => {
  if (!Array.isArray(memberRows) || memberRows.length === 0) return '';
  return [
    'Record\tCoordinates (+/-)\tProtein ID\tRole\tConfidence\tAssignment reason\tProduct / note',
    ...memberRows.map((row) => [
      row.record,
      row.coordinates,
      row.proteinId,
      row.role,
      row.confidence,
      row.assignmentReason,
      row.productOrNote
    ].join('\t'))
  ].join('\n');
};

const orthogroupMemberSequenceText = (memberRows, sequenceKind) => {
  const key = sequenceKindLabel(sequenceKind) === 'aa' ? 'aaFasta' : 'ntFasta';
  return (Array.isArray(memberRows) ? memberRows : [])
    .map((row) => String(row?.[key] ?? ''))
    .filter((text) => text.trim())
    .join('');
};

const orthogroupMemberSectionExtras = (memberRows, orthogroupId, displayName) => {
  const ntFasta = orthogroupMemberSequenceText(memberRows, 'nt');
  const aaFasta = orthogroupMemberSequenceText(memberRows, 'aa');
  return {
    memberRows,
    memberCopyText: orthogroupMemberCopyText(memberRows),
    memberNtFasta: ntFasta,
    memberAaFasta: aaFasta,
    memberNtCount: memberRows.filter((row) => normalizeText(row.ntFasta)).length,
    memberAaCount: memberRows.filter((row) => normalizeText(row.aaFasta)).length,
    memberNtFilename: orthogroupSequenceFilename(orthogroupId, displayName, 'nt'),
    memberAaFilename: orthogroupSequenceFilename(orthogroupId, displayName, 'aa')
  };
};

const resolveBlockMemberLabels = ({
  group,
  featureSvgIds,
  featureLookup
}) => {
  if (!group) return '';
  const values = [];
  uniqueMetadataValues(featureSvgIds).forEach((featureSvgId) => {
    const feature = featureLookup?.get?.(featureSvgId) || null;
    const member = group ? getGroupMemberForFeatureSvgId(group, featureSvgId) : null;
    if (!member) return;
    addUniqueDisplayText(values, resolveDisplayProteinId(feature, member, ''));
  });
  return values.join('; ');
};

const buildOrthogroupDetailRows = ({
  orthogroupId,
  idLabel = 'Orthogroup ID',
  displayName,
  description,
  scopeLabel,
  memberCount,
  recordCoverage,
  rbhOrthogroups,
  orthologPathCount,
  relatedEdgeCount
}) => {
  const rows = [];
  addRow(rows, idLabel, orthogroupId);
  addRow(rows, 'Display name', displayName);
  addRow(rows, 'Description', description);
  addRow(rows, 'Scope', scopeLabel);
  addRow(rows, 'Members', memberCount);
  addRow(rows, 'Record coverage', recordCoverage);
  addRow(rows, 'RBH seeds', Array.isArray(rbhOrthogroups) ? rbhOrthogroups.join('; ') : rbhOrthogroups);
  addRow(rows, 'Ortholog paths', orthologPathCount);
  addRow(rows, 'Related edges', relatedEdgeCount);
  return rows;
};

const buildBlockOrthogroups = ({
  orthogroupIds,
  orthogroups,
  featureLookup,
  queryFeatureSvgId,
  subjectFeatureSvgId,
  orthogroupNameOverrides,
  orthogroupDescriptionOverrides,
  groupScope
}) => orthogroupIds.map((orthogroupId) => {
  const group = getOrthogroupById(orthogroups, orthogroupId);
  const normalizedGroupScope = normalizeGroupMetadataScope(group?.scope || groupScope);
  const displayName = firstText(
    overrideValue(orthogroupNameOverrides, orthogroupId),
    group?.displayName,
    group?.display_name,
    group?.name
  );
  const description = firstText(
    overrideValue(orthogroupDescriptionOverrides, orthogroupId),
    group?.description
  );
  const memberCount = firstText(
    group?.member_count,
    group?.memberCount,
    Array.isArray(group?.members) ? group.members.length : ''
  );
  const scopeLabel = groupMetadataScopeLabel(normalizedGroupScope);
  const idLabel = normalizedGroupScope === 'adjacent_local' ? 'Collinear group ID' : 'Orthogroup ID';
  const recordCoverage = firstText(group?.record_coverage_count, group?.recordCoverage);
  const rbhOrthogroups = Array.isArray(group?.rbhOrthogroupIds) ? group.rbhOrthogroupIds : [];
  const orthologPathCount = Array.isArray(group?.orthologPaths) ? String(group.orthologPaths.length) : '';
  const relatedEdgeCount = Array.isArray(group?.relatedEdges) ? String(group.relatedEdges.length) : '';
  const memberRows = buildOrthogroupMemberRows(group, featureLookup, orthogroupId);
  return {
    id: orthogroupId,
    displayName,
    description,
    memberCount,
    recordCoverage,
    queryMember: resolveBlockMemberLabels({ group, featureSvgIds: queryFeatureSvgId, featureLookup }),
    subjectMember: resolveBlockMemberLabels({ group, featureSvgIds: subjectFeatureSvgId, featureLookup }),
    detailRows: buildOrthogroupDetailRows({
      orthogroupId,
      idLabel,
      displayName,
      description,
      scopeLabel,
      memberCount,
      recordCoverage,
      rbhOrthogroups,
      orthologPathCount,
      relatedEdgeCount
    }),
    ...orthogroupMemberSectionExtras(memberRows, orthogroupId, displayName)
  };
});

const buildFeatureRows = ({
  title,
  feature,
  recordId,
  interval,
  proteinId,
  locusId,
  displayName,
  featureSvgIds,
  featureLookup,
  group
}) => {
  const rows = [];
  const featureRows = buildFeatureListRows({
    featureSvgIds,
    featureLookup,
    group,
    recordId,
    interval,
    proteinId,
    locusId,
    displayName
  });
  const displayProteinIds = resolveFeatureSectionProteinIds({
    feature,
    featureSvgIds,
    featureLookup,
    group,
    fallbackProteinIds: proteinId,
    locusId,
    displayName
  });
  addRow(rows, 'Record', firstText(feature?.record_id, recordId));
  addRow(rows, 'Location', firstText(featureLocationText(feature), interval));
  addRow(rows, displayProteinIds.includes(';') ? 'Protein IDs' : 'Protein ID', displayProteinIds);
  addRow(rows, 'Gene', firstText(feature?.gene, qualifierFirstValue(feature, 'gene')));
  addRow(rows, 'Locus tag', firstText(feature?.locus_tag, feature?.locusTag, qualifierFirstValue(feature, 'locus_tag')));
  addRow(rows, 'Product', firstText(feature?.product, qualifierFirstValue(feature, 'product')));
  addRow(rows, 'Locus ID', locusId);
  addRow(rows, 'Display name', displayName);
  return section(title, rows, { featureRows });
};

export const buildPairwiseMatchPayload = (
  element,
  {
    featureLookup = new Map(),
    orthogroups = [],
    orthogroupNameOverrides = null,
    orthogroupDescriptionOverrides = null
  } = {}
) => {
  if (!element) return null;
  const matchKind = normalizeMatchKind(attr(element, 'data-match-kind'), element);
  const orthogroupId = attr(element, 'data-orthogroup-id');
  const orthogroupIds = uniqueMetadataValues(orthogroupId);
  const collinearityBlockId = attr(element, 'data-collinearity-block-id');
  const groupScope = firstText(attr(element, 'data-collinear-group-scope'), attr(element, 'data-group-scope'));
  const normalizedGroupScope = normalizeGroupMetadataScope(groupScope);
  const queryFeatureSvgId = attr(element, 'data-query-feature-svg-id');
  const subjectFeatureSvgId = attr(element, 'data-subject-feature-svg-id');
  const queryFeature = featureLookup.get?.(queryFeatureSvgId) || null;
  const subjectFeature = featureLookup.get?.(subjectFeatureSvgId) || null;
  const group = matchKind === 'collinear'
    ? null
    : getOrthogroupForMatch(orthogroups, {
      orthogroupId,
      queryFeatureSvgId,
      subjectFeatureSvgId
    }) || buildFallbackOrthogroup({
      orthogroupId,
      queryFeature,
      subjectFeature,
      featureLookup
    });
  const blockOrthogroups = matchKind === 'collinear'
    ? buildBlockOrthogroups({
      orthogroupIds,
      orthogroups,
      featureLookup,
      queryFeatureSvgId,
      subjectFeatureSvgId,
      orthogroupNameOverrides,
      orthogroupDescriptionOverrides,
      groupScope
    })
    : [];
  const displayName = firstText(
    overrideValue(orthogroupNameOverrides, orthogroupId),
    group?.displayName,
    group?.display_name,
    group?.name,
    featureProduct(queryFeature),
    featureProduct(subjectFeature)
  );
  const description = firstText(
    overrideValue(orthogroupDescriptionOverrides, orthogroupId),
    group?.description
  );
  const qInterval = intervalText(attr(element, 'data-qstart'), attr(element, 'data-qend'));
  const sInterval = intervalText(attr(element, 'data-sstart'), attr(element, 'data-send'));
  const title = MATCH_KIND_TITLES[matchKind] || MATCH_KIND_TITLES.pairwise;
  const subtitle = firstText(collinearityBlockId, orthogroupId, attr(element, 'data-gbdraw-pairwise-match-id'));

  const summaryRows = [];
  addRow(summaryRows, 'Query record', attr(element, 'data-query-record-id'));
  addRow(summaryRows, 'Subject record', attr(element, 'data-subject-record-id'));
  addRow(summaryRows, 'Query interval', qInterval);
  addRow(summaryRows, 'Subject interval', sInterval);
  addRow(summaryRows, 'Orientation', firstText(
    attr(element, 'data-collinearity-orientation'),
    attr(element, 'data-orientation')
  ));

  const alignmentRows = [];
  addRow(alignmentRows, 'Identity', attr(element, 'data-identity'));
  addRow(alignmentRows, 'Alignment length', attr(element, 'data-alignment-length'));
  addRow(alignmentRows, 'E-value', attr(element, 'data-evalue'));
  addRow(alignmentRows, 'Bit score', attr(element, 'data-bitscore'));
  addRow(alignmentRows, 'Mismatches', attr(element, 'data-mismatches'));
  addRow(alignmentRows, 'Gap opens', attr(element, 'data-gap-opens'));
  addRow(alignmentRows, 'Edge kind', attr(element, 'data-edge-kind'));
  addRow(alignmentRows, 'Render role', attr(element, 'data-render-role'));
  addRow(alignmentRows, 'RBH seed', attr(element, 'data-rbh-orthogroup-id'));
  addRow(alignmentRows, 'Path ID', attr(element, 'data-ortholog-path-id'));
  addRow(alignmentRows, 'Query members', attr(element, 'data-query-orthogroup-member-count'));
  addRow(alignmentRows, 'Subject members', attr(element, 'data-subject-orthogroup-member-count'));
  addRow(alignmentRows, 'Query member role', attr(element, 'data-query-orthogroup-role'));
  addRow(alignmentRows, 'Subject member role', attr(element, 'data-subject-orthogroup-role'));
  addRow(alignmentRows, 'Query confidence', attr(element, 'data-query-orthogroup-confidence'));
  addRow(alignmentRows, 'Subject confidence', attr(element, 'data-subject-orthogroup-confidence'));
  addRow(alignmentRows, 'Query assignment', attr(element, 'data-query-orthogroup-assignment-reason'));
  addRow(alignmentRows, 'Subject assignment', attr(element, 'data-subject-orthogroup-assignment-reason'));

  const blockRows = [];
  addRow(blockRows, 'Block ID', collinearityBlockId);
  addRow(blockRows, 'Kind', attr(element, 'data-collinearity-block-kind'));
  addRow(blockRows, 'Orientation', attr(element, 'data-collinearity-orientation'));
  addRow(blockRows, 'Color mode', attr(element, 'data-collinearity-color-mode'));
  if (matchKind === 'collinear') {
    addRow(blockRows, 'Average identity', attr(element, 'data-identity'));
    addRow(blockRows, 'Aligned length', attr(element, 'data-alignment-length'));
  }
  addRow(blockRows, 'Block score', attr(element, 'data-collinearity-block-score'));
  addRow(blockRows, 'Block e-value', attr(element, 'data-collinearity-block-evalue'));
  addRow(blockRows, 'Anchor', [
    attr(element, 'data-collinearity-anchor-index'),
    attr(element, 'data-collinearity-anchor-count')
  ].filter(Boolean).join(' / '));

  const orthogroupRows = [];
  if (matchKind !== 'collinear') {
    addRow(orthogroupRows, 'Orthogroup ID', orthogroupId);
    addRow(orthogroupRows, 'Display name', displayName);
    addRow(orthogroupRows, 'Description', description);
    addRow(orthogroupRows, 'Members', firstText(group?.member_count, group?.memberCount));
    addRow(orthogroupRows, 'Record coverage', firstText(group?.record_coverage_count, group?.recordCoverage));
    addRow(orthogroupRows, 'RBH seeds', Array.isArray(group?.rbhOrthogroupIds) ? group.rbhOrthogroupIds.join('; ') : '');
    addRow(orthogroupRows, 'Ortholog paths', Array.isArray(group?.orthologPaths) ? String(group.orthologPaths.length) : '');
    addRow(orthogroupRows, 'Related edges', Array.isArray(group?.relatedEdges) ? String(group.relatedEdges.length) : '');
  }
  const orthogroupMemberRows = buildOrthogroupMemberRows(group, featureLookup, orthogroupId);

  if (matchKind === 'orthogroup') {
    const summarySection = section(
      'Summary',
      orthogroupRows,
      orthogroupMemberSectionExtras(orthogroupMemberRows, orthogroupId, displayName)
    );
    return {
      id: firstText(attr(element, 'data-gbdraw-pairwise-match-id'), orthogroupId),
      title: orthogroupTitle(orthogroupId, displayName),
      subtitle: '',
      matchKind,
      orthogroupId,
      collinearityBlockId,
      queryFeatureSvgId,
      subjectFeatureSvgId,
      fill: firstText(element.getAttribute('fill'), element.style?.fill, '#94a3b8'),
      sections: summarySection.rows.length > 0 || summarySection.memberRows.length > 0
        ? [summarySection]
        : []
    };
  }

  const sections = [section('Summary', summaryRows)];
  if (matchKind !== 'collinear') {
    sections.push(section('Alignment', alignmentRows));
  }
  if (matchKind === 'orthogroup' || orthogroupRows.length > 0) {
    sections.push(section(
      'Orthogroup',
      orthogroupRows,
      orthogroupMemberSectionExtras(orthogroupMemberRows, orthogroupId, displayName)
    ));
  }
  if (matchKind === 'collinear') {
    const blockOrthogroupRows = [];
    const localGroups = normalizedGroupScope === 'adjacent_local';
    addRow(
      blockOrthogroupRows,
      localGroups ? 'Number of local collinear groups' : 'Number of orthogroups covered',
      String(orthogroupIds.length)
    );
    sections.push(section(localGroups ? 'Local collinear groups' : 'Orthogroups covered', blockOrthogroupRows, { blockOrthogroups }));
  }
  if (matchKind === 'collinear' || blockRows.length > 0) {
    sections.push(section('Collinearity', blockRows));
  }
  sections.push(buildFeatureRows({
    title: 'Query',
    feature: queryFeature,
    recordId: attr(element, 'data-query-record-id'),
    interval: qInterval,
    proteinId: attr(element, 'data-query-protein-id'),
    locusId: attr(element, 'data-query-locus-id'),
    displayName: attr(element, 'data-query-display-name'),
    featureSvgIds: queryFeatureSvgId,
    featureLookup,
    group
  }));
  sections.push(buildFeatureRows({
    title: 'Subject',
    feature: subjectFeature,
    recordId: attr(element, 'data-subject-record-id'),
    interval: sInterval,
    proteinId: attr(element, 'data-subject-protein-id'),
    locusId: attr(element, 'data-subject-locus-id'),
    displayName: attr(element, 'data-subject-display-name'),
    featureSvgIds: subjectFeatureSvgId,
    featureLookup,
    group
  }));

  return {
    id: firstText(attr(element, 'data-gbdraw-pairwise-match-id'), collinearityBlockId, orthogroupId),
    title,
    subtitle,
    matchKind,
    orthogroupId,
    groupScope: normalizedGroupScope,
    collinearityBlockId,
    queryFeatureSvgId,
    subjectFeatureSvgId,
    blockOrthogroupCount: blockOrthogroups.length || (matchKind === 'collinear' ? orthogroupIds.length : 0),
    blockOrthogroups,
    fill: firstText(element.getAttribute('fill'), element.style?.fill, '#94a3b8'),
    sections: sections.filter((entry) => (
      entry.rows.length > 0 ||
      (Array.isArray(entry.featureRows) && entry.featureRows.length > 0)
    ))
  };
};

export const buildPairwiseMatchHoverRows = (payload) => {
  const rows = [];
  const addFirst = (label, value) => addRow(rows, label, value);
  if (!payload) return rows;
  const summary = payload.sections.find((entry) => entry.title === 'Summary');
  const findValue = (sectionEntry, label) =>
    sectionEntry?.rows.find((row) => row.label === label)?.value || '';
  if (payload.matchKind === 'orthogroup') {
    addFirst('Kind', payload.matchKind);
    addFirst('Orthogroup', findValue(summary, 'Orthogroup ID') || payload.orthogroupId);
    addFirst('Display name', findValue(summary, 'Display name'));
    addFirst('Members', findValue(summary, 'Members'));
    return rows.slice(0, 6);
  }
  const alignment = payload.sections.find((entry) => entry.title === 'Alignment');
  const block = payload.sections.find((entry) => entry.title === 'Collinearity');
  const orthogroup = payload.sections.find((entry) => entry.title === 'Orthogroup');
  addFirst('Kind', payload.matchKind);
  addFirst('Identity', findValue(alignment, 'Identity') || findValue(block, 'Average identity'));
  addFirst('Query', findValue(summary, 'Query interval'));
  addFirst('Subject', findValue(summary, 'Subject interval'));
  if (payload.matchKind === 'collinear') {
    addFirst(payload.groupScope === 'adjacent_local' ? 'Collinear groups' : 'Orthogroups', String(payload.blockOrthogroupCount ?? ''));
  } else {
    addFirst('Orthogroup', findValue(orthogroup, 'Orthogroup ID'));
  }
  addFirst('Block', findValue(block, 'Block ID'));
  return rows.slice(0, 6);
};
