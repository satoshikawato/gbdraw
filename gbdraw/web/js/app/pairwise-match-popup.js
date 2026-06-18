import { resolveDisplayProteinId } from './feature-utils.js';
import { buildFeatureSequenceFastas } from './feature-sequence-fasta.js';

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
  const direct = featureLookup.get?.(svgId) || null;
  if (direct && (recordIndex === null || featureRecordIndex(direct) === recordIndex)) return direct;
  return getFeatureLookupValues(featureLookup).find((feature) => {
    if (normalizeText(feature?.svg_id || feature?.svgId) !== svgId) return false;
    return recordIndex === null || featureRecordIndex(feature) === recordIndex;
  }) || direct || null;
};

const groupHasFeatureSvgId = (group, featureSvgId) => {
  const id = normalizeText(featureSvgId);
  if (!id) return false;
  return (Array.isArray(group?.members) ? group.members : [])
    .some((member) => memberFeatureSvgId(member) === id);
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
        productOrNote: firstText(member?.product, member?.note),
        ntFasta: memberFastaText(member, feature, nucleotideSequence, aminoAcidSequence, 'nt'),
        aaFasta: memberFastaText(member, feature, nucleotideSequence, aminoAcidSequence, 'aa'),
        ntFilename: orthogroupMemberSequenceFilename(member, orthogroupId, 'nt'),
        aaFilename: orthogroupMemberSequenceFilename(member, orthogroupId, 'aa')
      };
    })
    .filter((row) => row.record || row.coordinates || row.proteinId || row.productOrNote || row.ntFasta || row.aaFasta);
};

const orthogroupMemberCopyText = (memberRows) => {
  if (!Array.isArray(memberRows) || memberRows.length === 0) return '';
  return [
    'Record\tCoordinates (+/-)\tProtein ID\tProduct / note',
    ...memberRows.map((row) => [row.record, row.coordinates, row.proteinId, row.productOrNote].join('\t'))
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

const buildFeatureRows = ({
  title,
  feature,
  recordId,
  interval,
  proteinId,
  unitId,
  locusId,
  displayName,
  svgId
}) => {
  const rows = [];
  addRow(rows, 'Record', firstText(feature?.record_id, recordId));
  addRow(rows, 'Location', firstText(featureLocationText(feature), interval));
  addRow(rows, 'Protein ID', resolveDisplayProteinId(feature, null, proteinId));
  addRow(rows, 'Gene', firstText(feature?.gene, qualifierFirstValue(feature, 'gene')));
  addRow(rows, 'Locus tag', firstText(feature?.locus_tag, feature?.locusTag, qualifierFirstValue(feature, 'locus_tag')));
  addRow(rows, 'Product', firstText(feature?.product, qualifierFirstValue(feature, 'product')));
  addRow(rows, 'Unit ID', unitId);
  addRow(rows, 'Locus ID', locusId);
  addRow(rows, 'Display name', displayName);
  return section(title, rows);
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
  const collinearityBlockId = attr(element, 'data-collinearity-block-id');
  const queryFeatureSvgId = attr(element, 'data-query-feature-svg-id');
  const subjectFeatureSvgId = attr(element, 'data-subject-feature-svg-id');
  const queryFeature = featureLookup.get?.(queryFeatureSvgId) || null;
  const subjectFeature = featureLookup.get?.(subjectFeatureSvgId) || null;
  const group = getOrthogroupForMatch(orthogroups, {
    orthogroupId,
    queryFeatureSvgId,
    subjectFeatureSvgId
  }) || buildFallbackOrthogroup({
    orthogroupId,
    queryFeature,
    subjectFeature,
    featureLookup
  });
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
  addRow(summaryRows, 'Match ID', attr(element, 'data-gbdraw-pairwise-match-id'));
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

  const blockRows = [];
  addRow(blockRows, 'Block ID', collinearityBlockId);
  addRow(blockRows, 'Kind', attr(element, 'data-collinearity-block-kind'));
  addRow(blockRows, 'Orientation', attr(element, 'data-collinearity-orientation'));
  addRow(blockRows, 'Color mode', attr(element, 'data-collinearity-color-mode'));
  addRow(blockRows, 'Block score', attr(element, 'data-collinearity-block-score'));
  addRow(blockRows, 'Block e-value', attr(element, 'data-collinearity-block-evalue'));
  addRow(blockRows, 'Anchor', [
    attr(element, 'data-collinearity-anchor-index'),
    attr(element, 'data-collinearity-anchor-count')
  ].filter(Boolean).join(' / '));
  addRow(blockRows, 'Query unit', attr(element, 'data-query-unit-id'));
  addRow(blockRows, 'Subject unit', attr(element, 'data-subject-unit-id'));
  addRow(blockRows, 'Query display', attr(element, 'data-query-display-name'));
  addRow(blockRows, 'Subject display', attr(element, 'data-subject-display-name'));

  const orthogroupRows = [];
  addRow(orthogroupRows, 'Orthogroup ID', orthogroupId);
  addRow(orthogroupRows, 'Display name', displayName);
  addRow(orthogroupRows, 'Description', description);
  addRow(orthogroupRows, 'Members', firstText(group?.member_count, group?.memberCount));
  addRow(orthogroupRows, 'Record coverage', firstText(group?.record_coverage_count, group?.recordCoverage));
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

  const sections = [
    section('Summary', summaryRows),
    section('Alignment', alignmentRows)
  ];
  if (matchKind === 'orthogroup' || orthogroupRows.length > 0) {
    sections.push(section(
      'Orthogroup',
      orthogroupRows,
      orthogroupMemberSectionExtras(orthogroupMemberRows, orthogroupId, displayName)
    ));
  }
  if (matchKind === 'collinear' || blockRows.length > 0) {
    sections.push(section('Collinearity', blockRows));
  }
  sections.push(buildFeatureRows({
    title: 'Query feature',
    feature: queryFeature,
    recordId: attr(element, 'data-query-record-id'),
    interval: qInterval,
    proteinId: attr(element, 'data-query-protein-id'),
    unitId: attr(element, 'data-query-unit-id'),
    locusId: attr(element, 'data-query-locus-id'),
    displayName: attr(element, 'data-query-display-name'),
    svgId: queryFeatureSvgId
  }));
  sections.push(buildFeatureRows({
    title: 'Subject feature',
    feature: subjectFeature,
    recordId: attr(element, 'data-subject-record-id'),
    interval: sInterval,
    proteinId: attr(element, 'data-subject-protein-id'),
    unitId: attr(element, 'data-subject-unit-id'),
    locusId: attr(element, 'data-subject-locus-id'),
    displayName: attr(element, 'data-subject-display-name'),
    svgId: subjectFeatureSvgId
  }));

  return {
    id: firstText(attr(element, 'data-gbdraw-pairwise-match-id'), collinearityBlockId, orthogroupId),
    title,
    subtitle,
    matchKind,
    orthogroupId,
    collinearityBlockId,
    queryFeatureSvgId,
    subjectFeatureSvgId,
    fill: firstText(element.getAttribute('fill'), element.style?.fill, '#94a3b8'),
    sections: sections.filter((entry) => entry.rows.length > 0)
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
  addFirst('Identity', findValue(alignment, 'Identity'));
  addFirst('Query', findValue(summary, 'Query interval'));
  addFirst('Subject', findValue(summary, 'Subject interval'));
  addFirst('Orthogroup', findValue(orthogroup, 'Orthogroup ID'));
  addFirst('Block', findValue(block, 'Block ID'));
  return rows.slice(0, 6);
};
