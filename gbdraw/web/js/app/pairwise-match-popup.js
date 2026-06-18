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

const attr = (element, name) => normalizeText(element?.getAttribute?.(name));

const firstText = (...values) => {
  for (const value of values) {
    const text = normalizeText(value);
    if (text) return text;
  }
  return '';
};

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

const proteinIdText = (feature, fallback) => firstText(
  fallback,
  feature?.sourceProteinId,
  feature?.source_protein_id,
  feature?.proteinId,
  feature?.protein_id
);

const getOrthogroupById = (orthogroups, orthogroupId) => {
  const id = normalizeText(orthogroupId);
  if (!id) return null;
  return (Array.isArray(orthogroups) ? orthogroups : [])
    .find((entry) => normalizeText(entry?.id || entry?.orthogroupId || entry?.orthogroup_id) === id) || null;
};

const overrideValue = (overrides, key) => {
  const normalizedKey = normalizeText(key);
  if (!normalizedKey || !overrides) return '';
  if (overrides instanceof Map) return normalizeText(overrides.get(normalizedKey));
  return normalizeText(overrides[normalizedKey]);
};

const section = (title, rows) => ({ title, rows: rows.filter((row) => normalizeText(row.value)) });

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
  addRow(rows, 'Protein ID', proteinIdText(feature, proteinId));
  addRow(rows, 'Source protein ID', firstText(feature?.sourceProteinId, feature?.source_protein_id));
  addRow(rows, 'Gene', firstText(feature?.gene, qualifierFirstValue(feature, 'gene')));
  addRow(rows, 'Locus tag', firstText(feature?.locus_tag, feature?.locusTag, qualifierFirstValue(feature, 'locus_tag')));
  addRow(rows, 'Product', firstText(feature?.product, qualifierFirstValue(feature, 'product')));
  addRow(rows, 'Unit ID', unitId);
  addRow(rows, 'Locus ID', locusId);
  addRow(rows, 'Display name', displayName);
  addRow(rows, 'Feature SVG ID', svgId);
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
  const group = getOrthogroupById(orthogroups, orthogroupId);
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
  addRow(summaryRows, 'Match style', attr(element, 'data-pairwise-match-style'));
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

  const sections = [
    section('Summary', summaryRows),
    section('Alignment', alignmentRows)
  ];
  if (matchKind === 'orthogroup' || orthogroupRows.length > 0) {
    sections.push(section('Orthogroup', orthogroupRows));
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
  const alignment = payload.sections.find((entry) => entry.title === 'Alignment');
  const block = payload.sections.find((entry) => entry.title === 'Collinearity');
  const orthogroup = payload.sections.find((entry) => entry.title === 'Orthogroup');
  const findValue = (sectionEntry, label) =>
    sectionEntry?.rows.find((row) => row.label === label)?.value || '';
  addFirst('Kind', payload.matchKind);
  addFirst('Identity', findValue(alignment, 'Identity'));
  addFirst('Query', findValue(summary, 'Query interval'));
  addFirst('Subject', findValue(summary, 'Subject interval'));
  addFirst('Orthogroup', findValue(orthogroup, 'Orthogroup ID'));
  addFirst('Block', findValue(block, 'Block ID'));
  return rows.slice(0, 6);
};
