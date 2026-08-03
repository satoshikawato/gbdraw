import { normalizeStringArray } from '../app/feature-utils.js';
import { STANDALONE_INTERACTIVE_SCRIPT, STANDALONE_INTERACTIVE_STYLE } from './standalone-interactivity-assets.js';
import { ensureSvgDefs } from './svg-serialization.js';

const SVG_NS = 'http://www.w3.org/2000/svg';
const FEATURE_ID_ATTRIBUTE = 'data-gbdraw-feature-id';
const RENDERED_FEATURE_ID_ATTRIBUTE = 'data-gbdraw-rendered-feature-id';
const FEATURE_SELECTOR = [
  `path[${FEATURE_ID_ATTRIBUTE}]`,
  `polygon[${FEATURE_ID_ATTRIBUTE}]`,
  `rect[${FEATURE_ID_ATTRIBUTE}]`,
  'path[id^="f"]',
  'polygon[id^="f"]',
  'rect[id^="f"]'
].join(', ');
const STANDALONE_MATCH_SELECTOR = [
  'path[data-gbdraw-match-id]',
  'path[data-gbdraw-pairwise-match-id]',
  'path[data-match-kind]',
  'path[data-pairwise-match-style]'
].join(', ');
const STANDALONE_ANNOTATION_SELECTOR = '[data-gbdraw-annotation-id]';
const INTERACTIVE_METADATA_ID = 'gbdraw-interactive-feature-metadata';
const INTERACTIVE_STYLE_ID = 'gbdraw-interactive-feature-style';
const INTERACTIVE_SCRIPT_ID = 'gbdraw-interactive-feature-script';
const INTERACTIVE_GLOW_FILTER_ID = 'gbdraw-interactive-feature-glow';
const INTERACTIVE_MATCH_GLOW_FILTER_ID = 'gbdraw-interactive-feature-match-glow';
const FEATURE_PART_SUFFIX_RE = /__part\d+$/;
const FEATURE_RECORD_SUFFIX_RE = /_record_\d+$/;
const INTERACTIVE_SCHEMA = 'gbdraw-interactive-feature-popup-v2';
const INTERACTIVE_CATALOG_SCHEMA = 3;

const compactWireValue = (value) => {
  if (Array.isArray(value)) {
    const items = value.map(compactWireValue).filter((item) => item !== undefined);
    return items.length ? items : undefined;
  }
  if (value && typeof value === 'object') {
    const compact = {};
    Object.entries(value).forEach(([key, entry]) => {
      const normalized = compactWireValue(entry);
      if (normalized !== undefined) compact[key] = normalized;
    });
    return Object.keys(compact).length ? compact : undefined;
  }
  if (value === null || value === undefined || value === '' || value === false) return undefined;
  return value;
};

const normalizeFeatureElementId = (value) =>
  String(value || '').trim().replace(FEATURE_PART_SUFFIX_RE, '');

export const stripEditorOnlyCursorStyles = (svg) => {
  if (!svg) return;
  svg.querySelectorAll('[style]').forEach((element) => {
    const style = element.getAttribute('style');
    if (!style || !/\bcursor\s*:/i.test(style)) return;
    element.style.removeProperty('cursor');
    if (!element.getAttribute('style')?.trim()) {
      element.removeAttribute('style');
    }
  });
};

const getElementFeatureId = (element) =>
  normalizeFeatureElementId(
    element?.getAttribute?.(RENDERED_FEATURE_ID_ATTRIBUTE) ||
    element?.getAttribute?.(FEATURE_ID_ATTRIBUTE) ||
    element?.getAttribute?.('id') ||
    ''
  );



const normalizeQualifierMap = (qualifiers) => {
  const normalized = {};
  if (!qualifiers || typeof qualifiers !== 'object' || Array.isArray(qualifiers)) {
    return normalized;
  }
  Object.entries(qualifiers).forEach(([key, value]) => {
    const normalizedKey = String(key || '').trim();
    const values = normalizeStringArray(value);
    if (!normalizedKey || values.length === 0) return;
    normalized[normalizedKey] = values;
  });
  return normalized;
};

const normalizeLocationParts = (parts) => {
  if (!Array.isArray(parts)) return [];
  return parts
    .map((part) => {
      const start = Number(part?.start);
      const end = Number(part?.end);
      const display = String(part?.display || '').trim() ||
        `${Number.isFinite(start) ? start + 1 : ''}..${Number.isFinite(end) ? end : ''}`;
      return {
        start: Number.isFinite(start) ? start : null,
        end: Number.isFinite(end) ? end : null,
        strand: String(part?.strand || '').trim(),
        display
      };
    })
    .filter((part) => part.display && part.display !== '..');
};

const buildStandaloneFeatureLocation = (feature) => {
  const start = Number(feature?.start);
  const end = Number(feature?.end);
  const startText = Number.isFinite(start) ? String(start + 1) : String(feature?.start ?? '');
  const endText = Number.isFinite(end) ? String(end) : String(feature?.end ?? '');
  const strand = String(feature?.strand || '').trim();
  const range = `${startText}..${endText}`;
  return strand ? `${range} (${strand})` : range;
};

const firstQualifierValue = (feature, key) => {
  const qualifiers = feature?.qualifiers && typeof feature.qualifiers === 'object'
    ? feature.qualifiers
    : {};
  const values = normalizeStringArray(qualifiers[key]);
  return values.find((value) => value.trim()) || '';
};

const getStandaloneFeatureLabel = (feature) => {
  const candidates = [
    feature?.label,
    feature?.gene,
    feature?.locus_tag,
    firstQualifierValue(feature, 'gene'),
    firstQualifierValue(feature, 'locus_tag'),
    firstQualifierValue(feature, 'product'),
    feature?.product,
    feature?.type,
    feature?.svg_id
  ];
  for (const candidate of candidates) {
    const label = String(candidate || '').trim();
    if (label) return label;
  }
  return 'Feature';
};

const normalizeFeatureIdKey = (value) => String(value || '').trim().toLowerCase();

const standaloneConsistentTextAlias = (payload, keys) => {
  const values = new Set(keys
    .filter((key) => Object.prototype.hasOwnProperty.call(payload || {}, key))
    .map((key) => String(payload?.[key] ?? '').trim())
    .filter(Boolean));
  return {
    valid: values.size <= 1,
    value: values.size === 1 ? values.values().next().value : ''
  };
};

const standaloneIntegerAlias = (candidates) => {
  const indexes = new Set();
  for (const candidate of candidates) {
    if (candidate === null || candidate === undefined || candidate === '') continue;
    const text = String(candidate).trim();
    if (typeof candidate === 'boolean' || !/^\d+$/.test(text)) {
      return { valid: false, value: null };
    }
    const recordIndex = Number(text);
    if (!Number.isSafeInteger(recordIndex)) {
      return { valid: false, value: null };
    }
    indexes.add(recordIndex);
  }
  return {
    valid: indexes.size <= 1,
    value: indexes.size === 1 ? indexes.values().next().value : null
  };
};

const standaloneRecordIndex = (candidates) => {
  const status = standaloneIntegerAlias(candidates);
  return status.valid ? status.value : null;
};

const standaloneFeatureRecordIndexStatus = (feature) => (
  standaloneIntegerAlias([
    feature?.fileIdx,
    feature?.file_idx,
    feature?.recordIndex,
    feature?.record_index,
    feature?.record_idx
  ])
);

const standaloneFeatureRecordIndex = (feature) => {
  const status = standaloneFeatureRecordIndexStatus(feature);
  return status.valid ? status.value : null;
};

const standaloneSourceIdentityStatus = (payload) => {
  const stable = standaloneConsistentTextAlias(payload, [
    'stable_svg_id', 'stableSvgId', 'stableFeatureSvgId',
    'stable_feature_svg_id', 'stable_feature_id', 'stableFeatureId'
  ]);
  const legacy = standaloneConsistentTextAlias(payload, ['featureSvgId', 'feature_svg_id']);
  return {
    valid: stable.valid && legacy.valid && (
      !stable.value || !legacy.value || stable.value === legacy.value
    ),
    value: stable.value || legacy.value
  };
};

const standaloneRenderedIdentityStatus = (payload) => standaloneConsistentTextAlias(
  payload,
  [
    'rendered_svg_id', 'renderedSvgId', 'renderedFeatureSvgId',
    'rendered_feature_svg_id'
  ]
);

const standaloneFeatureIdentityStatus = (feature) => {
  const source = standaloneSourceIdentityStatus(feature);
  const rendered = standaloneRenderedIdentityStatus(feature);
  const svg = standaloneConsistentTextAlias(feature, ['svg_id', 'svgId']);
  if (!source.valid || !rendered.valid || !svg.valid) {
    return { valid: false, source, rendered };
  }
  if (svg.value && rendered.value && svg.value !== rendered.value) {
    return { valid: false, source, rendered };
  }
  const renderedExplicit = Boolean(rendered.value);
  if (svg.value && !rendered.value) rendered.value = svg.value;
  if (svg.value && !source.value) {
    source.value = svg.value.replace(FEATURE_RECORD_SUFFIX_RE, '');
  }
  return { valid: true, source, rendered, renderedExplicit };
};

const standaloneFeatureStableSvgId = (feature) => {
  const status = standaloneFeatureIdentityStatus(feature);
  return status.valid ? status.source.value : '';
};

const standaloneMemberRecordIndexStatus = (member) => (
  standaloneIntegerAlias([
    member?.recordIndex,
    member?.record_index,
    member?.record_idx,
    member?.fileIdx,
    member?.file_idx
  ])
);

const standaloneMemberRecordIndex = (member) => {
  const status = standaloneMemberRecordIndexStatus(member);
  return status.valid ? status.value : null;
};

const standaloneMemberStableSvgId = (member) => {
  const stable = standaloneConsistentTextAlias(member, [
    'stableFeatureSvgId', 'stable_feature_svg_id',
    'stableFeatureId', 'stable_feature_id'
  ]);
  const legacy = standaloneConsistentTextAlias(
    member,
    ['featureSvgId', 'feature_svg_id']
  );
  if (
    !stable.valid
    || !legacy.valid
    || (stable.value && legacy.value && stable.value !== legacy.value)
  ) return '';
  return stable.value || legacy.value;
};

const standaloneMemberRenderedSvgId = (member) => {
  const status = standaloneConsistentTextAlias(member, [
    'renderedFeatureSvgId', 'rendered_feature_svg_id',
    'renderedSvgId', 'rendered_svg_id'
  ]);
  return status.valid ? status.value : '';
};

const standaloneFeatureIndexStatus = (payload) => standaloneIntegerAlias([
  payload?.featureIndex,
  payload?.feature_index,
  payload?.sourceFeatureIndex,
  payload?.source_feature_index
]);

const standaloneCanonicalIdentityStatus = (payload) => {
  const record = standaloneConsistentTextAlias(payload, ['recordKey', 'record_key']);
  const feature = standaloneConsistentTextAlias(
    payload,
    ['biologicalFeatureId', 'biological_feature_id']
  );
  return {
    valid: record.valid && feature.valid && Boolean(record.value) === Boolean(feature.value),
    value: record.value && feature.value ? [record.value, feature.value] : null
  };
};

const normalizeStandaloneContext = (options = {}) => ({
  popupMode: normalizeStandalonePopupMode(options.popupMode),
  featureCatalog:
    options.featureCatalog && typeof options.featureCatalog === 'object'
      ? options.featureCatalog
      : null,
  catalogResultIndex: options.catalogResultIndex ?? null,
  catalogResultName: String(options.catalogResultName || '').trim(),
  requireFeatureCatalog: options.requireFeatureCatalog === true,
  features: Array.isArray(options.features) ? options.features : [],
  editableLabels: Array.isArray(options.editableLabels) ? options.editableLabels : [],
  labelTextFeatureOverrides:
    options.labelTextFeatureOverrides && typeof options.labelTextFeatureOverrides === 'object'
      ? options.labelTextFeatureOverrides
      : {},
  labelTextBulkOverrides:
    options.labelTextBulkOverrides && typeof options.labelTextBulkOverrides === 'object'
      ? options.labelTextBulkOverrides
      : {},
  featureOrthogroupIndex: options.featureOrthogroupIndex instanceof Map
    ? options.featureOrthogroupIndex
    : new Map(),
  orthogroups: Array.isArray(options.orthogroups) ? options.orthogroups : [],
  orthogroupNameOverrides:
    options.orthogroupNameOverrides && typeof options.orthogroupNameOverrides === 'object'
      ? options.orthogroupNameOverrides
      : {},
  orthogroupDescriptionOverrides:
    options.orthogroupDescriptionOverrides && typeof options.orthogroupDescriptionOverrides === 'object'
      ? options.orthogroupDescriptionOverrides
      : {},
  legendEntries: Array.isArray(options.legendEntries) ? options.legendEntries : [],
  currentColors:
    options.currentColors && typeof options.currentColors === 'object'
      ? options.currentColors
      : {},
  sequenceSources: Array.isArray(options.sequenceSources) ? options.sequenceSources : []
});

const selectStandaloneCatalogItem = (context) => {
  const catalog = context?.featureCatalog;
  if (
    !catalog
    || catalog.schema !== INTERACTIVE_CATALOG_SCHEMA
    || !Array.isArray(catalog.items)
  ) {
    return null;
  }
  const resultIndex = standaloneIntegerAlias([context.catalogResultIndex]);
  const hasResultIndex = context.catalogResultIndex !== null
    && context.catalogResultIndex !== undefined
    && String(context.catalogResultIndex).trim() !== '';
  if (hasResultIndex && (!resultIndex.valid || resultIndex.value === null)) {
    return null;
  }
  const hasResultName = Boolean(context.catalogResultName);
  if (!hasResultIndex && !hasResultName) {
    return catalog.items.length === 1 ? catalog.items[0] : null;
  }
  const matches = catalog.items.filter((item) => (
    item
    && (!hasResultIndex || item.resultIndex === resultIndex.value)
    && (!hasResultName || String(item.resultName || '').trim() === context.catalogResultName)
  ));
  return matches.length === 1 ? matches[0] : null;
};

const cloneStandaloneCatalogItem = (item) => JSON.parse(JSON.stringify(item));

const standaloneCatalogFeatureKey = (recordKey, biologicalFeatureId) => (
  `${String(recordKey || '').trim()}\u0000${String(biologicalFeatureId || '').trim()}`
);

const getEditableLabelEntryForStandaloneFeature = (feature, context) => {
  const entries = Array.isArray(context?.editableLabels) ? context.editableLabels : [];
  const candidates = [
    normalizeFeatureIdKey(feature?.svg_id),
    normalizeFeatureIdKey(feature?.id)
  ].filter(Boolean);
  if (candidates.length === 0) return null;
  return entries.find((entry) => candidates.includes(normalizeFeatureIdKey(entry?.featureId))) || null;
};

const getLabelTextFeatureOverride = (feature, context) => {
  const overrides = context?.labelTextFeatureOverrides || {};
  const candidates = [
    String(feature?.svg_id || '').trim(),
    String(feature?.id || '').trim()
  ].filter(Boolean);
  for (const candidate of candidates) {
    if (Object.prototype.hasOwnProperty.call(overrides, candidate)) {
      const text = String(overrides[candidate] ?? '').trim();
      if (text) return text;
    }
  }
  return '';
};

const getStandaloneDisplayLabel = (feature, fallbackLabel, context) => {
  const editableEntry = getEditableLabelEntryForStandaloneFeature(feature, context);
  const editableText = String(editableEntry?.text || '').trim();
  if (editableText) return editableText;

  const directOverride = getLabelTextFeatureOverride(feature, context);
  if (directOverride) return directOverride;

  const sourceCandidates = [
    editableEntry?.sourceText,
    fallbackLabel,
    feature?.product,
    feature?.gene,
    feature?.locus_tag
  ].map((value) => String(value || '').trim()).filter(Boolean);
  for (const sourceText of sourceCandidates) {
    const bulkText = String(context?.labelTextBulkOverrides?.[sourceText] ?? '').trim();
    if (bulkText) return bulkText;
  }

  return fallbackLabel;
};

const getStandaloneSearchLabels = (feature, fallbackLabel, displayLabel) => {
  const labels = [
    displayLabel,
    fallbackLabel,
    feature?.label,
    feature?.gene,
    feature?.locus_tag,
    feature?.product,
    firstQualifierValue(feature, 'gene'),
    firstQualifierValue(feature, 'locus_tag'),
    firstQualifierValue(feature, 'product'),
    feature?.svg_id
  ];
  const seen = new Set();
  return labels
    .map((value) => String(value || '').trim())
    .filter((value) => {
      const key = value.toLowerCase();
      if (!value || seen.has(key)) return false;
      seen.add(key);
      return true;
    });
};

const buildStandaloneOrthogroupIndexKey = (recordIndex, svgId) => `${Number(recordIndex)}:${String(svgId || '').trim()}`;

const standaloneDirectOrthogroupMetadataMatches = (feature, member) => {
  if (!member || typeof member !== 'object' || Array.isArray(member)) return false;
  const featureRecord = standaloneFeatureRecordIndexStatus(feature);
  const memberRecord = standaloneMemberRecordIndexStatus(member);
  const featureIndex = standaloneFeatureIndexStatus(feature);
  const memberIndex = standaloneFeatureIndexStatus(member);
  const featureCanonical = standaloneCanonicalIdentityStatus(feature);
  const memberCanonical = standaloneCanonicalIdentityStatus(member);
  const featureIdentity = standaloneFeatureIdentityStatus(feature);
  const memberIdentity = standaloneFeatureIdentityStatus(member);
  if (
    !featureRecord.valid
    || !memberRecord.valid
    || featureRecord.value === null
    || featureRecord.value !== memberRecord.value
    || !featureIndex.valid
    || !memberIndex.valid
    || !featureCanonical.valid
    || !memberCanonical.valid
    || !featureIdentity.valid
    || !memberIdentity.valid
  ) return false;
  if (
    (featureIndex.value !== null || memberIndex.value !== null)
    && (
      featureIndex.value === null
      || memberIndex.value === null
      || featureIndex.value !== memberIndex.value
    )
  ) return false;
  if (
    featureCanonical.value
    || memberCanonical.value
  ) {
    if (
      !featureCanonical.value
      || !memberCanonical.value
      || featureCanonical.value[0] !== memberCanonical.value[0]
      || featureCanonical.value[1] !== memberCanonical.value[1]
    ) return false;
  }
  const featureStable = standaloneFeatureStableSvgId(feature);
  const memberStable = standaloneMemberStableSvgId(member);
  if (!featureStable || !memberStable || featureStable !== memberStable) return false;
  const featureRendered = featureIdentity.rendered;
  const memberRendered = memberIdentity.rendered;
  return (
    (!featureRendered.value && !memberRendered.value)
    || (
      Boolean(featureRendered.value)
      && Boolean(memberRendered.value)
      && featureRendered.value === memberRendered.value
    )
  );
};

const getStandaloneFeatureOrthogroupEntry = (feature, context) => {
  const directIdentity = standaloneConsistentTextAlias(
    feature,
    ['orthogroupId', 'orthogroup_id']
  );
  if (!directIdentity.valid) return null;
  const directOrthogroupId = directIdentity.value;
  if (directOrthogroupId) {
    const matchingGroups = (Array.isArray(context?.orthogroups) ? context.orthogroups : [])
      .filter((group) => {
        const identity = standaloneConsistentTextAlias(
          group,
          ['id', 'orthogroupId', 'orthogroup_id']
        );
        return identity.valid && identity.value === directOrthogroupId;
      });
    if (matchingGroups.length !== 1) return null;
    const members = Array.isArray(matchingGroups[0]?.members) ? matchingGroups[0].members : [];
    const matchingMembers = members.filter((member) => (
      standaloneDirectOrthogroupMetadataMatches(feature, member)
    ));
    if (matchingMembers.length !== 1) return null;
    const member = matchingMembers[0];
    return {
      orthogroupId: directOrthogroupId,
      orthogroupMemberCount: Number(feature?.orthogroupMemberCount || feature?.orthogroup_member_count || 0),
      orthogroupRecordCoverage: Number(feature?.orthogroupRecordCoverage || feature?.orthogroup_record_coverage || 0),
      proteinId: String(feature?.proteinId || feature?.protein_id || '').trim(),
      sourceProteinId: String(feature?.sourceProteinId || feature?.source_protein_id || '').trim(),
      orthogroupRepresentative: Boolean(feature?.orthogroupRepresentative || feature?.orthogroup_representative),
      orthogroupMember: member
    };
  }

  const index = context?.featureOrthogroupIndex instanceof Map ? context.featureOrthogroupIndex : new Map();
  const record = standaloneFeatureRecordIndexStatus(feature);
  const identities = standaloneFeatureIdentityStatus(feature);
  const featureIndex = standaloneFeatureIndexStatus(feature);
  const canonical = standaloneCanonicalIdentityStatus(feature);
  if (
    !record.valid
    || !identities.valid
    || !featureIndex.valid
    || !canonical.valid
    || (
      (identities.source.value || identities.rendered.value || featureIndex.value !== null)
      && record.value === null
    )
    || (
      !identities.source.value
      && !identities.rendered.value
      && featureIndex.value === null
      && !canonical.value
    )
    || index.size === 0
  ) return null;

  const matches = new Set();
  if (identities.source.value) {
    const candidates = new Set([
      index.get(buildStandaloneOrthogroupIndexKey(record.value, identities.source.value)),
      index.get(`record-id:${record.value}:${identities.source.value}`)
    ].filter(Boolean));
    if (candidates.size !== 1) return null;
    matches.add(candidates.values().next().value);
  }
  if (
    identities.rendered.value
    && (
      identities.renderedExplicit
      || index.recordsWithRenderedIds?.has(record.value)
    )
  ) {
    const entry = index.get(
      `record-rendered-id:${record.value}:${identities.rendered.value}`
    );
    if (!entry) return null;
    matches.add(entry);
  }
  if (featureIndex.value !== null) {
    const entry = index.get(
      `record-feature-index:${record.value}:${featureIndex.value}`
    );
    if (!entry) return null;
    matches.add(entry);
  }
  if (canonical.value) {
    const entry = index.get(`canonical:${canonical.value[0]}\u0000${canonical.value[1]}`);
    if (!entry) return null;
    matches.add(entry);
  }
  return matches.size === 1 ? matches.values().next().value : null;
};

const normalizeStandaloneOrthogroupMember = (member, renderedFeatureSvgId = '') => {
  if (!member || typeof member !== 'object' || Array.isArray(member)) return null;
  const start = Number(member.start);
  const end = Number(member.end);
  const record = standaloneMemberRecordIndexStatus(member);
  const recordIndex = record.valid ? record.value : null;
  const stableIdentity = standaloneConsistentTextAlias(member, [
    'stableFeatureSvgId', 'stable_feature_svg_id',
    'stableFeatureId', 'stable_feature_id'
  ]);
  const legacyIdentity = standaloneConsistentTextAlias(
    member,
    ['featureSvgId', 'feature_svg_id']
  );
  const renderedIdentity = standaloneConsistentTextAlias(member, [
    'renderedFeatureSvgId', 'rendered_feature_svg_id',
    'renderedSvgId', 'rendered_svg_id'
  ]);
  const stableFeatureSvgId = standaloneMemberStableSvgId(member);
  const featureIndex = standaloneFeatureIndexStatus(member);
  const canonical = standaloneCanonicalIdentityStatus(member);
  if (
    !record.valid
    || !stableIdentity.valid
    || !legacyIdentity.valid
    || (
      stableIdentity.value
      && legacyIdentity.value
      && stableIdentity.value !== legacyIdentity.value
    )
    || !renderedIdentity.valid
    || !featureIndex.valid
    || !canonical.valid
    || (
      (stableFeatureSvgId || featureIndex.value !== null)
      && recordIndex === null
    )
    || (!stableFeatureSvgId && featureIndex.value === null && !canonical.value)
  ) return null;
  return {
    orthogroup_id: String(member.orthogroupId || member.orthogroup_id || ''),
    protein_id: String(member.proteinId || member.protein_id || ''),
    source_protein_id: String(member.sourceProteinId || member.source_protein_id || ''),
    record_index: recordIndex,
    record_id: String(member.recordId || member.record_id || ''),
    feature_index: featureIndex.valid ? featureIndex.value : null,
    record_key: canonical.value?.[0] || '',
    recordKey: canonical.value?.[0] || '',
    biological_feature_id: canonical.value?.[1] || '',
    biologicalFeatureId: canonical.value?.[1] || '',
    label: String(member.label || ''),
    feature_svg_id: stableFeatureSvgId,
    stable_feature_svg_id: stableFeatureSvgId,
    rendered_feature_svg_id: String(renderedFeatureSvgId || ''),
    start: Number.isFinite(start) ? start : null,
    end: Number.isFinite(end) ? end : null,
    strand: String(member.strand || ''),
    representative: Boolean(member.representative),
    gene: String(member.gene || ''),
    locus_tag: String(member.locusTag || member.locus_tag || ''),
    gene_id: String(member.geneId || member.gene_id || ''),
    old_locus_tag: String(member.oldLocusTag || member.old_locus_tag || ''),
    product: String(member.product || ''),
    note: String(member.note || '')
  };
};

const normalizeStandaloneOrthogroupCandidate = (candidate) => {
  if (!candidate || typeof candidate !== 'object' || Array.isArray(candidate)) return null;
  return {
    text: String(candidate.text || ''),
    source: String(candidate.source || ''),
    member_count: Number.isFinite(Number(candidate.memberCount ?? candidate.member_count)) ? Number(candidate.memberCount ?? candidate.member_count) : 0,
    record_coverage_count: Number.isFinite(Number(candidate.recordCoverageCount ?? candidate.record_coverage_count)) ? Number(candidate.recordCoverageCount ?? candidate.record_coverage_count) : 0,
    representative_count: Number.isFinite(Number(candidate.representativeCount ?? candidate.representative_count)) ? Number(candidate.representativeCount ?? candidate.representative_count) : 0,
    score: Number.isFinite(Number(candidate.score)) ? Number(candidate.score) : 0
  };
};

const getStandaloneOrthogroupDisplayName = (group, context) => {
  const id = String(group?.id || '').trim();
  const override = id ? String(context?.orthogroupNameOverrides?.[id] || '').trim() : '';
  return override || String(group?.display_name || group?.displayName || group?.name || id).trim();
};

const getStandaloneOrthogroupDescription = (group, context) => {
  const id = String(group?.id || '').trim();
  const override = id ? String(context?.orthogroupDescriptionOverrides?.[id] || '').trim() : '';
  return override || String(group?.description || '').trim();
};

const catalogItemWithStandaloneOverrides = (sourceItem, context) => {
  const item = cloneStandaloneCatalogItem(sourceItem);
  const biologicalByKey = new Map();
  for (const feature of (Array.isArray(item.biologicalFeatures) ? item.biologicalFeatures : [])) {
    const recordKey = String(feature?.recordKey || '').trim();
    const biologicalFeatureId = String(feature?.biologicalFeatureId || '').trim();
    if (!recordKey || !biologicalFeatureId) return null;
    const key = standaloneCatalogFeatureKey(recordKey, biologicalFeatureId);
    if (biologicalByKey.has(key)) return null;
    biologicalByKey.set(key, feature);
  }
  const renderedKeys = new Set();
  const renderedIds = new Set();
  item.features = (Array.isArray(item.features) ? item.features : []).map((reference) => {
    const recordKey = String(reference?.recordKey || '').trim();
    const biologicalFeatureId = String(reference?.biologicalFeatureId || '').trim();
    if (!recordKey || !biologicalFeatureId) return null;
    const key = standaloneCatalogFeatureKey(recordKey, biologicalFeatureId);
    const renderedId = String(reference?.svgId || '').trim();
    if (!key || !renderedId || renderedKeys.has(key) || renderedIds.has(renderedId)) return null;
    const biological = biologicalByKey.get(key);
    if (!biological) return null;
    renderedKeys.add(key);
    renderedIds.add(renderedId);
    const feature = {
      ...biological,
      id: key,
      svg_id: renderedId
    };
    const fallbackLabel = getStandaloneFeatureLabel(feature);
    const displayLabel = getStandaloneDisplayLabel(feature, fallbackLabel, context);
    return displayLabel && displayLabel !== fallbackLabel
      ? { ...reference, displayLabel }
      : reference;
  });
  if (item.features.some((reference) => !reference)) return null;
  const groupIds = new Set();
  item.orthogroups = (Array.isArray(item.orthogroups) ? item.orthogroups : [])
    .map((group) => {
      const identity = standaloneConsistentTextAlias(
        group,
        ['id', 'orthogroupId', 'orthogroup_id']
      );
      const id = identity.valid ? identity.value : '';
      if (!id || groupIds.has(id)) return null;
      groupIds.add(id);
      const displayName = id
        ? String(context?.orthogroupNameOverrides?.[id] || '').trim()
        : '';
      const description = id
        ? String(context?.orthogroupDescriptionOverrides?.[id] || '').trim()
        : '';
      return {
        ...group,
        ...(displayName ? { display_name: displayName } : {}),
        ...(description ? { description } : {})
      };
    });
  if (item.orthogroups.some((group) => !group)) return null;
  return item;
};

const createStandaloneBiologicalFeatureResolver = (features) => {
  const byRecordAndId = new Map();
  const ambiguousRecordIds = new Set();
  const addUnique = (map, ambiguous, key, feature) => {
    if (!key || ambiguous.has(key)) return;
    const existing = map.get(key);
    if (existing && existing !== feature) {
      map.delete(key);
      ambiguous.add(key);
      return;
    }
    map.set(key, feature);
  };

  (Array.isArray(features) ? features : []).forEach((feature) => {
    const record = standaloneFeatureRecordIndexStatus(feature);
    const identities = standaloneFeatureIdentityStatus(feature);
    const featureIndex = standaloneFeatureIndexStatus(feature);
    const canonical = standaloneCanonicalIdentityStatus(feature);
    if (
      !record.valid
      || !identities.valid
      || !featureIndex.valid
      || !canonical.valid
      || (
        (identities.source.value || identities.rendered.value || featureIndex.value !== null)
        && record.value === null
      )
      || (
        !identities.source.value
        && featureIndex.value === null
        && !canonical.value
      )
    ) return;
    if (record.value !== null) {
      if (identities.source.value) {
        addUnique(
          byRecordAndId,
          ambiguousRecordIds,
          `record-id:${record.value}:${identities.source.value}`,
          feature
        );
      }
      if (identities.rendered.value) {
        addUnique(
          byRecordAndId,
          ambiguousRecordIds,
          `record-rendered-id:${record.value}:${identities.rendered.value}`,
          feature
        );
      }
      if (featureIndex.value !== null) {
        addUnique(
          byRecordAndId,
          ambiguousRecordIds,
          `record-feature-index:${record.value}:${featureIndex.value}`,
          feature
        );
      }
    }
    if (canonical.value) {
      addUnique(
        byRecordAndId,
        ambiguousRecordIds,
        `canonical:${canonical.value[0]}\u0000${canonical.value[1]}`,
        feature
      );
    }
  });

  return ({
    stableSvgId = '',
    renderedSvgId = '',
    recordIndex = null,
    featureIndex = null,
    recordKey = '',
    biologicalFeatureId = ''
  } = {}) => {
    const record = standaloneIntegerAlias([recordIndex]);
    const sourceIndex = standaloneIntegerAlias([featureIndex]);
    const stableId = String(stableSvgId || '').trim();
    const renderedId = String(renderedSvgId || '').trim();
    const canonical = [
      String(recordKey || '').trim(),
      String(biologicalFeatureId || '').trim()
    ];
    if (
      !record.valid
      || !sourceIndex.valid
      || (Boolean(canonical[0]) !== Boolean(canonical[1]))
      || ((stableId || renderedId || sourceIndex.value !== null) && record.value === null)
      || (!stableId && !renderedId && sourceIndex.value === null && !canonical[0])
    ) return null;
    const matches = new Set();
    if (stableId) {
      const feature = byRecordAndId.get(
        `record-id:${record.value}:${stableId}`
      );
      if (!feature) return null;
      matches.add(feature);
    }
    if (renderedId) {
      const feature = byRecordAndId.get(
        `record-rendered-id:${record.value}:${renderedId}`
      );
      if (!feature) return null;
      matches.add(feature);
    }
    if (sourceIndex.value !== null) {
      const feature = byRecordAndId.get(
        `record-feature-index:${record.value}:${sourceIndex.value}`
      );
      if (!feature) return null;
      matches.add(feature);
    }
    if (canonical[0]) {
      const feature = byRecordAndId.get(
        `canonical:${canonical[0]}\u0000${canonical[1]}`
      );
      if (!feature) return null;
      matches.add(feature);
    }
    return matches.size === 1 ? matches.values().next().value : null;
  };
};

const buildStandaloneOrthogroupPayloads = (features, context) => {
  const resolveRenderedFeature = createStandaloneBiologicalFeatureResolver(features);
  const resolveMemberRenderedFeatureId = (member) => {
    const record = standaloneMemberRecordIndexStatus(member);
    const featureIndex = standaloneFeatureIndexStatus(member);
    const canonical = standaloneCanonicalIdentityStatus(member);
    const stable = standaloneConsistentTextAlias(member, [
      'stableFeatureSvgId', 'stable_feature_svg_id',
      'stableFeatureId', 'stable_feature_id'
    ]);
    const legacy = standaloneConsistentTextAlias(member, ['featureSvgId', 'feature_svg_id']);
    const rendered = standaloneConsistentTextAlias(member, [
      'renderedFeatureSvgId', 'rendered_feature_svg_id',
      'renderedSvgId', 'rendered_svg_id'
    ]);
    if (
      !record.valid
      || !featureIndex.valid
      || !canonical.valid
      || !stable.valid
      || !legacy.valid
      || (stable.value && legacy.value && stable.value !== legacy.value)
      || !rendered.valid
    ) return '';
    const renderedFeature = resolveRenderedFeature({
      stableSvgId: stable.value || legacy.value,
      renderedSvgId: rendered.value,
      recordIndex: record.value,
      featureIndex: featureIndex.value,
      recordKey: canonical.value?.[0] || '',
      biologicalFeatureId: canonical.value?.[1] || ''
    });
    return String(renderedFeature?.svg_id || '').trim();
  };

  const groups = Array.isArray(context?.orthogroups) ? context.orthogroups : [];
  const groupIdCounts = new Map();
  groups.forEach((group) => {
    const identity = standaloneConsistentTextAlias(
      group,
      ['id', 'orthogroupId', 'orthogroup_id']
    );
    if (identity.valid && identity.value) {
      groupIdCounts.set(identity.value, (groupIdCounts.get(identity.value) || 0) + 1);
    }
  });
  return groups
    .map((group) => {
      const groupIdentity = standaloneConsistentTextAlias(
        group,
        ['id', 'orthogroupId', 'orthogroup_id']
      );
      if (
        !groupIdentity.valid
        || !groupIdentity.value
        || groupIdCounts.get(groupIdentity.value) !== 1
      ) return null;
      const members = Array.isArray(group?.members)
        ? group.members
          .map((member) => normalizeStandaloneOrthogroupMember(
            member,
            resolveMemberRenderedFeatureId(member)
          ))
          .filter(Boolean)
        : [];
      const memberCount = Number(group?.member_count || group?.memberCount || members.length || 0);
      const recordCoverageFallback = new Set(
        members
          .map((member) => Number(member.record_index))
          .filter((recordIndex) => Number.isInteger(recordIndex))
      ).size;
      const recordCoverageCount = Number(group?.record_coverage_count || group?.recordCoverageCount || recordCoverageFallback || 0);
      return {
        id: groupIdentity.value,
        name: String(group?.name || ''),
        display_name: getStandaloneOrthogroupDisplayName(group, context),
        description: getStandaloneOrthogroupDescription(group, context),
        name_confidence: String(group?.nameConfidence || group?.name_confidence || ''),
        name_candidates: Array.isArray(group?.nameCandidates || group?.name_candidates)
          ? (group.nameCandidates || group.name_candidates).map(normalizeStandaloneOrthogroupCandidate).filter(Boolean)
          : [],
        member_count: Number.isFinite(memberCount) ? memberCount : members.length,
        record_coverage_count: Number.isFinite(recordCoverageCount) ? recordCoverageCount : 0,
        members
      };
    })
    .filter(Boolean);
};

const collectRenderedFeatureIds = (svg) => {
  const ids = new Set();
  if (!svg) return ids;
  svg.querySelectorAll(FEATURE_SELECTOR).forEach((element) => {
    const id = getElementFeatureId(element);
    if (id) ids.add(id);
  });
  return ids;
};

const collectRenderedFeatureEntries = (svg) => {
  const entries = new Map();
  const ambiguousIds = new Set();
  if (!svg) return entries;
  svg.querySelectorAll(FEATURE_SELECTOR).forEach((element) => {
    const id = getElementFeatureId(element);
    if (!id || ambiguousIds.has(id)) return;
    const recordIndexText = element.getAttribute('data-gbdraw-record-index');
    const record = standaloneIntegerAlias([recordIndexText]);
    const stableId = String(
      element.getAttribute('data-gbdraw-stable-feature-id') || id
    ).trim();
    if (!record.valid) {
      entries.delete(id);
      ambiguousIds.add(id);
      return;
    }
    if (!entries.has(id)) {
      entries.set(id, {
        id,
        stableId,
        recordIndex: record.value,
        recordId: String(element.getAttribute('data-gbdraw-record-id') || '').trim(),
        elements: []
      });
    } else {
      const existing = entries.get(id);
      if (
        existing.stableId !== stableId
        || existing.recordIndex !== record.value
      ) {
        entries.delete(id);
        ambiguousIds.add(id);
        return;
      }
    }
    entries.get(id).elements.push(element);
  });
  return entries;
};

const normalizeColorKey = (value) => String(value || '').trim().toLowerCase();

const getRenderedElementFill = (element) => {
  const fill = String(element?.getAttribute?.('fill') || element?.style?.fill || '').trim();
  if (fill && fill.toLowerCase() !== 'none') return fill;
  return '';
};

const getRenderedFeatureFill = (entry) => {
  const elements = Array.isArray(entry?.elements) ? entry.elements : [];
  for (const element of elements) {
    const fill = getRenderedElementFill(element);
    if (fill) return fill;
  }
  return '';
};

const collectLegendCaptionsByColor = (svg, context) => {
  const captionsByColor = new Map();
  const addCaption = (color, caption) => {
    const key = normalizeColorKey(color);
    const text = String(caption || '').trim();
    if (key && text && !captionsByColor.has(key)) {
      captionsByColor.set(key, text);
    }
  };

  const legendEntries = Array.isArray(context?.legendEntries) ? context.legendEntries : [];
  legendEntries.forEach((entry) => addCaption(entry?.color, entry?.caption || entry?.originalCaption));

  const paletteColors = context?.currentColors && typeof context.currentColors === 'object'
    ? context.currentColors
    : {};
  Object.entries(paletteColors).forEach(([caption, color]) => addCaption(color, caption));

  if (!svg) return captionsByColor;
  const legendRoots = Array.from(svg.querySelectorAll('#legend, g[id*="legend"], [data-gbdraw-sticky-legend]'));
  legendRoots.forEach((root) => {
    root.querySelectorAll('rect, path, polygon').forEach((shape) => {
      const color = getRenderedElementFill(shape);
      if (!color) return;
      const group = shape.closest?.('g') || root;
      const groupText = group?.querySelector?.('text') || null;
      const siblingText = shape.nextElementSibling?.matches?.('text') ? shape.nextElementSibling : null;
      const textElement = groupText || siblingText;
      const caption = String(textElement?.textContent || '').trim();
      addCaption(color, caption);
    });
  });
  return captionsByColor;
};

const buildFallbackStandaloneFeaturePayload = (svgId, entry, captionsByColor) => {
  const fillColor = getRenderedFeatureFill(entry);
  const caption = captionsByColor.get(normalizeColorKey(fillColor)) || 'Feature';
  const label = caption === 'Feature' ? String(svgId) : caption;
  const searchLabels = Array.from(new Set([label, caption, svgId].map((value) => String(value || '').trim()).filter(Boolean)));
  return {
    svg_id: String(svgId || ''),
    stable_svg_id: String(entry?.stableId || svgId || ''),
    label,
    display_label: label,
    search_labels: searchLabels,
    record_id: String(entry?.recordId || ''),
    record_idx: Number.isInteger(entry?.recordIndex) ? entry.recordIndex : null,
    type: caption,
    start: null,
    end: null,
    strand: '',
    location: '',
    locus_tag: '',
    gene_id: '',
    old_locus_tag: '',
    fill_color: fillColor,
    orthogroup_id: '',
    orthogroup_member_count: 0,
    orthogroup_record_coverage: 0,
    protein_id: '',
    source_protein_id: '',
    orthogroup_representative: false,
    qualifiers: {},
    location_parts: [],
    nucleotide_sequence: '',
    amino_acid_sequence: '',
    sequence_warnings: [],
    nucleotide_fasta: '',
    amino_acid_fasta: '',
    orthogroup_member: null
  };
};

const normalizeStandalonePopupMode = (popupMode) => (
  popupMode === 'simple' ? 'simple' : 'rich'
);

const normalizeStandaloneBiologicalFeature = (feature, context) => {
  if (!feature || typeof feature !== 'object' || Array.isArray(feature)) return null;
  const record = standaloneFeatureRecordIndexStatus(feature);
  const featureIndex = standaloneFeatureIndexStatus(feature);
  const canonical = standaloneCanonicalIdentityStatus(feature);
  const identity = standaloneFeatureIdentityStatus(feature);
  if (
    !record.valid
    || !featureIndex.valid
    || !canonical.valid
    || !identity.valid
    || ((identity.source.value || featureIndex.value !== null) && record.value === null)
  ) return null;
  const stableSvgId = standaloneFeatureStableSvgId(feature);
  if (!stableSvgId) return null;
  const fallbackLabel = getStandaloneFeatureLabel(feature);
  const displayLabel = getStandaloneDisplayLabel(feature, fallbackLabel, context);
  const qualifiers = normalizeQualifierMap(feature.qualifiers);
  const translation = normalizeStringArray(qualifiers.translation).find((value) => value.trim()) || '';
  const orthogroupEntry = getStandaloneFeatureOrthogroupEntry(feature, context);
  return compactWireValue({
    svg_id: stableSvgId,
    stable_svg_id: stableSvgId,
    stable_feature_id: stableSvgId,
    record_idx: record.value,
    record_id: String(feature.record_id || feature.recordId || ''),
    feature_index: featureIndex.value,
    record_key: canonical.value?.[0] || '',
    biological_feature_id: canonical.value?.[1] || '',
    label: fallbackLabel,
    display_label: displayLabel,
    search_labels: getStandaloneSearchLabels(feature, fallbackLabel, displayLabel),
    type: String(feature.type || ''),
    start: Number.isFinite(Number(feature.start)) ? Number(feature.start) : null,
    end: Number.isFinite(Number(feature.end)) ? Number(feature.end) : null,
    strand: String(feature.strand || ''),
    location: buildStandaloneFeatureLocation(feature),
    locus_tag: String(feature.locus_tag || feature.locusTag || ''),
    gene_id: String(feature.gene_id || feature.geneId || ''),
    old_locus_tag: String(feature.old_locus_tag || feature.oldLocusTag || ''),
    gene: String(feature.gene || firstQualifierValue(feature, 'gene') || ''),
    product: String(feature.product || firstQualifierValue(feature, 'product') || ''),
    note: String(feature.note || firstQualifierValue(feature, 'note') || ''),
    orthogroup_id: String(orthogroupEntry?.orthogroupId || ''),
    orthogroup_member_count: Number.isFinite(Number(orthogroupEntry?.orthogroupMemberCount))
      ? Number(orthogroupEntry.orthogroupMemberCount)
      : 0,
    orthogroup_record_coverage: Number.isFinite(Number(orthogroupEntry?.orthogroupRecordCoverage))
      ? Number(orthogroupEntry.orthogroupRecordCoverage)
      : 0,
    protein_id: String(orthogroupEntry?.proteinId || feature.proteinId || feature.protein_id || ''),
    source_protein_id: String(
      orthogroupEntry?.sourceProteinId ||
      feature.sourceProteinId ||
      feature.source_protein_id ||
      firstQualifierValue(feature, 'protein_id') ||
      ''
    ),
    orthogroup_representative: Boolean(orthogroupEntry?.orthogroupRepresentative),
    qualifiers,
    location_parts: normalizeLocationParts(feature.location_parts),
    nucleotide_sequence: String(feature.nucleotide_sequence || feature.nucleotideSequence || ''),
    amino_acid_sequence: String(
      feature.amino_acid_sequence || feature.aminoAcidSequence || translation || ''
    ),
    sequence_warnings: normalizeStringArray(feature.sequence_warnings || feature.sequenceWarnings),
    orthogroup_member: normalizeStandaloneOrthogroupMember(
      orthogroupEntry?.orthogroupMember
    )
  });
};

const buildStandaloneBiologicalFeaturePayloads = (context, renderedFeatures) => {
  const candidates = [];
  const keyCounts = new Map();
  const resolveRenderedFeature = createStandaloneBiologicalFeatureResolver(renderedFeatures);
  (Array.isArray(context?.features) ? context.features : []).forEach((feature) => {
    const payload = normalizeStandaloneBiologicalFeature(feature, context);
    if (!payload) return;
    const recordIndex = Number(payload.record_idx);
    const canonicalKey = payload.record_key && payload.biological_feature_id
      ? standaloneCatalogFeatureKey(payload.record_key, payload.biological_feature_id)
      : '';
    const key = canonicalKey || (
      Number.isInteger(recordIndex)
        ? buildStandaloneOrthogroupIndexKey(recordIndex, payload.stable_svg_id)
        : `?:${payload.stable_svg_id}`
    );
    keyCounts.set(key, (keyCounts.get(key) || 0) + 1);
    candidates.push({ feature, payload, recordIndex, key });
  });
  const payloads = [];
  candidates.forEach(({ feature, payload, recordIndex, key }) => {
    if (keyCounts.get(key) !== 1) return;
    const renderedFeature = resolveRenderedFeature({
      stableSvgId: payload.stable_svg_id,
      renderedSvgId: String(feature?.svg_id || '').trim(),
      recordIndex: Number.isInteger(recordIndex) ? recordIndex : null
    });
    if (renderedFeature?.svg_id) {
      payload.rendered_svg_id = renderedFeature.svg_id;
    }
    payloads.push(payload);
  });
  return payloads;
};

const buildStandaloneFeaturePayloads = (svg, options = {}) => {
  const renderedEntries = collectRenderedFeatureEntries(svg);
  if (renderedEntries.size === 0) return [];

  const context = normalizeStandaloneContext(options);
  const normalizedPopupMode = normalizeStandalonePopupMode(context.popupMode);
  const resolveSourceFeature = createStandaloneBiologicalFeatureResolver(context.features);
  const payloads = [];
  const captionsByColor = collectLegendCaptionsByColor(svg, context);
  renderedEntries.forEach((entry, svgId) => {
    const feature = resolveSourceFeature({
      stableSvgId: entry.stableId,
      renderedSvgId: svgId,
      recordIndex: entry.recordIndex
    });
    if (!feature) {
      payloads.push(compactWireValue(buildFallbackStandaloneFeaturePayload(svgId, entry, captionsByColor)));
      return;
    }
    const fallbackLabel = getStandaloneFeatureLabel(feature);
    const displayLabel = getStandaloneDisplayLabel(feature, fallbackLabel, context);
    const orthogroupEntry = getStandaloneFeatureOrthogroupEntry(feature, context);
    const orthogroupMember = normalizeStandaloneOrthogroupMember(
      orthogroupEntry?.orthogroupMember,
      svgId
    );
    const qualifiers = normalizeQualifierMap(feature?.qualifiers);
    const aminoAcidSequence = String(feature?.amino_acid_sequence || '');
    const translation = normalizeStringArray(qualifiers.translation).find((value) => value.trim()) || '';
    const stableSvgId = String(entry.stableId || standaloneFeatureStableSvgId(feature) || svgId).trim();
    const recordIndex = Number.isInteger(entry.recordIndex)
      ? entry.recordIndex
      : standaloneFeatureRecordIndex(feature);
    const payload = {
      svg_id: svgId,
      stable_svg_id: stableSvgId,
      label: fallbackLabel,
      display_label: displayLabel,
      search_labels: Array.from(new Set([
        ...getStandaloneSearchLabels(feature, fallbackLabel, displayLabel),
        stableSvgId,
        svgId
      ].filter(Boolean))),
      record_id: String(entry.recordId || feature?.record_id || feature?.recordId || ''),
      record_idx: Number.isInteger(recordIndex) ? recordIndex : null,
      type: String(feature?.type || ''),
      start: Number.isFinite(Number(feature?.start)) ? Number(feature.start) : null,
      end: Number.isFinite(Number(feature?.end)) ? Number(feature.end) : null,
      strand: String(feature?.strand || ''),
      location: buildStandaloneFeatureLocation(feature),
      locus_tag: String(feature?.locus_tag || feature?.locusTag || ''),
      gene_id: String(feature?.gene_id || feature?.geneId || ''),
      old_locus_tag: String(feature?.old_locus_tag || feature?.oldLocusTag || ''),
      orthogroup_id: String(orthogroupEntry?.orthogroupId || ''),
      orthogroup_member_count: Number.isFinite(Number(orthogroupEntry?.orthogroupMemberCount))
        ? Number(orthogroupEntry.orthogroupMemberCount)
        : 0,
      orthogroup_record_coverage: Number.isFinite(Number(orthogroupEntry?.orthogroupRecordCoverage))
        ? Number(orthogroupEntry.orthogroupRecordCoverage)
        : 0,
      protein_id: String(orthogroupEntry?.proteinId || feature?.proteinId || feature?.protein_id || ''),
      source_protein_id: String(
        orthogroupEntry?.sourceProteinId ||
        feature?.sourceProteinId ||
        feature?.source_protein_id ||
        firstQualifierValue(feature, 'protein_id') ||
        ''
      ),
      orthogroup_representative: Boolean(orthogroupEntry?.orthogroupRepresentative)
    };
    if (normalizedPopupMode === 'rich') {
      Object.assign(payload, {
        qualifiers,
        location_parts: normalizeLocationParts(feature?.location_parts),
        nucleotide_sequence: String(feature?.nucleotide_sequence || ''),
        amino_acid_sequence: aminoAcidSequence && aminoAcidSequence !== translation
          ? aminoAcidSequence
          : '',
        sequence_warnings: normalizeStringArray(feature?.sequence_warnings),
        orthogroup_member: orthogroupMember
      });
    }
    payloads.push(compactWireValue(payload));
  });
  return payloads;
};

const standaloneAttr = (element, name) =>
  String(element?.getAttribute?.(name) || '').trim();

const standaloneElementTextAlias = (element, names) => {
  const values = new Set(names
    .map((name) => standaloneAttr(element, name))
    .filter(Boolean));
  return {
    valid: values.size <= 1,
    supplied: values.size > 0,
    value: values.size === 1 ? values.values().next().value : ''
  };
};

const standaloneElementMatchId = (element) => standaloneElementTextAlias(element, [
  'data-gbdraw-match-id',
  'data-gbdraw-pairwise-match-id'
]);

const firstStandaloneText = (...values) => {
  for (const value of values) {
    const text = String(value === null || value === undefined ? '' : value).trim();
    if (text) return text;
  }
  return '';
};

const splitStandaloneMetadataValues = (value) => firstStandaloneText(value)
  .split(';')
  .map((entry) => firstStandaloneText(entry))
  .filter(Boolean);

const uniqueStandaloneMetadataValues = (value) => {
  const seen = new Set();
  const values = [];
  splitStandaloneMetadataValues(value).forEach((entry) => {
    if (seen.has(entry)) return;
    seen.add(entry);
    values.push(entry);
  });
  return values;
};

const standaloneGeneratedProteinIdPattern =
  /^(?:h_[a-z2-7]{26}|f_[0-9a-f]{64}|gbd_r\d+_cds\d+|p_.+_\d+_\d+_-?\d+_[0-9a-f]{12}(?:_\d+)?)$/i;
// Keep unsupported historical shapes from leaking through display-only fallbacks.
const standaloneUnsupportedHistoricalProteinIdPattern = /@[^|]+\|.+~f_[0-9a-f]{64}$/i;
const standaloneGeneratedUnitIdPattern = /^gbd_r\d+_unit\d+$/i;

const isStandaloneInternalDisplayId = (value) => {
  const text = firstStandaloneText(value);
  return Boolean(text && (
    standaloneGeneratedProteinIdPattern.test(text) ||
    standaloneUnsupportedHistoricalProteinIdPattern.test(text) ||
    standaloneGeneratedUnitIdPattern.test(text) ||
    /^p_r_/i.test(text)
  ));
};

const addUniqueStandaloneDisplayText = (values, value) => {
  const text = firstStandaloneText(value);
  if (!text || isStandaloneInternalDisplayId(text) || values.includes(text)) return;
  values.push(text);
};

const firstStandaloneNonInternalDisplayText = (...values) => {
  for (const value of values) {
    if (Array.isArray(value)) {
      const nested = firstStandaloneNonInternalDisplayText(...value);
      if (nested) return nested;
      continue;
    }
    const text = firstStandaloneText(value);
    if (text && !isStandaloneInternalDisplayId(text)) return text;
  }
  return '';
};

const getStandaloneQualifierDisplayValue = (feature, key) => {
  const normalizedKey = String(key || '').trim().toLowerCase();
  const qualifiers = feature?.qualifiers && typeof feature.qualifiers === 'object' && !Array.isArray(feature.qualifiers)
    ? feature.qualifiers
    : {};
  if (!normalizedKey) return '';
  if (Object.prototype.hasOwnProperty.call(qualifiers, normalizedKey)) {
    return firstStandaloneNonInternalDisplayText(
      normalizeStringArray(qualifiers[normalizedKey])
    );
  }
  const matchingKey = Object.keys(qualifiers).find(
    (candidate) => candidate.toLowerCase() === normalizedKey
  );
  return matchingKey
    ? firstStandaloneNonInternalDisplayText(normalizeStringArray(qualifiers[matchingKey]))
    : '';
};

const resolveStandaloneDisplayProteinId = (feature, member = null, fallback = '') => {
  return firstStandaloneNonInternalDisplayText(
    feature?.displayProteinId,
    feature?.display_protein_id,
    member?.displayProteinId,
    member?.display_protein_id,
    feature?.sourceProteinId,
    feature?.source_protein_id,
    member?.sourceProteinId,
    member?.source_protein_id,
    getStandaloneQualifierDisplayValue(feature, 'protein_id'),
    feature?.locusTag,
    feature?.locus_tag,
    getStandaloneQualifierDisplayValue(feature, 'locus_tag'),
    member?.locusTag,
    member?.locus_tag,
    feature?.geneId,
    feature?.gene_id,
    getStandaloneQualifierDisplayValue(feature, 'gene_id'),
    member?.geneId,
    member?.gene_id,
    feature?.oldLocusTag,
    feature?.old_locus_tag,
    getStandaloneQualifierDisplayValue(feature, 'old_locus_tag'),
    member?.oldLocusTag,
    member?.old_locus_tag,
    feature?.ID,
    getStandaloneQualifierDisplayValue(feature, 'ID'),
    feature?.Name,
    getStandaloneQualifierDisplayValue(feature, 'Name'),
    feature?.Parent,
    getStandaloneQualifierDisplayValue(feature, 'Parent'),
    feature?.gene,
    getStandaloneQualifierDisplayValue(feature, 'gene'),
    member?.gene,
    member?.label,
    feature?.proteinId,
    feature?.protein_id,
    member?.proteinId,
    member?.protein_id,
    fallback
  );
};

const addStandaloneMatchRow = (rows, label, value) => {
  const text = String(value === null || value === undefined ? '' : value).trim();
  if (!text) return;
  rows.push([label, text]);
};

const standaloneIntervalText = (start, end) => {
  const startText = String(start || '').trim();
  const endText = String(end || '').trim();
  if (startText && endText) return `${startText}..${endText}`;
  return startText || endText;
};

const standaloneStrandText = (strand) => {
  const text = String(strand || '').trim();
  if (text === '1') return '+';
  if (text === '-1') return '-';
  return text;
};

const standaloneMemberLocationText = (member) => {
  if (!member || typeof member !== 'object') return '';
  const start = Number(member.start);
  const end = Number(member.end);
  const startText = Number.isFinite(start) ? String(start + 1) : String(member.start ?? '').trim();
  const endText = Number.isFinite(end) ? String(end) : String(member.end ?? '').trim();
  const range = startText && endText ? `${startText}..${endText}` : startText || endText;
  const strand = standaloneStrandText(member.strand);
  return range && strand ? `${range} (${strand})` : range;
};

const buildStandaloneMatchMemberRows = (orthogroup) => {
  const members = Array.isArray(orthogroup?.members) ? orthogroup.members : [];
  const orthogroupId = firstStandaloneText(orthogroup?.id, orthogroup?.orthogroupId, orthogroup?.orthogroup_id);
  const displayName = firstStandaloneText(orthogroup?.display_name, orthogroup?.displayName, orthogroup?.name);
  return members
    .map((member) => ({
      featureSvgId: standaloneMemberFeatureSvgId(member),
      feature_svg_id: standaloneMemberFeatureSvgId(member),
      stableFeatureSvgId: standaloneMemberFeatureSvgId(member),
      stable_feature_svg_id: standaloneMemberFeatureSvgId(member),
      renderedFeatureSvgId: standaloneMemberRenderedFeatureSvgId(member),
      rendered_feature_svg_id: standaloneMemberRenderedFeatureSvgId(member),
      recordIndex: standaloneMemberRecordIndex(member),
      record_index: standaloneMemberRecordIndex(member),
      orthogroupId,
      orthogroup_id: orthogroupId,
      displayName,
      display_name: displayName,
      record: firstStandaloneText(member?.recordId, member?.record_id),
      coordinates: standaloneMemberLocationText(member),
      proteinId: resolveStandaloneDisplayProteinId(null, member),
      productOrNote: firstStandaloneText(member?.product, member?.note)
    }))
    .filter((row) => row.record || row.coordinates || row.proteinId || row.productOrNote || row.featureSvgId);
};

const standaloneMemberCopyText = (memberRows) => {
  if (!Array.isArray(memberRows) || memberRows.length === 0) return '';
  return [
    'Record\tCoordinates (+/-)\tProtein ID\tProduct / note',
    ...memberRows.map((row) => [row.record, row.coordinates, row.proteinId, row.productOrNote].join('\t'))
  ].join('\n');
};

const standaloneMatchKind = (element) => {
  const explicit = standaloneAttr(element, 'data-match-kind').toLowerCase();
  if (['homology', 'pairwise', 'orthogroup', 'collinear'].includes(explicit)) return explicit;
  if (standaloneAttr(element, 'data-collinearity-block-id')) return 'collinear';
  if (standaloneAttr(element, 'data-orthogroup-id')) return 'orthogroup';
  return 'pairwise';
};

const STANDALONE_MATCH_TITLES = {
  homology: 'Homology ring match',
  pairwise: 'Pairwise match',
  orthogroup: 'Similarity-group match',
  collinear: 'Collinearity block'
};

const standaloneOrthogroupTitle = (orthogroupId, displayName) => {
  const id = String(orthogroupId || '').trim();
  const name = String(displayName || '').trim();
  if (id && name) return `${id}:${name}`;
  return id || name || STANDALONE_MATCH_TITLES.orthogroup;
};

const buildStandaloneMatchSection = (title, rows) => ({
  title,
  rows: rows.filter((row) => String(row?.[1] || '').trim())
});

const standaloneMemberFeatureSvgId = (member) => standaloneMemberStableSvgId(member);

const standaloneMemberRenderedFeatureSvgId = (member) => standaloneMemberRenderedSvgId(member);

const standaloneOrthogroupId = (entry) => {
  const identity = standaloneConsistentTextAlias(
    entry,
    ['id', 'orthogroupId', 'orthogroup_id']
  );
  return identity.valid ? identity.value : '';
};

const getStandaloneOrthogroupById = (orthogroups, orthogroupId) => {
  const id = firstStandaloneText(orthogroupId);
  if (!id) return null;
  const matches = (Array.isArray(orthogroups) ? orthogroups : [])
    .filter((entry) => standaloneOrthogroupId(entry) === id);
  return matches.length === 1 ? matches[0] : null;
};

const standaloneGroupHasFeatureSvgId = (group, featureSvgId, featuresById) => {
  const ids = splitStandaloneMetadataValues(featureSvgId);
  if (ids.length === 0) return false;
  return ids.some((id) => Boolean(
    getStandaloneGroupMemberForFeatureSvgId(group, id, featuresById)
  ));
};

const getStandaloneOrthogroupPayload = (
  orthogroups,
  orthogroupId,
  queryFeatureSvgId = '',
  subjectFeatureSvgId = '',
  featuresById = new Map()
) => {
  const id = String(orthogroupId || '').trim();
  const groups = Array.isArray(orthogroups) ? orthogroups : [];
  const directMatches = id
    ? groups.filter((entry) => standaloneOrthogroupId(entry) === id)
    : [];
  if (directMatches.length > 0) {
    return directMatches.length === 1 ? directMatches[0] : null;
  }
  const matches = groups.filter((entry) => (
    standaloneGroupHasFeatureSvgId(entry, queryFeatureSvgId, featuresById) ||
    standaloneGroupHasFeatureSvgId(entry, subjectFeatureSvgId, featuresById)
  ));
  return matches.length === 1 ? matches[0] : null;
};

const standaloneFeatureOrthogroupId = (feature) => {
  const ids = new Set([
    feature?.orthogroup_id,
    feature?.orthogroupId,
    feature?.orthogroup_member?.orthogroupId,
    feature?.orthogroup_member?.orthogroup_id,
    feature?.orthogroupMember?.orthogroupId,
    feature?.orthogroupMember?.orthogroup_id
  ].map((value) => String(value || '').trim()).filter(Boolean));
  return ids.size === 1 ? ids.values().next().value : '';
};

const standaloneFeatureProduct = (feature) => firstStandaloneNonInternalDisplayText(
  feature?.product,
  feature?.orthogroup_member?.product,
  feature?.orthogroupMember?.product,
  feature?.note,
  feature?.label,
  feature?.display_label,
  feature?.displayLabel
);

const standaloneFeatureToMember = (feature) => {
  if (!feature) return null;
  const member = feature?.orthogroup_member || feature?.orthogroupMember || {};
  return {
    recordId: firstStandaloneText(member?.recordId, member?.record_id, feature?.record_id),
    start: member?.start ?? feature?.start,
    end: member?.end ?? feature?.end,
    strand: firstStandaloneText(member?.strand, feature?.strand),
    proteinId: firstStandaloneText(member?.proteinId, member?.protein_id, feature?.protein_id, feature?.proteinId),
    sourceProteinId: firstStandaloneText(member?.sourceProteinId, member?.source_protein_id, feature?.source_protein_id, feature?.sourceProteinId),
    product: firstStandaloneText(member?.product, standaloneFeatureProduct(feature)),
    note: firstStandaloneText(member?.note, feature?.note),
    featureSvgId: firstStandaloneText(
      member?.stableFeatureSvgId,
      member?.stable_feature_svg_id,
      member?.featureSvgId,
      member?.feature_svg_id,
      feature?.stable_svg_id,
      feature?.svg_id
    ),
    stableFeatureSvgId: firstStandaloneText(
      member?.stableFeatureSvgId,
      member?.stable_feature_svg_id,
      feature?.stable_svg_id,
      member?.featureSvgId,
      member?.feature_svg_id,
      feature?.svg_id
    ),
    renderedFeatureSvgId: firstStandaloneText(
      member?.renderedFeatureSvgId,
      member?.rendered_feature_svg_id,
      feature?.svg_id
    )
  };
};

const buildStandaloneFallbackOrthogroup = ({ orthogroupId, queryFeature, subjectFeature, features }) => {
  const id = String(orthogroupId || '').trim();
  if (!id) return null;
  const matchingFeatures = (Array.isArray(features) ? features : [])
    .filter((feature) => standaloneFeatureOrthogroupId(feature) === id);
  const fallbackFeatures = matchingFeatures.length
    ? matchingFeatures
    : [queryFeature, subjectFeature].filter(Boolean);
  if (!fallbackFeatures.length) return null;
  const members = fallbackFeatures.map(standaloneFeatureToMember).filter(Boolean);
  const firstFeature = fallbackFeatures[0] || {};
  const recordCoverageFallback = new Set(
    fallbackFeatures
      .map((feature) => firstStandaloneText(feature?.record_id, feature?.recordId, feature?.record_idx, feature?.recordIndex))
      .filter(Boolean)
  ).size;
  return {
    id,
    name: standaloneFeatureProduct(firstFeature),
    member_count: firstStandaloneText(firstFeature?.orthogroup_member_count, firstFeature?.orthogroupMemberCount, members.length),
    record_coverage_count: firstStandaloneText(firstFeature?.orthogroup_record_coverage, firstFeature?.orthogroupRecordCoverage, recordCoverageFallback),
    members
  };
};

const getStandaloneGroupMemberForFeatureSvgId = (group, featureSvgId, featuresById) => {
  const id = firstStandaloneText(featureSvgId);
  const feature = id ? featuresById?.get?.(id) || null : null;
  if (!feature) return null;
  const members = Array.isArray(group?.members) ? group.members : [];
  const matches = members.filter((member) => (
    standaloneDirectOrthogroupMetadataMatches(feature, member)
  ));
  return matches.length === 1 ? matches[0] : null;
};

const resolveStandaloneFeatureSectionProteinIds = ({
  feature,
  featureSvgIds,
  featuresById,
  group,
  fallbackProteinIds,
  locusId,
  displayName
}) => {
  const values = [];
  const addFeatureProteinId = (candidateFeature, member = null) => {
    const text = resolveStandaloneDisplayProteinId(candidateFeature, member, '');
    addUniqueStandaloneDisplayText(values, text);
  };

  const ids = splitStandaloneMetadataValues(featureSvgIds);
  ids.forEach((featureSvgId) => {
    const candidateFeature = featuresById?.get?.(featureSvgId) || null;
    const member = getStandaloneGroupMemberForFeatureSvgId(group, featureSvgId, featuresById);
    addFeatureProteinId(candidateFeature, member);
  });
  if (feature) {
    addFeatureProteinId(
      feature,
      ids.length === 1
        ? getStandaloneGroupMemberForFeatureSvgId(group, ids[0], featuresById)
        : null
    );
  }
  if (values.length === 0) {
    splitStandaloneMetadataValues(locusId).forEach((value) => addUniqueStandaloneDisplayText(values, value));
  }
  if (values.length === 0) {
    splitStandaloneMetadataValues(displayName).forEach((value) => addUniqueStandaloneDisplayText(values, value));
  }
  if (values.length === 0) {
    splitStandaloneMetadataValues(fallbackProteinIds).forEach((value) => addUniqueStandaloneDisplayText(values, value));
  }
  return values.join('; ');
};

const standaloneFeatureLocationText = (feature) => {
  const direct = firstStandaloneText(feature?.location);
  if (direct && direct !== '..') return direct;
  const built = buildStandaloneFeatureLocation(feature);
  return built && built !== '..' ? built : '';
};

const standaloneFeatureRowLocusId = (feature, fallbackLocusId) => firstStandaloneNonInternalDisplayText(
  feature?.locusTag,
  feature?.locus_tag,
  getStandaloneQualifierDisplayValue(feature, 'locus_tag'),
  feature?.geneId,
  feature?.gene_id,
  getStandaloneQualifierDisplayValue(feature, 'gene_id'),
  fallbackLocusId
);

const standaloneFeatureRowDisplayName = (feature, fallbackDisplayName) => firstStandaloneNonInternalDisplayText(
  fallbackDisplayName,
  feature?.displayLabel,
  feature?.display_label,
  feature?.label,
  feature?.gene,
  getStandaloneQualifierDisplayValue(feature, 'gene'),
  feature?.locus_tag,
  feature?.locusTag,
  getStandaloneQualifierDisplayValue(feature, 'locus_tag'),
  standaloneFeatureProduct(feature)
);

const buildStandaloneFeatureListRows = ({
  featureSvgIds,
  featuresById,
  group,
  recordId,
  interval,
  proteinId,
  locusId,
  displayName
}) => {
  const featureIds = uniqueStandaloneMetadataValues(featureSvgIds);
  const proteinIds = splitStandaloneMetadataValues(proteinId);
  const locusIds = splitStandaloneMetadataValues(locusId);
  const displayNames = splitStandaloneMetadataValues(displayName);
  const count = Math.max(featureIds.length, proteinIds.length, locusIds.length, displayNames.length);
  if (count === 0) return [];

  return Array.from({ length: count }, (_unused, index) => {
    const svgId = featureIds[index] || '';
    const feature = svgId ? featuresById?.get?.(svgId) || null : null;
    const member = svgId
      ? getStandaloneGroupMemberForFeatureSvgId(group, svgId, featuresById)
      : null;
    const fallbackProteinId = firstStandaloneNonInternalDisplayText(
      locusIds[index],
      displayNames[index],
      proteinIds[index]
    );
    const resolvedProteinId = resolveStandaloneDisplayProteinId(feature, member, '');
    const displayProteinId = firstStandaloneText(
      isStandaloneInternalDisplayId(resolvedProteinId) ? '' : resolvedProteinId,
      fallbackProteinId,
      isStandaloneInternalDisplayId(proteinIds[index]) ? '' : proteinIds[index]
    );
    const rowRecord = firstStandaloneText(feature?.record_id, feature?.recordId, member?.recordId, member?.record_id, recordId);
    const rowLocation = firstStandaloneText(
      standaloneFeatureLocationText(feature),
      standaloneMemberLocationText(member),
      count === 1 ? interval : ''
    );
    const rowLocusId = standaloneFeatureRowLocusId(feature, locusIds[index]);
    const rowDisplayName = standaloneFeatureRowDisplayName(feature, displayNames[index]);
    const product = standaloneFeatureProduct(feature);
    const label = firstStandaloneNonInternalDisplayText(
      displayProteinId,
      rowDisplayName,
      rowLocusId,
      product,
      svgId,
      `Feature ${index + 1}`
    );
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
      svg_id: svgId,
      svgId,
      can_open: Boolean(feature?.svg_id),
      canOpen: Boolean(feature?.svg_id),
      label,
      record: rowRecord,
      location: rowLocation,
      protein_id: displayProteinId,
      proteinId: displayProteinId,
      locus_id: rowLocusId,
      locusId: rowLocusId,
      display_name: rowDisplayName,
      displayName: rowDisplayName,
      product,
      type: firstStandaloneText(feature?.type),
      copy_text: copyText,
      copyText
    };
  }).filter((row) => (
    row.svg_id ||
    row.record ||
    row.location ||
    row.protein_id ||
    row.locus_id ||
    row.display_name ||
    row.product
  ));
};

const resolveStandaloneBlockMemberLabels = ({
  group,
  featureSvgIds,
  featuresById
}) => {
  if (!group) return '';
  const values = [];
  uniqueStandaloneMetadataValues(featureSvgIds).forEach((featureSvgId) => {
    const feature = featuresById?.get?.(featureSvgId) || null;
    const member = group
      ? getStandaloneGroupMemberForFeatureSvgId(group, featureSvgId, featuresById)
      : null;
    if (!member) return;
    addUniqueStandaloneDisplayText(values, resolveStandaloneDisplayProteinId(feature, member, ''));
  });
  return values.join('; ');
};

const buildStandaloneOrthogroupDetailRows = ({
  orthogroupId,
  displayName,
  description,
  memberCount,
  recordCoverage
}) => {
  const rows = [];
  addStandaloneMatchRow(rows, 'Similarity group ID', orthogroupId);
  addStandaloneMatchRow(rows, 'Display name', displayName);
  addStandaloneMatchRow(rows, 'Description', description);
  addStandaloneMatchRow(rows, 'Members', memberCount);
  addStandaloneMatchRow(rows, 'Record coverage', recordCoverage);
  return rows;
};

const buildStandaloneBlockOrthogroups = ({
  orthogroupIds,
  orthogroups,
  featuresById,
  queryFeatureSvgId,
  subjectFeatureSvgId
}) => orthogroupIds.map((orthogroupId) => {
  const group = getStandaloneOrthogroupById(orthogroups, orthogroupId);
  const displayName = firstStandaloneText(group?.display_name, group?.displayName, group?.name);
  const description = firstStandaloneText(group?.description);
  const memberCount = firstStandaloneText(
    group?.member_count,
    group?.memberCount,
    Array.isArray(group?.members) ? group.members.length : ''
  );
  const recordCoverage = firstStandaloneText(group?.record_coverage_count, group?.recordCoverage);
  const memberRows = buildStandaloneMatchMemberRows(group);
  return {
    id: orthogroupId,
    display_name: displayName,
    displayName,
    description,
    member_count: memberCount,
    memberCount,
    record_coverage: recordCoverage,
    recordCoverage,
    query_member: resolveStandaloneBlockMemberLabels({ group, featureSvgIds: queryFeatureSvgId, featuresById }),
    subject_member: resolveStandaloneBlockMemberLabels({ group, featureSvgIds: subjectFeatureSvgId, featuresById }),
    detail_rows: buildStandaloneOrthogroupDetailRows({
      orthogroupId,
      displayName,
      description,
      memberCount,
      recordCoverage
    }),
    member_rows: memberRows,
    member_copy_text: standaloneMemberCopyText(memberRows)
  };
});

const buildStandaloneMatchFeatureSection = ({
  title,
  feature,
  recordId,
  interval,
  proteinId,
  locusId,
  displayName,
  featureSvgIds,
  featuresById,
  group
}) => {
  const rows = [];
  const featureRows = buildStandaloneFeatureListRows({
    featureSvgIds,
    featuresById,
    group,
    recordId,
    interval,
    proteinId,
    locusId,
    displayName
  });
  const displayProteinIds = resolveStandaloneFeatureSectionProteinIds({
    feature,
    featureSvgIds,
    featuresById,
    group,
    fallbackProteinIds: proteinId,
    locusId,
    displayName
  });
  addStandaloneMatchRow(rows, 'Record', firstStandaloneText(feature?.record_id, recordId));
  addStandaloneMatchRow(rows, 'Location', firstStandaloneText(feature?.location, interval));
  addStandaloneMatchRow(rows, displayProteinIds.includes(';') ? 'Protein IDs' : 'Protein ID', displayProteinIds);
  addStandaloneMatchRow(rows, 'Type', feature?.type);
  addStandaloneMatchRow(rows, 'Locus ID', locusId);
  addStandaloneMatchRow(rows, 'Display name', displayName);
  return {
    ...buildStandaloneMatchSection(title, rows),
    feature_rows: featureRows,
    featureRows
  };
};

const buildStandaloneMatchPayloadsV1 = (svg, { features = [], orthogroups = [] } = {}) => {
  if (!svg) return [];
  const featuresById = new Map();
  const ambiguousFeatureIds = new Set();
  (Array.isArray(features) ? features : []).forEach((feature) => {
    const svgId = String(feature?.svg_id || '').trim();
    if (!svgId || ambiguousFeatureIds.has(svgId)) return;
    if (featuresById.has(svgId) && featuresById.get(svgId) !== feature) {
      featuresById.delete(svgId);
      ambiguousFeatureIds.add(svgId);
      return;
    }
    featuresById.set(svgId, feature);
  });
  return Array.from(svg.querySelectorAll(STANDALONE_MATCH_SELECTOR)).map((element, index) => {
    const idStatus = standaloneElementMatchId(element);
    if (!idStatus.valid) return null;
    let id = idStatus.value;
    if (!id) {
      id = `pairwise_match_${index + 1}`;
      element.setAttribute('data-gbdraw-pairwise-match-id', id);
    }
    const matchKind = standaloneMatchKind(element);
    const orthogroupId = standaloneAttr(element, 'data-orthogroup-id');
    const orthogroupIds = uniqueStandaloneMetadataValues(orthogroupId);
    const blockId = standaloneAttr(element, 'data-collinearity-block-id');
    const queryFeatureSvgId = standaloneAttr(element, 'data-query-feature-svg-id');
    const subjectFeatureSvgId = standaloneAttr(element, 'data-subject-feature-svg-id');
    const queryFeature = featuresById.get(queryFeatureSvgId) || null;
    const subjectFeature = featuresById.get(subjectFeatureSvgId) || null;
    const orthogroup = matchKind === 'collinear'
      ? null
      : getStandaloneOrthogroupPayload(
        orthogroups,
        orthogroupId,
        queryFeatureSvgId,
        subjectFeatureSvgId,
        featuresById
      ) ||
        buildStandaloneFallbackOrthogroup({
          orthogroupId,
          queryFeature,
          subjectFeature,
          features
        });
    const blockOrthogroups = matchKind === 'collinear'
      ? buildStandaloneBlockOrthogroups({
        orthogroupIds,
        orthogroups,
        featuresById,
        queryFeatureSvgId,
        subjectFeatureSvgId
      })
      : [];
    const qInterval = standaloneIntervalText(standaloneAttr(element, 'data-qstart'), standaloneAttr(element, 'data-qend'));
    const sInterval = standaloneIntervalText(standaloneAttr(element, 'data-sstart'), standaloneAttr(element, 'data-send'));
    const summaryRows = [];
    addStandaloneMatchRow(summaryRows, 'Query record', firstStandaloneText(standaloneAttr(element, 'data-query-record-id'), standaloneAttr(element, 'data-query')));
    addStandaloneMatchRow(summaryRows, 'Subject record', firstStandaloneText(standaloneAttr(element, 'data-subject-record-id'), standaloneAttr(element, 'data-subject')));
    addStandaloneMatchRow(summaryRows, 'Query interval', qInterval);
    addStandaloneMatchRow(summaryRows, 'Subject interval', sInterval);
    addStandaloneMatchRow(summaryRows, 'Orientation', firstStandaloneText(standaloneAttr(element, 'data-collinearity-orientation'), standaloneAttr(element, 'data-orientation')));

    const alignmentRows = [];
    addStandaloneMatchRow(alignmentRows, 'Identity', standaloneAttr(element, 'data-identity'));
    addStandaloneMatchRow(alignmentRows, 'Alignment length', standaloneAttr(element, 'data-alignment-length'));
    addStandaloneMatchRow(alignmentRows, 'E-value', standaloneAttr(element, 'data-evalue'));
    addStandaloneMatchRow(alignmentRows, 'Bit score', standaloneAttr(element, 'data-bitscore'));
    addStandaloneMatchRow(alignmentRows, 'Mismatches', standaloneAttr(element, 'data-mismatches'));
    addStandaloneMatchRow(alignmentRows, 'Gap opens', standaloneAttr(element, 'data-gap-opens'));

    const orthogroupDisplayName = firstStandaloneText(
      orthogroup?.display_name,
      orthogroup?.displayName,
      orthogroup?.name,
      standaloneFeatureProduct(queryFeature),
      standaloneFeatureProduct(subjectFeature)
    );
    const orthogroupRows = [];
    if (matchKind !== 'collinear') {
      addStandaloneMatchRow(orthogroupRows, 'Similarity group ID', orthogroupId);
      addStandaloneMatchRow(orthogroupRows, 'Display name', orthogroupDisplayName);
      addStandaloneMatchRow(orthogroupRows, 'Description', orthogroup?.description);
      addStandaloneMatchRow(orthogroupRows, 'Members', firstStandaloneText(orthogroup?.member_count, orthogroup?.memberCount));
      addStandaloneMatchRow(orthogroupRows, 'Record coverage', firstStandaloneText(orthogroup?.record_coverage_count, orthogroup?.recordCoverage));
    }
    const orthogroupMemberRows = buildStandaloneMatchMemberRows(orthogroup);

    const blockRows = [];
    addStandaloneMatchRow(blockRows, 'Block ID', blockId);
    addStandaloneMatchRow(blockRows, 'Kind', standaloneAttr(element, 'data-collinearity-block-kind'));
    addStandaloneMatchRow(blockRows, 'Orientation', standaloneAttr(element, 'data-collinearity-orientation'));
    addStandaloneMatchRow(blockRows, 'Color mode', standaloneAttr(element, 'data-collinearity-color-mode'));
    if (matchKind === 'collinear') {
      addStandaloneMatchRow(blockRows, 'Average identity', standaloneAttr(element, 'data-identity'));
      addStandaloneMatchRow(blockRows, 'Aligned length', standaloneAttr(element, 'data-alignment-length'));
    }
    addStandaloneMatchRow(blockRows, 'Block score', standaloneAttr(element, 'data-collinearity-block-score'));
    addStandaloneMatchRow(blockRows, 'Block e-value', standaloneAttr(element, 'data-collinearity-block-evalue'));
    addStandaloneMatchRow(blockRows, 'Anchor', [
      standaloneAttr(element, 'data-collinearity-anchor-index'),
      standaloneAttr(element, 'data-collinearity-anchor-count')
    ].filter(Boolean).join(' / '));

    if (matchKind === 'orthogroup') {
      const summarySection = {
        ...buildStandaloneMatchSection('Summary', orthogroupRows),
        member_rows: orthogroupMemberRows,
        member_copy_text: standaloneMemberCopyText(orthogroupMemberRows)
      };
      const hoverRows = [];
      addStandaloneMatchRow(hoverRows, 'Kind', matchKind);
      addStandaloneMatchRow(hoverRows, 'Similarity group', orthogroupId);
      addStandaloneMatchRow(hoverRows, 'Display name', orthogroupDisplayName);
      addStandaloneMatchRow(hoverRows, 'Members', firstStandaloneText(orthogroup?.member_count, orthogroup?.memberCount));

      return {
        id,
        title: standaloneOrthogroupTitle(orthogroupId, orthogroupDisplayName),
        subtitle: '',
        match_kind: matchKind,
        orthogroup_id: orthogroupId,
        collinearity_block_id: blockId,
        fill: firstStandaloneText(element.getAttribute('fill'), '#94a3b8'),
        sections: (summarySection.rows.length || summarySection.member_rows.length) ? [summarySection] : [],
        hover_rows: hoverRows
      };
    }

    const sections = [buildStandaloneMatchSection('Summary', summaryRows)];
    if (matchKind !== 'collinear') {
      sections.push(buildStandaloneMatchSection('Alignment', alignmentRows));
    }
    if (orthogroupRows.length || matchKind === 'orthogroup') {
      sections.push({
        ...buildStandaloneMatchSection('Similarity group', orthogroupRows),
        member_rows: orthogroupMemberRows,
        member_copy_text: standaloneMemberCopyText(orthogroupMemberRows)
      });
    }
    if (matchKind === 'collinear') {
      const blockOrthogroupRows = [];
      addStandaloneMatchRow(blockOrthogroupRows, 'Number of similarity groups covered', String(orthogroupIds.length));
      sections.push({
        ...buildStandaloneMatchSection('Similarity groups covered', blockOrthogroupRows),
        block_orthogroups: blockOrthogroups
      });
    }
    if (blockRows.length || matchKind === 'collinear') {
      sections.push(buildStandaloneMatchSection('Collinearity', blockRows));
    }
    sections.push(buildStandaloneMatchFeatureSection({
      title: 'Query',
      feature: queryFeature,
      recordId: standaloneAttr(element, 'data-query-record-id'),
      interval: qInterval,
      proteinId: standaloneAttr(element, 'data-query-protein-id'),
      locusId: standaloneAttr(element, 'data-query-locus-id'),
      displayName: standaloneAttr(element, 'data-query-display-name'),
      featureSvgIds: queryFeatureSvgId,
      featuresById,
      group: orthogroup
    }));
    sections.push(buildStandaloneMatchFeatureSection({
      title: 'Subject',
      feature: subjectFeature,
      recordId: standaloneAttr(element, 'data-subject-record-id'),
      interval: sInterval,
      proteinId: standaloneAttr(element, 'data-subject-protein-id'),
      locusId: standaloneAttr(element, 'data-subject-locus-id'),
      displayName: standaloneAttr(element, 'data-subject-display-name'),
      featureSvgIds: subjectFeatureSvgId,
      featuresById,
      group: orthogroup
    }));

    const findRow = (sectionTitle, rowLabel) => {
      const section = sections.find((entry) => entry.title === sectionTitle);
      const row = section?.rows.find((entry) => entry[0] === rowLabel);
      return row?.[1] || '';
    };
    const hoverRows = [];
    addStandaloneMatchRow(hoverRows, 'Kind', matchKind);
    addStandaloneMatchRow(hoverRows, 'Identity', findRow('Alignment', 'Identity') || findRow('Collinearity', 'Average identity'));
    addStandaloneMatchRow(hoverRows, 'Query', findRow('Summary', 'Query interval'));
    addStandaloneMatchRow(hoverRows, 'Subject', findRow('Summary', 'Subject interval'));
    if (matchKind === 'collinear') {
      addStandaloneMatchRow(hoverRows, 'Similarity groups', String(blockOrthogroups.length || orthogroupIds.length));
    } else {
      addStandaloneMatchRow(hoverRows, 'Similarity group', orthogroupId);
    }
    addStandaloneMatchRow(hoverRows, 'Block', blockId);

    return {
      id,
      title: STANDALONE_MATCH_TITLES[matchKind] || STANDALONE_MATCH_TITLES.pairwise,
      subtitle: firstStandaloneText(blockId, orthogroupId, id),
      match_kind: matchKind,
      orthogroup_id: orthogroupId,
      collinearity_block_id: blockId,
      block_orthogroup_count: blockOrthogroups.length || (matchKind === 'collinear' ? orthogroupIds.length : 0),
      block_orthogroups: blockOrthogroups,
      fill: firstStandaloneText(element.getAttribute('fill'), '#94a3b8'),
      sections: sections.filter((section) => (
        section.rows.length > 0 ||
        (Array.isArray(section.feature_rows) && section.feature_rows.length > 0)
      )),
      hover_rows: hoverRows
    };
  }).filter(Boolean);
};

const buildStandaloneMatchPayloads = (svg) => {
  if (!svg) return [];
  const payloads = Array.from(svg.querySelectorAll(STANDALONE_MATCH_SELECTOR)).map((element, index) => {
    const idStatus = standaloneElementMatchId(element);
    if (!idStatus.valid) return null;
    let id = idStatus.value;
    if (!id) {
      id = `pairwise_match_${index + 1}`;
      element.setAttribute('data-gbdraw-match-id', id);
      element.setAttribute('data-gbdraw-pairwise-match-id', id);
    }
    return compactWireValue({
      id,
      match_kind: standaloneMatchKind(element),
      orthogroup_ids: uniqueStandaloneMetadataValues(standaloneAttr(element, 'data-orthogroup-id')),
      collinearity_block_id: standaloneAttr(element, 'data-collinearity-block-id'),
      fill: firstStandaloneText(element.getAttribute('fill'), '#94a3b8'),
      query_record_id: firstStandaloneText(
        standaloneAttr(element, 'data-query-record-id'),
        standaloneAttr(element, 'data-query')
      ),
      subject_record_id: firstStandaloneText(
        standaloneAttr(element, 'data-subject-record-id'),
        standaloneAttr(element, 'data-subject')
      ),
      query_record_index: standaloneAttr(element, 'data-query-record-index'),
      subject_record_index: standaloneAttr(element, 'data-subject-record-index'),
      source_index: standaloneAttr(element, 'data-source-index'),
      track_index: standaloneAttr(element, 'data-track-index'),
      track_label: standaloneAttr(element, 'data-track-label'),
      reference_side: standaloneAttr(element, 'data-reference-side'),
      reference_record_id: standaloneAttr(element, 'data-reference-record-id'),
      qstart: standaloneAttr(element, 'data-qstart'),
      qend: standaloneAttr(element, 'data-qend'),
      sstart: standaloneAttr(element, 'data-sstart'),
      send: standaloneAttr(element, 'data-send'),
      orientation: firstStandaloneText(
        standaloneAttr(element, 'data-collinearity-orientation'),
        standaloneAttr(element, 'data-orientation')
      ),
      identity: standaloneAttr(element, 'data-identity'),
      alignment_length: standaloneAttr(element, 'data-alignment-length'),
      evalue: standaloneAttr(element, 'data-evalue'),
      bitscore: standaloneAttr(element, 'data-bitscore'),
      mismatches: standaloneAttr(element, 'data-mismatches'),
      gap_opens: standaloneAttr(element, 'data-gap-opens'),
      block_kind: standaloneAttr(element, 'data-collinearity-block-kind'),
      collinear_group_scope: standaloneAttr(element, 'data-collinear-group-scope'),
      group_kind: standaloneAttr(element, 'data-group-kind'),
      block_color_mode: standaloneAttr(element, 'data-collinearity-color-mode'),
      block_score: standaloneAttr(element, 'data-collinearity-block-score'),
      block_evalue: standaloneAttr(element, 'data-collinearity-block-evalue'),
      anchor_index: standaloneAttr(element, 'data-collinearity-anchor-index'),
      anchor_count: standaloneAttr(element, 'data-collinearity-anchor-count'),
      query_feature_svg_id: standaloneAttr(element, 'data-query-feature-svg-id'),
      subject_feature_svg_id: standaloneAttr(element, 'data-subject-feature-svg-id'),
      query_stable_feature_svg_id: standaloneAttr(
        element,
        'data-query-stable-feature-svg-id'
      ),
      subject_stable_feature_svg_id: standaloneAttr(
        element,
        'data-subject-stable-feature-svg-id'
      ),
      query_feature_index: standaloneAttr(element, 'data-query-feature-index'),
      subject_feature_index: standaloneAttr(element, 'data-subject-feature-index'),
      query_protein_id: standaloneAttr(element, 'data-query-protein-id'),
      subject_protein_id: standaloneAttr(element, 'data-subject-protein-id'),
      query_locus_id: standaloneAttr(element, 'data-query-locus-id'),
      subject_locus_id: standaloneAttr(element, 'data-subject-locus-id'),
      query_display_name: standaloneAttr(element, 'data-query-display-name'),
      subject_display_name: standaloneAttr(element, 'data-subject-display-name')
    });
  }).filter(Boolean);
  const idCounts = new Map();
  payloads.forEach((payload) => {
    idCounts.set(payload.id, (idCounts.get(payload.id) || 0) + 1);
  });
  return payloads.filter((payload) => idCounts.get(payload.id) === 1);
};

const selectStandaloneSequenceSources = (matches, sources) => {
  if (!Array.isArray(matches) || matches.length === 0) return [];
  const linearReferences = [];
  const circularRecordIds = new Set();
  const comparisonReferences = [];
  matches.forEach((match) => {
    if (match?.match_kind === 'homology') {
      const referenceSide = String(match.reference_side || '').trim();
      const recordId = standaloneConsistentTextAlias(match, [
        `${referenceSide}_record_id`,
        `${referenceSide}RecordId`
      ]);
      if (!recordId.valid) return;
      if (recordId.value) circularRecordIds.add(recordId.value);
      const sourceIndex = standaloneIntegerAlias([match.source_index]);
      if (sourceIndex.valid && sourceIndex.value !== null) {
        const comparisonSide = referenceSide === 'query' ? 'subject' : 'query';
        const comparisonRecordId = standaloneConsistentTextAlias(match, [
          `${comparisonSide}_record_id`,
          `${comparisonSide}RecordId`
        ]);
        if (comparisonRecordId.valid) {
          comparisonReferences.push({
            sourceIndex: sourceIndex.value,
            recordId: comparisonRecordId.value
          });
        }
      }
      return;
    }
    ['query', 'subject'].forEach((role) => {
      const recordIndex = standaloneIntegerAlias([
        match?.[`${role}_record_index`]
      ]);
      const recordId = standaloneConsistentTextAlias(match, [
        `${role}_record_id`,
        `${role}RecordId`
      ]);
      if (recordIndex.valid && recordIndex.value !== null && recordId.valid) {
        linearReferences.push({ recordIndex: recordIndex.value, recordId: recordId.value });
      }
    });
  });
  return (Array.isArray(sources) ? sources : []).filter((source) => {
    const origin = String(source?.origin || '').trim();
    const sourceRecordId = standaloneConsistentTextAlias(source, ['recordId', 'record_id']);
    if (!sourceRecordId.valid) return false;
    const sourceAliases = new Set([
      sourceRecordId.value,
      ...(Array.isArray(source?.aliases) ? source.aliases : [])
    ].map((value) => String(value || '').trim()).filter(Boolean));
    if (origin === 'linear-record') {
      const recordIndex = standaloneIntegerAlias([
        source.recordIndex,
        source.record_index
      ]);
      return recordIndex.valid && recordIndex.value !== null && linearReferences.some((reference) => (
        reference.recordIndex === recordIndex.value
        && (!reference.recordId || sourceAliases.has(reference.recordId))
      ));
    }
    if (origin === 'homology-comparison') {
      const sourceIndex = standaloneIntegerAlias([
        source.sourceIndex,
        source.source_index
      ]);
      return sourceIndex.valid && sourceIndex.value !== null && comparisonReferences.some((reference) => (
        reference.sourceIndex === sourceIndex.value
        && (!reference.recordId || sourceAliases.has(reference.recordId))
      ));
    }
    if (origin === 'circular-reference') {
      return Array.from(sourceAliases).some((value) => circularRecordIds.has(value));
    }
    return false;
  });
};

const appendStandaloneFeatureGlowFilter = (
  defs,
  {
    id,
    color,
    opacity,
    blurStdDeviation,
    slope,
    extent
  }
) => {
  const existing = defs.querySelector(`#${CSS.escape(id)}`);
  if (existing?.parentNode) {
    existing.parentNode.removeChild(existing);
  }

  const filter = document.createElementNS(SVG_NS, 'filter');
  filter.setAttribute('id', id);
  filter.setAttribute('x', `-${extent}%`);
  filter.setAttribute('y', `-${extent}%`);
  filter.setAttribute('width', `${100 + (extent * 2)}%`);
  filter.setAttribute('height', `${100 + (extent * 2)}%`);
  filter.setAttribute('color-interpolation-filters', 'sRGB');

  const componentTransfer = document.createElementNS(SVG_NS, 'feComponentTransfer');
  componentTransfer.setAttribute('in', 'SourceGraphic');
  componentTransfer.setAttribute('result', 'gbdrawBrightenedFeature');
  ['R', 'G', 'B'].forEach((channel) => {
    const func = document.createElementNS(SVG_NS, `feFunc${channel}`);
    func.setAttribute('type', 'linear');
    func.setAttribute('slope', String(slope));
    componentTransfer.appendChild(func);
  });

  const blur = document.createElementNS(SVG_NS, 'feGaussianBlur');
  blur.setAttribute('in', 'SourceAlpha');
  blur.setAttribute('stdDeviation', String(blurStdDeviation));
  blur.setAttribute('result', 'gbdrawFeatureGlowBlur');

  const flood = document.createElementNS(SVG_NS, 'feFlood');
  flood.setAttribute('flood-color', color);
  flood.setAttribute('flood-opacity', String(opacity));
  flood.setAttribute('result', 'gbdrawFeatureGlowColor');

  const composite = document.createElementNS(SVG_NS, 'feComposite');
  composite.setAttribute('in', 'gbdrawFeatureGlowColor');
  composite.setAttribute('in2', 'gbdrawFeatureGlowBlur');
  composite.setAttribute('operator', 'in');
  composite.setAttribute('result', 'gbdrawFeatureGlow');

  const merge = document.createElementNS(SVG_NS, 'feMerge');
  ['gbdrawFeatureGlow', 'gbdrawBrightenedFeature'].forEach((resultName) => {
    const mergeNode = document.createElementNS(SVG_NS, 'feMergeNode');
    mergeNode.setAttribute('in', resultName);
    merge.appendChild(mergeNode);
  });

  filter.appendChild(componentTransfer);
  filter.appendChild(blur);
  filter.appendChild(flood);
  filter.appendChild(composite);
  filter.appendChild(merge);
  defs.appendChild(filter);
};

const ensureStandaloneFeatureGlowFilter = (svg) => {
  const defs = ensureSvgDefs(svg);
  appendStandaloneFeatureGlowFilter(defs, {
    id: INTERACTIVE_GLOW_FILTER_ID,
    color: '#2563eb',
    opacity: 0.85,
    blurStdDeviation: 3,
    slope: 1.2,
    extent: 35
  });
  appendStandaloneFeatureGlowFilter(defs, {
    id: INTERACTIVE_MATCH_GLOW_FILTER_ID,
    color: '#fbbf24',
    opacity: 0.32,
    blurStdDeviation: 1.5,
    slope: 1.04,
    extent: 25
  });
};

const removeExistingStandaloneInteractivityAssets = (svg) => {
  [
    INTERACTIVE_METADATA_ID,
    INTERACTIVE_STYLE_ID,
    INTERACTIVE_SCRIPT_ID,
    INTERACTIVE_GLOW_FILTER_ID,
    INTERACTIVE_MATCH_GLOW_FILTER_ID,
    'gbdraw-viewport-controls',
    'gbdraw-feature-search-controls',
    'gbdraw-feature-popup',
    'gbdraw-feature-hover-popup'
  ].forEach((id) => {
    const element = svg.querySelector(`#${CSS.escape(id)}`);
    if (element?.parentNode) {
      element.parentNode.removeChild(element);
    }
  });
};

const addClassToken = (element, token) => {
  const existing = String(element.getAttribute('class') || '').trim();
  const tokens = new Set(existing ? existing.split(/\s+/) : []);
  tokens.add(token);
  element.setAttribute('class', Array.from(tokens).join(' '));
};

const removeClassToken = (element, token) => {
  const existing = String(element.getAttribute('class') || '').trim();
  const tokens = existing ? existing.split(/\s+/).filter((entry) => entry !== token) : [];
  if (tokens.length) element.setAttribute('class', tokens.join(' '));
  else element.removeAttribute('class');
};

const parseSvgViewBox = (value) => {
  const parts = String(value || '')
    .trim()
    .split(/[\s,]+/)
    .map((part) => Number(part));
  if (parts.length < 4 || parts.some((part) => !Number.isFinite(part))) return null;
  if (parts[2] <= 0 || parts[3] <= 0) return null;
  return {
    x: parts[0],
    y: parts[1],
    width: parts[2],
    height: parts[3]
  };
};

const formatSvgViewBox = (rect) =>
  [rect.x, rect.y, rect.width, rect.height].map((value) => String(value)).join(' ');

const resolveStandaloneOriginalGeometry = (svg) => {
  const dataViewBox = parseSvgViewBox(svg.getAttribute('data-gbdraw-original-viewbox'));
  const currentViewBox = parseSvgViewBox(svg.getAttribute('viewBox'));
  const width = parseFloat(svg.getAttribute('width'));
  const height = parseFloat(svg.getAttribute('height'));
  const fallbackViewBox =
    Number.isFinite(width) && Number.isFinite(height) && width > 0 && height > 0
      ? { x: 0, y: 0, width, height }
      : null;
  const viewBox = dataViewBox || currentViewBox || fallbackViewBox || { x: 0, y: 0, width: 900, height: 650 };
  const originalWidth = svg.getAttribute('data-gbdraw-original-width') ||
    svg.getAttribute('width') ||
    `${viewBox.width}px`;
  const originalHeight = svg.getAttribute('data-gbdraw-original-height') ||
    svg.getAttribute('height') ||
    `${viewBox.height}px`;
  return {
    viewBox,
    originalWidth,
    originalHeight
  };
};

const applyStandaloneViewportRoot = (svg) => {
  svg.setAttribute('width', '100vw');
  svg.setAttribute('height', '100vh');
  svg.setAttribute('preserveAspectRatio', 'xMidYMid meet');
  svg.style.setProperty('width', '100vw');
  svg.style.setProperty('height', '100vh');
  svg.style.setProperty('display', 'block');
  svg.style.setProperty('background', '#ffffff');
};

export const enrichSvgWithStandaloneInteractivity = (svg, options = {}) => {
  if (!svg) return false;

  const context = normalizeStandaloneContext(options);
  const sourceCatalogItem = selectStandaloneCatalogItem(context);
  if ((context.featureCatalog || context.requireFeatureCatalog) && !sourceCatalogItem) {
    return false;
  }
  const catalogItem = sourceCatalogItem
    ? catalogItemWithStandaloneOverrides(sourceCatalogItem, context)
    : null;
  if (sourceCatalogItem && !catalogItem) return false;
  const originalGeometry = resolveStandaloneOriginalGeometry(svg);
  const normalizedPopupMode = normalizeStandalonePopupMode(context.popupMode);
  removeExistingStandaloneInteractivityAssets(svg);
  removeClassToken(svg, 'gbdraw-feature-search-active');
  removeClassToken(svg, 'gbdraw-feature-search-updating');
  svg.querySelectorAll('.gbdraw-interactive-feature--dimmed, .gbdraw-interactive-feature--match, .gbdraw-interactive-feature--active-match')
    .forEach((element) => {
      removeClassToken(element, 'gbdraw-interactive-feature--dimmed');
      removeClassToken(element, 'gbdraw-interactive-feature--match');
      removeClassToken(element, 'gbdraw-interactive-feature--active-match');
    });
  svg.querySelectorAll('[data-gbdraw-interactive-annotation]').forEach((element) => {
    element.removeAttribute('data-gbdraw-interactive-annotation');
    removeClassToken(element, 'gbdraw-interactive-annotation');
  });
  const features = catalogItem
    ? catalogItem.features
    : buildStandaloneFeaturePayloads(svg, {
      ...context,
      popupMode: normalizedPopupMode
    });
  const biologicalFeatures = catalogItem
    ? catalogItem.biologicalFeatures
    : buildStandaloneBiologicalFeaturePayloads(context, features);
  const orthogroups = catalogItem
    ? catalogItem.orthogroups
    : buildStandaloneOrthogroupPayloads(features, context);
  const matches = catalogItem
    ? catalogItem.comparisonMatches
    : buildStandaloneMatchPayloads(svg, { features, orthogroups });
  const sequenceSources = catalogItem
    ? (catalogItem.sequenceSources || [])
    : selectStandaloneSequenceSources(matches, context.sequenceSources);

  const featureIds = new Set(features.map((feature) => (
    catalogItem ? feature.svgId : feature.svg_id
  )));
  svg.querySelectorAll(FEATURE_SELECTOR).forEach((element) => {
    const id = getElementFeatureId(element);
    if (!featureIds.has(id)) return;
    element.setAttribute('data-gbdraw-interactive-feature', 'true');
    addClassToken(element, 'gbdraw-interactive-feature');
  });
  const matchIds = new Set(matches.map((match) => match.id));
  svg.querySelectorAll(STANDALONE_MATCH_SELECTOR).forEach((element) => {
    const idStatus = standaloneElementMatchId(element);
    const id = idStatus.valid ? idStatus.value : '';
    if (!matchIds.has(id)) return;
    element.setAttribute('data-gbdraw-interactive-match', 'true');
    addClassToken(element, 'gbdraw-interactive-pairwise-match');
  });
  const annotationDomIds = new Set(
    (Array.isArray(catalogItem?.annotations) ? catalogItem.annotations : [])
      .map((annotation) => String(annotation?.dom_id || annotation?.domId || '').trim())
      .filter(Boolean)
  );
  svg.querySelectorAll(STANDALONE_ANNOTATION_SELECTOR).forEach((element) => {
    if (!annotationDomIds.has(String(element.id || '').trim())) return;
    element.setAttribute('data-gbdraw-interactive-annotation', 'true');
    addClassToken(element, 'gbdraw-interactive-annotation');
  });

  const metadata = document.createElementNS(SVG_NS, 'metadata');
  metadata.setAttribute('id', INTERACTIVE_METADATA_ID);
  metadata.setAttribute('data-popup-mode', normalizedPopupMode);
  if (catalogItem) {
    metadata.setAttribute('data-schema', String(INTERACTIVE_CATALOG_SCHEMA));
    metadata.setAttribute('data-result-index', String(catalogItem.resultIndex));
    metadata.setAttribute('data-result-name', String(catalogItem.resultName || ''));
    metadata.textContent = JSON.stringify({
      schema: INTERACTIVE_CATALOG_SCHEMA,
      items: [catalogItem]
    });
  } else {
    metadata.setAttribute('data-schema', INTERACTIVE_SCHEMA);
    metadata.textContent = JSON.stringify(compactWireValue({
      schema: INTERACTIVE_SCHEMA,
      popup_mode: normalizedPopupMode,
      features,
      biological_features: biologicalFeatures,
      orthogroups,
      matches,
      sequence_sources: sequenceSources
    }));
  }

  const style = document.createElementNS(SVG_NS, 'style');
  style.setAttribute('id', INTERACTIVE_STYLE_ID);
  style.setAttribute('type', 'text/css');
  style.textContent = STANDALONE_INTERACTIVE_STYLE;

  const script = document.createElementNS(SVG_NS, 'script');
  script.setAttribute('id', INTERACTIVE_SCRIPT_ID);
  script.setAttribute('type', 'application/ecmascript');
  script.textContent = STANDALONE_INTERACTIVE_SCRIPT;

  ensureStandaloneFeatureGlowFilter(svg);
  svg.appendChild(metadata);
  svg.appendChild(style);
  svg.appendChild(script);
  svg.setAttribute('data-gbdraw-original-viewbox', formatSvgViewBox(originalGeometry.viewBox));
  svg.setAttribute('data-gbdraw-original-width', originalGeometry.originalWidth);
  svg.setAttribute('data-gbdraw-original-height', originalGeometry.originalHeight);
  svg.setAttribute('data-gbdraw-interactive-svg', 'true');
  applyStandaloneViewportRoot(svg);
  return true;
};
