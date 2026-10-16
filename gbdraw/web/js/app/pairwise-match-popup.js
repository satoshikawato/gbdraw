// @ts-check
import {
  formatFeatureLocation,
  isInternalProteinDisplayId,
  locationRangeText,
  resolveDisplayProteinId
} from '../services/feature-utils.js';
import { buildFeatureSequenceFastas } from '../services/feature-sequence-fasta.js';
import { buildMatchSequenceBundle } from '../services/match-sequences.js';
import { readRecordSourceSpan, recordSourceInterval } from '../services/record-source-coordinates.js';
import {
  groupMetadataScopeLabel,
  normalizeGroupMetadataScope
} from '../services/losat-normalization.js';
import {
  RENDERED_FEATURE_ID_KEYS,
  STABLE_FEATURE_ID_KEYS,
  featureIdentity,
  identityMatches,
  renderedFeatureIdentity,
  textAliasStatus,
  orthogroupIdStatus
} from '../services/feature-identity.js';

/** @import { AliasSource, FeatureIdentity } from '../services/feature-identity.js' */
/** @import { FeatureLike } from '../services/feature-utils.js' */
/** @import { SequenceSource, SequenceSourceResolution, SequenceSourceResolver } from '../services/match-sequences.js' */

/**
 * A rendered or source feature row of the catalogs, or a saved group member:
 * its fields are read by name, with snake_case aliases.
 * @typedef {FeatureLike & { orthogroupMember?: PopupMember | null, orthogroup_member?: PopupMember | null }} PopupFeature
 */

/**
 * A member record of a saved similarity group.
 * @typedef {FeatureLike} PopupMember
 */

/**
 * A saved similarity group: its fields are read by name, with snake_case aliases.
 * @typedef {AliasSource} PopupGroup
 */

/** @typedef {'pairwise' | 'orthogroup' | 'collinear' | 'homology'} MatchKind */
/** @typedef {'nt' | 'aa'} SequenceKind */
/** @typedef {ReturnType<typeof buildOrthogroupMemberRows>[number]} MemberRow One row of a group's member table. */
/** @typedef {ReturnType<typeof buildFeatureListRows>[number]} FeatureListRow One row of a feature list. */
/** @typedef {ReturnType<typeof buildGroupView>} GroupView The view of one similarity group. */

/**
 * One popup table row.
 * @typedef {{ label: string, value: string }} PopupRow
 */

/**
 * A popup section. `section` builds `title` and `rows`; the member and
 * feature sections add their list data and the Similarity group section
 * adds the copy and FASTA fields.
 * @typedef {object} PopupSection
 * @property {string} title
 * @property {PopupRow[]} rows
 * @property {FeatureListRow[]} [featureRows]
 * @property {MemberRow[]} [memberRows]
 * @property {string} [memberCopyText]
 * @property {string} [memberNtFasta]
 * @property {string} [memberAaFasta]
 * @property {number} [memberNtCount]
 * @property {number} [memberAaCount]
 * @property {string} [memberNtFilename]
 * @property {string} [memberAaFilename]
 * @property {GroupView[]} [blockOrthogroups]
 */

/**
 * What the popup and the hover summary read from one match element.
 * @typedef {object} PairwiseMatchDescriptor
 * @property {MatchKind} matchKind
 * @property {string} matchId
 * @property {string} orthogroupId
 * @property {string[]} orthogroupIds
 * @property {string} collinearityBlockId
 * @property {string} groupScope
 * @property {string} normalizedGroupScope
 * @property {string} queryFeatureSvgId
 * @property {string} subjectFeatureSvgId
 * @property {string} queryRecordId
 * @property {string} subjectRecordId
 * @property {string} queryInterval
 * @property {string} subjectInterval
 * @property {string} queryTableInterval
 * @property {string} subjectTableInterval
 * @property {string} identity
 * @property {string} orientation
 * @property {string} title
 * @property {string} subtitle
 * @property {string} fill
 * @property {{ query: ReturnType<typeof matchEndpointReferences>, subject: ReturnType<typeof matchEndpointReferences> }} endpoints
 */

/**
 * The values of one hover summary, taken from a descriptor or from a payload.
 * @typedef {object} PairwiseHoverRowValues
 * @property {MatchKind} [matchKind]
 * @property {string} [identity]
 * @property {string} [queryInterval]
 * @property {string} [subjectInterval]
 * @property {string} [orthogroupId]
 * @property {string} [displayName]
 * @property {string} [memberCount]
 * @property {string} [groupScope]
 * @property {number | string} [blockOrthogroupCount]
 * @property {string} [collinearityBlockId]
 */

/**
 * Overrides keyed by similarity-group ID, as a Map or a plain object.
 * @typedef {Map<string, string> | Record<string, string> | null} OrthogroupTextOverrides
 */

/**
 * Rendered features by SVG ID, as the SVG owner holds them.
 * @typedef {Map<string, PopupFeature> | Record<string, PopupFeature> | PopupFeature[] | null} FeatureLookup
 */

/**
 * What the match popup shows (`buildMatchPopupPayload`). A similarity-group
 * match has a `group`, a collinearity block a `block` and its groups, and the
 * other kinds the `sections` and the matched sequences.
 * @typedef {object} MatchPopupPayload
 * @property {string} id
 * @property {string} title
 * @property {string} subtitle
 * @property {MatchKind} matchKind
 * @property {string} orthogroupId
 * @property {string} collinearityBlockId
 * @property {string} queryFeatureSvgId
 * @property {string} subjectFeatureSvgId
 * @property {string} fill
 * @property {GroupView | null} [group]
 * @property {PopupSection[]} sections
 * @property {string} [groupScope]
 * @property {number} [blockOrthogroupCount]
 * @property {GroupView[]} [blockOrthogroups]
 * @property {ReturnType<typeof buildCollinearBlockView>} [block]
 * @property {ReturnType<typeof buildSequenceBundleForMatch>} [sequenceBundle]
 * @property {string} [sequenceTitle]
 */

/**
 * @typedef {object} MatchPopupPayloadOptions
 * @property {FeatureLookup} [featureLookup]
 * @property {PopupFeature[]} [sourceFeatures]
 * @property {PopupGroup[]} [orthogroups]
 * @property {OrthogroupTextOverrides} [orthogroupNameOverrides]
 * @property {OrthogroupTextOverrides} [orthogroupDescriptionOverrides]
 * @property {import('../services/match-sequences.js').SequenceSourceResolver | null} [resolveSequenceSource]
 */

/**
 * @typedef {object} MatchHoverSummaryOptions
 * @property {PopupGroup[] | (() => PopupGroup[])} [orthogroups]
 * @property {OrthogroupTextOverrides} [orthogroupNameOverrides]
 */

export const PAIRWISE_MATCH_SELECTOR = [
  'path[data-gbdraw-match-id]',
  'path[data-gbdraw-pairwise-match-id]',
  'path[data-match-kind]',
  'path[data-pairwise-match-style]'
].join(', ');

const MATCH_KIND_TITLES = {
  pairwise: 'Pairwise match',
  orthogroup: 'Similarity-group match',
  collinear: 'Collinearity block',
  homology: 'Homology ring match'
};

const MATCH_KIND_ALIASES = new Set(['pairwise', 'orthogroup', 'collinear', 'homology']);
/** @param {unknown} value */
const normalizeText = (value) => String(value === null || value === undefined ? '' : value).trim();
/** @param {unknown} value */
const normalizeSequence = (value) => String(value ?? '').replace(/\s+/g, '').toUpperCase();

/**
 * @param {Element | null | undefined} element
 * @param {string} name
 */
const attr = (element, name) => normalizeText(element?.getAttribute?.(name));
/**
 * @param {Element | null | undefined} element
 * @param {string} name
 */
const integerAttrStatus = (element, name) => {
  const value = attr(element, name);
  if (!value) return { valid: true, supplied: false, value: null };
  if (!/^\d+$/.test(value)) {
    return { valid: false, supplied: true, value: null };
  }
  const numeric = Number(value);
  const valid = Number.isSafeInteger(numeric) && numeric >= 0;
  return { valid, supplied: true, value: valid ? numeric : null };
};
/**
 * @param {Element | null | undefined} element
 * @param {string} name
 */
const integerAttr = (element, name) => integerAttrStatus(element, name).value;

/** @param {...unknown} values */
const firstText = (...values) => {
  for (const value of values) {
    const text = normalizeText(value);
    if (text) return text;
  }
  return '';
};

/** @param {unknown} value */
const splitMetadataValues = (value) => normalizeText(value)
  .split(';')
  .map((entry) => normalizeText(entry))
  .filter(Boolean);

/** @param {unknown} value */
const uniqueMetadataValues = (value) => {
  /** @type {Set<string>} */
  const seen = new Set();
  /** @type {string[]} */
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

/** @param {unknown} value */
const isInternalDisplayId = (value) => {
  const text = normalizeText(value);
  return Boolean(text && (
    isInternalProteinDisplayId(text) ||
    generatedProteinIdPattern.test(text) ||
    generatedUnitIdPattern.test(text)
  ));
};

/**
 * @param {string[]} values
 * @param {unknown} value
 */
const addUniqueDisplayText = (values, value) => {
  const text = normalizeText(value);
  if (!text || isInternalDisplayId(text) || values.includes(text)) return;
  values.push(text);
};

/** @param {...unknown} values */
const firstSequenceText = (...values) => {
  for (const value of values) {
    const sequence = normalizeSequence(value);
    if (sequence) return sequence;
  }
  return '';
};

/**
 * @param {unknown} value
 * @param {string} [fallback]
 */
const makeSafeFilename = (value, fallback = 'orthogroup') => {
  const cleaned = normalizeText(value).replace(/[^\w.-]+/g, '_').replace(/^_+|_+$/g, '');
  return cleaned || fallback;
};

/** @param {SequenceKind} sequenceKind */
const sequenceKindLabel = (sequenceKind) => (sequenceKind === 'aa' ? 'aa' : 'nt');
/** @param {SequenceKind} sequenceKind */
const sequenceExtension = (sequenceKind) => (sequenceKindLabel(sequenceKind) === 'aa' ? 'faa' : 'fna');

/**
 * @param {unknown} value
 * @param {SVGElement | null} [element]
 * @returns {MatchKind}
 */
const normalizeMatchKind = (value, element = null) => {
  const normalized = normalizeText(value).toLowerCase();
  // A member of `MATCH_KIND_ALIASES`, by the test.
  if (MATCH_KIND_ALIASES.has(normalized)) return /** @type {MatchKind} */ (normalized);
  if (attr(element, 'data-collinearity-block-id')) return 'collinear';
  if (attr(element, 'data-orthogroup-id')) return 'orthogroup';
  return 'pairwise';
};

/**
 * @param {PopupRow[]} rows
 * @param {string} label
 * @param {unknown} value
 */
const addRow = (rows, label, value) => {
  const text = normalizeText(value);
  if (!text) return;
  rows.push({ label, value: text });
};

/**
 * @param {unknown} start
 * @param {unknown} end
 */
const intervalText = (start, end) => {
  const startText = normalizeText(start);
  const endText = normalizeText(end);
  if (!startText && !endText) return '';
  if (startText && endText) return `${startText}..${endText}`;
  return startText || endText;
};

/**
 * @param {PopupFeature | null | undefined} feature
 * @param {string} key
 */
const qualifierFirstValue = (feature, key) => {
  const normalizedKey = normalizeText(key).toLowerCase();
  if (!feature || !normalizedKey) return '';
  const direct = feature[normalizedKey];
  if (direct !== null && direct !== undefined && direct !== '') return firstArrayValue(direct);
  // A plain object, by the test.
  const qualifiers = /** @type {Record<string, unknown>} */ (
    feature.qualifiers && typeof feature.qualifiers === 'object' && !Array.isArray(feature.qualifiers)
      ? feature.qualifiers
      : {}
  );
  if (qualifiers[normalizedKey] !== null && qualifiers[normalizedKey] !== undefined) {
    return firstArrayValue(qualifiers[normalizedKey]);
  }
  const matchingKey = Object.keys(qualifiers).find((candidate) => candidate.toLowerCase() === normalizedKey);
  return matchingKey ? firstArrayValue(qualifiers[matchingKey]) : '';
};

/** @param {unknown} value */
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

/** @param {PopupGroup | null | undefined} group */
const groupScopeStatus = (group) => {
  const presentationValues = [
    group?.presentationScope,
    group?.presentation_scope,
    group?.collinearGroupScope,
    group?.collinear_group_scope,
    group?.groupScope,
    group?.group_scope
  ]
    .map(normalizeText)
    .filter(Boolean)
    .map(normalizeGroupMetadataScope);
  if (normalizeText(group?.groupKind ?? group?.group_kind).toLowerCase() === 'collinear_gene_group') {
    presentationValues.push('adjacent_local');
  }
  const presentationScopes = new Set(presentationValues);
  if (presentationScopes.size > 1) {
    return { valid: false, supplied: true, value: '' };
  }
  if (presentationScopes.size === 1) {
    return { valid: true, supplied: true, value: presentationScopes.values().next().value };
  }
  const coreScope = normalizeText(group?.scope);
  return {
    valid: true,
    supplied: Boolean(coreScope),
    value: coreScope ? normalizeGroupMetadataScope(coreScope) : ''
  };
};

/** @param {PopupGroup | null | undefined} group */
const groupScopeValue = (group) => {
  const status = groupScopeStatus(group);
  return status.valid ? status.value : '';
};

/**
 * @template V
 * @param {Map<string, V[]>} buckets
 * @param {string} key
 * @param {V} value
 */
const addBucketValue = (buckets, key, value) => {
  if (!key) return;
  if (!buckets.has(key)) buckets.set(key, []);
  /** @type {V[]} */ (buckets.get(key)).push(value); // Set above when absent.
};

/**
 * The similarity groups by ID, built once per groups array (OV-930). A group
 * answers a lookup only when its ID and scope are valid; `identified` keeps,
 * in order, every group whose ID is valid and supplied, whatever its scope.
 * @typedef {object} OrthogroupIndex
 * @property {Map<string, PopupGroup[]>} byId
 * @property {Map<string, PopupGroup[]>} byScopedId Keyed by ID and normalized scope.
 * @property {PopupGroup[]} identified
 */

/** @type {WeakMap<PopupGroup[], OrthogroupIndex>} */
const orthogroupIndexes = new WeakMap();

/**
 * The index of a groups array. The groups are replaced as a whole array, never
 * edited in place, so the array's identity keys its index.
 * @param {PopupGroup[] | null | undefined} orthogroups
 * @returns {OrthogroupIndex}
 */
const orthogroupIndexFor = (orthogroups) => {
  const groups = Array.isArray(orthogroups) ? orthogroups : null;
  const cached = groups ? orthogroupIndexes.get(groups) : null;
  if (cached) return cached;
  /** @type {OrthogroupIndex} */
  const index = { byId: new Map(), byScopedId: new Map(), identified: [] };
  for (const group of groups || []) {
    const idStatus = orthogroupIdStatus(group);
    if (!idStatus.valid || !idStatus.value) continue;
    index.identified.push(group);
    const scopeStatus = groupScopeStatus(group);
    if (!scopeStatus.valid) continue;
    addBucketValue(index.byId, idStatus.value, group);
    addBucketValue(index.byScopedId, `${idStatus.value}\u001f${scopeStatus.value}`, group);
  }
  if (groups) orthogroupIndexes.set(groups, index);
  return index;
};

/**
 * The one group with this ID (and this scope, when one is asked for), or null
 * when none or several are.
 * @param {OrthogroupIndex} index
 * @param {unknown} orthogroupId
 * @param {unknown} [groupScope]
 * @returns {PopupGroup | null}
 */
const orthogroupInIndex = (index, orthogroupId, groupScope = '') => {
  const id = normalizeText(orthogroupId);
  if (!id) return null;
  const requestedScope = normalizeText(groupScope)
    ? normalizeGroupMetadataScope(groupScope)
    : '';
  const candidates = requestedScope
    ? index.byScopedId.get(`${id}\u001f${requestedScope}`) || []
    : index.byId.get(id) || [];
  return candidates.length === 1 ? candidates[0] : null;
};

/** @param {AliasSource | null | undefined} member */
const memberFeatureSvgId = (member) => {
  const status = textAliasStatus(member, STABLE_FEATURE_ID_KEYS);
  return status.valid ? status.value : '';
};

/** @param {AliasSource | null | undefined} feature */
const featureRenderedSvgId = (feature) => {
  const status = textAliasStatus(feature, [
    ...RENDERED_FEATURE_ID_KEYS,
    'svg_id',
    'svgId'
  ]);
  return status.valid ? status.value : '';
};

/**
 * @param {FeatureLookup} featureLookup
 * @returns {PopupFeature[]}
 */
const getFeatureLookupValues = (featureLookup) => {
  if (!featureLookup) return [];
  if (typeof featureLookup.values === 'function') return Array.from(featureLookup.values());
  if (Array.isArray(featureLookup)) return featureLookup;
  if (typeof featureLookup === 'object') return Object.values(featureLookup);
  return [];
};

/**
 * @param {Element} element
 * @param {'query' | 'subject'} role
 */
const matchEndpointReferences = (element, role) => {
  const renderedIds = splitMetadataValues(attr(element, `data-${role}-feature-svg-id`));
  const stableText = attr(element, `data-${role}-stable-feature-svg-id`);
  const sourceIndexText = attr(element, `data-${role}-feature-index`);
  const recordText = attr(element, `data-${role}-record-index`);
  const constrained = renderedIds.length > 0;
  if (!constrained) {
    return { constrained: false, valid: false, references: [] };
  }
  const stableIds = splitMetadataValues(stableText);
  const sourceIndexes = splitMetadataValues(sourceIndexText);
  if (
    (stableIds.length > 0 && stableIds.length !== renderedIds.length) ||
    (sourceIndexes.length > 0 && sourceIndexes.length !== renderedIds.length)
  ) return { constrained, valid: false, references: [] };

  const references = renderedIds.map((renderedId, index) => featureIdentity({
    recordIndex: recordText,
    featureIndex: sourceIndexes[index],
    stableFeatureSvgId: stableIds[index],
    renderedFeatureSvgId: renderedId
  }));
  return {
    constrained,
    valid: references.every((reference) => (
      reference.usable &&
      reference.recordIndex.supplied &&
      reference.renderedId.supplied &&
      (reference.stableId.supplied || reference.sourceIndex.supplied)
    )),
    references
  };
};

// Source coordinates first; the table (search-frame) interval only when it
// differs (PD-OI-076). Circular homology matches have no Linear record span.
/**
 * @param {Element} element
 * @param {MatchKind} matchKind
 * @param {'query' | 'subject'} role
 * @param {string} start
 * @param {string} end
 */
const matchIntervals = (element, matchKind, role, start, end) => {
  const source = matchKind === 'homology'
    ? null
    : recordSourceInterval(
      readRecordSourceSpan(element, attr(element, `data-${role}-record-index`)),
      start,
      end
    );
  return source
    ? { interval: intervalText(source.start, source.end), table: source.table ? intervalText(source.table.start, source.table.end) : '' }
    : { interval: intervalText(start, end), table: '' };
};

/**
 * @param {SVGElement} element
 * @returns {PairwiseMatchDescriptor}
 */
const readPairwiseMatchDescriptor = (element) => {
  const matchKind = normalizeMatchKind(attr(element, 'data-match-kind'), element);
  const orthogroupId = attr(element, 'data-orthogroup-id');
  const collinearityBlockId = attr(element, 'data-collinearity-block-id');
  const groupScope = firstText(
    attr(element, 'data-collinear-group-scope'),
    attr(element, 'data-group-scope')
  );
  const queryFeatureSvgId = attr(element, 'data-query-feature-svg-id');
  const subjectFeatureSvgId = attr(element, 'data-subject-feature-svg-id');
  const queryStart = attr(element, 'data-qstart');
  const queryEnd = attr(element, 'data-qend');
  const subjectStart = attr(element, 'data-sstart');
  const subjectEnd = attr(element, 'data-send');
  const queryIntervals = matchIntervals(element, matchKind, 'query', queryStart, queryEnd);
  const subjectIntervals = matchIntervals(element, matchKind, 'subject', subjectStart, subjectEnd);
  const matchId = firstText(
    attr(element, 'data-gbdraw-match-id'),
    attr(element, 'data-gbdraw-pairwise-match-id'),
    collinearityBlockId,
    orthogroupId
  );
  return {
    matchKind,
    matchId,
    orthogroupId,
    orthogroupIds: uniqueMetadataValues(orthogroupId),
    collinearityBlockId,
    groupScope,
    normalizedGroupScope: normalizeGroupMetadataScope(groupScope),
    queryFeatureSvgId,
    subjectFeatureSvgId,
    queryRecordId: attr(element, 'data-query-record-id'),
    subjectRecordId: attr(element, 'data-subject-record-id'),
    queryInterval: queryIntervals.interval,
    subjectInterval: subjectIntervals.interval,
    queryTableInterval: queryIntervals.table,
    subjectTableInterval: subjectIntervals.table,
    identity: attr(element, 'data-identity'),
    orientation: firstText(
      attr(element, 'data-collinearity-orientation'),
      attr(element, 'data-orientation')
    ),
    title: MATCH_KIND_TITLES[matchKind] || MATCH_KIND_TITLES.pairwise,
    subtitle: firstText(collinearityBlockId, orthogroupId, matchId),
    fill: firstText(attr(element, 'fill'), element?.style?.fill, '#94a3b8'),
    endpoints: {
      query: matchEndpointReferences(element, 'query'),
      subject: matchEndpointReferences(element, 'subject')
    }
  };
};

/** @param {PopupFeature | null | undefined} feature */
const featureOrthogroupIdStatus = (feature) => {
  const statuses = [
    orthogroupIdStatus(feature),
    orthogroupIdStatus(feature?.orthogroupMember),
    orthogroupIdStatus(feature?.orthogroup_member)
  ];
  if (statuses.some((status) => !status.valid)) {
    return { valid: false, supplied: true, value: '' };
  }
  const values = new Set(statuses.map((status) => status.value).filter(Boolean));
  return {
    valid: values.size <= 1,
    supplied: statuses.some((status) => status.supplied),
    value: values.size === 1 ? values.values().next().value : ''
  };
};

/** @param {PopupFeature | null | undefined} feature */
const featureProduct = (feature) => firstText(
  feature?.product,
  qualifierFirstValue(feature, 'product'),
  feature?.orthogroupMember?.product,
  feature?.orthogroup_member?.product,
  feature?.note,
  qualifierFirstValue(feature, 'note')
);

/**
 * @param {PopupFeature | null | undefined} feature
 * @param {FeatureIdentity | null} [resolvedIdentity]
 */
const featureToOrthogroupMember = (feature, resolvedIdentity = null) => {
  if (!feature) return null;
  const member = feature?.orthogroupMember || feature?.orthogroup_member || {};
  let identity = resolvedIdentity || featureIdentity(feature);
  if (!resolvedIdentity && !identity.usable) {
    identity = featureIdentity(feature, { allowLegacySvgStable: true });
  }
  return {
    recordIndex: identity.recordIndex.value,
    recordKey: firstText(member?.recordKey, member?.record_key, feature?.recordKey, feature?.record_key),
    biologicalFeatureId: firstText(
      member?.biologicalFeatureId,
      member?.biological_feature_id,
      feature?.biologicalFeatureId,
      feature?.biological_feature_id
    ),
    featureIndex: identity.sourceIndex.value,
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
    featureSvgId: firstText(
      memberFeatureSvgId(member),
      identity.stableId.value
    ),
    renderedFeatureSvgId: featureRenderedSvgId(feature)
  };
};

/** @param {FeatureIdentity | null | undefined} identity */
const identityIndexKeys = (identity) => {
  if (!identity?.usable) return [];
  const recordIndex = identity.recordIndex.supplied ? identity.recordIndex.value : null;
  /** @type {string[]} */
  const keys = [];
  if (identity.recordKey.supplied && identity.biologicalId.supplied) {
    keys.push(`canonical:${identity.recordKey.value}\u001f${identity.biologicalId.value}`);
  }
  if (recordIndex !== null && identity.sourceIndex.supplied) {
    keys.push(`source:${recordIndex}\u001f${identity.sourceIndex.value}`);
  }
  if (recordIndex !== null && identity.stableId.supplied) {
    keys.push(`stable:${recordIndex}\u001f${identity.stableId.value}`);
  }
  return keys;
};

/**
 * A rendered feature and its identity from its rendered ID.
 * @typedef {{ feature: PopupFeature, identity: FeatureIdentity }} RenderedRecord
 */

/**
 * A source feature, how a member reaches it, and the rendered ID it names.
 * @typedef {object} SourceRecord
 * @property {PopupFeature} feature
 * @property {FeatureIdentity | null | undefined} candidateIdentity
 * @property {FeatureIdentity | null} sourceIdentity
 * @property {ReturnType<typeof textAliasStatus>} sourceRendered
 * @property {ReturnType<typeof textAliasStatus>} legacyRendered
 */

/**
 * A group member and its identity.
 * @typedef {{ member: PopupMember, identity: FeatureIdentity }} MemberRecord
 */

/**
 * The lookups of one payload (`createPairwisePayloadContext`): the features,
 * groups, and members that the match's IDs name, each answered only when one
 * candidate matches.
 * @typedef {object} PairwisePayloadContext
 * @property {(orthogroupId: unknown, queryFeature: PopupFeature | null, subjectFeature: PopupFeature | null) => PopupFeature[]} fallbackFeaturesForOrthogroup
 * @property {(renderedSvgId: unknown) => PopupFeature | null} featureForLookupId
 * @property {(member: PopupMember) => PopupFeature | null} featureForMember
 * @property {(group: PopupGroup | null | undefined, featureSvgIds: unknown) => boolean} groupHasFeatureIds
 * @property {(group: PopupGroup | null | undefined, featureSvgId: unknown) => PopupMember | null} groupMemberForFeatureId
 * @property {(group: PopupGroup | null | undefined) => MemberRecord[]} groupMemberRecords
 * @property {(feature: PopupFeature) => FeatureIdentity} identityForFeature
 * @property {(orthogroupId: unknown, groupScope?: unknown) => PopupGroup | null} orthogroupById
 * @property {() => PopupGroup | null} orthogroupForMatch
 * @property {(member: PopupMember, feature: PopupFeature | null) => PopupFeature | null} renderedFeatureForMember
 */

/**
 * The lookups of a payload that depend only on the rendered and the source
 * features, and the member index of each group they meet. Built once per
 * (`featureLookup`, `sourceFeatures`) pair and reused by every later popup on
 * the same Result (OV-911).
 * @typedef {object} PayloadFeatureIndex
 * @property {(source: PopupFeature | PopupMember | null | undefined, mode: 'rendered' | 'standard' | 'legacy') => FeatureIdentity} cachedIdentity
 * @property {WeakMap<object, RenderedRecord>} renderedByObject
 * @property {Map<string, RenderedRecord[]>} renderedById
 * @property {Map<string, PopupFeature[]>} featuresByOrthogroupId
 * @property {WeakMap<object, ReturnType<typeof featureOrthogroupIdStatus>>} featureOrthogroupStatuses
 * @property {WeakMap<object, SourceRecord>} sourceByObject
 * @property {Map<string, SourceRecord[]>} sourceByIdentity
 * @property {(group: PopupGroup | null | undefined) => { records: MemberRecord[], byIdentity: Map<string, MemberRecord[]> }} memberIndexForGroup
 */

/**
 * @param {FeatureLookup} featureLookup
 * @param {PopupFeature[]} sourceFeatures
 * @returns {PayloadFeatureIndex}
 */
const buildPayloadFeatureIndex = (featureLookup, sourceFeatures) => {
  /** @type {WeakMap<object, Map<string, FeatureIdentity>>} */
  const identityCache = new WeakMap();
  /**
   * @param {PopupFeature | PopupMember | null | undefined} source
   * @param {'rendered' | 'standard' | 'legacy'} mode
   * @returns {FeatureIdentity}
   */
  const cachedIdentity = (source, mode) => {
    if (!source || (typeof source !== 'object' && typeof source !== 'function')) {
      if (mode === 'rendered') return renderedFeatureIdentity(source);
      return featureIdentity(source, { allowLegacySvgStable: mode === 'legacy' });
    }
    let modes = identityCache.get(source);
    if (!modes) {
      modes = new Map();
      identityCache.set(source, modes);
    }
    if (!modes.has(mode)) {
      modes.set(
        mode,
        mode === 'rendered'
          ? renderedFeatureIdentity(source)
          : featureIdentity(source, { allowLegacySvgStable: mode === 'legacy' })
      );
    }
    return /** @type {FeatureIdentity} */ (modes.get(mode)); // Set above when absent.
  };

  /** @type {RenderedRecord[]} */
  const renderedRecords = [];
  /** @type {WeakMap<object, RenderedRecord>} */
  const renderedByObject = new WeakMap();
  /** @type {Map<string, RenderedRecord[]>} */
  const renderedById = new Map();
  /** @type {Map<string, PopupFeature[]>} */
  const featuresByOrthogroupId = new Map();
  /** @type {WeakMap<object, ReturnType<typeof featureOrthogroupIdStatus>>} */
  const featureOrthogroupStatuses = new WeakMap();
  /** @type {Set<PopupFeature>} */
  const seenRenderedFeatures = new Set();
  for (const feature of getFeatureLookupValues(featureLookup)) {
    if (!feature) continue;
    const groupStatus = featureOrthogroupIdStatus(feature);
    featureOrthogroupStatuses.set(feature, groupStatus);
    if (groupStatus.valid && groupStatus.value) {
      addBucketValue(featuresByOrthogroupId, groupStatus.value, feature);
    }
    if (seenRenderedFeatures.has(feature)) continue;
    seenRenderedFeatures.add(feature);
    const identity = cachedIdentity(feature, 'rendered');
    const record = { feature, identity };
    renderedRecords.push(record);
    renderedByObject.set(feature, record);
    if (identity.usable && identity.renderedId.value) {
      addBucketValue(renderedById, identity.renderedId.value, record);
    }
  }

  const sourceIsBiological = Array.isArray(sourceFeatures) && sourceFeatures.length > 0;
  /** @type {WeakMap<object, SourceRecord>} */
  const sourceByObject = new WeakMap();
  /** @type {Map<string, SourceRecord[]>} */
  const sourceByIdentity = new Map();
  const sourceValues = sourceIsBiological
    ? sourceFeatures
    : renderedRecords.map((record) => record.feature);
  for (const feature of sourceValues) {
    if (!feature) continue;
    const standardIdentity = cachedIdentity(feature, 'standard');
    const legacyIdentity = sourceIsBiological
      ? cachedIdentity(feature, 'legacy')
      : null;
    const candidateIdentity = sourceIsBiological
      ? legacyIdentity
      : renderedByObject.get(feature)?.identity;
    /** @type {SourceRecord} */
    const record = {
      feature,
      candidateIdentity,
      sourceIdentity: standardIdentity.usable ? standardIdentity : legacyIdentity,
      sourceRendered: textAliasStatus(feature, RENDERED_FEATURE_ID_KEYS),
      legacyRendered: textAliasStatus(feature, ['svg_id', 'svgId'])
    };
    if (!sourceByObject.has(feature)) sourceByObject.set(feature, record);
    identityIndexKeys(candidateIdentity).forEach((key) => {
      addBucketValue(sourceByIdentity, key, record);
    });
  }

  /** @type {WeakMap<object, { records: MemberRecord[], byIdentity: Map<string, MemberRecord[]> }>} */
  const groupMemberIndexes = new WeakMap();
  /**
   * @param {PopupGroup | null | undefined} group
   * @returns {{ records: MemberRecord[], byIdentity: Map<string, MemberRecord[]> }}
   */
  const memberIndexForGroup = (group) => {
    if (!group || typeof group !== 'object') return { records: [], byIdentity: new Map() };
    const cached = groupMemberIndexes.get(group);
    if (cached) return cached;
    /** @type {{ records: MemberRecord[], byIdentity: Map<string, MemberRecord[]> }} */
    const index = { records: [], byIdentity: new Map() };
    for (const member of Array.isArray(group.members) ? group.members : []) {
      const record = { member, identity: cachedIdentity(member, 'standard') };
      index.records.push(record);
      identityIndexKeys(record.identity).forEach((key) => {
        addBucketValue(index.byIdentity, key, record);
      });
    }
    groupMemberIndexes.set(group, index);
    return index;
  };
  return {
    cachedIdentity,
    renderedByObject,
    renderedById,
    featuresByOrthogroupId,
    featureOrthogroupStatuses,
    sourceByObject,
    sourceByIdentity,
    memberIndexForGroup
  };
};

/** @type {WeakMap<object, WeakMap<PopupFeature[], PayloadFeatureIndex>>} */
const payloadFeatureIndexes = new WeakMap();

/**
 * The feature index of a lookup and a source list. The catalogs are replaced as
 * a whole, never edited in place, so their identities key the index.
 * @param {FeatureLookup} featureLookup
 * @param {PopupFeature[]} sourceFeatures
 * @returns {PayloadFeatureIndex}
 */
const payloadFeatureIndexFor = (featureLookup, sourceFeatures) => {
  if (!featureLookup || typeof featureLookup !== 'object' || !Array.isArray(sourceFeatures)) {
    return buildPayloadFeatureIndex(featureLookup, sourceFeatures);
  }
  let bySource = payloadFeatureIndexes.get(featureLookup);
  if (!bySource) {
    bySource = new WeakMap();
    payloadFeatureIndexes.set(featureLookup, bySource);
  }
  let index = bySource.get(sourceFeatures);
  if (!index) {
    index = buildPayloadFeatureIndex(featureLookup, sourceFeatures);
    bySource.set(sourceFeatures, index);
  }
  return index;
};

/**
 * @param {{
 *   featureLookup?: FeatureLookup,
 *   sourceFeatures?: PopupFeature[],
 *   orthogroups?: PopupGroup[],
 *   descriptor?: PairwiseMatchDescriptor | null
 * }} [options]
 * @returns {PairwisePayloadContext}
 */
const createPairwisePayloadContext = ({
  featureLookup = new Map(),
  sourceFeatures = [],
  orthogroups = [],
  descriptor = null
} = {}) => {
  const {
    cachedIdentity,
    renderedByObject,
    renderedById,
    featuresByOrthogroupId,
    featureOrthogroupStatuses,
    sourceByObject,
    sourceByIdentity,
    memberIndexForGroup
  } = payloadFeatureIndexFor(featureLookup, sourceFeatures);

  /** @type {Set<string>} */
  const rejectedRenderedIds = new Set();
  /** @param {unknown} renderedSvgId */
  const rawRenderedRecordForId = (renderedSvgId) => {
    const candidates = renderedById.get(normalizeText(renderedSvgId)) || [];
    return candidates.length === 1 ? candidates[0] : null;
  };
  for (const role of /** @type {const} */ (['query', 'subject'])) {
    const endpoint = descriptor?.endpoints?.[role];
    if (!endpoint?.constrained) continue;
    const renderedIds = splitMetadataValues(descriptor?.[`${role}FeatureSvgId`]);
    if (!endpoint.valid) {
      renderedIds.forEach((renderedId) => rejectedRenderedIds.add(renderedId));
      continue;
    }
    endpoint.references.forEach((reference) => {
      const endpointRecord = rawRenderedRecordForId(reference.renderedId.value);
      if (!endpointRecord || !identityMatches(
        reference,
        endpointRecord.identity,
        { includeRendered: true }
      )) rejectedRenderedIds.add(reference.renderedId.value);
    });
  }

  /** @param {unknown} renderedSvgId */
  const renderedRecordForId = (renderedSvgId) => {
    const id = normalizeText(renderedSvgId);
    if (!id || rejectedRenderedIds.has(id)) return null;
    return rawRenderedRecordForId(id);
  };
  /**
   * @param {unknown} renderedSvgId
   * @returns {PopupFeature | null}
   */
  const featureForLookupId = (renderedSvgId) => {
    const id = normalizeText(renderedSvgId);
    if (!id || rejectedRenderedIds.has(id)) return null;
    // Only a Map answers `get`; an array or plain object lookup yields no feature.
    const feature = /** @type {Partial<Map<string, PopupFeature>> | null} */ (featureLookup)?.get?.(id) || null;
    const record = feature ? renderedByObject.get(feature) : null;
    return record && rejectedRenderedIds.has(record.identity.renderedId.value)
      ? null
      : feature;
  };

  /**
   * @param {PopupGroup | null | undefined} group
   * @param {unknown} featureSvgId
   * @returns {PopupMember | null}
   */
  const groupMemberForFeatureId = (group, featureSvgId) => {
    const endpoint = renderedRecordForId(featureSvgId);
    if (!endpoint) return null;
    const memberIndex = memberIndexForGroup(group);
    /** @type {Set<MemberRecord>} */
    const seen = new Set();
    /** @type {MemberRecord[]} */
    const candidates = [];
    identityIndexKeys(endpoint.identity).forEach((key) => {
      for (const record of memberIndex.byIdentity.get(key) || []) {
        if (seen.has(record)) continue;
        seen.add(record);
        candidates.push(record);
      }
    });
    const matches = candidates.filter(({ identity }) => (
      identityMatches(identity, endpoint.identity, { includeRendered: true })
    ));
    return matches.length === 1 ? matches[0].member : null;
  };
  /**
   * @param {PopupGroup | null | undefined} group
   * @param {unknown} featureSvgIds
   */
  const groupHasFeatureIds = (group, featureSvgIds) => (
    splitMetadataValues(featureSvgIds)
      .some((featureSvgId) => Boolean(groupMemberForFeatureId(group, featureSvgId)))
  );

  const groupIndex = orthogroupIndexFor(orthogroups);
  const implicitGroupMatches = descriptor
    && !normalizeText(descriptor.orthogroupId)
    && descriptor.matchKind !== 'collinear'
    ? groupIndex.identified.filter((group) => (
      groupHasFeatureIds(group, descriptor.queryFeatureSvgId)
      || groupHasFeatureIds(group, descriptor.subjectFeatureSvgId)
    ))
    : [];

  /**
   * @param {unknown} orthogroupId
   * @param {unknown} [groupScope]
   * @returns {PopupGroup | null}
   */
  const orthogroupById = (orthogroupId, groupScope = '') => orthogroupInIndex(groupIndex, orthogroupId, groupScope);
  const orthogroupForMatch = () => {
    const direct = orthogroupById(descriptor?.orthogroupId, descriptor?.groupScope);
    if (direct) return direct;
    if (normalizeText(descriptor?.orthogroupId)) return null;
    return implicitGroupMatches.length === 1 ? implicitGroupMatches[0] : null;
  };

  /**
   * @param {PopupMember} member
   * @param {PopupFeature | null} feature
   * @returns {PopupFeature | null}
   */
  const renderedFeatureForMember = (member, feature) => {
    if (!feature) return null;
    const memberIdentity = cachedIdentity(member, 'standard');
    const sourceRecord = sourceByObject.get(feature);
    const sourceIdentity = sourceRecord?.sourceIdentity;
    if (
      !memberIdentity.usable
      || !sourceIdentity?.usable
      || !identityMatches(memberIdentity, sourceIdentity)
    ) return null;
    // `sourceIdentity?.usable` held above, so the source record exists.
    const memberRendered = memberIdentity.renderedId;
    const sourceRendered = /** @type {SourceRecord} */ (sourceRecord).sourceRendered;
    if (!memberRendered.valid || !sourceRendered.valid) return null;
    if (
      memberRendered.supplied
      && sourceRendered.supplied
      && memberRendered.value !== sourceRendered.value
    ) return null;
    let renderedSvgId = memberRendered.value || sourceRendered.value;
    if (!renderedSvgId) {
      if (!/** @type {SourceRecord} */ (sourceRecord).legacyRendered.valid) return null;
      renderedSvgId = /** @type {SourceRecord} */ (sourceRecord).legacyRendered.value;
    }
    const endpoint = renderedRecordForId(renderedSvgId);
    if (
      !endpoint
      || !identityMatches(memberIdentity, endpoint.identity, { includeRendered: true })
      || !identityMatches(sourceIdentity, endpoint.identity)
    ) return null;
    return endpoint.feature;
  };
  /**
   * @param {PopupMember} member
   * @returns {PopupFeature | null}
   */
  const featureForMember = (member) => {
    const memberIdentity = cachedIdentity(member, 'standard');
    const candidateKey = identityIndexKeys(memberIdentity)[0];
    if (!candidateKey) return null;
    const matches = (sourceByIdentity.get(candidateKey) || []).filter((record) => {
      if (!identityMatches(memberIdentity, record.candidateIdentity)) return false;
      // `identityMatches` is false for a missing or unusable candidate identity.
      if (!memberIdentity.renderedId.supplied && !(/** @type {FeatureIdentity} */ (record.candidateIdentity)).renderedId.supplied) {
        return true;
      }
      return Boolean(renderedFeatureForMember(member, record.feature));
    });
    return matches.length === 1 ? matches[0].feature : null;
  };

  /** @param {PopupFeature} feature */
  const identityForFeature = (feature) => {
    const standard = cachedIdentity(feature, 'standard');
    return standard.usable ? standard : cachedIdentity(feature, 'legacy');
  };
  /** @param {PopupFeature} feature */
  const orthogroupStatusForFeature = (feature) => (
    featureOrthogroupStatuses.get(feature) || featureOrthogroupIdStatus(feature)
  );
  /**
   * @param {unknown} orthogroupId
   * @param {PopupFeature | null} queryFeature
   * @param {PopupFeature | null} subjectFeature
   * @returns {PopupFeature[]}
   */
  const fallbackFeaturesForOrthogroup = (orthogroupId, queryFeature, subjectFeature) => {
    const id = normalizeText(orthogroupId);
    const matching = featuresByOrthogroupId.get(id) || [];
    if (matching.length) return matching;
    return [queryFeature, subjectFeature].filter(/** @param {PopupFeature | null} feature @returns {feature is PopupFeature} */ (feature) => {
      if (!feature) return false;
      const status = orthogroupStatusForFeature(feature);
      return status.valid && (!status.supplied || status.value === id);
    });
  };

  return {
    fallbackFeaturesForOrthogroup,
    featureForLookupId,
    featureForMember,
    groupHasFeatureIds,
    groupMemberForFeatureId,
    /** @param {PopupGroup | null | undefined} group */
    groupMemberRecords: (group) => memberIndexForGroup(group).records,
    identityForFeature,
    orthogroupById,
    orthogroupForMatch,
    renderedFeatureForMember
  };
};

/**
 * @param {{
 *   orthogroupId: unknown,
 *   queryFeature: PopupFeature | null,
 *   subjectFeature: PopupFeature | null,
 *   context: PairwisePayloadContext
 * }} input
 * @returns {PopupGroup | null}
 */
const buildFallbackOrthogroupWithContext = ({
  orthogroupId,
  queryFeature,
  subjectFeature,
  context
}) => {
  const id = normalizeText(orthogroupId);
  if (!id) return null;
  const fallbackFeatures = context.fallbackFeaturesForOrthogroup(
    id,
    queryFeature,
    subjectFeature
  );
  if (!fallbackFeatures.length) return null;
  const members = fallbackFeatures
    .map((feature) => featureToOrthogroupMember(
      feature,
      context.identityForFeature(feature)
    ))
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

/**
 * @param {OrthogroupTextOverrides | undefined} overrides
 * @param {unknown} key
 */
const overrideValue = (overrides, key) => {
  const normalizedKey = normalizeText(key);
  if (!normalizedKey || !overrides) return '';
  if (overrides instanceof Map) return normalizeText(overrides.get(normalizedKey));
  // Own entries only: `constructor` names an inherited member of a plain object, not an override.
  return Object.hasOwn(overrides, normalizedKey) ? normalizeText(overrides[normalizedKey]) : '';
};

/**
 * @param {string} title
 * @param {PopupRow[]} rows
 * @param {Partial<PopupSection>} [extras]
 * @returns {PopupSection}
 */
const section = (title, rows, extras = {}) => ({
  title,
  rows: rows.filter((row) => normalizeText(row.value)),
  ...extras
});

/**
 * @param {unknown} orthogroupId
 * @param {unknown} displayName
 */
const orthogroupTitle = (orthogroupId, displayName) => {
  const id = normalizeText(orthogroupId);
  const name = normalizeText(displayName);
  if (id && name) return `${id}:${name}`;
  return id || name || MATCH_KIND_TITLES.orthogroup;
};

/** @param {PopupMember | null | undefined} member */
const memberLocationText = (member) => {
  if (!member || typeof member !== 'object') return '';
  const range = locationRangeText(member);
  const strand = normalizeMemberStrand(member.strand);
  return range && strand ? `${range} (${strand})` : range;
};

/**
 * @param {PopupMember | null | undefined} member
 * @param {PopupFeature | null | undefined} feature
 * @param {SequenceKind} sequenceKind
 */
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

/** @param {unknown} strand */
const normalizeMemberStrand = (strand) => {
  if (strand === -1 || String(strand).trim() === '-1') return '-';
  if (strand === 1 || String(strand).trim() === '1') return '+';
  return normalizeText(strand);
};

/**
 * @param {PopupMember | null | undefined} member
 * @param {PopupFeature | null | undefined} feature
 * @param {string} nucleotideSequence
 * @param {string} aminoAcidSequence
 */
const buildMemberFeaturePayload = (member, feature, nucleotideSequence, aminoAcidSequence) => {
  const sourceFeature = feature && typeof feature === 'object' ? feature : {};
  const displayProteinId = resolveDisplayProteinId(
    sourceFeature,
    member,
    memberLocationText(member)
  );
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
    organism: sourceFeature.organism || member?.organism,
    nucleotide_sequence: nucleotideSequence,
    amino_acid_sequence: aminoAcidSequence
  };
};

/**
 * @param {PopupMember | null | undefined} member
 * @param {PopupFeature | null | undefined} feature
 * @param {string} nucleotideSequence
 * @param {string} aminoAcidSequence
 */
const buildMemberFastas = (member, feature, nucleotideSequence, aminoAcidSequence) => {
  const fastaFeature = buildMemberFeaturePayload(member, feature, nucleotideSequence, aminoAcidSequence);
  const fastas = buildFeatureSequenceFastas(fastaFeature, {
    nucleotideSequence,
    aminoAcidSequence
  });
  return {
    ntFasta: fastas.nucleotideFasta ? `${fastas.nucleotideFasta}\n` : '',
    aaFasta: fastas.aminoAcidFasta ? `${fastas.aminoAcidFasta}\n` : ''
  };
};

/**
 * @param {unknown} orthogroupId
 * @param {unknown} displayName
 * @param {SequenceKind} sequenceKind
 */
const orthogroupSequenceFilename = (orthogroupId, displayName, sequenceKind) => {
  const id = normalizeText(orthogroupId) || 'orthogroup';
  const name = makeSafeFilename(displayName || id, id);
  return `${makeSafeFilename(`${id}_${name}_${sequenceKindLabel(sequenceKind)}`)}.${sequenceExtension(sequenceKind)}`;
};

/**
 * @param {PopupMember} member
 * @param {unknown} orthogroupId
 * @param {SequenceKind} sequenceKind
 */
const orthogroupMemberSequenceFilename = (member, orthogroupId, sequenceKind) => {
  const id = normalizeText(orthogroupId) || 'orthogroup';
  const memberId = firstText(
    resolveDisplayProteinId(null, member, ''),
    memberFeatureSvgId(member),
    'member'
  );
  return `${makeSafeFilename(`${id}_${memberId}_${sequenceKindLabel(sequenceKind)}`)}.${sequenceExtension(sequenceKind)}`;
};

// Python's member defaults (OrthogroupMember in analysis/protein_colinearity.py).
// The member table names a role or confidence only when it differs (UX-B02).
const DEFAULT_MEMBER_ROLE = 'anchor';
const DEFAULT_MEMBER_CONFIDENCE = 'high';

/**
 * @param {string} role
 * @param {string} confidence
 */
const memberRoleNote = (role, confidence) => [
  role && role !== DEFAULT_MEMBER_ROLE ? role.replaceAll('_', ' ') : '',
  confidence && confidence !== DEFAULT_MEMBER_CONFIDENCE ? `${confidence} confidence` : ''
].filter(Boolean).join(' · ');

/**
 * @param {PopupGroup | null | undefined} group
 * @param {PairwisePayloadContext} context
 * @param {unknown} orthogroupId
 */
const buildOrthogroupMemberRows = (group, context, orthogroupId) => (
  context.groupMemberRecords(group)
    .map(({ member }) => {
      const feature = context.featureForMember(member);
      const renderedFeature = context.renderedFeatureForMember(member, feature);
      const nucleotideSequence = memberSequence(member, feature, 'nt');
      const aminoAcidSequence = memberSequence(member, feature, 'aa');
      const fastas = buildMemberFastas(
        member,
        feature,
        nucleotideSequence,
        aminoAcidSequence
      );
      const role = firstText(member?.role, member?.memberRole, member?.member_role);
      const confidence = firstText(member?.confidence, member?.memberConfidence, member?.member_confidence);
      return {
        feature: renderedFeature || feature,
        canOpen: Boolean(renderedFeature),
        record: firstText(member?.recordId, member?.record_id),
        coordinates: memberLocationText(member),
        proteinId: resolveDisplayProteinId(null, member),
        role,
        confidence,
        roleNote: memberRoleNote(role, confidence),
        assignmentReason: firstText(member?.assignmentReason, member?.assignment_reason),
        productOrNote: firstText(member?.product, member?.note),
        ntFasta: fastas.ntFasta,
        aaFasta: fastas.aaFasta,
        ntFilename: orthogroupMemberSequenceFilename(member, orthogroupId, 'nt'),
        aaFilename: orthogroupMemberSequenceFilename(member, orthogroupId, 'aa')
      };
    })
    .filter((row) => row.record || row.coordinates || row.proteinId || row.role || row.confidence || row.assignmentReason || row.productOrNote || row.ntFasta || row.aaFasta)
);

/**
 * @param {{
 *   feature: PopupFeature | null,
 *   featureSvgIds: unknown,
 *   context: PairwisePayloadContext,
 *   group: PopupGroup | null,
 *   fallbackProteinIds: unknown,
 *   locusId: unknown,
 *   displayName: unknown
 * }} input
 */
const resolveFeatureSectionProteinIds = ({
  feature,
  featureSvgIds,
  context,
  group,
  fallbackProteinIds,
  locusId,
  displayName
}) => {
  /** @type {string[]} */
  const values = [];
  /**
   * @param {PopupFeature | null} candidateFeature
   * @param {PopupMember | null} [member]
   */
  const addFeatureProteinId = (candidateFeature, member = null) => {
    const text = resolveDisplayProteinId(candidateFeature, member, '');
    addUniqueDisplayText(values, text);
  };

  const ids = splitMetadataValues(featureSvgIds);
  ids.forEach((featureSvgId) => {
    const candidateFeature = context.featureForLookupId(featureSvgId);
    const member = context.groupMemberForFeatureId(group, featureSvgId);
    addFeatureProteinId(candidateFeature, member);
  });
  if (feature) {
    addFeatureProteinId(
      feature,
      ids.length === 1 ? context.groupMemberForFeatureId(group, ids[0]) : null
    );
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

/** @param {...unknown} values */
const firstNonInternalDisplayText = (...values) => {
  for (const value of values) {
    const text = normalizeText(value);
    if (text && !isInternalDisplayId(text)) return text;
  }
  return '';
};

/**
 * @param {PopupFeature | null | undefined} feature
 * @param {unknown} fallbackDisplayName
 */
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

/**
 * @param {PopupFeature | null | undefined} feature
 * @param {unknown} fallbackLocusId
 */
const featureRowLocusId = (feature, fallbackLocusId) => firstText(
  feature?.locus_tag,
  feature?.locusTag,
  qualifierFirstValue(feature, 'locus_tag'),
  feature?.gene_id,
  feature?.geneId,
  qualifierFirstValue(feature, 'gene_id'),
  fallbackLocusId
);

/**
 * @param {unknown} primaryLabel
 * @param {...unknown} values
 */
const featureRowSubLabel = (primaryLabel, ...values) => {
  const seen = new Set(splitMetadataValues(primaryLabel));
  /** @type {string[]} */
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

/**
 * @param {{
 *   featureSvgIds: unknown,
 *   context: PairwisePayloadContext,
 *   group: PopupGroup | null,
 *   recordId: unknown,
 *   interval: unknown,
 *   proteinId: unknown,
 *   locusId: unknown,
 *   displayName: unknown
 * }} input
 */
const buildFeatureListRows = ({
  featureSvgIds,
  context,
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
    const feature = svgId ? context.featureForLookupId(svgId) : null;
    const member = svgId ? context.groupMemberForFeatureId(group, svgId) : null;
    const fallbackProteinId = firstNonInternalDisplayText(
      locusIds[index],
      displayNames[index],
      proteinIds[index]
    );
    const resolvedProteinId = resolveDisplayProteinId(feature, member, '');
    const displayProteinId = firstText(
      isInternalDisplayId(resolvedProteinId) ? '' : resolvedProteinId,
      fallbackProteinId,
      isInternalDisplayId(proteinIds[index]) ? '' : proteinIds[index]
    );
    const rowRecord = firstText(feature?.record_id, feature?.recordId, member?.recordId, member?.record_id, recordId);
    const rowLocation = firstText(
      formatFeatureLocation(feature),
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
      canOpen: Boolean(
        featureRenderedSvgId(feature) && (!group || member)
      ),
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

/** @param {MemberRow[]} memberRows */
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

/**
 * @param {MemberRow[]} memberRows
 * @param {SequenceKind} sequenceKind
 */
const orthogroupMemberSequenceText = (memberRows, sequenceKind) => {
  const key = sequenceKindLabel(sequenceKind) === 'aa' ? 'aaFasta' : 'ntFasta';
  return (Array.isArray(memberRows) ? memberRows : [])
    .map((row) => String(row?.[key] ?? ''))
    .filter((text) => text.trim())
    .join('');
};

/**
 * @param {MemberRow[]} memberRows
 * @param {unknown} orthogroupId
 * @param {unknown} displayName
 */
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

/**
 * "N genes in M records" (UX-B02): the one count line of a group.
 * @param {unknown} memberCount
 * @param {unknown} recordCoverage
 */
export const groupCountText = (memberCount, recordCoverage) => {
  const plural = (/** @type {unknown} */ value, /** @type {string} */ noun) => {
    const count = Number(normalizeText(value));
    return Number.isInteger(count) && count > 0 ? `${count} ${noun}${count === 1 ? '' : 's'}` : '';
  };
  const genes = plural(memberCount, 'gene');
  const records = plural(recordCoverage, 'record');
  return genes && records ? `${genes} in ${records}` : genes || records;
};

/**
 * The lookups of the payload context that the group and gene-pair views use.
 * @typedef {{
 *   orthogroupById: (orthogroupId: string, groupScope?: string) => PopupGroup | null,
 *   groupMemberForFeatureId: (group: PopupGroup | null, featureSvgId: string) => PopupMember | null
 * }} GroupLookupPort
 */

/**
 * The view of one similarity group in the match popup (UX-B02): the name with
 * the ID as a badge, one count line, the member table and its copy menu, and
 * the computed facts behind "Technical details". A description the user wrote
 * is shown under the name; a computed description is a technical detail.
 * @param {{
 *   orthogroupId: string,
 *   group: PopupGroup | null,
 *   context: PairwisePayloadContext,
 *   displayName: string,
 *   userDescription: string,
 *   groupScope: string,
 *   localGroups?: boolean
 * }} input
 */
const buildGroupView = ({
  orthogroupId,
  group,
  context,
  displayName,
  userDescription,
  groupScope,
  localGroups = false
}) => {
  const memberRows = buildOrthogroupMemberRows(group, context, orthogroupId);
  const memberCount = firstText(
    group?.member_count,
    group?.memberCount,
    Array.isArray(group?.members) ? group.members.length : ''
  );
  const recordCoverage = firstText(group?.record_coverage_count, group?.recordCoverage);
  /** @type {PopupRow[]} */
  const technicalRows = [];
  addRow(technicalRows, localGroups ? 'Collinear group ID' : 'Similarity group ID', orthogroupId);
  addRow(technicalRows, 'Description', userDescription ? '' : group?.description);
  addRow(technicalRows, 'Scope', groupMetadataScopeLabel(normalizeGroupMetadataScope(groupScopeValue(group) || groupScope)));
  addRow(technicalRows, 'RBH seeds', Array.isArray(group?.rbhOrthogroupIds) ? group.rbhOrthogroupIds.join('; ') : '');
  addRow(technicalRows, 'Group paths', firstText(
    group?.orthologPathCount,
    Array.isArray(group?.orthologPaths) ? String(group.orthologPaths.length) : ''
  ));
  addRow(technicalRows, 'Related edges', firstText(
    group?.relatedEdgeCount,
    Array.isArray(group?.relatedEdges) ? String(group.relatedEdges.length) : ''
  ));
  const id = normalizeText(orthogroupId);
  const name = normalizeText(displayName);
  return {
    id,
    displayName: name,
    title: name || id || MATCH_KIND_TITLES.orthogroup,
    badge: name && id && name !== id ? id : '',
    description: normalizeText(userDescription),
    memberCount,
    recordCoverage,
    countText: groupCountText(memberCount || memberRows.length, recordCoverage),
    technicalRows,
    ...orthogroupMemberSectionExtras(memberRows, orthogroupId, displayName)
  };
};

/**
 * @param {{
 *   orthogroupIds: string[],
 *   context: PairwisePayloadContext,
 *   orthogroupNameOverrides: OrthogroupTextOverrides | undefined,
 *   orthogroupDescriptionOverrides: OrthogroupTextOverrides | undefined,
 *   groupScope: string
 * }} input
 */
const buildBlockOrthogroups = ({
  orthogroupIds,
  context,
  orthogroupNameOverrides,
  orthogroupDescriptionOverrides,
  groupScope
}) => orthogroupIds.map((orthogroupId) => {
  const group = context.orthogroupById(orthogroupId, groupScope);
  return buildGroupView({
    orthogroupId,
    group,
    context,
    displayName: firstText(
      overrideValue(orthogroupNameOverrides, orthogroupId),
      group?.displayName,
      group?.display_name,
      group?.name
    ),
    userDescription: overrideValue(orthogroupDescriptionOverrides, orthogroupId),
    groupScope,
    localGroups: normalizeGroupMetadataScope(groupScopeValue(group) || groupScope) === 'adjacent_local'
  });
});

// Collinear block popup (UX-B01): one summary line, the gene pairs, and the
// block facts behind "Details".
const COLLINEAR_ORIENTATION_TEXT = Object.freeze({ plus: 'Same orientation', minus: 'Inverted' });
const COLLINEAR_BLOCK_NOTE = 'Block envelopes may include intergenic sequence and genes that are not anchors.';

/** @param {number} value */
const groupedInteger = (value) => value.toLocaleString('en-US');

/** @param {unknown} value */
const finiteNumber = (value) => {
  const text = normalizeText(value);
  const number = Number(text);
  return text && Number.isFinite(number) ? number : null;
};

/**
 * "AP027133.1 10,124..28,017 (17.9 kb)" for one side of a block.
 * @param {string} recordId
 * @param {string} interval
 * @param {string} tableInterval
 */
const blockSpanView = (recordId, interval, tableInterval) => {
  const match = /^(\d+)\.\.(\d+)$/.exec(normalizeText(interval));
  const start = match ? Number(match[1]) : 0;
  const end = match ? Number(match[2]) : 0;
  const length = match ? Math.abs(end - start) + 1 : 0;
  const lengthText = length >= 1000 ? `${(length / 1000).toFixed(1)} kb` : (length ? `${length} bp` : '');
  const intervalText = match ? `${groupedInteger(start)}..${groupedInteger(end)}` : normalizeText(interval);
  return {
    recordId: normalizeText(recordId),
    interval: normalizeText(interval),
    intervalText,
    lengthText,
    tableInterval: normalizeText(tableInterval),
    text: [normalizeText(recordId), intervalText, lengthText ? `(${lengthText})` : ''].filter(Boolean).join(' ')
  };
};

/**
 * @param {Element} element
 * @param {PairwiseMatchDescriptor} descriptor
 */
const collinearDetailRows = (element, descriptor) => {
  /** @type {PopupRow[]} */
  const rows = [];
  const identity = finiteNumber(descriptor.identity);
  const alignedLength = finiteNumber(attr(element, 'data-alignment-length'));
  addRow(rows, 'Block ID', descriptor.collinearityBlockId);
  addRow(rows, 'Kind', attr(element, 'data-collinearity-block-kind'));
  addRow(rows, 'Orientation', attr(element, 'data-collinearity-orientation'));
  addRow(rows, 'Color mode', attr(element, 'data-collinearity-color-mode'));
  if (descriptor.matchKind === 'collinear') {
    addRow(rows, 'Average identity', identity === null ? descriptor.identity : `${identity.toFixed(2)}%`);
    // LOSATP anchors: the summed protein alignment length.
    addRow(rows, 'Aligned length', alignedLength === null ? '' : `${groupedInteger(alignedLength)} aa`);
  }
  addRow(rows, 'Block score', attr(element, 'data-collinearity-block-score'));
  addRow(rows, 'Block e-value', attr(element, 'data-collinearity-block-evalue'));
  addRow(rows, 'Anchor', [
    attr(element, 'data-collinearity-anchor-index'),
    attr(element, 'data-collinearity-anchor-count')
  ].filter(Boolean).join(' / '));
  return rows;
};

/**
 * One row per gene pair (query gene | subject gene | group), in the order of
 * the block's groups; a gene in no group of the block gets a row of its own.
 * @param {{
 *   queryRows: FeatureListRow[],
 *   subjectRows: FeatureListRow[],
 *   blockOrthogroups: GroupView[],
 *   context: GroupLookupPort,
 *   groupScope: string
 * }} input
 */
const buildGenePairs = ({ queryRows, subjectRows, blockOrthogroups, context, groupScope }) => {
  /** @type {Set<FeatureListRow>} */
  const used = new Set();
  /**
   * @param {FeatureListRow[]} rows
   * @param {PopupGroup} group
   */
  const membersOf = (rows, group) => rows.filter((row) => (
    !used.has(row) && row.svgId && context.groupMemberForFeatureId(group, String(row.svgId))
  ));
  /** @type {{ key: string, query: FeatureListRow | null, subject: FeatureListRow | null, group: { id: string, displayName: string } | null }[]} */
  const pairs = [];
  blockOrthogroups.forEach((view) => {
    const group = context.orthogroupById(view.id, groupScope);
    const queries = group ? membersOf(queryRows, group) : [];
    const subjects = group ? membersOf(subjectRows, group) : [];
    [...queries, ...subjects].forEach((row) => used.add(row));
    const count = Math.max(queries.length, subjects.length, 1);
    for (let index = 0; index < count; index += 1) {
      pairs.push({
        key: `${view.id}-${index}`,
        query: queries[index] || null,
        subject: subjects[index] || null,
        group: { id: view.id, displayName: view.displayName }
      });
    }
  });
  [...queryRows, ...subjectRows].filter((row) => !used.has(row)).forEach((row) => {
    const isQuery = queryRows.includes(row);
    pairs.push({
      key: `${isQuery ? 'query' : 'subject'}-${row.key}`,
      query: isQuery ? row : null,
      subject: isQuery ? null : row,
      group: null
    });
  });
  return pairs;
};

/** @param {FeatureListRow | null} row */
const genePairCells = (row) => (row
  ? [row.record, row.location, row.proteinId, row.product]
  : ['', '', '', '']
).map((value) => normalizeText(value));

/**
 * @param {{
 *   element: Element,
 *   descriptor: PairwiseMatchDescriptor,
 *   genePairs: ReturnType<typeof buildGenePairs>,
 *   localGroups: boolean
 * }} input
 */
const buildCollinearBlockView = ({ element, descriptor, genePairs, localGroups }) => {
  const query = blockSpanView(descriptor.queryRecordId, descriptor.queryInterval, descriptor.queryTableInterval);
  const subject = blockSpanView(descriptor.subjectRecordId, descriptor.subjectInterval, descriptor.subjectTableInterval);
  const orientation = normalizeText(descriptor.orientation).toLowerCase();
  // Own keys only: `constructor` and `__proto__` name inherited members, not orientations (OV-322).
  const orientationText = Object.hasOwn(COLLINEAR_ORIENTATION_TEXT, orientation)
    ? COLLINEAR_ORIENTATION_TEXT[/** @type {keyof typeof COLLINEAR_ORIENTATION_TEXT} */ (orientation)]
    : orientation;
  const identity = finiteNumber(descriptor.identity);
  const identityText = identity === null ? '' : `${Math.round(identity)}% identity`;
  const pairCount = genePairs.filter((pair) => pair.query && pair.subject).length;
  const pairsLabel = `${pairCount} gene pair${pairCount === 1 ? '' : 's'}`;
  const facts = [orientationText, identityText, pairsLabel].filter(Boolean);
  const groupHeading = localGroups ? 'Collinear group' : 'Similarity group';
  return {
    query,
    subject,
    facts,
    identity: descriptor.identity,
    pairCount,
    pairsLabel,
    groupHeading,
    genePairs,
    summaryText: `Collinear block: ${query.text} vs ${subject.text}; ${facts.join('; ')}`,
    pairsCopyText: [
      ['Query record', 'Query location', 'Query protein ID', 'Query product',
        'Subject record', 'Subject location', 'Subject protein ID', 'Subject product',
        `${groupHeading} ID`, `${groupHeading} name`].join('\t'),
      ...genePairs.map((pair) => [
        ...genePairCells(pair.query),
        ...genePairCells(pair.subject),
        pair.group?.id || '',
        pair.group?.displayName || ''
      ].join('\t'))
    ].join('\n'),
    detailRows: collinearDetailRows(element, descriptor),
    note: COLLINEAR_BLOCK_NOTE
  };
};

/**
 * @param {{
 *   title: string,
 *   feature: PopupFeature | null,
 *   recordId: unknown,
 *   interval: unknown,
 *   proteinId: unknown,
 *   locusId: unknown,
 *   displayName: unknown,
 *   featureSvgIds: unknown,
 *   context: PairwisePayloadContext,
 *   group: PopupGroup | null
 * }} input
 */
const buildFeatureRows = ({
  title,
  feature,
  recordId,
  interval,
  proteinId,
  locusId,
  displayName,
  featureSvgIds,
  context,
  group
}) => {
  /** @type {PopupRow[]} */
  const rows = [];
  const featureRows = buildFeatureListRows({
    featureSvgIds,
    context,
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
    context,
    group,
    fallbackProteinIds: proteinId,
    locusId,
    displayName
  });
  addRow(rows, 'Record', firstText(feature?.record_id, recordId));
  addRow(rows, 'Location', firstText(formatFeatureLocation(feature), interval));
  addRow(rows, displayProteinIds.includes(';') ? 'Protein IDs' : 'Protein ID', displayProteinIds);
  addRow(rows, 'Gene', firstText(feature?.gene, qualifierFirstValue(feature, 'gene')));
  addRow(rows, 'Locus tag', firstText(feature?.locus_tag, feature?.locusTag, qualifierFirstValue(feature, 'locus_tag')));
  addRow(rows, 'Product', firstText(feature?.product, qualifierFirstValue(feature, 'product')));
  addRow(rows, 'Locus ID', locusId);
  addRow(rows, 'Display name', displayName);
  return section(title, rows, { featureRows });
};

/**
 * @param {Element} element
 * @param {MatchKind} matchKind
 */
const buildMatchSpans = (element, matchKind) => {
  if (matchKind === 'orthogroup') return [];
  const queryRecordIndex = integerAttrStatus(element, 'data-query-record-index');
  const subjectRecordIndex = integerAttrStatus(element, 'data-subject-record-index');
  const sourceIndex = integerAttrStatus(element, 'data-source-index');
  const referenceSide = attr(element, 'data-reference-side').toLowerCase();
  /** @param {'query' | 'subject'} role */
  const homologyDisplayRole = (role) => role === referenceSide ? 'Reference' : 'Comparison';
  return [
    {
      role: 'query',
      sourceKey: matchKind === 'homology' ? '' : (
        queryRecordIndex.value !== null ? `linear:record:${queryRecordIndex.value}` : ''
      ),
      recordId: firstText(attr(element, 'data-query-record-id'), attr(element, 'data-query')),
      start: attr(element, 'data-qstart'),
      end: attr(element, 'data-qend'),
      sourceSpan: matchKind === 'homology' ? null : readRecordSourceSpan(element, queryRecordIndex.value),
      displayRole: matchKind === 'homology' ? homologyDisplayRole('query') : 'Query',
      sourceIndex: sourceIndex.value,
      recordIndex: queryRecordIndex.value,
      identityValid: queryRecordIndex.valid && sourceIndex.valid,
      referenceSide
    },
    {
      role: 'subject',
      sourceKey: matchKind === 'homology' ? '' : (
        subjectRecordIndex.value !== null ? `linear:record:${subjectRecordIndex.value}` : ''
      ),
      recordId: firstText(attr(element, 'data-subject-record-id'), attr(element, 'data-subject')),
      start: attr(element, 'data-sstart'),
      end: attr(element, 'data-send'),
      sourceSpan: matchKind === 'homology' ? null : readRecordSourceSpan(element, subjectRecordIndex.value),
      displayRole: matchKind === 'homology' ? homologyDisplayRole('subject') : 'Subject',
      sourceIndex: sourceIndex.value,
      recordIndex: subjectRecordIndex.value,
      identityValid: subjectRecordIndex.valid && sourceIndex.valid,
      referenceSide
    }
  ];
};

/**
 * @param {Element} element
 * @param {MatchKind} matchKind
 * @param {string} matchId
 * @param {SequenceSourceResolver | null} resolveSequenceSource
 */
const buildSequenceBundleForMatch = (element, matchKind, matchId, resolveSequenceSource) => {
  const spans = buildMatchSpans(element, matchKind);
  if (!spans.length) return null;
  const bundle = buildMatchSequenceBundle(spans, {
    matchId,
    resolveSequenceSource: (sourceKey, recordId, context) => {
      /** @type {SequenceSourceResolution | (SequenceSource & { source?: undefined }) | null} */
      const resolved = typeof resolveSequenceSource === 'function'
        ? resolveSequenceSource(sourceKey, recordId, context)
        : null;
      const resolvedSource = resolved?.source !== undefined ? resolved.source : resolved;
      if (matchKind === 'homology' && context?.origin === 'homology-comparison' && !resolvedSource) {
        return {
          source: null,
          reason: 'Comparison sequence was not supplied for this BLAST source.'
        };
      }
      return resolved;
    },
    contextForSpan: (span) => {
      if (matchKind !== 'homology') {
        return { origin: 'linear-record', recordIndex: span.recordIndex };
      }
      const isReference = span.role === span.referenceSide;
      return isReference
        ? { origin: 'circular-reference', recordIndex: span.recordIndex }
        : { origin: 'homology-comparison', sourceIndex: span.sourceIndex };
    },
    unavailableReasonForSpan: (span) => (
      span.identityValid ? '' : 'The sequence identity metadata for this match is invalid.'
    )
  });
  // The span rows show input-file coordinates, as the Interval rows and the
  // FASTA headers do (R5, CI-05); the bases are read on the displayed record.
  bundle.entries = bundle.entries.map((entry, index) => {
    const source = recordSourceInterval(spans[index].sourceSpan, entry.span.start, entry.span.end);
    return source ? { ...entry, span: { ...entry.span, start: source.start, end: source.end } } : entry;
  });
  return bundle;
};

/**
 * @param {SVGElement} element
 * @param {MatchPopupPayloadOptions} [options]
 * @returns {MatchPopupPayload | null}
 */
export const buildMatchPopupPayload = (
  element,
  {
    featureLookup = new Map(),
    sourceFeatures = [],
    orthogroups = [],
    orthogroupNameOverrides = null,
    orthogroupDescriptionOverrides = null,
    resolveSequenceSource = null
  } = {}
) => {
  if (!element) return null;
  const descriptor = readPairwiseMatchDescriptor(element);
  const {
    matchKind,
    matchId,
    orthogroupId,
    orthogroupIds,
    collinearityBlockId,
    groupScope,
    normalizedGroupScope,
    queryFeatureSvgId,
    subjectFeatureSvgId,
    queryInterval: qInterval,
    subjectInterval: sInterval,
    title,
    subtitle
  } = descriptor;
  const context = createPairwisePayloadContext({
    featureLookup,
    sourceFeatures,
    orthogroups,
    descriptor
  });
  const queryFeature = context.featureForLookupId(queryFeatureSvgId);
  const subjectFeature = context.featureForLookupId(subjectFeatureSvgId);
  const group = matchKind === 'collinear'
    ? null
    : context.orthogroupForMatch() || buildFallbackOrthogroupWithContext({
      orthogroupId,
      queryFeature,
      subjectFeature,
      context
    });
  const userDescription = overrideValue(orthogroupDescriptionOverrides, orthogroupId);
  const groupView = matchKind !== 'collinear' && (normalizeText(orthogroupId) || group)
    ? buildGroupView({
      orthogroupId: firstText(orthogroupId, group ? orthogroupIdStatus(group).value : ''),
      group,
      context,
      displayName: firstText(
        overrideValue(orthogroupNameOverrides, orthogroupId),
        group?.displayName,
        group?.display_name,
        group?.name,
        featureProduct(queryFeature),
        featureProduct(subjectFeature)
      ),
      userDescription,
      groupScope
    })
    : null;
  const identity = {
    matchKind,
    orthogroupId,
    collinearityBlockId,
    queryFeatureSvgId,
    subjectFeatureSvgId,
    fill: descriptor.fill
  };

  if (matchKind === 'orthogroup') {
    return {
      id: firstText(matchId, orthogroupId),
      title: orthogroupTitle(orthogroupId, groupView?.displayName),
      subtitle: '',
      ...identity,
      group: groupView,
      sections: []
    };
  }

  const sequenceBundle = buildSequenceBundleForMatch(
    element,
    matchKind,
    matchId || 'match',
    resolveSequenceSource
  );

  /** @param {'query' | 'subject'} role */
  const featureListRows = (role) => buildFeatureListRows({
    featureSvgIds: role === 'query' ? queryFeatureSvgId : subjectFeatureSvgId,
    context,
    group,
    recordId: role === 'query' ? descriptor.queryRecordId : descriptor.subjectRecordId,
    interval: role === 'query' ? qInterval : sInterval,
    proteinId: attr(element, `data-${role}-protein-id`),
    locusId: attr(element, `data-${role}-locus-id`),
    displayName: attr(element, `data-${role}-display-name`)
  });

  if (matchKind === 'collinear') {
    const blockOrthogroups = buildBlockOrthogroups({
      orthogroupIds,
      context,
      orthogroupNameOverrides,
      orthogroupDescriptionOverrides,
      groupScope
    });
    return {
      id: matchId,
      title,
      subtitle,
      ...identity,
      groupScope: normalizedGroupScope,
      blockOrthogroupCount: blockOrthogroups.length || orthogroupIds.length,
      blockOrthogroups,
      block: buildCollinearBlockView({
        element,
        descriptor,
        genePairs: buildGenePairs({
          queryRows: featureListRows('query'),
          subjectRows: featureListRows('subject'),
          blockOrthogroups,
          context,
          groupScope
        }),
        localGroups: normalizedGroupScope === 'adjacent_local'
      }),
      sequenceBundle,
      sections: []
    };
  }

  /** @type {PopupRow[]} */
  const summaryRows = [];
  addRow(summaryRows, 'Query record', descriptor.queryRecordId);
  addRow(summaryRows, 'Subject record', descriptor.subjectRecordId);
  addRow(summaryRows, 'Query interval', qInterval);
  addRow(summaryRows, 'Query table interval', descriptor.queryTableInterval);
  addRow(summaryRows, 'Subject interval', sInterval);
  addRow(summaryRows, 'Subject table interval', descriptor.subjectTableInterval);
  addRow(summaryRows, 'Orientation', descriptor.orientation);
  if (matchKind === 'homology') {
    addRow(summaryRows, 'Ring label', attr(element, 'data-track-label'));
    const rawSourceIndex = integerAttr(element, 'data-source-index');
    addRow(summaryRows, 'Source index', rawSourceIndex !== null ? String(rawSourceIndex + 1) : '');
    addRow(summaryRows, 'Reference side', attr(element, 'data-reference-side'));
    addRow(summaryRows, 'Reference record', attr(element, 'data-reference-record-id'));
  }

  /** @type {PopupRow[]} */
  const alignmentRows = [];
  addRow(alignmentRows, 'Identity', descriptor.identity);
  addRow(alignmentRows, 'Alignment length', attr(element, 'data-alignment-length'));
  addRow(alignmentRows, 'E-value', attr(element, 'data-evalue'));
  addRow(alignmentRows, 'Bit score', attr(element, 'data-bitscore'));
  addRow(alignmentRows, 'Mismatches', attr(element, 'data-mismatches'));
  addRow(alignmentRows, 'Gap opens', attr(element, 'data-gap-opens'));
  addRow(alignmentRows, 'Edge kind', attr(element, 'data-edge-kind'));
  addRow(alignmentRows, 'Render role', attr(element, 'data-render-role'));
  addRow(alignmentRows, 'RBH seed', attr(element, 'data-rbh-orthogroup-id'));
  addRow(alignmentRows, 'Group path ID', attr(element, 'data-ortholog-path-id'));
  addRow(alignmentRows, 'Query group genes', attr(element, 'data-query-orthogroup-member-count'));
  addRow(alignmentRows, 'Subject group genes', attr(element, 'data-subject-orthogroup-member-count'));
  addRow(alignmentRows, 'Query member role', attr(element, 'data-query-orthogroup-role'));
  addRow(alignmentRows, 'Subject member role', attr(element, 'data-subject-orthogroup-role'));
  addRow(alignmentRows, 'Query confidence', attr(element, 'data-query-orthogroup-confidence'));
  addRow(alignmentRows, 'Subject confidence', attr(element, 'data-subject-orthogroup-confidence'));
  addRow(alignmentRows, 'Query assignment', attr(element, 'data-query-orthogroup-assignment-reason'));
  addRow(alignmentRows, 'Subject assignment', attr(element, 'data-subject-orthogroup-assignment-reason'));

  const sections = [section('Summary', summaryRows), section('Alignment', alignmentRows)];
  const blockRows = collinearDetailRows(element, descriptor);
  if (blockRows.length > 0) sections.push(section('Collinearity', blockRows));
  if (matchKind !== 'homology') {
    sections.push(buildFeatureRows({
      title: 'Query',
      feature: queryFeature,
      recordId: descriptor.queryRecordId,
      interval: qInterval,
      proteinId: attr(element, 'data-query-protein-id'),
      locusId: attr(element, 'data-query-locus-id'),
      displayName: attr(element, 'data-query-display-name'),
      featureSvgIds: queryFeatureSvgId,
      context,
      group
    }));
    sections.push(buildFeatureRows({
      title: 'Subject',
      feature: subjectFeature,
      recordId: descriptor.subjectRecordId,
      interval: sInterval,
      proteinId: attr(element, 'data-subject-protein-id'),
      locusId: attr(element, 'data-subject-locus-id'),
      displayName: attr(element, 'data-subject-display-name'),
      featureSvgIds: subjectFeatureSvgId,
      context,
      group
    }));
  }

  return {
    id: matchId,
    title,
    subtitle,
    ...identity,
    groupScope: normalizedGroupScope,
    group: groupView,
    sequenceTitle: 'Matched sequences',
    sequenceBundle,
    sections: sections.filter((entry) => (
      entry.rows.length > 0 ||
      (Array.isArray(entry.featureRows) && entry.featureRows.length > 0)
    ))
  };
};

/**
 * @param {PairwiseHoverRowValues} [values]
 */
const formatPairwiseMatchHoverRows = ({
  matchKind,
  identity,
  queryInterval,
  subjectInterval,
  orthogroupId,
  displayName,
  memberCount,
  groupScope,
  blockOrthogroupCount,
  collinearityBlockId
} = {}) => {
  /** @type {PopupRow[]} */
  const rows = [];
  /**
   * @param {string} label
   * @param {unknown} value
   */
  const addFirst = (label, value) => addRow(rows, label, value);
  if (!matchKind) return rows;
  if (matchKind === 'orthogroup') {
    addFirst('Kind', matchKind);
    addFirst('Similarity group', orthogroupId);
    addFirst('Display name', displayName);
    addFirst('Genes', memberCount);
    return rows.slice(0, 6);
  }
  addFirst('Kind', matchKind);
  addFirst('Identity', identity);
  addFirst('Query', queryInterval);
  addFirst('Subject', subjectInterval);
  if (matchKind === 'collinear') {
    addFirst(
      groupScope === 'adjacent_local'
        ? 'Collinear groups'
        : 'Similarity groups',
      String(blockOrthogroupCount ?? '')
    );
  } else {
    addFirst('Similarity group', orthogroupId);
  }
  addFirst('Block', collinearityBlockId);
  return rows.slice(0, 6);
};

/**
 * @param {SVGElement} element
 * @param {MatchHoverSummaryOptions} [options]
 */
export const buildPairwiseMatchHoverSummary = (
  element,
  {
    orthogroups = [],
    orthogroupNameOverrides = null
  } = {}
) => {
  if (!element) return null;
  const descriptor = readPairwiseMatchDescriptor(element);
  const groupValues = descriptor.matchKind === 'orthogroup'
    ? (typeof orthogroups === 'function' ? orthogroups() : orthogroups)
    : [];
  const group = descriptor.matchKind === 'orthogroup'
    ? orthogroupInIndex(
      orthogroupIndexFor(groupValues),
      descriptor.orthogroupId,
      descriptor.groupScope
    )
    : null;
  const displayName = descriptor.matchKind === 'orthogroup'
    ? firstText(
      overrideValue(orthogroupNameOverrides, descriptor.orthogroupId),
      group?.displayName,
      group?.display_name,
      group?.name
    )
    : '';
  const memberCount = descriptor.matchKind === 'orthogroup'
    ? firstText(group?.member_count, group?.memberCount)
    : '';
  return {
    id: descriptor.matchId,
    title: descriptor.matchKind === 'orthogroup'
      ? orthogroupTitle(descriptor.orthogroupId, displayName)
      : descriptor.title,
    subtitle: descriptor.matchKind === 'orthogroup' ? '' : descriptor.subtitle,
    fill: descriptor.fill,
    rows: formatPairwiseMatchHoverRows({
      matchKind: descriptor.matchKind,
      identity: descriptor.identity,
      queryInterval: descriptor.queryInterval,
      subjectInterval: descriptor.subjectInterval,
      orthogroupId: descriptor.orthogroupId,
      displayName,
      memberCount,
      groupScope: descriptor.normalizedGroupScope,
      blockOrthogroupCount: descriptor.orthogroupIds.length,
      collinearityBlockId: descriptor.collinearityBlockId
    })
  };
};

/** @param {MatchPopupPayload | null | undefined} payload */
export const buildPairwiseMatchHoverRows = (payload) => {
  if (!payload) return [];
  const summary = payload.sections.find((entry) => entry.title === 'Summary');
  const alignment = payload.sections.find((entry) => entry.title === 'Alignment');
  /**
   * @param {PopupSection | undefined} sectionEntry
   * @param {string} label
   */
  const findValue = (sectionEntry, label) => (
    sectionEntry?.rows.find((row) => row.label === label)?.value || ''
  );
  return formatPairwiseMatchHoverRows({
    matchKind: payload.matchKind,
    identity: payload.block ? payload.block.identity : findValue(alignment, 'Identity'),
    queryInterval: payload.block ? payload.block.query.interval : findValue(summary, 'Query interval'),
    subjectInterval: payload.block ? payload.block.subject.interval : findValue(summary, 'Subject interval'),
    orthogroupId: payload.orthogroupId,
    displayName: payload.matchKind === 'orthogroup' ? payload.group?.displayName : '',
    memberCount: payload.matchKind === 'orthogroup' ? payload.group?.memberCount : '',
    groupScope: payload.groupScope,
    blockOrthogroupCount: payload.blockOrthogroupCount,
    collinearityBlockId: payload.collinearityBlockId
  });
};
