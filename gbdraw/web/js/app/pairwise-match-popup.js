import {
  isInternalProteinDisplayId,
  resolveDisplayProteinId
} from './feature-utils.js';
import { buildFeatureSequenceFastas } from './feature-sequence-fasta.js';
import { buildMatchSequenceBundle } from './match-sequences.js';
import {
  groupMetadataScopeLabel,
  normalizeGroupMetadataScope
} from './losat-normalization.js';
import {
  RENDERED_FEATURE_ID_KEYS,
  STABLE_FEATURE_ID_KEYS,
  featureIdentity,
  identityMatches,
  renderedFeatureIdentity,
  textAliasStatus
} from '../services/feature-identity.js';

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
const normalizeText = (value) => String(value === null || value === undefined ? '' : value).trim();
const normalizeSequence = (value) => String(value ?? '').replace(/\s+/g, '').toUpperCase();

const attr = (element, name) => normalizeText(element?.getAttribute?.(name));
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
const integerAttr = (element, name) => integerAttrStatus(element, name).value;

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
    isInternalProteinDisplayId(text) ||
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

const ORTHOGROUP_ID_KEYS = Object.freeze(['id', 'orthogroupId', 'orthogroup_id']);
const orthogroupIdStatus = (source) => textAliasStatus(source, ORTHOGROUP_ID_KEYS);

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

const groupScopeValue = (group) => {
  const status = groupScopeStatus(group);
  return status.valid ? status.value : '';
};

const getOrthogroupById = (orthogroups, orthogroupId, groupScope = '') => {
  const id = normalizeText(orthogroupId);
  if (!id) return null;
  const requestedScope = normalizeText(groupScope)
    ? normalizeGroupMetadataScope(groupScope)
    : '';
  const matches = [];
  for (const entry of Array.isArray(orthogroups) ? orthogroups : []) {
    const status = orthogroupIdStatus(entry);
    const scopeStatus = groupScopeStatus(entry);
    if (
      status.valid
      && scopeStatus.valid
      && status.value === id
      && (!requestedScope || scopeStatus.value === requestedScope)
    ) matches.push(entry);
  }
  return matches.length === 1 ? matches[0] : null;
};

const memberFeatureSvgId = (member) => {
  const status = textAliasStatus(member, STABLE_FEATURE_ID_KEYS);
  return status.valid ? status.value : '';
};

const featureRenderedSvgId = (feature) => {
  const status = textAliasStatus(feature, [
    ...RENDERED_FEATURE_ID_KEYS,
    'svg_id',
    'svgId'
  ]);
  return status.valid ? status.value : '';
};

const getFeatureLookupValues = (featureLookup) => {
  if (!featureLookup) return [];
  if (typeof featureLookup.values === 'function') return Array.from(featureLookup.values());
  if (Array.isArray(featureLookup)) return featureLookup;
  if (typeof featureLookup === 'object') return Object.values(featureLookup);
  return [];
};

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
    queryInterval: intervalText(queryStart, queryEnd),
    subjectInterval: intervalText(subjectStart, subjectEnd),
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

const featureOrthogroupId = (feature) => {
  const status = featureOrthogroupIdStatus(feature);
  return status.valid ? status.value : '';
};

const featureProduct = (feature) => firstText(
  feature?.product,
  qualifierFirstValue(feature, 'product'),
  feature?.orthogroupMember?.product,
  feature?.orthogroup_member?.product,
  feature?.note,
  qualifierFirstValue(feature, 'note')
);

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

const addBucketValue = (buckets, key, value) => {
  if (!key) return;
  if (!buckets.has(key)) buckets.set(key, []);
  buckets.get(key).push(value);
};

const identityIndexKeys = (identity) => {
  if (!identity?.usable) return [];
  const recordIndex = identity.recordIndex.supplied ? identity.recordIndex.value : null;
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

const createPairwisePayloadContext = ({
  featureLookup = new Map(),
  sourceFeatures = [],
  orthogroups = [],
  descriptor = null
} = {}) => {
  const identityCache = new WeakMap();
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
    return modes.get(mode);
  };

  const renderedRecords = [];
  const renderedByObject = new WeakMap();
  const renderedById = new Map();
  const featuresByOrthogroupId = new Map();
  const featureOrthogroupStatuses = new WeakMap();
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
  const sourceByObject = new WeakMap();
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

  const rejectedRenderedIds = new Set();
  const rawRenderedRecordForId = (renderedSvgId) => {
    const candidates = renderedById.get(normalizeText(renderedSvgId)) || [];
    return candidates.length === 1 ? candidates[0] : null;
  };
  for (const role of ['query', 'subject']) {
    const endpoint = descriptor?.endpoints?.[role];
    if (!endpoint?.constrained) continue;
    const renderedIds = splitMetadataValues(descriptor[`${role}FeatureSvgId`]);
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

  const renderedRecordForId = (renderedSvgId) => {
    const id = normalizeText(renderedSvgId);
    if (!id || rejectedRenderedIds.has(id)) return null;
    return rawRenderedRecordForId(id);
  };
  const featureForLookupId = (renderedSvgId) => {
    const id = normalizeText(renderedSvgId);
    if (!id || rejectedRenderedIds.has(id)) return null;
    const feature = featureLookup?.get?.(id) || null;
    const record = feature ? renderedByObject.get(feature) : null;
    return record && rejectedRenderedIds.has(record.identity.renderedId.value)
      ? null
      : feature;
  };

  const groupMemberIndexes = new WeakMap();
  const memberIndexForGroup = (group) => {
    if (!group || typeof group !== 'object') return { records: [], byIdentity: new Map() };
    const cached = groupMemberIndexes.get(group);
    if (cached) return cached;
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
  const groupMemberForFeatureId = (group, featureSvgId) => {
    const endpoint = renderedRecordForId(featureSvgId);
    if (!endpoint) return null;
    const memberIndex = memberIndexForGroup(group);
    const seen = new Set();
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
  const groupHasFeatureIds = (group, featureSvgIds) => (
    splitMetadataValues(featureSvgIds)
      .some((featureSvgId) => Boolean(groupMemberForFeatureId(group, featureSvgId)))
  );

  const groupsById = new Map();
  const groupsByScopedId = new Map();
  const implicitGroupMatches = [];
  for (const group of Array.isArray(orthogroups) ? orthogroups : []) {
    const idStatus = orthogroupIdStatus(group);
    const scopeStatus = groupScopeStatus(group);
    if (!idStatus.valid || !idStatus.value) continue;
    if (scopeStatus.valid) {
      addBucketValue(groupsById, idStatus.value, group);
      addBucketValue(
        groupsByScopedId,
        `${idStatus.value}\u001f${scopeStatus.value}`,
        group
      );
    }
    if (
      descriptor
      && !normalizeText(descriptor.orthogroupId)
      && descriptor.matchKind !== 'collinear'
      && (
        groupHasFeatureIds(group, descriptor.queryFeatureSvgId)
        || groupHasFeatureIds(group, descriptor.subjectFeatureSvgId)
      )
    ) implicitGroupMatches.push(group);
  }

  const orthogroupById = (orthogroupId, groupScope = '') => {
    const id = normalizeText(orthogroupId);
    if (!id) return null;
    const requestedScope = normalizeText(groupScope)
      ? normalizeGroupMetadataScope(groupScope)
      : '';
    const candidates = requestedScope
      ? groupsByScopedId.get(`${id}\u001f${requestedScope}`) || []
      : groupsById.get(id) || [];
    return candidates.length === 1 ? candidates[0] : null;
  };
  const orthogroupForMatch = () => {
    const direct = orthogroupById(descriptor?.orthogroupId, descriptor?.groupScope);
    if (direct) return direct;
    if (normalizeText(descriptor?.orthogroupId)) return null;
    return implicitGroupMatches.length === 1 ? implicitGroupMatches[0] : null;
  };

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
    const memberRendered = memberIdentity.renderedId;
    const sourceRendered = sourceRecord.sourceRendered;
    if (!memberRendered.valid || !sourceRendered.valid) return null;
    if (
      memberRendered.supplied
      && sourceRendered.supplied
      && memberRendered.value !== sourceRendered.value
    ) return null;
    let renderedSvgId = memberRendered.value || sourceRendered.value;
    if (!renderedSvgId) {
      if (!sourceRecord.legacyRendered.valid) return null;
      renderedSvgId = sourceRecord.legacyRendered.value;
    }
    const endpoint = renderedRecordForId(renderedSvgId);
    if (
      !endpoint
      || !identityMatches(memberIdentity, endpoint.identity, { includeRendered: true })
      || !identityMatches(sourceIdentity, endpoint.identity)
    ) return null;
    return endpoint.feature;
  };
  const featureForMember = (member) => {
    const memberIdentity = cachedIdentity(member, 'standard');
    const candidateKey = identityIndexKeys(memberIdentity)[0];
    if (!candidateKey) return null;
    const matches = (sourceByIdentity.get(candidateKey) || []).filter((record) => {
      if (!identityMatches(memberIdentity, record.candidateIdentity)) return false;
      if (!memberIdentity.renderedId.supplied && !record.candidateIdentity.renderedId.supplied) {
        return true;
      }
      return Boolean(renderedFeatureForMember(member, record.feature));
    });
    return matches.length === 1 ? matches[0].feature : null;
  };

  const identityForFeature = (feature) => {
    const standard = cachedIdentity(feature, 'standard');
    return standard.usable ? standard : cachedIdentity(feature, 'legacy');
  };
  const orthogroupStatusForFeature = (feature) => (
    featureOrthogroupStatuses.get(feature) || featureOrthogroupIdStatus(feature)
  );
  const fallbackFeaturesForOrthogroup = (orthogroupId, queryFeature, subjectFeature) => {
    const id = normalizeText(orthogroupId);
    const matching = featuresByOrthogroupId.get(id) || [];
    if (matching.length) return matching;
    return [queryFeature, subjectFeature].filter((feature) => {
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
    groupMemberRecords: (group) => memberIndexForGroup(group).records,
    identityForFeature,
    orthogroupById,
    orthogroupForMatch,
    renderedFeatureForMember
  };
};

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

const getRenderedFeatureForMember = (member, feature, featureLookup) => (
  createPairwisePayloadContext({ featureLookup, sourceFeatures: [feature] })
    .renderedFeatureForMember(member, feature)
);

const getFeatureForMember = (member, featureLookup, sourceFeatures = []) => (
  createPairwisePayloadContext({ featureLookup, sourceFeatures }).featureForMember(member)
);

const getGroupMemberForFeatureSvgId = (group, featureSvgId, featureLookup = null) => (
  createPairwisePayloadContext({ featureLookup, orthogroups: [group] })
    .groupMemberForFeatureId(group, featureSvgId)
);

const getOrthogroupForMatch = (orthogroups, options = {}) => (
  createPairwisePayloadContext({
    featureLookup: options.featureLookup,
    orthogroups,
    descriptor: {
      matchKind: 'pairwise',
      orthogroupId: options.orthogroupId,
      groupScope: options.groupScope,
      queryFeatureSvgId: options.queryFeatureSvgId,
      subjectFeatureSvgId: options.subjectFeatureSvgId
    }
  }).orthogroupForMatch()
);

const buildFallbackOrthogroup = ({
  orthogroupId,
  queryFeature,
  subjectFeature,
  featureLookup
}) => buildFallbackOrthogroupWithContext({
  orthogroupId,
  queryFeature,
  subjectFeature,
  context: createPairwisePayloadContext({ featureLookup })
});

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
  const strand = normalizeMemberStrand(member.strand);
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

const orthogroupSequenceFilename = (orthogroupId, displayName, sequenceKind) => {
  const id = normalizeText(orthogroupId) || 'orthogroup';
  const name = makeSafeFilename(displayName || id, id);
  return `${makeSafeFilename(`${id}_${name}_${sequenceKindLabel(sequenceKind)}`)}.${sequenceExtension(sequenceKind)}`;
};

const orthogroupMemberSequenceFilename = (member, orthogroupId, sequenceKind) => {
  const id = normalizeText(orthogroupId) || 'orthogroup';
  const memberId = firstText(
    resolveDisplayProteinId(null, member, ''),
    memberFeatureSvgId(member),
    'member'
  );
  return `${makeSafeFilename(`${id}_${memberId}_${sequenceKindLabel(sequenceKind)}`)}.${sequenceExtension(sequenceKind)}`;
};

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
      return {
        feature: renderedFeature || feature,
        canOpen: Boolean(renderedFeature),
        record: firstText(member?.recordId, member?.record_id),
        coordinates: memberLocationText(member),
        proteinId: resolveDisplayProteinId(null, member),
        role: firstText(member?.role, member?.memberRole, member?.member_role),
        confidence: firstText(member?.confidence, member?.memberConfidence, member?.member_confidence),
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

const resolveFeatureSectionProteinIds = ({
  feature,
  featureSvgIds,
  context,
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
  context
}) => {
  if (!group) return '';
  const values = [];
  uniqueMetadataValues(featureSvgIds).forEach((featureSvgId) => {
    const feature = context.featureForLookupId(featureSvgId);
    const member = group ? context.groupMemberForFeatureId(group, featureSvgId) : null;
    if (!member) return;
    addUniqueDisplayText(values, resolveDisplayProteinId(feature, member, ''));
  });
  return values.join('; ');
};

const buildOrthogroupDetailRows = ({
  orthogroupId,
  idLabel = 'Similarity group ID',
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
  addRow(rows, 'Group paths', orthologPathCount);
  addRow(rows, 'Related edges', relatedEdgeCount);
  return rows;
};

const buildBlockOrthogroups = ({
  orthogroupIds,
  context,
  queryFeatureSvgId,
  subjectFeatureSvgId,
  orthogroupNameOverrides,
  orthogroupDescriptionOverrides,
  groupScope
}) => orthogroupIds.map((orthogroupId) => {
  const group = context.orthogroupById(orthogroupId, groupScope);
  const normalizedGroupScope = normalizeGroupMetadataScope(groupScopeValue(group) || groupScope);
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
  const idLabel = normalizedGroupScope === 'adjacent_local'
    ? 'Collinear group ID'
    : 'Similarity group ID';
  const recordCoverage = firstText(group?.record_coverage_count, group?.recordCoverage);
  const rbhOrthogroups = Array.isArray(group?.rbhOrthogroupIds) ? group.rbhOrthogroupIds : [];
  const orthologPathCount = firstText(
    group?.orthologPathCount,
    Array.isArray(group?.orthologPaths) ? String(group.orthologPaths.length) : ''
  );
  const relatedEdgeCount = firstText(
    group?.relatedEdgeCount,
    Array.isArray(group?.relatedEdges) ? String(group.relatedEdges.length) : ''
  );
  const memberRows = buildOrthogroupMemberRows(group, context, orthogroupId);
  return {
    id: orthogroupId,
    displayName,
    description,
    memberCount,
    recordCoverage,
    queryMember: resolveBlockMemberLabels({ group, featureSvgIds: queryFeatureSvgId, context }),
    subjectMember: resolveBlockMemberLabels({ group, featureSvgIds: subjectFeatureSvgId, context }),
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
  context,
  group
}) => {
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
  addRow(rows, 'Location', firstText(featureLocationText(feature), interval));
  addRow(rows, displayProteinIds.includes(';') ? 'Protein IDs' : 'Protein ID', displayProteinIds);
  addRow(rows, 'Gene', firstText(feature?.gene, qualifierFirstValue(feature, 'gene')));
  addRow(rows, 'Locus tag', firstText(feature?.locus_tag, feature?.locusTag, qualifierFirstValue(feature, 'locus_tag')));
  addRow(rows, 'Product', firstText(feature?.product, qualifierFirstValue(feature, 'product')));
  addRow(rows, 'Locus ID', locusId);
  addRow(rows, 'Display name', displayName);
  return section(title, rows, { featureRows });
};

const buildMatchSpans = (element, matchKind) => {
  if (matchKind === 'orthogroup') return [];
  const queryRecordIndex = integerAttrStatus(element, 'data-query-record-index');
  const subjectRecordIndex = integerAttrStatus(element, 'data-subject-record-index');
  const sourceIndex = integerAttrStatus(element, 'data-source-index');
  const referenceSide = attr(element, 'data-reference-side').toLowerCase();
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
      displayRole: matchKind === 'homology' ? homologyDisplayRole('subject') : 'Subject',
      sourceIndex: sourceIndex.value,
      recordIndex: subjectRecordIndex.value,
      identityValid: subjectRecordIndex.valid && sourceIndex.valid,
      referenceSide
    }
  ];
};

const buildSequenceBundleForMatch = (element, matchKind, matchId, resolveSequenceSource) => {
  const spans = buildMatchSpans(element, matchKind);
  if (!spans.length) return null;
  return buildMatchSequenceBundle(spans, {
    matchId,
    resolveSequenceSource: (sourceKey, recordId, context) => {
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
        ? { origin: 'circular-reference' }
        : { origin: 'homology-comparison', sourceIndex: span.sourceIndex };
    },
    unavailableReasonForSpan: (span) => (
      span.identityValid ? '' : 'The sequence identity metadata for this match is invalid.'
    )
  });
};

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
  const blockOrthogroups = matchKind === 'collinear'
    ? buildBlockOrthogroups({
      orthogroupIds,
      context,
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

  const summaryRows = [];
  addRow(summaryRows, 'Query record', descriptor.queryRecordId);
  addRow(summaryRows, 'Subject record', descriptor.subjectRecordId);
  addRow(summaryRows, 'Query interval', qInterval);
  addRow(summaryRows, 'Subject interval', sInterval);
  addRow(summaryRows, 'Orientation', descriptor.orientation);
  if (matchKind === 'homology') {
    addRow(summaryRows, 'Ring label', attr(element, 'data-track-label'));
    const rawSourceIndex = integerAttr(element, 'data-source-index');
    addRow(summaryRows, 'Source index', rawSourceIndex !== null ? String(rawSourceIndex + 1) : '');
    addRow(summaryRows, 'Reference side', attr(element, 'data-reference-side'));
    addRow(summaryRows, 'Reference record', attr(element, 'data-reference-record-id'));
  }

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
    addRow(blockRows, 'Average identity', descriptor.identity);
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
    addRow(orthogroupRows, 'Similarity group ID', orthogroupId);
    addRow(orthogroupRows, 'Display name', displayName);
    addRow(orthogroupRows, 'Description', description);
    addRow(orthogroupRows, 'Members', firstText(group?.member_count, group?.memberCount));
    addRow(orthogroupRows, 'Record coverage', firstText(group?.record_coverage_count, group?.recordCoverage));
    addRow(orthogroupRows, 'RBH seeds', Array.isArray(group?.rbhOrthogroupIds) ? group.rbhOrthogroupIds.join('; ') : '');
    addRow(orthogroupRows, 'Group paths', firstText(
      group?.orthologPathCount,
      Array.isArray(group?.orthologPaths) ? String(group.orthologPaths.length) : ''
    ));
    addRow(orthogroupRows, 'Related edges', firstText(
      group?.relatedEdgeCount,
      Array.isArray(group?.relatedEdges) ? String(group.relatedEdges.length) : ''
    ));
  }
  const orthogroupMemberRows = buildOrthogroupMemberRows(group, context, orthogroupId);

  if (matchKind === 'orthogroup') {
    const summarySection = section(
      'Summary',
      orthogroupRows,
      orthogroupMemberSectionExtras(orthogroupMemberRows, orthogroupId, displayName)
    );
    return {
      id: firstText(matchId, orthogroupId),
      title: orthogroupTitle(orthogroupId, displayName),
      subtitle: '',
      matchKind,
      orthogroupId,
      collinearityBlockId,
      queryFeatureSvgId,
      subjectFeatureSvgId,
      fill: descriptor.fill,
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
      'Similarity group',
      orthogroupRows,
      orthogroupMemberSectionExtras(orthogroupMemberRows, orthogroupId, displayName)
    ));
  }
  if (matchKind === 'collinear') {
    const blockOrthogroupRows = [];
    const localGroups = normalizedGroupScope === 'adjacent_local';
    addRow(
      blockOrthogroupRows,
      localGroups
        ? 'Number of local collinear groups'
        : 'Number of similarity groups covered',
      String(orthogroupIds.length)
    );
    sections.push(section(
      localGroups ? 'Local collinear groups' : 'Similarity groups covered',
      blockOrthogroupRows,
      { blockOrthogroups }
    ));
  }
  if (matchKind === 'collinear' || blockRows.length > 0) {
    sections.push(section('Collinearity', blockRows));
  }
  if (matchKind !== 'homology') sections.push(buildFeatureRows({
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
  if (matchKind !== 'homology') sections.push(buildFeatureRows({
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

  const sequenceBundle = buildSequenceBundleForMatch(
    element,
    matchKind,
    matchId || 'match',
    resolveSequenceSource
  );

  return {
    id: matchId,
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
    fill: descriptor.fill,
    sequenceTitle: matchKind === 'collinear' ? 'Collinear block spans' : 'Matched sequences',
    sequenceNote: matchKind === 'collinear'
      ? 'Block envelopes may include intergenic sequence and genes that are not anchors.'
      : '',
    sequenceBundle,
    sections: sections.filter((entry) => (
      entry.rows.length > 0 ||
      (Array.isArray(entry.featureRows) && entry.featureRows.length > 0)
    ))
  };
};

// Compatibility export retained for callers and older tests.
export const buildPairwiseMatchPayload = (element, options = {}) =>
  buildMatchPopupPayload(element, options);

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
  const rows = [];
  const addFirst = (label, value) => addRow(rows, label, value);
  if (!matchKind) return rows;
  if (matchKind === 'orthogroup') {
    addFirst('Kind', matchKind);
    addFirst('Similarity group', orthogroupId);
    addFirst('Display name', displayName);
    addFirst('Members', memberCount);
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
    ? getOrthogroupById(
      groupValues,
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

export const buildPairwiseMatchHoverRows = (payload) => {
  if (!payload) return [];
  const summary = payload.sections.find((entry) => entry.title === 'Summary');
  const alignment = payload.sections.find((entry) => entry.title === 'Alignment');
  const block = payload.sections.find((entry) => entry.title === 'Collinearity');
  const orthogroup = payload.sections.find(
    (entry) => entry.title === 'Similarity group'
  );
  const findValue = (sectionEntry, label) => (
    sectionEntry?.rows.find((row) => row.label === label)?.value || ''
  );
  return formatPairwiseMatchHoverRows({
    matchKind: payload.matchKind,
    identity: findValue(alignment, 'Identity') || findValue(block, 'Average identity'),
    queryInterval: findValue(summary, 'Query interval'),
    subjectInterval: findValue(summary, 'Subject interval'),
    orthogroupId: payload.matchKind === 'orthogroup'
      ? findValue(summary, 'Similarity group ID') || payload.orthogroupId
      : findValue(orthogroup, 'Similarity group ID'),
    displayName: findValue(summary, 'Display name'),
    memberCount: findValue(summary, 'Members'),
    groupScope: payload.groupScope,
    blockOrthogroupCount: payload.blockOrthogroupCount,
    collinearityBlockId: findValue(block, 'Block ID')
  });
};
