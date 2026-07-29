import { cloneJsonData } from './json-clone.js';

export const FEATURE_CATALOG_SCHEMA = 1;
export const FEATURE_CATALOG_ENCODING = 'biological-authority-v1';

const isObject = (value) => (
  value !== null && typeof value === 'object' && !Array.isArray(value)
);

const jsonEqual = (left, right) => {
  if (left === right) return true;
  if (Array.isArray(left) || Array.isArray(right)) {
    if (!Array.isArray(left) || !Array.isArray(right) || left.length !== right.length) {
      return false;
    }
    return left.every((value, index) => jsonEqual(value, right[index]));
  }
  if (!isObject(left) || !isObject(right)) return false;
  const leftKeys = Object.keys(left);
  const rightKeys = Object.keys(right);
  if (leftKeys.length !== rightKeys.length) return false;
  return leftKeys.every(
    (key) => Object.prototype.hasOwnProperty.call(right, key) && jsonEqual(left[key], right[key])
  );
};

const featureCatalogKey = (feature) => {
  if (!isObject(feature) || !Number.isInteger(feature.record_idx)) return null;
  const stableId = String(
    feature.stable_svg_id || feature.stable_feature_id || feature.svg_id || ''
  );
  return stableId ? `${feature.record_idx}\u0000${stableId}` : null;
};

const firstQualifier = (qualifiers, key) => {
  const values = isObject(qualifiers) && Array.isArray(qualifiers[key])
    ? qualifiers[key]
    : [];
  const value = values.find((item) => item !== '');
  return value === undefined ? '' : value;
};

const featureQualifiersAreStrings = (value) => (
  isObject(value)
  && Object.entries(value).every(
    ([key, items]) => (
      typeof key === 'string'
      && Array.isArray(items)
      && items.every((item) => typeof item === 'string')
    )
  )
);

const truncateCodePoints = (value, limit) => Array.from(value).slice(0, limit).join('');

const expandBiologicalFeature = (feature, index, profile) => {
  if (!isObject(feature) || !isObject(feature.qualifiers) || !isObject(feature.selector)) {
    throw new Error('Compact session biological features require qualifier and selector objects.');
  }
  if (
    !featureQualifiersAreStrings(feature.qualifiers)
    || (
      Object.prototype.hasOwnProperty.call(feature.selector, 'qualifiers')
      && !featureQualifiersAreStrings(feature.selector.qualifiers)
    )
  ) {
    throw new Error('Compact session feature qualifiers must contain string arrays.');
  }
  const expanded = cloneJsonData(feature);
  const svgId = typeof expanded.svg_id === 'string' ? expanded.svg_id : '';
  if (!svgId) throw new Error('Compact session biological features require svg_id.');

  if (!Object.prototype.hasOwnProperty.call(expanded, 'id')) expanded.id = `f${index}`;
  if (!Object.prototype.hasOwnProperty.call(expanded, 'stable_svg_id')) {
    expanded.stable_svg_id = svgId;
  }
  if (!Object.prototype.hasOwnProperty.call(expanded, 'stable_feature_id')) {
    expanded.stable_feature_id = svgId;
  }
  [
    'protein_id',
    'locus_tag',
    'gene_id',
    'old_locus_tag',
    'gene',
    'product'
  ].forEach((field) => {
    if (!Object.prototype.hasOwnProperty.call(expanded, field)) {
      expanded[field] = firstQualifier(expanded.qualifiers, field);
    }
  });
  if (!Object.prototype.hasOwnProperty.call(expanded, 'source_protein_id')) {
    expanded.source_protein_id = Object.prototype.hasOwnProperty.call(
      expanded,
      'protein_id'
    )
      ? expanded.protein_id
      : '';
  }
  if (!Object.prototype.hasOwnProperty.call(expanded, 'note')) {
    expanded.note = truncateCodePoints(firstQualifier(expanded.qualifiers, 'note'), 50);
  }
  if (!Object.prototype.hasOwnProperty.call(expanded, 'sequence_warnings')) {
    expanded.sequence_warnings = [];
  }
  if (!Object.prototype.hasOwnProperty.call(expanded.selector, 'qualifiers')) {
    expanded.selector.qualifiers = cloneJsonData(expanded.qualifiers);
  }
  if (!Object.prototype.hasOwnProperty.call(expanded.selector, 'hash')) {
    expanded.selector.hash = svgId;
  }
  if (profile === 'rich-v1') {
    if (!Object.prototype.hasOwnProperty.call(expanded, 'amino_acid_sequence')) {
      expanded.amino_acid_sequence = firstQualifier(expanded.qualifiers, 'translation');
    }
    if (typeof expanded.nucleotide_sequence !== 'string') {
      throw new Error('Compact rich biological features require nucleotide sequences.');
    }
  } else if (profile !== 'sanitized-v1') {
    throw new Error(`Unsupported compact feature catalog profile: ${profile}.`);
  }
  return expanded;
};

const compactBiologicalFeature = (feature, index, profile) => {
  if (!isObject(feature) || !isObject(feature.qualifiers) || !isObject(feature.selector)) {
    return null;
  }
  if (
    !featureQualifiersAreStrings(feature.qualifiers)
    || (
      Object.prototype.hasOwnProperty.call(feature.selector, 'qualifiers')
      && !featureQualifiersAreStrings(feature.selector.qualifiers)
    )
  ) {
    return null;
  }
  const compact = cloneJsonData(feature);
  const svgId = typeof compact.svg_id === 'string' ? compact.svg_id : '';
  if (!svgId) return null;

  if (compact.id === `f${index}`) delete compact.id;
  if (compact.stable_svg_id === svgId) delete compact.stable_svg_id;
  if (compact.stable_feature_id === svgId) delete compact.stable_feature_id;
  if (compact.source_protein_id === compact.protein_id) delete compact.source_protein_id;
  if (Array.isArray(compact.sequence_warnings) && compact.sequence_warnings.length === 0) {
    delete compact.sequence_warnings;
  }
  if (jsonEqual(compact.selector.qualifiers, compact.qualifiers)) {
    delete compact.selector.qualifiers;
  }
  if (compact.selector.hash === svgId) delete compact.selector.hash;

  [
    'protein_id',
    'locus_tag',
    'gene_id',
    'old_locus_tag',
    'gene',
    'product'
  ].forEach((field) => {
    if (compact[field] === firstQualifier(compact.qualifiers, field)) delete compact[field];
  });
  if (compact.note === truncateCodePoints(firstQualifier(compact.qualifiers, 'note'), 50)) {
    delete compact.note;
  }
  if (
    profile === 'rich-v1'
    && compact.amino_acid_sequence === firstQualifier(compact.qualifiers, 'translation')
  ) {
    delete compact.amino_acid_sequence;
  }

  try {
    return jsonEqual(expandBiologicalFeature(compact, index, profile), feature)
      ? compact
      : null;
  } catch (_error) {
    return null;
  }
};

const compactFeatureCatalog = (features) => {
  if (!isObject(features) || Object.prototype.hasOwnProperty.call(features, 'featureCatalog')) {
    return features;
  }
  const extracted = features.extractedFeatures;
  const biological = features.biologicalFeatures;
  if (
    !Array.isArray(extracted)
    || extracted.length === 0
    || !Array.isArray(biological)
    || biological.length === 0
    || !extracted.every(isObject)
    || !biological.every(isObject)
  ) {
    return features;
  }

  const rich = biological.every(
    (feature) => Object.prototype.hasOwnProperty.call(feature, 'nucleotide_sequence')
      && Object.prototype.hasOwnProperty.call(feature, 'amino_acid_sequence')
  );
  const sanitized = biological.every(
    (feature) => !Object.prototype.hasOwnProperty.call(feature, 'nucleotide_sequence')
      && !Object.prototype.hasOwnProperty.call(feature, 'amino_acid_sequence')
  );
  if (!rich && !sanitized) return features;
  const profile = rich ? 'rich-v1' : 'sanitized-v1';

  const biologicalIndexes = new Map();
  for (let index = 0; index < biological.length; index += 1) {
    const key = featureCatalogKey(biological[index]);
    if (key === null || biologicalIndexes.has(key)) return features;
    biologicalIndexes.set(key, index);
  }

  const compactBiological = [];
  for (let index = 0; index < biological.length; index += 1) {
    const compact = compactBiologicalFeature(biological[index], index, profile);
    if (compact === null) return features;
    compactBiological.push(compact);
  }

  const references = [];
  const referencedBiologicalIndexes = new Set();
  for (const feature of extracted) {
    const key = featureCatalogKey(feature);
    const biologicalIndex = key === null ? undefined : biologicalIndexes.get(key);
    if (biologicalIndex === undefined || typeof feature.id !== 'string') return features;
    if (referencedBiologicalIndexes.has(biologicalIndex)) return features;
    referencedBiologicalIndexes.add(biologicalIndex);
    const projected = cloneJsonData(biological[biologicalIndex]);
    delete projected.feature_index;
    projected.id = feature.id;
    const reference = [biologicalIndex, feature.id];
    if (typeof feature.rendered_feature_svg_id === 'string') {
      projected.rendered_feature_svg_id = feature.rendered_feature_svg_id;
      reference.push(feature.rendered_feature_svg_id);
    }
    if (!jsonEqual(projected, feature)) return features;
    references.push(reference);
  }

  const compactFeatures = { ...features };
  delete compactFeatures.extractedFeatures;
  compactFeatures.biologicalFeatures = compactBiological;
  compactFeatures.featureCatalog = {
    schema: FEATURE_CATALOG_SCHEMA,
    encoding: FEATURE_CATALOG_ENCODING,
    profile,
    extracted: references
  };
  return compactFeatures;
};

export const compactSessionFeatureCatalog = (session) => {
  if (
    !isObject(session) ||
    ![37, 38, 39].includes(session.version) ||
    !isObject(session.features)
  ) {
    return session;
  }
  const features = compactFeatureCatalog(session.features);
  return features === session.features ? session : { ...session, features };
};

export const expandSessionFeatureCatalog = (session) => {
  if (!isObject(session) || !isObject(session.features)) return session;
  const features = session.features;
  const catalog = features.featureCatalog;
  if (catalog === undefined) return session;
  if (
    !isObject(catalog)
    || catalog.schema !== FEATURE_CATALOG_SCHEMA
    || catalog.encoding !== FEATURE_CATALOG_ENCODING
    || !['rich-v1', 'sanitized-v1'].includes(catalog.profile)
    || !Array.isArray(catalog.extracted)
    || Object.prototype.hasOwnProperty.call(features, 'extractedFeatures')
    || !Array.isArray(features.biologicalFeatures)
    || !features.biologicalFeatures.every(isObject)
  ) {
    throw new Error('Invalid compact session feature catalog.');
  }

  if (catalog.extracted.length > features.biologicalFeatures.length) {
    throw new Error('Invalid compact extracted-feature reference.');
  }
  const referencedBiologicalIndexes = new Set();
  const validatedReferences = catalog.extracted.map((reference) => {
    if (
      !Array.isArray(reference)
      || ![2, 3].includes(reference.length)
      || !Number.isInteger(reference[0])
      || reference[0] < 0
      || reference[0] >= features.biologicalFeatures.length
      || typeof reference[1] !== 'string'
      || (reference.length === 3 && typeof reference[2] !== 'string')
    ) {
      throw new Error('Invalid compact extracted-feature reference.');
    }
    if (referencedBiologicalIndexes.has(reference[0])) {
      throw new Error('Invalid compact extracted-feature reference.');
    }
    referencedBiologicalIndexes.add(reference[0]);
    return reference;
  });
  const biologicalFeatures = features.biologicalFeatures.map(
    (feature, index) => expandBiologicalFeature(feature, index, catalog.profile)
  );
  const extractedFeatures = validatedReferences.map((reference) => {
    const projected = cloneJsonData(biologicalFeatures[reference[0]]);
    delete projected.feature_index;
    projected.id = reference[1];
    if (reference.length === 3) projected.rendered_feature_svg_id = reference[2];
    return projected;
  });

  const expandedFeatures = { ...features, biologicalFeatures, extractedFeatures };
  delete expandedFeatures.featureCatalog;
  return { ...session, features: expandedFeatures };
};
