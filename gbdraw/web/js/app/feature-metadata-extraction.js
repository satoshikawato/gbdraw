import { runFeatureExtraction } from '../services/diagram-generation.js';
import { cloneJsonValue } from '../services/json-clone.js';

const FEATURE_EXTRACTION_CACHE_LIMIT = 16;
const featureExtractionCache = new WeakMap();
const featureExtractionFileIds = new WeakMap();
let nextFeatureExtractionFileId = 1;

const getNow = () => (globalThis.performance?.now ? performance.now() : Date.now());

const normalizeRecordSelectorText = (value) => {
  const normalized = String(value ?? '').trim();
  if (!normalized || ['none', 'null', 'jsnull', 'undefined', 'jsundefined', '-'].includes(normalized.toLowerCase())) {
    return '';
  }
  return normalized;
};

const getFeatureExtractionFileId = (file) => {
  if (!file || (typeof file !== 'object' && typeof file !== 'function')) return '';
  if (!featureExtractionFileIds.has(file)) {
    featureExtractionFileIds.set(file, nextFeatureExtractionFileId);
    nextFeatureExtractionFileId += 1;
  }
  return featureExtractionFileIds.get(file);
};

export const cloneFeatureExtractionData = (data) => ({
  features: Array.isArray(data?.features)
    ? data.features.map((feature) => cloneJsonValue(
        feature,
        feature && typeof feature === 'object' ? { ...feature } : feature
      ))
    : [],
  biological_features: Array.isArray(data?.biological_features)
    ? data.biological_features.map((feature) => cloneJsonValue(
        feature,
        feature && typeof feature === 'object' ? { ...feature } : feature
      ))
    : [],
  record_ids: Array.isArray(data?.record_ids) ? [...data.record_ids] : [],
  selector_safety_scope: Array.isArray(data?.selector_safety_scope)
    ? cloneJsonValue(data.selector_safety_scope, [])
    : [],
  error: data?.error
});

export const buildFeatureExtractionCacheKey = ({
  regionSpec,
  recordSelector,
  reverseFlag,
  selectedFeatures,
  featureVisibility,
  format = 'genbank',
  mode = 'linear',
  fastaFile = null,
  includeBiologicalFeatures = false
}) =>
  JSON.stringify({
    regionSpec: String(regionSpec || ''),
    recordSelector: normalizeRecordSelectorText(recordSelector),
    reverseFlag: reverseFlag ? '1' : '0',
    selectedFeatures: Array.isArray(selectedFeatures) && selectedFeatures.length ? selectedFeatures : 'all',
    featureVisibility: String(featureVisibility || ''),
    format: String(format || 'genbank'),
    mode: String(mode || 'linear'),
    fastaFileId: getFeatureExtractionFileId(fastaFile),
    includeBiologicalFeatures: Boolean(includeBiologicalFeatures)
  });

export const getCachedFeatureExtraction = (file, key) => {
  if (!file) return null;
  const byKey = featureExtractionCache.get(file);
  const cached = byKey?.get(key) || null;
  return cached ? cloneFeatureExtractionData(cached) : null;
};

export const setCachedFeatureExtraction = (file, key, value) => {
  if (!file || value?.error) return;
  let byKey = featureExtractionCache.get(file);
  if (!byKey) {
    byKey = new Map();
    featureExtractionCache.set(file, byKey);
  }
  if (byKey.size >= FEATURE_EXTRACTION_CACHE_LIMIT) byKey.delete(byKey.keys().next().value);
  byKey.set(key, cloneFeatureExtractionData(value));
};

export const readFeatureExtractionData = async ({
  path,
  file,
  format = 'genbank',
  fastaPath = null,
  fastaFile = null,
  mode = 'linear',
  regionSpec,
  recordSelector,
  reverseFlag,
  selectedFeatures = null,
  includeBiologicalFeatures = false,
  featureVisibilityTablePath = null,
  featureVisibilityTsv = '',
  timingEntries = [],
  timingLabel = 'feature extraction',
  runFeatureExtractionImpl = runFeatureExtraction
}) => {
  const cacheKey = buildFeatureExtractionCacheKey({
    regionSpec,
    recordSelector,
    reverseFlag,
    selectedFeatures,
    featureVisibility: featureVisibilityTsv,
    format,
    mode,
    fastaFile,
    includeBiologicalFeatures
  });
  const cached = getCachedFeatureExtraction(file, cacheKey);
  if (cached) {
    if (Array.isArray(timingEntries)) timingEntries.push({ label: timingLabel, ms: 0, details: 'cache hit' });
    return cached;
  }

  if (!file) {
    throw new Error('Feature extraction input file is not available.');
  }
  if (String(format || '').toLowerCase() === 'gff' && (!fastaPath || !fastaFile)) {
    throw new Error('GFF3 feature extraction requires a paired FASTA file.');
  }

  const startedAt = getNow();
  const buffer = await file.arrayBuffer();
  const requestFiles = [{
    path,
    name: file.name || path.split('/').pop() || 'input.gb',
    bytes: buffer
  }];
  if (fastaPath && fastaFile) {
    requestFiles.push({
      path: fastaPath,
      name: fastaFile.name || fastaPath.split('/').pop() || 'input.fasta',
      bytes: await fastaFile.arrayBuffer()
    });
  }
  if (featureVisibilityTablePath && featureVisibilityTsv) {
    requestFiles.push({
      path: featureVisibilityTablePath,
      name: featureVisibilityTablePath.split('/').pop() || 'feature_visibility.tsv',
      bytes: new TextEncoder().encode(featureVisibilityTsv).buffer
    });
  }

  const response = await runFeatureExtractionImpl({
    path,
    format,
    fastaPath,
    mode,
    files: requestFiles,
    regionSpec: regionSpec || null,
    recordSelector: recordSelector || null,
    reverseFlag: Boolean(reverseFlag),
    selectedFeatures: Array.isArray(selectedFeatures) && selectedFeatures.length ? selectedFeatures : null,
    featureVisibilityTablePath: featureVisibilityTablePath || null,
    includeBiologicalFeatures: Boolean(includeBiologicalFeatures)
  });
  const featData = response?.result ?? response;
  if (Array.isArray(timingEntries)) {
    timingEntries.push({ label: timingLabel, ms: getNow() - startedAt, details: 'worker' });
  }
  setCachedFeatureExtraction(file, cacheKey, featData);
  return featData;
};

// Retained for compatibility with older callers. Source metadata must not use this
// prediction; actual rendered IDs are read back from the generated SVG.
export const makeLinearRenderedFeatureId = (stableSvgId, recordIndex, recordCount) => {
  const svgId = String(stableSvgId || '').trim();
  if (!svgId || Number(recordCount) <= 1) return svgId;
  return `${svgId}_record_${Number(recordIndex) + 1}`;
};

const stableFeatureId = (feature) => String(
  feature?.stable_feature_id ||
  feature?.stableFeatureId ||
  feature?.stable_svg_id ||
  feature?.stableFeatureSvgId ||
  feature?.svg_id ||
  feature?.svgId ||
  ''
).trim();

const normalizeBiologicalFeature = (feature, extra = {}) => {
  const stableId = stableFeatureId(feature);
  return {
    ...feature,
    ...extra,
    svg_id: stableId,
    stable_svg_id: stableId,
    stable_feature_id: stableId
  };
};

const biologicalFeaturesFromExtraction = (featData) => (
  Array.isArray(featData?.biological_features) && (
    featData.biological_features.length > 0 ||
    !Array.isArray(featData?.features) ||
    featData.features.length === 0
  )
    ? featData.biological_features
    : (Array.isArray(featData?.features) ? featData.features : [])
);

const hasRegionValue = (value) => value !== null && value !== undefined && value !== '';

export const buildLinearRegionExtractionContext = (linearSeqs = [], lInputType = 'gb') => {
  const regionSpecs = [];
  const recordSelectors = [];
  const reverseFlags = [];
  const seqs = Array.isArray(linearSeqs) ? linearSeqs : [];

  seqs.forEach((seq = {}, idx) => {
    const hasStart = hasRegionValue(seq.region_start);
    const hasEnd = hasRegionValue(seq.region_end);
    const recordIdRaw = seq.region_record_id ? String(seq.region_record_id).trim() : '';
    const wantsReverse = Boolean(seq.region_reverse);

    if (hasStart !== hasEnd) {
      throw new Error(`Sequence #${idx + 1}: Provide both Region start and end, or leave both empty.`);
    }

    recordSelectors.push(recordIdRaw || '');

    if (hasStart && hasEnd) {
      const start = Number(seq.region_start);
      const end = Number(seq.region_end);
      if (!Number.isFinite(start) || !Number.isFinite(end)) {
        throw new Error(`Sequence #${idx + 1}: Region start/end must be numbers.`);
      }
      if (!Number.isInteger(start) || !Number.isInteger(end)) {
        throw new Error(`Sequence #${idx + 1}: Region start/end must be integers.`);
      }
      if (start < 1 || end < 1) {
        throw new Error(`Sequence #${idx + 1}: Region start/end must be >= 1.`);
      }
      const canonicalStart = Math.min(start, end);
      const canonicalEnd = Math.max(start, end);
      const specBody = `${start}-${end}${wantsReverse ? ':rc' : ''}`;
      const canonicalSpecBody = `${canonicalStart}-${canonicalEnd}`;
      const fileLabel = String(lInputType || 'gb') === 'gb' ? `seq_${idx}.gb` : `seq_${idx}.gff`;
      const cliSpec = recordIdRaw ? `${fileLabel}:${recordIdRaw}:${specBody}` : `#${idx + 1}:${specBody}`;
      regionSpecs.push({ cli: cliSpec, file: canonicalSpecBody, displayFile: specBody });
      reverseFlags.push(false);
      return;
    }

    regionSpecs.push(null);
    reverseFlags.push(wantsReverse);
  });

  return { regionSpecs, recordSelectors, reverseFlags };
};

const extractionRegionSpec = (regionSpec) => {
  if (!regionSpec) return null;
  if (typeof regionSpec === 'string') return regionSpec || null;
  return regionSpec.displayFile || regionSpec.file || null;
};

export const extractFeatureMetadataForPreview = async ({
  mode,
  cInputType,
  lInputType,
  circularFile,
  circularFastaFile = null,
  linearSeqs = [],
  regionSpecs = null,
  recordSelectors = null,
  reverseFlags = null,
  featureVisibilityTablePath = null,
  featureVisibilityTsv = '',
  selectedFeatures = null,
  enrichFeature = (feature) => feature,
  readFeatureExtractionDataImpl = readFeatureExtractionData,
  timingEntries = []
}) => {
  const errors = [];

  if (mode === 'circular' && ['gb', 'gff'].includes(cInputType)) {
    const isGff = cInputType === 'gff';
    const featData = await readFeatureExtractionDataImpl({
      path: isGff ? '/input.gff' : '/input.gb',
      file: circularFile || null,
      format: isGff ? 'gff' : 'genbank',
      fastaPath: isGff ? '/input.fasta' : null,
      fastaFile: isGff ? circularFastaFile : null,
      mode: 'circular',
      regionSpec: null,
      recordSelector: null,
      reverseFlag: false,
      selectedFeatures,
      includeBiologicalFeatures: true,
      featureVisibilityTablePath,
      featureVisibilityTsv,
      timingEntries,
      timingLabel: 'feature extraction circular input'
    });
    if (featData?.error) {
      errors.push({ inputIndex: 0, message: String(featData.error) });
      return {
        extractedFeatures: [],
        biologicalFeatures: [],
        featureSelectorSafetyScope: [],
        featureRecordIds: [],
        selectedFeatureRecordIdx: 0,
        errors
      };
    }
    if (!Array.isArray(featData?.features)) {
      errors.push({ inputIndex: 0, message: 'Feature extraction returned an unexpected payload.' });
      return {
        extractedFeatures: [],
        biologicalFeatures: [],
        featureSelectorSafetyScope: [],
        featureRecordIds: [],
        selectedFeatureRecordIdx: 0,
        errors
      };
    }
    return {
      extractedFeatures: featData.features.map((feature) => normalizeBiologicalFeature(feature)),
      biologicalFeatures: biologicalFeaturesFromExtraction(featData)
        .map((feature) => normalizeBiologicalFeature(feature)),
      featureSelectorSafetyScope: Array.isArray(featData?.selector_safety_scope)
        ? featData.selector_safety_scope
        : [],
      featureRecordIds: featData?.record_ids || [],
      selectedFeatureRecordIdx: 0,
      errors
    };
  }

  if (mode === 'linear' && ['gb', 'gff'].includes(lInputType) && Array.isArray(linearSeqs) && linearSeqs.length > 0) {
    const isGff = lInputType === 'gff';
    const regionContext = regionSpecs && recordSelectors && reverseFlags
      ? { regionSpecs, recordSelectors, reverseFlags }
      : buildLinearRegionExtractionContext(linearSeqs, lInputType);
    let allFeatures = [];
    let allBiologicalFeatures = [];
    let allSelectorSafetyScope = [];
    const allRecordLabels = [];
    for (let i = 0; i < linearSeqs.length; i += 1) {
      const seq = linearSeqs[i] || {};
      const featData = await readFeatureExtractionDataImpl({
        path: isGff ? `/seq_${i}.gff` : `/seq_${i}.gb`,
        file: isGff ? seq.gff || null : seq.gb || null,
        format: isGff ? 'gff' : 'genbank',
        fastaPath: isGff ? `/seq_${i}.fasta` : null,
        fastaFile: isGff ? seq.fasta || null : null,
        mode: 'linear',
        regionSpec: extractionRegionSpec(regionContext.regionSpecs[i]),
        recordSelector: regionContext.recordSelectors[i] ?? '',
        reverseFlag: Boolean(regionContext.reverseFlags[i]),
        selectedFeatures,
        includeBiologicalFeatures: true,
        featureVisibilityTablePath,
        featureVisibilityTsv,
        timingEntries,
        timingLabel: `feature extraction linear input #${i + 1}`
      });

      if (featData?.error) {
        errors.push({ inputIndex: i, message: String(featData.error) });
        continue;
      }
      if (!Array.isArray(featData?.features)) {
        errors.push({ inputIndex: i, message: 'Feature extraction returned an unexpected payload.' });
        continue;
      }

      const normalizeLinearFeature = (feature) => {
        const enriched = enrichFeature(feature, i) || feature;
        return normalizeBiologicalFeature(enriched, {
          fileIdx: i,
          displayRecordId: `File ${i + 1}: ${feature.record_id}`,
          id: `file${i}_${feature.id}`
        });
      };
      const features = featData.features.map(normalizeLinearFeature);
      const biologicalFeatures = biologicalFeaturesFromExtraction(featData)
        .map(normalizeLinearFeature);
      allFeatures = allFeatures.concat(features);
      allBiologicalFeatures = allBiologicalFeatures.concat(biologicalFeatures);

      if (Array.isArray(featData.selector_safety_scope)) {
        allSelectorSafetyScope = allSelectorSafetyScope.concat(
          featData.selector_safety_scope.map((entry) => ({ ...entry, fileIdx: i }))
        );
      }
      (featData.record_ids || []).forEach((rid, ridx) => {
        allRecordLabels.push({ label: `File ${i + 1}: ${rid}`, fileIdx: i, recordIdx: ridx });
      });
    }

    return {
      extractedFeatures: allFeatures,
      biologicalFeatures: allBiologicalFeatures,
      featureSelectorSafetyScope: allSelectorSafetyScope,
      featureRecordIds: allRecordLabels.map((record) => record.label),
      selectedFeatureRecordIdx: 0,
      errors
    };
  }

  return {
    extractedFeatures: [],
    biologicalFeatures: [],
    featureSelectorSafetyScope: [],
    featureRecordIds: [],
    selectedFeatureRecordIdx: 0,
    errors
  };
};
