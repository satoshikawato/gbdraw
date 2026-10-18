// @ts-check
import { runFeatureExtraction } from './diagram-generation.js';
import { linearCardRegion } from './record-options.js';
import { cloneJsonValue } from './json-clone.js';
import { cloneFileBytesForTransfer } from './file-content-cache.js';

/** @import { FileLike } from './session-resource-backing.js' */

/**
 * A feature as the Python extraction returns it: JSON keyed by field name.
 * @typedef {Record<string, unknown>} ExtractedFeature
 */
/** @typedef {ExtractedFeature & { svg_id: string, stable_svg_id: string, stable_feature_id: string }} NormalizedFeature */

/**
 * What the Worker returns for one source, before the module checks it.
 * @typedef {{
 *   features?: ExtractedFeature[], biological_features?: ExtractedFeature[],
 *   record_ids?: string[], error?: unknown
 * }} FeatureExtractionPayload
 */

/**
 * The fields of a Linear record slot (`LinearSeq` in `state.js`) that the
 * region context reads, with the types `state.js` declares.
 * @typedef {{
 *   region_record_id?: string, region_start?: unknown, region_end?: unknown, region_reverse?: boolean
 * }} LinearRegionSeq
 */

/** @typedef {{ cli: string, file: string, displayFile: string }} LinearRegionSpec */

/** A drawn Linear record: its input and its record within that input. @typedef {{ fileIdx: number, recordIdx: number }} FeatureRecordInput */

const FEATURE_EXTRACTION_CACHE_LIMIT = 16;
/** @type {WeakMap<FileLike, Map<string, FeatureExtractionPayload>>} */
const featureExtractionCache = new WeakMap();
/** @type {WeakMap<FileLike, number>} */
const featureExtractionFileIds = new WeakMap();
let nextFeatureExtractionFileId = 1;

// Some hosts have no `performance.now`.
const getNow = () => (/** @type {{ now?: () => number } | undefined} */ (globalThis.performance)?.now
  ? performance.now() : Date.now());

/** @param {unknown} value */
const normalizeRecordSelectorText = (value) => {
  const normalized = String(value ?? '').trim();
  if (!normalized || ['none', 'null', 'jsnull', 'undefined', 'jsundefined', '-'].includes(normalized.toLowerCase())) {
    return '';
  }
  return normalized;
};

/** @param {FileLike | null | undefined} file */
const getFeatureExtractionFileId = (file) => {
  if (!file || (typeof file !== 'object' && typeof file !== 'function')) return '';
  if (!featureExtractionFileIds.has(file)) {
    featureExtractionFileIds.set(file, nextFeatureExtractionFileId);
    nextFeatureExtractionFileId += 1;
  }
  return featureExtractionFileIds.get(file);
};

/** @param {FeatureExtractionPayload | null | undefined} data */
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
  error: data?.error
});

/**
 * @param {{
 *   regionSpec: unknown, recordSelector: unknown, reverseFlag: unknown, selectedFeatures: unknown,
 *   featureVisibility: unknown, format?: string, mode?: string, fastaFile?: File | null,
 *   includeBiologicalFeatures?: boolean
 * }} options
 */
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

/** @param {FileLike | null | undefined} file @param {string} key */
export const getCachedFeatureExtraction = (file, key) => {
  if (!file) return null;
  const byKey = featureExtractionCache.get(file);
  const cached = byKey?.get(key) || null;
  return cached ? cloneFeatureExtractionData(cached) : null;
};

/** @param {FileLike | null | undefined} file @param {string} key @param {FeatureExtractionPayload} value */
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

/**
 * @typedef {object} FeatureExtractionReadOptions
 * @property {string} path
 * @property {File | null | undefined} file
 * @property {string} [format]
 * @property {string | null} [fastaPath]
 * @property {File | null} [fastaFile]
 * @property {string} [mode]
 * @property {string | null | undefined} [regionSpec]
 * @property {string | null | undefined} [recordSelector]
 * @property {boolean} [reverseFlag]
 * @property {string[] | null} [selectedFeatures]
 * @property {boolean} [includeBiologicalFeatures]
 * @property {string | null} [featureVisibilityTablePath]
 * @property {string} [featureVisibilityTsv]
 * @property {{ label: string, ms: number, details: string }[]} [timingEntries]
 * @property {string} [timingLabel]
 * @property {(request: Record<string, any>) => Promise<any>} [runFeatureExtractionImpl]
 */

/**
 * @param {FeatureExtractionReadOptions} options
 * @returns {Promise<FeatureExtractionPayload>}
 */
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
  const buffer = await cloneFileBytesForTransfer(file);
  const requestFiles = [{
    path,
    name: file.name || path.split('/').pop() || 'input.gb',
    bytes: buffer
  }];
  if (fastaPath && fastaFile) {
    requestFiles.push({
      path: fastaPath,
      name: fastaFile.name || fastaPath.split('/').pop() || 'input.fasta',
      bytes: await cloneFileBytesForTransfer(fastaFile)
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
/** @param {unknown} stableSvgId @param {unknown} recordIndex @param {unknown} recordCount */
export const makeLinearRenderedFeatureId = (stableSvgId, recordIndex, recordCount) => {
  const svgId = String(stableSvgId || '').trim();
  if (!svgId || Number(recordCount) <= 1) return svgId;
  return `${svgId}_record_${Number(recordIndex) + 1}`;
};

/** @param {ExtractedFeature | null | undefined} feature */
const stableFeatureId = (feature) => String(
  feature?.stable_feature_id ||
  feature?.stableFeatureId ||
  feature?.stable_svg_id ||
  feature?.stableFeatureSvgId ||
  feature?.svg_id ||
  feature?.svgId ||
  ''
).trim();

/**
 * @param {ExtractedFeature} feature
 * @param {Record<string, unknown>} [extra]
 * @returns {NormalizedFeature}
 */
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

/** @param {FeatureExtractionPayload | null | undefined} featData @returns {ExtractedFeature[]} */
const biologicalFeaturesFromExtraction = (featData) => (
  Array.isArray(featData?.biological_features) && (
    featData.biological_features.length > 0 ||
    !Array.isArray(featData?.features) ||
    featData.features.length === 0
  )
    ? featData.biological_features
    : (Array.isArray(featData?.features) ? featData.features : [])
);

/**
 * Linear card region endpoints: null when unset, else validated numbers
 * (`linearCardRegion` reads the fields); a region Generate cannot draw throws.
 * @param {Pick<LinearRegionSeq, 'region_start' | 'region_end'>} seq
 * @param {number} idx
 */
export const resolveLinearRegionBounds = (seq, idx) => {
  const region = linearCardRegion(seq, idx + 1);
  if (region && 'error' in region) throw region.error;
  return region;
};

/**
 * @param {LinearRegionSeq[] | null} [linearSeqs]
 * @param {string} [lInputType]
 */
export const buildLinearRegionExtractionContext = (linearSeqs = [], lInputType = 'gb') => {
  /** @type {(LinearRegionSpec | null)[]} */
  const regionSpecs = [];
  /** @type {string[]} */
  const recordSelectors = [];
  /** @type {boolean[]} */
  const reverseFlags = [];
  const seqs = Array.isArray(linearSeqs) ? linearSeqs : [];

  seqs.forEach((seq = {}, idx) => {
    const recordIdRaw = seq.region_record_id ? String(seq.region_record_id).trim() : '';
    const wantsReverse = Boolean(seq.region_reverse);
    const bounds = resolveLinearRegionBounds(seq, idx);

    recordSelectors.push(recordIdRaw || '');

    if (bounds) {
      const { start, end } = bounds;
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

/** @param {RegionSpecInput | null | undefined} regionSpec */
const extractionRegionSpec = (regionSpec) => {
  if (!regionSpec) return null;
  if (typeof regionSpec === 'string') return regionSpec || null;
  return regionSpec.displayFile || regionSpec.file || null;
};

/** @typedef {string | { cli?: string, file?: string, displayFile?: string }} RegionSpecInput */

/**
 * @typedef {object} ExtractFeatureMetadataOptions
 * @property {string} mode
 * @property {string} cInputType
 * @property {string} lInputType
 * @property {File | null} [circularFile]
 * @property {File | null} [circularFastaFile]
 * @property {Record<string, any>[]} [linearSeqs]
 * @property {Array<RegionSpecInput | null> | null} [regionSpecs]
 * @property {string[] | null} [recordSelectors]
 * @property {boolean[] | null} [reverseFlags]
 * @property {string | null} [featureVisibilityTablePath]
 * @property {string} [featureVisibilityTsv]
 * @property {string[] | null} [selectedFeatures]
 * @property {(feature: ExtractedFeature, inputIndex?: number) => ExtractedFeature | null | undefined} [enrichFeature]
 * @property {typeof readFeatureExtractionData} [readFeatureExtractionDataImpl]
 * @property {{ label: string, ms: number, details: string }[]} [timingEntries]
 */

/** @param {ExtractFeatureMetadataOptions} options */
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
  enrichFeature = /** @type {(feature: ExtractedFeature, inputIndex?: number) => ExtractedFeature | null | undefined} */ ((feature) => feature),
  readFeatureExtractionDataImpl = readFeatureExtractionData,
  timingEntries = []
}) => {
  /** @type {{ inputIndex: number, message: string }[]} */
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
        featureRecordIds: [],
        selectedFeatureRecordIdx: 0,
        errors
      };
    }
    return {
      extractedFeatures: featData.features.map((feature) => normalizeBiologicalFeature(feature)),
      biologicalFeatures: biologicalFeaturesFromExtraction(featData)
        .map((feature) => normalizeBiologicalFeature(feature)),
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
    /** @type {NormalizedFeature[]} */
    let allFeatures = [];
    /** @type {NormalizedFeature[]} */
    let allBiologicalFeatures = [];
    /** @type {{ label: string, fileIdx: number, recordIdx: number }[]} */
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

      /** @param {ExtractedFeature} feature */
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

      (featData.record_ids || []).forEach((rid, ridx) => {
        allRecordLabels.push({ label: `File ${i + 1}: ${rid}`, fileIdx: i, recordIdx: ridx });
      });
    }

    return {
      extractedFeatures: allFeatures,
      biologicalFeatures: allBiologicalFeatures,
      featureRecordIds: allRecordLabels.map((record) => record.label),
      /** @type {FeatureRecordInput[]} The drawn records in order. */
      featureRecordInputs: allRecordLabels.map(({ fileIdx, recordIdx }) => ({ fileIdx, recordIdx })),
      selectedFeatureRecordIdx: 0,
      errors
    };
  }

  return {
    extractedFeatures: [],
    biologicalFeatures: [],
    featureRecordIds: [],
    selectedFeatureRecordIdx: 0,
    errors
  };
};
