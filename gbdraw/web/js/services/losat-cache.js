// @ts-check
import {
  LEGACY_LOSAT_DERIVED_CACHE_SCHEMA,
  LEGACY_PROTEIN_CANDIDATE_SCHEMA,
  LOSAT_DERIVED_CACHE_SCHEMA,
  NUCLEOTIDE_LOSAT_CACHE_SCHEMA,
  PROTEIN_IDENTITY_MANIFEST_SCHEMA,
  PROTEIN_LOSAT_CACHE_SCHEMA
} from '../contracts/schema.generated.js';

// Membership is tested on unvalidated Session values, hence `ReadonlySet<unknown>`.
/** @type {ReadonlySet<unknown>} */
const LEGACY_CANDIDATE_STATES = new Set(['pending', 'promoted', 'rejected']);
const LOSAT_OUTFMT6_COLUMN_COUNT = 12;
const LOSAT_OUTFMT6_INTEGER_COLUMN_INDEXES = new Set([3, 4, 5, 6, 7, 8, 9]);
const STRICT_DECIMAL_INTEGER_RE = /^[+-]?[0-9]+$/;
const STRICT_DECIMAL_NUMBER_RE =
  /^[+-]?(?:[0-9]+(?:\.[0-9]*)?|\.[0-9]+)(?:[eE][+-]?[0-9]+)?$/;
const MANIFEST_PRESENTATION_IDENTITY_FRAGMENTS = [
  'viewfeaturesvgid',
  'viewfeaturehashparts',
  'renderedfeaturesvgid',
  'renderedsvgid'
];

/**
 * The runtime identity stored beside a raw entry (`losat_runtime_record` in Python). The Web's own writer gives
 * `kind: 'losat'`, `source: 'wasm'`, `version: null` and the program; a Session load stores the saved record
 * without checking its fields, so an entry's `runtime` is `unknown` until a reader has tested that it is an
 * object, and then each field of this record is `unknown`.
 * @typedef {{ kind?: unknown, source?: unknown, version?: unknown, path?: unknown, program?: unknown }} LosatCacheRuntime
 */

/**
 * A raw LOSAT cache entry that `classifyRawLosatCacheEntry` accepted as current. Its classification proves
 * `kind` and `text` (and the protein hashes of a protein entry). Generate writes the other fields with their
 * usual types (`program`, `outfmt`, `flow` strings, `args` a list); a Session load stores them unchecked, so
 * they are `unknown` here.
 * @typedef {{
 *   kind: 'raw-losat', text: string, schema?: unknown, identityKind?: unknown, idEncoding?: unknown,
 *   program?: unknown, outfmt?: unknown, args?: unknown, flow?: unknown, searchContext?: unknown,
 *   queryCanonicalHash?: unknown, subjectCanonicalHash?: unknown, queryProteinSetHash?: unknown,
 *   subjectProteinSetHash?: unknown, queryRuntimeBindingHash?: unknown, subjectRuntimeBindingHash?: unknown,
 *   queryRecordInstanceKey?: unknown, subjectRecordInstanceKey?: unknown, runtime?: unknown,
 *   [field: string]: unknown
 * }} RawLosatCacheEntry
 */

/**
 * A raw cache entry under the key it was found at. `key` is the `cacheKey` given: the caller's `Map` may use any key.
 * @typedef {{ key: unknown, entry: RawLosatCacheEntry }} RawLosatCacheMatch
 */

// The runtime identity the Web stores beside each raw entry it searches
// (design D9/D10, the CLI's losat_runtime_record). It is never part of a raw
// key. The shipped WASM has no build manifest, so its version is not recorded
// (null), as for a CLI runtime whose version probe fails.
/** @param {unknown} program */
export const webLosatRuntimeRecord = (program) => ({
  kind: 'losat',
  source: 'wasm',
  version: null,
  program: String(program || '')
});

// Default raw LOSAT TSV name of a Linear edge; the CLI port is
// gbdraw.comparisons.linear_losat.losat_edge_filename
// (tests/fixtures/losat_edge_filename_cases.json).
/** @param {unknown} label */
const losatFilenameLabel = (label) => String(label ?? '')
  .replace(/[\s/]+/gu, '.')
  .replace(/\.+/g, '.')
  .replace(/^\.+|\.+$/g, '')
  .replace(/[^\p{L}\p{N}_.-]+/gu, '_')
  .replace(/^_+|_+$/g, '');

/**
 * @param {unknown} leftLabel
 * @param {unknown} rightLabel
 * @param {string} suffix
 */
export const losatEdgeFilename = (leftLabel, rightLabel, suffix) => (
  `${losatFilenameLabel(leftLabel) || 'seq_1'}.${losatFilenameLabel(rightLabel) || 'seq_2'}.${suffix}.tsv`
);

/**
 * @param {unknown} value
 * @returns {value is Record<string, unknown>}
 */
export const isPlainObject = (value) => (
  Boolean(value) && typeof value === 'object' && !Array.isArray(value)
);

/**
 * @template T
 * @param {T} value A JSON-serializable value.
 * @returns {T}
 */
const cloneJson = (value) => JSON.parse(JSON.stringify(value));

/**
 * @param {unknown} args
 * @returns {string[]}
 */
export const normalizeLosatArgs = (args) => (
  Array.isArray(args) ? args.map((arg) => String(arg)) : []
);

/**
 * @param {unknown} left
 * @param {unknown} right
 * @returns {boolean}
 */
export const sameLosatArgs = (left, right) => {
  const a = normalizeLosatArgs(left);
  const b = normalizeLosatArgs(right);
  return a.length === b.length && a.every((value, index) => value === b[index]);
};

/**
 * @param {unknown} entry
 * @returns {entry is Record<string, unknown>}
 */
const hasRawShape = (entry) => (
  isPlainObject(entry) &&
  entry.kind === 'raw-losat' &&
  typeof entry.text === 'string' &&
  // `RegExp.test` coerces what it is given to text, so the type is tested first (OV-983).
  (entry.searchContext === undefined
    || (typeof entry.searchContext === 'string' && /^[0-9a-f]{64}$/.test(entry.searchContext)))
);

/**
 * @param {unknown} entry
 * @returns {boolean}
 */
export const isProteinRawLosatCacheEntry = (entry) => (
  hasRawShape(entry) &&
  entry.schema === PROTEIN_LOSAT_CACHE_SCHEMA &&
  entry.identityKind === 'protein' &&
  entry.idEncoding === 'runtime-handle-v1' &&
  String(entry.program || '').toLowerCase() === 'blastp' &&
  typeof entry.queryProteinSetHash === 'string' &&
  Boolean(entry.queryProteinSetHash) &&
  typeof entry.subjectProteinSetHash === 'string' &&
  Boolean(entry.subjectProteinSetHash) &&
  typeof entry.queryRuntimeBindingHash === 'string' &&
  Boolean(entry.queryRuntimeBindingHash) &&
  typeof entry.subjectRuntimeBindingHash === 'string' &&
  Boolean(entry.subjectRuntimeBindingHash) &&
  typeof entry.queryRecordInstanceKey === 'string' &&
  Boolean(entry.queryRecordInstanceKey) &&
  typeof entry.subjectRecordInstanceKey === 'string' &&
  Boolean(entry.subjectRecordInstanceKey)
);

/**
 * @param {unknown} entry
 * @returns {boolean}
 */
export const isLegacyProteinRawLosatCacheEntry = (entry) => (
  hasRawShape(entry) &&
  entry.schema === NUCLEOTIDE_LOSAT_CACHE_SCHEMA &&
  String(entry.program || '').toLowerCase() === 'blastp'
);

/**
 * @param {unknown} entry
 * @returns {boolean}
 */
export const isNucleotideRawLosatCacheEntry = (entry) => {
  if (!hasRawShape(entry) || entry.schema !== NUCLEOTIDE_LOSAT_CACHE_SCHEMA) return false;
  if (entry.identityKind && entry.identityKind !== 'nucleotide') return false;
  return String(entry.program || '').toLowerCase() !== 'blastp';
};

/**
 * @param {unknown} entry
 * @returns {'protein-current' | 'nucleotide-current' | 'protein-legacy' | 'invalid'}
 */
export const classifyRawLosatCacheEntry = (entry) => {
  if (isProteinRawLosatCacheEntry(entry)) return 'protein-current';
  if (isNucleotideRawLosatCacheEntry(entry)) return 'nucleotide-current';
  if (isLegacyProteinRawLosatCacheEntry(entry)) return 'protein-legacy';
  return 'invalid';
};

/**
 * @param {unknown} entry
 * @returns {boolean}
 */
export const isCurrentRawLosatCacheEntry = (entry) => {
  const classification = classifyRawLosatCacheEntry(entry);
  return classification === 'protein-current' || classification === 'nucleotide-current';
};

/**
 * A derived cache entry: what `isLosatDerivedCacheEntry` has checked. `mode` is stored unchecked.
 * @typedef {{ schema: number, kind: 'derived-losatp-payload', key: string, payload: Record<string, unknown>, idEncoding?: unknown, mode?: unknown }} LosatDerivedCacheEntry
 */

/**
 * @param {unknown} entry
 * @param {{ allowLegacy?: boolean }} [options]
 * @returns {entry is LosatDerivedCacheEntry}
 */
export const isLosatDerivedCacheEntry = (entry, { allowLegacy = true } = {}) => (
  isPlainObject(entry) &&
  (
    entry.schema === LOSAT_DERIVED_CACHE_SCHEMA ||
    (
      allowLegacy &&
      entry.schema === LEGACY_LOSAT_DERIVED_CACHE_SCHEMA
    )
  ) &&
  entry.kind === 'derived-losatp-payload' &&
  (
    entry.schema !== LOSAT_DERIVED_CACHE_SCHEMA ||
    entry.idEncoding === 'runtime-handle-v1'
  ) &&
  typeof entry.key === 'string' &&
  Boolean(entry.key) &&
  isPlainObject(entry.payload)
);

/**
 * @param {unknown} value
 * @returns {string}
 */
const stableObjectText = (value) => {
  if (Array.isArray(value)) return `[${value.map(stableObjectText).join(',')}]`;
  if (isPlainObject(value)) {
    return `{${Object.keys(value).sort().map((key) => (
      `${JSON.stringify(key)}:${stableObjectText(value[key])}`
    )).join(',')}}`;
  }
  return JSON.stringify(value);
};

/**
 * @param {Record<string, unknown>} target
 * @param {unknown} incoming
 * @param {string} owner
 */
const mergeManifestMap = (target, incoming, owner) => {
  Object.entries(incoming || {}).forEach(([key, value]) => {
    if (!Object.prototype.hasOwnProperty.call(target, key)) {
      target[key] = cloneJson(value);
      return;
    }
    if (stableObjectText(target[key]) !== stableObjectText(value)) {
      throw new Error(`Protein identity manifest has conflicting ${owner} '${key}'.`);
    }
  });
};

/**
 * @param {unknown} value
 * @param {WeakSet<object>} [visited]
 * @returns {boolean}
 */
const containsManifestPresentationIdentity = (value, visited = new WeakSet()) => {
  if (Array.isArray(value)) {
    return value.some((item) => containsManifestPresentationIdentity(item, visited));
  }
  if (!isPlainObject(value) || visited.has(value)) return false;
  visited.add(value);
  return Object.entries(value).some(([key, item]) => {
    const normalizedKey = String(key)
      .normalize('NFKC')
      .toLowerCase()
      .replace(/[^a-z0-9]+/g, '');
    return MANIFEST_PRESENTATION_IDENTITY_FRAGMENTS.some(
      (fragment) => normalizedKey.includes(fragment)
    ) || containsManifestPresentationIdentity(item, visited);
  });
};

/**
 * @typedef {object} ProteinIdentityInstance
 * @property {number} schema
 * @property {unknown} recordAnalysisId Coerces to an own key of `ProteinIdentityManifest.recordAnalyses`; its type is not checked.
 * @property {string} runtimeBindingHash
 * @property {string} displayBindingHash
 * @property {Record<string, string>} runtimeIds feature analysis ID -> runtime handle
 * @property {Record<string, { displayAlias: string, exportOrdinal?: number | null }>} featureMetadata
 */

/**
 * The protein identity manifest once `validateProteinIdentityManifest` accepts it.
 * @typedef {object} ProteinIdentityManifest
 * @property {number} schema
 * @property {Record<string, { schema: number, proteins: { featureAnalysisId: string }[] }>} proteinSets
 * @property {Record<string, { schema: number, proteinSetHash: string }>} recordAnalyses
 * @property {Record<string, ProteinIdentityInstance>} recordInstances
 */

/** @returns {ProteinIdentityManifest} */
export const emptyProteinIdentityManifest = () => ({
  schema: PROTEIN_IDENTITY_MANIFEST_SCHEMA,
  proteinSets: {},
  recordAnalyses: {},
  recordInstances: {}
});

/**
 * @param {unknown} manifest
 * @returns {manifest is ProteinIdentityManifest}
 */
export const validateProteinIdentityManifest = (manifest) => {
  if (!isPlainObject(manifest) || manifest.schema !== PROTEIN_IDENTITY_MANIFEST_SCHEMA) {
    return false;
  }
  if (containsManifestPresentationIdentity(manifest)) return false;
  if (
    !isPlainObject(manifest.proteinSets) ||
    !isPlainObject(manifest.recordAnalyses) ||
    !isPlainObject(manifest.recordInstances)
  ) return false;

  const proteinSetFeatureIds = new Map();
  for (const [proteinSetHash, proteinSet] of Object.entries(manifest.proteinSets)) {
    if (!proteinSetHash || !isPlainObject(proteinSet) || proteinSet.schema !== 1) return false;
    if (!Array.isArray(proteinSet.proteins)) return false;
    const featureIds = new Set();
    for (const protein of proteinSet.proteins) {
      const featureId = isPlainObject(protein) ? protein.featureAnalysisId : null;
      if (
        typeof featureId !== 'string' ||
        !/^f_[0-9a-f]{64}$/.test(featureId) ||
        featureIds.has(featureId)
      ) return false;
      featureIds.add(featureId);
    }
    proteinSetFeatureIds.set(proteinSetHash, featureIds);
  }

  for (const analysis of Object.values(manifest.recordAnalyses)) {
    if (!isPlainObject(analysis) || analysis.schema !== 1) return false;
    if (
      typeof analysis.proteinSetHash !== 'string' ||
      !proteinSetFeatureIds.has(analysis.proteinSetHash)
    ) return false;
  }

  const allRuntimeHandles = new Set();
  for (const [instanceKey, instance] of Object.entries(manifest.recordInstances)) {
    if (!instanceKey || !isPlainObject(instance) || instance.schema !== 2) return false;
    // `hasOwnProperty.call` and the lookup below coerce a non-string `recordAnalysisId` to a key.
    if (!Object.prototype.hasOwnProperty.call(manifest.recordAnalyses, /** @type {string} */ (instance.recordAnalysisId))) {
      return false;
    }
    if (
      typeof instance.runtimeBindingHash !== 'string' ||
      !instance.runtimeBindingHash ||
      typeof instance.displayBindingHash !== 'string' ||
      !instance.displayBindingHash ||
      !isPlainObject(instance.runtimeIds) ||
      !isPlainObject(instance.featureMetadata)
    ) return false;
    // The loop over `recordAnalyses` above proved each one is a plain object with a string `proteinSetHash`.
    const analysis = /** @type {{ proteinSetHash: string }} */ (
      manifest.recordAnalyses[/** @type {string} */ (instance.recordAnalysisId)]
    );
    const expectedFeatureIds = proteinSetFeatureIds.get(analysis.proteinSetHash);
    if (
      Object.keys(instance.runtimeIds).length !== expectedFeatureIds.size ||
      Object.keys(instance.featureMetadata).length !== expectedFeatureIds.size
    ) return false;

    const featuresByAlias = new Map();
    for (const [featureId, runtimeHandle] of Object.entries(instance.runtimeIds)) {
      if (!expectedFeatureIds.has(featureId)) return false;
      if (
        typeof runtimeHandle !== 'string' ||
        !/^h_[a-z2-7]{26}$/.test(runtimeHandle) ||
        allRuntimeHandles.has(runtimeHandle)
      ) return false;
      const metadata = instance.featureMetadata[featureId];
      if (
        !isPlainObject(metadata) ||
        typeof metadata.displayAlias !== 'string'
      ) return false;
      const alias = metadata.displayAlias.normalize('NFC').trim();
      if (!alias) return false;
      if (!featuresByAlias.has(alias)) featuresByAlias.set(alias, []);
      featuresByAlias.get(alias).push(featureId);
      allRuntimeHandles.add(runtimeHandle);
    }
    if (Object.keys(instance.featureMetadata).some((featureId) => !expectedFeatureIds.has(featureId))) {
      return false;
    }
    for (const featureIds of featuresByAlias.values()) {
      const ordered = [...featureIds].sort();
      for (let index = 0; index < ordered.length; index += 1) {
        const featureId = ordered[index];
        const expectedOrdinal = ordered.length > 1 ? index + 1 : null;
        // Every `featureId` here is a key of `runtimeIds`, whose metadata the loop above proved is a plain object.
        const rawOrdinal = /** @type {{ exportOrdinal?: unknown }} */ (instance.featureMetadata[featureId]).exportOrdinal;
        const actualOrdinal = rawOrdinal === undefined || rawOrdinal === null
          ? null
          : rawOrdinal;
        if (actualOrdinal !== expectedOrdinal) return false;
      }
    }
  }
  return true;
};

/** @type {WeakMap<object, { manifest: ProteinIdentityManifest, runtimeIdsByInstance: Map<string, Set<string>>, runtimeHandles: Set<string> | null }>} */
const validatedProteinIdentityIndexes = new WeakMap();

/**
 * @param {unknown} manifest
 * @returns {object | null}
 */
export const buildValidatedProteinIdentityIndex = (manifest) => {
  if (!validateProteinIdentityManifest(manifest)) return null;
  const runtimeIdsByInstance = new Map();
  Object.entries(manifest.recordInstances).forEach(([instanceKey, instance]) => {
    const runtimeIds = new Set(Object.values(instance.runtimeIds));
    runtimeIdsByInstance.set(instanceKey, runtimeIds);
  });
  const index = Object.freeze({});
  validatedProteinIdentityIndexes.set(index, {
    manifest,
    runtimeIdsByInstance,
    runtimeHandles: null
  });
  return index;
};

/** @param {object | null | undefined} index */
export const releaseValidatedProteinIdentityIndex = (index) => (
  // `WeakMap.delete` returns false for a non-object key, so a missing index is harmless.
  validatedProteinIdentityIndexes.delete(/** @type {object} */ (index))
);

/**
 * @param {unknown} manifest
 * @param {object | null | undefined} index
 */
const indexedProteinIdentity = (manifest, index) => {
  const validated = index && validatedProteinIdentityIndexes.get(index);
  return validated?.manifest === manifest ? validated : null;
};

/**
 * @param {unknown} manifests
 * @param {{ invalidInputMessage?: string | null }} [options]
 */
export const mergeProteinIdentityManifests = (manifests, { invalidInputMessage = null } = {}) => {
  const inputs = Array.isArray(manifests) ? manifests : [];
  // A caller-supplied input error takes precedence over any merge conflict.
  // Keep the default helper's interleaved failure order for other callers.
  const inputsValidated = invalidInputMessage !== null
    && inputs.every((manifest) => validateProteinIdentityManifest(manifest));
  if (invalidInputMessage !== null && !inputsValidated) {
    throw new Error(invalidInputMessage);
  }
  const merged = emptyProteinIdentityManifest();
  inputs.forEach((manifest) => {
    if (!inputsValidated && !validateProteinIdentityManifest(manifest)) {
      throw new Error('Cannot merge an invalid protein identity manifest.');
    }
    mergeManifestMap(merged.proteinSets, manifest.proteinSets, 'protein set');
    mergeManifestMap(merged.recordAnalyses, manifest.recordAnalyses, 'record analysis');
    mergeManifestMap(merged.recordInstances, manifest.recordInstances, 'record instance');
  });
  if (!validateProteinIdentityManifest(merged)) {
    throw new Error('Merged protein identity manifest is invalid.');
  }
  return merged;
};

/**
 * @param {any} entry The cache entry to check.
 * @param {any} manifest The protein identity manifest it must reference.
 * @param {{ identityIndex?: object | null }} [options] An index from `buildValidatedProteinIdentityIndex` for `manifest`.
 */
export const validateProteinRawEntryReferences = (
  entry,
  manifest,
  { identityIndex = null } = {}
) => {
  if (!isProteinRawLosatCacheEntry(entry)) {
    return false;
  }
  const indexedIdentity = indexedProteinIdentity(manifest, identityIndex);
  if (!indexedIdentity && !validateProteinIdentityManifest(manifest)) return false;
  const queryInstance = manifest.recordInstances[entry.queryRecordInstanceKey];
  const subjectInstance = manifest.recordInstances[entry.subjectRecordInstanceKey];
  if (!queryInstance || !subjectInstance) return false;
  if (queryInstance.runtimeBindingHash !== entry.queryRuntimeBindingHash) return false;
  if (subjectInstance.runtimeBindingHash !== entry.subjectRuntimeBindingHash) return false;
  const queryAnalysis = manifest.recordAnalyses[queryInstance.recordAnalysisId];
  const subjectAnalysis = manifest.recordAnalyses[subjectInstance.recordAnalysisId];
  if (
    queryAnalysis?.proteinSetHash !== entry.queryProteinSetHash ||
    subjectAnalysis?.proteinSetHash !== entry.subjectProteinSetHash
  ) return false;
  return rawProteinTextMatchesBindings(
    entry.text,
    indexedIdentity?.runtimeIdsByInstance.get(entry.queryRecordInstanceKey)
      || new Set(Object.values(queryInstance.runtimeIds)),
    indexedIdentity?.runtimeIdsByInstance.get(entry.subjectRecordInstanceKey)
      || new Set(Object.values(subjectInstance.runtimeIds))
  );
};

/**
 * @param {unknown} text
 * @param {unknown} queryIds
 * @param {unknown} subjectIds
 * @returns {boolean}
 */
export const rawProteinTextMatchesBindings = (text, queryIds, subjectIds) => {
  if (typeof text !== 'string' || !(queryIds instanceof Set) || !(subjectIds instanceof Set)) {
    return false;
  }
  for (const rawLine of text.split(/\r?\n/)) {
    const line = rawLine.trim();
    if (!line || line.startsWith('#')) continue;
    const columns = rawLine.split('\t');
    const hasInvalidNumericField = columns.slice(2).some((value, offset) => {
      const columnIndex = offset + 2;
      const pattern = LOSAT_OUTFMT6_INTEGER_COLUMN_INDEXES.has(columnIndex)
        ? STRICT_DECIMAL_INTEGER_RE
        : STRICT_DECIMAL_NUMBER_RE;
      return !pattern.test(value) || !Number.isFinite(Number(value));
    });
    if (
      columns.length !== LOSAT_OUTFMT6_COLUMN_COUNT ||
      !queryIds.has(columns[0]) ||
      !subjectIds.has(columns[1]) ||
      hasInvalidNumericField
    ) {
      return false;
    }
  }
  return true;
};

const DERIVED_SCALAR_PROTEIN_REFERENCE_KEYS = new Set([
  'proteinId',
  'queryProteinId',
  'subjectProteinId',
  'protein_id',
  'query_protein_id',
  'subject_protein_id'
]);
const DERIVED_UNIT_PROTEIN_REFERENCE_KEYS = new Set([
  'queryUnitId',
  'subjectUnitId',
  'query_unit_id',
  'subject_unit_id'
]);
const DERIVED_ARRAY_PROTEIN_REFERENCE_KEYS = new Set([
  'proteinIds',
  'sharedProteinIds',
  'protein_ids',
  'shared_protein_ids'
]);
const DERIVED_COMPOUND_EDGE_REFERENCE_KEYS = new Set([
  'supportingEdge',
  'supportingEdges',
  'supporting_edge',
  'supporting_edges',
  'edgeId',
  'edgeIds',
  'edge_id',
  'edge_ids'
]);
const RUNTIME_HANDLE_RE = /^h_[a-z2-7]{26}$/;
const FEATURE_ANALYSIS_ID_RE = /(?:^|[^A-Za-z0-9_])f_[0-9a-f]{64}(?:$|[^A-Za-z0-9_])/;
const LEGACY_PROTEIN_REFERENCE_RE = /(?:^|[^A-Za-z0-9._%+-])p_[A-Za-z0-9._%+-]+?_\d+_\d+_(?:-1|0|1)_[0-9a-f]{12}(?:_[2-9][0-9]*)?(?:$|[^A-Za-z0-9._%+-])/;
// Reject this unsupported historical shape defensively at the display boundary.
const LONG_TRANSPORT_ID_RE = /(?:^|[^A-Za-z0-9._%-])(?:[A-Za-z0-9._-]|%[0-9A-F]{2})+@(?:[A-Za-z0-9._-]|%[0-9A-F]{2})+\|(?:[A-Za-z0-9._-]|%[0-9A-F]{2})+~f_[0-9a-f]{64}(?:$|[^A-Za-z0-9._%-])/;
const COMPOUND_SUPPORTING_EDGE_RE = /^(h_[a-z2-7]{26})->(h_[a-z2-7]{26}):[A-Za-z][A-Za-z0-9._-]*$/;
const COMPOUND_PATH_EDGE_RE = /^[^:\s]+:\d+:(h_[a-z2-7]{26})->\d+:(h_[a-z2-7]{26}):[A-Za-z][A-Za-z0-9._-]*$/;

/**
 * @param {unknown} value
 * @returns {boolean}
 */
// `Number.isInteger` is true only for a number, which is what the cast to `number` states.
const isNonnegativeInteger = (value) => Number.isInteger(value) && /** @type {number} */ (value) >= 0;

/**
 * @param {LosatDerivedCacheEntry} entry A derived entry that `isLosatDerivedCacheEntry` accepted.
 * @returns {boolean}
 */
const isStrictEmptyDerivedResult = (entry) => {
  const mode = entry?.mode;
  if (!(/** @type {readonly unknown[]} */ (['orthogroup', 'collinear'])).includes(mode) || !isPlainObject(entry?.payload)) {
    return false;
  }
  const payload = entry.payload;
  const allowedKeys = new Set(['identity', 'provenance', 'pairs', 'orthogroups']);
  if (mode === 'collinear') {
    [
      'collinearGroups',
      'collinearGroupScope',
      'collinearityBlocks'
    ].forEach((key) => allowedKeys.add(key));
  }
  if (Object.keys(payload).some((key) => !allowedKeys.has(key))) return false;

  if (Object.prototype.hasOwnProperty.call(payload, 'identity')) {
    const identity = payload.identity;
    if (
      !isPlainObject(identity) ||
      identity.cacheSchema !== LOSAT_DERIVED_CACHE_SCHEMA ||
      identity.idEncoding !== 'runtime-handle-v1' ||
      identity.mode !== mode ||
      !Array.isArray(identity.rawCacheKeys) ||
      identity.rawCacheKeys.some((key) => typeof key !== 'string' || !key)
    ) return false;
  }

  if (
    !Array.isArray(payload.pairs) ||
    payload.pairs.length === 0 ||
    !Array.isArray(payload.orthogroups)
  ) {
    return false;
  }
  if (payload.orthogroups.length !== 0) return false;
  const allowedPairKeys = new Set([
    'pair_index',
    'query_index',
    'subject_index',
    'tsv',
    'rows',
    'hit_count'
  ]);
  const pairIndices = new Set();
  for (const pair of payload.pairs) {
    if (
      !isPlainObject(pair) ||
      Object.keys(pair).some((key) => !allowedPairKeys.has(key)) ||
      !isNonnegativeInteger(pair.pair_index) ||
      pairIndices.has(pair.pair_index) ||
      pair.tsv !== '' ||
      !Array.isArray(pair.rows) ||
      pair.rows.length !== 0 ||
      pair.hit_count !== 0
    ) return false;
    pairIndices.add(pair.pair_index);
    const hasQueryIndex = Object.prototype.hasOwnProperty.call(pair, 'query_index');
    const hasSubjectIndex = Object.prototype.hasOwnProperty.call(pair, 'subject_index');
    if (
      hasQueryIndex !== hasSubjectIndex ||
      (
        hasQueryIndex &&
        (
          !isNonnegativeInteger(pair.query_index) ||
          !isNonnegativeInteger(pair.subject_index)
        )
      )
    ) return false;
  }

  if (mode === 'collinear') {
    for (const key of ['collinearGroups', 'collinearityBlocks']) {
      if (
        Object.prototype.hasOwnProperty.call(payload, key) &&
        (!Array.isArray(payload[key]) || payload[key].length !== 0)
      ) return false;
    }
    if (
      Object.prototype.hasOwnProperty.call(payload, 'collinearGroupScope') &&
      !(/** @type {readonly unknown[]} */ (['adjacent_local', 'global_collinear'])).includes(payload.collinearGroupScope)
    ) return false;
  }
  return true;
};

/**
 * @param {any} entry The cache entry to check.
 * @param {any} manifest The protein identity manifest it must reference.
 * @param {{ identityIndex?: object | null }} [options] An index from `buildValidatedProteinIdentityIndex` for `manifest`.
 */
export const validateDerivedProteinReferences = (
  entry,
  manifest,
  { identityIndex = null } = {}
) => {
  if (!isLosatDerivedCacheEntry(entry, { allowLegacy: false })) return false;
  const indexedIdentity = indexedProteinIdentity(manifest, identityIndex);
  if (!indexedIdentity && !validateProteinIdentityManifest(manifest)) return false;
  let runtimeHandles = indexedIdentity?.runtimeHandles || null;
  if (!runtimeHandles) {
    runtimeHandles = new Set();
    Object.values(manifest.recordInstances).forEach((instance) => {
      Object.values(instance.runtimeIds || {}).forEach((handle) => /** @type {Set<string>} */ (runtimeHandles).add(handle)); // set just above, when it was missing
    });
    if (indexedIdentity) indexedIdentity.runtimeHandles = runtimeHandles;
  }
  /** @param {string} value */
  const hasForbiddenLegacyReference = (value) => (
    LEGACY_PROTEIN_REFERENCE_RE.test(value) ||
    LONG_TRANSPORT_ID_RE.test(value) ||
    FEATURE_ANALYSIS_ID_RE.test(value)
  );
  /** @param {string} value */
  const compoundEdgeReferences = (value) => {
    const match = COMPOUND_SUPPORTING_EDGE_RE.exec(value) || COMPOUND_PATH_EDGE_RE.exec(value);
    return match ? [match[1], match[2]] : null;
  };
  let sawReference = false;
  /**
   * @param {unknown} value
   * @param {string} [ownerKey]
   * @returns {boolean}
   */
  const visit = (value, ownerKey = '') => {
    if (
      DERIVED_COMPOUND_EDGE_REFERENCE_KEYS.has(ownerKey) &&
      typeof value !== 'string' &&
      !Array.isArray(value)
    ) return false;
    if (typeof value === 'string') {
      if (hasForbiddenLegacyReference(value)) return false;
      if (DERIVED_COMPOUND_EDGE_REFERENCE_KEYS.has(ownerKey)) {
        const references = compoundEdgeReferences(value);
        if (!references || references.some((reference) => !runtimeHandles.has(reference))) {
          return false;
        }
        sawReference = true;
        return true;
      }
      if (DERIVED_SCALAR_PROTEIN_REFERENCE_KEYS.has(ownerKey) && value) {
        sawReference = true;
        const references = value.split(';').map((reference) => reference.trim());
        return (
          references.length > 0 &&
          references.every((reference) => reference && runtimeHandles.has(reference))
        );
      }
      if (DERIVED_UNIT_PROTEIN_REFERENCE_KEYS.has(ownerKey) && value) {
        const references = value.split(';').map((reference) => reference.trim());
        const runtimeReferences = references.filter((reference) => reference.startsWith('h_'));
        if (runtimeReferences.length === 0) return true;
        sawReference = true;
        return (
          references.every(Boolean) &&
          runtimeReferences.every(
            (reference) => RUNTIME_HANDLE_RE.test(reference) && runtimeHandles.has(reference)
          )
        );
      }
      return true;
    }
    if (Array.isArray(value)) {
      if (DERIVED_COMPOUND_EDGE_REFERENCE_KEYS.has(ownerKey)) {
        return value.every((item) => typeof item === 'string' && visit(item, ownerKey));
      }
      if (DERIVED_ARRAY_PROTEIN_REFERENCE_KEYS.has(ownerKey)) {
        for (const item of value) {
          if (typeof item !== 'string' || !runtimeHandles.has(item)) return false;
          sawReference = true;
        }
        return true;
      }
      return value.every((item) => visit(item, ''));
    }
    if (!isPlainObject(value)) return true;
    return Object.entries(value).every(([key, item]) => {
      if (hasForbiddenLegacyReference(key)) return false;
      if (RUNTIME_HANDLE_RE.test(key)) {
        sawReference = true;
        if (!runtimeHandles.has(key)) return false;
      }
      return visit(item, key);
    });
  };
  return (
    visit(entry.payload) &&
    (
      sawReference ||
      Object.keys(entry.payload).length === 0 ||
      isStrictEmptyDerivedResult(entry)
    )
  );
};

/**
 * @param {unknown} cacheMap
 * @param {unknown} cacheKey
 * @param {Record<string, any>} [metadata]
 * @param {Record<string, any> | null} [manifest]
 * @param {{ identityIndex?: object | null }} [options]
 * @returns {RawLosatCacheMatch | null}
 */
export const getCurrentRawLosatCacheEntry = (
  cacheMap, cacheKey, metadata = {}, manifest = null, { identityIndex = null } = {}
) => {
  if (!(cacheMap instanceof Map) || !cacheKey) return null;
  const entry = cacheMap.get(cacheKey);
  if ((entry?.searchContext ?? null) !== (metadata.searchContext ?? null)) return null;
  const classification = classifyRawLosatCacheEntry(entry);
  if (classification === 'protein-current') {
    if (String(entry.program || '') !== String(metadata.program || 'blastp')) return null;
    if (String(entry.outfmt || '6') !== String(metadata.outfmt || '6')) return null;
    if (!sameLosatArgs(entry.args, metadata.args)) return null;
    if (
      metadata.queryRuntimeBindingHash &&
      entry.queryRuntimeBindingHash !== metadata.queryRuntimeBindingHash
    ) return null;
    if (
      metadata.subjectRuntimeBindingHash &&
      entry.subjectRuntimeBindingHash !== metadata.subjectRuntimeBindingHash
    ) return null;
    if (!validateProteinRawEntryReferences(entry, manifest, { identityIndex })) return null;
    return { key: cacheKey, entry };
  }
  if (classification !== 'nucleotide-current') return null;
  if (String(entry.program || '') !== String(metadata.program || '')) return null;
  if (String(entry.outfmt || '6') !== String(metadata.outfmt || '6')) return null;
  if (!sameLosatArgs(entry.args, metadata.args)) return null;
  if (String(entry.queryCanonicalHash || '') !== String(metadata.queryCanonicalHash || '')) return null;
  if (String(entry.subjectCanonicalHash || '') !== String(metadata.subjectCanonicalHash || '')) return null;
  if (String(entry.flow || '') !== String(metadata.flow || '')) return null;
  return { key: cacheKey, entry };
};

/**
 * A legacy protein cache entry held back until its references validate.
 * @typedef {object} LegacyProteinCandidate
 * @property {string} state 'pending', 'promoted', or 'rejected'
 * @property {Record<string, any>} originalEntry the raw cache entry, as imported
 * @property {string | null} rejectionReason
 */

/**
 * @typedef {object} LegacyProteinCandidateEnvelope
 * @property {number} schema
 * @property {LegacyProteinCandidate[]} entries
 */

/**
 * @param {Record<string, any>[]} entries the raw cache entries of an imported Session
 * @returns {LegacyProteinCandidateEnvelope}
 */
export const createLegacyProteinCandidateEnvelope = (entries) => ({
  schema: LEGACY_PROTEIN_CANDIDATE_SCHEMA,
  entries: (Array.isArray(entries) ? entries : [])
    .filter(isLegacyProteinRawLosatCacheEntry)
    .map((entry) => ({
      state: 'pending',
      originalEntry: cloneJson(entry),
      rejectionReason: null
    }))
});

/** @param {unknown} value */
export const normalizeLegacyProteinCandidateEnvelope = (value) => {
  if (!isPlainObject(value) || value.schema !== LEGACY_PROTEIN_CANDIDATE_SCHEMA) {
    return createLegacyProteinCandidateEnvelope([]);
  }
  const entries = (Array.isArray(value.entries) ? value.entries : [])
    .filter((candidate) => (
      isPlainObject(candidate) &&
      LEGACY_CANDIDATE_STATES.has(candidate.state) &&
      isLegacyProteinRawLosatCacheEntry(candidate.originalEntry)
    ))
    .map((candidate) => ({
      state: candidate.state,
      originalEntry: cloneJson(candidate.originalEntry),
      rejectionReason: candidate.rejectionReason == null
        ? null
        : String(candidate.rejectionReason)
    }));
  return { schema: LEGACY_PROTEIN_CANDIDATE_SCHEMA, entries };
};

/**
 * @param {unknown} envelope
 * @param {number} candidateIndex
 * @param {string} nextState
 * @param {string | null} [rejectionReason]
 */
export const transitionLegacyProteinCandidate = (
  envelope,
  candidateIndex,
  nextState,
  rejectionReason = null
) => {
  const normalized = normalizeLegacyProteinCandidateEnvelope(envelope);
  if (!LEGACY_CANDIDATE_STATES.has(nextState)) {
    throw new Error(`Unsupported legacy protein candidate state: ${nextState}`);
  }
  if (!Number.isInteger(candidateIndex) || candidateIndex < 0 || candidateIndex >= normalized.entries.length) {
    throw new Error('Legacy protein candidate index is out of range.');
  }
  const entries = normalized.entries.map((candidate, index) => (
    index === candidateIndex
      ? {
          ...candidate,
          state: nextState,
          rejectionReason: nextState === 'rejected'
            ? String(rejectionReason || 'Legacy cache validation failed.')
            : null
        }
      : candidate
  ));
  return { schema: LEGACY_PROTEIN_CANDIDATE_SCHEMA, entries };
};

/** @param {unknown} envelope */
export const serializableLegacyProteinCandidateEnvelope = (envelope) => {
  const normalized = normalizeLegacyProteinCandidateEnvelope(envelope);
  return {
    schema: LEGACY_PROTEIN_CANDIDATE_SCHEMA,
    entries: normalized.entries.filter((candidate) => candidate.state !== 'promoted')
  };
};

/**
 * @param {unknown} envelope
 * @param {{ states?: string[] }} [options]
 */
export const candidateOriginalEntries = (envelope, { states = ['pending'] } = {}) => {
  const acceptedStates = new Set(states);
  return normalizeLegacyProteinCandidateEnvelope(envelope).entries
    .map((candidate, index) => ({ candidate, index }))
    .filter(({ candidate }) => acceptedStates.has(candidate.state))
    .map(({ candidate, index }) => ({ index, entry: cloneJson(candidate.originalEntry) }));
};
