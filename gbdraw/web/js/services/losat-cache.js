// @ts-check
import {
  LEGACY_LOSAT_DERIVED_CACHE_SCHEMA,
  LEGACY_PROTEIN_CANDIDATE_SCHEMA,
  LOSAT_DERIVED_CACHE_SCHEMA,
  NUCLEOTIDE_LOSAT_CACHE_SCHEMA,
  PROTEIN_IDENTITY_MANIFEST_SCHEMA,
  PROTEIN_LOSAT_CACHE_SCHEMA
} from '../contracts/schema.generated.js';
import { LOSAT_WASM_BUILD } from './losat-build.generated.js';

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
// key. Its version is the shipped Wasm's build (gbdraw/web/wasm/losat/losat-build.json).
/** @param {unknown} program */
export const webLosatRuntimeRecord = (program) => ({
  kind: 'losat',
  source: 'wasm',
  version: LOSAT_WASM_BUILD.version,
  program: String(program || '')
});

/**
 * The LOSATN task a search's options name (`--task`, the raw-cache spelling; Generate always writes it), or null.
 * @param {unknown} args
 * @returns {string | null}
 */
export const losatSearchTask = (args) => {
  const words = Array.isArray(args) ? args.map(String) : [];
  const index = words.indexOf('--task');
  return index >= 0 && index + 1 < words.length ? words[index + 1] : null;
};

/**
 * Whether Generate may reuse a raw entry's rows, by the runtime that searched them. A CLI search and this LOSAT
 * Wasm build are always reused. Rows of the Wasm of gbdraw 0.14.4 and earlier (`version: null`) and rows with no
 * runtime record (Sessions before version 44, promoted legacy protein caches) are reused only for the searches
 * that Wasm gives the same rows as this build: LOSATP, and LOSATN task megablast (#1540's engine parity). Its
 * LOSATN task blastn gave other E-values and hits, its task dc-megablast ran blastn (OV-2377), and its TLOSATX
 * rows came in another order, so those are searched again once. So are rows of any other Wasm build.
 * @param {unknown} entry
 * @returns {boolean}
 */
export const isReusableLosatEntry = (entry) => {
  if (!isPlainObject(entry)) return false;
  const runtime = isPlainObject(entry.runtime) ? entry.runtime : null;
  if (runtime !== null && (runtime.source !== 'wasm' || runtime.version === LOSAT_WASM_BUILD.version)) return true;
  if (runtime !== null && runtime.version !== null) return false;
  return entry.program === 'blastp' || (entry.program === 'blastn' && losatSearchTask(entry.args) === 'megablast');
};

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
 * @param {unknown} incoming The map of an admitted manifest: its values are frozen, so the merge shares them.
 * @param {string} owner
 */
const mergeManifestMap = (target, incoming, owner) => {
  Object.entries(incoming || {}).forEach(([key, value]) => {
    if (!Object.prototype.hasOwnProperty.call(target, key)) {
      target[key] = value;
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

// The manifests the merge admitted (OV-1553): each passed `validateProteinIdentityManifest` once and is deeply
// frozen, so no reader can change it and its checks, its runtime ID sets and the raw tables verified against it
// stay valid for the life of the object.
/** @type {WeakSet<object>} */
const admittedManifests = new WeakSet();

/**
 * @param {unknown} manifest
 * @returns {manifest is ProteinIdentityManifest}
 */
const isAdmittedProteinIdentityManifest = (manifest) => (
  typeof manifest === 'object' && manifest !== null && admittedManifests.has(manifest)
);

/**
 * @param {unknown} value A JSON value.
 * @returns {void}
 */
const deepFreeze = (value) => {
  if (!value || typeof value !== 'object') return;
  Object.freeze(value);
  Object.values(value).forEach(deepFreeze);
};

/**
 * Validates `manifest` once and freezes it. The caller gives up the right to change it.
 * @param {unknown} manifest
 * @returns {manifest is ProteinIdentityManifest}
 */
const admitProteinIdentityManifest = (manifest) => {
  if (isAdmittedProteinIdentityManifest(manifest)) return true;
  if (!validateProteinIdentityManifest(manifest)) return false;
  deepFreeze(manifest);
  admittedManifests.add(manifest);
  return true;
};

/**
 * The runtime handles of one record instance of an admitted manifest, one Set per frozen instance object, so a
 * merge that shares the instance shares its Set.
 * @type {WeakMap<object, Set<string>>}
 */
const admittedInstanceRuntimeIds = new WeakMap();

/** @param {ProteinIdentityInstance} instance A record instance of an admitted manifest. */
const admittedRuntimeIdsOf = (instance) => {
  let runtimeIds = admittedInstanceRuntimeIds.get(instance);
  if (!runtimeIds) {
    runtimeIds = new Set(Object.values(instance.runtimeIds));
    admittedInstanceRuntimeIds.set(instance, runtimeIds);
  }
  return runtimeIds;
};

/**
 * @param {unknown} manifest
 * @returns {manifest is ProteinIdentityManifest}
 */
export const validateProteinIdentityManifest = (manifest) => {
  // An admitted manifest passed this check and is frozen, so the answer cannot have changed.
  if (isAdmittedProteinIdentityManifest(manifest)) return true;
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
 * The validated identity of an admitted manifest: built once per frozen manifest object.
 * @type {WeakMap<object, { manifest: ProteinIdentityManifest, runtimeIdsByInstance: Map<string, Set<string>>, runtimeHandles: Set<string> | null }>}
 */
const admittedIdentities = new WeakMap();

/**
 * The checked identity of `manifest`: the one `index` names, or the one of an admitted manifest.
 * @param {unknown} manifest
 * @param {object | null | undefined} index
 */
const indexedProteinIdentity = (manifest, index) => {
  const validated = index && validatedProteinIdentityIndexes.get(index);
  if (validated?.manifest === manifest) return validated;
  if (!isAdmittedProteinIdentityManifest(manifest)) return null;
  let identity = admittedIdentities.get(manifest);
  if (!identity) {
    identity = {
      manifest,
      runtimeIdsByInstance: new Map(Object.entries(manifest.recordInstances).map(
        ([instanceKey, instance]) => [instanceKey, admittedRuntimeIdsOf(instance)]
      )),
      runtimeHandles: null
    };
    admittedIdentities.set(manifest, identity);
  }
  return identity;
};

/**
 * The raw table each protein raw entry was last verified against (OV-1554): its text and the two runtime ID Sets
 * of an admitted manifest. The scan's answer depends on nothing else, and an admitted manifest's Sets never change.
 * @type {WeakMap<object, { text: string, queryIds: Set<string>, subjectIds: Set<string> }>}
 */
const verifiedRawProteinTexts = new WeakMap();

/** @type {{ inputs: unknown[], merged: ProteinIdentityManifest } | null} */
let lastMerge = null;

/**
 * Merges the manifests of the compared records into one admitted (frozen) manifest. Each input is admitted:
 * validated once in its life and frozen, so the caller passes manifests it owns (the extraction's) and does not
 * change them afterwards. The merged manifest shares the inputs' frozen parts. The same inputs, in the same
 * order, return the manifest of the previous merge (OV-1553); a failed merge keeps nothing.
 * @param {unknown} manifests
 * @param {{ invalidInputMessage?: string | null }} [options]
 * @returns {ProteinIdentityManifest}
 */
export const mergeProteinIdentityManifests = (manifests, { invalidInputMessage = null } = {}) => {
  const inputs = Array.isArray(manifests) ? manifests : [];
  if (
    lastMerge
    && lastMerge.inputs.length === inputs.length
    && lastMerge.inputs.every((input, index) => input === inputs[index])
  ) return lastMerge.merged;
  // A caller-supplied input error takes precedence over any merge conflict.
  // Keep the default helper's interleaved failure order for other callers.
  const inputsValidated = invalidInputMessage !== null
    && inputs.every((manifest) => admitProteinIdentityManifest(manifest));
  if (invalidInputMessage !== null && !inputsValidated) {
    throw new Error(invalidInputMessage);
  }
  const merged = emptyProteinIdentityManifest();
  inputs.forEach((manifest) => {
    if (!inputsValidated && !admitProteinIdentityManifest(manifest)) {
      throw new Error('Cannot merge an invalid protein identity manifest.');
    }
    mergeManifestMap(merged.proteinSets, manifest.proteinSets, 'protein set');
    mergeManifestMap(merged.recordAnalyses, manifest.recordAnalyses, 'record analysis');
    mergeManifestMap(merged.recordInstances, manifest.recordInstances, 'record instance');
  });
  if (!validateProteinIdentityManifest(merged)) {
    throw new Error('Merged protein identity manifest is invalid.');
  }
  // The values are the inputs' frozen ones; freezing the maps completes the merged manifest.
  [merged.proteinSets, merged.recordAnalyses, merged.recordInstances, merged].forEach((value) => Object.freeze(value));
  admittedManifests.add(merged);
  lastMerge = { inputs: [...inputs], merged };
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
  const queryIds = indexedIdentity?.runtimeIdsByInstance.get(entry.queryRecordInstanceKey)
    || new Set(Object.values(queryInstance.runtimeIds));
  const subjectIds = indexedIdentity?.runtimeIdsByInstance.get(entry.subjectRecordInstanceKey)
    || new Set(Object.values(subjectInstance.runtimeIds));
  const verified = verifiedRawProteinTexts.get(entry);
  if (
    verified
    && verified.text === entry.text
    && verified.queryIds === queryIds
    && verified.subjectIds === subjectIds
  ) return true;
  if (!rawProteinTextMatchesBindings(entry.text, queryIds, subjectIds)) return false;
  if (isAdmittedProteinIdentityManifest(manifest)) {
    verifiedRawProteinTexts.set(entry, { text: entry.text, queryIds, subjectIds });
  }
  return true;
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
 * A zero-hit pair holds exactly one empty table: the header-only `canonicalTsv` the current writer emits, or the
 * `rows: []` that derived entries carried from #312 until the canonical table replaced it. The Python twin is
 * `_has_one_empty_table_form` in `gbdraw/analysis/protein_artifacts.py`.
 * @param {Record<string, unknown>} pair
 */
const hasOneEmptyTableForm = (pair) => {
  const hasTable = Object.prototype.hasOwnProperty.call(pair, 'canonicalTsv');
  if (hasTable === Object.prototype.hasOwnProperty.call(pair, 'rows')) return false;
  if (!hasTable) return Array.isArray(pair.rows) && pair.rows.length === 0;
  return typeof pair.canonicalTsv === 'string' && /^[^\n]*\n$/.test(pair.canonicalTsv);
};

// The typed analysis result a derived payload of each mode carries: its payload member, the kind of its typed
// resource and the type of its value. The Python twin is `_DERIVED_RESULT_MEMBERS` in
// `gbdraw/analysis/protein_artifacts.py`.
/** @type {Readonly<Record<string, { member: string, kind: string, type: string }>>} */
const DERIVED_RESULT_MEMBERS = Object.freeze({
  orthogroup: { member: 'orthogroupResult', kind: 'orthogroupResult', type: 'OrthogroupGraphResult' },
  collinear: { member: 'collinearityResult', kind: 'result', type: 'CollinearityResult' }
});

/**
 * Whether `resource` is a typed analysis result of `mode` as the Python typed codec writes it: schema 3, the mode's
 * kind and value type, and an object of fields.
 * @param {unknown} mode
 * @param {unknown} resource
 * @returns {boolean}
 */
export const isDerivedAnalysisResult = (mode, resource) => {
  const expected = typeof mode === 'string' && Object.hasOwn(DERIVED_RESULT_MEMBERS, mode)
    ? DERIVED_RESULT_MEMBERS[mode]
    : null;
  if (!expected || !isPlainObject(resource)) return false;
  const value = resource.value;
  return resource.schema === 3
    && resource.kind === expected.kind
    && isPlainObject(value)
    && value.type === expected.type
    && isPlainObject(value.fields);
};

/**
 * Whether a JSON value names nothing: every array and object in it is empty. A typed object (`{ type, fields }`) is
 * read through its fields, so a typed result whose collections are all empty, nested results included, holds no
 * value. The Python twin is `_holds_no_value`.
 * @param {unknown} value
 * @returns {boolean}
 */
const holdsNoValue = (value) => {
  if (Array.isArray(value)) return value.length === 0;
  if (!isPlainObject(value)) return true;
  const keys = Object.keys(value);
  if (keys.length === 2 && Object.hasOwn(value, 'type') && isPlainObject(value.fields)) {
    return Object.values(value.fields).every(holdsNoValue);
  }
  return keys.length === 0;
};

/**
 * Whether a derived entry is the result of zero-hit searches, which names no protein (OV-1127): its mode computes a
 * result, its identity (if present) is current, every pair is a zero-hit pair, its typed result (if present) is the
 * mode's, and every other member holds no value: the typed result and the members of the form written before it
 * (`orthogroups: []`) are empty. The Python twin is `_is_zero_hit_result` in `gbdraw/analysis/protein_artifacts.py`;
 * `tests/fixtures/session_admission_vectors.json` pins both.
 * @param {LosatDerivedCacheEntry} entry A derived entry that `isLosatDerivedCacheEntry` accepted.
 * @returns {boolean}
 */
const isZeroHitDerivedResult = (entry) => {
  const mode = entry?.mode;
  const result = typeof mode === 'string' && Object.hasOwn(DERIVED_RESULT_MEMBERS, mode)
    ? DERIVED_RESULT_MEMBERS[mode]
    : null;
  if (!result || !isPlainObject(entry?.payload)) return false;
  const payload = entry.payload;

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

  if (!Array.isArray(payload.pairs) || payload.pairs.length === 0) return false;
  const allowedPairKeys = new Set([
    'pair_index',
    'query_index',
    'subject_index',
    'tsv',
    'canonicalTsv',
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
      !hasOneEmptyTableForm(pair) ||
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

  if (Object.hasOwn(payload, result.member) && !isDerivedAnalysisResult(mode, payload[result.member])) return false;
  return Object.entries(payload).every(([key, value]) => (
    ['identity', 'provenance', 'pairs'].includes(key)
    || holdsNoValue(key === result.member && isPlainObject(value) ? value.value : value)
  ));
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
      isZeroHitDerivedResult(entry)
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
  if (!isReusableLosatEntry(entry)) return null;
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

/**
 * The legacy protein artifacts a written Session keeps, from one write's attempts (OV-1132, OV-1143).
 * Python `written_legacy_protein_artifacts` is the other side; both read
 * tests/fixtures/legacy_protein_artifact_vectors.json. A promoted candidate leaves the envelope on
 * serialization; the evidence goes when current derived entries supersede it, or when this write
 * promoted the last pending candidate.
 * @template {{ entries: readonly unknown[] }} E
 * @param {{ candidates: unknown, evidence: E, promoted: readonly number[], hasCurrentDerived: boolean }} write
 * @returns {{ candidates: LegacyProteinCandidateEnvelope,
 *   evidence: E | { schema: number, entries: never[] } }}
 */
export const writtenLegacyProteinArtifacts = ({ candidates, evidence, promoted, hasCurrentDerived }) => {
  const written = promoted.reduce(
    (envelope, index) => transitionLegacyProteinCandidate(envelope, index, 'promoted'),
    normalizeLegacyProteinCandidateEnvelope(candidates)
  );
  const promotedTheLast = promoted.length > 0
    && !written.entries.some((candidate) => candidate.state === 'pending');
  return {
    candidates: written,
    evidence: hasCurrentDerived || promotedTheLast
      ? { schema: LEGACY_LOSAT_DERIVED_CACHE_SCHEMA, entries: [] }
      : evidence
  };
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
