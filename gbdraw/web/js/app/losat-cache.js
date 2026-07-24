export const PROTEIN_LOSAT_CACHE_SCHEMA = 3;
export const NUCLEOTIDE_LOSAT_CACHE_SCHEMA = 2;
export const LOSAT_DERIVED_CACHE_SCHEMA = 2;
export const LEGACY_LOSAT_DERIVED_CACHE_SCHEMA = 1;
export const PROTEIN_IDENTITY_MANIFEST_SCHEMA = 1;
export const LEGACY_PROTEIN_CANDIDATE_SCHEMA = 1;
export const SESSION_LOSAT_CACHE_BYTE_LIMIT = 109_051_904;

const LEGACY_CANDIDATE_STATES = new Set(['pending', 'promoted', 'rejected']);
const LOSAT_OUTFMT6_COLUMN_COUNT = 12;
const EMPTY_CACHE_ENVELOPE_BYTE_LENGTH = 14;

export const isPlainObject = (value) => (
  Boolean(value) && typeof value === 'object' && !Array.isArray(value)
);

const utf8ByteLength = (value) => {
  const text = String(value ?? '');
  let length = 0;
  for (let index = 0; index < text.length; index += 1) {
    const codePoint = text.charCodeAt(index);
    if (codePoint < 0x80) {
      length += 1;
    } else if (codePoint < 0x800) {
      length += 2;
    } else if (
      codePoint >= 0xd800 &&
      codePoint <= 0xdbff &&
      index + 1 < text.length &&
      text.charCodeAt(index + 1) >= 0xdc00 &&
      text.charCodeAt(index + 1) <= 0xdfff
    ) {
      length += 4;
      index += 1;
    } else {
      length += 3;
    }
  }
  return length;
};

const serializedJsonByteLength = (value) => {
  const serialized = JSON.stringify(value);
  return typeof serialized === 'string' ? utf8ByteLength(serialized) : 0;
};

export const serializedLosatArtifactByteLength = ({
  rawEntries = [],
  derivedEntries = []
} = {}) => (
  (2 * EMPTY_CACHE_ENVELOPE_BYTE_LENGTH) +
  rawEntries.reduce((sum, entry) => sum + serializedJsonByteLength(entry), 0) +
  derivedEntries.reduce((sum, entry) => sum + serializedJsonByteLength(entry), 0) +
  Math.max(0, rawEntries.length - 1) +
  Math.max(0, derivedEntries.length - 1)
);

export const pruneSerializedLosatArtifacts = ({
  rawEntries = [],
  derivedEntries = [],
  maxBytes = SESSION_LOSAT_CACHE_BYTE_LIMIT
} = {}) => {
  const raw = Array.isArray(rawEntries) ? rawEntries : [];
  const derived = Array.isArray(derivedEntries) ? derivedEntries : [];
  const limit = Number.isFinite(Number(maxBytes))
    ? Math.max(2 * EMPTY_CACHE_ENVELOPE_BYTE_LENGTH, Math.trunc(Number(maxBytes)))
    : SESSION_LOSAT_CACHE_BYTE_LIMIT;
  const rawSizes = raw.map((entry) => serializedJsonByteLength(entry));
  const derivedSizes = derived.map((entry) => serializedJsonByteLength(entry));
  const selectedRaw = new Set();
  const selectedDerived = new Set();
  let serializedBytes = 2 * EMPTY_CACHE_ENVELOPE_BYTE_LENGTH;

  const addEntry = (index, sizes, selected) => {
    if (selected.has(index)) return true;
    const separatorBytes = selected.size > 0 ? 1 : 0;
    const entryBytes = sizes[index] + separatorBytes;
    if (serializedBytes + entryBytes > limit) return false;
    selected.add(index);
    serializedBytes += entryBytes;
    return true;
  };

  // Reserve the last serialized raw result, then scan derived and remaining
  // raw entries in reverse serialized order without changing output order.
  const reservedRawIndex = raw.length - 1;
  if (reservedRawIndex >= 0) {
    addEntry(reservedRawIndex, rawSizes, selectedRaw);
  }

  for (let index = derived.length - 1; index >= 0; index -= 1) {
    addEntry(index, derivedSizes, selectedDerived);
  }

  for (let index = raw.length - 1; index >= 0; index -= 1) {
    addEntry(index, rawSizes, selectedRaw);
  }

  return {
    rawEntries: raw.filter((_entry, index) => selectedRaw.has(index)),
    derivedEntries: derived.filter((_entry, index) => selectedDerived.has(index)),
    serializedBytes
  };
};

const cloneJson = (value) => JSON.parse(JSON.stringify(value));

export const normalizeLosatArgs = (args) => (
  Array.isArray(args) ? args.map((arg) => String(arg)) : []
);

export const sameLosatArgs = (left, right) => {
  const a = normalizeLosatArgs(left);
  const b = normalizeLosatArgs(right);
  return a.length === b.length && a.every((value, index) => value === b[index]);
};

const hasRawShape = (entry) => (
  isPlainObject(entry) &&
  entry.kind === 'raw-losat' &&
  typeof entry.text === 'string'
);

export const isProteinRawLosatCacheEntry = (entry) => (
  hasRawShape(entry) &&
  entry.schema === PROTEIN_LOSAT_CACHE_SCHEMA &&
  entry.identityKind === 'protein' &&
  String(entry.program || '').toLowerCase() === 'blastp' &&
  typeof entry.queryProteinSetHash === 'string' &&
  Boolean(entry.queryProteinSetHash) &&
  typeof entry.subjectProteinSetHash === 'string' &&
  Boolean(entry.subjectProteinSetHash) &&
  typeof entry.queryBindingHash === 'string' &&
  Boolean(entry.queryBindingHash) &&
  typeof entry.subjectBindingHash === 'string' &&
  Boolean(entry.subjectBindingHash) &&
  typeof entry.queryRecordInstanceKey === 'string' &&
  Boolean(entry.queryRecordInstanceKey) &&
  typeof entry.subjectRecordInstanceKey === 'string' &&
  Boolean(entry.subjectRecordInstanceKey)
);

export const isLegacyProteinRawLosatCacheEntry = (entry) => (
  hasRawShape(entry) &&
  entry.schema === NUCLEOTIDE_LOSAT_CACHE_SCHEMA &&
  String(entry.program || '').toLowerCase() === 'blastp'
);

export const isNucleotideRawLosatCacheEntry = (entry) => {
  if (!hasRawShape(entry) || entry.schema !== NUCLEOTIDE_LOSAT_CACHE_SCHEMA) return false;
  if (entry.identityKind && entry.identityKind !== 'nucleotide') return false;
  return String(entry.program || '').toLowerCase() !== 'blastp';
};

export const classifyRawLosatCacheEntry = (entry) => {
  if (isProteinRawLosatCacheEntry(entry)) return 'protein-current';
  if (isNucleotideRawLosatCacheEntry(entry)) return 'nucleotide-current';
  if (isLegacyProteinRawLosatCacheEntry(entry)) return 'protein-legacy';
  return 'invalid';
};

export const isCurrentRawLosatCacheEntry = (entry) => {
  const classification = classifyRawLosatCacheEntry(entry);
  return classification === 'protein-current' || classification === 'nucleotide-current';
};

export const isLosatDerivedCacheEntry = (entry, { allowLegacy = true } = {}) => (
  isPlainObject(entry) &&
  (
    entry.schema === LOSAT_DERIVED_CACHE_SCHEMA ||
    (allowLegacy && entry.schema === LEGACY_LOSAT_DERIVED_CACHE_SCHEMA)
  ) &&
  entry.kind === 'derived-losatp-payload' &&
  typeof entry.key === 'string' &&
  Boolean(entry.key) &&
  isPlainObject(entry.payload)
);

const stableObjectText = (value) => {
  if (Array.isArray(value)) return `[${value.map(stableObjectText).join(',')}]`;
  if (isPlainObject(value)) {
    return `{${Object.keys(value).sort().map((key) => (
      `${JSON.stringify(key)}:${stableObjectText(value[key])}`
    )).join(',')}}`;
  }
  return JSON.stringify(value);
};

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

export const emptyProteinIdentityManifest = () => ({
  schema: PROTEIN_IDENTITY_MANIFEST_SCHEMA,
  proteinSets: {},
  recordAnalyses: {},
  recordInstances: {}
});

export const validateProteinIdentityManifest = (manifest) => {
  if (!isPlainObject(manifest) || manifest.schema !== PROTEIN_IDENTITY_MANIFEST_SCHEMA) {
    return false;
  }
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
        !featureId.startsWith('f_') ||
        featureIds.has(featureId)
      ) return false;
      featureIds.add(featureId);
    }
    proteinSetFeatureIds.set(proteinSetHash, featureIds);
  }

  const allTransportIds = new Set();
  for (const analysis of Object.values(manifest.recordAnalyses)) {
    if (!isPlainObject(analysis) || analysis.schema !== 1) return false;
    if (
      typeof analysis.proteinSetHash !== 'string' ||
      !proteinSetFeatureIds.has(analysis.proteinSetHash)
    ) {
      return false;
    }
  }
  for (const [instanceKey, instance] of Object.entries(manifest.recordInstances)) {
    if (!instanceKey || !isPlainObject(instance) || instance.schema !== 1) return false;
    if (!Object.prototype.hasOwnProperty.call(manifest.recordAnalyses, instance.recordAnalysisId)) {
      return false;
    }
    if (typeof instance.bindingHash !== 'string' || !instance.bindingHash) return false;
    if (!isPlainObject(instance.transportIds)) return false;
    const analysis = manifest.recordAnalyses[instance.recordAnalysisId];
    const expectedFeatureIds = proteinSetFeatureIds.get(analysis.proteinSetHash);
    if (Object.keys(instance.transportIds).length !== expectedFeatureIds.size) return false;
    for (const [featureId, transportId] of Object.entries(instance.transportIds)) {
      if (!expectedFeatureIds.has(featureId)) return false;
      if (
        typeof transportId !== 'string' ||
        !transportId ||
        /\s/.test(transportId) ||
        !transportId.endsWith(`~${featureId}`)
      ) return false;
      if (allTransportIds.has(transportId)) return false;
      allTransportIds.add(transportId);
    }
  }
  return true;
};

export const mergeProteinIdentityManifests = (manifests) => {
  const merged = emptyProteinIdentityManifest();
  (Array.isArray(manifests) ? manifests : []).forEach((manifest) => {
    if (!validateProteinIdentityManifest(manifest)) {
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

export const validateProteinRawEntryReferences = (entry, manifest) => {
  if (!isProteinRawLosatCacheEntry(entry) || !validateProteinIdentityManifest(manifest)) {
    return false;
  }
  const queryInstance = manifest.recordInstances[entry.queryRecordInstanceKey];
  const subjectInstance = manifest.recordInstances[entry.subjectRecordInstanceKey];
  if (!queryInstance || !subjectInstance) return false;
  if (queryInstance.bindingHash !== entry.queryBindingHash) return false;
  if (subjectInstance.bindingHash !== entry.subjectBindingHash) return false;
  const queryAnalysis = manifest.recordAnalyses[queryInstance.recordAnalysisId];
  const subjectAnalysis = manifest.recordAnalyses[subjectInstance.recordAnalysisId];
  return (
    queryAnalysis?.proteinSetHash === entry.queryProteinSetHash &&
    subjectAnalysis?.proteinSetHash === entry.subjectProteinSetHash
  );
};

export const proteinTransportIdSets = (manifest, queryInstanceKey, subjectInstanceKey) => {
  if (!validateProteinIdentityManifest(manifest)) return null;
  const query = manifest.recordInstances[queryInstanceKey];
  const subject = manifest.recordInstances[subjectInstanceKey];
  if (!query || !subject) return null;
  return {
    query: new Set(Object.values(query.transportIds)),
    subject: new Set(Object.values(subject.transportIds))
  };
};

export const rawProteinTextMatchesBindings = (text, queryIds, subjectIds) => {
  if (typeof text !== 'string' || !(queryIds instanceof Set) || !(subjectIds instanceof Set)) {
    return false;
  }
  for (const rawLine of text.split(/\r?\n/)) {
    const line = rawLine.trim();
    if (!line || line.startsWith('#')) continue;
    const columns = rawLine.split('\t');
    if (
      columns.length !== LOSAT_OUTFMT6_COLUMN_COUNT ||
      !queryIds.has(columns[0]) ||
      !subjectIds.has(columns[1])
    ) {
      return false;
    }
  }
  return true;
};

export const getCurrentRawLosatCacheEntry = (cacheMap, cacheKey, metadata = {}, manifest = null) => {
  if (!(cacheMap instanceof Map) || !cacheKey) return null;
  const entry = cacheMap.get(cacheKey);
  const classification = classifyRawLosatCacheEntry(entry);
  if (classification === 'protein-current') {
    if (String(entry.program || '') !== String(metadata.program || 'blastp')) return null;
    if (String(entry.outfmt || '6') !== String(metadata.outfmt || '6')) return null;
    if (!sameLosatArgs(entry.args, metadata.args)) return null;
    if (metadata.queryBindingHash && entry.queryBindingHash !== metadata.queryBindingHash) return null;
    if (metadata.subjectBindingHash && entry.subjectBindingHash !== metadata.subjectBindingHash) return null;
    if (!validateProteinRawEntryReferences(entry, manifest)) return null;
    const ids = proteinTransportIdSets(
      manifest,
      entry.queryRecordInstanceKey,
      entry.subjectRecordInstanceKey
    );
    if (!ids || !rawProteinTextMatchesBindings(entry.text, ids.query, ids.subject)) return null;
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

export const serializableLegacyProteinCandidateEnvelope = (envelope) => {
  const normalized = normalizeLegacyProteinCandidateEnvelope(envelope);
  return {
    schema: LEGACY_PROTEIN_CANDIDATE_SCHEMA,
    entries: normalized.entries.filter((candidate) => candidate.state !== 'promoted')
  };
};

export const candidateOriginalEntries = (envelope, { states = ['pending'] } = {}) => {
  const acceptedStates = new Set(states);
  return normalizeLegacyProteinCandidateEnvelope(envelope).entries
    .map((candidate, index) => ({ candidate, index }))
    .filter(({ candidate }) => acceptedStates.has(candidate.state))
    .map(({ candidate, index }) => ({ index, entry: cloneJson(candidate.originalEntry) }));
};
