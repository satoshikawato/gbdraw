import { formatFastaEntry } from './feature-sequence-fasta.js';
import {
  orderedConservationSources,
  orderedOptionalConservationFiles
} from './conservation-series.js';
import { readFileText } from '../services/file-content-cache.js';

const IUPAC_COMPLEMENT = Object.freeze({
  A: 'T', C: 'G', G: 'C', T: 'A', U: 'A',
  R: 'Y', Y: 'R', S: 'S', W: 'W', K: 'M', M: 'K',
  B: 'V', D: 'H', H: 'D', V: 'B', N: 'N',
  '-': '-'
});

const text = (value) => String(value ?? '').trim();
const normalizedSequence = (value) => String(value ?? '').replace(/\s+/g, '').toUpperCase();
const optionalNonnegativeIntegerStatus = (value) => {
  if (value === null || value === undefined || value === '') {
    return { valid: true, supplied: false, value: null };
  }
  const numeric = typeof value === 'number'
    ? value
    : (typeof value === 'string' && /^\d+$/.test(value.trim())
      ? Number(value.trim())
      : Number.NaN);
  const valid = Number.isSafeInteger(numeric) && numeric >= 0;
  return { valid, supplied: true, value: valid ? numeric : null };
};
const safePart = (value, fallback = 'sequence') => {
  const cleaned = text(value).replace(/[^A-Za-z0-9_.-]+/g, '_').replace(/^_+|_+$/g, '');
  return cleaned || fallback;
};

const aliasesForSource = (source) => {
  const values = [source?.recordId, ...(Array.isArray(source?.aliases) ? source.aliases : [])];
  const aliases = new Set();
  values.forEach((value) => {
    const full = text(value);
    if (!full) return;
    aliases.add(full);
    const first = full.split(/\s+/)[0];
    if (first) aliases.add(first);
  });
  return aliases;
};

const sequenceSourcesAgree = (left, right) => (
  left.key === right.key
  && left.recordId === right.recordId
  && left.sequence === right.sequence
  && left.origin === right.origin
  && left.recordIndex === right.recordIndex
  && left.sourceIndex === right.sourceIndex
  && (
    left.aliases.length === right.aliases.length
    && left.aliases.every((alias) => right.aliases.includes(alias))
  )
);

export const reverseComplementNucleotide = (sequence) => {
  const input = normalizedSequence(sequence);
  let output = '';
  for (let index = input.length - 1; index >= 0; index -= 1) {
    output += IUPAC_COMPLEMENT[input[index]] || 'N';
  }
  return output;
};

export const validateMatchCoordinates = (startRaw, endRaw, sequenceLength) => {
  const start = Number(startRaw);
  const end = Number(endRaw);
  const length = Number(sequenceLength);
  if (!Number.isInteger(start) || !Number.isInteger(end)) {
    return { valid: false, reason: 'Coordinates must be whole 1-based values.' };
  }
  if (start < 1 || end < 1) {
    return { valid: false, reason: 'Coordinates must be greater than or equal to 1.' };
  }
  if (!Number.isInteger(length) || length < 1) {
    return { valid: false, reason: 'The sequence source is empty.' };
  }
  if (start > length || end > length) {
    return { valid: false, reason: `Coordinates exceed the sequence length (${length.toLocaleString()} bp).` };
  }
  return { valid: true, start, end, orientation: start <= end ? '+' : '-' };
};

export const extractMatchedSpan = (sequence, startRaw, endRaw) => {
  const normalized = normalizedSequence(sequence);
  const validation = validateMatchCoordinates(startRaw, endRaw, normalized.length);
  if (!validation.valid) return { ...validation, sequence: '' };
  const low = Math.min(validation.start, validation.end);
  const high = Math.max(validation.start, validation.end);
  const sliced = normalized.slice(low - 1, high);
  return {
    ...validation,
    sequence: validation.orientation === '-' ? reverseComplementNucleotide(sliced) : sliced,
    sequenceLength: high - low + 1
  };
};

export const createSequenceSourceRegistry = (initialSources = []) => {
  const sources = new Map();
  const ambiguousKeys = new Set();

  const register = (source) => {
    const key = text(source?.key);
    const sequence = normalizedSequence(source?.sequence);
    const recordIndex = optionalNonnegativeIntegerStatus(source?.recordIndex);
    const sourceIndex = optionalNonnegativeIntegerStatus(source?.sourceIndex);
    if (!key || !sequence || !recordIndex.valid || !sourceIndex.valid) return null;
    const entry = {
      key,
      recordId: text(source?.recordId) || key,
      aliases: Array.from(aliasesForSource(source)),
      sequence,
      origin: text(source?.origin),
      recordIndex: recordIndex.value,
      sourceIndex: sourceIndex.value
    };
    if (ambiguousKeys.has(key)) return null;
    const existing = sources.get(key);
    if (existing) {
      if (sequenceSourcesAgree(existing, entry)) return existing;
      sources.delete(key);
      ambiguousKeys.add(key);
      return null;
    }
    sources.set(key, entry);
    return entry;
  };

  const reset = (nextSources = []) => {
    sources.clear();
    ambiguousKeys.clear();
    (Array.isArray(nextSources) ? nextSources : []).forEach(register);
  };

  // Artifact History only supplies entries previously produced by this registry.
  // Re-adopt those immutable entries without normalizing or copying large sequence
  // strings a second time.
  const resetTrusted = (nextSources = []) => {
    sources.clear();
    ambiguousKeys.clear();
    (Array.isArray(nextSources) ? nextSources : []).forEach((source) => {
      const key = text(source?.key);
      if (key) sources.set(key, source);
    });
  };

  const resolve = (sourceKey, recordId, context = {}) => {
    const expectedSourceIndex = optionalNonnegativeIntegerStatus(context.sourceIndex);
    const expectedRecordIndex = optionalNonnegativeIntegerStatus(context.recordIndex);
    if (!expectedSourceIndex.valid || !expectedRecordIndex.valid) {
      return { source: null, reason: 'The sequence source identity for this match is invalid.' };
    }
    const expectedOrigin = text(context.origin);
    const wanted = text(recordId);
    const agrees = (source) => (
      (!expectedOrigin || source.origin === expectedOrigin)
      && (!expectedSourceIndex.supplied || source.sourceIndex === expectedSourceIndex.value)
      && (!expectedRecordIndex.supplied || source.recordIndex === expectedRecordIndex.value)
      && (!wanted || source.recordId === wanted || aliasesForSource(source).has(wanted))
    );
    const explicitKey = text(sourceKey);
    if (explicitKey) {
      if (ambiguousKeys.has(explicitKey)) {
        return { source: null, reason: 'The sequence source key for this match is ambiguous.' };
      }
      const direct = sources.get(explicitKey);
      if (!direct) {
        return { source: null, reason: 'The sequence source for this match is unavailable.' };
      }
      return agrees(direct)
        ? { source: direct, reason: '' }
        : { source: null, reason: 'The sequence source key conflicts with the match identity.' };
    }

    let candidates = Array.from(sources.values()).filter((source) => {
      if (expectedOrigin && source.origin !== expectedOrigin) return false;
      if (expectedSourceIndex.supplied && source.sourceIndex !== expectedSourceIndex.value) return false;
      if (expectedRecordIndex.supplied && source.recordIndex !== expectedRecordIndex.value) return false;
      return true;
    });

    if (expectedRecordIndex.supplied && candidates.length === 1 && !wanted) {
      return { source: candidates[0], reason: '' };
    }
    if (!wanted) return { source: null, reason: 'The match does not identify a sequence record.' };
    const exact = candidates.filter((source) => source.recordId === wanted);
    if (exact.length === 1) return { source: exact[0], reason: '' };
    if (exact.length > 1) {
      return { source: null, reason: `Record ID “${wanted}” is ambiguous in the available sequence sources.` };
    }
    const alias = candidates.filter((source) => aliasesForSource(source).has(wanted));
    if (alias.length === 1) return { source: alias[0], reason: '' };
    if (alias.length > 1) {
      return { source: null, reason: `Record alias “${wanted}” is ambiguous in the available sequence sources.` };
    }
    return { source: null, reason: `No sequence source matched record “${wanted}”.` };
  };

  (Array.isArray(initialSources) ? initialSources : []).forEach(register);
  return {
    sources,
    register,
    reset,
    resetTrusted,
    resolve,
    values: () => Array.from(sources.values())
  };
};

const sourceIdentity = (source, overrides = {}) => {
  const recordIndex = optionalNonnegativeIntegerStatus(
    overrides.recordIndex ?? source?.recordIndex
  );
  const sourceIndex = optionalNonnegativeIntegerStatus(
    overrides.sourceIndex ?? source?.sourceIndex
  );
  return {
    key: text(overrides.key ?? source?.key),
    recordId: text(overrides.recordId ?? source?.recordId),
    origin: text(overrides.origin ?? source?.origin),
    recordIndex: recordIndex.valid && recordIndex.supplied ? recordIndex.value : null,
    sourceIndex: sourceIndex.valid && sourceIndex.supplied ? sourceIndex.value : null
  };
};

const sourceIdentityKey = (identity) => JSON.stringify([
  identity.origin,
  identity.recordIndex,
  identity.sourceIndex,
  identity.recordId,
  identity.key
]);

const consistentTextAlias = (value, keys, fallback = '') => {
  const values = new Set(
    keys
      .filter((key) => Object.prototype.hasOwnProperty.call(value || {}, key))
      .map((key) => text(value?.[key]))
      .filter(Boolean)
  );
  return {
    valid: values.size <= 1,
    value: values.size ? Array.from(values)[0] : fallback
  };
};

const consistentIntegerAlias = (value, keys) => {
  const statuses = keys
    .filter((key) => Object.prototype.hasOwnProperty.call(value || {}, key))
    .map((key) => optionalNonnegativeIntegerStatus(value?.[key]));
  const supplied = statuses.filter((status) => status.supplied);
  const values = new Set(supplied.filter((status) => status.valid).map((status) => status.value));
  return {
    valid: statuses.every((status) => status.valid) && values.size <= 1,
    supplied: supplied.length > 0,
    value: values.size === 1 ? Array.from(values)[0] : null
  };
};

const sourceResolutionContext = (identity) => ({
  origin: identity.origin,
  ...(identity.recordIndex === null ? {} : { recordIndex: identity.recordIndex }),
  ...(identity.sourceIndex === null ? {} : { sourceIndex: identity.sourceIndex })
});

const invalidCatalogSourceReason = (source, identity) => {
  if (!source || typeof source !== 'object' || Array.isArray(source)) {
    return 'The catalog sequence source is not an object.';
  }
  if (!identity.key) return 'The catalog sequence source key is empty.';
  if (!identity.recordId) return 'The catalog sequence source record ID is empty.';
  if (
    typeof source.sequence !== 'string'
    || !source.sequence
    || /\s/.test(source.sequence)
  ) {
    return 'The catalog sequence source has no canonical nucleotide sequence.';
  }
  if (!['linear-record', 'circular-reference', 'homology-comparison'].includes(identity.origin)) {
    return 'The catalog sequence source origin is invalid.';
  }
  const recordIndex = optionalNonnegativeIntegerStatus(source.recordIndex);
  const sourceIndex = optionalNonnegativeIntegerStatus(source.sourceIndex);
  if (!recordIndex.valid || !sourceIndex.valid) {
    return 'The catalog sequence source index identity is invalid.';
  }
  if (
    ['linear-record', 'circular-reference'].includes(identity.origin)
    && (!recordIndex.supplied || identity.recordIndex === null)
  ) {
    return 'The catalog primary sequence source has no record index.';
  }
  if (
    identity.origin === 'homology-comparison'
    && (!sourceIndex.supplied || identity.sourceIndex === null)
  ) {
    return 'The catalog comparison sequence source has no source index.';
  }
  return '';
};

export const resolveCircularComparisonSequenceAvailability = ({
  files,
  circularConservation
} = {}) => {
  const sourceMode = text(circularConservation?.source).toLowerCase();
  if (!['upload', 'losat'].includes(sourceMode)) return [];

  const comparisonFiles = sourceMode === 'upload'
    ? (files?.c_conservation_sequence_sources || [])
    : (files?.c_conservation_fastas || []);
  const orderedComparisonFiles = sourceMode === 'upload'
    ? orderedConservationSources(
        files?.c_conservation_blasts || [],
        circularConservation
      ).map((entry) => comparisonFiles[entry.sourceIndex] || null)
    : orderedOptionalConservationFiles(
        comparisonFiles,
        circularConservation
      );
  return orderedComparisonFiles.map((file, sourceIndex) => ({
    sourceIndex,
    file: file || null,
    availability: file ? 'recoverable' : 'never-supplied'
  }));
};

const comparisonSourceAvailabilityAt = (inventory, sourceIndex) => {
  if (!Array.isArray(inventory)) {
    return {
      valid: false,
      reason: 'The Circular comparison source availability inventory is missing.'
    };
  }
  const matches = inventory.filter((entry) => (
    entry
    && typeof entry === 'object'
    && !Array.isArray(entry)
    && entry.sourceIndex === sourceIndex
  ));
  if (matches.length === 0) {
    return {
      valid: false,
      reason: `The Circular comparison source availability is missing for source index ${sourceIndex}.`
    };
  }
  if (matches.length !== 1) {
    return {
      valid: false,
      reason: `The Circular comparison source availability conflicts for source index ${sourceIndex}.`
    };
  }

  const entry = matches[0];
  const availability = text(entry.availability).toLowerCase();
  if (!['recoverable', 'never-supplied'].includes(availability)) {
    return {
      valid: false,
      reason: `The Circular comparison source availability is unknown for source index ${sourceIndex}.`
    };
  }
  const hasFile = Boolean(entry.file);
  if (
    (availability === 'recoverable' && !hasFile)
    || (availability === 'never-supplied' && hasFile)
  ) {
    return {
      valid: false,
      reason: `The Circular comparison source availability conflicts with its FASTA binding at source index ${sourceIndex}.`
    };
  }
  return { valid: true, availability };
};

/**
 * Decide whether a current feature catalog already contains every sequence
 * source that its feature and match-popup consumers require.
 */
export const analyzeCatalogSequenceSourceCoverage = ({
  mode,
  catalogFeatureState,
  renderRequest,
  comparisonSourceAvailability
} = {}) => {
  const normalizedMode = text(mode).toLowerCase();
  const records = Array.isArray(renderRequest?.records) ? renderRequest.records : [];
  const displayedRecords = records.map((record, recordIndex) => ({
    recordIndex,
    recordKey: text(record?.recordKey)
  }));
  const displayedIndexByKey = new Map(
    displayedRecords.map((record) => [record.recordKey, record.recordIndex])
  );
  const items = Array.isArray(catalogFeatureState?.items)
    ? catalogFeatureState.items
    : null;
  const invalidCatalogSources = [];
  if (!['linear', 'circular'].includes(normalizedMode)) {
    invalidCatalogSources.push({ reason: 'The diagram mode is invalid.' });
  }
  if (!items) {
    invalidCatalogSources.push({ reason: 'The current feature catalog has no item list.' });
  }

  const catalogItems = items || [];
  const sourceEntries = [];
  catalogItems.forEach((item, itemIndex) => {
    const sources = Array.isArray(item?.sequenceSources) ? item.sequenceSources : [];
    sources.forEach((source, sequenceSourceIndex) => {
      sourceEntries.push({ itemIndex, sequenceSourceIndex, source });
    });
  });
  const registry = createSequenceSourceRegistry(
    sourceEntries.map((entry) => entry.source)
  );
  sourceEntries.forEach(({ itemIndex, sequenceSourceIndex, source }) => {
    const identity = sourceIdentity(source);
    let reason = invalidCatalogSourceReason(source, identity);
    if (!reason) {
      const resolution = registry.resolve(
        identity.key,
        identity.recordId,
        sourceResolutionContext(identity)
      );
      if (!resolution.source) reason = resolution.reason;
    }
    if (reason) {
      invalidCatalogSources.push({
        itemIndex,
        sequenceSourceIndex,
        source: identity,
        reason
      });
    }
  });

  const requirements = new Map();
  const failedRequirementKeys = new Set();
  const failureReasons = new Map();
  const optionalUnavailableRequirementKeys = new Set();
  const optionalUnavailableReasons = new Map();
  const usedDisplayedRecordKeys = new Set();
  const consumerCounts = { biologicalFeatures: 0, matchEndpoints: 0 };

  const addConsumer = ({ kind, id, recordKey = '', expectedSource, resolution }) => {
    const key = sourceIdentityKey(expectedSource);
    if (!requirements.has(key)) {
      requirements.set(key, {
        expectedSource,
        consumerCount: 0,
        biologicalFeatureCount: 0,
        matchEndpointCount: 0,
        exampleConsumers: []
      });
    }
    const requirement = requirements.get(key);
    requirement.consumerCount += 1;
    if (kind === 'biological-feature') {
      requirement.biologicalFeatureCount += 1;
      consumerCounts.biologicalFeatures += 1;
    } else {
      requirement.matchEndpointCount += 1;
      consumerCounts.matchEndpoints += 1;
    }
    if (requirement.exampleConsumers.length < 3) {
      requirement.exampleConsumers.push({ kind, id });
    }
    if (recordKey) usedDisplayedRecordKeys.add(recordKey);
    if (resolution?.source) return;
    if (resolution?.optionalUnavailable) {
      optionalUnavailableRequirementKeys.add(key);
      if (!optionalUnavailableReasons.has(key)) {
        optionalUnavailableReasons.set(key, new Set());
      }
      optionalUnavailableReasons.get(key).add(
        text(resolution.reason) || 'The optional comparison FASTA was never supplied.'
      );
      return;
    }
    failedRequirementKeys.add(key);
    if (!failureReasons.has(key)) failureReasons.set(key, new Set());
    failureReasons.get(key).add(
      text(resolution?.reason) || 'The required sequence source is unavailable.'
    );
  };

  catalogItems.forEach((item, itemIndex) => {
    const recordKeys = Array.isArray(item?.recordKeys)
      ? item.recordKeys.map(text)
      : [];
    const itemSources = Array.isArray(item?.sequenceSources) ? item.sequenceSources : [];
    const expectedPrimaryOrigin = normalizedMode === 'linear'
      ? 'linear-record'
      : 'circular-reference';
    (Array.isArray(item?.biologicalFeatures) ? item.biologicalFeatures : []).forEach(
      (feature) => {
        if (!Object.prototype.hasOwnProperty.call(feature || {}, 'sequenceSourceIndex')) return;
        const sequenceSourceIndex = feature.sequenceSourceIndex;
        const recordKey = text(feature?.recordKey);
        const itemRecordIndex = recordKeys.indexOf(recordKey);
        const source = Number.isInteger(sequenceSourceIndex)
          ? itemSources[sequenceSourceIndex]
          : null;
        const expectedSource = sourceIdentity(source, {
          key: source?.key || (
            itemRecordIndex >= 0
              ? `${normalizedMode === 'linear' ? 'linear' : 'circular'}:record:${itemRecordIndex}`
              : ''
          ),
          recordId: source?.recordId || feature?.record_id || feature?.recordId,
          origin: expectedPrimaryOrigin,
          recordIndex: itemRecordIndex
        });
        let resolution = { source: null, reason: '' };
        if (
          !Number.isInteger(sequenceSourceIndex)
          || sequenceSourceIndex < 0
          || !source
        ) {
          resolution.reason = 'The biological feature references a missing catalog source index.';
        } else if (itemRecordIndex < 0 || !displayedIndexByKey.has(recordKey)) {
          resolution.reason = 'The biological feature record is not displayed by the render request.';
        } else if (
          text(source.origin) !== expectedPrimaryOrigin
          || source.recordIndex !== itemRecordIndex
        ) {
          resolution.reason = 'The biological feature source conflicts with its record identity.';
        } else {
          resolution = registry.resolve(
            expectedSource.key,
            expectedSource.recordId,
            sourceResolutionContext(expectedSource)
          );
        }
        addConsumer({
          kind: 'biological-feature',
          id: `${recordKey}/${text(feature?.biologicalFeatureId)}`,
          recordKey,
          expectedSource,
          resolution
        });
      }
    );

    (Array.isArray(item?.comparisonMatches) ? item.comparisonMatches : []).forEach(
      (match, matchIndex) => {
        const matchId = consistentTextAlias(match, ['id', 'matchId', 'match_id'], `match-${matchIndex + 1}`);
        const matchKind = consistentTextAlias(match, ['match_kind', 'matchKind'], 'pairwise');
        const normalizedMatchKind = matchKind.value.toLowerCase();
        if (matchKind.valid && normalizedMatchKind === 'orthogroup') return;
        const referenceSide = consistentTextAlias(match, ['reference_side', 'referenceSide']);
        const normalizedReferenceSide = referenceSide.value.toLowerCase();
        const sourceIndex = consistentIntegerAlias(match, ['source_index', 'sourceIndex']);
        ['query', 'subject'].forEach((role) => {
          const recordId = consistentTextAlias(
            match,
            [`${role}_record_id`, `${role}RecordId`]
          );
          const recordIndex = consistentIntegerAlias(
            match,
            [`${role}_record_index`, `${role}RecordIndex`]
          );
          const homology = normalizedMatchKind === 'homology';
          const referenceIdentityValid = ['query', 'subject'].includes(normalizedReferenceSide);
          const expectedOrigin = homology
            ? (role === normalizedReferenceSide ? 'circular-reference' : 'homology-comparison')
            : 'linear-record';
          const displayedRecord = recordIndex.value === null
            ? null
            : displayedRecords[recordIndex.value] || null;
          const expectedSource = sourceIdentity(null, {
            key: expectedOrigin === 'linear-record' && recordIndex.value !== null
              ? `linear:record:${recordIndex.value}`
              : '',
            recordId: recordId.value,
            origin: expectedOrigin,
            recordIndex: expectedOrigin === 'linear-record' ? recordIndex.value : null,
            sourceIndex: expectedOrigin === 'homology-comparison' ? sourceIndex.value : null
          });
          let resolution = { source: null, reason: '' };
          if (!matchId.valid || !matchKind.valid || !recordId.valid) {
            resolution.reason = 'The match endpoint has conflicting sequence identity aliases.';
          } else if (homology && !referenceIdentityValid) {
            resolution.reason = 'The homology match reference side is invalid.';
          } else if (
            expectedOrigin === 'linear-record'
            && (!recordIndex.valid || !recordIndex.supplied || !displayedRecord)
          ) {
            resolution.reason = 'The match endpoint record index is invalid or not displayed.';
          } else if (
            expectedOrigin === 'homology-comparison'
            && (!sourceIndex.valid || !sourceIndex.supplied)
          ) {
            resolution.reason = 'The homology comparison source index is invalid.';
          } else {
            resolution = registry.resolve(
              expectedSource.key,
              expectedSource.recordId,
              sourceResolutionContext(expectedSource)
            );
            if (
              !resolution.source
              && expectedOrigin === 'homology-comparison'
              && expectedSource.recordId
            ) {
              const availability = comparisonSourceAvailabilityAt(
                comparisonSourceAvailability,
                expectedSource.sourceIndex
              );
              if (!availability.valid) {
                resolution = { source: null, reason: availability.reason };
              } else if (availability.availability === 'never-supplied') {
                resolution = {
                  source: null,
                  optionalUnavailable: true,
                  reason: 'The optional comparison FASTA was never supplied.'
                };
              }
            }
          }
          addConsumer({
            kind: 'match-endpoint',
            id: `${matchId.value}/${role}`,
            recordKey: displayedRecord?.recordKey || '',
            expectedSource,
            resolution
          });
        });
      }
    );
  });

  const requiredConsumers = Array.from(requirements.entries()).map(([key, value]) => ({
    key,
    ...value
  }));
  const missingConsumers = requiredConsumers
    .filter(({ key }) => failedRequirementKeys.has(key))
    .map((requirement) => ({
      ...requirement,
      reasons: Array.from(failureReasons.get(requirement.key) || [])
    }));
  const optionalUnavailableConsumers = requiredConsumers
    .filter(({ key }) => (
      optionalUnavailableRequirementKeys.has(key)
      && !failedRequirementKeys.has(key)
    ))
    .map((requirement) => ({
      ...requirement,
      reasons: Array.from(optionalUnavailableReasons.get(requirement.key) || [])
    }));
  const resolvedConsumers = requiredConsumers.filter(
    ({ key }) => (
      !failedRequirementKeys.has(key)
      && !optionalUnavailableRequirementKeys.has(key)
    )
  );
  const displayedRecordsWithoutConsumers = displayedRecords.filter(
    ({ recordKey }) => !usedDisplayedRecordKeys.has(recordKey)
  );
  return {
    complete: invalidCatalogSources.length === 0 && missingConsumers.length === 0,
    requiredConsumers,
    resolvedConsumers,
    missingConsumers,
    optionalUnavailableConsumers,
    invalidCatalogSources,
    consumerCounts,
    displayedRecordsWithoutConsumers
  };
};

const normalizeResolution = (resolved) => {
  if (!resolved) return { source: null, reason: '' };
  if (resolved.source !== undefined) return resolved;
  return { source: resolved, reason: '' };
};

export const buildMatchSequenceEntry = (
  span,
  {
    matchId = 'match',
    resolveSequenceSource,
    context = {},
    unavailableReason = ''
  } = {}
) => {
  const normalizedSpan = {
    role: text(span?.role) === 'subject' ? 'subject' : 'query',
    sourceKey: text(span?.sourceKey),
    recordId: text(span?.recordId),
    start: span?.start,
    end: span?.end,
    molecule: 'nucleotide',
    displayRole: text(span?.displayRole) || (text(span?.role) === 'subject' ? 'Subject' : 'Query')
  };
  const unavailable = (reason) => ({
    span: normalizedSpan,
    orientation: Number(normalizedSpan.start) <= Number(normalizedSpan.end) ? '+' : '-',
    sequenceLength: 0,
    fasta: '',
    filename: '',
    available: false,
    unavailableReason: reason || unavailableReason || 'Sequence is unavailable.'
  });
  if (unavailableReason) return unavailable(unavailableReason);
  if (typeof resolveSequenceSource !== 'function') return unavailable('Sequence source resolver is unavailable.');
  const resolution = normalizeResolution(resolveSequenceSource(
    normalizedSpan.sourceKey,
    normalizedSpan.recordId,
    context
  ));
  if (!resolution.source) return unavailable(resolution.reason);
  const extracted = extractMatchedSpan(
    resolution.source.sequence,
    normalizedSpan.start,
    normalizedSpan.end
  );
  if (!extracted.valid) return unavailable(extracted.reason);
  const recordId = normalizedSpan.recordId || resolution.source.recordId;
  const header = [
    `${safePart(matchId, 'match')}_${normalizedSpan.role}`,
    `record=${recordId}`,
    `coords=${normalizedSpan.start}..${normalizedSpan.end}`,
    `strand=${extracted.orientation}`
  ].join('|');
  const fasta = `${formatFastaEntry({ id: header, sequence: extracted.sequence })}\n`;
  return {
    span: { ...normalizedSpan, recordId },
    orientation: extracted.orientation,
    sequenceLength: extracted.sequenceLength,
    fasta,
    filename: `${safePart(matchId)}_${normalizedSpan.role}_${safePart(recordId)}_${normalizedSpan.start}-${normalizedSpan.end}.fna`,
    available: true,
    unavailableReason: ''
  };
};

export const buildMatchSequenceBundle = (
  spans,
  {
    matchId = 'match',
    resolveSequenceSource,
    contextForSpan = () => ({}),
    unavailableReasonForSpan = () => ''
  } = {}
) => {
  const entries = (Array.isArray(spans) ? spans : []).map((span) => buildMatchSequenceEntry(span, {
    matchId,
    resolveSequenceSource,
    context: contextForSpan(span),
    unavailableReason: unavailableReasonForSpan(span)
  }));
  const availableEntries = entries.filter((entry) => entry.available);
  return {
    entries,
    combinedFasta: availableEntries.length === entries.length && entries.length > 1
      ? availableEntries.map((entry) => entry.fasta).join('')
      : '',
    combinedFilename: `${safePart(matchId)}_both.fna`
  };
};

export const parseSequenceRecords = (value) => {
  const input = String(value ?? '');
  const records = [];
  let current = null;
  input.split(/\r?\n/).forEach((line) => {
    if (line.startsWith('>')) {
      if (current) records.push({ ...current, sequence: current.parts.join('').toUpperCase() });
      const header = line.slice(1).trim();
      current = { recordId: header.split(/\s+/)[0] || `record_${records.length + 1}`, header, parts: [] };
    } else if (current && line.trim()) {
      current.parts.push(line.replace(/\s+/g, ''));
    }
  });
  if (current) records.push({ ...current, sequence: current.parts.join('').toUpperCase() });
  return records.filter((record) => record.sequence);
};

const parseGenbankSequenceRecords = (value) => {
  const records = [];
  String(value ?? '').split(/^\/\/\s*$/m).forEach((chunk) => {
    const originMatch = chunk.match(/\nORIGIN\b([\s\S]*)$/i);
    if (!originMatch) return;
    const locus = text(chunk.match(/^LOCUS\s+(\S+)/m)?.[1]);
    const accession = text(chunk.match(/^ACCESSION\s+(\S+)/m)?.[1]);
    const version = text(chunk.match(/^VERSION\s+(\S+)/m)?.[1]);
    const sequence = normalizedSequence(originMatch[1].replace(/[^A-Za-z]/g, ''));
    if (!sequence) return;
    const recordId = version || accession || locus || `record_${records.length + 1}`;
    records.push({
      recordId,
      header: recordId,
      aliases: Array.from(new Set([recordId, version, accession, locus].filter(Boolean))),
      sequence
    });
  });
  return records;
};

const parseInputSequenceRecords = (fileText, inputType) => (
  text(inputType).toLowerCase() === 'gb'
    ? parseGenbankSequenceRecords(fileText)
    : parseSequenceRecords(fileText).map((record) => ({
        ...record,
        aliases: Array.from(new Set([record.recordId, record.header].filter(Boolean)))
      }))
);

const selectInputSequenceRecord = (records, selectorValue) => {
  if (!records.length) return null;
  const selector = text(selectorValue);
  if (!selector) return records[0];
  if (selector.startsWith('#')) {
    const recordIndex = Number(selector.slice(1)) - 1;
    return Number.isInteger(recordIndex) ? records[recordIndex] || null : null;
  }
  const matches = records.filter((record) => (
    record.recordId === selector ||
    (Array.isArray(record.aliases) && record.aliases.includes(selector))
  ));
  return matches.length === 1 ? matches[0] : null;
};

const materializeLinearSequenceRecord = (record, sequenceState) => {
  if (!record) return null;
  const startRaw = sequenceState?.region_start;
  const endRaw = sequenceState?.region_end;
  const hasStart = startRaw !== null && startRaw !== undefined && startRaw !== '';
  const hasEnd = endRaw !== null && endRaw !== undefined && endRaw !== '';
  if (hasStart !== hasEnd) return null;

  let sequence = normalizedSequence(record.sequence);
  let reverse = Boolean(sequenceState?.region_reverse);
  if (hasStart) {
    const start = Number(startRaw);
    const end = Number(endRaw);
    if (
      !Number.isInteger(start) ||
      !Number.isInteger(end) ||
      start < 1 ||
      end < 1 ||
      start > sequence.length ||
      end > sequence.length
    ) return null;
    sequence = sequence.slice(Math.min(start, end) - 1, Math.max(start, end));
    reverse = reverse || start > end;
  }
  return {
    ...record,
    sequence: reverse ? reverseComplementNucleotide(sequence) : sequence
  };
};

const readInputSequenceRecords = async (file, inputType) => {
  if (!file) return [];
  return parseInputSequenceRecords(await readFileText(file), inputType);
};

/**
 * Rebuild the transient match-sequence registry after loading a web session.
 *
 * Sessions already contain the source files, so sequence strings do not need a
 * second persisted copy. Sources that were never supplied (for example, a
 * BLAST upload without its optional comparison FASTA) remain unavailable.
 */
export const buildRestoredMatchSequenceSources = async ({
  mode,
  cInputType,
  lInputType,
  files,
  linearSeqs,
  circularConservation
} = {}) => {
  const sources = [];
  if (text(mode).toLowerCase() === 'linear') {
    const inputType = text(lInputType).toLowerCase() === 'gff' ? 'fasta' : 'gb';
    for (let recordIndex = 0; recordIndex < (linearSeqs || []).length; recordIndex += 1) {
      const sequenceState = linearSeqs[recordIndex] || {};
      const file = inputType === 'gb' ? sequenceState.gb : sequenceState.fasta;
      const records = await readInputSequenceRecords(file, inputType);
      const selected = selectInputSequenceRecord(records, sequenceState.region_record_id);
      const record = materializeLinearSequenceRecord(selected, sequenceState);
      if (!record?.sequence) continue;
      sources.push({
        key: `linear:record:${recordIndex}`,
        recordId: record.recordId,
        aliases: record.aliases,
        sequence: record.sequence,
        origin: 'linear-record',
        recordIndex
      });
    }
    return sources;
  }

  const circularInputType = text(cInputType).toLowerCase() === 'gff' ? 'fasta' : 'gb';
  const referenceFile = circularInputType === 'gb' ? files?.c_gb : files?.c_fasta;
  const referenceRecords = await readInputSequenceRecords(referenceFile, circularInputType);
  referenceRecords.forEach((record, recordIndex) => {
    sources.push({
      key: `circular:record:${recordIndex}`,
      recordId: record.recordId,
      aliases: record.aliases,
      sequence: record.sequence,
      origin: 'circular-reference',
      recordIndex
    });
  });

  const comparisonSourceAvailability = resolveCircularComparisonSequenceAvailability({
    files,
    circularConservation
  });
  for (const { sourceIndex: displayIndex, file } of comparisonSourceAvailability) {
    const records = await readInputSequenceRecords(
      file,
      'fasta'
    );
    records.forEach((record) => {
      sources.push({
        key: `homology:comparison:${displayIndex}:${record.recordId}`,
        recordId: record.recordId,
        aliases: record.aliases,
        sequence: record.sequence,
        origin: 'homology-comparison',
        sourceIndex: displayIndex
      });
    });
  }
  return sources;
};
