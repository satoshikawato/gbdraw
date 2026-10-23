// @ts-check
import { formatFastaEntry } from './feature-sequence-fasta.js';
import { recordSourceInterval } from './record-source-coordinates.js';
import {
  orderedConservationSources,
  orderedOptionalConservationFiles
} from './conservation-series.js';
import {
  gffFastaReaderSelector,
  readAllSequenceRecords,
  readComparisonSequenceRecords
} from './record-discovery.js';
import { linearCardRegion } from './record-options.js';
import { t } from '../utils/i18n.js';

/** @import { FeatureCatalogItem, FeatureCatalogRow } from './feature-catalog.js' */
/** @import { FileLike } from './session-resource-backing.js' */

/**
 * A registered sequence source: one record whose bases a match span reads.
 * @typedef {object} SequenceSource
 * @property {string} key
 * @property {string} recordId
 * @property {string[]} aliases
 * @property {string} sequence
 * @property {string} origin
 * @property {number | null} recordIndex
 * @property {number | null} sourceIndex
 */

/**
 * A source as callers supply it for registration: a row this module built, or
 * a sequence source of a feature catalog (a Session's is checked only for its
 * `sequence`). The registry reads every field as unvalidated and normalizes
 * it into a `SequenceSource`, or registers nothing.
 * @typedef {object} SequenceSourceInput
 * @property {unknown} [key]
 * @property {unknown} [recordId]
 * @property {unknown} [aliases]
 * @property {unknown} [sequence]
 * @property {unknown} [origin]
 * @property {unknown} [recordIndex]
 * @property {unknown} [sourceIndex]
 */

/**
 * The identity a match carries to pick its sequence source.
 * @typedef {object} SequenceSourceContext
 * @property {string} [origin]
 * @property {number | string | null} [recordIndex]
 * @property {number | string | null} [sourceIndex]
 */

/**
 * @typedef {object} SequenceSourceResolution
 * @property {SequenceSource | null} source
 * @property {string} reason
 * @property {boolean} [optionalUnavailable] The source was never supplied and is optional.
 */

/**
 * What a resolver returns, as `normalizeResolution` reads it: a resolution, or a bare source.
 * @typedef {SequenceSourceResolution | (SequenceSource & { source?: undefined }) | null | undefined} ResolverResult
 */

/**
 * The registry's sources by key, and the keys that two different sources claimed.
 * @typedef {Readonly<{ sources: Map<string, SequenceSource>, ambiguousKeys: Set<string> }>} SequenceSourceOwner
 */

/**
 * The identity fields a source is read for; every field may be absent or unnormalized.
 * @typedef {{ key?: unknown, recordId?: unknown, origin?: unknown, recordIndex?: unknown, sourceIndex?: unknown }} SourceIdentityFields
 */

/**
 * @typedef {object} SourceIdentity
 * @property {string} key
 * @property {string} recordId
 * @property {string} origin
 * @property {number | null} recordIndex
 * @property {number | null} sourceIndex
 */

/**
 * A record of a sequence file, in file order (`ParsedSequenceRecord`): `sequence` is null for a record without one.
 * @typedef {{ recordId: string, aliases: readonly string[], sequence: string | null }} InputSequenceRecord
 */

/**
 * Port that resolves a match span to its sequence source. A resolver may
 * return the bare source or a `{ source, reason }` resolution.
 * @typedef {(
 *   sourceKey: string,
 *   recordId: string,
 *   context: SequenceSourceContext
 * ) => SequenceSourceResolution | SequenceSource | null} SequenceSourceResolver
 */

/**
 * @typedef {object} MatchSpanExtraction
 * @property {boolean} valid
 * @property {string} [reason]
 * @property {number} [start]
 * @property {number} [end]
 * @property {string} [orientation]
 * @property {string} sequence
 * @property {number} [sequenceLength]
 */

/**
 * @typedef {object} MatchSequenceEntryOptions
 * @property {string} [matchId]
 * @property {SequenceSourceResolver} [resolveSequenceSource]
 * @property {SequenceSourceContext} [context]
 * @property {string} [unavailableReason]
 */

/**
 * @typedef {object} MatchSequenceBundleOptions
 * @property {string} [matchId]
 * @property {SequenceSourceResolver} [resolveSequenceSource]
 * @property {(span: any) => SequenceSourceContext} [contextForSpan]
 * @property {(span: any) => string} [unavailableReasonForSpan]
 */

/**
 * The slice of the session files that the sequence-source readers consult.
 * @typedef {Record<string, any>} MatchSequenceFiles
 */

/** @type {Readonly<Record<string, string>>} */
const IUPAC_COMPLEMENT = Object.freeze({
  A: 'T', C: 'G', G: 'C', T: 'A', U: 'A',
  R: 'Y', Y: 'R', S: 'S', W: 'W', K: 'M', M: 'K',
  B: 'V', D: 'H', H: 'D', V: 'B', N: 'N',
  '-': '-'
});

/** @param {unknown} value */
const text = (value) => String(value ?? '').trim();
/** @param {unknown} value */
const normalizedSequence = (value) => String(value ?? '').replace(/\s+/g, '').toUpperCase();
/**
 * @param {unknown} value
 * @returns {{ valid: boolean, supplied: boolean, value: number | null }}
 */
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
/**
 * @param {unknown} value
 * @param {string} [fallback]
 */
const safePart = (value, fallback = 'sequence') => {
  const cleaned = text(value).replace(/[^A-Za-z0-9_.-]+/g, '_').replace(/^_+|_+$/g, '');
  return cleaned || fallback;
};

/** @param {{ recordId?: unknown, aliases?: unknown } | null | undefined} source */
const aliasesForSource = (source) => {
  const values = [source?.recordId, ...(Array.isArray(source?.aliases) ? source.aliases : [])];
  /** @type {Set<string>} */
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

/**
 * @param {SequenceSource} left
 * @param {SequenceSource} right
 */
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

/** @param {unknown} sequence */
export const reverseComplementNucleotide = (sequence) => {
  const input = normalizedSequence(sequence);
  let output = '';
  for (let index = input.length - 1; index >= 0; index -= 1) {
    output += IUPAC_COMPLEMENT[input[index]] || 'N';
  }
  return output;
};

/**
 * @param {unknown} startRaw
 * @param {unknown} endRaw
 * @param {unknown} sequenceLength
 * @returns {{ valid: false, reason: string } | { valid: true, start: number, end: number, orientation: string }}
 */
export const validateMatchCoordinates = (startRaw, endRaw, sequenceLength) => {
  const start = Number(startRaw);
  const end = Number(endRaw);
  const length = Number(sequenceLength);
  if (!Number.isInteger(start) || !Number.isInteger(end)) {
    return { valid: false, reason: t('Coordinates must be whole 1-based values.') };
  }
  if (start < 1 || end < 1) {
    return { valid: false, reason: t('Coordinates must be greater than or equal to 1.') };
  }
  if (!Number.isInteger(length) || length < 1) {
    return { valid: false, reason: t('The sequence source is empty.') };
  }
  if (start > length || end > length) {
    return { valid: false, reason: t('Coordinates exceed the sequence length ({length} bp).', { length: length.toLocaleString() }) };
  }
  return { valid: true, start, end, orientation: start <= end ? '+' : '-' };
};

/**
 * @param {string} sequence
 * @param {unknown} startRaw
 * @param {unknown} endRaw
 * @returns {MatchSpanExtraction}
 */
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

/**
 * @param {SequenceSourceInput[]} [initialSources]
 */
export const createSequenceSourceRegistry = (initialSources = []) => {
  /** @type {SequenceSourceOwner} */
  let owner = Object.freeze({
    sources: new Map(),
    ambiguousKeys: new Set()
  });

  /**
   * @param {SequenceSourceOwner} targetOwner
   * @param {SequenceSourceInput | null | undefined} source
   * @param {{ trusted?: boolean }} [options]
   * @returns {SequenceSource | null}
   */
  const registerInto = (targetOwner, source, { trusted = false } = {}) => {
    const key = text(source?.key);
    if (trusted) {
      // A trusted source is an entry this registry produced (Artifact History).
      if (key) targetOwner.sources.set(key, /** @type {SequenceSource} */ (source));
      return key ? /** @type {SequenceSource} */ (source) : null;
    }
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
    if (targetOwner.ambiguousKeys.has(key)) return null;
    const existing = targetOwner.sources.get(key);
    if (existing) {
      if (sequenceSourcesAgree(existing, entry)) return existing;
      targetOwner.sources.delete(key);
      targetOwner.ambiguousKeys.add(key);
      return null;
    }
    targetOwner.sources.set(key, entry);
    return entry;
  };

  /** @param {SequenceSourceInput | null | undefined} source */
  const register = (source) => {
    const nextOwner = Object.freeze({
      sources: new Map(owner.sources),
      ambiguousKeys: new Set(owner.ambiguousKeys)
    });
    const registered = registerInto(nextOwner, source);
    owner = nextOwner;
    return registered;
  };

  /**
   * @param {readonly (SequenceSourceInput | null | undefined)[] | null} [nextSources]
   * @param {{ trusted?: boolean }} [options]
   * @returns {SequenceSourceOwner}
   */
  const buildOwner = (nextSources = [], { trusted = false } = {}) => {
    const nextOwner = Object.freeze({
      sources: new Map(),
      ambiguousKeys: new Set()
    });
    (Array.isArray(nextSources) ? nextSources : []).forEach((source) => {
      registerInto(nextOwner, source, { trusted });
    });
    return nextOwner;
  };

  /** @param {readonly (SequenceSourceInput | null | undefined)[] | null} [nextSources] */
  const reset = (nextSources = []) => {
    owner = buildOwner(nextSources);
  };

  // Artifact History only supplies entries previously produced by this registry.
  // Re-adopt those immutable entries without normalizing or copying large sequence
  // strings a second time.
  /** @param {readonly (SequenceSourceInput | null | undefined)[] | null} [nextSources] */
  const resetTrusted = (nextSources = []) => {
    owner = buildOwner(nextSources, { trusted: true });
  };

  /** @param {SequenceSourceOwner | null | undefined} nextOwner */
  const replaceTrustedOwner = (nextOwner) => {
    if (!(nextOwner?.sources instanceof Map) || !(nextOwner?.ambiguousKeys instanceof Set)) {
      throw new Error('A trusted sequence source owner is required.');
    }
    owner = nextOwner;
  };

  /**
   * @param {unknown} sourceKey
   * @param {unknown} recordId
   * @param {SequenceSourceContext} [context]
   * @returns {SequenceSourceResolution}
   */
  const resolve = (sourceKey, recordId, context = {}) => {
    const expectedSourceIndex = optionalNonnegativeIntegerStatus(context.sourceIndex);
    const expectedRecordIndex = optionalNonnegativeIntegerStatus(context.recordIndex);
    if (!expectedSourceIndex.valid || !expectedRecordIndex.valid) {
      return { source: null, reason: t('The sequence source identity for this match is invalid.') };
    }
    const expectedOrigin = text(context.origin);
    const wanted = text(recordId);
    /** @param {SequenceSource} source */
    const agrees = (source) => (
      (!expectedOrigin || source.origin === expectedOrigin)
      && (!expectedSourceIndex.supplied || source.sourceIndex === expectedSourceIndex.value)
      && (!expectedRecordIndex.supplied || source.recordIndex === expectedRecordIndex.value)
      && (!wanted || source.recordId === wanted || aliasesForSource(source).has(wanted))
    );
    const explicitKey = text(sourceKey);
    if (explicitKey) {
      if (owner.ambiguousKeys.has(explicitKey)) {
        return { source: null, reason: t('The sequence source key for this match is ambiguous.') };
      }
      const direct = owner.sources.get(explicitKey);
      if (!direct) {
        return { source: null, reason: t('The sequence source for this match is unavailable.') };
      }
      return agrees(direct)
        ? { source: direct, reason: '' }
        : { source: null, reason: t('The sequence source key conflicts with the match identity.') };
    }

    let candidates = Array.from(owner.sources.values()).filter((source) => {
      if (expectedOrigin && source.origin !== expectedOrigin) return false;
      if (expectedSourceIndex.supplied && source.sourceIndex !== expectedSourceIndex.value) return false;
      if (expectedRecordIndex.supplied && source.recordIndex !== expectedRecordIndex.value) return false;
      return true;
    });

    if (expectedRecordIndex.supplied && candidates.length === 1 && !wanted) {
      return { source: candidates[0], reason: '' };
    }
    if (!wanted) return { source: null, reason: t('The match does not identify a sequence record.') };
    const exact = candidates.filter((source) => source.recordId === wanted);
    if (exact.length === 1) return { source: exact[0], reason: '' };
    if (exact.length > 1) {
      return { source: null, reason: t('Record ID “{id}” is ambiguous in the available sequence sources.', { id: wanted }) };
    }
    const alias = candidates.filter((source) => aliasesForSource(source).has(wanted));
    if (alias.length === 1) return { source: alias[0], reason: '' };
    if (alias.length > 1) {
      return { source: null, reason: t('Record alias “{id}” is ambiguous in the available sequence sources.', { id: wanted }) };
    }
    return { source: null, reason: t('No sequence source matched record “{id}”.', { id: wanted }) };
  };

  owner = buildOwner(initialSources);
  return {
    get sources() { return owner.sources; },
    buildTrustedOwner: (/** @type {readonly (SequenceSourceInput | null | undefined)[] | null} */ nextSources = []) => buildOwner(nextSources, { trusted: true }),
    // E1: the owner `reset` would install, for a mode's artifact slot.
    buildSourceOwner: (/** @type {readonly (SequenceSourceInput | null | undefined)[]} */ nextSources) => buildOwner(nextSources),
    captureTrustedOwner: () => owner,
    register,
    replaceTrustedOwner,
    reset,
    resetTrusted,
    resolve,
    values: () => Array.from(owner.sources.values())
  };
};

/**
 * @param {SourceIdentityFields | null | undefined} source
 * @param {SourceIdentityFields} [overrides]
 * @returns {SourceIdentity}
 */
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

/** @param {SourceIdentity} identity */
const sourceIdentityKey = (identity) => JSON.stringify([
  identity.origin,
  identity.recordIndex,
  identity.sourceIndex,
  identity.recordId,
  identity.key
]);

/**
 * @param {Readonly<Record<string, unknown>> | null | undefined} value
 * @param {readonly string[]} keys
 * @param {string} [fallback]
 */
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

/**
 * @param {Readonly<Record<string, unknown>> | null | undefined} value
 * @param {readonly string[]} keys
 */
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

/**
 * @param {SourceIdentity} identity
 * @returns {SequenceSourceContext}
 */
const sourceResolutionContext = (identity) => ({
  origin: identity.origin,
  ...(identity.recordIndex === null ? {} : { recordIndex: identity.recordIndex }),
  ...(identity.sourceIndex === null ? {} : { sourceIndex: identity.sourceIndex })
});

/**
 * @param {{ sequence?: unknown, recordIndex?: unknown, sourceIndex?: unknown } | null | undefined} source
 * @param {SourceIdentity} identity
 */
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

/**
 * @param {{ files?: MatchSequenceFiles, circularConservation?: Record<string, any> }} [options]
 */
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

/**
 * @param {unknown} inventory
 * @param {number | null} sourceIndex
 * @returns {{ valid: false, reason: string } | { valid: true, availability: string }}
 */
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
/**
 * @param {{
 *   mode?: string,
 *   catalogFeatureState?: Record<string, any>,
 *   renderRequest?: Record<string, any>,
 *   comparisonSourceAvailability?: Record<string, any>[]
 * }} [options]
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
  /** @type {FeatureCatalogItem[] | null} */
  const items = Array.isArray(catalogFeatureState?.items)
    ? catalogFeatureState.items
    : null;
  /** @type {{ reason: string, itemIndex?: number, sequenceSourceIndex?: number, source?: SourceIdentity }[]} */
  const invalidCatalogSources = [];
  if (!['linear', 'circular'].includes(normalizedMode)) {
    invalidCatalogSources.push({ reason: 'The diagram mode is invalid.' });
  }
  if (!items) {
    invalidCatalogSources.push({ reason: 'The current feature catalog has no item list.' });
  }

  const catalogItems = items || [];
  /** @type {{ itemIndex: number, sequenceSourceIndex: number, source: FeatureCatalogRow }[]} */
  const sourceEntries = [];
  catalogItems.forEach((item, itemIndex) => {
    const sources = Array.isArray(item?.sequenceSources) ? item.sequenceSources : [];
    sources.forEach((source, sequenceSourceIndex) => {
      sourceEntries.push({ itemIndex, sequenceSourceIndex, source });
    });
  });
  // Catalog rows are untrusted: the registry reads a row with optional chaining and
  // skips one it cannot register, and the loop below reports each such row.
  const registry = createSequenceSourceRegistry(
    /** @type {SequenceSourceInput[]} */ (sourceEntries.map((entry) => entry.source))
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

  /**
   * @type {Map<string, {
   *   expectedSource: SourceIdentity, consumerCount: number, biologicalFeatureCount: number,
   *   matchEndpointCount: number, exampleConsumers: { kind: string, id: string }[]
   * }>}
   */
  const requirements = new Map();
  /** @type {Set<string>} */
  const failedRequirementKeys = new Set();
  /** @type {Map<string, Set<string>>} */
  const failureReasons = new Map();
  /** @type {Set<string>} */
  const optionalUnavailableRequirementKeys = new Set();
  /** @type {Map<string, Set<string>>} */
  const optionalUnavailableReasons = new Map();
  /** @type {Set<string>} */
  const usedDisplayedRecordKeys = new Set();
  const consumerCounts = { biologicalFeatures: 0, matchEndpoints: 0 };

  /**
   * @param {{
   *   kind: 'biological-feature' | 'match-endpoint', id: string, recordKey?: string,
   *   expectedSource: SourceIdentity, resolution?: SequenceSourceResolution
   * }} consumer
   */
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
    // The block above makes sure the key has its requirement.
    const requirement = /** @type {NonNullable<ReturnType<typeof requirements.get>>} */ (requirements.get(key));
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
      /** @type {Set<string>} */ (optionalUnavailableReasons.get(key)).add(
        text(resolution.reason) || t('The optional comparison FASTA was never supplied.')
      );
      return;
    }
    failedRequirementKeys.add(key);
    if (!failureReasons.has(key)) failureReasons.set(key, new Set());
    /** @type {Set<string>} */ (failureReasons.get(key)).add(
      text(resolution?.reason) || t('The required sequence source is unavailable.')
    );
  };

  catalogItems.forEach((item) => {
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
        /** @type {SequenceSourceResolution} */
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
          /** @type {SequenceSourceResolution} */
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
                // comparisonSourceAvailabilityAt gives every invalid result a reason.
                resolution = { source: null, reason: /** @type {string} */ (availability.reason) };
              } else if (availability.availability === 'never-supplied') {
                resolution = {
                  source: null,
                  optionalUnavailable: true,
                  reason: t('The optional comparison FASTA was never supplied.')
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

/**
 * @param {ResolverResult} resolved
 * @returns {SequenceSourceResolution}
 */
const normalizeResolution = (resolved) => {
  if (!resolved) return { source: null, reason: '' };
  if (resolved.source !== undefined) return resolved;
  return { source: resolved, reason: '' };
};

/**
 * @param {any} span
 * @param {MatchSequenceEntryOptions} [options]
 */
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
    displayRole: text(span?.displayRole) || (text(span?.role) === 'subject' ? t('Subject') : t('Query'))
  };
  /** @param {string} [reason] */
  const unavailable = (reason) => ({
    span: normalizedSpan,
    orientation: Number(normalizedSpan.start) <= Number(normalizedSpan.end) ? '+' : '-',
    sequenceLength: 0,
    fasta: '',
    filename: '',
    available: false,
    unavailableReason: reason || unavailableReason || t('Sequence is unavailable.')
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
  // The sequence is read on the displayed record; the header reports the
  // same bases in input-file coordinates (PD-OI-076).
  const source = recordSourceInterval(span?.sourceSpan, normalizedSpan.start, normalizedSpan.end);
  const coords = source || normalizedSpan;
  const strand = source ? (source.start <= source.end ? '+' : '-') : extracted.orientation;
  const header = [
    `${safePart(matchId, 'match')}_${normalizedSpan.role}`,
    `record=${recordId}`,
    `coords=${coords.start}..${coords.end}`,
    `strand=${strand}`
  ].join('|');
  const fasta = `${formatFastaEntry({ id: header, sequence: extracted.sequence })}\n`;
  return {
    span: { ...normalizedSpan, recordId },
    orientation: strand,
    sequenceLength: extracted.sequenceLength,
    fasta,
    filename: `${safePart(matchId)}_${normalizedSpan.role}_${safePart(recordId)}_${coords.start}-${coords.end}.fna`,
    available: true,
    unavailableReason: ''
  };
};

/**
 * @param {any[]} spans
 * @param {MatchSequenceBundleOptions} [options]
 */
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

/**
 * @param {readonly InputSequenceRecord[]} records
 * @param {unknown} selectorValue
 */
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

/**
 * @param {InputSequenceRecord | null} record
 * @param {{ region_start?: unknown, region_end?: unknown, region_reverse?: unknown } | null | undefined} sequenceState
 */
const materializeLinearSequenceRecord = (record, sequenceState) => {
  if (!record?.sequence) return null;
  let sequence = normalizedSequence(record.sequence);
  let reverse = Boolean(sequenceState?.region_reverse);
  // The region Generate draws: an end past the record end is drawn to it.
  const region = linearCardRegion(sequenceState ?? {}, 1, sequence.length);
  if (region && ('error' in region || !region.drawn)) return null;
  if (region?.drawn) {
    sequence = sequence.slice(region.drawn.start - 1, region.drawn.drawnEnd);
    reverse = reverse || region.start > region.end;
  }
  return {
    ...record,
    sequence: reverse ? reverseComplementNucleotide(sequence) : sequence
  };
};

/**
 * A Linear record file's or the Circular reference's records as the CLI loader reads them, once per File (OV-1680):
 * the browser reader's, or the Worker's for a file that reader declines (one that starts with a UTF-8 byte order
 * mark), so `#N` and the first record name the records the drawing reads (OV-1599).
 * @param {FileLike | null | undefined} file
 * @param {unknown} inputType
 * @returns {Promise<InputSequenceRecord[]>}
 */
const readInputSequenceRecords = async (file, inputType) => {
  if (!file) return [];
  const fmt = text(inputType).toLowerCase() === 'gb' ? 'genbank' : 'fasta';
  return (await readAllSequenceRecords(file, fmt))
    .map(({ id, aliases, sequence }) => ({ recordId: id, aliases, sequence }));
};

/**
 * A comparison genome's records as the ring search read them (`readComparisonSequenceRecords`, once per File): the
 * browser's for a FASTA file it reads exactly, the comparison reader's for the others. A file that cannot be read
 * gives no sources, and the other sources stay.
 * @param {FileLike | null} file null for a comparison file that was never supplied.
 * @returns {Promise<InputSequenceRecord[]>}
 */
const readComparisonRecords = async (file) => {
  if (!file) return [];
  try {
    return (await readComparisonSequenceRecords(file))
      .map(({ id, aliases, sequence }) => ({ recordId: id, aliases, sequence }));
  } catch {
    return [];
  }
};

/**
 * Rebuild the transient match-sequence registry after loading a web session.
 *
 * Sessions already contain the source files, so sequence strings do not need a
 * second persisted copy. Sources that were never supplied (for example, a
 * BLAST upload without its optional comparison FASTA) remain unavailable.
 */
/**
 * @param {{
 *   mode?: string,
 *   cInputType?: string,
 *   lInputType?: string,
 *   files?: MatchSequenceFiles,
 *   linearSeqs?: Record<string, any>[],
 *   circularConservation?: Record<string, any>
 * }} [options]
 */
export const buildRestoredMatchSequenceSources = async ({
  mode,
  cInputType,
  lInputType,
  files,
  linearSeqs = [],
  circularConservation
} = {}) => {
  /** @type {SequenceSourceInput[]} */
  const sources = [];
  if (text(mode).toLowerCase() === 'linear') {
    const inputType = text(lInputType).toLowerCase() === 'gff' ? 'fasta' : 'gb';
    for (let recordIndex = 0; recordIndex < (linearSeqs || []).length; recordIndex += 1) {
      const sequenceState = linearSeqs[recordIndex] || {};
      const file = inputType === 'gb' ? sequenceState.gb : sequenceState.fasta;
      const records = await readInputSequenceRecords(file, inputType);
      const selector = text(sequenceState.region_record_id);
      const selected = selectInputSequenceRecord(records, inputType === 'gb' ? selector
        : await gffFastaReaderSelector({ gffFile: sequenceState.gff, fastaFile: file, selector }));
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
    if (!record.sequence) return;
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
    const records = await readComparisonRecords(file);
    records.forEach((record) => {
      if (!record.sequence) return;
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
