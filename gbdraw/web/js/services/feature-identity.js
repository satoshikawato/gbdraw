// @ts-check
/** @typedef {Record<string, unknown>} AliasSource Any object whose alias keys are read. */
/** @typedef {{ valid: boolean, supplied: boolean, value: string }} TextAliasStatus */
/** @typedef {{ valid: boolean, supplied: boolean, value: number | null }} IntegerAliasStatus */
/**
 * @typedef {{
 *   valid: boolean,
 *   usable: boolean,
 *   recordIndex: IntegerAliasStatus,
 *   sourceIndex: IntegerAliasStatus,
 *   recordKey: TextAliasStatus,
 *   biologicalId: TextAliasStatus,
 *   renderedId: TextAliasStatus,
 *   stableId: TextAliasStatus
 * }} FeatureIdentity
 */
/**
 * A selector of a record in the committed request: by index or by record ID.
 * @typedef {{ kind?: string, index: number, value?: string }} CommittedRecordSelector
 */
/**
 * A record of the committed request, as the pure checks below, the record
 * display and the alignment read it. The committed request is of the current
 * schema, but the alignment owner's port (`getCommittedSession`) is not yet
 * typed, so this names only the fields they read rather than the generated
 * `CanonicalRecord`.
 * @typedef {{
 *   recordKey: string,
 *   cardinality?: string,
 *   selector?: CommittedRecordSelector | null,
 *   region?: { selector?: CommittedRecordSelector | null, reverseComplement?: unknown } | null,
 *   presentation?: { reverseComplement?: unknown } | null,
 *   display?: { isCircular: boolean | null, startCoordinate: number | null },
 *   source: { resourceId?: string, gffResourceId?: string, fastaResourceId?: string }
 * }} CommittedRecord
 */
/** @typedef {{ mode?: string, records: CommittedRecord[] }} CommittedRecordRequest */
/**
 * An input source of the record display (the File card and its optional pair).
 * @typedef {{ scope: string, sourceUid: string, source: unknown, paired?: unknown }} RecordDisplaySource
 */
/**
 * A feature as `isCurrentFeature` reads it. A catalog feature's `record_key` is text
 * (services/feature-catalog.js); a feature row of a Session before 40 is stored unchecked.
 * @typedef {{ record_key?: unknown }} RecordKeyedFeature
 */

/** @param {unknown} value */
const normalizeText = (value) => String(value ?? '').trim();

export const RECORD_INDEX_KEYS = Object.freeze([
  'recordIndex', 'record_index', 'recordIdx', 'record_idx', 'fileIdx', 'file_idx'
]);
export const SOURCE_FEATURE_INDEX_KEYS = Object.freeze([
  'featureIndex', 'feature_index', 'sourceFeatureIndex', 'source_feature_index'
]);
export const STABLE_FEATURE_ID_KEYS = Object.freeze([
  'stableFeatureSvgId', 'stable_feature_svg_id', 'stableFeatureId',
  'stable_feature_id', 'stableSvgId', 'stable_svg_id', 'featureSvgId',
  'feature_svg_id'
]);
export const RENDERED_FEATURE_ID_KEYS = Object.freeze([
  'renderedFeatureSvgId', 'rendered_feature_svg_id', 'renderedSvgId',
  'rendered_svg_id'
]);
export const RECORD_KEY_KEYS = Object.freeze(['recordKey', 'record_key']);
export const BIOLOGICAL_FEATURE_ID_KEYS = Object.freeze([
  'biologicalFeatureId', 'biological_feature_id'
]);
export const ORTHOGROUP_ID_KEYS = Object.freeze([
  'id', 'orthogroupId', 'orthogroup_id'
]);

/**
 * @param {AliasSource | null | undefined} source
 * @param {readonly string[]} keys
 * @returns {TextAliasStatus}
 */
export const textAliasStatus = (source, keys) => {
  let supplied = false;
  let value = '';
  for (const key of keys) {
    const normalized = normalizeText(source?.[key]);
    if (!normalized) continue;
    supplied = true;
    if (value && value !== normalized) return { valid: false, supplied, value: '' };
    value = normalized;
  }
  return { valid: true, supplied, value };
};

/**
 * @param {AliasSource | null | undefined} source
 * @param {readonly string[]} keys
 * @returns {IntegerAliasStatus}
 */
export const nonnegativeIntegerAliasStatus = (source, keys) => {
  let supplied = false;
  /** @type {number | null} */
  let value = null;
  for (const key of keys) {
    const raw = source?.[key];
    // A blank or whitespace-only value is not supplied, as in textAliasStatus
    // and Python's _consistent_nonnegative_integer_alias (OV-274).
    if (raw === null || raw === undefined || (typeof raw === 'string' && raw.trim() === '')) continue;
    supplied = true;
    const normalized = typeof raw === 'number'
      ? raw
      : (typeof raw === 'string' && /^\d+$/.test(raw.trim()) ? Number(raw.trim()) : Number.NaN);
    if (!Number.isSafeInteger(normalized) || normalized < 0) {
      return { valid: false, supplied, value: null };
    }
    if (value !== null && value !== normalized) {
      return { valid: false, supplied, value: null };
    }
    value = normalized;
  }
  return { valid: true, supplied, value };
};

/** @param {AliasSource | null | undefined} source */
export const orthogroupIdStatus = (source) => (
  textAliasStatus(source, ORTHOGROUP_ID_KEYS)
);

/**
 * @template G
 * @param {G[] | null | undefined} orthogroups
 * @returns {{ group: G, id: string }[]}
 */
export const uniqueOrthogroupEntries = (orthogroups) => {
  /** @type {Map<string, number>} */
  const counts = new Map();
  /** @type {{ group: G, id: string }[]} */
  const entries = [];
  (Array.isArray(orthogroups) ? orthogroups : []).forEach((group) => {
    // A group is a record of alias fields; callers hold it under their own shape.
    const status = orthogroupIdStatus(/** @type {AliasSource | null | undefined} */ (group));
    if (!status.valid || !status.supplied) return;
    counts.set(status.value, (counts.get(status.value) || 0) + 1);
    entries.push({ group, id: status.value });
  });
  return entries.filter(({ id }) => counts.get(id) === 1);
};

/**
 * @param {AliasSource | null | undefined} source
 * @param {{ allowLegacySvgStable?: boolean }} [options]
 * @returns {FeatureIdentity}
 */
export const featureIdentity = (source, { allowLegacySvgStable = false } = {}) => {
  const recordIndex = nonnegativeIntegerAliasStatus(source, RECORD_INDEX_KEYS);
  const sourceIndex = nonnegativeIntegerAliasStatus(source, SOURCE_FEATURE_INDEX_KEYS);
  const recordKey = textAliasStatus(source, RECORD_KEY_KEYS);
  const biologicalId = textAliasStatus(source, BIOLOGICAL_FEATURE_ID_KEYS);
  const renderedId = textAliasStatus(source, RENDERED_FEATURE_ID_KEYS);
  let stableId = textAliasStatus(source, STABLE_FEATURE_ID_KEYS);
  if (allowLegacySvgStable && !renderedId.supplied) {
    const legacySvgId = textAliasStatus(source, ['svg_id', 'svgId']);
    if (
      !legacySvgId.valid ||
      (stableId.supplied && legacySvgId.supplied && stableId.value !== legacySvgId.value)
    ) {
      stableId = { valid: false, supplied: true, value: '' };
    } else if (!stableId.supplied) {
      stableId = legacySvgId;
    }
  }
  const valid = [recordIndex, sourceIndex, recordKey, biologicalId, renderedId, stableId]
    .every((status) => status.valid)
    && recordKey.supplied === biologicalId.supplied;
  const scoped = recordIndex.supplied || recordKey.supplied;
  const identified = stableId.supplied || sourceIndex.supplied || biologicalId.supplied;
  return {
    valid,
    usable: valid && scoped && identified,
    recordIndex,
    sourceIndex,
    recordKey,
    biologicalId,
    renderedId,
    stableId
  };
};

/**
 * @param {AliasSource | null | undefined} feature
 * @returns {FeatureIdentity}
 */
export const renderedFeatureIdentity = (feature) => {
  const identity = featureIdentity(feature);
  const renderedId = textAliasStatus(feature, [
    ...RENDERED_FEATURE_ID_KEYS,
    'svg_id',
    'svgId'
  ]);
  return {
    ...identity,
    valid: identity.valid && renderedId.valid,
    usable: identity.usable && renderedId.valid && renderedId.supplied,
    renderedId
  };
};

/**
 * @param {FeatureIdentity | null | undefined} reference
 * @param {FeatureIdentity | null | undefined} candidate
 * @param {{ includeRendered?: boolean }} [options]
 */
export const identityMatches = (reference, candidate, { includeRendered = false } = {}) => {
  if (!reference?.usable || !candidate?.usable) return false;
  /** @type {('recordIndex' | 'sourceIndex' | 'recordKey' | 'biologicalId' | 'stableId' | 'renderedId')[]} */
  const fields = ['recordIndex', 'sourceIndex', 'recordKey', 'biologicalId', 'stableId'];
  if (includeRendered) fields.push('renderedId');
  return fields.every((field) => (
    !reference[field].supplied || (
      candidate[field].supplied && candidate[field].value === reference[field].value
    )
  ));
};

/** @typedef {'recordIndex' | 'sourceIndex' | 'recordKey' | 'biologicalId' | 'stableId' | 'renderedId'} IdentityField */

/**
 * Index entries by the fields their identities supply, and return a lookup of
 * the entries a reference identity matches (`identityMatches`, which still owns
 * the matching rule; the index only narrows the candidates). Entries whose
 * identity is not usable never match.
 * @template T
 * @param {readonly T[]} entries
 * @param {(entry: T) => FeatureIdentity} identityOf
 * @param {{ includeRendered?: boolean }} [options]
 * @returns {(reference: FeatureIdentity) => T[]}
 */
export const indexFeatureIdentities = (entries, identityOf, { includeRendered = true } = {}) => {
  /** @type {IdentityField[]} */
  const fields = ['recordIndex', 'sourceIndex', 'recordKey', 'biologicalId', 'stableId'];
  if (includeRendered) fields.push('renderedId');
  /** @type {Map<IdentityField, Map<string | number | null, T[]>>} */
  const indexes = new Map(fields.map((field) => [field, new Map()]));
  for (const entry of entries) {
    const identity = identityOf(entry);
    if (!identity.usable) continue;
    for (const field of fields) {
      const status = identity[field];
      if (!status.supplied) continue;
      // `fields` lists every key of `indexes`, so the lookup always finds a map.
      const index = /** @type {Map<string | number | null, T[]>} */ (indexes.get(field));
      const bucket = index.get(status.value);
      if (bucket) bucket.push(entry);
      else index.set(status.value, [entry]);
    }
  }
  return (reference) => {
    if (!reference.usable) return [];
    let candidates = entries;
    for (const field of fields) {
      const status = reference[field];
      if (!status.supplied) continue;
      const matches = /** @type {Map<string | number | null, T[]>} */ (indexes.get(field)).get(status.value) || [];
      if (matches.length < candidates.length) candidates = matches;
    }
    return candidates.filter((entry) => identityMatches(reference, identityOf(entry), { includeRendered }));
  };
};

/**
 * The lookup of the one member of `members` that is a feature: both identities
 * supply the same fields with the same values, the rendered ID included, and no
 * other member does. The members are indexed once, so a caller that resolves
 * many features against one group builds the lookup once (OV-914).
 * @template M
 * @param {M[] | null | undefined} members
 * @returns {(feature: AliasSource | null | undefined) => M | null}
 */
export const uniqueOrthogroupMemberResolver = (members) => {
  // A member is a record of alias fields; callers hold it under their own shape.
  const entries = (Array.isArray(members) ? members : []).map((member) => ({
    member,
    identity: featureIdentity(/** @type {AliasSource | null | undefined} */ (member))
  }));
  const candidatesFor = indexFeatureIdentities(entries, (entry) => entry.identity);
  return (feature) => {
    const reference = renderedFeatureIdentity(feature);
    if (!reference.usable) return null;
    const matches = candidatesFor(reference)
      .filter(({ identity }) => identityMatches(identity, reference, { includeRendered: true }));
    return matches.length === 1 ? matches[0].member : null;
  };
};

/**
 * @template M
 * @param {AliasSource | null | undefined} feature
 * @param {M[] | null | undefined} [members]
 * @returns {M | null}
 */
export const resolveUniqueOrthogroupMemberForFeature = (feature, members = []) => (
  uniqueOrthogroupMemberResolver(members)(feature)
);

// The committed request record that draws a record key: the record itself, or
// one record of a source drawn whole (`cardinality: 'all'`, `<recordKey>:<n>`).
/**
 * @param {CommittedRecordRequest | null | undefined} request
 * @param {string} key
 */
export const committedRecordForKey = (request, key) => request?.records.find((record) => record.recordKey === key
  || (record.cardinality === 'all' && key.startsWith(`${record.recordKey}:`)
    && /^[1-9]\d*$/.test(key.slice(record.recordKey.length + 1))));

// Whether an input source of the record display feeds a committed record.
/**
 * @param {CommittedRecordRequest | null | undefined} request
 * @param {CommittedRecord} record
 * @param {RecordDisplaySource} source
 */
export const committedRecordUsesSource = (request, record, source) => source.scope === request?.mode
  && (source.scope === 'circular' || record.recordKey === source.sourceUid
    || (record.selector?.kind === 'recordIndex'
      && record.recordKey === `${source.sourceUid}:${record.selector.index + 1}`));

// Whether a feature's record still reads the input its committed Result was
// drawn from. A pure check (R13): the composition root passes the record
// display's source binding, refreshed for the committed request: `request`,
// the current input `sources`, the `boundSources` captured when that request
// was committed, and `matchesSavedSource(file, resourceId)` for its Session.
/**
 * @param {RecordKeyedFeature} feature
 * @param {{
 *   request: CommittedRecordRequest | null | undefined,
 *   sources: RecordDisplaySource[],
 *   boundSources: RecordDisplaySource[],
 *   matchesSavedSource: (file: unknown, resourceId: string | undefined) => boolean
 * }} binding
 */
export const isCurrentFeature = (feature, { request, sources, boundSources, matchesSavedSource }) => {
  // A record key that is not text names no committed record (OV-618).
  const key = feature.record_key;
  const record = typeof key === 'string' ? committedRecordForKey(request, key) : undefined;
  if (!record) return false;
  return sources.some((source) => source.source && committedRecordUsesSource(request, record, source)
    && boundSources.some((bound) => bound.scope === source.scope && bound.sourceUid === source.sourceUid
      && bound.source === source.source && bound.paired === source.paired)
    && matchesSavedSource(source.source, record.source.resourceId || record.source.gffResourceId)
    && (!source.paired || matchesSavedSource(source.paired, record.source.fastaResourceId)));
};
