// @ts-check
import {
  formatFeatureLocation,
  isInternalProteinDisplayId,
  resolveDisplayProteinId,
  resolveInternalProteinId
} from '../../services/feature-utils.js';
import { uniqueOrthogroupMemberResolver } from '../../services/feature-identity.js';
import { msg, t, tm } from '../../utils/i18n.js';
/** @import { FeatureLike } from '../../services/feature-utils.js' */

/**
 * A similarity group as the search reads it, or an unvalidated row: its ID under
 * either name, the display name under any of three, and its member rows.
 * @typedef {object} SearchOrthogroup
 * @property {unknown} [id]
 * @property {unknown} [orthogroupId]
 * @property {unknown} [orthogroup_id]
 * @property {unknown} [display_name]
 * @property {unknown} [displayName]
 * @property {unknown} [name]
 * @property {unknown} [description]
 * @property {FeatureLike[] | null} [members]
 */

/**
 * A value the search can match, under the label that names where it came from.
 * @typedef {object} SearchItem
 * @property {string} label
 * @property {string} value
 * @property {string} alphabet `nucleotide` or `amino-acid` for a sequence, else empty.
 */

/**
 * @typedef {object} PreparedSearchItem
 * @property {string} label
 * @property {string} value
 * @property {string} normalizedValue Empty for a sequence, which a pattern matches.
 * @property {string} alphabet
 */

/**
 * @typedef {object} SearchMatchDetail
 * @property {string} label
 * @property {string} value
 * @property {string} match
 */

/**
 * @typedef {object} FeatureSearchMatcher
 * @property {boolean} active
 * @property {string} error
 * @property {(value: unknown, alphabet: string) => string} match
 * @property {(item: PreparedSearchItem) => string} [matchPrepared]
 * @property {(values: unknown[]) => boolean} test
 */

const RICH_FEATURE_SEARCH_FIELD_IDS = Object.freeze([
  'qualifier-key',
  'qualifier-value',
  'nucleotide',
  'amino-acid'
]);

const FEATURE_SEARCH_FIELD_DEFINITIONS = Object.freeze([
  { id: 'all', label: msg('All') },
  { id: 'label', label: msg('Label') },
  { id: 'type', label: msg('Feature type') },
  { id: 'record-id', label: msg('Record ID') },
  { id: 'location', label: msg('Location') },
  { id: 'strand', label: msg('Strand') },
  { id: 'orthogroup', label: msg('Similarity group') },
  { id: 'qualifier-key', label: msg('Qualifier key'), richOnly: true },
  { id: 'qualifier-value', label: msg('Qualifier value'), richOnly: true },
  { id: 'nucleotide', label: msg('Nucleotide'), richOnly: true },
  { id: 'amino-acid', label: msg('Amino acid'), richOnly: true }
]);

const FIELD_IDS = new Set(FEATURE_SEARCH_FIELD_DEFINITIONS.map((field) => field.id));
const RICH_FIELD_IDS = new Set(RICH_FEATURE_SEARCH_FIELD_IDS);

/** @param {unknown} field */
export const isRichFeatureSearchField = (field) => RICH_FIELD_IDS.has(String(field || ''));

/** @param {unknown} popupMode */
export const normalizeFeatureSearchPopupMode = (popupMode) => (
  popupMode === 'simple' ? 'simple' : 'rich'
);

/**
 * @param {unknown} field
 * @param {{ popupMode?: unknown }} [options]
 */
export const normalizeFeatureSearchField = (field, { popupMode = 'rich' } = {}) => {
  const selected = String(field || 'all').trim() || 'all';
  if (!FIELD_IDS.has(selected)) return 'all';
  if (normalizeFeatureSearchPopupMode(popupMode) === 'simple' && isRichFeatureSearchField(selected)) {
    return 'all';
  }
  return selected;
};

/** @param {{ popupMode?: unknown }} [options] */
export const getFeatureSearchFieldOptions = ({ popupMode = 'rich' } = {}) => {
  const normalizedPopupMode = normalizeFeatureSearchPopupMode(popupMode);
  return FEATURE_SEARCH_FIELD_DEFINITIONS.filter((field) => (
    normalizedPopupMode !== 'simple' || !field.richOnly
  )).map((field) => ({ ...field, label: tm(field.label) }));
};

/** @param {unknown} value */
const normalizeSearchText = (value) => String(value == null ? '' : value).toLowerCase();

/**
 * @param {{ label?: unknown, value?: unknown, alphabet?: unknown } | null | undefined} item
 * @returns {PreparedSearchItem}
 */
const prepareSearchItem = (item) => ({
  label: String(item?.label || ''),
  value: String(item?.value || ''),
  normalizedValue: item?.alphabet ? '' : normalizeSearchText(item?.value),
  alphabet: String(item?.alphabet || '')
});

/** @type {Readonly<Record<string, string>>} */
const NUCLEOTIDE_IUPAC = Object.freeze({
  A: 'A',
  C: 'C',
  G: 'G',
  T: 'TU',
  U: 'TU',
  R: 'AG',
  Y: 'CTU',
  S: 'GC',
  W: 'ATU',
  K: 'GTU',
  M: 'AC',
  B: 'CGTU',
  D: 'AGTU',
  H: 'ACTU',
  V: 'ACG',
  N: 'ACGTU',
  '-': '-'
});

/** @type {Readonly<Record<string, string>>} */
const AMINO_ACID_IUPAC = Object.freeze({
  A: 'A',
  C: 'C',
  D: 'D',
  E: 'E',
  F: 'F',
  G: 'G',
  H: 'H',
  I: 'I',
  K: 'K',
  L: 'L',
  M: 'M',
  N: 'N',
  P: 'P',
  Q: 'Q',
  R: 'R',
  S: 'S',
  T: 'T',
  V: 'V',
  W: 'W',
  Y: 'Y',
  U: 'U',
  O: 'O',
  B: 'DN',
  Z: 'EQ',
  J: 'IL',
  X: 'ACDEFGHIKLMNPQRSTVWYUO',
  '*': '*',
  '-': '-'
});

/** @param {unknown} alphabet */
const getIupacAlphabet = (alphabet) => {
  if (alphabet === 'nucleotide') return NUCLEOTIDE_IUPAC;
  if (alphabet === 'amino-acid') return AMINO_ACID_IUPAC;
  return null;
};

/** @param {unknown} value */
const uniqueCharacters = (value) => {
  /** @type {Set<string>} */
  const seen = new Set();
  /** @type {string[]} */
  const chars = [];
  String(value || '').split('').forEach((char) => {
    if (!char || seen.has(char)) return;
    seen.add(char);
    chars.push(char);
  });
  return chars;
};

/** @param {unknown} left @param {unknown} right */
const charSetsIntersect = (left, right) => {
  const rightSet = new Set(uniqueCharacters(right));
  return uniqueCharacters(left).some((char) => rightSet.has(char));
};

/** @param {unknown} value */
const escapeRegExpText = (value) => String(value).replace(/[.*+?^${}()|[\]\\]/g, '\\$&');
/** @param {unknown} value */
const escapeRegExpClassChar = (value) => String(value).replace(/[\\\]\[\^-]/g, '\\$&');

/** @param {unknown} query @param {unknown} alphabet */
const buildIupacQueryPattern = (query, alphabet) => {
  const map = getIupacAlphabet(alphabet);
  const normalized = String(query || '').replace(/\s+/g, '').toUpperCase();
  if (!map || !normalized) return null;
  const targetChars = Object.keys(map);
  /** @type {string[]} */
  const parts = [];
  normalized.split('').forEach((queryChar) => {
    const querySet = map[queryChar] || queryChar;
    const matchingTargets = targetChars.filter((targetChar) => (
      charSetsIntersect(querySet, map[targetChar] || targetChar)
    ));
    if (!matchingTargets.length) {
      parts.push(escapeRegExpText(queryChar));
    } else if (matchingTargets.length === 1) {
      parts.push(escapeRegExpText(matchingTargets[0]));
    } else {
      parts.push(`[${matchingTargets.map(escapeRegExpClassChar).join('')}]`);
    }
  });
  return parts.join('');
};

/**
 * @param {string[]} values
 * @param {unknown} value
 */
const appendSearchValues = (values, value) => {
  if (Array.isArray(value)) {
    value.forEach((entry) => appendSearchValues(values, entry));
    return;
  }
  if (value === null || value === undefined) return;
  const text = String(value).trim();
  if (text) values.push(text);
};

/**
 * @param {SearchItem[]} items
 * @param {string} label
 * @param {unknown} value
 * @param {{ alphabet?: unknown }} [options]
 */
const appendSearchItems = (items, label, value, options = {}) => {
  if (Array.isArray(value)) {
    value.forEach((entry) => appendSearchItems(items, label, entry, options));
    return;
  }
  if (value === null || value === undefined) return;
  const text = String(value).trim();
  if (!text || isInternalProteinDisplayId(text)) return;
  items.push({
    label,
    value: text,
    alphabet: options.alphabet ? String(options.alphabet) : ''
  });
};

/**
 * @param {FeatureLike | null | undefined} feature
 * @returns {Record<string, unknown>}
 */
const getFeatureQualifiers = (feature) => (
  feature?.qualifiers && typeof feature.qualifiers === 'object' && !Array.isArray(feature.qualifiers)
    // Checked to be a non-array object just above.
    ? /** @type {Record<string, unknown>} */ (feature.qualifiers)
    : {}
);

/** @param {FeatureLike | null | undefined} feature */
const resolveFeatureAminoAcidSequence = (feature) => {
  const direct = String(feature?.amino_acid_sequence || feature?.aminoAcidSequence || '').trim();
  if (direct) return direct;
  const qualifiers = getFeatureQualifiers(feature);
  /** @type {string[]} */
  const values = [];
  appendSearchValues(values, qualifiers.translation);
  return values.find((value) => value.trim()) || '';
};

/** @param {FeatureLike | null | undefined} feature @param {unknown} qualifierKey */
const getQualifierSearchItems = (feature, qualifierKey) => {
  const target = normalizeSearchText(qualifierKey).trim();
  const qualifiers = getFeatureQualifiers(feature);
  /** @type {SearchItem[]} */
  const items = [];
  Object.keys(qualifiers).sort().forEach((key) => {
    if (target && normalizeSearchText(key) !== target) return;
    appendSearchItems(items, t('Qualifier {key}', { key }), qualifiers[key]);
  });
  return items;
};

/** @param {FeatureLike | null | undefined} feature @param {string} key */
const firstQualifierValue = (feature, key) => {
  const qualifiers = getFeatureQualifiers(feature);
  /** @type {string[]} */
  const values = [];
  appendSearchValues(values, qualifiers[key]);
  return values.find((value) => value.trim()) || '';
};

/** @param {FeatureLike | null | undefined} feature */
const getFeatureLabelCandidates = (feature) => [
  feature?.display_label,
  feature?.displayLabel,
  feature?.label,
  feature?.gene,
  feature?.locus_tag,
  feature?.locusTag,
  feature?.product,
  firstQualifierValue(feature, 'gene'),
  firstQualifierValue(feature, 'locus_tag'),
  firstQualifierValue(feature, 'product'),
  feature?.search_labels,
  feature?.searchLabels,
  feature?.svg_id
];

/** @param {FeatureLike | null | undefined} feature */
const getLabelSearchItems = (feature) => {
  /** @type {SearchItem[]} */
  const items = [];
  /** @type {Set<string>} */
  const seen = new Set();
  getFeatureLabelCandidates(feature).forEach((value) => {
    /** @type {string[]} */
    const values = [];
    appendSearchValues(values, value);
    values.forEach((entry) => {
      const key = normalizeSearchText(entry);
      if (!key || seen.has(key)) return;
      seen.add(key);
      appendSearchItems(items, key === normalizeSearchText(feature?.svg_id) ? t('SVG ID') : t('Label'), entry);
    });
  });
  return items;
};

/** @param {{ id?: unknown, orthogroupId?: unknown, orthogroup_id?: unknown } | null | undefined} source */
const orthogroupIdFor = (source) => {
  const ids = new Set([
    source?.id,
    source?.orthogroupId,
    source?.orthogroup_id
  ].map((value) => String(value || '').trim()).filter(Boolean));
  return ids.size === 1 ? ids.values().next().value : '';
};

/** @param {unknown} orthogroups Rows of similarity groups; anything but an array holds none. */
const buildOrthogroupMap = (orthogroups) => {
  /** @type {Map<string, SearchOrthogroup[]>} */
  const candidates = new Map();
  (Array.isArray(orthogroups) ? orthogroups : []).forEach((/** @type {SearchOrthogroup} */ group) => {
    const id = orthogroupIdFor(group);
    if (!id) return;
    const matches = candidates.get(id) || [];
    matches.push(group);
    candidates.set(id, matches);
  });
  /** @type {Map<string, SearchOrthogroup>} */
  const groups = new Map();
  candidates.forEach((matches, id) => {
    if (matches.length === 1) groups.set(id, matches[0]);
  });
  return groups;
};

/** @param {FeatureLike | null | undefined} feature */
const getOrthogroupId = (feature) => orthogroupIdFor({
  orthogroupId: feature?.orthogroupId,
  orthogroup_id: feature?.orthogroup_id
});

// Each group's members are indexed once per group map, that is once per search
// index build, not once per feature and field (OV-914).
/** @type {WeakMap<Map<string, SearchOrthogroup>, Map<SearchOrthogroup, (feature: FeatureLike | null | undefined) => FeatureLike | null>>} */
const memberResolversByGroupMap = new WeakMap();

/**
 * @param {FeatureLike | null | undefined} feature
 * @param {SearchOrthogroup | null} group
 * @param {Map<string, SearchOrthogroup>} orthogroupsById
 */
const getFeatureOrthogroupMember = (feature, group, orthogroupsById) => {
  if (!group) return null;
  let resolvers = memberResolversByGroupMap.get(orthogroupsById);
  if (!resolvers) {
    resolvers = new Map();
    memberResolversByGroupMap.set(orthogroupsById, resolvers);
  }
  let resolve = resolvers.get(group);
  if (!resolve) {
    resolve = uniqueOrthogroupMemberResolver(group.members);
    resolvers.set(group, resolve);
  }
  return resolve(feature);
};

/** @param {FeatureLike | null | undefined} feature @param {Map<string, SearchOrthogroup>} orthogroupsById */
const getOrthogroupSearchItems = (feature, orthogroupsById) => {
  const orthogroupId = getOrthogroupId(feature);
  const candidateGroup = orthogroupsById.get(orthogroupId) || null;
  const group = orthogroupIdFor(candidateGroup) === orthogroupId ? candidateGroup : null;
  const member = getFeatureOrthogroupMember(feature, group, orthogroupsById);
  const proteinId = resolveDisplayProteinId(feature, member);
  const internalId = resolveInternalProteinId(feature, member);
  /** @type {SearchItem[]} */
  const items = [];
  appendSearchItems(items, t('Similarity group ID'), orthogroupId);
  appendSearchItems(items, t('Similarity group name'), group?.display_name || group?.displayName || group?.name);
  appendSearchItems(items, t('Similarity group description'), group?.description);
  appendSearchItems(items, t('Protein ID'), proteinId);
  if (
    internalId &&
    internalId !== proteinId &&
    !isInternalProteinDisplayId(internalId)
  ) {
    appendSearchItems(items, t('Internal protein ID'), internalId);
  }
  appendSearchItems(items, t('Similarity-group member'), member?.label);
  appendSearchItems(items, t('Similarity-group member gene'), member?.gene);
  appendSearchItems(items, t('Similarity-group member product'), member?.product);
  appendSearchItems(items, t('Similarity-group member note'), member?.note);
  appendSearchItems(items, t('Similarity-group member protein ID'), resolveDisplayProteinId(null, member));
  return items;
};

/**
 * @param {FeatureLike | null | undefined} feature
 * @param {unknown} field
 * @param {unknown} qualifierKey
 * @param {{ popupMode?: unknown, orthogroupsById?: Map<string, SearchOrthogroup> }} [options]
 * @returns {SearchItem[]}
 */
const featureSearchItems = (
  feature,
  field,
  qualifierKey,
  {
    popupMode = 'rich',
    orthogroupsById = new Map()
  } = {}
) => {
  /** @type {SearchItem[]} */
  const items = [];
  const selectedField = normalizeFeatureSearchField(field, { popupMode });
  const qualifiers = getFeatureQualifiers(feature);

  if (selectedField === 'label') return getLabelSearchItems(feature);
  if (selectedField === 'type') {
    appendSearchItems(items, t('Feature type'), feature?.type);
    return items;
  }
  if (selectedField === 'record-id') {
    appendSearchItems(items, t('Record ID'), feature?.record_id);
    appendSearchItems(items, t('Record ID'), feature?.recordId);
    appendSearchItems(items, t('Record ID'), feature?.displayRecordId);
    return items;
  }
  if (selectedField === 'location') {
    appendSearchItems(items, t('Location'), formatFeatureLocation(feature));
    return items;
  }
  if (selectedField === 'strand') {
    appendSearchItems(items, t('Strand'), feature?.strand);
    return items;
  }
  if (selectedField === 'orthogroup') {
    return getOrthogroupSearchItems(feature, orthogroupsById);
  }
  if (selectedField === 'qualifier-key') {
    appendSearchItems(items, t('Qualifier key'), Object.keys(qualifiers));
    return items;
  }
  if (selectedField === 'qualifier-value') {
    return getQualifierSearchItems(feature, qualifierKey);
  }
  if (selectedField === 'nucleotide') {
    appendSearchItems(items, t('Nucleotide sequence'), feature?.nucleotide_sequence || feature?.nucleotideSequence, { alphabet: 'nucleotide' });
    return items;
  }
  if (selectedField === 'amino-acid') {
    appendSearchItems(items, t('Amino acid sequence'), resolveFeatureAminoAcidSequence(feature), { alphabet: 'amino-acid' });
    return items;
  }

  items.push(...getLabelSearchItems(feature));
  appendSearchItems(items, t('Record ID'), feature?.record_id);
  appendSearchItems(items, t('Record ID'), feature?.recordId);
  appendSearchItems(items, t('Record ID'), feature?.displayRecordId);
  appendSearchItems(items, t('Feature type'), feature?.type);
  appendSearchItems(items, t('Location'), formatFeatureLocation(feature));
  appendSearchItems(items, t('Strand'), feature?.strand);
  items.push(...getOrthogroupSearchItems(feature, orthogroupsById));
  if (normalizeFeatureSearchPopupMode(popupMode) !== 'simple') {
    // All excludes sequence content; the Nucleotide and Amino acid fields search it (D-01).
    appendSearchItems(items, t('Qualifier key'), Object.keys(qualifiers));
    Object.keys(qualifiers).sort().forEach((key) => {
      if (key.toLowerCase() !== 'translation') appendSearchItems(items, t('Qualifier {key}', { key }), qualifiers[key]);
    });
  }
  return items;
};

/**
 * @param {unknown} query
 * @param {unknown} useRegex
 * @returns {FeatureSearchMatcher}
 */
const compileFeatureSearchMatcher = (query, useRegex) => {
  const trimmedQuery = String(query || '').trim();
  if (!trimmedQuery) {
    return { active: false, error: '', match: () => '', test: () => false };
  }
  if (useRegex) {
    try {
      const regex = new RegExp(trimmedQuery, 'i');
      return {
        active: true,
        error: '',
        match: (value) => {
          regex.lastIndex = 0;
          const match = String(value == null ? '' : value).match(regex);
          return match ? String(match[0] || '') : '';
        },
        matchPrepared: (item) => {
          regex.lastIndex = 0;
          const match = String(item?.value || '').match(regex);
          return match ? String(match[0] || '') : '';
        },
        test: (values) => values.some((value) => {
          regex.lastIndex = 0;
          return regex.test(String(value == null ? '' : value));
        })
      };
    } catch {
      return { active: true, error: 'Invalid JavaScript regular expression. Turn off Regular expression to return to word search.', match: () => '', test: () => false };
    }
  }

  const needle = normalizeSearchText(trimmedQuery);
  /** @type {Record<string, RegExp | null>} */
  const sequenceMatchers = {};
  /** @param {string} alphabet */
  const getSequenceRegex = (alphabet) => {
    if (!getIupacAlphabet(alphabet)) return null;
    if (!Object.prototype.hasOwnProperty.call(sequenceMatchers, alphabet)) {
      const pattern = buildIupacQueryPattern(trimmedQuery, alphabet);
      sequenceMatchers[alphabet] = pattern ? new RegExp(pattern, 'i') : null;
    }
    return sequenceMatchers[alphabet];
  };
  return {
    active: true,
    error: '',
    match: (value, alphabet) => {
      const text = String(value == null ? '' : value);
      const sequenceRegex = getSequenceRegex(alphabet);
      if (sequenceRegex) {
        sequenceRegex.lastIndex = 0;
        const sequenceMatch = text.match(sequenceRegex);
        return sequenceMatch ? String(sequenceMatch[0] || '') : '';
      }
      const index = normalizeSearchText(text).indexOf(needle);
      return index === -1 ? '' : text.slice(index, index + trimmedQuery.length);
    },
    matchPrepared: (item) => {
      const text = String(item?.value || '');
      const sequenceRegex = getSequenceRegex(item?.alphabet);
      if (sequenceRegex) {
        sequenceRegex.lastIndex = 0;
        const sequenceMatch = text.match(sequenceRegex);
        return sequenceMatch ? String(sequenceMatch[0] || '') : '';
      }
      const index = String(item?.normalizedValue || '').indexOf(needle);
      return index === -1 ? '' : text.slice(index, index + trimmedQuery.length);
    },
    test: (values) => values.some((value) => normalizeSearchText(value).includes(needle))
  };
};

/** @param {FeatureLike | null | undefined} feature @param {string} field @param {unknown} popupMode @param {Map<string, SearchOrthogroup>} orthogroupsById */
const preparedItemsForField = (feature, field, popupMode, orthogroupsById) => (
  featureSearchItems(feature, field, '', { popupMode, orthogroupsById }).map(prepareSearchItem)
);

/**
 * The prepared index that `buildFeatureSearchIndex` returns and
 * `runFeatureSearch` reads.
 * @typedef {object} FeatureSearchIndex
 * @property {string} popupMode
 * @property {string[]} featureOrder
 * @property {Map<string, FeatureSearchDocument>} byId
 * @property {Map<string, string[]>} qualifierFeatureIdsByKey
 * @property {Map<string, SearchOrthogroup>} orthogroupsById
 */

/**
 * One searchable feature: its prepared items by search field, and its qualifier values by key.
 * @typedef {object} FeatureSearchDocument
 * @property {FeatureLike} feature
 * @property {Map<string, PreparedSearchItem[]>} itemsByField
 * @property {Map<string, PreparedSearchItem[]>} qualifierValuesByKey
 */

/**
 * @typedef {object} BuildFeatureSearchIndexOptions
 * @property {FeatureLike[]} [features]
 * @property {string} [popupMode]
 * @property {SearchOrthogroup[]} [orthogroups]
 */

/**
 * @param {BuildFeatureSearchIndexOptions} [options]
 * @returns {FeatureSearchIndex}
 */
export const buildFeatureSearchIndex = ({
  features = [],
  popupMode = 'rich',
  orthogroups = []
} = {}) => {
  const normalizedPopupMode = normalizeFeatureSearchPopupMode(popupMode);
  const orthogroupsById = buildOrthogroupMap(orthogroups);
  /** @type {string[]} */
  const featureOrder = [];
  /** @type {Map<string, FeatureSearchDocument>} */
  const byId = new Map();
  /** @type {Map<string, string[]>} */
  const qualifierFeatureIdsByKey = new Map();

  (Array.isArray(features) ? features : []).forEach((feature) => {
    // A feature the Result does not draw is found by its `search_id` (R-5).
    const svgId = String(feature?.search_id || feature?.svg_id || '').trim();
    if (!svgId || byId.has(svgId)) return;
    const qualifiers = getFeatureQualifiers(feature);
    /** @type {Map<string, PreparedSearchItem[]>} */
    const qualifierValuesByKey = new Map();
    /** @type {PreparedSearchItem[]} */
    const qualifierKeys = [];
    Object.keys(qualifiers).sort().forEach((key) => {
      const normalizedKey = normalizeSearchText(key).trim();
      if (!normalizedKey) return;
      qualifierKeys.push(prepareSearchItem({ label: t('Qualifier key'), value: key }));
      /** @type {SearchItem[]} */
      const values = [];
      appendSearchItems(values, t('Qualifier {key}', { key }), qualifiers[key]);
      const preparedValues = values.map(prepareSearchItem);
      qualifierValuesByKey.set(normalizedKey, preparedValues);
      if (!qualifierFeatureIdsByKey.has(normalizedKey)) qualifierFeatureIdsByKey.set(normalizedKey, []);
      // Set just above when it was missing.
      /** @type {string[]} */ (qualifierFeatureIdsByKey.get(normalizedKey)).push(svgId);
    });

    /** @type {Map<string, PreparedSearchItem[]>} */
    const itemsByField = new Map();
    ['label', 'type', 'record-id', 'location', 'strand', 'orthogroup', 'nucleotide', 'amino-acid']
      .forEach((field) => {
        if (normalizedPopupMode === 'simple' && isRichFeatureSearchField(field)) return;
        itemsByField.set(field, preparedItemsForField(feature, field, normalizedPopupMode, orthogroupsById));
      });
    itemsByField.set('qualifier-key', qualifierKeys);
    itemsByField.set('qualifier-value', Array.from(qualifierValuesByKey.values()).flat());
    itemsByField.set('all', preparedItemsForField(feature, 'all', normalizedPopupMode, orthogroupsById));

    featureOrder.push(svgId);
    byId.set(svgId, { feature, itemsByField, qualifierValuesByKey });
  });

  return {
    popupMode: normalizedPopupMode,
    featureOrder,
    byId,
    qualifierFeatureIdsByKey,
    orthogroupsById
  };
};

/**
 * @param {FeatureSearchDocument | undefined} document
 * @param {FeatureSearchMatcher} matcher
 * @param {string} field
 * @param {unknown} qualifierKey
 */
const preparedFeatureSearchMatches = (document, matcher, field, qualifierKey) => {
  if (!document || !matcher?.active || matcher.error) return [];
  const normalizedQualifierKey = normalizeSearchText(qualifierKey).trim();
  const items = field === 'qualifier-value' && normalizedQualifierKey
    ? document.qualifierValuesByKey.get(normalizedQualifierKey) || []
    : document.itemsByField.get(field) || [];
  /** @type {SearchMatchDetail[]} */
  const details = [];
  items.forEach((item) => {
    const matchedText = matcher.matchPrepared
      ? matcher.matchPrepared(item)
      : matcher.match(item.value, item.alphabet);
    if (!matchedText) return;
    details.push({ label: item.label, value: item.value, match: matchedText });
  });
  return details;
};

/**
 * @typedef {object} RunFeatureSearchOptions
 * @property {FeatureLike[]} [features]
 * @property {Iterable<string> | null} [renderedFeatureIds]
 * @property {string} [query]
 * @property {string} [field]
 * @property {string} [qualifierKey]
 * @property {boolean} [useRegex]
 * @property {string} [popupMode]
 * @property {SearchOrthogroup[]} [orthogroups]
 * @property {FeatureSearchIndex | null} [searchIndex]
 * @property {string} [previousActiveId]
 */

/**
 * @param {RunFeatureSearchOptions} [options]
 */
export const runFeatureSearch = ({
  features,
  renderedFeatureIds,
  query,
  field,
  qualifierKey,
  useRegex,
  popupMode = 'rich',
  orthogroups = [],
  searchIndex = null,
  previousActiveId = ''
} = {}) => {
  const normalizedPopupMode = normalizeFeatureSearchPopupMode(popupMode);
  const normalizedField = normalizeFeatureSearchField(field, { popupMode: normalizedPopupMode });
  const matcher = compileFeatureSearchMatcher(query, useRegex);
  const renderedIds = renderedFeatureIds instanceof Set
    ? renderedFeatureIds
    : new Set(Array.from(renderedFeatureIds || []).map((id) => String(id || '').trim()).filter(Boolean));
  const renderedFeatureCount = renderedIds.size;

  if (!matcher.active || matcher.error) {
    return {
      field: normalizedField,
      error: matcher.error || '',
      matches: [],
      matchDetails: {},
      activeIndex: -1,
      renderedFeatureCount
    };
  }

  const preparedIndex = searchIndex?.byId instanceof Map
    ? searchIndex
    : buildFeatureSearchIndex({ features, popupMode: normalizedPopupMode, orthogroups });
  /** @type {string[]} */
  const matches = [];
  /** @type {Record<string, SearchMatchDetail[]>} */
  const matchDetails = {};
  let candidateIds = preparedIndex.featureOrder;
  const normalizedQualifierKey = normalizeSearchText(qualifierKey).trim();
  if (normalizedField === 'qualifier-value' && normalizedQualifierKey) {
    candidateIds = preparedIndex.qualifierFeatureIdsByKey.get(normalizedQualifierKey) || [];
  }
  candidateIds.forEach((svgId) => {
    if (!svgId || !renderedIds.has(svgId)) return;
    const details = preparedFeatureSearchMatches(
      preparedIndex.byId.get(svgId), matcher, normalizedField, qualifierKey
    );
    if (!details.length) return;
    matches.push(svgId);
    matchDetails[svgId] = details;
  });

  let activeIndex = -1;
  if (matches.length) {
    const previousIndex = previousActiveId ? matches.indexOf(previousActiveId) : -1;
    activeIndex = previousIndex === -1 ? 0 : previousIndex;
  }

  return {
    field: normalizedField,
    error: '',
    matches,
    matchDetails,
    activeIndex,
    renderedFeatureCount
  };
};

/** @param {unknown} value */
const collapseWhitespace = (value) => String(value == null ? '' : value).replace(/\s+/g, ' ').trim();

/** @param {unknown} value @param {unknown} matchText */
const searchMatchSnippet = (value, matchText) => {
  const text = collapseWhitespace(value);
  const match = collapseWhitespace(matchText);
  if (!text) return '';
  if (!match) return text.length > 80 ? `${text.slice(0, 77)}...` : text;
  const index = text.toLowerCase().indexOf(match.toLowerCase());
  if (index === -1) return text.length > 80 ? `${text.slice(0, 77)}...` : text;
  const start = Math.max(0, index - 24);
  const end = Math.min(text.length, index + match.length + 24);
  return `${start > 0 ? '...' : ''}${text.slice(start, end)}${end < text.length ? '...' : ''}`;
};

/** @param {SearchMatchDetail | null | undefined} detail */
const formatSearchMatchDetail = (detail) => {
  if (!detail) return '';
  const label = collapseWhitespace(detail.label);
  const snippet = searchMatchSnippet(detail.value, detail.match);
  if (!label && !snippet) return '';
  if (!snippet) return label;
  return label ? `${label}: ${snippet}` : snippet;
};

// The line under the search row. The record comes first so that a truncated
// line still says which record the active match is on (OV-211).
/**
 * @param {FeatureLike | null | undefined} feature
 * @param {SearchMatchDetail[] | null | undefined} details
 * @returns {string}
 */
export const formatSearchMatchLine = (feature, details) => {
  const detailText = details?.length ? formatSearchMatchDetail(details[0]) : '';
  if (!detailText) return '';
  const recordId = String(feature?.record_id || '');
  return `${recordId ? `${t('Record {id}', { id: recordId })} · ` : ''}${t('Matched {detail}', { detail: detailText })}`;
};
