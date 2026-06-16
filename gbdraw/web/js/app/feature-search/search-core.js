export const RICH_FEATURE_SEARCH_FIELD_IDS = Object.freeze([
  'qualifier-key',
  'qualifier-value',
  'nucleotide',
  'amino-acid'
]);

export const FEATURE_SEARCH_FIELD_DEFINITIONS = Object.freeze([
  { id: 'all', label: 'All' },
  { id: 'label', label: 'Label' },
  { id: 'type', label: 'Feature type' },
  { id: 'record-id', label: 'Record ID' },
  { id: 'location', label: 'Location' },
  { id: 'strand', label: 'Strand' },
  { id: 'orthogroup', label: 'Orthogroup' },
  { id: 'qualifier-key', label: 'Qualifier key', richOnly: true },
  { id: 'qualifier-value', label: 'Qualifier value', richOnly: true },
  { id: 'nucleotide', label: 'Nucleotide', richOnly: true },
  { id: 'amino-acid', label: 'Amino acid', richOnly: true }
]);

const FIELD_IDS = new Set(FEATURE_SEARCH_FIELD_DEFINITIONS.map((field) => field.id));
const RICH_FIELD_IDS = new Set(RICH_FEATURE_SEARCH_FIELD_IDS);

export const isRichFeatureSearchField = (field) => RICH_FIELD_IDS.has(String(field || ''));

export const normalizeFeatureSearchPopupMode = (popupMode) => (
  popupMode === 'simple' ? 'simple' : 'rich'
);

export const normalizeFeatureSearchField = (field, { popupMode = 'rich' } = {}) => {
  const selected = String(field || 'all').trim() || 'all';
  if (!FIELD_IDS.has(selected)) return 'all';
  if (normalizeFeatureSearchPopupMode(popupMode) === 'simple' && isRichFeatureSearchField(selected)) {
    return 'all';
  }
  return selected;
};

export const getFeatureSearchFieldOptions = ({ popupMode = 'rich' } = {}) => {
  const normalizedPopupMode = normalizeFeatureSearchPopupMode(popupMode);
  return FEATURE_SEARCH_FIELD_DEFINITIONS.filter((field) => (
    normalizedPopupMode !== 'simple' || !field.richOnly
  ));
};

const normalizeSearchText = (value) => String(value == null ? '' : value).toLowerCase();

export const NUCLEOTIDE_IUPAC = Object.freeze({
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

export const AMINO_ACID_IUPAC = Object.freeze({
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

const getIupacAlphabet = (alphabet) => {
  if (alphabet === 'nucleotide') return NUCLEOTIDE_IUPAC;
  if (alphabet === 'amino-acid') return AMINO_ACID_IUPAC;
  return null;
};

const uniqueCharacters = (value) => {
  const seen = new Set();
  const chars = [];
  String(value || '').split('').forEach((char) => {
    if (!char || seen.has(char)) return;
    seen.add(char);
    chars.push(char);
  });
  return chars;
};

const charSetsIntersect = (left, right) => {
  const rightSet = new Set(uniqueCharacters(right));
  return uniqueCharacters(left).some((char) => rightSet.has(char));
};

const escapeRegExpText = (value) => String(value).replace(/[.*+?^${}()|[\]\\]/g, '\\$&');
const escapeRegExpClassChar = (value) => String(value).replace(/[\\\]\[\^-]/g, '\\$&');

export const buildIupacQueryPattern = (query, alphabet) => {
  const map = getIupacAlphabet(alphabet);
  const normalized = String(query || '').replace(/\s+/g, '').toUpperCase();
  if (!map || !normalized) return null;
  const targetChars = Object.keys(map);
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

const appendSearchValues = (values, value) => {
  if (Array.isArray(value)) {
    value.forEach((entry) => appendSearchValues(values, entry));
    return;
  }
  if (value === null || value === undefined) return;
  const text = String(value).trim();
  if (text) values.push(text);
};

const appendSearchItems = (items, label, value, options = {}) => {
  if (Array.isArray(value)) {
    value.forEach((entry) => appendSearchItems(items, label, entry, options));
    return;
  }
  if (value === null || value === undefined) return;
  const text = String(value).trim();
  if (!text) return;
  items.push({
    label,
    value: text,
    alphabet: options.alphabet ? String(options.alphabet) : ''
  });
};

export const getFeatureQualifiers = (feature) => (
  feature?.qualifiers && typeof feature.qualifiers === 'object' && !Array.isArray(feature.qualifiers)
    ? feature.qualifiers
    : {}
);

const getQualifierSearchItems = (feature, qualifierKey) => {
  const target = normalizeSearchText(qualifierKey).trim();
  const qualifiers = getFeatureQualifiers(feature);
  const items = [];
  Object.keys(qualifiers).sort().forEach((key) => {
    if (target && normalizeSearchText(key) !== target) return;
    appendSearchItems(items, `Qualifier ${key}`, qualifiers[key]);
  });
  return items;
};

const firstQualifierValue = (feature, key) => {
  const qualifiers = getFeatureQualifiers(feature);
  const values = [];
  appendSearchValues(values, qualifiers[key]);
  return values.find((value) => value.trim()) || '';
};

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

const getLabelSearchItems = (feature) => {
  const items = [];
  const seen = new Set();
  getFeatureLabelCandidates(feature).forEach((value) => {
    const values = [];
    appendSearchValues(values, value);
    values.forEach((entry) => {
      const key = normalizeSearchText(entry);
      if (!key || seen.has(key)) return;
      seen.add(key);
      appendSearchItems(items, key === normalizeSearchText(feature?.svg_id) ? 'SVG ID' : 'Label', entry);
    });
  });
  return items;
};

const buildFeatureLocation = (feature) => {
  const direct = String(feature?.location || '').trim();
  if (direct) return direct;
  const parts = Array.isArray(feature?.location_parts) ? feature.location_parts : [];
  const partText = parts
    .map((part) => String(part?.display || '').trim())
    .filter(Boolean)
    .join(', ');
  if (partText) return partText;
  const start = Number(feature?.start);
  const end = Number(feature?.end);
  const startText = Number.isFinite(start) ? String(start + 1) : String(feature?.start ?? '');
  const endText = Number.isFinite(end) ? String(end) : String(feature?.end ?? '');
  const range = `${startText}..${endText}`;
  const strand = String(feature?.strand || '').trim();
  return strand ? `${range} (${strand})` : range;
};

const buildOrthogroupMap = (orthogroups) => {
  const groups = new Map();
  (Array.isArray(orthogroups) ? orthogroups : []).forEach((group) => {
    const id = String(group?.id || group?.orthogroupId || group?.orthogroup_id || '').trim();
    if (id && !groups.has(id)) groups.set(id, group);
  });
  return groups;
};

const getOrthogroupId = (feature) => String(feature?.orthogroup_id || feature?.orthogroupId || '').trim();

const getFeatureOrthogroupMember = (feature, group) => {
  if (feature?.orthogroup_member && typeof feature.orthogroup_member === 'object') return feature.orthogroup_member;
  if (feature?.orthogroupMember && typeof feature.orthogroupMember === 'object') return feature.orthogroupMember;
  const members = Array.isArray(group?.members) ? group.members : [];
  const svgId = String(feature?.svg_id || '').trim();
  const recordIndex = Number(feature?.record_idx ?? feature?.recordIndex ?? feature?.fileIdx);
  return members.find((member) => {
    const memberSvgId = String(member?.featureSvgId || member?.feature_svg_id || '').trim();
    if (memberSvgId !== svgId) return false;
    if (!Number.isInteger(recordIndex)) return true;
    return Number(member?.recordIndex ?? member?.record_index) === recordIndex;
  }) || null;
};

const displayProteinId = (feature, member) => String(
  feature?.source_protein_id ||
  feature?.sourceProteinId ||
  member?.sourceProteinId ||
  member?.source_protein_id ||
  feature?.protein_id ||
  feature?.proteinId ||
  member?.proteinId ||
  member?.protein_id ||
  ''
).trim();

const internalProteinId = (feature, member) => String(
  feature?.protein_id ||
  feature?.proteinId ||
  member?.proteinId ||
  member?.protein_id ||
  ''
).trim();

const getOrthogroupSearchItems = (feature, orthogroupsById) => {
  const orthogroupId = getOrthogroupId(feature);
  const group = orthogroupsById.get(orthogroupId) || null;
  const member = getFeatureOrthogroupMember(feature, group);
  const proteinId = displayProteinId(feature, member);
  const internalId = internalProteinId(feature, member);
  const items = [];
  appendSearchItems(items, 'Orthogroup ID', orthogroupId);
  appendSearchItems(items, 'Orthogroup name', group?.display_name || group?.displayName || group?.name);
  appendSearchItems(items, 'Orthogroup description', group?.description);
  appendSearchItems(items, 'Protein ID', proteinId);
  if (internalId && internalId !== proteinId) {
    appendSearchItems(items, 'Internal protein ID', internalId);
  }
  appendSearchItems(items, 'Orthogroup member', member?.label);
  appendSearchItems(items, 'Orthogroup member gene', member?.gene);
  appendSearchItems(items, 'Orthogroup member product', member?.product);
  appendSearchItems(items, 'Orthogroup member note', member?.note);
  appendSearchItems(items, 'Orthogroup member protein ID', member?.sourceProteinId || member?.source_protein_id || member?.proteinId || member?.protein_id);
  return items;
};

export const featureSearchItems = (
  feature,
  field,
  qualifierKey,
  {
    popupMode = 'rich',
    orthogroupsById = new Map()
  } = {}
) => {
  const items = [];
  const selectedField = normalizeFeatureSearchField(field, { popupMode });
  const qualifiers = getFeatureQualifiers(feature);

  if (selectedField === 'label') return getLabelSearchItems(feature);
  if (selectedField === 'type') {
    appendSearchItems(items, 'Feature type', feature?.type);
    return items;
  }
  if (selectedField === 'record-id') {
    appendSearchItems(items, 'Record ID', feature?.record_id);
    appendSearchItems(items, 'Record ID', feature?.recordId);
    appendSearchItems(items, 'Record ID', feature?.displayRecordId);
    return items;
  }
  if (selectedField === 'location') {
    appendSearchItems(items, 'Location', buildFeatureLocation(feature));
    appendSearchItems(items, 'Start', feature?.start);
    appendSearchItems(items, 'End', feature?.end);
    return items;
  }
  if (selectedField === 'strand') {
    appendSearchItems(items, 'Strand', feature?.strand);
    return items;
  }
  if (selectedField === 'orthogroup') {
    return getOrthogroupSearchItems(feature, orthogroupsById);
  }
  if (selectedField === 'qualifier-key') {
    appendSearchItems(items, 'Qualifier key', Object.keys(qualifiers));
    return items;
  }
  if (selectedField === 'qualifier-value') {
    return getQualifierSearchItems(feature, qualifierKey);
  }
  if (selectedField === 'nucleotide') {
    appendSearchItems(items, 'Nucleotide sequence', feature?.nucleotide_sequence || feature?.nucleotideSequence, { alphabet: 'nucleotide' });
    return items;
  }
  if (selectedField === 'amino-acid') {
    appendSearchItems(items, 'Amino acid sequence', feature?.amino_acid_sequence || feature?.aminoAcidSequence, { alphabet: 'amino-acid' });
    return items;
  }

  items.push(...getLabelSearchItems(feature));
  appendSearchItems(items, 'Record ID', feature?.record_id);
  appendSearchItems(items, 'Record ID', feature?.recordId);
  appendSearchItems(items, 'Record ID', feature?.displayRecordId);
  appendSearchItems(items, 'Feature type', feature?.type);
  appendSearchItems(items, 'Location', buildFeatureLocation(feature));
  appendSearchItems(items, 'Strand', feature?.strand);
  items.push(...getOrthogroupSearchItems(feature, orthogroupsById));
  if (normalizeFeatureSearchPopupMode(popupMode) !== 'simple') {
    appendSearchItems(items, 'Qualifier key', Object.keys(qualifiers));
    Object.keys(qualifiers).sort().forEach((key) => {
      appendSearchItems(items, `Qualifier ${key}`, qualifiers[key]);
    });
    appendSearchItems(items, 'Nucleotide sequence', feature?.nucleotide_sequence || feature?.nucleotideSequence, { alphabet: 'nucleotide' });
    appendSearchItems(items, 'Amino acid sequence', feature?.amino_acid_sequence || feature?.aminoAcidSequence, { alphabet: 'amino-acid' });
  }
  return items;
};

export const compileFeatureSearchMatcher = (query, useRegex) => {
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
        test: (values) => values.some((value) => {
          regex.lastIndex = 0;
          return regex.test(String(value == null ? '' : value));
        })
      };
    } catch {
      return { active: true, error: 'Invalid regex', match: () => '', test: () => false };
    }
  }

  const needle = normalizeSearchText(trimmedQuery);
  const sequenceMatchers = {};
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
    test: (values) => values.some((value) => normalizeSearchText(value).includes(needle))
  };
};

export const featureSearchMatches = (feature, matcher, field, qualifierKey, options = {}) => {
  if (!matcher?.active || matcher.error) return [];
  return featureSearchItems(feature, field, qualifierKey, options)
    .map((item) => {
      const matchedText = matcher.match ? matcher.match(item.value, item.alphabet) : '';
      if (!matchedText) return null;
      return {
        label: item.label,
        value: String(item.value),
        match: matchedText
      };
    })
    .filter(Boolean);
};

export const runFeatureSearch = ({
  features,
  renderedFeatureIds,
  query,
  field,
  qualifierKey,
  useRegex,
  popupMode = 'rich',
  orthogroups = [],
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

  const orthogroupsById = buildOrthogroupMap(orthogroups);
  const matches = [];
  const matchDetails = {};
  (Array.isArray(features) ? features : []).forEach((feature) => {
    const svgId = String(feature?.svg_id || '').trim();
    if (!svgId || !renderedIds.has(svgId)) return;
    const details = featureSearchMatches(feature, matcher, normalizedField, qualifierKey, {
      popupMode: normalizedPopupMode,
      orthogroupsById
    });
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

const collapseWhitespace = (value) => String(value == null ? '' : value).replace(/\s+/g, ' ').trim();

export const searchMatchSnippet = (value, matchText) => {
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

export const formatSearchMatchDetail = (detail) => {
  if (!detail) return '';
  const label = collapseWhitespace(detail.label);
  const snippet = searchMatchSnippet(detail.value, detail.match);
  if (!label && !snippet) return '';
  if (!snippet) return label;
  return label ? `${label}: ${snippet}` : snippet;
};
