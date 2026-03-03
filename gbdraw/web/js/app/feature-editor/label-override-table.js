const LABEL_OVERRIDE_COLUMN_COUNT = 5;
const PRIMARY_HEADER = ['record_id', 'feature_type', 'qualifier', 'value', 'label_text'];
const LEGACY_HEADER = ['record', 'feature_type', 'qualifier_key', 'qualifier_value_regex', 'label_text'];
const CANDIDATE_QUALIFIERS = ['protein_id', 'locus_tag', 'gene'];
const DEFAULT_LABEL_QUALIFIER_PRIORITY = ['product', 'gene', 'locus_tag', 'protein_id', 'old_locus_tag', 'note'];

const escapeRegexLiteral = (value) => String(value ?? '').replace(/[.*+?^${}()|[\]\\]/g, '\\$&');
const normalizeTsvCell = (value) => String(value ?? '').replace(/[\t\r\n]+/g, ' ').trim();
const toSortedKeys = (obj) =>
  Object.keys(obj || {}).sort((a, b) => String(a || '').localeCompare(String(b || '')));
const normalizeFeatureIdKey = (value) => String(value ?? '').trim().toLowerCase();

const normalizeStrandToken = (value) => {
  const token = String(value ?? '').trim().toLowerCase();
  if (token === '+' || token === 'positive' || token === 'forward' || token === '1') return '+';
  if (token === '-' || token === 'negative' || token === 'reverse' || token === '-1') return '-';
  return 'undefined';
};

const normalizeQualifierMap = (qualifiers) => {
  const normalized = {};
  if (!qualifiers || typeof qualifiers !== 'object') return normalized;
  Object.entries(qualifiers).forEach(([keyRaw, valuesRaw]) => {
    const key = String(keyRaw || '').trim().toLowerCase();
    if (!key) return;
    const values = Array.isArray(valuesRaw) ? valuesRaw : [valuesRaw];
    normalized[key] = values
      .filter((value) => value !== null && value !== undefined)
      .map((value) => String(value));
  });
  return normalized;
};

const isHeaderRow = (parts) => {
  if (!Array.isArray(parts) || parts.length !== LABEL_OVERRIDE_COLUMN_COUNT) return false;
  const normalized = parts.map((part) => String(part || '').trim().toLowerCase());
  const primary = PRIMARY_HEADER.every((value, idx) => normalized[idx] === value);
  if (primary) return true;
  return LEGACY_HEADER.every((value, idx) => normalized[idx] === value);
};

export const buildFeatureMetadataMap = (features) => {
  const metadataByFeatureId = new Map();
  if (!Array.isArray(features)) return metadataByFeatureId;

  features.forEach((feature) => {
    const featureId = String(feature?.svg_id || '').trim();
    const key = normalizeFeatureIdKey(featureId);
    if (!key || metadataByFeatureId.has(key)) return;

    const start = String(feature?.start ?? '').trim();
    const end = String(feature?.end ?? '').trim();
    const location = start && end ? `${start}..${end}` : '';
    const position = location ? `${location}:${normalizeStrandToken(feature?.strand)}` : '';
    const record = String(feature?.record_id || '').trim();
    const featureType = String(feature?.type || '').trim();

    metadataByFeatureId.set(key, {
      featureId,
      record,
      location,
      position,
      featureType,
      qualifiers: normalizeQualifierMap(feature?.qualifiers)
    });
  });

  return metadataByFeatureId;
};

const buildEditableLabelByFeatureId = (editableLabels) => {
  const labelsByFeatureId = new Map();
  if (!Array.isArray(editableLabels)) return labelsByFeatureId;

  editableLabels.forEach((entry) => {
    const featureIdKey = normalizeFeatureIdKey(entry?.featureId);
    if (!featureIdKey || labelsByFeatureId.has(featureIdKey)) return;
    labelsByFeatureId.set(featureIdKey, {
      text: String(entry?.text ?? ''),
      sourceText: String(entry?.sourceText ?? '')
    });
  });
  return labelsByFeatureId;
};

const normalizeVisibilityOverrides = (visibilityOverrides) => {
  const normalized = new Map();
  if (!visibilityOverrides || typeof visibilityOverrides !== 'object') return normalized;

  Object.entries(visibilityOverrides).forEach(([featureIdRaw, modeRaw]) => {
    const featureIdKey = normalizeFeatureIdKey(featureIdRaw);
    if (!featureIdKey) return;
    const mode = String(modeRaw || '').trim().toLowerCase();
    if (mode !== 'on' && mode !== 'off') return;
    normalized.set(featureIdKey, mode);
  });

  return normalized;
};

const getRegexForFeatureHash = (featureIdRaw) => `^${escapeRegexLiteral(String(featureIdRaw || '').trim())}$`;

const resolveDefaultLabelText = (metadata, editableLabelEntry = null) => {
  const editableText = String(editableLabelEntry?.text || '').trim();
  if (editableText) return editableText;
  const sourceText = String(editableLabelEntry?.sourceText || '').trim();
  if (sourceText) return sourceText;

  if (!metadata) return '';
  for (const qualifier of DEFAULT_LABEL_QUALIFIER_PRIORITY) {
    const values = Array.isArray(metadata?.qualifiers?.[qualifier]) ? metadata.qualifiers[qualifier] : [];
    const first = values.find((value) => String(value || '').trim() !== '');
    if (first) return String(first);
  }

  const featureType = String(metadata.featureType || '').trim();
  const location = String(metadata.location || '').trim();
  if (featureType && location) return `${featureType} ${location}`;
  if (featureType) return featureType;
  if (location) return location;
  return '';
};

const buildSourceTextByFeatureId = (editableLabels, featureOverrideSources = {}) => {
  const sourceByFeatureId = new Map();

  if (featureOverrideSources && typeof featureOverrideSources === 'object') {
    Object.entries(featureOverrideSources).forEach(([featureIdRaw, sourceTextRaw]) => {
      const key = normalizeFeatureIdKey(featureIdRaw);
      const sourceText = String(sourceTextRaw ?? '').trim();
      if (!key || !sourceText || sourceByFeatureId.has(key)) return;
      sourceByFeatureId.set(key, sourceText);
    });
  }

  if (!Array.isArray(editableLabels)) return sourceByFeatureId;

  editableLabels.forEach((entry) => {
    const featureId = String(entry?.featureId || '').trim();
    const sourceText = String(entry?.sourceText || '').trim();
    const key = normalizeFeatureIdKey(featureId);
    if (!key || sourceByFeatureId.has(key)) return;
    sourceByFeatureId.set(key, sourceText);
  });

  return sourceByFeatureId;
};

const buildFeatureIdsBySourceText = (editableLabels) => {
  const featureIdsBySourceText = new Map();
  if (!Array.isArray(editableLabels)) return featureIdsBySourceText;

  editableLabels.forEach((entry) => {
    const sourceText = String(entry?.sourceText || '');
    const featureIdKey = normalizeFeatureIdKey(entry?.featureId);
    if (!sourceText || !featureIdKey) return;
    if (!featureIdsBySourceText.has(sourceText)) {
      featureIdsBySourceText.set(sourceText, new Set());
    }
    featureIdsBySourceText.get(sourceText).add(featureIdKey);
  });

  return featureIdsBySourceText;
};

const buildQualifierUniquenessIndex = (metadataByFeatureId) => {
  const counts = new Map();
  metadataByFeatureId.forEach((metadata) => {
    const record = String(metadata?.record || '').trim();
    const featureType = String(metadata?.featureType || '').trim();
    if (!record || !featureType) return;

    CANDIDATE_QUALIFIERS.forEach((qualifier) => {
      const values = Array.isArray(metadata?.qualifiers?.[qualifier]) ? metadata.qualifiers[qualifier] : [];
      values.forEach((valueRaw) => {
        const value = String(valueRaw || '').trim();
        if (!value) return;
        const key = `${record}\u0000${featureType}\u0000${qualifier}\u0000${value}`;
        counts.set(key, (counts.get(key) || 0) + 1);
      });
    });
  });
  return counts;
};

const selectPreferredQualifier = (metadata, qualifierUniquenessIndex) => {
  if (!metadata) return null;
  const record = String(metadata?.record || '').trim();
  const featureType = String(metadata?.featureType || '').trim();
  if (!record || !featureType) return null;

  for (const qualifier of CANDIDATE_QUALIFIERS) {
    const values = Array.isArray(metadata?.qualifiers?.[qualifier]) ? metadata.qualifiers[qualifier] : [];
    for (const valueRaw of values) {
      const value = String(valueRaw || '').trim();
      if (!value) continue;
      const key = `${record}\u0000${featureType}\u0000${qualifier}\u0000${value}`;
      if ((qualifierUniquenessIndex.get(key) || 0) === 1) {
        return { qualifier, value };
      }
    }
  }
  return null;
};

export const buildLabelOverrideRows = (featureOverrides, bulkOverrides, options = {}) => {
  const rows = [];
  let skippedFeatureCount = 0;
  let skippedFeatureSourceCount = 0;
  let skippedMissingSourceCount = 0;
  const featureMetadataById = buildFeatureMetadataMap(options.extractedFeatures);
  const editableLabelByFeatureId = buildEditableLabelByFeatureId(options.editableLabels);
  const sourceTextByFeatureId = buildSourceTextByFeatureId(
    options.editableLabels,
    options.featureOverrideSources
  );
  const featureIdsBySourceText = buildFeatureIdsBySourceText(options.editableLabels);
  const qualifierUniquenessIndex = buildQualifierUniquenessIndex(featureMetadataById);
  const visibilityOverridesByFeatureId = normalizeVisibilityOverrides(options.visibilityOverrides);
  const featureOverrideKeyById = new Map();
  toSortedKeys(featureOverrides).forEach((featureIdRaw) => {
    const key = normalizeFeatureIdKey(featureIdRaw);
    if (!key || featureOverrideKeyById.has(key)) return;
    featureOverrideKeyById.set(key, featureIdRaw);
  });
  const consumedFeatureOverrideKeys = new Set();

  Array.from(visibilityOverridesByFeatureId.keys())
    .sort((a, b) => a.localeCompare(b))
    .forEach((featureIdKey) => {
      const visibilityMode = visibilityOverridesByFeatureId.get(featureIdKey);
      const metadata = featureMetadataById.get(featureIdKey);
      const recordId = normalizeTsvCell(metadata?.record || '*') || '*';
      const featureType = normalizeTsvCell(metadata?.featureType || '*') || '*';
      const featureIdRaw = featureOverrideKeyById.get(featureIdKey) || featureIdKey;
      const featureOverrideKey = featureOverrideKeyById.get(featureIdKey);
      if (featureOverrideKey) {
        consumedFeatureOverrideKeys.add(featureOverrideKey);
      }
      let nextText = '';
      if (visibilityMode === 'on') {
        const overrideText = featureOverrideKey
          ? String(featureOverrides[featureOverrideKey] ?? '')
          : '';
        const fallbackText = resolveDefaultLabelText(
          metadata,
          editableLabelByFeatureId.get(featureIdKey) || null
        );
        nextText = normalizeTsvCell(overrideText || fallbackText || featureIdRaw || featureIdKey);
      }
      rows.push(
        `${recordId}\t${featureType}\thash\t${getRegexForFeatureHash(featureIdRaw)}\t${nextText}`
      );
    });

  toSortedKeys(featureOverrides).forEach((featureIdRaw) => {
    if (consumedFeatureOverrideKeys.has(featureIdRaw)) return;
    const featureId = String(featureIdRaw || '').trim();
    const featureIdKey = normalizeFeatureIdKey(featureId);
    if (!featureId || !featureIdKey) {
      skippedFeatureCount += 1;
      return;
    }

    if (visibilityOverridesByFeatureId.has(featureIdKey)) return;
    const metadata = featureMetadataById.get(featureIdKey);
    const recordId = normalizeTsvCell(metadata?.record || '*') || '*';
    const featureType = normalizeTsvCell(metadata?.featureType || '*') || '*';
    const sourceText = sourceTextByFeatureId.get(featureIdKey) || '';
    const preferredQualifier = selectPreferredQualifier(metadata, qualifierUniquenessIndex);

    let qualifier = 'hash';
    let qualifierRegex = getRegexForFeatureHash(featureId);
    if (preferredQualifier?.qualifier && preferredQualifier?.value) {
      qualifier = preferredQualifier.qualifier;
      qualifierRegex = `^${escapeRegexLiteral(preferredQualifier.value)}$`;
    } else if (sourceText) {
      qualifier = 'label';
      qualifierRegex = `^${escapeRegexLiteral(sourceText)}$`;
    }
    const nextText = normalizeTsvCell(featureOverrides[featureIdRaw]);

    rows.push(`${recordId}\t${featureType}\t${qualifier}\t${qualifierRegex}\t${nextText}`);
  });

  toSortedKeys(bulkOverrides).forEach((sourceTextRaw) => {
    const sourceText = String(sourceTextRaw ?? '');
    if (!sourceText) return;
    const nextText = normalizeTsvCell(bulkOverrides[sourceTextRaw]);
    const featureIdKeys = Array.from(featureIdsBySourceText.get(sourceText) || []).sort((a, b) =>
      String(a || '').localeCompare(String(b || ''))
    );
    const expandedRows = [];
    const seen = new Set();

    featureIdKeys.forEach((featureIdKey) => {
      const metadata = featureMetadataById.get(featureIdKey);
      if (!metadata) return;
      const recordId = normalizeTsvCell(metadata?.record || '*') || '*';
      const featureType = normalizeTsvCell(metadata?.featureType || '*') || '*';
      const preferredQualifier = selectPreferredQualifier(metadata, qualifierUniquenessIndex);
      const qualifier = preferredQualifier?.qualifier || 'label';
      const qualifierValue = preferredQualifier?.value || sourceText;
      const qualifierRegex = `^${escapeRegexLiteral(qualifierValue)}$`;
      const row = `${recordId}\t${featureType}\t${qualifier}\t${qualifierRegex}\t${nextText}`;
      if (seen.has(row)) return;
      seen.add(row);
      expandedRows.push(row);
    });

    if (expandedRows.length > 0) {
      rows.push(...expandedRows);
      return;
    }

    rows.push(`*\t*\tlabel\t^${escapeRegexLiteral(sourceText)}$\t${nextText}`);
  });

  return { rows, skippedFeatureCount, skippedFeatureSourceCount, skippedMissingSourceCount };
};

export const buildLabelOverrideTsv = (featureOverrides, bulkOverrides, options = {}) => {
  const { rows, skippedFeatureCount, skippedFeatureSourceCount, skippedMissingSourceCount } = buildLabelOverrideRows(
    featureOverrides,
    bulkOverrides,
    options
  );

  return {
    tsv: rows.length > 0 ? `${rows.join('\n')}\n` : '',
    rows,
    skippedFeatureCount,
    skippedFeatureSourceCount,
    skippedMissingSourceCount
  };
};

export const parseLabelOverrideTsv = (text) => {
  const rows = [];
  const sourceText = String(text ?? '');
  const lines = sourceText.split(/\r?\n/);

  lines.forEach((lineRaw, idx) => {
    const lineNo = idx + 1;
    const trimmed = lineRaw.trim();
    if (!trimmed || trimmed.startsWith('#')) return;

    const parts = lineRaw.split('\t');
    if (parts.length !== LABEL_OVERRIDE_COLUMN_COUNT) {
      throw new Error(
        `Invalid label TSV at line ${lineNo}: expected ${LABEL_OVERRIDE_COLUMN_COUNT} columns, found ${parts.length}.`
      );
    }

    if (isHeaderRow(parts)) return;

    const [recordIdRaw, featureTypeRaw, qualifierRaw, valueRegexRaw, labelTextRaw] = parts;
    const recordId = String(recordIdRaw ?? '').trim();
    const featureType = String(featureTypeRaw ?? '').trim();
    const qualifier = String(qualifierRaw ?? '').trim();
    const valueRegex = String(valueRegexRaw ?? '').trim();
    const labelText = String(labelTextRaw ?? '').replace(/[\r\n]+/g, ' ');

    if (!recordId) {
      throw new Error(`Invalid label TSV at line ${lineNo}: column 1 (record_id) is required.`);
    }
    if (!featureType) {
      throw new Error(`Invalid label TSV at line ${lineNo}: column 2 (feature_type) is required.`);
    }
    if (!qualifier) {
      throw new Error(`Invalid label TSV at line ${lineNo}: column 3 (qualifier) is required.`);
    }
    if (!valueRegex) {
      throw new Error(`Invalid label TSV at line ${lineNo}: column 4 (value) is required.`);
    }

    let qualifierValuePattern;
    try {
      qualifierValuePattern = new RegExp(valueRegex, 'i');
    } catch (error) {
      throw new Error(
        `Invalid label TSV at line ${lineNo}: invalid regex '${valueRegex}' (${error.message}).`
      );
    }

    rows.push({
      lineNo,
      recordId,
      featureType,
      qualifier,
      valueRegex,
      qualifierValuePattern,
      labelText,
      isGlobalLabelRule: recordId === '*' && featureType === '*' && qualifier.toLowerCase() === 'label'
    });
  });

  return rows;
};
