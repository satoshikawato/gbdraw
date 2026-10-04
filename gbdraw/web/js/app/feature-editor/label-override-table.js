import { escapeRegexLiteral } from '../feature-selector.js';
import { recordStructuralMetric } from '../../services/runtime-test-hooks.js';

const LABEL_OVERRIDE_COLUMN_COUNT = 5;
const PRIMARY_HEADER = ['record_id', 'feature_type', 'qualifier', 'value', 'label_text'];
const LEGACY_HEADER = ['record', 'feature_type', 'qualifier_key', 'qualifier_value_regex', 'label_text'];

const normalizeTsvCell = (value) => String(value ?? '').replace(/[\t\r\n]+/g, ' ').trim();
const toSortedKeys = (obj) =>
  Object.keys(obj || {}).sort((a, b) => String(a || '').localeCompare(String(b || '')));

const isHeaderRow = (parts) => {
  if (!Array.isArray(parts) || parts.length !== LABEL_OVERRIDE_COLUMN_COUNT) return false;
  const normalized = parts.map((part) => String(part || '').trim().toLowerCase());
  const primary = PRIMARY_HEADER.every((value, idx) => normalized[idx] === value);
  if (primary) return true;
  return LEGACY_HEADER.every((value, idx) => normalized[idx] === value);
};

/**
 * Projects the bulk label edits (`labelTextBulkOverrides`, {sourceText: text})
 * for one request. Per-feature edits are identity rows (`featureOverrides`,
 * design Q4); a bulk edit is a source-text rule. It reaches every feature whose
 * label shows that source text, known from `labelTargets` (the displayed
 * Result's labels, [{identityKey, sourceText}]) and the label source text the
 * draft rows recorded (any Result, B6), as those features' `labelText`; a
 * feature's own text wins. A source text that no known label shows, and a
 * blank text, stay a `* * label ^text$` row of the label table.
 */
export const buildBulkLabelProjection = (bulkOverrides, options = {}) => {
  const bulkLabelText = {};
  const rows = [];
  const sourceTexts = toSortedKeys(bulkOverrides);
  if (sourceTexts.length === 0) return { bulkLabelText, rows };
  const labelTargets = Array.isArray(options.labelTargets) ? options.labelTargets : [];
  const featureOverrides = options.featureOverrides || {};
  recordStructuralMetric('labelOverrideTableBuildCount', 1, { featureCount: labelTargets.length });
  const identitiesBySourceText = new Map();
  const add = (sourceText, identityKey) => {
    if (!sourceText || !identityKey) return;
    if (!identitiesBySourceText.has(sourceText)) identitiesBySourceText.set(sourceText, new Set());
    identitiesBySourceText.get(sourceText).add(identityKey);
  };
  labelTargets.forEach((entry) => add(String(entry?.sourceText ?? ''), entry?.identityKey));
  Object.entries(featureOverrides || {}).forEach(([identityKey, row]) => add(row?.labelSourceText, identityKey));
  sourceTexts.forEach((sourceTextRaw) => {
    const sourceText = String(sourceTextRaw ?? '');
    if (!sourceText) return;
    const text = String(bulkOverrides[sourceTextRaw] ?? '');
    const identities = identitiesBySourceText.get(sourceText);
    // A blank text hides the matching labels; only the table rule says that.
    if (identities?.size && text.trim()) {
      identities.forEach((identityKey) => { bulkLabelText[identityKey] = text; });
      return;
    }
    rows.push(`*\t*\tlabel\t^${escapeRegexLiteral(sourceText)}$\t${normalizeTsvCell(text)}`);
  });
  return { bulkLabelText, rows };
};

export const serializeLabelOverrideRows = (rows) => {
  const serialized = (Array.isArray(rows) ? rows : [])
    .map((row) => [
      normalizeTsvCell(row?.recordId ?? row?.record_id),
      normalizeTsvCell(row?.featureType ?? row?.feature_type),
      normalizeTsvCell(row?.qualifier),
      normalizeTsvCell(row?.valueRegex ?? row?.value),
      normalizeTsvCell(row?.labelText ?? row?.label_text)
    ])
    .filter((fields) => fields.slice(0, 4).every(Boolean))
    .map((fields) => fields.join('\t'));
  return serialized.length > 0 ? `${serialized.join('\n')}\n` : '';
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

    rows.push({
      lineNo,
      recordId,
      featureType,
      qualifier,
      valueRegex,
      labelText,
      isGlobalLabelRule: recordId === '*' && featureType === '*' && qualifier.toLowerCase() === 'label'
    });
  });

  return rows;
};
