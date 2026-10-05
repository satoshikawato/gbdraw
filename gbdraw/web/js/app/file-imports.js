import { normalizeSpecificRuleColor, resolveColorToHex } from './color-utils.js';
import { diagnosticError } from '../services/error-normalization.js';
import { normalizeTsvCell } from '../utils/tsv-cell.js';

const SPECIFIC_RULE_COLUMNS = Object.freeze(['feature_type', 'qualifier', 'pattern', 'color']);

// Python reads these tables with a fixed column count and rejects any other row
// (gbdraw/io/table_text.py::read_literal_table), so the import does too (R4). The
// diagnostic is the one the Label override parser uses.
const requireColumnCount = (parts, columnCount, row) => {
  if (parts.length !== columnCount) {
    throw diagnosticError('TABLE_INVALID', { row, columnCount, reason: 'FIELDS' });
  }
};

// Sessions 31–39 stored these tables as their Web writer wrote them, before cell
// values were normalized (normalizeTsvCell): a tab in a value made extra cells and
// a line break a short row. For them (`repairs`), a row is read as the current
// writer writes it: the extra cells join the last field with one space, and a row
// without its required fields (`complete`) is dropped; `repairs` lists each one.
const tableCells = (line, columnCount, row, repairs, complete) => {
  const parts = line.split('\t');
  if (!repairs) {
    requireColumnCount(parts, columnCount, row);
    return parts.map((part) => part.trim());
  }
  const cells = [...parts.slice(0, columnCount - 1), parts.slice(columnCount - 1).join('\t')]
    .map(normalizeTsvCell);
  if (parts.length < columnCount || !complete(cells)) {
    repairs.push({ row, repair: 'dropped' });
    return null;
  }
  if (parts.length > columnCount) repairs.push({ row, repair: 'joined' });
  return cells;
};
const isTableColor = (color) => color.startsWith('#') || /^[a-z]+$/i.test(color);

export const parseColorTable = (text, { legacyRows = false } = {}) => {
  const colors = {};
  let count = 0;
  const repairs = legacyRows ? [] : null;
  const lines = text.split(/\r?\n/);

  for (const [index, line] of lines.entries()) {
    if (!line.trim() || line.trim().startsWith('#') || line.trim().startsWith('[')) continue;
    const cells = tableCells(line, 2, index + 1, repairs, ([key, color]) => key && isTableColor(color));
    if (!cells) continue;
    const [key, color] = cells;
    if (key.toLowerCase() === 'feature_type' && color.toLowerCase() === 'color') continue;
    if (key && isTableColor(color)) {
      colors[key] = resolveColorToHex(color);
      count++;
    }
  }

  return { colors, count, repairs: repairs || [] };
};

export const parseSpecificRules = (text) => {
  const rules = [];
  const rulesWithCaptions = [];
  const identities = new Set();
  const lines = String(text ?? '').split(/\r?\n/);

  for (let index = 0; index < lines.length; index += 1) {
    const line = lines[index];
    const lineNo = index + 1;
    if (!line.trim() || line.trim().startsWith('#')) continue;
    const parts = line.split('\t');
    if (
      parts.length >= 2 &&
      parts[0].trim().toLowerCase() === 'feature_type' &&
      parts[1].trim().toLowerCase() === 'qualifier_key'
    ) continue;

    if (parts.length < 4 || parts.length > 5) {
      throw diagnosticError('TABLE_INVALID', { row: lineNo, reason: 'SPECIFIC_COLUMNS' });
    }

    const [feat, qual, val, colorRaw, captionRaw = ''] = parts.map((part) => part.trim());
    const required = [feat, qual, val, colorRaw];
    const missingIndex = required.findIndex((value) => !value);
    if (missingIndex >= 0) {
      throw diagnosticError('TABLE_INVALID', { row: lineNo, field: SPECIFIC_RULE_COLUMNS[missingIndex], reason: 'REQUIRED' });
    }
    const color = normalizeSpecificRuleColor(colorRaw);
    if (!color) {
      throw diagnosticError('TABLE_INVALID', { row: lineNo, field: 'color', reason: 'COLOR' });
    }

    const rule = {
      feat,
      qual,
      val,
      color,
      cap: captionRaw,
      fromFile: true
    };
    const identity = JSON.stringify([rule.feat, rule.qual, rule.val, rule.color, rule.cap]);
    if (identities.has(identity)) continue;
    identities.add(identity);
    rules.push(rule);
    if (rule.cap) rulesWithCaptions.push(rule);
  }

  return { rules, rulesWithCaptions, count: rules.length };
};

export const serializeSpecificRules = (rules) => {
  const rows = (Array.isArray(rules) ? rules : [])
    .map((rule) => [
      normalizeTsvCell(rule?.feat),
      normalizeTsvCell(rule?.qual),
      normalizeTsvCell(rule?.val),
      normalizeTsvCell(rule?.color),
      normalizeTsvCell(rule?.cap)
    ])
    .filter((fields) => fields.slice(0, 4).every(Boolean))
    .map((fields) => fields.join('\t'));

  return rows.length > 0 ? `${rows.join('\n')}\n` : '';
};

export const serializeLabelWhitelistRules = (rules) => {
  const rows = (Array.isArray(rules) ? rules : [])
    .map((rule) => [rule?.feat, rule?.qual, rule?.key].map(normalizeTsvCell))
    .filter(([feat, qual]) => feat && qual)
    .map((fields) => fields.join('\t'));

  return rows.length > 0 ? `${rows.join('\n')}\n` : '';
};

export const serializeQualifierPriorityRules = (rules) => {
  const rows = (Array.isArray(rules) ? rules : [])
    .map((rule) => [rule?.feat, rule?.order].map(normalizeTsvCell))
    .filter((fields) => fields.every(Boolean))
    .map((fields) => fields.join('\t'));

  return rows.length > 0 ? `${rows.join('\n')}\n` : '';
};

export const parsePriorityRules = (text, { legacyRows = false } = {}) => {
  const rules = [];
  const repairs = legacyRows ? [] : null;
  const lines = text.split(/\r?\n/);

  for (const [index, line] of lines.entries()) {
    if (!line.trim() || line.trim().startsWith('#')) continue;
    const cells = tableCells(line, 2, index + 1, repairs, ([feat, order]) => feat && order);
    if (!cells) continue;
    const [feat, order] = cells;
    if (feat.toLowerCase() === 'feature_type' && order.toLowerCase() === 'priorities') continue;
    rules.push({ feat, order });
  }

  return { rules, count: rules.length, repairs: repairs || [] };
};

export const parseWhitelistRules = (text, { legacyRows = false } = {}) => {
  const rules = [];
  const repairs = legacyRows ? [] : null;
  const lines = text.split(/\r?\n/);

  for (const [index, line] of lines.entries()) {
    if (!line.trim() || line.trim().startsWith('#')) continue;
    const cells = tableCells(line, 3, index + 1, repairs, ([feat, qual]) => feat && qual);
    if (!cells) continue;
    const [feat, qual, key] = cells;
    rules.push({ feat, qual, key });
  }

  return { rules, count: rules.length, repairs: repairs || [] };
};

export const parseBlacklistWords = (text) => {
  const words = text
    .split(/[\r\n,]+/)
    .map((word) => word.trim())
    .filter((word) => word);

  return { words, count: words.length };
};
