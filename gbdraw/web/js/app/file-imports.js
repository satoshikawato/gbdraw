import { normalizeSpecificRuleColor, resolveColorToHex } from './color-utils.js';
import { diagnosticError } from '../services/error-normalization.js';

const SPECIFIC_RULE_COLUMNS = Object.freeze(['feature_type', 'qualifier', 'pattern', 'color']);

// Python reads these tables with a fixed column count and rejects any other row
// (gbdraw/io/table_text.py::read_literal_table), so the import does too (R4). The
// diagnostic is the one the Label override parser uses.
const requireColumnCount = (parts, columnCount, row) => {
  if (parts.length !== columnCount) {
    throw diagnosticError('TABLE_INVALID', { row, columnCount, reason: 'FIELDS' });
  }
};

export const parseColorTable = (text) => {
  const colors = {};
  let count = 0;
  const lines = text.split(/\r?\n/);

  for (const [index, line] of lines.entries()) {
    if (!line.trim() || line.trim().startsWith('#') || line.trim().startsWith('[')) continue;
    const parts = line.split('\t');
    requireColumnCount(parts, 2, index + 1);
    const key = parts[0].trim();
    const color = parts[1].trim();
    if (key.toLowerCase() === 'feature_type' && color.toLowerCase() === 'color') continue;
    if (key && (color.startsWith('#') || /^[a-z]+$/i.test(color))) {
      colors[key] = resolveColorToHex(color);
      count++;
    }
  }

  return { colors, count };
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

const normalizeTsvCell = (value) => String(value ?? '').replace(/[\t\r\n]+/g, ' ').trim();

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

export const parsePriorityRules = (text) => {
  const rules = [];
  const lines = text.split(/\r?\n/);

  for (const [index, line] of lines.entries()) {
    if (!line.trim() || line.trim().startsWith('#')) continue;
    const parts = line.split('\t');
    requireColumnCount(parts, 2, index + 1);
    if (
      parts[0].trim().toLowerCase() === 'feature_type' &&
      parts[1].trim().toLowerCase() === 'priorities'
    ) continue;
    rules.push({ feat: parts[0].trim(), order: parts[1].trim() });
  }

  return { rules, count: rules.length };
};

export const parseWhitelistRules = (text) => {
  const rules = [];
  const lines = text.split(/\r?\n/);

  for (const [index, line] of lines.entries()) {
    if (!line.trim() || line.trim().startsWith('#')) continue;
    const parts = line.split('\t');
    requireColumnCount(parts, 3, index + 1);
    rules.push({
      feat: parts[0].trim(),
      qual: parts[1].trim(),
      key: parts[2].trim()
    });
  }

  return { rules, count: rules.length };
};

export const parseBlacklistWords = (text) => {
  const words = text
    .split(/[\r\n,]+/)
    .map((word) => word.trim())
    .filter((word) => word);

  return { words, count: words.length };
};
