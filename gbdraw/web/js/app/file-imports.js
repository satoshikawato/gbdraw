import { resolveColorToHex } from './color-utils.js';

export const parseColorTable = (text) => {
  const colors = {};
  let count = 0;
  const lines = text.split(/\r?\n/);

  for (const line of lines) {
    if (!line.trim() || line.trim().startsWith('#') || line.trim().startsWith('[')) continue;
    const parts = line.split('\t');
    if (parts.length < 2) continue;
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
      throw new Error(`Invalid specific-color TSV at line ${lineNo}: expected 4 or 5 columns.`);
    }

    const [feat, qual, val, colorRaw, captionRaw = ''] = parts.map((part) => part.trim());
    const required = [feat, qual, val, colorRaw];
    const missingIndex = required.findIndex((value) => !value);
    if (missingIndex >= 0) {
      throw new Error(
        `Invalid specific-color TSV at line ${lineNo}: column ${missingIndex + 1} is required.`
      );
    }
    try {
      new RegExp(val);
    } catch (error) {
      throw new Error(`Invalid specific-color regex at line ${lineNo}: ${error.message}`);
    }
    const color = String(resolveColorToHex(colorRaw) || '').toLowerCase();
    if (!/^#(?:[0-9a-f]{3}|[0-9a-f]{4}|[0-9a-f]{6}|[0-9a-f]{8})$/.test(color)) {
      throw new Error(`Invalid specific-color value at line ${lineNo}: ${colorRaw}`);
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

  for (const line of lines) {
    if (!line.trim() || line.trim().startsWith('#')) continue;
    const parts = line.split('\t');
    if (parts.length < 2) continue;
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

  for (const line of lines) {
    if (!line.trim() || line.trim().startsWith('#')) continue;
    const parts = line.split('\t');
    if (parts.length < 3) continue;
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
