import { resolveColorToHex } from './color-utils.js';
import {
  FEATURE_INSTANCE_HASH_PATTERN,
  FEATURE_INSTANCE_HASH_QUALIFIER
} from '../services/feature-instance-identity.js';
import {
  FEATURE_SEMANTIC_SCOPE_QUALIFIER,
  isFeatureSemanticSelector
} from './feature-editor/semantic-fill-selectors.js';

export const INSTANCE_HASH_QUALIFIER = FEATURE_INSTANCE_HASH_QUALIFIER;
export const INSTANCE_HASH_PATTERN = FEATURE_INSTANCE_HASH_PATTERN;
export const CURRENT_SPECIFIC_RULE_SCHEMA = 6;

const normalizeRuleColor = (value) => {
  const raw = String(value ?? '').trim();
  if (raw.toLowerCase() === 'none') return 'none';
  return String(resolveColorToHex(raw) || '').toLowerCase();
};

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

export const parseSpecificRules = (
  text,
  { schema = CURRENT_SPECIFIC_RULE_SCHEMA } = {}
) => {
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
    const schema6ReservedLiteral = Number(schema) >= 6 && (
      qual === INSTANCE_HASH_QUALIFIER || qual === FEATURE_SEMANTIC_SCOPE_QUALIFIER
    );
    const instanceLiteral = schema6ReservedLiteral && qual === INSTANCE_HASH_QUALIFIER;
    const semanticLiteral = schema6ReservedLiteral && qual === FEATURE_SEMANTIC_SCOPE_QUALIFIER;
    if (instanceLiteral) {
      if (!INSTANCE_HASH_PATTERN.test(val)) {
        throw new Error(`Invalid feature-instance hash at line ${lineNo}: ${val}`);
      }
    } else if (semanticLiteral) {
      if (!isFeatureSemanticSelector(val)) {
        throw new Error(`Invalid Feature semantic selector at line ${lineNo}: ${val}`);
      }
    } else {
      try {
        new RegExp(val);
      } catch (error) {
        throw new Error(`Invalid specific-color regex at line ${lineNo}: ${error.message}`);
      }
    }
    const color = normalizeRuleColor(colorRaw);
    if (
      color !== 'none' &&
      !/^#(?:[0-9a-f]{3}|[0-9a-f]{4}|[0-9a-f]{6}|[0-9a-f]{8})$/.test(color)
    ) {
      throw new Error(`Invalid specific-color value at line ${lineNo}: ${colorRaw}`);
    }

    const rule = {
      feat,
      qual,
      val,
      color,
      cap: captionRaw,
      ...([INSTANCE_HASH_QUALIFIER, FEATURE_SEMANTIC_SCOPE_QUALIFIER].includes(qual)
        ? { match: schema6ReservedLiteral ? 'literal' : 'regex' }
        : {}),
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

export const serializeSpecificRules = (
  rules,
  { schema = CURRENT_SPECIFIC_RULE_SCHEMA } = {}
) => {
  const rows = (Array.isArray(rules) ? rules : [])
    .map((rule) => {
      const fields = [
        normalizeTsvCell(rule?.feat),
        normalizeTsvCell(rule?.qual),
        normalizeTsvCell(rule?.val),
        normalizeTsvCell(rule?.color),
        normalizeTsvCell(rule?.cap)
      ];
      if (fields[1] === INSTANCE_HASH_QUALIFIER) {
        if (Number(schema) < 6) {
          throw new Error('Feature-instance selectors require canonical request schema 6.');
        }
        if (!INSTANCE_HASH_PATTERN.test(fields[2])) {
          throw new Error(`Invalid feature-instance hash: ${fields[2]}`);
        }
      }
      if (fields[1] === FEATURE_SEMANTIC_SCOPE_QUALIFIER) {
        if (Number(schema) < 6) {
          throw new Error('Feature semantic selectors require canonical request schema 6.');
        }
        if (!isFeatureSemanticSelector(fields[2])) {
          throw new Error(`Invalid Feature semantic selector: ${fields[2]}`);
        }
      }
      return fields;
    })
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
