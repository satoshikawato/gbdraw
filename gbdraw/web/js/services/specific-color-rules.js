// @ts-check
import { resolveColorToHex } from '../utils/color-utils.js';
import { parseSpecificRules } from './file-imports.js';

export const SPECIFIC_COLOR_FILE_OWNER = 'specific-color-file';

const normalizeText = (value) => String(value ?? '').trim();
const normalizeColor = (value) => String(resolveColorToHex(normalizeText(value)) || '').toLowerCase();

export const normalizeSpecificRule = (rule, { fromFile = Boolean(rule?.fromFile) } = {}) => ({
  feat: normalizeText(rule?.feat),
  qual: normalizeText(rule?.qual),
  val: normalizeText(rule?.val),
  color: normalizeColor(rule?.color),
  cap: normalizeText(rule?.cap),
  ...(fromFile ? { fromFile: true } : {})
});

export const specificRuleIdentity = (rule) => {
  const normalized = normalizeSpecificRule(rule);
  return [
    normalized.feat,
    normalized.qual,
    normalized.val,
    normalized.color,
    normalized.cap
  ];
};

const identityKey = (rule) => JSON.stringify(specificRuleIdentity(rule));

export const applySpecificRuleProvenance = (canonicalRules, storedRules) => {
  const storedFileIdentities = new Set(
    (Array.isArray(storedRules) ? storedRules : [])
      .filter((rule) => rule?.fromFile)
      .map(identityKey)
  );
  return (Array.isArray(canonicalRules) ? canonicalRules : []).map((rule) => {
    const normalized = normalizeSpecificRule(rule, { fromFile: false });
    return storedFileIdentities.has(identityKey(normalized))
      ? { ...normalized, fromFile: true }
      : normalized;
  });
};

// N-06 (PD-OI-042): Python keeps a rule's caption, but draws a rule whose
// caption names a renderer row of another color (a feature type, "other ...",
// or a numeric track) as "<caption> [<hex>]", made unique like a legend key
// (gbdraw/legend/table.py). The renderer rows are the Generate rows of the
// current legend that no rule draws. The live legend and the Legend editor use
// this one allocation to address the row Generate draws for a rule.

// Python's _legend_fill_identity: lower-case #rrggbb, or none.
const fillIdentity = (value) => {
  const color = normalizeColor(value);
  return /^#[0-9a-f]{3}$/.test(color) ? `#${[...color.slice(1)].map((char) => char + char).join('')}` : color;
};
const allocatedRowBase = (rule) => `${rule.cap} [${fillIdentity(rule.color)}]`;
const isAllocatedRow = (caption, rule) => {
  const base = allocatedRowBase(rule);
  return caption === base || (caption.startsWith(`${base} (`) && /^\d+\)$/.test(caption.slice(base.length + 2)));
};
const drawsRow = (row, rule) => row.color === fillIdentity(rule.color)
  && (row.caption === rule.cap || isAllocatedRow(row.caption, rule));

export const rendererLegendRows = ({ legendEntries = [], originalLegendOrder = [], rules = [] } = {}) => {
  const generated = new Set((originalLegendOrder || []).map(normalizeText).filter(Boolean));
  const normalizedRules = (rules || []).map((rule) => normalizeSpecificRule(rule)).filter((rule) => rule.cap);
  return (legendEntries || [])
    .map((entry) => ({
      caption: normalizeText(entry?.caption),
      origin: normalizeText(entry?.originalCaption || entry?.caption),
      color: fillIdentity(entry?.color)
    }))
    .filter((row) => row.caption && generated.has(row.origin)
      && !normalizedRules.some((rule) => drawsRow(row, rule)))
    .map(({ caption, color }) => ({ caption, color }));
};

const uniqueLegendKey = (reserved, preferred) => {
  if (!reserved.has(preferred)) return preferred;
  let suffix = 2;
  while (reserved.has(`${preferred} (${suffix})`)) suffix += 1;
  return `${preferred} (${suffix})`;
};

/**
 * A Legend row of the renderer: its caption and its fill (`fillIdentity`).
 * @typedef {object} RendererLegendRow
 * @property {string} caption
 * @property {string} color
 */

/**
 * A specific color rule as `normalizeSpecificRule` returns it.
 * @typedef {object} SpecificColorRule
 * @property {string} feat
 * @property {string} qual
 * @property {string} val
 * @property {string} color
 * @property {string} cap
 * @property {boolean} [fromFile]
 */

/**
 * Returns rule -> the legend caption Generate draws for that rule.
 * @param {Partial<SpecificColorRule>[]} [rules]
 * @param {RendererLegendRow[]} [rendererRows]
 * @returns {(rule: Partial<SpecificColorRule> | null | undefined) => string}
 */
export const createRuleLegendCaptions = (rules = [], rendererRows = []) => {
  const rowColors = new Map((rendererRows || []).map((row) => [row.caption, row.color]));
  const normalizedRules = (rules || []).map((rule) => normalizeSpecificRule(rule));
  const reserved = new Set([...rowColors.keys(), ...normalizedRules.map((rule) => rule.cap).filter(Boolean)]);
  const allocated = new Map();
  const keyOf = (rule) => JSON.stringify([rule.cap, fillIdentity(rule.color)]);
  normalizedRules
    .filter((rule) => rule.cap && rowColors.has(rule.cap) && rowColors.get(rule.cap) !== fillIdentity(rule.color))
    .sort((left, right) => keyOf(left).localeCompare(keyOf(right)))
    .forEach((rule) => {
      if (allocated.has(keyOf(rule))) return;
      const caption = uniqueLegendKey(reserved, allocatedRowBase(rule));
      reserved.add(caption);
      allocated.set(keyOf(rule), caption);
    });
  return (rule) => {
    const normalized = normalizeSpecificRule(rule);
    return allocated.get(keyOf(normalized)) || normalized.cap;
  };
};

// Only a caption that names a current row of another color can be allocated.
const mayBeAllocated = (rule, legendEntries) => Boolean(rule.cap) && (legendEntries || []).some((entry) => (
  normalizeText(entry?.caption) === rule.cap && fillIdentity(entry?.color) !== fillIdentity(rule.color)
));

// The legend caption Generate draws for one rule of `context.rules`.
export const ruleLegendCaption = (rule, { rules = [], legendEntries = [], originalLegendOrder = [] } = {}) => {
  const normalized = normalizeSpecificRule(rule);
  if (!mayBeAllocated(normalized, legendEntries)) return normalized.cap;
  return createRuleLegendCaptions(rules, rendererLegendRows({ legendEntries, originalLegendOrder, rules }))(normalized);
};

export const buildLegendIntents = (rules, rendererRows = []) => {
  const legendCaption = createRuleLegendCaptions(rules, rendererRows);
  const byCaption = new Map();
  for (const rule of rules || []) {
    const normalized = normalizeSpecificRule(rule);
    const caption = legendCaption(normalized);
    if (caption && !byCaption.has(caption)) {
      byCaption.set(caption, { caption, color: normalized.color });
    }
  }
  return { intents: [...byCaption.values()] };
};

// The rules a legend row draws, by the allocation above (never by reading a
// suffix back): editing that row edits these rules.
export const legendRowRules = (caption, { rules = [], legendEntries = [], originalLegendOrder = [] } = {}) => {
  const target = normalizeText(caption);
  if (!target) return [];
  const legendCaption = createRuleLegendCaptions(
    rules, rendererLegendRows({ legendEntries, originalLegendOrder, rules })
  );
  return (rules || []).filter((rule) => legendCaption(rule) === target);
};

export const diffLegendIntents = (currentEntries, desiredIntents) => {
  const desired = new Map(
    (Array.isArray(desiredIntents) ? desiredIntents : [])
      .map((entry) => /** @type {[string, string]} */ ([normalizeText(entry?.caption), normalizeColor(entry?.color)]))
      .filter(([caption]) => caption)
  );
  const current = new Map();
  (Array.isArray(currentEntries) ? currentEntries : []).forEach((entry) => {
    const caption = normalizeText(entry?.caption);
    if (!caption || current.has(caption)) return;
    current.set(caption, { ...entry, caption, color: normalizeColor(entry?.color) });
  });

  const diff = { add: [], update: [], remove: [], unchanged: [] };
  current.forEach((entry, caption) => {
    if (!desired.has(caption)) {
      diff.remove.push(entry);
    } else if (desired.get(caption) === entry.color) {
      diff.unchanged.push(entry);
    } else {
      diff.update.push({ ...entry, color: desired.get(caption) });
    }
    desired.delete(caption);
  });
  desired.forEach((color, caption) => diff.add.push({ caption, color }));
  return diff;
};

export const prepareSpecificColorImport = (text, currentRules = []) => {
  const { rules } = parseSpecificRules(text);
  const retainedRules = (Array.isArray(currentRules) ? currentRules : [])
    .filter((rule) => !rule?.fromFile)
    .map((rule) => normalizeSpecificRule(rule));
  const nextRules = [...retainedRules, ...rules.map((rule) => normalizeSpecificRule(rule, { fromFile: true }))];
  return {
    nextRules,
    importedCount: rules.length
  };
};

