// @ts-check
import { resolveColorToHex } from '../utils/color-utils.js';
import { parseSpecificRules } from './file-imports.js';
/** @import { PythonLegendKey, PythonLegendRow } from './legend-svg.js' */

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

// Whether Python's caption normalization
// (gbdraw/features/colors.py::normalize_specific_color_captions) returns these
// rules as they are: it changes captions only when a non-blank caption has two
// or more colors, so a table whose every caption has one color string is its
// own normalization (equal strings normalize equally; one direction only, so
// `red` and `#ff0000` under one caption still ask Python). R4: the vectors in
// tests/fixtures/specific_color_caption_fixpoints.json hold both answers.
/** @param {Record<string, any>[]} rules */
export const ruleCaptionsAreNormalized = (rules) => {
  /** @type {Map<string, string>} */
  const colorOf = new Map();
  return rules.every((rule) => {
    const caption = String(rule?.cap ?? '');
    if (!caption.trim()) return true;
    const color = String(rule?.color ?? '').trim();
    if (!colorOf.has(caption)) colorOf.set(caption, color);
    return colorOf.get(caption) === color;
  });
};

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

/**
 * The Legend state the allocation reads.
 * @typedef {object} LegendRowContext
 * @property {Record<string, any>[]} [legendEntries]
 * @property {string[]} [originalLegendOrder]
 * @property {Partial<SpecificColorRule>[]} [rules]
 * @property {ReadonlyMap<PythonLegendKey, PythonLegendRow>} [pythonRows] Python's rows of the displayed
 *   Result. Given, a rule is compared with each row's key and the fill Python drew for it, as Python
 *   allocates (`generated_fills`, R15-4); without it, with the listed rows' captions and swatches.
 */

// The rows N-06 compares a rule with: Python's (key and drawn fill) when the
// caller gives them, else the listed entries (caption and swatch).
/** @param {LegendRowContext} context */
const comparedRows = ({ legendEntries = [], pythonRows }) => (pythonRows
  ? [...pythonRows.values()].flatMap((row) => (row.color
    ? [{ caption: String(row.key), origin: String(row.key), color: fillIdentity(row.color) }]
    : []))
  : (legendEntries || []).map((entry) => ({
    caption: normalizeText(entry?.caption),
    origin: normalizeText(entry?.originalCaption || entry?.caption),
    color: fillIdentity(entry?.color)
  })));

/**
 * @param {ReturnType<typeof comparedRows>} rows
 * @param {LegendRowContext} context
 * @returns {RendererLegendRow[]}
 */
const rendererRowsOf = (rows, { originalLegendOrder = [], rules = [] }) => {
  const generated = new Set((originalLegendOrder || []).map(normalizeText).filter(Boolean));
  const normalizedRules = (rules || []).map((rule) => normalizeSpecificRule(rule)).filter((rule) => rule.cap);
  return rows
    .filter((row) => row.caption && generated.has(row.origin)
      && !normalizedRules.some((rule) => drawsRow(row, rule)))
    .map(({ caption, color }) => ({ caption, color }));
};

/** @param {LegendRowContext} [context] */
export const rendererLegendRows = (context = {}) => rendererRowsOf(comparedRows(context), context);

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
/** @param {ReturnType<typeof normalizeSpecificRule>} rule @param {ReturnType<typeof comparedRows>} rows */
const mayBeAllocated = (rule, rows) => Boolean(rule.cap) && rows.some((row) => (
  row.caption === rule.cap && row.color !== fillIdentity(rule.color)
));

// The legend caption Generate draws for a rule of `context.rules`, read for
// many rules: the captions are allocated once, when a rule first needs it.
/**
 * @param {LegendRowContext} [context]
 * @returns {(rule: Partial<SpecificColorRule> | null | undefined) => string}
 */
export const ruleLegendCaptions = (context = {}) => {
  const rows = comparedRows(context);
  /** @type {((rule: Partial<SpecificColorRule> | null | undefined) => string) | null} */
  let allocated = null;
  return (rule) => {
    const normalized = normalizeSpecificRule(rule);
    if (!mayBeAllocated(normalized, rows)) return normalized.cap;
    allocated ||= createRuleLegendCaptions(context.rules, rendererRowsOf(rows, context));
    return allocated(normalized);
  };
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
/**
 * @param {string} caption
 * @param {LegendRowContext} [context]
 */
export const legendRowRules = (caption, context = {}) => {
  const target = normalizeText(caption);
  if (!target) return [];
  const legendCaption = createRuleLegendCaptions(context.rules, rendererLegendRows(context));
  return (context.rules || []).filter((rule) => legendCaption(rule) === target);
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

  /** @type {Record<'add' | 'update' | 'remove' | 'unchanged', Array<Record<string, any> & { caption: string, color: string }>>} */
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

