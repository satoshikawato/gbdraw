// @ts-check
import { resolveColorToHex } from '../utils/color-utils.js';
import { parseSpecificRules } from './file-imports.js';
import { ruleMatcher } from './rule-matchers.js';
/** @import { PythonLegendKey, PythonLegendRow, SwatchColor } from './legend-svg.js' */

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
/** @param {readonly { cap?: unknown, color?: unknown }[]} rules */
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
// (gbdraw/legend/table.py). The renderer rows are the rows Python drew for the
// displayed Result that no rule draws. The live legend, the stroke reach, the
// feature popup and the Legend editor use this one allocation to address the
// row Generate draws for a rule.

// Python's _legend_fill_identity: lower-case #rrggbb, or none.
const fillIdentity = (value) => {
  const color = normalizeColor(value);
  return /^#[0-9a-f]{3}$/.test(color) ? `#${[...color.slice(1)].map((char) => char + char).join('')}` : color;
};
/** @param {unknown} value @returns {SwatchColor} */
const fillOf = (value) => /** @type {SwatchColor} */ (fillIdentity(value));
// The key of the row Python draws under a caption: the caption itself (a
// feature type's row by the type, a rule's row by the rule's caption), where
// N-06 allocates no "<caption> [<hex>]" key.
/** @param {string} caption @returns {PythonLegendKey} */
const captionRowKey = (caption) => /** @type {PythonLegendKey} */ (caption);
const allocatedRowBase = (rule) => `${rule.cap} [${fillIdentity(rule.color)}]`;
const isAllocatedRow = (caption, rule) => {
  const base = allocatedRowBase(rule);
  return caption === base || (caption.startsWith(`${base} (`) && /^\d+\)$/.test(caption.slice(base.length + 2)));
};
const drawsRow = (row, rule) => row.color === fillIdentity(rule.color)
  && (row.caption === rule.cap || isAllocatedRow(row.caption, rule));

/**
 * A displayed feature as the allocation reads it: its type and the fill Python
 * drew it with (a catalog row's `fill_color`), with its known rule matches.
 * @typedef {{ type?: unknown, fill_color?: unknown }} LegendRowFeature
 */

/**
 * The Legend state the allocation reads.
 * @typedef {object} LegendRowContext
 * @property {ReadonlyMap<PythonLegendKey, PythonLegendRow>} pythonRows Python's rows of the displayed
 *   Result (`pythonLegendRows`; none before the first Generate). A rule is compared with each row's key
 *   and the fill Python drew for it, as Python allocates (`generated_fills`, R15-4), never with the
 *   color a row's swatch shows (OV-294).
 * @property {readonly string[]} originalLegendOrder
 * @property {readonly Partial<SpecificColorRule>[]} rules
 * @property {Iterable<LegendRowFeature>} features The features the displayed Result draws
 *   (feature-visibility.js `displayedDrawnFeatures`), with their known rule matches (`rowsLeftOut`).
 */

// Python's own Legend rows (gbdraw/legend/table.py::_generated_legend_fills),
// read from a key: a track's row ("<NT> content" and "<NT> skew", each also
// as "(+)" and "(-)", and "Depth"), with the palette key of its fill ('' for
// Depth, which no palette key colors), else the row of the feature type it
// names, the type's own row or its "other" row, which the type's palette key
// colors. tests/fixtures/legend_generated_rows.json holds Python's keys.
const TRACK_ROW = /^[ACGTU]{2} (content|skew)(?: \(([+-])\))?$/;
/**
 * @param {string} key
 * @returns {{ track: boolean, paletteKey: string }}
 */
export const generatedLegendRow = (key) => {
  if (key === 'Depth') return { track: true, paletteKey: '' };
  const track = TRACK_ROW.exec(key);
  if (track) return { track: true, paletteKey: track[1] === 'content' ? 'gc_content' : `skew_${track[2] === '-' ? 'low' : 'high'}` };
  if (!key.toLowerCase().startsWith('other ')) return { track: false, paletteKey: key };
  const type = key.slice(6).trim();
  if (type.toLowerCase() === 'proteins') return { track: false, paletteKey: 'CDS' };
  return { track: false, paletteKey: type.endsWith('s') ? type.slice(0, -1) : type };
};

// The Python rows N-06 does not compare a rule with, read from the displayed
// features and their known rule matches (a match not known yet colors no
// feature), only when a captioned rule may need them:
// - OV-306: Python draws no row under a feature type's name once a captioned
//   rule of that type colors a feature (`_generated_legend_fills`: the type
//   has a rule whose caption and color a feature took).
// - OV-307: a row Python drew for a rule is no row of the renderer, also when
//   a live edit of the rule's color leaves the row in the color Python drew it
//   in until Generate: a feature Python drew in the row's color takes a rule
//   of the row's caption, and the caption names no row Python generates
//   (`generatedLegendRow`: a track's row, or a displayed type's row).
/**
 * @param {LegendRowContext} context
 * @returns {Set<PythonLegendKey>}
 */
const rowsLeftOut = ({ rules = [], pythonRows, features = [] }) => {
  const captioned = rules.map((rule) => normalizeSpecificRule(rule)).filter((rule) => rule.cap);
  const typed = captioned.filter((rule) => pythonRows.has(captionRowKey(rule.feat)));
  /** @type {Map<string, string>} A caption naming a Python row of another color, and that row's fill. */
  const recolored = new Map();
  captioned.forEach((rule) => {
    const color = fillIdentity(pythonRows.get(captionRowKey(rule.cap))?.color);
    if (color && color !== fillIdentity(rule.color)) recolored.set(rule.cap, color);
  });
  if (typed.length === 0 && recolored.size === 0) return new Set();
  const recoloredTypes = new Set(captioned.filter((rule) => recolored.has(rule.cap)).map((rule) => rule.feat));
  const winnerOf = ruleMatcher([...rules]).firstIfKnown;
  /** @param {Partial<SpecificColorRule>} rule */
  const pair = (rule) => JSON.stringify([normalizeText(rule.cap), normalizeColor(rule.color)]);
  const taken = new Set();
  /** @type {Set<string>} */
  const drawnForRules = new Set();
  /** @type {Set<string>} */
  const types = new Set();
  for (const feature of features) {
    const type = normalizeText(feature.type);
    types.add(type);
    if (typed.length === 0 && !recoloredTypes.has(type) && !recoloredTypes.has('*')) continue;
    const rule = winnerOf(feature);
    const caption = normalizeText(rule?.cap);
    if (!rule || !caption) continue;
    taken.add(pair(rule));
    if (recolored.get(caption) === fillIdentity(feature.fill_color)) drawnForRules.add(caption);
  }
  /** @param {string} caption */
  const generated = (caption) => {
    const row = generatedLegendRow(caption);
    return row.track || types.has(row.paletteKey);
  };
  return new Set([
    ...typed.filter((rule) => taken.has(pair(rule))).map((rule) => rule.feat),
    ...[...drawnForRules].filter((caption) => !generated(caption))
  ].map(captionRowKey));
};

// Python's rows of the displayed Result, by key and drawn fill: those N-06
// compares a rule with, and those the rules take (`rowsLeftOut`).
/** @param {LegendRowContext} context */
const pythonRowsOf = (context) => {
  const leftOut = rowsLeftOut(context);
  /** @type {RendererLegendRow[]} */
  const compared = [];
  /** @type {Set<PythonLegendKey>} */
  const taken = new Set();
  for (const row of context.pythonRows.values()) {
    if (!row.color) continue;
    if (leftOut.has(row.key)) taken.add(row.key);
    else compared.push({ caption: row.key, color: fillOf(row.color) });
  }
  return { compared, taken };
};

/**
 * @param {RendererLegendRow[]} rows
 * @param {LegendRowContext} context
 * @returns {RendererLegendRow[]}
 */
const rendererRowsOf = (rows, { originalLegendOrder = [], rules = [] }) => {
  const generated = new Set(originalLegendOrder.map(normalizeText).filter(Boolean));
  const normalizedRules = rules.map((rule) => normalizeSpecificRule(rule)).filter((rule) => rule.cap);
  return rows.filter((row) => row.caption && generated.has(row.caption)
    && !normalizedRules.some((rule) => drawsRow(row, rule)));
};

// The rows a commit of `context.rules` reads: the renderer rows of the N-06
// allocation, and the keys of the rows Python drew that those rules take,
// whose listed rows the commit replaces (`rowsLeftOut`: the row of a feature
// type once a captioned rule of the type colors a feature the displayed
// Result draws, OV-306, and a row Python drew for a rule, OV-307).
/**
 * @param {LegendRowContext} context
 * @returns {{ rendererRows: RendererLegendRow[], takenKeys: Set<PythonLegendKey> }}
 */
export const ruleCommitLegendRows = (context) => {
  const { compared, taken } = pythonRowsOf(context);
  return { rendererRows: rendererRowsOf(compared, context), takenKeys: taken };
};

const uniqueLegendKey = (reserved, preferred) => {
  if (!reserved.has(preferred)) return preferred;
  let suffix = 2;
  while (reserved.has(`${preferred} (${suffix})`)) suffix += 1;
  return `${preferred} (${suffix})`;
};

/**
 * A Legend row of the renderer: Python's key and the fill Python drew it in
 * (`fillIdentity`), never the caption or the swatch a Legend edit shows.
 * @typedef {object} RendererLegendRow
 * @property {PythonLegendKey} caption
 * @property {SwatchColor} color
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
 * Returns rule -> the key of the Legend row Generate draws for that rule.
 * @param {Partial<SpecificColorRule>[]} [rules]
 * @param {RendererLegendRow[]} [rendererRows]
 * @returns {(rule: Partial<SpecificColorRule> | null | undefined) => PythonLegendKey}
 */
export const createRuleLegendCaptions = (rules = [], rendererRows = []) => {
  /** @type {Map<PythonLegendKey, SwatchColor>} */
  const rowColors = new Map((rendererRows || []).map((row) => [row.caption, row.color]));
  const normalizedRules = (rules || []).map((rule) => normalizeSpecificRule(rule));
  /** @type {Set<string>} */
  const reserved = new Set([...rowColors.keys(), ...normalizedRules.map((rule) => rule.cap).filter(Boolean)]);
  /** @type {Map<string, PythonLegendKey>} */
  const allocated = new Map();
  const keyOf = (rule) => JSON.stringify([rule.cap, fillIdentity(rule.color)]);
  normalizedRules
    .filter((rule) => {
      const rowColor = rule.cap ? rowColors.get(captionRowKey(rule.cap)) : undefined;
      return rowColor !== undefined && rowColor !== fillIdentity(rule.color);
    })
    .sort((left, right) => keyOf(left).localeCompare(keyOf(right)))
    .forEach((rule) => {
      if (allocated.has(keyOf(rule))) return;
      const caption = uniqueLegendKey(reserved, allocatedRowBase(rule));
      reserved.add(caption);
      allocated.set(keyOf(rule), captionRowKey(caption));
    });
  return (rule) => {
    const normalized = normalizeSpecificRule(rule);
    return allocated.get(keyOf(normalized)) || captionRowKey(normalized.cap);
  };
};

// Only a caption that names a current row of another color can be allocated.
/** @param {ReturnType<typeof normalizeSpecificRule>} rule @param {RendererLegendRow[]} rows */
const mayBeAllocated = (rule, rows) => Boolean(rule.cap) && rows.some((row) => (
  row.caption === captionRowKey(rule.cap) && row.color !== fillIdentity(rule.color)
));

// The key of the Legend row Generate draws for a rule of `context.rules`,
// read for many rules: the keys are allocated once, when a rule first needs it.
/**
 * @param {LegendRowContext} context
 * @returns {(rule: Partial<SpecificColorRule> | null | undefined) => PythonLegendKey}
 */
export const ruleLegendCaptions = (context) => {
  const rows = pythonRowsOf(context).compared;
  /** @type {((rule: Partial<SpecificColorRule> | null | undefined) => PythonLegendKey) | null} */
  let allocated = null;
  return (rule) => {
    const normalized = normalizeSpecificRule(rule);
    if (!mayBeAllocated(normalized, rows)) return captionRowKey(normalized.cap);
    allocated ||= createRuleLegendCaptions([...context.rules || []], rendererRowsOf(rows, context));
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
// suffix back): editing that row edits these rules. The row is Python's key
// (`legendEntryKey` of a listed row), never the caption a rename shows.
/**
 * @param {PythonLegendKey | undefined} key
 * @param {LegendRowContext} context
 * @returns {Partial<SpecificColorRule>[]}
 */
export const legendRowRules = (key, context) => {
  const target = normalizeText(key);
  if (!target) return [];
  const rules = context.rules || [];
  const legendCaption = createRuleLegendCaptions([...rules], rendererRowsOf(pythonRowsOf(context).compared, context));
  return rules.filter((rule) => legendCaption(rule) === target);
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

