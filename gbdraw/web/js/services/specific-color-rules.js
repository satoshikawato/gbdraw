// @ts-check
import { resolveColorToHex } from '../utils/color-utils.js';
import { markedRows, parseSpecificRules, withoutImportedRows } from './file-imports.js';
import { legendColorIdentity as fillIdentity } from './legend-row-names.js';
import { ruleMatcher } from './rule-matchers.js';
import { isPlainObject } from '../utils/contract-shape.js';
/** @import { PythonLegendKey, PythonLegendRow, SwatchColor } from './legend-svg.js' */

/**
 * A rule as the normalization reads it: a table row, a stored rule, or a rule built from a Legend row.
 * @typedef {Readonly<{ feat?: unknown, qual?: unknown, val?: unknown, color?: unknown, cap?: unknown, fromFile?: unknown }> | null | undefined} SpecificRuleInput
 */

/** @param {unknown} value */
const normalizeText = (value) => String(value ?? '').trim();
/** @param {unknown} value */
const normalizeColor = (value) => String(resolveColorToHex(normalizeText(value)) || '').toLowerCase();

/**
 * @param {SpecificRuleInput} rule
 * @param {{ fromFile?: boolean }} [options]
 */
export const normalizeSpecificRule = (rule, { fromFile = Boolean(rule?.fromFile) } = {}) => ({
  feat: normalizeText(rule?.feat),
  qual: normalizeText(rule?.qual),
  val: normalizeText(rule?.val),
  color: normalizeColor(rule?.color),
  cap: normalizeText(rule?.cap),
  ...(fromFile ? { fromFile: true } : {})
});

/** @param {SpecificRuleInput} rule */
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

/** @param {SpecificRuleInput} rule */
const identityKey = (rule) => JSON.stringify(specificRuleIdentity(rule));

/**
 * @param {unknown} canonicalRules
 * @param {unknown} storedRules
 */
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

// The key of the Legend row a rule draws (PD-OI-042 scenario revision 2):
// Python keeps a rule's caption and draws it, numbered by the naming rule
// (legend-row-names.js). The row's key is the caption, or "<caption> [<fill>]",
// made unique like a legend key, when the caption is written with several
// colors in the rule table or names a renderer row of another swatch (a
// feature type, "other ...", or a numeric track; another fill, or another
// line than a rule row's, D-13): gbdraw/legend/table.py::
// rule_legend_keys, with tests/fixtures/legend_rule_keys.json holding both
// answers. The renderer rows are the rows Python drew for the displayed
// Result that no rule draws. The live legend, the stroke reach, the feature
// popup and the Legend editor use this one allocation to address the row
// Generate draws for a rule; the key is an identity and is never shown.
/** @param {unknown} value @returns {SwatchColor} */
const fillOf = (value) => /** @type {SwatchColor} */ (fillIdentity(value));
// The key of the row Python draws under a caption: the caption itself (a
// feature type's row by the type, a rule's row by the rule's caption), where
// N-06 allocates no "<caption> [<hex>]" key.
/** @param {string} caption @returns {PythonLegendKey} */
const captionRowKey = (caption) => /** @type {PythonLegendKey} */ (caption);
/** @param {Pick<SpecificColorRule, 'cap' | 'color'>} rule */
const allocatedRowBase = (rule) => `${rule.cap} [${fillIdentity(rule.color)}]`;
// The caption of a rule as a saved Session or a color table file is read
// (PD-OI-042 scenario revision 2, accepted residual risk). Builds from
// first-parent `main` fe6861f0 until the naming rule saved the key above as
// the caption of a rule whose caption was written with several colors. A
// caption that ends in the color ID of its own rule is read without it, once;
// one that would be left blank is read as written. Python's reading is
// gbdraw/legend/table.py::saved_rule_caption; both answer
// tests/fixtures/saved_rule_captions.json.
/** @param {string} caption @param {unknown} color */
export const savedRuleCaption = (caption, color) => {
  const colorId = allocatedRowBase({ cap: '', color: String(color ?? '') });
  const name = caption.endsWith(colorId) ? caption.slice(0, -colorId.length) : caption;
  return name.trim() ? name : caption;
};
/**
 * Saved rules with each caption read by `savedRuleCaption`; an entry that is
 * no object stays (the rules' reader drops it, listed).
 * @template {{ readonly cap?: unknown, readonly color?: unknown }} Rule
 * @param {readonly Rule[]} rules
 * @returns {Rule[]}
 */
export const savedSpecificRules = (rules) => rules.map((rule) => {
  if (!rule || typeof rule !== 'object') return rule;
  const caption = String(rule.cap ?? '');
  const read = savedRuleCaption(caption, rule.color);
  return read === caption ? rule : { ...rule, cap: read };
});
/**
 * The color of a per-feature color edit: an edit object's `color`, or an
 * older bare color.
 * @param {unknown} edit
 * @returns {unknown}
 */
export const featureColorEditColor = (edit) => (
  edit && typeof edit === 'object' && Object.hasOwn(edit, 'color') ? /** @type {{ color?: unknown }} */ (edit).color : edit
);
// Per-feature color edits as a saved Session is read: an older build wrote a
// rule's rewritten caption into the edits of the features the rule colored,
// so an edit's caption is read as its rule's is, by the edit's own color.
/**
 * @param {unknown} edits Unvalidated: a map of no object is no edit, as its
 *   reader (`readFeatureColorEdits`) lists it (OV-866).
 * @returns {Record<string, unknown>}
 */
export const savedFeatureColorEdits = (edits) => Object.fromEntries(Object.entries(isPlainObject(edits) ? edits : {}).map(
  /** @param {[string, unknown]} entry @returns {[string, unknown]} */ ([key, edit]) => {
  if (!edit || typeof edit !== 'object' || typeof (/** @type {{ caption?: unknown }} */ (edit)).caption !== 'string') return [key, edit];
  const { caption, color } = /** @type {{ caption: string, color?: unknown }} */ (edit);
  const read = savedRuleCaption(caption, color);
  return [key, read === caption ? edit : { ...edit, caption: read }];
}));
/**
 * @param {string} caption
 * @param {Pick<SpecificColorRule, 'cap' | 'color'>} rule
 */
const isAllocatedRow = (caption, rule) => {
  const base = allocatedRowBase(rule);
  return caption === base || (caption.startsWith(`${base} (`) && /^\d+\)$/.test(caption.slice(base.length + 2)));
};
// Whether a row has the swatch a rule row of `color` is drawn with: that fill,
// and a rule row's line (Python's `legend_swatch` of both, D-13).
/** @param {RendererLegendRow} row @param {unknown} color */
const hasRuleSwatch = (row, color) => row.color === fillIdentity(color) && row.ruleLine !== false;
/**
 * @param {RendererLegendRow} row
 * @param {Pick<SpecificColorRule, 'cap' | 'color'>} rule
 */
const drawsRow = (row, rule) => hasRuleSwatch(row, rule.color)
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

// Python's own Legend rows (gbdraw/legend/table.py::_generated_legend_rows),
// read from a key: a default track's row ("<NT> content" and "<NT> skew", each
// also as "(+)" and "(-)", and "Depth"), else the row of the feature type it
// names, the type's own row or its "other" row, which the type's palette key
// colors. A track row has no palette key here: the track rows of the request
// that drew the Result give its fill (track-legend-rows.js, OV-484).
// tests/fixtures/legend_generated_rows.json holds Python's keys.
const TRACK_ROW = /^[ACGTU]{2} (content|skew)(?: \([+-]\))?$/;
/**
 * @param {string} key
 * @returns {{ track: boolean, paletteKey: string }}
 */
export const generatedLegendRow = (key) => {
  if (key === 'Depth' || TRACK_ROW.test(key)) return { track: true, paletteKey: '' };
  if (!key.toLowerCase().startsWith('other ')) return { track: false, paletteKey: key };
  const type = key.slice(6).trim();
  if (type.toLowerCase() === 'proteins') return { track: false, paletteKey: 'CDS' };
  return { track: false, paletteKey: type.endsWith('s') ? type.slice(0, -1) : type };
};

// The Python rows N-06 does not compare a rule with, read from the displayed
// features and their known rule matches (a match not known yet colors no
// feature), only when a captioned rule may need them:
// - OV-306: Python draws no row under a feature type's name once a captioned
//   rule of that type colors a feature (`_generated_legend_rows`: the type
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
    else compared.push({ caption: row.key, color: fillOf(row.color), ruleLine: row.ruleLine });
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

// Python's `_unique_legend_key` (gbdraw/legend/table.py).
/** @param {ReadonlySet<string>} reserved @param {string} preferred @returns {string} */
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
 * @property {boolean} [ruleLine] False when Python drew its swatch with another line than a rule row's
 *   (`PythonLegendRow.ruleLine`).
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
 * How a rule table write (rule-actions.js) builds the next table from the one at its turn (OV-444);
 * a row returned as it is or as `revise(row, patch)` keeps its identity. Null writes nothing.
 * @callback SpecificRuleBuild
 * @param {readonly SpecificColorRule[]} current
 * @param {(row: SpecificColorRule, patch: Partial<SpecificColorRule>) => SpecificColorRule} revise
 * @returns {Partial<SpecificColorRule>[] | null}
 */

// The rules of one caption and one fill draw one Legend row: that pair is the
// identity of their Legend group, for the row allocation and for the groups below.
/** @param {Pick<SpecificColorRule, 'cap' | 'color'>} rule */
const legendGroupKey = (rule) => JSON.stringify([rule.cap, fillIdentity(rule.color)]);

// The fills each non-blank caption of a rule table is written with.
/** @param {readonly Pick<SpecificColorRule, 'cap' | 'color'>[]} rules @returns {Map<string, Set<string>>} */
const captionFills = (rules) => {
  /** @type {Map<string, Set<string>>} */
  const fills = new Map();
  rules.forEach((rule) => {
    if (rule.cap) fills.set(rule.cap, (fills.get(rule.cap) || new Set()).add(fillIdentity(rule.color)));
  });
  return fills;
};

/**
 * The name of a color-rule Legend group: the caption its rules carry, never
 * the key of the row Generate draws for them (N-06).
 * @typedef {string & { readonly __brand: 'ColorRuleGroupName' }} ColorRuleGroupName
 */
/**
 * A color-rule Legend group (PD-OI-061 revision 4): the rules of one caption
 * and one color, which share one Legend row whether a feature takes them or
 * not. Its name and its color are its identity.
 * @typedef {object} ColorRuleLegendGroup
 * @property {ColorRuleGroupName} name
 * @property {SwatchColor} color
 * @property {readonly string[]} types The feature types of its rules, in rule order (`*`: every type).
 */

// A group takes a feature of a type one of its rules is for; a rule for every
// type (`*`) is for each type. The popup lists the groups that take the clicked
// feature, and a pick joins every group it lists (PD-OI-061 revision 4).
/** @param {ColorRuleLegendGroup} group @param {string} type */
export const colorRuleGroupTakes = (group, type) => group.types.includes(type) || group.types.includes('*');

// The color-rule Legend groups of a rule table, each once, in rule order; with
// `type`, those that take a feature of that type. The feature popup's list and
// the target of a Legend rename read this one derivation.
/**
 * @param {readonly SpecificRuleInput[]} rules
 * @param {string} [type]
 * @returns {ColorRuleLegendGroup[]}
 */
export const colorRuleLegendGroups = (rules, type) => {
  /** @type {Map<string, { name: ColorRuleGroupName, color: SwatchColor, types: string[] }>} */
  const groups = new Map();
  for (const rule of rules) {
    const normalized = normalizeSpecificRule(rule);
    if (!normalized.cap) continue;
    const key = legendGroupKey(normalized);
    const group = groups.get(key)
      || { name: /** @type {ColorRuleGroupName} */ (normalized.cap), color: fillOf(normalized.color), types: [] };
    if (!group.types.includes(normalized.feat)) group.types.push(normalized.feat);
    groups.set(key, group);
  }
  return [...groups.values()].filter((group) => type === undefined || colorRuleGroupTakes(group, type));
};

/**
 * Returns rule -> the key of the Legend row Generate draws for that rule.
 * @param {Partial<SpecificColorRule>[] | null} [rules]
 * @param {RendererLegendRow[]} [rendererRows]
 * @returns {(rule: SpecificRuleInput) => PythonLegendKey}
 */
export const createRuleLegendCaptions = (rules = [], rendererRows = []) => {
  /** @type {Map<string, RendererLegendRow>} */
  const rows = new Map((rendererRows || []).map((row) => [row.caption, row]));
  const normalizedRules = (rules || []).map((rule) => normalizeSpecificRule(rule)).filter((rule) => rule.cap);
  const fillsOf = captionFills(normalizedRules);
  /** @type {Set<string>} */
  const reserved = new Set([...rows.keys(), ...fillsOf.keys()]);
  /** @type {Map<string, PythonLegendKey>} */
  const allocated = new Map();
  // Python allocates in sorted order; the order never changes a key, since
  // each caption and fill has its own "<caption> [<fill>]".
  [...fillsOf.keys()].sort().forEach((cap) => {
    const fills = /** @type {Set<string>} */ (fillsOf.get(cap));
    const row = rows.get(cap);
    [...fills].sort().forEach((color) => {
      if (fills.size < 2 && (!row || hasRuleSwatch(row, color))) return;
      const key = uniqueLegendKey(reserved, allocatedRowBase({ cap, color }));
      reserved.add(key);
      allocated.set(legendGroupKey({ cap, color }), captionRowKey(key));
    });
  });
  return (rule) => {
    const normalized = normalizeSpecificRule(rule);
    return allocated.get(legendGroupKey(normalized)) || captionRowKey(normalized.cap);
  };
};

// The rule rows a change of the rule table re-keys (OV-645): a group of one
// caption and one fill that both tables hold is one Legend row, and its key
// changes when its caption starts or stops being written in several colors.
// Old key to new key; the Legend owner keeps such a row in its place.
/**
 * @param {readonly Partial<SpecificColorRule>[]} before
 * @param {readonly Partial<SpecificColorRule>[]} after
 * @param {RendererLegendRow[]} [rendererRows]
 * @returns {Map<PythonLegendKey, PythonLegendKey>}
 */
export const rekeyedRuleRows = (before, after, rendererRows = []) => {
  const kept = new Set(after.map((rule) => legendGroupKey(normalizeSpecificRule(rule))));
  const groups = before.map((rule) => normalizeSpecificRule(rule)).filter((group) => kept.has(legendGroupKey(group)));
  return changedRowKeys(groups.map((group) => [group, group]),
    createRuleLegendCaptions([...before], rendererRows), createRuleLegendCaptions([...after], rendererRows));
};

// The Legend rows a Remove of rules takes from a drawing not shown, and the
// rows it re-keys, read from the keys the drawing's Legend stores (OV-1280).
// Its Result is not mounted, so Python's fills and swatch lines are unknown,
// and no key is allocated again from the colors its rows list: a Legend
// recolor, or a type row drawn in a rule's fill with its own line (D-13),
// would read as a rule row. A rule group's row is its allocated key
// ("<cap> [<fill>]", or "<cap> [<fill>] (n)") when the drawing lists one,
// else its caption. `listedKeys` are the keys of every row the drawing lists,
// deleted and dormant rows too, which Python still draws. `gone`: the rows
// only removed groups drew; `rekeyed`: a kept group whose caption is now
// written in one color moves to the caption, unless a row that is no rule
// row has it (as `rekeyedRuleRows` allocates).
/**
 * @param {readonly Partial<SpecificColorRule>[]} before
 * @param {readonly Partial<SpecificColorRule>[]} after The rules of `before` that stay.
 * @param {Iterable<string>} listedKeys
 * @returns {{ gone: RuleRowColor[], rekeyed: Map<PythonLegendKey, PythonLegendKey> }}
 */
export const listedRuleRowRemoval = (before, after, listedKeys) => {
  const listed = [...listedKeys].map(captionRowKey);
  /** @param {readonly Partial<SpecificColorRule>[]} rules @returns {Map<string, SpecificColorRule>} */
  const groupsOf = (rules) => new Map(rules.map((rule) => normalizeSpecificRule(rule))
    .filter((rule) => rule.cap).map((rule) => [legendGroupKey(rule), rule]));
  const removed = groupsOf(before);
  const kept = groupsOf(after);
  kept.forEach((_group, id) => removed.delete(id));
  const captions = new Set([...removed.values(), ...kept.values()].map((group) => group.cap));
  /** @param {SpecificColorRule} group @returns {PythonLegendKey} */
  const storedKey = (group) => captionRowKey(listed.find((key) => !captions.has(key) && isAllocatedRow(key, group)) ?? group.cap);
  const ruleRows = new Set([...removed.values(), ...kept.values()].map(storedKey));
  const fills = captionFills([...kept.values()]);
  /** @type {Map<PythonLegendKey, PythonLegendKey>} */
  const rekeyed = new Map();
  kept.forEach((group) => {
    const from = storedKey(group);
    if (from !== group.cap && fills.get(group.cap)?.size === 1 && !listed.some((key) => key === group.cap && !ruleRows.has(key))) {
      rekeyed.set(from, captionRowKey(group.cap));
    }
  });
  const held = new Set([...[...kept.values()].map(storedKey), ...rekeyed.values()]);
  const gone = [...removed.values()].map((group) => ({ caption: storedKey(group), color: group.color }))
    .filter(({ caption }) => !held.has(caption));
  return { gone, rekeyed };
};

/**
 * A rule row: its Legend key and the color its rules draw it in.
 * @typedef {{ caption: PythonLegendKey, color: string }} RuleRowColor
 */

// OV-152: the keys of the Legend colors that copy a removed rule row's color,
// which go with the row: left behind, a copy would paint the row Generate
// draws for the features' palette color. A copy is stored under the row's key
// in the row's color (the popup stores the caption's color when it recolors a
// whole group, OV-61). A color set in the Legend editor is no copy: the editor
// recolors the rules of a row a rule draws (legend.js), so its colors sit on
// rows no rule draws and stay. Each removed row is read on its own, so two
// removed groups stored under one key each retire their own copy.
/**
 * @param {Readonly<Record<string, string>>} overrides
 * @param {readonly RuleRowColor[]} removedRows
 * @returns {PythonLegendKey[]}
 */
export const legendColorCopies = (overrides, removedRows) => {
  /** @type {PythonLegendKey[]} */
  const copies = [];
  for (const { caption, color } of removedRows) {
    if (Object.hasOwn(overrides, caption) && resolveColorToHex(overrides[caption]) === resolveColorToHex(color)) copies.push(caption);
  }
  return copies;
};

// The rule rows the reading of saved captions re-keys (F2): a caption saved
// with its own rule's color ID is read without it (`savedSpecificRules`), and
// when the read caption has one swatch, its row's key is the caption, no
// longer the text saved. `context.rules` are the rules as saved; each is
// paired with its own reading. A caption still written in several colors, or
// naming a row of another swatch, keeps its key.
/**
 * @param {LegendRowContext} context
 * @returns {Map<PythonLegendKey, PythonLegendKey>}
 */
export const savedRuleRowRekeys = (context) => {
  const saved = [...context.rules || []];
  const read = savedSpecificRules(saved);
  return changedRowKeys(saved.map((rule, index) => [rule, read[index]]),
    ruleLegendCaptions({ ...context, rules: saved }), ruleLegendCaptions({ ...context, rules: read }));
};

// Old key to new key of each pair of one rule before and after a change whose key differs.
/**
 * @param {readonly [Partial<SpecificColorRule>, Partial<SpecificColorRule>][]} pairs
 * @param {(rule: Partial<SpecificColorRule>) => PythonLegendKey} keyBefore
 * @param {(rule: Partial<SpecificColorRule>) => PythonLegendKey} keyAfter
 * @returns {Map<PythonLegendKey, PythonLegendKey>}
 */
const changedRowKeys = (pairs, keyBefore, keyAfter) => {
  /** @type {Map<PythonLegendKey, PythonLegendKey>} */
  const rekeyed = new Map();
  for (const [before, after] of pairs) {
    if (!normalizeText(before.cap) || !normalizeText(after.cap)) continue;
    const from = keyBefore(before);
    const to = keyAfter(after);
    if (from !== to) rekeyed.set(from, to);
  }
  return rekeyed;
};

// Only a caption written with several colors, or one that names a current row
// of another swatch, can be allocated.
/** @param {ReturnType<typeof normalizeSpecificRule>} rule @param {RendererLegendRow[]} rows */
const namesRowOfAnotherSwatch = (rule, rows) => Boolean(rule.cap) && rows.some((row) => (
  row.caption === captionRowKey(rule.cap) && !hasRuleSwatch(row, rule.color)
));

// The key of the Legend row Generate draws for a rule of `context.rules`,
// read for many rules: the keys are allocated once, when a rule first needs it.
/**
 * @param {LegendRowContext} context
 * @returns {(rule: Partial<SpecificColorRule> | null | undefined) => PythonLegendKey}
 */
export const ruleLegendCaptions = (context) => {
  const rows = pythonRowsOf(context).compared;
  /** @type {Map<string, Set<string>> | null} */
  let fillsOf = null;
  /** @type {((rule: Partial<SpecificColorRule> | null | undefined) => PythonLegendKey) | null} */
  let allocated = null;
  return (rule) => {
    const normalized = normalizeSpecificRule(rule);
    fillsOf ||= captionFills((context.rules || []).map((each) => normalizeSpecificRule(each)));
    if ((fillsOf.get(normalized.cap)?.size || 0) < 2 && !namesRowOfAnotherSwatch(normalized, rows)) return captionRowKey(normalized.cap);
    allocated ||= createRuleLegendCaptions([...context.rules || []], rendererRowsOf(rows, context));
    return allocated(normalized);
  };
};

// The key of the Legend row each color-rule group of `context.rules` is drawn
// in: a group is the rules of one caption and one color, which share that row.
/**
 * @param {LegendRowContext} context
 * @returns {(group: Pick<ColorRuleLegendGroup, 'name' | 'color'>) => PythonLegendKey}
 */
export const groupLegendRowKeys = (context) => {
  const rowOf = ruleLegendCaptions(context);
  return (group) => rowOf({ cap: group.name, color: group.color });
};

// The Legend rows of a rule table: one per row key, for the rules `draws`
// admits (a commit: those a displayed feature takes). The keys are allocated
// from the whole table, as Python allocates them, whichever rules draw a row.
/**
 * @param {Partial<SpecificColorRule>[] | null | undefined} rules
 * @param {RendererLegendRow[]} [rendererRows]
 * @param {(rule: Partial<SpecificColorRule>) => boolean} [draws]
 */
export const buildLegendIntents = (rules, rendererRows = [], draws = () => true) => {
  const legendCaption = createRuleLegendCaptions(rules, rendererRows);
  // An intent is a rule's Legend row: its key (`caption`), its color, and its
  // name (the rule's caption) where that is not the key, as `data-legend-name`.
  /** @type {Map<PythonLegendKey, { caption: PythonLegendKey, name?: string, color: string }>} */
  const byCaption = new Map();
  for (const rule of rules || []) {
    if (!draws(rule)) continue;
    const normalized = normalizeSpecificRule(rule);
    const caption = legendCaption(normalized);
    if (caption && !byCaption.has(caption)) {
      byCaption.set(caption, { caption, color: normalized.color, ...(normalized.cap === caption ? {} : { name: normalized.cap }) });
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

/** @typedef {Record<string, any> & { caption: string, color: string }} DiffEntry A Legend entry with its caption and color normalized; its other fields are kept. */
/**
 * @param {readonly { caption?: unknown, color?: unknown }[] | null | undefined} currentEntries
 * @param {readonly { caption?: unknown, color?: unknown }[] | null | undefined} desiredIntents
 */
export const diffLegendIntents = (currentEntries, desiredIntents) => {
  const desired = new Map(
    (Array.isArray(desiredIntents) ? desiredIntents : [])
      .map((entry) => /** @type {[string, string]} */ ([normalizeText(entry?.caption), normalizeColor(entry?.color)]))
      .filter(([caption]) => caption)
  );
  /** @type {Map<string, DiffEntry>} */
  const current = new Map();
  (Array.isArray(currentEntries) ? currentEntries : []).forEach((entry) => {
    const caption = normalizeText(entry?.caption);
    if (!caption || current.has(caption)) return;
    current.set(caption, { ...entry, caption, color: normalizeColor(entry?.color) });
  });

  /** @type {Record<'add' | 'update' | 'remove' | 'unchanged', DiffEntry[]>} */
  const diff = { add: [], update: [], remove: [], unchanged: [] };
  current.forEach((entry, caption) => {
    if (!desired.has(caption)) {
      diff.remove.push(entry);
    } else if (desired.get(caption) === entry.color) {
      diff.unchanged.push(entry);
    } else {
      diff.update.push({ ...entry, color: /** @type {string} */ (desired.get(caption)) }); // `desired.has(caption)` held above
    }
    desired.delete(caption);
  });
  desired.forEach((color, caption) => diff.add.push({ caption, color }));
  return diff;
};

/**
 * @param {string | null | undefined} text
 * @param {readonly SpecificRuleInput[]} [currentRules]
 */
export const prepareSpecificColorImport = (text, currentRules = []) => {
  const { rules } = parseSpecificRules(text);
  // The file's rows replace the rows of the card's earlier file (OV-1280).
  const retainedRules = withoutImportedRows((Array.isArray(currentRules) ? currentRules : [])
    .map((rule) => normalizeSpecificRule(rule)), 't_color');
  const nextRules = [...retainedRules, ...markedRows(savedSpecificRules(rules), 't_color').map((rule) => normalizeSpecificRule(rule))];
  return {
    nextRules,
    importedCount: rules.length
  };
};

