// @ts-check
import { currentPaletteName, isTypeColorValue, resolveColorToHex, ruleFill } from '../utils/color-utils.js';
import { TYPE_COLOR_VALUE } from '../contracts/schema.generated.js';
import { MODE_SLICE_DEFAULTS } from '../mode-slice-defaults.generated.js';
import { markedRows, parseSpecificRules, withoutImportedRows } from './file-imports.js';
import { legendColorIdentity as fillIdentity } from './legend-row-names.js';
import { ruleMatcher } from './rule-matchers.js';
import { isPlainObject } from '../utils/contract-shape.js';
/** @import { PythonLegendKey, PythonLegendRow, SwatchColor } from './legend-svg.js' */
/** @import { LegendOrderRuleRow } from '../contracts/rule-synthesis.generated.js' */

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

/**
 * OV-2435: whether a rule names a qualifier, the one check of the Add form, a row edit, and a row's hint. Python
 * files a rule under its qualifier and reads only the qualifiers a feature has, so a rule without one colors
 * nothing; the CLI color table refuses it (`read_color_table` in `gbdraw/io/colors.py`).
 * @param {SpecificRuleInput} rule
 * @returns {boolean}
 */
export const specificRuleHasQualifier = (rule) => normalizeText(rule?.qual) !== '';

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
// Whether a rule of a saved Session is one Reset fill color wrote before it
// wrote the type color value: a `hash` rule named like its type's "other"
// Legend row in the type's color (`ruleFill` of the type color value, so the
// fallback color for a type the colors do not hold) in one of the type colors
// of the Session's drawing (`typeColors()`, called only for a rule so named;
// `savedResetTypeColors`). The Web of first-parent `main` since
// a6aa233f3 (0.14.5) wrote it for a feature another rule still colored; "This
// feature only" under the scope dialog's default name writes that name in a
// picked color, so the color must match. Python's reading is
// gbdraw/legend/table.py::is_saved_reset_rule; both answer
// tests/fixtures/saved_rule_captions.json (`rules`).
/**
 * @param {{ readonly feat?: unknown, readonly qual?: unknown, readonly color?: unknown }} rule
 * @param {string} caption The caption as `savedRuleCaption` reads it.
 * @param {() => readonly FeatureTypeColors[]} typeColors
 */
const isSavedResetRule = (rule, caption, typeColors) => {
  const type = normalizeText(rule.feat);
  return normalizeText(rule.qual).toLowerCase() === 'hash' && caption === otherFeatureRowKey(type)
    && typeColors().some((colors) => fillIdentity(ruleFill(TYPE_COLOR_VALUE, type, colors)) === fillIdentity(rule.color));
};
// The type colors a reset rule of a Session drawn in `typeColors` with the
// palette `palette` may have been written in (`isSavedResetRule`):
// `typeColors`, and `typeColors` with the rows the palette held in the release
// that wrote such rules (`earlierPaletteRows`, OV-3259). Python's twin is
// gbdraw/session_drawings.py::saved_reset_type_colors; both answer
// tests/fixtures/saved_rule_captions.json (`rules`).
/**
 * @param {FeatureTypeColors} typeColors
 * @param {unknown} palette
 * @returns {FeatureTypeColors[]}
 */
export const savedResetTypeColors = (typeColors, palette) => {
  const name = currentPaletteName(String(palette || 'default'));
  const earlier = Object.hasOwn(MODE_SLICE_DEFAULTS.earlierPaletteRows.palettes, name)
    ? MODE_SLICE_DEFAULTS.earlierPaletteRows.palettes[name] : null;
  return earlier ? [typeColors, { ...typeColors, ...earlier }] : [typeColors];
};
/**
 * Saved rules with each caption read by `savedRuleCaption`, and, given the
 * type colors of the saved drawing (a Session; an imported file gives none),
 * each rule `isSavedResetRule` names read as the type color value without a
 * caption. An entry that is no object stays (the rules' reader drops it, listed).
 * @template {{ readonly feat?: unknown, readonly qual?: unknown, readonly cap?: unknown, readonly color?: unknown }} Rule
 * @param {readonly Rule[]} rules
 * @param {() => readonly FeatureTypeColors[]} [typeColors] The type colors a reset rule may have been written in
 *   (`savedResetTypeColors`).
 * @returns {Rule[]}
 */
export const savedSpecificRules = (rules, typeColors) => {
  /** @type {readonly FeatureTypeColors[] | null} */
  let colors = null;
  const colorsOnce = () => (colors ||= /** @type {() => readonly FeatureTypeColors[]} */ (typeColors)());
  return rules.map((rule) => {
    if (!rule || typeof rule !== 'object') return rule;
    const caption = String(rule.cap ?? '');
    const read = savedRuleCaption(caption, rule.color);
    if (typeColors && isSavedResetRule(rule, read, colorsOnce)) return { ...rule, cap: '', color: TYPE_COLOR_VALUE };
    return read === caption ? rule : { ...rule, cap: read };
  });
};
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
 * @property {FeatureTypeColors} palette The color each feature type takes without a rule, as the displayed
 *   Result applies it (`appliedFeatureColors`); the feature rows read it (`featureLegendRows`).
 */

/**
 * The colors of the feature types, by type and `default` (color-utils.js `appliedFeatureColors`).
 * @typedef {Readonly<Record<string, string>>} FeatureTypeColors
 */

/**
 * The color rules the drawn features take, read in one pass of `features`.
 * @typedef {object} DrawnRuleMatches
 * @property {boolean} known Whether every drawn feature's match is known. Each reader says what an
 *   unknown match means to it.
 * @property {ReadonlyMap<Partial<SpecificColorRule>, readonly LegendRowFeature[]>} takenBy Each captioned
 *   rule of `rules` that some drawn feature with a known match takes as its first rule, and those features.
 * @property {FeatureLegendRows} rows The feature rows Generate draws for them, a feature of an unknown
 *   match read as one no rule colors.
 */

/**
 * A drawn feature as the feature rows read it: its type, and the rule it takes (null: no rule; it is drawn
 * in its type's color).
 * @typedef {{ type: string, winner: SpecificRuleInput | null }} LegendRowUse
 */

/**
 * The row a feature is in: a type's own row or its "other" row, by its key; a captioned rule's row, as the
 * caption and color of its rules (its key: `ruleLegendCaptions`, tests/fixtures/legend_rule_keys.json); null
 * for none.
 * @typedef {PythonLegendKey | Pick<SpecificColorRule, 'cap' | 'color'> | null} FeatureLegendRow
 */

/**
 * The feature rows of a Legend (`featureLegendRows`).
 * @typedef {object} FeatureLegendRows
 * @property {readonly string[]} present The drawn types, in the order they first appear.
 * @property {(type: string) => boolean} rulesDraw Whether features of the type take a captioned rule of the
 *   type (or of every type): its features of no rule then form its "other" row.
 * @property {(type: string) => PythonLegendKey | null} typeRow The renderer's own row of a drawn type.
 * @property {(type: string, winner: SpecificRuleInput | null) => FeatureLegendRow} rowOf The row a feature of
 *   the type is in, by the rule it takes; a feature that is not drawn gets the row it would join.
 * @property {(rule: SpecificRuleInput) => boolean} drawsRow Whether the rules of a captioned rule's caption
 *   and color draw a row.
 * @property {(key: string) => string | null} otherRowType The drawn type whose "other" row (a row the drawing
 *   names, of its features of no rule) `key` is; null when `key` is no such row.
 */

/**
 * What a rule row's key reads of the Results it is drawn in: the features they draw and the type colors they draw
 * them in (`ruleFillsOf`). Every reader of a key passes it; a reader of no Result says so with no features.
 * @typedef {Pick<LegendRowContext, 'features' | 'palette'>} DrawnLegendState
 */

/**
 * What one Result of a drawing gives the keys of its rule rows: the features it draws, the type colors it draws them
 * in, and Python's rows of it (none before its first Generate, nor where a reader has none).
 * @typedef {Pick<LegendRowContext, 'features' | 'palette' | 'pythonRows'>} ResultLegendState
 */

// The key of a feature type's "other" row: Python's `_other_feature_legend_key`.
/** @param {string} type @returns {PythonLegendKey} */
export const otherFeatureRowKey = (type) => captionRowKey(type === 'CDS' ? 'other proteins' : `other ${type}s`);

// The type among `types` whose own row or "other" row `key` is.
/** @param {string} key @param {Iterable<string>} types @returns {string | undefined} */
const typeOfRow = (key, types) => [...types].find((type) => key === type || key === otherFeatureRowKey(type));
// The type whose row `key` is, among `types`; else the key itself, a key that
// names a type's row.
/** @param {PythonLegendKey} key @param {Iterable<string>} types @returns {string} */
export const legendRowType = (key, types) => typeOfRow(key, types) ?? key;

// THE feature rows of a Legend and the row each drawn feature is in, the twin
// of Python's `FeatureLegendRows` (gbdraw/legend/table.py), whose rows Generate
// draws: both answer tests/fixtures/legend_generated_rows.json (`feature_rows`).
// A rule for every type (`*`) is read as the same rule written for each type
// (OV-3460). The captioned rules the drawn features take, by the features'
// type, caption and color, are used (Python's `legend_rule_use`; OV-3465): a
// type whose features take a captioned rule draws its "other" row for its
// features of no rule, and no row for its features of an uncaptioned rule;
// any other type draws its own row (Legend contract 05A4-04 L-C2, L-C3). A
// rule of the type color value (`TYPE_COLOR_VALUE`) draws each type in its
// color (`ruleFill`); one without a caption draws its features as features of
// no rule (Python's `rule_use`).
/**
 * @param {readonly SpecificRuleInput[]} rules The rule table, in table order.
 * @param {Iterable<LegendRowUse>} drawn The drawn features.
 * @param {FeatureTypeColors} palette The type colors.
 * @returns {FeatureLegendRows}
 */
export const featureLegendRows = (rules, drawn, palette) => {
  const table = rules.map((rule) => normalizeSpecificRule(rule));
  /** @param {Pick<SpecificColorRule, 'cap' | 'color'>} rule */
  const use = (rule) => JSON.stringify([rule.cap, rule.color]);
  /** @param {string} type @param {Pick<SpecificColorRule, 'cap' | 'color'>} rule */
  const typeUse = (type, rule) => JSON.stringify([type, rule.cap, rule.color]);
  // A rule as the Legend reads a feature of the type that takes it: its
  // caption and the fill it draws the type in; null for an uncaptioned rule
  // of the type color value (Python's `rule_use`).
  /** @param {string} type @param {Pick<SpecificColorRule, 'cap' | 'color'>} rule @returns {Pick<SpecificColorRule, 'cap' | 'color'> | null} */
  const useOf = (type, rule) => (isTypeColorValue(rule.color) && !rule.cap
    ? null : { cap: rule.cap, color: normalizeColor(ruleFill(rule.color, type, palette)) });
  /** @type {Set<string>} The captioned rules the drawn features take, by `typeUse`. */
  const used = new Set();
  /** @type {Set<string>} */
  const defaultTypes = new Set();
  /** @type {Set<string>} */
  const present = new Set();
  for (const { type, winner } of drawn) {
    present.add(type);
    const rule = winner && useOf(type, normalizeSpecificRule(winner));
    if (!rule) defaultTypes.add(type);
    else if (rule.cap) used.add(typeUse(type, rule));
  }
  /** @param {string} type The captioned rules of a type, its own and those for every type, in the fill each draws it in. */
  const typeRules = (type) => table.filter((rule) => (rule.feat === type || rule.feat === '*') && rule.cap)
    .map((rule) => /** @type {Pick<SpecificColorRule, 'cap' | 'color'>} */ (useOf(type, rule)));
  /** @type {Map<string, boolean>} */
  const drawsByType = new Map();
  /** @param {string} type */
  const rulesDraw = (type) => {
    if (!drawsByType.has(type)) drawsByType.set(type, typeRules(type).some((rule) => used.has(typeUse(type, rule))));
    return /** @type {boolean} */ (drawsByType.get(type));
  };
  /** @param {string} type */
  const typeRow = (type) => {
    if (!table.some((rule) => rule.feat === type || rule.feat === '*') || !rulesDraw(type)) return captionRowKey(type);
    return defaultTypes.has(type) ? otherFeatureRowKey(type) : null;
  };
  return {
    present: [...present],
    rulesDraw,
    typeRow,
    rowOf: (type, winner) => {
      const rule = winner && useOf(type, normalizeSpecificRule(winner));
      if (rule?.cap) return rule;
      if (!rulesDraw(type)) return captionRowKey(type);
      return rule ? null : otherFeatureRowKey(type);
    },
    drawsRow: (rule) => {
      const normalized = normalizeSpecificRule(rule);
      return Boolean(normalized.cap) && [...present].some((type) => {
        const group = /** @type {Pick<SpecificColorRule, 'cap' | 'color'>} */ (useOf(type, normalized));
        return used.has(typeUse(type, group)) && typeRules(type).some((each) => use(each) === use(group));
      });
    },
    otherRowType: (key) => [...present].find((type) => key !== type && typeRow(type) === key) ?? null
  };
};

// Whether a key is the row of a default track: "<NT> content" and "<NT> skew",
// each also as "(+)" and "(-)", and "Depth" (gbdraw/legend/table.py::
// _generated_legend_rows; tests/fixtures/legend_generated_rows.json holds
// Python's keys). The track rows of the request that drew the Result give its
// fill (track-legend-rows.js, OV-484).
const TRACK_ROW = /^[ACGTU]{2} (content|skew)(?: \([+-]\))?$/;
/** @param {string} key */
export const isTrackLegendKey = (key) => key === 'Depth' || TRACK_ROW.test(key);

// The one owner of which rules the displayed Result's drawn features take as
// their first rule, whether every match is known (`ruleMatcher().firstIfKnown`),
// and the feature rows Generate draws for them (`featureLegendRows`).
/**
 * @param {Pick<LegendRowContext, 'rules' | 'features' | 'palette'>} context
 * @returns {DrawnRuleMatches}
 */
export const drawnRuleMatches = ({ rules = [], features = [], palette }) => {
  const winnerOf = ruleMatcher([...rules]).firstIfKnown;
  let known = true;
  /** @type {Map<Partial<SpecificColorRule>, LegendRowFeature[]>} */
  const takenBy = new Map();
  /** @type {LegendRowUse[]} */
  const drawn = [];
  for (const feature of features) {
    const rule = winnerOf(feature);
    if (rule === undefined) known = false;
    drawn.push({ type: normalizeText(feature.type), winner: rule ?? null });
    if (!rule || !normalizeText(rule.cap)) continue;
    const colored = takenBy.get(rule);
    if (colored) colored.push(feature);
    else takenBy.set(rule, [feature]);
  }
  return { known, takenBy, rows: featureLegendRows(rules, drawn, palette) };
};

/**
 * The type colors and the drawn types a rule of the type color value is read with (`ruleRowFills`).
 * @typedef {{ palette: FeatureTypeColors, types: readonly string[] }} RuleFills
 */

// The drawn types of a Legend row context, in the order they first appear.
/** @param {DrawnLegendState} context @returns {RuleFills} */
export const ruleFillsOf = ({ features = [], palette }) => ({
  palette, types: [...new Set(Array.from(features, (feature) => normalizeText(feature.type)))]
});

// A rule as the row it draws for a feature of `type` reads it, normalized:
// of the type color value, in the color of its own type, or, a rule for every
// type (`*`), of `type`; any other rule as written. A rule that comes resolved
// (a row of `ruleRowFills`) reads as it is. A rule for every type of the type
// color value draws a row per type (`typedRowFills`), so it has no row without
// `type`: reading it without one throws, rather than keying it in a color no
// row has (OV-3803).
/**
 * @param {SpecificRuleInput} rule
 * @param {RuleFills} fills
 * @param {unknown} [type] The type of the feature the row is for; required for a `*` rule of the type color value.
 * @returns {SpecificColorRule}
 */
const ruleInFills = (rule, fills, type) => {
  const normalized = normalizeSpecificRule(rule);
  if (!isTypeColorValue(normalized.color)) return normalized;
  const drawnType = normalized.feat === '*' ? normalizeText(type) : normalized.feat;
  if (normalized.feat === '*' && !drawnType) throw new Error('A color rule for every feature type in the type color is keyed per feature type: read its rows with captionedRuleRows.');
  return { ...normalized, color: normalizeColor(ruleFill(normalized.color, drawnType, fills.palette)) };
};
// The rows a rule can draw, each with the type it is drawn for (`*`: every
// type) and the rule in that row's fill (`ruleInFills`): one, or a rule for
// every type of the type color value once in the color of each drawn type
// (Python's `FeatureLegendRows.rule_fills`).
/**
 * @param {SpecificRuleInput} rule
 * @param {RuleFills} fills
 * @returns {{ type: string, row: SpecificColorRule }[]}
 */
const typedRowFills = (rule, fills) => {
  const normalized = normalizeSpecificRule(rule);
  if (normalized.feat !== '*' || !isTypeColorValue(normalized.color)) return [{ type: normalized.feat, row: ruleInFills(normalized, fills) }];
  return fills.types.map((type) => ({ type, row: ruleInFills(normalized, fills, type) }));
};
// Whether revising `rows` of `current` into `next` keeps each Legend row it
// touches whole: no other rule draws a row with a revised rule, before or
// after. Such a revision moves no row, since its rules keep their places
// (OV-158). One that takes part of a row, or meets another row's rules, can
// move a row: the row's first drawing rule changes (OV-2264). The rows are
// keyed in the type colors and drawn types of the shown drawing (`fills`).
/**
 * @param {readonly Partial<SpecificColorRule>[]} current
 * @param {readonly Partial<SpecificColorRule>[]} next
 * @param {readonly Partial<SpecificColorRule>[]} rows
 * @param {RuleFills} fills
 */
export const revisesWholeRuleRows = (current, next, rows, fills) => /** @type {const} */ ([
  [current, (/** @type {Partial<SpecificColorRule>} */ rule) => rows.includes(rule)],
  [next, (/** @type {Partial<SpecificColorRule>} */ rule) => !current.includes(rule)]
]).every(([table, revised]) => {
  const rowKey = createRuleLegendCaptions(table, fills);
  /** @param {Partial<SpecificColorRule>} rule */
  const keysOf = (rule) => ruleRowFills(rule, fills).map((row) => rowKey(row));
  const captioned = table.filter((rule) => String(rule.cap ?? '').trim());
  const touched = new Set(captioned.filter(revised).flatMap(keysOf));
  return !captioned.some((rule) => !revised(rule) && keysOf(rule).some((key) => touched.has(key)));
});

// The color a feature of `type` takes when it joins the row `rule` draws
// for it ("Use existing"): `shown`, the color of that row; `written`, the
// color its own rule is written in. A rule of the type color value keeps that
// value for a feature of its type (or of any type, a rule for every type), so
// the feature stays in the row after a palette change; a feature of another
// type takes the row's color, which it cannot follow. Null without a rule of a color.
/**
 * @param {SpecificRuleInput | null | undefined} rule
 * @param {unknown} type
 * @param {FeatureTypeColors} palette
 * @returns {{ shown: string, written: string } | null}
 */
export const joinedRuleColor = (rule, type, palette) => {
  const normalized = normalizeSpecificRule(rule ?? {});
  if (!normalized.color) return null;
  const shown = ruleInFills(normalized, { palette, types: [] }, type).color;
  const follows = isTypeColorValue(normalized.color) && (normalized.feat === '*' || normalized.feat === normalizeText(type));
  return { shown, written: follows ? TYPE_COLOR_VALUE : shown };
};
/** @param {SpecificRuleInput} rule @param {RuleFills} fills @returns {SpecificColorRule[]} */
export const ruleRowFills = (rule, fills) => typedRowFills(rule, fills).map(({ row }) => row);
/** @param {readonly SpecificRuleInput[]} rules @param {RuleFills} fills */
export const rulesInFills = (rules, fills) => rules.flatMap((rule) => ruleRowFills(rule, fills));
// The rows of `rules`, each as its key (`keyOf` of the rule in that row's
// fill) and the color it draws the row in (`ruleRowFills`).
/**
 * @param {readonly SpecificRuleInput[]} rules
 * @param {RuleFills} fills
 * @param {(rule: SpecificColorRule) => PythonLegendKey} keyOf
 * @returns {{ caption: PythonLegendKey, color: string }[]}
 */
export const ruleRowColors = (rules, fills, keyOf) => rulesInFills(rules, fills).map((row) => ({ caption: keyOf(row), color: row.color }));

// The Python rows N-06 does not compare a rule with, read from the displayed
// features and their known rule matches (`drawnRuleMatches`; a match not
// known yet colors no feature), only when a captioned rule may need them:
// - OV-306: Python draws no row under a feature type's name once a captioned
//   rule of that type, or of every type, colors a feature (`rulesDraw`).
// - OV-307: a row Python drew for a rule is no row of the renderer, also when
//   a live edit of the rule's color leaves the row in the color Python drew it
//   in until Generate: a feature Python drew in the row's color takes a rule
//   of the row's caption, and the caption names no row Python generates (a
//   track's row, or a displayed type's own or "other" row).
// `ruleMatches` is `drawnRuleMatches` of the same `context.rules` and
// `context.features`, when the caller read it already, so one compile reads
// the features once; it is an argument, not a context field, so a context
// copied with other rules never carries it.
/**
 * @param {LegendRowContext} context
 * @param {DrawnRuleMatches} [ruleMatches]
 * @returns {Set<PythonLegendKey>}
 */
const rowsLeftOut = (context, ruleMatches) => {
  const { rules = [], pythonRows, features = [] } = context;
  const captioned = rulesInFills(rules, ruleFillsOf(context)).filter((rule) => rule.cap);
  const typed = new Set(captioned.map((rule) => rule.feat).filter((feat) => feat === '*' || pythonRows.has(captionRowKey(feat))));
  /** @type {Map<string, string>} A caption naming a Python row of another color, and that row's fill. */
  const recolored = new Map();
  captioned.forEach((rule) => {
    const color = fillIdentity(pythonRows.get(captionRowKey(rule.cap))?.color);
    if (color && color !== fillIdentity(rule.color)) recolored.set(rule.cap, color);
  });
  if (typed.size === 0 && recolored.size === 0) return new Set();
  const { takenBy, rows } = ruleMatches ?? drawnRuleMatches({ rules, features, palette: context.palette });
  /** @type {Set<string>} */
  const drawnForRules = new Set();
  for (const [rule, colored] of takenBy) {
    const caption = normalizeText(rule.cap);
    for (const feature of colored) {
      if (recolored.get(caption) === fillIdentity(feature.fill_color)) drawnForRules.add(caption);
    }
  }
  /** @param {string} caption */
  const generated = (caption) => isTrackLegendKey(caption) || typeOfRow(caption, rows.present) !== undefined;
  const types = [...typed].flatMap((feat) => (feat === '*' ? rows.present : [feat]));
  return new Set([
    ...types.filter((type) => pythonRows.has(captionRowKey(type)) && rows.rulesDraw(type)),
    ...[...drawnForRules].filter((caption) => !generated(caption))
  ].map(captionRowKey));
};

// Python's rows of the displayed Result, by key and drawn fill: those N-06
// compares a rule with, and those the rules take (`rowsLeftOut`).
/**
 * @param {LegendRowContext} context
 * @param {DrawnRuleMatches} [ruleMatches] As for `rowsLeftOut`.
 */
const pythonRowsOf = (context, ruleMatches) => {
  const leftOut = rowsLeftOut(context, ruleMatches);
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
const rendererRowsOf = (rows, context) => {
  const { originalLegendOrder = [], rules = [] } = context;
  const generated = new Set(originalLegendOrder.map(normalizeText).filter(Boolean));
  const normalizedRules = rulesInFills(rules, ruleFillsOf(context)).filter((rule) => rule.cap);
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
 * A next rule table that Python answered with the Legend order Generate draws from it on the
 * displayed Result (`SynthesizedColorRules.legend`): a rule row, by the index in `rules` of a rule
 * that draws it and the feature type it is drawn for (`LegendOrderRuleRow`); a text is the key of the
 * renderer's own row of a feature type. Null when the build did not ask for it.
 * @typedef {{ rules: Partial<SpecificColorRule>[], legend: readonly (LegendOrderRuleRow | string)[] | null }} AnsweredRuleTable
 */

/**
 * How a rule table write (rule-actions.js) builds the next table from the one at its turn (OV-444);
 * a row returned as it is or as `revise(row, patch)` keeps its identity. Null writes nothing. A build
 * that Python answers (a group color edit) returns a promise; the write waits for it at its turn, and
 * the answer carries the Legend order of the features `drawnFeatures` gives (the features the
 * displayed Result draws, found by identity in the rule preparation's catalog).
 * @callback SpecificRuleBuild
 * @param {readonly SpecificColorRule[]} current
 * @param {(row: SpecificColorRule, patch: Partial<SpecificColorRule>) => SpecificColorRule} revise
 * @param {() => Iterable<object>} drawnFeatures
 * @param {() => LegendRowContext} legendRows The Legend state of the drawing whose table it builds (its displayed
 *   Result's rows, drawn features and type colors), by which it reads a row's rules and keys its rows.
 * @returns {Partial<SpecificColorRule>[] | AnsweredRuleTable | null
 *   | Promise<Partial<SpecificColorRule>[] | AnsweredRuleTable | null>}
 */

// The rules of one caption and one fill draw one Legend row: that pair is the
// identity of their Legend group, for the row allocation and for the groups below.
/** @param {Pick<SpecificColorRule, 'cap' | 'color'>} rule */
const legendGroupKey = (rule) => JSON.stringify([rule.cap, fillIdentity(rule.color)]);

// Whether `rule` is in the Legend group `group` names: the same caption and the
// same fill identity, a color name or spelling that is one fill being one fill.
/**
 * @param {Pick<SpecificColorRule, 'cap' | 'color'>} group
 * @param {SpecificRuleInput} rule
 */
export const inRuleGroup = (group, rule) => (
  legendGroupKey(normalizeSpecificRule(group)) === legendGroupKey(normalizeSpecificRule(rule))
);

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
// the target of a Legend rename read this one derivation. A group's color is
// the swatch its rows are drawn in: a rule of the type color value is in its
// type's color, and one for every type is in a group of each drawn type
// (`typedRowFills`), which takes features of that type only.
/**
 * @param {readonly SpecificRuleInput[]} rules
 * @param {RuleFills} fills
 * @param {string} [type]
 * @returns {ColorRuleLegendGroup[]}
 */
export const colorRuleLegendGroups = (rules, fills, type) => {
  /** @type {Map<string, { name: ColorRuleGroupName, color: SwatchColor, types: string[] }>} */
  const groups = new Map();
  const typed = type === undefined || fills.types.includes(type) ? fills : { ...fills, types: [...fills.types, type] };
  for (const { type: rowType, row } of rules.flatMap((rule) => typedRowFills(rule, typed))) {
    if (!row.cap) continue;
    const key = legendGroupKey(row);
    const group = groups.get(key)
      || { name: /** @type {ColorRuleGroupName} */ (row.cap), color: fillOf(row.color), types: [] };
    if (!group.types.includes(rowType)) group.types.push(rowType);
    groups.set(key, group);
  }
  return [...groups.values()].filter((group) => type === undefined || colorRuleGroupTakes(group, type));
};

/**
 * Returns (rule, type) -> the key of the Legend row Generate draws for that
 * rule, for a feature of `type`: THE key function of rule rows (Python's
 * `rule_legend_keys` over `FeatureLegendRows.rule_fills`; vector
 * tests/fixtures/legend_rule_keys.json). It reads every rule in the fill it
 * draws (`ruleInFills`, with `fills`), so a rule of the type color value is
 * keyed by its type's color; a rule for every type of that value draws a
 * row per drawn type, so it is keyed only with the type of the row's feature
 * (without one the read throws: `captionedRuleRows` lists every row). A rule
 * that comes resolved (`ruleRowFills`) is keyed as it is.
 * @param {readonly SpecificRuleInput[] | null | undefined} rules
 * @param {RuleFills} fills The type colors and the drawn types.
 * @param {RendererLegendRow[]} [rendererRows]
 * @returns {(rule: SpecificRuleInput, type?: unknown) => PythonLegendKey}
 */
export const createRuleLegendCaptions = (rules, fills, rendererRows = []) => {
  /** @type {Map<string, RendererLegendRow>} */
  const rows = new Map((rendererRows || []).map((row) => [row.caption, row]));
  const normalizedRules = rulesInFills(rules || [], fills).filter((rule) => rule.cap);
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
  return (rule, type) => {
    const normalized = ruleInFills(rule, fills, type);
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
 * @param {RuleFills} fills The type colors and drawn types the rows are drawn with (`ruleRowFills`).
 * @param {RendererLegendRow[]} [rendererRows]
 * @returns {Map<PythonLegendKey, PythonLegendKey>}
 */
export const rekeyedRuleRows = (before, after, fills, rendererRows = []) => {
  const kept = new Set(rulesInFills(after, fills).map(legendGroupKey));
  const groups = rulesInFills(before, fills).filter((group) => kept.has(legendGroupKey(group)));
  return changedRowKeys(groups.map((group) => [group, group]),
    createRuleLegendCaptions(before, fills, rendererRows), createRuleLegendCaptions(after, fills, rendererRows));
};

// The Legend rows a take-back of rules removes from a drawing (Delete
// settings: rows that name a feature instance of a turned-off record, which
// no Result draws), and the rows it re-keys, read from the keys the drawing's
// Legend stores (OV-1280). Python's fills and swatch lines of those rows are
// not read, and no key is allocated again from the colors its rows list: a Legend
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
 * @param {RuleFills} fills The type colors and drawn types the rows are drawn with (`ruleRowFills`).
 * @returns {{ gone: RuleRowColor[], rekeyed: Map<PythonLegendKey, PythonLegendKey> }}
 */
export const listedRuleRowRemoval = (before, after, listedKeys, fills) => {
  const listed = [...listedKeys].map(captionRowKey);
  /** @param {readonly Partial<SpecificColorRule>[]} rules @returns {Map<string, SpecificColorRule>} */
  const groupsOf = (rules) => new Map(rulesInFills(rules, fills)
    .filter((rule) => rule.cap).map((rule) => [legendGroupKey(rule), rule]));
  const removed = groupsOf(before);
  const kept = groupsOf(after);
  kept.forEach((_group, id) => removed.delete(id));
  const captions = new Set([...removed.values(), ...kept.values()].map((group) => group.cap));
  /** @param {SpecificColorRule} group @returns {PythonLegendKey} */
  const storedKey = (group) => captionRowKey(listed.find((key) => !captions.has(key) && isAllocatedRow(key, group)) ?? group.cap);
  const ruleRows = new Set([...removed.values(), ...kept.values()].map(storedKey));
  const keptFills = captionFills([...kept.values()]);
  /** @type {Map<PythonLegendKey, PythonLegendKey>} */
  const rekeyed = new Map();
  kept.forEach((group) => {
    const from = storedKey(group);
    if (from !== group.cap && keptFills.get(group.cap)?.size === 1 && !listed.some((key) => key === group.cap && !ruleRows.has(key))) {
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
// longer the text saved. `context.rules` are the rules as saved; each row a
// rule draws (`typedRowFills`: a rule for every type of the type color value,
// one per drawn type) is paired with the same row of its reading, which keeps
// the rule's type and color. A caption still written in several colors, or
// naming a row of another swatch, keeps its key.
/**
 * @param {LegendRowContext} context
 * @returns {Map<PythonLegendKey, PythonLegendKey>}
 */
export const savedRuleRowRekeys = (context) => {
  const saved = [...context.rules || []];
  const read = savedSpecificRules(saved);
  const fills = ruleFillsOf(context);
  /** @type {[SpecificColorRule, SpecificColorRule][]} */
  const pairs = saved.flatMap((rule, index) => {
    // The reading changes no rule's type and no rule for every type's color, so both draw the same rows.
    const readRows = typedRowFills(read[index], fills);
    return typedRowFills(rule, fills).map(({ row }, at) => /** @type {[SpecificColorRule, SpecificColorRule]} */ ([row, readRows[at].row]));
  });
  return changedRowKeys(pairs, ruleLegendCaptions({ ...context, rules: saved }), ruleLegendCaptions({ ...context, rules: read }));
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

// The key of the Legend row Generate draws for a rule of `context.rules`, for
// a feature of `type` (`createRuleLegendCaptions`, with the context's type
// colors and drawn types), read for many rules: the keys are allocated once,
// when a rule first needs it. A rule for every type of the type color value
// takes its row's type, as there.
/**
 * @param {LegendRowContext} context
 * @param {DrawnRuleMatches} [ruleMatches] As for `rowsLeftOut`.
 * @returns {(rule: Partial<SpecificColorRule> | null | undefined, type?: unknown) => PythonLegendKey}
 */
export const ruleLegendCaptions = (context, ruleMatches) => {
  const rows = pythonRowsOf(context, ruleMatches).compared;
  const fills = ruleFillsOf(context);
  /** @type {Map<string, Set<string>> | null} */
  let fillsOf = null;
  /** @type {((rule: SpecificRuleInput, type?: unknown) => PythonLegendKey) | null} */
  let allocated = null;
  return (rule, type) => {
    const normalized = ruleInFills(rule ?? {}, fills, type);
    fillsOf ||= captionFills(rulesInFills(context.rules || [], fills));
    if ((fillsOf.get(normalized.cap)?.size || 0) < 2 && !namesRowOfAnotherSwatch(normalized, rows)) return captionRowKey(normalized.cap);
    allocated ||= createRuleLegendCaptions(context.rules, fills, rendererRowsOf(rows, context));
    return allocated(normalized);
  };
};

// THE rows the captioned rules of `context.rules` draw: each rule once in
// every fill it draws (`typedRowFills`: a rule for every type of the type
// color value once per drawn type), with the key Generate draws the row under
// (`ruleLegendCaptions`) and its name, the rule's caption. A reader of a
// table's rule rows reads them here, never a rule's key without its type
// (OV-3803). `context.features` is read more than once.
/**
 * @param {LegendRowContext} context
 * @param {DrawnRuleMatches} [ruleMatches] As for `ruleLegendCaptions`.
 * @returns {{ key: PythonLegendKey, name: string }[]}
 */
export const captionedRuleRows = (context, ruleMatches) => {
  const captioned = (context.rules || []).filter((rule) => normalizeText(rule?.cap));
  if (captioned.length === 0) return [];
  const keyOf = ruleLegendCaptions(context, ruleMatches);
  const fills = ruleFillsOf(context);
  return captioned.flatMap((rule) => typedRowFills(rule, fills)).map(({ row }) => ({ key: keyOf(row), name: row.cap }));
};

// THE Legend states of the Results of a drawing that its rule rows are keyed
// in: the table (its rules and inventory) over each Result's own drawn state.
// Python keys the rule rows of each Result of a batch over the types that
// Result draws (`gbdraw circular` of several records): a `*` rule of the type
// color value in one type's color only where that Result draws no other type.
// A reader of a drawing's Results reads them one by one here, never as one
// state of every Result's features.
/**
 * @param {Pick<LegendRowContext, 'rules' | 'originalLegendOrder'>} table
 * @param {Iterable<ResultLegendState>} results
 * @returns {LegendRowContext[]}
 */
export const resultLegendRowContexts = (table, results) => Array.from(results, (result) => ({ ...table, ...result }));

// The keys of every rule row the Results of a drawing draw: the union of each
// Result's (`captionedRuleRows`), so a reader names the same rule rows
// whichever Result of a batch is shown.
/** @param {readonly LegendRowContext[]} contexts @returns {Set<PythonLegendKey>} */
export const resultsRuleRowKeys = (contexts) => new Set(contexts.flatMap((context) => captionedRuleRows(context).map((row) => row.key)));

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

/**
 * A rule's Legend row as a rule commit draws it: its key (`caption`), its name
 * (the rule's caption) where that is not the key, and its color.
 * @typedef {{ caption: PythonLegendKey, name?: string, color: string }} RuleLegendIntent
 */

// The Legend rows of a rule table: one per row key, for the rules `draws`
// admits (a commit: those a displayed feature takes). The keys are allocated
// from the whole table, as Python allocates them, whichever rules draw a row.
/**
 * @param {readonly Partial<SpecificColorRule>[] | null | undefined} rules
 * @param {RuleFills} fills The type colors and drawn types the rows are drawn with (`ruleRowFills`).
 * @param {RendererLegendRow[]} [rendererRows]
 * @param {(rule: Partial<SpecificColorRule>) => boolean} [draws]
 */
export const buildLegendIntents = (rules, fills, rendererRows = [], draws = () => true) => {
  const legendCaption = createRuleLegendCaptions(rules, fills, rendererRows);
  // An intent is a rule's Legend row: its key (`caption`), its color, and its
  // name (the rule's caption) where that is not the key, as `data-legend-name`.
  /** @type {Map<PythonLegendKey, RuleLegendIntent>} */
  const byCaption = new Map();
  for (const rule of rules || []) {
    if (!draws(rule)) continue;
    for (const normalized of ruleRowFills(rule, fills)) {
      const caption = legendCaption(normalized);
      if (caption && !byCaption.has(caption)) {
        byCaption.set(caption, { caption, color: normalized.color, ...(normalized.cap === caption ? {} : { name: normalized.cap }) });
      }
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
  const fills = ruleFillsOf(context);
  const legendCaption = createRuleLegendCaptions(rules, fills, rendererRowsOf(pythonRowsOf(context).compared, context));
  return rules.filter((rule) => ruleRowFills(rule, fills).some((each) => legendCaption(each) === target));
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

