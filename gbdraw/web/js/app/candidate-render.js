// @ts-check
/** @import { FeatureCatalogAdmission } from '../services/feature-catalog.js' */
/** @import { LegendStrokeOperation, SvgAdmissionRuntime, SvgResultTransform } from '../services/svg-result-ingestion.js' */
/** @import { PythonLegendKey, PythonLegendRow, RenderedFeatureId, SwatchColor } from '../services/legend-svg.js' */
/** @import { SpecificColorRule } from '../services/specific-color-rules.js' */
import { normalizeDefaultColor, resolveColorToHex } from '../utils/color-utils.js';
import { defaultLegendCaptionOrder, isLegendOrderEdited } from '../services/legend-svg.js';
import { cloneJsonValue } from '../services/json-clone.js';
import { biologicalFeatureKey } from '../services/feature-catalog.js';
import { ruleMatcher } from '../services/rule-matchers.js';
import { ruleLegendCaptions } from '../services/specific-color-rules.js';
import { recordStructuralMetric } from '../services/runtime-test-hooks.js';
import { resolveFeatureDrawn } from '../services/feature-visibility.js';
import {
  admitCurrentGeneratedResults
} from '../services/svg-result-ingestion.js';

/**
 * The committed editor state a direct mutation plan compiles. Python and the
 * editor owners define the row shapes (R7), so each map is a record here.
 * @typedef {object} EditorPlanOptions
 * @property {EditorAddressing} catalogAdmission The feature catalog's admission, or the addressing of
 *   a Result without a catalog (`displayedFeatureAddressing`).
 * @property {Record<string, any>} [featureColorOverrides]
 * @property {Record<string, any>} [featureStrokeOverrides]
 * @property {Record<string, any>} [featureOverrides]
 * @property {Record<string, any>[]} [legendEntries]
 * @property {Record<string, any>[]} [deletedLegendEntries]
 * @property {string[]} [originalLegendOrder]
 * @property {boolean} [sourceReplaced]
 * @property {Iterable<string>} [addedLegendCaptions] Captions the renderer drew without a manual rule.
 * @property {Iterable<string>} [unrequestedDepthCaptions] Captions of Depth series the request left out (Show Depth off); Python cannot name their rows (OV-81).
 * @property {Record<string, any>[]} [dormantLegendEntries] OV-120: renamed rows an earlier Generate did not draw. Their
 *   renames and styles apply where this Result draws the row again.
 * @property {Record<string, any>} [legendColorOverrides]
 * @property {Record<string, any>} [legendStrokeOverrides]
 * @property {Record<string, any>[]} [manualSpecificRules]
 * @property {string[] | null} [replayDefaultLegendOrder]
 * @property {SvgResultTransform[]} [resultTransforms] One transform per Result, by index.
 * @property {SvgResultTransform | null} [transformSvg] A transform applied to every Result.
 * @property {LivePreview | null} [livePreview] The draft's palette, specific-color rules, and Feature
 *   visibility as the displayed Result previews them. Generate leaves it out: Python drew them.
 */

/**
 * What Python draws from the draft at Generate, previewed live on the
 * displayed Result for the operation domains it shows (`domains`): the
 * applied palette, the specific-color rules through their prepared matches
 * (Python's, R4), and whether each rendered feature is drawn
 * (`resolveFeatureDrawn`). The compile runs only the stages those domains
 * need (`COMPILE_STAGES`).
 * @typedef {object} LivePreview
 * @property {readonly string[]} domains
 * @property {Record<string, string>} paletteColors
 * @property {ReturnType<typeof import('../services/feature-visibility.js').featureDrawnContext> | null} drawnContext
 *   `featureDrawnContext` of the drawing; read when `domains` has featureVisibility.
 * @property {RuleLegendRows | null} [ruleRows] The Legend rows a rule commit shows at once.
 * @property {ReadonlyMap<PythonLegendKey, PythonLegendRow>} [pythonRows] Python's Legend rows of the
 *   displayed Result (`pythonLegendRows`), by which the draft allocates a rule's row (N-06).
 */

/**
 * The Legend rows of a rule commit until Python draws them (`prepareFileLegendEntries`):
 * the rows its rules add, each before the row it replaces, and the rows they retire.
 * @typedef {{ add: Array<{ caption: string, color: string, before?: string }>, retire: string[] }} RuleLegendRows
 */

/**
 * What the compile reads of a feature catalog's admission.
 * @typedef {FeatureCatalogAdmission | ReturnType<typeof import('../services/feature-override-identity.js').displayedFeatureAddressing>} EditorAddressing
 */

/**
 * @typedef {EditorPlanOptions & SvgAdmissionRuntime & { generationResponse: any, catalogAdmission: FeatureCatalogAdmission }} CandidateCommitOptions
 */

const text = (value) => String(value ?? '').trim();
const hasOwn = (object, key) => Object.prototype.hasOwnProperty.call(object || {}, key);

// A committed paint is read with the Default colors domain (D-41).
const normalizePaint = (value, label) => {
  const raw = text(value);
  if (!raw) return '';
  const color = normalizeDefaultColor(raw);
  if (color) return color;
  throw new Error(`Invalid ${label} override in the committed editor state.`);
};

const normalizeStrokeWidth = (value) => {
  if (value === null || value === undefined || value === '') return null;
  const width = Number(value);
  if (!Number.isFinite(width) || width < 0) {
    throw new Error('Invalid stroke width override in the committed editor state.');
  }
  return width;
};

/**
 * The planner's mutable form of one Result's operations (the shape
 * services/svg-result-ingestion.js declares as `SvgMutationOperations`).
 * @typedef {Record<
 *   'featureFills' | 'featureStrokes' | 'featureVisibility' | 'labelText' | 'labelVisibility'
 *   | 'legendFills' | 'legendRenames' | 'legendDeletes' | 'legendAdds' | 'legendOrder',
 *   Record<string, any>[]
 * > & { legendStrokes: LegendStrokeOperation[], callerTransforms: SvgResultTransform[] }} MutableOperations
 */

/** @returns {MutableOperations} */
const emptyOperations = () => ({
  featureFills: [],
  featureStrokes: [],
  featureVisibility: [],
  labelText: [],
  labelVisibility: [],
  legendFills: [],
  legendStrokes: [],
  legendRenames: [],
  legendDeletes: [],
  legendAdds: [],
  legendOrder: [],
  callerTransforms: []
});

const operationCount = (operations) => Object.values(operations)
  .reduce((total, entries) => total + entries.length, 0);

const freezeOperations = (operations) => {
  Object.values(operations).forEach((entries) => {
    entries.forEach((entry) => {
      if (entry && typeof entry === 'object') Object.freeze(entry);
    });
    Object.freeze(entries);
  });
  return Object.freeze(operations);
};

const addToResults = (operationsByResult, resultIndexes, domain, operation) => {
  resultIndexes.forEach((resultIndex) => {
    const target = operationsByResult[resultIndex];
    if (target) target[domain].push(
      typeof operation === 'function' ? operation : { ...operation }
    );
  });
};

const matchingRuleDerivedFill = (override, manualSpecificRules) => {
  if (!override || typeof override !== 'object') return false;
  const caption = text(override.caption);
  const color = text(override.color).toLowerCase();
  if (!caption || !color) return false;
  return (Array.isArray(manualSpecificRules) ? manualSpecificRules : []).some((rule) => (
    text(rule?.cap) === caption && text(rule?.color).toLowerCase() === color
  ));
};

const samePaint = (left, right) => (
  text(resolveColorToHex(text(left))).toLowerCase() === text(resolveColorToHex(text(right))).toLowerCase()
);

// The palette color of a Legend row Python draws for a feature type or a
// track (`other <type>s` rows take their type's color).
/** @type {Record<string, string>} */
const LEGEND_PALETTE_KEYS = {
  'GC content': 'gc_content',
  'GC skew (+)': 'skew_high',
  'GC skew (-)': 'skew_low'
};
/** @param {string} caption @param {Record<string, string>} palette */
const paletteLegendColor = (caption, palette) => {
  const key = LEGEND_PALETTE_KEYS[caption] || caption;
  if (palette[key]) return palette[key];
  const lower = caption.toLowerCase();
  if (lower === 'other proteins') return palette.CDS || '';
  if (!lower.startsWith('other ')) return '';
  let type = caption.slice(6).trim();
  if (type.toLowerCase() === 'proteins') return palette.CDS || '';
  if (type.endsWith('s')) type = type.slice(0, -1);
  return palette[type] || '';
};

// The color the draft gives the features of a Legend row (its Python key)
// where a displayed Result predates the draft: its first rule's color, else
// its palette color, as Python draws them. A rule's row is the key Python
// draws for it (N-06: "<caption> [<hex>]" when the caption names a row of
// another color), allocated from the colors Python drew for its rows, never
// from their swatches (R15-4), nor a row its caption only spells. The live
// compile shows it on rows without a Legend color of their own, and a row's
// stroke reaches by it (`LegendRowReach.draftColor`); the feature popup reads
// it the same way.
/**
 * @param {{ rules: readonly Partial<SpecificColorRule>[], pythonRows?: ReadonlyMap<PythonLegendKey, PythonLegendRow>,
 *   originalLegendOrder: readonly string[], paletteColors: Record<string, string> }} draft
 * @returns {{ ruleCaptions: PythonLegendKey[], colorOf: (key: PythonLegendKey) => SwatchColor | null }}
 */
export const draftLegendRowColors = ({ rules, pythonRows = new Map(), originalLegendOrder, paletteColors }) => {
  const rowOf = ruleLegendCaptions({ rules: [...rules], pythonRows, originalLegendOrder: [...originalLegendOrder] });
  /** @type {Map<PythonLegendKey, string>} */
  const ruleColors = new Map();
  rules.forEach((rule) => {
    const key = /** @type {PythonLegendKey} */ (text(rule?.cap) ? text(rowOf(rule)) : '');
    if (key && !ruleColors.has(key)) ruleColors.set(key, text(rule.color));
  });
  return {
    ruleCaptions: [...ruleColors.keys()],
    colorOf: (key) => {
      const color = normalizePaint(ruleColors.get(key) || paletteLegendColor(key, paletteColors), 'legend color');
      return color ? /** @type {SwatchColor} */ (color) : null;
    }
  };
};

// The draft's Legend rows: the listed and the dormant entries, Python's keys
// of the generated and the deleted rows, and the row a styled caption (a key
// of `legendColorOverrides` or `legendStrokeOverrides`) addresses, by Python's
// key where Python draws it; null for a deleted row. The compile and the
// feature popup read the rows through it, so a deleted row's stroke reaches
// no feature in either.
/**
 * @param {Pick<EditorPlanOptions, 'legendEntries' | 'deletedLegendEntries' | 'dormantLegendEntries' | 'originalLegendOrder'>} legend
 */
export const draftLegendRows = ({
  legendEntries = [], deletedLegendEntries = [], dormantLegendEntries = [], originalLegendOrder = []
}) => {
  const currentEntries = normalizedLegendEntries(legendEntries);
  const originalCaptions = new Set(
    (Array.isArray(originalLegendOrder) ? originalLegendOrder : []).map(text).filter(Boolean)
  );
  const deletedCaptions = new Set(
    (Array.isArray(deletedLegendEntries) ? deletedLegendEntries : [])
      .map((entry) => text(entry?.originalCaption || entry?.caption))
      .filter((caption) => originalCaptions.has(caption))
  );
  // OV-120: a renamed row an earlier Generate hid (GC off, Show Depth off)
  // waits in the drawing; this Result may draw it again or still leave it out.
  const shownOriginals = new Set(currentEntries.map((entry) => entry.originalCaption));
  const dormantEntries = normalizedLegendEntries(dormantLegendEntries)
    .filter((entry) => entry.caption !== entry.originalCaption && !shownOriginals.has(entry.originalCaption)
      && !deletedCaptions.has(entry.originalCaption));
  const dormantOriginals = new Set(dormantEntries.map((entry) => entry.originalCaption));
  const entriesByCaption = new Map([...dormantEntries, ...currentEntries].map(entry => [entry.caption, entry]));
  /** @param {string} caption */
  const styledRow = (caption) => {
    const entry = entriesByCaption.get(caption);
    const originalCaption = entry?.originalCaption || caption;
    if (deletedCaptions.has(originalCaption)) return null;
    const dormant = dormantOriginals.has(originalCaption);
    const isOriginal = dormant || originalCaptions.has(originalCaption);
    return { entry, dormant, isOriginal, targetCaption: /** @type {PythonLegendKey} */ (isOriginal ? originalCaption : caption) };
  };
  return { currentEntries, originalCaptions, deletedCaptions, dormantEntries, styledRow };
};

// OV-276: the swatch color of each listed Legend panel row that the palette
// colors, by the row's caption: a row of Python's without a Legend color of
// its own that no rule draws takes the color the live compile's
// `legendFills` stage gives its Python key (`draftLegendRowColors`). The
// other rows keep their listed color. Nothing writes it into the rows, which
// are inputs of the rule preparation.
/**
 * @param {Pick<EditorPlanOptions, 'legendEntries' | 'deletedLegendEntries' | 'dormantLegendEntries' | 'originalLegendOrder'>
 *   & { legendColorOverrides: Readonly<Record<string, unknown>>, rules: readonly Partial<SpecificColorRule>[],
 *   pythonRows: ReadonlyMap<PythonLegendKey, PythonLegendRow>, paletteColors: Record<string, string> }} draft
 * @returns {Map<string, SwatchColor>}
 */
export const draftLegendPanelColors = ({ legendColorOverrides, rules, pythonRows, paletteColors, ...legend }) => {
  const rows = draftLegendRows(legend);
  const draft = draftLegendRowColors({ rules, pythonRows, originalLegendOrder: [...rows.originalCaptions], paletteColors });
  const ruleRows = new Set(draft.ruleCaptions);
  /** @type {Map<string, SwatchColor>} */
  const colors = new Map();
  rows.currentEntries.forEach(({ caption }) => {
    const row = hasOwn(legendColorOverrides, caption) ? null : rows.styledRow(caption);
    if (!row || !rows.originalCaptions.has(row.targetCaption) || ruleRows.has(row.targetCaption)) return;
    const color = draft.colorOf(row.targetCaption);
    if (color) colors.set(caption, color);
  });
  return colors;
};

// Python's fill of a feature type (gbdraw/features/colors.py
// `default_color_map.get(feature.type, "#d3d3d3")`): the palette has no
// fallback key.
/** @param {Record<string, string>} palette @param {unknown} type */
const paletteFeatureColor = (palette, type) => {
  const key = String(type ?? '');
  return (hasOwn(palette, key) && palette[key]) || '#d3d3d3';
};

// The stages of one compile and the operation domains each one writes. Generate
// runs every stage but `rules`; a live display runs the stages of the domains
// it shows. `rules` matches the specific-color rules for the live feature fills.
const COMPILE_STAGES = Object.freeze({
  fills: Object.freeze(['featureFills']),
  rules: Object.freeze(['featureFills']),
  visibility: Object.freeze(['featureVisibility']),
  labels: Object.freeze(['labelText', 'labelVisibility']),
  legend: Object.freeze(['legendRenames', 'legendDeletes', 'legendAdds', 'legendOrder']),
  legendFills: Object.freeze(['legendFills']),
  strokes: Object.freeze(['featureStrokes', 'legendStrokes'])
});

// The operation domains each kind of live edit shows on the displayed Result
// (app-setup): a Legend row stroke also strokes the row's features; the
// palette and the rules fill features and Legend rows; a Result display shows
// the Legend structure and the paint domains whose intent changed. An edit
// that returns Legend rows (Restore, a History step of the rows) shows their
// structure and fills: a returning row of Python's shows the color Generate
// gives it.
export const LIVE_EDIT_DOMAINS = Object.freeze({
  strokes: COMPILE_STAGES.strokes,
  legendFills: COMPILE_STAGES.legendFills,
  fills: Object.freeze(['featureFills', 'legendFills']),
  visibility: COMPILE_STAGES.visibility,
  legendStructure: COMPILE_STAGES.legend,
  legendRows: Object.freeze([...COMPILE_STAGES.legend, ...COMPILE_STAGES.legendFills]),
  // A deleted row strokes no feature, and a returned row strokes its own
  // again, as Generate draws them (OV-293).
  deletedRows: Object.freeze([...COMPILE_STAGES.legend, ...COMPILE_STAGES.legendFills, ...COMPILE_STAGES.strokes])
});
/** @type {ReadonlyArray<[string[], readonly string[]]>} */
const EDITOR_PAINT_PATHS = Object.freeze([
  [['editorState', 'featureStrokes'], LIVE_EDIT_DOMAINS.strokes],
  [['editorState', 'legend', 'strokeOverrides'], LIVE_EDIT_DOMAINS.strokes],
  [['editorState', 'legend', 'colorOverrides'], LIVE_EDIT_DOMAINS.legendFills],
  ...['entries', 'dormantEntries', 'addedCaptions']
    .map((/** @type {string} */ key) => /** @type {[string[], readonly string[]]} */ ([['editorState', 'legend', key], LIVE_EDIT_DOMAINS.legendRows])),
  [['editorState', 'legend', 'deletedEntries'], LIVE_EDIT_DOMAINS.deletedRows]
]);
// The domains an Undo or Redo shows, by the paths of the step's changes.
/** @param {unknown} changes A History step's change list. @returns {string[]} */
export const editorPaintDomains = (changes) => [...new Set((Array.isArray(changes) ? changes : []).flatMap(({ path } = {}) => (
  Array.isArray(path)
    ? EDITOR_PAINT_PATHS.filter(([paintPath]) => paintPath.every((key, index) => index >= path.length || path[index] === key))
      .flatMap(([, domains]) => domains)
    : []
)))];

const resolvedStableTargets = (catalogAdmission, key) => (
  catalogAdmission.renderedTargetsByOverrideKey.get(key) || []
);

const renderedResultIndexes = (catalogAdmission, renderedId) => (
  catalogAdmission.resultIndexesByRenderedId.get(renderedId) || new Set()
);

// Per-feature edits name a source identity in the drawing of the catalog's
// mode; each Result of that mode draws it with its own rendered ID (design Q4
// 6.2, R2).
const identityTargets = (catalogAdmission, row) => (
  resolvedStableTargets(catalogAdmission, biologicalFeatureKey(row?.recordKey, row?.biologicalFeatureId))
);

const normalizedLegendEntries = (entries) => (
  Array.isArray(entries)
    ? entries.map((entry) => ({
        caption: text(entry?.caption),
        originalCaption: text(entry?.originalCaption || entry?.caption),
        color: text(entry?.color),
        xPos: Number.isFinite(Number(entry?.xPos)) ? Number(entry.xPos) : null,
        yPos: Number.isFinite(Number(entry?.yPos)) ? Number(entry.yPos) : null,
        featureIds: Object.freeze([
          ...new Set((Array.isArray(entry?.featureIds) ? entry.featureIds : []).map(text).filter(Boolean))
        ])
      })).filter((entry) => entry.caption)
    : []
);

/** @param {EditorPlanOptions} options */
const compilePlanBundle = ({
  catalogAdmission,
  featureColorOverrides = {},
  featureStrokeOverrides = {},
  featureOverrides = {},
  legendEntries = [],
  deletedLegendEntries = [],
  originalLegendOrder = [],
  sourceReplaced = false,
  addedLegendCaptions = [],
  unrequestedDepthCaptions = [],
  dormantLegendEntries = [],
  legendColorOverrides = {},
  legendStrokeOverrides = {},
  manualSpecificRules = [],
  replayDefaultLegendOrder = null,
  resultTransforms = [],
  transformSvg = null,
  livePreview = null
}) => {
  if (!catalogAdmission || !Array.isArray(catalogAdmission.resultNames)) {
    throw new Error('Direct editor mutation planning requires an admitted feature catalog.');
  }
  const operationsByResult = catalogAdmission.resultNames.map(() => emptyOperations());
  const normalizedFeatureColorOverrides = {};
  const drawnFillById = 'fillByRenderedId' in catalogAdmission ? catalogAdmission.fillByRenderedId : null;
  // The fill Python drew a rendered feature with: its catalog row's, else the
  // displayed Result's (a Result without a catalog).
  /** @param {string} renderedId @param {Record<string, any> | undefined} feature */
  const pythonFill = (renderedId, feature) => text(feature?.fill_color) || text(drawnFillById?.get(renderedId));
  const shownDomains = new Set(livePreview?.domains || []);
  /** @param {string} domain */
  const shows = (domain) => !livePreview || shownDomains.has(domain);
  // Every stage runs through `runStage`, which records it; one structural
  // metric per compile names the stages that ran (the work guard reads it).
  /** @type {Set<keyof typeof COMPILE_STAGES>} */
  const stagesRun = new Set();
  /**
   * @template T
   * @param {keyof typeof COMPILE_STAGES} name
   * @param {() => T} run
   * @returns {T | undefined}
   */
  const runStage = (name, run) => {
    if (name !== 'rules' && !COMPILE_STAGES[name].some(shows)) return undefined;
    stagesRun.add(name);
    return run();
  };
  const rules = Array.isArray(manualSpecificRules) ? manualSpecificRules : [];

  // Fill precedence, live as at Generate: a feature's fill edit, else its
  // first matching rule, else its palette color, else what Python drew. Python
  // drew the rules and the palette of the request, so only the live preview
  // emits them, where they differ from Python's fill. A feature whose rule
  // match is not known yet keeps the fill the Result shows (`color: null`).
  runStage('fills', () => {
    Object.entries(featureColorOverrides || {}).forEach(([key, rawOverride]) => {
      const color = normalizePaint(
        rawOverride && typeof rawOverride === 'object' && hasOwn(rawOverride, 'color')
          ? rawOverride.color
          : rawOverride,
        'feature fill'
      );
      const targets = resolvedStableTargets(catalogAdmission, key);
      if (!color || targets.length === 0) return;
      normalizedFeatureColorOverrides[key] = rawOverride && typeof rawOverride === 'object'
        ? { ...cloneJsonValue(rawOverride, {}), color }
        : color;
      if (matchingRuleDerivedFill(rawOverride, manualSpecificRules)) return;
      targets.forEach(({ resultIndex, renderedId }) => {
        operationsByResult[resultIndex].featureFills.push({ renderedId, color });
      });
    });
    if (!livePreview) return;
    const paletteColors = livePreview.paletteColors || {};
    const winnerOf = rules.length > 0 ? ruleMatcher(rules).firstIfKnown : null;
    operationsByResult.forEach((operations, resultIndex) => {
      const edited = new Set(operations.featureFills.map(({ renderedId }) => renderedId));
      /** @type {Map<string, Record<string, any>>} */
      const rendered = catalogAdmission.renderedFeaturesByResult?.[resultIndex] || new Map();
      const unedited = [...rendered].filter(([renderedId]) => !edited.has(renderedId));
      const winners = winnerOf
        ? runStage('rules', () => new Map(unedited.map(([renderedId, feature]) => [renderedId, winnerOf(feature)])))
        : null;
      unedited.forEach(([renderedId, feature]) => {
        const rule = winners ? winners.get(renderedId) : null;
        if (rule === undefined) {
          operations.featureFills.push({ renderedId, color: null });
          return;
        }
        const color = normalizePaint(rule ? rule.color : paletteFeatureColor(paletteColors, feature.type), 'feature fill');
        if (color && !samePaint(color, pythonFill(renderedId, feature))) operations.featureFills.push({ renderedId, color });
      });
    });
  });

  // Python draws every rendered feature, so the live preview names the
  // rendered features the resolver hides; Generate replays the On and Off
  // edits of the identity rows.
  const hiddenRenderedIds = new Set();
  runStage('visibility', () => {
    if (livePreview) {
      const { drawnContext } = livePreview;
      if (!drawnContext) return;
      operationsByResult.forEach((operations, resultIndex) => {
        /** @type {Map<string, Record<string, any>>} */
        const rendered = catalogAdmission.renderedFeaturesByResult?.[resultIndex] || new Map();
        rendered.forEach((feature, renderedId) => {
          if (resolveFeatureDrawn(feature, drawnContext) === false) operations.featureVisibility.push({ renderedId, mode: 'off' });
        });
      });
      return;
    }
    Object.values(featureOverrides || {}).forEach((row) => {
      const featureMode = text(row?.featureVisibility).toLowerCase();
      if (featureMode !== 'on' && featureMode !== 'off') return;
      identityTargets(catalogAdmission, row).forEach(({ resultIndex, renderedId }) => {
        const operations = operationsByResult[resultIndex];
        if (!operations) return;
        operations.featureVisibility.push({ renderedId, mode: featureMode });
        if (featureMode === 'off') hiddenRenderedIds.add(renderedId);
      });
    });
  });
  // The preview binder shows labels live; Generate replays them.
  runStage('labels', () => {
    Object.values(featureOverrides || {}).forEach((row) => {
      const labelMode = text(row?.labelVisibility).toLowerCase();
      identityTargets(catalogAdmission, row).forEach(({ resultIndex, renderedId }) => {
        const operations = operationsByResult[resultIndex];
        if (!operations) return;
        if (typeof row?.labelText === 'string') {
          operations.labelText.push({ renderedId, value: row.labelText });
        }
        if (labelMode === 'on' || labelMode === 'off') {
          operations.labelVisibility.push({ renderedId, mode: labelMode });
        }
      });
    });
  });

  const legendRows = draftLegendRows({ legendEntries, deletedLegendEntries, dormantLegendEntries, originalLegendOrder });
  const { currentEntries, originalCaptions, deletedCaptions, dormantEntries } = legendRows;
  const manualCaptions = new Set(
    (Array.isArray(manualSpecificRules) ? manualSpecificRules : [])
      .map((rule) => text(rule?.cap))
      .filter(Boolean)
  );
  const rendererDerivedCaptions = new Set(
    Array.from(addedLegendCaptions || []).map(text).filter(Boolean)
  );
  // A Depth row the draft keeps but Show Depth hides is excused like a switched-off
  // GC row. Python reports the GC rows; it is sent no Depth source, so the draft says.
  // The rule reads the caption an operation addresses in Python's output, so a
  // rename of the row and the styles under its new name are excused too (OV-88).
  const unrequestedDepth = new Set(Array.from(unrequestedDepthCaptions || []).map(text).filter(Boolean));
  const allResultIndexes = operationsByResult.map((_, index) => index);

  runStage('legend', () => {
    // D-08 (PD-OI-063): an edited Legend order is replayed over the renderer's
    // slots. The renderer places generated entries in their generated order and
    // direct additions after them; only a different order emits an operation.
    // A renamed entry keeps its row, and the Legend layout places it (OV-156,
    // zero shift), so a rename carries no position. A displayed batch
    // Result that may still show an earlier edited order also receives its own
    // default order, the generated order of that Result (D-07, B20, OV-47).
    const legendOrderChanged = isLegendOrderEdited(currentEntries, [...originalCaptions]);
    const defaultOrder = Array.isArray(replayDefaultLegendOrder) ? replayDefaultLegendOrder.map(text).filter(Boolean) : [];
    const replayedLegendCaptions = legendOrderChanged
      ? currentEntries.map((entry) => entry.caption)
      : (defaultOrder.length > 0 ? defaultLegendCaptionOrder(currentEntries, defaultOrder) : null);
    if (replayedLegendCaptions) {
      addToResults(operationsByResult, allResultIndexes, 'legendOrder', {
        captions: Object.freeze(replayedLegendCaptions)
      });
    }

    currentEntries.forEach((entry) => {
      const isOriginal = originalCaptions.has(entry.originalCaption);
      if (
        isOriginal
        && entry.caption !== entry.originalCaption
      ) {
        // A rename onto a color rule's caption is explicit and wins (PV-02). A
        // rule that replaced the source row leaves nothing to rename.
        addToResults(operationsByResult, allResultIndexes, 'legendRenames', {
          from: entry.originalCaption,
          to: entry.caption,
          allowMissing: sourceReplaced || manualCaptions.has(entry.caption) || unrequestedDepth.has(entry.originalCaption)
        });
      }
      if (
        originalCaptions.size > 0
        && !isOriginal
        && !manualCaptions.has(entry.caption)
        && !rendererDerivedCaptions.has(entry.caption)
      ) {
        const color = normalizePaint(entry.color, 'added legend color');
        if (color) {
          addToResults(operationsByResult, allResultIndexes, 'legendAdds', {
            caption: entry.caption,
            color,
            xPos: entry.xPos,
            yPos: entry.yPos
          });
        }
      }
    });

    dormantEntries.forEach((entry) => {
      addToResults(operationsByResult, allResultIndexes, 'legendRenames', {
        from: entry.originalCaption, to: entry.caption, allowMissing: true
      });
    });

    deletedCaptions.forEach((caption) => {
      // Explicit deletions remain valid while their category is absent, including
      // subsequent Generate calls after a source replacement.
      addToResults(operationsByResult, allResultIndexes, 'legendDeletes', { caption, allowMissing: true });
    });
    // The rows of a rule commit, shown at once on the displayed Result until
    // Python draws them (gaps 1 and 2): the rule's row where the Result lacks
    // it, before the row it replaces, and the retired row of a removed or
    // renamed rule. Generate compiles neither: Python draws the rules' rows.
    const ruleRows = livePreview?.ruleRows;
    (ruleRows?.add || []).forEach(({ caption, color, before = '' }) => {
      const paint = normalizePaint(color, 'rule legend color');
      if (caption && paint) addToResults(operationsByResult, allResultIndexes, 'legendAdds', {
        caption, color: paint, xPos: null, yPos: null, ifAbsent: true, before
      });
    });
    (ruleRows?.retire || []).forEach((caption) => {
      addToResults(operationsByResult, allResultIndexes, 'legendDeletes', { caption, allowMissing: true, retire: true });
    });
  });

  // Category style preferences outlive the current generated entry projection.
  // Apply a returning category's preference without synthesizing a manual row.
  const renderedIdsByDirectCaption = new Map();
  Object.entries(featureColorOverrides || {}).forEach(([key, override]) => {
    const caption = text(override?.caption);
    if (!caption) return;
    const renderedIds = renderedIdsByDirectCaption.get(caption) || new Set();
    resolvedStableTargets(catalogAdmission, key).forEach(({ renderedId }) => {
      if (renderedId) renderedIds.add(renderedId);
    });
    renderedIdsByDirectCaption.set(caption, renderedIds);
  });
  // The row a styled caption addresses in Python's output, and where a Result
  // may lack it; null for a deleted category.
  /** @param {string} caption */
  const styledRow = (caption) => {
    const row = legendRows.styledRow(caption);
    if (!row) return null;
    const { entry, dormant, isOriginal, targetCaption } = row;
    const namedIds = [...(renderedIdsByDirectCaption.get(caption) || [])];
    const legendRenderedIds = entry && entry.featureIds.length > 0 ? entry.featureIds : namedIds;
    const allowMissing = !entry || (sourceReplaced && isOriginal) || unrequestedDepth.has(targetCaption)
      || dormant || (rendererDerivedCaptions.has(caption)
      && (legendRenderedIds.length === 0 || legendRenderedIds.every(id => hiddenRenderedIds.has(id))));
    // Each Result styles only the category features it renders. A batch Result
    // that renders none of them draws no row for the category (OV-45).
    // A row of unknown features may be absent from a Result that does not draw it
    // or whose draft removed it; admission decides that from Python's Legend row
    // facts (`metadata.legendRows`, OV-46, OV-63), not here.
    /** @param {number} resultIndex */
    const allowMissingIn = (resultIndex) => allowMissing
      || (legendRenderedIds.length > 0
        && legendRenderedIds.every((id) => !renderedResultIndexes(catalogAdmission, id).has(resultIndex)));
    return { entry, targetCaption, namedIds, allowMissingIn };
  };
  const styledCaptions = new Set([...Object.keys(legendColorOverrides), ...Object.keys(legendStrokeOverrides)]);
  // Live only: at Generate Python's row has the color (OV-288).
  const draftRowColors = livePreview && draftLegendRowColors({
    rules, pythonRows: livePreview.pythonRows, originalLegendOrder: [...originalCaptions], paletteColors: livePreview.paletteColors
  });
  /** @param {PythonLegendKey} key */
  const draftRowColor = (key) => (draftRowColors ? draftRowColors.colorOf(key) : null);

  runStage('legendFills', () => {
    /** @type {Set<string>} */
    const filledCaptions = new Set();
    styledCaptions.forEach((caption) => {
      if (!hasOwn(legendColorOverrides, caption)) return;
      const row = styledRow(caption);
      if (!row) return;
      const color = normalizePaint(legendColorOverrides[caption], 'legend color');
      if (!color) return;
      filledCaptions.add(row.targetCaption);
      operationsByResult.forEach((operations, resultIndex) => {
        operations.legendFills.push({ caption: row.targetCaption, color, allowMissing: row.allowMissingIn(resultIndex) });
      });
    });
    // A row without a Legend color of its own shows the color the draft gives
    // its features (OV-146: every swatch fill is an operation, so a reconcile
    // never returns a rule or palette row to an older fill). A Result may not
    // draw the row.
    if (!draftRowColors) return;
    new Set([.../** @type {Set<PythonLegendKey>} */ (originalCaptions), ...draftRowColors.ruleCaptions]).forEach((caption) => {
      if (filledCaptions.has(caption) || deletedCaptions.has(caption)) return;
      const color = draftRowColor(caption);
      if (color) addToResults(operationsByResult, allResultIndexes, 'legendFills', { caption, color, allowMissing: true });
    });
  });

  runStage('strokes', () => {
    Object.entries(featureStrokeOverrides || {}).forEach(([key, rawOverride]) => {
      if (!rawOverride || typeof rawOverride !== 'object') return;
      const strokeColor = hasOwn(rawOverride, 'strokeColor')
        ? normalizePaint(rawOverride.strokeColor, 'feature stroke color')
        : '';
      const strokeWidth = hasOwn(rawOverride, 'strokeWidth')
        ? normalizeStrokeWidth(rawOverride.strokeWidth)
        : null;
      const targets = resolvedStableTargets(catalogAdmission, key);
      if ((!strokeColor && strokeWidth === null) || targets.length === 0) return;
      targets.forEach(({ resultIndex, renderedId }) => {
        operationsByResult[resultIndex].featureStrokes.push({
          renderedId,
          strokeColor,
          strokeWidth
        });
      });
    });
    // A Legend row's stroke reaches the features the executor finds in each
    // Result from Python's row there (`legendRowFeatureIds`, OV-123, OV-288);
    // the compile gives what the editor intent says of the row.
    /** @type {Map<number, RenderedFeatureId[]>} */
    const ownStrokesByResult = new Map();
    /** @param {number} resultIndex */
    const ownStrokeIdsIn = (resultIndex) => {
      const known = ownStrokesByResult.get(resultIndex);
      if (known) return known;
      const own = operationsByResult[resultIndex].featureStrokes.map(({ renderedId }) => /** @type {RenderedFeatureId} */ (renderedId));
      ownStrokesByResult.set(resultIndex, own);
      return own;
    };
    styledCaptions.forEach((caption) => {
      const stroke = legendStrokeOverrides[caption];
      if (!stroke || typeof stroke !== 'object') return;
      const row = styledRow(caption);
      if (!row) return;
      const strokeColor = hasOwn(stroke, 'strokeColor') ? normalizePaint(stroke.strokeColor, 'legend stroke color') : '';
      const strokeWidth = hasOwn(stroke, 'strokeWidth') ? normalizeStrokeWidth(stroke.strokeWidth) : null;
      if (!strokeColor && strokeWidth === null) return;
      const listedIds = /** @type {RenderedFeatureId[]} */ (row.entry ? row.entry.featureIds : []);
      const namedIds = /** @type {RenderedFeatureId[]} */ (row.namedIds);
      const draftColor = draftRowColor(row.targetCaption);
      operationsByResult.forEach((operations, resultIndex) => {
        operations.legendStrokes.push({
          caption: row.targetCaption, strokeColor, strokeWidth,
          allowMissing: row.allowMissingIn(resultIndex),
          reach: { listedIds, namedIds, ownStrokeIds: ownStrokeIdsIn(resultIndex), draftColor }
        });
      });
    });
  });

  resultTransforms.forEach((transform, index) => {
    if (typeof transform === 'function') operationsByResult[index].callerTransforms.push(transform);
  });

  if (typeof transformSvg === 'function') {
    addToResults(operationsByResult, allResultIndexes, 'callerTransforms', transformSvg);
  }

  recordStructuralMetric('editorPlanCompile', 1, {
    stages: [...stagesRun],
    domains: livePreview ? [...shownDomains] : null
  });
  const frozenOperations = Object.freeze(operationsByResult.map(freezeOperations));
  const kind = frozenOperations.some((operations) => operationCount(operations) > 0)
    ? 'MUTATING'
    : 'EMPTY';
  return {
    plan: Object.freeze({ kind, operationsByResult: frozenOperations }),
    normalizedFeatureColorOverrides
  };
};

/**
 * Compile direct live-editor deltas without enumerating admitted Features.
 * @param {EditorPlanOptions} [options]
 */
export const compileDirectEditorMutationPlan = (options = /** @type {EditorPlanOptions} */ ({})) => (
  compilePlanBundle(options).plan
);

/** @param {CandidateCommitOptions} options */
export const prepareCandidateRenderCommit = ({
  generationResponse,
  catalogAdmission,
  selectedFeatureTypes = null,
  sanitizer = globalThis.DOMPurify || globalThis.window?.DOMPurify,
  parser = globalThis.DOMParser || globalThis.window?.DOMParser,
  ...editorState
}) => {
  const bundle = compilePlanBundle({ catalogAdmission, ...editorState });
  return {
    results: admitCurrentGeneratedResults(generationResponse, {
      catalogAdmission,
      mutationPlan: bundle.plan,
      sanitizer,
      parser,
      selectedFeatureTypes
    }),
    featureState: catalogAdmission.featureState,
    featureColorOverrides: bundle.normalizedFeatureColorOverrides,
    mutationPlan: bundle.plan
  };
};

export const prepareReflowResultCommit = prepareCandidateRenderCommit;
