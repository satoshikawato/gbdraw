// @ts-check
/** @import { FeatureCatalogAdmission } from '../services/feature-catalog.js' */
/** @import { SvgAdmissionRuntime, SvgResultTransform } from '../services/svg-result-ingestion.js' */
import { normalizeDefaultColor } from '../utils/color-utils.js';
import { defaultLegendCaptionOrder, isLegendOrderEdited, legendRowFeatureIds } from '../services/legend-svg.js';
import { cloneJsonValue } from '../services/json-clone.js';
import { biologicalFeatureKey } from '../services/feature-catalog.js';
import { ruleMatcher } from '../services/rule-matchers.js';
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
 * @property {(() => Map<string, string>) | null} [shownFills] The fill each rendered feature of the
 *   displayed Result shows: what a Legend row stroke reaches by color when the compile does not
 *   preview the feature fills.
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
 *   | 'legendFills' | 'legendStrokes' | 'legendRenames' | 'legendDeletes' | 'legendAdds' | 'legendOrder',
 *   Record<string, any>[]
 * > & { callerTransforms: SvgResultTransform[] }} MutableOperations
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
  fills: ['featureFills'],
  rules: ['featureFills'],
  visibility: ['featureVisibility'],
  labels: ['labelText', 'labelVisibility'],
  legend: ['legendRenames', 'legendDeletes', 'legendAdds', 'legendOrder'],
  legendFills: ['legendFills'],
  strokes: ['featureStrokes', 'legendStrokes']
});

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

  const currentEntries = normalizedLegendEntries(legendEntries);
  const originalCaptions = new Set(
    (Array.isArray(originalLegendOrder) ? originalLegendOrder : []).map(text).filter(Boolean)
  );
  const deletedCaptions = new Set(
    (Array.isArray(deletedLegendEntries) ? deletedLegendEntries : [])
      .map((entry) => text(entry?.originalCaption || entry?.caption))
      .filter((caption) => originalCaptions.has(caption))
  );
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
  // OV-120: a renamed row an earlier Generate hid (GC off, Show Depth off)
  // waits in the drawing; this Result may draw it again or still leave it out.
  const shownOriginals = new Set(currentEntries.map((entry) => entry.originalCaption));
  const dormantEntries = normalizedLegendEntries(dormantLegendEntries)
    .filter((entry) => entry.caption !== entry.originalCaption && !shownOriginals.has(entry.originalCaption)
      && !deletedCaptions.has(entry.originalCaption));
  const dormantOriginals = new Set(dormantEntries.map((entry) => entry.originalCaption));
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
  });

  // Category style preferences outlive the current generated entry projection.
  // Apply a returning category's preference without synthesizing a manual row.
  const entriesByCaption = new Map([...dormantEntries, ...currentEntries].map(entry => [entry.caption, entry]));
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
  // Live, a row renamed in the Legend editor shows under its new key on the
  // displayed Result, which the executor tries when Python's key is absent.
  /** @param {string} targetCaption @param {string | undefined} caption */
  const shownRow = (targetCaption, caption) => (
    livePreview && caption && caption !== targetCaption ? { renamedCaption: caption } : {}
  );
  // The row a styled caption addresses in Python's output, and where a Result
  // may lack it; null for a deleted category.
  /** @param {string} caption */
  const styledRow = (caption) => {
    const entry = entriesByCaption.get(caption);
    const originalCaption = entry?.originalCaption || caption;
    if (deletedCaptions.has(originalCaption)) return null;
    const dormant = dormantOriginals.has(originalCaption);
    const isOriginal = dormant || originalCaptions.has(originalCaption);
    const targetCaption = isOriginal ? originalCaption : caption;
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
    return { entry, targetCaption, namedIds, allowMissingIn, shown: shownRow(targetCaption, entry?.caption) };
  };
  const styledCaptions = new Set([...Object.keys(legendColorOverrides), ...Object.keys(legendStrokeOverrides)]);

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
        operations.legendFills.push({ caption: row.targetCaption, ...row.shown, color, allowMissing: row.allowMissingIn(resultIndex) });
      });
    });
    // A row without a Legend color of its own shows its rule's color, else its
    // palette color, as Python draws them (OV-146: every swatch fill is an
    // operation, so a reconcile never returns a rule or palette row to an older
    // fill). A Result may not draw the row.
    if (!livePreview) return;
    /** @type {Map<string, string>} */
    const ruleColors = new Map();
    rules.forEach((rule) => {
      const caption = text(rule?.cap);
      if (caption && !ruleColors.has(caption)) ruleColors.set(caption, text(rule.color));
    });
    const renamedFrom = new Map([...dormantEntries, ...currentEntries].map((entry) => [entry.originalCaption, entry.caption]));
    new Set([...originalCaptions, ...ruleColors.keys()]).forEach((caption) => {
      if (filledCaptions.has(caption) || deletedCaptions.has(caption)) return;
      const color = normalizePaint(
        ruleColors.get(caption) || paletteLegendColor(caption, livePreview.paletteColors || {}), 'legend color'
      );
      if (color) addToResults(operationsByResult, allResultIndexes, 'legendFills', {
        caption, ...shownRow(caption, renamedFrom.get(caption)), color, allowMissing: true
      });
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
    // The fill each feature of a Result is drawn with and the features with a
    // stroke edit of their own: what a Legend row's stroke reads to reach the
    // row's features, as the live stroke does on the mounted Result
    // (`legendRowFeatureIds`, OV-123). Live, a feature without a fill
    // operation shows the fill the displayed Result shows, unless this compile
    // previews the fills (then Python's fill is the one it keeps).
    const previewsFills = Boolean(livePreview) && shownDomains.has('featureFills');
    /** @type {Map<string, string> | null} */
    let shown = null;
    /** @param {string} renderedId */
    const shownFill = (renderedId) => {
      shown ??= livePreview?.shownFills?.() || new Map();
      return text(shown.get(renderedId));
    };
    /** @type {Map<number, { drawnFills: Array<[string, string]>, ownStrokeIds: string[] }>} */
    const drawnFeaturesByResult = new Map();
    /** @param {number} resultIndex */
    const drawnFeaturesIn = (resultIndex) => {
      const known = drawnFeaturesByResult.get(resultIndex);
      if (known) return known;
      const operations = operationsByResult[resultIndex];
      const edited = new Map(operations.featureFills.map(({ renderedId, color }) => [renderedId, color]));
      /** @type {Map<string, Record<string, any>>} */
      const rendered = catalogAdmission.renderedFeaturesByResult?.[resultIndex] || new Map();
      /** @param {string} renderedId @param {Record<string, any>} feature */
      const drawnFill = (renderedId, feature) => {
        const own = edited.get(renderedId);
        if (own) return own;
        const keeps = livePreview && (own === null || !previewsFills);
        return (keeps ? shownFill(renderedId) : '') || pythonFill(renderedId, feature);
      };
      const drawn = {
        /** @type {Array<[string, string]>} */
        drawnFills: [...rendered].map(([renderedId, feature]) => [renderedId, drawnFill(renderedId, feature)]),
        ownStrokeIds: operations.featureStrokes.map(({ renderedId }) => renderedId)
      };
      drawnFeaturesByResult.set(resultIndex, drawn);
      return drawn;
    };
    styledCaptions.forEach((caption) => {
      const stroke = legendStrokeOverrides[caption];
      if (!stroke || typeof stroke !== 'object') return;
      const row = styledRow(caption);
      if (!row) return;
      const strokeColor = hasOwn(stroke, 'strokeColor') ? normalizePaint(stroke.strokeColor, 'legend stroke color') : '';
      const strokeWidth = hasOwn(stroke, 'strokeWidth') ? normalizeStrokeWidth(stroke.strokeWidth) : null;
      if (strokeColor || strokeWidth !== null) operationsByResult.forEach((operations, resultIndex) => {
        operations.legendStrokes.push({
          caption: row.targetCaption, ...row.shown, strokeColor, strokeWidth, allowMissing: row.allowMissingIn(resultIndex),
          renderedIds: legendRowFeatureIds(row.entry, { ...drawnFeaturesIn(resultIndex), namedIds: row.namedIds })
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
