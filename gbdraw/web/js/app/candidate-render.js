// @ts-check
/** @import { FeatureCatalogAdmission } from '../services/feature-catalog.js' */
/** @import { SvgAdmissionRuntime, SvgResultTransform } from '../services/svg-result-ingestion.js' */
import { resolveColorToHex } from './color-utils.js';
import { defaultLegendCaptionOrder, isLegendOrderEdited } from './legend/utils.js';
import { cloneJsonValue } from '../services/json-clone.js';
import { biologicalFeatureKey } from '../services/feature-catalog.js';
import {
  admitCurrentGeneratedResults
} from '../services/svg-result-ingestion.js';

/**
 * The committed editor state a direct mutation plan compiles. Python and the
 * editor owners define the row shapes (R7), so each map is a record here.
 * @typedef {object} EditorPlanOptions
 * @property {FeatureCatalogAdmission} catalogAdmission
 * @property {Record<string, any>} [featureColorOverrides]
 * @property {Record<string, any>} [featureStrokeOverrides]
 * @property {Record<string, any>} [featureOverrides]
 * @property {Record<string, any>[]} [legendEntries]
 * @property {Record<string, any>[]} [deletedLegendEntries]
 * @property {string[]} [originalLegendOrder]
 * @property {boolean} [sourceReplaced]
 * @property {Iterable<string>} [addedLegendCaptions] Captions the renderer drew without a manual rule.
 * @property {Iterable<string>} [unrequestedDepthCaptions] Captions of Depth series the request left out (Show Depth off); Python cannot name their rows (OV-81).
 * @property {Record<string, any>} [legendColorOverrides]
 * @property {Record<string, any>} [legendStrokeOverrides]
 * @property {Record<string, any>[]} [manualSpecificRules]
 * @property {string[] | null} [replayDefaultLegendOrder]
 * @property {SvgResultTransform[]} [resultTransforms] One transform per Result, by index.
 * @property {SvgResultTransform | null} [transformSvg] A transform applied to every Result.
 */

/**
 * @typedef {EditorPlanOptions & SvgAdmissionRuntime & { generationResponse: any }} CandidateCommitOptions
 */

const text = (value) => String(value ?? '').trim();
const hasOwn = (object, key) => Object.prototype.hasOwnProperty.call(object || {}, key);

const normalizePaint = (value, label) => {
  const raw = text(value);
  if (!raw) return '';
  const resolved = text(resolveColorToHex(raw));
  if (
    /^(?:none|transparent)$/i.test(resolved)
    || /^#(?:[0-9a-f]{3}|[0-9a-f]{4}|[0-9a-f]{6}|[0-9a-f]{8})$/i.test(resolved)
    || /^rgba?\(\s*[-+.\d%]+(?:\s*[,/]\s*|\s+)[-+.\d%]+(?:\s*[,/]\s*|\s+)[-+.\d%]+(?:\s*[,/]\s*[-+.\d%]+)?\s*\)$/i.test(resolved)
    || /^hsla?\(\s*[-+.\d]+(?:deg|grad|rad|turn)?(?:\s*[,/]\s*|\s+)[-+.\d%]+(?:\s*[,/]\s*|\s+)[-+.\d%]+(?:\s*[,/]\s*[-+.\d%]+)?\s*\)$/i.test(resolved)
  ) {
    return resolved;
  }
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

const resolvedStableTargets = (catalogAdmission, key) => (
  catalogAdmission.renderedTargetsByOverrideKey.get(key) || []
);

const renderedResultIndexes = (catalogAdmission, renderedId) => (
  catalogAdmission.resultIndexesByRenderedId.get(renderedId) || new Set()
);

// Per-feature edits name a source identity in one mode; each Result of that
// mode draws it with its own rendered ID (design Q4 6.2, R2).
const identityTargets = (catalogAdmission, row) => (row?.scope === catalogAdmission.mode
  ? resolvedStableTargets(catalogAdmission, biologicalFeatureKey(row?.recordKey, row?.biologicalFeatureId))
  : []);

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
  legendColorOverrides = {},
  legendStrokeOverrides = {},
  manualSpecificRules = [],
  replayDefaultLegendOrder = null,
  resultTransforms = [],
  transformSvg = null
}) => {
  if (!catalogAdmission || !Array.isArray(catalogAdmission.resultNames)) {
    throw new Error('Direct editor mutation planning requires an admitted feature catalog.');
  }
  const operationsByResult = catalogAdmission.resultNames.map(() => emptyOperations());
  const normalizedFeatureColorOverrides = {};
  const normalizedFeatureStrokeOverrides = {};

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
    normalizedFeatureStrokeOverrides[key] = {
      ...cloneJsonValue(rawOverride, {}),
      ...(strokeColor ? { strokeColor } : {}),
      ...(strokeWidth !== null ? { strokeWidth } : {})
    };
  });

  const hiddenRenderedIds = new Set();
  Object.values(featureOverrides || {}).forEach((row) => {
    const targets = identityTargets(catalogAdmission, row);
    const featureMode = text(row?.featureVisibility).toLowerCase();
    const labelMode = text(row?.labelVisibility).toLowerCase();
    targets.forEach(({ resultIndex, renderedId }) => {
      const operations = operationsByResult[resultIndex];
      if (!operations) return;
      if (featureMode === 'on' || featureMode === 'off') {
        operations.featureVisibility.push({ renderedId, mode: featureMode });
        if (featureMode === 'off') hiddenRenderedIds.add(renderedId);
      }
      if (typeof row?.labelText === 'string') {
        operations.labelText.push({ renderedId, value: row.labelText });
      }
      if (labelMode === 'on' || labelMode === 'off') {
        operations.labelVisibility.push({ renderedId, mode: labelMode });
      }
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
  const unrequestedDepth = new Set(Array.from(unrequestedDepthCaptions || []).map(text).filter(Boolean));
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
  const allResultIndexes = operationsByResult.map((_, index) => index);

  // D-08 (PD-OI-063): an edited Legend order is replayed over the renderer's
  // slots. The renderer places generated entries in their generated order and
  // direct additions after them; only a different order emits an operation,
  // and a renamed entry then takes its slot from that order. A displayed batch
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
        allowMissing: sourceReplaced || manualCaptions.has(entry.caption),
        xPos: legendOrderChanged ? null : entry.xPos,
        yPos: legendOrderChanged ? null : entry.yPos
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

  // Category style preferences outlive the current generated entry projection.
  // Apply a returning category's preference without synthesizing a manual row.
  const entriesByCaption = new Map(currentEntries.map(entry => [entry.caption, entry]));
  const styledCaptions = new Set([...Object.keys(legendColorOverrides), ...Object.keys(legendStrokeOverrides)]);
  styledCaptions.forEach(caption => {
    const entry = entriesByCaption.get(caption);
    const originalCaption = entry?.originalCaption || caption;
    if (deletedCaptions.has(originalCaption)) return;
    const isOriginal = originalCaptions.has(originalCaption);
    const targetCaption = isOriginal ? originalCaption : caption;
    const legendRenderedIds = entry && entry.featureIds.length > 0
      ? entry.featureIds : [...(renderedIdsByDirectCaption.get(caption) || [])];
    const allowMissing = !entry || (sourceReplaced && isOriginal) || unrequestedDepth.has(caption)
      || (rendererDerivedCaptions.has(caption)
      && (legendRenderedIds.length === 0 || legendRenderedIds.every(id => hiddenRenderedIds.has(id))));
    // Each Result styles only the category features it renders. A batch Result
    // that renders none of them draws no row for the category (OV-45).
    const renderedIdsIn = (resultIndex) => legendRenderedIds
      .filter(id => renderedResultIndexes(catalogAdmission, id).has(resultIndex));
    // A row of unknown features may be absent from a Result that does not draw it
    // or whose draft removed it; admission decides that from Python's Legend row
    // facts (`metadata.legendRows`, OV-46, OV-63), not here.
    const allowMissingIn = (resultIndex) => allowMissing
      || (legendRenderedIds.length > 0 && renderedIdsIn(resultIndex).length === 0);
    if (hasOwn(legendColorOverrides, caption)) {
      const color = normalizePaint(legendColorOverrides[caption], 'legend color');
      if (color) operationsByResult.forEach((operations, resultIndex) => {
        operations.legendFills.push({ caption: targetCaption, color, allowMissing: allowMissingIn(resultIndex) });
      });
    }
    const stroke = legendStrokeOverrides[caption];
    if (stroke && typeof stroke === 'object') {
      const strokeColor = hasOwn(stroke, 'strokeColor') ? normalizePaint(stroke.strokeColor, 'legend stroke color') : '';
      const strokeWidth = hasOwn(stroke, 'strokeWidth') ? normalizeStrokeWidth(stroke.strokeWidth) : null;
      if (strokeColor || strokeWidth !== null) operationsByResult.forEach((operations, resultIndex) => {
        operations.legendStrokes.push({
          caption: targetCaption, strokeColor, strokeWidth, allowMissing: allowMissingIn(resultIndex),
          renderedIds: renderedIdsIn(resultIndex)
        });
      });
    }
  });

  deletedCaptions.forEach((caption) => {
    // Explicit deletions remain valid while their category is absent, including
    // subsequent Generate calls after a source replacement.
    addToResults(operationsByResult, allResultIndexes, 'legendDeletes', { caption, allowMissing: true });
  });

  resultTransforms.forEach((transform, index) => {
    if (typeof transform === 'function') operationsByResult[index].callerTransforms.push(transform);
  });

  if (typeof transformSvg === 'function') {
    addToResults(operationsByResult, allResultIndexes, 'callerTransforms', transformSvg);
  }

  const frozenOperations = Object.freeze(operationsByResult.map(freezeOperations));
  const kind = frozenOperations.some((operations) => operationCount(operations) > 0)
    ? 'MUTATING'
    : 'EMPTY';
  return {
    plan: Object.freeze({ kind, operationsByResult: frozenOperations }),
    normalizedFeatureColorOverrides,
    normalizedFeatureStrokeOverrides
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
    featureStrokeOverrides: bundle.normalizedFeatureStrokeOverrides,
    mutationPlan: bundle.plan
  };
};

export const prepareReflowResultCommit = prepareCandidateRenderCommit;
