import { resolveColorToHex } from './color-utils.js';
import { cloneJsonValue } from '../services/json-clone.js';
import {
  admitCurrentGeneratedResults
} from '../services/svg-result-ingestion.js';

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

const compilePlanBundle = ({
  catalogAdmission,
  featureColorOverrides = {},
  featureStrokeOverrides = {},
  featureVisibilityOverrides = {},
  labelTextFeatureOverrides = {},
  labelVisibilityOverrides = {},
  legendEntries = [],
  deletedLegendEntries = [],
  originalLegendOrder = [],
  addedLegendCaptions = [],
  legendColorOverrides = {},
  legendStrokeOverrides = {},
  manualSpecificRules = [],
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

  Object.entries(featureVisibilityOverrides || {}).forEach(([renderedId, rawMode]) => {
    const mode = text(rawMode).toLowerCase();
    const indexes = renderedResultIndexes(catalogAdmission, renderedId);
    if ((mode !== 'on' && mode !== 'off') || indexes.size === 0) return;
    addToResults(operationsByResult, indexes, 'featureVisibility', { renderedId, mode });
  });

  Object.entries(labelTextFeatureOverrides || {}).forEach(([renderedId, value]) => {
    const indexes = renderedResultIndexes(catalogAdmission, renderedId);
    if (indexes.size === 0) return;
    addToResults(operationsByResult, indexes, 'labelText', {
      renderedId,
      value: String(value ?? '')
    });
  });

  Object.entries(labelVisibilityOverrides || {}).forEach(([renderedId, rawMode]) => {
    const mode = text(rawMode).toLowerCase();
    const indexes = renderedResultIndexes(catalogAdmission, renderedId);
    if ((mode !== 'on' && mode !== 'off') || indexes.size === 0) return;
    addToResults(operationsByResult, indexes, 'labelVisibility', { renderedId, mode });
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
  const hiddenRenderedIds = new Set(
    Object.entries(featureVisibilityOverrides || {})
      .filter(([, mode]) => text(mode).toLowerCase() === 'off')
      .map(([renderedId]) => text(renderedId))
      .filter(Boolean)
  );
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

  currentEntries.forEach((entry) => {
    const isOriginal = originalCaptions.has(entry.originalCaption);
    if (
      isOriginal
      && entry.caption !== entry.originalCaption
      && !manualCaptions.has(entry.caption)
    ) {
      addToResults(operationsByResult, allResultIndexes, 'legendRenames', {
        from: entry.originalCaption,
        to: entry.caption,
        xPos: entry.xPos,
        yPos: entry.yPos
      });
    }
    const targetCaption = isOriginal ? entry.originalCaption : entry.caption;
    const legendRenderedIds = entry.featureIds.length > 0
      ? entry.featureIds
      : [...(renderedIdsByDirectCaption.get(entry.caption) || [])];
    const allowMissing = rendererDerivedCaptions.has(entry.caption)
      && (
        legendRenderedIds.length === 0
        || legendRenderedIds.every((renderedId) => hiddenRenderedIds.has(renderedId))
      );
    if (
      hasOwn(legendColorOverrides, entry.caption)
      && !deletedCaptions.has(entry.originalCaption)
    ) {
      const color = normalizePaint(legendColorOverrides[entry.caption], 'legend color');
      if (color) {
        addToResults(operationsByResult, allResultIndexes, 'legendFills', {
          caption: targetCaption,
          color,
          allowMissing
        });
      }
    }
    const stroke = legendStrokeOverrides?.[entry.caption];
    if (stroke && typeof stroke === 'object' && !deletedCaptions.has(entry.originalCaption)) {
      const strokeColor = hasOwn(stroke, 'strokeColor')
        ? normalizePaint(stroke.strokeColor, 'legend stroke color')
        : '';
      const strokeWidth = hasOwn(stroke, 'strokeWidth')
        ? normalizeStrokeWidth(stroke.strokeWidth)
        : null;
      if (strokeColor || strokeWidth !== null) {
        addToResults(operationsByResult, allResultIndexes, 'legendStrokes', {
          caption: targetCaption,
          strokeColor,
          strokeWidth,
          allowMissing,
          renderedIds: legendRenderedIds.filter((renderedId) => (
            renderedResultIndexes(catalogAdmission, renderedId).size > 0
          ))
        });
      }
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

  deletedCaptions.forEach((caption) => {
    addToResults(operationsByResult, allResultIndexes, 'legendDeletes', { caption });
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

/** Compile direct live-editor deltas without enumerating admitted Features. */
export const compileDirectEditorMutationPlan = (options = {}) => (
  compilePlanBundle(options).plan
);

export const prepareCandidateRenderCommit = ({
  generationResponse,
  catalogAdmission,
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
      parser
    }),
    featureState: catalogAdmission.featureState,
    featureColorOverrides: bundle.normalizedFeatureColorOverrides,
    featureStrokeOverrides: bundle.normalizedFeatureStrokeOverrides,
    mutationPlan: bundle.plan
  };
};

export const prepareReflowResultCommit = prepareCandidateRenderCommit;
