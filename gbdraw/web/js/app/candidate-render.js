import { resolveColorToHex } from './color-utils.js';
import {
  FEATURE_SELECTOR,
  filterFeatureFillTargets,
  getFeatureIdentity
} from './feature-dom.js';
import {
  biologicalFeatureKey,
  featureStateFromCatalog
} from '../services/feature-catalog.js';
import {
  migrateLegacyFeatureOverrides
} from '../services/feature-override-identity.js';
import { cloneJsonValue } from '../services/json-clone.js';
import { ingestSvgResults } from '../services/svg-result-ingestion.js';
import { PAIRWISE_LEGEND_SELECTOR } from './legend/utils.js';
import { applyStrokeOverridesToSvg } from './legend/stroke-actions.js';
import { applyLegendOrderToSvg } from './legend/manual-intent-command.js';
import { resolveOrderedSpecificRule } from './feature-utils.js';

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
    throw new Error('Invalid feature stroke width override in the committed editor state.');
  }
  return width;
};

const featureElementsById = (svg) => {
  const index = new Map();
  Array.from(svg.querySelectorAll(FEATURE_SELECTOR)).forEach((element) => {
    const id = getFeatureIdentity(element);
    if (!id) return;
    if (!index.has(id)) index.set(id, []);
    index.get(id).push(element);
  });
  return index;
};

const biologicalFeatureIndex = (item) => new Map(
  (Array.isArray(item?.biologicalFeatures) ? item.biologicalFeatures : []).map((feature) => [
    biologicalFeatureKey(feature?.recordKey, feature?.biologicalFeatureId),
    feature
  ])
);

const observeItemFeatureFills = ({
  svg,
  item,
  observedFillOverrides,
  specificRules,
  strokeOverrides
}) => {
  let changed = false;
  const elementsById = featureElementsById(svg);
  const biologicalByKey = biologicalFeatureIndex(item);
  item.features.forEach((rendered) => {
    const renderedId = text(rendered.svgId);
    const stableKey = biologicalFeatureKey(
      rendered.recordKey,
      rendered.biologicalFeatureId
    );
    const elements = elementsById.get(renderedId) || [];
    if (elements.length === 0) {
      throw new Error('Sanitized SVG content is missing a rendered feature binding.');
    }
    const biological = biologicalByKey.get(stableKey);
    if (!biological) {
      throw new Error('Feature catalogue is missing a biological Feature binding.');
    }
    const fillTargets = filterFeatureFillTargets(elements);
    if (fillTargets.length === 0) {
      throw new Error('Sanitized SVG content is missing a Feature fill target.');
    }
    const observedFills = new Set(fillTargets.map((element) => (
      normalizePaint(element.getAttribute('fill'), 'rendered feature fill').toLowerCase()
    )));
    if (observedFills.size !== 1 || observedFills.has('')) {
      throw new Error('Rendered Feature parts disagree on their committed fill.');
    }
    const [observedFill] = observedFills;
    const catalogFill = normalizePaint(rendered.fillColor, 'catalogued feature fill').toLowerCase();
    if (catalogFill && observedFill !== catalogFill) {
      throw new Error('Rendered Feature fill conflicts with admitted Feature metadata.');
    }
    const effectiveRule = resolveOrderedSpecificRule(biological, specificRules)?.rule || null;
    const caption = text(effectiveRule?.cap) || text(biological.type);
    const existingFill = observedFillOverrides[stableKey];
    if (
      existingFill
      && (
        text(existingFill.color).toLowerCase() !== observedFill
        || text(existingFill.caption) !== caption
      )
    ) {
      throw new Error('Rendered Feature copies disagree on their committed fill intent.');
    }
    observedFillOverrides[stableKey] = { color: observedFill, caption };

    const strokeOverride = strokeOverrides[stableKey];
    if (!strokeOverride || typeof strokeOverride !== 'object') return;
    const strokeValue = hasOwn(strokeOverride, 'strokeColor')
      ? normalizePaint(strokeOverride.strokeColor, 'feature stroke color')
      : '';
    const strokeWidth = hasOwn(strokeOverride, 'strokeWidth')
      ? normalizeStrokeWidth(strokeOverride.strokeWidth)
      : null;
    if (!strokeValue && strokeWidth === null) return;
    elements.forEach((element) => {
      if (strokeValue && element.getAttribute('stroke') !== strokeValue) {
        element.setAttribute('stroke', strokeValue);
        changed = true;
      }
      if (
        strokeWidth !== null
        && element.getAttribute('stroke-width') !== String(strokeWidth)
      ) {
        element.setAttribute('stroke-width', strokeWidth);
        changed = true;
      }
    });
  });
  return changed;
};

const removePairwiseIdentityLegend = (svg) => {
  let changed = false;
  svg.querySelectorAll(PAIRWISE_LEGEND_SELECTOR).forEach((legend) => {
    legend.remove();
    changed = true;
  });
  return changed;
};

export const prepareReflowResultCommit = ({
  results,
  suppressPairwiseIdentityLegend = false,
  features = [],
  featureStrokeOverrides = {},
  legendStrokeOverrides = {},
  legendOrderIntent = [],
  sanitizer = globalThis.DOMPurify || globalThis.window?.DOMPurify
}) => ({
  results: ingestSvgResults(results, {
    sanitizer,
    transformSvg: (svg) => {
      const legendChanged = suppressPairwiseIdentityLegend
        ? removePairwiseIdentityLegend(svg)
        : false;
      const strokeChanged = applyStrokeOverridesToSvg({
        svg,
        features,
        featureStrokeOverrides,
        legendStrokeOverrides
      }) > 0;
      const orderChanged = applyLegendOrderToSvg(svg, legendOrderIntent).movedEntries > 0;
      return legendChanged || strokeChanged || orderChanged;
    }
  })
});

export const prepareCandidateRenderCommit = ({
  results,
  catalog,
  mode = '',
  featureColorOverrides: derivedFillCache = {},
  featureStrokeOverrides,
  legendStrokeOverrides = {},
  legendOrderIntent = [],
  manualSpecificRules = [],
  legacyFeatures = [],
  preparedFeatureState = null,
  suppressPairwiseIdentityLegend = false,
  transformSvg = null,
  sanitizer = globalThis.DOMPurify || globalThis.window?.DOMPurify
}) => {
  const featureState = preparedFeatureState || featureStateFromCatalog(catalog, { mode });
  // Feature fills in a generated or saved Result are renderer-owned. The Web
  // cache is rebuilt from admitted SVG bytes below; it must never replay over
  // Python output or turn inherited palette/rule values into local intent.
  void derivedFillCache;
  const fillOverrides = {};
  const strokeOverrides = cloneJsonValue(featureStrokeOverrides, {});
  const migrationOptions = { legacyFeatures };

  migrateLegacyFeatureOverrides(
    strokeOverrides,
    featureState.extractedFeatures,
    migrationOptions
  );

  const itemsByIndex = new Map(
    catalog.items.map((item) => [item.resultIndex, item])
  );
  const processedResults = ingestSvgResults(results, {
    sanitizer,
    transformSvg: (svg, { resultIndex }) => {
      const item = itemsByIndex.get(resultIndex);
      if (!item) {
        throw new Error('The diagram engine returned incomplete feature metadata.');
      }
      const legendChanged = suppressPairwiseIdentityLegend
        ? removePairwiseIdentityLegend(svg)
        : false;
      const overridesChanged = observeItemFeatureFills({
        svg,
        item,
        observedFillOverrides: fillOverrides,
        specificRules: manualSpecificRules,
        strokeOverrides
      });
      const legendStrokeChanged = applyStrokeOverridesToSvg({
        svg,
        features: featureState.extractedFeatures,
        legendStrokeOverrides,
        featureStrokeOverrides: {}
      }) > 0;
      const callerChanged = typeof transformSvg === 'function'
        ? Boolean(transformSvg(svg, { resultIndex }))
        : false;
      const orderChanged = applyLegendOrderToSvg(svg, legendOrderIntent).movedEntries > 0;
      return legendChanged || overridesChanged || legendStrokeChanged || callerChanged || orderChanged;
    }
  });

  return {
    results: processedResults,
    featureState,
    featureColorOverrides: fillOverrides,
    featureStrokeOverrides: strokeOverrides
  };
};
