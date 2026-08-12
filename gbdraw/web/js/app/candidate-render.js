import { resolveColorToHex } from './color-utils.js';
import {
  FEATURE_SELECTOR,
  filterFeatureFillTargets,
  getFeatureIdentity
} from './feature-dom.js';
import { ruleMatchesFeature } from './feature-utils.js';
import {
  biologicalFeatureKey,
  featureStateFromCatalog
} from '../services/feature-catalog.js';
import {
  featureOverrideKey,
  migrateLegacyFeatureOverrides
} from '../services/feature-override-identity.js';
import { cloneJsonValue } from '../services/json-clone.js';
import { ingestSvgResults } from '../services/svg-result-ingestion.js';
import { PAIRWISE_LEGEND_SELECTOR } from './legend/utils.js';
import { applyStrokeOverridesToSvg } from './legend/stroke-actions.js';

const text = (value) => String(value ?? '').trim();
const hasOwn = (object, key) => Object.prototype.hasOwnProperty.call(object || {}, key);

const replaceRuleDerivedFillOverrides = (overrides, features, rules) => {
  const candidates = Array.isArray(features) ? features : [];
  const specificRules = Array.isArray(rules) ? rules : [];
  candidates.forEach((feature) => {
    for (const rule of specificRules) {
      if (!ruleMatchesFeature(feature, rule)) continue;
      const key = featureOverrideKey(feature);
      if (key) {
        overrides[key] = {
          color: rule.color,
          caption: rule.cap
        };
      }
      break;
    }
  });
};

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

const applyItemFeatureOverrides = ({
  svg,
  item,
  fillOverrides,
  strokeOverrides
}) => {
  let changed = false;
  const elementsById = featureElementsById(svg);
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
    const fillOverride = fillOverrides[stableKey];
    const fillValue = normalizePaint(
      fillOverride && typeof fillOverride === 'object' && hasOwn(fillOverride, 'color')
        ? fillOverride.color
        : fillOverride,
      'feature fill'
    );
    if (fillValue) {
      filterFeatureFillTargets(elements).forEach((element) => {
        if (element.getAttribute('fill') !== fillValue) {
          element.setAttribute('fill', fillValue);
          changed = true;
        }
      });
    }

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
      return legendChanged || strokeChanged;
    }
  })
});

export const prepareCandidateRenderCommit = ({
  results,
  catalog,
  mode = '',
  featureColorOverrides,
  featureStrokeOverrides,
  legendStrokeOverrides = {},
  manualSpecificRules = [],
  legacyFeatures = [],
  preparedFeatureState = null,
  suppressPairwiseIdentityLegend = false,
  transformSvg = null,
  sanitizer = globalThis.DOMPurify || globalThis.window?.DOMPurify
}) => {
  const featureState = preparedFeatureState || featureStateFromCatalog(catalog, { mode });
  const fillOverrides = cloneJsonValue(featureColorOverrides, {});
  const strokeOverrides = cloneJsonValue(featureStrokeOverrides, {});
  const migrationOptions = { legacyFeatures };

  migrateLegacyFeatureOverrides(
    fillOverrides,
    featureState.extractedFeatures,
    migrationOptions
  );
  migrateLegacyFeatureOverrides(
    strokeOverrides,
    featureState.extractedFeatures,
    migrationOptions
  );
  replaceRuleDerivedFillOverrides(
    fillOverrides,
    featureState.extractedFeatures,
    manualSpecificRules
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
      const overridesChanged = applyItemFeatureOverrides({
        svg,
        item,
        fillOverrides,
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
      return legendChanged || overridesChanged || legendStrokeChanged || callerChanged;
    }
  });

  return {
    results: processedResults,
    featureState,
    featureColorOverrides: fillOverrides,
    featureStrokeOverrides: strokeOverrides
  };
};
