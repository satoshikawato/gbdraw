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
import { serializeCleanSvg } from '../services/svg-serialization.js';
import { sanitizeSvgContent } from '../services/svg-sanitization.js';

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

const parseSanitizedSvg = (content, sanitizer) => {
  const sanitized = sanitizeSvgContent(content, sanitizer);
  if (typeof DOMParser !== 'function') {
    throw new Error('SVG parsing is unavailable.');
  }
  const document = new DOMParser().parseFromString(sanitized, 'image/svg+xml');
  if (
    document.querySelector('parsererror')
    || String(document.documentElement?.localName || '').toLowerCase() !== 'svg'
  ) {
    throw new Error('The diagram engine returned malformed SVG content.');
  }
  return document.documentElement;
};

const applyItemFeatureOverrides = ({
  svg,
  item,
  fillOverrides,
  strokeOverrides
}) => {
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
        element.setAttribute('fill', fillValue);
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
      if (strokeValue) element.setAttribute('stroke', strokeValue);
      if (strokeWidth !== null) element.setAttribute('stroke-width', strokeWidth);
    });
  });
};

export const prepareCandidateRenderCommit = ({
  results,
  catalog,
  mode = '',
  featureColorOverrides,
  featureStrokeOverrides,
  manualSpecificRules = [],
  legacyFeatures = [],
  sanitizer = globalThis.DOMPurify
}) => {
  const featureState = featureStateFromCatalog(catalog, { mode });
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
  const processedResults = results.map((result, resultIndex) => {
    const item = itemsByIndex.get(resultIndex);
    if (!item) {
      throw new Error('The diagram engine returned incomplete feature metadata.');
    }
    const svg = parseSanitizedSvg(result?.content, sanitizer);
    applyItemFeatureOverrides({
      svg,
      item,
      fillOverrides,
      strokeOverrides
    });
    return {
      ...result,
      content: serializeCleanSvg(svg)
    };
  });

  return {
    results: processedResults,
    featureState,
    featureColorOverrides: fillOverrides,
    featureStrokeOverrides: strokeOverrides
  };
};
