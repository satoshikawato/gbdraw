import { resolveColorToHex } from './color-utils.js';
import {
  filterFeatureFillTargets,
  getFeatureElementIndex,
  getFeatureElements
} from './feature-dom.js';
import { getFeatureCaption, ruleMatchesFeature } from './feature-utils.js';
import { getAllFeatureLegendGroups, PAIRWISE_LEGEND_SELECTOR } from './legend/utils.js';
import { biologicalFeatureKey } from '../services/feature-catalog.js';
import {
  featureOverrideKey,
  getFeatureOverride
} from '../services/feature-override-identity.js';
import { cloneJsonValue } from '../services/json-clone.js';

export const EDITABLE_LABEL_SELECTOR = 'text[data-label-editable="true"]';
export const LABEL_VISIBILITY_PREVIEW_ATTRIBUTE = 'data-gbdraw-label-visibility-preview';

const text = (value) => String(value ?? '').trim();
const hasOwn = (object, key) => Object.prototype.hasOwnProperty.call(object || {}, key);

/**
 * @typedef {Object} EditorSvgProjectionResult
 * @property {boolean} changed Whether any supported SVG presentation changed.
 * @property {number} featureFillCount Number of feature fill elements changed.
 * @property {number} featureStrokeCount Number of feature stroke elements changed.
 * @property {number} featureVisibilityCount Number of feature visibility elements changed.
 * @property {number} labelTextCount Number of label text elements changed.
 * @property {number} labelVisibilityCount Number of label visibility elements changed.
 * @property {number} legendCount Number of legend elements changed.
 * @property {number} suppressionCount Number of suppressed legend groups removed.
 * @property {string[]} unresolvedLabelFeatureIds Ambiguous or unavailable label bindings skipped safely.
 */

/**
 * @typedef {Object} EditorSvgProjectionInput
 * @property {Array<object>} [features]
 * @property {Object<string, *>} [featureColorOverrides]
 * @property {Object<string, *>} [featureStrokeOverrides]
 * @property {Object<string, string>} [featureVisibilityOverrides]
 * @property {Object<string, string>} [labelTextFeatureOverrides]
 * @property {Object<string, string>} [labelTextBulkOverrides]
 * @property {Object<string, string>} [labelTextFeatureOverrideSources]
 * @property {Object<string, string>} [labelVisibilityOverrides]
 * @property {Object<string, string>} [legendColorOverrides]
 * @property {Object<string, object>} [legendStrokeOverrides]
 * @property {Array<object>} [legendEntries]
 * @property {Array<object>} [manualSpecificRules]
 * @property {boolean} [suppressPairwiseIdentityLegend]
 */

export const normalizeEditorPaint = (value, label) => {
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

export const normalizeEditorStrokeWidth = (value, label = 'feature stroke width') => {
  if (value === null || value === undefined || value === '') return null;
  const width = Number(value);
  if (!Number.isFinite(width) || width < 0) {
    throw new Error(`Invalid ${label} override in the committed editor state.`);
  }
  return width;
};

export const getLabelText = (textElement) => {
  const textPath = textElement?.querySelector?.('textPath');
  return textPath ? (textPath.textContent || '') : (textElement?.textContent || '');
};

export const setLabelText = (textElement, value) => {
  if (!textElement) return false;
  const nextText = String(value ?? '');
  if (getLabelText(textElement) === nextText) return false;
  const textPath = textElement.querySelector?.('textPath');
  if (textPath) textPath.textContent = nextText;
  else textElement.textContent = nextText;
  return true;
};

export const normalizeLabelVisibilityMode = (value) => {
  const normalized = text(value).toLowerCase();
  return normalized === 'on' || normalized === 'off' ? normalized : 'default';
};

export const applyLabelVisibility = (textElement, modeRaw, {
  markPreview = false
} = {}) => {
  if (!textElement) return false;
  const mode = normalizeLabelVisibilityMode(modeRaw);
  const previewHidden = textElement.hasAttribute?.(LABEL_VISIBILITY_PREVIEW_ATTRIBUTE);
  let changed = false;
  if (mode === 'off') {
    if (markPreview && !previewHidden) {
      textElement.setAttribute(LABEL_VISIBILITY_PREVIEW_ATTRIBUTE, 'off');
      changed = true;
    }
    if (textElement.getAttribute?.('display') !== 'none') {
      textElement.setAttribute('display', 'none');
      changed = true;
    }
    return changed;
  }
  if (previewHidden) {
    textElement.removeAttribute(LABEL_VISIBILITY_PREVIEW_ATTRIBUTE);
    changed = true;
  }
  if (textElement.getAttribute?.('display') === 'none') {
    textElement.removeAttribute('display');
    changed = true;
  }
  return changed;
};

export const resetLabelsToSourceText = (svg) => {
  let changedCount = 0;
  svg?.querySelectorAll?.(EDITABLE_LABEL_SELECTOR)?.forEach((textElement) => {
    const sourceText = textElement.getAttribute('data-label-source-text');
    if (sourceText !== null && setLabelText(textElement, sourceText)) changedCount += 1;
  });
  return changedCount;
};

export const getLegendSwatches = (svg, caption) => {
  const swatches = [];
  getAllFeatureLegendGroups(svg).forEach((targetGroup) => {
    const entryGroup = Array.from(targetGroup.querySelectorAll('g[data-legend-key]'))
      .find((entry) => entry.getAttribute('data-legend-key') === caption);
    const swatch = Array.from(entryGroup?.querySelectorAll?.('path') || []).find((path) => {
      const fill = path.getAttribute('fill');
      return fill && fill !== 'none' && !fill.startsWith('url(');
    });
    if (swatch) swatches.push(swatch);
  });
  return swatches;
};

export const applyLegendColorOverridesToSvg = ({
  svg,
  legendColorOverrides = {}
} = {}) => {
  if (!svg) return 0;
  let changedCount = 0;
  Object.entries(legendColorOverrides || {}).forEach(([caption, color]) => {
    const normalized = normalizeEditorPaint(color, 'legend color');
    if (!normalized) return;
    getLegendSwatches(svg, caption).forEach((swatch) => {
      if (swatch.getAttribute('fill') === normalized) return;
      swatch.setAttribute('fill', normalized);
      changedCount += 1;
    });
  });
  return changedCount;
};

const applyStrokeAttributes = (element, overrides, label) => {
  if (!overrides || typeof overrides !== 'object') return false;
  const strokeColor = hasOwn(overrides, 'strokeColor')
    ? normalizeEditorPaint(overrides.strokeColor, `${label} color`)
    : '';
  const strokeWidth = hasOwn(overrides, 'strokeWidth')
    ? normalizeEditorStrokeWidth(overrides.strokeWidth, `${label} width`)
    : null;
  let changed = false;
  if (strokeColor && element.getAttribute('stroke') !== strokeColor) {
    element.setAttribute('stroke', strokeColor);
    changed = true;
  }
  if (strokeWidth !== null && element.getAttribute('stroke-width') !== String(strokeWidth)) {
    element.setAttribute('stroke-width', strokeWidth);
    changed = true;
  }
  return changed;
};

export const applyStrokeOverridesToSvg = ({
  svg,
  features = [],
  legendStrokeOverrides = {},
  featureStrokeOverrides = {}
} = {}) => {
  if (!svg) return 0;
  const featureIndex = getFeatureElementIndex(svg);
  let changedCount = 0;
  const applyToFeature = (feature, overrides) => {
    if (!overrides || typeof overrides !== 'object') return;
    const featureId = text(
      feature?.rendered_svg_id
      || feature?.renderedSvgId
      || feature?.rendered_feature_svg_id
      || feature?.renderedFeatureSvgId
      || feature?.svg_id
      || feature?.svgId
    );
    getFeatureElements(svg, featureId, featureIndex).forEach((element) => {
      if (applyStrokeAttributes(element, overrides, 'feature stroke')) changedCount += 1;
    });
  };

  Object.entries(legendStrokeOverrides || {}).forEach(([caption, overrides]) => {
    if (!overrides || typeof overrides !== 'object') return;
    (Array.isArray(features) ? features : [])
      .filter((feature) => getFeatureCaption(feature) === caption)
      .forEach((feature) => applyToFeature(feature, overrides));
    getLegendSwatches(svg, caption).forEach((swatch) => {
      if (applyStrokeAttributes(swatch, overrides, 'legend stroke')) changedCount += 1;
    });
  });

  (Array.isArray(features) ? features : []).forEach((feature) => {
    applyToFeature(feature, getFeatureOverride(featureStrokeOverrides, feature));
  });
  return changedCount;
};

export const applyLegendEntryStateToSvg = ({ svg, legendEntries = [] } = {}) => {
  if (!svg) return 0;
  let changedCount = 0;
  (Array.isArray(legendEntries) ? legendEntries : []).forEach((entry) => {
    const caption = text(entry?.caption);
    const sourceCaption = text(entry?.originalCaption) || caption;
    if (!caption) return;
    getAllFeatureLegendGroups(svg).forEach((targetGroup) => {
      const entryGroup = Array.from(targetGroup.querySelectorAll('g[data-legend-key]'))
        .find((candidate) => {
          const key = text(candidate.getAttribute('data-legend-key'));
          return key === caption || key === sourceCaption;
        });
      if (!entryGroup) return;
      if (text(entryGroup.getAttribute('data-legend-key')) !== caption) {
        entryGroup.setAttribute('data-legend-key', caption);
        changedCount += 1;
      }
      const label = entryGroup.querySelector('text');
      if (label && String(label.textContent || '') !== caption) {
        label.textContent = caption;
        changedCount += 1;
      }
      const color = normalizeEditorPaint(entry?.color, 'legend color');
      if (!color) return;
      const swatch = Array.from(entryGroup.querySelectorAll('path')).find((path) => {
        const fill = path.getAttribute('fill');
        return fill && fill !== 'none' && !fill.startsWith('url(');
      });
      if (swatch && swatch.getAttribute('fill') !== color) {
        swatch.setAttribute('fill', color);
        changedCount += 1;
      }
    });
  });
  return changedCount;
};

const replaceRuleDerivedFillOverrides = (overrides, features, rules) => {
  const hashRules = [];
  const qualifierRules = [];
  (Array.isArray(rules) ? rules : []).forEach((rule) => {
    (text(rule?.qual).toLowerCase() === 'hash' ? hashRules : qualifierRules).push(rule);
  });
  (Array.isArray(features) ? features : []).forEach((feature) => {
    const rule = hashRules.find((candidate) => ruleMatchesFeature(feature, candidate))
      || qualifierRules.find((candidate) => ruleMatchesFeature(feature, candidate));
    const key = rule ? featureOverrideKey(feature) : '';
    if (key) overrides[key] = { color: rule.color, caption: rule.cap };
  });
};

const featureBindings = (features, item) => {
  if (Array.isArray(item?.features)) {
    const byStableKey = new Map((Array.isArray(features) ? features : []).map((feature) => [
      featureOverrideKey(feature),
      feature
    ]));
    return item.features.map((rendered) => {
      const stableKey = biologicalFeatureKey(rendered.recordKey, rendered.biologicalFeatureId);
      return {
        feature: byStableKey.get(stableKey) || null,
        overrideKey: stableKey,
        renderedId: text(rendered.svgId)
      };
    });
  }
  return (Array.isArray(features) ? features : []).map((feature) => ({
    feature,
    overrideKey: featureOverrideKey(feature),
    renderedId: text(
      feature?.rendered_svg_id
      || feature?.renderedSvgId
      || feature?.rendered_feature_svg_id
      || feature?.renderedFeatureSvgId
      || feature?.svg_id
      || feature?.svgId
    )
  })).filter((binding) => binding.renderedId);
};

const labelProjection = (svg, state, { resetLabelState = false } = {}) => {
  const labels = Array.from(svg.querySelectorAll?.(
    `${EDITABLE_LABEL_SELECTOR}, text[data-label-feature-id], text[data-label-source-text]`
  ) || []);
  const byFeatureId = new Map();
  labels.forEach((label) => {
    const featureId = text(label.getAttribute('data-label-feature-id'));
    if (!featureId) return;
    if (!byFeatureId.has(featureId)) byFeatureId.set(featureId, []);
    byFeatureId.get(featureId).push(label);
  });
  const claimed = new Set();
  const unresolved = new Set();
  let textCount = 0;
  let visibilityCount = 0;

  if (resetLabelState) {
    labels.forEach((label) => {
      const source = label.getAttribute('data-label-source-text');
      if (source !== null && setLabelText(label, source)) textCount += 1;
      if (applyLabelVisibility(label, 'default')) visibilityCount += 1;
    });
  }

  Object.entries(state.labelTextFeatureOverrides).forEach(([featureId, nextText]) => {
    let candidates = byFeatureId.get(featureId) || [];
    if (candidates.length === 0) {
      const source = String(state.labelTextFeatureOverrideSources[featureId] ?? '');
      if (source) {
        candidates = labels.filter((label) => (
          String(label.getAttribute('data-label-source-text') ?? getLabelText(label)) === source
        ));
      }
    }
    if (candidates.length !== 1) {
      unresolved.add(featureId);
      return;
    }
    claimed.add(candidates[0]);
    if (setLabelText(candidates[0], nextText)) textCount += 1;
  });

  Object.entries(state.labelTextBulkOverrides).forEach(([source, nextText]) => {
    labels.forEach((label) => {
      if (claimed.has(label)) return;
      const labelSource = String(label.getAttribute('data-label-source-text') ?? getLabelText(label));
      if (labelSource === source && setLabelText(label, nextText)) textCount += 1;
    });
  });

  Object.entries(state.labelVisibilityOverrides).forEach(([featureId, mode]) => {
    const candidates = byFeatureId.get(featureId) || [];
    if (candidates.length !== 1) {
      unresolved.add(featureId);
      return;
    }
    if (applyLabelVisibility(candidates[0], mode)) visibilityCount += 1;
  });
  return {
    textCount,
    visibilityCount,
    unresolvedLabelFeatureIds: Array.from(unresolved).sort()
  };
};

const removePairwiseIdentityLegend = (svg) => {
  let changedCount = 0;
  svg.querySelectorAll(PAIRWISE_LEGEND_SELECTOR).forEach((legend) => {
    legend.remove();
    changedCount += 1;
  });
  return changedCount;
};

/**
 * Build a deterministic projector from explicit editor state. It never reads UI
 * globals, dispatches user actions, invokes Python, or schedules diagram work.
 * Projection is synchronous: a thrown validation or binding error leaves result
 * admission to the caller, and a no-op returns changed=false.
 *
 * @param {EditorSvgProjectionInput} input
 */
export const createEditorSvgProjection = (input = {}) => {
  const state = {
    features: cloneJsonValue(input.features, []),
    featureColorOverrides: cloneJsonValue(input.featureColorOverrides, {}),
    featureStrokeOverrides: cloneJsonValue(input.featureStrokeOverrides, {}),
    featureVisibilityOverrides: cloneJsonValue(input.featureVisibilityOverrides, {}),
    labelTextFeatureOverrides: cloneJsonValue(input.labelTextFeatureOverrides, {}),
    labelTextBulkOverrides: cloneJsonValue(input.labelTextBulkOverrides, {}),
    labelTextFeatureOverrideSources: cloneJsonValue(input.labelTextFeatureOverrideSources, {}),
    labelVisibilityOverrides: cloneJsonValue(input.labelVisibilityOverrides, {}),
    legendColorOverrides: cloneJsonValue(input.legendColorOverrides, {}),
    legendStrokeOverrides: cloneJsonValue(input.legendStrokeOverrides, {}),
    legendEntries: cloneJsonValue(input.legendEntries, []),
    suppressPairwiseIdentityLegend: Boolean(input.suppressPairwiseIdentityLegend)
  };
  replaceRuleDerivedFillOverrides(
    state.featureColorOverrides,
    state.features,
    input.manualSpecificRules
  );

  const project = (svg, {
    item = null,
    requireFeatureBindings = false,
    resetFeatureVisibility = false,
    resetLabelState = false
  } = {}) => {
    if (!svg) throw new Error('Editor SVG projection requires an SVG root.');
    const elementsById = getFeatureElementIndex(svg);
    let featureFillCount = 0;
    let featureStrokeCount = 0;
    let featureVisibilityCount = 0;
    const bindings = featureBindings(state.features, item);
    bindings.forEach(({ feature, overrideKey, renderedId }) => {
      const elements = elementsById.get(renderedId) || [];
      if (requireFeatureBindings && elements.length === 0) {
        throw new Error('Sanitized SVG content is missing a rendered feature binding.');
      }
      const fillOverride = state.featureColorOverrides[overrideKey]
        ?? (feature ? getFeatureOverride(state.featureColorOverrides, feature) : undefined);
      const fill = normalizeEditorPaint(
        fillOverride && typeof fillOverride === 'object' && hasOwn(fillOverride, 'color')
          ? fillOverride.color
          : fillOverride,
        'feature fill'
      );
      if (fill) {
        filterFeatureFillTargets(elements).forEach((element) => {
          if (element.getAttribute('fill') === fill) return;
          element.setAttribute('fill', fill);
          featureFillCount += 1;
        });
      }
      const strokeOverride = state.featureStrokeOverrides[overrideKey]
        ?? (feature ? getFeatureOverride(state.featureStrokeOverrides, feature) : undefined);
      elements.forEach((element) => {
        if (applyStrokeAttributes(element, strokeOverride, 'feature stroke')) featureStrokeCount += 1;
      });
      const hasRenderedVisibility = hasOwn(state.featureVisibilityOverrides, renderedId);
      const hasStableVisibility = hasOwn(state.featureVisibilityOverrides, overrideKey);
      const hasVisibilityOverride = hasRenderedVisibility || hasStableVisibility;
      const visibility = String(hasRenderedVisibility
        ? state.featureVisibilityOverrides[renderedId]
        : state.featureVisibilityOverrides[overrideKey] ?? '').trim().toLowerCase();
      elements.forEach((element) => {
        if (visibility === 'off') {
          if (element.getAttribute('display') !== 'none') {
            element.setAttribute('display', 'none');
            featureVisibilityCount += 1;
          }
        } else if (
          (hasVisibilityOverride || resetFeatureVisibility)
          && element.getAttribute('display') === 'none'
        ) {
          element.removeAttribute('display');
          featureVisibilityCount += 1;
        }
      });
    });

    const labelResult = labelProjection(svg, state, { resetLabelState });
    const projectedLegendEntries = state.legendEntries.map((entry) => ({
      ...entry,
      color: hasOwn(state.legendColorOverrides, entry?.caption)
        ? state.legendColorOverrides[entry.caption]
        : null
    }));
    const legendCount = applyLegendEntryStateToSvg({ svg, legendEntries: projectedLegendEntries })
      + applyLegendColorOverridesToSvg({ svg, legendColorOverrides: state.legendColorOverrides })
      + applyStrokeOverridesToSvg({
        svg,
        features: state.features,
        legendStrokeOverrides: state.legendStrokeOverrides,
        featureStrokeOverrides: {}
      });
    const suppressionCount = state.suppressPairwiseIdentityLegend
      ? removePairwiseIdentityLegend(svg)
      : 0;
    const changed = featureFillCount + featureStrokeCount + featureVisibilityCount
      + labelResult.textCount + labelResult.visibilityCount + legendCount + suppressionCount > 0;
    /** @type {EditorSvgProjectionResult} */
    const result = {
      changed,
      featureFillCount,
      featureStrokeCount,
      featureVisibilityCount,
      labelTextCount: labelResult.textCount,
      labelVisibilityCount: labelResult.visibilityCount,
      legendCount,
      suppressionCount,
      unresolvedLabelFeatureIds: labelResult.unresolvedLabelFeatureIds
    };
    return result;
  };

  return {
    featureColorOverrides: state.featureColorOverrides,
    featureStrokeOverrides: state.featureStrokeOverrides,
    project
  };
};
