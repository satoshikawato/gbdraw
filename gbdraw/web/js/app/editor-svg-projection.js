import { resolveColorToHex } from './color-utils.js';
import {
  filterFeatureFillTargets,
  getFeatureElementIndex,
  getFeatureElements
} from './feature-dom.js';
import { getFeatureCaption, ruleMatchesFeature } from './feature-utils.js';
import {
  featureVisibilityRuleMatchesFeature,
  resolveEffectiveFeatureVisibility
} from './feature-visibility.js';
import {
  getAllFeatureLegendGroups,
  PAIRWISE_LEGEND_SELECTOR,
  parseTransformXY
} from './legend/utils.js';
import { biologicalFeatureKey } from '../services/feature-catalog.js';
import {
  featureOverrideKey,
  getFeatureOverride
} from '../services/feature-override-identity.js';
import { cloneJsonValue } from '../services/json-clone.js';
import { recordStructuralMetric } from '../services/runtime-test-hooks.js';
import { createDomMutationJournal } from './dom-mutation-journal.js';

export const EDITABLE_LABEL_SELECTOR = 'text[data-label-editable="true"]';
export const LABEL_VISIBILITY_PREVIEW_ATTRIBUTE = 'data-gbdraw-label-visibility-preview';
export const LABEL_RENDERER_DISPLAY_ATTRIBUTE = 'data-gbdraw-label-renderer-display';

const text = (value) => String(value ?? '').trim();
const hasOwn = (object, key) => Object.prototype.hasOwnProperty.call(object || {}, key);
const setAttribute = (mutation, element, name, value) => (
  mutation
    ? mutation.setAttribute(element, name, value)
    : (element.setAttribute(name, value), true)
);
const removeAttribute = (mutation, element, name) => (
  mutation
    ? mutation.removeAttribute(element, name)
    : (element.removeAttribute(name), true)
);
const setTextContent = (mutation, element, value) => (
  mutation
    ? mutation.setTextContent(element, value)
    : (element.textContent = String(value ?? ''), true)
);

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
 * @property {Array<object>} [deletedLegendEntries]
 * @property {Array<string>} [originalLegendOrder]
 * @property {Array<string>|Set<string>} [addedLegendCaptions]
 * @property {Array<object>} [manualSpecificRules]
 * @property {Array<object>} [featureVisibilityManualRules]
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

export const setLabelText = (textElement, value, { mutation = null } = {}) => {
  if (!textElement) return false;
  const nextText = String(value ?? '');
  if (getLabelText(textElement) === nextText) return false;
  const textPath = textElement.querySelector?.('textPath');
  if (textPath) setTextContent(mutation, textPath, nextText);
  else setTextContent(mutation, textElement, nextText);
  return true;
};

export const normalizeLabelVisibilityMode = (value) => {
  const normalized = text(value).toLowerCase();
  return normalized === 'on' || normalized === 'off' ? normalized : 'default';
};

export const applyLabelVisibility = (textElement, modeRaw, {
  markPreview = false,
  mutation = null
} = {}) => {
  if (!textElement) return false;
  const mode = normalizeLabelVisibilityMode(modeRaw);
  const previewMode = textElement.getAttribute?.(LABEL_VISIBILITY_PREVIEW_ATTRIBUTE);
  let changed = false;
  if (mode === 'off' || mode === 'on') {
    if (markPreview && previewMode === null) {
      const rendererDisplay = textElement.getAttribute?.('display');
      setAttribute(mutation, textElement, LABEL_RENDERER_DISPLAY_ATTRIBUTE, rendererDisplay ?? '');
      changed = true;
    }
    if (markPreview && previewMode !== mode) {
      setAttribute(mutation, textElement, LABEL_VISIBILITY_PREVIEW_ATTRIBUTE, mode);
      changed = true;
    }
    if (mode === 'off') {
      if (textElement.getAttribute?.('display') !== 'none') {
        setAttribute(mutation, textElement, 'display', 'none');
        changed = true;
      }
    } else if (textElement.getAttribute?.('display') !== null) {
      removeAttribute(mutation, textElement, 'display');
      changed = true;
    }
    return changed;
  }
  if (previewMode === null) return false;
  const rendererDisplay = textElement.getAttribute?.(LABEL_RENDERER_DISPLAY_ATTRIBUTE);
  if (rendererDisplay) setAttribute(mutation, textElement, 'display', rendererDisplay);
  else removeAttribute(mutation, textElement, 'display');
  removeAttribute(mutation, textElement, LABEL_VISIBILITY_PREVIEW_ATTRIBUTE);
  removeAttribute(mutation, textElement, LABEL_RENDERER_DISPLAY_ATTRIBUTE);
  changed = true;
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
  legendColorOverrides = {},
  mutation = null
} = {}) => {
  if (!svg) return 0;
  let changedCount = 0;
  Object.entries(legendColorOverrides || {}).forEach(([caption, color]) => {
    const normalized = normalizeEditorPaint(color, 'legend color');
    if (!normalized) return;
    getLegendSwatches(svg, caption).forEach((swatch) => {
      if (swatch.getAttribute('fill') === normalized) return;
      setAttribute(mutation, swatch, 'fill', normalized);
      changedCount += 1;
    });
  });
  return changedCount;
};

const applyStrokeAttributes = (element, overrides, label, mutation = null) => {
  if (!overrides || typeof overrides !== 'object') return false;
  const strokeColor = hasOwn(overrides, 'strokeColor')
    ? normalizeEditorPaint(overrides.strokeColor, `${label} color`)
    : '';
  const strokeWidth = hasOwn(overrides, 'strokeWidth')
    ? normalizeEditorStrokeWidth(overrides.strokeWidth, `${label} width`)
    : null;
  let changed = false;
  if (strokeColor && element.getAttribute('stroke') !== strokeColor) {
    setAttribute(mutation, element, 'stroke', strokeColor);
    changed = true;
  }
  if (strokeWidth !== null && element.getAttribute('stroke-width') !== String(strokeWidth)) {
    setAttribute(mutation, element, 'stroke-width', strokeWidth);
    changed = true;
  }
  return changed;
};

export const applyStrokeOverridesToSvg = ({
  svg,
  features = [],
  legendStrokeOverrides = {},
  featureStrokeOverrides = {},
  mutation = null
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
      if (applyStrokeAttributes(element, overrides, 'feature stroke', mutation)) changedCount += 1;
    });
  };

  Object.entries(legendStrokeOverrides || {}).forEach(([caption, overrides]) => {
    if (!overrides || typeof overrides !== 'object') return;
    (Array.isArray(features) ? features : [])
      .filter((feature) => getFeatureCaption(feature) === caption)
      .forEach((feature) => applyToFeature(feature, overrides));
    getLegendSwatches(svg, caption).forEach((swatch) => {
      if (applyStrokeAttributes(swatch, overrides, 'legend stroke', mutation)) changedCount += 1;
    });
  });

  (Array.isArray(features) ? features : []).forEach((feature) => {
    applyToFeature(feature, getFeatureOverride(featureStrokeOverrides, feature));
  });
  return changedCount;
};

export const applyLegendEntryStateToSvg = ({ svg, legendEntries = [], mutation = null } = {}) => {
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
        setAttribute(mutation, entryGroup, 'data-legend-key', caption);
        changedCount += 1;
      }
      const label = entryGroup.querySelector('text');
      if (label && String(label.textContent || '') !== caption) {
        setTextContent(mutation, label, caption);
        changedCount += 1;
      }
      const color = normalizeEditorPaint(entry?.color, 'legend color');
      if (!color) return;
      const swatch = Array.from(entryGroup.querySelectorAll('path')).find((path) => {
        const fill = path.getAttribute('fill');
        return fill && fill !== 'none' && !fill.startsWith('url(');
      });
      if (swatch && swatch.getAttribute('fill') !== color) {
        setAttribute(mutation, swatch, 'fill', color);
        changedCount += 1;
      }
    });
  });
  return changedCount;
};

const directLegendEntryGroups = (targetGroup) => {
  const direct = Array.from(targetGroup?.children || []).filter(
    (child) => child.tagName?.toLowerCase() === 'g' && child.hasAttribute?.('data-legend-key')
  );
  return direct.length > 0
    ? direct
    : Array.from(targetGroup?.querySelectorAll?.('g[data-legend-key]') || []);
};

const legendEntryAnchor = (entryGroup) => {
  const groupOffset = parseTransformXY(entryGroup?.getAttribute?.('transform'));
  const anchor = entryGroup?.querySelector?.('text') || Array.from(
    entryGroup?.querySelectorAll?.('path') || []
  ).find((path) => {
    const fill = path.getAttribute('fill');
    return fill && fill !== 'none' && !fill.startsWith('url(');
  });
  const anchorOffset = parseTransformXY(anchor?.getAttribute?.('transform'));
  return { x: groupOffset.x + anchorOffset.x, y: groupOffset.y + anchorOffset.y };
};

const moveLegendEntryToAnchor = (entryGroup, anchor, mutation) => {
  const current = legendEntryAnchor(entryGroup);
  const deltaX = anchor.x - current.x;
  const deltaY = anchor.y - current.y;
  if (Math.abs(deltaX) < 1e-6 && Math.abs(deltaY) < 1e-6) return 0;
  let changed = 0;
  entryGroup.querySelectorAll?.('[transform]').forEach((node) => {
    const position = parseTransformXY(node.getAttribute('transform'));
    if (mutation.setAttribute(
      node,
      'transform',
      `translate(${position.x + deltaX}, ${position.y + deltaY})`
    )) changed += 1;
  });
  return changed;
};

const projectedLegendEntryColor = (entry) => normalizeEditorPaint(entry?.color, 'legend color');

const applyLegendStructureToSvg = ({
  svg,
  legendEntries = [],
  deletedLegendEntries = [],
  originalLegendOrder = [],
  addedLegendCaptions = [],
  mutation
} = {}) => {
  const desiredEntries = (Array.isArray(legendEntries) ? legendEntries : [])
    .filter((entry) => text(entry?.caption));
  const originalCaptions = new Set(
    (Array.isArray(originalLegendOrder) ? originalLegendOrder : []).map(text).filter(Boolean)
  );
  const addedCaptions = new Set(
    (Array.isArray(addedLegendCaptions) ? addedLegendCaptions : []).map(text).filter(Boolean)
  );
  const deletedCaptions = new Set();
  (Array.isArray(deletedLegendEntries) ? deletedLegendEntries : []).forEach((entry) => {
    const caption = text(entry?.caption);
    const originalCaption = text(entry?.originalCaption);
    if (caption) deletedCaptions.add(caption);
    if (originalCaption) deletedCaptions.add(originalCaption);
  });
  const structureInitialized = originalCaptions.size > 0
    || addedCaptions.size > 0
    || deletedCaptions.size > 0;
  if (!structureInitialized) return 0;

  let changedCount = 0;
  getAllFeatureLegendGroups(svg).forEach((targetGroup) => {
    const initialGroups = directLegendEntryGroups(targetGroup);
    if (initialGroups.length === 0 && desiredEntries.length === 0) return;
    const groupsByCaption = new Map(
      initialGroups.map((group) => [text(group.getAttribute('data-legend-key')), group])
    );
    const used = new Set();
    const desiredGroups = [];

    desiredEntries.forEach((entry) => {
      const caption = text(entry.caption);
      const sourceCaption = text(entry.originalCaption) || caption;
      let group = groupsByCaption.get(caption) || groupsByCaption.get(sourceCaption) || null;
      if (group && used.has(group)) group = null;
      if (!group) {
        const isAdded = addedCaptions.has(caption)
          || (!originalCaptions.has(caption) && !originalCaptions.has(sourceCaption));
        if (isAdded) {
          const template = initialGroups[0] || desiredGroups[0] || null;
          if (!template?.cloneNode) {
            throw new Error(`Cannot replay added Legend entry "${caption}" without a renderer template.`);
          }
          group = template.cloneNode(true);
          mutation.setAttribute(group, 'data-legend-key', caption);
          mutation.removeAttribute(group, 'data-legend-owner');
          const label = group.querySelector('text');
          if (label) mutation.setTextContent(label, caption);
          const color = projectedLegendEntryColor(entry);
          if (color) {
            const swatch = Array.from(group.querySelectorAll('path')).find((path) => {
              const fill = path.getAttribute('fill');
              return fill && fill !== 'none' && !fill.startsWith('url(');
            });
            if (swatch) mutation.setAttribute(swatch, 'fill', color);
          }
          mutation.appendChild(targetGroup, group);
          changedCount += 1;
        }
      }
      if (!group) return;
      if (text(group.getAttribute('data-legend-key')) !== caption) {
        mutation.setAttribute(group, 'data-legend-key', caption);
        const label = group.querySelector('text');
        if (label) mutation.setTextContent(label, caption);
        changedCount += 1;
      }
      used.add(group);
      desiredGroups.push(group);
    });

    initialGroups.forEach((group) => {
      if (used.has(group)) return;
      if (mutation.removeNode(group)) changedCount += 1;
    });

    const slots = initialGroups
      .filter((group) => group.parentNode === targetGroup || group.parentElement === targetGroup)
      .map(legendEntryAnchor)
      .sort((left, right) => {
        const yDelta = left.y - right.y;
        return Math.abs(yDelta) < 1 ? left.x - right.x : yDelta;
      });
    while (slots.length < desiredGroups.length) {
      const last = slots.at(-1) || { x: 0, y: 0 };
      const previous = slots.at(-2) || null;
      const stepY = previous && Math.abs(last.y - previous.y) >= 1
        ? last.y - previous.y
        : 20;
      slots.push({ x: last.x, y: last.y + stepY });
    }
    desiredGroups.forEach((group, index) => {
      changedCount += moveLegendEntryToAnchor(group, slots[index], mutation);
    });
    const remaining = directLegendEntryGroups(targetGroup).filter((group) => !desiredGroups.includes(group));
    const nextOrder = [...desiredGroups, ...remaining];
    const currentOrder = directLegendEntryGroups(targetGroup);
    if (nextOrder.some((group, index) => currentOrder[index] !== group)) {
      nextOrder.forEach((group) => mutation.appendChild(targetGroup, group));
      changedCount += 1;
    }
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

const labelProjection = (svg, state, { resetLabelState = false, mutation = null } = {}) => {
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
      if (source !== null && setLabelText(label, source, { mutation })) textCount += 1;
      if (applyLabelVisibility(label, 'default', { mutation })) visibilityCount += 1;
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
    if (setLabelText(candidates[0], nextText, { mutation })) textCount += 1;
  });

  Object.entries(state.labelTextBulkOverrides).forEach(([source, nextText]) => {
    labels.forEach((label) => {
      if (claimed.has(label)) return;
      const labelSource = String(label.getAttribute('data-label-source-text') ?? getLabelText(label));
      if (labelSource === source && setLabelText(label, nextText, { mutation })) textCount += 1;
    });
  });

  Object.entries(state.labelVisibilityOverrides).forEach(([featureId, mode]) => {
    const candidates = byFeatureId.get(featureId) || [];
    if (candidates.length !== 1) {
      unresolved.add(featureId);
      return;
    }
    if (applyLabelVisibility(candidates[0], mode, { markPreview: true, mutation })) visibilityCount += 1;
  });
  return {
    textCount,
    visibilityCount,
    unresolvedLabelFeatureIds: Array.from(unresolved).sort()
  };
};

const removePairwiseIdentityLegend = (svg, mutation = null) => {
  let changedCount = 0;
  svg.querySelectorAll(PAIRWISE_LEGEND_SELECTOR).forEach((legend) => {
    if (mutation) mutation.removeNode(legend);
    else legend.remove();
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
  let applicationCount = 0;
  const state = {
    // Validated Feature metadata is a large borrowed owner. The projection
    // never mutates, clones, or serializes this graph.
    features: Array.isArray(input.features) ? input.features : [],
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
    deletedLegendEntries: cloneJsonValue(input.deletedLegendEntries, []),
    originalLegendOrder: cloneJsonValue(input.originalLegendOrder, []),
    addedLegendCaptions: cloneJsonValue(
      input.addedLegendCaptions instanceof Set
        ? Array.from(input.addedLegendCaptions)
        : input.addedLegendCaptions,
      []
    ),
    featureVisibilityManualRules: cloneJsonValue(input.featureVisibilityManualRules, []),
    suppressPairwiseIdentityLegend: Boolean(input.suppressPairwiseIdentityLegend)
  };
  const structuralMetrics = Object.freeze({
    editorProjectionFullFeatureCloneCount: 0,
    editorProjectionFullFeatureSerializationCount: 0,
    get editorProjectionUnusedBuildCount() {
      return applicationCount === 0 ? 1 : 0;
    }
  });
  recordStructuralMetric('editorProjectionBuildCount');
  recordStructuralMetric('editorProjectionBorrowedFeatureOwnerCount', 1, {
    featureOwner: state.features
  });
  recordStructuralMetric('editorProjectionFullFeatureCloneCount', 0);
  recordStructuralMetric('editorProjectionFullFeatureSerializationCount', 0);
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
    applicationCount += 1;
    recordStructuralMetric('editorProjectionApplicationCount');
    recordStructuralMetric('editorProjectionUnusedBuildCount', 0);
    const elementsById = getFeatureElementIndex(svg);
    const bindings = featureBindings(state.features, item);
    bindings.forEach(({ feature, overrideKey, renderedId }) => {
      if (requireFeatureBindings && (elementsById.get(renderedId) || []).length === 0) {
        throw new Error('Sanitized SVG content is missing a rendered feature binding.');
      }
      const fillOverride = state.featureColorOverrides[overrideKey]
        ?? (feature ? getFeatureOverride(state.featureColorOverrides, feature) : undefined);
      normalizeEditorPaint(
        fillOverride && typeof fillOverride === 'object' && hasOwn(fillOverride, 'color')
          ? fillOverride.color
          : fillOverride,
        'feature fill'
      );
      const strokeOverride = state.featureStrokeOverrides[overrideKey]
        ?? (feature ? getFeatureOverride(state.featureStrokeOverrides, feature) : undefined);
      if (strokeOverride && typeof strokeOverride === 'object') {
        if (hasOwn(strokeOverride, 'strokeColor')) {
          normalizeEditorPaint(strokeOverride.strokeColor, 'feature stroke color');
        }
        if (hasOwn(strokeOverride, 'strokeWidth')) {
          normalizeEditorStrokeWidth(strokeOverride.strokeWidth, 'feature stroke width');
        }
      }
    });
    Object.values(state.legendColorOverrides).forEach((color) => (
      normalizeEditorPaint(color, 'legend color')
    ));
    state.legendEntries.forEach((entry) => normalizeEditorPaint(entry?.color, 'legend color'));
    Object.values(state.legendStrokeOverrides).forEach((strokeOverride) => {
      if (!strokeOverride || typeof strokeOverride !== 'object') return;
      if (hasOwn(strokeOverride, 'strokeColor')) {
        normalizeEditorPaint(strokeOverride.strokeColor, 'legend stroke color');
      }
      if (hasOwn(strokeOverride, 'strokeWidth')) {
        normalizeEditorStrokeWidth(strokeOverride.strokeWidth, 'legend stroke width');
      }
    });

    const mutation = createDomMutationJournal();
    let featureFillCount = 0;
    let featureStrokeCount = 0;
    let featureVisibilityCount = 0;
    try {
      bindings.forEach(({ feature, overrideKey, renderedId }) => {
        const elements = elementsById.get(renderedId) || [];
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
          mutation.setAttribute(element, 'fill', fill);
          featureFillCount += 1;
        });
      }
      const strokeOverride = state.featureStrokeOverrides[overrideKey]
        ?? (feature ? getFeatureOverride(state.featureStrokeOverrides, feature) : undefined);
      elements.forEach((element) => {
        if (applyStrokeAttributes(element, strokeOverride, 'feature stroke', mutation)) featureStrokeCount += 1;
      });
      const hasRenderedVisibility = hasOwn(state.featureVisibilityOverrides, renderedId);
      const hasStableVisibility = hasOwn(state.featureVisibilityOverrides, overrideKey);
      const hasVisibilityOverride = hasRenderedVisibility || hasStableVisibility;
      const hasManualVisibility = Boolean(feature) && state.featureVisibilityManualRules.some(
        (rule) => featureVisibilityRuleMatchesFeature(feature, rule)
      );
      const visibility = resolveEffectiveFeatureVisibility(
        feature || renderedId,
        state.featureVisibilityOverrides,
        null,
        state.featureVisibilityManualRules,
        { overrideKeys: [renderedId, overrideKey] }
      );
      elements.forEach((element) => {
        if (visibility === 'off') {
          if (element.getAttribute('display') !== 'none') {
            mutation.setAttribute(element, 'display', 'none');
            featureVisibilityCount += 1;
          }
        } else if (
          (hasVisibilityOverride || hasManualVisibility || resetFeatureVisibility)
          && element.getAttribute('display') === 'none'
        ) {
          mutation.removeAttribute(element, 'display');
          featureVisibilityCount += 1;
        }
      });
      });

      const labelResult = labelProjection(svg, state, { resetLabelState, mutation });
      const projectedLegendEntries = state.legendEntries.map((entry) => ({
        ...entry,
        color: hasOwn(state.legendColorOverrides, entry?.caption)
          ? state.legendColorOverrides[entry.caption]
          : null
      }));
      const legendStructureCount = applyLegendStructureToSvg({
        svg,
        legendEntries: state.legendEntries,
        deletedLegendEntries: state.deletedLegendEntries,
        originalLegendOrder: state.originalLegendOrder,
        addedLegendCaptions: state.addedLegendCaptions,
        mutation
      });
      const legendCount = legendStructureCount + applyLegendEntryStateToSvg({
        svg,
        legendEntries: projectedLegendEntries,
        mutation
      }) + applyLegendColorOverridesToSvg({
        svg,
        legendColorOverrides: state.legendColorOverrides,
        mutation
      }) + applyStrokeOverridesToSvg({
        svg,
        features: state.features,
        legendStrokeOverrides: state.legendStrokeOverrides,
        featureStrokeOverrides: {},
        mutation
      });
      const suppressionCount = state.suppressPairwiseIdentityLegend
        ? removePairwiseIdentityLegend(svg, mutation)
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
      mutation.commit();
      return result;
    } catch (error) {
      mutation.rollback();
      throw error;
    }
  };

  return {
    featureColorOverrides: state.featureColorOverrides,
    featureStrokeOverrides: state.featureStrokeOverrides,
    structuralMetrics,
    project
  };
};
