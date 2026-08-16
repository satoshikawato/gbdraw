import { resolveColorToHex } from './color-utils.js';
import {
  filterFeatureFillTargets,
  getFeatureElementIndex,
  getFeatureElements
} from './feature-dom.js';
import { getFeatureCaption, ruleMatchesFeature } from './feature-utils.js';
import { resolveEffectiveFeatureVisibility } from './feature-visibility.js';
import {
  COMPARISON_LEGEND_SELECTOR,
  getAllFeatureLegendGroups,
  getComparisonLegendGroup,
  PAIRWISE_LEGEND_SELECTOR,
  parseTransformXY
} from './legend/utils.js';
import {
  applyCompositionEdit,
  parseCompositionMetadata
} from './legend-layout/composition-actions.js';
import { biologicalFeatureKey } from '../services/feature-catalog.js';
import {
  featureOverrideKey,
  getFeatureOverride
} from '../services/feature-override-identity.js';
import { cloneJsonValue } from '../services/json-clone.js';
import { recordStructuralMetric } from '../services/runtime-test-hooks.js';
import { createDomMutationJournal } from './dom-mutation-journal.js';
import { applyFeatureVisibility } from './feature-visibility-dom.js';

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

export const resetLabelsToSourceText = (svg, { mutation = null } = {}) => {
  let changedCount = 0;
  svg?.querySelectorAll?.(EDITABLE_LABEL_SELECTOR)?.forEach((textElement) => {
    const sourceText = textElement.getAttribute('data-label-source-text');
    if (sourceText !== null && setLabelText(textElement, sourceText, { mutation })) changedCount += 1;
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

const LEGEND_ROW_TOLERANCE = 1;

const legendOrientation = (targetGroup, anchors, compositionOrientation = null) => {
  let owner = targetGroup;
  while (owner) {
    const explicit = text(owner.getAttribute?.('data-gbdraw-orientation')).toLowerCase();
    if (explicit === 'horizontal' || explicit === 'h') return 'horizontal';
    if (explicit === 'vertical' || explicit === 'v') return 'vertical';
    const ownerId = text(owner.id).toLowerCase();
    if (/(?:^|_)horizontal(?:_|$)|_h$/.test(ownerId)) return 'horizontal';
    if (/(?:^|_)vertical(?:_|$)|_v$/.test(ownerId)) return 'vertical';
    owner = owner.parentElement || owner.parentNode || null;
  }

  if (compositionOrientation === 'horizontal' || compositionOrientation === 'vertical') {
    return compositionOrientation;
  }

  if (anchors.length > 1) {
    const xs = anchors.map(({ x }) => x);
    const ys = anchors.map(({ y }) => y);
    const xSpan = Math.max(...xs) - Math.min(...xs);
    const ySpan = Math.max(...ys) - Math.min(...ys);
    if (xSpan >= LEGEND_ROW_TOLERANCE) return 'horizontal';
    if (ySpan >= LEGEND_ROW_TOLERANCE) return 'vertical';
  }
  return 'vertical';
};

const sortedLegendRows = (anchors) => {
  const sorted = [...anchors].sort((left, right) => {
    const yDelta = left.y - right.y;
    return Math.abs(yDelta) < LEGEND_ROW_TOLERANCE ? left.x - right.x : yDelta;
  });
  const rows = [];
  sorted.forEach((anchor) => {
    const row = rows.at(-1);
    if (!row || Math.abs(row.y - anchor.y) >= LEGEND_ROW_TOLERANCE) {
      rows.push({ y: anchor.y, anchors: [anchor] });
    } else {
      row.anchors.push(anchor);
    }
  });
  rows.forEach((row) => row.anchors.sort((left, right) => left.x - right.x));
  return rows;
};

const lastPositiveStep = (values, fallback) => {
  for (let index = values.length - 1; index > 0; index -= 1) {
    const step = values[index] - values[index - 1];
    if (Math.abs(step) >= LEGEND_ROW_TOLERANCE) return step;
  }
  return fallback;
};

const resolveLegendSlots = ({
  targetGroup,
  initialGroups,
  desiredCount,
  compositionOrientation = null
}) => {
  const anchors = initialGroups.map(legendEntryAnchor);
  const orientation = legendOrientation(targetGroup, anchors, compositionOrientation);
  const rows = sortedLegendRows(anchors);
  let topology = orientation;
  const slots = rows.flatMap((row) => row.anchors);

  if (orientation === 'vertical') {
    const stepY = lastPositiveStep(slots.map(({ y }) => y), 20);
    const x = slots.at(-1)?.x ?? 0;
    while (slots.length < desiredCount) {
      const y = (slots.at(-1)?.y ?? -stepY) + stepY;
      slots.push({ x, y });
    }
  } else if (rows.length <= 1) {
    const stepX = lastPositiveStep(slots.map(({ x }) => x), 20);
    const y = slots.at(-1)?.y ?? 0;
    while (slots.length < desiredCount) {
      const x = (slots.at(-1)?.x ?? -stepX) + stepX;
      slots.push({ x, y });
    }
  } else {
    topology = 'horizontal-wrapped';
    const prototypeRow = rows.reduce((longest, row) => (
      row.anchors.length > longest.anchors.length ? row : longest
    ), rows[0]);
    const columns = prototypeRow.anchors.map(({ x }) => x);
    const columnCount = columns.length;
    const stepY = lastPositiveStep(rows.map(({ y }) => y), 20);
    let row = rows.at(-1);
    let columnIndex = row.anchors.length;
    while (slots.length < desiredCount) {
      if (columnIndex >= columnCount) {
        row = { y: row.y + stepY, anchors: [] };
        columnIndex = 0;
      }
      const anchor = { x: columns[columnIndex] ?? 0, y: row.y };
      row.anchors.push(anchor);
      slots.push(anchor);
      columnIndex += 1;
    }
  }

  recordStructuralMetric('legendTopologyResolutionCount', 1, {
    topology,
    orientation,
    initialCount: initialGroups.length,
    desiredCount,
    targetGroupOwner: targetGroup
  });
  return slots;
};

const projectedLegendEntryColor = (entry) => normalizeEditorPaint(entry?.color, 'legend color');

const rendererLegendLayout = (svg) => {
  const hasCompositionMetadata = svg?.getAttribute?.('data-gbdraw-composition') !== null
    || svg?.getAttribute?.('data-gbdraw-composition-schema') !== null;
  if (!hasCompositionMetadata) return null;
  const metadata = parseCompositionMetadata(svg);
  if (!metadata.legend || metadata.legendSide === 'none') return null;
  const orientation = metadata.legendSide === 'top' || metadata.legendSide === 'bottom'
    ? 'horizontal'
    : 'vertical';
  return {
    metadata,
    orientation
  };
};

const measuredLegendTextWidth = (label) => {
  const fontSize = Number.parseFloat(label?.getAttribute?.('font-size')) || 16;
  const content = String(label?.textContent || '');
  try {
    const canvas = globalThis.document?.createElement?.('canvas');
    const context = canvas?.getContext?.('2d');
    if (context) {
      const fontStyle = label.getAttribute?.('font-style') || 'normal';
      const fontWeight = label.getAttribute?.('font-weight') || 'normal';
      const fontFamily = label.getAttribute?.('font-family') || 'sans-serif';
      context.font = `${fontStyle} ${fontWeight} ${fontSize}px ${fontFamily}`;
      const width = Number(context.measureText(content).width);
      if (Number.isFinite(width) && width > 0) return width;
    }
  } catch (_error) {
    // Detached SVG admission can still use the deterministic fallback below.
  }
  const measured = Number(label?.getComputedTextLength?.());
  if (Number.isFinite(measured) && measured > 0) return measured;
  const bboxWidth = Number(label?.getBBox?.().width);
  if (Number.isFinite(bboxWidth) && bboxWidth > 0) return bboxWidth;
  return Math.max(fontSize * 0.6, Array.from(content).length * fontSize * 0.6);
};

const comparisonLegendWidth = (targetGroup) => {
  let owner = targetGroup;
  while (owner && text(owner.id) !== 'legend') owner = owner.parentElement || owner.parentNode || null;
  const comparison = owner?.querySelector?.(COMPARISON_LEGEND_SELECTOR) || null;
  if (!comparison) return 0;
  const bboxWidth = Number(comparison.getBBox?.().width);
  if (Number.isFinite(bboxWidth) && bboxWidth > 0) return bboxWidth;
  let width = 0;
  Array.from(comparison.querySelectorAll?.('path') || []).forEach((path) => {
    if (!String(path.getAttribute?.('fill') || '').startsWith('url(')) return;
    const numbers = String(path.getAttribute?.('d') || '').match(/[+-]?(?:\d+(?:\.\d*)?|\.\d+)/g)
      ?.map(Number) || [];
    const xs = numbers.filter((_value, index) => index % 2 === 0);
    if (xs.length > 0) width = Math.max(width, Math.max(...xs) - Math.min(...xs));
  });
  return width;
};

const rendererLegendSlots = ({
  targetGroup,
  initialGroups,
  desiredGroups,
  rendererLayout
}) => {
  const metrics = rendererLayout.metadata.legendReflow;
  const rectSize = metrics.colorRectSize;
  const lineHeight = metrics.lineHeight;
  const textXOffset = metrics.textXOffset;
  const orientation = legendOrientation(targetGroup, [], rendererLayout.orientation);
  const initialEntries = initialGroups.map((group) => {
    const anchor = legendEntryAnchor(group);
    const label = group.querySelector?.('text');
    const textWidth = measuredLegendTextWidth(label);
    return {
      anchor,
      group,
      left: anchor.x - textXOffset,
      right: anchor.x + rectSize + textWidth + textXOffset,
      width: rectSize + textWidth + textXOffset * 2
    };
  });
  const initialRows = sortedLegendRows(initialEntries.map(({ anchor }) => anchor));
  const slots = initialRows.flatMap((row) => row.anchors).slice(0, desiredGroups.length);
  if (orientation === 'vertical') {
    const x = slots[0]?.x ?? textXOffset;
    while (slots.length < desiredGroups.length) {
      slots.push({
        x,
        y: (slots.at(-1)?.y ?? (rectSize / 2 - lineHeight)) + lineHeight
      });
    }
    return { orientation, topology: 'vertical', slots };
  }

  let availableWidth = rendererLayout.metadata.primary.finalBounds.width;
  const reservedWidth = comparisonLegendWidth(targetGroup);
  if (reservedWidth > 0) {
    availableWidth = Math.max(
      rectSize + textXOffset * 2,
      availableWidth - reservedWidth - textXOffset
    );
  }
  const visualEntries = [...initialEntries].sort((left, right) => {
    const yDelta = left.anchor.y - right.anchor.y;
    return Math.abs(yDelta) < LEGEND_ROW_TOLERANCE
      ? left.anchor.x - right.anchor.x
      : yDelta;
  });
  const rows = [];
  visualEntries.forEach((entry) => {
    const row = rows.at(-1);
    if (!row || Math.abs(row.y - entry.anchor.y) >= LEGEND_ROW_TOLERANCE) {
      rows.push({ y: entry.anchor.y, entries: [entry], right: entry.right });
    } else {
      row.entries.push(entry);
      row.right = Math.max(row.right, entry.right);
    }
  });

  let nextGroupIndex = slots.length;
  let activeRow = rows.at(-1) || {
    y: rectSize / 2,
    entries: [],
    right: 0
  };
  while (nextGroupIndex < desiredGroups.length) {
    const group = desiredGroups[nextGroupIndex];
    const width = rectSize
      + measuredLegendTextWidth(group.querySelector?.('text'))
      + textXOffset * 2;
    const left = activeRow.entries.length > 0 ? activeRow.right : 0;
    if (activeRow.entries.length > 0 && left + width > availableWidth) break;
    const anchor = { x: left + textXOffset, y: activeRow.y };
    slots.push(anchor);
    activeRow.entries.push({ anchor, group, left, right: left + width, width });
    activeRow.right = left + width;
    nextGroupIndex += 1;
  }

  const pendingRows = [];
  let pendingRow = { entries: [], width: 0 };
  for (; nextGroupIndex < desiredGroups.length; nextGroupIndex += 1) {
    const group = desiredGroups[nextGroupIndex];
    const width = rectSize
      + measuredLegendTextWidth(group.querySelector?.('text'))
      + textXOffset * 2;
    if (pendingRow.entries.length > 0 && pendingRow.width + width > availableWidth) {
      pendingRows.push(pendingRow);
      pendingRow = { entries: [], width: 0 };
    }
    pendingRow.entries.push({ group, width });
    pendingRow.width += width;
  }
  if (pendingRow.entries.length > 0) pendingRows.push(pendingRow);
  const lastY = rows.at(-1)?.y ?? (rectSize / 2 - lineHeight);
  pendingRows.forEach((row, rowIndex) => {
    let cursorX = Math.max(0, (availableWidth - row.width) * 0.5);
    const y = lastY + (rowIndex + 1) * lineHeight;
    row.entries.forEach((entry) => {
      slots.push({ x: cursorX + textXOffset, y });
      cursorX += entry.width;
    });
  });
  return {
    orientation,
    topology: rows.length + pendingRows.length > 1 ? 'horizontal-wrapped' : 'horizontal',
    slots
  };
};

const unionMeasuredBounds = (bounds) => {
  const usable = bounds.filter(Boolean);
  if (usable.length === 0) return null;
  const x = Math.min(...usable.map((item) => item.x));
  const y = Math.min(...usable.map((item) => item.y));
  const right = Math.max(...usable.map((item) => item.x + item.width));
  const bottom = Math.max(...usable.map((item) => item.y + item.height));
  return { x, y, width: right - x, height: bottom - y };
};

const numericAttribute = (element, name, fallback = 0) => {
  const value = Number.parseFloat(element?.getAttribute?.(name));
  return Number.isFinite(value) ? value : fallback;
};

const translatedBounds = (bounds, x, y) => bounds ? ({
  x: bounds.x + x,
  y: bounds.y + y,
  width: bounds.width,
  height: bounds.height
}) : null;

const expandedForStroke = (element, bounds) => {
  if (!bounds || text(element?.getAttribute?.('stroke')).toLowerCase() === 'none') return bounds;
  const halfWidth = Math.max(0, numericAttribute(element, 'stroke-width')) * 0.5;
  return {
    x: bounds.x - halfWidth,
    y: bounds.y - halfWidth,
    width: bounds.width + halfWidth * 2,
    height: bounds.height + halfWidth * 2
  };
};

const pathBounds = (path) => {
  const numbers = String(path?.getAttribute?.('d') || '')
    .match(/[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?/g)
    ?.map(Number) || [];
  if (numbers.length < 2) return null;
  const xs = numbers.filter((_value, index) => index % 2 === 0);
  const ys = numbers.filter((_value, index) => index % 2 === 1);
  if (xs.length === 0 || ys.length === 0) return null;
  const x = Math.min(...xs);
  const y = Math.min(...ys);
  return {
    x,
    y,
    width: Math.max(...xs) - x,
    height: Math.max(...ys) - y
  };
};

const ownDetachedBounds = (element) => {
  const tagName = text(element?.tagName).toLowerCase().replace(/^.*:/, '');
  if (tagName === 'path') {
    const fill = text(element.getAttribute?.('fill')).toLowerCase();
    const stroke = text(element.getAttribute?.('stroke')).toLowerCase();
    if ((fill === '' || fill === 'none') && (stroke === '' || stroke === 'none')) return null;
    return expandedForStroke(element, pathBounds(element));
  }
  if (tagName === 'rect') {
    return expandedForStroke(element, {
      x: numericAttribute(element, 'x'),
      y: numericAttribute(element, 'y'),
      width: Math.max(0, numericAttribute(element, 'width')),
      height: Math.max(0, numericAttribute(element, 'height'))
    });
  }
  if (tagName === 'line') {
    const x1 = numericAttribute(element, 'x1');
    const x2 = numericAttribute(element, 'x2');
    const y1 = numericAttribute(element, 'y1');
    const y2 = numericAttribute(element, 'y2');
    return expandedForStroke(element, {
      x: Math.min(x1, x2),
      y: Math.min(y1, y2),
      width: Math.abs(x2 - x1),
      height: Math.abs(y2 - y1)
    });
  }
  if (tagName === 'circle') {
    const radius = Math.max(0, numericAttribute(element, 'r'));
    return expandedForStroke(element, {
      x: numericAttribute(element, 'cx') - radius,
      y: numericAttribute(element, 'cy') - radius,
      width: radius * 2,
      height: radius * 2
    });
  }
  if (tagName !== 'text') return null;
  const fontSize = Math.max(1, numericAttribute(element, 'font-size', 16));
  const width = measuredLegendTextWidth(element);
  const anchor = text(element.getAttribute?.('text-anchor')).toLowerCase();
  const baseline = text(element.getAttribute?.('dominant-baseline')).toLowerCase();
  const x = numericAttribute(element, 'x');
  const y = numericAttribute(element, 'y');
  const left = anchor === 'middle' ? x - width / 2 : anchor === 'end' ? x - width : x;
  const top = baseline === 'hanging'
    ? y
    : (baseline === 'central' || baseline === 'middle') ? y - fontSize / 2 : y - fontSize;
  return expandedForStroke(element, { x: left, y: top, width, height: fontSize });
};

const isComparisonLegend = (element) => (
  element?.getAttribute?.('data-gbdraw-role') === 'comparison-legend'
  || /^pairwise_legend(?:_[hv])?$/.test(text(element?.id))
  || text(element?.id) === 'conservation_identity_legend'
);

const detachedLegendBounds = (
  root,
  { excludeComparison = false, includeRootTransform = false } = {}
) => {
  const visit = (element, offsetX, offsetY, includeTransform = true) => {
    if (!element) return null;
    if (
      element.getAttribute?.('display') === 'none'
      || element.style?.display === 'none'
      || (excludeComparison && isComparisonLegend(element))
    ) return null;
    const tagName = text(element.tagName).toLowerCase().replace(/^.*:/, '');
    if (tagName === 'defs' || tagName === 'lineargradient' || tagName === 'stop') return null;
    const transform = includeTransform
      ? parseTransformXY(element.getAttribute?.('transform'))
      : { x: 0, y: 0 };
    const nextX = offsetX + transform.x;
    const nextY = offsetY + transform.y;
    const own = translatedBounds(ownDetachedBounds(element), nextX, nextY);
    const children = Array.from(element.children || []).map(
      (child) => visit(child, nextX, nextY)
    );
    return unionMeasuredBounds([own, ...children]);
  };
  return visit(root, 0, 0, includeRootTransform);
};

const parsedViewBox = (value) => {
  const parts = String(value || '').trim().split(/[\s,]+/).map(Number);
  return parts.length === 4 && parts.every(Number.isFinite) ? parts : null;
};

const projectedCanvasPadding = (svg) => {
  const original = parsedViewBox(svg?.getAttribute?.('data-original-view-box'));
  const current = parsedViewBox(svg?.getAttribute?.('viewBox'));
  if (!original || !current) return null;
  const [baseX, baseY, baseWidth, baseHeight] = original;
  const [currentX, currentY, currentWidth, currentHeight] = current;
  return {
    top: Math.max(0, baseY - currentY),
    right: Math.max(0, currentX + currentWidth - (baseX + baseWidth)),
    bottom: Math.max(0, currentY + currentHeight - (baseY + baseHeight)),
    left: Math.max(0, baseX - currentX)
  };
};

const boundsDiffer = (left, right) => (
  !left
  || !right
  || ['x', 'y', 'width', 'height'].some(
    (field) => Math.abs(Number(left[field]) - Number(right[field])) > 1e-6
  )
);

const reconcileProjectedLegendComposition = (svg, mutation) => {
  const rendererLayout = rendererLegendLayout(svg);
  if (!rendererLayout) return 0;
  const legend = svg.getElementById?.('legend');
  if (!legend) return 0;
  const { lineHeight, textXOffset } = rendererLayout.metadata.legendReflow;
  let changed = 0;

  getAllFeatureLegendGroups(svg).forEach((targetGroup) => {
    const orientation = legendOrientation(targetGroup, [], rendererLayout.orientation);
    const owner = targetGroup === legend ? legend : targetGroup.parentElement;
    const comparison = getComparisonLegendGroup(owner);
    if (!comparison) return;
    const featureBounds = detachedLegendBounds(targetGroup, {
      excludeComparison: true,
      includeRootTransform: targetGroup !== legend
    });
    const comparisonBounds = detachedLegendBounds(comparison);
    if (!featureBounds || !comparisonBounds) return;
    const comparisonX = orientation === 'horizontal'
      ? featureBounds.x + featureBounds.width + textXOffset - comparisonBounds.x
      : featureBounds.x + featureBounds.width / 2
        - (comparisonBounds.x + comparisonBounds.width / 2);
    const comparisonY = orientation === 'horizontal'
      ? featureBounds.y + featureBounds.height / 2
        - (comparisonBounds.y + comparisonBounds.height / 2)
      : featureBounds.y + featureBounds.height + lineHeight / 2 - comparisonBounds.y;
    if (mutation.setAttribute(
      comparison,
      'transform',
      `translate(${comparisonX}, ${comparisonY})`
    )) changed += 1;
  });

  const localBounds = detachedLegendBounds(legend);
  if (!localBounds) return changed;
  const envelope = Array.from(legend.children || []).find((child) => (
    text(child.tagName).toLowerCase().replace(/^.*:/, '') === 'path'
    && ['none', ''].includes(text(child.getAttribute?.('fill')).toLowerCase())
    && ['none', ''].includes(text(child.getAttribute?.('stroke')).toLowerCase())
  ));
  if (envelope) {
    const right = localBounds.x + localBounds.width;
    const bottom = localBounds.y + localBounds.height;
    if (mutation.setAttribute(
      envelope,
      'd',
      `M ${localBounds.x},${localBounds.y} L ${right},${localBounds.y} L ${right},${bottom} L ${localBounds.x},${bottom} z`
    )) changed += 1;
  }
  if (!boundsDiffer(localBounds, rendererLayout.metadata.legend.localBounds)) return changed;
  applyCompositionEdit(svg, {
    legendLocalBounds: localBounds,
    titleLocalBounds: rendererLayout.metadata.title?.localBounds || null,
    primaryFinalBounds: rendererLayout.metadata.primary.finalBounds,
    canvasPadding: projectedCanvasPadding(svg),
    mutation
  });
  return changed + 1;
};

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
  const rendererLayout = rendererLegendLayout(svg);
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

    const rendererSlots = rendererLayout
      ? rendererLegendSlots({
          targetGroup,
          initialGroups,
          desiredGroups,
          rendererLayout
        })
      : null;
    const slots = rendererSlots?.slots || resolveLegendSlots({
        targetGroup,
        // Removed entries still own renderer-authored slots. Resolve topology from
        // the complete admitted Legend before projecting the desired structure.
        initialGroups,
        desiredCount: desiredGroups.length
      });
    if (rendererSlots) {
      recordStructuralMetric('legendTopologyResolutionCount', 1, {
        topology: rendererSlots.topology,
        orientation: rendererSlots.orientation,
        initialCount: initialGroups.length,
        desiredCount: desiredGroups.length,
        targetGroupOwner: targetGroup
      });
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
  recordStructuralMetric('editorProjectionBuildCount');
  recordStructuralMetric('editorProjectionBorrowedFeatureOwnerCount', 1, {
    featureOwner: state.features
  });
  replaceRuleDerivedFillOverrides(
    state.featureColorOverrides,
    state.features,
    input.manualSpecificRules
  );

  const project = (svg, {
    item = null,
    requireFeatureBindings = false,
    resetFeatureVisibility = false,
    resetLabelState = false,
    mutation: suppliedMutation = null
  } = {}) => {
    if (!svg) throw new Error('Editor SVG projection requires an SVG root.');
    recordStructuralMetric('editorProjectionApplicationCount');
    const elementsById = getFeatureElementIndex(svg);
    const bindings = featureBindings(state.features, item);
    recordStructuralMetric('editorProjectionFeatureBindingAccessCount', bindings.length, {
      featureOwner: state.features,
      itemOwner: item
    });
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

    const mutation = suppliedMutation || createDomMutationJournal();
    const ownsMutation = !suppliedMutation;
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
      const visibility = resolveEffectiveFeatureVisibility(
        feature || renderedId,
        state.featureVisibilityOverrides,
        null,
        state.featureVisibilityManualRules,
        { overrideKeys: [renderedId, overrideKey] }
      );
      elements.forEach((element) => {
        let changed = false;
        if (resetFeatureVisibility) {
          changed = applyFeatureVisibility(element, 'default', {
            markPreview: true,
            mutation,
            reason: 'editor-projection-reset'
          });
        }
        changed = applyFeatureVisibility(element, visibility, {
          markPreview: true,
          mutation,
          reason: 'editor-projection'
        }) || changed;
        if (changed) featureVisibilityCount += 1;
      });
      });

      const labelResult = labelProjection(svg, state, { resetLabelState, mutation });
      const projectedLegendEntries = state.legendEntries.map((entry) => ({
        ...entry,
        color: hasOwn(state.legendColorOverrides, entry?.caption)
          ? state.legendColorOverrides[entry.caption]
          : null
      }));
      const suppressionCount = state.suppressPairwiseIdentityLegend
        ? removePairwiseIdentityLegend(svg, mutation)
        : 0;
      const legendStructureCount = applyLegendStructureToSvg({
        svg,
        legendEntries: state.legendEntries,
        deletedLegendEntries: state.deletedLegendEntries,
        originalLegendOrder: state.originalLegendOrder,
        addedLegendCaptions: state.addedLegendCaptions,
        mutation
      });
      const legendStyleCount = applyLegendEntryStateToSvg({
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
      let legendCount = legendStructureCount + legendStyleCount;
      if (legendStructureCount > 0 || suppressionCount > 0) {
        legendCount += reconcileProjectedLegendComposition(svg, mutation);
      }
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
      if (ownsMutation) mutation.commit();
      return result;
    } catch (error) {
      if (ownsMutation) mutation.rollback();
      throw error;
    }
  };

  return {
    featureColorOverrides: state.featureColorOverrides,
    featureStrokeOverrides: state.featureStrokeOverrides,
    project
  };
};
