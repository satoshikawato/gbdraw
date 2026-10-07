// @ts-check
import {
  FEATURE_SELECTOR,
  getFeatureElementIndex,
  getFeatureElements,
  getFeatureIdentity
} from '../feature-editor/svg-actions.js';
import { getAllFeatureLegendGroups, mountedLegendRowFeatureIds, setsFeatureStroke } from './utils.js';
import {
  featureOverrideKey,
  migrateLegacyFeatureOverrides
} from '../../services/feature-override-identity.js';

const hasOwn = (value, key) => Object.prototype.hasOwnProperty.call(value || {}, key);

const applyStrokeAttributes = (element, strokeColor, strokeWidth) => {
  let changed = false;
  const color = String(strokeColor || '').trim();
  if (color && element.getAttribute('stroke') !== color) {
    element.setAttribute('stroke', color);
    changed = true;
  }
  if (strokeWidth !== undefined && strokeWidth !== null && strokeWidth !== '') {
    const numericWidth = Number(strokeWidth);
    if (!Number.isFinite(numericWidth) || numericWidth < 0) return changed;
    const width = String(numericWidth);
    if (element.getAttribute('stroke-width') !== width) {
      element.setAttribute('stroke-width', width);
      changed = true;
    }
  }
  return changed;
};

const restoreStrokeAttributes = (element, originalColor, originalWidth) => {
  let changed = false;
  const nextColor = originalColor === null ? null : String(originalColor);
  const nextWidth = originalWidth === null ? null : String(originalWidth);
  if (nextColor === null) {
    if (element.hasAttribute('stroke')) {
      element.removeAttribute('stroke');
      changed = true;
    }
  } else if (element.getAttribute('stroke') !== nextColor) {
    element.setAttribute('stroke', nextColor);
    changed = true;
  }
  if (nextWidth === null) {
    if (element.hasAttribute('stroke-width')) {
      element.removeAttribute('stroke-width');
      changed = true;
    }
  } else if (element.getAttribute('stroke-width') !== nextWidth) {
    element.setAttribute('stroke-width', nextWidth);
    changed = true;
  }
  return changed;
};

const getLegendSwatches = (svg, caption) => {
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

// The feature edits `legendRowFeatureIds` reads besides the drawn fills, by
// rendered ID: the features a feature color edit names into each Legend row,
// and the features with a stroke edit of their own (OV-123).
/**
 * @param {Array<Record<string, any>>} features
 * @param {{ featureColorOverrides?: Record<string, any>, featureStrokeOverrides?: Record<string, any> }} overrides
 */
const featureEditFacts = (features, { featureColorOverrides = {}, featureStrokeOverrides = {} }) => {
  /** @type {Map<string, string[]>} */
  const namedIdsByCaption = new Map();
  /** @type {string[]} */
  const ownStrokeIds = [];
  (Array.isArray(features) ? features : []).forEach((feature) => {
    const key = featureOverrideKey(feature);
    const svgId = String(feature?.svg_id || '').trim();
    if (!key || !svgId) return;
    const caption = String(featureColorOverrides?.[key]?.caption || '').trim();
    if (caption) namedIdsByCaption.set(caption, [...(namedIdsByCaption.get(caption) || []), svgId]);
    if (setsFeatureStroke(featureStrokeOverrides?.[key])) ownStrokeIds.push(svgId);
  });
  /** @param {string} caption */
  const ofRow = (caption) => ({ namedIds: namedIdsByCaption.get(caption) || [], ownStrokeIds });
  return { ownStrokeIds, ofRow };
};

/** @param {{ svg?: Element | null, legendColorOverrides?: Record<string, string> }} [options] */
export const applyLegendColorOverridesToSvg = ({
  svg,
  legendColorOverrides = {}
} = {}) => {
  if (!svg) return 0;
  let changedCount = 0;
  Object.entries(legendColorOverrides || {}).forEach(([caption, color]) => {
    const normalized = String(color || '').trim();
    if (!normalized) return;
    getLegendSwatches(svg, caption).forEach((swatch) => {
      if (swatch.getAttribute('fill') === normalized) return;
      swatch.setAttribute('fill', normalized);
      changedCount += 1;
    });
  });
  return changedCount;
};

/**
 * @param {{
 *   svg?: Element | null,
 *   features?: Array<Record<string, any>>,
 *   legendEntries?: Array<Record<string, any>>,
 *   legendStrokeOverrides?: Record<string, any>,
 *   featureColorOverrides?: Record<string, any>,
 *   featureStrokeOverrides?: Record<string, any>
 * }} [options]
 */
export const applyStrokeOverridesToSvg = ({
  svg,
  features = [],
  legendEntries = [],
  legendStrokeOverrides = {},
  featureColorOverrides = {},
  featureStrokeOverrides = {}
} = {}) => {
  if (!svg) return 0;
  const featureIndex = getFeatureElementIndex(svg);
  const edits = featureEditFacts(features, { featureColorOverrides, featureStrokeOverrides });
  let changedCount = 0;
  /** @param {string} svgId @param {Record<string, any>} overrides */
  const applyToFeature = (svgId, overrides) => {
    if (!overrides || typeof overrides !== 'object') return;
    getFeatureElements(svg, svgId, featureIndex).forEach((element) => {
      if (applyStrokeAttributes(element, overrides.strokeColor, overrides.strokeWidth)) {
        changedCount += 1;
      }
    });
  };

  Object.entries(legendStrokeOverrides || {}).forEach(([caption, overrides]) => {
    if (!overrides || typeof overrides !== 'object') return;
    mountedLegendRowFeatureIds(svg, caption, legendEntries, edits.ofRow(caption))
      .forEach((svgId) => applyToFeature(svgId, overrides));
    getLegendSwatches(svg, caption).forEach((swatch) => {
      if (applyStrokeAttributes(swatch, overrides.strokeColor, overrides.strokeWidth)) {
        changedCount += 1;
      }
    });
  });

  (Array.isArray(features) ? features : []).forEach((feature) => {
    const key = featureOverrideKey(feature);
    if (key) applyToFeature(feature?.svg_id, featureStrokeOverrides?.[key]);
  });
  return changedCount;
};

/**
 * @typedef {object} LegendStrokeActionsOptions
 * @property {Record<string, any>} state App state (state.js; not yet typed).
 * @property {((reason: string) => boolean) | null} [commitActiveResultEdit]
 *   The preview owner's commit of an edit to the displayed Result (R1, R13).
 */

/** @param {LegendStrokeActionsOptions} options */
export const createLegendStrokeActions = ({ state, commitActiveResultEdit = null }) => {
  const {
    extractedFeatures,
    legendEntries,
    legendStrokeOverrides,
    featureColorOverrides,
    featureStrokeOverrides,
    originalSvgStroke,
    svgContainer,
    legendStrokeOptionsOpen
  } = state;

  // OV-157: a row's Stroke options button shows or hides its stroke controls.
  // That is view state, kept out of the Legend entries, so the click records
  // no History step and the Session does not save it.
  /** @param {string} caption */
  const isLegendStrokeOptionsOpen = (caption) => Boolean(legendStrokeOptionsOpen?.has(caption));
  /** @param {string} caption */
  const toggleLegendStrokeOptions = (caption) => {
    if (!legendStrokeOptionsOpen?.delete(caption)) legendStrokeOptionsOpen?.add(caption);
  };
  const closeLegendStrokeOptions = () => legendStrokeOptionsOpen?.clear();

  const captureLegendSwatchStroke = (caption) => {
    const svg = svgContainer.value?.querySelector?.('svg');
    const swatch = svg ? getLegendSwatches(svg, caption)[0] : null;
    const widthValue = swatch?.getAttribute?.('stroke-width');
    const width = widthValue === null || widthValue === undefined || widthValue === ''
      ? null
      : Number(widthValue);
    return {
      originalStrokeColor: swatch?.getAttribute?.('stroke') ?? null,
      originalStrokeWidth: width !== null && Number.isFinite(width) ? width : null
    };
  };

  const persistStrokeEdit = (reason) => commitActiveResultEdit?.(reason);

  const liveFeatureEdits = () => featureEditFacts(extractedFeatures.value, { featureColorOverrides, featureStrokeOverrides });
  // The features of the mounted Result that a stroke on the Legend row
  // `caption` reaches (`legendRowFeatureIds`, OV-123).
  /** @param {Element} svg @param {string} caption */
  const rowFeatureIds = (svg, caption) => mountedLegendRowFeatureIds(
    svg, caption, legendEntries.value, liveFeatureEdits().ofRow(caption)
  );

  // The stroke the renderer gives a feature block, read from the first feature
  // path. Generate applies the stroke edits before the Result is mounted, so a
  // path that an edit reached gives the stroke the edit recorded it replaced:
  // its Legend row's (the row's swatch as drawn), else its own (OV-123).
  const captureOriginalStroke = () => {
    const svg = svgContainer.value?.querySelector?.('svg');
    const firstFeaturePath = svg?.querySelector?.('path[id^="f"]');
    if (!svg || !firstFeaturePath) return;
    const svgId = getFeatureIdentity(firstFeaturePath);
    const edits = liveFeatureEdits();
    const recorded = [
      ...Object.entries(legendStrokeOverrides)
        .filter(([caption]) => mountedLegendRowFeatureIds(
          svg, caption, legendEntries.value, { namedIds: edits.ofRow(caption).namedIds }
        ).includes(svgId))
        .map(([, override]) => override),
      ...extractedFeatures.value
        .filter((feature) => String(feature?.svg_id || '').trim() === svgId)
        .map((feature) => featureStrokeOverrides[featureOverrideKey(feature)])
    ].find((override) => hasOwn(override, 'originalStrokeColor'));
    const widthValue = recorded && hasOwn(recorded, 'originalStrokeWidth')
      ? recorded.originalStrokeWidth
      : firstFeaturePath.getAttribute('stroke-width');
    const strokeWidth = Number.parseFloat(String(widthValue ?? ''));
    originalSvgStroke.value = {
      color: recorded ? recorded.originalStrokeColor : firstFeaturePath.getAttribute('stroke'),
      width: Number.isFinite(strokeWidth) ? strokeWidth : null
    };
  };

  const getLegendEntryStrokeColor = (idx) => {
    const entry = legendEntries.value[idx];
    if (!entry) return '';
    const override = legendStrokeOverrides[entry.caption];
    if (override && override.strokeColor !== undefined) return override.strokeColor;
    return '';
  };

  const getLegendEntryStrokeWidth = (idx) => {
    const entry = legendEntries.value[idx];
    if (!entry) return '';
    const override = legendStrokeOverrides[entry.caption];
    if (override && override.strokeWidth !== undefined) return override.strokeWidth;
    return '';
  };

  const updateLegendEntryStrokeColor = (idx, color) => {
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    const entry = legendEntries.value[idx];
    if (!entry) return false;
    const normalized = String(color || '').trim();
    if (String(legendStrokeOverrides[entry.caption]?.strokeColor || '').trim() === normalized) {
      return false;
    }

    if (!legendStrokeOverrides[entry.caption]) {
      legendStrokeOverrides[entry.caption] = {
        ...captureLegendSwatchStroke(entry.caption),
        strokeColor: normalized,
        strokeWidth: getLegendEntryStrokeWidth(idx)
      };
    }
    legendStrokeOverrides[entry.caption].strokeColor = normalized;

    applyStrokeToFeaturesByCaption(entry.caption, normalized, null);
    return true;
  };

  const updateLegendEntryStrokeWidth = (idx, width) => {
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    const entry = legendEntries.value[idx];
    if (!entry) return false;

    const widthVal = parseFloat(width);
    if (isNaN(widthVal)) return false;
    if (Number(legendStrokeOverrides[entry.caption]?.strokeWidth) === widthVal) return false;

    if (!legendStrokeOverrides[entry.caption]) {
      legendStrokeOverrides[entry.caption] = {
        ...captureLegendSwatchStroke(entry.caption),
        strokeColor: getLegendEntryStrokeColor(idx),
        strokeWidth: widthVal
      };
    }
    legendStrokeOverrides[entry.caption].strokeWidth = widthVal;

    applyStrokeToFeaturesByCaption(entry.caption, null, widthVal);
    return true;
  };

  const setLegendEntryStrokeColorValue = (idx, value) => {
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    const entry = legendEntries.value[idx];
    if (!entry) return;
    if (value !== null) {
      return updateLegendEntryStrokeColor(idx, String(value || '').trim());
    }
    const override = legendStrokeOverrides[entry.caption];
    if (!override || !Object.prototype.hasOwnProperty.call(override, 'strokeColor')) return false;
    if (override) {
      delete override.strokeColor;
      if (override.strokeWidth === undefined || override.strokeWidth === '') {
        delete legendStrokeOverrides[entry.caption];
      }
    }
    const inheritedColor = originalSvgStroke.value.color;
    applyStrokeToFeaturesByCaption(entry.caption, inheritedColor, null, {
      removeStroke: inheritedColor === null
    });
    return true;
  };

  const resetLegendEntryStroke = (idx) => {
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    const entry = legendEntries.value[idx];
    if (!entry) return false;
    if (!svgContainer.value) return false;

    const svg = svgContainer.value.querySelector('svg');
    if (!svg) return false;

    const override = legendStrokeOverrides[entry.caption];
    const originalColor = originalSvgStroke.value.color;
    const originalWidth = originalSvgStroke.value.width;
    const originalSwatchColor = hasOwn(override, 'originalStrokeColor')
      ? override.originalStrokeColor
      : originalColor;
    const originalSwatchWidth = hasOwn(override, 'originalStrokeWidth')
      ? override.originalStrokeWidth
      : originalWidth;
    let updatedCount = 0;

    const featureIndex = getFeatureElementIndex(svg);
    rowFeatureIds(svg, entry.caption).forEach((svgId) => {
      getFeatureElements(svg, svgId, featureIndex).forEach((el) => {
        if (restoreStrokeAttributes(el, originalColor, originalWidth)) updatedCount++;
      });
    });

    getLegendSwatches(svg, entry.caption).forEach((swatch) => {
      if (restoreStrokeAttributes(swatch, originalSwatchColor, originalSwatchWidth)) {
        updatedCount++;
      }
    });

    const overrideRemoved = Object.prototype.hasOwnProperty.call(
      legendStrokeOverrides,
      entry.caption
    );
    if (overrideRemoved) delete legendStrokeOverrides[entry.caption];

    if (updatedCount > 0) {
      persistStrokeEdit('reset-legend-stroke');
    }
    return updatedCount > 0 || overrideRemoved;
  };

  const resetAllStrokes = () => {
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    if (!svgContainer.value) return false;
    const svg = svgContainer.value.querySelector('svg');
    if (!svg) return false;

    const originalColor = originalSvgStroke.value.color;
    const originalWidth = originalSvgStroke.value.width;

    const featurePaths = svg.querySelectorAll(FEATURE_SELECTOR);
    let updatedCount = 0;
    featurePaths.forEach((path) => {
      if (restoreStrokeAttributes(path, originalColor, originalWidth)) updatedCount++;
    });

    const legendGroups = getAllFeatureLegendGroups(svg);
    for (const targetGroup of legendGroups) {
      const paths = targetGroup.querySelectorAll('path');
      paths.forEach((p) => {
        const fill = p.getAttribute('fill');
        if (fill && fill !== 'none' && !fill.startsWith('url(')) {
          if (restoreStrokeAttributes(p, originalColor, originalWidth)) updatedCount++;
        }
      });
    }

    const overridesRemoved =
      Object.keys(legendStrokeOverrides).length > 0 ||
      Object.keys(featureStrokeOverrides).length > 0;
    Object.keys(legendStrokeOverrides).forEach((key) => delete legendStrokeOverrides[key]);
    Object.keys(featureStrokeOverrides).forEach((key) => delete featureStrokeOverrides[key]);

    if (updatedCount > 0) {
      persistStrokeEdit('reset-all-strokes');
      console.log(
        `Reset all strokes: updated ${updatedCount} elements to original (color=${originalColor}, width=${originalWidth})`
      );
    }
    return updatedCount > 0 || overridesRemoved;
  };

  const applyStrokeToFeaturesByCaption = (
    caption,
    strokeColor,
    strokeWidth,
    { removeStroke = false } = {}
  ) => {
    if (!svgContainer.value) return;
    const svg = svgContainer.value.querySelector('svg');
    if (!svg) return;

    let updatedCount = 0;

    const featureIndex = getFeatureElementIndex(svg);
    const reached = rowFeatureIds(svg, caption);
    reached.forEach((svgId) => {
      getFeatureElements(svg, svgId, featureIndex).forEach((el) => {
        if (removeStroke) {
          if (el.getAttribute('stroke') !== null) {
            el.removeAttribute('stroke');
            updatedCount++;
          }
        } else if (strokeColor !== null) {
          if (el.getAttribute('stroke') !== String(strokeColor)) {
            el.setAttribute('stroke', strokeColor);
            updatedCount++;
          }
        }
        if (strokeWidth !== null && el.getAttribute('stroke-width') !== String(strokeWidth)) {
          el.setAttribute('stroke-width', strokeWidth);
          updatedCount++;
        }
      });
    });
    console.log(`Applied stroke to ${reached.length} features for "${caption}"`);

    const legendGroups = getAllFeatureLegendGroups(svg);
    for (const targetGroup of legendGroups) {
      const entryGroup = targetGroup.querySelector(`g[data-legend-key="${CSS.escape(caption)}"]`);
      if (entryGroup) {
        const paths = entryGroup.querySelectorAll('path');
        for (const path of paths) {
          const fill = path.getAttribute('fill');
          if (fill && fill !== 'none' && !fill.startsWith('url(')) {
            if (removeStroke) {
              if (path.getAttribute('stroke') !== null) {
                path.removeAttribute('stroke');
                updatedCount++;
              }
            } else if (strokeColor !== null) {
              if (path.getAttribute('stroke') !== String(strokeColor)) {
                path.setAttribute('stroke', strokeColor);
                updatedCount++;
              }
            }
            if (strokeWidth !== null && path.getAttribute('stroke-width') !== String(strokeWidth)) {
              path.setAttribute('stroke-width', strokeWidth);
              updatedCount++;
            }
            break;
          }
        }
      }
    }

    if (updatedCount > 0) {
      persistStrokeEdit('legend-stroke');
      console.log(`Applied stroke to ${updatedCount} elements for caption "${caption}"`);
    }
  };

  const reapplyStrokeOverrides = () => {
    if (!svgContainer.value) return false;
    const svg = svgContainer.value.querySelector('svg');
    if (!svg) return false;

    if (
      Object.keys(legendStrokeOverrides).length === 0 &&
      Object.keys(featureStrokeOverrides).length === 0
    ) return false;
    migrateLegacyFeatureOverrides(featureStrokeOverrides, extractedFeatures.value);
    const totalUpdated = applyStrokeOverridesToSvg({
      svg,
      features: extractedFeatures.value,
      legendEntries: legendEntries.value,
      legendStrokeOverrides,
      featureColorOverrides,
      featureStrokeOverrides
    });

    if (totalUpdated > 0) {
      persistStrokeEdit('reapply-stroke-overrides');
    }
    return totalUpdated > 0;
  };

  /** @param {{ changes?: unknown }} [options] `changes` is a History change list; anything else is ignored. */
  const reconcileStrokeOverrides = ({ changes = null } = {}) => {
    const svg = svgContainer.value?.querySelector?.('svg');
    if (!svg) return false;
    const originalColor = originalSvgStroke.value.color;
    const originalWidth = originalSvgStroke.value.width;
    let changed = false;
    const historyChanges = Array.isArray(changes) ? changes : null;

    if (historyChanges) {
      const featureBaselines = new Map();
      const legendBaselines = new Map();
      const collectBaseline = (target, key, value) => {
        if (!key) return;
        const baseline = target.get(key) || {};
        if (hasOwn(value, 'originalStrokeColor')) {
          baseline.originalStrokeColor = value.originalStrokeColor;
        }
        if (hasOwn(value, 'originalStrokeWidth')) {
          baseline.originalStrokeWidth = value.originalStrokeWidth;
        }
        target.set(key, baseline);
      };

      historyChanges.forEach((change) => {
        const path = Array.isArray(change?.path) ? change.path : [];
        const isFeatureStroke = path[0] === 'editorState'
          && path[1] === 'featureStrokes'
          && path[2] === 'overrides';
        const isLegendStroke = path[0] === 'editorState'
          && path[1] === 'legend'
          && path[2] === 'strokeOverrides';
        if (!isFeatureStroke && !isLegendStroke) return;
        const key = String(path[3] || '').trim();
        const target = isFeatureStroke ? featureBaselines : legendBaselines;
        collectBaseline(target, key, change.before);
        collectBaseline(target, key, change.after);
        collectBaseline(
          target,
          key,
          isFeatureStroke ? featureStrokeOverrides[key] : legendStrokeOverrides[key]
        );
      });

      if (featureBaselines.size === 0 && legendBaselines.size === 0) {
        return reapplyStrokeOverrides();
      }

      legendBaselines.forEach((baseline, caption) => {
        const baselineColor = hasOwn(baseline, 'originalStrokeColor')
          ? baseline.originalStrokeColor
          : originalColor;
        const baselineWidth = hasOwn(baseline, 'originalStrokeWidth')
          ? baseline.originalStrokeWidth
          : originalWidth;
        const featureIndex = getFeatureElementIndex(svg);
        rowFeatureIds(svg, caption).forEach((featureId) => {
          getFeatureElements(svg, featureId, featureIndex).forEach((element) => {
            if (restoreStrokeAttributes(element, originalColor, originalWidth)) changed = true;
          });
        });
        getLegendSwatches(svg, caption).forEach((swatch) => {
          if (restoreStrokeAttributes(swatch, baselineColor, baselineWidth)) changed = true;
        });
      });

      const featuresByOverrideKey = new Map();
      extractedFeatures.value.forEach((feature) => {
        const key = featureOverrideKey(feature);
        if (key) featuresByOverrideKey.set(key, feature);
        const svgId = String(feature?.svg_id || '').trim();
        if (svgId) featuresByOverrideKey.set(svgId, feature);
      });
      featureBaselines.forEach((baseline, key) => {
        const feature = featuresByOverrideKey.get(key);
        if (!feature) return;
        const baselineColor = hasOwn(baseline, 'originalStrokeColor')
          ? baseline.originalStrokeColor
          : originalColor;
        const baselineWidth = hasOwn(baseline, 'originalStrokeWidth')
          ? baseline.originalStrokeWidth
          : originalWidth;
        getFeatureElements(svg, feature.svg_id).forEach((element) => {
          if (restoreStrokeAttributes(element, baselineColor, baselineWidth)) changed = true;
        });
      });
    } else {
      svg.querySelectorAll(FEATURE_SELECTOR).forEach((path) => {
        if (restoreStrokeAttributes(path, originalColor, originalWidth)) changed = true;
      });
    }

    const reapplied = reapplyStrokeOverrides();
    if (changed && !reapplied) persistStrokeEdit('history-stroke-reconcile');
    return changed || reapplied;
  };

  return {
    applyStrokeToFeaturesByCaption,
    captureOriginalStroke,
    closeLegendStrokeOptions,
    getLegendEntryStrokeColor,
    getLegendEntryStrokeWidth,
    isLegendStrokeOptionsOpen,
    reconcileStrokeOverrides,
    reapplyStrokeOverrides,
    resetAllStrokes,
    resetLegendEntryStroke,
    setLegendEntryStrokeColorValue,
    toggleLegendStrokeOptions,
    updateLegendEntryStrokeColor,
    updateLegendEntryStrokeWidth
  };
};
