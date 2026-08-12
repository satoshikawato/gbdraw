import { getFeatureCaption, ruleMatchesFeature } from '../feature-utils.js';
import {
  FEATURE_SELECTOR,
  getFeatureElementIndex,
  getFeatureElements
} from '../feature-editor/svg-actions.js';
import { getAllFeatureLegendGroups } from './utils.js';
import { serializeCleanSvg } from '../../services/svg-serialization.js';
import {
  featureOverrideKey,
  migrateLegacyFeatureOverrides
} from '../../services/feature-override-identity.js';

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
    getFeatureElements(svg, feature?.svg_id, featureIndex).forEach((element) => {
      if (applyStrokeAttributes(element, overrides.strokeColor, overrides.strokeWidth)) {
        changedCount += 1;
      }
    });
  };

  Object.entries(legendStrokeOverrides || {}).forEach(([caption, overrides]) => {
    if (!overrides || typeof overrides !== 'object') return;
    (Array.isArray(features) ? features : [])
      .filter((feature) => getFeatureCaption(feature) === caption)
      .forEach((feature) => applyToFeature(feature, overrides));
    getAllFeatureLegendGroups(svg).forEach((targetGroup) => {
      const entryGroup = Array.from(targetGroup.querySelectorAll('g[data-legend-key]'))
        .find((entry) => entry.getAttribute('data-legend-key') === caption);
      const swatch = Array.from(entryGroup?.querySelectorAll?.('path') || []).find((path) => {
        const fill = path.getAttribute('fill');
        return fill && fill !== 'none' && !fill.startsWith('url(');
      });
      if (swatch && applyStrokeAttributes(swatch, overrides.strokeColor, overrides.strokeWidth)) {
        changedCount += 1;
      }
    });
  });

  (Array.isArray(features) ? features : []).forEach((feature) => {
    const key = featureOverrideKey(feature);
    if (key) applyToFeature(feature, featureStrokeOverrides?.[key]);
  });
  return changedCount;
};

export const createLegendStrokeActions = ({ state, previewRuntime = null }) => {
  const {
    extractedFeatures,
    legendEntries,
    legendStrokeOverrides,
    featureStrokeOverrides,
    originalSvgStroke,
    results,
    selectedResultIndex,
    skipCaptureBaseConfig,
    svgContainer,
    manualSpecificRules
  } = state;

  const persistStrokeEdit = (svg, reason = 'feature-stroke') => {
    skipCaptureBaseConfig.value = true;
    if (previewRuntime?.markActiveResultDirty?.(reason)) {
      previewRuntime.flushActiveResult?.();
      return;
    }
    const resultIdx = selectedResultIndex.value;
    if (resultIdx >= 0 && results.value.length > resultIdx) {
      results.value[resultIdx] = { ...results.value[resultIdx], content: serializeCleanSvg(svg) };
    }
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
    const entry = legendEntries.value[idx];
    if (!entry) return false;
    const normalized = String(color || '').trim();
    if (String(legendStrokeOverrides[entry.caption]?.strokeColor || '').trim() === normalized) {
      return false;
    }

    if (!legendStrokeOverrides[entry.caption]) {
      legendStrokeOverrides[entry.caption] = {
        strokeColor: normalized,
        strokeWidth: getLegendEntryStrokeWidth(idx)
      };
    }
    legendStrokeOverrides[entry.caption].strokeColor = normalized;

    applyStrokeToFeaturesByCaption(entry.caption, normalized, null);
    return true;
  };

  const updateLegendEntryStrokeWidth = (idx, width) => {
    const entry = legendEntries.value[idx];
    if (!entry) return false;

    const widthVal = parseFloat(width);
    if (isNaN(widthVal)) return false;
    if (Number(legendStrokeOverrides[entry.caption]?.strokeWidth) === widthVal) return false;

    if (!legendStrokeOverrides[entry.caption]) {
      legendStrokeOverrides[entry.caption] = {
        strokeColor: getLegendEntryStrokeColor(idx),
        strokeWidth: widthVal
      };
    }
    legendStrokeOverrides[entry.caption].strokeWidth = widthVal;

    applyStrokeToFeaturesByCaption(entry.caption, null, widthVal);
    return true;
  };

  const setLegendEntryStrokeColorValue = (idx, value) => {
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
    const entry = legendEntries.value[idx];
    if (!entry) return false;
    if (!svgContainer.value) return false;

    const svg = svgContainer.value.querySelector('svg');
    if (!svg) return false;

    const originalColor = originalSvgStroke.value.color;
    const originalWidth = originalSvgStroke.value.width;
    let updatedCount = 0;

    if (entry.featureIds && entry.featureIds.length > 0) {
      for (const svgId of entry.featureIds) {
        const elements = getFeatureElements(svg, svgId);
        elements.forEach((el) => {
          if (restoreStrokeAttributes(el, originalColor, originalWidth)) updatedCount++;
        });
      }
    } else {
      const legendFillColor = entry.color;
      if (legendFillColor) {
        const normalizedColor = legendFillColor.toLowerCase();
        const featurePaths = svg.querySelectorAll(FEATURE_SELECTOR);
        featurePaths.forEach((path) => {
          const fill = path.getAttribute('fill');
          if (fill && fill.toLowerCase() === normalizedColor) {
            if (restoreStrokeAttributes(path, originalColor, originalWidth)) updatedCount++;
          }
        });
      }
    }

    getAllFeatureLegendGroups(svg).forEach((targetGroup) => {
      const entryGroup = Array.from(targetGroup.querySelectorAll('g[data-legend-key]'))
        .find((candidate) => candidate.getAttribute('data-legend-key') === entry.caption);
      const swatch = Array.from(entryGroup?.querySelectorAll?.('path') || []).find((path) => {
        const fill = path.getAttribute('fill');
        return fill && fill !== 'none' && !fill.startsWith('url(');
      });
      if (swatch && restoreStrokeAttributes(swatch, originalColor, originalWidth)) {
        updatedCount++;
      }
    });

    const overrideRemoved = Object.prototype.hasOwnProperty.call(
      legendStrokeOverrides,
      entry.caption
    );
    if (overrideRemoved) delete legendStrokeOverrides[entry.caption];

    if (updatedCount > 0) {
      persistStrokeEdit(svg, 'reset-legend-stroke');
    }
    return updatedCount > 0 || overrideRemoved;
  };

  const resetAllStrokes = () => {
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
      persistStrokeEdit(svg, 'reset-all-strokes');
      console.log(
        `Reset all strokes: updated ${updatedCount} elements to original (color=${originalColor}, width=${originalWidth})`
      );
    }
    return updatedCount > 0 || overridesRemoved;
  };

  const captureOriginalStrokeValues = (caption) => {
    if (!svgContainer.value) return { strokeColor: '#000000', strokeWidth: 0.5 };
    const svg = svgContainer.value.querySelector('svg');
    if (!svg) return { strokeColor: '#000000', strokeWidth: 0.5 };

    const legendEntry = legendEntries.value.find((e) => e.caption === caption);
    const legendFillColor = legendEntry?.color;

    const matchingFeatures = extractedFeatures.value.filter((f) => getFeatureCaption(f) === caption);

    const ruleMatch = manualSpecificRules.find((r) => r.cap === caption);
    if (ruleMatch) {
      const ruleFeatures = extractedFeatures.value.filter((f) => ruleMatchesFeature(f, ruleMatch));
      matchingFeatures.push(...ruleFeatures);
    }

    for (const feat of matchingFeatures) {
      const el = getFeatureElements(svg, feat.svg_id)[0] || null;
      if (el) {
        return {
          strokeColor: el.getAttribute('stroke') || '#000000',
          strokeWidth: parseFloat(el.getAttribute('stroke-width')) || 0.5
        };
      }
    }

    if (legendFillColor) {
      const normalizedLegendColor = legendFillColor.toLowerCase();
      const featurePaths = svg.querySelectorAll(FEATURE_SELECTOR);
      for (const path of featurePaths) {
        const fill = path.getAttribute('fill');
        if (fill && fill.toLowerCase() === normalizedLegendColor) {
          return {
            strokeColor: path.getAttribute('stroke') || '#000000',
            strokeWidth: parseFloat(path.getAttribute('stroke-width')) || 0.5
          };
        }
      }
    }

    return { strokeColor: '#000000', strokeWidth: 0.5 };
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
    const processedIds = new Set();

    const legendEntry = legendEntries.value.find((e) => e.caption === caption);

    if (legendEntry && legendEntry.featureIds && legendEntry.featureIds.length > 0) {
      for (const svgId of legendEntry.featureIds) {
        if (processedIds.has(svgId)) continue;
        processedIds.add(svgId);
        const elements = getFeatureElements(svg, svgId);
        elements.forEach((el) => {
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
      }
      console.log(`Applied stroke to ${legendEntry.featureIds.length} features via featureIds for "${caption}"`);
    } else {
      const legendFillColor = legendEntry?.color;
      if (legendFillColor) {
        const normalizedLegendColor = legendFillColor.toLowerCase();
        const featurePaths = svg.querySelectorAll(FEATURE_SELECTOR);
        featurePaths.forEach((path) => {
          const fill = path.getAttribute('fill');
          if (fill && fill.toLowerCase() === normalizedLegendColor) {
            const pathId = path.getAttribute('id');
            if (pathId && !processedIds.has(pathId)) {
              processedIds.add(pathId);
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
            }
          }
        });
      }
      console.log(`Fallback: Applied stroke via fill color matching for "${caption}"`);
    }

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
      persistStrokeEdit(svg, 'legend-stroke');
      console.log(`Applied stroke to ${updatedCount} elements for caption "${caption}"`);
    }
  };

  const reapplyStrokeOverrides = () => {
    if (!svgContainer.value) return;
    const svg = svgContainer.value.querySelector('svg');
    if (!svg) return;

    if (
      Object.keys(legendStrokeOverrides).length === 0 &&
      Object.keys(featureStrokeOverrides).length === 0
    ) return;
    migrateLegacyFeatureOverrides(featureStrokeOverrides, extractedFeatures.value);
    const totalUpdated = applyStrokeOverridesToSvg({
      svg,
      features: extractedFeatures.value,
      legendStrokeOverrides,
      featureStrokeOverrides
    });

    if (totalUpdated > 0) {
      persistStrokeEdit(svg, 'reapply-stroke-overrides');
    }
  };

  const reconcileStrokeOverrides = () => {
    const svg = svgContainer.value?.querySelector?.('svg');
    if (!svg) return false;
    const originalColor = originalSvgStroke.value.color;
    const originalWidth = originalSvgStroke.value.width;
    let changed = false;
    svg.querySelectorAll(FEATURE_SELECTOR).forEach((path) => {
      const nextColor = originalColor === null ? null : String(originalColor);
      const nextWidth = originalWidth === null ? null : String(originalWidth);
      if (nextColor === null) {
        if (path.hasAttribute('stroke')) changed = true;
        path.removeAttribute('stroke');
      } else if (path.getAttribute('stroke') !== nextColor) {
        path.setAttribute('stroke', nextColor);
        changed = true;
      }
      if (nextWidth === null) {
        if (path.hasAttribute('stroke-width')) changed = true;
        path.removeAttribute('stroke-width');
      } else if (path.getAttribute('stroke-width') !== nextWidth) {
        path.setAttribute('stroke-width', nextWidth);
        changed = true;
      }
    });
    reapplyStrokeOverrides();
    if (changed) persistStrokeEdit(svg, 'history-stroke-reconcile');
    return changed;
  };

  return {
    applyStrokeToFeaturesByCaption,
    captureOriginalStrokeValues,
    getLegendEntryStrokeColor,
    getLegendEntryStrokeWidth,
    reconcileStrokeOverrides,
    reapplyStrokeOverrides,
    resetAllStrokes,
    resetLegendEntryStroke,
    setLegendEntryStrokeColorValue,
    updateLegendEntryStrokeColor,
    updateLegendEntryStrokeWidth
  };
};
