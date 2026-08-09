import { getFeatureCaption, ruleMatchesFeature } from '../feature-utils.js';
import { FEATURE_SELECTOR, getFeatureElements } from '../feature-editor/svg-actions.js';
import { getAllFeatureLegendGroups } from './utils.js';
import { serializeCleanSvg } from '../../services/svg-serialization.js';
import {
  featureOverrideKey,
  migrateLegacyFeatureOverrides
} from '../../services/feature-override-identity.js';

export const createLegendStrokeActions = ({ state }) => {
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
    if (!entry) return;

    if (!legendStrokeOverrides[entry.caption]) {
      legendStrokeOverrides[entry.caption] = {
        strokeColor: color,
        strokeWidth: getLegendEntryStrokeWidth(idx)
      };
    }
    legendStrokeOverrides[entry.caption].strokeColor = color;

    applyStrokeToFeaturesByCaption(entry.caption, color, null);
  };

  const updateLegendEntryStrokeWidth = (idx, width) => {
    const entry = legendEntries.value[idx];
    if (!entry) return;

    const widthVal = parseFloat(width);
    if (isNaN(widthVal)) return;

    if (!legendStrokeOverrides[entry.caption]) {
      legendStrokeOverrides[entry.caption] = {
        strokeColor: getLegendEntryStrokeColor(idx),
        strokeWidth: widthVal
      };
    }
    legendStrokeOverrides[entry.caption].strokeWidth = widthVal;

    applyStrokeToFeaturesByCaption(entry.caption, null, widthVal);
  };

  const setLegendEntryStrokeColorValue = (idx, value) => {
    const entry = legendEntries.value[idx];
    if (!entry) return;
    if (value !== null) {
      updateLegendEntryStrokeColor(idx, String(value || '').trim());
      return;
    }
    const override = legendStrokeOverrides[entry.caption];
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
  };

  const resetLegendEntryStroke = (idx) => {
    const entry = legendEntries.value[idx];
    if (!entry) return;
    if (!svgContainer.value) return;

    const svg = svgContainer.value.querySelector('svg');
    if (!svg) return;

    const originalColor = originalSvgStroke.value.color;
    const originalWidth = originalSvgStroke.value.width;
    let updatedCount = 0;

    if (entry.featureIds && entry.featureIds.length > 0) {
      for (const svgId of entry.featureIds) {
        const elements = getFeatureElements(svg, svgId);
        elements.forEach((el) => {
          if (originalColor === null) {
            el.removeAttribute('stroke');
          } else {
            el.setAttribute('stroke', originalColor);
          }
          if (originalWidth === null) {
            el.removeAttribute('stroke-width');
          } else {
            el.setAttribute('stroke-width', originalWidth);
          }
          updatedCount++;
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
            if (originalColor === null) {
              path.removeAttribute('stroke');
            } else {
              path.setAttribute('stroke', originalColor);
            }
            if (originalWidth === null) {
              path.removeAttribute('stroke-width');
            } else {
              path.setAttribute('stroke-width', originalWidth);
            }
            updatedCount++;
          }
        });
      }
    }

    delete legendStrokeOverrides[entry.caption];

    if (updatedCount > 0) {
      skipCaptureBaseConfig.value = true;
      const resultIdx = selectedResultIndex.value;
      if (resultIdx >= 0 && results.value.length > resultIdx) {
        results.value[resultIdx] = { ...results.value[resultIdx], content: serializeCleanSvg(svg) };
      }
    }
  };

  const resetAllStrokes = () => {
    if (!svgContainer.value) return;
    const svg = svgContainer.value.querySelector('svg');
    if (!svg) return;

    const originalColor = originalSvgStroke.value.color;
    const originalWidth = originalSvgStroke.value.width;

    const featurePaths = svg.querySelectorAll(FEATURE_SELECTOR);
    let updatedCount = 0;
    featurePaths.forEach((path) => {
      if (originalColor === null) {
        path.removeAttribute('stroke');
      } else {
        path.setAttribute('stroke', originalColor);
      }
      if (originalWidth === null) {
        path.removeAttribute('stroke-width');
      } else {
        path.setAttribute('stroke-width', originalWidth);
      }
      updatedCount++;
    });

    const legendGroups = getAllFeatureLegendGroups(svg);
    for (const targetGroup of legendGroups) {
      const paths = targetGroup.querySelectorAll('path');
      paths.forEach((p) => {
        const fill = p.getAttribute('fill');
        if (fill && fill !== 'none' && !fill.startsWith('url(')) {
          if (originalColor === null) {
            p.removeAttribute('stroke');
          } else {
            p.setAttribute('stroke', originalColor);
          }
          if (originalWidth === null) {
            p.removeAttribute('stroke-width');
          } else {
            p.setAttribute('stroke-width', originalWidth);
          }
        }
      });
    }

    Object.keys(legendStrokeOverrides).forEach((key) => delete legendStrokeOverrides[key]);
    Object.keys(featureStrokeOverrides).forEach((key) => delete featureStrokeOverrides[key]);

    if (updatedCount > 0) {
      skipCaptureBaseConfig.value = true;
      const resultIdx = selectedResultIndex.value;
      if (resultIdx >= 0 && results.value.length > resultIdx) {
        results.value[resultIdx] = { ...results.value[resultIdx], content: serializeCleanSvg(svg) };
      }
      console.log(
        `Reset all strokes: updated ${updatedCount} elements to original (color=${originalColor}, width=${originalWidth})`
      );
    }
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
            el.removeAttribute('stroke');
            updatedCount++;
          } else if (strokeColor !== null) {
            el.setAttribute('stroke', strokeColor);
            updatedCount++;
          }
          if (strokeWidth !== null) {
            el.setAttribute('stroke-width', strokeWidth);
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
                path.removeAttribute('stroke');
                updatedCount++;
              } else if (strokeColor !== null) {
                path.setAttribute('stroke', strokeColor);
                updatedCount++;
              }
              if (strokeWidth !== null) {
                path.setAttribute('stroke-width', strokeWidth);
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
              path.removeAttribute('stroke');
            } else if (strokeColor !== null) {
              path.setAttribute('stroke', strokeColor);
            }
            if (strokeWidth !== null) {
              path.setAttribute('stroke-width', strokeWidth);
            }
            break;
          }
        }
      }
    }

    if (updatedCount > 0) {
      skipCaptureBaseConfig.value = true;
      const resultIdx = selectedResultIndex.value;
      if (resultIdx >= 0 && results.value.length > resultIdx) {
        results.value[resultIdx] = { ...results.value[resultIdx], content: serializeCleanSvg(svg) };
      }
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

    let totalUpdated = 0;

    for (const [caption, overrides] of Object.entries(legendStrokeOverrides)) {
      const { strokeColor, strokeWidth } = overrides;
      if (!strokeColor && strokeWidth === undefined) continue;

      const matchingFeatures = extractedFeatures.value.filter((f) => getFeatureCaption(f) === caption);

      for (const feat of matchingFeatures) {
        const elements = getFeatureElements(svg, feat.svg_id);
        elements.forEach((el) => {
          if (strokeColor) {
            el.setAttribute('stroke', strokeColor);
          }
          if (strokeWidth !== undefined && strokeWidth !== null) {
            el.setAttribute('stroke-width', strokeWidth);
          }
          totalUpdated++;
        });
      }
    }

    for (const feature of extractedFeatures.value) {
      const featureKey = featureOverrideKey(feature);
      const overrides = featureKey ? featureStrokeOverrides[featureKey] : null;
      const { strokeColor, strokeWidth } = overrides || {};
      if (!strokeColor && strokeWidth === undefined) continue;
      const elements = getFeatureElements(svg, feature.svg_id);
      elements.forEach((el) => {
        if (strokeColor) {
          el.setAttribute('stroke', strokeColor);
        }
        if (strokeWidth !== undefined && strokeWidth !== null) {
          el.setAttribute('stroke-width', strokeWidth);
        }
        totalUpdated++;
      });
    }

    if (totalUpdated > 0) {
      skipCaptureBaseConfig.value = true;
      const resultIdx = selectedResultIndex.value;
      if (resultIdx >= 0 && results.value.length > resultIdx) {
        results.value[resultIdx] = { ...results.value[resultIdx], content: serializeCleanSvg(svg) };
      }
    }
  };

  return {
    applyStrokeToFeaturesByCaption,
    captureOriginalStrokeValues,
    getLegendEntryStrokeColor,
    getLegendEntryStrokeWidth,
    reapplyStrokeOverrides,
    resetAllStrokes,
    resetLegendEntryStroke,
    setLegendEntryStrokeColorValue,
    updateLegendEntryStrokeColor,
    updateLegendEntryStrokeWidth
  };
};
