import { getFeatureCaption, ruleMatchesFeature } from '../feature-utils.js';
import {
  FEATURE_SELECTOR,
  getFeatureElements
} from '../feature-editor/svg-actions.js';
import { getAllFeatureLegendGroups } from './utils.js';
import {
  applyLegendColorOverridesToSvg,
  applyStrokeOverridesToSvg,
  getLegendSwatches
} from '../editor-svg-projection.js';
import {
  featureOverrideKey,
  migrateLegacyFeatureOverrides
} from '../../services/feature-override-identity.js';

export { applyLegendColorOverridesToSvg, applyStrokeOverridesToSvg };

const hasOwn = (value, key) => Object.prototype.hasOwnProperty.call(value || {}, key);

const restoreStrokeAttributes = (element, originalColor, originalWidth, mutation) => {
  let changed = false;
  const nextColor = originalColor === null ? null : String(originalColor);
  const nextWidth = originalWidth === null ? null : String(originalWidth);
  if (nextColor === null) {
    if (element.hasAttribute('stroke')) {
      mutation.removeAttribute(element, 'stroke');
      changed = true;
    }
  } else if (element.getAttribute('stroke') !== nextColor) {
    mutation.setAttribute(element, 'stroke', nextColor);
    changed = true;
  }
  if (nextWidth === null) {
    if (element.hasAttribute('stroke-width')) {
      mutation.removeAttribute(element, 'stroke-width');
      changed = true;
    }
  } else if (element.getAttribute('stroke-width') !== nextWidth) {
    mutation.setAttribute(element, 'stroke-width', nextWidth);
    changed = true;
  }
  return changed;
};

export const createLegendStrokeActions = ({ state, previewRuntime = null }) => {
  if (
    typeof previewRuntime?.commitDomEdit !== 'function'
    || typeof previewRuntime?.runDomEditSync !== 'function'
  ) {
    throw new Error('createLegendStrokeActions requires the preview runtime edit protocol.');
  }
  const {
    extractedFeatures,
    legendEntries,
    legendStrokeOverrides,
    featureStrokeOverrides,
    originalSvgStroke,
    svgContainer,
    manualSpecificRules
  } = state;

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

  const commitStrokeMutation = (reason, mutate, { lease = null } = {}) => (
    lease
      ? lease.mutate(mutate)
      : previewRuntime.commitDomEdit({
          reason,
          invalidateIndexes: ['features', 'legend'],
          mutate
        }).changed
  );

  const canonicalStrokeState = () => [
    { target: legendStrokeOverrides },
    { target: featureStrokeOverrides },
    { target: legendEntries, key: 'value', deep: true },
    { target: originalSvgStroke, key: 'value', deep: true },
    { target: manualSpecificRules }
  ];
  const strokeAction = (reason, action) => (...args) => {
    const options = args.at(-1);
    if (options?.lease) return action(...args);
    return previewRuntime.runDomEditSync({
      reason,
      canonicalState: canonicalStrokeState(),
      action: () => action(...args)
    });
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
    if (!svgContainer.value?.querySelector?.('svg')) return false;

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
    commitStrokeMutation('reset-legend-stroke', ({ svg, mutation }) => {
      if (entry.featureIds && entry.featureIds.length > 0) {
        for (const svgId of entry.featureIds) {
          getFeatureElements(svg, svgId).forEach((element) => {
            if (restoreStrokeAttributes(element, originalColor, originalWidth, mutation)) {
              updatedCount++;
            }
          });
        }
      } else {
        const legendFillColor = entry.color;
        if (legendFillColor) {
          const normalizedColor = legendFillColor.toLowerCase();
          svg.querySelectorAll(FEATURE_SELECTOR).forEach((path) => {
            const fill = path.getAttribute('fill');
            if (
              fill
              && fill.toLowerCase() === normalizedColor
              && restoreStrokeAttributes(path, originalColor, originalWidth, mutation)
            ) {
              updatedCount++;
            }
          });
        }
      }

      getLegendSwatches(svg, entry.caption).forEach((swatch) => {
        if (restoreStrokeAttributes(swatch, originalSwatchColor, originalSwatchWidth, mutation)) {
          updatedCount++;
        }
      });
      return updatedCount;
    });

    const overrideRemoved = Object.prototype.hasOwnProperty.call(
      legendStrokeOverrides,
      entry.caption
    );
    if (overrideRemoved) delete legendStrokeOverrides[entry.caption];

    return updatedCount > 0 || overrideRemoved;
  };

  const resetAllStrokes = () => {
    if (!svgContainer.value?.querySelector?.('svg')) return false;
    const originalColor = originalSvgStroke.value.color;
    const originalWidth = originalSvgStroke.value.width;
    let updatedCount = 0;
    commitStrokeMutation('reset-all-strokes', ({ svg, mutation }) => {
      svg.querySelectorAll(FEATURE_SELECTOR).forEach((path) => {
        if (restoreStrokeAttributes(path, originalColor, originalWidth, mutation)) updatedCount++;
      });

      for (const targetGroup of getAllFeatureLegendGroups(svg)) {
        targetGroup.querySelectorAll('path').forEach((path) => {
          const fill = path.getAttribute('fill');
          if (
            fill
            && fill !== 'none'
            && !fill.startsWith('url(')
            && restoreStrokeAttributes(path, originalColor, originalWidth, mutation)
          ) {
            updatedCount++;
          }
        });
      }
      return updatedCount;
    });

    const overridesRemoved =
      Object.keys(legendStrokeOverrides).length > 0 ||
      Object.keys(featureStrokeOverrides).length > 0;
    Object.keys(legendStrokeOverrides).forEach((key) => delete legendStrokeOverrides[key]);
    Object.keys(featureStrokeOverrides).forEach((key) => delete featureStrokeOverrides[key]);

    if (updatedCount > 0) {
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
    let updatedCount = 0;
    const legendEntry = legendEntries.value.find((e) => e.caption === caption);
    commitStrokeMutation('legend-stroke', ({ svg, mutation }) => {
      const processedIds = new Set();
      const applyStroke = (element) => {
        if (removeStroke) {
          if (element.getAttribute('stroke') !== null) {
            mutation.removeAttribute(element, 'stroke');
            updatedCount++;
          }
        } else if (strokeColor !== null && element.getAttribute('stroke') !== String(strokeColor)) {
          mutation.setAttribute(element, 'stroke', strokeColor);
          updatedCount++;
        }
        if (
          strokeWidth !== null
          && element.getAttribute('stroke-width') !== String(strokeWidth)
        ) {
          mutation.setAttribute(element, 'stroke-width', strokeWidth);
          updatedCount++;
        }
      };

      if (legendEntry?.featureIds?.length > 0) {
        for (const svgId of legendEntry.featureIds) {
          if (processedIds.has(svgId)) continue;
          processedIds.add(svgId);
          getFeatureElements(svg, svgId).forEach(applyStroke);
        }
      } else {
        const normalizedLegendColor = String(legendEntry?.color || '').toLowerCase();
        if (normalizedLegendColor) {
          svg.querySelectorAll(FEATURE_SELECTOR).forEach((path) => {
            if (String(path.getAttribute('fill') || '').toLowerCase() !== normalizedLegendColor) return;
            const pathId = path.getAttribute('id');
            if (pathId && !processedIds.has(pathId)) {
              processedIds.add(pathId);
              applyStroke(path);
            }
          });
        }
      }

      for (const targetGroup of getAllFeatureLegendGroups(svg)) {
        const entryGroup = targetGroup.querySelector(
          `g[data-legend-key="${CSS.escape(caption)}"]`
        );
        if (!entryGroup) continue;
        for (const path of entryGroup.querySelectorAll('path')) {
          const fill = path.getAttribute('fill');
          if (fill && fill !== 'none' && !fill.startsWith('url(')) {
            applyStroke(path);
            break;
          }
        }
      }
      return updatedCount;
    });

    if (updatedCount > 0) {
      console.log(`Applied stroke to ${updatedCount} elements for caption "${caption}"`);
    }
    return updatedCount > 0;
  };

  const reapplyStrokeOverrides = () => {
    if (
      Object.keys(legendStrokeOverrides).length === 0 &&
      Object.keys(featureStrokeOverrides).length === 0
    ) return false;
    migrateLegacyFeatureOverrides(featureStrokeOverrides, extractedFeatures.value);
    let totalUpdated = 0;
    commitStrokeMutation('reapply-stroke-overrides', ({ svg, mutation }) => {
      totalUpdated = applyStrokeOverridesToSvg({
        svg,
        features: extractedFeatures.value,
        legendStrokeOverrides,
        featureStrokeOverrides,
        mutation
      });
      return totalUpdated;
    });
    return totalUpdated > 0;
  };

  const reconcileStrokeOverrides = ({ changes = null, lease = null } = {}) => {
    const originalColor = originalSvgStroke.value.color;
    const originalWidth = originalSvgStroke.value.width;
    const historyChanges = Array.isArray(changes) ? changes : null;
    const featureBaselines = new Map();
    const legendBaselines = new Map();

    if (historyChanges) {
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
    }

    migrateLegacyFeatureOverrides(featureStrokeOverrides, extractedFeatures.value);
    let updatedCount = 0;
    commitStrokeMutation('history-stroke-reconcile', ({ svg, mutation }) => {
      mutation.captureState(legendStrokeOverrides);
      mutation.captureState(featureStrokeOverrides);
      if (historyChanges) {
        legendBaselines.forEach((baseline, caption) => {
          const baselineColor = hasOwn(baseline, 'originalStrokeColor')
            ? baseline.originalStrokeColor
            : originalColor;
          const baselineWidth = hasOwn(baseline, 'originalStrokeWidth')
            ? baseline.originalStrokeWidth
            : originalWidth;
          const featureIds = new Set();
          const entry = legendEntries.value.find((candidate) => candidate.caption === caption);
          (entry?.featureIds || []).forEach((featureId) => {
            featureIds.add(String(featureId || '').trim());
          });
          extractedFeatures.value
            .filter((feature) => getFeatureCaption(feature) === caption)
            .forEach((feature) => featureIds.add(String(feature?.svg_id || '').trim()));
          featureIds.forEach((featureId) => {
            getFeatureElements(svg, featureId).forEach((element) => {
              if (restoreStrokeAttributes(element, originalColor, originalWidth, mutation)) {
                updatedCount++;
              }
            });
          });
          getLegendSwatches(svg, caption).forEach((swatch) => {
            if (restoreStrokeAttributes(swatch, baselineColor, baselineWidth, mutation)) {
              updatedCount++;
            }
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
            if (restoreStrokeAttributes(element, baselineColor, baselineWidth, mutation)) {
              updatedCount++;
            }
          });
        });
      } else {
        svg.querySelectorAll(FEATURE_SELECTOR).forEach((path) => {
          if (restoreStrokeAttributes(path, originalColor, originalWidth, mutation)) {
            updatedCount++;
          }
        });
      }

      updatedCount += applyStrokeOverridesToSvg({
        svg,
        features: extractedFeatures.value,
        legendStrokeOverrides,
        featureStrokeOverrides,
        mutation
      });
      return updatedCount;
    }, { lease });
    return updatedCount > 0;
  };

  return {
    applyStrokeToFeaturesByCaption: strokeAction('legend-stroke', applyStrokeToFeaturesByCaption),
    captureOriginalStrokeValues,
    getLegendEntryStrokeColor,
    getLegendEntryStrokeWidth,
    reconcileStrokeOverrides: strokeAction('legend-stroke-reconcile', reconcileStrokeOverrides),
    reapplyStrokeOverrides: strokeAction('legend-stroke-reapply', reapplyStrokeOverrides),
    resetAllStrokes: strokeAction('legend-stroke-reset-all', resetAllStrokes),
    resetLegendEntryStroke: strokeAction('legend-stroke-reset', resetLegendEntryStroke),
    setLegendEntryStrokeColorValue: strokeAction(
      'legend-stroke-color',
      setLegendEntryStrokeColorValue
    ),
    updateLegendEntryStrokeColor: strokeAction(
      'legend-stroke-color',
      updateLegendEntryStrokeColor
    ),
    updateLegendEntryStrokeWidth: strokeAction(
      'legend-stroke-width',
      updateLegendEntryStrokeWidth
    )
  };
};
