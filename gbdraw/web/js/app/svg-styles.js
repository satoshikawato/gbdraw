import {
  interpolateColor,
  resolveCollinearMatchColor,
  resolvePairwiseLegendGradientColorKeys
} from './color-utils.js';
import {
  FEATURE_SELECTOR,
  getFeatureElementIndex,
  getFeatureFillElements
} from './feature-editor/svg-actions.js';
import {
  FEATURE_PART_CONNECTOR,
  getFeaturePart,
  isFeatureFillTarget
} from './feature-dom.js';
import { ruleMatchesFeature } from './feature-utils.js';
import {
  getAllFeatureLegendGroups,
  PAIRWISE_LEGEND_SELECTOR,
  parseTransformXY
} from './legend/utils.js';
import { serializeCleanSvg } from '../services/svg-serialization.js';
import { getGroupsByBaseIds } from '../services/svg-result-normalization.js';

export { getGroupsByBaseIds } from '../services/svg-result-normalization.js';

const normalizeComparableColor = (value) => String(value || '').trim().toLowerCase();

const setColorAttributeIfChanged = (element, attribute, value) => {
  if (normalizeComparableColor(element.getAttribute(attribute)) === normalizeComparableColor(value)) {
    return false;
  }
  element.setAttribute(attribute, value);
  return true;
};

export const NON_FEATURE_PALETTE_COLOR_KEYS = Object.freeze([
  'gc_content',
  'skew_high',
  'skew_low',
  'pairwise_match',
  'pairwise_match_min',
  'pairwise_match_max',
  'collinear_block_plus_min',
  'collinear_block_plus',
  'collinear_block_minus_min',
  'collinear_block_minus'
]);

const BUILT_IN_LEGEND_PALETTE_KEYS = Object.freeze({
  'gc content': 'gc_content',
  'gc skew (+)': 'skew_high',
  'gc skew (-)': 'skew_low'
});

/** Return the applied-palette owner for an editable built-in legend swatch. */
export const builtInLegendPaletteColorKey = (caption) => (
  BUILT_IN_LEGEND_PALETTE_KEYS[String(caption ?? '').trim().toLowerCase()] || ''
);

const NON_FEATURE_PALETTE_COLOR_KEY_SET = new Set(NON_FEATURE_PALETTE_COLOR_KEYS);
const text = (value) => String(value ?? '').trim();

export const changedNonFeaturePaletteColorKeys = (beforeColors = {}, afterColors = {}) => (
  NON_FEATURE_PALETTE_COLOR_KEYS.filter(
    (key) => normalizeComparableColor(beforeColors?.[key]) !== normalizeComparableColor(afterColors?.[key])
  )
);

const normalizedDinucleotide = (value, fallback = 'GC') => {
  const normalized = text(value || fallback).toUpperCase();
  return /^[ACGTU]{2}$/.test(normalized) ? normalized : text(fallback).toUpperCase() || 'GC';
};

/**
 * Freeze the track-slot ownership needed to distinguish palette-derived skew
 * colors from explicit per-track colors during a detached Result projection.
 */
export const buildNonFeaturePaletteProjectionContext = ({
  trackSlots = [],
  defaultDinucleotide = 'GC'
} = {}) => Object.freeze({
  skewSlots: Object.freeze(Object.fromEntries(
    (Array.isArray(trackSlots) ? trackSlots : [])
      .filter((slot) => text(slot?.renderer) === 'dinucleotide_skew' && text(slot?.id))
      .map((slot) => {
        const params = slot?.params && typeof slot.params === 'object' ? slot.params : {};
        const dinucleotide = normalizedDinucleotide(
          params.nt ?? params.dinucleotide,
          defaultDinucleotide
        );
        return [text(slot.id), Object.freeze({
          highPaletteOwned: !text(params.positive_color ?? params.high_color),
          lowPaletteOwned: !text(params.negative_color ?? params.low_color),
          legendLabel: text(params.legend_label) || `${dinucleotide} skew`
        })];
      })
  )),
  contentSlots: Object.freeze(Object.fromEntries(
    (Array.isArray(trackSlots) ? trackSlots : [])
      .filter((slot) => text(slot?.renderer) === 'dinucleotide_content' && text(slot?.id))
      .map((slot) => {
        const params = slot?.params && typeof slot.params === 'object' ? slot.params : {};
        const dinucleotide = normalizedDinucleotide(
          params.nt ?? params.dinucleotide,
          defaultDinucleotide
        );
        return [text(slot.id), Object.freeze({
          legendLabel: text(params.legend_label) || `${dinucleotide} content`
        })];
      })
  ))
});

const setProjectedColor = (element, attribute, value, counters, counterKey) => {
  if (!element || !text(value)) return false;
  if (!setColorAttributeIfChanged(element, attribute, value)) return false;
  counters[counterKey] += 1;
  return true;
};

const visibleFilledPaths = (group) => Array.from(group?.querySelectorAll?.('path') || []).filter((path) => {
  const fill = normalizeComparableColor(path.getAttribute?.('fill'));
  return fill && fill !== 'white' && fill !== '#ffffff' && fill !== 'none' && !fill.startsWith('url(');
});

const groupSlotId = (group) => text(
  group?.getAttribute?.('data-gbdraw-slot-id')
  || group?.getAttribute?.('id')
);

const defaultSkewSlot = (slotId) => (
  /^(?:gc_)?skew(?:_|$)/i.test(slotId)
  || /_gc_skew(?:_|$)/i.test(slotId)
);

const legendSwatchPath = (entryGroup) => Array.from(
  entryGroup?.querySelectorAll?.('path') || []
).find((path) => {
  const fill = text(path.getAttribute?.('fill'));
  return fill && fill.toLowerCase() !== 'none' && !fill.startsWith('url(');
}) || null;

const updateLegacyLegendSwatches = (legendGroup, colorsByCaption, counters) => {
  const texts = legendGroup?.querySelectorAll?.('text') || [];
  const paths = legendGroup?.querySelectorAll?.('path') || [];
  texts.forEach((textElement) => {
    const caption = text(textElement.textContent);
    const color = colorsByCaption.get(caption);
    if (!color) return;
    const textPosition = parseTransformXY(textElement.getAttribute?.('transform'));
    let bestPath = null;
    let bestX = -Infinity;
    paths.forEach((path) => {
      const pathPosition = parseTransformXY(path.getAttribute?.('transform'));
      const fill = text(path.getAttribute?.('fill'));
      if (
        Math.abs(pathPosition.y - textPosition.y) < 2
        && pathPosition.x < textPosition.x
        && fill
        && fill.toLowerCase() !== 'none'
        && !fill.startsWith('url(')
        && pathPosition.x > bestX
      ) {
        bestX = pathPosition.x;
        bestPath = path;
      }
    });
    setProjectedColor(bestPath, 'fill', color, counters, 'legendSwatches');
  });
};

const updateNonFeatureLegendSwatches = (svg, colorsByCaption, counters) => {
  getAllFeatureLegendGroups(svg).forEach((legendGroup) => {
    const keyedEntries = legendGroup.querySelectorAll?.('g[data-legend-key]') || [];
    if (keyedEntries.length === 0) {
      updateLegacyLegendSwatches(legendGroup, colorsByCaption, counters);
      return;
    }
    keyedEntries.forEach((entryGroup) => {
      const color = colorsByCaption.get(text(entryGroup.getAttribute?.('data-legend-key')));
      if (!color) return;
      setProjectedColor(
        legendSwatchPath(entryGroup),
        'fill',
        color,
        counters,
        'legendSwatches'
      );
    });
  });
};

const updatePairwiseLegendGradientStops = (pairwiseLegend, colors, counters = null) => {
  let updated = false;
  pairwiseLegend.querySelectorAll('linearGradient').forEach((gradient) => {
    const legendKey = gradient.closest('g[data-legend-key]')?.getAttribute('data-legend-key') || '';
    const { minKey, maxKey } = resolvePairwiseLegendGradientColorKeys(legendKey);
    const minColor = colors[minKey];
    const maxColor = colors[maxKey];
    if (!minColor || !maxColor) return;
    const stops = gradient.querySelectorAll('stop');
    if (stops.length < 2) return;
    const minChanged = setColorAttributeIfChanged(stops[0], 'stop-color', minColor);
    const maxChanged = setColorAttributeIfChanged(stops[1], 'stop-color', maxColor);
    if (minChanged || maxChanged) {
      updated = true;
      if (counters) counters.gradientStops += Number(minChanged) + Number(maxChanged);
    }
  });
  return updated;
};

const projectionError = (unsupported) => {
  const error = new Error(
    `Non-Feature palette projection requires regeneration metadata (${unsupported.length} unresolved target(s)).`
  );
  error.code = 'NON_FEATURE_PALETTE_PROJECTION_UNSUPPORTED';
  error.details = Object.freeze(unsupported.map((entry) => Object.freeze({ ...entry })));
  return error;
};

/**
 * Prepare GC/skew/comparison palette consumers on a detached clone. The input
 * SVG is never changed. Feature fills and Feature-derived legends intentionally
 * remain owned by the Feature bulk-style projection.
 */
export const prepareNonFeaturePaletteProjection = ({
  svg,
  beforeColors = {},
  afterColors = {},
  context = null,
  strict = true
} = {}) => {
  if (!svg?.cloneNode) throw new Error('Non-Feature palette projection requires an SVG root.');
  const projectedSvg = svg.cloneNode(true);
  const changedKeys = changedNonFeaturePaletteColorKeys(beforeColors, afterColors);
  const changedKeySet = new Set(changedKeys);
  const counters = {
    gcPaths: 0,
    skewPaths: 0,
    comparisonPaths: 0,
    legendSwatches: 0,
    gradientStops: 0
  };
  const unsupported = [];
  const colorsByCaption = new Map();

  if (changedKeySet.has('gc_content') && text(afterColors.gc_content)) {
    getGroupsByBaseIds(
      projectedSvg,
      ['gc_content'],
      ['dinucleotide_content']
    ).forEach((group) => {
      visibleFilledPaths(group).forEach((path) => {
        setProjectedColor(path, 'fill', afterColors.gc_content, counters, 'gcPaths');
      });
      const slotId = groupSlotId(group);
      const label = text(context?.contentSlots?.[slotId]?.legendLabel);
      if (label) colorsByCaption.set(label, afterColors.gc_content);
    });
    colorsByCaption.set('GC content', afterColors.gc_content);
  }

  if (changedKeySet.has('skew_high') || changedKeySet.has('skew_low')) {
    getGroupsByBaseIds(
      projectedSvg,
      ['skew', 'gc_skew'],
      ['dinucleotide_skew']
    ).forEach((group) => {
      const slotId = groupSlotId(group);
      const ownership = context?.skewSlots?.[slotId] || null;
      const isDefault = defaultSkewSlot(slotId);
      if (!ownership && !isDefault) {
        unsupported.push({ kind: 'skew-slot-ownership', slotId: slotId || '(missing)' });
        return;
      }
      const highPaletteOwned = ownership ? ownership.highPaletteOwned : true;
      const lowPaletteOwned = ownership ? ownership.lowPaletteOwned : true;
      const paths = visibleFilledPaths(group);
      if (changedKeySet.has('skew_high') && highPaletteOwned && paths[0]) {
        setProjectedColor(paths[0], 'fill', afterColors.skew_high, counters, 'skewPaths');
      }
      if (changedKeySet.has('skew_low') && lowPaletteOwned && paths[1]) {
        setProjectedColor(paths[1], 'fill', afterColors.skew_low, counters, 'skewPaths');
      }
      const label = text(ownership?.legendLabel) || (isDefault ? 'GC skew' : '');
      if (label && highPaletteOwned && text(afterColors.skew_high)) {
        colorsByCaption.set(`${label} (+)`, afterColors.skew_high);
      }
      if (label && lowPaletteOwned && text(afterColors.skew_low)) {
        colorsByCaption.set(`${label} (-)`, afterColors.skew_low);
      }
    });
  }

  const comparisonKeysChanged = changedKeys.some((key) => (
    key.startsWith('pairwise_match') || key.startsWith('collinear_block_')
  ));
  if (comparisonKeysChanged) {
    projectedSvg.querySelectorAll(
      'path[data-gbdraw-pairwise-match-id], path[data-gbdraw-match-id]'
    ).forEach((path) => {
      const colorMode = text(path.getAttribute('data-collinearity-color-mode'))
        .toLowerCase()
        .replace(/-/g, '_');
      const orientation = text(path.getAttribute('data-collinearity-orientation')).toLowerCase();
      const blockId = text(path.getAttribute('data-collinearity-block-id'));
      const factorText = text(path.getAttribute('data-identity-factor'));
      const factor = factorText ? Number(factorText) : NaN;
      let relevantKeys = ['pairwise_match_min', 'pairwise_match_max'];
      if (blockId && colorMode === 'orientation' && ['plus', 'minus'].includes(orientation)) {
        relevantKeys = [`collinear_block_${orientation}`];
      } else if (
        blockId
        && colorMode === 'orientation_identity'
        && ['plus', 'minus'].includes(orientation)
      ) {
        relevantKeys = [
          `collinear_block_${orientation}_min`,
          `collinear_block_${orientation}`
        ];
      }
      if (!relevantKeys.some((key) => changedKeySet.has(key))) return;
      if (!Number.isFinite(factor) && relevantKeys.length > 1) {
        unsupported.push({
          kind: 'comparison-identity-factor',
          matchId: text(
            path.getAttribute('data-gbdraw-pairwise-match-id')
            || path.getAttribute('data-gbdraw-match-id')
          ) || '(missing)'
        });
        return;
      }
      let color = resolveCollinearMatchColor({
        blockId,
        colorMode,
        orientation,
        identityFactor: Number.isFinite(factor) ? factor : null,
        colors: afterColors
      });
      if (!color) {
        const minColor = text(afterColors.pairwise_match_min);
        const maxColor = text(afterColors.pairwise_match_max);
        if (!minColor || !maxColor || !Number.isFinite(factor)) {
          unsupported.push({
            kind: 'comparison-palette-metadata',
            matchId: text(
              path.getAttribute('data-gbdraw-pairwise-match-id')
              || path.getAttribute('data-gbdraw-match-id')
            ) || '(missing)'
          });
          return;
        }
        color = interpolateColor(minColor, maxColor, factor);
      }
      setProjectedColor(path, 'fill', color, counters, 'comparisonPaths');
    });
    projectedSvg.querySelectorAll(PAIRWISE_LEGEND_SELECTOR).forEach((legend) => {
      updatePairwiseLegendGradientStops(legend, afterColors, counters);
    });
  }

  updateNonFeatureLegendSwatches(projectedSvg, colorsByCaption, counters);
  if (strict && unsupported.length > 0) throw projectionError(unsupported);
  return Object.freeze({
    svg: projectedSvg,
    changed: Object.values(counters).some((count) => count > 0),
    changedKeys: Object.freeze(changedKeys),
    unsupported: Object.freeze(unsupported.map((entry) => Object.freeze({ ...entry }))),
    counters: Object.freeze({ ...counters })
  });
};

export const createSvgStyles = ({ state, watch, legendActions, previewRuntime = null }) => {
  const {
    svgContent,
    extractedFeatures,
    appliedPaletteColors,
    manualSpecificRules,
    results,
    selectedResultIndex,
    skipCaptureBaseConfig,
    svgContainer,
    adv,
    mode,
    form
  } = state;
  void legendActions;

  const persistSvgEdit = (svg, reason) => {
    skipCaptureBaseConfig.value = true;
    if (previewRuntime?.markActiveResultDirty?.(reason)) return;
    const idx = selectedResultIndex.value;
    if (idx >= 0 && results.value.length > idx) {
      results.value[idx] = { ...results.value[idx], content: serializeCleanSvg(svg) };
    }
  };

  // Compatibility adapter for hydration callers. It deliberately prepares a
  // detached projection and never writes the mounted Result. Live acceptance
  // is owned by the all-Result bulk style command.
  const applyPaletteToSvg = ({
    beforeColors = appliedPaletteColors.value,
    afterColors = appliedPaletteColors.value,
    context = null,
    strict = false
  } = {}) => {
    const svg = svgContainer.value?.querySelector?.('svg');
    if (!svg) return null;
    return prepareNonFeaturePaletteProjection({
      svg,
      beforeColors,
      afterColors,
      context,
      strict
    });
  };

  const applySpecificRulesToSvg = () => {
    if (!svgContent.value || !extractedFeatures.value.length) return;
    if (!manualSpecificRules.length) return;
    if (!svgContainer.value) return;

    const svg = svgContainer.value.querySelector('svg');
    if (!svg) return;

    const featureElementIndex = getFeatureElementIndex(svg);
    let updatedCount = 0;

    extractedFeatures.value.forEach((feat) => {
      if (!feat.svg_id) return;

      let matchingRule = null;

      for (const rule of manualSpecificRules) {
        if ((rule.qual || '').toLowerCase() !== 'hash') continue;
        if (ruleMatchesFeature(feat, rule)) {
          matchingRule = rule;
          break;
        }
      }

      if (!matchingRule) {
        for (const rule of manualSpecificRules) {
          if ((rule.qual || '').toLowerCase() === 'hash') continue;
          if (ruleMatchesFeature(feat, rule)) {
            matchingRule = rule;
            break;
          }
        }
      }

      const elements = getFeatureFillElements(svg, feat.svg_id, featureElementIndex);
      if (elements.length > 0) {
        const newColor = matchingRule ? matchingRule.color : appliedPaletteColors.value[feat.type] || '#cccccc';
        elements.forEach((el) => {
          if (el.getAttribute('fill') !== newColor) {
            el.setAttribute('fill', newColor);
            updatedCount++;
          }
        });
      }
    });

    if (updatedCount > 0) {
      persistSvgEdit(svg, 'specific-rule-style');
      console.log(`Applied specific rules: updated ${updatedCount} elements`);
    }
  };

  const applyStylesToSvg = () => {
    if (!svgContent.value) return;
    if (!svgContainer.value) return;

    const svg = svgContainer.value.querySelector('svg');
    if (!svg) return;

    let updatedCount = 0;

    if (adv.block_stroke_color || adv.block_stroke_width !== null) {
      const featurePaths = Array.from(svg.querySelectorAll(FEATURE_SELECTOR)).filter(isFeatureFillTarget);
      featurePaths.forEach((path) => {
        if (adv.block_stroke_color) {
          path.setAttribute('stroke', adv.block_stroke_color);
          updatedCount++;
        }
        if (adv.block_stroke_width !== null) {
          path.setAttribute('stroke-width', adv.block_stroke_width);
          updatedCount++;
        }
      });
    }

    const axisGroups = getGroupsByBaseIds(svg, ['Axis']);
    axisGroups.forEach((axisGroup) => {
      const axisElements = axisGroup.querySelectorAll('path, line, circle');
      axisElements.forEach((el) => {
        if (adv.axis_stroke_color) {
          el.setAttribute('stroke', adv.axis_stroke_color);
          updatedCount++;
        }
        if (adv.axis_stroke_width !== null) {
          el.setAttribute('stroke-width', adv.axis_stroke_width);
          updatedCount++;
        }
      });
    });

    const tickGroups = getGroupsByBaseIds(svg, ['tick'], ['ticks']);
    tickGroups.forEach((tickGroup) => {
      const tickElements = tickGroup.querySelectorAll('path, line');
      tickElements.forEach((el) => {
        if (adv.axis_stroke_color) {
          el.setAttribute('stroke', adv.axis_stroke_color);
          updatedCount++;
        }
        if (adv.axis_stroke_width !== null) {
          el.setAttribute('stroke-width', adv.axis_stroke_width);
          updatedCount++;
        }
      });
    });

    if (adv.line_stroke_color || adv.line_stroke_width !== null) {
      const connectorPaths = Array.from(svg.querySelectorAll(FEATURE_SELECTOR)).filter(
        (path) => getFeaturePart(path) === FEATURE_PART_CONNECTOR
      );
      connectorPaths.forEach((path) => {
        if (adv.line_stroke_color) {
          path.setAttribute('stroke', adv.line_stroke_color);
          updatedCount++;
        }
        if (adv.line_stroke_width !== null) {
          path.setAttribute('stroke-width', adv.line_stroke_width);
          updatedCount++;
        }
      });
    }

    const lengthBarGroup = svg.getElementById('length_bar');
    if (lengthBarGroup) {
      const scaleElements = lengthBarGroup.querySelectorAll('line, path');
      scaleElements.forEach((el) => {
        if (adv.scale_stroke_color) {
          el.setAttribute('stroke', adv.scale_stroke_color);
          updatedCount++;
        }
        if (adv.scale_stroke_width !== null) {
          el.setAttribute('stroke-width', adv.scale_stroke_width);
          updatedCount++;
        }
      });
    }

    if (adv.block_stroke_color || adv.block_stroke_width !== null) {
      const legendGroups = [];
      const mainLegend = svg.getElementById('legend');
      if (mainLegend) {
        const featureLegend = mainLegend.querySelector('#feature_legend');
        if (featureLegend) legendGroups.push(featureLegend);
        else legendGroups.push(mainLegend);
      }
      const hLegend = svg.getElementById('legend_horizontal');
      if (hLegend) {
        const hFeatureLegend = hLegend.querySelector('#feature_legend_h');
        if (hFeatureLegend) legendGroups.push(hFeatureLegend);
      }
      const vLegend = svg.getElementById('legend_vertical');
      if (vLegend) {
        const vFeatureLegend = vLegend.querySelector('#feature_legend_v');
        if (vFeatureLegend) legendGroups.push(vFeatureLegend);
      }

      legendGroups.forEach((legendGroup) => {
        const paths = legendGroup.querySelectorAll('path');
        paths.forEach((path) => {
          const fill = path.getAttribute('fill');
          const stroke = path.getAttribute('stroke');
          const d = path.getAttribute('d') || '';
          if (
            fill &&
            fill !== 'none' &&
            fill !== 'white' &&
            fill !== '#ffffff' &&
            !fill.startsWith('url(') &&
            stroke &&
            d.includes('z') &&
            d.split(' ').length < 20
          ) {
            if (adv.block_stroke_color) {
              path.setAttribute('stroke', adv.block_stroke_color);
              updatedCount++;
            }
            if (adv.block_stroke_width !== null) {
              path.setAttribute('stroke-width', adv.block_stroke_width);
              updatedCount++;
            }
          }
        });
      });
    }

    if (updatedCount > 0) {
      persistSvgEdit(svg, 'diagram-style');
      console.log(`Applied styles: updated ${updatedCount} elements`);
    }
  };

  const applyTrackVisibility = () => {
    if (!svgContent.value) return;
    if (!svgContainer.value) return;

    const svg = svgContainer.value.querySelector('svg');
    if (!svg) return;

    let updated = false;

    const gcContentGroups = getGroupsByBaseIds(
      svg,
      ['gc_content'],
      ['dinucleotide_content']
    );
    if (gcContentGroups.length > 0) {
      const shouldHide = mode.value === 'circular' ? form.suppress_gc : !form.show_gc;
      gcContentGroups.forEach((gcContentGroup) => {
        const currentDisplay = gcContentGroup.getAttribute('display');
        if (shouldHide && currentDisplay !== 'none') {
          gcContentGroup.setAttribute('display', 'none');
          updated = true;
        } else if (!shouldHide && currentDisplay === 'none') {
          gcContentGroup.removeAttribute('display');
          updated = true;
        }
      });
    }

    const skewGroups = getGroupsByBaseIds(
      svg,
      ['skew', 'gc_skew'],
      ['dinucleotide_skew']
    );
    if (skewGroups.length > 0) {
      const shouldHide = mode.value === 'circular' ? form.suppress_skew : !form.show_skew;
      skewGroups.forEach((skewGroup) => {
        const currentDisplay = skewGroup.getAttribute('display');
        if (shouldHide && currentDisplay !== 'none') {
          skewGroup.setAttribute('display', 'none');
          updated = true;
        } else if (!shouldHide && currentDisplay === 'none') {
          skewGroup.removeAttribute('display');
          updated = true;
        }
      });
    }

    const depthGroups = getGroupsByBaseIds(svg, ['depth'], ['depth']);
    if (depthGroups.length > 0) {
      const shouldHide = !form.show_depth;
      depthGroups.forEach((depthGroup) => {
        const currentDisplay = depthGroup.getAttribute('display');
        if (shouldHide && currentDisplay !== 'none') {
          depthGroup.setAttribute('display', 'none');
          updated = true;
        } else if (!shouldHide && currentDisplay === 'none') {
          depthGroup.removeAttribute('display');
          updated = true;
        }
      });
    }

    if (updated) {
      persistSvgEdit(svg, 'track-visibility');
      console.log('Track visibility updated');
    }
  };

  watch(
    () => [
      adv.block_stroke_color,
      adv.block_stroke_width,
      adv.line_stroke_color,
      adv.line_stroke_width,
      adv.axis_stroke_color,
      adv.axis_stroke_width,
      adv.scale_stroke_color,
      adv.scale_stroke_width
    ],
    () => {
      applyStylesToSvg();
    }
  );

  watch(
    () => [form.suppress_gc, form.suppress_skew, form.show_gc, form.show_skew, form.show_depth],
    () => {
      applyTrackVisibility();
    }
  );

  return {
    applyPaletteToSvg,
    applySpecificRulesToSvg,
    applyStylesToSvg,
    applyTrackVisibility
  };
};
