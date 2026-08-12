import {
  estimateColorFactor,
  interpolateColor,
  resolveCollinearMatchColor,
  resolvePairwiseLegendGradientColorKeys
} from './color-utils.js';
import {
  FEATURE_SELECTOR,
  getFeatureElementIndex,
  getFeatureFillElements,
  getFeatureIdentity
} from './feature-editor/svg-actions.js';
import {
  FEATURE_PART_CONNECTOR,
  getFeaturePart,
  isFeatureFillTarget
} from './feature-dom.js';
import { ruleMatchesFeature } from './feature-utils.js';
import { PAIRWISE_LEGEND_SELECTOR, parseTransformXY } from './legend/utils.js';
import { serializeCleanSvg } from '../services/svg-serialization.js';
import { getFeatureOverride } from '../services/feature-override-identity.js';
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

const paletteColorsEqual = (left, right) => {
  const keys = new Set([...Object.keys(left || {}), ...Object.keys(right || {})]);
  return Array.from(keys).every(
    (key) => normalizeComparableColor(left?.[key]) === normalizeComparableColor(right?.[key])
  );
};

const paletteColorKeysEqual = (left, right, keys) => keys.every(
  (key) => normalizeComparableColor(left?.[key]) === normalizeComparableColor(right?.[key])
);

export const createSvgStyles = ({ state, watch, nextTick, legendActions, previewRuntime = null }) => {
  const {
    svgContent,
    extractedFeatures,
    featuresBySvgId,
    appliedPaletteColors,
    manualSpecificRules,
    featureColorOverrides,
    legendColorOverrides,
    pairwiseMatchFactors,
    results,
    selectedResultIndex,
    skipCaptureBaseConfig,
    svgContainer,
    adv,
    mode,
    form
  } = state;

  const { getAllFeatureLegendGroups } = legendActions;

  const persistSvgEdit = (svg, reason) => {
    skipCaptureBaseConfig.value = true;
    if (previewRuntime?.markActiveResultDirty?.(reason)) return;
    const idx = selectedResultIndex.value;
    if (idx >= 0 && results.value.length > idx) {
      results.value[idx] = { ...results.value[idx], content: serializeCleanSvg(svg) };
    }
  };

  const updatePairwiseLegendGradientStops = (pairwiseLegend, colors) => {
    let updated = false;
    pairwiseLegend.querySelectorAll('linearGradient').forEach((gradient) => {
      const legendKey = gradient.closest('g[data-legend-key]')?.getAttribute('data-legend-key') || '';
      const { minKey, maxKey } = resolvePairwiseLegendGradientColorKeys(legendKey);
      const minColor = colors[minKey];
      const maxColor = colors[maxKey];
      if (!minColor || !maxColor) return;
      const stops = gradient.querySelectorAll('stop');
      if (stops.length >= 2) {
        updated = setColorAttributeIfChanged(stops[0], 'stop-color', minColor) || updated;
        updated = setColorAttributeIfChanged(stops[1], 'stop-color', maxColor) || updated;
      }
    });
    return updated;
  };

  const applyPaletteToSvg = ({
    recolorPairwise = false,
    recolorCollinear = false
  } = {}) => {
    if (!svgContent.value || !extractedFeatures.value.length) return;
    if (!svgContainer.value) return;

    const svg = svgContainer.value.querySelector('svg');
    if (!svg) return;

    const colors = appliedPaletteColors.value;
    const featurePaths = Array.from(getFeatureElementIndex(svg).values()).flat();
    const featureLookup = featuresBySvgId?.value || new Map();
    let updatedCount = 0;

    featurePaths.forEach((path) => {
      if (!isFeatureFillTarget(path)) return;
      const svgId = getFeatureIdentity(path);
      const feat = featureLookup.get(svgId);
      if (!feat) return;

      const paletteColor = colors[feat.type];
      if (!paletteColor) return;

      const hasSpecificRule = manualSpecificRules.some((rule) => ruleMatchesFeature(feat, rule));

      if (!hasSpecificRule && !getFeatureOverride(featureColorOverrides, feat)) {
        const currentFill = path.getAttribute('fill');
        if (currentFill !== paletteColor) {
          path.setAttribute('fill', paletteColor);
          updatedCount++;
        }
      }
    });

    const gcContentGroups = getGroupsByBaseIds(
      svg,
      ['gc_content'],
      ['dinucleotide_content']
    );
    if (gcContentGroups.length > 0 && colors.gc_content) {
      gcContentGroups.forEach((gcContentGroup) => {
        const gcPaths = gcContentGroup.querySelectorAll('path');
        gcPaths.forEach((path) => {
          if (setColorAttributeIfChanged(path, 'fill', colors.gc_content)) updatedCount++;
        });
      });
    }

    const skewGroups = getGroupsByBaseIds(
      svg,
      ['skew', 'gc_skew'],
      ['dinucleotide_skew']
    );
    if (skewGroups.length > 0) {
      skewGroups.forEach((skewGroup) => {
        const skewPaths = skewGroup.querySelectorAll('path');
        let pathIndex = 0;
        skewPaths.forEach((path) => {
          const fill = path.getAttribute('fill');
          if (fill && fill !== 'white' && fill !== 'none') {
            if (pathIndex === 0 && colors.skew_high) {
              if (setColorAttributeIfChanged(path, 'fill', colors.skew_high)) updatedCount++;
            } else if (pathIndex === 1 && colors.skew_low) {
              if (setColorAttributeIfChanged(path, 'fill', colors.skew_low)) updatedCount++;
            }
            pathIndex++;
          }
        });
      });
    }

    if (
      (recolorPairwise || recolorCollinear)
      && colors.pairwise_match_min
      && colors.pairwise_match_max
    ) {
      let compIdx = 1;
      let compGroup = svg.getElementById(`comparison${compIdx}`);
      while (compGroup) {
        const matchPaths = compGroup.querySelectorAll('path');
        matchPaths.forEach((path, pathIdx) => {
          const pathKey = `comp${compIdx}_path${pathIdx}`;
          const currentFill = path.getAttribute('fill');
          if (currentFill) {
            const collinearityBlockId = path.getAttribute('data-collinearity-block-id') || '';
            const collinearityColorMode = path.getAttribute('data-collinearity-color-mode') || '';
            const metadataText = path.getAttribute('data-identity-factor');
            const metadataFactor = metadataText === null || metadataText === '' ? NaN : Number(metadataText);
            const collinearColor = resolveCollinearMatchColor({
              blockId: collinearityBlockId,
              colorMode: collinearityColorMode,
              orientation: path.getAttribute('data-collinearity-orientation') || '',
              identityFactor: Number.isFinite(metadataFactor) ? metadataFactor : null,
              colors
            });
            if (collinearColor) {
              if (
                recolorCollinear
                && setColorAttributeIfChanged(path, 'fill', collinearColor)
              ) {
                updatedCount++;
              }
              return;
            }
            if (collinearityBlockId && !collinearityColorMode) return;
            if (!recolorPairwise) return;

            let factor;
            if (Number.isFinite(metadataFactor)) {
              factor = metadataFactor;
              pairwiseMatchFactors.value[pathKey] = factor;
            } else if (pairwiseMatchFactors.value[pathKey] !== undefined) {
              factor = pairwiseMatchFactors.value[pathKey];
            } else {
              const origMin = window._origPairwiseMin || '#FFE7E7';
              const origMax = window._origPairwiseMax || '#FF7272';
              factor = estimateColorFactor(currentFill, origMin, origMax);
              pairwiseMatchFactors.value[pathKey] = factor;
            }
            const newColor = interpolateColor(colors.pairwise_match_min, colors.pairwise_match_max, factor);
            if (setColorAttributeIfChanged(path, 'fill', newColor)) updatedCount++;
          }
        });
        compIdx++;
        compGroup = svg.getElementById(`comparison${compIdx}`);
      }
    }

    const featureLegendGroups = getAllFeatureLegendGroups(svg);
    const keyToColorKey = {
      CDS: 'CDS',
      'D-loop': 'D-loop',
      repeat_region: 'repeat_region',
      tmRNA: 'tmRNA',
      tRNA: 'tRNA',
      rRNA: 'rRNA',
      ncRNA: 'ncRNA',
      misc_feature: 'misc_feature',
      mobile_element: 'mobile_element',
      'GC content': 'gc_content',
      'GC skew (+)': 'skew_high',
      'GC skew (-)': 'skew_low'
    };
    const resolveOtherLegendColor = (legendKey, palette) => {
      if (!legendKey) return null;
      const lowerKey = legendKey.toLowerCase();
      if (lowerKey === 'other proteins') return palette.CDS || null;
      if (!lowerKey.startsWith('other ')) return null;
      let raw = legendKey.slice(6).trim();
      if (!raw) return null;
      if (raw.toLowerCase() === 'proteins') return palette.CDS || null;
      if (raw.endsWith('s')) raw = raw.slice(0, -1);
      return palette[raw] || null;
    };
    const resolveLegendColor = (legendKey, palette) => {
      if (!legendKey) return null;
      const colorKey = keyToColorKey[legendKey];
      if (colorKey && palette[colorKey]) return palette[colorKey];
      if (palette[legendKey]) return palette[legendKey];
      return resolveOtherLegendColor(legendKey, palette);
    };
    featureLegendGroups.forEach((featureLegendGroup) => {
      if (!featureLegendGroup) return;

      const entryGroups = featureLegendGroup.querySelectorAll('g[data-legend-key]');

      if (entryGroups.length > 0) {
        entryGroups.forEach((entryGroup) => {
          const legendKey = entryGroup.getAttribute('data-legend-key');
          if (!legendKey) return;
          if (legendColorOverrides[legendKey]) return;

          const newColor = resolveLegendColor(legendKey, colors);
          if (!newColor) return;

          const paths = entryGroup.querySelectorAll('path');
          for (const path of paths) {
            const fill = path.getAttribute('fill');
            if (fill && fill !== 'none' && !fill.startsWith('url(')) {
              if (setColorAttributeIfChanged(path, 'fill', newColor)) updatedCount++;
              break;
            }
          }
        });
      } else {
        const texts = featureLegendGroup.querySelectorAll('text');
        const allPaths = featureLegendGroup.querySelectorAll('path');
        texts.forEach((textEl) => {
          const textContent = textEl.textContent?.trim();
          if (!textContent) return;
          if (legendColorOverrides[textContent]) return;

          const newColor = resolveLegendColor(textContent, colors);
          if (!newColor) return;

          const textPos = parseTransformXY(textEl.getAttribute('transform'));
          let bestPath = null;
          let bestX = -Infinity;
          for (const path of allPaths) {
            const pathPos = parseTransformXY(path.getAttribute('transform'));
            const fill = path.getAttribute('fill');
            if (
              Math.abs(pathPos.y - textPos.y) < 2 &&
              pathPos.x < textPos.x &&
              fill &&
              fill !== 'none' &&
              !fill.startsWith('url(')
            ) {
              if (pathPos.x > bestX) {
                bestX = pathPos.x;
                bestPath = path;
              }
            }
          }
          if (bestPath) {
            if (setColorAttributeIfChanged(bestPath, 'fill', newColor)) updatedCount++;
          }
        });
      }
    });

    if (colors.pairwise_match_min && colors.pairwise_match_max) {
      const allPairwiseLegends = svg.querySelectorAll(PAIRWISE_LEGEND_SELECTOR);
      allPairwiseLegends.forEach((pairwiseLegend) => {
        if (updatePairwiseLegendGradientStops(pairwiseLegend, colors)) updatedCount++;
      });
    }

    if (updatedCount > 0) {
      persistSvgEdit(svg, 'palette-style');
    }
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
    appliedPaletteColors,
    (colors, previousColors) => {
      if (paletteColorsEqual(colors, previousColors)) return;
      // Inferred comparison factors are lossy, so only re-interpolate a family
      // when one of its palette endpoints actually changed.
      const recolorPairwise = !paletteColorKeysEqual(
        colors,
        previousColors,
        ['pairwise_match_min', 'pairwise_match_max']
      );
      const recolorCollinear = !paletteColorKeysEqual(
        colors,
        previousColors,
        [
          'collinear_block_plus_min',
          'collinear_block_plus',
          'collinear_block_minus_min',
          'collinear_block_minus'
        ]
      );
      nextTick(() => {
        applyPaletteToSvg({ recolorPairwise, recolorCollinear });
        applySpecificRulesToSvg();
      });
    },
    { deep: true }
  );

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
