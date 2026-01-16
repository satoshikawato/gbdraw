import { estimateColorFactor, interpolateColor } from './color-utils.js';
import { ruleMatchesFeature } from './feature-utils.js';

export const createSvgStyles = ({ state, watch, legendActions }) => {
  const {
    svgContent,
    extractedFeatures,
    currentColors,
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

  const ensureUniquePairwiseGradientIds = (svg) => {
    if (!svg) return;

    const legendGroup = svg.getElementById('legend');
    if (!legendGroup) return;

    const horizontalLegend = legendGroup.querySelector('#legend_horizontal');
    const verticalLegend = legendGroup.querySelector('#legend_vertical');

    if (!horizontalLegend || !verticalLegend) return;

    const fixGradientId = (legend, suffix) => {
      const pairwiseLegend = legend.querySelector('#pairwise_legend');
      if (!pairwiseLegend) return;

      const gradient = pairwiseLegend.querySelector('linearGradient');
      if (!gradient) return;

      const currentId = gradient.id;
      if (currentId.endsWith(`_${suffix}`)) return;

      const baseId = currentId.replace(/_[hv]$/, '');
      const newId = `${baseId}_${suffix}`;
      gradient.setAttribute('id', newId);

      const paths = pairwiseLegend.querySelectorAll('path');
      paths.forEach((path) => {
        const fill = path.getAttribute('fill');
        if (fill && fill.includes('url(#')) {
          path.setAttribute('fill', `url(#${newId})`);
        }
      });
    };

    fixGradientId(horizontalLegend, 'h');
    fixGradientId(verticalLegend, 'v');
  };

  const applyPaletteToSvg = () => {
    if (!svgContent.value || !extractedFeatures.value.length) return;
    if (!svgContainer.value) return;

    const svg = svgContainer.value.querySelector('svg');
    if (!svg) return;

    ensureUniquePairwiseGradientIds(svg);

    const colors = currentColors.value;
    const featurePaths = svg.querySelectorAll('path[id^="f"]');
    let updatedCount = 0;

    featurePaths.forEach((path) => {
      const svgId = path.getAttribute('id');
      const feat = extractedFeatures.value.find((f) => f.svg_id === svgId);
      if (!feat) return;

      const paletteColor = colors[feat.type];
      if (!paletteColor) return;

      const hasSpecificRule = manualSpecificRules.some((rule) => ruleMatchesFeature(feat, rule));

      if (!hasSpecificRule && !featureColorOverrides[feat.id]) {
        const currentFill = path.getAttribute('fill');
        if (currentFill !== paletteColor) {
          path.setAttribute('fill', paletteColor);
          updatedCount++;
        }
      }
    });

    const gcContentGroup = svg.getElementById('gc_content');
    if (gcContentGroup && colors.gc_content) {
      const gcPaths = gcContentGroup.querySelectorAll('path');
      gcPaths.forEach((path) => {
        path.setAttribute('fill', colors.gc_content);
        updatedCount++;
      });
    }

    const skewGroup = svg.getElementById('skew') || svg.getElementById('gc_skew');
    if (skewGroup) {
      const skewPaths = skewGroup.querySelectorAll('path');
      let pathIndex = 0;
      skewPaths.forEach((path) => {
        const fill = path.getAttribute('fill');
        if (fill && fill !== 'white' && fill !== 'none') {
          if (pathIndex === 0 && colors.skew_high) {
            path.setAttribute('fill', colors.skew_high);
            updatedCount++;
          } else if (pathIndex === 1 && colors.skew_low) {
            path.setAttribute('fill', colors.skew_low);
            updatedCount++;
          }
          pathIndex++;
        }
      });
    }

    if (colors.pairwise_match_min && colors.pairwise_match_max) {
      let compIdx = 1;
      let compGroup = svg.getElementById(`comparison${compIdx}`);
      while (compGroup) {
        const matchPaths = compGroup.querySelectorAll('path');
        matchPaths.forEach((path, pathIdx) => {
          const pathKey = `comp${compIdx}_path${pathIdx}`;
          const currentFill = path.getAttribute('fill');
          if (currentFill) {
            let factor;
            if (pairwiseMatchFactors.value[pathKey] !== undefined) {
              factor = pairwiseMatchFactors.value[pathKey];
            } else {
              const origMin = window._origPairwiseMin || '#FFE7E7';
              const origMax = window._origPairwiseMax || '#FF7272';
              factor = estimateColorFactor(currentFill, origMin, origMax);
              pairwiseMatchFactors.value[pathKey] = factor;
            }
            const newColor = interpolateColor(colors.pairwise_match_min, colors.pairwise_match_max, factor);
            path.setAttribute('fill', newColor);
            updatedCount++;
          }
        });
        compIdx++;
        compGroup = svg.getElementById(`comparison${compIdx}`);
      }
    }

    const featureLegendGroups = getAllFeatureLegendGroups(svg);
    const keyToColorKey = {
      CDS: 'CDS',
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
              path.setAttribute('fill', newColor);
              break;
            }
          }
        });
      } else {
        const texts = featureLegendGroup.querySelectorAll('text');
        const allPaths = featureLegendGroup.querySelectorAll('path');
        const parseTransform = (transform) => {
          if (!transform) return { x: 0, y: 0 };
          const match = transform.match(/translate\(\s*([\d.-]+)\s*,\s*([\d.-]+)\s*\)/);
          return match ? { x: parseFloat(match[1]), y: parseFloat(match[2]) } : { x: 0, y: 0 };
        };
        texts.forEach((textEl) => {
          const textContent = textEl.textContent?.trim();
          if (!textContent) return;
          if (legendColorOverrides[textContent]) return;

          const newColor = resolveLegendColor(textContent, colors);
          if (!newColor) return;

          const textPos = parseTransform(textEl.getAttribute('transform'));
          let bestPath = null;
          let bestX = -Infinity;
          for (const path of allPaths) {
            const pathPos = parseTransform(path.getAttribute('transform'));
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
            bestPath.setAttribute('fill', newColor);
            updatedCount++;
          }
        });
      }
    });

    if (colors.pairwise_match_min && colors.pairwise_match_max) {
      const allPairwiseLegends = svg.querySelectorAll('#pairwise_legend, [id="pairwise_legend"]');
      allPairwiseLegends.forEach((pairwiseLegend) => {
        const gradient = pairwiseLegend.querySelector('linearGradient');
        if (gradient) {
          const stops = gradient.querySelectorAll('stop');
          if (stops.length >= 2) {
            stops[0].setAttribute('stop-color', colors.pairwise_match_min);
            stops[1].setAttribute('stop-color', colors.pairwise_match_max);
            updatedCount++;
          }
        }
      });
    }

    if (updatedCount > 0) {
      skipCaptureBaseConfig.value = true;
      const idx = selectedResultIndex.value;
      if (idx >= 0 && results.value.length > idx) {
        const serializer = new XMLSerializer();
        results.value[idx] = { ...results.value[idx], content: serializer.serializeToString(svg) };
      }
    }
  };

  const applySpecificRulesToSvg = () => {
    if (!svgContent.value || !extractedFeatures.value.length) return;
    if (!svgContainer.value) return;

    const svg = svgContainer.value.querySelector('svg');
    if (!svg) return;

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

      const elements = svg.querySelectorAll(`#${CSS.escape(feat.svg_id)}`);
      if (elements.length > 0) {
        const newColor = matchingRule ? matchingRule.color : currentColors.value[feat.type] || '#cccccc';
        elements.forEach((el) => {
          if (el.getAttribute('fill') !== newColor) {
            el.setAttribute('fill', newColor);
            updatedCount++;
          }
        });
      }
    });

    if (updatedCount > 0) {
      skipCaptureBaseConfig.value = true;
      const idx = selectedResultIndex.value;
      if (idx >= 0 && results.value.length > idx) {
        const serializer = new XMLSerializer();
        results.value[idx] = { ...results.value[idx], content: serializer.serializeToString(svg) };
      }
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
      const featurePaths = svg.querySelectorAll('path[id^="f"]');
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

    const axisGroup = svg.getElementById('Axis');
    if (axisGroup) {
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
    }

    const tickGroup = svg.getElementById('tick');
    if (tickGroup) {
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
    }

    if (adv.line_stroke_color || adv.line_stroke_width !== null) {
      const allPaths = svg.querySelectorAll('path');
      allPaths.forEach((path) => {
        const fill = path.getAttribute('fill');
        const stroke = path.getAttribute('stroke');
        const id = path.getAttribute('id') || '';
        if ((fill === 'none' || !fill) && stroke && !id.startsWith('f')) {
          if (adv.line_stroke_color) {
            path.setAttribute('stroke', adv.line_stroke_color);
            updatedCount++;
          }
          if (adv.line_stroke_width !== null) {
            path.setAttribute('stroke-width', adv.line_stroke_width);
            updatedCount++;
          }
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
      skipCaptureBaseConfig.value = true;
      const idx = selectedResultIndex.value;
      if (idx >= 0 && results.value.length > idx) {
        const serializer = new XMLSerializer();
        results.value[idx] = { ...results.value[idx], content: serializer.serializeToString(svg) };
      }
      console.log(`Applied styles: updated ${updatedCount} elements`);
    }
  };

  const applyTrackVisibility = () => {
    if (!svgContent.value) return;
    if (!svgContainer.value) return;

    const svg = svgContainer.value.querySelector('svg');
    if (!svg) return;

    let updated = false;

    const gcContentGroups = svg.querySelectorAll('#gc_content, [id="gc_content"]');
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

    const skewGroups = svg.querySelectorAll('#skew, [id="skew"], #gc_skew, [id="gc_skew"]');
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

    if (updated) {
      skipCaptureBaseConfig.value = true;
      const idx = selectedResultIndex.value;
      if (idx >= 0 && results.value.length > idx) {
        const serializer = new XMLSerializer();
        results.value[idx] = { ...results.value[idx], content: serializer.serializeToString(svg) };
      }
      console.log('Track visibility updated');
    }
  };

  watch(
    currentColors,
    () => {
      applyPaletteToSvg();
      applySpecificRulesToSvg();
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
    () => [form.suppress_gc, form.suppress_skew, form.show_gc, form.show_skew],
    () => {
      applyTrackVisibility();
    }
  );

  return {
    ensureUniquePairwiseGradientIds,
    applyPaletteToSvg,
    applySpecificRulesToSvg,
    applyStylesToSvg,
    applyTrackVisibility
  };
};
