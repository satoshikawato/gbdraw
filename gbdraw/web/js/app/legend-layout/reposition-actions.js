import { estimateColorFactor, interpolateColor } from '../color-utils.js';
import { getElementsBounds, getTransformedBBox } from './transform-utils.js';

export const createLegendRepositionActions = ({
  state,
  debugLog,
  legendActions,
  svgActions,
  diagramActions
}) => {
  const {
    svgContent,
    svgContainer,
    mode,
    generatedLegendPosition,
    circularBaseConfig,
    linearBaseConfig,
    diagramElements,
    diagramElementOriginalTransforms,
    diagramElementBaseTransforms,
    diagramOffset,
    legendInitialTransform,
    legendCurrentOffset,
    pairwiseMatchFactors,
    currentColors,
    legendColorOverrides,
    selectedResultIndex,
    results,
    skipCaptureBaseConfig,
    skipPositionReapply,
    form
  } = state;

  const { getAllFeatureLegendGroups, reflowDualLegendLayout, reflowSingleLegendLayout } = legendActions;
  const { ensureUniquePairwiseGradientIds } = svgActions;
  const { applyDiagramShift } = diagramActions;

  const getCircularAbsoluteConfig = (position, baseConfig) => {
    const { viewBoxWidth, viewBoxHeight, legendWidth, legendHeight } = baseConfig;
    const genVbW = baseConfig.generatedViewBoxWidth || viewBoxWidth;
    const genVbH = baseConfig.generatedViewBoxHeight || viewBoxHeight;
    console.log(
      `[DEBUG] getCircularAbsoluteConfig called with position="${position}" (type: ${typeof position}, length: ${position.length})`
    );
    console.log(`[DEBUG] position === 'left': ${position === 'left'}, position === 'right': ${position === 'right'}`);
    console.log(`[DEBUG] viewBoxWidth=${viewBoxWidth}, genVbW=${genVbW}`);

    switch (position) {
      case 'left':
        debugLog('Matched case: left');
        return {
          vbWidth: viewBoxWidth + legendWidth * 1.1,
          vbHeight: viewBoxHeight,
          diagramShiftX: legendWidth * 0.55,
          diagramShiftY: 0,
          legendX: 10,
          legendY: (viewBoxHeight - legendHeight) / 2
        };
      case 'right':
        debugLog('Matched case: right');
        return {
          vbWidth: viewBoxWidth + legendWidth * 1.1,
          vbHeight: viewBoxHeight,
          diagramShiftX: 0,
          diagramShiftY: 0,
          legendX: viewBoxWidth + legendWidth * 0.05,
          legendY: (viewBoxHeight - legendHeight) / 2
        };
      case 'upper_left':
        return {
          vbWidth: genVbW,
          vbHeight: genVbH,
          diagramShiftX: 0,
          diagramShiftY: 0,
          legendX: 0.025 * genVbW,
          legendY: 0.05 * genVbH
        };
      case 'upper_right':
        return {
          vbWidth: genVbW,
          vbHeight: genVbH,
          diagramShiftX: 0,
          diagramShiftY: 0,
          legendX: 0.85 * genVbW,
          legendY: 0.05 * genVbH
        };
      case 'lower_left':
        return {
          vbWidth: genVbW,
          vbHeight: genVbH,
          diagramShiftX: 0,
          diagramShiftY: 0,
          legendX: 0.025 * genVbW,
          legendY: 0.78 * genVbH
        };
      case 'lower_right':
        return {
          vbWidth: genVbW,
          vbHeight: genVbH,
          diagramShiftX: 0,
          diagramShiftY: 0,
          legendX: 0.875 * genVbW,
          legendY: 0.75 * genVbH
        };
      case 'none':
      default:
        console.log(`[DEBUG] Matched case: default (position was "${position}")`);
        return {
          vbWidth: genVbW,
          vbHeight: genVbH,
          diagramShiftX: 0,
          diagramShiftY: 0,
          legendX: 0,
          legendY: 0
        };
    }
  };

  const repositionForLegendChange = (newPosition, oldPosition) => {
    if (!svgContainer.value || !svgContent.value) return;
    const svg = svgContainer.value.querySelector('svg');
    if (!svg) return;

    const legendGroup = svg.getElementById('legend');
    const isLinear = mode.value === 'linear';

    let legendWidth = 120;
    let legendHeight = 150;
    if (legendGroup) {
      const bbox = legendGroup.getBBox();
      legendWidth = bbox.width || legendWidth;
      legendHeight = bbox.height || legendHeight;
    }

    if (!isLinear) {
      debugLog('Circular mode reposition:', { newPosition, oldPosition });
      debugLog('circularBaseConfig:', circularBaseConfig.value);

      if (legendGroup) {
        if (newPosition === 'none') {
          legendGroup.style.display = 'none';
        } else {
          legendGroup.style.display = '';
        }
      }

      if (newPosition === 'none') {
        generatedLegendPosition.value = newPosition;
        return;
      }

      const baseConfig = {
        ...circularBaseConfig.value,
        legendWidth: legendWidth,
        legendHeight: legendHeight
      };
      const origLegendWidth = circularBaseConfig.value.legendWidth || 0;
      const origLegendHeight = circularBaseConfig.value.legendHeight || 0;
      if (Math.abs(legendWidth - origLegendWidth) > 1 || Math.abs(legendHeight - origLegendHeight) > 1) {
        console.log(
          `[DEBUG] WARN Legend size CHANGED: ${origLegendWidth.toFixed(1)}x${origLegendHeight.toFixed(
            1
          )} -> ${legendWidth.toFixed(1)}x${legendHeight.toFixed(1)}`
        );
        console.log(
          `[DEBUG] WARN Shift calculation affected: old=${(origLegendWidth * 0.55).toFixed(1)}, new=${(
            legendWidth * 0.55
          ).toFixed(1)}`
        );
      }
      debugLog('baseConfig:', baseConfig);

      const config = getCircularAbsoluteConfig(newPosition, baseConfig);
      debugLog('calculated config:', config);

      svg.setAttribute('viewBox', `0 0 ${config.vbWidth} ${config.vbHeight}`);

      if (legendGroup) {
        legendGroup.setAttribute('transform', `translate(${config.legendX}, ${config.legendY})`);
        legendInitialTransform.value = { x: config.legendX, y: config.legendY };
        legendCurrentOffset.x = 0;
        legendCurrentOffset.y = 0;
      }

      debugLog('Circular: diagramElementBaseTransforms size:', diagramElementBaseTransforms.value.size);
      debugLog('Circular: diagramElements count:', diagramElements.value.length);

      debugLog('Base transforms at reposition start:');
      for (const [el, transform] of diagramElementBaseTransforms.value) {
        console.log(`  ${el.id}: {x:${transform.x}, y:${transform.y}}`);
      }

      if (diagramElements.value.length > 0) {
        diagramElements.value.forEach((el, idx) => {
          const currentTransform = el.getAttribute('transform');
          const baseTransform = diagramElementBaseTransforms.value.get(el) || { x: 0, y: 0 };
          const foundInMap = diagramElementBaseTransforms.value.has(el);
          const newX = baseTransform.x + config.diagramShiftX;
          const newY = baseTransform.y + config.diagramShiftY;
          console.log(
            `[DEBUG] Element ${idx} (${el.id}): currentTransform="${currentTransform}", foundInMap=${foundInMap}, baseTransform={x:${
              baseTransform.x
            }, y:${baseTransform.y}}, newTransform={x:${newX}, y:${newY}}`
          );
          el.setAttribute('transform', `translate(${newX}, ${newY})`);

          diagramElementOriginalTransforms.value.set(el, { x: newX, y: newY });
          console.log(`[DEBUG] repositionForLegendChange UPDATED originalTransforms for ${el.id}: (${newX}, ${newY})`);
        });
        diagramOffset.x = 0;
        diagramOffset.y = 0;
        console.log(`[DEBUG] repositionForLegendChange FINISHED updating ${diagramElements.value.length} elements`);
      }
    } else {
      debugLog('Linear mode reposition:', { newPosition, oldPosition });

      let viewBox = svg.getAttribute('viewBox');
      if (!viewBox) return;
      const parts = viewBox.split(/\s+/).map(parseFloat);
      if (parts.length !== 4) return;
      let [vbX, vbY, vbW, vbH] = parts;
      debugLog('viewBox:', { vbX, vbY, vbW, vbH });

      const isHorizontalLayout = (pos) => pos === 'top' || pos === 'bottom';

      const generatedPos = linearBaseConfig.value.generatedPosition || 'right';
      const wasHorizontal = isHorizontalLayout(generatedPos);
      const nowHorizontal = isHorizontalLayout(newPosition);
      debugLog('layout:', { generatedPos, wasHorizontal, nowHorizontal });
      const prevLegendPos = oldPosition || generatedLegendPosition.value || generatedPos;
      const wasHorizontalPos = isHorizontalLayout(prevLegendPos);
      const prevHorizontalLegendHeight = linearBaseConfig.value.horizontalLegendHeight || legendHeight;

      if (legendGroup) {
        if (newPosition === 'none') {
          legendGroup.setAttribute('display', 'none');
        } else {
          legendGroup.removeAttribute('display');

          const horizontalLegend = legendGroup.querySelector('#legend_horizontal');
          const verticalLegend = legendGroup.querySelector('#legend_vertical');
          const hasDualLegends = !!(horizontalLegend && verticalLegend);
          debugLog('dual legends found:', { horizontal: !!horizontalLegend, vertical: !!verticalLegend });

          if (!nowHorizontal && wasHorizontalPos) {
            legendGroup.setAttribute('transform', 'translate(0, 0)');
          }

          if (hasDualLegends) {
            debugLog('Horizontal legend innerHTML (first 500 chars):', horizontalLegend.innerHTML.substring(0, 500));
            debugLog('Vertical legend innerHTML (first 500 chars):', verticalLegend.innerHTML.substring(0, 500));
            debugLog('Before toggle - horizontal display:', horizontalLegend.getAttribute('display'));
            debugLog('Before toggle - vertical display:', verticalLegend.getAttribute('display'));

            horizontalLegend.removeAttribute('display');
            verticalLegend.removeAttribute('display');

            const hBbox = horizontalLegend.getBBox();
            const vBbox = verticalLegend.getBBox();
            debugLog('Horizontal legend bbox:', { width: hBbox.width, height: hBbox.height });
            debugLog('Vertical legend bbox:', { width: vBbox.width, height: vBbox.height });
            const hLegendWidth = hBbox.width || linearBaseConfig.value.horizontalLegendWidth || legendWidth;
            const hLegendHeight = hBbox.height || linearBaseConfig.value.horizontalLegendHeight || legendHeight;
            const vLegendWidth = vBbox.width || linearBaseConfig.value.verticalLegendWidth || legendWidth;
            const vLegendHeight = vBbox.height || linearBaseConfig.value.verticalLegendHeight || legendHeight;
            linearBaseConfig.value.horizontalLegendWidth = hLegendWidth;
            linearBaseConfig.value.horizontalLegendHeight = hLegendHeight;
            linearBaseConfig.value.verticalLegendWidth = vLegendWidth;
            linearBaseConfig.value.verticalLegendHeight = vLegendHeight;

            if (nowHorizontal) {
              verticalLegend.setAttribute('display', 'none');
              legendWidth = hLegendWidth;
              legendHeight = hLegendHeight;
              debugLog('Showing HORIZONTAL legend, using dimensions:', { legendWidth, legendHeight });
            } else {
              horizontalLegend.setAttribute('display', 'none');
              legendWidth = vLegendWidth;
              legendHeight = vLegendHeight;
              debugLog('Showing VERTICAL legend, using dimensions:', { legendWidth, legendHeight });
            }

            debugLog('After toggle - horizontal display:', horizontalLegend.getAttribute('display'));
            debugLog('After toggle - vertical display:', verticalLegend.getAttribute('display'));

            ensureUniquePairwiseGradientIds(svg);

            if (currentColors.value.pairwise_match_min && currentColors.value.pairwise_match_max) {
              [horizontalLegend, verticalLegend].forEach((legend, idx) => {
                if (!legend) return;
                const pairwiseLegend = legend.querySelector('#pairwise_legend');
                if (pairwiseLegend) {
                  const gradient = pairwiseLegend.querySelector('linearGradient');
                  if (gradient) {
                    const stops = gradient.querySelectorAll('stop');
                    if (stops.length >= 2) {
                      stops[0].setAttribute('stop-color', currentColors.value.pairwise_match_min);
                      stops[1].setAttribute('stop-color', currentColors.value.pairwise_match_max);
                      console.log(
                        `[DEBUG] Updated ${idx === 0 ? 'horizontal' : 'vertical'} pairwise gradient:`,
                        currentColors.value.pairwise_match_min,
                        '->',
                        currentColors.value.pairwise_match_max
                      );
                    }
                  }
                }
              });
            }
          }

          const vbAttr = nowHorizontal
            ? svg.getAttribute('data-horizontal-viewbox')
            : svg.getAttribute('data-vertical-viewbox');
          debugLog('Applying viewBox attribute:', vbAttr);
          if (vbAttr) {
            const vbParts = vbAttr.split(/\s+/).map(parseFloat);
            if (vbParts.length === 4) {
              [vbX, vbY, vbW, vbH] = vbParts;
              svg.setAttribute('viewBox', `${vbX} ${vbY} ${vbW} ${vbH}`);
              debugLog('ViewBox set to:', { vbX, vbY, vbW, vbH });
            }
          } else {
            const configVb = nowHorizontal
              ? linearBaseConfig.value.horizontalViewBox
              : linearBaseConfig.value.verticalViewBox;
            debugLog('Using config viewBox:', { nowHorizontal, configVb });
            if (configVb && configVb.w > 0 && configVb.h > 0) {
              vbX = configVb.x;
              vbY = configVb.y;
              vbW = configVb.w;
              vbH = configVb.h;
              svg.setAttribute('viewBox', `${vbX} ${vbY} ${vbW} ${vbH}`);
              debugLog('ViewBox set from config:', { vbX, vbY, vbW, vbH });
            }
          }

          if (hasDualLegends) {
            reflowDualLegendLayout(svg);

            const hDisplay = horizontalLegend.getAttribute('display');
            const vDisplay = verticalLegend.getAttribute('display');
            horizontalLegend.removeAttribute('display');
            verticalLegend.removeAttribute('display');

            const hBbox = horizontalLegend.getBBox();
            const vBbox = verticalLegend.getBBox();

            if (hDisplay !== null) {
              horizontalLegend.setAttribute('display', hDisplay);
            } else {
              horizontalLegend.removeAttribute('display');
            }
            if (vDisplay !== null) {
              verticalLegend.setAttribute('display', vDisplay);
            } else {
              verticalLegend.removeAttribute('display');
            }

            const hLegendWidth = hBbox.width || linearBaseConfig.value.horizontalLegendWidth || legendWidth;
            const hLegendHeight = hBbox.height || linearBaseConfig.value.horizontalLegendHeight || legendHeight;
            const vLegendWidth = vBbox.width || linearBaseConfig.value.verticalLegendWidth || legendWidth;
            const vLegendHeight = vBbox.height || linearBaseConfig.value.verticalLegendHeight || legendHeight;

            linearBaseConfig.value.horizontalLegendWidth = hLegendWidth;
            linearBaseConfig.value.horizontalLegendHeight = hLegendHeight;
            linearBaseConfig.value.verticalLegendWidth = vLegendWidth;
            linearBaseConfig.value.verticalLegendHeight = vLegendHeight;

            if (nowHorizontal) {
              legendWidth = hLegendWidth;
              legendHeight = hLegendHeight;
            } else {
              legendWidth = vLegendWidth;
              legendHeight = vLegendHeight;
            }

            const updatedViewBox = svg.getAttribute('viewBox');
            if (updatedViewBox) {
              const updatedParts = updatedViewBox.split(/\s+/).map(parseFloat);
              if (updatedParts.length === 4) {
                [vbX, vbY, vbW, vbH] = updatedParts;
              }
            }
          }

          if (!hasDualLegends) {
            const layout = nowHorizontal ? 'horizontal' : 'vertical';
            let maxWidthOverride = vbW;
            if (layout === 'horizontal') {
              const wrapAttr = svg.getAttribute('data-horizontal-wrap-width');
              if (wrapAttr) {
                const wrapWidth = parseFloat(wrapAttr);
                if (Number.isFinite(wrapWidth) && wrapWidth > 0) {
                  maxWidthOverride = wrapWidth;
                }
              }
            }
            const reflowResult = reflowSingleLegendLayout(svg, layout, maxWidthOverride);
            if (reflowResult) {
              legendWidth = reflowResult.legendWidth || legendWidth;
              legendHeight = reflowResult.legendHeight || legendHeight;
            }
          }

          const padding = 20;

          if (nowHorizontal) {
            linearBaseConfig.value.horizontalLegendWidth = legendWidth;
            linearBaseConfig.value.horizontalLegendHeight = legendHeight;

            const heightDelta = legendHeight - prevHorizontalLegendHeight;
            if (heightDelta > 0) {
              vbH += heightDelta;
            }

            svg.setAttribute('viewBox', `${vbX} ${vbY} ${vbW} ${vbH}`);
            svg.setAttribute('data-horizontal-viewbox', `${vbX} ${vbY} ${vbW} ${vbH}`);
            linearBaseConfig.value.horizontalViewBox = { x: vbX, y: vbY, w: vbW, h: vbH };
          } else {
            linearBaseConfig.value.verticalLegendWidth = legendWidth;
            linearBaseConfig.value.verticalLegendHeight = legendHeight;
          }

          const legendPadding = 40;
          const requiredHeight = legendHeight + legendPadding;
          if (!nowHorizontal && requiredHeight > vbH) {
            const newVbH = requiredHeight;
            debugLog('Expanding canvas height:', { oldHeight: vbH, newHeight: newVbH });

            vbH = newVbH;
            svg.setAttribute('viewBox', `${vbX} ${vbY} ${vbW} ${vbH}`);
            svg.setAttribute('data-vertical-viewbox', `${vbX} ${vbY} ${vbW} ${vbH}`);
            debugLog('ViewBox after height expansion:', { vbX, vbY, vbW, vbH });
            linearBaseConfig.value.verticalViewBox = { x: vbX, y: vbY, w: vbW, h: vbH };
          } else if (!nowHorizontal) {
            svg.setAttribute('data-vertical-viewbox', `${vbX} ${vbY} ${vbW} ${vbH}`);
            linearBaseConfig.value.verticalViewBox = { x: vbX, y: vbY, w: vbW, h: vbH };
          }

          const computeVerticalLegendY = () => {
            const ignoredIds = new Set(['length_bar', 'tick', 'labels', 'Axis']);
            const legendElements = diagramElements.value.filter((el) => !ignoredIds.has(el.id));
            const bounds = getElementsBounds(legendElements.length > 0 ? legendElements : diagramElements.value);
            if (!bounds) {
              return vbY + (vbH - legendHeight) / 2;
            }

            const centerY = bounds.y + bounds.height / 2;
            let legendY = centerY - legendHeight / 2;
            const viewportBottom = vbY + vbH;

            if (legendY < vbY || legendY + legendHeight > viewportBottom) {
              legendY = Math.max(vbY + (vbH - legendHeight) / 2, vbY + padding);
            }

            return legendY;
          };

          let finalX;
          let finalY;

          switch (newPosition) {
            case 'top':
              finalX = vbX + (vbW - legendWidth) / 2;
              finalY = vbY + padding;
              break;
            case 'bottom':
              finalX = vbX + (vbW - legendWidth) / 2;
              finalY = vbY + vbH - legendHeight - padding;
              break;
            case 'left':
              finalX = vbX + padding;
              finalY = computeVerticalLegendY();
              break;
            case 'right':
            default:
              finalX = vbX + vbW - legendWidth - padding;
              finalY = computeVerticalLegendY();
              break;
          }

          legendGroup.setAttribute('transform', `translate(${finalX}, ${finalY})`);
          legendInitialTransform.value = { x: finalX, y: finalY };
          legendCurrentOffset.x = 0;
          legendCurrentOffset.y = 0;
        }
      }

      const hVbAttr = svg.getAttribute('data-horizontal-viewbox');
      const vVbAttr = svg.getAttribute('data-vertical-viewbox');
      const hVb = hVbAttr
        ? hVbAttr.split(/\s+/).map(parseFloat)
        : linearBaseConfig.value.horizontalViewBox?.w > 0
          ? [
              linearBaseConfig.value.horizontalViewBox.x,
              linearBaseConfig.value.horizontalViewBox.y,
              linearBaseConfig.value.horizontalViewBox.w,
              linearBaseConfig.value.horizontalViewBox.h
            ]
          : null;
      const vVb = vVbAttr
        ? vVbAttr.split(/\s+/).map(parseFloat)
        : linearBaseConfig.value.verticalViewBox?.w > 0
          ? [
              linearBaseConfig.value.verticalViewBox.x,
              linearBaseConfig.value.verticalViewBox.y,
              linearBaseConfig.value.verticalViewBox.w,
              linearBaseConfig.value.verticalViewBox.h
            ]
          : null;
      const canRepositionDiagram = diagramElements.value.length > 0 && diagramElementBaseTransforms.value.size > 0;

      if (canRepositionDiagram) {

        let shiftX = 0;
        let shiftY = 0;

        if (hVb && vVb) {
          const vLegendWidth = linearBaseConfig.value.verticalLegendWidth || legendWidth;
          const hLegendHeight = linearBaseConfig.value.horizontalLegendHeight || legendHeight;

          const generatedWithLeftLegend = generatedPos === 'left';
          const needsLeftLegend = newPosition === 'left';

          if (!generatedWithLeftLegend && needsLeftLegend) {
            shiftX = vLegendWidth;
          } else if (generatedWithLeftLegend && !needsLeftLegend) {
            shiftX = -vLegendWidth;
          }

          const generatedWithTopLegend = generatedPos === 'top';
          const needsTopLegend = newPosition === 'top';

          if (!generatedWithTopLegend && needsTopLegend) {
            shiftY = hLegendHeight;
          } else if (generatedWithTopLegend && !needsTopLegend) {
            shiftY = -hLegendHeight;
          }
        }

        debugLog('Linear: diagramElementBaseTransforms size:', diagramElementBaseTransforms.value.size);
        debugLog('Linear: diagramElements count:', diagramElements.value.length);
        diagramElements.value.forEach((el, idx) => {
          const baseTransform = diagramElementBaseTransforms.value.get(el) || { x: 0, y: 0 };
          const foundInMap = diagramElementBaseTransforms.value.has(el);
          console.log(`[DEBUG] Linear Element ${idx} (${el.id}): foundInMap=${foundInMap}, baseTransform=`, baseTransform);
          const newX = baseTransform.x + shiftX;
          const newY = baseTransform.y + shiftY;
          el.setAttribute('transform', `translate(${newX}, ${newY})`);
          diagramElementOriginalTransforms.value.set(el, { x: newX, y: newY });
        });
        diagramOffset.x = 0;
        diagramOffset.y = 0;

        if (nowHorizontal && legendGroup) {
          const legendBounds = getTransformedBBox(legendGroup);
          const diagramBounds = getElementsBounds(diagramElements.value);
          const overlapPadding = 20;

          if (legendBounds && diagramBounds) {
            if (newPosition === 'top') {
              const overlap = legendBounds.y + legendBounds.height + overlapPadding - diagramBounds.y;
              if (overlap > 0) {
                applyDiagramShift(0, overlap);
                vbH += overlap;
                svg.setAttribute('viewBox', `${vbX} ${vbY} ${vbW} ${vbH}`);
                svg.setAttribute('data-horizontal-viewbox', `${vbX} ${vbY} ${vbW} ${vbH}`);
                linearBaseConfig.value.horizontalViewBox = { x: vbX, y: vbY, w: vbW, h: vbH };
              }
            } else if (newPosition === 'bottom') {
              const overlap = diagramBounds.y + diagramBounds.height + overlapPadding - legendBounds.y;
              if (overlap > 0) {
                vbH += overlap;
                svg.setAttribute('viewBox', `${vbX} ${vbY} ${vbW} ${vbH}`);
                svg.setAttribute('data-horizontal-viewbox', `${vbX} ${vbY} ${vbW} ${vbH}`);
                linearBaseConfig.value.horizontalViewBox = { x: vbX, y: vbY, w: vbW, h: vbH };
                const finalX = vbX + (vbW - legendWidth) / 2;
                const finalY = vbY + vbH - legendHeight - overlapPadding;
                legendGroup.setAttribute('transform', `translate(${finalX}, ${finalY})`);
                legendInitialTransform.value = { x: finalX, y: finalY };
                legendCurrentOffset.x = 0;
                legendCurrentOffset.y = 0;
              }
            }
          }
        }
      }
    }

    generatedLegendPosition.value = newPosition;

    const colors = currentColors.value;

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
    featureLegendGroups.forEach((featureLegendGroup) => {
      if (!featureLegendGroup) return;

      const entryGroups = featureLegendGroup.querySelectorAll('g[data-legend-key]');

      if (entryGroups.length > 0) {
        entryGroups.forEach((entryGroup) => {
          const legendKey = entryGroup.getAttribute('data-legend-key');
          if (!legendKey) return;
          if (legendColorOverrides[legendKey]) return;

          const colorKey = keyToColorKey[legendKey];
          const newColor = colorKey ? colors[colorKey] : colors[legendKey];
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
        const parseXY = (transform) => {
          if (!transform) return { x: 0, y: 0 };
          const match = transform.match(/translate\(\s*([\d.-]+)\s*,\s*([\d.-]+)\s*\)/);
          return match ? { x: parseFloat(match[1]), y: parseFloat(match[2]) } : { x: 0, y: 0 };
        };
        texts.forEach((textEl) => {
          const textContent = textEl.textContent?.trim();
          if (!textContent) return;
          if (legendColorOverrides[textContent]) return;
          const colorKey = keyToColorKey[textContent];
          const newColor = colorKey ? colors[colorKey] : colors[textContent];
          if (!newColor) return;
          const textPos = parseXY(textEl.getAttribute('transform'));
          let bestPath = null;
          let bestX = -Infinity;
          for (const path of allPaths) {
            const pathPos = parseXY(path.getAttribute('transform'));
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
          }
        });
      }
    });

    skipCaptureBaseConfig.value = true;
    skipPositionReapply.value = true;
    const idx = selectedResultIndex.value;
    if (idx >= 0 && results.value.length > idx) {
      const serializer = new XMLSerializer();
      const serialized = serializer.serializeToString(svg);

      if (!isLinear && diagramElements.value.length > 0) {
        const firstEl = diagramElements.value[0];
        const domTransform = firstEl.getAttribute('transform');
        const transformMatch = serialized.match(new RegExp(`id="${firstEl.id}"[^>]*transform="([^"]+)"`));
        const serializedTransform = transformMatch ? transformMatch[1] : 'NOT FOUND';
        console.log(`[DEBUG] Before save - DOM transform: ${domTransform}`);
        console.log(`[DEBUG] Before save - Serialized transform: ${serializedTransform}`);
      }

      results.value[idx] = { ...results.value[idx], content: serialized };
    }

    console.log(`Legend repositioned (${isLinear ? 'linear' : 'circular'}): ${oldPosition} -> ${newPosition}`);
  };

  return {
    repositionForLegendChange
  };
};
