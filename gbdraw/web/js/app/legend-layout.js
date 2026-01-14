import { estimateColorFactor, interpolateColor } from './color-utils.js';

export const createLegendLayout = ({ state, debugLog, legendActions, svgActions }) => {
  const {
    svgContent,
    svgContainer,
    mode,
    generatedLegendPosition,
    circularBaseConfig,
    linearBaseConfig,
    diagramElements,
    diagramElementIds,
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
    canvasPadding,
    originalSvgStroke,
    zoom,
    diagramDragging,
    diagramDragStart,
    form
  } = state;

  const { resetLegendPositionOnly, getAllFeatureLegendGroups, reflowSingleLegendLayout } = legendActions;
  const { ensureUniquePairwiseGradientIds } = svgActions;

  const parseTransform = (transformStr) => {
    if (!transformStr) return { x: 0, y: 0 };
    const match = transformStr.match(/translate\(\s*([-\d.]+)\s*,?\s*([-\d.]+)?\s*\)/);
    if (match) {
      return { x: parseFloat(match[1]) || 0, y: parseFloat(match[2]) || 0 };
    }
    return { x: 0, y: 0 };
  };

  const getTransformedBBox = (el) => {
    if (!el) return null;
    const bbox = el.getBBox();
    const offset = parseTransform(el.getAttribute('transform'));
    return { x: bbox.x + offset.x, y: bbox.y + offset.y, width: bbox.width, height: bbox.height };
  };

  const getElementsBounds = (elements) => {
    let minX = Infinity;
    let minY = Infinity;
    let maxX = -Infinity;
    let maxY = -Infinity;
    elements.forEach((el) => {
      if (!el) return;
      const bounds = getTransformedBBox(el);
      if (!bounds) return;
      minX = Math.min(minX, bounds.x);
      minY = Math.min(minY, bounds.y);
      maxX = Math.max(maxX, bounds.x + bounds.width);
      maxY = Math.max(maxY, bounds.y + bounds.height);
    });
    if (minX === Infinity) return null;
    return { x: minX, y: minY, width: maxX - minX, height: maxY - minY };
  };

  const applyDiagramShift = (deltaX, deltaY) => {
    if (diagramElements.value.length === 0) return;
    diagramElements.value.forEach((el) => {
      const current = parseTransform(el.getAttribute('transform'));
      const newX = current.x + deltaX;
      const newY = current.y + deltaY;
      el.setAttribute('transform', `translate(${newX}, ${newY})`);
      diagramElementOriginalTransforms.value.set(el, { x: newX, y: newY });
    });
    diagramOffset.x = 0;
    diagramOffset.y = 0;
  };

  const startDiagramDrag = (e) => {
    if (
      e.target.closest('#legend') ||
      e.target.closest('path[data-feature-id]') ||
      e.target.closest('rect[data-feature-id]')
    ) {
      return;
    }
    if (!svgContainer.value) return;
    const svg = svgContainer.value.querySelector('svg');
    if (!svg) return;

    const clickedGroup = e.target.closest('g[id]');
    if (!clickedGroup) return;

    const clickedId = clickedGroup.id;
    if (
      clickedId === 'legend' ||
      clickedId === 'feature_legend' ||
      clickedId === 'pairwise_legend' ||
      clickedId === 'horizontal_legend' ||
      clickedId === 'vertical_legend'
    ) {
      return;
    }

    e.preventDefault();
    diagramDragging.value = true;
    diagramDragStart.x = e.clientX;
    diagramDragStart.y = e.clientY;

    diagramElements.value.forEach((el) => {
      el.style.opacity = '0.8';
    });

    document.addEventListener('mousemove', onDiagramDrag);
    document.addEventListener('mouseup', endDiagramDrag);
  };

  const onDiagramDrag = (e) => {
    if (!diagramDragging.value) return;

    const deltaX = (e.clientX - diagramDragStart.x) / zoom.value;
    const deltaY = (e.clientY - diagramDragStart.y) / zoom.value;

    diagramElements.value.forEach((el) => {
      const original = diagramElementOriginalTransforms.value.get(el) || { x: 0, y: 0 };
      const newX = original.x + diagramOffset.x + deltaX;
      const newY = original.y + diagramOffset.y + deltaY;
      el.setAttribute('transform', `translate(${newX}, ${newY})`);
    });
  };

  const endDiagramDrag = (e) => {
    if (!diagramDragging.value) return;

    const deltaX = (e.clientX - diagramDragStart.x) / zoom.value;
    const deltaY = (e.clientY - diagramDragStart.y) / zoom.value;
    diagramOffset.x += deltaX;
    diagramOffset.y += deltaY;

    diagramDragging.value = false;
    document.removeEventListener('mousemove', onDiagramDrag);
    document.removeEventListener('mouseup', endDiagramDrag);

    diagramElements.value.forEach((el) => {
      el.style.opacity = '1';
    });
  };

  const resetDiagramPosition = () => {
    diagramOffset.x = 0;
    diagramOffset.y = 0;

    console.log('[DEBUG] resetDiagramPosition called');
    console.log('[DEBUG] diagramElements count:', diagramElements.value.length);
    console.log('[DEBUG] originalTransforms size:', diagramElementOriginalTransforms.value.size);
    console.log('[DEBUG] All originalTransforms entries:');
    diagramElementOriginalTransforms.value.forEach((transform, mapEl) => {
      console.log(`[DEBUG]   Map entry: ${mapEl.id} -> (${transform.x}, ${transform.y})`);
    });

    diagramElements.value.forEach((el, idx) => {
      const original = diagramElementOriginalTransforms.value.get(el);
      const currentTransform = el.getAttribute('transform');
      let foundInMap = false;
      diagramElementOriginalTransforms.value.forEach((_, mapEl) => {
        if (mapEl === el) foundInMap = true;
      });
      console.log(
        `[DEBUG] Reset element ${idx} (${el.id}): current="${currentTransform}", foundInMap=${foundInMap}, original=`,
        original
      );
      if (original && (original.x !== 0 || original.y !== 0)) {
        el.setAttribute('transform', `translate(${original.x}, ${original.y})`);
      } else {
        el.removeAttribute('transform');
      }
    });
  };

  const resetAllPositions = () => {
    resetDiagramPosition();
    resetLegendPositionOnly();
  };

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
          `[DEBUG] ⚠️ Legend size CHANGED: ${origLegendWidth.toFixed(1)}x${origLegendHeight.toFixed(
            1
          )} -> ${legendWidth.toFixed(1)}x${legendHeight.toFixed(1)}`
        );
        console.log(
          `[DEBUG] ⚠️ Shift calculation affected: old=${(origLegendWidth * 0.55).toFixed(1)}, new=${(
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

          if (!hasDualLegends) {
            const layout = nowHorizontal ? 'horizontal' : 'vertical';
            const reflowResult = reflowSingleLegendLayout(svg, layout, vbW);
            if (reflowResult) {
              legendWidth = reflowResult.legendWidth || legendWidth;
              legendHeight = reflowResult.legendHeight || legendHeight;
            }
          }

          const padding = 20;

          if (nowHorizontal) {
            linearBaseConfig.value.horizontalLegendWidth = legendWidth;
            linearBaseConfig.value.horizontalLegendHeight = legendHeight;

            if (legendWidth + padding * 2 > vbW) {
              vbW = legendWidth + padding * 2;
            }

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

          let finalX;
          let finalY;

          switch (newPosition) {
            case 'top':
              finalX = (vbW - legendWidth) / 2;
              finalY = padding;
              break;
            case 'bottom':
              finalX = (vbW - legendWidth) / 2;
              finalY = vbH - legendHeight - padding;
              break;
            case 'left':
              finalX = padding;
              finalY = (vbH - legendHeight) / 2;
              break;
            case 'right':
            default:
              finalX = vbW - legendWidth - padding;
              finalY = (vbH - legendHeight) / 2;
              break;
          }

          legendGroup.setAttribute('transform', `translate(${finalX}, ${finalY})`);
          legendInitialTransform.value = { x: finalX, y: finalY };
          legendCurrentOffset.x = 0;
          legendCurrentOffset.y = 0;
        }
      }

      const hasStoredViewBox =
        svg.hasAttribute('data-horizontal-viewbox') && svg.hasAttribute('data-vertical-viewbox');

      if (hasStoredViewBox && diagramElements.value.length > 0) {
        const hVbAttr = svg.getAttribute('data-horizontal-viewbox');
        const vVbAttr = svg.getAttribute('data-vertical-viewbox');
        const hVb = hVbAttr ? hVbAttr.split(/\s+/).map(parseFloat) : null;
        const vVb = vVbAttr ? vVbAttr.split(/\s+/).map(parseFloat) : null;

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
                const finalX = (vbW - legendWidth) / 2;
                const finalY = vbH - legendHeight - overlapPadding;
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

  let setupDiagramDragCallCount = 0;

  const setupDiagramDrag = (preserveOffset = false) => {
    setupDiagramDragCallCount++;
    if (!svgContainer.value) return;
    const svg = svgContainer.value.querySelector('svg');
    if (!svg) return;

    if (!preserveOffset) {
      diagramOffset.x = 0;
      diagramOffset.y = 0;
    }

    const knownIds = ['tick', 'labels', 'Axis', 'gc_content', 'skew', 'gc_skew', 'length_bar'];
    const foundElements = [];
    const foundIds = [];

    knownIds.forEach((id) => {
      const el = svg.getElementById(id);
      if (el) {
        foundElements.push(el);
        foundIds.push(id);
      }
    });

    const allGroups = svg.querySelectorAll('g[id]');
    allGroups.forEach((group) => {
      const id = group.id;
      if (!id) return;

      if (id === 'legend' || id === 'feature_legend' || id === 'pairwise_legend') return;

      if (foundElements.includes(group)) return;

      const isAccession = id.match(/^[A-Z]{2}_?\d+/) || id.match(/^[A-Z]+\d+\.\d+$/);
      const isKnown = knownIds.includes(id);
      const isDynamic =
        id.startsWith('record_') ||
        id.startsWith('definition_') ||
        id.startsWith('seq_') ||
        id.startsWith('track_') ||
        id.startsWith('match_') ||
        id.startsWith('comparison');

      if (isAccession || isKnown || isDynamic) {
        foundElements.push(group);
        if (!foundIds.includes(id)) {
          foundIds.push(id);
        }
      }
    });

    diagramElements.value = foundElements;
    diagramElementIds.value = foundIds;

    const originalTransforms = new Map();
    console.log(`[DEBUG] ========== setupDiagramDrag CALL #${setupDiagramDragCallCount} ==========`);
    console.log(
      `[DEBUG] setupDiagramDrag: preserveOffset=${preserveOffset}, offset=(${diagramOffset.x}, ${diagramOffset.y}), generatedLegendPosition=${generatedLegendPosition.value}`
    );
    foundElements.forEach((el, idx) => {
      const transform = parseTransform(el.getAttribute('transform'));
      console.log(`[DEBUG] setupDiagramDrag element ${idx} (${el.id}): DOM transform=(${transform.x}, ${transform.y})`);
      if (preserveOffset && (diagramOffset.x !== 0 || diagramOffset.y !== 0)) {
        const adjusted = {
          x: transform.x - diagramOffset.x,
          y: transform.y - diagramOffset.y
        };
        console.log(`[DEBUG] setupDiagramDrag element ${idx}: adjusted to (${adjusted.x}, ${adjusted.y})`);
        originalTransforms.set(el, adjusted);
      } else {
        originalTransforms.set(el, transform);
      }
    });
    diagramElementOriginalTransforms.value = originalTransforms;
    console.log(`[DEBUG] setupDiagramDrag FINISHED: Set ${originalTransforms.size} original transforms`);
    originalTransforms.forEach((transform, el) => {
      console.log(`[DEBUG]   -> ${el.id}: (${transform.x}, ${transform.y})`);
    });

    foundElements.forEach((el) => {
      el.style.cursor = 'grab';
    });

    svg.removeEventListener('mousedown', startDiagramDrag);
    svg.addEventListener('mousedown', startDiagramDrag);
  };

  const applyCanvasPadding = () => {
    if (!svgContainer.value) return;
    const svg = svgContainer.value.querySelector('svg');
    if (!svg) return;

    let viewBox = svg.getAttribute('viewBox');
    if (!viewBox) {
      const width = parseFloat(svg.getAttribute('width')) || 800;
      const height = parseFloat(svg.getAttribute('height')) || 600;
      viewBox = `0 0 ${width} ${height}`;
    }

    const parts = viewBox.split(/\s+/).map(parseFloat);
    if (parts.length !== 4) return;

    if (!svg.dataset.originalViewBox) {
      svg.dataset.originalViewBox = viewBox;
    }

    const [origX, origY, origW, origH] = svg.dataset.originalViewBox.split(/\s+/).map(parseFloat);

    const newX = origX - canvasPadding.left;
    const newY = origY - canvasPadding.top;
    const newW = origW + canvasPadding.left + canvasPadding.right;
    const newH = origH + canvasPadding.top + canvasPadding.bottom;

    svg.setAttribute('viewBox', `${newX} ${newY} ${newW} ${newH}`);

    const origWidth = parseFloat(svg.dataset.originalWidth || svg.getAttribute('width')) || origW;
    const origHeight = parseFloat(svg.dataset.originalHeight || svg.getAttribute('height')) || origH;
    if (!svg.dataset.originalWidth) {
      svg.dataset.originalWidth = origWidth;
      svg.dataset.originalHeight = origHeight;
    }

    const scaleX = newW / origW;
    const scaleY = newH / origH;
    svg.setAttribute('width', origWidth * scaleX);
    svg.setAttribute('height', origHeight * scaleY);
  };

  const resetCanvasPadding = () => {
    canvasPadding.top = 0;
    canvasPadding.right = 0;
    canvasPadding.bottom = 0;
    canvasPadding.left = 0;

    if (!svgContainer.value) return;
    const svg = svgContainer.value.querySelector('svg');
    if (!svg) return;

    if (svg.dataset.originalViewBox) {
      svg.setAttribute('viewBox', svg.dataset.originalViewBox);
    }
    if (svg.dataset.originalWidth) {
      svg.setAttribute('width', svg.dataset.originalWidth);
      svg.setAttribute('height', svg.dataset.originalHeight);
    }
  };

  const captureOriginalStroke = () => {
    if (!svgContainer.value) return;
    const svg = svgContainer.value.querySelector('svg');
    if (!svg) return;

    const firstFeaturePath = svg.querySelector('path[id^="f"]');
    if (firstFeaturePath) {
      const strokeColor = firstFeaturePath.getAttribute('stroke');
      const strokeWidthAttr = firstFeaturePath.getAttribute('stroke-width');
      const strokeWidth = strokeWidthAttr !== null ? parseFloat(strokeWidthAttr) : null;

      originalSvgStroke.value = { color: strokeColor, width: strokeWidth };
      console.log(`Captured original stroke: color=${strokeColor}, width=${strokeWidth}`);
    }
  };

  const captureBaseConfig = () => {
    if (!svgContainer.value) return;
    const svg = svgContainer.value.querySelector('svg');
    if (!svg) return;

    const viewBox = svg.getAttribute('viewBox');
    if (!viewBox) return;
    const parts = viewBox.split(/\s+/).map(parseFloat);
    if (parts.length !== 4) return;
    const [vbX, vbY, vbW, vbH] = parts;

    const legendGroup = svg.getElementById('legend');
    let legendWidth = 120;
    let legendHeight = 150;
    if (legendGroup) {
      const bbox = legendGroup.getBBox();
      legendWidth = bbox.width || legendWidth;
      legendHeight = bbox.height || legendHeight;
    }

    const isLinear = mode.value === 'linear';

    if (isLinear) {
      const verticalVb = svg.getAttribute('data-vertical-viewbox');
      const horizontalVb = svg.getAttribute('data-horizontal-viewbox');

      if (verticalVb && horizontalVb) {
        const vParts = verticalVb.split(/\s+/).map(parseFloat);
        const hParts = horizontalVb.split(/\s+/).map(parseFloat);
        linearBaseConfig.value.verticalViewBox = { x: vParts[0], y: vParts[1], w: vParts[2], h: vParts[3] };
        linearBaseConfig.value.horizontalViewBox = { x: hParts[0], y: hParts[1], w: hParts[2], h: hParts[3] };
      } else {
        let vLegendW = legendWidth;
        let vLegendH = legendHeight;
        let hLegendW = legendWidth;
        let hLegendH = legendHeight;

        if (legendGroup) {
          const horizontalLegend = legendGroup.querySelector('#legend_horizontal');
          const verticalLegend = legendGroup.querySelector('#legend_vertical');
          if (horizontalLegend && verticalLegend) {
            const hDisplay = horizontalLegend.getAttribute('display');
            const vDisplay = verticalLegend.getAttribute('display');
            horizontalLegend.removeAttribute('display');
            verticalLegend.removeAttribute('display');

            const hBbox = horizontalLegend.getBBox();
            const vBbox = verticalLegend.getBBox();
            hLegendW = hBbox.width || legendWidth;
            hLegendH = hBbox.height || legendHeight;
            vLegendW = vBbox.width || legendWidth;
            vLegendH = vBbox.height || legendHeight;

            if (hDisplay) horizontalLegend.setAttribute('display', hDisplay);
            else horizontalLegend.removeAttribute('display');
            if (vDisplay) verticalLegend.setAttribute('display', vDisplay);
            else verticalLegend.removeAttribute('display');
          }
        }

        const genPos = form.legend;
        const isGeneratedVertical = genPos === 'left' || genPos === 'right';

        if (isGeneratedVertical) {
          linearBaseConfig.value.verticalViewBox = { x: vbX, y: vbY, w: vbW, h: vbH };
          linearBaseConfig.value.horizontalViewBox = {
            x: vbX,
            y: vbY,
            w: vbW - vLegendW,
            h: vbH + hLegendH
          };
        } else {
          linearBaseConfig.value.horizontalViewBox = { x: vbX, y: vbY, w: vbW, h: vbH };
          linearBaseConfig.value.verticalViewBox = {
            x: vbX,
            y: vbY,
            w: vbW + vLegendW,
            h: vbH - hLegendH
          };
        }
      }

      if (legendGroup) {
        const horizontalLegend = legendGroup.querySelector('#legend_horizontal');
        const verticalLegend = legendGroup.querySelector('#legend_vertical');

        if (horizontalLegend && verticalLegend) {
          const hDisplay = horizontalLegend.getAttribute('display');
          const vDisplay = verticalLegend.getAttribute('display');
          horizontalLegend.removeAttribute('display');
          verticalLegend.removeAttribute('display');

          const hBbox = horizontalLegend.getBBox();
          const vBbox = verticalLegend.getBBox();

          linearBaseConfig.value.horizontalLegendWidth = hBbox.width || legendWidth;
          linearBaseConfig.value.horizontalLegendHeight = hBbox.height || legendHeight;
          linearBaseConfig.value.verticalLegendWidth = vBbox.width || legendWidth;
          linearBaseConfig.value.verticalLegendHeight = vBbox.height || legendHeight;

          if (hDisplay) horizontalLegend.setAttribute('display', hDisplay);
          if (vDisplay) verticalLegend.setAttribute('display', vDisplay);
        } else {
          linearBaseConfig.value.verticalLegendWidth = legendWidth;
          linearBaseConfig.value.verticalLegendHeight = legendHeight;
          linearBaseConfig.value.horizontalLegendWidth = legendWidth;
          linearBaseConfig.value.horizontalLegendHeight = legendHeight;
        }
      }

      linearBaseConfig.value.generatedPosition = form.legend;
      linearBaseConfig.value.diagramBaseTransforms = new Map(diagramElementOriginalTransforms.value);
    } else {
      const legendPos = form.legend;
      let baseVbW = vbW;

      if (legendPos === 'left' || legendPos === 'right') {
        baseVbW = vbW - legendWidth * 1.1;
      }

      const diagramCenterX = baseVbW / 2;
      const diagramCenterY = vbH / 2;

      circularBaseConfig.value = {
        viewBoxWidth: baseVbW,
        viewBoxHeight: vbH,
        generatedViewBoxWidth: vbW,
        generatedViewBoxHeight: vbH,
        diagramCenterX: diagramCenterX,
        diagramCenterY: diagramCenterY,
        legendWidth: legendWidth,
        legendHeight: legendHeight,
        generatedPosition: legendPos
      };
    }

    diagramElementBaseTransforms.value = new Map(diagramElementOriginalTransforms.value);
  };

  return {
    applyCanvasPadding,
    resetCanvasPadding,
    captureBaseConfig,
    captureOriginalStroke,
    repositionForLegendChange,
    setupDiagramDrag,
    resetAllPositions
  };
};
