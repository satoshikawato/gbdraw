import { getFeatureCaption, ruleMatchesFeature } from './feature-utils.js';

export const createLegendManager = ({ state, getPyodide, debugLog }) => {
  const {
    pyodideReady,
    results,
    selectedResultIndex,
    svgContainer,
    adv,
    mode,
    legendEntries,
    deletedLegendEntries,
    originalLegendOrder,
    originalLegendColors,
    newLegendCaption,
    newLegendColor,
    legendStrokeOverrides,
    legendColorOverrides,
    originalSvgStroke,
    legendDragging,
    legendDragStart,
    legendOriginalTransform,
    legendInitialTransform,
    legendCurrentOffset,
    zoom,
    extractedFeatures,
    manualSpecificRules,
    skipCaptureBaseConfig
  } = state;

  const addLegendEntry = async (caption, color) => {
    console.log(`addLegendEntry called with caption="${caption}", color="${color}"`);
    if (!svgContainer.value || !pyodideReady.value) {
      console.log(
        `addLegendEntry early return: svgContainer=${!!svgContainer.value}, pyodideReady=${pyodideReady.value}`
      );
      return false;
    }

    const svg = svgContainer.value.querySelector('svg');
    if (!svg) return false;

    const legendGroup = svg.getElementById('legend');
    if (!legendGroup) {
      console.log('No legend group found');
      return false;
    }

    const allTargetGroups = getAllFeatureLegendGroups(svg);
    if (allTargetGroups.length === 0) {
      console.log('No feature legend groups found');
      return false;
    }

    const targetGroup = allTargetGroups[0];

    const existingTexts = targetGroup.querySelectorAll('text');
    let finalCaption = caption.trim();

    for (const t of existingTexts) {
      const existingCaption = t.textContent?.trim();
      if (existingCaption === finalCaption) {
        const transform = t.getAttribute('transform');
        if (transform) {
          const match = transform.match(/translate\([^,]+,\s*([\d.]+)\)/);
          if (match) {
            const y = match[1];
            const rects = targetGroup.querySelectorAll('path');
            for (const r of rects) {
              const rt = r.getAttribute('transform');
              if (rt && rt.includes(`, ${y})`)) {
                const existingColor = r.getAttribute('fill');
                if (existingColor === color) {
                  console.log(`addLegendEntry: same caption and color already exist, returning: ${finalCaption}`);
                  return finalCaption;
                }
                break;
              }
            }
          }
        }
      }
    }

    const baseCaption = finalCaption.replace(/\s*\(\d+\)$/, '');
    let counter = 1;
    const existingCaptions = new Set();
    existingTexts.forEach((t) => existingCaptions.add(t.textContent?.trim()));

    while (existingCaptions.has(finalCaption)) {
      finalCaption = `${baseCaption} (${counter})`;
      counter++;
    }

    caption = finalCaption;

    let rectSize = 14;
    const firstColorRect = targetGroup.querySelector('path[fill]:not([fill="none"]):not([fill^="url("])');
    if (firstColorRect) {
      const d = firstColorRect.getAttribute('d');
      if (d) {
        const lMatch = d.match(/L\s+([\d.]+),/);
        if (lMatch) {
          rectSize = parseFloat(lMatch[1]);
        }
      }
    }

    let fontSize = 14;
    let fontFamily = 'Arial';
    const firstText = targetGroup.querySelector('text');
    if (firstText) {
      const fs = firstText.getAttribute('font-size');
      if (fs) fontSize = parseFloat(fs);
      const ff = firstText.getAttribute('font-family');
      if (ff) fontFamily = ff;
    }

    let strokeColor = 'black';
    let strokeWidth = 0.5;
    if (adv.block_stroke_color) {
      strokeColor = adv.block_stroke_color;
    }
    if (adv.block_stroke_width !== null && adv.block_stroke_width !== undefined) {
      strokeWidth = adv.block_stroke_width;
    }
    if (!adv.block_stroke_color && firstColorRect) {
      const existingStroke = firstColorRect.getAttribute('stroke');
      if (existingStroke && existingStroke !== 'none') {
        strokeColor = existingStroke;
      }
    }
    if ((adv.block_stroke_width === null || adv.block_stroke_width === undefined) && firstColorRect) {
      const existingStrokeWidth = firstColorRect.getAttribute('stroke-width');
      if (existingStrokeWidth) {
        strokeWidth = parseFloat(existingStrokeWidth);
      }
    }

    const lineMargin = (24 / 14) * rectSize;
    const xMargin = (22 / 14) * rectSize;

    let maxY = -lineMargin;
    const textElements = targetGroup.querySelectorAll('text');
    textElements.forEach((el) => {
      const transform = el.getAttribute('transform');
      if (transform) {
        const match = transform.match(/translate\([^,]+,\s*([\d.-]+)\)/);
        if (match) {
          const y = parseFloat(match[1]);
          if (y > maxY) maxY = y;
        }
      }
    });

    if (textElements.length === 0) {
      const colorRects = targetGroup.querySelectorAll('path');
      colorRects.forEach((el) => {
        const fill = el.getAttribute('fill');
        if (fill && fill !== 'none' && !fill.startsWith('url(')) {
          const transform = el.getAttribute('transform');
          if (transform) {
            const match = transform.match(/translate\([^,]+,\s*([\d.-]+)\)/);
            if (match) {
              const y = parseFloat(match[1]);
              if (y > maxY) maxY = y;
            }
          }
        }
      });
    }

    const newY = maxY + lineMargin;

    try {
      const escapedCaption = caption.replace(/\\/g, '\\\\').replace(/"/g, '\\"');
      const escapedFontFamily = fontFamily.replace(/\\/g, '\\\\').replace(/"/g, '\\"');
      const parser = new DOMParser();

      const parseTransformXY = (transform) => {
        if (!transform) return { x: 0, y: 0 };
        const match = transform.match(/translate\(\s*([\d.-]+)\s*,\s*([\d.-]+)\s*\)/);
        return match ? { x: parseFloat(match[1]), y: parseFloat(match[2]) } : { x: 0, y: 0 };
      };

      let entryWidth = rectSize + xMargin + 100;
      try {
        const pyodide = getPyodide();
        if (!pyodide) throw new Error('Pyodide not ready');
        const widthResult = pyodide.runPython(`
import json
from gbdraw.core.text import calculate_bbox_dimensions
width, _ = calculate_bbox_dimensions("${escapedCaption}", "${escapedFontFamily}", ${fontSize}, 72)
json.dumps({"width": width})
`);
        const widthData = JSON.parse(widthResult);
        entryWidth = rectSize + xMargin + widthData.width + xMargin;
      } catch (e) {
        console.warn('Could not calculate text width, using estimate');
      }

      const viewBox = svg.getAttribute('viewBox');
      let canvasWidth = 800;
      if (viewBox) {
        const parts = viewBox.split(/\s+/);
        if (parts.length >= 4) canvasWidth = parseFloat(parts[2]);
      }

      for (const group of allTargetGroups) {
        const parentId = group.parentElement?.id || '';
        const isHorizontalGroup = parentId === 'legend_horizontal';

        let newX = 0,
          newY = 0;

        if (isHorizontalGroup) {
          let featureLegendMaxWidth = canvasWidth;
          const horizontalLegend = legendGroup.querySelector('#legend_horizontal');
          if (horizontalLegend) {
            const hPairwiseLegend = horizontalLegend.querySelector('#pairwise_legend');
            if (hPairwiseLegend) {
              const pairwiseTransform = hPairwiseLegend.getAttribute('transform');
              if (pairwiseTransform) {
                const match = pairwiseTransform.match(/translate\(\s*([\d.-]+)\s*,/);
                if (match) {
                  featureLegendMaxWidth = parseFloat(match[1]) - xMargin;
                }
              }
            }
          }

          let maxY = rectSize / 2;
          let maxXOnMaxY = 0;
          let lastEntryRightEdge = 0;

          const groupTextElements = group.querySelectorAll('text');
          groupTextElements.forEach((el) => {
            const pos = parseTransformXY(el.getAttribute('transform'));
            if (pos.y > maxY) {
              maxY = pos.y;
              maxXOnMaxY = pos.x;
              const textBBox = el.getBBox();
              lastEntryRightEdge = pos.x + textBBox.width + xMargin;
            } else if (Math.abs(pos.y - maxY) < 1) {
              if (pos.x > maxXOnMaxY) {
                maxXOnMaxY = pos.x;
                const textBBox = el.getBBox();
                lastEntryRightEdge = pos.x + textBBox.width + xMargin;
              }
            }
          });

          if (groupTextElements.length === 0) {
            const colorRects = group.querySelectorAll('path');
            colorRects.forEach((el) => {
              const fill = el.getAttribute('fill');
              if (fill && fill !== 'none' && !fill.startsWith('url(')) {
                const pos = parseTransformXY(el.getAttribute('transform'));
                if (pos.y > maxY) {
                  maxY = pos.y;
                  maxXOnMaxY = pos.x;
                  lastEntryRightEdge = pos.x + rectSize + xMargin;
                } else if (Math.abs(pos.y - maxY) < 1) {
                  if (pos.x > maxXOnMaxY) {
                    maxXOnMaxY = pos.x + rectSize;
                    lastEntryRightEdge = pos.x + rectSize + xMargin;
                  }
                }
              }
            });
          }

          let nextX = groupTextElements.length > 0 ? lastEntryRightEdge : 0;

          if (nextX + entryWidth > featureLegendMaxWidth && nextX > 0) {
            newX = 0;
            newY = maxY + lineMargin;
          } else {
            newX = nextX;
            newY = maxY;
          }
        } else {
          let groupMaxY = -lineMargin;
          const groupTextElements = group.querySelectorAll('text');
          groupTextElements.forEach((el) => {
            const pos = parseTransformXY(el.getAttribute('transform'));
            if (pos.y > groupMaxY) groupMaxY = pos.y;
          });

          if (groupTextElements.length === 0) {
            const colorRects = group.querySelectorAll('path');
            colorRects.forEach((el) => {
              const fill = el.getAttribute('fill');
              if (fill && fill !== 'none' && !fill.startsWith('url(')) {
                const pos = parseTransformXY(el.getAttribute('transform'));
                if (pos.y > groupMaxY) groupMaxY = pos.y;
              }
            });
          }

          newX = 0;
          newY = groupMaxY + lineMargin;
        }

        const pyodide = getPyodide();
        if (!pyodide) throw new Error('Pyodide not ready');
        const resultJson = pyodide.runPython(
          `generate_legend_entry_svg("${escapedCaption}", "${color}", ${newY}, ${rectSize}, ${fontSize}, "${escapedFontFamily}", ${newX}, "${strokeColor}", ${strokeWidth})`
        );
        const result = JSON.parse(resultJson);

        const entryGroup = document.createElementNS('http://www.w3.org/2000/svg', 'g');
        entryGroup.setAttribute('data-legend-key', caption);

        const rectDoc = parser.parseFromString(
          `<svg xmlns="http://www.w3.org/2000/svg">${result.rect}</svg>`,
          'image/svg+xml'
        );
        const rectEl = rectDoc.querySelector('path');
        if (rectEl) {
          entryGroup.appendChild(document.importNode(rectEl, true));
        }

        const textDoc = parser.parseFromString(
          `<svg xmlns="http://www.w3.org/2000/svg">${result.text}</svg>`,
          'image/svg+xml'
        );
        const textEl = textDoc.querySelector('text');
        if (textEl) {
          entryGroup.appendChild(document.importNode(textEl, true));
        }

        group.appendChild(entryGroup);
      }

      const hasDualLegends =
        !!legendGroup.querySelector('#legend_horizontal') && !!legendGroup.querySelector('#legend_vertical');
      if (hasDualLegends) {
        reflowDualLegendLayout(svg);
      } else {
        updatePairwiseLegendPositions(svg);
      }

      skipCaptureBaseConfig.value = true;
      const idx = selectedResultIndex.value;
      if (idx >= 0 && results.value.length > idx) {
        const serializer = new XMLSerializer();
        results.value[idx] = { ...results.value[idx], content: serializer.serializeToString(svg) };
      }

      return caption;
    } catch (e) {
      console.error('Failed to add legend entry:', e);
      return false;
    }
  };

  const getAllFeatureLegendGroups = (svg) => {
    if (!svg) return [];

    const legendGroup = svg.getElementById('legend');
    if (!legendGroup) return [];

    const horizontalLegend = legendGroup.querySelector('#legend_horizontal');
    const verticalLegend = legendGroup.querySelector('#legend_vertical');

    const groups = [];

    if (horizontalLegend && verticalLegend) {
      const hFeatureLegend = horizontalLegend.querySelector('#feature_legend_h');
      const vFeatureLegend = verticalLegend.querySelector('#feature_legend_v');
      if (hFeatureLegend) groups.push(hFeatureLegend);
      if (vFeatureLegend) groups.push(vFeatureLegend);
      return groups;
    }

    const featureLegendGroup = legendGroup.querySelector('#feature_legend');
    return featureLegendGroup ? [featureLegendGroup] : [legendGroup];
  };

  const getVisibleFeatureLegendGroup = (svg) => {
    if (!svg) return null;

    const legendGroup = svg.getElementById('legend');
    if (!legendGroup) return null;

    const horizontalLegend = legendGroup.querySelector('#legend_horizontal');
    const verticalLegend = legendGroup.querySelector('#legend_vertical');

    if (horizontalLegend && verticalLegend) {
      const hVisible =
        horizontalLegend.style.display !== 'none' &&
        horizontalLegend.getAttribute('display') !== 'none';
      const targetLegend = hVisible ? horizontalLegend : verticalLegend;
      const featureGroup = hVisible
        ? targetLegend.querySelector('#feature_legend_h')
        : targetLegend.querySelector('#feature_legend_v');
      return featureGroup || targetLegend;
    }

    const featureLegendGroup = legendGroup.querySelector('#feature_legend');
    return featureLegendGroup || legendGroup;
  };

  const isCurrentLegendHorizontal = (svg) => {
    const legendGroup = svg.getElementById('legend');
    if (!legendGroup) return false;

    const horizontalLegend = legendGroup.querySelector('#legend_horizontal');
    const verticalLegend = legendGroup.querySelector('#legend_vertical');

    if (horizontalLegend && verticalLegend) {
      const hVisible =
        horizontalLegend.style.display !== 'none' &&
        horizontalLegend.getAttribute('display') !== 'none';
      return hVisible;
    }

    return false;
  };

  const updateLegendEntryColorByCaption = (caption, color) => {
    if (!svgContainer.value) return false;

    const svg = svgContainer.value.querySelector('svg');
    if (!svg) return false;

    const targetGroups = getAllFeatureLegendGroups(svg);
    if (targetGroups.length === 0) return false;

    let updated = false;

    for (const targetGroup of targetGroups) {
      const entryGroup = targetGroup.querySelector(`g[data-legend-key="${CSS.escape(caption)}"]`);
      if (entryGroup) {
        const paths = entryGroup.querySelectorAll('path');
        for (const path of paths) {
          const fill = path.getAttribute('fill');
          if (fill && fill !== 'none' && !fill.startsWith('url(')) {
            path.setAttribute('fill', color);
            updated = true;
            break;
          }
        }
      }
    }

    if (updated) {
      skipCaptureBaseConfig.value = true;
      const idx = selectedResultIndex.value;
      if (idx >= 0 && results.value.length > idx) {
        const serializer = new XMLSerializer();
        results.value[idx] = { ...results.value[idx], content: serializer.serializeToString(svg) };
      }
      console.log(`Updated legend entry color: "${caption}" to ${color}`);
    }

    return updated;
  };

  const legendEntryExists = (caption) => {
    if (!svgContainer.value) return false;

    const svg = svgContainer.value.querySelector('svg');
    if (!svg) return false;

    const targetGroups = getAllFeatureLegendGroups(svg);
    if (targetGroups.length === 0) return false;

    const targetGroup = targetGroups[0];

    const entryGroup = targetGroup.querySelector(`g[data-legend-key="${CSS.escape(caption)}"]`);
    return entryGroup !== null;
  };

  const parseTransformXY = (transform) => {
    if (!transform) return { x: 0, y: 0 };
    const match = transform.match(/translate\(\s*([\d.-]+)\s*,\s*([\d.-]+)\s*\)/);
    return match ? { x: parseFloat(match[1]), y: parseFloat(match[2]) } : { x: 0, y: 0 };
  };

  const removeLegendEntry = (caption) => {
    if (!svgContainer.value) return false;

    const svg = svgContainer.value.querySelector('svg');
    if (!svg) return false;

    const targetGroups = getAllFeatureLegendGroups(svg);
    if (targetGroups.length === 0) return false;

    let removed = false;

    for (const targetGroup of targetGroups) {
      const entryGroup = targetGroup.querySelector(`g[data-legend-key="${CSS.escape(caption)}"]`);
      if (entryGroup) {
        entryGroup.remove();
        removed = true;
      }
    }

    if (removed) {
      compactLegendEntries(svg);

      skipCaptureBaseConfig.value = true;
      const idx = selectedResultIndex.value;
      if (idx >= 0 && results.value.length > idx) {
        const serializer = new XMLSerializer();
        results.value[idx] = { ...results.value[idx], content: serializer.serializeToString(svg) };
      }
      console.log(`Removed legend entry: "${caption}"`);
    }

    return removed;
  };

  const expandCanvasForVerticalLegend = (svg) => {
    if (mode.value !== 'linear') return;
    if (isCurrentLegendHorizontal(svg)) return;

    const legendGroup = svg.getElementById('legend');
    if (!legendGroup) return;

    const verticalLegend = legendGroup.querySelector('#legend_vertical');
    if (!verticalLegend || verticalLegend.getAttribute('display') === 'none') return;

    let legendGroupY = 0;
    const legendGroupTransform = legendGroup.getAttribute('transform');
    if (legendGroupTransform) {
      const match = legendGroupTransform.match(/translate\(\s*([\d.-]+)\s*,\s*([\d.-]+)\s*\)/);
      if (match) legendGroupY = parseFloat(match[2]);
    }

    const featureLegend = verticalLegend.querySelector('#feature_legend_v');
    const pairwiseLegend = verticalLegend.querySelector('#pairwise_legend');

    let legendContentHeight = 0;
    if (featureLegend) {
      const featureBBox = featureLegend.getBBox();
      legendContentHeight = featureBBox.y + featureBBox.height;
    }

    if (pairwiseLegend) {
      const pairwiseTransform = pairwiseLegend.getAttribute('transform');
      let pairwiseY = 0;
      if (pairwiseTransform) {
        const match = pairwiseTransform.match(/translate\(\s*([\d.-]+)\s*,\s*([\d.-]+)\s*\)/);
        if (match) pairwiseY = parseFloat(match[2]);
      }
      const pairwiseBBox = pairwiseLegend.getBBox();
      const pairwiseBottom = pairwiseY + pairwiseBBox.y + pairwiseBBox.height;
      legendContentHeight = Math.max(legendContentHeight, pairwiseBottom);
    }

    const viewBox = svg.getAttribute('viewBox');
    if (!viewBox) return;
    const parts = viewBox.split(/\s+/).map(parseFloat);
    if (parts.length !== 4) return;
    let [vbX, vbY, vbW, vbH] = parts;

    const legendBottomEdge = legendGroupY + legendContentHeight;
    const bottomPadding = 20;

    console.log(
      `Canvas expansion check: legendGroupY=${legendGroupY.toFixed(1)}, legendContentHeight=${legendContentHeight.toFixed(
        1
      )}, legendBottomEdge=${legendBottomEdge.toFixed(1)}, vbH=${vbH.toFixed(1)}`
    );

    if (legendBottomEdge + bottomPadding > vbH) {
      const newVbH = legendBottomEdge + bottomPadding;
      console.log(`Expanding canvas for vertical legend: ${vbH.toFixed(1)} -> ${newVbH.toFixed(1)}`);

      svg.setAttribute('viewBox', `${vbX} ${vbY} ${vbW} ${newVbH}`);
    }
  };

  const updatePairwiseLegendPositions = (svg) => {
    const legendGroup = svg.getElementById('legend');
    if (!legendGroup) return;

    let rectSize = 14;
    const firstColorRect = legendGroup.querySelector('path[fill]:not([fill="none"]):not([fill^="url("])');
    if (firstColorRect) {
      const d = firstColorRect.getAttribute('d');
      if (d) {
        const lMatch = d.match(/L\s+([\d.]+),/);
        if (lMatch) rectSize = parseFloat(lMatch[1]);
      }
    }
    const lineMargin = (24 / 14) * rectSize;

    const verticalLegend = legendGroup.querySelector('#legend_vertical');
    if (verticalLegend) {
      const vFeatureLegend = verticalLegend.querySelector('#feature_legend_v');
      const vPairwiseLegend = verticalLegend.querySelector('#pairwise_legend');

      if (vFeatureLegend && vPairwiseLegend) {
        let maxFeatureY = 0;
        const featureTexts = vFeatureLegend.querySelectorAll('text');
        featureTexts.forEach((el) => {
          const pos = parseTransformXY(el.getAttribute('transform'));
          if (pos.y > maxFeatureY) maxFeatureY = pos.y;
        });

        const newPairwiseY = maxFeatureY + lineMargin + lineMargin / 2;

        const currentTransform = vPairwiseLegend.getAttribute('transform');
        let pairwiseX = 0;
        if (currentTransform) {
          const match = currentTransform.match(/translate\(\s*([\d.-]+)\s*,/);
          if (match) pairwiseX = parseFloat(match[1]);
        }

        vPairwiseLegend.setAttribute('transform', `translate(${pairwiseX}, ${newPairwiseY})`);
        console.log(`Repositioned vertical pairwise legend to y=${newPairwiseY}`);
      }
    }

    const horizontalLegend = legendGroup.querySelector('#legend_horizontal');
    if (horizontalLegend) {
      const hFeatureLegend = horizontalLegend.querySelector('#feature_legend_h');
      const hPairwiseLegend = horizontalLegend.querySelector('#pairwise_legend');

      if (hFeatureLegend && hPairwiseLegend) {
        let minFeatureY = Infinity,
          maxFeatureY = -Infinity;
        const featureTexts = hFeatureLegend.querySelectorAll('text');
        featureTexts.forEach((el) => {
          const pos = parseTransformXY(el.getAttribute('transform'));
          if (pos.y < minFeatureY) minFeatureY = pos.y;
          if (pos.y > maxFeatureY) maxFeatureY = pos.y;
        });

        if (minFeatureY !== Infinity && maxFeatureY !== -Infinity) {
          const featureHeight = maxFeatureY - minFeatureY + lineMargin;

          const pairwiseBBox = hPairwiseLegend.getBBox();
          const pairwiseHeight = pairwiseBBox.height;

          const currentTransform = hPairwiseLegend.getAttribute('transform');
          let pairwiseX = 0,
            currentPairwiseY = 0;
          if (currentTransform) {
            const match = currentTransform.match(/translate\(\s*([\d.-]+)\s*,\s*([\d.-]+)\s*\)/);
            if (match) {
              pairwiseX = parseFloat(match[1]);
              currentPairwiseY = parseFloat(match[2]);
            }
          }

          let newPairwiseY = currentPairwiseY;
          if (featureHeight > pairwiseHeight) {
            newPairwiseY = (featureHeight - pairwiseHeight) / 2;
          }

          hPairwiseLegend.setAttribute('transform', `translate(${pairwiseX}, ${newPairwiseY})`);
          console.log(`Repositioned horizontal pairwise legend to y=${newPairwiseY}`);
        }
      }
    }

    expandCanvasForVerticalLegend(svg);
  };

  const reflowDualLegendLayout = (svg) => {
    const legendGroup = svg.getElementById('legend');
    if (!legendGroup) return;

    const horizontalLegend = legendGroup.querySelector('#legend_horizontal');
    const verticalLegend = legendGroup.querySelector('#legend_vertical');
    if (!horizontalLegend || !verticalLegend) return;

    const viewBox = svg.getAttribute('viewBox');
    let baseWidth = 800;
    if (viewBox) {
      const parts = viewBox.split(/\s+/).map(parseFloat);
      if (parts.length === 4) baseWidth = parts[2];
    }

    const horizontalViewBox = svg.getAttribute('data-horizontal-viewbox');
    let horizontalWidth = baseWidth;
    if (horizontalViewBox) {
      const parts = horizontalViewBox.split(/\s+/).map(parseFloat);
      if (parts.length === 4) horizontalWidth = parts[2];
    }

    const layoutLegendGroup = (legend, layout, maxWidth) => {
      const featureGroup =
        layout === 'horizontal'
          ? legend.querySelector('#feature_legend_h')
          : legend.querySelector('#feature_legend_v');
      if (!featureGroup) return;

      const pairwiseLegend = legend.querySelector('#pairwise_legend');

      featureGroup.setAttribute('transform', 'translate(0, 0)');
      if (pairwiseLegend) {
        pairwiseLegend.setAttribute('transform', 'translate(0, 0)');
      }

      let rectSize = 14;
      const firstColorRect = featureGroup.querySelector('path[fill]:not([fill="none"]):not([fill^="url("])');
      if (firstColorRect) {
        const d = firstColorRect.getAttribute('d');
        if (d) {
          const lMatch = d.match(/L\s+([\d.]+),/);
          if (lMatch) rectSize = parseFloat(lMatch[1]);
        }
      }

      const lineHeight = (24 / 14) * rectSize;
      const textXOffset = (22 / 14) * rectSize;

      const texts = Array.from(featureGroup.querySelectorAll('text'));
      const rects = Array.from(featureGroup.querySelectorAll('path')).filter((r) => {
        const fill = r.getAttribute('fill');
        return fill && fill !== 'none' && !fill.startsWith('url(');
      });

      const entries = texts.map((t) => {
        const pos = parseTransformXY(t.getAttribute('transform'));
        const expectedRectX = pos.x - textXOffset;
        let matchedRect = null;
        let bestDistance = Infinity;
        for (const r of rects) {
          const rectPos = parseTransformXY(r.getAttribute('transform'));
          if (Math.abs(rectPos.y - pos.y) < 1) {
            const distance = Math.abs(rectPos.x - expectedRectX);
            if (distance < textXOffset && distance < bestDistance) {
              bestDistance = distance;
              matchedRect = r;
            }
          }
        }
        return { text: t, rect: matchedRect, x: pos.x, y: pos.y };
      });

      if (layout === 'horizontal') {
        entries.sort((a, b) => {
          if (Math.abs(a.y - b.y) < 1) return a.x - b.x;
          return a.y - b.y;
        });

        let newX = textXOffset;
        let newY = rectSize / 2;
        const wrapWidth = maxWidth || baseWidth;

        entries.forEach((entry) => {
          const textBBox = entry.text.getBBox();
          const entryWidth = rectSize + textXOffset + textBBox.width + textXOffset;

          if (newX + entryWidth > wrapWidth && newX > textXOffset) {
            newX = textXOffset;
            newY += lineHeight;
          }

          entry.text.setAttribute('transform', `translate(${newX}, ${newY})`);
          if (entry.rect) {
            entry.rect.setAttribute('transform', `translate(${newX - textXOffset}, ${newY})`);
          }

          newX += entryWidth;
        });
      } else {
        entries.sort((a, b) => a.y - b.y);

        let newY = rectSize / 2;
        entries.forEach((entry) => {
          entry.text.setAttribute('transform', `translate(${textXOffset}, ${newY})`);
          if (entry.rect) {
            entry.rect.setAttribute('transform', `translate(0, ${newY})`);
          }
          newY += lineHeight;
        });
      }

      const featureBBox = featureGroup.getBBox();
      const featureWidth = featureBBox.width || 0;
      let featureHeight = featureBBox.height || 0;
      if (featureHeight === 0) {
        featureHeight = layout === 'horizontal' ? lineHeight : rectSize / 2;
      }

      if (pairwiseLegend) {
        const pairwiseBBox = pairwiseLegend.getBBox();
        const pairwiseWidth = pairwiseBBox.width || 0;
        const pairwiseHeight = pairwiseBBox.height || 0;

        if (layout === 'horizontal') {
          const featureYOffset = pairwiseHeight > featureHeight ? (pairwiseHeight - featureHeight) / 2 : 0;
          const pairwiseYOffset = featureHeight > pairwiseHeight ? (featureHeight - pairwiseHeight) / 2 : 0;
          if (featureYOffset > 0) {
            featureGroup.setAttribute('transform', `translate(0, ${featureYOffset})`);
          }
          const pairwiseX = featureWidth + textXOffset;
          pairwiseLegend.setAttribute('transform', `translate(${pairwiseX}, ${pairwiseYOffset})`);
        } else {
          const featureXOffset = pairwiseWidth > featureWidth ? (pairwiseWidth - featureWidth) / 2 : 0;
          const pairwiseXOffset = featureWidth > pairwiseWidth ? (featureWidth - pairwiseWidth) / 2 : 0;
          if (featureXOffset > 0) {
            featureGroup.setAttribute('transform', `translate(${featureXOffset}, 0)`);
          }
          const pairwiseY = featureHeight + lineHeight / 2;
          pairwiseLegend.setAttribute('transform', `translate(${pairwiseXOffset}, ${pairwiseY})`);
        }
      }
    };

    const horizontalDisplay = horizontalLegend.getAttribute('display');
    const verticalDisplay = verticalLegend.getAttribute('display');
    horizontalLegend.removeAttribute('display');
    verticalLegend.removeAttribute('display');

    layoutLegendGroup(horizontalLegend, 'horizontal', horizontalWidth);
    layoutLegendGroup(verticalLegend, 'vertical', baseWidth);

    if (horizontalDisplay !== null) {
      horizontalLegend.setAttribute('display', horizontalDisplay);
    } else {
      horizontalLegend.removeAttribute('display');
    }

    if (verticalDisplay !== null) {
      verticalLegend.setAttribute('display', verticalDisplay);
    } else {
      verticalLegend.removeAttribute('display');
    }

    expandCanvasForVerticalLegend(svg);
  };

  const compactLegendEntries = (svg) => {
    const legendGroup = svg.getElementById('legend');
    const hasDualLegends =
      !!legendGroup?.querySelector('#legend_horizontal') && !!legendGroup?.querySelector('#legend_vertical');
    if (hasDualLegends) {
      reflowDualLegendLayout(svg);
      return;
    }

    const targetGroups = getAllFeatureLegendGroups(svg);
    if (targetGroups.length === 0) return;

    for (const targetGroup of targetGroups) {
      let rectSize = 14;
      const firstColorRect = targetGroup.querySelector(
        'path[fill]:not([fill="none"]):not([fill^="url("])'
      );
      if (firstColorRect) {
        const d = firstColorRect.getAttribute('d');
        if (d) {
          const lMatch = d.match(/L\s+([\d.]+),/);
          if (lMatch) rectSize = parseFloat(lMatch[1]);
        }
      }
      const lineHeight = (24 / 14) * rectSize;
      const textXOffset = (22 / 14) * rectSize;

      const texts = Array.from(targetGroup.querySelectorAll('text'));
      if (texts.length === 0) continue;

      const allRects = Array.from(targetGroup.querySelectorAll('path')).filter((r) => {
        const fill = r.getAttribute('fill');
        return fill && fill !== 'none' && !fill.startsWith('url(');
      });

      const entries = texts.map((t) => {
        const pos = parseTransformXY(t.getAttribute('transform'));
        const expectedRectX = pos.x - textXOffset;
        let matchedRect = null;
        let bestDistance = Infinity;
        for (const r of allRects) {
          const rectPos = parseTransformXY(r.getAttribute('transform'));
          if (Math.abs(rectPos.y - pos.y) < 1) {
            const distance = Math.abs(rectPos.x - expectedRectX);
            if (distance < textXOffset && distance < bestDistance) {
              bestDistance = distance;
              matchedRect = r;
            }
          }
        }
        return { text: t, rect: matchedRect, x: pos.x, y: pos.y };
      });

      const uniqueYs = [...new Set(entries.map((e) => Math.round(e.y)))];
      const isHorizontal = uniqueYs.length < entries.length;

      if (isHorizontal) {
        entries.sort((a, b) => {
          if (Math.abs(a.y - b.y) < 1) return a.x - b.x;
          return a.y - b.y;
        });

        let maxWidth = 800;
        const legendGroup = svg.getElementById('legend');
        if (legendGroup) {
          const horizontalLegend = legendGroup.querySelector('#legend_horizontal');
          if (horizontalLegend) {
            const hPairwiseLegend = horizontalLegend.querySelector('#pairwise_legend');
            if (hPairwiseLegend) {
              const pairwiseTransform = hPairwiseLegend.getAttribute('transform');
              if (pairwiseTransform) {
                const match = pairwiseTransform.match(/translate\(\s*([\d.-]+)\s*,/);
                if (match) {
                  maxWidth = parseFloat(match[1]) - textXOffset;
                }
              }
            } else {
              const viewBox = svg.getAttribute('viewBox');
              if (viewBox) {
                const parts = viewBox.split(/\s+/).map(parseFloat);
                if (parts.length === 4) maxWidth = parts[2] - textXOffset;
              }
            }
          }
        }

        let newX = textXOffset;
        let newY = rectSize / 2;

        entries.forEach((entry) => {
          const textBBox = entry.text.getBBox();
          const entryWidth = rectSize + textXOffset + textBBox.width + textXOffset;

          if (newX + entryWidth > maxWidth && newX > textXOffset) {
            newX = textXOffset;
            newY += lineHeight;
          }

          entry.text.setAttribute('transform', `translate(${newX}, ${newY})`);

          if (entry.rect) {
            const expectedRectX = newX - textXOffset;
            entry.rect.setAttribute('transform', `translate(${expectedRectX}, ${newY})`);
          }

          newX += entryWidth;
        });
      } else {
        entries.sort((a, b) => a.y - b.y);

        let newY = rectSize / 2;
        entries.forEach((entry) => {
          entry.text.setAttribute('transform', `translate(${textXOffset}, ${newY})`);

          if (entry.rect) {
            entry.rect.setAttribute('transform', `translate(0, ${newY})`);
          }

          newY += lineHeight;
        });
      }
    }

    updatePairwiseLegendPositions(svg);
  };

  const reflowSingleLegendLayout = (svg, layout, maxWidthOverride = null) => {
    const legendGroup = svg.getElementById('legend');
    if (!legendGroup) return null;

    const featureLegendGroup = legendGroup.querySelector('#feature_legend') || legendGroup;
    const pairwiseLegend = legendGroup.querySelector('#pairwise_legend');

    const textElements = Array.from(featureLegendGroup.querySelectorAll('text')).filter((el) => {
      if (!pairwiseLegend) return true;
      return !el.closest('#pairwise_legend');
    });
    if (textElements.length === 0) return null;

    let rectSize = 14;
    const colorRects = Array.from(featureLegendGroup.querySelectorAll('path')).filter((r) => {
      const fill = r.getAttribute('fill');
      if (!fill || fill === 'none' || fill.startsWith('url(')) return false;
      if (pairwiseLegend && r.closest('#pairwise_legend')) return false;
      return true;
    });
    if (colorRects.length > 0) {
      const d = colorRects[0].getAttribute('d');
      if (d) {
        const lMatch = d.match(/L\s+([\d.]+),/);
        if (lMatch) rectSize = parseFloat(lMatch[1]);
      }
    }

    const lineHeight = (24 / 14) * rectSize;
    const textXOffset = (22 / 14) * rectSize;
    const featureOffset = parseTransformXY(featureLegendGroup.getAttribute('transform'));
    const pairwiseBBox = pairwiseLegend ? pairwiseLegend.getBBox() : null;
    const pairwiseWidth = pairwiseBBox ? pairwiseBBox.width : 0;

    const entries = textElements.map((t) => {
      const pos = parseTransformXY(t.getAttribute('transform'));
      return { text: t, rect: null, x: pos.x, y: pos.y };
    });
    entries.forEach((entry) => {
      const expectedRectX = entry.x - textXOffset;
      let matchedRect = null;
      let bestDistance = Infinity;
      for (const r of colorRects) {
        const rectPos = parseTransformXY(r.getAttribute('transform'));
        if (Math.abs(rectPos.y - entry.y) < 1) {
          const distance = Math.abs(rectPos.x - expectedRectX);
          if (distance < textXOffset && distance < bestDistance) {
            bestDistance = distance;
            matchedRect = r;
          }
        }
      }
      entry.rect = matchedRect;
    });

    const computeFeatureBounds = () => {
      let minX = Infinity,
        minY = Infinity,
        maxX = -Infinity,
        maxY = -Infinity;
      entries.forEach((entry) => {
        if (entry.newX === undefined || entry.newY === undefined) return;
        const textWidth = entry.textWidth || 0;
        const rectX = entry.newX - textXOffset;
        const rectY = entry.newY - rectSize / 2;
        const entryRight = entry.newX + textWidth + textXOffset;
        const entryBottom = entry.newY + rectSize / 2;
        minX = Math.min(minX, rectX);
        minY = Math.min(minY, rectY);
        maxX = Math.max(maxX, entryRight);
        maxY = Math.max(maxY, entryBottom);
      });
      if (minX === Infinity) return null;
      return { x: minX, y: minY, width: maxX - minX, height: maxY - minY };
    };

    if (layout === 'horizontal') {
      let maxWidth = maxWidthOverride || 800;
      if (!maxWidthOverride) {
        const viewBox = svg.getAttribute('viewBox');
        if (viewBox) {
          const parts = viewBox.split(/\s+/).map(parseFloat);
          if (parts.length === 4) maxWidth = parts[2];
        }
      }

      const reservedOffset = Math.max(featureOffset.x, 0);
      let maxFeatureWidth = maxWidth - reservedOffset;
      if (pairwiseLegend && pairwiseWidth > 0) {
        maxFeatureWidth = maxWidth - reservedOffset - pairwiseWidth - textXOffset;
        if (maxFeatureWidth < rectSize + textXOffset * 2) {
          maxFeatureWidth = maxWidth - reservedOffset;
        }
      }

      entries.sort((a, b) => {
        if (Math.abs(a.y - b.y) < 1) return a.x - b.x;
        return a.y - b.y;
      });

      let newX = textXOffset;
      let newY = rectSize / 2;

      entries.forEach((entry) => {
        const textBBox = entry.text.getBBox();
        const entryWidth = rectSize + textXOffset + textBBox.width + textXOffset;

        if (newX + entryWidth > maxFeatureWidth && newX > textXOffset) {
          newX = textXOffset;
          newY += lineHeight;
        }

        entry.text.setAttribute('transform', `translate(${newX}, ${newY})`);
        if (entry.rect) {
          entry.rect.setAttribute('transform', `translate(${newX - textXOffset}, ${newY})`);
        }
        entry.newX = newX;
        entry.newY = newY;
        entry.textWidth = textBBox.width;

        newX += entryWidth;
      });
    } else {
      entries.sort((a, b) => a.y - b.y);
      let newY = rectSize / 2;
      entries.forEach((entry) => {
        const textBBox = entry.text.getBBox();
        entry.text.setAttribute('transform', `translate(${textXOffset}, ${newY})`);
        if (entry.rect) {
          entry.rect.setAttribute('transform', `translate(0, ${newY})`);
        }
        entry.newX = textXOffset;
        entry.newY = newY;
        entry.textWidth = textBBox.width;
        newY += lineHeight;
      });
    }

    if (pairwiseLegend) {
      const featureBounds = computeFeatureBounds();
      const fallbackBBox = featureLegendGroup.getBBox();
      const featureBBox = featureBounds || fallbackBBox;

      const featureX = featureOffset.x + (featureBounds ? featureBounds.x : fallbackBBox.x);
      const featureY = featureOffset.y + (featureBounds ? featureBounds.y : fallbackBBox.y);
      const featureWidth = featureBBox.width;
      const featureHeight = featureBBox.height;
      const effectivePairwiseBBox = pairwiseBBox || pairwiseLegend.getBBox();

      let pairwiseX = featureX;
      let pairwiseY = featureY;

      if (layout === 'horizontal') {
        const yOffset = Math.max((featureHeight - effectivePairwiseBBox.height) / 2, 0);
        pairwiseX = featureX + featureWidth + textXOffset;
        pairwiseY = featureY + yOffset;
      } else {
        const xOffset = Math.max((featureWidth - effectivePairwiseBBox.width) / 2, 0);
        pairwiseX = featureX + xOffset;
        pairwiseY = featureY + featureHeight + lineHeight / 2;
      }

      pairwiseLegend.setAttribute('transform', `translate(${pairwiseX}, ${pairwiseY})`);
    }

    const bbox = legendGroup.getBBox();
    return { legendWidth: bbox.width, legendHeight: bbox.height };
  };

  const extractLegendEntries = () => {
    if (!svgContainer.value) {
      legendEntries.value = [];
      return;
    }

    const svg = svgContainer.value.querySelector('svg');
    if (!svg) {
      legendEntries.value = [];
      return;
    }

    const targetGroup = getVisibleFeatureLegendGroup(svg);
    if (!targetGroup) {
      legendEntries.value = [];
      return;
    }

    const entries = [];

    const parseTransform = (transform) => {
      if (!transform) return { x: 0, y: 0 };
      const match = transform.match(/translate\(\s*([\d.-]+)\s*,\s*([\d.-]+)\s*\)/);
      return match ? { x: parseFloat(match[1]), y: parseFloat(match[2]) } : { x: 0, y: 0 };
    };

    let entryGroups = targetGroup.querySelectorAll('g[data-legend-key]');
    if (entryGroups.length === 0) {
      const texts = Array.from(targetGroup.querySelectorAll('text'));
      const allPaths = Array.from(targetGroup.querySelectorAll('path'));
      const groupsToAdd = [];

      texts.forEach((textEl) => {
        const caption = textEl.textContent?.trim();
        if (!caption) return;

        const textPos = parseTransform(textEl.getAttribute('transform'));
        let bestPath = null;
        let bestX = -Infinity;

        for (const path of allPaths) {
          const fill = path.getAttribute('fill');
          if (!fill || fill === 'none' || fill.startsWith('url(')) continue;
          const pathPos = parseTransform(path.getAttribute('transform'));
          if (Math.abs(pathPos.y - textPos.y) < 2 && pathPos.x < textPos.x) {
            if (pathPos.x > bestX) {
              bestX = pathPos.x;
              bestPath = path;
            }
          }
        }

        if (bestPath) {
          const entryGroup = document.createElementNS('http://www.w3.org/2000/svg', 'g');
          entryGroup.setAttribute('data-legend-key', caption);
          entryGroup.appendChild(bestPath);
          entryGroup.appendChild(textEl);
          groupsToAdd.push(entryGroup);
        }
      });

      if (groupsToAdd.length > 0) {
        groupsToAdd.forEach((group) => targetGroup.appendChild(group));
        skipCaptureBaseConfig.value = true;
        const resultIdx = selectedResultIndex.value;
        if (resultIdx >= 0 && results.value.length > resultIdx) {
          const serializer = new XMLSerializer();
          results.value[resultIdx] = { ...results.value[resultIdx], content: serializer.serializeToString(svg) };
        }
        entryGroups = targetGroup.querySelectorAll('g[data-legend-key]');
      }
    }

    entryGroups.forEach((entryGroup) => {
      const caption = entryGroup.getAttribute('data-legend-key');
      if (!caption) return;

      let color = '#cccccc';
      const paths = entryGroup.querySelectorAll('path');
      for (const path of paths) {
        const fill = path.getAttribute('fill');
        if (fill && fill !== 'none' && !fill.startsWith('url(')) {
          color = fill;
          break;
        }
      }

      let xPos = 0,
        yPos = 0;
      const groupTransform = parseTransform(entryGroup.getAttribute('transform'));
      if (groupTransform.x !== 0 || groupTransform.y !== 0) {
        xPos = groupTransform.x;
        yPos = groupTransform.y;
      } else {
        const textEl = entryGroup.querySelector('text');
        if (textEl) {
          const textTransform = parseTransform(textEl.getAttribute('transform'));
          xPos = textTransform.x;
          yPos = textTransform.y;
        }
      }

      const existingEntry = legendEntries.value.find((e) => e.caption === caption);
      const showStroke = existingEntry?.showStroke || false;
      const existingFeatureIds = existingEntry?.featureIds || [];
      const originalCaption = existingEntry?.originalCaption || caption;

      entries.push({
        caption,
        originalCaption,
        color,
        xPos,
        yPos,
        showStroke,
        featureIds: existingFeatureIds
      });
    });

    legendEntries.value = entries;

    if (originalLegendOrder.value.length === 0 && entries.length > 0) {
      originalLegendOrder.value = entries.map((e) => e.caption);
    }

    if (Object.keys(originalLegendColors.value).length === 0 && entries.length > 0) {
      entries.forEach((entry) => {
        originalLegendColors.value[entry.caption] = entry.color;
      });
    }
  };

  const updateLegendEntryColor = (idx, newColor) => {
    if (!svgContainer.value) return;
    const svg = svgContainer.value.querySelector('svg');
    if (!svg) return;

    const entry = legendEntries.value[idx];
    if (!entry) return;

    const caption = entry.caption;
    if (!caption) return;

    const targetGroups = getAllFeatureLegendGroups(svg);
    if (targetGroups.length === 0) return;

    for (const targetGroup of targetGroups) {
      const entryGroup = targetGroup.querySelector(`g[data-legend-key="${CSS.escape(caption)}"]`);
      if (entryGroup) {
        const paths = entryGroup.querySelectorAll('path');
        for (const path of paths) {
          const fill = path.getAttribute('fill');
          if (fill && fill !== 'none' && !fill.startsWith('url(')) {
            path.setAttribute('fill', newColor);
            break;
          }
        }
      }
    }

    entry.color = newColor;
    legendColorOverrides[caption] = newColor;

    skipCaptureBaseConfig.value = true;
    const resultIdx = selectedResultIndex.value;
    if (resultIdx >= 0 && results.value.length > resultIdx) {
      const serializer = new XMLSerializer();
      results.value[resultIdx] = { ...results.value[resultIdx], content: serializer.serializeToString(svg) };
    }
  };

  const updateLegendEntryCaption = (idx, newCaption) => {
    if (!svgContainer.value) return;
    const svg = svgContainer.value.querySelector('svg');
    if (!svg) return;

    const entry = legendEntries.value[idx];
    if (!entry) return;

    const oldCaption = entry.caption;
    const caption = newCaption.trim();
    if (!caption) return;

    const targetGroups = getAllFeatureLegendGroups(svg);
    if (targetGroups.length === 0) return;

    for (const targetGroup of targetGroups) {
      const entryGroup = targetGroup.querySelector(`g[data-legend-key="${CSS.escape(oldCaption)}"]`);
      if (entryGroup) {
        entryGroup.setAttribute('data-legend-key', caption);
        const textEl = entryGroup.querySelector('text');
        if (textEl) {
          textEl.textContent = caption;
        }
      }
    }

    if (legendColorOverrides[oldCaption]) {
      legendColorOverrides[caption] = legendColorOverrides[oldCaption];
      delete legendColorOverrides[oldCaption];
    }

    if (legendStrokeOverrides[oldCaption]) {
      legendStrokeOverrides[caption] = legendStrokeOverrides[oldCaption];
      delete legendStrokeOverrides[oldCaption];
    }

    const ruleMatches = manualSpecificRules.filter((r) => r.cap === oldCaption);
    for (const rule of ruleMatches) {
      rule.cap = caption;
    }

    entry.caption = caption;

    skipCaptureBaseConfig.value = true;
    const resultIdx = selectedResultIndex.value;
    if (resultIdx >= 0 && results.value.length > resultIdx) {
      const serializer = new XMLSerializer();
      results.value[resultIdx] = { ...results.value[resultIdx], content: serializer.serializeToString(svg) };
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
        const elements = svg.querySelectorAll(`#${CSS.escape(svgId)}`);
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
        const featurePaths = svg.querySelectorAll('path[id^="f"]');
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
        const serializer = new XMLSerializer();
        results.value[resultIdx] = { ...results.value[resultIdx], content: serializer.serializeToString(svg) };
      }
    }
  };

  const resetAllStrokes = () => {
    if (!svgContainer.value) return;
    const svg = svgContainer.value.querySelector('svg');
    if (!svg) return;

    const originalColor = originalSvgStroke.value.color;
    const originalWidth = originalSvgStroke.value.width;

    const featurePaths = svg.querySelectorAll('path[id^="f"]');
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

    if (updatedCount > 0) {
      skipCaptureBaseConfig.value = true;
      const resultIdx = selectedResultIndex.value;
      if (resultIdx >= 0 && results.value.length > resultIdx) {
        const serializer = new XMLSerializer();
        results.value[resultIdx] = { ...results.value[resultIdx], content: serializer.serializeToString(svg) };
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
      const el = svg.querySelector(`#${CSS.escape(feat.svg_id)}`);
      if (el) {
        return {
          strokeColor: el.getAttribute('stroke') || '#000000',
          strokeWidth: parseFloat(el.getAttribute('stroke-width')) || 0.5
        };
      }
    }

    if (legendFillColor) {
      const normalizedLegendColor = legendFillColor.toLowerCase();
      const featurePaths = svg.querySelectorAll('path[id^="f"]');
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

  const applyStrokeToFeaturesByCaption = (caption, strokeColor, strokeWidth) => {
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
        const elements = svg.querySelectorAll(`#${CSS.escape(svgId)}`);
        elements.forEach((el) => {
          if (strokeColor !== null) {
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
        const featurePaths = svg.querySelectorAll('path[id^="f"]');
        featurePaths.forEach((path) => {
          const fill = path.getAttribute('fill');
          if (fill && fill.toLowerCase() === normalizedLegendColor) {
            const pathId = path.getAttribute('id');
            if (pathId && !processedIds.has(pathId)) {
              processedIds.add(pathId);
              if (strokeColor !== null) {
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
            if (strokeColor !== null) {
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
        const serializer = new XMLSerializer();
        results.value[resultIdx] = { ...results.value[resultIdx], content: serializer.serializeToString(svg) };
      }
      console.log(`Applied stroke to ${updatedCount} elements for caption "${caption}"`);
    }
  };

  const deleteLegendEntry = (idx) => {
    const entry = legendEntries.value[idx];
    if (!entry) return;

    deletedLegendEntries.value.push({ ...entry });

    removeLegendEntry(entry.caption);
    extractLegendEntries();
  };

  const restoreDeletedLegendEntries = () => {
    if (deletedLegendEntries.value.length === 0) return;

    for (const entry of deletedLegendEntries.value) {
      addLegendEntry(entry.caption, entry.color);
    }
    deletedLegendEntries.value = [];
    extractLegendEntries();
  };

  const moveLegendEntryUp = (idx) => {
    if (idx <= 0) return;
    swapLegendEntries(idx, idx - 1);
  };

  const moveLegendEntryDown = (idx) => {
    if (idx >= legendEntries.value.length - 1) return;
    swapLegendEntries(idx, idx + 1);
  };

  const sortLegendEntries = (direction = 'asc') => {
    if (!svgContainer.value) return;
    const svg = svgContainer.value.querySelector('svg');
    if (!svg) return;

    const entries = [...legendEntries.value];
    if (entries.length < 2) return;

    const yPositions = entries.map((e) => e.yPos).sort((a, b) => a - b);

    const sortedEntries = [...entries].sort((a, b) => {
      const cmp = a.caption.localeCompare(b.caption, undefined, { sensitivity: 'base' });
      return direction === 'asc' ? cmp : -cmp;
    });

    const targetGroups = getAllFeatureLegendGroups(svg);
    if (targetGroups.length === 0) return;

    const yMapping = new Map();
    entries.forEach((entry) => {
      const sortedIdx = sortedEntries.findIndex((e) => e.caption === entry.caption && e.color === entry.color);
      if (sortedIdx !== -1) {
        yMapping.set(entry.yPos, yPositions[sortedIdx]);
      }
    });

    for (const targetGroup of targetGroups) {
      const allElements = targetGroup.querySelectorAll('path, text');
      allElements.forEach((el) => {
        const transform = el.getAttribute('transform');
        if (!transform) return;

        const match = transform.match(/translate\(([^,]+),\s*([\d.-]+)\)/);
        if (!match) return;

        const x = match[1];
        const y = parseFloat(match[2]);

        for (const [oldY, newY] of yMapping) {
          if (Math.abs(y - oldY) < 1) {
            el.setAttribute('transform', `translate(${x}, ${newY})`);
            break;
          }
        }
      });
    }

    skipCaptureBaseConfig.value = true;
    const resultIdx = selectedResultIndex.value;
    if (resultIdx >= 0 && results.value.length > resultIdx) {
      const serializer = new XMLSerializer();
      results.value[resultIdx] = { ...results.value[resultIdx], content: serializer.serializeToString(svg) };
    }

    extractLegendEntries();
  };

  const sortLegendEntriesByDefault = () => {
    if (!svgContainer.value) return;
    const svg = svgContainer.value.querySelector('svg');
    if (!svg) return;

    const entries = [...legendEntries.value];
    if (entries.length < 2) return;
    if (originalLegendOrder.value.length === 0) return;

    const yPositions = entries.map((e) => e.yPos).sort((a, b) => a - b);

    const sortedEntries = [...entries].sort((a, b) => {
      const aOrigIdx = originalLegendOrder.value.indexOf(a.caption);
      const bOrigIdx = originalLegendOrder.value.indexOf(b.caption);

      if (aOrigIdx !== -1 && bOrigIdx !== -1) {
        return aOrigIdx - bOrigIdx;
      }
      if (aOrigIdx !== -1) return -1;
      if (bOrigIdx !== -1) return 1;
      return a.caption.localeCompare(b.caption, undefined, { sensitivity: 'base' });
    });

    const targetGroups = getAllFeatureLegendGroups(svg);
    if (targetGroups.length === 0) return;

    const yMapping = new Map();
    entries.forEach((entry) => {
      const sortedIdx = sortedEntries.findIndex((e) => e.caption === entry.caption && e.color === entry.color);
      if (sortedIdx !== -1) {
        yMapping.set(entry.yPos, yPositions[sortedIdx]);
      }
    });

    for (const targetGroup of targetGroups) {
      const allElements = targetGroup.querySelectorAll('path, text');
      allElements.forEach((el) => {
        const transform = el.getAttribute('transform');
        if (!transform) return;

        const match = transform.match(/translate\(([^,]+),\s*([\d.-]+)\)/);
        if (!match) return;

        const x = match[1];
        const y = parseFloat(match[2]);

        for (const [oldY, newY] of yMapping) {
          if (Math.abs(y - oldY) < 1) {
            el.setAttribute('transform', `translate(${x}, ${newY})`);
            break;
          }
        }
      });
    }

    skipCaptureBaseConfig.value = true;
    const resultIdx = selectedResultIndex.value;
    if (resultIdx >= 0 && results.value.length > resultIdx) {
      const serializer = new XMLSerializer();
      results.value[resultIdx] = { ...results.value[resultIdx], content: serializer.serializeToString(svg) };
    }

    extractLegendEntries();
  };

  const swapLegendEntries = (idx1, idx2) => {
    if (!svgContainer.value) return;
    const svg = svgContainer.value.querySelector('svg');
    if (!svg) return;

    const targetGroups = getAllFeatureLegendGroups(svg);
    if (targetGroups.length === 0) return;

    const entry1 = legendEntries.value[idx1];
    const entry2 = legendEntries.value[idx2];
    if (!entry1 || !entry2) return;

    const y1 = entry1.yPos;
    const y2 = entry2.yPos;

    for (const targetGroup of targetGroups) {
      const allElements = targetGroup.querySelectorAll('path, text');
      allElements.forEach((el) => {
        const transform = el.getAttribute('transform');
        if (!transform) return;

        const match = transform.match(/translate\(([^,]+),\s*([\d.-]+)\)/);
        if (!match) return;

        const x = match[1];
        const y = parseFloat(match[2]);

        if (Math.abs(y - y1) < 1) {
          el.setAttribute('transform', `translate(${x}, ${y2})`);
        } else if (Math.abs(y - y2) < 1) {
          el.setAttribute('transform', `translate(${x}, ${y1})`);
        }
      });
    }

    legendEntries.value[idx1].yPos = y2;
    legendEntries.value[idx2].yPos = y1;

    const temp = { ...legendEntries.value[idx1] };
    legendEntries.value[idx1] = { ...legendEntries.value[idx2], yPos: y1 };
    legendEntries.value[idx2] = { ...temp, yPos: y2 };

    skipCaptureBaseConfig.value = true;
    const resultIdx = selectedResultIndex.value;
    if (resultIdx >= 0 && results.value.length > resultIdx) {
      const serializer = new XMLSerializer();
      results.value[resultIdx] = { ...results.value[resultIdx], content: serializer.serializeToString(svg) };
    }

    extractLegendEntries();
  };

  const addNewLegendEntry = async () => {
    if (!newLegendCaption.value.trim()) return;

    const added = await addLegendEntry(newLegendCaption.value.trim(), newLegendColor.value);
    if (added) {
      newLegendCaption.value = '';
      newLegendColor.value = '#808080';
      setTimeout(() => extractLegendEntries(), 100);
    }
  };

  const parseTransform = (transformStr) => {
    if (!transformStr) return { x: 0, y: 0 };
    const match = transformStr.match(/translate\(\s*([-\d.]+)\s*,?\s*([-\d.]+)?\s*\)/);
    if (match) {
      return { x: parseFloat(match[1]) || 0, y: parseFloat(match[2]) || 0 };
    }
    return { x: 0, y: 0 };
  };

  const startLegendDrag = (e) => {
    if (!svgContainer.value) return;
    const svg = svgContainer.value.querySelector('svg');
    if (!svg) return;
    const legendGroup = svg.getElementById('legend');
    if (!legendGroup) return;

    e.preventDefault();
    e.stopPropagation();

    legendDragging.value = true;
    legendDragStart.x = e.clientX;
    legendDragStart.y = e.clientY;

    const currentTransform = parseTransform(legendGroup.getAttribute('transform'));
    legendOriginalTransform.value = { ...currentTransform };
  };

  const onLegendDrag = (e) => {
    if (!legendDragging.value) return;
    if (!svgContainer.value) return;

    const svg = svgContainer.value.querySelector('svg');
    if (!svg) return;
    const legendGroup = svg.getElementById('legend');
    if (!legendGroup) return;

    const deltaX = (e.clientX - legendDragStart.x) / zoom.value;
    const deltaY = (e.clientY - legendDragStart.y) / zoom.value;

    const newX = legendOriginalTransform.value.x + deltaX;
    const newY = legendOriginalTransform.value.y + deltaY;

    legendGroup.setAttribute('transform', `translate(${newX}, ${newY})`);
    legendCurrentOffset.x = deltaX;
    legendCurrentOffset.y = deltaY;
  };

  const endLegendDrag = () => {
    if (!legendDragging.value) return;

    legendDragging.value = false;
  };

  const resetLegendPositionOnly = () => {
    if (!svgContainer.value) return;
    const svg = svgContainer.value.querySelector('svg');
    if (!svg) return;
    const legendGroup = svg.getElementById('legend');
    if (!legendGroup) return;

    const initial = legendInitialTransform.value;
    legendGroup.setAttribute('transform', `translate(${initial.x}, ${initial.y})`);
    legendCurrentOffset.x = 0;
    legendCurrentOffset.y = 0;

    skipCaptureBaseConfig.value = true;
    const idx = selectedResultIndex.value;
    if (idx >= 0 && results.value.length > idx) {
      const serializer = new XMLSerializer();
      results.value[idx] = { ...results.value[idx], content: serializer.serializeToString(svg) };
    }
  };

  const resetLegendPosition = () => {
    resetLegendPositionOnly();
    extractLegendEntries();
  };

  const setupLegendDrag = () => {
    if (!svgContainer.value) return;
    const svg = svgContainer.value.querySelector('svg');
    if (!svg) return;
    const legendGroup = svg.getElementById('legend');
    if (!legendGroup) return;

    const initialTransform = parseTransform(legendGroup.getAttribute('transform'));
    legendInitialTransform.value = { ...initialTransform };

    legendGroup.style.cursor = 'move';
    legendGroup.onmousedown = startLegendDrag;

    svg.onmousemove = onLegendDrag;
    svg.onmouseup = endLegendDrag;
    svg.onmouseleave = endLegendDrag;
  };

  const reapplyStrokeOverrides = () => {
    if (!svgContainer.value) return;
    const svg = svgContainer.value.querySelector('svg');
    if (!svg) return;

    const overrideCount = Object.keys(legendStrokeOverrides).length;
    if (overrideCount === 0) {
      debugLog('No stroke overrides to reapply');
      return;
    }

    console.log(`[DEBUG] Reapplying ${overrideCount} stroke overrides to new SVG`);

    let totalUpdated = 0;

    for (const [caption, overrides] of Object.entries(legendStrokeOverrides)) {
      const { strokeColor, strokeWidth } = overrides;
      if (!strokeColor && strokeWidth === undefined) continue;

      const matchingFeatures = extractedFeatures.value.filter((f) => getFeatureCaption(f) === caption);

      for (const feat of matchingFeatures) {
        const elements = svg.querySelectorAll(`#${CSS.escape(feat.svg_id)}`);
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

    if (totalUpdated > 0) {
      skipCaptureBaseConfig.value = true;
      const resultIdx = selectedResultIndex.value;
      if (resultIdx >= 0 && results.value.length > resultIdx) {
        const serializer = new XMLSerializer();
        results.value[resultIdx] = { ...results.value[resultIdx], content: serializer.serializeToString(svg) };
      }
      console.log(`[DEBUG] Reapplied stroke overrides to ${totalUpdated} elements`);
    }
  };

  return {
    addLegendEntry,
    addNewLegendEntry,
    getAllFeatureLegendGroups,
    getVisibleFeatureLegendGroup,
    isCurrentLegendHorizontal,
    updateLegendEntryColorByCaption,
    legendEntryExists,
    removeLegendEntry,
    expandCanvasForVerticalLegend,
    updatePairwiseLegendPositions,
    compactLegendEntries,
    reflowSingleLegendLayout,
    extractLegendEntries,
    updateLegendEntryColor,
    updateLegendEntryCaption,
    getLegendEntryStrokeColor,
    getLegendEntryStrokeWidth,
    updateLegendEntryStrokeColor,
    updateLegendEntryStrokeWidth,
    resetLegendEntryStroke,
    resetAllStrokes,
    captureOriginalStrokeValues,
    applyStrokeToFeaturesByCaption,
    deleteLegendEntry,
    restoreDeletedLegendEntries,
    moveLegendEntryUp,
    moveLegendEntryDown,
    sortLegendEntries,
    sortLegendEntriesByDefault,
    swapLegendEntries,
    startLegendDrag,
    onLegendDrag,
    endLegendDrag,
    resetLegendPositionOnly,
    resetLegendPosition,
    setupLegendDrag,
    reapplyStrokeOverrides
  };
};
