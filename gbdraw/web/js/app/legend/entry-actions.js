import {
  getAllFeatureLegendGroups,
  getLegendChildById,
  getVisibleFeatureLegendGroup,
  parseTransformXY
} from './utils.js';

export const createLegendEntryActions = ({ state, getPyodide, layoutActions }) => {
  const {
    pyodideReady,
    results,
    selectedResultIndex,
    svgContainer,
    adv,
    legendEntries,
    deletedLegendEntries,
    originalLegendOrder,
    originalLegendColors,
    newLegendCaption,
    newLegendColor,
    legendStrokeOverrides,
    legendColorOverrides,
    manualSpecificRules,
    skipCaptureBaseConfig
  } = state;

  const { updatePairwiseLegendPositions, reflowDualLegendLayout, compactLegendEntries } = layoutActions;

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
            const hPairwiseLegend = getLegendChildById(horizontalLegend, 'pairwise_legend');
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

  const addNewLegendEntry = async () => {
    if (!newLegendCaption.value.trim()) return;

    const added = await addLegendEntry(newLegendCaption.value.trim(), newLegendColor.value);
    if (added) {
      newLegendCaption.value = '';
      newLegendColor.value = '#808080';
      setTimeout(() => extractLegendEntries(), 100);
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

  return {
    addLegendEntry,
    addNewLegendEntry,
    deleteLegendEntry,
    extractLegendEntries,
    legendEntryExists,
    removeLegendEntry,
    restoreDeletedLegendEntries,
    updateLegendEntryCaption,
    updateLegendEntryColor,
    updateLegendEntryColorByCaption
  };
};
