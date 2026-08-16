import {
  getAllFeatureLegendGroups,
  getComparisonLegendGroup,
  getLegendChildById,
  isInsideComparisonLegend,
  parseTransformXY
} from './utils.js';
import { parseCompositionMetadata } from '../legend-layout/composition-actions.js';

export const createLegendLayoutActions = () => {
  const setAttribute = (element, name, value, mutation = null) => (
    mutation
      ? mutation.setAttribute(element, name, value)
      : (element.setAttribute(name, value), true)
  );
  const removeAttribute = (element, name, mutation = null) => (
    mutation
      ? mutation.removeAttribute(element, name)
      : (element.removeAttribute(name), true)
  );
  const getHorizontalWrapWidth = (svg) => {
    if (!svg) return null;
    return parseCompositionMetadata(svg).primary.finalBounds.width;
  };

  const getLegendReflowMetrics = (svg) => {
    const metrics = parseCompositionMetadata(svg).legendReflow;
    if (!metrics) {
      throw new Error('This diagram has no legend reflow metadata. Regenerate it before editing the legend.');
    }
    return metrics;
  };

  const centerHorizontalRows = (entries, textXOffset, availableWidth, mutation = null) => {
    if (!Number.isFinite(availableWidth) || availableWidth <= 0 || !Array.isArray(entries) || entries.length === 0) {
      return;
    }

    const rows = new Map();
    entries.forEach((entry) => {
      if (entry.newX === undefined || entry.newY === undefined) return;
      const rowKey = Number(entry.newY).toFixed(3);
      if (!rows.has(rowKey)) rows.set(rowKey, []);
      rows.get(rowKey).push(entry);
    });

    rows.forEach((rowEntries) => {
      let rowLeft = Infinity;
      let rowRight = -Infinity;

      rowEntries.forEach((entry) => {
        const measuredWidth =
          entry.textWidth !== undefined && Number.isFinite(entry.textWidth)
            ? entry.textWidth
            : entry.text.getBBox().width;
        const left = Number(entry.newX) - textXOffset;
        const right = Number(entry.newX) + measuredWidth + textXOffset;
        rowLeft = Math.min(rowLeft, left);
        rowRight = Math.max(rowRight, right);
      });

      if (!Number.isFinite(rowLeft) || !Number.isFinite(rowRight) || rowRight <= rowLeft) return;

      const rowWidth = rowRight - rowLeft;
      const targetLeft = Math.max(0, (availableWidth - rowWidth) * 0.5);
      const shiftX = targetLeft - rowLeft;
      if (Math.abs(shiftX) < 1e-6) return;

      rowEntries.forEach((entry) => {
        const shiftedX = Number(entry.newX) + shiftX;
        const y = Number(entry.newY);
        entry.newX = shiftedX;
        setAttribute(entry.text, 'transform', `translate(${shiftedX}, ${y})`, mutation);
        if (entry.rect) {
          setAttribute(
            entry.rect,
            'transform',
            `translate(${shiftedX - textXOffset}, ${y})`,
            mutation
          );
        }
      });
    });
  };

  const getSingleComparisonLegendParts = (pairwiseLegend) => {
    const texts = Array.from(pairwiseLegend?.querySelectorAll?.('text') || []);
    const titleTexts = texts.filter((text) => {
      const label = String(text.textContent || '').trim();
      return label && !/^\d+(?:\.\d+)?%$/.test(label);
    });
    const barPaths = Array.from(pairwiseLegend?.querySelectorAll?.('path') || []).filter((path) =>
      String(path.getAttribute('fill') || '').startsWith('url(')
    );
    if (titleTexts.length !== 1 || barPaths.length !== 1) return null;

    const title = titleTexts[0];
    const bar = barPaths[0];
    const d = bar.getAttribute('d') || '';
    const widthMatch = d.match(/L\s*([+-]?\d+(?:\.\d+)?),/);
    const parsedBarWidth = widthMatch ? parseFloat(widthMatch[1]) : NaN;
    const barWidth = Number.isFinite(parsedBarWidth) && parsedBarWidth > 0 ? parsedBarWidth : bar.getBBox().width;
    const titleWidth = title.getBBox().width;
    const percentLabels = texts.filter((text) => /^\d+(?:\.\d+)?%$/.test(String(text.textContent || '').trim()));
    return { title, titleWidth, bar, barWidth, percentLabels };
  };

  const getComparisonLegendAlignmentWidth = (pairwiseLegend, fallbackWidth = 0) => {
    const parts = getSingleComparisonLegendParts(pairwiseLegend);
    if (parts) {
      return parts.barWidth;
    }

    return fallbackWidth;
  };

  const centerSingleComparisonLegendParts = (pairwiseLegend, mutation = null) => {
    const parts = getSingleComparisonLegendParts(pairwiseLegend);
    if (!parts) return;

    const { title, bar, barWidth, percentLabels } = parts;
    const barX = 0;
    const titleX = barWidth / 2;
    const titlePos = parseTransformXY(title.getAttribute('transform'));
    const barPos = parseTransformXY(bar.getAttribute('transform'));

    setAttribute(title, 'text-anchor', 'middle', mutation);
    setAttribute(title, 'transform', `translate(${titleX}, ${titlePos.y})`, mutation);
    setAttribute(bar, 'transform', `translate(${barX}, ${barPos.y})`, mutation);

    percentLabels.forEach((label) => {
      const labelPos = parseTransformXY(label.getAttribute('transform'));
      const text = String(label.textContent || '').trim();
      const x = text === '100%' ? barX + barWidth : barX;
      setAttribute(label, 'transform', `translate(${x}, ${labelPos.y})`, mutation);
    });
  };

  const updatePairwiseLegendPositions = (svg, mutation = null) => {
    const legendGroup = svg.getElementById('legend');
    if (!legendGroup) return;

    const horizontalLegend = legendGroup.querySelector('#legend_horizontal');
    const verticalLegend = legendGroup.querySelector('#legend_vertical');
    const hasDualLegends = !!(horizontalLegend && verticalLegend);

    if (!hasDualLegends) {
      const legendSide = parseCompositionMetadata(svg).legendSide;
      const layout = legendSide === 'top' || legendSide === 'bottom' ? 'horizontal' : 'vertical';
      const maxWidth = layout === 'horizontal' ? getHorizontalWrapWidth(svg) : null;
      reflowSingleLegendLayout(svg, layout, maxWidth, mutation);
      return;
    }

    const { lineHeight: lineMargin } = getLegendReflowMetrics(svg);

    if (verticalLegend) {
      const vFeatureLegend = verticalLegend.querySelector('#feature_legend_v');
      const vPairwiseLegend = getLegendChildById(verticalLegend, 'pairwise_legend');

      if (vFeatureLegend && vPairwiseLegend) {
        centerSingleComparisonLegendParts(vPairwiseLegend, mutation);

        let maxFeatureY = 0;
        const featureTexts = vFeatureLegend.querySelectorAll('text');
        featureTexts.forEach((el) => {
          const pos = parseTransformXY(el.getAttribute('transform'));
          if (pos.y > maxFeatureY) maxFeatureY = pos.y;
        });

        const newPairwiseY = maxFeatureY + lineMargin + lineMargin / 2;

        setAttribute(vPairwiseLegend, 'transform', `translate(0, ${newPairwiseY})`, mutation);
        console.log(`Repositioned vertical pairwise legend to y=${newPairwiseY}`);
      }
    }

    if (horizontalLegend) {
      const hFeatureLegend = horizontalLegend.querySelector('#feature_legend_h');
      const hPairwiseLegend = getLegendChildById(horizontalLegend, 'pairwise_legend');

      if (hFeatureLegend && hPairwiseLegend) {
        centerSingleComparisonLegendParts(hPairwiseLegend, mutation);

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

          setAttribute(
            hPairwiseLegend,
            'transform',
            `translate(${pairwiseX}, ${newPairwiseY})`,
            mutation
          );
          console.log(`Repositioned horizontal pairwise legend to y=${newPairwiseY}`);
        }
      }
    }

  };

  const reflowDualLegendLayout = (svg, mutation = null) => {
    const legendGroup = svg.getElementById('legend');
    if (!legendGroup) return;

    const horizontalLegend = legendGroup.querySelector('#legend_horizontal');
    const verticalLegend = legendGroup.querySelector('#legend_vertical');
    if (!horizontalLegend || !verticalLegend) return;

    const horizontalWidth = getHorizontalWrapWidth(svg);
    const reflowMetrics = getLegendReflowMetrics(svg);

    const layoutLegendGroup = (legend, layout, maxWidth) => {
      const featureGroup =
        layout === 'horizontal'
          ? legend.querySelector('#feature_legend_h')
          : legend.querySelector('#feature_legend_v');
      if (!featureGroup) return;

      const pairwiseLegend = getLegendChildById(legend, 'pairwise_legend');

      setAttribute(featureGroup, 'transform', 'translate(0, 0)', mutation);
      if (pairwiseLegend) {
        setAttribute(pairwiseLegend, 'transform', 'translate(0, 0)', mutation);
        centerSingleComparisonLegendParts(pairwiseLegend, mutation);
      }

      const {
        colorRectSize: rectSize,
        lineHeight,
        textXOffset
      } = reflowMetrics;

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

        let cursorX = 0;
        let newY = rectSize / 2;
        let wrapWidth = Math.max(maxWidth, rectSize + textXOffset * 2);
        if (pairwiseLegend) {
          const pairwiseBBox = pairwiseLegend.getBBox();
          const pairwiseWidth = pairwiseBBox.width || 0;
          if (pairwiseWidth > 0) {
            wrapWidth = Math.max(
              wrapWidth - pairwiseWidth - textXOffset,
              rectSize + textXOffset * 2
            );
          }
        }

        entries.forEach((entry) => {
          const textBBox = entry.text.getBBox();
          const entryWidth = rectSize + textXOffset + textBBox.width + textXOffset;

          if (cursorX + entryWidth > wrapWidth && cursorX > 0) {
            cursorX = 0;
            newY += lineHeight;
          }
          const newX = cursorX + textXOffset;

          setAttribute(entry.text, 'transform', `translate(${newX}, ${newY})`, mutation);
          if (entry.rect) {
            setAttribute(
              entry.rect,
              'transform',
              `translate(${newX - textXOffset}, ${newY})`,
              mutation
            );
          }
          entry.newX = newX;
          entry.newY = newY;
          entry.textWidth = textBBox.width;

          cursorX += entryWidth;
        });
        centerHorizontalRows(entries, textXOffset, wrapWidth, mutation);
      } else {
        entries.sort((a, b) => a.y - b.y);

        let newY = rectSize / 2;
        entries.forEach((entry) => {
          setAttribute(entry.text, 'transform', `translate(${textXOffset}, ${newY})`, mutation);
          if (entry.rect) {
            setAttribute(entry.rect, 'transform', `translate(0, ${newY})`, mutation);
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
        const pairwiseAlignmentWidth = getComparisonLegendAlignmentWidth(pairwiseLegend, pairwiseWidth);

        if (layout === 'horizontal') {
          const featureYOffset = pairwiseHeight > featureHeight ? (pairwiseHeight - featureHeight) / 2 : 0;
          const pairwiseYOffset = featureHeight > pairwiseHeight ? (featureHeight - pairwiseHeight) / 2 : 0;
          if (featureYOffset > 0) {
            setAttribute(featureGroup, 'transform', `translate(0, ${featureYOffset})`, mutation);
          }
          const pairwiseX = featureBBox.x + featureWidth + textXOffset;
          setAttribute(
            pairwiseLegend,
            'transform',
            `translate(${pairwiseX}, ${pairwiseYOffset})`,
            mutation
          );
        } else {
          const featureXOffset =
            pairwiseAlignmentWidth > featureWidth ? (pairwiseAlignmentWidth - featureWidth) / 2 : 0;
          setAttribute(featureGroup, 'transform', `translate(${featureXOffset}, 0)`, mutation);
          const pairwiseY = featureHeight + lineHeight / 2;
          setAttribute(pairwiseLegend, 'transform', `translate(0, ${pairwiseY})`, mutation);
        }
      }
    };

    const horizontalDisplay = horizontalLegend.getAttribute('display');
    const verticalDisplay = verticalLegend.getAttribute('display');
    removeAttribute(horizontalLegend, 'display', mutation);
    removeAttribute(verticalLegend, 'display', mutation);

    layoutLegendGroup(horizontalLegend, 'horizontal', horizontalWidth);
    layoutLegendGroup(verticalLegend, 'vertical', horizontalWidth);

    if (horizontalDisplay !== null) {
      setAttribute(horizontalLegend, 'display', horizontalDisplay, mutation);
    } else {
      removeAttribute(horizontalLegend, 'display', mutation);
    }

    if (verticalDisplay !== null) {
      setAttribute(verticalLegend, 'display', verticalDisplay, mutation);
    } else {
      removeAttribute(verticalLegend, 'display', mutation);
    }

  };

  const compactLegendEntries = (svg, mutation = null) => {
    const legendGroup = svg.getElementById('legend');
    const hasDualLegends =
      !!legendGroup?.querySelector('#legend_horizontal') && !!legendGroup?.querySelector('#legend_vertical');
    if (hasDualLegends) {
      reflowDualLegendLayout(svg, mutation);
      return;
    }

    const targetGroups = getAllFeatureLegendGroups(svg);
    if (targetGroups.length === 0) return;

    for (const targetGroup of targetGroups) {
      const comparisonLegend = getComparisonLegendGroup(targetGroup);
      const {
        colorRectSize: rectSize,
        lineHeight,
        textXOffset
      } = getLegendReflowMetrics(svg);

      const texts = Array.from(targetGroup.querySelectorAll('text')).filter((el) => {
        if (!comparisonLegend) return true;
        return !comparisonLegend.contains(el);
      });
      if (texts.length === 0) continue;

      const allRects = Array.from(targetGroup.querySelectorAll('path')).filter((r) => {
        const fill = r.getAttribute('fill');
        if (comparisonLegend && comparisonLegend.contains(r)) return false;
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

        const maxWidth = getHorizontalWrapWidth(svg);

        let cursorX = 0;
        let newY = rectSize / 2;

        entries.forEach((entry) => {
          const textBBox = entry.text.getBBox();
          const entryWidth = rectSize + textXOffset + textBBox.width + textXOffset;

          if (cursorX + entryWidth > maxWidth && cursorX > 0) {
            cursorX = 0;
            newY += lineHeight;
          }
          const newX = cursorX + textXOffset;

          setAttribute(entry.text, 'transform', `translate(${newX}, ${newY})`, mutation);
          entry.newX = newX;
          entry.newY = newY;
          entry.textWidth = textBBox.width;

          if (entry.rect) {
            const expectedRectX = newX - textXOffset;
            setAttribute(entry.rect, 'transform', `translate(${expectedRectX}, ${newY})`, mutation);
          }

          cursorX += entryWidth;
        });
        centerHorizontalRows(entries, textXOffset, maxWidth, mutation);
      } else {
        entries.sort((a, b) => a.y - b.y);

        let newY = rectSize / 2;
        entries.forEach((entry) => {
          setAttribute(entry.text, 'transform', `translate(${textXOffset}, ${newY})`, mutation);

          if (entry.rect) {
            setAttribute(entry.rect, 'transform', `translate(0, ${newY})`, mutation);
          }

          newY += lineHeight;
        });
      }
    }

    updatePairwiseLegendPositions(svg, mutation);
  };

  const reflowSingleLegendLayout = (svg, layout, maxWidthOverride = null, mutation = null) => {
    const legendGroup = svg.getElementById('legend');
    if (!legendGroup) return null;

    const featureLegendGroup = legendGroup.querySelector('#feature_legend') || legendGroup;
    const isRootLegendGroup = featureLegendGroup === legendGroup;
    const pairwiseLegend = getComparisonLegendGroup(legendGroup);

    const textElements = Array.from(featureLegendGroup.querySelectorAll('text')).filter((el) => {
      return !isInsideComparisonLegend(el);
    });
    if (textElements.length === 0) return null;

    const {
      colorRectSize: rectSize,
      lineHeight,
      textXOffset
    } = getLegendReflowMetrics(svg);
    const colorRects = Array.from(featureLegendGroup.querySelectorAll('path')).filter((r) => {
      const fill = r.getAttribute('fill');
      if (!fill || fill === 'none' || fill.startsWith('url(')) return false;
      if (isInsideComparisonLegend(r)) return false;
      return true;
    });
    if (!isRootLegendGroup) {
      setAttribute(featureLegendGroup, 'transform', 'translate(0, 0)', mutation);
    }
    if (pairwiseLegend) {
      setAttribute(pairwiseLegend, 'transform', 'translate(0, 0)', mutation);
      centerSingleComparisonLegendParts(pairwiseLegend, mutation);
    }

    const featureOffset = { x: 0, y: 0 };
    const pairwiseBBox = pairwiseLegend ? pairwiseLegend.getBBox() : null;
    const pairwiseContentOffsetY = pairwiseBBox ? pairwiseBBox.y : 0;
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
      const maxWidth = Number(maxWidthOverride);
      if (!Number.isFinite(maxWidth) || maxWidth <= 0) {
        throw new Error('Horizontal legend reflow requires a positive metadata-owned wrap width.');
      }

      const availableWidth = Math.max(maxWidth, rectSize + textXOffset * 2);
      const reservedOffset = Math.max(isRootLegendGroup ? 0 : featureOffset.x, 0);
      let maxFeatureWidth = availableWidth - reservedOffset;
      if (pairwiseLegend && pairwiseWidth > 0) {
        maxFeatureWidth = availableWidth - reservedOffset - pairwiseWidth - textXOffset;
        if (maxFeatureWidth < rectSize + textXOffset * 2) {
          maxFeatureWidth = availableWidth - reservedOffset;
        }
      }

      entries.sort((a, b) => {
        if (Math.abs(a.y - b.y) < 1) return a.x - b.x;
        return a.y - b.y;
      });

      let cursorX = 0;
      let newY = rectSize / 2;

      entries.forEach((entry) => {
        const textBBox = entry.text.getBBox();
        const entryWidth = rectSize + textXOffset + textBBox.width + textXOffset;

        if (cursorX + entryWidth > maxFeatureWidth && cursorX > 0) {
          cursorX = 0;
          newY += lineHeight;
        }
        const newX = cursorX + textXOffset;

        setAttribute(entry.text, 'transform', `translate(${newX}, ${newY})`, mutation);
        if (entry.rect) {
          setAttribute(
            entry.rect,
            'transform',
            `translate(${newX - textXOffset}, ${newY})`,
            mutation
          );
        }
        entry.newX = newX;
        entry.newY = newY;
        entry.textWidth = textBBox.width;

        cursorX += entryWidth;
      });
      centerHorizontalRows(entries, textXOffset, maxFeatureWidth, mutation);
    } else {
      entries.sort((a, b) => a.y - b.y);
      let newY = rectSize / 2;
      entries.forEach((entry) => {
        const textBBox = entry.text.getBBox();
        setAttribute(entry.text, 'transform', `translate(${textXOffset}, ${newY})`, mutation);
        if (entry.rect) {
          setAttribute(entry.rect, 'transform', `translate(0, ${newY})`, mutation);
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
      const pairwiseAlignmentWidth = getComparisonLegendAlignmentWidth(pairwiseLegend, effectivePairwiseBBox.width);

      let pairwiseX = featureX;
      let pairwiseY = featureY;
      let featureShiftX = 0;
      let featureShiftY = 0;

      if (layout === 'horizontal') {
        const heightDiff = effectivePairwiseBBox.height - featureHeight;
        if (heightDiff > 0) {
          featureShiftY = heightDiff / 2;
        } else if (heightDiff < 0) {
          pairwiseY += (-heightDiff) / 2;
        }
        pairwiseX = featureX + featureWidth + textXOffset;
      } else {
        const widthDiff = pairwiseAlignmentWidth - featureWidth;
        if (widthDiff > 0) {
          featureShiftX = widthDiff / 2;
        }
        pairwiseY = featureY + featureHeight + lineHeight / 2;
      }

      if (!isRootLegendGroup) {
        setAttribute(
          featureLegendGroup,
          'transform',
          `translate(${featureShiftX}, ${featureShiftY})`,
          mutation
        );
      } else if (featureShiftX !== 0 || featureShiftY !== 0) {
        entries.forEach((entry) => {
          const shiftedX = Number(entry.newX) + featureShiftX;
          const shiftedY = Number(entry.newY) + featureShiftY;
          entry.newX = shiftedX;
          entry.newY = shiftedY;
          setAttribute(entry.text, 'transform', `translate(${shiftedX}, ${shiftedY})`, mutation);
          if (entry.rect) {
            setAttribute(
              entry.rect,
              'transform',
              `translate(${shiftedX - textXOffset}, ${shiftedY})`,
              mutation
            );
          }
        });
      }
      const adjustedPairwiseY = pairwiseY - pairwiseContentOffsetY;
      setAttribute(
        pairwiseLegend,
        'transform',
        `translate(${pairwiseX}, ${adjustedPairwiseY})`,
        mutation
      );
    }

    const bbox = legendGroup.getBBox();
    return { legendWidth: bbox.width, legendHeight: bbox.height };
  };

  return {
    compactLegendEntries,
    reflowDualLegendLayout,
    reflowSingleLegendLayout,
    updatePairwiseLegendPositions
  };
};
