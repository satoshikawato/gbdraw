import {
  getAllFeatureLegendGroups,
  getLegendChildById,
  isCurrentLegendHorizontal,
  parseTransform,
  parseTransformXY
} from './utils.js';

export const createLegendLayoutActions = ({ state }) => {
  const { mode, form, linearBaseConfig } = state;

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
    const pairwiseLegend = getLegendChildById(verticalLegend, 'pairwise_legend');

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
      svg.setAttribute('data-vertical-viewbox', `${vbX} ${vbY} ${vbW} ${newVbH}`);
      linearBaseConfig.value.verticalViewBox = { x: vbX, y: vbY, w: vbW, h: newVbH };
    }
  };

  const expandCanvasForHorizontalLegend = (svg) => {
    if (mode.value !== 'linear') return;
    if (!isCurrentLegendHorizontal(svg)) return;

    const legendGroup = svg.getElementById('legend');
    if (!legendGroup) return;

    const viewBox = svg.getAttribute('viewBox');
    if (!viewBox) return;
    const parts = viewBox.split(/\s+/).map(parseFloat);
    if (parts.length !== 4) return;
    const [vbX, vbY, vbW, vbH] = parts;

    const legendOffset = parseTransform(legendGroup.getAttribute('transform'));
    const legendBox = legendGroup.getBBox();
    const legendBottomEdge = legendOffset.y + legendBox.y + legendBox.height;
    const bottomPadding = 20;

    if (legendBottomEdge + bottomPadding > vbY + vbH) {
      const newVbH = legendBottomEdge + bottomPadding - vbY;
      console.log(`Expanding canvas for horizontal legend: ${vbH.toFixed(1)} -> ${newVbH.toFixed(1)}`);
      svg.setAttribute('viewBox', `${vbX} ${vbY} ${vbW} ${newVbH}`);
      svg.setAttribute('data-horizontal-viewbox', `${vbX} ${vbY} ${vbW} ${newVbH}`);
      linearBaseConfig.value.horizontalViewBox = { x: vbX, y: vbY, w: vbW, h: newVbH };
    }
  };

  const updatePairwiseLegendPositions = (svg) => {
    const legendGroup = svg.getElementById('legend');
    if (!legendGroup) return;

    const horizontalLegend = legendGroup.querySelector('#legend_horizontal');
    const verticalLegend = legendGroup.querySelector('#legend_vertical');
    const hasDualLegends = !!(horizontalLegend && verticalLegend);

    if (!hasDualLegends) {
      const layout = form.legend === 'top' || form.legend === 'bottom' ? 'horizontal' : 'vertical';
      reflowSingleLegendLayout(svg, layout);
      return;
    }

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

    if (verticalLegend) {
      const vFeatureLegend = verticalLegend.querySelector('#feature_legend_v');
      const vPairwiseLegend = getLegendChildById(verticalLegend, 'pairwise_legend');

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

    if (horizontalLegend) {
      const hFeatureLegend = horizontalLegend.querySelector('#feature_legend_h');
      const hPairwiseLegend = getLegendChildById(horizontalLegend, 'pairwise_legend');

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
    expandCanvasForHorizontalLegend(svg);
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

      const pairwiseLegend = getLegendChildById(legend, 'pairwise_legend');

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
    expandCanvasForHorizontalLegend(svg);
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
            const hPairwiseLegend = getLegendChildById(horizontalLegend, 'pairwise_legend');
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
    const isRootLegendGroup = featureLegendGroup === legendGroup;
    const pairwiseLegend = legendGroup.querySelector('#pairwise_legend');
    const pairwiseTransform = pairwiseLegend
      ? parseTransform(pairwiseLegend.getAttribute('transform'))
      : { x: 0, y: 0 };

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
    const featureOffset = isRootLegendGroup
      ? { x: 0, y: 0 }
      : parseTransformXY(featureLegendGroup.getAttribute('transform'));
    const pairwiseBBox = pairwiseLegend ? pairwiseLegend.getBBox() : null;
    const pairwiseContentOffsetY = pairwiseBBox ? pairwiseBBox.y - pairwiseTransform.y : 0;
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

      const reservedOffset = Math.max(isRootLegendGroup ? 0 : featureOffset.x, 0);
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

      const adjustedPairwiseY = pairwiseY - pairwiseContentOffsetY;
      pairwiseLegend.setAttribute('transform', `translate(${pairwiseX}, ${adjustedPairwiseY})`);
    }

    const bbox = legendGroup.getBBox();
    return { legendWidth: bbox.width, legendHeight: bbox.height };
  };

  return {
    compactLegendEntries,
    expandCanvasForHorizontalLegend,
    expandCanvasForVerticalLegend,
    reflowDualLegendLayout,
    reflowSingleLegendLayout,
    updatePairwiseLegendPositions
  };
};
