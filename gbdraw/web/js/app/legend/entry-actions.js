import {
  getAllFeatureLegendGroups,
  getVisibleFeatureLegendGroup,
  parseTransformXY
} from './utils.js';
import { parseCompositionMetadata } from '../legend-layout/composition-actions.js';
import {
  diffLegendIntents,
  SPECIFIC_COLOR_FILE_OWNER
} from '../specific-color-rules.js';

const normalizedColor = (value) => String(value || '').trim().toLowerCase();

const legendEntryColor = (entryGroup) => {
  for (const path of entryGroup?.querySelectorAll?.('path') || []) {
    const fill = path.getAttribute('fill');
    if (fill && fill !== 'none' && !fill.startsWith('url(')) return normalizedColor(fill);
  }
  return '';
};

const findLegendEntryGroup = (targetGroup, caption) => (
  Array.from(targetGroup?.querySelectorAll?.('g[data-legend-key]') || [])
    .find((entry) => entry.getAttribute('data-legend-key') === caption) || null
);

const directLegendEntryGroups = (targetGroup) => {
  const direct = Array.from(targetGroup?.children || []).filter(
    (child) => child.tagName?.toLowerCase() === 'g' && child.hasAttribute?.('data-legend-key')
  );
  return direct.length > 0
    ? direct
    : Array.from(targetGroup?.querySelectorAll?.('g[data-legend-key]') || []);
};

const legendEntryAnchor = (entryGroup) => {
  const groupOffset = parseTransformXY(entryGroup?.getAttribute?.('transform'));
  const target = entryGroup?.querySelector?.('text') || Array.from(
    entryGroup?.querySelectorAll?.('path') || []
  ).find((path) => {
    const fill = path.getAttribute('fill');
    return fill && fill !== 'none' && !fill.startsWith('url(');
  });
  const targetOffset = parseTransformXY(target?.getAttribute?.('transform'));
  return { x: groupOffset.x + targetOffset.x, y: groupOffset.y + targetOffset.y };
};

const moveLegendEntryToAnchor = (entryGroup, anchor, mutation = null) => {
  const current = legendEntryAnchor(entryGroup);
  const deltaX = anchor.x - current.x;
  const deltaY = anchor.y - current.y;
  if (Math.abs(deltaX) < 1e-6 && Math.abs(deltaY) < 1e-6) return false;
  entryGroup.querySelectorAll?.('[transform]').forEach((node) => {
    const position = parseTransformXY(node.getAttribute('transform'));
    const value = `translate(${position.x + deltaX}, ${position.y + deltaY})`;
    if (mutation) mutation.setAttribute(node, 'transform', value);
    else node.setAttribute('transform', value);
  });
  return true;
};

const setLegendEntryColor = (entryGroup, color, mutation = null) => {
  const path = Array.from(entryGroup?.querySelectorAll?.('path') || []).find((candidate) => {
    const fill = candidate.getAttribute('fill');
    return fill && fill !== 'none' && !fill.startsWith('url(');
  });
  if (!path || normalizedColor(path.getAttribute('fill')) === normalizedColor(color)) return false;
  if (mutation) mutation.setAttribute(path, 'fill', color);
  else path.setAttribute('fill', color);
  return true;
};

export const createLegendEntryActions = ({
  state,
  layoutActions,
  previewRuntime = null,
  legendHelperOperations = null
}) => {
  if (
    typeof previewRuntime?.commitDomEdit !== 'function'
    || typeof previewRuntime?.beginDomEditLease !== 'function'
  ) {
    throw new Error('createLegendEntryActions requires the preview runtime edit protocol.');
  }
  if (
    typeof legendHelperOperations?.measureText !== 'function'
    || typeof legendHelperOperations?.generateEntrySvg !== 'function'
  ) {
    throw new Error('createLegendEntryActions requires legend helper operations.');
  }
  const {
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
    manualSpecificRules
  } = state;

  const { updatePairwiseLegendPositions, reflowDualLegendLayout, compactLegendEntries } = layoutActions;
  let legendGeometryChangedHandler = null;
  const retiredEntryTemplates = new Map();

  const targetGroupKey = (targetGroup, index) => (
    String(targetGroup?.id || targetGroup?.parentElement?.id || `legend-target-${index}`)
  );

  const rememberRetiredEntry = (caption, targetGroup, targetIndex, entryGroup) => {
    const normalizedCaption = String(caption || '').trim();
    if (!normalizedCaption || !entryGroup?.cloneNode) return;
    if (!retiredEntryTemplates.has(normalizedCaption)) {
      retiredEntryTemplates.set(normalizedCaption, new Map());
    }
    retiredEntryTemplates.get(normalizedCaption).set(
      targetGroupKey(targetGroup, targetIndex),
      entryGroup.cloneNode(true)
    );
  };

  const restoredEntryTemplate = (caption, targetGroup, targetIndex) => {
    const templates = retiredEntryTemplates.get(String(caption || '').trim());
    if (!templates) return null;
    return templates.get(targetGroupKey(targetGroup, targetIndex)) || templates.values().next().value || null;
  };

  const applyLegendMutation = (reason, mutate) => {
    return previewRuntime.commitDomEdit({
      reason,
      invalidateIndexes: ['legend'],
      mutate
    }).changed;
  };

  const setLegendGeometryChangedHandler = (handler) => {
    legendGeometryChangedHandler = typeof handler === 'function' ? handler : null;
  };

  const onLegendGeometryChanged = () => legendGeometryChangedHandler?.();

  const addLegendEntry = async (caption, color, options = {}) => {
    const owner = String(options.owner || '').trim();
    const conflictPolicy = options.conflictPolicy || 'suffix';
    const shouldCommit = options.commit !== false;
    const shouldReflow = options.reflow !== false;
    console.log(`addLegendEntry called with caption="${caption}", color="${color}"`);
    if (!svgContainer.value) return false;

    const svg = svgContainer.value.querySelector('svg');
    if (!svg) return false;
    const composition = parseCompositionMetadata(svg);
    const reflowMetrics = composition.legendReflow;
    if (!reflowMetrics) {
      throw new Error('This diagram has no legend reflow metadata. Regenerate it before editing the legend.');
    }

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

    let finalCaption = caption.trim();

    const keyedEntry = findLegendEntryGroup(targetGroup, finalCaption);
    if (keyedEntry) {
      if (legendEntryColor(keyedEntry) === normalizedColor(color)) return finalCaption;
      if (conflictPolicy === 'error') {
        throw new Error(`Legend entry "${finalCaption}" already exists with a different color.`);
      }
    } else {
      const legacyText = Array.from(targetGroup.querySelectorAll('text'))
        .find((text) => text.textContent?.trim() === finalCaption);
      if (legacyText) {
        const textPosition = parseTransformXY(legacyText.getAttribute('transform'));
        const legacyPath = Array.from(targetGroup.querySelectorAll('path')).find((path) => {
          const fill = path.getAttribute('fill');
          if (!fill || fill === 'none' || fill.startsWith('url(')) return false;
          const pathPosition = parseTransformXY(path.getAttribute('transform'));
          return Math.abs(pathPosition.y - textPosition.y) < 2 && pathPosition.x < textPosition.x;
        });
        if (legacyPath && normalizedColor(legacyPath.getAttribute('fill')) === normalizedColor(color)) {
          return finalCaption;
        }
        if (legacyPath && conflictPolicy === 'error') {
          throw new Error(`Legend entry "${finalCaption}" already exists with a different color.`);
        }
      }
    }

    const baseCaption = finalCaption.replace(/\s*\(\d+\)$/, '');
    let counter = 1;
    const existingCaptions = new Set();
    targetGroup.querySelectorAll('text').forEach((t) => existingCaptions.add(t.textContent?.trim()));

    while (existingCaptions.has(finalCaption)) {
      finalCaption = `${baseCaption} (${counter})`;
      counter++;
    }

    caption = finalCaption;

    const {
      colorRectSize: rectSize,
      lineHeight: lineMargin,
      textXOffset: xMargin
    } = reflowMetrics;
    const firstColorRect = targetGroup.querySelector('path[fill]:not([fill="none"]):not([fill^="url("])');

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

    let maxY = -lineMargin;
    const textElements = targetGroup.querySelectorAll('text');
    textElements.forEach((el) => {
      const transform = el.getAttribute('transform');
      if (transform) {
        const { y } = parseTransformXY(transform);
        if (y > maxY) maxY = y;
      }
    });

    if (textElements.length === 0) {
      const colorRects = targetGroup.querySelectorAll('path');
      colorRects.forEach((el) => {
        const fill = el.getAttribute('fill');
        if (fill && fill !== 'none' && !fill.startsWith('url(')) {
          const transform = el.getAttribute('transform');
          if (transform) {
            const { y } = parseTransformXY(transform);
            if (y > maxY) maxY = y;
          }
        }
      });
    }

    const newY = maxY + lineMargin;
    const sharedLease = options.lease || null;
    const editLease = sharedLease || previewRuntime.beginDomEditLease({
      reason: 'legend-add',
      invalidateIndexes: ['legend']
    });
    if (!editLease?.current) return false;
    if (!sharedLease && !shouldCommit) {
      editLease.cancel();
      throw new Error('A non-committing Legend add requires an explicit shared edit lease.');
    }

    try {
      const parser = new DOMParser();
      const preparedEntries = [];

      const widthResponse = await legendHelperOperations.measureText({
        caption,
        fontFamily,
        fontSize
      });
      if (!editLease.current) throw new Error('The mounted Legend changed while measuring its entry.');
      if (widthResponse.result?.error) throw new Error(widthResponse.result.error);
      const measuredWidth = Number(widthResponse.result?.width);
      if (!Number.isFinite(measuredWidth) || measuredWidth < 0) {
        throw new Error('Python returned an invalid legend text width.');
      }
      const entryWidth = rectSize + xMargin + measuredWidth + xMargin;
      const canvasWidth = composition.primary.finalBounds.width;

      for (const group of allTargetGroups) {
        const parentId = group.parentElement?.id || '';
        const isHorizontalGroup = parentId === 'legend_horizontal';

        let newX = 0,
          newY = 0;

        if (isHorizontalGroup) {
          const featureLegendMaxWidth = canvasWidth;

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

        const entryResponse = await legendHelperOperations.generateEntrySvg({
            caption,
            color,
            yOffset: newY,
            rectSize,
            fontSize,
            fontFamily,
            xOffset: newX,
            strokeColor,
            strokeWidth
          });
        if (!editLease.current) throw new Error('The mounted Legend changed while generating its entry.');
        const result = entryResponse.result;
        if (result?.error) throw new Error(result.error);

        const entryGroup = document.createElementNS('http://www.w3.org/2000/svg', 'g');
        entryGroup.setAttribute('data-legend-key', caption);
        if (owner) entryGroup.setAttribute('data-legend-owner', owner);

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

        preparedEntries.push({ group, entryGroup });
      }

      const applied = editLease.mutate(({ svg: targetSvg, mutation }) => {
        if (targetSvg !== svg) return false;
        const currentTargets = getAllFeatureLegendGroups(targetSvg);
        if (
          currentTargets.length !== allTargetGroups.length
          || currentTargets.some((group, index) => group !== allTargetGroups[index])
        ) return false;
        preparedEntries.forEach(({ group, entryGroup }) => mutation.appendChild(group, entryGroup));
        if (shouldReflow) {
          const currentLegendGroup = targetSvg.getElementById('legend');
          const hasDualLegends = Boolean(
            currentLegendGroup?.querySelector('#legend_horizontal')
            && currentLegendGroup?.querySelector('#legend_vertical')
          );
          if (hasDualLegends) reflowDualLegendLayout(targetSvg);
          else updatePairwiseLegendPositions(targetSvg);
        }
        return preparedEntries.length;
      });
      if (!applied) {
        if (!sharedLease) editLease.cancel();
        return false;
      }
      if (!sharedLease) editLease.commit();
      if (shouldReflow) onLegendGeometryChanged();

      return caption;
    } catch (e) {
      if (!sharedLease) editLease.cancel();
      console.error('Failed to add legend entry:', e);
      if (options.throwOnError) throw e;
      return false;
    }
  };

  const updateLegendEntryColorByCaption = (caption, color) => {
    if (!svgContainer.value) return false;

    const svg = svgContainer.value.querySelector('svg');
    if (!svg) return false;

    const updated = applyLegendMutation('legend-fill', ({ svg: targetSvg, mutation }) => {
      let changed = false;
      for (const targetGroup of getAllFeatureLegendGroups(targetSvg)) {
        const entryGroup = targetGroup.querySelector(`g[data-legend-key="${CSS.escape(caption)}"]`);
        if (entryGroup) {
          const paths = entryGroup.querySelectorAll('path');
          for (const path of paths) {
            const fill = path.getAttribute('fill');
            if (fill && fill !== 'none' && !fill.startsWith('url(')) {
              if (normalizedColor(fill) === normalizedColor(color)) break;
              mutation.setAttribute(path, 'fill', color);
              changed = true;
              break;
            }
          }
        }
      }
      return changed;
    });

    if (updated) {
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

    const removed = applyLegendMutation('legend-remove', ({ svg: targetSvg, mutation }) => {
      let changed = false;
      getAllFeatureLegendGroups(targetSvg).forEach((targetGroup, targetIndex) => {
        const entryGroup = targetGroup.querySelector(`g[data-legend-key="${CSS.escape(caption)}"]`);
        if (entryGroup) {
          rememberRetiredEntry(caption, targetGroup, targetIndex, entryGroup);
          mutation.removeNode(entryGroup);
          changed = true;
        }
      });
      if (changed) compactLegendEntries(targetSvg);
      return changed;
    });

    if (removed) {
      onLegendGeometryChanged();
      console.log(`Removed legend entry: "${caption}"`);
    }

    return removed;
  };

  const reconcileLegendEntries = ({ restoreColorState = false } = {}) => {
    const svg = svgContainer.value?.querySelector?.('svg');
    if (!svg) return false;
    const targetGroups = getAllFeatureLegendGroups(svg);
    if (targetGroups.length === 0) return false;

    const desiredEntries = [];
    let entryColorStateChanged = false;
    const seenCaptions = new Set();
    (Array.isArray(legendEntries.value) ? legendEntries.value : []).forEach((entry) => {
      const caption = String(entry?.caption || '').trim();
      if (!caption || seenCaptions.has(caption)) return;
      seenCaptions.add(caption);
      const entryColor = String(entry?.color || '#cccccc');
      const color = restoreColorState
        ? String(
            legendColorOverrides[caption]
            || originalLegendColors.value?.[caption]
            || entryColor
          )
        : entryColor;
      if (normalizedColor(entry?.color) !== normalizedColor(color)) {
        entryColorStateChanged = true;
      }
      desiredEntries.push({ ...entry, caption, color });
    });

    const restoredCaptions = new Set();
    const changed = applyLegendMutation('legend-reconcile', ({ svg: targetSvg, mutation }) => {
      let domChanged = false;
      getAllFeatureLegendGroups(targetSvg).forEach((targetGroup, targetIndex) => {
      const initialGroups = directLegendEntryGroups(targetGroup);
      const groupsByCaption = new Map(
        initialGroups.map((group) => [String(group.getAttribute('data-legend-key') || '').trim(), group])
      );
      const assigned = new Map();
      const usedGroups = new Set();

      desiredEntries.forEach((entry) => {
        const exact = groupsByCaption.get(entry.caption);
        if (!exact || usedGroups.has(exact)) return;
        assigned.set(entry.caption, exact);
        usedGroups.add(exact);
      });

      const missing = desiredEntries.filter((entry) => !assigned.has(entry.caption));
      const extras = initialGroups.filter((group) => !usedGroups.has(group));
      const pairedCount = Math.min(missing.length, extras.length);

      for (let index = 0; index < pairedCount; index += 1) {
        const entry = missing[index];
        const group = extras[index];
        const previousCaption = String(group.getAttribute('data-legend-key') || '').trim();
        if (previousCaption !== entry.caption) {
          mutation.setAttribute(group, 'data-legend-key', entry.caption);
          const text = group.querySelector('text');
          if (text && String(text.textContent || '') !== entry.caption) {
            mutation.setTextContent(text, entry.caption);
          }
          domChanged = true;
        }
        assigned.set(entry.caption, group);
        usedGroups.add(group);
      }

      extras.slice(pairedCount).forEach((group) => {
        const caption = String(group.getAttribute('data-legend-key') || '').trim();
        rememberRetiredEntry(caption, targetGroup, targetIndex, group);
        mutation.removeNode(group);
        domChanged = true;
      });

      missing.slice(pairedCount).forEach((entry) => {
        const cached = restoredEntryTemplate(entry.caption, targetGroup, targetIndex);
        const fallback = initialGroups[0] || assigned.values().next().value || null;
        const group = (cached || fallback)?.cloneNode?.(true) || null;
        if (!group) return;
        mutation.setAttribute(group, 'data-legend-key', entry.caption);
        if (!cached) mutation.removeAttribute(group, 'data-legend-owner');
        const text = group.querySelector('text');
        if (text) mutation.setTextContent(text, entry.caption);
        mutation.appendChild(targetGroup, group);
        assigned.set(entry.caption, group);
        restoredCaptions.add(entry.caption);
        domChanged = true;
      });

      const orderedGroups = desiredEntries
        .map((entry) => assigned.get(entry.caption))
        .filter(Boolean);
      const slots = orderedGroups
        .map((group) => legendEntryAnchor(group))
        .sort((left, right) => {
          const yDelta = left.y - right.y;
          return Math.abs(yDelta) < 1 ? left.x - right.x : yDelta;
        });

      desiredEntries.forEach((entry, index) => {
        const group = assigned.get(entry.caption);
        if (!group) return;
        if (setLegendEntryColor(group, entry.color, mutation)) domChanged = true;
        const text = group.querySelector('text');
        if (text && String(text.textContent || '') !== entry.caption) {
          mutation.setTextContent(text, entry.caption);
          domChanged = true;
        }
        if (slots[index] && moveLegendEntryToAnchor(group, slots[index], mutation)) domChanged = true;
      });
      const currentOrder = directLegendEntryGroups(targetGroup);
      const orderChanged = currentOrder.length !== orderedGroups.length ||
        orderedGroups.some((group, index) => currentOrder[index] !== group);
      if (orderChanged) {
        orderedGroups.forEach((group) => mutation.appendChild(targetGroup, group));
        domChanged = true;
      }
      });
      if (!domChanged) return false;
      const legendGroup = targetSvg.getElementById('legend');
      const hasDualLegends = Boolean(
        legendGroup?.querySelector('#legend_horizontal') && legendGroup?.querySelector('#legend_vertical')
      );
      if (hasDualLegends) reflowDualLegendLayout(targetSvg);
      else updatePairwiseLegendPositions(targetSvg);
      return true;
    });

    if (restoreColorState && entryColorStateChanged) {
      legendEntries.value = desiredEntries;
    }
    if (!changed) return restoreColorState && entryColorStateChanged;
    restoredCaptions.forEach((caption) => retiredEntryTemplates.delete(caption));
    onLegendGeometryChanged();
    return true;
  };

  const syncFileLegendEntries = async (intents, { previousFileIntents = [] } = {}) => {
    const emptyDiff = { add: [], update: [], remove: [], unchanged: [] };
    if (!svgContainer.value) return emptyDiff;
    const svg = svgContainer.value.querySelector('svg');
    if (!svg || getAllFeatureLegendGroups(svg).length === 0) return emptyDiff;
    const lease = previewRuntime.beginDomEditLease({
      reason: 'legend-file-sync',
      invalidateIndexes: ['legend']
    });
    if (!lease) return emptyDiff;
    const editorSnapshot = legendEntries.value.map((entry) => ({ ...entry }));
    const provenance = new Map(
      previousFileIntents.map((entry) => [String(entry?.caption || '').trim(), normalizedColor(entry?.color)])
    );

    try {
      const desiredByCaption = new Map(
        intents.map((intent) => [intent.caption, normalizedColor(intent.color)])
      );
      lease.mutate(({ svg: targetSvg, mutation }) => {
        let changed = false;
        getAllFeatureLegendGroups(targetSvg).forEach((group) => {
          Array.from(group.querySelectorAll('g[data-legend-key]')).forEach((entry) => {
            const caption = entry.getAttribute('data-legend-key') || '';
            if (
              !entry.hasAttribute('data-legend-owner')
              && provenance.get(caption) === legendEntryColor(entry)
            ) {
              changed = mutation.setAttribute(
                entry,
                'data-legend-owner',
                SPECIFIC_COLOR_FILE_OWNER
              ) || changed;
            }
          });
        });
        return changed;
      });

      const targetGroups = getAllFeatureLegendGroups(svg);
      const primaryEntries = Array.from(targetGroups[0].querySelectorAll('g[data-legend-key]'));
      const ownedEntries = primaryEntries
        .filter((entry) => entry.getAttribute('data-legend-owner') === SPECIFIC_COLOR_FILE_OWNER)
        .map((entry) => ({
          caption: entry.getAttribute('data-legend-key') || '',
          color: legendEntryColor(entry)
        }));
      const reusableEntries = primaryEntries
        .filter((entry) => entry.getAttribute('data-legend-owner') !== SPECIFIC_COLOR_FILE_OWNER)
        .map((entry) => ({
          caption: entry.getAttribute('data-legend-key') || '',
          color: legendEntryColor(entry)
        }))
        .filter((entry) => desiredByCaption.get(entry.caption) === entry.color);

      for (const intent of intents) {
        for (const group of targetGroups) {
          const existing = findLegendEntryGroup(group, intent.caption);
          if (!existing || existing.getAttribute('data-legend-owner') === SPECIFIC_COLOR_FILE_OWNER) continue;
          if (legendEntryColor(existing) !== normalizedColor(intent.color)) {
            throw new Error(`Legend entry "${intent.caption}" already exists with a different color.`);
          }
        }
      }

      const diff = diffLegendIntents([...ownedEntries, ...reusableEntries], intents);
      lease.mutate(({ svg: targetSvg, mutation }) => {
        let changed = false;
        const currentTargets = getAllFeatureLegendGroups(targetSvg);
        diff.remove.forEach((entry) => {
          currentTargets.forEach((group) => {
            const target = findLegendEntryGroup(group, entry.caption);
            if (target?.getAttribute('data-legend-owner') === SPECIFIC_COLOR_FILE_OWNER) {
              changed = mutation.removeNode(target) || changed;
            }
          });
        });
        diff.update.forEach((entry) => {
          currentTargets.forEach((group) => {
            const target = findLegendEntryGroup(group, entry.caption);
            if (target?.getAttribute('data-legend-owner') !== SPECIFIC_COLOR_FILE_OWNER) return;
            const path = Array.from(target.querySelectorAll('path')).find((candidate) => {
              const fill = candidate.getAttribute('fill');
              return fill && fill !== 'none' && !fill.startsWith('url(');
            });
            if (path) changed = mutation.setAttribute(path, 'fill', entry.color) || changed;
          });
        });
        return changed;
      });
      for (const entry of diff.add) {
        const added = await addLegendEntry(entry.caption, entry.color, {
          owner: SPECIFIC_COLOR_FILE_OWNER,
          conflictPolicy: 'error',
          commit: false,
          lease,
          reflow: false,
          throwOnError: true
        });
        if (!added) throw new Error(`Failed to add Legend entry "${entry.caption}".`);
      }
      if (!lease.current) throw new Error('The mounted Legend changed during file synchronization.');
      lease.mutate(({ svg: targetSvg }) => {
        if (diff.add.length + diff.update.length + diff.remove.length === 0) return false;
        const legendGroup = targetSvg.getElementById('legend');
        const hasDualLegends = Boolean(
          legendGroup?.querySelector('#legend_horizontal')
          && legendGroup?.querySelector('#legend_vertical')
        );
        if (hasDualLegends) reflowDualLegendLayout(targetSvg);
        else updatePairwiseLegendPositions(targetSvg);
        return true;
      });
      lease.commit();
      onLegendGeometryChanged();
      extractLegendEntries();
      return diff;
    } catch (error) {
      lease.cancel();
      legendEntries.value = editorSnapshot;
      throw error;
    }
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

    const entryGroups = targetGroup.querySelectorAll('g[data-legend-key]');

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
      const groupTransform = parseTransformXY(entryGroup.getAttribute('transform'));
      const textEl = entryGroup.querySelector('text');
      if (textEl) {
        const textTransform = parseTransformXY(textEl.getAttribute('transform'));
        xPos = groupTransform.x + textTransform.x;
        yPos = groupTransform.y + textTransform.y;
      } else if (groupTransform.x !== 0 || groupTransform.y !== 0) {
        xPos = groupTransform.x;
        yPos = groupTransform.y;
      }

      const existingEntry = legendEntries.value.find((entry) => (
        entry.caption === caption
        && normalizedColor(entry.color) === normalizedColor(color)
      ));
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

    const visuallySortedEntries = [...entries].sort((a, b) => {
      const yDelta = a.yPos - b.yPos;
      if (Math.abs(yDelta) < 1) {
        const xDelta = a.xPos - b.xPos;
        if (Math.abs(xDelta) < 1) {
          return a.caption.localeCompare(b.caption, undefined, { sensitivity: 'base' });
        }
        return xDelta;
      }
      return yDelta;
    });

    legendEntries.value = visuallySortedEntries;

    if (originalLegendOrder.value.length === 0 && visuallySortedEntries.length > 0) {
      originalLegendOrder.value = visuallySortedEntries.map((e) => e.caption);
    }

    if (Object.keys(originalLegendColors.value).length === 0 && visuallySortedEntries.length > 0) {
      visuallySortedEntries.forEach((entry) => {
        originalLegendColors.value[entry.caption] = entry.color;
      });
    }
  };

  const updateLegendEntryColor = (idx, newColor) => {
    if (!svgContainer.value) return false;
    const svg = svgContainer.value.querySelector('svg');
    if (!svg) return false;

    const entry = legendEntries.value[idx];
    if (!entry) return false;

    const caption = entry.caption;
    if (!caption) return false;

    const changed = updateLegendEntryColorByCaption(caption, newColor);
    const stateChanged = normalizedColor(entry.color) !== normalizedColor(newColor);
    if (!changed && !stateChanged) return false;
    entry.color = newColor;
    legendColorOverrides[caption] = newColor;
    return true;
  };

  const updateLegendEntryCaption = (idx, newCaption) => {
    if (!svgContainer.value) return false;
    const svg = svgContainer.value.querySelector('svg');
    if (!svg) return false;

    const entry = legendEntries.value[idx];
    if (!entry) return false;

    const oldCaption = entry.caption;
    const caption = newCaption.trim();
    if (!caption || caption === oldCaption) return false;

    if (getAllFeatureLegendGroups(svg).length === 0) return false;
    const changed = applyLegendMutation('legend-rename', ({ svg: targetSvg, mutation }) => {
      let updated = false;
      for (const targetGroup of getAllFeatureLegendGroups(targetSvg)) {
        const entryGroup = targetGroup.querySelector(`g[data-legend-key="${CSS.escape(oldCaption)}"]`);
        if (!entryGroup) continue;
        mutation.setAttribute(entryGroup, 'data-legend-key', caption);
        const textEl = entryGroup.querySelector('text');
        if (textEl) mutation.setTextContent(textEl, caption);
        updated = true;
      }
      if (!updated) return false;
      const legendGroup = targetSvg.getElementById('legend');
      const hasDualLegends =
        !!legendGroup?.querySelector('#legend_horizontal') && !!legendGroup?.querySelector('#legend_vertical');
      if (hasDualLegends) reflowDualLegendLayout(targetSvg);
      else updatePairwiseLegendPositions(targetSvg);
      return true;
    });
    if (!changed) return false;

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
    onLegendGeometryChanged();
    return true;
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
    onLegendGeometryChanged,
    removeLegendEntry,
    reconcileLegendEntries,
    restoreDeletedLegendEntries,
    setLegendGeometryChangedHandler,
    syncFileLegendEntries,
    updateLegendEntryCaption,
    updateLegendEntryColor,
    updateLegendEntryColorByCaption
  };
};
