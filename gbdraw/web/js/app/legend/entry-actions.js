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

const moveLegendEntryToAnchor = (entryGroup, anchor) => {
  const current = legendEntryAnchor(entryGroup);
  const deltaX = anchor.x - current.x;
  const deltaY = anchor.y - current.y;
  if (Math.abs(deltaX) < 1e-6 && Math.abs(deltaY) < 1e-6) return false;
  entryGroup.querySelectorAll?.('[transform]').forEach((node) => {
    const position = parseTransformXY(node.getAttribute('transform'));
    node.setAttribute('transform', `translate(${position.x + deltaX}, ${position.y + deltaY})`);
  });
  return true;
};

const setLegendEntryColor = (entryGroup, color) => {
  const path = Array.from(entryGroup?.querySelectorAll?.('path') || []).find((candidate) => {
    const fill = candidate.getAttribute('fill');
    return fill && fill !== 'none' && !fill.startsWith('url(');
  });
  if (!path || normalizedColor(path.getAttribute('fill')) === normalizedColor(color)) return false;
  path.setAttribute('fill', color);
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
    || typeof previewRuntime?.runDomEdit !== 'function'
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

  const persistLegendReconciliation = (svg, reason = 'legend-reconcile') => (
    svg
      ? previewRuntime.commitDomEdit({
          reason,
          invalidateIndexes: ['legend'],
          mutate: () => true
        }).changed
      : false
  );

  const applyLegendMutation = (reason, mutate, { commit = true } = {}) => {
    if (!commit) return Boolean(mutate());
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

    try {
      const parser = new DOMParser();

      const widthResponse = await legendHelperOperations.measureText({
        caption,
        fontFamily,
        fontSize
      });
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

        group.appendChild(entryGroup);
      }

      if (shouldReflow) {
        const hasDualLegends =
          !!legendGroup.querySelector('#legend_horizontal') && !!legendGroup.querySelector('#legend_vertical');
        if (hasDualLegends) {
          reflowDualLegendLayout(svg);
        } else {
          updatePairwiseLegendPositions(svg);
        }
        onLegendGeometryChanged();
      }

      if (shouldCommit) {
        persistLegendReconciliation(svg);
      }

      return caption;
    } catch (e) {
      console.error('Failed to add legend entry:', e);
      if (options.throwOnError) throw e;
      return false;
    }
  };

  const updateLegendEntryColorByCaption = (caption, color, { commit = true } = {}) => {
    if (!svgContainer.value) return false;

    const svg = svgContainer.value.querySelector('svg');
    if (!svg) return false;

    const targetGroups = getAllFeatureLegendGroups(svg);
    if (targetGroups.length === 0) return false;

    const updated = applyLegendMutation('legend-fill', () => {
      let changed = false;
      for (const targetGroup of targetGroups) {
        const entryGroup = targetGroup.querySelector(`g[data-legend-key="${CSS.escape(caption)}"]`);
        if (entryGroup) {
          const paths = entryGroup.querySelectorAll('path');
          for (const path of paths) {
            const fill = path.getAttribute('fill');
            if (fill && fill !== 'none' && !fill.startsWith('url(')) {
              if (normalizedColor(fill) === normalizedColor(color)) break;
              path.setAttribute('fill', color);
              changed = true;
              break;
            }
          }
        }
      }
      return changed;
    }, { commit });

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

  const removeLegendEntry = (caption, { commit = true } = {}) => {
    if (!svgContainer.value) return false;

    const svg = svgContainer.value.querySelector('svg');
    if (!svg) return false;

    const targetGroups = getAllFeatureLegendGroups(svg);
    if (targetGroups.length === 0) return false;

    const removed = applyLegendMutation('legend-remove', () => {
      let changed = false;
      targetGroups.forEach((targetGroup, targetIndex) => {
        const entryGroup = targetGroup.querySelector(`g[data-legend-key="${CSS.escape(caption)}"]`);
        if (entryGroup) {
          rememberRetiredEntry(caption, targetGroup, targetIndex, entryGroup);
          entryGroup.remove();
          changed = true;
        }
      });
      if (changed) compactLegendEntries(svg);
      return changed;
    }, { commit });

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

    let changed = false;
    const restoredCaptions = new Set();

    targetGroups.forEach((targetGroup, targetIndex) => {
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
          group.setAttribute('data-legend-key', entry.caption);
          const text = group.querySelector('text');
          if (text && String(text.textContent || '') !== entry.caption) text.textContent = entry.caption;
          changed = true;
        }
        assigned.set(entry.caption, group);
        usedGroups.add(group);
      }

      extras.slice(pairedCount).forEach((group) => {
        const caption = String(group.getAttribute('data-legend-key') || '').trim();
        rememberRetiredEntry(caption, targetGroup, targetIndex, group);
        group.remove();
        changed = true;
      });

      missing.slice(pairedCount).forEach((entry) => {
        const cached = restoredEntryTemplate(entry.caption, targetGroup, targetIndex);
        const fallback = initialGroups[0] || assigned.values().next().value || null;
        const group = (cached || fallback)?.cloneNode?.(true) || null;
        if (!group) return;
        group.setAttribute('data-legend-key', entry.caption);
        if (!cached) group.removeAttribute('data-legend-owner');
        const text = group.querySelector('text');
        if (text) text.textContent = entry.caption;
        targetGroup.appendChild(group);
        assigned.set(entry.caption, group);
        restoredCaptions.add(entry.caption);
        changed = true;
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
        if (setLegendEntryColor(group, entry.color)) changed = true;
        const text = group.querySelector('text');
        if (text && String(text.textContent || '') !== entry.caption) {
          text.textContent = entry.caption;
          changed = true;
        }
        if (slots[index] && moveLegendEntryToAnchor(group, slots[index])) changed = true;
      });
      const currentOrder = directLegendEntryGroups(targetGroup);
      const orderChanged = currentOrder.length !== orderedGroups.length ||
        orderedGroups.some((group, index) => currentOrder[index] !== group);
      if (orderChanged) {
        orderedGroups.forEach((group) => targetGroup.appendChild(group));
        changed = true;
      }
    });

    if (restoreColorState && entryColorStateChanged) {
      legendEntries.value = desiredEntries;
    }
    if (!changed) return restoreColorState && entryColorStateChanged;
    restoredCaptions.forEach((caption) => retiredEntryTemplates.delete(caption));
    const legendGroup = svg.getElementById('legend');
    const hasDualLegends = Boolean(
      legendGroup?.querySelector('#legend_horizontal') && legendGroup?.querySelector('#legend_vertical')
    );
    if (hasDualLegends) reflowDualLegendLayout(svg);
    else updatePairwiseLegendPositions(svg);
    persistLegendReconciliation(svg);
    return true;
  };

  const syncFileLegendEntries = (intents, { previousFileIntents = [] } = {}) => (
    previewRuntime.runDomEdit({
      reason: 'legend-file-sync',
      action: async () => {
    if (!svgContainer.value) {
      return { add: [], update: [], remove: [], unchanged: [] };
    }
    const svg = svgContainer.value.querySelector('svg');
    if (!svg) return { add: [], update: [], remove: [], unchanged: [] };
    const targetGroups = getAllFeatureLegendGroups(svg);
    if (targetGroups.length === 0) return { add: [], update: [], remove: [], unchanged: [] };

    const svgSnapshot = svg.cloneNode(true);
    const resultIndex = previewRuntime.getActiveRuntime?.()?.resultIndex ?? 0;
    const editorSnapshot = legendEntries.value.map((entry) => ({ ...entry }));
    const provenance = new Map(
      previousFileIntents.map((entry) => [String(entry?.caption || '').trim(), normalizedColor(entry?.color)])
    );

    try {
      const desiredByCaption = new Map(intents.map((intent) => [intent.caption, normalizedColor(intent.color)]));
      targetGroups.forEach((group) => {
        Array.from(group.querySelectorAll('g[data-legend-key]')).forEach((entry) => {
          const caption = entry.getAttribute('data-legend-key') || '';
          if (
            !entry.hasAttribute('data-legend-owner') &&
            provenance.get(caption) === legendEntryColor(entry)
          ) {
            entry.setAttribute('data-legend-owner', SPECIFIC_COLOR_FILE_OWNER);
          }
        });
      });

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
      for (const entry of diff.remove) {
        targetGroups.forEach((group) => {
          const target = findLegendEntryGroup(group, entry.caption);
          if (target?.getAttribute('data-legend-owner') === SPECIFIC_COLOR_FILE_OWNER) target.remove();
        });
      }
      for (const entry of diff.update) {
        targetGroups.forEach((group) => {
          const target = findLegendEntryGroup(group, entry.caption);
          if (target?.getAttribute('data-legend-owner') !== SPECIFIC_COLOR_FILE_OWNER) return;
          const path = Array.from(target.querySelectorAll('path')).find((candidate) => {
            const fill = candidate.getAttribute('fill');
            return fill && fill !== 'none' && !fill.startsWith('url(');
          });
          if (path) path.setAttribute('fill', entry.color);
        });
      }
      for (const entry of diff.add) {
        await addLegendEntry(entry.caption, entry.color, {
          owner: SPECIFIC_COLOR_FILE_OWNER,
          conflictPolicy: 'error',
          commit: false,
          reflow: false,
          throwOnError: true
        });
      }

      const legendGroup = svg.getElementById('legend');
      const hasDualLegends =
        !!legendGroup?.querySelector('#legend_horizontal') && !!legendGroup?.querySelector('#legend_vertical');
      if (hasDualLegends) reflowDualLegendLayout(svg);
      else updatePairwiseLegendPositions(svg);
      onLegendGeometryChanged();
      persistLegendReconciliation(svg, 'legend-file-sync');
      extractLegendEntries();
      return diff;
    } catch (error) {
      svg.replaceWith(svgSnapshot);
      previewRuntime.mountResultSvg?.(resultIndex, svgSnapshot);
      legendEntries.value = editorSnapshot;
      throw error;
    }
      }
    })
  );

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

    const targetGroups = getAllFeatureLegendGroups(svg);
    if (targetGroups.length === 0) return false;
    const changed = applyLegendMutation('legend-rename', () => {
      let updated = false;
      for (const targetGroup of targetGroups) {
        const entryGroup = targetGroup.querySelector(`g[data-legend-key="${CSS.escape(oldCaption)}"]`);
        if (!entryGroup) continue;
        entryGroup.setAttribute('data-legend-key', caption);
        const textEl = entryGroup.querySelector('text');
        if (textEl) textEl.textContent = caption;
        updated = true;
      }
      if (!updated) return false;
      const legendGroup = svg.getElementById('legend');
      const hasDualLegends =
        !!legendGroup?.querySelector('#legend_horizontal') && !!legendGroup?.querySelector('#legend_vertical');
      if (hasDualLegends) reflowDualLegendLayout(svg);
      else updatePairwiseLegendPositions(svg);
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
