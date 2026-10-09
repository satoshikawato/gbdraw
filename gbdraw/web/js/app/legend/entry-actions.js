// @ts-check
/** @import { DrawingState } from '../../state.js' */
import { resolveColorToHex, toNativeColorInputValue } from '../../utils/color-utils.js';
import {
  defaultLegendCaptionOrder,
  getAllFeatureLegendGroups,
  getLegendEntrySwatch,
  getVisibleFeatureLegendGroup,
  isLegendOrderEdited,
  mountedLegendRowFeatureIds,
  orderLegendEntries,
  parseTransformXY,
  setsFeatureStroke
} from '../../services/legend-svg.js';
import { parseCompositionMetadata } from '../legend-layout/composition-actions.js';
import {
  diffLegendIntents,
  legendRowRules,
  SPECIFIC_COLOR_FILE_OWNER
} from '../../services/specific-color-rules.js';
import { getFeatureElementIndex, getFeatureFillElements, isAutoFeatureUnderlay } from '../../services/feature-dom.js';
import { featureOverrideKey } from '../../services/feature-override-identity.js';

const normalizedColor = (value) => {
  const resolved = String(resolveColorToHex(String(value || '').trim()) || value || '').trim().toLowerCase();
  return resolved.startsWith('#') ? toNativeColorInputValue(resolved) : resolved;
};

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

const setLegendEntryColor = (entryGroup, color) => {
  const path = Array.from(entryGroup?.querySelectorAll?.('path') || []).find((candidate) => {
    const fill = candidate.getAttribute('fill');
    return fill && fill !== 'none' && !fill.startsWith('url(');
  });
  if (!path || normalizedColor(path.getAttribute('fill')) === normalizedColor(color)) return false;
  path.setAttribute('fill', color);
  return true;
};

/**
 * @typedef {object} LegendEntryActionsOptions
 * @property {Record<string, any>} state App state (state.js; not yet typed).
 * @property {((reason: string) => boolean) | null} [commitActiveResultEdit]
 *   The preview owner's commit of an edit to the displayed Result (R1, R13).
 * @property {(() => string | undefined) | null} [readActiveResultIdentity]
 *   The preview owner's runtime identity of the mounted Result.
 */

/** @param {LegendEntryActionsOptions} options */
export const createLegendEntryActions = ({
  state,
  commitActiveResultEdit = null,
  readActiveResultIdentity = null
}) => {
  const {
    results,
    svgContainer,
    originalLegendOrder,
    originalLegendColors,
    newLegendCaption,
    newLegendColor,
    originalSvgStroke
  } = state;

  /** @type {((options?: { commit?: boolean }) => unknown) | null} */
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

  // A batch Result keeps the Legend entries that a deletion removed from it,
  // so an Undo or restore of that deletion returns them when the Result is
  // displayed again (D-07). Keyed by the mounted Result's runtime identity.
  const retiredEntriesByResult = new Map();

  // The generated inventory of Legend captions in Python's order, per Result
  // (OV-47): each Result has its own order, so the displayed Result's order is
  // `originalLegendOrder`, and the others wait here under their runtime
  // identity until they are displayed. `replayedInventoryResults` are the
  // Results a draw replayed an edited order onto: their drawing no longer shows
  // Python's order.
  /** @type {Map<string, string[]>} */
  const inventoryByResult = new Map();
  /** @type {Set<string>} */
  const replayedInventoryResults = new Set();

  const rememberResultEntry = (
    caption, targetGroup, targetIndex, entryGroup,
    identity = readActiveResultIdentity?.()
  ) => {
    if (!identity || !entryGroup?.cloneNode) return;
    const entries = retiredEntriesByResult.get(identity) || new Map();
    const targetKey = targetGroupKey(targetGroup, targetIndex);
    entries.set(`${targetKey}\u0000${caption}`, { caption, targetKey, node: entryGroup.cloneNode(true) });
    retiredEntriesByResult.set(identity, entries);
  };

  /** @param {{ resultIdentity?: string, deletedCaptions?: string[] }} [options] */
  const hasRetiredResultLegend = ({ resultIdentity, deletedCaptions = [] } = {}) => {
    const deleted = new Set(deletedCaptions.map((caption) => String(caption || '').trim()));
    return [...(retiredEntriesByResult.get(resultIdentity)?.values() || [])]
      .some((entry) => !deleted.has(entry.caption));
  };

  // Before the displayed Result receives the diagram-wide Legend operations:
  // return its entries whose deletion is no longer intended, and keep the
  // entries that the operations delete. Returns whether its Legend changes.
  /**
   * @param {SVGSVGElement} svg
   * @param {{ resultIdentity?: string, liveResultIdentities?: string[], deletedCaptions?: string[] }} [options]
   */
  const prepareDisplayedResultLegend = (svg, {
    resultIdentity: identity,
    liveResultIdentities = [],
    deletedCaptions = []
  } = {}) => {
    const liveIdentities = new Set(liveResultIdentities);
    [...retiredEntriesByResult.keys()].forEach((key) => {
      if (!liveIdentities.has(key)) retiredEntriesByResult.delete(key);
    });
    const deleted = new Set(deletedCaptions.map((caption) => String(caption || '').trim()));
    const retired = retiredEntriesByResult.get(identity);
    let changed = false;
    getAllFeatureLegendGroups(svg).forEach((targetGroup, targetIndex) => {
      const targetKey = targetGroupKey(targetGroup, targetIndex);
      retired?.forEach((entry, key) => {
        if (entry.targetKey !== targetKey || deleted.has(entry.caption)) return;
        if (!findLegendEntryGroup(targetGroup, entry.caption)) {
          targetGroup.appendChild(entry.node.cloneNode(true));
          changed = true;
        }
        retired.delete(key);
      });
      deleted.forEach((caption) => {
        const entryGroup = findLegendEntryGroup(targetGroup, caption);
        if (!entryGroup) return;
        rememberResultEntry(caption, targetGroup, targetIndex, entryGroup, identity);
        changed = true;
      });
    });
    return changed;
  };

  const restoredEntryTemplate = (caption, targetGroup, targetIndex) => {
    const templates = retiredEntryTemplates.get(String(caption || '').trim());
    if (!templates) return null;
    return templates.get(targetGroupKey(targetGroup, targetIndex)) || templates.values().next().value || null;
  };

  const persistLegendReconciliation = () => commitActiveResultEdit?.('history-legend-reconcile');

  const setLegendGeometryChangedHandler = (handler) => {
    legendGeometryChangedHandler = typeof handler === 'function' ? handler : null;
  };

  /** @param {{ commit?: boolean }} [options] */
  const onLegendGeometryChanged = (options) => legendGeometryChangedHandler?.(options);

  // The stroke Generate gives a row added here: the renderer's first Legend row
  // as it drew it, before that row's own stroke edit, whose captured original
  // Reset Stroke restores (OV-121, PD-OI-066). The draft's block stroke applies
  // on Generate, so it is not read here (R1).
  /**
   * @param {DrawingState} drawing
   * @param {Element} targetGroup
   * @returns {{ color: string | null, width: string | number | null }}
   */
  const templateRowStroke = (drawing, targetGroup) => {
    const template = targetGroup.querySelector('g[data-legend-key]');
    return rendererRowStroke(drawing, template, getLegendEntrySwatch(template)
      || targetGroup.querySelector('path[fill]:not([fill="none"]):not([fill^="url("])'));
  };

  // The stroke the renderer drew on Legend row `row` (its swatch), before the
  // row's own stroke edit, whose captured original Reset Stroke restores.
  /**
   * @param {DrawingState} drawing
   * @param {Element | null} row
   * @param {Element | null} [swatch]
   * @returns {{ color: string | null, width: string | number | null }}
   */
  const rendererRowStroke = (drawing, row, swatch = getLegendEntrySwatch(row)) => {
    const override = row ? drawing.legendStrokeOverrides[String(row.getAttribute('data-legend-key') || '').trim()] : null;
    /**
     * @param {unknown} value The row's own stroke edit, if any.
     * @param {string} original The override field that holds the drawn value.
     * @param {string} attribute
     * @param {() => string | number | null} inherited
     */
    const drawn = (value, original, attribute, inherited) => {
      if (value === undefined || value === null || value === '') return swatch?.getAttribute(attribute) ?? null;
      return Object.prototype.hasOwnProperty.call(override, original) ? override[original] : inherited();
    };
    return {
      color: drawn(override?.strokeColor, 'originalStrokeColor', 'stroke', () => originalSvgStroke.value.color),
      width: drawn(override?.strokeWidth, 'originalStrokeWidth', 'stroke-width', () => originalSvgStroke.value.width)
    };
  };

  const addLegendEntry = async (caption, color, options = {}) => {
    const drawing = state.activeDrawing();
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    const owner = String(options.owner || '').trim();
    const conflictPolicy = options.conflictPolicy || 'suffix';
    const shouldCommit = options.commit !== false;
    const shouldReflow = options.reflow !== false;
    console.log(`addLegendEntry called with caption="${caption}", color="${color}"`);
    if (!svgContainer.value) return false;

    const svg = options.svg || svgContainer.value.querySelector('svg');
    if (!svg) return false;
    if (!parseCompositionMetadata(svg).legendReflow) {
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

    // The row is a copy of the Legend's first row as the renderer drew it, as
    // Generate copies it (services/svg-result-ingestion.js, OV-86, OV-121),
    // appended last; the layout owner then places every row as Python lays
    // out these rows (zero shift).
    const templateStroke = templateRowStroke(drawing, targetGroup);
    /** @type {Element[]} */
    const added = [];
    for (const group of allTargetGroups) {
      const template = directLegendEntryGroups(group)[0];
      const entryGroup = /** @type {Element | null} */ (template?.cloneNode?.(true) || null);
      const swatch = getLegendEntrySwatch(entryGroup);
      if (!entryGroup || !swatch) {
        added.forEach((entry) => entry.remove());
        console.error('Failed to add legend entry: the Legend has no row to copy.');
        if (options.throwOnError) throw new Error('The Legend has no row to copy for a new entry.');
        return false;
      }
      entryGroup.removeAttribute('transform');
      entryGroup.setAttribute('data-legend-key', caption);
      if (owner) entryGroup.setAttribute('data-legend-owner', owner);
      else entryGroup.removeAttribute('data-legend-owner');
      const label = entryGroup.querySelector('text');
      if (label) label.textContent = caption;
      swatch.setAttribute('fill', color);
      if (templateStroke.color === null) swatch.removeAttribute('stroke');
      else swatch.setAttribute('stroke', String(templateStroke.color));
      if (templateStroke.width === null) swatch.removeAttribute('stroke-width');
      else swatch.setAttribute('stroke-width', String(templateStroke.width));
      group.appendChild(entryGroup);
      added.push(entryGroup);
    }

    if (shouldReflow) onLegendGeometryChanged();

    if (shouldCommit) {
      persistLegendReconciliation();
    }

    return caption;
  };

  const updateLegendEntryColorByCaption = (caption, color, { commit = true } = {}) => {
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
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
            if (normalizedColor(fill) === normalizedColor(color)) break;
            path.setAttribute('fill', color);
            updated = true;
            break;
          }
        }
      }
    }

    if (updated) {
      if (commit) persistLegendReconciliation();
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
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    if (!svgContainer.value) return false;

    const svg = svgContainer.value.querySelector('svg');
    if (!svg) return false;

    const targetGroups = getAllFeatureLegendGroups(svg);
    if (targetGroups.length === 0) return false;

    let removed = false;

    targetGroups.forEach((targetGroup, targetIndex) => {
      const entryGroup = targetGroup.querySelector(`g[data-legend-key="${CSS.escape(caption)}"]`);
      if (entryGroup) {
        rememberRetiredEntry(caption, targetGroup, targetIndex, entryGroup);
        rememberResultEntry(caption, targetGroup, targetIndex, entryGroup);
        entryGroup.remove();
        removed = true;
      }
    });

    if (removed) {
      // The layout owner lays the remaining rows out as Python would (zero
      // shift; OV-124) and docks the Legend.
      onLegendGeometryChanged();

      if (commit) persistLegendReconciliation();
      console.log(`Removed legend entry: "${caption}"`);
    }

    return removed;
  };

  const captureLegendEntryOwners = () => {
    const svg = svgContainer.value?.querySelector?.('svg');
    return svg ? getAllFeatureLegendGroups(svg).map((group, index) => ({
      target: targetGroupKey(group, index),
      entries: directLegendEntryGroups(group).map(entry => ({
        caption: String(entry.getAttribute('data-legend-key') || '').trim(),
        owner: String(entry.getAttribute('data-legend-owner') || '')
      }))
    })) : [];
  };

  // The one ordering of the mounted Legend (R3): each feature Legend group
  // takes the caption order through `orderLegendEntries`. History restore, the
  // B19 projection, and live Sort and Move (the port app/legend.js gives
  // sort-actions) share it. Null when no Legend is mounted.
  // `layOut: false` leaves the layout to a caller that lays the Legend out
  // once after more changes.
  const orderMountedLegend = (captionOrder, { keepFollowed = false, layOut = true } = {}) => {
    const targetGroups = getAllFeatureLegendGroups(svgContainer.value?.querySelector?.('svg'));
    if (targetGroups.length === 0) return null;
    const changed = targetGroups.reduce((moved, targetGroup) => (
      orderLegendEntries(targetGroup, captionOrder, { keepFollowed }) || moved
    ), false);
    // The rows in their new order take the places Python gives that order
    // (zero shift), not the slots of the old order.
    if (changed && layOut) onLegendGeometryChanged();
    return changed;
  };

  const legendCaption = (entry) => String(entry?.caption || '').trim();
  const generatedCaption = (entry) => String(entry?.originalCaption || entry?.caption || '').trim();

  // Whether the Legend editor changed the rows Python drew: an added row, a
  // deleted row, a renamed row, or another order.
  /** @param {DrawingState} drawing @param {SVGSVGElement} svg */
  const hasLegendRowEdits = (drawing, svg) => (
    drawing.deletedLegendEntries.value.length > 0
    || getAllFeatureLegendGroups(svg).some((group) => group.querySelector('g[data-legend-owner="direct-editor"]'))
    || (drawing.legendEntries.value || []).some((entry) => entry?.originalCaption && entry.caption !== entry.originalCaption)
    || isLegendOrderEdited(drawing.legendEntries.value || [], originalLegendOrder.value || [])
  );

  // Python lays out only the rows it draws. Generate, a rerender, and the
  // display of a batch Result apply the rows the editor added, deleted,
  // renamed, and ordered (services/svg-result-ingestion.js), so a mounted
  // Result with such edits is laid out by the layout owner as Python lays the
  // edited rows out, as the live edit was (zero shift; OV-122, OV-124,
  // OV-126, OV-127, PD-OI-066). Returns whether it was laid out again.
  /** @param {SVGSVGElement} svg */
  // The mount binder commits the Result once its later steps have bound it.
  const layOutMountedLegendEdits = (svg) => {
    if (getAllFeatureLegendGroups(svg).length === 0 || !hasLegendRowEdits(state.activeDrawing(), svg)) return false;
    onLegendGeometryChanged({ commit: false });
    return true;
  };

  // Whether a captured Legend list lists exactly the entries the mounted
  // Result draws, so that it describes this Result.
  const describesMountedLegend = (svg, entries) => {
    const drawn = Array.from(getVisibleFeatureLegendGroup(svg)?.querySelectorAll?.('g[data-legend-key]') || [])
      .map((group) => String(group.getAttribute('data-legend-key') || '').trim());
    const listed = new Set(entries.map(legendCaption));
    return drawn.length === listed.size && drawn.every((caption) => listed.has(caption));
  };

  // B19 (D-07, D-08): a History step made on another batch Result. The
  // mounted Result receives the step's shared Legend intent, from its other
  // side `from` to the restored list: captions, colors, deletions, direct
  // additions, and the order. An entry that only another Result draws is never
  // copied here and an entry only this Result draws is never removed; a
  // returning entry comes from this Result's retired entries or is a direct
  // addition. The Legend panel then lists this Result.
  /** @param {DrawingState} drawing */
  const projectLegendChange = (drawing, targetGroups, restored, from, entryOwners) => {
    const captionMap = (list, key) => new Map(list.filter(legendCaption).map((entry) => [key(entry), entry]));
    const before = captionMap(from, legendCaption);
    const beforeByGenerated = captionMap(from, generatedCaption);
    const after = captionMap(restored, legendCaption);
    const afterByGenerated = captionMap(restored, generatedCaption);
    const renames = new Map();
    before.forEach((entry, caption) => {
      const next = afterByGenerated.get(generatedCaption(entry));
      if (next && !after.has(caption) && !before.has(legendCaption(next))) renames.set(caption, legendCaption(next));
    });
    const renamed = new Set(renames.values());
    const removed = [...before.keys()].filter((caption) => !after.has(caption) && !renames.has(caption));
    const added = [...after.keys()].filter((caption) => !before.has(caption) && !renamed.has(caption));
    const orderChanged = [...before.keys()].join('\u0000') !== [...after.keys()].join('\u0000');
    const identity = readActiveResultIdentity?.();
    const retired = retiredEntriesByResult.get(identity);
    let changed = false;
    let removedEntry = false;
    let returnedEntry = false;

    targetGroups.forEach((targetGroup, targetIndex) => {
      const targetKey = targetGroupKey(targetGroup, targetIndex);
      const owners = new Map((entryOwners.find(group => group.target === targetKey)?.entries || [])
        .map(entry => [entry.caption, entry.owner]));
      const setCaption = (group, caption) => {
        group.setAttribute('data-legend-key', caption);
        const text = group.querySelector('text');
        if (text) text.textContent = caption;
      };
      renames.forEach((caption, previous) => {
        const group = findLegendEntryGroup(targetGroup, previous);
        if (!group || findLegendEntryGroup(targetGroup, caption)) return;
        setCaption(group, caption);
        changed = true;
      });
      removed.forEach((caption) => {
        const group = findLegendEntryGroup(targetGroup, caption);
        if (!group) return;
        rememberRetiredEntry(caption, targetGroup, targetIndex, group);
        rememberResultEntry(caption, targetGroup, targetIndex, group, identity);
        group.remove();
        changed = removedEntry = true;
      });
      added.forEach((caption) => {
        const key = `${targetKey}\u0000${caption}`;
        const direct = owners.get(caption) === 'direct-editor';
        const template = retired?.get(key)?.node || (direct
          ? restoredEntryTemplate(caption, targetGroup, targetIndex) || directLegendEntryGroups(targetGroup)[0]
          : null);
        if (!template || findLegendEntryGroup(targetGroup, caption)) return;
        const group = targetGroup.appendChild(template.cloneNode(true));
        setCaption(group, caption);
        if (direct) group.setAttribute('data-legend-owner', 'direct-editor');
        retired?.delete(key);
        changed = returnedEntry = true;
      });
      after.forEach((entry, caption) => {
        const group = findLegendEntryGroup(targetGroup, caption);
        const previous = before.get(caption) || beforeByGenerated.get(generatedCaption(entry));
        const stepColor = !previous || normalizedColor(previous.color) !== normalizedColor(entry.color);
        const color = drawing.legendColorOverrides[caption] || (stepColor ? entry.color : '');
        if (group && color && setLegendEntryColor(group, String(color))) changed = true;
      });
    });
    // A Result already in the order keeps its own entries' places (B18); a
    // returned entry takes its place in the order, and own entries follow. A
    // step that leaves the default order of the Result it was made on gives
    // this Result its own default order (OV-47).
    // No inventory is stored under an empty identity, so the guard only states what get() returned.
    const ownInventory = identity ? inventoryByResult.get(identity) : undefined;
    const directCaptions = new Set(entryOwners.flatMap((group) => group.entries)
      .filter((entry) => entry.owner === 'direct-editor').map((entry) => entry.caption));
    const generatedEntries = restored.filter((entry) => !directCaptions.has(legendCaption(entry)));
    const restoresDefault = [...inventoryByResult.values()].some((inventory) => (
      generatedEntries.every((entry) => inventory.includes(generatedCaption(entry)))
      && !isLegendOrderEdited(restored, inventory)
    ));
    const order = ownInventory?.length && restoresDefault
      ? defaultLegendCaptionOrder(restored, ownInventory)
      : [...after.keys()];
    if ((orderChanged || returnedEntry) && orderMountedLegend(order, { keepFollowed: !returnedEntry, layOut: false })) {
      changed = true;
    }
    if (changed || removedEntry || returnedEntry) {
      onLegendGeometryChanged();
      persistLegendReconciliation();
    }
    extractLegendEntries();
    return changed;
  };

  // History restore of the Legend (R3, R11), from the Legend side of the
  // History step: its entry owners and `from`, the list the step leaves. The
  // restored list installs as is on the Result it describes: outside a batch,
  // or when `from` lists the entries the mounted Result draws. Otherwise the
  // step was made on another batch Result and only its shared intent is
  // projected (B19).
  const reconcileLegendEntries = ({ entryOwners = [], from = null } = {}) => {
    const drawing = state.activeDrawing();
    const svg = svgContainer.value?.querySelector?.('svg');
    if (!svg) return false;
    const targetGroups = getAllFeatureLegendGroups(svg);
    if (targetGroups.length === 0) return false;
    if (results.value.length > 1 && Array.isArray(from) && !describesMountedLegend(svg, from)) {
      return projectLegendChange(drawing, targetGroups, drawing.legendEntries.value || [], from, entryOwners);
    }

    const desiredEntries = [];
    let entryColorStateChanged = false;
    const seenCaptions = new Set();
    (Array.isArray(drawing.legendEntries.value) ? drawing.legendEntries.value : []).forEach((entry) => {
      const caption = String(entry?.caption || '').trim();
      if (!caption || seenCaptions.has(caption)) return;
      seenCaptions.add(caption);
      const color = String(drawing.legendColorOverrides[caption] || entry?.color || '#cccccc');
      if (normalizedColor(entry?.color) !== normalizedColor(color)) {
        entryColorStateChanged = true;
      }
      desiredEntries.push({ ...entry, caption, color });
    });

    let changed = false;
    const restoredCaptions = new Set();

    targetGroups.forEach((targetGroup, targetIndex) => {
      const owners = new Map((entryOwners.find(group => group.target === targetGroupKey(targetGroup, targetIndex))?.entries || [])
        .map(entry => [entry.caption, entry.owner]));
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

      desiredEntries.forEach((entry) => {
        const group = assigned.get(entry.caption);
        if (!group) return;
        if (owners.has(entry.caption)) {
          const owner = owners.get(entry.caption);
          if (String(group.getAttribute('data-legend-owner') || '') !== owner) {
            if (owner) group.setAttribute('data-legend-owner', owner);
            else group.removeAttribute('data-legend-owner');
            changed = true;
          }
        }
        if (setLegendEntryColor(group, entry.color)) changed = true;
        const text = group.querySelector('text');
        if (text && String(text.textContent || '') !== entry.caption) {
          text.textContent = entry.caption;
          changed = true;
        }
      });
    });
    // The entries take the restored order in the Legend's slots, as live Sort
    // and Move place them.
    if (orderMountedLegend(desiredEntries.map((entry) => entry.caption), { layOut: false })) changed = true;

    if (entryColorStateChanged) drawing.legendEntries.value = desiredEntries;
    if (!changed) return entryColorStateChanged;
    restoredCaptions.forEach((caption) => retiredEntryTemplates.delete(caption));
    onLegendGeometryChanged();
    persistLegendReconciliation();
    return true;
  };

  // The legend rows the specific-color rules draw, prepared on a disposable copy
  // of the mounted legend (R13): the root (its composition metadata) and the
  // `#legend` group only, the one part this reads and replaces. The rule owner
  // runs its transition and then `apply` in one History step, while `isCurrent`
  // holds; no rule mutation runs here. Resolves to false when the Result or the
  // caller became stale.
  /**
   * @param {Record<string, any>[]} intents
   * @param {{ drawing?: DrawingState, previousFileIntents?: Record<string, any>[], isCurrent?: () => boolean, placement?: { caption: string, at: string } | null }} [options]
   */
  const prepareFileLegendEntries = async (intents, {
    drawing = state.activeDrawing(), previousFileIntents = [], isCurrent = () => true, placement = null
  } = {}) => {
    const mountedSvg = svgContainer.value?.querySelector('svg');
    const mountedLegendGroup = mountedSvg?.getElementById('legend');
    const svg = mountedSvg && mountedLegendGroup ? /** @type {SVGSVGElement} */ (mountedSvg.cloneNode(false)) : null;
    if (svg && mountedLegendGroup) svg.appendChild(mountedLegendGroup.cloneNode(true));
    const targetGroups = svg ? getAllFeatureLegendGroups(svg) : [];
    if (!svg || targetGroups.length === 0) {
      return { diff: { add: [], update: [], remove: [], unchanged: [] }, isCurrent: () => true, apply: () => {} };
    }

    const provenance = new Map();
    for (const entry of previousFileIntents) {
      const caption = String(entry?.caption || '').trim();
      if (!provenance.has(caption)) provenance.set(caption, new Set());
      provenance.get(caption).add(normalizedColor(entry?.color));
    }

    const desiredByCaption = new Map(intents.map((intent) => [intent.caption, normalizedColor(intent.color)]));
    targetGroups.forEach((group) => {
      Array.from(group.querySelectorAll('g[data-legend-key]')).forEach((entry) => {
        const caption = entry.getAttribute('data-legend-key') || '';
        if (
          !entry.hasAttribute('data-legend-owner') &&
          provenance.get(caption)?.has(legendEntryColor(entry))
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
        svg,
        owner: SPECIFIC_COLOR_FILE_OWNER,
        conflictPolicy: 'error',
        commit: false,
        reflow: false,
        throwOnError: true
      });
    }

    // OV-158 (Owner decision 2026-10-07): the row a rename draws takes the
    // place of the row it renames, through the one ordering of a Legend, so
    // the editor list records the edited order Generate replays (PD-OI-063).
    if (placement && intents.some((intent) => intent.caption === placement.caption)) {
      const order = (drawing.legendEntries.value || []).map(legendCaption)
        .map((/** @type {string} */ caption) => (caption === placement.at ? placement.caption : caption));
      targetGroups.forEach((group) => orderLegendEntries(group, order));
    }

    if (!isCurrent() || svgContainer.value.querySelector('svg') !== mountedSvg) return false;
    // The copy's rows are laid out once they are mounted: `apply` hands the
    // Legend to the layout owner (zero shift).
    return {
      diff,
      isCurrent: () => svgContainer.value.querySelector('svg') === mountedSvg,
      // Mounted geometry and the Result commit synchronously inside History.
      // An unchanged Legend (an empty diff) is neither laid out nor committed
      // again; a changed one is laid out and committed once.
      apply: () => {
        const mountedLegend = mountedSvg.getElementById('legend');
        const candidateLegend = svg.getElementById('legend');
        if (!mountedLegend || !candidateLegend || !mountedLegend.isEqualNode(candidateLegend)) {
          if (mountedLegend && candidateLegend) mountedLegend.replaceWith(candidateLegend);
          onLegendGeometryChanged({ commit: false });
          commitActiveResultEdit?.('legend-file-sync');
        }
        extractLegendEntries();
      }
    };
  };

  /**
   * Read the Legend of a mounted Result: its entries in visual order and the
   * captions the renderer generated (not the editor's direct entries). The one
   * reader of a mounted Legend; `previousEntries` keep stroke, feature ids,
   * and the generated caption of a renamed entry; a renamed row an earlier
   * Generate hid (`dormantEntries`, OV-120) keeps its generated caption too.
   * @param {SVGSVGElement} svg
   * @param {Record<string, any>[]} previousEntries
   * @param {Record<string, any>[]} [dormantEntries]
   */
  const readMountedLegend = (svg, previousEntries, dormantEntries = []) => {
    const targetGroup = getVisibleFeatureLegendGroup(svg);
    if (!targetGroup) return null;

    const entries = [];
    const generatedCaptions = new Set();

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

      const existingEntry = previousEntries.find((entry) => (
        entry.caption === caption
        && normalizedColor(entry.color) === normalizedColor(color)
      )) || dormantEntries.find((entry) => entry.caption === caption);
      const existingFeatureIds = existingEntry?.featureIds || [];
      const originalCaption = existingEntry?.originalCaption || caption;
      if (entryGroup.getAttribute('data-legend-owner') !== 'direct-editor') {
        generatedCaptions.add(originalCaption);
      }

      entries.push({
        caption,
        originalCaption,
        color,
        xPos,
        yPos,
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

    return { entries: visuallySortedEntries, generatedCaptions };
  };

  // The inventory after a draw: while an edited order is replayed, the
  // surviving categories keep their default order and new ones follow it;
  // otherwise the Result shows the renderer's order, which becomes the default
  // order. Deleted entries keep their place.
  /**
   * @param {DrawingState} drawing
   * @param {{ inventory: string[], orderEdited: boolean, generatedCaptions: Set<string>, rendered: string[] }} drawn
   */
  const drawnInventory = (drawing, { inventory, orderEdited, generatedCaptions, rendered }) => {
    const retainedCaptions = new Set([
      ...generatedCaptions,
      ...drawing.deletedLegendEntries.value.map(entry => entry.originalCaption || entry.caption)
    ]);
    const surviving = inventory.filter(caption => retainedCaptions.has(caption));
    if (orderEdited) return [...new Set([...surviving, ...rendered])];
    // A deleted entry the renderer no longer shows goes before the first
    // rendered entry that followed it, so a Restore returns it there (OV-154).
    const order = [...new Set(rendered)];
    surviving.filter(caption => !order.includes(caption)).forEach((caption) => {
      const next = inventory.slice(inventory.indexOf(caption) + 1).find(later => order.includes(later));
      order.splice(next === undefined ? order.length : order.indexOf(next), 0, caption);
    });
    return order;
  };
  const renderedCaptions = ({ entries, generatedCaptions }) => entries
    .map(entry => entry.originalCaption)
    .filter(caption => generatedCaptions.has(caption));

  /**
   * The renamed rows a draw leaves out (OV-120): each keeps its generated
   * caption and its rename, unless the row is deleted.
   * @param {DrawingState} drawing
   * @param {{ previous: Record<string, any>[], drawn: Record<string, any>[] }} rows
   * @returns {Record<string, any>[]}
   */
  const dormantAfterDraw = (drawing, { previous, drawn }) => {
    const drawnOriginals = new Set(drawn.map((entry) => String(entry?.originalCaption || entry?.caption || '').trim()));
    const deleted = new Set(drawing.deletedLegendEntries.value
      .map((entry) => String(entry?.originalCaption || entry?.caption || '').trim()));
    /** @type {Set<string>} */
    const seen = new Set();
    return previous.flatMap((entry) => {
      const caption = String(entry?.caption || '').trim();
      const originalCaption = String(entry?.originalCaption || '').trim();
      if (!caption || !originalCaption || caption === originalCaption || drawnOriginals.has(originalCaption)
        || deleted.has(originalCaption) || seen.has(originalCaption)) return [];
      seen.add(originalCaption);
      return [{ caption, originalCaption, color: entry.color, featureIds: [] }];
    });
  };

  /** @param {string[]} liveResultIdentities */
  const pruneResultInventories = (liveResultIdentities) => {
    const live = new Set(liveResultIdentities);
    [...inventoryByResult.keys()].forEach((key) => {
      if (!live.has(key)) inventoryByResult.delete(key);
    });
    [...replayedInventoryResults].forEach((key) => {
      if (!live.has(key)) replayedInventoryResults.delete(key);
    });
  };

  // A Result about to be displayed has no inventory yet: read it from the
  // mounted Legend as the renderer drew it. The stored SVG of a Result changes
  // only through edits made while it is displayed, so this runs before the
  // editor intent is projected onto it. The displayed Result keeps its own
  // inventory in `originalLegendOrder` until `adoptResultInventory`.
  /**
   * @param {SVGSVGElement} svg
   * @param {{ resultIdentity?: string, liveResultIdentities?: string[] }} [options]
   * @returns {string[]} The Result's inventory.
   */
  const captureResultInventory = (svg, { resultIdentity: identity = '', liveResultIdentities = [] } = {}) => {
    const drawing = state.activeDrawing();
    pruneResultInventories(liveResultIdentities);
    const stored = inventoryByResult.get(identity);
    if (!identity || stored) return stored || [];
    const mounted = readMountedLegend(svg, drawing.legendEntries.value || [], drawing.dormantLegendEntries.value || []);
    const inventory = mounted
      ? drawnInventory(drawing, {
        inventory: originalLegendOrder.value,
        orderEdited: replayedInventoryResults.has(identity),
        generatedCaptions: mounted.generatedCaptions,
        rendered: renderedCaptions(mounted)
      })
      : [];
    inventoryByResult.set(identity, inventory);
    replayedInventoryResults.delete(identity);
    return inventory;
  };

  // E1: a Result of the other diagram mode keeps the inventory it was drawn
  // with while that mode is not shown; a Result loaded from a Session's
  // `otherModeResult` arrives with its saved one. Kept only for a Result with
  // no stored inventory; `adoptResultInventory` adopts it when it is displayed.
  /** @param {string | undefined} identity @param {string[]} inventory */
  const rememberResultInventory = (identity, inventory) => {
    if (!identity || inventoryByResult.has(identity) || !Array.isArray(inventory)) return;
    inventoryByResult.set(identity, [...inventory]);
  };

  // `originalLegendOrder` is the displayed Result's inventory: after a Result
  // is displayed and after a History restore wrote another Result's value.
  // `restored` marks a restore (Session Load, checkpoint) that installed the
  // displayed Result's own inventory.
  /** @param {string | undefined} identity @param {{ restored?: boolean }} [options] */
  const adoptResultInventory = (identity, { restored = false } = {}) => {
    if (!identity) return;
    const stored = inventoryByResult.get(identity);
    if (stored && !restored) originalLegendOrder.value = [...stored];
    else inventoryByResult.set(identity, [...originalLegendOrder.value]);
  };

  /** @param {{ replaceGeneratedInventory?: boolean, liveResultIdentities?: string[] }} [options] */
  const extractLegendEntries = ({ replaceGeneratedInventory = false, liveResultIdentities = [] } = {}) => {
    const drawing = state.activeDrawing();
    if (!svgContainer.value) {
      drawing.legendEntries.value = [];
      return;
    }

    const svg = svgContainer.value.querySelector('svg');
    if (!svg) {
      drawing.legendEntries.value = [];
      return;
    }

    const previousEntries = drawing.legendEntries.value || [];
    const previousDormant = drawing.dormantLegendEntries.value || [];
    const mounted = readMountedLegend(svg, previousEntries, previousDormant);
    if (!mounted) {
      drawing.legendEntries.value = [];
      return;
    }
    const { entries: visuallySortedEntries, generatedCaptions } = mounted;

    // The entries shown before this extraction decide, as in Generate, whether
    // the displayed Result replays an edited order (D-08).
    const orderEdited = isLegendOrderEdited(drawing.legendEntries.value, originalLegendOrder.value);
    drawing.legendEntries.value = visuallySortedEntries;
    // OV-120: a draw keeps the rename of each row it does not draw; the rename
    // applies again when a later Generate draws the row.
    if (replaceGeneratedInventory) {
      drawing.dormantLegendEntries.value = dormantAfterDraw(drawing, {
        previous: [...previousEntries, ...previousDormant],
        drawn: visuallySortedEntries
      });
    }

    const identity = readActiveResultIdentity?.() || '';
    if (replaceGeneratedInventory || originalLegendOrder.value.length === 0) {
      // Keep deletion intent; live editor extraction alone must not advance
      // the accepted generated inventory.
      originalLegendOrder.value = drawnInventory(drawing, {
        inventory: originalLegendOrder.value,
        orderEdited,
        generatedCaptions,
        rendered: renderedCaptions(mounted)
      });
      if (identity) {
        // A draw that replays an edited order shows it on every Result.
        if (replaceGeneratedInventory && orderEdited) {
          liveResultIdentities.filter((key) => key !== identity).forEach((key) => {
            inventoryByResult.delete(key);
            replayedInventoryResults.add(key);
          });
        }
        inventoryByResult.set(identity, [...originalLegendOrder.value]);
      }
    } else if (identity && !inventoryByResult.has(identity)) {
      inventoryByResult.set(identity, [...originalLegendOrder.value]);
    }

    if (Object.keys(originalLegendColors.value).length === 0 && visuallySortedEntries.length > 0) {
      visuallySortedEntries.forEach((entry) => {
        originalLegendColors.value[entry.caption] = entry.color;
      });
    }
  };

  // OV-276: the Legend panel rows the palette colors take the color the
  // palette repaint gave their swatches. True when a row changed.
  /** @param {unknown} colorsByCaption */
  const setPaletteLegendEntryColors = (colorsByCaption) => {
    if (!(colorsByCaption instanceof Map) || colorsByCaption.size === 0) return false;
    const drawing = state.activeDrawing();
    let changed = false;
    const entries = (drawing.legendEntries.value || []).map((entry) => {
      const color = colorsByCaption.get(entry?.caption);
      if (!color || normalizedColor(color) === normalizedColor(entry.color)) return entry;
      changed = true;
      return { ...entry, color };
    });
    if (changed) drawing.legendEntries.value = entries;
    return changed;
  };

  const updateLegendEntryColor = (idx, newColor) => {
    const drawing = state.activeDrawing();
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    if (!svgContainer.value) return false;
    const svg = svgContainer.value.querySelector('svg');
    if (!svg) return false;

    const entry = drawing.legendEntries.value[idx];
    if (!entry) return false;

    const caption = entry.caption;
    if (!caption) return false;

    const targetGroups = getAllFeatureLegendGroups(svg);
    if (targetGroups.length === 0) return false;

    let changed = false;

    for (const targetGroup of targetGroups) {
      const entryGroup = targetGroup.querySelector(`g[data-legend-key="${CSS.escape(caption)}"]`);
      if (entryGroup) {
        const paths = entryGroup.querySelectorAll('path');
        for (const path of paths) {
          const fill = path.getAttribute('fill');
          if (fill && fill !== 'none' && !fill.startsWith('url(')) {
            if (normalizedColor(fill) !== normalizedColor(newColor)) {
              path.setAttribute('fill', newColor);
              changed = true;
            }
            break;
          }
        }
      }
    }

    const stateChanged = normalizedColor(entry.color) !== normalizedColor(newColor);
    if (!changed && !stateChanged) return false;
    entry.color = newColor;
    drawing.legendColorOverrides[caption] = newColor;
    if (changed) persistLegendReconciliation();
    return true;
  };

  const updateLegendEntryCaption = (idx, newCaption) => {
    const drawing = state.activeDrawing();
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    if (!svgContainer.value) return false;
    const svg = svgContainer.value.querySelector('svg');
    if (!svg) return false;

    const entry = drawing.legendEntries.value[idx];
    if (!entry) return false;

    const oldCaption = entry.caption;
    const caption = newCaption.trim();
    if (!caption || caption === oldCaption) return false;

    const targetGroups = getAllFeatureLegendGroups(svg);
    if (targetGroups.length === 0) return false;

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

    if (drawing.legendColorOverrides[oldCaption]) {
      drawing.legendColorOverrides[caption] = drawing.legendColorOverrides[oldCaption];
      delete drawing.legendColorOverrides[oldCaption];
    }

    if (drawing.legendStrokeOverrides[oldCaption]) {
      drawing.legendStrokeOverrides[caption] = drawing.legendStrokeOverrides[oldCaption];
      delete drawing.legendStrokeOverrides[oldCaption];
    }

    entry.caption = caption;
    onLegendGeometryChanged();

    persistLegendReconciliation();
    return true;
  };

  const addNewLegendEntry = async () => {
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    if (!newLegendCaption.value.trim()) return;

    const added = await addLegendEntry(newLegendCaption.value.trim(), newLegendColor.value, { owner: 'direct-editor' });
    if (added) {
      newLegendCaption.value = '';
      newLegendColor.value = '#808080';
      // The editor lists the row within the add, so the add's History step
      // holds it (OV-125).
      extractLegendEntries();
    }
  };

  const deleteLegendEntry = (idx) => {
    const drawing = state.activeDrawing();
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    const entry = drawing.legendEntries.value[idx];
    if (!entry) return;

    drawing.deletedLegendEntries.value.push({ ...entry });

    removeLegendEntry(entry.caption);
    extractLegendEntries();
  };

  // The editor list with `entry` at its place in the default order: before the
  // first entry the generated inventory ranks after it, else after the last
  // ranked entry; an entry the inventory does not rank goes last, as the
  // default order puts the editor's own rows (`defaultLegendEntryOrder`).
  /**
   * @param {Record<string, any>[]} entries
   * @param {Record<string, any>} entry
   * @param {string[]} inventory
   */
  const withEntryAtDefaultPlace = (entries, entry, inventory) => {
    /** @param {Record<string, any>} item */
    const rank = (item) => inventory.indexOf(generatedCaption(item));
    const own = rank(entry);
    if (own < 0) return [...entries, entry];
    const next = entries.findIndex((item) => rank(item) > own);
    const at = next >= 0 ? next : entries.reduce((last, item, index) => (rank(item) >= 0 ? index + 1 : last), 0);
    return [...entries.slice(0, at), entry, ...entries.slice(at)];
  };

  // OV-167: the stroke Generate draws on a restored row this page holds no
  // removed copy of (after a Session load, or a row removed before the page
  // opened). Python strokes a row by its kind: the rows of drawn features
  // with the feature block stroke, the GC and Depth rows with their track's
  // stroke. A drawn row of the same kind shows it; else, for a feature row, a
  // block of a feature the row colors shows the block stroke. Null keeps the
  // add path's stroke, the first row's (OV-121).
  /**
   * @param {DrawingState} drawing
   * @param {SVGSVGElement} svg
   * @param {Element} group
   * @param {string} caption
   * @param {Set<string>} returning The captions this Restore returns.
   * @returns {{ color: string | number | null, width: string | number | null } | null}
   */
  const restoredRowStroke = (drawing, svg, group, caption, returning) => {
    const entries = drawing.legendEntries.value || [];
    /** @param {string} key */
    const drawsFeatures = (key) => mountedLegendRowFeatureIds(svg, key, entries).length > 0;
    const featureRow = drawsFeatures(caption);
    const peer = directLegendEntryGroups(group).find((entry) => {
      const key = String(entry.getAttribute('data-legend-key') || '').trim();
      return key !== caption && !returning.has(key) && drawsFeatures(key) === featureRow;
    });
    if (peer) return rendererRowStroke(drawing, peer);
    return featureRow ? featureBlockStroke(drawing, svg, caption) : null;
  };

  // Python's feature block stroke, as a block of a feature row `caption`
  // colors draws it: not an automatic underlay (drawn without a stroke), and
  // not a feature or row with a stroke edit, unless the edit kept the drawn
  // original. Null when no such block is drawn.
  /**
   * @param {DrawingState} drawing
   * @param {SVGSVGElement} svg
   * @param {string} caption
   * @returns {{ color: string | number | null, width: string | number | null } | null}
   */
  const featureBlockStroke = (drawing, svg, caption) => {
    const rowEdit = drawing.legendStrokeOverrides[caption];
    if (rowEdit && Object.prototype.hasOwnProperty.call(rowEdit, 'originalStrokeColor')) {
      return { color: rowEdit.originalStrokeColor, width: rowEdit.originalStrokeWidth ?? null };
    }
    if (setsFeatureStroke(rowEdit)) return null;
    const featureEdits = drawing.featureStrokeOverrides;
    const ownStrokeIds = (state.extractedFeatures?.value || [])
      .filter((/** @type {Record<string, any>} */ feature) => setsFeatureStroke(featureEdits[featureOverrideKey(feature)]))
      .map((/** @type {Record<string, any>} */ feature) => String(feature?.svg_id || '').trim());
    const featureIndex = getFeatureElementIndex(svg);
    for (const id of mountedLegendRowFeatureIds(svg, caption, drawing.legendEntries.value || [], { ownStrokeIds })) {
      const block = getFeatureFillElements(svg, id, featureIndex)
        .find((element) => !isAutoFeatureUnderlay(element));
      if (block) return { color: block.getAttribute('stroke'), width: block.getAttribute('stroke-width') };
    }
    return null;
  };

  // OV-154: Restore and Restore all in the Legend editor. A restored row leaves
  // the deleted list, so Generate draws it again, and the displayed Result
  // draws it now: the add path draws the row in the color Generate gives it,
  // with the swatch stroke the row had when it was removed (the stroke
  // Generate draws for its kind when this page did not remove it, OV-167),
  // the editor list takes it at its place in the default order (PD-OI-063),
  // the one ordering of the mounted Legend moves it there, and the layout
  // owner lays the Legend out once, as Python lays these rows out (zero
  // shift). Resolves to whether a row returned.
  /** @param {number[] | null} [indexes] Indexes into the deleted list; every row when omitted. */
  const restoreDeletedLegendEntries = async (indexes = null) => {
    const drawing = state.activeDrawing();
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    const svg = svgContainer.value?.querySelector('svg');
    const deleted = drawing.deletedLegendEntries.value || [];
    const picked = new Set(indexes ?? deleted.map((/** @type {unknown} */ _, /** @type {number} */ index) => index));
    const returning = deleted.filter((/** @type {unknown} */ _, /** @type {number} */ index) => picked.has(index));
    if (!svg || returning.length === 0) return false;
    const inventory = (originalLegendOrder.value || []).map((/** @type {unknown} */ caption) => String(caption || '').trim());
    const context = {
      rules: drawing.manualSpecificRules,
      legendEntries: drawing.legendEntries.value || [],
      originalLegendOrder: inventory
    };
    let entries = [...(drawing.legendEntries.value || [])];
    const returningCaptions = new Set(returning.map(legendCaption));
    for (const entry of returning) {
      const caption = legendCaption(entry);
      const rowRules = legendRowRules(caption, context);
      const owner = inventory.includes(generatedCaption(entry)) ? ''
        : rowRules.length ? SPECIFIC_COLOR_FILE_OWNER : 'direct-editor';
      const color = String(drawing.legendColorOverrides[caption] || rowRules[0]?.color || entry.color);
      const added = await addLegendEntry(caption, color, { owner, conflictPolicy: 'error', commit: false, reflow: false });
      if (added !== caption) return false;
      getAllFeatureLegendGroups(svg).forEach((group, index) => {
        const removed = getLegendEntrySwatch(restoredEntryTemplate(caption, group, index));
        const swatch = getLegendEntrySwatch(findLegendEntryGroup(group, caption));
        if (!swatch) return;
        const stroke = removed
          ? { color: removed.getAttribute('stroke'), width: removed.getAttribute('stroke-width') }
          : restoredRowStroke(drawing, svg, group, caption, returningCaptions);
        if (!stroke) return;
        if (stroke.color === null) swatch.removeAttribute('stroke');
        else swatch.setAttribute('stroke', String(stroke.color));
        if (stroke.width === null) swatch.removeAttribute('stroke-width');
        else swatch.setAttribute('stroke-width', String(stroke.width));
      });
      entries = withEntryAtDefaultPlace(entries, { ...entry, color }, inventory);
    }
    drawing.deletedLegendEntries.value = deleted.filter((/** @type {unknown} */ _, /** @type {number} */ index) => !picked.has(index));
    drawing.legendEntries.value = entries;
    orderMountedLegend(entries.map(legendCaption), { layOut: false });
    onLegendGeometryChanged();
    persistLegendReconciliation();
    extractLegendEntries();
    return true;
  };

  return {
    addLegendEntry,
    captureLegendEntryOwners,
    adoptResultInventory,
    addNewLegendEntry,
    captureResultInventory,
    deleteLegendEntry,
    extractLegendEntries,
    legendEntryExists,
    hasRetiredResultLegend,
    layOutMountedLegendEdits,
    onLegendGeometryChanged,
    orderMountedLegend,
    prepareDisplayedResultLegend,
    prepareFileLegendEntries,
    rememberResultInventory,
    removeLegendEntry,
    reconcileLegendEntries,
    restoreDeletedLegendEntries,
    setLegendGeometryChangedHandler,
    setPaletteLegendEntryColors,
    updateLegendEntryCaption,
    updateLegendEntryColor,
    updateLegendEntryColorByCaption
  };
};
