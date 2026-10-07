// @ts-check
import { resolveColorToHex, toNativeColorInputValue } from '../../utils/color-utils.js';
import {
  defaultLegendCaptionOrder,
  getAllFeatureLegendGroups,
  getLegendEntrySwatch,
  getVisibleFeatureLegendGroup,
  isLegendOrderEdited,
  orderLegendEntries,
  parseTransformXY
} from './utils.js';
import { parseCompositionMetadata } from '../legend-layout/composition-actions.js';
import {
  diffLegendIntents,
  SPECIFIC_COLOR_FILE_OWNER
} from '../../services/specific-color-rules.js';

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
    selectedResultIndex,
    svgContainer,
    legendEntries,
    deletedLegendEntries,
    originalLegendOrder,
    originalLegendColors,
    newLegendCaption,
    newLegendColor,
    legendStrokeOverrides,
    legendColorOverrides,
    originalSvgStroke
  } = state;

  /** @type {(() => void) | null} */
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

  const onLegendGeometryChanged = () => legendGeometryChangedHandler?.();

  // The stroke Generate gives a row added here: the renderer's first Legend row
  // as it drew it, before that row's own stroke edit, whose captured original
  // Reset Stroke restores (OV-121, PD-OI-066). The draft's block stroke applies
  // on Generate, so it is not read here (R1).
  /**
   * @param {Element} targetGroup
   * @returns {{ color: string | null, width: string | number | null }}
   */
  const templateRowStroke = (targetGroup) => {
    const template = targetGroup.querySelector('g[data-legend-key]');
    const swatch = getLegendEntrySwatch(template)
      || targetGroup.querySelector('path[fill]:not([fill="none"]):not([fill^="url("])');
    const override = template ? legendStrokeOverrides[String(template.getAttribute('data-legend-key') || '').trim()] : null;
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
    const templateStroke = templateRowStroke(targetGroup);
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
  /** @param {SVGSVGElement} svg */
  const hasLegendRowEdits = (svg) => (
    deletedLegendEntries.value.length > 0
    || getAllFeatureLegendGroups(svg).some((group) => group.querySelector('g[data-legend-owner="direct-editor"]'))
    || (legendEntries.value || []).some((entry) => entry?.originalCaption && entry.caption !== entry.originalCaption)
    || isLegendOrderEdited(legendEntries.value || [], originalLegendOrder.value || [])
  );

  // Python lays out only the rows it draws. Generate, a rerender, and the
  // display of a batch Result apply the rows the editor added, deleted,
  // renamed, and ordered (services/svg-result-ingestion.js), so a mounted
  // Result with such edits is laid out by the layout owner as Python lays the
  // edited rows out, as the live edit was (zero shift; OV-122, OV-124,
  // OV-126, OV-127, PD-OI-066). Returns whether it was laid out again.
  /** @param {SVGSVGElement} svg */
  const layOutMountedLegendEdits = (svg) => {
    if (getAllFeatureLegendGroups(svg).length === 0 || !hasLegendRowEdits(svg)) return false;
    onLegendGeometryChanged();
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
  const projectLegendChange = (targetGroups, restored, from, entryOwners) => {
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
        const color = legendColorOverrides[caption] || (stepColor ? entry.color : '');
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
    const svg = svgContainer.value?.querySelector?.('svg');
    if (!svg) return false;
    const targetGroups = getAllFeatureLegendGroups(svg);
    if (targetGroups.length === 0) return false;
    if (results.value.length > 1 && Array.isArray(from) && !describesMountedLegend(svg, from)) {
      return projectLegendChange(targetGroups, legendEntries.value || [], from, entryOwners);
    }

    const desiredEntries = [];
    let entryColorStateChanged = false;
    const seenCaptions = new Set();
    (Array.isArray(legendEntries.value) ? legendEntries.value : []).forEach((entry) => {
      const caption = String(entry?.caption || '').trim();
      if (!caption || seenCaptions.has(caption)) return;
      seenCaptions.add(caption);
      const color = String(legendColorOverrides[caption] || entry?.color || '#cccccc');
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

    if (entryColorStateChanged) legendEntries.value = desiredEntries;
    if (!changed) return entryColorStateChanged;
    restoredCaptions.forEach((caption) => retiredEntryTemplates.delete(caption));
    onLegendGeometryChanged();
    persistLegendReconciliation();
    return true;
  };

  // The legend rows the specific-color rules draw, prepared on a disposable copy
  // of the mounted legend (R13). The rule owner runs its transition and then
  // `apply` in one History step, while `isCurrent` holds; no rule mutation runs
  // here. Resolves to false when the Result or the caller became stale.
  const prepareFileLegendEntries = async (intents, { previousFileIntents = [], isCurrent = () => true } = {}) => {
    const mountedSvg = svgContainer.value?.querySelector('svg');
    const svg = mountedSvg?.cloneNode(true);
    const targetGroups = svg ? getAllFeatureLegendGroups(svg) : [];
    if (targetGroups.length === 0) {
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

    if (!isCurrent() || svgContainer.value.querySelector('svg') !== mountedSvg) return false;
    // The copy's rows are laid out once they are mounted: `apply` hands the
    // Legend to the layout owner (zero shift).
    return {
      diff,
      isCurrent: () => svgContainer.value.querySelector('svg') === mountedSvg,
      // Mounted geometry and the Result commit synchronously inside History.
      apply: () => {
        const mountedLegend = mountedSvg.getElementById('legend');
        const candidateLegend = svg.getElementById('legend');
        if (mountedLegend && candidateLegend) mountedLegend.replaceWith(candidateLegend);
        onLegendGeometryChanged();
        commitActiveResultEdit?.('legend-file-sync');
        extractLegendEntries();
      }
    };
  };

  /**
   * Read the Legend of a mounted Result: its entries in visual order and the
   * captions the renderer generated (not the editor's direct entries). The one
   * reader of a mounted Legend; `previousEntries` keep stroke, feature ids,
   * and the generated caption of a renamed entry.
   * @param {SVGSVGElement} svg
   * @param {any[]} previousEntries
   */
  const readMountedLegend = (svg, previousEntries) => {
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
      ));
      const showStroke = existingEntry?.showStroke || false;
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

    return { entries: visuallySortedEntries, generatedCaptions };
  };

  // The inventory after a draw: while an edited order is replayed, the
  // surviving categories keep their default order and new ones follow it;
  // otherwise the Result shows the renderer's order, which becomes the default
  // order. Deleted entries keep their place.
  /**
   * @param {{ inventory: string[], orderEdited: boolean, generatedCaptions: Set<string>, rendered: string[] }} drawn
   */
  const drawnInventory = ({ inventory, orderEdited, generatedCaptions, rendered }) => {
    const retainedCaptions = new Set([
      ...generatedCaptions,
      ...deletedLegendEntries.value.map(entry => entry.originalCaption || entry.caption)
    ]);
    const surviving = inventory.filter(caption => retainedCaptions.has(caption));
    return [...new Set(orderEdited ? [...surviving, ...rendered] : [...rendered, ...surviving])];
  };
  const renderedCaptions = ({ entries, generatedCaptions }) => entries
    .map(entry => entry.originalCaption)
    .filter(caption => generatedCaptions.has(caption));

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
    pruneResultInventories(liveResultIdentities);
    const stored = inventoryByResult.get(identity);
    if (!identity || stored) return stored || [];
    const mounted = readMountedLegend(svg, legendEntries.value || []);
    const inventory = mounted
      ? drawnInventory({
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
    if (!svgContainer.value) {
      legendEntries.value = [];
      return;
    }

    const svg = svgContainer.value.querySelector('svg');
    if (!svg) {
      legendEntries.value = [];
      return;
    }

    const mounted = readMountedLegend(svg, legendEntries.value || []);
    if (!mounted) {
      legendEntries.value = [];
      return;
    }
    const { entries: visuallySortedEntries, generatedCaptions } = mounted;

    // The entries shown before this extraction decide, as in Generate, whether
    // the displayed Result replays an edited order (D-08).
    const orderEdited = isLegendOrderEdited(legendEntries.value, originalLegendOrder.value);
    legendEntries.value = visuallySortedEntries;

    const identity = readActiveResultIdentity?.() || '';
    if (replaceGeneratedInventory || originalLegendOrder.value.length === 0) {
      // Keep deletion intent; live editor extraction alone must not advance
      // the accepted generated inventory.
      originalLegendOrder.value = drawnInventory({
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

  const updateLegendEntryColor = (idx, newColor) => {
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    if (!svgContainer.value) return false;
    const svg = svgContainer.value.querySelector('svg');
    if (!svg) return false;

    const entry = legendEntries.value[idx];
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
    legendColorOverrides[caption] = newColor;
    if (changed) persistLegendReconciliation();
    return true;
  };

  const updateLegendEntryCaption = (idx, newCaption) => {
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
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
      const result = results.value[selectedResultIndex.value];
      setTimeout(() => {
        if (!state.sessionOperationAvailability?.()
          && results.value[selectedResultIndex.value] === result) extractLegendEntries();
      }, 100);
    }
  };

  const deleteLegendEntry = (idx) => {
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    const entry = legendEntries.value[idx];
    if (!entry) return;

    deletedLegendEntries.value.push({ ...entry });

    removeLegendEntry(entry.caption);
    extractLegendEntries();
  };

  const restoreDeletedLegendEntries = () => {
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    if (deletedLegendEntries.value.length === 0) return;

    for (const entry of deletedLegendEntries.value) {
      addLegendEntry(entry.caption, entry.color);
    }
    deletedLegendEntries.value = [];
    extractLegendEntries();
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
    removeLegendEntry,
    reconcileLegendEntries,
    restoreDeletedLegendEntries,
    setLegendGeometryChangedHandler,
    updateLegendEntryCaption,
    updateLegendEntryColor,
    updateLegendEntryColorByCaption
  };
};
