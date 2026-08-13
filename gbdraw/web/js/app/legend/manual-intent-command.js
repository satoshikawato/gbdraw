import {
  admitProjectedSvgResult,
  getCommittedSvgContent,
  getCommittedSvgResultMetadata,
  parseCommittedSvgResultRoot
} from '../../services/svg-result-ingestion.js';
import { catalogResultKey } from '../../services/feature-catalog.js';
import { styleFingerprint } from '../../services/style-revision.js';
import { cloneJsonData as cloneJson } from '../../services/json-clone.js';
import { normalizeRuleColor } from '../feature-editor/fill-view-model.js';
import {
  normalizeFeatureStrokeColor,
  normalizeFeatureStrokeWidth
} from '../feature-editor/stroke-view-model.js';
import {
  buildStyleSnapshotCommand,
  exactJsonByteLength,
  immutableStyleCommandSnapshot
} from '../feature-editor/style-snapshot-command.js';
import {
  deriveUsedFeatureFillGroupsByResult,
  prepareFeatureFillLegendProjections
} from './feature-fill-projection.js';
import { getAllFeatureLegendGroups, parseTransformXY } from './utils.js';

const MANUAL_OWNER = 'manual';
const ORDER_OWNER = 'legend-order';
const text = (value) => String(value ?? '').trim();
const captionKey = (value) => text(value).toLocaleLowerCase();

const readRef = (value) => (
  value && typeof value === 'object' && 'value' in value ? value.value : value
);

const writeRef = (target, value, label) => {
  if (!target || typeof target !== 'object' || !('value' in target)) {
    throw new Error(`Manual legend command requires a writable ${label} ref.`);
  }
  target.value = value;
};

const isPromiseLike = (value) => Boolean(value && typeof value.then === 'function');

const featureOwned = (entry) => (
  (Array.isArray(entry?.featureIds) && entry.featureIds.length > 0)
  || text(entry?.owner).toLowerCase().startsWith('feature')
);

const normalizeManualEntry = (entry, index) => {
  if (!entry || typeof entry !== 'object' || Array.isArray(entry)) {
    throw new Error(`Manual legend entry ${index + 1} is invalid.`);
  }
  if (featureOwned(entry)) {
    throw new Error('Feature-derived legend entries must use the Feature fill command.');
  }
  const caption = text(entry.caption);
  const color = normalizeRuleColor(entry.color);
  const rawStrokeColor = entry.strokeColor ?? entry.stroke_color;
  const rawStrokeWidth = entry.strokeWidth ?? entry.stroke_width;
  const strokeColor = rawStrokeColor === null || rawStrokeColor === undefined || rawStrokeColor === ''
    ? null
    : normalizeFeatureStrokeColor(rawStrokeColor);
  const strokeWidth = rawStrokeWidth === null || rawStrokeWidth === undefined || rawStrokeWidth === ''
    ? null
    : normalizeFeatureStrokeWidth(rawStrokeWidth);
  if (!caption) throw new Error(`Manual legend entry ${index + 1} has no caption.`);
  if (!color || color === 'none') {
    throw new Error(`Manual legend entry "${caption}" has an invalid color.`);
  }
  if (rawStrokeColor !== null && rawStrokeColor !== undefined && rawStrokeColor !== '' && !strokeColor) {
    throw new Error(`Manual legend entry "${caption}" has an invalid stroke color.`);
  }
  if (rawStrokeWidth !== null && rawStrokeWidth !== undefined && rawStrokeWidth !== '' && strokeWidth === null) {
    throw new Error(`Manual legend entry "${caption}" has an invalid stroke width.`);
  }
  return {
    ...cloneJson(entry),
    caption,
    originalCaption: text(entry.originalCaption) || caption,
    color,
    strokeColor,
    strokeWidth,
    owner: MANUAL_OWNER,
    featureIds: []
  };
};

export const normalizeManualLegendEntries = (entries = []) => {
  if (!Array.isArray(entries)) throw new Error('Manual legend intent must be an array.');
  const seen = new Set();
  return entries.map((entry, index) => {
    const normalized = normalizeManualEntry(entry, index);
    const key = captionKey(normalized.caption);
    if (seen.has(key)) {
      throw new Error(`Manual legend caption "${normalized.caption}" is duplicated.`);
    }
    seen.add(key);
    return normalized;
  });
};

/** Normalize the document-global caption precedence without inventing Result-local rows. */
export const normalizeLegendOrderIntent = (order = []) => {
  if (!Array.isArray(order)) throw new Error('Legend order intent must be an array.');
  const seen = new Set();
  const normalized = [];
  order.forEach((value) => {
    const caption = text(value);
    const key = captionKey(caption);
    if (!caption || seen.has(key)) return;
    seen.add(key);
    normalized.push(caption);
  });
  return immutableStyleCommandSnapshot(normalized);
};

const nextLegendOrderIntent = (currentOrder, intent) => {
  const current = [...normalizeLegendOrderIntent(currentOrder)];
  const kind = normalizedIntentKind(intent);
  if (kind === 'order') return normalizeLegendOrderIntent(intent?.order ?? intent?.captionOrder);
  if (current.length === 0) return immutableStyleCommandSnapshot(current);
  const caption = text(intent?.caption ?? intent?.oldCaption ?? intent?.entry?.caption);
  const key = captionKey(caption);
  if (kind === 'add') {
    const added = text(intent?.caption ?? intent?.entry?.caption);
    if (added && !current.some((entry) => captionKey(entry) === captionKey(added))) current.push(added);
  } else if (kind === 'remove') {
    return immutableStyleCommandSnapshot(current.filter((entry) => captionKey(entry) !== key));
  } else if (kind === 'rename') {
    const renamed = text(intent?.newCaption ?? intent?.nextCaption);
    return immutableStyleCommandSnapshot(current.map((entry) => (
      captionKey(entry) === key ? renamed : entry
    )).filter(Boolean));
  }
  return normalizeLegendOrderIntent(current);
};

const normalizedIntentKind = (intent) => text(intent?.kind ?? intent?.action).toLowerCase();

const assertedManualTarget = ({ entries, caption, selectedEntries, featureCaptionKeys }) => {
  const key = captionKey(caption);
  const selected = (Array.isArray(selectedEntries) ? selectedEntries : [])
    .find((entry) => captionKey(entry?.caption) === key);
  if (featureOwned(selected) || featureCaptionKeys.has(key)) {
    throw new Error('Feature-derived legend entries must use the Feature fill command.');
  }
  const index = entries.findIndex((entry) => captionKey(entry.caption) === key);
  if (index < 0) throw new Error(`Manual legend entry "${text(caption)}" is unavailable.`);
  return index;
};

/** Apply one document-global manual-only legend intent without mutating its inputs. */
export const applyManualLegendIntent = ({
  entries = [],
  intent,
  selectedEntries = [],
  featureCaptionKeys = new Set()
} = {}) => {
  if (featureOwned(intent?.entry) || featureOwned(intent?.targetEntry)) {
    throw new Error('Feature-derived legend entries must use the Feature fill command.');
  }
  const current = normalizeManualLegendEntries(entries);
  const next = current.map((entry) => cloneJson(entry));
  const usedFeatureCaptions = featureCaptionKeys instanceof Set
    ? featureCaptionKeys
    : new Set(Array.from(featureCaptionKeys || [], captionKey));
  const kind = normalizedIntentKind(intent);
  const caption = text(intent?.caption ?? intent?.oldCaption ?? intent?.entry?.caption);

  if (kind === 'add') {
    const rawEntry = {
      ...(intent?.entry && typeof intent.entry === 'object' ? cloneJson(intent.entry) : {}),
      caption,
      color: intent?.color ?? intent?.entry?.color
    };
    const added = normalizeManualEntry(rawEntry, next.length);
    const key = captionKey(added.caption);
    if (usedFeatureCaptions.has(key)) {
      throw new Error('Feature-derived legend entries must use the Feature fill command.');
    }
    if (next.some((entry) => captionKey(entry.caption) === key)) {
      throw new Error(`Manual legend caption "${added.caption}" already exists.`);
    }
    next.push(added);
    return immutableStyleCommandSnapshot(next);
  }

  if (!['remove', 'rename', 'color', 'stroke'].includes(kind)) {
    throw new Error(`Unsupported manual legend intent: ${kind || '(missing)'}`);
  }
  const index = assertedManualTarget({
    entries: next,
    caption,
    selectedEntries,
    featureCaptionKeys: usedFeatureCaptions
  });

  if (kind === 'remove') {
    next.splice(index, 1);
    return immutableStyleCommandSnapshot(next);
  }

  if (kind === 'rename') {
    const newCaption = text(intent?.newCaption ?? intent?.nextCaption);
    if (!newCaption) throw new Error('A renamed manual legend entry requires a caption.');
    const newKey = captionKey(newCaption);
    if (usedFeatureCaptions.has(newKey)) {
      throw new Error('Feature-derived legend entries must use the Feature fill command.');
    }
    if (next.some((entry, candidateIndex) => (
      candidateIndex !== index && captionKey(entry.caption) === newKey
    ))) {
      throw new Error(`Manual legend caption "${newCaption}" already exists.`);
    }
    next[index] = { ...next[index], caption: newCaption };
    return immutableStyleCommandSnapshot(next);
  }

  if (kind === 'color') {
    const color = normalizeRuleColor(intent?.color ?? intent?.newColor);
    if (!color || color === 'none') {
      throw new Error(`Manual legend entry "${caption}" has an invalid color.`);
    }
    next[index] = { ...next[index], color };
    return immutableStyleCommandSnapshot(next);
  }

  const requested = intent?.value && typeof intent.value === 'object' ? intent.value : intent;
  const inherit = requested?.kind === 'inherit' || requested?.inherit === true;
  if (inherit) {
    next[index] = { ...next[index], strokeColor: null, strokeWidth: null };
    return immutableStyleCommandSnapshot(next);
  }
  const hasColor = requested?.strokeColor !== undefined || requested?.color !== undefined;
  const hasWidth = requested?.strokeWidth !== undefined || requested?.width !== undefined;
  const strokeColor = hasColor
    ? normalizeFeatureStrokeColor(requested.strokeColor ?? requested.color)
    : next[index].strokeColor;
  const strokeWidth = hasWidth
    ? normalizeFeatureStrokeWidth(requested.strokeWidth ?? requested.width)
    : next[index].strokeWidth;
  if ((hasColor && !strokeColor) || (hasWidth && strokeWidth === null)) {
    throw new Error(`Manual legend entry "${caption}" has an invalid stroke.`);
  }
  next[index] = { ...next[index], strokeColor, strokeWidth };
  return immutableStyleCommandSnapshot(next);
};

const resultBindings = (results, catalog) => {
  if (!Array.isArray(results) || !Array.isArray(catalog?.items)) {
    throw new Error('Manual legend projection requires Results and a feature catalogue.');
  }
  if (results.length !== catalog.items.length || results.length === 0) {
    throw new Error('Manual legend Result topology does not match the feature catalogue.');
  }
  const keys = new Set();
  return catalog.items.map((item, resultIndex) => {
    const resultKey = catalogResultKey(item);
    const result = results[resultIndex];
    if (
      !resultKey
      || keys.has(resultKey)
      || !result
      || (item?.resultIndex !== undefined && Number(item.resultIndex) !== resultIndex)
      || (text(item?.resultName) && text(item.resultName) !== text(result?.name))
    ) {
      throw new Error(`Manual legend Result identity is invalid at index ${resultIndex}.`);
    }
    keys.add(resultKey);
    return { item, result, resultIndex, resultKey };
  });
};

const existingEntriesMap = ({ catalog, selectedResultIndex, selectedLegendEntries, supplied }) => {
  if (supplied instanceof Map) return supplied;
  if (supplied && typeof supplied === 'object') return supplied;
  const entries = {};
  (Array.isArray(catalog?.items) ? catalog.items : []).forEach((item, resultIndex) => {
    const resultKey = catalogResultKey(item);
    const catalogEntries = item?.legendEntries ?? item?.derivedLegendEntries;
    entries[resultKey] = resultIndex === selectedResultIndex
      ? cloneJson(selectedLegendEntries || [])
      : cloneJson(Array.isArray(catalogEntries) ? catalogEntries : []);
  });
  return entries;
};

const normalizedMounted = (mounted) => {
  const rawIndex = mounted?.resultIndex;
  const resultIndex = rawIndex !== null
    && rawIndex !== undefined
    && Number.isInteger(Number(rawIndex))
    && Number(rawIndex) >= 0
    ? Number(rawIndex)
    : null;
  return {
    ...(mounted || {}),
    resultIndex,
    resultKey: text(mounted?.resultKey),
    svg: mounted?.svg || null
  };
};

const sameMountedOwnership = (left, right) => (
  (left?.resultIndex ?? null) === (right?.resultIndex ?? null)
  && text(left?.resultKey) === text(right?.resultKey)
  && (left?.svg || null) === (right?.svg || null)
);

const currentMounted = (getMountedContext) => normalizedMounted(
  typeof getMountedContext === 'function' ? getMountedContext() : null
);

const mountedIsValid = (mounted, bindings) => {
  if (mounted.resultIndex === null) return !mounted.resultKey && !mounted.svg;
  const binding = bindings[mounted.resultIndex];
  return Boolean(
    binding
    && binding.resultKey === mounted.resultKey
    && mounted.svg?.cloneNode
  );
};

const resultContent = (result, committedContent = getCommittedSvgContent) => {
  const committed = committedContent(result);
  return committed === null || committed === undefined
    ? String(result?.content ?? '')
    : String(committed);
};

const directLegendEntryGroups = (targetGroup) => {
  const direct = Array.from(targetGroup?.children || []).filter((child) => (
    text(child?.localName ?? child?.tagName).toLowerCase() === 'g'
    && child.hasAttribute?.('data-legend-key')
  ));
  return direct.length > 0
    ? direct
    : Array.from(targetGroup?.querySelectorAll?.('g[data-legend-key]') || []);
};

const legendEntrySwatch = (group) => Array.from(group?.querySelectorAll?.('path') || [])
  .find((path) => {
    const fill = text(path?.getAttribute?.('fill')).toLowerCase();
    return fill && fill !== 'none' && !fill.startsWith('url(');
  }) || null;

const legendEntryAnchor = (group) => {
  const groupOffset = parseTransformXY(group?.getAttribute?.('transform'));
  const child = group?.querySelector?.('text') || legendEntrySwatch(group);
  const childOffset = parseTransformXY(child?.getAttribute?.('transform'));
  return { x: groupOffset.x + childOffset.x, y: groupOffset.y + childOffset.y };
};

const moveLegendEntryToAnchor = (group, anchor) => {
  const current = legendEntryAnchor(group);
  const deltaX = anchor.x - current.x;
  const deltaY = anchor.y - current.y;
  if (Math.abs(deltaX) < 1e-9 && Math.abs(deltaY) < 1e-9) return false;
  const transformed = group?.getAttribute?.('transform')
    ? [group]
    : [legendEntrySwatch(group), group?.querySelector?.('text')]
      .filter((node, index, values) => (
        values.indexOf(node) === index && node?.getAttribute?.('transform')
      ));
  if (transformed.length === 0) {
    group.setAttribute?.('transform', `translate(${deltaX}, ${deltaY})`);
  } else {
    transformed.forEach((node) => {
      const position = parseTransformXY(node.getAttribute('transform'));
      node.setAttribute('transform', `translate(${position.x + deltaX}, ${position.y + deltaY})`);
    });
  }
  return true;
};

/** Apply global caption precedence to only the captions present in one Result. */
export const orderedLegendCaptions = (currentCaptions = [], orderIntent = []) => {
  const current = normalizeLegendOrderIntent(currentCaptions);
  const byKey = new Map(current.map((caption) => [captionKey(caption), caption]));
  const ordered = [];
  normalizeLegendOrderIntent(orderIntent).forEach((caption) => {
    const present = byKey.get(captionKey(caption));
    if (!present) return;
    ordered.push(present);
    byKey.delete(captionKey(caption));
  });
  current.forEach((caption) => {
    const key = captionKey(caption);
    if (!byKey.has(key)) return;
    ordered.push(byKey.get(key));
    byKey.delete(key);
  });
  return immutableStyleCommandSnapshot(ordered);
};

const visualLegendEntryGroups = (targetGroup) => directLegendEntryGroups(targetGroup)
  .map((group) => ({ group, anchor: legendEntryAnchor(group) }))
  .sort((left, right) => left.anchor.y - right.anchor.y || left.anchor.x - right.anchor.x);

/** Mutate a detached SVG clone while preserving that Result's existing legend slots. */
export const applyLegendOrderToSvg = (svg, orderIntent = []) => {
  if (!svg?.cloneNode) throw new Error('Legend order projection requires an SVG root.');
  const counters = { targetGroups: 0, movedEntries: 0 };
  getAllFeatureLegendGroups(svg).forEach((targetGroup) => {
    const visual = visualLegendEntryGroups(targetGroup);
    if (visual.length < 2) return;
    const byCaption = new Map();
    visual.forEach(({ group }) => {
      const caption = text(group.getAttribute?.('data-legend-key'));
      const key = captionKey(caption);
      if (!key || byCaption.has(key)) {
        throw new Error('Legend order projection found a missing or duplicate caption.');
      }
      byCaption.set(key, group);
    });
    const desiredCaptions = orderedLegendCaptions(
      visual.map(({ group }) => group.getAttribute?.('data-legend-key')),
      orderIntent
    );
    const slots = visual.map(({ anchor }) => anchor);
    desiredCaptions.forEach((caption, index) => {
      const group = byCaption.get(captionKey(caption));
      if (!group) throw new Error(`Legend entry "${caption}" could not be ordered.`);
      if (moveLegendEntryToAnchor(group, slots[index])) counters.movedEntries += 1;
      targetGroup.appendChild?.(group);
    });
    counters.targetGroups += 1;
  });
  return Object.freeze(counters);
};

const orderForResult = (orders, resultKey, fallback) => {
  if (orders instanceof Map && orders.has(resultKey)) return orders.get(resultKey);
  if (orders && typeof orders === 'object' && Array.isArray(orders[resultKey])) {
    return orders[resultKey];
  }
  return fallback;
};

const orderSelectedLegendEntries = (entries, orderIntent) => {
  const current = Array.isArray(entries) ? entries : [];
  const byCaption = new Map(current.map((entry) => [captionKey(entry?.caption), entry]));
  const ordered = orderedLegendCaptions(current.map((entry) => entry?.caption), orderIntent);
  const slots = [...current]
    .map((entry) => ({ xPos: entry?.xPos, yPos: entry?.yPos }))
    .sort((left, right) => (
      (Number(left.yPos) || 0) - (Number(right.yPos) || 0)
      || (Number(left.xPos) || 0) - (Number(right.xPos) || 0)
    ));
  return ordered.map((caption, index) => ({
    ...cloneJson(byCaption.get(captionKey(caption))),
    ...(slots[index] || {})
  }));
};

/** Stage one document-global order intent across every admitted Result. */
export const prepareLegendOrderResultProjection = async ({
  results,
  catalog,
  legendOrderIntent = [],
  orderByResult = null,
  selectedResultIndex = null,
  selectedLegendEntries = [],
  mounted = null,
  parser = globalThis.DOMParser || globalThis.window?.DOMParser,
  parseResultRoot = parseCommittedSvgResultRoot,
  admitResult = admitProjectedSvgResult,
  committedContent = getCommittedSvgContent,
  committedMetadata = getCommittedSvgResultMetadata
} = {}) => {
  const bindings = resultBindings(results, catalog);
  const active = normalizedMounted(mounted);
  if (!mountedIsValid(active, bindings)) {
    throw new Error('The mounted Result changed before legend order projection.');
  }
  const selectedIndex = Number.isInteger(Number(selectedResultIndex))
    && Number(selectedResultIndex) >= 0
    && Number(selectedResultIndex) < bindings.length
    ? Number(selectedResultIndex)
    : null;
  const candidates = new Map();
  const preparedSvgByResultKey = new Map();
  const admissionMetadataByResultKey = {};
  const sourceOrderByResult = {};
  let mountedSources = 0;
  let detachedPasses = 0;
  let changedResults = 0;
  let movedEntries = 0;
  let targetGroups = 0;

  bindings.forEach((binding) => {
    const source = binding.resultIndex === active.resultIndex
      ? active.svg
      : parseResultRoot(binding.result, parser);
    if (binding.resultIndex === active.resultIndex) mountedSources += 1;
    else detachedPasses += 1;
    const primary = getAllFeatureLegendGroups(source)[0] || null;
    sourceOrderByResult[binding.resultKey] = primary
      ? visualLegendEntryGroups(primary).map(({ group }) => text(group.getAttribute?.('data-legend-key')))
      : [];
    const svg = source.cloneNode(true);
    const requestedOrder = orderForResult(orderByResult, binding.resultKey, legendOrderIntent);
    const applied = applyLegendOrderToSvg(svg, requestedOrder);
    movedEntries += applied.movedEntries;
    targetGroups += applied.targetGroups;
    const candidate = admitResult(binding.result, svg);
    if (!candidate || typeof candidate !== 'object') {
      throw new Error(`Legend order projection failed admission for Result ${binding.resultKey}.`);
    }
    const selected = resultContent(candidate, committedContent) === resultContent(binding.result, committedContent)
      ? binding.result
      : candidate;
    if (selected !== binding.result) changedResults += 1;
    candidates.set(binding.resultIndex, selected);
    preparedSvgByResultKey.set(binding.resultKey, svg);
    admissionMetadataByResultKey[binding.resultKey] = Object.freeze({
      before: committedMetadata(binding.result) || null,
      after: committedMetadata(selected) || null
    });
  });

  const selectedBinding = selectedIndex === null ? null : bindings[selectedIndex];
  const selectedOrder = selectedBinding
    ? orderForResult(orderByResult, selectedBinding.resultKey, legendOrderIntent)
    : legendOrderIntent;
  return Object.freeze({
    previousResults: results,
    nextResults: results.map((result, resultIndex) => candidates.get(resultIndex) || result),
    affectedResultKeys: Object.freeze(bindings.map((binding) => binding.resultKey)),
    admissionMetadataByResultKey: Object.freeze(admissionMetadataByResultKey),
    preparedSvgByResultKey,
    preparedMountedSvg: active.resultIndex === null
      ? null
      : preparedSvgByResultKey.get(active.resultKey) || null,
    selectedLegendEntries: orderSelectedLegendEntries(selectedLegendEntries, selectedOrder),
    sourceOrderByResult: Object.freeze(cloneJson(sourceOrderByResult)),
    counters: Object.freeze({
      affectedResults: bindings.length,
      mountedSerializations: mountedSources,
      detachedPasses,
      changedResults,
      movedEntries,
      targetGroups,
      measurements: 0,
      additions: 0,
      renames: 0,
      removals: 0,
      colorUpdates: 0
    })
  });
};

/**
 * Stage the manual intent into every Result. The default path uses only browser
 * SVG parsing, cloning, measurement, and admitted-result projection.
 */
export const prepareManualLegendResultProjection = async ({
  results,
  catalog,
  rules = [],
  paletteColors = {},
  manualLegendEntries = [],
  legendOrderIntent = [],
  selectedResultIndex = null,
  selectedLegendEntries = [],
  existingEntriesByResult = null,
  mounted = null,
  svgDefaultColor = '#cccccc',
  parser = globalThis.DOMParser || globalThis.window?.DOMParser,
  legendProjectionOptions = {},
  parseResultRoot = parseCommittedSvgResultRoot,
  admitResult = admitProjectedSvgResult,
  committedContent = getCommittedSvgContent,
  committedMetadata = getCommittedSvgResultMetadata,
  prepareLegends = prepareFeatureFillLegendProjections
} = {}) => {
  const bindings = resultBindings(results, catalog);
  const active = normalizedMounted(mounted);
  if (!mountedIsValid(active, bindings)) {
    throw new Error('The mounted Result changed before manual legend projection.');
  }
  const selectedIndex = Number.isInteger(Number(selectedResultIndex))
    && Number(selectedResultIndex) >= 0
    && Number(selectedResultIndex) < bindings.length
    ? Number(selectedResultIndex)
    : null;
  const sourcesByResult = new Map();
  let mountedSources = 0;
  let detachedPasses = 0;
  bindings.forEach((binding) => {
    if (binding.resultIndex === active.resultIndex) {
      sourcesByResult.set(binding.resultKey, active.svg);
      mountedSources += 1;
      return;
    }
    sourcesByResult.set(binding.resultKey, parseResultRoot(binding.result, parser));
    detachedPasses += 1;
  });

  const existing = existingEntriesMap({
    catalog,
    selectedResultIndex: selectedIndex,
    selectedLegendEntries,
    supplied: existingEntriesByResult
  });
  const derived = deriveUsedFeatureFillGroupsByResult({
    catalog,
    rules,
    paletteColors,
    manualLegendEntries,
    legendOrder: legendOrderIntent,
    existingEntriesByResult: existing,
    svgDefaultColor
  });
  if (derived.projections.length !== bindings.length) {
    throw new Error('Manual legend derivation did not cover every Result.');
  }
  const prepared = await prepareLegends({
    sourcesByResult,
    projections: derived.projections,
    ...legendProjectionOptions
  });
  const staged = prepared?.staged;
  if (!(staged instanceof Map)) {
    throw new Error('Manual legend preparation did not return Result-local staged roots.');
  }

  const candidates = new Map();
  const admissionMetadataByResultKey = {};
  const preparedSvgByResultKey = new Map();
  let changedResults = 0;
  bindings.forEach((binding) => {
    const stagedEntry = staged.get(binding.resultKey);
    const svg = stagedEntry?.svg || stagedEntry;
    if (!svg?.cloneNode) {
      throw new Error(`Manual legend projection is incomplete for Result ${binding.resultKey}.`);
    }
    const candidate = admitResult(binding.result, svg);
    if (!candidate || typeof candidate !== 'object') {
      throw new Error(`Manual legend projection failed admission for Result ${binding.resultKey}.`);
    }
    const selected = resultContent(candidate, committedContent) === resultContent(binding.result, committedContent)
      ? binding.result
      : candidate;
    if (selected !== binding.result) changedResults += 1;
    candidates.set(binding.resultIndex, selected);
    preparedSvgByResultKey.set(binding.resultKey, svg);
    admissionMetadataByResultKey[binding.resultKey] = Object.freeze({
      before: committedMetadata(binding.result) || null,
      after: committedMetadata(selected) || null
    });
  });

  const selectedProjection = selectedIndex === null
    ? null
    : derived.projections.find((entry) => entry.resultIndex === selectedIndex);
  return Object.freeze({
    previousResults: results,
    nextResults: results.map((result, resultIndex) => candidates.get(resultIndex) || result),
    affectedResultKeys: Object.freeze(bindings.map((binding) => binding.resultKey)),
    admissionMetadataByResultKey: Object.freeze(admissionMetadataByResultKey),
    preparedSvgByResultKey,
    preparedMountedSvg: active.resultIndex === null
      ? null
      : preparedSvgByResultKey.get(active.resultKey) || null,
    selectedLegendEntries: cloneJson(selectedProjection?.entries ?? selectedLegendEntries ?? []),
    counters: Object.freeze({
      affectedResults: bindings.length,
      mountedSerializations: mountedSources,
      detachedPasses,
      changedResults,
      measurements: Number(prepared?.counters?.measurements || 0),
      additions: Number(prepared?.counters?.additions || 0),
      renames: Number(prepared?.counters?.renames || 0),
      removals: Number(prepared?.counters?.removals || 0),
      colorUpdates: Number(prepared?.counters?.colorUpdates || 0)
    })
  });
};

const featureCaptionKeys = ({ catalog, rules, paletteColors, svgDefaultColor }) => {
  const projection = deriveUsedFeatureFillGroupsByResult({
    catalog,
    rules,
    paletteColors,
    manualLegendEntries: [],
    svgDefaultColor
  });
  const keys = new Set();
  projection.projections.forEach((result) => result.entries.forEach((entry) => {
    if (Array.isArray(entry.featureIds) && entry.featureIds.length > 0) {
      keys.add(captionKey(entry.caption));
    }
  }));
  return keys;
};

const paletteSnapshot = (state) => ({
  appliedPaletteName: text(readRef(state?.appliedPaletteName)) || 'default',
  appliedPaletteColors: cloneJson(readRef(state?.appliedPaletteColors) || {})
});

const currentStyleFingerprint = (state) => styleFingerprint({
  rules: state?.manualSpecificRules || [],
  ...paletteSnapshot(state)
});

const ledgerIsCurrent = (ledger, resultKeys, fingerprint) => resultKeys.every(
  (resultKey) => text(ledger?.[resultKey]) === fingerprint
);

const manualSignature = (entries) => JSON.stringify(normalizeManualLegendEntries(entries));

const resultSignature = (result) => String(result?.content ?? '');

const projectedResultsMatch = (published, projected) => (
  Array.isArray(published)
  && Array.isArray(projected)
  && published.length === projected.length
  && published.every((result, index) => (
    text(result?.name) === text(projected[index]?.name)
    && resultSignature(result) === resultSignature(projected[index])
  ))
);

const copyAttributes = (target, source) => {
  const targetNames = Array.from(target?.attributes || []).map((attribute) => (
    Array.isArray(attribute) ? attribute[0] : attribute.name
  )).filter(Boolean);
  targetNames.forEach((name) => target.removeAttribute?.(name));
  Array.from(source?.attributes || []).forEach((attribute) => {
    const name = Array.isArray(attribute) ? attribute[0] : attribute.name;
    const value = Array.isArray(attribute) ? attribute[1] : attribute.value;
    if (name) target.setAttribute?.(name, value);
  });
};

const replaceMountedRoot = (target, source) => {
  if (!target?.cloneNode || !source?.cloneNode) {
    throw new Error('Manual legend mounted projection requires SVG roots.');
  }
  copyAttributes(target, source);
  const children = Array.from(source.childNodes || source.children || [])
    .map((child) => child.cloneNode(true));
  if (typeof target.replaceChildren === 'function') {
    target.replaceChildren(...children);
  } else {
    Array.from(target.childNodes || target.children || []).forEach((child) => {
      child.remove?.();
      if (child.parentElement === target) target.removeChild?.(child);
    });
    children.forEach((child) => target.appendChild?.(child));
  }
  return true;
};

const captureExactState = ({ state, mounted, captureRuntimeState }) => {
  const runtime = typeof captureRuntimeState === 'function' ? captureRuntimeState() : null;
  if (isPromiseLike(runtime)) throw new Error('Manual legend runtime capture must be synchronous.');
  return {
    manualEntries: readRef(state.manualLegendEntries),
    orderIntent: readRef(state.legendOrderIntent),
    results: readRef(state.results),
    legendEntries: readRef(state.legendEntries),
    ledger: readRef(state.validatedStyleFingerprintByResultKey),
    revision: Number(readRef(state.semanticStyleRevision) || 0),
    fingerprint: text(readRef(state.semanticStyleFingerprint)),
    selectedResultIndex: readRef(state.selectedResultIndex),
    mounted,
    mountedSvg: mounted.svg?.cloneNode?.(true) || null,
    runtime
  };
};

const restoreExactState = ({
  state,
  snapshot,
  restoreMountedProjection,
  restoreRuntimeState,
  context,
  intent
}) => {
  if (snapshot.mounted?.svg && snapshot.mountedSvg) {
    const restored = typeof restoreMountedProjection === 'function'
      ? restoreMountedProjection({ snapshot, context, intent })
      : replaceMountedRoot(snapshot.mounted.svg, snapshot.mountedSvg);
    if (isPromiseLike(restored)) throw new Error('Manual legend mounted rollback must be synchronous.');
    if (restored === false) return false;
  }
  writeRef(state.manualLegendEntries, snapshot.manualEntries, 'manualLegendEntries');
  if (state.legendOrderIntent && snapshot.orderIntent !== undefined) {
    writeRef(state.legendOrderIntent, snapshot.orderIntent, 'legendOrderIntent');
  }
  writeRef(state.results, snapshot.results, 'results');
  writeRef(state.legendEntries, snapshot.legendEntries, 'legendEntries');
  writeRef(
    state.validatedStyleFingerprintByResultKey,
    snapshot.ledger,
    'validatedStyleFingerprintByResultKey'
  );
  writeRef(state.semanticStyleRevision, snapshot.revision, 'semanticStyleRevision');
  writeRef(state.semanticStyleFingerprint, snapshot.fingerprint, 'semanticStyleFingerprint');
  if (state.selectedResultIndex && snapshot.selectedResultIndex !== undefined) {
    writeRef(state.selectedResultIndex, snapshot.selectedResultIndex, 'selectedResultIndex');
  }
  if (typeof restoreRuntimeState === 'function') {
    const restored = restoreRuntimeState(snapshot.runtime);
    if (isPromiseLike(restored)) throw new Error('Manual legend runtime rollback must be synchronous.');
    if (restored === false) return false;
  }
  return true;
};

const requiredToken = (intent, key, fallback) => {
  const value = intent?.[key] ?? fallback;
  const number = Number(value);
  if (!Number.isSafeInteger(number) || number < 0) {
    throw new Error(`Manual legend intent requires ${key}.`);
  }
  return number;
};

/**
 * Build one History-compatible document-global manual legend command.
 *
 * Stable wiring surface: pass the intent tokens captured by the UI, state,
 * catalogue, and mounted/runtime adapters. All async Result/legend work finishes
 * before the synchronous commit begins.
 */
export const buildManualLegendIntentCommand = async ({
  intent,
  state,
  catalog,
  getMountedContext = null,
  prepareProjection = prepareManualLegendResultProjection,
  prepareOrderProjection = prepareLegendOrderResultProjection,
  commitMountedProjection = null,
  restoreMountedProjection = null,
  captureRuntimeState = null,
  restoreRuntimeState = null,
  getExistingEntriesByResult = null,
  reconcile = null,
  refreshPresentation = null,
  legendProjectionOptions = {},
  svgDefaultColor = '#cccccc',
  label = 'Change manual legend'
} = {}) => {
  if (
    !state
    || !catalog
    || typeof prepareProjection !== 'function'
    || typeof prepareOrderProjection !== 'function'
  ) {
    throw new Error('Manual legend command requires state, catalogue, and a projection owner.');
  }
  const intentKind = normalizedIntentKind(intent);
  const isOrderIntent = intentKind === 'order';
  if (isOrderIntent && !state.legendOrderIntent) {
    throw new Error('Legend order editing requires a writable document-global order ref.');
  }
  const expectedEpoch = requiredToken(intent, 'documentEpoch');
  const expectedGeneration = requiredToken(intent, 'resultGenerationKey');
  const expectedRevision = requiredToken(intent, 'semanticStyleRevision');
  const bindings = resultBindings(readRef(state.results), catalog);
  const resultKeys = bindings.map((binding) => binding.resultKey);
  const fingerprint = currentStyleFingerprint(state);
  const expectedFingerprint = text(intent?.styleFingerprint) || fingerprint;
  const initialMounted = currentMounted(getMountedContext);
  if (
    Number(readRef(state.documentEpoch) || 0) !== expectedEpoch
    || Number(readRef(state.resultGenerationKey) || 0) !== expectedGeneration
    || Number(readRef(state.semanticStyleRevision) || 0) !== expectedRevision
    || text(readRef(state.semanticStyleFingerprint)) !== expectedFingerprint
    || fingerprint !== expectedFingerprint
    || !ledgerIsCurrent(
      readRef(state.validatedStyleFingerprintByResultKey) || {},
      resultKeys,
      fingerprint
    )
    || !mountedIsValid(initialMounted, bindings)
    || (readRef(state.featureCatalog) && readRef(state.featureCatalog) !== catalog)
  ) {
    throw new Error('The manual legend intent became stale before command preparation.');
  }
  if (
    intent?.mountedResultKey !== undefined
    && text(intent.mountedResultKey) !== initialMounted.resultKey
  ) {
    throw new Error('The manual legend intent belongs to a different mounted Result.');
  }

  const selectedEntries = readRef(state.legendEntries) || [];
  const palette = paletteSnapshot(state);
  const usedFeatureCaptions = featureCaptionKeys({
    catalog,
    rules: state.manualSpecificRules || [],
    paletteColors: palette.appliedPaletteColors,
    svgDefaultColor
  });
  const beforeEntries = immutableStyleCommandSnapshot(
    normalizeManualLegendEntries(readRef(state.manualLegendEntries) || [])
  );
  const afterEntries = isOrderIntent
    ? beforeEntries
    : applyManualLegendIntent({
        entries: beforeEntries,
        intent,
        selectedEntries,
        featureCaptionKeys: usedFeatureCaptions
      });
  const beforeOrder = normalizeLegendOrderIntent(readRef(state.legendOrderIntent) || []);
  const afterOrder = nextLegendOrderIntent(beforeOrder, intent);
  const ledger = cloneJson(readRef(state.validatedStyleFingerprintByResultKey) || {});
  const before = immutableStyleCommandSnapshot({
    manualEntries: beforeEntries,
    legendOrder: beforeOrder,
    fingerprint,
    resultFingerprintByKey: ledger
  });
  const after = immutableStyleCommandSnapshot({
    manualEntries: afterEntries,
    legendOrder: afterOrder,
    fingerprint,
    resultFingerprintByKey: ledger
  });
  let transitionRevision = expectedRevision;
  let firstApplyPending = true;
  let implicitBeforeOrderByResult = null;

  const validateCurrent = ({ from, prepared = null, direction }) => {
    let currentBindings;
    try {
      currentBindings = resultBindings(readRef(state.results), catalog);
    } catch {
      return false;
    }
    if (
      Number(readRef(state.documentEpoch) || 0) !== expectedEpoch
      || Number(readRef(state.resultGenerationKey) || 0) !== expectedGeneration
      || Number(readRef(state.semanticStyleRevision) || 0) !== transitionRevision
      || currentStyleFingerprint(state) !== from.fingerprint
      || text(readRef(state.semanticStyleFingerprint)) !== from.fingerprint
      || manualSignature(readRef(state.manualLegendEntries) || []) !== manualSignature(from.manualEntries)
      || JSON.stringify(normalizeLegendOrderIntent(readRef(state.legendOrderIntent) || []))
        !== JSON.stringify(from.legendOrder)
      || !ledgerIsCurrent(
        readRef(state.validatedStyleFingerprintByResultKey) || {},
        resultKeys,
        from.fingerprint
      )
      || (readRef(state.featureCatalog) && readRef(state.featureCatalog) !== catalog)
    ) return false;
    const mounted = currentMounted(getMountedContext);
    if (!mountedIsValid(mounted, currentBindings)) return false;
    if (
      firstApplyPending
      && direction !== 'undo'
      && intent?.mountedResultKey !== undefined
      && text(intent.mountedResultKey) !== mounted.resultKey
    ) return false;
    if (prepared && (
      readRef(state.results) !== prepared.sourceResults
      || readRef(state.legendEntries) !== prepared.sourceLegendEntries
      || readRef(state.manualLegendEntries) !== prepared.sourceManualEntries
      || readRef(state.legendOrderIntent) !== prepared.sourceOrderIntent
      || Number(readRef(state.selectedResultIndex)) !== prepared.selectedResultIndex
      || !sameMountedOwnership(mounted, prepared.preflightMounted)
      || prepared.sourceResultRecords.some(({ resultIndex, result, content }) => (
        readRef(state.results)?.[resultIndex] !== result
        || resultSignature(result) !== content
      ))
    )) return false;
    return true;
  };

  const prepareTransition = async ({ from, to, direction }) => {
    const sourceResults = readRef(state.results);
    const sourceLegendEntries = readRef(state.legendEntries);
    const sourceManualEntries = readRef(state.manualLegendEntries);
    const sourceOrderIntent = readRef(state.legendOrderIntent);
    const selectedResultIndex = Number(readRef(state.selectedResultIndex));
    const mounted = currentMounted(getMountedContext);
    const suppliedExisting = typeof getExistingEntriesByResult === 'function'
      ? getExistingEntriesByResult({
          results: sourceResults,
          catalog,
          mounted,
          selectedResultIndex,
          selectedLegendEntries: sourceLegendEntries
        })
      : null;
    if (isPromiseLike(suppliedExisting)) {
      throw new Error('Manual legend existing-entry lookup must be synchronous.');
    }
    const projectionOwner = isOrderIntent ? prepareOrderProjection : prepareProjection;
    const projection = await projectionOwner({
      direction,
      intent,
      from,
      to,
      results: sourceResults,
      catalog,
      rules: state.manualSpecificRules || [],
      paletteColors: palette.appliedPaletteColors,
      manualLegendEntries: to.manualEntries,
      legendOrderIntent: to.legendOrder,
      orderByResult: isOrderIntent && to.legendOrder.length === 0
        ? implicitBeforeOrderByResult
        : null,
      selectedResultIndex,
      selectedLegendEntries: sourceLegendEntries,
      existingEntriesByResult: suppliedExisting,
      mounted,
      svgDefaultColor,
      legendProjectionOptions
    });
    if (
      isOrderIntent
      && direction !== 'undo'
      && !implicitBeforeOrderByResult
      && projection?.sourceOrderByResult
    ) {
      implicitBeforeOrderByResult = cloneJson(projection.sourceOrderByResult);
    }
    if (!projection || !Array.isArray(projection.nextResults)) {
      throw new Error('Manual legend projection did not prepare a complete Result array.');
    }
    if (projection.nextResults.length !== sourceResults.length) {
      throw new Error('Manual legend projection changed Result topology.');
    }
    projection.nextResults.forEach((result, resultIndex) => {
      if (!result || text(result.name) !== text(sourceResults[resultIndex]?.name)) {
        throw new Error(`Manual legend projection changed Result identity at index ${resultIndex}.`);
      }
    });
    const projectedKeys = Array.isArray(projection.affectedResultKeys)
      ? projection.affectedResultKeys.map(text)
      : resultKeys;
    if (
      projectedKeys.length !== resultKeys.length
      || resultKeys.some((key) => !projectedKeys.includes(key))
    ) {
      throw new Error('Manual legend projection did not cover every Result.');
    }
    const sourceResultRecords = sourceResults.map((result, resultIndex) => ({
      resultIndex,
      result,
      content: resultSignature(result)
    }));
    const retainedForHistory = {
      beforeResults: sourceResults,
      afterResults: projection.nextResults,
      admissionMetadataByResultKey: projection.admissionMetadataByResultKey || null
    };
    return {
      projection,
      sourceResults,
      sourceResultRecords,
      sourceLegendEntries,
      sourceManualEntries,
      sourceOrderIntent,
      selectedResultIndex,
      preflightMounted: mounted,
      selectedLegendEntries: cloneJson(
        projection.selectedLegendEntries ?? sourceLegendEntries ?? []
      ),
      retainedForHistory,
      retainedBytes: exactJsonByteLength(retainedForHistory),
      counters: {
        affectedResults: resultKeys.length,
        projectedResults: Number(projection.counters?.affectedResults || resultKeys.length),
        mountedSerializations: Number(projection.counters?.mountedSerializations || 0),
        detachedPasses: Number(projection.counters?.detachedPasses || 0),
        changedResults: Number(projection.counters?.changedResults || 0),
        resultArraySwaps: 1,
        legendMeasurements: Number(projection.counters?.measurements || 0),
        legendAdditions: Number(projection.counters?.additions || 0),
        legendRenames: Number(projection.counters?.renames || 0),
        legendRemovals: Number(projection.counters?.removals || 0),
        legendColorUpdates: Number(projection.counters?.colorUpdates || 0),
        legendOrderMoves: Number(projection.counters?.movedEntries || 0)
      }
    };
  };

  const commitTransition = ({ to, direction, prepared }) => {
    writeRef(
      state.manualLegendEntries,
      cloneJson(to.manualEntries),
      'manualLegendEntries'
    );
    if (state.legendOrderIntent) {
      writeRef(state.legendOrderIntent, cloneJson(to.legendOrder), 'legendOrderIntent');
    }
    if (prepared.preflightMounted.resultIndex !== null) {
      const committed = typeof commitMountedProjection === 'function'
        ? commitMountedProjection({
            direction,
            intent,
            mounted: prepared.preflightMounted,
            prepared: prepared.projection,
            to
          })
        : replaceMountedRoot(
            prepared.preflightMounted.svg,
            prepared.projection.preparedMountedSvg
          );
      if (isPromiseLike(committed)) {
        throw new Error('Manual legend mounted commit must be synchronous.');
      }
      if (committed === false) return false;
    }
    writeRef(state.legendEntries, cloneJson(prepared.selectedLegendEntries), 'legendEntries');
    writeRef(state.results, prepared.projection.nextResults, 'results');
    if (!projectedResultsMatch(readRef(state.results), prepared.projection.nextResults)) {
      throw new Error('Manual legend Result publication did not preserve the staged projection.');
    }
    prepared.publishedResults = readRef(state.results);
    writeRef(
      state.validatedStyleFingerprintByResultKey,
      Object.freeze(cloneJson(to.resultFingerprintByKey)),
      'validatedStyleFingerprintByResultKey'
    );
    const reconciled = typeof reconcile === 'function'
      ? reconcile({ direction, intent, prepared, to })
      : undefined;
    if (isPromiseLike(reconciled)) throw new Error('Manual legend reconciliation must be synchronous.');
    if (reconciled === false) return false;
    const refreshed = typeof refreshPresentation === 'function'
      ? refreshPresentation({ direction, intent, prepared, to })
      : undefined;
    if (isPromiseLike(refreshed)) throw new Error('Manual legend presentation refresh must be synchronous.');
    if (refreshed === false) return false;
    writeRef(
      state.semanticStyleRevision,
      Number(readRef(state.semanticStyleRevision) || 0) + 1,
      'semanticStyleRevision'
    );
    writeRef(state.semanticStyleFingerprint, to.fingerprint, 'semanticStyleFingerprint');
    return true;
  };

  const baseCommand = await buildStyleSnapshotCommand({
    label,
    before,
    after,
    metadata: {
      intentKind,
      caption: text(intent?.caption ?? intent?.oldCaption ?? intent?.entry?.caption),
      documentEpoch: expectedEpoch,
      resultGenerationKey: expectedGeneration,
      affectedResultKeys: resultKeys,
      resultExtent: 'all-affected-results',
      owner: isOrderIntent ? ORDER_OWNER : MANUAL_OWNER
    },
    counters: { affectedResults: resultKeys.length },
    isNoop: (left, right) => (
      manualSignature(left.manualEntries) === manualSignature(right.manualEntries)
      && JSON.stringify(left.legendOrder) === JSON.stringify(right.legendOrder)
    ),
    validateCurrent,
    prepareTransition,
    commitTransition,
    captureExactState: () => captureExactState({
      state,
      mounted: currentMounted(getMountedContext),
      captureRuntimeState
    }),
    restoreExactState: (snapshot, context) => restoreExactState({
      state,
      snapshot,
      restoreMountedProjection,
      restoreRuntimeState,
      context,
      intent
    }),
    assertCommitted: ({ to, prepared }) => (
      manualSignature(readRef(state.manualLegendEntries) || []) === manualSignature(to.manualEntries)
      && JSON.stringify(normalizeLegendOrderIntent(readRef(state.legendOrderIntent) || []))
        === JSON.stringify(to.legendOrder)
      && currentStyleFingerprint(state) === to.fingerprint
      && text(readRef(state.semanticStyleFingerprint)) === to.fingerprint
      && ledgerIsCurrent(
        readRef(state.validatedStyleFingerprintByResultKey) || {},
        resultKeys,
        to.fingerprint
      )
      && readRef(state.results) === prepared.publishedResults
      && projectedResultsMatch(readRef(state.results), prepared.projection.nextResults)
    )
  });

  if (baseCommand.noop) return baseCommand;
  const run = async (method, options) => {
    if (Number(readRef(state.semanticStyleRevision) || 0) !== transitionRevision) return false;
    const result = await method(options);
    if (result !== false) {
      transitionRevision = Number(readRef(state.semanticStyleRevision) || 0);
      firstApplyPending = false;
    }
    return result;
  };
  return Object.freeze({
    ...baseCommand,
    apply: (options = {}) => run(baseCommand.apply, options),
    revert: (options = {}) => run(baseCommand.revert, options),
    compensate: (options = {}) => {
      const restored = baseCommand.compensate(options);
      if (!isPromiseLike(restored) && restored !== false) {
        transitionRevision = Number(options?.snapshot?.revision ?? readRef(state.semanticStyleRevision) ?? 0);
      }
      return restored;
    }
  });
};
