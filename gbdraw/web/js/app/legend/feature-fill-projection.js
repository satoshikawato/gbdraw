import {
  normalizeRuleColor,
  resolveOrderedSpecificRule
} from '../feature-editor/fill-view-model.js';
import {
  biologicalFeatureKey,
  catalogResultKey
} from '../../services/feature-catalog.js';
import { parseCompositionMetadata } from '../legend-layout/composition-actions.js';
import { getAllFeatureLegendGroups, parseTransformXY } from './utils.js';
import { cloneJsonData as cloneJson } from '../../services/json-clone.js';

const SVG_NS = 'http://www.w3.org/2000/svg';
const text = (value) => String(value ?? '').trim();
const captionKey = (value) => text(value).toLowerCase();

const deepFreeze = (value) => {
  if (!value || typeof value !== 'object' || Object.isFrozen(value)) return value;
  Object.values(value).forEach(deepFreeze);
  return Object.freeze(value);
};

const resultEntries = (value, resultKey) => {
  if (value instanceof Map) return value.get(resultKey) || [];
  if (Array.isArray(value)) {
    const projection = value.find((entry) => text(entry?.resultKey) === resultKey);
    return projection?.entries || [];
  }
  return value && typeof value === 'object' ? value[resultKey] || [] : [];
};

const renderedFeaturesForItem = (item) => {
  const renderedByFeature = new Map();
  (Array.isArray(item?.features) ? item.features : []).forEach((feature) => {
    const key = biologicalFeatureKey(feature?.recordKey, feature?.biologicalFeatureId);
    if (!key) return;
    if (!renderedByFeature.has(key)) renderedByFeature.set(key, []);
    renderedByFeature.get(key).push(feature);
  });
  const biological = Array.isArray(item?.biologicalFeatures)
    ? item.biologicalFeatures
    : [];
  if (biological.length === 0) {
    return (Array.isArray(item?.features) ? item.features : []).map((feature) => ({
      feature,
      rendered: [feature]
    }));
  }
  return biological.map((feature) => ({
    feature,
    rendered: renderedByFeature.get(
      biologicalFeatureKey(feature?.recordKey, feature?.biologicalFeatureId)
    ) || []
  }));
};

const featureCaption = (feature, rule) => text(
  rule?.cap
  || feature?.type
  || feature?.effectiveLegendCaption
  || feature?.effective_legend_caption
  || feature?.legendCaption
  || feature?.legend_caption
);

const featureColor = (feature, rule, paletteColors, svgDefaultColor) => (
  normalizeRuleColor(rule?.color)
  || normalizeRuleColor(paletteColors?.[feature?.type])
  || normalizeRuleColor(feature?.fill ?? feature?.color)
  || normalizeRuleColor(svgDefaultColor)
  || '#cccccc'
);

const preservedEntry = (existingEntries, caption) => (
  existingEntries.find((entry) => captionKey(entry?.caption) === captionKey(caption)) || null
);

const mergeEntryStyle = (entry, previous) => ({
  ...(previous ? cloneJson(previous) : {}),
  ...entry,
  caption: entry.caption,
  originalCaption: text(previous?.originalCaption) || entry.caption,
  color: entry.color,
  featureIds: [...entry.featureIds]
});

const orderEntries = (entries, existingEntries, legendOrder = []) => {
  const preferredOrder = new Map();
  (Array.isArray(legendOrder) ? legendOrder : []).forEach((caption, index) => {
    const key = captionKey(caption);
    if (key && !preferredOrder.has(key)) preferredOrder.set(key, index);
  });
  const existingOrder = new Map(
    existingEntries.map((entry, index) => [captionKey(entry?.caption), index])
  );
  return [...entries].sort((left, right) => {
    const leftPreferred = preferredOrder.get(captionKey(left.caption));
    const rightPreferred = preferredOrder.get(captionKey(right.caption));
    if (leftPreferred !== undefined || rightPreferred !== undefined) {
      if (leftPreferred === undefined) return 1;
      if (rightPreferred === undefined) return -1;
      if (leftPreferred !== rightPreferred) return leftPreferred - rightPreferred;
    }
    const leftIndex = existingOrder.get(captionKey(left.caption));
    const rightIndex = existingOrder.get(captionKey(right.caption));
    if (leftIndex !== undefined || rightIndex !== undefined) {
      if (leftIndex === undefined) return 1;
      if (rightIndex === undefined) return -1;
      if (leftIndex !== rightIndex) return leftIndex - rightIndex;
    }
    return left.firstUse - right.firstUse || left.caption.localeCompare(right.caption);
  });
};

/**
 * Resolve actually used Feature legend groups independently for every Result.
 * Specific rules and palette colors remain the semantic owners; the returned
 * entries are immutable projection data.
 */
export const deriveUsedFeatureFillGroupsByResult = ({
  catalog,
  rules = [],
  paletteColors = {},
  manualLegendEntries = [],
  legendOrder = [],
  existingEntriesByResult = {},
  svgDefaultColor = '#cccccc'
} = {}) => {
  const items = Array.isArray(catalog?.items) ? catalog.items : [];
  const normalizedRules = Array.isArray(rules) ? rules : [];
  const projections = [];
  const counters = {
    results: items.length,
    renderedFeatures: 0,
    ruleResolutionUpperBound: 0,
    usedFeatureGroups: 0,
    manualGroups: 0
  };

  items.forEach((item, resultIndex) => {
    const resultKey = catalogResultKey(item);
    if (!resultKey) throw new Error(`Legend projection Result ${resultIndex} has no durable key.`);
    const existingEntries = resultEntries(existingEntriesByResult, resultKey);
    const groups = new Map();
    let firstUse = 0;

    renderedFeaturesForItem(item).forEach(({ feature, rendered }) => {
      if (rendered.length === 0 || feature?.hidden === true || feature?.rendered === false) return;
      counters.renderedFeatures += 1;
      counters.ruleResolutionUpperBound += normalizedRules.length;
      const resolved = resolveOrderedSpecificRule(feature, normalizedRules);
      const color = featureColor(
        feature,
        resolved?.rule,
        paletteColors,
        svgDefaultColor
      );
      if (color === 'none') return;
      const caption = featureCaption(feature, resolved?.rule);
      if (!caption) return;
      const key = captionKey(caption);
      const previous = groups.get(key);
      if (previous && previous.color !== color) {
        throw new Error(`Legend caption "${caption}" resolves to more than one color.`);
      }
      const featureIds = rendered.map((entry) => text(entry?.svgId ?? entry?.svg_id)).filter(Boolean);
      if (previous) {
        featureIds.forEach((featureId) => previous.featureIds.add(featureId));
        return;
      }
      groups.set(key, {
        caption,
        color,
        featureIds: new Set(featureIds),
        owner: 'feature',
        firstUse: firstUse++
      });
    });

    (Array.isArray(manualLegendEntries) ? manualLegendEntries : []).forEach((rawEntry) => {
      const caption = text(rawEntry?.caption);
      const color = normalizeRuleColor(rawEntry?.color);
      if (!caption || !color || color === 'none') return;
      const key = captionKey(caption);
      const featureGroup = groups.get(key);
      if (featureGroup) {
        if (featureGroup.color !== color) {
          throw new Error(`Manual legend caption "${caption}" conflicts with a used Feature color.`);
        }
        Object.assign(featureGroup, {
          ...cloneJson(rawEntry),
          caption: featureGroup.caption,
          color: featureGroup.color,
          featureIds: featureGroup.featureIds,
          owner: 'feature+manual'
        });
        return;
      }
      groups.set(key, {
        ...cloneJson(rawEntry),
        caption,
        color,
        featureIds: new Set(),
        owner: 'manual',
        firstUse: firstUse++
      });
      counters.manualGroups += 1;
    });

    const entries = orderEntries([...groups.values()], existingEntries, legendOrder).map((entry) => {
      const normalized = mergeEntryStyle({
        ...entry,
        featureIds: [...entry.featureIds].sort()
      }, preservedEntry(existingEntries, entry.caption));
      delete normalized.firstUse;
      return normalized;
    });
    counters.usedFeatureGroups += entries.filter((entry) => entry.featureIds.length > 0).length;
    projections.push({
      resultKey,
      resultIndex,
      resultName: text(item?.resultName),
      entries
    });
  });

  return deepFreeze({ projections, counters });
};

/**
 * Keep Result-owned non-Feature legend rows while replacing the semantic
 * Feature/manual projection. `beforeProjection` identifies rows owned by the
 * semantic Feature/manual side; every other existing row remains SVG-owned.
 */
export const preserveResultLocalNonFeatureLegendEntries = ({
  beforeProjection = null,
  afterProjection = null,
  existingEntries = []
} = {}) => {
  if (!afterProjection) throw new Error('A target Feature legend projection is required.');
  const beforeCaptions = new Set(
    (beforeProjection?.entries || []).map((entry) => captionKey(entry?.caption)).filter(Boolean)
  );
  const desired = new Map(
    (afterProjection?.entries || []).map((entry) => [captionKey(entry?.caption), cloneJson(entry)])
  );
  (Array.isArray(existingEntries) ? existingEntries : []).forEach((entry) => {
    const key = captionKey(entry?.caption);
    if (!key || beforeCaptions.has(key) || desired.has(key)) return;
    desired.set(key, cloneJson(entry));
  });
  const entries = [];
  (Array.isArray(existingEntries) ? existingEntries : []).forEach((entry) => {
    const key = captionKey(entry?.caption);
    if (!desired.has(key)) return;
    entries.push(desired.get(key));
    desired.delete(key);
  });
  entries.push(...desired.values());
  return deepFreeze({ ...cloneJson(afterProjection), entries });
};

const directEntryGroups = (targetGroup) => {
  const direct = Array.from(targetGroup?.children || []).filter((child) => (
    text(child?.localName ?? child?.tagName).toLowerCase() === 'g'
    && child.hasAttribute?.('data-legend-key')
  ));
  return direct.length > 0
    ? direct
    : Array.from(targetGroup?.querySelectorAll?.('g[data-legend-key]') || []);
};

const entryCaption = (group) => text(group?.getAttribute?.('data-legend-key'));

const entrySwatch = (group) => Array.from(group?.querySelectorAll?.('path') || []).find((path) => {
  const fill = text(path?.getAttribute?.('fill')).toLowerCase();
  return fill && fill !== 'none' && !fill.startsWith('url(');
}) || null;

const removeNode = (node) => {
  if (typeof node?.remove === 'function') node.remove();
  else if (node?.parentElement?.removeChild) node.parentElement.removeChild(node);
  else if (Array.isArray(node?.parentElement?.children)) {
    node.parentElement.children = node.parentElement.children.filter((child) => child !== node);
  }
};

const anchorFor = (group) => {
  const groupPosition = parseTransformXY(group?.getAttribute?.('transform'));
  const child = group?.querySelector?.('text') || entrySwatch(group);
  const childPosition = parseTransformXY(child?.getAttribute?.('transform'));
  return { x: groupPosition.x + childPosition.x, y: groupPosition.y + childPosition.y };
};

const moveToAnchor = (group, anchor) => {
  const current = anchorFor(group);
  const dx = anchor.x - current.x;
  const dy = anchor.y - current.y;
  if (Math.abs(dx) < 1e-9 && Math.abs(dy) < 1e-9) return;
  const groupHasTransform = Boolean(group?.getAttribute?.('transform'));
  const transformed = (groupHasTransform
    ? [group]
    : [entrySwatch(group), group?.querySelector?.('text')])
    .filter((node, index, values) => values.indexOf(node) === index && node?.getAttribute?.('transform'));
  if (transformed.length === 0) {
    group.setAttribute?.('transform', `translate(${dx}, ${dy})`);
    return;
  }
  transformed.forEach((node) => {
    const currentPosition = parseTransformXY(node.getAttribute('transform'));
    node.setAttribute('transform', `translate(${currentPosition.x + dx}, ${currentPosition.y + dy})`);
  });
};

const legendIsHorizontal = (targetGroup) => {
  let current = targetGroup;
  while (current) {
    const id = text(current?.id ?? current?.getAttribute?.('id')).toLowerCase();
    const orientation = text(current?.getAttribute?.('data-gbdraw-orientation')).toLowerCase();
    if (orientation === 'h' || id.includes('horizontal') || id.endsWith('_h')) return true;
    if (orientation === 'v' || id.includes('vertical') || id.endsWith('_v')) return false;
    current = current.parentElement || null;
  }
  return false;
};

const copyComputedFont = (source, target, view) => {
  const style = source && typeof view?.getComputedStyle === 'function'
    ? view.getComputedStyle(source)
    : null;
  for (const property of [
    'font-family', 'font-size', 'font-style', 'font-weight',
    'font-stretch', 'font-variant', 'letter-spacing', 'word-spacing'
  ]) {
    const value = style?.getPropertyValue?.(property) || source?.getAttribute?.(property);
    if (text(value)) target.setAttribute(property, value);
  }
};

/** Browser-owned SVG text measurement with no rendering-runtime dependency. */
export const measureLegendTextInBrowser = async ({ caption, textElement, svg } = {}) => {
  const document = svg?.ownerDocument || textElement?.ownerDocument || globalThis.document;
  if (!document?.createElementNS) throw new Error('Browser SVG text measurement is unavailable.');
  if (document.fonts?.ready) await document.fonts.ready;
  const host = document.createElementNS(SVG_NS, 'svg');
  host.setAttribute('aria-hidden', 'true');
  host.setAttribute('width', '0');
  host.setAttribute('height', '0');
  host.setAttribute('style', 'position:absolute;left:-100000px;top:-100000px;overflow:visible');
  const measuredText = document.createElementNS(SVG_NS, 'text');
  measuredText.textContent = text(caption);
  copyComputedFont(textElement, measuredText, document.defaultView || globalThis.window);
  host.appendChild(measuredText);
  const owner = document.body || document.documentElement;
  if (!owner?.appendChild) throw new Error('Browser SVG text measurement host is unavailable.');
  owner.appendChild(host);
  try {
    const measuredLength = Number(measuredText.getComputedTextLength?.());
    const bounds = measuredText.getBBox?.();
    const width = Number.isFinite(measuredLength) && measuredLength >= 0
      ? measuredLength
      : Number(bounds?.width);
    const height = Number(bounds?.height);
    if (!Number.isFinite(width) || width <= 0) {
      throw new Error(`Browser SVG text measurement failed for "${text(caption)}".`);
    }
    return { width, height: Number.isFinite(height) && height >= 0 ? height : 0 };
  } finally {
    removeNode(host);
  }
};

const validateDesiredEntries = (entries) => {
  const normalized = [];
  const byCaption = new Map();
  (Array.isArray(entries) ? entries : []).forEach((entry) => {
    const caption = text(entry?.caption);
    const color = normalizeRuleColor(entry?.color);
    if (!caption || !color || color === 'none') return;
    const key = captionKey(caption);
    if (byCaption.has(key) && byCaption.get(key).color !== color) {
      throw new Error(`Legend caption "${caption}" has conflicting colors.`);
    }
    if (byCaption.has(key)) return;
    const next = { ...cloneJson(entry), caption, color };
    byCaption.set(key, next);
    normalized.push(next);
  });
  return normalized;
};

const defaultLayout = ({ assigned, entries, slots, metrics, measurements, targetGroup }) => {
  const horizontal = legendIsHorizontal(targetGroup);
  let last = slots.at(-1) || { x: 0, y: 0 };
  entries.forEach((entry, index) => {
    let slot = slots[index];
    if (!slot) {
      const previous = entries[index - 1];
      const previousWidth = measurements.get(captionKey(previous?.caption))?.width || 0;
      slot = horizontal
          ? {
            x: last.x + (2 * metrics.textXOffset) + previousWidth,
            y: last.y
          }
        : { x: last.x, y: last.y + metrics.lineHeight };
    }
    moveToAnchor(assigned.get(captionKey(entry.caption)), slot);
    last = slot;
  });
};

/**
 * Clone and stage one Result-local legend. Missing metadata, templates, fonts,
 * or measurements reject before the source root is mutated.
 */
export const prepareFeatureFillLegendProjection = async ({
  sourceSvg,
  entries,
  allowAbsentLegend = false,
  measureText = measureLegendTextInBrowser,
  parseMetadata = parseCompositionMetadata,
  getTargetGroups = getAllFeatureLegendGroups,
  templateForEntry = null,
  layoutEntries = defaultLayout
} = {}) => {
  if (!sourceSvg?.cloneNode) throw new Error('A detached-capable SVG root is required for legend projection.');
  const desired = validateDesiredEntries(entries);
  const sourceTargetGroups = getTargetGroups(sourceSvg);
  if (!Array.isArray(sourceTargetGroups) || sourceTargetGroups.length === 0) {
    if (!allowAbsentLegend) throw new Error('This diagram has no Feature legend target.');
    const svg = sourceSvg.cloneNode(true);
    if (!svg || svg === sourceSvg) {
      throw new Error('Legend projection requires an independent SVG clone.');
    }
    return Object.freeze({
      svg,
      entries: deepFreeze([]),
      metadata: null,
      counters: deepFreeze({
        targetGroups: 0,
        colorUpdates: 0,
        strokeUpdates: 0,
        renames: 0,
        additions: 0,
        removals: 0,
        measurements: 0
      })
    });
  }
  const metadata = parseMetadata(sourceSvg);
  const metrics = metadata?.legendReflow;
  if (!metrics) throw new Error('This diagram has no legend reflow metadata. Regenerate it before editing the legend.');
  const svg = sourceSvg.cloneNode(true);
  if (!svg || svg === sourceSvg) {
    throw new Error('Legend projection requires an independent SVG clone.');
  }
  const targetGroups = getTargetGroups(svg);
  if (!Array.isArray(targetGroups) || targetGroups.length === 0) {
    throw new Error('This diagram has no Feature legend target.');
  }
  const counters = {
    targetGroups: targetGroups.length,
    colorUpdates: 0,
    strokeUpdates: 0,
    renames: 0,
    additions: 0,
    removals: 0,
    measurements: 0
  };

  for (const targetGroup of targetGroups) {
    const initial = directEntryGroups(targetGroup);
    const slots = initial.map(anchorFor).sort((left, right) => left.y - right.y || left.x - right.x);
    const unused = new Set(initial);
    const assigned = new Map();
    desired.forEach((entry) => {
      const exact = initial.find((group) => (
        unused.has(group) && entryCaption(group) === entry.caption
      ));
      if (!exact) return;
      unused.delete(exact);
      assigned.set(captionKey(entry.caption), exact);
    });
    const missing = desired.filter((entry) => !assigned.has(captionKey(entry.caption)));
    const reusable = [...unused];
    const measurements = new Map();

    for (let index = 0; index < missing.length; index += 1) {
      const entry = missing[index];
      let group = reusable[index] || null;
      if (group) {
        unused.delete(group);
        counters.renames += 1;
      } else {
        const supplied = typeof templateForEntry === 'function'
          ? templateForEntry({ entry, targetGroup, sourceSvg: svg, initial })
          : null;
        const template = supplied || initial[0] || assigned.values().next().value || null;
        group = template?.cloneNode?.(true) || null;
        if (!group || group === template) {
          throw new Error(`No independent legend entry template is available for "${entry.caption}".`);
        }
        targetGroup.appendChild(group);
        counters.additions += 1;
      }
      const label = group.querySelector?.('text');
      if (!label) throw new Error(`Legend entry "${entry.caption}" has no text template.`);
      const measurement = await measureText({ caption: entry.caption, textElement: label, svg });
      const width = Number(measurement?.width ?? measurement);
      if (!Number.isFinite(width) || width <= 0) {
        throw new Error(`Legend text measurement failed for "${entry.caption}".`);
      }
      measurements.set(captionKey(entry.caption), {
        width,
        height: Number(measurement?.height) || 0
      });
      counters.measurements += 1;
      group.setAttribute('data-legend-key', entry.caption);
      label.textContent = entry.caption;
      assigned.set(captionKey(entry.caption), group);
    }

    if (legendIsHorizontal(targetGroup) && desired.length > slots.length && slots.length > 0) {
      const previous = desired[slots.length - 1];
      const previousKey = captionKey(previous?.caption);
      if (!measurements.has(previousKey)) {
        const previousLabel = assigned.get(previousKey)?.querySelector?.('text');
        if (!previousLabel) {
          throw new Error(`Legend entry "${previous?.caption}" has no text template.`);
        }
        const measurement = await measureText({
          caption: previous.caption,
          textElement: previousLabel,
          svg
        });
        const width = Number(measurement?.width ?? measurement);
        if (!Number.isFinite(width) || width <= 0) {
          throw new Error(`Legend text measurement failed for "${previous.caption}".`);
        }
        measurements.set(previousKey, {
          width,
          height: Number(measurement?.height) || 0
        });
        counters.measurements += 1;
      }
    }

    [...unused].forEach((group) => {
      removeNode(group);
      counters.removals += 1;
    });
    desired.forEach((entry) => {
      const group = assigned.get(captionKey(entry.caption));
      if (!group) throw new Error(`Legend entry "${entry.caption}" could not be staged.`);
      const swatch = entrySwatch(group);
      if (!swatch) throw new Error(`Legend entry "${entry.caption}" has no color swatch template.`);
      if (text(swatch.getAttribute('fill')).toLowerCase() !== entry.color) {
        swatch.setAttribute('fill', entry.color);
        counters.colorUpdates += 1;
      }
      const strokeColor = text(entry.strokeColor);
      const currentStroke = swatch.getAttribute('stroke');
      if (strokeColor ? currentStroke !== strokeColor : currentStroke !== null) {
        if (strokeColor) swatch.setAttribute('stroke', strokeColor);
        else swatch.removeAttribute('stroke');
        counters.strokeUpdates += 1;
      }
      const strokeWidth = entry.strokeWidth;
      const currentStrokeWidth = swatch.getAttribute('stroke-width');
      if (strokeWidth === null || strokeWidth === undefined || strokeWidth === '') {
        if (currentStrokeWidth !== null) {
          swatch.removeAttribute('stroke-width');
          counters.strokeUpdates += 1;
        }
      } else if (currentStrokeWidth !== String(Number(strokeWidth))) {
        swatch.setAttribute('stroke-width', String(Number(strokeWidth)));
        counters.strokeUpdates += 1;
      }
      targetGroup.appendChild(group);
    });
    const result = layoutEntries({
      assigned,
      entries: desired,
      slots,
      metrics,
      measurements,
      targetGroup,
      svg
    });
    if (result && typeof result.then === 'function') {
      throw new Error('Legend layout commit must be synchronous after measurement.');
    }
  }

  return Object.freeze({
    svg,
    entries: deepFreeze(cloneJson(desired)),
    metadata: deepFreeze({ legendReflow: cloneJson(metrics) }),
    counters: deepFreeze(counters)
  });
};

/** Prepare every Result-local legend before publishing any of them. */
export const prepareFeatureFillLegendProjections = async ({
  sourcesByResult,
  projections,
  ...options
} = {}) => {
  const staged = new Map();
  const counters = {
    results: 0,
    measurements: 0,
    additions: 0,
    renames: 0,
    removals: 0,
    colorUpdates: 0,
    strokeUpdates: 0
  };
  for (const projection of Array.isArray(projections) ? projections : []) {
    const resultKey = text(projection?.resultKey);
    const sourceSvg = sourcesByResult instanceof Map
      ? sourcesByResult.get(resultKey)
      : sourcesByResult?.[resultKey];
    if (!sourceSvg) throw new Error(`Legend projection source is unavailable for Result ${resultKey}.`);
    const prepared = await prepareFeatureFillLegendProjection({
      ...options,
      sourceSvg,
      entries: projection.entries
    });
    staged.set(resultKey, prepared);
    counters.results += 1;
    for (const key of [
      'measurements', 'additions', 'renames', 'removals', 'colorUpdates', 'strokeUpdates'
    ]) {
      counters[key] += prepared.counters[key];
    }
  }
  return Object.freeze({ staged, counters: deepFreeze(counters) });
};
