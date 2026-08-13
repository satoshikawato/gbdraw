import {
  admitProjectedSvgResult,
  getCommittedSvgContent,
  getCommittedSvgResultMetadata,
  parseCommittedSvgResultRoot,
  projectMountedSvgResult
} from '../../services/svg-result-ingestion.js';
import { catalogResultKey } from '../../services/feature-catalog.js';
import { getFeatureElementIndex } from '../feature-dom.js';
import { getAllFeatureLegendGroups } from '../legend/utils.js';
import {
  normalizeFeatureStrokeColor,
  normalizeFeatureStrokeWidth
} from './stroke-view-model.js';

const text = (value) => String(value ?? '').trim();
const hasOwn = (value, key) => Object.prototype.hasOwnProperty.call(value || {}, key);

const normalizedStroke = (value) => {
  if (!value || typeof value !== 'object') {
    throw new Error('Feature stroke projection is missing an effective stroke.');
  }
  const color = hasOwn(value, 'strokeColor')
    ? normalizeFeatureStrokeColor(value.strokeColor)
    : null;
  const width = hasOwn(value, 'strokeWidth')
    ? normalizeFeatureStrokeWidth(value.strokeWidth)
    : null;
  if ((hasOwn(value, 'strokeColor') && !color) || (hasOwn(value, 'strokeWidth') && width === null)) {
    throw new Error('Feature stroke projection contains an invalid effective stroke.');
  }
  return { strokeColor: color, strokeWidth: width };
};

const applyStroke = (element, stroke) => {
  let changed = false;
  if (stroke.strokeColor === null) {
    if (element.hasAttribute('stroke')) {
      element.removeAttribute('stroke');
      changed = true;
    }
  } else if (element.getAttribute('stroke') !== stroke.strokeColor) {
    element.setAttribute('stroke', stroke.strokeColor);
    changed = true;
  }
  if (stroke.strokeWidth === null) {
    if (element.hasAttribute('stroke-width')) {
      element.removeAttribute('stroke-width');
      changed = true;
    }
  } else if (element.getAttribute('stroke-width') !== String(stroke.strokeWidth)) {
    element.setAttribute('stroke-width', stroke.strokeWidth);
    changed = true;
  }
  return changed;
};

const targetPair = (featureKey, expectedResultKey) => {
  let values;
  try {
    values = JSON.parse(featureKey);
  } catch {
    throw new Error(`Feature stroke target identity is malformed: ${text(featureKey)}`);
  }
  if (
    !Array.isArray(values)
    || values.length !== 3
    || text(values[0]) !== expectedResultKey
    || !text(values[1])
    || !text(values[2])
  ) {
    throw new Error(`Feature stroke target identity is cross-Result or incomplete: ${text(featureKey)}`);
  }
  return JSON.stringify([text(values[1]), text(values[2])]);
};

const applyItemStrokes = ({
  svg,
  item,
  resultKey,
  strokesByTargetKey,
  targetFeatureKeys
}) => {
  const index = getFeatureElementIndex(svg);
  const renderedByPair = new Map();
  (Array.isArray(item?.features) ? item.features : []).forEach((feature) => {
    const pairKey = JSON.stringify([
      text(feature?.recordKey),
      text(feature?.biologicalFeatureId)
    ]);
    if (!renderedByPair.has(pairKey)) renderedByPair.set(pairKey, []);
    renderedByPair.get(pairKey).push(feature);
  });
  const biologicalPairs = new Set(
    (Array.isArray(item?.biologicalFeatures) ? item.biologicalFeatures : []).map((feature) => (
      JSON.stringify([
        text(feature?.recordKey),
        text(feature?.biologicalFeatureId)
      ])
    ))
  );
  let changed = 0;
  const covered = new Set();
  targetFeatureKeys.forEach((featureKey) => {
    const pairKey = targetPair(featureKey, resultKey);
    if (!biologicalPairs.has(pairKey)) {
      throw new Error(`Feature stroke target identity is unavailable: ${text(featureKey)}`);
    }
    const rendered = renderedByPair.get(pairKey) || [];
    if (rendered.length === 0) {
      throw new Error(`Feature stroke target has no rendered binding: ${text(featureKey)}`);
    }
    const value = strokesByTargetKey instanceof Map
      ? strokesByTargetKey.get(featureKey)
      : strokesByTargetKey?.[featureKey];
    const stroke = normalizedStroke(value);
    let targetElements = 0;
    rendered.forEach((binding) => {
      const svgId = text(binding?.svgId ?? binding?.svg_id);
      const elements = index.get(svgId) || [];
      if (elements.length === 0) {
        throw new Error(`Feature stroke SVG target is unavailable: ${svgId}`);
      }
      targetElements += elements.length;
      elements.forEach((element) => {
        if (applyStroke(element, stroke)) changed += 1;
      });
    });
    if (targetElements > 0) covered.add(featureKey);
  });
  if (covered.size !== targetFeatureKeys.length) {
    throw new Error(
      `Feature stroke target coverage is incomplete (${covered.size}/${targetFeatureKeys.length}).`
    );
  }
  return changed;
};

const legendSwatch = (entry) => Array.from(entry?.querySelectorAll?.('path') || [])
  .find((path) => {
    const fill = text(path.getAttribute?.('fill'));
    return fill && fill !== 'none' && !fill.startsWith('url(');
  }) || null;

const applyLegendStroke = (svg, caption, value) => {
  if (!caption) return 0;
  const stroke = normalizedStroke(value);
  const groups = getAllFeatureLegendGroups(svg);
  let covered = 0;
  let changed = 0;
  groups.forEach((group) => {
    const entry = Array.from(group.querySelectorAll('g[data-legend-key]'))
      .find((candidate) => text(candidate.getAttribute('data-legend-key')) === caption);
    const swatch = legendSwatch(entry);
    if (!swatch) return;
    covered += 1;
    if (applyStroke(swatch, stroke)) changed += 1;
  });
  if (covered === 0) {
    throw new Error(`Feature stroke legend target is unavailable: ${caption}`);
  }
  return changed;
};

const resultBindings = (catalog) => {
  const byKey = new Map();
  (Array.isArray(catalog?.items) ? catalog.items : []).forEach((item, resultIndex) => {
    const resultKey = catalogResultKey(item);
    if (!resultKey || byKey.has(resultKey)) {
      throw new Error(`Feature stroke Result identity is missing or duplicate at index ${resultIndex}.`);
    }
    byKey.set(resultKey, { item, resultIndex });
  });
  return byKey;
};

/** Prepare every affected Result before the caller publishes one array swap. */
export const prepareFeatureStrokeResultProjection = ({
  results,
  catalog,
  strokesByTargetKey = {},
  affectedResultKeys = [],
  targetFeatureKeysByResult = null,
  legendCaption = '',
  legendStroke = null,
  mountedResultIndex = null,
  mountedSvg = null,
  mounted = null,
  preparedSvgByResultKey = null,
  parser = globalThis.DOMParser || globalThis.window?.DOMParser
} = {}) => {
  if (!Array.isArray(results)) throw new Error('Feature stroke projection requires Results.');
  const byResultKey = resultBindings(catalog);
  const mountedIndexValue = mountedResultIndex ?? mounted?.resultIndex;
  const activeResultIndex = Number.isInteger(Number(mountedIndexValue))
    && Number(mountedIndexValue) >= 0
    ? Number(mountedIndexValue)
    : null;
  const activeSvg = mountedSvg || mounted?.svg || null;
  if (activeResultIndex !== null && activeResultIndex >= results.length) {
    throw new Error(`Mounted Feature stroke Result index is unavailable: ${activeResultIndex}`);
  }
  if (
    text(mounted?.resultKey)
    && catalogResultKey(catalog?.items?.[activeResultIndex]) !== text(mounted.resultKey)
  ) {
    throw new Error('Mounted Feature stroke Result identity is stale.');
  }
  const orderedAffected = (Array.isArray(affectedResultKeys) ? affectedResultKeys : [])
    .map(text).filter(Boolean);
  const affected = new Set(orderedAffected);
  if (affected.size === 0 || affected.size !== orderedAffected.length) {
    throw new Error('Feature stroke projection requires unique affected Result identities.');
  }

  const candidates = new Map();
  const stagedRoots = new Map();
  const admissionMetadataByResultKey = {};
  const counters = {
    affectedResults: affected.size,
    mountedSerializations: 0,
    detachedPasses: 0,
    changedResults: 0,
    featureElementsChanged: 0,
    legendSwatchesChanged: 0
  };

  affected.forEach((resultKey) => {
    const binding = byResultKey.get(resultKey);
    if (!binding) throw new Error(`Feature stroke target Result is unavailable: ${resultKey}`);
    const { item, resultIndex } = binding;
    const result = results[resultIndex];
    if (!result) throw new Error(`Feature stroke target Result index is unavailable: ${resultIndex}`);
    const featureKeys = targetFeatureKeysByResult instanceof Map
      ? targetFeatureKeysByResult.get(resultKey)
      : targetFeatureKeysByResult?.[resultKey];
    if (!Array.isArray(featureKeys) || featureKeys.length === 0) {
      throw new Error(`Feature stroke target coverage is empty for Result ${resultKey}.`);
    }
    const transform = (svg) => {
      counters.featureElementsChanged += applyItemStrokes({
        svg,
        item,
        resultKey,
        strokesByTargetKey,
        targetFeatureKeys: featureKeys
      });
      if (legendCaption) {
        counters.legendSwatchesChanged += applyLegendStroke(svg, legendCaption, legendStroke);
      }
      stagedRoots.set(resultKey, svg);
    };

    const preparedRoot = preparedSvgByResultKey instanceof Map
      ? preparedSvgByResultKey.get(resultKey)
      : preparedSvgByResultKey?.[resultKey];
    let candidate;
    if (preparedRoot) {
      if (!preparedRoot?.cloneNode) {
        throw new Error(`Prepared Feature stroke SVG is invalid for Result ${resultKey}.`);
      }
      transform(preparedRoot);
      candidate = admitProjectedSvgResult(result, preparedRoot);
      if (resultIndex === activeResultIndex) counters.mountedSerializations += 1;
      else counters.detachedPasses += 1;
    } else if (resultIndex === activeResultIndex) {
      if (!activeSvg) throw new Error('The affected mounted Result has no mounted SVG root.');
      candidate = projectMountedSvgResult(result, activeSvg, { resultIndex, transformSvg: transform });
      counters.mountedSerializations += 1;
    } else {
      const root = parseCommittedSvgResultRoot(result, parser);
      transform(root);
      candidate = admitProjectedSvgResult(result, root);
      counters.detachedPasses += 1;
    }
    const selected = getCommittedSvgContent(result) === candidate.content ? result : candidate;
    admissionMetadataByResultKey[resultKey] = Object.freeze({
      before: getCommittedSvgResultMetadata(result) || null,
      after: getCommittedSvgResultMetadata(selected) || null
    });
    candidates.set(resultIndex, selected);
    if (selected !== result) counters.changedResults += 1;
  });

  const nextResults = results.map((result, resultIndex) => (
    candidates.has(resultIndex) ? candidates.get(resultIndex) : result
  ));
  const mountedResultKey = activeResultIndex === null
    ? ''
    : catalogResultKey(catalog?.items?.[activeResultIndex]);
  return Object.freeze({
    nextResults,
    previousResults: results,
    affectedResultKeys: Object.freeze([...affected]),
    admissionMetadataByResultKey: Object.freeze(admissionMetadataByResultKey),
    preparedSvgByResultKey: stagedRoots,
    mountedResultIndex: activeResultIndex,
    mountedResultKey,
    preparedMountedSvg: mountedResultKey ? stagedRoots.get(mountedResultKey) || null : null,
    counters: Object.freeze(counters)
  });
};

export const commitFeatureStrokeResultProjection = (state, projection) => {
  if (!projection || state.results.value !== projection.previousResults) {
    throw new Error('Feature stroke Result projection is stale.');
  }
  state.results.value = projection.nextResults;
  return true;
};
