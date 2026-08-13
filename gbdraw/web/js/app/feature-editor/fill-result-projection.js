import {
  admitProjectedSvgResult,
  getCommittedSvgContent,
  getCommittedSvgResultMetadata,
  parseCommittedSvgResultRoot,
  projectMountedSvgResult
} from '../../services/svg-result-ingestion.js';
import {
  filterFeatureFillTargets,
  getFeatureElementIndex
} from '../feature-dom.js';
import { catalogResultKey } from '../../services/feature-catalog.js';
import { normalizeRuleColor } from './fill-view-model.js';

const text = (value) => String(value ?? '').trim();

const itemResultKey = (item) => (
  catalogResultKey(item) || ''
);

const fillValueForFeature = (feature, fillByFeatureKey) => {
  const pairKey = JSON.stringify([
    text(feature?.recordKey),
    text(feature?.biologicalFeatureId)
  ]);
  const value = fillByFeatureKey.get(pairKey);
  if (value === null || value === undefined) return null;
  return normalizeRuleColor(value && typeof value === 'object' ? value.color : value);
};

const applyItemFills = (svg, item, fillByFeatureKey, requiredTargetKeys = null) => {
  const index = getFeatureElementIndex(svg);
  let changed = 0;
  const coveredPairs = new Set();
  (Array.isArray(item?.features) ? item.features : []).forEach((feature) => {
    const pairKey = JSON.stringify([
      text(feature?.recordKey),
      text(feature?.biologicalFeatureId)
    ]);
    if (requiredTargetKeys && !requiredTargetKeys.has(pairKey)) return;
    const fill = fillValueForFeature(feature, fillByFeatureKey);
    if (!fill) {
      if (requiredTargetKeys) {
        throw new Error(`Feature fill value is unavailable for rendered target ${pairKey}.`);
      }
      return;
    }
    const svgId = text(feature?.svgId ?? feature?.svg_id);
    const elements = index.get(svgId) || [];
    const targets = filterFeatureFillTargets(elements);
    if (targets.length === 0) {
      throw new Error(`Feature fill SVG target is unavailable: ${svgId}`);
    }
    coveredPairs.add(pairKey);
    targets.forEach((element) => {
      if (element.getAttribute('fill') === fill) return;
      element.setAttribute('fill', fill);
      changed += 1;
    });
  });
  if (requiredTargetKeys && coveredPairs.size !== requiredTargetKeys.size) {
    throw new Error(
      `Feature fill target coverage is incomplete (${coveredPairs.size}/${requiredTargetKeys.size}).`
    );
  }
  return { changed, covered: coveredPairs.size };
};

const resultKeyMap = (catalog) => {
  const bindings = new Map();
  (Array.isArray(catalog?.items) ? catalog.items : []).forEach((item, resultIndex) => {
    const resultKey = itemResultKey(item);
    if (!resultKey || bindings.has(resultKey)) {
      throw new Error(`Feature fill Result identity is missing or duplicate at index ${resultIndex}.`);
    }
    bindings.set(resultKey, { item, resultIndex });
  });
  return bindings;
};

/**
 * Prepare all affected Result replacements before publishing any state.
 * The mounted Result is cloned in memory; inactive Results use detached admitted roots.
 */
export const prepareFeatureFillResultProjection = ({
  results,
  catalog,
  fillsByTargetKey = {},
  affectedResultKeys = [],
  mountedResultIndex = null,
  mountedSvg = null,
  mounted = null,
  preparedSvgByResultKey = null,
  targetFeatureKeysByResult = null,
  parser = globalThis.DOMParser || globalThis.window?.DOMParser
} = {}) => {
  if (!Array.isArray(results)) throw new Error('Feature fill projection requires Results.');
  const byResultKey = resultKeyMap(catalog);
  const mountedIndexValue = mountedResultIndex ?? mounted?.resultIndex;
  const activeResultIndex = Number.isInteger(Number(mountedIndexValue))
    && Number(mountedIndexValue) >= 0
    ? Number(mountedIndexValue)
    : null;
  const activeSvg = mountedSvg || mounted?.svg || null;
  if (activeResultIndex !== null && activeResultIndex >= results.length) {
    throw new Error(`Mounted Feature fill Result index is unavailable: ${activeResultIndex}`);
  }
  if (
    text(mounted?.resultKey)
    && itemResultKey(catalog?.items?.[activeResultIndex]) !== text(mounted.resultKey)
  ) {
    throw new Error('Mounted Feature fill Result identity is stale.');
  }
  const normalizedTargetKeys = (Array.isArray(affectedResultKeys) ? affectedResultKeys : [])
    .map(text)
    .filter(Boolean);
  const targets = new Set(normalizedTargetKeys);
  if (targets.size === 0 || targets.size !== normalizedTargetKeys.length) {
    throw new Error('Feature fill projection requires unique affected Result identities.');
  }
  const rawFillByFeatureKey = fillsByTargetKey instanceof Map
    ? new Map(fillsByTargetKey)
    : new Map(Object.entries(fillsByTargetKey || {}));
  const fillByResultKey = new Map();
  (Array.isArray(catalog?.items) ? catalog.items : []).forEach((item, resultIndex) => {
    const resultKey = itemResultKey(item);
    if (!targets.has(resultKey)) return;
    const byPair = new Map();
    (Array.isArray(item?.biologicalFeatures) ? item.biologicalFeatures : []).forEach((feature) => {
      const targetKey = JSON.stringify([
        resultKey,
        text(feature?.recordKey),
        text(feature?.biologicalFeatureId)
      ]);
      const pairKey = JSON.stringify([
        text(feature?.recordKey),
        text(feature?.biologicalFeatureId)
      ]);
      const value = rawFillByFeatureKey.has(targetKey)
        ? rawFillByFeatureKey.get(targetKey)
        : rawFillByFeatureKey.get(pairKey);
      if (value !== undefined && value !== null) byPair.set(pairKey, value);
    });
    fillByResultKey.set(resultKey, byPair);
  });
  const candidates = new Map();
  const stagedRoots = new Map();
  const admissionMetadataByResultKey = {};
  const counters = {
    affectedResults: targets.size,
    mountedSerializations: 0,
    detachedPasses: 0,
    changedResults: 0
  };

  for (const resultKey of targets) {
    const binding = byResultKey.get(resultKey);
    if (!binding) throw new Error(`Feature fill target Result is unavailable: ${resultKey}`);
    const { item, resultIndex } = binding;
    const fillByFeatureKey = fillByResultKey.get(resultKey) || new Map();
    const result = results[resultIndex];
    if (!result) throw new Error(`Feature fill target Result index is unavailable: ${resultIndex}`);
    const preparedRoot = preparedSvgByResultKey instanceof Map
      ? preparedSvgByResultKey.get(resultKey)
      : preparedSvgByResultKey?.[resultKey];
    const plannedFeatureKeys = targetFeatureKeysByResult instanceof Map
      ? targetFeatureKeysByResult.get(resultKey)
      : targetFeatureKeysByResult?.[resultKey];
    const requiredPairKeys = Array.isArray(plannedFeatureKeys) ? new Set() : null;
    if (requiredPairKeys) {
      const biologicalPairKeys = new Set(
        (Array.isArray(item?.biologicalFeatures) ? item.biologicalFeatures : []).map((feature) => (
          JSON.stringify([
            text(feature?.recordKey),
            text(feature?.biologicalFeatureId)
          ])
        ))
      );
      const renderedPairKeys = new Set(
        (Array.isArray(item?.features) ? item.features : []).map((feature) => (
          JSON.stringify([
            text(feature?.recordKey),
            text(feature?.biologicalFeatureId)
          ])
        ))
      );
      const plannedPairKeys = new Set();
      plannedFeatureKeys.forEach((featureKey) => {
        let values;
        try {
          values = JSON.parse(featureKey);
        } catch {
          throw new Error(`Feature fill target identity is malformed: ${text(featureKey)}`);
        }
        if (
          !Array.isArray(values)
          || values.length !== 3
          || text(values[0]) !== resultKey
          || !text(values[1])
          || !text(values[2])
        ) {
          throw new Error(`Feature fill target identity is cross-Result or incomplete: ${text(featureKey)}`);
        }
        const pairKey = JSON.stringify([text(values[1]), text(values[2])]);
        if (!biologicalPairKeys.has(pairKey)) {
          throw new Error(`Feature fill target identity is unavailable: ${text(featureKey)}`);
        }
        if (plannedPairKeys.has(pairKey)) {
          throw new Error(`Feature fill target identity is duplicated: ${text(featureKey)}`);
        }
        plannedPairKeys.add(pairKey);
        if (renderedPairKeys.has(pairKey)) requiredPairKeys.add(pairKey);
      });
      if (plannedPairKeys.size === 0) {
        throw new Error(`Feature fill target coverage is empty for Result ${resultKey}.`);
      }
    }
    let candidate;
    if (preparedRoot) {
      if (!preparedRoot?.cloneNode) {
        throw new Error(`Prepared Feature fill SVG is invalid for Result ${resultKey}.`);
      }
      const svg = preparedRoot;
      applyItemFills(svg, item, fillByFeatureKey, requiredPairKeys);
      candidate = admitProjectedSvgResult(result, svg);
      stagedRoots.set(resultKey, svg);
      if (resultIndex === activeResultIndex) counters.mountedSerializations += 1;
      else counters.detachedPasses += 1;
    } else if (resultIndex === activeResultIndex) {
      if (!activeSvg) throw new Error('The affected mounted Result has no mounted SVG root.');
      candidate = projectMountedSvgResult(result, activeSvg, {
        resultIndex,
        transformSvg: (svg) => {
          applyItemFills(svg, item, fillByFeatureKey, requiredPairKeys);
          stagedRoots.set(resultKey, svg);
        }
      });
      counters.mountedSerializations += 1;
    } else {
      const root = parseCommittedSvgResultRoot(result, parser);
      applyItemFills(root, item, fillByFeatureKey, requiredPairKeys);
      stagedRoots.set(resultKey, root);
      candidate = admitProjectedSvgResult(result, root);
      counters.detachedPasses += 1;
    }
    const selectedResult = getCommittedSvgContent(result) === candidate.content
      ? result
      : candidate;
    admissionMetadataByResultKey[resultKey] = Object.freeze({
      before: getCommittedSvgResultMetadata(result) || null,
      after: getCommittedSvgResultMetadata(selectedResult) || null
    });
    if (selectedResult === result) {
      candidates.set(resultIndex, result);
    } else {
      candidates.set(resultIndex, selectedResult);
      counters.changedResults += 1;
    }
  }

  const nextResults = results.map((result, resultIndex) => (
    candidates.has(resultIndex) ? candidates.get(resultIndex) : result
  ));
  return Object.freeze({
    nextResults,
    previousResults: results,
    affectedResultKeys: Object.freeze([...targets]),
    admissionMetadataByResultKey: Object.freeze(admissionMetadataByResultKey),
    preparedSvgByResultKey: stagedRoots,
    mountedResultIndex: activeResultIndex,
    mountedResultKey: activeResultIndex === null
      ? ''
      : itemResultKey(catalog?.items?.[activeResultIndex]),
    preparedMountedSvg: activeResultIndex === null
      ? null
      : stagedRoots.get(itemResultKey(catalog?.items?.[activeResultIndex])) || null,
    counters: Object.freeze(counters)
  });
};

export const commitFeatureFillResultProjection = (state, projection) => {
  if (!projection || state.results.value !== projection.previousResults) {
    throw new Error('Feature fill Result projection is stale.');
  }
  state.results.value = projection.nextResults;
  return true;
};
