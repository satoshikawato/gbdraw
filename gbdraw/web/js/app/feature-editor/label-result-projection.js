import {
  admitProjectedSvgResult,
  getCommittedSvgContent,
  getCommittedSvgResultMetadata,
  parseCommittedSvgResultRoot,
  projectMountedSvgResult
} from '../../services/svg-result-ingestion.js';
import { catalogResultKey } from '../../services/feature-catalog.js';

const text = (value) => String(value ?? '').trim();

const EXCLUDED_LABEL_ANCESTOR = [
  '#legend',
  '#feature_legend',
  '#comparison_legend',
  '#horizontal_legend',
  '#vertical_legend',
  '#length_bar',
  'g[data-gbdraw-slot-renderer="ticks"]',
  'g[id="tick"]',
  'g[id^="tick_"]',
  'g[data-gbdraw-role="record-definition"]'
].join(', ');

const getLabelText = (element) => (
  element?.querySelector?.('textPath')?.textContent ?? element?.textContent ?? ''
);

const setLabelText = (element, value) => {
  const target = element?.querySelector?.('textPath') || element;
  if (!target) return false;
  const next = String(value ?? '');
  if (String(target.textContent ?? '') === next) return false;
  target.textContent = next;
  return true;
};

const labelElements = (svg) => {
  const labels = new Set();
  for (const selector of [
    'text[data-label-editable="true"]',
    'g[id="labels"] text, g[id^="labels_"] text',
    'g[id="label_text"] text, g[id^="label_text_"] text',
    'text[dominant-baseline="central"]'
  ]) {
    Array.from(svg?.querySelectorAll?.(selector) || []).forEach((element) => labels.add(element));
  }
  Array.from(svg?.querySelectorAll?.('text > textPath') || []).forEach((textPath) => {
    if (textPath?.parentElement) labels.add(textPath.parentElement);
  });
  return [...labels].filter((element) => !element.closest?.(EXCLUDED_LABEL_ANCESTOR));
};

const featureIdOf = (element) => text(
  element?.getAttribute?.('data-label-feature-id')
  || element?.getAttribute?.('data-gbdraw-label-feature-id')
  || element?.getAttribute?.('data-gbdraw-feature-id')
);

const sourceTextOf = (element) => {
  const explicit = element?.getAttribute?.('data-label-source-text');
  return explicit === null || explicit === undefined ? '' : String(explicit);
};

const normalizeTargetLabel = (value) => {
  if (!value || typeof value !== 'object') {
    throw new Error('Feature label projection is missing an effective label.');
  }
  if (typeof value.text !== 'string') {
    throw new Error('Feature label projection contains invalid label text.');
  }
  return {
    text: value.text,
    sourceText: String(value.sourceText ?? ''),
    fromText: String(value.fromText ?? ''),
    renderedSvgIds: [...new Set(
      (Array.isArray(value.renderedSvgIds) ? value.renderedSvgIds : [])
        .map(text).filter(Boolean)
    )]
  };
};

const applySourceGroup = (elements, sourceText, nextText) => {
  let matched = 0;
  let changed = 0;
  elements.forEach((element) => {
    const source = sourceTextOf(element);
    const current = String(getLabelText(element));
    if (source !== sourceText && (source || current !== sourceText)) return;
    matched += 1;
    if (!source && element?.setAttribute) {
      element.setAttribute('data-label-source-text', sourceText);
    }
    if (setLabelText(element, nextText)) changed += 1;
  });
  return { matched, changed };
};

const applyRenderedGroup = (elements, fromText, nextText) => {
  let matched = 0;
  let changed = 0;
  elements.forEach((element) => {
    if (String(getLabelText(element)) !== fromText) return;
    matched += 1;
    if (setLabelText(element, nextText)) changed += 1;
  });
  return { matched, changed };
};

const exactCandidate = (elements, value, used) => {
  const direct = elements.filter((element) => (
    !used.has(element) && value.renderedSvgIds.includes(featureIdOf(element))
  ));
  if (direct.length === 1) return direct[0];
  if (direct.length > 1) {
    throw new Error(`Feature label SVG identity is ambiguous: ${value.renderedSvgIds.join(', ')}`);
  }
  const byText = elements.filter((element) => {
    if (used.has(element)) return false;
    const current = String(getLabelText(element));
    const source = sourceTextOf(element);
    return (
      (value.fromText && current === value.fromText)
      || (value.sourceText && (source === value.sourceText || (!source && current === value.sourceText)))
      || current === value.text
    );
  });
  return byText.length === 1 ? byText[0] : null;
};

/** Apply one resolved label scope to a staged SVG root. */
export const applyFeatureLabelProjectionToSvg = ({
  svg,
  semanticScope,
  labelsByTargetKey = {},
  targetFeatureKeys = [],
  matchText = '',
  sourceText = ''
} = {}) => {
  if (!svg?.querySelectorAll) throw new Error('Feature label SVG projection is unavailable.');
  const elements = labelElements(svg);
  const keys = [...new Set((Array.isArray(targetFeatureKeys) ? targetFeatureKeys : [])
    .map(text).filter(Boolean))];
  if (keys.length === 0) throw new Error('Feature label projection has no stable targets.');
  const first = normalizeTargetLabel(
    labelsByTargetKey instanceof Map
      ? labelsByTargetKey.get(keys[0])
      : labelsByTargetKey?.[keys[0]]
  );
  if (semanticScope === 'source-annotation-label') {
    const outcome = applySourceGroup(elements, String(sourceText), first.text);
    return Object.freeze({
      targetFeatures: keys.length,
      coveredTargets: keys.length,
      matchedLabels: outcome.matched,
      changedLabels: outcome.changed
    });
  }
  if (semanticScope === 'rendered-label') {
    const outcome = applyRenderedGroup(elements, String(matchText), first.text);
    return Object.freeze({
      targetFeatures: keys.length,
      coveredTargets: keys.length,
      matchedLabels: outcome.matched,
      changedLabels: outcome.changed
    });
  }

  const used = new Set();
  let changed = 0;
  keys.forEach((key) => {
    const value = normalizeTargetLabel(
      labelsByTargetKey instanceof Map ? labelsByTargetKey.get(key) : labelsByTargetKey?.[key]
    );
    const element = exactCandidate(elements, value, used);
    if (!element) throw new Error(`Feature label SVG target is unavailable: ${key}`);
    used.add(element);
    if (!sourceTextOf(element) && value.sourceText && element?.setAttribute) {
      element.setAttribute('data-label-source-text', value.sourceText);
    }
    if (!featureIdOf(element) && value.renderedSvgIds[0] && element?.setAttribute) {
      element.setAttribute('data-label-feature-id', value.renderedSvgIds[0]);
    }
    if (setLabelText(element, value.text)) changed += 1;
  });
  return Object.freeze({
    targetFeatures: keys.length,
    coveredTargets: used.size,
    matchedLabels: used.size,
    changedLabels: changed
  });
};

const resultBindings = (catalog) => {
  const byKey = new Map();
  (Array.isArray(catalog?.items) ? catalog.items : []).forEach((item, resultIndex) => {
    const resultKey = catalogResultKey(item);
    if (!resultKey || byKey.has(resultKey)) {
      throw new Error(`Feature label Result identity is missing or duplicate at index ${resultIndex}.`);
    }
    byKey.set(resultKey, { item, resultIndex });
  });
  return byKey;
};

/** Prepare all affected Result strings before publishing one array swap. */
export const prepareFeatureLabelResultProjection = ({
  results,
  catalog,
  labelsByTargetKey = {},
  affectedResultKeys = [],
  targetFeatureKeysByResult = null,
  semanticScope = '',
  matchText = '',
  sourceText = '',
  mountedResultIndex = null,
  mountedSvg = null,
  mounted = null,
  preparedSvgByResultKey = null,
  parser = globalThis.DOMParser || globalThis.window?.DOMParser
} = {}) => {
  if (!Array.isArray(results)) throw new Error('Feature label projection requires Results.');
  const byResultKey = resultBindings(catalog);
  const indexValue = mountedResultIndex ?? mounted?.resultIndex;
  const activeResultIndex = Number.isInteger(Number(indexValue)) && Number(indexValue) >= 0
    ? Number(indexValue)
    : null;
  const activeSvg = mountedSvg || mounted?.svg || null;
  if (activeResultIndex !== null && activeResultIndex >= results.length) {
    throw new Error(`Mounted Feature label Result index is unavailable: ${activeResultIndex}`);
  }
  if (
    text(mounted?.resultKey)
    && catalogResultKey(catalog?.items?.[activeResultIndex]) !== text(mounted.resultKey)
  ) {
    throw new Error('Mounted Feature label Result identity is stale.');
  }
  const ordered = (Array.isArray(affectedResultKeys) ? affectedResultKeys : [])
    .map(text).filter(Boolean);
  const affected = new Set(ordered);
  if (affected.size === 0 || affected.size !== ordered.length) {
    throw new Error('Feature label projection requires unique affected Result identities.');
  }
  const candidates = new Map();
  const stagedRoots = new Map();
  const admissionMetadataByResultKey = {};
  const counters = {
    affectedResults: affected.size,
    mountedSerializations: 0,
    detachedPasses: 0,
    changedResults: 0,
    targetFeatures: 0,
    matchedLabels: 0,
    changedLabels: 0
  };

  affected.forEach((resultKey) => {
    const binding = byResultKey.get(resultKey);
    if (!binding) throw new Error(`Feature label target Result is unavailable: ${resultKey}`);
    const { resultIndex } = binding;
    const result = results[resultIndex];
    if (!result) throw new Error(`Feature label target Result index is unavailable: ${resultIndex}`);
    const featureKeys = targetFeatureKeysByResult instanceof Map
      ? targetFeatureKeysByResult.get(resultKey)
      : targetFeatureKeysByResult?.[resultKey];
    if (!Array.isArray(featureKeys) || featureKeys.length === 0) {
      throw new Error(`Feature label target coverage is empty for Result ${resultKey}.`);
    }
    const transform = (svg) => {
      const outcome = applyFeatureLabelProjectionToSvg({
        svg,
        semanticScope,
        labelsByTargetKey,
        targetFeatureKeys: featureKeys,
        matchText,
        sourceText
      });
      counters.targetFeatures += outcome.targetFeatures;
      counters.matchedLabels += outcome.matchedLabels;
      counters.changedLabels += outcome.changedLabels;
      stagedRoots.set(resultKey, svg);
    };
    const preparedRoot = preparedSvgByResultKey instanceof Map
      ? preparedSvgByResultKey.get(resultKey)
      : preparedSvgByResultKey?.[resultKey];
    let candidate;
    if (preparedRoot) {
      if (!preparedRoot?.cloneNode) {
        throw new Error(`Prepared Feature label SVG is invalid for Result ${resultKey}.`);
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
  const nextResults = results.map((result, index) => (
    candidates.has(index) ? candidates.get(index) : result
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

export const commitFeatureLabelResultProjection = (state, projection) => {
  if (!projection || state.results.value !== projection.previousResults) {
    throw new Error('Feature label Result projection is stale.');
  }
  state.results.value = projection.nextResults;
  return true;
};
