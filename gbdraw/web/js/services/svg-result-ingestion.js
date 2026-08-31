import { FEATURE_SELECTOR, filterFeatureFillTargets, getFeatureIdentity } from '../app/feature-dom.js';
import { getAllFeatureLegendGroups, parseTransformXY } from '../app/legend/utils.js';
import { isCurrentWorkerGenerationResponse } from './current-worker-result-source.js';
import { sanitizeSvgContent } from './svg-sanitization.js';
import { serializeCleanSvg } from './svg-serialization.js';
import { collectRenderedFeatureIdentitiesFromSvgRoot } from './session-feature-metadata.js';
import { normalizeSvgResultIds } from './svg-result-normalization.js';
import {
  recordSessionLifecycleEvent,
  recordStructuralMetric
} from './runtime-test-hooks.js';

// These symbols are intentionally module-private. A JSON/Session round trip can
// preserve Result data, but cannot reproduce runtime source or commit provenance.
const RESULT_SOURCE = Symbol('gbdraw.svgResultSource');
const COMMITTED_SVG_RESULT = Symbol('gbdraw.committedSvgResult');

const SVG_RESULT_SOURCE_CLASSES = Object.freeze({
  CURRENT_WORKER: 'current-worker',
  CURRENT_SESSION: 'current-session',
  LEGACY_IMPORT: 'legacy-import'
});

let nextResultIdentity = 1;

const text = (value) => String(value ?? '').trim();

const committedState = (result) => (
  result && typeof result === 'object' && result[COMMITTED_SVG_RESULT]
    ? result[COMMITTED_SVG_RESULT]
    : null
);

const requireResultList = (results) => {
  if (!Array.isArray(results)) {
    throw new Error('The diagram engine returned an invalid Result list.');
  }
  results.forEach((result) => {
    if (!result || typeof result !== 'object' || Array.isArray(result)) {
      throw new Error('The diagram engine returned an invalid SVG Result.');
    }
  });
  return results;
};

const requireAlignedCatalogAdmission = (admission, results) => {
  if (
    !admission
    || !Array.isArray(admission.resultNames)
    || !(admission.renderedTargetsByOverrideKey instanceof Map)
    || !(admission.resultIndexesByRenderedId instanceof Map)
    || !Array.isArray(admission.renderedIdentitiesByResult)
  ) {
    throw new Error('Current SVG Result admission requires an admitted feature catalog.');
  }
  const logicalResults = requireResultList(results);
  if (
    admission.resultNames.length !== logicalResults.length
    || admission.resultNames.some((name, index) => name !== text(logicalResults[index]?.name))
  ) {
    throw new Error('Current SVG Results do not align with the admitted feature catalog.');
  }
  return admission;
};

const createRuntimeSource = (sourceClass, results, catalogAdmission = null) => Object.freeze({
  [RESULT_SOURCE]: true,
  sourceClass,
  results,
  catalogAdmission
});

/**
 * Classify persisted current-session Results without granting current-worker
 * provenance. The catalog keeps this path parse-free and preserves lazy Worker use.
 */
export const createCurrentSessionResultSource = (results, catalogAdmission) => (
  createRuntimeSource(
    SVG_RESULT_SOURCE_CLASSES.CURRENT_SESSION,
    requireResultList(results),
    requireAlignedCatalogAdmission(catalogAdmission, results)
  )
);

/** Classify persisted historical input as the sole compatibility-normalization source. */
export const createLegacyImportResultSource = (results) => (
  createRuntimeSource(
    SVG_RESULT_SOURCE_CLASSES.LEGACY_IMPORT,
    requireResultList(results)
  )
);

const requireSource = (source, sourceClass) => {
  if (!source?.[RESULT_SOURCE] || source.sourceClass !== sourceClass) {
    throw new Error(`SVG Result source must be classified as ${sourceClass}.`);
  }
  return source;
};

const parseSanitizedSvg = (
  content,
  parser = globalThis.DOMParser || globalThis.window?.DOMParser
) => {
  if (typeof parser !== 'function') {
    throw new Error('SVG parsing is unavailable.');
  }
  const document = new parser().parseFromString(content, 'image/svg+xml');
  if (
    document?.querySelector?.('parsererror')
    || String(document?.documentElement?.localName || '').toLowerCase() !== 'svg'
  ) {
    throw new Error('The diagram engine returned malformed SVG content.');
  }
  return document.documentElement;
};

const markCommitted = (result, metadata) => {
  result[COMMITTED_SVG_RESULT] = {
    identity: nextResultIdentity++,
    mounted: false,
    mountedContent: result.content,
    metadata
  };
  return result;
};

const commitCatalogBackedResult = (result, metadata) => markCommitted(result, metadata);

const hasSanitizedSvgEnvelope = (content) => {
  const withoutDeclaration = String(content || '')
    .trim()
    .replace(/^<\?xml[^>]*>\s*/i, '');
  return /^<svg(?:\s|>)/i.test(withoutDeclaration);
};

const sanitizeSvgResultContent = (
  result,
  sanitizer,
  { phase, resultIndex }
) => {
  recordSessionLifecycleEvent('result-svg-characters', {
    phase,
    resultIndex,
    value: String(result.content || '').length
  });
  const content = sanitizeSvgContent(result.content, sanitizer);
  recordStructuralMetric('svgSanitizationCount', 1, { phase, resultIndex });
  return content;
};

const serializeAdmittedSvg = (svg, { phase, resultIndex }) => {
  const content = serializeCleanSvg(svg);
  recordStructuralMetric('svgSerializationCount', 1, { phase, resultIndex });
  return content;
};

const metadataFromCatalogAdmission = (catalogAdmission, resultIndex, sourceClass) => {
  const renderedFeatureIdentities = catalogAdmission.renderedIdentitiesByResult[resultIndex];
  if (!renderedFeatureIdentities) {
    throw new Error('The admitted feature catalog is missing Result identity metadata.');
  }
  return Object.freeze({ renderedFeatureIdentities, sourceClass });
};

export const isCommittedSvgResult = (result) => Boolean(committedState(result));

export const isCommittedSvgResultMounted = (result) => Boolean(
  committedState(result)?.mounted
);

export const getCommittedSvgContent = (result) => {
  const runtime = committedState(result);
  if (!runtime) return null;
  return runtime.mounted ? runtime.mountedContent : result.content;
};

export const getCommittedSvgResultMetadata = (result) => (
  committedState(result)?.metadata || null
);

export const markCommittedSvgResultMounted = (result) => {
  const runtime = committedState(result);
  if (!runtime) return false;
  runtime.mounted = true;
  runtime.mountedContent = result.content;
  return true;
};

export const markCommittedSvgResultUnmounted = (result) => {
  const runtime = committedState(result);
  if (!runtime) return false;
  runtime.mounted = false;
  runtime.mountedContent = result.content;
  return true;
};

const currentResultOperations = (plan, resultIndex) => (
  plan.operationsByResult[resultIndex] || null
);

const hasOperations = (operations) => Boolean(
  operations
  && (
    operations.featureFills.length
    || operations.featureStrokes.length
    || operations.featureVisibility.length
    || operations.labelText.length
    || operations.labelVisibility.length
    || operations.legendFills.length
    || operations.legendStrokes.length
    || operations.legendRenames.length
    || operations.legendDeletes.length
    || operations.legendAdds.length
    || operations.callerTransforms.length
  )
);

const requireCurrentMutationPlan = (plan, resultCount) => {
  if (
    !plan
    || (plan.kind !== 'EMPTY' && plan.kind !== 'MUTATING')
    || !Array.isArray(plan.operationsByResult)
    || plan.operationsByResult.length !== resultCount
  ) {
    throw new Error('Current SVG Result admission requires a valid mutation plan.');
  }
  const mutating = plan.operationsByResult.some(hasOperations);
  if ((plan.kind === 'MUTATING') !== mutating) {
    throw new Error('Current SVG Result mutation plan classification is inconsistent.');
  }
  return plan;
};

const freezeEmptyOperations = () => Object.freeze({
  featureFills: Object.freeze([]),
  featureStrokes: Object.freeze([]),
  featureVisibility: Object.freeze([]),
  labelText: Object.freeze([]),
  labelVisibility: Object.freeze([]),
  legendFills: Object.freeze([]),
  legendStrokes: Object.freeze([]),
  legendRenames: Object.freeze([]),
  legendDeletes: Object.freeze([]),
  legendAdds: Object.freeze([]),
  callerTransforms: Object.freeze([])
});

export const createEmptySvgMutationPlan = (resultCount) => {
  if (!Number.isSafeInteger(resultCount) || resultCount < 0) {
    throw new TypeError('An EMPTY SVG mutation plan requires a nonnegative Result count.');
  }
  return Object.freeze({
    kind: 'EMPTY',
    operationsByResult: Object.freeze(
      Array.from({ length: resultCount }, () => freezeEmptyOperations())
    )
  });
};

const setAttributeIfDifferent = (element, name, value) => {
  const normalized = String(value);
  if (element.getAttribute(name) === normalized) return false;
  element.setAttribute(name, normalized);
  return true;
};

const removeAttributeIfPresent = (element, name) => {
  if (!element.hasAttribute(name)) return false;
  element.removeAttribute(name);
  return true;
};

const legendSwatch = (entryGroup) => Array.from(
  entryGroup?.querySelectorAll?.('path') || []
).find((path) => {
  const fill = path.getAttribute('fill');
  return fill && fill !== 'none' && !fill.startsWith('url(');
}) || null;

const createLazyMutationIndex = (svg, { phase, resultIndex }) => {
  const state = {
    featureElements: null,
    labelElements: null,
    legendEntries: null,
    legendGroups: null
  };
  let announced = false;
  const announce = () => {
    if (announced) return;
    announced = true;
    recordStructuralMetric('svgMutationIndexBuildCount', 1, { phase, resultIndex });
  };
  return {
    features() {
      announce();
      if (state.featureElements) return state.featureElements;
      state.featureElements = new Map();
      Array.from(svg.querySelectorAll(FEATURE_SELECTOR)).forEach((element) => {
        const renderedId = getFeatureIdentity(element);
        if (!renderedId) return;
        if (!state.featureElements.has(renderedId)) state.featureElements.set(renderedId, []);
        state.featureElements.get(renderedId).push(element);
      });
      recordStructuralMetric('featureDomFullScanCount', 1, { phase, resultIndex });
      return state.featureElements;
    },
    labels() {
      announce();
      if (state.labelElements) return state.labelElements;
      state.labelElements = new Map();
      Array.from(svg.querySelectorAll('text[data-label-feature-id]')).forEach((element) => {
        const renderedId = text(element.getAttribute('data-label-feature-id'));
        if (!renderedId) return;
        if (!state.labelElements.has(renderedId)) state.labelElements.set(renderedId, []);
        state.labelElements.get(renderedId).push(element);
      });
      return state.labelElements;
    },
    legends() {
      announce();
      if (state.legendEntries) return {
        entries: state.legendEntries,
        groups: state.legendGroups
      };
      state.legendEntries = new Map();
      state.legendGroups = getAllFeatureLegendGroups(svg);
      state.legendGroups.forEach((group) => {
        const seen = new Set();
        Array.from(group.querySelectorAll('g[data-legend-key]')).forEach((entry) => {
          const caption = text(entry.getAttribute('data-legend-key'));
          if (!caption || seen.has(caption)) {
            throw new Error('Current SVG contains an ambiguous Legend binding.');
          }
          seen.add(caption);
          if (!state.legendEntries.has(caption)) state.legendEntries.set(caption, []);
          state.legendEntries.get(caption).push(entry);
        });
      });
      recordStructuralMetric('legendDomFullScanCount', 1, { phase, resultIndex });
      return { entries: state.legendEntries, groups: state.legendGroups };
    }
  };
};

const requireFeatureElements = (index, renderedId) => {
  const elements = index.features().get(renderedId) || [];
  if (elements.length === 0) {
    throw new Error('Sanitized SVG content is missing a rendered Feature binding.');
  }
  return elements;
};

const requireLabelElement = (index, renderedId) => {
  const elements = index.labels().get(renderedId) || [];
  if (elements.length !== 1) {
    throw new Error('Sanitized SVG content is missing or ambiguously binds an editable Label.');
  }
  return elements[0];
};

const requireLegendEntries = (index, caption) => {
  const entries = index.legends().entries.get(caption) || [];
  if (entries.length === 0) {
    throw new Error('Sanitized SVG content is missing a Legend binding.');
  }
  return entries;
};

const applyFeatureOperations = (index, operations) => {
  operations.featureFills.forEach(({ renderedId, color }) => {
    const targets = filterFeatureFillTargets(requireFeatureElements(index, renderedId));
    if (targets.length === 0) {
      throw new Error('Sanitized SVG content is missing a rendered Feature fill target.');
    }
    targets.forEach((element) => setAttributeIfDifferent(element, 'fill', color));
  });
  operations.featureStrokes.forEach(({ renderedId, strokeColor, strokeWidth }) => {
    requireFeatureElements(index, renderedId).forEach((element) => {
      if (strokeColor) setAttributeIfDifferent(element, 'stroke', strokeColor);
      if (strokeWidth !== null) setAttributeIfDifferent(element, 'stroke-width', strokeWidth);
    });
  });
  operations.featureVisibility.forEach(({ renderedId, mode }) => {
    requireFeatureElements(index, renderedId).forEach((element) => {
      if (mode === 'off') setAttributeIfDifferent(element, 'display', 'none');
      else removeAttributeIfPresent(element, 'display');
    });
  });
};

const setLabelText = (element, value) => {
  const textPath = element.querySelector?.('textPath');
  if (textPath) textPath.textContent = value;
  else element.textContent = value;
};

const applyLabelOperations = (index, operations) => {
  operations.labelText.forEach(({ renderedId, value }) => {
    setLabelText(requireLabelElement(index, renderedId), value);
  });
  operations.labelVisibility.forEach(({ renderedId, mode }) => {
    const element = requireLabelElement(index, renderedId);
    if (mode === 'off') {
      setAttributeIfDifferent(element, 'data-gbdraw-label-visibility-preview', 'off');
      setAttributeIfDifferent(element, 'display', 'none');
    } else {
      removeAttributeIfPresent(element, 'data-gbdraw-label-visibility-preview');
      removeAttributeIfPresent(element, 'display');
    }
  });
};

const updateLegendCaption = (entry, caption) => {
  entry.setAttribute('data-legend-key', caption);
  const label = entry.querySelector('text');
  if (label) label.textContent = caption;
};

const legendEntryAnchor = (entry) => {
  const groupOffset = parseTransformXY(entry?.getAttribute?.('transform'));
  const target = entry?.querySelector?.('text') || legendSwatch(entry);
  const targetOffset = parseTransformXY(target?.getAttribute?.('transform'));
  return { x: groupOffset.x + targetOffset.x, y: groupOffset.y + targetOffset.y };
};

const moveLegendEntryToAnchor = (entry, xPos, yPos) => {
  if (!Number.isFinite(xPos) || !Number.isFinite(yPos)) return;
  const current = legendEntryAnchor(entry);
  const deltaX = xPos - current.x;
  const deltaY = yPos - current.y;
  if (Math.abs(deltaX) < 1e-6 && Math.abs(deltaY) < 1e-6) return;
  if (entry.hasAttribute?.('transform')) {
    const groupOffset = parseTransformXY(entry.getAttribute('transform'));
    entry.setAttribute(
      'transform',
      `translate(${groupOffset.x + deltaX}, ${groupOffset.y + deltaY})`
    );
    return;
  }
  const transformedChildren = Array.from(entry.querySelectorAll?.('[transform]') || []);
  if (transformedChildren.length === 0) {
    entry.setAttribute('transform', `translate(${deltaX}, ${deltaY})`);
    return;
  }
  transformedChildren.forEach((node) => {
    const position = parseTransformXY(node.getAttribute('transform'));
    node.setAttribute('transform', `translate(${position.x + deltaX}, ${position.y + deltaY})`);
  });
};

const applyLegendOperations = (index, operations) => {
  operations.legendFills.forEach(({ caption, color }) => {
    requireLegendEntries(index, caption).forEach((entry) => {
      const swatch = legendSwatch(entry);
      if (!swatch) throw new Error('Sanitized SVG content is missing a Legend swatch.');
      setAttributeIfDifferent(swatch, 'fill', color);
    });
  });
  operations.legendStrokes.forEach(({ caption, strokeColor, strokeWidth, renderedIds }) => {
    (Array.isArray(renderedIds) ? renderedIds : []).forEach((renderedId) => {
      requireFeatureElements(index, renderedId).forEach((element) => {
        if (strokeColor) setAttributeIfDifferent(element, 'stroke', strokeColor);
        if (strokeWidth !== null) setAttributeIfDifferent(element, 'stroke-width', strokeWidth);
      });
    });
    requireLegendEntries(index, caption).forEach((entry) => {
      const swatch = legendSwatch(entry);
      if (!swatch) throw new Error('Sanitized SVG content is missing a Legend swatch.');
      if (strokeColor) setAttributeIfDifferent(swatch, 'stroke', strokeColor);
      if (strokeWidth !== null) setAttributeIfDifferent(swatch, 'stroke-width', strokeWidth);
    });
  });
  operations.legendRenames.forEach(({ from, to, xPos, yPos }) => {
    requireLegendEntries(index, from).forEach((entry) => {
      updateLegendCaption(entry, to);
      moveLegendEntryToAnchor(entry, xPos, yPos);
    });
  });
  operations.legendDeletes.forEach(({ caption }) => {
    requireLegendEntries(index, caption).forEach((entry) => entry.remove());
  });
  operations.legendAdds.forEach(({ caption, color, xPos, yPos }) => {
    const { entries, groups } = index.legends();
    if (entries.has(caption) || groups.length === 0) {
      throw new Error('Current SVG cannot admit the requested Legend addition.');
    }
    groups.forEach((group) => {
      const template = group.querySelector('g[data-legend-key]');
      const added = template?.cloneNode?.(true) || null;
      if (!added) throw new Error('Current SVG has no Legend entry template.');
      updateLegendCaption(added, caption);
      const swatch = legendSwatch(added);
      if (!swatch) throw new Error('Current SVG has no Legend swatch template.');
      swatch.setAttribute('fill', color);
      added.setAttribute('data-legend-owner', 'direct-editor');
      moveLegendEntryToAnchor(added, xPos, yPos);
      group.appendChild(added);
    });
  });
};

const admitCurrentResult = (
  result,
  metadata,
  operations,
  {
    sanitizer,
    parser,
    resultIndex,
    sourceClass
  }
) => {
  const phase = sourceClass;
  const sanitized = sanitizeSvgResultContent(result, sanitizer, { phase, resultIndex });
  if (!hasSanitizedSvgEnvelope(sanitized)) {
    throw new Error('The diagram engine returned malformed SVG content.');
  }
  if (!hasOperations(operations)) {
    return commitCatalogBackedResult({ ...result, content: sanitized }, metadata);
  }

  const svg = parseSanitizedSvg(sanitized, parser);
  recordStructuralMetric('applicationSvgParseCount', 1, { phase, resultIndex });
  const index = createLazyMutationIndex(svg, { phase, resultIndex });
  applyFeatureOperations(index, operations);
  applyLabelOperations(index, operations);
  applyLegendOperations(index, operations);
  operations.callerTransforms.forEach((transform) => transform(svg, { result, resultIndex }));
  const content = serializeAdmittedSvg(svg, { phase, resultIndex });
  return commitCatalogBackedResult({ ...result, content }, metadata);
};

const admitCatalogBackedResults = (
  source,
  mutationPlan,
  {
    sanitizer = globalThis.DOMPurify || globalThis.window?.DOMPurify,
    parser = globalThis.DOMParser || globalThis.window?.DOMParser
  } = {}
) => {
  const { results, catalogAdmission, sourceClass } = source;
  requireAlignedCatalogAdmission(catalogAdmission, results);
  const plan = requireCurrentMutationPlan(mutationPlan, results.length);
  recordSessionLifecycleEvent('svg.admission-started', {
    phase: sourceClass,
    resultCount: results.length,
    mutationKind: plan.kind
  });
  recordStructuralMetric('currentLegacyNormalizationCount', 0, { phase: sourceClass });
  recordStructuralMetric('legacyOverrideMigrationCount', 0, { phase: sourceClass });
  recordStructuralMetric('manualRuleFeatureMatchCount', 0, { phase: sourceClass });
  const admitted = results.map((result, resultIndex) => admitCurrentResult(
    result,
    metadataFromCatalogAdmission(catalogAdmission, resultIndex, sourceClass),
    currentResultOperations(plan, resultIndex),
    { sanitizer, parser, resultIndex, sourceClass }
  ));
  recordSessionLifecycleEvent('svg.admission-completed', {
    phase: sourceClass,
    resultCount: admitted.length,
    mutationKind: plan.kind
  });
  recordSessionLifecycleEvent('artifact.candidate-completed', {
    phase: sourceClass,
    resultCount: admitted.length,
    catalogFootprint: catalogAdmission.scalarMetrics
  });
  return admitted;
};

/**
 * Admit only a freshly decoded current Worker response. The runtime token minted
 * by diagram-generation.js and exact catalog object alignment are both required.
 */
export const admitCurrentGeneratedResults = (
  generationResponse,
  {
    catalogAdmission,
    mutationPlan,
    sanitizer,
    parser
  } = {}
) => {
  if (!isCurrentWorkerGenerationResponse(generationResponse)) {
    throw new Error('Current SVG Results require runtime Worker provenance.');
  }
  if (generationResponse.metadata?.featureCatalog !== catalogAdmission?.catalog) {
    throw new Error('Current Worker SVG Results do not own the admitted feature catalog.');
  }
  const source = createRuntimeSource(
    SVG_RESULT_SOURCE_CLASSES.CURRENT_WORKER,
    requireResultList(generationResponse.results),
    requireAlignedCatalogAdmission(catalogAdmission, generationResponse.results)
  );
  return admitCatalogBackedResults(source, mutationPlan, { sanitizer, parser });
};

export const admitCurrentSessionResults = (
  source,
  {
    mutationPlan,
    sanitizer,
    parser
  } = {}
) => admitCatalogBackedResults(
  requireSource(source, SVG_RESULT_SOURCE_CLASSES.CURRENT_SESSION),
  mutationPlan,
  { sanitizer, parser }
);

const ingestSvgResult = (
  result,
  {
    sanitizer,
    parser,
    transformSvg,
    sourceClass,
    resultIndex
  }
) => {
  if (isCommittedSvgResult(result)) return result;
  const content = sanitizeSvgResultContent(result, sanitizer, {
    phase: sourceClass,
    resultIndex
  });
  const svg = parseSanitizedSvg(content, parser);
  recordStructuralMetric('applicationSvgParseCount', 1, {
    phase: sourceClass,
    resultIndex
  });
  if (typeof transformSvg === 'function') transformSvg(svg, { result, resultIndex });
  normalizeSvgResultIds(svg);
  const renderedFeatureIdentities = collectRenderedFeatureIdentitiesFromSvgRoot(svg);
  recordStructuralMetric('svgIdentityScanCount', 1, { phase: sourceClass, resultIndex });
  const serialized = serializeAdmittedSvg(svg, { phase: sourceClass, resultIndex });
  return markCommitted(
    { ...result, content: serialized },
    Object.freeze({ renderedFeatureIdentities, sourceClass })
  );
};

/**
 * Admit historical/unclassified persisted SVG through compatibility repair.
 * Current Worker and current-session sources cannot enter this function.
 */
export const admitLegacyImportedResults = (
  source,
  {
    sanitizer = globalThis.DOMPurify || globalThis.window?.DOMPurify,
    parser = globalThis.DOMParser || globalThis.window?.DOMParser,
    transformSvg = null
  } = {}
) => {
  const { results, sourceClass } = requireSource(
    source,
    SVG_RESULT_SOURCE_CLASSES.LEGACY_IMPORT
  );
  recordSessionLifecycleEvent('svg.admission-started', {
    phase: sourceClass,
    resultCount: results.length
  });
  const admitted = results.map((result, resultIndex) => ingestSvgResult(result, {
    sanitizer,
    parser,
    transformSvg,
    sourceClass,
    resultIndex
  }));
  recordSessionLifecycleEvent('svg.admission-completed', {
    phase: sourceClass,
    resultCount: admitted.length
  });
  return admitted;
};
