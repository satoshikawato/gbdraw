// @ts-check
/** @import { FeatureCatalogAdmission } from './feature-catalog.js' */
import {
  AUTO_FEATURE_UNDERLAY_STROKE,
  FEATURE_PART_CONNECTOR,
  FEATURE_SELECTOR,
  filterFeatureFillTargets,
  getFeatureIdentity,
  getFeaturePart,
  isAutoFeatureUnderlay
} from './feature-dom.js';
import {
  getAllFeatureLegendGroups,
  getLegendEntrySwatch as legendSwatch,
  legendRowFeatureIds,
  moveLegendEntryToAnchor,
  orderLegendEntries,
  pythonLegendKey,
  setsFeatureStroke,
  SPECIFIC_COLOR_FILE_OWNER
} from './legend-svg.js';
import { isCurrentWorkerGenerationResponse } from './current-worker-result-source.js';
import { diagnosticError } from '../utils/error-normalization.js';
import { sanitizeSvgContent } from './svg-sanitization.js';
import { LEGEND_ORDER_RECORD, RESULT_BASE_SELECTOR, resultBaseAttribute } from './result-paint-bases.js';
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

/**
 * A caller's rewrite of one parsed Result SVG. The return value is ignored.
 * @typedef {(svg: Element, context: { result: any, resultIndex: number }) => unknown} SvgResultTransform
 */

/**
 * One Result's compiled editor operations. The planner is app/candidate-render.js;
 * this module declares the shape it applies.
 * A live feature fill of `color: null` keeps the fill the Result shows.
 * @typedef {Record<
 *   'featureFills' | 'featureStrokes' | 'featureVisibility' | 'labelText' | 'labelVisibility'
 *   | 'legendFills' | 'legendStrokes' | 'legendRenames' | 'legendDeletes' | 'legendAdds' | 'legendOrder',
 *   readonly Record<string, any>[]
 * > & { callerTransforms: readonly SvgResultTransform[] }} SvgMutationOperations
 */

/**
 * @typedef {object} SvgMutationPlan
 * @property {'EMPTY' | 'MUTATING'} kind
 * @property {readonly SvgMutationOperations[]} operationsByResult One entry per Result.
 * @property {number} [legacyNormalizationCount] The Results whose transform is a
 *   legacy normalization (`createSavedResultPlan`).
 */

/**
 * A classified Result list: the output of `createCurrentSessionResultSource` or
 * `createLegacyImportResultSource`, the only input of the session and legacy
 * admissions.
 * @typedef {object} SvgResultSource
 * @property {string} sourceClass
 * @property {Record<string, any>[]} results
 * @property {FeatureCatalogAdmission | null} catalogAdmission Set for a current source only.
 * @property {readonly SvgLegendRowFacts[] | null} [legendRows] Python's Legend row facts, one per Result; set for a current generated source only.
 */

/**
 * Python's Legend row facts for one Result: the row keys its Legend drew, and the
 * feature rows the draft removed there.
 * @typedef {object} SvgLegendRowFacts
 * @property {ReadonlySet<string>} drawn
 * @property {ReadonlySet<string>} suppressed
 */

/**
 * The DOM services and request facts an admission reads. A missing sanitizer
 * or parser falls back to the page's DOMPurify and DOMParser.
 * @typedef {object} SvgAdmissionRuntime
 * @property {any} [sanitizer] DOMPurify, or an object with `sanitize(svg, options)`.
 * @property {any} [parser] A DOMParser constructor.
 * @property {readonly string[] | null} [selectedFeatureTypes] The feature types of the request that drew the Results.
 */

/**
 * @typedef {SvgAdmissionRuntime & {
 *   catalogAdmission: FeatureCatalogAdmission,
 *   mutationPlan: SvgMutationPlan
 * }} CurrentGeneratedAdmissionOptions
 */

/** @typedef {SvgAdmissionRuntime & { mutationPlan: SvgMutationPlan }} CurrentSessionAdmissionOptions */

/**
 * @typedef {SvgAdmissionRuntime & {
 *   transformSvg?: SvgResultTransform | null
 * }} LegacyImportAdmissionOptions
 */

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

/**
 * @param {string} sourceClass
 * @param {Record<string, any>[]} results
 * @param {FeatureCatalogAdmission | null} [catalogAdmission]
 * @param {readonly SvgLegendRowFacts[] | null} [legendRows]
 * @returns {SvgResultSource}
 */
const createRuntimeSource = (
  sourceClass,
  results,
  catalogAdmission = null,
  legendRows = null
) => Object.freeze({
  [RESULT_SOURCE]: true,
  sourceClass,
  results,
  catalogAdmission,
  legendRows
});

const invalidResult = () => diagnosticError('RESULT_INVALID', {}, { stage: 'result-admission' });

// Python's Legend row facts for each Result (`metadata.legendRows`): the row keys
// its Legend drew, and the feature rows its records can name that the draft did not
// draw (hidden or recaptioned features). They are read for one admission and never
// persisted. Without them (a replayed Session, a stub) a required row stays required.
const normalizeLegendRowFacts = (legendRows, results) => {
  if (legendRows === undefined) return null;
  const keys = (list) => (
    Array.isArray(list) && list.every((key) => typeof key === 'string') ? new Set(list) : null
  );
  if (!Array.isArray(legendRows) || legendRows.length !== results.length) throw invalidResult();
  return Object.freeze(legendRows.map((row, index) => {
    const drawn = keys(row?.drawn);
    const suppressed = keys(row?.suppressed);
    if (!drawn || !suppressed || row.resultIndex !== index || row.resultName !== text(results[index]?.name)) {
      throw invalidResult();
    }
    return Object.freeze({ drawn, suppressed });
  }));
};

// A required Legend row that Result `resultIndex` does not contain may be absent when
// Python reports it as removed by the draft there (suppressed), or as belonging to
// another Result of the batch (drawn or suppressed there). A row Python drew in this
// Result must be present, and a key no Result reports is stale (OV-46, OV-63).
const legendRowMayBeAbsent = (legendRows, resultIndex, caption) => {
  if (!legendRows) return false;
  const own = legendRows[resultIndex];
  if (own.drawn.has(caption)) return false;
  if (own.suppressed.has(caption)) return true;
  return legendRows.some((other, index) => index !== resultIndex
    && (other.drawn.has(caption) || other.suppressed.has(caption)));
};

/**
 * Classify persisted current-session Results without granting current-worker
 * provenance. The catalog keeps this path parse-free and preserves lazy Worker use.
 */
/**
 * @param {Record<string, any>[]} results
 * @param {FeatureCatalogAdmission} catalogAdmission
 * @returns {SvgResultSource}
 */
export const createCurrentSessionResultSource = (results, catalogAdmission) => (
  createRuntimeSource(
    SVG_RESULT_SOURCE_CLASSES.CURRENT_SESSION,
    requireResultList(results),
    requireAlignedCatalogAdmission(catalogAdmission, results)
  )
);

/** Classify persisted historical input as the sole compatibility-normalization source. */
/**
 * @param {Record<string, any>[]} results
 * @returns {SvgResultSource}
 */
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

// The feature types of the request that drew the Results (its
// `diagramOptions.selectedFeaturesSet`), or null when it names none. The
// Features list reads them from the displayed Result's metadata (R13).
const resultFeatureTypes = (selectedFeatureTypes) => (
  Array.isArray(selectedFeatureTypes) ? Object.freeze(selectedFeatureTypes.map(String)) : null
);

const metadataFromCatalogAdmission = (catalogAdmission, resultIndex, sourceClass, selectedFeatureTypes) => {
  const renderedFeatureIdentities = catalogAdmission.renderedIdentitiesByResult[resultIndex];
  if (!renderedFeatureIdentities) {
    throw new Error('The admitted feature catalog is missing Result identity metadata.');
  }
  return Object.freeze({ renderedFeatureIdentities, sourceClass, selectedFeatureTypes });
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

export const getCommittedSvgResultRuntimeIdentity = (result) => (
  committedState(result)?.identity ?? null
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
    || operations.legendOrder.length
    || operations.callerTransforms.length
  )
);

const hasDetachedOperations = (operations) => Boolean(
  operations
  && (
    operations.featureFills.length
    || operations.featureStrokes.length
    || operations.featureVisibility.length
    || operations.legendFills.length
    || operations.legendStrokes.length
    || operations.legendRenames.length
    || operations.legendDeletes.length
    || operations.legendAdds.length
    || operations.legendOrder.length
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
  legendOrder: Object.freeze([]),
  callerTransforms: Object.freeze([])
});

/**
 * @param {number} resultCount
 * @returns {SvgMutationPlan}
 */
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

// An operation that keeps what the Result shows (a live value not known yet)
// still claims the attribute, so a reconcile leaves it as it is.
/**
 * @param {{ painted: Map<Element, Set<string>> }} index
 * @param {Element} element
 * @param {string} name
 */
const keepPaintAttribute = (index, element, name) => {
  index.painted.set(element, (index.painted.get(element) || new Set()).add(name));
};

// A paint attribute the executor changes keeps Python's value beside it
// (`resultBaseAttribute`), recorded on the first change. The index notes
// each attribute an operation sets, so a reconcile leaves it alone.
/**
 * @param {{ painted: Map<Element, Set<string>> }} index
 * @param {Element} element
 * @param {string} name
 * @param {string | number | null} value null removes the attribute.
 */
const setPaintAttribute = (index, element, name, value) => {
  keepPaintAttribute(index, element, name);
  const current = element.getAttribute(name);
  const next = value === null ? null : String(value);
  if (current === next) return false;
  const base = resultBaseAttribute(name);
  if (!element.hasAttribute(base)) element.setAttribute(base, current ?? '');
  if (next === null) element.removeAttribute(name);
  else element.setAttribute(name, next);
  return true;
};

// The domains a reconcile returns to Python's values, and the attributes
// each one owns on feature elements, on Legend swatches, and on Legend rows
// (a rename's key and a delete's hiding). A Legend row stroke also strokes
// the row's features.
const PAINT_DOMAIN_ATTRIBUTES = Object.freeze({
  featureFills: { feature: ['fill'], swatch: [], row: [] },
  featureStrokes: { feature: ['stroke', 'stroke-width'], swatch: [], row: [] },
  featureVisibility: { feature: ['display'], swatch: [], row: [] },
  legendFills: { feature: [], swatch: ['fill'], row: [] },
  legendStrokes: { feature: ['stroke', 'stroke-width'], swatch: ['stroke', 'stroke-width'], row: [] },
  legendRenames: { feature: [], swatch: [], row: ['data-legend-key'] },
  legendDeletes: { feature: [], swatch: [], row: ['display'] }
});
const RESULT_PAINT_DOMAINS = Object.freeze(
  /** @type {Array<keyof typeof PAINT_DOMAIN_ATTRIBUTES>} */ (Object.keys(PAINT_DOMAIN_ATTRIBUTES))
);

/** @param {Element} element @returns {'feature' | 'swatch' | 'row'} */
const paintElementKind = (element) => {
  if (element.hasAttribute('data-legend-key')) return 'row';
  for (let node = element.parentElement; node; node = node.parentElement) {
    if (node.hasAttribute('data-legend-key')) return 'swatch';
  }
  return 'feature';
};

/** @param {Element} entry @param {string} caption */
const updateLegendCaption = (entry, caption) => {
  entry.setAttribute('data-legend-key', caption);
  const label = entry.querySelector('text');
  if (label) label.textContent = caption;
};

/**
 * Return every attribute the executor changed in `domains` and no operation
 * of this pass set to the value Python drew (a Legend row's key with its
 * text). Returns whether a Legend row changed.
 * @param {Element} root
 * @param {readonly string[]} domains
 * @param {Map<Element, Set<string>>} painted
 */
const restorePaintBases = (root, domains, painted) => {
  const owned = { feature: new Set(), swatch: new Set(), row: new Set() };
  domains.forEach((domain) => {
    const attributes = PAINT_DOMAIN_ATTRIBUTES[/** @type {keyof typeof PAINT_DOMAIN_ATTRIBUTES} */ (domain)];
    attributes?.feature.forEach((name) => owned.feature.add(name));
    attributes?.swatch.forEach((name) => owned.swatch.add(name));
    attributes?.row.forEach((name) => owned.row.add(name));
  });
  let rowsChanged = false;
  if (owned.feature.size === 0 && owned.swatch.size === 0 && owned.row.size === 0) return rowsChanged;
  [root, ...Array.from(root.querySelectorAll(RESULT_BASE_SELECTOR))].forEach((element) => {
    const kind = paintElementKind(element);
    owned[kind].forEach((name) => {
      if (painted.get(element)?.has(name)) return;
      const base = resultBaseAttribute(name);
      const value = element.getAttribute(base);
      if (value === null) return;
      if (value === '') element.removeAttribute(name);
      else if (name === 'data-legend-key') updateLegendCaption(element, value);
      else element.setAttribute(name, value);
      element.removeAttribute(base);
      if (kind === 'row') rowsChanged = true;
    });
  });
  return rowsChanged;
};

const createLazyMutationIndex = (svg, { phase, resultIndex }) => {
  /** @type {{ featureElements: Map<string, Element[]> | null, legendEntries: Map<string, Element[]> | null, legendGroups: Element[] | null }} */
  const built = {
    featureElements: null,
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
    /** @type {Map<Element, Set<string>>} */
    painted: new Map(),
    features() {
      announce();
      if (built.featureElements) return built.featureElements;
      const featureElements = new Map();
      built.featureElements = featureElements;
      Array.from(svg.querySelectorAll(FEATURE_SELECTOR)).forEach((element) => {
        const renderedId = getFeatureIdentity(element);
        if (!renderedId) return;
        if (!featureElements.has(renderedId)) featureElements.set(renderedId, []);
        // `get` holds the array set on the line above.
        /** @type {Element[]} */ (featureElements.get(renderedId)).push(element);
      });
      recordStructuralMetric('featureDomFullScanCount', 1, { phase, resultIndex });
      return built.featureElements;
    },
    legends() {
      announce();
      if (built.legendEntries) return {
        entries: built.legendEntries,
        // `legendGroups` is set together with `legendEntries` below.
        groups: /** @type {Element[]} */ (built.legendGroups)
      };
      const legendEntries = new Map();
      built.legendEntries = legendEntries;
      built.legendGroups = getAllFeatureLegendGroups(svg);
      // A row of Python's is found by Python's key, also once renamed or
      // hidden; an editor row by its own key, unless a row of Python's has it.
      built.legendGroups.forEach((group) => {
        /** @type {Map<string, Element>} */
        const seen = new Map();
        const rows = Array.from(group.querySelectorAll('g[data-legend-key]'));
        [...rows.filter((entry) => !isEditorRow(entry)), ...rows.filter(isEditorRow)].forEach((entry) => {
          const caption = isEditorRow(entry) ? text(entry.getAttribute('data-legend-key')) : pythonLegendKey(entry);
          const other = seen.get(caption);
          if (other && !isEditorRow(other) && isEditorRow(entry)) return;
          if (!caption || other) {
            throw new Error('Current SVG contains an ambiguous Legend binding.');
          }
          seen.set(caption, entry);
          if (!legendEntries.has(caption)) legendEntries.set(caption, []);
          // `get` holds the array set on the line above.
          /** @type {Element[]} */ (legendEntries.get(caption)).push(entry);
        });
      });
      recordStructuralMetric('legendDomFullScanCount', 1, { phase, resultIndex });
      return { entries: built.legendEntries, groups: built.legendGroups };
    }
  };
};

/** @param {Element} entry */
const isEditorRow = (entry) => entry.getAttribute('data-legend-owner') === 'direct-editor';

const requireFeatureElements = (index, renderedId) => {
  const elements = index.features().get(renderedId) || [];
  if (elements.length === 0) {
    throw new Error('Sanitized SVG content is missing a rendered Feature binding.');
  }
  return elements;
};

/**
 * @param {any} index
 * @param {string} caption
 * @param {{ allowMissing?: boolean }} [options]
 * @param {(caption: string) => boolean} [mayBeAbsent]
 */
const requireLegendEntries = (index, caption, { allowMissing = false } = {}, mayBeAbsent = () => false) => {
  const entries = index.legends().entries.get(caption) || [];
  if (entries.length === 0 && !allowMissing && !mayBeAbsent(caption)) throw invalidResult();
  return entries;
};

const applyFeatureOperations = (index, operations) => {
  operations.featureFills.forEach(({ renderedId, color }) => {
    const targets = filterFeatureFillTargets(requireFeatureElements(index, renderedId));
    if (targets.length === 0) {
      throw new Error('Sanitized SVG content is missing a rendered Feature fill target.');
    }
    targets.forEach((element) => (color === null
      ? keepPaintAttribute(index, element, 'fill')
      : setPaintAttribute(index, element, 'fill', color)));
  });
  operations.featureStrokes.forEach(({ renderedId, strokeColor, strokeWidth }) => {
    requireFeatureElements(index, renderedId).forEach((element) => {
      if (strokeColor) setPaintAttribute(index, element, 'stroke', strokeColor);
      if (strokeWidth !== null) setPaintAttribute(index, element, 'stroke-width', strokeWidth);
    });
  });
  operations.featureVisibility.forEach(({ renderedId, mode }) => {
    requireFeatureElements(index, renderedId).forEach((element) => {
      setPaintAttribute(index, element, 'display', mode === 'off' ? 'none' : null);
    });
  });
};

// `mayBeAbsent(caption)` says whether Python's Legend row facts let a required row
// be missing from this Result (OV-63). Fills, strokes, renames, and deletes
// address a row of Python's by Python's key; an addition and its styles address
// the editor row by its caption; the order addresses the keys the rows show.
// Returns whether a Legend row was added, renamed, hidden, or moved.
const applyLegendOperations = (index, operations, { displayed = false, mayBeAbsent = /** @type {((caption: string) => boolean) | undefined} */ (undefined) } = {}) => {
  const requireRow = (operation) => requireLegendEntries(index, operation.caption, operation, mayBeAbsent);
  let changed = false;
  // Python never draws a row the Legend editor added, so the row is added first and
  // a fill or stroke on it then finds it like a generated row (OV-86). An added row
  // copies Python's first row as Python drew it, before this pass styles that row.
  // A rule commit's row (`ifAbsent`) is added only where the Result lacks it,
  // before the row it replaces, and owned by the rules; a row of that caption
  // the Result has stays as it is, shown again if a commit retired it.
  // An editor row goes before the rows of the later additions, as Generate
  // appends them in order, so a renamed editor row keeps its place among them.
  operations.legendAdds.forEach(({ caption, color, xPos, yPos, ifAbsent = false, before = '' }, addIndex) => {
    const { entries, groups } = index.legends();
    const existingEntries = entries.get(caption) || [];
    if (ifAbsent && existingEntries.length > 0) {
      existingEntries.forEach((entry) => {
        if (!entry.hasAttribute(resultBaseAttribute('display'))) changed = removeAttributeIfPresent(entry, 'display') || changed;
      });
      return;
    }
    if (existingEntries.length > 0) {
      existingEntries.forEach((entry) => {
        const swatch = legendSwatch(entry);
        if (!swatch) throw new Error('Current SVG has no Legend swatch template.');
        // The editor row's own color is its drawn fill; Python drew none.
        setAttributeIfDifferent(swatch, 'fill', color);
        removeAttributeIfPresent(swatch, resultBaseAttribute('fill'));
        entry.setAttribute('data-legend-owner', 'direct-editor');
        changed = moveLegendEntryToAnchor(entry, xPos, yPos) || changed;
      });
      return;
    }
    if (groups.length === 0) {
      throw new Error('Current SVG cannot admit the requested Legend addition.');
    }
    entries.set(caption, groups.map((group) => {
      const rows = Array.from(group.querySelectorAll('g[data-legend-key]'));
      const template = rows.find((entry) => !isEditorRow(entry));
      const row = /** @type {Element | null} */ (template?.cloneNode?.(true) || null);
      if (!row) throw new Error('Current SVG has no Legend entry template.');
      // The copy is of Python's row as drawn, not of that row's edits (OV-121);
      // Python draws no hidden row.
      restorePaintBases(row, RESULT_PAINT_DOMAINS, new Map());
      row.removeAttribute('display');
      updateLegendCaption(row, caption);
      const swatch = legendSwatch(row);
      if (!swatch) throw new Error('Current SVG has no Legend swatch template.');
      swatch.setAttribute('fill', color);
      row.setAttribute('data-legend-owner', ifAbsent ? SPECIFIC_COLOR_FILE_OWNER : 'direct-editor');
      moveLegendEntryToAnchor(row, xPos, yPos);
      const later = new Set(operations.legendAdds.slice(addIndex + 1).map((add) => add.caption));
      const next = rows.find((entry) => (before
        ? text(entry.getAttribute('data-legend-key')) === before
        : !ifAbsent && isEditorRow(entry) && later.has(text(entry.getAttribute('data-legend-key')))));
      if (next?.parentElement) next.parentElement.insertBefore(row, next);
      else group.appendChild(row);
      return row;
    }));
    changed = true;
  });
  operations.legendFills.forEach((operation) => {
    const { color } = operation;
    requireRow(operation).forEach((entry) => {
      const swatch = legendSwatch(entry);
      if (!swatch) throw new Error('Sanitized SVG content is missing a Legend swatch.');
      setPaintAttribute(index, swatch, 'fill', color);
    });
  });
  operations.legendStrokes.forEach((operation) => {
    const { strokeColor, strokeWidth, renderedIds } = operation;
    (Array.isArray(renderedIds) ? renderedIds : []).forEach((renderedId) => {
      requireFeatureElements(index, renderedId).forEach((element) => {
        if (strokeColor) setPaintAttribute(index, element, 'stroke', strokeColor);
        if (strokeWidth !== null) setPaintAttribute(index, element, 'stroke-width', strokeWidth);
      });
    });
    requireRow(operation).forEach((entry) => {
      const swatch = legendSwatch(entry);
      if (!swatch) throw new Error('Sanitized SVG content is missing a Legend swatch.');
      if (strokeColor) setPaintAttribute(index, swatch, 'stroke', strokeColor);
      if (strokeWidth !== null) setPaintAttribute(index, swatch, 'stroke-width', strokeWidth);
    });
  });
  // A rename keeps the row in place and Python's key as a record; the Legend
  // layout then places every row in the Legend's order, as Python does (OV-156).
  operations.legendRenames.forEach(({ from, to, allowMissing }) => {
    requireLegendEntries(index, from, { allowMissing }, mayBeAbsent).forEach((entry) => {
      if (!setPaintAttribute(index, entry, 'data-legend-key', to)) return;
      updateLegendCaption(entry, to);
      changed = true;
    });
  });
  // A delete hides the row and keeps it with Python's key, so a reconcile
  // without the delete shows it in its place.
  // A rule commit's retired row (`retire`) stays hidden, with no record, until
  // Python draws the Legend again.
  operations.legendDeletes.forEach(({ caption, allowMissing, retire = false }) => {
    requireLegendEntries(index, caption, { allowMissing }).forEach((entry) => {
      changed = (retire ? setAttributeIfDifferent(entry, 'display', 'none') : setPaintAttribute(index, entry, 'display', 'none'))
        || changed;
    });
  });
  // The edited Legend order is replayed last, over the renderer's slots (D-08).
  // A displayed batch Result that already follows it keeps its order, so the
  // entries only that Result draws keep their places (B18).
  operations.legendOrder.forEach(({ captions }) => {
    index.legends().groups.forEach((group) => { changed = orderLegendGroup(group, captions, displayed) || changed; });
  });
  return changed;
};

// The executor's one ordering of a feature Legend group (R3). The first change
// of a group's order records Python's order of its rows (their Python keys).
/** @param {Element} group @param {readonly string[]} captions @param {boolean} keepFollowed */
const orderLegendGroup = (group, captions, keepFollowed) => {
  const pythonOrder = JSON.stringify(Array.from(group.querySelectorAll('g[data-legend-key]'))
    .filter((entry) => !isEditorRow(entry)).map(pythonLegendKey));
  const moved = orderLegendEntries(group, captions, { keepFollowed });
  if (moved && !group.hasAttribute(LEGEND_ORDER_RECORD)) group.setAttribute(LEGEND_ORDER_RECORD, pythonOrder);
  return moved;
};

// The Legend structure a reconcile shows returns to Python's where no
// operation keeps it: editor rows no addition names are removed, and a group
// whose order no order operation sets takes Python's recorded order again.
// Returns whether a row was removed or moved.
/**
 * @param {{ legends: () => { groups: Element[] } }} index
 * @param {Record<string, any>} operations
 * @param {readonly string[]} domains
 */
const restoreLegendStructure = (index, operations, domains) => {
  let changed = false;
  const { groups } = index.legends();
  if (domains.includes('legendAdds')) {
    const added = new Set(operations.legendAdds.map(({ caption }) => caption));
    groups.forEach((group) => Array.from(group.querySelectorAll('g[data-legend-key]')).forEach((entry) => {
      if (!isEditorRow(entry) || added.has(text(entry.getAttribute('data-legend-key')))) return;
      entry.remove();
      changed = true;
    }));
  }
  if (domains.includes('legendOrder') && operations.legendOrder.length === 0) {
    groups.forEach((group) => {
      const recorded = group.getAttribute(LEGEND_ORDER_RECORD);
      if (recorded === null) return;
      /** @type {unknown} */
      let order = null;
      try { order = JSON.parse(recorded); } catch { order = null; }
      const shownKeys = new Map(Array.from(group.querySelectorAll('g[data-legend-key]'))
        .filter((entry) => !isEditorRow(entry))
        .map((entry) => [pythonLegendKey(entry), text(entry.getAttribute('data-legend-key'))]));
      if (Array.isArray(order)) {
        changed = orderLegendGroup(group, order.map((key) => shownKeys.get(text(key)) || ''), false) || changed;
      }
      group.removeAttribute(LEGEND_ORDER_RECORD);
    });
  }
  return changed;
};

/**
 * Reconcile one Result's mounted SVG with its compiled editor operations,
 * using the executor that Generate admission uses (D-07, PD-OI-062). Every
 * operation given is applied; then each attribute the executor changed
 * earlier in one of `domains` that no operation set returns to the value
 * Python drew, and the Legend structure of `domains` returns to Python's where
 * no operation keeps it: a row's key and text (`legendRenames`), a hidden row
 * (`legendDeletes`), the editor rows (`legendAdds`), and the order
 * (`legendOrder`). A second call changes nothing.
 * The preview binder owns label DOM identity, so label operations stay with
 * it. Legend operations are diagram-wide and a batch Result shows only its own
 * categories and features, so an absent caption or feature is skipped.
 * Returns whether the Legend's rows changed, so the caller lays it out once.
 * @param {Element} svg
 * @param {Record<string, any>} operations
 * @param {{ resultIndex?: number, domains?: readonly string[] }} [options]
 * @returns {boolean}
 */
export const reconcileMountedResult = (svg, operations, { resultIndex = 0, domains = RESULT_PAINT_DOMAINS } = {}) => {
  const index = createLazyMutationIndex(svg, { phase: 'result-selection', resultIndex });
  const present = ({ renderedId }) => (index.features().get(renderedId) || []).length > 0;
  applyFeatureOperations(index, {
    featureFills: operations.featureFills.filter(present),
    featureStrokes: operations.featureStrokes.filter(present),
    featureVisibility: operations.featureVisibility.filter(present)
  });
  let legendChanged = false;
  if (index.legends().groups.length > 0) {
    // Every row operation addresses Python's key (the index keeps it for a
    // renamed row), or an editor row's own key, and may miss its row here.
    const allowMissing = (operation) => ({ ...operation, allowMissing: true });
    legendChanged = applyLegendOperations(index, {
      legendFills: operations.legendFills.map(allowMissing),
      legendStrokes: operations.legendStrokes.map((operation) => ({
        ...allowMissing(operation),
        renderedIds: (operation.renderedIds || []).filter((renderedId) => present({ renderedId }))
      })),
      legendRenames: operations.legendRenames.map(allowMissing),
      legendDeletes: operations.legendDeletes.map(allowMissing),
      legendAdds: operations.legendAdds,
      legendOrder: operations.legendOrder
    }, { displayed: true });
    legendChanged = restoreLegendStructure(index, operations, domains) || legendChanged;
  }
  return restorePaintBases(svg, domains, index.painted) || legendChanged;
};

/**
 * The stroke and color edits a saved Result shows, as Load reads them for its
 * mode, and the values the Session recorded as Python's.
 * @typedef {object} SavedResultEdits
 * @property {Record<string, any>} featureColorOverrides
 * @property {Record<string, any>} featureStrokeOverrides
 * @property {Record<string, any>[]} legendEntries
 * @property {Record<string, string>} legendColorOverrides
 * @property {Record<string, any>} legendStrokeOverrides
 * @property {Record<string, string>} originalLegendColors
 * @property {'circular' | 'linear'} mode
 * @property {{ color: string | null, width: number | null } | null} blockStroke Python's block stroke of
 *   the set's Result as a Session 46 kept it (`originalSvgStroke`); null when the Session's is not known to be.
 */

/** @param {string} name @param {unknown} left @param {unknown} right */
const samePaint = (name, left, right) => {
  const a = text(left).toLowerCase();
  const b = text(right).toLowerCase();
  if (name !== 'stroke-width' || !a || !b) return a === b;
  return Number(a) === Number(b);
};

// Python strokes the feature parts of one kind alike within one size class
// (a Circular record; a whole Linear canvas): a block and its outline with
// the block stroke, a block drawn under an outline with none, a connector
// with the line stroke. A Legend swatch of a feature row has the block stroke
// of the first record.
/** @param {Element} element @param {readonly Element[]} parts @param {boolean} perRecord */
const strokeKind = (element, parts, perRecord) => {
  const record = perRecord ? text(element.getAttribute('data-gbdraw-record-index')) : '';
  if (getFeaturePart(element) === FEATURE_PART_CONNECTOR) return `${record}|connector`;
  const id = text(element.getAttribute('id'));
  const outlined = Boolean(id) && parts.some((part) => part.getAttribute('id') === `${id}__outline`);
  return `${record}|${outlined ? 'outlined' : 'block'}`;
};

// Session 46 (0.14.0) and older current Sessions saved a Result with its
// stroke and color edits drawn in but without the records of Python's values.
// Load records them once, on each element that shows its edit, so Reset and
// Undo return it to Python's value: fills from the catalog, swatch fills
// from the Legend's original colors, and strokes from the parts of the same
// kind that no saved stroke edit reached (a stroke an edit or the Session
// kept may be an edited one, or another Result's). A Session 46 with one
// Result kept that Result's block stroke. A kind no part of which escaped the
// edits, or whose parts disagree, is not known: its parts get no record and
// show the saved stroke until Generate.
/**
 * @param {Element} svg
 * @param {{ resultIndex: number, catalogAdmission: FeatureCatalogAdmission, edits: SavedResultEdits,
 *   blockStroke: SavedResultEdits['blockStroke'] }} saved
 */
const recordSavedEditBases = (svg, { resultIndex, catalogAdmission, edits, blockStroke }) => {
  const index = createLazyMutationIndex(svg, { phase: 'session-load', resultIndex });
  /** @param {Element} element @param {string} name @param {unknown} edited @param {unknown} original */
  const record = (element, name, edited, original) => {
    const current = element.getAttribute(name);
    if (edited === null || edited === undefined || edited === '') return;
    // A catalog row without Python's fill cannot say what Python drew: no record.
    if (name === 'fill' && !text(original)) return;
    // A value read from the drawing returns as written; one the Session kept
    // as a number is as good as any spelling of it.
    const drawn = typeof original === 'number' ? samePaint(name, current, original) : text(current) === text(original);
    if (drawn || !samePaint(name, current, edited)) return;
    if (!element.hasAttribute(resultBaseAttribute(name))) element.setAttribute(resultBaseAttribute(name), text(original));
  };
  /** @param {string} key */
  const renderedIdsOf = (key) => (catalogAdmission.renderedTargetsByOverrideKey.get(key) || [])
    .filter((target) => target.resultIndex === resultIndex).map((target) => target.renderedId);
  const elementsOf = (/** @type {string} */ renderedId) => index.features().get(renderedId) || [];

  /** @type {Map<string, string>} */
  const editedFills = new Map();
  /** @type {Map<string, string[]>} */
  const namedIdsByCaption = new Map();
  /** @type {Map<string, Record<string, any>>} */
  const renderedFeatures = catalogAdmission.renderedFeaturesByResult?.[resultIndex] || new Map();
  Object.entries(edits.featureColorOverrides).forEach(([key, edit]) => {
    const color = text(edit && typeof edit === 'object' ? edit.color : edit);
    const caption = text(edit?.caption);
    renderedIdsOf(key).forEach((renderedId) => {
      if (color) editedFills.set(renderedId, color);
      if (caption) namedIdsByCaption.set(caption, [...(namedIdsByCaption.get(caption) || []), renderedId]);
      filterFeatureFillTargets(elementsOf(renderedId)).forEach((element) => (
        record(element, 'fill', color, renderedFeatures.get(renderedId)?.fill_color)
      ));
    });
  });
  /** @type {Array<[string, Record<string, any>]>} The rendered features each stroke edit reaches. */
  const stroked = [];
  Object.entries(edits.featureStrokeOverrides).forEach(([key, edit]) => {
    if (setsFeatureStroke(edit)) renderedIdsOf(key).forEach((renderedId) => stroked.push([renderedId, edit]));
  });
  const ownStrokeIds = stroked.map(([renderedId]) => renderedId);
  const drawnFills = [...renderedFeatures].map(([renderedId, feature]) => (
    /** @type {[string, string]} */ ([renderedId, editedFills.get(renderedId) ?? text(feature?.fill_color)])
  ));
  /** @type {Array<[string, Record<string, any>]>} */
  const strokedRows = Object.entries(edits.legendStrokeOverrides).filter(([, edit]) => setsFeatureStroke(edit));
  strokedRows.forEach(([caption, edit]) => {
    const entry = edits.legendEntries.find((row) => text(row?.caption) === caption);
    legendRowFeatureIds(entry, { drawnFills, namedIds: namedIdsByCaption.get(caption) || [], ownStrokeIds })
      .forEach((renderedId) => stroked.push([renderedId, edit]));
  });

  // Python's stroke of each kind: the one every unedited part of it shows.
  const perRecord = edits.mode === 'circular';
  const reached = new Set(stroked.map(([renderedId]) => renderedId));
  /** @type {Map<string, { color: unknown, width: unknown } | null>} null: the parts disagree. */
  const drawnStrokes = new Map();
  index.features().forEach((parts, renderedId) => {
    if (reached.has(renderedId)) return;
    parts.filter((element) => !isAutoFeatureUnderlay(element)).forEach((element) => {
      const kind = strokeKind(element, parts, perRecord);
      const drawn = { color: element.getAttribute('stroke'), width: element.getAttribute('stroke-width') };
      const known = drawnStrokes.get(kind);
      if (known === undefined) drawnStrokes.set(kind, drawn);
      else if (known && !(samePaint('stroke', known.color, drawn.color) && samePaint('stroke-width', known.width, drawn.width))) {
        drawnStrokes.set(kind, null);
      }
    });
  });
  const firstPart = [...index.features().values()].flat().find((element) => !isAutoFeatureUnderlay(element));
  const firstBlocks = `${perRecord ? text(firstPart?.getAttribute('data-gbdraw-record-index')) : ''}|block`;
  if (blockStroke && !drawnStrokes.has(firstBlocks)) drawnStrokes.set(firstBlocks, blockStroke);
  /** @param {Element} element @param {{ color: unknown, width: unknown } | null | undefined} drawn @param {Record<string, any>} edit */
  const recordStroke = (element, drawn, edit) => {
    if (!drawn) return;
    record(element, 'stroke', edit.strokeColor, drawn.color);
    record(element, 'stroke-width', edit.strokeWidth, drawn.width);
  };
  stroked.forEach(([renderedId, edit]) => {
    const parts = elementsOf(renderedId);
    parts.forEach((element) => recordStroke(
      element,
      isAutoFeatureUnderlay(element) ? AUTO_FEATURE_UNDERLAY_STROKE : drawnStrokes.get(strokeKind(element, parts, perRecord)),
      edit
    ));
  });

  const rows = index.legends().entries;
  strokedRows.forEach(([caption, edit]) => {
    (rows.get(caption) || []).forEach((row) => {
      const swatch = legendSwatch(row);
      if (swatch) recordStroke(swatch, drawnStrokes.get(firstBlocks), edit);
    });
  });
  Object.entries(edits.legendColorOverrides).forEach(([caption, color]) => {
    const entry = edits.legendEntries.find((row) => text(row?.caption) === caption);
    const original = edits.originalLegendColors[text(entry?.originalCaption) || caption];
    if (original === undefined) return;
    (rows.get(caption) || []).forEach((row) => {
      const swatch = legendSwatch(row);
      if (swatch) record(swatch, 'fill', color, original);
    });
  });
};

/**
 * A legacy normalization of a Result saved without composition metadata (a
 * Session 40 Result of main 8228ffab or 7aad9e3e, OV-273): `applies` tells it
 * from the Result's text.
 * @typedef {{ transform: SvgResultTransform, applies: (content: unknown) => boolean }} LegacyResultNormalization
 */

/**
 * The plan a current Session's Results are admitted with at Load: none, or,
 * for each Result the legacy normalization applies to, that normalization,
 * and for each Result saved without records of Python's paint, the records of
 * the edits it shows (`recordSavedEditBases`).
 * @param {readonly Record<string, any>[]} results
 * @param {FeatureCatalogAdmission} catalogAdmission
 * @param {SavedResultEdits | null} edits
 * @param {LegacyResultNormalization | null} [legacy]
 * @returns {SvgMutationPlan}
 */
export const createSavedResultPlan = (results, catalogAdmission, edits, legacy = null) => {
  const editsShown = edits && [
    edits.featureColorOverrides, edits.featureStrokeOverrides, edits.legendColorOverrides, edits.legendStrokeOverrides
  ].some((overrides) => Object.keys(overrides || {}).length > 0);
  const needsRecords = (/** @type {Record<string, any>} */ result) => (
    Boolean(editsShown) && String(result?.content || '').indexOf('data-gbdraw-base-') < 0
  );
  const normalizes = (/** @type {Record<string, any>} */ result) => Boolean(legacy?.applies(result?.content));
  const legacyNormalizationCount = results.filter(normalizes).length;
  if (legacyNormalizationCount === 0 && !results.some(needsRecords)) return createEmptySvgMutationPlan(results.length);
  const savedEdits = /** @type {SavedResultEdits} */ (edits);
  // The Session's block stroke is that of one Result; a batch's Results may
  // differ in size class.
  const kept = edits?.blockStroke;
  const blockStroke = results.length === 1 && text(kept?.color) && kept?.width !== null && kept?.width !== undefined
    ? kept : null;
  return Object.freeze({
    kind: 'MUTATING',
    legacyNormalizationCount,
    operationsByResult: Object.freeze(results.map((result, resultIndex) => {
      // The records find a row by the key the normalized Legend shows, so the
      // normalization comes first.
      const transforms = [
        ...(legacy && normalizes(result) ? [legacy.transform] : []),
        ...(needsRecords(result)
          ? [(/** @type {Element} */ svg) => recordSavedEditBases(svg, { resultIndex, catalogAdmission, edits: savedEdits, blockStroke })]
          : [])
      ];
      return transforms.length > 0
        ? Object.freeze({ ...freezeEmptyOperations(), callerTransforms: Object.freeze(transforms) })
        : freezeEmptyOperations();
    }))
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
    sourceClass,
    legendRows
  }
) => {
  const phase = sourceClass;
  const sanitized = sanitizeSvgResultContent(result, sanitizer, { phase, resultIndex });
  if (!hasSanitizedSvgEnvelope(sanitized)) {
    throw new Error('The diagram engine returned malformed SVG content.');
  }
  // Label DOM identity is reconstructed by the one mounted-preview Label binder.
  // Keep its intent in the plan, but do not require transient mounted hooks on a
  // detached Worker SVG before PreviewRuntime has installed the candidate.
  if (!hasDetachedOperations(operations)) {
    return commitCatalogBackedResult({ ...result, content: sanitized }, metadata);
  }

  const svg = parseSanitizedSvg(sanitized, parser);
  recordStructuralMetric('applicationSvgParseCount', 1, { phase, resultIndex });
  const index = createLazyMutationIndex(svg, { phase, resultIndex });
  applyFeatureOperations(index, operations);
  applyLegendOperations(index, operations, {
    mayBeAbsent: (caption) => legendRowMayBeAbsent(legendRows, resultIndex, caption)
  });
  operations.callerTransforms.forEach((transform) => transform(svg, { result, resultIndex }));
  const content = serializeAdmittedSvg(svg, { phase, resultIndex });
  return commitCatalogBackedResult({ ...result, content }, metadata);
};

/**
 * @param {SvgResultSource} source
 * @param {SvgMutationPlan | undefined} mutationPlan `requireCurrentMutationPlan` throws when it is missing.
 * @param {SvgAdmissionRuntime} [options]
 */
const admitCatalogBackedResults = (
  source,
  mutationPlan,
  {
    sanitizer = globalThis.DOMPurify || globalThis.window?.DOMPurify,
    parser = globalThis.DOMParser || globalThis.window?.DOMParser,
    selectedFeatureTypes = null
  } = {}
) => {
  const { results, catalogAdmission, sourceClass, legendRows } = source;
  requireAlignedCatalogAdmission(catalogAdmission, results);
  const plan = requireCurrentMutationPlan(mutationPlan, results.length);
  const featureTypes = resultFeatureTypes(selectedFeatureTypes);
  recordSessionLifecycleEvent('svg.admission-started', {
    phase: sourceClass,
    resultCount: results.length,
    mutationKind: plan.kind
  });
  recordStructuralMetric('currentLegacyNormalizationCount', plan.legacyNormalizationCount || 0, { phase: sourceClass });
  recordStructuralMetric('legacyOverrideMigrationCount', 0, { phase: sourceClass });
  recordStructuralMetric('manualRuleFeatureMatchCount', 0, { phase: sourceClass });
  const admitted = results.map((result, resultIndex) => admitCurrentResult(
    result,
    metadataFromCatalogAdmission(catalogAdmission, resultIndex, sourceClass, featureTypes),
    currentResultOperations(plan, resultIndex),
    { sanitizer, parser, resultIndex, sourceClass, legendRows }
  ));
  recordSessionLifecycleEvent('svg.admission-completed', {
    phase: sourceClass,
    resultCount: admitted.length,
    mutationKind: plan.kind
  });
  recordSessionLifecycleEvent('artifact.candidate-completed', {
    phase: sourceClass,
    resultCount: admitted.length,
    // `requireAlignedCatalogAdmission` above throws for a null admission.
    catalogFootprint: /** @type {FeatureCatalogAdmission} */ (catalogAdmission).scalarMetrics
  });
  return admitted;
};

/**
 * Admit only a freshly decoded current Worker response. The runtime token minted
 * by diagram-generation.js and exact catalog object alignment are both required.
 */
/**
 * @param {any} generationResponse A decoded Worker response.
 * @param {Partial<CurrentGeneratedAdmissionOptions>} [options] Admission throws when the catalog admission or plan is missing.
 */
export const admitCurrentGeneratedResults = (
  generationResponse,
  {
    catalogAdmission,
    mutationPlan,
    sanitizer,
    parser,
    selectedFeatureTypes
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
    requireAlignedCatalogAdmission(catalogAdmission, generationResponse.results),
    normalizeLegendRowFacts(generationResponse.metadata?.legendRows, generationResponse.results)
  );
  return admitCatalogBackedResults(source, mutationPlan, { sanitizer, parser, selectedFeatureTypes });
};

/**
 * @param {SvgResultSource} source
 * @param {Partial<CurrentSessionAdmissionOptions>} [options] Admission throws when the plan is missing.
 */
export const admitCurrentSessionResults = (
  source,
  {
    mutationPlan,
    sanitizer,
    parser,
    selectedFeatureTypes
  } = {}
) => admitCatalogBackedResults(
  requireSource(source, SVG_RESULT_SOURCE_CLASSES.CURRENT_SESSION),
  mutationPlan,
  { sanitizer, parser, selectedFeatureTypes }
);

const ingestSvgResult = (
  result,
  {
    sanitizer,
    parser,
    transformSvg,
    sourceClass,
    resultIndex,
    selectedFeatureTypes
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
    Object.freeze({ renderedFeatureIdentities, sourceClass, selectedFeatureTypes })
  );
};

/**
 * Admit historical/unclassified persisted SVG through compatibility repair.
 * Current Worker and current-session sources cannot enter this function.
 */
/**
 * @param {SvgResultSource} source
 * @param {LegacyImportAdmissionOptions} [options]
 */
export const admitLegacyImportedResults = (
  source,
  {
    sanitizer = globalThis.DOMPurify || globalThis.window?.DOMPurify,
    parser = globalThis.DOMParser || globalThis.window?.DOMParser,
    transformSvg = null,
    selectedFeatureTypes = null
  } = {}
) => {
  const { results, sourceClass } = requireSource(
    source,
    SVG_RESULT_SOURCE_CLASSES.LEGACY_IMPORT
  );
  const featureTypes = resultFeatureTypes(selectedFeatureTypes);
  recordSessionLifecycleEvent('svg.admission-started', {
    phase: sourceClass,
    resultCount: results.length
  });
  const admitted = results.map((result, resultIndex) => ingestSvgResult(result, {
    sanitizer,
    parser,
    transformSvg,
    sourceClass,
    resultIndex,
    selectedFeatureTypes: featureTypes
  }));
  recordSessionLifecycleEvent('svg.admission-completed', {
    phase: sourceClass,
    resultCount: admitted.length
  });
  return admitted;
};
