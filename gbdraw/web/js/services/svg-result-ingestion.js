// @ts-check
/** @import { FeatureCatalogAdmission } from './feature-catalog.js' */
/** @import { LegendRowReach, PythonLegendKey, RenderedFeatureId } from './legend-svg.js' */
/** @import { RenderedIdentityCollection } from './session-feature-metadata.js' */
/** @import { SpecificColorRule } from './specific-color-rules.js' */
/** @import { SvgSanitizer } from './svg-sanitization.js' */
/** @import { DecodedSvgDocument, HeldSegment } from './current-worker-result-source.js' */
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
  legendRowName,
  legendRowShown,
  legendRowSwatch,
  LEGEND_MERGED_ROW,
  namedLegends,
  moveLegendEntryToAnchor,
  orderLegendEntries,
  pythonLegendKey,
  pythonLegendRow,
  pythonLegendRows,
  pythonOrderedLegendRows,
  recordedGroupOrder,
  setsFeatureStroke,
  SPECIFIC_COLOR_FILE_OWNER
} from './legend-svg.js';
import { legendEntryKey } from './drawing/paint.js';
import { isCurrentWorkerGenerationResponse, SEGMENTS } from './current-worker-result-source.js';
import { diagnosticError } from '../utils/error-normalization.js';
import { drawnLegendNames } from './legend-row-names.js';
import { ruleLegendCaptions } from './specific-color-rules.js';
import { sanitizeSvgContent, unwrapSvgFragment, wrapSvgFragment } from './svg-sanitization.js';
import { LEGEND_ORDER_RECORD, RESULT_BASE_SELECTOR, pythonDrawnAttribute, resultBaseAttribute } from './result-paint-bases.js';
import { serializeCleanSvg } from './svg-serialization.js';
import { ownEntry } from '../utils/own-entry.js';
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

/** @param {unknown} value */
const text = (value) => String(value ?? '').trim();

/**
 * The root of a parsed Result SVG: `parseSanitizedSvg` admits only a document
 * whose root element is an `<svg>`.
 * @typedef {SVGSVGElement} SvgResultRoot
 */

/**
 * One Result of a Worker response or a Session: its `name` and its SVG text
 * (`content`), or for a decoded Generate reply its `document`; the other
 * members pass through unread. A committed copy holds its runtime state under
 * `COMMITTED_SVG_RESULT`, and the segments of its document under `SEGMENTS`.
 * @typedef {Record<string, any> & {
 *   [COMMITTED_SVG_RESULT]?: CommittedSvgState,
 *   [SEGMENTS]?: ReadonlyMap<string, Readonly<HeldSegment>>
 * }} SvgResultEntry
 */

/**
 * The preview owner's reader of the Results (`currentResults` of
 * app/preview-runtime.js): the displayed Result's edits are written into its
 * text before the Results are returned (R1, R13).
 * @typedef {() => SvgResultEntry[]} CurrentResultsPort
 */

/**
 * What a Result's metadata holds once it is committed (`metadataFromCatalogAdmission`,
 * `ingestSvgResult`). A legacy import has no `legendRowFacts`.
 * @typedef {object} CommittedSvgResultMetadata
 * @property {RenderedIdentityCollection} renderedFeatureIdentities
 * @property {string} sourceClass
 * @property {readonly string[] | null | undefined} selectedFeatureTypes
 * @property {SvgLegendRowFacts | null} [legendRowFacts]
 */

/**
 * The runtime state of a committed Result, never persisted.
 * @typedef {object} CommittedSvgState
 * @property {number} identity
 * @property {boolean} mounted
 * @property {string} mountedContent The markup the mounted preview was built from; never read as the Result's text.
 * @property {CommittedSvgResultMetadata} metadata
 */

/**
 * Where an admission is: the phase it reports and the Result it handles.
 * @typedef {{ phase: string, resultIndex: number }} SvgAdmissionSite
 */

/**
 * A caller's rewrite of one parsed Result SVG. The return value is ignored.
 * A transform that returns false changed nothing; any other value counts as a change (OV-345).
 * @typedef {(svg: SvgResultRoot, context: { result: SvgResultEntry, resultIndex: number }) => unknown} SvgResultTransform
 */

/**
 * A stroke on a Legend row: its swatch, and the features the row reaches in
 * the Result, which the executor reads from Python's row there (OV-288).
 * @typedef {object} LegendStrokeOperation
 * @property {PythonLegendKey} caption The row: Python's key, or an editor row's own key.
 * @property {string} strokeColor Empty: the operation sets no color.
 * @property {number | null} strokeWidth Null: the operation sets no width.
 * @property {boolean} [allowMissing]
 * @property {LegendRowReach} reach
 */

/**
 * One Result's compiled editor operations. The planner is app/candidate-render.js;
 * this module declares the shape it applies.
 * A live feature fill of `color: null` keeps the fill the Result shows.
 * @typedef {Record<
 *   'featureFills' | 'featureStrokes' | 'featureVisibility' | 'labelText' | 'labelVisibility'
 *   | 'legendFills' | 'legendRenames' | 'legendDeletes' | 'legendAdds' | 'legendOrder',
 *   readonly Record<string, any>[]
 * > & { legendStrokes: readonly LegendStrokeOperation[], callerTransforms: readonly SvgResultTransform[] }} SvgMutationOperations
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
 * @property {SvgResultEntry[]} results
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
 * @property {SvgSanitizer} [sanitizer] DOMPurify, or an object with `sanitize(svg, options)`.
 * @property {typeof DOMParser} [parser] A DOMParser constructor.
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

/**
 * @param {SvgResultEntry | null | undefined} result
 * @returns {CommittedSvgState | null}
 */
const committedState = (result) => (
  result && typeof result === 'object' && result[COMMITTED_SVG_RESULT]
    ? result[COMMITTED_SVG_RESULT]
    : null
);

/**
 * @param {unknown} results
 * @returns {SvgResultEntry[]}
 */
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

/**
 * @param {FeatureCatalogAdmission | null | undefined} admission
 * @param {unknown} results
 * @returns {FeatureCatalogAdmission}
 */
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
 * @param {SvgResultEntry[]} results
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
/**
 * @param {unknown} legendRows
 * @param {SvgResultEntry[]} results
 * @returns {readonly SvgLegendRowFacts[] | null}
 */
const normalizeLegendRowFacts = (legendRows, results) => {
  if (legendRows === undefined) return null;
  /** @param {unknown} list */
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
/**
 * @param {readonly SvgLegendRowFacts[] | null | undefined} legendRows
 * @param {number} resultIndex
 * @param {string} caption
 * @returns {boolean}
 */
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
 * @param {SvgResultEntry[]} results
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
 * @param {SvgResultEntry[]} results
 * @returns {SvgResultSource}
 */
export const createLegacyImportResultSource = (results) => (
  createRuntimeSource(
    SVG_RESULT_SOURCE_CLASSES.LEGACY_IMPORT,
    requireResultList(results)
  )
);

/**
 * @param {(SvgResultSource & { [RESULT_SOURCE]?: true }) | null | undefined} source
 * @param {string} sourceClass
 * @returns {SvgResultSource}
 */
const requireSource = (source, sourceClass) => {
  if (!source?.[RESULT_SOURCE] || source.sourceClass !== sourceClass) {
    throw new Error(`SVG Result source must be classified as ${sourceClass}.`);
  }
  return source;
};

/**
 * @param {string} content
 * @param {typeof DOMParser} [parser]
 * @returns {SvgResultRoot}
 */
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
  // `localName` is `svg` above; the parser types the root of any document as an HTMLElement.
  return /** @type {SvgResultRoot} */ (/** @type {unknown} */ (document.documentElement));
};

/**
 * @param {SvgResultEntry} result
 * @param {CommittedSvgResultMetadata} metadata
 */
const markCommitted = (result, metadata) => {
  result[COMMITTED_SVG_RESULT] = {
    identity: nextResultIdentity++,
    mounted: false,
    mountedContent: result.content,
    metadata
  };
  return result;
};

/**
 * @param {SvgResultEntry} result
 * @param {CommittedSvgResultMetadata} metadata
 */
const commitCatalogBackedResult = (result, metadata) => markCommitted(result, metadata);

/** @param {unknown} content */
const hasSanitizedSvgEnvelope = (content) => {
  const withoutDeclaration = String(content || '')
    .trim()
    .replace(/^<\?xml[^>]*>\s*/i, '');
  return /^<svg(?:\s|>)/i.test(withoutDeclaration);
};

/**
 * Every Result text and every document piece is sanitized here, and only here.
 * @param {unknown} source
 * @param {SvgSanitizer | undefined} sanitizer
 */
const sanitizeText = (source, sanitizer) => sanitizeSvgContent(source, sanitizer);

/**
 * The sanitized start of a wrapped fragment, learned once per sanitizer, so a
 * test's fake sanitizer never lends its prefix to DOMPurify.
 * @type {WeakMap<SvgSanitizer, string>}
 */
const wrapPrefixes = new WeakMap();

/** @param {SvgSanitizer | undefined} sanitizer */
const wrapPrefix = (sanitizer) => {
  if (!sanitizer) throw new Error('SVG sanitization is unavailable.');
  let prefix = wrapPrefixes.get(sanitizer);
  if (prefix === undefined) {
    prefix = withoutSuffix(sanitizeText(wrapSvgFragment(''), sanitizer), '</svg>');
    wrapPrefixes.set(sanitizer, prefix);
  }
  return prefix;
};

/**
 * @param {string} fragment Whole SVG elements.
 * @param {SvgSanitizer | undefined} sanitizer
 */
const sanitizeFragment = (fragment, sanitizer) => {
  const sanitized = unwrapSvgFragment(sanitizeText(wrapSvgFragment(fragment), sanitizer), wrapPrefix(sanitizer));
  if (sanitized === null) throw invalidResult();
  return sanitized;
};

/**
 * @param {string} text
 * @param {string} suffix
 */
const withoutSuffix = (text, suffix) => {
  if (!text.endsWith(suffix)) throw invalidResult();
  return text.slice(0, text.length - suffix.length);
};

/**
 * Sanitize a Generate document piece by piece: the root start tag, each start
 * tag and raw piece, and each children text once per revision; a revision the
 * candidate holds keeps its sanitized text. Joined in order, the pieces are the
 * sanitization of the whole document for the pieces Python writes (gbdraw/web/CLAUDE.md,
 * Runtime data flow).
 *
 * The reuse (CW-06): scope, the committed Results of the shown drawing; input,
 * the children text (its revision) and the sanitizer, and also the page:
 * DOMPurify's SANITIZE_DOM drops an `id` or `name` value that names a property
 * of `document`, and the page's `embed`, `form`, `iframe`, `img` and `object`
 * elements add such names. A held text is the sanitization against the page of
 * the Generate that first admitted it. The app page has none of those named
 * elements, so the page part does not change while the page lives (pinned by
 * the stage S browser case); lifetime, the Result that keeps the text.
 * @param {Readonly<DecodedSvgDocument>} document
 * @param {ReadonlyMap<string, Readonly<HeldSegment>>} held
 * @param {SvgSanitizer | undefined} sanitizer
 * @param {SvgAdmissionSite} site
 * @returns {{ content: string, segments: Map<string, Readonly<HeldSegment>> }}
 */
const sanitizeDocumentPieces = (document, held, sanitizer, { phase, resultIndex }) => {
  /** @type {Map<string, Readonly<HeldSegment>>} */
  const segments = new Map();
  const parts = [withoutSuffix(sanitizeText(`${document.head}</svg>`, sanitizer), '</svg>')];
  let sanitizedPieces = 0;
  let memoHits = 0;
  for (const piece of document.pieces) {
    if (piece.kind === 'raw') {
      parts.push(sanitizeFragment(piece.text, sanitizer));
      continue;
    }
    parts.push(withoutSuffix(sanitizeFragment(`${piece.start}</g>`, sanitizer), '</g>'));
    let segment = segments.get(piece.revision) ?? held.get(piece.revision);
    if (segment) {
      memoHits += 1;
    } else {
      if (typeof piece.children !== 'string') throw invalidResult();
      segment = Object.freeze({
        revision: piece.revision,
        bytes: piece.bytes,
        sanitized: sanitizeFragment(piece.children, sanitizer)
      });
      sanitizedPieces += 1;
    }
    segments.set(piece.revision, segment);
    parts.push(segment.sanitized, '</g>');
  }
  parts.push('</svg>');
  recordStructuralMetric('svgSanitizedPieceCount', sanitizedPieces, { phase, resultIndex });
  recordStructuralMetric('svgMemoHitCount', memoHits, { phase, resultIndex });
  return { content: parts.join(''), segments };
};

/**
 * Sanitize one Result: a current Worker document piece by piece, any other
 * Result (a saved Session, a legacy import, a test) as its whole `content`.
 * @param {SvgResultEntry} result
 * @param {SvgSanitizer | undefined} sanitizer
 * @param {SvgAdmissionSite} site
 * @returns {{ content: string, segments: Map<string, Readonly<HeldSegment>> | null }}
 */
const sanitizeAdmittedSvg = (result, sanitizer, site) => {
  const { phase, resultIndex } = site;
  /** @type {Readonly<DecodedSvgDocument> | undefined} */
  const document = phase === SVG_RESULT_SOURCE_CLASSES.CURRENT_WORKER ? result.document : undefined;
  const admitted = document
    ? sanitizeDocumentPieces(document, result[SEGMENTS] ?? new Map(), sanitizer, site)
    : { content: sanitizeText(result.content, sanitizer), segments: null };
  recordSessionLifecycleEvent('result-svg-characters', {
    phase,
    resultIndex,
    value: document ? admitted.content.length : String(result.content || '').length
  });
  recordStructuralMetric('svgSanitizationCount', 1, { phase, resultIndex });
  return admitted;
};

/**
 * The Result a sanitized candidate commits as: its members without the decoded
 * document, the sanitized `content`, and under `SEGMENTS` the sanitized
 * segment of each revision of its document (none for a whole-text Result).
 * @param {SvgResultEntry} result
 * @param {{ content: string, segments: Map<string, Readonly<HeldSegment>> | null }} admitted
 * @returns {SvgResultEntry}
 */
const admittedResult = (result, { content, segments }) => {
  const { document: _document, [SEGMENTS]: _held, ...rest } = result;
  return segments ? { ...rest, content, [SEGMENTS]: segments } : { ...rest, content };
};

/**
 * @param {Element} svg
 * @param {SvgAdmissionSite} site
 */
const serializeAdmittedSvg = (svg, { phase, resultIndex }) => {
  const content = serializeCleanSvg(svg);
  recordStructuralMetric('svgSerializationCount', 1, { phase, resultIndex });
  return content;
};

// The feature types of the request that drew the Results (its
// `diagramOptions.selectedFeaturesSet`), or null when it names none. The
// Features list reads them from the displayed Result's metadata (R13).
/** @param {readonly string[] | null | undefined} selectedFeatureTypes */
const resultFeatureTypes = (selectedFeatureTypes) => (
  Array.isArray(selectedFeatureTypes) ? Object.freeze(selectedFeatureTypes.map(String)) : null
);

// `legendRowFacts`: Python's Legend row facts of the Result (the rows it drew and the
// keys no row draws), kept in memory with the committed Result; null for a Result Load admitted.
/**
 * @param {FeatureCatalogAdmission} catalogAdmission @param {number} resultIndex @param {string} sourceClass
 * @param {readonly string[] | null} selectedFeatureTypes
 * @param {readonly SvgLegendRowFacts[] | null | undefined} [legendRows]
 * @returns {CommittedSvgResultMetadata}
 */
const metadataFromCatalogAdmission = (catalogAdmission, resultIndex, sourceClass, selectedFeatureTypes, legendRows = null) => {
  const renderedFeatureIdentities = catalogAdmission.renderedIdentitiesByResult[resultIndex];
  if (!renderedFeatureIdentities) {
    throw new Error('The admitted feature catalog is missing Result identity metadata.');
  }
  return Object.freeze({
    renderedFeatureIdentities, sourceClass, selectedFeatureTypes, legendRowFacts: legendRows?.[resultIndex] ?? null
  });
};

/** @param {SvgResultEntry | null | undefined} result */
export const isCommittedSvgResult = (result) => Boolean(committedState(result));

/** @param {SvgResultEntry | null | undefined} result */
export const isCommittedSvgResultMounted = (result) => Boolean(
  committedState(result)?.mounted
);

/** @param {SvgResultEntry} result */
export const getCommittedSvgContent = (result) => {
  const runtime = committedState(result);
  if (!runtime) return null;
  return runtime.mounted ? runtime.mountedContent : result.content;
};

/** @param {SvgResultEntry | null | undefined} result */
export const getCommittedSvgResultMetadata = (result) => (
  committedState(result)?.metadata || null
);

/** @param {SvgResultEntry | null | undefined} result */
export const getCommittedSvgResultRuntimeIdentity = (result) => (
  committedState(result)?.identity ?? null
);

/** @param {SvgResultEntry} result */
export const markCommittedSvgResultMounted = (result) => {
  const runtime = committedState(result);
  if (!runtime) return false;
  runtime.mounted = true;
  runtime.mountedContent = result.content;
  return true;
};

/** @param {SvgResultEntry} result */
export const markCommittedSvgResultUnmounted = (result) => {
  const runtime = committedState(result);
  if (!runtime) return false;
  runtime.mounted = false;
  runtime.mountedContent = result.content;
  return true;
};

/**
 * @param {SvgMutationPlan} plan
 * @param {number} resultIndex
 * @returns {SvgMutationOperations | null}
 */
const currentResultOperations = (plan, resultIndex) => (
  plan.operationsByResult[resultIndex] || null
);

/** @param {SvgMutationOperations | null | undefined} operations */
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

/**
 * @param {SvgMutationOperations | null | undefined} operations
 * @returns {operations is SvgMutationOperations}
 */
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

/**
 * @param {SvgMutationPlan | undefined} plan
 * @param {number} resultCount
 * @returns {SvgMutationPlan}
 */
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

/**
 * @param {Element} element
 * @param {string} name
 * @param {unknown} value
 */
const setAttributeIfDifferent = (element, name, value) => {
  const normalized = String(value);
  if (element.getAttribute(name) === normalized) return false;
  element.setAttribute(name, normalized);
  return true;
};

/**
 * @param {Element} element
 * @param {string} name
 */
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
// each attribute an operation sets, so a reconcile leaves it alone, and
// whether the pass changed the Result.
/**
 * @param {{ painted: Map<Element, Set<string>>, changed: boolean }} index
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
  index.changed = true;
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

// THE writer of the names a Result's Legend draws (D-09, PD-OI-042 scenario
// revision 2): after the Legend operations of a pass, each shown row's text is
// the naming rule's answer (legend-row-names.js) for the rows the Legend shows,
// by their names and the swatches they show. No operation writes a row's text;
// Python's own texts are this rule's answer for the rows it drew. Rows of one
// name and one swatch are one row (D-28): the later of them is hidden under
// `LEGEND_MERGED_ROW`, which holds the key of the earlier row (D-35,
// `legendRowKeys`), and shows again on the pass that finds its name or
// swatch different (a delete or a retirement of a merged row first shows it,
// `unmergeLegendRow`). Returns whether a text or a merge changed.
/** @param {Parameters<typeof namedLegends>[0]} svg */
const nameLegendRows = (svg) => {
  let changed = false;
  const displayRecord = resultBaseAttribute('display');
  namedLegends(svg).forEach((legend) => {
    const all = Array.from(legend.querySelectorAll('g[data-legend-key]'));
    const merged = new Map(all.filter((row) => row.hasAttribute(LEGEND_MERGED_ROW))
      .map((row) => [row, row.getAttribute(LEGEND_MERGED_ROW)]));
    merged.forEach((_, row) => {
      row.removeAttribute(LEGEND_MERGED_ROW);
      if (!row.hasAttribute(displayRecord)) row.removeAttribute('display');
    });
    /** @type {Map<string, string>} A name and swatch, and the key of the row that draws them. */
    const seen = new Map();
    const rows = all.filter(legendRowShown).filter((row) => {
      const identity = JSON.stringify([legendRowName(row), legendRowSwatch(row)]);
      const into = seen.get(identity);
      if (into === undefined) return Boolean(seen.set(identity, pythonLegendKey(row)));
      row.setAttribute('display', 'none');
      row.setAttribute(LEGEND_MERGED_ROW, into);
      if (merged.get(row) !== into) changed = true;
      merged.delete(row);
      return false;
    });
    if (merged.size > 0) changed = true;
    const drawn = drawnLegendNames(rows.map((row) => ({ name: legendRowName(row), swatch: legendRowSwatch(row) })));
    rows.forEach((row, at) => {
      const label = row.querySelector('text');
      if (!label || label.textContent === drawn[at]) return;
      label.textContent = drawn[at];
      changed = true;
    });
  });
  return changed;
};

// A row the naming pass merged shows again, so an operation that hides it
// hides a shown row (D-28).
/** @param {Element} row */
const unmergeLegendRow = (row) => {
  if (!row.hasAttribute(LEGEND_MERGED_ROW)) return;
  row.removeAttribute(LEGEND_MERGED_ROW);
  row.removeAttribute('display');
};

// An added row's key and name: `data-legend-name` only where the name is not
// the key, as Python writes it (a rule's row, whose name is the rule's caption).
/** @param {Element} row @param {string} key @param {string | undefined} name */
const setAddedRowIdentity = (row, key, name) => {
  row.setAttribute('data-legend-key', key);
  if (name && name !== key) row.setAttribute('data-legend-name', name);
  else row.removeAttribute('data-legend-name');
};

/**
 * Return every attribute the executor changed in `domains` and no operation
 * of this pass set to the value Python drew (a Legend row's key; its text
 * follows in `nameLegendRows`). Returns whether an attribute returned, and whether a Legend row did.
 * @param {Element} root
 * @param {readonly string[]} domains
 * @param {Map<Element, Set<string>>} painted
 * @returns {{ changed: boolean, rowsChanged: boolean }}
 */
const restorePaintBases = (root, domains, painted) => {
  const owned = { feature: new Set(), swatch: new Set(), row: new Set() };
  domains.forEach((domain) => {
    const attributes = PAINT_DOMAIN_ATTRIBUTES[/** @type {keyof typeof PAINT_DOMAIN_ATTRIBUTES} */ (domain)];
    attributes?.feature.forEach((name) => owned.feature.add(name));
    attributes?.swatch.forEach((name) => owned.swatch.add(name));
    attributes?.row.forEach((name) => owned.row.add(name));
  });
  let changed = false;
  let rowsChanged = false;
  if (owned.feature.size === 0 && owned.swatch.size === 0 && owned.row.size === 0) return { changed, rowsChanged };
  [root, ...Array.from(root.querySelectorAll(RESULT_BASE_SELECTOR))].forEach((element) => {
    const kind = paintElementKind(element);
    owned[kind].forEach((name) => {
      if (painted.get(element)?.has(name)) return;
      const base = resultBaseAttribute(name);
      const value = element.getAttribute(base);
      if (value === null) return;
      if (value === '') element.removeAttribute(name);
      else element.setAttribute(name, value);
      element.removeAttribute(base);
      changed = true;
      if (kind === 'row') rowsChanged = true;
    });
  });
  return { changed, rowsChanged };
};

/**
 * The lookups one admission or reconcile builds from the Result on first use,
 * and what the operations paint.
 * @typedef {object} LazyMutationIndex
 * @property {Map<Element, Set<string>>} painted The attributes an operation set, by element.
 * @property {boolean} changed Whether a paint attribute changed.
 * @property {() => Map<string, Element[]>} features The elements of each rendered Feature.
 * @property {() => { entries: Map<string, Element[]>, groups: Element[] }} legends The Legend rows by key, and the feature Legend groups.
 */

/**
 * @param {SvgResultRoot} svg
 * @param {SvgAdmissionSite} site
 * @returns {LazyMutationIndex}
 */
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
    changed: false,
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
      const rowsByKey = new Map();
      built.legendEntries = rowsByKey;
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
          if (!rowsByKey.has(caption)) rowsByKey.set(caption, []);
          // `get` holds the array set on the line above.
          /** @type {Element[]} */ (rowsByKey.get(caption)).push(entry);
        });
      });
      recordStructuralMetric('legendDomFullScanCount', 1, { phase, resultIndex });
      return { entries: built.legendEntries, groups: built.legendGroups };
    }
  };
};

/** @param {Element} entry */
const isEditorRow = (entry) => entry.getAttribute('data-legend-owner') === 'direct-editor';

/**
 * @param {LazyMutationIndex} index
 * @param {string} renderedId
 */
const requireFeatureElements = (index, renderedId) => {
  const elements = index.features().get(renderedId) || [];
  if (elements.length === 0) {
    throw new Error('Sanitized SVG content is missing a rendered Feature binding.');
  }
  return elements;
};

/**
 * @param {LazyMutationIndex} index
 * @param {string} caption
 * @param {{ allowMissing?: boolean }} [options]
 * @param {(caption: string) => boolean} [mayBeAbsent]
 */
const requireLegendEntries = (index, caption, { allowMissing = false } = {}, mayBeAbsent = () => false) => {
  const entries = index.legends().entries.get(caption) || [];
  if (entries.length === 0 && !allowMissing && !mayBeAbsent(caption)) throw invalidResult();
  return entries;
};

/**
 * @param {LazyMutationIndex} index
 * @param {Pick<SvgMutationOperations, 'featureFills' | 'featureStrokes' | 'featureVisibility'>} operations
 */
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
// `restoresFills`: the pass returns the feature fills no operation sets to
// Python's (a reconcile that shows them).
// Returns whether a Legend row was added, renamed, hidden, or moved.
/**
 * @param {LazyMutationIndex} index
 * @param {Pick<SvgMutationOperations, 'legendFills' | 'legendStrokes' | 'legendRenames' | 'legendDeletes' | 'legendAdds' | 'legendOrder'>} operations
 * @param {{ displayed?: boolean, restoresFills?: boolean, mayBeAbsent?: (caption: string) => boolean }} [options]
 */
const applyLegendOperations = (index, operations, {
  displayed = false, restoresFills = false, mayBeAbsent = /** @type {((caption: string) => boolean) | undefined} */ (undefined)
} = {}) => {
  /** @param {SvgMutationOperations['legendFills'][number] | LegendStrokeOperation} operation A fill, a stroke, or the like: each names its row by `caption`. */
  const requireRow = (operation) => requireLegendEntries(index, operation.caption, operation, mayBeAbsent);
  // The fill each drawn feature shows once this pass ends, read once per pass
  // for the Legend row strokes (perf-014x S5).
  /** @type {Array<[RenderedFeatureId, string]> | null} */
  let drawnFills = null;
  const featureFills = () => {
    drawnFills ??= [...index.features()].flatMap(([renderedId, elements]) => {
      const element = filterFeatureFillTargets(elements)[0];
      if (!element) return [];
      const fill = restoresFills && !index.painted.get(element)?.has('fill')
        ? pythonDrawnAttribute(element, 'fill')
        : element.getAttribute('fill');
      return [/** @type {[RenderedFeatureId, string]} */ ([renderedId, text(fill)])];
    });
    return drawnFills;
  };
  let changed = false;
  // Python never draws a row the Legend editor added, so the row is added first and
  // a fill or stroke on it then finds it like a generated row (OV-86). An added row
  // copies Python's first row as Python drew it, before this pass styles that row,
  // whatever row a user order shows first (OV-2188).
  // A rule commit's row (`ifAbsent`) is added only where the Result lacks it,
  // before the row it replaces, and owned by the rules; a row of that caption
  // the Result has stays as it is, shown again if a commit retired it.
  // An editor row goes before the rows of the later additions, as Generate
  // appends them in order, so a renamed editor row keeps its place among them.
  operations.legendAdds.forEach(({ caption, name, color, xPos, yPos, ifAbsent = false, before = '' }, addIndex) => {
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
        const recolored = setAttributeIfDifferent(swatch, 'fill', color);
        const unrecorded = removeAttributeIfPresent(swatch, resultBaseAttribute('fill'));
        setAddedRowIdentity(entry, caption, name);
        if (setAttributeIfDifferent(entry, 'data-legend-owner', 'direct-editor') || recolored || unrecorded) index.changed = true;
        changed = moveLegendEntryToAnchor(entry, xPos, yPos) || changed;
      });
      return;
    }
    if (groups.length === 0) {
      throw new Error('Current SVG cannot admit the requested Legend addition.');
    }
    entries.set(caption, groups.map((group) => {
      const rows = Array.from(group.querySelectorAll('g[data-legend-key]'));
      // The template is Python's first row, whatever row a user order shows
      // first (OV-2188): its own row when it has one, since a rule commit's row
      // an older build saved may carry a copied GC row's "no line".
      const pythonRows = pythonOrderedLegendRows(group);
      const template = pythonRows.find((pythonRow) => !pythonRow.hasAttribute('data-legend-owner')) ?? pythonRows[0];
      const row = /** @type {Element | null} */ (template?.cloneNode(true) || null);
      if (!row) throw new Error('Current SVG has no Legend entry template.');
      // The copy is of Python's row as drawn, not of that row's edits (OV-121);
      // Python draws no hidden row.
      restorePaintBases(row, RESULT_PAINT_DOMAINS, new Map());
      row.removeAttribute('display');
      setAddedRowIdentity(row, caption, name);
      const swatch = legendSwatch(row);
      if (!swatch) throw new Error('Current SVG has no Legend swatch template.');
      swatch.setAttribute('fill', color);
      row.setAttribute('data-legend-owner', ifAbsent ? SPECIFIC_COLOR_FILE_OWNER : 'direct-editor');
      moveLegendEntryToAnchor(row, xPos, yPos);
      const later = new Set(operations.legendAdds.slice(addIndex + 1).map((add) => add.caption));
      const next = rows.find((entry) => (before
        ? pythonLegendKey(entry) === before
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
  operations.legendStrokes.forEach((/** @type {LegendStrokeOperation} */ operation) => {
    const { strokeColor, strokeWidth, reach } = operation;
    const rows = requireRow(operation);
    const pythonColor = rows.length > 0 ? pythonLegendRow(rows[0]).color : null;
    legendRowFeatureIds(reach, pythonColor, featureFills()).forEach((renderedId) => {
      requireFeatureElements(index, renderedId).forEach((element) => {
        if (strokeColor) setPaintAttribute(index, element, 'stroke', strokeColor);
        if (strokeWidth !== null) setPaintAttribute(index, element, 'stroke-width', strokeWidth);
      });
    });
    rows.forEach((entry) => {
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
      changed = setPaintAttribute(index, entry, 'data-legend-key', to) || changed;
    });
  });
  // A delete hides the row and keeps it with Python's key, so a reconcile
  // without the delete shows it in its place.
  // A rule commit's retired row (`retire`) stays hidden, with no record, until
  // Python draws the Legend again.
  operations.legendDeletes.forEach(({ caption, allowMissing, retire = false }) => {
    requireLegendEntries(index, caption, { allowMissing }).forEach((entry) => {
      unmergeLegendRow(entry);
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
  const pythonOrder = JSON.stringify(pythonOrderedLegendRows(group).map(pythonLegendKey));
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
 * @param {SvgMutationOperations} operations
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
      if (!group.hasAttribute(LEGEND_ORDER_RECORD)) return;
      const order = recordedGroupOrder(group);
      if (order) changed = orderLegendGroup(group, order.map(text), false) || changed;
      // The record's removal is a change, so the stored Result loses it too.
      group.removeAttribute(LEGEND_ORDER_RECORD);
      changed = true;
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
 * Returns whether the reconcile changed the Result, so a caller commits only
 * a change (OV-345), and whether the Legend's rows changed, so the caller
 * lays it out once.
 * @param {SvgResultRoot} svg
 * @param {SvgMutationOperations} operations
 * @param {{ resultIndex?: number, domains?: readonly string[] }} [options]
 * @returns {{ changed: boolean, legendChanged: boolean }}
 */
export const reconcileMountedResult = (svg, operations, { resultIndex = 0, domains = RESULT_PAINT_DOMAINS } = {}) => {
  const index = createLazyMutationIndex(svg, { phase: 'result-selection', resultIndex });
  const present = (/** @type {Readonly<Record<string, unknown>>} */ operation) => (
    index.features().get(String(operation.renderedId)) || []
  ).length > 0;
  applyFeatureOperations(index, {
    featureFills: operations.featureFills.filter(present),
    featureStrokes: operations.featureStrokes.filter(present),
    featureVisibility: operations.featureVisibility.filter(present)
  });
  let legendChanged = false;
  if (index.legends().groups.length > 0) {
    // Every row operation addresses Python's key (the index keeps it for a
    // renamed row), or an editor row's own key, and may miss its row here.
    /**
     * @template {Readonly<Record<string, unknown>>} T
     * @param {T} operation
     */
    const allowMissing = (operation) => ({ ...operation, allowMissing: true });
    legendChanged = applyLegendOperations(index, {
      legendFills: operations.legendFills.map(allowMissing),
      legendStrokes: operations.legendStrokes.map(allowMissing),
      legendRenames: operations.legendRenames.map(allowMissing),
      legendDeletes: operations.legendDeletes.map(allowMissing),
      legendAdds: operations.legendAdds,
      legendOrder: operations.legendOrder
    }, { displayed: true, restoresFills: domains.includes('featureFills') });
    legendChanged = restoreLegendStructure(index, operations, domains) || legendChanged;
  }
  const restored = restorePaintBases(svg, domains, index.painted);
  const named = nameLegendRows(svg);
  return {
    changed: index.changed || legendChanged || restored.changed || named,
    legendChanged: restored.rowsChanged || legendChanged || named
  };
};

/**
 * A saved feature color edit: a color, or a color with the Legend caption it
 * names (older Sessions wrote the color alone).
 * @typedef {string | { color?: unknown, caption?: unknown } | null} SavedFeatureColorEdit
 */
/**
 * A saved stroke edit of a feature or a Legend row.
 * @typedef {{ strokeColor?: unknown, strokeWidth?: unknown }} SavedStrokeEdit
 */
/**
 * A saved Legend row as Load reads it: its caption, Python's key, and the
 * features older Sessions listed.
 * @typedef {{ caption?: unknown, originalCaption?: unknown, featureIds?: unknown }} SavedLegendEntry
 */
/**
 * The stroke and color edits a saved Result shows, as Load reads them for its
 * mode, and the values the Session recorded as Python's.
 * @typedef {object} SavedResultEdits
 * @property {Record<string, unknown>} featureColorOverrides Unvalidated: each a `SavedFeatureColorEdit` or anything else.
 * @property {Record<string, SavedStrokeEdit>} featureStrokeOverrides
 * @property {SavedLegendEntry[]} legendEntries
 * @property {Record<string, string>} legendColorOverrides
 * @property {Record<string, SavedStrokeEdit>} legendStrokeOverrides
 * @property {Record<string, string>} originalLegendColors
 * @property {string[]} originalLegendOrder The captions of the rows Python drew (the Legend inventory).
 * @property {readonly Partial<SpecificColorRule>[]} rules The color rules the mode's drawing loads, read as saved rules.
 * @property {ReadonlyMap<PythonLegendKey, PythonLegendKey>} [legendRekeys] The rule rows whose key that reading
 *   changes (`savedRuleRowRekeys`, F2), old key to new; the other fields are already under the new keys.
 * @property {'circular' | 'linear'} mode
 * @property {boolean} resultsRecordPythonPaint The Session's Results carry records of Python's paint
 *   (`SourceFacts`, read from the file's own version: a Session 27-44 Python converted saved none).
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

// Sessions older than 46 (main writes 44) saved a Result with its stroke and
// color edits drawn in but without the records of Python's values. Load
// records them once, on each element that shows its edit, so Reset and Undo
// return it to Python's value: fills from the catalog, swatch fills from the
// Legend's original colors, and strokes from the parts of the same kind that
// no saved stroke edit reached (a stroke an edit or the Session kept may be an
// edited one, or another Result's). A kind no part of which escaped the
// edits, or whose parts disagree, is not known: its parts get no record and
// show the saved stroke until Generate. Returns whether it recorded a value.
/**
 * @param {SvgResultRoot} svg
 * @param {{ resultIndex: number, catalogAdmission: FeatureCatalogAdmission, edits: SavedResultEdits }} saved
 * @returns {boolean}
 */
const recordSavedEditBases = (svg, { resultIndex, catalogAdmission, edits }) => {
  const index = createLazyMutationIndex(svg, { phase: 'session-load', resultIndex });
  let recorded = false;
  /** @param {Element} element @param {string} name @param {unknown} edited @param {unknown} original */
  const record = (element, name, edited, original) => {
    const current = element.getAttribute(name);
    if (edited === null || edited === undefined || edited === '') return;
    // A catalog row without Python's fill cannot say what Python drew: no record.
    if (name === 'fill' && !text(original)) return;
    // A value read from the drawing returns as written; one the Session kept
    // as a number is as good as any spelling of it.
    const drawn = typeof original === 'number' ? samePaint(name, current, original) : text(current) === text(original);
    if (drawn || !samePaint(name, current, edited) || element.hasAttribute(resultBaseAttribute(name))) return;
    element.setAttribute(resultBaseAttribute(name), text(original));
    recorded = true;
  };
  /**
   * @param {string} key
   * @returns {RenderedFeatureId[]}
   */
  const renderedIdsOf = (key) => (catalogAdmission.renderedTargetsByOverrideKey.get(key) || [])
    .filter((/** @type {{ resultIndex: number }} */ target) => target.resultIndex === resultIndex).map((/** @type {{ renderedId: RenderedFeatureId }} */ target) => target.renderedId);
  const elementsOf = (/** @type {string} */ renderedId) => index.features().get(renderedId) || [];

  /** @type {Map<string, string>} */
  const editedFills = new Map();
  /** @type {Map<string, RenderedFeatureId[]>} */
  const namedIdsByCaption = new Map();
  /** @type {Map<string, { fill_color?: unknown }>} */
  const renderedFeatures = catalogAdmission.renderedFeaturesByResult?.[resultIndex] || new Map();
  Object.entries(edits.featureColorOverrides).forEach(([key, edit]) => {
    // The cast: an object (the test), whose fields are unchecked.
    const fields = edit && typeof edit === 'object' ? /** @type {{ color?: unknown, caption?: unknown }} */ (edit) : null;
    const color = text(fields ? fields.color : edit);
    const caption = text(fields ? fields.caption : '');
    renderedIdsOf(key).forEach((renderedId) => {
      if (color) editedFills.set(renderedId, color);
      if (caption) namedIdsByCaption.set(caption, [...(namedIdsByCaption.get(caption) || []), /** @type {RenderedFeatureId} */ (renderedId)]);
      filterFeatureFillTargets(elementsOf(renderedId)).forEach((element) => (
        record(element, 'fill', color, renderedFeatures.get(renderedId)?.fill_color)
      ));
    });
  });
  /** @type {Array<[string, SavedStrokeEdit]>} The rendered features each stroke edit reaches. */
  const stroked = [];
  Object.entries(edits.featureStrokeOverrides).forEach(([key, edit]) => {
    if (setsFeatureStroke(edit)) renderedIdsOf(key).forEach((renderedId) => stroked.push([renderedId, edit]));
  });
  const ownStrokeIds = /** @type {RenderedFeatureId[]} */ (stroked.map(([renderedId]) => renderedId));
  const drawnFills = [...renderedFeatures].map(([renderedId, feature]) => (
    /** @type {[RenderedFeatureId, string]} */ ([renderedId, editedFills.get(renderedId) ?? text(feature?.fill_color)])
  ));
  const rows = index.legends().entries;
  // A swatch shows its Legend color; its record keeps Python's, which is the
  // color the row's stroke reaches by (OV-288), so the records come first.
  // Both stores are keyed by the row's key.
  Object.entries(edits.legendColorOverrides).forEach(([key, color]) => {
    // A row key is data: only an own entry is a saved color (OV-982).
    const original = ownEntry(edits.originalLegendColors, key);
    if (original === undefined) return;
    (rows.get(key) || []).forEach((row) => {
      const swatch = legendSwatch(row);
      if (swatch) record(swatch, 'fill', color, original);
    });
  });
  /** @type {Array<[string, SavedStrokeEdit]>} */
  const strokedRows = Object.entries(edits.legendStrokeOverrides).filter(([, edit]) => setsFeatureStroke(edit));
  strokedRows.forEach(([key, edit]) => {
    const entry = edits.legendEntries.find((row) => legendEntryKey(row) === key);
    const row = rows.get(key)?.[0];
    /** @type {LegendRowReach} */
    const reach = {
      listedIds: /** @type {RenderedFeatureId[]} */ ((Array.isArray(entry?.featureIds) ? entry.featureIds : []).map(text).filter(Boolean)),
      namedIds: namedIdsByCaption.get(key) || [],
      ownStrokeIds,
      draftColor: null
    };
    legendRowFeatureIds(reach, row ? pythonLegendRow(row).color : null, drawnFills)
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
  /** @param {Element} element @param {{ color: unknown, width: unknown } | null | undefined} drawn @param {SavedStrokeEdit} edit */
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

  strokedRows.forEach(([caption, edit]) => {
    (rows.get(caption) || []).forEach((row) => {
      const swatch = legendSwatch(row);
      if (swatch) recordStroke(swatch, drawnStrokes.get(firstBlocks), edit);
    });
  });
  return recorded;
};

// The renamed generated rows a Session lists, by their caption: main's writer
// (Sessions 44 to 46) renamed a row of Python's by rewriting its key and kept
// Python's caption only as the entry's `originalCaption`. Only a caption no
// other listed entry has, that is not itself a generated caption (Python may
// draw a row of it), names one such row.
/** @param {SavedResultEdits} edits @returns {Map<string, string>} */
const savedLegendRenames = (edits) => {
  const generated = new Set(edits.originalLegendOrder.map(text));
  /** @type {Map<string, number>} */
  const listed = new Map();
  edits.legendEntries.forEach((entry) => listed.set(text(entry?.caption), (listed.get(text(entry?.caption)) || 0) + 1));
  return new Map(edits.legendEntries.flatMap((entry) => {
    const caption = text(entry?.caption);
    const original = text(entry?.originalCaption);
    return caption && original && caption !== original && generated.has(original) && !generated.has(caption)
      && listed.get(caption) === 1 ? [[caption, original]] : [];
  }));
};

// A Result saved before the executor recorded Python's Legend key shows such
// a rename without the record, so a live edit addressed by Python's key misses
// the row. Load records it once, on the one row of Python's with the renamed
// caption, in a Result that has no row of Python's caption (U3a H1).
// Returns whether it recorded a key.
/** @param {SvgResultRoot} svg @param {Map<string, string>} renames @returns {boolean} */
const recordSavedLegendKeys = (svg, renames) => {
  let recorded = false;
  const keyRecord = resultBaseAttribute('data-legend-key');
  const groups = getAllFeatureLegendGroups(svg);
  const rows = groups.flatMap((group) => Array.from(group.querySelectorAll('g[data-legend-key]')));
  const drawnKeys = new Set(rows.flatMap((row) => [text(row.getAttribute('data-legend-key')), text(row.getAttribute(keyRecord))]));
  groups.forEach((group) => {
    const groupRows = Array.from(group.querySelectorAll('g[data-legend-key]'));
    renames.forEach((original, caption) => {
      if (drawnKeys.has(original)) return;
      const named = groupRows.filter((row) => text(row.getAttribute('data-legend-key')) === caption);
      const [row] = named;
      if (named.length === 1 && !row.hasAttribute('data-legend-owner') && !row.hasAttribute(keyRecord)) {
        row.setAttribute(keyRecord, original);
        recorded = true;
      }
    });
  });
  return recorded;
};

// F2: a rule row whose key the reading of the saved captions changes (a
// caption saved with its own color ID whose rule is the only rule of that
// caption) shows its old key in the saved figure. Load moves the row to its
// new key, as it moves the edits stored under it (`rekeyLegendRecords`), and
// draws the names. A rename's record of the old key moves with it. Returns
// whether it moved a row.
/** @param {SvgResultRoot} svg @param {ReadonlyMap<PythonLegendKey, PythonLegendKey>} rekeyed @returns {boolean} */
const rekeySavedRuleRows = (svg, rekeyed) => {
  const lookup = /** @type {ReadonlyMap<string, string>} */ (rekeyed);
  let moved = false;
  getAllFeatureLegendGroups(svg).flatMap((group) => Array.from(group.querySelectorAll('g[data-legend-key]'))).forEach((row) => {
    ['data-legend-key', resultBaseAttribute('data-legend-key')].forEach((attribute) => {
      const to = lookup.get(text(row.getAttribute(attribute)));
      if (to === undefined) return;
      row.setAttribute(attribute, to);
      moved = true;
    });
  });
  return moved && (nameLegendRows(svg), true);
};

// D-19: a build of first-parent main (fe6861f0 to this naming rule) drew a
// color rule's row under its key, `Complex IV [#ffe700]` or `tRNA [#ff0000]`,
// as the row's text and wrote no `data-legend-name`. Load names such a row
// once, as Python names it now: the row whose key the key mirror allocates to
// a loaded rule takes that rule's caption as its name, and the naming pass
// draws the names. The key stays, so the edits saved under it stay on the row.
// A row that has a name, or shows a rename, keeps it. Returns whether it
// named a row.
const SAVED_COLOR_ID = / \[(?:#[0-9a-f]{6}|none)\]/;
/** @param {SvgResultRoot} svg @param {SavedResultEdits} edits @returns {boolean} */
const nameSavedRuleRows = (svg, edits) => {
  const keyOf = ruleLegendCaptions({
    rules: edits.rules, pythonRows: pythonLegendRows(svg), features: [], originalLegendOrder: edits.originalLegendOrder
  });
  /** @type {Map<string, string>} */
  const names = new Map(edits.rules.filter((rule) => text(rule.cap)).map((rule) => [keyOf(rule), text(rule.cap)]));
  const keyRecord = resultBaseAttribute('data-legend-key');
  let named = false;
  namedLegends(svg).flatMap((legend) => Array.from(legend.querySelectorAll('g[data-legend-key]'))).forEach((row) => {
    const key = text(row.getAttribute('data-legend-key'));
    const name = names.get(key);
    if (!name || name === key || row.hasAttribute('data-legend-name') || row.hasAttribute(keyRecord)) return;
    row.setAttribute('data-legend-name', name);
    named = true;
  });
  return named && (nameLegendRows(svg), true);
};

/**
 * A legacy normalization of a Result saved without composition metadata (a
 * Session 40 Result of main 8228ffab or 7aad9e3e, OV-273): `applies` tells it
 * from the Result's text.
 * @typedef {{ transform: SvgResultTransform, applies: (content: unknown) => boolean }} LegacyResultNormalization
 */

/**
 * The plan a current Session's Results are admitted with at Load: none, or,
 * for each Result the legacy normalization applies to, that normalization;
 * for each Result of a Session that saved no records of Python's paint (a
 * Session 27-44, `SourceFacts`), the records of the edits it shows
 * (`recordSavedEditBases`);
 * for each Result the records of Python's key on the rows a Session
 * renamed without one (`recordSavedLegendKeys`); and for each Result whose
 * text shows a color ID while the Session has captioned color rules, the names
 * of its rule rows (`nameSavedRuleRows`).
 * @param {readonly { content?: unknown }[]} results
 * @param {FeatureCatalogAdmission} catalogAdmission
 * @param {SavedResultEdits | null} edits
 * @param {LegacyResultNormalization | null} [legacy]
 * @returns {SvgMutationPlan}
 */
export const createSavedResultPlan = (results, catalogAdmission, edits, legacy = null) => {
  const editsShown = edits && [
    edits.featureColorOverrides, edits.featureStrokeOverrides, edits.legendColorOverrides, edits.legendStrokeOverrides
  ].some((overrides) => Object.keys(overrides || {}).length > 0);
  const needsRecords = Boolean(edits && editsShown && !edits.resultsRecordPythonPaint);
  const normalizes = (/** @type {{ content?: unknown }} */ result) => Boolean(legacy?.applies(result?.content));
  const legacyNormalizationCount = results.filter(normalizes).length;
  const renames = edits ? savedLegendRenames(edits) : new Map();
  const captioned = Boolean(edits?.rules.some((rule) => text(rule.cap)));
  const namesRuleRows = (/** @type {{ content?: unknown }} */ result) => captioned && SAVED_COLOR_ID.test(String(result?.content ?? ''));
  const rekeyed = edits?.legendRekeys ?? new Map();
  if (legacyNormalizationCount === 0 && renames.size === 0 && !needsRecords && rekeyed.size === 0 && !results.some(namesRuleRows)) {
    return createEmptySvgMutationPlan(results.length);
  }
  const savedEdits = /** @type {SavedResultEdits} */ (edits);
  return Object.freeze({
    kind: 'MUTATING',
    legacyNormalizationCount,
    operationsByResult: Object.freeze(results.map((result, resultIndex) => {
      // The records find a row by the key the normalized Legend shows, so the
      // normalization comes first, then the re-key of the rows whose stored
      // edits Load moved; the paint records find a renamed row by the key it
      // shows, so they come before the key records. The names read Python's
      // fills and the renamed rows, so they come last.
      const transforms = [
        ...(legacy && normalizes(result) ? [legacy.transform] : []),
        ...(rekeyed.size > 0 ? [(/** @type {SvgResultRoot} */ svg) => rekeySavedRuleRows(svg, rekeyed)] : []),
        ...(needsRecords
          ? [(/** @type {SvgResultRoot} */ svg) => recordSavedEditBases(svg, { resultIndex, catalogAdmission, edits: savedEdits })]
          : []),
        ...(renames.size > 0 ? [(/** @type {SvgResultRoot} */ svg) => recordSavedLegendKeys(svg, renames)] : []),
        ...(namesRuleRows(result) ? [(/** @type {SvgResultRoot} */ svg) => nameSavedRuleRows(svg, savedEdits)] : [])
      ];
      return transforms.length > 0
        ? Object.freeze({ ...freezeEmptyOperations(), callerTransforms: Object.freeze(transforms) })
        : freezeEmptyOperations();
    }))
  });
};

/**
 * @param {SvgResultEntry} result
 * @param {CommittedSvgResultMetadata} metadata
 * @param {SvgMutationOperations | null} operations
 * @param {{ sanitizer?: SvgSanitizer, parser?: typeof DOMParser, resultIndex: number, sourceClass: string, legendRows?: readonly SvgLegendRowFacts[] | null }} options
 */
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
  const admitted = admittedResult(result, sanitizeAdmittedSvg(result, sanitizer, { phase, resultIndex }));
  const sanitized = admitted.content;
  if (!hasSanitizedSvgEnvelope(sanitized)) {
    throw new Error('The diagram engine returned malformed SVG content.');
  }
  // Label DOM identity is reconstructed by the one mounted-preview Label binder.
  // Keep its intent in the plan, but do not require transient mounted hooks on a
  // detached Worker SVG before PreviewRuntime has installed the candidate.
  if (!hasDetachedOperations(operations)) {
    return commitCatalogBackedResult(admitted, metadata);
  }

  const svg = parseSanitizedSvg(sanitized, parser);
  recordStructuralMetric('applicationSvgParseCount', 1, { phase, resultIndex });
  const index = createLazyMutationIndex(svg, { phase, resultIndex });
  applyFeatureOperations(index, operations);
  const rowsChanged = [applyLegendOperations(index, operations, {
    mayBeAbsent: (caption) => legendRowMayBeAbsent(legendRows, resultIndex, caption)
  }), nameLegendRows(svg)].includes(true);
  // Every transform runs; one that returns false changed nothing.
  const transformed = operations.callerTransforms
    .map((transform) => transform(svg, { result: admitted, resultIndex }) !== false).includes(true);
  // A Result the plan leaves as it is keeps its bytes (OV-345).
  if (!index.changed && !rowsChanged && !transformed) {
    return commitCatalogBackedResult(admitted, metadata);
  }
  const content = serializeAdmittedSvg(svg, { phase, resultIndex });
  return commitCatalogBackedResult({ ...admitted, content }, metadata);
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
    // `requireAlignedCatalogAdmission` above throws for a null admission.
    metadataFromCatalogAdmission(/** @type {FeatureCatalogAdmission} */ (catalogAdmission), resultIndex, sourceClass, featureTypes, legendRows),
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

/**
 * @param {SvgResultEntry} result
 * @param {{ sanitizer?: SvgSanitizer, parser?: typeof DOMParser, transformSvg?: SvgResultTransform | null, sourceClass: string, resultIndex: number, selectedFeatureTypes: readonly string[] | null }} options
 */
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
  const { content } = sanitizeAdmittedSvg(result, sanitizer, {
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
