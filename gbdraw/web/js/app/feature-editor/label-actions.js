// @ts-check
/** @import { NoticeHandle, ShowNotice } from '../app-messages.js' */
/** @import { DrawingState, DiagramMode } from '../../state.js' */
import { RULE_INPUT_CHANGE_WORDS, diagnosticError, normalizeUserFacingError } from '../../utils/error-normalization.js';
import { DRAWN_SELECTOR_QUALIFIERS, drawnSelectorUnknown, ruleFeaturePayload } from '../rule-matching.js';
import { featureDrawnContext, featureDrawnInResult, getFeatureVisibilityOverride } from '../../services/feature-visibility.js';
import { parseLabelOverrideTsv, serializeLabelOverrideRows } from '../../services/label-override-table.js';
import { replaceBulkLabelEdits } from '../../services/drawing/label-tables.js';
import { escapeRegexLiteral } from '../../services/feature-selector.js';
import {
  featureIdentityKeyOf,
  normalizeFeatureOverrideLabelText,
  updateFeatureOverride
} from '../../services/feature-placement.js';
import { resultCatalogFeatures, resultRenderedFeatures } from '../../services/feature-catalog.js';
import { FEATURE_SELECTOR, getFeatureIdentity } from './svg-actions.js';
import { downloadTextFile } from '../../services/text-download.js';
import { editorTableRef } from '../../services/file-imports.js';
import { defaultFeatureRendering } from '../../utils/feature-rendering.js';
import { readFileText } from '../../services/file-content-cache.js';
import { ownEntry } from '../../utils/own-entry.js';
import { normalizeTsvCell } from '../../utils/tsv-cell.js';
import { COMPARISON_LEGEND_SELECTOR } from '../../services/legend-svg.js';
import {
  ANNOTATION_SELECTOR, DEPTH_GROUPS, DINUCLEOTIDE_CONTENT_GROUPS, LENGTH_BAR_ID, REGION_LABEL_FRONT_SELECTOR,
  TICK_GROUPS, trackGroupSelector
} from '../../services/drawn-groups.js';
import { joinSentences, msg, t, tm, tn } from '../../utils/i18n.js';

/**
 * A feature row of the displayed Result's catalog or of the extracted list.
 * Python owns the fields (R7), so the editor reads them as `unknown`.
 * @typedef {Record<string, unknown>} LabelFeature
 */
/**
 * What the label code reads of a render request's `diagramOptions`.
 * @typedef {{ configOverrides?: Readonly<Record<string, unknown>> | null, featureShapes?: Readonly<Record<string, unknown>> | null }} LabelDiagramOptions
 */
/** @typedef {{ x: number, y: number }} LabelPoint */
/** @typedef {{ id: string, x: number, y: number, groupId: string }} FeatureCentroid */
/** @typedef {{ all: FeatureCentroid[], byGroup: Map<string, FeatureCentroid[]> }} FeatureGeometry */
/** @typedef {{ key: string, idx: number, text: string, sourceText: string, featureId: string, kind: string, draftText: string }} EditableLabelEntry */
/**
 * The popup's label fields as one Apply reads them.
 * @typedef {{ text: string, sourceText: string, visibility: string, hasEditableLabel: boolean }} LabelPopupEdit
 */
/** @typedef {(options?: { forceReflow?: boolean }) => unknown} LabelApply */
/** @typedef {{ textEl: SVGTextElement, referencePoint: LabelPoint, candidates: FeatureCentroid[], candidateById: Map<string, FeatureCentroid>, kind: string }} LabelMeta */
/**
 * The input a Label TSV change reads: a file input, or the object the retry
 * passes in its place.
 * @typedef {{ files: ArrayLike<File> | null, value: string }} LabelTableInput
 */

// The groups whose text is never a feature label. A Linear Result draws a
// feature label with the central baseline, and so does it draw a region label
// (in its annotation's group, and its copy in front of the leader lines, which
// names that group) and the depth and GC content axis ticks (OV-1291).
export const EXCLUDED_GROUP_SELECTOR = [
  '#legend',
  '#feature_legend',
  COMPARISON_LEGEND_SELECTOR,
  '#horizontal_legend',
  '#vertical_legend',
  `#${LENGTH_BAR_ID}`,
  trackGroupSelector(TICK_GROUPS),
  ANNOTATION_SELECTOR,
  REGION_LABEL_FRONT_SELECTOR,
  trackGroupSelector(DEPTH_GROUPS),
  trackGroupSelector(DINUCLEOTIDE_CONTENT_GROUPS)
].join(', ');
const EDITABLE_LABEL_SELECTOR = 'text[data-label-editable="true"]';
const LABEL_BINDING_SCHEMA_ATTRIBUTE = 'data-gbdraw-label-binding-schema';
const LABEL_BINDING_SCHEMA = '1';
const LABEL_FEATURE_ID_ATTRIBUTE = 'data-label-feature-id';
const LABEL_VISIBILITY_PREVIEW_ATTRIBUTE = 'data-gbdraw-label-visibility-preview';
// The labels whose visibility a projection sets: the labels the editor bound,
// and the labels Python bound to their features, which the editor binds only
// when it needs them, so a feature a rule hides hides its label (OV-35).
const VISIBILITY_LABEL_SELECTOR = `${EDITABLE_LABEL_SELECTOR}, `
  + `text[${LABEL_BINDING_SCHEMA_ATTRIBUTE}="${LABEL_BINDING_SCHEMA}"][${LABEL_FEATURE_ID_ATTRIBUTE}]`;

// The locator of a feature whose label binding failed: its rendered ID, type,
// and one-based span from the displayed feature catalog (R6).
/**
 * @param {unknown} featureId
 * @param {(LabelFeature | null | undefined)[] | null | undefined} features
 * @param {number} count
 */
const labelBindingLocator = (featureId, features, count) => {
  const key = normalizeKeyToken(featureId);
  const feature = (Array.isArray(features) ? features : [])
    .find((item) => normalizeKeyToken(item?.svg_id) === key);
  const start = Number(feature?.start);
  const end = Number(feature?.end);
  return {
    featureId: String(feature?.svg_id || featureId).trim(),
    ...(feature?.type ? { featureType: String(feature.type) } : {}),
    ...(Number.isSafeInteger(start) && Number.isSafeInteger(end) ? { featureStart: start + 1, featureEnd: end } : {}),
    ...(count > 1 ? { featureCount: count } : {})
  };
};

// A required feature binds exactly one label, an optional one at most one.
// A missing forced label is a setting the user can change (LABEL_NOT_DRAWN);
// two labels for one feature repeat for the same inputs (RENDER_FAILED).
/**
 * @param {Iterable<Element> | null | undefined} labelElements
 * @param {Iterable<unknown> | null | undefined} requiredFeatureIds
 * @param {{ allowMissing?: boolean, features?: (LabelFeature | null | undefined)[] | null }} [options]
 */
export const requireUniqueEditableLabelBindings = (
  labelElements,
  requiredFeatureIds,
  { allowMissing = false, features = [] } = {}
) => {
  /** @type {Map<string, string>} */
  const required = new Map();
  Array.from(requiredFeatureIds || []).forEach((featureId) => {
    const key = normalizeKeyToken(featureId);
    if (key && !required.has(key)) required.set(key, String(featureId).trim());
  });
  if (required.size === 0) return;
  const counts = new Map(Array.from(required.keys(), (key) => [key, 0]));
  Array.from(labelElements || []).forEach((element) => {
    const key = normalizeKeyToken(element?.getAttribute?.(LABEL_FEATURE_ID_ATTRIBUTE));
    // `counts.has(key)` was just checked, so `get` returns a number.
    if (counts.has(key)) counts.set(key, /** @type {number} */ (counts.get(key)) + 1);
  });
  /** @param {(count: number) => boolean} predicate */
  const failing = (predicate) => Array.from(counts).filter(([, count]) => predicate(count))
    .map(([key]) => required.get(key));
  const missing = allowMissing ? [] : failing((count) => count === 0);
  const [code, featureIds] = missing.length
    ? ['LABEL_NOT_DRAWN', missing]
    : ['RENDER_FAILED', failing((count) => count > 1)];
  if (featureIds.length === 0) return;
  throw diagnosticError(code, {
    ...(code === 'LABEL_NOT_DRAWN' ? { reason: 'FORCED_LABEL' } : {}),
    ...labelBindingLocator(featureIds[0], features, featureIds.length)
  }, /** @type {{ stage?: string, operation?: string }} */ ({ stage: 'render', operation: 'generate' }));
};

// Whether a label is drawn only when it fits inside its feature, which only
// Python measures: a label edit then asks for the label rerender.
/** @param {LabelDiagramOptions | null | undefined} diagramOptions */
const labelsEmbeddedOnly = (diagramOptions) => diagramOptions?.configOverrides?.['labels.rendering'] === 'embedded_only';

// A feature drawn as underlay has no label (gbdraw/features/factory.py), and
// with Label Rendering = Embedded Only no label is drawn that does not fit
// inside its feature (gbdraw/labels/). Returns 'underlay', 'embedded_only', or
// '' for a render request's diagram options. Generate and the label rerender
// ask this of the features their own request draws.
/**
 * @param {LabelFeature | null | undefined} feature
 * @param {LabelDiagramOptions | null | undefined} diagramOptions
 */
export const labelDrawingBlocker = (feature, diagramOptions) => {
  const featureType = String(feature?.type || '').trim();
  if (featureType && (
    diagramOptions?.featureShapes?.[featureType] || defaultFeatureRendering(featureType)
  ) === 'underlay') return 'underlay';
  return labelsEmbeddedOnly(diagramOptions) ? 'embedded_only' : '';
};
// Why the diagram draws no label for a feature (`labelAbsenceReason`): the one
// table of reasons, as the sentence the Label Not Shown dialog renders and the
// popup note sentence that names it. The reason is read as `string` (it is ''
// when none applies), so a lookup is `undefined` for what the table lacks.
// The texts are marked with `msg` and shown with `tm` (OV-1679).
/** @type {Readonly<Record<string, { sentence: string, note: string } | undefined>>} */
export const LABEL_ABSENCE_REASONS = Object.freeze({
  hidden: { sentence: msg('The feature is hidden.'), note: msg('No label shown (the feature is hidden).') },
  underlay: { sentence: msg('Labels are not drawn for features drawn as "Underlay".'), note: msg('No label shown (drawn as Underlay).') },
  scope_none: {
    sentence: msg('Labels are set to "None" ("Show Labels" or "Label Mode"), so only a feature with Label visibility "On" has one.'),
    note: msg('No label shown (Show Labels = None).')
  },
  scope_first: {
    sentence: msg('"Show Labels" is "First Record Only", so labels are drawn only in the first record unless a feature has Label visibility "On".'),
    note: msg('No label shown (Show Labels = First Record Only).')
  },
  scope_orthogroup_top: {
    sentence: msg('"Show Labels" is "Top Similarity Group Record", so labels are drawn only for the features that setting selects, unless a feature has Label visibility "On".'),
    note: msg('No label shown (Show Labels = Top Similarity Group Record).')
  },
  whitelist: {
    sentence: msg('A label whitelist is set, and labels are drawn only for the features it lists.'),
    note: msg('No label shown (not in the label whitelist).')
  },
  blacklist: {
    sentence: msg('A label blacklist is set, and a label whose text contains one of its keywords is not drawn.'),
    note: msg('No label shown (matches the label blacklist).')
  },
  embedded_only: {
    sentence: msg('With "Label Rendering" = "Embedded Only", a label is drawn only when it fits inside its feature.'),
    note: msg('No label shown (does not fit; Label Rendering = Embedded Only).')
  }
});
// No reason found: the Show Labels scope and the label filters decide the rest.
const LABEL_ABSENCE_UNKNOWN = msg('The current Show Labels and label filter settings draw none.');

// The label display scope and filters of a render request that leave a feature
// with Default label visibility unlabeled: '' or a LABEL_ABSENCE_REASONS key.
// Python decides which feature a filter removes; a set filter is named as the
// reason once no other reason applies.
/**
 * @param {LabelFeature | null | undefined} feature
 * @param {LabelDiagramOptions | null | undefined} diagramOptions
 */
const labelScopeFilterReason = (feature, diagramOptions) => {
  const overrides = diagramOptions?.configOverrides || {};
  const scope = overrides['labels.circular.scope'] ?? overrides['labels.linear.scope'];
  if (scope === 'none') return 'scope_none';
  if (scope === 'first' && Number(feature?.record_idx) > 0) return 'scope_first';
  // Python picks the selected features (orthogroup_label_eligibility); the Web
  // names the setting and does not repeat that rule.
  if (scope === 'orthogroup_top') return 'scope_orthogroup_top';
  if (editorTableRef(diagramOptions, 'whitelist')) return 'whitelist';
  const blacklist = overrides['labels.filtering.blacklist_keywords'];
  return Array.isArray(blacklist) && blacklist.length > 0 ? 'blacklist' : '';
};

/** @param {string | null} value @param {number} [fallback] */
const toNumber = (value, fallback = 0) => {
  const parsed = Number.parseFloat(String(value));
  return Number.isFinite(parsed) ? parsed : fallback;
};

/** @param {unknown} value */
const normalizeKeyToken = (value) => String(value ?? '').trim().toLowerCase();
/** @param {unknown} name @param {string} [fallback] */
const makeSafeFilename = (name, fallback = 'gbdraw') => {
  const cleaned = String(name || '')
    .replace(/[^\w.-]+/g, '_')
    .replace(/^_+|_+$/g, '');
  return cleaned || fallback;
};
/**
 * @param {SVGSVGElement | null | undefined} svg
 * @param {SVGGraphicsElement | null | undefined} element
 * @param {number} x
 * @param {number} y
 * @returns {LabelPoint}
 */
const getSvgPoint = (svg, element, x, y) => {
  if (!svg || !element) return { x, y };
  const point = svg.createSVGPoint();
  point.x = x;
  point.y = y;
  const ctm = element.getCTM();
  if (!ctm) return { x, y };
  const transformed = point.matrixTransform(ctm);
  return { x: transformed.x, y: transformed.y };
};

/**
 * @param {LabelPoint | null | undefined} point
 * @returns {point is LabelPoint}
 */
const isFinitePoint = (point) =>
  !!point && Number.isFinite(point.x) && Number.isFinite(point.y);

/** @param {Element | null | undefined} textPathEl */
const getTextPathHref = (textPathEl) =>
  String(textPathEl?.getAttribute('href') || textPathEl?.getAttribute('xlink:href') || '').trim();

/**
 * @param {SVGSVGElement} svg
 * @param {Element} textEl
 * @returns {SVGPathElement | null}
 */
const resolveEmbeddedLabelPathElement = (svg, textEl) => {
  const textPathEl = textEl?.querySelector?.('textPath');
  if (!textPathEl || !svg) return null;
  const href = getTextPathHref(textPathEl);
  if (href.startsWith('#')) {
    const linkedEl = svg.getElementById(href.slice(1));
    // The tag name `path` makes it an SVGPathElement.
    if (linkedEl?.tagName?.toLowerCase() === 'path') return /** @type {SVGPathElement} */ (linkedEl);
  }
  const prev = textEl.previousElementSibling;
  if (prev?.tagName?.toLowerCase() === 'path') return /** @type {SVGPathElement} */ (prev);
  return null;
};

/** @param {SVGSVGElement} svg @param {Element} textEl */
const getEmbeddedLabelAnchor = (svg, textEl) => {
  const pathEl = resolveEmbeddedLabelPathElement(svg, textEl);
  if (!pathEl) return null;
  try {
    const totalLength = pathEl.getTotalLength();
    if (Number.isFinite(totalLength) && totalLength > 0) {
      const midpoint = pathEl.getPointAtLength(totalLength / 2);
      return getSvgPoint(svg, pathEl, midpoint.x, midpoint.y);
    }
  } catch {
    // Fall back to path bbox center when path length APIs are unavailable.
  }
  try {
    const bbox = pathEl.getBBox();
    return getSvgPoint(svg, pathEl, bbox.x + bbox.width / 2, bbox.y + bbox.height / 2);
  } catch {
    return null;
  }
};

/** @param {SVGSVGElement} svg @param {SVGGraphicsElement} element @returns {LabelPoint} */
const getElementCenter = (svg, element) => {
  try {
    const bbox = element.getBBox();
    return getSvgPoint(svg, element, bbox.x + bbox.width / 2, bbox.y + bbox.height / 2);
  } catch {
    if (element?.tagName?.toLowerCase() === 'text') {
      const anchor = getEmbeddedLabelAnchor(svg, element);
      // isFinitePoint is true only for a non-null point.
      if (isFinitePoint(anchor)) return anchor;
    }
    return { x: 0, y: 0 };
  }
};

/** @param {Element} textEl */
const getLabelText = (textEl) => {
  const textPath = textEl.querySelector('textPath');
  if (textPath) return textPath.textContent || '';
  return textEl.textContent || '';
};

/** @param {Element} textEl @param {unknown} value */
const setLabelText = (textEl, value) => {
  const nextText = String(value ?? '');
  const textPath = textEl.querySelector('textPath');
  if (textPath) {
    textPath.textContent = nextText;
    return;
  }
  textEl.textContent = nextText;
};

/** @param {Element} textEl */
const hasExcludedAncestor = (textEl) => Boolean(textEl.closest(EXCLUDED_GROUP_SELECTOR));

/** @param {Element} svg @param {SVGTextElement} textEl @returns {SVGLineElement | null} */
const getPhasedCircularFeatureLine = (svg, textEl) => {
  const textGroup = textEl.closest('g[id="label_text"], g[id^="label_text_"]');
  if (!textGroup) return null;
  const textGroups = Array.from(svg.querySelectorAll('g[id="label_text"], g[id^="label_text_"]'));
  const leaderGroups = Array.from(svg.querySelectorAll('g[id="label_leaders"], g[id^="label_leaders_"]'));
  const groupIndex = textGroups.indexOf(textGroup);
  if (groupIndex < 0 || groupIndex >= leaderGroups.length) return null;
  const featureId = String(textEl.getAttribute(LABEL_FEATURE_ID_ATTRIBUTE) || '').trim();
  if (!featureId) return null;
  // The leader layer draws two parts per label, both bound to its feature:
  // bend to letters, then bend to feature.
  const leaderLines = leaderGroups[groupIndex].querySelectorAll(`line[${LABEL_FEATURE_ID_ATTRIBUTE}="${CSS.escape(featureId)}"]`);
  return /** @type {SVGLineElement | undefined} */ (leaderLines[1]) || null;
};

/** @param {SVGSVGElement} svg @param {SVGTextElement} textEl */
const getCircularFeatureAnchor = (svg, textEl) => {
  const inLabelsGroup = Boolean(textEl.closest('g[id="labels"], g[id^="labels_"]'));
  /** @type {SVGGraphicsElement | null} */
  let featureLine = null;
  if (inLabelsGroup) {
    const prev = textEl.previousElementSibling;
    const prev2 = prev ? prev.previousElementSibling : null;
    const lines = [prev, prev2].filter(
      (candidate) => candidate && candidate.tagName && candidate.tagName.toLowerCase() === 'line'
    );
    // The filter above keeps only `line` elements.
    featureLine = /** @type {SVGLineElement | undefined} */ (lines[0]) || null;
  } else {
    featureLine = getPhasedCircularFeatureLine(svg, textEl);
  }
  if (!featureLine) return null;

  const x2 = toNumber(featureLine.getAttribute('x2'), NaN);
  const y2 = toNumber(featureLine.getAttribute('y2'), NaN);
  if (!Number.isFinite(x2) || !Number.isFinite(y2)) return null;
  return getSvgPoint(svg, featureLine, x2, y2);
};

/** @param {SVGSVGElement} svg @param {SVGTextElement} textEl @param {DiagramMode} mode @returns {LabelPoint | null} */
const getLabelReferencePoint = (svg, textEl, mode) => {
  if (mode === 'circular') {
    const circularAnchor = getCircularFeatureAnchor(svg, textEl);
    if (isFinitePoint(circularAnchor)) return circularAnchor;
  }
  const embeddedAnchor = getEmbeddedLabelAnchor(svg, textEl);
  if (isFinitePoint(embeddedAnchor)) return embeddedAnchor;
  const center = getElementCenter(svg, textEl);
  return isFinitePoint(center) ? center : null;
};

/** @param {SVGSVGElement} svg @param {DiagramMode} mode @returns {SVGTextElement[]} */
const collectEditableLabelElements = (svg, mode) => {
  /** @type {Set<SVGTextElement>} */
  const labels = new Set();
  // Every selector below names `text` elements (`text > textPath` takes the
  // parent), which are SVGTextElement.

  if (mode === 'circular') {
    /** @type {NodeListOf<SVGTextElement>} */ (svg.querySelectorAll('g[id="labels"] text, g[id^="labels_"] text')).forEach((textEl) => {
      if (hasExcludedAncestor(textEl)) return;
      labels.add(textEl);
    });
    /** @type {NodeListOf<SVGTextElement>} */ (svg.querySelectorAll('g[id="label_text"] text, g[id^="label_text_"] text')).forEach((textEl) => {
      if (hasExcludedAncestor(textEl)) return;
      labels.add(textEl);
    });
    svg.querySelectorAll('text > textPath').forEach((textPathEl) => {
      const parentText = textPathEl.parentElement;
      if (!parentText) return;
      if (hasExcludedAncestor(parentText)) return;
      // The selector `text > textPath` makes the parent a `text` element; the
      // DOM typing of `parentElement` is HTMLElement.
      labels.add(/** @type {SVGTextElement} */ (/** @type {unknown} */ (parentText)));
    });
  } else {
    /** @type {NodeListOf<SVGTextElement>} */ (svg.querySelectorAll('text[dominant-baseline="central"]')).forEach((textEl) => {
      if (hasExcludedAncestor(textEl)) return;
      labels.add(textEl);
    });
  }

  return Array.from(labels).sort((a, b) => {
    const aCenter = getLabelReferencePoint(svg, a, mode) || { x: 0, y: 0 };
    const bCenter = getLabelReferencePoint(svg, b, mode) || { x: 0, y: 0 };
    if (Math.abs(aCenter.y - bCenter.y) > 1) return aCenter.y - bCenter.y;
    return aCenter.x - bCenter.x;
  });
};

// OV-923: what the label binding reads of a mounted root besides its text
// and children: feature identities, group ids, and geometry. A hover
// highlight (style) or a selection (class) changes none of them.
const LABEL_BINDING_INPUT_ATTRIBUTES = Object.freeze([
  'id', 'data-gbdraw-feature-id', 'data-gbdraw-rendered-feature-id', 'display', 'transform', 'viewBox',
  'preserveAspectRatio', 'd', 'points', 'x', 'y', 'dx', 'dy', 'x1', 'y1', 'x2', 'y2', 'cx', 'cy', 'r', 'rx',
  'ry', 'width', 'height', 'href', 'startOffset', 'side', 'textLength', 'lengthAdjust', 'font-size',
  'font-family', 'font-weight', 'font-style', 'letter-spacing', 'text-anchor', 'dominant-baseline'
]);

/**
 * The one mounted root whose label binding is current: bound, and nothing the
 * binding reads changed on it since. The binding measures every feature
 * (getBBox), so opening the editor or a feature popup reuses a current one.
 * Without MutationObserver no binding counts as current.
 */
const createLabelBindingMarker = () => {
  /** @type {{ root: Element, viewport: string, stale: boolean, observer: MutationObserver } | null} */
  let current = null;
  // The root's CSS box sets its viewBox scale, which the binding's coordinates include.
  /** @param {Element} root */
  const viewportOf = (root) => `${root.clientWidth}x${root.clientHeight}`;
  return {
    /** @param {Element} root */
    isCurrent: (root) => {
      if (!current || current.root !== root) return false;
      if (current.observer.takeRecords().length > 0) current.stale = true;
      return !current.stale && current.viewport === viewportOf(root);
    },
    /** @param {Element} root */
    markCurrent: (root) => {
      current?.observer.disconnect();
      current = null;
      if (typeof MutationObserver !== 'function') return;
      const observer = new MutationObserver(() => {
        if (entry) entry.stale = true;
        observer.disconnect();
      });
      const entry = { root, viewport: viewportOf(root), stale: false, observer };
      observer.observe(root, {
        subtree: true, childList: true, characterData: true, attributes: true,
        attributeFilter: [...LABEL_BINDING_INPUT_ATTRIBUTES]
      });
      current = entry;
    }
  };
};

/** @param {SVGSVGElement} svg @returns {FeatureGeometry} */
const collectFeatureGeometry = (svg) => {
  /** @type {Map<string, { id: string, x: number, y: number, n: number, groupId: string }>} */
  const grouped = new Map();
  // FEATURE_SELECTOR names path, polygon, and rect elements, which are SVGGraphicsElement.
  /** @type {NodeListOf<SVGGraphicsElement>} */ (svg.querySelectorAll(FEATURE_SELECTOR)).forEach((el) => {
    const id = getFeatureIdentity(el);
    if (!id) return;
    const center = getElementCenter(svg, el);
    const groupId = el.closest('g[id]')?.id || '';
    if (!grouped.has(id)) {
      grouped.set(id, { id, x: 0, y: 0, n: 0, groupId });
    }
    // Set just above when absent.
    const item = /** @type {{ id: string, x: number, y: number, n: number, groupId: string }} */ (grouped.get(id));
    item.x += center.x;
    item.y += center.y;
    item.n += 1;
  });

  /** @type {FeatureCentroid[]} */
  const all = [];
  /** @type {Map<string, FeatureCentroid[]>} */
  const byGroup = new Map();
  grouped.forEach((item) => {
    if (item.n <= 0) return;
    const centroid = {
      id: item.id,
      x: item.x / item.n,
      y: item.y / item.n,
      groupId: item.groupId
    };
    all.push(centroid);
    if (!byGroup.has(centroid.groupId)) {
      byGroup.set(centroid.groupId, []);
    }
    // Set just above when absent.
    /** @type {FeatureCentroid[]} */ (byGroup.get(centroid.groupId)).push(centroid);
  });

  return { all, byGroup };
};

/** @param {FeatureGeometry | null | undefined} featureGeometry @param {Element} textEl @param {DiagramMode} mode */
const getFeatureCandidatesForLabel = (featureGeometry, textEl, mode) => {
  if (!featureGeometry || featureGeometry.all.length === 0) return [];
  const labelGroupId = textEl.closest('g[id]')?.id || '';
  let candidates = featureGeometry.all;
  if (mode === 'linear' && labelGroupId) {
    const grouped = featureGeometry.byGroup.get(labelGroupId);
    if (grouped && grouped.length > 0) {
      candidates = grouped;
    }
  }
  return candidates;
};

/** @param {DiagramMode} mode @param {string} kind */
const getDistanceThreshold = (mode, kind) => {
  if (mode === 'linear') return kind === 'embedded' ? 700 : 540;
  return kind === 'embedded' ? 520 : 420;
};

/** @param {LabelPoint} referencePoint @param {LabelPoint} candidate @param {DiagramMode} mode */
const computeCandidateDistance = (referencePoint, candidate, mode) => {
  const dx = candidate.x - referencePoint.x;
  const dy = candidate.y - referencePoint.y;
  return mode === 'linear' ? Math.abs(dx) + Math.abs(dy) * 0.6 : Math.hypot(dx, dy);
};

/**
 * @param {SVGSVGElement} svg
 * @param {SVGTextElement[]} labelElements
 * @param {FeatureGeometry | null | undefined} featureGeometry
 * @param {DiagramMode} mode
 */
const assignFeatureIdsToLabels = (svg, labelElements, featureGeometry, mode) => {
  /** @type {Map<SVGTextElement, string>} */
  const assignments = new Map();
  if (!featureGeometry || featureGeometry.all.length === 0) return assignments;

  const featureIds = new Set(featureGeometry.all.map((feature) => feature.id));
  /** @type {Set<string>} */
  const usedFeatureIds = new Set();
  // `filter(Boolean)` drops the nulls the callback returns.
  const labelMeta = /** @type {LabelMeta[]} */ (labelElements
    .map((textEl) => {
      const referencePoint = getLabelReferencePoint(svg, textEl, mode);
      if (!isFinitePoint(referencePoint)) return null;
      const candidates = getFeatureCandidatesForLabel(featureGeometry, textEl, mode);
      return {
        textEl,
        referencePoint,
        candidates,
        candidateById: new Map(candidates.map((candidate) => [candidate.id, candidate])),
        kind: textEl.querySelector('textPath') ? 'embedded' : 'regular'
      };
    })
    .filter(Boolean));

  labelMeta.forEach((meta) => {
    const existingId = String(meta.textEl.getAttribute('data-label-feature-id') || '').trim();
    if (!existingId || !featureIds.has(existingId)) return;
    if (!meta.candidateById.has(existingId)) return;
    if (usedFeatureIds.has(existingId)) return;
    assignments.set(meta.textEl, existingId);
    usedFeatureIds.add(existingId);
  });

  /** @type {{ textEl: SVGTextElement, featureId: string, distance: number }[]} */
  const edges = [];
  labelMeta.forEach((meta) => {
    if (assignments.has(meta.textEl)) return;
    const threshold = getDistanceThreshold(mode, meta.kind);
    meta.candidates.forEach((candidate) => {
      if (usedFeatureIds.has(candidate.id)) return;
      const distance = computeCandidateDistance(meta.referencePoint, candidate, mode);
      if (!Number.isFinite(distance) || distance > threshold) return;
      edges.push({ textEl: meta.textEl, featureId: candidate.id, distance });
    });
  });
  edges.sort((a, b) => a.distance - b.distance);

  const assignedLabels = new Set(assignments.keys());
  edges.forEach((edge) => {
    if (assignedLabels.has(edge.textEl)) return;
    if (usedFeatureIds.has(edge.featureId)) return;
    assignments.set(edge.textEl, edge.featureId);
    assignedLabels.add(edge.textEl);
    usedFeatureIds.add(edge.featureId);
  });

  return assignments;
};

/**
 * @typedef {object} FeatureLabelActionsOptions
 * @property {Record<string, any>} state App state (state.js; not yet typed).
 * @property {((reason: string) => boolean) | null} [commitActiveResultEdit]
 *   The preview owner's commit of an edit to the displayed Result (R1, R13).
 * @property {import('../rule-matching.js').RulePreparation} rulePreparation
 * @property {(value?: any) => { value: any }} ref Vue `ref`
 * @property {<T>(getter: () => T) => { value: T }} computed Vue `computed`
 * @property {((source: any, callback: (...args: any[]) => void, options?: Record<string, any>) => () => void) | null} [watch] Vue `watch`
 * @property {() => Promise<any>} [nextTick] Vue `nextTick`
 * @property {() => ({ diagramOptions?: Record<string, any> } | null)} [getCommittedRequest]
 *   The committed canonical request (Python owns the option fields, R7).
 * @property {((feature: Record<string, any>, mode: string, options?: Record<string, any>) => any) | null} [setFeatureVisibility]
 *   The visibility owner's transition: Show feature and label sets Feature visibility through it.
 * @property {() => boolean | Promise<boolean | { error: any }>} [prepareDrawnFeatureMatches]
 *   The visibility owner's preparation of the rule matches `resolveFeatureDrawn` reads (Label On asks the blocker of it).
 * @property {(payload: Record<string, any>, options?: Record<string, any>) => Promise<Record<string, any>>} evaluateLabelRules
 *   The root's Python evaluation of Label TSV rows against the displayed labels (R7).
 * @property {ShowNotice} [showNotice] The root's notice banner: a table import or export reports there.
 */

/** @param {FeatureLabelActionsOptions} options */
export const createFeatureLabelActions = ({
  state,
  commitActiveResultEdit = null,
  rulePreparation,
  ref,
  computed,
  watch = null,
  nextTick = () => Promise.resolve(),
  getCommittedRequest = () => null,
  setFeatureVisibility = null,
  prepareDrawnFeatureMatches = () => true,
  evaluateLabelRules,
  showNotice = () => {}
}) => {
  const {
    mode,
    results,
    selectedResultIndex,
    svgContainer,
    editableLabels,
    extractedFeatures,
    clickedFeature,
    labelTextScopeDialog,
    hiddenLabelTextDialog,
    labelOnDialog,
    autoLabelReflowEnabled,
    labelReflowRequestSeq,
    labelReflowForceRequestSeq,
    labelReflowLastError,
    labelReflowProcessing
  } = state;

  // The label edits of the identity-keyed draft: Label visibility, label text,
  // and the label's source text (design Q4); Feature visibility stays. The
  // bulk label edits become `bulk`.
  /** @type {readonly ('labelVisibility' | 'labelText' | 'labelSourceText')[]} */
  const LABEL_FIELDS = Object.freeze(['labelVisibility', 'labelText', 'labelSourceText']);
  /** @param {DrawingState} drawing @param {Record<string, string>} [bulk] */
  const clearOverrides = (drawing, bulk = {}) => {
    Object.values(drawing.featureOverrides).forEach((row) => {
      updateFeatureOverride(drawing.slices.featureEdits, row, Object.fromEntries(LABEL_FIELDS.map((field) => [field, null])));
    });
    replaceBulkLabelEdits(drawing.slices.labelTables, bulk);
  };

  // A mounted label names its feature by rendered ID; the displayed Result's
  // catalog binds that ID to the feature's source identity (R3). `drawing`
  // owns the displayed Result; its mode names the catalog's records.
  /** @param {DrawingState} drawing */
  const displayedFeatures = (drawing) => resultRenderedFeatures(state, drawing.mode);
  /**
   * @param {DrawingState} drawing
   * @param {unknown} featureIdRaw
   * @param {Map<string, LabelFeature> | null} [displayed]
   * @returns {LabelFeature | null}
   */
  const labelFeature = (drawing, featureIdRaw, displayed = displayedFeatures(drawing)) => {
    const featureId = String(featureIdRaw || '').trim();
    if (!featureId) return null;
    if (displayed) return displayed.get(featureId) || null;
    const key = normalizeKeyToken(featureId);
    return (extractedFeatures.value || []).find((/** @type {LabelFeature | null | undefined} */ feature) => normalizeKeyToken(feature?.svg_id) === key) || null;
  };
  /** @param {DrawingState} drawing @param {LabelFeature | null | undefined} feature */
  const rowOf = (drawing, feature) => {
    const key = featureIdentityKeyOf(feature);
    return key ? drawing.featureOverrides[key] || null : null;
  };
  /** @param {DrawingState} drawing @param {unknown} featureId @param {Map<string, LabelFeature> | null} [displayed] */
  const labelRow = (drawing, featureId, displayed) => rowOf(drawing, labelFeature(drawing, featureId, displayed));

  // Whether the label rerender after this edit draws the feature: the last
  // Generate's request with the current edits and rules (F-3, Owner Q2). The
  // feature visibility owner's resolver answers as Generate does.
  /** @param {DrawingState} drawing */
  const drawnContext = (drawing) => featureDrawnContext(drawing, { diagramOptions: getCommittedRequest()?.diagramOptions });
  /**
   * @param {DrawingState} drawing
   * @param {LabelFeature | null | undefined} feature
   * @param {ReturnType<typeof drawnContext>} [context]
   */
  const featureHidden = (drawing, feature, context = drawnContext(drawing)) => (
    featureDrawnInResult(feature, context, resultCatalogFeatures(state, drawing.mode)) === false
  );
  // Owner decisions Q1 and Q2 (2026-10-04): Label visibility On takes effect
  // only when the diagram can draw the label. Returns 'hidden', 'underlay', or
  // 'embedded_only', or '' when the label is drawn.
  /**
   * @param {DrawingState} drawing
   * @param {LabelFeature | null | undefined} feature
   * @param {LabelDiagramOptions | null | undefined} diagramOptions
   */
  const labelOnBlocker = (drawing, feature, diagramOptions) => (
    featureHidden(drawing, feature) ? 'hidden' : labelDrawingBlocker(feature, diagramOptions)
  );
  // Why the displayed Result has no label for the feature: the blocker, unless
  // a scope or filter reason applies to Default label visibility first.
  // '' or a LABEL_ABSENCE_REASONS key.
  /**
   * @param {DrawingState} drawing
   * @param {LabelFeature | null | undefined} feature
   * @param {LabelDiagramOptions | null | undefined} diagramOptions
   */
  const labelAbsenceReason = (drawing, feature, diagramOptions) => {
    const blocker = labelOnBlocker(drawing, feature, diagramOptions);
    if (blocker === 'hidden' || blocker === 'underlay'
      || normalizeVisibilityMode(rowOf(drawing, feature)?.labelVisibility) === 'on') return blocker;
    return labelScopeFilterReason(feature, diagramOptions) || blocker;
  };

  const commitLabelEdit = () => commitActiveResultEdit?.('feature-label');

  // R13: Generate and the label rerender clear the notices of the previous
  // label build through this port when they start; a rerender also clears its
  // last failure, as a queued request does, and so does a Generate that commits
  // its Results (OV-36). The rerender still reports its own failure (R1(c)).
  const clearLabelBuildNotices = ({ rerender = false } = {}) => {
    if (rerender) labelReflowLastError.value = null;
  };

  const queueLabelReflow = (force = false) => {
    labelReflowLastError.value = null;
    if (force) {
      labelReflowForceRequestSeq.value += 1;
      return;
    }
    if (!autoLabelReflowEnabled.value) return;
    labelReflowRequestSeq.value += 1;
  };

  /** @param {unknown} value */
  const normalizeVisibilityMode = (value) => {
    const normalized = String(value || '').trim().toLowerCase();
    return normalized === 'on' || normalized === 'off' ? normalized : 'default';
  };

  /** @param {DrawingState} drawing @param {Element} svg */
  const resetLabelsToSourceText = (drawing, svg) => {
    let changed = false;
    const displayed = displayedFeatures(drawing);
    svg.querySelectorAll(EDITABLE_LABEL_SELECTOR).forEach((textEl) => {
      const sourceText = labelRow(drawing, textEl.getAttribute('data-label-feature-id'), displayed)?.labelSourceText
        ?? textEl.getAttribute('data-label-source-text');
      if (sourceText === null) return;
      if (getLabelText(textEl) === sourceText) return;
      setLabelText(textEl, sourceText);
      textEl.setAttribute('data-label-source-text', sourceText);
      changed = true;
    });
    return changed;
  };

  // The label text intent on one SVG: each editable label shows its feature
  // override, else the bulk override of its source text, else that source
  // text. A label edited, imported, undone, or reset while another Result was
  // displayed therefore shows the current intent here (R3).
  /** @param {DrawingState} drawing @param {Element} svg */
  const projectLabelTextIntent = (drawing, svg) => {
    let changed = false;
    const displayed = displayedFeatures(drawing);
    svg.querySelectorAll(EDITABLE_LABEL_SELECTOR).forEach((textEl) => {
      const row = labelRow(drawing, textEl.getAttribute(LABEL_FEATURE_ID_ATTRIBUTE), displayed);
      const sourceText = row?.labelSourceText
        ?? textEl.getAttribute('data-label-source-text') ?? getLabelText(textEl);
      textEl.setAttribute('data-label-source-text', sourceText);
      // The source text is free text, so the bulk table is read by own key.
      const bulkText = ownEntry(drawing.labelTextBulkOverrides, sourceText);
      const desiredText = row?.labelText ?? bulkText ?? sourceText;
      if (getLabelText(textEl) === desiredText) return;
      setLabelText(textEl, desiredText);
      changed = true;
    });
    return changed;
  };

  /** @param {Element | null | undefined} svg @param {Element | null | undefined} textEl */
  const resolveCompleteLabelVisualUnit = (svg, textEl) => {
    if (!svg || !textEl) return null;
    const featureId = String(textEl.getAttribute(LABEL_FEATURE_ID_ATTRIBUTE) || '').trim();
    if (!featureId) return null;
    if (textEl.getAttribute(LABEL_BINDING_SCHEMA_ATTRIBUTE) !== LABEL_BINDING_SCHEMA) {
      return null;
    }
    const parts = Array.from(svg.querySelectorAll(`[${LABEL_FEATURE_ID_ATTRIBUTE}]`))
      .filter((element) => (
        String(element.getAttribute(LABEL_FEATURE_ID_ATTRIBUTE) || '').trim() === featureId
      ));
    if (!parts.includes(textEl)) return null;
    return { featureId, parts };
  };

  /** @param {Element | null | undefined} svg @param {Element | null | undefined} textEl @param {unknown} modeRaw */
  const applyLabelVisibilityPreview = (svg, textEl, modeRaw) => {
    const visualUnit = resolveCompleteLabelVisualUnit(svg, textEl);
    if (!visualUnit) return { available: false, changed: false };
    const visibilityMode = normalizeVisibilityMode(modeRaw);
    let changed = false;
    if (visibilityMode === 'off') {
      visualUnit.parts.forEach((part) => {
        if (part.getAttribute(LABEL_VISIBILITY_PREVIEW_ATTRIBUTE) !== 'off') {
          part.setAttribute(LABEL_VISIBILITY_PREVIEW_ATTRIBUTE, 'off');
          changed = true;
        }
        if (part.getAttribute('display') !== 'none') {
          part.setAttribute('display', 'none');
          changed = true;
        }
      });
      return { available: true, changed };
    }
    visualUnit.parts.forEach((part) => {
      if (!part.hasAttribute(LABEL_VISIBILITY_PREVIEW_ATTRIBUTE)) return;
      part.removeAttribute(LABEL_VISIBILITY_PREVIEW_ATTRIBUTE);
      part.removeAttribute('display');
      changed = true;
    });
    return { available: true, changed };
  };

  // The Label visibility On and Off each mounted Result was drawn with, by
  // feature identity key: the intent of the Result's first projection. Only
  // Python decides whether Default draws a label, so a feature whose On or Off
  // returns to Default asks for the rerender, as a missing On label does
  // (UJ-01, R3).
  /** @type {WeakMap<Element, Map<string, string>>} */
  const drawnLabelVisibility = new WeakMap();
  /** @param {DrawingState} drawing */
  const labelVisibilityIntent = (drawing) => new Map(Object.entries(drawing.featureOverrides)
    .map(([key, row]) => /** @type {[string, string]} */ ([key, normalizeVisibilityMode(row?.labelVisibility)]))
    .filter(([, visibility]) => visibility !== 'default'));
  /** @param {DrawingState} drawing @param {Element} svg */
  const labelVisibilityReturnsToDefault = (drawing, svg) => {
    const intent = labelVisibilityIntent(drawing);
    const drawn = drawnLabelVisibility.get(svg);
    if (!drawn) drawnLabelVisibility.set(svg, intent);
    return Array.from(drawn?.keys() || []).some((key) => !intent.has(key));
  };

  /** @param {DrawingState} drawing @param {Element} svg */
  const applyStoredVisibilityOverridesToSvg = (drawing, svg) => {
    let changed = false;
    let unavailableOverride = labelVisibilityReturnsToDefault(drawing, svg);
    const displayed = displayedFeatures(drawing);
    const context = drawnContext(drawing);
    svg.querySelectorAll(VISIBILITY_LABEL_SELECTOR).forEach((textEl) => {
      const featureId = String(
        textEl.getAttribute(LABEL_FEATURE_ID_ATTRIBUTE) || ''
      ).trim();
      const feature = featureId ? labelFeature(drawing, featureId, displayed) : null;
      // A hidden feature hides its label through this projection too (F-3).
      const visibilityMode = !featureId
        ? 'default'
        : (featureHidden(drawing, feature, context) ? 'off' : rowOf(drawing, feature)?.labelVisibility);
      const projection = applyLabelVisibilityPreview(svg, textEl, visibilityMode);
      changed = projection.changed || changed;
      if (normalizeVisibilityMode(visibilityMode) !== 'default' && !projection.available) {
        unavailableOverride = true;
      }
    });
    return { changed, unavailableOverride };
  };

  /** @param {Element} svg */
  const refreshEditableList = (svg) => {
    /** @type {EditableLabelEntry[]} */
    const nextEntries = [];
    svg.querySelectorAll(EDITABLE_LABEL_SELECTOR).forEach((textEl, index) => {
      const key = textEl.getAttribute('data-label-key');
      if (!key) return;
      const text = getLabelText(textEl);
      const sourceText = textEl.getAttribute('data-label-source-text') || text;
      const featureId = textEl.getAttribute('data-label-feature-id') || '';
      const kind = textEl.querySelector('textPath') ? 'embedded' : 'regular';
      nextEntries.push({
        key,
        idx: index + 1,
        text,
        sourceText,
        featureId,
        kind,
        draftText: text
      });
    });
    editableLabels.value = nextEntries;
  };

  const clickedFeatureId = () => String(clickedFeature.value?.svg_id || clickedFeature.value?.id || '').trim();
  /** @param {DrawingState} drawing @param {string} featureId @returns {LabelFeature | null} */
  const featureById = (drawing, featureId) => (
    (clickedFeatureId() === featureId ? clickedFeature.value?.feat : null)
    || labelFeature(drawing, featureId)
  );

  /** @param {unknown} featureId @returns {EditableLabelEntry | null} */
  const getEditableLabelByFeatureId = (featureId) => {
    const target = normalizeKeyToken(featureId);
    if (!target) return null;
    return (
      editableLabels.value.find(
        (/** @type {EditableLabelEntry | null | undefined} */ entry) => normalizeKeyToken(entry?.featureId) === target
      ) || null
    );
  };

  const syncClickedFeatureLabelState = () => {
    const drawing = state.activeDrawing();
    if (!clickedFeature.value) return;
    const featureId = String(clickedFeature.value.svg_id || clickedFeature.value.id || '').trim();
    const entry = getEditableLabelByFeatureId(featureId);
    const row = rowOf(drawing, clickedFeature.value.feat);
    const fallbackText = row?.labelText || clickedFeature.value.labelText || clickedFeature.value.label || '';
    const fallbackSource = row?.labelSourceText || clickedFeature.value.labelSourceText
      || clickedFeature.value.label || '';
    // Per-feature edits name the feature by source identity, which a Result
    // committed before feature catalogs existed does not carry.
    const editable = Boolean(featureId && featureIdentityKeyOf(clickedFeature.value.feat));
    clickedFeature.value.labelKey = entry?.key || '';
    clickedFeature.value.labelText = entry?.text ?? fallbackText;
    clickedFeature.value.labelSourceText = entry?.sourceText ?? fallbackSource;
    clickedFeature.value.labelVisibility = normalizeVisibilityMode(row?.labelVisibility);
    clickedFeature.value.hasEditableLabel = editable;
    clickedFeature.value.labelUnavailableReason = editable
      ? ''
      : (featureId
        ? t('Generate the diagram again to edit this feature\'s label.')
        : t('No editable feature label for this feature in current diagram.'));
  };

  /**
   * @param {{
   *   requiredFeatureIds?: readonly string[],
   *   optionalFeatureIds?: readonly string[],
   *   reportedLabelBinding?: { featureIds: readonly string[], report: (error: unknown) => void } | null,
   *   queueIncompleteVisibility?: boolean,
   *   reuseCurrentBinding?: boolean
   * }} [options] `reuseCurrentBinding`: opening the editor or a feature
   *   popup, which checks no feature ids, keeps a current binding (OV-923).
   * @returns {{ changed: boolean, rerender: boolean, textChanged: boolean } | undefined} The projection; undefined when no Result is mounted.
   */
  const syncLabelEditor = ({
    requiredFeatureIds = [],
    optionalFeatureIds = [],
    reportedLabelBinding = null,
    queueIncompleteVisibility = true,
    reuseCurrentBinding = false
  } = {}) => {
    const drawing = state.activeDrawing();
    if (!svgContainer.value) return;
    const svg = svgContainer.value.querySelector('svg');
    if (!svg) return;
    boundLabelRoots.add(svg);
    const reuse = reuseCurrentBinding && !requiredFeatureIds.length && !optionalFeatureIds.length
      && !reportedLabelBinding && labelBinding.isCurrent(svg);
    // A reused binding is one made with no feature ids to check (`reuse`).
    const labelElements = reuse ? [] : bindLabels(svg, requiredFeatureIds, optionalFeatureIds);

    const projection = projectLabelIntent(drawing, svg, { queueIncompleteVisibility });
    // A label reflow keeps its Result: the same check reports to the reflow
    // after the binding completes instead of failing it.
    if (reportedLabelBinding) {
      try {
        requireUniqueEditableLabelBindings(labelElements, reportedLabelBinding.featureIds, { features: extractedFeatures.value });
      } catch (error) {
        reportedLabelBinding.report(error);
      }
    }
    return projection;
  };

  const labelBinding = createLabelBindingMarker();
  // Binds each editable label of the root to its feature, by geometry.
  /**
   * @param {SVGSVGElement} svg
   * @param {readonly string[]} requiredFeatureIds
   * @param {readonly string[]} optionalFeatureIds
   * @returns {SVGTextElement[]} the editable labels, in key order
   */
  const bindLabels = (svg, requiredFeatureIds, optionalFeatureIds) => {
    // Label intent is keyed by feature identity, not by the mounted view: a
    // Result switch, record selection, hide, or reflow changes which labels
    // this SVG shows, never which overrides exist (FE-01, R2).
    const featureGeometry = collectFeatureGeometry(svg);
    const labelElements = collectEditableLabelElements(svg, mode.value);
    const features = extractedFeatures.value;
    requireUniqueEditableLabelBindings(
      labelElements,
      [...requiredFeatureIds, ...optionalFeatureIds],
      { allowMissing: true, features }
    );
    const featureAssignments = assignFeatureIdsToLabels(svg, labelElements, featureGeometry, mode.value);
    labelElements.forEach((textEl, index) => {
      textEl.style.cursor = 'text';
      textEl.setAttribute('data-label-editable', 'true');
      textEl.setAttribute('data-label-key', `label-${index + 1}`);
      const currentText = getLabelText(textEl);
      if (!textEl.hasAttribute('data-label-source-text')) {
        textEl.setAttribute('data-label-source-text', currentText);
      }
      const featureId = featureAssignments.get(textEl);
      if (featureId) {
        textEl.setAttribute('data-label-feature-id', featureId);
      } else {
        textEl.removeAttribute('data-label-feature-id');
      }
    });
    requireUniqueEditableLabelBindings(labelElements, requiredFeatureIds, { features });
    requireUniqueEditableLabelBindings(
      labelElements,
      optionalFeatureIds,
      { allowMissing: true, features }
    );
    // Marked after the binding's own writes and before the projection, so a
    // text the projection changes counts as a change since the binding.
    labelBinding.markCurrent(svg);
    return labelElements;
  };

  // One projection of the label intent onto the mounted Result, shared by a
  // live edit, a Label TSV import, History apply, and the display of a Result.
  // `rerender`: only the label rerender can show the intent; a caller that
  // queues its own reflow (`queueIncompleteVisibility: false`) forces it then.
  // `textChanged`: a label shows another text, which only Python fits.
  /** @param {DrawingState} drawing @param {Element} svg @param {{ queueIncompleteVisibility?: boolean }} [options] */
  const projectLabelIntent = (drawing, svg, { queueIncompleteVisibility = true } = {}) => {
    const textChanged = projectLabelTextIntent(drawing, svg);
    const visibilityProjection = applyStoredVisibilityOverridesToSvg(drawing, svg);
    refreshEditableList(svg);
    syncClickedFeatureLabelState();
    const changed = textChanged || visibilityProjection.changed;
    if (changed) commitLabelEdit();
    if (queueIncompleteVisibility && visibilityProjection.unavailableOverride) {
      queueLabelReflow(true);
    }
    return { changed, rerender: visibilityProjection.unavailableOverride, textChanged };
  };

  // The mounted roots whose labels a popup or the binder bound.
  /** @type {WeakSet<Element>} */
  const boundLabelRoots = new WeakSet();
  // The label intent on the mounted Result, and the one label reflow request
  // of the edit that projects it: `reflow` places the labels (Auto Reflow),
  // `rerender` asks for the automatic rerender. The labels of a root no popup
  // or binder has bound are bound first, once (OV-237). OV-200: with Embedded
  // Only a changed label text asks for the rerender, which draws the label
  // only when it fits, as Generate does; a label hidden, or shown as Python
  // drew it, needs no such decision.
  const reconcileLabelOverrides = ({ reflow = false, rerender = false } = {}) => {
    const drawing = state.activeDrawing();
    const svg = svgContainer.value?.querySelector?.('svg');
    if (!svg) return false;
    const projection = boundLabelRoots.has(svg) || svg.querySelector(EDITABLE_LABEL_SELECTOR)
      ? projectLabelIntent(drawing, svg, { queueIncompleteVisibility: false })
      : syncLabelEditor({ queueIncompleteVisibility: false });
    const force = rerender || Boolean(projection?.rerender)
      || (Boolean(projection?.textChanged) && labelsEmbeddedOnly(getCommittedRequest()?.diagramOptions));
    if (reflow || force) queueLabelReflow(force);
    return Boolean(projection?.changed);
  };

  // A feature visibility edit shows or hides the feature's label in the same
  // action, then Auto Reflow places the labels as Generate does (F-3). An edit
  // that draws a feature the Result does not draw (`rerender`) always reruns
  // the rerender, which draws it and its label (R-5). A caller that also
  // projects the label intent (`labels`) does so here, so the edit queues one
  // request.
  const applyFeatureVisibilityToLabels = ({ reflow = true, rerender = false, labels = false } = {}) => {
    const drawing = state.activeDrawing();
    if (labels) return reconcileLabelOverrides({ reflow, rerender });
    const svg = svgContainer.value?.querySelector?.('svg');
    if (!svg) return false;
    const projection = applyStoredVisibilityOverridesToSvg(drawing, svg);
    if (projection.changed) commitLabelEdit();
    if (reflow || rerender) queueLabelReflow(rerender || projection.unavailableOverride);
    return projection.changed;
  };

  // An edit whose Result Python must draw again asks for the automatic
  // rerender, Auto Reflow on or off: a specific-color rule edit or a History
  // step that changes a Legend source (Owner decision 2026-10-06, OV-43). A
  // request while one is pending joins it. D-36: an accepted request leaves a
  // redraw owed by the shown drawing (each mode has its own, PD-OI-086) until
  // the composition root shows a rerender of that drawing drawn from its
  // current paint (`redrawShown`, OV-346); while one is owed, the shown
  // figure can lack rows Python draws, so a History step does not record it
  // as drawn (D-24 (OV-731)).
  /** @type {WeakSet<DrawingState>} */
  const redrawOwed = new WeakSet();
  const requestAutomaticRerender = () => {
    if (!svgContainer.value?.querySelector?.('svg')) return false;
    redrawOwed.add(state.activeDrawing());
    queueLabelReflow(true);
    return true;
  };
  /** @param {DrawingState} drawing */
  const isRedrawOwed = (drawing) => redrawOwed.has(drawing);
  /** @param {DrawingState} drawing */
  const redrawShown = (drawing) => { redrawOwed.delete(drawing); };

  /** @param {string | null | undefined} labelKey @param {unknown} nextTextRaw */
  const requestLabelTextChangeByKey = (labelKey, nextTextRaw) => {
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    if (!labelKey) return;
    const entry = editableLabels.value.find((/** @type {EditableLabelEntry} */ candidate) => candidate.key === labelKey);
    if (!entry) return;
    const nextText = String(nextTextRaw ?? '');
    if (entry.text === nextText) return;

    labelTextScopeDialog.show = true;
    labelTextScopeDialog.labelKey = entry.key;
    labelTextScopeDialog.newText = nextText;
    labelTextScopeDialog.sourceText = entry.sourceText || entry.text;
    labelTextScopeDialog.featureId = entry.featureId || '';
    labelTextScopeDialog.matchingCount = editableLabels.value.filter(
      (/** @type {EditableLabelEntry} */ candidate) => candidate.sourceText === (entry.sourceText || entry.text)
    ).length;
    if (clickedFeature.value?.labelKey === entry.key) {
      clickedFeature.value.labelText = nextText;
    }
  };

  const closeLabelTextScopeDialog = () => {
    labelTextScopeDialog.show = false;
    labelTextScopeDialog.labelKey = '';
    labelTextScopeDialog.newText = '';
    labelTextScopeDialog.sourceText = '';
    labelTextScopeDialog.featureId = '';
    labelTextScopeDialog.matchingCount = 0;
  };

  /** @param {unknown} featureId @param {unknown} nextTextRaw */
  const requestLabelTextChangeByFeatureId = (featureId, nextTextRaw) => {
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    const entry = getEditableLabelByFeatureId(featureId);
    if (!entry) return false;
    requestLabelTextChangeByKey(entry.key, nextTextRaw);
    return true;
  };

  // Stores one feature's Label visibility and shows it in its open popup.
  /** @param {DrawingState} drawing @param {LabelFeature | null | undefined} feature @param {unknown} modeRaw */
  const setLabelVisibilityOverride = (drawing, feature, modeRaw) => {
    const nextMode = normalizeVisibilityMode(modeRaw);
    const previousMode = normalizeVisibilityMode(rowOf(drawing, feature)?.labelVisibility);
    updateFeatureOverride(drawing.slices.featureEdits, feature, { labelVisibility: nextMode === 'default' ? null : nextMode });
    if (featureIdentityKeyOf(clickedFeature.value?.feat) === featureIdentityKeyOf(feature)) {
      clickedFeature.value.labelVisibility = nextMode;
    }
    return previousMode !== nextMode;
  };

  // One feature's label text edit: the text equal to the label's own text
  // removes the edit; any other text is one line, kept with its source text.
  /**
   * @param {DrawingState} drawing
   * @param {LabelFeature | null | undefined} feature
   * @param {unknown} labelTextRaw
   * @param {unknown} sourceTextRaw
   * @param {unknown} baselineTextRaw
   */
  const applyDirectFeatureLabelOverride = (drawing, feature, labelTextRaw, sourceTextRaw, baselineTextRaw) => {
    if (!featureIdentityKeyOf(feature)) return false;
    const nextText = String(labelTextRaw ?? '');
    const sourceText = String(sourceTextRaw ?? '').trim();
    if (nextText === String(baselineTextRaw ?? '')) {
      return updateFeatureOverride(drawing.slices.featureEdits, feature, { labelText: null, labelSourceText: null });
    }
    return updateFeatureOverride(drawing.slices.featureEdits, feature, {
      labelText: normalizeFeatureOverrideLabelText(nextText),
      ...(sourceText ? { labelSourceText: sourceText } : {})
    });
  };

  /** @param {unknown} featureId @param {string} nextText */
  const applyDirectTextToCurrentSvg = (featureId, nextText) => {
    if (!svgContainer.value) return { svg: null, changed: false };
    const svg = svgContainer.value.querySelector('svg');
    if (!svg) return { svg: null, changed: false };
    const entry = getEditableLabelByFeatureId(featureId);
    if (!entry?.key) return { svg, changed: false };
    const targetEl = svg.querySelector(`text[data-label-key="${CSS.escape(entry.key)}"]`);
    if (!targetEl) return { svg, changed: false };
    const currentText = getLabelText(targetEl);
    if (currentText === nextText) return { svg, changed: false };
    setLabelText(targetEl, nextText);
    return { svg, changed: true };
  };

  /**
   * @param {unknown} featureId
   * @param {unknown} visibilityMode
   * @returns {{ available: boolean, changed: boolean, svg?: Element }}
   */
  const applyDirectVisibilityToCurrentSvg = (featureId, visibilityMode) => {
    if (!svgContainer.value) return { available: false, changed: false };
    const svg = svgContainer.value.querySelector('svg');
    if (!svg) return { available: false, changed: false };
    const entry = getEditableLabelByFeatureId(featureId);
    if (!entry?.key) return { available: false, changed: false };
    const targetEl = svg.querySelector(`text[data-label-key="${CSS.escape(entry.key)}"]`);
    const projection = applyLabelVisibilityPreview(svg, targetEl, visibilityMode);
    return { ...projection, svg };
  };

  // Label text alone keeps Default visibility, which follows Show Labels and
  // the label filters. When they leave this feature unlabeled, Label Not Shown
  // asks before anything is written, as applying On does (Owner Q1, Q2): the
  // caller's History step stays open until the choice, so Show and Keep hidden
  // are one step each and Cancel writes nothing and records none.
  /** @type {((choice: string) => unknown) | null} */
  let answerHiddenLabelTextDialog = null;
  // The one close of the dialog. A dialog hidden without a choice (Reset label
  // edits, an applied Load Label TSV, a mode switch) is answered as Cancel:
  // the text is not written and the History step of its Apply ends (OV-462).
  const closeHiddenLabelTextDialog = (choice = 'cancel') => {
    const answer = answerHiddenLabelTextDialog;
    answerHiddenLabelTextDialog = null;
    hiddenLabelTextDialog.show = false;
    hiddenLabelTextDialog.featureId = '';
    hiddenLabelTextDialog.reason = '';
    return answer?.(choice);
  };

  // One popup Label edit: the overrides, then the displayed Result.
  /**
   * @param {DrawingState} drawing
   * @param {LabelFeature} feature
   * @param {string} featureId
   * @param {LabelPopupEdit} edit
   * @param {{ forceReflow?: boolean }} [options]
   */
  const applyPopupLabelEdit = (drawing, feature, featureId, edit, { forceReflow = false } = {}) => {
    const row = rowOf(drawing, feature);
    if (edit.text.trim() && edit.visibility === 'default'
      && normalizeVisibilityMode(row?.labelVisibility) === 'default'
      && normalizeFeatureOverrideLabelText(edit.text) !== normalizeFeatureOverrideLabelText(String(row?.labelText ?? edit.sourceText))
      && !getEditableLabelByFeatureId(featureId)) {
      closeHiddenLabelTextDialog();
      hiddenLabelTextDialog.featureId = featureId;
      hiddenLabelTextDialog.reason = labelAbsenceReason(drawing, feature, getCommittedRequest()?.diagramOptions);
      hiddenLabelTextDialog.show = true;
      const writeText = () => applyDirectFeatureLabelOverride(drawing, feature, edit.text, edit.sourceText, edit.sourceText);
      return new Promise((resolve) => {
        answerHiddenLabelTextDialog = (choice) => {
          let written;
          if (choice === 'text_only') written = writeText();
          else if (choice === 'show') {
            written = applyLabelOn(drawing, featureId, ({ forceReflow: force = false } = {}) => {
              writeText();
              setLabelVisibilityOverride(drawing, feature, 'on');
              const projection = applyDirectVisibilityToCurrentSvg(featureId, 'on');
              queueLabelReflow(force || !projection.available);
            });
          } else written = syncClickedFeatureLabelState();
          resolve(written);
          return written;
        };
      });
    }

    const visibilityChanged = setLabelVisibilityOverride(drawing, feature, edit.visibility);
    const textChanged = applyDirectFeatureLabelOverride(drawing, feature, edit.text, edit.sourceText, edit.sourceText);

    const textProjection = edit.hasEditableLabel
      ? applyDirectTextToCurrentSvg(featureId, edit.text)
      : { svg: null, changed: false };
    const visibilityProjection = visibilityChanged
      ? applyDirectVisibilityToCurrentSvg(featureId, edit.visibility)
      : { available: true, changed: false, svg: null };
    const mutatedSvg = visibilityProjection.svg || textProjection.svg;
    if (mutatedSvg && (textProjection.changed || visibilityProjection.changed)) {
      commitLabelEdit();
      syncLabelEditor({ queueIncompleteVisibility: false });
    }

    if (visibilityChanged || (!edit.hasEditableLabel && textChanged)) {
      const svg = svgContainer.value?.querySelector?.('svg');
      queueLabelReflow(forceReflow || !visibilityProjection.available
        || Boolean(svg && labelVisibilityReturnsToDefault(drawing, svg)));
      return;
    }

    if (textChanged) {
      queueLabelReflow(labelsEmbeddedOnly(getCommittedRequest()?.diagramOptions));
    }
  };

  const updateClickedFeatureLabelText = async () => {
    const drawing = state.activeDrawing();
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    if (!clickedFeature.value) return;
    const featureId = clickedFeatureId();
    if (!featureId) return;
    const feature = clickedFeature.value.feat;
    if (!featureIdentityKeyOf(feature)) return;
    const sourceText = String(clickedFeature.value.labelSourceText || clickedFeature.value.label || '');
    const edit = {
      text: String(clickedFeature.value.labelText ?? ''),
      sourceText,
      visibility: normalizeVisibilityMode(clickedFeature.value.labelVisibility),
      hasEditableLabel: Boolean(clickedFeature.value.hasEditableLabel)
    };
    // A label edit carries one line of text; clearing the text hides the
    // label, as an empty label-table text does, unless Label visibility is On.
    if (!edit.text.trim() && edit.text !== sourceText) {
      edit.text = sourceText;
      if (edit.visibility !== 'on') edit.visibility = 'off';
      clickedFeature.value.labelText = sourceText;
      clickedFeature.value.labelVisibility = edit.visibility;
    }
    /** @type {LabelApply} */
    const apply = (options) => applyPopupLabelEdit(drawing, feature, featureId, edit, options);
    if (edit.visibility === 'on' && normalizeVisibilityMode(rowOf(drawing, feature)?.labelVisibility) !== 'on') {
      return applyLabelOn(drawing, featureId, apply);
    }
    return apply();
  };

  // The Label Not Shown dialog answers the Apply that opened it.
  /** @param {string} choice */
  const handleHiddenLabelTextChoice = (choice) => {
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    return closeHiddenLabelTextDialog(choice);
  };

  // Owner decisions Q1 and Q2 (2026-10-04): applying Label visibility On asks
  // first when the diagram cannot draw the label: a hidden feature (Show
  // feature and label, or Keep feature hidden), a feature drawn as underlay,
  // and, after the label reflow, a label that does not fit with Embedded Only
  // (Keep without label). The caller's History step stays open until the
  // choice, so a choice is one step and Cancel, which leaves the label intent
  // as it was, records none. Global settings never change.
  /** @param {DrawingState} drawing @param {string} featureId @param {LabelApply} apply */
  const applyLabelOn = async (drawing, featureId, apply) => {
    // The blocker reads Python's rule matches; prepared ones answer at once.
    const prepared = /** @type {any} */ (prepareDrawnFeatureMatches());
    if (typeof prepared?.then === 'function') await prepared;
    const feature = featureById(drawing, featureId);
    const diagramOptions = getCommittedRequest()?.diagramOptions;
    const blocker = labelOnBlocker(drawing, feature, diagramOptions);
    if (!blocker) return apply();
    return confirmLabelOn(drawing, { featureId, feature, diagramOptions, blocker, apply });
  };

  /** @type {((choice: string) => void) | null} */
  let answerLabelOnDialog = null;
  /** @param {string} reason @param {LabelFeature | null} feature @returns {Promise<string>} */
  const askLabelOn = (reason, feature) => {
    answerLabelOnDialog?.('cancel');
    labelOnDialog.reason = reason;
    labelOnDialog.featureType = String(feature?.type || '');
    labelOnDialog.show = true;
    return new Promise((resolve) => { answerLabelOnDialog = resolve; });
  };

  // The Label visibility On dialog answers the Apply that opened it.
  /** @param {string} choice */
  const handleLabelOnChoice = (choice) => {
    const answer = answerLabelOnDialog;
    answerLabelOnDialog = null;
    labelOnDialog.show = false;
    answer?.(choice);
  };

  const cancelLabelOn = () => {
    syncClickedFeatureLabelState();
    return false;
  };

  /** @param {DrawingState} drawing @param {LabelFeature | null} feature */
  const captureFeatureLabelIntent = (drawing, feature) => {
    const row = rowOf(drawing, feature);
    return Object.fromEntries(LABEL_FIELDS.map((field) => [field, row?.[field] ?? null]));
  };
  /** @param {DrawingState} drawing @param {LabelFeature | null} feature @param {Record<string, string | null>} captured */
  const restoreFeatureLabelIntent = (drawing, feature, captured) => updateFeatureOverride(drawing.slices.featureEdits, feature, captured);

  // Whether the label reflow that this Apply queued drew the feature's label;
  // null when no reflow replaced the Result (it was skipped or failed, or a
  // Generate took over), so nothing is asked.
  /** @param {string} featureId */
  const reflowDrawsLabel = async (featureId) => {
    const generationKey = state.resultGenerationKey?.value;
    const before = results.value[selectedResultIndex.value];
    await nextTick();
    if (labelReflowProcessing?.value && typeof watch === 'function') {
      await /** @type {Promise<void>} */ (new Promise((resolve) => {
        const stop = watch(labelReflowProcessing, (busy) => {
          if (busy) return;
          stop();
          resolve();
        });
      }));
    }
    const result = results.value[selectedResultIndex.value];
    if (!result || result === before || labelReflowLastError.value || state.processing?.value
      || state.resultGenerationKey?.value !== generationKey) return null;
    const root = new DOMParser().parseFromString(String(result.content || ''), 'image/svg+xml').documentElement;
    const key = normalizeKeyToken(featureId);
    return Array.from(root.querySelectorAll(`text[${LABEL_FEATURE_ID_ATTRIBUTE}]`))
      .some((element) => normalizeKeyToken(element.getAttribute(LABEL_FEATURE_ID_ATTRIBUTE)) === key);
  };

  /**
   * @param {DrawingState} drawing
   * @param {{ featureId: string, feature: LabelFeature | null, diagramOptions: LabelDiagramOptions | null | undefined, blocker: string, apply: LabelApply }} choice
   */
  const confirmLabelOn = async (drawing, { featureId, feature, diagramOptions, blocker, apply }) => {
    const showFeature = blocker === 'hidden';
    if (showFeature) {
      const choice = await askLabelOn('hidden', feature);
      if (choice === 'keep') return apply();
      if (choice !== 'show') return cancelLabelOn();
    }
    // A feature shown here can still be drawn as underlay or with Embedded Only.
    const drawnAs = showFeature ? labelDrawingBlocker(feature, diagramOptions) : blocker;
    if (drawnAs === 'underlay' && await askLabelOn('underlay', feature) !== 'keep') return cancelLabelOn();
    const embeddedOnly = drawnAs === 'embedded_only';
    const labelIntent = captureFeatureLabelIntent(drawing, feature);
    const featureVisibility = getFeatureVisibilityOverride(drawing.featureOverrides, feature);
    // The port returns false for a feature without an identity key, which is what a null feature has.
    if (showFeature) setFeatureVisibility?.(/** @type {LabelFeature} */ (feature), 'on', { triggerReflow: false });
    apply({ forceReflow: embeddedOnly });
    if (!embeddedOnly) return true;
    const applied = JSON.stringify(captureFeatureLabelIntent(drawing, feature));
    const drawn = await reflowDrawsLabel(featureId);
    // A later edit of this label decides instead of this Apply.
    if (drawn !== false || JSON.stringify(captureFeatureLabelIntent(drawing, feature)) !== applied) return true;
    if (await askLabelOn('embedded_only', feature) === 'keep') return true;
    // Cancel restores the label intent and feature visibility this Apply
    // changed, and the label reflow draws the Result they describe again.
    restoreFeatureLabelIntent(drawing, feature, labelIntent);
    if (showFeature) setFeatureVisibility?.(/** @type {LabelFeature} */ (feature), featureVisibility, { triggerReflow: false });
    reconcileLabelOverrides({ rerender: true });
    return cancelLabelOn();
  };

  // The popup note when the displayed Result draws no label for the clicked
  // feature: why the diagram cannot draw it (Owner Q1, Q2), or that On shows it,
  // in one line.
  const clickedFeatureLabelHint = computed(() => {
    const drawing = state.activeDrawing();
    const clicked = clickedFeature.value;
    if (!clicked?.hasEditableLabel) return '';
    const draft = normalizeVisibilityMode(clicked.labelVisibility);
    if (draft === 'off') return '';
    void results.value; // The committed request changes with the displayed Results.
    void rulePreparation?.pending?.value; // Python's rule matches arrive when an evaluation ends.
    const featureId = clickedFeatureId();
    const feature = featureById(drawing, featureId);
    const diagramOptions = getCommittedRequest()?.diagramOptions;
    const blocker = labelOnBlocker(drawing, feature, diagramOptions);
    // A hidden feature's label is hidden with it (F-3).
    if (clicked.labelKey && blocker !== 'hidden') return '';
    /** @type {string[]} */
    const next = [];
    if (normalizeVisibilityMode(rowOf(drawing, clicked.feat)?.labelVisibility) === 'on') {
      next.push(blocker === 'hidden' ? t('On applies when the feature is shown.') : t('On applies when the label can be drawn.'));
    } else if (!blocker && draft === 'default') {
      next.push(t('Choose On to show it.'));
    } else if (!blocker) {
      return '';
    }
    const note = LABEL_ABSENCE_REASONS[labelAbsenceReason(drawing, feature, diagramOptions)]?.note;
    return joinSentences(note ? tm(note) : t('No label shown.'), ...next);
  });

  // The Label Not Shown dialog: the same reason sentence as the popup note.
  const hiddenLabelTextMessage = computed(() => joinSentences(
    t('This feature has no label in the current Result.'),
    tm(LABEL_ABSENCE_REASONS[hiddenLabelTextDialog.reason]?.sentence ?? LABEL_ABSENCE_UNKNOWN),
    t('The edited text will not appear unless you show this label.')
  ));

  // A label list edit: one line of text, or an empty text, which hides the
  // label (Label visibility Off) as an empty label-table text does.
  /** @param {DrawingState} drawing @param {LabelFeature | null} feature @param {unknown} text @param {string} sourceText */
  const applyLabelTextEdit = (drawing, feature, text, sourceText) => {
    const labelText = normalizeFeatureOverrideLabelText(String(text ?? ''));
    return updateFeatureOverride(drawing.slices.featureEdits, feature, labelText
      ? { labelText, labelSourceText: sourceText || null }
      : { labelText: null, labelVisibility: 'off', labelSourceText: sourceText || null });
  };

  /** @param {string} choice */
  const handleLabelTextScopeChoice = (choice) => {
    const drawing = state.activeDrawing();
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    if (choice === 'cancel' || !labelTextScopeDialog.show) {
      closeLabelTextScopeDialog();
      return;
    }
    if (!svgContainer.value) {
      closeLabelTextScopeDialog();
      return;
    }
    const svg = svgContainer.value.querySelector('svg');
    if (!svg) {
      closeLabelTextScopeDialog();
      return;
    }

    const targetKey = String(labelTextScopeDialog.labelKey || '');
    const sourceText = String(labelTextScopeDialog.sourceText || '');
    const newText = String(labelTextScopeDialog.newText ?? '');
    const featureId = String(labelTextScopeDialog.featureId || '');

    if (choice === 'all') {
      let matchedCount = 0;
      let hasUntrackableMatch = false;
      const displayed = displayedFeatures(drawing);
      svg.querySelectorAll(EDITABLE_LABEL_SELECTOR).forEach((/** @type {Element} */ textEl) => {
        const candidateSource = textEl.getAttribute('data-label-source-text') || '';
        if (candidateSource !== sourceText) return;
        matchedCount += 1;
        const candidate = labelFeature(drawing, textEl.getAttribute('data-label-feature-id'), displayed);
        if (featureIdentityKeyOf(candidate)) {
          applyLabelTextEdit(drawing, candidate, newText, candidateSource);
        } else {
          hasUntrackableMatch = true;
        }
        setLabelText(textEl, newText);
      });
      /** @type {Record<string, string>} */
      const bulk = { ...drawing.labelTextBulkOverrides };
      if (hasUntrackableMatch || matchedCount === 0) bulk[sourceText] = newText;
      else delete bulk[sourceText];
      replaceBulkLabelEdits(drawing.slices.labelTables, bulk);
    } else if (choice === 'single' && targetKey) {
      const targetEl = svg.querySelector(`text[data-label-key="${CSS.escape(targetKey)}"]`);
      if (targetEl) {
        setLabelText(targetEl, newText);
      }
      const feature = labelFeature(drawing, featureId);
      if (featureIdentityKeyOf(feature)) applyLabelTextEdit(drawing, feature, newText, sourceText);
    } else {
      closeLabelTextScopeDialog();
      return;
    }

    commitLabelEdit();
    closeLabelTextScopeDialog();
    syncLabelEditor();
    queueLabelReflow();
  };

  /** @type {{ drawing: DrawingState, notice: NoticeHandle } | null} The notice of the last Label TSV load, and the drawing it edited. */
  let labelImportNotice = null;
  const resetAllLabelTextOverrides = () => {
    const drawing = state.activeDrawing();
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    // D-06: the notice of the Label TSV edits this removes ends with them.
    const importNotice = labelImportNotice;
    if (importNotice && importNotice.drawing === drawing) {
      importNotice.notice.dismiss();
      labelImportNotice = null;
    }
    const svg = svgContainer.value?.querySelector('svg');
    if (svg) resetLabelsToSourceText(drawing, svg);
    clearOverrides(drawing);
    if (!svgContainer.value) {
      editableLabels.value = [];
      closeLabelTextScopeDialog();
      closeHiddenLabelTextDialog();
      return;
    }
    if (!svg) return;
    applyStoredVisibilityOverridesToSvg(drawing, svg);

    closeLabelTextScopeDialog();
    closeHiddenLabelTextDialog();
    commitLabelEdit();
    queueLabelReflow(Boolean(syncLabelEditor({ queueIncompleteVisibility: false })?.rerender));
  };

  // A Label TSV row selects labels in every Result of a batch. The displayed
  // Result binds its labels as the editor does; the committed content of
  // another Result carries the renderer's label bindings (B6).
  /** @param {DrawingState} drawing @param {Element} element @param {LabelFeature | null} feature */
  const labelSourceText = (drawing, element, feature) => (
    rowOf(drawing, feature)?.labelSourceText
    ?? element.getAttribute('data-label-source-text') ?? getLabelText(element)
  );
  /** @param {DrawingState} drawing @param {SVGSVGElement} svg */
  const labelImportTargets = (drawing, svg) => {
    const elements = collectEditableLabelElements(svg, mode.value);
    const assignments = assignFeatureIdsToLabels(svg, elements, collectFeatureGeometry(svg), mode.value);
    const displayed = displayedFeatures(drawing);
    const targets = elements.map((element) => {
      const featureId = assignments.get(element) || '';
      const feature = featureId ? labelFeature(drawing, featureId, displayed) : null;
      return { featureId, feature, sourceText: labelSourceText(drawing, element, feature) };
    });
    const seen = new Set(targets.map(({ feature }) => featureIdentityKeyOf(feature)).filter(Boolean));
    results.value.forEach((/** @type {{ content?: unknown } | null | undefined} */ result, /** @type {number} */ index) => {
      if (index === selectedResultIndex.value || typeof result?.content !== 'string') return;
      const features = resultRenderedFeatures(state, drawing.mode, index);
      const root = new DOMParser().parseFromString(result.content, 'image/svg+xml').documentElement;
      root.querySelectorAll(`text[${LABEL_FEATURE_ID_ATTRIBUTE}]`).forEach((element) => {
        const featureId = String(element.getAttribute(LABEL_FEATURE_ID_ATTRIBUTE) || '').trim();
        const feature = featureId ? features?.get(featureId) || null : null;
        const key = featureIdentityKeyOf(feature);
        if (!key || seen.has(key)) return;
        seen.add(key);
        targets.push({ featureId, feature, sourceText: labelSourceText(drawing, element, feature) });
      });
    });
    return targets;
  };

  let labelImportRevision = 0;
  const labelImportFailure = ref(null);
  const canRetryLabelImportFailure = computed(() => {
    const drawing = state.activeDrawing();
    return Boolean(labelImportFailure.value
      && state.errorLog?.value === labelImportFailure.value.error
      && labelImportFailure.value.revision === labelImportRevision
      && rulePreparation.isCurrent(labelImportFailure.value.snapshot)
      && labelImportFailure.value.intent === labelIntentSignature(drawing));
  });
  const retryLabelImportFailure = () => canRetryLabelImportFailure.value ? labelImportFailure.value.retry() : false;
  const editLabelImportFailure = () => {
    if (canRetryLabelImportFailure.value) labelImportFailure.value.input?.click?.();
  };
  /** @param {DrawingState} drawing */
  const labelIntentSignature = (drawing) => JSON.stringify([drawing.featureOverrides, drawing.labelTextBulkOverrides]);
  /** @param {{ target: LabelTableInput | null, sourceInput?: LabelTableInput }} event A file input's change event, or the retry's stand-in for it. */
  const loadLabelOverrideTable = async (event) => {
    const drawing = state.activeDrawing();
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    const input = event?.target;
    const file = input?.files?.[0];
    if (!file) return;
    const sourceSvg = svgContainer.value?.querySelector('svg') || null;
    const before = rulePreparation.snapshot();
    const labelIntent = labelIntentSignature(drawing);
    const revision = ++labelImportRevision;
    const previousError = state.errorLog?.value;

    try {
      const text = await readFileText(file);
      if (input.files?.[0] !== file || (svgContainer.value?.querySelector('svg') || null) !== sourceSvg) return;
      if (state.sessionOperationAvailability?.()) return state.sessionOperationAvailability();
      const rows = parseLabelOverrideTsv(text);
      const loadedRowsText = tn(rows.length, 'Loaded {count} row.', 'Loaded {count} rows.');

      const svg = sourceSvg;
      const labels = svg ? labelImportTargets(drawing, svg) : [];
      const evaluation = await evaluateLabelRules({
        kind: 'label', rules: rows,
        features: labels.map((entry) => ruleFeaturePayload(
          entry.feature || { type: '', svg_id: entry.featureId }, entry.sourceText
        ))
      });
      if (state.sessionOperationAvailability?.()) return state.sessionOperationAvailability();
      if (revision !== labelImportRevision || input.files?.[0] !== file) return;
      // OV-447, OV-449 (Owner decisions): the table was evaluated against a
      // state that changed meanwhile, so the load is declined and names the
      // change. An edit of the label intent stays; it is named first, since
      // it also redraws the diagram and makes the rule preparation stale.
      const changed = labelIntent !== labelIntentSignature(drawing)
        ? t('Not applied: label text, label visibility, or feature visibility changed while the table loaded. The changes were kept. Load the table again.')
        : (svgContainer.value?.querySelector('svg') || null) !== sourceSvg
          ? t('Not applied: the displayed diagram changed while the table loaded. Load the table again.')
          : !rulePreparation.isCurrent(before)
            ? t('Not applied: {inputs} changed while the table loaded. Load the table again.', { inputs: tm(RULE_INPUT_CHANGE_WORDS) }) : '';
      if (changed) {
        showNotice({ text: joinSentences(loadedRowsText, changed), tone: 'warning' });
        return;
      }
      if (!svg) {
        showNotice({ text: joinSentences(loadedRowsText, t('No diagram is shown, so nothing was applied.')), tone: 'warning' });
        return;
      }
      // A `location` or `record_location` row cannot be decided for a label
      // whose drawn values the Result does not record (feature catalog 3 or 4
      // on a cropped, reverse-complemented, or rotated record). The import is
      // declined whole, so the label edits it would replace stay (R4).
      // `drawnSelectorUnknown` holds only for a feature that has a `drawnSelector` key, so `feature` is set below.
      const undecided = labels.filter(({ feature }) => drawnSelectorUnknown(feature) && rows.some((row) => (
        DRAWN_SELECTOR_QUALIFIERS.has(row.qualifier.toLowerCase())
        && (row.recordId === '*' || row.recordId === String(/** @type {LabelFeature} */ (feature).record_id ?? ''))
        && (row.featureType === '*' || row.featureType === String(/** @type {LabelFeature} */ (feature).type ?? ''))
      ))).length;
      if (undecided > 0) {
        showNotice({ text: joinSentences(loadedRowsText, t('Not applied.')), items: [
          tn(undecided,
        'location and record_location rows cannot be matched to {count} label on a cropped, reverse-complemented, or rotated record until the diagram records where it drew it. Generate the diagram again, then load the table.',
        'location and record_location rows cannot be matched to {count} labels on cropped, reverse-complemented, or rotated records until the diagram records where it drew them. Generate the diagram again, then load the table.')
        ] });
        return;
      }

      // The import replaces the label intent once; the displayed Result shows
      // it now and every other Result when it is displayed (B6, R3). A table
      // that applies to no label is declined whole, so the label edits stay
      // and History records no step (OV-131).
      let skippedNonTrackableCount = 0;
      const applicable = labels.flatMap((entry, index) => {
        const matchedRow = rows[evaluation.winners[index]];
        if (!matchedRow) return [];
        const tracked = Boolean(featureIdentityKeyOf(entry.feature));
        // A global `label` row follows the source text; any other row its
        // feature, as that feature's label edit (design Q4).
        if (matchedRow.isGlobalLabelRule ? !entry.sourceText : !tracked) {
          skippedNonTrackableCount += 1;
          return [];
        }
        return [{ entry, matchedRow, tracked }];
      });
      if (applicable.length === 0) {
        showNotice({ text: joinSentences(loadedRowsText, t('Not applied; the existing label edits were kept.')), items: [
          skippedNonTrackableCount > 0
            ? tn(skippedNonTrackableCount, 'The {count} matched label lacked a feature key.', 'The {count} matched labels lacked a feature key.')
            : t('No row matched a label of the diagram.')
        ] });
        return;
      }
      clearOverrides(drawing, Object.fromEntries(applicable
        .filter(({ matchedRow }) => matchedRow.isGlobalLabelRule)
        .map(({ entry, matchedRow }) => [entry.sourceText, String(matchedRow.labelText ?? '')])));
      applicable.forEach(({ entry, matchedRow, tracked }) => {
        if (!matchedRow.isGlobalLabelRule) {
          applyLabelTextEdit(drawing, entry.feature, matchedRow.labelText, entry.sourceText);
        } else if (tracked) {
          updateFeatureOverride(drawing.slices.featureEdits, entry.feature, { labelSourceText: entry.sourceText || null });
        }
      });

      closeLabelTextScopeDialog();
      closeHiddenLabelTextDialog();
      queueLabelReflow(Boolean(syncLabelEditor({ queueIncompleteVisibility: false })?.rerender));

      const message = joinSentences(loadedRowsText, tn(applicable.length, 'Applied to {count} label.', 'Applied to {count} labels.'));
      const skipped = skippedNonTrackableCount > 0
        ? [tn(skippedNonTrackableCount, '{count} match lacked a feature key and was not applied.', '{count} matches lacked a feature key and were not applied.')]
        : [];
      if (state.errorLog?.value === labelImportFailure.value?.error) state.errorLog.value = null;
      labelImportFailure.value = null;
      const notice = showNotice({ text: message, items: skipped });
      labelImportNotice = notice ? { drawing, notice } : null;
    } catch (error) {
      if (revision !== labelImportRevision || input.files?.[0] !== file || !rulePreparation.isCurrent(before)
        || state.errorLog?.value !== previousError) return;
      const model = normalizeUserFacingError(error, { operation: 'evaluateRules', stage: 'resource-staging' });
      if (state.errorLog) state.errorLog.value = model;
      labelImportFailure.value = { error: model, snapshot: before, intent: labelIntent, revision,
        input: event.sourceInput || input,
        retry: () => loadLabelOverrideTable({ target: { files: [file], value: '' }, sourceInput: event.sourceInput || input }) };
    } finally {
      if (input?.files?.[0] === file) input.value = '';
    }
  };

  // Export Label TSV writes the label rules the next Generate sends: a saved
  // label table, or the bulk label edits as `* * label` rows. Per-feature label
  // edits are identity rows; Export Feature Edits TSV writes them (design Q4 6.4).
  const downloadLabelOverrideTable = () => {
    const drawing = state.activeDrawing();
    const savedTable = serializeLabelOverrideRows(drawing.canonicalLabelOverrideRows?.value);
    const rows = savedTable
      ? savedTable.trimEnd().split('\n')
      : Object.keys(drawing.labelTextBulkOverrides).sort((a, b) => a.localeCompare(b)).filter(Boolean)
        .map((sourceText) => `*\t*\tlabel\t^${escapeRegexLiteral(sourceText)}$\t${
          normalizeTsvCell(drawing.labelTextBulkOverrides[sourceText])}`);
    if (rows.length === 0) {
      showNotice({ text: joinSentences(t('There are no label rules to export.'), t('Export Feature Edits TSV writes per-feature label edits.')) });
      return;
    }
    const selectedIdx = selectedResultIndex.value;
    const resultName =
      selectedIdx >= 0 && selectedIdx < results.value.length
        ? String(results.value[selectedIdx]?.name || '')
        : '';
    const outputName = `${makeSafeFilename(resultName, 'gbdraw')}.label_table.tsv`;
    downloadTextFile(outputName, `${rows.join('\n')}\n`, 'text/tab-separated-values');
  };

  return {
    applyFeatureVisibilityToLabels,
    clearLabelBuildNotices,
    clickedFeatureLabelHint,
    closeHiddenLabelTextDialog,
    closeLabelTextScopeDialog,
    hiddenLabelTextMessage,
    downloadLabelOverrideTable,
    loadLabelOverrideTable, canRetryLabelImportFailure, retryLabelImportFailure, editLabelImportFailure,
    getEditableLabelByFeatureId,
    handleHiddenLabelTextChoice,
    handleLabelOnChoice,
    handleLabelTextScopeChoice,
    requestLabelTextChangeByFeatureId,
    requestLabelTextChangeByKey,
    reconcileLabelOverrides,
    requestAutomaticRerender,
    isRedrawOwed,
    redrawShown,
    resetAllLabelTextOverrides,
    syncClickedFeatureLabelState,
    syncLabelEditor,
    updateClickedFeatureLabelText
  };
};
