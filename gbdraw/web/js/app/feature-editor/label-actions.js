// @ts-check
import { diagnosticError, normalizeUserFacingError } from '../../utils/error-normalization.js';
import { DRAWN_SELECTOR_QUALIFIERS, drawnSelectorUnknown, ruleFeaturePayload } from '../rule-matching.js';
import { featureDrawnContext, featureDrawnInResult, getFeatureVisibilityOverride } from '../../services/feature-visibility.js';
import { parseLabelOverrideTsv, serializeLabelOverrideRows } from '../../services/label-override-table.js';
import { escapeRegexLiteral } from '../../services/feature-selector.js';
import {
  featureIdentityKeyOf,
  normalizeFeatureOverrideLabelText,
  updateFeatureOverride
} from '../../services/feature-placement.js';
import { resultCatalogFeatures, resultRenderedFeatures } from '../../services/feature-catalog.js';
import { FEATURE_SELECTOR, getFeatureIdentity } from './svg-actions.js';
import { downloadTextFile } from '../../services/text-download.js';
import { defaultFeatureRendering } from '../../utils/feature-rendering.js';
import { readFileText } from '../../services/file-content-cache.js';
import { normalizeTsvCell } from '../../utils/tsv-cell.js';
import { COMPARISON_LEGEND_SELECTOR } from '../../services/legend-svg.js';

export const EXCLUDED_GROUP_SELECTOR = [
  '#legend',
  '#feature_legend',
  COMPARISON_LEGEND_SELECTOR,
  '#horizontal_legend',
  '#vertical_legend',
  '#length_bar',
  'g[data-gbdraw-slot-renderer="ticks"]',
  'g[id="tick"]',
  'g[id^="tick_"]'
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
export const requireUniqueEditableLabelBindings = (
  labelElements,
  requiredFeatureIds,
  { allowMissing = false, features = [] } = {}
) => {
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

// A feature drawn as underlay has no label (gbdraw/features/factory.py), and
// with Label Rendering = Embedded Only no label is drawn that does not fit
// inside its feature (gbdraw/labels/). Returns 'underlay', 'embedded_only', or
// '' for a render request's diagram options. Generate and the label rerender
// ask this of the features their own request draws.
export const labelDrawingBlocker = (feature, diagramOptions) => {
  const featureType = String(feature?.type || '').trim();
  if (featureType && (
    diagramOptions?.featureShapes?.[featureType] || defaultFeatureRendering(featureType)
  ) === 'underlay') return 'underlay';
  return diagramOptions?.configOverrides?.['labels.rendering'] === 'embedded_only' ? 'embedded_only' : '';
};
// Why the diagram draws no label for a feature (`labelAbsenceReason`): the one
// table of sentences the popup note and the Label Not Shown dialog both render.
export const LABEL_ABSENCE_REASONS = Object.freeze({
  hidden: ' The feature is hidden.',
  underlay: ' Labels are not drawn for features drawn as "Underlay".',
  scope_none: ' Labels are set to "None" ("Show Labels" or "Label Mode"), so only a feature with Label visibility "On" has one.',
  scope_first: ' "Show Labels" is "First Record Only", so labels are drawn only in the first record unless a feature has Label visibility "On".',
  scope_orthogroup_top: ' "Show Labels" is "Top Similarity Group Record", so labels are drawn only for the features that setting selects, unless a feature has Label visibility "On".',
  whitelist: ' A label whitelist is set, and labels are drawn only for the features it lists.',
  blacklist: ' A label blacklist is set, and a label whose text contains one of its keywords is not drawn.',
  embedded_only: ' With "Label Rendering" = "Embedded Only", a label is drawn only when it fits inside its feature.'
});
const LABEL_ABSENT_PREFIX = 'This feature has no label in the current Result.';
// No reason found: the Show Labels scope and the label filters decide the rest.
const LABEL_ABSENCE_UNKNOWN = ' The current Show Labels and label filter settings draw none.';

// The label display scope and filters of a render request that leave a feature
// with Default label visibility unlabeled: '' or a LABEL_ABSENCE_REASONS key.
// Python decides which feature a filter removes; a set filter is named as the
// reason once no other reason applies.
const labelScopeFilterReason = (feature, diagramOptions) => {
  const overrides = diagramOptions?.configOverrides || {};
  const scope = overrides['labels.circular.scope'] ?? overrides['labels.linear.scope'];
  if (scope === 'none') return 'scope_none';
  if (scope === 'first' && Number(feature?.record_idx) > 0) return 'scope_first';
  // Python picks the selected features (orthogroup_label_eligibility); the Web
  // names the setting and does not repeat that rule.
  if (scope === 'orthogroup_top') return 'scope_orthogroup_top';
  if (diagramOptions?.labelWhitelistFile) return 'whitelist';
  const blacklist = overrides['labels.filtering.blacklist_keywords'];
  return Array.isArray(blacklist) && blacklist.length > 0 ? 'blacklist' : '';
};

const toNumber = (value, fallback = 0) => {
  const parsed = Number.parseFloat(value);
  return Number.isFinite(parsed) ? parsed : fallback;
};

const normalizeKeyToken = (value) => String(value ?? '').trim().toLowerCase();
const makeSafeFilename = (name, fallback = 'gbdraw') => {
  const cleaned = String(name || '')
    .replace(/[^\w.-]+/g, '_')
    .replace(/^_+|_+$/g, '');
  return cleaned || fallback;
};
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

const isFinitePoint = (point) =>
  Boolean(point) && Number.isFinite(point.x) && Number.isFinite(point.y);

const getTextPathHref = (textPathEl) =>
  String(textPathEl?.getAttribute('href') || textPathEl?.getAttribute('xlink:href') || '').trim();

const resolveEmbeddedLabelPathElement = (svg, textEl) => {
  const textPathEl = textEl?.querySelector?.('textPath');
  if (!textPathEl || !svg) return null;
  const href = getTextPathHref(textPathEl);
  if (href.startsWith('#')) {
    const linkedEl = svg.getElementById(href.slice(1));
    if (linkedEl?.tagName?.toLowerCase() === 'path') return linkedEl;
  }
  const prev = textEl.previousElementSibling;
  if (prev?.tagName?.toLowerCase() === 'path') return prev;
  return null;
};

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

const getElementCenter = (svg, element) => {
  try {
    const bbox = element.getBBox();
    return getSvgPoint(svg, element, bbox.x + bbox.width / 2, bbox.y + bbox.height / 2);
  } catch {
    if (element?.tagName?.toLowerCase() === 'text') {
      const anchor = getEmbeddedLabelAnchor(svg, element);
      // isFinitePoint is true only for a non-null point.
      if (isFinitePoint(anchor)) return /** @type {{ x: number, y: number }} */ (anchor);
    }
    return { x: 0, y: 0 };
  }
};

const getLabelText = (textEl) => {
  const textPath = textEl.querySelector('textPath');
  if (textPath) return textPath.textContent || '';
  return textEl.textContent || '';
};

const setLabelText = (textEl, value) => {
  const nextText = String(value ?? '');
  const textPath = textEl.querySelector('textPath');
  if (textPath) {
    textPath.textContent = nextText;
    return;
  }
  textEl.textContent = nextText;
};

const hasExcludedAncestor = (textEl) => Boolean(textEl.closest(EXCLUDED_GROUP_SELECTOR));

const getPhasedCircularFeatureLine = (svg, textEl) => {
  const textGroup = textEl.closest('g[id="label_text"], g[id^="label_text_"]');
  if (!textGroup) return null;
  const textGroups = Array.from(svg.querySelectorAll('g[id="label_text"], g[id^="label_text_"]'));
  const leaderGroups = Array.from(svg.querySelectorAll('g[id="label_leaders"], g[id^="label_leaders_"]'));
  const groupIndex = textGroups.indexOf(textGroup);
  if (groupIndex < 0 || groupIndex >= leaderGroups.length) return null;
  const textIndex = Array.from(textGroup.querySelectorAll('text')).indexOf(textEl);
  if (textIndex < 0) return null;
  const leaderLines = Array.from(leaderGroups[groupIndex].querySelectorAll('line'));
  return leaderLines[(2 * textIndex) + 1] || null;
};

const getCircularFeatureAnchor = (svg, textEl) => {
  const inLabelsGroup = Boolean(textEl.closest('g[id="labels"], g[id^="labels_"]'));
  /** @type {Element | null} */
  let featureLine = null;
  if (inLabelsGroup) {
    const prev = textEl.previousElementSibling;
    const prev2 = prev ? prev.previousElementSibling : null;
    const lines = [prev, prev2].filter(
      (candidate) => candidate && candidate.tagName && candidate.tagName.toLowerCase() === 'line'
    );
    featureLine = lines[0] || null;
  } else {
    featureLine = getPhasedCircularFeatureLine(svg, textEl);
  }
  if (!featureLine) return null;

  const x2 = toNumber(featureLine.getAttribute('x2'), NaN);
  const y2 = toNumber(featureLine.getAttribute('y2'), NaN);
  if (!Number.isFinite(x2) || !Number.isFinite(y2)) return null;
  return getSvgPoint(svg, featureLine, x2, y2);
};

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

const collectEditableLabelElements = (svg, mode) => {
  const labels = new Set();

  if (mode === 'circular') {
    svg.querySelectorAll('g[id="labels"] text, g[id^="labels_"] text').forEach((textEl) => {
      if (hasExcludedAncestor(textEl)) return;
      labels.add(textEl);
    });
    svg.querySelectorAll('g[id="label_text"] text, g[id^="label_text_"] text').forEach((textEl) => {
      if (hasExcludedAncestor(textEl)) return;
      labels.add(textEl);
    });
    svg.querySelectorAll('text > textPath').forEach((textPathEl) => {
      const parentText = textPathEl.parentElement;
      if (!parentText) return;
      if (hasExcludedAncestor(parentText)) return;
      labels.add(parentText);
    });
  } else {
    svg.querySelectorAll('text[dominant-baseline="central"]').forEach((textEl) => {
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

const collectFeatureGeometry = (svg) => {
  const grouped = new Map();
  svg.querySelectorAll(FEATURE_SELECTOR).forEach((el) => {
    const id = getFeatureIdentity(el);
    if (!id) return;
    const center = getElementCenter(svg, el);
    const groupId = el.closest('g[id]')?.id || '';
    if (!grouped.has(id)) {
      grouped.set(id, { id, x: 0, y: 0, n: 0, groupId });
    }
    const item = grouped.get(id);
    item.x += center.x;
    item.y += center.y;
    item.n += 1;
  });

  const all = [];
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
    byGroup.get(centroid.groupId).push(centroid);
  });

  return { all, byGroup };
};

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

const getDistanceThreshold = (mode, kind) => {
  if (mode === 'linear') return kind === 'embedded' ? 700 : 540;
  return kind === 'embedded' ? 520 : 420;
};

const computeCandidateDistance = (referencePoint, candidate, mode) => {
  const dx = candidate.x - referencePoint.x;
  const dy = candidate.y - referencePoint.y;
  return mode === 'linear' ? Math.abs(dx) + Math.abs(dy) * 0.6 : Math.hypot(dx, dy);
};

const assignFeatureIdsToLabels = (svg, labelElements, featureGeometry, mode) => {
  const assignments = new Map();
  if (!featureGeometry || featureGeometry.all.length === 0) return assignments;

  const featureIds = new Set(featureGeometry.all.map((feature) => feature.id));
  const usedFeatureIds = new Set();
  const labelMeta = labelElements
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
    .filter(Boolean);

  labelMeta.forEach((meta) => {
    const existingId = String(meta.textEl.getAttribute('data-label-feature-id') || '').trim();
    if (!existingId || !featureIds.has(existingId)) return;
    if (!meta.candidateById.has(existingId)) return;
    if (usedFeatureIds.has(existingId)) return;
    assignments.set(meta.textEl, existingId);
    usedFeatureIds.add(existingId);
  });

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
  evaluateLabelRules
}) => {
  const {
    mode,
    generatedMode,
    results,
    selectedResultIndex,
    svgContainer,
    editableLabels,
    extractedFeatures,
    clickedFeature,
    labelTextScopeDialog,
    hiddenLabelTextDialog,
    labelOnDialog,
    featureOverrides,
    labelTextBulkOverrides,
    autoLabelReflowEnabled,
    labelReflowRequestSeq,
    labelReflowForceRequestSeq,
    labelReflowLastError,
    labelReflowProcessing
  } = state;

  // The label edits of the identity-keyed draft: Label visibility, label text,
  // and the label's source text (design Q4); Feature visibility stays.
  const LABEL_FIELDS = Object.freeze(['labelVisibility', 'labelText', 'labelSourceText']);
  const clearOverrides = () => {
    Object.values(featureOverrides).forEach((row) => {
      updateFeatureOverride(featureOverrides, row, Object.fromEntries(LABEL_FIELDS.map((field) => [field, null])));
    });
    Object.keys(labelTextBulkOverrides).forEach((key) => delete labelTextBulkOverrides[key]);
  };

  // A mounted label names its feature by rendered ID; the displayed Result's
  // catalog binds that ID to the feature's source identity (R3).
  const displayedFeatures = () => resultRenderedFeatures(state);
  const labelFeature = (featureIdRaw, displayed = displayedFeatures()) => {
    const featureId = String(featureIdRaw || '').trim();
    if (!featureId) return null;
    if (displayed) return displayed.get(featureId) || null;
    const key = normalizeKeyToken(featureId);
    return (extractedFeatures.value || []).find((feature) => normalizeKeyToken(feature?.svg_id) === key) || null;
  };
  const rowOf = (feature) => {
    const key = featureIdentityKeyOf(feature);
    return key ? featureOverrides[key] || null : null;
  };
  const labelRow = (featureId, displayed) => rowOf(labelFeature(featureId, displayed));

  // Whether the label rerender after this edit draws the feature: the last
  // Generate's request with the current edits and rules (F-3, Owner Q2). The
  // feature visibility owner's resolver answers as Generate does.
  const drawnContext = () => featureDrawnContext(state, { diagramOptions: getCommittedRequest()?.diagramOptions });
  const featureHidden = (feature, context = drawnContext()) => (
    featureDrawnInResult(feature, context, resultCatalogFeatures(state)) === false
  );
  // Owner decisions Q1 and Q2 (2026-10-04): Label visibility On takes effect
  // only when the diagram can draw the label. Returns 'hidden', 'underlay', or
  // 'embedded_only', or '' when the label is drawn.
  const labelOnBlocker = (feature, diagramOptions) => (
    featureHidden(feature) ? 'hidden' : labelDrawingBlocker(feature, diagramOptions)
  );
  // Why the displayed Result has no label for the feature: the blocker, unless
  // a scope or filter reason applies to Default label visibility first.
  // '' or a LABEL_ABSENCE_REASONS key.
  const labelAbsenceReason = (feature, diagramOptions) => {
    const blocker = labelOnBlocker(feature, diagramOptions);
    if (blocker === 'hidden' || blocker === 'underlay'
      || normalizeVisibilityMode(rowOf(feature)?.labelVisibility) === 'on') return blocker;
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

  const normalizeVisibilityMode = (value) => {
    const normalized = String(value || '').trim().toLowerCase();
    return normalized === 'on' || normalized === 'off' ? normalized : 'default';
  };

  const resetLabelsToSourceText = (svg) => {
    let changed = false;
    const displayed = displayedFeatures();
    svg.querySelectorAll(EDITABLE_LABEL_SELECTOR).forEach((textEl) => {
      const sourceText = labelRow(textEl.getAttribute('data-label-feature-id'), displayed)?.labelSourceText
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
  const projectLabelTextIntent = (svg) => {
    let changed = false;
    const displayed = displayedFeatures();
    svg.querySelectorAll(EDITABLE_LABEL_SELECTOR).forEach((textEl) => {
      const row = labelRow(textEl.getAttribute(LABEL_FEATURE_ID_ATTRIBUTE), displayed);
      const sourceText = row?.labelSourceText
        ?? textEl.getAttribute('data-label-source-text') ?? getLabelText(textEl);
      textEl.setAttribute('data-label-source-text', sourceText);
      const desiredText = row?.labelText ?? labelTextBulkOverrides[sourceText] ?? sourceText;
      if (getLabelText(textEl) === desiredText) return;
      setLabelText(textEl, desiredText);
      changed = true;
    });
    return changed;
  };

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

  const applyStoredVisibilityOverridesToSvg = (svg) => {
    let changed = false;
    let unavailableOverride = false;
    const displayed = displayedFeatures();
    const context = drawnContext();
    svg.querySelectorAll(VISIBILITY_LABEL_SELECTOR).forEach((textEl) => {
      const featureId = String(
        textEl.getAttribute(LABEL_FEATURE_ID_ATTRIBUTE) || ''
      ).trim();
      const feature = featureId ? labelFeature(featureId, displayed) : null;
      // A hidden feature hides its label through this projection too (F-3).
      const visibilityMode = !featureId
        ? 'default'
        : (featureHidden(feature, context) ? 'off' : rowOf(feature)?.labelVisibility);
      const projection = applyLabelVisibilityPreview(svg, textEl, visibilityMode);
      changed = projection.changed || changed;
      if (normalizeVisibilityMode(visibilityMode) !== 'default' && !projection.available) {
        unavailableOverride = true;
      }
    });
    return { changed, unavailableOverride };
  };

  const refreshEditableList = (svg) => {
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
  const featureById = (featureId) => (
    (clickedFeatureId() === featureId ? clickedFeature.value?.feat : null)
    || labelFeature(featureId)
  );

  const getEditableLabelByFeatureId = (featureId) => {
    const target = normalizeKeyToken(featureId);
    if (!target) return null;
    return (
      editableLabels.value.find(
        (entry) => normalizeKeyToken(entry?.featureId) === target
      ) || null
    );
  };

  const syncClickedFeatureLabelState = () => {
    if (!clickedFeature.value) return;
    const featureId = String(clickedFeature.value.svg_id || clickedFeature.value.id || '').trim();
    const entry = getEditableLabelByFeatureId(featureId);
    const row = rowOf(clickedFeature.value.feat);
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
        ? 'Generate the diagram again to edit this feature\'s label.'
        : 'No editable feature label for this feature in current diagram.');
  };

  /**
   * @param {{
   *   requiredFeatureIds?: readonly string[],
   *   optionalFeatureIds?: readonly string[],
   *   reportedLabelBinding?: { featureIds: readonly string[], report: (error: unknown) => void } | null,
   *   queueIncompleteVisibility?: boolean
   * }} [options]
   */
  const syncLabelEditor = ({
    requiredFeatureIds = [],
    optionalFeatureIds = [],
    reportedLabelBinding = null,
    queueIncompleteVisibility = true
  } = {}) => {
    // The retained Result can outlive its active mode. Keep its label intent
    // dormant until its own mode can project and reconcile those identities.
    if (generatedMode.value !== mode.value) return;
    if (!svgContainer.value) return;
    const svg = svgContainer.value.querySelector('svg');
    if (!svg) return;

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

    projectLabelIntent(svg, { queueIncompleteVisibility });
    // A label reflow keeps its Result: the same check reports to the reflow
    // after the binding completes instead of failing it.
    if (reportedLabelBinding) {
      try {
        requireUniqueEditableLabelBindings(labelElements, reportedLabelBinding.featureIds, { features });
      } catch (error) {
        reportedLabelBinding.report(error);
      }
    }
  };

  // One projection of the label intent onto the mounted Result, shared by a
  // live edit, a Label TSV import, History apply, and the display of a Result.
  const projectLabelIntent = (svg, { queueIncompleteVisibility = true } = {}) => {
    const textChanged = projectLabelTextIntent(svg);
    const visibilityProjection = applyStoredVisibilityOverridesToSvg(svg);
    refreshEditableList(svg);
    syncClickedFeatureLabelState();
    const changed = textChanged || visibilityProjection.changed;
    if (changed) commitLabelEdit();
    if (queueIncompleteVisibility && visibilityProjection.unavailableOverride) {
      queueLabelReflow(true);
    }
    return changed;
  };

  const reconcileLabelOverrides = () => {
    const svg = svgContainer.value?.querySelector?.('svg');
    return svg ? projectLabelIntent(svg) : false;
  };

  // A feature visibility edit shows or hides the feature's label in the same
  // action, then Auto Reflow places the labels as Generate does (F-3). An edit
  // that draws a feature the Result does not draw (`rerender`) always reruns
  // the rerender, which draws it and its label (R-5).
  const applyFeatureVisibilityToLabels = ({ reflow = true, rerender = false } = {}) => {
    if (generatedMode.value !== mode.value) return false;
    const svg = svgContainer.value?.querySelector?.('svg');
    if (!svg) return false;
    const projection = applyStoredVisibilityOverridesToSvg(svg);
    if (projection.changed) commitLabelEdit();
    if (reflow || rerender) queueLabelReflow(rerender || projection.unavailableOverride);
    return projection.changed;
  };

  // An edit whose Result Python must draw again asks for the automatic
  // rerender, Auto Reflow on or off: a specific-color rule edit or a History
  // step that changes a Legend source (Owner decision 2026-10-06, OV-43). A
  // request while one is pending joins it.
  const requestAutomaticRerender = () => {
    if (generatedMode.value !== mode.value || !svgContainer.value?.querySelector?.('svg')) return false;
    queueLabelReflow(true);
    return true;
  };

  const requestLabelTextChangeByKey = (labelKey, nextTextRaw) => {
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    if (!labelKey) return;
    const entry = editableLabels.value.find((candidate) => candidate.key === labelKey);
    if (!entry) return;
    const nextText = String(nextTextRaw ?? '');
    if (entry.text === nextText) return;

    labelTextScopeDialog.show = true;
    labelTextScopeDialog.labelKey = entry.key;
    labelTextScopeDialog.newText = nextText;
    labelTextScopeDialog.sourceText = entry.sourceText || entry.text;
    labelTextScopeDialog.featureId = entry.featureId || '';
    labelTextScopeDialog.matchingCount = editableLabels.value.filter(
      (candidate) => candidate.sourceText === (entry.sourceText || entry.text)
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

  const requestLabelTextChangeByFeatureId = (featureId, nextTextRaw) => {
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    const entry = getEditableLabelByFeatureId(featureId);
    if (!entry) return false;
    requestLabelTextChangeByKey(entry.key, nextTextRaw);
    return true;
  };

  // Stores one feature's Label visibility and shows it in its open popup.
  const setLabelVisibilityOverride = (feature, modeRaw) => {
    const nextMode = normalizeVisibilityMode(modeRaw);
    const previousMode = normalizeVisibilityMode(rowOf(feature)?.labelVisibility);
    updateFeatureOverride(featureOverrides, feature, { labelVisibility: nextMode === 'default' ? null : nextMode });
    if (featureIdentityKeyOf(clickedFeature.value?.feat) === featureIdentityKeyOf(feature)) {
      clickedFeature.value.labelVisibility = nextMode;
    }
    return previousMode !== nextMode;
  };

  // One feature's label text edit: the text equal to the label's own text
  // removes the edit; any other text is one line, kept with its source text.
  const applyDirectFeatureLabelOverride = (feature, labelTextRaw, sourceTextRaw, baselineTextRaw) => {
    if (!featureIdentityKeyOf(feature)) return false;
    const nextText = String(labelTextRaw ?? '');
    const sourceText = String(sourceTextRaw ?? '').trim();
    if (nextText === String(baselineTextRaw ?? '')) {
      return updateFeatureOverride(featureOverrides, feature, { labelText: null, labelSourceText: null });
    }
    return updateFeatureOverride(featureOverrides, feature, {
      labelText: normalizeFeatureOverrideLabelText(nextText),
      ...(sourceText ? { labelSourceText: sourceText } : {})
    });
  };

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

  // One popup Label edit: the overrides, then the displayed Result.
  const applyPopupLabelEdit = (feature, featureId, edit, { forceReflow = false } = {}) => {
    const visibilityChanged = setLabelVisibilityOverride(feature, edit.visibility);
    const textChanged = applyDirectFeatureLabelOverride(feature, edit.text, edit.sourceText, edit.sourceText);

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

    // Label text alone keeps Default visibility, which follows Show Labels and
    // the label filters. When they leave this feature unlabeled, ask whether
    // to show the label (On) or keep only the text.
    if (textChanged && !visibilityChanged && edit.text.trim()
      && edit.visibility === 'default'
      && !getEditableLabelByFeatureId(featureId)) {
      hiddenLabelTextDialog.featureId = featureId;
      hiddenLabelTextDialog.reason = labelAbsenceReason(feature, getCommittedRequest()?.diagramOptions);
      hiddenLabelTextDialog.show = true;
      return;
    }

    if (visibilityChanged || (!edit.hasEditableLabel && textChanged)) {
      queueLabelReflow(forceReflow || !visibilityProjection.available);
      return;
    }

    if (textChanged) {
      queueLabelReflow();
    }
  };

  const updateClickedFeatureLabelText = async () => {
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
    const apply = (options) => applyPopupLabelEdit(feature, featureId, edit, options);
    if (edit.visibility === 'on' && normalizeVisibilityMode(rowOf(feature)?.labelVisibility) !== 'on') {
      return applyLabelOn(featureId, apply);
    }
    return apply();
  };

  const closeHiddenLabelTextDialog = () => {
    hiddenLabelTextDialog.show = false;
    hiddenLabelTextDialog.featureId = '';
    hiddenLabelTextDialog.reason = '';
  };

  const handleHiddenLabelTextChoice = (choice) => {
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    if (!hiddenLabelTextDialog.show) return;
    const featureId = hiddenLabelTextDialog.featureId;
    closeHiddenLabelTextDialog();
    if (choice !== 'show' || !featureId) return;
    const feature = featureById(featureId);
    if (!featureIdentityKeyOf(feature)) return;
    return applyLabelOn(featureId, ({ forceReflow = false } = {}) => {
      setLabelVisibilityOverride(feature, 'on');
      const projection = applyDirectVisibilityToCurrentSvg(featureId, 'on');
      queueLabelReflow(forceReflow || !projection.available);
    });
  };

  // Owner decisions Q1 and Q2 (2026-10-04): applying Label visibility On asks
  // first when the diagram cannot draw the label: a hidden feature (Show
  // feature and label, or Keep feature hidden), a feature drawn as underlay,
  // and, after the label reflow, a label that does not fit with Embedded Only
  // (Keep without label). The caller's History step stays open until the
  // choice, so a choice is one step and Cancel, which leaves the label intent
  // as it was, records none. Global settings never change.
  const applyLabelOn = async (featureId, apply) => {
    // The blocker reads Python's rule matches; prepared ones answer at once.
    const prepared = /** @type {any} */ (prepareDrawnFeatureMatches());
    if (typeof prepared?.then === 'function') await prepared;
    const feature = featureById(featureId);
    const diagramOptions = getCommittedRequest()?.diagramOptions;
    const blocker = labelOnBlocker(feature, diagramOptions);
    if (!blocker) return apply();
    return confirmLabelOn({ featureId, feature, diagramOptions, blocker, apply });
  };

  /** @type {((choice: string) => void) | null} */
  let answerLabelOnDialog = null;
  const askLabelOn = (reason, feature) => {
    answerLabelOnDialog?.('cancel');
    labelOnDialog.reason = reason;
    labelOnDialog.featureType = String(feature?.type || '');
    labelOnDialog.show = true;
    return new Promise((resolve) => { answerLabelOnDialog = resolve; });
  };

  // The Label visibility On dialog answers the Apply that opened it.
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

  const captureFeatureLabelIntent = (feature) => {
    const row = rowOf(feature);
    return Object.fromEntries(LABEL_FIELDS.map((field) => [field, row?.[field] ?? null]));
  };
  const restoreFeatureLabelIntent = (feature, captured) => updateFeatureOverride(featureOverrides, feature, captured);

  // Whether the label reflow that this Apply queued drew the feature's label;
  // null when no reflow replaced the Result (it was skipped or failed, or a
  // Generate took over), so nothing is asked.
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

  const confirmLabelOn = async ({ featureId, feature, diagramOptions, blocker, apply }) => {
    const showFeature = blocker === 'hidden';
    if (showFeature) {
      const choice = await askLabelOn('hidden', feature);
      if (choice === 'keep') return apply();
      if (choice !== 'show') return cancelLabelOn();
    }
    // A feature shown here can still be drawn as underlay or with Embedded Only.
    const drawing = showFeature ? labelDrawingBlocker(feature, diagramOptions) : blocker;
    if (drawing === 'underlay' && await askLabelOn('underlay', feature) !== 'keep') return cancelLabelOn();
    const embeddedOnly = drawing === 'embedded_only';
    const labelIntent = captureFeatureLabelIntent(feature);
    const featureVisibility = getFeatureVisibilityOverride(featureOverrides, feature);
    if (showFeature) setFeatureVisibility?.(feature, 'on', { triggerReflow: false });
    apply({ forceReflow: embeddedOnly });
    if (!embeddedOnly) return true;
    const applied = JSON.stringify(captureFeatureLabelIntent(feature));
    const drawn = await reflowDrawsLabel(featureId);
    // A later edit of this label decides instead of this Apply.
    if (drawn !== false || JSON.stringify(captureFeatureLabelIntent(feature)) !== applied) return true;
    if (await askLabelOn('embedded_only', feature) === 'keep') return true;
    // Cancel restores the label intent and feature visibility this Apply
    // changed, and the label reflow draws the Result they describe again.
    restoreFeatureLabelIntent(feature, labelIntent);
    if (showFeature) setFeatureVisibility?.(feature, featureVisibility, { triggerReflow: false });
    reconcileLabelOverrides();
    queueLabelReflow(true);
    return cancelLabelOn();
  };

  // The popup note when the displayed Result draws no label for the clicked
  // feature: why the diagram cannot draw it (Owner Q1, Q2), or that On shows it.
  const clickedFeatureLabelHint = computed(() => {
    const clicked = clickedFeature.value;
    if (!clicked?.hasEditableLabel) return '';
    const draft = normalizeVisibilityMode(clicked.labelVisibility);
    if (draft === 'off') return '';
    void results.value; // The committed request changes with the displayed Results.
    void rulePreparation?.pending?.value; // Python's rule matches arrive when an evaluation ends.
    const featureId = clickedFeatureId();
    const feature = featureById(featureId);
    const diagramOptions = getCommittedRequest()?.diagramOptions;
    const blocker = labelOnBlocker(feature, diagramOptions);
    // A hidden feature's label is hidden with it (F-3).
    if (clicked.labelKey && blocker !== 'hidden') return '';
    let next = '';
    if (normalizeVisibilityMode(rowOf(clicked.feat)?.labelVisibility) === 'on') {
      next = ` Its Label visibility "On" applies when ${blocker === 'hidden' ? 'the feature is shown' : 'the label can be drawn'}.`;
    } else if (!blocker && draft === 'default') {
      next = ' Choose On to show it.';
    } else if (!blocker) {
      return '';
    }
    const reason = labelAbsenceReason(feature, diagramOptions);
    return `${LABEL_ABSENT_PREFIX}${LABEL_ABSENCE_REASONS[reason] || ''}${next}`;
  });

  // The Label Not Shown dialog: the same reason sentence as the popup note.
  const hiddenLabelTextMessage = computed(() => (
    `${LABEL_ABSENT_PREFIX}${LABEL_ABSENCE_REASONS[hiddenLabelTextDialog.reason] || LABEL_ABSENCE_UNKNOWN}`
    + ' The edited text will not appear unless you show this label.'
  ));

  // A label list edit: one line of text, or an empty text, which hides the
  // label (Label visibility Off) as an empty label-table text does.
  const applyLabelTextEdit = (feature, text, sourceText) => {
    const labelText = normalizeFeatureOverrideLabelText(String(text ?? ''));
    return updateFeatureOverride(featureOverrides, feature, labelText
      ? { labelText, labelSourceText: sourceText || null }
      : { labelText: null, labelVisibility: 'off', labelSourceText: sourceText || null });
  };

  const handleLabelTextScopeChoice = (choice) => {
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
      const displayed = displayedFeatures();
      svg.querySelectorAll(EDITABLE_LABEL_SELECTOR).forEach((textEl) => {
        const candidateSource = textEl.getAttribute('data-label-source-text') || '';
        if (candidateSource !== sourceText) return;
        matchedCount += 1;
        const candidate = labelFeature(textEl.getAttribute('data-label-feature-id'), displayed);
        if (featureIdentityKeyOf(candidate)) {
          applyLabelTextEdit(candidate, newText, candidateSource);
        } else {
          hasUntrackableMatch = true;
        }
        setLabelText(textEl, newText);
      });
      if (hasUntrackableMatch || matchedCount === 0) {
        labelTextBulkOverrides[sourceText] = newText;
      } else if (Object.prototype.hasOwnProperty.call(labelTextBulkOverrides, sourceText)) {
        delete labelTextBulkOverrides[sourceText];
      }
    } else if (choice === 'single' && targetKey) {
      const targetEl = svg.querySelector(`text[data-label-key="${CSS.escape(targetKey)}"]`);
      if (targetEl) {
        setLabelText(targetEl, newText);
      }
      const feature = labelFeature(featureId);
      if (featureIdentityKeyOf(feature)) applyLabelTextEdit(feature, newText, sourceText);
    } else {
      closeLabelTextScopeDialog();
      return;
    }

    commitLabelEdit();
    closeLabelTextScopeDialog();
    syncLabelEditor();
    queueLabelReflow();
  };

  const resetAllLabelTextOverrides = () => {
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    const svg = svgContainer.value?.querySelector('svg');
    if (svg) resetLabelsToSourceText(svg);
    clearOverrides();
    if (!svgContainer.value) {
      editableLabels.value = [];
      closeLabelTextScopeDialog();
      closeHiddenLabelTextDialog();
      return;
    }
    if (!svg) return;
    applyStoredVisibilityOverridesToSvg(svg);

    closeLabelTextScopeDialog();
    closeHiddenLabelTextDialog();
    commitLabelEdit();
    syncLabelEditor();
    queueLabelReflow();
  };

  // A Label TSV row selects labels in every Result of a batch. The displayed
  // Result binds its labels as the editor does; the committed content of
  // another Result carries the renderer's label bindings (B6).
  const labelSourceText = (element, feature) => (
    rowOf(feature)?.labelSourceText
    ?? element.getAttribute('data-label-source-text') ?? getLabelText(element)
  );
  const labelImportTargets = (svg) => {
    const elements = collectEditableLabelElements(svg, mode.value);
    const assignments = assignFeatureIdsToLabels(svg, elements, collectFeatureGeometry(svg), mode.value);
    const displayed = displayedFeatures();
    const targets = elements.map((element) => {
      const featureId = assignments.get(element) || '';
      const feature = featureId ? labelFeature(featureId, displayed) : null;
      return { featureId, feature, sourceText: labelSourceText(element, feature) };
    });
    const seen = new Set(targets.map(({ feature }) => featureIdentityKeyOf(feature)).filter(Boolean));
    results.value.forEach((result, index) => {
      if (index === selectedResultIndex.value || typeof result?.content !== 'string') return;
      const features = resultRenderedFeatures(state, index);
      const root = new DOMParser().parseFromString(result.content, 'image/svg+xml').documentElement;
      root.querySelectorAll(`text[${LABEL_FEATURE_ID_ATTRIBUTE}]`).forEach((element) => {
        const featureId = String(element.getAttribute(LABEL_FEATURE_ID_ATTRIBUTE) || '').trim();
        const feature = featureId ? features?.get(featureId) || null : null;
        const key = featureIdentityKeyOf(feature);
        if (!key || seen.has(key)) return;
        seen.add(key);
        targets.push({ featureId, feature, sourceText: labelSourceText(element, feature) });
      });
    });
    return targets;
  };

  let labelImportRevision = 0;
  const labelImportFailure = ref(null);
  const canRetryLabelImportFailure = computed(() => Boolean(labelImportFailure.value
    && state.errorLog?.value === labelImportFailure.value.error
    && labelImportFailure.value.revision === labelImportRevision
    && rulePreparation.isCurrent(labelImportFailure.value.snapshot)
    && labelImportFailure.value.intent === labelIntentSignature()));
  const retryLabelImportFailure = () => canRetryLabelImportFailure.value ? labelImportFailure.value.retry() : false;
  const editLabelImportFailure = () => {
    if (canRetryLabelImportFailure.value) labelImportFailure.value.input?.click?.();
  };
  const labelIntentSignature = () => JSON.stringify([featureOverrides, labelTextBulkOverrides]);
  const loadLabelOverrideTable = async (event) => {
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    const input = event?.target;
    const file = input?.files?.[0];
    if (!file) return;
    const sourceSvg = svgContainer.value?.querySelector('svg') || null;
    const before = rulePreparation.snapshot();
    const labelIntent = labelIntentSignature();
    const revision = ++labelImportRevision;
    const previousError = state.errorLog?.value;

    try {
      const text = await readFileText(file);
      if (input.files?.[0] !== file || (svgContainer.value?.querySelector('svg') || null) !== sourceSvg) return;
      if (state.sessionOperationAvailability?.()) return state.sessionOperationAvailability();
      const rows = parseLabelOverrideTsv(text);

      const svg = sourceSvg;
      const labels = svg ? labelImportTargets(svg) : [];
      const evaluation = await evaluateLabelRules({
        kind: 'label', rules: rows,
        features: labels.map((entry) => ruleFeaturePayload(
          entry.feature || { type: '', svg_id: entry.featureId }, entry.sourceText
        ))
      });
      if (state.sessionOperationAvailability?.()) return state.sessionOperationAvailability();
      if (revision !== labelImportRevision || labelIntent !== labelIntentSignature()
        || input.files?.[0] !== file || !rulePreparation.isCurrent(before)
        || (svgContainer.value?.querySelector('svg') || null) !== sourceSvg) return;
      if (!svg) {
        window.alert(`Loaded ${rows.length} row(s). No diagram is currently displayed.`);
        return;
      }
      // A `location` or `record_location` row cannot be decided for a label
      // whose drawn values the Result does not record (feature catalog 3 or 4
      // on a cropped, reverse-complemented, or rotated record). The import is
      // declined whole, so the label edits it would replace stay (R4).
      const undecided = labels.filter(({ feature }) => drawnSelectorUnknown(feature) && rows.some((row) => (
        DRAWN_SELECTOR_QUALIFIERS.has(row.qualifier.toLowerCase())
        && (row.recordId === '*' || row.recordId === String(feature.record_id ?? ''))
        && (row.featureType === '*' || row.featureType === String(feature.type ?? ''))
      ))).length;
      if (undecided > 0) {
        window.alert(`Loaded ${rows.length} row(s). Not applied: location and record_location rows cannot be `
          + `matched to ${undecided} label(s) on cropped, reverse-complemented, or rotated records until the `
          + 'diagram records where it drew them. Generate the diagram again, then load the table.');
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
        window.alert(`Loaded ${rows.length} row(s). Not applied: ${skippedNonTrackableCount > 0
          ? `the ${skippedNonTrackableCount} matched label(s) lacked a feature key`
          : 'no row matched a label of the diagram'}. The existing label edits were kept.`);
        return;
      }
      clearOverrides();
      applicable.forEach(({ entry, matchedRow, tracked }) => {
        if (matchedRow.isGlobalLabelRule) {
          labelTextBulkOverrides[entry.sourceText] = String(matchedRow.labelText ?? '');
          if (tracked) updateFeatureOverride(featureOverrides, entry.feature, { labelSourceText: entry.sourceText || null });
        } else {
          applyLabelTextEdit(entry.feature, matchedRow.labelText, entry.sourceText);
        }
      });

      closeLabelTextScopeDialog();
      closeHiddenLabelTextDialog();
      syncLabelEditor();
      queueLabelReflow();

      let message = `Loaded ${rows.length} row(s). Applied to ${applicable.length} label(s).`;
      if (skippedNonTrackableCount > 0) {
        message += ` ${skippedNonTrackableCount} match(es) lacked a feature key and were not applied.`;
      }
      if (state.errorLog?.value === labelImportFailure.value?.error) state.errorLog.value = null;
      labelImportFailure.value = null;
      window.alert(message);
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
    const savedTable = serializeLabelOverrideRows(state.canonicalLabelOverrideRows?.value);
    const rows = savedTable
      ? savedTable.trimEnd().split('\n')
      : Object.keys(labelTextBulkOverrides).sort((a, b) => a.localeCompare(b)).filter(Boolean)
        .map((sourceText) => `*\t*\tlabel\t^${escapeRegexLiteral(sourceText)}$\t${
          normalizeTsvCell(labelTextBulkOverrides[sourceText])}`);
    if (rows.length === 0) {
      window.alert('No label rules to export. Export Feature Edits TSV writes per-feature label edits.');
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
    resetAllLabelTextOverrides,
    syncClickedFeatureLabelState,
    syncLabelEditor,
    updateClickedFeatureLabelText
  };
};
