// @ts-check
import { resolveColorToHex } from '../../utils/color-utils.js';
import { reportRuleRunFailure, runWhenPrepared } from '../rule-matching.js';
import {
  formatFeatureLength,
  formatFeatureLocation,
  getFeatureCaption,
  normalizeStringArray,
  resolveFeatureProteinId
} from '../../services/feature-utils.js';
import {
  PAIRWISE_MATCH_SELECTOR,
  buildPairwiseMatchHoverSummary,
  buildMatchPopupPayload
} from '../pairwise-match-popup.js';
import { buildFeatureSequenceFastas } from '../../services/feature-sequence-fasta.js';
import { getFeatureOverride } from '../../services/feature-override-identity.js';
import { featureIdentityKeyOf, featureOverrideValue } from '../../services/feature-placement.js';
import { resultCatalogFeatures, stableFeatureOverrideKey } from '../../services/feature-catalog.js';
import { COMPARISON_LEGEND_SELECTOR } from '../legend/utils.js';
import { recordStructuralMetric } from '../../services/runtime-test-hooks.js';
import {
  featureIdentity,
  identityMatches,
  renderedFeatureIdentity
} from '../../services/feature-identity.js';
import {
  FEATURE_ID_ATTRIBUTE,
  FEATURE_SELECTOR,
  buildFeatureElementIndex,
  clearFeatureElementIndex,
  getFeatureElementIndex,
  getFeatureElements,
  getFeatureFillElements,
  getFeatureIdentity,
  normalizeFeatureIdentity
} from '../../services/feature-dom.js';

export {
  FEATURE_ID_ATTRIBUTE,
  FEATURE_SELECTOR,
  buildFeatureElementIndex,
  clearFeatureElementIndex,
  getFeatureElementIndex,
  getFeatureElements,
  getFeatureFillElements,
  getFeatureIdentity,
  normalizeFeatureIdentity
};

/**
 * The feature selection owner's reactions to a click or a drag on the preview.
 * @typedef {object} FeatureSelectionPort
 * @property {() => boolean} [consumeSuppressNextClick]
 * @property {(event: Event, svg: Element) => Element | null} [getSelectableFeatureTarget]
 * @property {(selectionId: string) => any} [toggleFeatureSelection]
 * @property {(selectionId: string, options?: { additive?: boolean }) => any} [selectFeatureRange]
 * @property {(options?: { clearStatus?: boolean }) => any} [clearFeatureSelection]
 * @property {(svgId: string) => any} [markPlainFeatureClick]
 * @property {(event: Event, svg: Element) => boolean} [startMarqueePointer]
 * @property {(event: Event) => boolean} [moveMarqueePointer]
 * @property {(event: Event) => boolean} [commitMarqueePointer]
 * @property {() => any} [cancelMarqueePointer]
 */

/**
 * The preview owner's pan and zoom interaction, which pauses feature hover.
 * @typedef {object} PreviewTransformInteractionPort
 * @property {() => boolean} [isActive]
 * @property {(listener: (change: { active: boolean, kind: string, event: Event, reconcile: boolean }) => void) => (() => void) | void} [subscribe]
 */

/**
 * The visibility owner's preparation of the rule matches `resolveFeatureDrawn`
 * reads. Strict, it rejects when the color preparation fails and resolves to
 * false when it is stale.
 * @typedef {(options?: { strict?: boolean }) => boolean | Promise<boolean | { error: any }>} PrepareDrawnFeatureMatchesPort
 */

/**
 * @typedef {object} FeatureSvgActionsOptions
 * @property {Record<string, any>} state App state (state.js; not yet typed).
 * @property {(feature: Record<string, any>) => string} getFeatureColor The rule owner's color of a feature.
 * @property {(feature: Record<string, any>) => string} [getEffectiveLegendCaption] The rule owner's Legend caption of a feature.
 * @property {(() => void) | null} [onFeaturePopupOpened] The label owner's sync of its editor when the popup opens.
 * @property {PrepareDrawnFeatureMatchesPort} [prepareDrawnFeatureMatches]
 *   The visibility owner's preparation of the rule matches `resolveFeatureDrawn` reads.
 * @property {FeatureSelectionPort | null} [featureSelection]
 * @property {((changes: { featureId: string, mode: string }[], options?: { reason?: string }) => boolean) | null} [applyFeatureVisibilityChanges]
 *   The preview owner's projection of feature visibility changes.
 * @property {PreviewTransformInteractionPort | null} [previewTransformInteraction]
 */

/**
 * The state of the feature and match hover summary. The tooltip element is
 * created on the first hover; `lastEvent` holds the last pointer position seen.
 * @typedef {object} HoverSummaryState
 * @property {HTMLElement | null} element
 * @property {number | null} timer
 * @property {number | null} frame
 * @property {boolean} visible
 * @property {string} activeSvgId
 * @property {{ clientX?: number, clientY?: number } | null} lastEvent
 */

/**
 * The Similarity alignment overlay over the preview. `observer` is set once the
 * overlay is built and stays null without ResizeObserver.
 * @typedef {object} AlignmentOverlay
 * @property {HTMLElement} element
 * @property {HTMLElement} scrollRoot
 * @property {number} frame
 * @property {ResizeObserver | null} observer
 * @property {() => void} schedule
 * @property {(recordKey: string, anchor: any) => void} onSelect
 * @property {(recordKey: string, candidateKey: string) => void} onHover
 */

/**
 * @typedef {object} AlignmentOverlayRequest
 * @property {any} reference
 * @property {any[]} [ambiguities]
 * @property {(recordKey: string, anchor: any) => void} onSelect
 * @property {(recordKey: string, candidateKey: string) => void} onHover
 */

/**
 * The delegated handlers bound to one Result SVG. The lookup maps are built on
 * first use, and the methods are set while the handlers are attached.
 * @typedef {object} DelegatedFeatureHandlers
 * @property {Element} svg
 * @property {Map<string, Element[]> | null} pathsByIdMap
 * @property {Map<string, Record<string, any>> | null} featureLookup
 * @property {Map<string, Set<string>> | null} featureIdsByOrthogroupId
 * @property {Map<string, Element[]> | null} comparisonElementsByOrthogroupId
 * @property {Map<string, Element[]> | null} comparisonElementsByCollinearityBlockId
 * @property {Map<string, Element[]> | null} comparisonElementsByMatchId
 * @property {boolean} pairwiseAffordancesPrepared
 * @property {string | null} activeHoverSvgId
 * @property {string} activeHoverKey
 * @property {Element | null} activeMatchHoverElement
 * @property {string} activeMatchHoverKey
 * @property {Element | null} pendingMatchElement
 * @property {Set<Element>} activeHoverElements
 * @property {{ clientX: number, clientY: number } | null} transformPointer
 * @property {boolean} reconcileHoverAfterTransform
 * @property {number | null} hoverReconcileFrame
 * @property {string} alignmentCandidateSvgId
 * @property {{ featureSvgId: string | null, matchElement: Element | null } | null} alignmentCandidateRestore
 * @property {AlignmentOverlay | null} alignmentOverlay
 * @property {Map<string, any>} alignmentCandidatesBySvgId
 * @property {((request: AlignmentOverlayRequest) => void) | null} showAlignmentOverlay
 * @property {(() => void) | null} clearAlignmentOverlay
 * @property {((anchor: any) => boolean) | null} previewAlignmentCandidate
 * @property {((options?: { restore?: boolean }) => void) | null} clearAlignmentCandidatePreview
 * @property {((eventLike: any, kind?: string) => void) | null} beginPreviewTransformInteraction
 * @property {((options?: { reconcile?: boolean }) => void) | null} endPreviewTransformInteraction
 * @property {(() => void) | null} cleanup
 */

/** @param {FeatureSvgActionsOptions} options */
export const createFeatureSvgActions = ({
  state,
  getFeatureColor,
  getEffectiveLegendCaption,
  onFeaturePopupOpened = null,
  // R13: the popup opens once the rule matches `resolveFeatureDrawn` reads are
  // prepared; the composition root injects the visibility owner's preparation.
  prepareDrawnFeatureMatches = () => true,
  featureSelection = null,
  // R13: the preview owner's projection of feature visibility changes.
  applyFeatureVisibilityChanges = null,
  previewTransformInteraction = null
}) => {
  const {
    orthogroups,
    collinearGroups,
    orthogroupNameOverrides,
    orthogroupDescriptionOverrides,
    extractedFeatures,
    biologicalFeatures,
    featuresBySvgId,
    featureColorOverrides,
    featureOverrides,
    svgContainer,
    clickedFeature,
    clickedFeaturePos,
    clickedPairwiseMatch,
    clickedPairwiseMatchPos,
    matchSequenceRegistry,
    selectedAnnotation,
    featurePopupSize,
    featureSelectionDrag,
    adv
  } = state;
  /** @type {DelegatedFeatureHandlers | null} */
  let delegatedFeatureHandlers = null;
  const isPreviewTransformInteractionActive = () => Boolean(
    previewTransformInteraction?.isActive?.()
  );
  /** @type {HoverSummaryState} */
  const hoverSummaryState = {
    element: null,
    timer: null,
    frame: null,
    visible: false,
    activeSvgId: '',
    lastEvent: null
  };

  const getOrthogroupIds = (value) =>
    Array.from(new Set(
      String(value || '')
        .split(';')
        .map((entry) => entry.trim())
        .filter(Boolean)
    ));

  const renderedFeatureSvgId = (feature) => String(
    feature?.rendered_svg_id ||
    feature?.renderedSvgId ||
    feature?.rendered_feature_svg_id ||
    feature?.renderedFeatureSvgId ||
    feature?.svg_id ||
    ''
  ).trim();

  const normalizeVisibilityMode = (value) => {
    const normalized = String(value || '').trim().toLowerCase();
    if (normalized === 'suppress') return 'exclude_matching';
    return ['on', 'off', 'exclude_matching'].includes(normalized) ? normalized : 'default';
  };

  const getPopupPosition = (eventLike, popupWidth = 720, popupHeight = 520) => {
    const margin = 12;
    const fallbackX = window.innerWidth / 2;
    const fallbackY = window.innerHeight / 2;
    const resolvedPopupWidth = Math.min(popupWidth, Math.max(0, window.innerWidth - (2 * margin)));
    const resolvedPopupHeight = Math.min(popupHeight, Math.max(0, window.innerHeight - (2 * margin)));
    const rawX = Number.isFinite(eventLike?.clientX) ? eventLike.clientX + 10 : fallbackX;
    const rawY = Number.isFinite(eventLike?.clientY) ? eventLike.clientY + 10 : fallbackY;
    const maxX = Math.max(margin, window.innerWidth - resolvedPopupWidth - margin);
    const maxY = Math.max(margin, window.innerHeight - resolvedPopupHeight - margin);
    return {
      x: Math.min(Math.max(rawX, margin), maxX),
      y: Math.min(Math.max(rawY, margin), maxY)
    };
  };

  const normalizeQualifierRows = (qualifiers) => {
    if (!qualifiers || typeof qualifiers !== 'object' || Array.isArray(qualifiers)) return [];
    return Object.entries(qualifiers)
      .map(([key, value]) => {
        const values = normalizeStringArray(value);
        return {
          key: String(key || ''),
          values,
          copyText: values.join('\n'),
          displayValue: values.join('\n')
        };
      })
      .filter((row) => row.key && row.values.length > 0)
      .sort((left, right) => left.key.localeCompare(right.key));
  };

  const getQualifierFirstValue = (feat, key) => {
    const normalizedKey = String(key || '').trim().toLowerCase();
    if (!normalizedKey) return '';
    const directValue = feat?.[normalizedKey];
    const qualifierValue = feat?.qualifiers && typeof feat.qualifiers === 'object'
      ? feat.qualifiers[normalizedKey]
      : null;
    const values = normalizeStringArray(directValue || qualifierValue)
      .map((value) => String(value || '').trim())
      .filter(Boolean);
    return values[0] || '';
  };

  const getHoverSummaryPrimaryLabel = (feat) => (
    getQualifierFirstValue(feat, 'gene') ||
    getQualifierFirstValue(feat, 'locus_tag') ||
    getQualifierFirstValue(feat, 'product') ||
    getFeatureCaption(feat) ||
    ''
  );

  const createHoverSummaryElement = (tagName, className = '', text = '') => {
    const element = document.createElement(tagName);
    if (className) element.className = className;
    if (text !== '') element.textContent = text;
    return element;
  };

  const addHoverSummaryRow = (container, label, value, { clamp = false } = {}) => {
    const normalizedValue = String(value === null || value === undefined ? '' : value).trim();
    if (!normalizedValue) return;
    const row = createHoverSummaryElement('div', 'feature-hover-summary-row');
    row.appendChild(createHoverSummaryElement('div', 'feature-hover-summary-label', label));
    row.appendChild(createHoverSummaryElement(
      'div',
      `feature-hover-summary-value${clamp ? ' is-clamped' : ''}`,
      normalizedValue
    ));
    container.appendChild(row);
  };

  const buildHoverSummaryRows = (feat, primaryLabel) => {
    const product = getQualifierFirstValue(feat, 'product');
    const gene = getQualifierFirstValue(feat, 'gene');
    const locusTag = getQualifierFirstValue(feat, 'locus_tag');
    const note = getQualifierFirstValue(feat, 'note');
    const locationText = formatFeatureLocation(feat);
    const effectiveCaption = String(getEffectiveLegendCaption?.(feat) || '').trim();
    const rows = [];

    if (gene && gene !== primaryLabel) rows.push(['Gene', gene]);
    if (locusTag && locusTag !== primaryLabel) rows.push(['Locus', locusTag]);
    if (product && product !== primaryLabel) rows.push(['Product', product, true]);
    if (note && note !== primaryLabel && note !== product) rows.push(['Note', note, true]);
    rows.push(['Length', formatFeatureLength(feat)]);
    rows.push(['Location', locationText]);
    rows.push(['Record', feat?.record_id || '']);
    if (feat?.orthogroupId) rows.push(['Similarity group', feat.orthogroupId]);
    if (effectiveCaption && effectiveCaption !== primaryLabel) rows.push(['Legend', effectiveCaption]);
    return rows;
  };

  const buildOrthogroupDetailRows = (feat) => {
    const member = feat?.orthogroupMember || feat?.orthogroup_member || null;
    const proteinId = resolveFeatureProteinId(feat, member);
    const rows = [
      { key: 'orthogroup_id', label: 'Similarity group ID', value: feat?.orthogroupId || feat?.orthogroup_id },
      { key: 'orthogroup_members', label: 'Members', value: feat?.orthogroupMemberCount || feat?.orthogroup_member_count },
      { key: 'orthogroup_coverage', label: 'Record coverage', value: feat?.orthogroupRecordCoverage || feat?.orthogroup_record_coverage },
      { key: 'protein_id', label: 'Protein ID', value: proteinId }
    ];
    return rows.filter((row) => String(row.value === null || row.value === undefined ? '' : row.value) !== '');
  };

  const buildDetailRows = ({ defaultLabel, feat, locationText }) => {
    const rows = [
      { key: 'label', label: 'Label', value: defaultLabel },
      { key: 'record_id', label: 'Record ID', value: feat.record_id },
      { key: 'type', label: 'Feature type', value: feat.type },
      { key: 'location', label: 'Location', value: locationText }
    ];
    rows.push(...buildOrthogroupDetailRows(feat));
    return rows
      .map((row) => ({ ...row, value: row.value === null || row.value === undefined ? '' : String(row.value) }))
      .filter((row) => row.value !== '');
  };

  /** @returns {Map<string, Record<string, any>>} */
  const buildFeatureLookup = () => {
    if (featuresBySvgId?.value instanceof Map) return featuresBySvgId.value;
    const indexed = new Map();
    const features = Array.isArray(extractedFeatures.value) ? extractedFeatures.value : [];
    for (const feat of features) {
      const svgId = renderedFeatureSvgId(feat);
      if (!svgId || indexed.has(svgId)) continue;
      indexed.set(svgId, feat);
    }
    return indexed;
  };

  const getFeatureTarget = (target, svg) => {
    if (!target || typeof target.closest !== 'function') return null;
    const matchEl = target.closest(PAIRWISE_MATCH_SELECTOR);
    if (matchEl && svg.contains(matchEl)) return null;
    const featureEl = target.closest(FEATURE_SELECTOR);
    if (!featureEl || !svg.contains(featureEl)) return null;
    return featureEl;
  };

  const getPairwiseMatchTarget = (target, svg) => {
    if (!target || typeof target.closest !== 'function') return null;
    const matchEl = target.closest(PAIRWISE_MATCH_SELECTOR);
    if (!matchEl || !svg.contains(matchEl)) return null;
    return matchEl;
  };

  const getTopmostSvgTarget = (eventLike, svg, selector) => {
    if (!svg || !selector || !Number.isFinite(eventLike?.clientX) || !Number.isFinite(eventLike?.clientY)) {
      return null;
    }
    const stack = typeof document.elementsFromPoint === 'function'
      ? document.elementsFromPoint(eventLike.clientX, eventLike.clientY)
      : [];
    for (const element of stack) {
      if (!element || !svg.contains(element)) continue;
      const target = element.matches?.(selector)
        ? element
        : element.closest?.(selector);
      if (target && svg.contains(target)) return target;
    }
    return null;
  };

  const getFeatureClickTarget = (eventLike, svg) =>
    getTopmostSvgTarget(eventLike, svg, FEATURE_SELECTOR) || getFeatureTarget(eventLike?.target, svg);

  const getPairwiseMatchClickTarget = (eventLike, svg) =>
    getTopmostSvgTarget(eventLike, svg, PAIRWISE_MATCH_SELECTOR) || getPairwiseMatchTarget(eventLike?.target, svg);

  const isBackgroundPreviewClick = (target, svg) => {
    if (!target || !svg?.contains?.(target)) return false;
    if (target === svg) return true;
    if (
      target.closest?.(
        [
          FEATURE_SELECTOR,
          PAIRWISE_MATCH_SELECTOR,
          'text[data-label-editable="true"]',
          '[data-label-key]',
          '[data-label-feature-id]',
          '#legend',
          '#feature_legend',
          COMPARISON_LEGEND_SELECTOR,
          '#horizontal_legend',
          '#vertical_legend',
          '[data-legend-key]'
        ].join(', ')
      )
    ) {
      return false;
    }
    if (state.layoutRepositionMode?.value && target.closest?.('g[id]')) return false;
    return true;
  };

  const cleanupDelegatedFeatureHandlers = () => {
    if (!delegatedFeatureHandlers?.cleanup) return;
    delegatedFeatureHandlers.cleanup();
    delegatedFeatureHandlers = null;
  };

  // `renderedSvgId` is '' for a feature the displayed Result does not draw.
  /**
   * @param {Record<string, any>} feat
   * @param {Element | null} [featureElement]
   * @param {string} [renderedSvgId]
   */
  const buildClickedFeaturePayload = (feat, featureElement = null, renderedSvgId = undefined) => {
    const defaultLabel = getFeatureCaption(feat);
    const existingOverride = getFeatureOverride(featureColorOverrides, feat);
    const effectiveCaption = String(getEffectiveLegendCaption?.(feat) || existingOverride?.caption || defaultLabel || '').trim();
    const locationText = formatFeatureLocation(feat);
    const locationParts = Array.isArray(feat.location_parts) ? feat.location_parts : [];
    const qualifierRows = normalizeQualifierRows(feat.qualifiers);
    const sequenceWarnings = normalizeStringArray(feat.sequence_warnings);
    const nucleotideSequence = String(feat.nucleotide_sequence || '');
    const aminoAcidSequence = String(feat.amino_acid_sequence || '');
    const { nucleotideFasta, aminoAcidFasta } = buildFeatureSequenceFastas(feat, {
      nucleotideSequence,
      aminoAcidSequence
    });

    const currentColor = resolveColorToHex(
      featureElement?.getAttribute('fill') || getFeatureColor(feat)
    );
    const currentStrokeColor = featureElement?.getAttribute('stroke') || '#000000';
    const currentStrokeWidth = parseFloat(featureElement?.getAttribute('stroke-width') ?? '') || 0.5;

    const actualSvgId = String(renderedSvgId ?? renderedFeatureSvgId(feat)).trim();
    const visibilityMode = normalizeVisibilityMode(featureOverrideValue(featureOverrides, feat, 'featureVisibility'));

    return {
      id: feat.id,
      svg_id: actualSvgId,
      rendered_svg_id: actualSvgId,
      stable_svg_id: feat.stable_svg_id || feat.stable_feature_id || feat.svg_id || '',
      label: defaultLabel,
      location: locationText,
      locationParts,
      color: currentColor,
      feat,
      activeTab: 'edit',
      recordId: String(feat.record_id || ''),
      gene: getQualifierFirstValue(feat, 'gene'),
      recordIdx: Number.isInteger(Number(feat.record_idx)) ? Number(feat.record_idx) : null,
      featureType: String(feat.type || ''),
      start: Number.isFinite(Number(feat.start)) ? Number(feat.start) : null,
      end: Number.isFinite(Number(feat.end)) ? Number(feat.end) : null,
      strand: String(feat.strand || ''),
      qualifiers: feat.qualifiers && typeof feat.qualifiers === 'object' ? feat.qualifiers : {},
      qualifierRows,
      sequenceWarnings,
      nucleotideSequence,
      aminoAcidSequence,
      nucleotideFasta,
      aminoAcidFasta,
      detailRows: buildDetailRows({ defaultLabel, feat, locationText }),
      legendName: effectiveCaption,
      appliedLegendName: effectiveCaption,
      strokeColor: currentStrokeColor,
      strokeWidth: currentStrokeWidth,
      originalStrokeColor: currentStrokeColor,
      originalStrokeWidth: currentStrokeWidth,
      labelKey: '',
      labelText: '',
      labelSourceText: '',
      labelVisibility: 'default',
      featureVisibility: visibilityMode,
      // Per-feature edits name the feature by source identity (design Q4).
      identityEditable: Boolean(featureIdentityKeyOf(feat)),
      proteinId: feat.proteinId || feat.protein_id || '',
      sourceProteinId: feat.sourceProteinId || feat.source_protein_id || '',
      orthogroupId: feat.orthogroupId || '',
      orthogroupMemberCount: feat.orthogroupMemberCount || 0,
      orthogroupRecordCoverage: feat.orthogroupRecordCoverage || 0,
      orthogroupRepresentative: Boolean(feat.orthogroupRepresentative),
      orthogroupMember: feat.orthogroupMember || null,
      hasEditableLabel: false,
      labelUnavailableReason: 'No editable feature label for this feature in current diagram.'
    };
  };

  /**
   * @param {Record<string, any>} feat
   * @param {{ clientX: number, clientY: number } | null} [eventLike]
   */
  const openPreparedFeatureEditor = (feat, eventLike = null) => {
    if (!feat) return null;
    if (!svgContainer.value) return null;
    const svg = svgContainer.value.querySelector('svg');
    if (!svg) return null;

    // A catalog feature opens as the displayed Result draws it; one that the
    // Result lists but does not draw (a hidden feature, R-5) opens without an
    // element, so its popup can show it again.
    let target = feat;
    let renderedSvgId = renderedFeatureSvgId(feat);
    const catalogFeatures = resultCatalogFeatures(state);
    const key = stableFeatureOverrideKey(feat);
    if (catalogFeatures && key) {
      target = catalogFeatures.renderedByIdentity.get(key)
        || catalogFeatures.biological.find((candidate) => stableFeatureOverrideKey(candidate) === key);
      if (!target) return null;
      renderedSvgId = catalogFeatures.renderedByIdentity.has(key) ? renderedFeatureSvgId(target) : '';
    } else if (!renderedSvgId) {
      return null;
    }
    const featureElements = renderedSvgId ? getFeatureElements(svg, renderedSvgId) : [];
    if (renderedSvgId && featureElements.length === 0) return null;

    hideHoverSummary();
    if (clickedPairwiseMatch) clickedPairwiseMatch.value = null;
    const featureElement = renderedSvgId
      ? getFeatureFillElements(svg, renderedSvgId)[0] || featureElements[0] || null
      : null;
    clickedFeature.value = buildClickedFeaturePayload(target, featureElement, renderedSvgId);
    if (featurePopupSize) {
      featurePopupSize.width = 0;
      featurePopupSize.height = 0;
    }

    const popupPosition = getPopupPosition(eventLike, adv?.rich_feature_popup === false ? 440 : 720);
    clickedFeaturePos.x = popupPosition.x;
    clickedFeaturePos.y = popupPosition.y;
    if (typeof onFeaturePopupOpened === 'function') {
      onFeaturePopupOpened();
    }
    return clickedFeature.value;
  };

  // The popup states whether the feature is drawn (resolveFeatureDrawn).
  /**
   * @param {Record<string, any>} feat
   * @param {{ clientX: number, clientY: number } | null} [eventLike]
   */
  const openFeatureEditorForFeature = (feat, eventLike = null) => reportRuleRunFailure(
    state, 'feature-extraction', () => runWhenPrepared(
      state, () => [prepareDrawnFeatureMatches({ strict: true })], () => openPreparedFeatureEditor(feat, eventLike)
    )
  );

  const hoverSummaryIsAllowed = () => {
    if (clickedFeature.value) return false;
    if (clickedPairwiseMatch?.value) return false;
    if (isPreviewTransformInteractionActive()) return false;
    if (window.matchMedia && !window.matchMedia('(hover: hover) and (pointer: fine)').matches) {
      return false;
    }
    return true;
  };

  const ensureHoverSummaryElement = () => {
    if (hoverSummaryState.element?.isConnected) return hoverSummaryState.element;
    const element = createHoverSummaryElement('div', 'feature-hover-summary');
    element.hidden = true;
    element.setAttribute('role', 'tooltip');
    document.body.appendChild(element);
    hoverSummaryState.element = element;
    return element;
  };

  const positionHoverSummary = (eventLike = hoverSummaryState.lastEvent) => {
    const element = hoverSummaryState.element;
    if (!element || element.hidden || !eventLike) return;
    const margin = 12;
    const offset = 14;
    // Number.isFinite is true only for a number, so a finite clientX is a number.
    const clientX = Number.isFinite(eventLike.clientX) ? /** @type {number} */ (eventLike.clientX) : window.innerWidth / 2;
    const clientY = Number.isFinite(eventLike.clientY) ? /** @type {number} */ (eventLike.clientY) : window.innerHeight / 2;
    const rect = element.getBoundingClientRect();
    let x = clientX + offset;
    let y = clientY + offset;

    if (x + rect.width + margin > window.innerWidth) x = clientX - rect.width - offset;
    if (y + rect.height + margin > window.innerHeight) y = clientY - rect.height - offset;
    x = Math.min(Math.max(x, margin), Math.max(margin, window.innerWidth - rect.width - margin));
    y = Math.min(Math.max(y, margin), Math.max(margin, window.innerHeight - rect.height - margin));
    element.style.left = `${Math.round(x)}px`;
    element.style.top = `${Math.round(y)}px`;
  };

  const scheduleHoverSummaryPosition = (eventLike) => {
    hoverSummaryState.lastEvent = {
      clientX: Number.isFinite(eventLike?.clientX) ? eventLike.clientX : hoverSummaryState.lastEvent?.clientX,
      clientY: Number.isFinite(eventLike?.clientY) ? eventLike.clientY : hoverSummaryState.lastEvent?.clientY
    };
    if (!hoverSummaryState.visible || hoverSummaryState.frame) return;
    hoverSummaryState.frame = window.requestAnimationFrame(() => {
      hoverSummaryState.frame = null;
      positionHoverSummary();
    });
  };

  const renderHoverSummary = (feat, featureElement, eventLike) => {
    if (!feat || !hoverSummaryIsAllowed()) {
      hideHoverSummary();
      return;
    }
    const element = ensureHoverSummaryElement();
    const primaryLabel = getHoverSummaryPrimaryLabel(feat);
    const featureType = String(feat?.type || 'Feature').trim() || 'Feature';
    const titleText = primaryLabel ? `${featureType}: ${primaryLabel}` : featureType;
    const locationText = formatFeatureLocation(feat);
    const color = resolveColorToHex(
      featureElement?.getAttribute?.('fill') || getFeatureColor(feat) || '#94a3b8'
    ) || '#94a3b8';

    element.replaceChildren();
    const title = createHoverSummaryElement('div', 'feature-hover-summary-title');
    const swatch = createHoverSummaryElement('div', 'feature-hover-summary-swatch');
    swatch.style.backgroundColor = color;
    const titleTextWrap = createHoverSummaryElement('div', 'feature-hover-summary-text');
    titleTextWrap.appendChild(createHoverSummaryElement('div', 'feature-hover-summary-heading', titleText));
    titleTextWrap.appendChild(createHoverSummaryElement('div', 'feature-hover-summary-subtitle', locationText));
    title.appendChild(swatch);
    title.appendChild(titleTextWrap);
    element.appendChild(title);

    buildHoverSummaryRows(feat, primaryLabel).forEach(([label, value, clamp]) => {
      addHoverSummaryRow(element, label, value, { clamp: Boolean(clamp) });
    });

    element.hidden = false;
    hoverSummaryState.visible = true;
    hoverSummaryState.activeSvgId = renderedFeatureSvgId(feat);
    scheduleHoverSummaryPosition(eventLike);
    positionHoverSummary(eventLike);
  };

  const scheduleHoverSummary = (feat, featureElement, eventLike) => {
    if (hoverSummaryState.timer) {
      window.clearTimeout(hoverSummaryState.timer);
      hoverSummaryState.timer = null;
    }
    if (!feat || !hoverSummaryIsAllowed()) {
      hideHoverSummary();
      return;
    }
    scheduleHoverSummaryPosition(eventLike);
    const show = () => {
      hoverSummaryState.timer = null;
      renderHoverSummary(feat, featureElement, eventLike);
    };
    if (hoverSummaryState.visible) {
      show();
      return;
    }
    hoverSummaryState.timer = window.setTimeout(show, 180);
  };

  const renderMatchHoverSummary = (summary, eventLike) => {
    if (!summary || !hoverSummaryIsAllowed()) {
      hideHoverSummary();
      return;
    }
    const element = ensureHoverSummaryElement();
    const color = resolveColorToHex(summary.fill || '#94a3b8') || '#94a3b8';

    element.replaceChildren();
    const title = createHoverSummaryElement('div', 'feature-hover-summary-title');
    const swatch = createHoverSummaryElement('div', 'feature-hover-summary-swatch');
    swatch.style.backgroundColor = color;
    const titleTextWrap = createHoverSummaryElement('div', 'feature-hover-summary-text');
    titleTextWrap.appendChild(createHoverSummaryElement('div', 'feature-hover-summary-heading', summary.title));
    titleTextWrap.appendChild(createHoverSummaryElement('div', 'feature-hover-summary-subtitle', summary.subtitle || summary.id || ''));
    title.appendChild(swatch);
    title.appendChild(titleTextWrap);
    element.appendChild(title);

    summary.rows.forEach((row) => {
      addHoverSummaryRow(element, row.label, row.value, { clamp: row.label === 'Query' || row.label === 'Subject' });
    });

    element.hidden = false;
    hoverSummaryState.visible = true;
    hoverSummaryState.activeSvgId = String(summary.id || '').trim();
    scheduleHoverSummaryPosition(eventLike);
    positionHoverSummary(eventLike);
  };

  const scheduleMatchHoverSummary = (summary, eventLike) => {
    if (hoverSummaryState.timer) {
      window.clearTimeout(hoverSummaryState.timer);
      hoverSummaryState.timer = null;
    }
    if (!summary || !hoverSummaryIsAllowed()) {
      hideHoverSummary();
      return;
    }
    scheduleHoverSummaryPosition(eventLike);
    const show = () => {
      hoverSummaryState.timer = null;
      renderMatchHoverSummary(summary, eventLike);
    };
    if (hoverSummaryState.visible) {
      show();
      return;
    }
    hoverSummaryState.timer = window.setTimeout(show, 180);
  };

  function hideHoverSummary() {
    if (hoverSummaryState.timer) {
      window.clearTimeout(hoverSummaryState.timer);
      hoverSummaryState.timer = null;
    }
    if (hoverSummaryState.frame) {
      window.cancelAnimationFrame(hoverSummaryState.frame);
      hoverSummaryState.frame = null;
    }
    if (hoverSummaryState.element) {
      hoverSummaryState.element.hidden = true;
    }
    hoverSummaryState.visible = false;
    hoverSummaryState.activeSvgId = '';
    hoverSummaryState.lastEvent = null;
  }

  const groupsForMatch = (matchElement) => (
    matchElement.getAttribute('data-match-kind') === 'collinear'
      ? collinearGroups?.value || []
      : orthogroups?.value || []
  );

  const buildMatchPayload = (matchElement, featureLookup) => buildMatchPopupPayload(matchElement, {
    featureLookup,
    sourceFeatures: Array.isArray(biologicalFeatures?.value) && biologicalFeatures.value.length > 0
      ? biologicalFeatures.value
      : (Array.isArray(extractedFeatures.value) ? extractedFeatures.value : []),
    orthogroups: groupsForMatch(matchElement),
    orthogroupNameOverrides,
    orthogroupDescriptionOverrides,
    resolveSequenceSource: matchSequenceRegistry?.resolve
  });

  // The summary builder reads a function as well as a list.
  const buildMatchHoverSummary = (matchElement) => buildPairwiseMatchHoverSummary(matchElement, /** @type {any} */ ({
    orthogroups: () => groupsForMatch(matchElement),
    orthogroupNameOverrides
  }));

  const openPairwiseMatchPopup = (matchElement, eventLike, featureLookup) => {
    if (!matchElement || !clickedPairwiseMatch || !clickedPairwiseMatchPos) return null;
    const payload = buildMatchPayload(matchElement, featureLookup);
    if (!payload) return null;
    hideHoverSummary();
    clickedFeature.value = null;
    clickedPairwiseMatch.value = payload;
    const popupPosition = getPopupPosition(eventLike, 460, 520);
    clickedPairwiseMatchPos.x = popupPosition.x;
    clickedPairwiseMatchPos.y = popupPosition.y;
    return payload;
  };

  const applyVisibilityPreviewChanges = (changes, { reason = 'feature-visibility' } = {}) => {
    const normalizedChanges = (Array.isArray(changes) ? changes : [])
      .map((change) => ({
        featureId: String(change?.featureId || change?.svgId || change?.id || '').trim(),
        mode: normalizeVisibilityMode(change?.mode)
      }))
      .filter((change) => change.featureId);
    if (normalizedChanges.length === 0) return false;
    return applyFeatureVisibilityChanges?.(normalizedChanges, { reason }) === true;
  };

  const applyVisibilityPreviewBySvgId = (svgId, modeRaw) => (
    applyVisibilityPreviewChanges([{ featureId: svgId, mode: modeRaw }])
  );

  /** @param {{ root?: Element | null, phase?: string, rootGeneration?: number }} [options] */
  const attachSvgFeatureHandlers = ({
    root = null,
    phase = 'preview-bind',
    rootGeneration = 0
  } = {}) => {
    const svg = root || svgContainer.value?.querySelector?.('svg') || null;
    if (!svg) return false;

    if (delegatedFeatureHandlers && delegatedFeatureHandlers.svg !== svg) {
      cleanupDelegatedFeatureHandlers();
    }

    const buildFeatureIdsByOrthogroupId = (featureLookup) => {
      const featureIdsByOrthogroupId = new Map();
      featureLookup.forEach((feat, svgId) => {
        if (!svgId) return;
        getOrthogroupIds(feat?.orthogroupId).forEach((orthogroupId) => {
          if (!featureIdsByOrthogroupId.has(orthogroupId)) {
            featureIdsByOrthogroupId.set(orthogroupId, new Set());
          }
          featureIdsByOrthogroupId.get(orthogroupId).add(svgId);
        });
      });
      return featureIdsByOrthogroupId;
    };

    if (delegatedFeatureHandlers?.svg === svg) {
      return false;
    }

    /** @type {DelegatedFeatureHandlers} */
    const handlerState = {
      svg,
      pathsByIdMap: null,
      featureLookup: null,
      featureIdsByOrthogroupId: null,
      comparisonElementsByOrthogroupId: null,
      comparisonElementsByCollinearityBlockId: null,
      comparisonElementsByMatchId: null,
      pairwiseAffordancesPrepared: false,
      activeHoverSvgId: null,
      activeHoverKey: '',
      activeMatchHoverElement: null,
      activeMatchHoverKey: '',
      pendingMatchElement: null,
      activeHoverElements: new Set(),
      transformPointer: null,
      reconcileHoverAfterTransform: false,
      hoverReconcileFrame: null,
      alignmentCandidateSvgId: '',
      alignmentCandidateRestore: null,
      alignmentOverlay: null,
      alignmentCandidatesBySvgId: new Map(),
      showAlignmentOverlay: null,
      clearAlignmentOverlay: null,
      previewAlignmentCandidate: null,
      clearAlignmentCandidatePreview: null,
      beginPreviewTransformInteraction: null,
      endPreviewTransformInteraction: null,
      cleanup: null
    };

    const ensureFeatureLookup = () => {
      if (!handlerState.featureLookup) handlerState.featureLookup = buildFeatureLookup();
      return handlerState.featureLookup;
    };

    const ensureFeaturePaths = () => {
      if (!handlerState.pathsByIdMap) {
        recordStructuralMetric('featureDomFullScanCount', 1, { phase: 'interaction' });
        // getFeatureElementIndex always returns a Map (the cached index or a new one).
        handlerState.pathsByIdMap = /** @type {Map<string, Element[]>} */ (getFeatureElementIndex(svg));
      }
      return handlerState.pathsByIdMap;
    };

    const ensureFeatureOrthogroupIndex = () => {
      if (!handlerState.featureIdsByOrthogroupId) {
        handlerState.featureIdsByOrthogroupId = buildFeatureIdsByOrthogroupId(
          ensureFeatureLookup()
        );
      }
      return handlerState.featureIdsByOrthogroupId;
    };

    const ensureComparisonIndexes = () => {
      if (
        handlerState.comparisonElementsByOrthogroupId
        && handlerState.comparisonElementsByCollinearityBlockId
      ) {
        return;
      }
      const byOrthogroup = new Map();
      const byBlock = new Map();
      const byMatch = new Map();
      svg.querySelectorAll('[data-orthogroup-id]').forEach((element) => {
        if (element.matches?.(FEATURE_SELECTOR)) return;
        getOrthogroupIds(element.getAttribute('data-orthogroup-id')).forEach((orthogroupId) => {
          if (!byOrthogroup.has(orthogroupId)) byOrthogroup.set(orthogroupId, []);
          byOrthogroup.get(orthogroupId).push(element);
        });
      });
      svg.querySelectorAll(PAIRWISE_MATCH_SELECTOR).forEach((element) => {
        const matchId = String(element.getAttribute('data-gbdraw-match-id') || element.getAttribute('data-gbdraw-pairwise-match-id') || '').trim();
        if (matchId) {
          if (!byMatch.has(matchId)) byMatch.set(matchId, []);
          byMatch.get(matchId).push(element);
        }
        const blockId = String(element.getAttribute('data-collinearity-block-id') || '').trim();
        if (!blockId) return;
        if (!byBlock.has(blockId)) byBlock.set(blockId, []);
        byBlock.get(blockId).push(element);
      });
      handlerState.comparisonElementsByOrthogroupId = byOrthogroup;
      handlerState.comparisonElementsByCollinearityBlockId = byBlock;
      handlerState.comparisonElementsByMatchId = byMatch;
    };

      const setHoverStyle = (element, highlight) => {
        if (!element?.style) return;
        if (highlight) {
          if (!element.hasAttribute('data-gbdraw-hover-opacity')) {
            element.setAttribute('data-gbdraw-hover-opacity', element.style.opacity || '');
            element.setAttribute('data-gbdraw-hover-filter', element.style.filter || '');
          }
          element.style.opacity = '0.7';
          element.style.filter = 'brightness(1.2)';
          handlerState.activeHoverElements.add(element);
          return;
        }
        if (element.hasAttribute('data-gbdraw-hover-opacity')) {
          element.style.opacity = element.getAttribute('data-gbdraw-hover-opacity') || '';
          element.style.filter = element.getAttribute('data-gbdraw-hover-filter') || '';
          element.removeAttribute('data-gbdraw-hover-opacity');
          element.removeAttribute('data-gbdraw-hover-filter');
        }
        handlerState.activeHoverElements.delete(element);
      };

      const setFeatureHover = (svgId) => {
        (ensureFeaturePaths().get(svgId) || []).forEach((element) => {
          setHoverStyle(element, true);
        });
      };

      const getFeatureHoverKey = (svgId) => {
        const feat = ensureFeatureLookup().get(svgId);
        const orthogroupId = String(feat?.orthogroupId || '').trim();
        return orthogroupId ? `orthogroup:${orthogroupId}` : `feature:${svgId}`;
      };

      const setOrthogroupHover = (orthogroupId) => {
        const id = String(orthogroupId || '').trim();
        if (!id) return;
        ensureComparisonIndexes();
        (ensureFeatureOrthogroupIndex().get(id) || new Set()).forEach((featureId) => {
          setFeatureHover(featureId);
        });
        // ensureComparisonIndexes above has set all three comparison maps.
        (/** @type {Map<string, Element[]>} */ (handlerState.comparisonElementsByOrthogroupId).get(id) || []).forEach((element) => {
          setHoverStyle(element, true);
        });
      };

      const setCollinearityBlockHover = (blockId) => {
        const id = String(blockId || '').trim();
        if (!id) return;
        ensureComparisonIndexes();
        // ensureComparisonIndexes above has set all three comparison maps.
        (/** @type {Map<string, Element[]>} */ (handlerState.comparisonElementsByCollinearityBlockId).get(id) || []).forEach((element) => {
          setHoverStyle(element, true);
        });
      };

      const setHoverHighlight = (svgId) => {
        const feat = ensureFeatureLookup().get(svgId);
        const orthogroupId = String(feat?.orthogroupId || '').trim();
        if (orthogroupId) {
          setOrthogroupHover(orthogroupId);
          return;
        }
        setFeatureHover(svgId);
      };

      const matchAttr = (element, name) => String(element?.getAttribute?.(name) || '').trim();

      const getMatchHoverKey = (matchElement) => {
        const blockId = matchAttr(matchElement, 'data-collinearity-block-id');
        if (blockId) return `collinearity:${blockId}`;
        const orthogroupId = matchAttr(matchElement, 'data-orthogroup-id');
        if (orthogroupId) return `orthogroup:${orthogroupId}`;
        return `match:${matchAttr(matchElement, 'data-gbdraw-match-id') || matchAttr(matchElement, 'data-gbdraw-pairwise-match-id') || matchAttr(matchElement, 'd')}`;
      };

      const matchFragments = (element) => {
        ensureComparisonIndexes();
        const id = matchAttr(element, 'data-gbdraw-match-id') || matchAttr(element, 'data-gbdraw-pairwise-match-id');
        // ensureComparisonIndexes above has set all three comparison maps.
        return /** @type {Map<string, Element[]>} */ (handlerState.comparisonElementsByMatchId).get(id) || [element];
      };

      const setMatchHover = (matchElement) => {
        if (!matchElement) return;
        const blockId = matchAttr(matchElement, 'data-collinearity-block-id');
        const orthogroupId = matchAttr(matchElement, 'data-orthogroup-id');
        matchFragments(matchElement).forEach((element) => setHoverStyle(element, true));
        if (blockId) {
          setCollinearityBlockHover(blockId);
          return;
        }
        if (orthogroupId) {
          setOrthogroupHover(orthogroupId);
        }
      };

      const clearTrackedHoverStyles = () => {
        [...handlerState.activeHoverElements].forEach((element) => {
          setHoverStyle(element, false);
        });
      };

      const clearActiveFeatureHover = () => {
        if (!handlerState.activeHoverSvgId) return;
        clearTrackedHoverStyles();
        handlerState.activeHoverSvgId = null;
        handlerState.activeHoverKey = '';
      };

      const clearActiveMatchHover = () => {
        if (!handlerState.activeMatchHoverElement) return;
        clearTrackedHoverStyles();
        handlerState.activeMatchHoverElement = null;
        handlerState.activeMatchHoverKey = '';
      };

      const clearAlignmentCandidatePreview = ({ restore = true } = {}) => {
        if (!handlerState.alignmentCandidateSvgId && !handlerState.alignmentCandidateRestore) return;
        clearTrackedHoverStyles();
        handlerState.alignmentCandidateSvgId = '';
        const previous = handlerState.alignmentCandidateRestore;
        handlerState.alignmentCandidateRestore = null;
        if (!restore || !previous) return;
        if (previous.featureSvgId && ensureFeatureLookup().has(previous.featureSvgId)) {
          setHoverHighlight(previous.featureSvgId);
          handlerState.activeHoverSvgId = previous.featureSvgId;
          handlerState.activeHoverKey = getFeatureHoverKey(previous.featureSvgId);
          return;
        }
        if (previous.matchElement?.isConnected) {
          setMatchHover(previous.matchElement);
          handlerState.activeMatchHoverElement = previous.matchElement;
          handlerState.activeMatchHoverKey = getMatchHoverKey(previous.matchElement);
        }
      };

      const previewAlignmentCandidate = (anchor) => {
        const requested = featureIdentity(anchor);
        if (!requested.usable) return false;
        const matches = Array.from(ensureFeatureLookup().entries())
          .filter(([, feature]) => {
            const candidate = renderedFeatureIdentity(feature);
            return candidate.usable
              && identityMatches(requested, candidate);
          });
        if (matches.length !== 1) return false;
        const restore = handlerState.alignmentCandidateRestore || {
          featureSvgId: handlerState.activeHoverSvgId,
          matchElement: handlerState.activeMatchHoverElement
        };
        if (handlerState.alignmentCandidateSvgId) {
          clearTrackedHoverStyles();
          handlerState.alignmentCandidateSvgId = '';
        }
        clearActiveFeatureHover();
        clearActiveMatchHover();
        hideHoverSummary();
        handlerState.alignmentCandidateRestore = restore;
        const [svgId] = matches[0];
        setFeatureHover(svgId);
        handlerState.alignmentCandidateSvgId = svgId;
        return true;
      };

      handlerState.previewAlignmentCandidate = previewAlignmentCandidate;
      handlerState.clearAlignmentCandidatePreview = clearAlignmentCandidatePreview;

      const uniqueRenderedFeature = (anchor) => {
        const requested = featureIdentity(anchor);
        if (!requested.usable) return null;
        const matches = Array.from(ensureFeatureLookup().entries()).filter(([, feature]) => (
          identityMatches(requested, renderedFeatureIdentity(feature))
        ));
        if (matches.length !== 1) return null;
        const [svgId] = matches[0];
        const blocks = getFeatureFillElements(svg, svgId, ensureFeaturePaths());
        return blocks.length === 1 ? { svgId, element: blocks[0] } : null;
      };

      const clearAlignmentOverlay = () => {
        const overlay = handlerState.alignmentOverlay;
        if (!overlay) return;
        if (overlay.frame) window.cancelAnimationFrame(overlay.frame);
        overlay.scrollRoot.removeEventListener('scroll', overlay.schedule);
        window.removeEventListener('resize', overlay.schedule);
        overlay.observer?.disconnect();
        overlay.element.remove();
        handlerState.alignmentOverlay = null;
        handlerState.alignmentCandidatesBySvgId.clear();
        clearAlignmentCandidatePreview();
      };

      /** @param {AlignmentOverlayRequest} request */
      const showAlignmentOverlay = ({ reference, ambiguities, onSelect, onHover }) => {
        clearAlignmentOverlay();
        const scrollRoot = svgContainer.value?.parentElement;
        const viewport = scrollRoot?.parentElement;
        if (!scrollRoot || !viewport || !svg.isConnected) return;
        const referenceFeature = uniqueRenderedFeature(reference);
        const candidates = [];
        (ambiguities || []).forEach((record) => record.candidates.forEach((candidate, index) => {
          const rendered = uniqueRenderedFeature(candidate.anchor);
          if (rendered) candidates.push({
            ...rendered, recordKey: record.recordKey, anchor: candidate.anchor,
            key: candidate.key, number: index + 1,
            accessibleName: 'Select ' + record.recordLabel + ' candidate ' + (index + 1)
              + ': ' + candidate.label + ', ' + candidate.coordinates
          });
        }));
        const counts = new Map();
        candidates.forEach(({ svgId }) => counts.set(svgId, (counts.get(svgId) || 0) + 1));
        const visibleCandidates = candidates.filter(({ svgId }) => counts.get(svgId) === 1);
        visibleCandidates.forEach((candidate) => {
          handlerState.alignmentCandidatesBySvgId.set(candidate.svgId, candidate);
        });

        const element = document.createElement('div');
        element.className = 'gbdraw-alignment-overlay';
        element.setAttribute('data-similarity-alignment-canvas', '');
        const guide = document.createElement('div');
        guide.className = 'gbdraw-alignment-guide';
        guide.setAttribute('aria-hidden', 'true');
        element.appendChild(guide);
        const badges = visibleCandidates.map((candidate) => {
          const button = document.createElement('button');
          button.type = 'button';
          button.className = 'gbdraw-alignment-badge';
          button.textContent = String(candidate.number);
          button.setAttribute('data-alignment-record-key', candidate.recordKey);
          button.setAttribute('data-alignment-candidate-key', candidate.key);
          button.setAttribute('aria-label', candidate.accessibleName);
          button.addEventListener('click', (event) => {
            event.preventDefault();
            event.stopPropagation();
            onSelect(candidate.recordKey, candidate.anchor);
          });
          button.addEventListener('mouseenter', () => onHover(candidate.recordKey, candidate.key));
          button.addEventListener('mouseleave', () => onHover('', ''));
          button.addEventListener('focus', () => onHover(candidate.recordKey, candidate.key));
          button.addEventListener('blur', () => onHover('', ''));
          element.appendChild(button);
          return { candidate, button };
        });
        viewport.appendChild(element);

        const geometry = (target, viewportRect) => {
          if (!target?.isConnected || target.getClientRects().length !== 1
            || window.getComputedStyle(target).visibility === 'hidden') return null;
          const rect = target.getBoundingClientRect();
          if (!rect.width || !rect.height || rect.right <= viewportRect.left
            || rect.left >= viewportRect.right || rect.bottom <= viewportRect.top) return null;
          return {
            x: (rect.left + rect.right) / 2 - viewportRect.left,
            y: (rect.top + rect.bottom) / 2 - viewportRect.top
          };
        };
        /** @type {AlignmentOverlay} */
        const overlay = {
          element,
          scrollRoot,
          frame: 0,
          observer: null,
          schedule: () => {
            if (!overlay.frame) overlay.frame = window.requestAnimationFrame(update);
          },
          onSelect,
          onHover
        };
        const update = () => {
          overlay.frame = 0;
          if (!element.isConnected || delegatedFeatureHandlers !== handlerState) return;
          const rect = viewport.getBoundingClientRect();
          const referencePoint = geometry(referenceFeature?.element, rect);
          guide.hidden = !referencePoint || referencePoint.x < 0 || referencePoint.x > rect.width;
          if (referencePoint && !guide.hidden) guide.style.left = referencePoint.x + 'px';
          badges.forEach(({ candidate, button }) => {
            const point = geometry(candidate.element, rect);
            button.hidden = !point || point.x < 0 || point.x > rect.width
              || point.y < 0 || point.y > rect.height;
            if (point && !button.hidden) {
              button.style.left = point.x + 'px';
              button.style.top = point.y + 'px';
            }
          });
        };
        scrollRoot.addEventListener('scroll', overlay.schedule, { passive: true });
        window.addEventListener('resize', overlay.schedule);
        if (window.ResizeObserver) {
          overlay.observer = new ResizeObserver(overlay.schedule);
          overlay.observer.observe(viewport);
          overlay.observer.observe(svgContainer.value);
        }
        handlerState.alignmentOverlay = overlay;
        overlay.schedule();
      };

      handlerState.showAlignmentOverlay = showAlignmentOverlay;
      handlerState.clearAlignmentOverlay = clearAlignmentOverlay;

      const clearPendingMatch = () => {
        if (!handlerState.pendingMatchElement) return;
        matchFragments(handlerState.pendingMatchElement).forEach((element) => element.classList.remove('gbdraw-match-pending'));
        handlerState.pendingMatchElement = null;
      };

      const setPendingMatch = (matchElement) => {
        if (handlerState.pendingMatchElement === matchElement) return;
        clearPendingMatch();
        handlerState.pendingMatchElement = matchElement;
        matchFragments(matchElement).forEach((element) => element.classList.add('gbdraw-match-pending'));
      };

      const rememberTransformPointer = (eventLike) => {
        if (!Number.isFinite(eventLike?.clientX) || !Number.isFinite(eventLike?.clientY)) return;
        handlerState.transformPointer = {
          clientX: eventLike.clientX,
          clientY: eventLike.clientY
        };
        handlerState.reconcileHoverAfterTransform = true;
      };

      const activateFeatureHover = (featureEl, eventLike) => {
        clearActiveMatchHover();
        const svgId = getFeatureIdentity(featureEl);
        if (!svgId) return false;
        const hoverKey = getFeatureHoverKey(svgId);
        if (handlerState.activeHoverKey !== hoverKey && handlerState.activeHoverSvgId) {
          clearActiveFeatureHover();
        }
        if (handlerState.activeHoverKey !== hoverKey) {
          setHoverHighlight(svgId);
        }
        handlerState.activeHoverSvgId = svgId;
        handlerState.activeHoverKey = hoverKey;
        const candidate = handlerState.alignmentCandidatesBySvgId.get(svgId);
        handlerState.alignmentOverlay?.onHover(candidate?.recordKey || '', candidate?.key || '');
        scheduleHoverSummary(ensureFeatureLookup().get(svgId), featureEl, eventLike);
        return true;
      };

      const activateMatchHover = (matchEl, eventLike) => {
        handlerState.alignmentOverlay?.onHover('', '');
        clearActiveFeatureHover();
        const matchKey = getMatchHoverKey(matchEl);
        if (handlerState.activeMatchHoverKey !== matchKey && handlerState.activeMatchHoverElement) {
          clearActiveMatchHover();
        }
        if (handlerState.activeMatchHoverKey !== matchKey) {
          setMatchHover(matchEl);
        }
        handlerState.activeMatchHoverElement = matchEl;
        handlerState.activeMatchHoverKey = matchKey;
        scheduleMatchHoverSummary(buildMatchHoverSummary(matchEl), eventLike);
        return true;
      };

      handlerState.beginPreviewTransformInteraction = (eventLike, kind = '') => {
        if (handlerState.alignmentOverlay) handlerState.alignmentOverlay.element.hidden = true;
        handlerState.alignmentOverlay?.onHover('', '');
        if (handlerState.hoverReconcileFrame !== null) {
          window.cancelAnimationFrame(handlerState.hoverReconcileFrame);
          handlerState.hoverReconcileFrame = null;
        }
        handlerState.reconcileHoverAfterTransform = Boolean(
          handlerState.activeHoverSvgId || handlerState.activeMatchHoverElement
        );
        rememberTransformPointer(eventLike);
        clearActiveFeatureHover();
        clearActiveMatchHover();
        hideHoverSummary();
        recordStructuralMetric('previewTransformHoverCleanupCount', 1, { kind });
      };

      handlerState.endPreviewTransformInteraction = ({ reconcile = true } = {}) => {
        if (handlerState.alignmentOverlay) {
          handlerState.alignmentOverlay.element.hidden = false;
          handlerState.alignmentOverlay.schedule();
        }
        if (!reconcile || !handlerState.reconcileHoverAfterTransform || !handlerState.transformPointer) {
          handlerState.reconcileHoverAfterTransform = false;
          handlerState.transformPointer = null;
          return;
        }
        if (handlerState.hoverReconcileFrame !== null) {
          window.cancelAnimationFrame(handlerState.hoverReconcileFrame);
        }
        handlerState.hoverReconcileFrame = window.requestAnimationFrame(() => {
          handlerState.hoverReconcileFrame = null;
          if (isPreviewTransformInteractionActive() || delegatedFeatureHandlers !== handlerState) return;
          const pointer = handlerState.transformPointer;
          handlerState.transformPointer = null;
          handlerState.reconcileHoverAfterTransform = false;
          recordStructuralMetric('previewTransformHoverReconcileCount', 1);
          const target = getTopmostSvgTarget(
            pointer,
            svg,
            `${PAIRWISE_MATCH_SELECTOR}, ${FEATURE_SELECTOR}`
          );
          const matchEl = getPairwiseMatchTarget(target, svg);
          if (matchEl) {
            activateMatchHover(matchEl, pointer);
            return;
          }
          const featureEl = getFeatureTarget(target, svg);
          if (featureEl) activateFeatureHover(featureEl, pointer);
        });
      };

      const handleMouseOver = (e) => {
        if (isPreviewTransformInteractionActive()) {
          rememberTransformPointer(e);
          return;
        }
        const featureEl = getFeatureTarget(e.target, svg);
        if (featureEl) {
          activateFeatureHover(featureEl, e);
          return;
        }
        const matchEl = getPairwiseMatchTarget(e.target, svg);
        if (!matchEl) return;
        activateMatchHover(matchEl, e);
      };

      const handleMouseMove = (e) => {
        if (isPreviewTransformInteractionActive()) {
          rememberTransformPointer(e);
          return;
        }
        if (
          clickedFeature.value ||
          clickedPairwiseMatch?.value ||
          featureSelectionDrag?.active
        ) {
          hideHoverSummary();
          return;
        }
        if (hoverSummaryState.visible || hoverSummaryState.timer) {
          scheduleHoverSummaryPosition(e);
        }
      };

      const handleMouseOut = (e) => {
        if (isPreviewTransformInteractionActive()) {
          rememberTransformPointer(e);
          return;
        }
        const featureEl = getFeatureTarget(e.target, svg);
        if (featureEl) {
          const svgId = getFeatureIdentity(featureEl);
          if (!svgId || handlerState.activeHoverSvgId !== svgId) return;
          const relatedFeature = getFeatureTarget(e.relatedTarget, svg);
          if (relatedFeature && getFeatureHoverKey(getFeatureIdentity(relatedFeature)) === handlerState.activeHoverKey) return;
          clearActiveFeatureHover();
          handlerState.alignmentOverlay?.onHover('', '');
          hideHoverSummary();
          return;
        }
        const matchEl = getPairwiseMatchTarget(e.target, svg);
        if (!matchEl || handlerState.activeMatchHoverElement !== matchEl) return;
        const relatedMatch = getPairwiseMatchTarget(e.relatedTarget, svg);
        if (relatedMatch && getMatchHoverKey(relatedMatch) === handlerState.activeMatchHoverKey) return;
        clearActiveMatchHover();
        hideHoverSummary();
      };

      const handleClick = (e) => {
        clearPendingMatch();
        if (featureSelection?.consumeSuppressNextClick?.()) {
          e.preventDefault();
          e.stopPropagation();
          return;
        }
        const selectableFeatureEl = featureSelection?.getSelectableFeatureTarget?.(e, svg) || null;
        const modifierSelection = Boolean(e.ctrlKey || e.metaKey);
        if (modifierSelection || e.shiftKey) {
          if (selectableFeatureEl) {
            const selectionId = getFeatureIdentity(selectableFeatureEl);
            e.preventDefault();
            e.stopPropagation();
            hideHoverSummary();
            if (modifierSelection && !e.shiftKey) {
              featureSelection?.toggleFeatureSelection?.(selectionId);
            } else {
              featureSelection?.selectFeatureRange?.(selectionId, { additive: modifierSelection });
            }
          } else if (isBackgroundPreviewClick(e.target, svg)) {
            featureSelection?.clearFeatureSelection?.({ clearStatus: true });
          }
          return;
        }
        const annotationEl = e.target?.closest?.('[data-gbdraw-annotation-id]');
        if (annotationEl && svg.contains(annotationEl)) {
          e.preventDefault();
          e.stopPropagation();
          if (selectedAnnotation) {
            selectedAnnotation.value = {
              id: annotationEl.getAttribute('data-gbdraw-annotation-id') || '',
              setId: annotationEl.getAttribute('data-gbdraw-annotation-set-id') || '',
              trackId: annotationEl.getAttribute('data-gbdraw-annotation-track-id') || ''
            };
          }
          hideHoverSummary();
          return;
        }

        const featureEl = getFeatureClickTarget(e, svg);
        if (featureEl) {
          const svgId = getFeatureIdentity(featureEl);
          if (!svgId) return;
          const candidate = handlerState.alignmentCandidatesBySvgId.get(svgId);
          if (candidate) {
            e.preventDefault();
            e.stopPropagation();
            hideHoverSummary();
            handlerState.alignmentOverlay?.onSelect(candidate.recordKey, candidate.anchor);
            return;
          }
          e.stopPropagation();
          hideHoverSummary();
          featureSelection?.markPlainFeatureClick?.(svgId);
          const feat = ensureFeatureLookup().get(svgId);
          if (feat) {
            openFeatureEditorForFeature(feat, e);
          } else {
            console.log(`No feature found for svg_id: ${svgId}`);
          }
          return;
        }
        const matchEl = getPairwiseMatchClickTarget(e, svg);
        if (matchEl) {
          e.stopPropagation();
          e.preventDefault();
          openPairwiseMatchPopup(matchEl, e, ensureFeatureLookup());
          return;
        }
        if (isBackgroundPreviewClick(e.target, svg)) {
          featureSelection?.clearFeatureSelection?.({ clearStatus: true });
          if (clickedFeature.value) clickedFeature.value = null;
          if (clickedPairwiseMatch?.value) clickedPairwiseMatch.value = null;
        }
      };

      const handleKeyDown = (e) => {
        if (e.key !== 'Enter' && e.key !== ' ') return;
        const matchEl = getPairwiseMatchTarget(e.target, svg);
        if (!matchEl) return;
        e.stopPropagation();
        e.preventDefault();
        openPairwiseMatchPopup(matchEl, e, ensureFeatureLookup());
      };

      const handlePointerDown = (e) => {
        if (e.button === 0 && e.isPrimary !== false) {
          const matchEl = getPairwiseMatchClickTarget(e, svg);
          if (matchEl) {
            setPendingMatch(matchEl);
            return;
          }
        }
        clearPendingMatch();
        if (!featureSelection?.startMarqueePointer?.(e, svg)) return;
        e.stopPropagation();
      };

      const handlePointerOut = (e) => {
        if (!handlerState.pendingMatchElement) return;
        const matchEl = getPairwiseMatchTarget(e.target, svg);
        if (matchEl !== handlerState.pendingMatchElement) return;
        if (getPairwiseMatchTarget(e.relatedTarget, svg) === matchEl) return;
        clearPendingMatch();
      };

      const handlePointerMove = (e) => {
        if (!featureSelection?.moveMarqueePointer?.(e)) return;
        e.preventDefault();
        e.stopPropagation();
        hideHoverSummary();
      };

      const handlePointerUp = (e) => {
        clearPendingMatch();
        if (!featureSelection?.commitMarqueePointer?.(e)) return;
        e.preventDefault();
        e.stopPropagation();
        hideHoverSummary();
      };

      const handlePointerCancel = () => {
        clearPendingMatch();
        featureSelection?.cancelMarqueePointer?.();
      };

      svg.addEventListener('mouseover', handleMouseOver);
      svg.addEventListener('mousemove', handleMouseMove);
      svg.addEventListener('mouseout', handleMouseOut);
      svg.addEventListener('click', handleClick);
      svg.addEventListener('keydown', handleKeyDown);
      svg.addEventListener('pointerdown', handlePointerDown, true);
      svg.addEventListener('pointerout', handlePointerOut, true);
      svg.addEventListener('pointermove', handlePointerMove, true);
      svg.addEventListener('pointerup', handlePointerUp, true);
      svg.addEventListener('pointercancel', handlePointerCancel, true);
      window.addEventListener?.('pointerup', clearPendingMatch, true);
      window.addEventListener?.('pointercancel', clearPendingMatch, true);
      window.addEventListener?.('blur', clearPendingMatch);
      handlerState.cleanup = () => {
        svg.removeEventListener('mouseover', handleMouseOver);
        svg.removeEventListener('mousemove', handleMouseMove);
        svg.removeEventListener('mouseout', handleMouseOut);
        svg.removeEventListener('click', handleClick);
        svg.removeEventListener('keydown', handleKeyDown);
        svg.removeEventListener('pointerdown', handlePointerDown, true);
        svg.removeEventListener('pointerout', handlePointerOut, true);
        svg.removeEventListener('pointermove', handlePointerMove, true);
        svg.removeEventListener('pointerup', handlePointerUp, true);
        svg.removeEventListener('pointercancel', handlePointerCancel, true);
        window.removeEventListener?.('pointerup', clearPendingMatch, true);
        window.removeEventListener?.('pointercancel', clearPendingMatch, true);
        window.removeEventListener?.('blur', clearPendingMatch);
        if (handlerState.hoverReconcileFrame !== null) {
          window.cancelAnimationFrame(handlerState.hoverReconcileFrame);
          handlerState.hoverReconcileFrame = null;
        }
        handlerState.reconcileHoverAfterTransform = false;
        handlerState.transformPointer = null;
        clearPendingMatch();
        clearAlignmentOverlay();
        clearAlignmentCandidatePreview({ restore: false });
        clearActiveFeatureHover();
        clearActiveMatchHover();
        hideHoverSummary();
      };
      delegatedFeatureHandlers = handlerState;
    recordStructuralMetric('perFeatureListenerRegistrationCount', 0, {
      phase,
      rootGeneration
    });
    return true;
  };

  const unsubscribePreviewTransformInteraction = previewTransformInteraction?.subscribe?.(({
    active,
    kind,
    event,
    reconcile
  }) => {
    if (active) {
      delegatedFeatureHandlers?.beginPreviewTransformInteraction?.(event, kind);
      return;
    }
    delegatedFeatureHandlers?.endPreviewTransformInteraction?.({ reconcile });
  }) || null;

  const dispose = () => {
    if (typeof unsubscribePreviewTransformInteraction === 'function') {
      unsubscribePreviewTransformInteraction();
    }
    cleanupDelegatedFeatureHandlers();
    if (hoverSummaryState.element?.isConnected) hoverSummaryState.element.remove();
    hoverSummaryState.element = null;
  };

  const previewAlignmentCandidate = (anchor) => (
    delegatedFeatureHandlers?.previewAlignmentCandidate?.(anchor) || false
  );

  const clearAlignmentCandidatePreview = () => {
    delegatedFeatureHandlers?.clearAlignmentCandidatePreview?.();
  };

  const showAlignmentOverlay = (request) => {
    delegatedFeatureHandlers?.showAlignmentOverlay?.(request);
  };

  const clearAlignmentOverlay = () => {
    delegatedFeatureHandlers?.clearAlignmentOverlay?.();
  };

  /** @param {{ root?: Element | null, phase?: string, rootGeneration?: number }} [options] */
  const preparePairwiseInteractionAffordances = ({
    root = null,
    phase = 'preview-bind',
    rootGeneration = 0
  } = {}) => {
    const svg = root || delegatedFeatureHandlers?.svg || null;
    if (!svg || delegatedFeatureHandlers?.svg !== svg) return false;
    if (delegatedFeatureHandlers.pairwiseAffordancesPrepared) return false;
    recordStructuralMetric('comparisonDomFullScanCount', 1, { phase, rootGeneration });
    // The pairwise match selector matches SVG path elements, which have a style.
    Array.from(/** @type {NodeListOf<SVGElement>} */ (svg.querySelectorAll(PAIRWISE_MATCH_SELECTOR))).forEach((element, index) => {
      if (element?.style) element.style.cursor = 'pointer';
      element.setAttribute('role', 'button');
      element.setAttribute('tabindex', '0');
      element.setAttribute('aria-label', `Pairwise match ${index + 1}`);
    });
    delegatedFeatureHandlers.pairwiseAffordancesPrepared = true;
    return true;
  };

  return {
    applyVisibilityPreviewBySvgId,
    applyVisibilityPreviewChanges,
    attachSvgFeatureHandlers,
    getFeatureElements,
    getFeatureFillElements,
    openFeatureEditorForFeature,
    previewAlignmentCandidate,
    clearAlignmentCandidatePreview,
    showAlignmentOverlay,
    clearAlignmentOverlay,
    preparePairwiseInteractionAffordances,
    dispose
  };
};
