// @ts-check
/** @import { DrawingState, PopupPosition } from '../../state.js' */
/** @import { PopupOrigin } from '../../utils/popup-bounds.js' */
/** @import { AliasSource, FeatureIdentity } from '../../services/feature-identity.js' */
/** @import { FeatureLike } from '../../services/feature-utils.js' */
/** @import { AlignmentFigureGeometry, AlignmentOverlayCandidate, AlignmentOverlayRecord, FigureSpan } from '../../services/alignment-review-view.js' */
import { resolveColorToHex } from '../../utils/color-utils.js';
import { placePopupNearPoint, popupViewportBounds } from '../../utils/popup-bounds.js';
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
import { COMPARISON_LEGEND_SELECTOR, legendEntryKey } from '../../services/legend-svg.js';
import { PREVIEW_ITEM_SELECTOR, previewItemTarget } from '../preview-item-popup.js';
import { recordStructuralMetric } from '../../services/runtime-test-hooks.js';
import {
  featureIdentity,
  indexFeatureIdentities,
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

// The hover highlight of a feature or band is written once the pointer has
// rested on it this long (OV-900, D-01): a pointer passing over bands writes
// none, so the browser re-rasters no band area. It comes before the hover
// summary (180 ms).
export const HOVER_HIGHLIGHT_DWELL_MS = 70;

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
 * A pointer position, as an event or a plain point carries it.
 * @typedef {{ clientX?: number, clientY?: number }} HoverPoint
 */
/**
 * A position that names the element under it, which a mouse or pointer event does.
 * @typedef {{ clientX: number, clientY: number, target?: EventTarget | null }} PointerPosition
 */
/**
 * An event target that may be an Element (the window and a document are not).
 * @typedef {EventTarget & { closest?: Element['closest'] }} MaybeElementTarget
 */

/**
 * The hover summary of a pairwise match, as the match popup builds it.
 * @typedef {{ id: string, title: string, subtitle: string, fill: string, rows: Array<{ label: string, value: string }> }} MatchHoverSummary
 */

/**
 * The feature selection owner's reactions to a click or a drag on the preview.
 * @typedef {object} FeatureSelectionPort
 * @property {() => boolean} [consumeSuppressNextClick]
 * @property {(event: Event, svg: Element) => Element | null} [getSelectableFeatureTarget]
 * @property {(selectionId: string) => any} [toggleFeatureSelection]
 * @property {(selectionId: string, options?: { additive?: boolean }) => any} [selectFeatureRange]
 * @property {(options?: { clearStatus?: boolean }) => any} [clearFeatureSelection]
 * @property {(svgId: string) => any} [markPlainFeatureClick]
 * @property {(event: PointerEvent, svg: Element) => boolean} [startMarqueePointer]
 * @property {(event: PointerEvent) => boolean} [moveMarqueePointer]
 * @property {(event: PointerEvent) => boolean} [commitMarqueePointer]
 * @property {() => any} [cancelMarqueePointer]
 */

/**
 * The preview owner's pan and zoom interaction, which pauses feature hover.
 * @typedef {object} PreviewTransformInteractionPort
 * @property {() => boolean} [isActive]
 * @property {(listener: (change: { active: boolean, kind: string, event: Event | null, reconcile: boolean }) => void) => (() => void) | void} [subscribe]
 */

/**
 * The visibility owner's preparation of the rule matches `resolveFeatureDrawn`
 * reads. Strict, it rejects when the color preparation fails and resolves to
 * false when it is stale.
 * @typedef {(options?: { strict?: boolean, view?: boolean }) => boolean | Promise<boolean | { error: any }>} PrepareDrawnFeatureMatchesPort
 */

/**
 * The name the Legend draws for the row a caption is the key of (the Legend owner's `shownLegendName`).
 * @typedef {(caption: string) => string} ShownLegendNamePort
 */

/**
 * @typedef {object} FeatureSvgActionsOptions
 * @property {Record<string, any>} state App state (state.js; not yet typed).
 * @property {(feature: Record<string, any>) => string} getFeatureColor The rule owner's color of a feature.
 * @property {(feature: Record<string, any>) => string} [getEffectiveLegendCaption] The rule owner's Legend caption of a feature.
 * @property {ShownLegendNamePort} [shownLegendName]
 * @property {(() => void) | null} [onFeaturePopupOpened] The label owner's sync of its editor when the popup opens.
 * @property {PrepareDrawnFeatureMatchesPort} [prepareDrawnFeatureMatches]
 *   The visibility owner's preparation of the rule matches `resolveFeatureDrawn` reads.
 * @property {FeatureSelectionPort | null} [featureSelection]
 *   The preview owner's projection of feature visibility changes.
 * @property {PreviewTransformInteractionPort | null} [previewTransformInteraction]
 * @property {() => boolean} [previewPopupsHeld] While it holds (the Review alignment panel is open, the popups
 *   hidden), no gesture on the figure opens, retargets or closes a popup.
 */

/**
 * One row of the hover summary: its label and value text nodes.
 * @typedef {{ row: HTMLElement, label: Text, valueElement: HTMLElement, value: Text }} HoverSummaryRowView
 */

/**
 * The hover summary's DOM, built once (OV-933): a show rewrites its text nodes,
 * the swatch colour and which rows are shown, so it adds or removes no node and
 * changes a class only when a row's clamping changes. The vendored Tailwind
 * runtime rescans every class on a node or class change.
 * @typedef {object} HoverSummaryView
 * @property {HTMLElement} element
 * @property {HTMLElement} swatch
 * @property {Text} heading
 * @property {Text} subtitle
 * @property {HoverSummaryRowView[]} rows
 */

/**
 * What one hover summary shows: a feature's or a match's.
 * @typedef {object} HoverSummaryContent
 * @property {string} color
 * @property {string} heading
 * @property {string} subtitle
 * @property {Array<{ label: string, value: unknown, clamp: boolean }>} rows
 */

/**
 * The state of the feature and match hover summary. The tooltip is built on the
 * first show; `lastEvent` holds the last pointer position seen.
 * @typedef {object} HoverSummaryState
 * @property {HoverSummaryView | null} view
 * @property {number | null} timer
 * @property {number | null} frame
 * @property {boolean} visible
 * @property {{ clientX?: number, clientY?: number } | null} lastEvent
 */

/**
 * @typedef {object} AlignmentOverlayRequest
 * @property {{ key: string, recordKey: string, anchor: AliasSource, name: string }} reference
 * @property {AlignmentOverlayRecord[]} records
 * @property {(recordKey: string, anchor: AliasSource) => void} onSelect
 * @property {(recordKey: string, candidateKey: string) => void} onHover
 */

/**
 * A candidate the Result draws, with its blocks and its ring button.
 * @typedef {AlignmentOverlayCandidate & { record: AlignmentOverlayRecord, svgId: string, blocks: Element[], button: HTMLButtonElement }} AlignmentOverlayRing
 */

/**
 * The Similarity alignment overlay over the preview viewport. `resizeObserver` stays null
 * without ResizeObserver.
 * @typedef {object} AlignmentOverlay
 * @property {HTMLElement} element
 * @property {HTMLElement} scrollRoot
 * @property {number} frame
 * @property {ResizeObserver | null} resizeObserver
 * @property {MutationObserver} transformObserver
 * @property {AlignmentOverlayRing | null} hover
 * @property {() => void} schedule
 * @property {(ring: AlignmentOverlayRing | null) => void} setHover
 * @property {(recordKey: string, anchor: AliasSource) => void} onSelect
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
 * @property {((reference: FeatureIdentity) => RenderedIdentityEntry[]) | null} renderedFeatureIdentities
 * @property {Set<Element>} selectedMatchElements The match elements marked as the open match popup's.
 * @property {boolean} pairwiseAffordancesPrepared
 * @property {string | null} activeHoverSvgId
 * @property {string} activeHoverKey
 * @property {Element | null} activeMatchHoverElement
 * @property {string} activeMatchHoverKey
 * @property {Element | null} pendingMatchElement
 * @property {Set<Element>} activeHoverElements
 * @property {string} hoverHighlightKey The hover whose highlight is written ('' for none).
 * @property {string} pendingHoverHighlightKey The hover whose highlight waits for the dwell.
 * @property {number | null} hoverHighlightTimer
 * @property {{ clientX: number, clientY: number } | null} transformPointer
 * @property {boolean} reconcileHoverAfterTransform
 * @property {number | null} hoverReconcileFrame
 * @property {string} alignmentCandidateSvgId
 * @property {{ featureSvgId: string | null, matchElement: Element | null } | null} alignmentCandidateRestore
 * @property {AlignmentOverlay | null} alignmentOverlay
 * @property {Map<string, AlignmentOverlayRing>} alignmentCandidatesBySvgId
 * @property {((request: AlignmentOverlayRequest) => AlignmentFigureGeometry | null) | null} showAlignmentOverlay
 * @property {(() => void) | null} clearAlignmentOverlay
 * @property {((anchor: AliasSource) => boolean) | null} previewAlignmentCandidate
 * @property {((options?: { restore?: boolean }) => void) | null} clearAlignmentCandidatePreview
 * @property {((eventLike: any, kind?: string) => void) | null} beginPreviewTransformInteraction
 * @property {((options?: { reconcile?: boolean }) => void) | null} endPreviewTransformInteraction
 * @property {((match: SelectedPairwiseMatch | null) => void) | null} markSelectedMatch
 * @property {(() => void) | null} cleanup
 */

/** @typedef {{ svgId: string, identity: FeatureIdentity }} RenderedIdentityEntry */

/**
 * The open match popup's payload, as far as it names the bands it shows.
 * @typedef {{ id?: unknown, matchKind?: unknown, orthogroupId?: unknown }} SelectedPairwiseMatch
 */

/** @param {FeatureSvgActionsOptions} options */
export const createFeatureSvgActions = ({
  state,
  getFeatureColor,
  getEffectiveLegendCaption,
  shownLegendName = (caption) => caption,
  onFeaturePopupOpened = null,
  // R13: the popup opens once the rule matches `resolveFeatureDrawn` reads are
  // prepared; the composition root injects the visibility owner's preparation.
  prepareDrawnFeatureMatches = () => true,
  featureSelection = null,
  previewTransformInteraction = null,
  previewPopupsHeld = () => false
}) => {
  const {
    orthogroups,
    extractedFeatures,
    biologicalFeatures,
    featuresBySvgId,
    svgContainer,
    clickedFeature,
    clickedFeaturePos,
    clickedPairwiseMatch,
    clickedPairwiseMatchPos,
    matchSequenceRegistry,
    selectedAnnotation,
    featurePopupSize,
    featureSelectionDrag
  } = state;
  /** @type {DelegatedFeatureHandlers | null} */
  let delegatedFeatureHandlers = null;
  const isPreviewTransformInteractionActive = () => Boolean(
    previewTransformInteraction?.isActive?.()
  );
  /** @type {HoverSummaryState} */
  const hoverSummaryState = {
    view: null,
    timer: null,
    frame: null,
    visible: false,
    lastEvent: null
  };
  // The fine-pointer query, asked once; its `matches` follows the device (OV-938).
  /** @type {MediaQueryList | null | undefined} */
  let finePointerQuery;

  /** @param {unknown} value */
  const getOrthogroupIds = (value) =>
    Array.from(new Set(
      String(value || '')
        .split(';')
        .map((entry) => entry.trim())
        .filter(Boolean)
    ));

  /** @param {FeatureLike | null | undefined} feature */
  const renderedFeatureSvgId = (feature) => String(
    feature?.rendered_svg_id ||
    feature?.renderedSvgId ||
    feature?.rendered_feature_svg_id ||
    feature?.renderedFeatureSvgId ||
    feature?.svg_id ||
    ''
  ).trim();

  /** @param {unknown} value */
  const normalizeVisibilityMode = (value) => {
    const normalized = String(value || '').trim().toLowerCase();
    if (normalized === 'suppress') return 'exclude_matching';
    return ['on', 'off', 'exclude_matching'].includes(normalized) ? normalized : 'default';
  };

  // OV-371: a popup opened at a clicked point is placed beside it with an
  // expected size; app-setup places it again with its rendered size. `avoid`
  // is the clicked item's rectangle, which the popup leaves uncovered; the
  // origin's `panel` is the panel it opens beside (OV-503).
  /**
   * @param {PopupPosition} pos
   * @param {PopupOrigin | null | undefined} eventLike
   * @param {number} width
   * @param {DOMRect | null} [avoid]
   */
  const placePopupAt = (pos, eventLike, width, avoid = null) => {
    const anchor = Number.isFinite(eventLike?.clientX) && Number.isFinite(eventLike?.clientY)
      ? {
          x: Number(eventLike?.clientX), y: Number(eventLike?.clientY),
          avoid: avoid && { left: avoid.left, top: avoid.top, right: avoid.right, bottom: avoid.bottom },
          panel: eventLike?.panel
        }
      : null;
    Object.assign(pos, placePopupNearPoint({ point: anchor, width, height: 520, bounds: popupViewportBounds() }), { anchor });
  };

  /** @param {unknown} qualifiers */
  const normalizeQualifierRows = (qualifiers) => {
    if (!qualifiers || typeof qualifiers !== 'object' || Array.isArray(qualifiers)) return [];
    return Object.entries(qualifiers)
      .map(([key, value]) => {
        const values = normalizeStringArray(value);
        return { key: String(key || ''), values };
      })
      .filter((row) => row.key && row.values.length > 0)
      .sort((left, right) => left.key.localeCompare(right.key));
  };

  /** @param {FeatureLike | null | undefined} feat @param {unknown} key */
  const getQualifierFirstValue = (feat, key) => {
    const normalizedKey = String(key || '').trim().toLowerCase();
    if (!normalizedKey) return '';
    const directValue = feat?.[normalizedKey];
    // A non-null object (checked): the qualifier table is read by a fixed key.
    const qualifierValue = feat?.qualifiers && typeof feat.qualifiers === 'object'
      ? /** @type {Record<string, unknown>} */ (feat.qualifiers)[normalizedKey]
      : null;
    const values = normalizeStringArray(directValue || qualifierValue)
      .map((value) => String(value || '').trim())
      .filter(Boolean);
    return values[0] || '';
  };

  /** @param {FeatureLike | null | undefined} feat */
  const getHoverSummaryPrimaryLabel = (feat) => (
    getQualifierFirstValue(feat, 'gene') ||
    getQualifierFirstValue(feat, 'locus_tag') ||
    getQualifierFirstValue(feat, 'product') ||
    getFeatureCaption(feat) ||
    ''
  );

  /** @param {string} tagName @param {string} [className] */
  const createHoverSummaryElement = (tagName, className = '') => {
    const element = document.createElement(tagName);
    if (className) element.className = className;
    return element;
  };

  /** @param {HTMLElement} parent @param {string} className */
  const appendHoverSummaryText = (parent, className) => {
    const element = createHoverSummaryElement('div', className);
    const text = document.createTextNode('');
    element.appendChild(text);
    parent.appendChild(element);
    return { element, text };
  };

  /** @param {FeatureLike} feat @param {string} primaryLabel */
  const buildHoverSummaryRows = (feat, primaryLabel) => {
    const product = getQualifierFirstValue(feat, 'product');
    const gene = getQualifierFirstValue(feat, 'gene');
    const locusTag = getQualifierFirstValue(feat, 'locus_tag');
    const note = getQualifierFirstValue(feat, 'note');
    const locationText = formatFeatureLocation(feat);
    // The feature's row as the Legend draws it; its caption is the row's key, which is never shown.
    const legendName = shownLegendName(String(getEffectiveLegendCaption?.(feat) || '').trim());
    /** @type {Array<[string, unknown, boolean?]>} */
    const rows = [];

    if (gene && gene !== primaryLabel) rows.push(['Gene', gene]);
    if (locusTag && locusTag !== primaryLabel) rows.push(['Locus', locusTag]);
    if (product && product !== primaryLabel) rows.push(['Product', product, true]);
    if (note && note !== primaryLabel && note !== product) rows.push(['Note', note, true]);
    rows.push(['Length', formatFeatureLength(feat)]);
    rows.push(['Location', locationText]);
    rows.push(['Record', feat?.record_id || '']);
    if (feat?.orthogroupId) rows.push(['Similarity group', feat.orthogroupId]);
    if (legendName && legendName !== primaryLabel) rows.push(['Legend', legendName]);
    return rows;
  };

  // Details lists what the popup header and the Edit tab do not show:
  // the header names the feature, its product, record, and location, and the
  // Similarity group section its group.
  /** @param {FeatureLike} feat */
  const buildDetailRows = (feat) => [
    { key: 'type', label: 'Feature type', value: feat.type },
    { key: 'length', label: 'Length', value: formatFeatureLength(feat) },
    { key: 'protein_id', label: 'Protein ID', value: resolveFeatureProteinId(feat, feat?.orthogroupMember || feat?.orthogroup_member || null) }
  ]
    .map((row) => ({ ...row, value: row.value === null || row.value === undefined ? '' : String(row.value) }))
    .filter((row) => row.value !== '');

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

  /** @param {MaybeElementTarget | null | undefined} target @param {Element} svg */
  const getFeatureTarget = (target, svg) => {
    if (!target || typeof target.closest !== 'function') return null;
    const matchEl = target.closest(PAIRWISE_MATCH_SELECTOR);
    if (matchEl && svg.contains(matchEl)) return null;
    const featureEl = target.closest(FEATURE_SELECTOR);
    if (!featureEl || !svg.contains(featureEl)) return null;
    return featureEl;
  };

  /** @param {MaybeElementTarget | null | undefined} target @param {Element} svg */
  const getPairwiseMatchTarget = (target, svg) => {
    if (!target || typeof target.closest !== 'function') return null;
    const matchEl = target.closest(PAIRWISE_MATCH_SELECTOR);
    if (!matchEl || !svg.contains(matchEl)) return null;
    return matchEl;
  };

  /** @param {PointerPosition | null | undefined} eventLike @param {Element | null | undefined} svg @param {string} selector */
  const getTopmostSvgTarget = (eventLike, svg, selector) => {
    if (!svg || !selector || !Number.isFinite(eventLike?.clientX) || !Number.isFinite(eventLike?.clientY)) {
      return null;
    }
    const stack = typeof document.elementsFromPoint === 'function'
      // `Number.isFinite(eventLike?.clientX)` holds only for a present eventLike.
      ? document.elementsFromPoint(/** @type {PointerPosition} */ (eventLike).clientX, /** @type {PointerPosition} */ (eventLike).clientY)
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

  /** @param {PointerPosition} eventLike @param {Element} svg */
  const getFeatureClickTarget = (eventLike, svg) =>
    getTopmostSvgTarget(eventLike, svg, FEATURE_SELECTOR) || getFeatureTarget(eventLike?.target, svg);

  /** @param {PointerPosition} eventLike @param {Element} svg */
  const getPairwiseMatchClickTarget = (eventLike, svg) =>
    getTopmostSvgTarget(eventLike, svg, PAIRWISE_MATCH_SELECTOR) || getPairwiseMatchTarget(eventLike?.target, svg);

  /** @param {Element | null | undefined} target @param {Element | null | undefined} svg */
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
   * @param {DrawingState} drawing
   * @param {Record<string, any>} feat
   * @param {Element | null} [featureElement]
   * @param {string} [renderedSvgId]
   */
  const buildClickedFeaturePayload = (drawing, feat, featureElement = null, renderedSvgId = undefined) => {
    const defaultLabel = getFeatureCaption(feat);
    const existingOverride = getFeatureOverride(drawing.featureColorOverrides, feat);
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
    const visibilityMode = normalizeVisibilityMode(featureOverrideValue(drawing.featureOverrides, feat, 'featureVisibility'));

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
      product: getQualifierFirstValue(feat, 'product'),
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
      detailRows: buildDetailRows(feat),
      legendName: effectiveCaption,
      appliedLegendName: effectiveCaption,
      strokeColor: currentStrokeColor,
      strokeWidth: currentStrokeWidth,
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
   * @param {DrawingState} drawing
   * @param {Record<string, any>} feat
   * @param {PopupOrigin | null} [eventLike]
   */
  const openPreparedFeatureEditor = (drawing, feat, eventLike = null) => {
    if (!feat) return null;
    if (!svgContainer.value) return null;
    const svg = svgContainer.value.querySelector('svg');
    if (!svg) return null;

    // A catalog feature opens as the displayed Result draws it; one that the
    // Result lists but does not draw (a hidden feature, R-5) opens without an
    // element, so its popup can show it again.
    // The catalog lookup below can miss; `if (!target) return null` follows it.
    /** @type {typeof feat | undefined} */
    let target = feat;
    let renderedSvgId = renderedFeatureSvgId(feat);
    const catalogFeatures = resultCatalogFeatures(state);
    const key = stableFeatureOverrideKey(feat);
    if (catalogFeatures && key) {
      target = catalogFeatures.renderedByIdentity.get(key)
        || catalogFeatures.biological.find((/** @type {FeatureLike} */ candidate) => stableFeatureOverrideKey(candidate) === key);
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
    clickedFeature.value = buildClickedFeaturePayload(drawing, target, featureElement, renderedSvgId);
    if (featurePopupSize) {
      featurePopupSize.width = 0;
      featurePopupSize.height = 0;
    }

    placePopupAt(clickedFeaturePos, eventLike, state.richFeaturePopup.value === false ? 440 : 720);
    if (typeof onFeaturePopupOpened === 'function') {
      onFeaturePopupOpened();
    }
    return clickedFeature.value;
  };

  // The popup states whether the feature is drawn (resolveFeatureDrawn).
  /**
   * @param {Record<string, any>} feat
   * @param {PopupOrigin | null} [eventLike]
   */
  const openFeatureEditorForFeature = (feat, eventLike = null) => {
    const drawing = state.activeDrawing();
    return reportRuleRunFailure(
      state, 'feature-extraction', () => runWhenPrepared(
        state, () => [prepareDrawnFeatureMatches({ strict: true, view: true })], () => openPreparedFeatureEditor(drawing, feat, eventLike)
      )
    );
  };

  // A popup retargeted while it was hidden was placed for the view of that moment. Once the
  // Result rests in the view it returned to, the popup opens again beside its feature.
  const showFeaturePopupAtFeature = async () => {
    const shown = clickedFeature.value;
    const surface = /** @type {HTMLElement | null} */ (svgContainer.value);
    if (!shown || !surface) return;
    await new Promise((resolve) => { window.requestAnimationFrame(resolve); });
    await Promise.allSettled(surface.getAnimations().map(({ finished }) => finished));
    if (clickedFeature.value !== shown) return;
    const rects = getFeatureFillElements(surface.querySelector('svg'), shown.svg_id)
      .map((block) => block.getBoundingClientRect());
    if (!rects.length) return;
    openFeatureEditorForFeature(shown.feat, {
      clientX: (Math.min(...rects.map(({ left }) => left)) + Math.max(...rects.map(({ right }) => right))) / 2,
      clientY: (Math.min(...rects.map(({ top }) => top)) + Math.max(...rects.map(({ bottom }) => bottom))) / 2
    });
  };

  const hoverSummaryIsAllowed = () => {
    if (clickedFeature.value) return false;
    if (clickedPairwiseMatch?.value) return false;
    if (isPreviewTransformInteractionActive()) return false;
    if (finePointerQuery === undefined) {
      finePointerQuery = window.matchMedia ? window.matchMedia('(hover: hover) and (pointer: fine)') : null;
    }
    if (finePointerQuery && !finePointerQuery.matches) {
      return false;
    }
    return true;
  };

  /** @returns {HoverSummaryView} */
  const ensureHoverSummaryView = () => {
    if (hoverSummaryState.view?.element.isConnected) return hoverSummaryState.view;
    const element = createHoverSummaryElement('div', 'feature-hover-summary');
    element.hidden = true;
    element.setAttribute('role', 'tooltip');
    const title = createHoverSummaryElement('div', 'feature-hover-summary-title');
    const swatch = createHoverSummaryElement('div', 'feature-hover-summary-swatch');
    const titleTextWrap = createHoverSummaryElement('div', 'feature-hover-summary-text');
    const heading = appendHoverSummaryText(titleTextWrap, 'feature-hover-summary-heading').text;
    const subtitle = appendHoverSummaryText(titleTextWrap, 'feature-hover-summary-subtitle').text;
    title.appendChild(swatch);
    title.appendChild(titleTextWrap);
    element.appendChild(title);
    document.body.appendChild(element);
    hoverSummaryState.view = { element, swatch, heading, subtitle, rows: [] };
    return hoverSummaryState.view;
  };

  /** @param {HTMLElement} container @returns {HoverSummaryRowView} */
  const appendHoverSummaryRow = (container) => {
    const row = createHoverSummaryElement('div', 'feature-hover-summary-row');
    const label = appendHoverSummaryText(row, 'feature-hover-summary-label').text;
    const { element: valueElement, text: value } = appendHoverSummaryText(row, 'feature-hover-summary-value');
    container.appendChild(row);
    return { row, label, valueElement, value };
  };

  const positionHoverSummary = (eventLike = hoverSummaryState.lastEvent) => {
    const element = hoverSummaryState.view?.element;
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

  /** @param {HoverPoint | null | undefined} eventLike */
  const scheduleHoverSummaryPosition = (eventLike) => {
    // `Number.isFinite(eventLike?.clientX)` holds only for a present eventLike.
    hoverSummaryState.lastEvent = {
      clientX: Number.isFinite(eventLike?.clientX) ? /** @type {HoverPoint} */ (eventLike).clientX : hoverSummaryState.lastEvent?.clientX,
      clientY: Number.isFinite(eventLike?.clientY) ? /** @type {HoverPoint} */ (eventLike).clientY : hoverSummaryState.lastEvent?.clientY
    };
    if (!hoverSummaryState.visible || hoverSummaryState.frame) return;
    hoverSummaryState.frame = window.requestAnimationFrame(() => {
      hoverSummaryState.frame = null;
      positionHoverSummary();
    });
  };

  /**
   * @param {FeatureLike} feat
   * @param {Element | null | undefined} featureElement
   * @returns {HoverSummaryContent}
   */
  const featureHoverContent = (feat, featureElement) => {
    const primaryLabel = getHoverSummaryPrimaryLabel(feat);
    const featureType = String(feat?.type || 'Feature').trim() || 'Feature';
    return {
      color: resolveColorToHex(
        featureElement?.getAttribute?.('fill') || getFeatureColor(feat) || '#94a3b8'
      ) || '#94a3b8',
      heading: primaryLabel ? `${featureType}: ${primaryLabel}` : featureType,
      subtitle: formatFeatureLocation(feat),
      rows: buildHoverSummaryRows(feat, primaryLabel)
        .map(([label, value, clamp]) => ({ label, value, clamp: Boolean(clamp) }))
    };
  };

  /**
   * @param {MatchHoverSummary | null} summary
   * @returns {HoverSummaryContent | null}
   */
  const matchHoverContent = (summary) => summary && {
    color: resolveColorToHex(summary.fill || '#94a3b8') || '#94a3b8',
    heading: summary.title,
    subtitle: summary.subtitle || summary.id || '',
    rows: summary.rows.map((row) => ({ ...row, clamp: row.label === 'Query' || row.label === 'Subject' }))
  };

  // The rows without a value are not shown; the view keeps its other rows
  // empty and hidden.
  /** @param {HoverSummaryContent | null} content @param {HoverPoint | null | undefined} eventLike */
  const renderHoverSummary = (content, eventLike) => {
    if (!content || !hoverSummaryIsAllowed()) {
      hideHoverSummary();
      return;
    }
    const view = ensureHoverSummaryView();
    view.swatch.style.backgroundColor = content.color;
    view.heading.data = content.heading;
    view.subtitle.data = content.subtitle;
    const rows = content.rows
      .map((row) => ({ ...row, value: String(row.value === null || row.value === undefined ? '' : row.value).trim() }))
      .filter((row) => row.value);
    while (view.rows.length < rows.length) view.rows.push(appendHoverSummaryRow(view.element));
    view.rows.forEach((rowView, index) => {
      const row = rows[index];
      rowView.row.style.display = row ? '' : 'none';
      rowView.label.data = row?.label || '';
      rowView.value.data = row?.value || '';
      if (row) rowView.valueElement.classList.toggle('is-clamped', row.clamp);
    });
    view.element.hidden = false;
    hoverSummaryState.visible = true;
    // Placed in the next frame, before it paints (OV-938): measuring here would
    // force a style pass right after the highlight writes.
    scheduleHoverSummaryPosition(eventLike);
  };

  // The summary is built when it shows, 180 ms after the pointer rests or at
  // once while one is shown, not on every pointer entry (OV-915).
  /** @param {(() => HoverSummaryContent | null) | null} buildContent @param {HoverPoint | null | undefined} eventLike */
  const scheduleHoverSummary = (buildContent, eventLike) => {
    if (hoverSummaryState.timer) {
      window.clearTimeout(hoverSummaryState.timer);
      hoverSummaryState.timer = null;
    }
    if (!buildContent || !hoverSummaryIsAllowed()) {
      hideHoverSummary();
      return;
    }
    scheduleHoverSummaryPosition(eventLike);
    const show = () => {
      hoverSummaryState.timer = null;
      renderHoverSummary(buildContent(), eventLike);
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
    if (hoverSummaryState.view) {
      hoverSummaryState.view.element.hidden = true;
    }
    hoverSummaryState.visible = false;
    hoverSummaryState.lastEvent = null;
  }

  // Every catalog group, collinear presentation groups included, is in
  // `orthogroups`; the popup's scoped lookup picks the one a ribbon names. A
  // collinear ribbon without a collinear group scope names ordinary similarity
  // groups, which the presentation subset lacks (OV-370).
  const groupsForMatch = () => orthogroups?.value || [];

  /** @param {DrawingState} drawing @param {Element} matchElement @param {Map<string, FeatureLike>} featureLookup */
  const buildMatchPayload = (drawing, matchElement, featureLookup) => buildMatchPopupPayload(/** @type {SVGElement} */ (matchElement), {
    featureLookup,
    sourceFeatures: Array.isArray(biologicalFeatures?.value) && biologicalFeatures.value.length > 0
      ? biologicalFeatures.value
      : (Array.isArray(extractedFeatures.value) ? extractedFeatures.value : []),
    orthogroups: groupsForMatch(),
    orthogroupNameOverrides: drawing.orthogroupNameOverrides,
    orthogroupDescriptionOverrides: drawing.orthogroupDescriptionOverrides,
    resolveSequenceSource: matchSequenceRegistry?.resolve
  });

  // The summary builder reads a function as well as a list.
  /** @param {Element} matchElement */
  const buildMatchHoverSummary = (matchElement) => buildPairwiseMatchHoverSummary(/** @type {SVGElement} */ (matchElement), /** @type {any} */ ({
    orthogroups: groupsForMatch,
    orthogroupNameOverrides: state.activeDrawing().orthogroupNameOverrides
  }));

  /** @param {Element | null} matchElement @param {PopupOrigin | null} eventLike @param {Map<string, FeatureLike>} featureLookup */
  const openPairwiseMatchPopup = (matchElement, eventLike, featureLookup) => {
    const drawing = state.activeDrawing();
    if (!matchElement || !clickedPairwiseMatch || !clickedPairwiseMatchPos || previewPopupsHeld()) return null;
    const payload = buildMatchPayload(drawing, matchElement, featureLookup);
    if (!payload) return null;
    hideHoverSummary();
    clickedFeature.value = null;
    clickedPairwiseMatch.value = payload;
    placePopupAt(clickedPairwiseMatchPos, eventLike, 460);
    return payload;
  };

  // A Legend row or a record definition opens its popup in the match popup
  // frame while Layout edit is off; with it on, a press drags them. A Legend
  // row opens only when the drawing lists it (not a track-data row).
  /** @param {{ element: Element, item: import('../preview-item-popup.js').PreviewItem }} target @param {{ clientX?: number, clientY?: number }} eventLike */
  const openPreviewItemPopup = ({ element, item }, eventLike) => {
    if (!clickedPairwiseMatch || !clickedPairwiseMatchPos || state.layoutRepositionMode?.value || previewPopupsHeld()) return false;
    if (item.kind === 'legend-entry') {
      /** @type {import('../../state.js').LegendEntry[]} */
      const entries = state.activeDrawing().legendEntries.value || [];
      const { key } = item;
      const entry = entries.find((row) => legendEntryKey(row) === key || String(row.caption) === key);
      if (!entry) return false;
      // The popup follows the row by its key, which a rename keeps.
      item = { kind: 'legend-entry', key: legendEntryKey(entry) };
    }
    const rect = element.getBoundingClientRect();
    const anchor = Number.isFinite(eventLike?.clientX)
      ? eventLike
      : { clientX: rect.left + rect.width / 2, clientY: rect.bottom };
    hideHoverSummary();
    clickedFeature.value = null;
    // A Legend row's feature count reads the rule matches: the popup opens at
    // once and counts when they are prepared.
    const opened = { previewItem: item, matchesReady: item.kind !== 'legend-entry' };
    clickedPairwiseMatch.value = opened;
    // The preview-item popup is 400 px wide (`.preview-item-popup`).
    placePopupAt(clickedPairwiseMatchPos, anchor, 400, rect);
    if (item.kind === 'legend-entry') {
      const openedKey = item.key;
      reportRuleRunFailure(state, 'feature-extraction', () => runWhenPrepared(
        state, () => [prepareDrawnFeatureMatches({ strict: true, view: true })], () => {
          const current = clickedPairwiseMatch.value;
          if (current?.previewItem?.kind === 'legend-entry' && current.previewItem.key === openedKey) current.matchesReady = true;
        }
      ));
    }
    return true;
  };

  /** @param {{ root?: SVGSVGElement | null, phase?: string, rootGeneration?: number }} [options] */
  const attachSvgFeatureHandlers = ({
    root = null,
    phase = 'preview-bind',
    rootGeneration = 0
  } = {}) => {
    /** @type {SVGSVGElement | null} */
    const svg = root || svgContainer.value?.querySelector?.('svg') || null;
    if (!svg) return false;

    if (delegatedFeatureHandlers && delegatedFeatureHandlers.svg !== svg) {
      cleanupDelegatedFeatureHandlers();
    }

    /** @param {Map<string, FeatureLike>} featureLookup */
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
      renderedFeatureIdentities: null,
      selectedMatchElements: new Set(),
      pairwiseAffordancesPrepared: false,
      activeHoverSvgId: null,
      activeHoverKey: '',
      activeMatchHoverElement: null,
      activeMatchHoverKey: '',
      pendingMatchElement: null,
      activeHoverElements: new Set(),
      hoverHighlightKey: '',
      pendingHoverHighlightKey: '',
      hoverHighlightTimer: null,
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
      markSelectedMatch: null,
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

      /** @param {(Element & { style?: CSSStyleDeclaration }) | null | undefined} element @param {boolean} highlight */
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

      /** @param {string} svgId */
      const setFeatureHover = (svgId) => {
        (ensureFeaturePaths().get(svgId) || []).forEach((element) => {
          setHoverStyle(element, true);
        });
      };

      /** @param {string} svgId */
      const getFeatureHoverKey = (svgId) => {
        const feat = ensureFeatureLookup().get(svgId);
        const orthogroupId = String(feat?.orthogroupId || '').trim();
        return orthogroupId ? `orthogroup:${orthogroupId}` : `feature:${svgId}`;
      };

      /** @param {unknown} orthogroupId */
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

      /** @param {unknown} blockId */
      const setCollinearityBlockHover = (blockId) => {
        const id = String(blockId || '').trim();
        if (!id) return;
        ensureComparisonIndexes();
        // ensureComparisonIndexes above has set all three comparison maps.
        (/** @type {Map<string, Element[]>} */ (handlerState.comparisonElementsByCollinearityBlockId).get(id) || []).forEach((element) => {
          setHoverStyle(element, true);
        });
      };

      /** @param {string} svgId */
      const setHoverHighlight = (svgId) => {
        const feat = ensureFeatureLookup().get(svgId);
        const orthogroupId = String(feat?.orthogroupId || '').trim();
        if (orthogroupId) {
          setOrthogroupHover(orthogroupId);
          return;
        }
        setFeatureHover(svgId);
      };

      /** @param {Element | null | undefined} element @param {string} name */
      const matchAttr = (element, name) => String(element?.getAttribute?.(name) || '').trim();

      // A similarity-group ribbon, as opposed to a pairwise match or a collinear
      // block that names groups.
      /** @param {Element} element */
      const isOrthogroupRibbon = (element) => {
        const kind = matchAttr(element, 'data-match-kind').toLowerCase();
        return kind === 'orthogroup' || (
          !kind && !element.hasAttribute('data-collinearity-block-id') && Boolean(matchAttr(element, 'data-orthogroup-id'))
        );
      };

      /** @param {Element | null | undefined} matchElement */
      const getMatchHoverKey = (matchElement) => {
        const blockId = matchAttr(matchElement, 'data-collinearity-block-id');
        if (blockId) return `collinearity:${blockId}`;
        const orthogroupId = matchAttr(matchElement, 'data-orthogroup-id');
        if (orthogroupId) return `orthogroup:${orthogroupId}`;
        return `match:${matchAttr(matchElement, 'data-gbdraw-match-id') || matchAttr(matchElement, 'data-gbdraw-pairwise-match-id') || matchAttr(matchElement, 'd')}`;
      };

      /** @param {Element} element */
      const matchFragments = (element) => {
        ensureComparisonIndexes();
        const id = matchAttr(element, 'data-gbdraw-match-id') || matchAttr(element, 'data-gbdraw-pairwise-match-id');
        // ensureComparisonIndexes above has set all three comparison maps.
        return /** @type {Map<string, Element[]>} */ (handlerState.comparisonElementsByMatchId).get(id) || [element];
      };

      /** @param {Element | null | undefined} matchElement */
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

      // The hover highlight (OV-900): one hover's styles are written at a time,
      // once the pointer has rested on it for HOVER_HIGHLIGHT_DWELL_MS. Until
      // then the written highlight stays, so a move from one target straight
      // to the next swaps the two in one step and never shows neither.
      const cancelPendingHoverHighlight = () => {
        if (handlerState.hoverHighlightTimer !== null) window.clearTimeout(handlerState.hoverHighlightTimer);
        handlerState.hoverHighlightTimer = null;
        handlerState.pendingHoverHighlightKey = '';
      };

      /** @param {string} key @param {() => void} write */
      const showHoverHighlight = (key, write) => {
        cancelPendingHoverHighlight();
        if (handlerState.hoverHighlightKey) clearTrackedHoverStyles();
        write();
        handlerState.hoverHighlightKey = key;
      };

      /** @param {string} key @param {() => void} write */
      const requestHoverHighlight = (key, write) => {
        if (key === handlerState.pendingHoverHighlightKey) return;
        cancelPendingHoverHighlight();
        if (key === handlerState.hoverHighlightKey) return;
        handlerState.pendingHoverHighlightKey = key;
        handlerState.hoverHighlightTimer = window.setTimeout(() => {
          handlerState.hoverHighlightTimer = null;
          showHoverHighlight(key, write);
        }, HOVER_HIGHLIGHT_DWELL_MS);
      };

      const clearHoverHighlight = () => {
        cancelPendingHoverHighlight();
        clearTrackedHoverStyles();
        handlerState.hoverHighlightKey = '';
      };

      // The hovered target, for the summary, the overlay, and the restore after
      // a Review alignment preview; the highlight follows it after the dwell.
      const forgetFeatureHover = () => {
        handlerState.activeHoverSvgId = null;
        handlerState.activeHoverKey = '';
      };

      const forgetMatchHover = () => {
        handlerState.activeMatchHoverElement = null;
        handlerState.activeMatchHoverKey = '';
      };

      const clearHover = () => {
        clearHoverHighlight();
        forgetFeatureHover();
        forgetMatchHover();
      };

      /** @param {string} svgId */
      const featureHighlightKey = (svgId) => `feature|${getFeatureHoverKey(svgId)}`;
      /** @param {Element} matchElement */
      const matchHighlightKey = (matchElement) => `match|${getMatchHoverKey(matchElement)}`;

      const clearAlignmentCandidatePreview = ({ restore = true } = {}) => {
        if (!handlerState.alignmentCandidateSvgId && !handlerState.alignmentCandidateRestore) return;
        clearHoverHighlight();
        handlerState.alignmentCandidateSvgId = '';
        handlerState.alignmentOverlay?.setHover(null);
        const previous = handlerState.alignmentCandidateRestore;
        handlerState.alignmentCandidateRestore = null;
        if (!restore || !previous) return;
        // The hover under the pointer comes back at once: the pointer has rested there.
        const featureSvgId = previous.featureSvgId;
        if (featureSvgId && ensureFeatureLookup().has(featureSvgId)) {
          showHoverHighlight(featureHighlightKey(featureSvgId), () => setHoverHighlight(featureSvgId));
          handlerState.activeHoverSvgId = featureSvgId;
          handlerState.activeHoverKey = getFeatureHoverKey(featureSvgId);
          return;
        }
        const matchElement = previous.matchElement;
        if (matchElement?.isConnected) {
          showHoverHighlight(matchHighlightKey(matchElement), () => setMatchHover(matchElement));
          handlerState.activeMatchHoverElement = matchElement;
          handlerState.activeMatchHoverKey = getMatchHoverKey(matchElement);
        }
      };

      // The rendered features by identity, indexed once per mounted root: the
      // Review alignment asks for one feature per candidate (OV-935).
      const ensureRenderedFeatureIdentities = () => {
        if (!handlerState.renderedFeatureIdentities) {
          /** @type {RenderedIdentityEntry[]} */
          const entries = Array.from(ensureFeatureLookup(), ([svgId, feature]) => ({
            svgId, identity: renderedFeatureIdentity(feature)
          }));
          handlerState.renderedFeatureIdentities = indexFeatureIdentities(
            entries, (entry) => entry.identity, { includeRendered: false }
          );
        }
        return handlerState.renderedFeatureIdentities;
      };

      // The rendered feature an anchor names, when exactly one is.
      /** @param {AliasSource} anchor */
      const uniqueRenderedSvgId = (anchor) => {
        const requested = featureIdentity(anchor);
        if (!requested.usable) return '';
        const matches = ensureRenderedFeatureIdentities()(requested);
        return matches.length === 1 ? matches[0].svgId : '';
      };

      /** @param {AliasSource} anchor */
      const previewAlignmentCandidate = (anchor) => {
        const svgId = uniqueRenderedSvgId(anchor);
        if (!svgId) return false;
        const restore = handlerState.alignmentCandidateRestore || {
          featureSvgId: handlerState.activeHoverSvgId,
          matchElement: handlerState.activeMatchHoverElement
        };
        handlerState.alignmentCandidateSvgId = '';
        clearHover();
        hideHoverSummary();
        handlerState.alignmentCandidateRestore = restore;
        setFeatureHover(svgId);
        handlerState.alignmentCandidateSvgId = svgId;
        handlerState.alignmentOverlay?.setHover(handlerState.alignmentCandidatesBySvgId.get(svgId) ?? null);
        return true;
      };

      handlerState.previewAlignmentCandidate = previewAlignmentCandidate;
      handlerState.clearAlignmentCandidatePreview = clearAlignmentCandidatePreview;

      // A block the Result shows: hidden ones get no ring, no strip box and no place in a strip.
      /** @param {Element} block */
      const isDrawnBlock = (block) => block.isConnected && block.getClientRects().length === 1
        && window.getComputedStyle(block).visibility !== 'hidden';
      // A feature drawn in several parts (exons) is one span from its first block to its last:
      // the alignment centres that whole span.
      /** @param {AliasSource} anchor */
      const uniqueRenderedFeature = (anchor) => {
        const svgId = uniqueRenderedSvgId(anchor);
        if (!svgId) return null;
        const blocks = getFeatureFillElements(svg, svgId, ensureFeaturePaths());
        return blocks.length && blocks.every(isDrawnBlock) ? { svgId, blocks } : null;
      };

      const clearAlignmentOverlay = () => {
        const overlay = handlerState.alignmentOverlay;
        if (!overlay) return;
        if (overlay.frame) window.cancelAnimationFrame(overlay.frame);
        overlay.scrollRoot.removeEventListener('scroll', overlay.schedule);
        window.removeEventListener('resize', overlay.schedule);
        overlay.resizeObserver?.disconnect();
        overlay.transformObserver.disconnect();
        overlay.element.remove();
        handlerState.alignmentOverlay = null;
        handlerState.alignmentCandidatesBySvgId.clear();
        clearAlignmentCandidatePreview();
      };

      // The figure as the strips of the Review alignment panel need it: spans in SVG user
      // units from the Result's left edge, which the preview zoom does not change.
      /**
       * @param {AlignmentOverlayRequest['reference']} reference
       * @param {Element[]} referenceBlocks
       * @param {AlignmentOverlayRing[]} rings
       * @returns {AlignmentFigureGeometry | null}
       */
      const measureAlignmentGeometry = (reference, referenceBlocks, rings) => {
        const frame = svg.getBoundingClientRect();
        const unit = frame.width / (/** @type {SVGSVGElement} */ (svg).viewBox?.baseVal?.width || 0);
        if (!(unit > 0) || !Number.isFinite(unit)) return null;
        /** @param {Element} element @returns {FigureSpan} */
        const spanOf = (element) => {
          const rect = element.getBoundingClientRect();
          return { left: (rect.left - frame.left) / unit, right: (rect.right - frame.left) / unit };
        };
        /** @param {FigureSpan[]} spans @returns {FigureSpan} */
        const around = (spans) => ({
          left: Math.min(...spans.map(({ left }) => left)), right: Math.max(...spans.map(({ right }) => right))
        });
        /** @param {Element} element */
        const recordKeyOf = (element) => element.closest('[data-record-key]')?.getAttribute('data-record-key') || '';
        /** @type {Map<string, FigureSpan[]>} */
        const featuresByRecord = new Map();
        const featurePaths = ensureFeaturePaths();
        featurePaths.forEach((_elements, svgId) => getFeatureFillElements(svg, svgId, featurePaths).forEach((block) => {
          const recordKey = recordKeyOf(block);
          if (!recordKey || !isDrawnBlock(block)) return;
          if (!featuresByRecord.has(recordKey)) featuresByRecord.set(recordKey, []);
          featuresByRecord.get(recordKey)?.push(spanOf(block));
        }));
        const referenceSpan = around(referenceBlocks.map(spanOf));
        const drawnCandidates = [
          { recordKey: reference.recordKey, key: reference.key, span: referenceSpan },
          ...rings.map(({ record, key, blocks }) => ({ recordKey: record.recordKey, key, span: around(blocks.map(spanOf)) }))
        ];
        const records = Array.from(svg.querySelectorAll('g[data-record-key]')).flatMap((group) => {
          const recordKey = group.getAttribute('data-record-key') || '';
          const features = featuresByRecord.get(recordKey) || [];
          // The record's axis is its own line; a label's leader line is not.
          const axis = Array.from(group.querySelectorAll(':scope > line:not([data-label-feature-id])')).map(spanOf);
          const spans = [...axis, ...features];
          if (!spans.length) return [];
          return [{
            recordKey, features, extent: around(spans),
            candidates: drawnCandidates.filter((candidate) => candidate.recordKey === recordKey)
              .map(({ key, span }) => ({ key, span }))
          }];
        });
        return { referenceCenter: (referenceSpan.left + referenceSpan.right) / 2, records };
      };

      /** @param {AlignmentOverlayRequest} request */
      const showAlignmentOverlay = ({ reference, records, onSelect, onHover }) => {
        // A choice rebuilds the overlay; the ring that had the focus keeps it.
        const focused = document.activeElement?.closest?.('.gbdraw-alignment-ring');
        const focusedKeys = focused instanceof HTMLElement && handlerState.alignmentOverlay?.element.contains(focused)
          ? { ...focused.dataset } : null;
        clearAlignmentOverlay();
        const scrollRoot = svgContainer.value?.parentElement;
        const viewport = scrollRoot?.parentElement;
        if (!scrollRoot || !viewport || !svg.isConnected) return null;
        const referenceFeature = uniqueRenderedFeature(reference.anchor);

        const element = document.createElement('div');
        element.className = 'gbdraw-alignment-overlay';
        element.setAttribute('data-similarity-alignment-canvas', '');
        const drawing = document.createElementNS('http://www.w3.org/2000/svg', 'svg');
        drawing.setAttribute('aria-hidden', 'true');
        const referenceRing = document.createElement('span');
        referenceRing.className = 'gbdraw-alignment-ring is-reference is-chosen';
        const flag = document.createElement('span');
        flag.className = 'gbdraw-alignment-flag';
        const flagName = document.createElement('span');
        flagName.textContent = reference.name;
        flag.append(flagName, ' · reference');
        const tag = document.createElement('span');
        tag.className = 'gbdraw-alignment-tag';
        [referenceRing, flag, tag].forEach((part) => part.setAttribute('aria-hidden', 'true'));
        element.append(drawing, referenceRing);

        /** @type {AlignmentOverlay} */
        const overlay = {
          element,
          scrollRoot,
          frame: 0,
          resizeObserver: null,
          transformObserver: new MutationObserver(() => overlay.schedule()),
          hover: null,
          schedule: () => {
            if (!overlay.frame) overlay.frame = window.requestAnimationFrame(update);
          },
          setHover: (ring) => {
            if (overlay.hover === ring) return;
            overlay.hover = ring;
            onHover(ring?.record.recordKey || '', ring?.key || '');
            overlay.schedule();
          },
          onSelect
        };

        /** @type {AlignmentOverlayRing[]} */
        const candidates = records.flatMap((record) => record.candidates.flatMap((candidate) => {
          const rendered = uniqueRenderedFeature(candidate.anchor);
          return rendered ? [{ ...candidate, ...rendered, record, button: document.createElement('button') }] : [];
        }));
        const counts = new Map();
        candidates.forEach(({ svgId }) => counts.set(svgId, (counts.get(svgId) || 0) + 1));
        const rings = candidates.filter(({ svgId }) => counts.get(svgId) === 1);
        rings.forEach((ring) => {
          const { button, record } = ring;
          const chosen = !record.skipped && record.chosenKey === ring.key;
          handlerState.alignmentCandidatesBySvgId.set(ring.svgId, ring);
          button.type = 'button';
          button.className = 'gbdraw-alignment-ring';
          button.classList.toggle('is-chosen', chosen);
          button.classList.toggle('is-need', record.needsChoice);
          button.classList.toggle('is-skipped', record.skipped);
          button.dataset.alignmentRecordKey = record.recordKey;
          button.dataset.alignmentCandidateKey = ring.key;
          button.setAttribute('aria-label', ring.accessibleName);
          button.setAttribute('aria-pressed', String(chosen));
          button.addEventListener('click', (event) => {
            event.preventDefault();
            event.stopPropagation();
            onSelect(record.recordKey, ring.anchor);
          });
          button.addEventListener('mouseenter', () => overlay.setHover(ring));
          button.addEventListener('mouseleave', () => overlay.setHover(null));
          button.addEventListener('focus', () => overlay.setHover(ring));
          button.addEventListener('blur', () => overlay.setHover(null));
          element.appendChild(button);
        });
        element.append(flag, tag);
        viewport.appendChild(element);

        /**
         * The box around a feature's drawn blocks in the viewport, or null when its centre is
         * outside the viewport.
         * @param {Element[] | null | undefined} blocks
         * @param {DOMRect} frame
         */
        const boxIn = (blocks, frame) => {
          if (!blocks?.length) return null;
          const rects = blocks.map((block) => block.getBoundingClientRect());
          const box = {
            left: Math.min(...rects.map((rect) => rect.left)) - frame.left,
            right: Math.max(...rects.map((rect) => rect.right)) - frame.left,
            top: Math.min(...rects.map((rect) => rect.top)) - frame.top,
            bottom: Math.max(...rects.map((rect) => rect.bottom)) - frame.top
          };
          if (!(box.right > box.left) || !(box.bottom > box.top)) return null;
          const x = (box.left + box.right) / 2;
          const y = (box.top + box.bottom) / 2;
          return x < 0 || x > frame.width || y < 0 || y > frame.height ? null : { ...box, x, y };
        };
        /**
         * @param {HTMLElement} ring
         * @param {NonNullable<ReturnType<typeof boxIn>> | null} box
         */
        const placeRing = (ring, box) => {
          ring.hidden = !box;
          if (!box) return;
          const width = Math.max(box.right - box.left, 8) + 6;
          Object.assign(ring.style, { left: (box.x - width / 2) + 'px', top: (box.top - 3) + 'px',
            width: width + 'px', height: (box.bottom - box.top + 6) + 'px' });
        };
        /**
         * A label centred on `x`, kept inside the viewport.
         * @param {HTMLElement} label
         * @param {number} x
         * @param {number} top
         * @param {number} frameWidth
         */
        const placeLabel = (label, x, top, frameWidth) => {
          label.hidden = false;
          const width = label.offsetWidth;
          label.style.left = Math.max(4, Math.min(x - width / 2, frameWidth - width - 4)) + 'px';
          label.style.top = top + 'px';
        };
        // Down from the candidate, then along to the reference line.
        /**
         * @param {NonNullable<ReturnType<typeof boxIn>>} box
         * @param {number} lineX
         * @param {string} tone
         */
        const arrowTo = (box, lineX, tone) => {
          if (Math.abs(box.x - lineX) < 8) return '';
          const side = box.x > lineX ? 1 : -1;
          const y = box.bottom + 8;
          const shaft = `M${box.x.toFixed(1)},${(box.bottom + 4).toFixed(1)} V${y.toFixed(1)} H${(lineX + side * 7).toFixed(1)}`;
          const head = `M${(lineX + side * 1.5).toFixed(1)},${y.toFixed(1)} l${side * 7},-4 v8 z`;
          return `<g class="gbdraw-alignment-arrow ${tone}"><path class="casing" d="${shaft}"/><path class="casing head" d="${head}"/>`
            + `<path class="shaft" d="${shaft}"/><path class="head" d="${head}"/></g>`;
        };
        let settled = '';
        const update = () => {
          overlay.frame = 0;
          if (!element.isConnected || delegatedFeatureHandlers !== handlerState) return;
          // Hidden while the preview pans or zooms; the end of the gesture places it again (OV-936).
          if (isPreviewTransformInteractionActive()) return;
          const frame = viewport.getBoundingClientRect();
          const paper = svgContainer.value.getBoundingClientRect();
          const referenceBox = boxIn(referenceFeature?.blocks, frame);
          const lineX = referenceBox?.x ?? null;
          placeRing(referenceRing, referenceBox);
          let shapes = '';
          flag.hidden = lineX === null;
          if (lineX !== null) {
            const top = Math.max(paper.top - frame.top - 8, 22);
            const bottom = Math.min(paper.bottom - frame.top, frame.height);
            const line = `x1="${lineX.toFixed(1)}" x2="${lineX.toFixed(1)}" y1="${top.toFixed(1)}" y2="${bottom.toFixed(1)}"`;
            shapes += `<line class="gbdraw-alignment-line casing" ${line}/><line class="gbdraw-alignment-line" ${line}/>`;
            placeLabel(flag, lineX, top - 19, frame.width);
          }
          const hover = overlay.hover;
          tag.hidden = true;
          rings.forEach((ring) => {
            const box = boxIn(ring.blocks, frame);
            const hovered = hover === ring;
            placeRing(ring.button, box);
            ring.button.classList.toggle('is-hover', hovered);
            if (!box) return;
            if (hovered) {
              tag.textContent = ring.name;
              placeLabel(tag, box.x, box.top - 25, frame.width);
            }
            // The chosen candidate's arrow gives way to the hovered candidate of its record.
            const chosen = ring.button.classList.contains('is-chosen');
            if (lineX === null || !(hovered || (chosen && hover?.record !== ring.record))) return;
            shapes += arrowTo(box, lineX, [
              ring.record.needsChoice ? 'is-need' : ring.record.skipped ? 'is-skipped' : '', chosen ? '' : 'is-other'
            ].join(' '));
          });
          drawing.innerHTML = shapes;
          // A zoom or a fit moves the Result through a transition: follow it until it rests.
          const position = [paper.left, paper.top, paper.width, frame.width, frame.height].join(',');
          if (position !== settled || svgContainer.value.getAnimations().length > 0) {
            settled = position;
            overlay.schedule();
          }
        };
        scrollRoot.addEventListener('scroll', overlay.schedule, { passive: true });
        window.addEventListener('resize', overlay.schedule);
        if (window.ResizeObserver) {
          overlay.resizeObserver = new ResizeObserver(overlay.schedule);
          overlay.resizeObserver.observe(viewport);
          overlay.resizeObserver.observe(svgContainer.value);
        }
        // Zoom and pan are the surface's transform.
        overlay.transformObserver.observe(svgContainer.value, { attributes: true, attributeFilter: ['style'] });
        handlerState.alignmentOverlay = overlay;
        overlay.schedule();
        if (focusedKeys) {
          rings.find(({ button }) => button.dataset.alignmentRecordKey === focusedKeys.alignmentRecordKey
            && button.dataset.alignmentCandidateKey === focusedKeys.alignmentCandidateKey)?.button.focus({ preventScroll: true });
        }
        return referenceFeature ? measureAlignmentGeometry(reference, referenceFeature.blocks, rings) : null;
      };

      handlerState.showAlignmentOverlay = showAlignmentOverlay;
      handlerState.clearAlignmentOverlay = clearAlignmentOverlay;

      // The bands of the open match popup: a similarity group's ribbons, or one
      // match's fragments. Only the bands that change are written (OV-937).
      /** @param {SelectedPairwiseMatch | null} match */
      handlerState.markSelectedMatch = (match) => {
        const matchId = String(match?.id || '').trim();
        const orthogroupId = String(match?.matchKind === 'orthogroup' ? match?.orthogroupId : '').trim();
        /** @type {Element[]} */
        let selected = [];
        if (orthogroupId || matchId) {
          ensureComparisonIndexes();
          // ensureComparisonIndexes above has set all three comparison maps.
          selected = orthogroupId
            ? (/** @type {Map<string, Element[]>} */ (handlerState.comparisonElementsByOrthogroupId).get(orthogroupId) || [])
              .filter((element) => element.matches(PAIRWISE_MATCH_SELECTOR) && isOrthogroupRibbon(element))
            : /** @type {Map<string, Element[]>} */ (handlerState.comparisonElementsByMatchId).get(matchId) || [];
        }
        const next = new Set(selected);
        handlerState.selectedMatchElements.forEach((element) => {
          if (!next.has(element)) element.classList.toggle('gbdraw-match-selected', false);
        });
        next.forEach((element) => {
          if (!handlerState.selectedMatchElements.has(element)) element.classList.toggle('gbdraw-match-selected', true);
        });
        handlerState.selectedMatchElements = next;
      };

      const clearPendingMatch = () => {
        if (!handlerState.pendingMatchElement) return;
        matchFragments(handlerState.pendingMatchElement).forEach((element) => element.classList.remove('gbdraw-match-pending'));
        handlerState.pendingMatchElement = null;
      };

      /** @param {Element} matchElement */
      const setPendingMatch = (matchElement) => {
        if (handlerState.pendingMatchElement === matchElement) return;
        clearPendingMatch();
        handlerState.pendingMatchElement = matchElement;
        matchFragments(matchElement).forEach((element) => element.classList.add('gbdraw-match-pending'));
      };

      /** @param {PointerPosition | null | undefined} eventLike */
      const rememberTransformPointer = (eventLike) => {
        if (!Number.isFinite(eventLike?.clientX) || !Number.isFinite(eventLike?.clientY)) return;
        // `Number.isFinite(eventLike?.clientX)` holds only for a present eventLike.
        handlerState.transformPointer = {
          clientX: /** @type {PointerPosition} */ (eventLike).clientX,
          clientY: /** @type {PointerPosition} */ (eventLike).clientY
        };
        handlerState.reconcileHoverAfterTransform = true;
      };

      /** @param {Element} featureEl @param {HoverPoint | null | undefined} eventLike */
      const activateFeatureHover = (featureEl, eventLike) => {
        forgetMatchHover();
        const svgId = getFeatureIdentity(featureEl);
        if (!svgId) {
          clearHover();
          return false;
        }
        const hoverKey = getFeatureHoverKey(svgId);
        requestHoverHighlight(featureHighlightKey(svgId), () => setHoverHighlight(svgId));
        handlerState.activeHoverSvgId = svgId;
        handlerState.activeHoverKey = hoverKey;
        handlerState.alignmentOverlay?.setHover(handlerState.alignmentCandidatesBySvgId.get(svgId) ?? null);
        const feat = ensureFeatureLookup().get(svgId);
        scheduleHoverSummary(feat ? () => featureHoverContent(feat, featureEl) : null, eventLike);
        return true;
      };

      /** @param {Element} matchEl @param {HoverPoint | null | undefined} eventLike */
      const activateMatchHover = (matchEl, eventLike) => {
        handlerState.alignmentOverlay?.setHover(null);
        forgetFeatureHover();
        const matchKey = getMatchHoverKey(matchEl);
        requestHoverHighlight(matchHighlightKey(matchEl), () => setMatchHover(matchEl));
        handlerState.activeMatchHoverElement = matchEl;
        handlerState.activeMatchHoverKey = matchKey;
        scheduleHoverSummary(() => matchHoverContent(buildMatchHoverSummary(matchEl)), eventLike);
        return true;
      };

      handlerState.beginPreviewTransformInteraction = (eventLike, kind = '') => {
        if (handlerState.alignmentOverlay) handlerState.alignmentOverlay.element.hidden = true;
        handlerState.alignmentOverlay?.setHover(null);
        if (handlerState.hoverReconcileFrame !== null) {
          window.cancelAnimationFrame(handlerState.hoverReconcileFrame);
          handlerState.hoverReconcileFrame = null;
        }
        handlerState.reconcileHoverAfterTransform = Boolean(
          handlerState.activeHoverSvgId || handlerState.activeMatchHoverElement
        );
        rememberTransformPointer(eventLike);
        clearHover();
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

      /** @param {MouseEvent} e */
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

      /** @param {MouseEvent} e */
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

      /** @param {MouseEvent} e */
      const handleMouseOut = (e) => {
        if (isPreviewTransformInteractionActive()) {
          rememberTransformPointer(e);
          return;
        }
        // A move straight onto another target leaves the highlight to that
        // target's hover; a move off every target clears it at once.
        const relatedFeature = getFeatureTarget(e.relatedTarget, svg);
        const relatedMatch = relatedFeature ? null : getPairwiseMatchTarget(e.relatedTarget, svg);
        const leave = () => {
          if (!relatedFeature && !relatedMatch) clearHoverHighlight();
          hideHoverSummary();
        };
        const featureEl = getFeatureTarget(e.target, svg);
        if (featureEl) {
          const svgId = getFeatureIdentity(featureEl);
          if (!svgId || handlerState.activeHoverSvgId !== svgId) return;
          if (relatedFeature && getFeatureHoverKey(getFeatureIdentity(relatedFeature)) === handlerState.activeHoverKey) return;
          forgetFeatureHover();
          handlerState.alignmentOverlay?.setHover(null);
          leave();
          return;
        }
        const matchEl = getPairwiseMatchTarget(e.target, svg);
        if (!matchEl || handlerState.activeMatchHoverElement !== matchEl) return;
        if (relatedMatch && getMatchHoverKey(relatedMatch) === handlerState.activeMatchHoverKey) return;
        forgetMatchHover();
        leave();
      };

      /** @param {MouseEvent} e */
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
          } else if (isBackgroundPreviewClick(/** @type {Element} */ (e.target), svg)) { // a click on the Result targets an element
            featureSelection?.clearFeatureSelection?.({ clearStatus: true });
          }
          return;
        }
        const annotationEl = /** @type {MaybeElementTarget | null} */ (e.target)?.closest?.('[data-gbdraw-annotation-id]');
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
            handlerState.alignmentOverlay?.onSelect(candidate.record.recordKey, candidate.anchor);
            return;
          }
          e.stopPropagation();
          hideHoverSummary();
          featureSelection?.markPlainFeatureClick?.(svgId);
          if (previewPopupsHeld()) return;
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
        const previewItem = previewItemTarget(e.target, svg);
        if (previewItem && openPreviewItemPopup(previewItem, e)) {
          e.stopPropagation();
          e.preventDefault();
          return;
        }
        if (isBackgroundPreviewClick(/** @type {Element} */ (e.target), svg)) { // a click on the Result targets an element
          featureSelection?.clearFeatureSelection?.({ clearStatus: true });
          if (previewPopupsHeld()) return;
          if (clickedFeature.value) clickedFeature.value = null;
          if (clickedPairwiseMatch?.value) clickedPairwiseMatch.value = null;
        }
      };

      /** @param {KeyboardEvent} e */
      const handleKeyDown = (e) => {
        if (e.key !== 'Enter' && e.key !== ' ') return;
        const matchEl = getPairwiseMatchTarget(e.target, svg);
        if (matchEl) {
          e.stopPropagation();
          e.preventDefault();
          openPairwiseMatchPopup(matchEl, {}, ensureFeatureLookup());
          return;
        }
        const previewItem = previewItemTarget(e.target, svg);
        if (previewItem && openPreviewItemPopup(previewItem, {})) {
          e.stopPropagation();
          e.preventDefault();
        }
      };

      /** @param {PointerEvent} e */
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

      /** @param {PointerEvent} e */
      const handlePointerOut = (e) => {
        if (!handlerState.pendingMatchElement) return;
        const matchEl = getPairwiseMatchTarget(e.target, svg);
        if (matchEl !== handlerState.pendingMatchElement) return;
        if (getPairwiseMatchTarget(e.relatedTarget, svg) === matchEl) return;
        clearPendingMatch();
      };

      /** @param {PointerEvent} e */
      const handlePointerMove = (e) => {
        if (!featureSelection?.moveMarqueePointer?.(e)) return;
        e.preventDefault();
        e.stopPropagation();
        hideHoverSummary();
      };

      /** @param {PointerEvent} e */
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
        clearHover();
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
      // A start (`active`) carries the event that began it; only the end carries null.
      delegatedFeatureHandlers?.beginPreviewTransformInteraction?.(/** @type {Event} */ (event), kind);
      return;
    }
    delegatedFeatureHandlers?.endPreviewTransformInteraction?.({ reconcile });
  }) || null;

  const dispose = () => {
    if (typeof unsubscribePreviewTransformInteraction === 'function') {
      unsubscribePreviewTransformInteraction();
    }
    cleanupDelegatedFeatureHandlers();
    if (hoverSummaryState.view?.element.isConnected) hoverSummaryState.view.element.remove();
    hoverSummaryState.view = null;
  };

  /** @param {AliasSource} anchor */
  const previewAlignmentCandidate = (anchor) => (
    delegatedFeatureHandlers?.previewAlignmentCandidate?.(anchor) || false
  );

  const clearAlignmentCandidatePreview = () => {
    delegatedFeatureHandlers?.clearAlignmentCandidatePreview?.();
  };

  // Shows the overlay and returns the figure geometry the panel's strips are laid out from.
  /** @param {AlignmentOverlayRequest} request */
  const showAlignmentOverlay = (request) => delegatedFeatureHandlers?.showAlignmentOverlay?.(request) ?? null;

  const clearAlignmentOverlay = () => {
    delegatedFeatureHandlers?.clearAlignmentOverlay?.();
  };

  // Marks the bands of the open match popup on the mounted Result.
  /** @param {SelectedPairwiseMatch | null} match */
  const markSelectedPairwiseMatch = (match) => {
    delegatedFeatureHandlers?.markSelectedMatch?.(match);
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

  // Legend rows and record definitions are buttons that open their popup;
  // an export strips these attributes (stripTransientPreviewState).
  /** @param {{ root?: Element | null }} [options] */
  const preparePreviewItemAffordances = ({ root = null } = {}) => {
    const svg = root || delegatedFeatureHandlers?.svg || null;
    if (!svg) return false;
    svg.querySelectorAll(PREVIEW_ITEM_SELECTOR).forEach((element) => {
      if (element.closest(COMPARISON_LEGEND_SELECTOR)) return;
      element.setAttribute('data-gbdraw-preview-item', '');
      element.setAttribute('role', 'button');
      element.setAttribute('tabindex', '0');
      element.setAttribute('aria-haspopup', 'dialog');
    });
    return true;
  };

  return {
    attachSvgFeatureHandlers,
    getFeatureElements,
    getFeatureFillElements,
    openFeatureEditorForFeature,
    previewAlignmentCandidate,
    clearAlignmentCandidatePreview,
    showAlignmentOverlay,
    clearAlignmentOverlay,
    markSelectedPairwiseMatch,
    showFeaturePopupAtFeature,
    preparePairwiseInteractionAffordances,
    preparePreviewItemAffordances,
    dispose
  };
};
