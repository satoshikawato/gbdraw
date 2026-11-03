// @ts-check
/** @import { DrawingState } from '../../state.js' */
/** @import { DrawingId } from '../../services/drawing-view.js' */
import {
  buildFeatureSearchIndex,
  formatSearchMatchLine,
  getFeatureSearchFieldOptions,
  isRichFeatureSearchField,
  normalizeFeatureSearchField,
  runFeatureSearch
} from './search-core.js';
import {
  applyPreviewActiveSearchMatch,
  centerPreviewFeature,
  createPreviewFeatureSearchDomState,
  disposePreviewFeatureSearchDomState,
  getFeatureScreenCenter,
  getPreviewFeatureElementIndex,
  resolvePreviewSvg,
  schedulePreviewFeatureSearchClasses
} from './preview-svg.js';
import { recordStructuralMetric } from '../../services/runtime-test-hooks.js';
import { featureOverrideValue } from '../../services/feature-placement.js';
import { ownEntry } from '../../utils/own-entry.js';
import { t, tn } from '../../utils/i18n.js';
/** @import { FeatureSearchIndex } from './search-core.js' */
/** @import { FeatureLike } from '../../services/feature-utils.js' */
/** @import { FeatureElementIndex } from './preview-svg.js' */
/** @import { EditableLabelEntry } from '../feature-editor/label-actions.js' */

/**
 * @typedef {object} PreviewFeatureSearchOptions
 * @property {Record<string, any>} state the Web state; its shape belongs to `state.js`
 * @property {(source: any, callback: (...args: any[]) => void, options?: Record<string, any>) => any} watch Vue `watch`
 * @property {(callback?: () => void) => Promise<void>} nextTick Vue `nextTick`
 * @property {<T>(getter: () => T) => { value: T }} computed Vue `computed`
 * @property {(feature: Record<string, any>, eventLike?: { clientX: number, clientY: number } | null) => any} openFeatureEditorForFeature
 *   feature editor port: opens the popup for a feature at a screen point
 * @property {() => Record<string, any>[]} [resolveOrthogroups] orthogroups named for the search (default: `state.orthogroups`)
 * @property {(() => boolean) | null} [isActiveResultReady] preview runtime port
 * @property {() => void} focusSearchQuery focuses the search bar's query input
 * @property {() => void} focusSearchToggle focuses the toolbar button that shows and hides the search bar
 */

/**
 * @param {PreviewFeatureSearchOptions} options
 */
export const createPreviewFeatureSearch = ({
  state,
  watch,
  nextTick,
  computed,
  openFeatureEditorForFeature,
  resolveOrthogroups = () => state.orthogroups.value,
  isActiveResultReady = null,
  focusSearchQuery,
  focusSearchToggle
}) => {
  const {
    svgContainer,
    canvasContainerRef,
    canvasPan,
    selectedResultIndex,
    svgContent,
    featureList,
    featureListState,
    orthogroups,
    previewFeatureSearchInput,
    previewFeatureSearchQuery,
    previewFeatureSearchField,
    previewFeatureSearchQualifierKey,
    previewFeatureSearchUseRegex,
    previewFeatureSearchMatches,
    previewFeatureSearchMatchDetails,
    previewFeatureSearchActiveIndex,
    previewFeatureSearchError,
    previewFeatureSearchRenderedCount,
    previewFeatureSearchVisible,
    clickedFeature
  } = state;

  // The rich popup is an app preference (`state.richFeaturePopup`), not a drawing setting.
  const getPopupMode = () => (state.richFeaturePopup.value === false ? 'simple' : 'rich');
  let refreshRequestId = 0;
  let appliedSearchField = normalizeFeatureSearchField(previewFeatureSearchField.value, { popupMode: getPopupMode() });
  let appliedQualifierKey = String(previewFeatureSearchQualifierKey.value || '');
  let appliedUseRegex = Boolean(previewFeatureSearchUseRegex.value);
  /** @type {FeatureSearchIndex | null} */
  let searchIndex = null;
  /** @type {FeatureElementIndex | null} */
  let featureElementIndex = null;
  /** @type {SVGSVGElement | null} */
  let featureElementIndexSvg = null;
  const appliedSearchDomState = createPreviewFeatureSearchDomState();
  const getSvg = () => resolvePreviewSvg(svgContainer.value);
  const getActiveMatchId = () => (
    previewFeatureSearchActiveIndex.value >= 0
      ? String(previewFeatureSearchMatches.value?.[previewFeatureSearchActiveIndex.value] || '')
      : ''
  );
  const queryIsActive = () => Boolean(String(previewFeatureSearchQuery.value || '').trim()) && !previewFeatureSearchError.value;
  const invalidateSearchIndex = () => {
    searchIndex = null;
  };
  const invalidateFeatureElementIndex = () => {
    featureElementIndex = null;
    featureElementIndexSvg = null;
  };
  // Search features searches the Features list of the displayed Result, so a
  // hidden feature is found too (R-5). A drawn feature is found by its
  // rendered ID; one the Result does not draw by its source identity.
  /** @type {Map<string, FeatureLike>} */
  let featuresBySearchId = new Map();
  /** @param {DrawingState} drawing */
  const ensureSearchIndex = (drawing) => {
    if (!searchIndex) {
      /** @type {Map<string, EditableLabelEntry>} */
      const labels = new Map((state.editableLabels?.value || []).map((/** @type {EditableLabelEntry} */ entry) => [entry.featureId, entry]));
      featuresBySearchId = new Map();
      searchIndex = buildFeatureSearchIndex({
        features: featureList.value.rows.map((/** @type {FeatureLike} */ feature) => {
          const renderedId = featureListState(feature).renderedId;
          const searchId = renderedId || String(feature.id || '');
          featuresBySearchId.set(searchId, feature);
          const entry = renderedId ? labels.get(renderedId) : null;
          const sources = [entry?.sourceText, feature.label, feature.product, feature.gene, feature.locus_tag];
          return {
            ...feature,
            search_id: searchId,
            search_labels: [
              feature.search_labels,
              entry?.text,
              featureOverrideValue(drawing.featureOverrides, feature, 'labelText'),
              // Only a non-empty text names a bulk override; a feature without that source reads none (OV-577).
              ...sources.flatMap((source) => (
                typeof source === 'string' && source !== '' ? [ownEntry(drawing.labelTextBulkOverrides, source)] : []
              ))
            ]
          };
        }),
        popupMode: getPopupMode(),
        orthogroups: resolveOrthogroups()
      });
      recordStructuralMetric('featureSearchIndexBuildCount', 1, { phase: 'feature-search' });
    }
    return searchIndex;
  };
  const updateRenderedCount = () => {
    previewFeatureSearchRenderedCount.value =
      state.featureCatalog?.value?.items?.[selectedResultIndex.value]?.features?.length
      ?? state.extractedFeatures.value?.length ?? 0;
  };
  const ensureFeatureElementIndex = (svg = getSvg()) => {
    if (featureElementIndexSvg !== svg || !featureElementIndex) {
      featureElementIndexSvg = svg;
      featureElementIndex = getPreviewFeatureElementIndex(svg);
      recordStructuralMetric('featureDomFullScanCount', 1, { phase: 'feature-search' });
    }
    // The check above assigns the index unless it is already set.
    return /** @type {FeatureElementIndex} */ (featureElementIndex);
  };
  const previewFeatureSearchFieldOptions = computed(() => getFeatureSearchFieldOptions({ popupMode: getPopupMode() }));
  const previewFeatureSearchQualifierEnabled = computed(() => (
    getPopupMode() !== 'simple' && previewFeatureSearchField.value === 'qualifier-value'
  ));
  const previewFeatureSearchHasMatches = computed(() => previewFeatureSearchMatches.value.length > 0);
  const previewFeatureSearchCanOpenActive = computed(() => (
    previewFeatureSearchHasMatches.value && previewFeatureSearchActiveIndex.value >= 0
  ));
  const previewFeatureSearchCanSearch = computed(() => (
    Boolean(String(previewFeatureSearchInput.value || '').trim()) ||
    Boolean(String(previewFeatureSearchQuery.value || '').trim())
  ));
  const previewFeatureSearchStatusText = computed(() => {
    if (previewFeatureSearchError.value) return previewFeatureSearchError.value;
    // No count before a search; after one, the active match of all.
    if (!String(previewFeatureSearchQuery.value || '').trim()) return '';
    const total = previewFeatureSearchMatches.value.length;
    if (!total) return t('No matches');
    const current = previewFeatureSearchActiveIndex.value >= 0
      ? previewFeatureSearchActiveIndex.value + 1
      : 0;
    return tn(total, '{current} of {count} match', '{current} of {count} matches', { current });
  });
  const previewFeatureSearchActiveDetail = computed(() => {
    const activeId = getActiveMatchId();
    const details = activeId ? previewFeatureSearchMatchDetails.value?.[activeId] || [] : [];
    return formatSearchMatchLine(getActiveMatchFeature(), details);
  });

  const clearPreviewClasses = () => {
    disposePreviewFeatureSearchDomState(appliedSearchDomState);
  };

  /** @param {DrawingState} drawing */
  const refreshSearchNow = (drawing, { preserveActive = true, center = false } = {}) => {
    const popupMode = getPopupMode();
    const normalizedField = normalizeFeatureSearchField(appliedSearchField, { popupMode });
    appliedSearchField = normalizedField;

    const previousActiveId = preserveActive ? getActiveMatchId() : '';
    const svg = getSvg();
    if (!String(previewFeatureSearchQuery.value || '').trim()) {
      updateRenderedCount();
      previewFeatureSearchError.value = '';
      previewFeatureSearchMatches.value = [];
      previewFeatureSearchMatchDetails.value = {};
      previewFeatureSearchActiveIndex.value = -1;
      clearPreviewClasses();
      return;
    }
    const featureIndex = ensureFeatureElementIndex(svg);
    const index = ensureSearchIndex(drawing);
    const searchResult = runFeatureSearch({
      features: featureList.value.rows,
      renderedFeatureIds: new Set(index.featureOrder),
      query: previewFeatureSearchQuery.value,
      field: normalizedField,
      qualifierKey: appliedQualifierKey,
      useRegex: appliedUseRegex,
      popupMode,
      orthogroups: orthogroups.value,
      searchIndex: index,
      previousActiveId
    });

    previewFeatureSearchRenderedCount.value = searchResult.renderedFeatureCount;
    previewFeatureSearchError.value = searchResult.error;
    previewFeatureSearchMatches.value = searchResult.matches;
    previewFeatureSearchMatchDetails.value = searchResult.matchDetails;
    previewFeatureSearchActiveIndex.value = searchResult.activeIndex;
    schedulePreviewFeatureSearchClasses({
      svg,
      matches: searchResult.matches,
      activeId: getActiveMatchId(),
      queryActive: queryIsActive(),
      featureIndex,
      appliedState: appliedSearchDomState
    });
    if (center && getActiveMatchId()) {
      centerPreviewFeature({
        svg,
        featureId: getActiveMatchId(),
        featureIndex,
        canvasContainer: canvasContainerRef.value,
        canvasPan
      });
    }
  };

  const scheduleRefreshSearch = (options = {}) => {
    const drawing = state.activeDrawing();
    const requestId = ++refreshRequestId;
    nextTick(() => {
      if (requestId !== refreshRequestId) return;
      refreshSearchNow(drawing, options);
    });
  };

  /** @param {unknown} value */
  const setQuery = (value) => {
    previewFeatureSearchInput.value = String(value || '');
  };

  /** @param {unknown} field */
  const setField = (field) => {
    previewFeatureSearchField.value = normalizeFeatureSearchField(field, { popupMode: getPopupMode() });
  };

  const applySearch = () => {
    const popupMode = getPopupMode();
    const normalizedField = normalizeFeatureSearchField(previewFeatureSearchField.value, { popupMode });
    if (normalizedField !== previewFeatureSearchField.value) {
      previewFeatureSearchField.value = normalizedField;
    }
    appliedSearchField = normalizedField;
    appliedQualifierKey = String(previewFeatureSearchQualifierKey.value || '');
    appliedUseRegex = Boolean(previewFeatureSearchUseRegex.value);
    previewFeatureSearchQuery.value = String(previewFeatureSearchInput.value || '');
    scheduleRefreshSearch({ preserveActive: false });
  };

  const getActiveMatchFeature = () => {
    const activeId = getActiveMatchId();
    return activeId ? featuresBySearchId.get(activeId) || null : null;
  };

  const openActiveMatch = ({ center = true } = {}) => {
    if (!previewFeatureSearchMatches.value.length) return;
    if (previewFeatureSearchActiveIndex.value < 0) {
      previewFeatureSearchActiveIndex.value = 0;
    }
    const activeId = getActiveMatchId();
    const feature = getActiveMatchFeature();
    if (!feature) return;

    const svg = getSvg();
    const featureIndex = ensureFeatureElementIndex(svg);
    if (center) {
      centerPreviewFeature({
        svg,
        featureId: activeId,
        featureIndex,
        canvasContainer: canvasContainerRef.value,
        canvasPan
      });
    }
    applyPreviewActiveSearchMatch({
      featureIndex,
      appliedState: appliedSearchDomState,
      activeId
    });
    nextTick(() => {
      const centerPoint = getFeatureScreenCenter(getSvg(), activeId, featureIndex);
      openFeatureEditorForFeature(feature, centerPoint);
    });
  };

  /**
   * @param {unknown} index
   * @param {{ center?: boolean }} [options]
   */
  const goToMatch = (index, { center = true } = {}) => {
    const count = previewFeatureSearchMatches.value.length;
    if (!count) {
      const svg = getSvg();
      const featureIndex = ensureFeatureElementIndex(svg);
      previewFeatureSearchActiveIndex.value = -1;
      applyPreviewActiveSearchMatch({ featureIndex, appliedState: appliedSearchDomState, activeId: '' });
      return;
    }
    const previousActiveId = getActiveMatchId();
    previewFeatureSearchActiveIndex.value = ((Number(index) || 0) % count + count) % count;
    const activeId = getActiveMatchId();
    const svg = getSvg();
    const featureIndex = ensureFeatureElementIndex(svg);
    if (!appliedSearchDomState.queryActive || !appliedSearchDomState.matchedIds.has(activeId)) {
      schedulePreviewFeatureSearchClasses({
        svg,
        matches: previewFeatureSearchMatches.value,
        activeId,
        queryActive: queryIsActive(),
        featureIndex,
        appliedState: appliedSearchDomState
      });
    } else {
      appliedSearchDomState.activeId = previousActiveId;
      applyPreviewActiveSearchMatch({ featureIndex, appliedState: appliedSearchDomState, activeId });
    }
    if (center && activeId) {
      centerPreviewFeature({
        svg,
        featureId: activeId,
        featureIndex,
        canvasContainer: canvasContainerRef.value,
        canvasPan
      });
    }
    if (clickedFeature?.value) {
      openActiveMatch({ center: false });
    }
  };

  const goToNext = () => goToMatch(previewFeatureSearchActiveIndex.value + 1);
  const goToPrevious = () => goToMatch(previewFeatureSearchActiveIndex.value - 1);

  /**
   * @typedef {object} SearchDraft
   * @property {string} input
   * @property {string} query
   * @property {string} field
   * @property {string} qualifierKey
   * @property {boolean} useRegex
   * @property {string} appliedField
   * @property {string} appliedQualifierKey
   * @property {boolean} appliedUseRegex
   */
  /** @type {Readonly<SearchDraft>} */
  const EMPTY_SEARCH = Object.freeze({
    input: '', query: '', field: 'all', qualifierKey: '', useRegex: false,
    appliedField: 'all', appliedQualifierKey: '', appliedUseRegex: false
  });
  /** @returns {SearchDraft} */
  const captureSearch = () => ({
    input: String(previewFeatureSearchInput.value || ''),
    query: String(previewFeatureSearchQuery.value || ''),
    field: String(previewFeatureSearchField.value || 'all'),
    qualifierKey: String(previewFeatureSearchQualifierKey.value || ''),
    useRegex: Boolean(previewFeatureSearchUseRegex.value),
    appliedField: appliedSearchField,
    appliedQualifierKey,
    appliedUseRegex
  });
  /** @param {Readonly<SearchDraft>} draft */
  const installSearch = (draft) => {
    previewFeatureSearchInput.value = draft.input;
    previewFeatureSearchQuery.value = draft.query;
    previewFeatureSearchField.value = draft.field;
    previewFeatureSearchQualifierKey.value = draft.qualifierKey;
    previewFeatureSearchUseRegex.value = draft.useRegex;
    appliedSearchField = draft.appliedField;
    appliedQualifierKey = draft.appliedQualifierKey;
    appliedUseRegex = draft.appliedUseRegex;
    previewFeatureSearchMatches.value = [];
    previewFeatureSearchMatchDetails.value = {};
    previewFeatureSearchActiveIndex.value = -1;
    previewFeatureSearchError.value = '';
    clearPreviewClasses();
  };

  const clearSearch = () => {
    installSearch(EMPTY_SEARCH);
    scheduleRefreshSearch({ preserveActive: false });
  };

  // A search belongs to the drawing it was typed in: the drawing transition
  // keeps the departing drawing's search and installs the arriving one's,
  // which the arriving Result's mount applies (`handleMountedResultReady`).
  /** @type {Map<DrawingId, SearchDraft>} */
  const searchByDrawing = new Map();
  /** @param {DrawingId} departing @param {DrawingId} arriving */
  const switchDrawingSearch = (departing, arriving) => {
    searchByDrawing.set(departing, captureSearch());
    installSearch(searchByDrawing.get(arriving) || EMPTY_SEARCH);
    searchByDrawing.delete(arriving);
  };
  // Hiding the search bar changes its visibility only: the search, its
  // matches, and its highlights stay. While an active search dims the diagram
  // (`queryIsActive`, the same test that sets the SVG root class), the toolbar
  // button names its match count, 0 included; otherwise the count is null.
  // Showing the bar from the toolbar focuses the query; hiding it returns
  // focus to the toolbar button.
  const hideSearchBar = () => {
    previewFeatureSearchVisible.value = false;
    nextTick(focusSearchToggle);
  };
  const toggleSearchBar = () => {
    previewFeatureSearchVisible.value = !previewFeatureSearchVisible.value;
    if (previewFeatureSearchVisible.value) nextTick(focusSearchQuery);
  };
  const previewFeatureSearchHiddenMatchCount = computed(() => (
    previewFeatureSearchVisible.value || !queryIsActive() ? null : previewFeatureSearchMatches.value.length
  ));

  // A loaded Session replaces the Results of both modes that a search named,
  // and shows the search bar again.
  const resetModeSearches = () => {
    searchByDrawing.clear();
    clearSearch();
    previewFeatureSearchVisible.value = true;
  };

  watch(
    () => getPopupMode(),
    () => {
      if (getPopupMode() === 'simple' && isRichFeatureSearchField(previewFeatureSearchField.value)) {
        previewFeatureSearchField.value = 'all';
      }
      if (getPopupMode() === 'simple' && isRichFeatureSearchField(appliedSearchField)) {
        appliedSearchField = 'all';
      }
      invalidateSearchIndex();
      if (queryIsActive() && isActiveResultReady?.()) scheduleRefreshSearch();
    }
  );
  watch([featureList, orthogroups], () => {
    updateRenderedCount();
    invalidateSearchIndex();
    if (queryIsActive() && isActiveResultReady?.()) scheduleRefreshSearch();
  });
  watch([selectedResultIndex, svgContent], () => {
    invalidateFeatureElementIndex();
  });
  watch([
    () => (state.editableLabels?.value || []).map((/** @type {EditableLabelEntry} */ entry) => [entry.featureId, entry.sourceText, entry.text]),
    () => Object.values(state.activeDrawing().featureOverrides || {}).map((row) => [row.recordKey, row.biologicalFeatureId, row.labelText]),
    () => state.activeDrawing().labelTextBulkOverrides,
    () => state.activeDrawing().orthogroupNameOverrides,
    () => state.activeDrawing().orthogroupDescriptionOverrides
  ], () => {
    invalidateSearchIndex();
    if (queryIsActive() && isActiveResultReady?.()) scheduleRefreshSearch();
  }, { deep: true });

  updateRenderedCount();

  const handleMountedResultReady = () => {
    updateRenderedCount();
    invalidateFeatureElementIndex();
    clearPreviewClasses();
    if (queryIsActive()) scheduleRefreshSearch({ preserveActive: false });
  };

  const dispose = () => {
    refreshRequestId += 1;
    disposePreviewFeatureSearchDomState(appliedSearchDomState);
  };

  return {
    previewFeatureSearchFieldOptions,
    previewFeatureSearchQualifierEnabled,
    previewFeatureSearchHasMatches,
    previewFeatureSearchCanOpenActive,
    previewFeatureSearchCanSearch,
    previewFeatureSearchStatusText,
    previewFeatureSearchActiveDetail,
    previewFeatureSearchHiddenMatchCount,
    hideSearchBar,
    toggleSearchBar,
    setQuery,
    setField,
    applySearch,
    goToNext,
    goToPrevious,
    clearSearch,
    switchDrawingSearch,
    resetModeSearches,
    handleMountedResultReady,
    openActiveMatch,
    refreshSearch: scheduleRefreshSearch,
    dispose
  };
};
