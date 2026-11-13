// @ts-check
/** @import { DrawingState } from '../../state.js' */
/** @import { DrawingId } from '../../services/drawing-view.js' */
import {
  buildFeatureSearchIndex,
  formatSearchMatchLine,
  getFeatureSearchFieldOptions,
  isRichFeatureSearchField,
  isSequenceFeatureSearchField,
  normalizeFeatureSearchField,
  runFeatureSearch
} from './search-core.js';
import { featureSequenceReads } from '../../services/catalog-sequences.js';
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
import { normalizeCaughtError } from '../../utils/error-normalization.js';
import { reactive } from '../../utils/reactivity.js';
/** @import { FeatureSearchIndex, FeatureSearchSequences } from './search-core.js' */
/** @import { FeatureSequenceRead } from '../../services/catalog-sequences.js' */
/** @import { UserFacingError } from '../../utils/error-normalization.js' */
/** @import { FeatureLike } from '../../services/feature-utils.js' */
/** @import { FeatureElementIndex } from './preview-svg.js' */
/** @import { EditableLabelEntry } from '../../services/label-override-table.js' */

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
 * @property {(features: readonly FeatureSequenceRead[]) => Promise<FeatureSearchSequences | null>} readFeatureSequences
 *   reads the features' sequences from Python for the shown Result (the catalog-sequence door); null when the shown
 *   Result changed before they arrived; throws on failure
 */

/**
 * The sequences the Nucleotide and Amino acid fields search (D5): Python reads them for the shown Result's features
 * once, and they are kept only with that Result's features (`drop` when the Features list changes). While they are
 * read the panel says so; a failure keeps its error for Retry. A read is current while its number is the last one,
 * so a reply to a read that ended, or a stale one (null: the shown Result changed), writes nothing.
 * @param {PreviewFeatureSearchOptions['readFeatureSequences']} readFeatureSequences
 * @param {PreviewFeatureSearchOptions['computed']} computed
 */
const createSearchSequences = (readFeatureSequences, computed) => {
  /** @type {FeatureSearchSequences | null} */
  let sequences = null;
  let readId = 0;
  const status = reactive({ reading: false, error: /** @type {UserFacingError | null} */ (null) });
  /** Ends the read: a reply still on its way writes nothing. */
  const end = () => {
    readId += 1;
    if (status.reading) status.reading = false;
    if (status.error) status.error = null;
  };
  /**
   * Reads the features' sequences unless a read runs, then calls `done` with them kept.
   * @param {readonly FeatureLike[]} features @param {() => void} done
   */
  const read = async (features, done) => {
    if (status.reading) return;
    const id = ++readId;
    Object.assign(status, { reading: true, error: null });
    try {
      const reply = await readFeatureSequences(featureSequenceReads(features));
      if (id !== readId) return;
      status.reading = false;
      if (!reply) return;
      sequences = reply;
      done();
    } catch (error) {
      if (id !== readId) return;
      Object.assign(status, {
        reading: false, error: normalizeCaughtError(error, { operation: 'readFeatureSequences', stage: 'helper' })
      });
    }
  };
  /**
   * Whether a search of `field` can run now. Another field can, and ends a read; a sequence field can once the
   * sequences are read, and else starts the read, which schedules the search again with `options` when they are.
   * @param {string} field @param {readonly FeatureLike[]} features
   * @param {(options: { preserveActive: boolean, center: boolean }) => void} schedule
   * @param {{ preserveActive: boolean, center: boolean }} options
   */
  const ready = (field, features, schedule, options) => {
    if (!isSequenceFeatureSearchField(field)) end();
    else if (!sequences) void read(features, () => schedule(options));
    return !isSequenceFeatureSearchField(field) || sequences !== null;
  };
  return {
    status,
    ready,
    end,
    /** The sequences read for the shown Result, or null before a sequence field searched it. */
    current: () => sequences,
    drop: () => {
      sequences = null;
      end();
    },
    /** The status line of a read that runs or failed, or null. */
    statusText: () => (status.reading ? t('Reading the sequence…') : status.error ? '' : null),
    /**
     * A failed read while its search is shown, with its Retry: that search runs, and so reads, again.
     * @param {() => boolean} active @param {(options: { preserveActive: boolean }) => void} schedule
     */
    failure: (active, schedule) => computed(() => (active() && status.error ? { error: status.error, retry: () => schedule({ preserveActive: false }) } : null))
  };
};

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
  focusSearchToggle,
  readFeatureSequences
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
  const searchSequences = createSearchSequences(readFeatureSequences, computed);
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
    if (!searchIndex || searchIndex.sequencesRead !== (searchSequences.current() !== null)) {
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
        orthogroups: resolveOrthogroups(),
        sequences: searchSequences.current()
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
    if (searchSequences.statusText() !== null) return searchSequences.statusText();
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
    const normalizedField = (appliedSearchField = normalizeFeatureSearchField(appliedSearchField, { popupMode }));

    const previousActiveId = preserveActive ? getActiveMatchId() : '';
    const svg = getSvg();
    // A search of a sequence field waits for the sequences Python reads, then runs again (D5).
    if (!String(previewFeatureSearchQuery.value || '').trim() || !searchSequences.ready(normalizedField,
      featureList.value.rows, scheduleRefreshSearch, { preserveActive, center })) {
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

  // What a change of the index's inputs does: the index is built again, and an active search runs again.
  const reindexSearch = () => {
    invalidateSearchIndex();
    if (queryIsActive() && isActiveResultReady?.()) scheduleRefreshSearch();
  };
  watch(
    getPopupMode,
    () => {
      if (getPopupMode() === 'simple' && isRichFeatureSearchField(previewFeatureSearchField.value)) {
        previewFeatureSearchField.value = 'all';
      }
      if (getPopupMode() === 'simple' && isRichFeatureSearchField(appliedSearchField)) {
        appliedSearchField = 'all';
      }
      reindexSearch();
    }
  );
  watch([featureList, orthogroups], () => {
    updateRenderedCount();
    searchSequences.drop();
    reindexSearch();
  });
  watch([selectedResultIndex, svgContent], invalidateFeatureElementIndex);
  watch([
    () => (state.editableLabels?.value || []).map((/** @type {EditableLabelEntry} */ entry) => [entry.featureId, entry.sourceText, entry.text]),
    () => Object.values(state.activeDrawing().featureOverrides || {}).map((row) => [row.recordKey, row.biologicalFeatureId, row.labelText]),
    () => state.activeDrawing().labelTextBulkOverrides,
    () => state.activeDrawing().orthogroupNameOverrides,
    () => state.activeDrawing().orthogroupDescriptionOverrides
  ], reindexSearch, { deep: true });

  updateRenderedCount();

  const handleMountedResultReady = () => {
    updateRenderedCount();
    invalidateFeatureElementIndex();
    clearPreviewClasses();
    if (queryIsActive()) scheduleRefreshSearch({ preserveActive: false });
  };

  const dispose = () => {
    refreshRequestId += 1;
    searchSequences.end();
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
    previewFeatureSearchSequenceFailure: searchSequences.failure(queryIsActive, scheduleRefreshSearch),
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
