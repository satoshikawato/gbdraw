import {
  formatSearchMatchDetail,
  getFeatureSearchFieldOptions,
  isRichFeatureSearchField,
  normalizeFeatureSearchField,
  runFeatureSearch
} from './search-core.js';
import {
  applyPreviewFeatureSearchClasses,
  centerPreviewFeature,
  getFeatureScreenCenter,
  getPreviewFeatureElementIndex,
  resolvePreviewSvg,
  stripPreviewFeatureSearchClasses
} from './preview-svg.js';

export const createPreviewFeatureSearch = ({
  state,
  watch,
  nextTick,
  computed,
  openFeatureEditorForFeature
}) => {
  const {
    svgContainer,
    canvasContainerRef,
    canvasPan,
    selectedResultIndex,
    svgContent,
    extractedFeatures,
    featuresBySvgId,
    orthogroups,
    adv,
    previewFeatureSearchInput,
    previewFeatureSearchQuery,
    previewFeatureSearchField,
    previewFeatureSearchQualifierKey,
    previewFeatureSearchUseRegex,
    previewFeatureSearchMatches,
    previewFeatureSearchMatchDetails,
    previewFeatureSearchActiveIndex,
    previewFeatureSearchError,
    previewFeatureSearchRenderedCount
  } = state;

  const getPopupMode = () => (adv?.rich_feature_popup === false ? 'simple' : 'rich');
  let refreshRequestId = 0;
  let appliedSearchField = normalizeFeatureSearchField(previewFeatureSearchField.value, { popupMode: getPopupMode() });
  let appliedQualifierKey = String(previewFeatureSearchQualifierKey.value || '');
  let appliedUseRegex = Boolean(previewFeatureSearchUseRegex.value);
  const getSvg = () => resolvePreviewSvg(svgContainer.value);
  const getActiveMatchId = () => (
    previewFeatureSearchActiveIndex.value >= 0
      ? String(previewFeatureSearchMatches.value?.[previewFeatureSearchActiveIndex.value] || '')
      : ''
  );
  const queryIsActive = () => Boolean(String(previewFeatureSearchQuery.value || '').trim()) && !previewFeatureSearchError.value;

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
    if (!String(previewFeatureSearchQuery.value || '').trim()) {
      return `0 / ${previewFeatureSearchRenderedCount.value} features`;
    }
    const current = previewFeatureSearchActiveIndex.value >= 0
      ? previewFeatureSearchActiveIndex.value + 1
      : 0;
    return `${current} / ${previewFeatureSearchMatches.value.length} features`;
  });
  const previewFeatureSearchActiveDetail = computed(() => {
    const activeId = getActiveMatchId();
    const details = activeId ? previewFeatureSearchMatchDetails.value?.[activeId] || [] : [];
    const detailText = details.length ? formatSearchMatchDetail(details[0]) : '';
    return detailText ? `Matched ${detailText}` : '';
  });

  const syncPreviewClasses = ({ center = false } = {}) => {
    const svg = getSvg();
    if (!svg) return;
    const activeId = getActiveMatchId();
    const featureIndex = getPreviewFeatureElementIndex(svg);
    applyPreviewFeatureSearchClasses({
      svg,
      matches: previewFeatureSearchMatches.value,
      activeId,
      queryActive: queryIsActive(),
      featureIndex
    });
    if (center && activeId) {
      centerPreviewFeature({
        svg,
        featureId: activeId,
        featureIndex,
        canvasContainer: canvasContainerRef.value,
        canvasPan
      });
    }
  };

  const clearPreviewClasses = () => {
    stripPreviewFeatureSearchClasses(getSvg());
  };

  const refreshSearchNow = ({ preserveActive = true, center = false } = {}) => {
    const popupMode = getPopupMode();
    const normalizedField = normalizeFeatureSearchField(appliedSearchField, { popupMode });
    appliedSearchField = normalizedField;

    const previousActiveId = preserveActive ? getActiveMatchId() : '';
    const svg = getSvg();
    const featureIndex = getPreviewFeatureElementIndex(svg);
    const renderedFeatureIds = new Set(featureIndex.keys());
    const searchResult = runFeatureSearch({
      features: extractedFeatures.value,
      renderedFeatureIds,
      query: previewFeatureSearchQuery.value,
      field: normalizedField,
      qualifierKey: appliedQualifierKey,
      useRegex: appliedUseRegex,
      popupMode,
      orthogroups: orthogroups.value,
      previousActiveId
    });

    previewFeatureSearchRenderedCount.value = searchResult.renderedFeatureCount;
    previewFeatureSearchError.value = searchResult.error;
    previewFeatureSearchMatches.value = searchResult.matches;
    previewFeatureSearchMatchDetails.value = searchResult.matchDetails;
    previewFeatureSearchActiveIndex.value = searchResult.activeIndex;
    applyPreviewFeatureSearchClasses({
      svg,
      matches: searchResult.matches,
      activeId: getActiveMatchId(),
      queryActive: queryIsActive(),
      featureIndex
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
    const requestId = ++refreshRequestId;
    nextTick(() => {
      if (requestId !== refreshRequestId) return;
      refreshSearchNow(options);
    });
  };

  const setQuery = (value) => {
    previewFeatureSearchInput.value = String(value || '');
  };

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

  const goToMatch = (index, { center = true } = {}) => {
    const count = previewFeatureSearchMatches.value.length;
    if (!count) {
      previewFeatureSearchActiveIndex.value = -1;
      syncPreviewClasses();
      return;
    }
    previewFeatureSearchActiveIndex.value = ((Number(index) || 0) % count + count) % count;
    syncPreviewClasses({ center });
  };

  const goToNext = () => goToMatch(previewFeatureSearchActiveIndex.value + 1);
  const goToPrevious = () => goToMatch(previewFeatureSearchActiveIndex.value - 1);

  const clearSearch = () => {
    previewFeatureSearchInput.value = '';
    previewFeatureSearchQuery.value = '';
    previewFeatureSearchField.value = 'all';
    previewFeatureSearchQualifierKey.value = '';
    previewFeatureSearchUseRegex.value = false;
    appliedSearchField = 'all';
    appliedQualifierKey = '';
    appliedUseRegex = false;
    previewFeatureSearchMatches.value = [];
    previewFeatureSearchMatchDetails.value = {};
    previewFeatureSearchActiveIndex.value = -1;
    previewFeatureSearchError.value = '';
    clearPreviewClasses();
    scheduleRefreshSearch({ preserveActive: false });
  };

  const openActiveMatch = () => {
    if (!previewFeatureSearchMatches.value.length) return;
    if (previewFeatureSearchActiveIndex.value < 0) {
      previewFeatureSearchActiveIndex.value = 0;
    }
    const activeId = getActiveMatchId();
    const feature = featuresBySvgId.value?.get?.(activeId) ||
      (Array.isArray(extractedFeatures.value)
        ? extractedFeatures.value.find((candidate) => String(candidate?.svg_id || '').trim() === activeId)
        : null);
    if (!feature) return;

    const svg = getSvg();
    const featureIndex = getPreviewFeatureElementIndex(svg);
    centerPreviewFeature({
      svg,
      featureId: activeId,
      featureIndex,
      canvasContainer: canvasContainerRef.value,
      canvasPan
    });
    syncPreviewClasses();
    nextTick(() => {
      const center = getFeatureScreenCenter(getSvg(), activeId, featureIndex);
      openFeatureEditorForFeature(feature, center);
    });
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
      scheduleRefreshSearch();
    }
  );
  watch(
    [
      selectedResultIndex,
      extractedFeatures,
      orthogroups,
      svgContent
    ],
    () => scheduleRefreshSearch()
  );

  scheduleRefreshSearch({ preserveActive: false });

  const dispose = () => {
    refreshRequestId += 1;
  };

  return {
    previewFeatureSearchFieldOptions,
    previewFeatureSearchQualifierEnabled,
    previewFeatureSearchHasMatches,
    previewFeatureSearchCanOpenActive,
    previewFeatureSearchCanSearch,
    previewFeatureSearchStatusText,
    previewFeatureSearchActiveDetail,
    setQuery,
    setField,
    applySearch,
    goToNext,
    goToPrevious,
    clearSearch,
    openActiveMatch,
    refreshSearch: scheduleRefreshSearch,
    dispose
  };
};
