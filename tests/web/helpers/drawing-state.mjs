// The drawing context of the Web state (gbdraw/web/js/state.js): every `state`
// key in exactly one class. A drawing holds one diagram mode's settings, its
// editor edits, and the values derived from them, under their `state` names and
// kinds; services read them from the drawing they are given
// (tests/web/drawing-context.test.mjs). Every other key is a generated
// artifact, a project input or cache, an app setting or catalog, or transient.
// Both modes share one drawing until the per-mode settings change.

export const DRAWING_DRAFT_KEYS = Object.freeze([
  'form', 'adv', 'losat', 'losatProgram', 'circularConservation', 'linearComparisonPlan',
  'linearRecordLayoutEnabled', 'linearRecordGap', 'linearRecordRows', 'recordDisplayDrafts',
  'unmanagedConfigOverrides', 'importedComparisonIntent', 'layoutPreferences', 'linearTypographyLinked',
  'modeProfileStateManager'
]);

export const DRAWING_EDITOR_KEYS = Object.freeze([
  'featureOverrides', 'featurePlacementOverrides', 'featureColorOverrides', 'featureStrokeOverrides',
  'featureVisibilityManualRules', 'labelTextBulkOverrides', 'canonicalLabelOverrideRows', 'legendEntries',
  'deletedLegendEntries', 'legendColorOverrides', 'legendStrokeOverrides', 'addedLegendCaptions',
  'fileLegendCaptions', 'manualSpecificRules', 'manualPriorityRules', 'filterMode', 'manualBlacklist',
  'manualWhitelist', 'selectedPalette', 'currentColors', 'pendingPaletteName', 'pendingPaletteColors',
  'annotationSets', 'orthogroupNameOverrides', 'orthogroupDescriptionOverrides', 'orthogroupDormantOverrides',
  'canvasPadding'
]);

// Computed from the members above (and project inputs); each drawing derives its own.
export const DRAWING_DERIVED_KEYS = Object.freeze([
  'activeLayoutPreferences', 'linearComparisonResolution', 'hasLinearComparisonIntent',
  'hasActiveLinearLosatIntent', 'hasActiveLinearUploadIntent', 'featureVisibilityRules', 'hasPendingPaletteDraft'
]);

export const DRAWING_KEYS = Object.freeze([...DRAWING_DRAFT_KEYS, ...DRAWING_EDITOR_KEYS, ...DRAWING_DERIVED_KEYS]);

// What a Generate produces and shows: the Results and their catalog, run facts,
// generated Legend inventory and colors, and the user's offsets on the Result.
export const ARTIFACT_KEYS = Object.freeze([
  'results', 'selectedResultIndex', 'lastRunInfo', 'trackSlotResolvedGeometry', 'annotationWarnings',
  'featureIdentityNotices', 'featureEditRemovalCount', 'comparisonWarnings', 'pairwiseMatchFactors',
  'orthogroups', 'collinearGroups', 'featureOrthogroupIndex', 'similarityAlignmentPlan',
  'similarityAlignmentResetReceipt', 'linearRecordTranslations', 'legacySimilarityAlignment',
  'appliedPaletteName', 'appliedPaletteColors', 'extractedFeatures', 'biologicalFeatures', 'featureCatalog',
  'featureRecordIds', 'editableLabels', 'originalLegendOrder', 'originalLegendColors', 'originalSvgStroke',
  'legendCurrentOffset', 'diagramOffset', 'lengthBarUserOffset', 'plotTitleUserOffset',
  'generatedLegendPosition', 'generatedMode', 'generatedMultiRecordCanvas', 'generatedCircularPlotTitlePosition'
]);

// Inputs and caches that every drawing shares.
export const PROJECT_KEYS = Object.freeze([
  'sessionTitle', 'files', 'linearSeqs', 'cInputType', 'lInputType', 'circularRecordList',
  'circularRecordDiscovery', 'losatCache', 'losatDerivedCache', 'losatCacheInfo', 'proteinIdentityManifest',
  'legacyProteinRawCandidates', 'legacyProteinDerivedEvidence'
]);

// App settings and catalogs outside any drawing, and the shown mode.
export const UI_KEYS = Object.freeze([
  'mode', 'downloadDpi', 'autoLabelReflowEnabled', 'paletteInstantPreviewEnabled', 'paletteDefinitions',
  'paletteNames', 'specificRulePresets', 'featureKeys', 'defaultColorKeys'
]);

// Never saved: progress and operation flags, input drafts, selection, search,
// dialogs, popups, drag and view state, and views derived from the others.
export const TRANSIENT_KEYS = Object.freeze([
  'processing', 'processingStatus', 'sessionSavePending', 'sessionImportPending', 'sessionOperationAvailability',
  'generationCancelRequested', 'errorLog', 'semanticFileWatchersSuppressed', 'sessionResourceDiscoveryDeferred',
  'sessionImportRollbackInProgress', 'failedGeneratePreservedResult', 'generationFailureRecovery',
  'resultPanelTab', 'matchSequenceRegistry', 'svgContent', 'svgResultIdentity', 'zoom', 'layoutRepositionMode',
  'isPanning', 'panStart', 'canvasPan', 'canvasContainerRef', 'suppressCircularMultiRecordDefaults',
  'selectedAnnotation', 'losatThreadingStatus', 'selectedOrthogroupAlignmentFeature', 'selectedOrthogroupId',
  'orthogroupSearch', 'orthogroupSortMode', 'showRightDrawer', 'rightDrawerTab', 'linearReorderNotice',
  'newSpecRule', 'specificRuleQualifierSuggestions', 'selectedSpecificPreset', 'specificRulePresetLoading',
  'featuresBySvgId', 'selectedFeatureIds', 'selectedFeatureAnchorId', 'featureSelectionStatus',
  'featureSelectionSuppressNextClick', 'featureSelectionDrag', 'selectedFeatureCount', 'selectedFeatures',
  'hasFeatureSelection', 'featureEditorStatus', 'featureEditorStatusText', 'featureExtractionPending',
  'featureExtractionError', 'selectedFeatureRecordIdx', 'resultGenerationKey', 'featurePanelTab',
  'featureSearchInput', 'featureSearch', 'previewFeatureSearchInput', 'previewFeatureSearchQuery',
  'previewFeatureSearchField', 'previewFeatureSearchQualifierKey', 'previewFeatureSearchUseRegex',
  'previewFeatureSearchMatches', 'previewFeatureSearchMatchDetails', 'previewFeatureSearchActiveIndex',
  'previewFeatureSearchError', 'previewFeatureSearchRenderedCount', 'featureListScrollTop',
  'featureListViewportHeight', 'isFeatureDrawerMounted', 'visibleFeatureRows', 'featureRecordPickerVisible',
  'featureListTopSpacerPx', 'featureListBottomSpacerPx', 'labelSearch', 'labelReflowProcessing',
  'labelReflowRequestSeq', 'labelReflowForceRequestSeq', 'labelReflowLastError', 'svgContainer',
  'clickedFeature', 'clickedFeaturePos', 'clickedPairwiseMatch', 'clickedPairwiseMatchPos',
  'pairwiseMatchPopupRef', 'pairwiseMatchPopupDrag', 'pairwiseMatchPopupSize', 'pairwiseMatchPopupResize',
  'featurePopupRef', 'featurePopupDrag', 'featurePopupSize', 'featurePopupResize', 'clickedLabel',
  'clickedLabelPos', 'featureStyleScopeDialog', 'resetColorDialog', 'legendRenameDialog',
  'labelTextScopeDialog', 'featureVisibilityScopeDialog', 'hiddenLabelTextDialog', 'labelOnDialog',
  'sidebarWidth', 'isResizing', 'newLegendCaption', 'newLegendColor', 'legendDragging', 'legendDragStart',
  'legendOriginalTransform', 'legendInitialTransform', 'diagramDragging', 'diagramDragStart',
  'diagramElementIds', 'diagramElementOriginalTransforms', 'diagramElements', 'lengthBarElement',
  'lengthBarOriginalTransform', 'plotTitleElement', 'plotTitleDragging', 'plotTitleDragStart',
  'plotTitleAutoTransform', 'showCanvasControls', 'shouldDeferCircularPreviewUpdates', 'skipCaptureBaseConfig',
  'skipExtractOnSvgChange', 'trustedArtifactRestoreInProgress', 'newColorFeat', 'newColorVal',
  'newPriorityRule', 'newFeatureToAdd', 'featureList', 'featureListState', 'filteredFeatures',
  'filteredEditableLabels'
]);

// The drawing store itself.
export const STORE_KEYS = Object.freeze(['drawings', 'activeDrawing']);

// A `state`-shaped node fixture's drawing: its own drawing members, by reference.
export const drawingOf = (fixture) => Object.freeze(Object.fromEntries(
  DRAWING_KEYS.filter((key) => Object.hasOwn(fixture, key)).map((key) => [key, fixture[key]])
));

// The fixture with the drawing store that services and owners read
// (`state.drawings`, `state.activeDrawing()`): the fixture is the one drawing
// of both modes, as in state.js, so a member the test replaces later is read.
// The store is not enumerable, so a test that clones or serializes the
// fixture sees only its members.
export const withDrawings = (fixture) => Object.defineProperties(fixture, {
  drawings: { value: Object.freeze({ circular: fixture, linear: fixture }), configurable: true },
  activeDrawing: { value: () => fixture, configurable: true }
});
