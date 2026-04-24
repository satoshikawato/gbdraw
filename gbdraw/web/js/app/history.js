const { ref, computed } = window.Vue;

const HISTORY_LIMIT = 20;

const clonePlain = (value, fallback) => {
  if (value === undefined) return fallback;
  return JSON.parse(JSON.stringify(value));
};

const replaceReactiveObject = (target, source = {}) => {
  Object.keys(target).forEach((key) => delete target[key]);
  Object.entries(source || {}).forEach(([key, value]) => {
    target[key] = value;
  });
};

const replaceReactiveArray = (target, items = []) => {
  target.splice(0, target.length, ...(Array.isArray(items) ? items : []));
};

const replaceRefArray = (targetRef, items = []) => {
  targetRef.value = Array.isArray(items) ? items : [];
};

const serializeCurrentResults = (state) => {
  const currentSvg = (() => {
    if (!state.svgContainer.value) return null;
    const svg = state.svgContainer.value.querySelector('svg');
    if (!svg) return null;
    const serializer = new XMLSerializer();
    return serializer.serializeToString(svg);
  })();

  return state.results.value.map((result, index) => ({
    name: String(result?.name || `Result ${index + 1}`),
    content:
      index === state.selectedResultIndex.value && typeof currentSvg === 'string'
        ? currentSvg
        : String(result?.content || '')
  }));
};

const buildHistorySnapshot = (state) => ({
  results: serializeCurrentResults(state),
  selectedResultIndex: state.selectedResultIndex.value,
  manualSpecificRules: clonePlain(state.manualSpecificRules, []),
  featureColorOverrides: clonePlain(state.featureColorOverrides, {}),
  featureVisibilityOverrides: clonePlain(state.featureVisibilityOverrides, {}),
  labelTextFeatureOverrides: clonePlain(state.labelTextFeatureOverrides, {}),
  labelTextBulkOverrides: clonePlain(state.labelTextBulkOverrides, {}),
  labelTextFeatureOverrideSources: clonePlain(state.labelTextFeatureOverrideSources, {}),
  labelVisibilityOverrides: clonePlain(state.labelVisibilityOverrides, {}),
  labelOverrideContextKey: String(state.labelOverrideContextKey.value || ''),
  deletedLegendEntries: clonePlain(state.deletedLegendEntries.value, []),
  originalLegendOrder: clonePlain(state.originalLegendOrder.value, []),
  originalLegendColors: clonePlain(state.originalLegendColors.value, {}),
  legendColorOverrides: clonePlain(state.legendColorOverrides, {}),
  legendStrokeOverrides: clonePlain(state.legendStrokeOverrides, {}),
  originalSvgStroke: clonePlain(state.originalSvgStroke.value, { color: null, width: null }),
  addedLegendCaptions: Array.from(state.addedLegendCaptions.value || []),
  fileLegendCaptions: Array.from(state.fileLegendCaptions.value || []),
  canvasPadding: {
    top: Number(state.canvasPadding.top) || 0,
    right: Number(state.canvasPadding.right) || 0,
    bottom: Number(state.canvasPadding.bottom) || 0,
    left: Number(state.canvasPadding.left) || 0
  },
  legendCurrentOffset: {
    x: Number(state.legendCurrentOffset.x) || 0,
    y: Number(state.legendCurrentOffset.y) || 0
  },
  legendInitialTransform: clonePlain(state.legendInitialTransform.value, { x: 0, y: 0 }),
  diagramOffset: {
    x: Number(state.diagramOffset.x) || 0,
    y: Number(state.diagramOffset.y) || 0
  },
  plotTitleAutoTransform: clonePlain(state.plotTitleAutoTransform.value, { x: 0, y: 0 }),
  plotTitleUserOffset: {
    x: Number(state.plotTitleUserOffset.x) || 0,
    y: Number(state.plotTitleUserOffset.y) || 0
  },
  generatedLegendPosition: String(state.generatedLegendPosition.value || ''),
  generatedMode: String(state.generatedMode.value || ''),
  generatedMultiRecordCanvas: Boolean(state.generatedMultiRecordCanvas.value),
  generatedCircularPlotTitlePosition: String(state.generatedCircularPlotTitlePosition.value || '')
});

const normalizeResultIndex = (results, requestedIndex) => {
  const count = Array.isArray(results) ? results.length : 0;
  if (count <= 0) return 0;
  const numeric = Number(requestedIndex);
  if (!Number.isInteger(numeric) || numeric < 0) return 0;
  return Math.min(numeric, count - 1);
};

export const createHistoryManager = ({ state, nextTick }) => {
  const undoStack = ref([]);
  const redoStack = ref([]);

  const setRestoring = (value) => {
    state.historyRestoring.value = Boolean(value);
  };

  const canUndo = computed(() => undoStack.value.length > 1);
  const canRedo = computed(() => redoStack.value.length > 0);

  const clearHistory = () => {
    undoStack.value = [];
    redoStack.value = [];
  };

  const buildEntry = (label = '') => {
    const snapshot = buildHistorySnapshot(state);
    return {
      label: String(label || ''),
      signature: JSON.stringify(snapshot),
      snapshot
    };
  };

  const waitForRender = async () => {
    await nextTick();
    await Promise.resolve();
    await nextTick();
  };

  const pushCheckpoint = async (label = '') => {
    if (state.historyRestoring.value) return false;
    if (!Array.isArray(state.results.value) || state.results.value.length === 0) return false;

    await waitForRender();

    const entry = buildEntry(label);
    const previous = undoStack.value[undoStack.value.length - 1];
    if (previous && previous.signature === entry.signature) {
      return false;
    }

    undoStack.value.push(entry);
    if (undoStack.value.length > HISTORY_LIMIT) {
      undoStack.value.splice(1, undoStack.value.length - HISTORY_LIMIT);
    }
    redoStack.value = [];
    return true;
  };

  const resetHistory = async (label = '') => {
    clearHistory();
    if (!Array.isArray(state.results.value) || state.results.value.length === 0) return false;
    await waitForRender();
    undoStack.value = [buildEntry(label)];
    return true;
  };

  const restoreSnapshot = async (snapshot) => {
    if (!snapshot) return false;

    setRestoring(true);
    try {
      state.clickedFeature.value = null;
      state.clickedLabel.value = null;
      state.skipCaptureBaseConfig.value = true;
      state.skipPositionReapply.value = true;

      replaceReactiveArray(state.manualSpecificRules, clonePlain(snapshot.manualSpecificRules, []));
      replaceReactiveObject(state.featureColorOverrides, clonePlain(snapshot.featureColorOverrides, {}));
      replaceReactiveObject(state.featureVisibilityOverrides, clonePlain(snapshot.featureVisibilityOverrides, {}));
      replaceReactiveObject(state.labelTextFeatureOverrides, clonePlain(snapshot.labelTextFeatureOverrides, {}));
      replaceReactiveObject(state.labelTextBulkOverrides, clonePlain(snapshot.labelTextBulkOverrides, {}));
      replaceReactiveObject(
        state.labelTextFeatureOverrideSources,
        clonePlain(snapshot.labelTextFeatureOverrideSources, {})
      );
      replaceReactiveObject(state.labelVisibilityOverrides, clonePlain(snapshot.labelVisibilityOverrides, {}));
      state.labelOverrideContextKey.value = String(snapshot.labelOverrideContextKey || '');

      replaceRefArray(state.deletedLegendEntries, clonePlain(snapshot.deletedLegendEntries, []));
      replaceRefArray(state.originalLegendOrder, clonePlain(snapshot.originalLegendOrder, []));
      state.originalLegendColors.value = clonePlain(snapshot.originalLegendColors, {});
      replaceReactiveObject(state.legendColorOverrides, clonePlain(snapshot.legendColorOverrides, {}));
      replaceReactiveObject(state.legendStrokeOverrides, clonePlain(snapshot.legendStrokeOverrides, {}));
      state.originalSvgStroke.value = clonePlain(snapshot.originalSvgStroke, { color: null, width: null });
      state.addedLegendCaptions.value = new Set(snapshot.addedLegendCaptions || []);
      state.fileLegendCaptions.value = new Set(snapshot.fileLegendCaptions || []);

      state.legendCurrentOffset.x = Number(snapshot.legendCurrentOffset?.x) || 0;
      state.legendCurrentOffset.y = Number(snapshot.legendCurrentOffset?.y) || 0;
      state.diagramOffset.x = Number(snapshot.diagramOffset?.x) || 0;
      state.diagramOffset.y = Number(snapshot.diagramOffset?.y) || 0;
      state.plotTitleUserOffset.x = Number(snapshot.plotTitleUserOffset?.x) || 0;
      state.plotTitleUserOffset.y = Number(snapshot.plotTitleUserOffset?.y) || 0;
      state.plotTitleAutoTransform.value = clonePlain(snapshot.plotTitleAutoTransform, { x: 0, y: 0 });

      state.generatedLegendPosition.value = String(snapshot.generatedLegendPosition || state.generatedLegendPosition.value);
      state.generatedMode.value = String(snapshot.generatedMode || state.generatedMode.value);
      state.generatedMultiRecordCanvas.value = Boolean(snapshot.generatedMultiRecordCanvas);
      state.generatedCircularPlotTitlePosition.value = String(
        snapshot.generatedCircularPlotTitlePosition || state.generatedCircularPlotTitlePosition.value
      );

      state.results.value = Array.isArray(snapshot.results)
        ? snapshot.results.map((result, index) => ({
            name: String(result?.name || `Result ${index + 1}`),
            content: String(result?.content || '')
          }))
        : [];
      state.selectedResultIndex.value = normalizeResultIndex(state.results.value, snapshot.selectedResultIndex);

      state.canvasPadding.top = Number(snapshot.canvasPadding?.top) || 0;
      state.canvasPadding.right = Number(snapshot.canvasPadding?.right) || 0;
      state.canvasPadding.bottom = Number(snapshot.canvasPadding?.bottom) || 0;
      state.canvasPadding.left = Number(snapshot.canvasPadding?.left) || 0;

      await waitForRender();
      state.legendInitialTransform.value = clonePlain(snapshot.legendInitialTransform, { x: 0, y: 0 });
      return true;
    } finally {
      await waitForRender();
      setRestoring(false);
    }
  };

  const undo = async () => {
    if (!canUndo.value || state.historyRestoring.value) return false;
    const current = undoStack.value.pop();
    if (!current) return false;
    redoStack.value.push(current);
    const previous = undoStack.value[undoStack.value.length - 1];
    return restoreSnapshot(previous?.snapshot);
  };

  const redo = async () => {
    if (!canRedo.value || state.historyRestoring.value) return false;
    const nextEntry = redoStack.value.pop();
    if (!nextEntry) return false;
    undoStack.value.push(nextEntry);
    return restoreSnapshot(nextEntry.snapshot);
  };

  return {
    canUndo,
    canRedo,
    clearHistory,
    pushCheckpoint,
    redo,
    resetHistory,
    undo
  };
};
