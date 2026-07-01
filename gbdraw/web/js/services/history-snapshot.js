import {
  buildFeatureVisibilityOverrideCache,
  featureVisibilityRulesFromOverrideCache,
  normalizeFeatureVisibilityRule
} from '../app/feature-visibility.js';

const cloneJsonData = (value) => {
  if (value === null || value === undefined) return value;
  return JSON.parse(JSON.stringify(value));
};

const clonePlainObject = (value) => {
  if (!value || typeof value !== 'object' || Array.isArray(value)) return {};
  return cloneJsonData(value);
};

const replacePlainObject = (target, source) => {
  if (!target || typeof target !== 'object') return;
  Object.keys(target).forEach((key) => delete target[key]);
  Object.entries(source || {}).forEach(([key, value]) => {
    target[key] = value;
  });
};

const cloneFeatureVisibilityRules = (rules) => (
  Array.isArray(rules) ? rules.map((rule) => normalizeFeatureVisibilityRule(rule)) : []
);

const replaceFeatureVisibilityRules = (state, rules) => {
  const normalizedRules = cloneFeatureVisibilityRules(rules);
  if (Array.isArray(state.featureVisibilityRules)) {
    state.featureVisibilityRules.splice(0, state.featureVisibilityRules.length, ...normalizedRules);
  }
  replacePlainObject(state.featureVisibilityOverrides, buildFeatureVisibilityOverrideCache(normalizedRules));
};

const setRef = (target, value) => {
  if (target && typeof target === 'object' && Object.prototype.hasOwnProperty.call(target, 'value')) {
    target.value = value;
  }
};

const getRef = (target, fallback = null) => (
  target && typeof target === 'object' && Object.prototype.hasOwnProperty.call(target, 'value')
    ? target.value
    : fallback
);

const nextFrame = () => new Promise((resolve) => {
  if (typeof window !== 'undefined' && typeof window.requestAnimationFrame === 'function') {
    window.requestAnimationFrame(() => resolve());
  } else {
    setTimeout(resolve, 0);
  }
});

const closeTransientState = (state) => {
  setRef(state.clickedFeature, null);
  setRef(state.clickedPairwiseMatch, null);
  setRef(state.clickedLabel, null);
  if (state.colorScopeDialog) state.colorScopeDialog.show = false;
  if (state.resetColorDialog) state.resetColorDialog.show = false;
  if (state.legendRenameDialog) state.legendRenameDialog.show = false;
  if (state.labelTextScopeDialog) state.labelTextScopeDialog.show = false;
  if (state.featureVisibilityScopeDialog) state.featureVisibilityScopeDialog.show = false;
  if (state.globalLabelModeDialog) state.globalLabelModeDialog.show = false;
  if (state.featurePopupDrag) state.featurePopupDrag.active = false;
  if (state.featurePopupResize) state.featurePopupResize.active = false;
  if (state.pairwiseMatchPopupDrag) state.pairwiseMatchPopupDrag.active = false;
  if (state.pairwiseMatchPopupResize) state.pairwiseMatchPopupResize.active = false;
};

const buildFallbackUiStateData = (state) => ({
  title: getRef(state.sessionTitle, ''),
  mode: getRef(state.mode, 'circular'),
  cInputType: getRef(state.cInputType, 'gb'),
  lInputType: getRef(state.lInputType, 'gb'),
  blastSource: getRef(state.blastSource, 'losat'),
  losatProgram: getRef(state.losatProgram, 'blastn'),
  selectedResultIndex: getRef(state.selectedResultIndex, 0),
  downloadDpi: getRef(state.downloadDpi, 300),
  canvasPadding: { ...(state.canvasPadding || {}) },
  generatedLegendPosition: getRef(state.generatedLegendPosition, 'left'),
  generatedMode: getRef(state.generatedMode, 'circular'),
  generatedMultiRecordCanvas: Boolean(getRef(state.generatedMultiRecordCanvas, false)),
  generatedCircularPlotTitlePosition: getRef(state.generatedCircularPlotTitlePosition, 'none'),
  circularLegendPosition: getRef(state.circularLegendPosition, 'left'),
  linearLegendPosition: getRef(state.linearLegendPosition, 'bottom'),
  circularPlotTitlePosition: getRef(state.circularPlotTitlePosition, 'none'),
  linearPlotTitlePosition: getRef(state.linearPlotTitlePosition, 'bottom'),
  circularSingleRecordLegendPosition: getRef(state.circularSingleRecordLegendPosition, 'left'),
  circularSingleRecordPlotTitlePosition: getRef(state.circularSingleRecordPlotTitlePosition, 'none'),
  circularMultiRecordLegendPosition: getRef(state.circularMultiRecordLegendPosition, null),
  circularMultiRecordPlotTitlePosition: getRef(state.circularMultiRecordPlotTitlePosition, null),
  autoLabelReflow: Boolean(getRef(state.autoLabelReflowEnabled, false)),
  paletteInstantPreviewEnabled: Boolean(getRef(state.paletteInstantPreviewEnabled, false)),
  appliedPaletteName: getRef(state.appliedPaletteName, 'default'),
  appliedPaletteColors: clonePlainObject(getRef(state.appliedPaletteColors, {})),
  pendingPaletteName: getRef(state.pendingPaletteName, ''),
  pendingPaletteColors: clonePlainObject(getRef(state.pendingPaletteColors, {})),
  circularBaseConfig: clonePlainObject(getRef(state.circularBaseConfig, {})),
  linearBaseConfig: {
    ...clonePlainObject(getRef(state.linearBaseConfig, {})),
    diagramBaseTransforms: []
  },
  legendCurrentOffset: { ...(state.legendCurrentOffset || {}) },
  diagramOffset: { ...(state.diagramOffset || {}) },
  lengthBarUserOffset: { ...(state.lengthBarUserOffset || {}) },
  plotTitleUserOffset: { ...(state.plotTitleUserOffset || {}) }
});

const applyFallbackUiStateData = (state, ui = {}) => {
  if (typeof ui.title === 'string') setRef(state.sessionTitle, ui.title);
  if (ui.mode) setRef(state.mode, ui.mode === 'linear' ? 'linear' : 'circular');
  if (ui.cInputType) setRef(state.cInputType, ui.cInputType);
  if (ui.lInputType) setRef(state.lInputType, ui.lInputType);
  if (ui.blastSource) setRef(state.blastSource, ui.blastSource);
  if (ui.losatProgram) setRef(state.losatProgram, ui.losatProgram);
  if (ui.downloadDpi) setRef(state.downloadDpi, ui.downloadDpi);
  if (ui.generatedLegendPosition) setRef(state.generatedLegendPosition, ui.generatedLegendPosition);
  if (ui.generatedMode) setRef(state.generatedMode, ui.generatedMode);
  if (Object.prototype.hasOwnProperty.call(ui, 'generatedMultiRecordCanvas')) {
    setRef(state.generatedMultiRecordCanvas, Boolean(ui.generatedMultiRecordCanvas));
  }
  if (ui.generatedCircularPlotTitlePosition) {
    setRef(state.generatedCircularPlotTitlePosition, ui.generatedCircularPlotTitlePosition);
  }
  if (ui.circularLegendPosition) setRef(state.circularLegendPosition, ui.circularLegendPosition);
  if (ui.linearLegendPosition) setRef(state.linearLegendPosition, ui.linearLegendPosition);
  if (ui.circularPlotTitlePosition) setRef(state.circularPlotTitlePosition, ui.circularPlotTitlePosition);
  if (ui.linearPlotTitlePosition) setRef(state.linearPlotTitlePosition, ui.linearPlotTitlePosition);
  if (Object.prototype.hasOwnProperty.call(ui, 'circularSingleRecordLegendPosition')) {
    setRef(state.circularSingleRecordLegendPosition, ui.circularSingleRecordLegendPosition);
  }
  if (Object.prototype.hasOwnProperty.call(ui, 'circularSingleRecordPlotTitlePosition')) {
    setRef(state.circularSingleRecordPlotTitlePosition, ui.circularSingleRecordPlotTitlePosition);
  }
  if (Object.prototype.hasOwnProperty.call(ui, 'circularMultiRecordLegendPosition')) {
    setRef(state.circularMultiRecordLegendPosition, ui.circularMultiRecordLegendPosition);
  }
  if (Object.prototype.hasOwnProperty.call(ui, 'circularMultiRecordPlotTitlePosition')) {
    setRef(state.circularMultiRecordPlotTitlePosition, ui.circularMultiRecordPlotTitlePosition);
  }
  if (state.canvasPadding && ui.canvasPadding) {
    state.canvasPadding.top = Number(ui.canvasPadding.top) || 0;
    state.canvasPadding.right = Number(ui.canvasPadding.right) || 0;
    state.canvasPadding.bottom = Number(ui.canvasPadding.bottom) || 0;
    state.canvasPadding.left = Number(ui.canvasPadding.left) || 0;
  }
  setRef(state.autoLabelReflowEnabled, Boolean(ui.autoLabelReflow));
  setRef(state.paletteInstantPreviewEnabled, Boolean(ui.paletteInstantPreviewEnabled));
  if (ui.appliedPaletteName !== undefined) setRef(state.appliedPaletteName, String(ui.appliedPaletteName || 'default'));
  if (ui.appliedPaletteColors) setRef(state.appliedPaletteColors, clonePlainObject(ui.appliedPaletteColors));
  if (ui.pendingPaletteName !== undefined) setRef(state.pendingPaletteName, String(ui.pendingPaletteName || ''));
  if (ui.pendingPaletteColors) setRef(state.pendingPaletteColors, clonePlainObject(ui.pendingPaletteColors));
  if (ui.circularBaseConfig) setRef(state.circularBaseConfig, clonePlainObject(ui.circularBaseConfig));
  if (ui.linearBaseConfig) {
    setRef(state.linearBaseConfig, {
      ...clonePlainObject(ui.linearBaseConfig),
      diagramBaseTransforms: new Map()
    });
  }
  if (state.legendCurrentOffset && ui.legendCurrentOffset) {
    state.legendCurrentOffset.x = Number(ui.legendCurrentOffset.x) || 0;
    state.legendCurrentOffset.y = Number(ui.legendCurrentOffset.y) || 0;
  }
  if (state.diagramOffset && ui.diagramOffset) {
    state.diagramOffset.x = Number(ui.diagramOffset.x) || 0;
    state.diagramOffset.y = Number(ui.diagramOffset.y) || 0;
  }
  if (state.lengthBarUserOffset && ui.lengthBarUserOffset) {
    state.lengthBarUserOffset.x = Number(ui.lengthBarUserOffset.x) || 0;
    state.lengthBarUserOffset.y = Number(ui.lengthBarUserOffset.y) || 0;
  }
  if (state.plotTitleUserOffset && ui.plotTitleUserOffset) {
    state.plotTitleUserOffset.x = Number(ui.plotTitleUserOffset.x) || 0;
    state.plotTitleUserOffset.y = Number(ui.plotTitleUserOffset.y) || 0;
  }
  if (Number.isInteger(ui.selectedResultIndex)) {
    const count = Array.isArray(getRef(state.results, [])) ? getRef(state.results, []).length : 0;
    setRef(state.selectedResultIndex, count > 0 ? Math.max(0, Math.min(ui.selectedResultIndex, count - 1)) : 0);
  }
};

const buildFallbackFeatureStateData = (state) => ({
  extractedFeatures: cloneJsonData(getRef(state.extractedFeatures, [])) || [],
  featureSelectorSafetyScope: cloneJsonData(getRef(state.featureSelectorSafetyScope, [])) || [],
  featureRecordIds: cloneJsonData(getRef(state.featureRecordIds, [])) || [],
  selectedFeatureRecordIdx: getRef(state.selectedFeatureRecordIdx, 0),
  featureColorOverrides: clonePlainObject(state.featureColorOverrides),
  featureVisibilityRules: cloneFeatureVisibilityRules(state.featureVisibilityRules),
  labelTextFeatureOverrides: clonePlainObject(state.labelTextFeatureOverrides),
  labelTextBulkOverrides: clonePlainObject(state.labelTextBulkOverrides),
  labelTextFeatureOverrideSources: clonePlainObject(state.labelTextFeatureOverrideSources),
  labelVisibilityOverrides: clonePlainObject(state.labelVisibilityOverrides),
  labelOverrideContextKey: getRef(state.labelOverrideContextKey, '')
});

const applyFallbackFeatureStateData = (state, features = {}) => {
  setRef(state.extractedFeatures, cloneJsonData(features.extractedFeatures) || []);
  setRef(state.featureSelectorSafetyScope, cloneJsonData(features.featureSelectorSafetyScope) || []);
  setRef(state.featureRecordIds, cloneJsonData(features.featureRecordIds) || []);
  setRef(
    state.selectedFeatureRecordIdx,
    Number.isInteger(features.selectedFeatureRecordIdx) ? features.selectedFeatureRecordIdx : 0
  );
  replacePlainObject(state.featureColorOverrides, clonePlainObject(features.featureColorOverrides));
  replaceFeatureVisibilityRules(
    state,
    Array.isArray(features.featureVisibilityRules)
      ? features.featureVisibilityRules
      : featureVisibilityRulesFromOverrideCache(features.featureVisibilityOverrides)
  );
  replacePlainObject(state.labelTextFeatureOverrides, clonePlainObject(features.labelTextFeatureOverrides));
  replacePlainObject(state.labelTextBulkOverrides, clonePlainObject(features.labelTextBulkOverrides));
  replacePlainObject(
    state.labelTextFeatureOverrideSources,
    clonePlainObject(features.labelTextFeatureOverrideSources)
  );
  replacePlainObject(state.labelVisibilityOverrides, clonePlainObject(features.labelVisibilityOverrides));
  setRef(state.labelOverrideContextKey, String(features.labelOverrideContextKey || ''));
};

const buildFallbackOrthogroupStateData = (state) => ({
  groups: cloneJsonData(getRef(state.orthogroups, [])) || [],
  selectedOrthogroupId: getRef(state.selectedOrthogroupId, ''),
  selectedOrthogroupAlignmentFeature: getRef(state.selectedOrthogroupAlignmentFeature, ''),
  orthogroupNameOverrides: clonePlainObject(state.orthogroupNameOverrides),
  orthogroupDescriptionOverrides: clonePlainObject(state.orthogroupDescriptionOverrides)
});

const applyFallbackOrthogroupStateData = (state, data = {}) => {
  setRef(state.orthogroups, cloneJsonData(data.groups) || []);
  setRef(state.selectedOrthogroupId, String(data.selectedOrthogroupId || ''));
  setRef(state.selectedOrthogroupAlignmentFeature, String(data.selectedOrthogroupAlignmentFeature || ''));
  replacePlainObject(state.orthogroupNameOverrides, clonePlainObject(data.orthogroupNameOverrides));
  replacePlainObject(state.orthogroupDescriptionOverrides, clonePlainObject(data.orthogroupDescriptionOverrides));
};

const buildFallbackResultsData = (state) => {
  const currentSvg = (() => {
    const svg = state.svgContainer?.value?.querySelector?.('svg');
    if (!svg || typeof XMLSerializer === 'undefined') return null;
    return new XMLSerializer().serializeToString(svg);
  })();
  const selected = getRef(state.selectedResultIndex, 0);
  return (getRef(state.results, []) || []).map((result, index) => ({
    name: result?.name || `Result ${index + 1}`,
    content: index === selected && currentSvg ? currentSvg : String(result?.content || '')
  }));
};

const applyFallbackResultsData = (state, results = []) => {
  setRef(
    state.results,
    Array.isArray(results)
      ? results.map((result, index) => ({
          name: result?.name || `Result ${index + 1}`,
          content: String(result?.content || '')
        }))
      : []
  );
};

const buildFilesData = (state, fileStore) => ({
  c_gb: fileStore.describeValue(state.files?.c_gb),
  c_gff: fileStore.describeValue(state.files?.c_gff),
  c_fasta: fileStore.describeValue(state.files?.c_fasta),
  c_depth: fileStore.describeValue(state.files?.c_depth),
  c_conservation_blasts: fileStore.describeValue(state.files?.c_conservation_blasts || []),
  c_conservation_fastas: fileStore.describeValue(state.files?.c_conservation_fastas || []),
  d_color: fileStore.describeValue(state.files?.d_color),
  t_color: fileStore.describeValue(state.files?.t_color),
  blacklist: fileStore.describeValue(state.files?.blacklist),
  whitelist: fileStore.describeValue(state.files?.whitelist),
  qualifier_priority: fileStore.describeValue(state.files?.qualifier_priority),
  linearSeqs: Array.from(state.linearSeqs || []).map((seq) => ({
    uid: seq.uid,
    gb: fileStore.describeValue(seq.gb),
    gff: fileStore.describeValue(seq.gff),
    fasta: fileStore.describeValue(seq.fasta),
    depth: fileStore.describeValue(seq.depth),
    blast: fileStore.describeValue(seq.blast),
    losat_gencode: seq.losat_gencode ?? 1,
    losat_filename: seq.losat_filename ?? '',
    definition: seq.definition ?? '',
    record_subtitle: seq.record_subtitle ?? '',
    region_record_id: seq.region_record_id ?? '',
    region_start: seq.region_start ?? null,
    region_end: seq.region_end ?? null,
    region_reverse: Boolean(seq.region_reverse)
  }))
});

const applyFilesData = (state, filesData, fileStore, normalizeLinearSeqList = null) => {
  if (!state.files) return;
  const restore = (value) => fileStore.restoreValue(value);
  state.files.c_gb = restore(filesData?.c_gb);
  state.files.c_gff = restore(filesData?.c_gff);
  state.files.c_fasta = restore(filesData?.c_fasta);
  state.files.c_depth = restore(filesData?.c_depth);
  state.files.c_conservation_blasts = Array.isArray(filesData?.c_conservation_blasts)
    ? restore(filesData.c_conservation_blasts).filter(Boolean)
    : [];
  state.files.c_conservation_fastas = Array.isArray(filesData?.c_conservation_fastas)
    ? restore(filesData.c_conservation_fastas).filter(Boolean)
    : [];
  state.files.d_color = restore(filesData?.d_color);
  state.files.t_color = restore(filesData?.t_color);
  state.files.blacklist = restore(filesData?.blacklist);
  state.files.whitelist = restore(filesData?.whitelist);
  state.files.qualifier_priority = restore(filesData?.qualifier_priority);

  if (!state.linearSeqs || typeof state.linearSeqs.splice !== 'function') return;
  const rows = Array.isArray(filesData?.linearSeqs)
    ? filesData.linearSeqs.map((seq) => ({
        uid: seq.uid,
        gb: restore(seq.gb),
        gff: restore(seq.gff),
        fasta: restore(seq.fasta),
        depth: restore(seq.depth),
        blast: restore(seq.blast),
        losat_gencode: seq.losat_gencode ?? 1,
        losat_filename: seq.losat_filename ?? '',
        definition: seq.definition ?? '',
        record_subtitle: seq.record_subtitle ?? '',
        region_record_id: seq.region_record_id ?? '',
        region_start: seq.region_start ?? null,
        region_end: seq.region_end ?? null,
        region_reverse: Boolean(seq.region_reverse)
      }))
    : [];
  const normalized = typeof normalizeLinearSeqList === 'function'
    ? normalizeLinearSeqList(rows)
    : rows;
  state.linearSeqs.splice(0, state.linearSeqs.length, ...normalized);
};

export const createHistorySnapshotService = ({
  state,
  fileStore,
  nextTick = async () => {},
  normalizeLinearSeqList = null,
  buildConfigData = null,
  applyConfigData = null,
  buildUiStateData = null,
  applyUiStateData = null,
  buildFeatureStateData = null,
  applyFeatureStateData = null,
  buildEditorStateData = null,
  applyEditorStateData = null,
  buildOrthogroupStateData = null,
  applyOrthogroupStateData = null,
  serializeResults = null,
  applyResultsData = null,
  buildRunStateData = null,
  applyRunStateData = null
}) => {
  if (!state || !fileStore) {
    throw new Error('createHistorySnapshotService requires state and fileStore.');
  }

  const buildHistorySnapshot = async () => {
    const config = typeof buildConfigData === 'function'
      ? buildConfigData()
      : { form: state.form, adv: state.adv };
    const ui = typeof buildUiStateData === 'function'
      ? buildUiStateData({ includePreviewNavigation: false })
      : buildFallbackUiStateData(state);
    const features = typeof buildFeatureStateData === 'function'
      ? buildFeatureStateData()
      : buildFallbackFeatureStateData(state);
    const editorState = typeof buildEditorStateData === 'function'
      ? buildEditorStateData()
      : {};
    const orthogroupState = typeof buildOrthogroupStateData === 'function'
      ? buildOrthogroupStateData()
      : buildFallbackOrthogroupStateData(state);
    const results = typeof serializeResults === 'function'
      ? serializeResults()
      : buildFallbackResultsData(state);
    const runState = typeof buildRunStateData === 'function'
      ? buildRunStateData()
      : {
          lastRunInfo: cloneJsonData(getRef(state.lastRunInfo, null)),
          pairwiseMatchFactors: clonePlainObject(getRef(state.pairwiseMatchFactors, {}))
        };

    return cloneJsonData({
      config,
      ui,
      files: buildFilesData(state, fileStore),
      results,
      features,
      editorState,
      orthogroupState,
      runState
    });
  };

  const applyHistorySnapshot = async (snapshot) => {
    if (!snapshot || typeof snapshot !== 'object') return;
    closeTransientState(state);

    const ui = snapshot.ui || {};
    if (ui.mode) setRef(state.mode, ui.mode === 'linear' ? 'linear' : 'circular');
    if (ui.cInputType) setRef(state.cInputType, ui.cInputType);
    if (ui.lInputType) setRef(state.lInputType, ui.lInputType);
    // The mode watcher clears generated metadata. Let that reset finish before
    // restoring snapshot-owned feature, label, and orthogroup state.
    await nextTick();

    if (typeof applyConfigData === 'function' && snapshot.config) {
      applyConfigData(snapshot.config);
    }

    if (typeof applyUiStateData === 'function') {
      applyUiStateData(ui, { restorePreviewNavigation: false });
    } else {
      applyFallbackUiStateData(state, ui);
    }
    await nextTick();

    applyFilesData(state, snapshot.files || {}, fileStore, normalizeLinearSeqList);

    if (state.skipCaptureBaseConfig) state.skipCaptureBaseConfig.value = true;
    if (state.skipPositionReapply) state.skipPositionReapply.value = true;
    if (state.skipExtractOnSvgChange) state.skipExtractOnSvgChange.value = false;

    if (typeof applyResultsData === 'function') {
      applyResultsData(snapshot.results || [], ui);
    } else {
      applyFallbackResultsData(state, snapshot.results || []);
      applyFallbackUiStateData(state, { selectedResultIndex: ui.selectedResultIndex });
    }

    if (typeof applyFeatureStateData === 'function') {
      applyFeatureStateData(snapshot.features || {});
    } else {
      applyFallbackFeatureStateData(state, snapshot.features || {});
    }

    if (typeof applyOrthogroupStateData === 'function') {
      applyOrthogroupStateData(snapshot.orthogroupState || {});
    } else {
      applyFallbackOrthogroupStateData(state, snapshot.orthogroupState || {});
    }

    if (typeof applyEditorStateData === 'function') {
      applyEditorStateData(snapshot.editorState || {});
    }

    if (typeof applyRunStateData === 'function') {
      applyRunStateData(snapshot.runState || {});
    } else {
      setRef(state.lastRunInfo, cloneJsonData(snapshot.runState?.lastRunInfo) || null);
      setRef(state.pairwiseMatchFactors, clonePlainObject(snapshot.runState?.pairwiseMatchFactors));
    }

    await nextTick();
    await nextFrame();
    if (typeof applyUiStateData === 'function') {
      applyUiStateData(ui, { restorePreviewNavigation: false });
    } else {
      applyFallbackUiStateData(state, ui);
    }
  };

  const snapshotSignature = (snapshot) => JSON.stringify(snapshot);

  return {
    applyHistorySnapshot,
    buildHistorySnapshot,
    snapshotSignature
  };
};
