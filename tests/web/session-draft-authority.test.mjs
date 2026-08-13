import assert from 'node:assert/strict';
import { webcrypto } from 'node:crypto';
import { readFile } from 'node:fs/promises';
import { installFakeSvgDom } from './fake-svg-dom.mjs';

globalThis.window = {
  Vue: {
    ref: (value) => ({ value }),
    reactive: (value) => value,
    computed: (getter) => ({ get value() { return getter(); } }),
    nextTick: async () => {}
  },
  DOMPurify: { sanitize: (value) => value }
};
globalThis.document = {};
globalThis.crypto = webcrypto;
installFakeSvgDom();
globalThis.File = class File extends Blob {
  constructor(parts, name, options = {}) {
    super(parts, options);
    this.name = String(name || 'file');
    this.lastModified = options.lastModified ?? Date.now();
  }
};
const alerts = [];
globalThis.alert = (message) => alerts.push(String(message));

const {
  applyEditorStateData,
  applyFeatureStateData,
  buildConfigData,
  buildEditorStateData,
  buildFeatureStateData,
  buildOrthogroupStateData,
  buildRunStateData,
  buildUiStateData,
  exportSession,
  importSession,
  overlayCurrentWriterDraftConfig,
  validateCurrentWriterActiveDraft
} = await import('../../gbdraw/web/js/services/config.js');
const { catalogResultKey } = await import(
  '../../gbdraw/web/js/services/feature-catalog.js'
);
const { captureStyleSnapshot } = await import(
  '../../gbdraw/web/js/services/style-revision.js'
);
const { state } = await import('../../gbdraw/web/js/state.js');

state.legendOrderIntent.value = ['Beta', 'Alpha'];
const legendOrderEditorSnapshot = buildEditorStateData();
assert.deepEqual(legendOrderEditorSnapshot.legend.orderIntent, ['Beta', 'Alpha']);
state.legendOrderIntent.value = [];
applyEditorStateData(legendOrderEditorSnapshot, { trusted: true });
assert.deepEqual(state.legendOrderIntent.value, ['Beta', 'Alpha']);
const legacyLegendEditorSnapshot = structuredClone(legendOrderEditorSnapshot);
delete legacyLegendEditorSnapshot.legend.orderIntent;
applyEditorStateData(legacyLegendEditorSnapshot, { trusted: true });
assert.deepEqual(
  state.legendOrderIntent.value,
  [],
  'sessions written before document-global order intent hydrate compatibly'
);
const {
  COMPOSITION_METADATA_ATTRIBUTE,
  COMPOSITION_SCHEMA_ATTRIBUTE
} = await import('../../gbdraw/web/js/app/legend-layout/composition-actions.js');
const { createLegendCanvasActions } = await import(
  '../../gbdraw/web/js/app/legend-layout/canvas-actions.js'
);
const { createSessionImportRollbackState } = await import(
  '../../gbdraw/web/js/app/app-setup.js'
);
const { runRecordDiscoveryWatcher } = await import(
  '../../gbdraw/web/js/app/watchers.js'
);

const canonicalFeature = {
  id: 'features',
  renderer: 'features',
  enabled: true,
  side: 'inside',
  width: null,
  radius: null,
  inner_gap_px: null,
  outer_gap_px: null,
  z: 0,
  params: { lane_direction: 'split' }
};
const disabledDraft = {
  id: 'disabled-draft',
  renderer: 'depth',
  enabled: false,
  side: 'outside',
  width: '27px',
  radius: '1.2',
  inner_gap_px: '4px',
  outer_gap_px: '5px',
  z: 3,
  params: { track_index: 99, nested: { keep: true } }
};
const projectedConfig = {
  form: { track_type: 'tuckin' },
  adv: {
    nt: 'GC',
    circular_track_slots_enabled: true,
    circular_track_slots: [canonicalFeature],
    circular_track_slots_axis_index: 1,
    linear_track_slots_enabled: false,
    linear_track_slots: [],
    linear_track_slots_axis_index: null
  }
};
const storedConfig = {
  form: { track_type: 'tuckin' },
  adv: {
    ...projectedConfig.adv,
    circular_track_slots: [disabledDraft, canonicalFeature],
    circular_track_slots_axis_index: 2,
    linear_track_slots_enabled: false,
    linear_track_slots: [{
      id: 'inactive-linear',
      renderer: 'spacer',
      enabled: false,
      side: 'below',
      height: '8px',
      spacing: '2px',
      z: 0,
      params: {}
    }],
    linear_track_slots_axis_index: 1,
    feature_width_circular: 19,
    depth_width_circular: 23
  }
};

assert.doesNotThrow(() => validateCurrentWriterActiveDraft({
  mode: 'circular',
  projectedConfig,
  storedConfig
}));
const restored = overlayCurrentWriterDraftConfig(projectedConfig, storedConfig);
assert.deepEqual(restored.adv.circular_track_slots, storedConfig.adv.circular_track_slots);
assert.deepEqual(restored.adv.linear_track_slots, storedConfig.adv.linear_track_slots);
assert.equal(restored.adv.circular_track_slots_axis_index, 2);
assert.equal(restored.adv.feature_width_circular, 19);
assert.equal(restored.adv.depth_width_circular, 23);

const mismatched = structuredClone(storedConfig);
mismatched.adv.circular_track_slots[1].width = '12px';
assert.throws(
  () => validateCurrentWriterActiveDraft({
    mode: 'circular',
    projectedConfig,
    storedConfig: mismatched
  }),
  /does not match the committed render request/
);

const inactive = structuredClone(mismatched);
inactive.adv.circular_track_slots_enabled = false;
assert.doesNotThrow(() => validateCurrentWriterActiveDraft({
  mode: 'circular',
  projectedConfig,
  storedConfig: inactive
}));

const divergentSession = JSON.parse(await readFile(
  'gbdraw/web/gallery/sessions/BGC0000708-BGC0000713.gbdraw-session.json',
  'utf8'
));
const committedComparisonCount = divergentSession.renderRequest.comparisons.length;
assert.ok(committedComparisonCount > 0);
divergentSession.config.linearComparisonPlan = {
  mode: 'none',
  defaultSource: 'losat',
  edges: []
};
divergentSession.config.losatProgram = 'tblastx';
Object.assign(divergentSession.config.adv, {
  comparison_height: 73,
  min_bitscore: 123,
  evalue: '1e-37',
  identity: 88,
  alignment_length: 456,
  pairwise_match_style: 'ribbon'
});
Object.assign(divergentSession.config.losat.blastp, {
  mode: 'collinear',
  maxHits: 17,
  collinearMinAnchors: 4,
  collinearMaxUnitGap: 9,
  collinearSearchScope: 'adjacent',
  collinearColorMode: 'orientation_identity'
});
const dormantComparisonDraft = structuredClone({
  losatProgram: divergentSession.config.losatProgram,
  adv: {
    comparison_height: divergentSession.config.adv.comparison_height,
    min_bitscore: divergentSession.config.adv.min_bitscore,
    evalue: divergentSession.config.adv.evalue,
    identity: divergentSession.config.adv.identity,
    alignment_length: divergentSession.config.adv.alignment_length,
    pairwise_match_style: divergentSession.config.adv.pairwise_match_style
  },
  blastp: {
    mode: divergentSession.config.losat.blastp.mode,
    maxHits: divergentSession.config.losat.blastp.maxHits,
    collinearMinAnchors:
      divergentSession.config.losat.blastp.collinearMinAnchors,
    collinearMaxUnitGap:
      divergentSession.config.losat.blastp.collinearMaxUnitGap,
    collinearSearchScope:
      divergentSession.config.losat.blastp.collinearSearchScope,
    collinearColorMode:
      divergentSession.config.losat.blastp.collinearColorMode
  }
});
divergentSession.renderRequest.diagramOptions.output.legend = 'right';
divergentSession.config.form.legend = 'right';
divergentSession.ui.generatedLegendPosition = 'right';
divergentSession.ui.layoutPreferences = {
  circular: {
    single: { legend: 'left', plotTitlePosition: 'none' },
    multi: { legend: null, plotTitlePosition: null }
  },
  linear: { legend: 'bottom', plotTitlePosition: 'bottom' }
};
divergentSession.ui.paletteInstantPreviewEnabled = false;
divergentSession.ui.pendingPaletteName = 'mint';
divergentSession.ui.pendingPaletteColors = { CDS: '#abcdef' };
const importEvent = {
  target: {
    files: [new Blob([JSON.stringify(divergentSession)], { type: 'application/json' })],
    value: 'selected'
  }
};

const imported = await importSession(importEvent);

assert.equal(imported.status, 'ok');
const currentHashBearingSession = structuredClone(imported.data);
assert.equal(imported.data.version, 41);
assert.equal(imported.data.renderRequest.schema, 6);
assert.deepEqual(
  imported.data.renderRequest.diagramOptions.colors.defaultColorsFile,
  { resourceId: 'colors-default-colors-file', representation: 'file' }
);
assert.equal(imported.data.renderRequest.diagramOptions.colors.defaultColors, null);
assert.ok(
  imported.data.editorState.featureCatalog.items.every((item) => (
    item.biologicalFeatures.every((feature) => /^fi1_[a-z2-7]{26}$/.test(feature.instanceHash))
  )),
  'safe v40 canonical Feature pairs must be upgraded before current-session hydration'
);
assert.ok(
  Object.keys(imported.data.features.featureColorOverrides).every((key) => key.includes('\u0000')),
  'recovered rule-derived overrides must use canonical Feature pairs'
);
assert.equal(state.featureExactScopeAvailable.value, true);
assert.equal(state.featureExactScopeDiagnostic.value, '');
const mismatchedPairSession = structuredClone(currentHashBearingSession);
const mismatchedPairFeature = mismatchedPairSession.editorState.featureCatalog.items
  .flatMap((item) => item.biologicalFeatures)
  .find((feature) => feature.instanceHash);
assert.ok(mismatchedPairFeature);
mismatchedPairFeature.instanceHash = mismatchedPairFeature.instanceHash ===
  'fi1_aaaaaaaaaaaaaaaaaaaaaaaaaa'
  ? 'fi1_bbbbbbbbbbbbbbbbbbbbbbbbbb'
  : 'fi1_aaaaaaaaaaaaaaaaaaaaaaaaaa';
const epochBeforeMismatchedPairImport = state.documentEpoch.value;
const alertCountBeforeMismatchedPairImport = alerts.length;
const mismatchedPairEvent = {
  target: {
    files: [new Blob([JSON.stringify(mismatchedPairSession)], { type: 'application/json' })],
    value: 'selected'
  }
};
const consoleErrorBeforeMismatchedPairImport = console.error;
console.error = () => {};
let mismatchedPairImport;
try {
  mismatchedPairImport = await importSession(mismatchedPairEvent);
} finally {
  console.error = consoleErrorBeforeMismatchedPairImport;
}
assert.equal(mismatchedPairImport.status, 'error');
assert.match(
  mismatchedPairImport.error.message,
  /does not match its canonical record and Feature identity/
);
assert.equal(state.documentEpoch.value, epochBeforeMismatchedPairImport);
assert.equal(alerts.length, alertCountBeforeMismatchedPairImport + 1);
assert.match(
  alerts.at(-1),
  /Failed to load session: Feature instance hash does not match its canonical record and Feature identity/
);
assert.equal(mismatchedPairEvent.target.value, '');
const exactCapabilitySnapshot = buildFeatureStateData();
state.featureExactScopeAvailable.value = false;
state.featureExactScopeDiagnostic.value = 'temporary diagnostic';
applyFeatureStateData(exactCapabilitySnapshot);
assert.equal(state.featureExactScopeAvailable.value, true);
assert.equal(state.featureExactScopeDiagnostic.value, '');
assert.equal(state.linearComparisonPlan.mode, 'none');
assert.deepEqual(state.linearComparisonPlan.edges, []);
assert.deepEqual({
  losatProgram: state.losatProgram.value,
  adv: {
    comparison_height: state.adv.comparison_height,
    min_bitscore: state.adv.min_bitscore,
    evalue: state.adv.evalue,
    identity: state.adv.identity,
    alignment_length: state.adv.alignment_length,
    pairwise_match_style: state.adv.pairwise_match_style
  },
  blastp: {
    mode: state.losat.blastp.mode,
    maxHits: state.losat.blastp.maxHits,
    collinearMinAnchors: state.losat.blastp.collinearMinAnchors,
    collinearMaxUnitGap: state.losat.blastp.collinearMaxUnitGap,
    collinearSearchScope: state.losat.blastp.collinearSearchScope,
    collinearColorMode: state.losat.blastp.collinearColorMode
  }
}, dormantComparisonDraft);
assert.equal(
  state.form.legend,
  'bottom',
  'the active editor preference must override the last generated request position'
);
assert.equal(state.adv.plot_title_position, 'bottom');
assert.equal(state.generatedLegendPosition.value, 'right');
assert.equal(state.appliedPaletteName.value, 'orange');
assert.equal(state.appliedPaletteColors.value.CDS, '#dddddd');
assert.equal(state.pendingPaletteName.value, 'mint');
assert.equal(state.pendingPaletteColors.value.CDS, '#abcdef');
assert.deepEqual(state.layoutPreferences.linear, {
  legend: 'bottom',
  plotTitlePosition: 'bottom'
});
assert.equal(imported.data.renderRequest.diagramOptions.output.legend, 'right');
assert.equal(
  imported.data.renderRequest.comparisons.length,
  committedComparisonCount,
  'the editable comparison draft must not replace the last committed render request'
);

const legacyLayoutFields = [
  'legend',
  'circularLegendPosition',
  'linearLegendPosition',
  'circularPlotTitlePosition',
  'linearPlotTitlePosition',
  'circularSingleRecordLegendPosition',
  'circularSingleRecordPlotTitlePosition',
  'circularMultiRecordLegendPosition',
  'circularMultiRecordPlotTitlePosition'
];
const withoutStoredLayoutPreferences = () => {
  const payload = structuredClone(divergentSession);
  delete payload.ui.layoutPreferences;
  legacyLayoutFields.forEach((field) => delete payload.ui[field]);
  return payload;
};
const importPayload = async (payload) => importSession({
  target: {
    files: [new Blob([JSON.stringify(payload)], { type: 'application/json' })],
    value: 'selected'
  }
});

for (const savedPlan of [{
  mode: 'adjacent',
  defaultSource: 'upload',
  edges: []
}, {
  mode: 'selected',
  defaultSource: 'losat',
  edges: [{
    id: 'saved-selected-edge',
    queryUid: 'record-1',
    subjectUid: 'record-2',
    included: true,
    fileActive: false,
    losatFilenameActive: true,
    source: 'losat',
    losatFilename: 'saved-selected.raw.tsv'
  }]
}]) {
  const savedPlanSession = structuredClone(divergentSession);
  savedPlanSession.config.linearComparisonPlan = structuredClone(savedPlan);
  const savedPlanImport = await importPayload(savedPlanSession);
  assert.equal(savedPlanImport.status, 'ok');
  assert.deepEqual(state.linearComparisonPlan, {
    ...savedPlan,
    edges: savedPlan.edges.map((edge) => ({ ...edge, file: null }))
  });
  assert.equal(
    savedPlanImport.data.renderRequest.comparisons.length,
    committedComparisonCount
  );
}

const absentLayoutSession = withoutStoredLayoutPreferences();
absentLayoutSession.ui.generatedLegendPosition = 'right';
const absentLayoutImport = await importPayload(absentLayoutSession);
assert.equal(absentLayoutImport.status, 'ok');
assert.equal(state.form.legend, 'right');
assert.equal(state.generatedLegendPosition.value, 'right');
assert.equal(state.layoutPreferences.linear.legend, 'right');

const legacyLayoutSession = withoutStoredLayoutPreferences();
legacyLayoutSession.ui.legend = 'bottom';
legacyLayoutSession.ui.linearPlotTitlePosition = 'bottom';
legacyLayoutSession.ui.generatedLegendPosition = 'right';
const legacyLayoutImport = await importPayload(legacyLayoutSession);
assert.equal(legacyLayoutImport.status, 'ok');
assert.equal(state.form.legend, 'bottom');
assert.equal(state.generatedLegendPosition.value, 'right');
assert.equal(state.layoutPreferences.linear.legend, 'bottom');

const partialLayoutSession = withoutStoredLayoutPreferences();
partialLayoutSession.ui.layoutPreferences = {
  linear: { legend: 'top' }
};
partialLayoutSession.ui.generatedLegendPosition = 'right';
const partialLayoutImport = await importPayload(partialLayoutSession);
assert.equal(partialLayoutImport.status, 'ok');
assert.equal(state.form.legend, 'top');
assert.equal(state.adv.plot_title_position, 'bottom');
assert.equal(state.generatedLegendPosition.value, 'right');
assert.deepEqual(state.layoutPreferences, {
  circular: {
    single: { legend: 'left', plotTitlePosition: 'none' },
    multi: { legend: null, plotTitlePosition: null }
  },
  linear: { legend: 'top', plotTitlePosition: 'bottom' }
});

state.orthogroups.value = [{ id: 'stale-drawer-group', members: [] }];
state.rightDrawerTab.value = 'orthogroups';
state.showRightDrawer.value = false;
const groupFreeDrawerImport = await importPayload(partialLayoutSession);
assert.equal(groupFreeDrawerImport.status, 'ok');
assert.equal(state.showRightDrawer.value, false);
assert.equal(state.rightDrawerTab.value, 'features');

const modeBeforeRollbackSetup = state.mode.value;
state.mode.value = 'circular';
state.modeProfileStateManager.transition(
  state.adv,
  modeBeforeRollbackSetup,
  state.mode.value
);
state.form.multi_record_canvas = false;
state.form.legend = 'left';
state.adv.plot_title_position = 'none';
state.generatedLegendPosition.value = 'left';
state.sessionTitle.value = 'keep-after-normalization-error';
Object.assign(state.canvasPadding, { top: 7, right: 8, bottom: 9, left: 10 });
Object.assign(state.canvasPan, { x: 27, y: 28 });
state.zoom.value = 1.25;
state.skipCaptureBaseConfig.value = false;
state.skipPositionReapply.value = false;
state.suppressCircularMultiRecordDefaults.value = true;
state.linearReorderNotice.value = 'keep reorder notice';
state.showRightDrawer.value = true;
state.rightDrawerTab.value = 'legend';
state.showCanvasControls.value = true;
state.isPanning.value = true;
Object.assign(state.panStart, { x: 1, y: 2, panX: 3, panY: 4 });
const selectedAnnotation = { id: 'keep-annotation' };
state.selectedAnnotation.value = selectedAnnotation;
state.selectedSpecificPreset.value = 'keep-preset';
state.specificRulePresetLoading.value = true;
Object.keys(state.newSpecRule).forEach((key) => delete state.newSpecRule[key]);
Object.assign(state.newSpecRule, { qualifier: 'gene', value: 'keep-rule' });
Object.keys(state.newPriorityRule).forEach((key) => delete state.newPriorityRule[key]);
Object.assign(state.newPriorityRule, { qualifier: 'product', priority: 7 });
state.newColorFeat.value = 'repeat_region';
state.newColorVal.value = '#123456';
state.newFeatureToAdd.value = 'misc_feature';
state.newLegendCaption.value = 'Keep legend caption';
state.newLegendColor.value = '#654321';
state.fileLegendCaptions.value = new Set(['Keep file legend']);
state.featureSearch.value = 'keep feature search';
state.labelSearch.value = 'keep label search';
Object.assign(state.featureVisibilitySelectorCache, {
  'keep-feature': { featureType: 'CDS', recordId: 'KEEP' }
});
state.selectedFeatureIds.value = new Set(['keep-feature-a', 'keep-feature-b']);
state.selectedFeatureAnchorId.value = 'keep-feature-a';
state.featureSelectionStatus.value = 'Keep selection';
state.featureSelectionSuppressNextClick.value = true;
Object.assign(state.featureSelectionDrag, {
  active: true,
  committed: true,
  startX: 31,
  startY: 32,
  currentX: 33,
  currentY: 34
});
state.labelReflowLastError.value = { summary: 'keep reflow error' };
state.labelOverrideBuildWarning.value = 'keep override warning';
state.labelLayoutDirtyReason.value = 'keep dirty reason';
const clickedFeature = { id: 'keep-feature' };
const clickedPairwiseMatch = { id: 'keep-match' };
const clickedLabel = { key: 'keep-label' };
state.clickedFeature.value = clickedFeature;
state.clickedPairwiseMatch.value = clickedPairwiseMatch;
state.clickedLabel.value = clickedLabel;
state.featureExtractionPending.value = true;
state.featureExtractionError.value = { summary: 'keep extraction error' };
Object.assign(state.featureEditorStatus, {
  status: 'summary-ready',
  generationId: 'keep-generation',
  error: null,
  summaryCount: 3,
  detailsCacheSize: 2
});
state.matchSequenceRegistry.reset([{
  key: 'keep-sequence',
  recordId: 'KEEP.1',
  aliases: ['KEEP'],
  sequence: 'AACCGGTT',
  origin: 'linear-record',
  recordIndex: 0,
  sourceIndex: null
}]);
state.circularRecordList.value = [{
  selector: '#1',
  record_id: 'KEEP.1',
  record_length: 8
}];
Object.assign(state.circularRecordDiscovery, {
  status: 'ready',
  error: '',
  inputType: 'gb',
  primaryFile: state.files.c_gb,
  pairedFile: null
});
const diagramElement = { id: 'keep-diagram-element' };
const lengthBarElement = { id: 'keep-length-bar' };
const plotTitleElement = { id: 'keep-plot-title' };
state.diagramElements.value = [diagramElement];
state.diagramElementIds.value = [diagramElement.id];
state.diagramElementOriginalTransforms.value = new Map([
  [diagramElement, { x: 11, y: 12 }]
]);
state.legendDragging.value = true;
Object.assign(state.legendDragStart, { x: 19, y: 20 });
state.legendOriginalTransform.value = { x: 21, y: 22 };
state.legendInitialTransform.value = { x: 13, y: 14 };
state.diagramDragging.value = true;
Object.assign(state.diagramDragStart, { x: 23, y: 24 });
state.lengthBarElement.value = lengthBarElement;
state.lengthBarOriginalTransform.value = { x: 15, y: 16 };
state.plotTitleElement.value = plotTitleElement;
state.plotTitleDragging.value = true;
Object.assign(state.plotTitleDragStart, { x: 25, y: 26 });
state.plotTitleAutoTransform.value = { x: 17, y: 18 };
state.featureListScrollTop.value = 144;
state.semanticFileWatchersSuppressed.value = true;

const depthTrackUiCounts = { circular: 2 };
const depthTracks = [
  { label: 'Keep depth 1' },
  { label: 'Keep depth 2' }
];
const featureListScrollRef = { value: { scrollTop: 144 } };
const selectedPairwiseBlockOrthogroupId = { value: 'keep-orthogroup' };
const sessionImportRollbackState = createSessionImportRollbackState({
  depthTrackUiCounts,
  depthTracks,
  featureListScrollTop: state.featureListScrollTop,
  featureListScrollRef,
  selectedPairwiseBlockOrthogroupId
});

const jsonClone = (value) => JSON.parse(JSON.stringify(value));
const rollbackState = () => ({
  config: jsonClone(buildConfigData()),
  ui: jsonClone(buildUiStateData()),
  features: jsonClone(buildFeatureStateData()),
  editorState: jsonClone(buildEditorStateData()),
  orthogroupState: jsonClone(buildOrthogroupStateData()),
  runState: jsonClone(buildRunStateData()),
  results: jsonClone(state.results.value),
  mode: state.mode.value,
  title: state.sessionTitle.value,
  legend: state.form.legend,
  plotTitlePosition: state.adv.plot_title_position,
  generatedLegendPosition: state.generatedLegendPosition.value,
  documentEpoch: state.documentEpoch.value,
  semanticStyleRevision: state.semanticStyleRevision.value,
  semanticStyleFingerprint: state.semanticStyleFingerprint.value,
  validatedStyleFingerprintByResultKey: jsonClone(
    state.validatedStyleFingerprintByResultKey.value
  ),
  semanticFileWatchersSuppressed: state.semanticFileWatchersSuppressed.value,
  sessionImportRollbackInProgress: state.sessionImportRollbackInProgress.value,
  skipCaptureBaseConfig: state.skipCaptureBaseConfig.value,
  skipPositionReapply: state.skipPositionReapply.value,
  suppressCircularMultiRecordDefaults: state.suppressCircularMultiRecordDefaults.value,
  linearReorderNotice: state.linearReorderNotice.value,
  showRightDrawer: state.showRightDrawer.value,
  rightDrawerTab: state.rightDrawerTab.value,
  showCanvasControls: state.showCanvasControls.value,
  isPanning: state.isPanning.value,
  panStart: { ...state.panStart },
  selectedAnnotation: state.selectedAnnotation.value,
  selectedSpecificPreset: state.selectedSpecificPreset.value,
  specificRulePresetLoading: state.specificRulePresetLoading.value,
  newSpecRule: structuredClone(state.newSpecRule),
  newPriorityRule: structuredClone(state.newPriorityRule),
  newColorFeat: state.newColorFeat.value,
  newColorVal: state.newColorVal.value,
  newFeatureToAdd: state.newFeatureToAdd.value,
  newLegendCaption: state.newLegendCaption.value,
  newLegendColor: state.newLegendColor.value,
  fileLegendCaptions: [...state.fileLegendCaptions.value],
  featureSearch: state.featureSearch.value,
  labelSearch: state.labelSearch.value,
  featureVisibilitySelectorCache: structuredClone(state.featureVisibilitySelectorCache),
  selectedFeatureIds: [...state.selectedFeatureIds.value],
  selectedFeatureAnchorId: state.selectedFeatureAnchorId.value,
  featureSelectionStatus: state.featureSelectionStatus.value,
  featureSelectionSuppressNextClick: state.featureSelectionSuppressNextClick.value,
  featureSelectionDrag: structuredClone(state.featureSelectionDrag),
  labelReflowLastError: state.labelReflowLastError.value,
  labelOverrideBuildWarning: state.labelOverrideBuildWarning.value,
  labelLayoutDirtyReason: state.labelLayoutDirtyReason.value,
  clickedFeature: state.clickedFeature.value,
  clickedPairwiseMatch: state.clickedPairwiseMatch.value,
  clickedLabel: state.clickedLabel.value,
  featureExtractionPending: state.featureExtractionPending.value,
  featureExtractionError: state.featureExtractionError.value,
  featureEditorStatus: structuredClone(state.featureEditorStatus),
  matchSequenceSources: structuredClone(state.matchSequenceRegistry.values()),
  circularRecordList: structuredClone(state.circularRecordList.value),
  circularRecordDiscovery: { ...state.circularRecordDiscovery },
  diagramElements: [...state.diagramElements.value],
  diagramElementIds: [...state.diagramElementIds.value],
  diagramElementOriginalTransforms: [...state.diagramElementOriginalTransforms.value],
  legendDragging: state.legendDragging.value,
  legendDragStart: { ...state.legendDragStart },
  legendOriginalTransform: structuredClone(state.legendOriginalTransform.value),
  legendInitialTransform: structuredClone(state.legendInitialTransform.value),
  diagramDragging: state.diagramDragging.value,
  diagramDragStart: { ...state.diagramDragStart },
  lengthBarElement: state.lengthBarElement.value,
  lengthBarOriginalTransform: structuredClone(state.lengthBarOriginalTransform.value),
  plotTitleElement: state.plotTitleElement.value,
  plotTitleDragging: state.plotTitleDragging.value,
  plotTitleDragStart: { ...state.plotTitleDragStart },
  plotTitleAutoTransform: structuredClone(state.plotTitleAutoTransform.value),
  circularDepthTrackUiCount: depthTrackUiCounts.circular,
  depthTracks: structuredClone(depthTracks),
  featureListScrollTop: state.featureListScrollTop.value,
  featureListElementScrollTop: featureListScrollRef.value.scrollTop,
  selectedPairwiseBlockOrthogroupId: selectedPairwiseBlockOrthogroupId.value
});

const recordDiscoveryDerivedState = {
  multiRecordPositions: [{ selector: 'KEEP.1', row: 2 }],
  circularRecords: ['KEEP.1'],
  linearSelectorCache: ['KEEP.1']
};
const expectedRecordDiscoveryDerivedState = structuredClone(
  recordDiscoveryDerivedState
);
state.sessionImportRollbackInProgress.value = true;
let guardedRecordDiscoveryResult;
let recordDiscoveryGeneration = 0;
try {
  guardedRecordDiscoveryResult = await runRecordDiscoveryWatcher({
    rollbackInProgress: state.sessionImportRollbackInProgress,
    semanticWatchersSuppressed: state.semanticFileWatchersSuppressed,
    refresh: async ({ suppress }) => {
      recordDiscoveryGeneration += 1;
      if (suppress) return;
      recordDiscoveryDerivedState.multiRecordPositions = [];
      recordDiscoveryDerivedState.circularRecords = [];
      recordDiscoveryDerivedState.linearSelectorCache = [];
      throw new Error('injected restored-file reparse failure');
    }
  });
} finally {
  state.sessionImportRollbackInProgress.value = false;
}
assert.equal(guardedRecordDiscoveryResult, false);
assert.equal(recordDiscoveryGeneration, 1);
assert.deepEqual(
  recordDiscoveryDerivedState,
  expectedRecordDiscoveryDerivedState
);
const stateBeforeFailedImport = rollbackState();

const malformedCompositionSvg = {
  getAttribute: (name) => {
    if (name === COMPOSITION_SCHEMA_ATTRIBUTE) return '1';
    if (name === COMPOSITION_METADATA_ATTRIBUTE) return '{broken';
    return null;
  }
};
const malformedCompositionContainer = {
  querySelector: (selector) => selector === 'svg' ? malformedCompositionSvg : null
};
const legendCanvasActions = createLegendCanvasActions({ state });
alerts.length = 0;
const originalConsoleError = console.error;
console.error = () => {};
let failedImport;
const failedImportEvent = {
  target: {
    files: [new Blob([JSON.stringify(divergentSession)], { type: 'application/json' })],
    value: 'selected'
  }
};
try {
  failedImport = await importSession(
    failedImportEvent,
    {
      rollbackState: sessionImportRollbackState,
      afterLoad: async () => {
        state.documentEpoch.value += 11;
        state.semanticStyleRevision.value += 7;
        state.semanticStyleFingerprint.value = 'sf1_interrupted-import';
        state.validatedStyleFingerprintByResultKey.value = Object.freeze({
          'interrupted-result': 'sf1_interrupted-import'
        });
        depthTrackUiCounts.circular = 5;
        state.featureListScrollTop.value = 0;
        featureListScrollRef.value.scrollTop = 0;
        selectedPairwiseBlockOrthogroupId.value = '';
        depthTracks.splice(
          0,
          depthTracks.length,
          ...Array.from({ length: 5 }, (_, index) => ({
            label: `Imported depth ${index + 1}`
          }))
        );
        const originalSvgContainer = state.svgContainer.value;
        state.svgContainer.value = malformedCompositionContainer;
        try {
          legendCanvasActions.captureBaseConfig();
        } finally {
          state.svgContainer.value = originalSvgContainer;
        }
      }
    }
  );
} finally {
  console.error = originalConsoleError;
}

assert.equal(failedImport.status, 'error');
assert.match(failedImport.error?.message || '', /composition metadata is not valid JSON/);
assert.deepEqual(rollbackState(), stateBeforeFailedImport);
assert.equal(alerts.length, 1);
assert.match(alerts[0], /^Failed to load session: .*composition metadata is not valid JSON/);
assert.equal(failedImportEvent.target.value, '');

const epochBeforeCurrentImport = state.documentEpoch.value;
const currentImport = await importPayload(currentHashBearingSession);
assert.equal(currentImport.status, 'ok');
assert.equal(state.documentEpoch.value, epochBeforeCurrentImport + 1);
assert.equal(state.semanticStyleRevision.value, 0);
const importedStyleSnapshot = captureStyleSnapshot(state);
assert.equal(state.semanticStyleFingerprint.value, importedStyleSnapshot.fingerprint);
const importedResultKeys = state.featureCatalog.value.items.map(catalogResultKey);
assert.equal(importedResultKeys.length, state.results.value.length);
assert.deepEqual(
  state.validatedStyleFingerprintByResultKey.value,
  Object.fromEntries(importedResultKeys.map((resultKey) => [
    resultKey,
    importedStyleSnapshot.fingerprint
  ]))
);

const exportPairFeature = state.featureCatalog.value.items
  .flatMap((item) => item.biologicalFeatures)
  .find((feature) => feature.instanceHash);
assert.ok(exportPairFeature);
const exportPairHash = exportPairFeature.instanceHash;
exportPairFeature.instanceHash = exportPairHash === 'fi1_aaaaaaaaaaaaaaaaaaaaaaaaaa'
  ? 'fi1_bbbbbbbbbbbbbbbbbbbbbbbbbb'
  : 'fi1_aaaaaaaaaaaaaaaaaaaaaaaaaa';
const consoleWarnBeforeMismatchedPairExport = console.warn;
console.warn = () => {};
try {
  await assert.rejects(
    exportSession('mismatched-feature-pair'),
    /Generate again before using Save Session/
  );
} finally {
  console.warn = consoleWarnBeforeMismatchedPairExport;
  exportPairFeature.instanceHash = exportPairHash;
}

const originalCreateElement = document.createElement;
document.createElement = () => ({
  addEventListener: () => {},
  click: () => {},
  parentNode: null
});
try {
  const saved = await exportSession('current-hash-bearing-roundtrip');
  assert.equal(saved.status, 'saved');
} finally {
  if (originalCreateElement === undefined) delete document.createElement;
  else document.createElement = originalCreateElement;
}
