import { installSessionImportWorker } from './helpers/session-import-node.mjs';
import assert from 'node:assert/strict';
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

installSessionImportWorker();

const {
  buildConfigData,
  buildEditorStateData,
  buildFeatureStateData,
  buildOrthogroupStateData,
  buildRunStateData,
  buildUiStateData,
  importSession,
  restoreCurrentWriterActiveConfig,
  serializeActiveRenderFiles
} = await import('../../gbdraw/web/js/services/config.js');
const {
  CURRENT_WRITER_ACTIVE_CONFIG_DOMAINS,
  validateCurrentWriterActiveConfig
} = await import('../../gbdraw/web/js/services/session-active-config-contract.js');
const { state } = await import('../../gbdraw/web/js/state.js');
const {
  buildCanonicalRenderRequest,
  projectCanonicalSessionRequest
} = await import('../../gbdraw/web/js/services/session-request.js');
state.selectedOrthogroupAlignmentFeature.value = 'legacy-selection';
assert.equal(
  Object.hasOwn(buildOrthogroupStateData(state.activeDrawing()), 'selectedOrthogroupAlignmentFeature'),
  false
);
state.selectedOrthogroupAlignmentFeature.value = '';
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

const compactActiveIntentSnapshot = () => ({
  modeAndInput: {
    mode: state.mode.value,
    inputType: state.mode.value === 'circular'
      ? state.cInputType.value
      : state.lInputType.value
  },
  palette: {
    selected: state.activeDrawing().selectedPalette.value,
    currentColors: {
      CDS: state.activeDrawing().currentColors.value.CDS,
      tRNA: state.activeDrawing().currentColors.value.tRNA
    },
    appliedName: state.appliedPaletteName.value,
    appliedColors: {
      CDS: state.appliedPaletteColors.value.CDS,
      tRNA: state.appliedPaletteColors.value.tRNA
    },
    pendingName: state.activeDrawing().pendingPaletteName.value,
    pendingColors: {
      CDS: state.activeDrawing().pendingPaletteColors.value.CDS,
      tRNA: state.activeDrawing().pendingPaletteColors.value.tRNA
    },
    instantPreview: state.paletteInstantPreviewEnabled.value
  },
  specificRules: structuredClone(state.activeDrawing().manualSpecificRules),
  qualifierPriorityRules: structuredClone(state.activeDrawing().manualPriorityRules),
  filters: {
    mode: state.activeDrawing().filterMode.value,
    whitelist: structuredClone(state.activeDrawing().manualWhitelist),
    blacklistText: state.activeDrawing().manualBlacklist.value
  },
  form: {
    plot_title: state.activeDrawing().form.plot_title,
    labels_mode: state.activeDrawing().form.labels_mode,
    show_scale: state.activeDrawing().form.show_scale,
    legend: state.activeDrawing().form.legend
  },
  adv: {
    axis_stroke_width: state.activeDrawing().adv.axis_stroke_width,
    label_font_size: state.activeDrawing().adv.label_font_size,
    feature_width_circular: state.activeDrawing().adv.feature_width_circular
  },
  annotationSets: state.activeDrawing().annotationSets.map((set) => ({
    id: set.id,
    annotations: set.annotations.map((annotation) => ({
      id: annotation.id,
      mark: annotation.mark
    }))
  })),
  trackSlots: {
    enabled: state.activeDrawing().adv.circular_track_slots_enabled,
    axisIndex: state.activeDrawing().adv.circular_track_slots_axis_index,
    slots: state.activeDrawing().adv.circular_track_slots.map((slot) => ({
      id: slot.id,
      renderer: slot.renderer,
      enabled: slot.enabled
    }))
  },
  layout: {
    preferences: structuredClone(state.activeDrawing().layoutPreferences),
    linearRecordLayout: {
      enabled: state.activeDrawing().linearRecordLayoutEnabled.value,
      recordGap: state.activeDrawing().linearRecordGap.value,
      rows: state.activeDrawing().linearRecordLayoutEnabled.value
        ? structuredClone(state.activeDrawing().linearRecordRows)
        : []
    }
  },
  editorOverrides: {
    fills: structuredClone(state.activeDrawing().featureColorOverrides),
    strokes: structuredClone(state.activeDrawing().featureStrokeOverrides),
    featureOverrides: structuredClone(state.activeDrawing().featureOverrides),
    legendColors: structuredClone(state.activeDrawing().legendColorOverrides),
    legendStrokes: structuredClone(state.activeDrawing().legendStrokeOverrides)
  }
});

const activeIntentDomainMismatches = (expected, actual) => {
  const mismatches = [];
  const visit = (expectedValue, actualValue, path, depth) => {
    if (
      depth < 2 &&
      expectedValue &&
      typeof expectedValue === 'object' &&
      !Array.isArray(expectedValue)
    ) {
      Object.entries(expectedValue).forEach(([key, value]) => {
        visit(value, actualValue?.[key], [...path, key], depth + 1);
      });
      return;
    }
    try {
      assert.deepEqual(actualValue, expectedValue);
    } catch {
      mismatches.push(path.join('.'));
    }
  };
  visit(expected, actual, [], 0);
  return mismatches;
};

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
  params: { lane_direction: 'inside' }
};
const disabledDraft = {
  id: 'disabled-draft',
  renderer: 'depth',
  enabled: false,
  side: 'outside',
  width: 27,
  radius: 1.2,
  inner_gap_px: 4,
  outer_gap_px: 5,
  z: 3,
  params: { track_index: 99, nested: { keep: true } }
};
const projectedConfig = {
  form: { track_type: 'tuckin' },
  adv: {
    nt: 'GC',
    ruler_label_font_size: 13,
    circular_track_slots_enabled: true,
    circular_track_slots_schema_version: 4,
    circular_track_slots: [canonicalFeature],
    circular_track_slots_axis_index: 1,
    linear_track_slots_enabled: false,
    linear_track_slots_schema_version: 2,
    linear_track_slots: [],
    linear_track_slots_axis_index: null
  }
};
const storedConfig = {
  form: { track_type: 'tuckin' },
  unmanagedConfigOverrides: {
    'objects.gc_content.percent_background_opacity': 0.42
  },
  adv: {
    ...projectedConfig.adv,
    linear_accession_visibility: 'auto',
    linear_length_visibility: 'auto',
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

assert.deepEqual(CURRENT_WRITER_ACTIVE_CONFIG_DOMAINS, [
  'form',
  'adv',
  'losat',
  'cliOptions',
  'colors',
  'palette',
  'paletteInstantPreviewEnabled',
  'rules',
  'qualifierPriorityRules',
  'filterMode',
  'whitelist',
  'blacklistText',
  'losatProgram',
  'circularConservation',
  'annotationSets',
  'recordDisplayDrafts',
  'featurePlacementOverrides',
  'modeProfiles',
  'unmanagedConfigOverrides',
  'linearRecordLayout',
  'linearComparisonPlan',
  'importedComparisonResolution',
  'webEdits'
]);
assert.doesNotThrow(() => validateCurrentWriterActiveConfig({
  mode: 'circular',
  storedConfig
}));
assert.throws(
  () => validateCurrentWriterActiveConfig({
    mode: 'linear',
    storedConfig: {
      ...storedConfig,
      importedComparisonResolution: { action: 'AUTOMATIC' }
    }
  }),
  /importedComparisonResolution\.action is invalid/
);
Object.assign(state.activeDrawing().importedComparisonIntent, {
  disposition: 'PRESERVED_READ_ONLY',
  action: 'INHERIT',
  message: 'Preserved for the current Session.',
  hasCommittedComparison: true
});
assert.deepEqual(buildConfigData(state.activeDrawing()).importedComparisonResolution, { action: 'INHERIT' });
Object.assign(state.activeDrawing().importedComparisonIntent, {
  disposition: 'EDITABLE',
  action: null,
  message: '',
  hasCommittedComparison: false
});
const restored = restoreCurrentWriterActiveConfig({
  mode: 'circular',
  projectedConfig,
  storedConfig
});
assert.deepEqual(restored.adv.circular_track_slots, storedConfig.adv.circular_track_slots);
assert.deepEqual(restored.adv.linear_track_slots, storedConfig.adv.linear_track_slots);
assert.equal(restored.adv.circular_track_slots_axis_index, 2);
assert.equal(restored.adv.feature_width_circular, 19);
assert.equal(restored.adv.depth_width_circular, 23);
assert.deepEqual(restored.unmanagedConfigOverrides, {
  'objects.gc_content.percent_background_opacity': 0.42
});

const storedBeforeIndependentRulerFont = structuredClone(storedConfig);
delete storedBeforeIndependentRulerFont.adv.ruler_label_font_size;
assert.equal(
  restoreCurrentWriterActiveConfig({
    mode: 'linear',
    projectedConfig,
    storedConfig: storedBeforeIndependentRulerFont
  }).adv.ruler_label_font_size,
  13
);

const galleryCompatibilityConfig = structuredClone(storedConfig);
galleryCompatibilityConfig.colors = { CDS: '#123456' };
galleryCompatibilityConfig.colorsAreOverrides = true;
galleryCompatibilityConfig.adv.losatProgram = 'blastp';
assert.doesNotThrow(() => validateCurrentWriterActiveConfig({
  mode: 'circular',
  storedConfig: galleryCompatibilityConfig
}));
const restoredGalleryCompatibilityConfig = restoreCurrentWriterActiveConfig({
  mode: 'circular',
  projectedConfig,
  storedConfig: galleryCompatibilityConfig
});
assert.equal(
  Object.prototype.hasOwnProperty.call(restoredGalleryCompatibilityConfig, 'colorsAreOverrides'),
  false
);
assert.equal(
  Object.prototype.hasOwnProperty.call(restoredGalleryCompatibilityConfig.adv, 'losatProgram'),
  false
);

const mismatched = structuredClone(storedConfig);
mismatched.adv.circular_track_slots[1].width = 12;
assert.doesNotThrow(() => validateCurrentWriterActiveConfig({
  mode: 'circular',
  storedConfig: mismatched
}));
assert.equal(
  restoreCurrentWriterActiveConfig({
    mode: 'circular',
    projectedConfig,
    storedConfig: mismatched
  }).adv.circular_track_slots[1].width,
  12
);

const inactive = structuredClone(mismatched);
inactive.adv.circular_track_slots_enabled = false;
assert.doesNotThrow(() => validateCurrentWriterActiveConfig({
  mode: 'circular',
  storedConfig: inactive
}));

const structurallyInvalid = structuredClone(storedConfig);
structurallyInvalid.adv.circular_track_slots[1].side = 'left';
assert.throws(
  () => validateCurrentWriterActiveConfig({
    mode: 'circular',
    storedConfig: structurallyInvalid
  }),
  /unsupported side/
);

const divergentSession = JSON.parse(await readFile(
  'gbdraw/web/gallery/sessions/BGC0000708-BGC0000713.gbdraw-session.json',
  'utf8'
));
// The Linear drawing's draft is its Session 46 slice.
const divergentSlice = divergentSession.modes.linear;
divergentSlice.ui ||= {};
// This test isolates comparison draft behavior from the unreleased schema-1 Gallery plan.
divergentSession.renderRequest.layout.similarityAlignment = null;
const committedComparisonCount = divergentSession.renderRequest.comparisons.length;
assert.ok(committedComparisonCount > 0);
divergentSlice.config.linearComparisonPlan = {
  mode: 'none',
  defaultSource: 'losat',
  edges: []
};
divergentSlice.config.losatProgram = 'tblastx';
Object.assign(divergentSlice.config.adv, {
  comparison_height: 73,
  min_bitscore: 123,
  evalue: '1e-37',
  identity: 88,
  alignment_length: 456,
  pairwise_match_style: 'ribbon'
});
Object.assign(divergentSlice.config.losat.blastp, {
  mode: 'collinear',
  maxHits: 17,
  candidateLimit: 11,
  orthogroupMemberMaxHits: 19,
  collinearMinAnchors: 4,
  collinearMaxUnitGap: 9,
  collinearUnitMode: 'locus',
  collinearAnchorMode: 'all',
  collinearMergeOrientation: 'strand',
  collinearSearchScope: 'adjacent',
  collinearColorMode: 'orientation_identity'
});
const dormantComparisonDraft = structuredClone({
  losatProgram: divergentSlice.config.losatProgram,
  adv: {
    comparison_height: divergentSlice.config.adv.comparison_height,
    min_bitscore: divergentSlice.config.adv.min_bitscore,
    evalue: divergentSlice.config.adv.evalue,
    identity: divergentSlice.config.adv.identity,
    alignment_length: divergentSlice.config.adv.alignment_length,
    pairwise_match_style: divergentSlice.config.adv.pairwise_match_style
  },
  blastp: {
    mode: divergentSlice.config.losat.blastp.mode,
    maxHits: divergentSlice.config.losat.blastp.maxHits,
    candidateLimit: divergentSlice.config.losat.blastp.candidateLimit,
    orthogroupMemberMaxHits:
      divergentSlice.config.losat.blastp.orthogroupMemberMaxHits,
    collinearMinAnchors:
      divergentSlice.config.losat.blastp.collinearMinAnchors,
    collinearMaxUnitGap:
      divergentSlice.config.losat.blastp.collinearMaxUnitGap,
    collinearUnitMode:
      divergentSlice.config.losat.blastp.collinearUnitMode,
    collinearAnchorMode:
      divergentSlice.config.losat.blastp.collinearAnchorMode,
    collinearMergeOrientation:
      divergentSlice.config.losat.blastp.collinearMergeOrientation,
    collinearSearchScope:
      divergentSlice.config.losat.blastp.collinearSearchScope,
    collinearColorMode:
      divergentSlice.config.losat.blastp.collinearColorMode
  }
});
divergentSession.renderRequest.diagramOptions.output.legend = 'right';
divergentSession.ui.generatedLegendPosition = 'right';
// The Linear slot of the layout (the slice holds its own mode's slot).
divergentSlice.ui.layoutPreferences = { legend: 'bottom', plotTitlePosition: 'bottom' };
divergentSession.ui.paletteInstantPreviewEnabled = false;
divergentSession.ui.appliedPaletteName = 'orchid';
divergentSession.ui.appliedPaletteColors = { CDS: '#123456' };
divergentSlice.ui.pendingPaletteName = 'mint';
divergentSlice.ui.pendingPaletteColors = { CDS: '#abcdef' };
const importEvent = {
  target: {
    files: [new Blob([JSON.stringify(divergentSession)], { type: 'application/json' })],
    value: 'selected'
  }
};

const imported = await importSession(importEvent);

assert.equal(imported.status, 'ok', imported.error?.message);
assert.equal(imported.comparisonDisposition, 'EDITABLE');
assert.deepEqual(state.activeDrawing().importedComparisonIntent, {
  disposition: 'EDITABLE',
  action: null,
  message: 'The saved comparison is represented by the current controls.',
  hasCommittedComparison: true
});
assert.equal(state.activeDrawing().linearComparisonPlan.mode, 'none');
assert.deepEqual(state.activeDrawing().linearComparisonPlan.edges, []);
assert.deepEqual({
  losatProgram: state.activeDrawing().losatProgram.value,
  adv: {
    comparison_height: state.activeDrawing().adv.comparison_height,
    min_bitscore: state.activeDrawing().adv.min_bitscore,
    evalue: state.activeDrawing().adv.evalue,
    identity: state.activeDrawing().adv.identity,
    alignment_length: state.activeDrawing().adv.alignment_length,
    pairwise_match_style: state.activeDrawing().adv.pairwise_match_style
  },
  blastp: {
    mode: state.activeDrawing().losat.blastp.mode,
    maxHits: state.activeDrawing().losat.blastp.maxHits,
    candidateLimit: state.activeDrawing().losat.blastp.candidateLimit,
    orthogroupMemberMaxHits: state.activeDrawing().losat.blastp.orthogroupMemberMaxHits,
    collinearMinAnchors: state.activeDrawing().losat.blastp.collinearMinAnchors,
    collinearMaxUnitGap: state.activeDrawing().losat.blastp.collinearMaxUnitGap,
    collinearUnitMode: state.activeDrawing().losat.blastp.collinearUnitMode,
    collinearAnchorMode: state.activeDrawing().losat.blastp.collinearAnchorMode,
    collinearMergeOrientation: state.activeDrawing().losat.blastp.collinearMergeOrientation,
    collinearSearchScope: state.activeDrawing().losat.blastp.collinearSearchScope,
    collinearColorMode: state.activeDrawing().losat.blastp.collinearColorMode
  }
}, dormantComparisonDraft);
assert.deepEqual(
  JSON.parse(JSON.stringify(buildConfigData(state.activeDrawing()).losat.blastp)),
  JSON.parse(JSON.stringify(state.activeDrawing().losat.blastp)),
  'the current Session writer must retain every editable LOSATP value'
);
assert.equal(
  state.activeDrawing().form.legend,
  'bottom',
  'the active editor preference must override the last generated request position'
);
assert.equal(state.activeDrawing().adv.plot_title_position, 'bottom');
assert.equal(state.generatedLegendPosition.value, 'right');
assert.equal(state.appliedPaletteName.value, 'orchid');
assert.equal(state.appliedPaletteColors.value.CDS, '#123456');
assert.equal(state.activeDrawing().pendingPaletteName.value, 'mint');
assert.equal(state.activeDrawing().pendingPaletteColors.value.CDS, '#abcdef');
assert.deepEqual(state.activeDrawing().layoutPreferences.linear, {
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
  delete payload.modes.linear.ui.layoutPreferences;
  legacyLayoutFields.forEach((field) => delete payload.ui[field]);
  return payload;
};
const importPayload = async (payload) => importSession({
  target: {
    files: [new Blob([JSON.stringify(payload)], { type: 'application/json' })],
    value: 'selected'
  }
});

const decisionRequiredSession = structuredClone(divergentSession);
decisionRequiredSession.renderRequest.comparisons = [{
  kind: 'nucleotideBlast',
  resourceId: 'missing-comparison-resource',
  queryRecordIndex: 0,
  subjectRecordIndex: 1
}];
decisionRequiredSession.losatCache = { entries: [] };
decisionRequiredSession.losatDerivedCache = { entries: [] };
const decisionRequiredImport = await importPayload(decisionRequiredSession);
assert.equal(
  decisionRequiredImport.status,
  'ok',
  decisionRequiredImport.error?.message
);
assert.equal(decisionRequiredImport.comparisonDisposition, 'DECISION_REQUIRED');
assert.deepEqual(state.activeDrawing().importedComparisonIntent, {
  disposition: 'DECISION_REQUIRED',
  action: null,
  message: 'The saved comparison is missing a required resource.',
  hasCommittedComparison: true
});
assert.equal(state.results.value.length, divergentSession.results.length);
assert.deepEqual(state.files.linearCanonicalComparisons, []);

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
  savedPlanSession.modes.linear.config.linearComparisonPlan = structuredClone(savedPlan);
  const savedPlanImport = await importPayload(savedPlanSession);
  assert.equal(savedPlanImport.status, 'ok');
  assert.deepEqual(state.activeDrawing().linearComparisonPlan, {
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
assert.equal(state.activeDrawing().form.legend, 'right');
assert.equal(state.generatedLegendPosition.value, 'right');
assert.equal(state.activeDrawing().layoutPreferences.linear.legend, 'right');

// Session 46 reads the layout from the slice only: the `ui` fields of the
// layout before Session 44 are not read, so the committed request's layout
// fills the slot.
const legacyLayoutSession = withoutStoredLayoutPreferences();
legacyLayoutSession.ui.legend = 'bottom';
legacyLayoutSession.ui.linearPlotTitlePosition = 'top';
legacyLayoutSession.ui.generatedLegendPosition = 'right';
const legacyLayoutImport = await importPayload(legacyLayoutSession);
assert.equal(legacyLayoutImport.status, 'ok');
assert.equal(state.activeDrawing().form.legend, 'right');
assert.equal(state.generatedLegendPosition.value, 'right');
assert.equal(state.activeDrawing().layoutPreferences.linear.legend, 'right');

const partialLayoutSession = withoutStoredLayoutPreferences();
partialLayoutSession.modes.linear.ui.layoutPreferences = { legend: 'top' };
partialLayoutSession.ui.generatedLegendPosition = 'right';
const partialLayoutImport = await importPayload(partialLayoutSession);
assert.equal(partialLayoutImport.status, 'ok');
assert.equal(state.activeDrawing().form.legend, 'top');
assert.equal(state.activeDrawing().adv.plot_title_position, 'bottom');
assert.equal(state.generatedLegendPosition.value, 'right');
assert.deepEqual(state.activeDrawing().layoutPreferences, {
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

const activeIntentSession = JSON.parse(await readFile(
  'gbdraw/web/gallery/sessions/HmmtDNA_basic_circular.gbdraw-session.json',
  'utf8'
));
const activeFeature = activeIntentSession.editorState.featureCatalog.items[0].features[0];
const activeFeatureId = activeFeature.svgId;
// The Circular drawing's draft is its Session 46 slice; per-feature edits are
// keyed by the feature's source identity (R2).
const activeSlice = activeIntentSession.modes.circular;
activeSlice.ui ||= {};
const activeFeatureIdentity = JSON.stringify([activeFeature.recordKey, activeFeature.biologicalFeatureId]);
Object.assign(activeSlice.config.form, {
  plot_title: 'Saved active draft',
  labels_mode: 'both',
  show_scale: false
});
Object.assign(activeSlice.config.adv, {
  axis_stroke_width: 7,
  label_font_size: 31,
  feature_width_circular: 19,
  circular_track_slots_enabled: false,
  circular_track_slots_schema_version: 4,
  circular_track_slots_axis_index: 0,
  circular_track_slots: [{
    ...disabledDraft,
    width: 27,
    radius: 1.2,
    inner_gap_px: 4,
    outer_gap_px: 5
  }]
});
Object.assign(activeSlice.config, {
  colors: { CDS: '#123456', tRNA: '#abcdef' },
  palette: 'orange',
  rules: [{
    feat: 'CDS',
    qual: 'gene',
    val: '^ND1$',
    color: '#f97316',
    cap: 'Saved rule',
    fromFile: false
  }],
  qualifierPriorityRules: [{ feat: 'CDS', order: 'gene,product' }],
  filterMode: 'Whitelist',
  whitelist: [{ feat: 'CDS', qual: 'gene', key: 'ND1' }],
  blacklistText: 'hypothetical, draft-only',
  annotationSets: [{
    id: 'saved-annotations',
    annotations: [{
      id: 'saved-window',
      target: { kind: 'coordinateSpan', start: 10, end: 30 },
      label: 'Saved window',
      mark: 'band'
    }]
  }],
  linearRecordLayout: {
    enabled: false,
    recordGap: 41,
    rows: []
  },
  linearComparisonPlan: { mode: 'none', defaultSource: 'losat', edges: [] }
});
activeIntentSession.ui.paletteInstantPreviewEnabled = false;
activeIntentSession.ui.appliedPaletteName = 'default';
activeIntentSession.ui.appliedPaletteColors = { CDS: '#e8b441', tRNA: '#71ee7d' };
activeSlice.ui.pendingPaletteName = 'orange';
activeSlice.ui.pendingPaletteColors = { CDS: '#123456', tRNA: '#abcdef' };
activeSlice.ui.selectedFeatureRecordIdx = 0;
activeSlice.ui.layoutPreferences = {
  single: { legend: 'left', plotTitlePosition: 'top' },
  multi: { legend: 'right', plotTitlePosition: 'bottom' }
};
activeSlice.features = {
  featureColorOverrides: {
    [activeFeatureId]: { color: '#334455', caption: 'Saved feature fill' }
  },
  featureVisibilityManualRules: [],
  featureOverrides: {
    [activeFeatureIdentity]: {
      recordKey: activeFeature.recordKey,
      biologicalFeatureId: activeFeature.biologicalFeatureId,
      featureVisibility: 'off',
      labelVisibility: 'off',
      labelText: 'Saved feature label',
      labelSourceText: 'Original label'
    }
  },
  labelOverrideRows: [],
  labelTextBulkOverrides: {}
};
activeSlice.editorState = {
  legend: {
    entries: [],
    deletedEntries: [],
    colorOverrides: { 'Saved feature fill': '#334455' },
    strokeOverrides: {
      'Saved feature fill': { strokeColor: '#112233', strokeWidth: 2 }
    },
    addedCaptions: []
  },
  featureStrokes: {
    overrides: {
      [activeFeatureId]: { strokeColor: '#654321', strokeWidth: 3 }
    }
  }
};
activeIntentSession.editorState = {
  ...activeIntentSession.editorState,
  legend: { originalOrder: [], originalColors: {} },
  originalSvgStroke: { color: '#000000', width: 1 }
};

const expectedActiveIntent = {
  modeAndInput: { mode: 'circular', inputType: 'gb' },
  palette: {
    selected: 'orange',
    currentColors: { CDS: '#123456', tRNA: '#abcdef' },
    appliedName: 'default',
    appliedColors: { CDS: '#e8b441', tRNA: '#71ee7d' },
    pendingName: 'orange',
    pendingColors: { CDS: '#123456', tRNA: '#abcdef' },
    instantPreview: false
  },
  specificRules: structuredClone(activeSlice.config.rules),
  qualifierPriorityRules: structuredClone(activeSlice.config.qualifierPriorityRules),
  filters: {
    mode: 'Whitelist',
    whitelist: structuredClone(activeSlice.config.whitelist),
    blacklistText: 'hypothetical, draft-only'
  },
  form: {
    plot_title: 'Saved active draft',
    labels_mode: 'both',
    show_scale: false,
    legend: 'left'
  },
  adv: {
    axis_stroke_width: 7,
    label_font_size: 31,
    feature_width_circular: 19
  },
  annotationSets: [{
    id: 'saved-annotations',
    annotations: [{ id: 'saved-window', mark: 'band' }]
  }],
  trackSlots: {
    enabled: false,
    axisIndex: 0,
    slots: [{ id: 'disabled-draft', renderer: 'depth', enabled: false }]
  },
  layout: {
    // The Circular drawing's slot; its Linear slot keeps the defaults.
    preferences: {
      circular: structuredClone(activeSlice.ui.layoutPreferences),
      linear: { legend: 'bottom', plotTitlePosition: 'bottom' }
    },
    linearRecordLayout: structuredClone(activeSlice.config.linearRecordLayout)
  },
  editorOverrides: {
    fills: structuredClone(activeSlice.features.featureColorOverrides),
    strokes: structuredClone(activeSlice.editorState.featureStrokes.overrides),
    featureOverrides: {
      [activeFeatureIdentity]: {
        recordKey: activeFeature.recordKey,
        biologicalFeatureId: activeFeature.biologicalFeatureId,
        featureVisibility: 'off',
        labelVisibility: 'off',
        labelText: 'Saved feature label',
        labelSourceText: 'Original label'
      }
    },
    legendColors: structuredClone(activeSlice.editorState.legend.colorOverrides),
    legendStrokes: structuredClone(activeSlice.editorState.legend.strokeOverrides)
  }
};

const activeIntentImport = await importPayload(activeIntentSession);
assert.equal(
  activeIntentImport.status,
  'ok',
  activeIntentImport.error?.message || 'active-intent session import failed'
);
const afterSessionLoadIntent = compactActiveIntentSnapshot();
const immediatelyBeforeGenerateIntent = compactActiveIntentSnapshot();
const activeFiles = await serializeActiveRenderFiles(state.mode.value, state, state.activeDrawing());
const firstGeneratedCanonical = buildCanonicalRenderRequest({ state, drawing: state.activeDrawing(), filesData: activeFiles });
const firstGeneratedProjection = projectCanonicalSessionRequest(firstGeneratedCanonical);

const afterLoadMismatches = activeIntentDomainMismatches(
  expectedActiveIntent,
  afterSessionLoadIntent
);
const beforeGenerateMismatches = activeIntentDomainMismatches(
  expectedActiveIntent,
  immediatelyBeforeGenerateIntent
);
assert.deepEqual(
  afterLoadMismatches,
  [],
  `active-intent domains lost after session load: ${afterLoadMismatches.join(', ')}; ` +
    `expected=${JSON.stringify(expectedActiveIntent.layout)} ` +
    `actual=${JSON.stringify(afterSessionLoadIntent.layout)}`
);
assert.deepEqual(
  beforeGenerateMismatches,
  [],
  `active-intent domains lost before first Generate: ${beforeGenerateMismatches.join(', ')}`
);
assert.equal(firstGeneratedProjection.mode, 'circular');
assert.equal(firstGeneratedProjection.inputType, 'gb');
assert.equal(firstGeneratedProjection.config.palette, 'orange');
assert.equal(firstGeneratedProjection.config.colors.CDS, '#123456');
assert.equal(firstGeneratedProjection.config.colors.tRNA, '#abcdef');
assert.deepEqual(
  firstGeneratedProjection.config.rules,
  activeSlice.config.rules.map(({ fromFile: _fromFile, ...rule }) => rule)
);
assert.deepEqual(
  firstGeneratedProjection.config.qualifierPriorityRules,
  activeSlice.config.qualifierPriorityRules
);
assert.equal(firstGeneratedProjection.config.filterMode, 'Whitelist');
assert.deepEqual(firstGeneratedProjection.config.whitelist, activeSlice.config.whitelist);
assert.equal(firstGeneratedProjection.config.form.plot_title, 'Saved active draft');
assert.equal(firstGeneratedProjection.config.form.labels_mode, 'both');
assert.equal(firstGeneratedProjection.config.form.show_scale, false);
assert.equal(firstGeneratedProjection.config.adv.axis_stroke_width, 7);
assert.equal(firstGeneratedProjection.config.adv.label_font_size, 31);
assert.deepEqual(
  firstGeneratedCanonical.renderRequest.diagramOptions.tracks.circularTrackSlots.find(
    (slot) => slot.renderer === 'features'
  )?.width,
  { value: 19, unit: 'px' }
);
assert.deepEqual(
  firstGeneratedProjection.config.annotationSets.map((set) => set.id),
  ['saved-annotations']
);

const stateBeforeInvalidActiveConfig = {
  activeIntent: compactActiveIntentSnapshot(),
  results: state.results.value,
  featureCatalog: state.featureCatalog.value,
  primaryFile: state.files.c_gb,
  svgContainer: state.svgContainer.value,
  sessionTitle: state.sessionTitle.value
};
const invalidActiveConfigSession = structuredClone(activeIntentSession);
invalidActiveConfigSession.modes.circular.config.form.linear_track_layout = 'spreadout';
const invalidActiveConfigEvent = {
  target: {
    files: [new Blob([JSON.stringify(invalidActiveConfigSession)], { type: 'application/json' })],
    value: 'selected'
  }
};
alerts.length = 0;
const consoleErrorBeforeInvalidActiveConfig = console.error;
console.error = () => {};
let invalidActiveConfigImport;
try {
  invalidActiveConfigImport = await importSession(invalidActiveConfigEvent);
} finally {
  console.error = consoleErrorBeforeInvalidActiveConfig;
}
assert.equal(invalidActiveConfigImport.status, 'error');
assert.equal(invalidActiveConfigImport.error.code, 'INPUT_INVALID');
assert.deepEqual(invalidActiveConfigImport.error.context, {field:'linear_track_layout',reason:'LINEAR_TRACK_LAYOUT'});
assert.deepEqual(compactActiveIntentSnapshot(), stateBeforeInvalidActiveConfig.activeIntent);
assert.strictEqual(state.results.value, stateBeforeInvalidActiveConfig.results);
assert.strictEqual(state.featureCatalog.value, stateBeforeInvalidActiveConfig.featureCatalog);
assert.strictEqual(state.files.c_gb, stateBeforeInvalidActiveConfig.primaryFile);
assert.strictEqual(state.svgContainer.value, stateBeforeInvalidActiveConfig.svgContainer);
assert.equal(state.sessionTitle.value, stateBeforeInvalidActiveConfig.sessionTitle);
assert.equal(invalidActiveConfigEvent.target.value, '');
assert.deepEqual(state.errorLog.value, invalidActiveConfigImport.error);
assert.equal(alerts.length, 0);

const legacyActiveIntent = {
  form: {
    plot_title: 'Legacy JSON next Generate',
    labels_mode: 'both',
    show_scale: false
  },
  adv: {
    axis_stroke_width: 9,
    label_font_size: 29,
    feature_width_circular: 23
  },
  colors: { CDS: '#0b4f6c', tRNA: '#f59e0b' },
  palette: 'orange',
  paletteInstantPreviewEnabled: true,
  rules: [{
    feat: 'CDS',
    qual: 'gene',
    val: '^ND2$',
    color: '#dc2626',
    cap: 'Legacy ND2',
    fromFile: false
  }],
  qualifierPriorityRules: [{ feat: 'CDS', order: 'product,gene' }],
  filterMode: 'Blacklist',
  whitelist: [{ feat: 'CDS', qual: 'gene', key: 'ND2' }],
  blacklistText: 'legacy-hidden'
};
alerts.length = 0;
const legacyActiveIntentImport = await importPayload(legacyActiveIntent);
assert.equal(legacyActiveIntentImport.status, 'legacy');
assert.equal(state.activeDrawing().selectedPalette.value, 'orange');
assert.equal(state.activeDrawing().currentColors.value.CDS, '#0b4f6c');
assert.equal(state.activeDrawing().currentColors.value.tRNA, '#f59e0b');
assert.equal(state.appliedPaletteName.value, 'orange');
assert.equal(state.appliedPaletteColors.value.CDS, '#0b4f6c');
assert.equal(state.activeDrawing().pendingPaletteName.value, '');
assert.deepEqual(state.activeDrawing().manualSpecificRules, legacyActiveIntent.rules);
assert.deepEqual(state.activeDrawing().manualPriorityRules, legacyActiveIntent.qualifierPriorityRules);
assert.equal(state.activeDrawing().filterMode.value, 'Blacklist');
assert.deepEqual(state.activeDrawing().manualWhitelist, legacyActiveIntent.whitelist);
assert.equal(state.activeDrawing().manualBlacklist.value, 'legacy-hidden');
assert.equal(state.activeDrawing().form.plot_title, 'Legacy JSON next Generate');
assert.equal(state.activeDrawing().form.show_scale, false);
assert.equal(state.activeDrawing().adv.axis_stroke_width, 9);
assert.equal(state.activeDrawing().adv.label_font_size, 29);
assert.equal(state.activeDrawing().adv.feature_width_circular, 23);
const legacyActiveFiles = await serializeActiveRenderFiles(state.mode.value, state, state.activeDrawing());
const legacyGeneratedCanonical = buildCanonicalRenderRequest({
  state,
  drawing: state.activeDrawing(),
  filesData: legacyActiveFiles
});
const legacyGeneratedProjection = projectCanonicalSessionRequest(legacyGeneratedCanonical);
assert.equal(legacyGeneratedProjection.config.palette, 'orange');
assert.equal(legacyGeneratedProjection.config.colors.CDS, '#0b4f6c');
assert.equal(legacyGeneratedProjection.config.colors.tRNA, '#f59e0b');
assert.deepEqual(legacyGeneratedProjection.config.rules, legacyActiveIntent.rules.map(
  ({ fromFile: _fromFile, ...rule }) => rule
));
assert.deepEqual(
  legacyGeneratedProjection.config.qualifierPriorityRules,
  legacyActiveIntent.qualifierPriorityRules
);
assert.equal(legacyGeneratedProjection.config.filterMode, 'Blacklist');
assert.equal(legacyGeneratedProjection.config.blacklistText, 'legacy-hidden');
assert.equal(legacyGeneratedProjection.config.form.plot_title, 'Legacy JSON next Generate');
assert.equal(legacyGeneratedProjection.config.form.show_scale, false);
assert.equal(legacyGeneratedProjection.config.adv.axis_stroke_width, 9);
assert.equal(legacyGeneratedProjection.config.adv.label_font_size, 29);
assert.deepEqual(
  legacyGeneratedCanonical.renderRequest.diagramOptions.tracks.circularTrackSlots.find(
    (slot) => slot.renderer === 'features'
  )?.width,
  { value: 23, unit: 'px' }
);
assert.deepEqual(alerts, [
  'Legacy configuration loaded. Save as a session to use the current format.'
]);

state.mode.value = 'circular';
state.activeDrawing().form.multi_record_canvas = false;
state.activeDrawing().form.legend = 'left';
state.activeDrawing().adv.plot_title_position = 'none';
state.generatedLegendPosition.value = 'left';
state.sessionTitle.value = 'keep-after-normalization-error';
Object.assign(state.activeDrawing().canvasPadding, { top: 7, right: 8, bottom: 9, left: 10 });
Object.assign(state.canvasPan, { x: 27, y: 28 });
state.zoom.value = 1.25;
state.skipCaptureBaseConfig.value = false;
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
state.activeDrawing().fileLegendCaptions.value = new Set(['Keep file legend']);
state.featureSearch.value = 'keep feature search';
state.labelSearch.value = 'keep label search';
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
// E1: the other mode's stashed artifact comes back by reference.
const stashedLinearArtifact = Object.freeze({ mode: 'linear' });
const artifactSlots = { circular: null, linear: stashedLinearArtifact };
const displayedArtifact = Object.freeze({ mode: 'circular' });
const installedArtifacts = [];
const sessionImportRollbackState = createSessionImportRollbackState({
  artifactSlots,
  captureDisplayedArtifact: () => displayedArtifact,
  installDisplayedArtifact: (slot) => installedArtifacts.push(slot),
  depthTrackUiCounts,
  depthTracks,
  featureListScrollTop: state.featureListScrollTop,
  featureListScrollRef,
  selectedPairwiseBlockOrthogroupId
});

const jsonClone = (value) => JSON.parse(JSON.stringify(value));
const rollbackState = () => ({
  config: jsonClone(buildConfigData(state.activeDrawing())),
  ui: jsonClone(buildUiStateData(state.activeDrawing())),
  features: jsonClone(buildFeatureStateData(state.activeDrawing())),
  editorState: jsonClone(buildEditorStateData(state.activeDrawing())),
  orthogroupState: jsonClone(buildOrthogroupStateData(state.activeDrawing())),
  runState: jsonClone(buildRunStateData()),
  results: jsonClone(state.results.value),
  mode: state.mode.value,
  title: state.sessionTitle.value,
  legend: state.activeDrawing().form.legend,
  plotTitlePosition: state.activeDrawing().adv.plot_title_position,
  generatedLegendPosition: state.generatedLegendPosition.value,
  semanticFileWatchersSuppressed: state.semanticFileWatchersSuppressed.value,
  sessionImportRollbackInProgress: state.sessionImportRollbackInProgress.value,
  skipCaptureBaseConfig: state.skipCaptureBaseConfig.value,
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
  fileLegendCaptions: [...state.activeDrawing().fileLegendCaptions.value],
  featureSearch: state.featureSearch.value,
  labelSearch: state.labelSearch.value,
  selectedFeatureIds: [...state.selectedFeatureIds.value],
  selectedFeatureAnchorId: state.selectedFeatureAnchorId.value,
  featureSelectionStatus: state.featureSelectionStatus.value,
  featureSelectionSuppressNextClick: state.featureSelectionSuppressNextClick.value,
  featureSelectionDrag: structuredClone(state.featureSelectionDrag),
  labelReflowLastError: state.labelReflowLastError.value,
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
// F-1 (R11): the rollback installs the captured state as it was, so unset
// slot sides, lane directions, and axis indexes stay unset and the Result's
// named stroke color stays named.
state.activeDrawing().adv.circular_track_slots.forEach((slot) => {
  slot.side = null;
  delete slot.params.lane_direction;
});
state.activeDrawing().adv.circular_track_slots_axis_index = null;
state.activeDrawing().adv.linear_track_slots_axis_index = null;
state.originalSvgStroke.value = { color: 'gray', width: 1 };
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
        artifactSlots.linear = null;
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
assert.equal(failedImport.error.code, 'INPUT_INVALID');
assert.deepEqual(failedImport.error.context, {field:'schema',reason:'JSON_FORMAT'});
assert.deepEqual(rollbackState(), stateBeforeFailedImport);
assert.strictEqual(artifactSlots.linear, stashedLinearArtifact);
assert.deepEqual(installedArtifacts, [displayedArtifact]);
assert.equal(artifactSlots.circular, null);
assert.equal(alerts.length, 0);
assert.deepEqual(state.errorLog.value, failedImport.error);
assert.equal(failedImportEvent.target.value, '');


// A late failed Session read cannot replace a notification from a later action.
const { normalizeUserFacingError } = await import('../../gbdraw/web/js/utils/error-normalization.js');
const beforeLateRead = rollbackState();
const lateFile = new File(['{}'], 'PRIVATE_LATE_SESSION.json');
let failLateRead;
// File methods do not cross structured clone. Hold the actual Worker reply
// boundary so this fixture exercises a late transport failure, not an ignored
// main-thread stream override.
const SessionWorker = globalThis.Worker;
globalThis.Worker = class extends SessionWorker {
  postMessage(message) {
    if (message.file === lateFile) {
      failLateRead = () => this.emit('error', { message: 'PRIVATE_LATE_READ_SENTINEL' });
    } else super.postMessage(message);
  }
};
const lateRead = importSession({target:{files:[lateFile],value:'selected'}});
while (!failLateRead) await new Promise(resolve=>setTimeout(resolve,0));
globalThis.Worker = SessionWorker;
const laterAlert = normalizeUserFacingError({code:'PDF_LIBRARY',operation:'export-pdf',stage:'initialization'});
state.errorLog.value = laterAlert;
failLateRead();
assert.deepEqual(await lateRead,{status:'stale'});
assert.equal(state.errorLog.value,laterAlert);
assert.deepEqual(rollbackState(),beforeLateRead);

// The same sole Session owner keeps a late Save error from replacing a newer alert.
const { exportSession } = await import('../../gbdraw/web/js/services/config.js');
let failLateSave;
const lateSave = exportSession('late save', { beforeExport: () => new Promise((_, reject) => {
  failLateSave = reject;
}), onError: () => assert.fail('A stale Save must not publish its error') });
while (!failLateSave) await new Promise(resolve => setTimeout(resolve, 0));
const saveLaterAlert = normalizeUserFacingError({ code: 'PDF_LIBRARY', operation: 'export-pdf', stage: 'initialization' });
state.errorLog.value = saveLaterAlert;
failLateSave(new Error('PRIVATE_LATE_SAVE_SENTINEL'));
assert.deepEqual(await lateSave, { status: 'stale' });
assert.equal(state.errorLog.value, saveLaterAlert);
assert.equal(state.sessionSavePending.value, false);
