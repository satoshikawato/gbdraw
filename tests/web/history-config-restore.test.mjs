import assert from 'node:assert/strict';

globalThis.window = { Vue: {
  ref: value => ({ value }), reactive: value => value,
  computed: getter => ({ get value() { return getter(); } }), nextTick: async () => {}
} };
const { state } = await import('../../gbdraw/web/js/state.js');
const { buildConfigData, applyConfigData } = await import('../../gbdraw/web/js/services/config.js');
const { createHistoryManager } = await import('../../gbdraw/web/js/services/history.js');
const { createHistoryFileStore } = await import('../../gbdraw/web/js/services/history-files.js');
const { createHistorySnapshotService } = await import('../../gbdraw/web/js/services/history-snapshot.js');
const snapshots = createHistorySnapshotService({
  state, fileStore: createHistoryFileStore(), buildConfigData, applyConfigData
});
const history = createHistoryManager({
  buildIntent: snapshots.buildHistoryIntent, applyIntent: snapshots.applyHistoryIntent,
  buildCheckpoint: () => assert.fail('Form edits must use intent History'),
  applyCheckpoint: () => assert.fail('Form edits must use intent History')
});

for (const [domain, key, first, second] of [
  ['form', 'circular_region_start', 1000, 2000],
  ['form', 'circular_region_end', 500, 1500],
  ['adv', 'window_size', 100, 200],
  ['adv', 'block_stroke_color', '#123456', '#abcdef']
]) {
  state.activeDrawing()[domain][key] = null;
  await history.initializeIntentBaseline();
  await history.runUndoable('Set explicit value', () => { state.activeDrawing()[domain][key] = first; });
  await history.runUndoable('Change explicit value', () => { state.activeDrawing()[domain][key] = second; });
  await history.undo();
  assert.equal(state.activeDrawing()[domain][key], first, `${domain}.${key}: explicit to explicit`);
  await history.undo();
  assert.equal(state.activeDrawing()[domain][key], null, `${domain}.${key}: explicit to Auto`);
  await history.redo();
  assert.equal(state.activeDrawing()[domain][key], first, `${domain}.${key}: Redo explicit value`);
}

// Record selection: a toggle is one intent step, and Undo and Redo restore
// the OFF records of the drawing it was made in.
{
  const drawing = state.activeDrawing();
  drawing.recordsOff.splice(0);
  await history.initializeIntentBaseline();
  await history.runUndoable('Turn a record OFF', () => { drawing.recordsOff.push('#2'); });
  await history.runUndoable('Turn another record OFF', () => { drawing.recordsOff.push('#3'); });
  await history.undo();
  assert.deepEqual([...drawing.recordsOff], ['#2']);
  await history.undo();
  assert.deepEqual([...drawing.recordsOff], []);
  await history.redo();
  assert.deepEqual([...drawing.recordsOff], ['#2']);
  drawing.recordsOff.splice(0);
}

// Delete settings (record selection D-08) is one step: Undo restores the
// card settings and every edit of the record it deleted.
{
  const config = await import('../../gbdraw/web/js/services/config.js');
  const { RECORD_SETTINGS_DEFAULTS, removeRecordEdits } = await import('../../gbdraw/web/js/services/record-draw-selection.js');
  const editorSnapshots = createHistorySnapshotService({
    state, fileStore: createHistoryFileStore(), buildConfigData, applyConfigData,
    buildFeatureStateData: config.buildFeatureStateData, applyFeatureStateData: config.applyFeatureStateData,
    buildEditorStateData: config.buildEditorStateData, applyEditorStateData: config.applyEditorStateData
  });
  const editorHistory = createHistoryManager({
    buildIntent: editorSnapshots.buildHistoryIntent, applyIntent: editorSnapshots.applyHistoryIntent,
    buildCheckpoint: () => assert.fail('Delete settings must use intent History'),
    applyCheckpoint: () => assert.fail('Delete settings must use intent History')
  });
  const drawing = state.drawings.linear;
  const card = state.linearSeqs[0];
  const uid = card.uid;
  card.definition = 'Kept organism';
  card.losat_gencode = 11;
  drawing.recordsOff.splice(0, Infinity, uid);
  const key = JSON.stringify([uid, 'cds-1']);
  drawing.featureOverrides[key] = {
    recordKey: uid, biologicalFeatureId: 'cds-1', featureVisibility: 'off', labelVisibility: null, labelText: null, labelSourceText: null
  };
  drawing.annotationSets.splice(0, Infinity, { id: 'set', label: 'Set', annotations: [
    { id: 'feature_1', target: { kind: 'featureIdentity', recordKey: uid, biologicalFeatureId: 'cds-1' },
      label: '', mark: 'highlight', lane: null, style: null, legendLabel: null, metadata: {} }
  ] });
  await editorHistory.initializeIntentBaseline();
  await editorHistory.runUndoable('Delete record settings', () => {
    removeRecordEdits([{ key: uid, requestKeys: [uid], ownsExpansions: true, bindingKeys: [], displaySource: uid, displaySelector: null }],
      { featureOverrides: drawing.featureOverrides, annotationSets: drawing.annotationSets, annotationBindingField: 'binding' });
    Object.assign(card, RECORD_SETTINGS_DEFAULTS);
  });
  assert.equal(card.definition, '');
  assert.equal(drawing.featureOverrides[key], undefined);
  assert.equal(drawing.annotationSets[0].annotations.length, 0);
  await editorHistory.undo();
  const restored = state.linearSeqs.find((sequence) => sequence.uid === uid);
  assert.equal(restored.definition, 'Kept organism');
  assert.equal(restored.losat_gencode, 11);
  assert.equal(drawing.featureOverrides[key]?.featureVisibility, 'off');
  assert.deepEqual(drawing.annotationSets[0].annotations.map((item) => item.id), ['feature_1']);
  assert.deepEqual([...drawing.recordsOff], [uid]);
  drawing.recordsOff.splice(0);
  drawing.annotationSets.splice(0);
  delete drawing.featureOverrides[key];
}

applyConfigData(state.activeDrawing(), { form: JSON.parse('{"unknown":1,"__proto__":{"polluted":true}}') });
assert.equal(Object.hasOwn(state.activeDrawing().form, 'unknown'), false);
assert.equal({}.polluted, undefined);

state.similarityAlignmentPlan.value = {
  schema: 1,
  mode: 'position',
  groupId: 'og-history',
  reference: {
    recordKey: 'record-a', biologicalFeatureId: 'feature-a',
    sourceFeatureIndex: 0, stableFeatureSvgId: 'feature-a'
  },
  records: []
};
state.linearRecordTranslations.value = [{ recordKey: 'record-a', x: 12, y: -3 }];
await history.initializeIntentBaseline();
await history.runUndoable('Clear alignment', () => {
  state.similarityAlignmentPlan.value = null;
  state.linearRecordTranslations.value = [{ recordKey: 'record-a', x: 21, y: 4 }];
});
await history.undo();
assert.equal(state.similarityAlignmentPlan.value.groupId, 'og-history');
assert.deepEqual(state.linearRecordTranslations.value, [
  { recordKey: 'record-a', x: 12, y: -3 }
]);
await history.redo();
assert.equal(state.similarityAlignmentPlan.value, null);
assert.deepEqual(state.linearRecordTranslations.value, [
  { recordKey: 'record-a', x: 21, y: 4 }
]);
// OV-383: Undo and Redo of a Linear setting keep the active alignment. The
// config a History step restores carries no alignment, so restoring it must
// not clear the plan or the record translations.
state.mode.value = 'linear';
const keptPlan = { schema: 1, mode: 'position', groupId: 'og-ov383',
  reference: { recordKey: 'record-a', biologicalFeatureId: 'feature-a', sourceFeatureIndex: 0, stableFeatureSvgId: 'feature-a' },
  records: [] };
const keptTranslations = [{ recordKey: 'record-a', x: 12, y: -3 }];
state.similarityAlignmentPlan.value = keptPlan;
state.linearRecordTranslations.value = keptTranslations;
const shaftBefore = state.drawings.linear.adv.arrow_shaft_width_ratio;
await history.initializeIntentBaseline();
await history.runUndoable('Arrow shaft', () => { state.drawings.linear.adv.arrow_shaft_width_ratio = 0.6; });
await history.undo();
assert.equal(state.drawings.linear.adv.arrow_shaft_width_ratio, shaftBefore);
assert.equal(state.similarityAlignmentPlan.value?.groupId, 'og-ov383', 'Undo keeps the alignment plan');
assert.deepEqual(state.linearRecordTranslations.value, keptTranslations, 'Undo keeps the record translations');
await history.redo();
assert.equal(state.drawings.linear.adv.arrow_shaft_width_ratio, 0.6);
assert.equal(state.similarityAlignmentPlan.value?.groupId, 'og-ov383', 'Redo keeps the alignment plan');
assert.deepEqual(state.linearRecordTranslations.value, keptTranslations, 'Redo keeps the record translations');
console.log('History restores nullable config values and preserves key guards.');

// SE-01, N-19, N-20: History checkpoints and the Session rollback hold the
// admitted feature catalog by reference; state admits no other catalog.
{
  const { FEATURE_CATALOG_SCHEMA, admitFeatureCatalog, featureStateFromCatalog } = await import('../../gbdraw/web/js/services/feature-catalog.js');
  const { buildEditorStateData, applyEditorStateData } = await import('../../gbdraw/web/js/services/config.js');
  const marker = 'se01-catalog-payload';
  const catalog = admitFeatureCatalog({
    schema: FEATURE_CATALOG_SCHEMA,
    items: [{
      resultIndex: 0, resultName: 'diagram.svg', recordKeys: ['record-a'],
      features: [{
        svgId: 'f0001', recordKey: 'record-a', biologicalFeatureId: 'feature-a', fillColor: '#abcdef',
        drawnSelector: { hash: 'stable-a', location: '1..6', recordLocation: 'record-a:1..6:+' }
      }],
      biologicalFeatures: [{
        recordKey: 'record-a', biologicalFeatureId: 'feature-a', stableFeatureId: 'stable-a', record_idx: 0,
        sourceFeatureIndex: 0, record_id: 'record-a', type: 'CDS', start: 1, end: 6, strand: 1,
        anchorProfile: { precision: 'exact', operator: 'single', partOrder: 'biological', strand: '+' },
        qualifiers: { note: [marker] }
      }],
      orthogroups: [], annotations: [], comparisonMatches: [], sequenceSources: []
    }]
  }, [{ name: 'diagram.svg', content: '<svg />' }], { adopt: true, mode: 'circular' }).catalog;
  state.mode.value = 'circular';
  state.results.value = [{ name: 'diagram.svg', content: '<svg />' }];
  state.featureCatalog.value = catalog;
  const liveCatalog = () => window.Vue.toRaw?.(state.featureCatalog.value) ?? state.featureCatalog.value;
  assert.strictEqual(liveCatalog(), catalog);

  const artifactSnapshots = createHistorySnapshotService({
    state, fileStore: createHistoryFileStore(), buildConfigData, applyConfigData,
    buildEditorStateData, applyEditorStateData
  });
  const artifactHistory = createHistoryManager({
    buildIntent: artifactSnapshots.buildHistoryIntent,
    applyIntent: artifactSnapshots.applyHistoryIntent,
    buildCheckpoint: artifactSnapshots.buildArtifactCheckpoint,
    applyCheckpoint: artifactSnapshots.applyArtifactCheckpoint,
    signatureFor: artifactSnapshots.snapshotSignature
  });
  await artifactHistory.captureBaseline();
  assert.equal(JSON.stringify(artifactHistory.getCurrentCheckpoint()).includes(marker), false);
  await artifactHistory.runUndoableCheckpoint('Change legend', () => {
    state.activeDrawing().legendEntries.value = [{ caption: 'SE-01', color: '#123456' }];
  });
  for (const direction of ['undo', 'redo', 'undo']) {
    await artifactHistory[direction]();
    assert.strictEqual(liveCatalog(), catalog, `${direction} keeps the admitted catalog object`);
    assert.equal(featureStateFromCatalog(liveCatalog(), { mode: 'circular' }).extractedFeatures.length, 1);
  }

  applyEditorStateData(state.activeDrawing(), buildEditorStateData(state.activeDrawing()));
  assert.strictEqual(liveCatalog(), catalog, 'Session rollback restores the admitted catalog');
  applyEditorStateData(state.activeDrawing(), { featureCatalog: structuredClone(catalog) });
  assert.equal(state.featureCatalog.value, null, 'an unadmitted catalog never enters state');
  state.featureCatalog.value = null;
  state.results.value = [];
  state.activeDrawing().legendEntries.value = [];
  console.log('History and Session rollback keep the admitted feature catalog by reference.');
}

// F-1 (G-H, R11): an Undo or Redo installs the captured state exactly. Unset
// slot sides, the Features lane direction, both axis indexes, and the Circular
// multi-record layout stay unset, and the Result's named stroke color stays
// named, so the next Generate and Session Save see the same state as before.
{
  const {
    buildEditorStateData, applyEditorStateData, buildUiStateData, applyUiStateData
  } = await import('../../gbdraw/web/js/services/config.js');
  state.mode.value = 'circular';
  state.activeDrawing().adv.circular_track_slots.forEach((slot) => {
    slot.side = null;
    delete slot.params.lane_direction;
  });
  state.activeDrawing().adv.circular_track_slots_axis_index = null;
  state.activeDrawing().adv.linear_track_slots_axis_index = null;
  Object.assign(state.activeDrawing().layoutPreferences.circular.multi, { legend: null, plotTitlePosition: null });
  state.originalSvgStroke.value = { color: 'gray', width: 1 };
  const unsetValues = () => ({
    sides: state.activeDrawing().adv.circular_track_slots.map((slot) => slot.side),
    laneDirection: state.activeDrawing().adv.circular_track_slots.map((slot) => slot.params.lane_direction ?? '<unset>'),
    circularAxis: state.activeDrawing().adv.circular_track_slots_axis_index,
    linearAxis: state.activeDrawing().adv.linear_track_slots_axis_index,
    multiLayout: { ...state.activeDrawing().layoutPreferences.circular.multi },
    originalSvgStroke: { ...state.originalSvgStroke.value }
  });
  const unset = unsetValues();
  assert.deepEqual(unset, {
    sides: [null, null, null, null],
    laneDirection: ['<unset>', '<unset>', '<unset>', '<unset>'],
    circularAxis: null,
    linearAxis: null,
    multiLayout: { legend: null, plotTitlePosition: null },
    originalSvgStroke: { color: 'gray', width: 1 }
  });
  const settings = () => JSON.stringify({
    config: buildConfigData(state.activeDrawing()), layoutPreferences: buildUiStateData(state.activeDrawing()).layoutPreferences
  });
  const restoreSnapshots = createHistorySnapshotService({
    state, fileStore: createHistoryFileStore(), buildConfigData, applyConfigData,
    buildUiStateData, applyUiStateData, buildEditorStateData, applyEditorStateData
  });
  const restoreHistory = createHistoryManager({
    buildIntent: restoreSnapshots.buildHistoryIntent,
    applyIntent: restoreSnapshots.applyHistoryIntent,
    buildCheckpoint: restoreSnapshots.buildArtifactCheckpoint,
    applyCheckpoint: restoreSnapshots.applyArtifactCheckpoint,
    signatureFor: restoreSnapshots.snapshotSignature
  });
  await restoreHistory.captureBaseline();
  await restoreHistory.initializeIntentBaseline();
  await restoreHistory.runUndoable('Rich Feature Popup', () => {
    state.richFeaturePopup.value = !state.richFeaturePopup.value;
  });
  await restoreHistory.runUndoable('Label Mode', () => { state.activeDrawing().form.labels_mode = 'both'; });
  const edited = settings();
  for (const direction of ['undo', 'redo', 'undo', 'undo', 'redo', 'redo']) {
    await restoreHistory[direction]();
    assert.deepEqual(unsetValues(), unset, `intent ${direction} keeps unset values unset`);
  }
  assert.equal(settings(), edited, 'Undo and Redo return to the edited settings');

  await restoreHistory.runUndoableCheckpoint('Change legend', () => {
    state.activeDrawing().legendEntries.value = [{ caption: 'F-1', color: '#123456' }];
  });
  const checkpointed = settings();
  for (const direction of ['undo', 'redo']) {
    await restoreHistory[direction]();
    assert.deepEqual(unsetValues(), unset, `checkpoint ${direction} keeps unset values unset`);
  }
  assert.equal(settings(), checkpointed, 'checkpoint Undo and Redo return to the same settings');

  // A restored stack belongs to state: an edit after the checkpoint Redo is
  // its own step and never writes into the checkpoint it was restored from.
  const steps = restoreHistory.getUndoCount();
  await restoreHistory.runUndoable('Move features outside', () => {
    state.activeDrawing().adv.circular_track_slots[0].side = 'outside';
  });
  assert.equal(restoreHistory.getUndoCount(), steps + 1, 'an edit after a restore adds one step');
  for (const direction of ['undo', 'undo', 'redo']) {
    await restoreHistory[direction]();
    assert.deepEqual(unsetValues(), unset, `${direction} after a later edit keeps the checkpoint unset`);
  }
  state.activeDrawing().legendEntries.value = [];
  state.originalSvgStroke.value = { color: null, width: null };
  console.log('History Undo and Redo keep unset settings unset.');
}
