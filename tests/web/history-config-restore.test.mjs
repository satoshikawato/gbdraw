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
  state[domain][key] = null;
  await history.initializeIntentBaseline();
  await history.runUndoable('Set explicit value', () => { state[domain][key] = first; });
  await history.runUndoable('Change explicit value', () => { state[domain][key] = second; });
  await history.undo();
  assert.equal(state[domain][key], first, `${domain}.${key}: explicit to explicit`);
  await history.undo();
  assert.equal(state[domain][key], null, `${domain}.${key}: explicit to Auto`);
  await history.redo();
  assert.equal(state[domain][key], first, `${domain}.${key}: Redo explicit value`);
}

applyConfigData({ form: JSON.parse('{"unknown":1,"__proto__":{"polluted":true}}') });
assert.equal(Object.hasOwn(state.form, 'unknown'), false);
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
console.log('History restores nullable config values and preserves key guards.');

// SE-01, N-19, N-20: History checkpoints and the Session rollback hold the
// admitted feature catalog by reference; state admits no other catalog.
{
  const { admitFeatureCatalog, featureStateFromCatalog } = await import('../../gbdraw/web/js/services/feature-catalog.js');
  const { buildEditorStateData, applyEditorStateData } = await import('../../gbdraw/web/js/services/config.js');
  const marker = 'se01-catalog-payload';
  const catalog = admitFeatureCatalog({
    schema: 4,
    items: [{
      resultIndex: 0, resultName: 'diagram.svg', recordKeys: ['record-a'],
      features: [{ svgId: 'f0001', recordKey: 'record-a', biologicalFeatureId: 'feature-a', fillColor: '#abcdef' }],
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
    state.legendEntries.value = [{ caption: 'SE-01', color: '#123456' }];
  });
  for (const direction of ['undo', 'redo', 'undo']) {
    await artifactHistory[direction]();
    assert.strictEqual(liveCatalog(), catalog, `${direction} keeps the admitted catalog object`);
    assert.equal(featureStateFromCatalog(liveCatalog(), { mode: 'circular' }).extractedFeatures.length, 1);
  }

  applyEditorStateData(buildEditorStateData({ preserveAdoptedCatalog: true }));
  assert.strictEqual(liveCatalog(), catalog, 'Session rollback restores the admitted catalog');
  applyEditorStateData({ featureCatalog: structuredClone(catalog) });
  assert.equal(state.featureCatalog.value, null, 'an unadmitted catalog never enters state');
  state.featureCatalog.value = null;
  state.results.value = [];
  state.legendEntries.value = [];
  console.log('History and Session rollback keep the admitted feature catalog by reference.');
}
