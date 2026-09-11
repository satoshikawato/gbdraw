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
console.log('History restores nullable config values and preserves key guards.');
