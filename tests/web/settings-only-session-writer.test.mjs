import assert from 'node:assert/strict';
import { test } from 'node:test';
import { gunzipSync } from 'node:zlib';

// Exercise the existing coordinator and gzip writer with browser I/O stubbed.
globalThis.window = {
  Vue: {
    ref: value => ({ value }), reactive: value => value,
    computed: getter => ({ get value() { return getter(); } }),
    nextTick: async () => {}
  },
  DOMPurify: { sanitize: value => value }
};
globalThis.document = {
  body: { appendChild: () => {} },
  createElement: () => ({ addEventListener: () => {}, click: () => {}, parentNode: null })
};
const { exportSession, buildConfigData } = await import('../../gbdraw/web/js/services/config.js');
const { state } = await import('../../gbdraw/web/js/state.js');
const { adoptCurrentSessionDocument } = await import('../../gbdraw/web/js/services/session-authority.js');

const savedDocument = async title => {
  const result = await exportSession(title);
  assert.equal(result.status, 'saved');
  return JSON.parse(gunzipSync(Buffer.from(await result.blob.arrayBuffer())));
};

test('source-free Save emits no render metadata and preserves valid raw scalar drafts', async () => {
  const freshConfig = buildConfigData();
  const fresh = await savedDocument('fresh settings');
  assert.equal(fresh.renderRequest, null);
  assert.deepEqual(fresh.runMetadata, {});
  assert.deepEqual(fresh.config, JSON.parse(JSON.stringify(freshConfig)));
  assert.deepEqual(buildConfigData(), freshConfig);
  assert.equal(adoptCurrentSessionDocument(fresh, 44).canonical, null);

  state.adv.circular_track_slots_enabled = true;
  state.adv.circular_track_slots.splice(0, state.adv.circular_track_slots.length, {
    id: 'features', renderer: 'features', enabled: true, side: 'inside',
    width: { value: '1.', unit: 'px' }, radius: { value: '1e-3', unit: 'factor' },
    inner_gap_px: null, outer_gap_px: null, z: 0, params: { lane_direction: 'inside' }
  });
  const before = buildConfigData();
  const saved = await savedDocument('typed text settings');
  assert.deepEqual(saved.config, JSON.parse(JSON.stringify(before)));
  assert.deepEqual(buildConfigData(), before);
  assert.deepEqual(saved.results, []);
  assert.equal(saved.editorState.featureCatalog, null);
  assert.deepEqual(saved.runMetadata, {});
  assert.equal(adoptCurrentSessionDocument(saved, 44).canonical, null);

  const forbiddenMetadata = structuredClone(saved);
  forbiddenMetadata.runMetadata.annotationWarnings = [];
  assert.throws(() => adoptCurrentSessionDocument(forbiddenMetadata, 44), /committed render artifacts/);
  state.adv.circular_track_slots[0].width = { value: '1e', unit: 'px' };
  await assert.rejects(exportSession('unfinished settings'), /width|positive|scalar/i);
  assert.deepEqual(state.adv.circular_track_slots[0].width, { value: '1e', unit: 'px' });
});
