import assert from 'node:assert/strict';
import { test } from 'node:test';
import { gunzipSync } from 'node:zlib';
import { readFile } from 'node:fs/promises';

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

const { writeCircularMeasureValue, changeCircularMeasureUnit } = await import(
  '../../gbdraw/web/js/app/circular-track-slots/measure-editor.js'
);
const { buildCircularTrackSlotPayload } = await import('../../gbdraw/web/js/app/circular-track-slots.js');
const { projectSettingsOnlySession } = await import('../../gbdraw/web/js/services/session-request.js');
const scalarFixtures = JSON.parse(await readFile(new URL(
  '../../docs/internal/issue-619-implementation-plan-20260927/SESSION_RESULTS/scalar-fixtures.json', import.meta.url
)));

test('current gzip writer/admission keeps S00 valid scalars and codec drafts, including disabled/inactive rows', async () => {
  const cases = scalarFixtures.filter(fixture => fixture.valid).map(fixture => [fixture.input, fixture.canonical]);
  cases.push(
    [writeCircularMeasureValue('1.', 'px'), { value: 1, unit: 'px' }],
    [changeCircularMeasureUnit({ value: '1e-3', unit: 'factor' }, 'px'), { value: 0.001, unit: 'px' }],
    [writeCircularMeasureValue('65%', 'px'), { value: 0.65, unit: 'factor' }]
  );
  for (const [mode, enabled] of [['circular', true], ['circular', false], ['linear', false]]) {
    state.mode.value = mode;
    const slot = state.adv.circular_track_slots[0];
    slot.enabled = enabled;
    for (const [index, [scalar, expected]] of cases.entries()) {
      slot.width = structuredClone(scalar);
      slot.radius = structuredClone(scalar);
      const before = buildConfigData();
      const saved = await savedDocument(`scalar settings ${mode} ${enabled} ${index}`);
      const adopted = adoptCurrentSessionDocument(saved, 44);
      assert.equal(adopted.canonical, null);
      const restored = projectSettingsOnlySession(adopted.document);
      assert.deepEqual(restored.config, JSON.parse(JSON.stringify(before)));
      assert.deepEqual(buildConfigData(), before);
      const restoredSlot = restored.config.adv.circular_track_slots[0];
      assert.deepEqual(restoredSlot.width, scalar);
      assert.deepEqual(restoredSlot.radius, scalar);
      assert.deepEqual(buildCircularTrackSlotPayload(restoredSlot).width, expected);
      assert.deepEqual(buildCircularTrackSlotPayload(restoredSlot).radius, expected);
      assert.deepEqual(saved.runMetadata, {});
      assert.deepEqual(saved.results, []);
    }
  }
});

test('writer and current admission reject invalid scalar drafts without altering them', async () => {
  state.mode.value = 'circular';
  const slot = state.adv.circular_track_slots[0];
  slot.enabled = true;
  slot.width = writeCircularMeasureValue('1.', 'px');
  slot.radius = null;
  const valid = await savedDocument('valid scalar settings');
  for (const scalar of [
    writeCircularMeasureValue('1e', 'px'), writeCircularMeasureValue('0', 'factor'),
    writeCircularMeasureValue('-1', 'px'), writeCircularMeasureValue('Infinity', 'factor'),
    { value: '', unit: 'px' }, { value: 1, unit: 'em' }, { value: true, unit: 'px' },
    { value: Infinity, unit: 'px' }, Infinity, NaN
  ]) {
    slot.width = scalar;
    const before = buildConfigData();
    await assert.rejects(exportSession('invalid scalar settings'), /width|positive|scalar/i);
    assert.deepEqual(buildConfigData(), before);
    const invalid = structuredClone(valid);
    invalid.config.adv.circular_track_slots[0].width = scalar;
    const invalidBefore = structuredClone(invalid);
    assert.throws(() => adoptCurrentSessionDocument(invalid, 44), /width|positive|scalar/i);
    assert.deepEqual(invalid, invalidBefore);
  }
});
