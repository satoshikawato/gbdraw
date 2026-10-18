import assert from 'node:assert/strict';
import { test } from 'node:test';
import { gunzipSync } from 'node:zlib';
import { readFile } from 'node:fs/promises';
import { installSessionImportWorker } from './helpers/session-import-node.mjs';
import { callerSession } from './helpers/session-caller.mjs';

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
const sessionConfig = await import('../../gbdraw/web/js/services/config.js');
const { SESSION_VERSION, buildConfigData, installShownDrawing } = sessionConfig;
const { state } = await import('../../gbdraw/web/js/state.js');
const { exportSession, importSession } = callerSession(sessionConfig, state);
globalThis.alert = () => {};
installSessionImportWorker();
const { adoptCurrentSessionDocument, validateSessionAuthorityInventory } = await import('../../gbdraw/web/js/services/session-authority.js');

const savedDocument = async title => {
  const result = await exportSession(title);
  assert.equal(result.status, 'saved');
  return JSON.parse(gunzipSync(Buffer.from(await result.blob.arrayBuffer())));
};

// These files are refused before the Load dialog; a dialog would cancel the
// Load and fail the check.
const refuseLoadDialog = async () => false;

// E-OD02, E-OD06: Save refuses a settings-only Session 46 without the
// settings of its mode; Load reads a missing `form` or `adv` as its default.
test('OV-591: Save writes form and adv in the settings of its mode, as Python requires; Load reads them as defaults', async () => {
  // gbdraw/session_io.py `_validate_settings_only_session`: "requires an active Web configuration and mode".
  const fresh = await savedDocument('settings of its mode');
  for (const [name, mutate] of [
    ['config {}', d => { d.modes.circular.config = {}; }],
    ['no form', d => { delete d.modes.circular.config.form; }],
    ['no adv', d => { delete d.modes.circular.config.adv; }]
  ]) {
    const document = structuredClone(fresh);
    mutate(document);
    assert.throws(() => validateSessionAuthorityInventory(structuredClone(document), SESSION_VERSION), /settings of its mode/, name);
    const loaded = await importSession({
      target: { files: [new Blob([JSON.stringify(document)], { type: 'application/json' })], value: 'selected' }
    }, { confirmLoadRepairs: refuseLoadDialog });
    assert.equal(loaded.status, 'ok', name);
    assert.deepEqual(loaded.repairs ?? [], [], name);
  }
  // The other mode's slice may hold less: Load fills it with that mode's defaults.
  const other = structuredClone(fresh);
  other.modes.linear.config = {};
  assert.equal(adoptCurrentSessionDocument(other, SESSION_VERSION).canonical, null);
  validateSessionAuthorityInventory(other, SESSION_VERSION);
});

test('OV-612: Save names the mode slice of a settings-only Session 46 by a string that is a key of modes', async () => {
  // `modes[ui.mode]` converted `["circular"]` to the key "circular"; Python requires `ui.mode` to be a diagram mode.
  const fresh = await savedDocument('mode key');
  for (const mode of [['circular'], [['linear']], 'constructor', '__proto__', 0]) {
    const document = structuredClone(fresh);
    document.ui.mode = mode;
    assert.throws(() => validateSessionAuthorityInventory(document, SESSION_VERSION), /settings of its mode/, JSON.stringify(mode));
  }
});

test('source-free Save emits no render metadata and preserves valid raw scalar drafts', async () => {
  const freshConfig = buildConfigData(state.activeDrawing());
  const fresh = await savedDocument('fresh settings');
  assert.equal(fresh.renderRequest, null);
  assert.equal(fresh.runMetadata, undefined);
  // The shown mode's drawing is its slice (Session 46).
  assert.deepEqual(fresh.modes.circular.config, JSON.parse(JSON.stringify(freshConfig)));
  assert.deepEqual(buildConfigData(state.activeDrawing()), freshConfig);
  assert.equal(adoptCurrentSessionDocument(fresh, SESSION_VERSION).canonical, null);

  state.activeDrawing().adv.circular_track_slots_enabled = true;
  state.activeDrawing().adv.circular_track_slots.splice(0, state.activeDrawing().adv.circular_track_slots.length, {
    id: 'features', renderer: 'features', enabled: true, side: 'inside',
    width: { value: '1.', unit: 'px' }, radius: { value: '1e-3', unit: 'factor' },
    inner_gap_px: null, outer_gap_px: null, z: 0, params: { lane_direction: 'inside' }
  });
  const before = buildConfigData(state.activeDrawing());
  const saved = await savedDocument('typed text settings');
  assert.deepEqual(saved.modes.circular.config, JSON.parse(JSON.stringify(before)));
  assert.deepEqual(buildConfigData(state.activeDrawing()), before);
  assert.deepEqual(saved.results, []);
  assert.equal(saved.editorState.featureCatalog, null);
  assert.equal(saved.runMetadata, undefined);
  assert.equal(adoptCurrentSessionDocument(saved, SESSION_VERSION).canonical, null);

  const forbiddenMetadata = structuredClone(saved);
  forbiddenMetadata.runMetadata = { annotationWarnings: [] };
  assert.throws(() => adoptCurrentSessionDocument(forbiddenMetadata, SESSION_VERSION), /committed render artifacts/);
  state.activeDrawing().adv.circular_track_slots[0].width = { value: '1e', unit: 'px' };
  await assert.rejects(exportSession('unfinished settings'), /width|positive|scalar/i);
  assert.deepEqual(state.activeDrawing().adv.circular_track_slots[0].width, { value: '1e', unit: 'px' });
});

const { writeCircularMeasureValue, changeCircularMeasureUnit } = await import(
  '../../gbdraw/web/js/services/circular-track-measure.js'
);
const { buildCircularTrackSlotPayload } = await import('../../gbdraw/web/js/services/circular-track-slot-model.js');
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
    installShownDrawing(mode);
    const slot = state.activeDrawing().adv.circular_track_slots[0];
    slot.enabled = enabled;
    for (const [index, [scalar, expected]] of cases.entries()) {
      slot.width = structuredClone(scalar);
      slot.radius = structuredClone(scalar);
      const before = buildConfigData(state.activeDrawing());
      const saved = await savedDocument(`scalar settings ${mode} ${enabled} ${index}`);
      const adopted = adoptCurrentSessionDocument(saved, SESSION_VERSION);
      assert.equal(adopted.canonical, null);
      const restored = projectSettingsOnlySession(adopted.document);
      assert.deepEqual(restored.config, JSON.parse(JSON.stringify(before)));
      assert.deepEqual(buildConfigData(state.activeDrawing()), before);
      const restoredSlot = restored.config.adv.circular_track_slots[0];
      assert.deepEqual(restoredSlot.width, scalar);
      assert.deepEqual(restoredSlot.radius, scalar);
      assert.deepEqual(buildCircularTrackSlotPayload(restoredSlot).width, expected);
      assert.deepEqual(buildCircularTrackSlotPayload(restoredSlot).radius, expected);
      assert.equal(saved.runMetadata, undefined);
      assert.deepEqual(saved.results, []);
    }
  }
});

test('writer and current admission reject invalid scalar drafts without altering them', async () => {
  installShownDrawing('circular');
  const slot = state.activeDrawing().adv.circular_track_slots[0];
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
    // A rejected Load restores the drawing's slots, so the slot is read again.
    state.activeDrawing().adv.circular_track_slots[0].width = scalar;
    const before = buildConfigData(state.activeDrawing());
    await assert.rejects(exportSession('invalid scalar settings'), /width|positive|scalar/i);
    assert.deepEqual(buildConfigData(state.activeDrawing()), before);
    // Session 46 Load lists a slot value Save refuses (OV-758; a value JSON
    // cannot write is the writer's check only).
    const text = JSON.stringify(scalar);
    if (text === undefined || text === 'null' || /null/.test(text)) continue;
    const invalid = structuredClone(valid);
    invalid.modes.circular.config.adv.circular_track_slots[0].width = scalar;
    const invalidBefore = structuredClone(invalid);
    /** @type {unknown} */
    let listed = null;
    const loaded = await importSession({
      target: { files: [new Blob([JSON.stringify(invalid)], { type: 'application/json' })], value: 'selected' }
    }, { confirmLoadRepairs: async (dialog) => { listed = dialog; return false; } });
    assert.equal(loaded.status, 'canceled', text);
    assert.match(JSON.stringify(listed), /circular_track_slots\[0\]\.width/, text);
    assert.deepEqual(invalid, invalidBefore);
    // A refused Load leaves the drawing as it was.
    assert.deepEqual(buildConfigData(state.activeDrawing()), before, `drawing after the refused Load of ${text}`);
  }
});
