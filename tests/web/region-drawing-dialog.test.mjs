// The New drawing from a region dialog's owner: open, derive (newest wins),
// Keep my values, refusal, Create as one step, Cancel records nothing.
import assert from 'node:assert/strict';
import test from 'node:test';

const { createRegionDrawing, regionSourceRecords } = await import('../../gbdraw/web/js/app/region-drawing.js');

const RECORDS = [
  { recordKey: 'record-1', caption: 'NC_1', length: 1200, features: [{ featureId: 'f1', label: 'dnaA' }] },
  { recordKey: 'record-2', caption: 'NC_2', length: 900, features: null }
];
const sent = { envelope: { renderRequest: { mode: 'circular', records: [] }, resources: {} }, filesByResource: new Map() };
const derivedReply = (adaptSizes) => ({
  status: 'derived', mode: 'linear', renderRequest: { mode: 'linear' }, resourceSources: {},
  returnedToAuto: adaptSizes ? [{ setting: 'objects.scale.interval', explicit: 25000, auto: 1000 },
    { setting: 'objects.features.block_stroke_width', explicit: 1234.5678, auto: null }] : [],
  notCarried: [{ kind: 'OTHER_MODE_SETTINGS', count: 0, items: ['canvas.circular.track_type'] }, { kind: 'COMPARISON_TABLES', count: 1, items: [] }]
});

const harness = ({ reply = null, refusal = 'ORIGIN_CROSSING' } = {}) => {
  const calls = { derive: [], created: [] };
  /** @type {Array<() => void>} */
  const pending = [];
  const owner = createRegionDrawing({
    shownDrawingId: () => 'circular',
    projectShownRequest: async () => sent,
    deriveRegion: async (payload) => {
      calls.derive.push(payload);
      if (reply === 'hold') await new Promise((resolve) => pending.push(resolve));
      const first = payload.regions[0];
      if ('start' in first && first.start > first.end) {
        return { result: { status: 'refused', reason: refusal, recordKey: first.recordKey, start: first.start, end: first.end } };
      }
      return { result: derivedReply(payload.adaptSizes) };
    },
    newDrawingFromRegion: async (sourceId, derivation) => {
      calls.created.push({ sourceId, derivation });
      return { drawingId: 'linear' };
    },
    shownRecords: async () => RECORDS
  });
  return { owner, calls, pending };
};

test('open lists the shown records by key, unchecked, and offers features only with a Result', async () => {
  const { owner } = harness();
  await owner.open();
  const { dialog } = owner;
  assert.equal(dialog.show, true);
  assert.deepEqual(dialog.rows.map((row) => [row.recordKey, row.checked, row.kind, row.start, row.end]),
    [['record-1', false, 'whole', '1', '1200'], ['record-2', false, 'whole', '1', '900']]);
  assert.equal(dialog.rows[1].features, null);
  assert.equal(dialog.status, 'idle');
});

test('a feature can be chosen only in a record whose Result drew features; each other row names why', async () => {
  const { owner } = harness();
  await owner.open();
  assert.deepEqual(owner.dialog.rows.map((row) => row.featureReason),
    ['', 'Generate the drawing to choose a feature.']);
  const ownerOfEmpty = createRegionDrawing({
    shownDrawingId: () => 'circular', projectShownRequest: async () => sent, deriveRegion: async () => ({ result: derivedReply(false) }),
    newDrawingFromRegion: async () => ({ drawingId: 'linear' }),
    shownRecords: async () => [{ recordKey: 'record-1', caption: 'NC_1', length: 1200, features: [] }]
  });
  await ownerOfEmpty.open();
  assert.deepEqual(ownerOfEmpty.dialog.rows.map((row) => row.featureReason), ['This record has no features to choose.']);
});

test('a derived reply shows the returned sizes and the not-carried lines; Keep my values derives without adapting', async () => {
  const { owner, calls } = harness();
  await owner.open();
  Object.assign(owner.dialog.rows[0], { checked: true, kind: 'coordinates', start: '101', end: '400' });
  await owner.refresh();
  assert.deepEqual(calls.derive.at(-1).regions, [{ recordKey: 'record-1', start: 101, end: 400 }]);
  assert.equal(calls.derive.at(-1).adaptSizes, true);
  assert.equal(owner.dialog.status, 'derived');
  // Values are grouped as the row's length is ("280,452 bp"), with every digit kept.
  assert.deepEqual(owner.dialog.returnedToAuto, [
    { setting: 'objects.scale.interval', label: 'Linear ruler interval', explicit: '25,000', auto: '1,000' },
    { setting: 'objects.features.block_stroke_width', label: 'Block Stroke Width', explicit: '1,234.5678', auto: 'Auto' }
  ]);
  // The other mode's settings are named by the mode they apply to (the new drawing is Linear).
  assert.deepEqual(owner.dialog.notCarried, ['Settings that apply only to Circular drawings', 'The uploaded comparison tables']);
  owner.dialog.keepValues = true;
  await owner.refresh();
  assert.equal(calls.derive.at(-1).adaptSizes, false);
  assert.deepEqual(owner.dialog.returnedToAuto, []);
});

test('a feature row sends the feature identity; an empty feature or a bad number derives nothing', async () => {
  const { owner, calls } = harness();
  await owner.open();
  Object.assign(owner.dialog.rows[0], { checked: true, kind: 'feature' });
  await owner.refresh();
  assert.equal(owner.dialog.status, 'invalid');
  owner.dialog.rows[0].featureId = 'f1';
  await owner.refresh();
  assert.deepEqual(calls.derive.at(-1).regions, [{ recordKey: 'record-1', featureId: 'f1' }]);
  Object.assign(owner.dialog.rows[1], { checked: true, kind: 'coordinates', start: '1.5', end: '10' });
  await owner.refresh();
  assert.equal(owner.dialog.status, 'invalid');
  assert.equal(calls.derive.length, 1);
});

test('the newest request wins', async () => {
  const { owner, pending } = harness({ reply: 'hold' });
  await owner.open();
  Object.assign(owner.dialog.rows[0], { checked: true, kind: 'coordinates', start: '900', end: '100' });
  const first = owner.refresh();
  Object.assign(owner.dialog.rows[0], { start: '100', end: '900' });
  const second = owner.refresh();
  while (pending.length < 2) await new Promise((resolve) => setImmediate(resolve));
  // The newer reply lands first; the older (a refusal) must not replace it.
  pending[1]();
  await second;
  pending[0]();
  await first;
  assert.equal(owner.dialog.status, 'derived');
});

test('a refusal shows its reason alone, worded for the topology, and Create does nothing', async () => {
  for (const [refusal, message] of [
    ['ORIGIN_CROSSING', 'Choose a start before the end: a region cannot cross the origin of a circular record.'],
    ['START_AFTER_END', 'Choose a start before the end.']
  ]) {
    const { owner, calls } = harness({ refusal });
    await owner.open();
    Object.assign(owner.dialog.rows[0], { checked: true, kind: 'coordinates', start: '900', end: '100' });
    await owner.refresh();
    assert.equal(owner.dialog.status, 'refused');
    assert.equal(owner.dialog.message, message, 'the reason only: no field name or invalid-region preface');
    assert.equal(await owner.confirm(), null);
    assert.equal(calls.created.length, 0);
  }
});

test('Create runs one list step and closes; Cancel records nothing', async () => {
  const { owner, calls } = harness();
  await owner.open();
  owner.dialog.rows[0].checked = true;
  await owner.refresh();
  assert.deepEqual(await owner.confirm(), { drawingId: 'linear' });
  assert.equal(calls.created.length, 1);
  assert.equal(calls.created[0].sourceId, 'circular');
  assert.equal(owner.dialog.show, false);
  await owner.open();
  owner.dialog.rows[1].checked = true;
  await owner.refresh();
  owner.close();
  assert.equal(owner.dialog.show, false);
  assert.equal(await owner.confirm(), null);
  assert.equal(calls.created.length, 1);
});

test('a close while the records are read keeps the dialog closed', async () => {
  let release = () => {};
  const owner = createRegionDrawing({
    shownDrawingId: () => 'circular',
    projectShownRequest: async () => sent,
    deriveRegion: async () => ({ result: null }),
    newDrawingFromRegion: async () => ({ drawingId: 'linear' }),
    shownRecords: () => new Promise((resolve) => { release = () => resolve(RECORDS); })
  });
  const opening = owner.open();
  assert.equal(owner.dialog.show, true);
  owner.close();
  release();
  await opening;
  assert.equal(owner.dialog.show, false);
  assert.deepEqual(owner.dialog.rows, []);
});

test('without a Result the records come from the request and discovery, with no feature choice', () => {
  // A Multi-Record Canvas of three discovered records, the third of unknown length, the second OFF.
  const sources = {
    mode: 'circular',
    circular: [
      { sourceIndex: 0, selector: '#1', recordId: 'NC_1', recordLength: 1200 },
      { sourceIndex: 1, selector: '#2', recordId: 'NC_2', recordLength: 900 },
      { sourceIndex: 2, selector: '#3', recordId: 'NC_3', recordLength: null }
    ]
  };
  const requestRecords = [
    { recordKey: 'record-1', selector: { kind: 'recordIndex', index: 0 } },
    { recordKey: 'record-3', selector: { kind: 'recordIndex', index: 2 } }
  ];
  const asked = [];
  const records = regionSourceRecords(requestRecords, sources, (recordKey) => { asked.push(recordKey); return null; });
  assert.deepEqual(records, [{ recordKey: 'record-1', caption: 'NC_1', length: 1200, features: null }]);
  assert.deepEqual(asked, ['record-1']);
});

test('every setting that can return to Auto is named by an English UI text, never its ID', async () => {
  const { adaptableSettingLabel } = await import('../../gbdraw/web/js/app/region-drawing.js');
  const { MODE_PROFILE_DATA } = await import('../../gbdraw/web/js/mode-profiles.generated.js');
  const { sourceCalls } = await import('../../tools/i18n/source-keys.mjs');
  const english = new Set(sourceCalls().map((call) => call.key));
  for (const [mode, settings] of Object.entries(MODE_PROFILE_DATA.autoSizes.adaptableSettings)) {
    for (const setting of settings) {
      const label = adaptableSettingLabel(/** @type {'circular' | 'linear'} */ (mode), setting);
      assert.ok(english.has(label), `${mode} ${setting}: "${label}" is no English UI text`);
      assert.notEqual(label, setting, `${mode} ${setting} is shown by its ID`);
    }
  }
  assert.equal(adaptableSettingLabel('linear', 'objects.unknown.size'), 'objects.unknown.size');
});

test('a setting that a section names in the pane names its track or section in the dialog list', async () => {
  const { adaptableSettingLabel } = await import('../../gbdraw/web/js/app/region-drawing.js');
  assert.deepEqual(['window', 'step', 'depth_window', 'depth_step', 'objects.gc_content.tick_font_size',
    'objects.definition.linear.font_size'].map((setting) => adaptableSettingLabel('linear', setting)), [
    'GC content/skew window', 'GC content/skew step', 'Depth Window', 'Depth Step', 'GC content tick font size',
    'Definition font size'
  ]);
});
