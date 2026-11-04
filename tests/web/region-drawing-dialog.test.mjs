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
  returnedToAuto: adaptSizes ? [{ setting: 'objects.scale.interval', explicit: 500, auto: 100 }] : [],
  notCarried: [{ kind: 'COMPARISON_TABLES', count: 1, items: [] }]
});

const harness = ({ reply = null } = {}) => {
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
        return { result: { status: 'refused', reason: 'ORIGIN_CROSSING', recordKey: first.recordKey, start: first.start, end: first.end } };
      }
      return { result: derivedReply(payload.adaptSizes) };
    },
    newDrawingFromRegion: async (sourceId, derivation) => {
      calls.created.push({ sourceId, derivation });
      return { drawingId: 'linear' };
    },
    shownRecords: async () => RECORDS,
    settingLabel: (setting) => (setting === 'objects.scale.interval' ? 'Tick interval' : setting)
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

test('a derived reply shows the returned sizes and the not-carried lines; Keep my values derives without adapting', async () => {
  const { owner, calls } = harness();
  await owner.open();
  Object.assign(owner.dialog.rows[0], { checked: true, kind: 'coordinates', start: '101', end: '400' });
  await owner.refresh();
  assert.deepEqual(calls.derive.at(-1).regions, [{ recordKey: 'record-1', start: 101, end: 400 }]);
  assert.equal(calls.derive.at(-1).adaptSizes, true);
  assert.equal(owner.dialog.status, 'derived');
  assert.deepEqual(owner.dialog.returnedToAuto, [{ setting: 'objects.scale.interval', label: 'Tick interval', explicit: '500', auto: '100' }]);
  assert.deepEqual(owner.dialog.notCarried, ['The uploaded comparison tables']);
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

test('a refusal shows the reason and Create does nothing', async () => {
  const { owner, calls } = harness();
  await owner.open();
  Object.assign(owner.dialog.rows[0], { checked: true, kind: 'coordinates', start: '900', end: '100' });
  await owner.refresh();
  assert.equal(owner.dialog.status, 'refused');
  assert.match(owner.dialog.message, /cross the origin/);
  assert.equal(await owner.confirm(), null);
  assert.equal(calls.created.length, 0);
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
    shownRecords: () => new Promise((resolve) => { release = () => resolve(RECORDS); }),
    settingLabel: (setting) => setting
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
