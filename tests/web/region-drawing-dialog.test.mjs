// The New drawing from a region dialog's owner: open, derive (newest wins),
// Keep my values, refusal, Create as one step, Cancel records nothing.
import assert from 'node:assert/strict';
import test from 'node:test';

const { createRegionDrawing, regionSourceRecords } = await import('../../gbdraw/web/js/app/region-drawing.js');
const { diagnosticError } = await import('../../gbdraw/web/js/utils/error-normalization.js');

const RECORDS = [
  { recordKey: 'record-1', caption: 'NC_1', length: 1200, features: [{ featureId: 'f1', label: 'dnaA' }] },
  { recordKey: 'record-2', caption: 'NC_2', length: 900, features: null }
];
/** The shown drawing's request, discovered records and Result features that give `records`. */
const recordPorts = (records) => {
  const sources = { mode: 'circular', circular: records.map((record, index) => (
    { sourceIndex: index, selector: `#${index + 1}`, recordId: record.caption, recordLength: record.length }
  )) };
  const requestRecords = records.map((record, index) => ({ recordKey: record.recordKey, selector: { kind: 'recordIndex', index } }));
  const features = new Map(records.map((record) => [record.recordKey, record.features]));
  const sent = { envelope: { renderRequest: { mode: 'circular', records: requestRecords }, resources: {} }, filesByResource: new Map() };
  return {
    projectShownRequest: async () => sent,
    readRecordSources: async () => sources,
    resultFeatures: (recordKey) => features.get(recordKey) ?? null
  };
};
const { projectShownRequest, readRecordSources, resultFeatures } = recordPorts(RECORDS);
const derivedReply = (adaptSizes) => ({
  status: 'derived', mode: 'linear', renderRequest: { mode: 'linear' }, resourceSources: {},
  returnedToAuto: adaptSizes ? [{ setting: 'objects.scale.interval', explicit: 25000, auto: 1000 },
    { setting: 'objects.features.block_stroke_width', explicit: 1234.5678, auto: null }] : [],
  notCarried: [{ kind: 'OTHER_MODE_SETTINGS', count: 0, items: ['canvas.circular.track_type'] }, { kind: 'COMPARISON_TABLES', count: 1, items: [] }]
});

const harness = ({ reply = null, refusal = 'START_AFTER_END' } = {}) => {
  const calls = { derive: [], created: [] };
  /** @type {Array<() => void>} */
  const pending = [];
  const owner = createRegionDrawing({
    shownDrawingId: () => 'circular',
    projectShownRequest,
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
    readRecordSources,
    resultFeatures
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
    shownDrawingId: () => 'circular', deriveRegion: async () => ({ result: derivedReply(false) }),
    newDrawingFromRegion: async () => ({ drawingId: 'linear' }),
    ...recordPorts([{ recordKey: 'record-1', caption: 'NC_1', length: 1200, features: [] }])
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

test('a refusal shows its reason alone and Create does nothing', async () => {
  for (const [refusal, message] of [
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
    projectShownRequest,
    deriveRegion: async () => ({ result: null }),
    newDrawingFromRegion: async () => ({ drawingId: 'linear' }),
    readRecordSources: () => new Promise((resolve) => { release = () => resolve(readRecordSources()); }),
    resultFeatures
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
  assert.equal(asked[0], 'record-1');
});

// The shown Result keyed its records as the committed request did; the rows
// are keyed as the draft is. The join is by the record's identity: its
// request key, then the keys it has under the other canvas (aliases).
test('a record drawn alone after a Multi-Record Canvas Generate offers the features the Result drew of it', async () => {
  const { circularSingleRecordKey } = await import('../../gbdraw/web/js/services/record-draw-selection.js');
  const entry = { sourceIndex: 1, selector: '#2', recordId: 'NC_2', recordLength: 900 };
  const sources = { mode: 'circular', circular: [{ sourceIndex: 0, selector: '#1', recordId: 'NC_1', recordLength: 1200 }, entry] };
  const draftKey = circularSingleRecordKey(entry);
  assert.notEqual(draftKey, 'record-2');
  // The Generated Multi-Record Canvas drew both records as record-1 and record-2.
  const catalog = new Map([['record-1', [{ featureId: 'f1', label: 'dnaA' }]], ['record-2', [{ featureId: 'f2', label: 'recA' }]]]);
  const records = regionSourceRecords([{ recordKey: draftKey, selector: { kind: 'recordIndex', index: 1 } }], sources,
    (recordKey) => catalog.get(recordKey) ?? null);
  assert.deepEqual(records, [{ recordKey: draftKey, caption: 'NC_2', length: 900, features: [{ featureId: 'f2', label: 'recA' }] }]);
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

// The read of a deferred list empties it and fills it when the discovery
// returns; the rows come from the request projected after that read, so every
// record of a Multi-Record Canvas is offered.
test('the rows come from the request projected after the records are read', async () => {
  const discovered = [{ sourceIndex: 0, selector: '#1', recordId: 'NC_1', recordLength: 1200 },
    { sourceIndex: 1, selector: '#2', recordId: 'NC_2', recordLength: 900 }];
  const list = [];
  const owner = createRegionDrawing({
    shownDrawingId: () => 'circular',
    // The projection keys the records of the list as it is when projected.
    projectShownRequest: async () => ({ envelope: { renderRequest: { mode: 'circular', records: list.map((entry) => (
      { recordKey: `record-${entry.sourceIndex + 1}`, selector: { kind: 'recordIndex', index: entry.sourceIndex } }
    )) }, resources: {} }, filesByResource: new Map() }),
    deriveRegion: async () => ({ result: null }),
    newDrawingFromRegion: async () => ({ drawingId: 'linear' }),
    readRecordSources: async () => {
      list.splice(0);
      await new Promise((resolve) => setTimeout(resolve, 5));
      list.push(...discovered);
      return { mode: 'circular', circular: [...list] };
    },
    resultFeatures: () => null
  });
  await owner.open();
  assert.deepEqual(owner.dialog.rows.map((row) => [row.recordKey, row.caption, row.length]),
    [['record-1', 'NC_1', 1200], ['record-2', 'NC_2', 900]]);
});

// A second Create while the first one runs adds nothing: the dialog is
// 'creating' and Create is disabled until the drawing is added.
test('a second Create while one runs adds one drawing', async () => {
  const { owner, calls } = harness();
  await owner.open();
  owner.dialog.rows[0].checked = true;
  await owner.refresh();
  const first = owner.confirm();
  assert.equal(owner.dialog.status, 'creating', 'Create is disabled while the drawing is added');
  assert.equal(owner.dialog.returnedToAuto.length, 2, 'the derived view stays shown');
  const second = owner.confirm();
  assert.deepEqual(await Promise.all([first, second]), [{ drawingId: 'linear' }, null]);
  assert.equal(calls.created.length, 1);
  assert.equal(owner.dialog.show, false);
});

// A Create whose prepare fails shows the reason in the dialog; the dialog
// stays open and the failure reaches no unhandled rejection.
test('a Create that fails shows its reason and stays open', async () => {
  const owner = createRegionDrawing({
    shownDrawingId: () => 'circular',
    deriveRegion: async () => ({ result: derivedReply(true) }),
    newDrawingFromRegion: async () => { throw diagnosticError('RESULT_INVALID'); },
    ...recordPorts(RECORDS)
  });
  await owner.open();
  owner.dialog.rows[0].checked = true;
  await owner.refresh();
  assert.equal(await owner.confirm(), null);
  assert.equal(owner.dialog.show, true);
  assert.equal(owner.dialog.status, 'error');
  assert.notEqual(owner.dialog.message, '');
  assert.equal(await owner.confirm(), null, 'Create stays disabled until the region derives again');
});

// The owner holds every dialog transition: the template calls editRow and
// setKeepValues and assigns no dialog state.
test('an edit through the owner writes the row found by its record key and derives again', async () => {
  const { owner, calls } = harness();
  await owner.open();
  await owner.editRow('record-2', { checked: true, kind: 'coordinates', start: '11', end: '20' });
  assert.deepEqual(calls.derive.at(-1).regions, [{ recordKey: 'record-2', start: 11, end: 20 }]);
  assert.deepEqual([owner.dialog.rows[0].checked, owner.dialog.rows[1].start], [false, '11']);
  assert.equal(owner.dialog.status, 'derived');
  await owner.setKeepValues(true);
  assert.equal(owner.dialog.keepValues, true);
  assert.equal(calls.derive.at(-1).adaptSizes, false);
});

test('an edit of an unknown record key changes nothing and derives nothing', async () => {
  const { owner, calls } = harness();
  await owner.open();
  await owner.editRow('record-2', { checked: true });
  const before = JSON.stringify(owner.dialog);
  const derives = calls.derive.length;
  await owner.editRow('record-9', { checked: false, start: '5' });
  assert.equal(JSON.stringify(owner.dialog), before);
  assert.equal(calls.derive.length, derives);
});

// The reply's vocabularies come from Python (mode-profiles.generated.js).
test('every kind of a part not carried that Python names has a text in each mode', async () => {
  const { notCarriedLine } = await import('../../gbdraw/web/js/app/region-drawing.js');
  const { MODE_PROFILE_DATA } = await import('../../gbdraw/web/js/mode-profiles.generated.js');
  const { sourceCalls } = await import('../../tools/i18n/source-keys.mjs');
  const english = new Set(sourceCalls().map((call) => call.key));
  assert.ok(MODE_PROFILE_DATA.regionDrawing.notCarriedKinds.length > 0);
  for (const kind of MODE_PROFILE_DATA.regionDrawing.notCarriedKinds) {
    for (const mode of ['circular', 'linear']) {
      const line = notCarriedLine(mode, kind);
      assert.ok(english.has(line), `${mode} ${kind}: "${line}" is no English UI text`);
    }
  }
});

test('a refusal is read for each reason Python names and for no other', async () => {
  for (const refusal of [...(await import('../../gbdraw/web/js/mode-profiles.generated.js')).MODE_PROFILE_DATA.regionDrawing.refusalReasons, 'NOT_A_REASON']) {
    const { owner } = harness({ refusal });
    await owner.open();
    await owner.editRow('record-1', { checked: true, kind: 'coordinates', start: '900', end: '100' });
    assert.equal(owner.dialog.status, refusal === 'NOT_A_REASON' ? 'error' : 'refused', refusal);
  }
});

// The dialog names a setting as the settings pane's control names it: each
// setting ID is matched to the form field the pane control writes (never to
// a caption), and the control's label key is read from index.html.
const PANE_FIELD = Object.freeze({
  'objects.features.block_stroke_width': 'adv.block_stroke_width',
  'objects.features.line_stroke_width': 'adv.line_stroke_width',
  'objects.legends.color_rect_size': 'adv.legend_box_size',
  'objects.legends.font_size': 'adv.legend_font_size',
  'labels.font_size': 'adv.label_font_size',
  'labels.font_size.linear': 'adv.label_font_size',
  'objects.axis.circular.stroke_width': 'adv.axis_stroke_width',
  'objects.axis.linear.stroke_width': 'adv.axis_stroke_width',
  'canvas.linear.default_cds_height': 'adv.feature_height',
  'objects.scale.font_size': 'adv.scale_font_size',
  'objects.scale.ruler_label_font_size': 'adv.ruler_label_font_size',
  depth_window: 'adv.depth_window_size',
  depth_step: 'adv.depth_step_size',
  'circularSlot.features.width': 'adv.feature_width_circular',
  'circularSlot.gc_content.width': 'adv.gc_content_width_circular',
  'circularSlot.gc_skew.width': 'adv.gc_skew_width_circular',
  'circularSlot.depth.width': 'adv.depth_width_circular',
  'objects.scale.interval': 'adv.scale_interval',
  'canvas.linear.default_gc_height': 'adv.gc_height',
  'canvas.linear.depth_height': 'adv.depth_height',
  'objects.depth.tick_font_size': 'track.config.tick_font_size'
});
// A pane control whose label is read with its section heading; the dialog's
// flat list names the section too.
const SECTION_NAMED = Object.freeze({
  window: ['adv.window_size', 'Window'],
  step: ['adv.step_size', 'Step'],
  'objects.gc_content.tick_font_size': ['adv.gc_content_tick_font_size', 'Tick Font Size'],
  'objects.definition.linear.font_size': ['adv.def_font_size', 'Default font size']
});
// No settings-pane control writes these.
const NO_PANE_CONTROL = Object.freeze(['canvas.linear.arrow_length_parameter']);

/** The label keys of the pane controls that write `field` (v-model or :value). */
const paneLabelKeys = async (field) => {
  const { readFileSync } = await import('node:fs');
  const html = readFileSync(new URL('../../gbdraw/web/index.html', import.meta.url), 'utf8');
  const escaped = field.replace(/[.*+?^${}()|[\]\\]/g, '\\$&');
  const tag = new RegExp(`<(?:input|select)\\b(?:[^>"]|"[^"]*")*?(?:v-model(?:\\.number)?="${escaped}"|:value="[^"]*\\b${escaped}\\b[^"]*")(?:[^>"]|"[^"]*")*>`, 'g');
  return [...html.matchAll(tag)].map((match) => /:aria-label="[^"]*?t\('([^']*)'/.exec(match[0])?.[1] ?? null);
};

test('each setting is named by the label of the pane control that writes it', async () => {
  const { adaptableSettingLabel } = await import('../../gbdraw/web/js/app/region-drawing.js');
  const { MODE_PROFILE_DATA } = await import('../../gbdraw/web/js/mode-profiles.generated.js');
  for (const [mode, settings] of Object.entries(MODE_PROFILE_DATA.autoSizes.adaptableSettings)) {
    for (const setting of settings) {
      const label = adaptableSettingLabel(/** @type {'circular' | 'linear'} */ (mode), setting);
      if (Object.hasOwn(PANE_FIELD, setting)) {
        const keys = await paneLabelKeys(PANE_FIELD[setting]);
        assert.ok(keys.length > 0, `${setting}: no pane control writes ${PANE_FIELD[setting]}`);
        assert.ok(keys.includes(label), `${mode} ${setting}: "${label}" is not the pane label (${keys.join(', ')})`);
      } else if (Object.hasOwn(SECTION_NAMED, setting)) {
        const [field, paneKey] = SECTION_NAMED[setting];
        assert.deepEqual([...new Set(await paneLabelKeys(field))], [paneKey], `${setting}: the pane label changed`);
        assert.notEqual(label, paneKey, `${setting}: the dialog needs the section in its name`);
      } else {
        assert.ok(NO_PANE_CONTROL.includes(setting), `${mode} ${setting} is in no group of this test`);
        assert.deepEqual(await paneLabelKeys(`adv.${setting.split('.').at(-1)}`), []);
      }
    }
  }
});
