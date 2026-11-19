// The question a Circular File card's file control asks before a Replace that
// removes Depth files (OV-2318), read from the dialog it opens: which records'
// Depth files leave, and what that does as Remove Depth file does it (OV-3955):
// a Depth series left with no file while others keep one makes Generate ask for
// a file or the series' removal (OV-3954, D-26); the last Depth file leaving
// turns Show Depth off, so nothing is asked; a row that reaches a row that
// leaves turns off too. The answer runs that Depth source change after the
// card files are written, in the same step.
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { test } from 'node:test';

globalThis.window = {
  Vue: {
    ref: (value) => ({ value }), reactive: (value) => value,
    computed: (getter) => ({ get value() { return getter(); } }), nextTick: async () => {}
  },
  DOMPurify: { sanitize: (value) => value }
};
globalThis.document = { body: { appendChild: () => {} }, createElement: () => ({ addEventListener: () => {} }) };

const { createCircularSourceReplacementQuestion } = await import('../../gbdraw/web/js/app/circular-record-controls.js');
const { readCircularSourceRecords } = await import('../../gbdraw/web/js/services/record-discovery.js');
const { circularSourceInputs } = await import('../../gbdraw/web/js/services/input-types.js');
const { activeDepthTrackIndices, depthSourceLoss } = await import('../../gbdraw/web/js/services/depth-track-state.js');
const { draftDrawsDepth } = await import('../../gbdraw/web/js/services/track-legend-rows.js');

const genbankA = readFileSync('tests/test_inputs/BGC0000708.gbk', 'utf8');
const genbankB = readFileSync('tests/test_inputs/BGC0000709.gbk', 'utf8');
const genbank = (name, ...texts) => new File([texts.join('\n')], name);

/**
 * One shown Circular drawing with one File card holding `file`, read, and the
 * record-major Depth matrix `depth` (one row per record, one column per
 * series, file names or null).
 */
const drawingWith = async (file, depth, { showDepth = true, id = 'drawing-1', slots = [] } = {}) => {
  const card = { uid: 'card-1', gb: file, gff: null, fasta: null };
  const inputs = circularSourceInputs('gb', [card]);
  const [{ records }] = await readCircularSourceRecords(inputs, (value) => value);
  return {
    id,
    inputType: { value: 'gb' },
    circularSources: [card],
    circularRecordDiscovery: {
      inputType: 'gb', status: 'ready', error: '',
      sources: inputs.map(({ uid, primaryFile, pairedFile }) => ({ uid, primaryFile, pairedFile, error: '' }))
    },
    circularRecordList: { value: records },
    files: { c_depth: depth.map((row) => row.map((name) => (name ? { name } : null))), c_depth_rows: null },
    form: { show_depth: showDepth },
    adv: { circular_track_slots_enabled: slots.length > 0, circular_track_slots: slots },
    committed: new Map([[card.uid, records.map((record) => record.record_id)]])
  };
};

/**
 * The dialog a Replace of the card file of `drawings` (the first one shown) with `value` opens, or null when it
 * opens none, with the Depth source changes its answer runs (`changes`: the drawing and the change of each).
 * The stack's loss is the track editor's (`circularDepthSourceLoss`: the series with a file before and after).
 */
const askedOnReplace = async (drawings, value) => {
  const list = [drawings].flat();
  /** @type {unknown[]} */
  const asked = [];
  const changes = [];
  const ask = createCircularSourceReplacementQuestion({
    state: { drawingOf: (id) => list.find((drawing) => drawing.id === id), drawingList: () => list, inputFiles: { fileOf: () => null } },
    askChoice: async (request) => { asked.push(request); return 'replace'; },
    showNotice: () => {},
    barTabs: () => list.map((drawing, index) => ({ kind: 'drawing', id: drawing.id, label: `Drawing ${index + 1}` })),
    committedCardRecords: async (drawing) => drawing.committed,
    circularDepthSourceLoss: (drawing, next) => depthSourceLoss(drawing.adv.circular_track_slots,
      activeDepthTrackIndices(drawing.files.c_depth), activeDepthTrackIndices(next),
      (rows) => draftDrawsDepth('circular', { form: drawing.form, adv: { ...drawing.adv, circular_track_slots: rows } })),
    changeCircularDepthSources: (drawing, _mutate, change) => { changes.push({ drawing: drawing.id, change }); }
  });
  const answer = await ask(list.map((drawing) => ({ drawingId: drawing.id, slot: 'source.gb', rowUid: 'card-1' })), value, () => true);
  assert.equal(typeof answer, 'function', 'the answer writes after the card files');
  assert.ok(asked.length <= 1);
  answer();
  return asked[0] ? { ...asked[0], changes } : null;
};
const SAME_SERIES = { setsFile: false, series: 'same' };

const LEAVING = 'The new file has no record for these records of the File card. Replacing the file removes their Depth files.';

test('OV-3955: a Replace that removes the last Depth file runs the Depth source change and asks for no file', async () => {
  // One record replaced by one record with another ID: its Depth file is the drawing's only one, so Show Depth
  // turns off (Remove Depth file's change) and Generate draws without Depth.
  const single = await askedOnReplace(await drawingWith(genbank('a.gbk', genbankA), [['sample-a.tsv']]), genbank('b.gbk', genbankB));
  assert.deepEqual(single.items, ['#1 BGC0000708: sample-a.tsv']);
  assert.equal(single.description, LEAVING);
  assert.deepEqual(single.changes, [{ drawing: 'drawing-1', change: SAME_SERIES }]);
});

test('OV-3954: a Replace that leaves a Depth series with no file says Generate will ask for one', async () => {
  // Two series: the record that leaves holds the only file of series #2; series #1 keeps a file of the record that stays.
  const two = await askedOnReplace(
    await drawingWith(genbank('ab.gbk', genbankA, genbankB), [['sample-a.tsv', null], ['sample-b.tsv', 'sample-b2.tsv']]),
    genbank('a.gbk', genbankA)
  );
  // A record's files are one list, joined by the shown language's list pattern (`listText`).
  assert.deepEqual(two.items, ['#2 BGC0000709: sample-b.tsv and sample-b2.tsv']);
  assert.equal(two.description, `${LEAVING} Depth series #2 stays empty, and Generate will ask for a file or its removal.`);
  assert.deepEqual(two.changes, [{ drawing: 'drawing-1', change: SAME_SERIES }]);
});

test('OV-3955: a Replace that removes a Depth row another row reaches names that row, which turns off too', async () => {
  const row = (id, renderer, params = {}, side = 'inside') => ({ id, renderer, enabled: true, side, z: 0, params });
  const slots = [
    row('features', 'features', { lane_direction: 'outside' }, 'outside'),
    row('depth', 'depth', { track_index: 0 }),
    { ...row('mark', 'annotations', { set_id: 'tint_set', anchor_slot: 'depth', layer: 'foreground' }, 'overlay'), z: 1 }
  ];
  const asked = await askedOnReplace(await drawingWith(genbank('a.gbk', genbankA), [['sample-a.tsv']], { slots }), genbank('b.gbk', genbankB));
  assert.deepEqual(asked.items, ['#1 BGC0000708: sample-a.tsv', 'Track row 3 (mark)']);
  assert.equal(asked.description, `${LEAVING} This row starts at or extends to a track that this change removes or turns off. It is turned off too, as it cannot be drawn without that track.`);
});

test('OV-3954: a Replace that leaves every Depth series a file, or a drawing that draws no Depth, adds no sentence', async () => {
  const kept = await askedOnReplace(
    await drawingWith(genbank('ab.gbk', genbankA, genbankB), [['sample-a.tsv'], ['sample-b.tsv']]),
    genbank('a.gbk', genbankA)
  );
  assert.deepEqual(kept.items, ['#2 BGC0000709: sample-b.tsv']);
  assert.equal(kept.description, LEAVING);

  // Show Depth off: Generate draws no Depth, so it asks for no file.
  const hidden = await askedOnReplace(
    await drawingWith(genbank('a.gbk', genbankA), [['sample-a.tsv']], { showDepth: false }),
    genbank('b.gbk', genbankB)
  );
  assert.equal(hidden.description, LEAVING);
});

test('OV-3954: a Replace in several drawings names each series left with no file with its drawing', async () => {
  const depth = [['sample-a.tsv', null], ['sample-b.tsv', 'sample-b2.tsv']];
  const asked = await askedOnReplace([
    await drawingWith(genbank('ab.gbk', genbankA, genbankB), depth),
    await drawingWith(genbank('ab.gbk', genbankA, genbankB), depth, { id: 'drawing-2' })
  ], genbank('a.gbk', genbankA));
  assert.deepEqual(asked.items, ['#2 BGC0000709 in Drawing 1: sample-b.tsv and sample-b2.tsv', '#2 BGC0000709 in Drawing 2: sample-b.tsv and sample-b2.tsv']);
  assert.equal(asked.description,
    `${LEAVING} Depth series #2 in Drawing 1 and #2 in Drawing 2 stay empty, and Generate will ask for a file or their removal.`);
  assert.deepEqual(asked.changes, [{ drawing: 'drawing-1', change: SAME_SERIES }, { drawing: 'drawing-2', change: SAME_SERIES }]);
});
