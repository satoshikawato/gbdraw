// The question a Circular File card's file control asks before a Replace that
// removes Depth files (OV-2318), read from the dialog it opens: which records'
// Depth files leave, and whether a Depth series is then left with no file, so
// that Generate asks for a file or the series' removal (OV-3954, D-26).
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

const genbankA = readFileSync('tests/test_inputs/BGC0000708.gbk', 'utf8');
const genbankB = readFileSync('tests/test_inputs/BGC0000709.gbk', 'utf8');
const genbank = (name, ...texts) => new File([texts.join('\n')], name);

/**
 * One shown Circular drawing with one File card holding `file`, read, and the
 * record-major Depth matrix `depth` (one row per record, one column per
 * series, file names or null).
 */
const drawingWith = async (file, depth, { showDepth = true, id = 'drawing-1' } = {}) => {
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
    adv: { circular_track_slots_enabled: false, circular_track_slots: [] },
    committed: new Map([[card.uid, records.map((record) => record.record_id)]])
  };
};

/**
 * The dialog a Replace of the card file of `drawings` (the first one shown) with `value` opens, or null when it
 * opens none.
 */
const askedOnReplace = async (drawings, value) => {
  const list = [drawings].flat();
  /** @type {unknown[]} */
  const asked = [];
  const ask = createCircularSourceReplacementQuestion({
    state: { drawingOf: (id) => list.find((drawing) => drawing.id === id), drawingList: () => list, inputFiles: { fileOf: () => null } },
    askChoice: async (request) => { asked.push(request); return 'replace'; },
    showNotice: () => {},
    barTabs: () => list.map((drawing, index) => ({ kind: 'drawing', id: drawing.id, label: `Drawing ${index + 1}` })),
    committedCardRecords: async (drawing) => drawing.committed
  });
  const answer = await ask(list.map((drawing) => ({ drawingId: drawing.id, slot: 'source.gb', rowUid: 'card-1' })), value, () => true);
  assert.equal(answer, true);
  assert.ok(asked.length <= 1);
  return asked[0] ?? null;
};

const LEAVING = 'The new file has no record for these records of the File card. Replacing the file removes their Depth files.';

test('OV-3954: a Replace that leaves a Depth series with no file says Generate will ask for one', async () => {
  // One record replaced by one record with another ID: its Depth file is the series' only one.
  const single = await askedOnReplace(await drawingWith(genbank('a.gbk', genbankA), [['sample-a.tsv']]), genbank('b.gbk', genbankB));
  assert.deepEqual(single.items, ['#1 BGC0000708: sample-a.tsv']);
  assert.equal(single.description, `${LEAVING} Depth series #1 stays empty, and Generate will ask for a file or its removal.`);

  // Two series: the record that leaves holds the only file of series #2; series #1 keeps a file of the record that stays.
  const two = await askedOnReplace(
    await drawingWith(genbank('ab.gbk', genbankA, genbankB), [['sample-a.tsv', null], ['sample-b.tsv', 'sample-b2.tsv']]),
    genbank('a.gbk', genbankA)
  );
  // A record's files are one list, joined by the shown language's list pattern (`listText`).
  assert.deepEqual(two.items, ['#2 BGC0000709: sample-b.tsv and sample-b2.tsv']);
  assert.equal(two.description, `${LEAVING} Depth series #2 stays empty, and Generate will ask for a file or its removal.`);
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
  const asked = await askedOnReplace([
    await drawingWith(genbank('a.gbk', genbankA), [['sample-a.tsv']]),
    await drawingWith(genbank('a.gbk', genbankA), [['sample-a.tsv']], { id: 'drawing-2' })
  ], genbank('b.gbk', genbankB));
  assert.deepEqual(asked.items, ['#1 BGC0000708 in Drawing 1: sample-a.tsv', '#1 BGC0000708 in Drawing 2: sample-a.tsv']);
  assert.equal(asked.description,
    `${LEAVING} Depth series #1 in Drawing 1 and #1 in Drawing 2 stay empty, and Generate will ask for a file or their removal.`);
});
