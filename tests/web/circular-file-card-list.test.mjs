// The Circular File card list (app/circular-file-cards.js, circular-multi-file
// D-01, D-10): one card per source file in drawing order, each card's header
// line, state, records line and fold; Add files and the File removal dialog
// call the source card writers (fakes here); a record's Topology card names its
// file. One pick or drop of GFF3 and FASTA files pairs them by file-name stem
// (`pickedFileCards`, circular-multi-file D-07). The shell's header line is the one the Linear File cards show
// (app/file-card.js), and the Choose records list of several cards groups its
// rows by file (app/record-selection.js). tests/web/circular-file-card-list.playwright.spec.js
// runs the list in the browser.
import assert from 'node:assert/strict';
import test from 'node:test';
import { loadPageVue } from './helpers/vue-runtime.mjs';

const Vue = loadPageVue();
globalThis.window = { Vue };
const { reactive } = Vue;
const { createCircularFileCards } = await import('../../gbdraw/web/js/app/circular-file-cards.js');
const { fileCardLine } = await import('../../gbdraw/web/js/app/file-card.js');
const { pickedFileCards } = await import('../../gbdraw/web/js/services/input-types.js');
const { createRecordSelection } = await import('../../gbdraw/web/js/app/record-selection.js');

const file = (name) => ({ fileId: `id-${name}`, name });
const record = (uid, k, recordId, length, identity = uid) => ({
  uid, k, identity, selector: `#${k}`, record_id: recordId, record_length: length, detectedTopology: 'circular'
});

const makeCards = ({ inputType = 'gb', sources, records = [], status = 'ready', failedSources = [], error = '' }) => {
  const drawing = reactive({
    inputType: { value: inputType },
    circularSources: sources,
    circularRecordList: { value: records }
  });
  const discovery = reactive({ status, failedSources, error });
  const calls = [];
  const cards = createCircularFileCards({
    shownDrawingOfMode: (mode) => (mode === 'circular' ? drawing : null),
    circularRecordDiscoveryState: { get value() { return discovery; } },
    addCircularSourceFiles: async (files) => {
      calls.push(['add', files.map((entry) => entry.name)]);
      pickedFileCards(drawing.inputType.value, files).forEach((slots, index) => (
        drawing.circularSources.push({ uid: `new-${index}`, gb: null, gff: null, fasta: null, ...slots })));
    },
    removeCircularSourceCard: async (uid) => {
      calls.push(['remove', uid]);
      drawing.circularSources.splice(drawing.circularSources.findIndex((source) => source.uid === uid), 1);
    }
  });
  return { cards, drawing, discovery, calls };
};

const SAKAI = () => [
  { uid: 'record-1', gb: file('BA000007.3.gb'), gff: null, fasta: null },
  { uid: 'record-2', gb: file('AB011549.2.gb'), gff: null, fasta: null },
  { uid: 'record-3', gb: file('AB011548.2.gb'), gff: null, fasta: null }
];
const SAKAI_RECORDS = [
  record('record-1', 1, 'BA000007.3', 5498578),
  record('record-2', 1, 'AB011549.2', 92721),
  record('record-3', 1, 'AB011548.2', 3306)
];

test('the header line names the files of the slots the input type reads, as Linear says it', () => {
  assert.equal(fileCardLine('gb', { gb: file('a.gb'), gff: null, fasta: null }), 'a.gb');
  assert.equal(fileCardLine('gff', { gb: null, gff: file('a.gff3'), fasta: null }), 'a.gff3 + no FASTA');
  assert.equal(fileCardLine('gff', { gb: null, gff: file('a.gff3'), fasta: file('a.fasta') }), 'a.gff3 + a.fasta');
  assert.equal(fileCardLine('fasta', { gb: file('a.gb'), gff: null, fasta: null }), 'No file');
});

test('one card per file in drawing order, each with its own records line, folded at first', () => {
  const { cards } = makeCards({ sources: SAKAI(), records: SAKAI_RECORDS });
  assert.deepEqual(cards.cards.value.map((card) => [card.index, card.uid, card.line, card.ready, card.failed, card.recordsLine]), [
    [0, 'record-1', 'BA000007.3.gb', true, false, '1 record · 5,498,578 bp'],
    [1, 'record-2', 'AB011549.2.gb', true, false, '1 record · 92,721 bp'],
    [2, 'record-3', 'AB011548.2.gb', true, false, '1 record · 3,306 bp']
  ]);
  assert.deepEqual(cards.cards.value.map((card) => cards.isFolded(card.uid)), [true, true, true]);
  cards.setFolded('record-2', false);
  assert.deepEqual(cards.cards.value.map((card) => cards.isFolded(card.uid)), [true, false, true]);
});

test('a card whose read failed shows its warning and no records; the other cards still list', () => {
  const { cards } = makeCards({ sources: SAKAI(), records: SAKAI_RECORDS, status: 'error', failedSources: ['record-2'] });
  assert.deepEqual(cards.cards.value.map((card) => [card.failed, card.recordsLine]), [
    [false, '1 record · 5,498,578 bp'], [true, ''], [false, '1 record · 3,306 bp']
  ]);
});

test('the read error names the failed card by its files when the list holds several cards', () => {
  const error = { summaryParts: [{ text: 'No records were found.' }, { text: 'This file looks like FASTA.' }] };
  const several = makeCards({ sources: SAKAI(), records: SAKAI_RECORDS, status: 'error', failedSources: ['record-2'], error });
  assert.equal(several.cards.readError.value, 'AB011549.2.gb: No records were found. This file looks like FASTA.');
  const one = makeCards({ sources: SAKAI().slice(1, 2), status: 'error', failedSources: ['record-2'], error });
  assert.equal(one.cards.readError.value, 'No records were found. This file looks like FASTA.');
  const unknown = makeCards({ sources: SAKAI().slice(1, 2), status: 'error', failedSources: ['record-2'] });
  assert.equal(unknown.cards.readError.value, 'Source inspection failed. Retry or replace the input.');
  assert.equal(makeCards({ sources: SAKAI(), records: SAKAI_RECORDS }).cards.readError.value, '');
});

test('records not read yet give no records line', () => {
  const { cards } = makeCards({ sources: SAKAI(), records: SAKAI_RECORDS, status: 'deferred' });
  assert.deepEqual(cards.cards.value.map((card) => card.recordsLine), ['', '', '']);
});

test('a file of several records: the records line counts them and a Topology card names its #k', () => {
  const { cards } = makeCards({
    sources: [{ uid: 'record-1', gb: file('two.gb'), gff: null, fasta: null }, { uid: 'record-2', gb: file('mito.gb'), gff: null, fasta: null }],
    records: [record('record-1', 1, 'A', 100, 'record-1:1'), record('record-1', 2, 'B', 50, 'record-1:2'), record('record-2', 1, 'M', 16569)]
  });
  assert.equal(cards.cards.value[0].recordsLine, '2 records · 150 bp');
  assert.equal(cards.recordPlace({ sourceUid: 'record-1', selector: '#2' }), 'two.gb #2');
  assert.equal(cards.recordPlace({ sourceUid: 'record-2', selector: '#1' }), 'mito.gb');
});

test('one pick or drop of GenBank or FASTA files is one card per file, in the order given, whatever the names say', () => {
  const names = (cards) => cards.map((slots) => Object.fromEntries(Object.entries(slots).map(([slot, entry]) => [slot, entry.name])));
  const picked = [file('BA000007.3.gb'), file('AB011549.2.gb'), file('pOSAK1.fasta')];
  assert.deepEqual(names(pickedFileCards('gb', picked)), [{ gb: 'BA000007.3.gb' }, { gb: 'AB011549.2.gb' }, { gb: 'pOSAK1.fasta' }]);
  assert.deepEqual(names(pickedFileCards('fasta', picked)), [{ fasta: 'BA000007.3.gb' }, { fasta: 'AB011549.2.gb' }, { fasta: 'pOSAK1.fasta' }]);
  assert.deepEqual(pickedFileCards('gb', []), []);
});

test('GFF3 and FASTA files of one pick pair by file-name stem; an unpaired file makes a card with the other slot empty', () => {
  const names = (cards) => cards.map((slots) => [slots.gff?.name ?? null, slots.fasta?.name ?? null]);
  // The picker's order: the card stands where its first file stood.
  assert.deepEqual(names(pickedFileCards('gff', [
    file('NC_013668.fasta'), file('NC_013668.gff3'), file('NC_013668-reannotated.gff3'), file('pX.fa')
  ])), [['NC_013668.gff3', 'NC_013668.fasta'], ['NC_013668-reannotated.gff3', null], [null, 'pX.fa']]);
  // The stem ignores letter case and a .gz ending; a name that is not .gff or .gff3 is the FASTA.
  assert.deepEqual(names(pickedFileCards('gff', [file('Chr1.GFF.gz'), file('chr1.fna.gz'), file('chr2.txt'), file('chr2.gff3')])),
    [['Chr1.GFF.gz', 'chr1.fna.gz'], ['chr2.gff3', 'chr2.txt']]);
  // Two GFF3 files of one stem: the first takes the FASTA, the second is a card of its own.
  assert.deepEqual(names(pickedFileCards('gff', [file('a.gff'), file('a.gff3'), file('a.fasta')])),
    [['a.gff', 'a.fasta'], ['a.gff3', null]]);
});

test('Add files: one call of the writer for every file of a pick; a card still waiting for its FASTA shows open', async () => {
  const gb = makeCards({ sources: SAKAI().slice(0, 1) });
  await gb.cards.addFiles([file('AB011549.2.gb'), file('AB011548.2.gb')]);
  assert.deepEqual(gb.calls, [['add', ['AB011549.2.gb', 'AB011548.2.gb']]]);
  assert.deepEqual(gb.cards.cards.value.map((card) => card.line), ['BA000007.3.gb', 'AB011549.2.gb', 'AB011548.2.gb']);
  assert.equal(gb.cards.isFolded('new-0'), true);

  const gff = makeCards({ inputType: 'gff', sources: [] });
  await gff.cards.addFiles([file('a.fasta'), file('a.gff3'), file('b.gff3')]);
  assert.deepEqual(gff.cards.cards.value.map((card) => [card.line, card.ready, gff.cards.isFolded(card.uid)]), [
    ['a.gff3 + a.fasta', true, true], ['b.gff3 + no FASTA', false, false]
  ]);
  await gff.cards.addFiles([]);
  assert.equal(gff.calls.length, 1);
});

test('the Add files zone names, accepts and offers the files of the input type', () => {
  const zoneOf = (inputType) => {
    const { cards } = makeCards({ inputType, sources: [] });
    const { label, fileKind, accept } = cards.addZone.value;
    return [label, fileKind, accept.split(',').filter((entry) => entry.startsWith('.')).join(' ')];
  };
  assert.deepEqual(zoneOf('gb'), ['GenBank/DDBJ files', 'genbank', '.gb .gbk .gbff .genbank .ddbj .embl .emb .txt .gz']);
  assert.deepEqual(zoneOf('gff'), ['GFF3 and FASTA files', 'gff', '.gff .gff3 .fa .fas .fasta .fna .ffn .faa .txt .gz']);
  assert.deepEqual(zoneOf('fasta'), ['FASTA files', 'fasta', '.fa .fas .fasta .fna .ffn .txt .gz']);
  // Once a card shows, the zone offers no Session file: the cards' own uploaders do.
  assert.equal(makeCards({ sources: SAKAI() }).cards.addZone.value.fileKind, undefined);
});

test('Remove asks first: the dialog names the card and its records, Cancel removes nothing, Remove file removes that card', async () => {
  const { cards, calls } = makeCards({ sources: SAKAI(), records: SAKAI_RECORDS });
  cards.requestRemoval('record-2');
  assert.equal(cards.removal.open, true);
  assert.equal(cards.removalTarget.value.line, 'AB011549.2.gb');
  assert.equal(cards.removalTarget.value.records.length, 1);
  cards.cancelRemoval();
  assert.equal(cards.removal.open, false);
  assert.deepEqual(calls, []);
  cards.requestRemoval('record-2');
  await cards.applyRemoval();
  assert.deepEqual(calls, [['remove', 'record-2']]);
  assert.equal(cards.removal.open, false);
  assert.deepEqual(cards.cards.value.map((card) => card.uid), ['record-1', 'record-3']);
});

test('the Choose records list of several cards groups its rows by file, sorted and searched within each file', async () => {
  const source = {
    key: 'circular', name: 'a.gb + b.gb',
    records: [
      { key: 'u1:1', recordId: 'chr', length: 300 }, { key: 'u1:2', recordId: 'plasmid', length: 20 },
      { key: 'u2', recordId: 'phage', length: 50 }
    ],
    groups: [{ key: 'u1', name: 'a.gb', keys: ['u1:1', 'u1:2'] }, { key: 'u2', name: 'b.gb', keys: ['u2'] }]
  };
  const owner = createRecordSelection({
    reactive, computed: Vue.computed,
    recordsOff: () => [], source: (mode, key) => (mode === 'circular' && key === 'circular' ? source : null),
    runUndoable: async (_label, change) => change(), withDialogChoice: (_label, handler) => handler,
    closeAfterDialogChoice: (close) => close(), afterRecordSetChange: () => {}, removeSource: () => true,
    askRemoveSource: () => null, deleteRecordSettings: () => true, prepareRecords: async () => {},
    autoOpenRequests: () => []
  });
  await owner.openRecordList('circular', 'circular');
  const view = () => owner.recordListView.value.groups.map((group) => [group.name, group.rows.map((row) => row.recordId)]);
  assert.deepEqual(view(), [['a.gb', ['chr', 'plasmid']], ['b.gb', ['phage']]]);
  owner.setRecordListSort('length-asc');
  assert.deepEqual(view(), [['a.gb', ['plasmid', 'chr']], ['b.gb', ['phage']]]);
  owner.setRecordListQuery('pha');
  assert.deepEqual(view(), [['b.gb', ['phage']]]);
  assert.deepEqual(owner.recordListView.value.rows.map((row) => row.key), ['u2']);
});

test('a list of one file has one group without a name, as before', async () => {
  const source = { key: 'f', name: 'f.gb', records: [{ key: 'f-1', recordId: 'x', length: 1 }] };
  const owner = createRecordSelection({
    reactive, computed: Vue.computed,
    recordsOff: () => [], source: (_mode, key) => (key === 'f' ? source : null),
    runUndoable: async (_label, change) => change(), withDialogChoice: (_label, handler) => handler,
    closeAfterDialogChoice: (close) => close(), afterRecordSetChange: () => {}, removeSource: () => true,
    askRemoveSource: () => null, deleteRecordSettings: () => true, prepareRecords: async () => {},
    autoOpenRequests: () => []
  });
  await owner.openRecordList('linear', 'f');
  assert.deepEqual(owner.recordListView.value.groups.map((group) => [group.name, group.rows.length]), [['', 1]]);
});
