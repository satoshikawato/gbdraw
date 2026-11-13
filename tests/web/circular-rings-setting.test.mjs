// OV-2818 (Owner 2026-11-09 「Q3=A」): a Circular drawing's Show Comparison
// Rings follows its ring files as Show Depth follows the Depth files: the
// first ring file turns it on, adding a file turns it on, removing the last
// turns it off, the source radio follows the files it switches to, and Reset
// Settings turns it off while the ring rows and their files stay. Session Load
// reads a Session 51's stored setting and turns an older Session's rings on by
// its ring files.
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { test } from 'node:test';
import { modeDrawing } from './helpers/drawing-state.mjs';

globalThis.window = {
  Vue: {
    ref: (value) => ({ value }), reactive: (value) => value,
    computed: (getter) => ({ get value() { return getter(); } }),
    nextTick: async () => {}
  },
  DOMPurify: { sanitize: (value) => value }
};

const { state } = await import('../../gbdraw/web/js/state.js');
await (await import('./helpers/drawing-state.mjs')).addStartLinearDrawing();
const {
  createCircularConservationSources, createCircularConservationSeriesPanel, createCircularConservationComparisonFiles
} = await import('../../gbdraw/web/js/app/circular-conservation.js');
const { createDrawingReset } = await import('../../gbdraw/web/js/app/reset-settings.js');
const { bindDrawingInputFileList, bindDrawingInputFileTree, setConservationBlastsSource } = await import('../../gbdraw/web/js/services/drawing/inputs.js');
const { sessionSourceFactsOf } = await import('../../gbdraw/web/js/services/session-drawings-upgrade.js');
const { circularRingSources, orderedConservationSources } = await import('../../gbdraw/web/js/services/conservation-series.js');

const circular = () => modeDrawing(state, 'circular');
// `arrayBuffer` lets Generate's serializer read it (`serializeActiveRenderFiles`).
const textFile = (name, text = '>ring\nACGT\n') => ({
  name, type: 'text/plain', size: text.length, lastModified: 0, data: btoa(text),
  arrayBuffer: async () => new TextEncoder().encode(text).buffer
});
const shownAction = (mode, action) => action(modeDrawing(state, mode));
const askingStackPorts = () => ({
  askChoice: async () => null, runUndoable: (_label, fn) => fn(), busy: () => false, rows: () => [], mode: 'circular'
});

const sources = createCircularConservationSources({
  state,
  editLosatRows: () => {},
  syncCircularConservationSlots: () => {},
  rekeyCircularConservationRows: () => {},
  circularConservationRowsLeaving: () => new Set(),
  askingStackPorts
});
const panel = createCircularConservationSeriesPanel({
  state, shownAction, shownDrawingOfMode: (mode) => modeDrawing(state, mode), editLosatRows: () => {},
  circularConservationFastaInput: { value: null },
  reconcileCircularConservation: sources.reconcileCircularConservation,
  refreshCircularConservationSlots: sources.refreshCircularConservationSlots,
  askRingSourceRemoval: sources.askRingSourceRemoval
});
const comparisonFiles = createCircularConservationComparisonFiles({
  state, shownAction, editLosatRows: () => {},
  // The step's commit waits for the record label read, which no Worker answers here.
  runUndoable: async () => undefined,
  setRingGenomeLists: sources.setRingGenomeLists,
  setCircularConservationSourceFiles: sources.setCircularConservationSourceFiles,
  reconcileCircularConservation: sources.reconcileCircularConservation,
  refreshCircularConservationSlots: sources.refreshCircularConservationSlots,
  askRingSourceRemoval: sources.askRingSourceRemoval
});
const { resetDrawing } = createDrawingReset({
  state,
  circularTrackNewRenderer: { value: '' },
  linearTrackNewRenderer: { value: '' },
  rederiveCircularConservation: sources.rederiveCircularConservation,
  applyLinearSeqMutation: () => {},
  invalidateLinearComparisonArtifacts: () => {},
  writeArtifactValues: () => {}
});

const addSeq = (name) => comparisonFiles.addCircularConservationComparisonFile(
  /** @type {any} */ ({ target: { files: [textFile(name)], value: name } })
);
const ringNames = () => circular().circularConservation.series.map((entry) => entry.fileName);

test('the first comparison file turns the rings on and removing the last turns them off', async () => {
  const drawing = circular();
  assert.equal(drawing.circularConservation.source, 'losat');
  assert.equal(drawing.circularConservation.enabled, false);
  await addSeq('zebrafish.fasta');
  assert.equal(drawing.circularConservation.enabled, true);
  await addSeq('fly.fasta');
  assert.deepEqual(ringNames(), ['zebrafish.fasta', 'fly.fasta']);

  panel.setCircularConservationShown(false);
  assert.equal(drawing.circularConservation.enabled, false, 'the checkbox turns the rings off');
  comparisonFiles.removeCircularConservationSource(1);
  assert.equal(drawing.circularConservation.enabled, false, 'removing a file while others stay keeps the setting');
  await addSeq('worm.fasta');
  assert.equal(drawing.circularConservation.enabled, true, 'adding a file turns the rings on');

  comparisonFiles.removeCircularConservationSource(1);
  assert.equal(drawing.circularConservation.enabled, true);
  comparisonFiles.removeCircularConservationSource(0);
  assert.deepEqual(ringNames(), []);
  assert.equal(drawing.circularConservation.enabled, false, 'removing the last file turns the rings off');
});

test('the source radio follows the ring files of the source it switches to', async () => {
  const drawing = circular();
  await addSeq('zebrafish.fasta');
  assert.equal(drawing.circularConservation.enabled, true);
  panel.changeCircularConservationSource('upload');
  assert.equal(drawing.circularConservation.enabled, false, 'Upload has no BLAST table');
  panel.setCircularConservationUploadFiles([textFile('zebrafish.blast.tsv', 'q\ts\t100\t8\t0\t0\t1\t8\t1\t8\t1e-9\t20\n')]);
  assert.equal(drawing.circularConservation.enabled, true, 'the first BLAST table turns the rings on');
  panel.setCircularConservationShown(false);
  panel.changeCircularConservationSource('losat');
  assert.equal(drawing.circularConservation.enabled, false, 'a switch between two sourced sources keeps the setting');
  panel.changeCircularConservationSource('upload');
  panel.setCircularConservationUploadFiles([]);
  panel.changeCircularConservationSource('losat');
  bindDrawingInputFileTree(state.inputFiles, drawing.slices.inputs, 'c_conservation_fastas', []);
  sources.reconcileCircularConservation(drawing);
});

test('Reset Settings turns the rings off and keeps the ring rows and their files', async () => {
  const drawing = circular();
  await addSeq('zebrafish.fasta');
  await addSeq('fly.fasta');
  drawing.circularConservation.series[0].label = 'Danio';
  assert.equal(drawing.circularConservation.enabled, true);

  resetDrawing(drawing);

  assert.equal(drawing.circularConservation.enabled, false);
  assert.deepEqual(drawing.files.c_conservation_fastas.map((files) => files.map((file) => file.name)), [['zebrafish.fasta'], ['fly.fasta']]);
  assert.deepEqual(ringNames(), ['zebrafish.fasta', 'fly.fasta']);
  assert.notEqual(drawing.circularConservation.series[0].label, 'Danio', 'the ring rows return to their defaults');
  // Turning the setting on again draws the kept rows.
  panel.setCircularConservationShown(true);
  assert.equal(drawing.circularConservation.enabled, true);
});

test('Load keeps a Session 51 rings setting and reads an older Session by its ring files', () => {
  const drawing = circular();
  assert.ok(ringNames().length > 0);
  assert.equal(sessionSourceFactsOf(51, false).ringSettingStored, true);
  assert.equal(sessionSourceFactsOf(50, true).ringSettingStored, false);
  assert.equal(sessionSourceFactsOf(46, true).ringSettingStored, false);

  drawing.circularConservation.enabled = false;
  sources.reconcileLoadedCircularConservation(drawing, true);
  assert.equal(drawing.circularConservation.enabled, false, 'a Session 51 keeps the rings off');
  sources.reconcileLoadedCircularConservation(drawing, false);
  assert.equal(drawing.circularConservation.enabled, true, 'an older Session with ring files shows its rings');
  bindDrawingInputFileTree(state.inputFiles, drawing.slices.inputs, 'c_conservation_fastas', []);
  sources.reconcileLoadedCircularConservation(drawing, false);
  assert.equal(drawing.circularConservation.enabled, false, 'an older Session without ring files has none');
});

// PR-B (D-02, Owner 2026-11-08 「すべて推奨通り」): a ring compares one genome,
// which may be several files. Its row adds files to that genome and removes
// them, and keeps its name and place; a LOSAT ring keeps its last file (Remove
// series removes the ring). The files keep the order they were added in: their
// order changes no hit drawn, and a new order would only search again (D-106).
test('a ring row adds and removes the files of its genome, in the order added, and keeps its name', async () => {
  const drawing = circular();
  await addSeq('chr1.fasta');
  await addSeq('fly.fasta');
  drawing.circularConservation.series[0].label = 'Two chromosomes';
  const rows = () => panel.circularConservationSeriesRows.value;
  const boundNames = () => drawing.files.c_conservation_fastas.map((files) => files.map((file) => file.name));
  comparisonFiles.addRingRowFiles(0, /** @type {any} */ ({ target: { files: [textFile('chr2.fasta'), textFile('chr3.fasta')], value: 'x' } }));
  assert.deepEqual(rows().map((row) => row.files), [['chr1.fasta', 'chr2.fasta', 'chr3.fasta'], ['fly.fasta']]);
  assert.deepEqual(boundNames(), [['chr1.fasta', 'chr2.fasta', 'chr3.fasta'], ['fly.fasta']]);
  assert.equal(drawing.circularConservation.series.length, 2, 'added files join the ring, not new rings');
  assert.equal(drawing.circularConservation.series[0].label, 'Two chromosomes');
  assert.match(String(drawing.circularConservation.series[0].sourceKey), /chr1\.fasta.* \+ .*chr2\.fasta.* \+ .*chr3\.fasta/);

  assert.equal('moveRingRowFile' in comparisonFiles, false, 'a ring row has no file order edit');
  const ringRowFileControls = readFileSync(new URL('../../gbdraw/web/index.html', import.meta.url), 'utf8')
    .split('data-ring-row-files')[1].split('data-ring-row-file-add')[0];
  assert.doesNotMatch(ringRowFileControls, /ph-arrow-(up|down)|Move File/, 'a ring row file has no move buttons');
  comparisonFiles.removeRingRowFile(0, 0);
  assert.deepEqual(rows()[0].files, ['chr2.fasta', 'chr3.fasta']);
  assert.equal(drawing.circularConservation.series[0].label, 'Two chromosomes');
  assert.deepEqual(rows().map((row) => row.canRemoveFile), [true, false]);
  comparisonFiles.removeRingRowFile(1, 0);
  assert.deepEqual(boundNames(), [['chr2.fasta', 'chr3.fasta'], ['fly.fasta']], 'a LOSAT ring keeps its last file');

  comparisonFiles.removeCircularConservationSource(1);
  comparisonFiles.removeCircularConservationSource(0);
  assert.deepEqual(ringNames(), []);
});

// Review 1 of #1777, item 1: a ring row reads and writes its files by the ring
// identity Generate pairs a name and color with (the series entry's source
// key), so two first files of one name and a moved series keep each genome
// with its own name.
test('a ring row reads and writes its files by the identity Generate uses, also for first files of one name', async () => {
  const drawing = circular();
  const add = (file) => comparisonFiles.addCircularConservationComparisonFile(/** @type {any} */ ({ target: { files: [file], value: 'x' } }));
  await add(textFile('sequence.fasta', '>A\nACGT\n'));
  await add(textFile('sequence.fasta', '>B\nACGTACGTACGT\n'));
  const [first, second] = drawing.circularConservation.series;
  first.label = 'Genome A';
  second.label = 'Genome B';
  panel.moveCircularConservationSeries(1, -1);
  const generated = () => orderedConservationSources(circularRingSources(drawing), drawing.circularConservation)
    .map((source) => [source.label, source.files.map((file) => file.size)]);
  const shown = () => panel.circularConservationSeriesRows.value
    .map((row) => [drawing.circularConservation.series[row.index].label, row.files]);
  assert.deepEqual(generated(), [['Genome B', [16]], ['Genome A', [8]]]);
  assert.deepEqual(shown(), [['Genome B', ['sequence.fasta']], ['Genome A', ['sequence.fasta']]]);

  comparisonFiles.addRingRowFiles(0, /** @type {any} */ ({ target: { files: [textFile('plasmid.fasta', '>P\nAC\n')], value: 'x' } }));
  assert.deepEqual(generated(), [['Genome B', [16, 6]], ['Genome A', [8]]], 'the plasmid joins the genome of its row');
  assert.deepEqual(shown(), [['Genome B', ['sequence.fasta', 'plasmid.fasta']], ['Genome A', ['sequence.fasta']]]);
  comparisonFiles.removeRingRowFile(0, 0);
  assert.deepEqual(generated(), [['Genome B', [6]], ['Genome A', [8]]], 'the file removed is the one of its row');

  comparisonFiles.removeCircularConservationSource(1);
  comparisonFiles.removeCircularConservationSource(0);
  assert.deepEqual(ringNames(), []);
});

// PR-B review 2, HIGH 1 (a): Remove series on the first of two rings whose first files share a name leaves the
// second ring with its own name, colour and genome.
test('Remove series leaves the other ring of the same file name with its own name, colour and genome', async () => {
  const drawing = circular();
  const add = (file) => comparisonFiles.addCircularConservationComparisonFile(/** @type {any} */ ({ target: { files: [file], value: 'x' } }));
  await add(textFile('sequence.fasta', '>A\nACGT\n'));
  await add(textFile('sequence.fasta', '>B\nACGTACGTACGT\n'));
  const [first, second] = drawing.circularConservation.series;
  first.label = 'Genome A';
  first.color = '#111111';
  second.label = 'Genome B';
  const colorB = second.color;
  comparisonFiles.removeCircularConservationSource(0);
  assert.deepEqual(drawing.circularConservation.series.map((entry) => [entry.label, entry.color]), [['Genome B', colorB]]);
  assert.deepEqual(orderedConservationSources(circularRingSources(drawing), drawing.circularConservation)
    .map((source) => [source.label, source.files.map((file) => file.size)]), [['Genome B', [16]]]);
  comparisonFiles.removeCircularConservationSource(0);
  assert.deepEqual(ringNames(), []);
});

// #1777 review 3, finding 1: two rings of one table name (a LOSAT-cache replay names each table after its first FASTA
// file, and NCBI names every download `sequence.fasta`) keep their own tables and genomes through a move and Remove
// series; their stored places tell them apart.
test('a moved replay ring of a shared table name keeps its table and genome, also after Remove series', () => {
  const drawing = circular();
  const shared = 'sequence.circular_conservation.losatn.tsv';
  const tables = [textFile(shared, 'b'.repeat(60)), textFile(shared, 'c'.repeat(63)), textFile('d.circular_conservation.losatn.tsv', 'd'.repeat(64))];
  const genomes = ['b', 'c', 'd'].map((name) => [textFile(`${name}.fasta`)]);
  bindDrawingInputFileTree(state.inputFiles, drawing.slices.inputs, 'c_conservation_fastas', genomes);
  bindDrawingInputFileList(state.inputFiles, drawing.slices.inputs, 'c_conservation_blasts', tables);
  setConservationBlastsSource(drawing.slices.inputs, 'losat-cache');
  drawing.circularConservation.series.splice(0, drawing.circularConservation.series.length, ...['b', 'c', 'd'].map((name, index) => ({
    sourceKey: `${name}.fasta|10|0|0`, fileName: tables[index].name, sourceIndex: index, label: `Genome ${name}`, color: `#00000${index}`, losat_gencode: 1
  })));
  const drawn = () => orderedConservationSources(drawing.files.c_conservation_blasts, drawing.circularConservation)
    .map((source) => [source.label, source.file.size]);
  const shown = () => panel.circularConservationSeriesRows.value.map((row) => [drawing.circularConservation.series[row.index].label, row.files]);
  panel.moveCircularConservationSeries(0, 1);
  assert.deepEqual(drawn(), [['Genome c', 63], ['Genome b', 60], ['Genome d', 64]]);
  assert.deepEqual(shown(), [['Genome c', ['c.fasta']], ['Genome b', ['b.fasta']], ['Genome d', ['d.fasta']]]);
  comparisonFiles.removeCircularConservationSource(1);
  assert.deepEqual(drawn(), [['Genome c', 63], ['Genome d', 64]]);
  assert.deepEqual(shown(), [['Genome c', ['c.fasta']], ['Genome d', ['d.fasta']]]);
  comparisonFiles.removeCircularConservationSource(1);
  comparisonFiles.removeCircularConservationSource(0);
  assert.deepEqual(ringNames(), []);
});

// The same for LOSAT rings whose first files share a name and whose source keys no longer match (keys of other dates):
// a move and Remove series leave each kept ring its own genome, name and colour.
test('moved LOSAT rings of one first file name keep their genomes through Remove series when their keys are stale', async () => {
  const drawing = circular();
  const add = (file) => comparisonFiles.addCircularConservationComparisonFile(/** @type {any} */ ({ target: { files: [file], value: 'x' } }));
  for (const bases of ['ACGT', 'ACGTACGT', 'ACGTACGTACGT']) await add(textFile('sequence.fasta', `>r\n${bases}\n`));
  drawing.circularConservation.series.forEach((entry, index) => {
    Object.assign(entry, { sourceKey: `stale-${index}`, label: `Genome ${'ABC'[index]}`, color: `#11111${index}` });
  });
  const generated = () => orderedConservationSources(circularRingSources(drawing), drawing.circularConservation)
    .map((source) => [source.label, source.color, source.files.map((file) => file.size)]);
  panel.moveCircularConservationSeries(0, 1);
  assert.deepEqual(generated(), [['Genome B', '#111111', [12]], ['Genome A', '#111110', [8]], ['Genome C', '#111112', [16]]]);
  comparisonFiles.removeCircularConservationSource(1);
  assert.deepEqual(generated(), [['Genome B', '#111111', [12]], ['Genome C', '#111112', [16]]]);
  comparisonFiles.removeCircularConservationSource(1);
  comparisonFiles.removeCircularConservationSource(0);
  assert.deepEqual(ringNames(), []);
});

// An uploaded BLAST table's ring takes its comparison genome's files the same way; they are optional.
test('an uploaded table ring row takes the files of its comparison genome, and may have none', () => {
  const drawing = circular();
  panel.changeCircularConservationSource('upload');
  panel.setCircularConservationUploadFiles([textFile('two.blast.tsv', 'q\ts\t100\t8\t0\t0\t1\t8\t1\t8\t1e-9\t20\n')]);
  const rows = () => panel.circularConservationSeriesRows.value;
  assert.deepEqual(rows().map((row) => [row.files, row.canRemoveFile]), [[[], false]]);
  comparisonFiles.addRingRowFiles(0, /** @type {any} */ ({ target: { files: [textFile('chr1.fasta'), textFile('chr2.fasta')], value: 'x' } }));
  assert.deepEqual(rows()[0].files, ['chr1.fasta', 'chr2.fasta']);
  assert.deepEqual(drawing.files.c_conservation_sequence_sources.map((files) => files.map((file) => file.name)), [['chr1.fasta', 'chr2.fasta']]);
  comparisonFiles.removeRingRowFile(0, 1);
  comparisonFiles.removeRingRowFile(0, 0);
  assert.deepEqual(rows()[0].files, []);
  assert.equal(drawing.circularConservation.series.length, 1, 'the table keeps its ring');
  panel.setCircularConservationUploadFiles([]);
  panel.changeCircularConservationSource('losat');
});

// #1777 review 4: three LOSAT rings whose files share name, size and date (NCBI's `sequence.fasta`, 8 B, no date) keep
// their own genomes through Remove series of the first or the middle one, and through a Session 51 Load that lost one
// of their files (its place stays empty, R1-B1): the ring that lost its file is the one that goes, every other ring
// keeps its genome, and a later Remove series or Add Seq keeps them so.
test('LOSAT rings of files of one name, size and date keep their genomes through Remove series and a Load that lost one', async () => {
  const drawing = circular();
  const bytesOf = (ref) => atob(state.inputFiles.fileOf(ref).data).split('\n')[1];
  const generated = () => orderedConservationSources(circularRingSources(drawing), drawing.circularConservation)
    .map((source) => [source.label, bytesOf(source.file)]);
  const shown = () => panel.circularConservationSeriesRows.value.map((row) => [drawing.circularConservation.series[row.index].label, row.files]);
  const sequenceFasta = (bases) => textFile('sequence.fasta', `>r\n${bases}\n`);
  const add = (bases) => comparisonFiles.addCircularConservationComparisonFile(/** @type {any} */ ({ target: { files: [sequenceFasta(bases)], value: 'x' } }));
  const clear = () => {
    while (drawing.circularConservation.series.length > 0) comparisonFiles.removeCircularConservationSource(0);
    bindDrawingInputFileTree(state.inputFiles, drawing.slices.inputs, 'c_conservation_fastas', []);
    assert.deepEqual(ringNames(), []);
  };
  for (const removed of [0, 1]) {
    for (const bases of ['ACGT', 'TTTT', 'GGGG']) await add(bases);
    drawing.circularConservation.series.forEach((entry, index) => Object.assign(entry, { label: `Genome ${'ABC'[index]}` }));
    comparisonFiles.removeCircularConservationSource(removed);
    assert.deepEqual(generated(), [['Genome A', 'ACGT'], ['Genome B', 'TTTT'], ['Genome C', 'GGGG']].filter((_, index) => index !== removed),
      `Remove series of ring ${removed}`);
    clear();
  }
  for (const lost of [0, 1]) {
    // What Load leaves of a Session 51 whose genome binding `lost` it could not read, then its reconcile (`afterImport`).
    const files = ['AAAA', 'CCCC', 'GGGG'].map((bases) => sequenceFasta(bases));
    bindDrawingInputFileTree(state.inputFiles, drawing.slices.inputs, 'c_conservation_fastas', files.map((file, index) => (index === lost ? [] : [file])));
    drawing.circularConservation.series.splice(0, drawing.circularConservation.series.length, ...['A', 'C', 'G'].map((name, index) => ({
      sourceKey: `sequence.fasta|8|0|${index}`, fileName: 'sequence.fasta', sourceIndex: index, label: `ring ${name}`, color: `#22222${index}`, losat_gencode: 1
    })));
    sources.reconcileLoadedCircularConservation(drawing, true);
    const kept = [['ring A', 'AAAA'], ['ring C', 'CCCC'], ['ring G', 'GGGG']].filter((_, index) => index !== lost);
    assert.deepEqual(generated(), kept, `Load without genome ${lost}`);
    assert.deepEqual(shown(), kept.map(([label]) => [label, ['sequence.fasta']]), `rows after Load without genome ${lost}`);
    sources.reconcileCircularConservation(drawing);
    assert.deepEqual(generated(), kept, `a later reconcile after Load without genome ${lost}`);
    await add('TTTT');
    assert.deepEqual(generated(), [...kept, ['sequence', 'TTTT']], `Add Seq after Load without genome ${lost}`);
    comparisonFiles.removeCircularConservationSource(0);
    assert.deepEqual(generated(), [kept[1], ['sequence', 'TTTT']], `Remove series after Load without genome ${lost}`);
    clear();
  }
});

// #1777 review 5: with the track stack on, each managed ring row stays with its ring when the source keys shift (three
// rings of NCBI `sequence.fasta` files of one size and date): through Remove series, whose reach dialog takes out the
// removed ring's row, through a Session 51 Load that lost one ring's file, whose row is the one that goes, and through
// a ring row's file edit. The stack's ring rows read the stack editor itself here.
const { createCircularTrackSlotEditor } = await import('../../gbdraw/web/js/app/circular-track-slots.js');
const stackEditor = createCircularTrackSlotEditor({ state });
/** @type {string[][]} */
const leavingAsked = [];
const stackSources = createCircularConservationSources({
  state,
  editLosatRows: () => {},
  syncCircularConservationSlots: (drawing) => stackEditor.syncCircularConservationSlots(drawing),
  rekeyCircularConservationRows: (drawing, keyChanges) => stackEditor.rekeyCircularConservationRows(drawing, keyChanges),
  circularConservationRowsLeaving: (drawing, entries, keyChanges) => {
    const leaving = stackEditor.circularConservationRowsLeaving(drawing, entries, keyChanges);
    leavingAsked.push([...leaving]);
    return leaving;
  },
  askingStackPorts
});
const stackFiles = createCircularConservationComparisonFiles({
  state, shownAction, editLosatRows: () => {},
  runUndoable: async () => undefined,
  setRingGenomeLists: stackSources.setRingGenomeLists,
  setCircularConservationSourceFiles: stackSources.setCircularConservationSourceFiles,
  reconcileCircularConservation: stackSources.reconcileCircularConservation,
  refreshCircularConservationSlots: stackSources.refreshCircularConservationSlots,
  askRingSourceRemoval: stackSources.askRingSourceRemoval
});
const stackPanel = createCircularConservationSeriesPanel({
  state, shownAction, shownDrawingOfMode: (mode) => modeDrawing(state, mode), editLosatRows: () => {},
  circularConservationFastaInput: { value: null },
  reconcileCircularConservation: stackSources.reconcileCircularConservation,
  refreshCircularConservationSlots: stackSources.refreshCircularConservationSlots,
  askRingSourceRemoval: stackSources.askRingSourceRemoval
});
const stackRingRows = () => circular().adv.circular_track_slots.filter((slot) => slot.renderer === 'sequence_conservation')
  .map((slot) => [slot.id, slot.width, slot.params.label]);
const stackGenerated = () => {
  const drawing = circular();
  return orderedConservationSources(circularRingSources(drawing), drawing.circularConservation)
    .map((source) => [source.label, atob(state.inputFiles.fileOf(source.file).data).split('\n')[1]]);
};
const withStack = async (run) => {
  const drawing = circular();
  const before = { enabled: drawing.adv.circular_track_slots_enabled, slots: [...drawing.adv.circular_track_slots], shown: drawing.circularConservation.enabled };
  drawing.adv.circular_track_slots_enabled = true;
  drawing.circularConservation.enabled = true;
  try {
    await run(drawing);
  } finally {
    while (drawing.circularConservation.series.length > 0) stackFiles.removeCircularConservationSource(0);
    bindDrawingInputFileTree(state.inputFiles, drawing.slices.inputs, 'c_conservation_fastas', []);
    drawing.adv.circular_track_slots.splice(0, drawing.adv.circular_track_slots.length, ...before.slots);
    drawing.adv.circular_track_slots_enabled = before.enabled;
    drawing.circularConservation.enabled = before.shown;
    leavingAsked.splice(0);
  }
};
const sequenceFastaOf = (bases) => textFile('sequence.fasta', `>r\n${bases}\n`);
/** Each ring row's own edit: its width. */
const widenRingRows = () => circular().adv.circular_track_slots.filter((slot) => slot.renderer === 'sequence_conservation')
  .forEach((slot, index) => { slot.width = `${10 + index}px`; });

test('the stack keeps each ring row with its ring through Remove series of rings of one file name, size and date', async () => {
  for (const removed of [0, 1, 2]) {
    await withStack(async (drawing) => {
      for (const bases of ['ACGT', 'TTTT', 'GGGG']) {
        await stackFiles.addCircularConservationComparisonFile(/** @type {any} */ ({ target: { files: [sequenceFastaOf(bases)], value: 'x' } }));
      }
      drawing.circularConservation.series.forEach((entry, index) => Object.assign(entry, { label: `G${'ABC'[index]}` }));
      stackSources.refreshCircularConservationSlots(drawing);
      widenRingRows();
      const rows = stackRingRows();
      assert.deepEqual(rows.map(([, width, label]) => [width, label]), [['10px', 'GA'], ['11px', 'GB'], ['12px', 'GC']]);
      leavingAsked.splice(0);
      stackFiles.removeCircularConservationSource(removed);
      assert.deepEqual(leavingAsked, [[rows[removed][0]]], `the reach dialog of Remove series ${removed} takes out the removed ring's row`);
      assert.deepEqual(stackRingRows(), rows.filter((_, index) => index !== removed), `rows after Remove series ${removed}`);
      assert.deepEqual(stackGenerated(), [['GA', 'ACGT'], ['GB', 'TTTT'], ['GC', 'GGGG']].filter((_, index) => index !== removed));
    });
  }
});

test('the stack keeps each ring row with its ring through a Session 51 Load that lost one ring file of one name, size and date', async () => {
  for (const lost of [0, 1, 2]) {
    await withStack(async (drawing) => {
      const files = ['AAAA', 'CCCC', 'GGGG'].map((bases) => sequenceFastaOf(bases));
      bindDrawingInputFileTree(state.inputFiles, drawing.slices.inputs, 'c_conservation_fastas', files.map((file) => [file]));
      drawing.circularConservation.series.splice(0, drawing.circularConservation.series.length, ...['A', 'C', 'G'].map((name, index) => ({
        sourceKey: `sequence.fasta|8|0|${index}`, fileName: 'sequence.fasta', sourceIndex: index, label: `ring ${name}`, color: `#22222${index}`, losat_gencode: 1
      })));
      stackSources.refreshCircularConservationSlots(drawing);
      widenRingRows();
      const rows = stackRingRows();
      // What Load leaves of that Session without genome `lost`: its place empty, the series and rows as saved; then `afterImport`.
      bindDrawingInputFileTree(state.inputFiles, drawing.slices.inputs, 'c_conservation_fastas', files.map((file, index) => (index === lost ? [] : [file])));
      stackSources.reconcileLoadedCircularConservation(drawing, true);
      assert.deepEqual(stackRingRows(), rows.filter((_, index) => index !== lost), `rows after Load without genome ${lost}`);
      assert.deepEqual(stackGenerated(), [['ring A', 'AAAA'], ['ring C', 'CCCC'], ['ring G', 'GGGG']].filter((_, index) => index !== lost));
      await stackFiles.addCircularConservationComparisonFile(/** @type {any} */ ({ target: { files: [sequenceFastaOf('TTTT')], value: 'x' } }));
      assert.deepEqual(stackRingRows().slice(0, 2), rows.filter((_, index) => index !== lost), `rows after Add Seq after Load without genome ${lost}`);
    });
  }
});

test('the stack keeps each ring row with its ring when a ring row of one file name, size and date gets another file', async () => {
  await withStack(async (drawing) => {
    for (const bases of ['ACGT', 'TTTT', 'GGGG']) {
      await stackFiles.addCircularConservationComparisonFile(/** @type {any} */ ({ target: { files: [sequenceFastaOf(bases)], value: 'x' } }));
    }
    drawing.circularConservation.series.forEach((entry, index) => Object.assign(entry, { label: `G${'ABC'[index]}` }));
    stackSources.refreshCircularConservationSlots(drawing);
    widenRingRows();
    const rows = stackRingRows();
    stackFiles.addRingRowFiles(0, /** @type {any} */ ({ target: { files: [textFile('plasmid.fasta')], value: 'x' } }));
    assert.deepEqual(stackRingRows(), rows, 'Add file on the first ring row');
    stackFiles.removeRingRowFile(0, 1);
    assert.deepEqual(stackRingRows(), rows, 'its remove button');
  });
});

test('the stack keeps each ring row with its ring through Remove series of replayed rings whose tables share name, size and date', async () => {
  for (const removed of [0, 1]) {
    await withStack(async (drawing) => {
      const shared = 'sequence.circular_conservation.losatn.tsv';
      const tables = ['b', 'c', 'd'].map((letter) => textFile(shared, letter.repeat(60)));
      bindDrawingInputFileTree(state.inputFiles, drawing.slices.inputs, 'c_conservation_fastas', ['b', 'c', 'd'].map((name) => [textFile(`${name}.fasta`)]));
      bindDrawingInputFileList(state.inputFiles, drawing.slices.inputs, 'c_conservation_blasts', tables);
      setConservationBlastsSource(drawing.slices.inputs, 'losat-cache');
      drawing.circularConservation.series.splice(0, drawing.circularConservation.series.length, ...['b', 'c', 'd'].map((name, index) => ({
        sourceKey: `${name}.fasta|10|0|0`, fileName: shared, sourceIndex: index, label: `Genome ${name}`, color: `#00000${index}`, losat_gencode: 1
      })));
      stackSources.refreshCircularConservationSlots(drawing);
      widenRingRows();
      const rows = stackRingRows();
      assert.deepEqual(rows.map(([, width, label]) => [width, label]), [['10px', 'Genome b'], ['11px', 'Genome c'], ['12px', 'Genome d']]);
      leavingAsked.splice(0);
      stackFiles.removeCircularConservationSource(removed);
      assert.deepEqual(leavingAsked, [[rows[removed][0]]], `the reach dialog of Remove series ${removed} takes out the removed ring's row`);
      assert.deepEqual(stackRingRows(), rows.filter((_, index) => index !== removed), `rows after Remove series ${removed}`);
      const drawn = orderedConservationSources(drawing.files.c_conservation_blasts, drawing.circularConservation)
        .map((source) => [source.label, atob(state.inputFiles.fileOf(source.file).data)[0]]);
      assert.deepEqual(drawn, [['Genome b', 'b'], ['Genome c', 'c'], ['Genome d', 'd']].filter((_, index) => index !== removed));
    });
  }
});

// #1777 review 6: after a LOSAT Generate's adoption (its tables at the rings' places, the rows moved to their keys),
// every later ring transition still keeps each ring's label, colour and stack row (id, width), and each ring keeps one
// identity (`ringIdentityDefects`): Add Seq, a ring row's file edit, Remove series, and the source radio. The
// adoption runs the real path: Generate's canonical request (`buildCanonicalRenderRequest`) adopted by
// `adoptCanonicalRenderArtifacts`, through the rings' adoption transition of the composition root.
const sessionConfig = await import('../../gbdraw/web/js/services/config.js');
const { buildCanonicalRenderRequest } = await import('../../gbdraw/web/js/services/session-request.js');
const { ringGenomeFileLists } = await import('../../gbdraw/web/js/services/drawing/inputs.js');
const { ringGenomeFilesInRingOrder } = await import('../../gbdraw/web/js/services/conservation-series.js');
const { NO_CONSERVATION_PAIRS, losatRequestRings, searchedRingSources } = await import('../../gbdraw/web/js/services/circular-conservation-keys.js');
const { readFileText } = await import('../../gbdraw/web/js/services/file-content-cache.js');
const { ringIdentityDefects } = await import('./helpers/ring-identity.mjs');
const REFERENCE_GENBANK = 'LOCUS       REFERENCE                   8 bp    DNA     linear   UNK 01-JAN-1980\nFEATURES             Location/Qualifiers\nORIGIN\n        1 aaccggtt\n//\n';
/** The reference record a Generate's request draws (as `session-losat-cache-validation.test.mjs` sets it). */
const placeReference = (drawing) => {
  drawing.circularSources.splice(0, drawing.circularSources.length,
    { uid: 'record-1', gb: state.inputFiles.registerFile(textFile('reference.gb', REFERENCE_GENBANK)), gff: null, fasta: null });
  drawing.circularRecordList.value = [{ uid: 'record-1', k: 1, identity: 'record-1', selector: '#1', record_id: 'REFERENCE',
    record_length: 8, detectedTopology: 'linear' }];
};
/**
 * A LOSAT Generate of the draft's rings as production builds it (the rings its ring keys read, `losatRequestRings`,
 * the request builder), adopted with the rings it searched: each ring's table is named after its first genome file,
 * as Generate names it, so three rings of `sequence.fasta` files give three tables of one name.
 * @param {import('../../gbdraw/web/js/state.js').DrawingState} drawing
 */
const generateRings = async (drawing) => {
  const searched = searchedRingSources(drawing);
  const filesData = await sessionConfig.serializeActiveRenderFiles(state, drawing);
  const tables = searched.map((ring, index) => ({
    name: `${ring.fileName.replace(/\.[^.]+$/, '')}.circular_conservation.losatn.tsv`,
    text: `q\ts\t100\t8\t0\t0\t1\t8\t1\t8\t1e-9\t2${index}\n`
  }));
  const canonical = buildCanonicalRenderRequest({
    state, drawing, filesData,
    resolvedCircularConservation: losatRequestRings(searched, tables, ringGenomeFileLists(filesData.c_conservation_fastas))
  });
  stackSources.adoptCircularRingResult(drawing, (installRings) => sessionConfig.adoptCanonicalRenderArtifacts(canonical, {
    adoptOwnedRequest: true, drawing, conservationPairs: NO_CONSERVATION_PAIRS, installRings
  }), { searched, builtFromDraft: true });
  assert.equal(drawing.files.c_conservation_blasts_source, 'losat-cache', 'the Generate adopted a LOSAT-cache result');
  return canonical;
};
/**
 * A rerender of the committed request (a Rotate): it replaces only the committed Session (R1 (c)).
 * @param {import('../../gbdraw/web/js/state.js').DrawingState} drawing
 * @param {Parameters<typeof sessionConfig.adoptCanonicalRenderArtifacts>[0]} canonical
 */
const rerenderCommitted = (drawing, canonical) => {
  stackSources.adoptCircularRingResult(drawing, (installRings) => sessionConfig.adoptCanonicalRenderArtifacts(canonical, {
    adoptOwnedRequest: true, drawing, conservationPairs: 'committed', installRings
  }), null);
};
const ringLooks = () => circular().circularConservation.series.map((entry) => [entry.label, entry.color]);
/** Each drawn ring's label and bases; a Generate binds its own views (Files) of the ring files. */
const drawnRings = async () => {
  const drawing = circular();
  const bases = async (file) => (await readFileText(state.inputFiles.fileOf(file))).split('\n')[1];
  const genomes = ringGenomeFilesInRingOrder(drawing);
  return Promise.all(orderedConservationSources(circularRingSources(drawing), drawing.circularConservation)
    .map(async (source, index) => [source.label, genomes[index].length > 0 ? await bases(genomes[index][0]) : null]));
};
/**
 * Three rings of NCBI `sequence.fasta` files of one size and date (or, with `distinct`, of files of their own names),
 * named GA, GB, GC, their rows 10, 11 and 12 px wide.
 */
const threeSequenceRings = async (drawing, { distinct = false } = {}) => {
  placeReference(drawing);
  for (const bases of ['ACGT', 'TTTT', 'GGGG']) {
    const file = distinct ? textFile(`genome_${bases}.fasta`, `>r\n${bases}\n`) : sequenceFastaOf(bases);
    await stackFiles.addCircularConservationComparisonFile(/** @type {any} */ ({ target: { files: [file], value: 'x' } }));
  }
  // Colours of the user's own, so a ring that loses its entry shows (a default label may come back from the stored
  // label list, a colour does not).
  drawing.circularConservation.series.forEach((entry, index) => Object.assign(entry, { label: `G${'ABC'[index]}`, color: `#12345${index}` }));
  stackSources.refreshCircularConservationSlots(drawing);
  widenRingRows();
  return { rows: stackRingRows(), looks: ringLooks() };
};

test('after a LOSAT Generate, Add Seq keeps each ring\'s label, colour and stack row (#1777 review 6)', async () => {
  await withStack(async (drawing) => {
    const { rows, looks } = await threeSequenceRings(drawing);
    await generateRings(drawing);
    assert.deepEqual(stackRingRows(), rows, 'the Generate keeps the rows');
    assert.deepEqual(ringIdentityDefects(drawing), [], 'one identity per ring after the Generate');
    await stackFiles.addCircularConservationComparisonFile(/** @type {any} */ ({ target: { files: [sequenceFastaOf('CCCC')], value: 'x' } }));
    assert.deepEqual(ringLooks().slice(0, 3), looks, 'labels and colours after Add Seq');
    assert.deepEqual(stackRingRows().slice(0, 3), rows, 'rows after Add Seq');
    assert.deepEqual(await drawnRings(), [['GA', 'ACGT'], ['GB', 'TTTT'], ['GC', 'GGGG'], ['sequence', 'CCCC']]);
    assert.deepEqual(ringIdentityDefects(drawing), []);
  });
});

test('after a LOSAT Generate, a ring row\'s file edit keeps each ring\'s stack row (#1777 review 6)', async () => {
  await withStack(async (drawing) => {
    const { rows, looks } = await threeSequenceRings(drawing);
    await generateRings(drawing);
    stackFiles.addRingRowFiles(0, /** @type {any} */ ({ target: { files: [textFile('plasmid.fasta')], value: 'x' } }));
    assert.deepEqual(stackRingRows(), rows, 'Add file on the first ring row');
    assert.deepEqual(ringLooks(), looks);
    assert.deepEqual(ringIdentityDefects(drawing), []);
    stackFiles.removeRingRowFile(0, 1);
    assert.deepEqual(stackRingRows(), rows, 'its remove button');
    assert.deepEqual(await drawnRings(), [['GA', 'ACGT'], ['GB', 'TTTT'], ['GC', 'GGGG']]);
    assert.deepEqual(ringIdentityDefects(drawing), []);
  });
});

test('after a LOSAT Generate, Remove series keeps the other rings\' labels and stack rows (#1777 review 6)', async () => {
  await withStack(async (drawing) => {
    const { rows, looks } = await threeSequenceRings(drawing);
    await generateRings(drawing);
    leavingAsked.splice(0);
    stackFiles.removeCircularConservationSource(1);
    assert.deepEqual(leavingAsked, [[rows[1][0]]], 'the reach dialog takes out the removed ring\'s row');
    assert.deepEqual(stackRingRows(), [rows[0], rows[2]]);
    assert.deepEqual(ringLooks(), [looks[0], looks[2]]);
    assert.deepEqual(ringIdentityDefects(drawing), []);
    await stackFiles.addCircularConservationComparisonFile(/** @type {any} */ ({ target: { files: [sequenceFastaOf('CCCC')], value: 'x' } }));
    assert.deepEqual(stackRingRows().slice(0, 2), [rows[0], rows[2]], 'then Add Seq');
    assert.deepEqual(await drawnRings(), [['GA', 'ACGT'], ['GC', 'GGGG'], ['sequence', 'CCCC']]);
    assert.deepEqual(ringIdentityDefects(drawing), []);
  });
});

test('after a LOSAT Generate, the source radio upload, LOSAT, upload keeps each ring\'s label and stack row (#1777 review 6)', async () => {
  await withStack(async (drawing) => {
    const { rows, looks } = await threeSequenceRings(drawing);
    await generateRings(drawing);
    for (const source of /** @type {const} */ (['upload', 'losat', 'upload'])) {
      stackPanel.changeCircularConservationSource(source);
      assert.deepEqual(stackRingRows(), rows, `rows after the radio's ${source}`);
      assert.deepEqual(ringLooks(), looks, `labels and colours after the radio's ${source}`);
      assert.deepEqual(ringIdentityDefects(drawing), [], `one identity per ring after the radio's ${source}`);
    }
    stackPanel.changeCircularConservationSource('losat');
  });
});

// #1777 review 7 (OV-3194): a rerender of the committed request (R1 (c)) leaves the draft's moved rings as the move
// left them; the next Generate draws them in that order (tests/web/ring-draft-transitions.test.mjs walks every edit).
test('a rerender of the committed request after a move leaves each ring where the move put it, with its own label, colour and stack row (#1777 review 7)', async () => {
  for (const distinct of [false, true]) {
    await withStack(async (drawing) => {
      const { rows, looks } = await threeSequenceRings(drawing, { distinct });
      const canonical = await generateRings(drawing);
      stackPanel.moveCircularConservationSeries(0, 1);
      rerenderCommitted(drawing, canonical);
      const files = distinct ? 'files of their own names' : 'sequence.fasta files';
      assert.deepEqual(ringLooks(), [looks[1], looks[0], looks[2]], `${files}: the moved rings, each with its own label and colour`);
      assert.deepEqual(stackRingRows(), rows, `${files}: each ring in its own row`);
      assert.deepEqual(await drawnRings(), [['GB', 'TTTT'], ['GA', 'ACGT'], ['GC', 'GGGG']], files);
      assert.deepEqual(ringIdentityDefects(drawing), [], files);
      await generateRings(drawing);
      assert.deepEqual(ringLooks(), [looks[1], looks[0], looks[2]], `${files}: the next Generate`);
      assert.deepEqual(stackRingRows(), rows, `${files}: the next Generate keeps each row`);
      assert.deepEqual(await drawnRings(), [['GB', 'TTTT'], ['GA', 'ACGT'], ['GC', 'GGGG']], `${files}: the next Generate`);
    });
  }
});

// OV-3178: the request carries no ring's LOSAT genetic code. A rerender of the committed request leaves the draft's codes
// (OV-3194), and a Generate gives each ring the code it searched it with (OV-3195), never the draft's entry at the same
// index.
test('a rerender of the committed rings after a move keeps each ring\'s own genetic code for the next Generate (#1777 round 9, OV-3178)', async () => {
  for (const distinct of [false, true]) {
    await withStack(async (drawing) => {
      await threeSequenceRings(drawing, { distinct });
      drawing.circularConservation.series.forEach((entry, index) => { entry.losat_gencode = [11, 4, 2][index]; });
      const files = distinct ? 'files of their own names' : 'sequence.fasta files';
      const codes = () => drawing.circularConservation.series.map((entry) => [entry.label, entry.losat_gencode]);
      // The rings and codes the next Generate searches, as it reads them (`circularConservationKeys`).
      const searched = () => Promise.all(orderedConservationSources(ringGenomeFileLists(drawing.files.c_conservation_fastas), drawing.circularConservation)
        .map(async (source) => [(await readFileText(state.inputFiles.fileOf(source.file))).split('\n')[1], source.losat_gencode]));
      const canonical = await generateRings(drawing);
      assert.deepEqual(codes(), [['GA', 11], ['GB', 4], ['GC', 2]], `${files}: the Generate keeps each ring's code`);
      stackPanel.moveCircularConservationSeries(0, 1);
      assert.deepEqual(await searched(), [['TTTT', 4], ['ACGT', 11], ['GGGG', 2]], `${files}: the moved rings`);
      rerenderCommitted(drawing, canonical);
      assert.deepEqual(codes(), [['GB', 4], ['GA', 11], ['GC', 2]], `${files}: a rerender of the committed request after a move`);
      assert.deepEqual(await searched(), [['TTTT', 4], ['ACGT', 11], ['GGGG', 2]], `${files}: the next Generate searches each ring with its code`);
      stackPanel.moveCircularConservationSeries(1, 1);
      await generateRings(drawing);
      assert.deepEqual(codes(), [['GB', 4], ['GC', 2], ['GA', 11]], `${files}: a Generate after a move`);
      assert.deepEqual(ringIdentityDefects(drawing), [], files);
    });
  }
});

test('replacing an uploaded BLAST table for every drawing keeps its ring\'s label and colour (#1777 round 8)', async () => {
  const drawing = circular();
  const table = (name) => textFile(name, 'q\ts\t100\t8\t0\t0\t1\t8\t1\t8\t1e-9\t20\n');
  panel.changeCircularConservationSource('upload');
  panel.setCircularConservationUploadFiles([table('zebrafish.blast.tsv'), table('fly.blast.tsv')]);
  Object.assign(drawing.circularConservation.series[0], { label: 'Danio', color: '#123456' });
  // The replacement of an input file that the drawing uses (`createInputReplacement`) writes the table at its place.
  sources.replaceCircularConservationBlasts(drawing, [table('zebrafish.v2.blast.tsv'), drawing.files.c_conservation_blasts[1]]);
  assert.deepEqual(drawing.circularConservation.series.map((entry) => [entry.fileName, entry.label, entry.color]),
    [['zebrafish.v2.blast.tsv', 'Danio', '#123456'], ['fly.blast.tsv', 'fly.blast', drawing.circularConservation.series[1].color]]);
  assert.deepEqual(ringIdentityDefects(drawing), []);
  panel.setCircularConservationUploadFiles([]);
  panel.changeCircularConservationSource('losat');
});

test('a new ring is named after its file, not after a ring that went (#1777 round 8, OV-3166)', async () => {
  await withStack(async (drawing) => {
    await threeSequenceRings(drawing);
    await generateRings(drawing);
    await stackFiles.addCircularConservationComparisonFile(/** @type {any} */ ({ target: { files: [sequenceFastaOf('CCCC')], value: 'x' } }));
    while (drawing.circularConservation.series.length > 0) stackFiles.removeCircularConservationSource(0);
    await stackFiles.addCircularConservationComparisonFile(/** @type {any} */ ({ target: { files: [textFile('genome_x.fasta')], value: 'x' } }));
    assert.deepEqual(ringLooks().map(([label]) => label), ['genome_x']);
    assert.deepEqual(ringIdentityDefects(drawing), []);
  });
});

test('re-picked tables of other names are named after their files, not after the rings that went (#1777 round 8, OV-3166; D-144)', () => {
  const drawing = circular();
  const table = (name) => textFile(name, 'q\ts\t100\t8\t0\t0\t1\t8\t1\t8\t1e-9\t20\n');
  panel.changeCircularConservationSource('upload');
  panel.setCircularConservationUploadFiles([table('a.tsv'), table('b.tsv')]);
  drawing.circularConservation.series.forEach((entry, index) => { entry.label = `U${index}`; });
  comparisonFiles.removeCircularConservationSource(1);
  assert.equal(drawing.circularConservation.labels, '', 'the series owns the labels; no ring edit writes the stored label list (OV-3193)');
  // D-144: a re-picked set of other names gets rings of its own; the list's U0 is the ring that went, not c.tsv.
  panel.setCircularConservationUploadFiles([table('c.tsv'), table('d.tsv')]);
  assert.deepEqual(drawing.circularConservation.series.map((entry) => [entry.fileName, entry.label]), [['c.tsv', 'c'], ['d.tsv', 'd']]);
  panel.setCircularConservationUploadFiles([]);
  panel.changeCircularConservationSource('losat');
});

// PD-OI-110 (D-21): under Same color the re-assignment controls are disabled, and disabled has no effect. A move of
// the series arrows with 「トラックを並び替えたら色を振り直す」 on re-assigns under 10 colors, not under Same color, so
// the indexes are those before the move when 10 colors comes back.
test('a series move re-assigns the scheme colours with the checkbox on, except under Same color', async () => {
  const drawing = circular();
  for (const name of ['a.fasta', 'b.fasta', 'c.fasta']) await addSeq(name);
  const indexes = () => drawing.circularConservation.series.map(({ fileName, schemeIndex }) => [fileName, schemeIndex]);
  assert.deepEqual(indexes(), [['a.fasta', 0], ['b.fasta', 1], ['c.fasta', 2]]);
  panel.setCircularConservationColorScheme('ten');
  panel.setCircularConservationReassignOnReorder(true);
  panel.moveCircularConservationSeries(2, -1);
  assert.deepEqual(indexes(), [['a.fasta', 0], ['c.fasta', 1], ['b.fasta', 2]], '10 colors re-assigns in the new order');

  panel.setCircularConservationColorScheme('same');
  panel.moveCircularConservationSeries(2, -1);
  assert.deepEqual(ringNames(), ['a.fasta', 'b.fasta', 'c.fasta']);
  panel.setCircularConservationColorScheme('ten');
  assert.deepEqual(indexes(), [['a.fasta', 0], ['b.fasta', 2], ['c.fasta', 1]], 'Same color kept the indexes of before the move');
  panel.reassignCircularConservationColors();
  assert.deepEqual(indexes(), [['a.fasta', 0], ['b.fasta', 1], ['c.fasta', 2]]);
  panel.setCircularConservationColorScheme('same');
  panel.moveCircularConservationSeries(1, 1);
  panel.reassignCircularConservationColors();
  assert.deepEqual(indexes(), [['a.fasta', 0], ['c.fasta', 2], ['b.fasta', 1]], 'the button does nothing under Same color');

  panel.setCircularConservationReassignOnReorder(false);
  while (drawing.circularConservation.series.length > 0) comparisonFiles.removeCircularConservationSource(0);
});

// Fix round 13: a re-assignment gives every ring its drawn position as its scheme index, a ring with its own colour
// too (it keeps its colour), so no two rings share an index and 「スキームの色を使う」 then shows the colour of the
// ring's place, never another ring's.
test('a re-assignment leaves no two rings one scheme index, so Use the scheme color shows the ring\'s own place', async () => {
  const drawing = circular();
  const TEN = ['#81a5cc', '#e5a972', '#8fc4c0'];
  const drawnColors = () => orderedConservationSources(circularRingSources(drawing), drawing.circularConservation).map(({ color }) => color);
  const assertDistinct = (where) => {
    const indexes = drawing.circularConservation.series.map(({ schemeIndex }) => schemeIndex);
    assert.equal(new Set(indexes).size, indexes.length, `no two rings share a scheme index ${where}: ${JSON.stringify(indexes)}`);
  };
  for (const route of ['checkbox', 'button']) {
    for (const name of ['a.fasta', 'b.fasta', 'c.fasta']) await addSeq(name);
    panel.setCircularConservationColorScheme('ten');
    panel.setCircularConservationRingColor(0, '#123456');
    panel.setCircularConservationReassignOnReorder(route === 'checkbox');
    panel.moveCircularConservationSeries(2, -1);
    panel.moveCircularConservationSeries(1, -1);
    assert.deepEqual(ringNames(), ['c.fasta', 'a.fasta', 'b.fasta']);
    if (route === 'button') panel.reassignCircularConservationColors();
    assertDistinct(`after the ${route} re-assignment`);
    assert.deepEqual(drawnColors(), [TEN[0], '#123456', TEN[2]], `the hand-set colour stays (${route})`);
    panel.useCircularConservationSchemeColor(1);
    assertDistinct(`after Use the scheme color (${route})`);
    assert.deepEqual(drawnColors(), TEN, `each ring draws the colour of its drawn position (${route})`);
    panel.setCircularConservationReassignOnReorder(false);
    panel.setCircularConservationColorScheme('same');
    while (drawing.circularConservation.series.length > 0) comparisonFiles.removeCircularConservationSource(0);
  }
});
