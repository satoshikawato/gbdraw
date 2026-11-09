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
const { bindDrawingInputFileList, bindDrawingInputFileTree } = await import('../../gbdraw/web/js/services/drawing/inputs.js');
const { sessionSourceFactsOf } = await import('../../gbdraw/web/js/services/session-drawings-upgrade.js');
const { circularRingSources, orderedConservationSources } = await import('../../gbdraw/web/js/services/conservation-series.js');

const circular = () => modeDrawing(state, 'circular');
const textFile = (name, text = '>ring\nACGT\n') => ({ name, type: 'text/plain', size: text.length, lastModified: 0, data: btoa(text) });
const shownAction = (mode, action) => action(modeDrawing(state, mode));
const askingStackPorts = () => ({
  askChoice: async () => null, runUndoable: (_label, fn) => fn(), busy: () => false, rows: () => [], mode: 'circular'
});

const sources = createCircularConservationSources({
  state,
  writeLosatCacheInfo: () => {},
  syncCircularConservationSlots: () => {},
  circularConservationRowsLeaving: () => new Set(),
  askingStackPorts
});
const panel = createCircularConservationSeriesPanel({
  state, shownAction, shownDrawingOfMode: (mode) => modeDrawing(state, mode), writeLosatCacheInfo: () => {},
  circularConservationFastaInput: { value: null },
  reconcileCircularConservation: sources.reconcileCircularConservation,
  refreshCircularConservationSlots: sources.refreshCircularConservationSlots,
  askRingSourceRemoval: sources.askRingSourceRemoval
});
const comparisonFiles = createCircularConservationComparisonFiles({
  state, shownAction, writeLosatCacheInfo: () => {},
  // The step's commit waits for the record label read, which no Worker answers here.
  runUndoable: async () => undefined,
  clearDerivedCircularConservationBlasts: sources.clearDerivedCircularConservationBlasts,
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
