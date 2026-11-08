// OV-2818 (Owner 2026-11-09 「Q3=A」): a Circular drawing's Show Comparison
// Rings follows its ring files as Show Depth follows the Depth files: the
// first ring file turns it on, adding a file turns it on, removing the last
// turns it off, the source radio follows the files it switches to, and Reset
// Settings turns it off while the ring rows and their files stay. Session Load
// reads a Session 51's stored setting and turns an older Session's rings on by
// its ring files.
import assert from 'node:assert/strict';
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
const { bindDrawingInputFileList } = await import('../../gbdraw/web/js/services/drawing/inputs.js');
const { sessionSourceFactsOf } = await import('../../gbdraw/web/js/services/session-drawings-upgrade.js');

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
  bindDrawingInputFileList(state.inputFiles, drawing.slices.inputs, 'c_conservation_fastas', []);
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
  assert.deepEqual(drawing.files.c_conservation_fastas.map((file) => file.name), ['zebrafish.fasta', 'fly.fasta']);
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
  bindDrawingInputFileList(state.inputFiles, drawing.slices.inputs, 'c_conservation_fastas', []);
  sources.reconcileLoadedCircularConservation(drawing, false);
  assert.equal(drawing.circularConservation.enabled, false, 'an older Session without ring files has none');
});
