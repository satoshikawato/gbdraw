// OV-2822: Reset Settings keeps a drawing's Depth files, and each kept file
// names its series again, as its upload named it; the series is not left
// without a label, which the next Generate draws as "Depth".
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
const { createDrawingReset } = await import('../../gbdraw/web/js/app/reset-settings.js');
const { bindDrawingInputFileTree, createLinearSeq } = await import('../../gbdraw/web/js/services/drawing/inputs.js');

const depthFile = (name) => {
  const text = 'record\tposition\tdepth\nr1\t1\t10\n';
  return { name, type: 'text/plain', size: text.length, lastModified: 0, data: btoa(text) };
};

const { resetDrawing } = createDrawingReset({
  state,
  circularTrackNewRenderer: { value: '' },
  linearTrackNewRenderer: { value: '' },
  rederiveCircularConservation: () => {},
  applyLinearSeqMutation: (drawing, sequences) => { drawing.linearSeqs.splice(0, drawing.linearSeqs.length, ...sequences); },
  invalidateLinearComparisonArtifacts: () => {},
  writeArtifactValues: () => {}
});

test('Reset names a kept Circular Depth file\'s series after the file', () => {
  const drawing = modeDrawing(state, 'circular');
  const kept = [state.inputFiles.registerFile(depthFile('reads-1kb.depth.tsv')), state.inputFiles.registerFile(depthFile('second.tsv'))];
  bindDrawingInputFileTree(state.inputFiles, drawing.slices.inputs, 'c_depth', [kept]);
  drawing.adv.depth_tracks.splice(0, drawing.adv.depth_tracks.length, { label: 'my title' }, { label: 'other' });

  resetDrawing(drawing);

  assert.deepEqual(drawing.adv.depth_tracks.map((track) => track.label), ['reads-1kb.depth', 'second']);
});

test('Reset names a kept Linear Depth file\'s series after the file', () => {
  const drawing = modeDrawing(state, 'linear');
  const seq = createLinearSeq();
  seq.depth = [state.inputFiles.registerFile(depthFile('linear-reads.tsv'))];
  drawing.linearSeqs.splice(0, drawing.linearSeqs.length, seq);
  drawing.adv.depth_tracks.splice(0, drawing.adv.depth_tracks.length, { label: 'my title' });

  resetDrawing(drawing);

  assert.deepEqual(drawing.adv.depth_tracks.map((track) => track.label), ['linear-reads']);
});
