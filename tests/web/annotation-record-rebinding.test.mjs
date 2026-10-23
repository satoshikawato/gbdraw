// Annotation record rebinding at the Linear discovery end (R10, D-15-46): the
// composition root (app/app-setup.js) rebinds the shown drawing's annotations
// to the records the catalog resolves in each transition that changes the
// catalog. A discovery that expands a multi-record file gives records 2..n
// new card keys, so the rebind runs after the expansion; a refresh that
// leaves the catalog as it was (after Save, at the end of Load) rebinds
// nothing, as the former catalog watcher did.
import assert from 'node:assert/strict';
import { File } from 'node:buffer';
import test from 'node:test';
import { installFakeSvgDom } from './fake-svg-dom.mjs';
import { modeDrawing } from './helpers/drawing-state.mjs';

globalThis.window = {
  Vue: {
    ref: (value) => ({ value }), shallowRef: (value) => ({ value }), reactive: (value) => value,
    isRef: (value) => Boolean(value && typeof value === 'object' && 'value' in value),
    computed: (getter) => ({ get value() { return getter(); } }), nextTick: async () => {},
    watch: () => () => {}, onMounted: () => {}, onUnmounted: () => {}, provide: () => {}
  },
  DOMPurify: { sanitize: (value) => value }
};
globalThis.document = { createElement: () => ({ click: () => {} }), querySelector: () => null };
installFakeSvgDom();
globalThis.File = File;

const { state } = await import('../../gbdraw/web/js/state.js');
const { createAppSetup } = await import('../../gbdraw/web/js/app/app-setup.js');
const { createAnnotationSet } = await import('../../gbdraw/web/js/services/annotation-state.js');
const { coordinateTarget } = await import('../../gbdraw/web/js/app/annotations/target-actions.js');
const { annotationRecordBinding } = await import('../../gbdraw/web/js/app/annotations/record-selector.js');
const { addAnnotationSet, addAnnotations } = await import('../../gbdraw/web/js/services/drawing/annotation-sets.js');

const genbank = (...records) => new File([records.map(([id, length]) => (
  `LOCUS       ${id} ${length} bp DNA\nACCESSION   ${id}\nORIGIN\n        1 ${'a'.repeat(length)}\n//\n`
)).join('')], 'records.gb');
const band = (id, recordId) => ({ id, target: coordinateTarget({ start: 1, end: 5, recordId }), mark: 'band' });

test('a Linear multi-record expansion binds an annotation to its record card, and a no-op refresh binds nothing', async () => {
  const app = createAppSetup();
  app.selectDrawing('linear');
  const drawing = modeDrawing(state, 'linear');
  app.setLinearSeqPrimaryFile(0, 'gb', genbank(['chr0', 40]));
  await app.refreshLinearRecordSelectors();
  addAnnotationSet(drawing.slices.annotationSets, createAnnotationSet({ id: 'set', annotations: [band('on-chr2', 'chr2')] }));
  app.addLinearSeq();
  app.setLinearSeqPrimaryFile(1, 'gb', genbank(['chr1', 50], ['chr2', 60]));
  await app.refreshLinearRecordSelectors();
  const card = state.linearSeqs.find((seq) => seq.region_record_id === 'chr2');
  assert.ok(card, 'the file expands into one card per record');
  const [set] = drawing.annotationSets;
  const bound = annotationRecordBinding(set.annotations[0]);
  assert.ok(bound.includes(JSON.stringify(card.uid)), `the annotation binds to the card of chr2, not ${bound}`);

  // The catalog is as the last rebind saw it: a refresh rebinds nothing.
  addAnnotations(drawing.slices.annotationSets, set, [band('on-chr0', 'chr0')]);
  await app.refreshLinearRecordSelectors();
  assert.equal(annotationRecordBinding(drawing.annotationSets[0].annotations[1]), '');
  app.selectDrawing('circular');
});
