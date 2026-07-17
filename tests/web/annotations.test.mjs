import assert from 'node:assert/strict';
import { cp, mkdtemp, writeFile } from 'node:fs/promises';
import { tmpdir } from 'node:os';
import { join } from 'node:path';
import { pathToFileURL } from 'node:url';

const tempRoot = await mkdtemp(join(tmpdir(), 'gbdraw-annotations-'));
await cp(join(process.cwd(), 'gbdraw', 'web', 'js'), join(tempRoot, 'js'), { recursive: true });
await writeFile(join(tempRoot, 'package.json'), '{"type":"module"}', 'utf8');
const load = (path) => import(pathToFileURL(join(tempRoot, 'js', path)));

const { createAnnotationSet, normalizeAnnotationSets } = await load('app/annotations/state.js');
const { coordinateTarget, featureTarget } = await load('app/annotations/target-actions.js');
const { encodeAnnotationTable, parseAnnotationTable } = await load('app/annotations/table-codec.js');
const { createAnnotationEditor } = await load('app/annotations.js');
const { buildLinearTrackSlotSpec, normalizeLinearTrackSlots } = await load('app/linear-track-slots.js');
const { buildCircularTrackSlotSpec, normalizeCircularTrackSlots } = await load('app/circular-track-slots.js');

const set = createAnnotationSet({
  id: 'review',
  annotations: [
    { id: 'coords', target: coordinateTarget({ start: 10, end: 20 }), label: 'Window', mark: 'band' },
    { id: 'gene', target: featureTarget({ selector: 'locus_tag=ABC_1' }), label: 'Gene', mark: 'bracket' }
  ]
});
const normalized = normalizeAnnotationSets([set]);
assert.equal(normalized[0].annotations[0].target.coordinateSpace, 'source');
assert.deepEqual(normalized[0].annotations[1].target.selectors[0], { key: 'locus_tag', value: 'ABC_1' });

const table = encodeAnnotationTable(normalized);
const restored = parseAnnotationTable(table);
assert.equal(restored[0].id, 'review');
assert.equal(restored[0].annotations[0].target.start, 10);
assert.equal(restored[0].annotations[1].target.selectors[0].value, 'ABC_1');

const state = {
  annotationSets: [],
  selectedFeatures: { value: [] },
  adv: { circular_track_slots: [], linear_track_slots: [] }
};
const editor = createAnnotationEditor({ state });
const created = editor.addAnnotationSet('review');
editor.addCoordinateAnnotation(created, { start: 5, end: 8 });
const copied = editor.duplicateAnnotationSet(created);
assert.equal(state.annotationSets.length, 2);
assert.notEqual(copied.id, created.id);
editor.removeAnnotationSet(copied);
assert.equal(state.annotationSets.length, 1);

const linearSlot = normalizeLinearTrackSlots([{
  id: 'notes', renderer: 'annotations', side: 'overlay',
  params: { set_id: 'review', anchor_slot: 'features', layer: 'underlay', show_labels: false }
}])[0];
assert.equal(linearSlot.side, 'overlay');
assert.match(buildLinearTrackSlotSpec(linearSlot), /set_id=review/);
assert.match(buildLinearTrackSlotSpec(linearSlot), /anchor_slot=features/);

const circularSlot = normalizeCircularTrackSlots([{
  id: 'notes', renderer: 'annotations', side: 'outside',
  params: { set_id: 'review', overflow: 'compress' }
}])[0];
assert.equal(circularSlot.side, 'outside');
assert.match(buildCircularTrackSlotSpec(circularSlot), /set_id=review/);
assert.match(buildCircularTrackSlotSpec(circularSlot), /overflow=compress/);

