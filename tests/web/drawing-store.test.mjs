// The drawing store writers (state.js; DRAWINGS-B design §3): `drawingTabs`
// holds the tab order and names, `drawingList()` follows it, and each writer
// refuses an unknown or duplicate drawing ID. `removeDrawing` refuses the last
// drawing and the shown one and drops the removed drawing's artifact slot.
import assert from 'node:assert/strict';
import { test } from 'node:test';

globalThis.window = { Vue: {
  ref: (value) => ({ value }), reactive: (value) => value, toRaw: (value) => value,
  computed: (getter) => ({ get value() { return getter(); } }), nextTick: async () => {}
} };
const { state, createDrawingState } = await import('../../gbdraw/web/js/state.js');
const { createArtifactSlotStore } = await import('../../gbdraw/web/js/services/artifact-slot.js');

const ids = () => state.drawingList().map((drawing) => drawing.id);
const refusal = (reason) => (error) => error?.code === 'INPUT_INVALID' && error.context?.reason === reason;

test('the store writers keep drawingList() in tab order and refuse bad IDs', () => {
  const original = state.drawingList();
  const originalTabs = state.drawingTabs.value;
  try {
    const tabsBefore = state.drawingTabs.value;
    assert.ok(Object.isFrozen(tabsBefore));
    assert.deepEqual(ids(), ['circular', 'linear']);

    const second = createDrawingState('linear-2', 'linear');
    state.addDrawing(second, 'Linear 2', 1);
    assert.deepEqual(ids(), ['circular', 'linear-2', 'linear']);
    assert.notEqual(state.drawingTabs.value, tabsBefore, 'each write replaces the frozen tab list');
    assert.equal(state.drawingOf('linear-2'), second);
    assert.throws(() => state.addDrawing(createDrawingState('linear-2', 'linear'), 'again', 0),
      refusal('DUPLICATE_DRAWING'));

    state.moveDrawing('linear-2', 99);
    assert.deepEqual(ids(), ['circular', 'linear', 'linear-2']);
    state.moveDrawing('linear-2', 0);
    assert.deepEqual(ids(), ['linear-2', 'circular', 'linear']);
    assert.throws(() => state.moveDrawing('nope', 0), refusal('UNKNOWN_DRAWING'));

    state.renameDrawing('linear-2', 'Phages');
    assert.deepEqual(state.drawingTabs.value.find((tab) => tab.id === 'linear-2'), { id: 'linear-2', name: 'Phages' });
    assert.throws(() => state.renameDrawing('nope', 'x'), refusal('UNKNOWN_DRAWING'));

    const slots = createArtifactSlotStore();
    slots.writeValues(second, { appliedPaletteName: 'kept' });
    assert.equal(slots.kept().length, 1);
    assert.throws(() => state.removeDrawing(state.activeDrawingId.value, slots), refusal('SHOWN_DRAWING'));
    assert.throws(() => state.removeDrawing('nope', slots), refusal('UNKNOWN_DRAWING'));
    assert.equal(state.removeDrawing('linear-2', slots)?.drawingId, 'linear-2', 'it hands its slot to the caller');
    assert.deepEqual(ids(), ['circular', 'linear']);
    assert.equal(slots.kept().length, 0, 'the removed drawing loses its slot');
    assert.throws(() => state.drawingOf('linear-2'));

    state.replaceDrawings([original[0]], [{ id: 'circular', name: 'Only' }]);
    assert.deepEqual(ids(), ['circular']);
    assert.throws(() => state.removeDrawing('circular', slots), refusal('LAST_DRAWING'));
    assert.throws(() => state.replaceDrawings(original, [{ id: 'circular', name: 'a' }, { id: 'circular', name: 'b' }]),
      refusal('DUPLICATE_DRAWING'));
    assert.throws(() => state.replaceDrawings(original, [{ id: 'circular', name: 'a' }]), refusal('UNKNOWN_DRAWING'));
    assert.throws(() => state.replaceDrawings([original[0]], [{ id: 'linear', name: 'a' }]), refusal('UNKNOWN_DRAWING'));
    assert.deepEqual(ids(), ['circular'], 'a refused write changes nothing');
  } finally {
    state.replaceDrawings(original, originalTabs);
  }
  assert.deepEqual(ids(), ['circular', 'linear']);
});
