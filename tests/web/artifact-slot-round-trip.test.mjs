// A mode switch keeps the displayed artifact under the drawing it belongs to:
// the drawing of its mode (`generatedMode`), which may differ from the drawing
// shown. A Load of a Session with no Result whose request mode differs from
// `ui.mode` leaves that state; a round trip through the other mode must give
// each drawing its own artifact back (A2 I1 review, Medium 1). This checks the
// shown artifact and the slots after each switch, not what Save writes.
import assert from 'node:assert/strict';
import { File } from 'node:buffer';
import test from 'node:test';
import { installFakeSvgDom } from './fake-svg-dom.mjs';

globalThis.window = {
  Vue: {
    ref: (value) => ({ value }), shallowRef: (value) => ({ value }), reactive: (value) => value,
    isRef: (value) => Boolean(value && typeof value === 'object' && 'value' in value),
    computed: (getter) => ({ get value() { return getter(); } }), nextTick: async () => {},
    watch: () => () => {}, onMounted: () => {}, onUnmounted: () => {}, provide: () => {}
  },
  DOMPurify: { sanitize: (value) => value }
};
globalThis.document = { createElement: () => ({ click: () => {} }) };
installFakeSvgDom();
globalThis.File = File;

const { state } = await import('../../gbdraw/web/js/state.js');
const { createAppSetup } = await import('../../gbdraw/web/js/app/app-setup.js');

test('a round trip gives the artifact of the other mode back to its own drawing', async () => {
  const app = createAppSetup();
  // As such a Load leaves it: Circular shown, its displayed artifact the Linear one, with no Result.
  state.generatedMode.value = 'linear';
  state.generatedLegendPosition.value = 'top';

  app.setDiagramMode('linear');
  assert.equal(state.generatedMode.value, 'linear', 'the Linear drawing shows its own artifact');
  assert.equal(state.generatedLegendPosition.value, 'top');
  // A Linear Generate.
  const linearResult = { name: 'linear.svg', content: '<svg id="linear"/>' };
  state.results.value = [linearResult];

  app.setDiagramMode('circular');
  assert.equal(state.generatedMode.value, 'circular', 'the Circular drawing comes back with its own empty artifact');
  assert.deepEqual(state.results.value, []);

  app.setDiagramMode('linear');
  assert.deepEqual(state.results.value, [linearResult], 'the Linear drawing keeps its Result');
});
