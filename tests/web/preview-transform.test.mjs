// OV-912: a wheel tick or a zoom button changes the preview zoom. The app
// template does not read the zoom, the pan, or the panning flag, so these
// render only the zoom label and write only the surface transform; History
// and the Session still read `state.zoom` and `state.canvasPan`.
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import test from 'node:test';
import { loadPageVue } from './helpers/vue-runtime.mjs';

const Vue = loadPageVue();
globalThis.window = { Vue, setTimeout, clearTimeout };
const { nextTick } = Vue;
const { state } = await import('../../gbdraw/web/js/state.js');
const { previewTransform, PreviewZoomPercent } = await import('../../gbdraw/web/js/components.js');

test('the app template reads no preview zoom, pan, or panning flag (OV-912)', () => {
  const html = readFileSync(new URL('../../gbdraw/web/index.html', import.meta.url), 'utf8');
  const surface = html.slice(html.lastIndexOf('<div', html.indexOf('ref="svgContainer"')), html.indexOf('>', html.indexOf('ref="svgContainer"')) + 1);
  assert.match(surface, /\sv-preview-transform[\s>]/);
  assert.doesNotMatch(surface, /:style=/);
  // Outside event handlers, no binding reads them.
  const bindings = [...html.matchAll(/(?:\s:[\w.-]+|\sv-(?:if|show|for|bind|html|text|model)(?::[\w.-]+)?)="([^"]*)"|\{\{([^}]*)\}\}/g)]
    // String literals in an expression are text, not reads.
    .map((match) => (match[1] ?? match[2]).replace(/'[^']*'|`[^`$]*`/g, "''"));
  for (const name of ['zoom', 'canvasPan', 'isPanning']) {
    const readers = bindings.filter((expression) => new RegExp(`(^|[^\\w.])${name}\\b`).test(expression));
    assert.deepEqual(readers, [], `${name} is read by ${readers.join(' | ')}`);
  }
  assert.match(html, /<preview-zoom-percent><\/preview-zoom-percent>/);
});

test('the surface follows the zoom, the pan, and panning until it unmounts (OV-912)', async () => {
  const element = { style: {} };
  state.zoom.value = 1;
  state.canvasPan.x = 0;
  state.canvasPan.y = 0;
  state.isPanning.value = false;
  previewTransform.beforeMount(element);
  assert.deepEqual({ ...element.style }, {
    transform: 'translate(0px, 0px) scale(1)', transformOrigin: 'top center', transition: 'transform 0.2s'
  });
  state.zoom.value = 1.5;
  state.canvasPan.x = 12;
  await nextTick();
  assert.equal(element.style.transform, 'translate(12px, 0px) scale(1.5)');
  state.isPanning.value = true;
  await nextTick();
  assert.equal(element.style.transition, 'none');
  state.isPanning.value = false;
  previewTransform.unmounted(element);
  state.zoom.value = 2;
  await nextTick();
  assert.equal(element.style.transform, 'translate(12px, 0px) scale(1.5)', 'an unmounted surface is not written');
  state.zoom.value = 1;
  state.canvasPan.x = 0;
  assert.equal(PreviewZoomPercent.template.includes('zoom'), true);
  assert.equal(PreviewZoomPercent.setup().zoom, state.zoom);
});
