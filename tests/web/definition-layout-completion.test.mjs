import assert from 'node:assert/strict';
import { test } from 'node:test';
import { readFile } from 'node:fs/promises';
import { setClassToken } from '../../gbdraw/web/js/services/svg-serialization.js';
import * as transforms from '../../gbdraw/web/js/app/legend-layout/transform-utils.js';

// Load the existing owners with a controlled scheduler and Worker response.
// All production control flow remains intact; only imported collaborators vary.
const loadOwner = async (path, name, dependencies) => {
  const source = (await readFile(new URL(`../../gbdraw/web/js/app/${path}`, import.meta.url), 'utf8'))
    .replace(/import\s+\{([\s\S]*?)\}\s+from\s+['"][^'"]+['"];?/g, 'const {$1} = dependencies;')
    .replace(/export const /g, 'const ');
  return new Function('dependencies', 'setTimeout', 'clearTimeout', 'window', `${source}\nreturn ${name};`)(
    dependencies, dependencies.setTimeout, dependencies.clearTimeout, { Vue: { reactive: value => value } });
};
const ref = value => ({ value });
const group = (id, transform) => {
  const attrs = { id, transform };
  return { id, attrs, style: { removeProperty() {} }, getAttribute: key => attrs[key] ?? null,
    setAttribute: (key, value) => { attrs[key] = value; }, addEventListener() {}, removeEventListener() {} };
};

test('incremental rebind adopts the replacement scale geometry, retaining same-root user offsets', async () => {
  const bar = group('length_bar', 'translate(0,100)');
  let activeBar = bar;
  const binding = { metadata: { primary: { automaticTranslation: [0, 0] }, title: null }, primary: { targets: [] }, title: { targets: [] } };
  const root = { getAttribute: () => '1', getElementById: id => id === 'length_bar' ? activeBar : null,
    addEventListener() {}, removeEventListener() {} };
  const create = await loadOwner('legend-layout/diagram-drag.js', 'createDiagramDragActions', {
    ...transforms, setClassToken, COMPOSITION_SCHEMA_ATTRIBUTE: 'schema', bindCompositionMetadata: () => binding,
    compositionUserDeltas: () => ({ primary: [] })
  });
  const state = Object.fromEntries(['results','selectedResultIndex','diagramElements','diagramElementIds',
    'diagramElementOriginalTransforms','diagramDragging','lengthBarElement','lengthBarOriginalTransform',
    'plotTitleElement','plotTitleDragging','plotTitleAutoTransform','layoutRepositionMode','zoom','skipCaptureBaseConfig'].map(key => [key, ref(null)]));
  for (const key of ['diagramOffset','diagramDragStart','lengthBarUserOffset','plotTitleDragStart','plotTitleUserOffset']) state[key] = { x: 0, y: 0 };
  state.svgContainer = ref({ querySelector: () => root });
  const actions = create({ state });
  actions.setupDiagramDrag(false);
  state.lengthBarUserOffset.y = 7;
  bar.setAttribute('transform', 'translate(0,107)');
  actions.setupDiagramDrag(true);
  assert.equal(bar.getAttribute('transform'), 'translate(0,107)');
  // A renderer-produced replacement already carries its own authoritative
  // translation. Merely binding it cannot apply a prior root's baseline.
  activeBar = group('length_bar', 'translate(0,144.5)');
  state.lengthBarUserOffset.y = 0;
  actions.setupDiagramDrag(true);
  assert.equal(activeBar.getAttribute('transform'), 'translate(0,144.5)');
});
