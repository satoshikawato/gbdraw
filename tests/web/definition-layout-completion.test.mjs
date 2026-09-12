import assert from 'node:assert/strict';
import { test } from 'node:test';
import { readFile } from 'node:fs/promises';
import * as transforms from '../../gbdraw/web/js/app/legend-layout/transform-utils.js';

// Load the existing owners with a controlled scheduler and Worker response.
// All production control flow remains intact; only imported collaborators vary.
const loadOwner = async (path, name, dependencies) => {
  const source = (await readFile(new URL(`../../gbdraw/web/js/app/${path}`, import.meta.url), 'utf8'))
    .replace(/import\s+\{([\s\S]*?)\}\s+from\s+['"][^'"]+['"];?/g, 'const {$1} = dependencies;')
    .replace(/export const /g, 'const ');
  return new Function('dependencies', 'setTimeout', 'clearTimeout', `${source}\nreturn ${name};`)(
    dependencies, dependencies.setTimeout, dependencies.clearTimeout);
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
    ...transforms, COMPOSITION_SCHEMA_ATTRIBUTE: 'schema', bindCompositionMetadata: () => binding,
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

test('controlled definition callbacks complete for the active root and ignore superseded Results', async () => {
  let callback;
  let resolveHelper;
  let mutations = 0;
  const root = { getAttribute: () => null, querySelectorAll: () => [], getElementById: () => null };
  let mounted = root;
  const state = { svgContent: ref('<svg/>'), svgResultIdentity: ref(1), mode: ref('circular'), shouldDeferCircularPreviewUpdates: ref(false),
    svgContainer: ref({ querySelector: () => mounted }), cInputType: ref('gb'), files: { c_gb: {} },
    linearSeqs: [], form: {}, adv: {}, selectedResultIndex: ref(0), results: ref([{ content: 'old' }]), skipCaptureBaseConfig: ref(false) };
  const create = await loadOwner('results.js', 'createResultsManager', {
    setTimeout(fn) { callback = fn; return 1; }, clearTimeout() {},
    isMultiRecordCanvasSvg: () => false, cloneFileBytesForTransfer: async () => new ArrayBuffer(0),
    DIAGRAM_HELPER_OPERATIONS: { REGENERATE_DEFINITION_SVGS: 'definitions' },
    runDiagramHelperOperation: () => new Promise(resolve => { resolveHelper = resolve; })
  });
  const owner = create({ state, legendLayout: { refreshCompositionGeometry() { mutations++; } } });
  owner.scheduleDefinitionUpdate();
  assert.equal(mutations, 0, 'fast path before the scheduled callback');
  callback();
  await new Promise(setImmediate);
  assert.equal(typeof resolveHelper, 'function');
  resolveHelper({ result: { definitions: [{ definition_group_id: 'unused' }] } });
  await new Promise(setImmediate);
  assert.equal(mutations, 1, 'the current root receives its legitimate delayed layout');
  owner.scheduleDefinitionUpdate();
  callback();
  await new Promise(setImmediate);
  mounted = { ...root };
  state.results.value = [{ content: 'replacement' }];
  state.svgResultIdentity.value = 2;
  resolveHelper({ result: { definitions: [{ definition_group_id: 'unused' }] } });
  await new Promise(setImmediate);
  assert.equal(mutations, 1, 'superseded callback must not reflow the current root');
  assert.equal(state.results.value[0].content, 'replacement');
});
