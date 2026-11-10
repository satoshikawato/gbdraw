// Issue #1427 (PD-OI-114 scenario revision 2): the Web sends the Legend row edits and the canvas padding to Python
// as config leaves, from one projection (`legendEditConfigLeaves`); a Result Python drew with padding is not padded
// again (PD-OI-064); a Legend whose rows are all hidden takes no room live, as Python draws it (D-36).
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import test from 'node:test';

import { installFakeSvgDom } from './fake-svg-dom.mjs';
import { savedDrawingView } from './helpers/saved-drawing-view.mjs';
import { gallerySessionDocument } from './helpers/session-documents.mjs';

installFakeSvgDom();
globalThis.CSS ||= { escape: (value) => String(value) };

const { legendEditConfigLeaves, pythonLegendRowKeys, readLegendEditLeaves } = await import('../../gbdraw/web/js/services/legend-draft-rows.js');
const { managedConfigOverridePathsForMode, requestLegendEditLeaves } = await import('../../gbdraw/web/js/services/session-request.js');
const { applyCanvasToSvg } = await import('../../gbdraw/web/js/app/legend-layout/composition-actions.js');
const { createLegendLayoutActions } = await import('../../gbdraw/web/js/app/legend/layout-actions.js');
const { CANVAS_SLICE, setCanvasPaddingSide } = await import('../../gbdraw/web/js/services/drawing/canvas.js');
const { readCanvasPadding } = await import('../../gbdraw/web/js/services/config-field-readers.js');

const row = (key, caption = key) => ({ caption, originalCaption: key, color: '#cccccc', xPos: 0, yPos: 0, featureIds: [] });
const NO_PADDING = { top: 0, right: 0, bottom: 0, left: 0 };
const leavesOf = (draft) => legendEditConfigLeaves({
  legendEntries: [], deletedLegendEntries: [], originalLegendOrder: [], rules: [], canvasPadding: NO_PADDING,
  pythonRowKeys: draft.originalLegendOrder ?? [], ...draft
});

test('a drawing without Legend edits or padding writes no leaf', () => {
  assert.deepEqual(leavesOf({ legendEntries: [row('CDS'), row('tRNA')], originalLegendOrder: ['CDS', 'tRNA'] }), {});
  assert.deepEqual(leavesOf({}), {});
});

test('a rename of a row Python draws is a row_names leaf by its key; an editor row is not', () => {
  const leaves = leavesOf({
    legendEntries: [row('CDS', 'Genes'), row('tRNA'), row('My row', 'Renamed editor row')],
    originalLegendOrder: ['CDS', 'tRNA']
  });
  assert.deepEqual(leaves, { 'objects.legends.row_names': [['CDS', 'Genes']] });
});

test('a rule row is not renamed through the leaves (its rule caption carries the name)', () => {
  const leaves = leavesOf({
    legendEntries: [row('Kinase', 'Kinases'), row('CDS', 'Genes')],
    originalLegendOrder: ['Kinase', 'CDS'],
    rules: [{ feat: 'CDS', qual: 'product', val: 'kinase', color: '#ff0000', cap: 'Kinase' }]
  });
  assert.deepEqual(leaves, { 'objects.legends.row_names': [['CDS', 'Genes']] });
});

test('a deleted row Python draws is a hidden_rows leaf; one the inventory lacks is not', () => {
  const leaves = leavesOf({
    legendEntries: [row('CDS')],
    deletedLegendEntries: [row('tRNA'), row('Editor row')],
    originalLegendOrder: ['CDS', 'tRNA']
  });
  assert.deepEqual(leaves, { 'objects.legends.hidden_rows': ['tRNA'] });
});

test('an edited order is a row_order leaf of Python keys; the default order writes none', () => {
  const inventory = ['CDS', 'tRNA', 'rRNA'];
  assert.deepEqual(leavesOf({ legendEntries: [row('rRNA'), row('CDS'), row('tRNA')], originalLegendOrder: inventory }),
    { 'objects.legends.row_order': ['rRNA', 'CDS', 'tRNA'] });
  assert.deepEqual(leavesOf({ legendEntries: inventory.map((key) => row(key)), originalLegendOrder: inventory }), {});
});

test('canvas padding is a canvas.padding leaf (top, right, bottom, left) only when a side is not 0', () => {
  assert.deepEqual(leavesOf({ canvasPadding: { top: 10, right: 0, bottom: 30, left: 40 } }),
    { 'canvas.padding': [10, 0, 30, 40] });
  // R7 (review S2, D-59): the projection never converts; the slice's writer keeps only drawable sides.
  assert.deepEqual(leavesOf({ canvasPadding: { top: -5, right: 12.5, bottom: 0, left: 0 } }),
    { 'canvas.padding': [-5, 12.5, 0, 0] });
});

// Review S2 (D-59): the Canvas padding writer stores a side the canvas can draw (a negative or non-finite side is
// 0), Load installs through the same writer, and the request carries the slice as it is.
test('D-59: a negative padding side is stored and drawn as 0; Load installs a saved negative side as 0', () => {
  const { handle, view } = CANVAS_SLICE.create();
  setCanvasPaddingSide(handle, 'left', '-5');
  assert.deepEqual({ ...view.canvasPadding }, NO_PADDING, 'typing -5 stores 0');
  setCanvasPaddingSide(handle, 'left', 12.5);
  setCanvasPaddingSide(handle, 'top', '');
  assert.deepEqual({ ...view.canvasPadding }, { ...NO_PADDING, left: 12.5 }, 'a field without a number leaves the side');
  setCanvasPaddingSide(handle, 'left', -5);
  assert.equal(view.canvasPadding.left, 0);
  CANVAS_SLICE.install(handle, { canvasPadding: { top: -5, right: Number.POSITIVE_INFINITY, bottom: 7, left: 40 } });
  assert.deepEqual({ ...view.canvasPadding }, { top: 0, right: 0, bottom: 7, left: 40 }, 'Load and History install drawable sides');
  assert.deepEqual(leavesOf({ canvasPadding: view.canvasPadding }), { 'canvas.padding': [0, 0, 7, 40] });
  // Load reads a stored side the canvas cannot draw as 0 and lists it in the Load repair dialog (G2).
  assert.deepEqual(readCanvasPadding({ top: -7.5, right: 12.5, bottom: 0, left: 0 }),
    { value: { top: 0, right: 12.5, bottom: 0, left: 0 }, lost: true, leaves: [{ at: '.top', stored: -7.5, used: 0 }] });
  const index = readFileSync(new URL('../../gbdraw/web/index.html', import.meta.url), 'utf8');
  for (const side of ['top', 'right', 'bottom', 'left']) {
    assert.match(index, new RegExp(`<input id="canvas-padding-${side}" type="number" min="0" aria-describedby="canvas-padding-help"`),
      `the ${side} field offers no negative side, and the panel's help tip describes it (review N2, R12)`);
  }
});

test('the four leaves are GUI-managed paths, so a loaded request never keeps them as preserved overrides', () => {
  for (const mode of ['circular', 'linear']) {
    const managed = managedConfigOverridePathsForMode(mode);
    for (const path of ['objects.legends.row_names', 'objects.legends.hidden_rows', 'objects.legends.row_order', 'canvas.padding']) {
      assert.ok(managed.includes(path), `${mode}: ${path}`);
    }
  }
});

// A Gallery Result as the app mounts it (the fake DOM keeps attribute text as written).
const COMPOSITION = 'data-gbdraw-composition';
const decodeAttribute = (value) => value.replace(/&quot;/g, '"').replace(/&lt;/g, '<').replace(/&gt;/g, '>')
  .replace(/&amp;/g, '&');
const galleryDocuments = { HmmtDNA_basic_circular: await gallerySessionDocument('HmmtDNA_basic_circular') };
const galleryRoot = (base) => {
  const saved = galleryDocuments[base];
  const content = savedDrawingView(saved, saved.ui.activeDrawingId).results[0].content;
  const svg = new DOMParser().parseFromString(content, 'image/svg+xml').documentElement;
  svg.setAttribute(COMPOSITION, decodeAttribute(svg.getAttribute(COMPOSITION)));
  return svg;
};
// The Result as Python draws it with `canvas.padding`: its canvas grown by the padding, recorded in the metadata.
const drawnWithPadding = (svg, padding) => {
  const metadata = JSON.parse(svg.getAttribute(COMPOSITION));
  metadata.spacing.canvasPadding = padding;
  svg.setAttribute(COMPOSITION, JSON.stringify(metadata));
  const [x, y, width, height] = svg.getAttribute('viewBox').split(/\s+/).map(Number);
  svg.setAttribute('viewBox', [x, y, width + padding.left + padding.right, height + padding.top + padding.bottom].join(' '));
  return svg;
};

test('PD-OI-064: a Result Python drew with the drawing padding is not padded again', () => {
  const padding = { top: 10, right: 20, bottom: 30, left: 40 };
  const svg = drawnWithPadding(galleryRoot('HmmtDNA_basic_circular'), padding);
  const viewBox = svg.getAttribute('viewBox');
  assert.equal(applyCanvasToSvg(svg, padding, { refit: false }), false);
  assert.equal(svg.getAttribute('viewBox'), viewBox);
});

test('PD-OI-064: a padding changed after Python drew it adds only the difference', () => {
  const svg = drawnWithPadding(galleryRoot('HmmtDNA_basic_circular'), { top: 10, right: 0, bottom: 0, left: 0 });
  const [x, y, width, height] = svg.getAttribute('viewBox').split(/\s+/).map(Number);
  assert.equal(applyCanvasToSvg(svg, { top: 30, right: 0, bottom: 0, left: 5 }, { refit: false }), true);
  assert.deepEqual(svg.getAttribute('viewBox').split(/\s+/).map(Number), [x - 5, y - 20, width + 5, height + 20]);
});

test('PD-OI-064: a Result that records no padding is padded as before', () => {
  const svg = galleryRoot('HmmtDNA_basic_circular');
  const [x, y, width, height] = svg.getAttribute('viewBox').split(/\s+/).map(Number);
  assert.equal(applyCanvasToSvg(svg, { top: 10, right: 20, bottom: 30, left: 40 }, { refit: false }), true);
  assert.deepEqual(svg.getAttribute('viewBox').split(/\s+/).map(Number), [x - 40, y - 10, width + 60, height + 40]);
});

test('D-36: a Legend whose rows are all hidden lays out as an empty box, not null', async () => {
  const svg = galleryRoot('HmmtDNA_basic_circular');
  const rows = Array.from(svg.getElementById('legend').querySelectorAll('g[data-legend-key]'));
  assert.ok(rows.length > 0);
  rows.forEach((entry) => entry.setAttribute('display', 'none'));
  const actions = createLegendLayoutActions({ readRegisteredFontFaces: () => [] });
  await actions.prepareLegendLayout();
  assert.deepEqual(actions.layOutLegend(svg), { minX: 0, minY: 0, maxX: 0, maxY: 0 });
});

// Review S1 (D-58): the leaves do not depend on which Result of a Circular batch is shown. Result B draws tRNA,
// Result A does not; tRNA deleted while B was shown stays hidden in the request built while A is shown, from the
// rows Python drew in every Result. An editor row, which no Result's Python rows hold, stays the Web's.
test('D-58: a row deleted on one Result of a batch stays hidden in the request while another Result is shown', () => {
  const resultB = galleryRoot('HmmtDNA_basic_circular');
  const resultA = galleryRoot('HmmtDNA_basic_circular');
  // A draws no tRNA row: the editor's own row stands where B's tRNA row is.
  resultA.querySelectorAll('g[data-legend-key="tRNA"]').forEach((element) => {
    element.setAttribute('data-legend-key', 'Added row');
    element.setAttribute('data-legend-owner', 'direct-editor');
  });
  const contents = [resultA, resultB].map((svg) => new XMLSerializer().serializeToString(svg));
  const keysOfA = ['CDS', 'rRNA', 'GC content', 'GC skew (+)', 'GC skew (-)'];
  const state = { originalLegendOrder: { value: keysOfA }, results: { value: contents.map((content) => ({ content })) } };
  const drawing = {
    slices: {},
    legendEntries: { value: keysOfA.map((key) => row(key)) },
    deletedLegendEntries: { value: [row('tRNA'), row('Added row')] },
    manualSpecificRules: [],
    canvasPadding: NO_PADDING
  };
  assert.deepEqual(requestLegendEditLeaves(state, drawing), { 'objects.legends.hidden_rows': ['tRNA'] });
  // The keys Python drew in each Result, the editor's rows left out, with the shown Result's inventory.
  const results = contents.map((content) => ({ content }));
  assert.deepEqual(pythonLegendRowKeys(keysOfA, results).sort(), ['CDS', 'GC content', 'GC skew (+)', 'GC skew (-)', 'rRNA', 'tRNA']);
  // A Result is parsed once while its text stays: a parser that refuses to run proves the second call reads no SVG.
  const refuse = class { parseFromString() { throw new Error('parsed again'); } };
  assert.deepEqual(pythonLegendRowKeys(keysOfA, results, /** @type {typeof DOMParser} */ (/** @type {unknown} */ (refuse))).sort(),
    ['CDS', 'GC content', 'GC skew (+)', 'GC skew (-)', 'rRNA', 'tRNA']);
});

// Load reads the leaves of a stored request once (`readLegendEditLeaves`): a Web request's overrides, else a CLI or
// Python Session's full config; a leaf Python would refuse, and a padding of zeros, is no edit.
test('a stored request gives its Legend edit leaves from its overrides, else from its full config', () => {
  const leaves = {
    'objects.legends.row_names': [['CDS', 'Genes']],
    'objects.legends.hidden_rows': ['tRNA'],
    'objects.legends.row_order': ['rRNA', 'CDS'],
    'canvas.padding': [10, 20, 30, 40]
  };
  assert.deepEqual(readLegendEditLeaves({ diagramOptions: { configOverrides: leaves } }), leaves);
  const config = {
    canvas: { padding: [10, 20, 30, 40] },
    objects: { legends: { row_names: [['CDS', 'Genes']], hidden_rows: ['tRNA'], row_order: ['rRNA', 'CDS'] } }
  };
  assert.deepEqual(readLegendEditLeaves({ diagramOptions: { config } }), leaves);
  // The override wins over the full config.
  assert.deepEqual(readLegendEditLeaves({ diagramOptions: { config, configOverrides: { 'objects.legends.hidden_rows': ['rRNA'] } } })[
    'objects.legends.hidden_rows'], ['rRNA']);
  const none = { canvas: { padding: [0, 0, 0, 0] }, objects: { legends: { row_names: [], hidden_rows: null } } };
  assert.deepEqual(readLegendEditLeaves({ diagramOptions: { config: none } }), {});
  assert.deepEqual(readLegendEditLeaves({ diagramOptions: { configOverrides: {
    'objects.legends.row_names': [['CDS']], 'objects.legends.hidden_rows': [' '], 'canvas.padding': [1, 2, -3, 4]
  } } }), {});
  assert.deepEqual(readLegendEditLeaves(null), {});
});
