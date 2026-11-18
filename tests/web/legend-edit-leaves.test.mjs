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
  legendEntries: [], deletedLegendEntries: [], originalLegendOrder: [], rules: [], results: { renderedFeaturesByResult: [[]], palette: {} }, canvasPadding: NO_PADDING,
  pythonRowKeys: draft.originalLegendOrder ?? [], ...draft
});

// The keys Python gives the rows of a rule for every type of the type color value over each drawn type set
// (tests/fixtures/legend_rule_keys.json, which tests/test_specific_color_captions.py checks against Python).
const RULE_KEY_CASES = JSON.parse(readFileSync(new URL('../fixtures/legend_rule_keys.json', import.meta.url), 'utf8')).cases;
const ruleKeyCase = (name) => {
  const found = RULE_KEY_CASES.find((each) => each.name === name);
  assert.ok(found, name);
  return found;
};
const TWO_TYPES = ruleKeyCase('a palette rule for every type over two drawn types, a key per type');
const ONE_TYPE = ruleKeyCase('a palette rule for every type over one drawn type, keyed by its caption');

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

// OV-3803: the rows of a rule of the type color value are keyed in the type
// colors and drawn types of the shown Result, as Python keys them
// (tests/fixtures/legend_rule_keys.json), so their renames are carried by the
// rules too: a rule for every type draws a row per drawn type, and a caption
// also written in a literal color is keyed in the type's color.
test('the rows of a rule of the type color value are rule rows by Python\'s keys (OV-3803)', () => {
  const drawn = { renderedFeaturesByResult: [[{ type: 'rRNA' }, { type: 'CDS' }]], palette: { rRNA: '#71ee7d', CDS: '#54bcf8', default: '#d3d3d3' } };
  for (const [shape, rules, keys] of [
    ['a rule for every type', [{ feat: '*', qual: 'product', val: 'x', color: 'palette', cap: 'Ribosomal' }],
      ['Ribosomal [#71ee7d]', 'Ribosomal [#54bcf8]']],
    ['a caption also written in a literal color', [
      { feat: 'rRNA', qual: 'product', val: 'x', color: 'palette', cap: 'small rRNA' },
      { feat: 'rRNA', qual: 'product', val: 'y', color: '#ff00ff', cap: 'small rRNA' }
    ], ['small rRNA [#71ee7d]', 'small rRNA [#ff00ff]']]
  ]) {
    const leaves = leavesOf({
      legendEntries: [...keys.map((key) => row(key, 'Renamed')), row('CDS', 'Genes')],
      originalLegendOrder: [...keys, 'CDS'],
      rules,
      results: drawn
    });
    assert.deepEqual(leaves, { 'objects.legends.row_names': [['CDS', 'Genes']] }, shape);
  }
});

// A batch's rule rows are keyed in each Result alone, as Python keys them
// (`gbdraw circular` of HmmtDNA, NC_001416 and MjeNMV with
// `* translation . palette Ribosomal`: `Ribosomal [#54bcf8]` in the Result that
// draws CDS, rRNA and tRNA; `Ribosomal` in the Results that draw CDS only;
// tests/fixtures/legend_rule_keys.json, a rule for every type over two drawn
// types and over one), and the leaves name the rule rows of any Result,
// whichever Result is shown.
test('the rule rows of a batch are the rows of each Result keyed alone, whichever Result is shown', () => {
  const results = {
    renderedFeaturesByResult: [TWO_TYPES.types.map((type) => ({ type })), ONE_TYPE.types.map((type) => ({ type }))],
    palette: TWO_TYPES.palette
  };
  const rules = [{ feat: '*', qual: 'translation', val: '.', color: 'palette', cap: 'Ribosomal' }];
  const keys = [...TWO_TYPES.keys[0], ...ONE_TYPE.keys[0]];
  for (const shown of [TWO_TYPES.keys[0], ONE_TYPE.keys[0]]) {
    const leaves = leavesOf({
      legendEntries: [...shown.map((key) => row(key, 'Renamed')), row('CDS', 'Genes')],
      originalLegendOrder: [...shown, 'CDS'],
      pythonRowKeys: [...keys, 'CDS'],
      rules,
      results
    });
    assert.deepEqual(leaves, { 'objects.legends.row_names': [['CDS', 'Genes']] }, `shown: ${shown.join(', ')}`);
  }
});

// OV-3805: a rule captioned like a renderer row of another swatch draws its
// own row, `<caption> [<fill>]` (tests/fixtures/legend_rule_keys.json), which
// the leaves read from Python's rows of the Result: the renderer row keeps its
// rename. HmmtDNA, rule `CDS product . #ff0000 "tRNA"`, the tRNA type row
// renamed "Transfer RNA".
test('a rename of a renderer row a rule caption names in another color is sent (OV-3805)', () => {
  const pythonRows = new Map([
    ['tRNA', { key: 'tRNA', color: '#e8b441' }], ['tRNA [#ff0000]', { key: 'tRNA [#ff0000]', color: '#ff0000' }]
  ]);
  const leaves = leavesOf({
    legendEntries: [row('tRNA [#ff0000]', 'tRNA'), row('tRNA', 'Transfer RNA')],
    originalLegendOrder: ['tRNA [#ff0000]', 'tRNA'],
    rules: [{ feat: 'CDS', qual: 'product', val: '.', color: '#ff0000', cap: 'tRNA' }],
    results: { renderedFeaturesByResult: [[{ type: 'CDS' }, { type: 'tRNA' }]], palette: { CDS: '#54bcf8', tRNA: '#e8b441' }, pythonRowsByResult: [pythonRows] }
  });
  assert.deepEqual(leaves, { 'objects.legends.row_names': [['tRNA', 'Transfer RNA']] });
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
const galleryResultText = (base) => {
  const saved = galleryDocuments[base];
  return savedDrawingView(saved, saved.ui.activeDrawingId).results[0].content;
};
const galleryRoot = (base) => {
  const content = galleryResultText(base);
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
    mode: 'circular',
    slices: {},
    legendEntries: { value: keysOfA.map((key) => row(key)) },
    deletedLegendEntries: { value: [row('tRNA'), row('Added row')] },
    dormantLegendEntries: { value: [] },
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

// A two-Result batch's catalog (Session 40+): Result 0 draws a feature of each of `TWO_TYPES.types`, Result 1 one of
// `ONE_TYPE.types`.
const catalogFeature = (recordKey, id, type, start) => ({
  recordKey, biologicalFeatureId: id, record_id: recordKey, type, start, end: start + 30, strand: 1,
  anchorProfile: { precision: 'exact', operator: 'single', partOrder: 'biological', strand: '+' }, qualifiers: { locus_tag: [id] }
});
const catalogItem = (resultIndex, recordKey, types) => {
  const biologicalFeatures = types.map((type, index) => catalogFeature(recordKey, `${recordKey}-${index}`, type, index * 100));
  return {
    resultIndex, resultName: `result-${resultIndex}.svg`, recordKeys: [recordKey], biologicalFeatures,
    features: biologicalFeatures.map(({ biologicalFeatureId }) => ({
      svgId: `svg-${biologicalFeatureId}`, recordKey, biologicalFeatureId, fillColor: '#000000',
      drawnSelector: { hash: `svg-${biologicalFeatureId}`, location: null, recordLocation: null }
    })),
    orthogroups: [], annotations: [], comparisonMatches: []
  };
};
const batchCatalog = () => ({ schema: 5, items: [catalogItem(0, 'REC1', TWO_TYPES.types), catalogItem(1, 'REC2', ONE_TYPE.types)] });
const ref = (value) => ({ value });
const batchState = (catalog, shown, extra = {}) => ({
  originalLegendOrder: ref([]), results: ref([0, 1].map((index) => ({ name: `result-${index}.svg`, content: '<svg/>' }))),
  selectedResultIndex: ref(shown), featureCatalog: catalog, extractedFeatures: ref([]), displayedResultMetadata: () => null,
  paletteDefinitions: ref({ default: TWO_TYPES.palette }), appliedPaletteName: ref('default'), appliedPaletteColors: ref({}),
  ...extra
});
const ruleDrawing = (rules, entries) => ({
  mode: 'circular', slices: {}, legendEntries: ref(entries), deletedLegendEntries: ref([]), dormantLegendEntries: ref([]),
  manualSpecificRules: rules, featureOverrides: {}, featureVisibilityManualRules: [], canvasPadding: NO_PADDING
});

// The Save and Generate leaves key a batch's rule rows in each Result's own features of the catalog
// (`requestLegendEditLeaves`), as Python keys them: every row of the rule is a rule row, whichever Result is shown.
test('the request\'s leaves key the rule rows of each Result of a batch catalog alone, whichever Result is shown', () => {
  const rules = [{ feat: '*', qual: 'translation', val: '.', color: 'palette', cap: 'Ribosomal' }];
  const keys = [...TWO_TYPES.keys[0], ...ONE_TYPE.keys[0]];
  for (const shown of [0, 1]) {
    // The shown Result's inventory lists every row, so each is Python's; the displayed Result's features are none.
    const state = batchState(ref(batchCatalog()), shown, { originalLegendOrder: ref([...keys, 'CDS']) });
    const drawing = ruleDrawing(rules, [...keys.map((key) => row(key, 'Renamed')), row('CDS', 'Genes')]);
    // The rule rows' renames are carried by the rule; the CDS type row's is sent.
    assert.deepEqual(requestLegendEditLeaves(state, drawing), { 'objects.legends.row_names': [['CDS', 'Genes']] }, `shown: Result ${shown}`);
  }
});

// Without a catalog (a Session 27-39 batch before its first Generate), the leaves key the rule rows in the shown
// Result alone: its drawn features and Python's rows of that Result. HmmtDNA batch, rule `CDS product . #ff0000
// "tRNA"`: Result 1 (shown) draws tRNA, so the rule's row there is `tRNA [#ff0000]`; Result 0 draws no tRNA, and
// Python keys the rule's row `tRNA`. The tRNA type row renamed on Result 1 is sent.
test('without a catalog the leaves read the shown Result\'s features with that Result\'s rows', () => {
  // The Result's text with its Legend row of `key` keyed `as`, its swatch drawn in `fill` (as Python writes it).
  const base = galleryResultText('HmmtDNA_basic_circular');
  const redrawn = (key, fill, as = key) => {
    const row = new RegExp(`<g data-legend-key="${key}"><path ([^>]*?)fill="#[0-9a-f]{6}"`);
    assert.match(base, row, key);
    return base.replace(row, `<g data-legend-key="${as}"><path $1fill="${fill}"`);
  };
  const contents = [redrawn('tRNA', '#ff0000'), redrawn('rRNA', '#ff0000', 'tRNA [#ff0000]')];
  const inventory = ['tRNA [#ff0000]', 'tRNA', 'CDS'];
  const state = batchState(null, 1, {
    originalLegendOrder: ref(inventory), results: ref(contents.map((content, index) => ({ name: `result-${index}.svg`, content }))),
    extractedFeatures: ref([{ svg_id: 'c1', type: 'CDS' }, { svg_id: 't1', type: 'tRNA' }])
  });
  const drawing = ruleDrawing([{ feat: 'CDS', qual: 'product', val: '.', color: '#ff0000', cap: 'tRNA' }],
    [row('tRNA [#ff0000]', 'tRNA'), row('tRNA', 'Transfer RNA'), row('CDS')]);
  assert.deepEqual(requestLegendEditLeaves(state, drawing), { 'objects.legends.row_names': [['tRNA', 'Transfer RNA']] });
});

// A request whose rules have no caption reads no Result's catalog features: the leaves admit the catalog only
// when a rule row's key needs them.
test('the request\'s leaves admit no catalog without a captioned rule', () => {
  // An admission reads the catalog's schema first (`admitFeatureCatalog`).
  let admissions = 0;
  const catalog = batchCatalog();
  const { schema } = catalog;
  Object.defineProperty(catalog, 'schema', { get: () => { admissions += 1; return schema; } });
  const drawing = ruleDrawing([{ feat: 'CDS', qual: 'product', val: '.', color: '#ff0000', cap: '' }], [row('CDS', 'Genes')]);
  const state = batchState(ref(catalog), 0, { originalLegendOrder: ref(['CDS']) });
  assert.deepEqual(requestLegendEditLeaves(state, drawing), { 'objects.legends.row_names': [['CDS', 'Genes']] });
  assert.equal(admissions, 0, 'no admission without a captioned rule');
  drawing.manualSpecificRules[0].cap = 'Red';
  requestLegendEditLeaves(state, drawing);
  assert.ok(admissions > 0, 'a captioned rule admits the catalog');
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
