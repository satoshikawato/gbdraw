// OV-282 (D-26): a Legend row of Python's holds no color but the user's
// Legend color; Python's swatch is the Result's (`originalLegendColors`), and
// the palette, the rules and the slots derive what a row shows. The writer
// (`listLegendRows` through `extractLegendEntries`), Save (`captureLegendPaint`),
// Load (`normalizeEditorStateData`) and the reader of a drawing that is not
// shown (`keptLegendRowColors`) keep to that one owner per value.
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import test from 'node:test';
import { File } from 'node:buffer';
import { installFakeSvgDom } from './fake-svg-dom.mjs';
import { installSessionImportWorker } from './helpers/session-import-node.mjs';
import { withDrawings } from './helpers/drawing-state.mjs';
import { Node, row, svgOf } from './helpers/legend-result-svg.mjs';

globalThis.window = {
  Vue: {
    ref: (value) => ({ value }), reactive: (value) => value,
    computed: (getter) => ({ get value() { return getter(); } }),
    nextTick: async () => {}
  },
  DOMPurify: { sanitize: (value) => value }
};
globalThis.document = { createElement: () => ({ click: () => {} }) };
installFakeSvgDom();
globalThis.File = File;
globalThis.CSS = { escape: (value) => String(value) };
installSessionImportWorker();

const { normalizeEditorStateData } = await import('../../gbdraw/web/js/services/config.js');
const { captureLegendPaint } = await import('../../gbdraw/web/js/services/drawing/paint.js');
const { createLegendEntryActions } = await import('../../gbdraw/web/js/app/legend/entry-actions.js');
const { displayedLegendRowColors, keptLegendRowColors } = await import('../../gbdraw/web/js/app/candidate-render.js');

const ref = (value) => ({ value });

const drawnResult = () => {
  const legend = new Node('g', { id: 'legend' }, [new Node('g', { id: 'feature_legend' }, [
    row('CDS', '#84b9ec', 0), row('tRNA', '#ddce76', 20), row('Mine', '#123456', 40, { 'data-legend-owner': 'direct-editor' })
  ])]);
  return new Node('svg', {}, [legend]);
};

test('a Generate lists the rows of Python\'s with no color but the user\'s, and Save writes no swatch', () => {
  const svg = drawnResult();
  const state = withDrawings({
    results: ref([{ name: 'r.svg', content: 'unchanged' }]), selectedResultIndex: ref(0),
    svgContainer: ref({ querySelector: (selector) => (selector === 'svg' ? svg : null) }),
    displayedResultMetadata: () => null, extractedFeatures: ref([]), adv: {},
    legendEntries: ref([{ caption: 'Mine', originalCaption: 'Mine', color: '#123456', featureIds: [] }]),
    deletedLegendEntries: ref([]), dormantLegendEntries: ref([]),
    originalLegendOrder: ref([]), originalLegendColors: ref({}),
    legendStrokeOverrides: {}, legendColorOverrides: { tRNA: '#ff8800' }, manualSpecificRules: [],
    skipCaptureBaseConfig: ref(false), svgContent: ref(null)
  });
  const actions = createLegendEntryActions({
    state, readLegendRowColors: (drawing) => displayedLegendRowColors(state, drawing, [])
  });
  actions.extractLegendEntries({ replaceGeneratedInventory: true });
  const saved = captureLegendPaint(state.activeDrawing().slices.paint).legendEntries;
  assert.deepEqual(saved.map((entry) => [entry.originalCaption, Object.hasOwn(entry, 'color') ? entry.color : null]), [
    ['CDS', null], ['tRNA', '#ff8800'], ['Mine', '#123456']
  ], 'the CDS row holds no swatch; the user-colored row its Legend color; the editor row its own');
  assert.deepEqual(state.originalLegendColors.value, { CDS: '#84b9ec', tRNA: '#ddce76' }, 'Python\'s swatches are the Result\'s');
  const shown = displayedLegendRowColors(state, state.activeDrawing(), []);
  assert.deepEqual(state.activeDrawing().legendEntries.value.map((entry) => shown.listed(entry)), ['#84b9ec', '#ff8800', '#123456']);
});

// OV-1725: each extraction writes Python's swatch of every row the Result
// draws, so a second Generate's swatches replace the first's, and a row with
// no derived color shows its swatch of the newest Result.
test('a second Generate\'s swatches replace the first\'s, and the rows show them', () => {
  let svg = svgOf([row('CDS', '#84b9ec', 0), row('tRNA', '#ddce76', 20)]);
  const state = withDrawings({
    results: ref([{ name: 'r.svg', content: 'unchanged' }]), selectedResultIndex: ref(0),
    svgContainer: ref({ querySelector: (selector) => (selector === 'svg' ? svg : null) }),
    displayedResultMetadata: () => null, extractedFeatures: ref([]), adv: {},
    legendEntries: ref([]), deletedLegendEntries: ref([]), dormantLegendEntries: ref([]),
    originalLegendOrder: ref([]), originalLegendColors: ref({}),
    legendStrokeOverrides: {}, legendColorOverrides: {}, manualSpecificRules: [],
    skipCaptureBaseConfig: ref(false), svgContent: ref(null)
  });
  const actions = createLegendEntryActions({ state, readLegendRowColors: (drawing) => displayedLegendRowColors(state, drawing, []) });
  actions.extractLegendEntries({ replaceGeneratedInventory: true });
  const first = state.originalLegendColors.value;
  svg = svgOf([row('CDS', '#2e7d32', 0), row('tRNA', '#ddce76', 20), row('My regions', '#abcdef', 40)]);
  actions.extractLegendEntries({ replaceGeneratedInventory: true });
  assert.deepEqual(state.originalLegendColors.value, { CDS: '#2e7d32', tRNA: '#ddce76', 'My regions': '#abcdef' });
  assert.deepEqual(first, { CDS: '#84b9ec', tRNA: '#ddce76' }, 'the first Result\'s colors, which a held figure keeps, stay as they were');
  const shown = displayedLegendRowColors(state, state.activeDrawing(), []);
  assert.deepEqual(state.activeDrawing().legendEntries.value.map((entry) => [entry.originalCaption, shown.listed(entry)]),
    [['CDS', '#2e7d32'], ['tRNA', '#ddce76'], ['My regions', '#abcdef']]);
  // A row the user colored shows the edit; its record keeps Python's swatch, which the extraction reads.
  svg = svgOf([row('CDS', '#ff0000', 0, {}, { 'data-gbdraw-base-fill': '#2e7d32' }), row('tRNA', '#ddce76', 20)]);
  state.activeDrawing().legendColorOverrides.CDS = '#ff0000';
  actions.extractLegendEntries();
  assert.equal(state.originalLegendColors.value.CDS, '#2e7d32');
  assert.equal(state.originalLegendColors.value['My regions'], '#abcdef', 'a row this Result does not draw keeps its swatch');
});

// A Session saved before this change holds the swatch of each row of
// Python's; Load reads only the colors a row owns.
const storedLegend = () => ({
  entries: [
    { caption: 'CDS', originalCaption: 'CDS', color: '#84b9ec', featureIds: [] },
    { caption: 'tRNA', originalCaption: 'tRNA', color: '#ff8800', featureIds: [] },
    { caption: 'Mine', originalCaption: 'Mine', color: '#123456', featureIds: [] },
    { caption: 'Plain', originalCaption: 'Plain', featureIds: [] }
  ],
  deletedEntries: [{ caption: 'rRNA', originalCaption: 'rRNA', color: '#7cecd5', featureIds: [] }],
  originalOrder: ['CDS', 'tRNA', 'rRNA', 'Plain'],
  originalColors: { CDS: '#84b9ec', tRNA: '#ddce76', rRNA: '#7cecd5', Plain: '#eeeeee' },
  colorOverrides: { tRNA: '#ff8800' }
});

test('Load keeps an entry\'s color only for an editor row and a user color, and reads an entry without one', () => {
  const repairs = [];
  const { legend } = normalizeEditorStateData({ legend: storedLegend() }, { repairs });
  const colors = (list) => list.map((entry) => [entry.originalCaption, Object.hasOwn(entry, 'color') ? entry.color : null]);
  assert.deepEqual(colors(legend.entries), [['CDS', null], ['tRNA', '#ff8800'], ['Mine', '#123456'], ['Plain', null]]);
  assert.deepEqual(colors(legend.deletedEntries), [['rRNA', null]]);
  assert.deepEqual(repairs, [], 'a recorded swatch is no loss');
  // An editor row without a color it can draw is still dropped and listed.
  const lost = [];
  const read = normalizeEditorStateData({ legend: { ...storedLegend(), entries: [{ caption: 'Bad', originalCaption: 'Bad', color: 'url(x)' }] } },
    { repairs: lost });
  assert.deepEqual(read.legend.entries, []);
  assert.deepEqual(lost.map((repair) => [repair.code, repair.label]), [['legend-color-dropped', 'Bad']]);
  // A dormant row is a renamed row of Python's (OV-120): its copy is discarded too. A swatch copy Load cannot read
  // is no loss; a stored color the user's Legend color replaces, which Save writes again, is listed (G2).
  const kept = [];
  const copies = normalizeEditorStateData({ legend: { ...storedLegend(),
    entries: [{ caption: 'C', originalCaption: 'CDS', color: 'url(x)' }, { caption: 'Mine', originalCaption: 'Mine', color: 'url(x)' }],
    dormantEntries: [{ caption: 'Renamed', originalCaption: 'Gone', color: '#999999' }],
    colorOverrides: { Mine: '#00ff00' } } }, { repairs: kept });
  assert.deepEqual(colors(copies.legend.entries), [['CDS', null], ['Mine', '#00ff00']]);
  assert.deepEqual(colors(copies.legend.dormantEntries), [['Gone', null]]);
  assert.deepEqual(kept.map((repair) => [repair.code, repair.label]), [['legend-color-dropped', 'Mine']]);
});

test('a Session that recorded swatches shows the same Legend colors after Load', () => {
  const session = JSON.parse(readFileSync(new URL('../fixtures/sessions/HmmtDNA_ATskew.v40-schema5.json', import.meta.url), 'utf8'));
  const stored = session.editorState.legend;
  const { legend } = normalizeEditorStateData(session.editorState, { repairs: [], sessionVersion: session.version });
  const shown = keptLegendRowColors({ legendColorOverrides: legend.colorOverrides }, legend);
  assert.ok(legend.entries.length > 0);
  assert.equal(legend.entries.some((entry) => Object.hasOwn(entry, 'color')), false, 'no row of Python\'s holds its swatch');
  assert.deepEqual(legend.entries.map((entry) => shown.listed(entry)), stored.entries.map((entry) => entry.color));
});

// A drawing that is not shown reads its kept Result's swatches: after a
// palette switch and a Generate in another drawing, its rows show Python's
// swatch of its own Result, not a copy an earlier extraction recorded.
test('a drawing that is not shown shows its kept Result\'s swatches, the user\'s colors and its editor rows\' own', () => {
  const drawing = { legendColorOverrides: { tRNA: '#ff8800' } };
  const shown = keptLegendRowColors(drawing, { originalOrder: ['CDS', 'tRNA'], originalColors: { CDS: '#2e7d32', tRNA: '#ddce76' } });
  assert.equal(shown.listed({ caption: 'CDS', originalCaption: 'CDS', color: '#84b9ec' }), '#2e7d32', 'not a stale copy');
  assert.equal(shown.listed({ caption: 'CDS', originalCaption: 'CDS' }), '#2e7d32');
  assert.equal(shown.deleted({ caption: 'tRNA', originalCaption: 'tRNA' }), '#ff8800');
  assert.equal(shown.listed({ caption: 'Mine', originalCaption: 'Mine', color: '#123456' }), '#123456');
});
