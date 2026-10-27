// PD-OI-089, OIC-029 (Owner 2026-10-09): the palette is the base layer and the user default
// colors (values that differ from the selected palette's color, as the `-d`
// table over `-p`) win over it. A palette switch keeps them, after a dialog when
// any exist; Default colors Reset asks before it discards them. A dialog
// choice is one History step, and Cancel records none (PD-OI-088).
import assert from 'node:assert/strict';
import { test } from 'node:test';
import { createResultsManager } from '../../gbdraw/web/js/app/results.js';
import { createPaletteLoader } from '../../gbdraw/web/js/app/palettes.js';
import { createHistoryManager } from '../../gbdraw/web/js/services/history.js';
import { createDialogChoice } from '../../gbdraw/web/js/app/history-inputs.js';
import { normalizePaletteColors, normalizePaletteDefinitions } from '../../gbdraw/web/js/utils/color-utils.js';
import { drawingStore } from './helpers/drawing-state.mjs';
import { PALETTE_NAME_ALIASES } from '../../gbdraw/web/js/contracts/schema.generated.js';
import {
  PALETTE_SLICE, clearQueuedPalette, queuePalette, replaceDefaultColors, restoreArtifactQueuedPalette, selectPalette, setDefaultColors
} from '../../gbdraw/web/js/services/drawing/palette.js';

const ref = (value) => ({ value });
const PALETTES = {
  default: { CDS: '#AABBCC', tRNA: '#111111', rRNA: '#121212', default: '#999999' },
  forest: { CDS: '#228B22', tRNA: '#222222', rRNA: '#232323', default: '#888888' }
};
const paletteOf = (name) => normalizePaletteColors({ ...PALETTES[name] });
// CDS `#abc` is the palette's `#AABBCC`, rRNA is Auto, and the comparison
// colors are the palette's defaults: three user colors (tRNA, misc_feature, gc_content).
const editedColors = () => ({
  ...paletteOf('default'), CDS: '#abc', tRNA: '#333333', rRNA: null, misc_feature: '#444444', gc_content: 'none'
});

// A drawing's palette slice, with its handle under `slices`.
const drawingOf = (members) => {
  const { handle, view } = PALETTE_SLICE.create(members);
  return { ...view, fileImports: {}, slices: { palette: handle } };
};

const setup = ({ colors = editedColors(), instant = true } = {}) => {
  const drawing = drawingOf({ currentColors: ref(colors) });
  const paletteColorsDialog = {
    show: false, kind: 'switch', fromPalette: '', toPalette: '', droppedPalette: '', count: 0, keysText: ''
  };
  const state = {
    paletteDefinitions: ref(PALETTES), paletteInstantPreviewEnabled: ref(instant), appliedPaletteName: ref('default'),
    appliedPaletteColors: ref({ ...colors }), activeDrawing: () => drawing, ...drawingStore(drawing, drawing), paletteColorsDialog
  };
  /** Whether the dialog was open at each History capture. */
  const openAtCapture = [];
  const history = createHistoryManager({
    buildIntent: () => (openAtCapture.push(paletteColorsDialog.show), {
      palette: drawing.selectedPalette.value, colors: { ...drawing.currentColors.value },
      applied: { ...state.appliedPaletteColors.value }, pending: drawing.pendingPaletteName.value
    }),
    signatureFor: (value) => JSON.stringify(value),
    applyIntent: () => {}, buildCheckpoint: () => ({}), applyCheckpoint: () => {}
  });
  // app-setup.js wires the dialog's choices the same way.
  const dialogChoice = createDialogChoice({ mutationPending: history.mutationPending, runUndoable: history.runUndoable, ref });
  const manager = createResultsManager({
    state, closeAfterDialogChoice: dialogChoice.closeAfterChoice, writeArtifactValues: () => { throw new Error('no hidden drawing'); }
  });
  const choose = dialogChoice.withHistory(
    () => 'Change setting', manager.handlePaletteColorsChoice, manager.cancelPaletteColorsDialog
  );
  return { drawing, state, manager, history, choose, paletteColorsDialog, pending: dialogChoice.pending, openAtCapture };
};

test('a user default color differs from the selected palette (#rgb = #rrggbb); Auto is none', () => {
  const { drawing, manager } = setup();
  assert.equal(manager.readUserDefaultColor(drawing, 'CDS'), null);
  assert.equal(manager.readUserDefaultColor(drawing, 'rRNA'), null);
  assert.equal(manager.readUserDefaultColor(drawing, 'pairwise_match'), null);
  assert.equal(manager.readUserDefaultColor(drawing, 'tRNA'), '#333333');
  assert.equal(manager.readUserDefaultColor(drawing, 'misc_feature'), '#444444');
  assert.equal(manager.readUserDefaultColor(drawing, 'gc_content'), 'none');
});

test('a palette switch with user colors opens the dialog and changes nothing until a choice', () => {
  const { drawing, state, manager, paletteColorsDialog } = setup();
  const before = JSON.stringify(drawing.currentColors.value);
  const select = { value: 'forest' };
  manager.requestPaletteChange({ target: select });
  assert.equal(select.value, 'default');
  assert.deepEqual({ ...paletteColorsDialog }, {
    show: true, kind: 'switch', fromPalette: 'default', toPalette: 'forest', droppedPalette: '', count: 3,
    keysText: 'tRNA, misc_feature, and gc_content'
  });
  assert.equal(drawing.selectedPalette.value, 'default');
  assert.equal(JSON.stringify(drawing.currentColors.value), before);
  assert.equal(state.appliedPaletteName.value, 'default');
});

test('Keep my N colors applies the new palette under the user colors as one History step', async () => {
  const { drawing, state, manager, history, choose, paletteColorsDialog } = setup();
  manager.selectPalette('forest');
  await choose('keep');
  const expected = { ...paletteOf('forest'), tRNA: '#333333', misc_feature: '#444444', gc_content: 'none' };
  assert.deepEqual(drawing.currentColors.value, expected);
  assert.equal(drawing.selectedPalette.value, 'forest');
  assert.equal(state.appliedPaletteName.value, 'forest');
  assert.deepEqual(state.appliedPaletteColors.value, expected);
  assert.equal(paletteColorsDialog.show, false);
  assert.equal(history.getUndoCount(), 1);
});

test("Use the palette's colors applies the new palette alone as one History step", async () => {
  const { drawing, manager, history, choose } = setup();
  manager.selectPalette('forest');
  await choose('palette');
  assert.deepEqual(drawing.currentColors.value, paletteOf('forest'));
  assert.equal(drawing.selectedPalette.value, 'forest');
  assert.equal(history.getUndoCount(), 1);
});

test('Cancel of the palette dialog changes nothing and records no History step', async () => {
  const { drawing, manager, history, choose, paletteColorsDialog } = setup();
  const before = JSON.stringify(drawing.currentColors.value);
  manager.selectPalette('forest');
  await choose('cancel');
  assert.equal(paletteColorsDialog.show, false);
  assert.equal(drawing.selectedPalette.value, 'default');
  assert.equal(JSON.stringify(drawing.currentColors.value), before);
  assert.equal(history.getUndoCount(), 0);
});

test('a palette switch without user colors applies at once, without a dialog', () => {
  const { drawing, manager, paletteColorsDialog } = setup({ colors: { ...paletteOf('default'), CDS: '#aabbcc', rRNA: null } });
  manager.requestPaletteChange({ target: { value: 'forest' } });
  assert.equal(paletteColorsDialog.show, false);
  assert.equal(drawing.selectedPalette.value, 'forest');
  assert.deepEqual(drawing.currentColors.value, paletteOf('forest'));
});

test('a queued switch (Instant Preview off) keeps the applied palette and queues the kept colors', async () => {
  const { drawing, state, manager, choose } = setup({ instant: false });
  const applied = JSON.stringify(state.appliedPaletteColors.value);
  manager.selectPalette('forest');
  await choose('keep');
  assert.equal(drawing.pendingPaletteName.value, 'forest');
  assert.equal(drawing.pendingPaletteColors.value.tRNA, '#333333');
  assert.equal(drawing.pendingPaletteColors.value.CDS, '#228B22');
  assert.equal(state.appliedPaletteName.value, 'default');
  assert.equal(JSON.stringify(state.appliedPaletteColors.value), applied);
});

// D-15-46 (M3): turning Instant Preview on applies the queued palette in the
// toggle's own transition, not in a watcher; turning it off queues nothing.
test('turning Instant Preview on applies the queued palette; turning it off keeps the applied one', async () => {
  const { drawing, state, manager, choose } = setup({ instant: false });
  manager.selectPalette('forest');
  await choose('keep');
  const colors = { ...drawing.currentColors.value };
  manager.setPaletteInstantPreview(true);
  assert.equal(state.paletteInstantPreviewEnabled.value, true);
  assert.equal(drawing.pendingPaletteName.value, '');
  assert.deepEqual(drawing.pendingPaletteColors.value, {});
  assert.equal(state.appliedPaletteName.value, 'forest');
  assert.deepEqual(state.appliedPaletteColors.value, colors);
  manager.setPaletteInstantPreview(false);
  assert.equal(state.paletteInstantPreviewEnabled.value, false);
  assert.equal(drawing.pendingPaletteName.value, '');
  assert.equal(state.appliedPaletteName.value, 'forest');
});

// Q5 A (Owner 2026-10-09): switching back to the applied palette while
// another is queued asks like any switch when user colors exist; the queue then
// clears and the colors apply live.
test('switching back to the applied palette while one is queued asks like any switch, then applies live', async () => {
  const { drawing, state, manager, choose, paletteColorsDialog, history } = setup({ instant: false });
  manager.selectPalette('forest');
  await choose('palette');
  drawing.currentColors.value = { ...drawing.currentColors.value, tRNA: '#555555' };
  const steps = history.getUndoCount();
  manager.selectPalette('default');
  assert.deepEqual({ ...paletteColorsDialog }, {
    show: true, kind: 'switch', fromPalette: 'forest', toPalette: 'default', droppedPalette: 'forest', count: 1,
    keysText: 'tRNA'
  });
  assert.equal(drawing.selectedPalette.value, 'forest');
  await choose('keep');
  const expected = { ...paletteOf('default'), tRNA: '#555555' };
  assert.equal(drawing.selectedPalette.value, 'default');
  assert.deepEqual(drawing.currentColors.value, expected);
  assert.equal(drawing.pendingPaletteName.value, '');
  assert.equal(state.appliedPaletteName.value, 'default');
  assert.deepEqual(state.appliedPaletteColors.value, expected);
  assert.equal(history.getUndoCount(), steps + 1);
});

test('switching back to the applied palette without user colors applies its colors at once', async () => {
  const { drawing, state, manager, choose, paletteColorsDialog } = setup({ instant: false });
  manager.selectPalette('forest');
  await choose('palette');
  manager.selectPalette('default');
  assert.equal(paletteColorsDialog.show, false);
  assert.deepEqual(drawing.currentColors.value, paletteOf('default'));
  assert.equal(drawing.pendingPaletteName.value, '');
  assert.deepEqual(state.appliedPaletteColors.value, paletteOf('default'));
});

test('Default colors Reset with user colors asks first; Reset discards them as one History step', async () => {
  const { drawing, manager, history, choose, paletteColorsDialog } = setup();
  manager.requestResetColors();
  assert.deepEqual({ ...paletteColorsDialog }, {
    show: true, kind: 'reset', fromPalette: 'default', toPalette: 'default', droppedPalette: '', count: 3,
    keysText: 'tRNA, misc_feature, and gc_content'
  });
  assert.equal(drawing.currentColors.value.tRNA, '#333333');
  await choose('cancel');
  assert.equal(drawing.currentColors.value.tRNA, '#333333');
  assert.equal(history.getUndoCount(), 0);
  manager.requestResetColors();
  await choose('reset');
  assert.deepEqual(drawing.currentColors.value, paletteOf('default'));
  assert.equal(paletteColorsDialog.show, false);
  assert.equal(history.getUndoCount(), 1);
});

test('Default colors Reset without user colors resets at once (Auto returns to the palette color)', () => {
  const { drawing, manager, paletteColorsDialog } = setup({ colors: { ...paletteOf('default'), rRNA: null } });
  manager.requestResetColors();
  assert.equal(paletteColorsDialog.show, false);
  assert.deepEqual(drawing.currentColors.value, paletteOf('default'));
});

test('the dialog names three keys and counts the rest', () => {
  const { manager, paletteColorsDialog } = setup({ colors: { ...editedColors(), CDS: '#000000', rRNA: '#010101' } });
  manager.selectPalette('forest');
  assert.equal(paletteColorsDialog.count, 5);
  assert.equal(paletteColorsDialog.keysText, 'CDS, tRNA, rRNA, and 2 more');
});

// A field edit in the Default colors list sets the color in the table it
// edits: the table stays the same object, so the list's grouping (which reads
// the table's keys) never moves the edited row to the other grid. A new key or
// a whole-table change gives a new table, which the grouping follows.
test('a Default colors field edit keeps the table object; a new key replaces it', () => {
  const { drawing, manager } = setup();
  const table = drawing.currentColors.value;
  manager.editDefaultColor('tRNA', 'none');
  assert.equal(drawing.currentColors.value, table);
  assert.equal(drawing.currentColors.value.tRNA, 'none');
  manager.editDefaultColor('tmRNA', '#010203');
  assert.notEqual(drawing.currentColors.value, table);
  assert.equal(drawing.currentColors.value.tmRNA, '#010203');
});

test('setDefaultColor writes the default color and applies it to the Result', () => {
  const { drawing, state, manager } = setup();
  manager.setDefaultColor(drawing, 'CDS', '#123456');
  assert.equal(drawing.currentColors.value.CDS, '#123456');
  assert.equal(state.appliedPaletteColors.value.CDS, '#123456');
  assert.equal(manager.readUserDefaultColor(drawing, 'CDS'), '#123456');
});

// Q1 B (Owner 2026-10-09): while a palette is queued, the color also
// reaches the shown Result now; a user color wins over any palette, so the next
// Generate draws the same color for the key.
test('setDefaultColor while a palette is queued writes the queued color and the shown Result color', async () => {
  const { drawing, state, manager, choose } = setup({ instant: false });
  manager.selectPalette('forest');
  await choose('palette');
  const applied = { ...state.appliedPaletteColors.value };
  manager.setDefaultColor(drawing, 'CDS', '#123456');
  assert.equal(drawing.pendingPaletteName.value, 'forest');
  assert.equal(drawing.pendingPaletteColors.value.CDS, '#123456');
  assert.equal(drawing.pendingPaletteColors.value.tRNA, '#222222');
  assert.equal(state.appliedPaletteName.value, 'default');
  assert.deepEqual(state.appliedPaletteColors.value, { ...applied, CDS: '#123456' });
});

// OV-262, OV-263: a Result restyle reads the applied colors as a live edit,
// the Generate commit, Load, and a History restore leave them: an Auto key is
// empty there. It shows the applied palette's color for the key, as Generate
// draws it, not the palette's `default` color. Edit-unify: the feature fills
// are the live compile's, given `appliedFeatureColors` (app-setup.js).
test('a restyle paints an Auto key with the applied palette color', async () => {
  const { compileDirectEditorMutationPlan } = await import('../../gbdraw/web/js/app/candidate-render.js');
  const { displayedFeatureAddressing } = await import('../../gbdraw/web/js/services/feature-override-identity.js');
  const { appliedFeatureColors } = await import('../../gbdraw/web/js/utils/color-utils.js');
  const feature = { id: 'f1', svg_id: 'f1', type: 'CDS' };
  const paletteColors = appliedFeatureColors({
    paletteDefinitions: ref(PALETTES), appliedPaletteName: ref('default'),
    appliedPaletteColors: ref({ ...paletteOf('default'), CDS: null })
  });
  const plan = compileDirectEditorMutationPlan({
    catalogAdmission: displayedFeatureAddressing([feature], ['diagram.svg'], 0), manualSpecificRules: [],
    livePreview: { drawnTracks: [], domains: ['featureFills'], paletteColors, drawnContext: null }
  });
  assert.deepEqual(plan.operationsByResult[0].featureFills, [{ renderedId: 'f1', color: '#AABBCC' }]);
});

// PD-OI-088, OIC-028: the dialog stays open and busy until its choice's
// History step ends, as the popup choice dialogs do (`createDialogChoice`).
test('the palette dialog closes after the History step of its choice', async () => {
  const { manager, choose, paletteColorsDialog, pending, openAtCapture } = setup();
  manager.selectPalette('forest');
  assert.equal(paletteColorsDialog.show, true);
  const choice = choose('keep');
  assert.equal(pending.value, true);
  openAtCapture.length = 0;
  await choice;
  assert.ok(openAtCapture.length > 0);
  assert.ok(openAtCapture.every(Boolean), 'open while the step records');
  assert.equal(paletteColorsDialog.show, false);
  assert.equal(pending.value, false);
});

// OV-276: the Legend panel shows each row the palette colors in the color the
// live compile gives its swatch (rows with a Legend color keep theirs). The
// color is derived: nothing writes it into the Legend rows, which are an input
// of the rule preparation (a write during a Generate made it stale).
// Edit-unify: the rows a rule draws are Python's (`draftLegendRowColors`).
const legendPanelSetup = async ({ rules = [], cdsDrawn = '#AABBCC' } = {}) => {
  const { draftLegendPanelColors, draftLegendRowColors } = await import('../../gbdraw/web/js/app/candidate-render.js');
  const entries = [
    { caption: 'CDS', color: '#AABBCC' }, { caption: 'other tRNAs', color: '#111111' },
    { caption: 'Manual', color: '#010101' }, { caption: 'rRNA', color: '#121212' }
  ];
  const originalLegendOrder = ['CDS', 'other tRNAs', 'rRNA'];
  const pythonRows = new Map([
    ['CDS', { key: 'CDS', color: cdsDrawn }], ['other tRNAs', { key: 'other tRNAs', color: '#111111' }],
    ['rRNA', { key: 'rRNA', color: '#121212' }]
  ]);
  const paletteColors = { ...paletteOf('default'), CDS: '#123456', tRNA: '#333333', rRNA: '#444444' };
  const panel = draftLegendPanelColors({
    legendEntries: entries, originalLegendOrder, legendColorOverrides: { rRNA: '#121212' }, rules, pythonRows, paletteColors, drawnTracks: []
  });
  const compiled = draftLegendRowColors({ rules, pythonRows, originalLegendOrder, paletteColors, drawnTracks: [] });
  return { panel, compiled, entries };
};

test('the Legend panel shows the palette color of the rows the compile colors', async () => {
  const { panel, compiled, entries } = await legendPanelSetup();
  assert.equal(compiled.colorOf('CDS'), '#123456');
  assert.deepEqual([...panel.listed], [['CDS', '#123456'], ['other tRNAs', '#333333']]);
  assert.deepEqual(entries[0], { caption: 'CDS', color: '#AABBCC' }, 'the derivation writes no Legend row');
});

// OV-282 (D-26): a deleted row the palette colors shows the palette's color
// in the Deleted items list, by its Python key; a deleted row with a Legend
// color of its own and an editor row keep the color they record.
test('the Deleted items list shows the palette color of a deleted palette row', async () => {
  const { draftLegendPanelColors } = await import('../../gbdraw/web/js/app/candidate-render.js');
  const panel = draftLegendPanelColors({
    legendEntries: [{ caption: 'CDS', color: '#AABBCC' }],
    deletedLegendEntries: [
      { caption: 'tRNAs', originalCaption: 'other tRNAs', color: '#111111' },
      { caption: 'rRNA', color: '#121212' }, { caption: 'Manual', color: '#010101' }
    ],
    originalLegendOrder: ['CDS', 'other tRNAs', 'rRNA'],
    legendColorOverrides: { rRNA: '#121212' }, rules: [], pythonRows: new Map(),
    paletteColors: { ...paletteOf('default'), CDS: '#123456', tRNA: '#333333', rRNA: '#444444' }, drawnTracks: []
  });
  assert.deepEqual([...panel.listed], [['CDS', '#123456']]);
  assert.deepEqual([...panel.deleted], [['other tRNAs', '#333333']]);
});

// OV-309: the default track rows take the palette color Python draws them
// in, under the request's dinucleotide and as one row or two
// (`_generated_legend_rows`); Depth and a row of another dinucleotide, which
// the request draws no track for, take no palette color (OV-484).
test('the compile gives every track row Python draws its palette color', async () => {
  const { draftLegendRowColors } = await import('../../gbdraw/web/js/app/candidate-render.js');
  const { drawnTracks } = await import('../../gbdraw/web/js/services/track-legend-rows.js');
  const paletteColors = { gc_content: '#a1a1a1', skew_high: '#b2b2b2', skew_low: '#c3c3c3' };
  const tracks = drawnTracks({ mode: 'circular', diagramOptions: { dinucleotide: 'AT', tracks: { circularTrackSlots: null } } });
  // The rows of one Result: Python draws a track's row once, or as (+) and (-).
  for (const [keys, colors] of [
    [['AT content', 'AT skew (+)', 'AT skew (-)', 'GC skew (+)', 'Depth'], ['#a1a1a1', '#b2b2b2', '#c3c3c3', null, null]],
    [['AT content (+)', 'AT content (-)', 'AT skew'], ['#a1a1a1', '#a1a1a1', '#b2b2b2']]
  ]) {
    const pythonRows = new Map(keys.map((key) => [key, { key, color: null }]));
    const { colorOf } = draftLegendRowColors({ rules: [], pythonRows, originalLegendOrder: keys, paletteColors, drawnTracks: tracks });
    assert.deepEqual(keys.map((key) => colorOf(key)), colors);
  }
});

// Rules of one caption in two colors draw a row per color, each under its own
// key (tests/fixtures/legend_rule_keys.json), so a rule's row has one color
// whichever rule a drawn feature takes (OV-347 was one live row for both).
test('rules of one caption in two colors color a row each', async () => {
  const { draftLegendRowColors } = await import('../../gbdraw/web/js/app/candidate-render.js');
  const psaA = { feat: 'CDS', qual: 'gene', val: 'psaA', color: '#ff0000', cap: 'photosystem I' };
  const psaB = { feat: 'CDS', qual: 'gene', val: 'psaB', color: '#00662c', cap: 'photosystem I' };
  const draft = draftLegendRowColors({
    rules: [psaA, psaB],
    pythonRows: new Map([['photosystem I [#00662c]', { key: 'photosystem I [#00662c]', color: '#00662c' }]]),
    features: [{ type: 'CDS', gene: 'psaA', fill_color: '#54bcf8' }, { type: 'CDS', gene: 'psaB', fill_color: '#00662c' }],
    originalLegendOrder: ['photosystem I [#00662c]', 'other proteins'],
    paletteColors: { CDS: '#54bcf8' }, drawnTracks: []
  });
  assert.deepEqual(draft.ruleCaptions, ['photosystem I [#ff0000]', 'photosystem I [#00662c]']);
  assert.deepEqual(draft.ruleCaptions.map(draft.colorOf), ['#ff0000', '#00662c']);
  assert.equal(draft.colorOf('photosystem I'), null);
});

// Review 2: a Legend row that a Specific color rule draws keeps the rule's
// color, also when its caption names a palette key (a rule captioned CDS
// whose color Python drew in the CDS row): the compile gives it the rule's
// color, and the panel does not report it.
test('a palette change leaves a Legend row that a rule draws', async () => {
  const rule = { feat: 'CDS', qual: 'product', val: '.', color: '#00aa00', cap: 'CDS' };
  const { panel, compiled } = await legendPanelSetup({ rules: [rule], cdsDrawn: '#00aa00' });
  assert.deepEqual([...panel.listed], [['other tRNAs', '#333333']]);
  assert.equal(compiled.colorOf('CDS'), '#00aa00');
});

// OV-527: a hand-edited Session may name a palette `constructor` or `__proto__`;
// the table lookup must read own keys only, so the drawing falls back to `default`.
test('a palette named like an inherited member falls back to the default palette', () => {
  for (const name of ['constructor', '__proto__', 'toString']) {
    const circular = drawingOf({ selectedPalette: ref(name) });
    const linear = drawingOf();
    const state = {
      paletteDefinitions: ref({}), paletteNames: ref([]), appliedPaletteName: ref(''), appliedPaletteColors: ref({}),
      ...drawingStore(circular, linear), activeDrawing: () => circular
    };
    assert.equal(createPaletteLoader({ state }).applyPalettes(PALETTES), true);
    assert.equal(circular.selectedPalette.value, 'default', name);
  }
});

// `herbal_garden` was `hearbal_garden` through 0.14.2 (PALETTE_NAME_ALIASES of
// gbdraw/io/colors.py). The palette slice stores the current name of every
// palette it is given (Load, History, a queued palette); a lookup by the former
// name finds the palette, and the list shows only the current name.
test('a former palette name is stored, looked up, and listed by its current name', () => {
  assert.deepEqual(PALETTE_NAME_ALIASES, { hearbal_garden: 'herbal_garden' });
  const drawing = drawingOf({});
  PALETTE_SLICE.install(drawing.slices.palette, { selectedPalette: 'hearbal_garden', pendingPaletteName: 'hearbal_garden' });
  assert.equal(drawing.selectedPalette.value, 'herbal_garden');
  assert.equal(drawing.pendingPaletteName.value, 'herbal_garden');
  queuePalette(drawing.slices.palette, 'hearbal_garden', { CDS: '#000000' });
  assert.equal(drawing.pendingPaletteName.value, 'herbal_garden');
  restoreArtifactQueuedPalette(drawing.slices.palette, 'hearbal_garden', { CDS: '#000000' });
  assert.equal(drawing.pendingPaletteName.value, 'herbal_garden');
  selectPalette(drawing.slices.palette, 'forest');
  assert.equal(drawing.selectedPalette.value, 'forest');

  const circular = drawingOf();
  const state = {
    paletteDefinitions: ref({}), paletteNames: ref([]), appliedPaletteName: ref(''), appliedPaletteColors: ref({}),
    ...drawingStore(circular, drawingOf()), activeDrawing: () => circular
  };
  const herbal = { CDS: '#addc7f', default: '#d3d3d3' };
  assert.equal(createPaletteLoader({ state }).applyPalettes({ ...PALETTES, herbal_garden: herbal }), true);
  assert.deepEqual(state.paletteNames.value, ['default', 'forest', 'herbal_garden']);
  assert.deepEqual(Object.keys(state.paletteDefinitions.value), ['default', 'forest', 'herbal_garden']);
  assert.equal(state.paletteDefinitions.value.hearbal_garden, state.paletteDefinitions.value.herbal_garden);
  assert.deepEqual(state.paletteDefinitions.value.herbal_garden, normalizePaletteColors(herbal));
});

// Load copies a saved `ui.appliedPaletteName` as written (services/config.js),
// so a 0.14.2 Result drawn with `hearbal_garden` is the applied palette under
// its former name: choosing `herbal_garden` is switching back to it, which
// applies live and drops the queued palette (Q5 A).
test('the current name of a Result drawn with a former palette name is the applied palette', () => {
  const { drawing, state, manager } = setup({ colors: paletteOf('default'), instant: false });
  const herbal = { CDS: '#addc7f', default: '#d3d3d3' };
  state.paletteDefinitions.value = normalizePaletteDefinitions({ ...PALETTES, herbal_garden: herbal });
  state.appliedPaletteName.value = 'hearbal_garden';
  queuePalette(drawing.slices.palette, 'forest', paletteOf('forest'));
  manager.selectPalette('herbal_garden');
  assert.equal(drawing.selectedPalette.value, 'herbal_garden');
  assert.equal(drawing.pendingPaletteName.value, '');
  assert.equal(state.appliedPaletteName.value, 'herbal_garden');
  assert.deepEqual(state.appliedPaletteColors.value, normalizePaletteColors(herbal));
});

// An artifact owner set holds its queued palette colors by reference, and a
// History restore gives the drawing that object (`restoreArtifactQueuedPalette`).
// Every palette writer replaces the colors object, so none writes the owner
// set's: a write into the frozen object would throw.
test('a queued palette restored by reference is never written through', () => {
  const owned = Object.freeze({ CDS: '#123456' });
  const drawing = drawingOf({});
  restoreArtifactQueuedPalette(drawing.slices.palette, 'forest', owned);
  assert.equal(drawing.pendingPaletteColors.value, owned);
  const writers = [
    () => queuePalette(drawing.slices.palette, 'default', { CDS: '#000000' }),
    () => restoreArtifactQueuedPalette(drawing.slices.palette, 'forest', owned),
    () => clearQueuedPalette(drawing.slices.palette),
    () => restoreArtifactQueuedPalette(drawing.slices.palette, 'forest', owned),
    () => selectPalette(drawing.slices.palette, 'forest', { CDS: '#111111' }),
    () => setDefaultColors(drawing.slices.palette, { tRNA: '#222222' }),
    () => replaceDefaultColors(drawing.slices.palette, { CDS: '#333333' }),
    () => PALETTE_SLICE.install(drawing.slices.palette, { pendingPaletteColors: { CDS: '#444444' } }),
    () => restoreArtifactQueuedPalette(drawing.slices.palette, 'forest', owned),
    () => PALETTE_SLICE.reset(drawing.slices.palette)
  ];
  writers.forEach((write) => write());
  assert.deepEqual(owned, { CDS: '#123456' });
  assert.notEqual(drawing.pendingPaletteColors.value, owned);
});

// M3.12: a Default colors edit syncs the palette draft itself (no watcher
// writes a drawing): live, it applies at once; with a queued palette, it joins
// the queue and the applied palette stays.
test('a Default colors edit applies live, or joins the queued palette, without a watcher', async () => {
  const live = setup({ instant: true });
  live.manager.editDefaultColor('tRNA', '#565656');
  assert.equal(live.drawing.currentColors.value.tRNA, '#565656');
  assert.equal(live.state.appliedPaletteColors.value.tRNA, '#565656');

  const queued = setup({ instant: false });
  queued.manager.selectPalette('forest');
  await queued.choose('keep');
  const applied = JSON.stringify(queued.state.appliedPaletteColors.value);
  queued.manager.editDefaultColor('tRNA', '#676767');
  assert.equal(queued.drawing.pendingPaletteName.value, 'forest');
  assert.equal(queued.drawing.pendingPaletteColors.value.tRNA, '#676767');
  assert.equal(JSON.stringify(queued.state.appliedPaletteColors.value), applied);
});
