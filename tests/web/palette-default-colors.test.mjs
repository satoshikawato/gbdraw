// D-15 (Owner 2026-10-09): the palette is the base layer and the user default
// colors (values that differ from the selected palette's color, as the `-d`
// table over `-p`) win over it. A palette switch keeps them, after a dialog when
// any exist; Default colors Reset asks before it discards them. A dialog
// choice is one History step, and Cancel records none (PD-OI-088).
import assert from 'node:assert/strict';
import { test } from 'node:test';
import { createResultsManager } from '../../gbdraw/web/js/app/results.js';
import { createHistoryManager } from '../../gbdraw/web/js/services/history.js';
import { createDialogChoice } from '../../gbdraw/web/js/app/history-inputs.js';
import { normalizePaletteColors } from '../../gbdraw/web/js/utils/color-utils.js';

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

const setup = ({ colors = editedColors(), instant = true } = {}) => {
  const drawing = {
    selectedPalette: ref('default'), currentColors: ref(colors), pendingPaletteName: ref(''), pendingPaletteColors: ref({})
  };
  const paletteColorsDialog = { show: false, kind: 'switch', fromPalette: '', toPalette: '', count: 0, keysText: '' };
  const state = {
    paletteDefinitions: ref(PALETTES), paletteInstantPreviewEnabled: ref(instant), appliedPaletteName: ref('default'),
    appliedPaletteColors: ref({ ...colors }), activeDrawing: () => drawing, paletteColorsDialog
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
  const manager = createResultsManager({ state, closeAfterDialogChoice: dialogChoice.closeAfterChoice });
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
    show: true, kind: 'switch', fromPalette: 'default', toPalette: 'forest', count: 3,
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
    show: true, kind: 'switch', fromPalette: 'forest', toPalette: 'default', count: 1, keysText: 'tRNA'
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
    show: true, kind: 'reset', fromPalette: 'default', toPalette: 'default', count: 3,
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
// draws it, not the palette's `default` color.
test('a restyle paints an Auto key with the applied palette color', async () => {
  const { createSvgStyles } = await import('../../gbdraw/web/js/app/svg-styles.js');
  const { withDrawings } = await import('./helpers/drawing-state.mjs');
  const feature = { id: 'f1', svg_id: 'f1', type: 'CDS' };
  const attributes = { id: 'f1', 'data-gbdraw-feature-id': 'f1', 'data-gbdraw-feature-part': 'block', fill: '#000000' };
  const element = { getAttribute: (key) => attributes[key] ?? null, setAttribute: (key, value) => { attributes[key] = value; } };
  const svg = {
    querySelectorAll: (selector) => (selector.includes('data-gbdraw-feature-id') ? [element] : []),
    getElementById: () => null
  };
  const state = withDrawings({
    mode: ref('circular'), svgContent: ref('<svg/>'), svgContainer: ref({ querySelector: () => svg }),
    extractedFeatures: ref([feature]), featuresBySvgId: ref(new Map([['f1', feature]])), manualSpecificRules: [],
    featureColorOverrides: {}, legendColorOverrides: {}, pairwiseMatchFactors: ref({}),
    paletteDefinitions: ref(PALETTES), appliedPaletteName: ref('default'),
    appliedPaletteColors: ref({ ...paletteOf('default'), CDS: null })
  });
  const styles = createSvgStyles({
    state, watch() {}, nextTick: (fn) => fn?.(), commitActiveResultEdit: () => true, projectPaletteAndRules: () => true
  });
  styles.applyPaletteToSvg();
  assert.equal(attributes.fill, '#AABBCC');
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

// OV-276: a palette repaint gives the composition root the color of each
// Legend row the palette colors (rows with a Legend color keep theirs), and
// the Legend owner shows those colors in the Legend panel rows.
test('a palette repaint reports its Legend row colors and the Legend panel rows take them', async () => {
  const { createSvgStyles } = await import('../../gbdraw/web/js/app/svg-styles.js');
  const { createLegendEntryActions } = await import('../../gbdraw/web/js/app/legend/entry-actions.js');
  const { withDrawings } = await import('./helpers/drawing-state.mjs');
  const element = (attributes) => ({
    getAttribute: (key) => attributes[key] ?? null, setAttribute: (key, value) => { attributes[key] = value; }
  });
  const row = (caption, fill) => ({
    ...element({ 'data-legend-key': caption }),
    querySelectorAll: (selector) => (selector === 'path' ? [element({ fill })] : [])
  });
  const rows = [row('CDS', '#AABBCC'), row('other tRNAs', '#111111'), row('Manual', '#010101'), row('rRNA', '#121212')];
  const featureLegend = { querySelectorAll: (selector) => (selector === 'g[data-legend-key]' ? rows : []) };
  const legend = { querySelector: (selector) => (selector === '#feature_legend' ? featureLegend : null) };
  const svg = { querySelectorAll: () => [], getElementById: (id) => (id === 'legend' ? legend : null) };
  const feature = { id: 'f1', svg_id: 'f1', type: 'CDS' };
  const entries = [
    { caption: 'CDS', color: '#AABBCC' }, { caption: 'other tRNAs', color: '#111111' },
    { caption: 'Manual', color: '#010101' }, { caption: 'rRNA', color: '#121212' }
  ];
  const state = withDrawings({
    mode: ref('circular'), svgContent: ref('<svg/>'), svgContainer: ref({ querySelector: () => svg }),
    extractedFeatures: ref([feature]), featuresBySvgId: ref(new Map()), manualSpecificRules: [],
    featureColorOverrides: {}, legendColorOverrides: { rRNA: '#121212' }, pairwiseMatchFactors: ref({}),
    paletteDefinitions: ref(PALETTES), appliedPaletteName: ref('default'),
    appliedPaletteColors: ref({ ...paletteOf('default'), CDS: '#123456', tRNA: '#333333', rRNA: '#444444' }),
    results: ref([]), legendEntries: ref(entries), originalLegendOrder: ref([]), originalLegendColors: ref({})
  });
  const styles = createSvgStyles({
    state, watch() {}, nextTick: (fn) => fn?.(), commitActiveResultEdit: () => true, projectPaletteAndRules: () => true
  });
  const legendRowColors = styles.applyPaletteToSvg();
  assert.deepEqual([...legendRowColors], [['CDS', '#123456'], ['other tRNAs', '#333333']]);
  const actions = createLegendEntryActions({ state });
  assert.equal(actions.setPaletteLegendEntryColors(legendRowColors), true);
  assert.deepEqual(state.activeDrawing().legendEntries.value, [
    { caption: 'CDS', color: '#123456' }, { caption: 'other tRNAs', color: '#333333' }, entries[2], entries[3]
  ]);
  assert.equal(state.activeDrawing().legendEntries.value[2], entries[2], 'an unchanged row keeps its object');
  assert.equal(actions.setPaletteLegendEntryColors(legendRowColors), false);
  assert.equal(actions.setPaletteLegendEntryColors(undefined), false);
});
