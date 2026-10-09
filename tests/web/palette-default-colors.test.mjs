// D-15 (Owner 2026-10-09): the palette is the base layer and the user default
// colors (values that differ from the selected palette's color, as the `-d`
// table over `-p`) win over it. A palette switch keeps them, after a dialog when
// any exist; Default colors Reset asks before it discards them. A dialog
// choice is one History step, and Cancel records none (PD-OI-088).
import assert from 'node:assert/strict';
import { test } from 'node:test';
import { createResultsManager } from '../../gbdraw/web/js/app/results.js';
import { createHistoryManager } from '../../gbdraw/web/js/services/history.js';
import { dialogChoiceWithHistory } from '../../gbdraw/web/js/app/history-inputs.js';
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
  const manager = createResultsManager({ state });
  const history = createHistoryManager({
    buildIntent: () => ({
      palette: drawing.selectedPalette.value, colors: { ...drawing.currentColors.value },
      applied: { ...state.appliedPaletteColors.value }, pending: drawing.pendingPaletteName.value
    }),
    signatureFor: (value) => JSON.stringify(value),
    applyIntent: () => {}, buildCheckpoint: () => ({}), applyCheckpoint: () => {}
  });
  // app-setup.js wires the dialog's choices the same way.
  const choose = dialogChoiceWithHistory(
    history, () => 'Change setting', manager.handlePaletteColorsChoice, manager.cancelPaletteColorsDialog
  );
  return { drawing, state, manager, history, choose, paletteColorsDialog };
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

// Q5 (Owner, pending): switching back to the applied palette while another is
// queued restores the applied colors without a dialog, as before D-15.
test('switching back to the applied palette while one is queued restores the applied colors', async () => {
  const { drawing, state, manager, choose, paletteColorsDialog } = setup({ instant: false });
  const applied = { ...state.appliedPaletteColors.value };
  manager.selectPalette('forest');
  await choose('palette');
  drawing.currentColors.value = { ...drawing.currentColors.value, tRNA: '#555555' };
  manager.selectPalette('default');
  assert.equal(paletteColorsDialog.show, false);
  assert.equal(drawing.selectedPalette.value, 'default');
  assert.deepEqual(drawing.currentColors.value, applied);
  assert.equal(drawing.pendingPaletteName.value, '');
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

test('the dialog names three keys and counts the rest (Q6 A)', () => {
  const { manager, paletteColorsDialog } = setup({ colors: { ...editedColors(), CDS: '#000000', rRNA: '#010101' } });
  manager.selectPalette('forest');
  assert.equal(paletteColorsDialog.count, 5);
  assert.equal(paletteColorsDialog.keysText, 'CDS, tRNA, rRNA, and 2 more');
});

// OV-262: Default colors Auto shows the palette's color live, as Generate
// draws it without a `-d` row, not the palette's `default` color.
test('an Auto default color shows the palette color on the Result', () => {
  const { drawing, state, manager } = setup({ colors: paletteOf('default') });
  drawing.currentColors.value = { ...drawing.currentColors.value, CDS: null };
  manager.syncPaletteDraftState();
  assert.equal(state.appliedPaletteColors.value.CDS, '#AABBCC');
});

test('setDefaultColor writes the default color and applies it to the Result', () => {
  const { drawing, state, manager } = setup();
  manager.setDefaultColor(drawing, 'CDS', '#123456');
  assert.equal(drawing.currentColors.value.CDS, '#123456');
  assert.equal(state.appliedPaletteColors.value.CDS, '#123456');
  assert.equal(manager.readUserDefaultColor(drawing, 'CDS'), '#123456');
});
