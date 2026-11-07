// Issue #1427 (PD-OI-114 scenario revision 2, Q-D1 A): on a figure drawn to a requested size, a Legend rename,
// deletion, Restore, Move, or Sort, or a Canvas padding change, shows at once and then asks the automatic rerender to
// fit the figure again; its Undo and Redo do too. A color or stroke edit, and an Auto figure, rerender nothing. While
// the rerender runs, the line under the Figure size field says the figure is being fitted.
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import test from 'node:test';

import { installFakeSvgDom } from './fake-svg-dom.mjs';
import { savedDrawingView } from './helpers/saved-drawing-view.mjs';

installFakeSvgDom();
globalThis.CSS ||= { escape: (value) => String(value) };

const { stepRefitsDrawnLegend } = await import('../../gbdraw/web/js/services/legend-draft-rows.js');
const { SETTING_APPLICATION, refitFieldsBeforeStep, settingApplicationHelp } = await import('../../gbdraw/web/js/services/setting-application.js');
const { requestedCanvasOf } = await import('../../gbdraw/web/js/app/legend-layout/composition-actions.js');
const { createFigureSizeField } = await import('../../gbdraw/web/js/app/figure-size-field.js');
const { createDefaultForm } = await import('../../gbdraw/web/js/services/session-active-config-contract.js');

const row = (key, caption = key, color = '#cccccc') => ({ caption, originalCaption: key, color, xPos: 0, yPos: 0, featureIds: [] });
const NO_PADDING = { top: 0, right: 0, bottom: 0, left: 0 };
const INVENTORY = ['CDS', 'tRNA', 'rRNA'];
const ENTRIES = [row('CDS'), row('tRNA'), row('rRNA')];
// The drawing's values as a step left them, with the shown Result's inventory and Python's rows of every Result.
const drawing = (values = {}) => ({
  legendEntries: ENTRIES, deletedLegendEntries: [], canvasPadding: NO_PADDING, originalLegendOrder: INVENTORY, ...values,
  // Python's rows hold the shown Result's inventory.
  pythonRowKeys: values.pythonRowKeys ?? values.originalLegendOrder ?? INVENTORY
});
const entriesChange = (before, after) => ({ path: ['editorState', 'legend', 'entries'], before, after, beforeHas: true, afterHas: true });

test('the Legend editor and Canvas padding rows are live, name the draft fields Python sizes, and say a sized figure is fitted again', () => {
  const rows = Object.fromEntries(SETTING_APPLICATION.map((entry) => [entry.id, entry]));
  assert.equal(rows['legend-editor'].application, 'live');
  assert.deepEqual(rows['legend-editor'].refitPaths, {
    'editorState.legend.entries': 'legendEntries', 'editorState.legend.deletedEntries': 'deletedLegendEntries', 'config.rules': 'rules'
  });
  assert.equal(rows['canvas-padding'].application, 'live');
  assert.deepEqual(rows['canvas-padding'].refitPaths, { 'ui.canvasPadding': 'canvasPadding' });
  assert.equal(rows['canvas-padding'].requestPaths, undefined, 'editor intent, not a draft setting');
  assert.equal(settingApplicationHelp('legend-editor'),
    'Live edit: legend text, order, stroke, and removal update the current Result. With a Figure size set, a rename, removal, or reorder then redraws the figure to fit it.');
  assert.equal(settingApplicationHelp('canvas-padding'),
    'Live edit: the padding updates the current Result. With a Figure size set, the figure is then redrawn to fit it.');
  assert.equal(settingApplicationHelp('features-editor').includes('Figure size'), false);
});

test('a rename, a deletion, a Restore, a Move, and a Sort change what Python draws; a color edit does not', () => {
  const renamed = [row('CDS', 'Genes'), row('tRNA'), row('rRNA')];
  assert.equal(stepRefitsDrawnLegend([entriesChange(ENTRIES, renamed)], undefined, drawing({ legendEntries: renamed })), true, 'rename');
  const moved = [row('tRNA'), row('CDS'), row('rRNA')];
  assert.equal(stepRefitsDrawnLegend([entriesChange(ENTRIES, moved)], undefined, drawing({ legendEntries: moved })), true, 'Move or Sort');
  const deletedRows = [row('tRNA')];
  const afterDelete = [row('CDS'), row('rRNA')];
  const deletion = [entriesChange(ENTRIES, afterDelete),
    { path: ['editorState', 'legend', 'deletedEntries'], before: [], after: deletedRows, beforeHas: true, afterHas: true }];
  assert.equal(stepRefitsDrawnLegend(deletion, undefined, drawing({ legendEntries: afterDelete, deletedLegendEntries: deletedRows })), true, 'deletion');
  // Restore is the deletion's reverse; its Undo puts the drawing back where the deletion left it.
  const restore = deletion.map((change) => ({ ...change, before: change.after, after: change.before }));
  assert.equal(stepRefitsDrawnLegend(restore, undefined, drawing()), true, 'Restore');
  // A Circular batch: Result B draws tRNA, the shown Result A does not. The Undo of tRNA's deletion made on B, while A
  // is shown, changes what the request hides, so the sized figure is fitted again: the check names Python's rows of
  // every Result, as the request does (`pythonLegendRowKeys`), not only the shown Result's.
  const shownA = { legendEntries: [row('CDS'), row('rRNA')], originalLegendOrder: ['CDS', 'rRNA'] };
  const deletedOnB = [{ path: ['editorState', 'legend', 'deletedEntries'], before: [], after: deletedRows, beforeHas: true, afterHas: true }];
  assert.equal(stepRefitsDrawnLegend(deletedOnB, 'undo', drawing({ ...shownA, pythonRowKeys: [...INVENTORY] })), true,
    'Undo, on Result A, of a deletion made on Result B');
  assert.equal(stepRefitsDrawnLegend(deletedOnB, 'undo', drawing({ ...shownA, pythonRowKeys: ['CDS', 'rRNA'] })), false,
    'a row no Result draws is not Python\'s');
  const recolored = [row('CDS', 'CDS', '#ff0000'), row('tRNA'), row('rRNA')];
  assert.equal(stepRefitsDrawnLegend([entriesChange(ENTRIES, recolored),
    { path: ['editorState', 'legend', 'colorOverrides', 'CDS'], before: undefined, after: '#ff0000', beforeHas: false, afterHas: true }],
  undefined, drawing({ legendEntries: recolored })), false, 'a Legend row color stays the Web\'s');
  assert.equal(stepRefitsDrawnLegend([{ path: ['editorState', 'legend', 'strokeOverrides', 'CDS'], before: undefined, after: { color: '#000000' }, beforeHas: false, afterHas: true }],
    undefined, drawing()), false, 'a Legend row stroke stays the Web\'s');
  assert.equal(stepRefitsDrawnLegend(null, undefined, drawing()), false, 'a checkpoint step names no change');
  assert.equal(stepRefitsDrawnLegend([{ path: ['config', 'form', 'plot_title'], before: '', after: 'T', beforeHas: true, afterHas: true }],
    undefined, drawing()), false, 'another setting');
});

test('a row only the Web draws is not Python\'s to fit; a rule row\'s rename is', () => {
  const added = [...ENTRIES, row('My row')];
  assert.equal(stepRefitsDrawnLegend([entriesChange(ENTRIES, added)], undefined, drawing({ legendEntries: added })), false);
  // A rule row is in the inventory (Python draws it from the color table): its rename reaches Python through its rule.
  const ruled = [...ENTRIES, row('Kinase')];
  const renamedRule = [...ENTRIES, row('Kinase', 'Kinases')];
  assert.equal(stepRefitsDrawnLegend([entriesChange(ruled, renamedRule)], undefined,
    drawing({ legendEntries: renamedRule, originalLegendOrder: [...INVENTORY, 'Kinase'] })), true);
});

// OV-2556: a Legend rename of a feature type row whose features carry other values writes one color rule per value
// with the new name; the row is then the rule row of that name, and the shown Result's inventory names it so. The
// rename is a checkpoint step, whose changes compare with the intent History last captured; Generate writes the Legend
// rows outside History, so that side may hold no rows. The request's leaves are then empty on both sides, but Python
// draws the row from the rules' caption, so the step fits the figure again.
test('a rename that writes color rules fits again; a rule color edit does not', () => {
  const rule = (val, cap, color = '#54bcf8') => ({ feat: 'CDS', qual: 'product', val, color, cap });
  const ruled = [row('Coding sequences'), row('tRNA'), row('rRNA')];
  const rules = [rule('^alpha$', 'Coding sequences'), rule('^beta$', 'Coding sequences')];
  const rename = [entriesChange([], ruled), { path: ['config', 'rules'], before: [], after: rules, beforeHas: true, afterHas: true }];
  const after = drawing({ legendEntries: ruled, originalLegendOrder: ['Coding sequences', 'tRNA', 'rRNA'], rules });
  assert.equal(stepRefitsDrawnLegend(rename, undefined, after), true, 'the rename');
  assert.equal(stepRefitsDrawnLegend(rename, 'undo', drawing({ rules: [] })), true, 'its Undo');
  const recolored = rules.map((entry) => ({ ...entry, color: '#ff0000' }));
  assert.equal(stepRefitsDrawnLegend([{ path: ['config', 'rules'], before: rules, after: recolored, beforeHas: true, afterHas: true }],
    undefined, { ...after, rules: recolored }), false, 'a rule color stays the size');
});

test('an Undo or a Redo of such a step fits again; the Undo reads the drawing as it was before the step', () => {
  const renamed = [row('CDS', 'Genes'), row('tRNA'), row('rRNA')];
  const change = entriesChange(ENTRIES, renamed);
  assert.equal(stepRefitsDrawnLegend([change], 'undo', drawing()), true);
  assert.equal(stepRefitsDrawnLegend([change], 'redo', drawing({ legendEntries: renamed })), true);
  // The other side of an Undo is the step's `after`.
  assert.deepEqual(refitFieldsBeforeStep([change], drawing(), true).legendEntries, renamed);
  assert.deepEqual(refitFieldsBeforeStep([change], drawing({ legendEntries: renamed }), false).legendEntries, ENTRIES);
});

test('a Canvas padding change fits again when it changes the padding Python draws', () => {
  const side = (before, after) => [{ path: ['ui', 'canvasPadding', 'top'], before, after, beforeHas: true, afterHas: true }];
  assert.equal(stepRefitsDrawnLegend(side(0, 10), undefined, drawing({ canvasPadding: { ...NO_PADDING, top: 10 } })), true);
  // A typed negative side is stored as 0 (D-59), so it makes no padding change.
  // A change of the whole padding object (Reset) holds every side.
  assert.equal(stepRefitsDrawnLegend([{ path: ['ui', 'canvasPadding'], before: { ...NO_PADDING, left: 20 }, after: NO_PADDING, beforeHas: true, afterHas: true }],
    undefined, drawing()), true);
});

// A Gallery Result as the app mounts it (the fake DOM keeps attribute text as written).
const COMPOSITION = 'data-gbdraw-composition';
const decodeAttribute = (value) => value.replace(/&quot;/g, '"').replace(/&lt;/g, '<').replace(/&gt;/g, '>')
  .replace(/&amp;/g, '&');
const galleryRoot = (base) => {
  const saved = JSON.parse(readFileSync(new URL(`../../gbdraw/web/gallery/sessions/${base}.gbdraw-session.json`, import.meta.url), 'utf8'));
  const content = savedDrawingView(saved, saved.ui.activeDrawingId).results[0].content;
  const svg = new DOMParser().parseFromString(content, 'image/svg+xml').documentElement;
  svg.setAttribute(COMPOSITION, decodeAttribute(svg.getAttribute(COMPOSITION)));
  return svg;
};
const sized = (svg, padding = null) => {
  const metadata = JSON.parse(svg.getAttribute(COMPOSITION));
  const [, , width, height] = svg.getAttribute('viewBox').split(/\s+/).map(Number);
  metadata.requestedCanvas = { widthPx: width, heightPx: null };
  if (padding) metadata.spacing.canvasPadding = padding;
  svg.setAttribute(COMPOSITION, JSON.stringify(metadata));
  return svg;
};

test('the size line counts only the padding the Web adds beyond the padding Python drew', () => {
  const padding = { top: 10, right: 0, bottom: 0, left: 0 };
  const drawnPadded = requestedCanvasOf(sized(galleryRoot('HmmtDNA_basic_circular'), padding));
  assert.equal(drawnPadded.addsPadding(padding), false, 'Python drew it: no cause once the refit is done');
  assert.equal(drawnPadded.addsPadding({ ...padding, top: 25 }), true, 'a change the refit has not drawn yet');
  assert.equal(drawnPadded.addsPadding(NO_PADDING), false, 'less padding than drawn adds no room');
  assert.equal(requestedCanvasOf(sized(galleryRoot('HmmtDNA_basic_circular'))).addsPadding(padding), true, 'a Result that records none');
});

test('while the rerender fits a sized figure, the line under the field says so instead of why it is larger', () => {
  const form = createDefaultForm();
  const shownDrawing = { id: 'circular', mode: 'circular', form };
  let shown = { first: '190 × 143.2 mm', second: '718 × 541 px', fitting: true,
    larger: { sides: [{ side: 'width', length: '180 mm' }], causes: ['legend-edits'] } };
  const field = createFigureSizeField({ state: { activeDrawing: () => shownDrawing }, watch: () => {}, shownSize: () => shown });
  assert.equal(field.figureSizeLine.value, 'Figure: 190 × 143.2 mm (718 × 541 px)');
  assert.equal(field.figureSizeFitting.value, 'Fitting the figure to the requested size.');
  assert.equal(field.figureSizeLarger.value, '');
  shown = { first: '180 × 143.2 mm', second: '680 × 541 px', fitting: false, larger: null };
  assert.equal(field.figureSizeFitting.value, '');
  assert.equal(field.figureSizeLarger.value, '');
});

// D-24 with the sized refit: a checkpoint step whose listener asked for the refit after its commit marks its figure
// undrawn, so History keeps the rerender's figure as the step's (Redo then shows the fitted figure without another
// Python run), and a restore of a figure that stayed undrawn (a refused refit) tells the listener so.
test('a refitted checkpoint step keeps the fitted figure; a restore of an unfitted one says it is undrawn', async () => {
  const { createHistoryManager } = await import('../../gbdraw/web/js/services/history.js');
  let rows = 'CDS';
  let figure = 'fitted-0';
  /** @type {unknown[][]} */
  const heard = [];
  const history = createHistoryManager({
    buildIntent: () => ({ rows }),
    applyIntent: (intent) => { rows = intent.rows; },
    buildCheckpoint: () => ({ rows, figure }),
    applyCheckpoint: (checkpoint) => { rows = checkpoint.rows; figure = checkpoint.figure; },
    figureDrawn: () => true
  });
  history.setSettingsStepListener((changes, direction, undrawn) => {
    heard.push([changes === null ? null : 'changes', direction, undrawn]);
    if (changes && !direction) history.markLastStepUndrawn();
  });
  await history.captureBaseline('initial');
  // The rename's step was drawn as far as its caller knew; the refit it asked for makes it undrawn until shown.
  await history.runUndoableCheckpoint('Rename legend item', () => { rows = 'Genes'; figure = 'live-1'; },
    { drawn: { change: true } });
  figure = 'fitted-1';
  assert.equal(await history.keepRedrawnResult(), true, 'the fitted figure is the step\'s');
  await history.undo();
  assert.deepEqual({ rows, figure }, { rows: 'CDS', figure: 'fitted-0' });
  await history.redo();
  assert.deepEqual({ rows, figure }, { rows: 'Genes', figure: 'fitted-1' }, 'Redo shows the fitted figure');
  assert.deepEqual(heard.slice(-2), [[null, undefined, false], [null, undefined, false]]);
  // A refit that never came (refused): the step keeps the live figure, and its Redo says it is undrawn.
  await history.runUndoableCheckpoint('Rename legend item', () => { rows = 'Long name'; figure = 'live-2'; },
    { drawn: { change: true } });
  await history.undo();
  assert.deepEqual(heard.at(-1), [null, undefined, false], 'the figure before it is drawn');
  await history.redo();
  assert.deepEqual({ rows, figure }, { rows: 'Long name', figure: 'live-2' });
  assert.deepEqual(heard.at(-1), [null, undefined, true], 'the listener refits it again');
  // A commit whose caller asks for the automatic rerender itself says so: the listener asks for no second one.
  await history.runUndoableCheckpoint('Add color rule', () => { rows = 'Alpha'; figure = 'live-3'; },
    { drawn: { change: false } });
  assert.deepEqual(heard.at(-1), ['changes', undefined, true]);
  // An intent step has no figure to mark.
  await history.runUndoable('Edit', () => { rows = 'Other'; });
  history.markLastStepUndrawn();
  assert.equal(history.stackEntries('undo').at(-1)?.type, 'intent');
});

// OV-2560 (#1634 residual): each side of every checkpoint step reads the drawing's one reader of whether its figure is
// the one its state draws, as History captures that side. A step whose caller says nothing of its figure (a file
// import, the Merge of a Legend rename onto a deleted row) that starts while the refit of the step before it is
// pending captures a figure Python has not fitted: its Undo shows that figure at once and tells the listener it is
// undrawn, so the automatic rerender fits it. A refit shown while the step is the last one is kept as its figure, so
// its Redo runs no Python (D-48). The reader is scoped to the History's own drawing.
test('a checkpoint step captured while its drawing owes a refit restores undrawn, whoever runs it', async () => {
  const { createHistoryRouter } = await import('../../gbdraw/web/js/services/history.js');
  const values = { circular: { colors: 'a', figure: 'fitted-0' }, other: { colors: 'x', figure: 'fitted-x' } };
  let shown = 'circular';
  /** @type {Set<string>} */
  const owed = new Set();
  /** @type {unknown[][]} */
  const heard = [];
  const history = createHistoryRouter({
    drawingIds: ['circular', 'other'],
    activeDrawingId: () => shown,
    sessionHistory: { buildIntent: () => ({}), applyIntent: () => {} },
    managerOptions: (id) => ({
      buildIntent: () => ({ colors: values[id].colors }),
      applyIntent: (intent) => { values[id].colors = intent.colors; },
      buildCheckpoint: () => ({ ...values[id] }),
      applyCheckpoint: (checkpoint) => { Object.assign(values[id], checkpoint); }
    })
  });
  history.setSettingsStepListener((changes, direction, undrawn) => heard.push([changes === null ? null : 'changes', direction, undrawn]));
  history.setFigureDrawnReader((drawingId) => !owed.has(drawingId));
  await history.captureBaseline('initial');
  // A Canvas padding refit is pending: the figure shows the padding as the Web drew it.
  owed.add('circular');
  values.circular.figure = 'padded-live';
  await history.runUndoableCheckpoint('Change uploaded file', () => { values.circular.colors = 'b'; });
  // The refit is shown while the import is the last step: History keeps it as the step's figure.
  owed.delete('circular');
  values.circular.figure = 'fitted-1';
  assert.equal(await history.keepRedrawnResult(), true, 'the fitted figure is the step\'s');
  await history.undo();
  assert.deepEqual(values.circular, { colors: 'a', figure: 'padded-live' }, 'Undo shows the figure before the step at once');
  assert.deepEqual(heard.at(-1), [null, undefined, true], 'and asks the listener to fit it');
  await history.redo();
  assert.deepEqual(values.circular, { colors: 'b', figure: 'fitted-1' });
  assert.deepEqual(heard.at(-1), [null, undefined, false], 'Redo shows the kept fitted figure with no Python run');
  // A refit owed by another drawing does not mark this drawing's figure: the step ends in the History that opened
  // it, which reads its own drawing, also when the other drawing is shown by then.
  owed.add('other');
  await history.runUndoableCheckpoint('Change uploaded file', () => { values.circular.colors = 'c'; shown = 'other'; });
  shown = 'circular';
  await history.undo();
  assert.deepEqual(heard.at(-1), [null, undefined, false], 'the figure before the step is drawn');
  await history.redo();
  assert.deepEqual(heard.at(-1), [null, undefined, false], 'the figure as the step ended is not marked undrawn');
  await history.undo();
  // A refit that never came (refused) leaves the step's figure undrawn: its Redo asks again (D-49).
  owed.delete('other');
  owed.add('circular');
  values.circular.figure = 'refused-live';
  await history.runUndoableCheckpoint('Change uploaded file', () => { values.circular.colors = 'd'; });
  await history.undo();
  assert.deepEqual(heard.at(-1), [null, undefined, true]);
  await history.redo();
  assert.deepEqual(heard.at(-1), [null, undefined, true], 'the Redo asks for the refit again');
});
