const { test, expect } = require('@playwright/test');
const { evaluateWithRetainedPromise, openApp, reveal } = require('./helpers/app-lifecycle.cjs');
const { readFileSync, writeFileSync } = require('node:fs');
const { gunzipSync } = require('node:zlib');
const { load, generate, switchMode, download } = require('./helpers/mode-transition.cjs');
const { expectLiveEqualsGenerate, settleLive } = require('./helpers/live-generate-parity.cjs');
const { RING_SESSION_FILE, ringSessionUpload } = require('./helpers/ring-session.cjs');

// R-3, Q3 (Owner, 2026-10-04): every edit that changes the direction of a
// draft feature slot runs through one transition, which asks before it leaves
// a lane Feature placement undrawable. Reset applies the edit and removes those
// rows as one History step; Cancel keeps the control's value and records none.
const HMMT = 'gbdraw/web/gallery/sessions/HmmtDNA_basic_circular.gbdraw-session.json';
const LAMBDA = 'gbdraw/web/gallery/sessions/lambda_basic_linear.gbdraw-session.json';

const placeLane = (page, side) => page.evaluate(async (lane) => {
  const app = window.__GBDRAW_APP__;
  const feature = app.extractedFeatures.find((item) => item.type === 'CDS');
  await app.featurePlacementActions.setPlacement([feature], lane);
  return JSON.stringify([feature.record_key, feature.biological_feature_id]);
}, side);
const layoutState = (page, key) => page.evaluate(async (overrideKey) => {
  const { state } = await import('./js/state.js');
  const lane = (mode) => {
    const slot = state.activeDrawing().adv[`${mode}_track_slots`].find((entry) => entry.renderer === 'features');
    return slot ? { enabled: slot.enabled !== false, side: slot.side, lane: slot.params?.lane_direction ?? null } : null;
  };
  return { trackType: state.activeDrawing().form.track_type, separate: state.activeDrawing().form.separate_strands,
    circularStack: state.activeDrawing().adv.circular_track_slots_enabled, circular: lane('circular'),
    linearStack: state.activeDrawing().adv.linear_track_slots_enabled, linear: lane('linear'),
    row: state.activeDrawing().featurePlacementOverrides[overrideKey]?.placement?.side || null,
    undo: window.__GBDRAW_HISTORY__.getUndoCount() };
}, key);
const dialog = (page) => page.getByRole('dialog', { name: 'Reset Feature placements?', exact: true });
const resetButton = (page) => dialog(page).getByRole('button', { name: 'Reset 1 placement to Auto', exact: true });
const cancelButton = (page) => dialog(page).getByRole('button', { name: 'Keep current layout', exact: true });
const openStack = async (page, mode) => {
  const toggle = await reveal(page.locator(`[aria-controls="${mode}-custom-track-slots-panel"]`));
  if (await page.locator(`#${mode}-custom-track-slots-panel`).count() === 0) await toggle.click();
  await expect(page.locator(`#${mode}-custom-track-slots-panel`)).toBeVisible();
};
const stackRow = (page, mode) => page.locator('.track-slot-row-head', {
  has: page.locator(`[aria-label="Enable ${mode} track slot features"]`)
});
const undo = (page) => page.getByRole('button', { name: /^Undo/ }).first().click();
// The Axis moves, Duplicate, and Remove/Delete of a stack row are items of its
// Row actions menu (Owner decision 2026-10-18); Move up and Move down stay buttons.
const rowMenu = (row) => row.locator('.track-slot-row-menu');
const menuToggle = (row) => rowMenu(row).locator(':scope > summary');
// The items of a closed menu are not in the accessibility tree; read them by class.
const menuItems = (row) => rowMenu(row).locator('.action-menu-item');
const menuItem = (row, label) => menuItems(row).filter({ hasText: new RegExp(`^\\s*${label}\\s*$`) });
const runRowAction = async (row, label) => {
  await menuToggle(row).click();
  await rowMenu(row).getByRole('button', { name: label, exact: true }).click();
};
const slotRow = (page, mode, id) => page.locator(`[data-capture="${mode}-track-slot-${id}"]`);
const stackOrder = (page, mode) => page.locator(`#${mode}-custom-track-slots-panel [data-capture^="${mode}-track-"]`)
  .evaluateAll((rows) => rows.map((row) => row.dataset.capture.replace(/^(?:circular|linear)-track-(?:slot-)?/, '')));

test('Circular custom stack Reset and row moves ask before dropping a lane placement (R-3)', async ({ browser }, testInfo) => {
  test.setTimeout(300000);
  const page = await load(browser, HMMT);
  try {
    await openStack(page, 'circular');
    await page.getByRole('checkbox', { name: 'Use custom stack', exact: true }).check();
    const key = await placeLane(page, 'outward');
    const placed = await layoutState(page, key);
    expect(placed).toMatchObject({ trackType: 'middle', circularStack: true,
      circular: { enabled: true, lane: 'split' }, row: 'outward' });

    // Preset Reset: the dialog asks before the stack changes; Cancel keeps it.
    const tuckin = page.getByRole('button', { name: 'Reset to Tuckin', exact: true });
    await tuckin.click();
    await expect(dialog(page)).toContainText('Choosing "Reset to Tuckin" removes the lane for 1 feature you placed manually.');
    await expect(resetButton(page)).toBeFocused();
    expect(await layoutState(page, key)).toEqual(placed);
    await page.screenshot({ path: testInfo.outputPath('preset-reset-dialog.png') });
    await cancelButton(page).click();
    await expect(dialog(page)).toBeHidden();
    await expect(tuckin).toBeFocused();
    expect(await layoutState(page, key)).toEqual(placed);

    // Row reorder: Reset applies the move and removes the row as one step.
    await runRowAction(stackRow(page, 'circular'), 'Move outside Axis');
    await expect(dialog(page)).toContainText('Choosing "Move outside Axis" removes the lane for 1 feature you placed manually.');
    expect(await layoutState(page, key)).toEqual(placed);
    // Keep current layout returns focus to the row's Row actions.
    await cancelButton(page).click();
    await expect(menuToggle(stackRow(page, 'circular'))).toBeFocused();
    expect(await layoutState(page, key)).toEqual(placed);
    await runRowAction(stackRow(page, 'circular'), 'Move outside Axis');
    await resetButton(page).click();
    await expect(dialog(page)).toBeHidden();
    expect(await layoutState(page, key)).toEqual({ ...placed,
      circular: { enabled: true, side: 'outside', lane: 'outside' }, row: null, undo: placed.undo + 1 });
    await undo(page);
    await expect.poll(() => layoutState(page, key)).toEqual(placed);

    await tuckin.click();
    await resetButton(page).click();
    await expect(dialog(page)).toBeHidden();
    expect(await layoutState(page, key)).toEqual({ ...placed, trackType: 'tuckin',
      circular: { enabled: true, side: 'inside', lane: 'inside' }, row: null, undo: placed.undo + 1 });
    await generate(page);
  } finally {
    await page.context().close();
  }
});

test('Linear Use custom stack asks before dropping a Linear lane, and Circular Separate Strands leaves it (R-3)', async ({ browser }, testInfo) => {
  test.setTimeout(300000);
  const page = await load(browser, LAMBDA);
  try {
    // A saved custom stack whose feature row sits above the Axis has no lanes.
    await openStack(page, 'linear');
    const custom = page.getByRole('checkbox', { name: 'Use custom stack', exact: true });
    await custom.check();
    await runRowAction(stackRow(page, 'linear'), 'Move above Axis');
    await custom.uncheck();
    const strands = await reveal(page.getByRole('checkbox', { name: 'Separate Strands', exact: true, includeHidden: true }));
    await strands.uncheck();
    const key = await placeLane(page, 'above');
    const placed = await layoutState(page, key);
    expect(placed).toMatchObject({ separate: false, linearStack: false, linear: { side: 'above' }, row: 'above' });

    // Use custom stack: Cancel keeps the checkbox off; Reset is one step.
    await custom.click();
    await expect(dialog(page)).toContainText('Changing "Use custom stack" to "On" removes the lane for 1 feature you placed manually.');
    await expect(custom).not.toBeChecked();
    expect(await layoutState(page, key)).toEqual(placed);
    await page.keyboard.press('Escape');
    await expect(dialog(page)).toBeHidden();
    await expect(custom).not.toBeChecked();
    await expect(custom).toBeFocused();
    expect(await layoutState(page, key)).toEqual(placed);
    await custom.click();
    await resetButton(page).click();
    await expect(custom).toBeChecked();
    expect(await layoutState(page, key)).toEqual({ ...placed, linearStack: true, row: null, undo: placed.undo + 1 });
    await undo(page);
    await expect.poll(() => layoutState(page, key)).toEqual(placed);
    await expect(custom).not.toBeChecked();

    // Each mode has its own drawing (PD-OI-086): the Circular panel's Separate
    // Strands is the Circular drawing's, so changing it asks nothing and leaves
    // the Linear Separate Strands and lane placement as they are.
    const linearDraft = () => page.evaluate(async (overrideKey) => {
      const { state } = await import('./js/state.js');
      const drawing = state.drawingOf('linear');
      return { separate: drawing.form.separate_strands,
        row: drawing.featurePlacementOverrides[overrideKey]?.placement?.side || null };
    }, key);
    expect(await linearDraft()).toEqual({ separate: false, row: 'above' });
    await switchMode(page, 'circular');
    const circularStrands = await reveal(page.locator('input[type=checkbox][aria-label="Separate Strands"]'));
    const circularBefore = await circularStrands.isChecked();
    await circularStrands.click();
    await expect(dialog(page)).toBeHidden();
    await expect(circularStrands).toBeChecked({ checked: !circularBefore });
    expect(await linearDraft()).toEqual({ separate: false, row: 'above' });
    await switchMode(page, 'linear');
    expect(await layoutState(page, key)).toMatchObject({ separate: false, linear: { side: 'above' }, row: 'above' });
    await generate(page);
  } finally {
    await page.context().close();
  }
});

test('Circular Row actions run from the menu and by keyboard, and focus stays with the row', async ({ browser }) => {
  test.setTimeout(300000);
  const page = await load(browser, 'gbdraw/web/gallery/sessions/HmmtDNA_ATskew.gbdraw-session.json');
  try {
    await openStack(page, 'circular');
    // The features row sits on the Axis, so the Axis marker is inside that row.
    expect(await stackOrder(page, 'circular')).toEqual(['features', 'gc_content', 'gc_skew', 'a_skew_2', 'ticks']);
    const gc = slotRow(page, 'circular', 'gc_content');
    // Move up and Move down stay the row's buttons; the other actions are menu items.
    await expect(gc.getByRole('button', { name: 'Circular track slot gc_content Move outward', exact: true })).toBeVisible();
    await expect(gc.getByRole('button', { name: 'Circular track slot gc_content Move inward', exact: true })).toBeVisible();
    await expect(menuItems(gc)).toHaveText(['Move outside Axis', 'Move inside Axis', 'Duplicate', 'Remove']);
    await expect(menuItem(gc, 'Move inside Axis')).toBeDisabled();

    await runRowAction(gc, 'Move outside Axis');
    await expect(rowMenu(gc)).not.toHaveAttribute('open', '');
    await expect(menuToggle(gc)).toBeFocused();
    expect(await stackOrder(page, 'circular')).toEqual(['gc_content', 'features', 'gc_skew', 'a_skew_2', 'ticks']);

    // Keyboard: Arrow Down opens the menu on its first enabled item, the
    // arrows move through the enabled items, Escape closes it, Enter runs one.
    await page.keyboard.press('ArrowDown');
    const inside = rowMenu(gc).getByRole('button', { name: 'Move inside Axis', exact: true });
    await expect(inside).toBeFocused();
    await expect(rowMenu(gc).getByRole('button', { name: 'Move outside Axis', exact: true })).toBeDisabled();
    await page.keyboard.press('ArrowDown');
    await expect(rowMenu(gc).getByRole('button', { name: 'Duplicate', exact: true })).toBeFocused();
    await page.keyboard.press('ArrowUp');
    await expect(inside).toBeFocused();
    await page.keyboard.press('Escape');
    await expect(rowMenu(gc)).not.toHaveAttribute('open', '');
    await expect(menuToggle(gc)).toBeFocused();
    await page.keyboard.press('ArrowDown');
    await expect(inside).toBeFocused();
    await page.keyboard.press('Enter');
    await expect(menuToggle(gc)).toBeFocused();
    expect(await stackOrder(page, 'circular')).toEqual(['features', 'gc_content', 'gc_skew', 'a_skew_2', 'ticks']);

    await runRowAction(gc, 'Duplicate');
    await expect(menuToggle(gc)).toBeFocused();
    expect(await stackOrder(page, 'circular')).toEqual(['features', 'gc_content', 'gc_content_2', 'gc_skew', 'a_skew_2', 'ticks']);
    // Remove: focus goes to the row now in its place.
    await runRowAction(slotRow(page, 'circular', 'gc_content_2'), 'Remove');
    await expect(menuToggle(slotRow(page, 'circular', 'gc_skew'))).toBeFocused();
    // The last row's Remove focuses the row before it.
    await runRowAction(slotRow(page, 'circular', 'ticks'), 'Remove');
    await expect(menuToggle(slotRow(page, 'circular', 'a_skew_2'))).toBeFocused();
    expect(await stackOrder(page, 'circular')).toEqual(['features', 'gc_content', 'gc_skew', 'a_skew_2']);

    // Track on Axis is offered where the renderer can sit on the Axis.
    const features = slotRow(page, 'circular', 'features');
    await expect(menuItem(features, 'Track on Axis')).toHaveCount(1);
    await expect(menuItem(gc, 'Track on Axis')).toHaveCount(0);
    await expect(menuItem(features, 'Remove')).toHaveClass(/action-menu-item-danger/);
  } finally {
    await page.context().close();
  }
});

test('Circular Track on Axis runs from the Row actions menu', async ({ browser }) => {
  test.setTimeout(300000);
  const page = await load(browser, 'gbdraw/web/gallery/sessions/HmmtDNA_ATskew.gbdraw-session.json');
  try {
    await openStack(page, 'circular');
    const features = slotRow(page, 'circular', 'features');
    const ticks = slotRow(page, 'circular', 'ticks');
    await expect(menuItem(features, 'Track on Axis')).toBeDisabled();
    await runRowAction(ticks, 'Track on Axis');
    await expect(menuToggle(ticks)).toBeFocused();
    await expect(menuItem(ticks, 'Track on Axis')).toBeDisabled();
    await expect(menuItem(features, 'Track on Axis')).toBeEnabled();
    expect(await page.evaluate(() => {
      const app = window.__GBDRAW_APP__;
      return app.adv.circular_track_slots[app.adv.circular_track_slots_axis_index]?.id;
    })).toBe('ticks');
  } finally {
    await page.context().close();
  }
});

test('Linear Row actions run from the menu, and Delete moves focus to the next row, then the row before', async ({ browser }) => {
  test.setTimeout(300000);
  const page = await load(browser, LAMBDA);
  try {
    await openStack(page, 'linear');
    await page.getByRole('checkbox', { name: 'Use custom stack', exact: true }).check();
    if (await page.locator('[data-capture^="linear-track-slot-"]').count() === 0) {
      await page.locator('button[title^="Replace the current custom stack from the simple controls"]').click();
    }
    expect(await stackOrder(page, 'linear')).toEqual(['features', 'gc_content', 'gc_skew']);
    const gc = slotRow(page, 'linear', 'gc_content');
    await expect(gc.getByRole('button', { name: 'Linear track slot gc_content Move up within side', exact: true })).toBeVisible();
    await expect(menuItems(gc)).toHaveText(['Move above Axis', 'Move below Axis', 'Duplicate', 'Delete']);
    // A disabled item names its reason in its title, as the row's button did.
    const duplicateFeatures = menuItem(slotRow(page, 'linear', 'features'), 'Duplicate');
    await expect(duplicateFeatures).toBeDisabled();
    await expect(duplicateFeatures).toHaveAttribute('title', 'Duplicating this row would make the custom stack invalid');
    await runRowAction(gc, 'Move above Axis');
    await expect(menuToggle(gc)).toBeFocused();
    expect(await stackOrder(page, 'linear')).toEqual(['gc_content', 'features', 'gc_skew']);
    await runRowAction(gc, 'Move below Axis');
    expect(await stackOrder(page, 'linear')).toEqual(['features', 'gc_content', 'gc_skew']);
    await runRowAction(gc, 'Duplicate');
    expect(await stackOrder(page, 'linear')).toEqual(['features', 'gc_content', 'gc_content_2', 'gc_skew']);
    await expect(menuItem(gc, 'Delete')).toHaveClass(/action-menu-item-danger/);
    await runRowAction(slotRow(page, 'linear', 'gc_content_2'), 'Delete');
    await expect(menuToggle(slotRow(page, 'linear', 'gc_skew'))).toBeFocused();
    await runRowAction(slotRow(page, 'linear', 'gc_skew'), 'Delete');
    await expect(menuToggle(gc)).toBeFocused();
    await runRowAction(gc, 'Delete');
    await expect(menuToggle(slotRow(page, 'linear', 'features'))).toBeFocused();
    // Deleting the last row refills the stack with the default rows, so a row
    // is always in its place; Add track takes focus only when no row is left.
    await runRowAction(slotRow(page, 'linear', 'features'), 'Delete');
    expect(await stackOrder(page, 'linear')).toEqual(['features', 'gc_content', 'gc_skew']);
    await expect(menuToggle(slotRow(page, 'linear', 'features'))).toBeFocused();
  } finally {
    await page.context().close();
  }
});

// R10: removing a Depth series is a reconcile that keeps the feature row on
// its side of the Axis, also when the removed row sat outside the Axis.
test('removing a Circular Depth series outside the Axis keeps the feature row inside', async ({ page }) => {
  await openApp(page, { waitForPalette: false });
  const stack = () => page.evaluate(() => {
    const app = window.__GBDRAW_APP__;
    const slots = app.adv.circular_track_slots;
    const axis = app.adv.circular_track_slots_axis_index;
    return slots.map((slot, index) => `${slot.id}:${index < axis ? 'outside' : 'inside'}:${slot.params?.lane_direction || ''}`);
  });
  await page.evaluate(() => {
    const app = window.__GBDRAW_APP__;
    const depth = (name) => new File(['position\tdepth\n1\t5\n'], name, { type: 'text/tab-separated-values' });
    app.files.c_depth = [[depth('a.tsv'), depth('b.tsv')]];
    Object.assign(app.form, { show_depth: true, track_type: 'tuckin' });
    app.resetCircularTrackSlotsFromSimpleControls();
    app.setCircularTrackSlotsEnabled(true);
    app.moveCircularTrackSlotOutside(app.adv.circular_track_slots.findIndex((slot) => slot.id === 'depth_1'));
  });
  const before = await stack();
  expect(before).toContain('depth_1:outside:');
  expect(before).toContain('features:inside:inside');
  await page.evaluate(() => window.__GBDRAW_APP__.removeCircularDepthTrack(0));
  const kept = (rows) => rows.filter((row) => !row.startsWith('depth'));
  expect(kept(await stack())).toEqual(kept(before));
});

// OV-1263 (audit G5b-03): a fresh Circular custom stack (every row inside the Axis) with an
// Annotations row set to Overlay. Generate keeps the Axis where the other rows put it; before,
// it took the Overlay Annotations row for the track on the Axis, and Outside axis then moved
// the Features, Ticks, and GC rows outside instead of the Annotations row.
test('Generate with an Overlay Annotations row keeps the Axis, and Outside axis moves only that row', async ({ page }) => {
  test.setTimeout(240000);
  const data = 'gbdraw/web/tutorial-data/tobacco-plastome-regions/';
  await openApp(page);
  await page.getByRole('button', { name: 'Circular', exact: true }).click();
  await page.getByLabel('GenBank/DDBJ File', { exact: true }).first().setInputFiles(`${data}NC_001879.gbk`);
  await page.locator('[data-setting="region-annotations"] input[type=file]').first()
    .setInputFiles(`${data}nicotiana-tabacum-regions.tsv`);
  await expect.poll(() => page.evaluate(async () => (await import('./js/state.js')).state.activeDrawing().annotationSets.length)).toBe(1);
  const stack = () => page.evaluate(async () => {
    const { state } = await import('./js/state.js');
    const { adv } = state.activeDrawing();
    return { axis: adv.circular_track_slots_axis_index, sides: Object.fromEntries(adv.circular_track_slots.map((slot) => [slot.id, slot.side])) };
  });
  const place = (side) => page.evaluate(async (placement) => {
    const { state } = await import('./js/state.js');
    const row = state.activeDrawing().adv.circular_track_slots.find((slot) => slot.renderer === 'annotations');
    return window.__GBDRAW_APP__.updateCircularTrackSlotPlacement(row, placement);
  }, side);
  await page.evaluate(() => {
    const app = window.__GBDRAW_APP__;
    app.setCircularTrackSlotsEnabled(true);
    app.addCircularTrackSlot('annotations');
  });
  await place('overlay');
  // An anchor under the row (the editor gives it the z its layer needs, OV-1297)
  // and Compress, so the Generate draws (two regions share a lane under Error).
  await page.evaluate(async () => {
    const { state } = await import('./js/state.js');
    const row = state.activeDrawing().adv.circular_track_slots.find((slot) => slot.renderer === 'annotations');
    window.__GBDRAW_APP__.setCircularOverlayAnnotationParam(row, 'anchor_slot', 'features');
    row.params.overflow = 'compress';
  });
  const before = await stack();
  expect(before.sides.annotations).toBe('overlay');
  await evaluateWithRetainedPromise(page, () => window.__GBDRAW_APP__.runAnalysis());
  expect(await page.evaluate(async () => {
    const { state } = await import('./js/state.js');
    return { error: state.errorLog.value?.code ?? null, results: state.results.value.length };
  })).toEqual({ error: null, results: 1 });
  expect((await stack()).axis).toBe(before.axis);
  await place('outside');
  expect((await stack()).sides).toEqual({ ...before.sides, annotations: 'outside' });
});

// PD-OI-104 scenario revision 7: removing a Depth file or series, or a stack
// row, that another row starts at or extends to, and Circular Hide GC Content
// with Use custom stack off, ask with one dialog. Cancel changes nothing and
// records no step; the choice is one Undo step, its Undo and Redo restore the
// source, the row, and the switches together, and Generate then draws.
const TOBACCO = 'gbdraw/web/gallery/sessions/tobacco-chloroplast.gbdraw-session.json';
const REGIONS_TSV = 'set_id\tid\tmark\trecord\tstart\tend\nmarks\tm1\thighlight\tNC_001416.1\t100\t900\n';
const depthTsv = (record) => [0, 1, 2, 3].map((index) => `${record}\t${index * 700 + 1}\t${10 + index}`).join('\n');
const rev7Depth = (mode, count) => ({ mode, count, tsv: depthTsv(mode === 'circular' ? 'NC_001879.2' : 'NC_001416.1') });

// A reach row of the shown stack: `anchor` a slot id; `to` an Extend to slot.
const addReachRow = (page, row) => evaluateWithRetainedPromise(page, async ({ id, anchor, to, setId, mode }) => {
  const { state } = await import('./js/state.js');
  await window.__GBDRAW_HISTORY__.runUndoable('Setup', () => {
    const rows = state.activeDrawing().adv[`${mode}_track_slots`];
    rows.push({ id, renderer: 'annotations', enabled: true, side: 'overlay', z: to ? -1 : 1, ...(mode === 'linear' ? { height: '', spacing: '' } : {}),
      params: { set_id: setId, anchor_slot: anchor, layer: to ? 'underlay' : 'foreground', marks: ['highlight'], ...(to ? { to_slot: to } : {}) } });
  });
}, row);
const addDepthSeries = (page, spec) => evaluateWithRetainedPromise(page, async ({ mode, count, tsv }) => {
  const app = window.__GBDRAW_APP__;
  const file = (name) => new File([tsv], name, { type: 'text/tab-separated-values' });
  await window.__GBDRAW_HISTORY__.runUndoable('Setup', () => {
    for (let index = 0; index < count; index += 1) {
      if (mode === 'circular') app.setCircularDepthFile(index, file(`s${index + 1}.tsv`));
      else app.setLinearDepthFile(app.linearSeqs[0], index, file(`s${index + 1}.tsv`));
    }
  });
}, spec);
const openReachFixture = async (browser, mode) => {
  const page = await load(browser, mode === 'circular' ? TOBACCO : LAMBDA);
  if (mode === 'linear') {
    await page.locator('[data-setting="region-annotations"] input[type=file]').first()
      .setInputFiles({ name: 'regions.tsv', mimeType: 'text/tab-separated-values', buffer: Buffer.from(REGIONS_TSV) });
    await expect.poll(() => page.evaluate(async () => (await import('./js/state.js')).state.activeDrawing().annotationSets.length)).toBe(1);
    await page.evaluate(() => window.__GBDRAW_APP__.setLinearTrackSlotsEnabled(true));
    await addReachRow(page, { id: 'tint', anchor: 'features', to: 'gc_content', setId: 'marks', mode });
  }
  return page;
};
const setIdOf = (mode) => (mode === 'circular' ? 'plastome_tint' : 'marks');
const stackState = (page, mode) => page.evaluate(async (stackMode) => {
  const { state } = await import('./js/state.js');
  const app = window.__GBDRAW_APP__;
  const { adv, form } = state.activeDrawing();
  const depth = stackMode === 'circular' ? [app.files.c_depth] : app.linearSeqs.map((seq) => seq.depth);
  return {
    rows: adv[`${stackMode}_track_slots`].map((row) => `${row.id}:${row.enabled !== false}`),
    depth: JSON.stringify(depth, (_, value) => (value instanceof File ? value.name : value)),
    series: adv.depth_tracks.length,
    showDepth: form.show_depth, suppressGc: form.suppress_gc ?? null
  };
}, mode);
const undoCount = (page) => page.evaluate(() => window.__GBDRAW_HISTORY__.getUndoCount());
const historyStep = async (page, step) => {
  await evaluateWithRetainedPromise(page, (name) => window.__GBDRAW_HISTORY__[name](), step);
  await settleLive(page);
};
const expectGenerates = async (page) => {
  await evaluateWithRetainedPromise(page, () => window.__GBDRAW_APP__.runAnalysis());
  expect(await page.evaluate(async () => {
    const { state } = await import('./js/state.js');
    return { error: state.errorLog.value?.code ?? null, results: state.results.value.length };
  })).toEqual({ error: null, results: 1 });
};
// Cancel, the choice, Undo, Redo, Generate. `trigger` makes the real click;
// `focused` is the control that holds the focus after Cancel.
const exerciseAsk = async (page, mode, { trigger, title, focused, expectAfter, beforeGenerate = async () => {}, live = false, cancel = 'Cancel', generates = true }) => {
  const prompt = page.locator('[data-choice-prompt]');
  const heading = prompt.getByRole('heading', { name: title, exact: true });
  const before = await stackState(page, mode);
  const steps = await undoCount(page);
  await trigger();
  await expect(heading).toBeVisible();
  await prompt.getByRole('button', { name: cancel, exact: true }).click();
  await expect(heading).toHaveCount(0);
  await settleLive(page);
  expect(await stackState(page, mode), 'Cancel changes nothing').toEqual(before);
  expect(await undoCount(page), 'Cancel records no step').toBe(steps);
  if (focused) await expect(focused()).toBeFocused();
  await trigger();
  await expect(heading).toBeVisible();
  await prompt.locator('[data-choice-prompt-action="hide"]').click();
  await expect(heading).toHaveCount(0);
  await settleLive(page);
  const after = await stackState(page, mode);
  expectAfter(after);
  expect(await undoCount(page), 'one Undo step').toBe(steps + 1);
  await historyStep(page, 'undo');
  expect(await stackState(page, mode), 'Undo restores it together').toEqual(before);
  await historyStep(page, 'redo');
  expect(await stackState(page, mode), 'Redo makes it again').toEqual(after);
  if (live) await expectLiveEqualsGenerate(page, { label: title });
  if (!generates) return;
  await beforeGenerate();
  await expectGenerates(page);
};
const rowNumber = (state, id) => state.rows.findIndex((row) => row.startsWith(`${id}:`)) + 1;
const openDepthSources = async (page, mode) => {
  const section = page.locator(mode === 'circular' ? '[data-setting="depth-sources"]' : '[data-linear-depth-settings]').first();
  await reveal(section);
  if (!(await section.evaluate((node) => node.open))) await section.locator(':scope > summary').click();
  return section;
};
const removeDepthFile = async (page, mode) => {
  if (mode === 'linear') {
    const source = page.locator('[data-linear-source-depth]').first();
    await reveal(source);
    if (!(await source.evaluate((node) => node.open))) await source.locator(':scope > summary').click();
    return source.getByRole('button', { name: /^Remove Depth TSV/ }).first().click();
  }
  const section = await openDepthSources(page, mode);
  return section.getByRole('button', { name: /^Remove Depth TSV/ }).first().click();
};
const expectDepthLeft = (after) => {
  expect(after.rows).toContain('mark:false');
  expect(after.rows.some((row) => row.startsWith('depth:'))).toBe(false);
  expect(after.depth).not.toContain('s1.tsv');
  expect(after.showDepth).toBe(false);
};
const reflowError = (page) => page.evaluate(async () => (await import('./js/state.js')).state.labelReflowLastError.value?.summary ?? null);
// A drawing with one Depth file, its row reached by `mark`, generated.
const generatedWithReach = async (browser, mode) => {
  const page = await openReachFixture(browser, mode);
  await addDepthSeries(page, rev7Depth(mode, 1));
  await addReachRow(page, { id: 'mark', anchor: 'depth', setId: setIdOf(mode), mode });
  await expectGenerates(page);
  await settleLive(page);
  return page;
};
// The Session a Save writes for that drawing, loaded: its Result draws Depth,
// and its Web input inventory names the Depth file.
const loadedWithReach = async (browser, mode, testInfo) => {
  const saved = await generatedWithReach(browser, mode);
  const file = testInfo.outputPath(`saved-${mode}.gbdraw-session.json`);
  await savedSession(saved, file);
  await saved.context().close();
  return load(browser, file);
};
const choose = async (page, title) => {
  const prompt = page.locator('[data-choice-prompt]');
  await expect(prompt.getByRole('heading', { name: title, exact: true })).toBeVisible();
  await prompt.locator('[data-choice-prompt-action="hide"]').click();
  await expect(prompt).toHaveCount(0);
  await settleLive(page);
};
const drawsDepth = (page) => page.evaluate(() => window.__GBDRAW_APP__.svgContainer.querySelector('svg').querySelectorAll('[id*="depth" i]').length > 0);

for (const mode of ['circular', 'linear']) {
  test(`rev 7: ${mode} Depth file removal after Generate asks, the row leaves, and the live Result equals Generate`, async ({ browser }) => {
    test.setTimeout(420000);
    const page = await generatedWithReach(browser, mode);
    try {
      const start = await stackState(page, mode);
      await exerciseAsk(page, mode, {
        trigger: () => removeDepthFile(page, mode),
        title: `Remove the Depth file and turn off Track row ${rowNumber(start, 'mark')}?`,
        live: true,
        expectAfter: expectDepthLeft
      });
    } finally { await page.context().close(); }
  });

  test(`rev 7: ${mode} Depth file removal after Load asks, removes the drawn track at once, and its Undo shows it again`, async ({ browser }, testInfo) => {
    test.setTimeout(420000);
    const page = await loadedWithReach(browser, mode, testInfo);
    try {
      const start = await stackState(page, mode);
      expect(await drawsDepth(page)).toBe(true);
      await removeDepthFile(page, mode);
      await choose(page, `Remove the Depth file and turn off Track row ${rowNumber(start, 'mark')}?`);
      expectDepthLeft(await stackState(page, mode));
      expect(await reflowError(page)).toBeNull();
      expect(await drawsDepth(page), 'the drawn Depth track left at once').toBe(false);
      await historyStep(page, 'undo');
      expect(await stackState(page, mode)).toEqual(start);
      expect(await reflowError(page)).toBeNull();
      expect(await drawsDepth(page), 'its Undo shows it again').toBe(true);
      await historyStep(page, 'redo');
      expect(await drawsDepth(page)).toBe(false);
      await expectLiveEqualsGenerate(page, { label: `${mode} file after Load` });
      await expectGenerates(page);
    } finally { await page.context().close(); }
  });

  test(`rev 7: ${mode} removing the Depth row a loaded Result drew removes the track at once, equal to Generate`, async ({ browser }, testInfo) => {
    test.setTimeout(420000);
    const page = await loadedWithReach(browser, mode, testInfo);
    try {
      await openStack(page, mode);
      const start = await stackState(page, mode);
      await runRowAction(slotRow(page, mode, 'depth'), mode === 'circular' ? 'Remove' : 'Delete');
      await choose(page, `Remove Track row ${rowNumber(start, 'depth')} and turn off Track row ${rowNumber(start, 'mark')}?`);
      expect(await reflowError(page)).toBeNull();
      expect(await drawsDepth(page)).toBe(false);
      await expectLiveEqualsGenerate(page, { label: `${mode} row after Load` });
    } finally { await page.context().close(); }
  });

  test(`rev 7: ${mode} Depth file removal writes only its own fields, with other pending values`, async ({ browser }) => {
    test.setTimeout(420000);
    const page = await openReachFixture(browser, mode);
    try {
      await addDepthSeries(page, rev7Depth(mode, 1));
      await addReachRow(page, { id: 'mark', anchor: 'depth', setId: setIdOf(mode), mode });
      const start = await stackState(page, mode);
      const pending = () => page.evaluate(async (stackMode) => {
        const { state } = await import('./js/state.js');
        const { form, adv } = state.activeDrawing();
        const { show_depth: _showDepth, ...otherForm } = form;
        const { depth_tracks: _series, [`${stackMode}_track_slots`]: _slots, [`${stackMode}_track_slots_axis_index`]: _axis, ...otherAdv } = adv;
        return JSON.stringify({ otherForm, otherAdv, untouched: adv[`${stackMode}_track_slots`].filter((row) => row.id !== 'depth' && row.id !== 'mark') });
      }, mode);
      await page.evaluate(async () => {
        const { state } = await import('./js/state.js');
        const drawing = state.activeDrawing();
        drawing.adv.label_font_size = 29;
        drawing.form.legend = 'bottom';
      });
      const before = await pending();
      await removeDepthFile(page, mode);
      await choose(page, `Remove the Depth file and turn off Track row ${rowNumber(start, 'mark')}?`);
      expect(await pending()).toBe(before);
    } finally { await page.context().close(); }
  });

  test(`rev 7: ${mode} removing the last Depth file asks, and its Undo restores the file and the row`, async ({ browser }) => {
    test.setTimeout(420000);
    const page = await openReachFixture(browser, mode);
    try {
      await addDepthSeries(page, rev7Depth(mode, 1));
      await addReachRow(page, { id: 'mark', anchor: 'depth', setId: setIdOf(mode), mode });
      const start = await stackState(page, mode);
      const name = `Track row ${rowNumber(start, 'mark')}`;
      await exerciseAsk(page, mode, {
        trigger: () => removeDepthFile(page, mode),
        title: `Remove the Depth file and turn off ${name}?`,
        expectAfter: (after) => {
          expect(after.rows).toContain('mark:false');
          expect(after.rows.some((row) => row.startsWith('depth:true'))).toBe(false);
          expect(after.depth).not.toContain('s1.tsv');
          expect(after.showDepth).toBe(false);
        }
      });
    } finally { await page.context().close(); }
  });

  test(`rev 7: ${mode} removing one Depth series of two asks, and keeps the other series`, async ({ browser }) => {
    test.setTimeout(420000);
    const page = await openReachFixture(browser, mode);
    try {
      await addDepthSeries(page, rev7Depth(mode, 2));
      const seriesRow = await page.evaluate(async (stackMode) => {
        const { state } = await import('./js/state.js');
        return state.activeDrawing().adv[`${stackMode}_track_slots`].find((row) => row.renderer === 'depth' && row.params.track_index === 0).id;
      }, mode);
      await addReachRow(page, { id: 'mark', anchor: seriesRow, setId: setIdOf(mode), mode });
      const start = await stackState(page, mode);
      const button = async () => {
        const section = await openDepthSources(page, mode);
        return section.getByRole('button', { name: mode === 'circular' ? 'Remove' : 'Remove series', exact: true }).first();
      };
      await exerciseAsk(page, mode, {
        trigger: async () => (await button()).click(),
        title: `Remove Depth series #1 and turn off Track row ${rowNumber(start, 'mark')}?`,
        expectAfter: (after) => {
          expect(after.rows).toContain('mark:false');
          expect(after.rows.some((row) => row.startsWith(`${seriesRow}:`))).toBe(false);
          expect(after.rows.filter((row) => row.startsWith('depth'))).toHaveLength(1);
          expect(after.depth).toContain('s2.tsv');
          expect(after.depth).not.toContain('s1.tsv');
          expect(after.showDepth).toBe(true);
        }
      });
    } finally { await page.context().close(); }
  });

  test(`rev 7: ${mode} removing a stack row another row reaches asks through the layout transition`, async ({ browser }) => {
    test.setTimeout(420000);
    const page = await openReachFixture(browser, mode);
    try {
      await openStack(page, mode);
      const reachId = mode === 'circular' ? 'plastome_tint' : 'tint';
      const start = await stackState(page, mode);
      await exerciseAsk(page, mode, {
        trigger: () => runRowAction(slotRow(page, mode, 'gc_content'), mode === 'circular' ? 'Remove' : 'Delete'),
        title: `Remove Track row ${rowNumber(start, 'gc_content')} and turn off Track row ${rowNumber(start, reachId)}?`,
        focused: () => menuToggle(slotRow(page, mode, 'gc_content')),
        live: mode === 'circular',
        expectAfter: (after) => {
          expect(after.rows.some((row) => row.startsWith('gc_content:'))).toBe(false);
          expect(after.rows).toContain(`${reachId}:false`);
        }
      });
    } finally { await page.context().close(); }
  });
}

test('rev 7: Circular Hide GC Content with Use custom stack off asks, and the stack turned on after the choice draws', async ({ browser }) => {
  test.setTimeout(420000);
  const page = await load(browser, TOBACCO);
  try {
    await page.evaluate(() => window.__GBDRAW_APP__.setCircularTrackSlotsEnabled(false));
    await openStack(page, 'circular');
    const box = await reveal(page.getByRole('checkbox', { name: 'Hide GC Content', exact: true }));
    const start = await stackState(page, 'circular');
    await exerciseAsk(page, 'circular', {
      trigger: () => box.click(),
      title: 'Hide GC content?',
      cancel: 'Keep GC content',
      focused: () => box,
      expectAfter: (after) => {
        expect(after.suppressGc).toBe(true);
        expect(after.rows).toContain('gc_content:false');
        expect(after.rows).toContain('plastome_tint:false');
      },
      beforeGenerate: () => page.evaluate(() => window.__GBDRAW_APP__.setCircularTrackSlotsEnabled(true))
    });
  } finally { await page.context().close(); }
});

// A request a Python caller or the CLI wrote may hold its slots as spec
// strings. Load makes the draft's rows from them, and the dialog names the row.
test('rev 7: a Session whose request holds string slot specs asks on Remove and on Hide GC Content, by slot id', async ({ browser }, testInfo) => {
  test.setTimeout(420000);
  const session = JSON.parse(readFileSync(TOBACCO, 'utf8'));
  const specOf = (slot) => {
    const params = Object.entries(slot.params || {}).filter(([key]) => key !== 'marks').map(([key, value]) => `${key}=${value}`);
    if (slot.side) params.push(`side=${slot.side}`);
    if (slot.z) params.push(`z=${slot.z}`);
    return `${slot.id}:${slot.renderer}${params.length ? `@${params.join(',')}` : ''}`;
  };
  const tracks = session.renderRequest.diagramOptions.tracks;
  tracks.circularTrackSlots = tracks.circularTrackSlots.map(specOf);
  delete session.modes.circular.config.adv.circular_track_slots;
  const file = testInfo.outputPath('string-slots.gbdraw-session.json');
  writeFileSync(file, JSON.stringify(session));
  const page = await load(browser, file);
  try {
    await openStack(page, 'circular');
    const start = await stackState(page, 'circular');
    expect(start.rows).toContain('plastome_tint:true');
    const prompt = page.locator('[data-choice-prompt]');
    const title = `Remove Track row ${rowNumber(start, 'gc_content')} and turn off Track row ${rowNumber(start, 'plastome_tint')}?`;
    await runRowAction(slotRow(page, 'circular', 'gc_content'), 'Remove');
    await expect(prompt.getByRole('heading', { name: title, exact: true })).toBeVisible();
    await prompt.locator('[data-choice-prompt-action="hide"]').click();
    await settleLive(page);
    expect((await stackState(page, 'circular')).rows).toContain('plastome_tint:false');
    await historyStep(page, 'undo');
    expect((await stackState(page, 'circular')).rows).toEqual(start.rows);
    const box = await reveal(page.getByRole('checkbox', { name: 'Hide GC Content', exact: true }));
    await box.click();
    await expect(prompt.getByRole('heading', { name: 'Hide GC content?', exact: true })).toBeVisible();
    await prompt.locator('[data-choice-prompt-action="hide"]').click();
    await settleLive(page);
    expect((await stackState(page, 'circular')).rows).toContain('plastome_tint:false');
    await expectGenerates(page);
  } finally { await page.context().close(); }
});

// PD-OI-058: a managed Depth row leaves with its series' last source, also
// after a Generate. No row reaches it here, so nothing asks.
const expectNoQuestion = async (page) => {
  await settleLive(page);
  await expect(page.locator('[data-choice-prompt]')).toHaveCount(0);
};
const linearDepthRows = (page) => page.evaluate(async () => (await import('./js/state.js')).state.activeDrawing().adv.linear_track_slots
  .filter((row) => row.renderer === 'depth').map((row) => ({ id: row.id, height: row.height })));
const savedSession = async (saved, file) => {
  const titled = await saved.evaluate(async () => Boolean((await import('./js/state.js')).state.sessionTitle?.value));
  await download(saved, 'Save Session', file, titled ? {} : { title: 'rev7' });
  return JSON.parse(gunzipSync(readFileSync(file)).toString('utf8'));
};
const heightsInSession = (session) => {
  const depthRows = (rows) => rows.filter((row) => row.renderer === 'depth').map((row) => row.height);
  return {
    draft: depthRows(session.modes.linear.config.adv.linear_track_slots),
    request: depthRows(session.renderRequest.diagramOptions.tracks.linearTrackSlots)
  };
};
const generatedLinear = async (browser, { reach = false } = {}) => {
  const page = await openReachFixture(browser, 'linear');
  await addDepthSeries(page, rev7Depth('linear', 1));
  if (reach) await addReachRow(page, { id: 'mark', anchor: 'depth', setId: 'marks', mode: 'linear' });
  await expectGenerates(page);
  await settleLive(page);
  return page;
};
// That drawing saved and loaded, its Depth row holding the height Generate
// writes, as `dev` saves it.
const loadedWithWrittenHeight = async (browser, testInfo, options) => {
  const saved = await generatedLinear(browser, options);
  const file = testInfo.outputPath('written-height.gbdraw-session.json');
  const session = await savedSession(saved, file);
  await saved.context().close();
  const writeHeight = (height) => (row) => (row.renderer === 'depth' ? { ...row, height } : row);
  const { adv } = session.modes.linear.config;
  adv.linear_track_slots = adv.linear_track_slots.map(writeHeight('10px'));
  const { tracks } = session.renderRequest.diagramOptions;
  tracks.linearTrackSlots = tracks.linearTrackSlots.map(writeHeight({ value: 10, unit: 'px' }));
  writeFileSync(file, JSON.stringify(session));
  const page = await load(browser, file);
  expect(await linearDepthRows(page)).toEqual([{ id: 'depth', height: '10px' }]);
  return page;
};
// What decides the drawing: the request's options, as a Save writes them, and
// the drawn size. A later Generate's SVG text is reserialized, so its bytes differ.
const drawnFrom = async (page, file) => ({
  options: (await savedSession(page, file)).renderRequest.diagramOptions,
  viewBox: await page.evaluate(() => window.__GBDRAW_APP__.svgContainer.querySelector('svg').getAttribute('viewBox').split(' ').map(Number))
});

for (const viaLoad of [false, true]) {
  const when = viaLoad ? 'after Load of a Session whose row holds the written height' : 'after Generate';
  const open = (browser, testInfo) => (viaLoad ? loadedWithWrittenHeight(browser, testInfo) : generatedLinear(browser));

  test(`rev 7: Linear Depth file removal with no row reaching it ${when} takes the row out without a question, and Generate draws`, async ({ browser }, testInfo) => {
    test.setTimeout(420000);
    const page = await open(browser, testInfo);
    try {
      await removeDepthFile(page, 'linear');
      await expectNoQuestion(page);
      expect((await stackState(page, 'linear')).rows.some((row) => row.startsWith('depth:'))).toBe(false);
      await expectGenerates(page);
    } finally { await page.context().close(); }
  });

  test(`rev 7: Linear Depth series removal ${when} takes the row out, and Generate draws`, async ({ browser }, testInfo) => {
    test.setTimeout(420000);
    const page = await open(browser, testInfo);
    try {
      const section = await openDepthSources(page, 'linear');
      await section.getByRole('button', { name: 'Remove series', exact: true }).first().click();
      await expectNoQuestion(page);
      expect(await linearDepthRows(page)).toEqual([]);
      await expectGenerates(page);
    } finally { await page.context().close(); }
  });
}

test('rev 7: Linear Depth file removal after Load of a Session whose row holds the written height asks, and its choice is one Undo step', async ({ browser }, testInfo) => {
  test.setTimeout(420000);
  const page = await loadedWithWrittenHeight(browser, testInfo, { reach: true });
  try {
    const start = await stackState(page, 'linear');
    await exerciseAsk(page, 'linear', {
      trigger: () => removeDepthFile(page, 'linear'),
      title: `Remove the Depth file and turn off Track row ${rowNumber(start, 'mark')}?`,
      expectAfter: expectDepthLeft
    });
  } finally { await page.context().close(); }
});

// A row's height is also written to its series, so it stays when the row leaves.
test('rev 7: a typed Linear Depth height: the row leaves with its last file, the series keeps the value, and the file coming back draws the same height from a managed row', async ({ browser }, testInfo) => {
  test.setTimeout(420000);
  const page = await openReachFixture(browser, 'linear');
  try {
    await addDepthSeries(page, rev7Depth('linear', 1));
    await page.evaluate(async () => {
      const { state } = await import('./js/state.js');
      const row = state.activeDrawing().adv.linear_track_slots.find((slot) => slot.id === 'depth');
      window.__GBDRAW_APP__.setLinearTrackSlotHeight(row, '40');
    });
    await expectGenerates(page);
    await settleLive(page);
    const seriesHeight = () => page.evaluate(async () => (await import('./js/state.js')).state.activeDrawing().adv.depth_tracks[0].height);
    const before = await drawnFrom(page, testInfo.outputPath('typed-before.gbdraw-session.json'));
    expect(before.options.depthTracks[0].height).toBe(40);
    expect(await seriesHeight()).toBe(40);
    await removeDepthFile(page, 'linear');
    await expectNoQuestion(page);
    expect(await linearDepthRows(page)).toEqual([]);
    expect(await seriesHeight()).toBe(40);
    await expectGenerates(page);
    await addDepthSeries(page, rev7Depth('linear', 1));
    expect(await linearDepthRows(page)).toEqual([{ id: 'depth', height: '' }]);
    expect(await seriesHeight()).toBe(40);
    await expectGenerates(page);
    expect(await drawnFrom(page, testInfo.outputPath('typed-after.gbdraw-session.json'))).toEqual(before);
    expect(await linearDepthRows(page)).toEqual([{ id: 'depth', height: '40px' }]);
  } finally { await page.context().close(); }
});

test('rev 7: a Linear Generate keeps the written height in the row, in the Session draft and its request', async ({ browser }, testInfo) => {
  test.setTimeout(420000);
  const page = await generatedLinear(browser);
  try {
    expect(await linearDepthRows(page)).toEqual([{ id: 'depth', height: '10px' }]);
    expect(heightsInSession(await savedSession(page, testInfo.outputPath('kept-height.gbdraw-session.json'))))
      .toEqual({ draft: ['10px'], request: [{ value: 10, unit: 'px' }] });
  } finally { await page.context().close(); }
});

// One File of two records: its Depth controls are the mixed row's Remove from
// all records and each record's own uploader. The records are stubs, so these
// stop before Generate.
const twoRecordFile = async (browser) => {
  const context = await browser.newContext({ viewport: { width: 1600, height: 1000 } });
  const page = await context.newPage();
  await openApp(page, { waitForPalette: false });
  await page.evaluate(() => {
    const app = window.__GBDRAW_APP__;
    app.selectDrawing('linear');
    app.addLinearSeq();
    const gb = new File(['LOCUS       MULTI\n'], 'multi.gbk', { type: 'text/plain' });
    app.linearSeqs[0].gb = gb;
    app.linearSeqs[1].gb = gb;
    app.setLinearTrackSlotsEnabled(true);
    app.setLinearDepthFile(app.linearSeqs[0], 0, new File(['a\t1\t1\n'], 's1.tsv'));
  });
  await addReachRow(page, { id: 'mark', anchor: 'depth', setId: 'marks', mode: 'linear' });
  return page;
};
test('rev 7: Linear Remove from all records of a File asks before it clears the last source of a series', async ({ browser }) => {
  test.setTimeout(420000);
  const page = await twoRecordFile(browser);
  try {
    const source = page.locator('[data-linear-source-depth]').first();
    await reveal(source);
    if (!(await source.evaluate((node) => node.open))) await source.locator(':scope > summary').click();
    const start = await stackState(page, 'linear');
    await exerciseAsk(page, 'linear', {
      trigger: () => source.getByRole('button', { name: 'Remove from all records', exact: true }).click(),
      title: `Remove the Depth file and turn off Track row ${rowNumber(start, 'mark')}?`,
      generates: false,
      expectAfter: expectDepthLeft
    });
  } finally { await page.context().close(); }
});
test('rev 7: Linear per-record uploader asks when it removes the last source of a series', async ({ browser }) => {
  test.setTimeout(420000);
  const page = await twoRecordFile(browser);
  try {
    const perRecord = page.locator('[data-linear-source-records] details[data-setting="depth-sources"]').first();
    await perRecord.evaluate((node) => { for (let open = node; open; open = open.parentElement.closest('details')) open.open = true; });
    await reveal(perRecord);
    const start = await stackState(page, 'linear');
    await exerciseAsk(page, 'linear', {
      trigger: () => perRecord.getByRole('button', { name: /^Remove Depth TSV/ }).first().click(),
      title: `Remove the Depth file and turn off Track row ${rowNumber(start, 'mark')}?`,
      generates: false,
      expectAfter: expectDepthLeft
    });
  } finally { await page.context().close(); }
});

test('rev 7: Circular Hide GC Skew with Use custom stack off asks when it turns off a row another row reaches', async ({ browser }) => {
  test.setTimeout(420000);
  const page = await load(browser, TOBACCO);
  try {
    await evaluateWithRetainedPromise(page, async () => {
      const { state } = await import('./js/state.js');
      await window.__GBDRAW_HISTORY__.runUndoable('Setup', () => {
        const drawing = state.activeDrawing();
        drawing.form.suppress_skew = false;
        // The rows as the editor's normalizer leaves them, so the question has nothing to normalize.
        const auto = { width: null, radius: null, inner_gap_px: null, outer_gap_px: null };
        drawing.adv.circular_track_slots.push(
          { id: 'gc_skew', renderer: 'dinucleotide_skew', enabled: true, side: 'inside', z: 0, params: {}, ...auto },
          { id: 'mark', renderer: 'annotations', enabled: true, side: 'overlay', z: 1,
            params: { set_id: 'plastome_tint', anchor_slot: 'gc_skew', layer: 'foreground', marks: ['highlight'], overflow: 'error', show_labels: true }, ...auto }
        );
      });
    });
    await page.evaluate(() => window.__GBDRAW_APP__.setCircularTrackSlotsEnabled(false));
    await openStack(page, 'circular');
    const box = await reveal(page.getByRole('checkbox', { name: 'Hide GC Skew', exact: true }));
    const start = await stackState(page, 'circular');
    await exerciseAsk(page, 'circular', {
      trigger: () => box.click(),
      title: 'Hide GC skew?',
      cancel: 'Keep GC skew',
      focused: () => box,
      expectAfter: (after) => {
        expect(after.rows).toContain('gc_skew:false');
        expect(after.rows).toContain('mark:false');
      },
      beforeGenerate: () => page.evaluate(() => window.__GBDRAW_APP__.setCircularTrackSlotsEnabled(true))
    });
  } finally { await page.context().close(); }
});

// Undo of a chosen removal returns the rows that reached it to the
// Result and to the committed request, and Redo takes them out again; a Save
// after the Undo holds them. (A row switch has its own cases in
// live-generate-parity-track-data.playwright.spec.js.)
const committedSlotsOf = (page, mode) => page.evaluate(async (stackMode) => {
  const { getCommittedCanonicalSession } = await import('./js/services/config.js');
  return JSON.parse(JSON.stringify(getCommittedCanonicalSession().renderRequest.diagramOptions.tracks[`${stackMode}TrackSlots`]));
}, mode);
const slotDrawnCount = (page, id) => page.evaluate((slotId) => window.__GBDRAW_APP__.svgContainer.querySelector('svg')
  .querySelectorAll(`[data-gbdraw-slot-id="${slotId}"] > *`).length, id);
const slotIdsOf = (slots) => slots.map((slot) => (typeof slot === 'string' ? slot.split(':')[0] : slot.id));
const chooseAsked = async (page) => {
  await page.locator('[data-choice-prompt] [data-choice-prompt-action="hide"]').click();
  await settleLive(page);
};
const REACH_UNDO_ROUTES = {
  'Depth file removal': (page, mode) => removeDepthFile(page, mode),
  'Depth series removal': async (page, mode) => {
    const section = await openDepthSources(page, mode);
    await section.getByRole('button', { name: mode === 'circular' ? 'Remove' : 'Remove series', exact: true }).first().click();
  },
  'row removal': async (page, mode) => {
    await openStack(page, mode);
    await runRowAction(slotRow(page, mode, 'depth'), mode === 'circular' ? 'Remove' : 'Delete');
  }
};
for (const mode of ['circular', 'linear']) {
  for (const [route, trigger] of Object.entries(REACH_UNDO_ROUTES)) {
    for (const after of ['Generate', 'Load']) {
      test(`rev 7: Undo and Redo of a chosen ${route} return the reaching row to the Result and the committed request (${mode}, after ${after})`, async ({ browser }, testInfo) => {
        test.setTimeout(420000);
        const page = after === 'Load' ? await loadedWithReach(browser, mode, testInfo) : await generatedWithReach(browser, mode);
        try {
          const drawnBefore = await slotDrawnCount(page, 'mark');
          const slotsBefore = await committedSlotsOf(page, mode);
          expect(drawnBefore).toBeGreaterThan(0);
          await trigger(page, mode);
          await chooseAsked(page);
          expect(await reflowError(page)).toBeNull();
          expect(await slotDrawnCount(page, 'mark'), 'the choice takes it out').toBe(0);
          const slotsAfter = await committedSlotsOf(page, mode);
          expect(slotIdsOf(slotsAfter)).not.toContain('mark');
          await historyStep(page, 'undo');
          expect(await reflowError(page)).toBeNull();
          expect(await slotDrawnCount(page, 'mark'), 'Undo draws it again').toBe(drawnBefore);
          expect(await committedSlotsOf(page, mode), 'Undo gives the committed request back').toEqual(slotsBefore);
          const saved = await savedSession(page, testInfo.outputPath('after-undo.gbdraw-session.json'));
          expect(slotIdsOf(saved.renderRequest.diagramOptions.tracks[`${mode}TrackSlots`]), 'a Save after the Undo holds it').toContain('mark');
          await historyStep(page, 'redo');
          expect(await slotDrawnCount(page, 'mark'), 'Redo takes it out again').toBe(0);
          expect(await committedSlotsOf(page, mode)).toEqual(slotsAfter);
          await historyStep(page, 'undo');
          await expectLiveEqualsGenerate(page, { label: `${mode} ${route} after ${after}` });
        } finally { await page.context().close(); }
      });
    }
  }
}

// PD-OI-104 scenario revision 9: a Depth series' last file.
const promptText = (page) => page.locator('[data-choice-prompt]').innerText();
const setRow = (page, mode, id, fields) => evaluateWithRetainedPromise(page, async ({ stackMode, rowId, values }) => {
  const { state } = await import('./js/state.js');
  await window.__GBDRAW_HISTORY__.runUndoable('Setup', () => {
    Object.assign(state.activeDrawing().adv[`${stackMode}_track_slots`].find((row) => row.id === rowId), values);
  });
}, { stackMode: mode, rowId: id, values: fields });
const addDepthFile = (page, mode, index = 0) => evaluateWithRetainedPromise(page, async ({ stackMode, series, tsv }) => {
  const app = window.__GBDRAW_APP__;
  const file = new File([tsv], `s${series + 1}.tsv`, { type: 'text/tab-separated-values' });
  await window.__GBDRAW_HISTORY__.runUndoable('Setup', () => (stackMode === 'circular' ? app.setCircularDepthFile(series, file) : app.setLinearDepthFile(app.linearSeqs[0], series, file)));
}, { stackMode: mode, series: index, tsv: rev7Depth(mode, 1).tsv });
const generateError = async (page) => {
  await evaluateWithRetainedPromise(page, () => window.__GBDRAW_APP__.runAnalysis());
  return page.evaluate(async () => (await import('./js/state.js')).state.errorLog.value?.code ?? null);
};
const depthRowsOf = (state) => state.rows.filter((row) => /^depth(_\d+)?:/.test(row));

for (const mode of ['circular', 'linear']) {
  for (const when of ['before the first Generate', 'after a Generate', 'after a Load']) {
    test(`rev 7: ${mode} a Depth row turned off stays off without a question when its last file goes and comes back, ${when}, and Generate draws`, async ({ browser }, testInfo) => {
      test.setTimeout(420000);
      let page = await openReachFixture(browser, mode);
      try {
        await addDepthSeries(page, rev7Depth(mode, 1));
        await page.evaluate(async (stackMode) => {
          const { state } = await import('./js/state.js');
          const app = window.__GBDRAW_APP__;
          const row = state.activeDrawing().adv[`${stackMode}_track_slots`].find((slot) => slot.id === 'depth');
          if (stackMode === 'circular') app.setCircularTrackSlotEnabled(row, false); else app.setLinearTrackSlotEnabled(row, false);
        }, mode);
        await settleLive(page);
        if (when !== 'before the first Generate') await expectGenerates(page);
        if (when === 'after a Load') {
          const file = testInfo.outputPath('off.gbdraw-session.json');
          await savedSession(page, file);
          await page.context().close();
          page = await load(browser, file);
          expect(depthRowsOf(await stackState(page, mode))).toEqual(['depth:false']);
        }
        const steps = await undoCount(page);
        const start = await stackState(page, mode);
        await removeDepthFile(page, mode);
        await expectNoQuestion(page);
        const removed = await stackState(page, mode);
        expect(depthRowsOf(removed), 'the row stays, off').toEqual(['depth:false']);
        expect(removed.showDepth).toBe(false);
        expect(await reflowError(page)).toBeNull();
        expect(await undoCount(page)).toBe(steps + 1);
        await historyStep(page, 'undo');
        expect(await stackState(page, mode)).toEqual(start);
        await historyStep(page, 'redo');
        await expectGenerates(page);
        expect(depthRowsOf(await stackState(page, mode))).toEqual(['depth:false']);
        const savedFile = testInfo.outputPath('off-removed.gbdraw-session.json');
        await savedSession(page, savedFile);
        await addDepthFile(page, mode);
        await settleLive(page);
        const back = await stackState(page, mode);
        expect(depthRowsOf(back), 'a re-added file adds no row and does not turn it on').toEqual(['depth:false']);
        expect(back.showDepth).toBe(true);
        await expectGenerates(page);
        expect(await drawsDepth(page), 'the row is off').toBe(false);
        await page.context().close();
        page = await load(browser, savedFile);
        expect(depthRowsOf(await stackState(page, mode)), 'a Session saved in that state keeps it').toEqual(['depth:false']);
        await addDepthFile(page, mode);
        await settleLive(page);
        expect(depthRowsOf(await stackState(page, mode))).toEqual(['depth:false']);
        await expectGenerates(page);
      } finally { await page.context().close(); }
    });
  }

  test(`rev 7: ${mode} a manual Depth row another row reaches asks, and the choice turns it off too`, async ({ browser }) => {
    test.setTimeout(420000);
    const page = await openReachFixture(browser, mode);
    try {
      await addDepthSeries(page, rev7Depth(mode, 1));
      await setRow(page, mode, 'depth', mode === 'circular' ? { id: 'cov' } : { id: 'cov', spacing: '6px' });
      await addReachRow(page, { id: 'mark', anchor: 'cov', setId: setIdOf(mode), mode });
      await expectGenerates(page);
      await settleLive(page);
      const start = await stackState(page, mode);
      const [cov, mark] = [rowNumber(start, 'cov'), rowNumber(start, 'mark')];
      const title = `Remove the Depth file and turn off Track rows ${cov} and ${mark}?`;
      await removeDepthFile(page, mode);
      await expect(page.locator('[data-choice-prompt]').getByRole('heading', { name: title, exact: true })).toBeVisible();
      const text = await promptText(page);
      expect(text).toContain(`Track row ${cov} (cov) has no Depth file left. It stays in the stack, turned off.`);
      expect(text).toContain(`Track row ${mark} (mark)`);
      expect(text).not.toContain('stays empty');
      await page.locator('[data-choice-prompt]').getByRole('button', { name: 'Cancel', exact: true }).click();
      await settleLive(page);
      await exerciseAsk(page, mode, {
        trigger: () => removeDepthFile(page, mode),
        title,
        live: true,
        expectAfter: (after) => {
          expect(after.rows).toEqual(expect.arrayContaining(['cov:false', 'mark:false']));
          expect(after.depth).not.toContain('s1.tsv');
          expect(after.showDepth).toBe(false);
        }
      });
      expect(await reflowError(page)).toBeNull();
    } finally { await page.context().close(); }
  });

  test(`rev 7: ${mode} a manual Depth row no row reaches stays on with no question, and Generate refuses it as before`, async ({ browser }) => {
    test.setTimeout(420000);
    const page = await openReachFixture(browser, mode);
    try {
      await addDepthSeries(page, rev7Depth(mode, 1));
      await setRow(page, mode, 'depth', { id: 'cov' });
      await removeDepthFile(page, mode);
      await expectNoQuestion(page);
      expect((await stackState(page, mode)).rows).toContain('cov:true');
      expect(await generateError(page)).toBe('TRACK_INVALID');
    } finally { await page.context().close(); }
  });

  test(`rev 7: ${mode} the last file of one of two series says the series stays empty, and Generate asks for a file or its removal`, async ({ browser }) => {
    test.setTimeout(420000);
    const page = await openReachFixture(browser, mode);
    try {
      await addDepthSeries(page, rev7Depth(mode, 2));
      const seriesRow = await page.evaluate(async (stackMode) => {
        const { state } = await import('./js/state.js');
        return state.activeDrawing().adv[`${stackMode}_track_slots`].find((row) => row.renderer === 'depth' && row.params.track_index === 0).id;
      }, mode);
      await addReachRow(page, { id: 'mark', anchor: seriesRow, setId: setIdOf(mode), mode });
      await expectGenerates(page);
      await settleLive(page);
      const start = await stackState(page, mode);
      const title = `Remove the Depth file and turn off Track row ${rowNumber(start, 'mark')}?`;
      await removeDepthFile(page, mode);
      await expect(page.locator('[data-choice-prompt]').getByRole('heading', { name: title, exact: true })).toBeVisible();
      expect(await promptText(page)).toContain('Depth series #1 stays empty, and Generate will ask for a file or its removal.');
      await page.locator('[data-choice-prompt]').getByRole('button', { name: 'Cancel', exact: true }).click();
      await settleLive(page);
      await exerciseAsk(page, mode, {
        trigger: () => removeDepthFile(page, mode),
        title,
        generates: false,
        expectAfter: (after) => {
          expect(after.rows).toContain('mark:false');
          expect(after.series, 'the series stays').toBe(2);
          expect(after.depth).not.toContain('s1.tsv');
          expect(after.depth).toContain('s2.tsv');
          expect(after.showDepth).toBe(true);
        }
      });
      expect(await generateError(page)).toBe('DEPTH_INVALID');
    } finally { await page.context().close(); }
  });

  test(`rev 7: ${mode} the sentence is left out when it is not true: a series removal, and the only series' last file`, async ({ browser }) => {
    test.setTimeout(420000);
    const page = await openReachFixture(browser, mode);
    try {
      await addDepthSeries(page, rev7Depth(mode, 2));
      const seriesRow = await page.evaluate(async (stackMode) => {
        const { state } = await import('./js/state.js');
        return state.activeDrawing().adv[`${stackMode}_track_slots`].find((row) => row.renderer === 'depth' && row.params.track_index === 0).id;
      }, mode);
      await addReachRow(page, { id: 'mark', anchor: seriesRow, setId: setIdOf(mode), mode });
      const section = await openDepthSources(page, mode);
      await section.getByRole('button', { name: mode === 'circular' ? 'Remove' : 'Remove series', exact: true }).first().click();
      await expect(page.locator('[data-choice-prompt]')).toBeVisible();
      expect(await promptText(page)).not.toContain('stays empty');
      await page.locator('[data-choice-prompt]').getByRole('button', { name: 'Cancel', exact: true }).click();
      await removeDepthFile(page, mode);
      await expect(page.locator('[data-choice-prompt]')).toBeVisible();
      expect(await promptText(page), 'a second series still has a file').toContain('stays empty');
      await chooseAsked(page);
      // The second series' last file, no row reaching it: nothing asks.
      await removeDepthFile(page, mode);
      await expectNoQuestion(page);
    } finally { await page.context().close(); }
  });
}

// Generate asks for a file only while it draws Depth: with the stack off, Show Depth.
test('rev 7: Linear with the stack and Show Depth off, the last file of one of two series asks without the empty-series sentence, and Generate draws', async ({ browser }) => {
  test.setTimeout(420000);
  const page = await openReachFixture(browser, 'linear');
  try {
    await addDepthSeries(page, rev7Depth('linear', 2));
    const seriesRow = await page.evaluate(async () => {
      const { state } = await import('./js/state.js');
      return state.activeDrawing().adv.linear_track_slots.find((row) => row.renderer === 'depth' && row.params.track_index === 0).id;
    });
    await addReachRow(page, { id: 'mark', anchor: seriesRow, setId: setIdOf('linear'), mode: 'linear' });
    await page.evaluate(async () => {
      window.__GBDRAW_APP__.setLinearTrackSlotsEnabled(false);
      (await import('./js/state.js')).state.activeDrawing().form.show_depth = false;
    });
    await settleLive(page);
    await expectGenerates(page);
    const start = await stackState(page, 'linear');
    await removeDepthFile(page, 'linear');
    await expect(page.locator('[data-choice-prompt]').getByRole('heading', { name: `Remove the Depth file and turn off Track row ${rowNumber(start, 'mark')}?`, exact: true })).toBeVisible();
    expect(await promptText(page)).not.toContain('stays empty');
    await chooseAsked(page);
    expect((await stackState(page, 'linear')).depth).toContain('s2.tsv');
    await expectGenerates(page);
  } finally { await page.context().close(); }
});

test('rev 7: Linear a moved Depth row leaves with its last file, and a re-added file puts it at the default place', async ({ browser }) => {
  test.setTimeout(420000);
  const page = await generatedLinear(browser);
  try {
    await page.evaluate(async () => {
      const { state } = await import('./js/state.js');
      const rows = state.activeDrawing().adv.linear_track_slots;
      window.__GBDRAW_APP__.moveLinearTrackSlot(rows.findIndex((row) => row.id === 'depth'), rows.findIndex((row) => row.id === 'gc_content'));
    });
    await settleLive(page);
    const order = async () => (await stackState(page, 'linear')).rows.map((row) => row.split(':')[0]);
    const moved = await order();
    await removeDepthFile(page, 'linear');
    await expectNoQuestion(page);
    expect(await order()).not.toContain('depth');
    await addDepthFile(page, 'linear');
    await settleLive(page);
    const back = await order();
    expect(back.at(-1), 'the default place').toBe('depth');
    expect(moved.at(-1)).not.toBe('depth');
    await expectGenerates(page);
  } finally { await page.context().close(); }
});

// A Linear File removal or clear holding the last file of a series (OV-2284).
const TESTA_GB = readFileSync('tests/fixtures/b_collide.gb', 'utf8');
const testaDepth = (series) => [0, 1, 2].map((index) => `TESTA\t${index * 900 + 1}\t${10 + index + series}`).join('\n');
const twoFiles = async (browser, { holding = [0], other = [], gb = TESTA_GB } = {}) => {
  const page = await openReachFixture(browser, 'linear');
  await evaluateWithRetainedPromise(page, async ({ text, on, off, tsvOf, lambdaTsv }) => {
    const app = window.__GBDRAW_APP__;
    app.addLinearSeq();
    app.setLinearSeqPrimaryFile(1, 'gb', new File([text], 'testa.gb', { type: 'text/plain', lastModified: 1 }));
    await app.refreshLinearRecordSelectors();
    await window.__GBDRAW_HISTORY__.runUndoable('Setup', () => {
      on.forEach((series) => app.setLinearDepthFile(app.linearSeqs[1], series, new File([tsvOf[series]], `s${series + 1}.tsv`)));
      off.forEach((series) => app.setLinearDepthFile(app.linearSeqs[0], series, new File([lambdaTsv], `s${series + 1}.tsv`)));
    });
  }, { text: gb, on: holding, off: other, tsvOf: [testaDepth(0), testaDepth(1)], lambdaTsv: rev7Depth('linear', 1).tsv });
  const ids = await page.evaluate(async () => (await import('./js/state.js')).state.activeDrawing().adv.linear_track_slots
    .filter((row) => row.renderer === 'depth').map((row) => [row.params.track_index, row.id]));
  for (const [series, id] of ids) {
    if (holding.includes(series)) await addReachRow(page, { id: `mark${series}`, anchor: id, setId: 'marks', mode: 'linear' });
  }
  return page;
};
const secondCard = (page) => page.locator('[data-linear-source-card]').nth(1);
const fileCardRemove = async (page, card = secondCard(page), action = 'Remove file') => {
  await reveal(card);
  await card.evaluate((node) => node.querySelectorAll('details').forEach((fold) => { fold.open = true; }));
  await card.getByRole('button', { name: /^Remove GenBank/ }).first().click();
  const dialog = page.locator('[data-linear-source-removal-dialog]');
  await expect(dialog).toBeVisible();
  await dialog.getByRole('button', { name: action, exact: true }).click();
  await expect(dialog).toHaveCount(0);
};
const LINEAR_FILE_ROUTES = {
  'File card Remove': (page) => fileCardRemove(page),
  'the record list Remove File': async (page) => {
    const list = page.locator('[data-record-list-dialog]');
    if (!(await list.count())) await secondCard(page).locator('[data-record-list-open]').click();
    await expect(list).toBeVisible();
    await list.getByRole('button', { name: 'Select none', exact: true }).click();
    const question = page.locator('[data-record-remove-file-dialog]');
    await expect(question).toBeVisible();
    await question.locator('[data-record-remove-file]').click();
    await expect(question).toHaveCount(0);
  },
  'the primary file set to none': (page) => page.evaluate(() => { window.__GBDRAW_APP__.setLinearSeqPrimaryFile(1, 'gb', null); })
};
for (const [route, trigger] of Object.entries(LINEAR_FILE_ROUTES)) {
  test(`rev 7: Linear ${route} holding a series' last file asks after its own question, as one Undo step`, async ({ browser }) => {
    test.setTimeout(420000);
    // The record list opens for a File with more than one record.
    const page = await twoFiles(browser, { gb: `${TESTA_GB}\n${TESTA_GB.replaceAll('TESTA', 'TESTB')}` });
    try {
      await expectGenerates(page);
      await settleLive(page);
      const start = await stackState(page, 'linear');
      await exerciseAsk(page, 'linear', {
        trigger: () => trigger(page),
        title: `Remove the File and turn off Track row ${rowNumber(start, 'mark0')}?`,
        beforeGenerate: async () => {
          const list = page.locator('[data-record-list-dialog]');
          if (await list.count()) await list.getByRole('button', { name: 'Close', exact: true }).click();
        },
        expectAfter: (after) => {
          expect(after.rows).toContain('mark0:false');
          expect(depthRowsOf(after), 'the row leaves').toEqual([]);
          expect(after.depth).not.toContain('s1.tsv');
          expect(after.showDepth).toBe(false);
        }
      });
    } finally { await page.context().close(); }
  });
}

test('rev 7: Linear File removal: Cancel on the second question keeps the File, returns the focus, and other pending values stay', async ({ browser }, testInfo) => {
  test.setTimeout(420000);
  const page = await twoFiles(browser);
  try {
    const pending = () => page.evaluate(async () => {
      const { state } = await import('./js/state.js');
      const { form, adv } = state.activeDrawing();
      const { show_depth: _showDepth, ...otherForm } = form;
      const { depth_tracks: _series, linear_track_slots: _slots, linear_track_slots_axis_index: _axis, ...otherAdv } = adv;
      return JSON.stringify({ otherForm, otherAdv, untouched: adv.linear_track_slots.filter((row) => row.id !== 'depth' && row.id !== 'mark0') });
    });
    await page.evaluate(async () => { (await import('./js/state.js')).state.activeDrawing().adv.label_font_size = 29; });
    const before = await pending();
    const start = await stackState(page, 'linear');
    const files = () => page.evaluate(() => window.__GBDRAW_APP__.linearSeqs.map((seq) => seq.gb?.name ?? null));
    const filesBefore = await files();
    const steps = await undoCount(page);
    await fileCardRemove(page);
    await expect(page.locator('[data-choice-prompt]').getByRole('button', { name: 'Cancel', exact: true }), 'the second question takes the focus').toBeFocused();
    await page.locator('[data-choice-prompt]').getByRole('button', { name: 'Cancel', exact: true }).click();
    await expect(page.locator('[data-choice-prompt]')).toHaveCount(0);
    await settleLive(page);
    expect(await files()).toEqual(filesBefore);
    expect(await stackState(page, 'linear')).toEqual(start);
    expect(await undoCount(page)).toBe(steps);
    await expect(secondCard(page).getByRole('button', { name: /^Remove GenBank/ }).first()).toBeFocused();
    await fileCardRemove(page);
    await chooseAsked(page);
    expect(await files()).toEqual([filesBefore[0]]);
    expect(await pending()).toBe(before);
    await historyStep(page, 'undo');
    expect(await files()).toEqual(filesBefore);
    expect(await stackState(page, 'linear')).toEqual(start);
    const saved = await savedSession(page, testInfo.outputPath('undo.gbdraw-session.json'));
    expect(saved.modes.linear.config.adv.linear_track_slots.find((row) => row.id === 'mark0').enabled).toBe(true);
    await expectGenerates(page);
  } finally { await page.context().close(); }
});

test('rev 7: Linear File removal losing two series asks once, with the union of the rows', async ({ browser }) => {
  test.setTimeout(420000);
  const page = await twoFiles(browser, { holding: [0, 1] });
  try {
    const start = await stackState(page, 'linear');
    await exerciseAsk(page, 'linear', {
      trigger: () => fileCardRemove(page),
      title: `Remove the File and turn off Track rows ${rowNumber(start, 'mark0')} and ${rowNumber(start, 'mark1')}?`,
      expectAfter: (after) => {
        expect(after.rows).toEqual(expect.arrayContaining(['mark0:false', 'mark1:false']));
        expect(depthRowsOf(after)).toEqual([]);
      }
    });
  } finally { await page.context().close(); }
});

test('rev 7: Linear File removal with another series left says the series stays empty', async ({ browser }) => {
  test.setTimeout(420000);
  const page = await twoFiles(browser, { holding: [0], other: [1] });
  try {
    const start = await stackState(page, 'linear');
    const title = `Remove the File and turn off Track row ${rowNumber(start, 'mark0')}?`;
    await fileCardRemove(page);
    await expect(page.locator('[data-choice-prompt]').getByRole('heading', { name: title, exact: true })).toBeVisible();
    expect(await promptText(page)).toContain('Depth series #1 stays empty, and Generate will ask for a file or its removal.');
    await chooseAsked(page);
    expect(await generateError(page)).toBe('DEPTH_INVALID');
  } finally { await page.context().close(); }
});

test('rev 7: Linear Empty this slot of the only File holding a series last file asks, as one Undo step', async ({ browser }) => {
  test.setTimeout(420000);
  const page = await openReachFixture(browser, 'linear');
  try {
    await addDepthSeries(page, rev7Depth('linear', 1));
    await addReachRow(page, { id: 'mark', anchor: 'depth', setId: 'marks', mode: 'linear' });
    const start = await stackState(page, 'linear');
    await exerciseAsk(page, 'linear', {
      trigger: () => fileCardRemove(page, page.locator('[data-linear-source-card]').first(), 'Empty this slot'),
      title: `Remove the File and turn off Track row ${rowNumber(start, 'mark')}?`,
      generates: false,
      expectAfter: (after) => {
        expect(after.rows).toContain('mark:false');
        expect(depthRowsOf(after)).toEqual([]);
      }
    });
  } finally { await page.context().close(); }
});

test('rev 7: Linear File removal with no row reaching what it takes asks only its own question', async ({ browser }) => {
  test.setTimeout(420000);
  const page = await twoFiles(browser);
  try {
    await setRow(page, 'linear', 'mark0', { enabled: false });
    await fileCardRemove(page);
    await expectNoQuestion(page);
    expect(depthRowsOf(await stackState(page, 'linear'))).toEqual([]);
    await expectGenerates(page);
  } finally { await page.context().close(); }
});

// A Circular managed comparison ring row that leaves with its source (OV-2289).
const loadRingSession = async (browser, { stack = true, session = ringSessionUpload() } = {}) => {
  const context = await browser.newContext({ viewport: { width: 1600, height: 1000 } });
  const page = await context.newPage();
  await page.addInitScript(() => { window.__LIFECYCLE__ = []; window.__GBDRAW_TEST_HOOKS__ = { onSessionLifecycleEvent: (event) => window.__LIFECYCLE__.push(event) }; });
  page.on('dialog', (dialog) => dialog.accept());
  await openApp(page);
  await page.locator('input[accept^=".json,"]').setInputFiles(session);
  await page.waitForFunction(() => window.__LIFECYCLE__.some((event) => event.name === 'interactiveReady') && !window.__GBDRAW_APP__.sessionImportPending);
  await page.locator('[data-setting="region-annotations"] input[type=file]').first().setInputFiles({
    name: 'regions.tsv', mimeType: 'text/tab-separated-values', buffer: Buffer.from('set_id\tid\tmark\trecord\tstart\tend\nmarks\tm1\thighlight\treference-a\t100\t900\n')
  });
  await expect.poll(() => page.evaluate(async () => (await import('./js/state.js')).state.activeDrawing().annotationSets.length)).toBeGreaterThan(0);
  await page.evaluate(() => window.__GBDRAW_APP__.setCircularTrackSlotsEnabled(true));
  await addReachRow(page, { id: 'mark', anchor: 'conservation_comparison_b', setId: 'marks', mode: 'circular' });
  if (!stack) await page.evaluate(() => window.__GBDRAW_APP__.setCircularTrackSlotsEnabled(false));
  return page;
};
const removeRing = async (page) => {
  const button = page.locator('button[title="Remove series"]').first();
  await reveal(button);
  await button.click();
};

test('rev 7: Circular removing a comparison source whose managed ring row another row reaches asks, as one Undo step, and Generate draws', async ({ browser }, testInfo) => {
  test.setTimeout(420000);
  const page = await loadRingSession(browser);
  try {
    const rings = () => page.evaluate(() => window.__GBDRAW_APP__.circularConservation.series.length);
    const undoSave = async () => {
      expect(await rings()).toBe(2);
      await historyStep(page, 'undo');
      expect(await rings()).toBe(3);
      const saved = await savedSession(page, testInfo.outputPath('ring-undo.gbdraw-session.json'));
      expect(saved.modes.circular.config.adv.circular_track_slots.find((row) => row.id === 'mark').enabled).not.toBe(false);
      await historyStep(page, 'redo');
    };
    await expectGenerates(page);
    await settleLive(page);
    const start = await stackState(page, 'circular');
    await exerciseAsk(page, 'circular', {
      trigger: () => removeRing(page),
      title: `Remove the comparison source and turn off Track row ${rowNumber(start, 'mark')}?`,
      beforeGenerate: undoSave,
      expectAfter: (after) => {
        expect(after.rows).toContain('mark:false');
        expect(after.rows.some((row) => row.startsWith('conservation_comparison_b:'))).toBe(false);
      }
    });
  } finally { await page.context().close(); }
});

// OV-2442: a search gives each ring another source file (other bytes, another key). The ring's row stays: its id,
// by which another row anchors to it and the next Generate draws it, and its edits. The Session's rows have no
// runtime record, so the first Generate searches; the stub gives other rows than the stored ones.
test('a ring search that writes other rows keeps the managed ring row, its edits, and the row that reaches it', async ({ browser }) => {
  test.setTimeout(420000);
  const page = await loadRingSession(browser, { session: RING_SESSION_FILE });
  try {
    const rows = () => page.evaluate(async () => {
      const { state } = await import('./js/state.js');
      return state.activeDrawing().adv.circular_track_slots.flatMap((row) => (
        row.renderer === 'sequence_conservation' ? [[row.id, row.width]]
          : row.renderer === 'annotations' ? [[row.id, row.params.anchor_slot]] : []
      ));
    });
    await page.evaluate(async () => {
      const { state } = await import('./js/state.js');
      state.losatExecution.executionMode = 'serial';
      window.__GBDRAW_LOSAT_EXECUTOR__ = async (jobs) => {
        window.__RING_SEARCHES__ = (window.__RING_SEARCHES__ || 0) + jobs.length;
        return jobs.map((job) => ({
          cacheKey: job.cacheKey,
          text: `${job.queryName || 'query'}\treference-a\t99\t300\t3\t0\t1\t300\t1\t300\t1e-20\t500\n`
        }));
      };
      const ring = state.activeDrawing().adv.circular_track_slots.find((row) => row.id === 'conservation_comparison_b');
      ring.width = '9px';
    });
    const before = await rows();
    expect(before).toEqual(expect.arrayContaining([['conservation_comparison_b', '9px'], ['mark', 'conservation_comparison_b']]));
    await expectGenerates(page);
    expect(await page.evaluate(() => window.__RING_SEARCHES__)).toBe(3);
    expect(await rows()).toEqual(before);
    const drawn = () => page.evaluate(() => {
      const app = window.__GBDRAW_APP__;
      return String(app.currentResults()[app.selectedResultIndex].content);
    });
    const first = await drawn();
    await expectGenerates(page);
    expect(await page.evaluate(() => window.__RING_SEARCHES__), 'the second Generate reuses the search').toBe(3);
    expect(await rows()).toEqual(before);
    expect(await drawn()).toBe(first);
  } finally { await page.context().close(); }
});

const addRing = (page, name) => evaluateWithRetainedPromise(page, async ({ name }) => {
  const file = new File([`>${name}\nACGTACGTACGTACGTACGTACGT\n`], `${name}.fasta`, { type: 'text/plain', lastModified: 5 });
  await window.__GBDRAW_APP__.addCircularConservationComparisonFile({ target: { files: [file], value: '' } });
}, { name });
// HmmtDNA, whose tRNA features give Rotate a target, with three rings from Add Seq on the custom stack.
const loadAddedRings = async (browser) => {
  const page = await load(browser, HMMT);
  await page.evaluate(async () => {
    const { state } = await import('./js/state.js');
    state.losatExecution.executionMode = 'serial';
    window.__GBDRAW_LOSAT_EXECUTOR__ = async (jobs) => jobs.map((job, index) => ({
      cacheKey: job.cacheKey,
      text: `${job.queryName || 'query'}\tNC_012920.1\t99\t${2100 + index}\t3\t0\t1\t${2100 + index}\t1\t${2100 + index}\t1e-20\t500\n`
    }));
  });
  for (const name of ['ringa', 'ringb', 'ringc']) await addRing(page, name);
  await page.locator('[data-setting="region-annotations"] input[type=file]').first().setInputFiles({
    name: 'regions.tsv', mimeType: 'text/tab-separated-values', buffer: Buffer.from('set_id\tid\tmark\trecord\tstart\tend\nmarks\tm1\thighlight\tNC_012920.1\t100\t900\n')
  });
  await expect.poll(() => page.evaluate(async () => (await import('./js/state.js')).state.activeDrawing().annotationSets.length)).toBeGreaterThan(0);
  await page.evaluate(() => window.__GBDRAW_APP__.setCircularTrackSlotsEnabled(true));
  return page;
};
const ringRows = (page) => page.evaluate(async () => {
  const { state } = await import('./js/state.js');
  return state.activeDrawing().adv.circular_track_slots.flatMap((row) => (
    row.renderer === 'sequence_conservation' ? [{ id: row.id, width: row.width, file: String(row.params.series_key).split('.')[0], label: row.params.label }]
      : []
  ));
});
const setRingWidth = (page, id, width) => page.evaluate(async ({ id, width }) => {
  const { state } = await import('./js/state.js');
  state.activeDrawing().adv.circular_track_slots.find((row) => row.id === id).width = width;
}, { id, width });
const rotateByTrna = async (page) => {
  await page.evaluate(() => { window.__GBDRAW_APP__.richFeaturePopup = true; });
  const search = page.getByRole('searchbox', { name: 'Search features', exact: true });
  await search.fill('tRNA');
  await search.press('Enter');
  await page.getByRole('button', { name: 'Open active feature', exact: true }).click();
  await page.getByRole('button', { name: 'Rotate record using this feature', exact: true }).click();
  await page.getByRole('region', { name: 'Rotate record using this feature', exact: true }).getByRole('button', { name: 'Apply now' }).click();
  await page.waitForFunction(() => !window.__GBDRAW_APP__.processing, null, { timeout: 120000 });
  await settleLive(page);
  expect(await page.evaluate(async () => (await import('./js/state.js')).state.errorLog.value?.code ?? null)).toBeNull();
};
// Each drawn ring by its row, with its colors and radii.
const drawnRings = (page) => page.evaluate(() => {
  const app = window.__GBDRAW_APP__;
  const doc = new DOMParser().parseFromString(String(app.currentResults()[app.selectedResultIndex].content), 'image/svg+xml');
  return [...doc.querySelectorAll('g[id*="_conservation_"]')].map((group) => {
    const paths = [...group.querySelectorAll('path')];
    return {
      row: group.id.replace(/^track_slot_[0-9a-f]+_/, ''),
      fills: [...new Set(paths.map((path) => path.getAttribute('fill')))],
      radii: [...new Set(paths.flatMap((path) => [...String(path.getAttribute('d')).matchAll(/A\s*([\d.]+),/g)].map((match) => match[1])))]
    };
  });
});

test('a Rotate after a ring reorder leaves each ring its own row, width, and anchor', async ({ browser }) => {
  test.setTimeout(420000);
  const page = await loadAddedRings(browser);
  try {
    await expectGenerates(page);
    await setRingWidth(page, 'conservation_ringa', '30px');
    await addReachRow(page, { id: 'mark', anchor: 'conservation_ringa', setId: 'marks', mode: 'circular' });
    await expectGenerates(page);
    await settleLive(page);
    const before = await ringRows(page);
    expect(before.map((row) => [row.id, row.file, row.width])).toEqual([
      ['conservation_ringa', 'ringa', '30px'], ['conservation_ringb', 'ringb', null], ['conservation_ringc', 'ringc', null]
    ]);
    const drawn = await drawnRings(page);
    expect(drawn.map((ring) => ring.row)).toEqual(before.map((row) => row.id));
    await page.evaluate(() => window.__GBDRAW_APP__.moveCircularConservationSeries(0, 1));
    await rotateByTrna(page);
    expect(await ringRows(page), 'the Rotate adopts the committed rings, and each row stays on its own').toEqual(before);
    await expectGenerates(page);
    await settleLive(page);
    const after = await ringRows(page);
    expect(after.map((row) => [row.id, row.label, row.width])).toEqual(before.map((row) => [row.id, row.label, row.width]));
    expect(await drawnRings(page), 'the next Generate draws each row in the color and at the radii it had before the Rotate').toEqual(drawn);
    expect(await page.evaluate(async () => (await import('./js/state.js')).state.activeDrawing().adv.circular_track_slots.find((row) => row.id === 'mark').params.anchor_slot))
      .toBe('conservation_ringa');
  } finally { await page.context().close(); }
});

test('a Rotate after a Remove and an Add Seq puts no edit on another ring', async ({ browser }) => {
  test.setTimeout(420000);
  const page = await loadAddedRings(browser);
  try {
    await expectGenerates(page);
    await settleLive(page);
    await page.evaluate(() => window.__GBDRAW_APP__.removeCircularConservationSource(0));
    await addRing(page, 'ringz');
    const added = (await ringRows(page)).find((row) => row.file === 'ringz');
    await setRingWidth(page, added.id, '7px');
    await rotateByTrna(page);
    expect((await ringRows(page)).map((row) => [row.file, row.label, row.width]), 'the edit of ringz is not on a committed ring').toEqual([
      ['ringa', 'ringa', null], ['ringb', 'ringb', null], ['ringc', 'ringc', null]
    ]);
  } finally { await page.context().close(); }
});

// While a Generate runs, every control that rewrites the ring series is disabled (Run LOSAT and Upload BLAST TSV,
// Move series, Remove series, Add Seq) and Reset Settings is refused, so a keyboard user behind the overlay cannot
// change the rings either. A row moves to its ring's new key only if the rings at adoption are the rings the request
// was built from; that compare stays for a change that does not come through these controls.
const revealRingControls = (page) => reveal(page.locator('button[aria-label="Add Seq"], button:has-text("Add Seq")').first());
const expectEachRowOnItsRing = async (page) => {
  const rows = await ringRows(page);
  expect(rows.map((row) => row.file).sort()).toEqual(['ringa', 'ringb', 'ringc']);
  expect(rows.map((row) => [row.id.replace(/^conservation_|_\d+$/g, ''), row.label]), 'each row\'s id and name are of the file it draws')
    .toEqual(rows.map((row) => [row.file, row.file]));
};
const ringLabels = (page) => page.evaluate(() => window.__GBDRAW_APP__.circularConservation.series.map((entry) => entry.label));
// Each control that changes the rings, as a keyboard user reaches it: Move series both ways, Remove series, Add Seq,
// Reset Settings, and the two source radios (last: Upload BLAST TSV would hide the others).
const ringControls = (page) => [
  page.getByRole('button', { name: 'Move series inward', exact: true }).first(),
  page.getByRole('button', { name: 'Move series outward', exact: true }).nth(1),
  page.getByRole('button', { name: 'Remove series', exact: true }).last(),
  page.getByRole('button', { name: 'Add Seq', exact: true }),
  page.getByRole('button', { name: 'Reset Settings', exact: true }),
  page.locator('[data-setting="pairwise-comparisons"]').getByRole('radio', { name: 'Upload BLAST TSV', exact: true }),
  page.locator('[data-setting="pairwise-comparisons"]').getByRole('radio', { name: 'Run LOSAT', exact: true })
];
// Focus the control by keyboard focus, not by pointer, and press the key; a disabled control takes no focus.
const pressOn = async (page, control, key) => {
  await page.evaluate(() => document.activeElement?.blur());
  await control.evaluate((element) => element.focus());
  await page.keyboard.press(key);
};
const holdNextSearch = (page, { fail = false } = {}) => page.evaluate((shouldFail) => {
  window.__RELEASE_SEARCH__ = null;
  window.__HOLD_SEARCH__ = true;
  window.__FAIL_SEARCH__ = shouldFail;
}, fail);
const releaseSearch = async (page) => {
  await page.evaluate(() => window.__RELEASE_SEARCH__());
  await page.waitForFunction(() => !window.__GBDRAW_APP__.processing, null, { timeout: 120000 });
  await settleLive(page);
};

test('the ring controls are disabled while a Generate runs, and work again after it ends or fails', async ({ browser }) => {
  test.setTimeout(420000);
  const page = await loadAddedRings(browser);
  try {
    await revealRingControls(page);
    await page.evaluate(() => {
      const search = window.__GBDRAW_LOSAT_EXECUTOR__;
      window.__GBDRAW_LOSAT_EXECUTOR__ = async (jobs) => {
        if (window.__HOLD_SEARCH__) {
          window.__HOLD_SEARCH__ = false;
          await new Promise((resolve) => { window.__RELEASE_SEARCH__ = resolve; });
          if (window.__FAIL_SEARCH__) throw new Error('stub search failure');
        }
        return search(jobs);
      };
    });
    let chosenFiles = 0;
    page.on('filechooser', async (chooser) => {
      chosenFiles += 1;
      await chooser.setFiles({ name: 'ringz.fasta', mimeType: 'text/plain', buffer: Buffer.from('>ringz\nACGTACGTACGTACGTACGTACGT\n') });
    });
    const generateButton = page.getByRole('button', { name: 'Generate Diagram', exact: true });
    // Out of file order, so a control that derives the series again from the files changes it.
    await ringControls(page)[0].click();
    expect(await ringLabels(page)).toEqual(['ringb', 'ringa', 'ringc']);
    for (const fail of [true, false]) {
      await holdNextSearch(page, { fail });
      await generateButton.click();
      await page.waitForFunction(() => Boolean(window.__RELEASE_SEARCH__), null, { timeout: 120000 });
      const [inward, outward, remove, addSeq, reset, upload, losat] = ringControls(page);
      await pressOn(page, inward, 'Enter');
      await pressOn(page, outward, ' ');
      await pressOn(page, remove, 'Enter');
      await pressOn(page, addSeq, ' ');
      await pressOn(page, reset, 'Enter');
      await expect(page.getByRole('dialog', { name: 'Reset all settings?', exact: true }), 'Reset Settings asks nothing while the search runs')
        .toHaveCount(0);
      expect(await page.evaluate(() => window.__GBDRAW_APP__.resetSettings()), 'the owner refuses Reset Settings while the search runs')
        .toEqual({ status: 'busy', reason: 'Generating diagram. Retry after generation finishes.' });
      await pressOn(page, upload, ' ');
      await pressOn(page, losat, ' ');
      expect(await ringLabels(page), 'no key press changes the rings while the search runs').toEqual(['ringb', 'ringa', 'ringc']);
      expect(chosenFiles, 'Add Seq opens no file chooser while the search runs').toBe(0);
      for (const control of ringControls(page)) await expect(control).toBeDisabled();
      await releaseSearch(page);
      expect(await page.evaluate(async () => (await import('./js/state.js')).state.errorLog.value === null)).toBe(!fail);
      for (const control of ringControls(page)) await expect(control).toBeEnabled();
    }
    await expectEachRowOnItsRing(page);
    await pressOn(page, ringControls(page)[0], 'Enter');
    expect(await ringLabels(page), 'Move series works again by keyboard after the Generate').toEqual(['ringa', 'ringb', 'ringc']);
  } finally { await page.context().close(); }
});

// A change through the hook, not a control: the compare at adoption replaces the rows instead of moving them to the
// other ring's key. This checks the state right after that Generate only; the series it adopts still takes the draft's
// entry at the same place (`adoptCanonicalRenderArtifacts`). No control in the Web UI rewrites the series while a
// Generate runs: the ring controls and the source radios are disabled, and Reset Settings is refused.
test('a ring reorder after the request is built puts no row on another ring', async ({ browser }) => {
  test.setTimeout(420000);
  const page = await loadAddedRings(browser);
  try {
    await revealRingControls(page);
    await page.evaluate(() => {
      const hooks = window.__GBDRAW_TEST_HOOKS__;
      const previous = hooks.onSessionLifecycleEvent;
      hooks.onSessionLifecycleEvent = (event) => {
        previous?.(event);
        if (event.name === 'canonical-request-built' && !window.__RING_MOVED__) {
          window.__GBDRAW_APP__.moveCircularConservationSeries(0, 1);
          window.__RING_MOVED__ = window.__GBDRAW_APP__.circularConservation.series.map((entry) => entry.label);
        }
      };
    });
    await expectGenerates(page);
    expect(await page.evaluate(() => window.__RING_MOVED__)).toEqual(['ringb', 'ringa', 'ringc']);
    await settleLive(page);
    await expectEachRowOnItsRing(page);
  } finally { await page.context().close(); }
});

test('rev 7: Circular removing a comparison source no row reaches asks nothing', async ({ browser }) => {
  test.setTimeout(420000);
  const page = await loadRingSession(browser);
  try {
    const remove = page.locator('button[title="Remove series"]').nth(1);
    await reveal(remove);
    await remove.click();
    await expectNoQuestion(page);
    expect((await stackState(page, 'circular')).rows.some((row) => row.startsWith('conservation_comparison_c:'))).toBe(false);
    await expectGenerates(page);
  } finally { await page.context().close(); }
});

test('rev 7: Circular clearing the uploaded comparison files asks, and Cancel keeps them', async ({ browser }) => {
  test.setTimeout(420000);
  const page = await load(browser, TOBACCO);
  try {
    const ringId = await evaluateWithRetainedPromise(page, async () => {
      const { state } = await import('./js/state.js');
      const app = window.__GBDRAW_APP__;
      const drawing = state.activeDrawing();
      drawing.circularConservation.source = 'upload';
      app.setCircularConservationUploadFiles([new File(['NC_001879.2\tq\t90\t100\t0\t0\t1\t100\t1\t100\t1e-5\t100\n'], 'ring.tsv')]);
      return drawing.adv.circular_track_slots.find((row) => row.renderer === 'sequence_conservation').id;
    });
    await addReachRow(page, { id: 'mark', anchor: ringId, setId: 'plastome_tint', mode: 'circular' });
    const start = await stackState(page, 'circular');
    const uploader = page.locator('button[aria-label="Remove BLAST outfmt 6/7 files"]');
    await reveal(uploader);
    await exerciseAsk(page, 'circular', {
      trigger: () => uploader.click(),
      title: `Remove the comparison source and turn off Track row ${rowNumber(start, 'mark')}?`,
      generates: false,
      expectAfter: (after) => {
        expect(after.rows).toContain('mark:false');
        expect(after.rows.some((row) => row.startsWith(`${ringId}:`))).toBe(false);
      }
    });
  } finally { await page.context().close(); }
});
