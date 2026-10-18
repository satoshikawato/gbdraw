const { test, expect } = require('@playwright/test');
const { openApp, reveal } = require('./helpers/app-lifecycle.cjs');
const { load, generate, switchMode } = require('./helpers/mode-transition.cjs');

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
    await expect(dialog(page)).toContainText('Reset to Tuckin removes the lane for 1 feature you placed manually.');
    await expect(resetButton(page)).toBeFocused();
    expect(await layoutState(page, key)).toEqual(placed);
    await page.screenshot({ path: testInfo.outputPath('preset-reset-dialog.png') });
    await cancelButton(page).click();
    await expect(dialog(page)).toBeHidden();
    await expect(tuckin).toBeFocused();
    expect(await layoutState(page, key)).toEqual(placed);

    // Row reorder: Reset applies the move and removes the row as one step.
    await runRowAction(stackRow(page, 'circular'), 'Move outside Axis');
    await expect(dialog(page)).toContainText('Move outside Axis removes the lane for 1 feature you placed manually.');
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
    await expect(dialog(page)).toContainText('Changing Use custom stack to On removes the lane for 1 feature you placed manually.');
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
