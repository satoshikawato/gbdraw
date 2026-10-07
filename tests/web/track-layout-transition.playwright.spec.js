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
  return JSON.stringify([feature.scope, feature.record_key, feature.biological_feature_id]);
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
const cancelButton = (page) => dialog(page).getByRole('button', { name: 'Cancel change', exact: true });
const openStack = async (page, mode) => {
  const toggle = await reveal(page.locator(`[aria-controls="${mode}-custom-track-slots-panel"]`));
  if (await page.locator(`#${mode}-custom-track-slots-panel`).count() === 0) await toggle.click();
  await expect(page.locator(`#${mode}-custom-track-slots-panel`)).toBeVisible();
};
const stackRow = (page, mode) => page.locator('.track-slot-row-head', {
  has: page.locator(`[aria-label="Enable ${mode} track slot features"]`)
});
const undo = (page) => page.getByRole('button', { name: /^Undo/ }).first().click();

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
    await expect(dialog(page)).toContainText('Reset to Tuckin leaves 1 Feature placement without its lane.');
    await expect(resetButton(page)).toBeFocused();
    expect(await layoutState(page, key)).toEqual(placed);
    await page.screenshot({ path: testInfo.outputPath('preset-reset-dialog.png') });
    await cancelButton(page).click();
    await expect(dialog(page)).toBeHidden();
    await expect(tuckin).toBeFocused();
    expect(await layoutState(page, key)).toEqual(placed);

    // Row reorder: Reset applies the move and removes the row as one step.
    await stackRow(page, 'circular').locator('button[title="Move outside Axis"]').click();
    await expect(dialog(page)).toContainText('Move outside Axis leaves 1 Feature placement without its lane.');
    expect(await layoutState(page, key)).toEqual(placed);
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

test('Circular Separate Strands and Linear Use custom stack ask before dropping a Linear lane (R-3)', async ({ browser }, testInfo) => {
  test.setTimeout(300000);
  const page = await load(browser, LAMBDA);
  try {
    // A saved custom stack whose feature row sits above the Axis has no lanes.
    await openStack(page, 'linear');
    const custom = page.getByRole('checkbox', { name: 'Use custom stack', exact: true });
    await custom.check();
    await stackRow(page, 'linear').locator('button[title="Move above Axis"]').click();
    await custom.uncheck();
    const strands = await reveal(page.getByRole('checkbox', { name: 'Separate Strands', exact: true, includeHidden: true }));
    await strands.uncheck();
    const key = await placeLane(page, 'above');
    const placed = await layoutState(page, key);
    expect(placed).toMatchObject({ separate: false, linearStack: false, linear: { side: 'above' }, row: 'above' });

    // Use custom stack: Cancel keeps the checkbox off; Reset is one step.
    await custom.click();
    await expect(dialog(page)).toContainText('Changing Use custom stack to On leaves 1 Feature placement without its lane.');
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

    // The Circular panel's Separate Strands is the same draft input as Linear's.
    await switchMode(page, 'circular');
    const circularStrands = await reveal(page.locator('input[type=checkbox][aria-label="Separate Strands"]'));
    const switched = await layoutState(page, key);
    await circularStrands.click();
    await expect(dialog(page)).toContainText('Changing Separate Strands to On leaves 1 Linear Feature placement without its lane.');
    await expect(circularStrands).not.toBeChecked();
    await page.screenshot({ path: testInfo.outputPath('circular-separate-strands-dialog.png') });
    await cancelButton(page).click();
    await expect(dialog(page)).toBeHidden();
    await expect(circularStrands).not.toBeChecked();
    expect(await layoutState(page, key)).toEqual(switched);
    await circularStrands.click();
    await resetButton(page).click();
    await expect(circularStrands).toBeChecked();
    expect(await layoutState(page, key)).toEqual({ ...switched, separate: true, row: null, undo: switched.undo + 1 });
    await switchMode(page, 'linear');
    await generate(page);
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
