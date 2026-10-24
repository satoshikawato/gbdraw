// The one place the specs switch the shown drawing (B4t): a click on the
// drawing's tab in the drawing bar, found by its drawing ID (never by the
// tab's text, which a rename or the UI language changes). `drawingIds` reads
// the store's order, so a spec never reads a tab to learn a drawing.
const { expect } = require('@playwright/test');

/** @param {import('@playwright/test').Page} page */
const shownDrawing = (page) => page.evaluate(async () => {
  const { state } = await import('./js/state.js');
  const drawing = state.activeDrawing();
  return { id: drawing.id, kind: drawing.mode };
});

/** @param {import('@playwright/test').Page} page @returns {Promise<string[]>} */
const drawingIds = (page) => page.evaluate(async () =>
  (await import('./js/state.js')).state.drawingList().map((drawing) => drawing.id));

/**
 * The tab of a drawing in the drawing bar.
 * @param {import('@playwright/test').Page} page @param {string} drawingId
 */
const drawingTab = (page, drawingId) => page.locator(`.drawing-bar [role="tab"][data-drawing-id="${drawingId}"]`);

/**
 * Show a drawing and wait until the app shows it. `kind` names the first
 * drawing of that kind in tab order, `name` the drawing with that stored name.
 * @param {import('@playwright/test').Page} page
 * @param {{ kind?: 'circular' | 'linear', name?: string }} target
 * @returns {Promise<string>} the shown DrawingId
 */
const showDrawing = async (page, { kind, name }) => {
  const drawingId = await page.evaluate(async ({ wantedKind, wantedName }) => {
    const { state } = await import('./js/state.js');
    const tab = state.drawingTabs.value.find((item) => (
      wantedName === undefined ? state.drawingOf(item.id).mode === wantedKind : item.name === wantedName));
    return tab ? tab.id : null;
  }, { wantedKind: kind, wantedName: name });
  if (!drawingId) throw new Error(`showDrawing: no drawing of kind ${kind} or named ${name}`);
  if (await drawingTab(page, drawingId).isVisible()) {
    await drawingTab(page, drawingId).click();
  } else {
    // At <= 768 px the tabs give way to the picker: its list holds every drawing in tab order.
    const picker = page.locator('.drawing-bar .drawing-picker');
    await picker.locator('summary').click();
    await picker.locator('.action-menu-item').nth((await drawingIds(page)).indexOf(drawingId)).click();
  }
  await expect.poll(async () => (await shownDrawing(page)).id).toBe(drawingId);
  return drawingId;
};

module.exports = { showDrawing, drawingTab, drawingIds, shownDrawing };
