// The one place the specs switch the shown drawing (B4t). Today the control is
// the mode button of the drawing's kind; the drawing bar replaces it by
// changing only this file. `drawingIds` reads the store's order, so a spec
// never reads a tab or a button to learn a drawing.
const { expect } = require('@playwright/test');

const KIND_LABEL = { circular: 'Circular', linear: 'Linear' };

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
 * Show a drawing and wait until the app shows one of that kind.
 * @param {import('@playwright/test').Page} page
 * @param {{ kind?: 'circular' | 'linear', name?: string }} target
 * @returns {Promise<string>} the shown DrawingId
 */
const showDrawing = async (page, { kind, name }) => {
  const wanted = kind || await page.evaluate(async (tabName) => {
    const { state } = await import('./js/state.js');
    const tab = state.drawingTabs.value.find((item) => item.name === tabName);
    return tab ? state.drawingOf(tab.id).mode : null;
  }, name);
  if (!KIND_LABEL[wanted]) throw new Error(`showDrawing: no drawing of kind ${kind} or named ${name}`);
  await page.getByRole('button', { name: KIND_LABEL[wanted], exact: true }).click();
  await expect.poll(() => page.evaluate(() => window.__GBDRAW_APP__.mode)).toBe(wanted);
  return (await shownDrawing(page)).id;
};

module.exports = { showDrawing, drawingIds, shownDrawing };
