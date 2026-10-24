// The drawing bar (PD-OI-090 〔drawing bar と切り替え〕): one tab per drawing,
// named by the drawing, under the header and above both panes. It works at
// 390 px with the compact Preview and Editor conditions of PD-OI-038 (OIC-025),
// and its tab menu renames a drawing as one History step.
const { test, expect } = require('@playwright/test');
const { openApp, evaluateWithRetainedPromise } = require('./helpers/app-lifecycle.cjs');
const { drawingIds, drawingTab, shownDrawing } = require('./helpers/drawing-bar.cjs');

const loadGallerySession = (page, filename) => evaluateWithRetainedPromise(page, async (path) => {
  const name = path.split('/').pop();
  const response = await fetch(path);
  if (!response.ok) throw new Error(`Could not load ${name}: ${response.status}`);
  const file = new File([await response.arrayBuffer()], name, { type: 'application/json' });
  return window.__GBDRAW_APP__.importSession({ target: { files: [file], value: 'selected' } });
}, `/gbdraw/web/gallery/sessions/${filename}`);

const tabName = (page, drawingId) => page.evaluate(async (id) => (
  (await import('./js/state.js')).state.drawingTabs.value.find((tab) => tab.id === id).name), drawingId);

test.describe('390 px', () => {
  test.use({ viewport: { width: 390, height: 844 }, hasTouch: true });

  test('the phone layout shows one picker row and one pane at a time, switched by reveal, Generate and the Editor', async ({ page }) => {
    test.setTimeout(240_000);
    await openApp(page);
    expect((await loadGallerySession(page, 'HmmtDNA_ATskew.gbdraw-session.json')).status).toBe('ok');
    const ids = await drawingIds(page);
    const shown = (await shownDrawing(page)).id;
    const other = ids.find((id) => id !== shown);
    const picker = page.locator('.drawing-bar .drawing-picker');
    const paneTab = (pane) => page.locator(`.phone-pane-tabs [role="tab"][data-phone-pane="${pane}"]`);
    const settings = page.locator('#settings-sections');
    const result = page.locator('#result-pane');
    await expect(page.locator('.drawing-bar-tabs')).toBeHidden();
    await expect(picker).toBeVisible();
    expect(await page.evaluate(() => document.documentElement.scrollWidth)).toBeLessThanOrEqual(390);

    // The picker names the shown drawing and lists every drawing in tab order; a tap shows another.
    await picker.locator('summary').tap();
    const items = picker.locator('.action-menu-item');
    await expect(items).toHaveCount(ids.length);
    await items.nth(ids.indexOf(other)).tap();
    await expect.poll(async () => (await shownDrawing(page)).id).toBe(other);
    // The keyboard reaches it too: Arrow Down opens the list on its first item.
    await picker.locator('summary').focus();
    await page.keyboard.press('ArrowDown');
    for (let index = 0; index < ids.indexOf(shown); index += 1) await page.keyboard.press('ArrowDown');
    await page.keyboard.press('Enter');
    await expect.poll(async () => (await shownDrawing(page)).id).toBe(shown);
    const menuBox = await page.locator('.drawing-bar-menu').evaluate(async (menu) => {
      menu.open = true;
      await new Promise(requestAnimationFrame);
      const box = menu.querySelector('.action-menu-items').getBoundingClientRect();
      menu.open = false;
      return { left: box.left, right: box.right };
    });
    expect(menuBox.left, 'the Drawing menu stays on screen').toBeGreaterThanOrEqual(0);
    expect(menuBox.right).toBeLessThanOrEqual(390);

    // Settings | Preview: one pane at a time, the ARIA tab pattern with manual activation.
    await expect(paneTab('settings')).toHaveAttribute('aria-selected', 'true');
    await expect(result).toBeHidden();
    await paneTab('preview').tap();
    await expect(settings).toBeHidden();
    await expect(result).toBeVisible();
    await paneTab('preview').focus();
    await page.keyboard.press('ArrowLeft');
    await expect(paneTab('settings')).toBeFocused();
    await page.keyboard.press('Enter');
    await expect(paneTab('settings')).toHaveAttribute('aria-selected', 'true');
    await expect(settings).toBeVisible();

    // A path that reveals a setting shows Settings; Generate and opening the Editor show the Preview.
    await paneTab('preview').tap();
    await page.evaluate(() => window.__GBDRAW_APP__.revealSettingsSection('input-genomes'));
    await expect(paneTab('settings')).toHaveAttribute('aria-selected', 'true');
    await expect(settings).toBeVisible();
    await page.locator('.generate-bar button').click();
    await expect(paneTab('preview')).toHaveAttribute('aria-selected', 'true');
    await page.waitForFunction(() => !window.__GBDRAW_APP__.processing, null, { timeout: 180_000 });
    await expect(page.getByRole('region', { name: 'Result Preview' })).toBeVisible();
    // The approved images (Owner R17-10 A): the Result comes in with its workspace at the pane's top, and the
    // canvas with the controls below it fills the pane; the Result header and the search lie above it.
    const workspaceFit = () => page.evaluate(() => {
      const box = (selector) => document.querySelector(selector).getBoundingClientRect();
      return { top: Math.round(box('.preview-workspace').top - box('#result-pane').top),
        end: Math.round(box('#result-pane').bottom - box('.preview-controls').bottom) };
    });
    await expect.poll(workspaceFit).toEqual({ top: 0, end: 0 });
    await paneTab('settings').tap();
    await page.evaluate(() => window.__GBDRAW_APP__.openRightDrawerTab('features'));
    await expect(paneTab('preview')).toHaveAttribute('aria-selected', 'true');
    expect(await page.evaluate(() => document.documentElement.scrollWidth)).toBeLessThanOrEqual(390);

    // PD-OI-038 at 390 x 844: with the Editor open, the canvas spans the Result card and shows at least
    // 200 px between the header and the Generate bar, the Editor docks under it, and its toggle takes a tap.
    await expect(page.locator('.right-drawer')).toBeVisible();
    const geometry = await page.getByRole('region', { name: 'Result Preview' }).evaluate(async (preview) => {
      await Promise.all(document.getAnimations().map((animation) => animation.finished.catch(() => {})));
      const canvas = preview.querySelector('.preview-canvas');
      const box = (element) => element.getBoundingClientRect();
      const style = getComputedStyle(preview);
      const generate = document.querySelector('.generate-bar');
      const bottom = getComputedStyle(generate).position === 'fixed' ? box(generate).top : innerHeight;
      const toggle = box(document.querySelector('.drawer-toggle'));
      return {
        canvasLeft: box(canvas).left - box(preview).left - parseFloat(style.borderLeftWidth),
        canvasRight: box(preview).right - parseFloat(style.borderRightWidth) - box(canvas).right,
        visibleHeight: Math.min(box(canvas).bottom, bottom) - Math.max(box(canvas).top, box(document.querySelector('.phone-pane-tabs')).bottom),
        drawerBelow: box(document.querySelector('.right-drawer')).top >= box(canvas).bottom - 0.5,
        toggleHit: document.querySelector('.drawer-toggle').contains(document.elementFromPoint(
          toggle.left + toggle.width / 2, toggle.top + toggle.height / 2)),
        pageWidth: document.documentElement.scrollWidth,
        // Opening the Editor brings the workspace to the pane's top: the canvas keeps 13rem, the whole dock
        // ends above the Generate bar, and its feature rows start in view.
        workspaceTop: Math.round(box(preview.querySelector('.preview-workspace')).top - box(document.querySelector('#result-pane')).top),
        canvasHeight: Math.round(box(canvas).height),
        dockAboveGenerate: box(document.querySelector('.right-drawer')).bottom <= bottom,
        rowsInView: [...document.querySelectorAll('.right-drawer [data-feature-list-row]')]
          .filter((row) => box(row).top >= box(canvas).bottom && box(row).bottom <= bottom).length
      };
    });
    expect(geometry).toMatchObject({ canvasLeft: 0, canvasRight: 0, drawerBelow: true, toggleHit: true,
      workspaceTop: 0, canvasHeight: 13 * 16, dockAboveGenerate: true });
    expect(geometry.rowsInView).toBeGreaterThanOrEqual(2);
    expect(geometry.visibleHeight).toBeGreaterThanOrEqual(200);
    expect(geometry.pageWidth).toBeLessThanOrEqual(390);
  });
});

test('Rename in the tab menu is one History step that Undo restores; a repeated name is refused with its reason', async ({ page }) => {
  test.setTimeout(180_000);
  await page.setViewportSize({ width: 1280, height: 900 });
  await openApp(page);
  const [first, second] = await drawingIds(page);
  await drawingTab(page, second).click();
  await expect.poll(async () => (await shownDrawing(page)).id).toBe(second);
  const name = await tabName(page, second);
  const firstName = await tabName(page, first);
  const menu = page.locator('.drawing-bar-menu');
  const dialog = page.getByRole('dialog', { name: 'Rename drawing' });
  const rename = async (value) => {
    await menu.locator('summary').click();
    await menu.getByRole('button', { name: 'Rename…', exact: true }).click();
    await expect(dialog.getByRole('textbox', { name: 'Drawing name' })).toHaveValue(name);
    await dialog.getByRole('textbox', { name: 'Drawing name' }).fill(value);
    await dialog.getByRole('button', { name: 'Rename', exact: true }).click();
  };

  await rename(firstName.toUpperCase());
  await expect(dialog).toContainText('Another drawing has this name.');
  await expect(dialog.getByRole('textbox', { name: 'Drawing name' })).toHaveValue(firstName.toUpperCase());
  await dialog.getByRole('button', { name: 'Cancel', exact: true }).click();
  await expect(dialog).toHaveCount(0);
  expect(await tabName(page, second), 'a refused rename changes nothing').toBe(name);

  await rename('Phage comparison');
  await expect(dialog).toHaveCount(0);
  await expect(drawingTab(page, second)).toHaveAccessibleName('Phage comparison');
  await expect(page.getByRole('tab', { name: 'Phage comparison', exact: true })).toHaveAttribute('aria-selected', 'true');
  await expect(menu.locator('summary')).toHaveAttribute('aria-label', 'Drawing menu for Phage comparison');
  expect(await page.evaluate(() => window.__GBDRAW_HISTORY__.undoLabel())).toBe('Rename drawing');

  await page.getByRole('button', { name: 'Undo', exact: true }).click();
  await expect(drawingTab(page, second)).toHaveAccessibleName(name);
  expect(await tabName(page, second)).toBe(name);
  expect((await shownDrawing(page)).id, 'Undo of a rename shows the same drawing').toBe(second);
});
