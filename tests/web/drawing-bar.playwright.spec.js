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

  test('the drawing bar reaches every tab without horizontal page scroll and keeps the compact Editor', async ({ page }) => {
    test.setTimeout(240_000);
    await openApp(page);
    expect((await loadGallerySession(page, 'HmmtDNA_ATskew.gbdraw-session.json')).status).toBe('ok');
    const ids = await drawingIds(page);
    const shown = (await shownDrawing(page)).id;
    const tabs = page.locator('.drawing-bar [role="tab"]');
    await expect(tabs).toHaveCount(ids.length);
    expect(await tabs.evaluateAll((elements) => elements.map((tab) => tab.dataset.drawingId)), 'tab order').toEqual(ids);
    expect(await page.evaluate(() => document.documentElement.scrollWidth)).toBeLessThanOrEqual(390);

    // Every tab is in view and takes a tap; the keyboard reaches each one (manual activation).
    for (const id of ids) {
      const box = await drawingTab(page, id).evaluate((tab) => {
        tab.scrollIntoView({ block: 'nearest', inline: 'nearest' });
        const rect = tab.getBoundingClientRect();
        return { left: rect.left, right: rect.right, hit: tab.contains(document.elementFromPoint(
          rect.left + rect.width / 2, rect.top + rect.height / 2)) };
      });
      expect(box.left, id).toBeGreaterThanOrEqual(0);
      expect(box.right, id).toBeLessThanOrEqual(390);
      expect(box.hit, `${id} takes a tap`).toBe(true);
    }
    const other = ids.find((id) => id !== shown);
    await drawingTab(page, shown).focus();
    await page.keyboard.press(ids.indexOf(other) > ids.indexOf(shown) ? 'ArrowRight' : 'ArrowLeft');
    await expect(drawingTab(page, other)).toBeFocused();
    expect((await shownDrawing(page)).id, 'moving the focus does not switch').toBe(shown);
    await page.keyboard.press('Enter');
    await expect.poll(async () => (await shownDrawing(page)).id).toBe(other);
    await expect(drawingTab(page, other)).toHaveAttribute('aria-selected', 'true');
    await drawingTab(page, shown).tap();
    await expect.poll(async () => (await shownDrawing(page)).id).toBe(shown);

    // PD-OI-038 at 390 x 844: with the Editor open, the canvas spans the Result card and shows at least
    // 200 px between the header and the Generate bar, the Editor docks under it, and its toggle takes a tap.
    await page.locator('.drawer-toggle').click();
    await expect(page.locator('.right-drawer')).toBeVisible();
    const geometry = await page.getByRole('region', { name: 'Result Preview' }).evaluate(async (preview) => {
      await Promise.all(document.getAnimations().map((animation) => animation.finished.catch(() => {})));
      const canvas = preview.querySelector('.preview-canvas');
      canvas.scrollIntoView({ block: 'start' });
      const header = document.querySelector('.app-header').getBoundingClientRect();
      window.scrollBy(0, canvas.getBoundingClientRect().top - header.bottom);
      const box = (element) => element.getBoundingClientRect();
      const style = getComputedStyle(preview);
      const generate = document.querySelector('.generate-bar');
      const bottom = getComputedStyle(generate).position === 'fixed' ? box(generate).top : innerHeight;
      const toggle = box(document.querySelector('.drawer-toggle'));
      return {
        canvasLeft: box(canvas).left - box(preview).left - parseFloat(style.borderLeftWidth),
        canvasRight: box(preview).right - parseFloat(style.borderRightWidth) - box(canvas).right,
        visibleHeight: Math.min(box(canvas).bottom, bottom) - Math.max(box(canvas).top, box(document.querySelector('.app-header')).bottom),
        drawerBelow: box(document.querySelector('.right-drawer')).top >= box(canvas).bottom - 0.5,
        toggleHit: document.querySelector('.drawer-toggle').contains(document.elementFromPoint(
          toggle.left + toggle.width / 2, toggle.top + toggle.height / 2)),
        pageWidth: document.documentElement.scrollWidth
      };
    });
    expect(geometry).toMatchObject({ canvasLeft: 0, canvasRight: 0, drawerBelow: true, toggleHit: true });
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
