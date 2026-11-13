// A right-to-left page before any right-to-left language: the case sets
// <html dir="rtl"> over the English text (pseudo right to left), as setUiLocale
// does for such a language. Checked here:
// - the figure (its markup and the box of every text, from the figure's corner)
//   and the downloaded SVG, PNG and PDF equal those of left to right: the
//   Preview and the PDF staging read left to right in every language;
// - at 1024x768, 1280x800 and 1440x900, with a Result, with the Editor open and
//   with the settings as the icon rail, the page does not scroll sideways;
// - the Editor and its toggle sit at the left edge of the Preview, the History
//   list starts under the right edge of its toggle, the drawing menu stays
//   inside the window, and the text of a rail button stands left of it.
// Not checked here: a right-to-left catalog (none exists yet).
const { test, expect } = require('@playwright/test');
const { readFileSync } = require('node:fs');
const { join } = require('node:path');
const { openApp } = require('./helpers/app-lifecycle.cjs');

test.use({ actionTimeout: 20_000, viewport: { width: 1280, height: 800 } });

const SESSION = join(process.cwd(), 'gbdraw/web/gallery/sessions/HmmtDNA_basic_circular.gbdraw');
const SIZES = [[1024, 768], [1280, 800], [1440, 900]];

const loadSession = async (page) => {
  await page.locator('input[accept*=".gbdraw"]').setInputFiles(SESSION);
  await page.waitForFunction(() => !window.__GBDRAW_APP__.sessionImportPending
    && window.__GBDRAW_APP__.results.length > 0, null, { timeout: 180_000 });
};
const setDirection = async (page, dir) => {
  await page.evaluate(async (value) => {
    document.documentElement.dir = value;
    await new Promise((resolve) => requestAnimationFrame(() => requestAnimationFrame(resolve)));
  }, dir);
};
const button = (page, name) => page.locator(`[aria-label="${name}"]:visible`).first();

// The figure's markup and each text's box from the figure's top left corner.
const figure = (page) => page.locator('.result-pane svg[data-gbdraw-composition]').first().evaluate((root) => {
  const corner = root.getBoundingClientRect();
  return {
    markup: root.outerHTML,
    texts: Array.from(root.querySelectorAll('text, tspan'), (element) => {
      const box = element.getBoundingClientRect();
      return [box.left - corner.left, box.top - corner.top, box.width, box.height].map((value) => Math.round(value * 100) / 100);
    })
  };
});
const download = async (page, testInfo, action, name) => {
  const pending = page.waitForEvent('download', { timeout: 120_000 });
  await page.evaluate((method) => window.__GBDRAW_APP__[method](), action);
  const path = testInfo.outputPath(name);
  await (await pending).saveAs(path);
  return readFileSync(path);
};
// The PDF without the moment it was written.
const pdfBody = (bytes) => bytes.toString('latin1').replace(/\/(?:CreationDate|ModDate) \(D:[^)]*\)/g, '').replace(/\/ID \[[^\]]*\]/g, '');
const exports = async (page, testInfo, dir) => ({
  svg: (await download(page, testInfo, 'downloadSVG', `${dir}.svg`)).toString('utf8'),
  png: (await download(page, testInfo, 'downloadPNG', `${dir}.png`)).toString('base64'),
  pdf: pdfBody(await download(page, testInfo, 'downloadPDF', `${dir}.pdf`))
});

test('in a right-to-left page the figure and its SVG, PNG and PDF equal those of left to right', async ({ page }, testInfo) => {
  test.setTimeout(300_000);
  await openApp(page);
  await loadSession(page);
  const ltr = { figure: await figure(page), exports: await exports(page, testInfo, 'ltr') };
  await setDirection(page, 'rtl');
  expect(await figure(page)).toEqual(ltr.figure);
  const rtl = await exports(page, testInfo, 'rtl');
  expect(rtl.svg).toBe(ltr.exports.svg);
  expect(rtl.png).toBe(ltr.exports.png);
  expect(rtl.pdf).toBe(ltr.exports.pdf);
});

const scrollsSideways = (page) => page.evaluate(() => document.scrollingElement.scrollWidth > document.scrollingElement.clientWidth);
const boxOf = (locator) => locator.evaluate((element) => {
  const { left, right, top, bottom } = element.getBoundingClientRect();
  return { left, right, top, bottom };
});
const insideWindow = (box, size) => box.left >= 0 && box.top >= 0 && box.right <= size[0] && box.bottom <= size[1];

test('a right-to-left page does not scroll sideways, and its drawer, lists and texts stand on its side, at three window sizes', async ({ page }) => {
  test.setTimeout(600_000);
  await openApp(page);
  await loadSession(page);
  await setDirection(page, 'rtl');
  // A step that the History list can show.
  const prefix = page.getByLabel('Output Prefix', { exact: true }).first();
  await prefix.fill('HmmtDNA_basic_circular_2');
  await prefix.press('Tab');
  for (const size of SIZES) {
    await page.setViewportSize({ width: size[0], height: size[1] });
    await setDirection(page, 'rtl');
    expect(await scrollsSideways(page), `${size} with a Result`).toBe(false);

    await button(page, 'History').click();
    const list = page.locator('#history-step-list');
    await expect(list).toBeVisible();
    const toggle = await boxOf(button(page, 'History'));
    const shown = await boxOf(list);
    expect(insideWindow(shown, size), `${size} History list inside the window`).toBe(true);
    // Its right edge under the toggle's, unless the window's 8 px margin moves it.
    const width = shown.right - shown.left;
    expect(Math.abs(shown.right - Math.max(8 + width, Math.min(toggle.right, size[0] - 8))), `${size} History list under its toggle's right edge`).toBeLessThan(1);
    await page.keyboard.press('Escape');
    await expect(list).toBeHidden();

    await button(page, 'Drawing menu for Circular').click();
    const menu = page.locator('.drawing-bar-menu[open] .action-menu-items');
    await expect(menu).toBeVisible();
    expect(insideWindow(await boxOf(menu), size), `${size} drawing menu inside the window`).toBe(true);
    await page.keyboard.press('Escape');

    await button(page, 'Editor').click();
    const drawer = page.locator('.right-drawer');
    await expect(drawer).toHaveAttribute('aria-hidden', 'false');
    await page.waitForTimeout(400);
    const workspace = await boxOf(page.locator('.preview-workspace'));
    const open = await boxOf(drawer);
    expect(Math.abs(open.left - workspace.left), `${size} the Editor at the left edge`).toBeLessThan(1);
    expect((await boxOf(button(page, 'Editor'))).left, `${size} the toggle beside the open Editor`).toBeGreaterThanOrEqual(open.right - 1);
    expect(await scrollsSideways(page), `${size} with the Editor open`).toBe(false);
    await button(page, 'Editor').click();
    await expect(drawer).toHaveAttribute('aria-hidden', 'true');

    await button(page, 'Collapse settings to icons').click();
    const entry = page.locator('[data-settings-display="rail"] .settings-index nav > button').first();
    await entry.hover();
    const tip = page.locator('.settings-index-tip');
    await expect(tip).toBeVisible();
    expect((await boxOf(tip)).right, `${size} a rail text left of its button`).toBeLessThanOrEqual((await boxOf(entry)).left);
    expect(await scrollsSideways(page), `${size} with the icon rail`).toBe(false);
    await page.mouse.move(size[0] / 2, size[1] - 2);
    await button(page, 'Expand settings').click();
  }
});
