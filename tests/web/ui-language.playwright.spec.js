// The header Language menu (PD-OI-106). A first visit shows English whatever
// the browser's language; a choice shows the UI in that language at once,
// without a reload, sets <html lang>, and is remembered per browser. It is no
// History step and changes no Result, setting, or Session: the drawn figure and
// a saved Session are the same bytes in every language. English comes back
// exactly. Storage that fails, or a catalog that fails to load, leaves English.
const { test, expect } = require('@playwright/test');
const { readFileSync } = require('node:fs');
const { join } = require('node:path');
const { pathToFileURL } = require('node:url');
const { gunzipSync } = require('node:zlib');
const { openApp } = require('./helpers/app-lifecycle.cjs');

test.use({ actionTimeout: 20_000, locale: 'ja-JP', viewport: { width: 1280, height: 800 } });

const LANGUAGE_KEY = 'gbdraw.ui-language';
const SESSION = 'HmmtDNA_basic_circular.gbdraw-session.json';
const CATALOG_ROOT = join(process.cwd(), 'gbdraw/web/js/utils/i18n-catalogs');
const catalogOf = async (tag) => (await import(pathToFileURL(join(CATALOG_ROOT, `${tag}.js`)).href)).default;

const menu = (page) => page.locator('header select[data-ui-language]');
const loadButton = (page) => page.locator('header [data-history-managed][aria-label]').filter({ has: page.locator('.ph-archive-tray') });
const pageLanguage = (page) => page.evaluate(() => document.documentElement.lang);
const stored = (page) => page.evaluate(() => Object.entries(localStorage));
const history = (page) => page.evaluate(() => ({
  undo: window.__GBDRAW_HISTORY__.getUndoCount(),
  redo: window.__GBDRAW_HISTORY__.getRedoCount(),
  undoEnabled: !document.querySelector('header .app-history-actions button').disabled
}));
// The figure as drawn: its markup and the computed style of its text (the
// pane around it may move, so its pixels are not compared).
const figure = async (page) => {
  const svg = page.locator('.result-pane svg[data-gbdraw-composition]').first();
  await expect(svg).toBeVisible();
  return svg.evaluate((root) => ({
    markup: root.outerHTML,
    text: Array.from(root.querySelectorAll('text, tspan'), (element) => {
      const style = getComputedStyle(element);
      return [style.fontFamily, style.fontSize, style.fontWeight, style.letterSpacing, style.textTransform, element.closest('[lang]')?.getAttribute('lang')].join('|');
    })
  }));
};
// The UI text a reader sees and hears: every text node and accessible attribute
// of the app. A notice keeps the language it was shown in (D-16), so it is
// dismissed first.
const uiText = async (page) => {
  const dismiss = page.locator('[data-app-notice] .app-notice-dismiss');
  if (await dismiss.count()) await dismiss.click();
  await expect(page.locator('[data-app-notice]')).toHaveCount(0);
  return page.evaluate(() => {
  const app = document.getElementById('app');
  const attributes = Array.from(app.querySelectorAll('[aria-label], [title], [placeholder], [alt]'), (element) => ['aria-label', 'title', 'placeholder', 'alt']
    .map((name) => element.getAttribute(name) ?? '').join('|'));
  return { text: app.innerText, attributes };
  });
};

const loadSession = async (page) => {
  await page.locator('input[accept^=".json,"]').setInputFiles(join(process.cwd(), 'gbdraw/web/gallery/sessions', SESSION));
  await page.waitForFunction(() => !window.__GBDRAW_APP__.sessionImportPending
    && window.__GBDRAW_APP__.results.length > 0, null, { timeout: 180_000 });
};

// Save Session's bytes without its time stamp. The loaded Session has a
// title, so Save asks for none; a second Save asks whether to download the
// same file again, and its primary action (in any language) does.
const savedSession = async (page, testInfo, name) => {
  const pending = page.waitForEvent('download', { timeout: 120_000 });
  await page.locator('header button[aria-describedby="session-save-help"]').click();
  const again = page.locator('[data-choice-prompt] button.btn-primary[data-choice-prompt-action]');
  const download = await Promise.race([pending, again.waitFor({ timeout: 120_000 }).then(async () => {
    await again.click();
    return pending;
  }, () => pending)]);
  const path = testInfo.outputPath(name);
  await download.saveAs(path);
  const bytes = readFileSync(path);
  const session = JSON.parse((bytes[0] === 0x1f && bytes[1] === 0x8b ? gunzipSync(bytes) : bytes).toString('utf8'));
  delete session.createdAt;
  return JSON.stringify(session);
};

test('a choice shows the UI in that language at once and changes no Result, History, or Session', async ({ page }, testInfo) => {
  test.setTimeout(420_000);
  const [ja, zhHans, zhHant] = await Promise.all(['ja', 'zh-Hans', 'zh-Hant'].map(catalogOf));
  const errors = [];
  page.on('pageerror', (error) => errors.push(String(error)));
  page.on('dialog', (dialog) => dialog.accept());

  // A first visit is English although the browser's language is Japanese.
  await openApp(page);
  expect(await page.evaluate(() => navigator.language)).toBe('ja-JP');
  await expect.poll(() => pageLanguage(page)).toBe('en');
  await expect(menu(page)).toHaveValue('en');
  await expect(menu(page)).toHaveAccessibleName('Language');
  await expect(menu(page).locator('option')).toHaveText(['English', '日本語', '简体中文', '繁體中文']);
  expect(await menu(page).locator('option').evaluateAll((options) => options.map((option) => option.lang)))
    .toEqual(['en', 'ja', 'zh-Hans', 'zh-Hant']);
  await expect(loadButton(page)).toHaveAccessibleName('Load Session');
  expect(await stored(page)).toEqual([]);

  await loadSession(page);
  const english = { figure: await figure(page), history: await history(page), ui: await uiText(page) };
  const englishSession = await savedSession(page, testInfo, 'en.json');

  await menu(page).selectOption('ja');
  await expect.poll(() => pageLanguage(page)).toBe('ja');
  await expect(loadButton(page)).toHaveAccessibleName(ja['Load Session']);
  await expect(loadButton(page)).toContainText(ja['Load Session']);
  await expect(menu(page)).toHaveAccessibleName(ja.Language);
  await expect(page.getByRole('navigation', { name: ja['Settings sections'] })
    .getByRole('button', { name: ja['Go to {section}'].replace('{section}', ja.Labels), exact: true })).toBeVisible();
  await expect(page.getByRole('button', { name: ja['Generate Diagram'], exact: true })).toBeVisible();
  // The diagram mode keeps the key 'Linear'; the Record Size Mode option has its own (OV-1387, i18n-catalogs.test.mjs).
  await expect(page.locator('header .mode-switch')).toContainText(ja.Linear);
  // No upper-casing or tracking in Japanese.
  expect(await page.locator('header .app-session-label').evaluate((element) => {
    const style = getComputedStyle(element);
    return [style.textTransform, style.letterSpacing];
  })).toEqual(['none', 'normal']);
  expect(await history(page)).toEqual(english.history);
  expect(await figure(page)).toEqual(english.figure);
  expect(await savedSession(page, testInfo, 'ja.json')).toBe(englishSession);
  expect(await stored(page)).toEqual([[LANGUAGE_KEY, 'ja']]);

  // A reload keeps Japanese, and the first frame the app paints is already Japanese.
  await page.addInitScript(() => {
    window.__firstNames = [];
    const read = () => {
      // The page holds the template until Vue mounts it (v-cloak).
      const mounted = window.__GBDRAW_APP__ && !document.getElementById('app')?.hasAttribute('v-cloak');
      const button = document.querySelector('header button[aria-describedby="session-save-help"]');
      if (mounted && button) window.__firstNames.push(button.getAttribute('aria-label'));
      if (window.__firstNames.length < 2) requestAnimationFrame(read);
    };
    requestAnimationFrame(read);
  });
  await openApp(page);
  await expect.poll(() => pageLanguage(page)).toBe('ja');
  await expect(menu(page)).toHaveValue('ja');
  expect(await page.evaluate(() => window.__firstNames)).toEqual([ja['Save Session'], ja['Save Session']]);
  await loadSession(page);
  const japanese = await history(page);

  for (const [tag, catalog] of [['zh-Hans', zhHans], ['zh-Hant', zhHant]]) {
    await menu(page).selectOption(tag);
    await expect.poll(() => pageLanguage(page)).toBe(tag);
    await expect(loadButton(page)).toHaveAccessibleName(catalog['Load Session']);
    await expect(page.getByRole('button', { name: catalog['Generate Diagram'], exact: true })).toBeVisible();
    expect(await page.locator('body').evaluate((body) => getComputedStyle(body).fontFamily)).not.toMatch(/^'?Inter'?, "?Noto Sans JP/);
    expect(await history(page)).toEqual(japanese);
    expect(await figure(page)).toEqual(english.figure);
  }

  // English comes back exactly.
  await menu(page).selectOption('en');
  await expect.poll(() => pageLanguage(page)).toBe('en');
  expect(await uiText(page)).toEqual(english.ui);
  expect(await figure(page)).toEqual(english.figure);
  expect(await stored(page)).toEqual([[LANGUAGE_KEY, 'en']]);
  expect(errors).toEqual([]);
});

test('storage that fails, or a catalog that cannot load, leaves English and the app usable', async ({ browser, baseURL }) => {
  test.setTimeout(300_000);
  const ja = await catalogOf('ja');
  const context = await browser.newContext({ baseURL, locale: 'ja-JP', viewport: { width: 1280, height: 800 } });
  try {
    // Storage that throws on every access: English, and a choice still applies to this page.
    const blocked = await context.newPage();
    const errors = [];
    blocked.on('pageerror', (error) => errors.push(String(error)));
    await blocked.addInitScript(() => {
      Object.defineProperty(window, 'localStorage', {
        configurable: true,
        get() { throw new DOMException('The operation is insecure.', 'SecurityError'); }
      });
    });
    await openApp(blocked);
    await expect.poll(() => pageLanguage(blocked)).toBe('en');
    await menu(blocked).selectOption('ja');
    await expect(loadButton(blocked)).toHaveAccessibleName(ja['Load Session']);
    await openApp(blocked);
    await expect.poll(() => pageLanguage(blocked)).toBe('en');
    await expect(loadButton(blocked)).toHaveAccessibleName('Load Session');
    expect(errors).toEqual([]);

    // A stored Japanese whose catalog fails to load: the app boots in English, and the menu says so.
    const offline = await context.newPage();
    offline.on('pageerror', (error) => errors.push(String(error)));
    await offline.route('**/js/utils/i18n-catalogs/*.js', (route) => route.abort());
    await offline.goto('/gbdraw/web/index.html', { waitUntil: 'domcontentloaded' });
    await offline.evaluate((key) => localStorage.setItem(key, 'ja'), LANGUAGE_KEY);
    await openApp(offline, { checkErrors: false });
    await expect.poll(() => pageLanguage(offline)).toBe('en');
    await expect(menu(offline)).toHaveValue('en');
    await expect(loadButton(offline)).toHaveAccessibleName('Load Session');
    // Choosing it again keeps English and the menu shows English.
    await menu(offline).selectOption('zh-Hans');
    await expect(menu(offline)).toHaveValue('en');
    await expect.poll(() => pageLanguage(offline)).toBe('en');
    expect(await offline.evaluate((key) => localStorage.getItem(key), LANGUAGE_KEY)).toBe('ja');
    expect(errors).toEqual([]);
  } finally {
    await context.close();
  }
});
