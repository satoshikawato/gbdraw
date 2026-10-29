// The header Language menu (PD-OI-106), with the browser's language set to
// Japanese. Checked here:
// - a first visit shows English; the menu lists each language in itself, with its lang;
// - a choice of 日本語 shows Japanese at once (header, settings index, Generate),
//   sets <html lang>, drops upper-casing and tracking, and is stored; a reload
//   keeps it and its first mounted frames are already Japanese;
// - in ja, zh-Hans and zh-Hant: the Undo and Redo counts, the preview figure
//   (markup, and the computed font and spacing of every text) and the
//   downloaded SVG export equal those of English; so does a saved Session
//   (createdAt and random row ids aside) in ja and zh-Hans; the PDF export measures the figure as
//   English (its staging element inherits lang="en"; checked in ja);
// - English comes back with the same UI text and accessible names;
// - a Generate in ja and in zh-Hans writes the same Source recipe and Exact
//   replay as one in English;
// - in fr (a Latin-script language) and ko: the page language, the shown
//   catalog, the plural form the locale's rules pick (fr: `one` for 0 and 1),
//   the Hangul font stack and line breaking of ko, and, in fr, the same saved
//   Session and Source recipe as English;
// - in th: the shown catalog in the header and the icon rail, the Thai font
//   stack, and the 12 px floor of the session label and the rail names (a rail
//   name may take two lines and stays inside its entry);
// - at 1024x768, 1280x800 and 1440x900, in all 11 languages: every rail name
//   stays inside the rail column (inside its border and padding) and no word of
//   it is split over two lines (ja フィーチャー, fr Disposition), and the footer holds its text without a page
//   scroll, with nothing drawn over it and the settings pane, its Generate
//   bar and the Result pane ending above it (stacked: en, ja; rail: all 11)
//   (OV-1506, OV-1427, OV-1509); the page does not scroll sideways,
//   the header's mode and Session names keep one line (OV-1515), and each
//   footer separator stays on the line of the link before it (OV-1517);
// - storage that throws, or a catalog that fails to load, leaves English.
// - in every language, the "Suggest a translation fix" link beside the menu is
//   named in that language (English where the catalog has no entry), opens a
//   new tab, and links the translation feedback form with only the shown
//   language and the gbdraw version, with the language tag in the issue title
//   (PD-OI-106 scenario revision 5);
//   at 1024, 1280 and 1440 px the header keeps the menu and the link in view;
// Not checked here: PNG and PDF bytes, the shown Run Info text, and the
// diagnostics text.
const { test, expect } = require('@playwright/test');
const { readFileSync } = require('node:fs');
const { join } = require('node:path');
const { pathToFileURL } = require('node:url');
const { gunzipSync } = require('node:zlib');
const { generateAndWaitForResult, openApp } = require('./helpers/app-lifecycle.cjs');

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

// The SVG export as downloaded.
const exportedSvg = async (page, testInfo, name) => {
  const pending = page.waitForEvent('download', { timeout: 120_000 });
  await page.evaluate(() => window.__GBDRAW_APP__.downloadSVG());
  const path = testInfo.outputPath(name);
  await (await pending).saveAs(path);
  return readFileSync(path, 'utf8');
};
// The lang that the PDF export's staging copy of the figure inherits.
const pdfStagingLang = async (page) => {
  await page.evaluate(() => {
    window.__pdfStaging = [];
    window.__pdfObserver = new MutationObserver((records) => records.forEach((record) => record.addedNodes.forEach((node) => {
      if (node.nodeType === 1 && node.querySelector('svg')) window.__pdfStaging.push(node.closest('[lang]')?.getAttribute('lang') ?? null);
    })));
    window.__pdfObserver.observe(document.body, { childList: true });
  });
  const pending = page.waitForEvent('download', { timeout: 120_000 });
  await page.evaluate(() => window.__GBDRAW_APP__.downloadPDF());
  await pending;
  return page.evaluate(() => { window.__pdfObserver.disconnect(); return window.__pdfStaging; });
};
// Generate, then the Run Info's Source recipe and Exact replay.
const generatedRecipes = async (page) => {
  await generateAndWaitForResult(page);
  return page.evaluate(() => {
    const { sourceRecipe, exactReplay } = window.__GBDRAW_APP__.lastRunInfo;
    return JSON.stringify({ sourceRecipe, exactReplay });
  });
};

const loadSession = async (page) => {
  await page.locator('input[accept^=".json,"]').setInputFiles(join(process.cwd(), 'gbdraw/web/gallery/sessions', SESSION));
  await page.waitForFunction(() => !window.__GBDRAW_APP__.sessionImportPending
    && window.__GBDRAW_APP__.results.length > 0, null, { timeout: 180_000 });
};

// Save Session's bytes without its time stamp and random ids. The loaded Session has a
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
  // A load makes new random row ids; the zh-Hans save follows a reload.
  return JSON.stringify(session).replace(/[0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12}/g, '<uuid>');
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
  await expect(menu(page).locator('option')).toHaveText([
    'English', 'Bahasa Indonesia', 'Español', 'Français', 'Português', 'Tiếng Việt', 'ไทย', '日本語', '简体中文', '繁體中文', '한국어'
  ]);
  expect(await menu(page).locator('option').evaluateAll((options) => options.map((option) => option.lang)))
    .toEqual(['en', 'id', 'es', 'fr', 'pt', 'vi', 'th', 'ja', 'zh-Hans', 'zh-Hant', 'ko']);
  await expect(loadButton(page)).toHaveAccessibleName('Load Session');
  expect(await stored(page)).toEqual([]);

  await loadSession(page);
  const english = { figure: await figure(page), history: await history(page), ui: await uiText(page) };
  const englishSession = await savedSession(page, testInfo, 'en.json');
  const englishSvg = await exportedSvg(page, testInfo, 'en.svg');
  expect(englishSvg).toContain('<svg');

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
  expect(await exportedSvg(page, testInfo, 'ja.svg')).toBe(englishSvg);
  expect(await pdfStagingLang(page)).toEqual(['en']);
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
    expect(await exportedSvg(page, testInfo, `${tag}.svg`)).toBe(englishSvg);
    if (tag === 'zh-Hans') expect(await savedSession(page, testInfo, `${tag}.json`)).toBe(englishSession);
  }

  // English comes back exactly.
  await menu(page).selectOption('en');
  await expect.poll(() => pageLanguage(page)).toBe('en');
  expect(await uiText(page)).toEqual(english.ui);
  expect(await figure(page)).toEqual(english.figure);
  expect(await stored(page)).toEqual([[LANGUAGE_KEY, 'en']]);

  // Generate writes the same Source recipe and Exact replay in every language.
  const englishRecipes = await generatedRecipes(page);
  expect(englishRecipes).toContain('gbdraw');
  for (const tag of ['ja', 'zh-Hans']) {
    await menu(page).selectOption(tag);
    await expect.poll(() => pageLanguage(page)).toBe(tag);
    expect(await generatedRecipes(page)).toBe(englishRecipes);
  }
  // OV-1678: the Run info tab shows its notes and the Reproducibility files
  // description in the UI language, while the Run info keeps them in English.
  const zh = await catalogOf('zh-Hans');
  await page.getByRole('button', { name: zh['Run info'], exact: true }).click();
  const panel = page.locator('.run-info-panel');
  const { notes, description } = await page.evaluate(() => window.__GBDRAW_APP__.lastRunInfo)
    .then((info) => ({ notes: info.reproducibility.notes, description: info.downloadBundle.description }));
  expect(notes).toContain('The Source recipe can be rerun with the uploaded file names shown here.');
  await expect(panel.getByText(zh['The Source recipe can be rerun with the uploaded file names shown here.'], { exact: true })).toBeVisible();
  const shown = await panel.locator('p').allTextContents();
  for (const english of [...notes, description]) expect(shown.map((text) => text.trim())).not.toContain(english);
  await panel.screenshot({ path: testInfo.outputPath('run-info-zh-Hans.png') });
  expect(errors).toEqual([]);
});

// The app's own translation module (one instance per page): the form a count shows.
const recordCounts = (page) => page.evaluate(async () => {
  const { tn } = await import('./js/utils/i18n.js');
  return [0, 1, 2].map((count) => tn(count, '{count} record', '{count} records'));
});

test('French and Korean show their catalogs and plural forms, and French changes no Session or Source recipe', async ({ page }, testInfo) => {
  test.setTimeout(420_000);
  const [fr, ko] = await Promise.all(['fr', 'ko'].map(catalogOf));
  const errors = [];
  page.on('pageerror', (error) => errors.push(String(error)));
  page.on('dialog', (dialog) => dialog.accept());
  await openApp(page);
  await loadSession(page);
  const englishSession = await savedSession(page, testInfo, 'en.json');
  expect(await recordCounts(page)).toEqual(['0 records', '1 record', '2 records']);

  await menu(page).selectOption('fr');
  await expect.poll(() => pageLanguage(page)).toBe('fr');
  await expect(loadButton(page)).toHaveAccessibleName(fr['Load Session']);
  await expect(page.getByRole('button', { name: fr['Generate Diagram'], exact: true })).toBeVisible();
  // French takes the `one` form for 0 and 1.
  expect(await recordCounts(page)).toEqual([
    fr['{count} record'].replace('{count}', '0'), fr['{count} record'].replace('{count}', '1'), fr['{count} records'].replace('{count}', '2')
  ]);
  expect(await savedSession(page, testInfo, 'fr.json')).toBe(englishSession);
  // A Generate in French writes the Source recipe and Exact replay of one in English.
  const frenchRecipes = await generatedRecipes(page);
  expect(frenchRecipes).toContain('gbdraw');
  // With the Result in French: the drawer's Legend toolbar stays inside the
  // drawer (OV-1520), and the feature popup's tabs show their whole names, also
  // at phone width (OV-1527).
  const drawer = page.locator('.right-drawer');
  await page.locator('.drawer-toggle').click();
  await drawer.locator('[aria-controls="right-drawer-panel-legend"]').click();
  expect(await drawer.evaluate((panel) => {
    const toolbar = panel.querySelector('[data-legend-tab-toolbar]');
    const right = panel.getBoundingClientRect().right;
    return {
      scrolls: toolbar.scrollWidth > toolbar.clientWidth || panel.scrollWidth > panel.clientWidth,
      past: [...toolbar.children].filter((child) => child.getBoundingClientRect().right > right + 0.5).map((child) => child.textContent.trim())
    };
  })).toEqual({ scrolls: false, past: [] });
  await drawer.locator('[aria-controls="right-drawer-panel-features"]').click();
  await drawer.getByRole('button', { name: fr.Edit, exact: true }).first().click();
  const featurePopup = page.locator('.feature-popup');
  await expect(featurePopup).toBeVisible();
  await page.setViewportSize({ width: 390, height: 844 });
  await expect.poll(() => featurePopup.evaluate((element) => {
    const strip = element.querySelector('button:has(> .ph-list-bullets)').parentElement;
    const box = strip.getBoundingClientRect();
    return [...strip.querySelectorAll('button')].filter((tab) => {
      const name = tab.querySelector('span');
      const rect = tab.getBoundingClientRect();
      return name.scrollWidth > name.clientWidth || rect.left < box.left - 0.5 || rect.right > box.right + 0.5;
    }).map((tab) => tab.textContent.trim());
  })).toEqual([]);
  await page.setViewportSize({ width: 1280, height: 800 });
  await page.getByRole('button', { name: fr['Close feature popup'], exact: true }).click();
  await page.locator('.drawer-toggle').click();
  await menu(page).selectOption('en');
  await expect.poll(() => pageLanguage(page)).toBe('en');
  expect(await generatedRecipes(page)).toBe(frenchRecipes);

  await menu(page).selectOption('ko');
  await expect.poll(() => pageLanguage(page)).toBe('ko');
  await expect(loadButton(page)).toHaveAccessibleName(ko['Load Session']);
  // Korean has one form for every count.
  expect(await recordCounts(page)).toEqual([0, 1, 2].map((count) => ko['{count} records'].replace('{count}', String(count))));
  expect(await page.locator('body').evaluate((body) => {
    const style = getComputedStyle(body);
    return [style.fontFamily.includes('Apple SD Gothic Neo'), style.wordBreak];
  })).toEqual([true, 'keep-all']);
  expect(await page.locator('header .app-session-label').evaluate((element) => {
    const style = getComputedStyle(element);
    return [style.textTransform, style.letterSpacing];
  })).toEqual(['none', 'normal']);
  expect(await stored(page)).toEqual([[LANGUAGE_KEY, 'ko']]);
  expect(errors).toEqual([]);
});

test('Thai shows its catalog in the header and the rail, in a Thai font and at 12 px or more', async ({ page }) => {
  const th = await catalogOf('th');
  const errors = [];
  page.on('pageerror', (error) => errors.push(String(error)));
  await openApp(page);
  await menu(page).selectOption('th');
  await expect.poll(() => pageLanguage(page)).toBe('th');
  await expect(loadButton(page)).toHaveAccessibleName(th['Load Session']);
  const session = page.locator('header .app-session-label');
  await expect(session).toHaveText(th.Session);
  expect(await session.evaluate((element) => {
    const style = getComputedStyle(element);
    return [style.fontSize, style.textTransform, style.letterSpacing];
  })).toEqual(['12px', 'none', 'normal']);
  expect(await page.locator('body').evaluate((body) => getComputedStyle(body).fontFamily)).toContain('Leelawadee UI');

  await page.getByRole('button', { name: th['Collapse settings to icons'], exact: true }).click();
  const font = page.locator('.settings-index nav > button[data-settings-entry="font"]');
  await expect(font.locator('span')).toHaveText(th.Font);
  expect(await font.evaluate((button) => {
    const name = button.querySelector('span');
    const outer = button.getBoundingClientRect();
    const inner = name.getBoundingClientRect();
    return {
      size: getComputedStyle(name).fontSize,
      inside: inner.left >= outer.left - 0.5 && inner.right <= outer.right + 0.5 && inner.top >= outer.top - 0.5 && inner.bottom <= outer.bottom + 0.5
    };
  })).toEqual({ size: '12px', inside: true });
  expect(errors).toEqual([]);
});

// The text of a rail name, as line boxes, against the rail column inside its
// border and padding, and each word of it (the words of the page language)
// against the lines it lies on; the footer's text lines against what is drawn
// over them, the bottoms of the settings pane, its Generate bar and the Result
// pane against the footer's top, and the page height against the window.
const railAndFooterFit = (page) => page.evaluate(() => {
  const names = [...document.querySelectorAll('.settings-index nav > button > span')].filter((span) => span.offsetParent);
  const outside = names.flatMap((span) => {
    const column = span.closest('.settings-index');
    const box = column.getBoundingClientRect();
    const style = getComputedStyle(column);
    const range = document.createRange();
    range.selectNodeContents(span);
    const lines = [...range.getClientRects()].filter((rect) => rect.width > 0);
    const left = Math.min(...lines.map((rect) => rect.left));
    const right = Math.max(...lines.map((rect) => rect.right));
    const innerLeft = box.left + parseFloat(style.borderLeftWidth) + parseFloat(style.paddingLeft);
    const innerRight = box.right - parseFloat(style.borderRightWidth) - parseFloat(style.paddingRight);
    return left < innerLeft - 0.5 || right > innerRight + 0.5 ? [span.textContent.trim()] : [];
  });
  const words = new Intl.Segmenter(document.documentElement.lang, { granularity: 'word' });
  const splitWords = names.flatMap((span) => {
    const text = span.firstChild;
    const range = document.createRange();
    const lineTop = (index) => {
      range.setStart(text, index);
      range.setEnd(text, index + 1);
      return [...range.getClientRects()].find((rect) => rect.width > 0)?.top;
    };
    return [...words.segment(text.textContent)].filter((word) => word.isWordLike).flatMap((word) => {
      const tops = [...word.segment].map((_, offset) => lineTop(word.index + offset)).filter((top) => top !== undefined);
      return tops.some((top) => Math.abs(top - tops[0]) > 3) ? [word.segment] : [];
    });
  });
  const footer = document.querySelector('[data-app-footer]');
  const box = footer.getBoundingClientRect();
  const content = [...footer.children].map((child) => child.getBoundingClientRect());
  const coveredLines = [...footer.querySelectorAll(':scope > span')].flatMap((span) => {
    const range = document.createRange();
    range.selectNodeContents(span);
    return [...range.getClientRects()].filter((rect) => rect.width > 1).flatMap((rect) => {
      const hit = document.elementFromPoint(rect.left + rect.width / 2, rect.top + rect.height / 2);
      return hit && footer.contains(hit) ? [] : [String(hit?.className || hit?.tagName || 'nothing').split(' ')[0]];
    });
  });
  const below = ['.settings-pane', '.generate-bar', '.result-pane'].filter((selector) => {
    const element = document.querySelector(selector);
    // A fixed element (the Generate bar) has no offsetParent: shown means it has a box.
    return element && element.getClientRects().length > 0 && element.getBoundingClientRect().bottom > box.top + 0.5;
  });
  // The header's mode and Session names keep one line (OV-1515), and a footer
  // separator stays on the line of the link before it (OV-1517).
  const lineTops = (element) => {
    const range = document.createRange();
    const walker = document.createTreeWalker(element, NodeFilter.SHOW_TEXT);
    const tops = [];
    for (let node = walker.nextNode(); node; node = walker.nextNode()) {
      range.selectNodeContents(node);
      tops.push(...[...range.getClientRects()].filter((rect) => rect.width > 0.5).map((rect) => rect.top));
    }
    return tops.filter((top, index) => tops.findIndex((other) => Math.abs(other - top) < 3) === index);
  };
  const wrappedHeaderNames = [...document.querySelectorAll('.mode-switch > button, .app-session-label')]
    .filter((element) => element.getClientRects().length > 0 && lineTops(element).length > 1)
    .map((element) => element.textContent.trim());
  const separatorsAlone = [...footer.querySelectorAll('a')].flatMap((link) => {
    const next = link.nextSibling;
    const at = next?.nodeType === Node.TEXT_NODE ? next.textContent.indexOf('|') : -1;
    if (at < 0) return [];
    const range = document.createRange();
    range.setStart(next, at);
    range.setEnd(next, at + 1);
    const separator = range.getBoundingClientRect();
    const linkLines = lineTops(link);
    // A link name keeps one line (Satoshi / Kawato split before OV-1517), and its separator follows on that line.
    return linkLines.length === 1 && Math.abs(separator.top - linkLines[0]) < 3 ? [] : [link.textContent.trim()];
  });
  return {
    names: names.length,
    outside,
    splitWords,
    wrappedHeaderNames,
    separatorsAlone,
    pageScrollsX: document.documentElement.scrollWidth > innerWidth,
    footerHoldsText: Math.min(...content.map((rect) => rect.top)) >= box.top - 0.5 && Math.max(...content.map((rect) => rect.bottom)) <= box.bottom + 0.5,
    coveredLines,
    below,
    pageScrolls: document.documentElement.scrollHeight > innerHeight
  };
});

// Control names and automatic values in the settings pane: the Colors Apply and
// Add buttons keep their name on one line (OV-1516); a Japanese control name
// breaks between phrases, but an element that asks for break-all keeps it
// (OV-1518, OV-1522); the automatic value shown in an empty field, in a number
// field or a track slot, is a block box that ends a long value with an
// ellipsis (OV-1519, OV-1521); a closed select does too (OV-1523); the
// Chinese and Japanese wrapping rule leaves an element that asks for one line
// on one line (OV-1524 uses the text-wrap-style longhand).
const controlFit = (page) => page.evaluate(() => {
  const colors = document.querySelector('[data-settings-section="colors"]');
  colors.querySelectorAll('details').forEach((details) => { details.open = true; });
  // The text lines of an element (its icon aside).
  const lines = (element) => {
    const range = document.createRange();
    const walker = document.createTreeWalker(element, NodeFilter.SHOW_TEXT);
    const tops = [];
    for (let node = walker.nextNode(); node; node = walker.nextNode()) {
      range.selectNodeContents(node);
      tops.push(...[...range.getClientRects()].filter((rect) => rect.width > 0.5).map((rect) => rect.top));
    }
    return tops.filter((top, index) => tops.findIndex((other) => Math.abs(other - top) < 3) === index).length;
  };
  const buttons = [...colors.querySelectorAll('button:has(> .ph-magic-wand), button:has(> .ph-plus)')].filter((button) => button.getClientRects().length > 0);
  const probe = (tag, className) => {
    const element = document.createElement(tag);
    element.className = className;
    element.textContent = 'x';
    document.getElementById('app').appendChild(element);
    const style = getComputedStyle(element);
    const read = { display: style.display, textOverflow: style.textOverflow, wordBreak: style.wordBreak, wrapMode: style.textWrapMode };
    element.remove();
    return read;
  };
  const placeholders = [...document.querySelectorAll('.auto-value-placeholder')].filter((span) => span.getClientRects().length > 0);
  const slot = probe('span', 'track-slot-auto-placeholder');
  const select = getComputedStyle(document.querySelector('header select[data-ui-language]'));
  return {
    buttons: buttons.length,
    wrappedButtons: buttons.filter((button) => lines(button) > 1).map((button) => button.textContent.trim()),
    placeholders: placeholders.length > 0,
    placeholderBoxes: [...new Set([...placeholders.map((span) => getComputedStyle(span)), slot].map((style) => `${style.display}/${style.textOverflow}`))],
    selectOverflow: select.textOverflow,
    buttonBreak: probe('button', '').wordBreak,
    breakAll: probe('h2', 'break-all').wordBreak,
    oneLine: [probe('p', 'truncate').wrapMode, probe('summary', 'whitespace-nowrap').wrapMode]
  };
});

test('rail names stay inside the rail and keep each word on one line, nothing covers the footer, and control names fit, in every language at three window sizes', async ({ page }) => {
  test.setTimeout(480_000);
  const errors = [];
  page.on('pageerror', (error) => errors.push(String(error)));
  await page.setViewportSize({ width: 1024, height: 768 });
  await openApp(page);
  for (const tag of ['zh-Hans', 'ja']) {
    await menu(page).selectOption(tag);
    await expect.poll(() => pageLanguage(page)).toBe(tag);
    await expect.poll(() => controlFit(page), { message: `${tag} at 1024x768` }).toEqual({
      buttons: 3,
      wrappedButtons: [],
      placeholders: true,
      placeholderBoxes: ['block/ellipsis'],
      selectOverflow: 'ellipsis',
      buttonBreak: tag === 'ja' ? 'auto-phrase' : 'normal',
      breakAll: 'break-all',
      oneLine: ['nowrap', 'nowrap']
    });
  }
  const viewports = [{ width: 1024, height: 768 }, { width: 1280, height: 800 }, { width: 1440, height: 900 }];
  const fit = { outside: [], splitWords: [], wrappedHeaderNames: [], separatorsAlone: [], pageScrollsX: false, footerHoldsText: true, coveredLines: [], below: [], pageScrolls: false };
  for (const display of ['stacked', 'rail']) {
    if (display === 'rail') await page.getByRole('button', { name: 'Collapse settings to icons', exact: true }).click();
    for (const tag of display === 'rail' ? ['en', 'ja', 'zh-Hans', 'zh-Hant', 'ko', 'es', 'fr', 'pt', 'id', 'vi', 'th'] : ['en', 'ja']) {
      await menu(page).selectOption(tag);
      await expect.poll(() => pageLanguage(page)).toBe(tag);
      for (const viewport of viewports) {
        await page.setViewportSize(viewport);
        await expect.poll(() => railAndFooterFit(page), { message: `${display}, ${tag} at ${viewport.width}x${viewport.height}` })
          .toEqual({ ...fit, names: display === 'rail' ? 13 : 0 });
      }
    }
    await menu(page).selectOption('en');
    await page.setViewportSize(viewports[0]);
  }
  expect(errors).toEqual([]);
});

test('in every language, the link beside the menu opens the feedback form with only the shown language and the gbdraw version', async ({ page }) => {
  test.setTimeout(180_000);
  const version = /^version\s*=\s*"([^"]+)"/m.exec(readFileSync(join(process.cwd(), 'pyproject.toml'), 'utf8'))[1];
  const errors = [];
  page.on('pageerror', (error) => errors.push(String(error)));
  await openApp(page);
  const link = page.locator('header .app-language a[data-translation-feedback]');
  // Beside the menu: the next element after it.
  expect(await menu(page).evaluate((select) => select.nextElementSibling?.matches('a[data-translation-feedback]'))).toBe(true);
  const tags = await menu(page).locator('option').evaluateAll((options) => options.map((option) => option.value));
  expect(tags).toHaveLength(11);
  const minTitle = await page.evaluate(async () => (await import('./js/app/header-fit.js')).SESSION_TITLE_MIN_PX);
  expect(minTitle).toBeGreaterThan(0);
  for (const tag of tags) {
    await menu(page).selectOption(tag);
    await expect.poll(() => pageLanguage(page)).toBe(tag);
    const catalog = tag === 'en' ? {} : await catalogOf(tag);
    const name = catalog['Suggest a translation fix'] ?? 'Suggest a translation fix';
    await expect(link).toBeVisible();
    await expect(link).toHaveAccessibleName(name);
    await expect(link).toHaveAttribute('title', name);
    await expect(link).toHaveAttribute('target', '_blank');
    await expect(link).toHaveAttribute('rel', 'noopener noreferrer');
    const href = new URL(await link.getAttribute('href'));
    expect(href.origin + href.pathname).toBe('https://github.com/satoshikawato/gbdraw/issues/new');
    expect([...href.searchParams]).toEqual([['template', 'translation-feedback.yml'], ['language', tag], ['version', version], ['title', `[Translation] ${tag}: `]]);
    // The header keeps one row inside the window, with the menu and the link in view (OV-1495: es, fr and pt at 1440 px).
    // OV-2177, OV-2178: from 769 px the header stays inside the window and keeps the Session title
    // SESSION_TITLE_MIN_PX wide or more, on one row or, where that row cannot, on two (en overflowed at 769 px and
    // left the title 10 px wide to 810 px). The last steps narrow the window, which the header follows too.
    for (const width of [1024, 1440, 1280, 769, 800, 850, 800]) {
      await page.setViewportSize({ width, height: 800 });
      // The header chooses its state in the page's resize event (app/header-fit.js), which can run after
      // setViewportSize returns: let it settle at each width (the same state and title width in two frames), so the
      // next, narrower width starts from that state.
      await page.evaluate((target) => new Promise((done) => {
        let last = '';
        let frames = 0;
        const step = () => {
          const header = document.querySelector('header');
          const now = `${innerWidth} ${header.getAttributeNames().join(' ')} ${header.querySelector('.app-session-button').getBoundingClientRect().width}`;
          frames = now === last && innerWidth === target ? frames + 1 : 0;
          last = now;
          if (frames >= 2) done(); else requestAnimationFrame(step);
        };
        requestAnimationFrame(step);
      }), width);
      await expect.poll(() => page.locator('header').evaluate((header, minTitle) => {
        const right = Math.max(...Array.from(header.querySelectorAll('.app-language, a[data-translation-feedback]'), (element) => element.getBoundingClientRect().right));
        // The Session caption and title button paint inside the Session block, never under the mode switch (OV-1530:
        // fr at 1280 px). A box clips what it hides: a part counts only up to the right edge of a clipping ancestor.
        const block = header.querySelector('.app-session-actions');
        const clipped = (node, edge) => {
          let shown = edge;
          for (let element = node.parentElement; element && element !== block; element = element.parentElement) {
            if (getComputedStyle(element).overflowX !== 'visible') shown = Math.min(shown, element.getBoundingClientRect().right);
          }
          return shown;
        };
        const range = document.createRange();
        const parts = [...block.querySelectorAll('*')].filter((element) => element.getClientRects().length > 0)
          .map((element) => clipped(element, element.getBoundingClientRect().right));
        const walker = document.createTreeWalker(block, NodeFilter.SHOW_TEXT);
        for (let node = walker.nextNode(); node; node = walker.nextNode()) {
          range.selectNodeContents(node);
          parts.push(...[...range.getClientRects()].filter((rect) => rect.width > 0.5).map((rect) => clipped(node, rect.right)));
        }
        // The mode switch bounds the block where it follows it on its row (the two-row header puts it first).
        const box = block.getBoundingClientRect();
        const mode = header.querySelector('.mode-switch').getBoundingClientRect();
        const follows = mode.left > box.left && mode.top < box.bottom && box.top < mode.bottom;
        const limit = Math.min(box.right, follows ? mode.left : Infinity);
        return {
          inside: right <= document.documentElement.clientWidth,
          scrolls: document.documentElement.scrollWidth > document.documentElement.clientWidth,
          sessionInside: Math.max(...parts) <= limit + 0.5,
          titleUsable: header.querySelector('.app-session-button').getBoundingClientRect().width >= minTitle
        };
      }, minTitle), { message: `${tag} at ${width} px` }).toEqual({ inside: true, scrolls: false, sessionInside: true, titleUsable: true });
    }
  }
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
