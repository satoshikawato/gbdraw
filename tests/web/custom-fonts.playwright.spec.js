const { test, expect } = require('@playwright/test');
const { execFileSync } = require('node:child_process');
const { resolve } = require('node:path');
const { readFileSync, writeFileSync } = require('node:fs');
const {
  openApp, reveal, generateAndWaitForResult, expectSessionLoadedWithoutWarnings, assertSessionLoadLeftWorkerIdle
} = require('./helpers/app-lifecycle.cjs');

// Custom fonts (D-01, D-06): a font file chosen in the Font section is
// registered in the page and in the diagram Worker, the drawing's font becomes
// its CSS list, and Generate draws and measures the diagram with it.

const SESSION = resolve('gbdraw/web/gallery/sessions/HmmtDNA_basic_circular.gbdraw-session.json');
const FAMILY = 'Gbdraw Test Serif';
const FONT_LIST = `'${FAMILY}', serif`;

// A bundled Liberation Serif face renamed, so the test needs no system font.
const renamedFont = (target) => {
  execFileSync('python', ['-c', `
import sys
from fontTools.ttLib import TTFont
font = TTFont("gbdraw/data/LiberationSerif-Regular.ttf")
names = font["name"]
names.names = [record for record in names.names if record.nameID not in {1, 2, 4, 6, 16, 17}]
for name_id, value in ((1, sys.argv[1]), (2, "Regular"), (4, sys.argv[1] + " Regular"), (6, sys.argv[1].replace(" ", "")), (16, sys.argv[1])):
    names.setName(value, name_id, 3, 1, 0x409)
font.save(sys.argv[2])
`, FAMILY, target], { cwd: resolve('.'), stdio: 'pipe' });
  return target;
};

const drawnFonts = (page) => page.evaluate(() => {
  const content = String(window.__GBDRAW_APP__.results[0]?.content || '');
  const label = document.querySelector('.gbdraw-preview-surface svg text[font-family]');
  return {
    families: [...new Set([...content.matchAll(/font-family="([^"]*)"/g)].map((match) => match[1]))],
    legendFace: (content.match(/fontFile(?:&quot;|")\s*:\s*(?:&quot;|")([^&"]*)/) || [])[1] || null,
    previewText: label ? getComputedStyle(label).fontFamily : null
  };
});

test('a font file chosen in Font is one History step and draws and measures the Generated diagram', async ({ page }, testInfo) => {
  test.setTimeout(300_000);
  const fontFile = renamedFont(testInfo.outputPath('gbdraw-test-serif.ttf'));
  await openApp(page);
  await page.locator('input[type="file"][accept^=".json,"]').setInputFiles(SESSION);
  await expectSessionLoadedWithoutWarnings(page, { timeout: 180_000 });

  const select = page.locator('[data-font-select]');
  await reveal(select);
  await expect(select).toHaveValue('default');
  await expect(page.locator('[data-font-custom-note]')).toHaveCount(0);
  // R8: the select's Load font file… opens the file picker.
  const chooser = page.waitForEvent('filechooser');
  await select.selectOption('load-file');
  await (await chooser).setFiles(fontFile);
  // Registration starts the diagram Worker's Python when it is not running yet.
  await expect(select).toHaveValue(FONT_LIST, { timeout: 180_000 });
  await expect(select.locator('option:checked')).toHaveText(FAMILY);
  await expect(page.locator('[data-font-custom-note]')).toBeVisible();
  await expect(page.locator('[data-font-load-error]')).toHaveCount(0);
  expect(await page.evaluate((family) => document.fonts.check(`12px "${family}"`), FAMILY)).toBe(true);

  // Loading and choosing the font is one History step; the font stays loaded.
  await page.getByRole('button', { name: 'Undo', exact: true }).click();
  await expect(select).toHaveValue('default');
  await expect(page.locator('[data-font-custom-note]')).toHaveCount(0);
  await expect(select.locator('option', { hasText: FAMILY })).toHaveCount(1);
  await page.getByRole('button', { name: 'Redo', exact: true }).click();
  await expect(select).toHaveValue(FONT_LIST);

  await generateAndWaitForResult(page);
  const drawn = await drawnFonts(page);
  expect(drawn.families).toEqual([FONT_LIST]);
  expect(drawn.legendFace).toMatch(/^user-[0-9a-f]{16}-Regular$/);
  expect(drawn.previewText).toContain(FAMILY);

  // R1: a new diagram Worker gets the page's fonts again, so a Generate after
  // the Worker is replaced still measures with the loaded face.
  await page.evaluate(async () => {
    const { disposeDiagramGenerationWorker } = await import('/gbdraw/web/js/services/diagram-generation.js');
    disposeDiagramGenerationWorker();
  });
  await generateAndWaitForResult(page);
  expect((await drawnFonts(page)).legendFace).toBe(drawn.legendFace);
  await expect(page.locator('[data-app-notice]', { hasText: 'could not be prepared again' })).toHaveCount(0);
});

// Receipt 〔読み込むとき〕 and 〔警告〕: a Session whose font is neither bundled
// nor loaded asks before it replaces the current work, without the diagram
// Worker. Cancel loads nothing; Load with <nearest> loads it, keeps the name in
// Font, draws the nearest font, and a Save says the Session keeps only the name.
test('a Session with a font that is not available asks first, and Load with the nearest font keeps its name', async ({ page }, testInfo) => {
  test.setTimeout(300_000);
  const MISSING = "'Garamond Premier Pro', serif";
  const session = JSON.parse(readFileSync(SESSION, 'utf8'));
  session.modes.circular.config.adv.font_family = MISSING;
  session.renderRequest.diagramOptions.configOverrides = {
    ...session.renderRequest.diagramOptions.configOverrides, 'objects.text.font_family': MISSING
  };
  const sessionFile = testInfo.outputPath('garamond.gbdraw-session.json');
  writeFileSync(sessionFile, JSON.stringify(session));
  await openApp(page);
  const sessionInput = page.locator('input[type="file"][accept^=".json,"]');
  const prompt = page.locator('[data-choice-prompt]');
  const select = page.locator('[data-font-select]');

  await sessionInput.setInputFiles(sessionFile);
  await expect(prompt.getByRole('heading', { name: 'This Session uses a font that is not available' })).toBeVisible({ timeout: 120_000 });
  await expect(prompt).toContainText('Its text is set in Garamond Premier Pro, which is not part of gbdraw and is not loaded in this page.');
  await expect(prompt).toContainText('Load with Liberation Serif: the Session opens and its text is drawn in Liberation Serif, the nearest available font. The Font setting keeps the name Garamond Premier Pro, so the font is used again once it is loaded. The saved preview keeps its look until you Generate.');
  await expect(prompt.locator('[data-choice-prompt-action="file"]')).toHaveText('Load font file…');
  await prompt.getByRole('button', { name: 'Cancel', exact: true }).click();
  await expect(prompt).toHaveCount(0);
  await expect(page.locator('[data-app-notice]', { hasText: 'Session loaded' })).toHaveCount(0);
  await reveal(select);
  await expect(select).toHaveValue('default');
  await assertSessionLoadLeftWorkerIdle(page);

  await sessionInput.setInputFiles(sessionFile);
  await prompt.locator('[data-choice-prompt-action="substitute"]').click();
  await expectSessionLoadedWithoutWarnings(page, { timeout: 180_000 });
  await expect(select).toHaveValue(MISSING);
  await expect(select.locator('option:checked')).toHaveText('Garamond Premier Pro — not available, drawn in Liberation Serif');
  await expect(page.locator('[data-font-custom-note]')).toBeVisible();
  await assertSessionLoadLeftWorkerIdle(page);

  await generateAndWaitForResult(page);
  expect((await drawnFonts(page)).legendFace).toBe('LiberationSerif-Regular');
  await expect(select).toHaveValue(MISSING);

  await page.getByRole('button', { name: 'Save Session', exact: true }).first().click();
  await expect(page.locator('[data-app-notice]')).toContainText('Session saved', { timeout: 120_000 });
  await expect(page.locator('[data-app-notice] li')).toHaveText([
    'The Session keeps the font name Garamond Premier Pro, not the font file. To open it elsewhere, the same font must be available there, or it is replaced with the nearest available font.'
  ]);
});
