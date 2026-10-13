const { test, expect } = require('@playwright/test');
const { execFileSync } = require('node:child_process');
const { resolve } = require('node:path');
const {
  openApp, reveal, generateAndWaitForResult, expectSessionLoadedWithoutWarnings
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
  await page.locator('[data-font-file-input]').setInputFiles(fontFile);
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
  const drawn = await page.evaluate(() => {
    const content = String(window.__GBDRAW_APP__.results[0]?.content || '');
    const label = document.querySelector('.gbdraw-preview-surface svg text[font-family]');
    return {
      families: [...new Set([...content.matchAll(/font-family="([^"]*)"/g)].map((match) => match[1]))],
      legendFace: (content.match(/fontFile(?:&quot;|")\s*:\s*(?:&quot;|")([^&"]*)/) || [])[1] || null,
      previewText: label ? getComputedStyle(label).fontFamily : null
    };
  });
  expect(drawn.families).toEqual([FONT_LIST]);
  expect(drawn.legendFace).toMatch(/^user-[0-9a-f]{16}-Regular$/);
  expect(drawn.previewText).toContain(FAMILY);

});
