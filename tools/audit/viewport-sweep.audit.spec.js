// Viewport-width sweep: the app and the Gallery page must not scroll horizontally at phone,
// tablet, laptop, and desktop widths. Each step records the elements that extend past the
// viewport and saves a screenshot. Override the widths with AUDIT_WIDTHS=390,768.
const { readFileSync } = require('node:fs');
const { join } = require('node:path');
const { test, expect } = require('@playwright/test');
const A = require('./helpers/audit-common.cjs');

const widths = (fallback) => (process.env.AUDIT_WIDTHS
  ? process.env.AUDIT_WIDTHS.split(',').map(Number)
  : fallback);

for (const width of widths([390, 768, 1280, 1920])) {
  test(`app at ${width} px`, async ({ page }) => {
    test.setTimeout(240_000);
    const outdir = A.outDir('viewport');
    const log = A.collect(page);
    page.on('dialog', (d) => d.accept());
    await page.setViewportSize({ width, height: 900 });
    const { openApp, generateAndWaitForResult } = A.helpers();
    await openApp(page);
    const initial = await A.measureOverflow(page);
    const chooser = page.waitForEvent('filechooser');
    await page.getByRole('button', { name: 'Choose GenBank/DDBJ File', exact: true }).click();
    await (await chooser).setFiles({
      name: 'HmmtDNA.gbk',
      mimeType: 'text/plain',
      buffer: readFileSync(join(A.REPO_ROOT, 'tests/test_inputs/HmmtDNA.gbk'))
    });
    await expect.poll(() => page.evaluate(() => window.__GBDRAW_APP__.circularRecordList.length), { timeout: 60_000 })
      .toBeGreaterThan(0);
    const afterUpload = await A.measureOverflow(page);
    const gen = await generateAndWaitForResult(page, { expectedStatus: null, requireCommittedResult: false });
    await page.waitForTimeout(800);
    const afterGenerate = await A.measureOverflow(page);
    await page.screenshot({ path: join(outdir, `app-${width}.png`) });
    A.writeEvidence(outdir, `app-${width}.json`, {
      initial, afterUpload, afterGenerate, generate: gen.result?.status, log
    });
    console.log(width, 'overflow', initial.overflowX, afterUpload.overflowX, afterGenerate.overflowX);
  });
}

for (const width of widths([390, 768, 1280])) {
  test(`Gallery at ${width} px`, async ({ page }) => {
    test.setTimeout(240_000);
    const outdir = A.outDir('viewport');
    const log = A.collect(page);
    await page.setViewportSize({ width, height: 900 });
    await page.goto('/gbdraw/web/gallery/', { waitUntil: 'domcontentloaded' });
    await page.waitForSelector('#sample-list button');
    await page.waitForTimeout(1500);
    const preview = await A.measureOverflow(page);
    await page.screenshot({ path: join(outdir, `gallery-${width}.png`) });
    const tabs = {};
    for (const tab of ['Tutorial', 'Command', 'Files']) {
      await page.getByRole('tab', { name: tab, exact: true }).click();
      await page.waitForTimeout(1200);
      tabs[tab] = await A.measureOverflow(page);
      if (tabs[tab].overflowX) await page.screenshot({ path: join(outdir, `gallery-${width}-${tab}.png`) });
    }
    A.writeEvidence(outdir, `gallery-${width}.json`, { preview, tabs, log });
    console.log(width, 'overflow', preview.overflowX,
      Object.entries(tabs).map(([k, v]) => `${k}:${v.overflowX}`).join(','));
  });
}
