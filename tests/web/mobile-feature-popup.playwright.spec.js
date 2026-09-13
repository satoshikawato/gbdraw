const { test, expect } = require('@playwright/test');
const { seeds, generate, popup, closeEditor, download } = require('./helpers/mode-transition.cjs');
const { load, agree, label } = require('./helpers/visual-state.cjs');

for (const mode of ['circular', 'linear']) {
  test(`B-03 ${mode} mobile search Edit Apply Label accepts pointer input`, async ({ browser }, info) => {
    test.setTimeout(300000);
    const page = await load(browser, mode === 'circular'
      ? 'gbdraw/web/gallery/sessions/HmmtDNA_ATskew.gbdraw-session.json' : seeds.linear, { width: 390, height: 844 });
    try {
      if (mode === 'circular') {
        await popup(page); await closeEditor(page);
        const search = page.getByRole('searchbox', { name: 'Search features', exact: true });
        await search.fill('tRNA'); await search.press('Enter');
        await page.getByRole('button', { name: 'Open active feature', exact: true }).click();
        await closeEditor(page);
      }
      await popup(page);
      await page.screenshot({ path: info.outputPath('popup-before.png') });
      await label(page, `MOBILE_${mode}`);
      for (const selector of ['[aria-label="Close feature popup"]', '[aria-label="Feature placement"]']) {
        const box = await page.locator('.feature-popup').locator(selector).boundingBox();
        expect(box.x).toBeGreaterThanOrEqual(0);
        expect(box.x + box.width).toBeLessThanOrEqual(390);
      }
      await agree(page, info, 'pointer-label-applied');
      await page.screenshot({ path: info.outputPath('popup-after.png') });
      await closeEditor(page);
      await generate(page);
      page.removeAllListeners('dialog');
      page.on('dialog', dialog => dialog.type() === 'prompt' ? dialog.accept('Mobile edit') : dialog.accept());
      const file = info.outputPath('edited.gbdraw-session.json.gz');
      await download(page, 'Save Session', file);
      const restored = await load(browser, file, { width: 390, height: 844 });
      try {
        expect(await restored.locator('.origin-top svg').textContent()).toContain(`MOBILE_${mode}`);
        await agree(restored, info, 'restored');
      } finally { await restored.context().close(); }
    } finally { await page.context().close(); }
  });
}
