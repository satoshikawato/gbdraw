const { test, expect } = require('@playwright/test');
const path = require('node:path');
const fs = require('node:fs/promises');
const { execFile } = require('node:child_process');
const { promisify } = require('node:util');
const { generate } = require('./helpers/mode-transition.cjs');
const { assertDiagramWorkerIdle, evaluateWithRetainedPromise, openApp } = require('./helpers/app-lifecycle.cjs');
const { capture, assertCoherent } = require('./helpers/visual-state.cjs');

for (const composite of [false, true]) {
  test(`B-05 current CLI ${composite ? 'two-source composite' : 'single'} replay applies a definition edit on Generate`, async ({ browser }, info) => {
    test.setTimeout(300000);
    const sources = ['tests/fixtures/sessions/cli-web-mito.gb', ...(composite ? ['tests/test_inputs/NC_001416.gb'] : [])].map(file => path.resolve(file));
    let file = info.outputPath('cli.gbdraw-session.json.gz');
    await promisify(execFile)('python', ['-m', 'gbdraw.cli', 'circular', '--gbk', ...sources,
      '--track_type', 'middle', '--labels', 'out', ...(composite ? ['--multi_record_canvas'] : []),
      '-o', info.outputPath('cli'), '--session_output', file], { env: { ...process.env, PYTHONPATH: process.cwd() } });
    for (const phase of ['cli', 'web']) {
      const context = await browser.newContext({ viewport: { width: 1600, height: 1000 } });
      const page = await context.newPage();
      await context.route('**/*', route => new URL(route.request().url()).hostname === '127.0.0.1' ? route.continue() : route.abort());
      page.on('dialog', dialog => dialog.type() === 'prompt' ? dialog.accept('CLI replay') : dialog.accept());
      const definitionFontSizes = () => page.evaluate(async () => (
        Array.from((await import('./js/state.js')).state.svgContainer.value
          .querySelectorAll('g[data-gbdraw-role="record-definition"]'))
          .map(group => Number(group.querySelector('text').getAttribute('font-size')))
      ));
      try {
        await openApp(page);
        await page.locator('input[accept^=".json,"]').setInputFiles(file);
        await page.waitForFunction(() => !window.__GBDRAW_APP__.sessionImportPending && window.__GBDRAW_APP__.results.length > 0);
        if (phase === 'web') {
          // PD-OI-059: the definition font applies on Generate. Editing it
          // after Load starts no Worker work and leaves the saved Result
          // unchanged until Generate.
          const saved = await page.evaluate(() => window.__GBDRAW_APP__.results[0].content);
          await page.getByLabel('Titles and Record Labels', { exact: true }).press('Enter');
          await page.getByLabel('Default font size', { exact: true }).fill('19');
          const rewritten = await page.waitForFunction(
            saved => window.__GBDRAW_APP__.results[0].content !== saved, saved, { timeout: 3000 }
          ).then(() => true, () => false);
          expect(rewritten).toBe(false);
          await assertDiagramWorkerIdle(page, 'A definition edit after Load must not start the diagram Worker');
          await assertCoherent(await capture(page, info, 'web-loaded-definition-draft'), 'loaded definition', true);
        }
        await generate(page);
        await page.evaluate(async () => { await window.Vue.nextTick(); await new Promise(requestAnimationFrame); });
        await assertCoherent(await capture(page, info, `${phase}-definition-completed`), phase, true);
        if (phase === 'web') expect(await definitionFontSizes()).toEqual(Array(composite ? 2 : 1).fill(19));
        if (phase === 'cli') {
          const pending = page.waitForEvent('download');
          expect(await evaluateWithRetainedPromise(page, () => window.__GBDRAW_APP__.saveSessionWithTitle())).toMatchObject({ status: 'saved' });
          file = info.outputPath('web.gbdraw-session.json.gz');
          await (await pending).saveAs(file);
        }
      } finally { await page.context().close(); }
    }
  });
}
