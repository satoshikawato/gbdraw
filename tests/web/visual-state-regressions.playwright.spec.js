const { test, expect } = require('@playwright/test');
const { seeds, generate, popup, closeEditor, switchMode } = require('./helpers/mode-transition.cjs');
const { openApp } = require('./helpers/app-lifecycle.cjs');
const { load, agree, reveal, check, color, label } = require('./helpers/visual-state.cjs');

for (const mode of ['circular', 'linear']) {
  test(`B-01 ${mode} color Redo publishes the completed SVG`, async ({ browser }, info) => {
    test.setTimeout(300000);
    const page = await load(browser, seeds[mode]);
    try {
      await generate(page);
      const target = await color(page);
      await agree(page, info, 'direct-color');
      await closeEditor(page);
      for (const action of ['Undo', 'Redo']) {
        await page.getByRole('button', { name: action, exact: true }).click();
        const observation = await agree(page, info, action);
        if (action === 'Redo') {
          const parts = observation.completed.mounted.filter(node => node.feature === target.featureId && node.part === 'block');
          expect(parts.length).toBeGreaterThan(0);
          expect(parts.every(node => node.attrs.fill === '#c83366')).toBe(true);
        }
      }
    } finally { await page.context().close(); }
  });
}

for (const mode of ['circular', 'linear']) {
  test(`B-04 ${mode} depth OFF and ON publish parent visibility`, async ({ browser }, info) => {
    test.setTimeout(300000);
    const page = await load(browser, seeds[mode]);
    try {
      const length = mode === 'circular' ? 16569 : 48502;
      const depth = Buffer.from(Array.from({ length: Math.ceil(length / 100) }, (_, i) => `${mode === 'circular' ? 'NC_012920.1' : 'NC_001416.1'}\t${i * 100 + 1}\t${20 + i % 10}\n`).join(''));
      await (await reveal(page.getByLabel('Depth TSV', { exact: true }))).setInputFiles({ name: 'depth.tsv', mimeType: 'text/plain', buffer: depth });
      await check(page, 'Show Depth', true); await generate(page);
      for (const enabled of [false, true]) {
        await check(page, 'Show Depth', enabled);
        await agree(page, info, `depth-${enabled}`);
        const displays = await page.locator('.origin-top svg g[id="depth"]').evaluateAll(groups => groups.map(group => group.getAttribute('display')));
        expect(displays.length).toBeGreaterThan(0);
        expect(displays.every(display => enabled ? display !== 'none' : display === 'none')).toBe(true);
      }
    } finally { await page.context().close(); }
  });
}

test('B-06 J36 rejected FASTA recovery and label regeneration preserve scale geometry', async ({ page }, info) => {
  test.setTimeout(300000);
  page.on('dialog', dialog => dialog.dismiss());
  await openApp(page);
  await switchMode(page, 'linear');
  await page.getByRole('radio', { name: 'GFF3 + FASTA', exact: true }).check();
  await page.getByLabel('GFF3', { exact: true }).setInputFiles('gbdraw/web/tutorial-data/lambda-gff3/NC_001416.gff3');
  const fasta = page.getByLabel('FASTA', { exact: true });
  const valid = 'gbdraw/web/tutorial-data/lambda-gff3/NC_001416.fna';
  await fasta.setInputFiles(valid); await generate(page);
  await fasta.setInputFiles({ name: 'mismatch.fasta', mimeType: 'text/plain', buffer: Buffer.from('>wrong_record_identity\n' + 'ATGC'.repeat(100) + '\n') });
  expect(await page.evaluate(() => window.__GBDRAW_APP__.runAnalysis())).toEqual({ status: 'error' });
  await agree(page, info, 'rejected-input');
  await fasta.setInputFiles(valid);
  await popup(page); await label(page, 'JOURNEY_J36');
  await agree(page, info, 'restored-label');
  await closeEditor(page); await generate(page);
  expect(await page.locator('.origin-top svg').textContent()).toContain('JOURNEY_J36');
  await agree(page, info, 'regenerated-scale');
});
