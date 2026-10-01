const { test, expect } = require('@playwright/test');
const { readFileSync } = require('node:fs');
const { join, resolve } = require('node:path');
const { generateAndWaitForResult, openApp, reveal } = require('./helpers/app-lifecycle.cjs');

// The Web Circular defaults use the Multi-Record Canvas. Its output must equal
// the single-record Circular path: no empty depth slot without depth input, a
// legend from the same slot-aware owner, and a center definition that fits.
const repoRoot = resolve(process.env.GBDRAW_REPO || process.cwd());
const hmmt = readFileSync(join(repoRoot, 'tests/test_inputs/HmmtDNA.gbk'), 'utf8');
const withOrganism = (organism) => hmmt.replace('/organism="Homo sapiens"', `/organism="${organism}"`);

const uploadGenbank = async (page, name, text) => {
  const chooser = page.waitForEvent('filechooser');
  await page.getByRole('button', { name: 'Choose GenBank/DDBJ File', exact: true }).click();
  await (await chooser).setFiles({ name, mimeType: 'text/plain', buffer: Buffer.from(text) });
  await expect.poll(
    () => page.evaluate(() => window.__GBDRAW_APP__.circularRecordList.length),
    { timeout: 60000 }
  ).toBeGreaterThan(0);
};

const inspectResult = (page) => page.evaluate(() => {
  const app = window.__GBDRAW_APP__;
  const doc = new DOMParser().parseFromString(String(app.results?.[0]?.content || ''), 'image/svg+xml');
  const definition = doc.querySelector('[data-gbdraw-role="record-definition"]');
  return {
    multiRecordCanvas: app.form.multi_record_canvas,
    slotIds: [...doc.querySelectorAll('[data-gbdraw-slot-id]')].map((node) => node.getAttribute('data-gbdraw-slot-id')),
    definitionLines: definition ? [...definition.querySelectorAll('text')].map((node) => node.textContent.trim()) : [],
    legendKeys: [...doc.querySelectorAll('#legend [data-legend-key]')].map((node) => node.getAttribute('data-legend-key'))
  };
});

for (const organism of [
  'Salmonella enterica subsp. enterica serovar Typhimurium',
  'Mycobacterium tuberculosis variant bovis BCG'
]) {
  test(`Web Circular defaults Generate ${organism} without an empty depth slot`, async ({ page }) => {
    test.setTimeout(240000);
    await openApp(page);
    await uploadGenbank(page, 'long-organism.gbk', withOrganism(organism));

    const outcome = await generateAndWaitForResult(page);
    const result = await inspectResult(page);

    expect(outcome.errorSummary).toBe('');
    expect(result.multiRecordCanvas).toBe(true);
    expect(result.slotIds).not.toContain('depth');
    expect(result.slotIds).toEqual(expect.arrayContaining(['gc_content', 'gc_skew']));
    expect(result.definitionLines.join(' ')).toContain(organism);
  });
}

test('Web Circular default legend uses custom slot labels and added skew slots', async ({ page }) => {
  test.setTimeout(240000);
  await openApp(page);
  await uploadGenbank(page, 'HmmtDNA.gbk', hmmt);

  const panel = page.locator('button[aria-controls="circular-custom-track-slots-panel"]');
  await reveal(panel);
  if (await panel.getAttribute('aria-expanded') !== 'true') await panel.click();
  await page.getByText('Use custom stack', { exact: true }).locator('input').check();
  const label = page.getByRole('group', { name: 'Circular track slot gc_content', exact: true })
    .getByRole('textbox', { name: 'Track legend label' });
  await label.fill('MY GC');
  await label.press('Tab');
  await page.getByRole('combobox', { name: 'New circular track renderer', exact: true }).selectOption('dinucleotide_skew');
  await page.getByRole('button', { name: /Add track/ }).click();
  const newId = await page.evaluate(() => window.__GBDRAW_APP__.adv.circular_track_slots.at(-1).id);
  const nucleotide = page.getByRole('group', { name: `Circular track slot ${newId}`, exact: true })
    .getByRole('textbox', { name: 'Track dinucleotide' });
  await nucleotide.fill('AT');
  await nucleotide.press('Tab');

  await generateAndWaitForResult(page);
  const result = await inspectResult(page);

  expect(result.multiRecordCanvas).toBe(true);
  expect(result.legendKeys).toContain('MY GC');
  expect(result.legendKeys).toEqual(expect.arrayContaining(['AT skew (+)', 'AT skew (-)']));
  expect(result.legendKeys).not.toContain('GC content');
});
