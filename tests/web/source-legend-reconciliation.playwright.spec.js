const { test, expect } = require('@playwright/test');
const fs = require('node:fs/promises');
const { load, generate, switchMode, download } = require('./helpers/mode-transition.cjs');

const seeds = {
  lambda: 'gbdraw/web/gallery/sessions/lambda_basic_linear.gbdraw-session.json',
  tobacco: 'gbdraw/web/gallery/sessions/tobacco-chloroplast.gbdraw-session.json'
};
const upload = async (page, name, linear = true) => {
  const session = JSON.parse(await fs.readFile(seeds[name], 'utf8'));
  const input = linear ? page.getByTestId('linear-genbank-1') : page.getByLabel('GenBank/DDBJ File', { exact: true });
  await input.setInputFiles({ name: `${name}.gb`, mimeType: 'text/plain', buffer: Buffer.from(
    session.resources[session.renderRequest.records[0].source.resourceId].data, 'base64') });
};
const inspect = page => page.evaluate(async () => {
  const { state: s } = await import('./js/state.js');
  const { getVisibleFeatureLegendGroup } = await import('./js/app/legend/utils.js');
  const captions = svg => [...(getVisibleFeatureLegendGroup(svg)?.querySelectorAll('g[data-legend-key]') || [])]
    .map(e => ({ caption: e.getAttribute('data-legend-key'), color: e.querySelector('path[fill]')?.getAttribute('fill') }))
    .sort((a, b) => a.caption.localeCompare(b.caption));
  const result = s.results.value[s.selectedResultIndex.value].content;
  const svg = new DOMParser().parseFromString(result, 'image/svg+xml').documentElement;
  return JSON.parse(JSON.stringify({
    entries: s.legendEntries.value.map(e => e.caption).sort(), original: s.originalLegendOrder.value,
    result: captions(svg), mounted: captions(s.svgContainer.value.querySelector('svg')),
    side: s.form.legend, fontSize: s.adv.legend_font_size, preferences: s.layoutPreferences.legend,
    palette: s.selectedPalette.value, colors: s.currentColors.value,
    featureIds: s.extractedFeatures.value.map(f => f.stable_feature_id),
    overrides: s.legendColorOverrides, deleted: s.deletedLegendEntries.value
  }));
});
const expectEntries = async (page, expected) => {
  const current = await inspect(page);
  expect(current.entries).toEqual([...expected].sort());
  expect(current.result.map(e => e.caption).sort()).toEqual([...expected].sort());
  expect(current.mounted).toEqual(current.result);
  return current;
};
const reveal = async locator => {
  for (const details of await locator.locator('xpath=ancestor::details').all()) {
    if (await details.getAttribute('open') === null) await details.locator(':scope > summary').click();
  }
  return locator;
};

test('L1-L8 generated legend categories reconcile while valid category and layout preferences survive', async ({ browser }, testInfo) => {
  test.setTimeout(1_800_000);
  const page = await load(browser, seeds.lambda);
  page.setDefaultTimeout(180_000);
  let fresh;
  try {
    await generate(page);
    const a = await expectEntries(page, ['CDS']);
    const palette = page.getByLabel('Palette', { exact: true });
    const alternatives = await palette.locator('option').evaluateAll(options => options.map(o => o.value));
    await (await reveal(palette)).selectOption(alternatives.find(value => value !== a.palette));
    await (await reveal(page.locator('#legend-position'))).selectOption('left');
    await page.evaluate(() => {
      const app = window.__GBDRAW_APP__;
      app.adv.legend_font_size = 18;
      app.updateLegendEntryColor(app.legendEntries.findIndex(e => e.caption === 'CDS'), '#123456');
    });
    await generate(page);
    const edited = await inspect(page);
    await upload(page, 'tobacco');
    await generate(page);
    const b = await expectEntries(page, ['CDS', 'repeat_region', 'rRNA', 'tRNA']);
    expect(b.featureIds.some(id => a.featureIds.includes(id))).toBe(false);
    expect(b.result.find(e => e.caption === 'CDS').color).toBe('#123456');
    for (const key of ['side', 'fontSize', 'preferences', 'palette', 'colors']) expect(b[key]).toEqual(edited[key]);
    await switchMode(page, 'circular');
    await switchMode(page, 'linear');
    expect((await inspect(page)).overrides).toEqual(b.overrides);
    await generate(page);
    await upload(page, 'lambda');
    await generate(page);
    const returned = await expectEntries(page, ['CDS']);
    for (const [action, captions] of [['Undo', b.entries], ['Redo', returned.entries]]) {
      await page.getByRole('button', { name: action, exact: true }).click();
      await expect.poll(() => page.evaluate(() => !window.__GBDRAW_HISTORY__.restoring.value
        && !window.__GBDRAW_HISTORY__.capturing.value)).toBe(true);
      await expectEntries(page, captions);
    }
    await generate(page);
    expect((await inspect(page)).result).toEqual(returned.result);
    const saved = testInfo.outputPath('legend-replacement.gbdraw-session.json.gz');
    await download(page, 'Save Session', saved);
    fresh = await load(browser, saved);
    await expectEntries(fresh, ['CDS']);
    await generate(fresh);
    expect((await expectEntries(fresh, ['CDS'])).result).toEqual(returned.result);
    expect(page.externalRequests).toEqual([]);
    expect(fresh.externalRequests).toEqual([]);
  } finally {
    await page.context().close();
    if (fresh) await fresh.context().close();
  }
});

test('manual legend rows survive replacement and absent customized categories do not reject a valid source', async ({ browser }, testInfo) => {
  test.setTimeout(1_800_000);
  const page = await load(browser);
  page.setDefaultTimeout(180_000);
  let fresh;
  try {
    await upload(page, 'tobacco', false);
    await generate(page);
    await page.evaluate(async () => {
      const app = window.__GBDRAW_APP__;
      app.newLegendCaption = 'Retained annotation';
      app.newLegendColor = '#884422';
      await app.addNewLegendEntry();
    });
    await expect.poll(() => page.evaluate(() => window.__GBDRAW_APP__.legendEntries.some(e => e.caption === 'Retained annotation'))).toBe(true);
    await page.evaluate(async () => {
      const app = window.__GBDRAW_APP__;
      app.updateLegendEntryColor(app.legendEntries.findIndex(e => e.caption === 'rRNA'), '#abcdef');
      await app.renameLegendEntry(app.legendEntries.findIndex(e => e.caption === 'rRNA'), 'Ribosomal RNA');
      app.deleteLegendEntry(app.legendEntries.findIndex(e => e.caption === 'tRNA'));
    });
    await generate(page);
    const a = await inspect(page);
    expect(a.entries).toContain('Retained annotation');
    expect(a.original).not.toContain('Retained annotation');
    const input = page.getByLabel('GenBank/DDBJ File', { exact: true });
    await input.setInputFiles({ name: 'corrupt.gb', mimeType: 'text/plain', buffer: Buffer.from('not a GenBank record') });
    expect(await page.evaluate(() => window.__GBDRAW_APP__.runAnalysis())).toEqual({ status: 'error' });
    expect(await inspect(page)).toEqual(a);
    await upload(page, 'lambda', false);
    await generate(page);
    const expected = ['CDS', 'GC content', 'GC skew (+)', 'GC skew (-)', 'Retained annotation'];
    await expectEntries(page, expected);
    await generate(page);
    await expectEntries(page, expected);
    const saved = testInfo.outputPath('manual-legend.gbdraw-session.json.gz');
    await download(page, 'Save Session', saved);
    fresh = await load(browser, saved);
    await generate(fresh);
    await expectEntries(fresh, expected);
    expect(page.externalRequests).toEqual([]);
    expect(fresh.externalRequests).toEqual([]);
  } finally {
    await page.context().close();
    if (fresh) await fresh.context().close();
  }
});
