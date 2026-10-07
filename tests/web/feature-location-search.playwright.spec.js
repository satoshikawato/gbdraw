const fs = require('node:fs');
const { execFileSync } = require('node:child_process');
const { join } = require('node:path');
const { pathToFileURL } = require('node:url');
const { test, expect } = require('@playwright/test');
const {
  evaluateWithRetainedPromise,
  generateAndWaitForResult,
  openApp,
  waitForAppShell
} = require('./helpers/app-lifecycle.cjs');

const FIXTURE = join(process.cwd(), 'tests/fixtures/feature_location_search.gb');
const SPLIT = '1101..1200, 1..200 (+)';
const MINUS = '901..1000, 701..800 (-)';

// A point inside the feature's filled shape, in viewport coordinates.
const featurePoint = (page, root, featureId) => page.evaluate(({ root, featureId }) => {
  const element = [...document.querySelectorAll(`${root} [data-gbdraw-feature-id="${CSS.escape(featureId)}"]`)]
    .find((node) => node.localName === 'path' && node.getAttribute('fill') !== 'none');
  const box = element.getBBox();
  for (let fy = 0.05; fy < 1; fy += 0.05) {
    for (let fx = 0.05; fx < 1; fx += 0.05) {
      const point = new DOMPoint(box.x + box.width * fx, box.y + box.height * fy);
      if (!element.isPointInFill(point)) continue;
      const screen = point.matrixTransform(element.getScreenCTM());
      if (document.elementFromPoint(screen.x, screen.y) === element) return { x: screen.x, y: screen.y };
    }
  }
  return null;
}, { root, featureId });

test('feature list, popup, hover, and Location search use 1-based INSDC locations', async ({ page }) => {
  test.setTimeout(240000);
  page.on('dialog', (dialog) => dialog.dismiss());
  await openApp(page);
  await page.getByLabel('GenBank/DDBJ File', { exact: true }).setInputFiles(FIXTURE);
  await generateAndWaitForResult(page);
  const ids = await page.evaluate(() => Object.fromEntries(window.__GBDRAW_APP__.extractedFeatures
    .filter((feature) => feature.locus_tag)
    .map((feature) => [feature.locus_tag, feature.svg_id])));

  const drawer = page.locator('.right-drawer');
  await page.locator('.drawer-toggle').click();
  await expect(drawer.getByText(/^Features \(\d+\)$/)).toBeVisible();
  for (const location of ['301..600 (+)', SPLIT, MINUS]) {
    await expect(drawer.getByTitle(location, { exact: true })).toHaveText(location);
  }
  await evaluateWithRetainedPromise(page, async (id) => {
    const app = window.__GBDRAW_APP__;
    await app.openFeatureEditorFromList(app.extractedFeatures.find((feature) => feature.svg_id === id), null);
  }, ids.LOC_0002);
  await expect(page.getByRole('dialog', { name: /^Feature details:/ }))
    .toContainText(`LOCTEST: ${SPLIT}`);
  expect(await page.evaluate(() => window.__GBDRAW_APP__.clickedFeature.detailRows
    .find((row) => row.key === 'location').value)).toBe(SPLIT);
  await page.getByRole('button', { name: 'Close feature popup', exact: true }).click();

  const point = await featurePoint(page, '.origin-top svg', ids.LOC_0002);
  expect(point).not.toBeNull();
  await page.mouse.move(point.x, point.y);
  const summary = page.getByRole('tooltip');
  await expect(summary).toContainText(SPLIT);
  await expect(summary).toContainText('300 bp');
  await page.mouse.move(0, 0);

  const search = page.getByRole('searchbox', { name: 'Search features', exact: true });
  const status = page.getByRole('status', { name: 'Feature search status' });
  const find = async (field, query) => {
    await page.getByLabel('Search field', { exact: true }).selectOption(field);
    await search.fill(query);
    await search.press('Enter');
  };
  await find('location', '300');
  await expect(status).toHaveText('0 / 0 features');
  await find('location', '1..200');
  await expect(status).toHaveText('1 / 1 features');
  await find('all', 'CYTB');
  await expect(status).toHaveText('1 / 1 features');
  expect(await page.evaluate(() => window.__GBDRAW_APP__.previewFeatureSearchMatches))
    .toEqual([ids.LOC_0003]);
  await find('nucleotide', 'CYTB');
  await expect.poll(() => page.evaluate(
    () => window.__GBDRAW_APP__.previewFeatureSearchMatches.length
  )).toBeGreaterThan(1);
});

test('Interactive SVG searches names in All and shows split locations and lengths', async ({ page }, testInfo) => {
  test.setTimeout(120000);
  const prefix = testInfo.outputPath('location-search');
  execFileSync('python', [
    '-m', 'gbdraw.cli', 'circular', '--gbk', FIXTURE, '-o', prefix, '-f', 'interactive_svg'
  ], { cwd: process.cwd(), stdio: 'ignore' });
  const svgPath = `${prefix}.interactive.svg`;
  const catalog = fs.readFileSync(svgPath, 'utf8');
  const featureId = (locusTag) => {
    const metadata = /<metadata id="gbdraw-interactive-feature-metadata"[^>]*>([\s\S]*?)<\/metadata>/.exec(catalog)[1]
      .replaceAll('&quot;', '"').replaceAll('&lt;', '<').replaceAll('&gt;', '>').replaceAll('&amp;', '&');
    const item = JSON.parse(metadata).items[0];
    const biological = item.biologicalFeatures.find((feature) => (
      [feature.locus_tag, ...(feature.qualifiers?.locus_tag || [])].includes(locusTag)
    ));
    return item.features.find((feature) => feature.biologicalFeatureId === biological.biologicalFeatureId).svgId;
  };
  await page.goto(pathToFileURL(svgPath).href);
  await page.getByRole('button', { name: 'Expand feature search' }).click();
  const count = page.locator('[data-search-count]');
  const find = async (field, query) => {
    await page.locator('[data-search-field]').selectOption(field);
    await page.locator('[data-search-query]').fill(query);
    await page.locator('[data-search-apply]').click();
  };
  await find('all', 'CYTB');
  await expect(count).toHaveText('1 / 1 features');
  await find('nucleotide', 'CYTB');
  await expect(count).not.toHaveText(/^\d+ \/ [01] features$/);
  await find('location', '300');
  await expect(count).toHaveText('0 / 0 features');
  await find('location', '1..200');
  await expect(count).toHaveText('1 / 1 features');
  await page.locator('[data-search-clear]').click();

  const point = await featurePoint(page, 'svg', featureId('LOC_0002'));
  expect(point).not.toBeNull();
  await page.mouse.move(point.x, point.y);
  const hover = page.locator('#gbdraw-feature-hover-popup');
  await expect(hover).toContainText(SPLIT);
  await expect(hover).toContainText('300 bp');
  await page.mouse.click(point.x, point.y);
  await expect(page.locator('#gbdraw-feature-popup .gfi-subtitle')).toHaveText(SPLIT);
});

test('Interactive SVG searches on Enter, fits the Qualifier key, and titles a feature as the app popup does', async ({ page }, testInfo) => {
  test.setTimeout(120000);
  const prefix = testInfo.outputPath('hmmt');
  execFileSync('python', [
    '-m', 'gbdraw.cli', 'circular', '--gbk', join(process.cwd(), 'tests/test_inputs/HmmtDNA.gbk'),
    '-o', prefix, '-f', 'interactive_svg'
  ], { cwd: process.cwd(), stdio: 'ignore' });
  const svgPath = `${prefix}.interactive.svg`;
  const metadata = /<metadata id="gbdraw-interactive-feature-metadata"[^>]*>([\s\S]*?)<\/metadata>/
    .exec(fs.readFileSync(svgPath, 'utf8'))[1]
    .replaceAll('&quot;', '"').replaceAll('&lt;', '<').replaceAll('&gt;', '>').replaceAll('&amp;', '&');
  const item = JSON.parse(metadata).items[0];
  const trnf = item.biologicalFeatures.find((feature) => (
    feature.type === 'tRNA' && (feature.qualifiers?.product || []).includes('tRNA-Phe')
  ));
  const trnfId = item.features.find((feature) => feature.biologicalFeatureId === trnf.biologicalFeatureId).svgId;
  await page.goto(pathToFileURL(svgPath).href);
  await page.getByRole('button', { name: 'Expand feature search' }).click();
  const query = page.locator('[data-search-query]');
  const count = page.locator('[data-search-count]');
  const popup = page.locator('#gbdraw-feature-popup');

  // Enter searches as the Search button does (FL-08); a second Enter opens the active match.
  await query.fill('CYTB');
  await query.press('Enter');
  await expect(count).toHaveText('1 / 1 features');
  await query.press('Enter');
  await expect(popup.locator('.gfi-title')).toHaveText('cytochrome b');
  await expect(popup.locator('.gfi-content')).toContainText('Protein IDYP_003024038.1');
  await popup.locator('[data-close]').click();
  // Enter on a focused button presses that button.
  await page.locator('[data-search-clear]').focus();
  await page.keyboard.press('Enter');
  // The bar's inputs are XHTML in an SVG document, which toHaveValue does not read.
  await expect.poll(() => query.evaluate((input) => input.value)).toBe('');
  await expect(count).toHaveText('0 / 37 features');
  await expect(popup).toBeHidden();

  // The Qualifier key input takes a row of the bar, and the bar shows all of its rows (FL-09).
  await page.locator('[data-search-field]').selectOption('qualifier-value');
  const qualifier = page.locator('[data-search-qualifier]');
  const fit = await qualifier.evaluate((input) => {
    const controls = input.closest('foreignObject');
    const bar = input.closest('.gfs');
    return {
      share: input.getBoundingClientRect().width / bar.getBoundingClientRect().width,
      overflow: bar.scrollHeight - Number(controls.getAttribute('height'))
    };
  });
  expect(fit.share).toBeGreaterThan(0.9);
  expect(fit.overflow).toBeLessThanOrEqual(0);
  await qualifier.fill('product');
  await query.fill('tRNA-Phe');
  await qualifier.press('Enter');
  await expect(count).toHaveText('1 / 1 features');
  await page.locator('[data-search-clear]').click();

  // The popup title follows the app popup: a tRNA is titled by its product, and
  // Details has no Protein ID row without a protein_id (GX-08).
  const point = await featurePoint(page, 'svg', trnfId);
  expect(point).not.toBeNull();
  await page.mouse.click(point.x, point.y);
  await expect(popup.locator('.gfi-title')).toHaveText('tRNA-Phe');
  await expect(popup.locator('.gfi-content')).toContainText('577..647 (+)');
  await expect(popup.locator('.gfi-content')).not.toContainText('Protein ID');
});

test('a none Specific Table color survives Generate, Save, and Load', async ({ page }, testInfo) => {
  test.setTimeout(360000);
  page.on('dialog', (dialog) => dialog.dismiss());
  await openApp(page);
  await page.getByLabel('GenBank/DDBJ File', { exact: true }).setInputFiles(FIXTURE);
  const table = testInfo.outputPath('none-color.tsv');
  fs.writeFileSync(table, 'CDS\tproduct\tduplicate\tNONE\thollow\n');
  await page.getByLabel('Specific Table (-t)', { exact: true }).setInputFiles(table);
  await page.evaluate(async () => {
    const app = window.__GBDRAW_APP__;
    await app.waitForAuxiliaryFileImport(app.files.t_color);
  });
  const rule = { feat: 'CDS', qual: 'product', val: 'duplicate', color: 'none', cap: 'hollow' };
  const rules = () => page.evaluate(() => window.__GBDRAW_APP__.manualSpecificRules
    .map(({ feat, qual, val, color, cap }) => ({ feat, qual, val, color, cap })));
  expect(await rules()).toEqual([rule]);
  expect(await page.evaluate(() => window.__GBDRAW_APP__.errorLog)).toBeNull();
  await generateAndWaitForResult(page);
  const fills = () => page.evaluate(() => {
    const app = window.__GBDRAW_APP__;
    const id = app.extractedFeatures.find((feature) => feature.locus_tag === 'LOC_0001').svg_id;
    return [...document.querySelectorAll(`.origin-top svg [data-gbdraw-feature-id="${CSS.escape(id)}"]`)]
      .map((element) => element.getAttribute('fill'))
      .filter((fill) => fill && fill !== 'none')
      .length;
  });
  expect(await fills()).toBe(0);
  await expect(page.locator('.origin-top svg')).toContainText('hollow');

  const pendingSave = page.waitForEvent('download');
  await evaluateWithRetainedPromise(page, async () => {
    window.__GBDRAW_APP__.sessionTitle = 'none-color';
    await window.__GBDRAW_APP__.saveSessionWithTitle();
  });
  const saved = await (await pendingSave).path();
  await page.reload({ waitUntil: 'domcontentloaded' });
  await waitForAppShell(page);
  await page.locator('input[accept^=".json,"]').first().setInputFiles(saved);
  await expect.poll(() => page.evaluate(() => window.__GBDRAW_APP__.sessionImportPending), {
    timeout: 180000
  }).toBe(false);
  expect(await page.evaluate(() => window.__GBDRAW_APP__.errorLog)).toBeNull();
  expect(await rules()).toEqual([rule]);
  await generateAndWaitForResult(page);
  expect(await page.evaluate(() => window.__GBDRAW_APP__.errorLog)).toBeNull();
  expect(await fills()).toBe(0);
});
