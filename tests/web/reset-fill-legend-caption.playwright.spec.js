// OV-3350 (Owner D-15 = B, 2026-11-11): after Reset fill color, a feature type
// that no other rule colors keeps its own Legend caption ("rRNA"), live and
// after Generate; while another rule still colors a feature of the type, the
// remainder row stays "other rRNAs". The HmmtDNA Session draws two rRNAs,
// s-rRNA and l-rRNA, in the palette color of rRNA.
const { test, expect } = require('@playwright/test');
const { evaluateWithRetainedPromise, generateAndWaitForResult } = require('./helpers/app-lifecycle.cjs');
const { HMMT_SESSION, settle, legendCaptions } = require('./helpers/audit-browser.cjs');
const { load } = require('./helpers/mode-transition.cjs');

/** @param {import('@playwright/test').Page} page @param {{ val: string, color: string, cap: string }} rule */
const addRrnaRule = (page, rule) => evaluateWithRetainedPromise(page, async (row) => {
  const app = window.__GBDRAW_APP__;
  Object.assign(app.newSpecRule, { feat: 'rRNA', qual: 'product', ...row });
  await app.addSpecificRule();
}, rule);

/** @param {import('@playwright/test').Page} page @param {string} product */
const colorThisFeature = (page, product) => evaluateWithRetainedPromise(page, async (name) => {
  const app = window.__GBDRAW_APP__;
  return app.setFeatureColorValue(app.extractedFeatures.find((item) => item.product === name), '#ff00ff');
}, product);

/** @param {import('@playwright/test').Page} page @param {string} product */
const resetFill = (page, product) => evaluateWithRetainedPromise(page, async (name) => {
  const app = window.__GBDRAW_APP__;
  await app.openFeatureEditorFromList(app.extractedFeatures.find((item) => item.product === name), null);
  await app.resetClickedFeatureFillColor();
}, product);

/** The fill of each rRNA on the mounted Result, by product. */
const rrnaFills = (page) => page.evaluate(() => {
  const app = window.__GBDRAW_APP__;
  const svg = app.svgContainer.querySelector('svg');
  return Object.fromEntries(app.extractedFeatures.filter((feature) => feature.type === 'rRNA').map((feature) => [
    feature.product,
    [...svg.querySelectorAll(`[data-gbdraw-feature-id="${CSS.escape(feature.svg_id)}"]`)]
      .map((element) => element.getAttribute('fill')).find((value) => value && value !== 'none')
  ]));
});

/** The rRNA rows of the Legend: those whose caption names rRNA. */
const rrnaRows = async (page, options) => (await legendCaptions(page, options)).filter((caption) => /rRNA/.test(caption));

const CASES = [
  {
    // The audit's scenario: a regex rule colors s-rRNA, Reset fill takes it back to the palette color.
    name: 'a regex rule that colors only the reset feature leaves the plain caption',
    setup: (page) => addRrnaRule(page, { val: '^s-rRNA$', color: '#ff00ff', cap: 'small rRNA' }),
    edited: ['small rRNA', 'other rRNAs'],
    reset: ['rRNA']
  },
  {
    // A This feature only color, then Reset fill: no rule of rRNA remains.
    name: 'a This feature only color taken back leaves the plain caption',
    setup: (page) => colorThisFeature(page, 's-rRNA'),
    edited: ['s-rRNA', 'other rRNAs'],
    reset: ['rRNA']
  },
  {
    // Another rule still colors l-rRNA: the reset feature joins the remainder row.
    name: 'a rule that still colors another feature of the type keeps the remainder row',
    setup: async (page) => {
      await addRrnaRule(page, { val: '^s-rRNA$', color: '#ff00ff', cap: 'small rRNA' });
      await settle(page);
      await addRrnaRule(page, { val: '^l-rRNA$', color: '#0055ff', cap: 'large rRNA' });
    },
    edited: ['small rRNA', 'large rRNA'],
    reset: ['other rRNAs', 'large rRNA']
  }
];

for (const CASE of CASES) {
  test(`Reset fill color: ${CASE.name} (OV-3350)`, async ({ browser }) => {
    test.setTimeout(240_000);
    const page = await load(browser, HMMT_SESSION);
    await settle(page);
    const palette = await page.evaluate(() => window.__GBDRAW_APP__.appliedPaletteColors.rRNA);
    expect(await rrnaRows(page)).toEqual(['rRNA']);
    await CASE.setup(page);
    await settle(page);
    expect((await rrnaRows(page)).sort()).toEqual([...CASE.edited].sort());
    expect((await rrnaFills(page))['s-rRNA']).toBe('#ff00ff');

    await resetFill(page, 's-rRNA');
    await settle(page);
    await expect.poll(async () => (await rrnaFills(page))['s-rRNA'], { timeout: 10_000 }).toBe(palette);
    const live = await rrnaRows(page);
    expect([...live].sort()).toEqual([...CASE.reset].sort());

    await page.evaluate(() => { window.__GBDRAW_APP__.clickedFeature = null; });
    await generateAndWaitForResult(page);
    await settle(page);
    // Live and Generate draw the same rows in the same order.
    expect(await rrnaRows(page, { source: 'result' })).toEqual(live);
    expect((await rrnaFills(page))['s-rRNA']).toBe(palette);
    await page.context().close();
  });
}
