// Figure size (Layout, issue #1427): a requested width and height apply on Generate, the line under the field shows
// the drawn size in the requested unit, the plain SVG states the size as the CLI writes it while the page keeps px,
// a Legend rename keeps the requested canvas, Save and Load keep the size, Undo restores Auto and undoes a unit
// change, and a Linear width the parts that do not shrink cannot fit is refused with the parts named and the previous
// figure kept.
const { test, expect } = require('@playwright/test');
const { readFileSync } = require('node:fs');
const { gunzipSync } = require('node:zlib');
const path = require('node:path');
const { evaluateWithRetainedPromise, generateAndWaitForResult } = require('./helpers/app-lifecycle.cjs');
const { showDrawing } = require('./helpers/drawing-bar.cjs');
const { HMMT, loadSessionFile, openFresh, openWithGenBank, settle } = require('./helpers/audit-browser.cjs');
const { download } = require('./helpers/mode-transition.cjs');

const root = path.resolve(process.env.GBDRAW_REPO || process.cwd());
const MM_180_PX = 180 * 96 / 25.4;

const openLayout = async (page) => {
  const details = page.locator('[data-setting="layout"] > details');
  // Opened directly: a click on the summary's middle lands on its help tip.
  await details.evaluate((element) => { element.open = true; });
  await page.locator('#figure-width').scrollIntoViewIfNeeded();
};
const setSide = async (page, side, amount, unit) => {
  await page.locator(`#figure-${side}-unit`).selectOption(unit);
  await page.locator(`#figure-${side}`).fill(String(amount));
  await page.locator(`#figure-${side}`).dispatchEvent('change');
  await settle(page);
};
const sizeState = (page) => page.evaluate(() => {
  const { form } = window.__GBDRAW_APP__;
  const svg = document.querySelector('.gbdraw-preview-surface svg');
  const text = (selector) => document.querySelector(selector)?.textContent?.trim() ?? null;
  return {
    width: form.figure_width,
    height: form.figure_height,
    amount: document.querySelector('#figure-width')?.value ?? null,
    unit: document.querySelector('#figure-width-unit')?.value ?? null,
    line: text('[data-figure-size-line]'),
    larger: text('[data-figure-size-larger]'),
    root: svg ? { width: svg.getAttribute('width'), height: svg.getAttribute('height'), viewBox: svg.getAttribute('viewBox') } : null
  };
});
const undo = async (page) => { await evaluateWithRetainedPromise(page, () => window.__GBDRAW_HISTORY__.undo()); await settle(page); };
const px = (value) => Number.parseFloat(String(value));

test('a Circular 180 x 180 mm figure: Undo, size line, SVG export units, Legend rename, Save and Load', async ({ page }, testInfo) => {
  test.setTimeout(420_000);
  await openWithGenBank(page, path.join(root, HMMT));
  await openLayout(page);
  // One edit is one History step: one Undo restores Auto.
  await setSide(page, 'width', 180, 'mm');
  expect((await sizeState(page)).width).toBe('180mm');
  await undo(page);
  expect(await sizeState(page)).toMatchObject({ width: null, amount: '' });
  await setSide(page, 'width', 180, 'mm');
  await setSide(page, 'height', 180, 'mm');
  expect(await sizeState(page)).toMatchObject({ width: '180mm', height: '180mm', amount: '180', unit: 'mm' });
  // A unit change with an amount set is one step too: 180 mm becomes 180 cm, and one Undo gives 180 mm back.
  await page.locator('#figure-width-unit').selectOption('cm');
  await settle(page);
  expect(await sizeState(page)).toMatchObject({ width: '180cm', amount: '180', unit: 'cm' });
  await undo(page);
  expect(await sizeState(page)).toMatchObject({ width: '180mm', height: '180mm', amount: '180', unit: 'mm' });
  await generateAndWaitForResult(page);
  await settle(page);
  const drawn = await sizeState(page);
  expect(drawn.line).toBe('Figure: 180 × 180 mm (680 × 680 px)');
  expect(drawn.larger).toBe(null);
  // The page's root is px (OV-2512); the viewBox is the requested box.
  expect(px(drawn.root.width)).toBeCloseTo(MM_180_PX, 2);
  expect(px(drawn.root.height)).toBeCloseTo(MM_180_PX, 2);
  expect(drawn.root.width).toMatch(/px$/);

  // The plain SVG states the size in mm, as the CLI writes it; the viewBox stays px.
  const svgFile = testInfo.outputPath('sized.svg');
  const exported = (await download(page, 'SVG', svgFile)).toString('utf8');
  const head = exported.slice(exported.indexOf('<svg'), exported.indexOf('>', exported.indexOf('<svg')));
  expect(head).toContain('width="180mm"');
  expect(head).toContain('height="180mm"');
  expect(head).toMatch(/viewBox="0(\.0)? 0(\.0)? 680\.31/);
  expect((await sizeState(page)).root.width, 'the page root stays px').toBe(drawn.root.width);
  // The PDF page has the requested physical size as it is (PD-OI-072: 0.75 pt per px), within 0.1 mm.
  const pdf = (await download(page, 'PDF', testInfo.outputPath('sized.pdf'))).toString('latin1');
  const boxes = [...pdf.matchAll(/MediaBox\s*\[([^\]]*)\]/g)];
  expect(boxes).toHaveLength(1);
  const [, , boxWidth, boxHeight] = boxes[0][1].trim().split(/\s+/).map(Number);
  const MM_PER_PT = 25.4 / 72;
  expect(Math.abs(boxWidth * MM_PER_PT - 180)).toBeLessThan(0.1);
  expect(Math.abs(boxHeight * MM_PER_PT - 180)).toBeLessThan(0.1);
  console.log(`figure-size: PDF page ${boxWidth * MM_PER_PT} x ${boxHeight * MM_PER_PT} mm`);

  // A Legend rename replans in the requested canvas (section 6): the box stays, the content centred in it.
  await evaluateWithRetainedPromise(page, async () => {
    const app = window.__GBDRAW_APP__;
    await app.renameLegendEntry(app.legendEntries.findIndex((entry) => entry.caption === 'CDS'), 'Genes');
  });
  await settle(page);
  const renamed = await sizeState(page);
  expect(px(renamed.root.width)).toBeCloseTo(MM_180_PX, 2);
  expect(px(renamed.root.height)).toBeCloseTo(MM_180_PX, 2);
  expect(renamed.line).toBe('Figure: 180 × 180 mm (680 × 680 px)');

  // Save and Load keep the size in the field (the Session is the positive fixture
  // tests/fixtures/sessions/HmmtDNA-figure-size-180mm.v51-web.gbdraw-session.json.gz).
  const saved = testInfo.outputPath('figure-size.gbdraw-session.json.gz');
  const session = JSON.parse(gunzipSync(await download(page, 'Save Session', saved, { title: '' })).toString('utf8'));
  expect(session.drawings[0].config.form).toMatchObject({ figure_width: '180mm', figure_height: '180mm' });
  expect(session.drawings[0].renderRequest.diagramOptions.configOverrides)
    .toMatchObject({ 'canvas.figure_width': '180mm', 'canvas.figure_height': '180mm' });
  await setSide(page, 'width', 90, 'cm');
  await loadSessionFile(page, saved);
  await openLayout(page);
  expect(await sizeState(page)).toMatchObject({ width: '180mm', height: '180mm', amount: '180', unit: 'mm' });
});

test('a Linear 85 mm width with a wide side Legend is refused, naming the parts, and the previous figure is kept', async ({ page }) => {
  test.setTimeout(300_000);
  await openFresh(page);
  await showDrawing(page, { kind: 'linear' });
  await page.evaluate(async (text) => {
    const app = window.__GBDRAW_APP__;
    app.setLinearSeqPrimaryFile(0, 'gb', new File([text], 'HmmtDNA.gbk', { type: 'text/plain' }));
    app.form.legend = 'right';
    app.adv.legend_font_size = 30;
    await window.Vue.nextTick();
  }, readFileSync(path.join(root, HMMT), 'utf8'));
  await settle(page);
  await generateAndWaitForResult(page);
  const shown = await page.evaluate(() => window.__GBDRAW_APP__.currentResults()[0].content);
  await openLayout(page);
  await setSide(page, 'width', 85, 'mm');
  const outcome = await generateAndWaitForResult(page, { expectedStatus: 'error' });
  expect(outcome.errorSummary).toContain('The figure cannot be 85 mm wide: the parts that do not shrink need');
  expect(outcome.errorSummary).toContain('place the Legend at the top or bottom');
  expect(outcome.errorSummary, 'the sentence names Figure width; no Setting line').not.toContain('Setting:');
  expect(await page.evaluate(() => window.__GBDRAW_APP__.currentResults()[0].content), 'the figure is kept').toBe(shown);
});

// PD-OI-114 rev 2 (Q-D1 A): on a figure drawn to a requested size, a Legend rename and a Canvas padding change show
// at once, the line under the field says the figure is being fitted, and the automatic rerender then draws it again
// (Python fits it with the edit); the edit and its rerender are one History step, and its Undo fits the figure again.
// A rename the box cannot hold is refused through the automatic rerender's failure note. On Auto nothing is redrawn.
test('a Legend rename and Canvas padding on a sized figure are fitted again by the automatic rerender; on Auto they are not', async ({ page }) => {
  test.setTimeout(420_000);
  await openWithGenBank(page, path.join(root, HMMT));
  await page.evaluate(() => {
    const hooks = window.__GBDRAW_TEST_HOOKS__ || {};
    const render = hooks.beforeDiagramGenerationResponse;
    window.__renders = 0;
    window.__holdRender = null;
    window.__GBDRAW_TEST_HOOKS__ = {
      ...hooks,
      beforeDiagramGenerationResponse: async (...args) => {
        window.__renders += 1;
        if (window.__holdRender) await window.__holdRender;
        return render?.(...args);
      }
    };
  });
  const renders = () => page.evaluate(() => window.__renders);
  const countFrom = () => page.evaluate(() => { window.__renders = 0; });
  const rename = (from, to) => evaluateWithRetainedPromise(page, async (names) => {
    const app = window.__GBDRAW_APP__;
    await app.renameLegendEntry(app.legendEntries.findIndex((entry) => entry.caption === names.from), names.to);
  }, { from, to }).then(() => settle(page));
  const fittingLine = () => page.evaluate(() => document.querySelector('[data-figure-size-fitting]')?.textContent?.trim() ?? null);

  await generateAndWaitForResult(page);
  await settle(page);
  await countFrom();
  await rename('CDS', 'Genes');
  expect(await page.evaluate(() => window.__GBDRAW_APP__.legendEntries.some((entry) => entry.caption === 'Genes'))).toBe(true);
  await undo(page);
  expect(await page.evaluate(() => window.__GBDRAW_APP__.legendEntries.some((entry) => entry.caption === 'CDS'))).toBe(true);
  expect(await renders(), 'an Auto figure is not redrawn').toBe(0);

  await openLayout(page);
  await setSide(page, 'width', 160, 'mm');
  await generateAndWaitForResult(page);
  await settle(page);
  const MM_160_PX = 160 * 96 / 25.4;
  expect(px((await sizeState(page)).root.width)).toBeCloseTo(MM_160_PX, 2);

  const undoCount = () => page.evaluate(() => window.__GBDRAW_HISTORY__.getUndoCount());
  const legendTexts = () => page.evaluate(() => [...document.querySelectorAll('.gbdraw-preview-surface svg #legend text')]
    .map((text) => text.textContent));
  const feedback = () => page.evaluate(() => document.querySelector('[data-live-application-feedback]')?.textContent?.trim() ?? null);

  // The rerender is held: the edit shows at once and the line says the figure is being fitted.
  const steps = await undoCount();
  await page.evaluate(() => { window.__renders = 0; window.__holdRender = new Promise((resolve) => { window.__releaseRender = resolve; }); });
  await evaluateWithRetainedPromise(page, async () => {
    const app = window.__GBDRAW_APP__;
    await app.renameLegendEntry(app.legendEntries.findIndex((entry) => entry.caption === 'CDS'), 'Genes');
  });
  await expect.poll(fittingLine).toBe('Fitting the figure to the requested size.');
  expect(await legendTexts(), 'the rename shows at once').toContain('Genes');
  expect(await page.evaluate(() => document.querySelector('[data-figure-size-larger]'))).toBe(null);
  await page.evaluate(() => { window.__holdRender = null; window.__releaseRender(); });
  await settle(page);
  await expect.poll(renders).toBe(1);
  await settle(page);
  const fitted = await sizeState(page);
  expect(await fittingLine()).toBe(null);
  expect(fitted.larger).toBe(null);
  expect(px(fitted.root.width)).toBeCloseTo(MM_160_PX, 2);
  expect(await legendTexts()).toContain('Genes');
  expect(await undoCount(), 'the rename and its rerender are one step').toBe(steps + 1);

  // One Undo reverts the rename with its rerender: History restores the fitted figure before it, so Python does not
  // run. Redo shows the fitted figure the rerender drew, again without a Python run (D-24).
  await countFrom();
  await undo(page);
  await settle(page);
  expect(px((await sizeState(page)).root.width)).toBeCloseTo(MM_160_PX, 2);
  expect(await page.evaluate(() => window.__GBDRAW_APP__.legendEntries.some((entry) => entry.caption === 'CDS'))).toBe(true);
  expect(await legendTexts()).not.toContain('Genes');
  expect(await undoCount()).toBe(steps);
  await evaluateWithRetainedPromise(page, () => window.__GBDRAW_HISTORY__.redo());
  await settle(page);
  const redone = await sizeState(page);
  expect(await legendTexts()).toContain('Genes');
  expect(px(redone.root.width)).toBeCloseTo(MM_160_PX, 2);
  expect(redone.larger).toBe(null);
  expect(await renders(), 'Undo and Redo of a fitted step run no Python').toBe(0);
  await undo(page);

  // A rename that makes the figure unable to fit: Python refuses the rerender, which reports it through the automatic
  // rerender's failure note; the edit stays as the Web drew it, and the size line says why the figure is larger.
  // One Undo reverts it and shows the fitted figure before it.
  const LONG = 'Protein-coding genes of the mitochondrion';
  await countFrom();
  await rename('CDS', LONG);
  await expect.poll(renders).toBe(1);
  await settle(page);
  const refused = await sizeState(page);
  const note = await feedback();
  console.log(`figure-size refit refused: note=${JSON.stringify(note)} line=${JSON.stringify(refused.line)} larger=${JSON.stringify(refused.larger)}`);
  expect(note).toContain('Live edit failed');
  expect(note).toContain('The figure cannot be 160 mm wide');
  expect(await legendTexts(), 'the edit is kept as the Web drew it').toContain(LONG);
  expect(px(refused.root.width)).toBeGreaterThan(MM_160_PX);
  expect(refused.larger).toBe('Larger than the requested 160 mm width: the edit could not be fitted (see the note above the preview).');
  expect(await undoCount()).toBe(steps + 1);
  await countFrom();
  await undo(page);
  await settle(page);
  const restored = await sizeState(page);
  console.log(`figure-size refused refit undone: note=${JSON.stringify(await feedback())} larger=${JSON.stringify(restored.larger)}`);
  expect(px(restored.root.width)).toBeCloseTo(MM_160_PX, 2);
  expect(restored.larger).toBe(null);
  expect(await legendTexts()).not.toContain(LONG);
  expect(await renders(), 'the figure before the rename is drawn: no Python run').toBe(0);

  // A Canvas padding change is fitted again too: the padding is inside the requested width.
  await countFrom();
  await page.evaluate(() => { window.__GBDRAW_APP__.showCanvasControls = true; });
  const left = page.locator('#canvas-padding-left');
  await left.fill('40');
  await left.press('Tab');
  await expect.poll(renders).toBe(1);
  await settle(page);
  const padded = await sizeState(page);
  expect(px(padded.root.width)).toBeCloseTo(MM_160_PX, 2);
  expect(padded.larger).toBe(null);
});
