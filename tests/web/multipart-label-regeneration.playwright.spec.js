const { test, expect } = require('@playwright/test');
const fs = require('node:fs/promises');
const { openApp } = require('./helpers/app-lifecycle.cjs');

const seed = 'gbdraw/web/gallery/sessions/tobacco-chloroplast.gbdraw-session.json';
const cases = [
  ['M1 single-part label', 'rbcL', 'label'],
  ['M2 multipart label', 'rps12', 'label'],
  ['M3 multipart placement', 'rps12', 'placement'],
  ['M4 multipart placement and label', 'rps12', 'both'],
  ['M5 multipart saved label', 'rps12', 'save'],
  ['M6 multipart mode round trip', 'rps12', 'mode'],
  ['M7 multipart label Undo and Redo', 'rps12', 'history']
];

for (const [name, gene, action] of cases) {
  test(`${name} remains bound through Generate`, async ({ browser }, testInfo) => {
    test.setTimeout(300000);
    let context;
    let page;
    const external = [];
    const errors = [];
    const load = async (file) => {
      if (context) await context.close();
      context = await browser.newContext({ viewport: { width: 1600, height: 1000 }, acceptDownloads: true });
      await context.route('**/*', (route) => {
        if (new URL(route.request().url()).hostname === '127.0.0.1') return route.continue();
        external.push(route.request().url());
        return route.abort();
      });
      page = await context.newPage();
      page.on('pageerror', (error) => errors.push(error.message));
      page.on('dialog', (dialog) => dialog.type() === 'prompt' ? dialog.accept(name) : dialog.dismiss());
      await page.addInitScript(() => {
        const post = Worker.prototype.postMessage;
        Worker.prototype.postMessage = function (message, ...rest) {
          if (message?.type === 'run') window.__LABEL_TEST_REQUEST__ = message.payload.request;
          return post.call(this, message, ...rest);
        };
      });
      await openApp(page);
      await page.locator('input[accept^=".json,"]').setInputFiles(file);
      await expect.poll(() => page.evaluate(() => !window.__GBDRAW_APP__.sessionImportPending
        && window.__GBDRAW_APP__.extractedFeatures.length > 0), { timeout: 180000 }).toBe(true);
    };
    const generate = async () => {
      expect(await page.evaluate(() => window.__GBDRAW_APP__.runAnalysis())).toEqual({ status: 'ok' });
      await expect.poll(() => page.evaluate(async () => {
        const { state } = await import('./js/state.js');
        const { isCommittedSvgResultMounted } = await import('./js/services/svg-result-ingestion.js');
        return isCommittedSvgResultMounted(state.results.value[state.selectedResultIndex.value]);
      })).toBe(true);
      expect(await page.evaluate(async () => (await import('./js/services/config.js'))
        .getCommittedCanonicalRenderRequest())).toEqual(await page.evaluate(() => window.__LABEL_TEST_REQUEST__));
    };
    const target = () => page.evaluate((gene) => {
      const feature = window.__GBDRAW_APP__.extractedFeatures.find((f) => f.type === 'CDS'
        && f.qualifiers.gene?.[0] === gene
        && (gene !== 'rps12' || f.qualifiers.locus_tag?.[0] === 'NitaCp049'));
      if (!feature) throw new Error('Retained biological target missing');
      return JSON.parse(JSON.stringify(feature));
    }, gene);
    const overrides = () => page.evaluate(async () => {
      const { state } = await import('./js/state.js');
      return JSON.parse(JSON.stringify({ labels: state.labelTextFeatureOverrides,
        placements: state.featurePlacementOverrides }));
    });
    const semantic = (content) => page.evaluate((content) => {
      const svg = new DOMParser().parseFromString(content, 'image/svg+xml');
      return [...svg.querySelectorAll('[data-gbdraw-feature-id], text')].map((node) => ({
        tag: node.localName, feature: node.getAttribute('data-gbdraw-feature-id'),
        label: node.getAttribute('data-label-feature-id'), text: node.textContent,
        part: node.getAttribute('data-gbdraw-feature-part'), d: node.getAttribute('d'),
        fill: node.getAttribute('fill'), display: node.getAttribute('display')
      }));
    }, content);
    const agree = async (stage, expectedText, feature) => {
      const artifact = await page.evaluate(async () => {
        const { state } = await import('./js/state.js');
        return { result: state.results.value[state.selectedResultIndex.value].content,
          mounted: state.svgContainer.value.querySelector('svg').outerHTML };
      });
      const result = await semantic(artifact.result);
      expect(await semantic(artifact.mounted)).toEqual(result);
      const download = page.waitForEvent('download');
      await page.getByRole('button', { name: 'SVG', exact: true }).click();
      const saved = testInfo.outputPath(`${stage}.svg`);
      await (await download).saveAs(saved);
      expect(await semantic(await fs.readFile(saved, 'utf8'))).toEqual(result);
      const parts = result.filter((node) => node.feature === feature.svg_id && node.part === 'block');
      expect(parts).toHaveLength(feature.location_parts.length);
      const labels = result.filter((node) => node.label === feature.svg_id);
      expect(labels.map((node) => node.text)).toEqual([expectedText]);
      return labels;
    };
    try {
      await load(seed);
      await generate();
      const feature = await target();
      expect(feature.location_parts).toHaveLength(gene === 'rps12' ? 3 : 1);
      const text = `edited ${name}`;
      await page.evaluate((id) => {
        const app = window.__GBDRAW_APP__;
        app.openFeatureEditorFromList(app.extractedFeatures.find((f) => f.id === id));
      }, feature.id);
      const label = page.locator('.feature-popup input[placeholder="Edit label text"]');
      await expect(label).toHaveValue(gene);
      if (['placement', 'both'].includes(action)) {
        const placement = page.getByLabel('Feature placement', { exact: true });
        await expect(placement.locator('option[value=outward]')).toHaveJSProperty('disabled', false);
        await placement.selectOption('outward');
      }
      if (action !== 'placement') {
        await label.fill(text);
        await page.getByRole('button', { name: 'Apply Label', exact: true }).click();
      }
      await page.getByRole('button', { name: 'Close feature popup', exact: true }).click();
      await expect.poll(() => page.evaluate(() => !window.__GBDRAW_HISTORY__.capturing.value)).toBe(true);
      const intent = await overrides();
      if (action !== 'placement') expect(intent.labels[feature.svg_id]).toBe(text);
      await agree('live', action === 'placement' ? gene : text, feature);
      if (action === 'save') {
        const pending = page.waitForEvent('download');
        await page.getByRole('button', { name: 'Save Session', exact: true }).click();
        const saved = testInfo.outputPath('edited.gbdraw-session.json.gz');
        await (await pending).saveAs(saved);
        await load(saved);
      } else if (action === 'mode') {
        for (const mode of ['Linear', 'Circular']) {
          await page.getByRole('button', { name: mode, exact: true }).click();
          await expect.poll(() => page.evaluate(() => window.__GBDRAW_APP__.mode)).toBe(mode.toLowerCase());
        }
      } else if (action === 'history') {
        await page.getByRole('button', { name: 'Undo', exact: true }).click();
        await expect.poll(overrides).toEqual({ labels: {}, placements: {} });
        await agree('undone', gene, feature);
        await page.getByRole('button', { name: 'Redo', exact: true }).click();
      }
      await expect.poll(overrides).toEqual(intent);
      await generate();
      expect(await overrides()).toEqual(intent);
      const regenerated = await target();
      expect(regenerated.id).toBe(feature.id);
      expect(regenerated.location_parts).toEqual(feature.location_parts);
      await agree('generated', action === 'placement' ? gene : text, regenerated);
      if (gene === 'rps12') {
        // The other trans-spliced copy keeps its own label and logical identity.
        const other = await page.evaluate(() => [...document.querySelectorAll('.origin-top text[data-label-feature-id]')]
          .filter((node) => node.textContent === 'rps12').map((node) => node.getAttribute('data-label-feature-id')));
        if (action !== 'placement') expect(other).not.toContain(feature.svg_id);
        expect(other).toHaveLength(action === 'placement' ? 2 : 1);
      }
      for (let step = 0; step < 3; step += 1) {
        await page.getByRole('button', { name: 'Zoom out', exact: true }).click();
      }
      await page.screenshot({ path: testInfo.outputPath('generated.png') });
      if (name.startsWith('M2 ')) {
        await page.setViewportSize({ width: 390, height: 844 });
        await expect(page.getByRole('button', { name: 'Generate Diagram', exact: true })).toBeEnabled();
        await page.screenshot({ path: testInfo.outputPath('mobile.png') });
      }
      expect(errors).toEqual([]);
      expect(external).toEqual([]);
    } finally {
      if (context) await context.close();
    }
  });
}
