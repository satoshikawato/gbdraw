const { test, expect } = require('@playwright/test');
const fs = require('node:fs/promises');
const { gunzipSync } = require('node:zlib');
const { openApp } = require('./helpers/app-lifecycle.cjs');

const seed = 'gbdraw/web/gallery/sessions/HmmtDNA_basic_circular.gbdraw-session.json';
const external = new WeakMap();
const load = async (browser, file = seed) => {
  const context = await browser.newContext({ viewport: { width: 1600, height: 1000 } });
  const attempted = [];
  await context.route('**/*', route => {
    if (new URL(route.request().url()).hostname === '127.0.0.1') return route.continue();
    attempted.push(route.request().url());
    return route.abort();
  });
  const page = await context.newPage();
  external.set(page, attempted);
  page.on('dialog', dialog => dialog.dismiss());
  await openApp(page);
  await page.locator('input[accept^=".json,"]').setInputFiles(file);
  await expect.poll(() => page.evaluate(() => !window.__GBDRAW_APP__.sessionImportPending
    && window.__GBDRAW_APP__.extractedFeatures.length > 0), { timeout: 180000 }).toBe(true);
  return page;
};
const generate = async page => {
  const key = await page.evaluate(async () => (await import('./js/state.js')).state.resultGenerationKey.value);
  await page.getByRole('button', { name: 'Generate Diagram', exact: true }).click();
  await expect.poll(() => page.evaluate(async () => {
    const { state } = await import('./js/state.js');
    return { key: state.resultGenerationKey.value, processing: state.processing.value, error: state.errorLog.value };
  }), { timeout: 180000 }).toEqual({ key: key + 1, processing: false, error: null });
};
const placement = async (page, value) => {
  await page.locator('.drawer-toggle').click();
  const edit = page.locator('.right-drawer').getByRole('button', { name: 'Edit', exact: true }).first();
  await expect(edit).toBeVisible();
  await edit.click();
  const control = page.getByRole('combobox', { name: 'Feature placement', exact: true });
  await expect(control.locator('option[value=outward]')).toHaveJSProperty('disabled', false);
  if (value) await control.selectOption(value);
  await expect(control).toHaveValue('outward');
  await page.getByRole('button', { name: 'Close feature popup', exact: true }).click();
  await page.locator('.drawer-toggle').click();
};
const snapshot = page => page.evaluate(async () => {
  const { state } = await import('./js/state.js');
  const { getCommittedCanonicalRenderRequest } = await import('./js/services/config.js');
  const feature = state.extractedFeatures.value.find(f => f.biological_feature_id === 'fb8ff22d9');
  return { mode: state.mode.value, records: state.circularRecordList.value,
    placements: state.featurePlacementOverrides,
    feature: feature && { recordKey: feature.record_key, biologicalFeatureId: feature.biological_feature_id },
    featureCount: state.extractedFeatures.value.length,
    committedPlacements: getCommittedCanonicalRenderRequest()?.diagramOptions.featurePlacements || [],
    results: state.results.value.map(result => result.content) };
});
const switchMode = async (page, mode) => {
  await page.getByRole('button', { name: mode === 'circular' ? 'Circular' : 'Linear', exact: true }).click();
  await expect.poll(() => page.evaluate(() => window.__GBDRAW_APP__.mode)).toBe(mode);
  if (mode === 'circular') await expect.poll(() => page.evaluate(async () =>
    (await import('./js/state.js')).state.circularRecordDiscovery.status)).not.toBe('loading');
};
const structure = (page, svg) => page.evaluate(svg => {
  const visit = node => node.nodeType === Node.ELEMENT_NODE
    ? [node.tagName, Object.fromEntries([...node.attributes].map(a => [a.name, a.value]).sort()), [...node.childNodes].map(visit)]
    : node.textContent;
  return visit(new DOMParser().parseFromString(svg, 'image/svg+xml').documentElement);
}, svg);

for (const committed of [false, true]) {
  test(`${committed ? 'committed' : 'ungenerated'} Circular placement retains source identity through mode history and Save/Load`, async ({ browser }, testInfo) => {
    test.setTimeout(360000);
    const pages = [];
    try {
      const reference = await load(browser);
      pages.push(reference);
      await generate(reference);
      await placement(reference, 'outward');
      await generate(reference);
      const unswitched = await snapshot(reference);
      const expectedSvg = await structure(reference, unswitched.results[0]);
      await reference.context().close();

      const page = await load(browser);
      pages.push(page);
      await generate(page);
      await placement(page, 'outward');
      if (committed) await generate(page);
      const before = await snapshot(page);
      expect(before.feature).toEqual({ recordKey: 'record-1', biologicalFeatureId: 'fb8ff22d9' });
      expect(Object.values(before.placements)).toEqual([{ ...before.feature,
        placement: { kind: 'lane', side: 'outward', level: 1 } }]);
      await page.evaluate(async () => { window.__MODE_SOURCE__ = (await import('./js/state.js')).state.files.c_gb; });
      const checkReturned = async () => {
        await expect.poll(async () => (await snapshot(page)).feature).toEqual(before.feature);
        const returned = await snapshot(page);
        expect(returned.records).toEqual(before.records);
        expect(returned.placements).toEqual(before.placements);
        expect(returned.featureCount).toBe(before.featureCount);
        expect(returned.results).toEqual(before.results);
        expect(await page.evaluate(async () => (await import('./js/state.js')).state.files.c_gb === window.__MODE_SOURCE__)).toBe(true);
        await placement(page);
      };
      await switchMode(page, 'linear');
      await switchMode(page, 'circular');
      await checkReturned();
      await page.screenshot({ path: testInfo.outputPath('returned-circular.png') });

      // Each mode control remains one undoable operation; the placement is not rewritten.
      for (const mode of ['linear', 'circular']) {
        await page.getByRole('button', { name: /^Undo/ }).first().click();
        await expect.poll(() => page.evaluate(() => window.__GBDRAW_APP__.mode)).toBe(mode);
      }
      await checkReturned();
      for (const mode of ['linear', 'circular']) {
        await page.getByRole('button', { name: /^Redo/ }).first().click();
        await expect.poll(() => page.evaluate(() => window.__GBDRAW_APP__.mode)).toBe(mode);
      }
      await checkReturned();

      const pending = page.waitForEvent('download');
      await page.getByRole('button', { name: 'Save Session', exact: true }).click();
      const download = await pending;
      const saved = testInfo.outputPath(download.suggestedFilename());
      await download.saveAs(saved);
      const bytes = await fs.readFile(saved);
      const savedSession = JSON.parse((bytes[0] === 0x1f ? gunzipSync(bytes) : bytes).toString());
      expect(savedSession.renderRequest.diagramOptions.featurePlacements || []).toEqual(before.committedPlacements);
      await generate(page);
      const generated = await snapshot(page);
      expect(generated.committedPlacements).toEqual(unswitched.committedPlacements);
      expect(await structure(page, generated.results[0])).toEqual(expectedSvg);
      await fs.writeFile(testInfo.outputPath('roundtrip.svg'), generated.results[0]);

      const restored = await load(browser, saved);
      pages.push(restored);
      const loaded = await snapshot(restored);
      expect(loaded.placements).toEqual(before.placements);
      expect(loaded.feature).toEqual(before.feature);
      // Save flushes the mounted SVG, including its editor metadata and styles.
      expect(loaded.results).toHaveLength(savedSession.results.length);
      for (const [index, result] of savedSession.results.entries()) {
        expect(await structure(restored, loaded.results[index])).toEqual(await structure(restored, result.content));
      }
      await placement(restored);
      await generate(restored);
      expect(await structure(restored, (await snapshot(restored)).results[0])).toEqual(expectedSvg);
      for (let step = 0; step < 3; step += 1) await restored.getByRole('button', { name: 'Zoom out', exact: true }).click();
      await restored.screenshot({ path: testInfo.outputPath('restored-generated.png') });
      for (const checked of pages) expect(external.get(checked)).toEqual([]);
    } finally {
      for (const page of pages) await page.context().close();
    }
  });
}

for (const operation of ['replacement', 'removal']) {
  test(`Circular source ${operation} invalidates old identity and placement after mode inactivity`, async ({ browser }) => {
    test.setTimeout(180000);
    const page = await load(browser);
    try {
      await generate(page);
      await placement(page, 'outward');
      const before = await snapshot(page);
      const source = await page.evaluate(async () => {
        const { state } = await import('./js/state.js');
        window.__MODE_SOURCE__ = state.files.c_gb;
        window.__MODE_FEATURE__ = state.extractedFeatures.value[0];
        return (await import('./js/services/file-content-cache.js')).readFileText(state.files.c_gb);
      });
      await switchMode(page, 'linear');
      await switchMode(page, 'circular');
      if (operation === 'replacement') {
        // Keep the biological record ID: file/source identity must still invalidate its binding.
        await page.getByLabel('GenBank/DDBJ File', { exact: true }).setInputFiles({ name: 'replacement.gbk',
          mimeType: 'text/plain', buffer: Buffer.from(source.replace(/Homo sapiens/g, 'Replacement source')) });
        await expect.poll(() => page.evaluate(async () => (await import('./js/state.js')).state.circularRecordDiscovery.status)).toBe('ready');
      } else {
        await page.getByLabel('GenBank/DDBJ File', { exact: true })
          .locator('xpath=ancestor::div[@role="group"]')
          .getByRole('button', { name: /Remove/ }).click();
      }
      await expect.poll(() => page.evaluate(async () => Object.keys((await import('./js/state.js')).state.featurePlacementOverrides))).toEqual([]);
      const invalidation = await page.evaluate(async () => {
        const { state } = await import('./js/state.js');
        return { sameSource: state.files.c_gb === window.__MODE_SOURCE__,
          identities: state.circularRecordDiscovery.canonicalRecordIdentities,
          choices: window.__GBDRAW_APP__.featurePlacementActions.choices([window.__MODE_FEATURE__]) };
      });
      expect(invalidation.sameSource).toBe(false);
      expect(invalidation.identities).toEqual([]);
      expect(invalidation.choices.every(choice => !choice.enabled)).toBe(true);
      const changed = await snapshot(page);
      expect(changed.records.some(record => record.recordKey === before.feature.recordKey)).toBe(false);
      if (operation === 'replacement') await generate(page);
      else expect(changed.records).toEqual([]);
      expect(external.get(page)).toEqual([]);
    } finally { await page.context().close(); }
  });
}

test('no-placement Circular and Linear mode round trips remain generatable', async ({ browser }) => {
  test.setTimeout(240000);
  for (const [file, original, temporary] of [[seed, 'circular', 'linear'],
    ['gbdraw/web/gallery/sessions/lambda_basic_linear.gbdraw-session.json', 'linear', 'circular']]) {
    const page = await load(browser, file);
    try {
      await generate(page);
      const before = await snapshot(page);
      await switchMode(page, temporary);
      await switchMode(page, original);
      await generate(page);
      expect((await snapshot(page)).placements).toEqual({});
      expect(await structure(page, (await snapshot(page)).results[0])).toEqual(await structure(page, before.results[0]));
      expect(external.get(page)).toEqual([]);
    } finally { await page.context().close(); }
  }
});
