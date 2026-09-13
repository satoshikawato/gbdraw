const { test, expect } = require('@playwright/test');
const { createHash } = require('node:crypto');
const { gunzipSync } = require('node:zlib');
const { load, generate, snapshot, download } = require('./helpers/mode-transition.cjs');

const hash = value => createHash('sha256').update(value).digest('hex');
const inspect = async (page, testInfo, name) => {
  const exported = await download(page, 'SVG', testInfo.outputPath(`${name}.svg`));
  const current = await snapshot(page);
  const editable = await page.evaluate(async () =>
    (await import('./js/services/config.js')).buildConfigData());
  const mounted = await page.evaluate(async () => {
    const { state } = await import('./js/state.js');
    const { serializeCleanSvg } = await import('./js/services/svg-serialization.js');
    return serializeCleanSvg(state.svgContainer.value.querySelector('svg'));
  });
  const result = {
    request: current.request, editable, markedMounted: current.markedMounted,
    selected: hash(current.result), mounted: hash(mounted), exported: hash(exported)
  };
  await testInfo.attach(name, { body: JSON.stringify(result), contentType: 'application/json' });
  return result;
};
const expectArtifact = (actual, expected) => {
  expect.soft(actual.selected, 'selected Result').toBe(expected.selected);
  expect.soft(actual.mounted, 'mounted SVG').toBe(expected.mounted);
  expect.soft(actual.exported, 'exported SVG').toBe(expected.exported);
  expect.soft(actual.markedMounted).toBe(true);
  expect.soft(actual.request, 'committed canonical request').toEqual(expected.request);
};

test('Undo Generate restores request A and Result A while preserving draft B through Save and fresh Load', async ({ browser }, testInfo) => {
  test.setTimeout(360000);
  const page = await load(browser);
  try {
    await generate(page);
    const a = await inspect(page, testInfo, 'a');
    await page.locator('summary[aria-label="Labels"]').click();
    const labels = page.locator('#circular-label-mode');
    await labels.focus();
    await labels.selectOption('none');
    await labels.press('Tab');
    await expect.poll(() => page.evaluate(() => window.__GBDRAW_HISTORY__.undoLabel())).toBe('Change setting');
    await generate(page);
    const b = await inspect(page, testInfo, 'b');
    expect(b.request).not.toEqual(a.request);
    expect(b.selected).not.toBe(a.selected);

    await page.evaluate(() => window.__GBDRAW_HISTORY__.undo());
    const restored = await inspect(page, testInfo, 'undo');
    expectArtifact(restored, a);
    expect(restored.editable).toEqual(b.editable);
    await expect(labels).toHaveValue('none');
    await page.screenshot({ path: testInfo.outputPath('undo.png') });

    const path = testInfo.outputPath('undo.gbdraw-session.json.gz');
    const saved = JSON.parse(gunzipSync(await download(page, 'Save Session', path)));
    expect.soft(saved.renderRequest).toEqual(a.request);
    expect(saved.config.form.labels_mode).toBe('none');
    const fresh = await load(browser, path);
    try {
      const loaded = await inspect(fresh, testInfo, 'fresh-load');
      expectArtifact(loaded, a);
      expect(loaded.editable).toEqual(b.editable);
      await generate(fresh);
      expectArtifact(await inspect(fresh, testInfo, 'fresh-generate'), b);
      expect(fresh.externalRequests).toEqual([]);
    } finally { await fresh.context().close(); }

    await page.evaluate(() => window.__GBDRAW_HISTORY__.redo());
    const redone = await inspect(page, testInfo, 'redo');
    expectArtifact(redone, b);
    expect(redone.editable).toEqual(b.editable);
    await page.evaluate(() => window.__GBDRAW_HISTORY__.undo());
    await page.evaluate(() => window.__GBDRAW_HISTORY__.undo());
    await expect(labels).toHaveValue('out');
    await generate(page);
    expectArtifact(await inspect(page, testInfo, 'undo-setting-generate'), a);
    expect(page.externalRequests).toEqual([]);
  } finally { await page.context().close(); }
});
