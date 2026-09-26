const { test, expect } = require('@playwright/test');
const { createHash } = require('node:crypto');
const { gunzipSync } = require('node:zlib');
const { load, generate, snapshot, download } = require('./helpers/mode-transition.cjs');

const status = page => page.evaluate(async () => {
  const { getGenerationApplicationStatus } = await import('./js/services/config.js');
  const calls = [], originals = [];
  const spy = (owner, key) => {
    const original = owner[key];
    originals.push(() => { owner[key] = original; });
    owner[key] = () => { calls.push(key); throw new Error(`Unexpected Status side effect: ${key}`); };
  };
  for (const key of ['Worker', 'atob', 'btoa']) spy(window, key);
  for (const key of ['arrayBuffer', 'text']) spy(Blob.prototype, key);
  spy(crypto.subtle, 'digest'); spy(SVGElement.prototype, 'cloneNode');
  try { return { ...getGenerationApplicationStatus(), calls }; }
  finally { originals.reverse().forEach(restore => restore()); }
});

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
    expect(await status(page)).toMatchObject({ status:'clean', calls:[] });
    await page.locator('summary[aria-label="Labels"]').click();
    const labels = page.locator('#circular-label-mode');
    await labels.focus();
    await labels.selectOption('none');
    await labels.press('Tab');
    await expect.poll(() => page.evaluate(() => window.__GBDRAW_HISTORY__.undoLabel())).toBe('Change setting');
    expect(await status(page)).toMatchObject({ status:'pending', calls:[] });
    await generate(page);
    const b = await inspect(page, testInfo, 'b');
    expect(await status(page)).toMatchObject({ status:'clean', calls:[] });
    expect(b.request).not.toEqual(a.request);
    expect(b.selected).not.toBe(a.selected);

    await page.evaluate(() => window.__GBDRAW_HISTORY__.undo());
    const restored = await inspect(page, testInfo, 'undo');
    expect(await status(page)).toMatchObject({ status:'pending', calls:[] });
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
      expect(await status(fresh)).toMatchObject({ status:'pending', calls:[] });
      await generate(fresh);
      expectArtifact(await inspect(fresh, testInfo, 'fresh-generate'), b);
      expect(await status(fresh)).toMatchObject({ status:'clean', calls:[] });
      expect(fresh.externalRequests).toEqual([]);
    } finally { await fresh.context().close(); }

    await page.evaluate(() => window.__GBDRAW_HISTORY__.redo());
    const redone = await inspect(page, testInfo, 'redo');
    expect(await status(page)).toMatchObject({ status:'clean', calls:[] });
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


test('Live palette and its History retain a separate scale Pending', async ({ browser }) => {
  test.setTimeout(180000);
  const page = await load(browser);
  try {
    await generate(page);
    const before = await page.evaluate(() => window.__GBDRAW_APP__.results[0].content);
    await page.evaluate(async () => {
      const { state } = await import('./js/state.js');
      await window.__GBDRAW_HISTORY__.runUndoable('Change scale', () => { state.adv.scale_interval = 12345; });
      state.paletteInstantPreviewEnabled.value = true;
      await window.__GBDRAW_HISTORY__.runUndoable('Live palette', async () => {
        state.currentColors.value = { ...state.currentColors.value, CDS:'#123456' };
        await window.Vue.nextTick();
      });
    });
    const current = await status(page);
    expect(current).toMatchObject({status:'pending',calls:[]});
    expect(current.differences.some(path=>path.includes('scale.interval'))).toBe(true);
    const after = await page.evaluate(() => window.__GBDRAW_APP__.results[0].content);
    expect(after).not.toBe(before);
    for (const action of ['undo','redo']) {
      await page.evaluate(action => window.__GBDRAW_HISTORY__[action](),action);
      expect(await status(page)).toMatchObject({status:'pending',calls:[]});
    }
    await page.evaluate(async () => {
      const { state } = await import('./js/state.js'); state.adv.scale_interval = null;
      await window.Vue.nextTick();
    });
    expect(await status(page)).toMatchObject({status:'clean',calls:[]});
  } finally { await page.context().close(); }
});


test('Live global stroke advances only its applied fields and preserves scale Pending', async ({ browser }) => {
  test.setTimeout(180000);
  const page = await load(browser);
  try {
    await generate(page);
    await page.evaluate(async () => {
      const { state } = await import('./js/state.js');
      state.adv.scale_interval = 12345;
      await window.__GBDRAW_HISTORY__.runUndoable('Live stroke', async () => {
        state.adv.block_stroke_width = 2;
        await window.Vue.nextTick();
      });
    });
    const current = await status(page);
    expect(current).toMatchObject({status:'pending',calls:[]});
    expect(current.differences.every(path=>path.includes('scale.interval'))).toBe(true);
    expect(await page.locator('.gbdraw-preview-surface svg path[data-gbdraw-feature-id][data-gbdraw-feature-part="block"]').first().getAttribute('stroke-width')).toBe('2');
    await page.evaluate(async () => {
      const { state } = await import('./js/state.js');state.adv.scale_interval = null;
      await window.Vue.nextTick();
    });
    expect(await status(page)).toMatchObject({status:'clean',calls:[]});
    await page.evaluate(() => window.__GBDRAW_HISTORY__.undo());
    expect((await status(page)).status).not.toBe('clean');
    await page.evaluate(() => window.__GBDRAW_HISTORY__.redo());
    expect(await status(page)).toMatchObject({status:'pending',calls:[]});
  } finally { await page.context().close(); }
});
