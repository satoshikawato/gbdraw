const { test, expect } = require('@playwright/test');
const { readFileSync, writeFileSync } = require('node:fs');
const { resolve } = require('node:path');
const { gunzipSync } = require('node:zlib');
const { execFileSync } = require('node:child_process');
const { openApp, reveal, getDiagramWorkerActivity } = require('./helpers/app-lifecycle.cjs');

const fixture = resolve('gbdraw/web/gallery/sessions/tobacco-chloroplast.gbdraw-session.json');
const sessionInput = 'input[type="file"][accept^=".json,"]';
const valueControl = (page, slot, field) => page.getByRole('textbox', { name: `Circular track slot ${slot} ${field} value`, exact: true });
const unitControl = (page, slot, field) => page.getByRole('combobox', { name: `Circular track slot ${slot} ${field} unit`, exact: true });
const dialogs = new WeakMap();

const prepare = async page => {
  const messages = [];
  dialogs.set(page, messages);
  page.on('dialog', async dialog => {
    messages.push(dialog.message());
    await dialog.accept(dialog.type() === 'prompt' ? '' : undefined);
  });
  await page.addInitScript(() => {
    window.__C619_EVENTS__ = [];
    window.__GBDRAW_TEST_HOOKS__ = {
      onSessionLifecycleEvent: event => window.__C619_EVENTS__.push(event)
    };
  });
  await openApp(page);
};
const load = async (page, file = fixture, accepted = true) => {
  const before = dialogs.get(page).length;
  const readyBefore = await page.evaluate(() => window.__C619_EVENTS__.filter(event => event.name === 'interactiveReady').length);
  const input = page.locator(sessionInput);
  await input.setInputFiles(file);
  if (accepted) {
    await expect.poll(() => dialogs.get(page).length, { timeout: 180_000 }).toBeGreaterThan(before);
  } else {
    // Rejected Load reports an operation alert and clears the input at the
    // import owner's finally boundary before the preserved state is compared.
    await expect.poll(() => input.inputValue()).toBe('');
    await expect(page.getByRole('alert', { name: 'Operation error' })).toBeVisible();
    const error = await page.evaluate(() => window.__GBDRAW_APP__.errorLog);
    expect(error.stage).toBe('request-validation');
    expect(JSON.stringify(error)).not.toMatch(/1e|gc_content|invalid\.gbdraw-session/);
  }
  await page.waitForFunction(({ readyBefore, accepted }) => (
    !window.__GBDRAW_APP__.sessionImportPending
    && (!accepted || window.__C619_EVENTS__.filter(event => event.name === 'interactiveReady').length > readyBefore)
  ), { readyBefore, accepted });
  const messages = dialogs.get(page).slice(before);
  expect(messages.some(message => message === 'Session loaded successfully!')).toBe(accepted);
};
const openPanel = async page => {
  const button = page.locator('button[aria-controls="circular-custom-track-slots-panel"]');
  await reveal(button);
  if (await button.getAttribute('aria-expanded') !== 'true') await button.click();
  await expect(valueControl(page, 'features', 'Width')).toBeVisible();
  return button;
};
const counts = page => page.evaluate(() => [window.__GBDRAW_HISTORY__.getUndoCount(), window.__GBDRAW_HISTORY__.getRedoCount()]);
const expectCounts = async (page, undo, redo = 0) => expect.poll(() => counts(page)).toEqual([undo, redo]);
const scalar = (page, slot = 'gc_content', field = 'width') => page.evaluate(({ slot, field }) => (
  window.__GBDRAW_APP__.adv.circular_track_slots.find(row => row.id === slot)[field]
), { slot, field });
const snapshot = page => page.evaluate(async () => {
  const config = await import('./js/services/config.js');
  const state = (await import('./js/state.js')).state;
  const digest = async text => Array.from(new Uint8Array(await crypto.subtle.digest('SHA-256', new TextEncoder().encode(text))))
    .map(value => value.toString(16).padStart(2, '0')).join('');
  return {
    config: config.buildConfigData(),
    request: config.getCommittedCanonicalRenderRequest(),
    resultHashes: await Promise.all(state.results.value.map(result => digest(result.content))),
    history: [window.__GBDRAW_HISTORY__.getUndoCount(), window.__GBDRAW_HISTORY__.getRedoCount()]
  };
});
const workerCounts = async page => {
  const { constructions, initializations, helpers, runs } = await getDiagramWorkerActivity(page);
  return { constructions, initializations, helpers, runs };
};
const edit = async (page, text, slot = 'gc_content', field = 'Width') => {
  const input = valueControl(page, slot, field);
  await input.focus();
  await input.fill(text);
  await input.press('Tab');
};
const paste = async (page, text) => {
  await page.context().grantPermissions(['clipboard-read', 'clipboard-write']);
  await page.evaluate(text => navigator.clipboard.writeText(text), text);
  const input = valueControl(page, 'gc_content', 'Width');
  await input.focus();
  await input.press('ControlOrMeta+A');
  await input.press('ControlOrMeta+V');
  await input.press('Tab');
};
const selectUnit = async (page, unit, slot = 'gc_content', field = 'Width') => {
  const select = unitControl(page, slot, field);
  await select.focus();
  await select.selectOption(unit);
  await select.press('Tab');
};
const undo = async page => page.getByRole('button', { name: 'Undo', exact: true }).click();
const redo = async page => page.getByRole('button', { name: 'Redo', exact: true }).click();

const canonicalSlots = request => request.diagramOptions.tracks.circularTrackSlots;
const geometry = page => page.evaluate(async () => (await import('./js/state.js')).state.trackSlotResolvedGeometry.value);
const generate = async (page, valid = true) => {
  const before = await getDiagramWorkerActivity(page);
  await page.getByRole('button', { name: 'Generate Diagram', exact: true }).click();
  if (valid) {
    await expect.poll(async () => (await getDiagramWorkerActivity(page)).settledRuns,
      { timeout: 180_000 }).toBe(before.settledRuns + 1);
  }
  await page.waitForFunction(() => !window.__GBDRAW_APP__.processing);
  const error = await page.evaluate(() => window.__GBDRAW_APP__.errorLog);
  if (valid) expect(error).toBeNull();
  else expect(error).not.toBeNull();
};
const save = async (page, testInfo, name) => {
  const pending = page.waitForEvent('download');
  await page.getByRole('button', { name: 'Save Session', exact: true }).click();
  const download = await pending;
  const path = testInfo.outputPath(`${name}.gbdraw-session.json.gz`);
  await download.saveAs(path);
  await page.waitForFunction(() => !window.__GBDRAW_APP__.sessionSavePending);
  return { path, document: JSON.parse(gunzipSync(readFileSync(path))) };
};
const freshLoad = async (browser, testInfo, file, verify) => {
  const fresh = await browser.newPage({ baseURL: testInfo.project.use.baseURL });
  try {
    await prepare(fresh);
    await load(fresh, file);
    await verify(fresh);
  } finally {
    await fresh.close();
  }
};
// Compare every drawing element independently of XML serialization and
// transient preview/export attributes, including text and transforms.
const svgGeometry = (page, svg) => page.evaluate(content => {
  const doc = new DOMParser().parseFromString(content, 'image/svg+xml');
  return [...doc.querySelectorAll('svg,g,path,circle,ellipse,rect,line,polyline,polygon,text,textPath')].map(node => ({
    tag: node.localName,
    attributes: [...node.attributes].filter(attr => [
      'viewBox', 'transform', 'd', 'cx', 'cy', 'r', 'rx', 'ry', 'x', 'y',
      'x1', 'x2', 'y1', 'y2', 'points', 'width', 'height', 'font-size', 'startOffset'
    ].includes(attr.name)).map(attr => [attr.name, attr.value]).sort(),
    text: ['text', 'textPath'].includes(node.localName) ? node.textContent : null
  }));
}, svg);

test.beforeEach(async ({ page }) => {
  test.setTimeout(180_000);
  await prepare(page);
  await load(page);
});

for (const width of [1440, 390]) {
  test(`saved scalar controls are readable, accessible, and read-only at ${width}px`, async ({ page }, testInfo) => {
    await page.setViewportSize({ width, height: 1000 });
    const before = await snapshot(page);
    const workers = await workerCounts(page);
    const button = await openPanel(page);
    for (const [slot, field, value, unit] of [
      ['plastome_regions', 'Width', '20', 'px'], ['plastome_regions', 'Radius', '0.65', 'factor'],
      ['gc_content', 'Width', '0.08', 'factor'], ['gc_content', 'Radius', '0.56', 'factor']
    ]) {
      const input = valueControl(page, slot, field);
      const select = unitControl(page, slot, field);
      await expect(input).toHaveValue(value);
      await expect(select).toHaveValue(unit);
      await expect(input).toHaveAttribute('type', 'text');
      await expect(input).toHaveAttribute('inputmode', 'decimal');
      await expect(input).toHaveAttribute('aria-describedby', 'circular-measure-help');
      await input.focus();
      await input.press('Tab');
      await expect(select).toBeFocused();
      await select.press('Tab');
      const numericBox = await input.boundingBox();
      const unitBox = await select.boundingBox();
      expect(numericBox.width).toBeGreaterThan(25);
      expect(unitBox.x).toBeGreaterThanOrEqual(numericBox.x + numericBox.width);
      expect(unitBox.x + unitBox.width).toBeLessThanOrEqual(width);
    }
    await expect(page.locator('#circular-custom-track-slots-panel')).not.toContainText('[object Object]');
    await expect(page.locator('#circular-measure-help')).toContainText('65% (0.65 ×R)');
    await valueControl(page, 'features', 'Width').scrollIntoViewIfNeeded();
    await page.screenshot({ path: testInfo.outputPath(`controls-${width}.png`) });
    await button.click();
    await openPanel(page);
    expect(await snapshot(page)).toEqual(before);
    expect(await workerCounts(page)).toEqual(workers);
    await testInfo.attach('read-only-state', { body: JSON.stringify({ width, workers, rawSlots: before.config.adv.circular_track_slots }), contentType: 'application/json' });
  });
}

test('numeric and unit edits each have one History step, no derived status, and exact restores', async ({ page }) => {
  await openPanel(page);
  const before = await snapshot(page);
  const workers = await workerCounts(page);
  const initial = await counts(page);
  await selectUnit(page, 'px');
  await expectCounts(page, initial[0] + 1);
  expect(await scalar(page)).toEqual({ value: 0.08, unit: 'px' });
  await expect(page.locator('[data-generation-application-summary]')).toHaveCount(0);
  await edit(page, ' 1. ');
  await expectCounts(page, initial[0] + 2);
  expect(await scalar(page)).toEqual({ value: '1.', unit: 'px' });
  await selectUnit(page, 'factor');
  await expectCounts(page, initial[0] + 3);
  expect(await scalar(page)).toEqual({ value: '1.', unit: 'factor' });
  await selectUnit(page, 'factor');
  await expectCounts(page, initial[0] + 3);
  expect(await workerCounts(page)).toEqual(workers);
  await undo(page);
  await expectCounts(page, initial[0] + 2, 1);
  await expect(valueControl(page, 'gc_content', 'Width')).toHaveValue('1.');
  await expect(unitControl(page, 'gc_content', 'Width')).toHaveValue('px');
  await undo(page);
  await expectCounts(page, initial[0] + 1, 2);
  await expect(valueControl(page, 'gc_content', 'Width')).toHaveValue('0.08');
  await expect(unitControl(page, 'gc_content', 'Width')).toHaveValue('px');
  await redo(page);
  await expectCounts(page, initial[0] + 2, 1);
  await redo(page);
  await expectCounts(page, initial[0] + 3);
  await expect(valueControl(page, 'gc_content', 'Width')).toHaveValue('1.');
  await expect(unitControl(page, 'gc_content', 'Width')).toHaveValue('factor');
  const after = await snapshot(page);
  expect(after.request).toEqual(before.request);
  expect(after.resultHashes).toEqual(before.resultHashes);
  const restoredActivity = await workerCounts(page);
  expect(restoredActivity.runs).toBe(workers.runs);
});

test('plain input and complete suffix paste share the codec; invalid text stays visible and undoable', async ({ page }) => {
  await openPanel(page);
  const before = await snapshot(page);
  await selectUnit(page, 'px');
  await edit(page, '1.5');
  expect(await scalar(page)).toEqual({ value: '1.5', unit: 'px' });
  await paste(page, '65%');
  await expect(valueControl(page, 'gc_content', 'Width')).toHaveValue('0.65');
  await expect(unitControl(page, 'gc_content', 'Width')).toHaveValue('factor');
  await paste(page, '20px');
  await expect(valueControl(page, 'gc_content', 'Width')).toHaveValue('20');
  await expect(unitControl(page, 'gc_content', 'Width')).toHaveValue('px');
  const input = valueControl(page, 'gc_content', 'Width');
  for (const text of ['bad', '20p', '20pxx', '1epx', '1e', '0', '-1', 'Infinity', 'NaN', '1e309']) {
    const previous = await scalar(page);
    const count = (await counts(page))[0];
    await edit(page, text);
    await expectCounts(page, count + 1);
    await expect(input).toHaveValue(text);
    await expect(input).toHaveAttribute('aria-invalid', 'true');
    expect(await scalar(page)).toEqual({ value: text, unit: 'px' });
    const errorId = (await input.getAttribute('aria-describedby')).split(' ').find(id => id.endsWith('-error'));
    await expect(page.locator(`#${errorId}`)).toContainText('positive finite');
    await selectUnit(page, 'factor');
    await expect(input).toHaveValue(text);
    expect(await scalar(page)).toEqual({ value: text, unit: 'factor' });
    await undo(page);
    await expect(unitControl(page, 'gc_content', 'Width')).toHaveValue('px');
    await undo(page);
    expect(await scalar(page)).toEqual(previous);
    await redo(page);
    await expect(input).toHaveValue(text);
    await undo(page);
    expect(await scalar(page)).toEqual(previous);
  }
  const after = await snapshot(page);
  expect(after.request).toEqual(before.request);
  expect(after.resultHashes).toEqual(before.resultHashes);
});

test('IME preserves uncommitted text and adapts only after compositionend, with trim and keyboard unit edit', async ({ page }) => {
  await openPanel(page);
  const input = valueControl(page, 'gc_content', 'Width');
  const before = await scalar(page);
  const count = (await counts(page))[0];
  await input.focus();
  await input.dispatchEvent('compositionstart');
  await input.evaluate(element => {
    element.value = '65%';
    element.dispatchEvent(new InputEvent('input', { bubbles: true, isComposing: true, inputType: 'insertCompositionText' }));
  });
  await expect(input).toHaveValue('65%');
  expect(await scalar(page)).toEqual(before);
  await expectCounts(page, count);
  await input.dispatchEvent('compositionend');
  await expect(input).toHaveValue('0.65');
  expect(await scalar(page)).toEqual({ value: '0.65', unit: 'factor' });
  await input.press('Tab');
  await expectCounts(page, count + 1);
  const select = unitControl(page, 'gc_content', 'Width');
  await expect(select).toBeFocused();
  await select.press('ArrowUp');
  await select.press('Tab');
  await expect(select).toHaveValue('px');
  await expectCounts(page, count + 2);
  await edit(page, ' 1e-3 ');
  expect(await scalar(page)).toEqual({ value: '1e-3', unit: 'px' });
  await expect(input).toHaveValue('1e-3');
  await expectCounts(page, count + 3);
});

test('Auto preference is transient across input, Clear, History, panel remount, Load, and Reset', async ({ page }) => {
  await openPanel(page);
  const before = await snapshot(page);
  const workers = await workerCounts(page);
  await expect(unitControl(page, 'features', 'Width')).toHaveValue('factor');
  await selectUnit(page, 'px', 'features');
  expect(await snapshot(page)).toEqual(before);
  const autoId = (await valueControl(page, 'features', 'Width').getAttribute('aria-describedby')).split(' ').find(id => id.endsWith('-auto'));
  await expect(page.locator(`#${autoId}`)).toContainText('px (auto)');
  await edit(page, '1.5', 'features');
  expect(await scalar(page, 'features')).toEqual({ value: '1.5', unit: 'px' });
  await expectCounts(page, before.history[0] + 1);
  await edit(page, '', 'features');
  expect(await scalar(page, 'features')).toBeNull();
  await expectCounts(page, before.history[0] + 2);
  await undo(page);
  await expect(valueControl(page, 'features', 'Width')).toHaveValue('1.5');
  await expect(unitControl(page, 'features', 'Width')).toHaveValue('px');
  await undo(page);
  expect(await scalar(page, 'features')).toBeNull();
  await expect(valueControl(page, 'features', 'Width')).toHaveValue('');
  await redo(page);
  await expect(unitControl(page, 'features', 'Width')).toHaveValue('px');
  await redo(page);
  expect(await scalar(page, 'features')).toBeNull();
  await selectUnit(page, 'px', 'features');
  const remountWorkers = await workerCounts(page);
  const remountBefore = await snapshot(page);
  await (await openPanel(page)).click();
  await openPanel(page);
  await expect(unitControl(page, 'features', 'Width')).toHaveValue('factor');
  expect(await snapshot(page)).toEqual(remountBefore);
  expect(await workerCounts(page)).toEqual(remountWorkers);
  await selectUnit(page, 'px', 'features');
  await load(page);
  await openPanel(page);
  await expect(unitControl(page, 'features', 'Width')).toHaveValue('factor');
  await selectUnit(page, 'px', 'features');
  await page.getByRole('button', { name: 'Reset Settings', exact: true }).click();
  await openPanel(page);
  await expect(unitControl(page, 'features', 'Width')).toHaveValue('factor');
  const finalActivity = await workerCounts(page);
  expect(finalActivity.runs).toBe(workers.runs);
});

test('DOM edits Save and fresh Load retain draft apart from committed Result; failed admission is atomic', async ({ page, browser }, testInfo) => {
  await openPanel(page);
  const before = await snapshot(page);
  await selectUnit(page, 'px');
  await edit(page, '1.');
  const edited = await snapshot(page);
  const downloadPromise = page.waitForEvent('download');
  await page.getByRole('button', { name: 'Save Session', exact: true }).click();
  const download = await downloadPromise;
  const path = testInfo.outputPath('pending-draft.gbdraw-session.json.gz');
  await download.saveAs(path);
  await page.waitForFunction(() => !window.__GBDRAW_APP__.sessionSavePending);
  expect(await snapshot(page)).toEqual(edited);
  const saved = JSON.parse(gunzipSync(readFileSync(path)));
  const row = saved.config.adv.circular_track_slots.find(slot => slot.id === 'gc_content');
  expect(row.width).toEqual({ value: '1.', unit: 'px' });
  const fresh = await browser.newPage({ baseURL: testInfo.project.use.baseURL });
  try {
    await prepare(fresh);
    await load(fresh, path);
    await openPanel(fresh);
    await expect(valueControl(fresh, 'gc_content', 'Width')).toHaveValue('1.');
    await expect(unitControl(fresh, 'gc_content', 'Width')).toHaveValue('px');
    const restored = await snapshot(fresh);
    await testInfo.attach('fresh-result.svg', { body: await fresh.evaluate(() => window.__GBDRAW_APP__.results[0].content), contentType: 'image/svg+xml' });
    await testInfo.attach('original-result.svg', { body: saved.results[0].content, contentType: 'image/svg+xml' });
    expect(restored.config.adv.circular_track_slots).toEqual(edited.config.adv.circular_track_slots);
    expect(restored.request).toEqual(before.request);
    expect(restored.resultHashes).toEqual(before.resultHashes);
    await expectCounts(fresh, 0);
    const rejected = structuredClone(saved);
    rejected.config.adv.circular_track_slots.find(slot => slot.id === 'gc_content').width = { value: '1e', unit: 'px' };
    await load(fresh, { name: 'invalid.gbdraw-session.json', mimeType: 'application/json', buffer: Buffer.from(JSON.stringify(rejected)) }, false);
    expect(await snapshot(fresh)).toEqual(restored);
    await expect(valueControl(fresh, 'gc_content', 'Width')).toHaveValue('1.');
  } finally {
    await fresh.close();
  }
});

test('typed and legacy scalars display without writes; decimal, exponent, and precision edits project exactly', async ({ page }, testInfo) => {
  await openPanel(page);
  const cases = [
    [{ value: 20, unit: 'px' }, '20', 'px', 20],
    [{ value: 0.08, unit: 'factor' }, '0.08', 'factor', 0.08],
    [1.5, '1.5', 'factor', 1.5], ['1.5', '1.5', 'factor', 1.5],
    ['20px', '20', 'px', 20], ['65%', '0.65', 'factor', 0.65],
    [{ value: '1.', unit: 'px' }, '1.', 'px', 1],
    [{ value: '1e-3', unit: 'factor' }, '1e-3', 'factor', 0.001],
    [{ value: '0.12345678901234567', unit: 'factor' }, '0.12345678901234567', 'factor', Number('0.12345678901234567')],
    [{ value: '1e-12', unit: 'px' }, '1e-12', 'px', 1e-12],
    [{ value: '1.2345678901234567e+20', unit: 'px' }, '1.2345678901234567e+20', 'px', Number('1.2345678901234567e+20')],
    [null, '', 'factor', null], ['', '', 'factor', null]
  ];
  const observations = [];
  for (const [raw, text, unit, value] of cases) {
    await page.evaluate(async raw => {
      document.activeElement?.blur();
      const app = window.__GBDRAW_APP__;
      await window.__GBDRAW_HISTORY__.runUndoable('Set scalar fixture', () => {
        app.updateCircularTrackSlotMeasure(app.adv.circular_track_slots.find(row => row.id === 'gc_content'), 'width', raw);
      });
    }, raw);
    const before = await snapshot(page);
    const workers = await workerCounts(page);
    await expect(valueControl(page, 'gc_content', 'Width')).toHaveValue(text);
    await expect(unitControl(page, 'gc_content', 'Width')).toHaveValue(unit);
    await valueControl(page, 'gc_content', 'Width').focus();
    await valueControl(page, 'gc_content', 'Width').press('Tab');
    expect(await snapshot(page)).toEqual(before);
    expect(await workerCounts(page)).toEqual(workers);
    await edit(page, text);
    const projected = await page.evaluate(async () => {
      const { buildCircularTrackSlotPayload } = await import('./js/app/circular-track-slots.js');
      return buildCircularTrackSlotPayload(window.__GBDRAW_APP__.adv.circular_track_slots.find(row => row.id === 'gc_content')).width;
    });
    expect(projected).toEqual(value === null ? null : { value, unit });
    observations.push({ raw, text, unit, projected });
  }
  for (const unit of ['px', 'factor']) {
    await selectUnit(page, unit);
    await edit(page, '1.5');
    expect(await scalar(page)).toEqual({ value: '1.5', unit });
  }
  await testInfo.attach('scalar-matrix', { body: JSON.stringify(observations), contentType: 'application/json' });
});

test('Generate preserves four canonical pairs and same-value geometry; unit edit affects only its track; SVG download matches Result', async ({ page, browser }, testInfo) => {
  test.setTimeout(1_800_000);
  await openPanel(page);
  const original = await snapshot(page);
  const originalGeometry = await geometry(page);
  await generate(page);
  const first = await snapshot(page);
  expect(canonicalSlots(first.request)).toEqual(canonicalSlots(original.request));
  const firstGeometry = await geometry(page);
  const withoutOutputName = value => ({ ...value, records: value.records.map(({ resultName, ...record }) => record) });
  expect(withoutOutputName(firstGeometry)).toEqual(withoutOutputName(originalGeometry));
  const firstSvg = await page.evaluate(() => window.__GBDRAW_APP__.results[0].content);
  for (const [slot, field, value] of [
    ['plastome_regions', 'Width', '20'], ['plastome_regions', 'Radius', '0.65'],
    ['gc_content', 'Width', '0.08'], ['gc_content', 'Radius', '0.56']
  ]) await edit(page, value, slot, field);
  await generate(page);
  const same = await snapshot(page);
  expect(canonicalSlots(same.request)).toEqual(canonicalSlots(first.request));
  expect(await geometry(page)).toEqual(firstGeometry);
  expect(await svgGeometry(page, await page.evaluate(() => window.__GBDRAW_APP__.results[0].content)))
    .toEqual(await svgGeometry(page, firstSvg));
  await selectUnit(page, 'px');
  const pending = await snapshot(page);
  expect(pending.request).toEqual(same.request);
  expect(pending.resultHashes).toEqual(same.resultHashes);
  await generate(page);
  const changed = await snapshot(page);
  const expected = structuredClone(canonicalSlots(same.request));
  expected.find(row => row.id === 'gc_content').width = { value: 0.08, unit: 'px' };
  expect(canonicalSlots(changed.request)).toEqual(expected);
  const changedGeometry = await geometry(page);
  const a = firstGeometry.records[0];
  const b = changedGeometry.records[0];
  expect(b.axisRadiusPx).toBe(a.axisRadiusPx);
  for (const slot of a.slots.filter(row => row.slotId !== 'gc_content')) {
    expect(b.slots.find(row => row.slotId === slot.slotId)).toEqual(slot);
  }
  expect(b.slots.find(row => row.slotId === 'gc_content').widthPx).toBe(0.08);
  expect(b.slots.find(row => row.slotId === 'gc_content')).toEqual({
    ...a.slots.find(row => row.slotId === 'gc_content'), widthPx: 0.08, widthFactor: null
  });
  const currentSvg = await page.evaluate(() => window.__GBDRAW_APP__.results[0].content);
  const pendingDownload = page.waitForEvent('download');
  await page.getByRole('button', { name: 'SVG', exact: true }).click();
  const downloaded = await pendingDownload;
  const svgPath = testInfo.outputPath('current-result.svg');
  await downloaded.saveAs(svgPath);
  expect(await svgGeometry(page, readFileSync(svgPath, 'utf8'))).toEqual(await svgGeometry(page, currentSvg));
  const saved = await save(page, testInfo, 'generated');
  await freshLoad(browser, testInfo, saved.path, async fresh => {
    const restored = await snapshot(fresh);
    expect(restored.request).toEqual(changed.request);
    expect(restored.resultHashes).toEqual(changed.resultHashes);
    expect(await geometry(fresh)).toEqual(changedGeometry);
  });
  writeFileSync(testInfo.outputPath('generate-geometry.json'), JSON.stringify({
    originalPairs: canonicalSlots(original.request), first: firstGeometry, same: firstGeometry, changed: changedGeometry,
    svgPath, sessionPath: saved.path
  }, null, 2));
  writeFileSync(testInfo.outputPath('result-source.svg'), currentSvg);
});

test('invalid Generate keeps draft, committed request and Result; correction regenerates successfully', async ({ page }, testInfo) => {
  test.setTimeout(1_800_000);
  await openPanel(page);
  const before = await snapshot(page);
  const observations = [];
  for (const text of ['bad', '1e', '0', '-1', 'Infinity', 'NaN', '1e309', '20pxx']) {
    const historyBeforeEdit = await counts(page);
    await edit(page, text);
    await expectCounts(page, historyBeforeEdit[0] + 1);
    const draft = await scalar(page);
    const beforeGenerate = await snapshot(page);
    const workers = await workerCounts(page);
    await generate(page, false);
    await expect(valueControl(page, 'gc_content', 'Width')).toHaveValue(text);
    await expect(valueControl(page, 'gc_content', 'Width')).toHaveAttribute('aria-invalid', 'true');
    expect(await scalar(page)).toEqual(draft);
    const failed = await snapshot(page);
    expect(failed.config.adv.circular_track_slots).toEqual(beforeGenerate.config.adv.circular_track_slots);
    expect(failed.history).toEqual(beforeGenerate.history);
    expect(failed.request).toEqual(before.request);
    expect(failed.resultHashes).toEqual(before.resultHashes);
    expect((await workerCounts(page)).runs).toBe(workers.runs);
    observations.push({ text, draft, resultPreserved: true, history: failed.history,
      fullConfigEqual: JSON.stringify(failed.config) === JSON.stringify(beforeGenerate.config) });
  }
  await page.evaluate(() => {
    const app = window.__GBDRAW_APP__;
    app.updateCircularTrackSlotMeasure(app.adv.circular_track_slots.find(row => row.id === 'gc_content'), 'width', { value: 1, unit: 'em' });
  });
  await generate(page, false);
  expect(await scalar(page)).toEqual({ value: 1, unit: 'em' });
  const failed = await snapshot(page);
  expect(failed.request).toEqual(before.request);
  expect(failed.resultHashes).toEqual(before.resultHashes);
  await selectUnit(page, 'factor');
  await edit(page, '0.09');
  await generate(page);
  expect(canonicalSlots((await snapshot(page)).request).find(row => row.id === 'gc_content').width)
    .toEqual({ value: 0.09, unit: 'factor' });
  expect((await snapshot(page)).resultHashes).not.toEqual(before.resultHashes);
  await testInfo.attach('failed-generate-recovery', { body: JSON.stringify(observations), contentType: 'application/json' });
});

for (const continuation of ['disabled', 'inactive biological']) {
  test(`actual Save and fresh Load preserve ${continuation} Circular scalar drafts`, async ({ page, browser }, testInfo) => {
    await openPanel(page);
    await edit(page, '0.12345678901234567');
    if (continuation === 'disabled') await page.evaluate(() => {
      window.__GBDRAW_APP__.adv.circular_track_slots.find(row => row.id === 'gc_content').enabled = false;
    });
    if (continuation === 'inactive biological') await page.getByRole('button', { name: 'Linear', exact: true }).click();
    const before = await snapshot(page);
    const saved = await save(page, testInfo, continuation.replaceAll(' ', '-'));
    expect(saved.document.config.adv.circular_track_slots).toEqual(before.config.adv.circular_track_slots);
    await freshLoad(browser, testInfo, saved.path, async fresh => {
      const restored = await snapshot(fresh);
      expect(restored.config.adv.circular_track_slots).toEqual(before.config.adv.circular_track_slots);
      expect(restored.request).toEqual(before.request);
      expect(restored.resultHashes).toEqual(before.resultHashes);
      if (continuation === 'inactive biological') {
        await expect(fresh.getByRole('button', { name: 'Linear', exact: true }))
          .toHaveAttribute('aria-pressed', 'true');
        expect(restored.config.modeProfiles.activeMode).toBe('linear');
        await fresh.getByRole('button', { name: 'Circular', exact: true }).click();
      }
      await openPanel(fresh);
      await expect(valueControl(fresh, 'gc_content', 'Width')).toHaveValue('0.12345678901234567');
      if (continuation === 'inactive biological') {
        expect(before.config.modeProfiles.activeMode).toBe('linear');
      }
    });
  });
}

test('source-free settings Save and fresh Load retain an inactive disabled scalar and Auto', async ({ browser }, testInfo) => {
  const page = await browser.newPage({ baseURL: testInfo.project.use.baseURL });
  try {
    await prepare(page);
    await page.evaluate(slots => {
      const app = window.__GBDRAW_APP__;
      app.adv.circular_track_slots_enabled = true;
      app.adv.circular_track_slots = slots.filter(row => row.renderer !== 'annotations');
    }, JSON.parse(readFileSync(fixture, 'utf8')).config.adv.circular_track_slots);
    await openPanel(page);
    await edit(page, '1e-12');
    await selectUnit(page, 'px');
    await page.evaluate(() => { window.__GBDRAW_APP__.adv.circular_track_slots.find(row => row.id === 'gc_content').enabled = false; });
    await selectUnit(page, 'px', 'features');
    await page.getByRole('button', { name: 'Linear', exact: true }).click();
    const before = await snapshot(page);
    expect(before.request).toBeNull();
    expect(before.resultHashes).toEqual([]);
    const saved = await save(page, testInfo, 'settings-only');
    expect(saved.document.runMetadata).toBeUndefined();
    await freshLoad(browser, testInfo, saved.path, async fresh => {
      const restored = await snapshot(fresh);
      expect(restored.config).toEqual(before.config);
      expect(restored.request).toBeNull();
      expect(restored.resultHashes).toEqual([]);
      await fresh.getByRole('button', { name: 'Circular', exact: true }).click();
      await openPanel(fresh);
      await expect(valueControl(fresh, 'gc_content', 'Width')).toHaveValue('1e-12');
      await expect(unitControl(fresh, 'gc_content', 'Width')).toHaveValue('px');
      await expect(valueControl(fresh, 'features', 'Width')).toHaveValue('');
      await expect(unitControl(fresh, 'features', 'Width')).toHaveValue('factor');
    });
  } finally {
    await page.close();
  }
});

test('CLI-origin Session continues through numeric edits, Generate, Save and fresh Load', async ({ page, browser }, testInfo) => {
  test.setTimeout(1_800_000);
  const prefix = testInfo.outputPath('cli-origin');
  const args = ['-m', 'gbdraw.cli', 'circular', '--gbk', resolve('tests/fixtures/sessions/cli-web-mito.gb'),
    '--output', prefix, '--format', 'svg', '--save_session'];
  execFileSync('python', args, { env: { ...process.env, PYTHONPATH: process.cwd() } });
  const path = `${prefix}.gbdraw-session.json`;
  const cli = JSON.parse(readFileSync(path, 'utf8'));
  expect(cli.config).toBeUndefined();
  await load(page, path);
  const before = await snapshot(page);
  expect(before.request).toEqual(cli.renderRequest);
  expect(before.resultHashes).toHaveLength(1);
  const panel = page.locator('button[aria-controls="circular-custom-track-slots-panel"]');
  await reveal(panel);
  if (await panel.getAttribute('aria-expanded') !== 'true') await panel.click();
  await page.getByRole('checkbox', { name: 'Use custom stack', exact: true }).check();
  await openPanel(page);
  await selectUnit(page, 'px', 'features');
  await edit(page, '20', 'features');
  await generate(page);
  const generated = await snapshot(page);
  expect(canonicalSlots(generated.request).find(row => row.id === 'features').width).toEqual({ value: 20, unit: 'px' });
  const saved = await save(page, testInfo, 'cli-continuation');
  await freshLoad(browser, testInfo, saved.path, async fresh => {
    const restored = await snapshot(fresh);
    expect(restored.request).toEqual(generated.request);
    expect(restored.resultHashes).toEqual(generated.resultHashes);
    await openPanel(fresh);
    await expect(valueControl(fresh, 'features', 'Width')).toHaveValue('20');
    await expect(unitControl(fresh, 'features', 'Width')).toHaveValue('px');
  });
  writeFileSync(testInfo.outputPath('cli-origin-command.json'), JSON.stringify({ command: ['python', ...args], exit: 0, input: path }, null, 2));
});

test('nonfinite numeric action writes are rejected before History; finite, Auto, and text drafts retain their behavior', async ({ page }) => {
  await openPanel(page);
  const baseline = await snapshot(page);
  const workers = await workerCounts(page);
  for (const kind of ['nan', 'infinity', 'negativeInfinity', 'typedNan', 'typedInfinity', 'typedNegativeInfinity']) {
    const rejected = await page.evaluate(async kind => {
      const app = window.__GBDRAW_APP__, history = window.__GBDRAW_HISTORY__;
      const row = app.adv.circular_track_slots.find(slot => slot.id === 'gc_content');
      const before = row.width;
      const values = { nan: NaN, infinity: Infinity, negativeInfinity: -Infinity,
        typedNan: { value: NaN, unit: 'px' }, typedInfinity: { value: Infinity, unit: 'factor' },
        typedNegativeInfinity: { value: -Infinity, unit: 'px' } };
      let error = null;
      try {
        await history.runUndoable('Rejected numeric measure', () => app.updateCircularTrackSlotMeasure(row, 'width', values[kind]));
      } catch (caught) { error = caught.message; }
      return { error, same: row.width === before };
    }, kind);
    expect(rejected.error).toMatch(/positive finite px or factor scalar/);
    expect(rejected.same).toBe(true);
    expect(await snapshot(page)).toEqual(baseline);
    expect(await workerCounts(page)).toEqual(workers);
  }
  await page.evaluate(async () => {
    const app = window.__GBDRAW_APP__, row = app.adv.circular_track_slots.find(slot => slot.id === 'gc_content');
    await window.__GBDRAW_HISTORY__.runUndoable('Same measure', () => app.updateCircularTrackSlotMeasure(row, 'width', row.width));
  });
  expect(await snapshot(page)).toEqual(baseline);
  for (const [kind, expected] of [
    ['bareFinite', 0.12], ['typedFinite', { value: 3, unit: 'px' }],
    ['auto', null], ['typedLexeme', { value: '1e-3', unit: 'factor' }],
    ['textInfinity', { value: 'Infinity', unit: 'px' }],
    ['textIncomplete', { value: '1e', unit: 'px' }],
    ['textSuffix', { value: '20em', unit: 'px' }]
  ]) {
    await page.evaluate(async kind => {
      const app = window.__GBDRAW_APP__, row = app.adv.circular_track_slots.find(slot => slot.id === 'gc_content');
      const values = { bareFinite: 0.12, typedFinite: { value: 3, unit: 'px' }, auto: null,
        typedLexeme: { value: '1e-3', unit: 'factor' },
        textInfinity: { value: 'Infinity', unit: 'px' },
        textIncomplete: { value: '1e', unit: 'px' },
        textSuffix: { value: '20em', unit: 'px' } };
      await window.__GBDRAW_HISTORY__.runUndoable('Accepted measure draft', () => app.updateCircularTrackSlotMeasure(row, 'width', values[kind]));
    }, kind);
    expect(await scalar(page)).toEqual(expected);
    const edited = await snapshot(page);
    await page.evaluate(() => window.__GBDRAW_HISTORY__.undo());
    await page.evaluate(() => window.__GBDRAW_HISTORY__.redo());
    expect(await scalar(page)).toEqual(expected);
    const restored = await snapshot(page);
    expect(restored).toEqual(edited);
    expect(restored.request).toEqual(baseline.request);
    expect(restored.resultHashes).toEqual(baseline.resultHashes);
  }
  await expect(valueControl(page, 'gc_content', 'Width')).toHaveValue('20em');
  const beforeSave = await snapshot(page);
  const downloads = [];
  page.on('download', download => downloads.push(download));
  await page.getByRole('button', { name: 'Save Session', exact: true }).click();
  await expect.poll(() => page.evaluate(() => window.__GBDRAW_APP__.errorLog?.stage)).not.toBeNull();
  await page.waitForFunction(() => !window.__GBDRAW_APP__.sessionSavePending);
  expect(downloads).toEqual([]);
  expect(await snapshot(page)).toEqual(beforeSave);
  await generate(page, false);
  await expect(valueControl(page, 'gc_content', 'Width')).toHaveValue('20em');
  const failed = await snapshot(page);
  expect(failed.history).toEqual(beforeSave.history);
  expect(failed.request).toEqual(baseline.request);
  expect(failed.resultHashes).toEqual(baseline.resultHashes);
});
