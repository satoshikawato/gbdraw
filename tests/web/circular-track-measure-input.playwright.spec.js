const { test, expect } = require('@playwright/test');
const { readFileSync } = require('node:fs');
const { resolve } = require('node:path');
const { gunzipSync } = require('node:zlib');
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
  await page.locator(sessionInput).setInputFiles(file);
  await expect.poll(() => dialogs.get(page).length, { timeout: 180_000 }).toBeGreaterThan(before);
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

test('numeric and unit edits each have one History step, Pending feedback, and exact restores', async ({ page }) => {
  await openPanel(page);
  const before = await snapshot(page);
  const workers = await workerCounts(page);
  const initial = await counts(page);
  await selectUnit(page, 'px');
  await expectCounts(page, initial[0] + 1);
  expect(await scalar(page)).toEqual({ value: 0.08, unit: 'px' });
  await expect(page.locator('[data-generation-application-summary]')).toContainText('Pending');
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
