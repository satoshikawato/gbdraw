const { test, expect } = require('@playwright/test');
const { resolve } = require('node:path');
const { openApp, getDiagramWorkerActivity } = require('./helpers/app-lifecycle.cjs');
const { expectLiveEqualsGenerate } = require('./helpers/live-generate-parity.cjs');
const { appAction, history } = require('./helpers/live-generate-parity-steps.cjs');

const fixture = resolve('gbdraw/web/gallery/sessions/tobacco-chloroplast.gbdraw-session.json');

const visibleRuleColor = (page, color = '#ff0000') => page.evaluate(color => {
  const svg = (window.__GBDRAW_APP__.results?.[0]?.content || '').toLowerCase();
  const mounted = [...document.querySelectorAll('svg [fill], svg [style]')]
    .filter(node => `${node.getAttribute('fill') || ''} ${node.getAttribute('style') || ''}`.toLowerCase().includes(color));
  return { result: svg.split(color).length - 1, mounted: mounted.length };
}, color);

const edit = (page, label, action) => page.evaluate(async ({ label, action }) => {
  const { state } = await import('./js/state.js');
  await window.__GBDRAW_HISTORY__.runUndoable(label, () => {
    if (action === 'scale') state.activeDrawing().adv.scale_interval = 12345;
    if (action === 'scaleAgain') state.activeDrawing().adv.scale_interval = 23456;
    if (action === 'color') state.activeDrawing().manualSpecificRules[0].color = '#ff0000';
    if (action === 'predicate') state.activeDrawing().manualSpecificRules[0].val = 'psaA_NO_MATCH';
  });
}, { label, action });

const restore = (page, direction) => page.evaluate(direction => window.__GBDRAW_HISTORY__[direction](), direction);

test('History prepares cold rules once, skips warm config preparation, and rechecks changed predicates', async ({ page }) => {
  test.setTimeout(180_000);
  page.on('dialog', dialog => dialog.accept());
  await openApp(page);
  await page.locator('input[accept^=".json,"]').setInputFiles(fixture);
  await page.waitForFunction(() => !window.__GBDRAW_APP__.sessionImportPending);
  const context = page.context();
  const cdp = await context.newCDPSession(page);
  await cdp.send('Profiler.enable');
  await cdp.send('Profiler.startPreciseCoverage', { callCount: true, detailed: true });
  const prepareCalls = async () => (await cdp.send('Profiler.takePreciseCoverage')).result
    .filter(entry => entry.url.endsWith('/app/rule-matching.js'))
    .flatMap(entry => entry.functions)
    .filter(entry => entry.functionName === 'prepare')
    .reduce((total, entry) => total + Math.max(...entry.ranges.map(range => range.count)), 0);
  try {
    const before = await getDiagramWorkerActivity(page);
    expect(await visibleRuleColor(page)).toEqual({ result: 0, mounted: 0 });
    const baseGreen = await visibleRuleColor(page, '#00662c');
    expect(baseGreen.result).toBeGreaterThan(0);
    expect(baseGreen.mounted).toBeGreaterThan(0);
    await edit(page, 'Unrelated scale', 'scale');
    expect(await prepareCalls()).toBe(0);
    await restore(page, 'undo');
    expect(await prepareCalls()).toBe(1);
    expect(await visibleRuleColor(page)).toEqual({ result: 0, mounted: 0 });
    expect(await visibleRuleColor(page, '#00662c')).toEqual(baseGreen);
    await restore(page, 'redo');
    expect(await prepareCalls()).toBe(0);
    expect(await visibleRuleColor(page, '#00662c')).toEqual(baseGreen);
    const warm = await getDiagramWorkerActivity(page);
    expect(warm.constructions).toBe(before.constructions);
    expect(warm.helpers - before.helpers).toBe(1);
    expect(warm.runs).toBe(before.runs);
    await edit(page, 'Rule color', 'color');
    expect(await prepareCalls()).toBe(0);
    await restore(page, 'undo');
    expect(await prepareCalls()).toBe(1);
    expect(await visibleRuleColor(page)).toEqual({ result: 0, mounted: 0 });
    expect(await visibleRuleColor(page, '#00662c')).toEqual(baseGreen);
    await restore(page, 'redo');
    expect(await prepareCalls()).toBe(1);
    const colored = await visibleRuleColor(page);
    expect(colored.result).toBeGreaterThan(0);
    expect(colored.mounted).toBeGreaterThan(0);
    await edit(page, 'Rule predicate', 'predicate');
    expect(await prepareCalls()).toBe(0);
    await restore(page, 'undo');
    expect(await prepareCalls()).toBe(1);
    expect(await visibleRuleColor(page)).toEqual(colored);
    await restore(page, 'redo');
    expect(await prepareCalls()).toBe(1);
    expect(await visibleRuleColor(page)).toEqual({ result: 0, mounted: 0 });
    await edit(page, 'Independent scale', 'scaleAgain');
    expect(await prepareCalls()).toBe(0);
    expect(await page.evaluate(() => window.__GBDRAW_APP__.manualSpecificRules[0].val)).toBe('psaA_NO_MATCH');
    expect(await visibleRuleColor(page)).toEqual({ result: 0, mounted: 0 });
    await restore(page, 'undo');
    expect(await prepareCalls()).toBe(0);
    await restore(page, 'redo');
    expect(await prepareCalls()).toBe(0);
    expect(await visibleRuleColor(page)).toEqual({ result: 0, mounted: 0 });
    const after = await getDiagramWorkerActivity(page);
    expect(after.constructions).toBe(warm.constructions);
    expect(after.helpers - warm.helpers).toBe(1);
    expect(after.runs).toBe(before.runs);
  } finally {
    await cdp.send('Profiler.stopPreciseCoverage');
  }
});

// OV-347 (D-15-6 (3)): Undo and Redo of a rule color and a rule predicate
// change, made through the rule owner, show the Legend rows and fills Generate
// draws for the restored rules.
test('Redo of a rule predicate change shows what Generate draws', async ({ page }) => {
  test.setTimeout(240_000);
  page.on('dialog', dialog => dialog.accept());
  await openApp(page);
  await page.locator('input[accept^=".json,"]').setInputFiles(fixture);
  await page.waitForFunction(() => !window.__GBDRAW_APP__.sessionImportPending);
  await appAction(page, 'setSpecificRuleField', 0, 'color', '#ff0000');
  await appAction(page, 'setSpecificRuleField', 0, 'val', 'psaA_NO_MATCH');
  await history(page, 'undo');
  await history(page, 'undo');
  await history(page, 'redo');
  await history(page, 'redo');
  await expectLiveEqualsGenerate(page, { label: 'Redo of a rule predicate change' });
});
