const { test, expect } = require('@playwright/test');
const { join } = require('node:path');
const { readFileSync, writeFileSync } = require('node:fs');
const { gunzipSync } = require('node:zlib');
const { evaluateWithRetainedPromise, openApp, reveal, generateAndWaitForResult, getDiagramWorkerActivity } = require('./helpers/app-lifecycle.cjs');
const { inspectSafeDetails } = require('./helpers/operation-error.cjs');

const fixture = mode => join(process.cwd(), 'gbdraw/web/gallery/sessions', mode === 'linear'
  ? 'lambda_basic_linear.gbdraw-session.json' : 'HmmtDNA_basic_circular.gbdraw-session.json');
const load = async (page, file) => {
  await page.locator('input[accept^=".json,"]').setInputFiles(file);
  await page.waitForFunction(() => !window.__GBDRAW_APP__.sessionImportPending);
};
const setup = async (page, mode = 'circular') => {
  page.on('dialog', dialog => dialog.accept());
  await openApp(page);
  await load(page, fixture(mode));
  await page.waitForFunction(() => window.__GBDRAW_APP__.extractedFeatures.length);
  await page.evaluate(async () => {
    const a = window.__GBDRAW_APP__;
    Object.assign(a.newSpecRule, { feat: 'CDS', qual: 'product', val: '(?i)protein', color: '#f01234', cap: '' });
    await a.addSpecificRule();
  });
  await reveal(page.getByLabel('Color rule 1 pattern', { exact: true }));
};
const snapshot = page => page.evaluate(async () => {
  const { serializeResults } = await import('/gbdraw/web/js/services/config.js');
  const xml = content => {
    const parsed = new DOMParser().parseFromString(content, 'image/svg+xml');
    if (parsed.querySelector('parsererror')) throw new Error('Invalid Result SVG');
    return new XMLSerializer().serializeToString(parsed);
  };
  const a = window.__GBDRAW_APP__, h = window.__GBDRAW_HISTORY__;
  // Observe the same accepted Result boundary as Save/Export. XML round trips
  // remove only lexical empty-tag spelling; every node/attribute stays exact.
  return { canonical: JSON.stringify(a.manualSpecificRules), svg: xml(a.svgContent),
    result: serializeResults().map(r => ({ name: r.name, content: xml(r.content) })), undo: h.getUndoCount(), redo: h.getRedoCount() };
});
const draft = page => page.evaluate(() => {
  const a = window.__GBDRAW_APP__;
  const row = a.manualSpecificRules[0];
  return row ? { text: a.specificRulePattern(row), draft: a.specificRulePatternDraft(row) } : null;
});
const reject = async (page, text = '[') => {
  const field = page.getByLabel('Color rule 1 pattern', { exact: true });
  await reveal(field); await field.fill(text); await field.press('Tab');
  await expect.poll(async () => (await draft(page))?.draft?.error?.code).toBe('REGEX_SYNTAX');
  return field;
};
const status = page => page.locator('[data-color-rule-pattern-status]').first();

for (const mode of ['circular', 'linear']) for (const width of [1600, 390]) {
  test(`rejected ${mode} pattern stays correctable with keyboard and safe field diagnostics at ${width}px`, async ({ page }, info) => {
    test.setTimeout(180000);
    await page.setViewportSize({ width, height: width === 390 ? 740 : 1000 });
    const logs = [];
    page.on('console', message => logs.push(message.text()));
    page.on('pageerror', error => logs.push(error.message));
    await setup(page, mode);
    const before = await snapshot(page);
    const activity = await getDiagramWorkerActivity(page);
    const field = await reject(page, '😀[PRIVATE_PATTERN_SENTINEL');
    await expect(field).toHaveValue('😀[PRIVATE_PATTERN_SENTINEL');
    await expect(field).toHaveAttribute('aria-invalid', 'true');
    const described = await field.getAttribute('aria-describedby');
    expect(described).toBe(await status(page).getAttribute('id'));
    await expect(status(page)).toContainText('Not applied');
    await expect(status(page)).toContainText('last accepted rule');
    await expect(status(page)).toContainText('current Result');
    expect(await snapshot(page)).toEqual(before);
    const error = (await draft(page)).draft.error;
    expect(error).toMatchObject({ code: 'REGEX_SYNTAX', stage: 'rule-validation', context: { position: 1, positionUnit: 'python-character' } });
    await inspectSafeDetails(page, status(page).getByRole('alert', { name: 'Rule error' }), error);
    await page.screenshot({ path: info.outputPath('rejected-field.png') });
    const retry = status(page).getByRole('button', { name: 'Retry', exact: true });
    await retry.scrollIntoViewIfNeeded(); await retry.focus(); await retry.press('Enter');
    await expect.poll(async () => (await draft(page))?.draft?.pending).toBe(false);
    expect(await snapshot(page)).toEqual(before);
    const afterRetry = await getDiagramWorkerActivity(page);
    expect(afterRetry.constructions).toBe(activity.constructions);
    expect(afterRetry.helpers - activity.helpers).toBe(4); // captions + one Python syntax/match operation per attempt
    await page.evaluate(() => { const a = window.__GBDRAW_APP__; a.openRightDrawerTab('features'); a.closeRightDrawer(); a.openRightDrawerTab('features'); });
    await expect(field).toHaveValue('😀[PRIVATE_PATTERN_SENTINEL');
    await page.evaluate(() => window.__GBDRAW_APP__.closeRightDrawer());
    const revert = status(page).getByRole('button', { name: 'Revert', exact: true });
    await revert.scrollIntoViewIfNeeded(); await revert.focus(); await revert.press('Enter');
    await expect(field).toBeFocused();
    await expect(field).toHaveValue('(?i)protein');
    await expect(status(page)).toHaveCount(0);
    expect(await snapshot(page)).toEqual(before);
    await reject(page);
    await field.fill('(?P<enzyme>protein)');
    await field.dispatchEvent('change'); // retain input focus while the real asynchronous live edit commits
    await expect.poll(async () => (await draft(page))?.draft).toBe(null);
    await expect(field).toBeFocused();
    expect((await snapshot(page)).undo).toBe(before.undo + 1);
    expect(await page.evaluate(() => window.__GBDRAW_APP__.manualSpecificRules[0].val)).toBe('(?P<enzyme>protein)');
    expect(await page.evaluate(() => window.__GBDRAW_APP__.undoHistory())).toBe(true);
    expect(await page.evaluate(() => window.__GBDRAW_APP__.manualSpecificRules[0].val)).toBe('(?i)protein');
    await page.evaluate(() => window.__GBDRAW_APP__.redoHistory());
    expect(await page.evaluate(() => window.__GBDRAW_APP__.manualSpecificRules[0].val)).toBe('(?P<enzyme>protein)');
    expect(logs.join('\n')).not.toContain('PRIVATE_');
    console.log(JSON.stringify({ mode, width, helpers: afterRetry.helpers - activity.helpers, constructions: afterRetry.constructions, acceptedHistoryDelta: 1 }));
  });
}

test('native runtime initialization and transport preparation failures retain their field drafts; Retry succeeds atomically', async ({ page }) => {
  test.setTimeout(180000); page.on('dialog', d => d.accept());
  await openApp(page); await load(page, fixture('circular'));
  // Current Session preview stays lazy. Add a persisted accepted row without
  // invoking Python so the first pattern change really exercises cold init.
  await page.evaluate(() => {
    const a = window.__GBDRAW_APP__;
    a.manualSpecificRules.push({ feat: 'CDS', qual: 'product', val: 'NADH', color: '#f01234', cap: '' });
    const send = Worker.prototype.postMessage;
    Worker.prototype.postMessage = function (message, ...args) {
      if (message.type === 'init') {
        Worker.prototype.postMessage = send;
        message = { ...message, pyodideModuleUrl: '' };
      }
      return send.call(this, message, ...args);
    };
  });
  const before = await snapshot(page);
  await page.evaluate(() => window.__GBDRAW_APP__.setSpecificRuleField(0, 'val', '(?P<enzyme>NADH)'));
  expect((await draft(page)).draft.error).toMatchObject({ code: 'WORKER_INIT', stage: 'initialization' });
  expect(await snapshot(page)).toEqual(before);
  await page.evaluate(() => window.__GBDRAW_APP__.retrySpecificRulePattern(window.__GBDRAW_APP__.manualSpecificRules[0]));
  expect((await draft(page)).draft).toBe(null);
  expect((await snapshot(page)).undo).toBe(before.undo + 1);
  const accepted = await snapshot(page);
  await page.evaluate(() => {
    const send = Worker.prototype.postMessage;
    Worker.prototype.postMessage = function (message, ...args) {
      if (message.operation === 'evaluateRules') {
        Worker.prototype.postMessage = send;
        throw { code: 'RESOURCE_INVALID', operation: 'evaluateRules', stage: 'resource-staging', message: 'PRIVATE_PREPARATION_SENTINEL' };
      }
      return send.call(this, message, ...args);
    };
  });
  await page.evaluate(() => window.__GBDRAW_APP__.setSpecificRuleField(0, 'val', '(?i)NADH'));
  expect((await draft(page)).draft.error).toMatchObject({ code: 'RESOURCE_INVALID', stage: 'resource-staging' });
  expect((await draft(page)).draft.error.summary).not.toMatch(/invalid.*regular|PRIVATE_|Invalid rule/i);
  expect(await snapshot(page)).toEqual(accepted);
  await page.evaluate(() => window.__GBDRAW_APP__.retrySpecificRulePattern(window.__GBDRAW_APP__.manualSpecificRules[0]));
  expect((await draft(page)).draft).toBe(null);
  expect((await snapshot(page)).undo).toBe(accepted.undo + 1);
});

test('real History, failed Session rollback, fresh Save/Load, Export and Generate preserve accepted boundaries', async ({ page, browser }, info) => {
  test.setTimeout(180000); await setup(page);
  await generateAndWaitForResult(page);
  await page.evaluate(() => window.__GBDRAW_APP__.setSpecificRuleField(0, 'val', '(?i)NADH'));
  const field = await reject(page, '[PRIVATE_DRAFT_SENTINEL');
  const before = await snapshot(page);
  // An unrelated checkpoint replaces the rules while restoring its own setting.
  await page.evaluate(() => window.__GBDRAW_HISTORY__.runUndoableCheckpoint('Unrelated font', () => { window.__GBDRAW_APP__.adv.legend_font_size = 17; }));
  await page.evaluate(() => window.__GBDRAW_APP__.undoHistory());
  await expect(field).toHaveValue('[PRIVATE_DRAFT_SENTINEL');
  await page.evaluate(() => window.__GBDRAW_APP__.redoHistory());
  await expect(field).toHaveValue('[PRIVATE_DRAFT_SENTINEL');
  // Save and fresh Load must carry only the accepted rule and current Result.
  const downloaded = page.waitForEvent('download');
  await evaluateWithRetainedPromise(page, () => window.__GBDRAW_APP__.saveSessionWithTitle());
  const sessionPath = await (await downloaded).path();
  const bytes = readFileSync(sessionPath);
  const document = JSON.parse((bytes[0] === 0x1f ? gunzipSync(bytes) : bytes).toString());
  expect(JSON.stringify(document)).not.toContain('PRIVATE_DRAFT_SENTINEL');
  const context = await browser.newContext(); const fresh = await context.newPage();
  fresh.on('dialog', d => d.accept()); await openApp(fresh); await load(fresh, sessionPath);
  expect((await getDiagramWorkerActivity(fresh)).constructions).toBe(0);
  expect((await draft(fresh)).text).toBe('(?i)NADH');
  expect((await draft(fresh)).draft).toBe(null);
  expect((await snapshot(fresh)).canonical).toBe((await snapshot(page)).canonical);
  writeFileSync(info.outputPath('save-fresh-result-boundary.json'), JSON.stringify({ fresh: await snapshot(fresh), current: await snapshot(page), document }));
  expect((await snapshot(fresh)).svg).toBe((await snapshot(page)).svg);
  await context.close();
  const stable = await snapshot(page);
  await load(page, { name: 'invalid.json', mimeType: 'application/json', buffer: Buffer.from('{invalid') });
  await expect(field).toHaveValue('[PRIVATE_DRAFT_SENTINEL');
  expect(await snapshot(page)).toEqual(stable);
  // Fail after reset/commit began: actual Session owner must restore source first.
  await page.evaluate(async () => {
    const { state } = await import('/gbdraw/web/js/state.js');
    // The first palette write of the reset throws once, then the ref behaves normally again.
    const accessor = Object.getOwnPropertyDescriptor(Object.getPrototypeOf(state.selectedPalette), 'value');
    Object.defineProperty(state.selectedPalette, 'value', {
      configurable: true,
      get() { return accessor.get.call(this); },
      set() { delete state.selectedPalette.value; throw new Error('PRIVATE_COMMIT_SENTINEL'); }
    });
  });
  await load(page, fixture('circular'));
  await expect(page.getByLabel('Color rule 1 pattern', { exact: true })).toHaveValue('[PRIVATE_DRAFT_SENTINEL');
  expect(await snapshot(page)).toEqual(stable);
  const exportDownload = page.waitForEvent('download');
  await page.evaluate(() => window.__GBDRAW_APP__.downloadSVG());
  const exported = readFileSync(await (await exportDownload).path(), 'utf8');
  expect(exported).not.toContain('PRIVATE_');
  expect(exported).toContain('#f01234');
  const exportXml = await page.evaluate(text => new XMLSerializer().serializeToString(new DOMParser().parseFromString(text, 'image/svg+xml')), exported);
  expect(exportXml).toBe((await snapshot(page)).result[0].content);
  expect((await draft(page)).text).toBe('[PRIVATE_DRAFT_SENTINEL');
  await generateAndWaitForResult(page);
  expect((await draft(page)).text).toBe('(?i)NADH');
  expect((await draft(page)).draft).toBe(null);
  expect((await snapshot(page)).canonical).toBe(before.canonical);
  await reject(page);
  await page.evaluate(() => window.__GBDRAW_APP__.undoHistory()); // Undo document replacement drops its field draft only if its target changed; existing checkpoint can retain unchanged rules.
  // Explicit target edit Undo/Redo releases the replaced row's draft.
  await page.evaluate(() => window.__GBDRAW_APP__.setSpecificRuleField(0, 'val', '(?P<enzyme>protein)'));
  await reject(page);
  await page.evaluate(() => window.__GBDRAW_APP__.undoHistory());
  expect((await draft(page)).draft).toBe(null);
  await page.evaluate(() => window.__GBDRAW_APP__.redoHistory());
  expect((await draft(page)).draft).toBe(null);
  await reject(page);
  await load(page, fixture('linear'));
  expect(await draft(page)).toBe(null);
  await page.evaluate(async () => { const a = window.__GBDRAW_APP__; Object.assign(a.newSpecRule, { feat: 'CDS', qual: 'product', val: 'protein', color: '#f01234', cap: '' }); await a.addSpecificRule(); });
  await reject(page);
  await page.evaluate(() => window.__GBDRAW_APP__.resetSettings());
  expect(await draft(page)).toBe(null);
});

for (const boundary of ['drawer', 'mode cycle', 'remove', 'reorder', 'Session', 'new edit', 'new keystroke', 'Worker cancellation']) {
  test(`held real Worker edit stays isolated across ${boundary}`, async ({ page }) => {
    test.setTimeout(180000); await setup(page);
    if (boundary === 'mode cycle') {
      await generateAndWaitForResult(page);
      // Settle the existing mode-keyed mount's stroke binding before freezing
      // the accepted boundary. The held edit must preserve every SVG attribute.
      await page.evaluate(() => { const a = window.__GBDRAW_APP__; a.setDiagramMode('linear'); a.setDiagramMode('circular'); });
      await snapshot(page);
    }
    await page.evaluate(() => {
      const send = Worker.prototype.postMessage;
      Worker.prototype.postMessage = function (message, ...args) {
        if (message.operation === 'evaluateRules') {
          Worker.prototype.postMessage = send;
          window.__releasePattern = () => send.call(this, message, ...args);
          return;
        }
        return send.call(this, message, ...args);
      };
      window.__pendingPattern = window.__GBDRAW_APP__.setSpecificRuleField(0, 'val', '(?P<enzyme>NADH)');
    });
    await page.waitForFunction(() => window.__releasePattern);
    const before = await snapshot(page);
    expect((await draft(page)).draft.pending).toBe(true);
    if (boundary === 'drawer') await page.evaluate(() => { const a = window.__GBDRAW_APP__; a.openRightDrawerTab('features'); a.closeRightDrawer(); a.openRightDrawerTab('features'); });
    if (boundary === 'mode cycle') await page.evaluate(() => { const a = window.__GBDRAW_APP__; a.setDiagramMode('linear'); a.setDiagramMode('circular'); });
    if (boundary === 'remove') await page.evaluate(() => window.__GBDRAW_APP__.removeSpecificRule(0));
    if (boundary === 'reorder') await page.evaluate(async () => { const a = window.__GBDRAW_APP__; Object.assign(a.newSpecRule, { feat: 'CDS', qual: 'product', val: 'other', color: '#abcdef', cap: '' }); await a.addSpecificRule(); await a.moveSpecificRuleDown(0); });
    if (boundary === 'Session') await load(page, fixture('circular'));
    if (boundary === 'new edit') await page.evaluate(() => window.__GBDRAW_APP__.setSpecificRuleField(0, 'val', '(?P<enzyme>protein)'));
    if (boundary === 'new keystroke') await page.getByLabel('Color rule 1 pattern', { exact: true }).fill('[new');
    if (boundary === 'Worker cancellation') await page.evaluate(async () => { const { cancelDiagramGeneration } = await import('/gbdraw/web/js/services/diagram-generation.js'); cancelDiagramGeneration(); await window.__pendingPattern; });
    const current = await snapshot(page);
    if (boundary === 'mode cycle') writeFileSync(test.info().outputPath('mode-boundary.json'), JSON.stringify({ before, current }));
    await page.evaluate(async () => { window.__releasePattern(); await window.__pendingPattern; });
    expect(await snapshot(page)).toEqual(current);
    if (['drawer', 'mode cycle'].includes(boundary)) { expect((await draft(page)).text).toBe('(?P<enzyme>NADH)'); expect(await snapshot(page)).toEqual(before); }
    if (boundary === 'reorder') expect(await page.evaluate(() => { const a = window.__GBDRAW_APP__; return a.specificRulePattern(a.manualSpecificRules[1]); })).toBe('(?P<enzyme>NADH)');
    if (['remove', 'Session'].includes(boundary)) expect(await draft(page)).toBe(null);
    if (boundary === 'Worker cancellation') { expect((await draft(page)).draft.error).toBe(null); expect(await snapshot(page)).toEqual(before); }
    if (boundary === 'new keystroke') expect((await draft(page)).text).toBe('[new');
  });
}


test('cold helper cancellation retains its draft and Retry starts the real runtime once', async ({ page }) => {
  test.setTimeout(180000);
  page.on('dialog', d => d.accept()); await openApp(page); await load(page, fixture('circular'));
  await page.evaluate(() => {
    const a = window.__GBDRAW_APP__;
    a.manualSpecificRules.push({ feat: 'CDS', qual: 'product', val: 'NADH', color: '#f01234', cap: '' });
    const send = Worker.prototype.postMessage;
    Worker.prototype.postMessage = function (message, ...args) {
      if (message.type === 'init') { Worker.prototype.postMessage = send; window.__coldInitHeld = true; return; }
      return send.call(this, message, ...args);
    };
    window.__pendingPattern = a.setSpecificRuleField(0, 'val', '(?P<enzyme>NADH)');
  });
  await page.waitForFunction(() => window.__coldInitHeld);
  const before = await snapshot(page);
  await page.evaluate(async () => { const { cancelDiagramGeneration } = await import('/gbdraw/web/js/services/diagram-generation.js'); cancelDiagramGeneration(); await window.__pendingPattern; });
  expect((await draft(page)).text).toBe('(?P<enzyme>NADH)');
  expect((await draft(page)).draft.error).toBe(null);
  expect((await draft(page)).draft.pending).toBe(false);
  expect(await snapshot(page)).toEqual(before);
  await page.evaluate(() => window.__GBDRAW_APP__.retrySpecificRulePattern(window.__GBDRAW_APP__.manualSpecificRules[0]));
  expect((await draft(page)).draft).toBe(null);
  expect((await snapshot(page)).undo).toBe(before.undo + 1);
  expect((await getDiagramWorkerActivity(page)).constructions).toBe(2);
});
