const { test, expect } = require('@playwright/test');
const { readFileSync } = require('node:fs');
const { gunzipSync } = require('node:zlib');
const { evaluateWithRetainedPromise, openApp } = require('./helpers/app-lifecycle.cjs');

const fixture = (name = 'HmmtDNA_basic_circular') => ({
  name: `${name}.gbdraw-session.json`, mimeType: 'application/json',
  buffer: readFileSync(`gbdraw/web/gallery/sessions/${name}.gbdraw-session.json`)
});
const setup = async (page) => {
  page.on('dialog', dialog => dialog.accept());
  const errors = [];
  page.on('pageerror', error => errors.push(error.message));
  page.on('console', message => { if (message.type() === 'error') errors.push(message.text()); });
  await openApp(page);
  await page.locator('input[type="file"][accept*="application/json"]').setInputFiles(fixture());
  await page.waitForFunction(() => !window.__GBDRAW_APP__.sessionImportPending);
  expect(errors).toEqual([]);
};
const capture = page => page.evaluate(async () => {
  const { state } = await import('/gbdraw/web/js/state.js');
  const config = await import('/gbdraw/web/js/services/config.js');
  const history = window.__GBDRAW_HISTORY__;
  window.beforeSession = {
    canonical: config.getCommittedCanonicalSession(),
    request: config.getCommittedCanonicalRenderRequest(),
    resources: config.getCommittedCanonicalSession()?.resources,
    results: state.results.value, file: state.files.c_gb,
    cache: state.losatCache.value, evidence: state.proteinIdentityManifest.value,
    undo: history.getUndoCount(), redo: history.getRedoCount(),
    config: JSON.stringify(config.buildConfigData())
  };
});
const intact = page => page.evaluate(async () => {
  const { state } = await import('/gbdraw/web/js/state.js');
  const config = await import('/gbdraw/web/js/services/config.js');
  const before = window.beforeSession;
  const history = window.__GBDRAW_HISTORY__;
  return {
    canonical: before.canonical === config.getCommittedCanonicalSession(),
    request: before.request === config.getCommittedCanonicalRenderRequest(),
    resources: before.resources === config.getCommittedCanonicalSession()?.resources,
    results: before.results === state.results.value,
    file: before.file === state.files.c_gb,
    cache: before.cache === state.losatCache.value,
    evidence: before.evidence === state.proteinIdentityManifest.value,
    history: before.undo === history.getUndoCount() && before.redo === history.getRedoCount(),
    config: before.config === JSON.stringify(config.buildConfigData())
  };
});
const expectIntact = async page => {
  const invariant = await intact(page);
  expect(invariant).toEqual(Object.fromEntries(Object.keys(invariant).map(key => [key, true])));
};

// Each case represents an independently reachable owner, rather than a disabled
// template assertion. Invalid arguments deliberately ensure busy wins admission.
const blockedActions = [
  ['selectResult', [0]], ['setOptionalNumberInputValue', [null, 'label_rotation', '45', true]],
  ['setMode', ['linear']], ['setCircularInputType', ['gff']],
  ['runAnalysis', []], ['cancelGeneration', []],
  ['setLabelFilterMode', ['Whitelist']], ['addWhitelistRule', []],
  ['removeWhitelistRule', [0]], ['removePriorityRule', [0]], ['resetSettings', []], ['resetLayout', []],
  ['clearLosatCache', []], ['setLosatPairFilename', ['missing', 'changed']],
  ['inspectCircularSourceRecords', []], ['refreshCircularRecordOrder', []],
  ['setCircularRecordPresentationSelector', ['#999']],
  ['setCircularRecordRow', [0, 4]], ['moveCircularRecordOrderUp', [0]],
  ['setLinearInputType', ['gff']], ['addLinearSeq', []],
  ['setLinearComparisonGlobalAction', ['none']], ['setLinearComparisonFile', ['missing', null]],
  ['setDepthTrackColor', [0, '#ffffff']], ['addCircularDepthTrack', []],
  ['addCircularTrackSlot', ['spacer']], ['addLinearTrackSlot', ['spacer']],
  ['addAnnotationSet', []], ['importAnnotationTableFile', []],
  ['addCustomColor', []], ['addFeature', []], ['setSpecificRuleField', [0, 'val', 'changed']],
  ['setFeatureVisibility', [null, 'hide']], ['resetAllLabelTextOverrides', []],
  ['addNewLegendEntry', []], ['sortLegendEntries', []], ['resetAllStrokes', []],
  ['resetCanvasPadding', []], ['editSessionTitle', []]
];

for (const operation of ['save', 'load']) {
  test(`${operation} occupies all semantic owners and keeps browsing available`, async ({ page }, info) => {
    await setup(page);
    await capture(page);
    await page.evaluate(async operation => {
      const service = await import('/gbdraw/web/js/services/config.js');
      let release;
      const gate = new Promise(resolve => { release = resolve; });
      window.releaseSession = release;
      const options = { beforeExport: () => gate, beforeImport: () => gate };
      window.sessionOperation = operation === 'save'
        ? service.exportSession('exclusive-save', options)
        : service.importSession({ target: { files: [new File(['{invalid'], 'invalid.json')], value: '' } }, options);
    }, operation);
    await page.waitForFunction(operation => window.__GBDRAW_APP__[
      operation === 'save' ? 'sessionSavePending' : 'sessionImportPending'
    ], operation);
    const outcomes = await evaluateWithRetainedPromise(page, async cases => {
      const app = window.__GBDRAW_APP__;
      const history = window.__GBDRAW_HISTORY__;
      const values = [];
      for (const [name, args] of cases) {
        if (typeof app[name] !== 'function') throw new Error(`Missing owner ${name}`);
        values.push([name, await app[name](...args)]);
      }
      values.push(['feature placement', await app.featurePlacementActions.setPlacement([], 'invalid')]);
      values.push(['undo', await history.undo()], ['redo', await history.redo()]);
      values.push(['history action', await history.runUndoable('forbidden', () => { throw new Error('ran'); })]);
      const service = await import('/gbdraw/web/js/services/config.js');
      values.push(['cross operation', app.sessionSavePending
        ? await service.importSession({ target: { files: [new File(['{}'], 'other.json')], value: '' } })
        : await service.exportSession('cross-save')]);
      return values;
    }, blockedActions);
    for (const [name, result] of outcomes) {
      expect(result?.status, name).toBe('busy');
      expect(result.reason, name).toMatch(/Retry after/);
    }
    await expect(page.getByRole('button', { name: 'Generate Diagram', exact: true })).toBeDisabled();
    await expect(page.getByRole('button', { name: 'Save Session', exact: true })).toBeDisabled();
    await expect(page.getByRole('button', { name: 'Load Session', exact: true })).toBeDisabled();
    await expect(page.locator('.upload-zone').first()).toHaveAttribute('tabindex', '-1');
    await expect(page.locator('.upload-zone').first()).toHaveAttribute('aria-disabled', 'true');
    const canvas = page.locator('.preview-canvas');
    const rect = await canvas.boundingBox();
    await page.mouse.move(rect.x + 8, rect.y + 8);
    await page.mouse.down();
    await page.mouse.move(rect.x + 28, rect.y + 28, { steps: 3 });
    await page.mouse.up();
    const beforeZoom = await page.evaluate(() => window.__GBDRAW_APP__.zoom);
    await page.mouse.wheel(0, -100);
    await expect.poll(() => page.evaluate(() => window.__GBDRAW_APP__.zoom)).toBeGreaterThan(beforeZoom);
    const search = page.getByRole('searchbox', { name: 'Search features', exact: true });
    await search.fill('gene');
    await search.press('ArrowLeft');
    await search.press('Enter');
    await expect(search).toHaveValue('gene');
    const browsedZoom = await page.evaluate(() => window.__GBDRAW_APP__.zoom);
    await expectIntact(page);
    const downloadPromise = operation === 'save' ? page.waitForEvent('download') : null;
    await page.evaluate(() => window.releaseSession());
    expect(await page.evaluate(async () => (await window.sessionOperation).status))
      .toBe(operation === 'save' ? 'saved' : 'error');
    if (downloadPromise) {
      const download = await downloadPromise;
      const path = info.outputPath('exclusive-session.json.gz');
      await download.saveAs(path);
      const saved = JSON.parse(gunzipSync(readFileSync(path)));
      expect(saved.ui.zoom).not.toBe(browsedZoom);
      expect(saved.config).toEqual(JSON.parse(await page.evaluate(() => window.beforeSession.config)));
    }
    await expectIntact(page);
    await expect(page.getByRole('button', { name: 'Load Session', exact: true })).toBeEnabled();
    await expect(page.locator('.upload-zone').first()).toHaveAttribute('tabindex', '0');
    await expect(page.locator('.upload-zone').first()).toHaveAttribute('aria-disabled', 'false');
  });
}

test('failed adoption restores source before transients and preserves all admitted owners', async ({ page }) => {
  await setup(page);
  await page.evaluate(async () => {
    await window.__GBDRAW_HISTORY__.runUndoable('retain edit', () => {
      window.__GBDRAW_APP__.adv.block_stroke_width = 2;
    });
    window.rollbackEvents = [];
    window.__GBDRAW_TEST_HOOKS__ = { onSessionLifecycleEvent: event => rollbackEvents.push(event.name) };
  });
  await capture(page);
  const replacement = fixture('lambda_basic_linear').buffer.toString('utf8');
  const outcome = await evaluateWithRetainedPromise(page, async text => {
    const service = await import('/gbdraw/web/js/services/config.js');
    return service.importSession({ target: { files: [new File([text], 'replacement.json')], value: '' } }, {
      beforePreviewMount: () => { throw new Error('controlled adoption failure'); }
    });
  }, replacement);
  expect(outcome.status).toBe('error');
  await expectIntact(page);
  const events = await page.evaluate(() => rollbackEvents);
  expect(events.indexOf('session-candidate-prepared')).toBeLessThan(events.indexOf('session-candidate-adopted'));
  expect(events.indexOf('session-rollback-source-restored')).toBeGreaterThan(-1);
  expect(events.indexOf('session-rollback-source-restored')).toBeLessThan(events.indexOf('session-rollback-transients-reconciled'));
  expect(await page.evaluate(() => window.__GBDRAW_APP__.sessionImportPending)).toBe(false);
});

test('Generate rejects Save and Load from processing publication through settlement', async ({ page }) => {
  test.setTimeout(180_000);
  await setup(page);
  await page.evaluate(() => {
    const history = window.__GBDRAW_HISTORY__;
    const original = history.runUndoableArtifactReplacement;
    let release;
    const gate = new Promise(resolve => { release = resolve; });
    window.releaseGenerate = release;
    history.runUndoableArtifactReplacement = async (...args) => {
      window.generateGateReached = true;
      await gate;
      return original(...args);
    };
    window.generationOperation = window.__GBDRAW_APP__.runAnalysis();
  });
  await page.waitForFunction(() => window.generateGateReached);
  const busy = await evaluateWithRetainedPromise(page, async () => {
    const service = await import('/gbdraw/web/js/services/config.js');
    return [await service.exportSession('blocked'), await service.importSession({
      target: { files: [new File(['{}'], 'blocked.json')], value: '' }
    })];
  });
  busy.forEach(result => {
    expect(result.status).toBe('busy');
    expect(result.reason).toMatch(/Generating diagram.*Retry/);
  });
  await expect(page.locator('[data-session-busy-reason]')).toContainText('Generating diagram');
  await page.evaluate(() => window.releaseGenerate());
  const generated = await evaluateWithRetainedPromise(page, async () => ({
    result: await window.generationOperation, error: window.__GBDRAW_APP__.errorLog
  }));
  expect(generated.result.status, JSON.stringify(generated.error)).toBe('ok');
  await expect(page.getByRole('button', { name: 'Save Session', exact: true })).toBeEnabled();
  await page.evaluate(async () => {
    const { state } = await import('/gbdraw/web/js/state.js');
    const send = Worker.prototype.postMessage;
    Worker.prototype.postMessage = function (message, ...args) {
      if (message.type === 'run' && state.labelReflowProcessing.value) {
        window.reflowRunHeld = true;
        window.releaseReflow = () => {
          Worker.prototype.postMessage = send;
          send.call(this, message, ...args);
        };
      } else send.call(this, message, ...args);
    };
    state.labelReflowForceRequestSeq.value += 1;
  });
  await page.waitForFunction(() => window.reflowRunHeld);
  const reflowBusy = await evaluateWithRetainedPromise(page, async () => {
    const service = await import('/gbdraw/web/js/services/config.js');
    return [await service.exportSession('blocked-reflow'), await service.importSession({
      target: { files: [new File(['{}'], 'blocked.json')], value: '' }
    })];
  });
  reflowBusy.forEach(result => {
    expect(result.status).toBe('busy');
    expect(result.reason).toMatch(/Updating diagram.*Retry/);
  });
  await expect(page.locator('[data-session-busy-reason]')).toContainText('Updating diagram');
  await page.evaluate(() => window.releaseReflow());
  await page.waitForFunction(() => !window.__GBDRAW_APP__.labelReflowProcessing);
  expect(await page.evaluate(() => window.__GBDRAW_APP__.labelReflowLastError)).toBeFalsy();
  await expect(page.getByRole('button', { name: 'Load Session', exact: true })).toBeEnabled();

});

test('teardown releases pending and canceled Save cannot download or clear a newer operation', async ({ page }) => {
  await setup(page);
  const outcome = await evaluateWithRetainedPromise(page, async () => {
    const service = await import('/gbdraw/web/js/services/config.js');
    const { state } = await import('/gbdraw/web/js/state.js');
    let release;
    const gate = new Promise(resolve => { release = resolve; });
    const saving = service.exportSession('disposed', { beforeExport: () => gate });
    await Promise.resolve();
    const published = state.sessionSavePending.value;
    service.disposeSessionOperations();
    const cleared = !state.sessionSavePending.value && !state.sessionImportPending.value;
    let releaseNew;
    const newGate = new Promise(resolve => { releaseNew = resolve; });
    const newer = service.importSession({ target: {
      files: [new File(['{invalid'], 'newer.json')], value: ''
    } }, { beforeImport: () => newGate });
    release();
    const canceled = await saving;
    const newerStillPending = state.sessionImportPending.value;
    releaseNew();
    const newStatus = (await newer).status;
    return { published, cleared, status: canceled.status, newerStillPending, newStatus };
  });
  expect(outcome).toEqual({ published: true, cleared: true, status: 'canceled', newerStillPending: true, newStatus: 'error' });
});

for (const failure of ['baseline-error', 'teardown-during-baseline']) {
  test(`${failure} rolls back an adopted candidate without losing History`, async ({ page }) => {
    await setup(page);
    await page.evaluate(async () => {
      await window.__GBDRAW_HISTORY__.runUndoable('retain edit', () => {
        window.__GBDRAW_APP__.adv.block_stroke_width = 2;
      });
    });
    await capture(page);
    const replacement = fixture('lambda_basic_linear').buffer.toString('utf8');
    const outcome = await evaluateWithRetainedPromise(page, async ({ text, failure }) => {
      const service = await import('/gbdraw/web/js/services/config.js');
      const history = window.__GBDRAW_HISTORY__;
      const initialize = history.initializeIntentBaseline;
      history.initializeIntentBaseline = async (...args) => {
        if (failure === 'baseline-error') throw new Error('controlled baseline failure');
        service.disposeSessionOperations();
        return initialize(...args);
      };
      try {
        return await window.__GBDRAW_APP__.importSession({ target: {
          files: [new File([text], 'replacement.json')], value: ''
        } });
      } finally {
        history.initializeIntentBaseline = initialize;
      }
    }, { text: replacement, failure });
    expect(outcome.status).toBe('error');
    await expectIntact(page);
    expect(await page.evaluate(() => window.__GBDRAW_APP__.sessionImportPending)).toBe(false);
  });
}

test('uncataloged multi-record draft Save prepares records privately without expanding live sources', async ({ page }, info) => {
  page.on('dialog', dialog => dialog.accept());
  await openApp(page);
  const text = readFileSync('tests/fixtures/sessions/cli-web-mito.gb', 'utf8');
  const download = page.waitForEvent('download');
  await page.evaluate(async text => {
    const app = window.__GBDRAW_APP__;
    app.setMode('linear');
    app.sessionTitle = 'Private draft';
    app.linearSeqs[0].gb = new File([text, text], 'two-records.gb');
    window.draftSource = app.linearSeqs[0];
    window.privateSaveEvents = [];
    window.__GBDRAW_TEST_HOOKS__ = { onSessionLifecycleEvent: event => {
      if (event.name === 'session-save-catalog-preparation-end') {
        privateSaveEvents.push({ count: app.linearSeqs.length, sameSource: app.linearSeqs[0] === draftSource });
      }
    } };
    window.draftSave = app.saveSessionWithTitle();
  }, text);
  expect((await evaluateWithRetainedPromise(page, async () => await window.draftSave)).status).toBe('saved');
  const file = info.outputPath('private-draft.json.gz');
  await (await download).saveAs(file);
  const saved = JSON.parse(gunzipSync(readFileSync(file)));
  expect(saved.renderRequest.records).toHaveLength(1);
  expect(saved.renderRequest.records[0]).toMatchObject({ cardinality: 'all', selector: null });
  const resourceId = saved.renderRequest.records[0].source.resourceId;
  expect(Buffer.from(saved.resources[resourceId].data, 'base64').toString('utf8')).toBe(text + text);
  expect(saved.results).toEqual([]);
  expect(await page.evaluate(() => privateSaveEvents)).toEqual([{ count: 1, sameSource: true }]);
});
