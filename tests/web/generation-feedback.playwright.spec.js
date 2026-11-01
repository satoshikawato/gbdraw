const { test, expect } = require('@playwright/test');
const { readFileSync, writeFileSync } = require('node:fs');
const { join } = require('node:path');
const { openApp } = require('./helpers/app-lifecycle.cjs');

const source = readFileSync(join(__dirname, '../test_inputs/HmmtDNA.gbk'), 'utf8');
// Applied authority is the committed canonical request; drafts never replace it.
const keepCommitted = (page, key) => page.evaluate(async key => {
  window[key] = (await import('./js/services/config.js')).getCommittedCanonicalRenderRequest();
  return Boolean(window[key]);
}, key);
const committedIs = (page, key) => page.evaluate(async key =>
  (await import('./js/services/config.js')).getCommittedCanonicalRenderRequest() === window[key], key);
const status = page => page.locator('[data-generation-progress]');
const generate = page => page.getByRole('button', { name: 'Generate Diagram', exact: true });
const prepare = async (page, baseURL) => {
  const externalRequests = [];
  await page.context().route(/^https?:\/\//, async route => {
    if (new URL(route.request().url()).origin !== new URL(baseURL).origin) {
      externalRequests.push(route.request().url());
      return route.abort('blockedbyclient');
    }
    return route.continue();
  });
  await page.addInitScript(() => {
    window.feedbackEvents = [];
    window.__GBDRAW_TEST_HOOKS__ = {
      onSessionLifecycleEvent(event) {
        window.feedbackEvents.push({ ...event, observedAt: performance.now() });
      }
    };
  });
  await openApp(page);
  await page.evaluate(() => {
    window.feedbackStatuses = [];
    window.Vue.watch(() => window.__GBDRAW_APP__.processingStatus, value => {
      window.feedbackStatuses.push({ value, at: performance.now() });
    }, { flush: 'sync' });
  });
  return externalRequests;
};
const resetObservations = page => page.evaluate(() => {
  window.feedbackStatuses = [];
  window.feedbackEvents = [];
});
const awaitSuccess = async page => {
  await expect.poll(() => page.evaluate(() => window.feedbackEvents.some(event => event.name === 'generate.completed')),
    { timeout: 180_000 }).toBe(true);
  await expect.poll(() => page.evaluate(() => {
    const app = window.__GBDRAW_APP__;
    return { processing: app.processing, error: app.errorLog?.summary || '', results: app.results.length };
  }), { timeout: 180_000 }).toEqual({ processing: false, error: '', results: 1 });
  await expect(status(page)).toHaveCount(0);
  await expect(generate(page)).toBeEnabled();
  await expect(page.locator('.gbdraw-preview-surface svg').first()).toBeVisible();
};
const record = async (page, testInfo, label) => {
  const evidence = await page.evaluate(() => ({
    statuses: window.feedbackStatuses,
    events: window.feedbackEvents,
    workers: window.__GBDRAW_DIAGRAM_WORKER_ACTIVITY__
  }));
  const path = testInfo.outputPath(`${label}.json`);
  writeFileSync(path, JSON.stringify(evidence, null, 2));
  await testInfo.attach(label, { path, contentType: 'application/json' });
  return evidence;
};
const screenshot = async (page, testInfo, label) => {
  const path = testInfo.outputPath(`${label}.png`);
  await page.screenshot({ path });
  await testInfo.attach(label, { path, contentType: 'image/png' });
};
const expectRenderStages = evidence => {
  const values = evidence.statuses.map(entry => entry.value);
  const expected = [
    'Preparing render inputs and session...',
    'Preparing diagram input resources...',
    'Rendering diagram...',
    'Finalizing diagram results...',
    'Preparing preview...',
    ''
  ];
  let previous = -1;
  for (const value of expected) {
    const index = values.indexOf(value, previous + 1);
    expect(index, JSON.stringify(values)).toBeGreaterThan(previous);
    previous = index;
  }
  const names = evidence.events.map(event => event.name);
  expect(names.indexOf('generate.completed')).toBeGreaterThan(names.indexOf('preview.post-bind-frame-completed'));
  expect(values.every(value => !value.includes('%'))).toBe(true);
};

test('Generate reports real cold/warm stages and preserves the successful Result through cancel, error, and retry', async ({ page, baseURL }, testInfo) => {
  test.setTimeout(180_000);
  const external = await prepare(page, baseURL);
  const setInput = text => page.evaluate(async content => {
    const app = window.__GBDRAW_APP__;
    app.files.c_gb = new File([content], 'HmmtDNA.gbk', { type: 'text/plain' });
    await window.Vue.nextTick();
  }, text);
  await setInput(source);
  await generate(page).click();
  // First use names the step the runtime is on (no percentages); the plain text shows only before the first step.
  const runtimeStatus = /^Preparing diagram runtime \(first use\)(\.\.\.|: step [1-5] of 5, [A-Za-z ]+\.\.\.)$/;
  await expect(status(page)).toHaveText(runtimeStatus);
  await expect(generate(page)).toBeDisabled();
  await screenshot(page, testInfo, 'cold-runtime-visible');
  await awaitSuccess(page);
  const cold = await record(page, testInfo, 'cold');
  expectRenderStages(cold);
  const coldRuntime = cold.statuses.map(entry => entry.value).filter(value => value.startsWith('Preparing diagram runtime'));
  expect(coldRuntime.length).toBeGreaterThan(0);
  expect(coldRuntime.every(value => runtimeStatus.test(value) && !value.includes('%'))).toBe(true);
  expect(coldRuntime.some(value => /: step \d of 5, /.test(value))).toBe(true);
  expect(cold.workers.constructions).toBe(1);
  expect(await keepCommitted(page, 'coldCommitted')).toBe(true);
  // Step 3a: the warm Generate reuses the catalog build and keeps this object.
  expect(cold.events.map(event => event.name)).toContain('catalog.admission-started');
  await page.evaluate(async () => {
    const { state } = await import('./js/state.js');
    window.coldCatalog = window.Vue.toRaw(state.featureCatalog.value);
  });

  await resetObservations(page);
  await generate(page).click();
  await awaitSuccess(page);
  const warm = await record(page, testInfo, 'warm');
  expectRenderStages(warm);
  expect(warm.statuses.map(entry => entry.value).filter(value => value.startsWith('Preparing diagram runtime'))).toEqual([]);
  expect(warm.workers.constructions).toBe(1);
  expect(warm.events.map(event => event.name)).toContain('catalog.reused');
  expect(warm.events.map(event => event.name)).not.toContain('catalog.admission-started');
  expect(await page.evaluate(async () => {
    const { state } = await import('./js/state.js');
    return Boolean(window.coldCatalog) && window.Vue.toRaw(state.featureCatalog.value) === window.coldCatalog;
  })).toBe(true);

  await page.evaluate(async () => {
    const app = window.__GBDRAW_APP__;
    window.successfulResults = app.results;
    // An edit since the shown Result, made as a History step like a user edit (UJ-10).
    await window.__GBDRAW_HISTORY__.runUndoable('Set scale interval', () => { app.adv.scale_interval = 12345; });
    window.feedbackHistory = [window.__GBDRAW_HISTORY__.getUndoCount(), window.__GBDRAW_HISTORY__.getRedoCount()];
    const config = await import('./js/services/config.js');
    window.appliedBeforeCancel = config.getCommittedCanonicalRenderRequest();
    window.__GBDRAW_TEST_HOOKS__.beforeDiagramGenerationResponse = () => new Promise(resolve => {
      window.releaseFeedbackResponse = resolve;
    });
  });
  await generate(page).click();
  await page.waitForFunction(() => Boolean(window.releaseFeedbackResponse));
  await expect(status(page)).toHaveText('Finalizing diagram results...');
  await expect(generate(page)).toBeDisabled();
  expect(await page.evaluate(() => window.__GBDRAW_APP__.results === window.successfulResults)).toBe(true);
  await screenshot(page, testInfo, 'finalization-visible');
  await page.getByRole('button', { name: /Cancel$/ }).click();
  await expect.poll(() => page.evaluate(() => window.__GBDRAW_APP__.processing)).toBe(false);
  await page.evaluate(() => {
    delete window.__GBDRAW_TEST_HOOKS__.beforeDiagramGenerationResponse;
    window.releaseFeedbackResponse();
  });
  expect(await page.evaluate(() => ({
    status: window.__GBDRAW_APP__.processingStatus,
    sameResult: window.__GBDRAW_APP__.results === window.successfulResults,
    preserved: window.__GBDRAW_APP__.failedGeneratePreservedResult
  }))).toEqual({ status: 'Canceled.', sameResult: true, preserved: true });
  await record(page, testInfo, 'canceled');
  await expect(page.locator('[data-generation-application-feedback]')).toHaveCount(0);
  expect(await page.evaluate(() => [window.__GBDRAW_HISTORY__.getUndoCount(), window.__GBDRAW_HISTORY__.getRedoCount()])).toEqual(await page.evaluate(() => window.feedbackHistory));
  expect(await committedIs(page, 'appliedBeforeCancel')).toBe(true);

  await setInput('invalid GenBank');
  await generate(page).click();
  await expect.poll(() => page.evaluate(() => window.__GBDRAW_APP__.processing), { timeout: 180_000 }).toBe(false);
  await expect(page.getByRole('alert', { name: 'Generation Error' })).toBeVisible();
  await expect(generate(page)).toBeEnabled();
  expect(await page.evaluate(() => window.__GBDRAW_APP__.results === window.successfulResults)).toBe(true);
  await record(page, testInfo, 'error');
  await expect(page.locator('[data-generation-application-feedback]')).toHaveCount(0);
  expect(await page.evaluate(() => [window.__GBDRAW_HISTORY__.getUndoCount(), window.__GBDRAW_HISTORY__.getRedoCount()])).toEqual(await page.evaluate(() => window.feedbackHistory));
  expect(await committedIs(page, 'appliedBeforeCancel')).toBe(true);

  await setInput(source);
  await page.setViewportSize({ width: 390, height: 844 });
  await resetObservations(page);
  await generate(page).click();
  await awaitSuccess(page);
  const retry = await record(page, testInfo, 'retry-mobile');
  expectRenderStages(retry);
  // Cancel ended the Worker: the retry's new Worker sends and admits the catalog.
  expect(retry.workers.constructions).toBe(2);
  expect(retry.events.map(event => event.name)).toContain('catalog.admission-started');
  expect(await committedIs(page, 'appliedBeforeCancel')).toBe(false);
  await screenshot(page, testInfo, 'ready-mobile');
  expect(external).toEqual([]);
});

// 3c B3 (PD-OI-101): after the first status a Generate yields a task, not
// frames. A Cancel that a user clicks while the main thread works at a later
// status is a queued input event, and Chromium dispatches it before the posted
// continuation. These cases click Cancel as a user does: the mouse is pressed
// on Cancel while the run waits for its first paint (frames held there), and
// released while the main thread is busy at the target status, so the mouseup
// and its click are queued there. Pins Chromium only.
const cancelByQueuedClick = async (page, testInfo, target) => {
  let mouseUp = null;
  const notifications = [];
  await page.route('**/__gbdraw-test/relay/**', async (route) => {
    // The release is sent while the page's task still runs and waits in the
    // renderer's input queue until that task ends.
    const busy = route.request().url().endsWith('/busy');
    notifications.push({ busy, at: Date.now() });
    if (busy && !mouseUp) mouseUp = page.mouse.up();
    await route.fulfill({ status: 204, body: '' });
  });
  await page.evaluate((targetStatus) => {
    const nativeFrame = window.requestAnimationFrame;
    const queued = { target: targetStatus, armed: false, holding: false, held: [], busy: null, input: [] };
    window.queuedCancel = queued;
    // A request made by the busy task itself reaches Playwright only when that
    // task ends (its Network event is sent then), so a Worker that is not busy
    // makes it at once.
    const relay = new Worker(URL.createObjectURL(new Blob([`onmessage = async ({ data }) => {
      await fetch(${JSON.stringify(`${location.origin}/__gbdraw-test/relay/`)} + data, { cache: 'no-store' }).catch(() => {});
      postMessage(data);
    };`], { type: 'text/javascript' })));
    for (const type of ['mousedown', 'mouseup', 'click']) {
      document.addEventListener(type, (event) => {
        if (!event.target.closest?.('button')?.textContent?.includes('Cancel')) return;
        queued.input.push({ type, timeStamp: event.timeStamp, dispatchedAt: performance.now() });
      }, { capture: true });
    }
    // afterPaint keeps the function it read, so the hold is a flag, not a swap.
    const heldFrame = (callback) => {
      if (!queued.holding) return nativeFrame(callback);
      queued.held.push(callback);
      return 0;
    };
    queued.release = () => {
      queued.holding = false;
      window.requestAnimationFrame = nativeFrame;
      for (const callback of queued.held.splice(0)) nativeFrame(callback);
    };
    const stopWatch = window.Vue.watch(() => window.__GBDRAW_APP__.processingStatus, (value) => {
      if (value === 'Preparing input files...' && !queued.armed) {
        queued.armed = true;
        queued.holding = true;
        window.requestAnimationFrame = heldFrame;
      }
      if (value !== queued.target || queued.busy) return;
      const start = performance.now();
      relay.postMessage('busy');
      while (performance.now() - start < 400) { /* the user's click arrives while this task runs */ }
      const end = performance.now();
      queued.busy = { start, end };
    }, { flush: 'sync' });
    queued.stop = () => { stopWatch(); relay.terminate(); };
    // The relay runs and has made one request before the Generate starts.
    return new Promise((resolve) => {
      relay.onmessage = () => { relay.onmessage = null; resolve(); };
      relay.postMessage('ready');
    });
  }, target);
  await generate(page).click();
  const cancel = page.getByRole('button', { name: /Cancel$/ });
  await expect(cancel).toBeVisible();
  const box = await cancel.boundingBox();
  await page.mouse.move(box.x + box.width / 2, box.y + box.height / 2);
  await page.mouse.down();
  await expect.poll(() => page.evaluate(() => window.queuedCancel.input.map(({ type }) => type))).toEqual(['mousedown']);
  await page.evaluate(() => window.queuedCancel.release());
  await expect.poll(() => page.evaluate(() => window.__GBDRAW_APP__.processing), { timeout: 60_000 }).toBe(false);
  await mouseUp;
  await page.unroute('**/__gbdraw-test/relay/**');
  const queued = await page.evaluate(() => {
    const { busy, input } = window.queuedCancel;
    window.queuedCancel.stop();
    return { busy, input };
  });
  expect(queued.busy).not.toBeNull();
  const path = testInfo.outputPath('queued-click.json');
  writeFileSync(path, JSON.stringify({ ...queued, notifications }, null, 2));
  await testInfo.attach('queued-click', { path, contentType: 'application/json' });
  const [, up, click] = queued.input;
  expect(queued.input.map(({ type }) => type)).toEqual(['mousedown', 'mouseup', 'click']);
  // The release reached the page during the busy task (an input event's
  // timeStamp is when it arrived; compared on the page's clock, not Node's) and
  // was dispatched after it.
  expect(up.timeStamp).toBeLessThan(queued.busy.end);
  expect(click.dispatchedAt).toBeGreaterThanOrEqual(queued.busy.end);
};
const warmReady = async (page, baseURL) => {
  const external = await prepare(page, baseURL);
  await page.evaluate(async content => {
    const app = window.__GBDRAW_APP__;
    app.files.c_gb = new File([content], 'HmmtDNA.gbk', { type: 'text/plain' });
    await window.Vue.nextTick();
  }, source);
  await generate(page).click();
  await awaitSuccess(page);
  await page.evaluate(() => { window.committedResults = window.__GBDRAW_APP__.results; });
  await resetObservations(page);
  return external;
};
const expectCanceledWarm = async (page, testInfo, label) => {
  const evidence = await record(page, testInfo, label);
  expect(await page.evaluate(() => ({
    status: window.__GBDRAW_APP__.processingStatus,
    sameResult: window.__GBDRAW_APP__.results === window.committedResults
  }))).toEqual({ status: 'Canceled.', sameResult: true });
  expect(evidence.statuses.filter(entry => entry.value === 'Canceled.')).toHaveLength(1);
  // GE-09: no Worker request was active, so the warm Worker stays.
  expect(evidence.workers.constructions).toBe(1);
  expect(evidence.workers.instances[0].terminated).toBe(false);
  return evidence;
};

test('A Cancel clicked during input resolution keeps the warm Worker for the next Generate', async ({ page, baseURL }, testInfo) => {
  test.setTimeout(240_000);
  const external = await warmReady(page, baseURL);
  const runsBefore = await page.evaluate(() => window.__GBDRAW_DIAGRAM_WORKER_ACTIVITY__.instances[0].runs.length);
  await cancelByQueuedClick(page, testInfo, 'Preparing render inputs and session...');
  const evidence = await expectCanceledWarm(page, testInfo, 'cancel-input-resolution');
  // Seen before the post: the run sent no request.
  expect(evidence.workers.instances[0].runs).toHaveLength(runsBefore);
  await resetObservations(page);
  await generate(page).click();
  await awaitSuccess(page);
  const next = await record(page, testInfo, 'after-cancel-warm');
  expect(next.statuses.some(entry => entry.value.startsWith('Preparing diagram runtime'))).toBe(false);
  expect(next.workers.constructions).toBe(1);
  expect(external).toEqual([]);
});

test('A Cancel clicked during SVG admission ends the Generate before the candidate mounts', async ({ page, baseURL }, testInfo) => {
  test.setTimeout(240_000);
  const external = await warmReady(page, baseURL);
  await cancelByQueuedClick(page, testInfo, 'Preparing preview...');
  const evidence = await expectCanceledWarm(page, testInfo, 'cancel-admission');
  expect(evidence.events.filter(event => event.name === 'preview.mount-observed')).toEqual([]);
  expect(external).toEqual([]);
});

test('Linear comparison preparation, real search counts, cache reuse, and no-comparison generation use the same status', async ({ page, baseURL }, testInfo) => {
  test.setTimeout(180_000);
  const external = await prepare(page, baseURL);
  await page.evaluate(async content => {
    const app = window.__GBDRAW_APP__;
    app.selectDrawing('linear');
    await window.Vue.nextTick();
    while (app.linearSeqs.length < 2) app.addLinearSeq();
    for (let i = 0; i < 2; i += 1) {
      app.setLinearSeqPrimaryFile(i, 'gb', new File([content], `record-${i}.gbk`, { type: 'text/plain' }));
    }
    await app.setLinearComparisonGlobalAction('losat');
    await window.Vue.nextTick();
  }, source);
  await page.getByRole('combobox', { name: 'LOSAT execution', exact: true }).selectOption('serial');
  await generate(page).click();
  await awaitSuccess(page);
  const search = await record(page, testInfo, 'linear-search');
  expectRenderStages(search);
  expect(search.statuses.map(entry => entry.value)).toEqual(expect.arrayContaining([
    'Preparing comparisons...', 'Preparing LOSAT jobs...',
    'Running LOSAT: 0/1 source job complete', 'Running LOSAT: 1/1 source job complete',
    'Preparing nucleotide comparison results...'
  ]));

  await resetObservations(page);
  await generate(page).click();
  await awaitSuccess(page);
  const cached = await record(page, testInfo, 'linear-cached');
  expectRenderStages(cached);
  expect(cached.statuses.map(entry => entry.value)).toContain('Using cached LOSAT results...');
  expect(cached.statuses.some(entry => entry.value.startsWith('Running LOSAT:'))).toBe(false);

  await page.evaluate(() => window.__GBDRAW_APP__.setLinearComparisonGlobalAction('none'));
  await resetObservations(page);
  await generate(page).click();
  await awaitSuccess(page);
  const none = await record(page, testInfo, 'linear-none');
  expectRenderStages(none);
  expect(none.statuses.some(entry => /comparison|LOSAT/.test(entry.value))).toBe(false);
  expect(external).toEqual([]);
});


test('S03 Linear LOSAT draft edits never replace the committed request before Generate', async ({ page, baseURL }) => {
  test.setTimeout(180000);
  await prepare(page,baseURL);
  await page.evaluate(async content => {
    const app=window.__GBDRAW_APP__;
    app.selectDrawing('linear');await window.Vue.nextTick();
    while(app.linearSeqs.length<2) app.addLinearSeq();
    for(let i=0;i<2;i++) app.setLinearSeqPrimaryFile(i,'gb',new File([content],`record-${i}.gbk`,{type:'text/plain'}));
    await app.setLinearComparisonGlobalAction('losat');await window.Vue.nextTick();
  },source);
  await page.getByRole('combobox', { name: 'LOSAT execution', exact: true }).selectOption('serial');
  await generate(page).click();
  await expect.poll(()=>page.evaluate(()=>window.feedbackEvents.some(event=>event.name==='generate.completed')),
    {timeout:180000}).toBe(true);
  expect(await page.evaluate(()=>window.__GBDRAW_APP__.errorLog?.summary||'')).toBe('');
  expect(await keepCommitted(page,'losatCommitted')).toBe(true);
  await page.evaluate(async()=>{
    const {state}=await import('./js/state.js');state.activeDrawing().losat.blastn.task='blastn';await window.Vue.nextTick();
  });
  expect(await committedIs(page,'losatCommitted')).toBe(true);
  await page.evaluate(async()=>{
    const {state}=await import('./js/state.js');state.activeDrawing().losat.blastn.task='megablast';await window.Vue.nextTick();
  });
  expect(await committedIs(page,'losatCommitted')).toBe(true);
});


test('S04 Linear Generate cancel, stale completion, failure and retry retain the committed Result and History', async ({ browser }, info) => {
  test.setTimeout(240000);
  const { load, generate, seeds } = require('./helpers/mode-transition.cjs');
  const page=await load(browser,seeds.linear);
  try {
    await generate(page);
    const before=await page.evaluate(async()=>{
      const app=window.__GBDRAW_APP__;
      app.adv.scale_font_size=19;
      window.s04PreviousResults=app.results;
      window.s04PreviousBasis=(await import('./js/services/config.js')).getCommittedCanonicalRenderRequest();
      return [window.__GBDRAW_HISTORY__.getUndoCount(),window.__GBDRAW_HISTORY__.getRedoCount()];
    });
    const preserved=async()=>{
      await expect(page.locator('[data-generation-application-feedback]')).toHaveCount(0);
      expect(await page.evaluate(()=>window.__GBDRAW_APP__.results===window.s04PreviousResults)).toBe(true);
      expect(await page.evaluate(()=>[window.__GBDRAW_HISTORY__.getUndoCount(),window.__GBDRAW_HISTORY__.getRedoCount()])).toEqual(before);
      expect(await committedIs(page,'s04PreviousBasis')).toBe(true);
    };
    await page.evaluate(()=>{
      window.__GBDRAW_TEST_HOOKS__.beforeDiagramGenerationResponse=()=>new Promise(resolve=>{window.s04ReleaseStale=resolve;});
    });
    await page.getByRole('button',{name:'Generate Diagram',exact:true}).click();
    await page.waitForFunction(()=>Boolean(window.s04ReleaseStale));
    await page.getByRole('button',{name:/Cancel$/}).click();
    await expect.poll(()=>page.evaluate(()=>window.__GBDRAW_APP__.processing)).toBe(false);
    await page.evaluate(()=>{delete window.__GBDRAW_TEST_HOOKS__.beforeDiagramGenerationResponse;window.s04ReleaseStale();});
    await page.evaluate(async()=>{await window.Vue.nextTick();});
    await preserved();
    await page.evaluate(()=>{
      window.__GBDRAW_TEST_HOOKS__.beforeDiagramGenerationResponse=()=>{throw new Error('S04 forced Generate failure');};
    });
    await page.getByRole('button',{name:'Generate Diagram',exact:true}).click();
    await expect.poll(()=>page.evaluate(()=>window.__GBDRAW_APP__.processing),{timeout:180000}).toBe(false);
    await expect(page.getByRole('alert',{name:'Generation Error'})).toContainText('without recognized diagnostic information');
    expect(await page.evaluate(()=>window.__GBDRAW_APP__.errorLog)).toMatchObject({
      code: 'UNKNOWN', operation: 'generate', stage: 'render'
    });
    await expect(page.getByRole('alert',{name:'Generation Error'})).not.toContainText('S04 forced Generate failure');
    await preserved();
    await page.screenshot({path:info.outputPath('linear-generate-failure-preserved.png')});
    await page.evaluate(()=>{delete window.__GBDRAW_TEST_HOOKS__.beforeDiagramGenerationResponse;});
    await generate(page);
    expect(await committedIs(page,'s04PreviousBasis')).toBe(false);
    expect(page.externalRequests).toEqual([]);
  } finally {await page.context().close();}
});


test('OV-16 a Generate that fails before it replaces the loaded Session Result settles and keeps that Result ready', async ({ browser }, info) => {
  test.setTimeout(240000);
  const { load, generate } = require('./helpers/mode-transition.cjs');
  const page = await load(browser);
  try {
    const before = await page.evaluate(() => {
      window.ov16LoadedResults = window.__GBDRAW_APP__.results;
      const observe = window.__GBDRAW_TEST_HOOKS__.onSessionLifecycleEvent;
      // The first Generate step throws once, before any candidate replaces the Result.
      window.__GBDRAW_TEST_HOOKS__.onSessionLifecycleEvent = (event) => {
        observe(event);
        if (event.name !== 'generation-input-resolution-start') return;
        window.__GBDRAW_TEST_HOOKS__.onSessionLifecycleEvent = observe;
        throw new Error('OV-16 forced early Generate failure');
      };
      return [window.__GBDRAW_HISTORY__.getUndoCount(), window.__GBDRAW_HISTORY__.getRedoCount()];
    });
    await page.getByRole('button', { name: 'Generate Diagram', exact: true }).click();
    await expect.poll(() => page.evaluate(() => window.__GBDRAW_APP__.processing), { timeout: 30000 }).toBe(false);
    await expect(page.getByRole('alert', { name: 'Generation Error' })).toBeVisible();
    await expect(page.getByRole('alert', { name: 'Generation Error' })).not.toContainText('OV-16 forced');
    expect(await page.evaluate(() => window.__GBDRAW_APP__.errorLog)).toMatchObject({ code: 'UNKNOWN', operation: 'generate' });
    expect(await page.evaluate(() => ({
      sameResult: window.__GBDRAW_APP__.results === window.ov16LoadedResults,
      preserved: window.__GBDRAW_APP__.failedGeneratePreservedResult,
      history: [window.__GBDRAW_HISTORY__.getUndoCount(), window.__GBDRAW_HISTORY__.getRedoCount()],
      restoring: window.__GBDRAW_HISTORY__.restoring.value,
      restoredReady: window.__MODE_EVENTS__.some(event => event.name === 'preview.restore-ready-receipt-accepted')
    }))).toEqual({ sameResult: true, preserved: true, history: before, restoring: false, restoredReady: true });
    await page.screenshot({ path: info.outputPath('early-generate-failure-settled.png') });
    await generate(page);
    expect(page.externalRequests).toEqual([]);
  } finally { await page.context().close(); }
});
