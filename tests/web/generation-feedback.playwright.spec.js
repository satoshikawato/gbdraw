const { test, expect } = require('@playwright/test');
const { readFileSync, writeFileSync } = require('node:fs');
const { join } = require('node:path');
const { openApp } = require('./helpers/app-lifecycle.cjs');

const source = readFileSync(join(__dirname, '../test_inputs/HmmtDNA.gbk'), 'utf8');
const applicationStatus = page => page.evaluate(async () =>
  (await import('./js/services/config.js')).getGenerationApplicationStatus());
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
  await expect(status(page)).toHaveText('Preparing diagram runtime (first use)...');
  await expect(generate(page)).toBeDisabled();
  await screenshot(page, testInfo, 'cold-runtime-visible');
  await awaitSuccess(page);
  const cold = await record(page, testInfo, 'cold');
  expectRenderStages(cold);
  expect(cold.workers.constructions).toBe(1);
  expect((await applicationStatus(page)).status).toBe('clean');

  await resetObservations(page);
  await generate(page).click();
  await awaitSuccess(page);
  const warm = await record(page, testInfo, 'warm');
  expectRenderStages(warm);
  expect(warm.statuses.map(entry => entry.value)).not.toContain('Preparing diagram runtime (first use)...');
  expect(warm.workers.constructions).toBe(1);

  await page.evaluate(async () => {
    const app = window.__GBDRAW_APP__;
    window.successfulResults = app.results;
    app.adv.scale_interval = 12345;
    window.feedbackHistory = [window.__GBDRAW_HISTORY__.getUndoCount(), window.__GBDRAW_HISTORY__.getRedoCount()];
    const config = await import('./js/services/config.js');
    window.appliedBeforeCancel = config.canonicalRenderArtifactOwner.capture().appliedGenerationIntent;
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
  expect((await applicationStatus(page)).status).toBe('pending');
  await expect(page.locator('[data-generation-application-feedback] strong')).toHaveText('Pending');
  expect(await page.evaluate(() => [window.__GBDRAW_HISTORY__.getUndoCount(), window.__GBDRAW_HISTORY__.getRedoCount()])).toEqual(await page.evaluate(() => window.feedbackHistory));
  expect(await page.evaluate(async () => window.appliedBeforeCancel ===
    (await import('./js/services/config.js')).canonicalRenderArtifactOwner.capture().appliedGenerationIntent)).toBe(true);

  await setInput('invalid GenBank');
  await generate(page).click();
  await expect.poll(() => page.evaluate(() => window.__GBDRAW_APP__.processing), { timeout: 180_000 }).toBe(false);
  await expect(page.getByRole('alert', { name: 'Generation Error' })).toBeVisible();
  await expect(generate(page)).toBeEnabled();
  expect(await page.evaluate(() => window.__GBDRAW_APP__.results === window.successfulResults)).toBe(true);
  await record(page, testInfo, 'error');
  expect((await applicationStatus(page)).status).toBe('pending');
  await expect(page.locator('[data-generation-application-feedback] strong')).toHaveText('Pending');
  expect(await page.evaluate(() => [window.__GBDRAW_HISTORY__.getUndoCount(), window.__GBDRAW_HISTORY__.getRedoCount()])).toEqual(await page.evaluate(() => window.feedbackHistory));
  expect(await page.evaluate(async () => window.appliedBeforeCancel ===
    (await import('./js/services/config.js')).canonicalRenderArtifactOwner.capture().appliedGenerationIntent)).toBe(true);

  await setInput(source);
  await page.setViewportSize({ width: 390, height: 844 });
  await resetObservations(page);
  await generate(page).click();
  await awaitSuccess(page);
  expectRenderStages(await record(page, testInfo, 'retry-mobile'));
  expect((await applicationStatus(page)).status).toBe('clean');
  await screenshot(page, testInfo, 'ready-mobile');
  expect(external).toEqual([]);
});

test('Linear comparison preparation, real search counts, cache reuse, and no-comparison generation use the same status', async ({ page, baseURL }, testInfo) => {
  test.setTimeout(180_000);
  const external = await prepare(page, baseURL);
  await page.evaluate(async content => {
    const app = window.__GBDRAW_APP__;
    app.mode = 'linear';
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
    'Running LOSAT: 0/1 source jobs complete', 'Running LOSAT: 1/1 source jobs complete',
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


test('S03 Linear LOSAT intent follows successful Result and task reversal', async ({ page, baseURL }) => {
  test.setTimeout(180000);
  await prepare(page,baseURL);
  await page.evaluate(async content => {
    const app=window.__GBDRAW_APP__;
    app.mode='linear';await window.Vue.nextTick();
    while(app.linearSeqs.length<2) app.addLinearSeq();
    for(let i=0;i<2;i++) app.setLinearSeqPrimaryFile(i,'gb',new File([content],`record-${i}.gbk`,{type:'text/plain'}));
    await app.setLinearComparisonGlobalAction('losat');await window.Vue.nextTick();
  },source);
  await page.getByRole('combobox', { name: 'LOSAT execution', exact: true }).selectOption('serial');
  await generate(page).click();
  await expect.poll(()=>page.evaluate(()=>window.feedbackEvents.some(event=>event.name==='generate.completed')),
    {timeout:180000}).toBe(true);
  expect(await page.evaluate(()=>window.__GBDRAW_APP__.errorLog?.summary||'')).toBe('');
  expect((await applicationStatus(page)).status).toBe('clean');
  await page.evaluate(async()=>{
    const {state}=await import('./js/state.js');state.losat.blastn.task='blastn';await window.Vue.nextTick();
  });
  expect((await applicationStatus(page)).status).toBe('pending');
  await page.evaluate(async()=>{
    const {state}=await import('./js/state.js');state.losat.blastn.task='megablast';await window.Vue.nextTick();
  });
  expect((await applicationStatus(page)).status).toBe('clean');
});


test('S04 settings-only and uncertain/invalid inputs never claim Applied; effective settings are observed without extra work', async ({ browser }, info) => {
  const { load, generate, seeds } = require('./helpers/mode-transition.cjs');
  const page = await browser.newPage();
  let linear;
  try {
    await openApp(page);
    const label=page.locator('[data-generation-application-feedback] strong');
    await expect(label).toHaveText('Invalid settings');
    await page.evaluate(async content => {
      window.__GBDRAW_APP__.files.c_gb=new File([content],'HmmtDNA.gbk',{type:'text/plain'});
      await window.Vue.nextTick();
    },source);
    await expect(label).toHaveText('Not generated');
    await expect(page.locator('[data-generation-application-feedback]')).toContainText('There is no Result yet');
    linear=await load(browser,seeds.linear);
    await expect(linear.locator('[data-generation-application-feedback] strong')).toHaveText('Pending');
    await generate(linear);
    const region=linear.locator('[data-generation-application-feedback]');
    await expect(region.locator('strong')).toHaveText('Applied');
    await linear.evaluate(async () => {
      const { state } = await import('./js/state.js');
      state.adv.linear_accession_visibility='show';state.adv.linear_length_visibility='show';
      // Changing a disabled slot and the inactive Circular profile is not a generation difference.
      const {createDefaultAdv}=await import('./js/services/session-active-config-contract.js');
      state.adv.circular_track_slots.push({...createDefaultAdv('linear').circular_track_slots[0],id:'inactive',enabled:false,height:'90px'});
      await window.Vue.nextTick();
    });
    await expect(region.locator('strong')).toHaveText('Applied');
    await linear.evaluate(async () => {
      const { state } = await import('./js/state.js');
      state.form.keep_definition_left_aligned='invalid';await window.Vue.nextTick();
    });
    await expect(region.locator('strong')).toHaveText('Invalid settings');
    await expect(region).toContainText('Check the generation inputs');
    await linear.evaluate(async () => {
      const { state } = await import('./js/state.js');state.form.keep_definition_left_aligned=true;
      const file=state.linearSeqs[0].gb;
      const {readFileBytes}=await import('./js/services/file-content-cache.js');
      state.linearSeqs[0].gb=new File([await readFileBytes(file)],file.name,{type:file.type});
      await window.Vue.nextTick();
    });
    await expect(region.locator('strong')).toHaveText('Pending');
    await expect.poll(() => applicationStatus(linear)).toMatchObject({status:'pending'});
    await linear.screenshot({path:info.outputPath('same-name-source-pending.png')});
  } finally { await page.close();if(linear) await linear.context().close(); }
});

test('S04 Linear Generate cancel, stale completion, failure and retry retain the Pending Result and History', async ({ browser }, info) => {
  test.setTimeout(240000);
  const { load, generate, seeds } = require('./helpers/mode-transition.cjs');
  const page=await load(browser,seeds.linear);
  try {
    await generate(page);
    const before=await page.evaluate(async()=>{
      const app=window.__GBDRAW_APP__;
      app.adv.scale_font_size=19;
      window.s04PreviousResults=app.results;
      window.s04PreviousBasis=(await import('./js/services/config.js')).canonicalRenderArtifactOwner.capture().appliedGenerationIntent;
      return [window.__GBDRAW_HISTORY__.getUndoCount(),window.__GBDRAW_HISTORY__.getRedoCount()];
    });
    const preserved=async()=>{
      await expect(page.locator('[data-generation-application-feedback] strong')).toHaveText('Pending');
      expect(await page.evaluate(()=>window.__GBDRAW_APP__.results===window.s04PreviousResults)).toBe(true);
      expect(await page.evaluate(()=>[window.__GBDRAW_HISTORY__.getUndoCount(),window.__GBDRAW_HISTORY__.getRedoCount()])).toEqual(before);
      expect(await page.evaluate(async()=>window.s04PreviousBasis===(await import('./js/services/config.js')).canonicalRenderArtifactOwner.capture().appliedGenerationIntent)).toBe(true);
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
    await expect(page.getByRole('alert',{name:'Generation Error'})).toContainText('S04 forced Generate failure');
    await preserved();
    await page.screenshot({path:info.outputPath('linear-generate-failure-pending.png')});
    await page.evaluate(()=>{delete window.__GBDRAW_TEST_HOOKS__.beforeDiagramGenerationResponse;});
    await generate(page);
    await expect(page.locator('[data-generation-application-feedback] strong')).toHaveText('Applied');
    expect(page.externalRequests).toEqual([]);
  } finally {await page.context().close();}
});
