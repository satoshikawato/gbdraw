const { test, expect } = require('@playwright/test');
const { openApp, getDiagramWorkerActivity, evaluateWithRetainedPromise } = require('./helpers/app-lifecycle.cjs');

test('real Python adapters retain causes through one lazy Worker and successful retry', async ({ page }) => {
  test.setTimeout(180000);
  const consoleMessages = [];
  page.on('console', message => consoleMessages.push(message.text()));
  page.on('pageerror', error => consoleMessages.push(error.message));
  await openApp(page);
  expect((await getDiagramWorkerActivity(page)).constructions).toBe(0);
  const observations = await page.evaluate(async () => {
    const service = await import('/gbdraw/web/js/services/diagram-generation.js');
    const { normalizeUserFacingError } = await import('/gbdraw/web/js/utils/error-normalization.js');
    const pattern = '😀[PRIVATE_PATTERN_SENTINEL';
    const failures = [];
    for (const kind of ['color', 'label']) {
      for (const features of [[], [{ type: 'CDS', qualifiers: { product: ['unrelated'] }, selector: {}, record: 'PRIVATE_RECORD_SENTINEL' }]]) {
        const rules = kind === 'color' ? [{ feat: 'CDS', qual: 'product', val: pattern }]
          : [{ recordId: '*', featureType: 'CDS', qualifier: 'product', valueRegex: pattern }];
        try {
          await service.runDiagramHelperOperation(service.DIAGRAM_HELPER_OPERATIONS.EVALUATE_RULES, { kind, features, rules });
          throw new Error('Syntax failure was accepted');
        } catch (error) {
          const model = normalizeUserFacingError(error);
          failures.push({ kind, model, again: normalizeUserFacingError(model) });
        }
      }
    }
    const session = await (await fetch('/gbdraw/web/gallery/sessions/HmmtDNA_basic_circular.gbdraw-session.json')).json();
    const request = structuredClone(session.renderRequest);
    request.diagramOptions.colors.colorTable = { resourceId: 'error-colors', representation: 'canonicalTsv' };
    const resource = (value) => {
      const bytes = new TextEncoder().encode(`feature_type\tqualifier_key\tvalue\tcolor\tcaption\nCDS\tproduct\t${value}\t#ff0000\t\n`);
      return { kind: 'canonical-tsv', name: 'PRIVATE_FILE_SENTINEL.tsv', encoding: 'base64',
        type: 'text/tab-separated-values', size: bytes.length, data: btoa(String.fromCharCode(...bytes)) };
    };
    const failed = await service.runDiagramGeneration({ request,
      resources: { ...session.resources, 'error-colors': resource(pattern) } });
    const render = normalizeUserFacingError(failed.results.error);
    const successful = await service.runDiagramGeneration({ request,
      resources: { ...session.resources, 'error-colors': resource('(?i)NADH') } });
    const helperRetry = await service.runDiagramHelperOperation(service.DIAGRAM_HELPER_OPERATIONS.EVALUATE_RULES, {
      kind: 'color', features: [], rules: [{ feat: 'CDS', qual: 'product', val: '(?P<enzyme>NADH)' }]
    });
    return { failures, render, success: successful.results.length,
      hasCatalog: Boolean(successful.metadata.featureCatalog), helperRetry: helperRetry.result };
  });
  for (const { kind, model, again } of observations.failures) {
    expect(model.code).toBe('REGEX_SYNTAX');
    expect(model.operation).toBe('evaluateRules');
    expect(model.stage).toBe('rule-validation');
    expect(model.context.position).toBe(1);
    expect(model.context.positionUnit).toBe('python-character');
    if (kind === 'label') expect(model.context.row).toBe(1);
    expect(again).toEqual(model);
  }
  expect(observations.render.code).toBe('REGEX_SYNTAX');
  expect(observations.render.operation).toBe('generate');
  expect(observations.render.stage).toBe('rule-validation');
  expect(observations.render.context.position).toBe(1);
  expect(observations.success).toBeGreaterThan(0);
  expect(observations.hasCatalog).toBe(true);
  expect(observations.helperRetry.winners).toEqual([]);
  expect(JSON.stringify(observations)).not.toContain('PRIVATE_');
  expect(consoleMessages.join('\n')).not.toContain('PRIVATE_');
  const activity = await getDiagramWorkerActivity(page);
  expect(activity.constructions).toBe(1);
  expect(activity.instances[0].initializations).toBe(1);
});

// UI-07: Generate with no input names the missing input before the diagram
// Worker (Pyodide) starts, in both modes, and Retry repeats the same check.
test('Generate without input fails before the diagram Worker starts', async ({ page }) => {
  await openApp(page);
  for (const mode of ['circular', 'linear']) {
    await page.evaluate((target) => window.__GBDRAW_APP__.setDiagramMode(target), mode);
    for (let attempt = 0; attempt < 2; attempt += 1) {
      const outcome = await evaluateWithRetainedPromise(page, () => window.__GBDRAW_APP__.runAnalysis());
      expect(outcome).toMatchObject({ status: 'error', error: { code: 'INPUT_REQUIRED', operation: 'generate' },
        recovery: 'no-result' });
    }
    await expect(page.getByRole('alert', { name: 'Generation Error' })).toContainText('Supply GenBank input');
  }
  expect((await getDiagramWorkerActivity(page)).constructions).toBe(0);
});

const { join } = require('node:path');
const { pathToFileURL } = require('node:url');
const { failNextNativeRender, inspectSafeDetails } = require('./helpers/operation-error.cjs');
const { generateAndWaitForResult } = require('./helpers/app-lifecycle.cjs');
for (const mode of ['circular', 'linear']) {
  for (const width of [1600, 390]) {
    test(`@pr-smoke native Generate failure keeps ${mode} Result and accessible diagnostics at ${width}px`, async ({ page }, info) => {
      test.setTimeout(180000);
      const logs = [];
      page.on('console', message => logs.push(message.text()));
      page.on('pageerror', error => logs.push(error.message));
      page.on('dialog', dialog => dialog.accept());
      await page.setViewportSize({ width, height: 900 });
      await openApp(page);
      if (mode === 'circular') {
        const first = await evaluateWithRetainedPromise(page, () => window.__GBDRAW_APP__.runAnalysis());
        expect(first).toMatchObject({status:'error',error:{code:'INPUT_REQUIRED',operation:'generate'},recovery:'no-result'});
        await expect(page.getByRole('alert', {name:'Generation Error'})).toContainText('No successful Result is available yet');
      }
      const session = mode === 'circular' ? 'HmmtDNA_basic_circular' : 'BGC0000708-BGC0000713';
      await page.locator('input[accept^=".json,"]').setInputFiles(join(process.cwd(), `gbdraw/web/gallery/sessions/${session}.gbdraw-session.json`));
      await page.waitForFunction(() => !window.__GBDRAW_APP__.sessionImportPending && window.__GBDRAW_APP__.results.length);
      const snapshot = () => page.evaluate(async () => {
        const { getCommittedCanonicalSession } = await import('./js/services/config.js');
        const a = window.__GBDRAW_APP__, h = window.__GBDRAW_HISTORY__;
        return { results: JSON.stringify(a.results), session: JSON.stringify(getCommittedCanonicalSession()),
          undo: h.getUndoCount(), redo: h.getRedoCount() };
      });
      const before = await snapshot();
      await failNextNativeRender(page);
      const failed = await evaluateWithRetainedPromise(page, () => window.__GBDRAW_APP__.runAnalysis());
      expect(failed).toMatchObject({status:'error',error:{code:'REGEX_SYNTAX',operation:'generate',stage:'rule-validation',
        context:{position:1,positionUnit:'python-character',reason:'UNTERMINATED_SET'}},recovery:'preserved'});
      expect(await snapshot()).toEqual(before);
      const alert = page.getByRole('alert', {name:'Generation Error'});
      await expect(alert).toContainText('last successful Result and committed request are unchanged');
      await inspectSafeDetails(page, alert, failed.error);
      await alert.screenshot({path:info.outputPath('operation-error-alert.png')});
      await page.locator('.preview-feature-search').screenshot({path:info.outputPath('search-controls.png')});
      await page.screenshot({path:info.outputPath('operation-error.png'),fullPage:true});
      expect(logs.join('\n')).not.toContain('PRIVATE_');
      const retry = alert.getByRole('button', {name:'Retry Generate',exact:true});
      await retry.scrollIntoViewIfNeeded(); await retry.focus(); await retry.press('Enter');
      await page.waitForFunction(() => !window.__GBDRAW_APP__.processing, null, {timeout:180000});
      await expect(alert).toHaveCount(0);
      expect((await snapshot()).undo).toBe(before.undo + 1);
      expect(logs.join('\n')).not.toContain('PRIVATE_');
    });
  }
}

// Replace the next Worker run's first record display start (an engine-side check).
const setNextDisplayStart = (page, startCoordinate) => page.evaluate((start) => {
  const send = Worker.prototype.postMessage;
  Worker.prototype.postMessage = function (message, ...args) {
    if (message.type !== 'run') return send.call(this, message, ...args);
    Worker.prototype.postMessage = send;
    const request = structuredClone(message.payload.request);
    request.records[0].display.startCoordinate = start;
    return send.call(this, { ...message, payload: { ...message.payload, request } }, ...args);
  };
}, startCoordinate);
const loadSession = async (page, name) => {
  await page.locator('input[accept^=".json,"]').setInputFiles(join(process.cwd(), `gbdraw/web/gallery/sessions/${name}.gbdraw-session.json`));
  await page.waitForFunction(() => !window.__GBDRAW_APP__.sessionImportPending && window.__GBDRAW_APP__.results.length);
};

// B9 (P07: offer only working actions): a render-stage failure that no
// producer classifies repeats for the same inputs, so the panel names the
// render failure and its exception class and offers Save Session, not Retry.
test('an unclassified render failure offers Save Session and no Retry Generate', async ({ page }) => {
  test.setTimeout(180000);
  await openApp(page);
  await loadSession(page, 'BGC0000708-BGC0000713');
  // A display start on a Linear record reaches the engine; its native check has no diagnostic.
  await setNextDisplayStart(page, 1);
  const failed = await evaluateWithRetainedPromise(page, () => window.__GBDRAW_APP__.runAnalysis());
  expect(failed).toMatchObject({ status: 'error', recovery: 'preserved', error: { code: 'RENDER_FAILED',
    operation: 'generate', stage: 'render', context: { exceptionType: 'ValidationError' }, actions: ['save-session'] } });
  const alert = page.getByRole('alert', { name: 'Generation Error' });
  await expect(alert).toContainText('The diagram engine failed while drawing this diagram.');
  await expect(alert).toContainText('Python exception: ValidationError.');
  await expect(alert).not.toContainText('Input validation failed');
  await expect(alert.getByRole('button', { name: 'Retry Generate', exact: true })).toHaveCount(0);
  await expect(alert.getByRole('button', { name: 'Save Session', exact: true })).toBeVisible();
  await inspectSafeDetails(page, alert, failed.error);
});

// OV-130 (R6): a preserved setting of the other diagram mode (OV-106) fails the
// engine's request decoding. The panel names the setting, the mode it belongs to
// and the two ways out, instead of the unclassified input-validation text.
test('a Circular setting in a Linear request names the setting and the ways out', async ({ page }) => {
  test.setTimeout(180000);
  await openApp(page);
  await loadSession(page, 'BGC0000708-BGC0000713');
  await page.evaluate(() => {
    const send = Worker.prototype.postMessage;
    Worker.prototype.postMessage = function (message, ...args) {
      if (message.type !== 'run') return send.call(this, message, ...args);
      Worker.prototype.postMessage = send;
      const request = structuredClone(message.payload.request);
      request.diagramOptions.configOverrides = { ...request.diagramOptions.configOverrides, 'objects.ticks.tick_width': 4 };
      return send.call(this, { ...message, payload: { ...message.payload, request } }, ...args);
    };
  });
  const failed = await evaluateWithRetainedPromise(page, () => window.__GBDRAW_APP__.runAnalysis());
  expect(failed).toMatchObject({ status: 'error', recovery: 'preserved', error: { code: 'MODE_SETTING',
    operation: 'generate', stage: 'request-validation', actions: ['edit-input'],
    context: { reason: 'CIRCULAR_SETTING', configPath: 'objects.ticks.tick_width' } } });
  const alert = page.getByRole('alert', { name: 'Generation Error' });
  await expect(alert).toContainText('A preserved session setting does not apply to this diagram mode.'
    + ' Setting: objects.ticks.tick_width. It applies only to Circular diagrams.'
    + ' Reset it under Preserved session settings, or switch to Circular.');
  await expect(alert).not.toContainText('Input validation failed');
  await expect(alert.getByRole('button', { name: 'Retry Generate', exact: true })).toHaveCount(0);
  await inspectSafeDetails(page, alert, failed.error);
});

// B10: the engine's display-start range check is a user-fixable input error.
test('a display start beyond the record is an input error with a working action', async ({ page }) => {
  test.setTimeout(180000);
  await openApp(page);
  await loadSession(page, 'HmmtDNA_basic_circular');
  await setNextDisplayStart(page, 99999999);
  const failed = await evaluateWithRetainedPromise(page, () => window.__GBDRAW_APP__.runAnalysis());
  expect(failed).toMatchObject({ status: 'error', recovery: 'preserved', error: { code: 'INPUT_INVALID',
    operation: 'generate', stage: 'render', context: { field: 'start', reason: 'DISPLAY_START_BOUNDS' } } });
  expect(failed.error.actions).toContain('retry');
  const alert = page.getByRole('alert', { name: 'Generation Error' });
  await expect(alert).toContainText('Use a display start between 1 and the record length.');
  await expect(alert).not.toContainText('The diagram engine failed while drawing');
  await inspectSafeDetails(page, alert, failed.error);
});

// B11: a live edit rerenders the committed Session with the current editor
// tables, so a failure that repeats for that request offers no Retry of the
// live edit; the note above the Result shows the cause and Generate.
test('a live edit that fails the same way for the same request offers no Retry', async ({ page }) => {
  test.setTimeout(180000);
  await openApp(page);
  await loadSession(page, 'BGC0000708-BGC0000713');
  const before = await page.evaluate(() => window.__GBDRAW_APP__.results[0].content);
  await setNextDisplayStart(page, 1);
  // The forced reflow request is the path a label visibility edit uses.
  await page.evaluate(async () => {
    const { state } = await import('/gbdraw/web/js/state.js');
    state.labelReflowForceRequestSeq.value += 1;
  });
  await expect.poll(() => page.evaluate(() => window.__GBDRAW_APP__.labelReflowLastError?.code ?? null),
    { timeout: 120000 }).toBe('RENDER_FAILED');
  const failed = await page.evaluate(() => JSON.parse(JSON.stringify(window.__GBDRAW_APP__.labelReflowLastError)));
  expect(failed).toMatchObject({ code: 'RENDER_FAILED', stage: 'render', context: { exceptionType: 'ValidationError' },
    actions: ['generate'] });
  const note = page.locator('[data-live-application-feedback]');
  await expect(note).toContainText('Live edit failed: direct edits already applied are kept');
  await expect(note).toContainText('The diagram engine failed while drawing this diagram.');
  await expect(note).toContainText('Change the edit, or change the settings and use Generate.');
  await expect(note).not.toContainText('Retry the live edit');
  expect(await page.evaluate(() => window.__GBDRAW_APP__.results[0].content)).toBe(before);
});

test('@pr-smoke live and downloaded standalone search retain JavaScript regex and word targets', async ({page}, info) => {
  test.setTimeout(180000);
  page.on('dialog', dialog=>dialog.accept());
  await openApp(page);
  await page.locator('input[accept^=".json,"]').setInputFiles(join(process.cwd(),'gbdraw/web/gallery/sessions/HmmtDNA_basic_circular.gbdraw-session.json'));
  await page.waitForFunction(()=>!window.__GBDRAW_APP__.sessionImportPending && window.__GBDRAW_APP__.extractedFeatures.length);
  const query=page.getByRole('searchbox',{name:'Search features',exact:true});
  await page.getByRole('combobox',{name:'Search field',exact:true}).selectOption('qualifier-value');
  await page.getByRole('textbox',{name:'Qualifier key',exact:true}).fill('product');
  const regex=page.getByRole('checkbox',{name:'Regex (JavaScript, i)',exact:true});
  await regex.check(); await query.fill('(?<enzyme>nadh)'); await query.press('Enter');
  const targets=await page.evaluate(()=>[...window.__GBDRAW_APP__.previewFeatureSearchMatches].sort());
  expect(targets).toHaveLength(7);
  await query.fill('(?P<enzyme>NADH)'); await query.press('Enter');
  await expect(page.getByRole('status',{name:'Feature search status'})).toContainText('Invalid JavaScript regular expression. Turn off Regex to return to word search.');
  await regex.uncheck(); await query.fill('NADH'); await query.press('Enter');
  await expect.poll(()=>page.evaluate(()=>[...window.__GBDRAW_APP__.previewFeatureSearchMatches].sort())).toEqual(targets);
  await page.evaluate(()=>{window.__GBDRAW_APP__.richFeaturePopup=true;});
  const pending=page.waitForEvent('download');
  await page.evaluate(()=>window.__GBDRAW_APP__.downloadInteractiveSVG());
  const download=await pending;
  const stream=await download.createReadStream();const chunks=[];
  for await(const chunk of stream)chunks.push(chunk);
  const source=Buffer.concat(chunks).toString('utf8');
  expect(source).toContain('Regex (JavaScript, i)');
  expect(source).toContain('Invalid JavaScript regular expression. Turn off Regex to return to word search.');
  const svgPath=info.outputPath('downloaded.interactive.svg');
  await download.saveAs(svgPath);
  await page.goto(pathToFileURL(svgPath).href);
  await page.getByRole('button',{name:'Expand feature search',exact:true}).click();
  const standaloneQuery=page.getByRole('searchbox',{name:'Search features',exact:true});
  await page.getByRole('combobox',{name:'Search field',exact:true}).selectOption('qualifier-value');
  await page.getByRole('textbox',{name:'Qualifier key for qualifier value search',exact:true}).fill('product');
  const standaloneRegex=page.getByRole('checkbox',{name:'Regex (JavaScript, i)',exact:true});
  await standaloneRegex.check();await standaloneQuery.fill('(?<enzyme>nadh)');await page.locator('[data-search-apply]').click();
  const selected=()=>page.locator('.gbdraw-interactive-feature--match[data-gbdraw-feature-id]').evaluateAll(elements=>[...new Set(elements.map(e=>e.getAttribute('data-gbdraw-feature-id')))].sort());
  await expect(page.locator('[data-search-count]')).toContainText('/ 7 features');
  await expect.poll(selected).toEqual(targets);
  await standaloneQuery.fill('(?P<enzyme>NADH)');await page.locator('[data-search-apply]').click();
  await expect(page.locator('[data-search-count]')).toContainText('Invalid JavaScript regular expression');
  await standaloneRegex.uncheck();await standaloneQuery.fill('NADH');await page.locator('[data-search-apply]').click();
  await expect.poll(selected).toEqual(targets);
  await page.screenshot({path:info.outputPath('standalone-search.png')});
});

test('table import failures retain exact cause and concrete retry without a draft lifecycle',async({page})=>{
  test.setTimeout(180000);page.on('dialog',dialog=>dialog.accept());
  await openApp(page);
  await page.locator('input[accept^=".json,"]').setInputFiles(join(process.cwd(),'gbdraw/web/gallery/sessions/HmmtDNA_basic_circular.gbdraw-session.json'));
  await page.waitForFunction(()=>!window.__GBDRAW_APP__.sessionImportPending && window.__GBDRAW_APP__.extractedFeatures.length);
  const before=await page.evaluate(()=>({svg:window.__GBDRAW_APP__.svgContent,rules:JSON.stringify(window.__GBDRAW_APP__.manualSpecificRules),history:window.__GBDRAW_HISTORY__.getUndoCount()}));
  const uploader=page.getByLabel('Specific Table (-t)',{exact:true});
  await uploader.setInputFiles({name:'PRIVATE_FILE_SENTINEL.tsv',mimeType:'text/plain',buffer:Buffer.from('CDS\tproduct\t(?<enzyme>NADH)\t#ff0000\t\n')});
  const alert=page.getByRole('alert',{name:'Rule error'});
  await expect(alert).toContainText('Python regular expression is invalid',{timeout:180000});
  await expect(alert.getByRole('button',{name:'Retry table import',exact:true})).toBeVisible();
  await page.evaluate(()=>{window.__S03_PREVIOUS_IMPORT_ERROR__=window.__GBDRAW_APP__.errorLog;});
  await alert.getByRole('button',{name:'Retry table import',exact:true}).click();
  await page.waitForFunction(()=>window.__GBDRAW_APP__.errorLog!==window.__S03_PREVIOUS_IMPORT_ERROR__ && !window.__GBDRAW_APP__.files.t_color);
  await expect(alert).toContainText('Python regular expression is invalid',{timeout:180000});
  expect(await page.evaluate(()=>window.__GBDRAW_APP__.svgContent)).toBe(before.svg);
  expect(await page.evaluate(()=>JSON.stringify(window.__GBDRAW_APP__.manualSpecificRules))).toBe(before.rules);
  expect(await page.evaluate(()=>window.__GBDRAW_HISTORY__.getUndoCount())).toBe(before.history);
  await uploader.setInputFiles({name:'corrected.tsv',mimeType:'text/plain',buffer:Buffer.from('CDS\tproduct\t(?P<enzyme>NADH)\t#ff0000\t\n')});
  await page.waitForFunction(()=>window.__GBDRAW_APP__.manualSpecificRules.length===1);
  await expect(alert).toHaveCount(0);
  expect(await page.evaluate(()=>window.__GBDRAW_HISTORY__.getUndoCount())).toBe(before.history+1);
});


test('Label TSV retry retains its real reselect input and accepted History', async ({page})=>{
  test.setTimeout(180000);page.on('dialog',dialog=>dialog.accept());
  await openApp(page);
  await page.locator('input[accept^=".json,"]').setInputFiles(join(process.cwd(),'gbdraw/web/gallery/sessions/HmmtDNA_basic_circular.gbdraw-session.json'));
  await page.waitForFunction(()=>!window.__GBDRAW_APP__.sessionImportPending && window.__GBDRAW_APP__.extractedFeatures.length);
  await page.evaluate(()=>window.__GBDRAW_APP__.openRightDrawerTab('features'));
  const before=await page.evaluate(()=>({svg:window.__GBDRAW_APP__.svgContent,history:window.__GBDRAW_HISTORY__.getUndoCount()}));
  // The Features list has two TSV file inputs (Label TSV, feature edits TSV); its button opens this one.
  const picker=page.waitForEvent('filechooser');
  await page.getByRole('button',{name:'Load Label TSV',exact:true}).click();
  await (await picker).setFiles({name:'PRIVATE_LABEL_SENTINEL.tsv',mimeType:'text/plain',buffer:Buffer.from('*\tCDS\tproduct\t(?<enzyme>NADH)\tPRIVATE_LABEL_SENTINEL\n')});
  const alert=page.getByRole('alert',{name:'Rule error'});
  await expect(alert).toContainText('Python regular expression is invalid',{timeout:180000});
  await page.evaluate(()=>window.__GBDRAW_APP__.retryLabelImportFailure());
  await expect(alert).toContainText('Python regular expression is invalid');
  expect(await page.evaluate(()=>window.__GBDRAW_APP__.svgContent)).toBe(before.svg);
  expect(await page.evaluate(()=>window.__GBDRAW_HISTORY__.getUndoCount())).toBe(before.history);
  const pending=page.waitForEvent('filechooser');
  await alert.getByRole('button',{name:'Reselect Label TSV',exact:true}).click();
  const chooser=await pending;
  await chooser.setFiles({name:'corrected-labels.tsv',mimeType:'text/plain',buffer:Buffer.from('*\tCDS\tproduct\t(?P<enzyme>NADH)\tCORRECTED\n')});
  await expect(page.locator('.origin-top svg text').filter({hasText:/^CORRECTED$/})).toHaveCount(7);
  await expect(alert).toHaveCount(0);
  expect(await page.evaluate(()=>window.__GBDRAW_HISTORY__.getUndoCount())).toBe(before.history+1);
});
