const { test, expect } = require('@playwright/test');
const { join } = require('node:path');
const { evaluateWithRetainedPromise, openApp, generateAndWaitForResult, getDiagramWorkerActivity } = require('./helpers/app-lifecycle.cjs');

const loadSession = async (page) => {
  await openApp(page);
  await page.locator('input[accept^=".json,"]').setInputFiles(join(process.cwd(), 'gbdraw/web/gallery/sessions/HmmtDNA_basic_circular.gbdraw-session.json'));
  await page.waitForFunction(() => !window.__GBDRAW_APP__.sessionImportPending && window.__GBDRAW_APP__.extractedFeatures.length);
};
// A rule that changes the Legend asks for the automatic rerender (OV-42, OV-43,
// #857); the Result and its Legend are Python's once it has finished.
const rerendered = (page) => page.waitForFunction(() => !window.__GBDRAW_APP__.processing
  && !window.__GBDRAW_APP__.labelReflowProcessing);
// The wait runs from a promise the page retains: V8 can collect a promise that
// page.evaluate awaits ("Promise was collected"), which Playwright reports as a
// destroyed execution context.
const waitForColorTableImport = (page) => evaluateWithRetainedPromise(page, async () => {
  const a = window.__GBDRAW_APP__;
  await a.waitForAuxiliaryFileImport(a.files.t_color);
});
const fills = (page) => page.evaluate(() => {
  const app = window.__GBDRAW_APP__;
  return app.extractedFeatures.filter((feature) => {
    const element = document.querySelector(`.origin-top svg [data-gbdraw-feature-id="${feature.svg_id}"][fill="#f01234"]`);
    return Boolean(element);
  }).map((feature) => feature.svg_id).sort();
});

test('Python-only color rules commit atomically, retain History and agree with Generate', async ({ page }) => {
  test.setTimeout(180000);
  const messages = [];
  page.on('dialog', async (dialog) => { messages.push(dialog.message()); await dialog.accept(); });
  await loadSession(page);
  expect((await getDiagramWorkerActivity(page)).constructions).toBe(0);
  messages.length = 0;
  const before = await page.evaluate(() => window.__GBDRAW_APP__.svgContent);
  await evaluateWithRetainedPromise(page, async () => {
    const a = window.__GBDRAW_APP__;
    Object.assign(a.newSpecRule, { feat: 'CDS', qual: 'product', val: '(?i)NADH', color: '#f01234', cap: '' });
    await a.addSpecificRule();
  });
  const matched = await fills(page);
  expect(matched).toHaveLength(7);
  expect(messages).toEqual([]);
  await page.evaluate(() => window.__GBDRAW_APP__.undoHistory());
  expect(await fills(page)).toEqual([]);
  await page.evaluate(() => window.__GBDRAW_APP__.redoHistory());
  expect(await fills(page)).toEqual(matched);
  const stable = await page.evaluate(() => ({ rules: JSON.stringify(window.__GBDRAW_APP__.manualSpecificRules), svg: window.__GBDRAW_APP__.svgContent }));
  await page.evaluate(() => window.__GBDRAW_APP__.setSpecificRuleField(0, 'val', '(?<enzyme>NADH)'));
  expect(await page.evaluate(() => { const a = window.__GBDRAW_APP__; return a.specificRulePatternDraft(a.manualSpecificRules[0]).error; })).toMatchObject({code:'REGEX_SYNTAX',context:{reason:'UNKNOWN_EXTENSION'}});
  expect(await page.evaluate(() => ({ rules: JSON.stringify(window.__GBDRAW_APP__.manualSpecificRules), svg: window.__GBDRAW_APP__.svgContent }))).toEqual(stable);
  await page.evaluate(() => window.__GBDRAW_APP__.setSpecificRuleField(0, 'val', '(?P<enzyme>NADH)'));
  expect(await fills(page)).toEqual(matched);
  await generateAndWaitForResult(page);
  expect(await fills(page)).toEqual(matched);
  expect(await page.evaluate(() => window.__GBDRAW_APP__.svgContent)).not.toBe(before);
});

test('Python label TSV syntax is validated before replacing live overrides', async ({ page }) => {
  test.setTimeout(180000);
  const messages = [];
  page.on('dialog', async (dialog) => { messages.push(dialog.message()); await dialog.accept(); });
  await loadSession(page);
  const importLabel = (pattern, label) => evaluateWithRetainedPromise(page, async ({ pattern, label }) => {
    const input = { files: [new File([`*\tCDS\tproduct\t${pattern}\t${label}\n`], 'labels.tsv')], value: 'labels.tsv' };
    await window.__GBDRAW_APP__.loadLabelOverrideTable({ target: input });
  }, { pattern, label });
  await importLabel('(?i)NADH', 'MATCHED');
  await expect(page.locator('[data-app-notice] p')).toHaveText('Loaded 1 row. Applied to 7 labels.');
  await page.evaluate(() => window.__GBDRAW_APP__.undoHistory());
  expect(await page.locator('.origin-top svg text').filter({ hasText: /^MATCHED$/ }).count()).toBe(0);
  await page.evaluate(() => window.__GBDRAW_APP__.redoHistory());
  expect(await page.locator('.origin-top svg text').filter({ hasText: /^MATCHED$/ }).count()).toBe(7);
  const snapshot = await page.evaluate(() => ({ svg: window.__GBDRAW_APP__.svgContent, overrides: JSON.stringify(window.__GBDRAW_APP__.featureOverrides) }));
  await importLabel('(?<enzyme>NADH)', 'INVALID');
  expect(await page.evaluate(() => window.__GBDRAW_APP__.errorLog)).toMatchObject({code:'REGEX_SYNTAX',context:{reason:'UNKNOWN_EXTENSION'}});
  expect(await page.evaluate(() => ({ svg: window.__GBDRAW_APP__.svgContent, overrides: JSON.stringify(window.__GBDRAW_APP__.featureOverrides) }))).toEqual(snapshot);
  await generateAndWaitForResult(page);
  expect(await page.locator('.origin-top svg text').filter({ hasText: /^MATCHED$/ }).count()).toBe(7);
});

test('late Python evaluation cannot replace a newer rule or a loaded Session', async ({ page }) => {
  test.setTimeout(180000);
  page.on('dialog', dialog => dialog.accept());
  await page.addInitScript(() => {
    const send = Worker.prototype.postMessage;
    Worker.prototype.postMessage = function (message, ...args) {
      if (window.__HOLD_RULES__ && message.operation === 'evaluateRules') {
        window.__HOLD_RULES__ = false;
        window.__RELEASE_RULES__ = () => send.call(this, message, ...args);
        return;
      }
      return send.call(this, message, ...args);
    };
  });
  await loadSession(page);
  await evaluateWithRetainedPromise(page, async () => {
    const a = window.__GBDRAW_APP__;
    Object.assign(a.newSpecRule, { feat: 'CDS', qual: 'product', val: '(?i)NADH', color: '#f01234', cap: '' });
    await a.addSpecificRule();
    window.__HOLD_RULES__ = true;
    window.__OLD_RULE_EDIT__ = a.setSpecificRuleField(0, 'val', '(?i)cytochrome');
  });
  await page.waitForFunction(() => window.__RELEASE_RULES__);
  await page.evaluate(() => window.__GBDRAW_APP__.setSpecificRuleField(0, 'val', '(?i)ATP'));
  const newer = await fills(page);
  expect(newer.length).toBeGreaterThan(0);
  await page.evaluate(async () => { window.__RELEASE_RULES__(); await window.__OLD_RULE_EDIT__; });
  expect(await fills(page)).toEqual(newer);
  expect(await page.evaluate(() => window.__GBDRAW_APP__.manualSpecificRules[0].val)).toBe('(?i)ATP');
  // No late evaluation crosses a Load: a rule edit's only request is its
  // match (a pattern not matched before, so the edit asks), and Load waits
  // for it (OV-377). A Load chosen while the match is held is refused and
  // changes nothing. Once the edit is applied, the loaded Session replaces
  // the rules and their fills. (The caption normalization request that a
  // Load used to take over is gone with the caption rewrite.)
  const sessionFile = join(process.cwd(), 'gbdraw/web/gallery/sessions/HmmtDNA_basic_circular.gbdraw-session.json');
  await evaluateWithRetainedPromise(page, async () => {
    const a = window.__GBDRAW_APP__;
    window.__RELEASE_RULES__ = null;
    window.__HOLD_RULES__ = true;
    window.__OLD_RULE_EDIT__ = a.setSpecificRuleField(0, 'val', '(?i)nadh');
  });
  await page.waitForFunction(() => window.__RELEASE_RULES__);
  const held = () => page.evaluate(() => {
    const a = window.__GBDRAW_APP__;
    return { available: a.sessionLoadAvailable, importing: a.sessionImportPending, rules: JSON.stringify(a.manualSpecificRules) };
  });
  const beforeLoad = await held();
  expect(beforeLoad).toMatchObject({ available: false, importing: false });
  await page.locator('input[accept^=".json,"]').setInputFiles(sessionFile);
  expect(await held()).toEqual(beforeLoad);
  await page.evaluate(async () => { window.__RELEASE_RULES__(); await window.__OLD_RULE_EDIT__; });
  expect(await page.evaluate(() => {
    const a = window.__GBDRAW_APP__;
    return [a.sessionLoadAvailable, a.manualSpecificRules.map((rule) => rule.val)];
  })).toEqual([true, ['(?i)nadh']]);
  expect((await fills(page)).length).toBeGreaterThan(0);
  await page.locator('input[accept^=".json,"]').setInputFiles(sessionFile);
  await page.waitForFunction(() => !window.__GBDRAW_APP__.sessionImportPending && !window.__GBDRAW_APP__.manualSpecificRules.length);
  expect(await fills(page)).toEqual([]);
  expect(await page.evaluate(() => window.__GBDRAW_APP__.manualSpecificRules)).toEqual([]);
});

test('Session preview stays lazy and saved Python rules remain editable and regenerable', async ({ page, browser }) => {
  test.setTimeout(180000);
  page.on('dialog', dialog => dialog.accept());
  await loadSession(page);
  await page.getByLabel('Specific rules file', { exact: true }).setInputFiles({ name: 'python.tsv', mimeType: 'text/plain', buffer: Buffer.from('CDS\tproduct\t(?P<enzyme>NADH)\t#f01234\t\n') });
  await page.evaluate(async () => { const a = window.__GBDRAW_APP__; await a.waitForAuxiliaryFileImport(a.files.t_color); });
  expect(await fills(page)).toHaveLength(7);
  const download = page.waitForEvent('download');
  await evaluateWithRetainedPromise(page, () => window.__GBDRAW_APP__.saveSessionWithTitle());
  const file = await (await download).path();
  const context = await browser.newContext();
  const fresh = await context.newPage();
  fresh.on('dialog', dialog => dialog.accept());
  await openApp(fresh);
  await fresh.locator('input[accept^=".json,"]').setInputFiles(file);
  await fresh.waitForFunction(() => !window.__GBDRAW_APP__.sessionImportPending && window.__GBDRAW_APP__.extractedFeatures.length);
  expect((await getDiagramWorkerActivity(fresh)).constructions).toBe(0);
  expect(await fills(fresh)).toHaveLength(7);
  await fresh.evaluate(() => window.__GBDRAW_APP__.setSpecificRuleField(0, 'color', '#abcdef'));
  await generateAndWaitForResult(fresh);
  expect(await fresh.locator('.origin-top svg [fill="#abcdef"]').count()).toBeGreaterThanOrEqual(7);
  await context.close();
});

test('Unicode regex corpus has identical live and generated color/label targets', async ({ page }) => {
  test.setTimeout(180000);
  page.on('dialog', dialog => dialog.accept());
  await openApp(page);
  await page.getByLabel('GenBank/DDBJ File', { exact: true }).setInputFiles(join(process.cwd(), 'tests/fixtures/regex_rules.gb'));
  await page.evaluate(() => { window.__GBDRAW_APP__.form.labels_mode = 'out'; window.__GBDRAW_APP__.autoLabelReflowEnabled = false; });
  await generateAndWaitForResult(page);
  for (const [pattern, expected] of [
    ['(?i)NADH', ['R0', 'R1']], ['(?P<enzyme>NADH)', ['R0', 'R1']],
    ['NADH\\Z', ['R0']], ['\\bβ', ['R2']], ['i', ['R3', 'R4', 'R5']]
  ]) {
    await evaluateWithRetainedPromise(page, async (pattern) => {
      const a = window.__GBDRAW_APP__;
      await a.clearAllSpecificRules();
      Object.assign(a.newSpecRule, { feat: 'CDS', qual: 'product', val: pattern, color: '#f01234', cap: '' });
      await a.addSpecificRule();
    }, pattern);
    const ids = await fills(page);
    const tags = await page.evaluate(ids => window.__GBDRAW_APP__.extractedFeatures.filter(f => ids.includes(f.svg_id)).map(f => f.locus_tag).sort(), ids);
    expect(tags).toEqual(expected);
    await evaluateWithRetainedPromise(page, async (pattern) => {
      await window.__GBDRAW_APP__.loadLabelOverrideTable({ target: { files: [new File([`*\tCDS\tproduct\t${pattern}\tMATCH\n`], 'labels.tsv')], value: 'labels.tsv' } });
    }, pattern);
    expect(await page.locator('.origin-top svg text').filter({ hasText: /^MATCH$/ }).count()).toBe(expected.length);
    await generateAndWaitForResult(page);
    expect(await fills(page)).toEqual(ids);
    expect(await page.locator('.origin-top svg text').filter({ hasText: /^MATCH$/ }).count()).toBe(expected.length);
  }
  await page.evaluate(() => window.__GBDRAW_APP__.clearAllSpecificRules());
  expect(await fills(page)).toEqual([]);
});

// The presets are the built-in color tables (gbdraw/data/color_tables/, `-t
// builtin:<id>` on the CLI): applying one writes its rows, and Undo, Redo and
// Generate draw what the rows color.
const presetFills = (page, color) => page.evaluate((fill) => window.__GBDRAW_APP__.extractedFeatures
  .filter((feature) => document.querySelector(`.origin-top svg [data-gbdraw-feature-id="${feature.svg_id}"][fill="${fill}"]`))
  .map((feature) => feature.svg_id).sort(), color);
const legendTexts = (page) => page.evaluate(() => [...document.querySelectorAll('.origin-top svg text')].map((node) => node.textContent));

test('the mitochondrion preset updates the legend; Undo, Redo and Generate agree', async ({ page }) => {
  test.setTimeout(180000);
  await loadSession(page);
  const preset = await page.evaluate(() => window.__GBDRAW_APP__.specificRulePresets.find((entry) => entry.id === 'mitochondrion'));
  const cox = preset.rules.find((rule) => rule.cap === 'complex IV (cytochrome c oxidase)');
  await evaluateWithRetainedPromise(page, async () => {
    const a = window.__GBDRAW_APP__;
    a.selectedSpecificPreset = 'mitochondrion';
    await a.applySpecificRulePreset();
  });
  await rerendered(page);
  expect(await page.evaluate(() => window.__GBDRAW_APP__.manualSpecificRules.length)).toBe(preset.rules.length);
  const matched = await presetFills(page, cox.color);
  expect(matched).toHaveLength(3);
  expect(await legendTexts(page)).toContain(cox.cap);
  await page.evaluate(() => window.__GBDRAW_APP__.undoHistory());
  expect(await presetFills(page, cox.color)).toEqual([]);
  await page.evaluate(() => window.__GBDRAW_APP__.redoHistory());
  await rerendered(page);
  expect(await presetFills(page, cox.color)).toEqual(matched);
  await generateAndWaitForResult(page);
  expect(await presetFills(page, cox.color)).toEqual(matched);
  expect(await legendTexts(page)).toContain(cox.cap);
});

test('the plastid preset replaces the tobacco Gallery rules through the Preset scheme select', async ({ page }) => {
  test.setTimeout(180000);
  await openApp(page);
  await page.locator('input[accept^=".json,"]').setInputFiles(join(process.cwd(), 'gbdraw/web/gallery/sessions/tobacco-chloroplast.gbdraw-session.json'));
  await page.waitForFunction(() => !window.__GBDRAW_APP__.sessionImportPending && window.__GBDRAW_APP__.extractedFeatures.length);
  const preset = await page.evaluate(() => window.__GBDRAW_APP__.specificRulePresets.find((entry) => entry.id === 'plastid'));
  const select = page.getByLabel('Preset scheme', { exact: true });
  await select.evaluate((element) => { for (let node = element.closest('details'); node; node = node.parentElement?.closest('details')) node.open = true; });
  // The Session already holds the plastid rows, so apply another preset first: applying the same rows again changes nothing to draw.
  for (const id of ['mitochondrion', 'plastid']) {
    const rows = await page.evaluate((presetId) => window.__GBDRAW_APP__.specificRulePresets.find((entry) => entry.id === presetId).rules, id);
    await select.selectOption(id);
    await select.locator('xpath=..').getByRole('button', { name: 'Apply' }).click();
    await page.waitForFunction((first) => window.__GBDRAW_APP__.manualSpecificRules[0]?.cap === first, rows[0].cap);
    await rerendered(page);
  }
  expect(await page.evaluate(() => window.__GBDRAW_APP__.manualSpecificRules.length)).toBe(preset.rules.length);
  await rerendered(page);
  const psi = preset.rules.find((rule) => rule.cap === 'photosystem I');
  const live = await presetFills(page, psi.color);
  expect(live.length).toBeGreaterThan(0);
  expect(await legendTexts(page)).toEqual(expect.arrayContaining(['photosystem I', 'transfer RNAs', 'origin of replication']));
  await generateAndWaitForResult(page);
  expect(await presetFills(page, psi.color)).toEqual(live);
});

test('invalid color TSV preserves the selected file, preview, rules and History', async ({ page }) => {
  test.setTimeout(180000);
  const messages = [];
  page.on('dialog', async dialog => { messages.push(dialog.message()); await dialog.accept(); });
  await loadSession(page);
  const upload = page.getByLabel('Specific rules file', { exact: true });
  await upload.setInputFiles({ name: 'valid.tsv', mimeType: 'text/plain', buffer: Buffer.from('CDS\tproduct\t(?i)NADH\t#f01234\t\n') });
  await waitForColorTableImport(page);
  await expect.poll(() => page.evaluate(() => window.__GBDRAW_HISTORY__.getUndoCount())).toBe(1);
  const inspect = () => page.evaluate(() => {
    const a = window.__GBDRAW_APP__;
    return { file: a.files.t_color?.name, rules: JSON.stringify(a.manualSpecificRules), svg: a.svgContent, history: window.__GBDRAW_HISTORY__.getUndoCount() };
  });
  const before = await inspect();
  await upload.setInputFiles({ name: 'invalid.tsv', mimeType: 'text/plain', buffer: Buffer.from('CDS\tproduct\t(?<enzyme>NADH)\t#abcdef\tInvalid\n') });
  await waitForColorTableImport(page);
  expect(await page.evaluate(() => window.__GBDRAW_APP__.errorLog)).toMatchObject({code:'REGEX_SYNTAX',context:{reason:'UNKNOWN_EXTENSION'}});
  await expect.poll(inspect).toEqual(before);
  await generateAndWaitForResult(page);
  expect(await fills(page)).toHaveLength(7);
});

test('25,000-feature Python preparation keeps the event loop responsive and reuses one worker', async ({ page }) => {
  test.setTimeout(180000);
  await openApp(page);
  expect((await getDiagramWorkerActivity(page)).constructions).toBe(0);
  const result = await page.evaluate(async () => {
    const { createRulePreparation } = await import('/gbdraw/web/js/app/rule-matching.js');
    const { ruleMatchesFeature } = await import('/gbdraw/web/js/services/rule-matchers.js');
    const { runDiagramHelperOperation, DIAGRAM_HELPER_OPERATIONS } = await import('/gbdraw/web/js/services/diagram-generation.js');
    const features = Array.from({ length: 25000 }, (_, i) => ({ type: 'CDS', svg_id: `f${i}`, qualifiers: { product: [i % 2 ? 'other' : 'β-lactamase'] } }));
    // The app's drawing store and card files, with the 25,000-feature catalog.
    const { state: appState } = await import('/gbdraw/web/js/state.js');
    const state = Object.create(appState, { extractedFeatures: { value: { value: features } } });
    let calls = 0;
    const preparation = createRulePreparation({ state, evaluate: async payload => {
      calls++;
      return (await runDiagramHelperOperation(DIAGRAM_HELPER_OPERATIONS.EVALUATE_RULES, payload)).result;
    } });
    const rules = [{ feat: 'CDS', qual: 'product', val: '\\bβ' }];
    let ticks = 0;
    const timer = setInterval(() => ticks++, 16);
    const start = performance.now();
    await preparation.prepare(rules);
    const preparedMs = performance.now() - start;
    clearInterval(timer);
    const cachedStart = performance.now();
    const synchronous = preparation.prepare(rules) === true;
    const count = features.filter(f => ruleMatchesFeature(f, rules[0])).length;
    return { count, calls, synchronous, ticks, preparedMs, reuseMs: performance.now() - cachedStart };
  });
  console.log(JSON.stringify({ pythonRulePerformance: result }));
  expect(result.count).toBe(12500);
  expect(result.calls).toBe(1);
  expect(result.synchronous).toBe(true);
  expect(result.ticks).toBeGreaterThan(0);
  expect((await getDiagramWorkerActivity(page)).constructions).toBe(1);
});
