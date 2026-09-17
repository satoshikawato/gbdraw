const { test, expect } = require('@playwright/test');
const { join } = require('node:path');
const { openApp, generateAndWaitForResult, getDiagramWorkerActivity } = require('./helpers/app-lifecycle.cjs');

const loadSession = async (page) => {
  await openApp(page);
  await page.locator('input[accept^=".json,"]').setInputFiles(join(process.cwd(), 'gbdraw/web/gallery/sessions/HmmtDNA_basic_circular.gbdraw-session.json'));
  await page.waitForFunction(() => !window.__GBDRAW_APP__.sessionImportPending && window.__GBDRAW_APP__.extractedFeatures.length);
};
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
  await page.evaluate(async () => {
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
  expect(messages.join('\n')).toMatch(/(unknown extension|Invalid rule)/);
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
  const importLabel = (pattern, label) => page.evaluate(async ({ pattern, label }) => {
    const input = { files: [new File([`*\tCDS\tproduct\t${pattern}\t${label}\n`], 'labels.tsv')], value: 'labels.tsv' };
    await window.__GBDRAW_APP__.loadLabelOverrideTable({ target: input });
  }, { pattern, label });
  await importLabel('(?i)NADH', 'MATCHED');
  expect(messages.at(-1)).toMatch(/Applied to 7 label/);
  await page.evaluate(() => window.__GBDRAW_APP__.undoHistory());
  expect(await page.locator('.origin-top svg text').filter({ hasText: /^MATCHED$/ }).count()).toBe(0);
  await page.evaluate(() => window.__GBDRAW_APP__.redoHistory());
  expect(await page.locator('.origin-top svg text').filter({ hasText: /^MATCHED$/ }).count()).toBe(7);
  const snapshot = await page.evaluate(() => ({ svg: window.__GBDRAW_APP__.svgContent, overrides: JSON.stringify(window.__GBDRAW_APP__.labelTextFeatureOverrides) }));
  await importLabel('(?<enzyme>NADH)', 'INVALID');
  expect(messages.at(-1)).toMatch(/(unknown extension|Failed to load label)/);
  expect(await page.evaluate(() => ({ svg: window.__GBDRAW_APP__.svgContent, overrides: JSON.stringify(window.__GBDRAW_APP__.labelTextFeatureOverrides) }))).toEqual(snapshot);
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
  await page.evaluate(async () => {
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
  await page.evaluate(() => {
    window.__RELEASE_RULES__ = null;
    window.__HOLD_RULES__ = true;
    window.__OLD_RULE_EDIT__ = window.__GBDRAW_APP__.setSpecificRuleField(0, 'val', '(?i)cytochrome');
  });
  await page.waitForFunction(() => window.__RELEASE_RULES__);
  await page.locator('input[accept^=".json,"]').setInputFiles(join(process.cwd(), 'gbdraw/web/gallery/sessions/HmmtDNA_basic_circular.gbdraw-session.json'));
  await page.waitForFunction(() => !window.__GBDRAW_APP__.sessionImportPending && !window.__GBDRAW_APP__.manualSpecificRules.length);
  await page.evaluate(async () => { window.__RELEASE_RULES__(); await window.__OLD_RULE_EDIT__; });
  expect(await fills(page)).toEqual([]);
  expect(await page.evaluate(() => window.__GBDRAW_APP__.manualSpecificRules)).toEqual([]);
});

test('Session preview stays lazy and saved Python rules remain editable and regenerable', async ({ page, browser }) => {
  test.setTimeout(180000);
  page.on('dialog', dialog => dialog.accept());
  await loadSession(page);
  await page.getByLabel('Specific Table (-t)', { exact: true }).setInputFiles({ name: 'python.tsv', mimeType: 'text/plain', buffer: Buffer.from('CDS\tproduct\t(?P<enzyme>NADH)\t#f01234\t\n') });
  await page.evaluate(async () => { const a = window.__GBDRAW_APP__; await a.waitForAuxiliaryFileImport(a.files.t_color); });
  expect(await fills(page)).toHaveLength(7);
  const download = page.waitForEvent('download');
  await page.evaluate(() => window.__GBDRAW_APP__.saveSessionWithTitle());
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
    await page.evaluate(async (pattern) => {
      const a = window.__GBDRAW_APP__;
      await a.clearAllSpecificRules();
      Object.assign(a.newSpecRule, { feat: 'CDS', qual: 'product', val: pattern, color: '#f01234', cap: '' });
      await a.addSpecificRule();
    }, pattern);
    const ids = await fills(page);
    const tags = await page.evaluate(ids => window.__GBDRAW_APP__.extractedFeatures.filter(f => ids.includes(f.svg_id)).map(f => f.locus_tag).sort(), ids);
    expect(tags).toEqual(expected);
    await page.evaluate(async (pattern) => {
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

test('Python preset rules update the legend and invalid replacements preserve the diagram', async ({ page }) => {
  test.setTimeout(180000);
  const messages = [];
  page.on('dialog', async dialog => { messages.push(dialog.message()); await dialog.accept(); });
  await loadSession(page);
  const preset = await page.evaluate(() => window.__GBDRAW_APP__.specificRulePresets[0]);
  let body = 'CDS\tproduct\t(?i)NADH\t#f01234\tPython preset\n';
  await page.route(`**/${preset.path.replace(/^\.\//, '')}`, route => route.fulfill({ body, contentType: 'text/plain' }));
  await page.evaluate(async id => {
    const a = window.__GBDRAW_APP__;
    a.selectedSpecificPreset = id;
    await a.applySpecificRulePreset();
  }, preset.id);
  expect(await fills(page)).toHaveLength(7);
  expect(await page.locator('.origin-top svg text').filter({ hasText: /^Python preset$/ }).count()).toBe(1);
  const snapshot = await page.evaluate(() => ({ rules: JSON.stringify(window.__GBDRAW_APP__.manualSpecificRules), svg: window.__GBDRAW_APP__.svgContent }));
  body = 'CDS\tproduct\t(?<enzyme>NADH)\t#abcdef\tInvalid\n';
  await page.evaluate(() => window.__GBDRAW_APP__.applySpecificRulePreset());
  expect(messages.at(-1)).toMatch(/Invalid rule/);
  expect(await page.evaluate(() => ({ rules: JSON.stringify(window.__GBDRAW_APP__.manualSpecificRules), svg: window.__GBDRAW_APP__.svgContent }))).toEqual(snapshot);
  await page.evaluate(() => window.__GBDRAW_APP__.undoHistory());
  expect(await fills(page)).toEqual([]);
  await page.evaluate(() => window.__GBDRAW_APP__.redoHistory());
  expect(await fills(page)).toHaveLength(7);
  await generateAndWaitForResult(page);
  expect(await fills(page)).toHaveLength(7);
});

test('invalid color TSV preserves the selected file, preview, rules and History', async ({ page }) => {
  test.setTimeout(180000);
  const messages = [];
  page.on('dialog', async dialog => { messages.push(dialog.message()); await dialog.accept(); });
  await loadSession(page);
  const upload = page.getByLabel('Specific Table (-t)', { exact: true });
  await upload.setInputFiles({ name: 'valid.tsv', mimeType: 'text/plain', buffer: Buffer.from('CDS\tproduct\t(?i)NADH\t#f01234\t\n') });
  await page.evaluate(async () => { const a = window.__GBDRAW_APP__; await a.waitForAuxiliaryFileImport(a.files.t_color); });
  await expect.poll(() => page.evaluate(() => window.__GBDRAW_HISTORY__.getUndoCount())).toBe(1);
  const inspect = () => page.evaluate(() => {
    const a = window.__GBDRAW_APP__;
    return { file: a.files.t_color?.name, rules: JSON.stringify(a.manualSpecificRules), svg: a.svgContent, history: window.__GBDRAW_HISTORY__.getUndoCount() };
  });
  const before = await inspect();
  await upload.setInputFiles({ name: 'invalid.tsv', mimeType: 'text/plain', buffer: Buffer.from('CDS\tproduct\t(?<enzyme>NADH)\t#abcdef\tInvalid\n') });
  await page.evaluate(async () => { const a = window.__GBDRAW_APP__; await a.waitForAuxiliaryFileImport(a.files.t_color); });
  await expect.poll(() => messages.at(-1)).toMatch(/Failed to load rules file/);
  await expect.poll(inspect).toEqual(before);
  await generateAndWaitForResult(page);
  expect(await fills(page)).toHaveLength(7);
});

test('25,000-feature Python preparation keeps the event loop responsive and reuses one worker', async ({ page }) => {
  test.setTimeout(180000);
  await openApp(page);
  expect((await getDiagramWorkerActivity(page)).constructions).toBe(0);
  const result = await page.evaluate(async () => {
    const { createRulePreparation, ruleMatchesFeature } = await import('/gbdraw/web/js/app/rule-matching.js');
    const { runDiagramHelperOperation, DIAGRAM_HELPER_OPERATIONS } = await import('/gbdraw/web/js/services/diagram-generation.js');
    const features = Array.from({ length: 25000 }, (_, i) => ({ type: 'CDS', svg_id: `f${i}`, qualifiers: { product: [i % 2 ? 'other' : 'β-lactamase'] } }));
    const state = { extractedFeatures: { value: features }, manualSpecificRules: [] };
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
