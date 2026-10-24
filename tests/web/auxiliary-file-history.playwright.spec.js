const { test, expect } = require('@playwright/test');
const { openApp, reveal, evaluateWithRetainedPromise } = require('./helpers/app-lifecycle.cjs');
const { readFileSync } = require('node:fs');
const { join } = require('node:path');
const { load, generate, switchMode } = require('./helpers/mode-transition.cjs');

const cases = [
  { label: 'Default colors file', key: 'd_color', text: 'CDS\t#123456\n', value: 'a.currentColors.CDS', expected: '#123456' },
  { label: 'Specific rules file', key: 't_color', text: 'CDS\tproduct\tNADH\t#123456\tImported\n', value: 'a.manualSpecificRules', expected: [{ feat: 'CDS', qual: 'product', val: 'NADH', color: '#123456', cap: 'Imported', fromFile: true }] },
  { label: 'Priority File (TSV)', key: 'qualifier_priority', text: 'CDS\tgene,product\n', value: 'a.manualPriorityRules', expected: [{ feat: 'CDS', order: 'gene,product', importedFrom: 'qualifier_priority' }] },
  { label: 'Whitelist File', key: 'whitelist', text: 'CDS\tgene\tmygene\n', value: 'a.manualWhitelist', expected: [{ feat: 'CDS', qual: 'gene', key: 'mygene', importedFrom: 'whitelist' }] },
  { label: 'Blacklist File', key: 'blacklist', text: 'mygene\n', value: 'a.manualBlacklist', expected: 'mygene' }
];

for (const sample of cases) {
  test(`${sample.label} imports and Undo/Redo restores its content with its selection`, async ({ page }) => {
    await openApp(page);
    await page.evaluate((key) => {
      const a = window.__GBDRAW_APP__;
      a.form.labels_mode = 'out';
      if (key === 'blacklist') a.setLabelFilterMode('Blacklist');
      if (key === 'whitelist') a.setLabelFilterMode('Whitelist');
    }, sample.key);
    const inspect = () => page.evaluate((key) => {
      const a = window.__GBDRAW_APP__;
      const values = { d_color: a.currentColors.CDS, t_color: a.manualSpecificRules,
        qualifier_priority: a.manualPriorityRules, whitelist: a.manualWhitelist,
        blacklist: a.manualBlacklist };
      return { file: a.files[key]?.name || null, value: JSON.parse(JSON.stringify(values[key])) };
    }, sample.key);
    const before = await inspect();
    const expected = sample.key === 'blacklist' ? `${before.value}, mygene` : sample.expected;
    await page.getByLabel(sample.label, { exact: true }).setInputFiles({ name: 'import.tsv', mimeType: 'text/plain', buffer: Buffer.from(sample.text) });
    await page.evaluate(async (key) => { const a = window.__GBDRAW_APP__; await a.waitForAuxiliaryFileImport(a.files[key]); }, sample.key);
    await expect.poll(inspect).toEqual({ file: 'import.tsv', value: expected });
    await expect.poll(() => page.evaluate(() => window.__GBDRAW_HISTORY__.getUndoCount())).toBe(1);
    await page.getByRole('button', { name: 'Undo', exact: true }).click();
    await expect.poll(inspect).toEqual(before);
    await page.getByRole('button', { name: 'Redo', exact: true }).click();
    await expect.poll(inspect).toEqual({ file: 'import.tsv', value: expected });
  });
}

// OV-1280 (D-19): Remove on a card leaves the settings as if its file had
// never been loaded, with every hand edit made before or after the load; one
// Undo brings back the card and what it loaded.
const removeCases = [
  { label: 'Default colors file', key: 'd_color', text: 'CDS\t#123456\ntRNA\t#654321\n',
    loaded: { CDS: '#123456', tRNA: '#abcdef' }, removed: { CDS: '#111111', tRNA: '#abcdef' } },
  { label: 'Specific rules file', key: 't_color', text: 'CDS\tproduct\talpha\t#123456\tAlpha\nCDS\tproduct\tbeta\t#222222\tBeta\n',
    loaded: [['hand', '#445566', null], ['alpha', '#123456', 't_color'], ['beta', '#999999', null]],
    removed: [['hand', '#445566', null], ['beta', '#999999', null]] },
  { label: 'Priority File (TSV)', key: 'qualifier_priority', text: 'CDS\tproduct,gene\nrRNA\tnote\n',
    loaded: [['CDS', 'product,gene', 'qualifier_priority'], ['tRNA', 'product', null], ['rRNA', 'locus_tag', null]],
    removed: [['CDS', 'gene', null], ['tRNA', 'product', null], ['rRNA', 'locus_tag', null]] },
  { label: 'Whitelist File', key: 'whitelist', text: 'CDS\tgene\talpha\nCDS\tgene\tbeta\n',
    loaded: [['mine', null], ['alpha', 'whitelist'], ['betas', null]], removed: [['mine', null], ['betas', null]] },
  { label: 'Blacklist File', key: 'blacklist', text: 'alpha\nbeta\n',
    loaded: 'hypothetical, alpha, beta, mine', removed: 'hypothetical, mine' }
];

for (const sample of removeCases) {
  test(`Remove on the ${sample.label} card takes back only what the file loaded, and Undo restores it (OV-1280)`, async ({ page }) => {
    await openApp(page);
    const inspect = () => page.evaluate((key) => {
      const a = window.__GBDRAW_APP__;
      const values = {
        d_color: () => ({ CDS: a.currentColors.CDS, tRNA: a.currentColors.tRNA }),
        t_color: () => a.manualSpecificRules.map((rule) => [rule.val, rule.color, rule.fromFile ? 't_color' : null]),
        qualifier_priority: () => a.manualPriorityRules.map((rule) => [rule.feat, rule.order, rule.importedFrom ?? null]),
        whitelist: () => a.manualWhitelist.map((rule) => [rule.key, rule.importedFrom ?? null]),
        blacklist: () => a.manualBlacklist
      };
      return { file: a.files[key]?.name || null, value: values[key]() };
    }, sample.key);
    // A hand edit before the load.
    await evaluateWithRetainedPromise(page, async (key) => {
      const a = window.__GBDRAW_APP__;
      a.form.labels_mode = 'out';
      if (key === 'd_color') a.currentColors.CDS = '#111111';
      if (key === 't_color') {
        Object.assign(a.newSpecRule, { feat: 'CDS', qual: 'product', val: 'hand', color: '#445566', cap: '' });
        await a.addSpecificRule();
      }
      if (key === 'qualifier_priority') {
        for (const [feat, order] of [['CDS', 'gene'], ['tRNA', 'product']]) {
          Object.assign(a.newPriorityRule, { feat, order });
          a.addPriorityRule();
        }
      }
      if (key === 'whitelist') {
        a.setLabelFilterMode('Whitelist');
        a.addWhitelistRule();
        a.setWhitelistRuleField(a.manualWhitelist[0], 'key', 'mine');
      }
      if (key === 'blacklist') {
        a.setLabelFilterMode('Blacklist');
        a.setBlacklistText('hypothetical');
      }
    }, sample.key);
    const before = await inspect();
    await page.getByLabel(sample.label, { exact: true }).setInputFiles({ name: 'import.tsv', mimeType: 'text/plain', buffer: Buffer.from(sample.text) });
    await page.evaluate(async (key) => { const a = window.__GBDRAW_APP__; await a.waitForAuxiliaryFileImport(a.files[key]); }, sample.key);
    // The card's section may be closed; its file input names the card.
    await reveal(page.getByLabel(sample.label, { exact: true }));
    // A hand edit after the load, of a value the file set where there is one.
    const field = { whitelist: ['input[placeholder="Pattern"]', 2, 'betas'],
      blacklist: ['textarea[placeholder="hypothetical, unknown"]', 0, 'hypothetical, alpha, beta, mine'] }[sample.key];
    if (field) {
      await page.locator(field[0]).nth(field[1]).fill(field[2]);
      // A text field's edit is its own History step once the field loses focus.
      await page.evaluate(() => document.activeElement?.blur());
    } else {
      await page.evaluate(async (key) => {
        const a = window.__GBDRAW_APP__;
        if (key === 'd_color') a.currentColors.tRNA = '#abcdef';
        if (key === 't_color') await a.setSpecificRuleField(2, 'color', '#999999');
        if (key === 'qualifier_priority') {
          Object.assign(a.newPriorityRule, { feat: 'rRNA', order: 'locus_tag' });
          a.addPriorityRule();
        }
      }, sample.key);
    }
    const loaded = { file: 'import.tsv', value: sample.loaded };
    await expect.poll(inspect).toEqual(loaded);
    const undoCount = await page.evaluate(() => window.__GBDRAW_HISTORY__.getUndoCount());
    await page.getByRole('button', { name: `Remove ${sample.label}`, exact: true }).click();
    await expect.poll(inspect).toEqual({ file: null, value: sample.removed });
    expect(await page.evaluate(() => window.__GBDRAW_HISTORY__.getUndoCount())).toBe(undoCount + 1);
    await page.getByRole('button', { name: 'Undo', exact: true }).click();
    await expect.poll(inspect).toEqual(loaded);
    expect(before.file).toBeNull();
  });
}

// OV-1280 (REVIEW-1 B1): Remove acts in both drawings. The drawing that is not
// shown takes back its rows from its own state: no rule of it is matched
// against the shown drawing's features, and the shown figure keeps its Legend.
for (const drawn of [false, true]) test(`Remove on the Specific rules card while the other mode is shown leaves that drawing as if the file had never been loaded${drawn ? ' (a Generate drew the file rows)' : ''}`, async ({ browser }) => {
  // A Session load, two Generates, and a Pyodide start come before the checks.
  test.setTimeout(300_000);
  const page = await load(browser);
  const circular = () => page.evaluate(async () => {
    const { state } = await import('./js/state.js');
    const drawing = state.drawingOfMode('circular');
    return {
      rules: drawing.manualSpecificRules.map((rule) => [rule.val, Boolean(rule.fromFile)]),
      legend: drawing.legendEntries.value.map((entry) => entry.caption).filter((caption) => ['Hand', 'Imported'].includes(caption))
    };
  });
  // The mounted figure's text of the two rows.
  const shownLegend = () => page.evaluate(async () => {
    const { state } = await import('./js/state.js');
    const root = state.svgContainer.value?.querySelector('svg');
    return [...(root?.querySelectorAll('text') ?? [])].map((node) => node.textContent).filter((text) => ['Hand', 'Imported'].includes(text));
  });
  await evaluateWithRetainedPromise(page, async () => {
    const a = window.__GBDRAW_APP__;
    Object.assign(a.newSpecRule, { feat: 'CDS', qual: 'product', val: 'cytochrome', color: '#445566', cap: 'Hand' });
    await a.addSpecificRule();
  });
  await page.getByLabel('Specific rules file', { exact: true }).setInputFiles({ name: 'rules.tsv', mimeType: 'text/plain', buffer: Buffer.from('CDS\tproduct\tNADH\t#123456\tImported\n') });
  await page.evaluate(async () => { const a = window.__GBDRAW_APP__; await a.waitForAuxiliaryFileImport(a.files.t_color); });
  const loaded = { rules: [['cytochrome', false], ['NADH', true]], legend: ['Hand', 'Imported'] };
  await expect.poll(circular).toEqual(loaded);
  if (drawn) await generate(page);
  await expect.poll(shownLegend).toEqual(['Hand', 'Imported']);
  await switchMode(page, 'linear');
  await page.evaluate(async (content) => {
    const app = window.__GBDRAW_APP__;
    app.setLinearInputType('gb');
    if (!app.linearSeqs.length) app.addLinearSeq();
    app.setLinearSeqPrimaryFile(0, 'gb', new File([content], 'lambda.gbk', { type: 'text/plain', lastModified: 1 }));
    await app.refreshLinearRecordSelectors();
  }, readFileSync(join(__dirname, '../test_inputs/NC_001416.gb'), 'utf8'));
  await generate(page);
  const linear = () => page.evaluate(async () => {
    const { state } = await import('./js/state.js');
    return {
      legend: state.drawingOfMode('linear').legendEntries.value.map((entry) => entry.caption),
      svg: state.svgContainer.value?.querySelector('svg')?.outerHTML ?? ''
    };
  });
  const linearBefore = await linear();
  await reveal(page.getByLabel('Specific rules file', { exact: true }));
  await page.getByRole('button', { name: 'Remove Specific rules file', exact: true }).click();
  // The rows and their Legend rows go at once, in the Remove's step (REVIEW-3 R3-B1).
  const removed = { rules: [['cytochrome', false]], legend: ['Hand'] };
  await expect.poll(circular).toEqual(removed);
  const undoCount = () => page.evaluate(() => window.__GBDRAW_HISTORY__.getUndoCount());
  const removedAt = await undoCount();
  expect(await linear()).toEqual(linearBefore);
  // One step: Undo brings back the card and the Circular rows.
  await page.getByRole('button', { name: 'Undo', exact: true }).click();
  await expect.poll(circular).toEqual(loaded);
  expect(await page.evaluate(() => window.__GBDRAW_APP__.files.t_color?.name ?? null)).toBe('rules.tsv');
  await page.getByRole('button', { name: 'Redo', exact: true }).click();
  await expect.poll(circular).toEqual(removed);
  // Shown again, the Circular figure is drawn without the file's row.
  await switchMode(page, 'circular');
  await expect.poll(shownLegend, { timeout: 120_000 }).toEqual(['Hand']);
  await expect.poll(circular).toEqual(removed);
  // REVIEW-4 N-1: back in Linear, one Undo of the Remove (a project step:
  // a switch is no step, D-15-51) brings the file's rows back while Circular
  // is not shown; shown again, its figure, last shown without them, is drawn
  // with them.
  await switchMode(page, 'linear');
  const stacks = () => page.evaluate(() => [window.__GBDRAW_HISTORY__.getUndoCount(), window.__GBDRAW_HISTORY__.getRedoCount()].join('/'));
  const beforeUndo = await stacks();
  // The click waits while a rerender keeps the Undo button disabled (OV-1339).
  await page.getByRole('button', { name: 'Undo', exact: true }).click();
  await expect.poll(stacks, { timeout: 60_000 }).not.toBe(beforeUndo);
  expect(await undoCount()).toBe(removedAt - 1);
  // The rules and their Legend rows come back now, while Circular is not
  // shown (OV-1655).
  await expect.poll(circular).toEqual(loaded);
  expect(await page.evaluate(() => window.__GBDRAW_APP__.mode)).toBe('linear');
  await switchMode(page, 'circular');
  await expect.poll(shownLegend, { timeout: 120_000 }).toEqual(['Hand', 'Imported']);
});

// OV-1267: Reset Settings removes the settings a table loaded, so the card
// that names the table goes too; one Undo brings back both.
test('Reset Settings removes each table card with the settings it loaded, and Undo restores both', async ({ page }) => {
  await openApp(page);
  await page.evaluate(() => { window.__GBDRAW_APP__.form.labels_mode = 'out'; });
  const inspect = () => page.evaluate(() => {
    const a = window.__GBDRAW_APP__;
    const plain = (value) => JSON.parse(JSON.stringify(value));
    return {
      files: ['d_color', 't_color', 'qualifier_priority', 'whitelist', 'blacklist'].map((key) => a.files[key]?.name || null),
      values: plain([a.currentColors.CDS, a.manualSpecificRules, a.manualPriorityRules, a.manualWhitelist, a.manualBlacklist])
    };
  });
  const defaults = await inspect();
  for (const sample of cases) {
    // The Whitelist and Blacklist cards show under their filter mode.
    await page.evaluate((key) => {
      const a = window.__GBDRAW_APP__;
      if (key === 'whitelist') a.setLabelFilterMode('Whitelist');
      if (key === 'blacklist') a.setLabelFilterMode('Blacklist');
    }, sample.key);
    await page.getByLabel(sample.label, { exact: true }).setInputFiles({ name: `${sample.key}.tsv`, mimeType: 'text/plain', buffer: Buffer.from(sample.text) });
    await page.evaluate(async (key) => { const a = window.__GBDRAW_APP__; await a.waitForAuxiliaryFileImport(a.files[key]); }, sample.key);
  }
  const loaded = await inspect();
  expect(loaded.files).toEqual(['d_color.tsv', 't_color.tsv', 'qualifier_priority.tsv', 'whitelist.tsv', 'blacklist.tsv']);
  expect(loaded.values).not.toEqual(defaults.values);
  await page.evaluate(() => window.__GBDRAW_APP__.resetSettings());
  await expect.poll(inspect).toEqual(defaults);
  await page.getByRole('button', { name: 'Undo', exact: true }).click();
  await expect.poll(inspect).toEqual(loaded);
});

test('a slow earlier specific table cannot overwrite a later uploaded table', async ({ page }) => {
  await page.addInitScript(() => {
    const read = File.prototype.arrayBuffer;
    File.prototype.arrayBuffer = async function () {
      const bytes = await read.call(this);
      if (this.name === 'older.tsv') {
        window.__HELD_FILE__ = this;
        await new Promise((resolve) => { window.__RELEASE_READ__ = resolve; });
      }
      return bytes;
    };
  });
  await openApp(page);
  const upload = page.getByLabel('Specific rules file', { exact: true });
  await upload.setInputFiles({ name: 'older.tsv', mimeType: 'text/plain', buffer: Buffer.from('CDS\tproduct\tNADH\t#ff0000\tOlder\n') });
  await page.waitForFunction(() => window.__RELEASE_READ__);
  await upload.setInputFiles({ name: 'newer.tsv', mimeType: 'text/plain', buffer: Buffer.from('CDS\tproduct\tNADH\t#0000ff\tNewer\n') });
  await page.evaluate(async () => { const a = window.__GBDRAW_APP__; await a.waitForAuxiliaryFileImport(a.files.t_color); });
  await expect.poll(() => page.evaluate(() => window.__GBDRAW_APP__.manualSpecificRules[0]?.cap)).toBe('Newer');
  await page.evaluate(async () => {
    window.__RELEASE_READ__();
    await window.__GBDRAW_APP__.waitForAuxiliaryFileImport(window.__HELD_FILE__);
  });
  expect(await page.evaluate(() => ({ file: window.__GBDRAW_APP__.files.t_color.name, captions: window.__GBDRAW_APP__.manualSpecificRules.map((rule) => rule.cap) }))).toEqual({ file: 'newer.tsv', captions: ['Newer'] });
});

// H-GUI-16: a Specific rules import still applying keeps History's checkpoint
// open, so a Priority File uploaded meanwhile is no History step; Generate
// pressed at once must still wait for both imports (OV-1193), and the Priority
// upload must not discard the Specific rules import (OV-1205).
test('Generate waits for the Specific rules and Priority File imports still applying (OV-1193)', async ({ browser }) => {
  // A Session load and a Pyodide start come before the checks.
  test.setTimeout(180_000);
  const page = await load(browser);
  await page.evaluate(() => {
    window.__RUN_RESOURCES__ = [];
    window.__RELEASE__ = {};
    const post = Worker.prototype.postMessage;
    Worker.prototype.postMessage = function trackedRun(message, ...rest) {
      if (message?.type === 'run') window.__RUN_RESOURCES__.push(message.payload.resourceManifest.map((entry) => entry.resourceId));
      return post.call(this, message, ...rest);
    };
    const read = File.prototype.arrayBuffer;
    File.prototype.arrayBuffer = async function heldRead() {
      const bytes = await read.call(this);
      if (this.name.startsWith('held-')) await new Promise((resolve) => { window.__RELEASE__[this.name] = resolve; });
      return bytes;
    };
  });
  const upload = async (label, name, text) => {
    await page.getByLabel(label, { exact: true }).setInputFiles({ name, mimeType: 'text/plain', buffer: Buffer.from(text) });
    await page.waitForFunction((held) => window.__RELEASE__[held], name);
  };
  await upload('Specific rules file', 'held-rules.tsv', 'CDS\tproduct\tNADH\t#123456\tImported\n');
  await upload('Priority File (TSV)', 'held-priority.tsv', 'CDS\tgene,product\n');
  const key = await page.evaluate(async () => (await import('./js/state.js')).state.resultGenerationKey.value);
  const generateStarts = () => page.evaluate(() => window.__MODE_EVENTS__.filter((event) => event.name === 'generate-start').length);
  const startsBefore = await generateStarts();
  await page.getByRole('button', { name: 'Generate Diagram', exact: true }).click();
  // The click is honoured: Generate shows its processing state while the imports apply.
  await expect.poll(() => page.evaluate(() => window.__GBDRAW_APP__.processing)).toBe(true);
  await page.evaluate(async () => {
    const a = window.__GBDRAW_APP__;
    window.__RELEASE__['held-rules.tsv']();
    await a.waitForAuxiliaryFileImport(a.files.t_color);
  });
  // Generate starts building its request only after the Priority File applies;
  // one that did not wait would start once the rules applied (`generate-start`).
  await page.waitForFunction((before) => window.__MODE_EVENTS__.filter((event) => event.name === 'generate-start').length > before,
    startsBefore, { timeout: 5000 }).catch(() => {});
  expect(await generateStarts()).toBe(startsBefore);
  await page.evaluate(() => window.__RELEASE__['held-priority.tsv']());
  await expect.poll(() => page.evaluate(async () => {
    const { state } = await import('./js/state.js');
    return { key: state.resultGenerationKey.value, processing: state.processing.value, error: state.errorLog.value };
  }), { timeout: 180000 }).toEqual({ key: key + 1, processing: false, error: null });
  // The first render request is the Generate's; a later rerender may follow it.
  expect(await page.evaluate(() => window.__RUN_RESOURCES__[0]))
    .toEqual(expect.arrayContaining(['colors-color-table-file', 'qualifier-priority-file']));
  expect(await page.evaluate(() => {
    const a = window.__GBDRAW_APP__;
    return { priority: JSON.parse(JSON.stringify(a.manualPriorityRules)), rules: a.manualSpecificRules.map((rule) => rule.cap) };
  })).toEqual({ priority: [{ feat: 'CDS', order: 'gene,product', importedFrom: 'qualifier_priority' }], rules: ['Imported'] });
});

// OV-1205: a Priority File uploaded while the Specific rules are
// being evaluated changes the rule inputs under that evaluation; the rule
// owner prepares the rules again, so they apply and the file stays selected.
test('a Priority File uploaded while the Specific rules are evaluated leaves them applied (OV-1205)', async ({ browser }) => {
  // A Session load and a Pyodide start come before the checks.
  test.setTimeout(180_000);
  const page = await load(browser);
  await page.evaluate(() => {
    window.__HOLD_EVALUATION__ = true;
    const post = Worker.prototype.postMessage;
    Worker.prototype.postMessage = function heldEvaluation(message, ...rest) {
      if (window.__HOLD_EVALUATION__ && message?.type === 'helper' && message.operation === 'evaluateRules'
        && JSON.stringify(message).includes('NADH')) {
        window.__HOLD_EVALUATION__ = false;
        window.__RELEASE_EVALUATION__ = () => post.call(this, message, ...rest);
        return undefined;
      }
      return post.call(this, message, ...rest);
    };
  });
  await page.getByLabel('Specific rules file', { exact: true }).setInputFiles({ name: 'rules.tsv', mimeType: 'text/plain', buffer: Buffer.from('CDS\tproduct\tNADH\t#123456\tImported\n') });
  await page.waitForFunction(() => window.__RELEASE_EVALUATION__);
  await page.getByLabel('Priority File (TSV)', { exact: true }).setInputFiles({ name: 'priority.tsv', mimeType: 'text/plain', buffer: Buffer.from('CDS\tgene,product\n') });
  const applied = await page.evaluate(async () => {
    const a = window.__GBDRAW_APP__;
    window.__RELEASE_EVALUATION__();
    return [await a.waitForAuxiliaryFileImport(a.files.t_color), await a.waitForAuxiliaryFileImport(a.files.qualifier_priority)];
  });
  expect(applied).toEqual([true, true]);
  expect(await page.evaluate(() => {
    const a = window.__GBDRAW_APP__;
    return { file: a.files.t_color?.name || null, rules: a.manualSpecificRules.map((rule) => rule.cap), error: a.errorLog };
  })).toEqual({ file: 'rules.tsv', rules: ['Imported'], error: null });
});
