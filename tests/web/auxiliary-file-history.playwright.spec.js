const { test, expect } = require('@playwright/test');
const { openApp } = require('./helpers/app-lifecycle.cjs');
const { load } = require('./helpers/mode-transition.cjs');

const cases = [
  { label: 'Default colors file', key: 'd_color', text: 'CDS\t#123456\n', value: 'a.currentColors.CDS', expected: '#123456' },
  { label: 'Specific rules file', key: 't_color', text: 'CDS\tproduct\tNADH\t#123456\tImported\n', value: 'a.manualSpecificRules', expected: [{ feat: 'CDS', qual: 'product', val: 'NADH', color: '#123456', cap: 'Imported', fromFile: true }] },
  { label: 'Priority File (TSV)', key: 'qualifier_priority', text: 'CDS\tgene,product\n', value: 'a.manualPriorityRules', expected: [{ feat: 'CDS', order: 'gene,product' }] },
  { label: 'Whitelist File', key: 'whitelist', text: 'CDS\tgene\tmygene\n', value: 'a.manualWhitelist', expected: [{ feat: 'CDS', qual: 'gene', key: 'mygene' }] },
  { label: 'Blacklist File', key: 'blacklist', text: 'mygene\n', value: 'a.manualBlacklist', expected: 'mygene' }
];

for (const sample of cases) {
  test(`${sample.label} imports and Undo/Redo restores its content with its selection`, async ({ page }) => {
    await openApp(page);
    await page.evaluate((key) => {
      const a = window.__GBDRAW_APP__;
      a.form.labels_mode = 'out';
      if (key === 'blacklist') a.filterMode = 'Blacklist';
      if (key === 'whitelist') a.filterMode = 'Whitelist';
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
  })).toEqual({ priority: [{ feat: 'CDS', order: 'gene,product' }], rules: ['Imported'] });
});
