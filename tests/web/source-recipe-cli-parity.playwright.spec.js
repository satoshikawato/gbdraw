// G-E (Web GUI audit 2026-09-30; GE-03, TR-07, CO-07 class): the Run Info
// Source recipe reproduces the Result. One probe per option family changes one
// setting and generates; the recipe then runs with the CLI in a fresh
// directory, and its SVG must equal the GUI Result
// (tests/utils/svg_compare.compare_svgs, interactive binding metadata ignored).
// The CLI SVG must also differ from the Result before the change, so a setting
// that does not reach the drawing, or a comparison that misses it, fails. A
// lossy recipe must be unavailable with a reason (P10) instead of running.
// The full sweep stays an audit tool (remediation/W8_verification_workflow.md).
const { test, expect } = require('@playwright/test');
const { spawnSync } = require('node:child_process');
const { readFileSync, writeFileSync } = require('node:fs');
const path = require('node:path');
const { evaluateWithRetainedPromise, generateAndWaitForResult, reveal } = require('./helpers/app-lifecycle.cjs');
const { BATCH_FIXTURE, openFresh, settle } = require('./helpers/audit-browser.cjs');

test.describe.configure({ retries: 0 });

const [TESTA, TESTB] = readFileSync(BATCH_FIXTURE, 'utf8').split(/^\/\/\s*$/m)
  .map((chunk) => chunk.trim()).filter((chunk) => chunk.startsWith('LOCUS')).map((chunk) => `${chunk}\n//\n`);
// TESTA and TESTB carry the same CDS at the same coordinates.
const BLAST = ['TESTA\tTESTB\t100.0\t300\t0\t0\t301\t600\t301\t600\t1e-80\t550',
  'TESTA\tTESTB\t95.0\t200\t10\t0\t2001\t2200\t2001\t2200\t1e-60\t300'].join('\n') + '\n';
const replayEnv = { ...process.env };
delete replayEnv.PYTHONPATH;
delete replayEnv.PYTHONHOME;

// A long page promise awaited inside page.evaluate can be collected by the inspector (OV-327).
const setApp = (page, body) => evaluateWithRetainedPromise(page, body);
const slotLabel = async (page, slot, text) => {
  const field = page.getByRole('group', { name: `Circular track slot ${slot}`, exact: true })
    .getByRole('textbox', { name: 'Track legend label' });
  await field.fill(text);
  await field.press('Tab');
};

// [name, change, { lossy }]. A lossy probe must yield an unavailable recipe.
const CIRCULAR = [
  ['definition font size', (page) => setApp(page, () => { window.__GBDRAW_APP__.adv.def_font_size = 30; })],
  ['track type', (page) => setApp(page, () => { window.__GBDRAW_APP__.form.track_type = 'middle'; })],
  ['GC window and step', (page) => setApp(page, () => {
    Object.assign(window.__GBDRAW_APP__.adv, { window_size: 200, step_size: 50 });
  })],
  ['legend position', (page) => setApp(page, () => { window.__GBDRAW_APP__.form.legend = 'upper_right'; })],
  ['label blacklist', (page) => setApp(page, () => {
    const app = window.__GBDRAW_APP__;
    app.setLabelFilterMode('Blacklist');
    app.manualBlacklist = 'duplicate';
  })],
  ['feature color rule', (page) => setApp(page, async () => {
    const app = window.__GBDRAW_APP__;
    Object.assign(app.newSpecRule, { feat: 'CDS', qual: 'product', val: 'gtg start', color: '#c83366', cap: 'GTG start' });
    await app.addSpecificRule();
  })],
  ['track legend label with "," and " #"', async (page) => {
    const button = page.locator('button[aria-controls="circular-custom-track-slots-panel"]');
    await reveal(button);
    if (await button.getAttribute('aria-expanded') !== 'true') await button.click();
    await page.getByText('Use custom stack', { exact: true }).locator('input').check();
    await slotLabel(page, 'gc_skew', 'GC skew (1 kb, AT-rich)');
    await slotLabel(page, 'gc_content', 'GC content #1');
  }, { lossy: true }]
];
const LINEAR = [
  ['uploaded BLAST table', async (page, testInfo) => {
    const table = testInfo.outputPath('testa-testb.tsv');
    writeFileSync(table, BLAST);
    await page.getByRole('button', { name: 'Use uploaded BLAST TSV for all adjacent pairs', exact: true }).click();
    await page.locator('input[type="file"][aria-label="BLAST TSV for #1 to #2"]').setInputFiles(table);
    return [table];
  }],
  ['record crop and reverse complement', (page) => setApp(page, () => {
    const app = window.__GBDRAW_APP__;
    app.setLinearRecordCrop(app.linearSeqs[0], 'region_start', 201);
    app.setLinearRecordCrop(app.linearSeqs[0], 'region_end', 3800);
    app.linearSeqs[1].region_reverse = true;
  })],
  // Without a ruler-label font the CLI ruler labels follow --scale_font_size.
  ['ruler scale font without a ruler-label font', (page) => setApp(page, () => {
    const app = window.__GBDRAW_APP__;
    app.form.scale_style = 'ruler';
    app.adv.scale_font_size = 10;
  }), { lossy: true }]
];

const MODES = [
  {
    name: 'Circular',
    probes: CIRCULAR,
    open: async (page, testInfo) => {
      const input = testInfo.outputPath('testa.gb');
      writeFileSync(input, TESTA);
      await openFresh(page);
      await page.getByLabel('GenBank/DDBJ File', { exact: true }).setInputFiles(input);
      await expect.poll(() => page.evaluate(() => window.__GBDRAW_APP__.circularRecordList.length)).toBe(1);
      // Labels on, so the label probes change the drawing.
      await page.evaluate(() => { window.__GBDRAW_APP__.form.labels_mode = 'out'; });
      return [input];
    }
  },
  {
    name: 'Linear',
    probes: LINEAR,
    open: async (page, testInfo) => {
      const inputs = [TESTA, TESTB].map((text, index) => {
        const file = testInfo.outputPath(`record-${index}.gb`);
        writeFileSync(file, text);
        return file;
      });
      await openFresh(page);
      await page.getByRole('button', { name: 'Linear', exact: true }).click();
      await page.waitForFunction(() => window.__GBDRAW_APP__?.mode === 'linear');
      await page.getByTestId('linear-genbank-1').setInputFiles(inputs[0]);
      await page.evaluate(() => window.__GBDRAW_APP__.addLinearSeq());
      await page.getByTestId('linear-genbank-2').setInputFiles(inputs[1]);
      await expect.poll(() => page.evaluate(() => window.__GBDRAW_APP__.linearSeqs.map((seq) => seq.gb?.name)))
        .toEqual(['record-0.gb', 'record-1.gb']);
      return inputs;
    }
  }
];

const generated = async (page, testInfo, name) => {
  await settle(page);
  await generateAndWaitForResult(page);
  await settle(page);
  const { svg, recipe } = await page.evaluate(() => {
    const app = window.__GBDRAW_APP__;
    const source = app.lastRunInfo?.sourceRecipe || {};
    return { svg: app.results[0].content, recipe: { available: source.available, command: source.command || '', reason: source.unavailableReason || '' } };
  });
  const file = testInfo.outputPath(`${name.replace(/\W+/g, '-')}.gui.svg`);
  writeFileSync(file, svg);
  return { file, recipe };
};

const bundle = async (page, testInfo, name) => {
  if (!await page.evaluate(() => window.__GBDRAW_APP__.runInfoHasCliHelperFiles)) return '-';
  const button = page.getByRole('button', { name: /Download reproducibility files/ });
  if (!await button.isVisible()) await page.getByRole('button', { name: /Run info/i }).click();
  const pending = page.waitForEvent('download');
  await button.click();
  const target = testInfo.outputPath(`${name.replace(/\W+/g, '-')}.bundle.zip`);
  await (await pending).saveAs(target);
  return target;
};

for (const mode of MODES) {
  test(`${mode.name}: each Source recipe probe replays in the CLI to the GUI Result`, async ({ page }, testInfo) => {
    test.setTimeout(900_000);
    const inputs = await mode.open(page, testInfo);
    let previous = await generated(page, testInfo, `${mode.name} baseline`);
    expect(previous.recipe.available, previous.recipe.reason).toBe(true);
    for (const [name, change, { lossy = false } = {}] of mode.probes) {
      await test.step(name, async () => {
        inputs.push(...(await change(page, testInfo) || []));
        const current = await generated(page, testInfo, `${mode.name} ${name}`);
        if (lossy) {
          expect.soft(current.recipe.available, `${name}: lossy recipe`).toBe(false);
          expect.soft(current.recipe.reason, `${name}: unavailable reason`).toMatch(/^Source recipe unavailable: /);
          previous = current;
          return;
        }
        expect(current.recipe.available, `${name}: ${current.recipe.reason}`).toBe(true);
        const archive = await bundle(page, testInfo, `${mode.name} ${name}`);
        const replay = spawnSync(process.env.GBDRAW_PYTHON || 'python', [
          path.resolve('tests/web/helpers/source-recipe-replay.py'),
          testInfo.outputPath(`${name.replace(/\W+/g, '-')}-cli`), current.recipe.command,
          current.file, previous.file, archive, ...new Set(inputs)
        ], { encoding: 'utf8', env: replayEnv });
        expect.soft(replay.status, `${mode.name} ${name}: ${current.recipe.command}\n${replay.stdout}${replay.stderr}`).toBe(0);
        previous = current;
      });
    }
  });
}
