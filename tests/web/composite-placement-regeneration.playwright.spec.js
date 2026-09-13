const { test, expect } = require('@playwright/test');
const fs = require('node:fs/promises');
const { execFile } = require('node:child_process');
const { promisify } = require('node:util');
const { openApp } = require('./helpers/app-lifecycle.cjs');

const singleSeed = 'gbdraw/web/gallery/sessions/HmmtDNA_basic_circular.gbdraw-session.json';
const allPlacements = ['auto', 'main', 'outward', 'inward'];

const compositeSeed = async (testInfo) => {
  const single = JSON.parse(await fs.readFile(singleSeed, 'utf8'));
  const lambda = JSON.parse(await fs.readFile('gbdraw/web/gallery/sessions/lambda_basic_linear.gbdraw-session.json', 'utf8'));
  const second = { ...lambda.renderRequest.records[0], recordKey: 'record-2',
    source: { kind: 'genbank', resourceId: 'record-2-genbank' } };
  const input = testInfo.outputPath('two-source-input.json');
  const output = testInfo.outputPath('two-source.gbdraw-session.json.gz');
  await fs.writeFile(input, JSON.stringify({ format: single.format, version: single.version,
    results: [], editorState: { featureCatalog: null }, ui: { ...single.ui, generatedMultiRecordCanvas: true },
    config: { ...single.config, form: { ...single.config.form, multi_record_canvas: true } },
    renderRequest: { ...single.renderRequest, grouping: 'grid',
      records: [...single.renderRequest.records, second], layout: { multiRecordSizeMode: 'auto',
        multiRecordMinRadiusRatio: 0.55, multiRecordColumnGapRatio: 0.1, multiRecordRowGapRatio: 0.05,
        multiRecordPositions: null } },
    resources: { ...single.resources, 'record-2-genbank': lambda.resources[lambda.renderRequest.records[0].source.resourceId] } }));
  // The normal writer creates matching SVG, feature metadata and current file bindings.
  const { stdout, stderr } = await promisify(execFile)('python', ['-m', 'gbdraw.cli', 'circular',
    '--session', input, '-o', testInfo.outputPath('two-source'), '--session_output', output],
  { cwd: testInfo.outputDir, env: { ...process.env, PYTHONPATH: process.cwd() },
    timeout: 1_800_000, maxBuffer: 1_000_000 });
  await fs.writeFile(testInfo.outputPath('seed-cli.log'), stdout + stderr);
  return output;
};

for (const composite of [false, true]) {
  test(`Circular ${composite ? 'composite P2-P8' : 'single P1'} placement capability survives generation`, async ({ browser }, testInfo) => {
    test.setTimeout(1_800_000);
    let context;
    let page;
    const external = [];
    const errors = [];
    const load = async (file) => {
      if (context) await context.close();
      context = await browser.newContext({ viewport: { width: 1600, height: 1000 } });
      await context.route('**/*', (route) => {
        if (new URL(route.request().url()).hostname === '127.0.0.1') return route.continue();
        external.push(route.request().url());
        return route.abort();
      });
      page = await context.newPage();
      page.on('pageerror', (error) => errors.push(error.message));
      page.on('dialog', (dialog) => dialog.type() === 'prompt' ? dialog.accept('Composite placement') : dialog.accept());
      await openApp(page);
      await page.locator('input[accept^=".json,"]').setInputFiles(file);
      await expect.poll(() => page.evaluate(() => !window.__GBDRAW_APP__.sessionImportPending
        && window.__GBDRAW_APP__.extractedFeatures.length > 0), { timeout: 180000 }).toBe(true);
    };
    const features = () => page.evaluate(() => {
      const items = window.__GBDRAW_APP__.extractedFeatures;
      return [...new Set(items.map((feature) => feature.record_key))]
        .map((key) => items.find((feature) => feature.record_key === key && feature.type === 'CDS'));
    });
    const generate = async () => {
      expect(await page.evaluate(() => window.__GBDRAW_APP__.runAnalysis())).toEqual({ status: 'ok' });
      await expect.poll(() => page.evaluate(async () => {
        const { state } = await import('./js/state.js');
        return (await import('./js/services/svg-result-ingestion.js'))
          .isCommittedSvgResultMounted(state.results.value[state.selectedResultIndex.value]);
      })).toBe(true);
    };
    const close = () => page.getByRole('button', { name: 'Close feature popup', exact: true }).click();
    const open = async (feature) => {
      await page.evaluate((feature) => window.__GBDRAW_APP__.openFeatureEditorFromList(feature), feature);
      return page.getByLabel('Feature placement', { exact: true });
    };
    const check = async (stage, enabled = allPlacements, targets = null) => {
      const selected = targets || await features();
      expect(selected).toHaveLength(composite ? 2 : 1);
      for (const feature of selected) {
        const control = await open(feature);
        const semantic = await page.evaluate((feature) => {
          const actions = window.__GBDRAW_APP__.featurePlacementActions;
          return { choices: actions.choices([feature]), value: actions.valueFor(feature) };
        }, feature);
        expect(semantic.choices.filter((choice) => choice.enabled).map((choice) => choice.value)).toEqual(enabled);
        expect(await control.evaluate((node) => [...node.options].filter((option) => !option.disabled)
          .map((option) => option.value))).toEqual(enabled);
        await expect(control).toHaveValue(semantic.value);
        for (const invalid of allPlacements.filter((value) => !enabled.includes(value))) {
          const rejected = await page.evaluate(({ feature, invalid }) => {
            try { window.__GBDRAW_APP__.featurePlacementActions.setPlacement([feature], invalid); return false; }
            catch (error) { return error.message.includes('Unavailable'); }
          }, { feature, invalid });
          expect(rejected).toBe(true);
        }
        await close();
      }
      const state = await page.evaluate(async () => {
        const { state } = await import('./js/state.js');
        const config = await import('./js/services/config.js');
        const { getSessionResourceSource } = await import('./js/services/file-content-cache.js');
        const session = config.getCommittedCanonicalSession();
        const binding = getSessionResourceSource(state.files.c_gb);
        return { records: state.circularRecordList.value, request: session.renderRequest,
          resources: Object.fromEntries(Object.entries(session.resources).map(([id, descriptor]) => [id, { size: descriptor.size, name: descriptor.name }])),
          source: { name: state.files.c_gb.name, size: state.files.c_gb.size,
            components: binding?.descriptors?.map((part) => part.resourceId) || [] },
          form: config.buildConfigData().form, geometry: state.trackSlotResolvedGeometry.value,
          catalog: state.featureCatalog.value, placements: state.featurePlacementOverrides };
      });
      await fs.writeFile(testInfo.outputPath(`${stage}.json`), JSON.stringify({ selected, enabled, ...state }, null, 2));
      return { ...state, recordOrder: selected.map((feature) => [feature.record_key, feature.record_id]) };
    };
    try {
      await load(composite ? await compositeSeed(testInfo) : singleSeed);
      const original = await check('P1-P2-before');
      if (composite) expect(original.source.components).toHaveLength(2);
      await generate();
      const generated = await check('P1-P2-after');
      expect(generated.recordOrder).toEqual(original.recordOrder);
      if (composite) {
        await generate();
        await check('P3-second-Generate');
        const target = (await features())[0];
        const placement = await open(target);
        await placement.selectOption('outward');
        await page.locator('.feature-popup input[placeholder="Edit label text"]').fill('COMPOSITE_RETAINED_LABEL');
        await page.getByRole('button', { name: 'Apply Label', exact: true }).click();
        await close();
        await generate();
        const edited = await check('P4-feature-edit');
        expect(Object.values(edited.placements)).toHaveLength(1);
        expect(await page.evaluate(() => window.__GBDRAW_APP__.results[0].content)).toContain('COMPOSITE_RETAINED_LABEL');
        const download = page.waitForEvent('download');
        await page.getByRole('button', { name: 'Save Session', exact: true }).click();
        const saved = testInfo.outputPath('P5-saved.gbdraw-session.json.gz');
        await (await download).saveAs(saved);
        await load(saved);
        expect((await check('P5-loaded')).placements).toEqual(edited.placements);
        await generate();
        expect((await check('P5-generated')).placements).toEqual(edited.placements);
        for (const mode of ['Linear', 'Circular']) {
          await page.getByRole('button', { name: mode, exact: true }).click();
          await expect.poll(() => page.evaluate(() => window.__GBDRAW_APP__.mode)).toBe(mode.toLowerCase());
        }
        expect((await check('P6-returned')).placements).toEqual(edited.placements);
        await generate();
        expect((await check('P6-generated')).placements).toEqual(edited.placements);

        // P8: capability follows draft track semantics, even while the selected value differs.
        for (const separate of [false, true]) {
          await page.getByRole('checkbox', { name: 'Separate Strands', exact: true }).setChecked(separate);
          for (const preset of ['tuckin', 'spreadout', 'middle']) {
            await page.locator('#circular-track-preset').selectOption(preset);
            const valid = preset === 'middle' ? allPlacements : ['auto', 'main'];
            await check(`P8-${separate}-${preset}-draft`, valid);
            // Clear the now-unsupported outward intent through the existing Auto action.
            const control = await open((await features())[0]);
            await control.selectOption('auto');
            await close();
            await generate();
            const rendered = await check(`P8-${separate}-${preset}-generated`, valid);
            expect(rendered.geometry.records.every((record) =>
              JSON.stringify(record.featurePlacementTargets.map((placement) => placement.side || placement.kind))
              === JSON.stringify(valid.filter((value) => value !== 'auto')))).toBe(true);
          }
        }
        // P7: an actual new File is a replacement even with identical name and bytes.
        const old = await features();
        const source = await page.evaluate(async () => {
          const { state } = await import('./js/state.js');
          return { name: state.files.c_gb.name, mimeType: state.files.c_gb.type,
            bytes: [...await (await import('./js/services/file-content-cache.js')).readFileBytes(state.files.c_gb)] };
        });
        await page.getByLabel('GenBank/DDBJ File', { exact: true }).setInputFiles({
          name: source.name, mimeType: source.mimeType, buffer: Buffer.from(source.bytes) });
        await expect.poll(() => page.evaluate(async () => (await import('./js/state.js'))
          .state.circularRecordDiscovery.status)).toBe('ready');
        expect((await check('P7-replaced-before-Generate', [], old)).placements).toEqual({});
        await generate();
        await check('P7-replaced-after-Generate');
      }
      expect(errors).toEqual([]);
      expect(external).toEqual([]);
    } finally {
      if (context) await context.close();
    }
  });
}
