const { test, expect } = require('@playwright/test');
const fs = require('node:fs/promises');
const { execFile } = require('node:child_process');
const { promisify } = require('node:util');
const { openApp, reveal, evaluateWithRetainedPromise } = require('./helpers/app-lifecycle.cjs');
const { showDrawing } = require('./helpers/drawing-bar.cjs');
const { shownSavedDrawing } = require('./helpers/saved-session.cjs');

const singleSeed = 'gbdraw/web/gallery/sessions/HmmtDNA_basic_circular.gbdraw-session.json';
const allPlacements = ['auto', 'main', 'outward', 'inward'];

// A Circular drawing of two GenBank files, one File card each: the CLI writes
// its Session from a request of two records with their own sources.
const twoCardSeed = async (testInfo) => {
  const single = JSON.parse(await fs.readFile(singleSeed, 'utf8'));
  const lambda = JSON.parse(await fs.readFile('gbdraw/web/gallery/sessions/lambda_basic_linear.gbdraw-session.json', 'utf8'));
  const singleDrawing = shownSavedDrawing(single);
  const lambdaRequest = shownSavedDrawing(lambda).renderRequest;
  const second = { ...lambdaRequest.records[0], recordKey: 'record-2',
    source: { kind: 'genbank', resourceId: 'record-2-genbank' } };
  const input = testInfo.outputPath('two-source-input.json');
  const output = testInfo.outputPath('two-source.gbdraw-session.json.gz');
  await fs.writeFile(input, JSON.stringify({ ...single,
    drawings: [{ ...singleDrawing, results: [], editorState: { featureCatalog: null },
      ui: { ...singleDrawing.ui, generatedMultiRecordCanvas: true },
      config: { ...singleDrawing.config, form: { ...singleDrawing.config.form, multi_record_canvas: true } },
      renderRequest: { ...singleDrawing.renderRequest, grouping: 'grid',
        records: [...singleDrawing.renderRequest.records, second], layout: { multiRecordSizeMode: 'auto',
          multiRecordMinRadiusRatio: 0.55, multiRecordColumnGapRatio: 0.1, multiRecordRowGapRatio: 0.05,
          multiRecordPositions: null } } }],
    resources: { ...single.resources, 'record-2-genbank': lambda.resources[lambdaRequest.records[0].source.resourceId] } }));
  // The normal writer creates matching SVG, feature metadata and current file bindings.
  const { stdout, stderr } = await promisify(execFile)('python', ['-m', 'gbdraw.cli', 'circular',
    '--session', input, '-o', testInfo.outputPath('two-source'), '--session_output', output],
  { cwd: testInfo.outputDir, env: { ...process.env, PYTHONPATH: process.cwd() },
    timeout: 1_800_000, maxBuffer: 1_000_000 });
  await fs.writeFile(testInfo.outputPath('seed-cli.log'), stdout + stderr);
  return output;
};

for (const twoCards of [false, true]) {
  test(`Circular ${twoCards ? 'two File cards P2-P8' : 'single P1'} placement capability survives generation`, async ({ browser }, testInfo) => {
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
      page.on('dialog', (dialog) => dialog.accept());
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
      expect(await evaluateWithRetainedPromise(page, () => window.__GBDRAW_APP__.runAnalysis())).toEqual({ status: 'ok' });
      await expect.poll(() => page.evaluate(async () => {
        const { state } = await import('./js/state.js');
        return (await import('./js/services/svg-result-ingestion.js'))
          .isCommittedSvgResultMounted(window.__GBDRAW_APP__.currentResults()[state.selectedResultIndex.value]);
      })).toBe(true);
    };
    const close = () => page.getByRole('button', { name: 'Close feature popup', exact: true }).click();
    const open = async (feature) => {
      await evaluateWithRetainedPromise(page, (feature) => window.__GBDRAW_APP__.openFeatureEditorFromList(feature), feature);
      return page.getByLabel('Feature placement', { exact: true });
    };
    const check = async (stage, enabled = allPlacements, targets = null) => {
      const selected = targets || await features();
      if (!targets) expect(selected).toHaveLength(twoCards ? 2 : 1);
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
            catch (error) { return /no feature track for this placement|not in the displayed Result/.test(error.message); }
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
        const drawing = state.activeDrawing();
        return { records: drawing.circularRecordList.value, request: session.renderRequest,
          resources: Object.fromEntries(Object.entries(session.resources).map(([id, descriptor]) => [id, { size: descriptor.size, name: descriptor.name }])),
          cards: drawing.circularSources.map((card) => {
            const file = state.inputFiles.fileOf(card.gb);
            return { uid: card.uid, name: file.name, size: file.size,
              resourceId: getSessionResourceSource(file)?.resourceId ?? null };
          }),
          form: config.buildConfigData(state.activeDrawing()).form, geometry: state.trackSlotResolvedGeometry.value,
          catalog: state.featureCatalog.value, placements: state.activeDrawing().featurePlacementOverrides };
      });
      await fs.writeFile(testInfo.outputPath(`${stage}.json`), JSON.stringify({ selected, enabled, ...state }, null, 2));
      return { ...state, recordOrder: selected.map((feature) => [feature.record_key, feature.record_id]) };
    };
    try {
      await load(twoCards ? await twoCardSeed(testInfo) : singleSeed);
      const original = await check('P1-P2-before');
      expect(original.cards).toHaveLength(twoCards ? 2 : 1);
      // Each card draws its own record, keyed by the card (a one-record file's identity is its uid).
      expect(original.recordOrder.map(([key]) => key)).toEqual(original.cards.map((card) => card.uid));
      await generate();
      const generated = await check('P1-P2-after');
      expect(generated.recordOrder).toEqual(original.recordOrder);
      if (twoCards) {
        await generate();
        await check('P3-second-Generate');
        const target = (await features())[0];
        const placement = await open(target);
        await placement.selectOption('outward');
        await page.locator('.feature-popup input[placeholder="Edit label text"]').fill('CARD_RETAINED_LABEL');
        await page.getByRole('button', { name: 'Apply Label', exact: true }).click();
        await close();
        await generate();
        const edited = await check('P4-feature-edit');
        expect(Object.values(edited.placements)).toHaveLength(1);
        // The edit stays on the first card's record.
        expect(Object.keys(edited.placements).map((key) => JSON.parse(key)[0])).toEqual([target.record_key]);
        expect(target.record_key).toBe(original.cards[0].uid);
        expect(await page.evaluate(() => window.__GBDRAW_APP__.currentResults()[0].content)).toContain('CARD_RETAINED_LABEL');
        const download = page.waitForEvent('download');
        await page.getByRole('button', { name: 'Save Session', exact: true }).click();
        const saved = testInfo.outputPath('P5-saved.gbdraw-session.json.gz');
        await (await download).saveAs(saved);
        await load(saved);
        expect((await check('P5-loaded')).placements).toEqual(edited.placements);
        await generate();
        expect((await check('P5-generated')).placements).toEqual(edited.placements);
        for (const kind of ['linear', 'circular']) await showDrawing(page, { kind });
        expect((await check('P6-returned')).placements).toEqual(edited.placements);
        await generate();
        expect((await check('P6-generated')).placements).toEqual(edited.placements);

        // P8: capability follows draft track semantics, even while the selected value differs.
        await reveal(page.locator('#circular-track-preset'));
        for (const separate of [false, true]) {
          await page.getByRole('checkbox', { name: 'Separate Strands', exact: true }).setChecked(separate);
          for (const preset of ['tuckin', 'spreadout', 'middle']) {
            await page.locator('#circular-track-preset').selectOption(preset);
            if (!separate && preset === 'tuckin') {
              // Q3: leaving Middle asks first, and Reset clears the outward intent (OV-09).
              await page.getByRole('dialog', { name: 'Reset Feature placements?', exact: true })
                .getByRole('button', { name: 'Reset 1 placement to Auto', exact: true }).click();
            }
            const valid = preset === 'middle' ? allPlacements : ['auto', 'main'];
            await check(`P8-${separate}-${preset}-draft`, valid);
            // The existing Auto action keeps a cleared intent at Auto.
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
        // P7: an actual new File in the first card is a replacement even with
        // identical name and bytes: its record's features lose their placement
        // until Generate, the other card's keep theirs.
        const [replacedFeature, keptFeature] = await features();
        const source = await page.evaluate(async () => {
          const { state } = await import('./js/state.js');
          const file = state.inputFiles.fileOf(state.activeDrawing().circularSources[0].gb);
          return { name: file.name, mimeType: file.type,
            bytes: [...await (await import('./js/services/file-content-cache.js')).readFileBytes(file)] };
        });
        const firstCard = page.locator('[data-circular-source-card]').first();
        await firstCard.locator('summary').click({ position: { x: 8, y: 8 } });
        await firstCard.locator('input[type=file]').setInputFiles({
          name: source.name, mimeType: source.mimeType, buffer: Buffer.from(source.bytes) });
        await expect.poll(() => page.evaluate(async () => (await import('./js/state.js'))
          .state.activeDrawing().circularRecordDiscovery.status)).toBe('ready');
        const replaced = await check('P7-replaced-before-Generate', [], [replacedFeature]);
        expect(replaced.cards.map((card) => card.uid)).toEqual(original.cards.map((card) => card.uid));
        expect(replaced.placements).toEqual({});
        await check('P7-kept-before-Generate', allPlacements, [keptFeature]);
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
