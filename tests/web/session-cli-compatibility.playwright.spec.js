const { test, expect } = require('@playwright/test');
const fs = require('node:fs/promises');
const path = require('node:path');
const { createHash } = require('node:crypto');
const { gunzipSync, gzipSync } = require('node:zlib');
const { execFile } = require('node:child_process');
const { promisify } = require('node:util');
const { CURRENT_REQUEST_SCHEMA, CURRENT_SESSION_VERSION, evaluateWithRetainedPromise, generateAndWaitForResult, openApp, expectSessionLoadedWithoutWarnings, settleLive } = require('./helpers/app-lifecycle.cjs');
const { shownSavedDrawing } = require('./helpers/saved-session.cjs');
const { drawingIds, shownDrawing } = require('./helpers/drawing-bar.cjs');
const { chooseCoordinates, openRegionDialog } = require('./helpers/region-drawing.cjs');

const root = process.cwd();
const mito = path.join(root, 'tests/fixtures/sessions/cli-web-mito.gb');
const lambda = path.join(root, 'tests/test_inputs/NC_001416.gb');
const gff = path.join(root, 'tests/test_inputs/NC_013668.gff3');
const fasta = path.join(root, 'tests/test_inputs/NC_013668.fasta');
const hash = bytes => createHash('sha256').update(bytes).digest('hex');
const readSession = async file => JSON.parse(gunzipSync(await fs.readFile(file)));
const cases = [
  { name: 'original single', mode: 'circular', args: ['--gbk', mito, '--labels', 'out'], sources: [mito] },
  { name: 'composite', mode: 'circular', args: ['--gbk', mito, lambda, '--multi_record_canvas'], sources: [mito, lambda] },
  { name: 'linear', mode: 'linear', args: ['--gbk', mito, lambda], sources: [mito, lambda] },
  {
    name: 'linear crop and reverse', mode: 'linear', sources: [mito, lambda],
    args: ['--gbk', mito, lambda, '--region', '#1:1000-8000', '--reverse_complement', '0', '--reverse_complement', '1']
  },
  { name: 'gff fasta', mode: 'circular', args: ['--gff', gff, '--fasta', fasta], sources: [gff, fasta] },
  // D-01, OV-221: the tables of a CLI Session reach the next Web Generate.
  ...[
    ['linear_tables', [lambda]], ['circular_records_tables', [path.join(root, 'tests/test_inputs/HmmtDNA.gbk')]],
    ['linear_records_blacklist', [path.join(root, 'tests/test_inputs/HmmtDNA.gbk'), lambda]],
    ['linear_blast', [path.join(root, 'tests/fixtures/web_comparison_shared_block.gb')]]
  ]
    .map(([id, sources]) => {
      const entry = require('../fixtures/cli_session_cross_surface/cases.json').cases.find(item => item.id === id);
      const args = entry.args.map(arg => (arg.startsWith('tests/') ? path.join(root, arg) : arg));
      return { name: id.replace(/_/g, ' '), mode: entry.mode, args, sources };
    })
];

const cli = async (mode, args, testInfo, label) => {
  const prefix = testInfo.outputPath(label);
  const file = `${prefix}.gbdraw-session.json.gz`;
  const command = ['-m', 'gbdraw.cli', mode, ...args, '-o', prefix, '--session_output', file];
  const { stdout, stderr } = await promisify(execFile)('python', command, {
    cwd: testInfo.outputDir, env: { ...process.env, PYTHONPATH: root },
    timeout: 1_800_000, maxBuffer: 1_000_000
  });
  await fs.writeFile(testInfo.outputPath(`${label}-cli.json`), JSON.stringify({ command, stdout, stderr }, null, 2));
  return file;
};

const svgSemantics = content => {
  const svg = new DOMParser().parseFromString(content, 'image/svg+xml').documentElement;
  const records = new Map([...svg.querySelectorAll('[data-gbdraw-record-index]')]
    .map(element => [Number(element.getAttribute('data-gbdraw-record-index')), element.getAttribute('data-gbdraw-record-id')]));
  // The CLI's Python and Pyodide differ in the last digits of path numbers;
  // compare them at the 4 decimals of tests/utils/svg_compare.py. The Web's
  // Result drops the trailing space of an open outline path (`L 2,-8 `), which
  // the path grammar ignores.
  const path = d => d && d.trim().replace(/-?\d*\.?\d+(?:e[-+]?\d+)?/gi, number => String(Math.round(Number(number) * 1e4) / 1e4));
  return {
    records: [...records.entries()].sort((a, b) => a[0] - b[0]),
    features: [...svg.querySelectorAll('[data-gbdraw-feature-id]')].map(element =>
      ['data-gbdraw-feature-id', 'data-gbdraw-record-id', 'd', 'fill', 'stroke', 'display', 'visibility']
        .map(key => key === 'd' ? path(element.getAttribute(key)) : element.getAttribute(key))),
    text: [...svg.querySelectorAll('text')].map(element => element.textContent)
  };
};

const snapshot = page => page.evaluate(async () => {
  const { state } = await import('./js/state.js');
  const { getCommittedCanonicalSession } = await import('./js/services/config.js');
  const { getSessionResourceSource, readFileBytes } = await import('./js/services/file-content-cache.js');
  const digest = async bytes => [...new Uint8Array(await crypto.subtle.digest('SHA-256', bytes))]
    .map(byte => byte.toString(16).padStart(2, '0')).join('');
  const drawing = state.activeDrawing();
  // The drawing names its files by ref; the store holds the files.
  const refs = drawing.mode === 'linear' ? drawing.linearSeqs.map(seq => seq.gb)
    : drawing.circularSources.flatMap(card => (drawing.inputType.value === 'gff' ? [card.gff, card.fasta] : [card.gb]));
  const files = refs.map(ref => state.inputFiles.fileOf(ref));
  const session = getCommittedCanonicalSession();
  return {
    mode: state.mode.value, request: session.renderRequest,
    resourceIds: Object.keys(session.resources),
    files: await Promise.all(files.map(async file => {
      const source = getSessionResourceSource(file);
      const parts = source.descriptors || [{
        ...source, name: file.name, type: file.type, lastModified: file.lastModified
      }];
      return {
        name: file.name, type: file.type, lastModified: file.lastModified, size: file.size,
        isArray: Array.isArray(file), sha256: await digest(await readFileBytes(file)),
        parts: await Promise.all(parts.map(async part => ({
          resourceId: part.resourceId, name: part.name, type: part.type, lastModified: part.lastModified,
          sha256: await digest(Uint8Array.from(atob(part.descriptor.data), char => char.charCodeAt(0)))
        })))
      };
    })),
    selectedIndex: state.selectedResultIndex.value,
    selected: window.__GBDRAW_APP__.currentResults()[state.selectedResultIndex.value]?.content || '',
    mounted: state.svgContainer.value?.querySelector('svg')?.outerHTML || '',
    error: state.errorLog.value
  };
});

const load = async (browser, file, viewport, requirePreview = true) => {
  const context = await browser.newContext({ viewport });
  const external = [], errors = [], dialogs = [];
  await context.route('**/*', route => {
    if (new URL(route.request().url()).hostname === '127.0.0.1') return route.continue();
    external.push(route.request().url());
    return route.abort();
  });
  const page = await context.newPage();
  page.on('pageerror', error => errors.push(String(error)));
  page.on('dialog', async dialog => { dialogs.push(dialog.message()); await dialog.accept(); });
  await openApp(page);
  await page.locator('input[accept^=".json,"]').setInputFiles(file);
  await expectSessionLoadedWithoutWarnings(page, { timeout: 180_000 });
  await expect.poll(() => page.evaluate(() => !window.__GBDRAW_APP__.sessionImportPending)).toBe(true);
  // On a phone a Load with a Result shows the Preview pane, as Generate does.
  if (requirePreview) await expect(page.locator('.origin-top svg')).toBeVisible();
  return { context, page, external, errors, dialogs };
};

for (const entry of cases) {
  test(`Session current CLI ${entry.name} survives Web Generate and bidirectional replay`, async ({ browser }, testInfo) => {
    test.setTimeout(1_800_000);
    expect(hash(await fs.readFile(mito))).toBe('f2e922c26561d37a8d3b410ba972526c709cad5635f598bbdcfc234f02642c92');
    const sourceHashes = await Promise.all(entry.sources.map(async file => hash(await fs.readFile(file))));
    let file = await cli(entry.mode, entry.args, testInfo, 'fresh-cli');
    let webFile;
    for (const phase of ['fresh-cli', 'cli-replay', 'web-cli-replay']) {
      if (phase === 'cli-replay') file = await cli(entry.mode, ['--session', file], testInfo, phase);
      if (phase === 'web-cli-replay') file = await cli(entry.mode, ['--session', webFile], testInfo, phase);
      const session = await readSession(file);
      expect(session.version).toBe(CURRENT_SESSION_VERSION);
      expect(shownSavedDrawing(session).renderRequest.schema).toBe(CURRENT_REQUEST_SCHEMA);
      expect(shownSavedDrawing(session).webFiles.bindings.schema).toBe(2);
      if (entry.name === 'composite') {
        // The two files of `--gbk a b` are two Circular File cards.
        expect(shownSavedDrawing(session).webFiles.bindings.circularSources.map(card => Boolean(card.gb))).toEqual([true, true]);
      }
      const run = await load(browser, file,
        entry.name === 'original single' ? { width: 390, height: 844 } : { width: 1600, height: 1000 });
      try {
        const before = await snapshot(run.page);
        expect(before.mode).toBe(entry.mode);
        expect(before.request).toEqual(shownSavedDrawing(session).renderRequest);
        // A CLI Session binds each input file once; the Web shows one row per
        // drawn record, so a two-record file (linear blast) fills two rows bound
        // to the same resource, and a Session the Web saved lists both rows. The
        // checks below read each resource once on either side.
        const resourceKey = parts => parts.map(part => part.resourceId).join('+');
        const distinct = (items, parts) => [...new Map(items.map(item => [resourceKey(parts(item)), item])).values()];
        const distinctFiles = files => distinct(files, file => file.parts);
        const distinctBindings = items => distinct(items, binding => binding.components || [binding]);
        expect(distinctFiles(before.files).flatMap(file => file.parts.map(part => part.sha256))).toEqual(sourceHashes);
        expect(before.files.every(file => !file.isArray)).toBe(true);
        // A records table binds its input through the request record, not the
        // Web inventory; the file bytes above identify it.
        if (!entry.args.includes('--records_table')) {
          const inventory = shownSavedDrawing(session).webFiles.bindings;
          const bindings = entry.mode === 'linear' ? inventory.linearSeqs.map(seq => seq.gb)
            : inventory.circularSources.flatMap(card => (entry.name === 'gff fasta' ? [card.gff, card.fasta] : [card.gb]));
          const metadata = file => ({ name: file.name, type: file.type, lastModified: file.lastModified });
          expect(distinctFiles(before.files).map(metadata)).toEqual(distinctBindings(bindings).map(metadata));
          expect(distinctFiles(before.files).flatMap(file => file.parts.map(part => ({ resourceId: part.resourceId, ...metadata(part) }))))
            .toEqual(distinctBindings(bindings).flatMap(binding => (binding.components || [binding]).map(part => ({ resourceId: part.resourceId, ...metadata(part) }))));
          // Every row that shares a resource shows the same file.
          for (const file of before.files) {
            expect(metadata(file)).toEqual(metadata(before.files.find(other => resourceKey(other.parts) === resourceKey(file.parts))));
          }
        }
        for (const part of before.files.flatMap(file => file.parts)) expect(before.resourceIds).toContain(part.resourceId);
        const expectedSvg = await run.page.evaluate(svgSemantics, before.selected);
        expect(expectedSvg.records.length).toBe(shownSavedDrawing(session).renderRequest.records.length);
        expect(await run.page.evaluate(svgSemantics, before.mounted)).toEqual(expectedSvg);
        // A read-only imported comparison (-b) is reused through Inherit (D-04):
        // the user picks it in the comparison resolution panel before Generate.
        if (entry.args.includes('-b')) {
          const resolution = run.page.locator('[data-imported-comparison-resolution]');
          await resolution.getByRole('button', { name: 'Inherit saved comparison' }).click();
          await expect(resolution).toContainText('INHERIT');
        }
        await generateAndWaitForResult(run.page);
        const after = await snapshot(run.page);
        const generatedSvg = await run.page.evaluate(svgSemantics, after.selected);
        // Generate draws the figure of the loaded Session (D-01).
        expect(generatedSvg).toEqual(expectedSvg);
        expect(after.request.records).toHaveLength(shownSavedDrawing(session).renderRequest.records.length);
        // Generate keeps each record's crop and orientation.
        const transform = record => [record.region?.start ?? null, record.region?.end ?? null,
          Boolean(record.region?.reverseComplement || record.presentation.reverseComplement)];
        expect(after.request.records.map(transform)).toEqual(shownSavedDrawing(session).renderRequest.records.map(transform));
        expect(after.files).toEqual(before.files);
        expect(distinctFiles(after.files).flatMap(file => file.parts.map(part => part.sha256))).toEqual(sourceHashes);
        expect(await run.page.evaluate(svgSemantics, after.mounted)).toEqual(generatedSvg);
        const download = run.page.waitForEvent('download');
        await run.page.getByRole('button', { name: 'SVG', exact: true }).click();
        const exported = testInfo.outputPath(`${phase}-export.svg`);
        await (await download).saveAs(exported);
        expect(await run.page.evaluate(svgSemantics, await fs.readFile(exported, 'utf8'))).toEqual(generatedSvg);
        if (phase === 'cli-replay') {
          const saved = run.page.waitForEvent('download');
          const outcome = await evaluateWithRetainedPromise(run.page, async () => await window.__GBDRAW_APP__.saveSessionWithTitle());
          expect(outcome.status).toBe('saved');
          webFile = testInfo.outputPath('web-current.gbdraw-session.json.gz');
          await (await saved).saveAs(webFile);
          const webSession = await readSession(webFile);
          expect(shownSavedDrawing(webSession).config.form).toBeTruthy();
        }
        await fs.writeFile(testInfo.outputPath(`${phase}-state.json`), JSON.stringify({
          before: { ...before, selected: undefined, mounted: undefined },
          after: { ...after, selected: undefined, mounted: undefined },
          records: generatedSvg.records, external: run.external, errors: run.errors, dialogs: run.dialogs
        }, null, 2));
        expect(run.external).toEqual([]);
        expect(run.errors).toEqual([]);
      } finally {
        await run.context.close();
      }
    }
  });
}

// OV-2995: a CLI Session that draws one record in two windows (`--region` twice) with a LOSAT comparison names
// drawn units, one per window, in its comparisons. The Web offers Inherit for it (its record indexes cannot hold
// window pairs), and Generate draws the saved figure with the saved unit pairs.
test('Session CLI windowed record with LOSAT regenerates through Inherit', async ({ browser }, testInfo) => {
  test.setTimeout(900_000);
  const bgc = id => path.join(root, `tests/test_inputs/${id}.gbk`);
  const file = await cli('linear', ['--gbk', bgc('BGC0000708'), bgc('BGC0000709'), '--region', 'BGC0000708:1-15000',
    '--region', 'BGC0000708:25000-40000', '--region', 'BGC0000709:1-30000', '--losat', 'losatp', '-f', 'svg'], testInfo, 'windowed');
  const saved = shownSavedDrawing(await readSession(file)).renderRequest;
  const unitPairs = request => request.comparisons.filter(comparison => comparison.kind === 'precomputedProteinComparison')
    .map(comparison => [comparison.queryRecordIndex, comparison.subjectRecordIndex]);
  expect(saved.records.map(record => record.windows.length)).toEqual([2, 1]);
  expect(unitPairs(saved)).toEqual([[0, 2], [1, 2]]);
  const run = await load(browser, file, { width: 1600, height: 1000 });
  try {
    const before = await snapshot(run.page);
    const expectedSvg = await run.page.evaluate(svgSemantics, before.selected);
    const resolution = run.page.locator('[data-imported-comparison-resolution]');
    await resolution.getByRole('button', { name: 'Inherit saved comparison' }).click();
    await expect(resolution).toContainText('INHERIT');
    await generateAndWaitForResult(run.page);
    const after = await snapshot(run.page);
    expect(unitPairs(after.request)).toEqual([[0, 2], [1, 2]]);
    expect(after.request.records.map(record => record.windows)).toEqual(saved.records.map(record => record.windows));
    expect(await run.page.evaluate(svgSemantics, after.selected)).toEqual(expectedSvg);
    expect(run.external).toEqual([]);
    expect(run.errors).toEqual([]);
  } finally {
    await run.context.close();
  }
});

// OV-3003: a Session whose LOSATP pipeline names the window pairs and holds no table and no LOSAT cache (a request
// with LOSATP pairs, as Python writes it). The browser runs LOSAT for exactly those pairs through Inherit and draws
// what the CLI drew from its own search of the same pairs.
test('Session with a windowed LOSATP pipeline and no results runs its pairs through Inherit', async ({ browser }, testInfo) => {
  test.setTimeout(900_000);
  const bgc = id => path.join(root, `tests/test_inputs/${id}.gbk`);
  const file = await cli('linear', ['--gbk', bgc('BGC0000708'), bgc('BGC0000709'), '--region', 'BGC0000708:1-15000',
    '--region', 'BGC0000708:25000-40000', '--region', 'BGC0000709:1-30000', '--losat', 'losatp', '--losatp_mode', 'pairwise',
    '-f', 'svg'], testInfo, 'windowed-pairwise');
  const session = await readSession(file);
  const drawing = shownSavedDrawing(session);
  const [pipeline] = drawing.renderRequest.comparisons;
  expect(drawing.renderRequest.comparisons).toEqual([expect.objectContaining({ kind: 'generatedProteinComparison', mode: 'pairwise', pairs: [] })]);
  // The CLI searched each window of BGC0000708 against BGC0000709: units 0 and 1 against unit 2.
  pipeline.pairs = [{ queryRecordIndex: 0, subjectRecordIndex: 2 }, { queryRecordIndex: 1, subjectRecordIndex: 2 }];
  drawing.losatCache = { entries: [] };
  drawing.losatDerivedCache = { entries: [] };
  const unrendered = testInfo.outputPath('windowed-pipeline.gbdraw-session.json.gz');
  await fs.writeFile(unrendered, gzipSync(JSON.stringify(session)));
  const cliRun = await load(browser, file, { width: 1600, height: 1000 });
  const expectedSvg = await cliRun.page.evaluate(svgSemantics, (await snapshot(cliRun.page)).selected);
  await cliRun.context.close();
  const run = await load(browser, unrendered, { width: 1600, height: 1000 });
  try {
    const resolution = run.page.locator('[data-imported-comparison-resolution]');
    await resolution.getByRole('button', { name: 'Inherit saved comparison' }).click();
    await expect(resolution).toContainText('INHERIT');
    // This spec's pages route every request, so they are not cross-origin isolated: LOSAT runs serially.
    await run.page.evaluate(async () => {
      (await import('./js/state.js')).state.losatExecution.executionMode = 'serial';
    });
    await generateAndWaitForResult(run.page);
    const after = await snapshot(run.page);
    expect(after.request.comparisons.filter(comparison => comparison.kind === 'precomputedProteinComparison')
      .map(comparison => [comparison.queryRecordIndex, comparison.subjectRecordIndex])).toEqual([[0, 2], [1, 2]]);
    expect(await run.page.evaluate(svgSemantics, after.selected)).toEqual(expectedSvg);
    expect(run.external).toEqual([]);
    expect(run.errors).toEqual([]);
  } finally {
    await run.context.close();
  }
});

// Where a figure places its legend and title: the size, each composition
// role's transform, and which legend layout it shows (at 2 decimals).
const placement = content => {
  const svg = new DOMParser().parseFromString(content, 'image/svg+xml').documentElement;
  const round = value => value && value.replace(/-?\d*\.?\d+(?:e[-+]?\d+)?/gi, number => String(Math.round(Number(number) * 100) / 100));
  return {
    size: ['width', 'height', 'viewBox'].map(name => round(svg.getAttribute(name))),
    roles: [...svg.querySelectorAll('[data-gbdraw-composition-role]')]
      .map(element => [element.getAttribute('data-gbdraw-composition-role'), round(element.getAttribute('transform'))]),
    legend: [...svg.querySelectorAll('#legend_horizontal, #legend_vertical')].map(element => [element.id, element.getAttribute('display')])
  };
};

// A Linear drawing made from a region of a Circular record ("New drawing from
// a region…"): `gbdraw render --drawing ID` draws the Web's Save of it as the
// Web's Generate drew it, with the legend and the title in the same place.
// tests/fixtures/sessions/web-region-circular-linear.v51 is this test's Save.
test('Session Linear region drawing of a Circular record draws the Web figure on the CLI', async ({ browser }, testInfo) => {
  test.setTimeout(900_000);
  const run = await load(browser, path.join(root, 'gbdraw/web/gallery/sessions/HmmtDNA_basic_circular.gbdraw-session.json'),
    { width: 1600, height: 1000 });
  try {
    await generateAndWaitForResult(run.page);
    const [source] = await drawingIds(run.page);
    const dialog = await openRegionDialog(run.page);
    await chooseCoordinates(dialog, 1001, 2000);
    await expect.poll(() => run.page.evaluate(() => window.__GBDRAW_APP__.regionDrawings.dialog.status), { timeout: 60_000 })
      .toBe('derived');
    await dialog.locator('[data-region-drawing-create]').click();
    await expect(dialog).toHaveCount(0);
    const [drawing] = (await drawingIds(run.page)).filter(id => id !== source);
    expect(await shownDrawing(run.page)).toEqual({ id: drawing, kind: 'linear' });
    await generateAndWaitForResult(run.page);
    const web = await run.page.evaluate(() => window.__GBDRAW_APP__.currentResults().map(result => String(result.content || '')));
    expect(web).toHaveLength(1);

    const saved = run.page.waitForEvent('download');
    const outcome = await evaluateWithRetainedPromise(run.page, async () => {
      window.__GBDRAW_APP__.sessionTitle = 'drawings-circular-region';
      return await window.__GBDRAW_APP__.saveSessionWithTitle();
    });
    expect(outcome.status).toBe('saved');
    const webFile = testInfo.outputPath('web-region.gbdraw-session.json.gz');
    await (await saved).saveAs(webFile);
    const session = await readSession(webFile);
    // The Web writes the legend and title placement of the new Linear drawing into its request.
    const request = session.drawings.find(item => item.id === drawing).renderRequest;
    expect(request.mode).toBe('linear');
    expect(request.diagramOptions.output).toEqual({ legend: 'bottom', plotTitlePosition: 'bottom' });

    const prefix = testInfo.outputPath('cli-region');
    await promisify(execFile)('python', ['-m', 'gbdraw.cli', 'render', '--session', webFile, '--drawing', drawing,
      '-o', prefix, '-f', 'svg'], {
      cwd: testInfo.outputDir, env: { ...process.env, PYTHONPATH: root }, timeout: 600_000, maxBuffer: 1_000_000
    });
    const cliSvg = await fs.readFile(`${prefix}.svg`, 'utf8');
    expect(await run.page.evaluate(placement, cliSvg)).toEqual(await run.page.evaluate(placement, web[0]));
    expect(await run.page.evaluate(svgSemantics, cliSvg)).toEqual(await run.page.evaluate(svgSemantics, web[0]));
    expect(run.external).toEqual([]);
    expect(run.errors).toEqual([]);
  } finally {
    await run.context.close();
  }
});

// OV-269, OV-273: a Session 40 the CLI of main 8228ffab wrote (no composition
// metadata in its Result) loads with its tables, and Generate draws the
// case's current CLI figure.
test('Session 40 early CLI linear tables survives Web Generate', async ({ browser }, testInfo) => {
  test.setTimeout(1_800_000);
  const entry = cases.find(item => item.name === 'linear tables');
  const prefix = testInfo.outputPath('cli');
  await promisify(execFile)('python', ['-m', 'gbdraw.cli', entry.mode, ...entry.args, '-o', prefix, '-f', 'svg'], {
    cwd: testInfo.outputDir, env: { ...process.env, PYTHONPATH: root }, timeout: 1_800_000, maxBuffer: 1_000_000
  });
  const file = testInfo.outputPath('v40-early.gbdraw-session.json');
  await fs.writeFile(file, gunzipSync(await fs.readFile(
    path.join(root, 'tests/fixtures/sessions/cli-linear-tables.v40-early.gbdraw-session.json.gz'))));
  const run = await load(browser, file, { width: 1600, height: 1000 });
  try {
    const tables = () => run.page.evaluate(async () => {
      const drawing = (await import('./js/state.js')).state.activeDrawing();
      return [drawing.manualSpecificRules.length, drawing.currentColors.value.CDS, drawing.featureVisibilityRules.value.length,
        drawing.canonicalLabelOverrideRows.value.length, drawing.filterMode.value, drawing.manualWhitelist.length,
        drawing.manualPriorityRules.length];
    });
    // The rows of colors.tsv, default_colors.tsv, visibility.tsv, labels.tsv, whitelist.tsv, priority.tsv.
    const loaded = [2, '#123456', 2, 2, 'Whitelist', 3, 1];
    expect(await tables()).toEqual(loaded);
    await generateAndWaitForResult(run.page);
    expect(await tables()).toEqual(loaded);
    const generated = (await snapshot(run.page)).selected;
    expect(await run.page.evaluate(svgSemantics, generated))
      .toEqual(await run.page.evaluate(svgSemantics, await fs.readFile(`${prefix}.svg`, 'utf8')));
    expect(run.errors).toEqual([]);
  } finally {
    await run.context.close();
  }
});

test('Session schema-1 Web binding control still loads', async ({ browser }, testInfo) => {
  const run = await load(browser, path.join(root, 'tests/fixtures/sessions/single.v41-bindings1.json'), { width: 1600, height: 1000 }, false);
  try {
    expect((await snapshot(run.page)).files).toHaveLength(1);
    expect(run.errors).toEqual([]);
  } finally {
    await run.context.close();
  }
});

// Issue #1427: a CLI Session whose request carries Legend row edits (renames, a hidden row, an edited order of the
// others) and canvas padding, and no Web editor state for them, opens with the Legend editor and Canvas padding
// showing them. An edit then keeps the rows Python drew, and a padding edit pads the canvas by its own change only.
test('Session CLI Legend row edits and canvas padding open in the Legend editor and Canvas padding', async ({ browser }, testInfo) => {
  test.setTimeout(600_000);
  const rows = testInfo.outputPath('legend-rows.tsv');
  await fs.writeFile(rows, 'key\tname\thidden\tposition\ntRNA\t\ttrue\t\nCDS\tGenes\t\t1\nrRNA\tRibosomal RNA\t\t2\n');
  const hmmt = path.join(root, 'tests/test_inputs/HmmtDNA.gbk');
  const file = await cli('circular', ['--gbk', hmmt, '--legend_rows', rows, '--canvas_padding', '10,20,30,40'], testInfo, 'legend-rows');
  const { context, page, errors } = await load(browser, file, { width: 1400, height: 1000 });
  const read = () => page.evaluate(async () => {
    const { state } = await import('./js/state.js');
    const { resultLegendRows } = await import('./js/services/legend-svg.js');
    const drawing = state.activeDrawing();
    const svg = state.svgContainer.value?.querySelector('svg');
    return {
      entries: drawing.legendEntries.value.map((entry) => [entry.originalCaption, entry.caption]),
      deleted: drawing.deletedLegendEntries.value.map((entry) => entry.originalCaption || entry.caption),
      padding: { ...drawing.canvasPadding },
      shown: (resultLegendRows(svg) || []).filter((row) => row.shown).map((row) => row.shownKey),
      size: ['width', 'height'].map((side) => Number.parseFloat(svg.getAttribute(side)))
    };
  });
  const loaded = await read();
  const drawnRows = ['Genes', 'Ribosomal RNA', 'GC content', 'GC skew (+)', 'GC skew (-)'];
  expect(loaded.entries.slice(0, 2)).toEqual([['CDS', 'Genes'], ['rRNA', 'Ribosomal RNA']]);
  expect(loaded.deleted).toEqual(['tRNA']);
  expect(loaded.padding).toEqual({ top: 10, right: 20, bottom: 30, left: 40 });
  expect(loaded.shown).toEqual(drawnRows);
  // An edit of another Legend row shows the drawing's edits again: the rows Python drew stay.
  await page.evaluate(() => {
    const app = window.__GBDRAW_APP__;
    return app.updateLegendEntryColor(app.legendEntries.findIndex((entry) => entry.caption === 'GC content'), '#7b2cbf');
  });
  await settleLive(page);
  const recolored = await read();
  expect(recolored.shown).toEqual(drawnRows);
  expect(recolored.deleted).toEqual(['tRNA']);
  expect(recolored.size).toEqual(loaded.size);
  // A padding edit grows the canvas by its own change, never by less the padding Python drew.
  await page.evaluate(() => window.__GBDRAW_APP__.setCanvasPaddingSide('top', 15));
  await settleLive(page);
  const padded = await read();
  expect(padded.padding.top).toBe(15);
  expect(padded.size[0]).toBeCloseTo(loaded.size[0], 3);
  expect(padded.size[1]).toBeCloseTo(loaded.size[1] + 5, 3);
  expect(errors).toEqual([]);
  await context.close();
});
