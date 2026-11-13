const { test, expect } = require('@playwright/test');
const fs = require('node:fs/promises');
const path = require('node:path');
const { createHash } = require('node:crypto');
const { gunzipSync } = require('node:zlib');
const { execFile } = require('node:child_process');
const { promisify } = require('node:util');
const { CURRENT_REQUEST_SCHEMA, CURRENT_SESSION_VERSION, evaluateWithRetainedPromise, generateAndWaitForResult, openApp, expectSessionLoadedWithoutWarnings, settleLive } = require('./helpers/app-lifecycle.cjs');
const { shownSavedDrawing } = require('./helpers/saved-session.cjs');
const { readSessionFile, readSessionFileJson, sessionFileBytesOf } = require('./helpers/session-file.cjs');
const { drawingIds, shownDrawing } = require('./helpers/drawing-bar.cjs');
const { chooseCoordinates, openRegionDialog } = require('./helpers/region-drawing.cjs');
const { popupEdit } = require('./helpers/live-generate-parity-steps.cjs');

const root = process.cwd();
const mito = path.join(root, 'tests/fixtures/sessions/cli-web-mito.gb');
const lambda = path.join(root, 'tests/test_inputs/NC_001416.gb');
const gff = path.join(root, 'tests/test_inputs/NC_013668.gff3');
const fasta = path.join(root, 'tests/test_inputs/NC_013668.fasta');
const hash = bytes => createHash('sha256').update(bytes).digest('hex');
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
  const file = `${prefix}.gbdraw`;
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
  const { sessionResourceStream } = await import('./js/services/session-resource-backing.js');
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
          // The bytes the descriptor holds: a loaded container's entry copy, or a base64 text.
          sha256: await digest(await new Response(await sessionResourceStream(part.descriptor)).arrayBuffer())
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
  await page.locator('input[accept*=".gbdraw"]').setInputFiles(file);
  await expectSessionLoadedWithoutWarnings(page, { timeout: 180_000 });
  await expect.poll(() => page.evaluate(() => !window.__GBDRAW_APP__.sessionImportPending)).toBe(true);
  // On a phone a Load with a Result shows the Preview pane, as Generate does.
  if (requirePreview) await expect(page.locator('.origin-top svg')).toBeVisible();
  return { context, page, external, errors, dialogs };
};
// A further Load into the page of `load`: it is done when the shown Results are no longer the ones before it.
const loadInto = async (run, file) => {
  const shown = () => run.page.evaluate(() => window.__GBDRAW_APP__.currentResults().map(result => String(result.content || '')).join('\n'));
  const before = await shown();
  await run.page.locator('input[accept*=".gbdraw"]').setInputFiles(file);
  await expect.poll(async () => !(await run.page.evaluate(() => window.__GBDRAW_APP__.sessionImportPending)) && (await shown()) !== before,
    { timeout: 180_000 }).toBe(true);
  await expectSessionLoadedWithoutWarnings(run.page, { timeout: 180_000 });
  await expect(run.page.locator('.origin-top svg')).toBeVisible();
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
      const session = await readSessionFile(file);
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
          webFile = testInfo.outputPath('web-current.gbdraw');
          await (await saved).saveAs(webFile);
          const webSession = await readSessionFile(webFile);
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

// The CLI arguments of one record (BGC0000708) drawn in two windows against BGC0000709, with a LOSATP comparison.
const windowedLosatpArgs = (mode) => {
  const bgc = id => path.join(root, `tests/test_inputs/${id}.gbk`);
  return ['--gbk', bgc('BGC0000708'), bgc('BGC0000709'), '--region', 'BGC0000708:1-15000', '--region',
    'BGC0000708:25000-40000', '--region', 'BGC0000709:1-30000', '--losat', 'losatp',
    ...(mode ? ['--losatp_mode', mode] : []), '-f', 'svg'];
};

// OV-2995: a CLI Session that draws one record in two windows (`--region` twice) with a LOSAT comparison names
// drawn units, one per window, in its comparisons. The Web offers Inherit for it (its record indexes cannot hold
// window pairs), and Generate draws the saved figure with the saved unit pairs. A Pairwise run stores the tables it
// drew as Similarity groups does, so Inherit draws them without LOSAT (OV-3037).
// Returns the figure the CLI drew, read from the Load (`svgSemantics`).
const windowedInherit = async (browser, file) => {
  const saved = shownSavedDrawing(await readSessionFile(file)).renderRequest;
  const unitPairs = request => request.comparisons.filter(comparison => comparison.kind === 'precomputedProteinComparison')
    .map(comparison => [comparison.queryRecordIndex, comparison.subjectRecordIndex]);
  expect(saved.records.map(record => record.windows.length)).toEqual([2, 1]);
  expect(unitPairs(saved)).toEqual([[0, 2], [1, 2]]);
  const run = await load(browser, file, { width: 1600, height: 1000 });
  try {
    const before = await snapshot(run.page);
    const expectedSvg = await run.page.evaluate(svgSemantics, before.selected);
    // The stored tables draw: no LOSAT search runs.
    await run.page.evaluate(() => {
      window.__GBDRAW_LOSAT_EXECUTOR__ = async () => { throw new Error('LOSAT ran'); };
    });
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
    return expectedSvg;
  } finally {
    await run.context.close();
  }
};
test('Session CLI windowed record with LOSAT regenerates through Inherit', async ({ browser }, testInfo) => {
  test.setTimeout(900_000);
  await windowedInherit(browser, await cli('linear', windowedLosatpArgs(null), testInfo, 'windowed'));
});

// OV-3047: a CLI Session of whole records stores the LOSATP results it drew and no Web controls. It loads with the
// LOSATP mode the CLI drew, and the first Generate draws the stored results, with no LOSAT run, as the CLI drew them.
for (const mode of ['pairwise', 'collinear']) test(`Session CLI LOSATP ${mode} loads in its mode and draws its stored results`, async ({ browser }, testInfo) => {
  test.setTimeout(900_000);
  const bgc = id => path.join(root, `tests/test_inputs/${id}.gbk`);
  const file = await cli('linear', ['--gbk', bgc('BGC0000708'), bgc('BGC0000709'), '--losat', 'losatp',
    '--losatp_mode', mode, '-f', 'svg'], testInfo, `whole-${mode}`);
  const run = await load(browser, file, { width: 1600, height: 1000 });
  try {
    const expectedSvg = await run.page.evaluate(svgSemantics, (await snapshot(run.page)).selected);
    expect(await run.page.evaluate(() => window.__GBDRAW_APP__.losat.blastp.mode)).toBe(mode);
    await run.page.evaluate(() => {
      window.__GBDRAW_LOSAT_EXECUTOR__ = async () => { throw new Error('LOSAT ran'); };
    });
    await generateAndWaitForResult(run.page);
    expect(await run.page.evaluate(svgSemantics, (await snapshot(run.page)).selected)).toEqual(expectedSvg);
    expect(run.external).toEqual([]);
    expect(run.errors).toEqual([]);
  } finally {
    await run.context.close();
  }
});

// OV-3003, OV-3037: a Session whose LOSATP Pairwise pipeline holds no table and no LOSAT cache: one that names the
// window pairs (as Python's locate writes it), and one that names none (as the CLI wrote a Pairwise run before it
// stored its tables). Inherit runs it in the browser, for the named pairs or else the records of adjacent rows, and
// draws what the CLI drew from its own search.
// `expectedSvg`: the figure of the CLI Session `file`, as its Load shows it.
const windowedPipeline = async (browser, testInfo, file, expectedSvg, pairs) => {
  const session = await readSessionFileJson(file);
  const drawing = shownSavedDrawing(session);
  const tables = drawing.renderRequest.comparisons.filter(comparison => comparison.kind === 'precomputedProteinComparison');
  const [marker] = drawing.renderRequest.comparisons.filter(comparison => comparison.kind === 'generatedProteinComparison');
  // The CLI searched each window of BGC0000708 against BGC0000709: units 0 and 1 against unit 2.
  expect(tables.map(table => [table.queryRecordIndex, table.subjectRecordIndex])).toEqual([[0, 2], [1, 2]]);
  expect(marker.mode).toBe('none');
  // The tables leave the document, so the written file holds no entry of them.
  tables.forEach(table => { delete session.resources[table.resourceId]; });
  drawing.renderRequest.comparisons = [{
    ...marker, mode: 'pairwise',
    pairs: pairs.map(([queryRecordIndex, subjectRecordIndex]) => ({ queryRecordIndex, subjectRecordIndex }))
  }];
  drawing.losatCache = { entries: [] };
  drawing.losatDerivedCache = { entries: [] };
  const unrendered = testInfo.outputPath(`windowed-pipeline-${pairs.length}.gbdraw`);
  await fs.writeFile(unrendered, await sessionFileBytesOf(session));
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
};
// One CLI run of the Pairwise windows serves the three: its Load is the Inherit case's, and the figure that Load
// shows is what the two pipeline Sessions made from it must draw.
test('Session CLI windowed record with LOSAT pairwise regenerates through Inherit, and windowed LOSATP pipelines that name their pairs or none run through Inherit', async ({ browser }, testInfo) => {
  test.setTimeout(900_000);
  const file = await cli('linear', windowedLosatpArgs('pairwise'), testInfo, 'windowed-pairwise');
  const expectedSvg = await test.step('Session CLI windowed record with LOSAT pairwise regenerates through Inherit',
    () => windowedInherit(browser, file));
  for (const [name, pairs] of [['names its pairs', [[0, 2], [1, 2]]], ['names no pairs', []]]) {
    await test.step(`Session with a windowed LOSATP pipeline that ${name} and no results runs it through Inherit`,
      () => windowedPipeline(browser, testInfo, file, expectedSvg, pairs));
  }
});

// OV-2732: a Session file container whose drawing the CLI wrote (no Web draft)
// with a Specific colors table and a Label whitelist loads: Load reads the
// tables from the input copies, the editors show their rows, and Generate draws
// the colors of the saved Result.
test('Session file container of a CLI drawing with a color table and a whitelist loads and draws its colors', async ({ browser }, testInfo) => {
  test.setTimeout(1_800_000);
  const tables = path.join(root, 'tests/fixtures/cli_session_cross_surface');
  const prefix = testInfo.outputPath('container-tables');
  await promisify(execFile)('python', ['-m', 'gbdraw.cli', 'linear', '--gbk', lambda, '--show_labels', 'all',
    '-t', path.join(tables, 'colors.tsv'), '--label_whitelist', path.join(tables, 'whitelist.tsv'),
    '-o', prefix, '--session_output', `${prefix}.gbdraw`], {
    cwd: testInfo.outputDir, env: { ...process.env, PYTHONPATH: root }, timeout: 1_800_000, maxBuffer: 1_000_000
  });
  const file = `${prefix}.gbdraw`;
  // The selected Result, the mounted figure and the error, read without the inputs' bytes.
  const figures = page => page.evaluate(async () => {
    const { state } = await import('./js/state.js');
    return { selected: window.__GBDRAW_APP__.currentResults()[state.selectedResultIndex.value]?.content || '',
      mounted: state.svgContainer.value?.querySelector('svg')?.outerHTML || '', error: state.errorLog.value };
  });
  const run = await load(browser, file, { width: 1600, height: 1000 });
  try {
    const before = await figures(run.page);
    expect(before.error).toBeNull();
    const saved = await run.page.evaluate(svgSemantics, before.selected);
    // The saved Result draws the portal protein in the table's color.
    expect(saved.features.some(([, , , fill]) => fill === '#ff0000')).toBe(true);
    expect(await run.page.evaluate(svgSemantics, before.mounted)).toEqual(saved);
    expect(await run.page.evaluate(async () => {
      const { state } = await import('./js/state.js');
      const drawing = state.activeDrawing();
      return { rules: drawing.manualSpecificRules.map(rule => [rule.feat, rule.qual, rule.val, rule.color]),
        whitelist: drawing.manualWhitelist.map(row => [row.feat, row.qual, row.key]) };
    })).toEqual({
      rules: [['CDS', 'product', 'portal', '#ff0000'], ['CDS', 'gene', '^COX', '#ff0000']],
      whitelist: [['CDS', 'product', 'terminase'], ['CDS', 'gene', 'COX'], ['CDS', 'gene', 'CYTB']]
    });
    await generateAndWaitForResult(run.page);
    const after = await figures(run.page);
    expect(after.error).toBeNull();
    expect(await run.page.evaluate(svgSemantics, after.selected)).toEqual(saved);
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
  const run = await load(browser, path.join(root, 'gbdraw/web/gallery/sessions/HmmtDNA_basic_circular.gbdraw'),
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
    const webFile = testInfo.outputPath('web-region.gbdraw');
    await (await saved).saveAs(webFile);
    const session = await readSessionFile(webFile);
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

// OV-3031: a CLI Session with two feature density tracks keeps them at Web Generate: the new Result draws the same
// density slot groups (paths and fills) as the CLI's replay of that Session. The Result the Session stores keeps the
// layout of the commit that wrote the fixture, so it gives only the slots and the fills.
test('Session CLI feature density tracks survive Web Generate', async ({ browser }, testInfo) => {
  test.setTimeout(600_000);
  const file = path.join(root, 'tests/fixtures/sessions/HmmtDNA-feature-density.v51-cli.gbdraw');
  const prefix = testInfo.outputPath('feature-density-replay');
  await promisify(execFile)('python', ['-m', 'gbdraw.cli', 'circular', '--session', file, '-o', prefix, '-f', 'svg'], {
    cwd: testInfo.outputDir, env: { ...process.env, PYTHONPATH: root }, timeout: 600_000, maxBuffer: 1_000_000
  });
  const run = await load(browser, file, { width: 1600, height: 1000 });
  try {
    const density = (svg) => run.page.evaluate((text) => {
      const root = new DOMParser().parseFromString(text, 'image/svg+xml');
      return [...root.querySelectorAll('[data-gbdraw-slot-renderer="feature_density"]')].map((group) => [
        group.getAttribute('data-gbdraw-slot-id'),
        [...group.querySelectorAll('path')].map((item) => `${item.getAttribute('fill')}|${item.getAttribute('d')}`).sort()
      ]);
    }, svg);
    const fills = (groups) => groups.map(([id, paths]) => [id, paths.map((item) => item.split('|')[0])]);
    const stored = await density((await readSessionFileJson(file)).drawings[0].results[0].content);
    expect(stored.map(([id]) => id)).toEqual(['density_genes', 'density_trna']);
    const replayed = await density(await fs.readFile(`${prefix}.svg`, 'utf8'));
    expect(fills(replayed)).toEqual(fills(stored));
    await generateAndWaitForResult(run.page);
    expect(await density((await snapshot(run.page)).selected)).toEqual(replayed);
    expect(run.errors).toEqual([]);
  } finally {
    await run.context.close();
  }
});

test('Session CLI one circle loads its layout, Generates the CLI figure, saves it, and Reset returns to the grid', async ({ browser }, testInfo) => {
  test.setTimeout(900_000);
  const file = await cli('circular', ['--gbk', mito, lambda, '--one_circle', '--one_circle_gap_degrees', '2.5',
    '--one_circle_no_record_band'], testInfo, 'one-circle');
  const session = await readSessionFile(file);
  expect(shownSavedDrawing(session).renderRequest.layout).toMatchObject({
    recordLayout: 'one_circle', oneCircleGapDegrees: 2.5, oneCircleRecordBand: false
  });
  const run = await load(browser, file, { width: 1600, height: 1000 });
  try {
    const settings = () => run.page.evaluate(() => {
      const app = window.__GBDRAW_APP__;
      const checked = (selector) => document.querySelector(selector)?.checked;
      return {
        canvas: app.form.multi_record_canvas, layout: app.adv.circular_record_layout, gap: app.adv.one_circle_gap_degrees,
        minLength: app.adv.one_circle_min_length, band: app.adv.one_circle_record_band,
        grid: checked('[data-circular-canvas-setting]'), oneCircle: checked('[data-circular-one-circle-setting]'),
        shown: Boolean(document.querySelector('[data-capture="one-circle-controls"]'))
      };
    });
    expect(await settings()).toEqual({
      canvas: true, layout: 'one_circle', gap: 2.5, minLength: 1000, band: false, grid: false, oneCircle: true, shown: true
    });
    // OV-3490: the One circle fields of a one-circle Result apply live, so their block shows no clock.
    await expect(run.page.locator('[data-capture="one-circle-controls"] [data-timing-clock]')).toHaveCount(0);
    const stored = shownSavedDrawing(session).results[0].content;
    await generateAndWaitForResult(run.page);
    const generated = (await snapshot(run.page)).selected;
    expect(await run.page.evaluate(svgSemantics, generated)).toEqual(await run.page.evaluate(svgSemantics, stored));

    const saved = run.page.waitForEvent('download');
    const outcome = await evaluateWithRetainedPromise(run.page, async () => {
      window.__GBDRAW_APP__.sessionTitle = 'one-circle';
      return await window.__GBDRAW_APP__.saveSessionWithTitle();
    });
    expect(outcome.status).toBe('saved');
    const webFile = testInfo.outputPath('web-one-circle.gbdraw');
    await (await saved).saveAs(webFile);
    const webDrawing = shownSavedDrawing(await readSessionFile(webFile));
    expect(webDrawing.renderRequest.layout).toMatchObject({
      recordLayout: 'one_circle', oneCircleGapDegrees: 2.5, oneCircleMinLength: 1000, oneCircleRecordBand: false
    });
    expect(webDrawing.config.adv).toMatchObject({
      circular_record_layout: 'one_circle', one_circle_gap_degrees: 2.5, one_circle_record_band: false
    });

    await evaluateWithRetainedPromise(run.page, async () => { await window.__GBDRAW_APP__.resetSettings(); });
    await settleLive(run.page);
    expect(await settings()).toEqual({
      canvas: true, layout: 'grid', gap: 1, minLength: 1000, band: true, grid: true, oneCircle: false, shown: false
    });
    expect(run.errors).toEqual([]);
  } finally {
    await run.context.close();
  }
});

// OV-3079, OV-3155: a CLI stack that names no Axis index draws each slot on its own side (a set's row without a side
// outside or above, a Linear features slot without a side below). The Session it saves, loaded and generated in the
// Web, draws the CLI's figure: the region set stays outside (above) the Axis, and Overlay set rows paint where Python
// paints them.
const NO_AXIS_STACKS = [
  ['circular', 'circular', ['features:features', 'regions:annotations@set_id=regions,side=outside'],
    1, [['regions', 'outside'], ['features', 'inside']]],
  ['linear', 'linear', ['features:features', 'regions:annotations@set_id=regions,side=above'],
    1, [['regions', 'above'], ['features', 'below']]],
  ['circular blank-side', 'circular', ['hl:annotations@set_id=hl,side=overlay,anchor_slot=features,layer=underlay,z=-1',
    'features:features', 'regions:annotations@set_id=regions'], 2, [['hl', 'overlay'], ['regions', 'outside'], ['features', 'inside']]],
  ['circular overlay tie', 'circular', ['features:features', 'a:annotations@set_id=a,side=overlay,anchor_slot=features,z=1',
    'gc:content', 'regions:annotations@set_id=regions,side=outside', 'b:annotations@set_id=b,side=overlay,anchor_slot=features,z=1'],
  2, [['a', 'overlay'], ['regions', 'outside'], ['features', 'inside'], ['gc', 'inside'], ['b', 'overlay']]],
  // A foreground overlay paints over every row of its z it reaches, wherever it was written (OV-3184, OV-3185), so it
  // keeps its place before the row after it. `a` spans features to `sk` and covers `gc` of its z.
  ['circular overlay span', 'circular', ['features:features', 'gc:content@z=1',
    'a:annotations@set_id=a,side=overlay,anchor_slot=features,to_slot=sk,z=1', 'regions:annotations@set_id=regions,side=outside,z=1',
    'sk:skew'], 2, [['a', 'overlay'], ['regions', 'outside'], ['features', 'inside'], ['gc', 'inside'], ['sk', 'inside']], ['gc']],
  // `a` with its own width reaches `sk` of its z past its anchor.
  ['circular overlay wide', 'circular', ['features:features', 'gc:content', 'sk:skew@z=1',
    'a:annotations@set_id=a,side=overlay,anchor_slot=gc,w=120px,z=1', 'regions:annotations@set_id=regions,side=outside,z=1'],
  2, [['a', 'overlay'], ['regions', 'outside'], ['features', 'inside'], ['gc', 'inside'], ['sk', 'inside']], ['sk']],
  // A span over the Axis covers `y` of its z inside, written before it, and `x` outside, written after it.
  ['circular overlay across the Axis', 'circular', ['features:features', 'y:content@z=1',
    'a:annotations@set_id=a,side=overlay,anchor_slot=oc,to_slot=t,z=1', 'oc:content@side=outside', 'x:skew@side=outside,z=1', 't:skew'],
  3, [['a', 'overlay'], ['oc', 'outside'], ['x', 'outside'], ['features', 'inside'], ['y', 'inside'], ['t', 'inside']], ['y', 'x']],
  // The span `a` keeps its place after `b`, written before it, although the row after `a` lies on the other side:
  // region labels keep off the labels of the slots before theirs.
  ['linear overlay after an Annotations row', 'linear', ['features:features', 'b:annotations@set_id=b,side=below,z=1',
    'a:annotations@set_id=a,side=overlay,anchor_slot=features,to_slot=sk,z=1', 'up:skew@side=above', 'sk:skew@side=below'],
  1, [['up', 'above'], ['features', 'below'], ['b', 'below'], ['a', 'overlay'], ['sk', 'below']], ['b']],
  ['circular overlay after an Annotations row', 'circular', ['features:features', 'b:annotations@set_id=b,side=inside,z=1',
    'a:annotations@set_id=a,side=overlay,anchor_slot=features,to_slot=sk,z=1', 'up:skew@side=outside', 'sk:skew'],
  1, [['up', 'outside'], ['features', 'inside'], ['b', 'inside'], ['a', 'overlay'], ['sk', 'inside']], ['b']],
  // Rows written out of side order: the Legend lists them by side in the CLI as in the Web (OV-3187).
  ['circular track rows out of side order', 'circular', ['features:features', 'gc:content', 'sk:skew@side=outside'],
    1, [['sk', 'outside'], ['features', 'inside'], ['gc', 'inside']]],
  ['linear Annotations rows out of side order', 'linear', ['features:features', 'b:annotations@set_id=b,side=below',
    'r:annotations@set_id=r,side=above'], 1, [['r', 'above'], ['features', 'below'], ['b', 'below']]]
];
// One stack: the CLI draws it, and its Session loads into `run` (a page of its own when `run.page` is unset).
const noAxisStack = async (browser, run, testInfo, [, mode, slots, axis, stack, coveredOfItsZ = []]) => {
  await fs.writeFile(testInfo.outputPath('regions.tsv'), [
    'set_id\tid\tmark\tstart\tend\tlabel\tlegend_label',
    'regions\thead\tline\t1\t9500\tHead genes\tRegion set', 'regions\ttail\tline\t9600\t19000\tTail genes\tRegion set',
    ...Object.entries({
      hl: 'hl\th1\thighlight\t20000\t26000\t\tHighlight set', a: 'a\ta1\thighlight\t1000\t9000\t\tAlpha set',
      b: 'b\tb1\thighlight\t5000\t12000\t\tBeta set', r: 'r\tr1\thighlight\t14000\t20000\t\tRho set'
    })
      .filter(([set]) => slots.some((slot) => slot.includes(`set_id=${set},`))).map(([, row]) => row)
  ].join('\n') + '\n');
  // The CLI writes no output over an earlier stack's: the steps of one case share its output folder.
  await Promise.all(['no-axis.gbdraw', 'no-axis.svg'].map((name) => fs.rm(testInfo.outputPath(name), { force: true })));
  const file = await cli(mode, ['--gbk', lambda, '--annotation_table', 'regions.tsv',
    ...slots.flatMap((slot) => [`--${mode}_track_slot`, slot]), '-f', 'svg'], testInfo, 'no-axis');
  const stored = shownSavedDrawing(await readSessionFile(file));
  const [slotsKey, axisKey] = [`${mode}TrackSlots`, `${mode}TrackAxisIndex`];
  expect(stored.renderRequest.diagramOptions.tracks[axisKey]).toBe(null);
  if (run.page) await loadInto(run, file);
  else Object.assign(run, await load(browser, file, { width: 1600, height: 1000 }));
  // The figure's size and each slot group's place and paths, by slot id (4 decimals, as `svgSemantics`).
  const slotGroups = (svg) => run.page.evaluate((text) => {
    const round = (value) => value && value.replace(/-?\d*\.?\d+(?:e[-+]?\d+)?/gi, (number) => String(Math.round(Number(number) * 1e4) / 1e4));
    const root = new DOMParser().parseFromString(text, 'image/svg+xml').documentElement;
    return {
      size: [round(root.getAttribute('width')), round(root.getAttribute('height'))],
      slots: [...root.querySelectorAll('[data-gbdraw-slot-id]')].map((group) => [
        group.getAttribute('data-gbdraw-slot-id'), round(group.getAttribute('transform')),
        [...group.querySelectorAll('path')].map((item) => round(item.getAttribute('d'))).sort()
      ]).sort((a, b) => a[0].localeCompare(b[0]))
    };
  }, svg);
  // The Web stack lists the rows outside (above) the Axis first, so its slot groups and their texts come in that order.
  const semantics = async (svg) => {
    const read = await run.page.evaluate(svgSemantics, svg);
    return { ...read, text: [...read.text].sort() };
  };
  // Overlay set rows overlap each other and the rows of their z they reach: whether each overlay's group comes after
  // each of those in the document is whether Python paints it over them.
  const overlays = stack.filter(([, drawnSide]) => drawnSide === 'overlay').map(([id]) => id);
  const painted = (svg) => run.page.evaluate(([text, overlayIds, otherIds]) => {
    const order = [...new DOMParser().parseFromString(text, 'image/svg+xml').querySelectorAll('[data-gbdraw-slot-id]')]
      .map((group) => group.getAttribute('data-gbdraw-slot-id'));
    return overlayIds.flatMap((overlay) => [...overlayIds, ...otherIds].filter((id) => id !== overlay)
      .map((id) => [overlay, id, order.indexOf(overlay) > order.indexOf(id)]));
  }, [svg, overlays, coveredOfItsZ]);
  // The Legend's entries, in the order it lists them.
  const legend = (svg) => run.page.evaluate((text) => [...new DOMParser().parseFromString(text, 'image/svg+xml')
    .querySelectorAll('#legend text')].map((item) => item.textContent), svg);
  const drawn = await slotGroups(stored.results[0].content);
  const before = await snapshot(run.page);
  expect(await slotGroups(before.selected)).toEqual(drawn);
  await generateAndWaitForResult(run.page);
  const after = await snapshot(run.page);
  expect(await slotGroups(after.selected)).toEqual(drawn);
  expect(await painted(after.selected)).toEqual(await painted(stored.results[0].content));
  expect(await legend(after.selected)).toEqual(await legend(stored.results[0].content));
  expect(await semantics(after.selected)).toEqual(await semantics(before.selected));
  const tracks = after.request.diagramOptions.tracks;
  expect({ axis: tracks[axisKey], slots: tracks[slotsKey].map(({ id, side: drawnSide }) => [id, drawnSide]) })
    .toEqual({ axis, slots: stack });
  expect(run.errors).toEqual([]);
};
// The stacks of one mode share a page: each Load replaces the drawing (a Load replaces the whole project), and each
// step names its stack.
for (const mode of ['circular', 'linear']) {
  test(`Session CLI ${mode} stacks without an Axis index draw their region sets where the CLI drew them`, async ({ browser }, testInfo) => {
    test.setTimeout(900_000);
    const run = {};
    try {
      for (const entry of NO_AXIS_STACKS.filter(([, entryMode]) => entryMode === mode)) {
        await test.step(`Session CLI ${entry[0]} stack without an Axis index draws its region set where the CLI drew it`,
          () => noAxisStack(browser, run, testInfo, entry));
      }
    } finally {
      await run.context?.close();
    }
  });
}

test('Session schema-1 Web binding control still loads', async ({ browser }, testInfo) => {
  const run = await load(browser, path.join(root, 'tests/fixtures/sessions/single.v41-bindings1.json'), { width: 1600, height: 1000 }, false);
  try {
    expect((await snapshot(run.page)).files).toHaveLength(1);
    expect(run.errors).toEqual([]);
  } finally {
    await run.context.close();
  }
});

// OV-3189: a CLI stack that puts the scale above the feature track and names no tick label layout draws its labels
// outside the ticks, as the Web orients such a row (both read the layout by place that Python owns). The Session it
// saves, loaded and generated in the Web, draws the CLI's figure, and the labels stand next to their ticks, clear of
// the feature track.
test('Session CLI scale above the features without a tick label layout draws the CLI figure', async ({ browser }, testInfo) => {
  test.setTimeout(600_000);
  const hmmt = path.join(root, 'tests/test_inputs/HmmtDNA.gbk');
  const file = await cli('circular', ['--gbk', hmmt, '-f', 'svg', '--circular_track_slot', 'ticks:ticks',
    '--circular_track_slot', 'features:features', '--circular_track_slot', 'gc_content:dinucleotide_content'], testInfo, 'ticks-first');
  const stored = shownSavedDrawing(await readSessionFile(file));
  expect(stored.renderRequest.diagramOptions.tracks.circularTrackSlots.find(({ id }) => id === 'ticks').params).toEqual({});
  const run = await load(browser, file, { width: 1600, height: 1000 });
  try {
    // Each slot group's place and paths, by slot id (4 decimals, as `svgSemantics`).
    const slotGroups = (svg) => run.page.evaluate((text) => {
      const round = (value) => value && value.replace(/-?\d*\.?\d+(?:e[-+]?\d+)?/gi, (number) => String(Math.round(Number(number) * 1e4) / 1e4));
      return [...new DOMParser().parseFromString(text, 'image/svg+xml').querySelectorAll('[data-gbdraw-slot-id]')]
        .map((group) => [group.getAttribute('data-gbdraw-slot-id'), round(group.getAttribute('transform')),
          [...group.querySelectorAll('path')].map((item) => round(item.getAttribute('d'))).sort()])
        .sort((a, b) => a[0].localeCompare(b[0]));
    }, svg);
    // Radii from the drawing: the tick marks' outer ends, the radius of the label paths (the labels stand outward
    // of it), and the feature track's outer edge (its largest arc).
    const radii = (svg) => run.page.evaluate((text) => {
      const root = new DOMParser().parseFromString(text, 'image/svg+xml');
      const numbers = (d) => (d.match(/-?\d*\.?\d+/g) || []).map(Number);
      const ticks = root.querySelector('[data-gbdraw-slot-id="ticks"]');
      const marks = [...ticks.querySelectorAll('path:not([id])')].map((item) => numbers(item.getAttribute('d')))
        .map(([x1, y1, x2, y2]) => Math.max(Math.hypot(x1, y1), Math.hypot(x2, y2)));
      const arcs = (items) => [...items].flatMap((item) => [...item.getAttribute('d').matchAll(/A\s*(-?[\d.]+)/g)].map((m) => Number(m[1])));
      return {
        tickOuter: Math.max(...marks),
        labelPaths: [...new Set(arcs(ticks.querySelectorAll('path[id^="circular_tick_label_path"]')))],
        featureOuter: Math.max(...arcs(root.querySelectorAll('[data-gbdraw-slot-id="features"] path')))
      };
    }, svg);
    const drawn = await slotGroups(stored.results[0].content);
    await generateAndWaitForResult(run.page);
    const after = await snapshot(run.page);
    expect(await slotGroups(after.selected)).toEqual(drawn);
    const { tickOuter, labelPaths, featureOuter } = await radii(after.selected);
    expect(tickOuter).toBeGreaterThan(featureOuter);
    expect(labelPaths.length).toBe(1);
    // The labels start within the label margin past their ticks, outside the feature track.
    expect(labelPaths[0] - tickOuter).toBeGreaterThan(0);
    expect(labelPaths[0] - tickOuter).toBeLessThan(5);
    expect(run.errors).toEqual([]);
  } finally {
    await run.context.close();
  }
});

// OV-3301: with the scale outside the feature track and its labels outside the ticks, a number a feature label's
// leader line runs through is drawn in front of it on a backing. The CLI Session, loaded and generated in the Web,
// draws the CLI's scale with those backings; a label hidden live takes its leader's backing away, as Generate does.
test('Session CLI scale outside the features backs the numbers a leader runs through, live and from Generate', async ({ browser }, testInfo) => {
  test.setTimeout(600_000);
  const hmmt = path.join(root, 'tests/test_inputs/HmmtDNA.gbk');
  const file = await cli('circular', ['--gbk', hmmt, '-f', 'svg', '--circular_track_order', 'ticks,features,gc_content,gc_skew',
    '--circular_track_axis_index', '1', '--labels', 'out', '--label_font_size', '10'], testInfo, 'scale-outside');
  const stored = shownSavedDrawing(await readSessionFile(file));
  const run = await load(browser, file, { width: 1600, height: 1000 });
  try {
    // The scale group's place and paths (4 decimals, as `svgSemantics`), and the features whose leader line enters
    // a backing (an annular sector drawn clockwise from its first point).
    const scale = (svg) => run.page.evaluate((text) => {
      const root = new DOMParser().parseFromString(text, 'image/svg+xml');
      const round = (value) => value && value.replace(/-?\d*\.?\d+(?:e[-+]?\d+)?/gi, (number) => String(Math.round(Number(number) * 1e4) / 1e4));
      const ticks = root.querySelector('[data-gbdraw-slot-id="ticks"]');
      const turn = (x, y) => (Math.atan2(y, x) + 2 * Math.PI) % (2 * Math.PI);
      const backings = [...ticks.querySelectorAll('path[fill-opacity]')].map((item) => {
        const [x1, y1, outer, , , , , x2, y2, , , inner] = (item.getAttribute('d').match(/-?\d*\.?\d+/g) || []).map(Number);
        return { inner, outer, start: turn(x1, y1), span: (turn(x2, y2) - turn(x1, y1) + 2 * Math.PI) % (2 * Math.PI) };
      });
      const inside = ([x, y]) => backings.some(({ inner, outer, start, span }) => Math.hypot(x, y) >= inner
        && Math.hypot(x, y) <= outer && (turn(x, y) - start + 2 * Math.PI) % (2 * Math.PI) <= span);
      const backed = [...root.querySelectorAll('[id="label_leaders"] line')].filter((line) => {
        const [x1, y1, x2, y2] = ['x1', 'y1', 'x2', 'y2'].map((key) => Number(line.getAttribute(key)));
        return Array.from({ length: 101 }, (_, i) => [x1 + (x2 - x1) * i / 100, y1 + (y2 - y1) * i / 100]).some(inside);
      }).map((line) => line.getAttribute('data-label-feature-id'));
      return {
        transform: round(ticks.getAttribute('transform')),
        paths: [...ticks.querySelectorAll('path')].map((item) => round(item.getAttribute('d'))).sort(),
        backed: [...new Set(backed)].sort()
      };
    }, svg);
    const drawn = await scale(stored.results[0].content);
    expect(drawn.backed.length).toBeGreaterThan(0);
    await generateAndWaitForResult(run.page);
    expect(await scale((await snapshot(run.page)).selected)).toEqual(drawn);
    await run.page.evaluate(() => { window.__GBDRAW_APP__.autoLabelReflowEnabled = true; });
    await popupEdit(run.page, { biological_feature_id: drawn.backed[0] }, { labelVisibility: 'off' });
    await settleLive(run.page);
    const live = await scale((await snapshot(run.page)).selected);
    expect(live.backed).not.toContain(drawn.backed[0]);
    await generateAndWaitForResult(run.page);
    expect(await scale((await snapshot(run.page)).selected)).toEqual(live);
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
