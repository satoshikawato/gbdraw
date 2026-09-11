const { test, expect } = require('@playwright/test');
const fs = require('node:fs/promises');
const path = require('node:path');
const { createHash } = require('node:crypto');
const { gunzipSync } = require('node:zlib');
const { execFile } = require('node:child_process');
const { promisify } = require('node:util');
const { openApp, generateAndWaitForResult } = require('./helpers/app-lifecycle.cjs');

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
  { name: 'gff fasta', mode: 'circular', args: ['--gff', gff, '--fasta', fasta], sources: [gff, fasta] }
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
  return {
    records: [...records.entries()].sort((a, b) => a[0] - b[0]),
    features: [...svg.querySelectorAll('[data-gbdraw-feature-id]')].map(element =>
      ['data-gbdraw-feature-id', 'data-gbdraw-record-id', 'd', 'fill', 'stroke', 'display', 'visibility'].map(key => element.getAttribute(key))),
    text: [...svg.querySelectorAll('text')].map(element => element.textContent)
  };
};

const snapshot = page => page.evaluate(async () => {
  const { state } = await import('./js/state.js');
  const { getCommittedCanonicalSession } = await import('./js/services/config.js');
  const { getSessionResourceSource, readFileBytes } = await import('./js/services/file-content-cache.js');
  const digest = async bytes => [...new Uint8Array(await crypto.subtle.digest('SHA-256', bytes))]
    .map(byte => byte.toString(16).padStart(2, '0')).join('');
  const files = state.mode.value === 'linear' ? state.linearSeqs.map(seq => seq.gb)
    : state.cInputType.value === 'gff' ? [state.files.c_gff, state.files.c_fasta] : [state.files.c_gb];
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
    selected: state.results.value[state.selectedResultIndex.value]?.content || '',
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
  const dialog = page.waitForEvent('dialog', { timeout: 180_000 });
  await page.locator('input[accept^=".json,"]').setInputFiles(file);
  expect((await dialog).message()).toBe('Session loaded successfully!');
  await expect.poll(() => page.evaluate(() => !window.__GBDRAW_APP__.sessionImportPending)).toBe(true);
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
      expect(session.version).toBe(41);
      expect(session.renderRequest.schema).toBe(7);
      expect(session.webFiles.bindings.schema).toBe(2);
      if (entry.name === 'composite') {
        expect(session.webFiles.bindings.c_gb.kind).toBe('composite');
        expect(session.webFiles.bindings.c_gb.components).toHaveLength(2);
      }
      const run = await load(browser, file,
        entry.name === 'original single' ? { width: 390, height: 844 } : { width: 1600, height: 1000 });
      try {
        const before = await snapshot(run.page);
        expect(before.mode).toBe(entry.mode);
        expect(before.request).toEqual(session.renderRequest);
        expect(before.files.flatMap(file => file.parts.map(part => part.sha256))).toEqual(sourceHashes);
        expect(before.files.every(file => !file.isArray)).toBe(true);
        const inventory = session.webFiles.bindings;
        const bindings = entry.mode === 'linear' ? inventory.linearSeqs.map(seq => seq.gb)
          : entry.name === 'gff fasta' ? [inventory.c_gff, inventory.c_fasta] : [inventory.c_gb];
        const metadata = file => ({ name: file.name, type: file.type, lastModified: file.lastModified });
        expect(before.files.map(metadata)).toEqual(bindings.map(metadata));
        expect(before.files.flatMap(file => file.parts.map(part => ({ resourceId: part.resourceId, ...metadata(part) }))))
          .toEqual(bindings.flatMap(binding => (binding.components || [binding]).map(part => ({ resourceId: part.resourceId, ...metadata(part) }))));
        for (const part of before.files.flatMap(file => file.parts)) expect(before.resourceIds).toContain(part.resourceId);
        const expectedSvg = await run.page.evaluate(svgSemantics, before.selected);
        expect(expectedSvg.records.length).toBe(session.renderRequest.records.length);
        expect(await run.page.evaluate(svgSemantics, before.mounted)).toEqual(expectedSvg);
        await generateAndWaitForResult(run.page);
        const after = await snapshot(run.page);
        const generatedSvg = await run.page.evaluate(svgSemantics, after.selected);
        expect(generatedSvg.records).toEqual(expectedSvg.records);
        expect(after.request.records).toHaveLength(session.renderRequest.records.length);
        expect(after.files).toEqual(before.files);
        expect(after.files.flatMap(file => file.parts.map(part => part.sha256))).toEqual(sourceHashes);
        expect(await run.page.evaluate(svgSemantics, after.mounted)).toEqual(generatedSvg);
        const download = run.page.waitForEvent('download');
        await run.page.getByRole('button', { name: 'SVG', exact: true }).click();
        const exported = testInfo.outputPath(`${phase}-export.svg`);
        await (await download).saveAs(exported);
        expect(await run.page.evaluate(svgSemantics, await fs.readFile(exported, 'utf8'))).toEqual(generatedSvg);
        if (phase === 'cli-replay') {
          const saved = run.page.waitForEvent('download');
          const outcome = await run.page.evaluate(async () => await window.__GBDRAW_APP__.saveSessionWithTitle());
          expect(outcome.status).toBe('saved');
          webFile = testInfo.outputPath('web-current.gbdraw-session.json.gz');
          await (await saved).saveAs(webFile);
          expect((await readSession(webFile)).config.form).toBeTruthy();
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

test('Session schema-1 Web binding control still loads', async ({ browser }, testInfo) => {
  const run = await load(browser, path.join(root, 'tests/fixtures/sessions/single.v41-bindings1.json'), { width: 1600, height: 1000 }, false);
  try {
    expect((await snapshot(run.page)).files).toHaveLength(1);
    expect(run.errors).toEqual([]);
  } finally {
    await run.context.close();
  }
});
