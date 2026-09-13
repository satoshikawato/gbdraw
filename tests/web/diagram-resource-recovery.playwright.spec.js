const { test, expect } = require('@playwright/test');
const fs = require('node:fs/promises');
const { createHash } = require('node:crypto');
const { openApp } = require('./helpers/app-lifecycle.cjs');
const { snapshot, download } = require('./helpers/mode-transition.cjs');

const gffPath = 'gbdraw/web/tutorial-data/lambda-gff3/NC_001416.gff3';
const fastaPath = 'gbdraw/web/tutorial-data/lambda-gff3/NC_001416.fna';
// Exact retained SESSION 05A4 M23 rejected input (424 bytes).
const mismatch = Buffer.from('>wrong_record_identity\n' + 'ATGC'.repeat(100) + '\n');
const semantics = (page, content) => page.evaluate(content => {
  const root = new DOMParser().parseFromString(content, 'image/svg+xml').documentElement;
  return {
    features: [...root.querySelectorAll('[data-gbdraw-feature-id]')].map(e =>
      Object.fromEntries(['data-gbdraw-feature-id', 'data-gbdraw-record-id', 'data-gbdraw-feature-part',
        'd', 'fill', 'stroke', 'display', 'visibility'].map(key => [key, e.getAttribute(key)]))),
    labels: [...root.querySelectorAll('text')].map(e => ({ text: e.textContent, id: e.id }))
  };
}, content);

const installEvidence = async page => {
  await page.addInitScript(() => {
    const evidence = window.__RESOURCE_RECOVERY__ = { events: [], metrics: [], workers: [], runs: [], rejections: [] };
    window.__MODE_EVENTS__ = [];
    window.__GBDRAW_TEST_HOOKS__ = {
      onStructuralMetric: event => evidence.metrics.push(event),
      onSessionLifecycleEvent: event => {
        evidence.events.push(event);
        if (event.name === 'result-admission-start' && evidence.rejectAdmission) {
          evidence.rejectAdmission = false;
          throw new Error('Injected Result admission rejection');
        }
        if (event.name === 'worker-resource-linking-end' && evidence.cancelAfterStaging) {
          evidence.cancelAfterStaging = false;
          window.__GBDRAW_APP__.cancelGeneration();
        }
      }
    };
    window.addEventListener('unhandledrejection', event => evidence.rejections.push(String(event.reason)));
    const NativeWorker = window.Worker;
    window.Worker = new Proxy(NativeWorker, { construct(target, args) {
      const worker = Reflect.construct(target, args, target);
      if (!String(args[0]).includes('diagram-generation-worker.js')) return worker;
      const entry = { initializations: 0, terminated: false, messages: [] };
      evidence.workers.push(entry);
      const post = worker.postMessage.bind(worker), terminate = worker.terminate.bind(worker);
      worker.postMessage = (message, transfer) => {
        if (message.type === 'init') entry.initializations++;
        if (message.type === 'run') evidence.runs.push({
          requestId: message.requestId, request: message.payload.request,
          manifest: structuredClone(message.payload.resourceManifest),
          staged: message.payload.stagedResources.map(r => ({
            resourceId: r.resourceId, cacheToken: r.cacheToken, bytes: r.bytes.byteLength
          }))
        });
        return transfer === undefined ? post(message) : post(message, transfer);
      };
      worker.terminate = () => { entry.terminated = true; return terminate(); };
      worker.addEventListener('message', ({ data }) => {
        if (data.type !== 'test-lifecycle') entry.messages.push({
          type: data.type, requestId: data.requestId, ok: data.ok,
          error: data.error || data.results?.error || null
        });
      });
      return worker;
    } });
  });
};

test('W1/W2/W4/W8/W9: retained M23 recovery succeeds once and resource ownership survives admission and cancellation', async ({ browser }, testInfo) => {
  test.setTimeout(1_800_000);
  const context = await browser.newContext({ viewport: { width: 1600, height: 1000 }, acceptDownloads: true });
  const external = [], errors = [], outcomes = [];
  await context.route('**/*', route => {
    if (new URL(route.request().url()).hostname === '127.0.0.1') return route.continue();
    external.push(route.request().url());
    return route.abort();
  });
  const page = await context.newPage();
  page.setDefaultTimeout(180_000);
  page.on('pageerror', error => errors.push(String(error)));
  await installEvidence(page);
  const generate = async expected => {
    const result = await page.evaluate(() => window.__GBDRAW_APP__.runAnalysis());
    outcomes.push({ result, error: await page.evaluate(() => window.__GBDRAW_APP__.errorLog) });
    expect(result.status).toBe(expected);
  };
  const inspect = () => page.evaluate(() => window.__RESOURCE_RECOVERY__);
  const agree = async name => {
    const state = await snapshot(page);
    expect(state.markedMounted).toBe(true);
    expect(state.result).toBe(state.payload);
    const expected = await semantics(page, state.result);
    expect(await semantics(page, state.mounted)).toEqual(expected);
    const exported = await download(page, 'SVG', testInfo.outputPath(`${name}.svg`));
    expect(await semantics(page, exported.toString())).toEqual(expected);
    return { request: state.request, semantics: expected };
  };
  try {
    const files = await Promise.all([fs.readFile(gffPath), fs.readFile(fastaPath)]);
    await fs.writeFile(testInfo.outputPath('inputs.json'), JSON.stringify(
      [gffPath, fastaPath, 'mismatch.fasta'].map((name, i) => ({ name, size: [...files, mismatch][i].length,
        sha256: createHash('sha256').update([...files, mismatch][i]).digest('hex') })), null, 2));
    await openApp(page);
    expect((await inspect()).workers).toEqual([]);
    await page.getByRole('button', { name: 'Linear', exact: true }).click();
    await page.getByRole('radio', { name: 'GFF3 + FASTA', exact: true }).check();
    await page.getByLabel('GFF3', { exact: true }).setInputFiles(gffPath);
    const fasta = page.getByLabel('FASTA', { exact: true });
    await fasta.setInputFiles(fastaPath);
    await generate('ok');
    const a = await agree('original-A');
    await fasta.setInputFiles({ name: 'mismatch.fasta', mimeType: 'text/plain', buffer: mismatch });
    await generate('error');
    expect(outcomes.at(-1).error.summary).toContain('No matching FASTA record');
    expect(await agree('rejected-B')).toEqual(a);
    const rejectedEvidence = await inspect();
    expect(rejectedEvidence.workers).toHaveLength(1);
    expect(rejectedEvidence.workers[0].terminated).toBe(false);
    expect(rejectedEvidence.workers[0].messages.find(m => m.requestId === rejectedEvidence.runs[1].requestId && m.type === 'run'))
      .toMatchObject({ ok: true, error: { type: 'ParseError' } });
    await fasta.setInputFiles(fastaPath);
    await generate('ok'); // Mandatory first recovery: no retry in the test.
    expect(await agree('first-recovery-A')).toEqual(a);
    const recovered = await inspect();
    expect(recovered.runs).toHaveLength(3);
    expect(recovered.runs[2].request).toEqual(a.request);
    expect(recovered.runs.map(run => run.staged.map(r => r.resourceId))).toEqual([
      ['record-1-gff3', 'record-1-fasta'], ['record-1-fasta'], ['record-1-fasta']
    ]);
    expect(recovered.runs.map(run => run.staged.reduce((sum, r) => sum + r.bytes, 0))).toEqual([86047, 424, 49253]);
    await generate('ok');
    expect((await inspect()).runs.at(-1).staged).toEqual([]);
    expect(await agree('warm-A')).toEqual(a);

    // A successful Worker render changes the same GFF resource, then Result
    // admission rejects it. The following A must restore its resource bytes.
    await page.getByLabel('GFF3', { exact: true }).setInputFiles({
      name: 'candidate.gff3', mimeType: 'text/plain', buffer: Buffer.concat([files[0], Buffer.from('\n# candidate\n')])
    });
    await page.evaluate(() => { window.__RESOURCE_RECOVERY__.rejectAdmission = true; });
    await generate('error');
    expect(outcomes.at(-1).error.summary).toContain('Injected Result admission rejection');
    expect(await agree('unadmitted-candidate')).toEqual(a);
    expect((await inspect()).workers).toHaveLength(1);
    await page.getByLabel('GFF3', { exact: true }).setInputFiles(gffPath);
    await generate('ok');
    expect((await inspect()).runs.at(-1).staged.map(r => r.resourceId)).toEqual(['record-1-gff3']);
    expect(await agree('admission-recovery-A')).toEqual(a);

    await page.getByLabel('GFF3', { exact: true }).setInputFiles({
      name: 'canceled.gff3', mimeType: 'text/plain', buffer: Buffer.concat([files[0], Buffer.from('\n# canceled candidate\n')])
    });
    await page.evaluate(() => { window.__RESOURCE_RECOVERY__.cancelAfterStaging = true; });
    await generate('canceled');
    expect(await agree('canceled')).toEqual(a);
    expect((await inspect()).workers[0].terminated).toBe(true);
    await page.getByLabel('GFF3', { exact: true }).setInputFiles(gffPath);
    await generate('ok');
    expect((await inspect()).runs.at(-1).staged.map(r => r.resourceId)).toEqual(['record-1-gff3', 'record-1-fasta']);
    expect(await agree('cancellation-recovery-A')).toEqual(a);
    const evidence = await inspect();
    expect(evidence.workers).toHaveLength(2);
    expect(evidence.workers.map(w => w.initializations)).toEqual([1, 1]);
    expect(evidence.workers[1].terminated).toBe(false);
    expect(JSON.stringify(evidence.workers)).not.toContain('Render resource cache miss');
    expect(evidence.rejections).toEqual([]);
    expect(errors).toEqual([]);
    expect(external).toEqual([]);
  } finally {
    await fs.writeFile(testInfo.outputPath('resource-recovery.json'), JSON.stringify({
      outcomes, errors, external, evidence: await inspect()
    }, null, 2));
    await context.close();
  }
});
