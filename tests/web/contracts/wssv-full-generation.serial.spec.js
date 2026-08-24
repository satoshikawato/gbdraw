const { test, expect } = require('@playwright/test');
const { spawnSync } = require('node:child_process');
const { writeFileSync } = require('node:fs');
const { resolve } = require('node:path');
const {
  getDiagramWorkerActivity,
  openApp
} = require('../helpers/app-lifecycle.cjs');

const SESSION_URL = '/gbdraw/web/gallery/sessions/WSSV_genome_comparison.gbdraw-session.json';
const EXPECTED_FILES = [
  'CN01.fasta',
  'WSSV-TW.fasta',
  'WSSV-CN.fasta',
  'WSSV-TH.fasta',
  'JP01A.fasta',
  'JP01B.fasta',
  'Pc2020.fasta',
  'E1.fasta',
  '0722-1.fasta',
  'CN03.fasta',
  'CN04.fasta',
  'WSSV-AU.fasta',
  'EU129.fa',
  'GCF7.fa',
  'MES-753.fa',
  'Shantou2019.fa',
  'POMZ1.fa',
  'POMZ4.fa',
  'MG18PR-0187-N40S.fa',
  'Angostura2013.fa'
];
const EXPECTED_LABELS = [
  'CN01',
  'WSSV-TW',
  'WSSV-CN',
  'WSSV-TH',
  'JP01A',
  'JP01B',
  'Pc2020',
  'E1',
  '0722-1',
  'CN03',
  'CN04',
  'WSSV-AU',
  'EU129',
  'GCF7',
  'MES-753',
  'Shantou2019',
  'POMZ1',
  'POMZ4',
  'MG18PR-0187-N40S',
  'Angostura2013'
];
const EXPECTED_COLORS = [
  '#6e91b7',
  '#f4a251',
  '#77b26f',
  '#e67577',
  '#8fc4c0',
  '#f0d369',
  '#be92b2',
  '#ffafb7',
  '#ae8e7c',
  '#c6bebb',
  '#6e91b7',
  '#f4a251',
  '#e67577',
  '#8fc4c0',
  '#bcb4ca',
  '#f0d369',
  '#be92b2',
  '#ffafb7',
  '#ae8e7c',
  '#c6bebb'
];
const repoRoot = resolve(process.env.GBDRAW_REPO || process.cwd());

const installWssvProbe = (page) => page.addInitScript(() => {
  const metrics = [];
  window.__GBDRAW_WSSV_PROBE__ = {
    metrics,
    losatCalls: 0,
    snapshot() {
      return {
        metrics: metrics.map((metric) => ({ ...metric })),
        losatCalls: this.losatCalls
      };
    }
  };
  window.__GBDRAW_TEST_HOOKS__ = {
    onStructuralMetric(metric) {
      metrics.push({ ...metric });
    }
  };
  window.__GBDRAW_LOSAT_EXECUTOR__ = async () => {
    window.__GBDRAW_WSSV_PROBE__.losatCalls += 1;
    throw new Error('WSSV replay must not execute LOSAT outside its verified cache.');
  };
});

test.describe.configure({ mode: 'serial' });

test('bundled WSSV session regenerates 20 conservation rings from lazy FASTAs', async ({
  browser
}, testInfo) => {
  const context = await browser.newContext();
  const page = await context.newPage();
  const externalRequests = [];
  page.on('request', (request) => {
    const url = new URL(request.url());
    if (!['127.0.0.1', 'localhost'].includes(url.hostname)) {
      externalRequests.push(request.url());
    }
  });
  page.on('dialog', (dialog) => dialog.accept());

  try {
    await installWssvProbe(page);
    await openApp(page);

    const imported = await page.evaluate(async (sessionUrl) => {
      const response = await fetch(sessionUrl);
      if (!response.ok) throw new Error(`Could not load WSSV session: ${response.status}`);
      const file = new File(
        [await response.text()],
        'WSSV_genome_comparison.gbdraw-session.json',
        { type: 'application/json' }
      );
      const result = await window.__GBDRAW_APP__.importSession({
        target: { files: [file], value: '' }
      });
      const { state } = await import('/gbdraw/web/js/state.js');
      const sources = state.files.c_conservation_fastas;
      const selected = window.__GBDRAW_APP__.results[
        window.__GBDRAW_APP__.selectedResultIndex
      ];
      return {
        result,
        content: String(selected?.content || ''),
        previewVisible: Boolean(document.querySelector('.shadow-xl.origin-top > svg')),
        resultCount: window.__GBDRAW_APP__.results.length,
        selectedResultIndex: window.__GBDRAW_APP__.selectedResultIndex,
        names: sources.map((source) => source.name),
        metadataOnly: sources.map((source) => ({
          frozen: Object.isFrozen(source),
          ownFields: ['text', 'arrayBuffer', 'data', 'resourceId']
            .filter((field) => Object.hasOwn(source, field))
        })),
        metrics: window.__GBDRAW_WSSV_PROBE__.snapshot().metrics
      };
    }, SESSION_URL);

    expect(imported.result?.status).toBe('ok');
    expect(imported.previewVisible).toBe(true);
    expect(imported.resultCount).toBe(1);
    expect(imported.selectedResultIndex).toBe(0);
    expect(imported.names).toEqual(EXPECTED_FILES);
    expect(imported.metadataOnly).toEqual(
      EXPECTED_FILES.map(() => ({ frozen: true, ownFields: [] }))
    );
    expect(imported.metrics.filter(({ name, resourceId }) => (
      name === 'base64DecodeCount'
      && (
        resourceId === 'record-1-genbank'
        || String(resourceId || '').startsWith('conservation-losat-fasta-files-')
      )
    ))).toHaveLength(21);
    expect(await getDiagramWorkerActivity(page)).toMatchObject({
      constructions: 0,
      initializations: 0,
      helpers: 0,
      runs: 0
    });

    const generated = await page.evaluate(async () => {
      const app = window.__GBDRAW_APP__;
      const startedAt = performance.now();
      const result = await app.runAnalysis();
      const elapsedMs = performance.now() - startedAt;
      const selected = app.results?.[app.selectedResultIndex];
      const content = String(selected?.content || '');
      const documentRoot = new DOMParser().parseFromString(content, 'image/svg+xml');
      const {
        isCommittedSvgResult,
        isCommittedSvgResultMounted
      } = await import('/gbdraw/web/js/services/svg-result-ingestion.js');
      const slots = [...documentRoot.querySelectorAll(
        '[data-gbdraw-slot-renderer="sequence_conservation"]'
      )];
      return {
        result,
        elapsedMs,
        errorLog: app.errorLog,
        resultCount: app.results.length,
        selectedResultIndex: app.selectedResultIndex,
        selectedName: String(selected?.name || ''),
        committedResult: isCommittedSvgResult(selected),
        mountedResult: isCommittedSvgResultMounted(selected),
        previewVisible: Boolean(document.querySelector('.shadow-xl.origin-top > svg')),
        slots: slots.map((slot) => ({
          id: slot.getAttribute('data-gbdraw-slot-id'),
          sourceIndex: Number(slot.getAttribute('data-source-index')),
          label: slot.getAttribute('data-track-label'),
          color: slot.getAttribute('data-track-color')
        })),
        renderedMatchCount: documentRoot.querySelectorAll(
          '[data-gbdraw-match-id]'
        ).length,
        comparisonLegendCount: documentRoot.querySelectorAll(
          '[data-gbdraw-role="comparison-legend"]'
        ).length,
        unsafeContent: /<script|foreignObject|\son\w+\s*=|javascript:/i.test(content),
        oldMissingFileMessage: content.includes(
          'Input file is not available for browser FASTA extraction.'
        ),
        content,
        probe: window.__GBDRAW_WSSV_PROBE__.snapshot()
      };
    });

    expect(generated.result).toEqual({ status: 'ok' });
    expect(generated.errorLog).toBeNull();
    expect(generated.resultCount).toBe(1);
    expect(generated.selectedResultIndex).toBe(0);
    expect(generated.selectedName).toBe('out.svg');
    expect(generated.committedResult).toBe(true);
    expect(generated.mountedResult).toBe(true);
    expect(generated.previewVisible).toBe(true);
    expect(generated.slots).toEqual(EXPECTED_LABELS.map((label, index) => ({
      id: `conservation_${index + 1}`,
      sourceIndex: index,
      label,
      color: EXPECTED_COLORS[index]
    })));
    expect(generated.renderedMatchCount).toBeGreaterThan(0);
    expect(generated.comparisonLegendCount).toBe(1);
    expect(generated.unsafeContent).toBe(false);
    expect(generated.oldMissingFileMessage).toBe(false);
    expect(generated.probe.losatCalls).toBe(0);
    expect(externalRequests).toEqual([]);

    const initialPath = testInfo.outputPath('wssv-initial.svg');
    const generatedPath = testInfo.outputPath('wssv-first-generate.svg');
    writeFileSync(initialPath, imported.content, 'utf8');
    writeFileSync(generatedPath, generated.content, 'utf8');
    const comparisonCommand = [
      'import sys',
      'from tests.utils.svg_compare import compare_svgs',
      'result = compare_svgs(sys.argv[1], sys.argv[2])',
      'print(result.message)',
      'print("\\n".join(result.differences))',
      'raise SystemExit(0 if result.equal else 1)'
    ].join(';');
    const semanticComparison = spawnSync(
      process.env.GBDRAW_PYTHON || 'python',
      ['-c', comparisonCommand, initialPath, generatedPath],
      { cwd: repoRoot, encoding: 'utf8' }
    );
    expect(
      semanticComparison.status,
      `${semanticComparison.stdout}\n${semanticComparison.stderr}`
    ).toBe(0);

    const fastaTextReads = generated.probe.metrics.filter(({ name, resourceId }) => (
      name === 'resourceTextReadCount'
      && String(resourceId || '').startsWith('conservation-losat-fasta-files-')
    ));
    expect(fastaTextReads.map(({ resourceId }) => resourceId)).toEqual(
      EXPECTED_FILES.map((_, index) => `conservation-losat-fasta-files-${index + 1}`)
    );
    expect(generated.probe.metrics.filter(({ name, resourceId }) => (
      name === 'base64DecodeCount'
      && (
        resourceId === 'record-1-genbank'
        || String(resourceId || '').startsWith('conservation-losat-fasta-files-')
      )
    ))).toHaveLength(21);

    const worker = await getDiagramWorkerActivity(page);
    expect(worker.constructions).toBe(1);
    expect(worker.initializations).toBe(1);
    expect(worker.runs).toBe(1);
    expect(worker.settledInitializations).toBe(1);
    expect(worker.settledRuns).toBe(1);
    expect(worker.instances).toHaveLength(1);
    expect(worker.instances[0].terminated).toBe(false);

    const evidence = {
      elapsedMs: generated.elapsedMs,
      resourceTextReads: fastaTextReads.length,
      ringCount: generated.slots.length,
      renderedMatchCount: generated.renderedMatchCount,
      losatCalls: generated.probe.losatCalls,
      worker
    };
    await testInfo.attach('wssv-full-generation-evidence.json', {
      body: Buffer.from(JSON.stringify(evidence, null, 2)),
      contentType: 'application/json'
    });
    console.log(`WSSV full-generation evidence: ${JSON.stringify(evidence)}`);
  } finally {
    await context.close();
  }
});
