const { test, expect } = require('@playwright/test');
const { spawnSync } = require('node:child_process');
const { readFileSync, writeFileSync } = require('node:fs');
const { basename, resolve } = require('node:path');
const { openApp } = require('../helpers/app-lifecycle.cjs');

const repoRoot = resolve(process.env.GBDRAW_REPO || process.cwd());
const examples = JSON.parse(readFileSync(
  resolve(repoRoot, 'gbdraw/web/gallery/examples.json'),
  'utf8'
));
const isolatedExamples = new Set([
  'vibrio-harveyi-group-collinear'
]);
const commonExamples = examples.filter(({ id }) => !isolatedExamples.has(id));
const compareCommand = [
  'import sys',
  'from tests.utils.svg_compare import compare_svgs',
  'result = compare_svgs(sys.argv[1], sys.argv[2])',
  'print(result.message)',
  'print("\\n".join(result.differences))',
  'raise SystemExit(0 if result.equal else 1)'
].join(';');

test.describe.configure({ mode: 'serial' });

for (const example of commonExamples) {
  test(`${example.id} preserves committed SVG semantics on first Generate`, async ({
    browser
  }, testInfo) => {
    test.setTimeout(600_000);
    const context = await browser.newContext();
    const page = await context.newPage();
    const externalRequests = [];
    const browserErrors = [];
    page.on('console', (message) => {
      if (message.type() === 'error') browserErrors.push(message.text());
    });
    page.on('pageerror', (error) => browserErrors.push(error.message));
    page.on('request', (request) => {
      const hostname = new URL(request.url()).hostname;
      if (!['127.0.0.1', 'localhost'].includes(hostname)) {
        externalRequests.push(request.url());
      }
    });
    page.on('dialog', (dialog) => dialog.accept());
    await page.addInitScript(() => {
      window.__GBDRAW_GALLERY_LOSAT_CALLS__ = 0;
      window.__GBDRAW_LOSAT_EXECUTOR__ = async () => {
        window.__GBDRAW_GALLERY_LOSAT_CALLS__ += 1;
        throw new Error('Published Gallery replay must use its committed cache.');
      };
    });

    try {
      await openApp(page, { waitForPalette: false });
      const loaded = await page.evaluate(async ({ session, id }) => {
        const response = await fetch(`/gbdraw/web/gallery/${session.replace(/^\.\//, '')}`);
        if (!response.ok) throw new Error(`Could not load ${id}: ${response.status}`);
        const bytes = await response.arrayBuffer();
        const file = new File(
          [bytes],
          session.split('/').pop(),
          { type: session.endsWith('.gz') ? 'application/gzip' : 'application/json' }
        );
        const result = await window.__GBDRAW_APP__.importSession({
          target: { files: [file], value: '' }
        });
        const app = window.__GBDRAW_APP__;
        const { state } = await import('/gbdraw/web/js/state.js');
        const selected = app.results?.[app.selectedResultIndex];
        return {
          status: result?.status,
          svg: String(selected?.content || ''),
          mode: state.mode.value,
          recordCount: state.mode.value === 'linear'
            ? state.linearSeqs.length
            : state.circularRecordList.value.length,
          slotIds: state.adv.circular_track_slots.map(({ id }) => id),
          annotationSetIds: state.annotationSets.map(({ id: setId }) => setId),
          legend: state.form.legend,
          arrowShaftWidthRatio: state.adv.arrow_shaft_width_ratio,
          ruleCaptions: state.manualSpecificRules.map(({ cap }) => cap).filter(Boolean)
        };
      }, example);
      expect(loaded.status).toBe('ok');
      expect(loaded.svg).toContain('<svg');

      if (example.id === 'HmmtDNA_ATskew') {
        expect(loaded.slotIds).toEqual([
          'features', 'gc_content', 'gc_skew', 'a_skew_2', 'ticks'
        ]);
      } else if (example.id === 'tobacco-chloroplast') {
        expect(loaded.annotationSetIds).toContain('plastome_regions');
        expect(loaded.slotIds).toHaveLength(3);
        expect(loaded.legend).toBe('upper_left');
      } else if (example.id === 'majanivirus_orthogroup') {
        expect(loaded.recordCount).toBe(9);
        expect(Number(loaded.arrowShaftWidthRatio)).toBe(1);
        expect(loaded.ruleCaptions).toHaveLength(3);
      }

      const generated = await page.evaluate(async () => {
        const app = window.__GBDRAW_APP__;
        const result = await app.runAnalysis();
        const { state } = await import('/gbdraw/web/js/state.js');
        const error = state.errorLog.value;
        const selected = app.results?.[app.selectedResultIndex];
        const svg = String(selected?.content || '');
        return {
          status: result?.status,
          error: result?.error?.message || result?.message || error?.summary
            || error?.message || (error ? JSON.stringify(error) : null),
          svg,
          matchCount: (svg.match(/data-gbdraw-pairwise-match-id=/g) || []).length,
          losatCalls: window.__GBDRAW_GALLERY_LOSAT_CALLS__
        };
      });
      expect(
        generated.status,
        [generated.error, ...browserErrors]
          .filter(Boolean).join('\n')
      ).toBe('ok');
      expect(generated.error).toBeNull();
      expect(generated.losatCalls).toBe(0);
      expect(externalRequests).toEqual([]);
      if (example.id === 'majanivirus_orthogroup') {
        expect(generated.matchCount).toBe(627);
      }

      const loadedPath = testInfo.outputPath(`${basename(example.id)}-loaded.svg`);
      const generatedPath = testInfo.outputPath(`${basename(example.id)}-generated.svg`);
      writeFileSync(loadedPath, loaded.svg, 'utf8');
      writeFileSync(generatedPath, generated.svg, 'utf8');
      const comparison = spawnSync(
        process.env.GBDRAW_PYTHON || 'python',
        ['-c', compareCommand, loadedPath, generatedPath],
        { cwd: repoRoot, encoding: 'utf8' }
      );
      expect(
        comparison.status,
        `${comparison.stdout}\n${comparison.stderr}`
      ).toBe(0);
    } finally {
      await context.close();
    }
  });
}
