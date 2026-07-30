const { test, expect } = require('@playwright/test');
const { createReadStream, existsSync } = require('node:fs');
const { createServer } = require('node:http');
const { extname, join, normalize, resolve, sep } = require('node:path');

const repoRoot = resolve(process.env.GBDRAW_REPO || process.cwd());
const contentTypes = {
  '.html': 'text/html; charset=utf-8',
  '.js': 'text/javascript; charset=utf-8',
  '.mjs': 'text/javascript; charset=utf-8',
  '.css': 'text/css; charset=utf-8',
  '.json': 'application/json; charset=utf-8',
  '.svg': 'image/svg+xml',
  '.wasm': 'application/wasm',
  '.whl': 'application/octet-stream'
};

let server;
let baseUrl;

test.beforeAll(async () => {
  await new Promise((resolveServer, rejectServer) => {
    server = createServer((request, response) => {
      const url = new URL(request.url || '/', 'http://127.0.0.1');
      const requestedPath = normalize(decodeURIComponent(url.pathname)).replace(/^(\.\.(?:\/|\\|$))+/, '');
      const filePath = resolve(repoRoot, requestedPath.replace(/^[/\\]+/, ''));
      if ((!filePath.startsWith(`${repoRoot}${sep}`) && filePath !== repoRoot) || !existsSync(filePath)) {
        response.writeHead(404);
        response.end('Not found');
        return;
      }
      response.writeHead(200, {
        'Content-Type': contentTypes[extname(filePath)] || 'application/octet-stream'
      });
      createReadStream(filePath).pipe(response);
    });
    server.once('error', rejectServer);
    server.listen(0, '127.0.0.1', () => {
      baseUrl = `http://127.0.0.1:${server.address().port}`;
      resolveServer();
    });
  });
});

test.afterAll(async () => {
  await new Promise((resolveClose) => server.close(resolveClose));
});

test('Linear record rows and N-to-M comparison batches remain keyed by sequence uid', async ({ page }) => {
  await page.goto(`${baseUrl}/gbdraw/web/index.html`, { waitUntil: 'domcontentloaded' });
  await page.waitForFunction(() => window.__GBDRAW_APP__);

  const state = await page.evaluate(() => {
    const app = window.__GBDRAW_APP__;
    app.mode = 'linear';
    app.addLinearSeq();
    app.addLinearSeq();
    app.addLinearSeq();
    app.linearRecordLayoutEnabled = true;
    app.syncLinearRecordLayout();
    app.setLinearRecordRow(app.linearSeqs[0].uid, 1);
    app.setLinearRecordRow(app.linearSeqs[1].uid, 1);
    app.setLinearRecordRow(app.linearSeqs[2].uid, 2);
    app.setLinearRecordRow(app.linearSeqs[3].uid, 2);
    app.addLinearComparisonBatch(true);
    return {
      tokens: app.linearLayoutTokens,
      comparisonCount: app.linearComparisons.length,
      endpoints: app.linearComparisons.map((item) => [item.queryUid, item.subjectUid]),
      uniqueUids: new Set(app.linearSeqs.map((item) => item.uid)).size
    };
  });

  expect(state.tokens).toEqual(['#1@1', '#2@1', '#3@2', '#4@2']);
  expect(state.comparisonCount).toBe(4);
  expect(state.uniqueUids).toBe(4);
  expect(new Set(state.endpoints.map((pair) => pair.join('->'))).size).toBe(4);
  await expect(page.locator('input[aria-label="Linear record row"]')).toHaveCount(4);
  await expect(page.getByText('All adjacent-row pairs', { exact: true })).toBeVisible();
});

test('Label-scoped feature colors and legends survive linear regeneration', async ({ page }) => {
  test.setTimeout(240000);
  const makeGenbank = (recordId) => {
    const sequence = 'atg'.repeat(100);
    const origin = sequence.match(/.{1,60}/g).map((chunk, index) => {
      const groups = chunk.match(/.{1,10}/g).join(' ');
      return `${String(index * 60 + 1).padStart(9)} ${groups}`;
    }).join('\n');
    return `LOCUS       ${recordId.padEnd(24)} 300 bp    DNA     linear   UNA 01-JAN-2000
DEFINITION  feature color regeneration test.
ACCESSION   ${recordId}
VERSION     ${recordId}
KEYWORDS    .
SOURCE      synthetic construct
  ORGANISM  synthetic construct
            .
FEATURES             Location/Qualifiers
     CDS             1..90
                     /product="wsv360-like protein"
ORIGIN
${origin}
//
`;
  };

  await page.goto(`${baseUrl}/gbdraw/web/index.html`, { waitUntil: 'domcontentloaded' });
  await page.waitForFunction(() => window.__GBDRAW_APP__);
  await page.waitForFunction(() => window.__GBDRAW_APP__?.pyodideReady === true, null, { timeout: 180000 });

  await page.evaluate(({ firstRecord, secondRecord }) => {
    const app = window.__GBDRAW_APP__;
    app.mode = 'linear';
    app.lInputType = 'gb';
    app.addLinearSeq();
    app.setLinearSeqPrimaryFile(0, 'gb', new File([firstRecord], 'record-a.gbk', {
      type: 'text/plain', lastModified: 1
    }));
    app.setLinearSeqPrimaryFile(1, 'gb', new File([secondRecord], 'record-b.gbk', {
      type: 'text/plain', lastModified: 2
    }));
    Object.assign(app.form, {
      legend: 'bottom',
      show_gc: false,
      show_skew: false,
      show_depth: false,
      show_labels_linear: 'none'
    });
  }, {
    firstRecord: makeGenbank('ColorRecA'),
    secondRecord: makeGenbank('ColorRecB')
  });

  const firstRun = await page.evaluate(async () => {
    const app = window.__GBDRAW_APP__;
    const result = await app.runAnalysis();
    return { result, errorLog: app.errorLog };
  });
  expect(firstRun).toEqual({ result: { status: 'ok' }, errorLog: null });

  const edited = await page.evaluate(async () => {
    const app = window.__GBDRAW_APP__;
    const target = app.extractedFeatures.find((feature) => feature.product === 'wsv360-like protein');
    if (!target) throw new Error('Expected product feature was not extracted.');
    await app.requestFeatureColorChange(target, '#8cf04f');
    const dialog = {
      displayLabel: app.colorScopeDialog.displayLabel,
      displayLabelSiblingCount: app.colorScopeDialog.displayLabelSiblingCount
    };
    await app.handleColorScopeChoice('displayLabel');
    return {
      dialog,
      rules: app.manualSpecificRules.map(({ feat, qual, val, color, cap }) => ({
        feat, qual, val, color, cap
      }))
    };
  });
  expect(edited.dialog).toEqual({
    displayLabel: 'wsv360-like protein',
    displayLabelSiblingCount: 1
  });
  expect(edited.rules).toEqual([{
    feat: 'CDS',
    qual: 'product',
    val: '^wsv360-like protein$',
    color: '#8cf04f',
    cap: 'wsv360-like protein'
  }]);

  const secondRun = await page.evaluate(async () => {
    const app = window.__GBDRAW_APP__;
    const result = await app.runAnalysis();
    return { result, errorLog: app.errorLog };
  });
  expect(secondRun).toEqual({ result: { status: 'ok' }, errorLog: null });

  const regenerated = await page.evaluate(() => {
    const app = window.__GBDRAW_APP__;
    const svgText = app.results?.[0]?.content || '';
    const svg = new DOMParser().parseFromString(svgText, 'image/svg+xml').documentElement;
    const featureIds = app.extractedFeatures
      .filter((feature) => feature.product === 'wsv360-like protein')
      .map((feature) => feature.svg_id);
    const featureFills = featureIds.map((featureId) => {
      const roots = [...svg.querySelectorAll('[data-gbdraw-feature-id]')]
        .filter((element) => element.getAttribute('data-gbdraw-feature-id') === featureId);
      return [...roots, ...roots.flatMap((root) => [...root.querySelectorAll('[fill]')])]
        .map((element) => String(element.getAttribute('fill') || '').toLowerCase())
        .filter((fill) => fill && fill !== 'none');
    });
    const legendEntries = [...svg.querySelectorAll('[data-legend-key]')]
      .filter((entry) => entry.getAttribute('data-legend-key') === 'wsv360-like protein');
    const legendFills = legendEntries.flatMap((entry) => [...entry.querySelectorAll('[fill]')])
      .map((element) => String(element.getAttribute('fill') || '').toLowerCase())
      .filter((fill) => fill && fill !== 'none');
    return {
      featureIds,
      featureFills,
      legendEntryCount: legendEntries.length,
      legendFills,
      rules: app.manualSpecificRules.map(({ feat, qual, val, color, cap }) => ({
        feat, qual, val, color, cap
      }))
    };
  });

  expect(regenerated.featureIds).toHaveLength(2);
  regenerated.featureFills.forEach((fills) => expect(fills).toContain('#8cf04f'));
  expect(regenerated.legendEntryCount).toBeGreaterThan(0);
  expect(regenerated.legendFills).toContain('#8cf04f');
  expect(regenerated.rules).toEqual(edited.rules);

  const reset = await page.evaluate(async () => {
    const app = window.__GBDRAW_APP__;
    const target = app.extractedFeatures.find((feature) => feature.product === 'wsv360-like protein');
    if (!target) throw new Error('Expected a feature to reset.');
    const defaultColor = app.appliedPaletteColors.CDS;
    app.clickedFeature = {
      svg_id: target.svg_id,
      feat: target,
      color: '#8cf04f'
    };
    app.resetColorDialog.defaultColor = defaultColor;
    app.resetColorDialog.caption = 'wsv360-like protein';
    await app.handleResetColorChoice('this');
    return {
      defaultColor,
      targetId: target.svg_id,
      rules: app.manualSpecificRules.map(({ feat, qual, val, color, cap }) => ({
        feat, qual, val, color, cap
      }))
    };
  });
  expect(reset.rules).toContainEqual({
    feat: 'CDS',
    qual: 'hash',
    val: reset.targetId.replace(/_record_\d+$/i, ''),
    color: reset.defaultColor,
    cap: 'other proteins'
  });
  expect(reset.rules).toContainEqual(edited.rules[0]);

  const thirdRun = await page.evaluate(async () => {
    const app = window.__GBDRAW_APP__;
    const result = await app.runAnalysis();
    const svg = new DOMParser().parseFromString(app.results?.[0]?.content || '', 'image/svg+xml').documentElement;
    const fillsById = Object.fromEntries(
      app.extractedFeatures
        .filter((feature) => feature.product === 'wsv360-like protein')
        .map((feature) => {
          const roots = [...svg.querySelectorAll('[data-gbdraw-feature-id]')]
            .filter((element) => element.getAttribute('data-gbdraw-feature-id') === feature.svg_id);
          const fills = [...roots, ...roots.flatMap((root) => [...root.querySelectorAll('[fill]')])]
            .map((element) => String(element.getAttribute('fill') || '').toLowerCase())
            .filter((fill) => fill && fill !== 'none');
          return [feature.svg_id, fills];
        })
    );
    return { result, errorLog: app.errorLog, fillsById };
  });
  expect(thirdRun.result).toEqual({ status: 'ok' });
  expect(thirdRun.errorLog).toBeNull();
  expect(thirdRun.fillsById[reset.targetId]).toContain(reset.defaultColor.toLowerCase());
  const siblingFills = Object.entries(thirdRun.fillsById)
    .filter(([featureId]) => featureId !== reset.targetId)
    .flatMap(([, fills]) => fills);
  expect(siblingFills).toContain('#8cf04f');
});

test('Region annotations expose and persist an explicit target-record selection', async ({ page }) => {
  await page.goto(`${baseUrl}/gbdraw/web/index.html`, { waitUntil: 'domcontentloaded' });
  await page.waitForFunction(() => window.__GBDRAW_APP__);

  const genbank = `LOCUS       RecA                      10 bp    DNA     linear   UNA 01-JAN-2000
DEFINITION  first.
ACCESSION   RecA
VERSION     RecA
KEYWORDS    .
SOURCE      synthetic construct
  ORGANISM  synthetic construct
            .
FEATURES             Location/Qualifiers
ORIGIN
        1 aaaaaaaaaa
//
LOCUS       RecB                      12 bp    DNA     linear   UNA 01-JAN-2000
DEFINITION  second.
ACCESSION   RecB
VERSION     RecB
KEYWORDS    .
SOURCE      synthetic construct
  ORGANISM  synthetic construct
            .
FEATURES             Location/Qualifiers
ORIGIN
        1 cccccccccccc
//
`;
  await page.evaluate((content) => {
    const app = window.__GBDRAW_APP__;
    app.mode = 'linear';
    app.lInputType = 'gb';
    app.setLinearSeqPrimaryFile(0, 'gb', new File([content], 'two-records.gb', {
      type: 'text/plain',
      lastModified: 1
    }));
    const set = app.addAnnotationSet('review');
    app.addCoordinateAnnotation(set, { start: 1, end: 10 });
  }, genbank);

  await page.getByText('Region Annotations', { exact: false }).click();
  const selector = page.getByLabel('Annotation target record');
  await expect(selector).toHaveCount(1);
  await expect(selector).toHaveValue('');
  await expect(selector.locator('option')).toHaveText([
    'Select target record',
    '#1 · RecA · 10 bp',
    '#2 · RecB · 12 bp'
  ], { timeout: 60000 });
  await expect(page.getByText('Choose the record that this annotation targets.')).toBeVisible();

  const rejected = await page.evaluate(() => window.__GBDRAW_APP__.runAnalysis());
  expect(rejected).toEqual({ status: 'error' });
  await expect(page.getByText('Choose a target record for region annotation review/region_1.')).toBeVisible();

  await selector.selectOption({ label: '#2 · RecB · 12 bp' });
  await expect(page.getByText('Choose the record that this annotation targets.')).toHaveCount(0);
  const target = await page.evaluate(() => (
    window.__GBDRAW_APP__.annotationSets[0].annotations[0].target.record
  ));
  expect(target).toEqual({ kind: 'recordId', value: 'RecB' });
});

test('Region annotation IDs accept continuous typing without losing focus', async ({ page }) => {
  await page.goto(`${baseUrl}/gbdraw/web/index.html`, { waitUntil: 'domcontentloaded' });
  await page.waitForFunction(() => window.__GBDRAW_APP__);

  await page.evaluate(() => {
    const app = window.__GBDRAW_APP__;
    const set = app.addAnnotationSet('review');
    app.addCoordinateAnnotation(set);
  });

  await page.getByText('Region Annotations', { exact: false }).click();
  const idInput = page.getByPlaceholder('annotation_id');
  await idInput.selectText();
  await idInput.pressSequentially('Repeat');

  await expect(idInput).toBeFocused();
  await expect(idInput).toHaveValue('Repeat');
  await idInput.press('Tab');
  await expect.poll(() => page.evaluate(() => (
    window.__GBDRAW_APP__.annotationSets[0].annotations[0].id
  ))).toBe('Repeat');
});

test('GFF annotation targets follow FASTA record order', async ({ page }) => {
  await page.goto(`${baseUrl}/gbdraw/web/index.html`, { waitUntil: 'domcontentloaded' });
  await page.waitForFunction(() => window.__GBDRAW_APP__);

  const gff = `##gff-version 3
##sequence-region RecB 1 12
RecB\ttest\tgene\t1\t3\t.\t+\t.\tID=gene_b
##sequence-region RecA 1 10
RecA\ttest\tgene\t2\t4\t.\t+\t.\tID=gene_a
`;
  const fasta = `>RecB
CCCCCCCCCCCC
>RecA
AAAAAAAAAA
`;
  await page.evaluate(({ gffText, fastaText }) => {
    const app = window.__GBDRAW_APP__;
    app.mode = 'linear';
    app.lInputType = 'gff';
    app.setLinearSeqPrimaryFile(0, 'gff', new File([gffText], 'records.gff3', {
      type: 'text/plain',
      lastModified: 2
    }));
    app.setLinearSeqPrimaryFile(0, 'fasta', new File([fastaText], 'records.fasta', {
      type: 'text/plain',
      lastModified: 3
    }));
    const set = app.addAnnotationSet('gff-review');
    app.addCoordinateAnnotation(set, { start: 1, end: 3 });
  }, { gffText: gff, fastaText: fasta });

  await page.getByText('Region Annotations', { exact: false }).click();
  await expect(page.getByLabel('Annotation target record').locator('option')).toHaveText([
    'Select target record',
    '#1 · RecB · 12 bp',
    '#2 · RecA · 10 bp'
  ], { timeout: 60000 });
});
