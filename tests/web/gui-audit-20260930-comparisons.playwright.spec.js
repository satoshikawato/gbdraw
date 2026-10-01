// Web GUI audit 2026-09-30: Linear comparisons and CLI comparison Sessions.
// A test marked test.fail(true, '<ID>') asserts the correct behavior of a
// current defect; the PR that fixes the audit ID removes the mark.
const { test, expect } = require('@playwright/test');
const { execFile } = require('node:child_process');
const { mkdirSync, readFileSync, writeFileSync } = require('node:fs');
const path = require('node:path');
const { promisify } = require('node:util');
const { evaluateWithRetainedPromise, generateAndWaitForResult } = require('./helpers/app-lifecycle.cjs');
const { BATCH_FIXTURE, openFresh, settle } = require('./helpers/audit-browser.cjs');

test.describe.configure({ retries: 0 });

const root = process.cwd();
const splitRecords = (file) => readFileSync(file, 'utf8').split(/^\/\/\s*$/m)
  .map((chunk) => chunk.trim()).filter((chunk) => chunk.startsWith('LOCUS')).map((chunk) => `${chunk}\n//\n`);
// R2c 2001..3000 and R3c 1..1000 are one shared 1 kb block; nothing else is shared.
const [R2C, R3C] = splitRecords('tests/fixtures/web_comparison_shared_block.gb');
const [TESTA, TESTB] = splitRecords(BATCH_FIXTURE);
const bgc = (id) => ({ name: `${id}.gbk`, text: readFileSync(path.join(root, 'tests/test_inputs', `${id}.gbk`), 'utf8') });

const openLinearWith = async (page, files) => {
  await openFresh(page);
  await page.getByRole('button', { name: 'Linear', exact: true }).click();
  await page.waitForFunction(() => window.__GBDRAW_APP__?.mode === 'linear');
  await page.evaluate(async (items) => {
    const app = window.__GBDRAW_APP__;
    while (app.linearSeqs.length < items.length) app.addLinearSeq();
    items.forEach((item, index) => app.setLinearSeqPrimaryFile(index, 'gb', new File(
      [item.text], item.name, { type: 'text/plain', lastModified: 1000 + index }
    )));
    await window.Vue.nextTick();
  }, files);
  await settle(page);
};

// Serial LOSAT: the Playwright server does not send COOP/COEP (CO-01).
const useLosat = async (page, { task, proteinMode = null, configure = null }) => {
  await page.evaluate(async (options) => {
    const app = window.__GBDRAW_APP__;
    await app.setLinearComparisonGlobalAction('losat');
    app.setLinearComparisonLosatMode(options.task);
    if (options.proteinMode) app.setLinearComparisonLosatpMode(options.proteinMode);
    app.losat.executionMode = 'serial';
    app.adv.pairwise_match_style = 'ribbon';
  }, { task, proteinMode });
  if (configure) await page.evaluate(configure);
  await settle(page);
};

const uploadComparison = async (page, name, text) => {
  await page.evaluate(async ({ fileName, table }) => {
    const app = window.__GBDRAW_APP__;
    await app.setLinearComparisonGlobalAction('upload');
    const edgeKey = app.linearComparisonTimeline.rows[0].boundaryAfter.pairs[0].edgeKey;
    app.setLinearComparisonCardFile(edgeKey, new File([table], fileName, { type: 'text/tab-separated-values' }));
    await window.Vue.nextTick();
  }, { fileName: name, table: text });
  await settle(page);
};

// Committed ribbons: data attributes and the x span of each endpoint side.
// A ribbon path runs query side (points 1-2) then subject side (points 3-4).
const committedRibbons = (page) => page.evaluate(() => {
  const app = window.__GBDRAW_APP__;
  const content = String(app.results?.[app.selectedResultIndex]?.content || '');
  const svg = new DOMParser().parseFromString(content, 'image/svg+xml');
  const numbers = (value) => (String(value || '').match(/-?\d+(?:\.\d+)?(?:e-?\d+)?/gi) || []).map(Number);
  const span = (values) => [Math.round(Math.min(...values)), Math.round(Math.max(...values))];
  return [...svg.querySelectorAll('[data-gbdraw-pairwise-match-id]')].map((element) => {
    const point = numbers(element.getAttribute('d'));
    const data = Object.fromEntries([...element.attributes]
      .filter(({ name }) => name.startsWith('data-'))
      .map(({ name, value }) => [name.slice(5), value]));
    return { data, queryX: span([point[0], point[2]]), subjectX: span([point[4], point[6]]) };
  });
});

const featureSpan = (page, locusTag) => page.evaluate((tag) => {
  const app = window.__GBDRAW_APP__;
  const feature = app.extractedFeatures.find((item) => item.locus_tag === tag);
  const content = String(app.results?.[app.selectedResultIndex]?.content || '');
  const svg = new DOMParser().parseFromString(content, 'image/svg+xml');
  const xs = [...svg.querySelectorAll(`[data-gbdraw-feature-id="${CSS.escape(feature.svg_id)}"]`)]
    .flatMap((element) => (String(element.getAttribute('d') || '').match(/-?\d+(?:\.\d+)?/g) || [])
      .map(Number).filter((_, index) => index % 2 === 0));
  return [Math.round(Math.min(...xs)), Math.round(Math.max(...xs))];
}, locusTag);

const contains = (outer, inner) => outer[0] <= inner[0] + 1 && outer[1] >= inner[1] - 1;

test('a LOSATP rerun after a Feature visibility change matches a fresh run with that rule', async ({ page, browser }) => {
  test.fail(true, 'CO-02');
  test.setTimeout(900_000);
  const files = [bgc('BGC0000708'), bgc('BGC0000709')];
  const hideProtein = async () => {
    const app = window.__GBDRAW_APP__;
    await app.addFeatureVisibilityRule();
    const index = app.featureVisibilityManualRules.length - 1;
    for (const [field, value] of Object.entries({
      recordId: '*', featureType: 'CDS', qualifier: 'protein_id', value: 'CAG38712.1', action: 'off'
    })) await app.setFeatureVisibilityRuleField(index, field, value);
  };
  const summary = async (target) => {
    const ribbons = await committedRibbons(target);
    return {
      ribbons: ribbons.length,
      groups: new Set(ribbons.map(({ data }) => data['orthogroup-id'])).size,
      unboundEndpoints: ribbons.filter(({ data }) => !data['query-feature-svg-id'] || !data['subject-feature-svg-id']).length
    };
  };
  await openLinearWith(page, files);
  await useLosat(page, { task: 'blastp', proteinMode: 'orthogroup' });
  await generateAndWaitForResult(page);
  await page.evaluate(hideProtein);
  await settle(page);
  expect(await page.evaluate(() => window.__GBDRAW_APP__.featureVisibilityManualRules.map(({ value }) => value)))
    .toEqual(['CAG38712.1']);
  await generateAndWaitForResult(page);
  const rerun = await summary(page);

  const freshPage = await browser.newPage();
  try {
    await openLinearWith(freshPage, files);
    await useLosat(freshPage, { task: 'blastp', proteinMode: 'orthogroup', configure: hideProtein });
    await generateAndWaitForResult(freshPage);
    const fresh = await summary(freshPage);
    expect(fresh.unboundEndpoints).toBe(0);
    expect(rerun).toEqual(fresh);
  } finally {
    await freshPage.close();
  }
});

test('Rotate with Orient feature forward keeps LOSATN ribbons on the homologous block', async ({ page }) => {
  test.fail(true, 'CO-03');
  test.setTimeout(600_000);
  await openLinearWith(page, [{ name: 'R2c.gb', text: R2C }, { name: 'R3c.gb', text: R3C }]);
  await useLosat(page, { task: 'blastn' });
  await generateAndWaitForResult(page);
  // R3_1 (101..400) lies inside the block shared with R2c (R3c 1..1000).
  const baseline = await committedRibbons(page);
  expect(baseline).toHaveLength(1);
  expect(contains(baseline[0].subjectX, await featureSpan(page, 'R3_1'))).toBe(true);

  const minusCds = await page.evaluate(() => window.__GBDRAW_APP__.extractedFeatures
    .find((feature) => feature.locus_tag === 'R3_2').svg_id);
  await page.locator(`svg [data-gbdraw-feature-id="${minusCds}"]`).first()
    .dispatchEvent('click', { clientX: 200, clientY: 200 });
  await page.waitForFunction(() => Boolean(window.__GBDRAW_APP__.clickedFeature?.feat), null, { timeout: 30_000 });
  await page.evaluate(() => window.__GBDRAW_APP__.setFeatureRecordRotationOrientForward(true));
  expect(await page.evaluate(() => window.__GBDRAW_APP__.featureRecordRotationDraft.canApply)).toBe(true);
  await evaluateWithRetainedPromise(page, async () => {
    await window.__GBDRAW_APP__.applyFeatureRecordRotation();
  });
  await settle(page);
  await generateAndWaitForResult(page);
  const rotated = await committedRibbons(page);
  expect(rotated).toHaveLength(1);
  const shared = await featureSpan(page, 'R3_1');
  expect(contains(rotated[0].subjectX, shared), JSON.stringify({ ribbon: rotated[0].subjectX, shared })).toBe(true);
});

test('Save Raw LOSAT TSV of a reversed record re-uploads to the same ribbons', async ({ page }) => {
  test.fail(true, 'CO-07');
  test.setTimeout(600_000);
  await openLinearWith(page, [{ name: 'R2c.gb', text: R2C }, { name: 'R3c.gb', text: R3C }]);
  await useLosat(page, {
    task: 'blastn',
    configure: () => { window.__GBDRAW_APP__.linearSeqs[1].region_reverse = true; }
  });
  await generateAndWaitForResult(page);
  const searched = (await committedRibbons(page)).map(({ queryX, subjectX }) => ({ queryX, subjectX }));
  expect(searched).toHaveLength(1);
  const edgeKey = await page.evaluate(() => window.__GBDRAW_APP__.linearComparisonResolution.edges[0].edgeKey);
  const pending = page.waitForEvent('download', { timeout: 60_000 });
  await page.evaluate((key) => window.__GBDRAW_APP__.downloadLosatPair(key, ''), edgeKey);
  const exported = readFileSync(await (await pending).path(), 'utf8');
  expect(exported.trim()).not.toBe('');
  await uploadComparison(page, 'R2c.R3c.losatn.tsv', exported);
  await generateAndWaitForResult(page);
  const uploaded = (await committedRibbons(page)).map(({ queryX, subjectX }) => ({ queryX, subjectX }));
  expect(uploaded).toEqual(searched);
});

test('a renamed Similarity group keeps its name only on the same members after regrouping', async ({ page }) => {
  test.fail(true, 'CO-08');
  test.setTimeout(900_000);
  await openLinearWith(page, [bgc('BGC0000708'), bgc('BGC0000709')]);
  await useLosat(page, { task: 'blastp', proteinMode: 'orthogroup' });
  await generateAndWaitForResult(page);
  const groups = () => page.evaluate(() => {
    const app = window.__GBDRAW_APP__;
    const member = (item) => String(item.sourceProteinId || item.source_protein_id || item.proteinId
      || `${item.recordIndex}:${item.start}-${item.end}`);
    return (app.orthogroups || []).map((group) => ({
      id: group.id,
      name: app.resolveOrthogroupName(group),
      members: (group.members || []).map(member).sort().join(',')
    }));
  });
  const before = await groups();
  const target = before.find(({ members }) => members.includes('CAG38711.1'));
  expect(target).toBeTruthy();
  await page.evaluate((id) => window.__GBDRAW_APP__.setOrthogroupNameOverride(id, 'LivB-renamed'), target.id);
  await settle(page);
  await page.evaluate(() => { window.__GBDRAW_APP__.adv.identity = 46; });
  await settle(page);
  await generateAndWaitForResult(page);
  const after = await groups();
  expect(after.length).not.toBe(before.length);
  expect(after.filter(({ name }) => name === 'LivB-renamed').map(({ members }) => members)
    .filter((members) => members !== target.members)).toEqual([]);
});

test('comparison filters outside their domain are rejected instead of replaced', async ({ page }) => {
  test.setTimeout(600_000);
  await openLinearWith(page, [{ name: 'TESTA.gb', text: TESTA }, { name: 'TESTB.gb', text: TESTB }]);
  await uploadComparison(page, 'TESTA.TESTB.tsv', [
    'TESTA\tTESTB\t95.0\t300\t10\t0\t301\t600\t301\t600\t1e-50\t300',
    'TESTA\tTESTB\t80.0\t100\t20\t0\t2001\t2100\t2100\t2001\t1e-20\t120'
  ].join('\n') + '\n');
  const defaults = await page.evaluate(() => {
    const { identity, min_bitscore: bitscore, alignment_length: length } = window.__GBDRAW_APP__.adv;
    return { identity, bitscore, length };
  });
  await generateAndWaitForResult(page);
  for (const [field, value] of [['identity', -5], ['min_bitscore', -3], ['alignment_length', 150.5]]) {
    await page.evaluate(({ name, invalid, initial }) => {
      Object.assign(window.__GBDRAW_APP__.adv, {
        identity: initial.identity, min_bitscore: initial.bitscore, alignment_length: initial.length, [name]: invalid
      });
    }, { name: field, invalid: value, initial: defaults });
    await settle(page);
    await generateAndWaitForResult(page, { expectedStatus: 'error' });
    expect(await page.evaluate((name) => window.__GBDRAW_APP__.adv[name], field)).toBe(value);
  }
});

// SE-06 (current CLI sidecar) and N-17 (schema-7 sidecar written by main; see
// tests/fixtures/sessions/se06-main-linear-blast-cli.provenance.json).
const cliBlastSessions = {
  'current CLI': async (directory) => {
    writeFileSync(path.join(directory, 'R2c.gb'), R2C);
    writeFileSync(path.join(directory, 'R3c.gb'), R3C);
    writeFileSync(path.join(directory, 'R2c_R3c.tsv'), 'R2c\tR3c\t100.000\t1000\t0\t0\t2001\t3000\t1\t1000\t0.0\t1847\n');
    const session = path.join(directory, 'cli-blast.gbdraw-session.json.gz');
    await promisify(execFile)('python', [
      '-m', 'gbdraw.cli', 'linear', '--gbk', 'R2c.gb', 'R3c.gb', '-b', 'R2c_R3c.tsv',
      '-o', 'cli-blast', '--session_output', session
    ], { cwd: directory, env: { ...process.env, PYTHONPATH: root }, timeout: 600_000, maxBuffer: 1_000_000 });
    return session;
  },
  'main v42 CLI': async () => path.join(root, 'tests/fixtures/sessions/se06-main-linear-blast-cli.v42.gbdraw-session.json.gz')
};
for (const [source, writeSession] of Object.entries(cliBlastSessions)) {
  test(`a ${source} Linear BLAST Session inherits its comparison and generates`, async ({ page }, testInfo) => {
    test.setTimeout(600_000);
    const directory = testInfo.outputDir;
    mkdirSync(directory, { recursive: true });
    const session = await writeSession(directory);
    await openFresh(page);
    await page.locator('input[accept^=".json,"]').setInputFiles(session);
    await page.waitForFunction(() => !window.__GBDRAW_APP__.sessionImportPending
      && window.__GBDRAW_APP__.results.length > 0, null, { timeout: 300_000 });
    await settle(page);
    const loaded = await committedRibbons(page);
    expect(loaded).toHaveLength(1);
    // D-36: the CLI comparison stays read-only and is reused through Inherit.
    expect(await page.evaluate(() => window.__GBDRAW_APP__.importedComparisonIntent.disposition))
      .toBe('PRESERVED_READ_ONLY');
    await evaluateWithRetainedPromise(page, async () => {
      await window.__GBDRAW_APP__.inheritImportedComparison();
    });
    await settle(page);
    await generateAndWaitForResult(page);
    expect((await committedRibbons(page)).map(({ queryX, subjectX }) => ({ queryX, subjectX })))
      .toEqual(loaded.map(({ queryX, subjectX }) => ({ queryX, subjectX })));
  });
}

test('the match popup and its FASTA header report source coordinates of a cropped record', async ({ page }) => {
  test.fail(true, 'CO-10');
  test.setTimeout(600_000);
  await openLinearWith(page, [{ name: 'R2c.gb', text: R2C }, { name: 'R3c.gb', text: R3C }]);
  await useLosat(page, {
    task: 'blastn',
    configure: () => {
      const app = window.__GBDRAW_APP__;
      app.setLinearRecordCrop(app.linearSeqs[0], 'region_start', 1001);
      app.setLinearRecordCrop(app.linearSeqs[0], 'region_end', 3000);
    }
  });
  await generateAndWaitForResult(page);
  expect(await committedRibbons(page)).toHaveLength(1);
  const popup = await evaluateWithRetainedPromise(page, async () => {
    const app = window.__GBDRAW_APP__;
    app.clickedPairwiseMatch = null;
    app.svgContainer.querySelector('svg [data-gbdraw-pairwise-match-id]')
      .dispatchEvent(new MouseEvent('click', { bubbles: true, clientX: 100, clientY: 100 }));
    for (let attempt = 0; attempt < 100 && !app.clickedPairwiseMatch; attempt += 1) {
      await new Promise((resolve) => setTimeout(resolve, 50));
    }
    const match = JSON.parse(JSON.stringify(app.clickedPairwiseMatch));
    return {
      query: JSON.stringify(match.sections.find(({ title }) => title === 'Query') || null),
      headers: (match.sequenceBundle?.entries || []).map(({ fasta }) => String(fasta || '').split('\n')[0])
    };
  });
  // The shared block is R2c source 2001..3000, which is 1001..2000 of the crop.
  const plain = (text) => text.replace(/(\d),(?=\d{3})/g, '$1');
  expect(plain(popup.query), popup.query).toMatch(/2001\.\.3000/);
  const queryHeader = popup.headers.find((header) => /_query\|/.test(header)) || '';
  expect(queryHeader, JSON.stringify(popup.headers)).toMatch(/coords=2001\.\.3000/);
});

