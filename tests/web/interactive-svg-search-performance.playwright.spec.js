const { execFileSync } = require('node:child_process');
const fs = require('node:fs');
const http = require('node:http');
const path = require('node:path');
const { pathToFileURL } = require('node:url');
const { test, expect } = require('@playwright/test');

let moduleServer;
let moduleOrigin;

test.beforeAll(async () => {
  moduleServer = http.createServer((request, response) => {
    if (request.url === '/blank.html') {
      response.writeHead(200, { 'Content-Type': 'text/html; charset=utf-8' });
      response.end('<!doctype html><html><body></body></html>');
      return;
    }
    const relativePath = String(request.url || '').replace(/^\/+/, '');
    const filePath = path.resolve(process.cwd(), relativePath);
    if (!filePath.startsWith(path.resolve(process.cwd()) + path.sep) || !fs.existsSync(filePath)) {
      response.writeHead(404);
      response.end('Not found');
      return;
    }
    response.writeHead(200, { 'Content-Type': 'text/javascript; charset=utf-8' });
    response.end(fs.readFileSync(filePath));
  });
  await new Promise((resolve) => moduleServer.listen(0, '127.0.0.1', resolve));
  moduleOrigin = `http://127.0.0.1:${moduleServer.address().port}`;
});

test.afterAll(async () => {
  if (moduleServer) await new Promise((resolve) => moduleServer.close(resolve));
});

test('standalone search updates only hits and active feature parts at 25k features', async ({ page }, testInfo) => {
  test.setTimeout(120_000);
  const svgPath = testInfo.outputPath('interactive-search-25000.svg');
  const generator = String.raw`
import sys
from gbdraw.render.interactive_svg import InteractiveSvgContext, enrich_svg

count = 25_000
hits = {17, 1017, 5017, 15017, 24017}
paths = []
features = []
for index in range(count):
    svg_id = f"f{index}"
    label = f"needle feature {index}" if index in hits else f"feature {index}"
    x = index % 500
    y = index // 500
    paths.append(
        f'<rect id="{svg_id}" data-gbdraw-feature-id="{svg_id}" '
        f'x="{x}" y="{y}" width="0.8" height="0.8" fill="#54bcf8" />'
    )
    features.append({"svg_id": svg_id, "label": label, "type": "CDS"})
source = '<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 500 50">' + ''.join(paths) + '</svg>'
with open(sys.argv[1], 'w', encoding='utf-8') as handle:
    handle.write(enrich_svg(source, InteractiveSvgContext(features=features, popup_mode='simple')))
`;
  execFileSync('python', ['-c', generator, svgPath], { cwd: process.cwd(), stdio: 'inherit' });

  await page.addInitScript(() => {
    window.__featureClassMutations = [];
    const originalToggle = DOMTokenList.prototype.toggle;
    DOMTokenList.prototype.toggle = function patchedToggle(token, force) {
      if (String(token).includes('feature--match') || String(token).includes('active-match') || String(token).includes('dimmed')) {
        window.__featureClassMutations.push({ token: String(token), force: Boolean(force) });
      }
      return originalToggle.call(this, token, force);
    };
  });
  await page.goto(pathToFileURL(svgPath).href);
  await page.getByRole('button', { name: 'Expand feature search' }).click();
  await page.evaluate(() => { window.__featureClassMutations = []; });

  await page.locator('[data-search-query]').fill('needle');
  const started = Date.now();
  await page.locator('[data-search-apply]').click();
  await expect(page.locator('svg')).toHaveClass(/gbdraw-feature-search-active/);
  await page.evaluate(() => new Promise((resolve) => requestAnimationFrame(() => requestAnimationFrame(resolve))));
  expect(Date.now() - started).toBeLessThan(1000);
  await expect(page.locator('.gbdraw-interactive-feature--match')).toHaveCount(5);
  await expect(page.locator('.gbdraw-interactive-feature--active-match')).toHaveCount(1);
  await expect(page.locator('.gbdraw-interactive-feature--dimmed')).toHaveCount(0);
  const resultMutations = await page.evaluate(() => window.__featureClassMutations);
  expect(resultMutations.length).toBeLessThanOrEqual(6);

  await page.evaluate(() => { window.__featureClassMutations = []; });
  const navigationStarted = Date.now();
  await page.locator('[data-search-next]').click();
  expect(Date.now() - navigationStarted).toBeLessThan(250);
  const navigationMutations = await page.evaluate(() => window.__featureClassMutations);
  expect(navigationMutations).toHaveLength(2);
  expect(navigationMutations.every((entry) => entry.token.includes('active-match'))).toBeTruthy();

  await page.evaluate(() => { window.__featureClassMutations = []; });
  await page.locator('[data-search-clear]').click();
  await page.evaluate(() => new Promise((resolve) => requestAnimationFrame(() => requestAnimationFrame(resolve))));
  await expect(page.locator('svg')).not.toHaveClass(/gbdraw-feature-search-active/);
  const clearMutations = await page.evaluate(() => window.__featureClassMutations);
  expect(clearMutations.length).toBeLessThanOrEqual(6);
});

test('v2 derives FASTA and amino-acid search from translation while the runtime accepts v1', async ({ page }, testInfo) => {
  const v2Path = testInfo.outputPath('interactive-v2.svg');
  const v1Path = testInfo.outputPath('interactive-v1.svg');
  const generator = String.raw`
import json
import sys
import xml.etree.ElementTree as ET
from gbdraw.render.interactive_svg import InteractiveSvgContext, enrich_svg

source = '<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 100 80"><rect id="fseq" data-gbdraw-feature-id="fseq" x="10" y="10" width="20" height="10" fill="#54bcf8" /></svg>'
context = InteractiveSvgContext(features=[{
    "svg_id": "fseq",
    "record_id": "rec1",
    "type": "CDS",
    "start": 0,
    "end": 9,
    "strand": "+",
    "qualifiers": {
        "product": ["protein A"],
        "protein_id": ["WP_000001.1"],
        "translation": ["MPEPTIDE"],
    },
    "nucleotide_sequence": "ATGAAATAA",
    "amino_acid_sequence": "MPEPTIDE",
}])
v2 = enrich_svg(source, context)
with open(sys.argv[1], 'w', encoding='utf-8') as handle:
    handle.write(v2)
root = ET.fromstring(v2)
metadata = next(item for item in root.iter() if item.tag.rsplit('}', 1)[-1] == 'metadata')
payload = json.loads(metadata.text)
payload["schema"] = "gbdraw-interactive-feature-popup-v1"
feature = payload["features"][0]
feature["amino_acid_sequence"] = "MPEPTIDE"
feature["nucleotide_fasta"] = ">rec1:1-9 protein A\nATGAAATAA"
feature["amino_acid_fasta"] = ">WP_000001.1 protein A\nMPEPTIDE"
metadata.set("data-schema", payload["schema"])
metadata.text = json.dumps(payload, separators=(',', ':'))
with open(sys.argv[2], 'w', encoding='utf-8') as handle:
    handle.write(ET.tostring(root, encoding='unicode'))
`;
  execFileSync('python', ['-c', generator, v2Path, v1Path], { cwd: process.cwd(), stdio: 'inherit' });

  for (const [schema, svgPath] of [['v2', v2Path], ['v1', v1Path]]) {
    await page.goto(pathToFileURL(svgPath).href);
    await page.getByRole('button', { name: 'Expand feature search' }).click();
    await page.locator('[data-search-field]').selectOption('amino-acid');
    await page.locator('[data-search-query]').fill('MPEPTIDE');
    await page.locator('[data-search-apply]').click();
    await expect(page.locator('.gbdraw-interactive-feature--match')).toHaveCount(1);
    await page.locator('[data-search-open]').click();
    await page.locator('[data-tab="sequence"]').click();
    const sequences = await page.locator('.gfi-pre').allTextContents();
    expect(sequences.join('\n'), `${schema} nucleotide FASTA`).toContain('ATGAAATAA');
    expect(sequences.join('\n'), `${schema} amino-acid FASTA`).toContain('MPEPTIDE');
    expect(sequences.join('\n'), `${schema} protein ID`).toContain('WP_000001.1');
  }
});

test('v2 raw match metadata is materialized lazily for the popup', async ({ page }, testInfo) => {
  const svgPath = testInfo.outputPath('interactive-v2-match.svg');
  const generator = String.raw`
import sys
from gbdraw.render.interactive_svg import InteractiveSvgContext, enrich_svg

source = '''<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 100 80">
<rect id="fq" data-gbdraw-feature-id="fq" x="10" y="10" width="20" height="10" fill="#54bcf8" />
<rect id="fs" data-gbdraw-feature-id="fs" x="10" y="50" width="20" height="10" fill="#54bcf8" />
<path data-gbdraw-pairwise-match-id="m1" data-match-kind="pairwise" data-orthogroup-id="og1"
 data-query-record-id="rec1" data-subject-record-id="rec2" data-qstart="1" data-qend="9"
 data-sstart="10" data-send="18" data-query-feature-svg-id="fq" data-subject-feature-svg-id="fs"
 data-identity="99.1" data-alignment-length="9" fill="#94a3b8" d="M 20 20 L 30 20 L 30 50 L 20 50 Z" />
<path data-gbdraw-pairwise-match-id="m2" data-match-kind="collinear" data-orthogroup-id="og1"
 data-collinearity-block-id="block1" data-collinearity-block-kind="syntenic"
 data-collinear-group-scope="adjacent_local" data-group-kind="collinear_gene_group"
 data-query-record-id="rec1" data-subject-record-id="rec2" data-qstart="1" data-qend="9"
 data-sstart="10" data-send="18" data-query-feature-svg-id="fq" data-subject-feature-svg-id="fs"
 data-identity="95.0" data-alignment-length="9" fill="#64748b" d="M 40 20 L 50 20 L 50 50 L 40 50 Z" />
</svg>'''
features = [
 {"svg_id": "fq", "record_id": "rec1", "type": "CDS", "start": 0, "end": 9, "orthogroup_id": "og1", "qualifiers": {"product": ["Protein A"], "protein_id": ["P1"]}},
 {"svg_id": "fs", "record_id": "rec2", "type": "CDS", "start": 9, "end": 18, "orthogroup_id": "og1", "qualifiers": {"product": ["Protein B"], "protein_id": ["P2"]}},
]
orthogroups = [{"id": "og1", "name": "Family", "member_count": 2, "record_coverage_count": 2, "members": [
 {"featureSvgId": "fq", "recordId": "rec1", "sourceProteinId": "P1", "product": "Protein A"},
 {"featureSvgId": "fs", "recordId": "rec2", "sourceProteinId": "P2", "product": "Protein B"},
]}]
with open(sys.argv[1], 'w', encoding='utf-8') as handle:
    handle.write(enrich_svg(source, InteractiveSvgContext(features=features, orthogroups=orthogroups)))
`;
  execFileSync('python', ['-c', generator, svgPath], { cwd: process.cwd(), stdio: 'inherit' });
  await page.goto(pathToFileURL(svgPath).href);
  await page.locator('[data-gbdraw-pairwise-match-id="m1"]').click();
  await expect(page.locator('.gfi-title')).toHaveText('Pairwise match');
  const sectionTitles = await page.locator('.gfi-block-title').allTextContents();
  expect(sectionTitles).toEqual(['Summary', 'Alignment', 'Orthogroup', 'Query', 'Subject']);
  await expect(page.locator('.gfi-content')).toContainText('99.1');
  await expect(page.locator('.gfi-content')).toContainText('Protein A');
  await expect(page.locator('.gfi-content')).toContainText('Protein B');
  await page.locator('[data-close]').click();
  await page.locator('[data-gbdraw-pairwise-match-id="m2"]').click();
  const collinearTitles = await page.locator('.gfi-block-title').allTextContents();
  expect(collinearTitles).toEqual(['Summary', 'Local collinear groups', 'Collinearity', 'Query', 'Subject']);
  await expect(page.locator('.gfi-content')).toContainText('Number of local collinear groups');
});

test('web exporter writes the same compact v2 feature and raw match contract', async ({ page }) => {
  await page.goto(`${moduleOrigin}/blank.html`);
  const payload = await page.evaluate(async ({ origin }) => {
    const { enrichSvgWithStandaloneInteractivity } = await import(
      `${origin}/gbdraw/web/js/services/standalone-interactivity.js`
    );
    const svg = document.createElementNS('http://www.w3.org/2000/svg', 'svg');
    svg.setAttribute('viewBox', '0 0 100 80');
    svg.innerHTML = `
      <rect id="fq" data-gbdraw-feature-id="fq" x="1" y="1" width="10" height="5" fill="#54bcf8" />
      <rect id="fs" data-gbdraw-feature-id="fs" x="1" y="20" width="10" height="5" fill="#54bcf8" />
      <path data-gbdraw-pairwise-match-id="m1" data-match-kind="pairwise"
        data-orthogroup-id="og1" data-query-record-id="rec1" data-subject-record-id="rec2"
        data-qstart="1" data-qend="9" data-sstart="10" data-send="18"
        data-query-feature-svg-id="fq" data-subject-feature-svg-id="fs"
        data-identity="99.1" data-alignment-length="9" d="M 1 2 L 2 3" />`;
    enrichSvgWithStandaloneInteractivity(svg, {
      popupMode: 'rich',
      features: [
        {
          svg_id: 'fq', record_id: 'rec1', type: 'CDS', start: 0, end: 9,
          qualifiers: { translation: ['MPEPTIDE'] },
          nucleotide_sequence: 'ATGAAATAA', amino_acid_sequence: 'MPEPTIDE'
        },
        { svg_id: 'fs', record_id: 'rec2', type: 'CDS', start: 9, end: 18 }
      ]
    });
    return JSON.parse(svg.querySelector('#gbdraw-interactive-feature-metadata').textContent);
  }, { origin: moduleOrigin });

  expect(payload.schema).toBe('gbdraw-interactive-feature-popup-v2');
  expect(payload.features[0].nucleotide_fasta).toBeUndefined();
  expect(payload.features[0].amino_acid_fasta).toBeUndefined();
  expect(payload.features[0].amino_acid_sequence).toBeUndefined();
  expect(payload.features[0].qualifiers.translation).toEqual(['MPEPTIDE']);
  expect(payload.matches[0]).toMatchObject({
    id: 'm1',
    match_kind: 'pairwise',
    orthogroup_ids: ['og1'],
    query_record_id: 'rec1',
    subject_record_id: 'rec2',
    identity: '99.1',
    alignment_length: '9'
  });
  expect(payload.matches[0].sections).toBeUndefined();
  expect(payload.matches[0].hover_rows).toBeUndefined();
});

test('preview search renderer applies result and active differences only', async ({ page }) => {
  await page.goto(`${moduleOrigin}/blank.html`);
  const result = await page.evaluate(async ({ origin }) => {
    const preview = await import(`${origin}/gbdraw/web/js/app/feature-search/preview-svg.js`);
    const svg = document.createElementNS('http://www.w3.org/2000/svg', 'svg');
    for (let index = 0; index < 1000; index += 1) {
      const rect = document.createElementNS('http://www.w3.org/2000/svg', 'rect');
      rect.setAttribute('id', `f${index}`);
      rect.setAttribute('data-gbdraw-feature-id', `f${index}`);
      svg.appendChild(rect);
    }
    document.body.appendChild(svg);
    const mutations = [];
    const originalToggle = DOMTokenList.prototype.toggle;
    DOMTokenList.prototype.toggle = function patchedToggle(token, force) {
      if (String(token).includes('feature-search-match') || String(token).includes('active-match') || String(token).includes('dimmed')) {
        mutations.push({ token: String(token), force: Boolean(force) });
      }
      return originalToggle.call(this, token, force);
    };
    const state = preview.createPreviewFeatureSearchDomState();
    const featureIndex = preview.getPreviewFeatureElementIndex(svg);
    const matches = ['f1', 'f21', 'f41', 'f61', 'f81'];
    preview.schedulePreviewFeatureSearchClasses({
      svg, matches, activeId: 'f1', queryActive: true, featureIndex, appliedState: state
    });
    await new Promise((resolve) => requestAnimationFrame(() => requestAnimationFrame(resolve)));
    const resultMutationCount = mutations.length;
    mutations.length = 0;
    preview.applyPreviewActiveSearchMatch({ featureIndex, appliedState: state, activeId: 'f21' });
    const navigationMutationCount = mutations.length;
    DOMTokenList.prototype.toggle = originalToggle;
    return {
      resultMutationCount,
      navigationMutationCount,
      matchCount: svg.querySelectorAll('.gbdraw-preview-feature-search-match').length,
      dimmedCount: svg.querySelectorAll('.gbdraw-preview-feature-search-dimmed').length,
      rootActive: svg.classList.contains('gbdraw-preview-feature-search-results-active')
    };
  }, { origin: moduleOrigin });

  expect(result).toEqual({
    resultMutationCount: 6,
    navigationMutationCount: 2,
    matchCount: 5,
    dimmedCount: 0,
    rootActive: true
  });
});
