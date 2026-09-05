const { execFileSync } = require('node:child_process');
const fs = require('node:fs');
const { pathToFileURL } = require('node:url');
const { test, expect } = require('@playwright/test');

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
    features.append({"svg_id": svg_id, "record_idx": 0, "label": label, "type": "CDS"})
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
  const applyElapsedMs = await page.locator('[data-search-apply]').evaluate(async (button) => {
    const started = performance.now();
    button.click();
    await new Promise((resolve) => requestAnimationFrame(() => requestAnimationFrame(resolve)));
    return performance.now() - started;
  });
  await expect(page.locator('svg')).toHaveClass(/gbdraw-feature-search-active/);
  expect(applyElapsedMs).toBeLessThan(1000);
  await expect(page.locator('.gbdraw-interactive-feature--match')).toHaveCount(5);
  await expect(page.locator('.gbdraw-interactive-feature--active-match')).toHaveCount(1);
  await expect(page.locator('.gbdraw-interactive-feature--dimmed')).toHaveCount(0);
  const resultMutations = await page.evaluate(() => window.__featureClassMutations);
  expect(resultMutations.length).toBeLessThanOrEqual(6);

  await page.evaluate(() => { window.__featureClassMutations = []; });
  const navigationElapsedMs = await page.locator('[data-search-next]').evaluate((button) => {
    const started = performance.now();
    button.click();
    return performance.now() - started;
  });
  expect(navigationElapsedMs).toBeLessThan(250);
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

test('v3 derives FASTA and the runtime still accepts v1/v2 metadata', async ({ page }, testInfo) => {
  const v3Path = testInfo.outputPath('interactive-v3.svg');
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
    "record_idx": 0,
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
v3 = enrich_svg(source, context)
with open(sys.argv[1], 'w', encoding='utf-8') as handle:
    handle.write(v3)
root = ET.fromstring(v3)
metadata = next(item for item in root.iter() if item.tag.rsplit('}', 1)[-1] == 'metadata')
feature = {
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
}
payload = {
    "schema": "gbdraw-interactive-feature-popup-v2",
    "popup_mode": "rich",
    "features": [feature],
    "biological_features": [feature],
    "orthogroups": [],
    "matches": [],
}
metadata.set("data-schema", payload["schema"])
metadata.text = json.dumps(payload, separators=(',', ':'))
with open(sys.argv[2], 'w', encoding='utf-8') as handle:
    handle.write(ET.tostring(root, encoding='unicode'))
payload["schema"] = "gbdraw-interactive-feature-popup-v1"
feature["amino_acid_sequence"] = "MPEPTIDE"
feature["nucleotide_fasta"] = ">rec1:1-9 protein A\nATGAAATAA"
feature["amino_acid_fasta"] = ">WP_000001.1 protein A\nMPEPTIDE"
metadata.set("data-schema", payload["schema"])
metadata.text = json.dumps(payload, separators=(',', ':'))
with open(sys.argv[3], 'w', encoding='utf-8') as handle:
    handle.write(ET.tostring(root, encoding='unicode'))
`;
  execFileSync('python', ['-c', generator, v3Path, v2Path, v1Path], { cwd: process.cwd(), stdio: 'inherit' });

  for (const [schema, svgPath] of [['v3', v3Path], ['v2', v2Path], ['v1', v1Path]]) {
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
 data-query-protein-id="h_aaaaaaaaaaaaaaaaaaaaaaaaaa"
 data-subject-protein-id="f_bbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb"
 data-identity="99.1" data-alignment-length="9" fill="#94a3b8" d="M 20 20 L 30 20 L 30 50 L 20 50 Z" />
<path data-gbdraw-pairwise-match-id="m2" data-match-kind="collinear" data-orthogroup-id="og1"
 data-collinearity-block-id="block1" data-collinearity-block-kind="syntenic"
 data-collinear-group-scope="adjacent_local" data-group-kind="collinear_gene_group"
 data-query-record-id="rec1" data-subject-record-id="rec2" data-qstart="1" data-qend="9"
 data-sstart="10" data-send="18" data-query-feature-svg-id="fq" data-subject-feature-svg-id="fs"
 data-identity="95.0" data-alignment-length="9" fill="#64748b" d="M 40 20 L 50 20 L 50 50 L 40 50 Z" />
<path data-gbdraw-pairwise-match-id="m3" data-match-kind="orthogroup" data-orthogroup-id="og1"
 data-query-record-id="rec1" data-subject-record-id="rec2"
 data-query-feature-svg-id="fq" data-subject-feature-svg-id="fs"
 fill="#94a3b8" d="M 60 20 L 66 20 L 66 50 L 60 50 Z" />
<path data-gbdraw-pairwise-match-id="m4" data-match-kind="orthogroup" data-orthogroup-id="og1"
 data-query-record-id="rec1" data-subject-record-id="rec2"
 data-query-feature-svg-id="fq" data-subject-feature-svg-id="fs"
 fill="#94a3b8" d="M 74 20 L 80 20 L 80 50 L 74 50 Z" />
</svg>'''
features = [
 {"svg_id": "fq", "record_idx": 0, "record_id": "rec1", "type": "CDS", "start": 0, "end": 9, "orthogroup_id": "og1",
  "proteinId": "h_aaaaaaaaaaaaaaaaaaaaaaaaaa",
  "displayProteinId": "f_bbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb",
  "sourceProteinId": "record@instance|alias~f_cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc",
  "locus_tag": "WP_QUERY_LOCUS.1",
  "qualifiers": {"product": ["Protein A"], "protein_id": ["h_aaaaaaaaaaaaaaaaaaaaaaaaaa", "WP_QUERY.1"], "translation": ["MPEPTIDE"]},
  "nucleotide_sequence": "ATGAAATAA", "amino_acid_sequence": "MPEPTIDE"},
 {"svg_id": "fs", "record_idx": 1, "record_id": "rec2", "type": "CDS", "start": 9, "end": 18, "orthogroup_id": "og1",
  "proteinId": "h_dddddddddddddddddddddddddd",
  "sourceProteinId": "f_eeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeee",
  "locus_tag": "WP_SUBJECT.1",
  "qualifiers": {"product": ["Protein B"], "protein_id": ["h_dddddddddddddddddddddddddd"], "translation": ["MSUBJECT"]},
  "nucleotide_sequence": "ATGCCCTAA", "amino_acid_sequence": "MSUBJECT"},
]
orthogroups = [{"id": "og1", "name": "Family", "member_count": 2, "record_coverage_count": 2, "members": [
 {"featureSvgId": "fq", "stableFeatureSvgId": "fq", "recordIndex": 0, "recordId": "rec1", "proteinId": "h_aaaaaaaaaaaaaaaaaaaaaaaaaa",
  "displayProteinId": "f_bbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb",
  "sourceProteinId": "record@instance|alias~f_cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc",
  "label": "h_aaaaaaaaaaaaaaaaaaaaaaaaaa", "locusTag": "WP_MEMBER_QUERY.1", "product": "Protein A"},
 {"featureSvgId": "fs", "stableFeatureSvgId": "fs", "recordIndex": 1, "recordId": "rec2", "proteinId": "h_dddddddddddddddddddddddddd",
  "sourceProteinId": "f_eeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeee",
  "label": "h_dddddddddddddddddddddddddd", "locusTag": "WP_SUBJECT.1", "product": "Protein B"},
]}]
with open(sys.argv[1], 'w', encoding='utf-8') as handle:
    handle.write(enrich_svg(source, InteractiveSvgContext(features=features, orthogroups=orthogroups)))
`;
  execFileSync('python', ['-c', generator, svgPath], { cwd: process.cwd(), stdio: 'inherit' });
  await page.addInitScript(() => {
    window.__copiedText = '';
    Object.defineProperty(navigator, 'clipboard', {
      configurable: true,
      value: {
        writeText: async (value) => {
          window.__copiedText = String(value);
        }
      }
    });
  });
  await page.goto(pathToFileURL(svgPath).href);
  const firstMatch = page.locator('[data-gbdraw-pairwise-match-id="m1"]');
  const secondMatch = page.locator('[data-gbdraw-pairwise-match-id="m2"]');
  await expect(firstMatch).toHaveAttribute('role', 'button');
  await expect(firstMatch).toHaveAttribute('tabindex', '0');

  await firstMatch.hover();
  await page.mouse.down();
  const pendingStyle = await firstMatch.evaluate((element) => {
    const style = getComputedStyle(element);
    return {
      focused: element === document.activeElement,
      focusVisible: element.matches(':focus-visible'),
      outlineStyle: style.outlineStyle,
      stroke: style.stroke
    };
  });
  expect(pendingStyle).toEqual({
    focused: true,
    focusVisible: false,
    outlineStyle: 'none',
    stroke: 'rgb(245, 158, 11)'
  });
  await expect(firstMatch).toHaveClass(/gbdraw-interactive-pairwise-match--pending/);
  await expect(firstMatch).not.toHaveClass(/gbdraw-interactive-pairwise-match--selected/);
  await expect(page.locator('.gfi-title')).toHaveCount(0);
  await firstMatch.dispatchEvent('pointercancel', {
    bubbles: true,
    button: 0,
    isPrimary: true,
    pointerId: 1,
    pointerType: 'mouse'
  });
  await expect(firstMatch).not.toHaveClass(/gbdraw-interactive-pairwise-match--pending/);
  await page.mouse.move(1, 1);
  await page.mouse.up();

  await firstMatch.dispatchEvent('pointerdown', {
    bubbles: true,
    button: 0,
    isPrimary: true,
    pointerId: 2,
    pointerType: 'mouse'
  });
  await secondMatch.dispatchEvent('pointerdown', {
    bubbles: true,
    button: 0,
    isPrimary: true,
    pointerId: 3,
    pointerType: 'mouse'
  });
  await expect(firstMatch).not.toHaveClass(/gbdraw-interactive-pairwise-match--pending/);
  await expect(secondMatch).toHaveClass(/gbdraw-interactive-pairwise-match--pending/);
  await secondMatch.dispatchEvent('pointercancel', {
    bubbles: true,
    button: 0,
    isPrimary: true,
    pointerId: 3,
    pointerType: 'mouse'
  });
  await expect(secondMatch).not.toHaveClass(/gbdraw-interactive-pairwise-match--pending/);

  await firstMatch.click();
  await expect(firstMatch)
    .toHaveClass(/gbdraw-interactive-pairwise-match--selected/);
  await expect(firstMatch).not.toHaveClass(/gbdraw-interactive-pairwise-match--pending/);
  await expect(page.locator('.gfi-title')).toHaveText('Pairwise match');
  const sectionTitles = await page.locator('.gfi-block-title').allTextContents();
  expect(sectionTitles).toEqual([
    'Matched sequences',
    'Query span',
    'Subject span',
    'Summary',
    'Alignment',
    'Similarity group',
    'Query',
    'Subject',
  ]);
  await expect(page.locator('.gfi-content')).toContainText('99.1');
  await expect(page.locator('.gfi-content')).toContainText('Protein A');
  await expect(page.locator('.gfi-content')).toContainText('Protein B');
  await expect(page.locator('.gfi-content')).toContainText('WP_QUERY.1');
  await expect(page.locator('.gfi-content')).toContainText('WP_SUBJECT.1');
  const popupText = await page.locator('.gfi-content').innerText();
  expect(popupText).not.toMatch(/h_[a-z2-7]{26}/i);
  expect(popupText).not.toMatch(/f_[0-9a-f]{64}/i);
  expect(popupText).not.toContain('record@instance|alias~');
  await page.locator('.gfi-match-feature-table').first().getByRole('button', { name: 'Copy' }).click();
  const copiedText = await page.evaluate(() => window.__copiedText);
  expect(copiedText).toContain('WP_QUERY.1');
  expect(copiedText).not.toMatch(/h_[a-z2-7]{26}/i);
  expect(copiedText).not.toMatch(/f_[0-9a-f]{64}/i);
  const queryMemberRow = page.locator('.gfi-og-members-table tbody tr').filter({ hasText: 'WP_QUERY.1' });
  const downloadPromise = page.waitForEvent('download');
  await queryMemberRow.getByRole('button', { name: 'DL aa' }).click();
  const download = await downloadPromise;
  expect(download.suggestedFilename()).toContain('WP_QUERY.1');
  expect(download.suggestedFilename()).not.toMatch(/h_[a-z2-7]{26}/i);
  expect(download.suggestedFilename()).not.toMatch(/f_[0-9a-f]{64}/i);
  await page.locator('[data-close]').click();
  await expect(page.locator('[data-gbdraw-pairwise-match-id="m1"]'))
    .not.toHaveClass(/gbdraw-interactive-pairwise-match--selected/);
  await secondMatch.click();
  const collinearTitles = await page.locator('.gfi-block-title').allTextContents();
  expect(collinearTitles).toEqual([
    'Collinear block spans',
    'Query span',
    'Subject span',
    'Summary',
    'Local collinear groups',
    'Collinearity',
    'Query',
    'Subject',
  ]);
  await expect(page.locator('.gfi-content')).toContainText('Number of local collinear groups');
  await page.locator('[data-close]').click();
  await page.locator('[data-gbdraw-pairwise-match-id="m3"]').click();
  await expect(page.locator('[data-gbdraw-pairwise-match-id="m3"]'))
    .toHaveClass(/gbdraw-interactive-pairwise-match--selected/);
  await expect(page.locator('[data-gbdraw-pairwise-match-id="m4"]'))
    .toHaveClass(/gbdraw-interactive-pairwise-match--selected/);
  await expect(page.locator('[data-gbdraw-pairwise-match-id="m1"]'))
    .not.toHaveClass(/gbdraw-interactive-pairwise-match--selected/);
  await expect.poll(() => page.locator('[data-gbdraw-pairwise-match-id="m3"]')
    .evaluate((element) => getComputedStyle(element).outlineStyle))
    .toBe('none');

  await page.locator('[data-close]').click();
  await firstMatch.focus();
  await page.keyboard.press('Shift+Tab');
  await page.keyboard.press('Tab');
  await expect(firstMatch).toBeFocused();
  const keyboardFocusStyle = await firstMatch.evaluate((element) => {
    const style = getComputedStyle(element);
    return {
      focusVisible: element.matches(':focus-visible'),
      outlineColor: style.outlineColor,
      outlineStyle: style.outlineStyle,
      outlineWidth: style.outlineWidth
    };
  });
  expect(keyboardFocusStyle).toEqual({
    focusVisible: true,
    outlineColor: 'rgb(37, 99, 235)',
    outlineStyle: 'solid',
    outlineWidth: '2px'
  });
  await firstMatch.press('Enter');
  await expect(page.locator('.gfi-title')).toHaveText('Pairwise match');
  await expect(firstMatch).toHaveClass(/gbdraw-interactive-pairwise-match--selected/);
  await expect.poll(() => firstMatch.evaluate((element) => getComputedStyle(element).outlineStyle))
    .toBe('solid');
  await page.locator('[data-close]').click();

  await secondMatch.focus();
  await secondMatch.press('Space');
  await expect(page.locator('.gfi-title')).toHaveText('Collinearity block');
  await expect(secondMatch).toHaveClass(/gbdraw-interactive-pairwise-match--selected/);
});

test('standalone homology popup exports exact spans and keeps missing comparison optional', async ({ page }, testInfo) => {
  const svgPath = testInfo.outputPath('interactive-homology-match.svg');
  const generator = String.raw`
import sys
from gbdraw.render.interactive_svg import InteractiveSvgContext, enrich_svg

source = '''<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 120 80">
<path data-ring-background="true" fill="#eeeeee" d="M 5 5 L 115 5 L 115 15 L 5 15 Z" />
<path data-gbdraw-match-id="homology_ring2_hit17" data-match-kind="homology"
 data-source-index="0" data-track-index="2" data-track-label="Comparison A"
 data-reference-side="query" data-reference-record-id="ref"
 data-query-record-id="ref" data-subject-record-id="cmp"
 data-qstart="2" data-qend="5" data-sstart="6" data-send="3"
 data-identity="98.5" data-alignment-length="4" data-evalue="1e-20"
 fill="#ef4444" d="M 10 25 L 110 25 L 110 35 L 10 35 Z" />
<path data-gbdraw-match-id="homology_ring3_hit1" data-match-kind="homology"
 data-source-index="1" data-track-index="3" data-track-label="BLAST only"
 data-reference-side="query" data-reference-record-id="ref"
 data-query-record-id="ref" data-subject-record-id="missing"
 data-qstart="1" data-qend="4" data-sstart="1" data-send="4"
 fill="#3b82f6" d="M 10 45 L 110 45 L 110 55 L 10 55 Z" />
</svg>'''
sources = [
 {"key": "circular:record:0", "recordId": "ref", "aliases": ["REF"],
  "sequence": "AACCGGTT", "origin": "circular-reference", "recordIndex": 0},
 {"key": "homology:comparison:0:cmp", "recordId": "cmp", "aliases": ["CMP"],
  "sequence": "TTGCAACC", "origin": "homology-comparison", "sourceIndex": 0},
]
with open(sys.argv[1], 'w', encoding='utf-8') as handle:
    handle.write(enrich_svg(source, InteractiveSvgContext(sequence_sources=sources)))
`;
  execFileSync('python', ['-c', generator, svgPath], { cwd: process.cwd(), stdio: 'inherit' });
  await page.addInitScript(() => {
    window.__copiedText = '';
    Object.defineProperty(navigator, 'clipboard', {
      configurable: true,
      value: { writeText: async (value) => { window.__copiedText = String(value); } },
    });
  });
  await page.goto(pathToFileURL(svgPath).href);

  await page.locator('[data-ring-background="true"]').click();
  await expect(page.locator('.gfi-title')).toHaveCount(0);
  await page.locator('[data-gbdraw-match-id="homology_ring2_hit17"]').click();
  await expect(page.locator('.gfi-title')).toHaveText('Homology ring match');
  await expect(page.locator('.gfi-content')).toContainText('Comparison A');
  await expect(page.locator('.gfi-content')).toContainText('Reference span');
  await expect(page.locator('.gfi-content')).toContainText('Comparison span');

  const referenceBlock = page.locator('.gfi-block').filter({ hasText: 'Reference span' }).last();
  await referenceBlock.getByRole('button', { name: 'Copy' }).click();
  await expect.poll(() => page.evaluate(() => window.__copiedText)).toBe(
    '>homology_ring2_hit17_query|record=ref|coords=2..5|strand=+\nACCG\n'
  );

  const downloadPromise = page.waitForEvent('download');
  await page.locator('.gfi-block-actions').getByRole('button', { name: 'FASTA' }).click();
  const download = await downloadPromise;
  expect(download.suggestedFilename()).toBe('homology_ring2_hit17_both.fna');
  expect(await fs.promises.readFile(await download.path(), 'utf8')).toBe(
    '>homology_ring2_hit17_query|record=ref|coords=2..5|strand=+\nACCG\n' +
    '>homology_ring2_hit17_subject|record=cmp|coords=6..3|strand=-\nTTGC\n'
  );

  await page.locator('[data-close]').click();
  await page.locator('[data-gbdraw-match-id="homology_ring3_hit1"]').click();
  await expect(page.locator('.gfi-content')).toContainText(
    'Comparison sequence was not supplied for this BLAST source.'
  );
  await expect(page.locator('.gfi-block-actions')).toHaveCount(0);
});

test('preview search renderer applies result and active differences only', async ({ page }) => {
  await page.goto('/');
  const origin = new URL(page.url()).origin;
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
    const repeatedFeatureIndex = preview.getPreviewFeatureElementIndex(svg);
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
      featureIndexReused: repeatedFeatureIndex === featureIndex,
      matchCount: svg.querySelectorAll('.gbdraw-preview-feature-search-match').length,
      dimmedCount: svg.querySelectorAll('.gbdraw-preview-feature-search-dimmed').length,
      rootActive: svg.classList.contains('gbdraw-preview-feature-search-results-active')
    };
  }, { origin });

  expect(result).toEqual({
    resultMutationCount: 6,
    navigationMutationCount: 2,
    featureIndexReused: true,
    matchCount: 5,
    dimmedCount: 0,
    rootActive: true
  });
});
