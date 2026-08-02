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

test('browser export embeds the exact selected schema-3 item and expands references', async ({
  page
}, testInfo) => {
  await page.goto(`${moduleOrigin}/blank.html`);
  const exported = await page.evaluate(async ({ origin }) => {
    const { enrichSvgWithStandaloneInteractivity } = await import(
      `${origin}/gbdraw/web/js/services/standalone-interactivity.js`
    );
    const fullNote = `${'x'.repeat(49)}😀tail`;
    const catalog = {
      schema: 3,
      items: [{
        resultIndex: 0,
        resultName: 'diagram.svg',
        recordKeys: ['record-key-a', 'record-key-b'],
        features: [{
          svgId: 'rendered-visible',
          recordKey: 'record-key-a',
          biologicalFeatureId: 'stable-visible',
          fillColor: '#54bcf8'
        }, {
          svgId: 'rendered-collision',
          recordKey: 'record-key-a',
          biologicalFeatureId: 'h_aaaaaaaaaaaaaaaaaaaaaaaaaa',
          fillColor: '#f59e0b'
        }],
        biologicalFeatures: [
          {
            recordKey: 'record-key-a',
            biologicalFeatureId: 'stable-visible',
            record_id: 'rec-visible',
            type: 'CDS',
            start: 0,
            end: 9,
            strand: '+',
            aminoAcidSequence: 'MVISIBLE',
            translationFromAminoAcidSequence: true,
            qualifiers: {
              locus_tag: ['VP_1'],
              protein_id: ['VP_1'],
              old_locus_tag: ['OLD_VP_1'],
              gene: ['visible_gene'],
              product: ['Visible protein'],
              note: [fullNote]
            },
            sequenceSourceIndex: 0
          },
          {
            recordKey: 'record-key-a',
            biologicalFeatureId: 'h_aaaaaaaaaaaaaaaaaaaaaaaaaa',
            type: 'CDS',
            start: 20,
            end: 29,
            strand: '+',
            qualifiers: {
              note: [fullNote]
            },
            sequenceSourceIndex: 0
          },
          {
            recordKey: 'record-key-b',
            biologicalFeatureId: 'stable-hidden',
            record_id: 'rec-hidden',
            type: 'CDS',
            start: 9,
            end: 18,
            strand: '+',
            location_parts: [
              { start: 9, end: 12, strand: '+' },
              { start: 12, end: 18, strand: '+' }
            ],
            product: 'Hidden override',
            amino_acid_sequence: 'MHIDDEN',
            translationFromAminoAcidSequence: true,
            qualifiers: {
              locus_tag: ['HP_1'],
              protein_id: ['HP_1'],
              product: ['Hidden protein']
            },
            sequenceSourceIndex: 1
          }
        ],
        orthogroups: [{
          id: 'og-hidden',
          name: 'hidden-test',
          description: 'Original group description',
          member_count: 2,
          record_coverage_count: 2,
          members: [
            {
              recordKey: 'record-key-a',
              biologicalFeatureId: 'stable-visible',
              representative: true
            },
            {
              recordKey: 'record-key-b',
              biologicalFeatureId: 'stable-hidden'
            }
          ]
        }],
        annotations: [{
          dom_id: 'annotation-review-window',
          id: 'review-window',
          set_id: 'review',
          track_id: 'annotations-1',
          record_id: 'rec-visible',
          record_index: 0,
          segments: [[2, 8]],
          label: 'Review window',
          mark: 'band',
          lane: 0,
          metadata: { reviewer: 'Ada' }
        }],
        comparisonMatches: [],
        sequenceSources: [{
          key: 'linear:record:0',
          origin: 'linear-record',
          recordIndex: 0,
          sequence: `ATGAAATAA${'N'.repeat(11)}ATGCCCTAA`
        }, {
          key: 'linear:record:1',
          origin: 'linear-record',
          recordIndex: 1,
          sequence: `${'N'.repeat(9)}ATGCCCTAA`
        }]
      }, {
        resultIndex: 1,
        resultName: 'other.svg',
        recordKeys: [],
        features: [],
        biologicalFeatures: [],
        orthogroups: [],
        annotations: [],
        comparisonMatches: []
      }]
    };
    for (let index = 0; index < 128; index += 1) {
      catalog.items[0].biologicalFeatures.push({
        recordKey: 'record-key-a',
        biologicalFeatureId: `bulk-feature-${index}`,
        type: 'CDS',
        start: 0,
        end: 3,
        strand: '+',
        qualifiers: {},
        sequenceSourceIndex: 0
      });
    }
    const svg = document.createElementNS('http://www.w3.org/2000/svg', 'svg');
    svg.setAttribute('xmlns', 'http://www.w3.org/2000/svg');
    svg.setAttribute('viewBox', '0 0 120 80');
    svg.innerHTML = `
      <rect id="rendered-visible" data-gbdraw-feature-id="shared-stable"
        data-gbdraw-rendered-feature-id="rendered-visible"
        x="5" y="5" width="25" height="12" fill="#54bcf8" />
      <rect id="rendered-collision" data-gbdraw-feature-id="shared-stable"
        data-gbdraw-rendered-feature-id="rendered-collision"
        x="35" y="5" width="25" height="12" fill="#f59e0b" />
      <g id="annotation-review-window" data-gbdraw-annotation-id="review-window"
        data-gbdraw-annotation-set-id="review"
        data-gbdraw-annotation-track-id="annotations-1"
        data-gbdraw-record-id="rec-visible" data-gbdraw-record-index="0"
        data-gbdraw-annotation-mark="band" data-gbdraw-annotation-label="Review window">
        <rect x="5" y="25" width="20" height="6" fill="#94a3b8" />
      </g>`;
    const enriched = enrichSvgWithStandaloneInteractivity(svg, {
      popupMode: 'rich',
      featureCatalog: catalog,
      catalogResultIndex: 0,
      catalogResultName: 'diagram.svg',
      requireFeatureCatalog: true,
      labelTextFeatureOverrides: {
        'rendered-visible': 'Edited visible label'
      },
      orthogroupNameOverrides: {
        'og-hidden': 'Edited similarity group'
      },
      orthogroupDescriptionOverrides: {
        'og-hidden': 'Edited group description'
      }
    });
    const metadata = svg.querySelector('#gbdraw-interactive-feature-metadata');
    return {
      enriched,
      catalog,
      embedded: JSON.parse(metadata.textContent),
      schema: metadata.getAttribute('data-schema'),
      resultIndex: metadata.getAttribute('data-result-index'),
      resultName: metadata.getAttribute('data-result-name'),
      sourceDisplayLabel: catalog.items[0].features[0].displayLabel,
      sourceGroup: catalog.items[0].orthogroups[0],
      svgText: new XMLSerializer().serializeToString(svg)
    };
  }, { origin: moduleOrigin });

  expect(exported.enriched).toBe(true);
  expect(exported.embedded.schema).toBe(3);
  expect(exported.embedded.items).toHaveLength(1);
  expect(exported.embedded.items[0].features[0].displayLabel)
    .toBe('Edited visible label');
  expect(exported.embedded.items[0].orthogroups[0].display_name)
    .toBe('Edited similarity group');
  expect(exported.embedded.items[0].orthogroups[0].description)
    .toBe('Edited group description');
  expect(exported.embedded.items[0].annotations[0].id).toBe('review-window');
  expect(exported.sourceDisplayLabel).toBeUndefined();
  expect(exported.sourceGroup.display_name).toBeUndefined();
  expect(exported.sourceGroup.description).toBe('Original group description');
  expect(exported.schema).toBe('3');
  expect(exported.resultIndex).toBe('0');
  expect(exported.resultName).toBe('diagram.svg');
  expect(exported.svgText.match(/data-gbdraw-interactive-feature="true"/g)).toHaveLength(2);

  const svgPath = testInfo.outputPath('interactive-v3.svg');
  fs.writeFileSync(svgPath, exported.svgText, 'utf8');
  await page.addInitScript(() => {
    window.__copiedText = '';
    window.__expandedCatalogFeatures = {};
    window.__sourceValidationScans = 0;
    window.__sourceValidationScanCounts = {};
    const sourceSequences = new Set([
      `ATGAAATAA${'N'.repeat(11)}ATGCCCTAA`,
      `${'N'.repeat(9)}ATGCCCTAA`
    ]);
    const nativeRegexTest = RegExp.prototype.test;
    RegExp.prototype.test = function (value) {
      if (this.source === '\\s' && sourceSequences.has(value)) {
        window.__sourceValidationScans += 1;
        window.__sourceValidationScanCounts[value] = (
          window.__sourceValidationScanCounts[value] || 0
        ) + 1;
      }
      return nativeRegexTest.call(this, value);
    };
    const nativeMapSet = Map.prototype.set;
    Map.prototype.set = function (key, value) {
      const biologicalFeatureId = String(
        value && (
          value.biologicalFeatureId || value.biological_feature_id
        ) || ''
      );
      if (
        biologicalFeatureId
        && value.qualifiers
        && Array.isArray(value.location_parts)
      ) {
        window.__expandedCatalogFeatures[biologicalFeatureId] = {
          aminoAcidSequence: (
            value.amino_acid_sequence || value.aminoAcidSequence
          ),
          gene: value.gene,
          locationParts: value.location_parts,
          locusTag: value.locus_tag,
          note: value.note,
          oldLocusTag: value.old_locus_tag,
          product: value.product,
          proteinId: value.protein_id,
          translation: value.qualifiers && value.qualifiers.translation,
          translationMarker: value.translationFromAminoAcidSequence
        };
      }
      return nativeMapSet.call(this, key, value);
    };
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

  const expandedFeatures = await page.evaluate(
    () => window.__expandedCatalogFeatures
  );
  expect(await page.evaluate(() => window.__sourceValidationScans)).toBe(2);
  expect(await page.evaluate(() => Object.values(
    window.__sourceValidationScanCounts
  ))).toEqual([1, 1]);
  expect(expandedFeatures['stable-visible']).toMatchObject({
    aminoAcidSequence: 'MVISIBLE',
    gene: 'visible_gene',
    locusTag: 'VP_1',
    oldLocusTag: 'OLD_VP_1',
    product: 'Visible protein',
    proteinId: 'VP_1',
    note: `${'x'.repeat(49)}😀`,
    translation: ['MVISIBLE'],
    locationParts: [{
      start: 0,
      end: 9,
      strand: '+',
      display: '1..9'
    }]
  });
  expect(expandedFeatures['stable-visible'].translationMarker).toBeUndefined();
  expect(expandedFeatures['stable-hidden'].product).toBe('Hidden override');

  await page.locator('[data-gbdraw-rendered-feature-id="rendered-visible"]').click();
  await expect(page.locator('.gfi-title')).toContainText('Edited visible label');
  await expect(page.locator('#gbdraw-feature-popup')).toContainText(
    'Edited similarity group'
  );
  const memberBlock = page.locator('.gfi-block').filter({
    hasText: 'Similarity-group members'
  }).last();
  await expect(memberBlock.locator('tbody tr')).toHaveCount(2);
  await expect(memberBlock).toContainText('Visible protein');
  await expect(memberBlock).toContainText('Hidden override');
  await memberBlock.locator('.gfi-block-actions').getByRole('button', {
    name: 'Copy aa (2)'
  }).click();
  await expect.poll(() => page.evaluate(() => window.__copiedText)).toContain('>VP_1');
  const copied = await page.evaluate(() => window.__copiedText);
  expect(copied).toContain('>HP_1');
  expect(copied).toContain('MHIDDEN');
  await memberBlock.locator('.gfi-block-actions').getByRole('button', {
    name: 'Copy nt (2)'
  }).click();
  const copiedNucleotide = await page.evaluate(() => window.__copiedText);
  expect(copiedNucleotide).toContain('ATGAAATAA');
  expect(copiedNucleotide).toContain('ATGCCCTAA');

  await page.locator('[data-close]').click();
  await page.locator('[data-gbdraw-rendered-feature-id="rendered-collision"]').click();
  await page.getByRole('button', { name: 'Sequence' }).click();
  const nucleotideBlock = page.locator('.gfi-block').filter({
    hasText: 'Nucleotide'
  }).last();
  await nucleotideBlock.getByRole('button', { name: 'Copy', exact: true }).click();
  const nucleotideFasta = await page.evaluate(() => window.__copiedText);
  expect(nucleotideFasta).toContain('>record:21..29');
  expect(nucleotideFasta).toContain('ATGCCCTAA');
  expect(nucleotideFasta).not.toMatch(/h_[a-z2-7]{26}/i);

  await page.locator('[data-close]').click();
  const annotation = page.locator('#annotation-review-window');
  await expect(annotation).toHaveAttribute('data-gbdraw-interactive-annotation', 'true');
  await annotation.click();
  await expect(page.locator('.gfi-title')).toContainText('Review window');
  await expect(page.locator('#gbdraw-feature-popup')).toContainText('3..8');
  const labelRow = page.locator('.gfi-row').filter({ hasText: 'Label' });
  await labelRow.getByRole('button', { name: 'Copy' }).click();
  await expect.poll(() => page.evaluate(() => window.__copiedText))
    .toBe('Review window');

  await page.locator('[data-close]').click();
  await page.getByRole('button', { name: 'Expand feature search' }).click();
  await page.locator('[data-search-query]').fill('Visible protein');
  await page.locator('[data-search-apply]').click();
  await expect(page.locator('.gbdraw-interactive-feature--match')).toHaveCount(1);
});

test('standalone rejects conflicting compact provenance markers', async ({
  page
}, testInfo) => {
  await page.goto(`${moduleOrigin}/blank.html`);
  const variants = await page.evaluate(async ({ origin }) => {
    const { enrichSvgWithStandaloneInteractivity } = await import(
      `${origin}/gbdraw/web/js/services/standalone-interactivity.js`
    );
    const makeCatalog = () => ({
      schema: 3,
      items: [{
        resultIndex: 0,
        resultName: 'invalid.svg',
        recordKeys: ['record-key'],
        features: [{
          svgId: 'rendered-invalid',
          recordKey: 'record-key',
          biologicalFeatureId: 'invalid-feature'
        }],
        biologicalFeatures: [{
          recordKey: 'record-key',
          biologicalFeatureId: 'invalid-feature',
          record_id: 'record-id',
          type: 'CDS',
          start: 0,
          end: 3,
          strand: '+',
          sequenceSourceIndex: 0,
          amino_acid_sequence: 'M',
          translationFromAminoAcidSequence: true,
          qualifiers: {}
        }],
        orthogroups: [],
        annotations: [],
        comparisonMatches: [],
        sequenceSources: [{
          origin: 'linear-record',
          recordIndex: 0,
          sequence: 'ATG'
        }]
      }]
    });
    const cases = [];
    for (const conflict of [
      'DIFFERENT', '', ' ', null, [], [null], 0, false
    ]) {
      const catalog = makeCatalog();
      catalog.items[0].biologicalFeatures[0]
        .qualifiers.translation = conflict;
      cases.push({ name: `translation-${cases.length}`, catalog });
    }
    for (const invalidAminoAcid of [0, false, {}, []]) {
      const catalog = makeCatalog();
      catalog.items[0].biologicalFeatures[0]
        .amino_acid_sequence = invalidAminoAcid;
      cases.push({ name: `amino-${cases.length}`, catalog });
    }
    for (const shadowingValue of [null, '']) {
      const catalog = makeCatalog();
      const feature = catalog.items[0].biologicalFeatures[0];
      feature.aminoAcidSequence = feature.amino_acid_sequence;
      feature.amino_acid_sequence = shadowingValue;
      cases.push({ name: `amino-alias-${cases.length}`, catalog });
    }
    for (const invalidSequence of [123, {}, 'AT G']) {
      const catalog = makeCatalog();
      catalog.items[0].sequenceSources[0].sequence = invalidSequence;
      cases.push({ name: `sequence-${cases.length}`, catalog });
    }
    const unreferencedInvalidSource = makeCatalog();
    unreferencedInvalidSource.items[0].sequenceSources.push({
      origin: 'linear-record',
      recordIndex: 0,
      sequence: 'AT G'
    });
    cases.push({
      name: 'sequence-unreferenced',
      catalog: unreferencedInvalidSource
    });
    const coexisting = makeCatalog();
    coexisting.items[0].biologicalFeatures[0].nucleotide_sequence = 'ATG';
    cases.push({ name: 'coexisting-sequence', catalog: coexisting });

    return cases.map(({ name, catalog }) => {
      const svg = document.createElementNS('http://www.w3.org/2000/svg', 'svg');
      svg.setAttribute('xmlns', 'http://www.w3.org/2000/svg');
      svg.setAttribute('viewBox', '0 0 100 50');
      svg.innerHTML = `
        <rect id="rendered-invalid"
          data-gbdraw-feature-id="invalid-feature"
          data-gbdraw-rendered-feature-id="rendered-invalid"
          x="5" y="5" width="20" height="10" />`;
      enrichSvgWithStandaloneInteractivity(svg, {
        popupMode: 'rich',
        featureCatalog: catalog,
        catalogResultIndex: 0,
        catalogResultName: 'invalid.svg',
        requireFeatureCatalog: true
      });
      return {
        name,
        svgText: new XMLSerializer().serializeToString(svg)
      };
    });
  }, { origin: moduleOrigin });

  await page.addInitScript(() => {
    window.__expandedInvalidFeature = false;
    const nativeMapSet = Map.prototype.set;
    Map.prototype.set = function (key, value) {
      if (
        value
        && (
          value.biologicalFeatureId === 'invalid-feature'
          || value.biological_feature_id === 'invalid-feature'
        )
      ) {
        window.__expandedInvalidFeature = true;
      }
      return nativeMapSet.call(this, key, value);
    };
  });

  for (const variant of variants) {
    const svgPath = testInfo.outputPath(`${variant.name}.svg`);
    fs.writeFileSync(svgPath, variant.svgText, 'utf8');
    await page.goto(pathToFileURL(svgPath).href);
    await expect(page.locator('#gbdraw-feature-search-controls')).toBeAttached();
    expect(await page.evaluate(() => window.__expandedInvalidFeature)).toBe(false);
  }
});

test('Download Interactive SVG forwards live editor overrides without mutating the catalog', async ({
  page
}) => {
  await page.goto(`${moduleOrigin}/blank.html`);
  await page.addScriptTag({
    url: `${moduleOrigin}/gbdraw/web/vendor/vue/vue.global.js`
  });
  await page.addScriptTag({
    url: `${moduleOrigin}/gbdraw/web/vendor/dompurify/purify.min.js`
  });

  const exported = await page.evaluate(async ({ origin }) => {
    const { state } = await import(`${origin}/gbdraw/web/js/state.js`);
    const { downloadInteractiveSVG } = await import(
      `${origin}/gbdraw/web/js/services/export.js`
    );
    const catalog = {
      schema: 3,
      items: [{
        resultIndex: 0,
        resultName: 'live.svg',
        recordKeys: ['record-a'],
        features: [{
          svgId: 'rendered-a',
          recordKey: 'record-a',
          biologicalFeatureId: 'biological-a'
        }],
        biologicalFeatures: [{
          recordKey: 'record-a',
          biologicalFeatureId: 'biological-a',
          type: 'CDS',
          start: 0,
          end: 9,
          product: 'Original feature label'
        }],
        orthogroups: [{
          id: 'group-a',
          name: 'Original group',
          description: 'Original description',
          members: [{
            recordKey: 'record-a',
            biologicalFeatureId: 'biological-a'
          }]
        }],
        annotations: [],
        comparisonMatches: []
      }]
    };
    const container = document.createElement('div');
    container.innerHTML = `
      <svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 100 50">
        <rect id="rendered-a" data-gbdraw-feature-id="biological-a"
          data-gbdraw-rendered-feature-id="rendered-a"
          x="5" y="5" width="20" height="10" />
      </svg>`;
    document.body.appendChild(container);

    state.results.value = [{ name: 'live.svg', content: container.innerHTML }];
    state.selectedResultIndex.value = 0;
    state.svgContainer.value = container;
    state.featureCatalog.value = catalog;
    state.labelTextFeatureOverrides['rendered-a'] = 'Live feature label';
    state.orthogroupNameOverrides['group-a'] = 'Live group name';
    state.orthogroupDescriptionOverrides['group-a'] = 'Live description';

    let downloadedBlob = null;
    const originalCreateObjectURL = URL.createObjectURL;
    const originalRevokeObjectURL = URL.revokeObjectURL;
    const originalClick = HTMLAnchorElement.prototype.click;
    URL.createObjectURL = (blob) => {
      downloadedBlob = blob;
      return 'blob:gbdraw-test';
    };
    URL.revokeObjectURL = () => {};
    HTMLAnchorElement.prototype.click = () => {};
    try {
      downloadInteractiveSVG();
      const svgText = await downloadedBlob.text();
      const doc = new DOMParser().parseFromString(svgText, 'image/svg+xml');
      const metadata = doc.querySelector('#gbdraw-interactive-feature-metadata');
      return {
        embedded: JSON.parse(metadata.textContent),
        sourceFeature: catalog.items[0].features[0],
        sourceGroup: catalog.items[0].orthogroups[0]
      };
    } finally {
      URL.createObjectURL = originalCreateObjectURL;
      URL.revokeObjectURL = originalRevokeObjectURL;
      HTMLAnchorElement.prototype.click = originalClick;
    }
  }, { origin: moduleOrigin });

  expect(exported.embedded.items[0].features[0].displayLabel)
    .toBe('Live feature label');
  expect(exported.embedded.items[0].orthogroups[0].display_name)
    .toBe('Live group name');
  expect(exported.embedded.items[0].orthogroups[0].description)
    .toBe('Live description');
  expect(exported.sourceFeature.displayLabel).toBeUndefined();
  expect(exported.sourceGroup.display_name).toBeUndefined();
  expect(exported.sourceGroup.description).toBe('Original description');
});
