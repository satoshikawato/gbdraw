const fs = require('node:fs');
const { join } = require('node:path');
const { pathToFileURL } = require('node:url');
const { gunzipSync } = require('node:zlib');
const { test, expect } = require('@playwright/test');
const {
  generateAndWaitForResult,
  openApp,
  waitForAppShell
} = require('./helpers/app-lifecycle.cjs');

const installDiagramRequestObserver = (page) => page.addInitScript(() => {
  window.__GBDRAW_DIAGRAM_RUNS__ = [];
  const NativeWorker = window.Worker;
  window.Worker = new Proxy(NativeWorker, {
    construct(target, args) {
      const worker = Reflect.construct(target, args, target);
      if (!String(args[0] || '').includes('diagram-generation-worker.js')) return worker;
      const nativePostMessage = worker.postMessage.bind(worker);
      worker.postMessage = (message, transfer) => {
        if (message?.type === 'run' && message?.payload?.request) {
          window.__GBDRAW_DIAGRAM_RUNS__.push(
            JSON.parse(JSON.stringify(message.payload.request))
          );
        }
        if (transfer === undefined) return nativePostMessage(message);
        return nativePostMessage(message, transfer);
      };
      return worker;
    }
  });
});

const makeCircularRecord = (recordId, gene, start, end) => `LOCUS       ${recordId.padEnd(24)} 360 bp    DNA     circular UNA 01-JAN-2000
DEFINITION  feature popup record rotation acceptance.
ACCESSION   ${recordId}
VERSION     ${recordId}
KEYWORDS    .
SOURCE      synthetic construct
  ORGANISM  synthetic construct
            .
FEATURES             Location/Qualifiers
     source          1..360
     CDS             ${start}..${end}
                     /gene="${gene}"
                     /locus_tag="${gene}"
                     /product="${gene} protein"
                     /translation="MKKKKKKKKK"
ORIGIN
        1 ${'acgt'.repeat(15)}
       61 ${'acgt'.repeat(15)}
      121 ${'acgt'.repeat(15)}
      181 ${'acgt'.repeat(15)}
      241 ${'acgt'.repeat(15)}
      301 ${'acgt'.repeat(15)}
//
`;

test('feature popup record rotation works by pointer and keyboard in rich and simple layouts', async ({
  page
}) => {
  test.setTimeout(240000);
  page.on('dialog', (dialog) => dialog.dismiss());
  await installDiagramRequestObserver(page);
  await openApp(page);
  await page.getByLabel('GenBank/DDBJ File', { exact: true }).setInputFiles(join(
    process.cwd(), 'tests/test_inputs/HmmtDNA.gbk'
  ));
  await generateAndWaitForResult(page);
  const before = await page.evaluate(async () => {
    const { state } = await import('/gbdraw/web/js/state.js');
    return {
      svg: window.__GBDRAW_APP__.svgContent,
      prefix: window.__GBDRAW_DIAGRAM_RUNS__.at(-1).output.prefix,
      history: window.__GBDRAW_HISTORY__.getUndoCount(),
      biological: state.biologicalFeatures.value.map((feature) => [
        feature.record_key,
        feature.biological_feature_id
      ])
    };
  });
  await page.evaluate(() => {
    const app = window.__GBDRAW_APP__;
    app.adv.rich_feature_popup = true;
    app.form.prefix = 'UNRELATED_PENDING_PREFIX';
  });

  const search = page.getByRole('searchbox', { name: 'Search features', exact: true });
  await search.fill('tRNA');
  await search.press('Enter');
  await page.getByRole('button', { name: 'Open active feature', exact: true }).click();
  const disclosure = page.getByRole('button', { name: /Record actions · Rotate record/ });
  await expect(disclosure).toHaveAttribute('aria-expanded', 'false');
  await disclosure.click();
  const actions = page.getByRole('region', { name: 'Record actions' });
  await expect(actions).toBeVisible();
  await expect(actions).toContainText('Coordinates refer to the original record.');
  await expect(actions.getByLabel('Record rotation anchor')).toHaveValue('five-prime');
  await expect(actions.getByLabel('Record rotation signed offset')).toHaveValue('0');

  await actions.getByLabel('Record rotation anchor').selectOption('midpoint');
  await actions.getByLabel('Record rotation signed offset').fill('2');
  await disclosure.click();
  await expect(actions).toBeHidden();
  await disclosure.click();
  await expect(actions.getByLabel('Record rotation anchor')).toHaveValue('midpoint');
  await expect(actions.getByLabel('Record rotation signed offset')).toHaveValue('2');
  const expectedStart = await page.evaluate(() => (
    window.__GBDRAW_APP__.featureRecordRotationDraft.startCoordinate
  ));
  expect(expectedStart).toBeGreaterThan(0);
  await actions.getByRole('button', { name: 'Apply and regenerate' }).click();
  await expect(actions.locator('[aria-live="polite"]')).toContainText('Regenerating');
  await page.waitForFunction(() => !window.__GBDRAW_APP__.processing, null, {
    timeout: 240000
  });
  await expect(actions.locator('[aria-live="polite"]')).toContainText(
    'Record rotation applied and regenerated.'
  );
  expect(await search.inputValue()).toBe('tRNA');
  expect(await page.evaluate(() => window.__GBDRAW_APP__.form.prefix))
    .toBe('UNRELATED_PENDING_PREFIX');
  expect(await page.evaluate(() => {
    const controls = window.__GBDRAW_APP__.recordDisplayControls;
    const rows = controls.rows?.value || controls.rows || [];
    const row = rows.find((entry) => entry.recordId === 'NC_012920.1') || rows[0];
    const draft = controls.draftFor(row);
    return {
      startCoordinate: draft.startCoordinate,
      anchorFeatureId: draft.anchorIntent?.biologicalFeatureId || ''
    };
  })).toMatchObject({
    startCoordinate: expectedStart
  });
  const committed = await page.evaluate(async () => {
    const app = window.__GBDRAW_APP__;
    const { state } = await import('/gbdraw/web/js/state.js');
    const request = window.__GBDRAW_DIAGRAM_RUNS__.at(-1);
    return {
      schema: request.schema,
      startCoordinate: request.records[0].display.startCoordinate,
      prefix: request.output.prefix,
      history: window.__GBDRAW_HISTORY__.getUndoCount(),
      biological: state.biologicalFeatures.value.map((feature) => [
        feature.record_key,
        feature.biological_feature_id
      ]),
      target: {
        recordKey: app.featureRecordRotationDraft.identity.recordKey,
        biologicalFeatureId: app.featureRecordRotationDraft.identity.biologicalFeatureId
      }
    };
  });
  expect(committed).toMatchObject({
    schema: 8,
    startCoordinate: expectedStart,
    prefix: before.prefix,
    history: before.history + 1
  });
  expect(committed.biological).toEqual(before.biological);
  expect(committed.target.recordKey).toBeTruthy();
  expect(committed.target.biologicalFeatureId).toBeTruthy();
  const transformedSvg = await page.evaluate(() => window.__GBDRAW_APP__.svgContent);
  expect(transformedSvg).not.toBe(before.svg);

  await page.getByRole('button', { name: 'Close feature popup', exact: true }).click();
  await page.evaluate(() => window.__GBDRAW_HISTORY__.undo());
  await expect.poll(() => page.evaluate(() => window.__GBDRAW_APP__.svgContent)).toBe(before.svg);
  expect(await page.evaluate(async () => {
    const { state } = await import('/gbdraw/web/js/state.js');
    return state.recordDisplayDrafts.length;
  })).toBe(0);
  await page.evaluate(() => window.__GBDRAW_HISTORY__.redo());
  await expect.poll(() => page.evaluate(() => window.__GBDRAW_APP__.svgContent)).toBe(transformedSvg);
  expect(await page.evaluate(async () => {
    const { state } = await import('/gbdraw/web/js/state.js');
    return state.recordDisplayDrafts[0];
  })).toMatchObject({
    startCoordinate: expectedStart,
    anchorIntent: {
      schema: 1,
      placement: 'anchor',
      anchor: 'midpoint',
      offsetBp: 2
    }
  });

  await page.setViewportSize({ width: 390, height: 844 });
  await page.evaluate(() => { window.__GBDRAW_APP__.adv.rich_feature_popup = false; });
  await page.getByRole('button', { name: 'Open active feature', exact: true }).click();
  const simplePopup = page.locator('.feature-popup--simple');
  await expect(simplePopup).toBeVisible();
  const bounds = await simplePopup.boundingBox();
  expect(bounds.x).toBeGreaterThanOrEqual(0);
  expect(bounds.x + bounds.width).toBeLessThanOrEqual(390);

  const simpleDisclosure = simplePopup.getByRole('button', { name: /Record actions · Rotate record/ });
  await expect(simpleDisclosure).toHaveAttribute('aria-expanded', 'false');
  await simpleDisclosure.focus();
  await page.keyboard.press('Enter');
  const simpleActions = simplePopup.getByRole('region', { name: 'Record actions' });
  const anchor = simpleActions.getByLabel('Record rotation anchor');
  await anchor.focus();
  await anchor.press('Enter');
  await anchor.press('ArrowDown');
  await anchor.press('Enter');
  await expect(anchor).toHaveValue('midpoint');
  const offset = simpleActions.getByLabel('Record rotation signed offset');
  await offset.focus();
  await page.keyboard.press('ControlOrMeta+A');
  await page.keyboard.type('-3');
  await expect(offset).toHaveValue('-3');
  const cancel = simpleActions.getByRole('button', { name: 'Cancel', exact: true });
  await cancel.focus();
  await page.keyboard.press('Enter');
  await expect(simplePopup).toBeVisible();
  await expect(simpleDisclosure).toHaveAttribute('aria-expanded', 'false');
  await expect(simpleActions).toBeHidden();
  expect(await page.evaluate(() => window.__GBDRAW_APP__.svgContent)).toBe(transformedSvg);
  await simpleDisclosure.click();
  const staleActions = simplePopup.getByRole('region', { name: 'Record actions' });
  const staleBefore = await page.evaluate(async () => {
    const { state } = await import('/gbdraw/web/js/state.js');
    const target = window.__GBDRAW_APP__.featureRecordRotationDraft.identity;
    const duplicate = (features) => {
      const feature = features.find((entry) => (
        entry.record_key === target.recordKey
        && entry.biological_feature_id === target.biologicalFeatureId
      ));
      features.push(JSON.parse(JSON.stringify(feature)));
    };
    duplicate(state.extractedFeatures.value);
    duplicate(state.biologicalFeatures.value);
    return {
      svg: window.__GBDRAW_APP__.svgContent,
      history: window.__GBDRAW_HISTORY__.getUndoCount(),
      drafts: JSON.stringify(state.recordDisplayDrafts)
    };
  });
  await staleActions.getByRole('button', { name: 'Apply and regenerate' }).click();
  await expect(staleActions.locator('[aria-live="polite"]')).toContainText(
    'no longer present'
  );
  expect(await page.evaluate(async () => {
    const { state } = await import('/gbdraw/web/js/state.js');
    return {
      svg: window.__GBDRAW_APP__.svgContent,
      history: window.__GBDRAW_HISTORY__.getUndoCount(),
      drafts: JSON.stringify(state.recordDisplayDrafts)
    };
  })).toEqual(staleBefore);
  await page.evaluate(async () => {
    const { state } = await import('/gbdraw/web/js/state.js');
    state.extractedFeatures.value.pop();
    state.biologicalFeatures.value.pop();
  });
  await page.getByRole('button', { name: 'Close feature popup', exact: true }).click();
  await page.setViewportSize({ width: 1280, height: 900 });

  const pendingSave = page.waitForEvent('download');
  await page.evaluate(async () => {
    window.__GBDRAW_APP__.sessionTitle = 'feature-popup-record-rotation';
    await window.__GBDRAW_APP__.saveSessionWithTitle();
  });
  const savedPath = await (await pendingSave).path();
  const saved = JSON.parse(gunzipSync(fs.readFileSync(savedPath)));
  expect(saved.version).toBe(44);
  expect(saved.renderRequest.schema).toBe(8);
  expect(saved.renderRequest.records[0].display.startCoordinate).toBe(expectedStart);
  expect(saved.config.recordDisplayDrafts[0].anchorIntent).toMatchObject({
    schema: 1,
    placement: 'anchor',
    anchor: 'midpoint',
    offsetBp: 2
  });
  await page.reload({ waitUntil: 'domcontentloaded' });
  await waitForAppShell(page);
  await page.locator('input[accept^=".json,"]').first().setInputFiles(savedPath);
  await expect.poll(() => page.evaluate(() => window.__GBDRAW_APP__.sessionImportPending), {
    timeout: 180000
  }).toBe(false);
  expect(await page.evaluate(() => window.__GBDRAW_APP__.svgContent)).toBe(transformedSvg);
  expect(await page.evaluate(async () => {
    const { state } = await import('/gbdraw/web/js/state.js');
    return state.recordDisplayDrafts[0];
  })).toMatchObject({
    startCoordinate: expectedStart,
    anchorIntent: {
      schema: 1,
      recordKey: committed.target.recordKey,
      biologicalFeatureId: committed.target.biologicalFeatureId
    }
  });
});

test('same-file record rotation keeps chromosome targets independent', async ({ page }) => {
  test.setTimeout(240000);
  await installDiagramRequestObserver(page);
  await openApp(page);
  await page.getByLabel('GenBank/DDBJ File', { exact: true }).setInputFiles({
    name: 'two-chromosomes.gbk',
    mimeType: 'text/plain',
    buffer: Buffer.from(
      makeCircularRecord('chromosome_I', 'dnaA', 21, 105)
      + makeCircularRecord('chromosome_II', 'parB', 151, 225)
    )
  });
  await generateAndWaitForResult(page);
  const baseline = await page.evaluate(async () => {
    const { state } = await import('/gbdraw/web/js/state.js');
    return {
      svg: window.__GBDRAW_APP__.svgContent,
      biological: state.biologicalFeatures.value.map((feature) => [
        feature.record_key,
        feature.biological_feature_id
      ])
    };
  });

  const rotateFromPopup = async (query, anchor, offset) => {
    const search = page.getByRole('searchbox', { name: 'Search features', exact: true });
    await search.fill(query);
    await search.press('Enter');
    await page.getByRole('button', { name: 'Open active feature', exact: true }).click();
    await page.getByRole('button', { name: /Record actions · Rotate record/ }).click();
    const actions = page.getByRole('region', { name: 'Record actions' });
    await expect(actions).toBeVisible();
    await actions.getByLabel('Record rotation anchor').selectOption(anchor);
    await actions.getByLabel('Record rotation signed offset').fill(String(offset));
    const draft = await page.evaluate(() => ({
      identity: { ...window.__GBDRAW_APP__.featureRecordRotationDraft.identity },
      startCoordinate: window.__GBDRAW_APP__.featureRecordRotationDraft.startCoordinate
    }));
    await actions.getByRole('button', { name: 'Apply and regenerate' }).click();
    await page.waitForFunction(() => !window.__GBDRAW_APP__.processing, null, {
      timeout: 240000
    });
    await expect(actions.locator('[aria-live="polite"]')).toContainText(
      'Record rotation applied and regenerated.'
    );
    const snapshot = await page.evaluate(async () => {
      const { state } = await import('/gbdraw/web/js/state.js');
      return {
        request: JSON.parse(JSON.stringify(window.__GBDRAW_DIAGRAM_RUNS__.at(-1))),
        svg: window.__GBDRAW_APP__.svgContent,
        biological: state.biologicalFeatures.value.map((feature) => [
          feature.record_key,
          feature.biological_feature_id
        ])
      };
    });
    await page.getByRole('button', { name: 'Close feature popup', exact: true }).click();
    return { draft, snapshot };
  };

  const chromosomeI = await rotateFromPopup('dnaA', 'five-prime', -5);
  expect(chromosomeI.snapshot.request.records).toHaveLength(2);
  const firstIndex = chromosomeI.snapshot.request.records.findIndex((record) => (
    record.recordKey === chromosomeI.draft.identity.recordKey
  ));
  expect(firstIndex).toBeGreaterThanOrEqual(0);
  const secondIndex = firstIndex === 0 ? 1 : 0;
  expect(chromosomeI.snapshot.request.records[firstIndex].display.startCoordinate)
    .toBe(chromosomeI.draft.startCoordinate);
  expect(chromosomeI.snapshot.request.records[secondIndex].display.startCoordinate).toBeNull();
  expect(chromosomeI.snapshot.svg).not.toBe(baseline.svg);
  expect(chromosomeI.snapshot.biological).toEqual(baseline.biological);

  const preservedFirstRecord = chromosomeI.snapshot.request.records[firstIndex];
  const chromosomeII = await rotateFromPopup('parB', 'midpoint', 7);
  const secondTargetIndex = chromosomeII.snapshot.request.records.findIndex((record) => (
    record.recordKey === chromosomeII.draft.identity.recordKey
  ));
  expect(secondTargetIndex).toBe(secondIndex);
  expect(chromosomeII.snapshot.request.records[firstIndex]).toEqual(preservedFirstRecord);
  expect(chromosomeII.snapshot.request.records[secondIndex].display.startCoordinate)
    .toBe(chromosomeII.draft.startCoordinate);
  expect(chromosomeII.snapshot.request.records.map((record) => record.presentation.gridRow))
    .toEqual(chromosomeI.snapshot.request.records.map((record) => record.presentation.gridRow));
  expect(chromosomeII.snapshot.request.tracks).toEqual(chromosomeI.snapshot.request.tracks);
  expect(chromosomeII.snapshot.request.comparisons)
    .toEqual(chromosomeI.snapshot.request.comparisons);
  expect(chromosomeII.snapshot.svg).not.toBe(chromosomeI.snapshot.svg);
  expect(chromosomeII.snapshot.biological).toEqual(baseline.biological);
});

test('both modes record rotation resolves the same circular source anchor', async ({ browser }) => {
  test.setTimeout(300000);
  const source = makeCircularRecord('shared_anchor', 'anchor_gene', 41, 125);
  const rotateInMode = async (mode) => {
    const context = await browser.newContext();
    const page = await context.newPage();
    await installDiagramRequestObserver(page);
    await openApp(page);
    if (mode === 'linear') {
      await page.getByRole('button', { name: 'Linear', exact: true }).click();
    }
    const upload = mode === 'linear'
      ? page.getByTestId('linear-genbank-1')
      : page.getByLabel('GenBank/DDBJ File', { exact: true });
    await upload.setInputFiles({
      name: 'shared-anchor.gbk',
      mimeType: 'text/plain',
      buffer: Buffer.from(source)
    });
    await generateAndWaitForResult(page);
    const beforeSvg = await page.evaluate(() => window.__GBDRAW_APP__.svgContent);
    const search = page.getByRole('searchbox', { name: 'Search features', exact: true });
    await search.fill('anchor_gene');
    await search.press('Enter');
    await page.getByRole('button', { name: 'Open active feature', exact: true }).click();
    await page.getByRole('button', { name: /Record actions · Rotate record/ }).click();
    const actions = page.getByRole('region', { name: 'Record actions' });
    await actions.getByLabel('Record rotation anchor').selectOption('midpoint');
    await actions.getByLabel('Record rotation signed offset').fill('-11');
    const expected = await page.evaluate(() => ({
      identity: { ...window.__GBDRAW_APP__.featureRecordRotationDraft.identity },
      startCoordinate: window.__GBDRAW_APP__.featureRecordRotationDraft.startCoordinate
    }));
    await actions.getByRole('button', { name: 'Apply and regenerate' }).click();
    await page.waitForFunction(() => !window.__GBDRAW_APP__.processing, null, {
      timeout: 240000
    });
    await expect(actions.locator('[aria-live="polite"]')).toContainText(
      'Record rotation applied and regenerated.'
    );
    const accepted = await page.evaluate(async (recordKey) => {
      const { state } = await import('/gbdraw/web/js/state.js');
      const request = window.__GBDRAW_DIAGRAM_RUNS__.at(-1);
      const record = request.records.find((entry) => entry.recordKey === recordKey);
      const draft = state.recordDisplayDrafts.find((entry) => (
        entry.anchorIntent?.recordKey === recordKey
      ));
      return {
        schema: request.schema,
        startCoordinate: record?.display?.startCoordinate,
        anchorIntent: draft?.anchorIntent,
        svg: window.__GBDRAW_APP__.svgContent
      };
    }, expected.identity.recordKey);
    expect(accepted.schema).toBe(8);
    expect(accepted.startCoordinate).toBe(expected.startCoordinate);
    expect(accepted.anchorIntent).toMatchObject({
      schema: 1,
      recordKey: expected.identity.recordKey,
      biologicalFeatureId: expected.identity.biologicalFeatureId,
      anchor: 'midpoint',
      offsetBp: -11
    });
    expect(accepted.svg).not.toBe(beforeSvg);
    await context.close();
    return accepted;
  };

  const circular = await rotateInMode('circular');
  const linear = await rotateInMode('linear');
  expect(linear.startCoordinate).toBe(circular.startCoordinate);
  expect(linear.anchorIntent.anchor).toBe(circular.anchorIntent.anchor);
  expect(linear.anchorIntent.offsetBp).toBe(circular.anchorIntent.offsetBp);
});

test('browser export embeds the exact selected schema-4 item and expands references', async ({
  page
}, testInfo) => {
  await page.goto('/');
  const origin = new URL(page.url()).origin;
  const exported = await page.evaluate(async ({ origin }) => {
    const { enrichSvgWithStandaloneInteractivity } = await import(
      `${origin}/gbdraw/web/js/services/standalone-interactivity.js`
    );
    const fullNote = `${'x'.repeat(49)}😀tail`;
    const catalog = {
      schema: 4,
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
  }, { origin });

  expect(exported.enriched).toBe(true);
  expect(exported.embedded.schema).toBe(4);
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
  expect(exported.schema).toBe('4');
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
  await page.clock.install();

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
  const groupCopyButtons = memberBlock.locator(
    '.gfi-block-actions [data-copy-feedback-key]'
  );
  const groupNtCopy = groupCopyButtons.nth(0);
  const groupAaCopy = groupCopyButtons.nth(1);
  await expect(groupNtCopy).toHaveText('Copy nt (2)');
  await expect(groupAaCopy).toHaveText('Copy aa (2)');
  await groupAaCopy.click();
  await expect.poll(() => page.evaluate(() => window.__copiedText)).toContain('>VP_1');
  const copied = await page.evaluate(() => window.__copiedText);
  expect(copied).toContain('>HP_1');
  expect(copied).toContain('MHIDDEN');
  await expect(groupAaCopy).toHaveText('Copied!');
  await expect(groupAaCopy).toHaveAccessibleName('Copied!');
  await expect(groupAaCopy).toHaveAttribute('aria-live', 'polite');
  await expect(groupAaCopy).toHaveAttribute('aria-atomic', 'true');
  await expect(groupNtCopy).toHaveText('Copy nt (2)');

  await page.getByRole('button', { name: 'Qualifiers', exact: true }).click();
  await page.getByRole('button', { name: 'Details', exact: true }).click();
  await expect(groupAaCopy).toHaveText('Copied!');
  await page.clock.fastForward(1500);
  await expect(groupAaCopy).toHaveText('Copy aa (2)');

  const firstMemberCopyButtons = memberBlock.locator(
    'tbody tr'
  ).first().locator('[data-copy-feedback-key]');
  const secondMemberAaCopy = memberBlock.locator(
    'tbody tr'
  ).nth(1).locator('[data-copy-feedback-key]').nth(1);
  const firstMemberAaCopy = firstMemberCopyButtons.nth(1);
  await firstMemberAaCopy.click();
  await expect(firstMemberAaCopy).toHaveText('Copied!');
  await expect(secondMemberAaCopy).toHaveText('Copy aa');
  expect(await page.evaluate(() => window.__copiedText)).toContain('>VP_1');
  expect(await page.evaluate(() => window.__copiedText)).not.toContain('>HP_1');
  await page.clock.fastForward(1500);
  await expect(firstMemberAaCopy).toHaveText('Copy aa');

  await groupNtCopy.click();
  await expect(groupNtCopy).toHaveText('Copied!');
  await expect(groupAaCopy).toHaveText('Copy aa (2)');
  const copiedNucleotide = await page.evaluate(() => window.__copiedText);
  expect(copiedNucleotide).toContain('ATGAAATAA');
  expect(copiedNucleotide).toContain('ATGCCCTAA');
  await page.clock.fastForward(1500);

  await page.evaluate(() => {
    window.__manualCopy = null;
    Object.defineProperty(navigator, 'clipboard', {
      configurable: true,
      value: { writeText: async () => { throw new Error('denied'); } }
    });
    window.prompt = (_title, value) => {
      window.__manualCopy = String(value);
      return null;
    };
  });
  await groupAaCopy.click();
  await expect(groupAaCopy).toHaveText('Copy manually');
  await expect(groupAaCopy).not.toHaveText('Copied!');
  expect(await page.evaluate(() => window.__manualCopy)).toBe(copied);
  await page.clock.fastForward(1500);
  await expect(groupAaCopy).toHaveText('Copy aa (2)');

  await page.evaluate(() => {
    window.__manualCopy = null;
    Object.defineProperty(navigator, 'clipboard', {
      configurable: true,
      value: {}
    });
  });
  await groupAaCopy.click();
  await expect(groupAaCopy).toHaveText('Copy manually');
  expect(await page.evaluate(() => window.__manualCopy)).toBe(copied);
  await page.clock.fastForward(1500);

  await page.evaluate(() => {
    Object.defineProperty(navigator, 'clipboard', {
      configurable: true,
      value: { writeText: async () => { throw new Error('denied'); } }
    });
    window.prompt = () => { throw new Error('prompt failed'); };
  });
  await groupAaCopy.click();
  await expect(groupAaCopy).toHaveText('Copy failed');
  await expect(groupAaCopy).not.toHaveText('Copied!');
  await page.clock.fastForward(1500);

  await page.evaluate(() => {
    Object.defineProperty(navigator, 'clipboard', {
      configurable: true,
      value: {
        writeText: async (value) => { window.__copiedText = String(value); }
      }
    });
  });
  await groupAaCopy.click();
  await page.clock.fastForward(1000);
  await groupAaCopy.click();
  await page.clock.fastForward(600);
  await expect(groupAaCopy).toHaveText('Copied!');
  await page.clock.fastForward(900);
  await expect(groupAaCopy).toHaveText('Copy aa (2)');

  await page.evaluate(() => {
    window.__copyResolvers = [];
    Object.defineProperty(navigator, 'clipboard', {
      configurable: true,
      value: {
        writeText: (value) => {
          window.__copiedText = String(value);
          return new Promise((resolve, reject) => {
            window.__copyResolvers.push({ reject, resolve });
          });
        }
      }
    });
  });
  await groupAaCopy.click();
  await groupAaCopy.click();
  await page.evaluate(() => window.__copyResolvers[1].resolve());
  await expect(groupAaCopy).toHaveText('Copied!');
  await page.evaluate(() => window.__copyResolvers[0].reject(new Error('stale denial')));
  await expect(groupAaCopy).toHaveText('Copied!');
  await page.clock.fastForward(1500);
  await expect(groupAaCopy).toHaveText('Copy aa (2)');

  await page.evaluate(() => {
    Object.defineProperty(navigator, 'clipboard', {
      configurable: true,
      value: {
        writeText: async (value) => { window.__copiedText = String(value); }
      }
    });
  });
  await groupAaCopy.click();
  await expect(groupAaCopy).toHaveText('Copied!');
  await page.locator('[data-close]').click();
  await page.locator('[data-gbdraw-rendered-feature-id="rendered-visible"]').click();
  const reopenedMemberBlock = page.locator('.gfi-block').filter({
    hasText: 'Similarity-group members'
  }).last();
  await expect(reopenedMemberBlock.locator(
    '.gfi-block-actions [data-copy-feedback-key]'
  ).nth(1)).toHaveText('Copy aa (2)');

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
  await page.goto('/');
  const origin = new URL(page.url()).origin;
  const variants = await page.evaluate(async ({ origin }) => {
    const { enrichSvgWithStandaloneInteractivity } = await import(
      `${origin}/gbdraw/web/js/services/standalone-interactivity.js`
    );
    const makeCatalog = () => ({
      schema: 4,
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
  }, { origin });

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
  await page.goto('/');
  const origin = new URL(page.url()).origin;
  await page.addScriptTag({
    url: '/gbdraw/web/vendor/vue/vue.global.js'
  });
  await page.addScriptTag({
    url: '/gbdraw/web/vendor/dompurify/purify.min.js'
  });

  const exported = await page.evaluate(async ({ origin }) => {
    const { state } = await import(`${origin}/gbdraw/web/js/state.js`);
    const { downloadInteractiveSVG } = await import(
      `${origin}/gbdraw/web/js/services/export.js`
    );
    const { captureSvgExport } = await import(
      `${origin}/gbdraw/web/js/services/svg-serialization.js`
    );
    const catalog = {
      schema: 4,
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
      await downloadInteractiveSVG(captureSvgExport(state, { interactive: true }));
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
  }, { origin });

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
