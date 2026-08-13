const { test, expect } = require('@playwright/test');
const { readFileSync } = require('node:fs');
const { join, resolve } = require('node:path');
const { gunzipSync } = require('node:zlib');

const repoRoot = resolve(process.env.GBDRAW_REPO || process.cwd());

const makeComparisonGenbank = (recordId, base = 'atg') => {
  const sequence = base.repeat(100);
  const origin = sequence.match(/.{1,60}/g).map((chunk, index) => {
    const groups = chunk.match(/.{1,10}/g).join(' ');
    return `${String(index * 60 + 1).padStart(9)} ${groups}`;
  }).join('\n');
  return `LOCUS       ${recordId.padEnd(24)} 300 bp    DNA     linear   UNA 01-JAN-2000
DEFINITION  linear comparison browser test.
ACCESSION   ${recordId}
VERSION     ${recordId}
KEYWORDS    .
SOURCE      synthetic construct
  ORGANISM  synthetic construct
            .
FEATURES             Location/Qualifiers
     CDS             1..90
                     /product="comparison test protein"
ORIGIN
${origin}
//
`;
};

const linearRecordCard = (page, uid) => (
  page.locator(`[data-linear-record-card="${uid}"]`)
);

const linearComparisonBoundary = (page, upperRow, lowerRow) => (
  page.locator(`[data-linear-comparison-boundary="${upperRow}->${lowerRow}"]`)
);

const linearComparisonPair = (page, edgeKey) => (
  page.locator(`fieldset[data-edge-key="${edgeKey}"]`)
);

const linearSelectedPairs = (page) => (
  page.locator('details[data-linear-comparison-disclosure="selected-pairs"]')
);

const openLinearSelectedPairs = async (page) => {
  const details = linearSelectedPairs(page);
  await expect(details).toHaveCount(1);
  if (!await details.evaluate((element) => element.open)) {
    await details.locator(':scope > summary').press('Enter');
  }
  await expect(details).toHaveAttribute('open', '');
  return details;
};

const openLinearAdvancedComparison = async (page) => {
  const details = page.locator(
    'details[data-linear-comparison-disclosure="advanced"]'
  );
  await expect(details).toHaveCount(1);
  if (!await details.evaluate((element) => element.open)) {
    await details.locator(':scope > summary').press('Enter');
  }
  await expect(details).toHaveAttribute('open', '');
  return details;
};

const expectExactEdgeKeys = async (container, expectedKeys) => {
  const pairs = container.locator('fieldset[data-edge-key]');
  await expect(pairs).toHaveCount(expectedKeys.length);
  const actualKeys = await pairs.evaluateAll((elements) => (
    elements.map((element) => element.dataset.edgeKey)
  ));
  expect(new Set(actualKeys).size).toBe(actualKeys.length);
  expect([...actualKeys].sort()).toEqual([...expectedKeys].sort());
};

const focusLastTabStop = async (container) => {
  await container.evaluate((element) => {
    const candidates = [...element.querySelectorAll(
      'a[href], button, input, select, textarea, [tabindex]'
    )].filter((candidate) => (
      !candidate.disabled &&
      candidate.tabIndex >= 0 &&
      candidate.getClientRects().length > 0
    ));
    const target = candidates.at(-1);
    if (!target) throw new Error('Expected at least one tab stop.');
    target.focus();
  });
};

const installDiagramRequestObserver = async (page) => {
  await page.addInitScript(() => {
    window.__GBDRAW_DIAGRAM_RUNS__ = [];
    window.__GBDRAW_RUNTIME_URLS__ = { fetches: [], workers: [] };
    const nativeFetch = window.fetch.bind(window);
    window.fetch = (...args) => {
      window.__GBDRAW_RUNTIME_URLS__.fetches.push(String(args[0]?.url || args[0] || ''));
      return nativeFetch(...args);
    };
    const NativeWorker = window.Worker;
    window.Worker = new Proxy(NativeWorker, {
      construct(target, args) {
        const worker = Reflect.construct(target, args, target);
        const workerUrl = String(args[0] || '');
        window.__GBDRAW_RUNTIME_URLS__.workers.push(workerUrl);
        if (workerUrl.includes('diagram-generation-worker.js')) {
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
        }
        return worker;
      }
    });
  });
};

test('Similarity-group popup selects every OG match without a focus outline', async ({ page }) => {
  await page.goto('/gbdraw/web/index.html', { waitUntil: 'domcontentloaded' });
  await page.waitForFunction(() => window.__GBDRAW_APP__);

  await page.evaluate(async () => {
    const app = window.__GBDRAW_APP__;
    const { ingestSvgResult } = await import('/gbdraw/web/js/services/svg-result-ingestion.js');
    app.mode = 'linear';
    app.results.splice(0, app.results.length, ingestSvgResult({
      name: 'pairwise-selection.svg',
      content: `<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 100 80">
        <path data-gbdraw-pairwise-match-id="match-1" data-match-kind="orthogroup"
          data-orthogroup-id="og-1"
          data-query-record-id="query" data-subject-record-id="subject"
          data-qstart="1" data-qend="20" data-sstart="5" data-send="24"
          fill="#94a3b8" d="M 5 20 L 25 20 L 25 60 L 5 60 Z" />
        <path data-gbdraw-pairwise-match-id="match-2" data-match-kind="orthogroup"
          data-orthogroup-id="og-1" fill="#94a3b8"
          d="M 35 20 L 50 20 L 50 60 L 35 60 Z" />
        <path data-gbdraw-pairwise-match-id="match-3" data-match-kind="orthogroup"
          data-orthogroup-id="og-2" fill="#94a3b8"
          d="M 60 20 L 75 20 L 75 60 L 60 60 Z" />
        <path data-gbdraw-pairwise-match-id="match-4" data-match-kind="pairwise"
          data-orthogroup-id="og-1"
          fill="#94a3b8" d="M 82 20 L 95 20 L 95 60 L 82 60 Z" />
      </svg>`
    }));
    app.selectedResultIndex = 0;
  });

  const firstOgMatch = page.getByRole('button', { name: 'Pairwise match 1', exact: true });
  const secondOgMatch = page.getByRole('button', { name: 'Pairwise match 2', exact: true });
  const otherOgMatch = page.getByRole('button', { name: 'Pairwise match 3', exact: true });
  const pairwiseMatch = page.getByRole('button', { name: 'Pairwise match 4', exact: true });
  await expect(firstOgMatch).toBeVisible();
  await firstOgMatch.press('Enter');

  await expect(page.getByRole('dialog', { name: 'Pairwise match details' })).toBeVisible();
  await expect(firstOgMatch).toHaveClass(/\bgbdraw-match-selected\b/);
  await expect(secondOgMatch).toHaveClass(/\bgbdraw-match-selected\b/);
  await expect(otherOgMatch).not.toHaveClass(/\bgbdraw-match-selected\b/);
  await expect(pairwiseMatch).not.toHaveClass(/\bgbdraw-match-selected\b/);
  await expect.poll(() => firstOgMatch.evaluate((element) => getComputedStyle(element).outlineStyle))
    .toBe('none');

  await page.getByRole('button', { name: 'Close match popup' }).click();
  await expect(page.getByRole('dialog', { name: 'Pairwise match details' })).toHaveCount(0);
  await expect(firstOgMatch).not.toHaveClass(/\bgbdraw-match-selected\b/);
  await expect(secondOgMatch).not.toHaveClass(/\bgbdraw-match-selected\b/);

  await pairwiseMatch.press('Enter');
  await expect(pairwiseMatch).toHaveClass(/\bgbdraw-match-selected\b/);
  await expect(firstOgMatch).not.toHaveClass(/\bgbdraw-match-selected\b/);
  await expect(secondOgMatch).not.toHaveClass(/\bgbdraw-match-selected\b/);
});

test('Linear record rows and N-to-M comparison batches remain keyed by sequence uid', async ({ page }) => {
  await page.goto('/gbdraw/web/index.html', { waitUntil: 'domcontentloaded' });
  await page.waitForFunction(() => window.__GBDRAW_APP__);

  const setup = await page.evaluate(() => {
    const app = window.__GBDRAW_APP__;
    app.mode = 'linear';
    app.addLinearSeq();
    app.addLinearSeq();
    app.addLinearSeq();
    app.files.linearCanonicalComparisons = [{ id: 'stale' }];
    app.losatCacheInfo = [
      { edgeKey: 'stale->edge', key: 'linear-cache' },
      { key: 'circular-cache' }
    ];
    app.setLinearRecordLayoutEnabled(true);
    app.setLinearRecordRow(app.linearSeqs[0].uid, 1);
    app.setLinearRecordRow(app.linearSeqs[1].uid, 1);
    app.setLinearRecordRow(app.linearSeqs[2].uid, 2);
    app.setLinearRecordRow(app.linearSeqs[3].uid, 2);
    app.setLinearComparisonGlobalAction('losat');
    return {
      uids: app.linearSeqs.map((item) => item.uid),
      tokens: app.linearLayoutTokens,
      comparisonMode: app.linearComparisonPlan.mode,
      uniqueUids: new Set(app.linearSeqs.map((item) => item.uid)).size,
      canonicalComparisonCount: app.files.linearCanonicalComparisons.length,
      cacheInfo: app.losatCacheInfo.map((entry) => entry.key)
    };
  });

  const [uidA, uidB, uidC, uidD] = setup.uids;
  const zippedEdgeKeys = [`${uidA}->${uidC}`, `${uidB}->${uidD}`];
  const crossProductEdgeKeys = [
    `${uidA}->${uidC}`, `${uidA}->${uidD}`,
    `${uidB}->${uidC}`, `${uidB}->${uidD}`
  ];
  const boundary = linearComparisonBoundary(page, 1, 2);

  expect(setup.tokens).toEqual(['#1@1', '#2@1', '#3@2', '#4@2']);
  expect(setup.comparisonMode).toBe('adjacent');
  expect(setup.uniqueUids).toBe(4);
  expect(setup.canonicalComparisonCount).toBe(0);
  expect(setup.cacheInfo).toEqual(['circular-cache']);
  const advancedComparison = page.locator(
    'details[data-linear-comparison-disclosure="advanced"]'
  );
  await advancedComparison.locator(':scope > summary').press('Enter');
  await expect(page.getByRole('spinbutton', {
    name: /^Linear record row for sequence \d+$/
  })).toHaveCount(4);
  await expect(page.locator('[data-linear-display-row]')).toHaveCount(2);
  await expect(page.locator('[data-linear-comparison-boundary]')).toHaveCount(1);
  await expectExactEdgeKeys(boundary, zippedEdgeKeys);

  await openLinearSelectedPairs(page);
  const editedPair = linearComparisonPair(page, `${uidA}->${uidC}`);
  await editedPair.getByRole('radio', { name: 'Upload BLAST TSV' }).check();
  const materialized = await page.evaluate(() => ({
    mode: window.__GBDRAW_APP__.linearComparisonPlan.mode,
    edges: window.__GBDRAW_APP__.linearComparisonPlan.edges
      .filter((edge) => edge.included)
      .map((edge) => [
        `${edge.queryUid}->${edge.subjectUid}`,
        edge.source
      ])
  }));
  expect(materialized.mode).toBe('selected');
  expect(Object.fromEntries(materialized.edges)).toEqual({
    [`${uidA}->${uidC}`]: 'upload',
    [`${uidB}->${uidD}`]: 'losat'
  });

  await page.getByRole('button', {
    name: 'All adjacent-row pairs (cross-product)', exact: true
  }).click();
  await expectExactEdgeKeys(boundary, crossProductEdgeKeys);
  const state = await page.evaluate(() => ({
    comparisonMode: window.__GBDRAW_APP__.linearComparisonPlan.mode,
    comparisonCount: window.__GBDRAW_APP__.linearComparisonPlan.edges
      .filter((item) => item.included).length,
    endpoints: window.__GBDRAW_APP__.linearComparisonPlan.edges
      .filter((item) => item.included)
      .map((item) => [item.queryUid, item.subjectUid])
  }));
  expect(state.comparisonMode).toBe('selected');
  expect(state.comparisonCount).toBe(4);
  expect(new Set(state.endpoints.map((pair) => pair.join('->'))).size).toBe(4);

  const globalComparisonControls = page.getByRole('group', {
    name: 'Set all adjacent comparisons'
  });
  await expect(globalComparisonControls.getByRole('button', {
    name: 'Set no comparison'
  })).toBeVisible();
  await expect(globalComparisonControls.getByRole('button', {
    name: 'Run LOSAT for all adjacent pairs'
  })).toBeVisible();
  await expect(globalComparisonControls.getByRole('button', {
    name: 'Use uploaded BLAST TSV for all adjacent pairs'
  })).toBeVisible();
  await globalComparisonControls.getByRole('button', {
    name: 'Set no comparison'
  }).click();
  const optedOut = await page.evaluate(() => ({
    mode: window.__GBDRAW_APP__.linearComparisonPlan.mode,
    retainedDrafts: window.__GBDRAW_APP__.linearComparisonPlan.edges.length,
    resolvedEdges: window.__GBDRAW_APP__.linearComparisonResolution.edges.length
  }));
  expect(optedOut).toEqual({ mode: 'none', retainedDrafts: 4, resolvedEdges: 0 });
  await expect(page.getByRole('group', { name: 'LOSAT Mode' })).toHaveCount(0);
  await expect(page.getByRole('combobox', { name: 'LOSATP mode' })).toHaveCount(0);
});

test('Linear records precede comparison pairs in DOM and keyboard order at narrow width', async ({ page }) => {
  await page.setViewportSize({ width: 390, height: 844 });
  await page.goto('/gbdraw/web/index.html', { waitUntil: 'domcontentloaded' });
  await page.waitForFunction(() => window.__GBDRAW_APP__);

  const uids = await page.evaluate(() => {
    const app = window.__GBDRAW_APP__;
    app.mode = 'linear';
    while (app.linearSeqs.length < 5) app.addLinearSeq();
    app.linearSeqs.forEach((sequence, index) => {
      sequence.definition = `Timeline record ${index + 1}`;
    });
    app.setLinearComparisonGlobalAction('losat');
    return app.linearSeqs.map((sequence) => sequence.uid);
  });
  const edgeKeys = uids.slice(0, -1).map((uid, index) => `${uid}->${uids[index + 1]}`);
  const expectedTimelineTokens = [
    ...uids.map((uid) => `record:${uid}`),
    ...uids.slice(0, -1).map((_, index) => `boundary:${index + 1}->${index + 2}`)
  ];

  const rows = page.locator('[data-linear-display-row]');
  const boundaries = page.locator('[data-linear-comparison-boundary]');
  await expect(rows).toHaveCount(5);
  await expect(boundaries).toHaveCount(4);
  expect(await rows.evaluateAll((elements) => (
    elements.map((element) => element.dataset.linearDisplayRow)
  ))).toEqual(['1', '2', '3', '4', '5']);
  expect(await boundaries.evaluateAll((elements) => (
    elements.map((element) => element.dataset.linearComparisonBoundary)
  ))).toEqual(['1->2', '2->3', '3->4', '4->5']);
  const timelineTokens = await page.locator(
    '[data-linear-record-card], [data-linear-comparison-boundary]'
  ).evaluateAll((elements) => elements.map((element) => (
    element.dataset.linearRecordCard
      ? `record:${element.dataset.linearRecordCard}`
      : `boundary:${element.dataset.linearComparisonBoundary}`
  )));
  expect(timelineTokens).toEqual(expectedTimelineTokens);
  for (let index = 0; index < edgeKeys.length; index += 1) {
    await expectExactEdgeKeys(linearComparisonBoundary(page, index + 1, index + 2), [edgeKeys[index]]);
  }

  const recordList = page.locator('[data-linear-record-list]');
  await expect(recordList.locator('[data-linear-comparison-boundary]')).toHaveCount(0);
  await expect(recordList.locator('[data-edge-key]')).toHaveCount(0);
  await expect(linearSelectedPairs(page)).not.toHaveAttribute('open', '');

  await focusLastTabStop(linearRecordCard(page, uids[0]));
  await page.keyboard.press('Tab');
  expect(await page.evaluate(() => ({
    edgeKey: document.activeElement?.closest('[data-edge-key]')?.dataset.edgeKey || '',
    recordUid: document.activeElement?.closest('[data-linear-record-card]')
      ?.dataset.linearRecordCard || ''
  }))).toMatchObject({ edgeKey: '' });

  await page.getByRole('button', { name: 'Add sequence' }).last().focus();
  const focusTrail = [];
  for (let step = 0; step < 5; step += 1) {
    await page.keyboard.press('Tab');
    focusTrail.push(await page.evaluate(() => ({
      edgeKey: document.activeElement?.closest('[data-edge-key]')?.dataset.edgeKey || '',
      command: document.activeElement?.getAttribute('aria-label') || ''
    })));
    if (focusTrail.at(-1).command === 'Set no comparison') break;
  }
  expect(focusTrail.some((entry) => entry.edgeKey)).toBe(false);
  expect(focusTrail.at(-1).command).toBe('Set no comparison');

  const overflow = await page.locator('.settings-pane').evaluate((settingsPane) => {
    const comparison = settingsPane.querySelector('[data-linear-comparison-card]');
    return {
      paneClientWidth: settingsPane.clientWidth,
      paneScrollWidth: settingsPane.scrollWidth,
      comparisonClientWidth: comparison?.clientWidth || 0,
      comparisonScrollWidth: comparison?.scrollWidth || 0
    };
  });
  expect(overflow.paneScrollWidth).toBeLessThanOrEqual(overflow.paneClientWidth + 1);
  expect(overflow.comparisonScrollWidth).toBeLessThanOrEqual(
    overflow.comparisonClientWidth + 1
  );
});

test('Linear region controls do not overlap at supported sidebar widths', async ({ page }) => {
  await page.goto('/gbdraw/web/index.html', { waitUntil: 'domcontentloaded' });
  await page.waitForFunction(() => window.__GBDRAW_APP__);
  await page.getByRole('button', { name: 'Linear', exact: true }).click();

  const recordOptions = page.locator('details[data-linear-record-options]').first();
  await expect(recordOptions).toHaveCount(1);
  if (!await recordOptions.evaluate((element) => element.open)) {
    await recordOptions.locator(':scope > summary').press('Enter');
  }
  await expect(recordOptions).toHaveAttribute('open', '');

  const recordSelector = page.getByRole('combobox', {
    name: 'Record selector for sequence 1', exact: true
  });
  const reverseComplement = page.getByRole('checkbox', {
    name: 'Reverse complement for sequence 1', exact: true
  });
  await expect(recordSelector).toBeVisible();
  await expect(reverseComplement).toBeVisible();

  for (const width of [320, 240]) {
    await page.evaluate((sidebarWidth) => {
      window.__GBDRAW_APP__.sidebarWidth = sidebarWidth;
    }, width);
    await expect.poll(
      () => page.evaluate(() => window.__GBDRAW_APP__.sidebarWidth)
    ).toBe(width);

    const [selectorGeometry, reverseGeometry] = await Promise.all([
      recordSelector.evaluate((select) => {
        const cell = select.parentElement;
        const grid = cell?.parentElement;
        const bounds = (element) => {
          const rect = element?.getBoundingClientRect();
          return rect ? {
            x: rect.x, y: rect.y, width: rect.width, height: rect.height
          } : null;
        };
        return { control: bounds(select), cell: bounds(cell), grid: bounds(grid) };
      }),
      reverseComplement.evaluate((checkbox) => {
        let label = checkbox.parentElement;
        while (label && String(label.tagName || '').toLowerCase() !== 'label') {
          label = label.parentElement;
        }
        const rect = label?.getBoundingClientRect();
        return {
          label: rect ? {
            x: rect.x, y: rect.y, width: rect.width, height: rect.height
          } : null,
          clientWidth: label?.clientWidth || 0,
          scrollWidth: label?.scrollWidth || 0
        };
      })
    ]);
    const selectorBox = selectorGeometry.control;
    const selectorCellBox = selectorGeometry.cell;
    const gridBox = selectorGeometry.grid;
    const reverseBox = reverseGeometry.label;
    expect(selectorBox).not.toBeNull();
    expect(selectorCellBox).not.toBeNull();
    expect(gridBox).not.toBeNull();
    expect(reverseBox).not.toBeNull();
    const intersects = (
      selectorBox.x < reverseBox.x + reverseBox.width
      && selectorBox.x + selectorBox.width > reverseBox.x
      && selectorBox.y < reverseBox.y + reverseBox.height
      && selectorBox.y + selectorBox.height > reverseBox.y
    );
    expect(intersects, `Region controls overlap at ${width}px`).toBe(false);
    expect(selectorBox.x).toBeGreaterThanOrEqual(selectorCellBox.x - 1);
    expect(selectorBox.x + selectorBox.width).toBeLessThanOrEqual(
      selectorCellBox.x + selectorCellBox.width + 1
    );
    expect(reverseGeometry.scrollWidth).toBeLessThanOrEqual(
      reverseGeometry.clientWidth + 1
    );
    expect(selectorCellBox.x).toBeGreaterThanOrEqual(gridBox.x - 1);
    expect(reverseBox.x + reverseBox.width).toBeLessThanOrEqual(
      gridBox.x + gridBox.width + 1
    );
  }
});

test('Selected pairs focuses Add and repairs an unplaced draft in its boundary', async ({ page }) => {
  await page.goto('/gbdraw/web/index.html', { waitUntil: 'domcontentloaded' });
  await page.waitForFunction(() => window.__GBDRAW_APP__);

  const uids = await page.evaluate(() => {
    const app = window.__GBDRAW_APP__;
    app.mode = 'linear';
    app.addLinearSeq();
    app.addLinearSeq();
    app.linearSeqs.forEach((sequence, index) => {
      sequence.definition = `Advanced record ${index + 1}`;
    });
    app.setLinearComparisonGlobalAction('losat');
    app.clearSelectedLinearComparisons();
    return app.linearSeqs.map((sequence) => sequence.uid);
  });
  const [uidA, uidB, uidC] = uids;
  await openLinearSelectedPairs(page);
  await page.getByRole('button', { name: 'Add', exact: true }).click();

  const addedPair = linearComparisonPair(page, `${uidA}->${uidB}`);
  await expect(addedPair).toBeVisible();
  expect(await addedPair.evaluate((element) => element.contains(document.activeElement))).toBe(true);
  expect(await page.evaluate(() => ({
    mode: window.__GBDRAW_APP__.linearComparisonPlan.mode,
    included: window.__GBDRAW_APP__.linearComparisonPlan.edges.filter((edge) => edge.included).length
  }))).toEqual({ mode: 'selected', included: 1 });

  await page.evaluate(({ queryUid, subjectUid }) => {
    const app = window.__GBDRAW_APP__;
    app.linearComparisonPlan.mode = 'selected';
    app.linearComparisonPlan.edges.splice(0, app.linearComparisonPlan.edges.length, {
      id: 'non-adjacent-browser-draft',
      queryUid,
      subjectUid,
      included: true,
      fileActive: false,
      losatFilenameActive: false,
      source: 'losat',
      file: null,
      losatFilename: ''
    });
  }, { queryUid: uidA, subjectUid: uidC });

  await expect(linearComparisonPair(page, `${uidA}->${uidC}`)).toHaveCount(0);
  const unplacedDraft = page.locator(
    '[data-linear-unplaced-draft="non-adjacent-browser-draft"]'
  );
  await expect(unplacedDraft).toBeVisible();
  await expect(unplacedDraft).toContainText(/adjacent (display )?rows/i);
  await unplacedDraft.getByRole('combobox', {
    name: 'To record for unplaced comparison non-adjacent-browser-draft'
  }).selectOption(uidB);

  await expect(unplacedDraft).toHaveCount(0);
  await expectExactEdgeKeys(linearComparisonBoundary(page, 1, 2), [`${uidA}->${uidB}`]);
  await expect(addedPair.getByRole('radio', { name: 'Run LOSAT' })).toBeChecked();
  expect(await page.evaluate(() => ({
    valid: window.__GBDRAW_APP__.linearComparisonResolution.valid,
    edgeKeys: window.__GBDRAW_APP__.linearComparisonResolution.edges.map((edge) => edge.edgeKey)
  }))).toEqual({ valid: true, edgeKeys: [`${uidA}->${uidB}`] });
});

test('Comparison card actions target the active owner of duplicate directional drafts', async ({ page }) => {
  await page.goto('/gbdraw/web/index.html', { waitUntil: 'domcontentloaded' });
  await page.waitForFunction(() => window.__GBDRAW_APP__);

  const [uidA, uidB] = await page.evaluate(() => {
    const app = window.__GBDRAW_APP__;
    app.mode = 'linear';
    app.addLinearSeq();
    const [first, second] = app.linearSeqs;
    app.linearComparisonPlan.mode = 'selected';
    app.linearComparisonPlan.defaultSource = 'losat';
    app.linearComparisonPlan.edges.splice(0, app.linearComparisonPlan.edges.length,
      {
        id: 'inactive-first-duplicate', queryUid: first.uid, subjectUid: second.uid,
        included: false, fileActive: false, losatFilenameActive: false,
        source: 'losat', file: null, losatFilename: 'inactive-first.tsv'
      },
      {
        id: 'active-second-duplicate', queryUid: first.uid, subjectUid: second.uid,
        included: true, fileActive: false, losatFilenameActive: true,
        source: 'losat', file: null, losatFilename: 'active-second.tsv'
      }
    );
    return [first.uid, second.uid];
  });
  const edgeKey = `${uidA}->${uidB}`;
  await openLinearSelectedPairs(page);
  const pair = linearComparisonPair(page, edgeKey);
  await expect(pair).toBeVisible();
  expect(await page.evaluate((key) => {
    const app = window.__GBDRAW_APP__;
    const owner = app.linearComparisonTimeline.rows
      .flatMap((row) => row.boundaryAfter?.pairs || [])
      .find((entry) => entry.edgeKey === key);
    return [
      owner?.edgeId,
      owner?.draft?.id,
      owner?.resolved?.id,
      app.linearComparisonTimeline.unplacedDrafts.map((entry) => entry.draft.id)
    ];
  }, edgeKey)).toEqual([
    'active-second-duplicate',
    'active-second-duplicate',
    'active-second-duplicate',
    ['inactive-first-duplicate']
  ]);

  await pair.getByRole('radio', { name: 'Upload BLAST TSV' }).check();
  expect(await page.evaluate(() => window.__GBDRAW_APP__.linearComparisonPlan.edges.map((edge) => [
    edge.id, edge.included, edge.source, edge.fileActive, edge.file?.name || '', edge.losatFilename
  ]))).toEqual([
    ['inactive-first-duplicate', false, 'losat', false, '', 'inactive-first.tsv'],
    ['active-second-duplicate', true, 'upload', false, '', 'active-second.tsv']
  ]);

  await pair.locator('input[type="file"][aria-label="BLAST TSV for #1 to #2"]').setInputFiles({
    name: 'active-second-upload.tsv',
    mimeType: 'text/tab-separated-values',
    buffer: Buffer.from('A\tB\t100\t30\t0\t0\t1\t30\t1\t30\t1e-20\t90\n')
  });
  expect(await page.evaluate(() => window.__GBDRAW_APP__.linearComparisonPlan.edges.map((edge) => [
    edge.id, edge.included, edge.source, edge.fileActive, edge.file?.name || '', edge.losatFilename
  ]))).toEqual([
    ['inactive-first-duplicate', false, 'losat', false, '', 'inactive-first.tsv'],
    ['active-second-duplicate', true, 'upload', true, 'active-second-upload.tsv', 'active-second.tsv']
  ]);
});

test('Normalize Record Lengths rejects a shared Linear row and remains recoverable', async ({ page }) => {
  test.setTimeout(300000);
  await installDiagramRequestObserver(page);
  await page.goto('/gbdraw/web/index.html', { waitUntil: 'domcontentloaded' });
  await page.waitForFunction(() => window.__GBDRAW_APP__);
  const sharedRowUids = await page.evaluate((records) => {
    const app = window.__GBDRAW_APP__;
    app.mode = 'linear';
    app.lInputType = 'gb';
    app.addLinearSeq();
    records.forEach((content, index) => app.setLinearSeqPrimaryFile(index, 'gb', new File(
      [content], `normalize-record-${index + 1}.gbk`, {
        type: 'text/plain',
        lastModified: index + 1
      }
    )));
    Object.assign(app.form, {
      legend: 'none',
      show_gc: false,
      show_skew: false,
      show_depth: false,
      show_labels_linear: 'none',
      normalize_length: true
    });
    app.setLinearComparisonGlobalAction('none');
    app.setLinearRecordLayoutEnabled(true);
    app.setLinearRecordRow(app.linearSeqs[0].uid, 1);
    app.setLinearRecordRow(app.linearSeqs[1].uid, 1);
    return app.linearSeqs.map((sequence) => sequence.uid);
  }, [
    makeComparisonGenbank('NormalizeRecA', 'atg'),
    makeComparisonGenbank('NormalizeRecB', 'gct')
  ]);

  await expect(page.locator('[data-linear-display-row="1"]')).toBeVisible();
  await expect(page.locator('[data-linear-display-row]')).toHaveCount(1);
  await expect(page.locator('[data-linear-comparison-boundary]')).toHaveCount(0);
  await expect(linearRecordCard(page, sharedRowUids[0])).toBeVisible();
  await expect(linearRecordCard(page, sharedRowUids[1])).toBeVisible();

  const invalid = await page.evaluate(async () => {
    const app = window.__GBDRAW_APP__;
    const result = await app.runAnalysis();
    return {
      result,
      errorSummary: String(app.errorLog?.summary || ''),
      errorDetails: Array.isArray(app.errorLog?.details)
        ? app.errorLog.details.map((detail) => String(detail))
        : [],
      dispatchedRequests: window.__GBDRAW_DIAGRAM_RUNS__.length,
      normalizeLength: app.form.normalize_length
    };
  });
  expect(invalid.result).toEqual({ status: 'error' });
  expect([invalid.errorSummary, ...invalid.errorDetails].join(' ')).toMatch(
    /Normalize Record Lengths.*same Linear row/i
  );
  expect(invalid.dispatchedRequests).toBe(0);
  expect(invalid.normalizeLength).toBe(true);

  const recovered = await page.evaluate(async () => {
    const app = window.__GBDRAW_APP__;
    app.form.normalize_length = false;
    const result = await app.runAnalysis();
    return {
      result,
      errorSummary: String(app.errorLog?.summary || ''),
      dispatchedRequests: window.__GBDRAW_DIAGRAM_RUNS__.length,
      normalizeLength: app.form.normalize_length
    };
  });
  expect(recovered).toEqual({
    result: { status: 'ok' },
    errorSummary: '',
    dispatchedRequests: 1,
    normalizeLength: false
  });
});

test('No comparison completes a real render without touching dormant comparison work', async ({ page }) => {
  test.setTimeout(300000);
  await installDiagramRequestObserver(page);
  await page.addInitScript(() => {
    window.__GBDRAW_LOSAT_EXECUTOR_CALLS__ = 0;
    window.__GBDRAW_LOSAT_EXECUTOR__ = async () => {
      window.__GBDRAW_LOSAT_EXECUTOR_CALLS__ += 1;
      throw new Error('No comparison must not execute LOSAT.');
    };
  });
  await page.goto('/gbdraw/web/index.html', { waitUntil: 'domcontentloaded' });
  await page.waitForFunction(() => window.__GBDRAW_APP__);
  const dormantEdgeKey = await page.evaluate(async (records) => {
    const app = window.__GBDRAW_APP__;
    const { state } = await import('./js/state.js');
    app.mode = 'linear';
    app.lInputType = 'gb';
    app.addLinearSeq();
    app.addLinearSeq();
    records.forEach((content, index) => app.setLinearSeqPrimaryFile(index, 'gb', new File(
      [content], `none-record-${index + 1}.gbk`, { type: 'text/plain', lastModified: index + 1 }
    )));
    Object.assign(app.form, {
      legend: 'bottom', show_gc: false, show_skew: false,
      show_depth: false, show_labels_linear: 'none'
    });
    const dormant = new File([
      'DormantA\tDormantB\t100\t30\t0\t0\t1\t30\t1\t30\t1e-20\t100\n'
    ], 'dormant-comparison.tsv', { type: 'text/tab-separated-values', lastModified: 99 });
    const nativeArrayBuffer = dormant.arrayBuffer.bind(dormant);
    window.__GBDRAW_DORMANT_READS__ = 0;
    dormant.arrayBuffer = () => {
      window.__GBDRAW_DORMANT_READS__ += 1;
      return nativeArrayBuffer();
    };
    const [first, second] = app.linearSeqs;
    app.linearComparisonPlan.edges.splice(0, app.linearComparisonPlan.edges.length, {
      id: 'dormant-none-edge', queryUid: first.uid, subjectUid: second.uid,
      included: true, fileActive: true, losatFilenameActive: true,
      source: 'upload', file: dormant, losatFilename: 'dormant-losat-name.tsv'
    });
    app.setLinearComparisonGlobalAction('none');
    const rawCache = new Map([['dormant-cache', { text: 'must remain unread' }]]);
    const nativeGet = rawCache.get.bind(rawCache);
    window.__GBDRAW_CACHE_LOOKUPS__ = 0;
    rawCache.get = (key) => {
      window.__GBDRAW_CACHE_LOOKUPS__ += 1;
      return nativeGet(key);
    };
    state.losatCache.value = window.Vue.markRaw(rawCache);
    return `${first.uid}->${second.uid}`;
  }, [
    makeComparisonGenbank('NoneRecA', 'atg'),
    makeComparisonGenbank('NoneRecB', 'gct'),
    makeComparisonGenbank('NoneRecC', 'tta')
  ]);

  const outcome = await page.evaluate(async () => {
    const app = window.__GBDRAW_APP__;
    const { state } = await import('./js/state.js');
    const result = await app.runAnalysis();
    const request = window.__GBDRAW_DIAGRAM_RUNS__.at(-1);
    const svg = new DOMParser().parseFromString(app.results[0].content, 'image/svg+xml');
    return {
      result,
      error: app.errorLog,
      request: [request.records.length, request.comparisons.length],
      planMode: app.linearComparisonPlan.mode,
      timelineSource: app.linearComparisonTimeline.rows[0].boundaryAfter.pairs[0].source,
      resolution: [Object.isFrozen(app.linearComparisonResolution), app.linearComparisonResolution.edges.length],
      draft: app.linearComparisonPlan.edges.map((edge) => [
        edge.id, edge.included, edge.fileActive, edge.losatFilenameActive,
        edge.file?.name, edge.losatFilename
      ]),
      skipped: [
        window.__GBDRAW_DORMANT_READS__, window.__GBDRAW_CACHE_LOOKUPS__,
        window.__GBDRAW_LOSAT_EXECUTOR_CALLS__
      ],
      cacheSize: state.losatCache.value.size,
      losatRuntimeUrls: [
        ...window.__GBDRAW_RUNTIME_URLS__.fetches,
        ...window.__GBDRAW_RUNTIME_URLS__.workers
      ].filter((url) => url.includes('losat')),
      comparisonSvgNodes: svg.querySelectorAll(
        '[data-query-row], [data-subject-row], [data-gbdraw-role="comparison-legend"], #pairwise_legend'
      ).length
    };
  });
  expect(outcome).toEqual({
    result: { status: 'ok' }, error: null, request: [3, 0], planMode: 'none', timelineSource: 'none',
    resolution: [true, 0],
    draft: [['dormant-none-edge', true, true, true, 'dormant-comparison.tsv', 'dormant-losat-name.tsv']],
    skipped: [0, 0, 0], cacheSize: 1, losatRuntimeUrls: [], comparisonSvgNodes: 0
  });
  await expect(page.getByText('Raw LOSAT results', { exact: true })).not.toBeVisible();

  await openLinearSelectedPairs(page);
  const dormantPair = linearComparisonPair(page, dormantEdgeKey);
  await expect(dormantPair).toBeVisible();
  const dormantUpload = dormantPair.locator(
    'input[type="file"][aria-label="BLAST TSV for #1 to #2"]'
  );
  await expect(dormantUpload).toHaveCount(0);
  await expect(dormantPair.getByRole('textbox', {
    name: 'Raw LOSAT filename for #1 to #2'
  })).toHaveCount(0);
  await expect(dormantPair.getByRole('button', {
    name: 'Save Raw LOSAT TSV for #1 to #2'
  })).toHaveCount(0);
  await expect(dormantPair).toContainText('Retained, inactive BLAST TSV: dormant-comparison.tsv');

  await openLinearAdvancedComparison(page);
  const dormantRawResult = page.locator(`[data-linear-raw-result="${dormantEdgeKey}"]`);
  await expect(dormantRawResult).toContainText(
    'Retained, inactive filename: dormant-losat-name.tsv'
  );

  await dormantPair.getByRole('button', { name: 'Reuse BLAST TSV for #1 to #2' }).click();
  await expect(dormantUpload).toHaveCount(1);
  await expect(dormantPair).toContainText('dormant-comparison.tsv');
  expect(await page.evaluate(() => {
    const edge = window.__GBDRAW_APP__.linearComparisonPlan.edges[0];
    return [
      window.__GBDRAW_APP__.linearComparisonPlan.mode,
      edge.included,
      edge.source,
      edge.fileActive,
      edge.file?.name
    ];
  })).toEqual(['selected', true, 'upload', true, 'dormant-comparison.tsv']);

  await dormantPair.getByRole('radio', { name: 'No comparison' }).check();
  await expect(dormantUpload).toHaveCount(0);
  await expect(dormantPair).toContainText('Retained, inactive BLAST TSV: dormant-comparison.tsv');
  expect(await page.evaluate(() => {
    const edge = window.__GBDRAW_APP__.linearComparisonPlan.edges[0];
    return [edge.included, edge.fileActive, edge.file?.name, edge.losatFilename];
  })).toEqual([false, true, 'dormant-comparison.tsv', 'dormant-losat-name.tsv']);

  await dormantRawResult.getByRole('button', {
    name: 'Reuse Raw LOSAT filename for #1 to #2'
  }).click();
  await expect(dormantRawResult.getByRole('textbox', {
    name: 'Raw LOSAT filename for #1 to #2'
  })).toHaveValue('dormant-losat-name.tsv');
  await expect(dormantRawResult.getByRole('button', {
    name: 'Save Raw LOSAT TSV for #1 to #2'
  })).toBeDisabled();
});

test('Automatic Linear renders every record from one GenBank source and survives session reload', async ({ page }) => {
  test.setTimeout(300000);
  await installDiagramRequestObserver(page);
  await page.goto('/gbdraw/web/index.html', { waitUntil: 'domcontentloaded' });
  await page.waitForFunction(() => window.__GBDRAW_APP__);
  await page.waitForFunction(
    () => Object.keys(window.__GBDRAW_APP__?.paletteDefinitions || {}).length > 0,
    null,
    { timeout: 180000 }
  );

  const sourceName = 'automatic-multi-record.gbk';
  await page.evaluate(({ content, name }) => {
    const app = window.__GBDRAW_APP__;
    app.mode = 'linear';
    app.lInputType = 'gb';
    const source = new File(
      [content], name, { type: 'text/plain', lastModified: 1 }
    );
    const nativeArrayBuffer = source.arrayBuffer.bind(source);
    window.__GBDRAW_SOURCE_READS__ = 0;
    source.arrayBuffer = () => {
      window.__GBDRAW_SOURCE_READS__ += 1;
      return nativeArrayBuffer();
    };
    app.setLinearSeqPrimaryFile(0, 'gb', source);
    Object.assign(app.form, {
      legend: 'none',
      show_gc: false,
      show_skew: false,
      show_depth: false,
      show_labels_linear: 'none'
    });
    app.setLinearComparisonGlobalAction('none');
    app.sessionTitle = 'automatic-multi-record';
  }, {
    content: makeComparisonGenbank('AutomaticA', 'atg') +
      makeComparisonGenbank('AutomaticB', 'gct'),
    name: sourceName
  });
  await expect.poll(() => page.evaluate(() => (
    window.__GBDRAW_APP__.linearRecordOptions(window.__GBDRAW_APP__.linearSeqs[0]).length
  ))).toBe(3);

  await page.evaluate(() => {
    window.__GBDRAW_AUTOMATIC_RUN__ = { done: false, result: null, error: '' };
    window.__GBDRAW_APP__.runAnalysis().then((result) => {
      Object.assign(window.__GBDRAW_AUTOMATIC_RUN__, { done: true, result });
    }).catch((error) => {
      Object.assign(window.__GBDRAW_AUTOMATIC_RUN__, {
        done: true,
        error: error?.message || String(error)
      });
    });
  });
  await page.waitForFunction(() => window.__GBDRAW_AUTOMATIC_RUN__?.done === true);
  const generated = await page.evaluate(() => {
    const app = window.__GBDRAW_APP__;
    const run = window.__GBDRAW_AUTOMATIC_RUN__;
    if (run.error) throw new Error(run.error);
    const request = window.__GBDRAW_DIAGRAM_RUNS__.at(-1);
    return {
      result: run.result,
      error: app.errorLog,
      sourceReads: window.__GBDRAW_SOURCE_READS__,
      cardCount: app.linearSeqs.length,
      selector: app.linearSeqs[0].region_record_id,
      grouping: request.grouping,
      selectors: request.records.map((record) => record.selector),
      sharedResourceCount: new Set(
        request.records.map((record) => record.source.resourceId)
      ).size
    };
  });
  expect(generated).toEqual({
    result: { status: 'ok' },
    error: null,
    sourceReads: 1,
    cardCount: 1,
    selector: '',
    grouping: 'single',
    selectors: [
      { kind: 'recordIndex', index: 0 },
      { kind: 'recordIndex', index: 1 }
    ],
    sharedResourceCount: 1
  });

  const sessionDownloadPromise = page.waitForEvent('download', { timeout: 120000 });
  expect((await page.evaluate(() => window.__GBDRAW_APP__.saveSessionWithTitle())).status)
    .toBe('saved');
  const sessionPath = await (await sessionDownloadPromise).path();
  const session = JSON.parse(gunzipSync(readFileSync(sessionPath)).toString('utf8'));
  expect(session.webFiles.bindings.linearSeqs).toHaveLength(1);
  expect(session.webFiles.bindings.linearSeqs[0].region_record_id).toBe('');
  expect(session.renderRequest.records.map((record) => record.selector)).toEqual([
    { kind: 'recordIndex', index: 0 },
    { kind: 'recordIndex', index: 1 }
  ]);
  expect(new Set(
    session.renderRequest.records.map((record) => record.source.resourceId)
  ).size).toBe(1);

  await page.addInitScript((recordSourceName) => {
    const nativeArrayBuffer = File.prototype.arrayBuffer;
    let releaseRecordRead;
    const recordReadGate = new Promise((resolve) => { releaseRecordRead = resolve; });
    window.__GBDRAW_RECORD_TEXT_READS__ = 0;
    window.__GBDRAW_RELEASE_RECORD_TEXT__ = releaseRecordRead;
    File.prototype.arrayBuffer = async function (...args) {
      if (this.name !== recordSourceName) return nativeArrayBuffer.apply(this, args);
      const { state } = await import('./js/state.js');
      if (state.semanticFileWatchersSuppressed.value) {
        return nativeArrayBuffer.apply(this, args);
      }
      window.__GBDRAW_RECORD_TEXT_READS__ += 1;
      await recordReadGate;
      return nativeArrayBuffer.apply(this, args);
    };
  }, sourceName);
  await page.reload({ waitUntil: 'domcontentloaded' });
  await page.waitForFunction(() => window.__GBDRAW_APP__);
  await page.waitForFunction(
    () => Object.keys(window.__GBDRAW_APP__?.paletteDefinitions || {}).length > 0,
    null,
    { timeout: 180000 }
  );
  const dialogPromise = page.waitForEvent('dialog', { timeout: 120000 });
  await page.locator('input[accept^=".json,"]').first().setInputFiles(sessionPath);
  const dialog = await dialogPromise;
  expect(dialog.message()).toBe('Session loaded successfully!');
  await dialog.accept();
  await page.waitForFunction(() => window.__GBDRAW_RECORD_TEXT_READS__ === 1);

  await page.evaluate(async () => {
    const { state } = await import('./js/state.js');
    state.labelReflowForceRequestReason.value = 'multi-record-session-reflow';
    state.labelReflowForceRequestSeq.value += 1;
    await window.Vue.nextTick();
  });
  await expect.poll(() => page.evaluate(() => window.__GBDRAW_RECORD_TEXT_READS__)).toBe(1);
  expect(await page.evaluate(async () => {
    const { state } = await import('./js/state.js');
    return [window.__GBDRAW_DIAGRAM_RUNS__.length, state.labelReflowLastError.value];
  })).toEqual([0, null]);
  await page.evaluate(() => window.__GBDRAW_RELEASE_RECORD_TEXT__());
  await page.waitForFunction(() => window.__GBDRAW_DIAGRAM_RUNS__.length === 1);

  const restored = await page.evaluate(async () => {
    const app = window.__GBDRAW_APP__;
    const { state } = await import('./js/state.js');
    const request = window.__GBDRAW_DIAGRAM_RUNS__.at(-1);
    return {
      error: state.labelReflowLastError.value,
      cardCount: app.linearSeqs.length,
      selector: app.linearSeqs[0].region_record_id,
      optionCount: app.linearRecordOptions(app.linearSeqs[0]).length,
      selectors: request.records.map((record) => record.selector),
      sharedResource: new Set(
        request.records.map((record) => record.source.resourceId)
      ).size
    };
  });
  expect(restored).toEqual({
    error: null,
    cardCount: 1,
    selector: '',
    optionCount: 3,
    selectors: [
      { kind: 'recordIndex', index: 0 },
      { kind: 'recordIndex', index: 1 }
    ],
    sharedResource: 1
  });
});

test('Sparse upload and mixed selected renders keep snapshots and raw cache identity stable', async ({ page }) => {
  test.setTimeout(420000);
  await installDiagramRequestObserver(page);
  await page.addInitScript(() => {
    window.__GBDRAW_LOSAT_EXECUTOR_CALLS__ = 0;
    window.__GBDRAW_LOSAT_EXECUTOR_JOBS__ = [];
    window.__GBDRAW_LOSAT_EXECUTOR__ = async (jobs) => {
      window.__GBDRAW_LOSAT_EXECUTOR_CALLS__ += 1;
      window.__GBDRAW_LOSAT_EXECUTOR_JOBS__.push(jobs.map((job) => [
        job.edgeKey, job.ordinal, job.queryIndex, job.subjectIndex, job.cacheKey
      ]));
      if (window.__GBDRAW_LOSAT_EXECUTOR_CALLS__ > 1) {
        throw new Error('A content-addressed cache hit must not execute another LOSAT job.');
      }
      await new Promise((resolve) => { window.__GBDRAW_RELEASE_LOSAT__ = resolve; });
      return jobs.map((job) => ({
        cacheKey: job.cacheKey,
        text: 'MixedRecB\tMixedRecC\t100\t60\t0\t0\t1\t60\t1\t60\t1e-30\t150\n'
      }));
    };
  });
  await page.goto('/gbdraw/web/index.html', { waitUntil: 'domcontentloaded' });
  await page.waitForFunction(() => window.__GBDRAW_APP__);
  const [uidA, uidB, uidC] = await page.evaluate((records) => {
    const app = window.__GBDRAW_APP__;
    app.mode = 'linear';
    app.lInputType = 'gb';
    app.addLinearSeq();
    app.addLinearSeq();
    records.forEach((content, index) => app.setLinearSeqPrimaryFile(index, 'gb', new File(
      [content], `mixed-record-${index + 1}.gbk`, { type: 'text/plain', lastModified: index + 10 }
    )));
    Object.assign(app.form, {
      legend: 'bottom', show_gc: false, show_skew: false,
      show_depth: false, show_labels_linear: 'none'
    });
    app.setLinearComparisonGlobalAction('losat');
    app.setLinearComparisonLosatMode('blastp');
    app.setLinearComparisonLosatpMode('collinear');
    app.setLinearComparisonLosatMode('blastn');
    app.losat.executionMode = 'serial';
    const [first, second, third] = app.linearSeqs;
    const blastRow = 'MixedRecA\tMixedRecB\t100\t60\t0\t0\t1\t60\t1\t60\t1e-30\t150\n';
    app.linearComparisonPlan.mode = 'adjacent';
    app.linearComparisonPlan.defaultSource = 'upload';
    app.linearComparisonPlan.edges.splice(0, app.linearComparisonPlan.edges.length,
      {
        id: 'upload-a-b', queryUid: first.uid, subjectUid: second.uid,
        included: false, fileActive: true, losatFilenameActive: false, source: 'upload',
        file: new File([blastRow], 'mixed-a-to-b.tsv'), losatFilename: ''
      },
      {
        id: 'retained-b-c', queryUid: second.uid, subjectUid: third.uid,
        included: false, fileActive: false, losatFilenameActive: false, source: 'upload',
        file: new File([blastRow], 'retained-b-to-c.tsv'), losatFilename: ''
      }
    );
    return [first.uid, second.uid, third.uid];
  }, [
    makeComparisonGenbank('MixedRecA', 'atg'),
    makeComparisonGenbank('MixedRecB', 'gct'),
    makeComparisonGenbank('MixedRecC', 'tta')
  ]);
  const requestPairs = (entries) => entries.map((entry) => [
    entry.kind, entry.queryRecordIndex, entry.subjectRecordIndex
  ]);
  const sparseUpload = await page.evaluate(async () => {
    const app = window.__GBDRAW_APP__;
    const result = await app.runAnalysis();
    return {
      result,
      comparisons: window.__GBDRAW_DIAGRAM_RUNS__.at(-1).comparisons,
      executorCalls: window.__GBDRAW_LOSAT_EXECUTOR_CALLS__
    };
  });
  expect(sparseUpload.result).toEqual({ status: 'ok' });
  expect(requestPairs(sparseUpload.comparisons)).toEqual([['nucleotideBlast', 0, 1]]);
  expect(sparseUpload.executorCalls).toBe(0);

  const independentReuse = await page.evaluate(() => {
    const app = window.__GBDRAW_APP__;
    app.setLinearComparisonLosatFilename('retained-b-c', 'retained-name.tsv');
    app.deactivateLinearComparisonLosatFilename('retained-b-c');
    app.reuseLinearComparisonFile('retained-b-c');
    const afterFile = app.linearComparisonPlan.edges.find((edge) => edge.id === 'retained-b-c');
    const fileOnly = [afterFile.fileActive, afterFile.losatFilenameActive, afterFile.source];
    app.deactivateLinearComparisonFile('retained-b-c');
    app.reuseLinearComparisonLosatFilename('retained-b-c');
    const afterName = app.linearComparisonPlan.edges.find((edge) => edge.id === 'retained-b-c');
    return [fileOnly, [afterName.fileActive, afterName.losatFilenameActive, afterName.source]];
  });
  expect(independentReuse).toEqual([[true, false, 'upload'], [false, true, 'losat']]);

  const selectedResolution = await page.evaluate(async () => {
    const app = window.__GBDRAW_APP__;
    const { buildPairwiseLosatJobSpecs, resolveLinearComparisonPlan } = await import(
      './js/app/linear-comparisons.js'
    );
    const [first, second, third] = app.linearSeqs;
    const upload = app.linearComparisonPlan.edges.find((edge) => edge.id === 'upload-a-b').file;
    app.linearComparisonPlan.mode = 'selected';
    app.linearComparisonPlan.defaultSource = 'losat';
    app.linearComparisonPlan.edges.splice(0, app.linearComparisonPlan.edges.length,
      {
        id: 'selected-upload-a-b', queryUid: first.uid, subjectUid: second.uid,
        included: true, fileActive: true, losatFilenameActive: false,
        source: 'upload', file: upload, losatFilename: ''
      },
      {
        id: 'selected-losat-b-c', queryUid: second.uid, subjectUid: third.uid,
        included: true, fileActive: false, losatFilenameActive: true,
        source: 'losat', file: null, losatFilename: 'selected-b-to-c.tsv'
      },
      {
        id: 'omitted-a-c', queryUid: first.uid, subjectUid: third.uid,
        included: false, fileActive: false, losatFilenameActive: false,
        source: 'upload', file: null, losatFilename: ''
      }
    );
    const programJobs = ['blastn', 'tblastx', 'blastp'].map((program) => {
      const resolution = resolveLinearComparisonPlan({
        plan: app.linearComparisonPlan,
        sequences: app.linearSeqs,
        losatProgram: program,
        blastpMode: 'pairwise'
      });
      return buildPairwiseLosatJobSpecs({ resolution, program, blastpMode: 'pairwise' })
        .map((job) => [job.program, job.edgeKey, job.queryIndex, job.subjectIndex]);
    });
    return [app.linearComparisonResolution.valid, Object.isFrozen(app.linearComparisonResolution),
      app.linearRecordLayoutEnabled,
      app.linearComparisonResolution.edges.map((edge) => [edge.edgeKey, edge.ordinal, edge.source]),
      programJobs];
  });
  expect(selectedResolution).toEqual([true, true, false, [
    [`${uidA}->${uidB}`, 0, 'upload'], [`${uidB}->${uidC}`, 1, 'losat']
  ], [
    [['blastn', `${uidB}->${uidC}`, 1, 2]],
    [['tblastx', `${uidB}->${uidC}`, 1, 2]],
    [['blastp', `${uidB}->${uidC}`, 1, 2]]
  ]]);

  await page.evaluate(() => { window.__GBDRAW_MIXED_RUN__ = window.__GBDRAW_APP__.runAnalysis(); });
  await page.waitForFunction(() => window.__GBDRAW_LOSAT_EXECUTOR_CALLS__ === 1);
  await page.evaluate(() => {
    window.__GBDRAW_APP__.setLinearComparisonGlobalAction('none');
    window.__GBDRAW_RELEASE_LOSAT__();
  });
  const firstMixed = await page.evaluate(async () => {
    const app = window.__GBDRAW_APP__;
    const result = await window.__GBDRAW_MIXED_RUN__;
    return {
      result,
      modeAndDrafts: [app.linearComparisonPlan.mode, app.linearComparisonPlan.edges.length],
      comparisons: window.__GBDRAW_DIAGRAM_RUNS__.at(-1).comparisons,
      cacheInfo: app.losatCacheInfo.map((entry) => [
        entry.edgeKey, entry.ordinal, entry.queryIndex, entry.subjectIndex, entry.filename
      ]),
      telemetry: app.lastRunInfo.losatTelemetry,
      jobs: window.__GBDRAW_LOSAT_EXECUTOR_JOBS__
    };
  });
  expect(firstMixed.result).toEqual({ status: 'ok' });
  expect(firstMixed.modeAndDrafts).toEqual(['none', 3]);
  expect(requestPairs(firstMixed.comparisons)).toEqual([
    ['nucleotideBlast', 0, 1], ['nucleotideBlast', 1, 2]
  ]);
  expect(firstMixed.cacheInfo).toEqual([[
    `${uidB}->${uidC}`, 1, 1, 2, 'selected-b-to-c.tsv'
  ]]);
  expect(firstMixed.telemetry).toMatchObject({ cacheHits: 0, cacheMisses: 1, uniqueJobs: 1 });
  expect(firstMixed.jobs[0][0].slice(0, 4)).toEqual([`${uidB}->${uidC}`, 1, 1, 2]);
  expect(firstMixed.jobs[0][0][4]).not.toBe('');

  const reordered = await page.evaluate(() => {
    const app = window.__GBDRAW_APP__;
    app.linearComparisonPlan.mode = 'selected';
    app.setLinearRecordLayoutEnabled(true);
    app.moveLinearSeqUp(2);
    app.moveLinearSeqUp(1);
    return [app.linearSeqs.map((record) => record.uid), app.linearComparisonResolution.edges.map(
      (edge) => [edge.edgeKey, edge.queryIndex, edge.subjectIndex]
    ), app.losatCacheInfo.map(
      (entry) => [entry.edgeKey, entry.ordinal, entry.queryIndex, entry.subjectIndex]
    )];
  });
  expect(reordered).toEqual([[uidC, uidA, uidB], [
    [`${uidA}->${uidB}`, 1, 2], [`${uidB}->${uidC}`, 2, 0]
  ], [
    [`${uidB}->${uidC}`, 1, 2, 0]
  ]]);
  await openLinearSelectedPairs(page);
  await openLinearAdvancedComparison(page);
  const uploadPair = linearComparisonPair(page, `${uidA}->${uidB}`);
  const losatPair = linearComparisonPair(page, `${uidB}->${uidC}`);
  const losatRawResult = page.locator(
    `[data-linear-raw-result="${uidB}->${uidC}"]`
  );
  await expect(uploadPair).toBeVisible();
  await expect(losatPair).toBeVisible();
  expect(await uploadPair.evaluate((element) => (
    element.closest('[data-linear-comparison-boundary]')
      ?.dataset.linearComparisonBoundary
  ))).toBe('1->2');
  expect(await losatPair.evaluate((element) => (
    element.closest('[data-linear-comparison-boundary]')
      ?.dataset.linearComparisonBoundary
  ))).toBe('2->3');
  await expect(uploadPair.locator(
    'input[type="file"][aria-label="BLAST TSV for #2 to #3"]'
  )).toHaveCount(1);
  await expect(uploadPair).toContainText('mixed-a-to-b.tsv');
  await expect(losatRawResult.getByRole('textbox', {
    name: 'Raw LOSAT filename for #3 to #1'
  })).toHaveValue('selected-b-to-c.tsv');
  await expect(losatRawResult.getByRole('button', {
    name: 'Save Raw LOSAT TSV for #3 to #1'
  })).toBeEnabled();

  await page.evaluate((uid) => window.__GBDRAW_APP__.setLinearRecordRow(uid, 1), uidC);
  expect(await losatPair.evaluate((element) => (
    element.closest('[data-linear-comparison-boundary]')
      ?.dataset.linearComparisonBoundary
  ))).toBe('1->2');
  await expect(losatRawResult.getByRole('textbox', {
    name: 'Raw LOSAT filename for #3 to #1'
  })).toHaveValue('selected-b-to-c.tsv');
  await expect(losatRawResult.getByRole('button', {
    name: 'Save Raw LOSAT TSV for #3 to #1'
  })).toBeEnabled();

  const movedWithinRow = await page.evaluate((uid) => {
    const app = window.__GBDRAW_APP__;
    app.moveLinearRecordWithinRow(uid, -1);
    return [
      app.linearSeqs.map((record) => record.uid),
      app.losatCacheInfo.map(
        (entry) => [entry.edgeKey, entry.ordinal, entry.queryIndex, entry.subjectIndex]
      )
    ];
  }, uidA);
  expect(movedWithinRow).toEqual([[uidA, uidC, uidB], [
    [`${uidB}->${uidC}`, 1, 2, 1]
  ]]);
  await expect(losatRawResult.getByRole('button', {
    name: 'Save Raw LOSAT TSV for #3 to #2'
  })).toBeEnabled();
  await page.evaluate((uid) => window.__GBDRAW_APP__.moveLinearRecordWithinRow(uid, 1), uidA);
  await expect(losatRawResult.getByRole('button', {
    name: 'Save Raw LOSAT TSV for #3 to #1'
  })).toBeEnabled();

  await page.evaluate((uid) => window.__GBDRAW_APP__.setLinearRecordRow(uid, 3), uidC);
  expect(await losatPair.evaluate((element) => (
    element.closest('[data-linear-comparison-boundary]')
      ?.dataset.linearComparisonBoundary
  ))).toBe('2->3');
  await expect(losatRawResult.getByRole('button', {
    name: 'Save Raw LOSAT TSV for #3 to #1'
  })).toBeEnabled();
  const rawBeforeRerunPromise = page.waitForEvent('download', { timeout: 120000 });
  await losatRawResult.getByRole('button', {
    name: 'Save Raw LOSAT TSV for #3 to #1'
  }).click();
  const rawBeforeRerun = await rawBeforeRerunPromise;
  expect(rawBeforeRerun.suggestedFilename()).toBe('selected-b-to-c.tsv');
  expect(readFileSync(await rawBeforeRerun.path(), 'utf8')).toBe(
    'MixedRecB\tMixedRecC\t100\t60\t0\t0\t1\t60\t1\t60\t1e-30\t150\n'
  );
  expect(await page.evaluate(() => window.__GBDRAW_LOSAT_EXECUTOR_CALLS__)).toBe(1);

  const cachedRun = await page.evaluate(async () => {
    const app = window.__GBDRAW_APP__;
    const result = await app.runAnalysis();
    return {
      result,
      comparisons: window.__GBDRAW_DIAGRAM_RUNS__.at(-1).comparisons,
      executorCalls: window.__GBDRAW_LOSAT_EXECUTOR_CALLS__,
      telemetry: app.lastRunInfo.losatTelemetry,
      cacheInfo: app.losatCacheInfo.map((entry) => [
        entry.edgeKey, entry.queryIndex, entry.subjectIndex, entry.filename
      ])
    };
  });
  expect(cachedRun.result).toEqual({ status: 'ok' });
  expect(cachedRun.executorCalls).toBe(1);
  expect(cachedRun.telemetry).toMatchObject({ cacheHits: 1, cacheMisses: 0, uniqueJobs: 0 });
  expect(requestPairs(cachedRun.comparisons)).toEqual([
    ['nucleotideBlast', 1, 2], ['nucleotideBlast', 2, 0]
  ]);
  expect(cachedRun.cacheInfo).toEqual([[
    `${uidB}->${uidC}`, 2, 0, 'selected-b-to-c.tsv'
  ]]);
  await expect(losatRawResult.getByRole('button', {
    name: 'Save Raw LOSAT TSV for #3 to #1'
  })).toBeEnabled();

  await page.evaluate(() => { window.__GBDRAW_APP__.sessionTitle = 'linear-comparison-browser-matrix'; });
  const sessionDownloadPromise = page.waitForEvent('download', { timeout: 120000 });
  expect((await page.evaluate(() => window.__GBDRAW_APP__.saveSessionWithTitle())).status).toBe('saved');
  const sessionPath = await (await sessionDownloadPromise).path();
  await page.reload({ waitUntil: 'domcontentloaded' });
  await page.waitForFunction(() => window.__GBDRAW_APP__);
  const dialogPromise = page.waitForEvent('dialog', { timeout: 120000 });
  await page.locator('input[accept^=".json,"]').setInputFiles(sessionPath);
  const dialog = await dialogPromise;
  expect(dialog.message()).toBe('Session loaded successfully!');
  await dialog.accept();
  await page.waitForFunction(() => window.__GBDRAW_APP__?.losatCacheInfo?.length === 1);
  const restored = await page.evaluate(async () => {
    const app = window.__GBDRAW_APP__;
    const { state } = await import('./js/state.js');
    return {
      mode: app.linearComparisonPlan.mode,
      resolution: app.linearComparisonResolution.edges.map(
        (edge) => [edge.edgeKey, edge.queryIndex, edge.subjectIndex]
      ),
      cache: app.losatCacheInfo.map((entry) => [
        entry.edgeKey, entry.ordinal, entry.queryUid, entry.subjectUid,
        entry.queryIndex, entry.subjectIndex, entry.filename
      ]),
      cacheSize: state.losatCache.value.size,
      ui: {
        intent: app.linearComparisonUi.intentKey,
        losatMode: app.linearComparisonUi.activeLosatModeKey,
        losatpMode: app.linearComparisonUi.activeLosatpModeKey,
        summary: app.linearComparisonUi.summaryText,
        settings: app.linearComparisonUi.sectionKeys.settings
      }
    };
  });
  expect(restored).toEqual({
    mode: 'selected',
    resolution: [[`${uidA}->${uidB}`, 1, 2], [`${uidB}->${uidC}`, 2, 0]],
    cache: [[`${uidB}->${uidC}`, 1, uidB, uidC, 2, 0, 'selected-b-to-c.tsv']],
    cacheSize: 1,
    ui: {
      intent: 'custom',
      losatMode: 'blastn',
      losatpMode: 'collinear',
      summary: expect.stringContaining('LOSATN · 2 selected pairs · 1 LOSAT, 1 upload'),
      settings: [
        'losat-mode', 'blastn-task', 'upload-readiness',
        'result-filters', 'comparison-appearance'
      ]
    }
  });
  const restoredComparisonCard = page.locator('[data-linear-comparison-card]');
  const restoredCommands = restoredComparisonCard.getByRole('group', {
    name: 'Set all adjacent comparisons'
  });
  await expect(restoredCommands.getByRole('button')).toHaveCount(3);
  await expect(restoredComparisonCard.getByRole('status')).toContainText(
    'Current: Selected pairs (2; 1 LOSAT, 1 upload)'
  );
  await expect(restoredComparisonCard.getByText('Custom', { exact: true })).toBeVisible();
  await expect(linearSelectedPairs(page).locator(':scope > summary')).toContainText(
    'Selected pairs (2)'
  );
  await openLinearSelectedPairs(page);
  await openLinearAdvancedComparison(page);
  const restoredUploadPair = linearComparisonPair(page, `${uidA}->${uidB}`);
  const restoredLosatRawResult = page.locator(
    `[data-linear-raw-result="${uidB}->${uidC}"]`
  );
  await expect(restoredUploadPair).toContainText('mixed-a-to-b.tsv');
  await expect(restoredLosatRawResult.getByRole('textbox', {
    name: 'Raw LOSAT filename for #3 to #1'
  })).toHaveValue('selected-b-to-c.tsv');
  const rawDownloadPromise = page.waitForEvent('download', { timeout: 120000 });
  await restoredLosatRawResult.getByRole('button', {
    name: 'Save Raw LOSAT TSV for #3 to #1'
  }).click();
  const rawDownload = await rawDownloadPromise;
  expect(rawDownload.suggestedFilename()).toBe('selected-b-to-c.tsv');
  expect(readFileSync(await rawDownload.path(), 'utf8')).toBe(
    'MixedRecB\tMixedRecC\t100\t60\t0\t0\t1\t60\t1\t60\t1e-30\t150\n'
  );
});

test('Existing Gallery session restores edge-owned raw LOSAT cache', async ({ page }) => {
  test.setTimeout(180000);
  page.on('dialog', (dialog) => dialog.accept());
  await page.goto('/gbdraw/web/index.html', { waitUntil: 'domcontentloaded' });
  await page.waitForFunction(() => window.__GBDRAW_APP__);

  const imported = await page.evaluate(async () => {
    const response = await fetch(
      '/gbdraw/web/gallery/sessions/BGC0000708-BGC0000713.gbdraw-session.json'
    );
    const text = await response.text();
    const session = JSON.parse(text);
    const sourceVisibleEntries = session.losatCache.entries.filter(
      (entry) => entry.display !== false
    );
    const file = new File(
      [text],
      'BGC0000708-BGC0000713.gbdraw-session.json',
      { type: 'application/json' }
    );
    const result = await window.__GBDRAW_APP__.importSession({
      target: { files: [file], value: '' }
    });
    return {
      status: result?.status,
      sourceHasNoEdgeIdentity: sourceVisibleEntries.every(
        (entry) => !entry.edgeKey && !entry.queryUid && !entry.subjectUid
      ),
      restoredEdgeKeys: window.__GBDRAW_APP__.losatCacheInfo.map(
        (entry) => entry.edgeKey
      )
    };
  });

  expect(imported).toEqual({
    status: 'ok',
    sourceHasNoEdgeIdentity: true,
    restoredEdgeKeys: [
      'record-1->record-2',
      'record-2->record-3',
      'record-3->record-4',
      'record-4->record-5'
    ]
  });
  await openLinearSelectedPairs(page);
  await openLinearAdvancedComparison(page);
  const firstRawResult = page.locator(
    '[data-linear-raw-result="record-1->record-2"]'
  );
  await expect(firstRawResult).toContainText('Raw result ready');
  await expect(firstRawResult.getByRole('button', {
    name: 'Save Raw LOSAT TSV for #1 to #2'
  })).toBeEnabled();
});

test('Candidate render post-processing sanitizes and reapplies stable styles before commit', async ({ page }) => {
  await page.goto('/gbdraw/web/index.html', { waitUntil: 'domcontentloaded' });
  await page.waitForFunction(() => window.__GBDRAW_APP__);

  const outcome = await page.evaluate(async () => {
    const { prepareCandidateRenderCommit } = await import('./js/app/candidate-render.js');
    const stableKey = `record-1\u0000feature-1`;
    const sourceResult = {
      name: 'candidate.svg',
      format: 'svg',
      content: [
        '<svg xmlns="http://www.w3.org/2000/svg">',
        '<script>throw new Error("unsafe")</script>',
        '<path id="rendered-1" data-gbdraw-feature-id="shared-feature"',
        ' data-gbdraw-rendered-feature-id="rendered-1"',
        ' data-gbdraw-feature-part="block" fill="#000000" stroke="#000000" stroke-width="1"/>',
        '<path id="rendered-2" data-gbdraw-feature-id="shared-feature"',
        ' data-gbdraw-rendered-feature-id="rendered-2"',
        ' data-gbdraw-feature-part="block" fill="#000000" stroke="#000000" stroke-width="1"/>',
        '</svg>'
      ].join('')
    };
    const catalog = {
      schema: 3,
      items: [{
        resultIndex: 0,
        resultName: 'candidate.svg',
        recordKeys: ['record-1'],
        features: [{
          svgId: 'rendered-1',
          recordKey: 'record-1',
          biologicalFeatureId: 'feature-1'
        }, {
          svgId: 'rendered-2',
          recordKey: 'record-1',
          biologicalFeatureId: 'feature-2'
        }],
        biologicalFeatures: [{
          recordKey: 'record-1',
          biologicalFeatureId: 'feature-1',
          stableFeatureId: 'stable-feature-1',
          type: 'CDS',
          record_idx: 0,
          record_id: 'record-1',
          start: 0,
          end: 10,
          strand: 1,
          qualifiers: {}
        }, {
          recordKey: 'record-1',
          biologicalFeatureId: 'feature-2',
          stableFeatureId: 'stable-feature-2',
          type: 'CDS',
          record_idx: 0,
          record_id: 'record-1',
          start: 20,
          end: 30,
          strand: 1,
          qualifiers: {}
        }],
        orthogroups: [],
        annotations: [],
        comparisonMatches: []
      }]
    };
    const prepared = prepareCandidateRenderCommit({
      results: [sourceResult],
      catalog,
      featureColorOverrides: {
        [stableKey]: { color: '#ff00ff', caption: 'Candidate style' }
      },
      featureStrokeOverrides: {
        [stableKey]: { strokeColor: '#ff00ff', strokeWidth: 5 }
      }
    });
    const svg = new DOMParser().parseFromString(
      prepared.results[0].content,
      'image/svg+xml'
    );
    const feature = svg.getElementById('rendered-1');
    const sibling = svg.getElementById('rendered-2');
    return {
      sourceUnchanged:
        sourceResult.content.includes('<script>')
        && sourceResult.content.includes('fill="#000000"'),
      scriptRemoved: !prepared.results[0].content.includes('<script>'),
      fill: feature?.getAttribute('fill'),
      stroke: feature?.getAttribute('stroke'),
      strokeWidth: feature?.getAttribute('stroke-width'),
      siblingFill: sibling?.getAttribute('fill'),
      renderedBindingPreserved:
        feature?.getAttribute('data-gbdraw-rendered-feature-id') === 'rendered-1',
      featureCount: prepared.featureState.extractedFeatures.length
    };
  });

  expect(outcome).toEqual({
    sourceUnchanged: true,
    scriptRemoved: true,
    fill: '#ff00ff',
    stroke: '#ff00ff',
    strokeWidth: '5',
    siblingFill: '#000000',
    renderedBindingPreserved: true,
    featureCount: 2
  });
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

  await page.goto('/gbdraw/web/index.html', { waitUntil: 'domcontentloaded' });
  await page.waitForFunction(() => window.__GBDRAW_APP__);

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
    await app.requestFeatureColorChange(target, '#ff00ff');
    const dialog = {
      displayLabel: app.featureStyleScopeDialog.displayLabel,
      displayLabelSiblingCount: app.featureStyleScopeDialog.displayLabelSiblingCount
    };
    await app.handleColorScopeChoice('displayLabel');
    app.extractedFeatures
      .filter((feature) => feature.product === 'wsv360-like protein')
      .forEach((feature) => {
        const key = String(feature.stable_override_key || feature.id || '');
        if (!key) throw new Error('Expected a stable feature override key.');
        app.featureStrokeOverrides[key] = {
          strokeColor: '#ff00ff',
          strokeWidth: 5
        };
      });
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
    color: '#ff00ff',
    cap: 'wsv360-like protein'
  }]);

  const secondRun = await page.evaluate(async () => {
    const app = window.__GBDRAW_APP__;
    app.editableLabels = [{
      key: 'stale-label',
      idx: 1,
      text: 'stale',
      sourceText: 'stale',
      featureId: 'stale-feature',
      kind: 'regular',
      draftText: 'stale'
    }];
    const result = await app.runAnalysis();
    await window.Vue.nextTick();
    await window.Vue.nextTick();
    return {
      result,
      errorLog: app.errorLog,
      hasStaleEditableLabel: app.editableLabels.some(
        (entry) => entry.key === 'stale-label'
      )
    };
  });
  expect(secondRun).toEqual({
    result: { status: 'ok' },
    errorLog: null,
    hasStaleEditableLabel: false
  });

  const regenerated = await page.evaluate(() => {
    const app = window.__GBDRAW_APP__;
    const svgText = app.results?.[0]?.content || '';
    const svg = new DOMParser().parseFromString(svgText, 'image/svg+xml').documentElement;
    const featureIds = app.extractedFeatures
      .filter((feature) => feature.product === 'wsv360-like protein')
      .map((feature) => feature.svg_id);
    const featureStyles = featureIds.map((featureId) => {
      const roots = [...svg.querySelectorAll('[data-gbdraw-feature-id]')]
        .filter((element) => element.getAttribute('data-gbdraw-feature-id') === featureId);
      const elements = [...roots, ...roots.flatMap((root) => [...root.querySelectorAll('*')])];
      return {
        fills: elements
          .map((element) => String(element.getAttribute('fill') || '').toLowerCase())
          .filter((fill) => fill && fill !== 'none'),
        strokeColors: elements
          .map((element) => String(element.getAttribute('stroke') || '').toLowerCase())
          .filter(Boolean),
        strokeWidths: elements
          .map((element) => String(element.getAttribute('stroke-width') || ''))
          .filter(Boolean)
      };
    });
    const legendEntries = [...svg.querySelectorAll('[data-legend-key]')]
      .filter((entry) => entry.getAttribute('data-legend-key') === 'wsv360-like protein');
    const legendFills = legendEntries.flatMap((entry) => [...entry.querySelectorAll('[fill]')])
      .map((element) => String(element.getAttribute('fill') || '').toLowerCase())
      .filter((fill) => fill && fill !== 'none');
    return {
      featureIds,
      featureStyles,
      legendEntryCount: legendEntries.length,
      legendFills,
      rules: app.manualSpecificRules.map(({ feat, qual, val, color, cap }) => ({
        feat, qual, val, color, cap
      }))
    };
  });

  expect(regenerated.featureIds).toHaveLength(2);
  regenerated.featureStyles.forEach(({ fills, strokeColors, strokeWidths }) => {
    expect(fills).toContain('#ff00ff');
    expect(strokeColors).toContain('#ff00ff');
    expect(strokeWidths).toContain('5');
  });
  expect(regenerated.legendEntryCount).toBeGreaterThan(0);
  expect(regenerated.legendFills).toContain('#ff00ff');
  expect(regenerated.rules).toEqual(edited.rules);

  const postprocessingFailure = await page.evaluate(async () => {
    const app = window.__GBDRAW_APP__;
    const targetId = String(app.extractedFeatures?.[0]?.svg_id || '');
    app.selectedFeatureIds = new Set(targetId ? [targetId] : []);
    app.selectedFeatureAnchorId = targetId;
    app.editableLabels = [{
      key: 'postprocess-rollback-label',
      idx: 1,
      text: 'Keep after post-processing failure',
      sourceText: 'Original rollback label',
      featureId: targetId,
      kind: 'regular',
      draftText: 'Keep after post-processing failure'
    }];
    const snapshot = () => JSON.stringify({
      results: app.results,
      resultGenerationKey: app.resultGenerationKey,
      selectedResultIndex: app.selectedResultIndex,
      selectedFeatureIds: [...app.selectedFeatureIds].sort(),
      selectedFeatureAnchorId: app.selectedFeatureAnchorId,
      featureColorOverrides: app.featureColorOverrides,
      featureStrokeOverrides: app.featureStrokeOverrides,
      extractedFeatures: app.extractedFeatures,
      featureCatalog: app.featureCatalog,
      legendEntries: app.legendEntries,
      lastRunInfo: app.lastRunInfo,
      appliedPaletteName: app.appliedPaletteName,
      appliedPaletteColors: app.appliedPaletteColors,
      pendingPaletteName: app.pendingPaletteName,
      pendingPaletteColors: app.pendingPaletteColors,
      editableLabels: app.editableLabels
    });
    const before = snapshot();
    const originalSanitize = window.DOMPurify.sanitize;
    let result;
    try {
      window.DOMPurify.sanitize = () => {
        throw new Error('Forced candidate post-processing failure.');
      };
      result = await app.runAnalysis();
    } finally {
      window.DOMPurify.sanitize = originalSanitize;
    }
    await window.Vue.nextTick();
    const beforeState = JSON.parse(before);
    const afterState = JSON.parse(snapshot());
    return {
      result,
      errorSummary: String(app.errorLog?.summary || ''),
      snapshotPreserved: JSON.stringify(afterState) === before,
      changedFields: Object.keys(beforeState).filter(
        (key) => JSON.stringify(beforeState[key]) !== JSON.stringify(afterState[key])
      )
    };
  });
  expect(postprocessingFailure).toEqual({
    result: { status: 'error' },
    errorSummary: 'Forced candidate post-processing failure.',
    snapshotPreserved: true,
    changedFields: []
  });

  const staleResponse = await page.evaluate(async () => {
    const app = window.__GBDRAW_APP__;
    app.errorLog = null;
    const snapshot = () => JSON.stringify({
      results: app.results,
      resultGenerationKey: app.resultGenerationKey,
      selectedResultIndex: app.selectedResultIndex,
      selectedFeatureIds: [...app.selectedFeatureIds].sort(),
      selectedFeatureAnchorId: app.selectedFeatureAnchorId,
      featureColorOverrides: app.featureColorOverrides,
      featureStrokeOverrides: app.featureStrokeOverrides,
      extractedFeatures: app.extractedFeatures,
      featureCatalog: app.featureCatalog,
      legendEntries: app.legendEntries,
      lastRunInfo: app.lastRunInfo,
      appliedPaletteName: app.appliedPaletteName,
      appliedPaletteColors: app.appliedPaletteColors,
      pendingPaletteName: app.pendingPaletteName,
      pendingPaletteColors: app.pendingPaletteColors,
      editableLabels: app.editableLabels
    });
    const before = snapshot();
    const originalPostMessage = Worker.prototype.postMessage;
    let delayedRunObserved = false;
    Worker.prototype.postMessage = function delayedFirstRun(message, transfer) {
      if (!delayedRunObserved && message?.type === 'run') {
        delayedRunObserved = true;
        window.setTimeout(() => {
          originalPostMessage.call(this, message, transfer);
        }, 1500);
        return;
      }
      originalPostMessage.call(this, message, transfer);
    };

    let firstResult;
    let secondResult;
    try {
      const firstRun = app.runAnalysis();
      for (let attempt = 0; attempt < 500 && !delayedRunObserved; attempt += 1) {
        await new Promise((resolve) => window.setTimeout(resolve, 2));
      }
      if (!delayedRunObserved) throw new Error('The first worker request was not observed.');
      const secondRun = app.runAnalysis();
      [firstResult, secondResult] = await Promise.all([firstRun, secondRun]);
    } finally {
      Worker.prototype.postMessage = originalPostMessage;
    }
    await window.Vue.nextTick();
    const beforeState = JSON.parse(before);
    const afterState = JSON.parse(snapshot());
    return {
      firstResult,
      secondResult,
      errorLog: app.errorLog,
      snapshotPreserved: JSON.stringify(afterState) === before,
      changedFields: Object.keys(beforeState).filter(
        (key) => JSON.stringify(beforeState[key]) !== JSON.stringify(afterState[key])
      )
    };
  });
  expect(staleResponse).toEqual({
    firstResult: { status: 'stale' },
    secondResult: { status: 'error' },
    errorLog: null,
    snapshotPreserved: true,
    changedFields: []
  });

  const reset = await page.evaluate(async () => {
    const app = window.__GBDRAW_APP__;
    const target = app.extractedFeatures.find((feature) => feature.product === 'wsv360-like protein');
    if (!target) throw new Error('Expected a feature to reset.');
    const defaultColor = app.appliedPaletteColors.CDS;
    app.clickedFeature = {
      svg_id: target.svg_id,
      feat: target,
      color: '#ff00ff'
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
  expect(siblingFills).toContain('#ff00ff');
});

test('Region annotations expose and persist an explicit target-record selection', async ({ page }) => {
  await page.goto('/gbdraw/web/index.html', { waitUntil: 'domcontentloaded' });
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
  await page.goto('/gbdraw/web/index.html', { waitUntil: 'domcontentloaded' });
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
  await page.goto('/gbdraw/web/index.html', { waitUntil: 'domcontentloaded' });
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
