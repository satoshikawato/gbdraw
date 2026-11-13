// Sequence search (PD-OI-117): a search adds one annotation set of its hits, drawn at once without Generate as one
// History step; Undo removes the hits at once and Redo draws them again. A hit across the origin of a circular record
// drawn whole is one hit in both modes. The chips and a row click outline hits outside the Result. A gradient and a
// search colour edit apply at once; a Legend colour edit stays. A protein query is searched with tblastn and its hits
// are drawn in genome coordinates. No hits, Cancel and a query that is no sequence add nothing.
// The live draw equals Generate in tests/web/live-generate-parity.playwright.spec.js.
const { test, expect } = require('@playwright/test');
const { readFileSync } = require('node:fs');
const { evaluateWithRetainedPromise, reveal, settleLive } = require('./helpers/app-lifecycle.cjs');
const { HMMT, loadSessionFile, openFresh } = require('./helpers/audit-browser.cjs');
const { semanticSnapshot } = require('./helpers/live-generate-parity.cjs');
const { download, load, seeds } = require('./helpers/mode-transition.cjs');
const { generate, legendRowColor, open } = require('./helpers/live-generate-parity-steps.cjs');
const {
  drawnHits, fasta, generationKey, genbankBases, reverseComplement, searchFromPanel, searchOutcome, span, undoLabels
} = require('./helpers/sequence-search.cjs');

test.describe.configure({ retries: 0 });

const OUTLINE = '[data-search-hit-outline] path.search-hit-outline';
const MT = genbankBases(HMMT).NC_012920;
// The last and first bases of the circular record (a hit across its origin), a short one whose halves fall under the
// thresholds alone, a plus-strand query and a minus-strand query.
const ORIGIN_QUERIES = [
  { name: 'origin600', bases: span(MT, MT.length - 299, 300) },
  { name: 'origin50', bases: span(MT, MT.length - 24, 25) },
  { name: 'plus3001', bases: span(MT, 3001, 3300) },
  { name: 'minus8001', bases: reverseComplement(span(MT, 8001, 8300)) }
];

// The annotation rows of a set in the shown drawing, as Region Annotations holds them.
const setRows = (page, setId) => page.evaluate(async (id) => {
  const { state } = await import('./js/state.js');
  const set = state.activeDrawing().annotationSets.find((entry) => entry.id === id);
  return set ? {
    keepHighlightOverlay: set.keepHighlightOverlay === true,
    rows: set.annotations.map(({ id: rowId, target, strand }) => ({
      id: rowId, start: target.start, end: target.end, wrapsOrigin: target.wrapsOrigin === true, strand: strand ?? null
    }))
  } : null;
}, setId);

// The outline overlay, whether the Result or its serialized SVG holds any of it, and the History depth.
const outlineState = (page) => page.evaluate(() => {
  const app = window.__GBDRAW_APP__;
  const svg = app.svgContainer.querySelector('svg');
  return {
    paths: document.querySelectorAll('[data-search-hit-outline] path.search-hit-outline').length,
    inResult: svg.querySelectorAll('.search-hit-outline, [data-search-hit-outline]').length
      + (new XMLSerializer().serializeToString(svg).includes('search-hit-outline') ? 1 : 0)
      + app.currentResults().filter((result) => String(result.content || '').includes('search-hit-outline')).length,
    undo: window.__GBDRAW_HISTORY__.getUndoCount()
  };
});

// The search of one mode: drawn at once, one History step, the hit across the origin one row, the chips and a row
// outline it outside the Result and its export, Undo and Redo apply at once. Reports the time from the History step
// to the drawn set.
const searchDrawsAtOnce = async (page, testInfo, mode) => {
  const generated = await generationKey(page);
  const steps = (await undoLabels(page)).length;
  // The time from the search's History step to its drawn set, measured in the page.
  await page.evaluate(() => {
    const start = window.__GBDRAW_HISTORY__.getUndoCount();
    window.__SEARCH_TIMING__ = {};
    const tick = () => {
      const now = performance.now();
      const timing = window.__SEARCH_TIMING__;
      if (timing.step === undefined && window.__GBDRAW_HISTORY__.getUndoCount() > start) timing.step = now;
      if (timing.step !== undefined && window.__GBDRAW_APP__.svgContainer.querySelector('svg [data-gbdraw-annotation-set-id="Search 1"]')) {
        timing.drawn = now;
        return;
      }
      requestAnimationFrame(tick);
    };
    requestAnimationFrame(tick);
  });
  await searchFromPanel(page, fasta(ORIGIN_QUERIES), { task: 'blastn' });
  const { outcome, searchError } = await searchOutcome(page);
  expect(searchError).toBeNull();
  expect(outcome.setId).toBe('Search 1');
  expect(await generationKey(page), 'drawn without Generate').toBe(generated);
  const labels = await undoLabels(page);
  // Newest first.
  expect([labels.length, labels[0]], 'one History step').toEqual([steps + 1, 'Sequence search']);
  const timing = await page.evaluate(() => window.__SEARCH_TIMING__);
  testInfo.annotations.push({ type: 'timing', description: `${mode}: search set drawn ${Math.round(timing.drawn - timing.step)} ms after its History step` });

  // One hit across the origin per crossing query, on the plus strand; the minus query hits the minus strand.
  const { keepHighlightOverlay, rows } = await setRows(page, 'Search 1');
  expect(keepHighlightOverlay).toBe(true);
  const byQuery = (name) => outcome.hits.filter((hit) => hit.query === name);
  expect(byQuery('origin600').map(({ start, end, strand }) => ({ start, end, strand })))
    .toEqual([{ start: MT.length - 299, end: 300, strand: '+' }]);
  expect(byQuery('origin50').map(({ start, end }) => ({ start, end }))).toEqual([{ start: MT.length - 24, end: 25 }]);
  expect(byQuery('plus3001').map(({ start, end, strand }) => ({ start, end, strand }))).toEqual([{ start: 3001, end: 3300, strand: '+' }]);
  expect(byQuery('minus8001').map(({ start, end, strand }) => ({ start, end, strand }))).toEqual([{ start: 8001, end: 8300, strand: '-' }]);
  const crossing = byQuery('origin600')[0].annotationId;
  expect(rows.find(({ id }) => id === crossing)).toEqual(expect.objectContaining({ wrapsOrigin: true, strand: '+' }));
  // Save Raw LOSAT TSV keeps the rows as LOSAT wrote them: the record and its origin window as subjects.
  const subjects = new Set(outcome.rawText.trim().split('\n').map((line) => line.split('\t')[1]));
  expect([...subjects].sort()).toEqual([`NC_012920.1:1-${MT.length}`, expect.stringMatching(/^NC_012920\.1:\d+-\d+$/)].sort());
  // The table lists the crossing hit once, across the origin.
  const section = page.locator('[data-settings-section="sequence-search"]');
  await expect(section.locator('[data-sequence-search-hit]')).toHaveCount(4);
  await expect(section.locator('[data-sequence-search-hit]', { hasText: 'origin600' })).toContainText(`${(MT.length - 299).toLocaleString('en-US')}–300`);
  // Every hit is drawn. The crossing one is one mark over the origin of a Circular drawing, and in Linear a mark at
  // each end of the row: a fill and the border that leaves its cut at the record's end open.
  const drawn = await drawnHits(page, 'Search 1');
  expect(Object.keys(drawn).sort()).toEqual(outcome.hits.map(({ annotationId }) => annotationId).sort());
  expect(drawn[crossing]).toBe(mode === 'circular' ? 1 : 4);

  // The chips and a row click outline hits outside the Result, without a History step; All outlines nothing.
  const depth = await page.evaluate(() => window.__GBDRAW_HISTORY__.getUndoCount());
  await expect(page.locator(OUTLINE)).toHaveCount(0);
  await section.getByRole('button', { name: 'origin600 (1)', exact: true }).click();
  await expect(page.locator(OUTLINE)).toHaveCount(drawn[crossing]);
  await expect(section.locator('[data-sequence-search-hit]')).toHaveCount(1);
  await section.getByRole('button', { name: 'All (4)', exact: true }).click();
  await expect(page.locator(OUTLINE)).toHaveCount(0);
  await section.locator('[data-sequence-search-hit]', { hasText: 'plus3001' }).click();
  await expect(page.locator(OUTLINE)).toHaveCount(1);
  expect(await outlineState(page)).toEqual({ paths: 1, inResult: 0, undo: depth });
  const exported = (await download(page, 'SVG', testInfo.outputPath(`${mode}-outlined.svg`))).toString('utf8');
  expect(exported).toContain('data-gbdraw-annotation-set-id="Search 1"');
  expect(exported).not.toContain('search-hit-outline');

  // Undo removes the hits at once, Redo draws them again; neither runs Generate.
  await page.getByRole('button', { name: 'Undo', exact: true }).click();
  await expect.poll(() => drawnHits(page, 'Search 1')).toEqual({});
  await settleLive(page);
  await expect(page.locator(OUTLINE)).toHaveCount(0);
  await page.getByRole('button', { name: 'Redo', exact: true }).click();
  await expect.poll(() => drawnHits(page, 'Search 1')).toEqual(drawn);
  await settleLive(page);
  expect(await generationKey(page)).toBe(generated);
};

test('Circular: a search draws its hits at once in one History step, across the origin too, and Undo and Redo apply at once', async ({ browser }, testInfo) => {
  test.setTimeout(600_000);
  const page = await load(browser, seeds.circular);
  try {
    await searchDrawsAtOnce(page, testInfo, 'circular');
    expect(page.externalRequests).toEqual([]);
  } finally {
    await page.context().close();
  }
});

// OV-2777: a number keeps the English digits and grouping in every language (PD-OI-106), whatever the browser's
// language: here a German and a Bengali browser (4.641.652 and ৪৬,৪১,৬৫২ before) with the UI in English. The query's
// length, a hit's Start–End and a feature's Length are three of the displays that read the browser's language before.
for (const locale of ['de-DE', 'bn-BD']) {
  test(`a ${locale} browser shows the query length, a hit's Start–End and a feature's Length with English digits and grouping`, async ({ browser }) => {
    test.setTimeout(600_000);
    const page = await load(browser, seeds.circular, undefined, { locale });
    try {
      expect(await page.evaluate(() => [navigator.language, (1300).toLocaleString()])).toEqual([locale, expect.not.stringMatching(/^1,300$/)]);
      const query = fasta([{ name: 'plus3001', bases: span(MT, 3001, 4300) }]);
      const section = page.locator('[data-settings-section="sequence-search"]');
      // The panel lists the read queries until a search ends.
      await page.getByRole('button', { name: 'Go to Sequence search', exact: true }).click();
      const unfold = page.locator('[data-sequence-search-unfold]');
      if (await unfold.count()) await unfold.click();
      await page.locator('#sequence-search-query').fill(query);
      await page.locator('#sequence-search-query').dispatchEvent('change');
      // Soft: each of the three displays reports its own text.
      await expect.soft(section.locator('[data-sequence-search-query]', { hasText: 'plus3001' })).toContainText('1,300 bp', { timeout: 60_000 });
      await searchFromPanel(page, query, { task: 'blastn' });
      expect((await searchOutcome(page)).searchError).toBeNull();
      await expect.soft(section.locator('[data-sequence-search-hit]', { hasText: 'plus3001' }).first()).toContainText('3,001–4,300');
      // COX1 is 5904..7445: 1,542 bp.
      await evaluateWithRetainedPromise(page, async () => {
        const app = window.__GBDRAW_APP__;
        await app.openFeatureEditorFromList(app.extractedFeatures.find((feature) => feature.type === 'CDS' && feature.gene === 'COX1'), null);
      });
      const details = page.getByRole('dialog', { name: /^Feature details:/ });
      await details.getByRole('button', { name: 'Details', exact: true }).click();
      await expect.soft(details.locator('[data-feature-popup-details] dt', { hasText: 'Length' }).locator('xpath=following-sibling::dd[1]')).toHaveText('1,542 bp');
    } finally {
      await page.context().close();
    }
  });
}

test('Linear: a search of a circular record drawn whole draws a hit across the origin as one hit at both ends', async ({ page }, testInfo) => {
  test.setTimeout(600_000);
  await open(page, { mode: 'linear', results: 'single', reflow: 'off', genbank: readFileSync(HMMT, 'utf8') });
  await searchDrawsAtOnce(page, testInfo, 'linear');
  // The crossing hit's marks lie at the two ends of the record's row, none in its middle.
  const ends = await page.evaluate(() => {
    const app = window.__GBDRAW_APP__;
    const svg = app.svgContainer.querySelector('svg');
    const hit = app.sequenceSearch.shownOutcome.hits.find((entry) => entry.query === 'origin600');
    const group = [...svg.querySelectorAll('g[data-gbdraw-annotation-id]')].find((entry) => entry.getAttribute('data-gbdraw-annotation-id') === hit.annotationId);
    const marks = [...group.querySelectorAll('path, rect')].map((mark) => mark.getBoundingClientRect());
    const record = [...svg.querySelectorAll('[data-gbdraw-feature-id]')].map((feature) => feature.getBoundingClientRect());
    const recordLeft = Math.min(...record.map((box) => box.left));
    const recordRight = Math.max(...record.map((box) => box.right));
    const middle = 0.5 * (recordLeft + recordRight);
    return {
      first: Math.min(...marks.map((box) => box.left)), last: Math.max(...marks.map((box) => box.right)), recordLeft, recordRight,
      overMiddle: marks.filter((box) => box.left < middle && box.right > middle).length,
    };
  });
  expect(ends.first).toBeLessThan(ends.recordLeft + 40);
  expect(ends.last).toBeGreaterThan(ends.recordRight - 40);
  expect(ends.overMiddle).toBe(0);
});

// With Use custom stack the search's set is drawn at once in front of the features, and the stack does not change.
test('Circular custom stack: a search set is drawn at once in front of the features and the stack stays as it was', async ({ browser }) => {
  test.setTimeout(600_000);
  const page = await load(browser, seeds.circular);
  try {
    await page.evaluate(() => window.__GBDRAW_APP__.setCircularTrackSlotsEnabled(true));
    await generate(page);
    const stack = () => page.evaluate(() => JSON.parse(JSON.stringify(window.__GBDRAW_APP__.adv.circular_track_slots)));
    const before = await stack();
    const generated = await generationKey(page);
    await searchFromPanel(page, fasta(ORIGIN_QUERIES.slice(2)));
    expect(await generationKey(page)).toBe(generated);
    expect(Object.keys(await drawnHits(page, 'Search 1'))).toHaveLength(2);
    expect(await stack()).toEqual(before);
    // The set's slot is drawn after (in front of) the feature track.
    const order = await page.evaluate(() => {
      const svg = window.__GBDRAW_APP__.svgContainer.querySelector('svg');
      const all = [...svg.querySelectorAll('*')];
      const hit = svg.querySelector('[data-gbdraw-annotation-set-id="Search 1"]');
      const feature = [...svg.querySelectorAll('[data-gbdraw-feature-id]')].pop();
      return { hit: all.indexOf(hit), feature: all.indexOf(feature) };
    });
    expect(order.hit).toBeGreaterThan(order.feature);
    await page.getByRole('button', { name: 'Undo', exact: true }).click();
    await expect.poll(() => drawnHits(page, 'Search 1')).toEqual({});
    expect(await stack()).toEqual(before);
  } finally {
    await page.context().close();
  }
});

// The gradient by query position: one Legend row "Query position"; a change of a hit's gradient end colour in Region
// Annotations applies at once.
test('Linear gradient: the Legend has one Query position row and a gradient colour edit applies at once', async ({ browser }) => {
  test.setTimeout(600_000);
  const page = await load(browser, seeds.linear);
  try {
    const lambda = genbankBases('tests/test_inputs/NC_001416.gb').NC_001416;
    await searchFromPanel(page, fasta([{ name: 'cosL', bases: span(lambda, 5001, 5400) }, { name: 'rev', bases: reverseComplement(span(lambda, 20001, 20300)) }]),
      { colorMode: 'gradient' });
    expect(Object.keys(await drawnHits(page, 'Search 1'))).toHaveLength(2);
    // The panel shows the Legend's bar while it is unfolded.
    await page.locator('[data-sequence-search-unfold]').click();
    await expect(page.locator('[data-sequence-search-gradient]')).toBeVisible();
    // The drawn rows, as the parity check reads them (a Linear figure holds its Legend in two orientations, one hidden).
    const { legend } = await semanticSnapshot(page);
    expect(legend.filter(({ caption }) => caption === 'Query position')).toHaveLength(1);
    // A palette change keeps the bar's colours: they are the set's, not the palette's identity colours (OV-2621).
    const queryBar = async () => (await semanticSnapshot(page)).legend.filter(({ caption }) => caption === 'Query position');
    const bar = await queryBar();
    await page.evaluate(async () => {
      const app = window.__GBDRAW_APP__;
      app.setPaletteInstantPreview(true);
      await app.selectPalette('emerald_twilight');
    });
    await settleLive(page);
    expect(await queryBar(), 'the Query position bar after a palette change').toEqual(bar);
    const generated = await generationKey(page);
    const steps = (await undoLabels(page)).length;
    // The first hit's drawn steps, from its start colour to its end colour.
    const steps0 = () => page.evaluate(() => {
      const app = window.__GBDRAW_APP__;
      const svg = app.svgContainer.querySelector('svg');
      const id = app.sequenceSearch.shownOutcome.hits[0].annotationId;
      const group = [...svg.querySelectorAll('g[data-gbdraw-annotation-set-id="Search 1"][data-gbdraw-annotation-id]')]
        .find((entry) => entry.getAttribute('data-gbdraw-annotation-id') === id);
      return [...group.querySelectorAll('[fill]')].map((mark) => mark.getAttribute('fill').toLowerCase()).filter((fill) => fill !== 'none');
    });
    const committedEnd = () => page.evaluate(async () => {
      const { getCommittedCanonicalRenderRequest } = await import('./js/services/config.js');
      const set = getCommittedCanonicalRenderRequest().diagramOptions.annotations.sets.find((entry) => entry.id === 'Search 1');
      return (set.annotations[0].style ?? set.defaultStyle).gradient.endColor;
    });
    const before = await steps0();
    const end = page.locator('input[data-annotation-gradient="end"]').first();
    await reveal(end);
    await end.scrollIntoViewIfNeeded();
    await end.fill('#123456');
    // The applied request carries the new end colour, and the hit is drawn again from it, without Generate.
    await expect.poll(committedEnd, { timeout: 120_000 }).toBe('#123456');
    await expect.poll(steps0, { timeout: 120_000 }).not.toEqual(before);
    await settleLive(page);
    expect(await generationKey(page), 'applied without Generate').toBe(generated);
    expect((await undoLabels(page)).length).toBe(steps + 1);
  } finally {
    await page.context().close();
  }
});

// A search's set takes a colour edit at once (its rows carry a Legend label), and a Legend colour edit made before
// stays on its row and its features.
test('Circular: a search colour edit applies at once and keeps a Legend colour edit', async ({ browser }) => {
  test.setTimeout(600_000);
  const page = await load(browser, seeds.circular);
  try {
    await legendRowColor(page, 'CDS', '#7b2cbf');
    await searchFromPanel(page, fasta(ORIGIN_QUERIES.slice(2)));
    const generated = await generationKey(page);
    const legend = () => page.evaluate(() => Object.fromEntries(window.__GBDRAW_APP__.legendEntries.map((entry) => [entry.caption, entry.color])));
    expect((await legend()).CDS).toBe('#7b2cbf');
    const section = page.locator('[data-setting="region-annotations"]');
    const fill = section.locator('input[aria-label="Annotation fill color"]').last();
    await reveal(fill);
    await fill.scrollIntoViewIfNeeded();
    await fill.fill('#00aa55');
    const hitFills = () => page.evaluate(() => {
      const svg = window.__GBDRAW_APP__.svgContainer.querySelector('svg');
      return [...svg.querySelectorAll('[data-gbdraw-annotation-set-id="Search 1"] [fill]')].map((mark) => mark.getAttribute('fill').toLowerCase());
    });
    await expect.poll(hitFills, { timeout: 120_000 }).toContain('#00aa55');
    await settleLive(page);
    expect(await generationKey(page), 'applied without Generate').toBe(generated);
    expect((await legend()).CDS, 'the Legend edit stays').toBe('#7b2cbf');
    const cdsFills = await page.evaluate(() => {
      const svg = window.__GBDRAW_APP__.svgContainer.querySelector('svg');
      const row = [...svg.querySelectorAll('g[data-legend-key]')].find((entry) => entry.textContent.trim() === 'CDS');
      return row ? [...row.querySelectorAll('[fill]')].map((mark) => mark.getAttribute('fill').toLowerCase()) : [];
    });
    expect(cdsFills).toContain('#7b2cbf');
  } finally {
    await page.context().close();
  }
});

// Before the first Generate Search says why it is off; a protein query is refused with its reason under the query box;
// a cancelled search, a failed search (dc-megablast on serial LOSAT) and a search without hits add nothing and keep the
// Result.
test('a query that is no sequence, Cancel, dc-megablast in Serial and a search without hits add nothing and keep the Result', async ({ page }) => {
  test.setTimeout(600_000);
  await openFresh(page);
  await page.getByRole('button', { name: 'Go to Sequence search', exact: true }).click();
  await expect(page.locator('#sequence-search-reason')).toHaveText('Generate to search');
  await expect(page.locator('[data-sequence-search-run]')).toBeDisabled();
  await loadSessionFile(page, seeds.linear);
  await settleLive(page);
  const result = () => page.evaluate(() => window.__GBDRAW_APP__.currentResults()[0].content);
  const before = { result: await result(), steps: (await undoLabels(page)).length, sets: await page.evaluate(() => window.__GBDRAW_APP__.annotationSets.length) };
  const unchanged = async () => {
    expect(await result()).toBe(before.result);
    expect((await undoLabels(page)).length).toBe(before.steps);
    expect(await page.evaluate(() => window.__GBDRAW_APP__.annotationSets.length)).toBe(before.sets);
  };

  // A query holding a character that is neither a nucleotide nor an amino acid code: the reason under the query box,
  // and Search off.
  await page.getByRole('button', { name: 'Go to Sequence search', exact: true }).click();
  const query = page.locator('#sequence-search-query');
  await query.fill('>prot\nMKVLAAGIVGLLL-AQSEWQRTYHPLIDDFEKMNQ\n');
  await query.press('Tab');
  // The first read starts Python's query reader in the diagram Worker.
  await expect(page.locator('#sequence-search-query-error')).toBeVisible({ timeout: 120_000 });
  await expect(page.locator('#sequence-search-query-error')).toContainText('First character that is neither a nucleotide nor an amino acid code: "-".');
  await expect(page.locator('#sequence-search-reason')).toHaveText('Correct the query to search');
  await expect(page.locator('[data-sequence-search-run]')).toBeDisabled();
  await unchanged();

  // Cancel while LOSAT loads: the full-window overlay shows the search, and Cancel adds nothing.
  const lambda = genbankBases('tests/test_inputs/NC_001416.gb').NC_001416;
  /** @type {() => void} */
  let release = () => {};
  const held = new Promise((resolve) => { release = resolve; });
  await page.route('**/wasm/losat/**', async (route) => { await held; await route.continue(); });
  await page.evaluate(() => { window.__GBDRAW_APP__.losatExecution.executionMode = 'serial'; });
  await query.fill(fasta([{ name: 'cosL', bases: span(lambda, 5001, 5400) }]));
  await page.locator('[data-sequence-search-run]').click();
  // Python reads the query and the record's bases first; the overlay shows once LOSAT starts.
  await expect.poll(() => page.evaluate(() => {
    const app = window.__GBDRAW_APP__;
    const { panel } = app.sequenceSearch;
    return { processing: app.processing, running: panel.running, error: panel.searchError?.summary ?? panel.queryError?.summary ?? null };
  }), { timeout: 120_000 }).toEqual({ processing: true, running: true, error: null });
  await expect(page.getByText('Searching sequences...', { exact: true })).toBeVisible();
  await page.getByRole('button', { name: /Cancel$/ }).click();
  release();
  await expect.poll(() => page.evaluate(() => window.__GBDRAW_APP__.sequenceSearch.panel.running), { timeout: 120_000 }).toBe(false);
  await settleLive(page);
  await page.unroute('**/wasm/losat/**');
  expect((await searchOutcome(page)).outcome).toBeNull();
  await unchanged();

  // dc-megablast needs threaded LOSAT, as in Pairwise Comparisons: a Serial search fails with that reason, before any
  // search, and adds nothing (a failed search, not a cancelled one).
  await searchFromPanel(page, fasta([{ name: 'cosL', bases: span(lambda, 5001, 5400) }]), { task: 'dc-megablast' });
  expect((await searchOutcome(page)).searchError?.summary).toContain(
    'LOSATN task dc-megablast runs only with threaded LOSAT. Select Threaded or Auto under Execution in a browser where Threaded is available, or choose task megablast or blastn.'
  );
  expect((await searchOutcome(page)).outcome).toBeNull();
  await unchanged();

  // No hit passes the thresholds: the panel says so, and nothing is added.
  const unrelated = Array.from({ length: 300 }, (_, index) => 'ACGT'[(index * 7 + Math.floor(index / 5)) % 4]).join('');
  await searchFromPanel(page, fasta([{ name: 'nohit', bases: unrelated }]), { task: 'blastn' });
  await expect(page.locator('[data-sequence-search-no-hits]')).toHaveText(/^No hits for nohit \(300 bp\) The thresholds removed \d+ alignments?\.$/);
  await unchanged();
});

// PD-OI-117: a protein query is searched with tblastn against the same records, beside a nucleotide query searched
// with blastn, and its hit is drawn in genome coordinates on the strand of its gene. Here the cI repressor of lambda,
// whose gene lies on the minus strand at 37227..37940 (the stop codon is not in the protein).
// The protein query is pasted first; blastn still runs first, so its rows lead the raw table and its name.
test('Linear: a protein query is searched with tblastn and its hit is drawn on its gene, beside a blastn hit', async ({ page }) => {
  test.setTimeout(600_000);
  await openFresh(page);
  await loadSessionFile(page, seeds.linear);
  await settleLive(page);
  const lambdaText = readFileSync('tests/test_inputs/NC_001416.gb', 'utf8');
  const cI = lambdaText.match(/\/gene="cI"\n(?:\s+\/[^\n]*\n)*?\s+\/translation="([^"]+)"/)[1].replace(/\s+/g, '');
  const lambda = genbankBases('tests/test_inputs/NC_001416.gb').NC_001416;
  await searchFromPanel(page, fasta([{ name: 'cI', bases: cI }, { name: 'cosL', bases: span(lambda, 5001, 5400) }]), { task: 'blastn' });
  const { outcome, searchError } = await searchOutcome(page);
  expect(searchError).toBeNull();
  expect(outcome.setId).toBe('Search 1');
  const section = page.locator('[data-settings-section="sequence-search"]');
  await section.locator('[data-sequence-search-unfold]').click();
  await expect(section.locator('[data-sequence-search-programs]')).toHaveText('blastn, tblastn');
  await expect(section.locator('[data-sequence-search-query]', { hasText: 'cI' })).toContainText(`${cI.length} aa`);
  const hit = (name) => outcome.hits.filter((entry) => entry.query === name).map(({ start, end, strand, qstart, qend }) => ({ start, end, strand, qstart, qend }));
  expect(hit('cI')).toEqual([{ start: 37230, end: 37940, strand: '-', qstart: 1, qend: cI.length }]);
  expect(hit('cosL')).toEqual([{ start: 5001, end: 5400, strand: '+', qstart: 1, qend: 400 }]);
  expect(outcome.rawText.split('\t', 1)[0]).toBe('cosL');
  const [proteinHit] = outcome.hits.filter((entry) => entry.query === 'cI');
  await expect(section.locator('[data-sequence-search-hit]', { hasText: 'cI' }).locator('td').nth(4)).toHaveText(`${proteinHit.alignmentLength} aa`);
  await expect(section.locator('[data-sequence-search-hit]', { hasText: 'cosL' }).locator('td').nth(4)).toHaveText('400');
  const [saved] = await Promise.all([page.waitForEvent('download'), section.getByRole('button', { name: 'Save Raw LOSAT TSV' }).click()]);
  expect(saved.suggestedFilename()).toBe('Search_1.losatn-tlosatn.tsv');
  const { rows } = await setRows(page, 'Search 1');
  expect(rows.map(({ start, end, strand }) => ({ start, end, strand }))).toEqual(expect.arrayContaining([
    { start: 37230, end: 37940, strand: '-' }, { start: 5001, end: 5400, strand: '+' }
  ]));
  const programs = await page.evaluate(() => Object.fromEntries(window.__GBDRAW_APP__.annotationSets
    .find((set) => set.id === 'Search 1').annotations.map((row) => [row.metadata.query, row.metadata.program])));
  expect(programs).toEqual({ cosL: 'blastn', cI: 'tblastn' });
  const drawn = await drawnHits(page, 'Search 1');
  expect(Object.keys(drawn).sort()).toEqual(outcome.hits.map(({ annotationId }) => annotationId).sort());
});

// Receipt: the panel's queries and settings are not saved; after a Load the panel is empty, and the hits of earlier
// searches are in Region Annotations (review B1). A Web Save then Load keeps the search set: its gradient, its rows'
// strands, and its overlay marker, and draws it again (review N7).
test('Circular: a Web Save then Load keeps a search set and empties the panel', async ({ browser }, testInfo) => {
  test.setTimeout(600_000);
  const page = await load(browser, seeds.circular);
  try {
    await searchFromPanel(page, fasta(ORIGIN_QUERIES.slice(2)), { colorMode: 'gradient' });
    await page.locator('[data-settings-section="sequence-search"] [data-sequence-search-hit]', { hasText: 'plus3001' }).click();
    await expect(page.locator(OUTLINE)).toHaveCount(1);
    const set = () => page.evaluate(() => JSON.parse(JSON.stringify(
      window.__GBDRAW_APP__.annotationSets.find((entry) => entry.id === 'Search 1')
    )));
    const saved = await set();
    expect(saved.keepHighlightOverlay).toBe(true);
    expect(saved.annotations.map(({ strand, style }) => [strand, Boolean(style?.gradient)])).toEqual([['+', true], ['-', true]]);
    const drawn = await drawnHits(page, 'Search 1');
    const path = testInfo.outputPath('search.gbdraw');
    await download(page, 'Save Session', path);
    await loadSessionFile(page, path);
    await settleLive(page);
    const panel = await page.evaluate(() => {
      const search = window.__GBDRAW_APP__.sequenceSearch;
      const { queryText, queries, outcome, folded, chip, selectedHit } = search.panel;
      return { queryText, queries: queries.length, outcome, folded, chip, selectedHit, rows: search.hitRows.length, chips: search.chips.length };
    });
    expect(panel).toEqual({ queryText: '', queries: 0, outcome: null, folded: false, chip: '', selectedHit: '', rows: 0, chips: 0 });
    await expect(page.locator('[data-settings-section="sequence-search"] [data-sequence-search-hit]')).toHaveCount(0);
    await expect(page.locator(OUTLINE)).toHaveCount(0);
    expect(await set()).toEqual(saved);
    expect(await drawnHits(page, 'Search 1')).toEqual(drawn);
  } finally {
    await page.context().close();
  }
});

// Ticks or unticks the record-list box of the Linear card at `index` (a record ON/OFF waits for Generate).
const drawCard = async (page, index, on) => {
  const uid = await page.evaluate((at) => window.__GBDRAW_APP__.linearSeqs[at].uid, index);
  const box = page.locator(`input[data-record-draw="${uid}"]`);
  await reveal(box);
  await box.scrollIntoViewIfNeeded();
  if (on) await box.check(); else await box.uncheck();
};

// The record index each drawn hit of a set sits on in the shown Result (Python's `data-gbdraw-record-index`).
const hitRecords = (page, setId) => page.evaluate((id) => [...window.__GBDRAW_APP__.svgContainer.querySelector('svg')
  .querySelectorAll('g[data-gbdraw-annotation-id]')]
  .filter((group) => group.getAttribute('data-gbdraw-annotation-set-id') === id)
  .map((group) => group.getAttribute('data-gbdraw-record-index')), setId);

// Review N1: a search's rows bind their record through the record catalog, as Add to annotation set does, so the hits
// of a record turned OFF wait with it: Generate draws the other record, and the hits come back with the record.
test('Linear: the hits of a record turned OFF wait with it and come back with it', async ({ page }) => {
  test.setTimeout(600_000);
  const lambdaText = readFileSync('tests/test_inputs/NC_001416.gb', 'utf8');
  await open(page, { mode: 'linear', results: 'single', reflow: 'off', genbank: `${readFileSync(HMMT, 'utf8')}${lambdaText}` });
  const lambda = genbankBases('tests/test_inputs/NC_001416.gb').NC_001416;
  await searchFromPanel(page, fasta([{ name: 'cosL', bases: span(lambda, 5001, 5400) }]));
  expect((await searchOutcome(page)).searchError).toBeNull();
  const drawn = await drawnHits(page, 'Search 1');
  expect(Object.keys(drawn)).toHaveLength(1);
  const draw = async (on) => drawCard(page, await page.evaluate(() => window.__GBDRAW_APP__.linearSeqs.findIndex((seq) => seq.region_record_id.startsWith('NC_001416'))), on);
  await draw(false);
  // Generate fails when a row names a record the request no longer draws (`generate` expects an ok run).
  await generate(page);
  expect(await drawnHits(page, 'Search 1')).toEqual({});
  await draw(true);
  await generate(page);
  expect(await drawnHits(page, 'Search 1')).toEqual(drawn);
});

// Review r2 B1, receipt: the search runs against the records the shown drawing draws, as last generated. A record
// turned OFF since the Generate (which waits for Generate) neither refuses the search nor moves its hits: they are
// drawn at once on the record the Result draws, and Generate then draws them on the draft's records.
test('Linear: a record turned OFF after Generate does not refuse a search of the other record', async ({ page }) => {
  test.setTimeout(600_000);
  const lambdaText = readFileSync('tests/test_inputs/NC_001416.gb', 'utf8');
  await open(page, { mode: 'linear', results: 'single', reflow: 'off', genbank: `${readFileSync(HMMT, 'utf8')}${lambdaText}` });
  const lambda = genbankBases('tests/test_inputs/NC_001416.gb').NC_001416;
  await drawCard(page, 0, false);
  await searchFromPanel(page, fasta([{ name: 'cosL', bases: span(lambda, 5001, 5400) }]));
  expect((await searchOutcome(page)).searchError).toBeNull();
  expect(await hitRecords(page, 'Search 1')).toEqual(['1']);
  await generate(page);
  expect(await hitRecords(page, 'Search 1')).toEqual(['0']);
});

// `count` records of one file that share the ID contig1 (lambda's bases, 4 kb each), the first turned OFF after
// Generate: a hit on the second is drawn on the second at once (the Result's `#2`, not the draft's), and on it again
// after Generate.
const genbankRecord = (id, bases) => {
  const lines = [];
  for (let at = 0; at < bases.length; at += 60) {
    lines.push(`${String(at + 1).padStart(9)} ${bases.slice(at, at + 60).toLowerCase().match(/.{1,10}/g).join(' ')}`);
  }
  return `LOCUS       ${id.padEnd(16)}${String(bases.length).padStart(12)} bp    DNA     linear   UNK 01-JAN-1980
DEFINITION  ${id}.
ACCESSION   ${id}
VERSION     ${id}
FEATURES             Location/Qualifiers
     CDS             1..300
                     /product="test protein"
ORIGIN
${lines.join('\n')}
//
`;
};
for (const count of [3, 2]) {
  test(`Linear: ${count} records named contig1, the first turned OFF after Generate: a hit on the second is drawn on it`, async ({ page }) => {
    test.setTimeout(600_000);
    const lambda = genbankBases('tests/test_inputs/NC_001416.gb').NC_001416;
    const contigs = Array.from({ length: count }, (_unused, index) => genbankRecord('contig1', span(lambda, index * 4000 + 1, (index + 1) * 4000)));
    await open(page, { mode: 'linear', results: 'single', reflow: 'off', genbank: contigs.join('') });
    expect(await page.evaluate(() => window.__GBDRAW_APP__.linearSeqs.length)).toBe(count);
    await drawCard(page, 0, false);
    await searchFromPanel(page, fasta([{ name: 'second', bases: span(lambda, 5001, 5400) }]));
    expect((await searchOutcome(page)).searchError).toBeNull();
    expect(await hitRecords(page, 'Search 1')).toEqual(['1']);
    await generate(page);
    expect(await hitRecords(page, 'Search 1')).toEqual(['0']);
  });
}

// OV-2631, D-57: three cards of records named contig1; after Generate the second card is deleted, so the draft's `#2`
// is the third record. A hit on the second record (drawn at once on it, as the Result draws it) is held as "target
// record unavailable": Generate refuses it instead of drawing it on the third record. Deleting the set lets Generate
// pass.
test('Linear: a hit on a record whose card was deleted after Generate is refused at Generate, not moved', async ({ page }) => {
  test.setTimeout(600_000);
  const lambda = genbankBases('tests/test_inputs/NC_001416.gb').NC_001416;
  const contigs = [0, 1, 2].map((index) => genbankRecord('contig1', span(lambda, index * 4000 + 1, (index + 1) * 4000)));
  await open(page, { mode: 'linear', results: 'single', reflow: 'off', genbank: contigs[0] });
  await page.evaluate(async (texts) => {
    const app = window.__GBDRAW_APP__;
    texts.forEach((text, index) => {
      app.addLinearSeq();
      app.setLinearSeqPrimaryFile(index + 1, 'gb', new File([text], `contig-${index + 2}.gb`, { type: 'text/plain', lastModified: 1000 }));
    });
    await window.Vue.nextTick();
  }, contigs.slice(1));
  await generate(page);
  expect(await page.evaluate(() => window.__GBDRAW_APP__.linearSeqs.length)).toBe(3);
  await page.evaluate(() => { const app = window.__GBDRAW_APP__; app.requestLinearSourceRemoval(app.linearSourceGroups[1]); });
  await page.locator('[data-linear-source-removal-dialog]').getByRole('button', { name: 'Remove file', exact: true }).click();
  await expect.poll(() => page.evaluate(() => window.__GBDRAW_APP__.linearSeqs.length)).toBe(2);
  await searchFromPanel(page, fasta([{ name: 'second', bases: span(lambda, 5001, 5400) }]));
  expect((await searchOutcome(page)).searchError).toBeNull();
  expect(await hitRecords(page, 'Search 1')).toEqual(['1']);
  const rejected = await evaluateWithRetainedPromise(page, () => window.__GBDRAW_APP__.runAnalysis());
  expect(rejected.status).toBe('error');
  expect(rejected.error).toMatchObject({ code: 'ANNOTATION_TARGET', context: { reason: 'TARGET_RECORD' } });
  await expect(page.getByRole('alert', { name: 'Generation Error' })).toContainText('Choose an available target record.');
  await page.evaluate(() => { const app = window.__GBDRAW_APP__; app.removeAnnotationSet(app.annotationSets.find(({ id }) => id === 'Search 1')); });
  await generate(page);
  expect(await hitRecords(page, 'Search 1')).toEqual([]);
});
