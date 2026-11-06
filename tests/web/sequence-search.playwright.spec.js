// Sequence search (PD-OI-117): a search adds one annotation set of its hits, drawn at once without Generate as one
// History step; Undo removes the hits at once and Redo draws them again. A hit across the origin of a circular record
// drawn whole is one hit in both modes. The chips and a row click outline hits outside the Result. A gradient and a
// search colour edit apply at once; a Legend colour edit stays. No hits, Cancel and a protein query add nothing.
// The live draw equals Generate in tests/web/live-generate-parity.playwright.spec.js.
const { test, expect } = require('@playwright/test');
const { readFileSync } = require('node:fs');
const { reveal, settleLive } = require('./helpers/app-lifecycle.cjs');
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
// a cancelled search and a search without hits add nothing and keep the Result.
test('a protein query, Cancel and a search without hits add nothing and keep the Result', async ({ page }) => {
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

  // A protein query: the reason under the query box, and Search off.
  await page.getByRole('button', { name: 'Go to Sequence search', exact: true }).click();
  const query = page.locator('#sequence-search-query');
  await query.fill('>prot\nMKVLAAGIVGLLLAQSEWQRTYHPLIDDFEKMNQ\n');
  await query.press('Tab');
  // The first read starts Python's query reader in the diagram Worker.
  await expect(page.locator('#sequence-search-query-error')).toBeVisible({ timeout: 120_000 });
  await expect(page.locator('#sequence-search-query-error')).not.toHaveText('');
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

  // No hit passes the thresholds: the panel says so, and nothing is added.
  const unrelated = Array.from({ length: 300 }, (_, index) => 'ACGT'[(index * 7 + Math.floor(index / 5)) % 4]).join('');
  await searchFromPanel(page, fasta([{ name: 'nohit', bases: unrelated }]), { task: 'blastn' });
  await expect(page.locator('[data-sequence-search-no-hits]')).toHaveText(/^No hits for nohit \(300 bp\) The thresholds removed \d+ alignments?\.$/);
  await unchanged();
});
