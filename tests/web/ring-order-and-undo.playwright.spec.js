// Circular comparison rings through History, a running Generate, and the ring order of the custom track stack.
// - OV-2445: Undo past a Track stack edit made before a Generate that searched the rings, then Generate: the
//   Generate's adoption of its ring result (the series keys, the managed stack rows, the ring files) is part of its
//   History step, so Undo, Redo, a second Undo, and a Save and Load after the Undo keep the stack rows on their rings.
// - OV-2317: while a Generate runs, the source radio's and Add comparison sequence's handlers refuse with the run's
//   busy answer, as their disabled controls show, and change no ring.
// - OV-3852 (Owner's choice A, 2026-11-18): with Use custom stack on, Move series and a move of a ring's Track stack
//   row change one ring order: the list, the stack rows and the drawn rings agree, and each move is one Undo step.
const { test, expect } = require('@playwright/test');
const path = require('node:path');
const { evaluateWithRetainedPromise, reveal } = require('./helpers/app-lifecycle.cjs');
const { load, download } = require('./helpers/mode-transition.cjs');
const { settleLive } = require('./helpers/live-generate-parity.cjs');

const root = path.resolve(process.env.GBDRAW_REPO || process.cwd());
const HMMT = path.join(root, 'gbdraw/web/gallery/sessions/HmmtDNA_basic_circular.gbdraw');
// Each case runs two or three Generates with their searches.
test.describe.configure({ timeout: 240000 });

// Add comparison sequence through its real file input.
const addRing = async (page, name) => {
  await reveal(page.getByRole('button', { name: 'Add comparison sequence', exact: true }));
  const input = page.locator('input[type=file][accept*=".fasta"][data-history-managed]').last();
  await input.setInputFiles({ name: `${name}.fasta`, mimeType: 'text/plain', buffer: Buffer.from(`>${name}\nACGTACGTACGTACGTACGTACGT\n`) });
  await page.evaluate(() => window.__GBDRAW_APP__.settleComparisonRecordLabels?.());
};
// A stub LOSAT search: each job gives one hit, so each ring draws.
const stubLosat = (page) => page.evaluate(async () => {
  const { state } = await import('./js/state.js');
  state.losatExecution.executionMode = 'serial';
  window.__GBDRAW_LOSAT_EXECUTOR__ = async (jobs) => jobs.map((job, index) => ({
    cacheKey: job.cacheKey,
    text: `${job.queryName || 'query'}\tNC_012920.1\t99\t${2100 + index}\t3\t0\t1\t${2100 + index}\t1\t${2100 + index}\t1e-20\t500\n`
  }));
});
/** HmmtDNA with rings a, b, c added by Add comparison sequence and Use custom stack on. */
const loadRings = async (browser) => {
  const page = await load(browser, HMMT);
  await stubLosat(page);
  for (const name of ['ringa', 'ringb', 'ringc']) await addRing(page, name);
  await page.evaluate(() => window.__GBDRAW_APP__.setCircularTrackSlotsEnabled(true));
  // The stack's step ends here, so Undo back to the start leaves the stack on.
  await page.evaluate(() => window.__GBDRAW_HISTORY__.settlePendingIntent());
  return page;
};
const generate = async (page) => {
  const key = await page.evaluate(async () => (await import('./js/state.js')).state.resultGenerationKey.value);
  await page.getByRole('button', { name: 'Generate Diagram', exact: true }).click();
  await page.waitForFunction(() => !window.__GBDRAW_APP__.processing, null, { timeout: 180000 });
  await settleLive(page);
  return page.evaluate(async (before) => {
    const { state } = await import('./js/state.js');
    return { generated: state.resultGenerationKey.value > before, error: state.errorLog.value?.code ?? null };
  }, key);
};
// The rings as the list, the stack rows and the drawn figure show them, and whether each ring row names a series.
const rings = (page) => page.evaluate(async () => {
  const { state } = await import('./js/state.js');
  const app = window.__GBDRAW_APP__;
  const drawing = state.activeDrawing();
  const series = drawing.circularConservation.series.map((entry) => ({ label: String(entry.label || ''), key: String(entry.sourceKey || '') }));
  const rows = drawing.adv.circular_track_slots.filter((slot) => slot.renderer === 'sequence_conservation')
    .map((slot) => ({ id: slot.id, key: String(slot.params.series_key || ''), width: slot.width ?? null }));
  const result = app.currentResults()[app.selectedResultIndex];
  const svg = result ? new DOMParser().parseFromString(String(result.content), 'image/svg+xml') : null;
  // Each drawn ring, outermost first, by the series place it draws.
  const groups = svg ? [...svg.querySelectorAll('[data-gbdraw-slot-renderer="sequence_conservation"]')] : [];
  const drawn = groups.map((group) => Number(group.getAttribute('data-source-index')));
  // Each drawn ring's stack row, in the order the figure draws them (its group id ends with the row id).
  const drawnRows = groups.map((group) => String(group.id).replace(/^track_slot_[0-9a-f]+_/, ''));
  const keys = new Set(series.map((entry) => entry.key));
  return { series, rows, drawn, drawnRows, rowsOnSeries: rows.every((row) => keys.has(row.key)) };
});
const undoCount = (page) => page.evaluate(() => window.__GBDRAW_HISTORY__.getUndoCount());
const clickHistory = async (page, name) => {
  await page.getByRole('button', { name }).first().click();
  await page.waitForFunction(() => !window.__GBDRAW_APP__.processing);
  await settleLive(page);
};
// A Track stack edit before any search: ring a's row 30 px wide (one History step).
const widenRingA = (page) => evaluateWithRetainedPromise(page, async () => {
  const { state } = await import('./js/state.js');
  await window.__GBDRAW_HISTORY__.runUndoable('Width', () => {
    state.activeDrawing().adv.circular_track_slots.find((slot) => slot.id === 'conservation_ringa').width = '30px';
  });
});

test.describe('OV-2445: Undo across a Generate that searched the rings', () => {
  for (const after of ['generate', 'redo', 'undo again', 'save and load']) {
    test(`Undo past a stack edit, then ${after}: the rows stay on their rings and Generate draws`, async ({ browser }, testInfo) => {
      let page = await loadRings(browser);
      const start = await undoCount(page);
      await widenRingA(page);
      expect(await generate(page)).toEqual({ generated: true, error: null });
      expect((await rings(page)).rowsOnSeries).toBe(true);
      while ((await undoCount(page)) > start) await clickHistory(page, /^Undo/);
      const undone = await rings(page);
      expect(undone.rowsOnSeries, JSON.stringify(undone)).toBe(true);
      expect(undone.rows.map((row) => row.width)).toEqual([null, null, null]);
      if (after === 'redo') {
        await clickHistory(page, /^Redo/);
        expect((await rings(page)).rows.map((row) => row.width)).toEqual(['30px', null, null]);
      }
      if (after === 'undo again') {
        await clickHistory(page, /^Redo/);
        await clickHistory(page, /^Undo/);
      }
      if (after === 'save and load') {
        const file = testInfo.outputPath('undone.gbdraw');
        await download(page, 'Save Session', file);
        await page.context().close();
        page = await load(browser, file);
        await stubLosat(page);
      }
      const state = await rings(page);
      expect(state.rowsOnSeries, JSON.stringify(state)).toBe(true);
      expect(await generate(page)).toEqual({ generated: true, error: null });
      const drawn = await rings(page);
      expect(drawn.rowsOnSeries, JSON.stringify(drawn)).toBe(true);
      expect(drawn.drawn).toHaveLength(3);
      await page.context().close();
    });
  }
});

test.describe('OV-2317: ring handlers during a running Generate', () => {
  for (const kind of ['source radio', 'Add comparison sequence']) {
    test(`the ${kind} handler refuses while a Generate runs, and the Generate draws its rings`, async ({ browser }) => {
      const page = await loadRings(browser);
      // Hold the Generate in its LOSAT search.
      await page.evaluate(() => {
        const search = window.__GBDRAW_LOSAT_EXECUTOR__;
        window.__GBDRAW_LOSAT_EXECUTOR__ = async (jobs) => {
          await new Promise((resolve) => { window.__RELEASE_SEARCH__ = resolve; });
          return search(jobs);
        };
      });
      const before = await rings(page);
      await page.evaluate(() => { window.__RUN__ = window.__GBDRAW_APP__.runAnalysis(); });
      await page.waitForFunction(() => window.__RELEASE_SEARCH__ && window.__GBDRAW_APP__.processing, null, { timeout: 120000 });
      const answer = await page.evaluate(async (which) => {
        const app = window.__GBDRAW_APP__;
        const file = new File(['>q\nACGTACGTACGTACGTACGTACGT\n'], 'ringq.fasta', { type: 'text/plain' });
        const result = which === 'source radio'
          ? await app.changeCircularConservationSource('upload')
          : await app.addCircularConservationComparisonFile({ target: { files: [file], value: '' } });
        return result ?? null;
      }, kind);
      expect(answer?.status).toBe('busy');
      expect(await rings(page)).toEqual(before);
      await page.evaluate(() => window.__RELEASE_SEARCH__());
      await page.evaluate(() => window.__RUN__);
      await page.waitForFunction(() => !window.__GBDRAW_APP__.processing, null, { timeout: 120000 });
      await settleLive(page);
      const after = await rings(page);
      expect(await page.evaluate(async () => (await import('./js/state.js')).state.errorLog.value?.code ?? null)).toBeNull();
      expect(after.rowsOnSeries, JSON.stringify(after)).toBe(true);
      expect(after.drawn).toHaveLength(3);
      await page.context().close();
    });
  }
});

// The list, the stack rows and the drawn rings, each as ring names (a, b, c), in their order.
const order = async (page, { drawn = false } = {}) => {
  const now = await rings(page);
  const name = (text) => text.replace(/^.*ring/, '');
  const byKey = new Map(now.series.map((entry) => [entry.key, name(entry.label)]));
  return {
    list: now.series.map((entry) => name(entry.label)),
    rows: now.rows.map((row) => byKey.get(row.key) ?? `?${row.id}`),
    ...(drawn ? { drawn: now.drawnRows.map(name) } : {})
  };
};
const one = (names) => ({ list: names, rows: names });
const settleStep = (page) => page.evaluate(() => window.__GBDRAW_HISTORY__.settlePendingIntent());
const click = async (page, label) => {
  // A stack row's buttons show with the Track stack open.
  const stack = page.locator('[data-track-slots-toggle]');
  if (label.startsWith('Circular track slot') && (await stack.getAttribute('aria-expanded')) !== 'true') {
    await stack.evaluate((element) => element.click());
  }
  const button = page.locator(`button[aria-label="${label}"]`).first();
  await expect(button).toBeEnabled();
  // The ring list and the stack rows sit in closed sections; the click is the button's own handler.
  await button.evaluate((element) => element.click());
  await settleStep(page);
};

test.describe('OV-3852: one ring order with Use custom stack on', () => {
  for (const [what, label] of [
    ['Move series', 'Move series inward'],
    ['a move of a ring\'s Track stack row', 'Circular track slot conservation_ringb Move outward']
  ]) {
    test(`${what} moves the list and the stack rows together, the figure draws that order, and Undo is one step`, async ({ browser }) => {
      const page = await loadRings(browser);
      expect(await order(page)).toEqual(one(['a', 'b', 'c']));
      const start = await undoCount(page);
      await click(page, label);
      expect(await order(page)).toEqual(one(['b', 'a', 'c']));
      expect(await undoCount(page)).toBe(start + 1);
      expect(await generate(page)).toEqual({ generated: true, error: null });
      expect(await order(page, { drawn: true })).toEqual({ ...one(['b', 'a', 'c']), drawn: ['b', 'a', 'c'] });
      await clickHistory(page, /^Undo/);
      await clickHistory(page, /^Undo/);
      expect(await order(page)).toEqual(one(['a', 'b', 'c']));
      await page.context().close();
    });
  }

  test('a Session whose stack drew its rings in another order than its list keeps the drawn order at Load', async ({ browser }, testInfo) => {
    let page = await loadRings(browser);
    // A Session an earlier version saved: the stack rows moved, its list kept (the stack drew the rows' order).
    await page.evaluate(async () => {
      const { state } = await import('./js/state.js');
      const slots = state.activeDrawing().adv.circular_track_slots;
      const at = (id) => slots.findIndex((slot) => slot.id === id);
      const [ringA] = slots.splice(at('conservation_ringa'), 1);
      slots.splice(at('conservation_ringc'), 0, ringA);
    });
    expect(await order(page)).toEqual({ list: ['a', 'b', 'c'], rows: ['b', 'a', 'c'] });
    const file = testInfo.outputPath('rows-moved.gbdraw');
    await download(page, 'Save Session', file);
    await page.context().close();
    page = await load(browser, file);
    await stubLosat(page);
    expect(await order(page)).toEqual(one(['b', 'a', 'c']));
    expect(await generate(page)).toEqual({ generated: true, error: null });
    expect(await order(page, { drawn: true })).toEqual({ ...one(['b', 'a', 'c']), drawn: ['b', 'a', 'c'] });
    await page.context().close();
  });
});
