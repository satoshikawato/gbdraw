// Sequence search (PD-OI-116) steps shared by the search specs, the live-vs-Generate parity cases, the R3 work
// allowlist and the random walk. Queries are cut from the GenBank file the drawing draws, so a test names the bases
// it expects to hit; the search itself runs LOSAT and Python as a user's search does. The test server is not
// cross-origin isolated, so LOSAT runs Serial (the default Threaded refuses there, as for Pairwise Comparisons).
const { expect } = require('@playwright/test');
const { readFileSync } = require('node:fs');
const { evaluateWithRetainedPromise, settleLive } = require('./app-lifecycle.cjs');

/**
 * The bases of each record of a GenBank file, by its LOCUS name, upper case.
 * @param {string} path
 * @returns {Record<string, string>}
 */
const genbankBases = (path) => Object.fromEntries(readFileSync(path, 'utf8').split(/^\/\/\s*$/m)
  .filter((entry) => /^LOCUS/m.test(entry))
  .map((entry) => [entry.match(/^LOCUS\s+(\S+)/m)[1], entry.split(/^ORIGIN.*$/m)[1].replace(/[^a-z]/gi, '').toUpperCase()]));

const COMPLEMENT = { A: 'T', C: 'G', G: 'C', T: 'A', N: 'N' };
/** @param {string} bases */
const reverseComplement = (bases) => [...bases].reverse().map((base) => COMPLEMENT[base] || 'N').join('');

/**
 * Bases `start`..`end` (1-based, inclusive) of a record; `start > end` reads across the origin.
 * @param {string} bases @param {number} start @param {number} end
 */
const span = (bases, start, end) => (start <= end ? bases.slice(start - 1, end) : bases.slice(start - 1) + bases.slice(0, end));

/** @param {{ name: string, bases: string }[]} queries */
const fasta = (queries) => queries.map(({ name, bases }) => `>${name}\n${bases}\n`).join('');

/**
 * Searches from the panel as a user does: the section opened from the section index, the query pasted, the task and
 * the hit colour chosen, then Search. Resolves when the search has ended and anything it drew has settled.
 * @param {import('@playwright/test').Page} page
 * @param {string} text
 * @param {{ task?: string, colorMode?: 'query' | 'gradient' }} [options]
 */
const searchFromPanel = async (page, text, { task = 'megablast', colorMode = 'query' } = {}) => {
  await page.evaluate(() => { window.__GBDRAW_APP__.losatExecution.executionMode = 'serial'; });
  await page.getByRole('button', { name: 'Go to Sequence search', exact: true }).click();
  const unfold = page.locator('[data-sequence-search-unfold]');
  if (await unfold.count()) await unfold.click();
  await page.locator('#sequence-search-query').fill(text);
  await page.locator('#sequence-search-task').selectOption(task);
  await page.locator('#sequence-search-color-mode').selectOption(colorMode);
  await page.locator('[data-sequence-search-run]').click();
  await expect.poll(() => page.evaluate(() => window.__GBDRAW_APP__.sequenceSearch.panel.running), { timeout: 300_000 }).toBe(false);
  await settleLive(page);
};

/**
 * Searches through the panel's controller (`app.sequenceSearch`), for the cases that drive the app by its actions.
 * @param {import('@playwright/test').Page} page
 * @param {string} text
 * @param {{ task?: string, colorMode?: 'query' | 'gradient' }} [options]
 */
const searchByAction = async (page, text, { task = 'megablast', colorMode = 'query' } = {}) => {
  await evaluateWithRetainedPromise(page, async ({ text, task, colorMode }) => {
    window.__GBDRAW_APP__.losatExecution.executionMode = 'serial';
    const search = window.__GBDRAW_APP__.sequenceSearch;
    search.panel.task = task;
    search.panel.colorMode = colorMode;
    await search.readQueries(text);
    await search.search();
  }, { text, task, colorMode });
  await settleLive(page);
};

/**
 * The latest search of the shown drawing (`shownOutcome`) without its raw text, and the panel's error, as data.
 * @param {import('@playwright/test').Page} page
 */
const searchOutcome = (page) => page.evaluate(() => {
  const search = window.__GBDRAW_APP__.sequenceSearch;
  const outcome = search.shownOutcome ? JSON.parse(JSON.stringify(search.shownOutcome)) : null;
  return {
    outcome,
    searchError: search.panel.searchError ? JSON.parse(JSON.stringify(search.panel.searchError)) : null,
    noHitText: search.noHitText
  };
});

/**
 * The drawn hits of a set on the shown Result: annotation id -> the number of its drawn marks (paths and rects, the
 * gradient steps left out, as the outline reads them).
 * @param {import('@playwright/test').Page} page @param {string} setId
 */
const drawnHits = (page, setId) => page.evaluate((id) => {
  const svg = window.__GBDRAW_APP__.svgContainer?.querySelector('svg');
  return Object.fromEntries([...(svg?.querySelectorAll('g[data-gbdraw-annotation-id]') || [])]
    .filter((group) => group.getAttribute('data-gbdraw-annotation-set-id') === id)
    .map((group) => [group.getAttribute('data-gbdraw-annotation-id'), group.querySelectorAll('path, rect').length]));
}, setId);

/**
 * Every drawn annotation mark of the shown Result, by set and annotation id: its slot, and each mark's shape, geometry
 * and paint (a gradient read as its stop colours). The semantic snapshot of the parity check leaves geometry out; a
 * search hit's lane and its gradient are geometry and paint of these marks.
 * @param {import('@playwright/test').Page} page
 */
const drawnAnnotationMarks = (page) => page.evaluate(() => {
  const svg = window.__GBDRAW_APP__.svgContainer.querySelector('svg');
  const paint = (value) => {
    const reference = String(value || '').match(/^url\(\s*['"]?#([^'")]+)['"]?\s*\)/);
    if (!reference) return String(value || '').toLowerCase();
    const target = svg.querySelector(`#${CSS.escape(reference[1])}`);
    return `gradient(${[...(target?.querySelectorAll('stop') || [])]
      .map((stop) => `${stop.getAttribute('offset')} ${String(stop.getAttribute('stop-color')).toLowerCase()}`).join(', ')})`;
  };
  const geometry = ['d', 'x', 'y', 'width', 'height', 'transform', 'opacity', 'fill-opacity', 'stroke', 'stroke-width'];
  return Object.fromEntries([...svg.querySelectorAll('g[data-gbdraw-annotation-id]')].map((group) => [
    `${group.getAttribute('data-gbdraw-annotation-set-id')} / ${group.getAttribute('data-gbdraw-annotation-id')}`,
    {
      slot: group.closest('[data-gbdraw-slot-id]')?.getAttribute('data-gbdraw-slot-id') ?? null,
      marks: [...group.querySelectorAll('path, rect, polygon, line')].map((mark) => ({
        shape: mark.localName,
        fill: paint(mark.getAttribute('fill')),
        ...Object.fromEntries(geometry.filter((name) => mark.hasAttribute(name)).map((name) => [name, mark.getAttribute(name)]))
      }))
    }
  ]));
});

/** The undo labels of the shown drawing's History. @param {import('@playwright/test').Page} page */
const undoLabels = (page) => page.evaluate(() => window.__GBDRAW_HISTORY__.timeline().undo.map((step) => step.label));

/** The Generate counter: a live draw leaves it as it was. @param {import('@playwright/test').Page} page */
const generationKey = (page) => page.evaluate(async () => (await import('./js/state.js')).state.resultGenerationKey.value);

module.exports = {
  drawnAnnotationMarks, drawnHits, fasta, generationKey, genbankBases, reverseComplement, searchByAction, searchFromPanel, searchOutcome, span,
  undoLabels
};
