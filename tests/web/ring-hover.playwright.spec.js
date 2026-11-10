// The hover on similarity rings in the preview (PD-OI-110 scenario revision 3, gbdraw/web/js/app/feature-editor/
// ring-hover.js): over a ring, that ring keeps its strength and the others are drawn at 35%, its Legend row is bold,
// a radial line marks the pointer, and a popup lists the position (bp) and every ring's identity there, outermost
// first; over a Legend name, the same without line and popup; leaving restores everything. On the Multi-Record Canvas
// the series keeps its strength in every record. Nothing reaches the Result: its SVG download and Save text are the
// same with and without a hover.
const { test, expect } = require('@playwright/test');
const { readFileSync } = require('node:fs');
const path = require('node:path');
const { generateAndWaitForResult } = require('./helpers/app-lifecycle.cjs');
const { openWithGenBank, settle } = require('./helpers/audit-browser.cjs');

const root = path.resolve(process.env.GBDRAW_REPO || process.cwd());

/** outfmt 6 rows of one ring: `[identity, start, end]` on each named record (query and subject spans alike). */
const ringTable = (records, hits) => records.flatMap((record) => hits.map(([identity, start, end]) => [
  'RING', record, identity.toFixed(1), String(end - start + 1), '0', '0', String(start), String(end), String(start), String(end), '1e-50', '200'
].join('\t'))).join('\n') + '\n';

const uploadRings = (page, tables) => page.evaluate(async (texts) => {
  const app = window.__GBDRAW_APP__;
  app.circularConservation.source = 'upload';
  app.setCircularConservationUploadFiles(texts.map(({ name, text }) => new File([text], name, { type: 'text/plain' })));
  await window.Vue.nextTick();
}, tables);

// The screen point at a display position of a ring (by series) of a record, in the middle of the ring's band.
const ringPoint = (page, { sourceIndex, recordIndex = 0, position }) => page.evaluate(async ({ sourceIndex: series, recordIndex: record, position: at }) => {
  const { state } = await import('./js/state.js');
  const { regionAxisFrame, shownAxisRecords } = await import('./js/services/region-axis.js');
  const svg = window.__GBDRAW_APP__.svgContainer.querySelector('svg');
  const { axis } = shownAxisRecords(state).find((entry) => entry.recordIndex === record);
  const frame = regionAxisFrame(svg, record);
  const ring = svg.querySelector(`g#record_${record} g[data-gbdraw-slot-renderer="sequence_conservation"][data-source-index="${series}"]`);
  const radius = (Number(ring.getAttribute('data-ring-inner-radius')) + Number(ring.getAttribute('data-ring-outer-radius'))) / 2;
  const angle = (axis.startAngleDeg + axis.degreesPerBp * at) * Math.PI / 180;
  const point = new DOMPoint(radius * Math.cos(angle), radius * Math.sin(angle)).matrixTransform(frame.getScreenCTM());
  const hit = document.elementFromPoint(point.x, point.y);
  return { x: point.x, y: point.y, onRing: Boolean(hit && ring.contains(hit)) };
}, { sourceIndex, recordIndex, position });

// What the hover shows: each ring group's opacity by record and series, the bold Legend rows, the overlay.
const hoverView = (page) => page.evaluate(() => {
  const svg = window.__GBDRAW_APP__.svgContainer.querySelector('svg');
  const rings = [...svg.querySelectorAll('g[data-gbdraw-slot-renderer="sequence_conservation"]')].map((group) => ({
    record: group.closest('[data-gbdraw-record-index]')?.getAttribute('data-gbdraw-record-index') ?? '0',
    series: group.getAttribute('data-source-index'),
    opacity: getComputedStyle(group).opacity
  }));
  const rows = [...svg.querySelectorAll('#conservation_identity_legend g[data-source-index]')].map((row) => ({
    series: row.getAttribute('data-source-index'),
    name: row.getAttribute('data-legend-key'),
    weight: getComputedStyle(row.querySelector('text')).fontWeight
  }));
  const overlay = document.querySelector('.gbdraw-ring-hover-overlay');
  const shown = Boolean(overlay && !overlay.hidden);
  return {
    rings,
    rows,
    shown,
    position: shown ? overlay.querySelector('.gbdraw-ring-hover-position').textContent : '',
    popup: shown ? [...overlay.querySelectorAll('.gbdraw-ring-hover-row')].map((row) => ({
      name: row.querySelector('.gbdraw-ring-hover-name').textContent,
      value: row.querySelector('.gbdraw-ring-hover-value').textContent,
      swatch: getComputedStyle(row.querySelector('.gbdraw-ring-hover-swatch')).backgroundColor,
      hovered: row.hasAttribute('data-hovered')
    })) : [],
    line: shown ? Number(overlay.querySelector('line').getAttribute('x2')) !== 0 || Number(overlay.querySelector('line').getAttribute('y2')) !== 0 : false
  };
});

// The SVG download (Export reads the mounted root) and the Results that Save writes.
const exported = async (page) => {
  const [download] = await Promise.all([page.waitForEvent('download'), page.evaluate(() => window.__GBDRAW_APP__.downloadSVG())]);
  const text = readFileSync(await download.path(), 'utf8');
  const results = await page.evaluate(() => window.__GBDRAW_APP__.currentResults().map((result) => String(result.content || '')));
  return { text, results };
};

// The base the popup names: the pointer lands on a whole pixel, a few bases from the asked position on a 3 kb ring.
const shownBase = (position) => Number(position.replace(/[^0-9]/g, ''));

// The whole figure in the preview, so every ring point is on screen.
const fitPreview = async (page) => {
  await page.getByRole('button', { name: 'Fit to preview' }).click();
  await settle(page);
};

const expectHover = async (page, series, { popup }) => {
  await expect.poll(async () => (await hoverView(page)).rings.every((ring) => ring.opacity === (ring.series === series ? '1' : '0.35'))).toBe(true);
  const view = await hoverView(page);
  for (const row of view.rows) expect(row.weight, `Legend row ${row.name}`).toBe(row.series === series ? '700' : '400');
  expect(view.shown).toBe(popup);
  return view;
};

const expectRestored = async (page) => {
  await expect.poll(async () => {
    const view = await hoverView(page);
    return view.rings.every((ring) => ring.opacity === '1') && view.rows.every((row) => row.weight === '400') && !view.shown;
  }).toBe(true);
};

test('a ring hover dims the other rings, bolds its Legend row, and lists every ring\'s identity; a Legend name does the same without popup', async ({ page }) => {
  test.setTimeout(420_000);
  await openWithGenBank(page, path.join(root, 'tests/fixtures/forced_label_underlay.gb'));
  await uploadRings(page, [
    { name: 'outer.blast', text: ringTable(['FORCEDLBL'], [[95, 1, 1000], [80, 1501, 2500]]) },
    { name: 'inner.blast', text: ringTable(['FORCEDLBL'], [[99, 1, 600], [90, 2001, 3000]]) }
  ]);
  await settle(page);
  await generateAndWaitForResult(page);
  await settle(page);
  await fitPreview(page);
  const before = await exported(page);

  // The inner ring at base 301 (display position 300): both rings have a hit there.
  const inner = await page.evaluate(() => [...window.__GBDRAW_APP__.svgContainer.querySelectorAll('g[data-gbdraw-slot-renderer="sequence_conservation"]')]
    .sort((a, b) => Number(a.getAttribute('data-ring-outer-radius')) - Number(b.getAttribute('data-ring-outer-radius')))[0].getAttribute('data-source-index'));
  const point = await ringPoint(page, { sourceIndex: inner, position: 300.5 });
  expect(point.onRing).toBe(true);
  await page.mouse.move(point.x, point.y);
  const view = await expectHover(page, inner, { popup: true });
  expect(view.position).toMatch(/^[0-9,]+ bp$/);
  expect(Math.abs(shownBase(view.position) - 301)).toBeLessThanOrEqual(5);
  expect(view.line).toBe(true);
  // Outermost first, named as the Legend names them (its rows run outer to inner).
  expect(view.popup.map((row) => row.name)).toEqual(view.rows.map((row) => row.name));
  expect(view.popup.map((row) => row.value)).toEqual(inner === '1' ? ['95.0%', '99.0%'] : ['99.0%', '95.0%']);
  expect(view.popup.filter((row) => row.hovered).length).toBe(1);
  // Each swatch is the colour the ring is drawn in at the position: a lower identity a lighter one.
  expect(view.popup.every((row) => row.swatch !== 'rgba(0, 0, 0, 0)')).toBe(true);
  expect(new Set(view.popup.map((row) => row.swatch)).size).toBe(2);
  // Nothing reaches the Result: the download and the Results Save writes are those without a hover.
  const during = await exported(page);
  expect(during.text).toBe(before.text);
  expect(during.results).toEqual(before.results);

  // Where the inner ring has a hit and the outer has none, the popup says so.
  const gap = await ringPoint(page, { sourceIndex: inner, position: 2200.5 });
  await page.mouse.move(gap.x, gap.y);
  await expect.poll(async () => Math.abs(shownBase((await hoverView(page)).position) - 2201)).toBeLessThanOrEqual(5);
  const values = Object.fromEntries((await hoverView(page)).popup.map((row) => [row.name, row.value]));
  expect(values).toEqual({ outer: '80.0%', inner: '90.0%' });

  // Leaving the drawing restores everything.
  await page.mouse.move(2, 2);
  await expectRestored(page);

  // A Legend name: the same rings and bold row, no popup and no line.
  const name = await page.evaluate((series) => {
    const row = window.__GBDRAW_APP__.svgContainer.querySelector(`#conservation_identity_legend g[data-source-index="${series}"] text`);
    const box = row.getBoundingClientRect();
    return { x: box.left + box.width / 2, y: box.top + box.height / 2 };
  }, inner);
  await page.mouse.move(name.x, name.y);
  await expectHover(page, inner, { popup: false });
  await page.mouse.move(2, 2);
  await expectRestored(page);
  expect((await exported(page)).text).toBe(before.text);
});

test('on the Multi-Record Canvas the hovered series keeps its strength in every record', async ({ page }) => {
  test.setTimeout(420_000);
  await openWithGenBank(page, path.join(root, 'tests/fixtures/web_batch_two_records.gb'), () => {
    window.__GBDRAW_APP__.form.multi_record_canvas = true;
  });
  await uploadRings(page, [
    { name: 'first.blast', text: ringTable(['TESTA', 'TESTB'], [[95, 1, 2000]]) },
    { name: 'second.blast', text: ringTable(['TESTB'], [[70, 1, 1000]]) }
  ]);
  await settle(page);
  await generateAndWaitForResult(page);
  await settle(page);
  await fitPreview(page);
  expect(await page.evaluate(() => window.__GBDRAW_APP__.svgContainer.querySelectorAll('svg [data-gbdraw-record-index="1"] g[data-gbdraw-slot-renderer="sequence_conservation"]').length)).toBe(2);

  const point = await ringPoint(page, { sourceIndex: '0', recordIndex: 1, position: 500.5 });
  expect(point.onRing).toBe(true);
  await page.mouse.move(point.x, point.y);
  const view = await expectHover(page, '0', { popup: true });
  // Both records keep series 0 and dim series 1.
  expect(view.rings.map(({ record, series, opacity }) => `${record}:${series}:${opacity}`).sort())
    .toEqual(['0:0:1', '0:1:0.35', '1:0:1', '1:1:0.35']);
  expect(Math.abs(shownBase(view.position) - 501)).toBeLessThanOrEqual(10);
  expect(Object.fromEntries(view.popup.map((row) => [row.name, row.value]))).toEqual({ first: '95.0%', second: '70.0%' });
  await page.mouse.move(2, 2);
  await expectRestored(page);
});
