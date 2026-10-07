// Zero shift (Owner, 2026-10-07): the Web lays a Legend out exactly as Python
// does. The Legend layout owner (app/legend/layout-actions.js, through the
// port in services/legend-layout.js) is run on the Gallery Results Python drew,
// unedited and with the rows of tests/fixtures/legend_layout_vectors.json
// edited into the SVG, and every position, group offset, and Legend size must
// equal Python's (the vectors) exactly. The app cases then check that a live
// Legend edit and the next Generate draw the same Legend and canvas.
const { test, expect } = require('@playwright/test');
const { readFileSync } = require('node:fs');
const { resolve } = require('node:path');
const { evaluateWithRetainedPromise, generateAndWaitForResult, openApp } = require('./helpers/app-lifecycle.cjs');
const { openWithGenBank, openFresh } = require('./helpers/audit-browser.cjs');
const { settleLive } = require('./helpers/live-generate-parity.cjs');

test.describe.configure({ retries: 0 });

const VECTORS = JSON.parse(readFileSync(resolve(__dirname, '../fixtures/legend_layout_vectors.json'), 'utf8'));
const FIXTURE = 'tests/fixtures/forced_label_underlay.gb';

test('the Legend layout owner equals Python on every Gallery Legend and its edits', async ({ page }) => {
  test.setTimeout(300_000);
  await openApp(page);
  const cases = VECTORS.layout.filter((vector) => !vector.name.startsWith('synthetic'));
  const report = await evaluateWithRetainedPromise(page, async (layoutCases) => {
    const { createLegendLayoutActions } = await import('./js/app/legend/layout-actions.js');
    const actions = createLegendLayoutActions();
    await actions.prepareLegendLayout();
    const sessionText = async (base) => {
      for (const suffix of ['.gbdraw-session.json', '.gbdraw-session.json.gz']) {
        const response = await fetch(`./gallery/sessions/${base}${suffix}`);
        if (!response.ok) continue;
        if (!suffix.endsWith('.gz')) return response.text();
        const stream = response.body.pipeThrough(new DecompressionStream('gzip'));
        return new Response(stream).text();
      }
      throw new Error(`No Gallery Session ${base}`);
    };
    const sessions = new Map();
    const translateOf = (element) => {
      let x = 0;
      let y = 0;
      const pattern = /translate\(\s*([-+.\deE]+)(?:\s*[, ]\s*([-+.\deE]+))?\s*\)/g;
      for (const match of String(element?.getAttribute('transform') || '').matchAll(pattern)) {
        x += Number(match[1]);
        y += Number(match[2] || 0);
      }
      return [x, y];
    };
    const entriesOf = (group) => Array.from(group?.children || [])
      .filter((child) => child.localName === 'g' && child.hasAttribute('data-legend-key'));
    const swatchOf = (entry) => Array.from(entry.querySelectorAll('path'))
      .find((path) => { const fill = path.getAttribute('fill'); return fill && fill !== 'none' && !fill.startsWith('url('); });
    // Edit the solid rows of one feature group into `keys`, as the Legend
    // editor does: a renamed row keeps its node, a deleted one is removed, an
    // added one copies the first row, and the order follows `keys`.
    const editRows = (group, keys) => {
      const original = entriesOf(group);
      const byKey = new Map(original.map((entry) => [entry.getAttribute('data-legend-key'), entry]));
      const wanted = new Set(keys);
      const used = new Set();
      const next = keys.map((key, index) => {
        let entry = byKey.get(key);
        if (!entry) {
          const atIndex = original[index];
          if (atIndex && !wanted.has(atIndex.getAttribute('data-legend-key')) && !used.has(atIndex)) entry = atIndex;
          else entry = original[0].cloneNode(true);
          entry.setAttribute('data-legend-key', key);
          entry.querySelector('text').textContent = key;
        }
        used.add(entry);
        return entry;
      });
      original.forEach((entry) => { if (!used.has(entry)) entry.remove(); });
      next.forEach((entry) => group.appendChild(entry));
    };
    const mismatches = [];
    let numbers = 0;
    const same = (actual, expected, path) => {
      if (Array.isArray(expected)) { expected.forEach((value, index) => same(actual?.[index], value, `${path}[${index}]`)); return; }
      numbers += 1;
      if (!(actual === expected)) mismatches.push(`${path}: Web ${actual}, Python ${expected}`);
    };
    for (const vector of layoutCases) {
      const [base, index] = vector.name.split('#');
      if (!sessions.has(base)) sessions.set(base, JSON.parse(await sessionText(base)));
      const content = sessions.get(base).results[Number(index)].content;
      const svg = new DOMParser().parseFromString(content, 'image/svg+xml').documentElement;
      const legend = svg.getElementById('legend');
      const keys = vector.rows.filter((row) => row.type === 'solid').map((row) => row.key);
      const label = `${vector.name} / ${vector.edit}`;
      const expected = vector.expected;
      if (vector.mode === 'linear') {
        editRows(legend.querySelector('#feature_legend_h'), keys);
        editRows(legend.querySelector('#feature_legend_v'), keys);
      } else {
        editRows(legend.querySelector('#feature_legend') || legend, keys);
      }
      const box = actions.layOutLegend(svg, { side: vector.options.side });
      same([box.minX, box.minY, box.maxX, box.maxY], expected.localBounds, `${label} localBounds`);
      const checkEntries = (group, expectedEntries, path) => {
        const entries = entriesOf(group);
        same(entries.map((entry) => entry.getAttribute('data-legend-key')).length, expectedEntries.length, `${path} rows`);
        expectedEntries.forEach((row, rowIndex) => {
          const entry = entries[rowIndex];
          if (entry?.getAttribute('data-legend-key') !== row.key) mismatches.push(`${path}[${rowIndex}] key ${entry?.getAttribute('data-legend-key')} vs ${row.key}`);
          same(translateOf(swatchOf(entry)), [row.rectX, row.rectY], `${path}[${rowIndex}] swatch`);
          same(translateOf(entry.querySelector('text')), [row.textX, row.textY], `${path}[${rowIndex}] caption`);
          same(translateOf(entry), [0, 0], `${path}[${rowIndex}] entry group`);
        });
      };
      if (vector.mode === 'linear') {
        for (const [name, suffix] of [['horizontal', 'h'], ['vertical', 'v']]) {
          const orientation = expected[name];
          const group = legend.querySelector(`#legend_${name}`);
          checkEntries(group.querySelector(`#feature_legend_${suffix}`), orientation.feature.entries, `${label} ${name}`);
          same(translateOf(group.querySelector(`#feature_legend_${suffix}`)), [orientation.featureX, orientation.featureY], `${label} ${name} feature group`);
          const comparison = group.querySelector('[data-gbdraw-role="comparison-legend"]');
          if (orientation.gradient) same(translateOf(comparison), [orientation.gradientX, orientation.gradientY], `${label} ${name} comparison group`);
        }
      } else {
        checkEntries(legend.querySelector('#feature_legend') || legend, expected.entries, label);
        const conservation = legend.querySelector('#conservation_identity_legend');
        if (expected.gradient) same(translateOf(conservation), [expected.gradientX, expected.gradientY], `${label} conservation group`);
      }
    }
    return { cases: layoutCases.length, numbers, mismatches: mismatches.slice(0, 20), mismatchCount: mismatches.length };
  }, cases);
  expect(report.cases).toBeGreaterThanOrEqual(90);
  expect(report.mismatches, `${report.mismatchCount} mismatches of ${report.numbers}`).toEqual([]);
});

// The rows of the mounted Legend, the Legend position, and the canvas, as
// numbers (the writers format transforms differently).
const legendGeometry = (page) => page.evaluate(() => {
  const numbers = (value) => (String(value || '').match(/[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?/g) || []).map(Number);
  const svg = window.__GBDRAW_APP__.svgContainer.querySelector('svg');
  const legend = svg.getElementById('legend');
  const shown = Array.from(legend.querySelectorAll('g[data-legend-key]'))
    .filter((entry) => !entry.closest('[display="none"]') && !entry.closest('[data-gbdraw-role="comparison-legend"]'));
  return {
    canvas: [...numbers(svg.getAttribute('width')), ...numbers(svg.getAttribute('height')), ...numbers(svg.getAttribute('viewBox'))],
    legend: numbers(legend.getAttribute('transform')),
    rows: shown.map((entry) => [
      entry.getAttribute('data-legend-key'),
      numbers(entry.querySelector('text')?.getAttribute('transform')),
      Array.from(entry.querySelectorAll('path')).flatMap((path) => numbers(path.getAttribute('transform')))
    ])
  };
});

const openLinear = async (page, fixture = FIXTURE, form = {}) => {
  await openFresh(page);
  await page.getByRole('button', { name: 'Linear', exact: true }).click();
  await page.waitForFunction(() => window.__GBDRAW_APP__?.mode === 'linear');
  await page.evaluate(async ({ text, name, values }) => {
    const app = window.__GBDRAW_APP__;
    Object.assign(app.form, values);
    app.setLinearSeqPrimaryFile(0, 'gb', new File([text], name, { type: 'text/plain', lastModified: 1000 }));
    await window.Vue.nextTick();
  }, { text: readFileSync(fixture, 'utf8'), name: fixture.split('/').pop(), values: form });
  await settleLive(page);
};

const generate = async (page) => { await generateAndWaitForResult(page); await settleLive(page); };
const rowIndex = (page, pattern) => page.evaluate((source) => (
  window.__GBDRAW_APP__.legendEntries.findIndex((entry) => new RegExp(source).test(entry.caption))
), pattern);

for (const mode of ['linear', 'circular']) {
  test(`live Legend edits equal the next Generate exactly (${mode})`, async ({ page }) => {
    test.setTimeout(300_000);
    if (mode === 'linear') await openLinear(page);
    else await openWithGenBank(page, FIXTURE, () => { window.__GBDRAW_APP__.form.labels_mode = 'out'; });
    await generate(page);
    const unedited = await legendGeometry(page);
    const steps = [
      ['add a row', async () => evaluateWithRetainedPromise(page, async () => {
        const app = window.__GBDRAW_APP__;
        app.newLegendCaption = 'Manual row with a long caption';
        app.newLegendColor = '#118833';
        await app.addNewLegendEntry();
      })],
      ['delete the first row', async () => {
        const index = await rowIndex(page, '.');
        await page.evaluate((i) => window.__GBDRAW_APP__.deleteLegendEntry(i), index);
      }],
      ['sort the rows', async () => page.evaluate(() => window.__GBDRAW_APP__.sortLegendEntries('desc'))]
    ];
    if (mode === 'circular') {
      steps.push(['rename GC content', async () => evaluateWithRetainedPromise(page, async () => {
        const app = window.__GBDRAW_APP__;
        const index = app.legendEntries.findIndex((entry) => /^GC content/.test(entry.caption));
        await app.renameLegendEntry(index, 'A much longer Legend caption for GC');
      })]);
    }
    for (const [name, run] of steps) {
      await run();
      await page.waitForTimeout(150);
      await settleLive(page);
      const live = await legendGeometry(page);
      await generate(page);
      const generated = await legendGeometry(page);
      expect(generated, `${mode}: ${name}`).toEqual(live);
    }
    expect(unedited.rows.length).toBeGreaterThan(0);
  });
}

// OV-156: in Linear, a row without features renamed in the Legend editor keeps
// its place in the Legend, live and at Generate, also when it is not the last
// row. HmmtDNA with GC content and GC skew; `GC skew (+)` is renamed (the M9
// rename), then renamed onto `GC skew (-)` with Suffix (the Linear half of the
// OV-151 M10 rename).
test('a renamed Linear row that is not last keeps its place live and at Generate (OV-156)', async ({ page }) => {
  test.setTimeout(300_000);
  await openLinear(page, 'tests/test_inputs/HmmtDNA.gbk', { show_gc: true, show_skew: true });
  await generate(page);
  const captions = () => page.evaluate(() => window.__GBDRAW_APP__.legendEntries.map((entry) => entry.caption));
  const drawnKeys = (geometry) => geometry.rows.map(([key]) => key);
  const before = await captions();
  expect(before).toContain('GC skew (+)');
  expect(before.indexOf('GC skew (+)'), 'GC skew (+) is not the last row').toBeLessThan(before.length - 1);
  const rename = (from, to) => evaluateWithRetainedPromise(page, async ({ from: source, to: target }) => {
    const app = window.__GBDRAW_APP__;
    await app.renameLegendEntry(app.legendEntries.findIndex((entry) => entry.caption === source), target);
  }, { from, to });
  const steps = [
    ['rename GC skew (+)', 'skew plus', async () => rename('GC skew (+)', 'skew plus')],
    ['rename it onto GC skew (-) with Suffix', 'GC skew (-) (1)', async () => {
      await rename('skew plus', 'GC skew (-)');
      expect(await page.evaluate(() => window.__GBDRAW_APP__.legendRenameDialog.show), 'the rename asks').toBe(true);
      await evaluateWithRetainedPromise(page, async () => { await window.__GBDRAW_APP__.handleLegendRenameChoice('suffix'); });
    }]
  ];
  for (const [name, caption, run] of steps) {
    await run();
    await page.waitForTimeout(150);
    await settleLive(page);
    const expected = before.map((entry) => (entry === 'GC skew (+)' ? caption : entry));
    expect(await captions(), `${name}: the editor keeps the order`).toEqual(expected);
    const live = await legendGeometry(page);
    await generate(page);
    expect(await captions(), `${name}: Generate keeps the editor order`).toEqual(expected);
    const generated = await legendGeometry(page);
    // The shown orientation draws each caption once, in the editor's order.
    expect(drawnKeys(generated), `${name}: drawn rows`).toEqual(expected.filter((entry) => drawnKeys(generated).includes(entry)));
    expect(new Set(drawnKeys(generated)).size, `${name}: one row per caption`).toBe(drawnKeys(generated).length);
    expect(generated, `${name}: live equals Generate`).toEqual(live);
  }
});

// OV-143 (E1 REVIEW-3, controls K1-K8): single-mode Legend edits whose screen
// differed from the next Generate on dev. Each shown Legend equals the next
// Generate exactly, and no two rows overlap.
const rowOverlaps = (page) => page.evaluate(() => {
  const svg = window.__GBDRAW_APP__.svgContainer.querySelector('svg');
  const shown = Array.from(svg.getElementById('legend').querySelectorAll('g[data-legend-key]'))
    .filter((entry) => !entry.closest('[display="none"]') && !entry.closest('[data-gbdraw-role="comparison-legend"]'))
    .map((entry) => [entry.getAttribute('data-legend-key'), entry.getBoundingClientRect()]);
  const overlaps = [];
  shown.forEach(([first, a], index) => shown.slice(index + 1).forEach(([second, b]) => {
    const width = Math.min(a.right, b.right) - Math.max(a.left, b.left);
    const height = Math.min(a.bottom, b.bottom) - Math.max(a.top, b.top);
    if (width > 0.5 && height > 0.5) overlaps.push(`${first} / ${second}`);
  }));
  return overlaps;
});
const expectShownEqualsGenerate = async (page, label) => {
  await settleLive(page);
  const shown = await legendGeometry(page);
  expect(await rowOverlaps(page), `${label}: shown rows overlap`).toEqual([]);
  await generate(page);
  expect(await rowOverlaps(page), `${label}: generated rows overlap`).toEqual([]);
  expect(await legendGeometry(page), `${label}: shown equals Generate`).toEqual(shown);
  return shown;
};
const openControlCircular = async (page) => {
  await openWithGenBank(page, FIXTURE, () => { window.__GBDRAW_APP__.autoLabelReflowEnabled = false; });
  await generate(page);
};
const renameRow = async (page, from, to) => {
  await evaluateWithRetainedPromise(page, async ({ source, target }) => {
    const app = window.__GBDRAW_APP__;
    await app.renameLegendEntry(app.legendEntries.findIndex((entry) => entry.caption === source), target);
  }, { source: from, target: to });
  await settleLive(page);
};
// A delete recorded as one History step, as the Legend editor records it.
const recordedDelete = async (page, caption) => {
  await evaluateWithRetainedPromise(page, async (target) => {
    const history = window.__GBDRAW_HISTORY__;
    const transaction = await history.begin('Change setting', { source: 'legend-zero-shift', owner: {} });
    const app = window.__GBDRAW_APP__;
    await app.deleteLegendEntry(app.legendEntries.findIndex((entry) => entry.caption === target));
    await window.Vue.nextTick();
    if (transaction) await history.commit(transaction);
  }, caption);
  await settleLive(page);
};
const historyStep = async (page, name) => {
  await evaluateWithRetainedPromise(page, async (step) => { await window.__GBDRAW_HISTORY__[step](); }, name);
  await settleLive(page);
};

test('OV-143 K1: a GC content row renamed in Circular stays in its slot', async ({ page }) => {
  test.setTimeout(300_000);
  await openControlCircular(page);
  const before = await legendGeometry(page);
  await renameRow(page, 'GC content', 'GC%');
  const shown = await expectShownEqualsGenerate(page, 'K1');
  const row = (geometry, key) => geometry.rows.find(([caption]) => caption === key);
  expect(row(shown, 'GC%')[1][1], 'K1: the row keeps its y').toBe(row(before, 'GC content')[1][1]);
});

test('OV-143 K4: an added row and a renamed CDS row in Circular do not overlap', async ({ page }) => {
  test.setTimeout(300_000);
  await openControlCircular(page);
  await evaluateWithRetainedPromise(page, async () => {
    const app = window.__GBDRAW_APP__;
    app.newLegendCaption = 'Extra';
    app.newLegendColor = '#7b2cbf';
    await app.addNewLegendEntry();
  });
  await page.waitForFunction(() => window.__GBDRAW_APP__.legendEntries.some((entry) => entry.caption === 'Extra'));
  await renameRow(page, 'CDS', 'Coding');
  await expectShownEqualsGenerate(page, 'K4');
});

test('OV-143 K6: a row moved up in Linear takes its place in the row', async ({ page }) => {
  test.setTimeout(300_000);
  await openLinear(page);
  await generate(page);
  await page.evaluate(() => {
    const app = window.__GBDRAW_APP__;
    app.moveLegendEntryUp(app.legendEntries.findIndex((entry) => entry.caption === 'repeat_region'));
  });
  await expectShownEqualsGenerate(page, 'K6');
});

test('OV-143 K7: two deletes in Circular are drawn compact live and at Generate', async ({ page }) => {
  test.setTimeout(300_000);
  await openControlCircular(page);
  await recordedDelete(page, 'repeat_region');
  await recordedDelete(page, 'GC content');
  await expectShownEqualsGenerate(page, 'K7');
});

test('OV-143 K8: Undo of a delete returns the row to its slot; Generate agrees', async ({ page }) => {
  test.setTimeout(300_000);
  await openControlCircular(page);
  await recordedDelete(page, 'repeat_region');
  await recordedDelete(page, 'GC content');
  await historyStep(page, 'undo');
  await expectShownEqualsGenerate(page, 'K8 Undo of the second delete');
});

test('OV-143 K8: Undo and Redo of two deletes walk the laid-out Legends', async ({ page }) => {
  test.setTimeout(300_000);
  await openControlCircular(page);
  const unedited = await legendGeometry(page);
  await recordedDelete(page, 'repeat_region');
  const oneDeleted = await legendGeometry(page);
  await recordedDelete(page, 'GC content');
  const twoDeleted = await legendGeometry(page);
  const walk = [['undo', oneDeleted], ['undo', unedited], ['redo', oneDeleted], ['redo', twoDeleted]];
  for (const [step, expected] of walk) {
    await historyStep(page, step);
    expect(await rowOverlaps(page), `K8 ${step}: rows overlap`).toEqual([]);
    expect(await legendGeometry(page), `K8 ${step}`).toEqual(expected);
  }
  await expectShownEqualsGenerate(page, 'K8 after Redo');
});
