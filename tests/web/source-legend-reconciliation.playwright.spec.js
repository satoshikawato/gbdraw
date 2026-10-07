const { test, expect } = require('@playwright/test');
const fs = require('node:fs/promises');
const { load, generate, switchMode, download } = require('./helpers/mode-transition.cjs');
const { evaluateWithRetainedPromise } = require('./helpers/app-lifecycle.cjs');

const seeds = {
  lambda: 'gbdraw/web/gallery/sessions/lambda_basic_linear.gbdraw-session.json',
  tobacco: 'gbdraw/web/gallery/sessions/tobacco-chloroplast.gbdraw-session.json'
};
const upload = async (page, name, linear = true) => {
  const session = JSON.parse(await fs.readFile(seeds[name], 'utf8'));
  const input = linear ? page.getByTestId('linear-genbank-1') : page.getByLabel('GenBank/DDBJ File', { exact: true });
  await input.setInputFiles({ name: `${name}.gb`, mimeType: 'text/plain', buffer: Buffer.from(
    session.resources[session.renderRequest.records[0].source.resourceId].data, 'base64') });
};
const inspect = page => page.evaluate(async () => {
  const { state: s } = await import('./js/state.js');
  const { getCommittedCanonicalRenderRequest } = await import('./js/services/config.js');
  const { getVisibleFeatureLegendGroup } = await import('./js/app/legend/utils.js');
  const digest = async value => [...new Uint8Array(await crypto.subtle.digest('SHA-256',
    new TextEncoder().encode(typeof value === 'string' ? value : JSON.stringify(value))))]
    .map(byte => byte.toString(16).padStart(2, '0')).join('');
  const captions = svg => [...(getVisibleFeatureLegendGroup(svg)?.querySelectorAll('g[data-legend-key]') || [])]
    .map(e => ({ caption: e.getAttribute('data-legend-key'), color: e.querySelector('path[fill]')?.getAttribute('fill') }))
    .sort((a, b) => a.caption.localeCompare(b.caption));
  const result = s.results.value[s.selectedResultIndex.value].content;
  const svg = new DOMParser().parseFromString(result, 'image/svg+xml').documentElement;
  return JSON.parse(JSON.stringify({
    entries: s.legendEntries.value.map(e => e.caption).sort(), original: s.originalLegendOrder.value,
    entryState: s.legendEntries.value, originalColors: s.originalLegendColors.value,
    manualRules: s.manualSpecificRules, strokes: s.legendStrokeOverrides,
    selectedResult: s.selectedResultIndex.value,
    resultsDigest: await digest(s.results.value),
    mountedDigest: await digest(s.svgContainer.value.querySelector('svg').outerHTML),
    requestDigest: await digest(getCommittedCanonicalRenderRequest()),
    result: captions(svg), mounted: captions(s.svgContainer.value.querySelector('svg')),
    side: s.form.legend, fontSize: s.adv.legend_font_size, preferences: s.layoutPreferences.legend,
    layoutPreferences: s.layoutPreferences,
    palette: s.selectedPalette.value, colors: s.currentColors.value,
    featureIds: s.extractedFeatures.value.map(f => f.stable_feature_id),
    overrides: s.legendColorOverrides, deleted: s.deletedLegendEntries.value
  }));
});
const expectEntries = async (page, expected) => {
  const current = await inspect(page);
  expect(current.entries).toEqual([...expected].sort());
  expect(current.result.map(e => e.caption).sort()).toEqual([...expected].sort());
  expect(current.mounted).toEqual(current.result);
  return current;
};
const reveal = async locator => {
  for (const details of await locator.locator('xpath=ancestor::details').all()) {
    if (await details.getAttribute('open') === null) await details.locator(':scope > summary').press('Enter');
  }
  return locator;
};

// The swatch of row `caption` in each Legend group of the selected Result and of the mounted SVG.
const rowStyle = (page, caption) => page.evaluate(async target => {
  const { state: s } = await import('./js/state.js');
  const { getAllFeatureLegendGroups, getLegendEntrySwatch } = await import('./js/app/legend/utils.js');
  const styles = svg => getAllFeatureLegendGroups(svg).map(group => {
    const swatch = getLegendEntrySwatch(group.querySelector(`g[data-legend-key="${CSS.escape(target)}"]`));
    return swatch && [swatch.getAttribute('fill'), swatch.getAttribute('stroke'), swatch.getAttribute('stroke-width')];
  });
  const content = s.results.value[s.selectedResultIndex.value].content;
  return {
    result: styles(new DOMParser().parseFromString(content, 'image/svg+xml').documentElement),
    mounted: styles(s.svgContainer.value.querySelector('svg'))
  };
}, caption);
const legendIndex = caption => window.__GBDRAW_APP__.legendEntries.findIndex(e => e.caption === caption);
const expectRow = async (page, caption, style) => {
  const current = await rowStyle(page, caption);
  expect(current.result.length).toBeGreaterThan(0);
  expect(current.result).toEqual(current.result.map(() => style));
  expect(current.mounted).toEqual(current.result);
};

// OV-86: a style on a row the Legend editor added failed every later Generate.
for (const mode of ['linear', 'circular']) {
  test(`M1 ${mode}: a style on a Legend editor added row survives Generate, rename, Session replay, and removal`, async ({ browser }, testInfo) => {
    test.setTimeout(900_000);
    const page = await load(browser);
    page.setDefaultTimeout(180_000);
    let fresh;
    try {
      if (mode === 'linear') await switchMode(page, 'linear');
      const input = mode === 'linear'
        ? page.getByTestId('linear-genbank-1')
        : page.getByLabel('GenBank/DDBJ File', { exact: true });
      await input.setInputFiles('tests/fixtures/forced_label_underlay.gb');
      await generate(page);
      await page.evaluate(async () => {
        const app = window.__GBDRAW_APP__;
        app.newLegendCaption = 'Manual row';
        app.newLegendColor = '#118833';
        await app.addNewLegendEntry();
      });
      await expect.poll(() => page.evaluate(legendIndex, 'Manual row')).toBeGreaterThanOrEqual(0);
      await generate(page);
      // Linear styles the fill first and Circular the stroke; each alone failed.
      const edits = {
        fill: () => page.evaluate(caption => {
          const app = window.__GBDRAW_APP__;
          return app.updateLegendEntryColor(app.legendEntries.findIndex(e => e.caption === caption), '#7b2cbf');
        }, 'Manual row'),
        stroke: () => page.evaluate(caption => {
          const app = window.__GBDRAW_APP__;
          const index = app.legendEntries.findIndex(e => e.caption === caption);
          return app.updateLegendEntryStrokeColor(index, '#e63946') && app.updateLegendEntryStrokeWidth(index, 2);
        }, 'Manual row')
      };
      for (const edit of mode === 'linear' ? ['fill', 'stroke'] : ['stroke', 'fill']) {
        expect(await edits[edit]()).toBe(true);
        await generate(page);
      }
      const styled = ['#7b2cbf', '#e63946', '2'];
      await expectRow(page, 'Manual row', styled);
      await evaluateWithRetainedPromise(page, async index => {
        await window.__GBDRAW_APP__.renameLegendEntry(index, 'Manual renamed');
      }, await page.evaluate(legendIndex, 'Manual row'));
      await expect.poll(() => page.evaluate(legendIndex, 'Manual renamed')).toBeGreaterThanOrEqual(0);
      await generate(page);
      await expectRow(page, 'Manual renamed', styled);
      await expectRow(page, 'Manual row', null);
      const saved = testInfo.outputPath(`added-row-style-${mode}.gbdraw-session.json.gz`);
      await download(page, 'Save Session', saved);
      fresh = await load(browser, saved);
      await expectRow(fresh, 'Manual renamed', styled);
      await generate(fresh);
      await expectRow(fresh, 'Manual renamed', styled);
      await fresh.evaluate(caption => {
        const app = window.__GBDRAW_APP__;
        app.deleteLegendEntry(app.legendEntries.findIndex(e => e.caption === caption));
      }, 'Manual renamed');
      await generate(fresh);
      await expectRow(fresh, 'Manual renamed', null);
      expect(page.externalRequests).toEqual([]);
      expect(fresh.externalRequests).toEqual([]);
    } finally {
      await page.context().close();
      if (fresh) await fresh.context().close();
    }
  });
}

// `forced_label_underlay.gb` generated in `mode`.
const loadGenerated = async (browser, mode) => {
  const page = await load(browser);
  page.setDefaultTimeout(180_000);
  if (mode === 'linear') await switchMode(page, 'linear');
  const input = mode === 'linear'
    ? page.getByTestId('linear-genbank-1')
    : page.getByLabel('GenBank/DDBJ File', { exact: true });
  await input.setInputFiles('tests/fixtures/forced_label_underlay.gb');
  await generate(page);
  return page;
};
const addLegendRow = async (page, caption, color) => {
  await evaluateWithRetainedPromise(page, async row => {
    const app = window.__GBDRAW_APP__;
    app.newLegendCaption = row.caption;
    app.newLegendColor = row.color;
    await app.addNewLegendEntry();
  }, { caption, color });
  await expect.poll(() => page.evaluate(legendIndex, caption)).toBeGreaterThanOrEqual(0);
};
// The swatch stroke of row `caption` (color normalized, width as a number) in each
// Legend group of the selected Result and of the mounted SVG.
const rowStroke = (page, caption) => page.evaluate(async target => {
  const { state: s } = await import('./js/state.js');
  const { getAllFeatureLegendGroups, getLegendEntrySwatch } = await import('./js/app/legend/utils.js');
  const paint = document.createElement('canvas').getContext('2d');
  const color = value => { paint.fillStyle = '#010203'; paint.fillStyle = String(value); return String(paint.fillStyle); };
  const strokes = svg => getAllFeatureLegendGroups(svg).map(group => {
    const swatch = getLegendEntrySwatch(group.querySelector(`g[data-legend-key="${CSS.escape(target)}"]`));
    return swatch && [color(swatch.getAttribute('stroke')), Number(swatch.getAttribute('stroke-width'))];
  });
  const content = s.results.value[s.selectedResultIndex.value].content;
  return {
    result: strokes(new DOMParser().parseFromString(content, 'image/svg+xml').documentElement),
    mounted: strokes(s.svgContainer.value.querySelector('svg'))
  };
}, caption);
const expectStroke = async (page, caption, stroke) => {
  const current = await rowStroke(page, caption);
  expect(current.result.length).toBeGreaterThan(0);
  expect(current.result).toEqual(current.result.map(() => stroke));
  expect(current.mounted).toEqual(current.result);
};

// OV-121 (PD-OI-066): a row added in the Legend editor without a stroke of its own
// takes the renderer's first row as Generate copies it, before that row's own
// stroke edit, live and at Generate; a later edit of the first row leaves it.
for (const mode of ['linear', 'circular']) {
  test(`M2 ${mode}: a Legend editor added row takes the first row's drawn stroke live and at Generate`, async ({ browser }) => {
    test.setTimeout(600_000);
    const page = await loadGenerated(browser, mode);
    try {
      const first = await page.evaluate(() => window.__GBDRAW_APP__.legendEntries[0].caption);
      const [drawn] = (await rowStroke(page, first)).result;
      const strokeFirstRow = (color, width) => page.evaluate(edit => {
        const app = window.__GBDRAW_APP__;
        return app.updateLegendEntryStrokeColor(0, edit.color) && app.updateLegendEntryStrokeWidth(0, edit.width);
      }, { color, width });
      expect(await strokeFirstRow('#e63946', 3)).toBe(true);
      await addLegendRow(page, 'Manual row', '#118833');
      await expectStroke(page, 'Manual row', drawn);
      await generate(page);
      await expectStroke(page, first, ['#e63946', 3]);
      await expectStroke(page, 'Manual row', drawn);
      expect(await strokeFirstRow('#2a9d8f', 1)).toBe(true);
      await expectStroke(page, 'Manual row', drawn);
      await generate(page);
      await expectStroke(page, first, ['#2a9d8f', 1]);
      await expectStroke(page, 'Manual row', drawn);
      expect(page.externalRequests).toEqual([]);
    } finally {
      await page.context().close();
    }
  });
}

// The caption box [left, top, right, bottom] of each drawn Legend row of the
// mounted SVG, and the canvas [width, height], in canvas units.
const legendCaptionBoxes = page => page.evaluate(async () => {
  const { state: s } = await import('./js/state.js');
  const svg = s.svgContainer.value.querySelector('svg');
  const width = svg.viewBox.baseVal?.width || svg.width.baseVal.value;
  const height = svg.viewBox.baseVal?.height || svg.height.baseVal.value;
  const frame = svg.getBoundingClientRect();
  const scale = frame.width / width;
  const rows = [...svg.querySelectorAll('#legend g[data-legend-key]')]
    .filter(entry => !entry.closest('[display="none"]'))
    .map(entry => {
      const box = entry.querySelector('text').getBoundingClientRect();
      const left = (box.left - frame.left) / scale;
      const top = (box.top - frame.top) / scale;
      return [entry.getAttribute('data-legend-key'), [left, top, left + box.width / scale, top + box.height / scale]];
    });
  return { canvas: [width, height], rows: Object.fromEntries(rows) };
});
const expectInsideCanvas = ({ canvas: [width, height], rows }) => {
  const outside = Object.entries(rows).filter(([, [left, top, right, bottom]]) => (
    left < 0 || top < 0 || right > width || bottom > height
  ));
  expect(outside).toEqual([]);
};
const expectSameBoxes = (actual, expected) => {
  const near = (left, right) => left.length === right.length && left.every((value, index) => Math.abs(value - right[index]) <= 1);
  expect(Object.keys(actual.rows).sort()).toEqual(Object.keys(expected.rows).sort());
  expect(near(actual.canvas, expected.canvas), `canvas ${actual.canvas} vs ${expected.canvas}`).toBe(true);
  const moved = Object.keys(expected.rows).filter(caption => !near(actual.rows[caption], expected.rows[caption]))
    .map(caption => `${caption}: ${actual.rows[caption].map(Math.round)} vs ${expected.rows[caption].map(Math.round)}`);
  expect(moved).toEqual([]);
};

// OV-122 (PD-OI-066, OIC-027): Generate draws a row added in the Legend editor,
// and the Legend around it, where the live add placed them, inside the canvas.
for (const mode of ['linear', 'circular']) {
  test(`M3 ${mode}: a Legend editor added row stays where the live add placed it, inside the canvas`, async ({ browser }) => {
    test.setTimeout(600_000);
    const page = await loadGenerated(browser, mode);
    try {
      await addLegendRow(page, 'Manual row', '#118833');
      const live = await legendCaptionBoxes(page);
      expect(Object.keys(live.rows)).toContain('Manual row');
      expectInsideCanvas(live);
      await generate(page);
      const generated = await legendCaptionBoxes(page);
      expectInsideCanvas(generated);
      expectSameBoxes(generated, live);
      expect(page.externalRequests).toEqual([]);
    } finally {
      await page.context().close();
    }
  });
}

test('L1-L8 generated legend categories reconcile while valid category and layout preferences survive', async ({ browser }, testInfo) => {
  test.setTimeout(1_800_000);
  const page = await load(browser, seeds.lambda);
  page.setDefaultTimeout(180_000);
  let fresh;
  try {
    await generate(page);
    const a = await expectEntries(page, ['CDS']);
    const palette = page.getByLabel('Palette', { exact: true });
    const alternatives = await palette.locator('option').evaluateAll(options => options.map(o => o.value));
    await (await reveal(palette)).selectOption(alternatives.find(value => value !== a.palette));
    await (await reveal(page.locator('#legend-position'))).selectOption('left');
    await page.evaluate(() => {
      const app = window.__GBDRAW_APP__;
      app.adv.legend_font_size = 18;
      app.updateLegendEntryColor(app.legendEntries.findIndex(e => e.caption === 'CDS'), '#123456');
    });
    await generate(page);
    const edited = await inspect(page);
    await upload(page, 'tobacco');
    await generate(page);
    const b = await expectEntries(page, ['CDS', 'repeat_region', 'rRNA', 'tRNA']);
    expect(b.featureIds.some(id => a.featureIds.includes(id))).toBe(false);
    expect(b.result.find(e => e.caption === 'CDS').color).toBe('#123456');
    for (const key of ['side', 'fontSize', 'preferences', 'palette', 'colors']) expect(b[key]).toEqual(edited[key]);
    await switchMode(page, 'circular');
    await switchMode(page, 'linear');
    expect((await inspect(page)).overrides).toEqual(b.overrides);
    await generate(page);
    await upload(page, 'lambda');
    await generate(page);
    const returned = await expectEntries(page, ['CDS']);
    for (const [action, captions] of [['Undo', b.entries], ['Redo', returned.entries]]) {
      await page.getByRole('button', { name: action, exact: true }).click();
      await expect.poll(() => page.evaluate(() => !window.__GBDRAW_HISTORY__.restoring.value
        && !window.__GBDRAW_HISTORY__.capturing.value)).toBe(true);
      await expectEntries(page, captions);
    }
    await generate(page);
    expect((await inspect(page)).result).toEqual(returned.result);
    const saved = testInfo.outputPath('legend-replacement.gbdraw-session.json.gz');
    await download(page, 'Save Session', saved);
    fresh = await load(browser, saved);
    await expectEntries(fresh, ['CDS']);
    await generate(fresh);
    expect((await expectEntries(fresh, ['CDS'])).result).toEqual(returned.result);
    expect(page.externalRequests).toEqual([]);
    expect(fresh.externalRequests).toEqual([]);
  } finally {
    await page.context().close();
    if (fresh) await fresh.context().close();
  }
});

test('G1-G3 rejected candidates preserve A and dormant category intent and manual rows survive A-B-A replay', async ({ browser }, testInfo) => {
  test.setTimeout(1_800_000);
  const page = await load(browser);
  page.setDefaultTimeout(180_000);
  let fresh;
  try {
    await upload(page, 'tobacco', false);
    await generate(page);
    await page.evaluate(async () => {
      const app = window.__GBDRAW_APP__;
      app.newLegendCaption = 'Retained annotation';
      app.newLegendColor = '#884422';
      await app.addNewLegendEntry();
    });
    await expect.poll(() => page.evaluate(() => window.__GBDRAW_APP__.legendEntries.some(e => e.caption === 'Retained annotation'))).toBe(true);
    await page.evaluate(async () => {
      const app = window.__GBDRAW_APP__;
      app.updateLegendEntryColor(app.legendEntries.findIndex(e => e.caption === 'repeat_region'), '#224466');
      app.updateLegendEntryColor(app.legendEntries.findIndex(e => e.caption === 'rRNA'), '#abcdef');
      await app.renameLegendEntry(app.legendEntries.findIndex(e => e.caption === 'rRNA'), 'Ribosomal RNA');
      app.deleteLegendEntry(app.legendEntries.findIndex(e => e.caption === 'tRNA'));
    });
    await generate(page);
    const a = await inspect(page);
    expect(a.entries).toContain('Retained annotation');
    expect(a.original).not.toContain('Retained annotation');
    const input = page.getByLabel('GenBank/DDBJ File', { exact: true });
    await input.setInputFiles({ name: 'corrupt.gb', mimeType: 'text/plain', buffer: Buffer.from('not a GenBank record') });
    expect((await evaluateWithRetainedPromise(page, () => window.__GBDRAW_APP__.runAnalysis())).status).toBe('error');
    const rejected = await inspect(page);
    expect(rejected).toEqual(a);
    await upload(page, 'lambda', false);
    await generate(page);
    const expected = ['CDS', 'GC content', 'GC skew (+)', 'GC skew (-)', 'Retained annotation'];
    const b = await expectEntries(page, expected);
    await generate(page);
    await expectEntries(page, expected);
    const saved = testInfo.outputPath('manual-legend.gbdraw-session.json.gz');
    await download(page, 'Save Session', saved);
    fresh = await load(browser, saved);
    await generate(fresh);
    await expectEntries(fresh, expected);
    await upload(fresh, 'tobacco', false);
    await generate(fresh);
    const returned = await inspect(fresh);
    expect(returned.entries).toContain('Retained annotation');
    expect(returned.entries).toContain('Ribosomal RNA');
    expect(returned.entries).not.toContain('rRNA');
    expect(returned.entries).not.toContain('tRNA');
    expect(returned.result.find(e => e.caption === 'repeat_region').color).toBe('#224466');
    expect(returned.result.find(e => e.caption === 'Ribosomal RNA').color).toBe('#abcdef');
    expect(returned.deleted).toEqual(a.deleted);
    await fs.writeFile(testInfo.outputPath('G1-G3.json'), JSON.stringify({ a, rejected, b, returned }, null, 2));
    expect(page.externalRequests).toEqual([]);
    expect(fresh.externalRequests).toEqual([]);
  } finally {
    await page.context().close();
    if (fresh) await fresh.context().close();
  }
});
