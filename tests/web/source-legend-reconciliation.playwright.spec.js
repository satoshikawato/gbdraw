const { test, expect } = require('@playwright/test');
const fs = require('node:fs/promises');
const { gunzipSync } = require('node:zlib');
const { load, generate, switchMode, download, loadEditorLegendRows } = require('./helpers/mode-transition.cjs');
const { evaluateWithRetainedPromise } = require('./helpers/app-lifecycle.cjs');
const { openBatch, openWithGenBank } = require('./helpers/audit-browser.cjs');
const { expectLiveEqualsGenerate, settleLive, showResult } = require('./helpers/live-generate-parity.cjs');

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
  const { getVisibleFeatureLegendGroup } = await import('./js/services/legend-svg.js');
  const digest = async value => [...new Uint8Array(await crypto.subtle.digest('SHA-256',
    new TextEncoder().encode(typeof value === 'string' ? value : JSON.stringify(value))))]
    .map(byte => byte.toString(16).padStart(2, '0')).join('');
  const captions = svg => [...(getVisibleFeatureLegendGroup(svg)?.querySelectorAll('g[data-legend-key]') || [])]
    .map(e => ({ caption: e.getAttribute('data-legend-key'), color: e.querySelector('path[fill]')?.getAttribute('fill') }))
    .sort((a, b) => a.caption.localeCompare(b.caption));
  const result = s.results.value[s.selectedResultIndex.value].content;
  const svg = new DOMParser().parseFromString(result, 'image/svg+xml').documentElement;
  return JSON.parse(JSON.stringify({
    entries: s.activeDrawing().legendEntries.value.map(e => e.caption).sort(), original: s.originalLegendOrder.value,
    entryState: s.activeDrawing().legendEntries.value, originalColors: s.originalLegendColors.value,
    manualRules: s.activeDrawing().manualSpecificRules, strokes: s.activeDrawing().legendStrokeOverrides,
    selectedResult: s.selectedResultIndex.value,
    resultsDigest: await digest(s.results.value),
    mountedDigest: await digest(s.svgContainer.value.querySelector('svg').outerHTML),
    requestDigest: await digest(getCommittedCanonicalRenderRequest()),
    result: captions(svg), mounted: captions(s.svgContainer.value.querySelector('svg')),
    side: s.activeDrawing().form.legend, fontSize: s.activeDrawing().adv.legend_font_size, preferences: s.activeDrawing().layoutPreferences.legend,
    layoutPreferences: s.activeDrawing().layoutPreferences,
    palette: s.activeDrawing().selectedPalette.value, colors: s.activeDrawing().currentColors.value,
    featureIds: s.extractedFeatures.value.map(f => f.stable_feature_id),
    overrides: s.activeDrawing().legendColorOverrides, deleted: s.activeDrawing().deletedLegendEntries.value
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
  const { getAllFeatureLegendGroups, getLegendEntrySwatch } = await import('./js/services/legend-svg.js');
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
      await addLegendRow(page, 'Manual row', '#118833');
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
// A Legend row the editor added, from a Session that holds it (R15-2 retired
// Add legend item): the page's Session saved with the row, loaded, generated.
const addLegendRow = async (page, caption, color) => {
  await loadEditorLegendRows(page, [[caption, color]]);
  await expect.poll(() => page.evaluate(legendIndex, caption)).toBeGreaterThanOrEqual(0);
};
// The swatch stroke of row `caption` (color normalized, width as a number) in each
// Legend group of the selected Result and of the mounted SVG.
const rowStroke = (page, caption) => page.evaluate(async target => {
  const { state: s } = await import('./js/state.js');
  const { getAllFeatureLegendGroups, getLegendEntrySwatch } = await import('./js/services/legend-svg.js');
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

const deleteLegendRow = async (page, caption) => {
  await page.evaluate(row => {
    const app = window.__GBDRAW_APP__;
    app.deleteLegendEntry(app.legendEntries.findIndex(e => e.caption === row));
  }, caption);
  await expect.poll(() => page.evaluate(legendIndex, caption)).toBe(-1);
};
// The Legend rows and the canvas after Generate are those of the live Result.
const expectGenerateKeepsLegend = async (page, step) => {
  const live = await legendCaptionBoxes(page);
  expectInsideCanvas(live);
  await generate(page);
  const generated = await legendCaptionBoxes(page);
  await test.step(step, () => expectSameBoxes(generated, live));
};

// OV-124 (PD-OI-066, OIC-027): Generate draws a Legend from which the Legend
// editor deleted a row as the live delete laid it out: the other rows close the
// gap, and the Legend and the canvas are fitted again.
for (const mode of ['linear', 'circular']) {
  test(`M4 ${mode}: a Legend without a deleted row is laid out live and at Generate alike`, async ({ browser }) => {
    test.setTimeout(600_000);
    const page = await loadGenerated(browser, mode);
    try {
      await deleteLegendRow(page, 'repeat_region');
      await expectGenerateKeepsLegend(page, 'delete a generated row');
      expect(page.externalRequests).toEqual([]);
    } finally {
      await page.context().close();
    }
  });
}

// OV-124 on a batch Result displayed after the delete (D-07, PD-OI-062): the
// display and the display after Generate lay out its Legend alike.
test('M6 circular batch: a Result displayed after a Legend row delete is laid out as after Generate', async ({ browser }) => {
  test.setTimeout(600_000);
  const context = await browser.newContext({ viewport: { width: 1600, height: 1000 } });
  const page = await context.newPage();
  try {
    await openBatch(page);
    await deleteLegendRow(page, 'GC content');
    await showResult(page, 1);
    const displayed = await legendCaptionBoxes(page);
    expect(Object.keys(displayed.rows)).not.toContain('GC content');
    expectInsideCanvas(displayed);
    await generate(page);
    await showResult(page, 1);
    expectSameBoxes(await legendCaptionBoxes(page), displayed);
  } finally {
    await context.close();
  }
});

// OV-124 with OV-122: a generated row and an editor-added row deleted in the
// Legend editor.
for (const mode of ['linear', 'circular']) {
  test(`M5 ${mode}: Legend rows added and deleted in the editor are laid out live and at Generate alike`, async ({ browser }) => {
    test.setTimeout(600_000);
    const page = await loadGenerated(browser, mode);
    try {
      await addLegendRow(page, 'Manual row', '#118833');
      await deleteLegendRow(page, 'repeat_region');
      await expectGenerateKeepsLegend(page, 'an added row, then delete a generated row');
      await deleteLegendRow(page, 'Manual row');
      await expectGenerateKeepsLegend(page, 'delete the added row');
      expect(page.externalRequests).toEqual([]);
    } finally {
      await page.context().close();
    }
  });
}

// The mounted Legend as the editor shows it: the caption boxes and the canvas
// (`legendCaptionBoxes`), the editor rows, and the strokes (color normalized,
// width as a number) of each row's swatch and of each feature.
const legendView = async page => ({
  ...await legendCaptionBoxes(page),
  ...await page.evaluate(async () => {
    const { state: s } = await import('./js/state.js');
    const { getAllFeatureLegendGroups, getLegendEntrySwatch } = await import('./js/services/legend-svg.js');
    const svg = s.svgContainer.value.querySelector('svg');
    const paint = document.createElement('canvas').getContext('2d');
    const stroke = element => {
      const color = element?.getAttribute('stroke') ?? 'none';
      paint.fillStyle = '#010203';
      paint.fillStyle = color;
      return `${color === 'none' ? color : paint.fillStyle} ${Number(element?.getAttribute('stroke-width'))}`;
    };
    const features = {};
    svg.querySelectorAll('[data-gbdraw-feature-id]').forEach(element => {
      const id = element.getAttribute('data-gbdraw-rendered-feature-id') || element.getAttribute('data-gbdraw-feature-id');
      (features[id] ||= new Set()).add(stroke(element));
    });
    return {
      entries: s.activeDrawing().legendEntries.value.map(entry => `${entry.caption} ${entry.color}`),
      swatches: getAllFeatureLegendGroups(svg).map(group => [...group.querySelectorAll('g[data-legend-key]')]
        .map(row => `${row.getAttribute('data-legend-key')}: ${stroke(getLegendEntrySwatch(row))}`).sort()),
      features: Object.fromEntries(Object.entries(features).map(([id, strokes]) => [id, [...strokes].sort().join(', ')]))
    };
  })
});
// How the Legend view `actual` differs from `expected`; places within 1 unit.
const legendViewDifferences = (actual, expected) => {
  const near = (left, right) => Boolean(left && right) && left.length === right.length
    && left.every((value, index) => Math.abs(value - right[index]) <= 1);
  const shown = value => (value ? value.map(Math.round).join(',') : 'none');
  const differences = [];
  if (!near(actual.canvas, expected.canvas)) differences.push(`canvas ${shown(actual.canvas)} vs ${shown(expected.canvas)}`);
  new Set([...Object.keys(actual.rows), ...Object.keys(expected.rows)]).forEach(caption => {
    if (!near(actual.rows[caption], expected.rows[caption])) {
      differences.push(`row ${caption}: ${shown(actual.rows[caption])} vs ${shown(expected.rows[caption])}`);
    }
  });
  ['entries', 'swatches'].forEach(key => {
    if (JSON.stringify(actual[key]) !== JSON.stringify(expected[key])) {
      differences.push(`${key}: ${JSON.stringify(actual[key])} vs ${JSON.stringify(expected[key])}`);
    }
  });
  new Set([...Object.keys(actual.features), ...Object.keys(expected.features)]).forEach(id => {
    if (actual.features[id] !== expected.features[id]) {
      differences.push(`feature ${id}: ${actual.features[id]} vs ${expected.features[id]}`);
    }
  });
  return differences;
};

// OV-125 (R11): each Legend editor action, through the editor's function or its
// control, is one History step. Undo restores the Legend before it (the rows,
// their order and places, the canvas, and the strokes of the swatches and the
// features), Redo the Legend after it, and Generate then draws what live shows.
for (const mode of ['linear', 'circular']) {
  test(`M7 ${mode}: each Legend editor action is one History step that Undo and Redo restore`, async ({ browser }) => {
    test.setTimeout(900_000);
    const page = await loadGenerated(browser, mode);
    const undoCount = () => page.evaluate(() => window.__GBDRAW_HISTORY__.getUndoCount());
    const step = async direction => {
      await evaluateWithRetainedPromise(page, name => window.__GBDRAW_HISTORY__[name](), direction);
      await settleLive(page);
    };
    // `actions` change the Legend one after another, each in one step: Undo
    // then returns each view before an action, and Redo each view after it.
    const expectSteps = (name, actions) => test.step(name, async () => {
      const views = [await legendView(page)];
      const start = await undoCount();
      for (const act of actions) {
        await act();
        await settleLive(page);
        views.push(await legendView(page));
        expect.soft(legendViewDifferences(views.at(-1), views.at(-2)), `${name}: the action changes the Legend`)
          .not.toEqual([]);
      }
      const recorded = await undoCount();
      expect.soft(recorded, `${name}: one Undo step per action`).toBe(start + actions.length);
      if (recorded !== start + actions.length) return;
      for (let index = actions.length - 1; index >= 0; index -= 1) {
        await step('undo');
        expect.soft(legendViewDifferences(await legendView(page), views[index]), `${name}: Undo ${index + 1}`).toEqual([]);
      }
      for (let index = 1; index <= actions.length; index += 1) {
        await step('redo');
        expect.soft(legendViewDifferences(await legendView(page), views[index]), `${name}: Redo ${index}`).toEqual([]);
      }
      expect.soft(await undoCount(), `${name}: Redo returns the steps`).toBe(start + actions.length);
    });
    const strokeWidth = (caption, width) => () => evaluateWithRetainedPromise(page, async row => {
      const app = window.__GBDRAW_APP__;
      await app.updateLegendEntryStrokeWidth(app.legendEntries.findIndex(e => e.caption === row.caption), row.width);
    }, { caption, width });
    // Reset Stroke and the History restore give a row's features the stroke of
    // the renderer's first feature, so each mode strokes a row whose features
    // the renderer draws with that stroke.
    const stroked = mode === 'linear' ? 'CDS' : 'repeat_region';
    try {
      // An editor-added row, from a Session that holds it (R15-2).
      await addLegendRow(page, 'Second row', '#7b2cbf');
      await page.locator('.drawer-toggle').click();
      await page.evaluate(() => window.__GBDRAW_APP__.openRightDrawerTab('legend'));
      await expectSteps('stroke width', [strokeWidth(stroked, 3)]);
      await expectSteps('Reset Stroke', [strokeWidth(stroked, 1), () => evaluateWithRetainedPromise(page, async caption => {
        const app = window.__GBDRAW_APP__;
        await app.resetLegendEntryStroke(app.legendEntries.findIndex(e => e.caption === caption));
      }, stroked)]);
      await expectSteps('delete a row', [() => deleteLegendRow(page, 'repeat_region')]);
      // The second delete uses the editor's Remove control of the row.
      await expectSteps('delete the added row', [
        async () => {
          await page.locator('.right-drawer').getByRole('button', { name: 'Remove Second row', exact: true }).click();
          await expect.poll(() => page.evaluate(legendIndex, 'Second row')).toBe(-1);
        }
      ]);
      // After an Undo, Generate draws the Legend live shows.
      await step('undo');
      const live = await legendView(page);
      expect.soft(Object.keys(live.rows)).toContain('Second row');
      expectInsideCanvas(live);
      await generate(page);
      expect(legendViewDifferences(await legendView(page), live), 'Generate after Undo').toEqual([]);
      expect(page.externalRequests).toEqual([]);
    } finally {
      await page.context().close();
    }
  });
}

// OV-126 (PD-OI-066): an automatic rerender draws the Legend editor's rows as
// Generate does. A Feature Visibility rule that hides FL1 changes a Legend
// source, so it draws the Result again (OV-42).
for (const mode of ['linear', 'circular']) {
  test(`M8 ${mode}: Legend rows added and deleted in the editor keep the live layout through an automatic rerender`, async ({ browser }) => {
    test.setTimeout(600_000);
    const page = await loadGenerated(browser, mode);
    try {
      await addLegendRow(page, 'Manual row', '#118833');
      await deleteLegendRow(page, 'repeat_region');
      await settleLive(page);
      const edited = await legendCaptionBoxes(page);
      await page.evaluate(() => {
        window.__renders = 0;
        window.__GBDRAW_TEST_HOOKS__.beforeDiagramGenerationResponse = () => { window.__renders += 1; };
      });
      await evaluateWithRetainedPromise(page, async () => {
        const app = window.__GBDRAW_APP__;
        await app.addFeatureVisibilityRule();
        const index = app.featureVisibilityManualRules.length - 1;
        const fields = { recordId: 'FORCEDLBL', featureType: 'CDS', qualifier: 'locus_tag', value: '^fl1$', action: 'off' };
        for (const [field, value] of Object.entries(fields)) await app.setFeatureVisibilityRuleField(index, field, value);
      });
      await settleLive(page);
      expect(await page.evaluate(() => window.__renders), 'one automatic rerender').toBe(1);
      const rerendered = await legendCaptionBoxes(page);
      expect(Object.keys(rerendered.rows).sort()).toEqual(Object.keys(edited.rows).sort());
      await expectGenerateKeepsLegend(page, 'Generate after the rerender');
      expect(page.externalRequests).toEqual([]);
    } finally {
      await page.context().close();
    }
  });
}

// OV-127 (PD-OI-066): a rename of a row without features lays the Legend out
// again live; Generate draws the renamed Legend with that layout.
for (const mode of ['linear', 'circular']) {
  test(`M9 ${mode}: a Legend row renamed in the editor is laid out live and at Generate alike`, async ({ browser }) => {
    test.setTimeout(600_000);
    const page = await loadGenerated(browser, mode);
    try {
      if (mode === 'linear') {
        await page.evaluate(() => { window.__GBDRAW_APP__.form.show_gc = true; });
        await generate(page);
      }
      await evaluateWithRetainedPromise(page, async () => {
        const app = window.__GBDRAW_APP__;
        await app.renameLegendEntry(app.legendEntries.findIndex(e => e.caption === 'GC content'),
          'A much longer Legend caption for this row');
      });
      await expect.poll(() => page.evaluate(legendIndex, 'A much longer Legend caption for this row'))
        .toBeGreaterThanOrEqual(0);
      await settleLive(page);
      await expectGenerateKeepsLegend(page, 'rename a row without features');
      expect(page.externalRequests).toEqual([]);
    } finally {
      await page.context().close();
    }
  });
}

// GUI audit FL-07, OV-127: HmmtDNA, Circular, Legend on the left; `GC skew (+)`
// renamed through the Legend editor's text field. The rows keep Python's pitch
// (41.14 units) live and after Generate.
test('M10 circular: a GC skew row renamed in the Legend editor keeps the row pitch live and at Generate', async ({ browser }) => {
  test.setTimeout(600_000);
  const context = await browser.newContext({ viewport: { width: 1600, height: 1000 } });
  const page = await context.newPage();
  const rowTops = () => page.evaluate(() => {
    const svg = window.__GBDRAW_APP__.svgContainer.querySelector('svg');
    return [...svg.querySelectorAll('#legend g[data-legend-key]')].filter(entry => !entry.closest('[display="none"]'))
      .map(entry => entry.querySelector('text').getCTM().f).sort((a, b) => a - b);
  });
  const expectPitch = (tops, step) => expect(tops.slice(1).map((top, index) => Math.round((top - tops[index]) * 10) / 10),
    step).toEqual(tops.slice(1).map(() => 41.1));
  try {
    await openWithGenBank(page, 'tests/test_inputs/HmmtDNA.gbk', () => { window.__GBDRAW_APP__.form.legend = 'left'; });
    await generate(page);
    expectPitch(await rowTops(), 'Generate');
    await page.locator('.drawer-toggle').click();
    await page.evaluate(() => window.__GBDRAW_APP__.openRightDrawerTab('legend'));
    const field = page.locator('.right-drawer input[type="text"]').nth(await page.evaluate(() => (
      [...document.querySelectorAll('.right-drawer input[type="text"]')].findIndex(input => input.value === 'GC skew (+)')
    )));
    await field.fill('skew plus');
    await field.press('Tab');
    await expect.poll(() => page.evaluate(legendIndex, 'skew plus')).toBeGreaterThanOrEqual(0);
    await settleLive(page);
    expectPitch(await rowTops(), 'live rename');
    await expectGenerateKeepsLegend(page, 'Generate after the rename');
    expectPitch(await rowTops(), 'Generate after the rename');
  } finally {
    await context.close();
  }
});

// OV-128: a row added in the Legend editor, renamed, removed, and returned by
// Undo keeps the first row's drawn stroke (OV-121), not that row's stroke edit,
// live and at Generate.
for (const mode of ['linear', 'circular']) {
  test(`M11 ${mode}: a removed added row returned by Undo keeps its stroke live and at Generate`, async ({ browser }) => {
    test.setTimeout(600_000);
    const page = await loadGenerated(browser, mode);
    try {
      const [drawn] = (await rowStroke(page, 'CDS')).result;
      await evaluateWithRetainedPromise(page, async () => {
        const app = window.__GBDRAW_APP__;
        const index = app.legendEntries.findIndex(e => e.caption === 'CDS');
        await app.setLegendEntryStrokeColorValue(index, '#e63946');
        await app.updateLegendEntryStrokeWidth(index, 3);
      });
      await addLegendRow(page, 'Manual row', '#118833');
      await evaluateWithRetainedPromise(page, async () => {
        const app = window.__GBDRAW_APP__;
        await window.__GBDRAW_HISTORY__.runUndoable('Rename legend item', () => (
          app.renameLegendEntry(app.legendEntries.findIndex(e => e.caption === 'Manual row'), 'Manual renamed')));
      });
      await expect.poll(() => page.evaluate(legendIndex, 'Manual renamed')).toBeGreaterThanOrEqual(0);
      await deleteLegendRow(page, 'Manual renamed');
      await evaluateWithRetainedPromise(page, () => window.__GBDRAW_HISTORY__.undo());
      await settleLive(page);
      const { mounted } = await rowStroke(page, 'Manual renamed');
      expect(mounted.length).toBeGreaterThan(0);
      expect(mounted).toEqual(mounted.map(() => drawn));
      await generate(page);
      await expectStroke(page, 'Manual renamed', drawn);
      expect(page.externalRequests).toEqual([]);
    } finally {
      await page.context().close();
    }
  });
}

// The editor captions, and the number of drawn rows of each caption the reader sees.
const legendRowsByCaption = page => page.evaluate(() => {
  const app = window.__GBDRAW_APP__;
  const svg = app.svgContainer.querySelector('svg');
  const drawn = {};
  for (const entry of svg.querySelectorAll('#legend g[data-legend-key]')) {
    if (entry.closest('[display="none"]')) continue;
    const caption = entry.querySelector('text')?.textContent.trim() || '';
    drawn[caption] = (drawn[caption] || 0) + 1;
  }
  return { editor: app.legendEntries.map(e => e.caption), drawn };
});
const oneRowEach = ({ editor }) => Object.fromEntries(editor.map(caption => [caption, 1]));

// OV-151 (GUI audit FL-04, PD-OI-061 revision 2): GC skew (+) and GC skew (-)
// draw no features, so renaming one onto the other offers Suffix and Cancel,
// never Merge, and a forced Merge changes nothing. After Suffix each caption has
// one drawn row, live and at Generate, and the editor order stays. The renamed
// row is the last row; in Linear a renamed row before it moves at Generate
// (OV-156, owned by the Legend layout port).
for (const mode of ['linear', 'circular']) {
  test(`M12 ${mode}: a GC skew row renamed onto the other GC skew row offers no Merge`, async ({ browser }) => {
    test.setTimeout(600_000);
    const page = await loadGenerated(browser, mode);
    try {
      if (mode === 'linear') {
        await page.evaluate(() => {
          const { form } = window.__GBDRAW_APP__;
          form.show_gc = true;
          form.show_skew = true;
        });
        await generate(page);
      }
      const before = await legendRowsByCaption(page);
      expect(before.editor.slice(-2)).toEqual(['GC skew (+)', 'GC skew (-)']);
      expect(before.drawn).toEqual(oneRowEach(before));
      const rename = () => evaluateWithRetainedPromise(page, async () => {
        const app = window.__GBDRAW_APP__;
        await app.renameLegendEntry(app.legendEntries.findIndex(e => e.caption === 'GC skew (-)'), 'GC skew (+)');
      });
      const dialog = page.locator('div.fixed', { has: page.getByRole('heading', { name: 'Legend Name Conflict' }) });

      await rename();
      await expect(dialog).toBeVisible();
      await expect(dialog.getByRole('button', { name: /^Merge into existing/ })).toHaveCount(0);
      await expect(dialog.getByRole('button', { name: 'Keep current color and add a suffix' })).toHaveCount(1);
      await expect(dialog.getByRole('button', { name: 'Cancel' })).toHaveCount(1);
      await evaluateWithRetainedPromise(page, () => window.__GBDRAW_APP__.handleLegendRenameChoice('merge'));
      await expect(dialog).toHaveCount(0);
      await settleLive(page);
      expect(await legendRowsByCaption(page), 'a forced Merge changes nothing').toEqual(before);

      await rename();
      await dialog.getByRole('button', { name: 'Keep current color and add a suffix' }).click();
      await expect.poll(() => page.evaluate(legendIndex, 'GC skew (+) (1)')).toBeGreaterThanOrEqual(0);
      await settleLive(page);
      const renamed = await legendRowsByCaption(page);
      expect(renamed.editor).toEqual([...before.editor.slice(0, -1), 'GC skew (+) (1)']);
      expect(renamed.drawn, 'one drawn row per caption, live').toEqual(oneRowEach(renamed));
      await expectGenerateKeepsLegend(page, 'Generate after the Suffix');
      expect(await legendRowsByCaption(page), 'one drawn row per caption at Generate, editor order kept').toEqual(renamed);
      expect(page.externalRequests).toEqual([]);
    } finally {
      await page.context().close();
    }
  });
}

// A saved Session file as JSON (gzip or plain).
const readSession = bytes => JSON.parse((bytes[0] === 0x1f && bytes[1] === 0x8b ? gunzipSync(bytes) : bytes).toString('utf8'));

// OV-157 (R11): a row's Stroke options button only shows or hides its stroke
// controls. It records no History step, and neither the Legend entries (which
// History and the Session hold) nor the saved Session carry it; a stroke edit
// made in it is still one step and is saved.
test('M13 circular: Stroke options is a disclosure without a History step or a saved field', async ({ browser }, testInfo) => {
  test.setTimeout(600_000);
  const page = await loadGenerated(browser, 'circular');
  try {
    const undoCount = () => page.evaluate(() => window.__GBDRAW_HISTORY__.getUndoCount());
    await page.locator('.drawer-toggle').click();
    await page.evaluate(() => window.__GBDRAW_APP__.openRightDrawerTab('legend'));
    const drawer = page.locator('.right-drawer');
    const caption = await page.evaluate(() => window.__GBDRAW_APP__.legendEntries[0].caption);
    const toggle = drawer.getByRole('button', { name: `Stroke options for ${caption}`, exact: true });
    const start = await undoCount();
    await toggle.click();
    await expect(toggle).toHaveAttribute('aria-expanded', 'true');
    await expect(drawer.getByLabel(`Legend stroke color for ${caption}`, { exact: true })).toBeVisible();
    await settleLive(page);
    expect(await undoCount(), 'opening the stroke options records no step').toBe(start);
    expect(await page.evaluate(() => window.__GBDRAW_APP__.legendEntries.filter(e => Object.hasOwn(e, 'showStroke')).length),
      'the Legend entries do not hold the disclosure').toBe(0);

    await evaluateWithRetainedPromise(page, () => window.__GBDRAW_APP__.updateLegendEntryStrokeWidth(0, 2));
    await settleLive(page);
    expect(await undoCount(), 'a stroke width edit in it is one step').toBe(start + 1);
    await toggle.click();
    await expect(toggle).toHaveAttribute('aria-expanded', 'false');
    await settleLive(page);
    expect(await undoCount(), 'closing it records no step').toBe(start + 1);

    const saved = readSession(await download(page, 'Save Session', testInfo.outputPath('stroke-options.gbdraw-session.json.gz')));
    const { legend } = saved.modes.circular.editorState;
    expect(legend.entries.filter(e => Object.hasOwn(e, 'showStroke')), 'the Session does not save the disclosure').toEqual([]);
    expect(Number(legend.strokeOverrides[caption]?.strokeWidth), 'the stroke edit is saved').toBe(2);
    expect(page.externalRequests).toEqual([]);
  } finally {
    await page.context().close();
  }
});

// The drawn Legend captions in reading order (top to bottom, then left to right).
const drawnLegendOrder = page => page.evaluate(() => {
  const svg = window.__GBDRAW_APP__.svgContainer.querySelector('svg');
  return [...svg.querySelectorAll('#legend g[data-legend-key]')]
    .filter(entry => !entry.closest('[display="none"]'))
    .map(entry => ({ caption: entry.querySelector('text')?.textContent.trim() || '', box: entry.querySelector('text').getBoundingClientRect() }))
    .sort((left, right) => (Math.abs(left.box.y - right.box.y) < 2 ? left.box.x - right.box.x : left.box.y - right.box.y))
    .map(({ caption }) => caption);
});
// The drawn Legend's numbers: the canvas, the Legend's offset, and each shown
// row's caption and swatch translations. Zero shift: a live Legend edit and
// the next Generate give the same numbers exactly.
const legendLayoutNumbers = page => page.evaluate(() => {
  const numbers = value => (String(value || '').match(/[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?/g) || []).map(Number);
  const svg = window.__GBDRAW_APP__.svgContainer.querySelector('svg');
  const legend = svg.getElementById('legend');
  return {
    canvas: [...numbers(svg.getAttribute('width')), ...numbers(svg.getAttribute('height')), ...numbers(svg.getAttribute('viewBox'))],
    legend: numbers(legend.getAttribute('transform')),
    rows: [...legend.querySelectorAll('g[data-legend-key]')]
      .filter(entry => !entry.closest('[display="none"]') && !entry.closest('[data-gbdraw-role="comparison-legend"]'))
      .map(entry => [
        entry.getAttribute('data-legend-key'),
        numbers(entry.querySelector('text')?.getAttribute('transform')),
        [...entry.querySelectorAll('path')].flatMap(path => numbers(path.getAttribute('transform')))
      ])
  };
});
// `loadGenerated` with GC content in Linear too, so each mode has a middle row.
const loadWithGc = async (browser, mode) => {
  const page = await loadGenerated(browser, mode);
  if (mode === 'linear') {
    await page.evaluate(() => { window.__GBDRAW_APP__.form.show_gc = true; });
    await generate(page);
  }
  return page;
};

// OV-154: the Legend editor lists the rows it deleted, each with Restore, and
// Restore all. Each click is one History step; a restored row returns at once
// where Generate draws it, and the Session saves the shorter deleted list.
for (const mode of ['linear', 'circular']) {
  test(`M14 ${mode}: deleted Legend rows return through Restore and Restore all`, async ({ browser }, testInfo) => {
    test.setTimeout(900_000);
    const page = await loadWithGc(browser, mode);
    let fresh;
    try {
      const undoCount = () => page.evaluate(() => window.__GBDRAW_HISTORY__.getUndoCount());
      const step = async direction => {
        await evaluateWithRetainedPromise(page, name => window.__GBDRAW_HISTORY__[name](), direction);
        await settleLive(page);
      };
      await page.locator('.drawer-toggle').click();
      await page.evaluate(() => window.__GBDRAW_APP__.openRightDrawerTab('legend'));
      const drawer = page.locator('.right-drawer');
      const deletedList = drawer.getByRole('list', { name: 'Deleted items' });
      const generated = await legendRowsByCaption(page);
      expect(generated.editor).toEqual(expect.arrayContaining(['CDS', 'repeat_region', 'GC content']));
      await expect(deletedList).toHaveCount(0);
      await deleteLegendRow(page, 'repeat_region');
      await deleteLegendRow(page, 'CDS');
      await settleLive(page);
      await expect(deletedList.getByRole('listitem')).toHaveCount(2);
      const without = caption => generated.editor.filter(entry => entry !== caption);

      let start = await undoCount();
      await deletedList.getByRole('button', { name: 'Restore repeat_region' }).click();
      await expect.poll(() => page.evaluate(legendIndex, 'repeat_region')).toBeGreaterThanOrEqual(0);
      await settleLive(page);
      expect(await undoCount(), 'Restore is one step').toBe(start + 1);
      await expect(deletedList.getByRole('listitem')).toHaveCount(1);
      const one = await legendRowsByCaption(page);
      expect(one.editor, 'the row returns at its place').toEqual(without('CDS'));
      expect(one.drawn).toEqual(oneRowEach(one));
      expect(await drawnLegendOrder(page)).toEqual(without('CDS'));
      await step('undo');
      expect((await legendRowsByCaption(page)).editor).toEqual(without('CDS').filter(entry => entry !== 'repeat_region'));
      await expect(deletedList.getByRole('listitem')).toHaveCount(2);
      await step('redo');
      expect(await legendRowsByCaption(page)).toEqual(one);
      await expect(deletedList.getByRole('listitem')).toHaveCount(1);

      const saved = testInfo.outputPath(`restore-${mode}.gbdraw-session.json.gz`);
      const session = readSession(await download(page, 'Save Session', saved));
      expect(session.modes[mode].editorState.legend.deletedEntries.map(entry => entry.caption), 'the Session saves the shorter list').toEqual(['CDS']);
      const restored = await legendLayoutNumbers(page);
      await expectLiveEqualsGenerate(page, { label: `${mode}: Restore` });
      expect(await legendLayoutNumbers(page), 'Restore: the screen is the Generate').toEqual(restored);
      expect((await legendRowsByCaption(page)).editor).toEqual(without('CDS'));

      start = await undoCount();
      await drawer.getByRole('button', { name: 'Restore all' }).click();
      await expect.poll(() => page.evaluate(legendIndex, 'CDS')).toBeGreaterThanOrEqual(0);
      await settleLive(page);
      expect(await undoCount(), 'Restore all is one step').toBe(start + 1);
      await expect(deletedList).toHaveCount(0);
      const all = await legendRowsByCaption(page);
      expect(all.editor, 'every row returns at its place').toEqual(generated.editor);
      expect(await drawnLegendOrder(page)).toEqual(generated.editor);
      const restoredAll = await legendLayoutNumbers(page);
      await expectLiveEqualsGenerate(page, { label: `${mode}: Restore all` });
      expect(await legendLayoutNumbers(page), 'Restore all: the screen is the Generate').toEqual(restoredAll);
      expect(await legendRowsByCaption(page)).toEqual(all);

      fresh = await load(browser, saved);
      fresh.setDefaultTimeout(180_000);
      await fresh.locator('.drawer-toggle').click();
      await fresh.evaluate(() => window.__GBDRAW_APP__.openRightDrawerTab('legend'));
      const loadedList = fresh.locator('.right-drawer').getByRole('list', { name: 'Deleted items' });
      await expect(loadedList.getByRole('listitem')).toHaveCount(1);
      await loadedList.getByRole('button', { name: 'Restore CDS' }).click();
      // The first Python helper after a Load starts the diagram Worker.
      await expect.poll(() => fresh.evaluate(legendIndex, 'CDS'), { timeout: 180_000 }).toBeGreaterThanOrEqual(0);
      await settleLive(fresh);
      expect((await legendRowsByCaption(fresh)).editor, 'a loaded Session restores its row').toEqual(generated.editor);
      expect(page.externalRequests).toEqual([]);
    } finally {
      await page.context().close();
      if (fresh) await fresh.context().close();
    }
  });
}

// OV-167: after a Session load there is no removed copy of a deleted row, so a
// restored row takes the stroke Generate draws for a row of its kind, also when
// the rows left to copy from are GC rows (drawn without a stroke).
for (const mode of ['linear', 'circular']) {
  test(`M16 ${mode}: a row restored after a Session load with only GC rows left takes Generate's stroke`, async ({ browser }, testInfo) => {
    test.setTimeout(900_000);
    const page = await loadWithGc(browser, mode);
    let fresh;
    try {
      const generated = await rowStroke(page, 'CDS');
      expect(generated.result.length).toBeGreaterThan(0);
      await deleteLegendRow(page, 'CDS');
      await deleteLegendRow(page, 'repeat_region');
      await settleLive(page);
      const left = (await legendRowsByCaption(page)).editor;
      expect(left.length, 'GC rows are left').toBeGreaterThan(0);
      expect(left.every(caption => /^GC /.test(caption)), `only GC rows are left: ${left}`).toBe(true);
      const saved = testInfo.outputPath(`ov167-${mode}.gbdraw-session.json.gz`);
      await download(page, 'Save Session', saved);

      fresh = await load(browser, saved);
      fresh.setDefaultTimeout(180_000);
      await fresh.locator('.drawer-toggle').click();
      await fresh.evaluate(() => window.__GBDRAW_APP__.openRightDrawerTab('legend'));
      await fresh.locator('.right-drawer').getByRole('list', { name: 'Deleted items' })
        .getByRole('button', { name: 'Restore CDS' }).click();
      // The first Python helper after a Load starts the diagram Worker.
      await expect.poll(() => fresh.evaluate(legendIndex, 'CDS'), { timeout: 180_000 }).toBeGreaterThanOrEqual(0);
      await settleLive(fresh);
      expect(await rowStroke(fresh, 'CDS'), 'the restored row live').toEqual(generated);
      await generate(fresh);
      expect(await rowStroke(fresh, 'CDS'), 'the row after Generate').toEqual(generated);
      expect(page.externalRequests).toEqual([]);
    } finally {
      await page.context().close();
      if (fresh) await fresh.context().close();
    }
  });
}

// OV-168: Reset Stroke and Reset all strokes give each part the stroke Generate
// draws. In Circular on forced_label_underlay the first feature path is the
// repeat_region underlay, which Python draws without a stroke; the reset
// stroke of a feature is the block stroke, and an underlay keeps none.
const drawnStrokes = page => page.evaluate(async () => {
  const { state: s } = await import('./js/state.js');
  const { getAllFeatureLegendGroups, getLegendEntrySwatch } = await import('./js/services/legend-svg.js');
  const paint = document.createElement('canvas').getContext('2d');
  const color = value => { paint.fillStyle = '#010203'; paint.fillStyle = String(value); return String(paint.fillStyle); };
  const stroke = element => [color(element.getAttribute('stroke')), Number(element.getAttribute('stroke-width'))];
  const read = svg => ({
    features: Object.fromEntries([...svg.querySelectorAll('path[data-gbdraw-feature-id], path[id^="f"]')]
      .map(path => [path.getAttribute('id'), stroke(path)])),
    swatches: Object.fromEntries(getAllFeatureLegendGroups(svg).flatMap((group, index) => [...group.querySelectorAll('g[data-legend-key]')]
      .map(entry => [`${index}:${entry.getAttribute('data-legend-key')}`, stroke(getLegendEntrySwatch(entry))])))
  });
  const content = s.results.value[s.selectedResultIndex.value].content;
  return {
    result: read(new DOMParser().parseFromString(content, 'image/svg+xml').documentElement),
    mounted: read(s.svgContainer.value.querySelector('svg'))
  };
});
test('M17 circular: Reset Stroke and Reset all strokes give the strokes Generate draws, next to an underlay', async ({ browser }) => {
  test.setTimeout(600_000);
  const page = await loadGenerated(browser, 'circular');
  try {
    const generated = await drawnStrokes(page);
    expect(Object.keys(generated.mounted.features).length).toBeGreaterThan(1);
    const strokeRows = async () => {
      for (const [caption, stroke, width] of [['CDS', '#e63946', 3], ['repeat_region', '#2a9d8f', 4]]) {
        await evaluateWithRetainedPromise(page, async row => {
          const app = window.__GBDRAW_APP__;
          const index = app.legendEntries.findIndex(e => e.caption === row.caption);
          await app.setLegendEntryStrokeColorValue(index, row.stroke);
          await app.updateLegendEntryStrokeWidth(index, row.width);
        }, { caption, stroke, width });
      }
      await settleLive(page);
      expect(await drawnStrokes(page), 'the rows are stroked').not.toEqual(generated);
    };
    const expectGenerated = async label => {
      await settleLive(page);
      expect(await drawnStrokes(page), `${label}: live`).toEqual(generated);
      await generate(page);
      expect(await drawnStrokes(page), `${label}: after Generate`).toEqual(generated);
    };
    await strokeRows();
    for (const caption of ['CDS', 'repeat_region']) {
      await evaluateWithRetainedPromise(page, async row => {
        const app = window.__GBDRAW_APP__;
        await app.resetLegendEntryStroke(app.legendEntries.findIndex(e => e.caption === row));
      }, caption);
    }
    await expectGenerated('Reset Stroke');
    await strokeRows();
    await evaluateWithRetainedPromise(page, async () => { await window.__GBDRAW_APP__.resetAllStrokes(); });
    await expectGenerated('Reset all strokes');
    expect(page.externalRequests).toEqual([]);
  } finally {
    await page.context().close();
  }
});

// OV-158 (Owner decision 2026-10-07): renaming a feature row in the Legend
// editor turns it into a rule row, which the live edit appended to the Legend
// and Generate kept last. The row keeps its place, live and at Generate,
// recorded as an edited order (PD-OI-063), and Sort by default keeps it there.
for (const mode of ['linear', 'circular']) {
  test(`M15 ${mode}: a renamed feature row keeps its place live and at Generate`, async ({ browser }) => {
    test.setTimeout(600_000);
    const page = await loadWithGc(browser, mode);
    try {
      const before = (await legendRowsByCaption(page)).editor;
      const index = before.indexOf('repeat_region');
      expect(index, 'a middle row').toBeGreaterThan(0);
      expect(index).toBeLessThan(before.length - 1);
      await evaluateWithRetainedPromise(page, async () => {
        const app = window.__GBDRAW_APP__;
        await app.renameLegendEntry(app.legendEntries.findIndex(e => e.caption === 'repeat_region'), 'Repeats');
      });
      await expect.poll(() => page.evaluate(legendIndex, 'Repeats')).toBeGreaterThanOrEqual(0);
      await settleLive(page);
      const renamed = before.map(caption => (caption === 'repeat_region' ? 'Repeats' : caption));
      expect((await legendRowsByCaption(page)).editor, 'live').toEqual(renamed);
      expect(await drawnLegendOrder(page), 'live, drawn').toEqual(renamed);
      await expectLiveEqualsGenerate(page, { label: `${mode}: rename` });
      expect((await legendRowsByCaption(page)).editor, 'after Generate').toEqual(renamed);
      expect(await drawnLegendOrder(page), 'after Generate, drawn').toEqual(renamed);

      await page.evaluate(() => window.__GBDRAW_APP__.sortLegendEntriesByDefault());
      await settleLive(page);
      expect((await legendRowsByCaption(page)).editor, 'Sort by default').toEqual(renamed);
      expect(await drawnLegendOrder(page), 'Sort by default, drawn').toEqual(renamed);
      await expectLiveEqualsGenerate(page, { label: `${mode}: Sort by default` });
      expect((await legendRowsByCaption(page)).editor).toEqual(renamed);
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
    await addLegendRow(page, 'Retained annotation', '#884422');
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
