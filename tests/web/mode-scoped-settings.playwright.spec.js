// PD-OI-086 (mode-scoped settings): Circular and Linear each keep their own
// drawing, with its settings and editor edits. A switch shows the other mode's
// drawing and writes nothing into either; an edit, a Generate, a Depth source
// change, Reset Settings, and History act on their own mode's drawing (Reset
// and History on both); a Session keeps both. The cases name the findings
// they guard (OV-80, OV-82, OV-84, OV-101, OV-106, OV-120, OV-142, OV-159)
// and the plan's History cases (H1, H2).
const { test, expect } = require('@playwright/test');
const { readFileSync } = require('node:fs');
const { gunzipSync } = require('node:zlib');
const { evaluateWithRetainedPromise, generateAndWaitForResult } = require('./helpers/app-lifecycle.cjs');
const { openFresh, openWithGenBank } = require('./helpers/audit-browser.cjs');
const { semanticSnapshot, settleLive } = require('./helpers/live-generate-parity.cjs');
const { download } = require('./helpers/mode-transition.cjs');

test.describe.configure({ retries: 0 });

const SINGLE = 'tests/fixtures/forced_label_underlay.gb';
const COLOR = '#7b2cbf';
const MODE_NAMES = { circular: 'Circular', linear: 'Linear' };
// A Depth TSV for the FORCEDLBL record.
const depthTsv = (base = 10) => Array.from({ length: 4 }, (_, index) => `FORCEDLBL\t${index * 700 + 1}\t${base + index}`).join('\n');

const switchMode = async (page, mode) => {
  await page.getByRole('button', { name: MODE_NAMES[mode], exact: true }).click();
  await page.waitForFunction((value) => window.__GBDRAW_APP__?.mode === value, mode);
  await expect.poll(() => page.evaluate(() => window.__GBDRAW_APP__.circularRecordDiscovery?.status)).not.toBe('loading');
  await settleLive(page);
};

const generate = async (page, options) => {
  const outcome = await generateAndWaitForResult(page, options);
  await settleLive(page);
  return outcome;
};

// A Generate whose status is the assertion: the status and the alert it left.
const generateOutcome = async (page) => {
  const outcome = await evaluateWithRetainedPromise(page, async () => {
    const app = window.__GBDRAW_APP__;
    const result = await app.runAnalysis();
    const error = app.errorLog;
    return { status: result?.status ?? null, error: error ? `${error.code}: ${error.summary}` : null };
  });
  await settleLive(page);
  return outcome;
};

const setLinearFile = async (page, path = SINGLE) => {
  await page.evaluate(async ({ text, name }) => {
    window.__GBDRAW_APP__.setLinearSeqPrimaryFile(0, 'gb', new File([text], name, { type: 'text/plain', lastModified: 1000 }));
    await window.Vue.nextTick();
  }, { text: readFileSync(path, 'utf8'), name: path.split('/').pop() });
  await settleLive(page);
};

// The Circular file loaded, Circular shown, labels outside, Auto Reflow off.
const openCircular = (page) => openWithGenBank(page, SINGLE, () => {
  const app = window.__GBDRAW_APP__;
  app.form.labels_mode = 'out';
  app.autoLabelReflowEnabled = false;
});

// The same file in both modes, Circular shown.
const openBoth = async (page) => {
  await openCircular(page);
  await switchMode(page, 'linear');
  await setLinearFile(page);
  await switchMode(page, 'circular');
};

// Both modes generated, Linear shown.
const bothResults = async (page) => {
  await openBoth(page);
  await generate(page);
  await switchMode(page, 'linear');
  await generate(page);
};

const setDepthFile = async (page, mode, index, name, text) => {
  await page.evaluate(({ target, slot, fileName, content }) => {
    const app = window.__GBDRAW_APP__;
    const file = content === null ? null : new File([content], fileName, { type: 'text/tab-separated-values' });
    if (target === 'circular') app.setCircularDepthFile(slot, file);
    else app.setLinearDepthFile(app.linearSeqs[0], slot, file);
  }, { target: mode, slot: index, fileName: name, content: text });
  await settleLive(page);
};

// Values of both drawings: `paths` maps a name to a dotted path in a drawing
// (a ref reads as its value).
const readDrawings = (page, paths) => page.evaluate(async (fieldPaths) => {
  const { state } = await import('/gbdraw/web/js/state.js');
  const unwrap = (value) => (value && typeof value === 'object' && value.__v_isRef ? value.value : value);
  const at = (drawing, path) => path.split('.').reduce((value, key) => unwrap(value)?.[key], drawing);
  const copy = (value) => {
    const plain = unwrap(value);
    return plain === undefined ? null : JSON.parse(JSON.stringify(plain));
  };
  const pick = (drawing) => Object.fromEntries(Object.entries(fieldPaths).map(([name, path]) => [name, copy(at(drawing, path))]));
  return { mode: state.mode.value, circular: pick(state.drawings.circular), linear: pick(state.drawings.linear) };
}, paths);

// One History step that assigns `values` (dotted paths) in the shown mode's drawing.
const editStep = async (page, label, values) => {
  await evaluateWithRetainedPromise(page, async ({ stepLabel, assignments }) => {
    const app = window.__GBDRAW_APP__;
    await window.__GBDRAW_HISTORY__.runUndoable(stepLabel, () => {
      for (const [path, value] of Object.entries(assignments)) {
        const keys = path.split('.');
        const last = keys.pop();
        keys.reduce((target, key) => target[key], app)[last] = value;
      }
    });
  }, { stepLabel: label, assignments: values });
  await settleLive(page);
};

const historyStep = async (page, step) => {
  await evaluateWithRetainedPromise(page, (name) => window.__GBDRAW_HISTORY__[name](), step);
  await settleLive(page);
};
const undoCount = (page) => page.evaluate(() => window.__GBDRAW_HISTORY__.getUndoCount());

// The Legend editor rows and the drawn Legend rows (caption and swatch fill).
const legendState = async (page) => ({
  rows: await page.evaluate(() => window.__GBDRAW_APP__.legendEntries.map((entry) => entry.caption)),
  drawn: (await semanticSnapshot(page)).legend.map(({ caption, fill }) => ({ caption, fill: String(fill).toLowerCase() }))
});
const drawnRow = async (page, caption) => (await legendState(page)).drawn.find((row) => row.caption === caption) || null;

const renameRow = async (page, from, to) => {
  await evaluateWithRetainedPromise(page, async ({ source, target }) => {
    const app = window.__GBDRAW_APP__;
    const index = app.legendEntries.findIndex((entry) => entry.caption === source);
    if (index < 0) throw new Error(`no Legend row "${source}": ${app.legendEntries.map((entry) => entry.caption).join(', ')}`);
    await app.renameLegendEntry(index, target);
  }, { source: from, target: to });
  await settleLive(page);
};

const colorRow = async (page, caption, color = COLOR) => {
  expect(await page.evaluate(({ target, value }) => {
    const app = window.__GBDRAW_APP__;
    const index = app.legendEntries.findIndex((entry) => entry.caption === target);
    return index >= 0 && app.updateLegendEntryColor(index, value);
  }, { target: caption, value: color }), `Legend row ${caption}`).toBeTruthy();
  await settleLive(page);
};

// Save Session, then a fresh page that loads the saved file.
const saveAndLoadFresh = async (page, browser, path) => {
  await download(page, 'Save Session', path);
  const context = await browser.newContext({ baseURL: new URL(page.url()).origin });
  const fresh = await context.newPage();
  await openFresh(fresh);
  await fresh.locator('input[accept^=".json,"]').setInputFiles(path);
  await fresh.waitForFunction(() => !window.__GBDRAW_APP__.sessionImportPending, null, { timeout: 180_000 });
  expect(await fresh.evaluate(() => window.__GBDRAW_APP__.errorLog), 'Load').toBeNull();
  await settleLive(fresh);
  return fresh;
};

// OV-82 G5: a switch to a mode without a Depth file and back leaves Circular
// Show Depth on; each switch is one History step.
test('OV-82 G5: a round trip through a mode without Depth keeps Circular Show Depth on, one History step per switch', async ({ page }) => {
  test.setTimeout(180_000);
  await openCircular(page);
  await setDepthFile(page, 'circular', 0, 'depth.tsv', depthTsv());
  expect((await readDrawings(page, { showDepth: 'form.show_depth' })).circular.showDepth, 'Circular Show Depth after its file').toBe(true);

  const before = await undoCount(page);
  await switchMode(page, 'linear');
  expect(await undoCount(page), 'the switch to Linear is one History step').toBe(before + 1);
  await switchMode(page, 'circular');
  expect(await undoCount(page), 'the switch back is one History step').toBe(before + 2);
  expect(await readDrawings(page, { showDepth: 'form.show_depth' }), 'after the round trip')
    .toMatchObject({ mode: 'circular', circular: { showDepth: true } });

  await historyStep(page, 'undo');
  expect(await readDrawings(page, { showDepth: 'form.show_depth' }), 'Undo of the last switch')
    .toMatchObject({ mode: 'linear', circular: { showDepth: true } });
});

// OV-82 G2: clearing the Linear Depth file leaves Circular Show Depth and its
// series label.
test('OV-82 G2: clearing the Linear Depth file keeps Circular Show Depth on and its series label', async ({ page }) => {
  test.setTimeout(180_000);
  await openBoth(page);
  await setDepthFile(page, 'circular', 0, 'circ-depth.tsv', depthTsv(10));
  const circular = (await readDrawings(page, { showDepth: 'form.show_depth', label: 'adv.depth_tracks.0.label' })).circular;
  expect(circular, 'Circular after its Depth file').toEqual({ showDepth: true, label: 'circ-depth' });
  await switchMode(page, 'linear');
  await setDepthFile(page, 'linear', 0, 'lin-depth.tsv', depthTsv(20));
  expect((await readDrawings(page, { showDepth: 'form.show_depth' })).linear.showDepth, 'Linear Show Depth after its file').toBe(true);

  await setDepthFile(page, 'linear', 0, 'lin-depth.tsv', null);
  await switchMode(page, 'circular');
  expect((await readDrawings(page, { showDepth: 'form.show_depth', label: 'adv.depth_tracks.0.label' })).circular,
    'Circular after the Linear Depth file is cleared').toEqual(circular);
});

// OV-101 G1: removing the Linear Depth track leaves Circular's series label and
// Show Depth.
test('OV-101 G1: removing the Linear Depth track keeps the Circular Depth series label and Show Depth', async ({ page }) => {
  test.setTimeout(180_000);
  await openBoth(page);
  await setDepthFile(page, 'circular', 0, 'circ-depth.tsv', depthTsv(10));
  const fields = { showDepth: 'form.show_depth', label: 'adv.depth_tracks.0.label' };
  const circular = (await readDrawings(page, fields)).circular;
  expect(circular, 'Circular after its Depth file').toEqual({ showDepth: true, label: 'circ-depth' });
  await switchMode(page, 'linear');
  await setDepthFile(page, 'linear', 0, 'lin-depth.tsv', depthTsv(20));

  await page.evaluate(() => window.__GBDRAW_APP__.removeLinearDepthTrack(0));
  await settleLive(page);
  expect((await readDrawings(page, fields)).circular, 'Circular after the Linear Depth track is removed').toEqual(circular);
  await switchMode(page, 'circular');
  expect((await readDrawings(page, fields)).circular, 'Circular shown again').toEqual(circular);
});

// OV-159: label_rendering belongs to its mode. A Circular Generate with Label
// mode None does not rewrite the Linear value, and a Linear Above Feature
// placement does not rewrite the Circular value on Load.
test('OV-159 label_rendering (a): a Circular Generate with Label mode None keeps Linear External Only', async ({ page }) => {
  test.setTimeout(180_000);
  await openCircular(page);
  await switchMode(page, 'linear');
  await editStep(page, 'Label rendering', { 'form.show_labels_linear': 'all', 'adv.label_rendering': 'external_only' });
  await switchMode(page, 'circular');
  await editStep(page, 'Label mode', { 'form.labels_mode': 'none' });
  await generate(page);
  await switchMode(page, 'linear');
  expect((await readDrawings(page, { rendering: 'adv.label_rendering' })).linear.rendering, 'Linear label_rendering')
    .toBe('external_only');
});

test('OV-159 label_rendering (b): Linear Above Feature does not reset Circular External Only through Save and Load', async ({ page, browser }, testInfo) => {
  test.setTimeout(300_000);
  await openBoth(page);
  await switchMode(page, 'linear');
  await editStep(page, 'Label placement', { 'adv.label_placement': 'above_feature' });
  await switchMode(page, 'circular');
  await editStep(page, 'Label rendering', { 'adv.label_rendering': 'external_only' });
  await generate(page);
  const fields = { rendering: 'adv.label_rendering', placement: 'adv.label_placement' };
  expect((await readDrawings(page, fields)).circular.rendering, 'Circular before Save').toBe('external_only');
  const fresh = await saveAndLoadFresh(page, browser, testInfo.outputPath('ov159.gbdraw-session.json'));
  try {
    const loaded = await readDrawings(fresh, fields);
    expect(loaded.circular.rendering, 'Circular label_rendering after Load').toBe('external_only');
    expect(loaded.linear.placement, 'Linear label_placement after Load').toBe('above_feature');
  } finally {
    await fresh.context().close();
  }
});

// OV-80: a Linear Depth row renamed `Coverage` stays Linear's; Circular (its own
// file, no Depth) generates, and the next Linear Generate draws `Coverage` with
// Show Depth on.
test('OV-80 Depth-row rename: a Linear Coverage rename survives a Circular Generate and returns with Show Depth on', async ({ page }) => {
  test.setTimeout(240_000);
  await openBoth(page);
  await switchMode(page, 'linear');
  await setDepthFile(page, 'linear', 0, 'depth.tsv', depthTsv());
  await generate(page);
  await renameRow(page, 'depth', 'Coverage');
  expect((await legendState(page)).drawn.map(({ caption }) => caption)).toContain('Coverage');

  await switchMode(page, 'circular');
  expect(await generateOutcome(page), 'the Circular Generate').toEqual({ status: 'ok', error: null });
  await switchMode(page, 'linear');
  expect((await legendState(page)).rows, 'the Linear Legend rows shown again').toContain('Coverage');
  expect((await readDrawings(page, { showDepth: 'form.show_depth' })).linear.showDepth, 'Linear Show Depth').toBe(true);
  await generate(page);
  expect((await legendState(page)).drawn.map(({ caption }) => caption), 'the next Linear Generate').toContain('Coverage');
});

// OV-80 F5-F7: a Linear Legend move, delete, or added row stays in Linear. The
// Circular Result, shown again and generated, keeps its own Legend; Linear,
// shown again and generated, keeps the edit.
const LINEAR_LEGEND_EDITS = [
  {
    name: 'a moved row (F5)',
    edit: (page) => page.evaluate(() => {
      const app = window.__GBDRAW_APP__;
      app.moveLegendEntryUp(app.legendEntries.findIndex((entry) => entry.caption === 'repeat_region'));
    }),
    rows: ['repeat_region', 'CDS']
  },
  {
    name: 'a deleted row (F6)',
    edit: (page) => evaluateWithRetainedPromise(page, async () => {
      const app = window.__GBDRAW_APP__;
      await app.deleteLegendEntry(app.legendEntries.findIndex((entry) => entry.caption === 'repeat_region'));
    }),
    rows: ['CDS']
  },
  {
    name: 'an added row (F7)',
    edit: (page) => evaluateWithRetainedPromise(page, async (color) => {
      const app = window.__GBDRAW_APP__;
      app.newLegendCaption = 'Extra';
      app.newLegendColor = color;
      await app.addNewLegendEntry();
    }, COLOR),
    rows: ['CDS', 'repeat_region', 'Extra']
  }
];
for (const { name, edit, rows } of LINEAR_LEGEND_EDITS) {
  test(`OV-80 ${name}: a Linear Legend edit stays in Linear across a Circular Generate`, async ({ page }) => {
    test.setTimeout(240_000);
    await openBoth(page);
    await generate(page);
    const circular = await legendState(page);
    await switchMode(page, 'linear');
    await generate(page);
    await edit(page);
    await expect.poll(() => page.evaluate(() => window.__GBDRAW_APP__.legendEntries.map((entry) => entry.caption)),
      { message: 'the Linear edit' }).toEqual(rows);
    await settleLive(page);
    const linear = await legendState(page);

    await switchMode(page, 'circular');
    expect(await legendState(page), 'Circular shown again').toEqual(circular);
    await generate(page);
    expect(await legendState(page), 'the Circular Generate').toEqual(circular);
    await switchMode(page, 'linear');
    expect(await legendState(page), 'Linear shown again').toEqual(linear);
    await generate(page);
    expect((await legendState(page)).rows, 'the next Linear Generate').toEqual(rows);
  });
}

// OV-120: a Legend rename and its color wait while a Generate hides the row,
// are saved in the Session, and draw again when the row returns.
const OV120_CASES = [
  {
    name: 'GC content renamed GC percent while GC content is off',
    row: 'GC content',
    renamed: 'GC percent',
    prepare: async () => {},
    // The Hide GC Content control.
    hide: (page, hidden) => page.evaluate((value) => window.__GBDRAW_APP__.setCircularGcSuppressed(value), hidden)
  },
  {
    name: 'depth renamed Coverage while Show Depth is off',
    row: 'depth',
    renamed: 'Coverage',
    prepare: (page) => setDepthFile(page, 'circular', 0, 'depth.tsv', depthTsv()),
    hide: (page, hidden) => editStep(page, 'Show Depth', { 'form.show_depth': !hidden })
  }
];
for (const { name, row, renamed, prepare, hide } of OV120_CASES) {
  test(`OV-120: ${name} draws again in its color after Save and Load`, async ({ page, browser }, testInfo) => {
    test.setTimeout(300_000);
    await openCircular(page);
    await prepare(page);
    await generate(page);
    await renameRow(page, row, renamed);
    await colorRow(page, renamed);
    expect(await drawnRow(page, renamed), 'the rename on the Result').toEqual({ caption: renamed, fill: COLOR });
    await hide(page, true);
    await settleLive(page);
    await generate(page);
    expect((await legendState(page)).drawn.map(({ caption }) => caption), 'the Generate without the row').not.toContain(renamed);

    const fresh = await saveAndLoadFresh(page, browser, testInfo.outputPath('ov120.gbdraw-session.json'));
    try {
      await hide(fresh, false);
      await settleLive(fresh);
      await generate(fresh);
      const captions = (await legendState(fresh)).drawn.map(({ caption }) => caption);
      expect(captions, 'the row returns under its rename').not.toContain(row);
      expect(await drawnRow(fresh, renamed), 'the row returns renamed, in its color').toEqual({ caption: renamed, fill: COLOR });
    } finally {
      await fresh.context().close();
    }
  });
}

// OV-142 N1/N5: a rename and a move in Circular stay in Circular; pure round
// trips change neither Legend.
test('OV-142 N1/N5: a Circular rename and move leave the Linear Legend, and round trips reorder nothing', async ({ page }) => {
  test.setTimeout(300_000);
  await bothResults(page);
  const linearBefore = await legendState(page);
  await switchMode(page, 'circular');
  await renameRow(page, 'CDS', 'Coding');
  await page.evaluate(() => {
    const app = window.__GBDRAW_APP__;
    app.moveLegendEntryUp(app.legendEntries.findIndex((entry) => entry.caption === 'GC content'));
  });
  await settleLive(page);
  const circularAfter = await legendState(page);
  expect(circularAfter.rows, 'the Circular edits').toContain('Coding');

  for (const trip of [1, 2]) {
    await switchMode(page, 'linear');
    expect(await legendState(page), `Linear on arrival ${trip}`).toEqual(linearBefore);
    await switchMode(page, 'circular');
    expect(await legendState(page), `Circular on return ${trip}`).toEqual(circularAfter);
  }
});

// OV-142 N4: a Circular Legend color on the row Linear renamed does not
// repaint Linear, shown again or generated; Circular draws its color.
test('OV-142 N4: a Circular color on the row Linear renamed leaves the Linear swatch', async ({ page }) => {
  test.setTimeout(300_000);
  await bothResults(page);
  await renameRow(page, 'CDS', 'Coding');
  const linearRenamed = await legendState(page);
  expect(linearRenamed.rows, 'the Linear rename').toContain('Coding');

  await switchMode(page, 'circular');
  expect((await legendState(page)).rows, 'Circular keeps its CDS row').toContain('CDS');
  await colorRow(page, 'CDS', '#654321');
  expect(await drawnRow(page, 'CDS'), 'the Circular swatch').toEqual({ caption: 'CDS', fill: '#654321' });

  await switchMode(page, 'linear');
  expect(await legendState(page), 'Linear shown again').toEqual(linearRenamed);
  await generate(page);
  expect(await legendState(page), 'the next Linear Generate').toEqual(linearRenamed);
  await switchMode(page, 'circular');
  expect(await drawnRow(page, 'CDS'), 'Circular CDS keeps its caption and color').toEqual({ caption: 'CDS', fill: '#654321' });
});

// OV-142 N3: leaving a mode before its Result is presented, then coming back,
// changes neither Legend and raises no error.
test('OV-142 N3: leaving a mode before its Result is presented changes no Legend', async ({ page }) => {
  test.setTimeout(300_000);
  await bothResults(page);
  await switchMode(page, 'circular');
  const circularBefore = await legendState(page);
  await switchMode(page, 'linear');
  await renameRow(page, 'CDS', 'Coding');
  const linearRenamed = await legendState(page);

  // Two switches in one task, then a Linear Generate left before its preview is ready.
  const switched = await evaluateWithRetainedPromise(page, async () => {
    const app = window.__GBDRAW_APP__;
    const there = app.setDiagramMode('circular');
    const back = app.setDiagramMode('linear');
    return { there: there?.status ?? there, back: back?.status ?? back };
  });
  expect(switched, 'both switches run').toEqual({ there: 'ok', back: 'ok' });
  await settleLive(page);
  const left = await evaluateWithRetainedPromise(page, async () => {
    const app = window.__GBDRAW_APP__;
    const result = await app.runAnalysis();
    let away = app.setDiagramMode('circular');
    for (let tries = 0; tries < 200 && away?.status !== 'ok'; tries += 1) {
      await new Promise((resolve) => requestAnimationFrame(resolve));
      away = app.setDiagramMode('circular');
    }
    return { generate: result?.status ?? null, away: away?.status ?? away };
  });
  expect(left, 'the Generate and the early switch').toEqual({ generate: 'ok', away: 'ok' });
  await settleLive(page);

  expect(await page.evaluate(() => window.__GBDRAW_APP__.errorLog), 'no error').toBeNull();
  expect(await legendState(page), 'Circular shown again').toEqual(circularBefore);
  await generate(page);
  expect(await legendState(page), 'the next Circular Generate').toEqual(circularBefore);
  await switchMode(page, 'linear');
  expect(await legendState(page), 'Linear shown again').toEqual(linearRenamed);
});

// OV-106: an unmanaged config override loaded with one mode's Session stays in
// that mode; the other mode generates and saves. Each case loads a Session 44
// fixture with one leaf the Web does not manage added to the saved request's
// config overrides, as the OV-106 probes did.
const OV106_CASES = [
  { fixture: 'tests/fixtures/sessions/HmmtDNA_basic_circular.v44-schema8.gbdraw-session.json.gz', from: 'circular', to: 'linear', path: 'objects.ticks.tick_width', value: 4 },
  { fixture: 'tests/fixtures/sessions/lambda_basic_linear.v44-schema8.gbdraw-session.json.gz', from: 'linear', to: 'circular', path: 'objects.blast_match.curve_tension', value: 0.3 }
];
const withConfigOverride = (fixture, path, value) => {
  const session = JSON.parse(gunzipSync(readFileSync(fixture)));
  expect(session.version, fixture).toBe(44);
  session.renderRequest.diagramOptions.configOverrides[path] = value;
  return Buffer.from(JSON.stringify(session));
};
for (const { fixture, from, to, path, value } of OV106_CASES) {
  test(`OV-106: a ${MODE_NAMES[from]} unmanaged override (${path}) does not fail the ${MODE_NAMES[to]} Generate and Save`, async ({ page }, testInfo) => {
    test.setTimeout(300_000);
    await openFresh(page);
    await page.locator('input[accept^=".json,"]').setInputFiles({
      name: `ov106-${from}.gbdraw-session.json`, mimeType: 'application/json', buffer: withConfigOverride(fixture, path, value)
    });
    await page.waitForFunction(() => !window.__GBDRAW_APP__.sessionImportPending
      && window.__GBDRAW_APP__.results.length > 0, null, { timeout: 180_000 });
    await settleLive(page);
    const overrides = { overrides: 'unmanagedConfigOverrides' };
    const loaded = await readDrawings(page, overrides);
    expect(loaded.mode, 'the Session opens on its mode').toBe(from);
    expect(Object.keys(loaded[from].overrides), `${from} override`).toContain(path);

    await switchMode(page, to);
    if (to === 'linear') await setLinearFile(page);
    else {
      await page.getByLabel('GenBank/DDBJ File', { exact: true }).setInputFiles(SINGLE);
      await expect.poll(() => page.evaluate(() => window.__GBDRAW_APP__.circularRecordList.length), { timeout: 60_000 }).toBeGreaterThan(0);
      await settleLive(page);
    }
    expect(await generateOutcome(page), `the ${to} Generate`).toEqual({ status: 'ok', error: null });
    await download(page, 'Save Session', testInfo.outputPath('ov106.gbdraw-session.json'));
    const after = await readDrawings(page, overrides);
    expect(Object.keys(after[from].overrides), `the ${from} override stays in ${from}`).toContain(path);
    expect(Object.keys(after[to].overrides), `${to} has no ${from} override`).not.toContain(path);
  });
}

// Depth series count: two Circular Depth series do not reach the Linear
// request, which has one Linear Depth file.
test('Depth series count: two Circular series and one Linear Depth file both generate', async ({ page }) => {
  test.setTimeout(240_000);
  await openCircular(page);
  await setDepthFile(page, 'circular', 0, 'depth-a.tsv', depthTsv(10));
  await setDepthFile(page, 'circular', 1, 'depth-b.tsv', depthTsv(30));
  expect(await generateOutcome(page), 'the Circular Generate with two series').toEqual({ status: 'ok', error: null });
  await switchMode(page, 'linear');
  await setLinearFile(page);
  await setDepthFile(page, 'linear', 0, 'depth-l.tsv', depthTsv(20));
  expect(await generateOutcome(page), 'the Linear Generate with one series').toEqual({ status: 'ok', error: null });
});

// OV-84: a replaced Circular file drops the strokes of the features it no
// longer has; the stroke of a feature it keeps stays and is drawn.
const STROKE = { strokeColor: '#ff0000', strokeWidth: 3 };
test('OV-84 source replacement: the stroke of a feature the new Circular file lacks is removed, the other stays', async ({ page }) => {
  test.setTimeout(240_000);
  await openCircular(page);
  await generate(page);
  const keys = await evaluateWithRetainedPromise(page, async ({ strokeColor, strokeWidth }) => {
    const app = window.__GBDRAW_APP__;
    const found = {};
    for (const tag of ['FL1', 'FL2']) {
      const feature = app.extractedFeatures.find((item) => item.locus_tag === tag);
      await app.openFeatureEditorFromList(feature, null);
      await window.Vue.nextTick();
      if (await app.updateClickedFeatureStroke(strokeColor, strokeWidth) !== true) throw new Error(`stroke not applied on ${tag}`);
      if (app.featureStyleScopeDialog.show) await app.handleFeatureStyleScopeChoice('single');
      app.clickedFeature = null;
      found[tag] = feature.stable_override_key;
    }
    return found;
  }, STROKE);
  await settleLive(page);
  const strokes = () => readDrawings(page, { strokes: 'featureStrokeOverrides' }).then((value) => value.circular.strokes);
  expect(Object.keys(await strokes()), 'both strokes stored').toEqual(expect.arrayContaining([keys.FL1, keys.FL2]));

  // The same file without the FL2 CDS.
  const withoutFl2 = readFileSync(SINGLE, 'utf8')
    .replace(/ {5}CDS {13}complement\(2001\.\.2600\)\n(?: {21}\/.*\n)+/, '');
  expect(withoutFl2).not.toContain('FL2');
  await page.getByLabel('GenBank/DDBJ File', { exact: true }).setInputFiles({
    name: 'forced_label_underlay_without_fl2.gb', mimeType: 'text/plain', buffer: Buffer.from(withoutFl2)
  });
  await expect.poll(() => page.evaluate(() => window.__GBDRAW_APP__.circularRecordList.length), { timeout: 60_000 }).toBeGreaterThan(0);
  await settleLive(page);
  await generate(page);
  expect(await page.evaluate(() => window.__GBDRAW_APP__.extractedFeatures.some((item) => item.locus_tag === 'FL2')), 'FL2 is gone').toBe(false);

  const stored = await strokes();
  expect(Object.keys(stored), 'the FL2 stroke is removed').not.toContain(keys.FL2);
  expect(stored[keys.FL1], 'the FL1 stroke stays').toMatchObject(STROKE);
  expect(await page.evaluate(async ({ key, strokeColor, strokeWidth }) => {
    const { getFeatureElements } = await import('/gbdraw/web/js/services/feature-dom.js');
    const app = window.__GBDRAW_APP__;
    const feature = app.extractedFeatures.find((item) => item.stable_override_key === key);
    const svg = new DOMParser().parseFromString(String(app.results[app.selectedResultIndex]?.content || ''), 'image/svg+xml').documentElement;
    const elements = feature ? getFeatureElements(svg, feature.svg_id) : [];
    return elements.length > 0 && elements.every((element) => element.getAttribute('stroke') === strokeColor
      && element.getAttribute('stroke-width') === String(strokeWidth));
  }, { key: keys.FL1, ...STROKE }), 'the FL1 stroke is drawn').toBe(true);
});

// Plan §8 samples: settings and editor edits of one mode do not reach the other.
const SAMPLE_FIELDS = {
  labelFontSize: 'adv.label_font_size',
  axisStrokeWidth: 'adv.axis_stroke_width',
  scaleInterval: 'adv.scale_interval',
  separateStrands: 'form.separate_strands',
  palette: 'selectedPalette',
  rules: 'manualSpecificRules',
  annotationSets: 'annotationSets',
  depthMin: 'adv.depth_min',
  gcContentMode: 'adv.gc_content_mode'
};
const SAMPLE_RULES = {
  circular: { feat: 'CDS', qual: 'locus_tag', val: '^FL1$', color: '#e63946', cap: 'alpha' },
  linear: { feat: 'CDS', qual: 'locus_tag', val: '^FL2$', color: '#2a9d8f', cap: 'beta' }
};
// The sample edits in the shown mode: settings, a palette, a color rule, and an annotation set.
const editSamples = async (page, { values, palette, rule, annotationSet }) => {
  await editStep(page, 'Sample settings', values);
  await page.evaluate((name) => {
    const app = window.__GBDRAW_APP__;
    app.selectedPalette = name;
    app.updatePalette();
  }, palette);
  await settleLive(page);
  await evaluateWithRetainedPromise(page, async (fields) => {
    const app = window.__GBDRAW_APP__;
    Object.assign(app.newSpecRule, fields);
    await app.addSpecificRule();
  }, rule);
  await settleLive(page);
  await page.evaluate((id) => { window.__GBDRAW_APP__.addAnnotationSet(id); }, annotationSet);
  await settleLive(page);
};
const ruleSummary = (rules) => (rules || []).map((rule) => [rule.feat, rule.qual, rule.val, rule.color, rule.cap]);
const sampleSummary = (drawing) => ({
  ...drawing,
  rules: ruleSummary(drawing.rules),
  annotationSets: (drawing.annotationSets || []).map((set) => set.id)
});

test('per-mode samples: settings, palette, color rule, annotation set, Depth min, and GC mode stay in their mode', async ({ page }) => {
  test.setTimeout(240_000);
  await openCircular(page);
  const read = async () => {
    const value = await readDrawings(page, SAMPLE_FIELDS);
    return { mode: value.mode, circular: sampleSummary(value.circular), linear: sampleSummary(value.linear) };
  };
  const initial = await read();
  const palettes = (await page.evaluate(() => [...window.__GBDRAW_APP__.paletteNames]))
    .filter((name) => name !== initial.circular.palette && name !== initial.linear.palette);
  expect(palettes.length, 'two other palettes').toBeGreaterThanOrEqual(2);

  await editSamples(page, {
    values: {
      'adv.label_font_size': 13, 'adv.axis_stroke_width': 2.5, 'adv.scale_interval': 500,
      'form.separate_strands': !initial.circular.separateStrands, 'adv.depth_min': 1, 'adv.gc_content_mode': 'percent'
    },
    palette: palettes[0],
    rule: SAMPLE_RULES.circular,
    annotationSet: 'circular-sample'
  });
  const circularEdited = (await read()).circular;
  expect(circularEdited, 'the Circular edits').toMatchObject({
    labelFontSize: 13, axisStrokeWidth: 2.5, scaleInterval: 500, separateStrands: !initial.circular.separateStrands,
    palette: palettes[0], depthMin: 1, gcContentMode: 'percent', annotationSets: ['circular-sample']
  });

  await switchMode(page, 'linear');
  expect((await read()).linear, 'Linear reads its own defaults').toEqual(initial.linear);

  await editSamples(page, {
    values: {
      'adv.label_font_size': 9, 'adv.axis_stroke_width': 1.5, 'adv.scale_interval': 200,
      'form.separate_strands': !initial.linear.separateStrands, 'adv.depth_min': 2, 'adv.gc_content_mode': 'deviation'
    },
    palette: palettes[1],
    rule: SAMPLE_RULES.linear,
    annotationSet: 'linear-sample'
  });
  const linearEdited = (await read()).linear;
  expect(linearEdited, 'the Linear edits').toMatchObject({
    labelFontSize: 9, axisStrokeWidth: 1.5, scaleInterval: 200, separateStrands: !initial.linear.separateStrands,
    palette: palettes[1], depthMin: 2, gcContentMode: 'deviation', annotationSets: ['linear-sample']
  });

  await switchMode(page, 'circular');
  expect(await read(), 'each mode keeps its own values').toEqual({ mode: 'circular', circular: circularEdited, linear: linearEdited });
});

// H1: Undo and Redo across a switch restore each step's mode and both drawings.
test('H1: Undo and Redo of an edit, a switch, and an edit in the other mode restore each step', async ({ page }) => {
  test.setTimeout(240_000);
  await openCircular(page);
  const read = async () => {
    const value = await readDrawings(page, { size: 'adv.label_font_size' });
    return { mode: value.mode, circular: value.circular.size, linear: value.linear.size };
  };
  const s0 = await read();
  const before = await undoCount(page);
  await editStep(page, 'Circular label font size', { 'adv.label_font_size': 13 });
  const s1 = { mode: 'circular', circular: 13, linear: s0.linear };
  expect(await read(), 'after the Circular edit').toEqual(s1);
  await switchMode(page, 'linear');
  const s2 = { ...s1, mode: 'linear' };
  expect(await read(), 'after the switch').toEqual(s2);
  await editStep(page, 'Linear label font size', { 'adv.label_font_size': 9 });
  const s3 = { ...s2, linear: 9 };
  expect(await read(), 'after the Linear edit').toEqual(s3);
  expect(await undoCount(page), 'three History steps').toBe(before + 3);

  for (const [step, expected] of [['undo', s2], ['undo', s1], ['undo', s0], ['redo', s1], ['redo', s2], ['redo', s3]]) {
    await historyStep(page, step);
    expect(await read(), `${step} to ${JSON.stringify(expected)}`).toEqual(expected);
  }
});

// H2: Reset Settings resets both drawings, and one Undo restores both modes' edits.
test('H2: one Undo of Reset Settings restores the edits of both modes', async ({ page }) => {
  test.setTimeout(240_000);
  await openCircular(page);
  const fields = { size: 'adv.label_font_size', width: 'adv.axis_stroke_width' };
  const initial = await readDrawings(page, fields);
  await editStep(page, 'Circular label font size', { 'adv.label_font_size': 13 });
  await switchMode(page, 'linear');
  await editStep(page, 'Linear axis stroke width', { 'adv.axis_stroke_width': 2 });
  const edited = {
    mode: 'linear',
    circular: { size: 13, width: initial.circular.width },
    linear: { size: initial.linear.size, width: 2 }
  };

  await page.evaluate(() => { window.__GBDRAW_APP__.resetSettings(); });
  await settleLive(page);
  expect(await readDrawings(page, fields), 'Reset Settings').toEqual({ ...initial, mode: 'linear' });
  await historyStep(page, 'undo');
  expect(await readDrawings(page, fields), 'one Undo of Reset Settings').toEqual(edited);
});
