// PD-OI-066 (R1, R3): Legend edits on the live Result equal the next Generate:
// the automatic rerender a Legend source asks for (OV-35, OV-42 to OV-44),
// Legend colors on rows that the draft removes or that only one Result or mode
// draws (OV-63, OV-80), and Legend renames (OV-62). The matrix of edit kinds is
// in live-generate-parity.playwright.spec.js.
const { expect } = require('@playwright/test');
const { execFileSync } = require('node:child_process');
const { readFileSync, writeFileSync } = require('node:fs');
const path = require('node:path');
const { evaluateWithRetainedPromise, generateAndWaitForResult, warmAppTest: test } = require('./helpers/app-lifecycle.cjs');
const {
  diffSemanticSnapshots, expectLiveEqualsGenerate, semanticSnapshot, settleLive, showResult
} = require('./helpers/live-generate-parity.cjs');
const { loadSessionFile, openFresh, openWithGenBank } = require('./helpers/audit-browser.cjs');
const { download, load } = require('./helpers/mode-transition.cjs');
const {
  SINGLE_FIXTURE, open, generate, popupEdit, addVisibilityRule, appAction, addColorRule, history,
  FL1_OFF, legendRowColor, FL1_ALPHA, DEPTH_TSV, colorLegendRow, openCanvas, renameRow,
  legendRowStrokeColor, editorLegendRows, switchPalette, deleteLegendRow
} = require('./helpers/live-generate-parity-steps.cjs');

test.describe.configure({ retries: 0 });

// Counts the diagram renders (Generate and the automatic rerender) from here on.
const countRenders = async (page) => {
  await page.evaluate(() => {
    window.__parityRenders = 0;
    window.__GBDRAW_TEST_HOOKS__ = {
      ...window.__GBDRAW_TEST_HOOKS__,
      beforeDiagramGenerationResponse: () => { window.__parityRenders += 1; }
    };
  });
  return () => page.evaluate(() => window.__parityRenders);
};

// With Auto Reflow off, an edit that changes no Legend source asks for no
// rerender: the other labels may wait for the reflow (Owner, OV-35), and the
// Result still shows what Generate draws apart from that placement.
test('an edit that changes no Legend source does not rerender with Auto Reflow off', async ({ page }) => {
  test.setTimeout(120_000);
  await open(page, { mode: 'circular', results: 'single', reflow: 'off' });
  await addColorRule(page, FL1_ALPHA);
  await generate(page);
  const renders = await countRenders(page);
  await popupEdit(page, 'FL2', { labelText: 'beta protein, edited live' });
  expect(await renders(), 'label text').toBe(0);
  await appAction(page, 'setSpecificRuleField', 0, 'color', '#7b2cbf');
  expect(await renders(), 'color of a rule a drawn feature uses').toBe(0);
  await popupEdit(page, 'FL2', { labelVisibility: 'off' });
  expect(await renders(), 'Label visibility Off').toBe(0);
  await popupEdit(page, 'FL2', { visibility: 'off' });
  expect(await renders(), 'Off for a CDS after the first').toBe(0);
  await expectLiveEqualsGenerate(page, { label: 'edits that change no Legend source' });
});

// An edit that changes a Legend source rerenders once, whether Auto Reflow
// would also place the labels or not.
test('an edit that changes a Legend source rerenders once', async ({ page }) => {
  test.setTimeout(120_000);
  await open(page, { mode: 'circular', results: 'single', reflow: 'off' });
  let renders = await countRenders(page);
  await addVisibilityRule(page, FL1_OFF);
  expect(await renders(), 'visibility rule, Auto Reflow off').toBe(1);
  await page.evaluate(() => { window.__GBDRAW_APP__.autoLabelReflowEnabled = true; });
  renders = await countRenders(page);
  await appAction(page, 'removeFeatureVisibilityRule', 0);
  expect(await renders(), 'visibility rule delete, Auto Reflow on').toBe(1);
  renders = await countRenders(page);
  await addColorRule(page, FL1_ALPHA);
  expect(await renders(), 'color rule, Auto Reflow on').toBe(1);
  await expectLiveEqualsGenerate(page, { label: 'edits that change a Legend source' });
});

// S4 (perf 0.14.x, counts not timings) at the edit port: a color rule that
// adds a Legend row shows it in the rule commit's one show (`showEditorIntent`),
// which lays the Legend out and serializes the Result once. Counted: the Result serializations of the preview owner's commit
// (`flushActiveResult`) until the automatic rerender responds.
const countResultSerializations = async (page) => {
  await page.evaluate(() => {
    window.__resultSerializations = 0;
    window.__countResultSerializations = true;
    if (window.__resultSerializationProbe) return;
    window.__resultSerializationProbe = true;
    const serialize = XMLSerializer.prototype.serializeToString;
    XMLSerializer.prototype.serializeToString = function (node) {
      if (window.__countResultSerializations && /flushActiveResult/.test(new Error().stack || '')) {
        window.__resultSerializations += 1;
      }
      return serialize.call(this, node);
    };
    const hooks = window.__GBDRAW_TEST_HOOKS__ || {};
    const render = hooks.beforeDiagramGenerationResponse;
    window.__GBDRAW_TEST_HOOKS__ = {
      ...hooks,
      beforeDiagramGenerationResponse: (...args) => { window.__countResultSerializations = false; return render?.(...args); }
    };
  });
  return () => page.evaluate(() => window.__resultSerializations);
};

test('a color rule that adds a Legend row serializes the Result once', async ({ page }) => {
  test.setTimeout(120_000);
  await open(page, { mode: 'circular', results: 'single', reflow: 'off' });
  const serializations = await countResultSerializations(page);
  await addColorRule(page, FL1_ALPHA);
  expect(await page.evaluate(() => window.__GBDRAW_APP__.legendEntries.map((entry) => entry.caption))).toContain('alpha');
  expect(await serializations(), 'color rule commit that adds a Legend row').toBe(1);
  await expectLiveEqualsGenerate(page, { label: 'Legend row added by a rule' });
});

// PD-OI-061 revision 4 (was OV-60): renaming a generated Legend row (GC
// content) onto the name of a color rule that draws no row asks. A row without
// features cannot merge into the rule's group, so the dialog offers Keep both
// and Cancel; Generate draws the kept name too.
// PD-OI-061 scenario revision 5: the dialog's button names the row Keep both
// draws, and its sentence the target row's new drawn name; Keep both draws
// exactly those names, live.
const keepBothButton = (page) => page.getByRole('button', { name: /^Keep both: name it "/ });
const chooseKeepBoth = async (page) => {
  const said = await page.evaluate(() => {
    const dialog = window.__GBDRAW_APP__.legendRenameDialog;
    return { name: dialog.keepBothName, after: dialog.targetNameAfter, color: String(dialog.currentColor).toLowerCase() };
  });
  await expect(keepBothButton(page)).toHaveText(`Keep both: name it "${said.name}"`);
  const sentence = page.locator('#legend-rename-description p', { hasText: /The existing "/ });
  await expect(sentence).toHaveCount(said.after ? 1 : 0);
  await evaluateWithRetainedPromise(page, async () => { await window.__GBDRAW_APP__.handleLegendRenameChoice('keep-both'); });
  await settleLive(page);
  const drawn = (await semanticSnapshot(page)).legend;
  expect(drawn.filter((row) => row.caption === said.name).map((row) => row.fill.toLowerCase()), 'the renamed row is drawn with the name the button said')
    .toContain(said.color);
  if (said.after) expect(drawn.map((row) => row.caption), 'the target row is drawn with the name the sentence said').toContain(said.after);
  return said;
};

test('a Legend rename onto the name of a color rule without a drawn row asks and keeps both', async ({ page }) => {
  test.setTimeout(120_000);
  await open(page, { mode: 'circular', results: 'single', reflow: 'on' });
  await addColorRule(page, { feat: 'CDS', qual: 'locus_tag', val: '^NOMATCH$', color: '#c83366', cap: 'Zeta' });
  await renameRow(page, 'GC content', 'Zeta');
  await expect(page.getByRole('heading', { name: 'Legend Name Conflict' })).toBeVisible();
  await expect(page.locator('#legend-rename-description p'), 'OV-586: the dialog says what the name is')
    .toHaveText('"Zeta" is the name of a color rule with another color.');
  await expect(page.getByRole('button', { name: /^Merge into "/ }), 'Merge is not offered').toHaveCount(0);
  await expect(keepBothButton(page), 'no other row draws the name: no number').toHaveText('Keep both: name it "Zeta"');
  await chooseKeepBoth(page);
  const entries = await page.evaluate(() => window.__GBDRAW_APP__.legendEntries.map((entry) => [entry.originalCaption, entry.caption]));
  expect(entries, 'the renamed row keeps the typed name').toContainEqual(['GC content', 'Zeta']);
  await expectLiveEqualsGenerate(page, { label: 'Legend rename onto the name of a color rule without a row, Keep both' });
});

// OV-46: a Legend-only color on a generated row that only one Result of a batch
// draws. The compiler cannot tell which Result draws the row, so each Result may
// miss it; the Generate succeeds when one Result draws it.
test('a Legend color on a row only one Result draws survives the next Generate', async ({ page }) => {
  test.setTimeout(150_000);
  await open(page, { mode: 'circular', results: 'batch', reflow: 'off' });
  await popupEdit(page, 'TESTA_0005', { fill: '#c83366' });
  await generate(page);
  const legendOf = async (index) => {
    await showResult(page, index);
    return (await semanticSnapshot(page)).legend;
  };
  const before = [await legendOf(0), await legendOf(1)];
  expect(before[0].map(({ caption }) => caption)).toContain('other proteins');
  expect(before[1].map(({ caption }) => caption)).not.toContain('other proteins');
  await showResult(page, 0);
  expect(await page.evaluate(() => {
    const app = window.__GBDRAW_APP__;
    return app.updateLegendEntryColor(app.legendEntries.findIndex((entry) => entry.caption === 'other proteins'), '#00aa00');
  })).toBeTruthy();
  await settleLive(page);
  await generate(page);
  const after = [await legendOf(0), await legendOf(1)];
  expect(after[0].find(({ caption }) => caption === 'other proteins').fill.toLowerCase()).toBe('#00aa00');
  expect(after[1]).toEqual(before[1]);
});

// OV-152 (GUI audit FL-02): a popup color for a whole Legend row writes rules
// and keeps a copy of their color as the row's Legend color. Removing those
// rules (Clear All, or deleting the last of them) retires the copy in the same
// History step, so the row takes the palette color its features return to,
// live and at Generate. Undo brings back the rules and the copy. A Legend color
// set on a row that no rule draws stays. "Apply to all" on a palette row writes
// no rule (D-15), so both cases start from a one-feature row.
const COPIED_LEGEND_COLOR_CASES = {
  'Clear All after a one-feature row edit': {
    states: { mode: 'circular', results: 'single', reflow: 'off' },
    color: (page) => popupEdit(page, { type: 'repeat_region' }, { fill: '#e63946' }),
    caption: 'repeat_region',
    direct: 'CDS',
    remove: (page) => appAction(page, 'clearAllSpecificRules')
  },
  'deleting the rule of a one-feature row': {
    states: { mode: 'linear', results: 'single', reflow: 'on' },
    color: (page) => popupEdit(page, { type: 'repeat_region' }, { fill: '#e63946' }),
    caption: 'repeat_region',
    direct: 'CDS',
    remove: (page) => appAction(page, 'removeSpecificRule', 0)
  }
};
for (const [name, { states, color, caption, direct, remove }] of Object.entries(COPIED_LEGEND_COLOR_CASES)) {
  test(`a Legend color copied from popup rules leaves with them: ${name} (${states.mode})`, async ({ page }) => {
    test.setTimeout(180_000);
    await open(page, states);
    const stored = () => page.evaluate(async () => ({ ...(await import('/gbdraw/web/js/state.js')).state.activeDrawing().legendColorOverrides }));
    const ruleCount = () => page.evaluate(() => window.__GBDRAW_APP__.manualSpecificRules.length);
    const liveFill = async () => (await semanticSnapshot(page)).legend.find((row) => row.caption === caption)?.fill;
    await legendRowColor(page, direct, '#7b2cbf');
    await color(page);
    const colored = await stored();
    expect(colored).toEqual({ [direct]: '#7b2cbf', [caption]: '#e63946' });
    const rules = await ruleCount();
    expect(rules).toBeGreaterThan(0);

    await remove(page);
    expect(await ruleCount()).toBe(0);
    expect(await stored(), 'the copy leaves; the Legend color of a row no rule draws stays')
      .toEqual({ [direct]: '#7b2cbf' });
    await history(page, 'undo');
    expect(await ruleCount(), 'Undo brings back the rules').toBe(rules);
    expect(await stored(), 'and the copy, in the same step').toEqual(colored);
    expect(await liveFill()).toBe('#e63946');
    await history(page, 'redo');
    expect(await stored()).toEqual({ [direct]: '#7b2cbf' });

    const palette = (await page.evaluate(async (type) => (
      String((await import('/gbdraw/web/js/state.js')).state.appliedPaletteColors.value[type])
    ), caption)).toLowerCase();
    expect(await liveFill(), 'the live row takes the palette color').toBe(palette);
    const { generated } = await expectLiveEqualsGenerate(page, { label: `${name}: after the removal` });
    expect(generated.legend.find((row) => row.caption === caption)?.fill, 'Generate draws the palette color').toBe(palette);
    expect(generated.legend.find((row) => row.caption === direct)?.fill).toBe('#7b2cbf');
  });
}

// D-15: "Apply to all" on a palette row (no Specific color rule draws it) sets
// the type's default color in one History step. The live Result equals the
// next Generate with one CDS row, and a CDS feature hidden at the time takes
// the color when it is shown again, with no second row.
test('Apply to all on a palette row sets the default color: one row, also for a feature shown later', async ({ page }) => {
  test.setTimeout(240_000);
  await open(page, { mode: 'circular', results: 'single', reflow: 'off' });
  await popupEdit(page, 'FL2', { visibility: 'off' });
  const undoCount = await page.evaluate(() => window.__GBDRAW_HISTORY__.getUndoCount());
  await popupEdit(page, 'FL1', { fill: '#e63946', scope: 'caption' });
  expect(await page.evaluate(() => ({
    color: window.__GBDRAW_APP__.currentColors.CDS,
    rules: window.__GBDRAW_APP__.manualSpecificRules.length,
    undo: window.__GBDRAW_HISTORY__.getUndoCount()
  }))).toEqual({ color: '#e63946', rules: 0, undo: undoCount + 1 });
  const cdsRows = (snapshot) => snapshot.legend.filter((row) => /^(CDS|other proteins)$/.test(row.caption));
  const { generated } = await expectLiveEqualsGenerate(page, { label: 'Apply to all on the CDS palette row' });
  expect(cdsRows(generated).map((row) => [row.caption, row.fill.toLowerCase()])).toEqual([['CDS', '#e63946']]);
  await popupEdit(page, 'FL2', { visibility: 'on' });
  const shown = await expectLiveEqualsGenerate(page, { label: 'the hidden CDS shown again' });
  expect(cdsRows(shown.generated).map((row) => [row.caption, row.fill.toLowerCase()])).toEqual([['CDS', '#e63946']]);
});

// OV-63: Python's Legend row facts excuse only a row the draft removed. A Legend
// style on a key no feature of the records can produce is a stale operation, and
// the Generate still fails at result admission.
test('a Legend style no feature can produce still fails the Generate', async ({ page }) => {
  test.setTimeout(90_000);
  await open(page, { mode: 'circular', results: 'single', reflow: 'off' });
  await page.evaluate(async () => {
    const app = window.__GBDRAW_APP__;
    const { state } = await import('/gbdraw/web/js/state.js');
    state.originalLegendOrder.value = [...state.originalLegendOrder.value, 'Ghost'];
    app.legendEntries.push({ caption: 'Ghost', originalCaption: 'Ghost', color: '#123456', yPos: 400 });
    state.activeDrawing().legendColorOverrides.Ghost = '#123456';
  });
  const outcome = await generateAndWaitForResult(page, { expectedStatus: 'error' });
  expect(outcome.errorSummary).toContain('could not be accepted');
  expect(await page.evaluate(() => JSON.stringify(window.__GBDRAW_APP__.errorLog))).toContain('RESULT_INVALID');
});

// OV-63 siblings: a Legend style on a row that the draft then removes by another
// route than hiding features (a track switched off, a feature type deselected, the
// Legend set to none). Each case colors the row, removes it, and requires that
// Generate succeeds. A track slot edit is live, so its Result must equal
// Generate's; the GC and skew switches, the Features selection, and the Legend
// position apply on Generate, so those cases check the Generate Result for the
// absent rows.
const SIBLINGS = [
  {
    name: 'GC content row after GC content is switched off',
    mode: 'circular',
    absent: ['GC content'],
    run: async (page) => {
      await colorLegendRow(page, 'GC content');
      await appAction(page, 'setCircularGcSuppressed', true);
    }
  },
  {
    name: 'feature type row after the type is deselected in the Features selection',
    mode: 'circular',
    absent: ['repeat_region'],
    run: async (page) => {
      await colorLegendRow(page, 'repeat_region');
      await page.evaluate(() => {
        const features = window.__GBDRAW_APP__.adv.features;
        features.splice(features.indexOf('repeat_region'), 1);
      });
      await settleLive(page);
    }
  },
  {
    name: 'GC skew row after its track is disabled',
    mode: 'circular',
    run: async (page) => {
      await colorLegendRow(page, 'GC skew (+)');
      await page.evaluate(() => {
        const app = window.__GBDRAW_APP__;
        app.setCircularTrackSlotEnabled(app.adv.circular_track_slots.find((slot) => slot.id === 'gc_skew'), false);
      });
      await settleLive(page);
    }
  },
  {
    name: 'GC content row after its track slot is removed',
    mode: 'circular',
    run: async (page) => {
      await colorLegendRow(page, 'GC content');
      await page.evaluate(() => {
        const app = window.__GBDRAW_APP__;
        app.removeCircularTrackSlot(app.adv.circular_track_slots.findIndex((slot) => slot.id === 'gc_content'));
      });
      await settleLive(page);
    }
  },
  {
    name: 'Depth row after Show Depth is switched off, Circular',
    mode: 'circular',
    absent: ['depth'],
    run: async (page) => {
      await page.evaluate((text) => {
        window.__GBDRAW_APP__.setCircularDepthFile(0, new File([text], 'depth.tsv', { type: 'text/tab-separated-values' }));
      }, DEPTH_TSV);
      await settleLive(page);
      await generate(page);
      await colorLegendRow(page, 'depth');
      await page.evaluate(() => { window.__GBDRAW_APP__.setTrackShown('show_depth', false); });
      await settleLive(page);
    }
  },
  {
    name: 'Depth row after Show Depth is switched off, Linear',
    mode: 'linear',
    absent: ['depth'],
    run: async (page) => {
      await page.evaluate((text) => {
        const app = window.__GBDRAW_APP__;
        app.setLinearDepthFile(app.linearSeqs[0], 0, new File([text], 'depth.tsv', { type: 'text/tab-separated-values' }));
      }, DEPTH_TSV);
      await settleLive(page);
      await generate(page);
      await colorLegendRow(page, 'depth');
      await page.evaluate(() => { window.__GBDRAW_APP__.setTrackShown('show_depth', false); });
      await settleLive(page);
    }
  },
  {
    name: 'GC skew row after GC skew is switched off',
    mode: 'circular',
    absent: ['GC skew (-)'],
    run: async (page) => {
      await colorLegendRow(page, 'GC skew (-)');
      await appAction(page, 'setCircularSkewSuppressed', true);
    }
  },
  {
    name: 'a feature row after the Legend is set to none, Circular one record',
    absent: ['CDS'],
    mode: 'circular',
    canvas: false,
    run: async (page) => {
      await colorLegendRow(page, 'CDS');
      await page.evaluate(() => { window.__GBDRAW_APP__.form.legend = 'none'; });
      await settleLive(page);
    }
  },
  {
    name: 'a feature row after the Legend is set to none, Linear',
    absent: ['CDS'],
    mode: 'linear',
    run: async (page) => {
      await colorLegendRow(page, 'CDS');
      await page.evaluate(() => { window.__GBDRAW_APP__.form.legend = 'none'; });
      await settleLive(page);
    }
  },
  {
    name: 'a feature row after the Legend is set to none, Circular Multi-Record Canvas',
    absent: ['CDS'],
    mode: 'circular',
    canvas: true,
    run: async (page) => {
      await colorLegendRow(page, 'CDS');
      await page.evaluate(() => { window.__GBDRAW_APP__.form.legend = 'none'; });
      await settleLive(page);
    }
  }
];

for (const { name, mode, canvas = null, absent = null, run } of SIBLINGS) {
  test(`a Legend color survives Generate: ${name}`, async ({ page }) => {
    test.setTimeout(120_000);
    await openCanvas(page, mode, canvas);
    await run(page);
    if (!absent) {
      await expectLiveEqualsGenerate(page, { label: name });
      return;
    }
    await generate(page);
    const captions = (await semanticSnapshot(page)).legend.map(({ caption }) => caption);
    for (const caption of absent) expect(captions, name).not.toContain(caption);
  });
}

// OV-81: the Legend color of a Depth row hidden by Show Depth is kept; showing
// Depth again draws the row in that color.
const DEPTH_ADDERS = [
  ['circular', (text) => { window.__GBDRAW_APP__.setCircularDepthFile(0, new File([text], 'depth.tsv', { type: 'text/tab-separated-values' })); }],
  ['linear', (text) => { const app = window.__GBDRAW_APP__; app.setLinearDepthFile(app.linearSeqs[0], 0, new File([text], 'depth.tsv', { type: 'text/tab-separated-values' })); }]
];
for (const [mode, addDepth] of DEPTH_ADDERS) {
  test(`a Depth Legend color returns when Show Depth is switched on again (${mode})`, async ({ page }) => {
    test.setTimeout(120_000);
    await openCanvas(page, mode, null);
    await page.evaluate(addDepth, DEPTH_TSV);
    await settleLive(page);
    await generate(page);
    await colorLegendRow(page, 'depth', '#7b2cbf');
    await page.evaluate(() => { window.__GBDRAW_APP__.setTrackShown('show_depth', false); });
    await settleLive(page);
    await generate(page);
    expect((await semanticSnapshot(page)).legend.map(({ caption }) => caption)).not.toContain('depth');
    await page.evaluate(() => { window.__GBDRAW_APP__.setTrackShown('show_depth', true); });
    await settleLive(page);
    await generate(page);
    const row = (await semanticSnapshot(page)).legend.find(({ caption }) => caption === 'depth');
    expect(row?.fill.toLowerCase()).toBe('#7b2cbf');
  });

  // OV-88, OV-120: a Depth row renamed in the Legend is excused like an
  // unrenamed one: its rename and the styles under the new name do not fail the
  // Generate. The Result without the row keeps the rename and its color dormant,
  // so the row returns as `Coverage` in its color.
  test(`a Depth row renamed in the Legend does not fail Generate while Show Depth hides it (${mode})`, async ({ page }) => {
    test.setTimeout(120_000);
    await openCanvas(page, mode, null);
    await page.evaluate(addDepth, DEPTH_TSV);
    await settleLive(page);
    await generate(page);
    await renameRow(page, 'depth', 'Coverage');
    await settleLive(page);
    await colorLegendRow(page, 'Coverage', '#7b2cbf');
    await page.evaluate(() => { window.__GBDRAW_APP__.setTrackShown('show_depth', false); });
    await settleLive(page);
    await generate(page);
    const hidden = (await semanticSnapshot(page)).legend.map(({ caption }) => caption);
    expect(hidden).not.toContain('Coverage');
    expect(hidden).not.toContain('depth');
    await page.evaluate(() => { window.__GBDRAW_APP__.setTrackShown('show_depth', true); });
    await settleLive(page);
    await generate(page);
    const returned = (await semanticSnapshot(page)).legend;
    expect(returned.map(({ caption }) => caption), 'the row returns under its rename').not.toContain('depth');
    expect(returned.find(({ caption }) => caption === 'Coverage')?.fill.toLowerCase(), 'Coverage in its color').toBe('#7b2cbf');
  });
}

// OV-80: Legend styles belong to the drawing of their mode (PD-OI-086). A color
// set on one mode's Result, on a row the other mode does not draw, stays in that
// mode's drawing only, does not fail the other mode's Generate, and returns with
// the row. Show Depth is per mode, so the Depth case returns with Depth on.
const LOCTEST_FIXTURE = 'tests/fixtures/feature_location_search.gb';
const switchMode = async (page, target) => {
  await page.getByRole('button', { name: target === 'linear' ? 'Linear' : 'Circular', exact: true }).click();
  await page.waitForFunction((value) => window.__GBDRAW_APP__?.mode === value, target);
  await settleLive(page);
};
const OTHER_MODE_ROWS = [
  { name: 'a feature type only the Linear file has', from: 'linear', linearFile: LOCTEST_FIXTURE, row: 'tRNA' },
  { name: 'a feature type only the Circular file has', from: 'circular', linearFile: LOCTEST_FIXTURE, row: 'repeat_region' },
  ...['linear', 'circular'].map((from) => ({
    name: `an annotation set of a selected ${from === 'linear' ? 'Linear' : 'Circular'} feature`,
    from,
    row: 'Region X',
    // A selected-feature target belongs to the mode of the Result it was picked on (R2).
    prepare: (page) => evaluateWithRetainedPromise(page, async () => {
      const app = window.__GBDRAW_APP__;
      app.selectedFeatureIds = [app.extractedFeatures.find((feature) => feature.type === 'CDS').svg_id];
      await window.Vue.nextTick();
      app.addAnnotationSet();
      const set = app.annotationSets[app.annotationSets.length - 1];
      await app.addSelectedFeatureAnnotations(set);
      app.setAnnotationSetLegendLabel(set, 'Region X');
    })
  })),
  {
    name: 'a Depth series only Linear has',
    from: 'linear',
    row: 'depth',
    prepare: (page) => page.evaluate((text) => {
      const app = window.__GBDRAW_APP__;
      app.setLinearDepthFile(app.linearSeqs[0], 0, new File([text], 'depth.tsv', { type: 'text/tab-separated-values' }));
    }, DEPTH_TSV)
  }
];
for (const { name, from, linearFile = SINGLE_FIXTURE, row, prepare = null } of OTHER_MODE_ROWS) {
  test(`a Legend color on ${name} survives the other mode's Generate (OV-80)`, async ({ page }) => {
    test.setTimeout(180_000);
    const other = from === 'linear' ? 'circular' : 'linear';
    await open(page, { mode: 'circular', results: 'single', reflow: 'off' });
    await switchMode(page, 'linear');
    await page.evaluate(async ({ text, fileName }) => {
      window.__GBDRAW_APP__.setLinearSeqPrimaryFile(0, 'gb', new File([text], fileName, { type: 'text/plain', lastModified: 1000 }));
      await window.Vue.nextTick();
    }, { text: readFileSync(linearFile, 'utf8'), fileName: linearFile.split('/').pop() });
    await settleLive(page);
    await switchMode(page, from);
    await generate(page);
    if (prepare) {
      await prepare(page);
      await settleLive(page);
      await generate(page);
    }
    await colorLegendRow(page, row);
    await switchMode(page, other);
    await generate(page);
    expect((await semanticSnapshot(page)).legend.map(({ caption }) => caption), `${other} draws no ${row} row`).not.toContain(row);
    expect(await page.evaluate(async ({ caption, own, shown }) => {
      const { state } = await import('/gbdraw/web/js/state.js');
      return {
        stored: state.drawingOf(own).legendColorOverrides[caption] ?? null,
        other: state.drawingOf(shown).legendColorOverrides[caption] ?? null
      };
    }, { caption: row, own: from, shown: other }), `the color stays stored in ${from} only`).toEqual({ stored: '#7b2cbf', other: null });
    await switchMode(page, from);
    await generate(page);
    const drawn = (await semanticSnapshot(page)).legend.find(({ caption }) => caption === row);
    expect(drawn?.fill.toLowerCase(), `${from} draws ${row} in the color`).toBe('#7b2cbf');
  });
}

// OV-62 (PD-OI-061 amended): two Legend rows merge only when both draw features
// of one same type. Anything else offers Keep both and Cancel only.
const legendCaptions = (page) => page.evaluate(() => window.__GBDRAW_APP__.legendEntries.map((entry) => entry.caption));
const choose = async (page, choice) => {
  await evaluateWithRetainedPromise(page, async (picked) => { await window.__GBDRAW_APP__.handleLegendRenameChoice(picked); }, choice);
  await settleLive(page);
};
const mergeButton = (page) => page.getByRole('button', { name: /^Merge into "/ });

test('a Legend rename of GC content onto CDS offers no Merge, and Keep both equals Generate', async ({ page }) => {
  test.setTimeout(120_000);
  await open(page, { mode: 'circular', results: 'single', reflow: 'on' });
  await addColorRule(page, { feat: 'CDS', qual: 'locus_tag', val: '^FL1$', color: '#e63946', cap: 'CDS' });
  const before = await legendCaptions(page);
  await renameRow(page, 'GC content', 'CDS');
  expect(await page.evaluate(() => window.__GBDRAW_APP__.legendRenameDialog.show), 'the rename asks').toBe(true);
  await expect(mergeButton(page), 'Merge is not offered').toHaveCount(0);
  await expect(keepBothButton(page), 'Keep both is offered').toHaveCount(1);
  await choose(page, 'merge');
  expect(await page.evaluate(() => window.__GBDRAW_APP__.legendRenameDialog.show), 'a refused Merge closes the dialog').toBe(false);
  expect(await legendCaptions(page), 'a refused Merge changes no row').toEqual(before);
  await renameRow(page, 'GC content', 'CDS');
  await chooseKeepBoth(page);
  expect(await page.evaluate(() => window.__GBDRAW_APP__.legendEntries.find((entry) => entry.originalCaption === 'GC content')?.caption),
    'the renamed row keeps the typed name').toBe('CDS');
  await expectLiveEqualsGenerate(page, { label: 'GC content renamed onto CDS with Keep both' });
});

// R15-3 (OV-285): a rename onto a deleted row's caption asks as onto a listed
// row. Rows that cannot merge (GC rows have no features) offer Keep both and
// Cancel only; rows of one feature type offer Restore and merge, which is one
// History step.
const legendState = (page) => page.evaluate(() => {
  const app = window.__GBDRAW_APP__;
  return { listed: app.legendEntries.map((entry) => entry.caption), deleted: app.deletedLegendEntries.map((entry) => entry.caption) };
});
const restoreMergeButton = (page) => page.getByRole('button', { name: /^Restore ".*" and merge into it/ });
const conflictAsked = async (page) => {
  expect(await page.evaluate(() => window.__GBDRAW_APP__.legendRenameDialog.show), 'the rename asks').toBe(true);
  await expect(page.getByRole('heading', { name: 'Legend Name Conflict' })).toBeVisible();
};

test('a Legend rename onto a deleted GC skew row asks with Keep both and Cancel only (OV-285)', async ({ page }) => {
  test.setTimeout(120_000);
  await open(page, { mode: 'circular', results: 'single', reflow: 'on' });
  await deleteLegendRow(page, 'GC skew (-)');
  const before = await legendState(page);
  const shown = await semanticSnapshot(page);
  await renameRow(page, 'GC skew (+)', 'GC skew (-)');
  await conflictAsked(page);
  await expect(restoreMergeButton(page), 'Restore and merge is not offered').toHaveCount(0);
  await expect(mergeButton(page), 'Merge is not offered').toHaveCount(0);
  await expect(keepBothButton(page), 'Keep both is offered').toHaveCount(1);
  await choose(page, 'cancel');
  expect(await page.evaluate(() => window.__GBDRAW_APP__.legendRenameDialog.show), 'Cancel closes the dialog').toBe(false);
  expect(await legendState(page), 'Cancel changes no row').toEqual(before);
  expect(diffSemanticSnapshots(shown, await semanticSnapshot(page)), 'Cancel leaves the Result as it was').toEqual([]);
  await renameRow(page, 'GC skew (+)', 'GC skew (-)');
  await conflictAsked(page);
  await chooseKeepBoth(page);
  expect(await legendState(page), 'Keep both keeps the typed name beside the deleted row').toEqual({
    listed: before.listed.map((caption) => (caption === 'GC skew (+)' ? 'GC skew (-)' : caption)),
    deleted: before.deleted
  });
  await expectLiveEqualsGenerate(page, { label: 'GC skew (+) renamed onto the deleted GC skew (-) with Keep both' });
});

// PR 3 (D-25, D-28): two rows of one name keep their own Legend colors under
// their keys; once their swatches match they are one row, live and at
// Generate. The dialog said the names Keep both draws. D-35: an edit of that
// one row edits both.
test('two rows of one name keep their own colors, and one swatch draws one row that edits as one, live = Generate (D-28, D-35)', async ({ page }) => {
  test.setTimeout(180_000);
  await open(page, { mode: 'circular', results: 'single', reflow: 'on' });
  const colorRow = async (key, color) => {
    await page.evaluate(({ row, value }) => {
      const app = window.__GBDRAW_APP__;
      return app.updateLegendEntryColor(app.legendEntries.findIndex((entry) => (entry.originalCaption || entry.caption) === row), value);
    }, { row: key, value: color });
    await settleLive(page);
  };
  const drawn = async () => (await semanticSnapshot(page)).legend
    .filter((row) => row.caption.startsWith('GC skew')).map((row) => [row.caption, row.fill.toLowerCase()]);
  await renameRow(page, 'GC skew (-)', 'GC skew (+)');
  await conflictAsked(page);
  const said = await chooseKeepBoth(page);
  expect([said.name, said.after], 'the button and the sentence').toEqual(['GC skew (+) (2)', 'GC skew (+) (1)']);
  await colorRow('GC skew (+)', '#e63946');
  await colorRow('GC skew (-)', '#2a9d8f');
  expect(await drawn(), 'each row keeps its own color').toEqual([['GC skew (+) (1)', '#e63946'], ['GC skew (+) (2)', '#2a9d8f']]);
  await expectLiveEqualsGenerate(page, { label: 'two rows of one name, each colored' });
  await colorRow('GC skew (-)', '#e63946');
  expect(await drawn(), 'one name and one swatch: one row').toEqual([['GC skew (+)', '#e63946']]);
  // D-35: a color or stroke edit of the one row edits every row it merges:
  // it stays one row, and the Legend tab lists one row, live and after Generate.
  const tabRows = () => page.evaluate(() => {
    const app = window.__GBDRAW_APP__;
    return app.legendEntries.filter((entry) => !app.legendEntryMerged(entry)).map((entry) => app.shownLegendEntryName(entry))
      .filter((name) => name.startsWith('GC skew'));
  });
  const listedKeys = () => page.evaluate(() => window.__GBDRAW_APP__.legendEntries
    .map((entry) => entry.originalCaption || entry.caption).filter((key) => key.startsWith('GC skew')));
  await colorRow('GC skew (+)', '#f4a261');
  expect(await drawn(), 'the recolored row stays one row').toEqual([['GC skew (+)', '#f4a261']]);
  expect(await tabRows(), 'the Legend tab lists one row').toEqual(['GC skew (+)']);
  await expectLiveEqualsGenerate(page, { label: 'two rows of one name and one swatch, recolored as one' });
  expect(await drawn(), 'one row after Generate').toEqual([['GC skew (+)', '#f4a261']]);
  expect(await listedKeys(), 'after Generate the merged row stays listed, not dormant').toEqual(['GC skew (+)', 'GC skew (-)']);
  expect(await tabRows(), 'the Legend tab lists one row after Generate').toEqual(['GC skew (+)']);
  await colorRow('GC skew (+)', '#264653');
  await evaluateWithRetainedPromise(page, async () => {
    const app = window.__GBDRAW_APP__;
    await app.setLegendEntryStrokeColorValue(app.legendEntries.findIndex((entry) => (entry.originalCaption || entry.caption) === 'GC skew (+)'), '#1d3557');
  });
  await settleLive(page);
  expect(await drawn(), 'recolored and stroked after Generate: one row').toEqual([['GC skew (+)', '#264653']]);
  expect(await tabRows(), 'the Legend tab still lists one row').toEqual(['GC skew (+)']);
  await expectLiveEqualsGenerate(page, { label: 'the one row recolored and stroked after Generate' });
  expect(await drawn(), 'one row after the next Generate').toEqual([['GC skew (+)', '#264653']]);
});

const RULE_ROWS = [
  { feat: 'CDS', qual: 'locus_tag', val: '^FL1$', color: '#e63946', cap: 'alpha' },
  { feat: 'CDS', qual: 'locus_tag', val: '^FL2$', color: '#2a9d8f', cap: 'beta' },
  { feat: 'repeat_region', qual: 'note', val: '^RPT_ONE$', color: '#7b2cbf', cap: 'rep' }
];

// Two color-rule rows of CDS: a listed target asks Merge or Keep both (below),
// a deleted one asks Restore and merge; either merge adopts the target's color.
test('a Legend rename onto a deleted row of one feature type restores and merges it in one History step (OV-285)', async ({ page }) => {
  test.setTimeout(180_000);
  await open(page, { mode: 'circular', results: 'single', reflow: 'on' });
  for (const rule of RULE_ROWS.slice(0, 2)) await addColorRule(page, rule);
  await deleteLegendRow(page, 'beta');
  const before = await legendState(page);
  const shown = await semanticSnapshot(page);
  await renameRow(page, 'alpha', 'beta');
  await conflictAsked(page);
  await expect(restoreMergeButton(page), 'Restore and merge is offered').toHaveCount(1);
  await expect(mergeButton(page), 'the listed-row Merge label is not shown').toHaveCount(0);
  await expect(keepBothButton(page), 'Keep both is offered').toHaveCount(1);
  await choose(page, 'merge');
  const merged = await legendState(page);
  expect(merged.listed.filter((caption) => caption === 'beta'), 'the restored row is listed once').toHaveLength(1);
  expect(merged.listed, 'alpha joined beta').not.toContain('alpha');
  expect(merged.deleted, 'beta left the deleted list').not.toContain('beta');
  expect(await page.evaluate(() => {
    const app = window.__GBDRAW_APP__;
    return {
      rules: app.manualSpecificRules.map(({ val, cap, color }) => [val, cap, color]),
      row: app.legendEntries.find((entry) => entry.caption === 'beta').color
    };
  }), 'the row returns in its own color, which the merged feature adopts').toEqual({
    rules: [['^FL1$', 'beta', '#2a9d8f'], ['^FL2$', 'beta', '#2a9d8f']], row: '#2a9d8f'
  });
  await history(page, 'undo');
  expect(await legendState(page), 'one Undo returns the rows before the rename').toEqual(before);
  expect(diffSemanticSnapshots(shown, await semanticSnapshot(page)), 'one Undo shows the Result before the rename').toEqual([]);
  await history(page, 'redo');
  expect(await legendState(page), 'Redo restores and merges again').toEqual(merged);
  await expectLiveEqualsGenerate(page, { label: 'alpha renamed onto the deleted beta row with Restore and merge' });
});

test('a Legend rename between rows of different feature types offers no Merge', async ({ page }) => {
  test.setTimeout(120_000);
  await open(page, { mode: 'circular', results: 'single', reflow: 'on' });
  for (const rule of RULE_ROWS) await addColorRule(page, rule);
  await renameRow(page, 'rep', 'alpha');
  expect(await page.evaluate(() => window.__GBDRAW_APP__.legendRenameDialog.show), 'the rename asks').toBe(true);
  await expect(mergeButton(page), 'Merge is not offered').toHaveCount(0);
  await chooseKeepBoth(page);
  expect((await semanticSnapshot(page)).legend.map((row) => row.caption).filter((caption) => caption.startsWith('alpha')),
    'both rows of the name are numbered in legend order, never drawn by key').toEqual(['alpha (1)', 'alpha (2)']);
  await expectLiveEqualsGenerate(page, { label: 'repeat_region row renamed onto a CDS row with Keep both' });
});

// PD-OI-061 revision 4 (OV-582): two color-rule rows of one feature type. The
// rename asks as for any other row, Merge leaves one row in the target's
// color, and Generate agrees.
// Every Legend key and text of the Result, rule caption, and list option that carries a color ID.
const colorIdNames = (page) => page.evaluate(() => {
  const app = window.__GBDRAW_APP__;
  const rows = [...app.svgContainer.querySelectorAll('g[data-legend-key]')];
  return [
    ...rows.map((row) => row.getAttribute('data-legend-key')),
    ...rows.map((row) => row.querySelector('text')?.textContent || ''),
    ...app.legendEntries.map((entry) => entry.caption),
    ...app.manualSpecificRules.map((rule) => rule.cap),
    ...[...document.querySelectorAll('select[aria-label="Legend names of color rules"] option')].map((option) => option.textContent.trim())
  ].filter((name) => name.includes(' [#'));
});
const ruleTable = (page) => page.evaluate(() => window.__GBDRAW_APP__.manualSpecificRules.map(({ val, cap, color }) => [val, cap, color]));
test('a Legend rename of a rule row onto another rule row of one feature type asks, and Merge equals Generate', async ({ page }) => {
  test.setTimeout(180_000);
  await open(page, { mode: 'circular', results: 'single', reflow: 'on' });
  for (const rule of RULE_ROWS) await addColorRule(page, rule);
  const before = { rows: await legendState(page), rules: await ruleTable(page) };
  const shown = await semanticSnapshot(page);
  await renameRow(page, 'beta', 'alpha');
  await conflictAsked(page);
  await expect(page.locator('#legend-rename-description p').first(), 'a Legend item keeps its sentence')
    .toHaveText('"alpha" already exists with another color.');
  await expect(keepBothButton(page), 'Keep both numbers both rows').toHaveText('Keep both: name it "alpha (2)"');
  await expect(page.locator('#legend-rename-description p').nth(1)).toHaveText('The existing "alpha" becomes "alpha (1)".');
  await page.getByRole('button', { name: 'Merge into "alpha"' }).click();
  await page.waitForFunction(() => !window.__GBDRAW_APP__.legendRenameDialog.show);
  await settleLive(page);
  const merged = await legendState(page);
  expect(merged.listed.filter((caption) => caption === 'alpha'), 'one row keeps the name').toHaveLength(1);
  expect(merged.listed, 'beta joined alpha').not.toContain('beta');
  expect(await ruleTable(page), 'the renamed rule adopts the target\'s name and color').toEqual([
    ['^FL1$', 'alpha', '#e63946'], ['^FL2$', 'alpha', '#e63946'], ['^RPT_ONE$', 'rep', '#7b2cbf']
  ]);
  expect(await colorIdNames(page), 'no name carries a color ID').toEqual([]);
  await history(page, 'undo');
  expect({ rows: await legendState(page), rules: await ruleTable(page) }, 'one Undo returns the rows and the rules').toEqual(before);
  expect(diffSemanticSnapshots(shown, await semanticSnapshot(page)), 'one Undo shows the Result before the rename').toEqual([]);
  await history(page, 'redo');
  expect(await legendState(page), 'Redo merges again').toEqual(merged);
  await expectLiveEqualsGenerate(page, { label: 'same-type rule row merged into a rule row' });
});

// PD-OI-061 revision 4 (OV-580, OV-581): the list under the popup's Legend
// name offers the color-rule groups of the clicked feature's type, and a pick
// joins the picked group: the feature takes its name and its color.
const UNDRAWN_RULE = { feat: 'CDS', qual: 'locus_tag', val: '^NOMATCH$', color: '#5a5a5a', cap: 'Connector' };
const openPopup = async (page, match) => {
  await evaluateWithRetainedPromise(page, async (fields) => {
    const app = window.__GBDRAW_APP__;
    await app.openFeatureEditorFromList(app.filteredFeatures.find((item) => (
      Object.entries(fields).every(([field, value]) => item[field] === value))), null);
    await window.Vue.nextTick();
  }, match);
  return page.getByLabel('Legend names of color rules');
};
const featureFill = async (page, match) => {
  const id = await page.evaluate((fields) => window.__GBDRAW_APP__.filteredFeatures.find((item) => (
    Object.entries(fields).every(([field, value]) => item[field] === value))).svg_id, match);
  return (await semanticSnapshot(page)).features[id]?.fill;
};
const legendFills = async (page) => Object.fromEntries(
  (await semanticSnapshot(page)).legend.map(({ caption, fill }) => [caption, String(fill).toLowerCase()])
);

test('a pick from the popup Legend name list joins the picked group with its color', async ({ page }) => {
  test.setTimeout(180_000);
  await open(page, { mode: 'circular', results: 'single', reflow: 'on' });
  for (const rule of [...RULE_ROWS, UNDRAWN_RULE]) await addColorRule(page, rule);
  const before = { rows: await legendState(page), rules: await ruleTable(page) };
  const shown = await semanticSnapshot(page);
  const list = await openPopup(page, { locus_tag: 'FL2' });
  expect(await list.locator('option').allTextContents().then((names) => names.map((name) => name.trim())),
    'the list offers the CDS groups, each once, drawn or not, and no group of another type')
    .toEqual(['Choose from color rules…', 'alpha', 'beta', 'Connector']);
  await list.selectOption({ label: 'alpha' });
  await page.waitForFunction(() => window.__GBDRAW_APP__.clickedFeature?.appliedLegendName === 'alpha', null, { timeout: 30_000 });
  await settleLive(page);
  expect(await page.evaluate(() => window.__GBDRAW_APP__.legendRenameDialog.show), 'a pick asks no conflict dialog').toBe(false);
  expect(await featureFill(page, { locus_tag: 'FL2' }), 'the feature takes the group\'s color').toEqual(['#e63946']);
  const joined = await legendState(page);
  expect(joined.listed.filter((caption) => caption === 'alpha'), 'one alpha row').toHaveLength(1);
  expect(joined.listed, 'the row the feature left is gone').not.toContain('beta');
  expect((await legendFills(page)).alpha, 'the row keeps its color').toBe('#e63946');
  expect(await colorIdNames(page), 'no name carries a color ID').toEqual([]);
  await page.evaluate(() => { window.__GBDRAW_APP__.clickedFeature = null; });
  await history(page, 'undo');
  expect({ rows: await legendState(page), rules: await ruleTable(page) }, 'one Undo returns the rows and the rules').toEqual(before);
  expect(diffSemanticSnapshots(shown, await semanticSnapshot(page)), 'one Undo shows the Result before the pick').toEqual([]);
  await history(page, 'redo');
  expect(await legendState(page), 'Redo joins again').toEqual(joined);
  await expectLiveEqualsGenerate(page, { label: 'FL2 joined alpha from the popup list' });
  expect(await colorIdNames(page), 'Generate adds no color ID').toEqual([]);
});

test('a pick of a color rule that draws no row asks the scope and draws the feature in the rule color', async ({ page }) => {
  test.setTimeout(180_000);
  await open(page, { mode: 'circular', results: 'single', reflow: 'on' });
  for (const rule of [RULE_ROWS[0], UNDRAWN_RULE]) await addColorRule(page, rule);
  const list = await openPopup(page, { locus_tag: 'FL2' });
  await list.selectOption({ label: 'Connector' });
  await expect(page.getByRole('heading', { name: 'Legend Name Scope' })).toBeVisible();
  await page.getByRole('button', { name: /This feature only/ }).click();
  await page.waitForFunction(() => window.__GBDRAW_APP__.clickedFeature?.appliedLegendName === 'Connector', null, { timeout: 30_000 });
  await settleLive(page);
  expect(await page.evaluate(() => window.__GBDRAW_APP__.legendRenameDialog.show), 'the scope answer is followed by no conflict dialog').toBe(false);
  expect(await featureFill(page, { locus_tag: 'FL2' }), 'the feature takes the rule\'s color').toEqual(['#5a5a5a']);
  expect(await featureFill(page, { product: 'dup alpha' }), 'another feature of the row keeps its color').not.toEqual(['#5a5a5a']);
  expect((await legendFills(page)).Connector, 'the rule\'s row is drawn in its color').toBe('#5a5a5a');
  expect((await ruleTable(page)).filter(([, cap]) => cap.startsWith('Connector')).map(([, cap, color]) => [cap, color]),
    'the rule and the feature share one name and one color').toEqual([['Connector', '#5a5a5a'], ['Connector', '#5a5a5a']]);
  expect(await colorIdNames(page), 'no name carries a color ID').toEqual([]);
  await page.evaluate(() => { window.__GBDRAW_APP__.clickedFeature = null; });
  await expectLiveEqualsGenerate(page, { label: 'FL2 joined the Connector rule from the popup list' });
  expect(await colorIdNames(page), 'Generate adds no color ID').toEqual([]);
});

// The reported step (OV-580): the clicked feature's row is a color rule's
// row, and the answer to the scope question is "All N". The whole row joins
// the picked group: its rule takes the group's name and color.
test('a pick with All N from a color rule\'s row moves the whole row into the picked group', async ({ page }) => {
  test.setTimeout(180_000);
  await open(page, { mode: 'circular', results: 'single', reflow: 'on' });
  await addColorRule(page, { feat: 'CDS', qual: 'product', val: 'protein$', color: '#2a9d8f', cap: 'beta' });
  await addColorRule(page, { feat: 'CDS', qual: 'product', val: '^dup', color: '#e63946', cap: 'alpha' });
  const before = { rows: await legendState(page), rules: await ruleTable(page) };
  const shown = await semanticSnapshot(page);
  expect([await featureFill(page, { locus_tag: 'FL1' }), await featureFill(page, { locus_tag: 'FL2' })],
    'the rule draws both features of the row').toEqual([['#2a9d8f'], ['#2a9d8f']]);
  const list = await openPopup(page, { locus_tag: 'FL2' });
  await list.selectOption({ label: 'alpha' });
  await expect(page.getByRole('heading', { name: 'Legend Name Scope' })).toBeVisible();
  await page.getByRole('button', { name: 'All 2 "beta" features' }).click();
  await page.waitForFunction(() => window.__GBDRAW_APP__.clickedFeature?.appliedLegendName === 'alpha', null, { timeout: 30_000 });
  await settleLive(page);
  expect(await page.evaluate(() => window.__GBDRAW_APP__.legendRenameDialog.show), 'the scope answer is followed by no conflict dialog').toBe(false);
  expect([await featureFill(page, { locus_tag: 'FL1' }), await featureFill(page, { locus_tag: 'FL2' })],
    'both features of the row take the group\'s color').toEqual([['#e63946'], ['#e63946']]);
  const joined = await legendState(page);
  expect(joined.listed.filter((caption) => caption === 'alpha'), 'one alpha row').toHaveLength(1);
  expect(joined.listed, 'the row that joined is gone').not.toContain('beta');
  expect((await legendFills(page)).alpha, 'the row keeps its color').toBe('#e63946');
  expect(await ruleTable(page), 'the row\'s rule takes the group\'s name and color')
    .toEqual([['protein$', 'alpha', '#e63946'], ['^dup', 'alpha', '#e63946']]);
  expect(await colorIdNames(page), 'no name carries a color ID').toEqual([]);
  await page.evaluate(() => { window.__GBDRAW_APP__.clickedFeature = null; });
  await history(page, 'undo');
  expect({ rows: await legendState(page), rules: await ruleTable(page) }, 'one Undo returns the rows and the rules').toEqual(before);
  expect(diffSemanticSnapshots(shown, await semanticSnapshot(page)), 'one Undo shows the Result before the pick').toEqual([]);
  await history(page, 'redo');
  expect(await legendState(page), 'Redo joins again').toEqual(joined);
  await expectLiveEqualsGenerate(page, { label: 'the beta rule row joined alpha from the popup list with All N' });
  expect(await colorIdNames(page), 'Generate adds no color ID').toEqual([]);
});

// D-07: a color rule for every feature type (`*`) is a rule for each type, so
// the list offers its group for any feature, and every group the list offers
// can be joined: also one whose row draws only features of another type.
test('the popup list offers a rule for every feature type, and a pick joins a group with rules of several types', async ({ page }) => {
  test.setTimeout(180_000);
  await open(page, { mode: 'circular', results: 'single', reflow: 'on' });
  for (const rule of [
    RULE_ROWS[2],
    { feat: 'CDS', qual: 'locus_tag', val: '^NOMATCH$', color: RULE_ROWS[2].color, cap: 'rep' },
    { feat: '*', qual: 'locus_tag', val: '^NOMATCH2$', color: '#123456', cap: 'Everything' }
  ]) await addColorRule(page, rule);
  const pick = async (locusTag, name) => {
    const list = await openPopup(page, { locus_tag: locusTag });
    expect(await list.locator('option').allTextContents().then((names) => names.map((text) => text.trim())),
      `${locusTag}: the list offers the every-type group and the group with a CDS rule`)
      .toEqual(['Choose from color rules…', 'Everything', 'rep']);
    await list.selectOption({ label: name });
    await expect(page.getByRole('heading', { name: 'Legend Name Scope' })).toBeVisible();
    await page.getByRole('button', { name: /This feature only/ }).click();
    await page.waitForFunction((caption) => window.__GBDRAW_APP__.clickedFeature?.appliedLegendName === caption, name, { timeout: 30_000 });
    await settleLive(page);
    expect(await page.evaluate(() => window.__GBDRAW_APP__.legendRenameDialog.show), `${name}: a pick is followed by no conflict dialog`).toBe(false);
  };
  await pick('FL1', 'rep');
  expect(await featureFill(page, { locus_tag: 'FL1' }), 'the CDS takes the color of the group its row shares with a repeat').toEqual(['#7b2cbf']);
  await pick('FL2', 'Everything');
  expect(await featureFill(page, { locus_tag: 'FL2' }), 'the CDS takes the color of the every-type rule').toEqual(['#123456']);
  const rows = await legendState(page);
  expect(rows.listed.filter((caption) => ['rep', 'Everything'].includes(caption)).sort(), 'one row per group').toEqual(['Everything', 'rep']);
  expect(await legendFills(page), 'each row in its group\'s color').toMatchObject({ rep: '#7b2cbf', Everything: '#123456' });
  expect(await colorIdNames(page), 'no name carries a color ID').toEqual([]);
  await page.evaluate(() => { window.__GBDRAW_APP__.clickedFeature = null; });
  await expectLiveEqualsGenerate(page, { label: 'FL1 joined a group of two feature types, FL2 the group of an every-type rule' });
  expect(await colorIdNames(page), 'Generate adds no color ID').toEqual([]);
});

// U3a (R14-8): Legend edits the Result executor shows live, without the
// automatic rerender, combined with the edits of other domains. Each case
// ends with live = Generate, in both modes. Auto Reflow decides only whether
// a label or feature visibility edit queues the label reflow
// (`queueLabelReflow` without force); these steps queue none or a forced
// one, and no request carries the setting, so a case runs the same code with
// Auto Reflow on and off. Each runs with it off; the cases with `reflowOn`
// also run with it on in that mode, so a Legend path that starts reading the
// setting is seen.
const ALPHA_ROW = { feat: 'CDS', qual: 'locus_tag', val: '^FL1$', color: '#e63946', cap: 'alpha' };
const renameAndSettle = async (page, caption, name) => { await renameRow(page, caption, name); await settleLive(page); };
const ruleIndex = (page, caption) => page.evaluate(
  (cap) => window.__GBDRAW_APP__.manualSpecificRules.findIndex((rule) => rule.cap === cap), caption
);
const U3A_PARITY = [
  {
    // Linear draws no GC rows here: its featureless row is an editor row,
    // which comes from a Session (R15-2).
    name: '(1) rename a row without features, then color it',
    run: async (page, mode) => {
      const [from, to] = mode === 'circular' ? ['GC content', 'GC %'] : ['Note', 'Note row'];
      if (mode !== 'circular') await editorLegendRows(page, [[from, '#123456']]);
      await renameAndSettle(page, from, to);
      await colorLegendRow(page, to, '#264653');
    }
  },
  {
    name: '(2) rename a palette row with features, then stroke it',
    run: async (page) => { await renameAndSettle(page, 'repeat_region', 'Repeats'); await legendRowStrokeColor(page, 'Repeats', '#e63946'); }
  },
  {
    // Editor rows come from a Session (R15-2).
    name: '(3) color an editor row, then change the palette',
    run: async (page) => {
      await editorLegendRows(page, [['Added', '#123456'], ['Second', '#654321']]);
      await colorLegendRow(page, 'Added', '#abcdef');
      await switchPalette(page, 'arctic');
    }
  },
  {
    // The row returns with the current palette fill and Python's stroke.
    name: '(4) delete a row, change the palette, Undo the delete',
    reflowOn: 'circular',
    run: async (page) => {
      await deleteLegendRow(page, 'repeat_region');
      await switchPalette(page, 'arctic');
      await history(page, 'undo');
    }
  },
  {
    name: '(5) rename a rule row, then color its rule in the Rules panel',
    setup: (page) => addColorRule(page, ALPHA_ROW),
    run: async (page) => {
      await renameAndSettle(page, 'alpha', 'Alpha row');
      await appAction(page, 'setSpecificRuleField', await ruleIndex(page, 'Alpha row'), 'color', '#2a9d8f');
    }
  },
  {
    name: '(6) sort, rename a rule row, delete the renamed row',
    setup: (page) => addColorRule(page, ALPHA_ROW),
    run: async (page, mode) => {
      await page.evaluate(() => window.__GBDRAW_APP__.sortLegendEntries('desc'));
      await settleLive(page);
      await renameAndSettle(page, 'alpha', 'Alpha row');
      if (mode === 'circular') await renameAndSettle(page, 'GC content', 'GC %');
      await deleteLegendRow(page, 'Alpha row');
    }
  },
  {
    // Two rows of one feature type and one color: the rename merges them.
    name: '(7) merge a rule row into another row of its type',
    setup: async (page) => {
      await addColorRule(page, ALPHA_ROW);
      await addColorRule(page, { ...ALPHA_ROW, val: '^FL2$', cap: 'beta' });
    },
    run: async (page) => {
      await renameAndSettle(page, 'beta', 'alpha');
      expect(await page.evaluate(() => window.__GBDRAW_APP__.legendEntries.filter((entry) => entry.caption === 'alpha').length),
        'one alpha row').toBe(1);
    }
  },
  {
    // U3a review H2: Python places the merged row at the first rule of its
    // caption (rule 1 here), before gamma, so the rename asks for the rerender.
    name: '(9) merge a rule row into a later row of its type',
    setup: async (page) => {
      await addColorRule(page, ALPHA_ROW);
      await addColorRule(page, { feat: 'CDS', qual: 'product', val: '^dup alpha$', color: '#2a9d8f', cap: 'gamma' });
      await addColorRule(page, { ...ALPHA_ROW, val: '^FL2$', cap: 'beta' });
    },
    run: async (page) => { await renameAndSettle(page, 'alpha', 'beta'); }
  },
  {
    // U3a review M1: the Rules panel relabels the whole row and asks for the
    // rerender. Redo restores the figure the step kept: the rerender's once it
    // is shown (D-24 (OV-731)), else the commit's, whose new row is appended,
    // so Redo asks for the rerender again.
    name: '(10) Undo and Redo of a Rules-panel caption change of a whole row',
    reflowOn: 'linear',
    setup: (page) => addColorRule(page, ALPHA_ROW),
    run: async (page) => {
      await appAction(page, 'setSpecificRuleField', await ruleIndex(page, 'alpha'), 'cap', 'Alpha row');
      await history(page, 'undo');
      await history(page, 'redo');
    }
  },
  {
    // U3a review L1 (PLAN section 4): Undo and Redo of the case (7) merge
    // restore Results that show the rows exactly, so neither asks Python.
    name: '(11) Undo and Redo of a merge within one type',
    setup: async (page) => {
      await addColorRule(page, ALPHA_ROW);
      await addColorRule(page, { ...ALPHA_ROW, val: '^FL2$', cap: 'beta' });
    },
    run: async (page) => {
      await renameAndSettle(page, 'beta', 'alpha');
      const renders = await countRenders(page);
      await history(page, 'undo');
      await history(page, 'redo');
      expect(await renders(), 'Undo and Redo of the merge ask for no rerender').toBe(0);
    }
  }
];
for (const mode of ['circular', 'linear']) {
  for (const { name, setup = null, run, reflowOn = null } of U3A_PARITY) {
    for (const reflow of reflowOn === mode ? ['off', 'on'] : ['off']) {
      test(`U3a Legend parity ${name} (${mode}, Auto Reflow ${reflow})`, async ({ page }) => {
        test.setTimeout(180_000);
        await open(page, { mode, results: 'single', reflow });
        if (setup) {
          await setup(page);
          await generate(page);
        }
        await run(page, mode);
        await expectLiveEqualsGenerate(page, { label: name });
      });
    }
  }
  if (mode === 'circular') {
    // A batch Result that never drew the delete is not laid out again; the
    // Restore on it and the display of Result 1 show what Generate draws.
    test('U3a Legend parity (8) delete on Result 1, Restore on Result 2 (circular batch, Auto Reflow off)', async ({ page }) => {
      test.setTimeout(180_000);
      await open(page, { mode, results: 'batch', reflow: 'off' });
      await deleteLegendRow(page, 'tRNA');
      await showResult(page, 1);
      await evaluateWithRetainedPromise(page, async () => { await window.__GBDRAW_APP__.restoreAllDeletedLegendEntries(); });
      await settleLive(page);
      await showResult(page, 0);
      await expectLiveEqualsGenerate(page, { label: 'Result 1 after the Restore on Result 2' });
    });
  }
  // U3a review M2: after the deleted type's features are hidden and
  // Generate runs, Python draws no row of it; its row facts say so, and the
  // Restore asks for no rerender.
  test(`U3a Legend parity (12) Restore of a row the Result never drew (${mode}, Auto Reflow off)`, async ({ page }) => {
    test.setTimeout(180_000);
    await open(page, { mode, results: 'single', reflow: 'off' });
    await deleteLegendRow(page, 'repeat_region');
    await addVisibilityRule(page, { recordId: 'FORCEDLBL', featureType: 'repeat_region', qualifier: 'note', value: '^RPT_ONE$', action: 'off' });
    await generate(page);
    const renders = await countRenders(page);
    await evaluateWithRetainedPromise(page, async () => { await window.__GBDRAW_APP__.restoreAllDeletedLegendEntries(); });
    await settleLive(page);
    expect(await renders(), 'the Restore asks for no rerender').toBe(0);
    await expectLiveEqualsGenerate(page, { label: 'Restore of a row the Result never drew' });
  });
}

// U3a compat: a Session 44 that main saved, whose Result has a renamed, a
// deleted, and an added Legend row without Legend structure records (D-15-34)
// (tests/fixtures/sessions/forced-label-underlay-legend-rows.provenance.json).
// Load shows the Legend it saved, Undo of a new edit returns that Result, the
// Restore of the deleted row asks for the automatic rerender (O-2), and live
// equals Generate.
const LEGEND_ROWS_SESSION = 'tests/fixtures/sessions/forced-label-underlay-legend-rows.v44.gbdraw-session.json.gz';
const LEGEND_ROWS_SAVED = JSON.parse(readFileSync('tests/fixtures/sessions/forced-label-underlay-legend-rows.provenance.json', 'utf8'))
  .sessions['forced-label-underlay-legend-rows.v44.gbdraw-session.json.gz'].legend;
const shownLegendRows = async (page) => (await semanticSnapshot(page)).legend.map((row) => row.caption);
const loadLegendRowsSession = async (page, { fresh = true } = {}) => {
  if (fresh) await openFresh(page);
  await loadSessionFile(page, LEGEND_ROWS_SESSION);
  await settleLive(page);
  expect(await legendCaptions(page), 'the saved Legend entries').toEqual(LEGEND_ROWS_SAVED.entries);
  expect(await page.evaluate(() => window.__GBDRAW_APP__.deletedLegendEntries.map((entry) => entry.caption)), 'the deleted rows')
    .toEqual(LEGEND_ROWS_SAVED.deleted);
  expect(await shownLegendRows(page), 'the rows the saved Result shows').toEqual(LEGEND_ROWS_SAVED.shownRows);
};
test('U3a compat: a Session 44 main saved loads its Legend rows, and live equals Generate', async ({ page }) => {
  test.setTimeout(240_000);
  await loadLegendRowsSession(page);
  await expectLiveEqualsGenerate(page, { label: 'Session 44 with Legend row edits, after Load' });
});
// U3a review H1: the saved Result renamed GC content to "GC %" without a
// record of Python's key; Load records it, so each live edit of the row
// reaches it. A Generate after each step would replace the loaded Result, so
// the color, stroke and rename are checked on the live row and then together
// against Generate; the delete starts from a fresh Load.
test('U3a compat: the row a Session 44 main saved renamed takes color, stroke, rename and delete live', async ({ page }) => {
  test.setTimeout(300_000);
  await loadLegendRowsSession(page);
  const liveRow = (caption) => page.evaluate((key) => {
    const row = [...window.__GBDRAW_APP__.svgContainer.querySelectorAll('g[data-legend-key]')]
      .find((element) => element.getAttribute('data-legend-key') === key && element.getAttribute('display') !== 'none');
    const swatch = row?.querySelector('path, rect');
    return row ? { text: row.querySelector('text')?.textContent?.trim(), fill: swatch?.getAttribute('fill'), stroke: swatch?.getAttribute('stroke') } : null;
  }, caption);
  await colorLegendRow(page, 'GC %', '#264653');
  expect((await liveRow('GC %'))?.fill, 'the color reaches the row').toBe('#264653');
  await legendRowStrokeColor(page, 'GC %', '#e63946');
  expect((await liveRow('GC %'))?.stroke, 'the stroke reaches the row').toBe('#e63946');
  await renameAndSettle(page, 'GC %', 'GC ratio');
  expect(await liveRow('GC %'), 'the old caption is gone').toBe(null);
  expect((await liveRow('GC ratio'))?.text, 'the rename reaches the row').toBe('GC ratio');
  await expectLiveEqualsGenerate(page, { label: 'color, stroke and rename of the renamed row after Load' });

  await loadLegendRowsSession(page, { fresh: false });
  await deleteLegendRow(page, 'GC %');
  expect(await liveRow('GC %'), 'the delete hides the row').toBe(null);
  await expectLiveEqualsGenerate(page, { label: 'delete of the renamed row after Load' });
});
test('U3a compat: Undo of a new edit and the Restore of a deleted row on a Session 44 main saved', async ({ page }) => {
  test.setTimeout(240_000);
  await loadLegendRowsSession(page);
  const loaded = await semanticSnapshot(page);
  await colorLegendRow(page, 'CDS', '#264653');
  await history(page, 'undo');
  expect(diffSemanticSnapshots(loaded, await semanticSnapshot(page)), 'Undo returns the loaded Result').toEqual([]);
  const renders = await countRenders(page);
  await evaluateWithRetainedPromise(page, async () => { await window.__GBDRAW_APP__.restoreAllDeletedLegendEntries(); });
  await settleLive(page);
  expect(await renders(), 'the Restore asks for the automatic rerender (O-2)').toBe(1);
  expect(await shownLegendRows(page), 'the restored row is shown').toContain('repeat_region');
  await expectLiveEqualsGenerate(page, { label: 'Restore after Load' });
});
// EU review A-M1: Reset Settings returns the deleted row the saved Result
// lacks as the Restore does (O-2), so the live Legend shows it as Generate
// draws it. Reset also returns labels_mode "out" to the default, a setting
// only Generate draws, so the case compares the Legend.
test('U3a compat: Reset Settings on a Session 44 main saved shows its deleted row as Generate draws it', async ({ page }) => {
  test.setTimeout(240_000);
  await loadLegendRowsSession(page);
  await evaluateWithRetainedPromise(page, async () => { await window.__GBDRAW_APP__.resetSettings(); });
  await settleLive(page);
  const live = await semanticSnapshot(page);
  expect(live.legend.map((row) => row.caption), 'the deleted row is shown').toContain('repeat_region');
  const listed = await legendCaptions(page);
  expect(listed, 'the Legend list shows the deleted row').toContain('repeat_region');
  await generateAndWaitForResult(page);
  await settleLive(page);
  expect(diffSemanticSnapshots(live, await semanticSnapshot(page)).filter(({ kind }) => kind === 'legend'),
    'the Legend after Reset Settings equals the Legend Generate draws').toEqual([]);
  expect(listed, 'the Legend list after Reset Settings').toEqual(await legendCaptions(page));
});

// R15-2 (OV-241 decision B) retired Add legend item. A Session main saved with
// a Legend row the editor added (MANUAL_ROW, no features, in its Linear Result;
// tests/fixtures/sessions/two-mode-project.provenance.json) loads and shows the
// row, and Generate keeps it.
test('R15-2: a Legend row the editor added in a Session main saved loads, shows, and survives Generate (linear)', async ({ page }) => {
  test.setTimeout(240_000);
  await openFresh(page);
  await page.locator('input[accept^=".json,"]').setInputFiles('tests/fixtures/sessions/two-mode-project.v44.gbdraw-session.json.gz');
  await page.waitForFunction(() => !window.__GBDRAW_APP__.sessionImportPending, null, { timeout: 180_000 });
  expect(await page.evaluate(() => window.__GBDRAW_APP__.errorLog)).toBeNull();
  if (await page.evaluate(() => window.__GBDRAW_APP__.mode) !== 'linear') {
    await page.getByRole('button', { name: 'Linear', exact: true }).click();
  }
  await settleLive(page);
  const manualRow = () => page.evaluate(() => {
    const app = window.__GBDRAW_APP__;
    const entry = app.legendEntries.find((item) => item.caption === 'MANUAL_ROW');
    const shown = [...app.svgContainer.querySelectorAll('g[data-legend-key="MANUAL_ROW"]')]
      .filter((row) => !row.closest('[display="none"]'));
    return {
      entry: entry ? { color: entry.color, featureIds: entry.featureIds } : null,
      fills: shown.map((row) => row.querySelector('path, rect')?.getAttribute('fill'))
    };
  });
  expect(await manualRow(), 'Load lists and shows the row').toEqual({ entry: { color: '#445566', featureIds: [] }, fills: ['#445566'] });
  await generate(page);
  expect(await manualRow(), 'Generate keeps the row').toEqual({ entry: { color: '#445566', featureIds: [] }, fills: ['#445566'] });
});


// OV-243 (U3b B2): the Legend list keys a row by Python's key, not by its
// caption and color, so a renamed row whose color a Generate changes stays the
// rename of Python's row (R2: no Generate prunes the rename), and a later
// Legend edit of the row reaches it.
test('OV-243: a renamed GC skew row keeps its key through a palette change and Generate, and its delete hides it (circular)', async ({ page }) => {
  test.setTimeout(240_000);
  await open(page, { mode: 'circular', results: 'single', reflow: 'off' });
  await renameRow(page, 'GC skew (+)', 'Skew+');
  await settleLive(page);
  await switchPalette(page, 'arctic');
  const shownRows = async () => (await semanticSnapshot(page)).legend.map((row) => row.caption);
  for (const run of ['first', 'second']) {
    await generate(page);
    expect(await shownRows(), `the ${run} Generate draws the rename`).toContain('Skew+');
  }
  await deleteLegendRow(page, 'Skew+');
  expect(await shownRows(), 'the delete hides the row live').not.toContain('Skew+');
  expect(await page.evaluate(() => window.__GBDRAW_APP__.deletedLegendEntries.map((entry) => [entry.caption, entry.originalCaption])),
    'the deleted row is the rename of Python\'s row').toEqual([['Skew+', 'GC skew (+)']]);
  await expectLiveEqualsGenerate(page, { label: 'delete of the renamed GC skew row' });
});

// OV-288 (R15-4): a Legend row's stroke reaches the features the renderer
// draws for the row, in Python's row color, whatever color the row's swatch
// takes. Both edit orders, with Undo and Redo of the second edit, and Save
// then Load, show what Generate draws: the row's features keep the stroke.
const ROW_STROKE = '#e63946';
// forced_label_underlay.gb: the repeat_region feature and the CDS FL2.
const REPEAT_REGION = 'f9cf91913';
const RULE_FEATURE = 'f841fb8a8';
const rowStrokedFeatures = (snapshot) => Object.entries(snapshot.features)
  .filter(([, feature]) => feature.stroke.includes(ROW_STROKE)).map(([id]) => id);
for (const order of ['stroke then color', 'color then stroke']) {
  test(`OV-288: a Legend row stroke and color on CDS (${order}) stroke the row's features, live = Generate (circular)`, async ({ page, browser }, testInfo) => {
    test.setTimeout(300_000);
    await open(page, { mode: 'circular', results: 'single', reflow: 'off' });
    const stroke = () => legendRowStrokeColor(page, 'CDS', ROW_STROKE);
    const color = () => legendRowColor(page, 'CDS', '#7b2cbf');
    if (order === 'stroke then color') { await stroke(); await color(); } else { await color(); await stroke(); }
    const stroked = rowStrokedFeatures(await semanticSnapshot(page));
    expect(stroked.length, 'the stroke reaches the CDS features live').toBeGreaterThan(0);
    await history(page, 'undo');
    await history(page, 'redo');
    expect(rowStrokedFeatures(await semanticSnapshot(page)), 'Undo and Redo of the second edit').toEqual(stroked);
    const { generated } = await expectLiveEqualsGenerate(page, { label: `${order}, after Undo and Redo` });
    expect(rowStrokedFeatures(generated), 'Generate strokes the same features').toEqual(stroked);
    if (order !== 'stroke then color') return;
    const saved = testInfo.outputPath('ov288.gbdraw-session.json');
    await download(page, 'Save Session', saved, { title: '' });
    const loaded = await load(browser, saved);
    expect(rowStrokedFeatures(await semanticSnapshot(loaded)), 'Load shows the stroke').toEqual(stroked);
    await expectLiveEqualsGenerate(loaded, { label: `${order}, after Save and Load` });
    await loaded.context().close();
  });
}

// OV-288, OV-291: a rule whose caption names a row of another color draws its
// own row "<caption> [<hex>]" (N-06). The draft's color of a row is the color
// of the rule Python draws that row for, so a palette change shown live keeps
// the row the caption names in its palette color, and that row's stroke
// reaches the features drawn in it, live as at Generate. The rule's feature
// is drawn in the rule's own row, which the stroke does not reach (OV-292 A,
// R15-6).
test('OV-291: a palette change and a stroke on a row a rule\'s caption names, live = Generate (circular)', async ({ page }) => {
  test.setTimeout(240_000);
  await open(page, { mode: 'circular', results: 'single', reflow: 'off' });
  await addColorRule(page, { feat: 'CDS', qual: 'locus_tag', val: '^FL2$', color: '#2266aa', cap: 'repeat_region' });
  await generate(page);
  await switchPalette(page, 'arctic');
  await legendRowStrokeColor(page, 'repeat_region', ROW_STROKE);
  const stroked = rowStrokedFeatures(await semanticSnapshot(page));
  expect(stroked, 'the repeat_region feature, not the rule\'s feature (OV-292)').toEqual([REPEAT_REGION]);
  const { generated } = await expectLiveEqualsGenerate(page, { label: 'palette change, then a stroke on the row a rule caption names' });
  expect(rowStrokedFeatures(generated), 'Generate strokes the same feature').toEqual(stroked);
});

// Review H1 of OV-288 (R15-4): N-06 allocates a rule's row by the colors
// Python gives its rows, not by their swatches. A Legend color equal to the
// rule's color on the row the rule's caption names leaves the rule in its own
// row, so that row's stroke still reaches the features drawn in it, live as
// at Generate.
test('OV-288: a Legend color equal to a rule\'s color on the row its caption names keeps the row\'s stroke reach, live = Generate (circular)', async ({ page }) => {
  test.setTimeout(240_000);
  await open(page, { mode: 'circular', results: 'single', reflow: 'off' });
  await addColorRule(page, { feat: 'CDS', qual: 'locus_tag', val: '^FL2$', color: '#2266aa', cap: 'repeat_region' });
  await generate(page);
  await legendRowColor(page, 'repeat_region', '#2266aa');
  await legendRowStrokeColor(page, 'repeat_region', ROW_STROKE);
  const stroked = rowStrokedFeatures(await semanticSnapshot(page));
  expect(stroked, 'the repeat_region feature, not the rule\'s feature (OV-292)').toEqual([REPEAT_REGION]);
  const { generated } = await expectLiveEqualsGenerate(page, { label: 'a Legend color equal to the rule color, then a stroke' });
  expect(rowStrokedFeatures(generated), 'Generate strokes the same features').toEqual(stroked);
});

// OV-294 A (R15-6): the Legend editor finds the rules a row draws by Python's
// rows, never by the row's swatch. A Legend color on the row a rule's caption
// names, equal to the rule's color, leaves the rule in its own row, so a second
// Legend color edit of that row recolors the row, not the rule's feature.
test('OV-294: a second Legend color on the row a rule\'s caption names recolors the row, not the rule, live = Generate (circular)', async ({ page }) => {
  test.setTimeout(240_000);
  await open(page, { mode: 'circular', results: 'single', reflow: 'off' });
  await addColorRule(page, { feat: 'CDS', qual: 'locus_tag', val: '^FL2$', color: '#2266aa', cap: 'repeat_region' });
  await generate(page);
  await legendRowColor(page, 'repeat_region', '#2266aa');
  await legendRowColor(page, 'repeat_region', '#ff8800');
  const shown = async (snapshot) => ({
    rules: await page.evaluate(() => window.__GBDRAW_APP__.manualSpecificRules.map((rule) => rule.color)),
    ruleFeature: snapshot.features[RULE_FEATURE].fill,
    rows: Object.fromEntries(snapshot.legend.filter((row) => row.key.startsWith('repeat_region')).map((row) => [row.key, row.fill]))
  });
  const expected = { rules: ['#2266aa'], ruleFeature: ['#2266aa'], rows: { repeat_region: '#ff8800', 'repeat_region [#2266aa]': '#2266aa' } };
  expect(await shown(await semanticSnapshot(page)), 'the second edit recolors the row').toEqual(expected);
  await history(page, 'undo');
  // The first Legend color gives the row the rule row's name and swatch, so
  // the two are one row again (D-28): the rule's row draws it.
  expect(await shown(await semanticSnapshot(page)), 'Undo returns the row\'s first Legend color').toEqual({
    ...expected, rows: { 'repeat_region [#2266aa]': '#2266aa' }
  });
  await history(page, 'redo');
  expect(await shown(await semanticSnapshot(page)), 'Redo recolors the row again').toEqual(expected);
  const { generated } = await expectLiveEqualsGenerate(page, { label: 'a second Legend color on the row a rule caption names, after Undo and Redo' });
  expect(await shown(generated), 'Generate draws the same').toEqual(expected);
});

// OV-307: a row Python drew for a rule keeps the rule's color as Python drew
// it until Generate, so the rule that a Legend color edit recolored live still
// draws that row: a second Legend color edit recolors the rule again, and the
// row's stroke reaches the rule's feature, live as at Generate.
test('OV-307: two Legend color edits and a stroke on a rule\'s own row recolor and stroke the rule\'s feature, live = Generate (circular)', async ({ page, browser }, testInfo) => {
  test.setTimeout(300_000);
  await open(page, { mode: 'circular', results: 'single', reflow: 'off' });
  await addColorRule(page, { feat: 'CDS', qual: 'locus_tag', val: '^FL2$', color: '#2266aa', cap: 'Group' });
  await generate(page);
  await legendRowColor(page, 'Group', '#ff8800');
  await legendRowColor(page, 'Group', '#11aa55');
  await legendRowStrokeColor(page, 'Group', ROW_STROKE);
  const shown = async (snapshot, target = page) => ({
    rules: await target.evaluate(() => window.__GBDRAW_APP__.manualSpecificRules.map((rule) => rule.color)),
    ruleFeature: snapshot.features[RULE_FEATURE].fill,
    stroked: rowStrokedFeatures(snapshot),
    row: snapshot.legend.filter((row) => row.key === 'Group').map((row) => row.fill)
  });
  const expected = { rules: ['#11aa55'], ruleFeature: ['#11aa55'], stroked: [RULE_FEATURE], row: ['#11aa55'] };
  expect(await shown(await semanticSnapshot(page)), 'both edits recolor the rule').toEqual(expected);
  await history(page, 'undo');
  await history(page, 'undo');
  expect(await shown(await semanticSnapshot(page)), 'two Undo steps return the first edit\'s rule color, unstroked').toEqual({
    rules: ['#ff8800'], ruleFeature: ['#ff8800'], stroked: [], row: ['#ff8800']
  });
  await history(page, 'redo');
  await history(page, 'redo');
  expect(await shown(await semanticSnapshot(page)), 'two Redo steps return both edits').toEqual(expected);
  const { generated } = await expectLiveEqualsGenerate(page, { label: 'two Legend colors and a stroke on a rule row, after Undo and Redo' });
  expect(await shown(generated), 'Generate draws the same').toEqual(expected);
  const saved = testInfo.outputPath('ov307.gbdraw-session.json');
  await download(page, 'Save Session', saved, { title: '' });
  const loaded = await load(browser, saved);
  expect(await shown(await semanticSnapshot(loaded), loaded), 'Load shows the same').toEqual(expected);
  await expectLiveEqualsGenerate(loaded, { label: 'two Legend colors and a stroke on a rule row, after Save and Load' });
  await loaded.context().close();
});

// OV-293 (PD-OI-066): Generate draws no stroke for a deleted Legend row, so
// its delete takes the stroke off its features live; Undo, Redo and Restore
// show the stroke as Generate draws it. A batch Result displayed after the
// delete shows it too.
test('OV-293: deleting a stroked Legend row takes its stroke off live, and Undo, Redo and Restore follow, live = Generate (circular)', async ({ page }) => {
  test.setTimeout(300_000);
  await open(page, { mode: 'circular', results: 'single', reflow: 'off' });
  await legendRowStrokeColor(page, 'CDS', ROW_STROKE);
  const stroked = rowStrokedFeatures(await semanticSnapshot(page));
  expect(stroked.length, 'the stroke reaches the CDS features live').toBeGreaterThan(0);
  await deleteLegendRow(page, 'CDS');
  expect(rowStrokedFeatures(await semanticSnapshot(page)), 'the delete takes the stroke off').toEqual([]);
  await history(page, 'undo');
  expect(rowStrokedFeatures(await semanticSnapshot(page)), 'Undo of the delete').toEqual(stroked);
  await history(page, 'redo');
  expect(rowStrokedFeatures(await semanticSnapshot(page)), 'Redo of the delete').toEqual([]);
  await expectLiveEqualsGenerate(page, { label: 'delete of a stroked row' });
  await appAction(page, 'restoreAllDeletedLegendEntries');
  expect(rowStrokedFeatures(await semanticSnapshot(page)), 'Restore returns the stroke').toEqual(stroked);
  await expectLiveEqualsGenerate(page, { label: 'Restore of a stroked row' });
});

test('OV-293: a batch Result displayed after a stroked Legend row is deleted shows no stroke, live = Generate (circular, two-Result batch)', async ({ page }) => {
  test.setTimeout(300_000);
  await open(page, { mode: 'circular', results: 'batch', reflow: 'off' });
  await legendRowStrokeColor(page, 'CDS', ROW_STROKE);
  await showResult(page, 1);
  expect(rowStrokedFeatures(await semanticSnapshot(page)).length, 'Result 2 shows the stroke').toBeGreaterThan(0);
  await showResult(page, 0);
  await deleteLegendRow(page, 'CDS');
  await showResult(page, 1);
  expect(rowStrokedFeatures(await semanticSnapshot(page)), 'Result 2 shows the delete').toEqual([]);
  await expectLiveEqualsGenerate(page, { label: 'a batch Result displayed after the delete of a stroked row' });
});

// One naming rule on every surface (PD-OI-042 scenario revision 2, PD-OI-095).
// After the Legend operations of a pass the executor writes every shown row's
// text from the Web mirror of the rule (services/legend-row-names.js over the
// names and swatches the rows show), also on rows no operation touched. So on a
// Result as Python drew it the mirror must answer what Python wrote: the same
// texts, in the order the Legend draws its rows (the solid rows, then the
// gradient rows). Each figure is drawn here by the repository's Python and read
// in a browser DOM the way the executor reads a Result.
const REPO_ROOT = path.resolve(__dirname, '..', '..');
const FORCED_TO_LOCTEST_BLAST = [['FORCEDLBL', 'LOCTEST', '85.0', '300', '0', '0', '101', '400', '301', '600', '1e-50', '200']];
const RING_BLAST = [['RING', 'FORCEDLBL', '85.0', '300', '0', '0', '1', '300', '101', '400', '1e-50', '200']];
const REGION_COLUMNS = ['set_id', 'id', 'mark', 'record', 'start', 'end', 'label', 'legend_label', 'fill', 'hatch_angle', 'hatch_spacing', 'hatch_color'];
const region = (id, start, legendLabel, fill, hatch = ['', '', '']) => (
  ['regions', id, 'band', 'FORCEDLBL', String(start), String(start + 500), id, legendLabel, fill, ...hatch]
);
const LINEAR_PAIR = [SINGLE_FIXTURE, 'tests/fixtures/feature_location_search.gb'];
// `rows`: each row Python drew, in document order, as [key, data-legend-name (null: the key is the name), text].
const PYTHON_NAMED_RESULTS = [{
  name: 'a hatch row and a plain row of one name',
  mode: 'circular',
  inputs: [SINGLE_FIXTURE],
  files: { 'regions.tsv': [REGION_COLUMNS, region('a', 100, 'Island', '#88aaff', ['45', '4', '#202020']), region('b', 1000, 'Island', '#88aaff'), region('c', 2000, 'Plain', '#ffaa88')] },
  args: ['--annotation_table', 'regions.tsv'],
  rows: [
    ['CDS', null, 'CDS'], ['repeat_region', null, 'repeat_region'],
    ['GC content', null, 'GC content'], ['GC skew (+)', null, 'GC skew (+)'], ['GC skew (-)', null, 'GC skew (-)'],
    ['Island', null, 'Island (1)'], ['Island (2)', 'Island', 'Island (2)'], ['Plain', null, 'Plain']
  ]
}, {
  name: 'a gradient row in the two copies of a Linear Legend',
  mode: 'linear',
  inputs: LINEAR_PAIR,
  files: { 'pair.blast': FORCED_TO_LOCTEST_BLAST },
  args: ['-b', 'pair.blast'],
  rows: [['CDS', null, 'CDS'], ['repeat_region', null, 'repeat_region'], ['tRNA', null, 'tRNA'], ['Pairwise match identity', null, 'Pairwise match identity']]
}, {
  name: 'a solid row and a gradient row of one name (Linear, a region annotation)',
  mode: 'linear',
  inputs: LINEAR_PAIR,
  files: { 'pair.blast': FORCED_TO_LOCTEST_BLAST, 'regions.tsv': [REGION_COLUMNS, region('a', 100, 'Pairwise match identity', '#88aaff')] },
  args: ['-b', 'pair.blast', '--annotation_table', 'regions.tsv'],
  // The table holds the gradient row first; the Legend draws it last.
  rows: [
    ['CDS', null, 'CDS'], ['repeat_region', null, 'repeat_region'], ['tRNA', null, 'tRNA'],
    ['Pairwise match identity (2)', 'Pairwise match identity', 'Pairwise match identity (1)'],
    ['Pairwise match identity', null, 'Pairwise match identity (2)']
  ]
}, {
  name: 'a solid row and a gradient row of one name (Circular, a color rule)',
  mode: 'circular',
  inputs: [SINGLE_FIXTURE],
  files: { 'ring.blast': RING_BLAST, 'rules.tsv': [['CDS', 'locus_tag', 'FL1', '#ff0000', 'Conservation identity']] },
  args: ['--conservation_blast', 'ring.blast', '-t', 'rules.tsv'],
  rows: [
    ['Conservation identity', null, 'Conservation identity (1)'], ['other proteins', null, 'other proteins'], ['repeat_region', null, 'repeat_region'],
    ['GC content', null, 'GC content'], ['GC skew (+)', null, 'GC skew (+)'], ['GC skew (-)', null, 'GC skew (-)'],
    ['Conservation identity (2)', 'Conservation identity', 'Conservation identity (2)']
  ]
}, {
  name: 'a caption in two colors in the two copies of a Linear Legend',
  mode: 'linear',
  inputs: LINEAR_PAIR,
  files: {
    'pair.blast': FORCED_TO_LOCTEST_BLAST,
    'rules.tsv': [['CDS', 'locus_tag', 'FL1', '#ff0000', 'Group'], ['CDS', 'locus_tag', 'LOC_0001', '#00ff00', 'Group']]
  },
  args: ['-b', 'pair.blast', '-t', 'rules.tsv'],
  rows: [
    ['Group [#ff0000]', 'Group', 'Group (1)'], ['Group [#00ff00]', 'Group', 'Group (2)'], ['other proteins', null, 'other proteins'],
    ['repeat_region', null, 'repeat_region'], ['tRNA', null, 'tRNA'], ['Pairwise match identity', null, 'Pairwise match identity']
  ]
}, {
  name: 'a caption in two colors and a rule named like a feature type row',
  mode: 'circular',
  inputs: ['tests/fixtures/feature_location_search.gb'],
  files: { 'rules.tsv': [
    ['CDS', 'locus_tag', 'LOC_0001', '#ff0000', 'Group'], ['CDS', 'locus_tag', 'LOC_0002', '#00ff00', 'Group'],
    ['CDS', 'locus_tag', 'LOC_0003', '#0000ff', 'tRNA']
  ] },
  args: ['-t', 'rules.tsv'],
  rows: [
    ['Group [#ff0000]', 'Group', 'Group (1)'], ['Group [#00ff00]', 'Group', 'Group (2)'],
    ['tRNA [#0000ff]', 'tRNA', 'tRNA (1)'], ['tRNA', null, 'tRNA (2)'],
    ['GC content', null, 'GC content'], ['GC skew (+)', null, 'GC skew (+)'], ['GC skew (-)', null, 'GC skew (-)']
  ]
}, {
  name: 'tracks of one label in one color and in two, and a track named like a feature type row',
  mode: 'circular',
  inputs: [SINGLE_FIXTURE],
  args: [
    'features:features', 'gc_content:dinucleotide_content', 'gc_content_2:dinucleotide_content', 'gc_skew:dinucleotide_skew',
    'gc_skew_2:dinucleotide_skew@positive_color=#00ff00', 'cds_named:dinucleotide_skew@legend_label=CDS,positive_color=#123456,negative_color=#123456'
  ].flatMap((slot) => ['--circular_track_slot', slot]),
  rows: [
    ['CDS', null, 'CDS (1)'], ['repeat_region', null, 'repeat_region'], ['GC content', null, 'GC content'],
    ['GC skew (+)', null, 'GC skew (+) (1)'], ['GC skew (-)', null, 'GC skew (-)'],
    ['GC skew (+) (2)', 'GC skew (+)', 'GC skew (+) (2)'], ['CDS (2)', 'CDS', 'CDS (2)']
  ]
}];

const drawWithPython = (testInfo, { mode, inputs, files = {}, args = [] }) => {
  Object.entries(files).forEach(([name, rows]) => {
    writeFileSync(testInfo.outputPath(name), `${rows.map((row) => row.join('\t')).join('\n')}\n`);
  });
  execFileSync(process.env.PYTHON || 'python', [
    '-m', 'gbdraw.cli', mode, '--gbk', ...inputs.map((input) => path.resolve(REPO_ROOT, input)),
    ...args.map((arg) => (arg in files ? testInfo.outputPath(arg) : arg)), '-o', testInfo.outputPath('figure'), '-f', 'svg'
  ], { cwd: REPO_ROOT, stdio: 'pipe' });
  return readFileSync(testInfo.outputPath('figure.svg'), 'utf8');
};

// A document of the app's origin with its sanitizer and without the app: the
// executor is a module of DOM functions.
const EXECUTOR_PAGE = '/gbdraw/web/legend-names-executor.html';
const openExecutorPage = async (page) => {
  await page.route(`**${EXECUTOR_PAGE}`, (route) => route.fulfill({
    contentType: 'text/html', body: '<!doctype html><script src="/gbdraw/web/vendor/dompurify/purify.min.js"></script>'
  }));
  await page.goto(EXECUTOR_PAGE);
};

// The rows of each Legend of a Result (a Linear Result holds a horizontal and
// a vertical copy) as Python wrote them, as the admitted Result holds them,
// and after the executor's reconcile without an operation.
const namedByExecutor = (page, content) => page.evaluate(async (svgText) => {
  const { createEmptySvgMutationPlan, reconcileMountedResult } = await import('/gbdraw/web/js/services/svg-result-ingestion.js');
  const { sanitizeSvgContent } = await import('/gbdraw/web/js/services/svg-sanitization.js');
  const parse = (text) => new DOMParser().parseFromString(text, 'image/svg+xml').documentElement;
  const legends = (svg) => {
    const copies = Array.from(svg.querySelectorAll('#legend_horizontal, #legend_vertical'));
    return (copies.length ? copies : [svg.querySelector('#legend')]).map((legend) => (
      Array.from(legend.querySelectorAll('g[data-legend-key]')).map((row) => [
        row.getAttribute('data-legend-key'), row.getAttribute('data-legend-name'), row.querySelector('text').textContent
      ])
    ));
  };
  const result = parse(sanitizeSvgContent(svgText));
  const admitted = legends(result);
  const outcome = reconcileMountedResult(result, createEmptySvgMutationPlan(1).operationsByResult[0]);
  return { python: legends(parse(svgText)), admitted, outcome, named: legends(result) };
}, content);

for (const figure of PYTHON_NAMED_RESULTS) {
  test(`the naming pass keeps the names Python drew: ${figure.name} (${figure.mode})`, async ({ page }, testInfo) => {
    const content = drawWithPython(testInfo, figure);
    await openExecutorPage(page);
    const { python, admitted, outcome, named } = await namedByExecutor(page, content);
    expect(python, 'what Python drew').toEqual(Array(figure.mode === 'linear' ? 2 : 1).fill(figure.rows));
    expect(admitted, 'the sanitized Result keeps each row\'s key, name, and text').toEqual(python);
    expect(named, 'the naming pass rewrites no row').toEqual(python);
    expect(outcome).toEqual({ changed: false, legendChanged: false });
    expect(python.flat().filter(([, , text]) => text.includes(' [#')), 'no drawn name is a row key').toEqual([]);
  });
}

// The same rule while the Result is edited: the names the live Result draws
// are the names the next Generate draws, the Legend tab shows a row by its
// drawn name, and a row key (`name [fill]`, `label (2)`) is never a name.
// Rows are read from each Legend the Result holds, in document order, as
// [key, data-legend-name, text].
const legendNameFacts = (page) => page.evaluate(() => {
  const app = window.__GBDRAW_APP__;
  const svg = app.svgContainer.querySelector('svg');
  const copies = [...svg.querySelectorAll('#legend_horizontal, #legend_vertical')];
  return {
    legends: (copies.length ? copies : [svg.querySelector('#legend')]).map((legend) => (
      [...legend.querySelectorAll('g[data-legend-key]')].filter((row) => row.getAttribute('display') !== 'none').map((row) => [
        row.getAttribute('data-legend-key'), row.getAttribute('data-legend-name'), row.querySelector('text').textContent
      ])
    )),
    shown: app.legendEntries.map((entry) => app.shownLegendName(entry.caption))
  };
});
const expectLegendNames = async (page, names, label) => {
  const { legends, shown } = await legendNameFacts(page);
  for (const rows of legends) expect(rows.map(([, , text]) => text), `${label}: the Legend draws`).toEqual(names);
  expect(shown, `${label}: the Legend tab shows each row by its drawn name`).toEqual(names);
  expect([...legends.flat().map(([, , text]) => text), ...shown].filter((name) => name.includes(' [#')), `${label}: no name is a row key`).toEqual([]);
  return legends[0];
};
const expectNamesLiveAndGenerated = async (page, names, label) => {
  await expectLegendNames(page, names, `${label}, live`);
  const { generated } = await expectLiveEqualsGenerate(page, { label });
  await expectLegendNames(page, names, `${label}, Generate`);
  return generated;
};
const setRuleField = (page, index, field, value) => appAction(page, 'setSpecificRuleField', index, field, value);

test('legend names: a color rule named like a feature type row is numbered with it, and shares its row in its color, live = Generate (circular)', async ({ page }) => {
  test.setTimeout(300_000);
  await openWithGenBank(page, LOCTEST_FIXTURE, () => { window.__GBDRAW_APP__.form.labels_mode = 'out'; });
  await generate(page);
  const tracks = ['GC content', 'GC skew (+)', 'GC skew (-)'];
  const typeFill = (await semanticSnapshot(page)).legend.find(({ key }) => key === 'tRNA').fill;
  await addColorRule(page, { feat: 'CDS', qual: 'locus_tag', val: '^LOC_0001$', color: '#0000ff', cap: 'tRNA' });
  await expectNamesLiveAndGenerated(page, ['tRNA (1)', 'other proteins', 'tRNA (2)', ...tracks], 'a rule named tRNA in another color');
  expect(await ruleTable(page), 'the caption stays as written').toEqual([['^LOC_0001$', 'tRNA', '#0000ff']]);
  await setRuleField(page, 0, 'color', typeFill);
  await expectNamesLiveAndGenerated(page, ['tRNA', 'other proteins', ...tracks], 'a rule named tRNA in the color of the tRNA row');
});

// A track of one color labelled `GC content`: a second source of the GC
// content row's label, in its own color.
const GC_CONTENT_LABELLED = { legend_label: 'GC content', positive_color: '#123456', negative_color: '#123456' };
const addCircularSlot = (page, renderer, params = {}) => page.evaluate(({ kind, values }) => {
  const app = window.__GBDRAW_APP__;
  app.setCircularTrackSlotsEnabled(true);
  const before = new Set(app.adv.circular_track_slots.map((slot) => slot.id));
  app.addCircularTrackSlot(kind);
  const slot = app.adv.circular_track_slots.find((entry) => !before.has(entry.id));
  Object.assign(slot.params, values);
  return slot.id;
}, { kind: renderer, values: params });
const legendKeys = async (page) => (await legendNameFacts(page)).legends[0].map(([key]) => key);

test('legend names: tracks of one label share a row in one color and are numbered in two, live = Generate (circular)', async ({ page }) => {
  test.setTimeout(300_000);
  await open(page, { mode: 'circular', results: 'single', reflow: 'off' });
  const names = (...gcContent) => ['CDS', 'repeat_region', gcContent[0], 'GC skew (+)', 'GC skew (-)', ...gcContent.slice(1)];
  // Two sources of one color: one row. A Legend edit of another row runs the naming pass on the live Result.
  const second = await addCircularSlot(page, 'dinucleotide_content');
  await generate(page);
  await colorLegendRow(page, 'CDS', '#264653');
  await expectNamesLiveAndGenerated(page, names('GC content'), 'two tracks of one label and one color');
  // Three sources, two of one color: the third keeps its key and is drawn as the second row of the name.
  await addCircularSlot(page, 'dinucleotide_skew', GC_CONTENT_LABELLED);
  await generate(page);
  expect(await legendKeys(page), 'a source that shares a row still takes its key').toEqual(names('GC content', 'GC content (3)'));
  await colorLegendRow(page, 'CDS', '#2a9d8f');
  await expectNamesLiveAndGenerated(page, names('GC content (1)', 'GC content (2)'), 'three sources of one label, two of one color');
  // Two sources in two colors.
  await page.evaluate((id) => {
    const app = window.__GBDRAW_APP__;
    app.setCircularTrackSlotEnabled(app.adv.circular_track_slots.find((slot) => slot.id === id), false);
  }, second);
  await generate(page);
  expect(await legendKeys(page)).toEqual(names('GC content', 'GC content (2)'));
  await colorLegendRow(page, 'GC content (2)', '#7b2cbf');
  await expectNamesLiveAndGenerated(page, names('GC content (1)', 'GC content (2)'), 'two tracks of one label in two colors');
});

// A color rule's caption written again in the rule list is kept as written
// (the Legend editor's rename onto a taken name asks instead, PD-OI-061).
test('legend names: a rule renamed onto the name of another rule numbers both rows, and renamed away leaves the bare name, live = Generate (linear)', async ({ page }) => {
  test.setTimeout(300_000);
  await open(page, { mode: 'linear', results: 'single', reflow: 'off' });
  for (const rule of RULE_ROWS.slice(0, 2)) await addColorRule(page, rule);
  const rest = ['other proteins', 'repeat_region'];
  await expectLegendNames(page, ['alpha', 'beta', ...rest], 'two rules of two names');
  await setRuleField(page, 1, 'cap', 'alpha');
  await expectNamesLiveAndGenerated(page, ['alpha (1)', 'alpha (2)', ...rest], 'a rule renamed onto the name of another rule');
  expect(await ruleTable(page), 'both captions stay as written').toEqual([['^FL1$', 'alpha', '#e63946'], ['^FL2$', 'alpha', '#2a9d8f']]);
  await setRuleField(page, 1, 'cap', 'gamma');
  await expectNamesLiveAndGenerated(page, ['alpha', 'gamma', ...rest], 'the rule renamed away');
});

// OV-645: the key of a rule's row changes when its name starts or stops being
// written in several colors (`alpha` and `alpha [#e63946]`). The row is the
// same row: it keeps its place when the Legend editor renames the other row of
// its name, and, in an order the user set, when the rule list makes the name
// shared and single again. A row the rule list renames is drawn anew, last.
test('legend names: a rule row keeps its place when its key changes with the sharing of its name, also in a user order, live = Generate (linear)', async ({ page }) => {
  test.setTimeout(300_000);
  await open(page, { mode: 'linear', results: 'single', reflow: 'off' });
  for (const rule of RULE_ROWS.slice(0, 2)) await addColorRule(page, rule);
  await setRuleField(page, 1, 'cap', 'alpha');
  await expectLegendNames(page, ['alpha (1)', 'alpha (2)', 'other proteins', 'repeat_region'], 'two rules of one name');
  await renameAndSettle(page, 'alpha [#2a9d8f]', 'gamma');
  await expectNamesLiveAndGenerated(page, ['alpha', 'gamma', 'other proteins', 'repeat_region'], 'the Legend editor renames one of two rows of a name');
  expect(await ruleTable(page)).toEqual([['^FL1$', 'alpha', '#e63946'], ['^FL2$', 'gamma', '#2a9d8f']]);

  await appAction(page, 'moveLegendEntry', 0, 2);
  await expectLegendNames(page, ['gamma', 'other proteins', 'alpha', 'repeat_region'], 'the user moves the alpha row');
  await setRuleField(page, 1, 'cap', 'alpha');
  await expectNamesLiveAndGenerated(page, ['other proteins', 'alpha (1)', 'repeat_region', 'alpha (2)'], 'the name becomes shared in a user order');
  await setRuleField(page, 1, 'cap', 'beta');
  await expectNamesLiveAndGenerated(page, ['other proteins', 'alpha', 'repeat_region', 'beta'], 'the name stops being shared in a user order');
});

// OV-1235: a rule color that writes a caption in a second color splits its
// row in two (`alpha [#e63946]`, `alpha [#2a9d8f]`), which Python places in
// rule order; the live edit asks Python, as for any row that loses part of its
// features (OV-43).
test('OV-1235: a rule color that splits the row of its caption draws the rows in rule order, live = Generate (linear)', async ({ page }) => {
  test.setTimeout(300_000);
  await open(page, { mode: 'linear', results: 'single', reflow: 'off' });
  await addColorRule(page, RULE_ROWS[0]);
  await addColorRule(page, { ...RULE_ROWS[1], cap: 'alpha', color: RULE_ROWS[0].color });
  await expectLegendNames(page, ['alpha', 'other proteins', 'repeat_region'], 'two rules of one name and one color');
  await setRuleField(page, 1, 'color', RULE_ROWS[1].color);
  await expectNamesLiveAndGenerated(page, ['alpha (1)', 'alpha (2)', 'other proteins', 'repeat_region'], 'the second rule in another color');
});

// OV-1236: Python draws no row for a color of a caption that no drawn feature
// uses, so hiding the only feature of `alpha [#2a9d8f]` leaves one row, named
// `alpha`; the live edit asks Python (OV-42).
test('OV-1236: hiding the only feature of one color of a shared caption leaves one row of the caption, live = Generate (linear)', async ({ page }) => {
  test.setTimeout(300_000);
  await open(page, { mode: 'linear', results: 'single', reflow: 'off' });
  await addColorRule(page, RULE_ROWS[0]);
  await addColorRule(page, { ...RULE_ROWS[1], cap: 'alpha' });
  await expectLegendNames(page, ['alpha (1)', 'alpha (2)', 'other proteins', 'repeat_region'], 'two rules of one name');
  await addVisibilityRule(page, { ...FL1_OFF, value: '^fl2$' });
  await expectNamesLiveAndGenerated(page, ['alpha', 'other proteins', 'repeat_region'], 'FL2 hidden');
});

// OV-1486 (PD-OI-063, D-08): Generate replays a sorted order after a color
// rule's removal brings back a row the order never held: the sorted rows keep
// their places and the returning row follows them, live and at Generate. The
// inventory a drawn Result leaves is Python's order of its rows, which the
// executor recorded before it replayed the sorted order, not the order shown.
test('OV-1486: a sorted Legend order keeps its places when removing a color rule brings back the CDS row, live = Generate (linear)', async ({ page }) => {
  test.setTimeout(300_000);
  await open(page, { mode: 'linear', results: 'single', reflow: 'off' });
  await addColorRule(page, RULE_ROWS[0]);
  await generate(page);
  await expectLegendNames(page, ['alpha', 'other proteins', 'repeat_region'], 'a rule takes one CDS');
  await appAction(page, 'sortLegendEntries', 'desc');
  await expectNamesLiveAndGenerated(page, ['repeat_region', 'other proteins', 'alpha'], 'the Legend sorted Z to A');
  await appAction(page, 'removeSpecificRule', 0);
  await expectNamesLiveAndGenerated(page, ['repeat_region', 'CDS'], 'the rule removed in the sorted order');
});

// OV-660: a stroke stored on a rule row follows the row when its key changes
// with the sharing of its name (`alpha [#e63946]` to `alpha` and back), live
// and at Generate.
test('legend names: a stroke stored on a rule row follows the row when its key changes, live = Generate (linear)', async ({ page }) => {
  test.setTimeout(300_000);
  await open(page, { mode: 'linear', results: 'single', reflow: 'off' });
  for (const rule of RULE_ROWS.slice(0, 2)) await addColorRule(page, rule);
  await setRuleField(page, 1, 'cap', 'alpha');
  await page.evaluate(() => {
    const app = window.__GBDRAW_APP__;
    const index = app.legendEntries.findIndex((entry) => entry.caption === 'alpha [#e63946]');
    return app.updateLegendEntryStrokeColor(index, '#111111') && app.updateLegendEntryStrokeWidth(index, 3);
  });
  await settleLive(page);
  const redRow = (snapshot) => snapshot.legend.filter(({ fill }) => fill.toLowerCase() === '#e63946')
    .map(({ key, stroke }) => [key, stroke.toLowerCase()]);
  const rest = ['other proteins', 'repeat_region'];
  await setRuleField(page, 1, 'cap', 'gamma');
  expect(redRow(await semanticSnapshot(page)), 'the other rule renamed away, live').toEqual([['alpha', '#111111']]);
  expect(redRow(await expectNamesLiveAndGenerated(page, ['alpha', 'gamma', ...rest], 'the other rule renamed away')))
    .toEqual([['alpha', '#111111']]);
  await setRuleField(page, 1, 'cap', 'alpha');
  expect(redRow(await semanticSnapshot(page)), 'the name shared again, live').toEqual([['alpha [#e63946]', '#111111']]);
  expect(redRow(await expectNamesLiveAndGenerated(page, ['alpha (1)', 'alpha (2)', ...rest], 'the name shared again')))
    .toEqual([['alpha [#e63946]', '#111111']]);
});

// A row key keeps the string 0.13.0 gave the second source of a label
// (`GC content (2)`), so the edits a Session stores under it stay on that
// track's row. Once both tracks have one swatch they share the first row, and
// the edits of the second are not used (PD-OI-095, accepted residual risk).
test('legend names: a color and a stroke stored under GC content (2) stay on the second track\'s row, and are not used once the tracks share a row (circular)', async ({ page }) => {
  test.setTimeout(300_000);
  await open(page, { mode: 'circular', results: 'single', reflow: 'off' });
  const names = (...gcContent) => ['CDS', 'repeat_region', gcContent[0], 'GC skew (+)', 'GC skew (-)', ...gcContent.slice(1)];
  const labelled = await addCircularSlot(page, 'dinucleotide_skew', GC_CONTENT_LABELLED);
  await generate(page);
  const firstFill = (await semanticSnapshot(page)).legend.find(({ key }) => key === 'GC content').fill;
  await colorLegendRow(page, 'GC content (2)', '#7b2cbf');
  await page.evaluate(() => {
    const app = window.__GBDRAW_APP__;
    const index = app.legendEntries.findIndex((entry) => entry.caption === 'GC content (2)');
    return app.updateLegendEntryStrokeColor(index, '#e63946') && app.updateLegendEntryStrokeWidth(index, 2);
  });
  await settleLive(page);
  expect(await page.evaluate(async () => {
    const drawing = (await import('/gbdraw/web/js/state.js')).state.activeDrawing();
    return [Object.keys(drawing.legendColorOverrides), Object.keys(drawing.legendStrokeOverrides)];
  }), 'the edits are stored under the row key').toEqual([['GC content (2)'], ['GC content (2)']]);
  const generated = await expectNamesLiveAndGenerated(page, names('GC content (1)', 'GC content (2)'), 'edits on the second track\'s row');
  const edited = (row) => [row.key, row.caption, row.fill.toLowerCase(), row.stroke.toLowerCase()];
  expect(generated.legend.filter(({ key }) => key.startsWith('GC content')).map(edited), 'Generate keeps the edits on the second track\'s row').toEqual([
    ['GC content', 'GC content (1)', firstFill.toLowerCase(), 'none'], ['GC content (2)', 'GC content (2)', '#7b2cbf', '#e63946']
  ]);
  // The second track in the color of the first: one row, drawn as Python draws the first.
  await page.evaluate(({ id, color }) => {
    const slot = window.__GBDRAW_APP__.adv.circular_track_slots.find((entry) => entry.id === id);
    Object.assign(slot.params, { positive_color: color, negative_color: color });
  }, { id: labelled, color: firstFill });
  await generate(page);
  await expectLegendNames(page, names('GC content'), 'both tracks in one color');
  expect((await semanticSnapshot(page)).legend.filter(({ key }) => key.startsWith('GC content')).map(edited), 'the shared row has no edit of the second')
    .toEqual([['GC content', 'GC content', firstFill.toLowerCase(), 'none']]);
});

// Undo and Redo of a color rule commit: the rule's row returns with its name
// beside its key, and the Legend draws the same names.
test('legend names: Undo and Redo of a color rule commit keep the name of the rule\'s row and the drawn names (circular)', async ({ page }) => {
  test.setTimeout(300_000);
  await open(page, { mode: 'circular', results: 'single', reflow: 'off' });
  const tracks = ['GC content', 'GC skew (+)', 'GC skew (-)'].map((name) => [name, null, name]);
  const rows = async () => (await legendNameFacts(page)).legends[0];
  const python = await rows();
  expect(python).toEqual([['CDS', null, 'CDS'], ['repeat_region', null, 'repeat_region'], ...tracks]);
  const one = [
    ['repeat_region [#e63946]', 'repeat_region', 'repeat_region (1)'], ['other proteins', null, 'other proteins'],
    ['repeat_region', null, 'repeat_region (2)'], ...tracks
  ];
  const two = [one[0], ['repeat_region [#2a9d8f]', 'repeat_region', 'repeat_region (2)'], one[1], ['repeat_region', null, 'repeat_region (3)'], ...tracks];
  await addColorRule(page, { ...RULE_ROWS[0], cap: 'repeat_region' });
  expect(await rows(), 'a rule named like the repeat_region row').toEqual(one);
  await history(page, 'undo');
  expect(await rows(), 'Undo of the commit').toEqual(python);
  await history(page, 'redo');
  expect(await rows(), 'Redo of the commit').toEqual(one);
  await addColorRule(page, { ...RULE_ROWS[1], cap: 'repeat_region' });
  expect(await rows(), 'a second rule of the name').toEqual(two);
  await history(page, 'undo');
  expect(await rows(), 'Undo of the second commit').toEqual(one);
  await history(page, 'redo');
  expect(await rows(), 'Redo of the second commit').toEqual(two);
  await expectNamesLiveAndGenerated(page, two.map(([, , text]) => text), 'after Undo and Redo');
  expect(await rows(), 'Generate writes the same keys and names').toEqual(two);
});

// OV-1985 (release-tier random walk, seed 40634318370): a caption-scope stroke
// on a rule row whose key a later rule edit changes. The fill of FL2 "Use
// existing color" writes a CDS rule in the palette color, and a second CDS rule
// in another color makes Python draw the row as `CDS [#54bcf8]`; the stroke is
// stored under that key. A caption edit of the second rule leaves CDS in one
// color, so the commit moves the row and its stroke to the key `CDS` (OV-660);
// after Generate, the caption written back moves them to `CDS [#54bcf8]`
// again. Each time the live Result shows the stroke on the moved row, as
// Generate draws it.
for (const { name, generateBetween } of [
  { name: 'the row leaves its color ID', generateBetween: false },
  { name: 'the row takes its color ID back after Generate', generateBetween: true }
]) {
  test(`a Legend row stroke moves with its row when a rule caption edit re-keys it, live = Generate: ${name} (OV-1985)`, async ({ page }) => {
    test.setTimeout(240_000);
    await open(page, { mode: 'linear', results: 'single', reflow: 'on' });
    await popupEdit(page, { type: 'CDS', start: 2000, end: 2600 }, { fill: '#e63946', scope: 'useExisting' });
    await addColorRule(page, { feat: 'CDS', qual: 'locus_tag', val: '^FL2$', color: '#e63946', cap: 'CDS' });
    await popupEdit(page, { type: 'CDS', start: 2850, end: 2900 }, { stroke: '#123456', scope: 'caption' });
    const setCaption = async (value) => {
      await evaluateWithRetainedPromise(page, (caption) => {
        const app = window.__GBDRAW_APP__;
        return app.setSpecificRuleField(app.manualSpecificRules.findIndex((rule) => rule.val === '^FL2$'), 'cap', caption);
      }, value);
      await settleLive(page);
    };
    await setCaption('transfer RNA');
    if (generateBetween) {
      await generate(page);
      await setCaption('CDS');
    }
    const cdsRow = (snapshot) => snapshot.legend.find((row) => row.caption === 'CDS');
    expect(cdsRow(await semanticSnapshot(page))?.stroke, 'the live CDS row shows the stroke').toBe('#123456');
    await expectLiveEqualsGenerate(page, { label: `a rule caption edit after a caption-scope stroke: ${name}` });
  });
}

// OV-1986: a Legend row's stroke reaches the features drawn in the row's color,
// so a rule commit that changes which features a stroked row reaches shows the
// strokes too, live and on Undo and Redo, as Generate draws them. (A) The rule
// of a stroked row is removed: its feature leaves the row and its stroke. (B)
// A pattern edit moves FL2 from the stroked row `other proteins` to the rule
// row `Other`; no Legend row is added or removed. On dev before the fix both
// commits already asked for the automatic rerender, which redrew the strokes,
// so these cases pin the result rather than fail first.
const ruleIndexOf = (page, cap) => page.evaluate((caption) => (
  window.__GBDRAW_APP__.manualSpecificRules.findIndex((rule) => rule.cap === caption)
), cap);
test('a stroked rule row whose rule is removed leaves its feature unstroked, live = Generate (OV-1986)', async ({ page }) => {
  test.setTimeout(240_000);
  await open(page, { mode: 'linear', results: 'single', reflow: 'off' });
  await addColorRule(page, { feat: 'CDS', qual: 'locus_tag', val: '^FL2$', color: '#2266aa', cap: 'Special' });
  await legendRowStrokeColor(page, 'Special', '#123456');
  const index = await ruleIndexOf(page, 'Special');
  await evaluateWithRetainedPromise(page, (rule) => window.__GBDRAW_APP__.removeSpecificRule(rule), index);
  await settleLive(page);
  await expectLiveEqualsGenerate(page, { label: 'the rule of a stroked row removed' });
});
for (const { name, history: steps } of [
  { name: 'the pattern edit', history: [] },
  { name: 'the pattern edit, then Undo and Redo', history: ['undo', 'redo'] }
]) {
  test(`a rule pattern edit that moves a feature out of a stroked row takes the row's stroke from it, live = Generate: ${name} (OV-1986)`, async ({ page }) => {
    test.setTimeout(240_000);
    await open(page, { mode: 'linear', results: 'single', reflow: 'off' });
    await addColorRule(page, { feat: 'CDS', qual: 'locus_tag', val: '^FL1$', color: '#2266aa', cap: 'Other' });
    await legendRowStrokeColor(page, 'other proteins', '#123456');
    const index = await ruleIndexOf(page, 'Other');
    await evaluateWithRetainedPromise(page, (rule) => window.__GBDRAW_APP__.setSpecificRuleField(rule, 'val', '^FL[12]$'), index);
    await settleLive(page);
    for (const step of steps) await history(page, step);
    await expectLiveEqualsGenerate(page, { label: `a pattern edit after a Legend row stroke: ${name}` });
  });
}

// OV-1987: a rule that takes a CDS feature under another caption makes Python
// name the CDS palette row "other proteins", so the stroke edit of row `CDS`
// names a row the Result no longer draws, and Generate strokes no feature with
// it. A later rule edit shows the strokes live; the row it names is absent, so
// it reaches no feature either, although the draft still gives the key `CDS`
// its palette color.
test('a stroke on a Legend row the Result no longer draws reaches no feature after a rule edit, live = Generate (OV-1987)', async ({ page }) => {
  test.setTimeout(240_000);
  await open(page, { mode: 'linear', results: 'single', reflow: 'off' });
  await legendRowStrokeColor(page, 'CDS', '#7b2cbf');
  await addColorRule(page, { feat: 'CDS', qual: 'locus_tag', val: 'FL2$', color: '#f4a261', cap: 'transfer RNA' });
  const index = await ruleIndexOf(page, 'transfer RNA');
  await evaluateWithRetainedPromise(page, (rule) => window.__GBDRAW_APP__.setSpecificRuleField(rule, 'color', '#123456'), index);
  await settleLive(page);
  await expectLiveEqualsGenerate(page, { label: 'a rule edit after the stroked row left the Result' });
});
