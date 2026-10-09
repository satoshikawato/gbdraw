// PD-OI-066 (R1, R3): Legend edits on the live Result equal the next Generate:
// the automatic rerender a Legend source asks for (OV-35, OV-42 to OV-44),
// Legend colors on rows that the draft removes or that only one Result or mode
// draws (OV-63, OV-80), and Legend renames (OV-62). The matrix of edit kinds is
// in live-generate-parity.playwright.spec.js.
const { test, expect } = require('@playwright/test');
const { readFileSync } = require('node:fs');
const { evaluateWithRetainedPromise, generateAndWaitForResult } = require('./helpers/app-lifecycle.cjs');
const {
  diffSemanticSnapshots, expectLiveEqualsGenerate, semanticSnapshot, settleLive, showResult
} = require('./helpers/live-generate-parity.cjs');
const { loadSessionFile, openFresh } = require('./helpers/audit-browser.cjs');
const {
  SINGLE_FIXTURE, open, generate, popupEdit, addVisibilityRule, appAction, addColorRule, history,
  FL1_OFF, legendRowColor, FL1_ALPHA, DEPTH_TSV, colorLegendRow, openCanvas, renameRow,
  legendRowStrokeColor, legendRowAdd, switchPalette, deleteLegendRow
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
// which lays the Legend out and serializes the Result once; so does Add legend
// item. Counted: the Result serializations of the preview owner's commit
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

test('a color rule or Add legend item that adds a Legend row serializes the Result once', async ({ page }) => {
  test.setTimeout(120_000);
  await open(page, { mode: 'circular', results: 'single', reflow: 'off' });
  let serializations = await countResultSerializations(page);
  await addColorRule(page, FL1_ALPHA);
  expect(await page.evaluate(() => window.__GBDRAW_APP__.legendEntries.map((entry) => entry.caption))).toContain('alpha');
  expect(await serializations(), 'color rule commit that adds a Legend row').toBe(1);
  serializations = await countResultSerializations(page);
  await legendRowAdd(page, 'Added', '#123456');
  expect(await serializations(), 'Add legend item').toBe(1);
  await expectLiveEqualsGenerate(page, { label: 'Legend rows added by a rule and by Add legend item' });
});

// OV-60: renaming a generated Legend row (GC content) onto the caption of a
// color rule that draws no row is an explicit rename; Generate draws it too.
test('a Legend rename onto the caption of a color rule without a drawn row equals Generate', async ({ page }) => {
  test.setTimeout(120_000);
  await open(page, { mode: 'circular', results: 'single', reflow: 'on' });
  await addColorRule(page, { feat: 'CDS', qual: 'locus_tag', val: '^NOMATCH$', color: '#c83366', cap: 'Zeta' });
  await evaluateWithRetainedPromise(page, async () => {
    const app = window.__GBDRAW_APP__;
    await app.renameLegendEntry(app.legendEntries.findIndex((entry) => entry.caption === 'GC content'), 'Zeta');
  });
  await settleLive(page);
  const captions = await page.evaluate(() => window.__GBDRAW_APP__.legendEntries.map((entry) => entry.caption));
  expect(captions, 'the renamed row is drawn live').toContain('Zeta');
  await expectLiveEqualsGenerate(page, { label: 'Legend rename onto the caption of a color rule without a row' });
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
      await page.evaluate(() => { window.__GBDRAW_APP__.form.show_depth = false; });
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
      await page.evaluate(() => { window.__GBDRAW_APP__.form.show_depth = false; });
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
    await page.evaluate(() => { window.__GBDRAW_APP__.form.show_depth = false; });
    await settleLive(page);
    await generate(page);
    expect((await semanticSnapshot(page)).legend.map(({ caption }) => caption)).not.toContain('depth');
    await page.evaluate(() => { window.__GBDRAW_APP__.form.show_depth = true; });
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
    await page.evaluate(() => { window.__GBDRAW_APP__.form.show_depth = false; });
    await settleLive(page);
    await generate(page);
    const hidden = (await semanticSnapshot(page)).legend.map(({ caption }) => caption);
    expect(hidden).not.toContain('Coverage');
    expect(hidden).not.toContain('depth');
    await page.evaluate(() => { window.__GBDRAW_APP__.form.show_depth = true; });
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
        stored: state.drawings[own].legendColorOverrides[caption] ?? null,
        other: state.drawings[shown].legendColorOverrides[caption] ?? null
      };
    }, { caption: row, own: from, shown: other }), `the color stays stored in ${from} only`).toEqual({ stored: '#7b2cbf', other: null });
    await switchMode(page, from);
    await generate(page);
    const drawn = (await semanticSnapshot(page)).legend.find(({ caption }) => caption === row);
    expect(drawn?.fill.toLowerCase(), `${from} draws ${row} in the color`).toBe('#7b2cbf');
  });
}

// OV-62 (PD-OI-061 amended): two Legend rows merge only when both draw features
// of one same type. Anything else offers Suffix and Cancel only.
const legendCaptions = (page) => page.evaluate(() => window.__GBDRAW_APP__.legendEntries.map((entry) => entry.caption));
const choose = async (page, choice) => {
  await evaluateWithRetainedPromise(page, async (picked) => { await window.__GBDRAW_APP__.handleLegendRenameChoice(picked); }, choice);
  await settleLive(page);
};
const mergeButton = (page) => page.getByRole('button', { name: /Merge into existing/ });
const suffixButton = (page) => page.getByRole('button', { name: /add a suffix/ });

test('a Legend rename of GC content onto CDS offers no Merge, and Suffix equals Generate', async ({ page }) => {
  test.setTimeout(120_000);
  await open(page, { mode: 'circular', results: 'single', reflow: 'on' });
  await addColorRule(page, { feat: 'CDS', qual: 'locus_tag', val: '^FL1$', color: '#e63946', cap: 'CDS' });
  const before = await legendCaptions(page);
  await renameRow(page, 'GC content', 'CDS');
  expect(await page.evaluate(() => window.__GBDRAW_APP__.legendRenameDialog.show), 'the rename asks').toBe(true);
  await expect(mergeButton(page), 'Merge is not offered').toHaveCount(0);
  await expect(suffixButton(page), 'Suffix is offered').toHaveCount(1);
  await choose(page, 'merge');
  expect(await page.evaluate(() => window.__GBDRAW_APP__.legendRenameDialog.show), 'a refused Merge closes the dialog').toBe(false);
  expect(await legendCaptions(page), 'a refused Merge changes no row').toEqual(before);
  await renameRow(page, 'GC content', 'CDS');
  await choose(page, 'suffix');
  const renamed = await legendCaptions(page);
  expect(renamed.filter((caption) => caption === 'CDS'), 'one row keeps the caption CDS').toHaveLength(1);
  expect(renamed, 'GC content is gone').not.toContain('GC content');
  expect(renamed.some((caption) => caption !== 'CDS' && caption.startsWith('CDS')), 'the plot has a suffixed caption').toBe(true);
  await expectLiveEqualsGenerate(page, { label: 'GC content renamed onto CDS with Suffix' });
});

const RULE_ROWS = [
  { feat: 'CDS', qual: 'locus_tag', val: '^FL1$', color: '#e63946', cap: 'alpha' },
  { feat: 'CDS', qual: 'locus_tag', val: '^FL2$', color: '#2a9d8f', cap: 'beta' },
  { feat: 'repeat_region', qual: 'note', val: '^RPT_ONE$', color: '#7b2cbf', cap: 'rep' }
];

test('a Legend rename between rows of different feature types offers no Merge', async ({ page }) => {
  test.setTimeout(120_000);
  await open(page, { mode: 'circular', results: 'single', reflow: 'on' });
  for (const rule of RULE_ROWS) await addColorRule(page, rule);
  await renameRow(page, 'rep', 'alpha');
  expect(await page.evaluate(() => window.__GBDRAW_APP__.legendRenameDialog.show), 'the rename asks').toBe(true);
  await expect(mergeButton(page), 'Merge is not offered').toHaveCount(0);
  await choose(page, 'suffix');
  const renamed = await legendCaptions(page);
  expect(renamed.filter((caption) => caption === 'alpha'), 'one row keeps the caption').toHaveLength(1);
  expect(renamed.some((caption) => caption !== 'alpha' && caption.startsWith('alpha')), 'the repeat row has a suffixed caption').toBe(true);
  await expectLiveEqualsGenerate(page, { label: 'repeat_region row renamed onto a CDS row with Suffix' });
});

// Two color-rule rows of one feature type: the target is owned by a rule, so the
// rename keeps its PD-OI-042 disambiguation instead of asking, and Generate agrees.
test('a Legend rename of a rule row onto another rule row of one feature type equals Generate', async ({ page }) => {
  test.setTimeout(120_000);
  await open(page, { mode: 'circular', results: 'single', reflow: 'on' });
  for (const rule of RULE_ROWS) await addColorRule(page, rule);
  await renameRow(page, 'beta', 'alpha');
  expect(await page.evaluate(() => window.__GBDRAW_APP__.legendRenameDialog.show), 'a rule-owned caption does not ask').toBe(false);
  await settleLive(page);
  await expectLiveEqualsGenerate(page, { label: 'same-type rule row renamed onto a rule row' });
});

// U3a (R14-8): Legend edits the Result executor shows live, without the
// automatic rerender, combined with the edits of other domains. Each case
// ends with live = Generate, in both modes, with Auto Reflow on and off.
const ALPHA_ROW = { feat: 'CDS', qual: 'locus_tag', val: '^FL1$', color: '#e63946', cap: 'alpha' };
const renameAndSettle = async (page, caption, name) => { await renameRow(page, caption, name); await settleLive(page); };
const ruleIndex = (page, caption) => page.evaluate(
  (cap) => window.__GBDRAW_APP__.manualSpecificRules.findIndex((rule) => rule.cap === cap), caption
);
const U3A_PARITY = [
  {
    // Linear draws no GC rows here: its featureless row is an added one.
    name: '(1) rename a row without features, then color it',
    run: async (page, mode) => {
      const [from, to] = mode === 'circular' ? ['GC content', 'GC %'] : ['Note', 'Note row'];
      if (mode !== 'circular') await legendRowAdd(page, from, '#123456');
      await renameAndSettle(page, from, to);
      await colorLegendRow(page, to, '#264653');
    }
  },
  {
    name: '(2) rename a palette row with features, then stroke it',
    run: async (page) => { await renameAndSettle(page, 'repeat_region', 'Repeats'); await legendRowStrokeColor(page, 'Repeats', '#e63946'); }
  },
  {
    name: '(3) add a row and color it, add a row and change the palette',
    run: async (page) => {
      await legendRowAdd(page, 'Added', '#123456');
      await colorLegendRow(page, 'Added', '#abcdef');
      await legendRowAdd(page, 'Second', '#654321');
      await switchPalette(page, 'arctic');
    }
  },
  {
    // The row returns with the current palette fill and Python's stroke.
    name: '(4) delete a row, change the palette, Undo the delete',
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
    // rerender; Redo restores the Result the commit left before it, whose
    // new row is appended, so Redo asks for it too.
    name: '(10) Undo and Redo of a Rules-panel caption change of a whole row',
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
  for (const reflow of ['off', 'on']) {
    for (const { name, setup = null, run } of U3A_PARITY) {
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
    for (const reflow of ['off', 'on']) {
      // A batch Result that never drew the delete is not laid out again; the
      // Restore on it and the display of Result 1 show what Generate draws.
      test(`U3a Legend parity (8) delete on Result 1, Restore on Result 2 (circular batch, Auto Reflow ${reflow})`, async ({ page }) => {
        test.setTimeout(180_000);
        await open(page, { mode, results: 'batch', reflow });
        await deleteLegendRow(page, 'tRNA');
        await showResult(page, 1);
        await evaluateWithRetainedPromise(page, async () => { await window.__GBDRAW_APP__.restoreAllDeletedLegendEntries(); });
        await settleLive(page);
        await showResult(page, 0);
        await expectLiveEqualsGenerate(page, { label: 'Result 1 after the Restore on Result 2' });
      });
    }
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

// U3a compat: a Session 46 saved before U3a, whose Result has a renamed, a
// deleted, and an added Legend row
// (tests/fixtures/sessions/forced-label-underlay-legend-rows.provenance.json).
// Load shows the Legend it saved, Undo of a new edit returns that Result, the
// Restore of the deleted row asks for the automatic rerender (O-2), and live
// equals Generate.
const LEGEND_ROWS_SESSION = 'tests/fixtures/sessions/forced-label-underlay-legend-rows.v46.gbdraw-session.json.gz';
const LEGEND_ROWS_SAVED = JSON.parse(readFileSync('tests/fixtures/sessions/forced-label-underlay-legend-rows.provenance.json', 'utf8'))
  .sessions['forced-label-underlay-legend-rows.v46.gbdraw-session.json.gz'].legend;
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
test('U3a compat: a Session 46 saved before U3a loads its Legend rows, and live equals Generate', async ({ page }) => {
  test.setTimeout(240_000);
  await loadLegendRowsSession(page);
  await expectLiveEqualsGenerate(page, { label: 'Session 46 with Legend row edits, after Load' });
});
// U3a review H1: the saved Result renamed GC content to "GC %" without a
// record of Python's key; Load records it, so each live edit of the row
// reaches it. A Generate after each step would replace the loaded Result, so
// the color, stroke and rename are checked on the live row and then together
// against Generate; the delete starts from a fresh Load.
test('U3a compat: the row a Session 46 saved before U3a renamed takes color, stroke, rename and delete live', async ({ page }) => {
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
test('U3a compat: Undo of a new edit and the Restore of a deleted row on a Session 46 saved before U3a', async ({ page }) => {
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

// OV-239 (Owner decision R15-1): Add legend item of the caption of a deleted
// row asks before anything changes. Restore is that row's Restore (it returns
// with its own color, not the entered one); Add as gives the new row the
// suffixed caption; Cancel changes nothing. Each choice is one History step,
// Cancel none, and live equals Generate.
const addConflictDialog = (page) => page.locator('div.fixed', {
  has: page.getByRole('heading', { name: 'A deleted Legend item named “repeat_region” exists' })
});
const undoCount = (page) => page.evaluate(() => window.__GBDRAW_HISTORY__.getUndoCount());
const requestDeletedCaptionAdd = async (page) => {
  await deleteLegendRow(page, 'repeat_region');
  const before = { steps: await undoCount(page), legend: await semanticSnapshot(page) };
  await legendRowAdd(page, 'repeat_region', '#123456');
  await expect(addConflictDialog(page)).toBeVisible();
  await expect(addConflictDialog(page).getByRole('button')).toHaveText(['Restore', 'Add as “repeat_region (1)”', 'Cancel']);
  expect(await undoCount(page), 'the request records no step').toBe(before.steps);
  expect(diffSemanticSnapshots(before.legend, await semanticSnapshot(page)), 'the request changes nothing').toEqual([]);
  return before;
};
const deletedCaptions = (page) => page.evaluate(() => window.__GBDRAW_APP__.deletedLegendEntries.map((entry) => entry.caption));
const ADD_CONFLICT_CHOICES = [
  {
    choice: 'Restore',
    check: async (page) => {
      expect(await deletedCaptions(page)).toEqual([]);
      expect(await legendCaptions(page)).toContain('repeat_region');
      expect(await page.evaluate(() => window.__GBDRAW_APP__.legendEntries.find((entry) => entry.caption === 'repeat_region').color),
        'the row returns with its own color').not.toBe('#123456');
    }
  },
  {
    choice: 'Add as “repeat_region (1)”',
    check: async (page) => {
      expect(await deletedCaptions(page)).toEqual(['repeat_region']);
      expect(await page.evaluate(() => window.__GBDRAW_APP__.legendEntries.find((entry) => entry.caption === 'repeat_region (1)')?.color))
        .toBe('#123456');
    }
  },
  { choice: 'Cancel', check: async (page) => { expect(await deletedCaptions(page)).toEqual(['repeat_region']); } }
];
for (const { choice, check } of ADD_CONFLICT_CHOICES) {
  test(`OV-239: Add legend item of a deleted row's caption, ${choice} (circular)`, async ({ page }) => {
    test.setTimeout(180_000);
    await open(page, { mode: 'circular', results: 'single', reflow: 'off' });
    const before = await requestDeletedCaptionAdd(page);
    await addConflictDialog(page).getByRole('button', { name: choice, exact: true }).click();
    await expect(addConflictDialog(page)).toHaveCount(0);
    await settleLive(page);
    await check(page);
    if (choice === 'Cancel') {
      expect(await undoCount(page), 'Cancel records no step').toBe(before.steps);
      expect(diffSemanticSnapshots(before.legend, await semanticSnapshot(page)), 'Cancel changes nothing').toEqual([]);
    } else {
      expect(await undoCount(page), 'the choice is one History step').toBe(before.steps + 1);
    }
    await expectLiveEqualsGenerate(page, { label: `Add of a deleted row's caption, ${choice}` });
  });
}
