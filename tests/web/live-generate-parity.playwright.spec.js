// PD-OI-066 scenario revision 2 (WAIT-FOR-G1-APPLIES; R1, R3): each live edit
// shows on the displayed Result what the next Generate from the same draft
// draws when no pending Applies on Generate value changes those elements;
// otherwise what a Generate of the committed Session plus the applied live
// values draws. These cases change no Applies on Generate value. The
// matrix pairs every edit kind with states of the editor; the check is
// expectLiveEqualsGenerate (tests/web/helpers/live-generate-parity.cjs), minus
// the differences tests/web/contracts/live-generate-parity-allowed.json allows.
// A case marked test.fail names the finding (OV-xx) of a confirmed mismatch;
// the PR that fixes it removes the mark. The Legend cases are in
// live-generate-parity-legend.playwright.spec.js, the track-data and mode cases
// in live-generate-parity-track-data.playwright.spec.js, and the steps they all
// use in tests/web/helpers/live-generate-parity-steps.cjs.
const { expect } = require('@playwright/test');
const { readFileSync } = require('node:fs');
const { answerSessionTitle, evaluateWithRetainedPromise, expectNoLoadRepairDialog, WARM_APP_PAGE, warmAppTest: test } = require('./helpers/app-lifecycle.cjs');
const { loadSessionFile, openFresh } = require('./helpers/audit-browser.cjs');
const { expectLiveEqualsGenerate, semanticSnapshot, settleLive, showResult } = require('./helpers/live-generate-parity.cjs');
const { download, load, loadEditorLegendRows } = require('./helpers/mode-transition.cjs');
const {
  SINGLE_FIXTURE, open, openBothDrawings, otherMode, switchDrawing, generate, popupEdit, addVisibilityRule, appAction, addColorRule, history, FL1_OFF,
  legendRowColor, FL1_ALPHA, renameRow, legendRowStrokeColor, editorLegendRows, switchPalette, deleteLegendRow,
  BATCH_0004_OFF, renameLegendRow, resetSettings, displayWhileSaveHeld, saveWhileDisplaying, loadFeatureEdits
} = require('./helpers/live-generate-parity-steps.cjs');
const { WORK_ALLOWLIST } = require('./helpers/live-work-allowlist.mjs');

test.describe.configure({ retries: 0 });

// Show Labels None (Circular Label Mode None), then Generate.
const showNoLabels = async (page) => {
  await page.evaluate(() => {
    const app = window.__GBDRAW_APP__;
    if (app.mode === 'linear') app.form.show_labels_linear = 'none';
    else app.form.labels_mode = 'none';
  });
  await generate(page);
};

const BATCH_CDS_OFF = { recordId: '*', featureType: 'CDS', qualifier: 'locus_tag', value: '.', action: 'off' };

// A stroke on a Legend row (the Legend editor's stroke controls).
const legendRowStroke = async (page, caption, color, width) => {
  await page.evaluate(({ row, value, size }) => {
    const app = window.__GBDRAW_APP__;
    const index = app.legendEntries.findIndex((entry) => entry.caption === row);
    return app.updateLegendEntryStrokeColor(index, value) && app.updateLegendEntryStrokeWidth(index, size);
  }, { row: caption, value: color, size: width });
  await settleLive(page);
};

// OV-484: custom track slots with two more skew slots, a GC skew slot without
// colors of its own and an AT skew slot with its own colors, then Generate.
// The second GC skew slot has the first slot's label and swatches, so its rows
// share the first slot's rows (Owner receipt legend.row-name-multiplicity:
// one label and one swatch draw one row; tests/fixtures/legend_track_slot_rows.json, "two content
// slots and two skew slots of one label and one swatch share their rows").
const addSkewSlots = async (page) => {
  await page.evaluate(() => {
    const app = window.__GBDRAW_APP__;
    const linear = app.mode === 'linear';
    if (linear) app.form.show_skew = true;
    const slots = () => (linear ? app.adv.linear_track_slots : app.adv.circular_track_slots);
    const add = () => {
      const before = new Set(slots().map((slot) => slot.id));
      if (linear) app.addLinearTrackSlot('dinucleotide_skew'); else app.addCircularTrackSlot('dinucleotide_skew');
      return slots().find((slot) => !before.has(slot.id));
    };
    if (linear) app.setLinearTrackSlotsEnabled(true); else app.setCircularTrackSlotsEnabled(true);
    if (!slots().some((slot) => slot.renderer === 'dinucleotide_skew')) add();
    Object.assign(add().params, { nt: 'GC' });
    Object.assign(add().params, { nt: 'AT', positive_color: '#deaf6e', negative_color: '#7294e3' });
  });
  await generate(page);
  const rows = await page.evaluate(() => window.__GBDRAW_APP__.legendEntries.map((entry) => [entry.caption, entry.color]));
  expect(rows).toEqual(expect.arrayContaining([['AT skew (+)', '#deaf6e'], ['AT skew (-)', '#7294e3']]));
  const captions = rows.map(([caption]) => caption);
  expect(captions.filter((caption) => caption.startsWith('GC skew'))).toEqual(['GC skew (+)', 'GC skew (-)']);
};

// The edit kinds. Each appears at least once in the matrix.
const KINDS = [
  'visibility rule add', 'visibility rule action', 'visibility rule delete', 'feature Off', 'feature On',
  'label text', 'label Off', 'label On', 'color rule add', 'color rule color', 'color rule delete', 'feature color',
  'undo', 'redo', 'Result switch', 'legend color', 'legend stroke', 'feature stroke',
  'feature legend name', 'feature color reset', 'label Default', 'palette'
];

// The matrix: one edit kind in one set of states, with the setup before the
// edit. The states cover every pair of values (checked below); a Linear
// diagram has one Result.
const CASES = [
  {
    kind: 'visibility rule add',
    edit: 'Feature Visibility rule add',
    states: { mode: 'circular', results: 'single', reflow: 'off', labels: 'unbound' },
    run: (page) => addVisibilityRule(page, FL1_OFF)
  },
  {
    kind: 'visibility rule action',
    edit: 'Feature Visibility rule action change',
    states: { mode: 'circular', results: 'batch', reflow: 'on', labels: 'unbound' },
    setup: async (page) => { await addVisibilityRule(page, BATCH_0004_OFF); await generate(page); },
    run: (page) => appAction(page, 'setFeatureVisibilityRuleField', 0, 'action', 'show')
  },
  {
    kind: 'visibility rule delete',
    edit: 'Feature Visibility rule delete',
    states: { mode: 'linear', results: 'single', reflow: 'off', labels: 'unbound' },
    setup: async (page) => { await addVisibilityRule(page, FL1_OFF); await generate(page); },
    run: (page) => appAction(page, 'removeFeatureVisibilityRule', 0)
  },
  {
    kind: 'feature Off',
    edit: 'Feature visibility Off (popup)',
    states: { mode: 'circular', results: 'batch', reflow: 'off', labels: 'bound' },
    run: (page) => popupEdit(page, 'TESTA_0004', { visibility: 'off' })
  },
  // OV-1485: of two CDS at the same location, the one left drawn keeps its
  // rendered ID (its instance suffix) in the Result Generate draws.
  {
    kind: 'feature Off',
    edit: 'Feature visibility Off (popup) of one of two identical features',
    states: { mode: 'linear', results: 'single', reflow: 'off', labels: 'unbound' },
    run: (page) => popupEdit(page, { type: 'CDS', product: 'dup alpha' }, { visibility: 'off' })
  },
  {
    kind: 'feature On',
    edit: 'Feature visibility On (popup) for a feature a rule hid at Generate',
    states: { mode: 'circular', results: 'single', reflow: 'on', labels: 'bound' },
    setup: async (page) => { await addVisibilityRule(page, FL1_OFF); await generate(page); },
    run: (page) => popupEdit(page, 'FL1', { visibility: 'on' })
  },
  {
    kind: 'label text',
    edit: 'Label text',
    states: { mode: 'linear', results: 'single', reflow: 'on', labels: 'bound' },
    run: (page) => popupEdit(page, 'FL1', { labelText: 'alpha protein, edited live' })
  },
  {
    kind: 'label Off',
    edit: 'Label visibility Off for a label with leaders',
    states: { mode: 'circular', results: 'batch', reflow: 'off', labels: 'bound' },
    run: (page) => popupEdit(page, { type: 'tRNA', product: 'tRNA-Leu' }, { labelVisibility: 'off' })
  },
  {
    kind: 'label text',
    edit: 'Label text',
    states: { mode: 'circular', results: 'single', reflow: 'off', labels: 'bound' },
    // Too long to stay inside FL1's arc: Generate draws it outside with leaders.
    run: (page) => popupEdit(page, 'FL1', { labelText: 'alpha protein, edited live with a text far too long to stay inside its own arc' })
  },
  {
    kind: 'label On',
    edit: 'Label visibility On for a label Off at Generate',
    states: { mode: 'linear', results: 'single', reflow: 'off', labels: 'bound' },
    setup: async (page) => { await popupEdit(page, 'FL2', { labelVisibility: 'off' }); await generate(page); },
    run: (page) => popupEdit(page, 'FL2', { labelVisibility: 'on' })
  },
  {
    kind: 'color rule add',
    edit: 'Specific color rule add',
    states: { mode: 'circular', results: 'batch', reflow: 'on', labels: 'bound' },
    run: (page) => addColorRule(page, { feat: 'CDS', qual: 'product', val: '^gtg start$', color: '#e63946', cap: 'GTG start' })
  },
  {
    kind: 'color rule color',
    edit: 'Specific color rule color change',
    states: { mode: 'linear', results: 'single', reflow: 'off', labels: 'unbound' },
    setup: async (page) => {
      await addColorRule(page, { feat: 'CDS', qual: 'locus_tag', val: '^FL2$', color: '#2a9d8f', cap: 'beta' });
      await generate(page);
    },
    run: (page) => appAction(page, 'setSpecificRuleField', 0, 'color', '#7b2cbf')
  },
  {
    kind: 'feature color',
    edit: 'Feature color (popup, this feature only)',
    states: { mode: 'circular', results: 'single', reflow: 'on', labels: 'bound' },
    run: (page) => popupEdit(page, 'FL1', { fill: '#c83366' })
  },
  {
    kind: 'feature color',
    edit: 'Feature color change of a feature colored at the last Generate',
    states: { mode: 'circular', results: 'batch', reflow: 'off', labels: 'unbound' },
    setup: async (page) => { await popupEdit(page, 'TESTA_0005', { fill: '#c83366' }); await generate(page); },
    run: (page) => popupEdit(page, 'TESTA_0005', { fill: '#123456' })
  },
  // OV-344: a This feature only color on a hidden feature of the second
  // Result, in the caption and color of a rule only the first Result draws.
  // Generate draws no row for it on the second Result.
  {
    kind: 'feature color',
    edit: 'Feature color (popup, this feature only) of a hidden feature, like a rule only another Result draws',
    states: { mode: 'circular', results: 'batch', reflow: 'off', labels: 'unbound' },
    setup: async (page) => {
      await popupEdit(page, 'TESTA_0002', { fill: '#7b2cbf' });
      await generate(page);
      await showResult(page, 1);
      await popupEdit(page, 'TESTB_0001', { visibility: 'off' });
    },
    run: (page) => popupEdit(page, 'TESTB_0001', { fill: '#7b2cbf' })
  },
  {
    kind: 'undo',
    edit: 'Undo of a Feature Visibility rule add',
    states: { mode: 'circular', results: 'batch', reflow: 'off', labels: 'unbound' },
    setup: (page) => addVisibilityRule(page, BATCH_0004_OFF),
    run: (page) => history(page, 'undo')
  },
  {
    kind: 'redo',
    edit: 'Redo of a specific color rule color change',
    states: { mode: 'linear', results: 'single', reflow: 'on', labels: 'unbound' },
    setup: async (page) => {
      await addColorRule(page, { feat: 'CDS', qual: 'locus_tag', val: '^FL1$', color: '#e63946', cap: 'alpha' });
      await generate(page);
      await appAction(page, 'setSpecificRuleField', 0, 'color', '#7b2cbf');
      await history(page, 'undo');
    },
    run: (page) => history(page, 'redo')
  },
  {
    kind: 'Result switch',
    edit: 'Result switch after a feature color with the same-product scope',
    states: { mode: 'circular', results: 'batch', reflow: 'on', labels: 'bound' },
    setup: (page) => popupEdit(page, 'TESTA_0001', { fill: '#c83366', scope: 'annotationLabel' }),
    run: (page) => showResult(page, 1)
  },
  {
    kind: 'Result switch',
    edit: 'Result switch after a Feature Visibility rule add',
    states: { mode: 'circular', results: 'batch', reflow: 'off', labels: 'unbound' },
    setup: (page) => addVisibilityRule(page, BATCH_0004_OFF),
    run: (page) => showResult(page, 1)
  },
  // OV-42, OV-43, OV-44 (Owner decision 2026-10-06, option A): an edit that
  // changes what a Legend derives from asks for the automatic rerender, also
  // with Auto Reflow off, and Python redraws the Legend of every Result.
  {
    kind: 'visibility rule add',
    edit: 'Feature Visibility rule add that hides every CDS',
    states: { mode: 'circular', results: 'batch', reflow: 'off', labels: 'unbound' },
    run: (page) => addVisibilityRule(page, BATCH_CDS_OFF)
  },
  {
    kind: 'feature Off',
    edit: 'Feature visibility Off (popup) for the first drawn CDS',
    states: { mode: 'linear', results: 'single', reflow: 'off', labels: 'bound' },
    run: (page) => popupEdit(page, 'FL1', { visibility: 'off' })
  },
  {
    kind: 'undo',
    edit: 'Undo of a Feature Visibility rule add that reorders the Legend',
    states: { mode: 'circular', results: 'single', reflow: 'off', labels: 'unbound' },
    setup: (page) => addVisibilityRule(page, FL1_OFF),
    run: (page) => history(page, 'undo')
  },
  {
    kind: 'redo',
    edit: 'Redo of a Feature Visibility rule add that reorders the Legend',
    states: { mode: 'circular', results: 'single', reflow: 'off', labels: 'unbound' },
    setup: async (page) => { await addVisibilityRule(page, FL1_OFF); await history(page, 'undo'); },
    run: (page) => history(page, 'redo')
  },
  {
    kind: 'color rule delete',
    edit: 'Specific color rule delete of a rule Generate drew',
    states: { mode: 'circular', results: 'single', reflow: 'off', labels: 'unbound' },
    setup: async (page) => { await addColorRule(page, FL1_ALPHA); await generate(page); },
    run: (page) => appAction(page, 'removeSpecificRule', 0)
  },
  {
    kind: 'undo',
    edit: 'Undo of a specific color rule add',
    states: { mode: 'circular', results: 'single', reflow: 'off', labels: 'unbound' },
    setup: (page) => addColorRule(page, FL1_ALPHA),
    run: (page) => history(page, 'undo')
  },
  {
    kind: 'redo',
    edit: 'Redo of a specific color rule add',
    states: { mode: 'circular', results: 'single', reflow: 'off', labels: 'unbound' },
    setup: async (page) => { await addColorRule(page, FL1_ALPHA); await history(page, 'undo'); },
    run: (page) => history(page, 'redo')
  },
  // OV-61: a rule row recolored after a popup edit colored its only feature.
  // The rule is the latest explicit color; the popup's Legend color does not
  // come back at Generate.
  {
    kind: 'legend color',
    edit: 'Legend row color of a feature colored in the popup',
    states: { mode: 'linear', results: 'single', reflow: 'off', labels: 'unbound' },
    setup: (page) => popupEdit(page, { type: 'repeat_region' }, { fill: '#e63946' }),
    run: (page) => legendRowColor(page, 'repeat_region', '#f4a261')
  },
  // OV-123: a stroke on a generated Legend row reaches the features drawn in
  // the row's color, live and at Generate; a batch Result shows it on display.
  {
    kind: 'legend stroke',
    edit: 'Legend row stroke on a generated row',
    states: { mode: 'linear', results: 'single', reflow: 'off', labels: 'unbound' },
    run: (page) => legendRowStroke(page, 'CDS', '#e63946', 3)
  },
  {
    kind: 'legend stroke',
    edit: 'Legend row stroke on a generated row',
    states: { mode: 'circular', results: 'single', reflow: 'on', labels: 'bound' },
    run: (page) => legendRowStroke(page, 'CDS', '#e63946', 3)
  },
  {
    kind: 'undo',
    edit: 'Undo of a Legend row stroke on a generated row',
    states: { mode: 'linear', results: 'single', reflow: 'on', labels: 'bound' },
    setup: (page) => legendRowStrokeColor(page, 'CDS', '#e63946'),
    run: (page) => history(page, 'undo')
  },
  {
    kind: 'redo',
    edit: 'Redo of a Legend row stroke on a generated row',
    states: { mode: 'circular', results: 'single', reflow: 'off', labels: 'unbound' },
    setup: async (page) => { await legendRowStrokeColor(page, 'CDS', '#e63946'); await history(page, 'undo'); },
    run: (page) => history(page, 'redo')
  },
  {
    kind: 'legend stroke',
    edit: 'Legend row stroke Reset of a stroke Generate drew',
    states: { mode: 'linear', results: 'single', reflow: 'on', labels: 'unbound' },
    setup: async (page) => { await legendRowStroke(page, 'CDS', '#e63946', 3); await generate(page); },
    run: (page) => appAction(page, 'resetLegendEntryStroke', 0)
  },
  // A feature's own stroke (the popup) wins over its Legend row's stroke, live
  // and at Generate; without it, the feature shows the row's stroke.
  {
    kind: 'legend stroke',
    edit: 'Legend row stroke on a row with a feature stroked in the popup',
    states: { mode: 'circular', results: 'single', reflow: 'off', labels: 'bound' },
    setup: (page) => popupEdit(page, 'FL1', { stroke: '#2a9d8f' }),
    run: (page) => legendRowStroke(page, 'CDS', '#e63946', 3)
  },
  {
    kind: 'feature stroke',
    edit: 'Feature stroke (popup, this feature only) on a feature of a stroked Legend row',
    states: { mode: 'linear', results: 'single', reflow: 'off', labels: 'bound' },
    setup: (page) => legendRowStroke(page, 'CDS', '#e63946', 3),
    run: (page) => popupEdit(page, 'FL1', { stroke: '#2a9d8f' })
  },
  {
    kind: 'feature stroke',
    edit: 'Feature stroke Reset (popup) of a feature in a stroked Legend row',
    states: { mode: 'circular', results: 'single', reflow: 'on', labels: 'bound' },
    setup: async (page) => {
      await legendRowStroke(page, 'CDS', '#e63946', 3);
      await popupEdit(page, 'FL1', { stroke: '#2a9d8f' });
      await generate(page);
    },
    run: (page) => popupEdit(page, 'FL1', { resetStroke: true })
  },
  {
    kind: 'Result switch',
    edit: 'Result switch after a Legend row stroke on a generated row',
    states: { mode: 'circular', results: 'batch', reflow: 'off', labels: 'unbound' },
    setup: (page) => legendRowStroke(page, 'CDS', '#e63946', 3),
    run: (page) => showResult(page, 1)
  },
  // OV-144 (D-07, PD-OI-062): a batch Result shown with a stroke edit that was
  // removed while another Result was shown returns to Python's stroke when it
  // is shown again.
  {
    kind: 'Result switch',
    edit: 'Result switch back to a Result shown with a Legend row stroke that was reset since',
    states: { mode: 'circular', results: 'batch', reflow: 'on', labels: 'unbound' },
    setup: async (page) => {
      await legendRowStroke(page, 'CDS', '#e63946', 3);
      await showResult(page, 1);
      await showResult(page, 0);
      await page.evaluate((row) => {
        const app = window.__GBDRAW_APP__;
        return app.resetLegendEntryStroke(app.legendEntries.findIndex((entry) => entry.caption === row));
      }, 'CDS');
      await settleLive(page);
    },
    run: (page) => showResult(page, 1)
  },
  {
    kind: 'color rule color',
    edit: 'Specific color rule color change of a feature colored in the popup',
    states: { mode: 'circular', results: 'single', reflow: 'on', labels: 'bound' },
    setup: (page) => popupEdit(page, { type: 'repeat_region' }, { fill: '#e63946' }),
    run: (page) => appAction(page, 'setSpecificRuleField', 0, 'color', '#f4a261')
  },
  {
    kind: 'Result switch',
    edit: 'Result switch after a feature color with the same-product scope',
    states: { mode: 'circular', results: 'batch', reflow: 'off', labels: 'bound' },
    setup: (page) => popupEdit(page, 'TESTA_0001', { fill: '#c83366', scope: 'annotationLabel' }),
    run: (page) => showResult(page, 1)
  },
  {
    kind: 'Result switch',
    edit: 'Result switch after a color change of a same-product rule Generate drew',
    states: { mode: 'circular', results: 'batch', reflow: 'off', labels: 'unbound' },
    setup: async (page) => {
      await popupEdit(page, 'TESTA_0001', { fill: '#c83366', scope: 'annotationLabel' });
      await generate(page);
      await appAction(page, 'setSpecificRuleField', 0, 'color', '#2a9d8f');
    },
    run: (page) => showResult(page, 1)
  },
  // OV-63: a Legend color on a row whose features are then all hidden. Generate
  // draws no row for them, and the stored color does not fail the Generate.
  {
    kind: 'feature Off',
    edit: 'Feature visibility Off (popup) for every feature of a Legend row with a Legend color',
    states: { mode: 'circular', results: 'single', reflow: 'off', labels: 'unbound' },
    setup: async (page) => {
      await page.evaluate(() => {
        const app = window.__GBDRAW_APP__;
        return app.updateLegendEntryColor(app.legendEntries.findIndex((entry) => entry.caption === 'repeat_region'), '#7b2cbf');
      });
      await settleLive(page);
    },
    run: (page) => popupEdit(page, { type: 'repeat_region' }, { visibility: 'off' })
  },
  {
    kind: 'feature Off',
    edit: 'Feature visibility Off (popup) for every feature of a Legend row renamed in the Legend',
    states: { mode: 'circular', results: 'single', reflow: 'off', labels: 'unbound' },
    setup: async (page) => {
      await page.evaluate(() => {
        const app = window.__GBDRAW_APP__;
        return app.renameLegendEntry(app.legendEntries.findIndex((entry) => entry.caption === 'repeat_region'), 'Repeats');
      });
      await settleLive(page);
    },
    run: (page) => popupEdit(page, { type: 'repeat_region' }, { visibility: 'off' })
  },
  {
    kind: 'visibility rule add',
    edit: 'Feature Visibility rule add that hides every feature of a Legend row with a Legend color',
    states: { mode: 'circular', results: 'single', reflow: 'off', labels: 'bound' },
    setup: async (page) => {
      await page.evaluate(() => {
        const app = window.__GBDRAW_APP__;
        return app.updateLegendEntryColor(app.legendEntries.findIndex((entry) => entry.caption === 'CDS'), '#7b2cbf');
      });
      await settleLive(page);
    },
    run: (page) => addVisibilityRule(page, BATCH_CDS_OFF)
  },
  {
    kind: 'color rule add',
    edit: 'Specific color rule add that recaptions every feature of a Legend row with a Legend color',
    states: { mode: 'circular', results: 'single', reflow: 'on', labels: 'unbound' },
    setup: async (page) => {
      await page.evaluate(() => {
        const app = window.__GBDRAW_APP__;
        return app.updateLegendEntryColor(app.legendEntries.findIndex((entry) => entry.caption === 'CDS'), '#7b2cbf');
      });
      await settleLive(page);
    },
    run: (page) => addColorRule(page, { feat: 'CDS', qual: 'locus_tag', val: '.', color: '#c83366', cap: 'Zeta' })
  },
  // UJ-03 (GUI journeys 2026-10-06, fixed by #857): a popup Legend name for
  // one feature of a shared row, and Reset fill color after a This feature
  // only color, change the rows Python derives (`other proteins`); the
  // automatic rerender shows them live.
  {
    kind: 'feature legend name',
    edit: 'Legend name (popup, this feature only) for one feature of a shared row',
    states: { mode: 'circular', results: 'single', reflow: 'on', labels: 'unbound' },
    run: (page) => popupEdit(page, 'FL1', { legendName: 'Complex IV' })
  },
  {
    kind: 'feature legend name',
    edit: 'Legend name (popup, this feature only) for one feature of a shared row',
    states: { mode: 'linear', results: 'single', reflow: 'off', labels: 'bound' },
    run: (page) => popupEdit(page, 'FL1', { legendName: 'Complex IV' })
  },
  {
    kind: 'feature color reset',
    edit: 'Reset fill color (popup) after a This feature only color',
    states: { mode: 'circular', results: 'single', reflow: 'off', labels: 'bound' },
    setup: (page) => popupEdit(page, 'FL1', { fill: '#c83366' }),
    run: (page) => popupEdit(page, 'FL1', { resetFill: true })
  },
  {
    kind: 'feature color reset',
    edit: 'Reset fill color (popup) of a This feature only color Generate drew',
    states: { mode: 'linear', results: 'single', reflow: 'on', labels: 'unbound' },
    setup: async (page) => { await popupEdit(page, 'FL1', { fill: '#c83366' }); await generate(page); },
    run: (page) => popupEdit(page, 'FL1', { resetFill: true })
  },
  // Review U2b #1 (cross-domain): a Legend structure edit, then a paint edit.
  // A row without features (GC content) is renamed live: the rename rewrites
  // the row's key on the Result, and a later palette or rule reconcile must
  // still find the row. (A feature row's rename is a rule and rerenders.)
  {
    kind: 'palette',
    edit: 'Palette change after a Legend row rename with a Legend color',
    states: { mode: 'circular', results: 'single', reflow: 'off', labels: 'unbound' },
    setup: async (page) => { await legendRowColor(page, 'GC content', '#ff0000'); await renameLegendRow(page, 'GC content', 'GC percent'); },
    run: (page) => switchPalette(page, 'alpine_retreat')
  },
  {
    kind: 'palette',
    edit: 'Palette change after a palette change and a Legend row rename',
    states: { mode: 'circular', results: 'single', reflow: 'on', labels: 'bound' },
    setup: async (page) => { await switchPalette(page, 'alpine_retreat'); await renameLegendRow(page, 'GC content', 'GC percent'); },
    run: (page) => switchPalette(page, 'arctic')
  },
  {
    kind: 'palette',
    edit: 'Palette change with a Legend editor row of a loaded Session',
    states: { mode: 'circular', results: 'batch', reflow: 'on', labels: 'unbound' },
    setup: (page) => editorLegendRows(page, [['Manual row', '#118833']]),
    run: (page) => switchPalette(page, 'alpine_retreat')
  },
  // OV-484: the Legend rows of a skew slot with its own colors keep them
  // through a palette change, and the rows a second slot of one label and
  // one swatch shares with the first slot follow the palette.
  ...['circular', 'linear'].map((mode) => ({
    kind: 'palette',
    edit: 'Palette change with skew slots of their own colors and of one label',
    states: { mode, results: 'single', reflow: 'off', labels: 'unbound' },
    setup: addSkewSlots,
    run: (page) => switchPalette(page, 'alpine_retreat')
  })),
  // Review U2BFIX #5: a Legend row delete (a writer not yet migrated to the
  // edit port) with a palette change, a Legend row stroke, and a Legend row
  // color shown through the port.
  {
    kind: 'palette',
    edit: 'Palette change after a Legend row delete',
    states: { mode: 'circular', results: 'batch', reflow: 'off', labels: 'unbound' },
    setup: (page) => deleteLegendRow(page, 'tRNA'),
    run: (page) => switchPalette(page, 'alpine_retreat')
  },
  {
    kind: 'legend stroke',
    edit: 'Legend row stroke after a Legend row delete',
    states: { mode: 'linear', results: 'single', reflow: 'off', labels: 'unbound' },
    setup: (page) => deleteLegendRow(page, 'repeat_region'),
    run: (page) => legendRowStroke(page, 'CDS', '#e63946', 3)
  },
  {
    kind: 'legend color',
    edit: 'Legend row color after a Legend row delete',
    states: { mode: 'circular', results: 'single', reflow: 'on', labels: 'bound' },
    setup: (page) => deleteLegendRow(page, 'repeat_region'),
    run: (page) => legendRowColor(page, 'CDS', '#f4a261')
  },
  {
    kind: 'color rule add',
    edit: 'Specific color rule add after a Legend row rename with a Legend color',
    states: { mode: 'circular', results: 'single', reflow: 'on', labels: 'bound' },
    setup: async (page) => { await legendRowColor(page, 'GC content', '#7b2cbf'); await renameLegendRow(page, 'GC content', 'GC percent'); },
    run: (page) => addColorRule(page, FL1_ALPHA)
  },
  // Review U2a #4 (OV-146) on a batch Result: Undo of a Legend color or a rule
  // commit after a palette switch shows the palette Generate draws.
  {
    kind: 'undo',
    edit: 'Undo of a Legend row color on the other batch Result after a palette switch',
    states: { mode: 'circular', results: 'batch', reflow: 'off', labels: 'unbound' },
    setup: async (page) => {
      await switchPalette(page, 'alpine_retreat');
      await legendRowColor(page, 'CDS', '#ff0000');
      await showResult(page, 1);
    },
    run: (page) => history(page, 'undo')
  },
  {
    kind: 'undo',
    edit: 'Undo of a color rule commit after a palette switch',
    states: { mode: 'circular', results: 'batch', reflow: 'on', labels: 'unbound' },
    setup: async (page) => {
      await switchPalette(page, 'alpine_retreat');
      await addColorRule(page, { feat: 'CDS', qual: 'locus_tag', val: '_0004$', color: '#c83366', cap: 'CDS' });
    },
    run: (page) => history(page, 'undo')
  },
  {
    kind: 'Result switch',
    edit: 'Result switch after deleting a same-product rule Generate drew',
    states: { mode: 'circular', results: 'batch', reflow: 'off', labels: 'unbound' },
    setup: async (page) => {
      await popupEdit(page, 'TESTA_0001', { fill: '#c83366', scope: 'annotationLabel' });
      await generate(page);
      await appAction(page, 'removeSpecificRule', 0);
    },
    run: (page) => showResult(page, 1)
  },
  // UJ-01: a label drawn only because of Label visibility On (Show Labels
  // None) leaves the Result when the intent returns to Default, by Undo or by
  // the popup; a label Off that the rerender left out comes back.
  {
    kind: 'undo',
    edit: 'Undo of a label text with Label visibility On under Show Labels None',
    states: { mode: 'circular', results: 'single', reflow: 'off', labels: 'unbound' },
    setup: async (page) => {
      await showNoLabels(page);
      await popupEdit(page, 'FL1', { labelText: 'FL1-X', labelVisibility: 'on' });
    },
    run: (page) => history(page, 'undo')
  },
  {
    kind: 'undo',
    edit: 'Undo of a label text with Label visibility On under Show Labels None',
    states: { mode: 'linear', results: 'single', reflow: 'on', labels: 'bound' },
    setup: async (page) => {
      await showNoLabels(page);
      await popupEdit(page, 'FL1', { labelText: 'FL1-X', labelVisibility: 'on' });
    },
    run: (page) => history(page, 'undo')
  },
  {
    kind: 'label Default',
    edit: 'Label visibility Default after On under Show Labels None',
    states: { mode: 'circular', results: 'batch', reflow: 'off', labels: 'bound' },
    setup: async (page) => {
      await showNoLabels(page);
      await popupEdit(page, 'TESTA_0002', { labelVisibility: 'on' });
    },
    run: (page) => popupEdit(page, 'TESTA_0002', { labelVisibility: 'default' })
  },
  {
    kind: 'label Default',
    edit: 'Label visibility Default after On under Show Labels None',
    states: { mode: 'linear', results: 'single', reflow: 'on', labels: 'unbound' },
    setup: async (page) => {
      await showNoLabels(page);
      await popupEdit(page, 'FL1', { labelVisibility: 'on' });
    },
    run: (page) => popupEdit(page, 'FL1', { labelVisibility: 'default' })
  },
  {
    kind: 'undo',
    edit: 'Undo of Label visibility Off that the rerender drew',
    states: { mode: 'circular', results: 'single', reflow: 'on', labels: 'unbound' },
    setup: (page) => popupEdit(page, 'FL2', { labelVisibility: 'off' }),
    run: (page) => history(page, 'undo')
  }
];

// The label binding state before the edit: a popup opened on another feature.
const BIND_TARGET = { single: 'FL2', batch: 'TESTA_0001' };

const statesName = ({ mode, results, reflow, labels }) => (
  `${mode}, ${results === 'batch' ? 'two-Result batch' : 'one Result'}, Auto Reflow ${reflow}, labels ${labels}`
);

test('the matrix covers every edit kind and every pair of states', () => {
  expect(KINDS.filter((kind) => !CASES.some((entry) => entry.kind === kind))).toEqual([]);
  expect(CASES.filter((entry) => !KINDS.includes(entry.kind)).map(({ edit }) => edit)).toEqual([]);
  const values = { mode: ['circular', 'linear'], results: ['single', 'batch'], reflow: ['on', 'off'], labels: ['bound', 'unbound'] };
  const factors = Object.keys(values);
  const missing = [];
  for (const [index, left] of factors.entries()) {
    for (const right of factors.slice(index + 1)) {
      for (const leftValue of values[left]) {
        for (const rightValue of values[right]) {
          if (left === 'mode' && leftValue === 'linear' && right === 'results' && rightValue === 'batch') continue;
          if (!CASES.some(({ states }) => states[left] === leftValue && states[right] === rightValue)) {
            missing.push(`${left}=${leftValue} with ${right}=${rightValue}`);
          }
        }
      }
    }
  }
  expect(missing).toEqual([]);
  expect(CASES.filter(({ states }) => states.mode === 'linear' && states.results === 'batch')).toEqual([]);
});

// An allowed difference applies to one kind and its attributes under a stated
// condition, and cites the decision or rule that allows it.
test('each allowed difference names its scope and cites its decision', () => {
  const { entries } = require('./contracts/live-generate-parity-allowed.json');
  for (const entry of entries) {
    expect(Object.keys(entry).sort(), entry.id).toEqual(['attributes', 'condition', 'decision', 'id', 'kind', 'reason']);
    expect(['feature', 'label', 'leader', 'legend', 'text'], entry.id).toContain(entry.kind);
    expect(entry.attributes.length, entry.id).toBeGreaterThan(0);
    expect(Object.keys(entry.condition).length, entry.id).toBeGreaterThan(0);
    expect(entry.reason.trim() && entry.decision.trim(), entry.id).toBeTruthy();
  }
  expect(new Set(entries.map(({ id }) => id)).size).toBe(entries.length);
});

// OV-287, OV-289, OV-290 (R-04, D-05): Reset Settings clears the editor
// edits, and the displayed Result shows the reset draft as the next Generate
// draws it, one case per edit domain. The Legend list (caption and the color
// the editor shows) equals the list after that Generate. A Reset before the
// first Generate puts the draft at the default settings, so the Reset under
// test changes no setting that only Generate draws.
const legendList = (page) => page.evaluate(() => {
  const app = window.__GBDRAW_APP__;
  return app.legendEntries.map((entry) => ({ caption: entry.caption, color: app.legendEntryColor(entry) }));
});
// A rerender held until released, so a case does not pass by timing.
const holdRerenders = (page) => page.evaluate(() => {
  window.__GBDRAW_TEST_HOOKS__ ||= {};
  window.__heldRerenders = [];
  window.__GBDRAW_TEST_HOOKS__.beforeDiagramGenerationResponse = () => new Promise((resolve) => { window.__heldRerenders.push(resolve); });
});
const releaseRerenders = async (page) => {
  await page.evaluate(() => {
    delete window.__GBDRAW_TEST_HOOKS__.beforeDiagramGenerationResponse;
    window.__heldRerenders.forEach((resolve) => resolve());
  });
  await settleLive(page);
};
const rerenderHeld = (page) => page.waitForFunction(() => window.__heldRerenders.length > 0);
const RESET_CASES = [
  {
    domain: 'a Legend row stroke (OV-287)',
    edit: (page) => legendRowStrokeColor(page, 'CDS', '#e63946'),
    // Undo shows the edit again and Redo the reset draft, on the Result.
    undoRedo: true
  },
  {
    domain: 'a Legend row color (U5 review L2)',
    edit: (page) => legendRowColor(page, 'CDS', '#7b2cbf')
  },
  {
    domain: 'a deleted Legend row (OV-289)',
    edit: (page) => deleteLegendRow(page, 'GC content')
  },
  {
    // Reset clears a row the editor added (R15-2: only a Session holds one,
    // `loadEditorLegendRows`); the deleted row stays gone.
    domain: 'a Legend row added in the editor and deleted',
    edit: async (page) => { await loadEditorLegendRows(page, [['Manual row', '#118833']]); await deleteLegendRow(page, 'Manual row'); },
    sameDrawing: true
  },
  {
    // The next two: Generate after Reset draws neither the added row nor the new caption.
    domain: 'a Legend row added in the editor',
    edit: (page) => loadEditorLegendRows(page, [['Manual row', '#118833']])
  },
  {
    // R15-14: Reset returns every Legend rename, a feature row's and a track
    // row's, to the caption Python draws; Undo of Reset restores them.
    domain: 'Legend rows renamed, a feature row and a track row (R15-14)',
    edit: async (page) => {
      await renameRow(page, 'CDS', 'Coding');
      await settleLive(page);
      await renameRow(page, 'GC content', 'GC %');
      await settleLive(page);
    },
    undoRedo: true,
    resetShowsDrawn: true
  },
  {
    // A popup "this feature only" color draws its own Legend row; Reset
    // removes its rule, and the rule owner asks for the rerender (OV-43).
    domain: 'a feature fill, this feature only, with its own Legend row (OV-290)',
    edit: (page) => popupEdit(page, 'FL1', { fill: '#2a9d8f' })
  },
  {
    domain: 'a feature hidden (OV-290)',
    edit: (page) => popupEdit(page, 'FL2', { visibility: 'off' }),
    undoRedo: true
  },
  {
    // The Result Generate drew lacks the feature, so Reset asks for the rerender.
    domain: 'a feature a visibility rule hid at Generate',
    edit: async (page) => { await addVisibilityRule(page, FL1_OFF); await generate(page); }
  },
  {
    domain: 'a label text with Label visibility On (OV-290)',
    edit: (page) => popupEdit(page, 'FL1', { labelText: 'FL1-X', labelVisibility: 'on' })
  }
];
for (const { domain, edit, undoRedo, sameDrawing, resetShowsDrawn } of RESET_CASES) {
  test(`Reset Settings after ${domain} shows what Generate draws`, async ({ page }) => {
    test.setTimeout(120_000);
    await open(page, { mode: 'circular', results: 'single', reflow: 'off' });
    await resetSettings(page);
    await generate(page);
    const drawn = await semanticSnapshot(page);
    await edit(page);
    const edited = await semanticSnapshot(page);
    if (sameDrawing) expect(edited, 'the edit leaves the drawing as Generate drew it').toEqual(drawn);
    else expect(edited, 'the edit shows on the Result').not.toEqual(drawn);
    await resetSettings(page);
    const reset = await semanticSnapshot(page);
    if (resetShowsDrawn) expect(reset, 'Reset Settings shows the drawing Generate drew').toEqual(drawn);
    if (undoRedo) {
      await history(page, 'undo');
      expect(await semanticSnapshot(page), 'Undo of Reset Settings').toEqual(edited);
      await history(page, 'redo');
      expect(await semanticSnapshot(page), 'Redo of Reset Settings').toEqual(reset);
    }
    const listed = await legendList(page);
    await expectLiveEqualsGenerate(page, { label: `Reset Settings after ${domain}` });
    expect(listed, 'the Legend list after Reset Settings').toEqual(await legendList(page));
  });
}

// OV-2293: Reset Settings while the other drawing is shown. Generate drew the
// edit into this drawing's Result, so its bytes lack what the reset draft draws:
// the arrival asks for the rerender a shown Reset asks for. The settings that
// apply on Generate are written back, so the comparison covers what the Reset
// shows live.
const LABEL_SETTINGS = Object.freeze([
  ['circular', 'form', 'labels_mode'], ['circular', 'form', 'multi_record_canvas'], ['circular', 'adv', 'circular_grouping_intent'],
  ['linear', 'form', 'show_labels_linear']
]);
const readLabelSettings = (page) => page.evaluate(async (settings) => {
  const { state } = await import('/gbdraw/web/js/state.js');
  return settings.map(([drawingMode, group, name]) => [drawingMode, group, name, window.Vue.unref(state.drawingOfMode(drawingMode)[group])[name]]);
}, LABEL_SETTINGS);
const writeLabelSettings = (page, values) => page.evaluate(async (written) => {
  const { state } = await import('/gbdraw/web/js/state.js');
  written.forEach(([drawingMode, group, name, value]) => { window.Vue.unref(state.drawingOfMode(drawingMode)[group])[name] = value; });
}, values);
// Reset writes Auto Reflow off and default label settings. A case keeps both, so
// its display runs with the state its title names.
const resetKeepingLabelSettings = async (page) => {
  const kept = await readLabelSettings(page);
  const autoReflow = await page.evaluate(() => window.__GBDRAW_APP__.autoLabelReflowEnabled);
  await resetSettings(page);
  await writeLabelSettings(page, kept);
  await setAutoReflow(page, autoReflow);
};
const resetFromOtherDrawing = async (page, mode) => {
  await switchDrawing(page, otherMode(mode));
  await resetKeepingLabelSettings(page);
  await switchDrawing(page, mode);
};
const OTHER_DRAWING_RESET_EDITS = {
  'a hidden label': (page) => popupEdit(page, 'FL2', { labelVisibility: 'off' }),
  'a hidden feature': (page) => popupEdit(page, 'FL2', { visibility: 'off' }),
  'an edited label text': (page) => popupEdit(page, 'FL1', { labelText: 'walk label' }),
  'a visibility rule': (page) => addVisibilityRule(page, { recordId: '*', featureType: 'CDS', qualifier: 'locus_tag', value: '^FL2$', action: 'off' })
};
for (const mode of ['linear', 'circular']) {
  for (const [name, edit] of Object.entries(OTHER_DRAWING_RESET_EDITS)) {
    test(`Reset Settings from the other drawing after ${name} and a Generate shows what Generate draws (${mode})`, async ({ page }) => {
      test.setTimeout(240_000);
      await openBothDrawings(page, { mode, results: 'single', reflow: 'off' });
      const drawn = await semanticSnapshot(page);
      await edit(page);
      await generate(page);
      expect(await semanticSnapshot(page), 'Generate draws the edit').not.toEqual(drawn);
      await resetFromOtherDrawing(page, mode);
      await expectLiveEqualsGenerate(page, { label: `Reset Settings from the other drawing after ${name}` });
    });
    test(`Reset Settings from the other drawing after ${name} that no Generate drew shows what Generate draws (${mode})`, async ({ page }) => {
      test.setTimeout(240_000);
      await openBothDrawings(page, { mode, results: 'single', reflow: 'off' });
      await edit(page);
      await resetFromOtherDrawing(page, mode);
      await expectLiveEqualsGenerate(page, { label: `Reset Settings from the other drawing after ${name}, not generated` });
    });
  }
}

// The shown Result itself, and another Result of the same drawing.
for (const [name, edit] of Object.entries(OTHER_DRAWING_RESET_EDITS)) {
  test(`Reset Settings on the Result that Generate drew after ${name} shows what Generate draws`, async ({ page }) => {
    test.setTimeout(240_000);
    await open(page, { mode: 'circular', results: 'single', reflow: 'off' });
    await edit(page);
    await generate(page);
    await resetKeepingLabelSettings(page);
    await expectLiveEqualsGenerate(page, { label: `Reset Settings on the Result that Generate drew after ${name}` });
  });
}
for (const [name, change] of Object.entries({
  'a hidden label': { labelVisibility: 'off' },
  'a hidden feature': { visibility: 'off' },
  'an edited label text': { labelText: 'walk label' }
})) {
  test(`Reset Settings while another batch Result is shown after ${name} and a Generate shows what Generate draws`, async ({ page }) => {
    test.setTimeout(240_000);
    await open(page, { mode: 'circular', results: 'batch', reflow: 'off' });
    await showResult(page, 1);
    await popupEdit(page, { biological_feature_id: 'f5aa08bf5', record_idx: 1 }, change);
    await generate(page);
    await showResult(page, 0);
    await resetKeepingLabelSettings(page);
    await showResult(page, 1);
    await expectLiveEqualsGenerate(page, { label: `Reset Settings on another batch Result after ${name}` });
  });
}

test('Undo of Reset Settings from the other drawing shows what Generate draws, and Redo shows the reset', async ({ page }) => {
  test.setTimeout(240_000);
  await openBothDrawings(page, { mode: 'linear', results: 'single', reflow: 'off' });
  await popupEdit(page, 'FL2', { visibility: 'off' });
  await generate(page);
  await resetFromOtherDrawing(page, 'linear');
  const reset = await semanticSnapshot(page);
  expect(Object.keys(reset.features), 'the rerender draws the feature').toContain('f841fb8a8');
  await history(page, 'undo');
  await settleLive(page);
  expect(Object.keys((await semanticSnapshot(page)).features), 'Undo hides it again').not.toContain('f841fb8a8');
  await history(page, 'redo');
  await settleLive(page);
  expect(Object.keys((await semanticSnapshot(page)).features), 'Redo draws it again').toContain('f841fb8a8');
});

// OV-2326: a Result restored by a Load and not displayed is taken as drawn with
// what its drawing holds at the Load: a switch alone asks for no rerender, and
// a Reset from the other drawing asks for it.
const loadedWithLinearShown = async (page, browser, testInfo) => {
  await openBothDrawings(page, { mode: 'circular', results: 'single', reflow: 'off' });
  await addColorRule(page, FL1_ALPHA);
  await popupEdit(page, 'FL2', { visibility: 'off' });
  await generate(page);
  await switchDrawing(page, 'linear');
  const saved = testInfo.outputPath('reset-after-load.gbdraw-session.json');
  await download(page, 'Save Session', saved, { title: '' });
  return load(browser, saved);
};
test('a switch to a Result a Load restored asks for no rerender', async ({ page, browser }, testInfo) => {
  test.setTimeout(240_000);
  const loaded = await loadedWithLinearShown(page, browser, testInfo);
  const rerenderRequests = () => loaded.evaluate(async () => (await import('/gbdraw/web/js/state.js')).state.labelReflowForceRequestSeq.value);
  const before = await rerenderRequests();
  await switchDrawing(loaded, 'circular');
  expect(await rerenderRequests()).toBe(before);
  await loaded.context().close();
});
test('Reset Settings from the other drawing after a Load shows the removed color rule and hidden feature as Generate draws', async ({ page, browser }, testInfo) => {
  test.setTimeout(240_000);
  const loaded = await loadedWithLinearShown(page, browser, testInfo);
  await resetKeepingLabelSettings(loaded);
  await switchDrawing(loaded, 'circular');
  await expectLiveEqualsGenerate(loaded, { label: 'Reset Settings from the other drawing after a Load' });
  await loaded.context().close();
});

// OV-2293: when a displayed Result asks for the rerender. A drawn label edit left
// its drawing while the Result was not shown, on a feature that Result draws, or
// a feature it lacks is drawn again; once, with Auto Reflow on or off. Nothing
// else: an edit the Result followed live, or that another Result draws, needs no
// Python run at the display.
test.describe('a Result displayed after its intent changed', () => {
  const BATCH_FEATURE = { biological_feature_id: 'f5aa08bf5', record_idx: 1 };
  const rerenderRequests = (page) => page.evaluate(async () => {
    const { state } = await import('/gbdraw/web/js/state.js');
    return state.labelReflowRequestSeq.value + state.labelReflowForceRequestSeq.value;
  });
  const expectNoRerenderOfDisplay = async (page, display) => {
    await settleLive(page);
    const before = await rerenderRequests(page);
    await display();
    await settleLive(page);
    expect(await rerenderRequests(page), 'the display asks for no rerender').toBe(before);
  };
  const failRerenders = (page) => page.evaluate(() => {
    window.__GBDRAW_TEST_HOOKS__ ||= {};
    window.__GBDRAW_TEST_HOOKS__.beforeDiagramGenerationResponse = () => { throw new Error('forced rerender failure'); };
  });
  const allowRerenders = (page) => page.evaluate(() => { delete window.__GBDRAW_TEST_HOOKS__.beforeDiagramGenerationResponse; });
  const saveAndLoad = async (page, browser, testInfo, name) => {
    const saved = testInfo.outputPath(`${name}.gbdraw-session.json`);
    await download(page, 'Save Session', saved, { title: '' });
    return load(browser, saved);
  };

  for (const reflow of ['off', 'on']) {
    const flag = `Auto Reflow ${reflow}`;
    test(`a label text cleared live after a Generate asks for no rerender at a Result switch (${flag})`, async ({ page }) => {
      test.setTimeout(240_000);
      await open(page, { mode: 'circular', results: 'batch', reflow });
      await showResult(page, 1);
      await popupEdit(page, BATCH_FEATURE, { labelText: 'walk label' });
      await generate(page);
      await showResult(page, 1);
      await popupEdit(page, BATCH_FEATURE, { labelText: '' });
      await expectNoRerenderOfDisplay(page, () => showResult(page, 0));
      await expectNoRerenderOfDisplay(page, () => showResult(page, 1));
      await expectLiveEqualsGenerate(page, { label: 'the text cleared on its own Result, then a Result switch and back' });
    });
    test(`a live visibility edit asks for no rerender at a Result switch (${flag})`, async ({ page }) => {
      test.setTimeout(240_000);
      await open(page, { mode: 'circular', results: 'batch', reflow });
      await showResult(page, 1);
      await popupEdit(page, BATCH_FEATURE, { visibility: 'off' });
      await expectNoRerenderOfDisplay(page, () => showResult(page, 0));
      await expectNoRerenderOfDisplay(page, () => showResult(page, 1));
    });
    test(`an Undo of a live label text and of its Generate asks for no rerender at a Result switch (${flag})`, async ({ page }) => {
      test.setTimeout(240_000);
      await open(page, { mode: 'circular', results: 'batch', reflow });
      await showResult(page, 1);
      await popupEdit(page, BATCH_FEATURE, { labelText: 'walk label' });
      await generate(page);
      await history(page, 'undo');
      await history(page, 'undo');
      await settleLive(page);
      const shown = await page.evaluate(() => window.__GBDRAW_APP__.selectedResultIndex);
      await expectNoRerenderOfDisplay(page, () => showResult(page, shown === 0 ? 1 : 0));
    });
    // With Auto Reflow off the bytes hold the source text beside the edit, so the
    // display projects it. With it on, Python drew the edit: the display draws again.
    test(`an Undo of a live label text while another Result is shown shows what Generate draws (${flag})`, async ({ page }) => {
      test.setTimeout(240_000);
      await open(page, { mode: 'circular', results: 'batch', reflow });
      await showResult(page, 1);
      await popupEdit(page, BATCH_FEATURE, { labelText: 'x' });
      await settleLive(page);
      await showResult(page, 0);
      await history(page, 'undo');
      if (reflow === 'off') await expectNoRerenderOfDisplay(page, () => showResult(page, 1));
      else await showResult(page, 1);
      await expectLiveEqualsGenerate(page, { label: 'an Undo of a live text while another Result was shown' });
    });
    for (const mode of ['linear', 'circular']) {
      const feature = 'FL1';
      const both = `${mode}, ${flag}`;
      test(`a Save and Load with a label text cleared after the Generate asks for no rerender at the switch (${both})`, async ({ page, browser }, testInfo) => {
        test.setTimeout(240_000);
        await openBothDrawings(page, { mode, results: 'single', reflow });
        await popupEdit(page, feature, { labelText: 'walk label' });
        await generate(page);
        await popupEdit(page, feature, { labelText: '' });
        await settleLive(page);
        await switchDrawing(page, otherMode(mode));
        const loaded = await saveAndLoad(page, browser, testInfo, 'cleared-label');
        await expectNoRerenderOfDisplay(loaded, () => switchDrawing(loaded, mode));
        await loaded.context().close();
      });
      test(`Save and Load after Reset Settings from the other drawing draws the hidden feature when its Result is shown (${both})`, async ({ page, browser }, testInfo) => {
        test.setTimeout(240_000);
        await openBothDrawings(page, { mode, results: 'single', reflow });
        await popupEdit(page, 'FL2', { visibility: 'off' });
        await generate(page);
        await switchDrawing(page, otherMode(mode));
        await resetKeepingLabelSettings(page);
        const loaded = await saveAndLoad(page, browser, testInfo, 'reset-saved');
        await switchDrawing(loaded, mode);
        await expectLiveEqualsGenerate(loaded, { label: 'Reset, Save, Load, then the other drawing' });
        await loaded.context().close();
      });
      // What stays until Generate: nothing in a loaded Session says what the bytes
      // of a Result hold, and a Load alone moves no label.
      test(`Save and Load after Reset Settings from the other drawing leave a label text until the next Generate (${both})`, async ({ page, browser }, testInfo) => {
        test.setTimeout(240_000);
        await openBothDrawings(page, { mode, results: 'single', reflow });
        await popupEdit(page, feature, { labelText: 'walk label' });
        await generate(page);
        await switchDrawing(page, otherMode(mode));
        await resetKeepingLabelSettings(page);
        const loaded = await saveAndLoad(page, browser, testInfo, 'reset-saved-label');
        await expectNoRerenderOfDisplay(loaded, () => switchDrawing(loaded, mode));
        const drawnTexts = () => loaded.evaluate(() => [...window.__GBDRAW_APP__.svgContainer.querySelectorAll('svg text[data-label-feature-id]')]
          .map((element) => element.textContent));
        expect(await drawnTexts(), 'the label text stays').toContain('walk label');
        await generate(loaded);
        expect(await drawnTexts(), 'Generate draws the reset text').not.toContain('walk label');
        await loaded.context().close();
      });
      test(`Reset Settings from the other drawing after a Load reaches a label text the Result draws (${both})`, async ({ page, browser }, testInfo) => {
        test.setTimeout(240_000);
        await openBothDrawings(page, { mode, results: 'single', reflow });
        await popupEdit(page, feature, { labelText: 'walk label' });
        await generate(page);
        await switchDrawing(page, otherMode(mode));
        const loaded = await saveAndLoad(page, browser, testInfo, 'reset-after-load-label');
        await resetKeepingLabelSettings(loaded);
        await switchDrawing(loaded, mode);
        await expectLiveEqualsGenerate(loaded, { label: 'a label text in a loaded Result, then Reset from the other drawing' });
        await loaded.context().close();
      });
      test(`an Undo of a color rule step keeps what a later Reset from the other drawing asks of the Result (${both})`, async ({ page }) => {
        test.setTimeout(240_000);
        await openBothDrawings(page, { mode, results: 'single', reflow });
        await popupEdit(page, feature, { labelText: 'walk label' });
        await generate(page);
        await addColorRule(page, { feat: 'CDS', qual: 'locus_tag', val: '^FL2$', color: '#2a9d8f', cap: 'beta' });
        await settleLive(page);
        await history(page, 'undo');
        await settleLive(page);
        await resetFromOtherDrawing(page, mode);
        await expectLiveEqualsGenerate(page, { label: 'a checkpoint Undo, then Reset from the other drawing' });
      });
      test(`an Undo of Reset Settings from the other drawing after its Result was shown draws the edit again (${both})`, async ({ page }) => {
        test.setTimeout(240_000);
        await openBothDrawings(page, { mode, results: 'single', reflow });
        await popupEdit(page, feature, { labelText: 'walk label' });
        await generate(page);
        await resetFromOtherDrawing(page, mode);
        await settleLive(page);
        await switchDrawing(page, otherMode(mode));
        await history(page, 'undo');
        await settleLive(page);
        await switchDrawing(page, mode);
        await expectLiveEqualsGenerate(page, { label: 'an Undo of the Reset after the Result was shown' });
      });
      test(`an Undo of Reset Settings before its Result was shown asks for no rerender (${both})`, async ({ page }) => {
        test.setTimeout(240_000);
        await openBothDrawings(page, { mode, results: 'single', reflow });
        await popupEdit(page, feature, { labelText: 'walk label' });
        await generate(page);
        await switchDrawing(page, otherMode(mode));
        await resetKeepingLabelSettings(page);
        await settleLive(page);
        await history(page, 'undo');
        await expectNoRerenderOfDisplay(page, () => switchDrawing(page, mode));
      });
      for (const [name, text] of [['blank', ''], ['edited', 'bulk label']]) {
        test(`Reset Settings from the other drawing after a Label TSV ${name} text and a Generate shows what Generate draws (${both})`, async ({ page }) => {
          test.setTimeout(240_000);
          await openBothDrawings(page, { mode, results: 'single', reflow });
          await evaluateWithRetainedPromise(page, (bulkText) => window.__GBDRAW_APP__.loadLabelOverrideTable({
            target: { files: [new File([`*\t*\tlabel\t^alpha protein$\t${bulkText}\n`], 'labels.tsv', { type: 'text/plain' })] }
          }), text);
          await settleLive(page);
          await generate(page);
          await resetFromOtherDrawing(page, mode);
          await expectLiveEqualsGenerate(page, { label: `a Label TSV ${name} text, then Reset from the other drawing` });
        });
      }
      for (const [name, edit] of [['hidden feature', { visibility: 'off' }], ['label text', { labelText: 'walk label' }]]) {
        test(`a rerender that failed on a display does not run again at the next display (${name}, ${both})`, async ({ page }) => {
          test.setTimeout(240_000);
          await openBothDrawings(page, { mode, results: 'single', reflow });
          await popupEdit(page, 'FL2', edit);
          await generate(page);
          await switchDrawing(page, otherMode(mode));
          await resetKeepingLabelSettings(page);
          await failRerenders(page);
          const before = await rerenderRequests(page);
          await switchDrawing(page, mode);
          await expect.poll(() => rerenderRequests(page), 'the display asks once').toBe(before + 1);
          await settleLive(page);
          await allowRerenders(page);
          expect(await page.evaluate(() => window.__GBDRAW_APP__.labelReflowLastError), 'the failure is reported').not.toBeNull();
          await switchDrawing(page, otherMode(mode));
          await expectNoRerenderOfDisplay(page, () => switchDrawing(page, mode));
        });
      }
    }
    test(`Reset Settings with a Result an Undo of a later Generate returned reaches what its bytes drew (${flag})`, async ({ page }) => {
      test.setTimeout(240_000);
      await open(page, { mode: 'circular', results: 'batch', reflow });
      await showResult(page, 1);
      await popupEdit(page, BATCH_FEATURE, { visibility: 'off' });
      await generate(page);
      await generate(page);
      await history(page, 'undo');
      await settleLive(page);
      await showResult(page, 0);
      await resetKeepingLabelSettings(page);
      await showResult(page, 1);
      await expectLiveEqualsGenerate(page, { label: 'a hide in a Result an Undo of a Generate returned' });
    });
    test(`Reset Settings with a label text in a Result an Undo of a later Generate returned shows what Generate draws (${flag})`, async ({ page }) => {
      test.setTimeout(240_000);
      await open(page, { mode: 'circular', results: 'batch', reflow });
      await showResult(page, 1);
      await popupEdit(page, BATCH_FEATURE, { labelText: 'walk label' });
      await generate(page);
      await generate(page);
      await history(page, 'undo');
      await settleLive(page);
      await showResult(page, 0);
      await resetKeepingLabelSettings(page);
      await showResult(page, 1);
      await expectLiveEqualsGenerate(page, { label: 'a label text in a Result an Undo of a Generate returned' });
    });
    // An Undo of a Generate says nothing of what the Result it returns was drawn
    // with beyond the state it restores, so this text stays until the next Generate.
    test(`an Undo of a Generate that replaced an unshown Result holding a label text shows what Generate draws (${flag})`, async ({ page }) => {
      test.fail(true, 'OV-2293: a label edit that left before the Generate an Undo takes back');
      test.setTimeout(240_000);
      await open(page, { mode: 'circular', results: 'batch', reflow });
      await showResult(page, 1);
      await popupEdit(page, BATCH_FEATURE, { labelText: 'walk label' });
      await generate(page);
      await showResult(page, 0);
      await resetKeepingLabelSettings(page);
      await settleLive(page);
      await generate(page);
      await history(page, 'undo');
      await settleLive(page);
      // That Undo returns the state Reset wrote.
      await setAutoReflow(page, reflow === 'on');
      await showResult(page, 1);
      await expectLiveEqualsGenerate(page, { label: 'a Result an Undo of a Generate returned, its edit gone before that Generate' });
    });
    test(`an Undo of a Generate that replaced an unshown Result without a hidden feature draws the feature (${flag})`, async ({ page }) => {
      test.setTimeout(240_000);
      await open(page, { mode: 'circular', results: 'batch', reflow });
      await showResult(page, 1);
      await popupEdit(page, BATCH_FEATURE, { visibility: 'off' });
      await generate(page);
      await showResult(page, 0);
      await resetKeepingLabelSettings(page);
      await settleLive(page);
      await generate(page);
      await history(page, 'undo');
      await settleLive(page);
      // That Undo returns the state Reset wrote.
      await setAutoReflow(page, reflow === 'on');
      await showResult(page, 1);
      await expectLiveEqualsGenerate(page, { label: 'a hide in a Result an Undo of a Generate returned, gone before that Generate' });
    });
    for (const [name, change] of [['hidden feature', { visibility: 'off' }], ['label text', { labelText: 'walk label' }]]) {
      test(`a rerender that failed on a Result display does not run again at the next display (${name}, ${flag})`, async ({ page }) => {
        test.setTimeout(240_000);
        await open(page, { mode: 'circular', results: 'batch', reflow });
        await showResult(page, 1);
        await popupEdit(page, BATCH_FEATURE, change);
        await generate(page);
        await showResult(page, 0);
        await resetKeepingLabelSettings(page);
        await settleLive(page);
        await failRerenders(page);
        const before = await rerenderRequests(page);
        await showResult(page, 1);
        await expect.poll(() => rerenderRequests(page), 'the display asks once').toBe(before + 1);
        await settleLive(page);
        await allowRerenders(page);
        expect(await page.evaluate(() => window.__GBDRAW_APP__.labelReflowLastError), 'the failure is reported').not.toBeNull();
        await expectNoRerenderOfDisplay(page, () => showResult(page, 0));
        await expectNoRerenderOfDisplay(page, () => showResult(page, 1));
      });
    }
    test(`a rerender that failed on an edit does not run again when another batch Result is displayed (${flag})`, async ({ page }) => {
      test.setTimeout(240_000);
      await open(page, { mode: 'circular', results: 'batch', reflow });
      await showResult(page, 1);
      await failRerenders(page);
      await popupEdit(page, BATCH_FEATURE, { labelVisibility: 'off' });
      await popupEdit(page, BATCH_FEATURE, { labelVisibility: 'default' });
      await settleLive(page);
      await allowRerenders(page);
      await expectNoRerenderOfDisplay(page, () => showResult(page, 0));
    });
  }
  // A Load alone moves no label of a cropped Result: its features outside the crop ask nothing.
  test('a switch to a cropped Result a Load restored asks for no rerender', async ({ page, browser }, testInfo) => {
    test.setTimeout(240_000);
    await open(page, { mode: 'linear', results: 'single', reflow: 'off' });
    await page.evaluate(() => { Object.assign(window.__GBDRAW_APP__.linearSeqs[0], { region_start: 1, region_end: 1000 }); });
    await generate(page);
    await switchDrawing(page, 'circular');
    await page.getByLabel('GenBank/DDBJ File', { exact: true }).setInputFiles(SINGLE_FIXTURE);
    await expect.poll(() => page.evaluate(() => window.__GBDRAW_APP__.circularRecordList.length), { timeout: 60_000 }).toBeGreaterThan(0);
    await settleLive(page);
    await generate(page);
    const loaded = await saveAndLoad(page, browser, testInfo, 'cropped');
    await expectNoRerenderOfDisplay(loaded, () => switchDrawing(loaded, 'linear'));
    await loaded.context().close();
  });
  // Undo or Redo back to Label visibility On, on a figure that lacks the label:
  // only a rerender can draw it, and History steps ask for none here (OV-2413).
  for (const mode of ['linear', 'circular']) {
    const labelsNoneOn = async (page) => {
      await open(page, { mode, results: 'single', reflow: 'off' });
      await showNoLabels(page);
      await popupEdit(page, 'FL1', { labelVisibility: 'on' });
      await settleLive(page);
    };
    test(`an Undo back to Label visibility On draws the label (${mode})`, async ({ page }) => {
      test.fail(true, 'OV-2413: an Undo back to Label visibility On on a figure without the label');
      test.setTimeout(240_000);
      await labelsNoneOn(page);
      await popupEdit(page, 'FL1', { labelVisibility: 'default' });
      await settleLive(page);
      await history(page, 'undo');
      await settleLive(page);
      await expectLiveEqualsGenerate(page, { label: 'an Undo back to On' });
    });
    test(`a Redo of Label visibility On draws the label (${mode})`, async ({ page }) => {
      test.fail(true, 'OV-2413: a Redo of Label visibility On on a figure without the label');
      test.setTimeout(240_000);
      await labelsNoneOn(page);
      await history(page, 'undo');
      await settleLive(page);
      await history(page, 'redo');
      await settleLive(page);
      await expectLiveEqualsGenerate(page, { label: 'a Redo of On' });
    });
  }
  // Text typed in the Features list stays through a rerender, as in the popup.
  test('text typed in the Features list label input survives a rerender', async ({ page }) => {
    test.setTimeout(180_000);
    await open(page, { mode: 'linear', results: 'single', reflow: 'on' });
    await page.getByRole('button', { name: 'Editor', exact: true }).click();
    await page.evaluate(() => {
      const app = window.__GBDRAW_APP__;
      app.featurePanelTab = 'features';
      app.featureListExpandedId = app.filteredFeatures.find((feature) => feature.locus_tag === 'FL2').id;
    });
    const input = page.locator('[data-feature-list-editor] input[type="text"]');
    await input.waitFor();
    await holdRerenders(page);
    const edit = popupEdit(page, 'FL1', { labelText: 'held FL1' });
    await rerenderHeld(page);
    await page.evaluate(() => {
      const app = window.__GBDRAW_APP__;
      app.featureListExpandedId = app.filteredFeatures.find((feature) => feature.locus_tag === 'FL2').id;
    });
    await input.click();
    await page.keyboard.press('Control+A');
    await page.keyboard.type('typed in the list');
    await releaseRerenders(page);
    await edit;
    await expect(input).toHaveValue('typed in the list');
    await input.blur();
    await page.evaluate(() => window.__GBDRAW_APP__.handleLabelTextScopeChoice('single'));
    await settleLive(page);
    expect(await page.evaluate(() => window.__GBDRAW_APP__.getEditableLabelByFeatureId(
      window.__GBDRAW_APP__.extractedFeatures.find((feature) => feature.locus_tag === 'FL2').svg_id).text)).toBe('typed in the list');
    await expectLiveEqualsGenerate(page, { label: 'the text typed during a rerender' });
  });
});

// Undo and Redo of a feature hide, a Label Off, or a label text apply at once with
// Auto Reflow off: when Undo or Redo returns, the figure shows the step, with any
// rerender held. Guards of how History steps work, which also keep the moved
// Legend and the displayed batch Result where they are. With Auto Reflow on, two
// of them are confirmed mismatches, marked with their findings.
test.describe('Undo and Redo', () => {
  const featureIdOf = (page, locusTag) => page.evaluate((tag) => {
    const app = window.__GBDRAW_APP__;
    return app.extractedFeatures.find((feature) => feature.locus_tag === tag)?.svg_id;
  }, locusTag);
  // The step, and what the figure shows of the feature as it returns.
  const stepAndRead = (page, name, id) => evaluateWithRetainedPromise(page, async ({ name: step, id: featureId }) => {
    await window.__GBDRAW_HISTORY__[step]();
    await window.Vue.nextTick();
    const root = window.__GBDRAW_APP__.svgContainer.querySelector('svg');
    const shown = (element) => !element.closest('[display="none"]');
    return {
      feature: [...root.querySelectorAll(`[data-gbdraw-feature-id="${CSS.escape(featureId)}"]`)].filter(shown).length > 0,
      labels: [...root.querySelectorAll(`text[data-label-feature-id="${CSS.escape(featureId)}"]`)].filter(shown).map((label) => label.textContent)
    };
  }, { name, id });
  const SHOWN = { feature: true, labels: ['beta protein'] };
  const AT_ONCE_CASES = [
    { domain: 'a feature hide', edit: (page) => popupEdit(page, 'FL2', { visibility: 'off' }), edited: { feature: false, labels: [] } },
    { domain: 'a Label Off', edit: (page) => popupEdit(page, 'FL2', { labelVisibility: 'off' }), edited: { feature: true, labels: [] } },
    { domain: 'a label text', edit: (page) => popupEdit(page, 'FL2', { labelText: 'at once text' }), edited: { feature: true, labels: ['at once text'] } },
    {
      domain: 'a feature hide by a visibility rule',
      edit: (page) => addVisibilityRule(page, { recordId: '*', featureType: 'CDS', qualifier: 'locus_tag', value: '^FL2$', action: 'off' }),
      edited: { feature: false, labels: [] }
    }
  ];
  const rerenderRequests = (page) => page.evaluate(async () => {
    const { state } = await import('/gbdraw/web/js/state.js');
    return state.labelReflowRequestSeq.value + state.labelReflowForceRequestSeq.value;
  });
  const movedLegend = (page) => page.evaluate(async () => {
    const { compositionUserDeltas } = await import('/gbdraw/web/js/app/legend-layout/composition-actions.js');
    return compositionUserDeltas(window.__GBDRAW_APP__.svgContainer.querySelector('svg'));
  });
  const dragLegend = (page) => page.evaluate(async () => {
    const app = window.__GBDRAW_APP__;
    app.layoutRepositionMode = true;
    await window.Vue.nextTick();
    const svg = app.svgContainer.querySelector('svg');
    const parts = [...svg.querySelectorAll('[data-gbdraw-composition-role="legend"]')];
    const target = parts.find((element) => element.id !== 'length_bar') || parts[0];
    const box = target.getBoundingClientRect();
    const [x, y] = [box.left + Math.min(box.width / 2, 8), box.top + Math.min(box.height / 2, 8)];
    const mouse = (type, clientX, clientY, buttons) => new MouseEvent(type, { bubbles: true, cancelable: true, clientX, clientY, buttons, view: window });
    const frame = () => new Promise((resolve) => { requestAnimationFrame(resolve); });
    target.dispatchEvent(mouse('mousedown', x, y, 1));
    await frame();
    svg.dispatchEvent(mouse('mousemove', x + 30, y + 20, 1));
    await frame();
    svg.dispatchEvent(mouse('mouseup', x + 30, y + 20, 0));
    await new Promise((resolve) => { setTimeout(resolve, 100); });
    app.layoutRepositionMode = false;
  });
  test('a moved Legend returns through Undo after a label edit, Auto Reflow on (linear)', async ({ page }) => {
    test.fail(true, 'OV-2437: the Undo of the drag after a label edit does nothing with Auto Reflow on');
    test.setTimeout(180_000);
    await open(page, { mode: 'linear', results: 'single', reflow: 'on' });
    const start = await movedLegend(page);
    await dragLegend(page);
    await settleLive(page);
    await popupEdit(page, 'FL2', { labelText: 'beside the Legend' });
    await settleLive(page);
    await history(page, 'undo');
    await settleLive(page);
    await history(page, 'undo');
    await settleLive(page);
    expect(await movedLegend(page), 'Undo of the drag returns it').toEqual(start);
  });
  for (const mode of ['linear', 'circular']) {
    for (const { domain, edit, edited } of AT_ONCE_CASES) {
      test(`Undo and Redo of ${domain} apply at once with no rerender, Auto Reflow off (${mode})`, async ({ page }) => {
        test.setTimeout(180_000);
        await open(page, { mode, results: 'single', reflow: 'off' });
        const id = await featureIdOf(page, 'FL2');
        await edit(page);
        await settleLive(page);
        const before = await rerenderRequests(page);
        await holdRerenders(page);
        expect(await stepAndRead(page, 'undo', id), 'Undo').toEqual(SHOWN);
        expect(await stepAndRead(page, 'redo', id), 'Redo').toEqual(edited);
        expect(await stepAndRead(page, 'undo', id), 'Undo again').toEqual(SHOWN);
        await releaseRerenders(page);
        expect(await rerenderRequests(page), 'the steps ask for no rerender').toBe(before);
        await expectLiveEqualsGenerate(page, { label: `Undo of ${domain}` });
      });
    }
    // With Auto Reflow on, Undo and Redo of a hide or a Label Off wait for the rerender.
    for (const { domain, edit, edited } of AT_ONCE_CASES.slice(0, 2)) {
      test(`Undo and Redo of ${domain} apply at once, Auto Reflow on (${mode})`, async ({ page }) => {
        test.fail(true, 'OV-2436: Undo and Redo of a hide or a Label Off wait for the rerender with Auto Reflow on');
        test.setTimeout(180_000);
        await open(page, { mode, results: 'single', reflow: 'on' });
        const id = await featureIdOf(page, 'FL2');
        await edit(page);
        await settleLive(page);
        await holdRerenders(page);
        const undo = await stepAndRead(page, 'undo', id);
        const redo = await stepAndRead(page, 'redo', id);
        await releaseRerenders(page);
        await settleLive(page);
        expect(undo, 'Undo').toEqual(SHOWN);
        expect(redo, 'Redo').toEqual(edited);
      });
    }
    test(`a moved Legend returns through Undo of a label edit, and comes back through Redo, Auto Reflow off (${mode})`, async ({ page }) => {
      test.setTimeout(180_000);
      await open(page, { mode, results: 'single', reflow: 'off' });
      const moved = () => movedLegend(page);
      const start = await moved();
      await dragLegend(page);
      await settleLive(page);
      const dragged = await moved();
      expect(dragged, 'the drag moves the Legend').not.toEqual(start);
      await popupEdit(page, 'FL2', { labelText: 'beside the Legend' });
      await settleLive(page);
      await history(page, 'undo');
      expect(await moved(), 'Undo of the label edit keeps the Legend').toEqual(dragged);
      await history(page, 'undo');
      expect(await moved(), 'Undo of the drag returns it').toEqual(start);
      await history(page, 'redo');
      expect(await moved(), 'Redo of the drag').toEqual(dragged);
    });
  }
  for (const [name, change] of [['a hide', { visibility: 'off' }], ['a label text', { labelText: 'on result 1' }]]) {
    test(`Undo and Redo of ${name} leave another batch Result displayed`, async ({ page }) => {
      test.setTimeout(240_000);
      await open(page, { mode: 'circular', results: 'batch', reflow: 'off' });
      const shownResult = () => page.evaluate(() => window.__GBDRAW_APP__.selectedResultIndex);
      await showResult(page, 1);
      await popupEdit(page, { biological_feature_id: 'f5aa08bf5', record_idx: 1 }, change);
      await settleLive(page);
      await showResult(page, 0);
      await history(page, 'undo');
      await settleLive(page);
      expect(await shownResult(), 'Undo').toBe(0);
      await history(page, 'redo');
      await settleLive(page);
      expect(await shownResult(), 'Redo').toBe(0);
    });
  }
});

// OV-2410 (Circular here; Linear below): Generate draws this label outside its arc with a leader for a text
// that does not fit, and inside it without one for a text that does. The label
// text edit places the label only with Auto Reflow on, as documented; with it
// off the label keeps its place until a reflow or Generate. That place also comes
// back with an Undo of the Generate while Auto Reflow is on, which the allowed
// difference does not cover.
const LEADER_FEATURE = { biological_feature_id: 'f5aa08bf5', record_idx: 1 };
const LEADER_TEXT = 'alpha protein, edited live';
const openLeaderCase = async (page, reflow) => {
  await open(page, { mode: 'circular', results: 'batch', reflow });
  await showResult(page, LEADER_FEATURE.record_idx);
};
const leaderDifferences = ({ tolerated }) => tolerated.map(({ kind, id, attribute, live, generated }) => ({ kind, id, attribute, live, generated }));
const LEADER_AWAITS_REFLOW = [{ kind: 'leader', id: 'f5aa08bf5', attribute: 'drawn', live: false, generated: true }];
const setAutoReflow = (page, enabled) => page.evaluate((value) => { window.__GBDRAW_APP__.autoLabelReflowEnabled = value; }, enabled);
const undoGenerate = async (page) => {
  await generate(page);
  await showResult(page, LEADER_FEATURE.record_idx);
  await history(page, 'undo');
};
test('a label text edit that moves the label outside its arc draws its leader with Auto Reflow on', async ({ page }) => {
  test.setTimeout(240_000);
  await openLeaderCase(page, 'on');
  await popupEdit(page, LEADER_FEATURE, { labelText: LEADER_TEXT });
  await expectLiveEqualsGenerate(page, { allowed: [], label: 'the edit with Auto Reflow on' });
  await undoGenerate(page);
  await expectLiveEqualsGenerate(page, { allowed: [], label: 'the Undo of the Generate with Auto Reflow on' });
});
test('a label text edit with Auto Reflow off leaves only its leader to the next placement', async ({ page }) => {
  test.setTimeout(240_000);
  await openLeaderCase(page, 'off');
  await popupEdit(page, LEADER_FEATURE, { labelText: LEADER_TEXT });
  expect(leaderDifferences(await expectLiveEqualsGenerate(page, { label: 'the edit with Auto Reflow off' }))).toEqual(LEADER_AWAITS_REFLOW);
});
test('an Undo of a Generate returns labels placed for the state before it when Auto Reflow was turned on meanwhile', async ({ page }) => {
  test.fail(true, 'OV-2410: the leader of a label text placed with Auto Reflow off stays through this Undo');
  test.setTimeout(240_000);
  await openLeaderCase(page, 'off');
  await popupEdit(page, LEADER_FEATURE, { labelText: LEADER_TEXT });
  await generate(page);
  await showResult(page, LEADER_FEATURE.record_idx);
  await setAutoReflow(page, true);
  await history(page, 'undo');
  await expectLiveEqualsGenerate(page, { label: 'the Undo with Auto Reflow on' });
});
test('a label text edit on a Linear diagram draws the leader Generate draws, Auto Reflow off', async ({ page }) => {
  test.fail(true, 'OV-2410: a label text edit with Auto Reflow off leaves the leader to the next placement');
  test.setTimeout(240_000);
  await open(page, { mode: 'linear', results: 'single', reflow: 'off' });
  await popupEdit(page, 'FL1', { labelText: 'alpha protein, edited live with a much longer text' });
  await expectLiveEqualsGenerate(page, { allowed: [], label: 'the edit on a Linear diagram' });
});
// The random walk's allowance for an Undo covers the edited label alone: a leader
// lost on any other label is a difference. The test removes one to show it.
test('a leader lost on another label fails where the edited label is left unplaced', async ({ page }) => {
  test.setTimeout(240_000);
  await openLeaderCase(page, 'off');
  await popupEdit(page, LEADER_FEATURE, { labelText: 'x' });
  await generate(page);
  await showResult(page, LEADER_FEATURE.record_idx);
  await setAutoReflow(page, true);
  await history(page, 'undo');
  await settleLive(page);
  const hidden = await page.evaluate((edited) => {
    const root = window.__GBDRAW_APP__.svgContainer.querySelector('svg');
    const others = [...root.querySelectorAll('[data-label-feature-id]')]
      .filter((element) => element.localName !== 'text' && element.getAttribute('data-label-feature-id') !== edited);
    others.forEach((element) => element.setAttribute('display', 'none'));
    return [...new Set(others.map((element) => element.getAttribute('data-label-feature-id')))];
  }, LEADER_FEATURE.biological_feature_id);
  expect(hidden, 'another label draws a leader').toHaveLength(1);
  const unexpected = await expectLiveEqualsGenerate(page, {
    label: 'the Undo with one other leader removed',
    unplaced: { id: LEADER_FEATURE.biological_feature_id, text: null }
  }).then(() => [], (error) => error?.matcherResult?.actual);
  expect(unexpected).toEqual([`leader ${hidden[0]} drawn: live false, generated true`]);
});
test('a label text that fits its arc draws no leader live or from Generate, Auto Reflow off', async ({ page }) => {
  test.setTimeout(240_000);
  await openLeaderCase(page, 'off');
  await popupEdit(page, LEADER_FEATURE, { labelText: 'x' });
  await expectLiveEqualsGenerate(page, { allowed: [], label: 'a text that fits' });
});

// A Reset that changes the palette: the palette watcher shows the default
// palette's fills, and Reset projects the other domains.
test('Reset Settings after a palette change and a Legend row stroke shows what Generate draws', async ({ page }) => {
  test.setTimeout(120_000);
  await open(page, { mode: 'circular', results: 'single', reflow: 'off' });
  await resetSettings(page);
  const palette = await page.evaluate(() => window.__GBDRAW_APP__.paletteNames
    .find((name) => name !== window.__GBDRAW_APP__.selectedPalette));
  await page.locator('summary[aria-label="Colors"]').click();
  await page.getByRole('combobox', { name: 'Palette', exact: true }).selectOption(palette);
  await settleLive(page);
  await generate(page);
  const drawn = await semanticSnapshot(page);
  await legendRowStrokeColor(page, 'CDS', '#e63946');
  await resetSettings(page);
  expect(await page.evaluate(() => window.__GBDRAW_APP__.selectedPalette), 'Reset returns to the default palette')
    .not.toBe(palette);
  const reset = await semanticSnapshot(page);
  expect(reset.features, 'the default palette shows on the features').not.toEqual(drawn.features);
  const listed = await legendList(page);
  await expectLiveEqualsGenerate(page, { label: 'Reset Settings after a palette change' });
  expect(listed, 'the Legend list after Reset Settings').toEqual(await legendList(page));
});

// OV-1385: a palette color of None draws the Legend row's swatch with fill
// none, as Python draws it, and the next live edits of that row still find
// the swatch: no page error, and each step shows what Generate draws (Linear,
// two genomes, as reported).
test('Palette color None and the live edits after it show the Legend row Generate draws', async ({ page }) => {
  test.setTimeout(240_000);
  /** @type {string[]} */
  const errors = [];
  page.on('pageerror', (error) => errors.push(error.message));
  await openFresh(page);
  await page.evaluate(async ({ first, second }) => {
    const app = window.__GBDRAW_APP__;
    app.selectDrawing('linear');
    await window.Vue.nextTick();
    if (app.linearSeqs.length < 2) app.addLinearSeq();
    app.setLinearSeqPrimaryFile(0, 'gb', new File([first], 'MjeNMV.gb', { type: 'text/plain' }));
    app.setLinearSeqPrimaryFile(1, 'gb', new File([second], 'MelaMJNV.gb', { type: 'text/plain' }));
    await window.Vue.nextTick();
  }, { first: readFileSync('examples/MjeNMV.gb', 'utf8'), second: readFileSync('examples/MelaMJNV.gb', 'utf8') });
  await generate(page);
  const swatchFills = (key) => page.evaluate((row) => [...window.__GBDRAW_APP__.svgContainer
    .querySelectorAll(`g[data-legend-key="${CSS.escape(row)}"]`)]
    .map((entry) => [...entry.children].find((child) => child.localName === 'path')?.getAttribute('fill')), key);
  await page.locator('summary[aria-label="Colors"]').click();
  const mode = page.getByRole('combobox', { name: 'CDS feature color mode', exact: true });
  const setMode = async (value) => { await mode.selectOption(value); await settleLive(page); };
  await setMode('none');
  expect(await swatchFills('CDS'), 'the swatch Python draws for None').toEqual(['none', 'none']);
  await setMode('auto');
  expect((await swatchFills('CDS')).every((fill) => fill && fill !== 'none'), 'Auto fills the swatch again').toBe(true);
  await setMode('none');
  await expectLiveEqualsGenerate(page, { label: 'CDS color None after Auto' });
  expect(await swatchFills('CDS'), 'Generate draws the swatch with fill none').toEqual(['none', 'none']);
  await setMode('auto');
  await expectLiveEqualsGenerate(page, { label: 'CDS color Auto after a Generate of None' });
  expect(errors, 'no live edit throws').toEqual([]);
});

for (const { edit, states, setup, run, knownMismatch } of CASES) {
  test(`${edit}: ${statesName(states)}`, async ({ page }) => {
    if (knownMismatch) test.fail(true, knownMismatch);
    test.setTimeout(90_000);
    await open(page, states);
    if (setup) await setup(page);
    if (states.labels === 'bound') await popupEdit(page, BIND_TARGET[states.results]);
    await run(page);
    await expectLiveEqualsGenerate(page, { label: `${edit} (${statesName(states)})` });
  });
}

// OV-484: a request without explicit track slots draws its skew track under
// Python's preset slot id. A palette change repaints that track as the next
// Generate draws it (the semantic snapshot does not read track fills).
for (const mode of ['circular', 'linear']) {
  test(`Palette change repaints the skew track of a request without track slots: ${mode}`, async ({ page }) => {
    test.setTimeout(90_000);
    await open(page, { mode, results: 'single', reflow: 'off' });
    if (mode === 'linear') {
      await page.evaluate(() => { window.__GBDRAW_APP__.form.show_skew = true; });
      await generate(page);
    }
    const skewFills = () => page.evaluate(() => [...document.querySelectorAll('svg g[data-gbdraw-slot-renderer="dinucleotide_skew"]')]
      .map((group) => [...group.querySelectorAll('path')].map((path) => String(path.getAttribute('fill') || '').toLowerCase())
        .filter((fill) => fill && fill !== 'none' && fill !== 'white')));
    const drawn = await skewFills();
    expect(drawn.flat().length, 'the Result draws a skew track').toBeGreaterThan(0);
    await switchPalette(page, 'alpine_retreat');
    const live = await skewFills();
    expect(live, 'the palette repaints the skew track').not.toEqual(drawn);
    await generate(page);
    expect(live, 'the skew track of the live edit is the one Generate draws').toEqual(await skewFills());
  });
}

// A Session 30 (gbdraw 0.13.0) has no request: the Web shows its saved Result,
// described by the track settings Load recorded. A palette change repaints its
// skew track and the track's Legend rows as Generate then draws them.
test('Palette change repaints the skew track and its Legend rows of a loaded Session 30 Result', async ({ page }) => {
  test.setTimeout(240_000);
  await openFresh(page);
  // A CLI Session of 0.13.0 lists no extracted features, so the Load is awaited by its Result.
  await page.locator('input[accept^=".json,"]').setInputFiles('tests/fixtures/sessions/scale-interval-zero-circular-cli.v30.gbdraw-session.json.gz');
  await expectNoLoadRepairDialog(page);
  await page.waitForFunction(() => !window.__GBDRAW_APP__.sessionImportPending
    && window.__GBDRAW_APP__.svgContainer?.querySelector('g#skew'), null, { timeout: 60_000 });
  await settleLive(page);
  expect(await page.evaluate(async () => (await import('/gbdraw/web/js/services/config.js')).getCommittedCanonicalRenderRequest()),
    'the loaded Result has no request').toBeNull();
  const skew = () => page.evaluate(() => {
    const fills = (selector) => [...window.__GBDRAW_APP__.svgContainer.querySelectorAll(selector)]
      .map((shape) => String(shape.getAttribute('fill') || '').toLowerCase()).filter((fill) => fill && fill !== 'none' && fill !== 'white');
    return {
      track: fills('g#skew > path, g#gc_skew > path, g[data-gbdraw-slot-renderer="dinucleotide_skew"] path'),
      rows: ['GC skew (+)', 'GC skew (-)'].map((key) => fills(`g[data-legend-key="${key}"] :is(rect, path)`)[0])
    };
  });
  const drawn = await skew();
  expect(drawn.track.length, 'the loaded Result draws a skew track').toBe(2);
  expect(drawn.rows, 'the track\'s Legend rows show its fills').toEqual(drawn.track);
  await switchPalette(page, 'alpine_retreat');
  const live = await skew();
  expect(live.track, 'the palette repaints the skew track').not.toEqual(drawn.track);
  expect(live.rows, 'the palette repaints the track\'s Legend rows').toEqual(live.track);
  await generate(page);
  expect(live, 'the live edit shows what Generate draws').toEqual(await skew());
});

// OV-144 (D-07, PD-OI-062): the executor's record of Python's strokes travels
// with a Result through Save Session and Load, so a stroke edit removed after
// Load leaves a batch Result that showed it when that Result is shown again.
test('a stroke removed after Load leaves a batch Result that showed it (circular, two-Result batch)', async ({ page, browser }, testInfo) => {
  test.setTimeout(240_000);
  await open(page, { mode: 'circular', results: 'batch', reflow: 'off' });
  await legendRowStroke(page, 'CDS', '#e63946', 3);
  await showResult(page, 1);
  await showResult(page, 0);
  const saved = testInfo.outputPath('stroked-batch.gbdraw-session.json');
  await download(page, 'Save Session', saved, { title: '' });
  const loaded = await load(browser, saved);
  if (await loaded.evaluate(() => window.__GBDRAW_APP__.selectedResultIndex) !== 0) await showResult(loaded, 0);
  await loaded.evaluate((row) => {
    const app = window.__GBDRAW_APP__;
    return app.resetLegendEntryStroke(app.legendEntries.findIndex((entry) => entry.caption === row));
  }, 'CDS');
  await settleLive(loaded);
  await showResult(loaded, 1);
  await expectLiveEqualsGenerate(loaded, { label: 'a stroke removed after Load' });
  await loaded.context().close();
});

// U2BFIX review #1: Save keeps each batch Result's bytes as stored, so a Result
// not displayed since the paint edits is saved without them. Its first
// display after Load shows every paint domain, as Generate draws it.
test('a batch Result saved before it showed the paint edits shows them after Load (circular, two-Result batch)', async ({ page, browser }, testInfo) => {
  test.setTimeout(240_000);
  await open(page, { mode: 'circular', results: 'batch', reflow: 'off' });
  await legendRowStroke(page, 'CDS', '#e63946', 3);
  await popupEdit(page, 'TESTA_0001', { stroke: '#2a9d8f' });
  await legendRowColor(page, 'tRNA', '#7b2cbf');
  const saved = testInfo.outputPath('stale-batch.gbdraw-session.json');
  await download(page, 'Save Session', saved, { title: '' });
  const loaded = await load(browser, saved);
  if (await loaded.evaluate(() => window.__GBDRAW_APP__.selectedResultIndex) !== 0) await showResult(loaded, 0);
  await showResult(loaded, 1);
  await expectLiveEqualsGenerate(loaded, { label: 'a batch Result saved without the paint edits, after Load' });
  await loaded.context().close();
});

// The same for an Undo of Generate, which restores the Results generated
// before it with their bytes as they were kept.
test('a batch Result kept without the paint edits shows them after an Undo of Generate (circular, two-Result batch)', async ({ page }) => {
  test.setTimeout(240_000);
  await open(page, { mode: 'circular', results: 'batch', reflow: 'off' });
  await legendRowColor(page, 'tRNA', '#7b2cbf');
  await legendRowStroke(page, 'CDS', '#e63946', 3);
  await generate(page);
  await history(page, 'undo');
  await showResult(page, 1);
  await expectLiveEqualsGenerate(page, { label: 'a batch Result kept without the paint edits, after an Undo of Generate' });
});

// U2BFIX2 review M3: a batch Result whose display a Save started meanwhile
// declines shows the palette and visibility edits made on another Result once
// Save ends, while it stays displayed, so a later Save keeps them.
test('a batch Result displayed while Save runs shows the edits once Save ends (circular, two-Result batch)', async ({ page }) => {
  test.setTimeout(240_000);
  await open(page, { mode: 'circular', results: 'batch', reflow: 'off' });
  await legendRowColor(page, 'tRNA', '#7b2cbf');
  await addVisibilityRule(page, BATCH_0004_OFF);
  await displayWhileSaveHeld(page, 1);
  await expectLiveEqualsGenerate(page, { label: 'a batch Result displayed while Save ran' });
});

// U2BFIX3 review M-A: Save Session started while a batch Result is being
// displayed writes that Result with the edits, so Load shows it as Generate
// draws it. The color rule's rerender draws both Results again and leaves
// the rule matches to prepare: the Legend row color after it is the edit the
// displayed Result lacks, and the Save prepares the matches it needs. (The
// case above holds its Save with its own `beforeExport`, so it never reaches
// the app's Save.)
test('a batch Result displayed while Save runs is saved with the edits (circular, two-Result batch)', async ({ page, browser }, testInfo) => {
  test.setTimeout(240_000);
  await open(page, { mode: 'circular', results: 'batch', reflow: 'off' });
  await addVisibilityRule(page, BATCH_0004_OFF);
  await addColorRule(page, { feat: 'CDS', qual: 'locus_tag', val: '_0002$', color: '#2266aa', cap: 'CDS' });
  await legendRowColor(page, 'tRNA', '#7b2cbf');
  const [file] = await Promise.all([page.waitForEvent('download'), saveWhileDisplaying(page, 1)]);
  expect(await page.evaluate(() => window.__GBDRAW_APP__.selectedResultIndex), 'the Save raced the display').toBe(1);
  const saved = testInfo.outputPath('declined-display.gbdraw-session.json');
  await file.saveAs(saved);
  await settleLive(page);
  await expectLiveEqualsGenerate(page, { label: 'a batch Result displayed while Save ran' });
  const loaded = await load(browser, saved);
  expect(await loaded.evaluate(() => window.__GBDRAW_APP__.selectedResultIndex), 'Load shows the saved Result').toBe(1);
  await expectLiveEqualsGenerate(loaded, { label: 'a batch Result saved while it was displayed, after Load' });
  await loaded.context().close();
});

// U2BFIX2 review L1: a label edit removed while another Result was displayed
// leaves a batch Result saved with it; its first display after Load shows the
// label as Generate draws it.
test('a batch Result saved with a label edit removed meanwhile shows the label after Load (circular, two-Result batch)', async ({ page, browser }, testInfo) => {
  test.setTimeout(240_000);
  await open(page, { mode: 'circular', results: 'batch', reflow: 'off' });
  await showResult(page, 1);
  await popupEdit(page, 'TESTB_0001', { labelText: 'edited on B' });
  await showResult(page, 0);
  await appAction(page, 'resetAllLabelTextOverrides');
  const saved = testInfo.outputPath('label-batch.gbdraw-session.json');
  await download(page, 'Save Session', saved, { title: '' });
  const loaded = await load(browser, saved);
  if (await loaded.evaluate(() => window.__GBDRAW_APP__.selectedResultIndex) !== 0) await showResult(loaded, 0);
  await showResult(loaded, 1);
  await expectLiveEqualsGenerate(loaded, { label: 'a label edit removed before Save, after Load' });
  await loaded.context().close();
});

// OV-129, OV-150: a stroke reset or Undo returns each feature part and Legend
// swatch to the stroke Python drew for it: a connector line its
// `line_stroke_*`, a block its `block_stroke_*` of the genome size class.
// Python's strokes are read before the first edit; the stroke width counts
// too, which the semantic snapshot leaves out.
const drawnStrokes = (page) => page.evaluate(() => {
  const root = window.__GBDRAW_APP__.svgContainer?.querySelector('svg');
  const strokes = (element) => [element.getAttribute('stroke'), element.getAttribute('stroke-width')];
  const features = [...root.querySelectorAll('[data-gbdraw-feature-id]')].map((element) => [
    element.getAttribute('data-gbdraw-rendered-feature-id') || element.getAttribute('data-gbdraw-feature-id'),
    element.getAttribute('data-gbdraw-feature-part') || element.localName,
    ...strokes(element)
  ]);
  const swatches = [...root.querySelectorAll('g[data-legend-key]')].map((entry) => {
    const swatch = [...entry.querySelectorAll('path')].find((path) => {
      const fill = path.getAttribute('fill');
      return fill && fill !== 'none' && !fill.startsWith('url(');
    });
    return [entry.getAttribute('data-legend-key'), 'swatch', ...(swatch ? strokes(swatch) : [])];
  });
  return [...features, ...swatches];
});
const STROKE_RESET_CASES = {
  'popup Reset Stroke of a spliced feature': {
    finding: 'OV-150: the connector of a spliced feature takes the block stroke on Reset Stroke',
    edit: (page) => popupEdit(page, 'TESTA_0004', { stroke: '#e63946' }),
    revert: (page) => popupEdit(page, 'TESTA_0004', { resetStroke: true })
  },
  'Reset all strokes after a Legend row stroke': {
    finding: 'OV-129: Reset all strokes gives connectors and every size class the first feature stroke',
    edit: (page) => legendRowStroke(page, 'CDS', '#e63946', 3),
    revert: (page) => appAction(page, 'resetAllStrokes')
  },
  'Undo of a Legend row stroke': {
    finding: 'OV-129: the History stroke reconcile gives connectors the first feature stroke',
    edit: (page) => legendRowStrokeColor(page, 'CDS', '#e63946'),
    revert: (page) => history(page, 'undo')
  }
};
// The findings (OV-150, OV-129) are fixed: every stroke edit is shown by the
// executor, which records Python's stroke (EU U2a).
for (const [name, { edit, revert }] of Object.entries(STROKE_RESET_CASES)) {
  test(`${name} returns each feature part to Python's stroke (circular, two-Result batch)`, async ({ page }) => {
    test.setTimeout(120_000);
    await open(page, { mode: 'circular', results: 'batch', reflow: 'off' });
    const drawn = await drawnStrokes(page);
    await edit(page);
    expect(await drawnStrokes(page), 'the edit changes strokes').not.toEqual(drawn);
    await revert(page);
    expect(await drawnStrokes(page), 'every part is back at the stroke Python drew').toEqual(drawn);
    await expectLiveEqualsGenerate(page, { label: name });
  });
}

// A Session 44 Result, which main saved without records of Python's paint,
// gets the records at Load (EU U2a), so Reset all strokes returns its feature
// parts and swatches to Python's strokes
// (tests/fixtures/sessions/forced-label-underlay-strokes.provenance.json).
test('Reset all strokes after loading a Session 44 Result saved without paint records matches Generate', async ({ page }) => {
  test.setTimeout(240_000);
  await openFresh(page);
  await loadSessionFile(page, 'tests/fixtures/sessions/forced-label-underlay-strokes.v44.gbdraw-session.json.gz');
  await appAction(page, 'resetAllStrokes');
  expect((await drawnStrokes(page)).some(([, , stroke]) => ['#2a9d8f', '#e63946'].includes(stroke)), 'no edited stroke is left')
    .toBe(false);
  await expectLiveEqualsGenerate(page, { label: 'Reset all strokes after Load' });
});

// OV-287: Reset Settings clears the stroke and Legend edits, and the displayed
// Result shows the reset draft: each feature part and swatch returns to
// Python's stroke, and the Legend shows Python's colors and its deleted row,
// which the Legend list shows again (OV-289).
// The draft starts at the default settings, so the second Reset changes no
// setting that only Generate draws.
test('Reset Settings after Legend and popup strokes, a Legend row color and a Legend delete shows what Generate draws', async ({ page }) => {
  test.setTimeout(180_000);
  await open(page, { mode: 'circular', results: 'single', reflow: 'off' });
  await resetSettings(page);
  await generate(page);
  const drawn = await drawnStrokes(page);
  await legendRowStrokeColor(page, 'CDS', '#e63946');
  await legendRowColor(page, 'CDS', '#7b2cbf');
  await popupEdit(page, 'FL2', { stroke: '#2a9d8f' });
  await deleteLegendRow(page, 'GC content');
  expect(await drawnStrokes(page), 'the edits change strokes').not.toEqual(drawn);
  await resetSettings(page);
  expect(await drawnStrokes(page), 'every part is back at the stroke Python drew').toEqual(drawn);
  const listed = () => page.evaluate(() => window.__GBDRAW_APP__.legendEntries.map((entry) => entry.caption));
  const listedLive = await listed();
  await expectLiveEqualsGenerate(page, { label: 'Reset Settings' });
  // OV-289: the deleted row is listed again, as Generate lists it.
  expect(listedLive, 'the Legend list after Reset').toEqual(await listed());
});
// EU U2a review #1: a Session older than 40 adopts no feature catalog. The
// editor port addresses the displayed features from the Session's extracted
// features and the fills it draws, so a Legend row stroke, its reset and a
// palette change reach the loaded 0.13.0 (Session 30) Result.
test('a stroke and a palette change reach a Session 30 Result loaded without a feature catalog', async ({ page }) => {
  test.setTimeout(240_000);
  await openFresh(page);
  await loadSessionFile(page, 'tests/fixtures/sessions/BGC0000708-BGC0000713.v30.gbdraw-session.json.gz');
  expect(await page.evaluate(async () => (await import('/gbdraw/web/js/state.js')).state.featureCatalog.value)).toBeNull();
  const drawn = await drawnStrokes(page);
  const row = await page.evaluate(() => window.__GBDRAW_APP__.legendEntries[0].caption);
  await legendRowStroke(page, row, '#e63946', 3);
  expect((await drawnStrokes(page)).some(([, , stroke]) => stroke === '#e63946'), 'the row stroke is shown').toBe(true);
  await appAction(page, 'resetAllStrokes');
  expect(await drawnStrokes(page), 'Reset all strokes returns Python\'s strokes').toEqual(drawn);
  const fills = () => page.evaluate(() => [...window.__GBDRAW_APP__.svgContainer.querySelectorAll('[data-gbdraw-feature-id]')]
    .map((element) => String(element.getAttribute('fill')).toLowerCase()));
  const before = await fills();
  await page.evaluate(() => {
    const app = window.__GBDRAW_APP__;
    app.setPaletteInstantPreview(true);
    // The fixture has a user default color, so the switch asks first (D-15);
    // take the palette's colors, as the dialog's "Change palette" does.
    return app.selectPalette(app.paletteNames.find((name) => name !== app.selectedPalette), 'palette');
  });
  await settleLive(page);
  // A 0.13.0 Result has no part or label records, so it is compared with
  // Generate by the fill of each feature. It also has no stroke-only outline
  // (fill none) on a partial feature, which Generate now draws.
  expect(await fills(), 'the palette change reaches the Result').not.toEqual(before);
  const featureFills = async () => Object.fromEntries(Object.entries((await semanticSnapshot(page)).features)
    .map(([id, { fill }]) => [id, fill.filter((value) => value !== 'none')]));
  const live = await featureFills();
  await generate(page);
  // OV-401: Generate writes the source-record hash into the IDs of the
  // reverse-complemented record 5, where the loaded 0.13.0 Result keeps the
  // drawn hash it was saved with. Each loaded ID is read as the ID Generate
  // gives the same feature: the catalog row's ID with its drawn hash.
  const generatedIds = await page.evaluate(async () => {
    const { state } = await import('/gbdraw/web/js/state.js');
    return Object.fromEntries((state.featureCatalog.value?.items || []).flatMap((item) => item.features)
      .filter((row) => row.drawnSelector?.hash)
      .map((row) => [
        row.svgId.replace(row.biologicalFeatureId.replace(/~\d+$/, ''), row.drawnSelector.hash), row.svgId
      ]));
  });
  expect(Object.fromEntries(Object.entries(live).map(([id, fill]) => [generatedIds[id] ?? id, fill])),
    'each feature is filled as Generate fills it').toEqual(await featureFills());
});

// The work guard: each row of WORK_ALLOWLIST (tests/web/helpers/live-work-allowlist.mjs),
// run in order on one two-Result batch page.
test('each edit kind runs only the compile stages and worker requests of its allowlist', async ({ page }) => {
  test.setTimeout(420_000);
  await open(page, { mode: 'circular', results: 'batch', reflow: 'off' });
  await addColorRule(page, { feat: 'CDS', qual: 'locus_tag', val: '_0001$', color: '#c83366', cap: 'CDS' });
  await generate(page);
  await page.evaluate(() => {
    window.__workLog = [];
    const hooks = window.__GBDRAW_TEST_HOOKS__ || {};
    const metric = hooks.onStructuralMetric;
    const render = hooks.beforeDiagramGenerationResponse;
    window.__GBDRAW_TEST_HOOKS__ = {
      ...hooks,
      onStructuralMetric: (event) => { metric?.(event); window.__workLog.push(event); },
      beforeDiagramGenerationResponse: (...args) => { window.__workLog.push({ name: 'render' }); return render?.(...args); }
    };
  });
  const observed = [];
  for (const { kind, run, before } of WORK_ALLOWLIST) {
    if (before) await before(page);
    await page.evaluate(() => { window.__workLog.length = 0; });
    await run(page);
    const log = await page.evaluate(() => window.__workLog.map(({ name, stages: ran, domains, operation }) => ({ name, ran, domains, operation })));
    // A rerender compiles Generate's plan (no domains) and shows its Result.
    const live = log.filter(({ name, domains }) => name === 'editorPlanCompile' && domains);
    const sent = log.flatMap(({ name, operation }) => (
      name === 'diagramHelperRequest' ? [operation] : (name === 'render' ? ['render'] : [])
    )).sort();
    observed.push({ kind, stages: [...new Set(live.flatMap(({ ran }) => ran))].sort(), compiles: live.length, requests: sent });
  }
  console.log(JSON.stringify(observed));
  expect(observed).toEqual(WORK_ALLOWLIST.map(({ kind, stages, compiles, requests }) => ({
    kind, stages: [...stages].sort(), compiles, requests: [...requests].sort()
  })));
});

// OV-200 (PD-OI-066): with Label Rendering = Embedded Only, Python draws no
// label that does not fit inside its feature. A Feature Edits TSV that turns
// such a label On, with Auto Reflow off, shows what Generate draws: no label.
// OV-237: a TSV label text reaches labels no popup has bound yet.
const LONG_LABEL = 'A_LABEL_TEXT_FAR_TOO_LONG_TO_FIT_INSIDE_ITS_FEATURE_'.repeat(4);
const TSV_LABEL_CASES = [
  { mode: 'circular', rendering: 'embedded_only', visibility: 'on', text: LONG_LABEL, name: 'Label visibility On for a label that does not fit, Embedded Only' },
  { mode: 'linear', rendering: 'embedded_only', visibility: 'on', text: LONG_LABEL, name: 'Label visibility On for a label that does not fit, Embedded Only' },
  { mode: 'circular', rendering: null, visibility: '', text: 'alpha edited', name: 'a label text' },
  // U2BFIX2 review M2: the table's label follow with Auto Reflow on.
  { mode: 'circular', rendering: 'embedded_only', visibility: '', text: LONG_LABEL, reflow: 'on', name: 'a label text that does not fit, Embedded Only, Auto Reflow on' }
];
for (const { mode, rendering, visibility, text: labelText, reflow = 'off', name } of TSV_LABEL_CASES) {
  test(`Load Feature Edits TSV with ${name} (${mode}, one Result, labels unbound)`, async ({ page }) => {
    test.setTimeout(180_000);
    await open(page, { mode, results: 'single', reflow });
    if (rendering) {
      await page.evaluate((value) => { window.__GBDRAW_APP__.adv.label_rendering = value; }, rendering);
      await generate(page);
    }
    await loadFeatureEdits(page, [['FL1', visibility, labelText]]);
    expect(await page.evaluate(() => Object.values(window.__GBDRAW_APP__.featureOverrides).some((row) => row.labelText)),
      'the table applied the row').toBe(true);
    await expectLiveEqualsGenerate(page, { label: `Feature Edits TSV, ${name} (${mode})` });
  });
}
// The same for a popup label text that no longer fits (Label visibility Default).
test('Popup label text that does not fit, Embedded Only (circular, one Result)', async ({ page }) => {
  test.setTimeout(180_000);
  await open(page, { mode: 'circular', results: 'single', reflow: 'off' });
  await page.evaluate(() => { window.__GBDRAW_APP__.adv.label_rendering = 'embedded_only'; });
  await generate(page);
  await popupEdit(page, 'FL1', { labelText: LONG_LABEL });
  await expectLiveEqualsGenerate(page, { label: 'popup label text that does not fit' });
});

// CD-05: a case reuses the warm page only when the restore brought it back to
// its fresh boot. The restore takes back inputs, settings and edits; what it
// cannot take back (a page global here) gives the next case a fresh page.
test.describe('the warm app page', () => {
  test.describe.configure({ mode: 'serial' });
  const warmPage = () => test.info().annotations.find(({ type }) => type === WARM_APP_PAGE)?.description;
  const freshBootView = (page) => page.evaluate(async () => {
    const { state } = await import('/gbdraw/web/js/state.js');
    const app = window.__GBDRAW_APP__;
    return {
      mode: app.mode,
      results: app.results.length,
      linearSources: state.linearSeqs.filter((seq) => seq.gb || seq.gff || seq.fasta).length,
      showLabelsLinear: state.drawingOf('linear').form.show_labels_linear,
      undo: window.__GBDRAW_HISTORY__.getUndoCount(),
      leftByCase: window.leftByCase ?? null
    };
  });

  test('a case leaves a Linear Result, a setting and a popup edit', async ({ page }) => {
    test.setTimeout(180_000);
    await open(page, { mode: 'linear', results: 'single', reflow: 'off' });
    await popupEdit(page, 'FL1', { fill: '#123456' });
    expect((await freshBootView(page)).results).toBe(1);
  });

  test('the next case reuses the page at its fresh boot', async ({ page }) => {
    expect(warmPage()).toBe('reused');
    await openFresh(page);
    expect(await freshBootView(page)).toEqual({
      mode: 'circular', results: 0, linearSources: 0, showLabelsLinear: 'none', undo: 0, leftByCase: null
    });
    await page.evaluate(() => { window.leftByCase = true; });
  });

  test('a case after one that left a page global gets a fresh page', async ({ page }) => {
    expect(warmPage()).toBe('fresh: page.globals.leftByCase: "<missing>" -> "boolean"');
    await openFresh(page);
    expect((await freshBootView(page)).leftByCase).toBe(null);
  });

  // OV-1001: the app asks before it downloads a Session file of a name it
  // downloaded before. The answer is the downloaded name, or the question.
  const saveSession = async (page) => {
    const downloaded = page.waitForEvent('download').then((file) => file.suggestedFilename());
    const question = page.locator('#choice-prompt-title');
    await page.getByRole('banner').getByRole('button', { name: 'Save Session', exact: true }).click();
    await answerSessionTitle(page, '');
    const asked = question.waitFor().then(() => question.textContent());
    return Promise.race([downloaded, asked]);
  };

  test('a case saves a Session', async ({ page }) => {
    await openFresh(page);
    expect(await saveSession(page)).toBe('gbdraw_session.json.gz');
  });

  test('the next case saves a Session of that name without the question', async ({ page }) => {
    await openFresh(page);
    expect(await saveSession(page)).toBe('gbdraw_session.json.gz');
    expect(warmPage()).toBe('fresh: savedSessions: [] -> ["gbdraw_session.json.gz"]');
  });
});
