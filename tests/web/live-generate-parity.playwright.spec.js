// PD-OI-066 (LIVE-EDIT-EQUALS-REGENERATION; R1, R3): each live edit shows on
// the displayed Result what the next Generate from the same draft draws. The
// matrix pairs every edit kind with states of the editor; the check is
// expectLiveEqualsGenerate (tests/web/helpers/live-generate-parity.cjs), minus
// the differences tests/web/contracts/live-generate-parity-allowed.json allows.
// A case marked test.fail names the finding (OV-xx) of a confirmed mismatch;
// the PR that fixes it removes the mark. The Legend cases are in
// live-generate-parity-legend.playwright.spec.js, the track-data and mode cases
// in live-generate-parity-track-data.playwright.spec.js, and the steps they all
// use in tests/web/helpers/live-generate-parity-steps.cjs.
const { test, expect } = require('@playwright/test');
const { readFileSync, writeFileSync } = require('node:fs');
const { gunzipSync } = require('node:zlib');
const { evaluateWithRetainedPromise } = require('./helpers/app-lifecycle.cjs');
const { BATCH_FIXTURE, loadSessionFile, openFresh } = require('./helpers/audit-browser.cjs');
const { expectLiveEqualsGenerate, semanticSnapshot, settleLive, showResult } = require('./helpers/live-generate-parity.cjs');
const { download, load, loadEditorLegendRows, switchMode } = require('./helpers/mode-transition.cjs');
const {
  open, generate, popupEdit, addVisibilityRule, appAction, addColorRule, history, FL1_OFF,
  legendRowColor, FL1_ALPHA, renameRow, legendRowStrokeColor, editorLegendRows, switchPalette, deleteLegendRow
} = require('./helpers/live-generate-parity-steps.cjs');

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

const BATCH_0004_OFF = { recordId: '*', featureType: 'CDS', qualifier: 'locus_tag', value: '_0004$', action: 'off' };
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

// A Legend row renamed in the Legend editor (`renameRow`), settled.
const renameLegendRow = async (page, caption, name) => {
  await renameRow(page, caption, name);
  await settleLive(page);
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
    kind: 'color rule add',
    edit: 'Specific color rule add',
    states: { mode: 'linear', results: 'single', reflow: 'off', labels: 'bound' },
    run: (page) => addColorRule(page, FL1_ALPHA)
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
const resetSettings = (page) => evaluateWithRetainedPromise(page, async () => { await window.__GBDRAW_APP__.resetSettings(); })
  .then(() => settleLive(page));
const legendList = (page) => page.evaluate(() => {
  const app = window.__GBDRAW_APP__;
  return app.legendEntries.map((entry) => ({ caption: entry.caption, color: app.legendEntryColor(entry) }));
});
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
    domain: 'a Legend row renamed',
    edit: async (page) => { await renameRow(page, 'CDS', 'Coding'); await settleLive(page); }
  },
  {
    // A popup "this feature only" color draws its own Legend row; Reset
    // removes its rule, and the rule owner asks for the rerender (OV-43).
    domain: 'a feature fill, this feature only, with its own Legend row (OV-290)',
    edit: (page) => popupEdit(page, 'FL1', { fill: '#2a9d8f' })
  },
  {
    domain: 'a feature hidden (OV-290)',
    edit: (page) => popupEdit(page, 'FL2', { visibility: 'off' })
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
for (const { domain, edit, undoRedo, sameDrawing } of RESET_CASES) {
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
    if (undoRedo) {
      const reset = await semanticSnapshot(page);
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
  await download(page, 'Save Session', saved);
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
  await download(page, 'Save Session', saved);
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
// The Save is held (its own `beforeExport`) until the Result it raced is mounted.
const displayWhileSaveHeld = async (page, resultIndex) => {
  await page.evaluate(async (index) => {
    const service = await import('/gbdraw/web/js/services/config.js');
    let release;
    const gate = new Promise((resolve) => { release = resolve; });
    window.releaseSave = release;
    window.__GBDRAW_APP__.selectResult(index);
    window.pendingSave = service.exportSession('declined-display', { beforeExport: () => gate });
  }, resultIndex);
  await expect.poll(() => page.evaluate(async (index) => {
    const { isCommittedSvgResultMounted } = await import('/gbdraw/web/js/services/svg-result-ingestion.js');
    const app = window.__GBDRAW_APP__;
    return app.sessionSavePending && app.selectedResultIndex === index && isCommittedSvgResultMounted(app.results[index]);
  }, resultIndex)).toBe(true);
  await settleLive(page);
  await evaluateWithRetainedPromise(page, async () => { window.releaseSave(); await window.pendingSave; });
  await settleLive(page);
};
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
const saveWhileDisplaying = (page, resultIndex) => evaluateWithRetainedPromise(page, async (index) => {
  const app = window.__GBDRAW_APP__;
  const saved = app.saveSessionWithTitle();
  app.selectResult(index);
  await saved;
}, resultIndex);
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
  await download(page, 'Save Session', saved);
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

// A Session 46 Result saved before the executor recorded Python's paint gets
// the records at Load (EU U2a), so Reset all strokes returns its feature
// parts and swatches to Python's strokes
// (tests/fixtures/sessions/forced-label-underlay-strokes.provenance.json).
test('Reset all strokes after loading a Session 46 Result saved without paint records matches Generate', async ({ page }) => {
  test.setTimeout(240_000);
  await openFresh(page);
  await loadSessionFile(page, 'tests/fixtures/sessions/forced-label-underlay-strokes.v46.gbdraw-session.json.gz');
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
// U2a review #3 (OV-195): Load reads the stroke Python drew on a connector
// from the connectors no stroke edit reached (Linear draws every record in one
// size class), not from the feature edit or the Session's block stroke, so
// Reset Stroke of a spliced feature on a Result saved before the executor's
// records returns its connector to the line stroke. The Session is saved here
// and its records are removed, as a Save before EU U1 wrote it.
test('Reset Stroke of a spliced feature after loading a Session saved without paint records matches Generate (linear, two records)', async ({ page, browser }, testInfo) => {
  test.setTimeout(240_000);
  await openFresh(page);
  await page.getByRole('button', { name: 'Linear', exact: true }).click();
  await page.waitForFunction(() => window.__GBDRAW_APP__?.mode === 'linear');
  await page.evaluate(async (text) => {
    window.__GBDRAW_APP__.setLinearSeqPrimaryFile(0, 'gb', new File([text], 'web_batch_two_records.gb', { type: 'text/plain', lastModified: 1000 }));
    await window.Vue.nextTick();
  }, readFileSync(BATCH_FIXTURE, 'utf8'));
  await settleLive(page);
  await generate(page);
  const drawn = await drawnStrokes(page);
  expect(drawn.filter(([, part]) => part === 'connector').length, 'both records draw a connector').toBe(2);
  await popupEdit(page, 'TESTA_0004', { stroke: '#e63946' });
  const saved = testInfo.outputPath('spliced-stroke.gbdraw-session.json');
  const bytes = await download(page, 'Save Session', saved);
  const session = JSON.parse((bytes[0] === 0x1f ? gunzipSync(bytes) : bytes).toString('utf8'));
  session.results.forEach((result) => { result.content = result.content.replace(/ data-gbdraw-base-[a-z-]+="[^"]*"/g, ''); });
  writeFileSync(saved, JSON.stringify(session));
  const loaded = await load(browser, saved);
  expect((await drawnStrokes(loaded)).some(([, , stroke]) => stroke === '#e63946'), 'the Result shows the saved stroke').toBe(true);
  await popupEdit(loaded, 'TESTA_0004', { resetStroke: true });
  expect(await drawnStrokes(loaded), 'every part is back at the stroke Python drew').toEqual(drawn);
  await expectLiveEqualsGenerate(loaded, { label: 'Reset Stroke after Load' });
  await loaded.context().close();
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
    app.paletteInstantPreviewEnabled = true;
    // The fixture has a user default color, so the switch asks first (D-15);
    // take the palette's colors, as the dialog's "Change palette" does.
    return app.selectPalette(app.paletteNames.find((name) => name !== app.selectedPalette), 'palette');
  });
  await settleLive(page);
  // A 0.13.0 Result has no part or label records, so it is compared with
  // Generate by the fill of each feature.
  expect(await fills(), 'the palette change reaches the Result').not.toEqual(before);
  const featureFills = async () => Object.fromEntries(Object.entries((await semanticSnapshot(page)).features)
    .map(([id, { fill }]) => [id, fill]));
  const live = await featureFills();
  await generate(page);
  expect(live, 'each feature is filled as Generate fills it').toEqual(await featureFills());
});

// Load Feature Edits TSV with one row per [locus_tag, label_visibility,
// label_text] of the first record.
const loadFeatureEdits = (page, rows) => evaluateWithRetainedPromise(page, async (edits) => {
  const app = window.__GBDRAW_APP__;
  const lines = edits.map(([tag, visibility, text]) => (
    `#1\thash=${app.extractedFeatures.find((item) => item.locus_tag === tag).svg_id}\t\t${visibility}\t${text}`
  ));
  const table = ['record\tfeature_selector\tfeature_visibility\tlabel_visibility\tlabel_text', ...lines, ''].join('\n');
  await app.loadFeatureEditTable({ target: { files: [new File([table], 'edits.tsv', { type: 'text/plain' })], value: '' } });
}, rows).then(() => settleLive(page));

// Label Rendering = Embedded Only, generated; the first Result is displayed.
// Records the first-record features whose labels it draws.
/** @type {string[]} */
let embeddedLabelFeatures = [];
const generateEmbeddedOnly = async (page) => {
  await page.evaluate(() => { window.__GBDRAW_APP__.adv.label_rendering = 'embedded_only'; });
  await generate(page);
  if (await page.evaluate(() => window.__GBDRAW_APP__.selectedResultIndex) !== 0) await showResult(page, 0);
  embeddedLabelFeatures = await page.evaluate(() => {
    const app = window.__GBDRAW_APP__;
    const svg = app.svgContainer.querySelector('svg');
    return app.extractedFeatures.filter((item) => String(item.locus_tag || '').startsWith('TESTA_')
      && svg.querySelector(`text[data-label-feature-id="${item.svg_id}"]`)).map((item) => item.locus_tag);
  });
  expect(embeddedLabelFeatures.length, 'Embedded Only draws two labels of the first record').toBeGreaterThan(1);
};

// D-15: the popup's Apply to all on the first feature of `type`, a palette row.
const applyToAllOnPaletteRow = (page, type, color) => evaluateWithRetainedPromise(page, async ({ featureType, value }) => {
  const app = window.__GBDRAW_APP__;
  await app.openFeatureEditorFromList(app.filteredFeatures.find((item) => item.type === featureType), null);
  await window.Vue.nextTick();
  await app.updateClickedFeatureColor(value);
  if (app.featureStyleScopeDialog.defaultColorType !== featureType) throw new Error(`${featureType} is no palette row`);
  await app.handleFeatureStyleScopeChoice('caption');
  app.clickedFeature = null;
}, { featureType: type, value: color }).then(() => settleLive(page));
// D-15: a switch to a palette not in `avoid`, answered with Keep my colors.
const switchPaletteKeeping = (page, avoid) => evaluateWithRetainedPromise(page, async (names) => {
  const app = window.__GBDRAW_APP__;
  app.selectPalette(app.paletteNames.find((name) => !names.includes(name)));
  if (!app.paletteColorsDialog.show) throw new Error('the palette switch did not ask');
  await app.handlePaletteColorsChoice('keep');
}, avoid).then(() => settleLive(page));

// The work guard (allowlist) at the app, on a two-Result batch: each edit kind
// runs exactly the compile stages of its entry (the compile's structural
// metric, live compiles only), a Result display or a History step compiles
// once, and an edit sends exactly the worker requests of its entry, each as
// often as listed (`render` is the automatic rerender). Stages are
// `COMPILE_STAGES` in app/candidate-render.js. Every path that shows edits
// on a Result has rows here: live edits, Undo, Redo, Result display, mode
// switch, Reset Settings, Session Load, and the rerender (`render`).
const WORK_ALLOWLIST = [
  {
    // A stroke action prepares only the saved rules (OV-198), whose matches
    // are known: no request.
    kind: 'feature stroke (popup)', stages: ['strokes'], compiles: 1, requests: [],
    run: (page) => evaluateWithRetainedPromise(page, async () => {
      const app = window.__GBDRAW_APP__;
      await app.setClickedFeatureStrokeColorValue('#2a9d8f');
      if (app.featureStyleScopeDialog.show) await app.handleFeatureStyleScopeChoice('single');
    }),
    before: (page) => evaluateWithRetainedPromise(page, async () => {
      const app = window.__GBDRAW_APP__;
      await app.openFeatureEditorFromList(app.filteredFeatures.find((item) => item.locus_tag === 'TESTA_0001'), null);
      await window.Vue.nextTick();
    }).then(() => settleLive(page))
  },
  { kind: 'Legend row stroke', stages: ['strokes'], compiles: 1, requests: [], run: (page) => legendRowStrokeColor(page, 'CDS', '#e63946') },
  { kind: 'Legend row color (no rule)', stages: ['legendFills'], compiles: 1, requests: [], run: (page) => legendRowColor(page, 'tRNA', '#7b2cbf') },
  // The step restores the editor's Legend rows (their colors among them): Legend rows and fills.
  { kind: 'History step (Undo of a Legend row color)', stages: ['legend', 'legendFills'], compiles: 1, requests: [], run: (page) => history(page, 'undo') },
  // The rules' matches are prepared: the palette sends no request.
  { kind: 'palette change', stages: ['fills', 'rules', 'legendFills'], compiles: 1, requests: [], run: (page) => switchPalette(page, 'arctic') },
  {
    // D-15: Apply to all on a palette row sets the type's default color and
    // adds no rule; the palette watcher shows it in one compile. The choice
    // first prepares the matches of the rules that could draw the row
    // (`withTargetRules`), one rule request.
    kind: 'Apply to all on a palette row', stages: ['fills', 'rules', 'legendFills'], compiles: 1, requests: ['evaluateRules'],
    run: (page) => applyToAllOnPaletteRow(page, 'tRNA', '#5e60ce')
  },
  {
    // D-15: a palette switch asks while a user default color exists; Keep is
    // one History step, which the palette watcher shows in one compile.
    kind: 'palette switch, Keep my colors', stages: ['fills', 'rules', 'legendFills'], compiles: 1, requests: [],
    run: (page) => switchPaletteKeeping(page, ['default', 'arctic'])
  },
  {
    kind: 'Result display after a palette change and a stroke', stages: ['legend', 'fills', 'rules', 'legendFills', 'strokes'],
    compiles: 1, requests: [], run: (page) => showResult(page, 1)
  },
  {
    // One compile for the field set that changes the rules Generate reads (the
    // value); the other steps of the add leave them as they were. The rule's
    // matches are one request; it hides features, so no rerender.
    kind: 'Feature visibility rule add', stages: ['visibility'], compiles: 1, requests: ['evaluateRules'],
    run: (page) => addVisibilityRule(page, BATCH_0004_OFF)
  },
  {
    // The rule's caption CDS gets a second color, so Python splits the CDS
    // row in two and the automatic rerender draws the Result again (OV-43);
    // its compile is Generate's plan, and the commit's own compile shows the
    // rows it adds. Requests: the commit's caption normalization and rule
    // matches; the rerender's rules are Python's captions already (OV-238).
    kind: 'color rule commit', stages: ['fills', 'rules', 'legendFills', 'legend'], compiles: 1,
    requests: ['evaluateRules', 'evaluateRules', 'render'],
    run: (page) => addColorRule(page, { feat: 'CDS', qual: 'locus_tag', val: '_0002$', color: '#2266aa', cap: 'CDS' })
  },
  {
    kind: 'Result display, nothing changed', stages: ['legend'], compiles: 1, requests: [],
    before: (page) => showResult(page, 0), run: (page) => showResult(page, 1)
  },
  {
    // U2BFIX3 review L-A: the display a Save declined shows the Legend; the
    // fills follow in one more compile once the Save ends. The color rule
    // commit's rerender left the rule matches to prepare: the display's
    // preparation is discarded while the Save runs, so the repaint prepares
    // them again.
    kind: 'Result display declined by Save, shown when Save ends', stages: ['legend', 'fills', 'rules', 'legendFills'],
    compiles: 2, requests: ['evaluateRules', 'evaluateRules'],
    before: async (page) => { await showResult(page, 0); await legendRowColor(page, 'tRNA', '#2a9d8f'); },
    run: (page) => displayWhileSaveHeld(page, 1)
  },
  {
    // U2BFIX3 review M-A: Save shows the fills the declined display left
    // behind before it writes the Result (the matches are prepared).
    kind: 'Save while the displayed Result lacks paint', stages: ['legend', 'fills', 'rules', 'legendFills'],
    compiles: 2, requests: [],
    before: async (page) => { await showResult(page, 0); await legendRowColor(page, 'tRNA', '#e9c46a'); },
    run: async (page) => {
      await Promise.all([page.waitForEvent('download'), saveWhileDisplaying(page, 1)]);
      await settleLive(page);
    }
  },
  {
    // OV-200 (U2BFIX2 review M1): only Python decides whether the label text
    // fits, so the table asks for the rerender; its one label follow places
    // the labels too, so one reflow runs (Auto Reflow on). Requests: the
    // table and the visibility rule's matches on the catalog Generate drew;
    // the rerender's rules are Python's captions already (OV-238).
    kind: 'Feature Edits TSV load (Embedded Only)', stages: ['visibility'], compiles: 1,
    requests: ['readFeatureOverrideTable', 'evaluateRules', 'render'],
    before: async (page) => {
      await generateEmbeddedOnly(page);
      await page.evaluate(() => { window.__GBDRAW_APP__.autoLabelReflowEnabled = true; });
    },
    run: (page) => loadFeatureEdits(page, [[embeddedLabelFeatures[0], '', 'dup edited']])
  },
  {
    // A label a table hides needs no fit decision: no rerender (Auto Reflow
    // off). The table keeps the label text it loaded before. Requests: the
    // table, and the color and visibility rule matches on the catalog the
    // rerender drew.
    kind: 'Feature Edits TSV load, label Off (Embedded Only)', stages: ['visibility'], compiles: 1,
    requests: ['readFeatureOverrideTable', 'evaluateRules', 'evaluateRules'],
    before: (page) => page.evaluate(() => { window.__GBDRAW_APP__.autoLabelReflowEnabled = false; }),
    run: (page) => loadFeatureEdits(page, [[embeddedLabelFeatures[0], '', 'dup edited'], [embeddedLabelFeatures[1], 'off', '']])
  },
  {
    // The rerender decides the fit; its rules are Python's captions already (OV-238).
    kind: 'popup label text (Embedded Only)', stages: [], compiles: 0, requests: ['render'],
    run: (page) => popupEdit(page, embeddedLabelFeatures[0], { labelText: 'dup' })
  },
  // U3a (R14-8): the Legend edits. A row's edits show live in one compile; a
  // rule change asks for the automatic rerender only when the displayed
  // Result cannot show the rows it regroups (a split, or a row another batch
  // Result draws). Captions with one color each send no caption request. An
  // edit after a rerender or a History restore first matches the rules on the
  // features of the Result it drew (one `evaluateRules`).
  {
    kind: 'Legend rename, row without features', stages: ['legend', 'legendFills', 'strokes'], compiles: 1, requests: ['evaluateRules'],
    run: (page) => renameLegendRow(page, 'GC content', 'GC %')
  },
  {
    // The rule takes one feature of the "other proteins" row of Result 1: a split.
    kind: 'color rule commit that splits a row', stages: ['fills', 'rules', 'legendFills', 'legend'], compiles: 1,
    requests: ['evaluateRules', 'render'],
    run: (page) => addColorRule(page, { feat: 'CDS', qual: 'locus_tag', val: '^TESTA_0005$', color: '#2a9d8f', cap: 'Five' })
  },
  {
    kind: 'Legend rename, rule row one Result draws (whole row)', stages: ['fills', 'rules', 'legendFills', 'legend'],
    compiles: 1, requests: ['evaluateRules'], run: (page) => renameLegendRow(page, 'Five', 'Fifth')
  },
  // The rename's History step holds the Result before and after it: Undo and
  // Redo restore the Result that shows the rows, so neither compiles nor asks.
  { kind: 'History step (Undo of the rule row rename)', stages: [], compiles: 0, requests: [], run: (page) => history(page, 'undo') },
  { kind: 'History step (Redo of the rule row rename)', stages: [], compiles: 0, requests: [], run: (page) => history(page, 'redo') },
  {
    kind: 'Legend row color, rule row', stages: ['fills', 'rules', 'legendFills'], compiles: 1, requests: ['evaluateRules'],
    run: (page) => legendRowColor(page, 'Fifth', '#264653')
  },
  {
    kind: 'Rules-panel color of a captioned rule', stages: ['fills', 'rules', 'legendFills'], compiles: 1, requests: [],
    run: (page) => page.evaluate(() => {
      const app = window.__GBDRAW_APP__;
      return app.setSpecificRuleField(app.manualSpecificRules.findIndex((rule) => rule.cap === 'Fifth'), 'color', '#e76f51');
    }).then(() => settleLive(page))
  },
  {
    // Result 2 draws the row too: Python draws it again there.
    kind: 'Legend rename, rule row both Results draw', stages: ['fills', 'rules', 'legendFills', 'legend'], compiles: 1,
    requests: ['render'], run: (page) => renameLegendRow(page, 'CDS [#2266aa]', 'Two')
  },
  {
    // P-1 (not decided): one rule per feature of the row, so their matches are
    // a request. The two records' tRNAs have one hash, so the rule also
    // renames Result 2's row, which Python draws again.
    kind: 'Legend rename, palette row with features', stages: ['fills', 'rules', 'legendFills', 'legend'], compiles: 1,
    requests: ['evaluateRules', 'evaluateRules', 'render'], run: (page) => renameLegendRow(page, 'tRNA', 'transfer RNA')
  },
  // A delete, its Undo and a Restore also show the strokes: a deleted row's
  // stroke leaves its features (OV-293).
  { kind: 'Legend delete', stages: ['legend', 'strokes'], compiles: 1, requests: [], run: (page) => deleteLegendRow(page, 'GC skew (+)') },
  {
    kind: 'History step (Undo of a Legend delete)', stages: ['legend', 'legendFills', 'strokes'], compiles: 1, requests: [],
    run: (page) => history(page, 'undo')
  },
  {
    kind: 'Legend Restore', stages: ['legend', 'legendFills', 'strokes'], compiles: 1, requests: [],
    before: (page) => deleteLegendRow(page, 'GC skew (+)'),
    run: (page) => page.evaluate(() => window.__GBDRAW_APP__.restoreAllDeletedLegendEntries()).then(() => settleLive(page))
  },
  {
    kind: 'Legend sort', stages: ['legend'], compiles: 1, requests: [],
    run: (page) => page.evaluate(() => window.__GBDRAW_APP__.sortLegendEntries('desc')).then(() => settleLive(page))
  },
  {
    // U3b review L2: a mode switch arrives at the Result its mode drew, which
    // already shows that drawing's intent (it was left unchanged), so the
    // display compiles nothing. OV-286: the arriving Result was applied with
    // another palette than the departing one (arctic, default), so the palette
    // watcher prepares the rules once more (one `evaluateRules`).
    kind: 'mode switch back to a Result', stages: [], compiles: 0, requests: ['evaluateRules'],
    before: async (page) => {
      await switchMode(page, 'linear');
      await page.evaluate(async (text) => {
        window.__GBDRAW_APP__.setLinearSeqPrimaryFile(0, 'gb', new File([text], 'batch.gb', { type: 'text/plain', lastModified: 1000 }));
        await window.Vue.nextTick();
      }, readFileSync(BATCH_FIXTURE, 'utf8'));
      await settleLive(page);
      await generate(page);
    },
    run: (page) => switchMode(page, 'circular').then(() => settleLive(page))
  },
  {
    // A Session this page saved, loaded into it: its Results already show
    // their drawing's intent, so the Load compiles nothing and asks nothing.
    kind: 'Session Load', stages: [], compiles: 0, requests: [],
    before: (page) => download(page, 'Save Session', test.info().outputPath('work-allowlist.gbdraw-session.json')),
    run: async (page) => {
      const shown = () => page.evaluate(async () => {
        const { getCommittedSvgResultRuntimeIdentity } = await import('/gbdraw/web/js/services/svg-result-ingestion.js');
        const app = window.__GBDRAW_APP__;
        return app.sessionImportPending ? null : getCommittedSvgResultRuntimeIdentity(app.results[app.selectedResultIndex]);
      });
      const before = await shown();
      await page.locator('input[accept^=".json,"]').setInputFiles(test.info().outputPath('work-allowlist.gbdraw-session.json'));
      await expect.poll(async () => { const now = await shown(); return now !== null && now !== before; }, { timeout: 180_000 }).toBe(true);
      await settleLive(page);
    }
  },
  {
    // R15-3 (OV-285): a rename onto a deleted row's caption asks; its Merge is
    // the Restore (which shows the row) and then the merge's rule commit, in
    // one History checkpoint. Each shows once: two compiles. The merge's
    // rule commit prepares its candidate rules (one request); the rows stay
    // on the displayed Result, so no rerender. It runs after the Load, so it
    // leaves the rule preparation of the rows above as they were (OV-286).
    kind: 'Legend rename onto a deleted row (Restore and merge)', stages: ['fills', 'legend', 'legendFills', 'rules', 'strokes'],
    compiles: 2, requests: ['evaluateRules'],
    before: (page) => deleteLegendRow(page, 'Two'),
    run: async (page) => {
      await renameRow(page, 'Fifth', 'Two');
      expect(await page.evaluate(() => window.__GBDRAW_APP__.legendRenameDialog.deletedTargetKey), 'the rename asks').not.toBe('');
      await evaluateWithRetainedPromise(page, () => window.__GBDRAW_APP__.handleLegendRenameChoice('merge'));
      await settleLive(page);
    }
  },
  {
    // Reset Settings as it is today (both drawings, one History checkpoint;
    // the drawing-scoped Reset of W3 replaces it). It shows the strokes, the
    // Legend rows and colors (OV-287) and the feature visibility (OV-290) it
    // clears in one compile; the palette watcher shows the default palette's
    // fills in a compile of its own, so a Reset that changes the palette
    // compiles twice. With nothing to reflow it sends no request.
    kind: 'Reset Settings', stages: ['fills', 'legend', 'legendFills', 'strokes', 'visibility'], compiles: 2, requests: [],
    run: resetSettings
  }
];

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
