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
const { evaluateWithRetainedPromise } = require('./helpers/app-lifecycle.cjs');
const { expectLiveEqualsGenerate, semanticSnapshot, settleLive, showResult } = require('./helpers/live-generate-parity.cjs');
const {
  open, generate, popupEdit, addVisibilityRule, appAction, addColorRule, history, FL1_OFF,
  legendRowColor, FL1_ALPHA, renameRow
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

// The Legend editor's stroke color control, one History step.
const legendRowStrokeColor = async (page, caption, color) => {
  await evaluateWithRetainedPromise(page, async ({ row, value }) => {
    const app = window.__GBDRAW_APP__;
    await app.setLegendEntryStrokeColorValue(app.legendEntries.findIndex((entry) => entry.caption === row), value);
  }, { row: caption, value: color });
  await settleLive(page);
};

// A row added in the Legend editor.
const legendRowAdd = async (page, caption, color) => {
  await evaluateWithRetainedPromise(page, async ({ row, value }) => {
    const app = window.__GBDRAW_APP__;
    app.newLegendCaption = row;
    app.newLegendColor = value;
    await app.addNewLegendEntry();
  }, { row: caption, value: color });
  await settleLive(page);
};

// The edit kinds. Each appears at least once in the matrix.
const KINDS = [
  'visibility rule add', 'visibility rule action', 'visibility rule delete', 'feature Off', 'feature On',
  'label text', 'label Off', 'label On', 'color rule add', 'color rule color', 'color rule delete', 'feature color',
  'undo', 'redo', 'Result switch', 'legend color', 'legend add', 'legend stroke', 'feature stroke',
  'feature legend name', 'feature color reset', 'label Default'
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
  // OV-121: a row added in the Legend editor takes the first row's stroke as
  // Generate copies it, not that row's stroke edit.
  {
    kind: 'legend add',
    edit: 'Legend editor row add after a stroke on the first row',
    states: { mode: 'linear', results: 'single', reflow: 'on', labels: 'bound' },
    setup: async (page) => { await legendRowStroke(page, 'CDS', '#e63946', 3); await generate(page); },
    run: (page) => legendRowAdd(page, 'Manual row', '#118833')
  },
  {
    kind: 'legend add',
    edit: 'Legend editor row add after a stroke on the first row',
    states: { mode: 'circular', results: 'single', reflow: 'off', labels: 'unbound' },
    setup: async (page) => { await legendRowStroke(page, 'CDS', '#e63946', 3); await generate(page); },
    run: (page) => legendRowAdd(page, 'Manual row', '#118833')
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
const deleteLegendRow = (page, caption) => evaluateWithRetainedPromise(page, async (row) => {
  const app = window.__GBDRAW_APP__;
  await app.deleteLegendEntry(app.legendEntries.findIndex((entry) => entry.caption === row));
}, caption).then(() => settleLive(page));
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
    // Generate does not draw a row the editor added; the deleted row stays gone.
    domain: 'a Legend row added in the editor and deleted',
    edit: async (page) => { await legendRowAdd(page, 'Manual row', '#118833'); await deleteLegendRow(page, 'Manual row'); },
    sameDrawing: true
  },
  {
    // The next two: Generate after Reset draws neither the added row nor the new caption.
    domain: 'a Legend row added in the editor',
    edit: (page) => legendRowAdd(page, 'Manual row', '#118833')
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
