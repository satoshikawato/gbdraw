// The steps the PD-OI-066 parity specs share: opening the fixtures in each
// state, Generate, popup edits, and the editor actions they drive
// (live-generate-parity.playwright.spec.js and its -legend and -track-data
// siblings). The comparison itself is tests/web/helpers/live-generate-parity.cjs.
const { expect } = require('@playwright/test');
const { readFileSync } = require('node:fs');
const { evaluateWithRetainedPromise, generateAndWaitForResult } = require('./app-lifecycle.cjs');
const { BATCH_FIXTURE, openFresh, openWithGenBank } = require('./audit-browser.cjs');
const { settleLive } = require('./live-generate-parity.cjs');
const { loadEditorLegendRows } = require('./mode-transition.cjs');

const SINGLE_FIXTURE = 'tests/fixtures/forced_label_underlay.gb';

// States: mode, results (one Result, or a two-Result Circular batch), Auto
// Reflow, and labels (bound: a feature popup was opened after the last
// Generate and before the edit; unbound: none was). A popup edit binds the
// labels itself.
const open = async (page, { mode, results, reflow }) => {
  if (mode === 'linear') {
    await openFresh(page);
    await page.getByRole('button', { name: 'Linear', exact: true }).click();
    await page.waitForFunction(() => window.__GBDRAW_APP__?.mode === 'linear');
    await page.evaluate(async (text) => {
      const app = window.__GBDRAW_APP__;
      app.setLinearSeqPrimaryFile(0, 'gb', new File([text], 'forced_label_underlay.gb', { type: 'text/plain', lastModified: 1000 }));
      await window.Vue.nextTick();
    }, readFileSync(SINGLE_FIXTURE, 'utf8'));
    await settleLive(page);
    await page.evaluate(() => { window.__GBDRAW_APP__.form.show_labels_linear = 'all'; });
  } else if (results === 'batch') {
    await openWithGenBank(page, BATCH_FIXTURE, () => {
      const app = window.__GBDRAW_APP__;
      app.form.labels_mode = 'out';
      app.form.multi_record_canvas = false;
      app.adv.circular_grouping_intent = 'batch';
    });
  } else {
    await openWithGenBank(page, SINGLE_FIXTURE, () => { window.__GBDRAW_APP__.form.labels_mode = 'out'; });
  }
  await page.evaluate((enabled) => { window.__GBDRAW_APP__.autoLabelReflowEnabled = enabled; }, reflow === 'on');
  await generate(page);
  expect(await page.evaluate(() => window.__GBDRAW_APP__.results.length)).toBe(results === 'batch' ? 2 : 1);
};

const generate = async (page) => {
  await generateAndWaitForResult(page);
  await settleLive(page);
};

// One popup edit through the actions its controls call. The popup opens from
// the Features list (the displayed Result's features, hidden ones included),
// as its Edit button does; `match` is a locus_tag or an object of feature
// fields. A scope dialog takes `scope` (default: this feature only). No case
// expects a Label visibility On dialog, so one fails the edit. An empty edit
// opens and closes the popup, as a reader looking at the feature does, which
// binds the labels.
const popupEdit = async (page, match, edit = {}) => {
  await evaluateWithRetainedPromise(page, async ({ target, change }) => {
    const app = window.__GBDRAW_APP__;
    const matches = (item) => (typeof target === 'string'
      ? item.locus_tag === target
      : Object.entries(target).every(([field, value]) => item[field] === value));
    await app.openFeatureEditorFromList(app.filteredFeatures.find(matches), null);
    await window.Vue.nextTick();
    if (!app.clickedFeature) throw new Error(`no popup for ${JSON.stringify(target)}`);
    if (change.fill) {
      await app.updateClickedFeatureColor(change.fill);
      if (app.featureStyleScopeDialog.show) await app.handleFeatureStyleScopeChoice(change.scope || 'single');
    }
    if (change.stroke) {
      await app.setClickedFeatureStrokeColorValue(change.stroke);
      if (app.featureStyleScopeDialog.show) await app.handleFeatureStyleScopeChoice(change.scope || 'single');
    }
    if (change.resetStroke) await app.resetClickedFeatureStroke();
    if (change.legendName !== undefined) {
      app.clickedFeature.legendName = change.legendName;
      await app.handleLegendNameCommit();
      if (app.legendRenameDialog.show && app.legendRenameDialog.mode === 'scope') {
        await app.handleLegendRenameChoice(change.scope || 'single');
      }
      if (app.legendRenameDialog.show) throw new Error(`unexpected Legend name dialog: ${app.legendRenameDialog.mode}`);
    }
    if (change.resetFill) {
      await app.resetClickedFeatureFillColor();
      if (app.resetColorDialog.show) await app.handleResetColorChoice(change.scope || 'this');
    }
    if (change.visibility) {
      app.clickedFeature.featureVisibility = change.visibility;
      await app.updateClickedFeatureVisibility(change.visibility);
      if (app.featureVisibilityScopeDialog.show) await app.handleFeatureVisibilityScopeChoice(change.scope || 'feature');
    }
    if (change.labelText !== undefined || change.labelVisibility) {
      if (change.labelText !== undefined) app.clickedFeature.labelText = change.labelText;
      if (change.labelVisibility) app.clickedFeature.labelVisibility = change.labelVisibility;
      const applied = Promise.resolve(app.updateClickedFeatureLabelText());
      const asked = new Promise((resolve) => {
        const poll = () => (app.labelOnDialog.show ? resolve('asked') : setTimeout(poll, 20));
        poll();
      });
      if (await Promise.race([applied.then(() => 'applied'), asked]) === 'asked') {
        app.handleLabelOnChoice('cancel');
        throw new Error(`unexpected Label visibility On dialog: ${app.labelOnDialog.reason}`);
      }
      if (app.labelTextScopeDialog.show) await app.handleLabelTextScopeChoice('single');
      if (app.hiddenLabelTextDialog.show) await app.handleHiddenLabelTextChoice('show');
    }
    app.clickedFeature = null;
  }, { target: match, change: edit });
  await settleLive(page);
};

const addVisibilityRule = async (page, fields) => {
  await evaluateWithRetainedPromise(page, async (ruleFields) => {
    const app = window.__GBDRAW_APP__;
    await app.addFeatureVisibilityRule();
    const index = app.featureVisibilityManualRules.length - 1;
    for (const [field, value] of Object.entries(ruleFields)) await app.setFeatureVisibilityRuleField(index, field, value);
  }, fields);
  await settleLive(page);
};

const appAction = async (page, name, ...args) => {
  await page.evaluate(({ action, values }) => window.__GBDRAW_APP__[action](...values), { action: name, values: args });
  await settleLive(page);
};

const addColorRule = async (page, rule) => {
  await page.evaluate((fields) => {
    const app = window.__GBDRAW_APP__;
    Object.assign(app.newSpecRule, fields);
    return app.addSpecificRule();
  }, rule);
  await settleLive(page);
};

const history = async (page, step) => {
  await page.evaluate((name) => window.__GBDRAW_HISTORY__[name](), step);
  await settleLive(page);
};

const FL1_OFF = { recordId: 'FORCEDLBL', featureType: 'CDS', qualifier: 'locus_tag', value: '^fl1$', action: 'off' };

// The Legend row of a type with one feature (the Legend editor's color control).
const legendRowColor = async (page, caption, color) => {
  await page.evaluate(({ row, value }) => {
    const app = window.__GBDRAW_APP__;
    return app.updateLegendEntryColor(app.legendEntries.findIndex((entry) => entry.caption === row), value);
  }, { row: caption, value: color });
  await settleLive(page);
};

const FL1_ALPHA = { feat: 'CDS', qual: 'locus_tag', val: '^FL1$', color: '#e63946', cap: 'alpha' };

// A Depth file of four windows on the record of SINGLE_FIXTURE.
const DEPTH_TSV = Array.from({ length: 4 }, (_, index) => `FORCEDLBL\t${index * 700 + 1}\t${10 + index}`).join('\n');

// Colors a Legend row through the Legend editor and checks that the row took it.
const colorLegendRow = async (page, caption, color = '#7b2cbf') => {
  expect(await page.evaluate(({ target, value }) => {
    const app = window.__GBDRAW_APP__;
    const index = app.legendEntries.findIndex((entry) => entry.caption === target);
    return index >= 0 && app.updateLegendEntryColor(index, value);
  }, { target: caption, value: color }), `Legend row ${caption}`).toBeTruthy();
  await settleLive(page);
};

// Opens one Result; a Circular diagram is generated again with the Multi-Record
// Canvas set to `canvas` unless it is null.
const openCanvas = async (page, mode, canvas) => {
  await open(page, { mode, results: 'single', reflow: 'off' });
  if (mode === 'circular' && canvas !== null) {
    await page.evaluate((value) => { window.__GBDRAW_APP__.form.multi_record_canvas = value; }, canvas);
    await generate(page);
  }
};

// Renames a Legend row through the Legend editor.
const renameRow = (page, from, to) => evaluateWithRetainedPromise(page, async ({ from, to }) => {
  const app = window.__GBDRAW_APP__;
  const index = app.legendEntries.findIndex((entry) => entry.caption === from);
  if (index < 0) throw new Error(`no Legend row "${from}": ${app.legendEntries.map((entry) => entry.caption).join(', ')}`);
  await app.renameLegendEntry(index, to);
}, { from, to });

// The Legend editor's stroke color control, one History step.
const legendRowStrokeColor = async (page, caption, color) => {
  await evaluateWithRetainedPromise(page, async ({ row, value }) => {
    const app = window.__GBDRAW_APP__;
    await app.setLegendEntryStrokeColorValue(app.legendEntries.findIndex((entry) => entry.caption === row), value);
  }, { row: caption, value: color });
  await settleLive(page);
};

// Legend editor rows ([caption, color]) from a Session that holds them
// (R15-2 retired Add legend item), loaded and generated.
const editorLegendRows = async (page, rows) => {
  await loadEditorLegendRows(page, rows);
  await settleLive(page);
};

// A palette change with instant preview (the palette menu).
const switchPalette = async (page, name) => {
  await page.evaluate((palette) => {
    const app = window.__GBDRAW_APP__;
    app.paletteInstantPreviewEnabled = true;
    return app.selectPalette(palette);
  }, name);
  await settleLive(page);
};

// A row deleted in the Legend editor, settled.
const deleteLegendRow = async (page, caption) => {
  await evaluateWithRetainedPromise(page, async (row) => {
    const app = window.__GBDRAW_APP__;
    await app.deleteLegendEntry(app.legendEntries.findIndex((entry) => entry.caption === row));
  }, caption);
  await settleLive(page);
};

module.exports = {
  SINGLE_FIXTURE,
  open,
  generate,
  popupEdit,
  addVisibilityRule,
  appAction,
  addColorRule,
  history,
  FL1_OFF,
  legendRowColor,
  FL1_ALPHA,
  DEPTH_TSV,
  colorLegendRow,
  openCanvas,
  renameRow,
  legendRowStrokeColor,
  editorLegendRows,
  switchPalette,
  deleteLegendRow
};
