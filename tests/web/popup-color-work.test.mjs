// OV-225 (R14-3, a work allowlist, not timings): from a popup color pick to
// the scope dialog, the app runs exactly the listed History stages and Worker
// requests. The dialog reads only the saved rules, so the rules the edit may
// add (the clicked feature's hash and label rules) are prepared when it
// commits. On a Session just loaded, that preparation started the Worker's
// Python runtime before the dialog opened (about 6 s on Vnig_TUMSAT-TG-2018).
import assert from 'node:assert/strict';
import { test } from 'node:test';
import { createRulePreparation, runWhenPrepared } from '../../gbdraw/web/js/app/rule-matching.js';
import { createFeatureColorActions } from '../../gbdraw/web/js/app/feature-editor/color-actions.js';
import { createHistoryManager } from '../../gbdraw/web/js/services/history.js';
import { createDialogChoice } from '../../gbdraw/web/js/app/history-inputs.js';
import { createResultsManager } from '../../gbdraw/web/js/app/results.js';
import { withDrawings } from './helpers/drawing-state.mjs';

const ref = (value) => ({ value });

const featuresOf = (products) => products.map((product, index) => ({
  id: `feature-${index}`, svg_id: `f${index}`, type: 'CDS', product, qualifiers: { product: [`p${index}`] }
}));
const savedRule = { feat: 'CDS', qual: 'product', val: '^p3$', color: '#111111', cap: 'p3' };

// Python's answer for `hash` and `product` rules (gbdraw/web_support/rule_matching.py).
const evaluateRules = async ({ features: payloads, rules, kind }) => {
  if (kind === 'color-captions') return { rules };
  const matches = payloads.map((payload) => rules.flatMap((rule, index) => {
    const values = rule.qual === 'hash' ? [payload.selector.hash] : payload.qualifiers[rule.qual] || [];
    return values.some((value) => new RegExp(rule.val).test(value)) ? [index] : [];
  }));
  return { matches, priorities: matches.map((row) => row.map(() => 0)) };
};

const setup = ({
  siblings, features = featuresOf(['p', 'p', 'p', 'p']), savedRules = [savedRule], commit = async () => true,
  userColors = {}, queued = false, appliedColors = { CDS: '#cccccc' }
}) => {
  /** @type {string[]} */
  const stages = [];
  /** @type {string[] | null} */
  let atDialog = null;
  /** @type {{ busy: boolean, pending: boolean }[]} */
  const opened = [];
  /** @type {{ pending: boolean }[]} */
  const closed = [];
  // A dialog records the stages run before it opened, whether it showed busy
  // when it opened, and whether History's step was still open when it closed.
  const dialog = (fields = {}) => ({
    ...fields,
    _show: false,
    get show() { return this._show; },
    set show(value) {
      if (value && !this._show) {
        atDialog = [...stages];
        opened.push({ busy: busy(), pending: history.mutationPending() });
      }
      if (!value && this._show) closed.push({ pending: history.mutationPending() });
      this._show = value;
    }
  });
  const featureStyleScopeDialog = dialog();
  const manualSpecificRules = [...savedRules];
  const state = withDrawings({
    extractedFeatures: ref(features), biologicalFeatures: ref(features), manualSpecificRules,
    svgResultIdentity: ref('one'), legendEntries: ref([{ caption: 'CDS', color: '#cccccc' }]),
    featureColorOverrides: {}, legendColorOverrides: {}, legendStrokeOverrides: {}, featureStrokeOverrides: {},
    addedLegendCaptions: ref(new Set()), results: ref([]), selectedResultIndex: ref(0), svgContainer: ref({ querySelector: () => null }),
    clickedFeature: ref({ feat: features[0], svg_id: features[0].svg_id, legendName: '' }), featureStyleScopeDialog,
    resetColorDialog: dialog(), legendRenameDialog: dialog(), originalLegendOrder: ref([]), originalLegendColors: ref({}),
    originalSvgStroke: ref({ color: null, width: null }), appliedPaletteColors: ref({ ...appliedColors }),
    skipCaptureBaseConfig: ref(false), skipExtractOnSvgChange: ref(false),
    // The palette members results.js reads: "default" is applied; "forest" is
    // queued while Instant Preview is off.
    paletteDefinitions: ref({ default: { CDS: '#cccccc' }, forest: { CDS: '#228b22' } }),
    paletteInstantPreviewEnabled: ref(false), appliedPaletteName: ref('default'),
    paletteColorsDialog: { show: false, kind: 'switch', fromPalette: '', toPalette: '', count: 0, keysText: '' },
    selectedPalette: ref(queued ? 'forest' : 'default'), pendingPaletteName: ref(queued ? 'forest' : ''),
    pendingPaletteColors: ref(queued ? { CDS: '#228b22' } : {}),
    currentColors: ref({ CDS: queued ? '#228b22' : '#cccccc', ...userColors }), hasPendingPaletteDraft: ref(queued)
  });
  const preparation = createRulePreparation({
    state,
    evaluate: (payload) => {
      stages.push(`worker:evaluateRules:${payload.kind}:${payload.rules.map((rule) => rule.qual).join('+')}`);
      return evaluateRules(payload);
    }
  });
  const history = createHistoryManager({
    buildIntent: () => {
      stages.push('history:buildIntent');
      return { rules: JSON.parse(JSON.stringify(manualSpecificRules)), colors: { ...state.currentColors.value } };
    },
    signatureFor: (value) => { stages.push('history:signature'); return JSON.stringify(value); },
    applyIntent: () => {}, buildCheckpoint: () => ({}), applyCheckpoint: () => {}
  });
  // app-setup.js: `createDialogChoice(...)`.
  const dialogChoice = createDialogChoice({ mutationPending: history.mutationPending, runUndoable: history.runUndoable, ref });
  // app-setup.js: the palette owner's ports (D-15).
  const palette = createResultsManager({ state, closeAfterDialogChoice: dialogChoice.closeAfterChoice });
  const actions = createFeatureColorActions({
    state, nextTick: async () => {}, onLegendGeometryChanged: () => {}, extractLegendEntries: () => {},
    getFeatureElements: () => [], getFeatureFillElements: () => [],
    closeAfterDialogChoice: dialogChoice.closeAfterChoice,
    readUserDefaultColor: palette.readUserDefaultColor,
    setDefaultColor: (drawing, key, color) => {
      stages.push('setDefaultColor');
      palette.setDefaultColor(drawing, key, color);
    },
    ruleActions: {
      runWithRuleMatches: (rules, commit) => runWhenPrepared(state, () => [preparation.prepare(rules)], commit),
      commitSpecificRules: async (rules) => {
        stages.push('commitSpecificRules');
        if (!await commit()) return false;
        manualSpecificRules.splice(0, manualSpecificRules.length, ...rules);
        return true;
      },
      countFeaturesMatchingRule: () => 0, findExistingColorForCaption: () => null,
      findFeaturesWithSameDisplayedLabel: () => [], findFeaturesWithSameIndividualLabel: () => [],
      findFeaturesWithSameLegendItem: (feature) => (siblings ? features.filter((other) => other !== feature) : []),
      findMatchingRegexRule: () => null, getDisplayedFeatureLabel: () => '',
      getIndividualFeatureLabel: (feature) => feature.qualifiers.product[0],
      effectiveLegendCaptions: () => () => 'CDS',
      getFeatureQualifier: (feature) => ({ qual: 'hash', val: feature.svg_id }),
      getLabelSpecificRule: (feature, label) => (label ? { feat: feature.type, qual: 'product', val: `^${label}$` } : null),
      getLegendRowRules: () => []
    }
  });
  // app-setup.js: `undoableAction(label, action)`.
  const pick = (color) => history.runUndoable('Change feature color', () => actions.updateClickedFeatureColor(color));
  const rename = (caption) => {
    state.clickedFeature.value.legendName = caption;
    return history.runUndoable('Rename legend item', () => actions.handleLegendNameCommit());
  };
  const reset = () => history.runUndoable('Reset feature color', () => actions.resetClickedFeatureFillColor());
  // app-setup.js: `scopeChoiceWithHistory(label, choice, cancel)`.
  const choose = (label, handler, cancel) => dialogChoice.withHistory(() => label, handler, cancel);
  // app-setup.js: `dialogChoicePending`, what a dialog's `busy` reads.
  const busy = () => dialogChoice.pending.value;
  const choices = {
    scope: choose('Change feature color', actions.handleFeatureStyleScopeChoice, actions.cancelFeatureStyleScope),
    rename: choose('Rename legend item', actions.handleLegendRenameChoice, actions.cancelLegendRename),
    reset: choose('Reset feature color', actions.handleResetColorChoice, actions.cancelResetColor)
  };
  return {
    stages, preparation, pick, rename, reset, choices, history, manualSpecificRules, dialogStages: () => atDialog,
    opened, closed, busy, state,
    featureStyleScopeDialog, legendRenameDialog: state.legendRenameDialog, resetColorDialog: state.resetColorDialog
  };
};

test('a popup color pick opens the scope dialog after the History capture alone', async () => {
  const { stages, preparation, pick, dialogStages, featureStyleScopeDialog } = setup({ siblings: true });
  assert.equal(await preparation.prepare([savedRule]), true);
  stages.length = 0;
  await pick('#123456');
  assert.equal(featureStyleScopeDialog.show, true);
  assert.deepEqual(dialogStages(), ['history:buildIntent', 'history:signature']);
});

test('a popup color pick without a dialog prepares the rules it commits', async () => {
  const { stages, preparation, pick, dialogStages } = setup({ siblings: false });
  assert.equal(await preparation.prepare([savedRule]), true);
  stages.length = 0;
  await pick('#123456');
  assert.equal(dialogStages(), null);
  assert.deepEqual(stages, [
    'history:buildIntent', 'history:signature', 'worker:evaluateRules:color:hash+product', 'commitSpecificRules',
    'history:buildIntent', 'history:signature'
  ]);
});

// The popup's Legend name and its fill reset open their dialogs from the saved
// rules too, and prepare the rules they may add when the choice commits.
for (const [name, run, dialogOf] of [
  ['a popup Legend rename', (setup_) => setup_.rename('Renamed'), (setup_) => setup_.legendRenameDialog],
  ['a popup fill reset', (setup_) => setup_.reset(), (setup_) => setup_.resetColorDialog]
]) {
  test(`${name} opens its dialog after the History capture alone`, async () => {
    const setup_ = setup({ siblings: true });
    assert.equal(await setup_.preparation.prepare([savedRule]), true);
    setup_.stages.length = 0;
    await run(setup_);
    assert.equal(dialogOf(setup_).show, true);
    assert.deepEqual(setup_.dialogStages(), ['history:buildIntent', 'history:signature']);
  });
}

// PD-OI-088 (OIC-028, D-12): a dialog opens with its choices ready, although
// it opens inside the step of the edit that asked for it. From a choice until
// its History step ends, the dialog stays open and busy; another choice or
// Cancel does nothing; the dialog closes once the step has ended.
const clickedHashRule = { feat: 'CDS', qual: 'hash', val: 'f0', color: '#222222', cap: 'p' };
for (const [name, open, dialogOf, choice, savedRules] of [
  ['scope', (setup_) => setup_.pick('#123456'), (setup_) => setup_.featureStyleScopeDialog, 'single', [savedRule]],
  ['rename', (setup_) => setup_.rename('Renamed'), (setup_) => setup_.legendRenameDialog, 'single', [savedRule]],
  ['reset', (setup_) => setup_.reset(), (setup_) => setup_.resetColorDialog, 'this', [savedRule, clickedHashRule]]
]) {
  test(`a popup ${name} dialog ignores a second choice and Cancel until its choice commits`, async () => {
    /** @type {(value: boolean) => void} */
    let release = () => {};
    let gated = false;
    const setup_ = setup({
      siblings: true, savedRules,
      commit: () => (gated ? new Promise((resolve) => { release = resolve; }) : Promise.resolve(true))
    });
    assert.equal(await setup_.preparation.prepare(savedRules), true);
    await open(setup_);
    assert.equal(dialogOf(setup_).show, true);
    const busyAtOpen = setup_.opened.map(({ busy }) => busy);
    const undoCount = setup_.history.getUndoCount();
    gated = true;
    setup_.stages.length = 0;
    const first = setup_.choices[name](choice);
    while (!setup_.stages.includes('commitSpecificRules')) await new Promise((resolve) => setImmediate(resolve));
    assert.equal(setup_.history.mutationPending(), true);
    const busyDuringChoice = setup_.busy();
    assert.equal(dialogOf(setup_).show, true);
    assert.equal(setup_.choices[name](choice), undefined);
    assert.equal(setup_.choices[name]('cancel'), undefined);
    assert.equal(dialogOf(setup_).show, true);
    release(true);
    await first;
    assert.equal(dialogOf(setup_).show, false);
    assert.deepEqual(
      { busyAtOpen, busyDuringChoice, busyAfter: setup_.busy(), stepOpenAtClose: setup_.closed.map(({ pending }) => pending) },
      { busyAtOpen: [false], busyDuringChoice: true, busyAfter: false, stepOpenAtClose: [false] }
    );
    assert.equal(setup_.history.mutationPending(), false);
    assert.equal(setup_.stages.filter((stage) => stage === 'commitSpecificRules').length, 1);
    assert.equal(setup_.history.getUndoCount(), undoCount + 1);
  });
}

// A choice whose commit fails leaves its dialog open and ready, so the user
// can choose again or cancel; it records no History step.
for (const [name, open, dialogOf, choice, savedRules] of [
  ['scope', (setup_) => setup_.pick('#123456'), (setup_) => setup_.featureStyleScopeDialog, 'single', [savedRule]],
  ['rename', (setup_) => setup_.rename('Renamed'), (setup_) => setup_.legendRenameDialog, 'single', [savedRule]],
  ['reset', (setup_) => setup_.reset(), (setup_) => setup_.resetColorDialog, 'this', [savedRule, clickedHashRule]]
]) {
  test(`a popup ${name} dialog stays open and ready when its choice fails`, async () => {
    let fail = false;
    const setup_ = setup({
      siblings: true, savedRules,
      commit: () => (fail ? Promise.reject(new Error('commit failed')) : Promise.resolve(true))
    });
    setup_.state.errorLog = { value: null };
    assert.equal(await setup_.preparation.prepare(savedRules), true);
    await open(setup_);
    const undoCount = setup_.history.getUndoCount();
    fail = true;
    await setup_.choices[name](choice);
    assert.equal(dialogOf(setup_).show, true);
    assert.deepEqual(setup_.closed, []);
    assert.equal(setup_.busy(), false);
    assert.equal(setup_.history.getUndoCount(), undoCount);
    assert.notEqual(setup_.state.errorLog.value, null);
  });
}

// D-15: "Apply to all" on a Legend row that no Specific color rule draws sets
// the feature type's default color as one History step; it writes no rule.
// Q1 B (Owner 2026-10-09): also while a palette is queued.
for (const queued of [false, true]) test(`Apply to all on a palette row sets the default color as one History step and commits no rule${queued ? ' (palette queued)' : ''}`, async () => {
  const setup_ = setup({ siblings: true, savedRules: [], queued });
  setup_.state.legendColorOverrides.CDS = '#999999';
  assert.equal(await setup_.preparation.prepare([]), true);
  await setup_.pick('#123456');
  assert.equal(setup_.featureStyleScopeDialog.defaultColorType, 'CDS');
  assert.equal(setup_.featureStyleScopeDialog.replacedDefaultColor, null);
  const undoCount = setup_.history.getUndoCount();
  setup_.stages.length = 0;
  await setup_.choices.scope('caption');
  assert.deepEqual(setup_.stages, [
    'history:buildIntent', 'history:signature', 'worker:evaluateRules:color:hash+product', 'setDefaultColor',
    'history:buildIntent', 'history:signature'
  ]);
  assert.equal(setup_.state.currentColors.value.CDS, '#123456');
  // The shown Result takes the color now; with a queued palette, the queued
  // colors take it too and the applied palette stays (Q1 B).
  assert.equal(setup_.state.appliedPaletteColors.value.CDS, '#123456');
  assert.equal(setup_.state.appliedPaletteName.value, 'default');
  assert.deepEqual(
    { pending: setup_.state.pendingPaletteName.value, queuedCDS: setup_.state.pendingPaletteColors.value.CDS },
    queued ? { pending: 'forest', queuedCDS: '#123456' } : { pending: '', queuedCDS: undefined }
  );
  assert.deepEqual(setup_.manualSpecificRules, []);
  assert.equal('CDS' in setup_.state.legendColorOverrides, false);
  assert.equal(setup_.featureStyleScopeDialog.show, false);
  assert.equal(setup_.history.getUndoCount(), undoCount + 1);
});

test('the scope dialog of a palette row names the user default color it replaces', async () => {
  const setup_ = setup({ siblings: true, savedRules: [], userColors: { CDS: '#aaaaaa' } });
  assert.equal(await setup_.preparation.prepare([]), true);
  await setup_.pick('#123456');
  assert.equal(setup_.featureStyleScopeDialog.defaultColorType, 'CDS');
  assert.equal(setup_.featureStyleScopeDialog.replacedDefaultColor, '#aaaaaa');
});

// A rule draws one of the row's features (p3): "Apply to all" writes
// per-feature rules, as before D-15.
test('Apply to all on a row a rule draws a feature of keeps writing rules', async () => {
  const setup_ = setup({ siblings: true, savedRules: [savedRule] });
  assert.equal(await setup_.preparation.prepare([savedRule]), true);
  await setup_.pick('#123456');
  assert.equal(setup_.featureStyleScopeDialog.defaultColorType, null);
  await setup_.choices.scope('caption');
  assert.equal(setup_.stages.includes('setDefaultColor'), false);
  assert.equal(setup_.stages.includes('commitSpecificRules'), true);
  assert.equal(setup_.state.currentColors.value.CDS, '#cccccc');
});

// Review 1: the choice does what the dialog showed. A rule draws one of the
// row's features when the dialog opens, so it shows no default-color line;
// when the saved rules change before the choice (here the rule goes away),
// "Apply to all" still writes rules instead of a default color.
test('Apply to all follows the dialog: no default-color line at open keeps the rule path', async () => {
  const setup_ = setup({ siblings: true, savedRules: [savedRule] });
  await setup_.pick('#123456');
  assert.equal(setup_.featureStyleScopeDialog.show, true);
  assert.equal(setup_.featureStyleScopeDialog.defaultColorType, null);
  setup_.manualSpecificRules.splice(0);
  setup_.stages.length = 0;
  await setup_.choices.scope('caption');
  assert.equal(setup_.stages.includes('setDefaultColor'), false);
  assert.equal(setup_.stages.includes('commitSpecificRules'), true);
  assert.equal(setup_.state.currentColors.value.CDS, '#cccccc');
});

// Review 3 (OV-262 in the popup): with the type's default color Auto, the
// popup's Reset fill color resets to the applied palette's color, as Generate
// draws it.
test('the popup Reset fill color after Auto resets to the applied palette color', async () => {
  const clickedRule = { feat: 'CDS', qual: 'hash', val: 'f0', color: '#222222', cap: 'p0' };
  const setup_ = setup({ siblings: true, savedRules: [clickedRule], appliedColors: { CDS: null } });
  setup_.state.currentColors.value = { CDS: null };
  setup_.state.svgContainer.value = { querySelector: () => null };
  assert.equal(await setup_.preparation.prepare([clickedRule]), true);
  await setup_.reset();
  assert.equal(setup_.resetColorDialog.show, true);
  setup_.stages.length = 0;
  await setup_.choices.reset('this_with_legend');
  assert.equal(setup_.stages.includes('commitSpecificRules'), true);
  assert.ok(setup_.manualSpecificRules.some((rule) => rule.qual === 'hash' && rule.val === 'f0' && rule.color === '#cccccc'),
    JSON.stringify(setup_.manualSpecificRules));
});
