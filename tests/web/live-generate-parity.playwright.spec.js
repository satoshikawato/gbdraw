// PD-OI-066 (LIVE-EDIT-EQUALS-REGENERATION; R1, R3): each live edit shows on
// the displayed Result what the next Generate from the same draft draws. The
// matrix pairs every edit kind with states of the editor; the check is
// expectLiveEqualsGenerate (tests/web/helpers/live-generate-parity.cjs), minus
// the differences tests/web/contracts/live-generate-parity-allowed.json allows.
// A case marked test.fail names the finding (OV-xx) of a confirmed mismatch;
// the PR that fixes it removes the mark.
const { test, expect } = require('@playwright/test');
const { readFileSync } = require('node:fs');
const { evaluateWithRetainedPromise, generateAndWaitForResult } = require('./helpers/app-lifecycle.cjs');
const { BATCH_FIXTURE, openFresh, openWithGenBank } = require('./helpers/audit-browser.cjs');
const { expectLiveEqualsGenerate, settleLive, showResult } = require('./helpers/live-generate-parity.cjs');

test.describe.configure({ retries: 0 });

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
const BATCH_0004_OFF = { recordId: '*', featureType: 'CDS', qualifier: 'locus_tag', value: '_0004$', action: 'off' };

// The edit kinds. Each appears at least once in the matrix.
const KINDS = [
  'visibility rule add', 'visibility rule action', 'visibility rule delete', 'feature Off', 'feature On',
  'label text', 'label Off', 'label On', 'color rule add', 'color rule color', 'feature color',
  'undo', 'redo', 'Result switch'
];

// The matrix: one edit kind in one set of states, with the setup before the
// edit. The states cover every pair of values (checked below); a Linear
// diagram has one Result.
const CASES = [
  {
    kind: 'visibility rule add',
    edit: 'Feature Visibility rule add',
    states: { mode: 'circular', results: 'single', reflow: 'off', labels: 'unbound' },
    run: (page) => addVisibilityRule(page, FL1_OFF),
    knownMismatch: 'OV-42: the legend keeps the order (and the rows) of the features drawn before the edit'
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
    run: (page) => addColorRule(page, { feat: 'CDS', qual: 'product', val: '^gtg start$', color: '#e63946', cap: 'GTG start' }),
    knownMismatch: 'OV-43: the legend keeps the CDS caption and appends the rule row last'
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
    run: (page) => popupEdit(page, 'FL1', { fill: '#c83366' }),
    knownMismatch: 'OV-43: the legend keeps the CDS caption and appends the rule row last'
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
    run: (page) => showResult(page, 1),
    knownMismatch: 'OV-43, OV-44: the other Result lacks the rule row and keeps the CDS caption'
  },
  {
    kind: 'Result switch',
    edit: 'Result switch after a Feature Visibility rule add',
    states: { mode: 'circular', results: 'batch', reflow: 'off', labels: 'unbound' },
    setup: (page) => addVisibilityRule(page, BATCH_0004_OFF),
    run: (page) => showResult(page, 1)
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
