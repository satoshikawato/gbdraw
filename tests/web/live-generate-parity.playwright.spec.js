// PD-OI-066 (LIVE-EDIT-EQUALS-REGENERATION; R1, R3): each live edit shows on
// the displayed Result what the next Generate from the same draft draws. The
// matrix pairs every edit kind with states of the editor; the check is
// expectLiveEqualsGenerate (tests/web/helpers/live-generate-parity.cjs), minus
// the differences tests/web/contracts/live-generate-parity-allowed.json allows.
// A case marked test.fail names the finding (OV-xx) of a confirmed mismatch;
// the PR that fixes it removes the mark. The cases after the matrix check what
// the Legend fixes of OV-42 to OV-44 ask of the automatic rerender.
const { test, expect } = require('@playwright/test');
const { readFileSync } = require('node:fs');
const { evaluateWithRetainedPromise, generateAndWaitForResult, reveal } = require('./helpers/app-lifecycle.cjs');
const { BATCH_FIXTURE, openFresh, openWithGenBank } = require('./helpers/audit-browser.cjs');
const { expectLiveEqualsGenerate, semanticSnapshot, settleLive, showResult } = require('./helpers/live-generate-parity.cjs');

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
const BATCH_CDS_OFF = { recordId: '*', featureType: 'CDS', qualifier: 'locus_tag', value: '.', action: 'off' };
// The Legend row of a type with one feature (the Legend editor's color control).
const legendRowColor = async (page, caption, color) => {
  await page.evaluate(({ row, value }) => {
    const app = window.__GBDRAW_APP__;
    return app.updateLegendEntryColor(app.legendEntries.findIndex((entry) => entry.caption === row), value);
  }, { row: caption, value: color });
  await settleLive(page);
};

const FL1_ALPHA = { feat: 'CDS', qual: 'locus_tag', val: '^FL1$', color: '#e63946', cap: 'alpha' };

// The edit kinds. Each appears at least once in the matrix.
const KINDS = [
  'visibility rule add', 'visibility rule action', 'visibility rule delete', 'feature Off', 'feature On',
  'label text', 'label Off', 'label On', 'color rule add', 'color rule color', 'color rule delete', 'feature color',
  'undo', 'redo', 'Result switch', 'legend color'
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

// OV-66: Undo of a Depth source removal brings back the Depth track and its
// tick text, and Redo hides them again; the live Result equals what Generate
// draws. A Depth source History step restores files, which suppresses the
// track-visibility watcher, so the step projects the visibility itself.
const DEPTH_TSV = Array.from({ length: 4 }, (_, index) => `FORCEDLBL\t${index * 700 + 1}\t${10 + index}`).join('\n');
const DEPTH_CASES = {
  'circular, clearing the Depth file': {
    mode: 'circular',
    add: (app, file) => app.setCircularDepthFile(0, file),
    remove: (app) => app.setCircularDepthFile(0, null)
  },
  'circular, removing the Depth track': {
    mode: 'circular',
    add: (app, file) => app.setCircularDepthFile(0, file),
    remove: (app) => app.removeCircularDepthTrack(0)
  },
  'linear, clearing the Depth file': {
    mode: 'linear',
    add: (app, file) => app.setLinearDepthFile(app.linearSeqs[0], 0, file),
    remove: (app) => app.setLinearDepthFile(app.linearSeqs[0], 0, null)
  }
};
for (const [name, { mode, add, remove }] of Object.entries(DEPTH_CASES)) {
  test(`Undo and Redo of a Depth source removal match Generate (${name})`, async ({ page }) => {
    test.setTimeout(180_000);
    await open(page, { mode, results: 'single', reflow: 'off' });
    const inStep = async (label, change, ...args) => {
      await page.evaluate(async ({ stepLabel, source, values }) => {
        const run = new Function('app', 'text', `return (${source})(app, text && new File([text], 'depth.tsv', { type: 'text/tab-separated-values' }));`);
        await window.__GBDRAW_HISTORY__.runUndoable(stepLabel, () => run(window.__GBDRAW_APP__, values[0]));
      }, { stepLabel: label, source: change.toString(), values: args });
      await settleLive(page);
    };
    await inStep('Change uploaded file', add, DEPTH_TSV);
    await generate(page);
    await inStep('Remove Depth', remove, null);
    await history(page, 'undo');
    await expectLiveEqualsGenerate(page, { label: `${name}: Undo` });
    await history(page, 'redo');
    await expectLiveEqualsGenerate(page, { label: `${name}: Redo` });
  });
}

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
    state.legendColorOverrides.Ghost = '#123456';
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
const colorLegendRow = async (page, caption, color = '#7b2cbf') => {
  expect(await page.evaluate(({ target, value }) => {
    const app = window.__GBDRAW_APP__;
    const index = app.legendEntries.findIndex((entry) => entry.caption === target);
    return index >= 0 && app.updateLegendEntryColor(index, value);
  }, { target: caption, value: color }), `Legend row ${caption}`).toBeTruthy();
  await settleLive(page);
};

const openCanvas = async (page, mode, canvas) => {
  await open(page, { mode, results: 'single', reflow: 'off' });
  if (mode === 'circular' && canvas !== null) {
    await page.evaluate((value) => { window.__GBDRAW_APP__.form.multi_record_canvas = value; }, canvas);
    await generate(page);
  }
};

// OV-65: a Legend color on a row named only by a track's data (an annotation set,
// a depth file) follows the caption: a data change retires the styles of the
// captions the data no longer names, in the History step of the change, so
// Undo brings back the data and the style (RETIRING_CASES, further below).
// A region annotation with a legend label draws a Legend row from its set; the
// slot of the set is added through the track slot control.
const addAnnotationRow = async (page, label) => {
  await page.evaluate((legendLabel) => {
    const app = window.__GBDRAW_APP__;
    const set = app.addAnnotationSet('regions');
    const annotation = app.addCoordinateAnnotation(set, { start: 100, end: 400 });
    annotation.legendLabel = legendLabel;
  }, label);
  await settleLive(page);
  await generate(page);
};

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
for (const [mode, addDepth] of [
  ['circular', (text) => { window.__GBDRAW_APP__.setCircularDepthFile(0, new File([text], 'depth.tsv', { type: 'text/tab-separated-values' })); }],
  ['linear', (text) => { const app = window.__GBDRAW_APP__; app.setLinearDepthFile(app.linearSeqs[0], 0, new File([text], 'depth.tsv', { type: 'text/tab-separated-values' })); }]
]) {
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
}

// OV-62 (PD-OI-061 amended): two Legend rows merge only when both draw features
// of one same type. Anything else offers Suffix and Cancel only.
const renameRow = (page, from, to) => evaluateWithRetainedPromise(page, async ({ from, to }) => {
  const app = window.__GBDRAW_APP__;
  const index = app.legendEntries.findIndex((entry) => entry.caption === from);
  if (index < 0) throw new Error(`no Legend row "${from}": ${app.legendEntries.map((entry) => entry.caption).join(', ')}`);
  await app.renameLegendEntry(index, to);
}, { from, to });
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

// OV-65: Legend styles follow the captions that track data names. Each data
// change below runs in one History step, as its control's does.
const inHistoryStep = (page, label, body, arg) => page.evaluate(
  async ({ stepLabel, source, value }) => {
    const change = new Function('app', 'value', `return (${source})(app, value);`);
    await window.__GBDRAW_HISTORY__.runUndoable(stepLabel, () => change(window.__GBDRAW_APP__, value));
  },
  { stepLabel: label, source: body.toString(), value: arg }
);

const legendStyleOf = (page, caption) => page.evaluate(async (target) => {
  const { state } = await import('/gbdraw/web/js/state.js');
  return {
    color: state.legendColorOverrides[target] ?? null,
    stroke: state.legendStrokeOverrides[target] ?? null
  };
}, caption);

const undo = async (page) => {
  await page.evaluate(() => window.__GBDRAW_APP__.undoHistory());
  await settleLive(page);
};

// Types a new text into a label field of an opened section, and leaves the
// field, as the reader does; the History step commits when the field loses focus.
const typeIntoLabel = async (page, section, label, text) => {
  await page.evaluate((name) => {
    document.querySelectorAll('details > summary').forEach((summary) => {
      if (summary.textContent.trim().startsWith(name) || summary.getAttribute('aria-label') === name) {
        summary.parentElement.open = true;
      }
    });
  }, section);
  const field = page.getByLabel(label, { exact: true });
  await field.fill(text);
  await field.blur();
  await settleLive(page);
};

const addDepthFile = async (page, name, mode = 'circular') => {
  await inHistoryStep(page, 'Change uploaded file', (app, { text, fileName, linear }) => {
    const file = new File([text], fileName, { type: 'text/tab-separated-values' });
    if (linear) app.setLinearDepthFile(app.linearSeqs[0], 0, file);
    else app.setCircularDepthFile(0, file);
  }, { text: DEPTH_TSV, fileName: name, linear: mode === 'linear' });
  await settleLive(page);
};
const depthCaption = (page) => page.evaluate(() => (
  window.__GBDRAW_APP__.legendEntries.find((entry) => /depth/i.test(entry.caption))?.caption ?? null
));
// The setup of a Depth case: a file `depth.tsv` drawn once, with a Legend color on its row.
const colorDepthRow = (mode = 'circular') => async (page) => {
  await addDepthFile(page, 'depth.tsv', mode);
  await generate(page);
  await colorLegendRow(page, await depthCaption(page));
};
const depthFileName = (page) => page.evaluate(() => {
  const app = window.__GBDRAW_APP__;
  return (app.mode === 'linear' ? app.linearSeqs[0]?.depth?.[0] : app.files.c_depth?.[0]?.[0])?.name ?? null;
});

const ANNOTATION_TSV = (legendLabel) => [
  'set_id\tid\tmark\tstart\tend\tlegend_label',
  `regions\tregion_1\thighlight\t100\t400\t${legendLabel}`
].join('\n');

const drawnLegendCaptions = async (page) => (await semanticSnapshot(page)).legend.map(({ caption }) => caption);
const NO_STYLE = { color: null, stroke: null };
const COLORED = '#7b2cbf';

// Each case runs a data change that removes a caption's rows as one History
// step. The style is retired in that step: Generate succeeds and draws no row
// of the data; Undo brings back the data and the style, and the live Result
// equals Generate.
const RETIRING_CASES = [
  {
    name: 'removing an annotation set',
    caption: async () => 'Region X',
    setup: async (page) => {
      await addAnnotationRow(page, 'Region X');
      await colorLegendRow(page, 'Region X');
    },
    change: (page) => inHistoryStep(page, 'Delete set', (app) => app.removeAnnotationSet(app.annotationSets[0])),
    restored: (page) => page.evaluate(() => window.__GBDRAW_APP__.annotationSets.length === 1)
  },
  {
    name: 'removing the depth file',
    caption: async (page) => depthCaption(page),
    setup: colorDepthRow(),
    change: (page) => inHistoryStep(page, 'Change uploaded file', (app) => app.setCircularDepthFile(0, null)),
    restored: async (page) => await depthFileName(page) === 'depth.tsv'
  },
  // A Depth track removal and a file of another name reach the same rule.
  {
    name: 'removing the Depth track',
    caption: async (page) => depthCaption(page),
    setup: colorDepthRow(),
    change: (page) => inHistoryStep(page, 'Remove Depth', (app) => app.removeCircularDepthTrack(0)),
    restored: async (page) => await depthFileName(page) === 'depth.tsv'
  },
  {
    name: 'replacing the depth file with a file of another name',
    caption: async (page) => depthCaption(page),
    setup: colorDepthRow(),
    change: (page) => addDepthFile(page, 'coverage.tsv'),
    restored: async (page) => await depthFileName(page) === 'depth.tsv',
    drawn: ['coverage']
  },
  {
    name: 'removing the Depth track in Linear',
    mode: 'linear',
    caption: async (page) => depthCaption(page),
    setup: colorDepthRow('linear'),
    change: (page) => inHistoryStep(page, 'Remove Depth', (app) => app.removeLinearDepthTrack(0)),
    restored: async (page) => await depthFileName(page) === 'depth.tsv'
  },
  {
    name: 'replacing the depth file with a file of another name in Linear',
    mode: 'linear',
    caption: async (page) => depthCaption(page),
    setup: colorDepthRow('linear'),
    change: (page) => addDepthFile(page, 'coverage.tsv', 'linear'),
    restored: async (page) => await depthFileName(page) === 'depth.tsv',
    drawn: ['coverage']
  },
  {
    name: 'replacing the annotation data with another legendLabel',
    caption: async () => 'Region X',
    setup: async (page) => {
      await addAnnotationRow(page, 'Region X');
      await colorLegendRow(page, 'Region X');
    },
    change: (page) => inHistoryStep(page, 'Import annotations', async (app, tsv) => {
      await app.importAnnotationTableFile({ target: { files: [new File([tsv], 'annotations.tsv')], value: '' } });
    }, ANNOTATION_TSV('Region Y')),
    restored: (page) => page.evaluate(() => window.__GBDRAW_APP__.annotationSets[0]?.annotations[0]?.legendLabel === 'Region X'),
    drawn: ['Region Y']
  },
  // OV-67: a label edit is a data change of the caption. The edits go through
  // the label fields a reader types into.
  {
    name: 'editing the legend label of an annotation set',
    caption: async () => 'Region X',
    setup: async (page) => {
      await page.evaluate(() => {
        const app = window.__GBDRAW_APP__;
        const set = app.addAnnotationSet('regions');
        app.addCoordinateAnnotation(set, { start: 100, end: 400 });
        set.legendLabel = 'Region X';
      });
      await settleLive(page);
      await generate(page);
      await colorLegendRow(page, 'Region X');
    },
    change: (page) => typeIntoLabel(page, 'Region Annotations', 'Set legend label', 'Region Z'),
    restored: (page) => page.evaluate(() => window.__GBDRAW_APP__.annotationSets[0]?.legendLabel === 'Region X'),
    drawn: ['Region Z']
  },
  {
    name: 'editing the legend title of a Depth series',
    caption: async (page) => depthCaption(page),
    setup: colorDepthRow(),
    change: (page) => typeIntoLabel(page, 'Depth TSV tracks', 'Depth legend title', 'Coverage'),
    restored: (page) => page.evaluate(() => window.__GBDRAW_APP__.adv.depth_tracks[0]?.label === 'depth'),
    drawn: ['Coverage']
  }
];

test.describe('OV-65 Legend styles follow the captions of track data', () => {
  test.beforeEach(() => { test.setTimeout(180_000); });

  for (const { name, mode = 'circular', caption, setup, change, restored, drawn = [] } of RETIRING_CASES) {
    test(`${name} retires the styles of its rows, and Generate succeeds`, async ({ page }) => {
      await openCanvas(page, mode, null);
      await setup(page);
      const row = await caption(page);
      expect(row, 'row caption').toBeTruthy();
      await change(page);
      await settleLive(page);
      expect(await legendStyleOf(page, row)).toEqual(NO_STYLE);
      await generate(page);
      const captions = await drawnLegendCaptions(page);
      expect(captions).not.toContain(row);
      for (const kept of drawn) expect(captions).toContain(kept);
    });

    test(`${name}: Undo restores the data and the style, and live equals Generate`, async ({ page }) => {
      await openCanvas(page, mode, null);
      await setup(page);
      const row = await caption(page);
      await change(page);
      await settleLive(page);
      await undo(page);
      expect(await restored(page), 'data restored').toBe(true);
      expect((await legendStyleOf(page, row)).color).toBe(COLORED);
      await expectLiveEqualsGenerate(page, { label: `${name}, Undo` });
    });
  }

  test('replacing the depth file with the label unchanged keeps the row style', async ({ page }) => {
    await openCanvas(page, 'circular', null);
    await addDepthFile(page, 'depth.tsv');
    await generate(page);
    const caption = await depthCaption(page);
    await colorLegendRow(page, caption);
    await addDepthFile(page, 'depth.tsv');
    expect((await legendStyleOf(page, caption)).color).toBe(COLORED);
    await expectLiveEqualsGenerate(page, { label: 'depth file replaced, same label' });
  });

  test('replacing the annotation data with the legendLabel unchanged keeps the row style', async ({ page }) => {
    await openCanvas(page, 'circular', null);
    await addAnnotationRow(page, 'Region X');
    await colorLegendRow(page, 'Region X');
    await inHistoryStep(page, 'Import annotations', async (app, tsv) => {
      await app.importAnnotationTableFile({ target: { files: [new File([tsv], 'annotations.tsv')], value: '' } });
    }, ANNOTATION_TSV('Region X'));
    await settleLive(page);
    expect((await legendStyleOf(page, 'Region X')).color).toBe(COLORED);
    await expectLiveEqualsGenerate(page, { label: 'annotation data replaced, same legendLabel' });
  });

  // OV-68: the generic Track legend label of a slot renames a GC row. Python
  // lists the default GC caption among the rows it can produce (OV-63), so the
  // Legend color on the old caption stays stored, as for a switched-off track,
  // and Generate succeeds. The case guards that.
  test('renaming the GC content row through its Track legend label keeps Generate working', async ({ page }) => {
    await openCanvas(page, 'circular', null);
    await colorLegendRow(page, 'GC content');
    const panel = page.locator('button[aria-controls="circular-custom-track-slots-panel"]');
    await reveal(panel);
    if (await panel.getAttribute('aria-expanded') !== 'true') await panel.click({ timeout: 15_000 });
    await page.getByText('Use custom stack', { exact: true }).locator('input').check({ timeout: 15_000 });
    await settleLive(page);
    const field = page.getByRole('group', { name: 'Circular track slot gc_content', exact: true })
      .getByRole('textbox', { name: 'Track legend label' });
    await field.fill('Custom GC', { timeout: 15_000 });
    await field.blur();
    await settleLive(page);
    await generate(page);
    const captions = await drawnLegendCaptions(page);
    expect(captions).toContain('Custom GC');
    expect(captions).not.toContain('GC content');
  });
});
