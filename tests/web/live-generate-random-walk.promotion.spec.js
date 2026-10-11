// Promotion check, not a PR or dev check (docs/internal/WEB_PERIODIC_AUDIT.md,
// Promotion PR checklist). A seeded random walk over the edit kinds of the
// PD-OI-066 parity matrix (tests/web/live-generate-parity.playwright.spec.js):
// after every step, the displayed Result must equal the Result the next
// Generate draws from the same draft (expectLiveEqualsGenerate, minus
// tests/web/contracts/live-generate-parity-allowed.json).
//
// Run (one browser, fixed seed, default budget):
//   GBDRAW_RANDOM_WALK_SEED=20261006 npx playwright test --config=tests/web/playwright/promotion.config.js --workers=1
// GBDRAW_RANDOM_WALK_SEED (default 20261006) seeds one PRNG per fixture;
// GBDRAW_RANDOM_WALK_STEPS (default 20) is the number of steps per fixture.
// A failure names the seed, the fixture, the step index, the steps so far, and
// the difference; the same seed and budget replay the same walk. The file
// suffix is `.promotion.spec.js` so that no PR or dev Playwright configuration
// (testMatch `.playwright.spec.js`) collects it.
const { test, expect } = require('@playwright/test');
const { readFileSync } = require('node:fs');
const { evaluateWithRetainedPromise, generateAndWaitForResult } = require('./helpers/app-lifecycle.cjs');
const { BATCH_FIXTURE, openFresh, openWithGenBank } = require('./helpers/audit-browser.cjs');
const { expectLiveEqualsGenerate, settleLive, showResult } = require('./helpers/live-generate-parity.cjs');

test.describe.configure({ retries: 0 });

const SEED = Number.parseInt(process.env.GBDRAW_RANDOM_WALK_SEED || '20261006', 10);
const STEPS = Number.parseInt(process.env.GBDRAW_RANDOM_WALK_STEPS || '20', 10);
if (!Number.isInteger(SEED) || SEED < 0) throw new Error('GBDRAW_RANDOM_WALK_SEED must be a non-negative integer');
if (!Number.isInteger(STEPS) || STEPS < 1) throw new Error('GBDRAW_RANDOM_WALK_STEPS must be a positive integer');
const SINGLE_FIXTURE = 'tests/fixtures/forced_label_underlay.gb';

// mulberry32: a 32-bit seeded generator; the same seed gives the same walk.
const createRandom = (seed) => {
  let state = seed >>> 0;
  const next = () => {
    state = (state + 0x6D2B79F5) >>> 0;
    let value = state;
    value = Math.imul(value ^ (value >>> 15), value | 1);
    value ^= value + Math.imul(value ^ (value >>> 7), value | 61);
    return ((value ^ (value >>> 14)) >>> 0) / 4294967296;
  };
  const int = (count) => Math.floor(next() * count);
  const pick = (items) => items[int(items.length)];
  const weighted = (entries) => {
    const total = entries.reduce((sum, [, weight]) => sum + weight, 0);
    let point = next() * total;
    for (const [value, weight] of entries) {
      point -= weight;
      if (point < 0) return value;
    }
    return entries[entries.length - 1][0];
  };
  return { next, int, pick, weighted };
};

const FIXTURES = [
  { name: 'Circular, one Result', mode: 'circular', results: 'single' },
  { name: 'Linear, one Result', mode: 'linear', results: 'single' },
  { name: 'Circular, two-Result batch', mode: 'circular', results: 'batch' }
];

const COLORS = ['#c83366', '#123456', '#2a9d8f', '#e63946', '#7b2cbf', '#f4a261', '#264653'];
const LABEL_TEXTS = [
  'alpha protein, edited live',
  'x',
  'walk label',
  'a label text far too long to stay inside its own arc, drawn outside with leaders'
];
const CAPTIONS = ['Walk caption', 'transfer RNA', 'Renamed group', 'CDS'];
const escapeRegex = (text) => text.replace(/[.*+?^${}()|[\]\\]/g, '\\$&');

// The setup of the parity spec: one fixture, Generate, the labels drawn
// outside, Auto Reflow as drawn by the walk.
const generate = async (page) => {
  await generateAndWaitForResult(page);
  await settleLive(page);
};

const openFixture = async (page, { mode, results }, reflow) => {
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
  await page.evaluate((enabled) => { window.__GBDRAW_APP__.autoLabelReflowEnabled = enabled; }, reflow);
  await generate(page);
};

// What the next step may use: the features of the displayed Result (sorted by
// a stable key), the editor lists, and the History and Result state. Undo is
// available above `undoFloor`, the History depth after the first Generate:
// Undo of that Generate leaves no Result to compare.
const readState = (page, undoFloor) => page.evaluate((floor) => {
  const app = window.__GBDRAW_APP__;
  const history = window.__GBDRAW_HISTORY__;
  const featureKey = (feature) => [feature.type, feature.start, feature.end, feature.locus_tag || '', feature.product || ''].join('|');
  const features = app.filteredFeatures.map((feature) => ({
    key: featureKey(feature), type: String(feature.type || ''), tag: String(feature.locus_tag || ''), product: String(feature.product || '')
  })).sort((left, right) => (left.key < right.key ? -1 : left.key > right.key ? 1 : 0));
  return {
    features,
    resultCount: app.results.length,
    resultIndex: app.selectedResultIndex,
    reflow: Boolean(app.autoLabelReflowEnabled),
    visibilityRules: app.featureVisibilityManualRules.length,
    colorRules: app.manualSpecificRules.length,
    legend: app.legendEntries.map((entry) => String(entry.caption || '')),
    undoDepth: history.getUndoCount(),
    canUndo: Boolean(history.canUndo()) && history.getUndoCount() > floor
  };
}, undoFloor);

// The next step, drawn from the generator and the state: { kind, ...arguments }.
// `picks` answers the dialogs the edit opens, in the order they open.
const nextStep = (random, state, { mode, results }) => {
  const features = state.features;
  const tagged = features.filter((feature) => feature.tag);
  const feature = () => random.pick(features);
  const rulePattern = () => {
    const target = random.pick(tagged);
    return { target, value: random.next() < 0.5 ? `^${escapeRegex(target.tag)}$` : `${escapeRegex(target.tag.slice(-3))}$` };
  };
  const index = (count) => random.int(count);
  const kinds = [
    ['fill', 5], ['featureVisibility', 2], ['labelText', 3], ['labelVisibility', 3],
    ['visibilityRuleAdd', 2], ['colorRuleAdd', 2], ['reflow', 1]
  ];
  if (state.visibilityRules) kinds.push(['visibilityRuleEdit', 2], ['visibilityRuleRemove', 2]);
  if (state.colorRules) kinds.push(['colorRuleEdit', 2], ['colorRuleRemove', 1]);
  if (state.legend.length) kinds.push(['legendRename', 2], ['legendColor', 2], ['legendSort', 2], ['legendSortDefault', 1]);
  // The check after each step runs Generate, which ends the Redo stack, so Redo
  // is exercised as one step that follows its own Undo.
  if (state.canUndo) kinds.push(['undo', 2], ['undoRedo', 2]);
  if (state.resultCount > 1) kinds.push(['resultSwitch', 2]);
  const kind = random.weighted(kinds);
  const picks = [random.int(1000), random.int(1000), random.int(1000)];
  switch (kind) {
    case 'fill': return { kind, feature: feature().key, color: random.pick(COLORS), picks };
    case 'featureVisibility': return { kind, feature: feature().key, value: random.pick(['on', 'off']), picks };
    case 'labelText': return { kind, feature: feature().key, text: random.pick(LABEL_TEXTS), picks };
    case 'labelVisibility': return { kind, feature: feature().key, value: random.pick(['on', 'off']), picks };
    case 'visibilityRuleAdd': {
      if (!tagged.length) return { kind: 'reflow', enabled: !state.reflow };
      const { target, value } = rulePattern();
      return {
        kind,
        fields: { recordId: '*', featureType: target.type, qualifier: 'locus_tag', value, action: random.pick(['off', 'show']) }
      };
    }
    case 'visibilityRuleEdit':
      return { kind, index: index(state.visibilityRules), field: 'action', value: random.pick(['off', 'show']) };
    case 'visibilityRuleRemove': return { kind, index: index(state.visibilityRules) };
    case 'colorRuleAdd': {
      if (!tagged.length) return { kind: 'reflow', enabled: !state.reflow };
      const { target, value } = rulePattern();
      return {
        kind,
        fields: { feat: target.type, qual: 'locus_tag', val: value, color: random.pick(COLORS), cap: random.pick(CAPTIONS) }
      };
    }
    case 'colorRuleEdit': {
      const field = random.pick(['color', 'cap']);
      return { kind, index: index(state.colorRules), field, value: field === 'color' ? random.pick(COLORS) : random.pick(CAPTIONS), picks };
    }
    case 'colorRuleRemove': return { kind, index: index(state.colorRules) };
    case 'legendRename': return { kind, index: index(state.legend.length), caption: random.pick(CAPTIONS), picks };
    case 'legendColor': return { kind, index: index(state.legend.length), color: random.pick(COLORS) };
    case 'legendSort': return { kind, direction: random.pick(['asc', 'desc']) };
    case 'resultSwitch': {
      const others = [...Array(state.resultCount).keys()].filter((value) => value !== state.resultIndex);
      return { kind, index: random.pick(others) };
    }
    case 'reflow': return { kind, enabled: !state.reflow };
    default: return { kind };
  }
};

// One step through the actions the editor controls call. Scope and Label
// visibility On dialogs that the edit opens are answered from `picks`.
const applyStep = (page, step) => evaluateWithRetainedPromise(page, async (current) => {
  const app = window.__GBDRAW_APP__;
  const history = window.__GBDRAW_HISTORY__;
  const notes = [];
  window.__GBDRAW_WALK_NOTES__ = notes;
  let picked = 0;
  const choose = (choices) => {
    const choice = choices[(current.picks?.[picked++ % current.picks.length] ?? 0) % choices.length];
    return choice;
  };
  const sleep = (ms) => new Promise((resolve) => setTimeout(resolve, ms));
  const featureKey = (feature) => [feature.type, feature.start, feature.end, feature.locus_tag || '', feature.product || ''].join('|');

  // Answers the dialogs the edit opens until the edit has finished and none is open.
  const answerDialogs = async (applied) => {
    // An answer can wait for a later dialog (a hidden label text answer opens
    // Label visibility On), so answers are started, not awaited, in the loop.
    const answers = [];
    let outstanding = 0;
    const answer = (result) => {
      outstanding += 1;
      const settled = Promise.resolve(result).finally(() => { outstanding -= 1; });
      answers.push(settled);
      settled.catch(() => {});
    };
    let finished = false;
    applied.then(() => { finished = true; }, () => { finished = true; });
    const deadline = Date.now() + 150_000;
    // A dialog is answered once, until it has closed.
    let answering = '';
    while (Date.now() < deadline) {
      const style = app.featureStyleScopeDialog;
      const visibility = app.featureVisibilityScopeDialog;
      const shown = ['featureStyleScopeDialog', 'featureVisibilityScopeDialog', 'labelOnDialog', 'labelTextScopeDialog',
        'hiddenLabelTextDialog', 'legendRenameDialog'].find((name) => app[name].show) || '';
      if (!shown) answering = '';
      if (shown && shown === answering) {
        await sleep(20);
        continue;
      }
      answering = shown;
      if (style.show) {
        const offered = ['single'];
        if (style.matchingRule) offered.push('rule');
        if (style.siblingCount > 0) offered.push('caption');
        const text = (value) => String(value || '').trim().toLowerCase();
        if (style.displayLabelSiblingCount > 0 && text(style.displayLabel) && text(style.displayLabel) !== text(style.legendName)) offered.push('displayLabel');
        if (style.annotationLabelSiblingCount > 0 && text(style.annotationLabel) && text(style.annotationLabel) !== text(style.legendName)
          && text(style.annotationLabel) !== text(style.displayLabel)) offered.push('annotationLabel');
        if (style.kind === 'fill' && style.existingCaptionColor) offered.push('useExisting');
        const choice = choose(offered);
        notes.push(`style scope ${choice} of [${offered}]`);
        answer(app.handleFeatureStyleScopeChoice(choice));
      } else if (visibility.show) {
        const offered = (visibility.scopes || []).map((scope) => scope.id);
        const choice = choose(offered);
        notes.push(`visibility scope ${choice} of [${offered}]`);
        answer(app.handleFeatureVisibilityScopeChoice(choice));
      } else if (app.labelOnDialog.show) {
        const reason = app.labelOnDialog.reason;
        const choice = reason === 'hidden' ? choose(['show', 'keep']) : choose(['keep', 'cancel']);
        notes.push(`label On ${reason}: ${choice}`);
        answer(app.handleLabelOnChoice(choice));
      } else if (app.labelTextScopeDialog.show) {
        const choice = choose(['single', 'all']);
        notes.push(`label text scope ${choice}`);
        answer(app.handleLabelTextScopeChoice(choice));
      } else if (app.hiddenLabelTextDialog.show) {
        const choice = choose(['show', 'text_only', 'cancel']);
        notes.push(`hidden label text ${choice}`);
        answer(app.handleHiddenLabelTextChoice(choice));
      } else if (app.legendRenameDialog.show) {
        const choice = choose(app.legendRenameDialog.mode === 'target' ? ['merge', 'suffix'] : ['single', 'group']);
        notes.push(`legend rename ${app.legendRenameDialog.mode} ${choice}`);
        answer(app.handleLegendRenameChoice(choice));
      } else if (finished && outstanding === 0) {
        break;
      } else {
        // Once the edit has settled, its promise wins the race at once: the
        // loop would then never yield a task, and an answer that waits for one
        // (a Worker reply) would never settle while the timers pile up until
        // the renderer runs out of memory (OV-330).
        await (finished ? sleep(20) : Promise.race([applied.catch(() => {}), sleep(20)]));
      }
    }
    await Promise.all([applied, ...answers]);
  };

  const withPopup = async (edit) => {
    const feature = app.filteredFeatures.find((item) => featureKey(item) === current.feature);
    if (!feature) throw new Error(`no feature ${current.feature}`);
    await app.openFeatureEditorFromList(feature, null);
    await window.Vue.nextTick();
    if (!app.clickedFeature) throw new Error(`no popup for ${current.feature}`);
    try {
      await edit();
    } finally {
      app.clickedFeature = null;
    }
  };

  switch (current.kind) {
    case 'fill':
      await withPopup(() => answerDialogs(Promise.resolve(app.updateClickedFeatureColor(current.color))));
      break;
    case 'featureVisibility':
      await withPopup(() => {
        app.clickedFeature.featureVisibility = current.value;
        return answerDialogs(Promise.resolve(app.updateClickedFeatureVisibility(current.value)));
      });
      break;
    case 'labelText':
      await withPopup(() => {
        app.clickedFeature.labelText = current.text;
        return answerDialogs(Promise.resolve(app.updateClickedFeatureLabelText()));
      });
      break;
    case 'labelVisibility':
      await withPopup(() => {
        app.clickedFeature.labelVisibility = current.value;
        return answerDialogs(Promise.resolve(app.updateClickedFeatureLabelText()));
      });
      break;
    case 'visibilityRuleAdd': {
      await app.addFeatureVisibilityRule();
      const index = app.featureVisibilityManualRules.length - 1;
      for (const [field, value] of Object.entries(current.fields)) await app.setFeatureVisibilityRuleField(index, field, value);
      break;
    }
    case 'visibilityRuleEdit':
      await app.setFeatureVisibilityRuleField(current.index, current.field, current.value);
      break;
    case 'visibilityRuleRemove':
      await app.removeFeatureVisibilityRule(current.index);
      break;
    case 'colorRuleAdd':
      Object.assign(app.newSpecRule, current.fields);
      await app.addSpecificRule();
      break;
    case 'colorRuleEdit':
      await answerDialogs(Promise.resolve(app.setSpecificRuleField(current.index, current.field, current.value)));
      break;
    case 'colorRuleRemove':
      await app.removeSpecificRule(current.index);
      break;
    case 'legendRename':
      await answerDialogs(Promise.resolve(app.renameLegendEntry(current.index, current.caption)));
      break;
    case 'legendColor':
      await app.updateLegendEntryColor(current.index, current.color);
      break;
    case 'legendSort':
      await app.sortLegendEntries(current.direction);
      break;
    case 'legendSortDefault':
      await app.sortLegendEntriesByDefault();
      break;
    case 'undo':
      await history.undo();
      break;
    case 'undoRedo':
      await history.undo();
      await history.redo();
      break;
    case 'reflow':
      app.autoLabelReflowEnabled = current.enabled;
      break;
    default:
      throw new Error(`unknown step kind ${current.kind}`);
  }
  return notes;
}, step);

// A step that never settles (for example, a dialog no answer covers) fails with
// the dialogs that are open, instead of waiting for the test deadline.
const STEP_TIMEOUT_MS = 180_000;
const applyStepWithin = async (page, step) => {
  let timer;
  const expired = new Promise((_, reject) => {
    timer = setTimeout(async () => {
      const open = await page.evaluate(() => Object.entries(window.__GBDRAW_APP__)
        .filter(([name, value]) => /Dialog$/.test(name) && value?.show).map(([name]) => name)).catch((error) => [`(page unreadable: ${error?.message || error})`]);
      const answered = await page.evaluate(() => window.__GBDRAW_WALK_NOTES__).catch(() => []);
      reject(new Error(`the step did not settle within ${STEP_TIMEOUT_MS / 1000} s; open dialogs: [${open}]; answered: [${answered}]`));
    }, STEP_TIMEOUT_MS);
  });
  try {
    return await Promise.race([applyStep(page, step), expired]);
  } finally {
    clearTimeout(timer);
  }
};

const describeStep = (index, step, notes = []) => {
  const { picks, ...rest } = step;
  return `#${index} ${JSON.stringify(rest)}${notes.length ? ` (${notes.join('; ')})` : ''}`;
};

test.beforeAll(() => {
  // eslint-disable-next-line no-console
  console.log(`live-generate random walk: seed ${SEED}, ${STEPS} steps per fixture`);
});

for (const [fixtureIndex, fixture] of FIXTURES.entries()) {
  test(`random walk (seed ${SEED}): ${fixture.name}`, async ({ page }) => {
    test.setTimeout(300_000 + STEPS * 240_000);
    // One generator per fixture, so a fixture replays alone with the same seed.
    const random = createRandom(SEED + Math.imul(fixtureIndex + 1, 0x9E3779B1));
    const sequence = [];
    const mismatches = [];
    const report = (index, difference) => (
      `seed ${SEED}, fixture "${fixture.name}", step ${index}\n`
      + `steps so far:\n${sequence.map((line) => `  ${line}`).join('\n')}\n${difference}`
    );

    const reflow = random.next() < 0.5;
    sequence.push(`#init Auto Reflow ${reflow ? 'on' : 'off'}`);
    // eslint-disable-next-line no-console
    console.log(`seed ${SEED}, fixture "${fixture.name}": Auto Reflow ${reflow ? 'on' : 'off'}`);
    await openFixture(page, fixture, reflow);
    expect(await page.evaluate(() => window.__GBDRAW_APP__.results.length)).toBe(fixture.results === 'batch' ? 2 : 1);
    const undoFloor = (await readState(page, 0)).undoDepth;

    for (let index = 1; index <= STEPS; index += 1) {
      const state = await readState(page, undoFloor);
      const step = nextStep(random, state, fixture);
      let notes = [];
      try {
        if (step.kind === 'resultSwitch') await showResult(page, step.index);
        else notes = await applyStepWithin(page, step);
        await settleLive(page);
      } catch (error) {
        sequence.push(describeStep(index, step));
        throw new Error(`${report(index, `the step failed: ${error?.message || error}`)}`
          + `${mismatches.length ? `\nmismatches before it:\n${mismatches.join('\n')}` : ''}`);
      }
      sequence.push(describeStep(index, step, notes));
      // eslint-disable-next-line no-console
      console.log(`  ${sequence[sequence.length - 1]}`);
      try {
        await expectLiveEqualsGenerate(page, { label: `step ${index} (${step.kind})` });
      } catch (error) {
        // expectLiveEqualsGenerate fails with the unexpected differences as `actual`.
        if (!Array.isArray(error?.matcherResult?.actual)) throw error;
        mismatches.push(report(index, `differences:\n${error.matcherResult.actual.map((line) => `  ${line}`).join('\n')}`));
        // Generate has drawn the draft; the walk goes on from that state.
      }
    }
    expect(mismatches.length, `${mismatches.length} of ${STEPS} steps differ from Generate:\n\n${mismatches.join('\n\n')}\n`).toBe(0);
  });
}
