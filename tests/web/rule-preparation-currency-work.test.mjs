// S7 (perf 0.14.x, counts not timings): the currency checks of the rule
// preparation (`snapshot`, `isCurrent`; a rule commit runs about ten) read the
// drawing inputs as JSON computed from the reactive drawing, so the N-sized
// feature color overrides are stringified once per change, not once per check;
// and every kind of write to an input makes an earlier snapshot stale.
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { test } from 'node:test';
import vm from 'node:vm';
import { withModeDrawings } from './helpers/drawing-state.mjs';

const context = vm.createContext({ console });
vm.runInContext(readFileSync(new URL('../../gbdraw/web/vendor/vue/vue.global.js', import.meta.url), 'utf8'), context);
globalThis.window = { Vue: context.Vue };
const { ref, reactive, toRaw } = context.Vue;
const { createRulePreparation } = await import('../../gbdraw/web/js/app/rule-matching.js');

const COUNT = 500;
const drawing = () => ({
  manualSpecificRules: reactive([{ feat: 'CDS', qual: 'product', val: 'p0', color: '#111111', cap: 'CDS' }]),
  legendEntries: ref([{ caption: 'CDS', color: '#111111', featureIds: ['f0'] }]),
  legendColorOverrides: reactive({ CDS: '#111111' }),
  legendStrokeOverrides: reactive({ CDS: { strokeColor: '#000000', strokeWidth: 1 } }),
  featureColorOverrides: reactive(Object.fromEntries(Array.from({ length: COUNT }, (_, index) => (
    [`f${index}`, { color: '#111111', caption: 'CDS' }]
  )))),
  featureOverrides: reactive({
    row: { recordKey: 'r', biologicalFeatureId: 'f0', featureVisibility: 'shown', labelText: 'a' }
  })
});
const setup = () => {
  const drawings = { circular: drawing(), linear: drawing() };
  const state = withModeDrawings({
    mode: ref('circular'), extractedFeatures: ref([]), biologicalFeatures: ref([]), svgResultIdentity: ref('one'),
    cInputType: ref('gb'), lInputType: ref('gb'), files: {}, linearSeqs: [], results: ref([]), selectedResultIndex: ref(0)
  }, drawings);
  return { drawings, state, preparation: createRulePreparation({ state, evaluate: async () => ({}) }) };
};

// JSON.stringify calls whose argument is the given reactive object.
const stringifications = (target, run) => {
  const raw = toRaw(target);
  const stringify = JSON.stringify;
  let count = 0;
  JSON.stringify = (value, ...rest) => {
    if (value && typeof value === 'object' && toRaw(value) === raw) count += 1;
    return stringify(value, ...rest);
  };
  try { run(); } finally { JSON.stringify = stringify; }
  return count;
};

test('the currency checks of one commit stringify the feature color overrides once per change', () => {
  const { drawings, preparation } = setup();
  const overrides = drawings.circular.featureColorOverrides;
  // A commit: its snapshots and checks before the write, then the check that reads the write.
  const before = stringifications(overrides, () => {
    const inputs = preparation.snapshot();
    const candidate = preparation.snapshot();
    for (let check = 0; check < 8; check += 1) assert.equal(preparation.isCurrent(candidate), true);
    assert.equal(preparation.isCurrent(inputs), true);
  });
  assert.equal(before, 1);
  const inputs = preparation.snapshot();
  const after = stringifications(overrides, () => {
    overrides.f1 = { color: '#222222', caption: 'CDS' };
    assert.equal(preparation.isCurrent(inputs), false);
    const next = preparation.snapshot();
    for (let check = 0; check < 8; check += 1) assert.equal(preparation.isCurrent(next), true);
  });
  assert.equal(after, 1);
});

// Each kind of write to a drawing input the matches depend on (the writers
// reach these values through the reactive drawing only).
const WRITES = {
  'rule added': (d) => { d.manualSpecificRules.push({ feat: 'CDS', qual: 'product', val: 'p1', color: '#222222', cap: 'B' }); },
  'rule table replaced in place': (d) => { d.manualSpecificRules.splice(0, 1, { ...d.manualSpecificRules[0], val: 'p9' }); },
  'rule field edited': (d) => { d.manualSpecificRules[0].color = '#333333'; },
  'Legend replaced': (d) => { d.legendEntries.value = [{ caption: 'B', color: '#111111', featureIds: [] }]; },
  'Legend row edited': (d) => { d.legendEntries.value[0].caption = 'B'; },
  'Legend row member added': (d) => { d.legendEntries.value[0].featureIds.push('f9'); },
  'Legend color set': (d) => { d.legendColorOverrides.CDS = '#444444'; },
  'Legend color removed': (d) => { delete d.legendColorOverrides.CDS; },
  'Legend stroke edited': (d) => { d.legendStrokeOverrides.CDS.strokeWidth = 2; },
  'feature color added': (d) => { d.featureColorOverrides.new = { color: '#555555', caption: 'CDS' }; },
  'feature color edited': (d) => { d.featureColorOverrides.f3.color = '#555555'; },
  'feature color removed': (d) => { delete d.featureColorOverrides.f4; },
  'feature hidden': (d) => { d.featureOverrides.row.featureVisibility = 'hidden'; },
  'feature edit added': (d) => {
    d.featureOverrides.other = { recordKey: 'r', biologicalFeatureId: 'f1', featureVisibility: 'hidden' };
  },
  'feature edit removed': (d) => { delete d.featureOverrides.row; }
};

for (const [name, write] of Object.entries(WRITES)) {
  test(`an earlier snapshot is stale after a write: ${name}`, () => {
    const { drawings, preparation } = setup();
    const before = preparation.snapshot();
    assert.equal(preparation.isCurrent(before), true);
    write(drawings.circular);
    assert.equal(preparation.isCurrent(before), false);
    assert.equal(preparation.isCurrent(preparation.snapshot()), true);
    // The other mode's drawing is not an input of the shown one.
    const shown = preparation.snapshot();
    write(drawings.linear);
    assert.equal(preparation.isCurrent(shown), true);
  });
}

test('a mode switch reads the inputs of the drawing it shows', () => {
  const { drawings, state, preparation } = setup();
  drawings.linear.featureColorOverrides.f2 = { color: '#666666', caption: 'CDS' };
  const before = preparation.snapshot();
  state.mode.value = 'linear';
  assert.equal(preparation.isCurrent(before), false);
  assert.equal(preparation.snapshot().featureColors, JSON.stringify(drawings.linear.featureColorOverrides));
});

test('an edit outside the inputs keeps a snapshot current', () => {
  const { drawings, preparation } = setup();
  const before = preparation.snapshot();
  drawings.circular.featureOverrides.row.labelText = 'b';
  assert.equal(preparation.isCurrent(before), true);
});
