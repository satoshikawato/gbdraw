import assert from 'node:assert/strict';
import { test } from 'node:test';
import { createRulePreparation, ruleMatchesFeature, firstMatchingRule } from '../../gbdraw/web/js/app/rule-matching.js';
import { evaluatePythonRules } from './helpers/python-rule-evaluator.mjs';

const setup = (evaluate = evaluatePythonRules) => {
  const features = ['NADH', 'β-lactamase', 'ı', 'other'].map((product, i) => ({
    type: 'CDS', svg_id: `f${i}`, qualifiers: { product: [product] }
  }));
  const state = { extractedFeatures: { value: features }, manualSpecificRules: [], svgResultIdentity: { value: 'one' } };
  return { state, features, preparation: createRulePreparation({ state, evaluate }) };
};
const rule = val => ({ feat: 'CDS', qual: 'product', val });

test('one prepared Python result supplies synchronous membership without JS regex translation', async () => {
  const { features, preparation } = setup();
  const rules = ['(?i)NADH', '(?P<enzyme>NADH)', 'NADH\\Z', '\\bβ', 'i'].map(rule);
  assert.equal(ruleMatchesFeature(features[0], rules[0]), null);
  assert.equal(await preparation.prepare(rules), true);
  assert.deepEqual(rules.map(r => features.map(f => ruleMatchesFeature(f, r))), [
    [true, false, false, false], [true, false, false, false], [true, false, false, false],
    [false, true, false, false], [false, false, true, false]
  ]);
  let committed = false;
  preparation.run(rules, () => { committed = true; });
  assert.equal(committed, true, 'prepared edits commit before returning');
});

test('invalid syntax with an empty or unrelated catalog cannot reach the commit', async () => {
  const { state, preparation } = setup();
  for (const features of [state.extractedFeatures.value, []]) {
    state.extractedFeatures.value = features;
    await assert.rejects(() => preparation.run([rule('(?<enzyme>NADH)')], () => assert.fail('invalid commit')));
  }
});

for (const replace of [s => { s.svgResultIdentity.value = 'two'; }, s => { s.manualSpecificRules.push(rule('other')); }, s => { s.extractedFeatures.value = []; }]) {
  test('late evaluation cannot commit to a replaced artifact, rules or catalog', async () => {
    let release;
    const { state, preparation } = setup(payload => new Promise(resolve => { release = () => evaluatePythonRules(payload).then(resolve); }));
    let commits = 0;
    const pending = preparation.run([rule('NADH')], () => { commits++; });
    replace(state);
    await release();
    await pending;
    assert.equal(commits, 0);
  });
}


test('prepared precedence follows Python even when qualifier and wildcard rules interleave', async () => {
  const { features, preparation } = setup();
  const rules = [
    { feat: '*', qual: 'hash', val: 'f0' },
    { feat: 'CDS', qual: 'product', val: 'NADH' },
    { feat: 'CDS', qual: 'hash', val: 'f0' }
  ];
  await preparation.prepare(rules);
  assert.equal(firstMatchingRule(features[0], rules), rules[2]);
  assert.equal(firstMatchingRule(features[0], rules.slice(0, 2)), rules[1]);
});

test('25,000-feature matching is prepared once and reused synchronously', async () => {
  let calls = 0;
  const { state, preparation } = setup(async payload => { calls++; return evaluatePythonRules(payload); });
  state.extractedFeatures.value = Array.from({ length: 25000 }, (_, i) => ({
    type: 'CDS', svg_id: `f${i}`, qualifiers: { product: [i % 2 ? 'other' : 'β-lactamase'] }
  }));
  const rules = [rule('\\bβ'), rule('(?i)lactamase')];
  const start = performance.now();
  await preparation.prepare(rules);
  const preparedMs = performance.now() - start;
  const matchStart = performance.now();
  assert.equal(state.extractedFeatures.value.filter(f => ruleMatchesFeature(f, rules[0])).length, 12500);
  assert.equal(preparation.prepare(rules), true);
  assert.equal(calls, 1);
  console.log(JSON.stringify({ featureCount: 25000, preparedMs, synchronousReuseMs: performance.now() - matchStart }));
});
