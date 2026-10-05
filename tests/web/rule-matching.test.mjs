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

test('full candidate normalizes captions on matching cache hits and retains provenance', async () => {
  const messages = [];
  const { state, features } = setup();
  const preparation = createRulePreparation({ state, evaluate: evaluatePythonRules, notify: message => messages.push(message) });
  const source = [
    { ...rule('NADH'), color: '#112233', cap: 'Shared', fromFile: true },
    { ...rule('other'), color: '#445566', cap: 'Shared' },
    { ...rule('absent'), color: '#778899', cap: 'Shared' }
  ];
  await preparation.prepare(source);
  const candidate = await preparation.prepareCandidate(source);
  assert.deepEqual(candidate.rules.map(r => r.cap), ['Shared [#112233]', 'Shared [#445566]', 'Shared [#778899]']);
  assert.equal(candidate.rules[0].fromFile, true);
  assert.equal(firstMatchingRule(features[0], candidate.rules).cap, 'Shared [#112233]');
  assert.equal(preparation.prepare(source), true, 'matching remains synchronous');
  preparation.notifyChanges(candidate);
  assert.match(messages[0], /Updated 3/);
  const edited = source.map((r, i) => i === 1 ? { ...r, color: '#112233' } : r);
  const next = await preparation.prepareCandidate(edited);
  assert.equal(next.rules[1].cap, 'Shared [#112233]');
  const repeated = await preparation.prepareCandidate(next.rules);
  assert.deepEqual(repeated.rules, next.rules);
  assert.deepEqual(repeated.changes, []);
});

for (const replace of [s => { s.files.t_color = {}; }, s => { s.results.value = [{ name: 'replacement', content: '<svg/>' }]; }, s => { s.legendEntries.value = [{caption:'Manual'}]; }, s => { s.linearSeqs.push({gb:{}}); }, s => { s.linearSeqs[0].gff = {}; }, s => { s.cInputType.value='gff'; }]) {
  test('caption preparation rejects replaced source, Result and legend snapshots', async () => {
    const { state } = setup();
    state.files = { t_color: null };
    state.linearSeqs = [{gb:null,gff:null,fasta:null}];
    state.cInputType = {value:'gbk'};
    state.results = { value: [] };
    state.legendEntries = { value: [] };
    let release;
    const preparation = createRulePreparation({ state, evaluate: payload => new Promise(resolve => { release = () => evaluatePythonRules(payload).then(resolve); }) });
    const pending = preparation.prepareCandidate([{...rule('NADH'), color:'#112233', cap:'Shared'}]);
    replace(state);
    await release();
    assert.equal(await pending, null);
  });
}


test('canceled caption helper cannot admit rules or emit a notification', async () => {
  const {state}=setup();
  const notices=[];
  const preparation=createRulePreparation({state,evaluate:async()=>{throw new Error('canceled');},notify:value=>notices.push(value)});
  await assert.rejects(()=>preparation.prepareCandidate([{...rule('NADH'),color:'#112233',cap:'Shared'}]),/canceled/);
  assert.deepEqual(state.manualSpecificRules,[]);
  assert.deepEqual(notices,[]);
});

test('historical rule-derived overrides rebind by source caption and color while direct overrides survive',async()=>{
  const {state,preparation}=setup();
  state.featureColorOverrides={a:{caption:'Shared',color:'#112233'},b:{caption:'Shared',color:'#445566'},direct:{caption:'Manual',color:'#abcdef'}};
  const candidate=await preparation.prepareCandidate([
    {...rule('NADH'),cap:'Shared',color:'#112233'}, {...rule('other'),cap:'Shared',color:'#445566'}
  ]);
  assert.deepEqual(candidate.featureColorOverrides,{a:{caption:'Shared [#112233]',color:'#112233'},b:{caption:'Shared [#445566]',color:'#445566'},direct:state.featureColorOverrides.direct});
  assert.equal(state.featureColorOverrides.a.caption,'Shared');
});


test('Generate preparation forwards its progress observer through caption and membership evaluation', async () => {
  const observations = [];
  const options = { onProgress: event => observations.push(event) };
  const kinds = [];
  const { preparation } = setup(async (payload, received) => {
    assert.equal(received, options);
    kinds.push(payload.kind);
    received.onProgress({ stage: 'preparing-runtime', requestId: 17 });
    return evaluatePythonRules(payload);
  });
  const candidate = await preparation.prepareCandidate([
    { ...rule('NADH'), color: '#112233', cap: 'NADH' }
  ], options);
  assert.deepEqual(kinds, ['color-captions', 'color']);
  assert.equal(candidate.rules[0].cap, 'NADH');
  assert.equal(observations.length, 2);
});


test('derived comparison cache invalidation keeps prepared colors current while source replacement rejects them', async () => {
  let release;
  const gate = new Promise(resolve => { release = resolve; });
  const { state, preparation } = setup(async payload => {
    if (payload.kind === 'color-captions') await gate;
    return evaluatePythonRules(payload);
  });
  state.files = { c_gb: {}, t_color: null, linearCanonicalComparisons: [] };
  const pending = preparation.prepareCandidate([{ ...rule('NADH'), color: '#112233', cap: 'NADH' }]);
  state.files.linearCanonicalComparisons = [];
  release();
  const candidate = await pending;
  assert.ok(candidate, 'invalidating request-owned comparison reuse cannot invalidate color rules');
  assert.equal(preparation.isCurrent(candidate.snapshot), true);
  state.files.c_gb = {};
  assert.equal(preparation.isCurrent(candidate.snapshot), false, 'replaced biological sources remain guarded');
});


for (const replace of [
  state => { state.files.c_gb = new File(['changed bases'], 'same.gb', { lastModified: 123 }); },
  state => { state.linearSeqs[0].gb = new File(['changed bases'], 'same.gb', { lastModified: 123 }); },
  state => { state.files.c_conservation_fastas[0] = new File(['changed bases'], 'same.gb', { lastModified: 123 }); }
]) {
  test('pending preparation rejects same-named biological replacement, including an in-place source list edit', async () => {
    let release;
    const gate = new Promise(resolve => { release = resolve; });
    const original = new File(['original bases'], 'same.gb', { lastModified: 123 });
    const { state, preparation } = setup(async payload => {
      await gate;
      return { rules: payload.rules };
    });
    state.files = { c_gb: original, c_conservation_fastas: [original], t_color: null };
    state.linearSeqs = [{ gb: original, gff: null, fasta: null }];
    const catalog = state.extractedFeatures.value;
    const pending = preparation.prepareCandidate([]);
    replace(state);
    release();
    assert.equal(await pending, null);
    assert.equal(state.extractedFeatures.value, catalog);
    assert.deepEqual(state.manualSpecificRules, []);
  });
}


for (const owner of ['extractedFeatures', 'biologicalFeatures']) {
  test('pending preparation rejects a replaced catalog owner even when its feature members are unchanged', async () => {
    let release;
    const gate = new Promise(resolve => { release = resolve; });
    const { state, preparation } = setup(async payload => {
      await gate;
      return { rules: payload.rules };
    });
    state.biologicalFeatures = { value: state.extractedFeatures.value };
    const pending = preparation.prepareCandidate([]);
    state[owner].value = [...state[owner].value];
    release();
    assert.equal(await pending, null);
  });
}

test('prepared-rule readiness distinguishes warm config, new predicates, and new feature objects', async () => {
  let calls = 0;
  const { state, preparation } = setup(async payload => { calls++; return evaluatePythonRules(payload); });
  const first = { ...rule('NADH'), color: '#00662c' };
  state.manualSpecificRules.push(first);
  assert.equal(preparation.isPrepared(), false);
  assert.equal(await preparation.prepare(), true);
  assert.equal(calls, 1);
  assert.equal(preparation.isPrepared(), true);
  first.color = '#ff0000';
  assert.equal(preparation.isPrepared(), true, 'color-only changes reuse the validated match');
  first.val = 'other';
  assert.equal(preparation.isPrepared(), false);
  assert.equal(await preparation.prepare(), true);
  assert.equal(calls, 2);
  state.extractedFeatures.value = state.extractedFeatures.value.map(feature => ({ ...feature }));
  assert.equal(preparation.isPrepared(), false);
  assert.equal(await preparation.prepare(), true);
  assert.equal(calls, 3);
});
