import assert from 'node:assert/strict';
import { test } from 'node:test';
import { createFeatureRuleActions } from '../../gbdraw/web/js/app/feature-editor/rule-actions.js';
import { createRulePreparation } from '../../gbdraw/web/js/app/rule-matching.js';
import { createLegendManager } from '../../gbdraw/web/js/app/legend.js';
import { diffLegendIntents } from '../../gbdraw/web/js/services/specific-color-rules.js';
import { evaluatePythonRules } from './helpers/python-rule-evaluator.mjs';

const setup = (evaluate = evaluatePythonRules) => {
  const state = {
    manualSpecificRules: [], extractedFeatures: { value: [
      { type: 'CDS', svg_id: 'a', qualifiers: { gene: ['a'] } },
      { type: 'CDS', svg_id: 'b', qualifiers: { gene: ['b'] } }
    ] },
    featureColorOverrides: {}, featureOverrides: {}, results: { value: [{name:'figure',content:'before'}] },
    svgResultIdentity: { value:'before' }, fileLegendCaptions: { value:new Set() }, addedLegendCaptions: { value:new Set() },
    legendEntries: { value:[] }, files: {t_color:null}, legendColorOverrides: {},
    newSpecRule: {feat:'CDS',qual:'gene',val:'a',color:'#112233',cap:'Shared'}
  };
  const notices = [], transactions = [], transactionScopes = [];
  const preparation = createRulePreparation({state, evaluate, notify:message=>notices.push(message)});
  let prepareLegend = async () => {};
  let previousIntents = [];
  let openTransaction = null;
  const legendApplies = [];
  let rerenders = 0;
  const transact = scope => async (label, commit) => {
    const before=JSON.stringify(state.manualSpecificRules);
    openTransaction = label;
    try { await commit(); } finally { openTransaction = null; }
    if (before!==JSON.stringify(state.manualSpecificRules)) {
      transactions.push(label);
      transactionScopes.push(scope);
    }
  };
  const actions = createFeatureRuleActions({ref:value=>({value}),computed:get=>({get value(){return get();}}),state, rulePreparation:preparation,
    runUndoableCheckpoint: transact('checkpoint'),
    runUndoable: transact('intent'),
    prepareFileLegendEntries: async (intents, {isCurrent,previousFileIntents}) => {
    previousIntents=previousFileIntents;
    await prepareLegend(intents);
    if (!isCurrent()) return false;
    return { diff: diffLegendIntents(state.legendEntries.value, intents), isCurrent: () => true, apply: () => {
      legendApplies.push({ transaction: openTransaction, rules: state.manualSpecificRules.map(rule => rule.cap) });
      state.legendEntries.value=intents;
      state.results.value=[{name:'figure',content:'after'}];
    } };
  }, projectPaletteAndRules:()=>true, ports:{requestAutomaticRerender:()=>{rerenders+=1;return true;}}, nextTick:async()=>{}});
  return {state,actions,preparation,notices,transactions,transactionScopes,legendApplies, setLegendPreparation: fn => {prepareLegend=fn;}, previousIntents:()=>previousIntents, rerenders:()=>rerenders};
};
const rules = [
  {feat:'CDS',qual:'gene',val:'a',color:'#112233',cap:'Shared',fromFile:true},
  {feat:'CDS',qual:'gene',val:'b',color:'#445566',cap:'Shared'}
];

test('rule action admits full canonical rules and legend in one transaction after preparation', async () => {
  const s=setup();
  let release, started;
  const ready=new Promise(resolve=>{started=resolve;});
  s.setLegendPreparation(()=>new Promise(resolve=>{release=resolve;started();}));
  const pending=s.actions.commitSpecificRules(rules);
  await ready;
  assert.deepEqual(s.state.manualSpecificRules,[]);
  assert.equal(s.state.results.value[0].content,'before');
  assert.deepEqual(s.state.legendEntries.value,[]);
  release();
  assert.equal(await pending,true);
  s.setLegendPreparation(async () => {});
  assert.deepEqual(s.state.manualSpecificRules.map(r=>r.cap),['Shared [#112233]','Shared [#445566]']);
  assert.equal(s.state.manualSpecificRules[0].fromFile,true);
  assert.deepEqual([...s.state.fileLegendCaptions.value],['Shared [#112233]']);
  assert.equal(s.transactions.length,1);
  assert.deepEqual(s.legendApplies,[{transaction:'Change specific color rules',rules:['Shared [#112233]','Shared [#445566]']}],
    'the legend rows apply after the rule transition, inside its History step (R13)');
  assert.match(s.notices[0],/Updated 2/);
  await s.actions.setSpecificRuleField(1,'cap','Renamed');
  assert.equal(s.state.manualSpecificRules[1].cap,'Renamed');
  assert.equal(s.state.manualSpecificRules[0].cap,'Shared [#112233]');
  await s.actions.moveSpecificRuleUp(1);
  assert.equal(s.state.manualSpecificRules[0].cap,'Renamed');
  s.state.addedLegendCaptions.value.add('Independent manual legend');
  await s.actions.removeSpecificRule(0);
  assert(s.state.addedLegendCaptions.value.has('Independent manual legend'));
  assert.equal(s.state.manualSpecificRules.length,1);
});

for (const outcome of ['stale','error']) {
  test(`legend ${outcome} candidate leaves canonical rules, provenance, Result and History unchanged`, async () => {
    const s=setup();
    await s.actions.commitSpecificRules(rules);
    s.state.addedLegendCaptions.value.add('Independent manual legend');
    const before=JSON.stringify(s.state.manualSpecificRules), result=s.state.results.value;
    const captions=[...s.state.fileLegendCaptions.value], count=s.transactions.length;
    s.setLegendPreparation(async()=>{
      if(outcome==='error') throw new Error('measurement failed');
      s.state.svgResultIdentity.value='replacement';
    });
    const candidate=rules.map(r=>({...r,cap:'Changed'}));
    if(outcome==='error') await assert.rejects(()=>s.actions.commitSpecificRules(candidate),/measurement failed/);
    else assert.equal(await s.actions.commitSpecificRules(candidate),false);
    assert.equal(JSON.stringify(s.state.manualSpecificRules),before);
    assert.equal(s.state.results.value,result);
    assert.deepEqual([...s.state.fileLegendCaptions.value],captions);
    assert.equal(s.transactions.length,count);
  });
}


test('historical caption ownership retains every source color before normalization',async()=>{
  const s=setup();s.state.manualSpecificRules.push(...rules.map(rule=>({...rule})));
  await s.actions.setSpecificRuleField(0,'val','a');
  assert.deepEqual(s.previousIntents(),[{caption:'Shared',color:'#112233'},{caption:'Shared',color:'#445566'}]);
  assert.deepEqual(s.state.manualSpecificRules.map(rule=>rule.cap),['Shared [#112233]','Shared [#445566]']);
});

test('complete caption recolor uses bounded intent while caption replacement retains checkpoint History', async () => {
  const s = setup();
  s.state.legendEntries.value = [{ caption: 'Shared', color: '#012345' }];
  const recolor = [{ feat: 'CDS', qual: 'gene', val: 'a', color: '#abcdef', cap: 'Shared' }];
  assert.equal(await s.actions.commitSpecificRules(recolor, 'Recolor', {
    previousLegendIntents: [{ caption: 'Shared', color: '#012345' }]
  }), true);
  assert.deepEqual(s.transactionScopes, ['intent']);
  assert.deepEqual(s.state.legendEntries.value, [{ caption: 'Shared', color: '#abcdef' }]);
  assert.deepEqual(s.state.manualSpecificRules.map(rule => rule.cap), ['Shared']);
  assert.equal(s.transactions.length, 1);
  assert.equal(await s.actions.commitSpecificRules(recolor.map(rule => ({ ...rule, cap: 'Renamed' }))), true);
  assert.deepEqual(s.transactionScopes, ['intent', 'checkpoint']);
  assert.deepEqual(s.state.legendEntries.value, [{ caption: 'Renamed', color: '#abcdef' }]);
});

// A rule commit draws a legend row only for a rule a shown, rendered feature
// matches: not for an unmatched rule, a feature only the biological catalog
// holds, or a hidden feature.
test('biological safety rows and hidden rendered features do not create unused legends', async () => {
  const s = setup();
  s.state.biologicalFeatures = { value: [{ type: 'CDS', svg_id: 'unrendered', qualifiers: { gene: ['c'] } }] };
  // Per-feature visibility is the feature's identity row (design Q4).
  Object.assign(s.state.extractedFeatures.value[1], { scope: 'circular', record_key: 'record-1', biological_feature_id: 'bio-b' });
  s.state.featureOverrides[JSON.stringify(['circular', 'record-1', 'bio-b'])] = {
    scope: 'circular', recordKey: 'record-1', biologicalFeatureId: 'bio-b',
    featureVisibility: 'off', labelVisibility: null, labelText: null, labelSourceText: null
  };
  let intents = null;
  s.setLegendPreparation(async (next) => { intents = next; });
  assert.equal(await s.actions.commitSpecificRules([
    { feat: 'CDS', qual: 'gene', val: 'a', color: '#112233', cap: 'Shared' },
    { feat: 'CDS', qual: 'gene', val: 'b', color: '#445566', cap: 'Shared' },
    { feat: 'CDS', qual: 'gene', val: 'c', color: '#778899', cap: 'Shared' },
    { feat: 'CDS', qual: 'gene', val: 'absent', color: '#aabbcc', cap: 'Shared' }
  ]), true);
  assert.deepEqual(s.state.manualSpecificRules.map(rule => rule.cap),
    ['Shared [#112233]', 'Shared [#445566]', 'Shared [#778899]', 'Shared [#aabbcc]']);
  assert.deepEqual(intents, [{ caption: 'Shared [#112233]', color: '#112233' }]);
});

// F-2 (D-14, PD-OI-069): Generate matches a hash rule against the drawn
// feature, whose hash on a cropped or reverse-complemented record differs from
// the source identity in `selector.hash`. The rule "This feature only" writes
// matches that feature live too.
test('a This feature only rule on a cropped record paints its feature live as Generate does', async () => {
  const s = setup();
  const cropped = {
    type: 'CDS', svg_id: 'fbd3d0b74_record_1', stable_feature_id: 'f3ccacda4', record_id: 'TESTA',
    selector: { type: 'CDS', start: 300, end: 600, strand: '+', hash: 'f3ccacda4', qualifiers: { locus_tag: ['TESTA_0001'] } },
    qualifiers: { locus_tag: ['TESTA_0001'] }
  };
  s.state.extractedFeatures.value = [cropped];
  const qualifier = s.actions.getFeatureQualifier(cropped);
  assert.deepEqual(qualifier, { qual: 'hash', val: 'fbd3d0b74' });
  let intents = null;
  s.setLegendPreparation(async (next) => { intents = next; });
  assert.equal(await s.actions.commitSpecificRules([
    { feat: 'CDS', ...qualifier, color: '#c83366', cap: 'duplicate protein' }
  ]), true);
  assert.deepEqual(Object.values(s.state.featureColorOverrides), [{ color: '#c83366', caption: 'duplicate protein' }]);
  assert.deepEqual(intents, [{ caption: 'duplicate protein', color: '#c83366' }]);
});

// N-06 (PD-OI-042): a rule captioned like a generated row of another color is
// drawn as "<caption> [<hex>]". The live commit adds that row, the row is tied to
// the rule, and the generated row is not.
test('a rule captioned like a generated row commits, and is edited through, its suffixed row', async () => {
  const s = setup();
  const trna = { type: 'tRNA', svg_id: 'trna', qualifiers: { product: ['tRNA-Phe'] } };
  s.state.extractedFeatures.value = [{ type: 'CDS', svg_id: 'cds', qualifiers: { gene: ['a'] } }, trna];
  s.state.originalLegendOrder = { value: ['CDS', 'tRNA'] };
  s.state.legendEntries.value = [
    { caption: 'CDS', originalCaption: 'CDS', color: '#54bcf8' },
    { caption: 'tRNA', originalCaption: 'tRNA', color: '#e8b441' }
  ];
  let intents = null;
  s.setLegendPreparation(async (next) => { intents = next; });
  const rule = { feat: 'tRNA', qual: 'product', val: '.*', color: '#ff0000', cap: 'CDS' };
  assert.equal(await s.actions.commitSpecificRules([rule]), true);
  assert.deepEqual(intents, [{ caption: 'CDS [#ff0000]', color: '#ff0000' }]);
  assert.deepEqual(s.state.manualSpecificRules.map((row) => row.cap), ['CDS'], 'the rule keeps its caption');

  // The Result Generate drew: the generated CDS row and the rule's row.
  s.state.originalLegendOrder.value = ['CDS', 'CDS [#ff0000]'];
  s.state.legendEntries.value = [
    { caption: 'CDS', originalCaption: 'CDS', color: '#54bcf8' },
    { caption: 'CDS [#ff0000]', originalCaption: 'CDS [#ff0000]', color: '#ff0000' }
  ];
  assert.equal(s.actions.getEffectiveLegendCaption(trna), 'CDS [#ff0000]');
  assert.deepEqual(s.actions.getLegendRowRules('CDS'), []);
  assert.deepEqual(s.actions.getLegendRowRules('CDS [#ff0000]'), s.state.manualSpecificRules);

  // Recoloring the rule retires the drawn row and adds the row for the new color.
  assert.equal(await s.actions.commitSpecificRules([{ ...rule, color: '#00ff00' }]), true);
  assert.deepEqual(s.previousIntents(), [{ caption: 'CDS [#ff0000]', color: '#ff0000' }]);
  assert.deepEqual(intents, [{ caption: 'CDS [#00ff00]', color: '#00ff00' }]);
});

test('the Legend editor recolors the rule of a suffixed row, and only that row', () => {
  const ref = (value) => ({ value });
  const rule = { feat: 'tRNA', qual: 'product', val: '.*', color: '#ff0000', cap: 'CDS' };
  const state = {
    manualSpecificRules: [rule], svgContainer: ref(null), results: ref([]), selectedResultIndex: ref(0),
    legendEntries: ref([
      { caption: 'CDS', originalCaption: 'CDS', color: '#54bcf8' },
      { caption: 'CDS [#ff0000]', originalCaption: 'CDS [#ff0000]', color: '#ff0000' }
    ]),
    originalLegendOrder: ref(['CDS', 'CDS [#ff0000]']), deletedLegendEntries: ref([]), originalLegendColors: ref({}),
    legendStrokeOverrides: {}, legendColorOverrides: {}, adv: {}
  };
  const committed = [];
  const legend = createLegendManager({ state, commitLegendRowRules: (next, label) => {
    committed.push({ rules: next, label });
    return true;
  } });
  assert.equal(legend.updateLegendEntryColor(1, '#00ff00'), true);
  assert.deepEqual(committed, [{ rules: [{ ...rule, color: '#00ff00' }], label: 'Change legend color' }]);
  assert.equal(legend.updateLegendEntryColor(0, '#123456'), false, 'the generated CDS row is no rule row');
  assert.equal(committed.length, 1);
});

test('a color action prepares the color matches of its rules without the caption evaluation, and a prepared table answers at once', async () => {
  const kinds = [];
  const s = setup(async (payload, options) => { kinds.push(payload.kind); return evaluatePythonRules(payload, options); });
  const commits = [];
  const pending = s.actions.runWithRuleMatches(rules, () => { commits.push('cold'); return 'cold'; });
  assert.ok(pending instanceof Promise, 'unprepared rules are prepared first');
  assert.deepEqual(commits, []);
  assert.equal(await pending, 'cold');
  assert.deepEqual(kinds, ['color'], 'the matches only; commitSpecificRules evaluates the captions');
  assert.equal(s.actions.runWithRuleMatches(rules, () => { commits.push('warm'); return 'warm'; }), 'warm');
  assert.deepEqual(commits, ['cold', 'warm']);
  assert.deepEqual(kinds, ['color']);
});

test('a color action does not commit when its rules were replaced during the preparation, and rejects an invalid pattern', async () => {
  let release;
  const s = setup((payload) => new Promise((resolve) => { release = () => evaluatePythonRules(payload).then(resolve); }));
  let commits = 0;
  const pending = s.actions.runWithRuleMatches(rules, () => { commits++; });
  s.state.manualSpecificRules.push({ ...rules[0], val: 'other' });
  await release();
  await pending;
  assert.equal(commits, 0);
  const invalid = setup();
  await assert.rejects(() => invalid.actions.runWithRuleMatches(
    [{ feat: 'CDS', qual: 'gene', val: '(?<enzyme>a)', color: '#112233', cap: '' }], () => assert.fail('invalid commit')
  ));
});

test('a rule commit retires the Legend color a popup edit left on a row it recolors (OV-61)', async () => {
  const s = setup();
  const rule = { feat: 'CDS', qual: 'gene', val: 'a', color: '#112233', cap: 'Shared' };
  s.state.manualSpecificRules.push({ ...rule });
  s.state.legendColorOverrides.Shared = '#112233';
  s.state.legendColorOverrides.Other = '#abcdef';
  assert.equal(await s.actions.commitSpecificRules([{ ...rule }]), true);
  assert.equal(s.state.legendColorOverrides.Shared, '#112233', 'the same color keeps the row\'s Legend color');
  assert.equal(await s.actions.commitSpecificRules([{ ...rule, color: '#445566' }]), true);
  assert.equal(Object.hasOwn(s.state.legendColorOverrides, 'Shared'), false, 'the recolored row drops the stale Legend color');
  assert.equal(s.state.legendColorOverrides.Other, '#abcdef', 'a row the commit does not draw keeps its Legend color');
});

test('a rule commit that removes a rule row retires the Legend color copied from its rule (OV-152)', async () => {
  const s = setup();
  const shared = { feat: 'CDS', qual: 'gene', val: 'a', color: '#112233', cap: 'Shared' };
  const kept = { feat: 'CDS', qual: 'gene', val: 'b', color: '#445566', cap: 'Kept' };
  const unmatched = { feat: 'CDS', qual: 'gene', val: 'none', color: '#778899', cap: 'Unmatched' };
  s.state.errorLog = { value: null };
  s.state.manualSpecificRules.push({ ...shared }, { ...kept }, { ...unmatched });
  // The popup's copy of the rule color, a Legend color set on a row no rule
  // draws, and copies on rows whose rules stay.
  Object.assign(s.state.legendColorOverrides, {
    Shared: '#112233', Direct: '#abcdef', Kept: '#445566', Unmatched: '#778899'
  });
  assert.equal(await s.actions.removeSpecificRule(0), true);
  assert.equal(Object.hasOwn(s.state.legendColorOverrides, 'Shared'), false, 'the removed rule row drops the copied color');
  assert.deepEqual(s.state.legendColorOverrides, { Direct: '#abcdef', Kept: '#445566', Unmatched: '#778899' },
    'rows no removed rule drew keep their Legend colors, also a rule that draws no feature now');
  assert.equal(s.transactions.length, 1, 'the retirement is part of the rule step');

  // A Legend color that is not the removed rule's color is no copy of it.
  s.state.legendColorOverrides.Kept = '#000000';
  assert.equal(await s.actions.clearAllSpecificRules(), true);
  assert.deepEqual(s.state.legendColorOverrides, { Direct: '#abcdef', Kept: '#000000' },
    'Clear All retires the copies of every rule it removes and keeps the other colors');

  // Rules that recolor a whole type keep its caption; removing the last of
  // them asks Python to draw the type's default row again.
  const typeRule = { feat: 'CDS', qual: 'gene', val: 'a', color: '#112233', cap: 'CDS' };
  assert.equal(await s.actions.commitSpecificRules([{ ...typeRule }, { ...typeRule, val: 'b' }]), true);
  const before = s.rerenders();
  assert.equal(await s.actions.removeSpecificRule(0), true);
  assert.equal(s.rerenders(), before, 'a rule of the type row remains');
  assert.equal(await s.actions.removeSpecificRule(0), true);
  assert.equal(s.rerenders(), before + 1, 'the last rule of the type row asks for the rerender');
});
