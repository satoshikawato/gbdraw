import assert from 'node:assert/strict';
import { test } from 'node:test';
import { createFeatureRuleActions } from '../../gbdraw/web/js/app/feature-editor/rule-actions.js';
import { createRulePreparation } from '../../gbdraw/web/js/app/rule-matching.js';
import { createLegendManager } from '../../gbdraw/web/js/app/legend.js';
import { readFileSync } from 'node:fs';
import {
  diffLegendIntents, generatedLegendRow, legendRowRules, ruleCommitLegendRows, ruleLegendCaptions
} from '../../gbdraw/web/js/services/specific-color-rules.js';
import { recordRuleMatches, ruleKey } from '../../gbdraw/web/js/services/rule-matchers.js';
import { setFeatureVisibilityOverride } from '../../gbdraw/web/js/services/feature-visibility.js';
import { evaluatePythonRules } from './helpers/python-rule-evaluator.mjs';
import { withDrawings } from './helpers/drawing-state.mjs';

const setup = (evaluate = evaluatePythonRules) => {
  const state = {
    manualSpecificRules: [], extractedFeatures: { value: [
      { type: 'CDS', svg_id: 'a', qualifiers: { gene: ['a'] } },
      { type: 'CDS', svg_id: 'b', qualifiers: { gene: ['b'] } }
    ] },
    featureColorOverrides: {}, featureOverrides: {}, results: { value: [{name:'figure',content:'before'}] },
    svgResultIdentity: { value:'before' }, fileLegendCaptions: { value:new Set() }, addedLegendCaptions: { value:new Set() },
    legendEntries: { value:[] }, files: {t_color:null}, legendColorOverrides: {},
    newSpecRule: {feat:'CDS',qual:'gene',val:'a',color:'#112233',cap:'Shared'}, labelReflowProcessing: { value:false }
  };
  // The automatic rerender's busy flag, with the Vue `watch` the owner waits on.
  const watchers = new Set();
  const watch = (source, callback) => { const entry = {source, callback}; watchers.add(entry); return () => watchers.delete(entry); };
  const setRerendering = (busy) => {
    state.labelReflowProcessing.value = busy;
    [...watchers].filter(entry => entry.source === state.labelReflowProcessing).forEach(entry => entry.callback(busy));
  };
  const notices = [], transactions = [], transactionScopes = [];
  const preparation = createRulePreparation({state: withDrawings(state), evaluate, notify:message=>notices.push(message)});
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
  const actions = createFeatureRuleActions({ref:value=>({value}),computed:get=>({get value(){return get();}}),state: withDrawings(state), rulePreparation:preparation,
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
  }, projectPaletteAndRules:()=>true, ports:{requestAutomaticRerender:()=>{rerenders+=1;return true;}}, nextTick:async()=>{}, watch});
  return {state,actions,preparation,notices,transactions,transactionScopes,legendApplies, setLegendPreparation: fn => {prepareLegend=fn;}, previousIntents:()=>previousIntents, rerenders:()=>rerenders, setRerendering};
};
// Review 3 (OV-262 in the Features drawer): with the type's default color
// Auto, the drawer's color input shows the applied palette's color, as
// Generate draws it.
test('getFeatureColor reads the applied palette color for a type set to Auto', () => {
  const s = setup();
  Object.assign(s.state, {
    paletteDefinitions: { value: { default: { CDS: '#5b8fd1' } } },
    appliedPaletteName: { value: 'default' }, appliedPaletteColors: { value: { CDS: null } }
  });
  assert.equal(s.actions.getFeatureColor(s.state.extractedFeatures.value[0]), '#5b8fd1');
});

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
    let replacements=0;
    s.setLegendPreparation(async()=>{
      if(outcome==='error') throw new Error('measurement failed');
      // Every preparation goes stale, so no attempt may commit (OV-166 retries).
      s.state.svgResultIdentity.value=`replacement-${++replacements}`;
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


// OV-166: a rule commit made during an automatic rerender waits for the
// rerender, but its preparation can still span the binding of the Result the
// rerender wrote: the binder reads that Result's Legend rows, which changes
// the inputs the candidate was prepared from while the Results stay the same.
// The commit prepares once more against the bound Result and applies; it is
// not dropped.
test('a rule commit whose preparation spans the rerendered Result binding prepares again and applies (OV-166)', async () => {
  const s=setup();
  const results=s.state.results.value;
  let preparations=0;
  s.setLegendPreparation(async()=>{
    preparations+=1;
    // The binder's extraction of the rerendered Result, once.
    if(preparations===1) s.state.legendEntries.value=[{caption:'CDS',originalCaption:'CDS',color:'#808080'}];
  });
  assert.equal(await s.actions.commitSpecificRules(rules),true);
  assert.equal(preparations,2,'one more preparation, against the bound Result');
  assert.deepEqual(s.state.manualSpecificRules.map(r=>r.cap),['Shared [#112233]','Shared [#445566]']);
  assert.equal(s.transactions.length,1);
  assert.notEqual(results,s.state.results.value,'the Legend rows applied');
});

// OV-280: a run whose rules Python must still match (a color choice's
// candidate rules) waits for the automatic rerender of an earlier edit. The
// rerender replaces the Result, and its binding the Legend rows (`other
// proteins`); a preparation that spanned them went stale and the run did
// nothing. A run whose rules are prepared still runs in the caller's tick
// (PD-OI-088: the dialog opens first).
test('a rule run waits for the automatic rerender in flight before it prepares (OV-280)', async () => {
  const s=setup();
  const [saved, candidate]=[{feat:'CDS',qual:'gene',val:'a',color:'#112233',cap:'CDS'},{feat:'CDS',qual:'gene',val:'b',color:'#445566',cap:'CDS'}];
  assert.equal(await s.preparation.prepare([saved]),true);
  s.setRerendering(true);
  assert.equal(s.actions.runWithRuleMatches([saved],()=>'opened'),'opened','prepared rules run at once');
  let runs=0;
  const run=s.actions.runWithRuleMatches([saved,candidate],()=>{runs+=1;return 'applied';});
  // The rerender lands: its Result, the bound Legend rows, then idle.
  s.state.svgResultIdentity.value='rerendered';
  s.state.legendEntries.value=[{caption:'CDS',originalCaption:'CDS',color:'#112233'},{caption:'other proteins',originalCaption:'other proteins',color:'#808080'}];
  assert.equal(runs,0,'not before the rerender is idle');
  s.setRerendering(false);
  assert.equal(await run,'applied');
  assert.equal(runs,1);
  assert(s.preparation.isPrepared([saved,candidate]));
});

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
  s.state.featureOverrides[JSON.stringify(['record-1', 'bio-b'])] = {
    recordKey: 'record-1', biologicalFeatureId: 'bio-b',
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
// The mounted Result as `pythonLegendRows` (services/legend-svg.js) reads it:
// Python's Legend rows, each its key and the fill Python drew its swatch in.
const resultSvg = (rows) => {
  const element = (attributes, children = []) => ({
    getAttribute: (name) => attributes[name] ?? null,
    querySelector: () => null,
    querySelectorAll: () => children
  });
  const legend = element({}, rows.map(([key, fill]) => element({ 'data-legend-key': key }, [element({ fill })])));
  return { getElementById: (id) => (id === 'legend' ? legend : null) };
};
const mounted = (rows) => {
  const svg = resultSvg(rows);
  return { value: { querySelector: () => svg } };
};

// OV-294 A (R15-6): the rules a Legend row draws are read from Python's rows
// of the displayed Result, never from the row's swatch. OV-307: a row Python
// drew for a rule keeps the fill Python drew it in after a live edit of the
// rule's color, and still is that rule's row, unless its caption names a
// feature type's row.
test('the rules a Legend row draws follow Python\'s rows, also after a live recolor of a rule', () => {
  const rows = (entries) => new Map(entries.map(([key, color]) => [key, { key, color }]));
  const fl2 = { type: 'CDS', svg_id: 'fl2', fill_color: '#2266aa', qualifiers: {} };
  const named = { feat: 'CDS', qual: 'locus_tag', val: '^FL2$', color: '#2266aa', cap: 'repeat_region' };
  const generated = {
    rules: [named], features: [fl2], originalLegendOrder: ['CDS', 'repeat_region', 'repeat_region [#2266aa]'],
    pythonRows: rows([['CDS', '#54bcf8'], ['repeat_region', '#d3d3d3'], ['repeat_region [#2266aa]', '#2266aa']])
  };
  recordRuleMatches([fl2], [ruleKey(named)], () => ({ matched: [0], priorities: [0], declined: [] }));
  assert.deepEqual(legendRowRules('repeat_region', generated), [], 'the row the caption names is no rule row');
  assert.deepEqual(legendRowRules('repeat_region [#2266aa]', generated), [named]);

  const recolored = { ...named, cap: 'Group', color: '#11aa55' };
  recordRuleMatches([fl2], [ruleKey(recolored)], () => ({ matched: [0], priorities: [0], declined: [] }));
  const group = {
    rules: [recolored], features: [fl2], originalLegendOrder: ['CDS', 'Group'],
    pythonRows: rows([['CDS', '#54bcf8'], ['Group', '#2266aa']])
  };
  assert.deepEqual(legendRowRules('Group', group), [recolored], 'Python drew the row for the rule');
  // A caption naming a feature type's row is compared with that row, whatever its fill.
  const other = { ...recolored, cap: 'other proteins' };
  recordRuleMatches([fl2], [ruleKey(other)], () => ({ matched: [0], priorities: [0], declined: [] }));
  const otherRows = { ...group, rules: [other], pythonRows: rows([['CDS', '#54bcf8'], ['other proteins', '#2266aa']]),
    originalLegendOrder: ['CDS', 'other proteins'] };
  assert.deepEqual(legendRowRules('other proteins', otherRows), []);
  assert.deepEqual(legendRowRules('other proteins [#11aa55]', otherRows), [other]);
});

test('a rule captioned like a generated row commits, and is edited through, its suffixed row', async () => {
  const s = setup();
  const trna = { type: 'tRNA', svg_id: 'trna', qualifiers: { product: ['tRNA-Phe'] } };
  s.state.extractedFeatures.value = [{ type: 'CDS', svg_id: 'cds', qualifiers: { gene: ['a'] } }, trna];
  s.state.originalLegendOrder = { value: ['CDS', 'tRNA'] };
  s.state.legendEntries.value = [
    { caption: 'CDS', originalCaption: 'CDS', color: '#54bcf8' },
    { caption: 'tRNA', originalCaption: 'tRNA', color: '#e8b441' }
  ];
  s.state.svgContainer = mounted([['CDS', '#54bcf8'], ['tRNA', '#e8b441']]);
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
  s.state.svgContainer = mounted([['CDS', '#54bcf8'], ['CDS [#ff0000]', '#ff0000']]);
  assert.equal(s.actions.getEffectiveLegendCaption(trna), 'CDS [#ff0000]');
  assert.deepEqual(s.actions.getLegendRowRules('CDS'), []);
  assert.deepEqual(s.actions.getLegendRowRules('CDS [#ff0000]'), s.state.manualSpecificRules);

  // Recoloring the rule retires the drawn row and adds the row for the new color.
  assert.equal(await s.actions.commitSpecificRules([{ ...rule, color: '#00ff00' }]), true);
  assert.deepEqual(s.previousIntents(), [{ caption: 'CDS [#ff0000]', color: '#ff0000' }]);
  assert.deepEqual(intents, [{ caption: 'CDS [#00ff00]', color: '#00ff00' }]);
});

// OV-308 (OV-306 on the commit path): Python draws no row for a feature type
// once a rule captioned with that type colors a feature, so the rule draws the
// row under the type's name in its own color. The commit adds that row, not
// "<type> [<hex>]", and replaces the type's row Python drew, which it owns.
test('a rule captioned with its own feature type commits the type\'s row in its color', async () => {
  const s = setup();
  s.state.originalLegendOrder = { value: ['CDS'] };
  s.state.legendEntries.value = [{ caption: 'CDS', originalCaption: 'CDS', color: '#54bcf8' }];
  s.state.svgContainer = mounted([['CDS', '#54bcf8']]);
  let intents = null;
  s.setLegendPreparation(async (next) => { intents = next; });
  const rule = { feat: 'CDS', qual: 'gene', val: '.', color: '#c83366', cap: 'CDS' };
  assert.equal(await s.actions.commitSpecificRules([rule]), true);
  assert.deepEqual(intents, [{ caption: 'CDS', color: '#c83366' }]);
  assert.deepEqual(s.previousIntents(), [{ caption: 'CDS', color: '#54bcf8' }], 'the commit owns the type\'s row');
});

// OV-310: the commit owns the listed row of a Python row its rules take by
// the row's key, and diffs it by the caption and the swatch the row shows: a
// Legend color on the type's row is no other row of that caption.
test('a rule captioned with its own feature type commits while the type\'s row shows a Legend color', async () => {
  const s = setup();
  const ref = (value) => ({ value });
  Object.assign(s.state, {
    originalLegendOrder: ref(['CDS']), deletedLegendEntries: ref([]), dormantLegendEntries: ref([]),
    originalLegendColors: ref({}), legendStrokeOverrides: {}, legendColorOverrides: { CDS: '#ff8800' }, adv: {},
    selectedResultIndex: ref(0), svgContainer: mounted([['CDS', '#54bcf8']])
  });
  s.state.legendEntries.value = [{ caption: 'CDS', originalCaption: 'CDS', color: '#ff8800' }];
  let intents = null;
  s.setLegendPreparation(async (next) => { intents = next; });
  const rule = { feat: 'CDS', qual: 'gene', val: '.', color: '#c83366', cap: 'CDS' };
  assert.equal(await s.actions.commitSpecificRules([rule]), true);
  assert.deepEqual(intents, [{ caption: 'CDS', color: '#c83366' }]);
  assert.deepEqual(s.previousIntents(), [{ caption: 'CDS', color: '#ff8800' }], 'the listed row, by its key');
  // The Legend preparation owns that row and updates it.
  s.state.legendEntries.value = [{ caption: 'CDS', originalCaption: 'CDS', color: '#ff8800' }];
  const legend = createLegendManager({ state: withDrawings(s.state), commitLegendRowRules: () => true });
  const prepared = await legend.prepareFileLegendEntries(intents, { previousFileIntents: s.previousIntents() });
  assert.deepEqual(prepared.diff.update.map(({ caption, color }) => ({ caption, color })), [{ caption: 'CDS', color: '#c83366' }]);
});

// L2 (OV-306 to OV-308 in batch): the rows the rules take are read from the
// features the displayed Result draws: a rule that colors a feature only in
// another batch Result, or a feature hidden since, takes no row of this one.
test('a rule commit owns no type row of the displayed Result for a feature it does not draw', async () => {
  const s = setup();
  const a1 = { type: 'CDS', svg_id: 'a1', fill_color: '#54bcf8', qualifiers: { locus_tag: ['A1'] } };
  const fl2 = { type: 'CDS', svg_id: 'fl2', fill_color: '#54bcf8', qualifiers: { locus_tag: ['FL2'] } };
  const hidden = { type: 'CDS', svg_id: 'a2', fill_color: '#54bcf8', qualifiers: { locus_tag: ['A2'] }, record_key: 'A', biological_feature_id: 'a2' };
  s.state.extractedFeatures.value = [a1, hidden, fl2];
  s.state.results.value = [{ name: 'A', content: 'a' }, { name: 'B', content: 'b' }];
  setFeatureVisibilityOverride(s.state.featureOverrides, hidden, 'off');
  Object.assign(s.state, {
    displayedResultMetadata: () => ({ renderedFeatureIdentities: { renderedIds: new Set(['a1', 'a2']) } }),
    originalLegendOrder: { value: ['CDS'] }, svgContainer: mounted([['CDS', '#54bcf8']])
  });
  s.state.legendEntries.value = [{ caption: 'CDS', originalCaption: 'CDS', color: '#54bcf8' }];
  s.setLegendPreparation(async () => {});
  const rule = { feat: 'CDS', qual: 'locus_tag', val: '^(FL2|A2)$', color: '#2266aa', cap: 'Group' };
  assert.equal(await s.actions.commitSpecificRules([rule]), true);
  assert.deepEqual(s.previousIntents(), [], 'Python keeps the CDS row of the displayed Result');
});

// L1: a rule captioned with a track's row (Python's own rows,
// `_generated_legend_fills`) is no rule of that row, also in the track's
// color: Python draws it in "<caption> [<hex>]" and keeps the track's row.
test('a rule captioned with a track row draws its own row, and the commit leaves the track row alone', () => {
  const rows = (entries) => new Map(entries.map(([key, color]) => [key, { key, color }]));
  const fl2 = { type: 'CDS', svg_id: 'fl2', fill_color: '#a0a0a0', qualifiers: {} };
  const live = { feat: 'CDS', qual: 'locus_tag', val: '^FL2$', color: '#ff0000', cap: 'GC content' };
  const next = { ...live, color: '#00ff00' };
  recordRuleMatches([fl2], [ruleKey(live), ruleKey(next)], () => ({ matched: [0, 1], priorities: [0, 0], declined: [] }));
  const context = {
    rules: [live], features: [fl2], originalLegendOrder: ['other proteins', 'GC content'],
    pythonRows: rows([['other proteins', '#54bcf8'], ['GC content', '#a0a0a0']])
  };
  assert.equal(ruleLegendCaptions(context)(live), 'GC content [#ff0000]', 'Python\'s key for the rule');
  assert.deepEqual(legendRowRules('GC content', context), [], 'the GC content row is the track\'s');
  assert.deepEqual([...ruleCommitLegendRows({ ...context, rules: [live, next] }).takenKeys], []);
});

// The Python rows the Web app reads from their keys, with the palette key of
// their fill, are those Python generates (tests/test_legend_row_facts.py).
test('the generated Legend rows read from their keys are Python\'s', () => {
  const { cases } = JSON.parse(readFileSync(new URL('../fixtures/legend_generated_rows.json', import.meta.url), 'utf8'));
  for (const { name, features_present: types, rows } of cases) {
    for (const [key, paletteKey] of rows) {
      const row = generatedLegendRow(key);
      assert.equal(row.paletteKey, paletteKey, `${name}: ${key}`);
      assert.ok(row.track || types.includes(row.paletteKey), `${name}: ${key} is a track row or a present type's row`);
    }
  }
  assert.deepEqual(generatedLegendRow('Group'), { track: false, paletteKey: 'Group' }, 'any other key names a type\'s row');
});

test('the Legend editor recolors the rule of a suffixed row, and only that row', () => {
  const ref = (value) => ({ value });
  const rule = { feat: 'tRNA', qual: 'product', val: '.*', color: '#ff0000', cap: 'CDS' };
  const state = {
    manualSpecificRules: [rule], svgContainer: mounted([['CDS', '#54bcf8'], ['CDS [#ff0000]', '#ff0000']]),
    results: ref([]), selectedResultIndex: ref(0),
    legendEntries: ref([
      { caption: 'CDS', originalCaption: 'CDS', color: '#54bcf8' },
      { caption: 'CDS [#ff0000]', originalCaption: 'CDS [#ff0000]', color: '#ff0000' }
    ]),
    originalLegendOrder: ref(['CDS', 'CDS [#ff0000]']), deletedLegendEntries: ref([]), originalLegendColors: ref({}),
    legendStrokeOverrides: {}, legendColorOverrides: {}, adv: {}
  };
  const committed = [];
  const legend = createLegendManager({ state: withDrawings(state), commitLegendRowRules: (next, label) => {
    committed.push({ rules: next, label });
    return true;
  }, readShownLegendColor: (entry) => entry?.color });
  assert.equal(legend.updateLegendEntryColor(1, '#00ff00'), true);
  assert.deepEqual(committed, [{ rules: [{ ...rule, color: '#00ff00' }], label: 'Change legend color' }]);
  assert.deepEqual([legend.legendRowHasRules(0), legend.legendRowHasRules(1)], [false, true]);
  assert.equal(legend.updateLegendEntryColor(0, '#123456'), true, 'the generated CDS row is no rule row');
  assert.equal(committed.length, 1);
  assert.deepEqual(state.legendColorOverrides, { CDS: '#123456' }, 'a Legend-only color');
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

// U3a 1d (R14-8): a rule commit asks for the automatic rerender only when the
// displayed Result cannot show the Legend rows it regroups (OV-43 narrowed): a
// whole-row rename shows its row at once (the row the commit adds, the old row
// it retires), a rename of part of a row asks Python. History restores the
// rows with the rules, so Undo and Redo of the rename ask nothing either.
test('a whole-row Legend rename asks no rerender, an appended or partial one asks one, Undo and Redo ask none', async () => {
  const { resultCatalogFeatures } = await import('../../gbdraw/web/js/services/feature-catalog.js');
  const anchorProfile = { precision: 'exact', operator: 'single', partOrder: 'biological', strand: '+' };
  const ids = ['A', 'B'];
  const catalog = { schema: 5, items: [{
    resultIndex: 0, resultName: 'result-0.svg', recordKeys: ['REC1'],
    biologicalFeatures: ids.map((id, index) => ({
      recordKey: 'REC1', biologicalFeatureId: id, record_id: 'REC1', type: 'CDS', start: index * 100, end: index * 100 + 30,
      strand: 1, anchorProfile, qualifiers: { locus_tag: [id] }
    })),
    features: ids.map((id) => ({ svgId: `svg-${id}`, recordKey: 'REC1', biologicalFeatureId: id, fillColor: '#000000',
      drawnSelector: { hash: `svg-${id}`, location: null, recordLocation: null } })),
    orthogroups: [], annotations: [], comparisonMatches: []
  }] };
  // The displayed Result's Legend: the keys of the rows it shows.
  let shownKeys = [];
  const row = (key) => ({ getAttribute: (name) => (name === 'data-legend-key' ? key : null) });
  const group = { querySelectorAll: () => shownKeys.map(row) };
  const svg = { getElementById: (id) => (id === 'legend' ? { querySelector: (selector) => (selector === '#feature_legend' ? group : null) } : null) };
  const state = {
    featureCatalog: { value: catalog }, generatedMode: { value: 'circular' }, selectedResultIndex: { value: 0 },
    results: { value: [{ name: 'result-0.svg', content: '<svg />' }] }, svgContainer: { value: { querySelector: () => svg } },
    manualSpecificRules: [], featureColorOverrides: {}, featureOverrides: {}, featureVisibilityManualRules: [],
    svgResultIdentity: { value: 'result' }, fileLegendCaptions: { value: new Set() }, addedLegendCaptions: { value: new Set() },
    legendEntries: { value: [] }, deletedLegendEntries: { value: [] }, files: { t_color: null }, legendColorOverrides: {},
    extractedFeatures: { value: [] }
  };
  state.extractedFeatures.value = [...resultCatalogFeatures(state).renderedByIdentity.values()];
  const drawingState = withDrawings(state);
  let rerenders = 0;
  const actions = createFeatureRuleActions({
    ref: value => ({ value }), computed: get => ({ get value() { return get(); } }), state: drawingState,
    rulePreparation: createRulePreparation({ state: drawingState, evaluate: evaluatePythonRules }),
    runUndoableCheckpoint: async (_label, commit) => commit(), runUndoable: async (_label, commit) => commit(),
    prepareFileLegendEntries: async (intents) => {
      const diff = diffLegendIntents(state.legendEntries.value, intents);
      return { diff, isCurrent: () => true, apply: () => {
        state.legendEntries.value = intents;
        return { add: diff.add.map(({ caption, color }) => ({ caption, color })), retire: diff.remove.map(({ caption }) => caption) };
      } };
    },
    projectPaletteAndRules: () => true, ports: { requestAutomaticRerender: () => { rerenders += 1; return true; } }
  });
  const tagRule = (pattern, cap) => ({ feat: 'CDS', qual: 'locus_tag', val: `^(${pattern})$`, color: '#112233', cap });
  assert.equal(await actions.commitSpecificRules([tagRule('A|B', 'alpha')]), true);
  shownKeys = ['alpha'];
  const first = rerenders;
  const alpha = state.manualSpecificRules.map(rule => ({ ...rule }));
  // A rename in the Rules panel appends the new row, where Python keeps the
  // row in the old row's place; the Legend rename places it (OV-158).
  const rulesPanel = state.manualSpecificRules.map(rule => ({ ...rule }));
  assert.equal(await actions.commitSpecificRules([tagRule('A|B', 'beta')]), true);
  assert.equal(rerenders, first + 1, 'an appended relabeled row asks Python');
  state.manualSpecificRules.splice(0, state.manualSpecificRules.length, ...rulesPanel);
  state.legendEntries.value = [{ caption: 'alpha', color: '#112233' }];
  assert.equal(await actions.commitSpecificRules([tagRule('A|B', 'beta')], 'Rename legend item',
    { legendPlacement: { caption: 'beta', at: 'alpha' } }), true);
  assert.equal(rerenders, first + 1, 'the renamed row is shown at once in its place');
  shownKeys = ['beta'];
  const beta = state.manualSpecificRules.map(rule => ({ ...rule }));
  // Undo restores the rules and the rows they drew, then Redo.
  state.manualSpecificRules.splice(0, state.manualSpecificRules.length, ...alpha.map(rule => ({ ...rule })));
  shownKeys = ['alpha'];
  actions.followRestoredRules(beta);
  state.manualSpecificRules.splice(0, state.manualSpecificRules.length, ...beta.map(rule => ({ ...rule })));
  shownKeys = ['beta'];
  actions.followRestoredRules(alpha);
  assert.equal(rerenders, first + 1, 'Undo and Redo of the rename ask no rerender');
  assert.equal(await actions.commitSpecificRules([tagRule('A', 'gamma'), tagRule('B', 'beta')]), true);
  assert.equal(rerenders, first + 2, 'a rename of part of the row asks Python');
});
