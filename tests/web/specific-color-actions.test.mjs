import assert from 'node:assert/strict';
import { test } from 'node:test';
import { createFeatureRuleActions } from '../../gbdraw/web/js/app/feature-editor/rule-actions.js';
import { createRulePreparation } from '../../gbdraw/web/js/app/rule-matching.js';
import { createLegendManager } from '../../gbdraw/web/js/app/legend.js';
import { createFeatureColorActions } from '../../gbdraw/web/js/app/feature-editor/color-actions.js';
import { createHistoryManager } from '../../gbdraw/web/js/services/history.js';
import { compileDirectEditorMutationPlan } from '../../gbdraw/web/js/app/candidate-render.js';
import { readFileSync } from 'node:fs';
import {
  createRuleLegendCaptions, diffLegendIntents, generatedLegendRow, legendRowRules, prepareSpecificColorImport,
  rekeyedRuleRows, ruleCommitLegendRows, ruleLegendCaptions, savedRuleCaption, savedSpecificRules
} from '../../gbdraw/web/js/services/specific-color-rules.js';
import { legendColorIdentity, sameLegendLine } from '../../gbdraw/web/js/services/legend-row-names.js';
import { recordRuleMatches, recordVisibilityMatch, ruleKey } from '../../gbdraw/web/js/services/rule-matchers.js';
import { setFeatureVisibilityOverride } from '../../gbdraw/web/js/services/feature-visibility.js';
import { evaluatePythonRules, synthesizePythonColorRules } from './helpers/python-evaluator.mjs';
import { withDrawings, withModeDrawings } from './helpers/drawing-state.mjs';
import { FEATURE_EDITS_SLICE } from '../../gbdraw/web/js/services/drawing/feature-edits.js';
// A feature-edits slice over a draft map, so a writer writes the map.
const editsOf = (featureOverrides, featurePlacementOverrides) => FEATURE_EDITS_SLICE.create({ featureOverrides, featurePlacementOverrides }).handle;

// The Legend owner's writer of a row's own color (app/feature-editor.js wires it); these tests record its calls.
const legendRowColorWrites = [];
/** @param {unknown[]} write */
const recordLegendRowColor = (...write) => { legendRowColorWrites.push(write); };

/**
 * @param {typeof evaluatePythonRules} [evaluate]
 * @param {{ createHistory?: ((state: Record<string, any>) => Record<string, any>) | null } & Record<string, any>} [options]
 *   `createHistory`: the History manager the rule owner commits through (a recorder of its steps by default);
 *   the other options go to the rule owner (`askChoice`).
 */
const setup = (evaluate = evaluatePythonRules, { createHistory = null, linear = null, ...options } = {}) => {
  const state = {
    manualSpecificRules: [], extractedFeatures: { value: [
      { type: 'CDS', svg_id: 'a', qualifiers: { gene: ['a'] } },
      { type: 'CDS', svg_id: 'b', qualifiers: { gene: ['b'] } }
    ] },
    featureColorOverrides: {}, featureOverrides: {}, results: { value: [{name:'figure',content:'before'}] },
    svgResultIdentity: { value:'before' }, addedLegendCaptions: { value:new Set() },
    legendEntries: { value:[] }, files: {t_color:null}, legendColorOverrides: {},
    newSpecRule: {feat:'CDS',qual:'gene',val:'a',color:'#112233',cap:'Shared'}, labelReflowProcessing: { value:false },
    editableLabels: { value: [] }, manualPriorityRules: [],
    deletedLegendEntries: { value: [] }, dormantLegendEntries: { value: [] }, legendStrokeOverrides: {},
    specificRulePresets: [], selectedSpecificPreset: { value: '' }, specificRulePresetLoading: { value: false }
  };
  // The automatic rerender's busy flag, with the Vue `watch` the owner waits on.
  const watchers = new Set();
  const watch = (source, callback) => { const entry = {source, callback}; watchers.add(entry); return () => watchers.delete(entry); };
  const setRerendering = (busy) => {
    state.labelReflowProcessing.value = busy;
    [...watchers].filter(entry => entry.source === state.labelReflowProcessing).forEach(entry => entry.callback(busy));
  };
  const transactions = [], transactionScopes = [];
  // `linear`: a second drawing; `state.mode.value` picks the shown one.
  const view = linear ? withModeDrawings(Object.assign(state, { mode: { value: 'circular' } }), { circular: state, linear }) : withDrawings(state);
  const preparation = createRulePreparation({state: view, evaluate, synthesize: synthesizePythonColorRules});
  let prepareLegend = async () => {};
  let previousIntents = [];
  let lastRekeyed = null;
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
  const history = createHistory?.(state);
  const actions = createFeatureRuleActions({ref:value=>({value}),computed:get=>({get value(){return get();}}),state: view, rulePreparation:preparation,
    runUndoableCheckpoint: history?.runUndoableCheckpoint || transact('checkpoint'),
    runUndoable: history?.runUndoable || transact('intent'),
    prepareFileLegendEntries: async (intents, {isCurrent,previousFileIntents,rekeyed}) => {
    previousIntents=previousFileIntents;
    lastRekeyed=rekeyed;
    await prepareLegend(intents);
    if (!isCurrent()) return false;
    return { diff: diffLegendIntents(state.legendEntries.value, intents), isCurrent: () => true, apply: () => {
      legendApplies.push({ transaction: openTransaction, rules: state.manualSpecificRules.map(rule => rule.cap) });
      state.legendEntries.value=intents;
      state.results.value=[{name:'figure',content:'after'}];
    } };
  },
  // The Legend owner's writer, which the take-back of a drawing not shown calls (OV-1280, AUTH-ORIGIN-1);
  // `state` already has its drawing store, which a test may have made per mode.
  moveLegendRecords: (...args) => createLegendManager({ state, commitLegendRowRules: () => true }).moveLegendRecords(...args),
  projectPaletteAndRules:()=>true, ports:{requestAutomaticRerender:()=>{rerenders+=1;return true;},isRedrawOwed:()=>false}, nextTick:async()=>{}, watch, ...options});
  return {state,actions,history,preparation,transactions,transactionScopes,legendApplies, setLegendPreparation: fn => {prepareLegend=fn;}, previousIntents:()=>previousIntents, rekeyed:()=>lastRekeyed, rerenders:()=>rerenders, setRerendering};
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
  const pending=s.actions.commitSpecificRules(() => rules);
  await ready;
  assert.deepEqual(s.state.manualSpecificRules,[]);
  assert.equal(s.state.results.value[0].content,'before');
  assert.deepEqual(s.state.legendEntries.value,[]);
  release();
  assert.equal(await pending,true);
  s.setLegendPreparation(async () => {});
  // A caption written in two colors stays as written (PD-OI-042 scenario revision 2).
  assert.deepEqual(s.state.manualSpecificRules.map(r=>r.cap),['Shared','Shared']);
  assert.equal(s.state.manualSpecificRules[0].fromFile,true);
  assert.equal(s.transactions.length,1);
  assert.deepEqual(s.legendApplies,[{transaction:'Change specific color rules',rules:['Shared','Shared']}],
    'the legend rows apply after the rule transition, inside its History step (R13)');
  await s.actions.setSpecificRuleField(1,'cap','Renamed');
  assert.equal(s.state.manualSpecificRules[1].cap,'Renamed');
  assert.equal(s.state.manualSpecificRules[0].cap,'Shared');
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
    await s.actions.commitSpecificRules(() => rules);
    s.state.addedLegendCaptions.value.add('Independent manual legend');
    const before=JSON.stringify(s.state.manualSpecificRules), result=s.state.results.value;
    const count=s.transactions.length;
    let replacements=0;
    s.setLegendPreparation(async()=>{
      if(outcome==='error') throw new Error('measurement failed');
      // Every preparation goes stale, so no attempt may commit (OV-166 retries).
      s.state.svgResultIdentity.value=`replacement-${++replacements}`;
    });
    const candidate=rules.map(r=>({...r,cap:'Changed'}));
    if(outcome==='error') await assert.rejects(()=>s.actions.commitSpecificRules(() => candidate),/measurement failed/);
    else assert.equal(await s.actions.commitSpecificRules(() => candidate),false);
    assert.equal(JSON.stringify(s.state.manualSpecificRules),before);
    assert.equal(s.state.results.value,result);
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
  assert.equal(await s.actions.commitSpecificRules(() => rules),true);
  assert.equal(preparations,2,'one more preparation, against the bound Result');
  assert.deepEqual(s.state.manualSpecificRules.map(r=>r.cap),['Shared','Shared']);
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

// OV-1133 (S4-Q5): a pattern edit made while the automatic rerender runs
// waits for it. The rerender replaces the Result the edit's turn began
// against, so the edit prepares again against the new Result and applies; it
// is not dropped as stale.
test('a pattern edit made during the automatic rerender waits for it, then applies (OV-1133)', async () => {
  const s=setup();
  await s.actions.commitSpecificRules(() => [{feat:'CDS',qual:'gene',val:'a',color:'#112233',cap:'Shared'}]);
  const row=s.state.manualSpecificRules[0];
  s.setRerendering(true);
  const pending=s.actions.setSpecificRuleField(0,'val','b');
  await new Promise(resolve => setImmediate(resolve));
  assert.equal(row.val,'a','not before the rerender is idle');
  // The rerender lands: its Result, then idle.
  s.state.svgResultIdentity.value='rerendered';
  s.state.results.value=[{name:'figure',content:'rerendered'}];
  s.setRerendering(false);
  assert.equal(await pending,true);
  assert.equal(row.val,'b');
  assert.equal(s.actions.specificRulePatternDraft(row),null,'no draft stays');
  assert.equal(s.transactions.length,2);
});

// The edit that waits is the edit of its Session: a Session operation made
// while it waits (a Load, Reset Settings, an Undo across a Load, or one still
// running when the rerender ends) supersedes it, and nothing is committed onto
// the other Session.
for (const operation of ['Load', 'Reset Settings', 'Undo across a Load', 'Session pending']) {
  test(`a pattern edit waiting for the automatic rerender is superseded by ${operation}`, async () => {
    const s=setup();
    await s.actions.commitSpecificRules(() => [{feat:'CDS',qual:'gene',val:'a',color:'#112233',cap:'Shared'}]);
    s.setRerendering(true);
    const pending=s.actions.setSpecificRuleField(0,'val','b');
    await new Promise(resolve => setImmediate(resolve));
    const replace = (rows) => s.state.manualSpecificRules.splice(0, s.state.manualSpecificRules.length, ...rows);
    if (operation === 'Load') {
      replace([{feat:'CDS',qual:'gene',val:'z',color:'#445566',cap:'Loaded'}]);
      s.actions.clearSpecificRulePatternDrafts();
    }
    if (operation === 'Reset Settings') { s.actions.clearSpecificRulePatternDrafts(); replace([]); }
    // History restores copies of the rows it saved.
    if (operation === 'Undo across a Load') replace(s.state.manualSpecificRules.map(rule => ({...rule})));
    if (operation === 'Session pending') s.state.sessionOperationAvailability = () => ({ status: 'busy' });
    const after=JSON.stringify(s.state.manualSpecificRules);
    s.state.svgResultIdentity.value='rerendered';
    s.state.results.value=[{name:'figure',content:'rerendered'}];
    s.setRerendering(false);
    const outcome=await pending;
    assert.notEqual(outcome,true);
    assert.equal(JSON.stringify(s.state.manualSpecificRules),after,'nothing is committed onto the other Session');
    assert.equal(s.transactions.length,1);
  });
}

test('a caption written in two colors keeps a row per color, each under its own key, and stays as written',async()=>{
  const s=setup();s.state.manualSpecificRules.push(...rules.map(rule=>({...rule})));
  await s.actions.setSpecificRuleField(0,'val','a');
  // The rows the rules had: one caption in two colors is two rows, each under its own key.
  assert.deepEqual(s.previousIntents(),[{caption:'Shared [#112233]',color:'#112233'},{caption:'Shared [#445566]',color:'#445566'}]);
  assert.deepEqual(s.state.manualSpecificRules.map(rule=>rule.cap),['Shared','Shared']);
});

test('complete caption recolor uses bounded intent while caption replacement retains checkpoint History', async () => {
  const s = setup();
  s.state.legendEntries.value = [{ caption: 'Shared', color: '#012345' }];
  const recolor = [{ feat: 'CDS', qual: 'gene', val: 'a', color: '#abcdef', cap: 'Shared' }];
  assert.equal(await s.actions.commitSpecificRules(() => recolor, 'Recolor', {
    previousLegendIntents: [{ caption: 'Shared', color: '#012345' }]
  }), true);
  assert.deepEqual(s.transactionScopes, ['intent']);
  assert.deepEqual(s.state.legendEntries.value, [{ caption: 'Shared', color: '#abcdef' }]);
  assert.deepEqual(s.state.manualSpecificRules.map(rule => rule.cap), ['Shared']);
  assert.equal(s.transactions.length, 1);
  assert.equal(await s.actions.commitSpecificRules(() => recolor.map(rule => ({ ...rule, cap: 'Renamed' }))), true);
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
  assert.equal(await s.actions.commitSpecificRules(() => [
    { feat: 'CDS', qual: 'gene', val: 'a', color: '#112233', cap: 'Shared' },
    { feat: 'CDS', qual: 'gene', val: 'b', color: '#445566', cap: 'Shared' },
    { feat: 'CDS', qual: 'gene', val: 'c', color: '#778899', cap: 'Shared' },
    { feat: 'CDS', qual: 'gene', val: 'absent', color: '#aabbcc', cap: 'Shared' }
  ]), true);
  assert.deepEqual(s.state.manualSpecificRules.map(rule => rule.cap), ['Shared', 'Shared', 'Shared', 'Shared']);
  assert.deepEqual(intents, [{ caption: 'Shared [#112233]', color: '#112233', name: 'Shared' }]);
});

// F-2 (release D-39): a hash rule names the feature's hash in its source
// record (`selector.hash`), which Generate matches on a cropped or
// reverse-complemented record too, and which its rendered ID carries. The rule
// "This feature only" writes matches that feature live too.
test('a This feature only rule on a cropped record paints its feature live as Generate does', async () => {
  const s = setup();
  const cropped = {
    type: 'CDS', svg_id: 'f3ccacda4_record_1', stable_feature_id: 'f3ccacda4', record_id: 'TESTA',
    selector: { type: 'CDS', start: 300, end: 600, strand: '+', hash: 'f3ccacda4', qualifiers: { locus_tag: ['TESTA_0001'] } },
    qualifiers: { locus_tag: ['TESTA_0001'] }
  };
  s.state.extractedFeatures.value = [cropped];
  const qualifier = s.actions.getFeatureQualifier(cropped);
  assert.deepEqual(qualifier, { qual: 'hash', val: 'f3ccacda4' });
  let intents = null;
  s.setLegendPreparation(async (next) => { intents = next; });
  assert.equal(await s.actions.commitSpecificRules(() => [
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
  const element = (attributes, children = [], localName = 'g') => ({
    localName,
    children,
    getAttribute: (name) => attributes[name] ?? null,
    hasAttribute: (name) => name in attributes,
    querySelector: () => null,
    querySelectorAll: () => children
  });
  const legend = element({}, rows.map(([key, fill, attributes = {}]) => element({ 'data-legend-key': key, ...attributes }, [element({ fill }, [], 'path')])));
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
  assert.equal(await s.actions.commitSpecificRules(() => [rule]), true);
  // The row's key is never drawn: the row is named by the rule's caption.
  assert.deepEqual(intents, [{ caption: 'CDS [#ff0000]', color: '#ff0000', name: 'CDS' }]);
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
  assert.equal(await s.actions.commitSpecificRules(() => [{ ...rule, color: '#00ff00' }]), true);
  assert.deepEqual(s.previousIntents(), [{ caption: 'CDS [#ff0000]', color: '#ff0000' }]);
  assert.deepEqual(intents, [{ caption: 'CDS [#00ff00]', color: '#00ff00', name: 'CDS' }]);
});

// A caption written in two colors is two rows, each under its own key. A
// feature belongs to the row of its rule, and a feature with only a color edit
// of that caption and color (its rule match not known) to the same row.
test('a feature of a caption written in two colors belongs to the row of its color', () => {
  const s = setup();
  s.state.manualSpecificRules.push(...rules.map((rule) => ({ ...rule })));
  const [a, b] = s.state.extractedFeatures.value;
  const edited = { type: 'CDS', svg_id: 'c', qualifiers: { gene: ['c'] } };
  s.state.extractedFeatures.value.push(edited);
  const keys = s.state.manualSpecificRules.map(ruleKey);
  const matches = [[0], [1], []];
  recordRuleMatches([a, b, edited], keys, (index) => ({ matched: matches[index], priorities: matches[index].map(() => 0), declined: [] }));
  s.state.featureColorOverrides.c = { color: '#445566', caption: 'Shared' };
  s.state.originalLegendOrder = { value: ['Shared [#112233]', 'Shared [#445566]'] };
  s.state.svgContainer = mounted([['Shared [#112233]', '#112233'], ['Shared [#445566]', '#445566']]);
  const captionOf = s.actions.effectiveLegendCaptions();
  assert.deepEqual([a, b, edited].map(captionOf), ['Shared [#112233]', 'Shared [#445566]', 'Shared [#445566]']);
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
  assert.equal(await s.actions.commitSpecificRules(() => [rule]), true);
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
  assert.equal(await s.actions.commitSpecificRules(() => [rule]), true);
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
  setFeatureVisibilityOverride(editsOf(s.state.featureOverrides), hidden, 'off');
  Object.assign(s.state, {
    displayedResultMetadata: () => ({ renderedFeatureIdentities: { renderedIds: new Set(['a1', 'a2']) } }),
    originalLegendOrder: { value: ['CDS'] }, svgContainer: mounted([['CDS', '#54bcf8']])
  });
  s.state.legendEntries.value = [{ caption: 'CDS', originalCaption: 'CDS', color: '#54bcf8' }];
  s.setLegendPreparation(async () => {});
  const rule = { feat: 'CDS', qual: 'locus_tag', val: '^(FL2|A2)$', color: '#2266aa', cap: 'Group' };
  assert.equal(await s.actions.commitSpecificRules(() => [rule]), true);
  assert.deepEqual(s.previousIntents(), [], 'Python keeps the CDS row of the displayed Result');
});

// OV-344: the rows a rule commit draws are those of the displayed Result. A
// rule whose features only another batch Result draws, or that a Feature
// visibility rule hides, draws no row on it, as Generate draws none.
test('a rule commit draws no row for a rule that only features the displayed Result does not draw match (OV-344)', async () => {
  const s = setup();
  const shown = { type: 'CDS', svg_id: 'b1', qualifiers: { locus_tag: ['B1'] } };
  const ruleHidden = { type: 'CDS', svg_id: 'b2', qualifiers: { locus_tag: ['B2'] } };
  const otherResult = { type: 'CDS', svg_id: 'a2', qualifiers: { locus_tag: ['A2'] } };
  s.state.extractedFeatures.value = [otherResult, shown, ruleHidden];
  s.state.results.value = [{ name: 'A', content: 'a' }, { name: 'B', content: 'b' }];
  const hideB2 = { recordId: '*', featureType: 'CDS', qualifier: 'locus_tag', value: '^B2$', action: 'off' };
  s.state.featureVisibilityManualRules = [hideB2];
  // Python's matches of the visibility rule (R4).
  [shown, ruleHidden, otherResult].forEach((feature) => recordVisibilityMatch(feature, hideB2, { matches: feature === ruleHidden }));
  s.state.displayedResultMetadata = () => ({ renderedFeatureIdentities: { renderedIds: new Set(['b1', 'b2']) } });
  let intents = null;
  s.setLegendPreparation(async (next) => { intents = next; });
  assert.equal(await s.actions.commitSpecificRules(() => [
    { feat: 'CDS', qual: 'locus_tag', val: '^A2$', color: '#7b2cbf', cap: 'other Result' },
    { feat: 'CDS', qual: 'locus_tag', val: '^B2$', color: '#264653', cap: 'hidden' },
    { feat: 'CDS', qual: 'locus_tag', val: '^B1$', color: '#c83366', cap: 'shown' }
  ]), true);
  assert.deepEqual(intents, [{ caption: 'shown', color: '#c83366' }]);
});

// L1: a rule captioned with a track's row (Python's own rows,
// `_generated_legend_rows`) is no rule of that row, also in the track's
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

// The Python rows the Web app reads from their keys are those Python
// generates (tests/test_legend_row_facts.py): a present type's row with the
// palette key of its fill, or a track's row, whose fill the track rows of the
// request give (tests/web/track-legend-row-colors.test.mjs).
test('the generated Legend rows read from their keys are Python\'s', () => {
  const { cases } = JSON.parse(readFileSync(new URL('../fixtures/legend_generated_rows.json', import.meta.url), 'utf8'));
  for (const { name, features_present: types, rows } of cases) {
    for (const [key, paletteKey] of rows) {
      const row = generatedLegendRow(key);
      assert.equal(row.track, !types.includes(paletteKey), `${name}: ${key} is a track row or a present type's row`);
      assert.equal(row.paletteKey, row.track ? '' : paletteKey, `${name}: ${key}`);
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
  const legend = createLegendManager({ state: withDrawings(state), commitLegendRowRules: (build, label) => {
    committed.push({ rules: build(state.manualSpecificRules, (row, patch) => ({ ...row, ...patch })), label });
    return true;
  }, readShownLegendColor: (entry) => entry?.color });
  assert.equal(legend.updateLegendEntryColor(1, '#00ff00'), true);
  assert.deepEqual(committed, [{ rules: [{ ...rule, color: '#00ff00' }], label: 'Change legend color' }]);
  assert.deepEqual([legend.legendRowHasRules(0), legend.legendRowHasRules(1)], [false, true]);
  assert.equal(legend.updateLegendEntryColor(0, '#123456'), true, 'the generated CDS row is no rule row');
  assert.equal(committed.length, 1);
  assert.deepEqual(state.legendColorOverrides, { CDS: '#123456' }, 'a Legend-only color');
});

// D-35: a Legend color edit of a row that merges a rule's row (one name, one
// swatch, D-28) revises that rule and writes the Legend color of the other
// keys in the rule commit's one step, so the row stays one row.
test('the Legend editor recolors every row a merged row stands for, a rule row through its rule, in one step (D-35)', () => {
  const ref = (value) => ({ value });
  const rule = { feat: 'tRNA', qual: 'product', val: '.*', color: '#ff0000', cap: 'CDS' };
  const state = {
    manualSpecificRules: [rule],
    // Python drew CDS in its palette color; its Legend color makes it the rule's swatch.
    svgContainer: mounted([['CDS', '#54bcf8'], ['CDS [#ff0000]', '#ff0000', { display: 'none', 'data-legend-merged': 'CDS' }]]),
    results: ref([]), selectedResultIndex: ref(0),
    legendEntries: ref([
      { caption: 'CDS', originalCaption: 'CDS', color: '#ff0000' },
      { caption: 'CDS [#ff0000]', originalCaption: 'CDS [#ff0000]', color: '#ff0000' }
    ]),
    originalLegendOrder: ref(['CDS', 'CDS [#ff0000]']), deletedLegendEntries: ref([]), originalLegendColors: ref({}),
    legendStrokeOverrides: {}, legendColorOverrides: { CDS: '#ff0000' }, adv: {}
  };
  const committed = [];
  const legend = createLegendManager({ state: withDrawings(state), commitLegendRowRules: (build, label, options) => {
    committed.push({ rules: build(state.manualSpecificRules, (row, patch) => ({ ...row, ...patch })), label });
    options?.afterCommit?.([]);
    return true;
  }, readShownLegendColor: (entry) => entry?.color });
  assert.equal(legend.legendRowHasRules(0), true, 'the shown row stands for a rule row');
  assert.equal(legend.updateLegendEntryColor(0, '#00ff00'), true);
  assert.deepEqual(committed, [{ rules: [{ ...rule, color: '#00ff00' }], label: 'Change legend color' }]);
  assert.deepEqual(state.legendColorOverrides, { CDS: '#00ff00' }, 'the other key\'s Legend color, in the rule commit');
  assert.equal(state.legendEntries.value[0].color, '#00ff00');
});

// Review N-1 of #1051 (D-35): the feature popup's "Legend item" scope edits the
// row the Legend draws, so on a merged row it reaches every row the row stands
// for, through the Legend owner's writers, and the row stays one row.
test('the popup\'s Legend item scope recolors and restrokes every row a merged row stands for (D-35, review N-1)', async () => {
  const ref = (value) => ({ value });
  const cds = { id: 'c1', svg_id: 'c1', type: 'CDS', qualifiers: {} };
  // CDS renamed `Shared` in the color of the tRNA rule row `Shared`: one row.
  const fixture = (defaultColorType) => {
    const rule = { feat: 'tRNA', qual: 'product', val: '.*', color: '#445566', cap: 'Shared' };
    const state = {
      manualSpecificRules: [rule],
      svgContainer: mounted([['CDS', '#445566'], ['Shared', '#445566', { display: 'none', 'data-legend-merged': 'CDS' }]]),
      results: ref([]), selectedResultIndex: ref(0), appliedPaletteColors: ref({ CDS: '#445566' }),
      extractedFeatures: ref([cds]), biologicalFeatures: ref([cds]), clickedFeature: ref(null),
      featureStyleScopeDialog: { show: true, kind: 'fill', feat: cds, color: '#ff0000', strokeColor: null, strokeWidth: null,
        legendName: 'CDS', matchingRule: null, siblingCount: 1, ruleMatchCount: 0, defaultColorType },
      resetColorDialog: {}, legendRenameDialog: {},
      legendEntries: ref([
        { caption: 'Shared', originalCaption: 'CDS', color: '#445566' },
        { caption: 'Shared', originalCaption: 'Shared', color: '#445566' }
      ]),
      originalLegendOrder: ref(['CDS', 'Shared']), deletedLegendEntries: ref([]), dormantLegendEntries: ref([]),
      originalLegendColors: ref({}), featureColorOverrides: {}, featureStrokeOverrides: {},
      legendStrokeOverrides: {}, legendColorOverrides: { CDS: '#445566' }, adv: {}
    };
    const committed = [];
    const defaults = [];
    const legend = createLegendManager({ state: withDrawings(state), commitLegendRowRules: () => assert.fail('no Legend editor commit'),
      readShownLegendColor: (entry) => entry?.color });
    const popup = createFeatureColorActions({
      state: withDrawings(state), readDrawnTracks: () => [], getFeatureElements: () => [], getFeatureFillElements: () => [],
      readUserDefaultColor: () => null, setDefaultColor: (_drawing, type, color) => { defaults.push([type, color]); },
      writeLegendRowColor: legend.writeLegendRowColor,
      ruleActions: {
        commitSpecificRules: async (build, label, options = {}) => {
          state.manualSpecificRules = build(state.manualSpecificRules, (row, patch) => ({ ...row, ...patch }));
          committed.push(label);
          options.afterCommit?.([]);
          return true;
        },
        getLegendRowRules: (key) => (key === 'Shared' ? state.manualSpecificRules : []),
        effectiveLegendCaptions: () => () => 'CDS',
        findFeaturesWithSameLegendItem: () => [], findFeaturesWithSameDisplayedLabel: () => [],
        findFeaturesWithSameIndividualLabel: () => [], getDisplayedFeatureLabel: () => '', getIndividualFeatureLabel: () => '',
        getFeatureQualifier: (feature) => ({ qual: 'hash', val: feature.svg_id }), groupColorRules: () => assert.fail('a merged rule row is recolored through its rule'),
        runWithRuleMatches: (_rules, commit) => commit()
      }
    });
    return { state, popup, committed, defaults };
  };
  for (const defaultColorType of [null, 'CDS']) {
    const { state, popup, committed, defaults } = fixture(defaultColorType);
    await popup.handleColorScopeChoice('caption');
    assert.deepEqual(state.manualSpecificRules.map(({ cap, color }) => [cap, color]), [['Shared', '#ff0000']],
      `the merged rule row's rule (${defaultColorType})`);
    assert.equal(committed.length, 1, 'one rule commit');
    if (defaultColorType) assert.deepEqual(defaults, [['CDS', '#ff0000']], 'the palette row\'s default color');
    else assert.deepEqual(state.legendColorOverrides, { CDS: '#ff0000' }, 'the other key\'s Legend color, in the rule commit');
  }
  // One route sets a type's default color: the Legend popup (`setTypeDefaultColor`) ends where the scope choice does,
  // the row's own Legend color released through the Legend owner.
  const viaScope = fixture('CDS');
  await viaScope.popup.handleColorScopeChoice('caption');
  const viaPopup = fixture('CDS');
  await viaPopup.popup.setTypeDefaultColor('CDS', '#ff0000');
  const outcome = ({ state, committed, defaults }) => ({
    rules: state.manualSpecificRules, overrides: { ...state.legendColorOverrides }, committed, defaults
  });
  assert.deepEqual(outcome(viaPopup), outcome(viaScope));
  assert.equal('CDS' in viaScope.state.legendColorOverrides, false, 'the row follows its type\'s default color again');
  const { state, popup } = fixture(null);
  Object.assign(state.featureStyleScopeDialog, { kind: 'stroke', color: null, strokeColor: '#00aa00', strokeWidth: 2 });
  assert.equal(await popup.handleFeatureStyleScopeChoice('caption'), true);
  assert.deepEqual(Object.keys(state.legendStrokeOverrides).sort(), ['CDS', 'Shared'], 'the stroke of every key');
  assert.deepEqual([state.legendStrokeOverrides.Shared.strokeColor, state.legendStrokeOverrides.Shared.strokeWidth], ['#00aa00', 2]);
});

test('a color action prepares the color matches of its rules without the caption evaluation, and a prepared table answers at once', async () => {
  const kinds = [];
  const s = setup(async (payload, options) => { kinds.push(payload.kind); return evaluatePythonRules(payload, options); });
  const commits = [];
  const pending = s.actions.runWithRuleMatches(rules, () => { commits.push('cold'); return 'cold'; });
  assert.ok(pending instanceof Promise, 'unprepared rules are prepared first');
  assert.deepEqual(commits, []);
  assert.equal(await pending, 'cold');
  assert.deepEqual(kinds, ['color'], 'the matches only');
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
  assert.equal(await s.actions.commitSpecificRules(() => [{ ...rule }]), true);
  assert.equal(s.state.legendColorOverrides.Shared, '#112233', 'the same color keeps the row\'s Legend color');
  assert.equal(await s.actions.commitSpecificRules(() => [{ ...rule, color: '#445566' }]), true);
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
  assert.equal(await s.actions.commitSpecificRules(() => [{ ...typeRule }, { ...typeRule, val: 'b' }]), true);
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
    svgResultIdentity: { value: 'result' }, addedLegendCaptions: { value: new Set() },
    legendEntries: { value: [] }, deletedLegendEntries: { value: [] }, files: { t_color: null }, legendColorOverrides: {},
    extractedFeatures: { value: [] }
  };
  state.extractedFeatures.value = [...resultCatalogFeatures(state).renderedByIdentity.values()];
  const drawingState = withDrawings(state);
  let rerenders = 0;
  const actions = createFeatureRuleActions({
    ref: value => ({ value }), computed: get => ({ get value() { return get(); } }), state: drawingState,
    rulePreparation: createRulePreparation({ state: drawingState, evaluate: evaluatePythonRules, synthesize: synthesizePythonColorRules }),
    runUndoableCheckpoint: async (_label, commit) => commit(), runUndoable: async (_label, commit) => commit(),
    prepareFileLegendEntries: async (intents) => {
      const diff = diffLegendIntents(state.legendEntries.value, intents);
      return { diff, isCurrent: () => true, apply: () => {
        state.legendEntries.value = intents;
        return { add: diff.add.map(({ caption, color }) => ({ caption, color })), retire: diff.remove.map(({ caption }) => caption) };
      } };
    },
    projectPaletteAndRules: () => true, ports: { requestAutomaticRerender: () => { rerenders += 1; return true; }, isRedrawOwed: () => false }
  });
  const tagRule = (pattern, cap) => ({ feat: 'CDS', qual: 'locus_tag', val: `^(${pattern})$`, color: '#112233', cap });
  assert.equal(await actions.commitSpecificRules(() => [tagRule('A|B', 'alpha')]), true);
  shownKeys = ['alpha'];
  const first = rerenders;
  const alpha = state.manualSpecificRules.map(rule => ({ ...rule }));
  // A rename in the Rules panel appends the new row, where Python keeps the
  // row in the old row's place; the Legend rename places it (OV-158).
  const rulesPanel = state.manualSpecificRules.map(rule => ({ ...rule }));
  assert.equal(await actions.commitSpecificRules(() => [tagRule('A|B', 'beta')]), true);
  assert.equal(rerenders, first + 1, 'an appended relabeled row asks Python');
  state.manualSpecificRules.splice(0, state.manualSpecificRules.length, ...rulesPanel);
  state.legendEntries.value = [{ caption: 'alpha', color: '#112233' }];
  assert.equal(await actions.commitSpecificRules(() => [tagRule('A|B', 'beta')], 'Rename legend item',
    { legendPlacement: { name: 'beta', color: '#112233', at: 'alpha', replaces: true } }), true);
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
  assert.equal(await actions.commitSpecificRules(() => [tagRule('A', 'gamma'), tagRule('B', 'beta')]), true);
  assert.equal(rerenders, first + 2, 'a rename of part of the row asks Python');
});

// OV-445, OV-443: a rule commit runs once its rules are prepared, and a newer
// commit supersedes a pending one. An Add is never dropped that way: a second
// Add, or an edit made while an Add is pending, waits for it and applies to
// the table that holds its rule; an Add waits for a pending edit the same way.
// Each write is its own History step (R11).
/** @param {Parameters<typeof setup>[1]} [options] */
const queuedRuleWrites = (options) => {
  const s=setup(undefined, options);
  s.state.errorLog={value:null};
  s.state.manualSpecificRules.push(
    {feat:'CDS',qual:'gene',val:'b',color:'#445566',cap:''},
    {feat:'CDS',qual:'gene',val:'zz',color:'#778899',cap:''});
  // Holds the first commit's Legend preparation; later commits run at once.
  let release, started;
  const ready=new Promise(resolve=>{started=resolve;});
  s.setLegendPreparation(()=>new Promise(resolve=>{release=resolve;started();}));
  const held=async(pending)=>{ await ready; s.setLegendPreparation(async()=>{}); return [pending,()=>release()]; };
  const table=()=>s.state.manualSpecificRules.map(rule=>`${rule.val} ${rule.color} ${rule.cap}`);
  return {s,held,table};
};
const withDocument=async(run)=>{
  const previousDocument=globalThis.document;
  globalThis.document={querySelector:()=>null};
  try { await run(); } finally { globalThis.document=previousDocument; }
};
const added='a #112233 Shared';
const ruleEdits={
  color:{run:s=>s.actions.setSpecificRuleField(0,'color','#000000'),label:'Edit specific color rule',
    table:['b #000000 ','zz #778899 ']},
  caption:{run:s=>s.actions.setSpecificRuleField(0,'cap','Edited'),label:'Edit specific color rule',
    table:['b #445566 Edited','zz #778899 ']},
  pattern:{run:s=>s.actions.setSpecificRuleField(0,'val','b|c'),label:'Edit specific color rule',
    table:['b|c #445566 ','zz #778899 ']},
  move:{run:s=>s.actions.moveSpecificRuleUp(1),label:'Move specific color rule',table:['zz #778899 ','b #445566 ']},
  remove:{run:s=>s.actions.removeSpecificRule(0),label:'Remove specific color rule',table:['zz #778899 ']}
};
for (const second of ['a', 'b']) {
  test(`a second Add while the first is pending adds both rules, one History step each (OV-445, pattern ${second})`, () => withDocument(async () => {
    const {s,held,table}=queuedRuleWrites();
    const [first,release]=await held(s.actions.addSpecificRule());
    s.state.newSpecRule.val=second;
    const next=s.actions.addSpecificRule();
    release();
    assert.deepEqual(await Promise.all([first,next]),[true,true]);
    assert.deepEqual(table(),['b #445566 ','zz #778899 ',added,`${second} #112233 Shared`]);
    assert.deepEqual(s.transactions,['Add specific color rule','Add specific color rule']);
    assert.equal(s.state.newSpecRule.val,'','the pattern of the last Add is cleared');
  }));
}
// An Add whose pattern may match slowly asks first. Its write takes its turn
// at the answer, so it also applies to the table an earlier Add wrote.
test('an Add that asks about a slow pattern applies after a pending Add, to the table with its rule', () => withDocument(async () => {
  let answer;
  const {s,held,table}=queuedRuleWrites({askChoice:()=>new Promise(resolve=>{answer=resolve;})});
  const [first,release]=await held(s.actions.addSpecificRule());
  s.state.newSpecRule.val='(a+)+';
  const asked=s.actions.addSpecificRule();
  assert.deepEqual(table(),['b #445566 ','zz #778899 '],'nothing is written before the answer');
  answer('add');
  release();
  assert.deepEqual(await Promise.all([first,asked]),[true,true]);
  assert.deepEqual(table(),['b #445566 ','zz #778899 ',added,'(a+)+ #112233 Shared']);
  assert.deepEqual(s.transactions,['Add specific color rule','Add specific color rule']);
}));
// The question is asked before the Add's History step (the real manager): the
// step is recorded once, at the answer, under the Add's name, and "Edit
// pattern" records none. An edit made while the question is asked is a step
// of its own, before the Add's.
test('an Add that asks about a slow pattern is one History step recorded at the answer; Edit pattern records none', () => withDocument(async () => {
  let answer;
  const s=setup(undefined,{askChoice:()=>new Promise(resolve=>{answer=resolve;}),createHistory:state=>createHistoryManager({
    buildIntent:()=>structuredClone({rules:state.manualSpecificRules}),
    applyIntent:({rules})=>{state.manualSpecificRules.splice(0,Infinity,...structuredClone(rules));},
    buildCheckpoint:()=>({}),applyCheckpoint:()=>{}})});
  s.state.errorLog={value:null};
  s.state.manualSpecificRules.push({feat:'CDS',qual:'gene',val:'b',color:'#445566',cap:''});
  Object.assign(s.state.newSpecRule,{val:'(a+)+',cap:''});
  const table=()=>s.state.manualSpecificRules.map(rule=>`${rule.val} ${rule.color}`);
  await s.history.initializeIntentBaseline();
  const ask=async()=>{
    const asked=s.actions.addSpecificRule();
    await new Promise(resolve=>setImmediate(resolve));
    assert.equal(s.history.mutationPending(),false,'no step is open while the question is asked');
    return [asked];
  };
  let [asked]=await ask();
  answer('cancel');
  await asked;
  assert.deepEqual([table(),s.history.getUndoCount()],[['b #445566'],0]);

  [asked]=await ask();
  assert.equal(await s.actions.setSpecificRuleField(0,'color','#000000'),true);
  assert.deepEqual([s.history.getUndoCount(),s.history.undoLabel()],[1,'Edit specific color rule']);
  answer('add');
  assert.equal(await asked,true);
  assert.deepEqual([table(),s.history.getUndoCount(),s.history.undoLabel()],[['b #000000','(a+)+ #112233'],2,'Add specific color rule']);
  await s.history.undo();
  assert.deepEqual([table(),s.history.undoLabel()],[['b #000000'],'Edit specific color rule']);
  await s.history.redo();
  assert.deepEqual(table(),['b #000000','(a+)+ #112233']);
}));
for (const [kind, edit] of Object.entries(ruleEdits)) {
  test(`a ${kind} edit made while an Add is pending applies after it, to the table with its rule (OV-443)`, () => withDocument(async () => {
    const {s,held,table}=queuedRuleWrites();
    const [add,release]=await held(s.actions.addSpecificRule());
    const edited=edit.run(s);
    release();
    assert.deepEqual(await Promise.all([add,edited]),[true,true]);
    assert.deepEqual(table(),[...edit.table,added]);
    assert.deepEqual(s.transactions,['Add specific color rule',edit.label]);
  }));
  if (kind === 'pattern') continue;
  test(`an Add made while a ${kind} edit is pending applies after it (OV-443)`, () => withDocument(async () => {
    const {s,held,table}=queuedRuleWrites();
    const [edited,release]=await held(edit.run(s));
    const add=s.actions.addSpecificRule();
    release();
    assert.deepEqual(await Promise.all([edited,add]),[true,true]);
    assert.deepEqual(table(),[...edit.table,added]);
    assert.deepEqual(s.transactions,[edit.label,'Add specific color rule']);
  }));
}
// OV-448: a pending pattern edit is the one write a newer write cancels. Its
// draft stays, and the newer write does not wait for it, however long its
// evaluation takes.
const afterHeldPattern={
  add:{run:s=>s.actions.addSpecificRule(),label:'Add specific color rule',table:['b #445566 ','zz #778899 ',added]},
  color:{run:s=>s.actions.setSpecificRuleField(1,'color','#000000'),label:'Edit specific color rule',table:['b #445566 ','zz #000000 ']},
  move:ruleEdits.move,
  remove:{run:s=>s.actions.removeSpecificRule(1),label:'Remove specific color rule',table:['b #445566 ']},
  clear:{run:s=>s.actions.clearAllSpecificRules(),label:'Clear specific color rules',table:[]}
};
for (const [kind, later] of Object.entries(afterHeldPattern)) {
  test(`a later ${kind} does not wait for a held pattern edit, which stays as a draft (OV-448)`, () => withDocument(async () => {
    const {s,held,table}=queuedRuleWrites();
    const row=s.state.manualSpecificRules[0];
    const [edited,release]=await held(ruleEdits.pattern.run(s));
    assert.equal(await later.run(s),true,'while the pattern edit is still held');
    assert.deepEqual(table(),later.table);
    release();
    assert.equal(await edited,false);
    assert.deepEqual(table(),later.table);
    assert.deepEqual(s.transactions,[later.label]);
    if (kind !== 'clear') assert.deepEqual([s.actions.specificRulePattern(row),s.actions.specificRulePatternDraft(row).pending],['b|c',false]);
  }));
}

// OV-444: every writer of the rule table, inside the rule owner or not,
// builds its rules from the table at its turn, so a write never drops a
// pending Add. A pattern keystroke or a drawer close cancels only pattern
// commits, never a pending color edit.
const popupFor=(s,ports={})=>{
  Object.assign(s.state,{clickedFeature:{value:null},svgContainer:{value:null},biologicalFeatures:{value:[]},
    featureStyleScopeDialog:{},resetColorDialog:{},legendRenameDialog:{},originalLegendOrder:{value:[]},
    originalLegendColors:{value:{}},deletedLegendEntries:{value:[]},legendStrokeOverrides:{},featureStrokeOverrides:{},
    selectedResultIndex:{value:0},appliedPaletteColors:{value:{CDS:'#cccccc'}}});
  return createFeatureColorActions({state:s.state,ruleActions:s.actions,getFeatureElements:()=>[],getFeatureFillElements:()=>[],
    readUserDefaultColor:()=>null,setDefaultColor:()=>{},writeLegendRowColor:recordLegendRowColor,...ports});
};
test('a popup recolor of a rule made while an Add is pending applies after it (OV-444)', () => withDocument(async () => {
  const {s,held,table}=queuedRuleWrites();
  const popup=popupFor(s);
  const [add,release]=await held(s.actions.addSpecificRule());
  Object.assign(s.state.featureStyleScopeDialog,{show:true,feat:s.state.extractedFeatures.value[1],color:'#000000',
    matchingRule:s.state.manualSpecificRules[0]});
  const recolored=popup.handleColorScopeChoice('rule');
  release();
  await Promise.all([add,recolored]);
  assert.deepEqual(table(),['b #000000 ','zz #778899 ',added]);
  assert.deepEqual(s.transactions,['Add specific color rule','Change specific color rule']);
}));
// OV-460: a color action's rule write is a second `runUndoable` call, made
// after the Python preparation, on the step its action runs on. It joins that
// step (`joinsOpenStep`, set by color-actions.js and passed to History by
// rule-actions.js): one step, under the name the step had when the action
// started, here the color input's.
test('a popup recolor of a rule and its rule write are one History step under the step\'s name', () => withDocument(async () => {
  const {s,held,table}=queuedRuleWrites({createHistory:state=>createHistoryManager({
    buildIntent:()=>JSON.parse(JSON.stringify(state.manualSpecificRules)),
    applyIntent:()=>{},buildCheckpoint:()=>({}),applyCheckpoint:()=>{}})});
  await s.history.initializeIntentBaseline();
  const popup=popupFor(s);
  Object.assign(s.state.featureStyleScopeDialog,{show:true,feat:s.state.extractedFeatures.value[1],color:'#000000',
    matchingRule:s.state.manualSpecificRules[0]});
  // app/history-inputs.js: the control begins the step. app-setup.js: the action runs on it.
  await s.history.begin('Change color',{source:'input-adapter',owner:{tagName:'INPUT'}});
  const [recolored,release]=await held(s.history.runUndoable('Change feature color',()=>popup.handleColorScopeChoice('rule')));
  assert.equal(s.history.mutationPending(),true,'the step is open while the rule write is prepared');
  release();
  await recolored;
  assert.deepEqual(table(),['b #000000 ','zz #778899 ']);
  assert.equal(s.history.getUndoCount(),1);
  assert.equal(s.history.undoLabel(),'Change color');
}));
test('a Legend recolor of a rule row made while an Add is pending applies after it (OV-444)', () => withDocument(async () => {
  const {s,held,table}=queuedRuleWrites();
  s.state.manualSpecificRules[0].cap='B';
  Object.assign(s.state,{svgContainer:mounted([['CDS','#54bcf8'],['B','#445566']]),selectedResultIndex:{value:0},
    originalLegendOrder:{value:['CDS','B']},deletedLegendEntries:{value:[]},originalLegendColors:{value:{}},
    legendStrokeOverrides:{},adv:{}});
  s.state.legendEntries.value=[{caption:'CDS',originalCaption:'CDS',color:'#54bcf8'},{caption:'B',originalCaption:'B',color:'#445566'}];
  assert.equal(await s.preparation.prepare(s.state.manualSpecificRules),true);
  const legend=createLegendManager({state:s.state,commitLegendRowRules:s.actions.commitSpecificRules,readShownLegendColor:entry=>entry?.color});
  assert.equal(legend.legendRowHasRules(1),true);
  const [add,release]=await held(s.actions.addSpecificRule());
  const recolored=legend.updateLegendEntryColor(1,'#000000');
  release();
  await Promise.all([add,recolored]);
  assert.deepEqual(table(),['b #000000 B','zz #778899 ',added]);
  assert.deepEqual(s.transactions,['Add specific color rule','Change legend color']);
}));
test('a rules TSV import made while an Add is pending applies after it (OV-444)', () => withDocument(async () => {
  const {s,held,table}=queuedRuleWrites();
  const [add,release]=await held(s.actions.addSpecificRule());
  const imported=s.actions.importSpecificRules('CDS\tgene\tb\t#000000\tFile\n');
  release();
  assert.deepEqual(await Promise.all([add,imported]),[true,1]);
  assert.deepEqual(table(),['b #445566 ','zz #778899 ',added,'b #000000 File']);
  assert.deepEqual(s.transactions,['Add specific color rule','Import specific color rules']);
}));
for (const [kind, cancel] of [
  ['a pattern keystroke', s => s.actions.editSpecificRulePattern(s.state.manualSpecificRules[1], 'typing')],
  ['a drawer close', s => s.actions.suspendSpecificRulePatternDrafts()]
]) {
  test(`${kind} cancels only pattern commits, not a pending color edit (OV-444)`, async () => {
    const {s,held,table}=queuedRuleWrites();
    const [edited,release]=await held(s.actions.setSpecificRuleField(0,'color','#000000'));
    cancel(s);
    release();
    assert.equal(await edited,true);
    assert.deepEqual(table(),['b #000000 ','zz #778899 ']);
    assert.deepEqual(s.transactions,['Edit specific color rule']);
  });
}

// The key of the row a rule draws (PD-OI-042 scenario revision 2): Python's answer is
// gbdraw/legend/table.py::rule_legend_keys (tests/test_specific_color_captions.py).
for (const { name, rules, renderer_rows: rendererRows, rule_stroke: ruleStroke = [], keys } of JSON.parse(
  readFileSync(new URL('../fixtures/legend_rule_keys.json', import.meta.url), 'utf8')
).cases) {
  test(`rule row keys: ${name}`, () => {
    // A renderer row is its key and fill, then its line; `rule_stroke` is the line of a rule row.
    const rows = rendererRows.map(([caption, color, stroke, width]) => ({
      caption, color: legendColorIdentity(color), ruleLine: sameLegendLine({ stroke, width }, { stroke: ruleStroke[0], width: ruleStroke[1] })
    }));
    const keyOf = createRuleLegendCaptions(rules, rows);
    assert.deepEqual(rules.map(keyOf), keys);
    assert.deepEqual(rules.map(createRuleLegendCaptions([...rules].reverse(), rows)), keys);
    // The lazy reader of a Result's rows allocates the same keys.
    const pythonRows = new Map(rows.map((row) => [row.caption, { key: row.caption, color: row.color, ruleLine: row.ruleLine }]));
    const lazy = ruleLegendCaptions({ rules, pythonRows, features: [], originalLegendOrder: rows.map((row) => row.caption) });
    assert.deepEqual(rules.map(lazy), keys);
  });
}

// OV-645: the rules of one caption and one fill draw one Legend row, and its key
// changes when the caption starts or stops being written in several colors.
// The Legend owner keeps such a row in its place, so the mirror names it: the
// old key and the new key of each row both rule tables hold.
test('rule row keys: a rule table change re-keys the rows whose caption starts or stops being shared (OV-645)', () => {
  const red = { feat: 'CDS', qual: 'gene', val: 'a', color: '#e63946', cap: 'alpha' };
  const teal = { feat: 'CDS', qual: 'gene', val: 'b', color: '#2a9d8f', cap: 'alpha' };
  const rekeyed = (before, after, rows) => [...rekeyedRuleRows(before, after, rows)];
  assert.deepEqual(rekeyed([red, teal], [red, { ...teal, cap: 'gamma' }]), [['alpha [#e63946]', 'alpha']], 'a name stops being shared');
  assert.deepEqual(rekeyed([red, { ...teal, cap: 'beta' }], [red, teal]), [['alpha', 'alpha [#e63946]']], 'a name becomes shared');
  assert.deepEqual(rekeyed([red, teal], [red]), [['alpha [#e63946]', 'alpha']], 'the other rule of the name is removed');
  assert.deepEqual(rekeyed([red, teal], [teal, { ...red, color: 'rgb(230, 57, 70)' }]), [], 'the same rows in another order and spelling');
  // A recolored or renamed rule draws another row; a row no table change touches keeps its key.
  assert.deepEqual(rekeyed([red, teal], [red, { ...teal, color: '#000000' }]), []);
  assert.deepEqual(rekeyed([red], [{ ...red, cap: 'beta' }]), []);
  const typeRow = [{ caption: 'alpha', color: '#54bcf8', ruleLine: true }];
  assert.deepEqual(rekeyed([red], [red, { ...teal, cap: 'beta' }], typeRow), [], 'a rule named like a renderer row keeps its allocated key');
});

// OV-640 (D-15): a feature's row is read by its key ("Shared [#112233]" for a
// caption written in two colors). A caption written from that reading is the
// row's name, taken through the Legend owner's port; no rule is captioned
// with a key, and "Reset all" resets the rules of the row, not of the name.
const twoColorPopup=()=>{
  const s=setup();
  s.state.errorLog={value:null};
  s.state.manualSpecificRules.push(...rules.map(rule=>({...rule})));
  const [a,b]=s.state.extractedFeatures.value;
  recordRuleMatches([a,b],s.state.manualSpecificRules.map(ruleKey),index=>({matched:[index],priorities:[0],declined:[]}));
  s.setLegendPreparation(async()=>{});
  const popup=popupFor(s,{legendRowName:caption=>(/^Shared \[#[0-9a-f]{6}\]$/.test(caption)?'Shared':caption)});
  const table=()=>s.state.manualSpecificRules.map(rule=>`${rule.qual} ${rule.color} ${rule.cap}`);
  return {s,a,b,popup,table};
};
test('a color for one feature of a keyed row captions its rule with a name, never the row key (OV-640)', () => withDocument(async () => {
  const {s,a,popup,table}=twoColorPopup();
  assert.equal(s.actions.effectiveLegendCaptions()(a),'Shared [#112233]');
  Object.assign(s.state.featureStyleScopeDialog,{show:true,kind:'fill',feat:a,color:'#c83366',legendName:'Shared [#112233]',
    matchingRule:s.state.manualSpecificRules[0],siblingCount:0,ruleMatchCount:2});
  await popup.handleColorScopeChoice('single');
  assert.deepEqual(table().filter(row=>row.includes(' [#')),[],'no caption is a row key');
  assert.deepEqual(table().filter(row=>row.startsWith('hash ')).map(row=>row.split(' ')[1]),['#c83366']);
  assert.equal(table().filter(row=>row.endsWith(' Shared')).length,2,'the feature left its rule\'s row: it is named by itself');
}));
test('the existing color of a keyed row and a hollow fill caption the feature with the row\'s name (OV-640)', () => withDocument(async () => {
  const {s,a,popup,table}=twoColorPopup();
  Object.assign(s.state.featureStyleScopeDialog,{show:true,kind:'fill',feat:a,color:'#c83366',legendName:'Shared [#445566]',
    existingCaptionColor:'#445566',matchingRule:null,siblingCount:0,ruleMatchCount:0});
  await popup.handleColorScopeChoice('useExisting');
  assert.deepEqual(table().filter(row=>row.startsWith('hash ')),['hash #445566 Shared']);
  await popup.setFeatureColorValue(a,'none');
  assert.deepEqual(table().filter(row=>row.startsWith('hash ')),['hash none Shared']);
}));
test('Reset all on a keyed row removes the rules of that row only (OV-640)', () => withDocument(async () => {
  const {s,a,popup,table}=twoColorPopup();
  s.state.clickedFeature.value={feat:a,svg_id:a.svg_id};
  await popup.handleResetColorChoice('all');
  assert.deepEqual(table(),['gene #445566 Shared']);
}));

// The compatibility reading of PD-OI-042 scenario revision 2 (accepted residual
// risk): a caption that ends in its own rule's color ID is read without it from
// a saved Session and from a color table file. Python's answer is
// gbdraw/legend/table.py::saved_rule_caption (tests/test_specific_color_captions.py).
const savedCaptions = JSON.parse(
  readFileSync(new URL('../fixtures/saved_rule_captions.json', import.meta.url), 'utf8')
).cases;
for (const { name, caption, color, read } of savedCaptions) {
  test(`saved rule caption: ${name}`, () => {
    assert.equal(savedRuleCaption(caption, color), read);
  });
}

test('saved rules and an imported color table are read without the color IDs of their own rules', () => {
  const saved = [
    { feat: 'CDS', qual: 'gene', val: '^cox1$', color: '#ffe700', cap: 'Complex IV [#ffe700]', fromFile: true },
    { feat: 'CDS', qual: 'gene', val: '^cox2$', color: '#1f77b4', cap: 'Complex IV [#1f77b4]' },
    { feat: 'CDS', qual: 'gene', val: '^nad1$', color: '#2ca02c', cap: 'Complex I [#ffe700]' }
  ];
  const read = savedSpecificRules(saved);
  assert.deepEqual(read, [
    { ...saved[0], cap: 'Complex IV' }, { ...saved[1], cap: 'Complex IV' }, saved[2]
  ]);
  // The rows are the keys the older build drew as text, so Legend edits saved under them stay on their rows.
  assert.deepEqual(read.map(createRuleLegendCaptions(read)), ['Complex IV [#ffe700]', 'Complex IV [#1f77b4]', 'Complex I [#ffe700]']);
  const text = saved.map((rule) => [rule.feat, rule.qual, rule.val, rule.color, rule.cap].join('\t')).join('\n');
  assert.deepEqual(prepareSpecificColorImport(text).nextRules.map((rule) => rule.cap), ['Complex IV', 'Complex IV', 'Complex I [#ffe700]']);
});

// OV-1100 (R15-19 B, OIC-030): a rule commit writes the drawing it was made
// in, also when the shown drawing changes while its rules are matched, and the
// owner reports the write as pending until it ends, so a switch waits for it.
test('a rule commit writes its own drawing and is pending until it ends', async () => {
  const linear = {
    manualSpecificRules: [], featureColorOverrides: { b: { color: '#abcdef', caption: 'Linear fill' } },
    featureOverrides: {}, addedLegendCaptions: { value: new Set() }, legendEntries: { value: [] },
    legendColorOverrides: {}, newSpecRule: { feat: 'CDS', qual: 'gene', val: 'a', color: '#000000', cap: '' }
  };
  const s = setup(evaluatePythonRules, { linear });
  let release, started;
  const ready = new Promise(resolve => { started = resolve; });
  s.setLegendPreparation(() => new Promise(resolve => { release = resolve; started(); }));
  assert.equal(s.actions.ruleWritePending(), false);
  const pending = s.actions.commitSpecificRules(() => [{ feat: 'CDS', qual: 'gene', val: 'b', color: '#112233', cap: 'B' }]);
  await ready;
  assert.equal(s.actions.ruleWritePending(), true, 'a switch waits for the rule being matched');
  s.state.mode.value = 'linear';
  s.setLegendPreparation(async () => {});
  release();
  await pending;
  assert.equal(s.actions.ruleWritePending(), false);
  assert.deepEqual(s.state.manualSpecificRules.map(rule => rule.cap), ['B'], 'the Circular rules');
  assert.equal(s.state.featureColorOverrides.b?.color, '#112233', 'the Circular fill');
  assert.deepEqual(linear.manualSpecificRules, []);
  assert.deepEqual(linear.featureColorOverrides, { b: { color: '#abcdef', caption: 'Linear fill' } }, 'the Linear fill is kept');
});

// REVIEW-3 R3-B1 (OV-1280): Remove on the Specific rules card while another
// mode is shown finishes that drawing's change at once, in the Remove's own
// step: its rules and its Legend records (rows, added captions, Legend colors,
// re-keyed rows), by the keys of its own stored Legend rows (a caption a hand
// rule shares in another color, a file rule captioned like a type's row).
// Nothing is owed: when the drawing is shown again, its display draws it again
// if its Legend cannot show the rules (app-setup.js).
const hiddenRemoval = (options = {}) => {
  const hand = { feat: 'CDS', qual: 'gene', val: 'a', color: '#445566', cap: 'Hypothetical' };
  const fileRules = [
    { feat: 'CDS', qual: 'gene', val: 'b', color: '#123456', cap: 'Hypothetical', fromFile: true },
    { feat: 'CDS', qual: 'gene', val: 'b', color: '#654321', cap: 'CDS', fromFile: true },
    { feat: 'CDS', qual: 'gene', val: 'b', color: '#abcdef', cap: 'Imported', fromFile: true }
  ];
  const drawn = [['CDS', '#54bcf8'], ['Hypothetical [#445566]', '#445566'], ['Hypothetical [#123456]', '#123456'],
    ['CDS [#654321]', '#654321'], ['Imported', '#abcdef']];
  const s = setup(undefined, options);
  s.state.manualSpecificRules.splice(0, s.state.manualSpecificRules.length, { ...hand }, ...fileRules.map((rule) => ({ ...rule })));
  // The user renamed the file rule's row before the Remove.
  s.state.legendEntries.value = drawn.map(([key, color]) => ({ caption: key === 'Imported' ? 'Imported protein' : key, originalCaption: key, color }));
  s.state.addedLegendCaptions.value = new Set(drawn.slice(1).map(([key]) => key));
  Object.assign(s.state.legendColorOverrides, { 'Imported': '#abcdef', 'CDS': '#000000' });
  s.state.mode = { value: 'linear' };
  withModeDrawings(s.state, { circular: s.state, linear: { manualSpecificRules: [] } });
  const legend = () => ({
    entries: s.state.legendEntries.value.map((entry) => [entry.originalCaption, entry.caption]),
    added: [...s.state.addedLegendCaptions.value],
    colors: { ...s.state.legendColorOverrides }
  });
  return { s, hand, drawn, legend };
};
const removedLegend = {
  entries: [['CDS', 'CDS'], ['Hypothetical', 'Hypothetical']],
  added: ['Hypothetical'],
  colors: { 'CDS': '#000000' }
};

test('a Specific rules removal in a drawing not shown takes back its Legend rows at once, so a Save before it is shown lists none', async () => {
  const { s, hand, legend } = hiddenRemoval();
  assert.equal(await s.actions.removeImportedSpecificRules(s.state), true);
  // What Save Session writes for the drawing now.
  assert.deepEqual(s.state.manualSpecificRules, [hand]);
  assert.deepEqual(legend(), removedLegend);
  assert.equal(s.rerenders(), 0, 'nothing is drawn while it is not shown');

  // The shown path allocates the same keys: its re-keyed row and its rows.
  const shown = setup();
  shown.state.manualSpecificRules.splice(0, 0, { ...hand }, { feat: 'CDS', qual: 'gene', val: 'b', color: '#123456', cap: 'Hypothetical', fromFile: true });
  shown.state.originalLegendOrder = { value: ['CDS', 'Hypothetical [#445566]', 'Hypothetical [#123456]'] };
  shown.state.svgContainer = mounted([['CDS', '#54bcf8'], ['Hypothetical [#445566]', '#445566'], ['Hypothetical [#123456]', '#123456']]);
  let intents = null;
  shown.setLegendPreparation(async (next) => { intents = next; });
  assert.equal(await shown.actions.removeImportedSpecificRules(shown.state), true);
  assert.deepEqual([...shown.rekeyed()], [['Hypothetical [#445566]', 'Hypothetical']]);
  assert.deepEqual(intents.map((intent) => intent.caption), removedLegend.added);
});

// R3-B1 route 4: a row Python no longer draws leaves no rename behind, so the
// next Generate of that drawing does not fail on it.
test('a renamed Legend row of a removed file rule leaves no rename for the next Generate of a drawing not shown', async () => {
  const { s, drawn } = hiddenRemoval();
  assert.equal(await s.actions.removeImportedSpecificRules(s.state), true);
  const plan = compileDirectEditorMutationPlan({
    catalogAdmission: { resultNames: ['figure'], renderedTargetsByOverrideKey: new Map(), resultIndexesByRenderedId: new Map() },
    legendEntries: s.state.legendEntries.value,
    originalLegendOrder: drawn.map(([key]) => key)
  });
  assert.deepEqual(plan.operationsByResult[0].legendRenames.filter(({ from }) => from === 'Imported'), []);
});

// R3-B1 route 2: the Remove's step holds the whole change, so Undo then Redo,
// also after the drawing was shown, gives the Remove again.
test('Undo then Redo of a Specific rules removal in a drawing not shown gives the removal again', () => withDocument(async () => {
  // The checkpoint holds what History's snapshot of the other mode holds.
  const { s, legend } = hiddenRemoval({ createHistory: (state) => createHistoryManager({
    buildIntent: () => ({}), applyIntent: () => {},
    buildCheckpoint: () => structuredClone({ rules: state.manualSpecificRules, entries: state.legendEntries.value,
      added: [...state.addedLegendCaptions.value], colors: { ...state.legendColorOverrides } }),
    applyCheckpoint: ({ rules, entries, added, colors }) => {
      state.manualSpecificRules.splice(0, Infinity, ...structuredClone(rules));
      state.legendEntries.value = structuredClone(entries);
      state.addedLegendCaptions.value = new Set(added);
      Object.keys(state.legendColorOverrides).forEach((key) => { delete state.legendColorOverrides[key]; });
      Object.assign(state.legendColorOverrides, colors);
    } }) });
  await s.history.initializeIntentBaseline();
  const before = legend();
  // components.js `clearFile`: the card's Remove is one checkpoint step.
  await s.history.runUndoableCheckpoint('Remove Specific rules file', () => s.actions.removeImportedSpecificRules(s.state));
  assert.deepEqual(legend(), removedLegend);
  s.state.mode.value = 'circular';
  await s.history.undo();
  assert.deepEqual(legend(), before);
  await s.history.redo();
  assert.deepEqual(legend(), removedLegend);
  assert.equal(s.state.manualSpecificRules.some((rule) => rule.fromFile), false);
}));

// R3-B1 route 3: nothing waits for the drawing to be shown, so a switch to it
// and back before anything was prepared leaves its Legend as the Remove wrote it.
test('a mode switch to a drawing whose file rules were removed while it was not shown, and back, changes nothing', async () => {
  const { s, legend } = hiddenRemoval();
  const removal = s.actions.removeImportedSpecificRules(s.state);
  s.state.mode.value = 'circular';
  s.state.mode.value = 'linear';
  assert.equal(await removal, true);
  assert.deepEqual(legend(), removedLegend);
});

// REVIEW-4 B1: the hidden take-back finds a removed rule group's Legend row by
// the key the drawing stores for it (its allocated key when one is listed,
// else its caption), never by keys allocated again from the colors the rows
// list, and a deleted row is one Python draws. So a type row that shares the
// caption keeps its rename, color, and Delete.
const hiddenRuleRow = ({ rule, rules = [rule], entries, deleted = [], colors = {}, options = {} }) => {
  const s = setup(undefined, options);
  s.state.manualSpecificRules.splice(0, Infinity, ...rules.map((one) => ({ ...one, fromFile: true })));
  s.state.legendEntries.value = entries.map(([key, caption, color]) => ({ caption, originalCaption: key, color }));
  s.state.deletedLegendEntries.value = deleted.map(([key, caption, color]) => ({ caption, originalCaption: key, color }));
  s.state.addedLegendCaptions.value = new Set(entries.map(([key]) => key).filter((key) => key.includes(' [')));
  Object.assign(s.state.legendColorOverrides, colors);
  s.state.mode = { value: 'linear' };
  withModeDrawings(s.state, { circular: s.state, linear: { manualSpecificRules: [] } });
  const plannedRenames = (order) => compileDirectEditorMutationPlan({
    catalogAdmission: { resultNames: ['figure'], renderedTargetsByOverrideKey: new Map(), resultIndexesByRenderedId: new Map() },
    legendEntries: s.state.legendEntries.value, originalLegendOrder: order
  }).operationsByResult[0].legendRenames;
  return { s, plannedRenames };
};
const listedRows = (entries) => entries.map((entry) => [entry.originalCaption, entry.caption, entry.color]);

test('a hidden Remove of a file rule drawn beside a deleted type row of its caption keeps the Delete (REVIEW-4 B1)', async () => {
  const { s, plannedRenames } = hiddenRuleRow({
    rule: { feat: 'CDS', qual: 'gene', val: 'b', color: '#654321', cap: 'CDS' },
    entries: [['CDS [#654321]', 'Imported', '#654321']],
    deleted: [['CDS', 'CDS', '#54bcf8']]
  });
  assert.equal(await s.actions.removeImportedSpecificRules(s.state), true);
  assert.deepEqual(listedRows(s.state.legendEntries.value), [], 'the file rule row leaves');
  assert.deepEqual(listedRows(s.state.deletedLegendEntries.value), [['CDS', 'CDS', '#54bcf8']], 'the type row stays deleted');
  assert.deepEqual([...s.state.addedLegendCaptions.value], []);
  assert.deepEqual(plannedRenames(['CDS', 'CDS [#654321]']), []);
});

for (const [variant, colors] of [['recolored in the Legend', { rRNA: '#ff0000' }], ['drawn in that fill with its own line (D-13)', {}]]) {
  test(`a hidden Remove of a file rule keeps the rename and color of a type row of its caption ${variant} (REVIEW-4 B1)`, async () => {
    const { s, plannedRenames } = hiddenRuleRow({
      rule: { feat: 'tRNA', qual: 'product', val: 'x', color: '#ff0000', cap: 'rRNA' },
      entries: [['rRNA', 'Ribosomal RNA', '#ff0000'], ['rRNA [#ff0000]', 'Imported', '#ff0000']],
      colors
    });
    assert.equal(await s.actions.removeImportedSpecificRules(s.state), true);
    assert.deepEqual(listedRows(s.state.legendEntries.value), [['rRNA', 'Ribosomal RNA', '#ff0000']]);
    assert.deepEqual({ ...s.state.legendColorOverrides }, colors, 'the type row keeps its color');
    assert.deepEqual([...s.state.addedLegendCaptions.value], []);
    assert.deepEqual(plannedRenames(['rRNA', 'rRNA [#ff0000]']), [{ from: 'rRNA', to: 'Ribosomal RNA', allowMissing: false }],
      'no rename is left for the row Python no longer draws');
  });
}

// REVIEW-5 B1: the Legend color copy of a removed file rule row (OV-61) is
// retired by the Legend owner's writer, not by the rule owner (AUTH-ORIGIN-1).
test('a hidden Remove retires the Legend color copy of a removed file rule row through the Legend owner (REVIEW-5 B1)', async () => {
  /** @type {Record<string, string>[]} */
  const seen = [];
  let writer = null;
  const { s } = hiddenRuleRow({
    rule: { feat: 'CDS', qual: 'gene', val: 'b', color: '#ff0000', cap: 'File' },
    entries: [['File', 'File', '#ff0000']],
    colors: { File: '#ff0000' },
    options: { moveLegendRecords: (...args) => { seen.push({ ...args[0].legendColorOverrides }); return writer(...args); } }
  });
  writer = createLegendManager({ state: s.state, commitLegendRowRules: () => true }).moveLegendRecords;
  assert.equal(await s.actions.removeImportedSpecificRules(s.state), true);
  assert.deepEqual(seen, [{ File: '#ff0000' }], 'the copy is still there when the Legend owner is called');
  assert.deepEqual({ ...s.state.legendColorOverrides }, {});
  assert.deepEqual(listedRows(s.state.legendEntries.value), []);
});

// Two removed rule groups stored under one key (no allocated key is listed)
// each retire their own Legend color copy (OV-152): the copy of either color goes.
for (const copy of ['#ff0000', '#0000ff']) {
  test(`a hidden Remove of two file rules stored under one key retires the Legend color copy of either (${copy})`, async () => {
    const { s } = hiddenRuleRow({
      rules: [{ feat: 'CDS', qual: 'gene', val: 'b', color: '#ff0000', cap: 'File' },
        { feat: 'CDS', qual: 'gene', val: 'c', color: '#0000ff', cap: 'File' }],
      entries: [],
      colors: { File: copy, Other: copy }
    });
    assert.equal(await s.actions.removeImportedSpecificRules(s.state), true);
    assert.deepEqual({ ...s.state.legendColorOverrides }, { Other: copy }, 'the copy goes; a color of another row stays');
  });
}

// D-21 (Q-ORIGIN-1 = A): preset rows are the user's, as hand-made rows are:
// Remove and a new rules file never take them, in any drawing.
const withPreset = (s, text) => {
  s.state.specificRulePresets.splice(0, Infinity, { id: 'preset', path: 'preset.tsv' });
  s.state.selectedSpecificPreset.value = 'preset';
  const fetched = globalThis.fetch;
  globalThis.fetch = async () => ({ ok: true, text: async () => text });
  return s.actions.applySpecificRulePreset().finally(() => { globalThis.fetch = fetched; });
};
const presetText = 'CDS\tgene\ta\t#0000ff\tPreset\n';
const fileText = 'CDS\tgene\tb\t#ff0000\tFile\n';
const rows = (s) => s.state.manualSpecificRules.map((rule) => [rule.cap, rule.fromFile === true]);

test('a preset applied after a rules file keeps its rows when the card is removed (D-21)', async () => {
  const s = setup();
  assert.equal(await s.actions.importSpecificRules(fileText), 1);
  assert.deepEqual(rows(s), [['File', true]]);
  await withPreset(s, presetText);
  assert.deepEqual(rows(s), [['Preset', false]]);
  assert.equal(await s.actions.removeImportedSpecificRules(s.state), true);
  assert.deepEqual(rows(s), [['Preset', false]]);
});

test('a rules file after a preset adds its rows after the preset rows, and Remove gives the preset rows back (D-21)', async () => {
  const s = setup();
  await withPreset(s, presetText);
  assert.equal(await s.actions.importSpecificRules(fileText), 1);
  assert.deepEqual(rows(s), [['Preset', false], ['File', true]]);
  assert.equal(await s.actions.removeImportedSpecificRules(s.state), true);
  assert.deepEqual(rows(s), [['Preset', false]]);
});

// watchers.js: a rules file loaded in one mode takes back the card's earlier
// rows in the other drawing; a preset there has none.
test('a rules file loaded in one mode leaves a preset of the drawing not shown as it is (D-21)', async () => {
  const s = setup();
  await withPreset(s, presetText);
  s.state.mode = { value: 'linear' };
  withModeDrawings(s.state, { circular: s.state, linear: { manualSpecificRules: [] } });
  const listed = JSON.stringify([s.state.manualSpecificRules, s.state.legendEntries.value]);
  assert.equal(await s.actions.removeImportedSpecificRules(s.state), true);
  assert.equal(JSON.stringify([s.state.manualSpecificRules, s.state.legendEntries.value]), listed);
});
