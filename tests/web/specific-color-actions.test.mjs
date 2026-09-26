import assert from 'node:assert/strict';
import { test } from 'node:test';
import { createFeatureRuleActions } from '../../gbdraw/web/js/app/feature-editor/rule-actions.js';
import { createRulePreparation } from '../../gbdraw/web/js/app/rule-matching.js';
import { evaluatePythonRules } from './helpers/python-rule-evaluator.mjs';

const setup = () => {
  const state = {
    manualSpecificRules: [], extractedFeatures: { value: [
      { type: 'CDS', svg_id: 'a', qualifiers: { gene: ['a'] } },
      { type: 'CDS', svg_id: 'b', qualifiers: { gene: ['b'] } }
    ] },
    featureColorOverrides: {}, results: { value: [{name:'figure',content:'before'}] },
    svgResultIdentity: { value:'before' }, fileLegendCaptions: { value:new Set() }, addedLegendCaptions: { value:new Set() },
    legendEntries: { value:[] }, files: {t_color:null},
    newSpecRule: {feat:'CDS',qual:'gene',val:'a',color:'#112233',cap:'Shared'}
  };
  const notices = [], transactions = [];
  const preparation = createRulePreparation({state, evaluate:evaluatePythonRules, notify:message=>notices.push(message)});
  let prepareLegend = async () => {};
  let previousIntents = [];
  const actions = createFeatureRuleActions({state, rulePreparation:preparation, history:{
    runUndoable: async (label, commit) => {
      const before=JSON.stringify(state.manualSpecificRules);
      await commit();
      if (before!==JSON.stringify(state.manualSpecificRules)) transactions.push(label);
    }
  }, legendActions: {
    syncFileLegendEntries: async (intents, {isCurrent,commit,previousFileIntents}) => {
      previousIntents=previousFileIntents;
      await prepareLegend(intents);
      if (!isCurrent()) return false;
      commit();
      state.legendEntries.value=intents;
      state.results.value=[{name:'figure',content:'after'}];
      return true;
    }
  }, svgActions:{applyPaletteToSvg(){},applySpecificRulesToSvg(){}}, nextTick:async()=>{}});
  return {state,actions,preparation,notices,transactions, setLegendPreparation: fn => {prepareLegend=fn;}, previousIntents:()=>previousIntents};
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
