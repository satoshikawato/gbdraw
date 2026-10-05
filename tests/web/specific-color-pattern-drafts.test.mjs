import assert from 'node:assert/strict';
import { test } from 'node:test';
import { execFileSync } from 'node:child_process';
import { createFeatureRuleActions } from '../../gbdraw/web/js/app/feature-editor/rule-actions.js';
import { DiagramGenerationCanceledError } from '../../gbdraw/web/js/services/diagram-generation.js';
import { createRulePreparation, ruleMatchesFeature } from '../../gbdraw/web/js/app/rule-matching.js';

const nativeEvaluate = async payload => {
  const response = JSON.parse(execFileSync(process.env.PYTHON || 'python', ['-c', `
import json,sys
from gbdraw.web_support.rule_matching import evaluate_rules_json
from gbdraw.web_support.error_adapter import serialize_web_error,private_web_execution
p=json.load(sys.stdin)
try:
    with private_web_execution():
        result=json.loads(evaluate_rules_json(json.dumps(p['features']),json.dumps(p['rules']),p['kind']))
    print(json.dumps({'result':result}))
except Exception as e:
    print(json.dumps({'error':serialize_web_error(e,operation='evaluateRules',stage='rule-validation')}))
`], { input: JSON.stringify(payload), encoding: 'utf8', maxBuffer: 32 * 1024 * 1024 }));
  if (response.error) throw response.error;
  return response.result;
};
const rule = (val = 'NADH') => ({ feat: 'CDS', qual: 'product', val, color: '#f01234', cap: '' });
const setup = ({ features = [{ type: 'CDS', svg_id: 'one', qualifiers: { product: ['NADH'] } }], evaluate = nativeEvaluate } = {}) => {
  const state = {
    manualSpecificRules: [rule()], extractedFeatures: { value: features },
    featureColorOverrides: {}, errorLog: { value: null },
    results: { value: [{ name: 'figure', content: '<svg>accepted</svg>' }] },
    svgResultIdentity: { value: 'accepted' }, mode: { value: 'circular' }, generatedMode: { value: 'circular' },
    fileLegendCaptions: { value: new Set() }, addedLegendCaptions: { value: new Set() },
    legendEntries: { value: [] }, files: { t_color: null }
  };
  const calls = [], history = [];
  let available = true;
  let finalize = async () => {};
  const preparation = createRulePreparation({ state, evaluate: payload => { calls.push(payload.kind); return evaluate(payload); } });
  const transact = async (label, commit) => {
    const before = JSON.stringify(state.manualSpecificRules);
    await commit();
    if (before !== JSON.stringify(state.manualSpecificRules)) history.push(label);
    await finalize();
  };
  const actions = createFeatureRuleActions({ state, ref: value => ({ value }), computed: get => ({ get value() { return get(); } }),
    rulePreparation: preparation, isPatternEditAvailable: () => available, nextTick: async () => {},
    history: { runUndoableCheckpoint: transact, runUndoable: transact },
    prepareFileLegendEntries: async (_, { isCurrent }) => isCurrent() && {
      diff: { add: [], remove: [] },
      isCurrent: () => true,
      apply: () => { state.results.value = [{ name: 'figure', content: JSON.stringify(state.manualSpecificRules) }]; }
    }, projectPaletteAndRules: () => true
  });
  const row = state.manualSpecificRules[0];
  const stable = () => ({ canonical: JSON.stringify(state.manualSpecificRules), result: state.results.value, history: history.length });
  return { state, row, actions, preparation, calls, history, stable, setAvailable: value => { available = value; }, setFinalize: fn => { finalize = fn; } };
};

for (const features of [[], [{ type: 'tRNA', svg_id: 'other', qualifiers: { product: ['unrelated'] } }]]) {
  test('native invalid syntax retains display text and classified field cause even without matching features', async () => {
    const s = setup({ features }); const before = s.stable();
    assert.equal(await s.actions.setSpecificRuleField(0, 'val', '😀PRIVATE_PATTERN_SENTINEL['), false);
    assert.equal(s.actions.specificRulePattern(s.row), '😀PRIVATE_PATTERN_SENTINEL[');
    const draft = s.actions.specificRulePatternDraft(s.row);
    assert.equal(draft.pending, false);
    assert.equal(draft.error.code, 'REGEX_SYNTAX');
    assert.equal(draft.error.context.positionUnit, 'python-character');
    assert.equal(draft.error.context.position, [...'😀PRIVATE_PATTERN_SENTINEL'].length);
    assert(!JSON.stringify(draft.error).includes('PRIVATE_PATTERN_SENTINEL'));
    assert.deepEqual(s.stable(), before);
    assert.equal(await s.actions.retrySpecificRulePattern(s.row), false);
    assert.deepEqual(s.stable(), before);
    s.actions.revertSpecificRulePattern(s.row);
    assert.equal(s.actions.specificRulePattern(s.row), 'NADH');
    assert.equal(s.actions.specificRulePatternDraft(s.row), null);
    assert.deepEqual(s.stable(), before);
    assert.deepEqual(s.calls, ['color-captions', 'color', 'color-captions', 'color']);
  });
}

for (const failure of [
  { code: 'WORKER_INIT', stage: 'initialization' },
  { code: 'RESOURCE_INVALID', stage: 'resource-staging', context: { reason: 'WORKSPACE' } }
]) {
  test(`${failure.stage} failure remains a correctable draft and Retry commits exactly once`, async () => {
    let fail = true;
    const s = setup({ evaluate: payload => fail ? Promise.reject({ ...failure, message: 'PRIVATE_EXCEPTION_SENTINEL', pattern: 'PRIVATE_PATTERN_SENTINEL' }) : nativeEvaluate(payload) });
    const before = s.stable();
    assert.equal(await s.actions.setSpecificRuleField(0, 'val', '(?P<enzyme>NADH)'), false);
    assert.deepEqual(s.stable(), before);
    const error = s.actions.specificRulePatternDraft(s.row).error;
    assert.equal(error.code, failure.code);
    assert.equal(error.stage, failure.stage);
    assert(!JSON.stringify(error).includes('PRIVATE_'));
    assert.equal(await s.actions.retrySpecificRulePattern(s.row), false);
    assert.deepEqual(s.stable(), before);
    fail = false;
    assert.equal(await s.actions.retrySpecificRulePattern(s.row), true);
    assert.equal(s.state.manualSpecificRules[0].val, '(?P<enzyme>NADH)');
    assert.equal(s.actions.specificRulePatternDraft(s.state.manualSpecificRules[0]), null);
    assert.notEqual(s.stable().result, before.result);
    assert.equal(s.history.length, 1);
    assert.deepEqual(s.calls, ['color-captions', 'color-captions', 'color-captions', 'color']);
  });
}

test('keystrokes never evaluate; failed correction preserves all accepted owners and valid Python correction commits once', async () => {
  const s = setup(); const before = s.stable();
  s.actions.editSpecificRulePattern(s.row, '[');
  assert.equal(s.calls.length, 0);
  assert.equal(s.actions.specificRulePatternDraft(s.row).pending, false);
  await s.actions.setSpecificRuleField(0, 'val', '[');
  await s.actions.setSpecificRuleField(0, 'val', '(?<js>NADH)');
  assert.equal(s.actions.specificRulePattern(s.row), '(?<js>NADH)');
  assert.deepEqual(s.stable(), before);
  assert.equal(await s.actions.setSpecificRuleField(0, 'val', '(?i)NADH'), true);
  assert.equal(s.state.manualSpecificRules[0].val, '(?i)NADH');
  assert.equal(s.history.length, 1);
});

const hold = () => {
  let release, entered;
  const ready = new Promise(resolve => { entered = resolve; });
  let held = true;
  const evaluate = async payload => {
    if (held) {
      held = false;
      await new Promise(resolve => { release = resolve; entered(); });
    }
    return nativeEvaluate(payload);
  };
  return { evaluate, ready, release: () => release() };
};
for (const change of ['new edit', 'new keystroke', 'remove', 'reorder', 'drawer', 'mode cycle', 'catalog', 'Result', 'document', 'Session pending']) {
  test(`held old edit cannot commit after ${change}`, async () => {
    const gate = hold(); const s = setup({ evaluate: gate.evaluate });
    s.state.manualSpecificRules.push(rule('other'));
    const before = s.stable();
    const pending = s.actions.setSpecificRuleField(0, 'val', '(?i)NADH');
    await gate.ready;
    assert.equal(s.actions.specificRulePatternDraft(s.row).pending, true);
    assert.deepEqual(s.stable(), before);
    if (change === 'new edit') await s.actions.setSpecificRuleField(0, 'val', '(?P<enzyme>NADH)');
    if (change === 'new keystroke') s.actions.editSpecificRulePattern(s.row, '[new');
    if (change === 'remove') await s.actions.removeSpecificRule(0);
    if (change === 'reorder') await s.actions.moveSpecificRuleDown(0);
    if (change === 'drawer') s.actions.suspendSpecificRulePatternDrafts();
    if (change === 'mode cycle') { s.actions.suspendSpecificRulePatternDrafts(); s.state.mode.value = 'linear'; s.state.mode.value = 'circular'; }
    if (change === 'catalog') s.state.extractedFeatures.value = [];
    if (change === 'Result') s.state.svgResultIdentity.value = 'replaced';
    if (change === 'document') s.actions.clearSpecificRulePatternDrafts();
    if (change === 'Session pending') { s.actions.suspendSpecificRulePatternDrafts(); s.setAvailable(false); }
    const current = s.stable();
    gate.release();
    assert.equal(await pending, false);
    assert.deepEqual(s.stable(), current);
    if (['remove', 'document'].includes(change)) assert.equal(s.actions.specificRulePatternDraft(s.row), null);
    if (change === 'reorder') {
      assert.equal(s.state.manualSpecificRules[1], s.row);
      assert.equal(s.actions.specificRulePattern(s.row), '(?i)NADH');
    }
    if (change === 'new edit') assert.equal(s.state.manualSpecificRules[0].val, '(?P<enzyme>NADH)');
    if (change === 'new keystroke') assert.equal(s.actions.specificRulePattern(s.row), '[new');
    if (change === 'drawer') assert.equal(s.actions.specificRulePatternDraft(s.row).pending, false);
  });
}

test('invisible mode cannot apply; unchanged bulk History and failed Session restore keep draft; target replacement and reset release it', async () => {
  const s = setup(); await s.actions.setSpecificRuleField(0, 'val', '[');
  const before = s.stable();
  s.actions.suspendSpecificRulePatternDrafts(); s.state.mode.value = 'linear';
  assert.equal(await s.actions.retrySpecificRulePattern(s.row), false);
  assert.deepEqual(s.stable(), before);
  s.state.mode.value = 'circular';
  for (const operation of ['unrelated Undo', 'unrelated Redo', 'failed Session rollback']) {
    const saved = s.actions.captureSpecificRulePatternDrafts();
    const id = s.actions.specificRulePatternFieldId(s.state.manualSpecificRules[0]);
    s.state.manualSpecificRules.splice(0, 1, { ...s.state.manualSpecificRules[0] });
    s.actions.restoreSpecificRulePatternDrafts(saved);
    const row = s.state.manualSpecificRules[0];
    assert.equal(s.actions.specificRulePattern(row), '[', operation);
    assert.equal(s.actions.specificRulePatternDraft(row).error.code, 'REGEX_SYNTAX');
    assert.equal(s.actions.specificRulePatternFieldId(row), id);
  }
  const saved = s.actions.captureSpecificRulePatternDrafts();
  s.state.manualSpecificRules[0] = rule('other');
  s.actions.restoreSpecificRulePatternDrafts(saved);
  assert.equal(s.actions.specificRulePatternDraft(s.state.manualSpecificRules[0]), null);
  s.actions.restoreSpecificRulePatternDrafts(saved); // Redo never revives a discarded draft.
  assert.equal(s.actions.specificRulePatternDraft(s.state.manualSpecificRules[0]), null);
  s.actions.editSpecificRulePattern(s.state.manualSpecificRules[0], '[');
  const oldDocument = s.actions.captureSpecificRulePatternDrafts();
  s.actions.clearSpecificRulePatternDrafts();
  s.actions.restoreSpecificRulePatternDrafts(oldDocument);
  assert.equal(s.actions.specificRulePatternDraft(s.state.manualSpecificRules[0]), null);
});

test('duplicate accepted patterns remain bound to the actual row through reorder and removal', async () => {
  const s = setup(); s.state.manualSpecificRules.push(rule());
  const second = s.state.manualSpecificRules[1];
  s.actions.editSpecificRulePattern(s.row, '[first');
  s.actions.editSpecificRulePattern(second, '[second');
  await s.actions.moveSpecificRuleDown(0);
  assert.equal(s.actions.specificRulePattern(s.state.manualSpecificRules[0]), '[second');
  assert.equal(s.actions.specificRulePattern(s.state.manualSpecificRules[1]), '[first');
  await s.actions.removeSpecificRule(0);
  assert.equal(s.actions.specificRulePatternDraft(second), null);
  assert.equal(s.actions.specificRulePattern(s.state.manualSpecificRules[0]), '[first');
});

test('25,000 features use the existing single matching evaluation, cached retry preparation and no keystroke helper', async () => {
  const features = Array.from({ length: 25000 }, (_, i) => ({ type: 'CDS', svg_id: `f${i}`, qualifiers: { product: [i % 2 ? 'other' : 'NADH'] } }));
  const s = setup({ features });
  const start = performance.now();
  await s.actions.setSpecificRuleField(0, 'val', '(?P<enzyme>NADH)');
  assert.equal(s.history.length, 1);
  assert.deepEqual(s.calls, ['color-captions', 'color']);
  assert.equal(features.filter(f => ruleMatchesFeature(f, s.state.manualSpecificRules[0])).length, 12500);
  assert.equal(s.preparation.prepare(), true);
  const row = s.state.manualSpecificRules[0];
  s.actions.editSpecificRulePattern(row, '[never evaluated');
  s.actions.revertSpecificRulePattern(row);
  assert.deepEqual(s.calls, ['color-captions', 'color']);
  console.log(JSON.stringify({ featureCount: features.length, preparationMs: performance.now() - start, matchingCalls: 1, captionCalls: 1, history: s.history.length }));
});


test('a committed edit finishing History cannot erase a later field draft', async () => {
  const s = setup();
  let release, started;
  const ready = new Promise(resolve => { started = resolve; });
  s.setFinalize(() => new Promise(resolve => { release = resolve; started(); }));
  const applied = s.actions.setSpecificRuleField(0, 'val', '(?P<enzyme>NADH)');
  await ready;
  assert.equal(s.state.manualSpecificRules[0], s.row);
  assert.equal(s.row.val, '(?P<enzyme>NADH)');
  s.actions.editSpecificRulePattern(s.row, '[new correction');
  release();
  assert.equal(await applied, true);
  assert.equal(s.actions.specificRulePattern(s.row), '[new correction');
  assert.equal(s.actions.specificRulePatternDraft(s.row).pending, false);
  assert.equal(s.history.length, 1);
});


test('transport cancellation keeps the field draft without publishing a failure or changing accepted owners', async () => {
  const s = setup({ evaluate: () => Promise.reject(new DiagramGenerationCanceledError()) });
  const before = s.stable();
  assert.equal(await s.actions.setSpecificRuleField(0, 'val', '(?P<enzyme>NADH)'), false);
  assert.equal(s.actions.specificRulePattern(s.row), '(?P<enzyme>NADH)');
  assert.equal(s.actions.specificRulePatternDraft(s.row).pending, false);
  assert.equal(s.actions.specificRulePatternDraft(s.row).error, null);
  assert.equal(s.state.errorLog.value, null);
  assert.deepEqual(s.stable(), before);
});

for (const operation of ['save', 'load']) {
  test(`Session ${operation} leaves corrective pattern drafts and accepted owners unchanged`, async () => {
    const s = setup();
    s.actions.editSpecificRulePattern(s.row, '[');
    const draft = s.actions.specificRulePatternDraft(s.row);
    const before = s.stable();
    const busy = { status: 'busy', operation };
    s.state.sessionOperationAvailability = () => busy;
    assert.deepEqual(s.actions.editSpecificRulePattern(s.row, 'changed'), busy);
    assert.deepEqual(s.actions.revertSpecificRulePattern(s.row), busy);
    assert.deepEqual(await s.actions.retrySpecificRulePattern(s.row), busy);
    assert.deepEqual(await s.actions.setSpecificRuleField(0, 'val', 'changed'), busy);
    assert.deepEqual(s.actions.retrySpecificRuleFailure(), busy);
    assert.equal(s.actions.specificRulePatternDraft(s.row), draft);
    assert.equal(draft.text, '[');
    assert.equal(draft.pending, false);
    assert.equal(s.calls.length, 0);
    assert.deepEqual(s.stable(), before);
    s.state.sessionOperationAvailability = () => null;
    s.actions.revertSpecificRulePattern(s.row);
    assert.equal(s.actions.specificRulePatternDraft(s.row), null);
  });
}
