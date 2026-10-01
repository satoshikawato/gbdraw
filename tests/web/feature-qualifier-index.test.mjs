import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { test } from 'node:test';
import vm from 'node:vm';
import { createFeatureRuleActions } from '../../gbdraw/web/js/app/feature-editor/rule-actions.js';
import { assertKnownDefect } from './helpers/known-defect.mjs';

const context = vm.createContext({ console });
vm.runInContext(readFileSync('gbdraw/web/vendor/vue/vue.global.js', 'utf8'), context);
const { ref, computed } = context.Vue;

const feature = (type, id) => ({ type, svg_id: id });
const create = (features) => {
  const extractedFeatures = ref(features);
  const state = {
    extractedFeatures,
    manualSpecificRules: [],
    featureColorOverrides: {},
    fileLegendCaptions: ref(new Set()),
    addedLegendCaptions: ref(new Set()),
    legendEntries: ref([]),
    results: ref([]),
    svgResultIdentity: ref(''),
    newSpecRule: {},
    newPriorityRule: {}
  };
  const actions = createFeatureRuleActions({
    state, ref, computed, nextTick: async () => {},
    rulePreparation: { isCurrent: () => true },
    history: {}, legendActions: {}, svgActions: {}
  });
  return { actions, extractedFeatures };
};

test('feature hash counts keep type, duplicate rendered identity and replacement semantics', () => {
  const a = feature('CDS', 'same_record_1');
  const b = feature('CDS', 'same_record_2');
  const otherType = feature('tRNA', 'same_record_3');
  const s = create([a, b, otherType]);
  assert.deepEqual(s.actions.getFeatureQualifier(a), { qual: 'hash', val: 'same_record_1' });
  assert.deepEqual(s.actions.getFeatureQualifier(otherType), { qual: 'hash', val: 'same' });
  s.extractedFeatures.value = [a, otherType];
  assert.deepEqual(s.actions.getFeatureQualifier(a), { qual: 'hash', val: 'same' });
  s.extractedFeatures.value.push(feature('CDS', 'same_record_4'));
  assert.deepEqual(s.actions.getFeatureQualifier(a), { qual: 'hash', val: 'same_record_1' });
});

test('qualifier preparation does not scan the full feature list for each target', () => {
  const features = Array.from({ length: 3000 }, (_, index) => feature('CDS', `shared_record_${index + 1}`));
  const s = create(features);
  assert.deepEqual(s.actions.getFeatureQualifier(features[0]), { qual: 'hash', val: 'shared_record_1' });
  assert.deepEqual(s.actions.getFeatureQualifier(features.at(-1)), { qual: 'hash', val: 'shared_record_3000' });
});

// FE-09 (D-14, PD-OI-069): Python matches a color rule only by the stable hash
// (gbdraw/features/selector_values.py), so "This feature only" on a record whose
// ID repeats must also write the stable hash. The first test above records the
// current rendered-instance fallback; the FE-09 fix rewrites it.
test('duplicate record IDs still write the stable hash Python matches (FE-09 known defect)', async () => {
  const first = feature('CDS', 'same_record_1');
  const second = feature('CDS', 'same_record_2');
  const s = create([first, second]);
  await assertKnownDefect('FE-09', () => {
    assert.deepEqual(s.actions.getFeatureQualifier(first), { qual: 'hash', val: 'same' });
    assert.deepEqual(s.actions.getFeatureQualifier(second), { qual: 'hash', val: 'same' });
  });
});
