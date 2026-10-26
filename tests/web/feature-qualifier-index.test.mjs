import assert from 'node:assert/strict';
import { test } from 'node:test';
import { createFeatureRuleActions } from '../../gbdraw/web/js/app/feature-editor/rule-actions.js';
import { loadPageVue } from './helpers/vue-runtime.mjs';

const { ref, computed } = loadPageVue();

const feature = (type, id) => ({ type, svg_id: id });
const create = (features) => {
  const extractedFeatures = ref(features);
  const state = {
    extractedFeatures,
    manualSpecificRules: [],
    featureColorOverrides: {},
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
    projectPaletteAndRules: () => true
  });
  return { actions, extractedFeatures };
};

// FE-09: a feature's single-feature condition is its stable hash, also when
// duplicate record IDs give several rendered instances one hash. It only
// prepares the matches of a color action; the rule "This feature only" writes
// is Python's and names the feature instance (PD-OI-069 revision 2,
// tests/test_web_rule_synthesis.py).
test('the single-feature condition is the stable hash, including for duplicate record IDs', () => {
  const first = feature('CDS', 'same_record_1');
  const second = feature('CDS', 'same_record_2');
  const otherType = feature('tRNA', 'same_record_3');
  const s = create([first, second, otherType]);
  assert.deepEqual(s.actions.getFeatureQualifier(first), { qual: 'hash', val: 'same' });
  assert.deepEqual(s.actions.getFeatureQualifier(second), { qual: 'hash', val: 'same' });
  assert.deepEqual(s.actions.getFeatureQualifier(otherType), { qual: 'hash', val: 'same' });
  s.extractedFeatures.value = [first];
  assert.deepEqual(s.actions.getFeatureQualifier(first), { qual: 'hash', val: 'same' });
});

test('a feature without a generation hash has no single-feature qualifier', () => {
  const s = create([]);
  assert.equal(s.actions.getFeatureQualifier({ type: 'CDS', svg_id: '' }), null);
});
