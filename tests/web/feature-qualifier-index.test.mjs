import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { test } from 'node:test';
import vm from 'node:vm';
import { createFeatureRuleActions } from '../../gbdraw/web/js/app/feature-editor/rule-actions.js';

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
    projectPaletteAndRules: () => true
  });
  return { actions, extractedFeatures };
};

// FE-09 (D-14, PD-OI-069): Python matches a color rule only by the stable hash
// (gbdraw/features/selector_values.py), so "This feature only" writes the
// stable hash even when duplicate record IDs give several rendered instances
// one hash. Those duplicates then share the rule.
test('This feature only writes the stable hash Python matches, including for duplicate record IDs', () => {
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
