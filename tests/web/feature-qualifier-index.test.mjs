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

// FE-09, PD-OI-069 revision 2: "This feature only" names the feature instance
// the catalog carries (Python's gbdraw/features/ids.py `feature_instance`)
// exactly, so a copy of the same record ID with the same stable hash keeps its
// color. A feature without an instance (one of a catalog written before
// catalogs carried it) is named by the stable hash those copies share.
test('This feature only names the feature instance, also for duplicate record IDs', () => {
  const instanceOf = (id, instance) => ({ ...feature('CDS', id), selector: { hash: 'same', instance } });
  const first = instanceOf('same_record_1', 'same@linear-seq-a');
  const second = instanceOf('same_record_2', 'same@linear-seq-b');
  const identical = instanceOf('same_record_3', 'same~1.5@record-1');
  const s = create([first, second, identical]);
  assert.deepEqual(s.actions.getFeatureQualifier(first), { qual: 'hash', val: '^same@linear-seq-a$' });
  assert.deepEqual(s.actions.getFeatureQualifier(second), { qual: 'hash', val: '^same@linear-seq-b$' });
  // The value is the exact pattern Python writes (`exact_value_pattern`).
  assert.deepEqual(s.actions.getFeatureQualifier(identical), { qual: 'hash', val: '^same~1\\.5@record-1$' });
  assert.deepEqual(s.actions.getFeatureQualifier(feature('CDS', 'same_record_1')), { qual: 'hash', val: 'same' });
});

test('a feature without a generation hash has no single-feature qualifier', () => {
  const s = create([]);
  assert.equal(s.actions.getFeatureQualifier({ type: 'CDS', svg_id: '' }), null);
});
