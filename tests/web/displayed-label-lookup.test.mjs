import assert from 'node:assert/strict';
import { test } from 'node:test';
import { createFeatureRuleActions } from '../../gbdraw/web/js/app/feature-editor/rule-actions.js';
import { withDrawings } from './helpers/drawing-state.mjs';

const setup = (features, entries) => {
  const state = {
    manualSpecificRules: [], featureColorOverrides: {}, featureOverrides: {}, labelTextBulkOverrides: {},
    extractedFeatures: { value: features }, editableLabels: { value: entries }
  };
  const never = () => { throw new Error('unused'); };
  return createFeatureRuleActions({
    ref: (value) => ({ value }), computed: (get) => ({ get value() { return get(); } }),
    state: withDrawings(state), rulePreparation: {}, runUndoable: never, runUndoableCheckpoint: never,
    prepareFileLegendEntries: never, projectPaletteAndRules: never, ports: {}, nextTick: async () => {}
  });
};

// An entry that counts every read of its featureId.
const countingEntry = (featureId, text, reads) => ({
  get featureId() { reads.count += 1; return featureId; },
  text
});

test('one same-displayed-label search reads each editable label entry once', () => {
  const N = 200;
  const E = 100;
  const features = Array.from({ length: N }, (_, i) => ({ id: `f${i}`, svg_id: `svg${i}`, type: 'CDS', gene: `g${i}`, start: i, end: i + 1 }));
  const reads = { count: 0 };
  const entries = Array.from({ length: E }, (_, i) => countingEntry(`svg${i * 2}`, i % 2 ? 'Shared' : `own ${i}`, reads));
  const actions = setup(features, entries);

  reads.count = 0;
  const found = actions.findFeaturesWithSameDisplayedLabel(features[2], 'Shared');
  assert.ok(reads.count <= E + N, `featureId reads: ${reads.count}`);
  assert.deepEqual(found.map((f) => f.svg_id), features.filter((_, i) => i % 4 === 2 && i !== 2).map((f) => f.svg_id));
});

test('the displayed label follows in-place edits and uses the first entry of a duplicated feature ID', () => {
  const features = [{ id: 'a', svg_id: 'svg-a', type: 'CDS', gene: 'geneA', start: 1, end: 2 }];
  const entries = [
    { featureId: ' SVG-A ', text: 'first' },
    { featureId: 'svg-a', text: 'second' }
  ];
  const actions = setup(features, entries);
  assert.equal(actions.getDisplayedFeatureLabel(features[0]), 'first');
  entries[0].text = 'edited';
  assert.equal(actions.getDisplayedFeatureLabel(features[0]), 'edited');
  entries[0].text = '';
  entries[0].featureId = 'other';
  assert.equal(actions.getDisplayedFeatureLabel(features[0]), 'second');
});
