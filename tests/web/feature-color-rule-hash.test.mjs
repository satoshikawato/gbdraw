import assert from 'node:assert/strict';
import { test } from 'node:test';
import {
  getFeatureColorRuleHash,
  getFeatureGenerationHash
} from '../../gbdraw/web/js/app/feature-utils.js';
import { ruleFeaturePayload } from '../../gbdraw/web/js/app/rule-matching.js';

const instanceId = 'ffa1f4c4a__instance_record_1_25cba05a594010cc';

test('the color rule hash drops record and instance suffixes (PD-OI-069)', () => {
  assert.equal(getFeatureColorRuleHash({ svg_id: 'ffa1f4c4a' }), 'ffa1f4c4a');
  assert.equal(getFeatureColorRuleHash({ svg_id: 'ffa1f4c4a_record_2' }), 'ffa1f4c4a');
  assert.equal(getFeatureColorRuleHash({ svg_id: instanceId }), 'ffa1f4c4a');
  assert.equal(
    getFeatureColorRuleHash({ rendered_feature_svg_id: instanceId, svg_id: 'other' }),
    'ffa1f4c4a'
  );
  assert.equal(getFeatureColorRuleHash(null), '');
  assert.equal(getFeatureColorRuleHash({}), '');
});

test('label and visibility callers keep the instance id', () => {
  assert.equal(getFeatureGenerationHash({ svg_id: instanceId }), instanceId);
});

test('the live rule payload uses the hash Python matches', () => {
  const payload = ruleFeaturePayload({ type: 'CDS', svg_id: instanceId, qualifiers: {} });
  assert.equal(payload.selector.hash, 'ffa1f4c4a');
});
