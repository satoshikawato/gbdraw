import assert from 'node:assert/strict';
import { test } from 'node:test';
import { getFeatureColorRuleHash } from '../../gbdraw/web/js/app/feature-utils.js';
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

test('the live rule payload uses the hash Python matches', () => {
  const payload = ruleFeaturePayload({ type: 'CDS', svg_id: instanceId, qualifiers: {} });
  assert.equal(payload.selector.hash, 'ffa1f4c4a');
});

// OV-02: catalog 5 gives the drawn record's selector values; a cropped record
// draws other coordinates than its source (design Q4 3.6).
test('the live rule payload sends the drawn selector values of catalog 5', () => {
  const payload = ruleFeaturePayload({
    type: 'misc_feature', svg_id: 'fb5977f81', record_id: 'TESTA', start: 2400, end: 2500, strand: '+',
    drawnSelector: { hash: 'f3b928d8c', location: '2200..2300', recordLocation: 'TESTA:2200..2300:+' },
    qualifiers: {}
  });
  assert.deepEqual(payload.selector, {
    hash: 'f3b928d8c', location: '2200..2300', record_location: 'TESTA:2200..2300:+'
  });
});

test('a catalog before schema 5 sends no coordinates for location rules (R4)', () => {
  const payload = ruleFeaturePayload({
    type: 'misc_feature', svg_id: 'f3b928d8c_record_1', record_id: 'TESTA', start: 2400, end: 2500,
    drawnSelector: null, qualifiers: {}
  });
  assert.deepEqual(payload.selector, { hash: 'f3b928d8c', location: null, record_location: null });
});

test('a catalog before schema 5 sends source coordinates where its rendered ID carries the source hash', () => {
  // The record was drawn with its source coordinates (no crop, reverse
  // complement, or rotation), so the source values are the drawn ones.
  const payload = ruleFeaturePayload({
    type: 'CDS', svg_id: 'ffa1f4c4a_record_1', record_id: 'TESTA', drawnSelector: null, qualifiers: {},
    selector: { hash: 'ffa1f4c4a', location: '300..600', record_location: 'TESTA:300..600:+', qualifiers: {} }
  });
  assert.deepEqual(payload.selector, { hash: 'ffa1f4c4a', location: '300..600', record_location: 'TESTA:300..600:+' });
});
