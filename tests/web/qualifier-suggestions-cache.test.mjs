// OV-924: the Specific color rule qualifier suggestions name every qualifier
// of the drawn features. They are read once per feature list, from the raw
// rows, not again on each rule edit or mode switch.
import assert from 'node:assert/strict';
import test from 'node:test';
import { loadPageVue } from './helpers/vue-runtime.mjs';

const Vue = loadPageVue();
globalThis.window = { Vue, setTimeout, clearTimeout };
const { state } = await import('../../gbdraw/web/js/state.js');

test('qualifier suggestions scan the features once per feature list (OV-924)', () => {
  let qualifierReads = 0;
  const feature = (qualifiers) => new Proxy({ type: 'CDS', qualifiers }, {
    get(target, key, receiver) {
      if (key === 'qualifiers') qualifierReads += 1;
      return Reflect.get(target, key, receiver);
    }
  });
  state.extractedFeatures.value = [feature({ custom_tag: ['a'] }), feature({ VendorKey: ['b'], product: ['c'] })];
  const first = state.specificRuleQualifierSuggestions.value;
  assert.equal(first.includes('custom_tag'), true);
  assert.equal(first.includes('VendorKey'), true);
  const afterScan = qualifierReads;
  const rules = state.activeDrawing().manualSpecificRules;
  rules.push({ qual: 'rule_only' });
  const withRule = state.specificRuleQualifierSuggestions.value;
  assert.equal(withRule.includes('rule_only'), true);
  assert.equal(withRule.includes('custom_tag'), true);
  rules.pop();
  assert.equal(state.specificRuleQualifierSuggestions.value.includes('rule_only'), false);
  assert.equal(qualifierReads, afterScan, 'a rule edit does not rescan the features');
  state.extractedFeatures.value = [feature({ other_tag: ['x'] })];
  const replaced = state.specificRuleQualifierSuggestions.value;
  assert.equal(replaced.includes('other_tag'), true);
  assert.equal(replaced.includes('custom_tag'), false);
  state.extractedFeatures.value = [];
});
