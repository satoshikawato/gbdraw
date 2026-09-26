// Synthetic baseline evidence. Does not initialize a Worker or mutate the app.
import assert from 'node:assert/strict';
import { STANDALONE_INTERACTIVE_SCRIPT } from '../../../../gbdraw/web/js/services/standalone-interactivity-assets.js';
import { runFeatureSearch } from '../../../../gbdraw/web/js/app/feature-search/search-core.js';
import { resolveEffectiveFeatureVisibility } from '../../../../gbdraw/web/js/app/feature-visibility.js';
import { createFeatureRuleActions } from '../../../../gbdraw/web/js/app/feature-editor/rule-actions.js';

const patterns = ['(?i)NADH', '(?P<enzyme>NADH)', '[', '(?<enzyme>NADH)', String.raw`NADH\Z`, String.raw`\bβ`, 'i', 'β'];
const corpus = ['NADH', 'nadh', 'β-lactamase', 'ı', 'i', 'İ', 'unrelated'];
// Execute the exact standalone regex branch, extracted by fixed source anchors.
// This is a function probe, not standalone document/browser initialization.
const start = STANDALONE_INTERACTIVE_SCRIPT.indexOf('  function compileSearchMatcher(query, useRegex) {');
const end = STANDALONE_INTERACTIVE_SCRIPT.indexOf('\n    var needle = normalizeSearchText(trimmedQuery);', start);
assert.ok(start >= 0 && end > start);
const standaloneCompile = new Function(`${STANDALONE_INTERACTIVE_SCRIPT.slice(start, end)}\n}\nreturn compileSearchMatcher;`)();
const rows = patterns.map((pattern) => {
  const standalone = standaloneCompile(pattern, true);
  const targets = corpus.flatMap((text, index) => standalone.test([text]) ? [index] : []);
  const search = runFeatureSearch({ features: [], query: pattern, useRegex: true });
  assert.equal(search.error, standalone.error);
  return { pattern, emptySearchError: search.error, standaloneError: standalone.error, standaloneCorpusTargets: targets };
});
const visibilityFallback = ['(?i)^NADH$', '(?P<enzyme>NADH)', '(?<enzyme>NADH)', '[', '^NADH$', '^nadh$'].map((pattern) => ({
  pattern,
  result: resolveEffectiveFeatureVisibility('NADH', {}, null, [{ recordId: '*', featureType: '*', qualifier: 'hash', value: pattern, action: 'hide' }])
}));
// Inject a non-syntax preparation rejection into the real field action.
// History/SVG are stubs; this does not demonstrate an actual Worker failure.
const accepted = { feat: 'CDS', qual: 'product', val: 'NADH', color: '#000000', cap: '' };
const state = { manualSpecificRules: [accepted], newSpecRule: {}, specificRulePresets: [], manualPriorityRules: [], adv: {} };
const alerts = [];
const previousAlert = globalThis.alert;
globalThis.alert = (message) => alerts.push(message);
const input = { isConnected: true, value: '[' };
try {
  const actions = createFeatureRuleActions({ state, nextTick: async () => {}, legendActions: {},
    rulePreparation: { snapshot: () => ({}), isCurrent: () => true, prepare: () => Promise.reject(new Error('Synthetic runtime unavailable')) },
    history: { runUndoable: () => { throw new Error('Unexpected commit'); } }, svgActions: {} });
  await actions.setSpecificRuleField(0, 'val', '[', input);
  assert.equal(state.manualSpecificRules[0], accepted);
  assert.equal(input.value, 'NADH');
  assert.deepEqual(alerts, ['Invalid rule: Synthetic runtime unavailable']);
} finally {
  if (previousAlert === undefined) delete globalThis.alert;
  else globalThis.alert = previousAlert;
}
console.log(JSON.stringify({ scope: 'function probes; empty app-search catalog, standalone regex branch, cache-miss hash fallback, stubbed rule preparation rejection', corpus, rows, visibilityFallback, rejectedField: { alerts, displayedValue: input.value, acceptedValue: accepted.val } }, null, 2));
