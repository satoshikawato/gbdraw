// OV-555: a table keyed by file text or by an id is read by own key, so that
// `constructor` names no entry. On a Vue reactive table the own-key check alone
// (`Object.hasOwn`, `hasOwnProperty.call`) is no dependency: a computed that
// found no entry must still run again when the entry is written first.
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import test from 'node:test';
import vm from 'node:vm';

const context = vm.createContext({ console });
vm.runInContext(readFileSync(new URL('../../gbdraw/web/vendor/vue/vue.global.js', import.meta.url), 'utf8'), context);
const { computed, reactive } = context.Vue;
const { ownEntry } = await import('../../gbdraw/web/js/utils/own-entry.js');
const { featureOverrideKey, getFeatureOverride } = await import('../../gbdraw/web/js/services/feature-override-identity.js');

// One reader over one reactive table: the value before the entry exists, after
// its first write, after a change and after its removal.
const followed = (read, key, first, second) => {
  const table = reactive({});
  const shown = computed(() => read(table));
  const seen = [shown.value];
  table[key] = first;
  seen.push(shown.value);
  table[key] = second;
  seen.push(shown.value);
  delete table[key];
  seen.push(shown.value);
  return seen;
};

test('ownEntry inside a computed follows the first write, a change and the removal of a key', () => {
  // The expression that `ensureSearchIndex` (app/feature-search/preview-actions.js) uses for a label's sources.
  const sources = ['hypothetical protein', 'constructor'];
  const labels = (table) => sources.map((source) => ownEntry(table, source)).join('|');
  assert.deepEqual(
    followed(labels, 'hypothetical protein', 'Transporter', 'Permease'),
    ['|', 'Transporter|', 'Permease|', '|']
  );
});

test('only an own key names an entry, on a plain and on a reactive table', () => {
  for (const table of [{ gene: 'own' }, reactive({ gene: 'own' })]) {
    assert.equal(ownEntry(table, 'gene'), 'own');
    for (const inherited of ['constructor', 'toString', 'hasOwnProperty', '__proto__']) {
      assert.equal(ownEntry(table, inherited), undefined, inherited);
    }
  }
  assert.equal(ownEntry(reactive({ constructor: 'own text' }), 'constructor'), 'own text');
  assert.equal(ownEntry(null, 'gene'), undefined);
  assert.equal(ownEntry(undefined, 'gene'), undefined);
});

test('a feature override read follows the first write of the feature\'s entry', () => {
  const feature = { id: 'record-1\u0000cds-1', svg_id: 'f1', type: 'CDS' };
  const key = featureOverrideKey(feature);
  assert.ok(key);
  assert.deepEqual(
    followed((table) => getFeatureOverride(table, feature)?.color ?? 'palette', key, { color: '#ff0000' }, { color: '#00ff00' }),
    ['palette', '#ff0000', '#00ff00', 'palette']
  );
  assert.equal(getFeatureOverride(null, feature), undefined);
  assert.equal(getFeatureOverride(reactive({}), { id: 'constructor' }), undefined);
});
