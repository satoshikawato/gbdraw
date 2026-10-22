// OV-555: a table keyed by file text or by an id is read by own key, so that
// `constructor` names no entry. On a Vue reactive table the own-key check alone
// (`Object.hasOwn`, `hasOwnProperty.call`) is no dependency: a computed that
// found no entry must still run again when the entry is written first.
import assert from 'node:assert/strict';
import test from 'node:test';
import { loadPageVue } from './helpers/vue-runtime.mjs';

const { computed, reactive } = loadPageVue();
const { hasOwnEntry, ownEntry } = await import('../../gbdraw/web/js/utils/own-entry.js');
const { featureOverrideKey, getFeatureOverride } = await import('../../gbdraw/web/js/services/feature-override-identity.js');
const { draftLegendPanelColors } = await import('../../gbdraw/web/js/app/candidate-render.js');
const { createFeatureRuleActions } = await import('../../gbdraw/web/js/app/feature-editor/rule-actions.js');
const { withDrawings } = await import('./helpers/drawing-state.mjs');

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

// OV-556: a reader that asks whether a key is there (and not for its value).
test('hasOwnEntry: only an own key is there, also one that holds undefined', () => {
  for (const table of [{ gene: 'own', empty: undefined }, reactive({ gene: 'own', empty: undefined })]) {
    assert.equal(hasOwnEntry(table, 'gene'), true);
    assert.equal(hasOwnEntry(table, 'empty'), true, 'an own key that holds undefined is there');
    for (const absent of ['constructor', 'toString', 'hasOwnProperty', '__proto__', 'product', '']) {
      assert.equal(hasOwnEntry(table, absent), false, absent);
    }
  }
  assert.equal(hasOwnEntry(reactive({ constructor: 'own text' }), 'constructor'), true);
  assert.equal(hasOwnEntry(null, 'gene'), false);
  assert.equal(hasOwnEntry(undefined, 'gene'), false);
});

test('hasOwnEntry inside a computed follows the first write, a change and the removal of a key', () => {
  assert.deepEqual(
    followed((table) => hasOwnEntry(table, 'constructor'), 'constructor', 'Renamed', undefined),
    [false, true, true, false]
  );
});

test('the Legend panel color of a row follows the first color edit of that row (OV-556)', () => {
  // A row the palette colors shows the palette's color until the row gets its own.
  const legend = {
    legendEntries: [{ caption: 'CDS', originalCaption: 'CDS', color: '#84b9ec' }], originalLegendOrder: ['CDS'],
    pythonRows: new Map([['CDS', { key: 'CDS', color: '#84b9ec' }]]), rules: [], paletteColors: { CDS: '#112233' },
    drawnTracks: []
  };
  const paletteColored = (legendColorOverrides) => (
    draftLegendPanelColors({ ...legend, legendColorOverrides }).listed.get('CDS') ?? 'the row\'s own color'
  );
  assert.deepEqual(
    followed(paletteColored, 'CDS', '#ff0000', '#00ff00'),
    ['#112233', 'the row\'s own color', 'the row\'s own color', '#112233']
  );
  assert.equal(paletteColored({}), '#112233');
  assert.equal(paletteColored(null), '#112233');
});

test('the rendering of a feature type follows the first write of its entry (OV-556)', () => {
  const adv = reactive({ feature_shapes: {} });
  const state = withDrawings({ adv, manualSpecificRules: [] });
  const { getFeatureShape } = createFeatureRuleActions({
    state, ref: (value) => ({ value }), computed: (get) => ({ get value() { return get(); } })
  });
  // `constructor` is a feature type like any other: the type's default, then its own entry.
  const shown = computed(() => ['repeat_region', 'constructor'].map(getFeatureShape));
  assert.deepEqual(shown.value, ['underlay', 'rectangle']);
  adv.feature_shapes.repeat_region = 'arrow';
  assert.deepEqual(shown.value, ['arrow', 'rectangle']);
  adv.feature_shapes.constructor = 'underlay';
  assert.deepEqual(shown.value, ['arrow', 'underlay']);
  delete adv.feature_shapes.repeat_region;
  assert.deepEqual(shown.value, ['underlay', 'underlay']);
  // An own key is read as it is: a value that is no rendering is an error, as before (also `undefined`).
  adv.feature_shapes.CDS = undefined;
  assert.throws(() => getFeatureShape('CDS'), /Unsupported feature rendering 'undefined'/);
  adv.feature_shapes = null;
  assert.equal(getFeatureShape('CDS'), 'arrow');
});
