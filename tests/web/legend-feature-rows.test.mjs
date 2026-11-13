// The feature rows of a Legend and the row of each drawn feature: the Web's
// twin (`featureLegendRows`) answers the vector Python's owner answers
// (gbdraw/legend/table.py::FeatureLegendRows, tests/test_legend_row_facts.py).
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import test from 'node:test';
import { featureLegendRows, isTrackLegendKey, legendRowType } from '../../gbdraw/web/js/services/specific-color-rules.js';

const fixture = JSON.parse(readFileSync(new URL('../fixtures/legend_generated_rows.json', import.meta.url), 'utf8'));

test('each drawn feature is in the row Python draws for it', () => {
  for (const { name, rules, drawn, feature_rows: featureRows, type_rows: typeRows } of fixture.feature_rows) {
    const table = rules.map(([feat, cap, color]) => ({ feat, qual: 'product', val: '.*', color, cap }));
    const uses = drawn.map(([type, winner]) => ({ type, winner: winner === null ? null : table[winner] }));
    const rows = featureLegendRows(table, uses);
    const answer = (/** @type {ReturnType<typeof rows.rowOf>} */ row) => (row && typeof row === 'object' ? [row.cap, row.color] : row);
    assert.deepEqual(uses.map(({ type, winner }) => answer(rows.rowOf(type, winner))), featureRows, name);
    assert.deepEqual(Object.fromEntries(rows.present.map((type) => [type, rows.typeRow(type)])), typeRows, name);
    // A key is a type's "other" row exactly when the type draws its features of no rule in it (D-23).
    for (const type of rows.present) {
      const other = type === 'CDS' ? 'other proteins' : `other ${type}s`;
      assert.equal(rows.otherRowType(other), typeRows[type] === other ? type : null, `${name}: ${other}`);
      assert.equal(rows.otherRowType(type), null, `${name}: ${type} is no "other" row`);
    }
    // A captioned rule draws a row exactly when a drawn feature is in it.
    const drawnRows = new Set(uses.map(({ type, winner }) => JSON.stringify(answer(rows.rowOf(type, winner)))));
    for (const rule of table.filter((each) => each.cap)) {
      assert.equal(rows.drawsRow(rule), drawnRows.has(JSON.stringify([rule.cap, rule.color])), `${name}: ${rule.cap}`);
    }
  }
});

// The default rows Python generates (tests/test_legend_row_facts.py): a track's
// row, whose fill the track rows of the request give
// (tests/web/track-legend-row-colors.test.mjs), or a present type's own or
// "other" row, whose type's palette color fills it.
test('the generated Legend rows read from their keys are Python\'s', () => {
  for (const { name, features_present: types, rows } of fixture.cases) {
    for (const [key, paletteKey] of rows) {
      const track = !types.includes(paletteKey);
      assert.equal(isTrackLegendKey(key), track, `${name}: ${key} is a track row or a present type's row`);
      if (!track) assert.equal(legendRowType(key, types), paletteKey, `${name}: ${key}`);
    }
  }
  assert.equal(legendRowType('Group', ['CDS']), 'Group', 'any other key names a type\'s row');
});
