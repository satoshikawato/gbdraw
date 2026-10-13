import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { test } from 'node:test';
import { drawnLegendNames, legendSwatch } from '../../gbdraw/web/js/services/legend-row-names.js';

// THE naming rule of the Legend (PD-OI-042 scenario revision 2, PD-OI-095):
// the Web app numbers the rows it displays as Python numbers the rows it
// draws (tests/test_legend_row_names.py answers the same cases).
const { cases } = JSON.parse(readFileSync(new URL('../fixtures/legend_names.json', import.meta.url), 'utf8'));

/** The fixture's swatch (Python's legend entry fields) as the Web app reads a row's paint. */
const paintOf = (swatch) => (swatch.type === 'gradient'
  ? { type: 'gradient', minColor: swatch.min_color, maxColor: swatch.max_color }
  : { type: 'solid', fill: swatch.fill, stroke: swatch.stroke, width: swatch.width });

for (const { name, rows, drawn } of cases) {
  test(`drawn Legend names: ${name}`, () => {
    assert.deepEqual(drawnLegendNames(rows.map((row) => ({ name: row.name, swatch: legendSwatch(paintOf(row.swatch)) }))), drawn);
  });
}

test('a stroke width read from an SVG attribute is a number', () => {
  assert.equal(legendSwatch({ fill: 'red', stroke: 'gray', width: '2.0' }), legendSwatch({ fill: '#FF0000', stroke: 'Gray', width: 2 }));
  assert.equal(legendSwatch({ fill: 'red', stroke: 'gray', width: '' }), legendSwatch({ fill: '#f00', stroke: null, width: null }));
});
