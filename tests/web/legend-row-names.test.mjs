import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { test } from 'node:test';
import { drawnLegendNames, keepBothNames, legendSwatch } from '../../gbdraw/web/js/services/legend-row-names.js';

// THE naming rule of the Legend (PD-OI-042 scenario revision 2, PD-OI-095):
// the Web app numbers the rows it displays as Python numbers the rows it
// draws (tests/test_legend_row_names.py answers the same cases).
const { cases } = JSON.parse(readFileSync(new URL('../fixtures/legend_names.json', import.meta.url), 'utf8'));

/** The fixture's swatch (Python's legend entry fields) as the Web app reads a row's paint. */
const paintOf = (swatch) => {
  if (swatch.type === 'ring') return { type: 'ring', trackIndex: Number(swatch.track_index) };
  return swatch.type === 'gradient'
    ? { type: 'gradient', minColor: swatch.min_color, maxColor: swatch.max_color }
    : { type: 'solid', fill: swatch.fill, stroke: swatch.stroke, width: swatch.width };
};

for (const { name, rows, drawn } of cases) {
  test(`drawn Legend names: ${name}`, () => {
    assert.deepEqual(drawnLegendNames(rows.map((row) => ({ name: row.name, swatch: legendSwatch(paintOf(row.swatch)) }))), drawn);
  });
}

test('a stroke width read from an SVG attribute is a number', () => {
  assert.equal(legendSwatch({ fill: 'red', stroke: 'gray', width: '2.0' }), legendSwatch({ fill: '#FF0000', stroke: 'Gray', width: 2 }));
  assert.equal(legendSwatch({ fill: 'red', stroke: 'gray', width: '' }), legendSwatch({ fill: '#f00', stroke: null, width: null }));
});

// PD-OI-061 scenario revision 5: Keep both keeps the typed name on the renamed
// row; every row of the name is numbered in legend order. The dialog names the
// renamed row's drawn name and, when the target's drawn name changes, the
// target's new drawn name.
const red = legendSwatch({ fill: '#ff0000' });
const blue = legendSwatch({ fill: '#0000ff' });
const green = legendSwatch({ fill: '#00ff00' });
const row = (key, name, swatch) => ({ key, name, swatch });

test('Keep both: the renamed row keeps its slot and both rows of the name are numbered', () => {
  const rows = [row('CDS', 'CDS', red), row('tRNA', 'tRNA', blue)];
  assert.deepEqual(keepBothNames(rows, { name: 'CDS', renamedKey: 'tRNA', swatch: blue, targetKey: 'CDS' }),
    { keepBothName: 'CDS (2)', targetNameAfter: 'CDS (1)' });
});

test('Keep both: a split-off row goes right after the target row (D-29)', () => {
  const rows = [row('CDS', 'CDS', red), row('tRNA', 'tRNA', blue)];
  assert.deepEqual(keepBothNames(rows, { name: 'CDS', renamedKey: null, swatch: green, targetKey: 'CDS' }),
    { keepBothName: 'CDS (2)', targetNameAfter: 'CDS (1)' });
  const before = [row('tRNA', 'tRNA', blue), row('CDS', 'CDS', red)];
  assert.deepEqual(keepBothNames(before, { name: 'CDS', renamedKey: null, swatch: green, targetKey: 'CDS' }),
    { keepBothName: 'CDS (2)', targetNameAfter: 'CDS (1)' });
});

test('Keep both: the renamed row before the target takes the first number', () => {
  const rows = [row('tRNA', 'tRNA', blue), row('CDS', 'CDS', red)];
  assert.deepEqual(keepBothNames(rows, { name: 'CDS', renamedKey: 'tRNA', swatch: blue, targetKey: 'CDS' }),
    { keepBothName: 'CDS (1)', targetNameAfter: 'CDS (2)' });
});

test('Keep both: a number skips a literal numbered name another row has', () => {
  const rows = [row('X (2)', 'X (2)', green), row('X', 'X', red), row('Y', 'Y', blue)];
  assert.deepEqual(keepBothNames(rows, { name: 'X', renamedKey: 'Y', swatch: blue, targetKey: 'X' }),
    { keepBothName: 'X (3)', targetNameAfter: 'X (1)' });
});

test('Keep both: a deleted or rule-only target draws no row, so no sentence', () => {
  const rows = [row('tRNA', 'tRNA', blue)];
  assert.deepEqual(keepBothNames(rows, { name: 'CDS', renamedKey: 'tRNA', swatch: blue, targetKey: null }),
    { keepBothName: 'CDS', targetNameAfter: null });
});

test('Keep both: an equal swatch is one row, drawn without a number (D-28)', () => {
  const rows = [row('CDS', 'CDS', red), row('tRNA', 'tRNA', red)];
  assert.deepEqual(keepBothNames(rows, { name: 'CDS', renamedKey: 'tRNA', swatch: red, targetKey: 'CDS' }),
    { keepBothName: 'CDS', targetNameAfter: null });
});

test('Keep both: names that differ only in case are drawn as written', () => {
  const rows = [row('CDS', 'CDS', red), row('tRNA', 'tRNA', blue)];
  assert.deepEqual(keepBothNames(rows, { name: 'cds', renamedKey: 'tRNA', swatch: blue, targetKey: 'CDS' }),
    { keepBothName: 'cds', targetNameAfter: null });
});

test('Keep both: a typed shown numbered name is kept as written (D-26)', () => {
  const rows = [row('CDS [red]', 'CDS', red), row('CDS [green]', 'CDS', green), row('tRNA', 'tRNA', blue)];
  assert.deepEqual(keepBothNames(rows, { name: 'CDS (2)', renamedKey: 'tRNA', swatch: blue, targetKey: 'CDS [green]' }),
    { keepBothName: 'CDS (2)', targetNameAfter: 'CDS (3)' });
});
