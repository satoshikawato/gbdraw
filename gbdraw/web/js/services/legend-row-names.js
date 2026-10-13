// @ts-check
import { resolveColorToHex } from '../utils/color-utils.js';

// THE naming rule of the Legend, for what the Web app displays (PD-OI-042
// scenario revision 2, PD-OI-095). A Legend row has three separate facts: its
// key (Python's identity of the row, `data-legend-key`, under which row edits
// and Sessions are stored; never shown), its name (what the user or the
// renderer wrote: `data-legend-name`, or a rename; never rewritten), and the
// name the Legend draws, which this module derives from the names and the
// swatches of the rows being drawn. Python owns the rule
// (gbdraw/legend/table.py::drawn_legend_names over ::legend_swatch); both
// answer tests/fixtures/legend_names.json.

/** @typedef {string & { readonly __brand: 'LegendRowName' }} LegendRowName The name of a Legend row as written; never a key, never numbered by gbdraw. */
/** @typedef {string & { readonly __brand: 'DrawnLegendName' }} DrawnLegendName The name the Legend draws for a row: its name, numbered when rows of other swatches share it. */
/** @typedef {string & { readonly __brand: 'LegendSwatchIdentity' }} LegendSwatchIdentity What two Legend rows compare to be the same swatch. */

/**
 * The paint of a Legend row's swatch as it is drawn: a solid swatch's fill
 * and stroke, or a gradient's two end colors.
 * @typedef {{ type?: 'solid', fill?: unknown, stroke?: unknown, width?: unknown }
 *   | { type: 'gradient', minColor?: unknown, maxColor?: unknown }} LegendSwatchPaint
 */

/**
 * A row as the naming rule reads it.
 * @typedef {object} NamedLegendRow
 * @property {LegendRowName} name
 * @property {LegendSwatchIdentity} swatch
 */

// Python's `_legend_fill_identity`: `none`, lower-case #rrggbb (a color name
// equals its hex), else the text in lower case (a hatch's pattern paint).
/** @param {unknown} value @returns {string} */
export const legendColorIdentity = (value) => {
  const color = String(resolveColorToHex(String(value ?? '').trim()) || '').toLowerCase();
  return /^#[0-9a-f]{3}$/.test(color) ? `#${[...color.slice(1)].map((char) => char + char).join('')}` : color;
};

// Python's `legend_swatch`: a solid row is its fill and its stroke (color and
// width; a stroke of no width or no color is no stroke), a hatch being part of
// the fill; a gradient row is its two end colors.
/** @param {LegendSwatchPaint} paint @returns {LegendSwatchIdentity} */
export const legendSwatch = (paint) => {
  if (paint.type === 'gradient') {
    return /** @type {LegendSwatchIdentity} */ (JSON.stringify(['gradient', legendColorIdentity(paint.minColor), legendColorIdentity(paint.maxColor)]));
  }
  const width = Number(paint.width);
  const stroke = legendColorIdentity(paint.stroke) || 'none';
  const drawn = Number.isFinite(width) && width > 0 && stroke !== 'none';
  return /** @type {LegendSwatchIdentity} */ (JSON.stringify(['solid', legendColorIdentity(paint.fill), drawn ? stroke : 'none', drawn ? width : 0]));
};

/**
 * The name each row is drawn with, for the rows being drawn in legend order.
 * Rows of one name and one swatch are one row and share a drawn name. A name
 * left with one row is drawn as written; a name left with two or more rows is
 * drawn `name (1)`, `name (2)`, ... in legend order, and a number skips a
 * drawn name that another row already has.
 * @param {readonly NamedLegendRow[]} rows
 * @returns {DrawnLegendName[]}
 */
export const drawnLegendNames = (rows) => {
  /** @type {Map<LegendRowName, LegendSwatchIdentity[]>} */
  const swatches = new Map();
  rows.forEach(({ name, swatch }) => {
    const group = swatches.get(name) || [];
    if (!group.includes(swatch)) group.push(swatch);
    swatches.set(name, group);
  });
  /** @type {Set<string>} */
  const taken = new Set([...swatches].filter(([, group]) => group.length === 1).map(([name]) => name));
  /** @type {Map<string, DrawnLegendName>} */
  const numbered = new Map();
  /** @param {NamedLegendRow} row */
  const rowKey = ({ name, swatch }) => JSON.stringify([name, swatch]);
  swatches.forEach((group, name) => {
    if (group.length === 1) return;
    let number = 0;
    group.forEach((swatch) => {
      number += 1;
      while (taken.has(`${name} (${number})`)) number += 1;
      const drawn = /** @type {DrawnLegendName} */ (`${name} (${number})`);
      numbered.set(rowKey({ name, swatch }), drawn);
      taken.add(drawn);
    });
  });
  return rows.map((row) => numbered.get(rowKey(row)) || /** @type {DrawnLegendName} */ (/** @type {string} */ (row.name)));
};
