// @ts-check
/** @import { LegendFontMetricsTable, LegendFontFaceTable } from '../utils/legend-font-metrics.generated.js' */
import { FALLBACK_FACE_PREFIX, bundledFacePrefix, comparableFamilyName, isGenericFamily, parseFontFamilyList } from './font-families.js';

// Python's Legend measurement and layout, ported operation for operation so
// that the Web places an edited Legend's rows, sizes the Legend and docks it
// exactly where Python would (zero shift). Sources:
// - text: gbdraw/core/text.py `drawn_text_box` and gbdraw/core/fonts.py
//   `resolve_font_face`;
// - rows: gbdraw/legend/metrics.py, linear_layout.py, circular_layout.py;
// - Legend size: gbdraw/configurators/legend.py `_linear_legend_local_bounds`,
//   `_circular_legend_local_bounds`;
// - docking and canvas: gbdraw/layout/composition.py `plan_composition`.
// The ports keep Python's min/max boxes and its order of operations, because
// doubles are equal only when computed the same way. The vectors in
// tests/fixtures/legend_layout_vectors.json hold Python's results; the node
// and Python tests require exact equality. A change to one side fails them.
// State-free: callers pass the font metrics (`loadLegendFontMetrics`, with
// the fonts registered in the page added by `withRegisteredFontFaces`).

const LEGEND_LINE_HEIGHT_RATIO = 24.0 / 14.0;
const LEGEND_TEXT_OFFSET_RATIO = 22.0 / 14.0;
const GRADIENT_BAR_WIDTH_RATIO = 10;
const GRADIENT_LABEL_GAP_RATIO = 0.2;
const SINGLE_GRADIENT_TRAILING_GAP_RATIO = 0.35;
// The similarity-ring block (gbdraw/legend/circular_layout.py, PD-OI-110).
const RING_OUTER_CAPTION = 'outer';
const RING_INNER_CAPTION = 'inner';
const RING_BAND_LABEL = 'identity';
const RING_BAND_LINE_GAP = 1.5;
const GLYPH_FIELDS = 5;
const KERNING_KEY_BASE = 65536;
const HORIZONTAL_SIDES = new Set(['top', 'bottom']);
const DOCK_SIDES = new Set(['left', 'right', 'top', 'bottom']);
const OVERLAY_SIDES = new Set(['upper_left', 'upper_right', 'lower_left', 'lower_right']);

/**
 * One measured face, indexed for lookups.
 * @typedef {object} LegendFontFace
 * @property {number} unitsPerEm
 * @property {ReadonlyMap<number, number>} glyphByCodePoint
 * @property {readonly number[]} glyphs Flat [advance, lsb, xMax, yMax, yMin, ...] by glyph index.
 * @property {readonly number[]} missingGlyph The facts of a character the face lacks: an em box
 *   ([unitsPerEm, 0, unitsPerEm, line ascent, -line descent]), as Python measures it.
 * @property {ReadonlyMap<number, number>} kerning Keyed `left * 65536 + right`.
 * @property {number} lineHeightEm The face's line box (line ascent + line descent) in em.
 */

/**
 * The measured faces, ready for `measureTextBox`: the bundled ones and those of
 * the fonts registered in the page.
 * @typedef {object} LegendFontMetrics
 * @property {ReadonlyMap<string, LegendFontFace>} faces By measurement name: a bundled face stem
 *   ("LiberationSans-Regular") or a registered face's `user-<hash>-<style>`.
 * @property {readonly string[]} bundledFaces
 * @property {ReadonlyMap<string, readonly RegisteredLegendFace[]>} registeredFamilies By comparable
 *   family name, in registration order.
 */

/**
 * A face of a font registered in the page, as Python registered it
 * (`register_font_file`, `font_metrics_table`).
 * @typedef {object} RegisteredLegendFace
 * @property {string} family
 * @property {number} weight The OS/2 weight class.
 * @property {boolean} italic
 * @property {string} measurementName
 * @property {LegendFontFaceTable} table
 */

/**
 * An axis-aligned box as Python's `Aabb` holds it.
 * @typedef {{ minX: number, minY: number, maxX: number, maxY: number }} LayoutBox
 */

/**
 * One Legend row of Python's legend table, in table order.
 * @typedef {object} LegendLayoutRow
 * @property {string} key The row key, which places the row's layout.
 * @property {string} name The name the row draws, which the layout measures (Python measures
 *   `drawn_names[key].drawn`); a renamed row keeps its key and draws its new name.
 * @property {'solid' | 'gradient' | 'ring'} type `ring`: a drawn similarity ring of the ring block.
 * @property {string} [stroke] The swatch stroke color, `none` for no stroke.
 * @property {number} [strokeWidth]
 * @property {number | string} [minValue] Gradient and ring rows: the lowest identity, which names the scale's
 *   first label.
 * @property {boolean} [swatch] Ring rows: whether the row shows its ring's color (the rings differ in color).
 */

/**
 * What Python's layout reads besides the rows (the `legendReflow` metadata of the Result).
 * @typedef {object} LegendLayoutOptions
 * @property {string} side The Legend side (`left`, `right`, `top`, `bottom`, or an overlay corner).
 * @property {number} wrapWidth The width Python passed to `measure_legend`.
 * @property {string} fontFile The measured face stem.
 * @property {number} fontSize
 * @property {number} colorRectSize
 */

/**
 * @typedef {object} LegendEntryLayout
 * @property {string} key
 * @property {LegendLayoutRow} row
 * @property {number} rectX
 * @property {number} rectY
 * @property {number} textX
 * @property {number} textY
 */

/** @param {unknown} value */
const text = (value) => String(value ?? '');

// ---- fonts ----

/**
 * Index a generated metrics table (`utils/legend-font-metrics.generated.js`).
 * @param {LegendFontMetricsTable} table
 * @returns {LegendFontMetrics}
 */
export const createLegendFontMetrics = (table) => {
  /** @type {Map<string, LegendFontFace>} */
  const faces = new Map();
  Object.entries(table.faces).forEach(([name, face]) => faces.set(name, indexFace(face)));
  return Object.freeze({ faces, bundledFaces: table.bundledFaces, registeredFamilies: new Map() });
};

/**
 * Index one face table.
 * @param {LegendFontFaceTable} face
 * @returns {LegendFontFace}
 */
const indexFace = (face) => {
  /** @type {Map<number, number>} */
  const glyphByCodePoint = new Map();
  for (let index = 0; index < face.cmap.length; index += 2) {
    glyphByCodePoint.set(face.cmap[index], face.cmap[index + 1]);
  }
  /** @type {Map<number, number>} */
  const kerning = new Map();
  for (let index = 0; index < face.kerning.length; index += 3) {
    kerning.set(face.kerning[index] * KERNING_KEY_BASE + face.kerning[index + 1], face.kerning[index + 2]);
  }
  const missingGlyph = Object.freeze([face.unitsPerEm, 0, face.unitsPerEm, face.lineAscent, -face.lineDescent]);
  // Python's `FontVerticalMetrics.ascent + .descent`, each divided by the em first.
  const lineHeightEm = (face.lineAscent / face.unitsPerEm) + (face.lineDescent / face.unitsPerEm);
  return Object.freeze({ unitsPerEm: face.unitsPerEm, glyphByCodePoint, glyphs: face.glyphs, kerning, missingGlyph, lineHeightEm });
};

/**
 * The bundled metrics with the faces of the fonts registered in the page, in
 * registration order (Python's registry order).
 * @param {LegendFontMetrics} metrics
 * @param {readonly RegisteredLegendFace[]} registered
 * @returns {LegendFontMetrics}
 */
export const withRegisteredFontFaces = (metrics, registered) => {
  if (registered.length === 0) return metrics;
  const faces = new Map(metrics.faces);
  /** @type {Map<string, RegisteredLegendFace[]>} */
  const registeredFamilies = new Map();
  registered.forEach((face) => {
    if (!faces.has(face.measurementName)) faces.set(face.measurementName, indexFace(face.table));
    // As Python's registry: a face of the same family, weight and style replaces the earlier one.
    const key = comparableFamilyName(face.family);
    const others = (registeredFamilies.get(key) || []).filter((other) => other.weight !== face.weight || other.italic !== face.italic);
    registeredFamilies.set(key, [...others, face]);
  });
  return Object.freeze({ faces, bundledFaces: metrics.bundledFaces, registeredFamilies });
};

/** @type {Promise<LegendFontMetrics> | null} */
let pendingFontMetrics = null;

/**
 * The bundled-font metrics, loaded on first use (about 320 KB of generated
 * module, 80 KB gzipped), so the app shell does not load them.
 * @returns {Promise<LegendFontMetrics>}
 */
export const loadLegendFontMetrics = () => {
  if (!pendingFontMetrics) {
    pendingFontMetrics = import('../utils/legend-font-metrics.generated.js')
      .then((module) => createLegendFontMetrics(module.LEGEND_FONT_METRICS))
      .catch((error) => {
        pendingFontMetrics = null;
        throw error;
      });
  }
  return pendingFontMetrics;
};

/**
 * Python's `collapse_svg_white_space` (gbdraw/core/text.py): a line break or
 * tab becomes a space and a run of spaces one space, with no trimming; NBSP and
 * the other Unicode spaces stay.
 * @param {string} value
 */
const collapseSvgWhiteSpace = (value) => value
  .replace(/\r\n/g, ' ')
  .replace(/[\t\n\r]/g, ' ')
  .replace(/ +/g, ' ');

/**
 * CSS font matching for weight 400 and style normal among a family's faces:
 * upright before italic; then 400 to 500 ascending, below 400 descending, above
 * 500 ascending. The Legend is measured with the normal weight and style.
 * @param {readonly RegisteredLegendFace[]} faces
 * @returns {RegisteredLegendFace | undefined}
 */
const normalFace = (faces) => {
  const upright = faces.filter((face) => !face.italic);
  const candidates = upright.length > 0 ? upright : faces;
  /** @param {number} weight */
  const rank = (weight) => (weight >= 400 && weight <= 500 ? [0, weight] : weight < 400 ? [1, -weight] : [2, weight]);
  return [...candidates].sort((left, right) => {
    const [leftClass, leftKey] = rank(left.weight);
    const [rightClass, rightKey] = rank(right.weight);
    return leftClass - rightClass || leftKey - rightKey;
  })[0];
};

/**
 * The face Python measures a font family with for the default weight and style
 * (`resolve_font_face` in a browser: the fonts registered in the page, then
 * bundled aliases; no fontconfig): the first family of the list that resolves,
 * else Liberation Sans.
 * @param {LegendFontMetrics} metrics
 * @param {string} fontFamily A CSS-style list ("'Liberation Sans', Arial, sans-serif").
 * @returns {string} A measurement name such as "LiberationSans-Regular".
 */
export const resolveLegendFontFace = (metrics, fontFamily) => {
  /** @param {string} prefix */
  const regular = (prefix) => {
    const stem = `${prefix}-Regular`;
    return metrics.bundledFaces.includes(stem) ? stem : null;
  };
  for (const family of parseFontFamilyList(fontFamily)) {
    // A family loaded in the page wins over the bundled font of the same alias name (OV-685).
    const registered = isGenericFamily(family) ? undefined : normalFace(metrics.registeredFamilies.get(comparableFamilyName(family)) || []);
    if (registered) return registered.measurementName;
    const prefix = bundledFacePrefix(family);
    const face = prefix ? regular(prefix) : null;
    if (face) return face;
  }
  return regular(FALLBACK_FACE_PREFIX) || metrics.bundledFaces[0];
};

/**
 * The face an edited Legend is measured with: the face Python recorded
 * (`legendReflow.fontFile`) when these metrics have it, else the face the
 * label's font family resolves to here (bundled or registered in this page).
 * A Result drawn on the CLI with a system font records a face the Web has no
 * metrics for (OV-683).
 * @param {LegendFontMetrics} metrics
 * @param {string | null | undefined} recordedFontFile
 * @param {string} fontFamily The label's `font-family`.
 * @returns {string}
 */
export const legendMeasurementFace = (metrics, recordedFontFile, fontFamily) => (
  recordedFontFile && metrics.faces.has(recordedFontFile)
    ? recordedFontFile
    : resolveLegendFontFace(metrics, fontFamily)
);

/**
 * The box text is drawn in, as gbdraw's Python layout measures it:
 * `drawn_text_box(text, family, size)` in gbdraw/core/text.py, with the face it
 * resolves to. Units are CSS px at `fontSize` px.
 *
 * The width is the advance sum with Python's kerning
 * (`calculate_svg_bbox_dimensions`, what a browser's `getComputedTextLength`
 * reports) or, when a glyph's ink ends past it, the right edge of the ink
 * (`text_ink`). The height is the face's line box, whatever the characters.
 * White space collapses as an SVG renderer draws it (`collapseSvgWhiteSpace`),
 * and a character the face lacks counts as an em box (OV-165).
 * @param {LegendFontMetrics} metrics From `loadLegendFontMetrics()` or `createLegendFontMetrics()`.
 * @param {{ text: string, fontFile: string, fontSize: number }} request
 *   `fontFile`: a measurement name in `metrics` (`resolveLegendFontFace` or `legendMeasurementFace`).
 * @returns {{ width: number, height: number }}
 */
export const measureTextBox = (metrics, { text: value, fontFile, fontSize }) => {
  const face = metrics.faces.get(fontFile);
  if (!face) throw new Error(`No font metrics for ${JSON.stringify(fontFile)}.`);
  const characters = Array.from(collapseSvgWhiteSpace(text(value)));
  if (characters.length === 0) return { width: 0.0, height: 0.0 };
  let origin = 0;
  let inkRight = -Infinity;
  /** @type {number | null} */
  let previous = null;
  characters.forEach((character) => {
    const glyph = face.glyphByCodePoint.get(/** @type {number} */ (character.codePointAt(0)));
    const offset = glyph === undefined ? -1 : glyph * GLYPH_FIELDS;
    const facts = offset < 0 ? face.missingGlyph : face.glyphs.slice(offset, offset + GLYPH_FIELDS);
    const [advance, inkLeft, xMax, yMax, yMin] = facts;
    if (previous !== null) {
      origin += (previous < 0 || glyph === undefined)
        ? 0
        : (face.kerning.get(previous * KERNING_KEY_BASE + glyph) || 0);
    }
    // Only a painted glyph has ink (`_face_text_ink`).
    if (xMax > inkLeft || yMax > yMin) inkRight = Math.max(inkRight, origin + xMax);
    origin += advance;
    previous = glyph === undefined ? -1 : glyph;
  });
  const size = Number(fontSize);
  const advanceWidth = origin * (size / face.unitsPerEm);
  const width = inkRight === -Infinity ? advanceWidth : Math.max(advanceWidth, (inkRight / face.unitsPerEm) * size);
  return { width, height: face.lineHeightEm * size };
};

// ---- shared ----

/**
 * CPython (3.12 and later) `sum()` of floats: Neumaier compensated summation
 * (Python/bltinmodule.c). Python's layout calls `sum`, so its port must too.
 * @param {readonly number[]} values
 * @returns {number}
 */
export const pythonFloatSum = (values) => {
  if (values.length === 0) return 0;
  let result = 0 + values[0];
  let compensation = 0.0;
  for (let index = 1; index < values.length; index += 1) {
    const value = values[index];
    const total = result + value;
    if (Math.abs(result) >= Math.abs(value)) compensation += (result - total) + value;
    else compensation += (value - total) + result;
    result = total;
  }
  if (compensation && Number.isFinite(compensation)) result += compensation;
  return result;
};

/** @param {number | string | undefined} value */
const minGradientLabelText = (value) => {
  const minIdentity = Number(value || 0);
  if (minIdentity === Math.trunc(minIdentity)) return `${Math.trunc(minIdentity)}%`;
  return `${String(minIdentity)}%`;
};

/**
 * @param {LegendFontMetrics} metrics
 * @param {LegendLayoutOptions} options
 */
const measurer = (metrics, { fontFile, fontSize }) => (
  /** @param {string} caption */ (caption) => measureTextBox(metrics, { text: caption, fontFile, fontSize })
);

/** @param {LayoutBox} box */
export const boxWidth = (box) => box.maxX - box.minX;
/** @param {LayoutBox} box */
export const boxHeight = (box) => box.maxY - box.minY;
/** @type {(minX: number, minY: number, maxX: number, maxY: number) => LayoutBox} */
const layoutBox = (minX, minY, maxX, maxY) => ({ minX, minY, maxX, maxY });
/** @type {(box: LayoutBox, dx: number, dy: number) => LayoutBox} */
const translatedBox = (box, dx, dy) => layoutBox(box.minX + dx, box.minY + dy, box.maxX + dx, box.maxY + dy);
/** @type {(box: LayoutBox, padding: number) => LayoutBox} */
const expandedBox = (box, padding) => layoutBox(box.minX - padding, box.minY - padding, box.maxX + padding, box.maxY + padding);
/** @param {readonly LayoutBox[]} boxes */
const unionBoxes = (boxes) => {
  let { minX, minY, maxX, maxY } = boxes[0];
  boxes.slice(1).forEach((box) => {
    minX = Math.min(minX, box.minX);
    minY = Math.min(minY, box.minY);
    maxX = Math.max(maxX, box.maxX);
    maxY = Math.max(maxY, box.maxY);
  });
  return layoutBox(minX, minY, maxX, maxY);
};

/** @param {LegendLayoutRow} row */
const halfStrokeWidth = (row) => {
  if (text(row.stroke ?? 'none').trim().toLowerCase() === 'none') return 0.0;
  const width = Number(row.strokeWidth || 0.0);
  if (!Number.isFinite(width) || width < 0) throw new Error('Legend stroke widths must be finite and non-negative.');
  return 0.5 * width;
};

// ---- Linear (gbdraw/legend/linear_layout.py) ----

/**
 * @typedef {object} LinearGradientEntryLayout
 * @property {string} key
 * @property {LegendLayoutRow} row
 * @property {number} titleX
 * @property {number} titleY
 * @property {number} barX
 * @property {number} barY
 */

/**
 * @typedef {object} LinearGradientLayout
 * @property {boolean} compact
 * @property {LinearGradientEntryLayout[]} entries
 * @property {number} width
 * @property {number} height
 * @property {number} barWidth
 * @property {string} minLabelText
 * @property {number} minLabelX
 * @property {number} maxLabelX
 * @property {number} scaleLabelY
 */

/**
 * @typedef {object} LinearOrientationLayout
 * @property {{ entries: LegendEntryLayout[], width: number, height: number, numLines: number }} feature
 * @property {LinearGradientLayout | null} gradient
 * @property {number} featureX
 * @property {number} featureY
 * @property {number} gradientX
 * @property {number} gradientY
 * @property {number} width
 * @property {number} height
 */

/**
 * @typedef {object} LinearLegendLayout
 * @property {LinearOrientationLayout} horizontal
 * @property {LinearOrientationLayout} vertical
 * @property {'horizontal' | 'vertical'} activeOrientation
 */

/**
 * Python's `single_gradient_block_x`: a one-gradient block is as wide as its
 * centred title, with the bar centred under it (OV-1261, OV-1640).
 * @param {number} barWidth
 * @param {number} titleWidth
 * @returns {{ width: number, barX: number }}
 */
const singleGradientBlockX = (barWidth, titleWidth) => {
  const width = Math.max(barWidth, titleWidth);
  return { width, barX: (width - barWidth) / 2.0 };
};

/**
 * @param {readonly LegendLayoutRow[]} rows
 * @param {(caption: string) => { width: number, height: number }} measure
 * @param {number} rect
 * @param {number} lineHeight
 * @returns {LinearGradientLayout | null}
 */
const linearGradientLayout = (rows, measure, rect, lineHeight) => {
  const gradients = rows.filter((row) => row.type === 'gradient');
  if (gradients.length === 0) return null;
  const barWidth = GRADIENT_BAR_WIDTH_RATIO * rect;
  const minLabelText = minGradientLabelText(gradients[0].minValue ?? 0);
  // The scale labels share one line box, whatever their digits.
  const scaleLabelHeight = measure('100%').height;
  if (gradients.length > 1) {
    const labelWidth = Math.max(...gradients.map((row) => measure(row.name).width));
    const barX = labelWidth + (GRADIENT_LABEL_GAP_RATIO * rect);
    const entries = gradients.map((row, index) => ({
      key: row.key,
      row,
      titleX: 0.0,
      titleY: (rect / 2.0) + (index * lineHeight),
      barX,
      barY: (rect / 2.0) + (index * lineHeight)
    }));
    const scaleLabelY = rect + ((gradients.length - 1) * lineHeight) + 2.0;
    return {
      compact: true, entries, width: barX + barWidth, height: scaleLabelY + scaleLabelHeight, barWidth,
      minLabelText, minLabelX: barX, maxLabelX: barX + barWidth, scaleLabelY
    };
  }
  const row = gradients[0];
  const title = measure(row.name);
  const { width, barX } = singleGradientBlockX(barWidth, title.width);
  const barY = title.height + (rect / 2.0);
  const scaleLabelY = barY + (rect / 2.0) + 2.0;
  return {
    compact: false,
    entries: [{ key: row.key, row, titleX: width / 2.0, titleY: 0.0, barX, barY }],
    width, height: scaleLabelY + scaleLabelHeight, barWidth,
    minLabelText, minLabelX: barX, maxLabelX: barX + barWidth, scaleLabelY
  };
};

/**
 * Python's Linear Legend layout, both orientations (`build_linear_legend_layout`).
 * @param {readonly LegendLayoutRow[]} rows
 * @param {LegendLayoutOptions} options
 * @param {LegendFontMetrics} metrics
 * @returns {LinearLegendLayout}
 */
export const buildLinearLegendLayout = (rows, options, metrics) => {
  const activeOrientation = HORIZONTAL_SIDES.has(options.side) ? 'horizontal' : 'vertical';
  const rect = options.colorRectSize;
  if (rows.length === 0) {
    /** @type {LinearOrientationLayout} */
    const empty = {
      feature: { entries: [], width: 0.0, height: 0.0, numLines: 0 }, gradient: null,
      featureX: 0.0, featureY: 0.0, gradientX: 0.0, gradientY: 0.0, width: 0.0, height: 0.0
    };
    return { horizontal: empty, vertical: empty, activeOrientation };
  }
  const lineHeight = LEGEND_LINE_HEIGHT_RATIO * rect;
  const offset = LEGEND_TEXT_OFFSET_RATIO * rect;
  const measure = measurer(metrics, options);
  const solids = rows.filter((row) => row.type === 'solid')
    .map((row) => ({ row, textWidth: measure(row.name).width }));
  const gradient = linearGradientLayout(rows, measure, rect, lineHeight);
  const gradientWidth = gradient ? gradient.width : 0.0;

  let wrapWidth = Number(options.wrapWidth);
  if (gradientWidth > 0.0 && wrapWidth > 0.0) {
    wrapWidth = Math.max(wrapWidth - gradientWidth - offset, rect + (2.0 * offset));
  }
  let yOffset = rect / 2.0;
  let currentX = 0.0;
  let currentRowWidth = 0.0;
  let maxRowWidth = 0.0;
  let height = lineHeight;
  let numLines = solids.length ? 1 : 0;
  /** @type {LegendEntryLayout[]} */
  const horizontalEntries = [];
  for (const { row, textWidth } of solids) {
    const entryWidth = rect + offset + textWidth + offset;
    if (wrapWidth > 0.0 && currentX + entryWidth > wrapWidth) {
      maxRowWidth = Math.max(maxRowWidth, currentRowWidth);
      currentX = 0.0;
      currentRowWidth = 0.0;
      yOffset += lineHeight;
      height += lineHeight;
      numLines += 1;
    }
    horizontalEntries.push({ key: row.key, row, rectX: currentX, rectY: yOffset, textX: currentX + offset, textY: yOffset });
    currentX += entryWidth;
    currentRowWidth += entryWidth;
  }
  const horizontalFeature = {
    entries: horizontalEntries, width: Math.max(maxRowWidth, currentRowWidth), height, numLines
  };

  let verticalY = rect / 2.0;
  let maxWidth = 0.0;
  /** @type {LegendEntryLayout[]} */
  const verticalEntries = [];
  for (const { row, textWidth } of solids) {
    verticalEntries.push({ key: row.key, row, rectX: 0.0, rectY: verticalY, textX: offset, textY: verticalY });
    maxWidth = Math.max(maxWidth, offset + textWidth);
    verticalY += lineHeight;
  }
  const verticalFeature = { entries: verticalEntries, width: maxWidth, height: verticalY, numLines: solids.length };

  if (!gradient) {
    return {
      horizontal: {
        feature: horizontalFeature, gradient: null, featureX: 0.0, featureY: 0.0, gradientX: 0.0, gradientY: 0.0,
        width: horizontalFeature.width, height: horizontalFeature.height
      },
      vertical: {
        feature: verticalFeature, gradient: null, featureX: 0.0, featureY: 0.0, gradientX: 0.0, gradientY: 0.0,
        width: verticalFeature.width, height: verticalFeature.height
      },
      activeOrientation
    };
  }
  return {
    horizontal: {
      feature: horizontalFeature,
      gradient,
      featureX: 0.0,
      featureY: Math.max(0.0, (gradient.height - horizontalFeature.height) / 2.0),
      gradientX: horizontalFeature.width + offset,
      gradientY: Math.max(0.0, (horizontalFeature.height - gradient.height) / 2.0),
      width: horizontalFeature.width + offset + gradient.width,
      height: Math.max(horizontalFeature.height, gradient.height)
    },
    vertical: {
      feature: verticalFeature,
      gradient,
      featureX: Math.max(0.0, (gradient.width - verticalFeature.width) / 2.0),
      featureY: 0.0,
      gradientX: 0.0,
      gradientY: verticalFeature.height + (lineHeight / 2.0),
      width: Math.max(verticalFeature.width, gradient.width),
      height: verticalFeature.height + (lineHeight / 2.0) + gradient.height
    },
    activeOrientation
  };
};

/**
 * The local bounds of the shown orientation (`_linear_legend_local_bounds`).
 * @param {LinearLegendLayout} layout
 * @param {number} colorRectSize
 * @returns {LayoutBox}
 */
export const linearLegendLocalBounds = (layout, colorRectSize) => {
  const active = layout[layout.activeOrientation];
  const rect = colorRectSize;
  const boxes = [layoutBox(0.0, 0.0, active.width, active.height)];
  for (const entry of active.feature.entries) {
    boxes.push(expandedBox(layoutBox(
      active.featureX + entry.rectX,
      active.featureY + entry.rectY - (0.5 * rect),
      active.featureX + entry.rectX + rect,
      active.featureY + entry.rectY + (0.5 * rect)
    ), halfStrokeWidth(entry.row)));
  }
  if (active.gradient) {
    for (const entry of active.gradient.entries) {
      boxes.push(expandedBox(layoutBox(
        active.gradientX + entry.barX,
        active.gradientY + entry.barY - (0.5 * rect),
        active.gradientX + entry.barX + active.gradient.barWidth,
        active.gradientY + entry.barY + (0.5 * rect)
      ), halfStrokeWidth(entry.row)));
    }
  }
  return unionBoxes(boxes);
};

// ---- Circular (gbdraw/legend/circular_layout.py) ----

/**
 * An entry of a compact block: its label, and its bar at the block's `barX`.
 * @typedef {object} CircularCompactGradientEntryLayout
 * @property {string} key
 * @property {LegendLayoutRow} row
 * @property {number} labelX
 * @property {number} labelY
 * @property {number} barY
 */

/**
 * The entry of a one-gradient block: its centred title and its bar.
 * @typedef {object} CircularSingleGradientEntryLayout
 * @property {string} key
 * @property {LegendLayoutRow} row
 * @property {number} titleX
 * @property {number} titleY
 * @property {number} barX
 * @property {number} barY
 */

/**
 * A gradient block: its entries, and the one scale (minimum and 100% labels)
 * below the bars that every form draws.
 * @typedef {object} CircularGradientLayout
 * @property {boolean} compact
 * @property {number} width
 * @property {number} height
 * @property {number} barWidth
 * @property {number} barX
 * @property {string} minLabelText
 * @property {number} minLabelX
 * @property {number} maxLabelX
 * @property {number} scaleLabelY
 * @property {CircularCompactGradientEntryLayout[]} compactEntries
 * @property {CircularSingleGradientEntryLayout[]} singleEntries
 */

/**
 * One drawn ring's row of the ring block: its swatch (when the block shows swatches) and name, centred on `y`.
 * @typedef {object} CircularRingRowLayout
 * @property {string} key
 * @property {LegendLayoutRow} row
 * @property {number} swatchX
 * @property {number} textX
 * @property {number} y
 */

/**
 * The similarity-ring block (Python's `CircularRingLegendLayout`): caption `outer`, the drawn rings outermost
 * first, caption `inner`, and the identity band labelled `identity` with its scale labels below it.
 * @typedef {object} CircularRingLayout
 * @property {number} width
 * @property {number} height
 * @property {boolean} showSwatches
 * @property {number} outerCaptionY
 * @property {number} innerCaptionY
 * @property {CircularRingRowLayout[]} rows
 * @property {number} barWidth
 * @property {number} bandY
 * @property {number} bandLabelX
 * @property {string} minLabelText
 * @property {number} minLabelX
 * @property {number} maxLabelX
 * @property {number} scaleLabelY
 */

/**
 * @typedef {object} CircularLegendLayout
 * @property {boolean} horizontal
 * @property {number} width
 * @property {number} height
 * @property {number} featureWidth
 * @property {number} featureHeight
 * @property {number} pairwiseLegendWidth
 * @property {number} lineMargin
 * @property {number} xMargin
 * @property {number} numLines
 * @property {number} numColumns
 * @property {number} numItemsPerLine
 * @property {LegendEntryLayout[]} entries
 * @property {CircularGradientLayout | null} gradient The gradient rows a saved Result drawn before the ring block
 *   still holds.
 * @property {CircularRingLayout | null} rings The similarity-ring block; it takes the gradient block's place.
 * @property {number} gradientX Where the gradient or ring block sits.
 * @property {number} gradientY
 */

/** @typedef {{ row: LegendLayoutRow, key: string, textWidth: number, entryWidth: number }} MeasuredCircularEntry */

/**
 * Python's `_build_ring_layout`.
 * @param {readonly MeasuredCircularEntry[]} rings
 * @param {(caption: string) => { width: number, height: number }} measure
 * @param {number} rect
 * @returns {CircularRingLayout | null}
 */
const circularRingLayout = (rings, measure, rect) => {
  if (rings.length === 0) return null;
  const line = LEGEND_LINE_HEIGHT_RATIO * rect;
  const xMargin = LEGEND_TEXT_OFFSET_RATIO * rect;
  const showSwatches = rings.some((entry) => Boolean(entry.row.swatch));
  const textX = showSwatches ? xMargin : 0.0;
  const barWidth = GRADIENT_BAR_WIDTH_RATIO * rect;
  const bandLabelX = barWidth + (xMargin - rect);
  const outerCaptionY = rect / 2;
  const rows = rings.map((entry, index) => ({
    key: entry.key, row: entry.row, swatchX: 0.0, textX, y: outerCaptionY + ((index + 1) * line)
  }));
  const innerCaptionY = outerCaptionY + ((rings.length + 1) * line);
  const bandY = innerCaptionY + (RING_BAND_LINE_GAP * line);
  const scaleLabelY = bandY + (rect / 2) + 2;
  const labelHeight = measure('100%').height;
  const captionWidth = Math.max(measure(RING_OUTER_CAPTION).width, measure(RING_INNER_CAPTION).width);
  const bandWidth = bandLabelX + measure(RING_BAND_LABEL).width;
  return {
    width: Math.max(captionWidth, bandWidth, ...rings.map((entry) => textX + entry.textWidth)),
    height: scaleLabelY + labelHeight + (line * SINGLE_GRADIENT_TRAILING_GAP_RATIO),
    showSwatches,
    outerCaptionY,
    innerCaptionY,
    rows,
    barWidth,
    bandY,
    bandLabelX,
    minLabelText: minGradientLabelText(rings[0].row.minValue ?? 0),
    minLabelX: 0.0,
    maxLabelX: barWidth,
    scaleLabelY
  };
};

/**
 * @param {readonly MeasuredCircularEntry[]} gradients
 * @param {(caption: string) => { width: number, height: number }} measure
 * @param {number} rect
 * @returns {CircularGradientLayout | null}
 */
const circularGradientLayout = (gradients, measure, rect) => {
  if (gradients.length === 0) return null;
  const barWidth = GRADIENT_BAR_WIDTH_RATIO * rect;
  const rowHeight = LEGEND_LINE_HEIGHT_RATIO * rect;
  const firstMinLabel = minGradientLabelText(gradients[0].row.minValue ?? 0);
  if (gradients.length > 1) {
    const labelWidth = Math.max(...gradients.map((entry) => entry.textWidth));
    const barX = labelWidth + rect;
    const scaleLabelY = rect + ((gradients.length - 1) * rowHeight) + 2;
    // The scale labels share one line box, whatever their digits.
    const labelHeight = measure('100%').height;
    return {
      compact: true,
      width: barX + barWidth,
      height: scaleLabelY + labelHeight,
      barWidth,
      barX,
      minLabelText: firstMinLabel,
      minLabelX: barX,
      maxLabelX: barX + barWidth,
      scaleLabelY,
      compactEntries: gradients.map((entry, index) => ({
        key: entry.key, row: entry.row, labelX: 0.0,
        labelY: (rect / 2) + (index * rowHeight), barY: (rect / 2) + (index * rowHeight)
      })),
      singleEntries: []
    };
  }
  const entry = gradients[0];
  const title = measure(entry.row.name);
  const barY = title.height + (rect / 2);
  const scaleLabelY = barY + (rect / 2) + 2;
  const labelHeight = measure('100%').height;
  const height = scaleLabelY + labelHeight + (rowHeight * SINGLE_GRADIENT_TRAILING_GAP_RATIO);
  const { width, barX } = singleGradientBlockX(barWidth, title.width);
  return {
    compact: false,
    width,
    height: Math.max(0.0, height),
    barWidth,
    barX,
    minLabelText: minGradientLabelText(entry.row.minValue ?? 0),
    minLabelX: barX,
    maxLabelX: barX + barWidth,
    scaleLabelY,
    compactEntries: [],
    singleEntries: [{ key: entry.key, row: entry.row, titleX: width / 2.0, titleY: 0.0, barX, barY }]
  };
};

/**
 * Python's Circular Legend layout (`build_circular_legend_layout`).
 * @param {readonly LegendLayoutRow[]} rows
 * @param {LegendLayoutOptions} options
 * @param {LegendFontMetrics} metrics
 * @returns {CircularLegendLayout}
 */
export const buildCircularLegendLayout = (rows, options, metrics) => {
  const rect = options.colorRectSize;
  const canvasWidth = Number(options.wrapWidth);
  const lineMargin = LEGEND_LINE_HEIGHT_RATIO * rect;
  const xMargin = LEGEND_TEXT_OFFSET_RATIO * rect;
  const measure = measurer(metrics, options);
  /** @type {MeasuredCircularEntry[]} */
  const solids = [];
  /** @type {MeasuredCircularEntry[]} */
  const gradients = [];
  /** @type {MeasuredCircularEntry[]} */
  const ringEntries = [];
  rows.forEach((row) => {
    if (row.type !== 'solid' && row.type !== 'gradient' && row.type !== 'ring') return;
    const textWidth = measure(row.name).width;
    const measured = { row, key: row.key, textWidth, entryWidth: textWidth + (2 * xMargin) };
    ({ solid: solids, gradient: gradients, ring: ringEntries })[row.type].push(measured);
  });
  const rings = circularRingLayout(ringEntries, measure, rect);
  const gradientRows = rings ? null : circularGradientLayout(gradients, measure, rect);
  // The block beside or below the solid rows: the ring block, or a saved Result's gradient rows.
  const gradient = rings || gradientRows;
  let pairwiseLegendWidth = 0.0;
  if (rings) {
    pairwiseLegendWidth = rings.width;
  } else if (gradients.length > 0) {
    const barWidth = GRADIENT_BAR_WIDTH_RATIO * rect;
    // One gradient reserves its block, title included (OV-1611).
    pairwiseLegendWidth = gradients.length > 1
      ? barWidth + (Math.max(...gradients.map((entry) => entry.textWidth)) + xMargin)
      : singleGradientBlockX(barWidth, gradients[0].textWidth).width;
  }
  const base = { pairwiseLegendWidth, lineMargin, xMargin, gradient: gradientRows, rings };

  if (HORIZONTAL_SIDES.has(options.side)) {
    let width = 0.0;
    const desiredSolid = pythonFloatSum(solids.map((entry) => entry.entryWidth)) + (solids.length ? 2 * xMargin : 0.0);
    let desired = desiredSolid;
    if (gradient) desired += Math.max(pairwiseLegendWidth, gradient.width);
    if (desired > 0) {
      let minWidth = solids.length ? Math.max(...solids.map((entry) => entry.entryWidth + (2 * xMargin))) : 0.0;
      if (gradient) {
        minWidth += Math.max(pairwiseLegendWidth, gradient.width);
        if (solids.length) minWidth += xMargin;
      }
      if (minWidth <= 0) minWidth = rect;
      width = Math.max(minWidth, Math.min(desired, canvasWidth));
    }
    let wrapWidth = width;
    if (gradient) wrapWidth = Math.max(rect, wrapWidth - pairwiseLegendWidth - xMargin);
    if (wrapWidth <= 0) wrapWidth = Infinity;
    /** @type {MeasuredCircularEntry[][]} */
    const rowsOfEntries = [];
    if (solids.length) {
      /** @type {MeasuredCircularEntry[]} */
      let current = [];
      let currentWidth = 0.0;
      const limit = wrapWidth <= 0 ? Infinity : wrapWidth;
      for (const entry of solids) {
        if (Number.isFinite(limit) && current.length && ((currentWidth + xMargin + entry.entryWidth) > limit)) {
          rowsOfEntries.push(current);
          current = [];
          currentWidth = 0.0;
        }
        current.push(entry);
        currentWidth += entry.entryWidth;
      }
      if (current.length) rowsOfEntries.push(current);
    }
    const featureHeight = rowsOfEntries.length ? rect + ((rowsOfEntries.length - 1) * lineMargin) : 0.0;
    const gradientHeight = gradient ? gradient.height : 0.0;
    const featureYOffset = Math.max(0.0, (gradientHeight - featureHeight) / 2.0);
    /** @type {LegendEntryLayout[]} */
    const entries = [];
    rowsOfEntries.forEach((rowEntries, rowIndex) => {
      const rowY = gradient ? featureYOffset + (rect / 2.0) + (rowIndex * lineMargin) : rowIndex * lineMargin;
      const rowWidth = pythonFloatSum(rowEntries.map((entry) => entry.entryWidth));
      const rowStartX = Number.isFinite(wrapWidth) ? xMargin + Math.max(0.0, (wrapWidth - rowWidth) * 0.5) : xMargin;
      let currentX = rowStartX;
      for (const entry of rowEntries) {
        entries.push({ key: entry.key, row: entry.row, rectX: currentX - xMargin, rectY: rowY, textX: currentX, textY: rowY });
        currentX += entry.entryWidth;
      }
    });
    const numLines = (solids.length || gradient) ? Math.max(1, rowsOfEntries.length) : 0;
    const numItemsPerLine = rowsOfEntries.reduce((most, rowEntries) => Math.max(most, rowEntries.length), 0);
    const featureWidth = Number.isFinite(wrapWidth) ? wrapWidth : width;
    let gradientX = 0.0;
    let gradientY = 0.0;
    let height = featureHeight;
    if (gradient) {
      gradientY = Math.max(0.0, (featureHeight - gradient.height) / 2.0);
      // Without feature rows the gradient block is the whole Legend (OV-1612).
      gradientX = Math.max(solids.length ? featureWidth + xMargin : 0.0, width - gradient.width);
      const contentBottom = Math.max(
        entries.length ? Math.max(...entries.map((entry) => entry.rectY)) + rect / 2.0 : 0.0,
        gradientY + gradient.height
      );
      height = contentBottom + (rect / 2.0);
    }
    return {
      ...base, horizontal: true, width, height, featureWidth, featureHeight, numLines,
      numColumns: numItemsPerLine + (gradient ? 1 : 0), numItemsPerLine, entries, gradientX, gradientY
    };
  }

  const featureBlockWidth = solids.length ? Math.max(...solids.map((entry) => xMargin + entry.textWidth)) : 0.0;
  const width = Math.max(featureBlockWidth, pairwiseLegendWidth, gradient ? gradient.width : 0.0);
  const alignmentWidth = Math.max(width, featureBlockWidth, gradient ? gradient.width : 0.0);
  const featureXOffset = gradient && featureBlockWidth > 0 ? Math.max(0.0, (alignmentWidth - featureBlockWidth) / 2.0) : 0.0;
  const entries = solids.map((entry, index) => ({
    key: entry.key, row: entry.row,
    rectX: featureXOffset, rectY: index * lineMargin, textX: featureXOffset + xMargin, textY: index * lineMargin
  }));
  const featureHeight = solids.length ? rect + ((solids.length - 1) * lineMargin) : 0.0;
  let gradientX = 0.0;
  let gradientY = 0.0;
  let height = featureHeight;
  if (gradient) {
    gradientX = Math.max(0.0, (alignmentWidth - gradient.width) / 2.0);
    gradientY = solids.length * lineMargin + (solids.length ? lineMargin * 0.5 : 0.0);
    height = Math.max(featureHeight, gradientY + gradient.height) + (rect / 2.0);
  }
  let numLines = solids.length;
  // outer, the rings, inner, the band and its scale labels
  if (rings) numLines += ringEntries.length + 4;
  else if (gradient) numLines += gradients.length + (gradients.length > 1 ? 1 : 2);
  return {
    ...base, horizontal: false, width, height, featureWidth: featureBlockWidth, featureHeight, numLines,
    numColumns: 1, numItemsPerLine: solids.length ? 1 : 0, entries, gradientX, gradientY
  };
};

/**
 * The Circular Legend's local bounds (`_circular_legend_local_bounds`).
 * @param {CircularLegendLayout} layout
 * @param {number} colorRectSize
 * @returns {LayoutBox}
 */
export const circularLegendLocalBounds = (layout, colorRectSize) => {
  const rect = colorRectSize;
  const boxes = [layoutBox(0.0, -0.5 * rect, layout.width, layout.height - (0.5 * rect))];
  for (const entry of layout.entries) {
    boxes.push(expandedBox(
      layoutBox(entry.rectX, entry.rectY - (0.5 * rect), entry.rectX + rect, entry.rectY + (0.5 * rect)),
      halfStrokeWidth(entry.row)
    ));
  }
  const rings = layout.rings;
  if (rings) {
    const x = layout.gradientX;
    const y = layout.gradientY;
    if (rings.showSwatches) {
      for (const row of rings.rows) {
        boxes.push(layoutBox(x + row.swatchX, y + row.y - (0.5 * rect), x + row.swatchX + rect, y + row.y + (0.5 * rect)));
      }
    }
    boxes.push(layoutBox(x, y + rings.bandY - (0.5 * rect), x + rings.barWidth, y + rings.bandY + (0.5 * rect)));
  }
  const gradient = layout.gradient;
  if (gradient) {
    for (const entry of [...gradient.compactEntries, ...gradient.singleEntries]) {
      boxes.push(expandedBox(layoutBox(
        layout.gradientX + gradient.barX,
        layout.gradientY + entry.barY - (0.5 * rect),
        layout.gradientX + gradient.barX + gradient.barWidth,
        layout.gradientY + entry.barY + (0.5 * rect)
      ), halfStrokeWidth(entry.row)));
    }
  }
  return unionBoxes(boxes);
};

// ---- composition (gbdraw/layout/composition.py) ----

/**
 * @typedef {object} CompositionSpacingPx
 * @property {number} edgePaddingPx
 * @property {number} dockGapPx
 * @property {number} titleGapPx
 * @property {number} stackGapPx
 * @property {number} overlayClearancePx
 */

/**
 * @typedef {object} CompositionOverlayPolicy
 * @property {readonly string[]} candidateScoreOrder
 * @property {readonly string[]} canvasGrowthCandidateOrder
 * @property {readonly string[]} canvasGrowthScoreOrder
 * @property {number} quadrantBoundaryRatio
 */

/**
 * @typedef {object} CompositionPlanRequest
 * @property {LayoutBox} primary
 * @property {LayoutBox | null} [legend]
 * @property {LayoutBox | null} [title]
 * @property {string} [legendSide]
 * @property {string} [titleSide]
 * @property {readonly LayoutBox[]} [overlayObstacles] In the primary item's coordinates.
 * @property {CompositionSpacingPx} spacing
 * @property {CompositionOverlayPolicy} overlayPolicy
 */

/** @typedef {{ role: string, dx: number, dy: number, local: LayoutBox, bounds: LayoutBox }} WorkingPlacement */

/**
 * @typedef {object} CompositionPlanResult
 * @property {LayoutBox} canvas From (0, 0).
 * @property {Array<{ role: string, translation: [number, number], finalBounds: LayoutBox }>} placements
 * @property {LayoutBox[]} overlayObstacles
 * @property {number[]} overlayConflictIndices
 * @property {string | null} overlayResolution `anchored`, `shifted`, `canvas_growth`, or null.
 */

/** @type {(role: string, local: LayoutBox, minX: number, minY: number) => WorkingPlacement} */
const alignMin = (role, local, minX, minY) => {
  const dx = minX - local.minX;
  const dy = minY - local.minY;
  return { role, dx, dy, local, bounds: translatedBox(local, dx, dy) };
};

/** @type {(box: LayoutBox, other: LayoutBox, clearance?: number) => boolean} */
const boxesIntersect = (box, other, clearance = 0.0) => (
  box.minX < other.maxX + clearance
  && box.maxX > other.minX - clearance
  && box.minY < other.maxY + clearance
  && box.maxY > other.minY - clearance
);

/** @param {readonly number[]} left @param {readonly number[]} right */
const compareTuples = (left, right) => {
  for (let index = 0; index < left.length; index += 1) {
    if (left[index] < right[index]) return -1;
    if (left[index] > right[index]) return 1;
  }
  return 0;
};

/** @type {(primaryMin: number, primaryMax: number, itemSize: number, nearMin: boolean, ratio: number) => [number, number] | null} */
const overlayAxisRange = (primaryMin, primaryMax, itemSize, nearMin, ratio) => {
  let minimum = primaryMin;
  let maximum = primaryMax - itemSize;
  const midpointStart = minimum + (ratio * (maximum - minimum));
  if (nearMin) maximum = Math.min(maximum, midpointStart);
  else minimum = Math.max(minimum, midpointStart);
  return minimum > maximum ? null : [minimum, maximum];
};

/** @type {(range: [number, number], itemSize: number, obstacles: readonly LayoutBox[], clearance: number, axis: 'x' | 'y') => number[]} */
const overlayAxisValues = ([minimum, maximum], itemSize, obstacles, clearance, axis) => {
  const values = new Set([minimum, maximum]);
  obstacles.forEach((obstacle) => {
    if (axis === 'x') {
      values.add(obstacle.minX - clearance - itemSize);
      values.add(obstacle.maxX + clearance);
    } else {
      values.add(obstacle.minY - clearance - itemSize);
      values.add(obstacle.maxY + clearance);
    }
  });
  return [...values].filter((value) => minimum <= value && value <= maximum).sort((a, b) => a - b);
};

/**
 * @param {LayoutBox} primary
 * @param {LayoutBox} legend
 * @param {string} side
 * @param {readonly LayoutBox[]} obstacles
 * @param {CompositionSpacingPx} spacing
 * @param {CompositionOverlayPolicy} policy
 * @returns {{ placement: WorkingPlacement, conflicts: number[], resolution: string }}
 */
const placeOverlayLegend = (primary, legend, side, obstacles, spacing, policy) => {
  const clearance = spacing.overlayClearancePx;
  const left = side === 'upper_left' || side === 'lower_left';
  const upper = side === 'upper_left' || side === 'upper_right';
  const anchorX = left ? primary.minX : primary.maxX - boxWidth(legend);
  const anchorY = upper ? primary.minY : primary.maxY - boxHeight(legend);
  const anchor = alignMin('legend', legend, anchorX, anchorY);
  /** @param {LayoutBox} bounds */
  const conflictsOf = (bounds) => obstacles
    .map((obstacle, index) => (boxesIntersect(bounds, obstacle, clearance) ? index : -1))
    .filter((index) => index >= 0);
  const initialConflicts = conflictsOf(anchor.bounds);
  const xRange = overlayAxisRange(primary.minX, primary.maxX, boxWidth(legend), left, policy.quadrantBoundaryRatio);
  const yRange = overlayAxisRange(primary.minY, primary.maxY, boxHeight(legend), upper, policy.quadrantBoundaryRatio);
  if (xRange && yRange) {
    if (initialConflicts.length === 0) return { placement: anchor, conflicts: [], resolution: 'anchored' };
    const xValues = overlayAxisValues(xRange, boxWidth(legend), obstacles, clearance, 'x');
    const yValues = overlayAxisValues(yRange, boxHeight(legend), obstacles, clearance, 'y');
    /** @param {[number, number]} coordinate */
    const score = ([x, y]) => {
      /** @type {Record<string, number>} */
      const metrics = {
        totalAnchorDistance: Math.abs(x - anchorX) + Math.abs(y - anchorY),
        xAnchorDistance: Math.abs(x - anchorX),
        yAnchorDistance: Math.abs(y - anchorY),
        nearEdgeX: left ? x : -x,
        nearEdgeY: upper ? y : -y
      };
      return policy.candidateScoreOrder.map((name) => metrics[name]);
    };
    /** @type {Array<[number, number]>} */
    const coordinates = [];
    xValues.forEach((x) => yValues.forEach((y) => coordinates.push([x, y])));
    const ordered = coordinates
      .map((coordinate) => ({ coordinate, key: score(coordinate) }))
      .sort((a, b) => compareTuples(a.key, b.key));
    for (const { coordinate } of ordered) {
      const candidate = alignMin('legend', legend, coordinate[0], coordinate[1]);
      if (conflictsOf(candidate.bounds).length === 0) {
        return { placement: candidate, conflicts: initialConflicts, resolution: 'shifted' };
      }
    }
  }
  const alignedX = left ? primary.minX : primary.maxX - boxWidth(legend);
  const alignedY = upper ? primary.minY : primary.maxY - boxHeight(legend);
  const horizontalX = left ? primary.minX - clearance - boxWidth(legend) : primary.maxX + clearance;
  const verticalY = upper ? primary.minY - clearance - boxHeight(legend) : primary.maxY + clearance;
  /** @type {Record<string, WorkingPlacement>} */
  const byName = {
    horizontal: alignMin('legend', legend, horizontalX, alignedY),
    vertical: alignMin('legend', legend, alignedX, verticalY)
  };
  const candidates = policy.canvasGrowthCandidateOrder.map((name) => byName[name]);
  const growthKey = (/** @type {WorkingPlacement} */ candidate, /** @type {number} */ index) => {
    const union = unionBoxes([primary, candidate.bounds]);
    /** @type {Record<string, number>} */
    const metrics = {
      addedArea: (boxWidth(union) * boxHeight(union)) - (boxWidth(primary) * boxHeight(primary)),
      addedExtent: (boxWidth(union) - boxWidth(primary)) + (boxHeight(union) - boxHeight(primary)),
      candidateOrder: index
    };
    return policy.canvasGrowthScoreOrder.map((name) => metrics[name]);
  };
  let best = 0;
  let bestKey = growthKey(candidates[0], 0);
  candidates.slice(1).forEach((candidate, offset) => {
    const key = growthKey(candidate, offset + 1);
    if (compareTuples(key, bestKey) < 0) {
      best = offset + 1;
      bestKey = key;
    }
  });
  return { placement: candidates[best], conflicts: initialConflicts, resolution: 'canvas_growth' };
};

/**
 * @param {LayoutBox} primary
 * @param {LayoutBox} title
 * @param {string} side
 * @param {CompositionSpacingPx} spacing
 * @param {WorkingPlacement | null} legend
 * @returns {WorkingPlacement}
 */
const placeTitle = (primary, title, side, spacing, legend) => {
  const centerX = 0.5 * (primary.minX + primary.maxX);
  const centerY = 0.5 * (primary.minY + primary.maxY);
  if (side === 'center') {
    return alignMin('title', title, centerX - (0.5 * boxWidth(title)), centerY - (0.5 * boxHeight(title)));
  }
  const titleMinX = centerX - (0.5 * boxWidth(title));
  let result;
  if (side === 'top') {
    result = alignMin('title', title, titleMinX, primary.minY - spacing.titleGapPx - boxHeight(title));
  } else if (side === 'bottom') {
    result = alignMin('title', title, titleMinX, primary.maxY + spacing.titleGapPx);
  } else {
    throw new Error(`Unsupported title placement ${JSON.stringify(side)}.`);
  }
  // Python's _place_title: a Legend the title meets, whatever its placement,
  // stays nearest the primary and the title takes the outer lane.
  if (legend && boxesIntersect(result.bounds, legend.bounds)) {
    result = alignMin(
      'title',
      title,
      titleMinX,
      side === 'top'
        ? legend.bounds.minY - spacing.stackGapPx - boxHeight(title)
        : legend.bounds.maxY + spacing.stackGapPx
    );
  }
  return result;
};

/** @param {LayoutBox | null | undefined} box */
const isEmptyBox = (box) => !box || boxWidth(box) <= 0.0 || boxHeight(box) <= 0.0;

/**
 * Python's `plan_composition`: dock or overlay the Legend, place the title,
 * pad the union, and move everything to start at (0, 0).
 * @param {CompositionPlanRequest} request
 * @returns {CompositionPlanResult}
 */
export const planLegendComposition = ({
  primary,
  legend = null,
  title = null,
  legendSide = 'none',
  titleSide = 'none',
  overlayObstacles = [],
  spacing,
  overlayPolicy
}) => {
  if (isEmptyBox(primary)) throw new Error('Primary bounds must have positive width and height.');
  /** @type {WorkingPlacement[]} */
  const working = [{ role: 'primary', dx: 0.0, dy: 0.0, local: primary, bounds: translatedBox(primary, 0.0, 0.0) }];
  /** @type {WorkingPlacement | null} */
  let legendPlacement = null;
  /** @type {number[]} */
  let overlayConflictIndices = [];
  /** @type {string | null} */
  let overlayResolution = null;
  if (legend && !isEmptyBox(legend) && legendSide !== 'none') {
    const primaryBounds = working[0].bounds;
    if (DOCK_SIDES.has(legendSide)) {
      const gap = spacing.dockGapPx;
      if (legendSide === 'left') {
        legendPlacement = alignMin('legend', legend, primaryBounds.minX - gap - boxWidth(legend),
          0.5 * (primaryBounds.minY + primaryBounds.maxY - boxHeight(legend)));
      } else if (legendSide === 'right') {
        legendPlacement = alignMin('legend', legend, primaryBounds.maxX + gap,
          0.5 * (primaryBounds.minY + primaryBounds.maxY - boxHeight(legend)));
      } else if (legendSide === 'top') {
        legendPlacement = alignMin('legend', legend, 0.5 * (primaryBounds.minX + primaryBounds.maxX - boxWidth(legend)),
          primaryBounds.minY - gap - boxHeight(legend));
      } else {
        legendPlacement = alignMin('legend', legend, 0.5 * (primaryBounds.minX + primaryBounds.maxX - boxWidth(legend)),
          primaryBounds.maxY + gap);
      }
    } else if (OVERLAY_SIDES.has(legendSide)) {
      const overlay = placeOverlayLegend(primaryBounds, legend, legendSide, overlayObstacles, spacing, overlayPolicy);
      legendPlacement = overlay.placement;
      overlayConflictIndices = overlay.conflicts;
      overlayResolution = overlay.resolution;
    } else {
      throw new Error(`Unsupported legend placement ${JSON.stringify(legendSide)}.`);
    }
    working.push(legendPlacement);
  }
  if (title && !isEmptyBox(title) && titleSide !== 'none') {
    working.push(placeTitle(working[0].bounds, title, titleSide, spacing, legendPlacement));
  }
  const padded = expandedBox(unionBoxes(working.map((placement) => placement.bounds)), spacing.edgePaddingPx);
  const outerDx = -padded.minX;
  const outerDy = -padded.minY;
  return {
    canvas: layoutBox(0.0, 0.0, boxWidth(padded), boxHeight(padded)),
    placements: working.map((placement) => ({
      role: placement.role,
      translation: /** @type {[number, number]} */ ([placement.dx + outerDx, placement.dy + outerDy]),
      finalBounds: translatedBox(placement.bounds, outerDx, outerDy)
    })),
    overlayObstacles: overlayObstacles.map((obstacle) => translatedBox(obstacle, outerDx, outerDy)),
    overlayConflictIndices,
    overlayResolution
  };
};
