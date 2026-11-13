// @ts-check
// Select region: from a point on a drawn record to the bases it shows. Python
// writes each record's axis into its track-slot geometry
// (`gbdraw/render/track_slot_metadata.py`, `axis`); this module only reads it.
// Points are in the frame of the record's drawn group (`regionAxisFrame`):
// the Circular centre is (0, 0), a Linear row starts at x 0 on its axis line.
import { ANNOTATION_SELECTOR } from './drawn-groups.js';
import { numberText, t } from '../utils/i18n.js';

/**
 * Which source base a drawn base shows. Display index `i` (0 to `drawnLength - 1` along the axis) shows source base
 * `(anchorBase - 1 + step * i) mod sourceLength + 1`; the record's own base `k` (its feature coordinates, 0-based)
 * is source base `(recordBase - 1 + step * k) mod sourceLength + 1`. `isCircular` says whether a span may cross the
 * record's origin (Python's `spans_may_cross_origin`: a circular record with no region applied, the annotation
 * resolver's rule; such a record draws every base); a drag crosses display index 0 only there.
 * @typedef {object} AxisSourceMap
 * @property {number} anchorBase
 * @property {number} recordBase
 * @property {1 | -1} step
 * @property {number} sourceLength
 * @property {number} drawnLength
 * @property {boolean} isCircular
 */

/**
 * @typedef {object} CircularAxisGeometry
 * @property {'circular'} kind
 * @property {number} radiusPx
 * @property {number} bandInnerPx
 * @property {number} bandOuterPx
 * @property {number} startAngleDeg The angle of display index 0, in SVG degrees (clockwise from +x).
 * @property {number} degreesPerBp
 * @property {AxisSourceMap} source
 */

/**
 * @typedef {object} LinearAxisRow
 * @property {number} xStartPx
 * @property {number} xEndPx
 * @property {number} yPx
 * @property {number} bandTopPx
 * @property {number} bandBottomPx
 * @property {number} displayStart
 * @property {number} displayEnd Exclusive.
 */

/**
 * @typedef {object} LinearAxisGeometry
 * @property {'linear'} kind
 * @property {LinearAxisRow[]} rows
 * @property {AxisSourceMap} source
 */

/** @typedef {CircularAxisGeometry | LinearAxisGeometry} RegionAxis */

/**
 * One drawn record that a region can be selected on: a record drawn whole or in one window, or one window of a record
 * drawn in several. Each window of such a record names the record's key and its window index (D-01), so a drawn unit
 * is found by `recordKey` and `windowIndex`, never by its place in the Result.
 * @typedef {object} RegionAxisRecord
 * @property {number} resultIndex
 * @property {number} recordIndex The drawn unit's index in the Result.
 * @property {string} recordKey The key of the request record the unit draws.
 * @property {number | null} windowIndex The unit's window of a record drawn in several; null for any other unit.
 * @property {string} recordId
 * @property {RegionAxis} axis
 */

/**
 * A run of drawn bases: `length` bases from display index `start`, forward along the axis (clockwise on a circle,
 * where it may pass display index 0 on a whole record).
 * @typedef {{ start: number, length: number }} DisplaySpan
 */

/**
 * A span in the input file's coordinates, 1-based with both ends included; `start > end` when it crosses the origin.
 * @typedef {{ start: number, end: number, length: number, wrapsOrigin: boolean }} SourceSpan
 */

/**
 * A selected span, as the Region popup holds it: which record of which Result (its window, of a record drawn in
 * several), the drawn bases, and their source coordinates.
 * @typedef {object} RegionSelection
 * @property {'region'} kind
 * @property {number} resultIndex
 * @property {number} recordIndex
 * @property {string} recordKey
 * @property {number | null} windowIndex
 * @property {string} recordId
 * @property {DisplaySpan} span
 * @property {SourceSpan} source
 */

/** @typedef {{ x: number, y: number }} FramePoint */
/** @typedef {{ start: number, end: number }} RecordPart A feature part in the input file's coordinates (the catalog's `location_parts`), 0-based, end excluded. */

/**
 * A drawn feature of the feature catalog as the converter reads it: its record key, the drawn unit it lies in, and its
 * parts in the input file's coordinates (fields of a row of `extractedFeatures`, an `ExtractedFeature`, read by name).
 * @typedef {object} RegionFeatureRow
 * @property {unknown} [svg_id]
 * @property {unknown} [record_key]
 * @property {unknown} [drawn_record_index] The unit's index in the Result, as `RegionAxisRecord.recordIndex`.
 * @property {unknown} [start]
 * @property {unknown} [end]
 * @property {unknown} [location_parts]
 */

/** @param {unknown} value @returns {value is Record<string, unknown>} */
const isObject = (value) => Boolean(value) && typeof value === 'object' && !Array.isArray(value);
/** @param {unknown} value @returns {value is number} */
const finite = (value) => typeof value === 'number' && Number.isFinite(value);
/** @param {unknown} value @returns {value is number} */
const whole = (value) => Number.isSafeInteger(value);
/** @param {number} value @param {number} modulus */
const mod = (value, modulus) => ((value % modulus) + modulus) % modulus;

/** @param {unknown} value @returns {AxisSourceMap | null} */
const readSource = (value) => {
  if (!isObject(value)) return null;
  const { anchorBase, recordBase, step, sourceLength, drawnLength, isCircular } = value;
  if (!whole(sourceLength) || sourceLength < 1 || !whole(drawnLength) || drawnLength < 1 || drawnLength > sourceLength) {
    return null;
  }
  if (!whole(anchorBase) || !whole(recordBase) || anchorBase < 1 || recordBase < 1) return null;
  if (anchorBase > sourceLength || recordBase > sourceLength || (step !== 1 && step !== -1)) return null;
  if (typeof isCircular !== 'boolean') return null;
  return { anchorBase, recordBase, step, sourceLength, drawnLength, isCircular };
};

/** @param {unknown} value @returns {LinearAxisRow | null} */
const readRow = (value) => {
  if (!isObject(value)) return null;
  const { xStartPx, xEndPx, yPx, bandTopPx, bandBottomPx, displayStart, displayEnd } = value;
  if (![xStartPx, xEndPx, yPx, bandTopPx, bandBottomPx].every(finite)) return null;
  if (!whole(displayStart) || !whole(displayEnd) || displayStart < 0 || displayEnd <= displayStart) return null;
  return /** @type {LinearAxisRow} */ ({ xStartPx, xEndPx, yPx, bandTopPx, bandBottomPx, displayStart, displayEnd });
};

/** @param {unknown} value @returns {RegionAxis | null} */
const readAxis = (value) => {
  if (!isObject(value)) return null;
  const source = readSource(value.source);
  if (!source) return null;
  if (value.kind === 'circular') {
    const { radiusPx, bandInnerPx, bandOuterPx, startAngleDeg, degreesPerBp } = value;
    if (![radiusPx, bandInnerPx, bandOuterPx, startAngleDeg, degreesPerBp].every(finite)) return null;
    return /** @type {CircularAxisGeometry} */ ({
      kind: 'circular', radiusPx, bandInnerPx, bandOuterPx, startAngleDeg, degreesPerBp, source
    });
  }
  if (value.kind !== 'linear' || !Array.isArray(value.rows) || value.rows.length === 0) return null;
  const rows = value.rows.map(readRow);
  if (rows.some((row) => row === null)) return null;
  return { kind: 'linear', rows: /** @type {LinearAxisRow[]} */ (rows), source };
};

/**
 * The records of one Result whose axis Python wrote; empty for a Result generated before it did (then Select region
 * cannot be used until the next Generate).
 * @param {unknown} trackSlotGeometry The run metadata `trackSlotGeometry` of the shown drawing.
 * @param {number} resultIndex
 * @returns {RegionAxisRecord[]}
 */
export const readRegionAxisRecords = (trackSlotGeometry, resultIndex) => {
  if (!isObject(trackSlotGeometry) || !Array.isArray(trackSlotGeometry.records)) return [];
  /** @type {RegionAxisRecord[]} */
  const records = [];
  for (const record of trackSlotGeometry.records) {
    if (!isObject(record) || record.resultIndex !== resultIndex) continue;
    const axis = readAxis(record.axis);
    const { recordIndex, recordKey, recordId, windowIndex } = record;
    if (!axis || !whole(recordIndex) || typeof recordKey !== 'string' || !recordKey) continue;
    // Python writes `windowIndex` only on the windows of a record drawn in several.
    if (windowIndex !== undefined && !(whole(windowIndex) && windowIndex >= 0)) continue;
    records.push({
      resultIndex, recordIndex, recordKey, windowIndex: windowIndex ?? null,
      recordId: typeof recordId === 'string' ? recordId : '', axis
    });
  }
  return records;
};

/**
 * The records of the shown Result whose axis Python wrote: the Result the preview shows of the shown drawing.
 * @param {{ trackSlotResolvedGeometry: { value: unknown }, selectedResultIndex: { value: unknown } }} shown
 * @returns {RegionAxisRecord[]}
 */
export const shownAxisRecords = ({ trackSlotResolvedGeometry, selectedResultIndex }) => (
  readRegionAxisRecords(trackSlotResolvedGeometry.value, Number(selectedResultIndex.value) || 0)
);

/**
 * The drawn group whose frame a record's axis is in: the innermost group of the record that is not one of its parts
 * (a definition, or a drawn annotation, OV-2482). Null when the Result does not draw the record.
 * @param {ParentNode} root The Result's SVG root.
 * @param {number} recordIndex
 * @returns {SVGGraphicsElement | null}
 */
export const regionAxisFrame = (root, recordIndex) => {
  const groups = Array.from(root.querySelectorAll(
    `g[data-gbdraw-record-index="${recordIndex}"]:not([data-gbdraw-role]):not(${ANNOTATION_SELECTOR})`
  ));
  const innermost = groups.filter((group) => !groups.some((other) => other !== group && group.contains(other)));
  return innermost.length === 1 ? /** @type {SVGGraphicsElement} */ (innermost[0]) : null;
};

/**
 * A screen point in a record frame's coordinates, and the frame's screen scale.
 * @param {SVGGraphicsElement} frame @param {number} clientX @param {number} clientY
 * @returns {{ point: FramePoint, scale: number } | null}
 */
export const framePoint = (frame, clientX, clientY) => {
  const matrix = frame.getScreenCTM();
  if (!matrix) return null;
  const point = new DOMPoint(clientX, clientY).matrixTransform(matrix.inverse());
  return { point: { x: point.x, y: point.y }, scale: Math.hypot(matrix.a, matrix.b) || 1 };
};

/** The SVG angle of a Circular display position (bases from display index 0). @param {CircularAxisGeometry} axis @param {number} position */
const circularAngleDeg = (axis, position) => axis.startAngleDeg + axis.degreesPerBp * position;

/**
 * Where a point lies along the axis, in bases from display index 0 (a fraction inside a base); `row` is the Linear
 * row. `inBand` says whether the point is on the record's feature tracks, widened by `slack` frame units.
 * @param {RegionAxis} axis
 * @param {FramePoint} point
 * @param {number} [slack]
 * @returns {{ position: number, row: number, inBand: boolean }}
 */
export const axisPosition = (axis, point, slack = 0) => {
  const { drawnLength } = axis.source;
  if (axis.kind === 'circular') {
    const radius = Math.hypot(point.x, point.y);
    const angle = Math.atan2(point.y, point.x) * 180 / Math.PI;
    const turn = 360 / Math.abs(axis.degreesPerBp);
    const position = mod((angle - axis.startAngleDeg) / axis.degreesPerBp, turn);
    return {
      position: Math.min(position, drawnLength - 1e-9),
      row: 0,
      inBand: radius >= axis.bandInnerPx - slack && radius <= axis.bandOuterPx + slack
    };
  }
  const nearest = axis.rows.reduce((best, row, index) => {
    const distance = Math.max(row.bandTopPx - point.y, point.y - row.bandBottomPx, 0);
    return distance < best.distance ? { index, distance } : best;
  }, { index: 0, distance: Infinity });
  return { ...linearRowPosition(axis, nearest.index, point.x), inBand: nearest.distance <= slack && linearOnRow(axis, nearest.index, point.x, slack) };
};

/** @param {LinearAxisGeometry} axis @param {number} row @param {number} x @param {number} slack */
const linearOnRow = (axis, row, x, slack) => {
  const { xStartPx, xEndPx } = axis.rows[row];
  return x >= Math.min(xStartPx, xEndPx) - slack && x <= Math.max(xStartPx, xEndPx) + slack;
};

/**
 * A Linear position on one row, clamped to the row's bases (a drag stays on the row it started on).
 * @param {LinearAxisGeometry} axis
 * @param {number} row
 * @param {number} x
 * @returns {{ position: number, row: number }}
 */
export const linearRowPosition = (axis, row, x) => {
  const { xStartPx, xEndPx, displayStart, displayEnd } = axis.rows[row];
  const fraction = xEndPx === xStartPx ? 0 : (x - xStartPx) / (xEndPx - xStartPx);
  const position = displayStart + fraction * (displayEnd - displayStart);
  return { position: Math.min(Math.max(position, displayStart), displayEnd - 1e-9), row };
};

/**
 * The span of a Circular drag that started at `start` and has turned `sweep` bases (negative counter-clockwise).
 * It crosses display index 0 only where Python says a span may cross the origin (`isCircular`); on any other record
 * it stops at the record's ends. A whole turn or more is the whole record.
 * @param {CircularAxisGeometry} axis
 * @param {number} start
 * @param {number} sweep
 * @returns {DisplaySpan}
 */
export const circularDragSpan = (axis, start, sweep) => {
  const { drawnLength } = axis.source;
  const from = Math.floor(start);
  if (Math.abs(sweep) >= drawnLength) return { start: 0, length: drawnLength };
  const to = Math.floor(start + sweep);
  if (!axis.source.isCircular) {
    const low = Math.max(0, Math.min(from, to));
    const high = Math.min(drawnLength - 1, Math.max(from, to));
    return { start: low, length: high - low + 1 };
  }
  return sweep >= 0
    ? { start: from, length: mod(to - from, drawnLength) + 1 }
    : { start: mod(to, drawnLength), length: mod(from - to, drawnLength) + 1 };
};

/**
 * The span between two positions on one Linear row.
 * @param {number} start
 * @param {number} end
 * @returns {DisplaySpan}
 */
export const linearDragSpan = (start, end) => {
  const low = Math.floor(Math.min(start, end));
  return { start: low, length: Math.floor(Math.max(start, end)) - low + 1 };
};

/** @param {AxisSourceMap} source @param {number} index */
const displaySourceBase = (source, index) => mod(source.anchorBase - 1 + source.step * index, source.sourceLength) + 1;

/**
 * The source base a display index shows.
 * @param {AxisSourceMap} source
 * @param {number} index
 */
export const sourceBaseAt = (source, index) => displaySourceBase(source, mod(index, source.sourceLength));

/**
 * A display span in the input file's coordinates: on the source strand, so a reverse-complemented record gives the
 * same span as the forward one.
 * @param {AxisSourceMap} source
 * @param {DisplaySpan} span
 * @returns {SourceSpan}
 */
export const sourceSpanOf = (source, span) => {
  if (span.length >= source.sourceLength) return { start: 1, end: source.sourceLength, length: source.sourceLength, wrapsOrigin: false };
  const first = sourceBaseAt(source, span.start);
  const last = sourceBaseAt(source, span.start + span.length - 1);
  const [start, end] = source.step === 1 ? [first, last] : [last, first];
  return { start, end, length: span.length, wrapsOrigin: start > end };
};

/**
 * The display index of the record's own base `offset` (0-based, as its feature coordinates count).
 * @param {AxisSourceMap} source
 * @param {number} offset
 */
export const displayIndexOfRecordOffset = (source, offset) => mod(
  offset + source.step * (source.recordBase - source.anchorBase),
  source.sourceLength
);

/**
 * The display index of the first drawn base of a feature part given in the input file's coordinates (OV-2718): the
 * part's first base, or its last on a reverse-complemented record, as the drawn record's own base.
 * @param {AxisSourceMap} source
 * @param {RecordPart} part
 */
const displayIndexOfPart = (source, { start, end }) => displayIndexOfRecordOffset(
  source, source.step === 1 ? mod(start - (source.recordBase - 1), source.sourceLength) : mod(source.recordBase - end, source.sourceLength)
);

/**
 * Whether every part of a feature lies inside a display span.
 * @param {AxisSourceMap} source
 * @param {DisplaySpan} span
 * @param {readonly RecordPart[]} parts
 */
export const partsInsideSpan = (source, span, parts) => parts.length > 0 && parts.every(({ start, end }) => {
  const length = end - start;
  if (!(length > 0)) return false;
  if (span.length >= source.sourceLength) return true;
  const offset = mod(displayIndexOfPart(source, { start, end }) - span.start, source.sourceLength);
  return offset + length <= span.length;
});

/**
 * The display position of the middle of a feature's parts (its first to its last part along the axis), in bases.
 * @param {AxisSourceMap} source
 * @param {readonly RecordPart[]} parts
 * @returns {number | null}
 */
export const partsMidpoint = (source, parts) => {
  if (parts.length === 0) return null;
  const starts = parts.map((part) => displayIndexOfPart(source, part));
  const first = starts[0];
  // Each part's extent along the axis measured from the first part's start.
  const reach = Math.max(...parts.map(({ start, end }, index) => mod(starts[index] - first, source.sourceLength) + (end - start)));
  return mod(first + reach / 2, source.sourceLength);
};

/**
 * The parts of a feature in the input file's coordinates.
 * @param {RegionFeatureRow} feature
 * @returns {RecordPart[]}
 */
export const featureParts = (feature) => {
  const listed = Array.isArray(feature.location_parts) && feature.location_parts.length ? feature.location_parts : [feature];
  return listed.flatMap((/** @type {{ start?: unknown, end?: unknown }} */ part) => {
    const start = Number(part?.start);
    const end = Number(part?.end);
    return Number.isSafeInteger(start) && Number.isSafeInteger(end) && end > start ? [{ start, end }] : [];
  });
};

/**
 * Whether display position `position` lies on the way from `from` to `to`: clockwise on a Circular record, as its
 * bases are drawn, and between the two on a Linear one.
 * @param {RegionAxis} axis
 * @param {number} from
 * @param {number} to
 * @param {number} position
 */
export const onTheWay = (axis, from, to, position) => (axis.kind === 'circular'
  ? mod(position - from, axis.source.sourceLength) <= mod(to - from, axis.source.sourceLength)
  : position >= Math.min(from, to) && position <= Math.max(from, to));

/** @param {number} value */
const fixed = (value) => Number(value.toFixed(3));

/**
 * The path of a selected span on the record's feature tracks, in the record's frame.
 * @param {RegionAxis} axis
 * @param {DisplaySpan} span
 * @returns {string}
 */
export const regionBandPath = (axis, span) => {
  if (axis.kind === 'circular') {
    const inner = axis.bandInnerPx;
    const outer = axis.bandOuterPx;
    if (span.length >= axis.source.drawnLength) {
      const ring = (/** @type {number} */ radius, /** @type {0 | 1} */ sweep) => (
        `M ${fixed(radius)} 0 A ${fixed(radius)} ${fixed(radius)} 0 1 ${sweep} ${fixed(-radius)} 0 `
        + `A ${fixed(radius)} ${fixed(radius)} 0 1 ${sweep} ${fixed(radius)} 0 Z`
      );
      return `${ring(outer, 1)} ${ring(inner, 0)}`;
    }
    const from = circularAngleDeg(axis, span.start) * Math.PI / 180;
    const to = circularAngleDeg(axis, span.start + span.length) * Math.PI / 180;
    const large = Math.abs(axis.degreesPerBp * span.length) > 180 ? 1 : 0;
    const sweep = axis.degreesPerBp > 0 ? 1 : 0;
    const point = (/** @type {number} */ radius, /** @type {number} */ angle) => `${fixed(radius * Math.cos(angle))} ${fixed(radius * Math.sin(angle))}`;
    return `M ${point(outer, from)} A ${fixed(outer)} ${fixed(outer)} 0 ${large} ${sweep} ${point(outer, to)} `
      + `L ${point(inner, to)} A ${fixed(inner)} ${fixed(inner)} 0 ${large} ${1 - sweep} ${point(inner, from)} Z`;
  }
  return axis.rows.flatMap((row) => {
    const start = Math.max(span.start, row.displayStart);
    const end = Math.min(span.start + span.length, row.displayEnd);
    if (end <= start) return [];
    const perBase = (row.xEndPx - row.xStartPx) / (row.displayEnd - row.displayStart);
    const x0 = row.xStartPx + perBase * (start - row.displayStart);
    const x1 = row.xStartPx + perBase * (end - row.displayStart);
    const left = Math.min(x0, x1);
    return [`M ${fixed(left)} ${fixed(row.bandTopPx)} h ${fixed(Math.abs(x1 - x0))} V ${fixed(row.bandBottomPx)} H ${fixed(left)} Z`];
  }).join(' ');
};

/**
 * The numbers of a source span as the UI shows them.
 * @param {SourceSpan} source
 */
export const regionSpanText = (source) => ({
  start: numberText(source.start), end: numberText(source.end), length: numberText(source.length)
});

/**
 * The Region popup's heading, from the selection alone: "Region", the record, and "start-end (N bp)".
 * @param {RegionSelection} selection
 */
export const regionSelectionHeading = (selection) => ({
  title: t('Region'),
  subtitle: selection.recordId || selection.recordKey,
  span: t('{start}-{end} ({length} bp)', regionSpanText(selection.source))
});
