// @ts-check
// The canvas slice of a drawing: its canvas padding (Session 46
// `modes.<m>.ui.canvasPadding`). A drawing shows the padding under its member
// name, read-only to the type checker; only this module writes it, through the
// slice's handle (`drawing.slices.canvas`). Load, a History restore, and Reset
// Layout write it with `install` and `reset`; Save and a History capture read
// it with `capture`.
import { reactive } from '../../utils/reactivity.js';
import { createSliceHandle, sliceMembers } from './slice-handle.js';
import { normalizeCanvasPadding } from '../canvas-padding.js';
/** @import { SliceHandle } from './slice-handle.js' */

/** @typedef {import('../canvas-padding.js').CanvasPadding} CanvasPadding */
/** @typedef {{ canvasPadding: CanvasPadding }} CanvasData What `capture` returns and `install` writes. */
/** @typedef {{ readonly canvasPadding: Readonly<CanvasPadding> }} CanvasView */
/** @typedef {SliceHandle<'canvas'>} CanvasSlice The handle a writer takes. */

const SIDES = /** @type {const} */ (['top', 'right', 'bottom', 'left']);

// A canvas handle exists only through `create`, which gives it these members.
/** @param {CanvasSlice} slice */
const membersOf = (slice) => /** @type {CanvasData} */ (sliceMembers(slice, 'canvas'));

/**
 * @param {CanvasSlice} slice
 * @param {(side: typeof SIDES[number]) => number} sideValue
 */
const writeSides = (slice, sideValue) => {
  const { canvasPadding } = membersOf(slice);
  SIDES.forEach((side) => { canvasPadding[side] = sideValue(side); });
};

/**
 * A new slice, or one over a drawing's existing members (a test fixture).
 * @param {Partial<CanvasData>} [members]
 * @returns {{ handle: CanvasSlice, view: CanvasView }}
 */
const create = (members = {}) => {
  const data = { canvasPadding: members.canvasPadding ?? reactive({ top: 0, right: 0, bottom: 0, left: 0 }) };
  return { handle: createSliceHandle('canvas', data), view: Object.freeze({ canvasPadding: data.canvasPadding }) };
};

/**
 * @param {CanvasSlice} slice
 * @returns {CanvasData}
 */
const capture = (slice) => ({ canvasPadding: { ...membersOf(slice).canvasPadding } });

// Copies the padding `data` holds into the slice's padding as the canvas draws
// it (`normalizeCanvasPadding`: a negative or non-finite side is 0, review S2,
// D-59); without one the slice keeps its padding (a History step or a Session
// that holds none).
/**
 * @param {CanvasSlice} slice
 * @param {Partial<CanvasData>} data
 */
const install = (slice, { canvasPadding }) => {
  if (!canvasPadding) return;
  const drawable = normalizeCanvasPadding(canvasPadding);
  writeSides(slice, (side) => drawable[side]);
};

/** @param {CanvasSlice} slice */
const reset = (slice) => writeSides(slice, () => 0);

/**
 * One side of the padding, as its number field gives it and the canvas draws
 * it: a negative side is 0 (review S2, D-59), so the slice holds only values
 * the request sends as they are (R7). A field that holds no number (empty, or
 * a value being typed) leaves the side.
 * @param {CanvasSlice} slice
 * @param {typeof SIDES[number]} side
 * @param {unknown} value
 */
export const setCanvasPaddingSide = (slice, side, value) => {
  const number = typeof value === 'number' ? value : Number.parseFloat(String(value));
  if (Number.isFinite(number)) membersOf(slice).canvasPadding[side] = normalizeCanvasPadding({ [side]: number })[side];
};

export const CANVAS_SLICE = Object.freeze({ members: Object.freeze(['canvasPadding']), create, capture, install, reset });
