// @ts-check
// The Canvas padding a drawing draws, px per side: the one normalizer of it, a
// leaf module so the canvas slice (services/drawing/canvas.js), the canvas fit,
// and decoration continuity share it.

/** @typedef {{ top: number, right: number, bottom: number, left: number }} CanvasPadding */

/**
 * Whether a padding side is one the canvas draws: a finite number of 0 or more. Load reads a stored side by it too
 * (`readCanvasPadding`), and lists a side it is not.
 * @param {unknown} number
 * @returns {number is number}
 */
export const isDrawablePaddingSide = (number) => typeof number === 'number' && Number.isFinite(number) && number >= 0;

/**
 * The padding a drawing draws, the one normalizer of it: each side as a
 * number, a side that holds no finite number of 0 or more being 0. The
 * slice's writers store the padding through it (review S2, D-59), so the
 * request projection sends the slice as it is (R7); the canvas fit reads a
 * padding it may not be given (none, or a partial one) through it.
 * @param {Readonly<Partial<Record<string, unknown>>> | null | undefined} value
 * @returns {CanvasPadding}
 */
export const normalizeCanvasPadding = (value) => {
  /** @param {'top' | 'right' | 'bottom' | 'left'} side */
  const sideOf = (side) => {
    const number = Number(value?.[side]);
    return isDrawablePaddingSide(number) ? number : 0;
  };
  return { top: sideOf('top'), right: sideOf('right'), bottom: sideOf('bottom'), left: sideOf('left') };
};

/** @param {Readonly<Partial<Record<string, unknown>>> | null | undefined} value */
export const hasCanvasPadding = (value) => Object.values(normalizeCanvasPadding(value)).some((number) => number !== 0);
