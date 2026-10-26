// @ts-check
// Which drawings a Session 49 Save writes and which one its Load opens
// (PD-OI-090 〔Save と Load〕, R16-8; DESIGN §4.4). State-free: the caller says
// which drawing is shown and which drawings are used.

/**
 * What Save writes: the used drawings in tab order, and the drawing Load opens.
 * @template {{ readonly id: string }} D
 * @typedef {{ drawings: D[], activeDrawingId: D['id'] }} SavedDrawings
 */

/**
 * A used drawing is written; the shown one opens when it is used. An unused
 * shown drawing is not written: the next written tab opens, else the
 * previous. When no drawing is used, the shown one is written (a Session holds
 * at least one drawing).
 * @template {{ readonly id: string }} D
 * @param {readonly D[]} drawings The drawings in tab order.
 * @param {D} shown The drawing shown at Save, found in `drawings` by its ID.
 * @param {(drawing: D) => boolean} isUsed
 * @returns {SavedDrawings<D>}
 */
export const savedDrawingsOf = (drawings, shown, isUsed) => {
  const used = drawings.map((drawing) => isUsed(drawing));
  const saved = drawings.filter((_drawing, index) => used[index]);
  if (!saved.length) return { drawings: [shown], activeDrawingId: shown.id };
  const shownIndex = drawings.findIndex((drawing) => drawing.id === shown.id);
  if (used[shownIndex]) return { drawings: saved, activeDrawingId: shown.id };
  const next = drawings.find((_drawing, index) => index > shownIndex && used[index]);
  const previous = drawings.findLast((_drawing, index) => index < shownIndex && used[index]);
  return { drawings: saved, activeDrawingId: /** @type {D} */ (next ?? previous).id };
};
