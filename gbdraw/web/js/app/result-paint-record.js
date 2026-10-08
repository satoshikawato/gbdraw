// @ts-check
// Which editor state each Result's bytes show (design §5.4). A Result
// displayed again receives a domain only when that state changed since it
// was last shown, so a Result without new edits gets no projection work and
// an Undo reaches a Result that is displayed again. The displayed Result
// follows every live edit; app-setup reads and records the states.

/**
 * @typedef {{
 *   colors: [unknown, string],
 *   visibility: string,
 *   strokes: string,
 *   labels: string,
 *   legendOrder: string
 * }} EditorPaintState
 */

// A paint state no editor state equals: a Result recorded with it shows
// every paint domain on its next display.
/** @param {EditorPaintState} current @returns {EditorPaintState} */
const unknownPaint = (current) => ({ ...current, colors: [null, ''], visibility: '', strokes: '' });

export const createResultPaintRecord = () => {
  /** @type {Map<string, EditorPaintState>} */
  const shownByResult = new Map();
  let bound = '';
  // The paint domains the displayed Result was not shown (a declined
  // display); recorded with it when it leaves.
  /** @type {Partial<EditorPaintState> | null} */
  let behind = null;
  /** @param {EditorPaintState} state */
  const recordBound = (state) => {
    if (bound && shownByResult.has(bound)) shownByResult.set(bound, { ...state, ...behind });
    behind = null;
  };
  return {
    /**
     * A committed artifact: each Result first seen shows `current`, and a
     * Result no longer live is forgotten. `displayed` is the Result mounted
     * by the commit. A `restored` artifact (Session Load, a History restore)
     * keeps each Result's bytes as they were saved or kept: a batch Result
     * not displayed since the last edits lacks them, so every Result but the
     * mounted one shows every paint domain on its first display.
     * @param {readonly string[]} identities
     * @param {string} displayed
     * @param {EditorPaintState} current
     * @param {{ restored?: boolean }} [options]
     */
    commit(identities, displayed, current, { restored = false } = {}) {
      const live = new Set(identities.filter(Boolean));
      live.forEach((identity) => {
        if (shownByResult.has(identity)) return;
        shownByResult.set(identity, restored && identity !== displayed ? unknownPaint(current) : current);
      });
      [...shownByResult.keys()].forEach((identity) => {
        if (!live.has(identity)) shownByResult.delete(identity);
      });
      bound = displayed;
      behind = null;
    },
    /**
     * The displayed Result leaves (a mode switch) having followed every live
     * edit. Returns its identity, or '' when none is recorded.
     * @param {EditorPaintState} current
     */
    depart(current) {
      const departed = bound && shownByResult.has(bound) ? bound : '';
      recordBound(current);
      bound = '';
      return departed;
    },
    /**
     * A Result is displayed: the Result shown until now followed every live
     * edit. Returns the state the arriving Result's bytes show.
     * @param {string} identity
     * @param {EditorPaintState} current
     * @returns {EditorPaintState}
     */
    display(identity, current) {
      if (bound !== identity) recordBound(current);
      bound = identity;
      return shownByResult.get(identity) || current;
    },
    /**
     * The displayed Result now shows `current`, or, when its palette and
     * rules projection `declined`, the fills and visibility of `previous`.
     * @param {string} identity
     * @param {EditorPaintState} current
     * @param {EditorPaintState} previous
     * @param {{ declined?: boolean }} [options]
     */
    shown(identity, current, previous, { declined = false } = {}) {
      behind = declined ? { colors: previous.colors, visibility: previous.visibility } : null;
      shownByResult.set(identity, { ...current, ...behind });
    },
    /** @param {string} identity The Result shown again with the intent it left with (E1). */
    rebind(identity) {
      bound = identity;
      behind = null;
    }
  };
};
