// @ts-check
// R10, Q3, R13: a track stack editor routes each of its edits of a feature-slot
// input through `changeTrackLayout`, the feature placement owner's transition
// that the composition root injects as a port. A DOM event the template passes
// last names the control for Cancel and focus (its listener's element, not an
// icon inside a button); the edit itself runs without the event.

/**
 * The edited checkbox, select, or button, with the members the placement owner reads.
 * @typedef {HTMLElement & {
 *   type?: string,
 *   checked?: boolean,
 *   selectedOptions?: ArrayLike<{ text: string }>,
 *   labels?: ArrayLike<{ textContent: string }> | null
 * }} LayoutControl
 */

/**
 * R10: the feature placement owner's transition for an edit of a draft
 * feature-slot input. It runs `apply` (and may ask before keeping it) and
 * returns what `apply` returns, a busy outcome, or `false` while the user is asked.
 * @typedef {(apply: () => any, control?: LayoutControl | null) => any} ChangeTrackLayout
 */

/**
 * The Legend port a stack edit runs in: a drawn track the edit hides or shows
 * again takes its Legend rows with it (R2, OV-1198).
 * @typedef {<T>(change: () => T) => T} ChangeTrackVisibility
 */

/**
 * The port a track stack editor passes to `featureSlotEdits`: each stack edit
 * (add, renderer change, enable, disable, remove, stack switch) runs in the
 * Legend port once. The port wraps the whole layout transition, so an edit
 * whose lane check restores the inputs and asks retires nothing; the edit that
 * the Reset choice of that dialog applies later runs in the port on its own,
 * inside the Reset's History step.
 * @param {ChangeTrackLayout} changeTrackLayout
 * @param {ChangeTrackVisibility} changeTrackVisibility
 * @returns {ChangeTrackLayout}
 */
export const trackLayoutWithVisibility = (changeTrackLayout, changeTrackVisibility) => (apply, control) => {
  let inTransition = true;
  try {
    return changeTrackVisibility(() => changeTrackLayout(
      () => (inTransition ? apply() : changeTrackVisibility(apply)), control
    ));
  } finally {
    inTransition = false;
  }
};

/**
 * @template {Record<string, (...args: any[]) => any>} Actions
 * @param {ChangeTrackLayout} changeTrackLayout
 * @param {Actions} actions
 * @returns {Actions} The same names; each takes the optional trailing DOM event as well.
 */
export const featureSlotEdits = (changeTrackLayout, actions) => /** @type {Actions} */ (Object.fromEntries(
  Object.entries(actions).map(([name, action]) => [name, (...args) => {
    const event = typeof Event === 'function' && args.at(-1) instanceof Event ? args.pop() : null;
    return changeTrackLayout(() => action(...args), event?.currentTarget);
  }])
));
