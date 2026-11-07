// @ts-check
// Keeps Tab inside the modal dialog of the overlay that handles the keydown.
/** @param {KeyboardEvent & { currentTarget: Element }} event Tab on the overlay that holds the dialog. */
export const trapDialogFocus = (event) => {
  const dialog = /** @type {HTMLElement | null} */ (event.currentTarget?.querySelector('[role="dialog"]'));
  const controls = /** @type {HTMLElement[]} */ (Array.from(dialog?.querySelectorAll('button:not(:disabled)') || []));
  if (!controls.length) {
    // A busy choice disables every button: Tab stays on the dialog panel.
    if (dialog) {
      event.preventDefault();
      dialog.focus({ preventScroll: true });
    }
    return;
  }
  // `controls` is not empty, so `dialog` is an element and `last` is a button.
  const first = controls[0];
  const last = /** @type {HTMLElement} */ (controls.at(-1));
  if (event.shiftKey && (document.activeElement === first || !(/** @type {HTMLElement} */ (dialog)).contains(document.activeElement))) {
    event.preventDefault();
    last.focus();
  } else if (!event.shiftKey && (document.activeElement === last || !(/** @type {HTMLElement} */ (dialog)).contains(document.activeElement))) {
    event.preventDefault();
    first.focus();
  }
};
