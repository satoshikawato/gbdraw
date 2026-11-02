// The timing notes `SettingNote` renders (gbdraw/web/js/components.js,
// PD-OI-037 rev 3): a change that applies on Generate shows a clock icon, a
// help tip named "Applies on Generate" whose description holds the note; a
// change that applies at once shows no mark.
const CLOCK = '.apply-note .help-tip-trigger';

/** The clocks inside `scope`, which may be a note itself. @param {import('@playwright/test').Locator} scope */
const clocksIn = (scope) => scope.locator(`:scope${CLOCK}, ${CLOCK}`);

/**
 * The notes of the clocks inside `scope`, in document order: the text of each
 * clock's description, which the page keeps in a hidden element of <body>.
 * @param {import('@playwright/test').Locator} scope
 * @returns {Promise<string[]>}
 */
const clockNotes = (scope) => clocksIn(scope).evaluateAll((buttons) => buttons.map((button) => {
  const description = document.getElementById(button.getAttribute('aria-describedby') || '');
  if (!description) throw new Error(`A clock without its description: ${button.outerHTML}`);
  return description.textContent.trim();
}));

module.exports = { CLOCK, clocksIn, clockNotes };
