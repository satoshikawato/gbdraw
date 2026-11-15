// The palette menu as a user works it (Colors > Palette), for the specs that
// change the palette with Instant Preview on: the Instant Preview checkbox,
// when it is off, and then a choice in the Palette select. The page imports
// this module (`/tests/web/helpers/palette-menu-page.mjs`).
//
// Each is its own History step, as in the app: turning Instant Preview on is
// one project step (OV-1160), and a choice in the select is one step of the
// shown drawing, which the input adapter records from the select's change
// event (OV-1071). A palette set while the Instant Preview step still captures
// its drawings would become part of that project step, which an Undo in any
// drawing takes, so an Undo in the other drawing would take the palette back
// (OV-3605).

/**
 * The members of the app (`window.__GBDRAW_APP__`) this module reads.
 * @typedef {{ paletteInstantPreviewEnabled: boolean, setPaletteInstantPreview: (enabled: boolean) => unknown }} PaletteMenuApp
 */

/**
 * Turns Instant Preview on when it is off, then chooses `name` in the Palette
 * select. A palette dialog the choice opens stays open for the caller to
 * answer.
 * @param {PaletteMenuApp} app
 * @param {string} name
 */
export const choosePaletteWithInstantPreview = async (app, name) => {
  if (!app.paletteInstantPreviewEnabled) await app.setPaletteInstantPreview(true);
  const select = document.querySelector('select[aria-describedby="palette-application-timing"]');
  if (!(select instanceof HTMLSelectElement)) throw new Error('no Palette select');
  if (select.disabled) throw new Error('the Palette select is disabled');
  if (![...select.options].some((option) => option.value === name)) throw new Error(`no palette ${name} in the Palette select`);
  select.value = name;
  select.dispatchEvent(new Event('change', { bubbles: true }));
};
