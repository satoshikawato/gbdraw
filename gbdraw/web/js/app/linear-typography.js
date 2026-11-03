// @ts-check
/** @import { SessionBusy } from '../state.js' */

/** @typedef {number | string | null} FontSizeInput A font size input: a number, the raw text when it is not numeric, or null when empty. */

/** @param {unknown} value */
const optionalNumberInputValue = (value) => {
  const text = String(value ?? '').trim();
  if (text === '') return null;
  const numeric = Number(text);
  return Number.isNaN(numeric) ? text : numeric;
};

/**
 * @typedef {Object} LinearTypographyDrawing
 * @property {{ scale_font_size: FontSizeInput, ruler_label_font_size: FontSizeInput }} adv The two font sizes of the advanced-options state owned by `state.js`.
 * @property {{ value: boolean }} linearTypographyLinked The ref that links the two font sizes.
 */

/**
 * @typedef {Object} LinearTypographyControllerOptions
 * @property {() => LinearTypographyDrawing | null} shownDrawing The shown drawing when it is the Linear one, else null (the handlers then do nothing); resolved at each call.
 * @property {() => SessionBusy | null} [mutationAvailability] Returns a busy outcome while a Session operation blocks edits.
 */

/** @param {LinearTypographyControllerOptions} options */
export const createLinearTypographyController = ({ shownDrawing, mutationAvailability = () => null }) => {
  /** @param {unknown} value */
  const setScaleFontSize = (value) => {
    const sessionBusy = mutationAvailability();
    if (sessionBusy) return sessionBusy;
    const drawing = shownDrawing();
    if (!drawing) return null;
    const nextValue = optionalNumberInputValue(value);
    drawing.adv.scale_font_size = nextValue;
    if (drawing.linearTypographyLinked.value) drawing.adv.ruler_label_font_size = nextValue;
    return nextValue;
  };

  /** @param {unknown} value */
  const setRulerLabelFontSize = (value) => {
    const sessionBusy = mutationAvailability();
    if (sessionBusy) return sessionBusy;
    const drawing = shownDrawing();
    if (!drawing || drawing.linearTypographyLinked.value) return false;
    drawing.adv.ruler_label_font_size = optionalNumberInputValue(value);
    return true;
  };

  /** @param {unknown} nextLinked */
  const setLinked = (nextLinked) => {
    const sessionBusy = mutationAvailability();
    if (sessionBusy) return sessionBusy;
    const drawing = shownDrawing();
    if (!drawing) return false;
    const nextValue = Boolean(nextLinked);
    if (drawing.linearTypographyLinked.value === nextValue) return false;
    drawing.linearTypographyLinked.value = nextValue;
    if (nextValue) drawing.adv.ruler_label_font_size = drawing.adv.scale_font_size;
    return true;
  };

  return {
    setLinked,
    setRulerLabelFontSize,
    setScaleFontSize
  };
};
