// @ts-check

/** @typedef {number | string | null} FontSizeInput A font size input: a number, the raw text when it is not numeric, or null when empty. */
/** @typedef {{ status: string, reason: string }} SessionBusyOutcome */

/** @param {unknown} value */
const optionalNumberInputValue = (value) => {
  const text = String(value ?? '').trim();
  if (text === '') return null;
  const numeric = Number(text);
  return Number.isNaN(numeric) ? text : numeric;
};

/**
 * @typedef {Object} LinearTypographyControllerOptions
 * @property {{ scale_font_size: FontSizeInput, ruler_label_font_size: FontSizeInput }} adv The two font sizes of the advanced-options state owned by `state.js`.
 * @property {{ value: boolean }} linked The ref that links the two font sizes.
 * @property {() => SessionBusyOutcome | null} [mutationAvailability] Returns a busy outcome while a Session operation blocks edits.
 */

/** @param {LinearTypographyControllerOptions} options */
export const createLinearTypographyController = ({ adv, linked, mutationAvailability = () => null }) => {
  /** @param {unknown} value */
  const setScaleFontSize = (value) => {
    const sessionBusy = mutationAvailability();
    if (sessionBusy) return sessionBusy;
    const nextValue = optionalNumberInputValue(value);
    adv.scale_font_size = nextValue;
    if (linked.value) adv.ruler_label_font_size = nextValue;
    return nextValue;
  };

  /** @param {unknown} value */
  const setRulerLabelFontSize = (value) => {
    const sessionBusy = mutationAvailability();
    if (sessionBusy) return sessionBusy;
    if (linked.value) return false;
    adv.ruler_label_font_size = optionalNumberInputValue(value);
    return true;
  };

  /** @param {unknown} nextLinked */
  const setLinked = (nextLinked) => {
    const sessionBusy = mutationAvailability();
    if (sessionBusy) return sessionBusy;
    const nextValue = Boolean(nextLinked);
    if (linked.value === nextValue) return false;
    linked.value = nextValue;
    if (nextValue) adv.ruler_label_font_size = adv.scale_font_size;
    return true;
  };

  return {
    setLinked,
    setRulerLabelFontSize,
    setScaleFontSize
  };
};
