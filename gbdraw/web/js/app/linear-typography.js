const optionalNumberInputValue = (value) => {
  const text = String(value ?? '').trim();
  if (text === '') return null;
  const numeric = Number(text);
  return Number.isNaN(numeric) ? text : numeric;
};

const linearTypographyValuesMatch = (adv = {}) => (
  Object.is(adv.scale_font_size, adv.ruler_label_font_size)
);

export const reconcileImportedLinearTypographyLink = ({ adv, linked, ui = {} }) => {
  if (!linked || typeof linked !== 'object' || !('value' in linked)) return false;
  // Omission takes the fresh linked default; unequal values still open unlinked.
  linked.value = (
    (ui.linearTypographyLinked ?? true) === true
    && linearTypographyValuesMatch(adv)
  );
  return linked.value;
};

export const createLinearTypographyController = ({ adv, linked, mutationAvailability = () => null }) => {
  const setScaleFontSize = (value) => {
    const sessionBusy = mutationAvailability();
    if (sessionBusy) return sessionBusy;
    const nextValue = optionalNumberInputValue(value);
    adv.scale_font_size = nextValue;
    if (linked.value) adv.ruler_label_font_size = nextValue;
    return nextValue;
  };

  const setRulerLabelFontSize = (value) => {
    const sessionBusy = mutationAvailability();
    if (sessionBusy) return sessionBusy;
    if (linked.value) return false;
    adv.ruler_label_font_size = optionalNumberInputValue(value);
    return true;
  };

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
