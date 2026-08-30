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
  linked.value = (
    ui.linearTypographyLinked === true
    && linearTypographyValuesMatch(adv)
  );
  return linked.value;
};

export const createLinearTypographyController = ({ adv, linked }) => {
  const setScaleFontSize = (value) => {
    const nextValue = optionalNumberInputValue(value);
    adv.scale_font_size = nextValue;
    if (linked.value) adv.ruler_label_font_size = nextValue;
    return nextValue;
  };

  const setRulerLabelFontSize = (value) => {
    if (linked.value) return false;
    adv.ruler_label_font_size = optionalNumberInputValue(value);
    return true;
  };

  const setLinked = (nextLinked) => {
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
