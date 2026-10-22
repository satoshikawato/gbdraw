// @ts-check
// Parse HTML or SVG markup without the page document (OV-1030, OV-1032).
// Parsing into an element of the page document costs, for each parsed node, an
// update of every live `labels` list of the page (Chromium keeps them once a
// script, such as Playwright's `getByLabel`, has read them): seconds for a
// large Result instead of tens of milliseconds. A <template>'s content belongs
// to an inert document, so the markup is parsed there and only the finished
// tree reaches the page, in one insertion (the preview) or one adoption (a
// node appended to the page).

/**
 * The parsed nodes of `markup`, owned by an inert document. Node methods that
 * need the page (layout, style, measurement) work on them only after they are
 * inserted into the page.
 * @param {string | null | undefined} markup
 * @returns {DocumentFragment}
 */
export const parseMarkupInertly = (markup) => {
  const template = document.createElement('template');
  template.innerHTML = String(markup ?? '');
  return template.content;
};
