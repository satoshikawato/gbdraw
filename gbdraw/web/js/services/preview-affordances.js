// @ts-check
// Marks that the preview adds to a mounted Result and that saving and export
// remove (services/svg-serialization.js stripTransientPreviewState), and the
// one writer of the Preview's classes there.

// A match band that the preview made a button (role, tabindex and the
// translated aria-label; app/feature-editor/svg-actions.js
// preparePairwiseInteractionAffordances). Saving finds it by this mark, not by
// its label text.
export const PREVIEW_MATCH_BUTTON_ATTRIBUTE = 'data-gbdraw-preview-match';

// The marks the label binder (app/feature-editor/label-actions.js bindLabels)
// writes on each label of a mounted Result it binds for editing: that it is
// editable, its key in the Features list, and the text it showed when bound.
// Python writes none of them. An export strips them, so its bytes do not
// depend on whether a label was clicked and a renamed label's old text does
// not leave the page (OV-3848); the Result keeps them for the editor.
export const LABEL_EDITABLE_ATTRIBUTE = 'data-label-editable';
export const LABEL_KEY_ATTRIBUTE = 'data-label-key';
export const LABEL_SOURCE_TEXT_ATTRIBUTE = 'data-label-source-text';
export const LABEL_BINDER_ATTRIBUTES = Object.freeze([
  LABEL_EDITABLE_ATTRIBUTE, LABEL_KEY_ATTRIBUTE, LABEL_SOURCE_TEXT_ATTRIBUTE
]);

/** @param {Element} root An export's copy of the Result. */
export const stripLabelBinderAttributes = (root) => {
  root.querySelectorAll(LABEL_BINDER_ATTRIBUTES.map((name) => `[${name}]`).join(', ')).forEach((element) => {
    LABEL_BINDER_ATTRIBUTES.forEach((name) => element.removeAttribute(name));
  });
};

/**
 * THE writer of a Preview class on an element of the Result: it adds or removes `token`, and removes the class
 * attribute when no class is left, since Python never draws an empty one and Export and Save read the Result
 * (OV-3241: a `classList.toggle` that removed the last class left `class=""` in the SVG download).
 * @param {Element | null | undefined} element
 * @param {string} token
 * @param {unknown} enabled
 */
export const setClassToken = (element, token, enabled) => {
  if (!element) return;
  if (element.classList?.toggle) {
    element.classList.toggle(token, Boolean(enabled));
    if (element.classList.length === 0) element.removeAttribute('class');
    return;
  }
  const existing = String(element.getAttribute('class') || '').trim();
  const tokens = existing ? existing.split(/\s+/).filter(Boolean) : [];
  const nextTokens = tokens.filter((entry) => entry !== token);
  if (enabled) nextTokens.push(token);
  if (nextTokens.length) {
    element.setAttribute('class', nextTokens.join(' '));
  } else {
    element.removeAttribute('class');
  }
};
