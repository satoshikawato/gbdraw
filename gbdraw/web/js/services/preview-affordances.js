// @ts-check
// Marks that the preview adds to a mounted Result and that saving and export
// remove (services/svg-serialization.js stripTransientPreviewState), and the
// one writer of the Preview's classes there.

// A match band that the preview made a button (role, tabindex and the
// translated aria-label; app/feature-editor/svg-actions.js
// preparePairwiseInteractionAffordances). Saving finds it by this mark, not by
// its label text.
export const PREVIEW_MATCH_BUTTON_ATTRIBUTE = 'data-gbdraw-preview-match';

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
