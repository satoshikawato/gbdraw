// @ts-check
// Marks that the preview adds to a mounted Result and that saving and export
// remove (services/svg-serialization.js stripTransientPreviewState), the
// label source-text record, and the one writer of the Preview's classes there.
import { FEATURE_LABEL_SELECTOR, LABEL_FEATURE_ID_ATTRIBUTE } from './drawn-groups.js';

// A match band that the preview made a button (role, tabindex and the
// translated aria-label; app/feature-editor/svg-actions.js
// preparePairwiseInteractionAffordances). Saving finds it by this mark, not by
// its label text.
export const PREVIEW_MATCH_BUTTON_ATTRIBUTE = 'data-gbdraw-preview-match';

// The marks the label binder (app/feature-editor/label-actions.js bindLabels)
// writes on each label of a mounted Result it binds for editing: that it is
// editable and its key in the Features list. Python writes neither, a
// Result's text (History, Sessions) holds neither, and Load drops them
// (services/svg-sanitization.js).
// The source-text record: a label's source text (the text it shows without
// label edits) where that differs from the text it shows. Like a paint base
// (result-paint-bases.js), it travels with the Result through History and
// Sessions, so an Undo of an edit returns to it while no feature row holds
// it; an export strips it. Its value has one rule, `labelSourceTextOf`, and
// the Result's text holds it independent of clicks: a Result entering the
// page (Generate, a rerender, Load) records it from the feature rows
// (`recordLabelSourceTexts`), the label projection writes it on each bound
// label of the mounted Result by the same rule, and serialization drops one
// equal to the shown text (`stripLabelBinderMarks`). So neither a
// saved Session nor an export depends on whether a label was clicked, and a
// renamed label's old text does not leave the page (OV-3848).
export const LABEL_EDITABLE_ATTRIBUTE = 'data-label-editable';
export const LABEL_KEY_ATTRIBUTE = 'data-label-key';
export const LABEL_SOURCE_TEXT_ATTRIBUTE = 'data-label-source-text';
// The labels the binder made editable.
export const EDITABLE_LABEL_SELECTOR = `text[${LABEL_EDITABLE_ATTRIBUTE}="true"]`;

/** @param {Element} textEl A label's text element. @returns {string} the text it shows */
export const labelTextOf = (textEl) => {
  const textPath = textEl.querySelector('textPath');
  if (textPath) return textPath.textContent || '';
  return textEl.textContent || '';
};

/**
 * A label's source text, by the one rule the label actions, the projection
 * and the admission share: its feature row's recorded source text, else the label's
 * record, else the text it shows.
 * @param {string | null | undefined} rowSourceText
 * @param {Element} textEl
 * @returns {string}
 */
export const labelSourceTextOf = (rowSourceText, textEl) => rowSourceText
  ?? textEl.getAttribute(LABEL_SOURCE_TEXT_ATTRIBUTE) ?? labelTextOf(textEl);

/**
 * Drops the binder's marks and every source-text record equal to the shown text.
 * @param {Element} root A copy of the mounted Result that saving or an export writes, or an admitted Result.
 * @returns {boolean} whether it changed the SVG
 */
export const stripLabelBinderMarks = (root) => {
  let changed = false;
  root.querySelectorAll(`[${LABEL_EDITABLE_ATTRIBUTE}], [${LABEL_KEY_ATTRIBUTE}]`).forEach((element) => {
    element.removeAttribute(LABEL_EDITABLE_ATTRIBUTE);
    element.removeAttribute(LABEL_KEY_ATTRIBUTE);
    changed = true;
  });
  root.querySelectorAll(`[${LABEL_SOURCE_TEXT_ATTRIBUTE}]`).forEach((element) => {
    if (element.getAttribute(LABEL_SOURCE_TEXT_ATTRIBUTE) === labelTextOf(element)) {
      element.removeAttribute(LABEL_SOURCE_TEXT_ATTRIBUTE);
      changed = true;
    }
  });
  return changed;
};

/**
 * The fields of a feature row (a draft row of the drawing or a Session) that
 * the source-text records read.
 * @typedef {object} LabelSourceRow
 * @property {string} recordKey
 * @property {string} biologicalFeatureId
 * @property {string | null} [labelSourceText]
 */

/**
 * The source text each feature row records, by Result and by the rendered
 * feature ID its labels name (`data-label-feature-id`); `targetsOf` gives a
 * row's rendered features (services/feature-catalog.js featureRowTargets).
 * @param {Readonly<Record<string, LabelSourceRow>> | null | undefined} rows
 * @param {(row: LabelSourceRow) => readonly { resultIndex: number, renderedId: string }[]} targetsOf
 * @returns {Map<number, Map<string, string>>}
 */
export const labelSourceTextsByResult = (rows, targetsOf) => {
  /** @type {Map<number, Map<string, string>>} */
  const byResult = new Map();
  Object.values(rows || {}).forEach((row) => {
    const sourceText = row.labelSourceText;
    if (typeof sourceText !== 'string' || !sourceText) return;
    targetsOf(row).forEach(({ resultIndex, renderedId }) => {
      const sources = byResult.get(resultIndex) || new Map();
      sources.set(renderedId, sourceText);
      byResult.set(resultIndex, sources);
    });
  });
  return byResult;
};

/**
 * THE writer of a Result's source-text records when it enters the page
 * (Generate, a rerender, Load): each label of a feature whose row records a
 * source text holds it by `labelSourceTextOf`, then the binder's marks and the
 * records equal to the shown text leave.
 * @param {Element} root The admitted Result.
 * @param {ReadonlyMap<string, string>} rowSources The rows' source texts by rendered feature ID.
 * @returns {boolean} whether it changed the SVG
 */
export const recordLabelSourceTexts = (root, rowSources) => {
  let changed = false;
  if (rowSources.size > 0) {
    root.querySelectorAll(FEATURE_LABEL_SELECTOR).forEach((textEl) => {
      const rowSource = rowSources.get(textEl.getAttribute(LABEL_FEATURE_ID_ATTRIBUTE) || '');
      if (rowSource === undefined) return;
      const sourceText = labelSourceTextOf(rowSource, textEl);
      if (sourceText === labelTextOf(textEl) || textEl.getAttribute(LABEL_SOURCE_TEXT_ATTRIBUTE) === sourceText) return;
      textEl.setAttribute(LABEL_SOURCE_TEXT_ATTRIBUTE, sourceText);
      changed = true;
    });
  }
  return stripLabelBinderMarks(root) || changed;
};

/** @param {Element} root An export's copy of the Result. */
export const stripLabelSourceRecords = (root) => {
  root.querySelectorAll(`[${LABEL_SOURCE_TEXT_ATTRIBUTE}]`).forEach((element) => {
    element.removeAttribute(LABEL_SOURCE_TEXT_ATTRIBUTE);
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
