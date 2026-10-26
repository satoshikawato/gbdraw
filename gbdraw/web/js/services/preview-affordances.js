// @ts-check
// Marks that the preview adds to a mounted Result and that saving and export
// remove (services/svg-serialization.js stripTransientPreviewState).

// A match band that the preview made a button (role, tabindex and the
// translated aria-label; app/feature-editor/svg-actions.js
// preparePairwiseInteractionAffordances). Saving finds it by this mark, not by
// its label text.
export const PREVIEW_MATCH_BUTTON_ATTRIBUTE = 'data-gbdraw-preview-match';
