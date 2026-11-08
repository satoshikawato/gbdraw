// @ts-check
/** @import { RegionSelection } from './region-axis.js' */
/** @import { MatchPopupPayload } from './pairwise-match-popup.js' */
// What the match popup frame (`state.clickedPairwiseMatch`) shows, and how a
// click in the preview finds the Legend row or record definition it opens.
// The frame shows a match (`buildMatchPopupPayload`) or a preview item
// (`app/preview-item-popup.js` derives that popup's view).
import { COMPARISON_LEGEND_SELECTOR, pythonLegendKey } from './legend-svg.js';

/**
 * The item a preview popup is about: a Legend row by its key, or the
 * definition of the Result's record at `recordIndex`.
 * @typedef {{ kind: 'legend-entry', key: string }
 *   | { kind: 'record-definition', recordIndex: number, recordId: string }
 *   | RegionSelection} PreviewItem
 * @typedef {{ previewItem: PreviewItem, matchesReady: boolean }} PreviewItemPopup
 *   A popup about a preview item; a Legend row's popup counts its features once `matchesReady`.
 * @typedef {MatchPopupPayload | PreviewItemPopup} MatchFramePopup
 *   What the match popup frame shows (`state.clickedPairwiseMatch`): a match, or a preview item.
 */

/**
 * The open match-frame popup when it is about a preview item, else null.
 * @param {MatchFramePopup | null} opened
 * @returns {PreviewItemPopup | null}
 */
export const previewItemPopupOf = (opened) => (opened && 'previewItem' in opened ? opened : null);

/**
 * The open match-frame popup when it shows a match, else null.
 * @param {MatchFramePopup | null} opened
 * @returns {MatchPopupPayload | null}
 */
export const matchPopupOf = (opened) => (opened && !('previewItem' in opened) ? opened : null);

export const PREVIEW_ITEM_SELECTOR =
  'g[data-legend-key], [data-gbdraw-role="record-definition"], [data-gbdraw-role="record-definition-row"]';

/**
 * The Legend row or record definition an event in the preview reached.
 * @param {EventTarget | null} target @param {Element} svg
 * @returns {{ element: Element, item: PreviewItem } | null}
 */
export const previewItemTarget = (target, svg) => {
  const element = target instanceof Element ? target.closest(PREVIEW_ITEM_SELECTOR) : null;
  if (!element || !svg.contains(element) || element.closest(COMPARISON_LEGEND_SELECTOR)) return null;
  if (element.hasAttribute('data-legend-key')) {
    // `pythonLegendKey` reads the key trimmed.
    const key = pythonLegendKey(element);
    return key ? { element, item: { kind: 'legend-entry', key } } : null;
  }
  const index = (element.getAttribute('data-gbdraw-record-index') ?? '').trim();
  const recordIndex = Number(index);
  if (!index || !Number.isInteger(recordIndex) || recordIndex < 0) return null;
  return { element, item: { kind: 'record-definition', recordIndex, recordId: (element.getAttribute('data-gbdraw-record-id') ?? '').trim() } };
};
