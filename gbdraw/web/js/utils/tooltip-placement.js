// @ts-check
// Where the house tooltip (`.ui-tooltip`, css/tailwind.input.css) sits: centred
// above its anchor at viewport coordinates, kept 10 px inside the viewport's
// sides, and below the anchor when the anchor is within 60 px of the top
// (header and browser frame) or the tooltip would not fit above (OV-2449: a
// long timing note). A help tip (HelpTip, components.js) and the whole hint of
// a cut field (app/field-hint.js) are placed here.

/**
 * The anchor's box in viewport coordinates (a DOMRect has these).
 * @typedef {{ top: number, bottom: number, left: number, width: number }} TooltipAnchor
 */
/**
 * The tooltip's position as CSS values of its `position: fixed` box.
 * @typedef {{ top: string, left: string, transform: string }} TooltipPlacement
 */

/** The tooltip's widest box (max-w-[250px]) with a margin. */
const TOOLTIP_MAX_WIDTH = 260;
/** The gap between the anchor and the tooltip. */
const GAP = 12;

/**
 * @param {TooltipAnchor} anchor
 * @param {number} tooltipHeight The shown tooltip's height, 0 before it is shown.
 * @param {number} viewportWidth
 * @returns {TooltipPlacement}
 */
export const placeTooltip = (anchor, tooltipHeight, viewportWidth) => {
  const halfWidth = TOOLTIP_MAX_WIDTH / 2;
  let left = anchor.left + anchor.width / 2;
  if (left < halfWidth + 10) left = halfWidth + 10;
  else if (left > viewportWidth - halfWidth - 10) left = viewportWidth - halfWidth - 10;
  const below = anchor.top < 60 || anchor.top - GAP - tooltipHeight < 8;
  return {
    top: `${below ? anchor.bottom + GAP : anchor.top - GAP}px`,
    left: `${left}px`,
    transform: below ? 'translate(-50%, 0)' : 'translate(-50%, -100%)'
  };
};
