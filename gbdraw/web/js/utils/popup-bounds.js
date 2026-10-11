// @ts-check

// The area a popup or floating panel may cover: below the app header
// (`header.app-header`, OV-384: its buttons stay reachable) and the notice
// banner (`[data-app-notice]`, so its text and Dismiss stay readable), and above
// the app footer (`[data-app-footer]`) when they are on screen, else the viewport. The
// feature and match popups, their drags and resizes, their first placement,
// and the alignment palette read these limits.
/** @returns {number} */
const popupViewportBottom = () => {
  const viewportBottom = Math.max(1, window.innerHeight || 1);
  const footerTop = document.querySelector('[data-app-footer]')?.getBoundingClientRect().top;
  return typeof footerTop === 'number' && footerTop > 0 && footerTop < viewportBottom
    ? footerTop
    : viewportBottom;
};

/** @returns {number} */
const popupViewportTop = () => {
  const bottoms = ['header.app-header', '[data-app-notice]']
    .map((selector) => document.querySelector(selector)?.getBoundingClientRect().bottom)
    .filter(/** @returns {bottom is number} */ (bottom) => typeof bottom === 'number' && bottom > 0);
  return bottoms.length ? Math.min(Math.max(...bottoms), popupViewportBottom() - 1) : 0;
};

/**
 * @typedef {{ top: number, bottom: number, left: number, right: number }} PopupBounds
 * @typedef {{ x: number, y: number, maxHeight: number }} PopupPlacement
 * A clicked point, and the rectangle of the clicked item when the popup must
 * not cover more than the point (a Circular center definition).
 * @typedef {{ x: number, y: number, avoid?: PopupBounds | null }} PopupAnchor
 */

/** @returns {PopupBounds} */
export const popupViewportBounds = () => ({
  top: popupViewportTop(),
  bottom: popupViewportBottom(),
  left: 0,
  right: Math.max(1, window.innerWidth || 1)
});

const clamp = (/** @type {number} */ value, /** @type {number} */ min, /** @type {number} */ max) => (
  Math.min(Math.max(value, min), Math.max(min, max))
);

/**
 * OV-371: the placement of a popup opened at a clicked point. The popup lies
 * fully inside the bounds and beside the point (right, then left), else below
 * or above it, so it does not cover what was clicked. When no side has room
 * for most of it, it covers the point. `height` is the popup's expected
 * height (its content, at most the bounds); `maxHeight` keeps a popup placed
 * above the point from growing over it. Without a point it is centred. With
 * `point.avoid`, the popup stays beside, below, or above that whole rectangle.
 * @param {{
 *   point: PopupAnchor | null,
 *   width: number,
 *   height: number,
 *   bounds: PopupBounds,
 *   minHeight?: number,
 *   margin?: number,
 *   gap?: number
 * }} request
 * @returns {PopupPlacement}
 */
export const placePopupNearPoint = ({ point, width, height, bounds, minHeight = 180, margin = 12, gap = 12 }) => {
  const top = bounds.top + margin;
  const bottom = bounds.bottom - margin;
  const left = bounds.left + margin;
  const right = bounds.right - margin;
  const w = Math.min(width, Math.max(0, right - left));
  const h = Math.min(height, Math.max(0, bottom - top));
  const at = (/** @type {number} */ x, /** @type {number} */ y, end = bottom) => ({ x, y, maxHeight: Math.max(0, end - y) });
  if (!point) return at(clamp((left + right - w) / 2, left, right - w), clamp((top + bottom - h) / 2, top, bottom - h));
  const { x: px, y: py } = point;
  const avoid = {
    left: Math.min(px, point.avoid?.left ?? px), right: Math.max(px, point.avoid?.right ?? px),
    top: Math.min(py, point.avoid?.top ?? py), bottom: Math.max(py, point.avoid?.bottom ?? py)
  };
  const besideY = clamp(py - h / 2, top, bottom - h);
  const rightX = Math.max(left, avoid.right + gap);
  if (rightX + w <= right) return at(rightX, besideY);
  const leftX = Math.min(right - w, avoid.left - gap - w);
  if (leftX >= left) return at(leftX, besideY);
  const x = clamp(px - w / 2, left, right - w);
  const belowY = Math.max(top, avoid.bottom + gap);
  const aboveEnd = Math.min(bottom, avoid.top - gap);
  const below = bottom - belowY;
  const above = aboveEnd - top;
  const useBelow = below >= h || (above < h && below >= above);
  // A side shorter than two thirds of the bounds would leave too little of the
  // popup: it then covers the point at full height (a phone screen).
  if ((useBelow ? below : above) < Math.min(h, Math.max(minHeight, (2 / 3) * (bottom - top)))) return at(x, besideY);
  if (useBelow) return at(x, belowY);
  return at(x, Math.max(top, aboveEnd - h), aboveEnd);
};
