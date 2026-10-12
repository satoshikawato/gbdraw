// @ts-check

// The area a popup or floating panel may cover: below the app header
// (`header.app-header`, OV-384: its buttons stay reachable) and the notice
// banner (`[data-app-notice]`, so its text and Dismiss stay readable), and above
// the app footer (`[data-app-footer]`) when they are on screen, else the viewport. The
// feature and match popups, their drags and resizes, their first placement,
// and the alignment palette read these limits.
// OV-496: the notice is a limit only while it ends in the upper half of the
// area under the header. On a one-column page it lies below the sidebar, lower
// in the window or far under it; a popup may then cover it.
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
  const bottom = popupViewportBottom();
  const bottomOf = (/** @type {string} */ selector) => document.querySelector(selector)?.getBoundingClientRect().bottom ?? 0;
  const header = Math.min(Math.max(0, bottomOf('header.app-header')), bottom - 1);
  const notice = bottomOf('[data-app-notice]');
  return notice > header && notice <= (header + bottom) / 2 ? notice : header;
};

/**
 * @typedef {{ top: number, bottom: number, left: number, right: number }} PopupBounds
 * @typedef {{ x: number, y: number, maxHeight: number }} PopupPlacement
 * A clicked point, and the rectangle of the clicked item when the popup must
 * not cover more than the point (a Circular center definition). `panel` is
 * the rectangle of the panel that holds the clicked control (the Editor
 * drawer, OV-503): the popup opens beside the whole panel when it can.
 * @typedef {{ x: number, y: number, avoid?: PopupBounds | null, panel?: PopupBounds | null }} PopupAnchor
 * Where a popup was opened from: the pointer event's point, and the panel of
 * a `PopupAnchor` when a control inside one opened the popup.
 * @typedef {{ clientX?: number, clientY?: number, panel?: PopupBounds | null }} PopupOrigin
 */

/** @returns {PopupBounds} */
export const popupViewportBounds = () => ({
  top: popupViewportTop(),
  bottom: popupViewportBottom(),
  left: 0,
  right: Math.max(1, window.innerWidth || 1)
});

// The popup bounds inside a rectangle of the page (the Result's canvas).
/** @param {PopupBounds} box @returns {PopupBounds} */
export const popupBoundsWithin = (box) => {
  const bounds = popupViewportBounds();
  return {
    top: Math.max(bounds.top, box.top), bottom: Math.min(bounds.bottom, box.bottom),
    left: Math.max(bounds.left, box.left), right: Math.min(bounds.right, box.right)
  };
};

// A press that moved farther than this before its release is a drag (a pan of
// the figure), not a click that opens or closes a popup of the Result.
const CLICK_SLOP_PX = 4;
/**
 * @typedef {{ clientX: number, clientY: number }} PressPoint Where a press began.
 * @param {PressPoint | null} press @param {PressPoint} release
 */
export const isClickAtPress = (press, release) =>
  !press || Math.hypot(release.clientX - press.clientX, release.clientY - press.clientY) <= CLICK_SLOP_PX;

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
 * With `point.panel`, it lies beside the panel when the bounds have room there;
 * without that room the panel changes nothing.
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
  for (const keep of point.panel ? [point.panel, avoid] : [avoid]) {
    const rightX = Math.max(left, keep.right + gap);
    if (rightX + w <= right) return at(rightX, besideY);
    const leftX = Math.min(right - w, keep.left - gap - w);
    if (leftX >= left) return at(leftX, besideY);
  }
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
