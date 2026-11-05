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
 * `obstacles` are rectangles the popup keeps clear of as well (the Legend and
 * the labels of the Result beside a selected region, D-28). `within` is a
 * rectangle the popup stays inside besides the bounds (a selected region's
 * popup stays above the preview toolbar, D-30).
 * @typedef {{
 *   x: number, y: number, avoid?: PopupBounds | null, panel?: PopupBounds | null,
 *   obstacles?: readonly PopupBounds[], within?: PopupBounds | null
 * }} PopupAnchor
 * Where a popup was opened from: the pointer event's point, and the panel,
 * obstacles, and limit of a `PopupAnchor`.
 * @typedef {{
 *   clientX?: number, clientY?: number, panel?: PopupBounds | null,
 *   obstacles?: readonly PopupBounds[], within?: PopupBounds | null
 * }} PopupOrigin
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
 * The placement nearest to (x, y) at which a popup of this size lies inside
 * the bounds: where a drag may leave it (OV-384: below the header), and where
 * a window that changed size leaves it. A popup larger than the bounds starts
 * at their top left. Its height then follows the bounds (`maxHeight` 0).
 * @param {{ x: number, y: number, width: number, height: number, bounds: PopupBounds, margin?: number }} request
 * @returns {PopupPlacement}
 */
export const placePopupWithin = ({ x, y, width, height, bounds, margin = 12 }) => ({
  x: clamp(x, bounds.left + margin, bounds.right - width - margin),
  y: clamp(y, bounds.top + margin, bounds.bottom - height - margin),
  maxHeight: 0
});

/**
 * OV-371: the placement of a popup opened at a clicked point. The popup lies
 * fully inside the bounds and beside the point (right, then left), else below
 * or above it, so it does not cover what was clicked. When no side has room
 * for most of it, it covers the point. `height` is the popup's expected
 * height (its content, at most the bounds); `maxHeight` keeps a popup placed
 * above the point from growing over it. Without a point it is centred. With
 * `point.avoid`, the popup stays beside, below, or above that whole rectangle.
 * With `point.panel`, it lies beside the panel when the bounds have room there;
 * without that room the panel changes nothing. With `point.obstacles`, that
 * placement stands only when it keeps `gap` clear of them and of the point
 * (`avoid`); else the popup takes the clear place nearest to the point, or,
 * when none is clear, the place that covers the least area of them. With
 * `point.within`, the bounds are cut to that rectangle first.
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
export const placePopupNearPoint = (request) => {
  const within = request.point?.within;
  const bounded = within ? { ...request, bounds: {
    top: Math.max(request.bounds.top, within.top), bottom: Math.min(request.bounds.bottom, within.bottom),
    left: Math.max(request.bounds.left, within.left), right: Math.min(request.bounds.right, within.right)
  } } : request;
  const placed = placeBesidePoint(bounded);
  const obstacles = bounded.point?.obstacles;
  return obstacles?.length ? placeClearOf(bounded, placed, obstacles) : placed;
};

/**
 * The placement beside the point that `placePopupNearPoint` describes, without its obstacles.
 * @param {Parameters<typeof placePopupNearPoint>[0]} request
 * @returns {PopupPlacement}
 */
const placeBesidePoint = ({ point, width, height, bounds, minHeight = 180, margin = 12, gap = 12 }) => {
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

/** @param {PopupBounds} a @param {PopupBounds} b */
const overlapArea = (a, b) => Math.max(0, Math.min(a.right, b.right) - Math.max(a.left, b.left))
  * Math.max(0, Math.min(a.bottom, b.bottom) - Math.max(a.top, b.top));

// The places tried: a fixed grid over the bounds, the place beside the point, and the edges beside the point (with
// its `avoid` box) and beside the boxes nearest to it. Their number never depends on how many boxes there are (D-30).
const GRID_STOPS = 9;
const NEAREST_BOXES = 4;
/** The most places one placement tries. */
export const MAX_POPUP_PLACES = (GRID_STOPS + 2 + 2 * (NEAREST_BOXES + 1)) ** 2;
// The boxes are indexed by the cells of a grid over the bounds, so a place is measured only against the boxes of
// the cells it covers.
const INDEX_CELLS = 8;

/**
 * The bounds a popup of the request's size may start in, and the boxes it keeps `gap` from.
 * @param {Parameters<typeof placePopupNearPoint>[0]} request
 * @param {readonly PopupBounds[]} obstacles
 */
const placementFrame = ({ point, width, height, bounds, margin = 12, gap = 12 }, obstacles) => {
  const top = bounds.top + margin;
  const bottom = bounds.bottom - margin;
  const left = bounds.left + margin;
  const right = bounds.right - margin;
  const w = Math.min(width, Math.max(0, right - left));
  const h = Math.min(height, Math.max(0, bottom - top));
  const px = point?.x ?? (left + right) / 2;
  const py = point?.y ?? (top + bottom) / 2;
  const avoid = {
    left: Math.min(px, point?.avoid?.left ?? px), right: Math.max(px, point?.avoid?.right ?? px),
    top: Math.min(py, point?.avoid?.top ?? py), bottom: Math.max(py, point?.avoid?.bottom ?? py)
  };
  const keep = [avoid, ...obstacles].map((box) => ({
    left: box.left - gap, right: box.right + gap, top: box.top - gap, bottom: box.bottom + gap
  }));
  return { top, bottom, left, right, w, h, px, py, keep };
};

/**
 * The places `placePopupNearPoint` tries for a popup that must keep clear of `obstacles`: at most `MAX_POPUP_PLACES`,
 * whatever the number of obstacles.
 * @param {Parameters<typeof placePopupNearPoint>[0]} request
 * @param {PopupPlacement} besidePoint
 * @param {readonly PopupBounds[]} obstacles
 * @returns {{ x: number, y: number }[]}
 */
export const popupPlaces = (request, besidePoint, obstacles) => {
  const { top, bottom, left, right, w, h, px, py, keep } = placementFrame(request, obstacles);
  const [around, ...boxes] = keep;
  const distanceTo = (/** @type {PopupBounds} */ box) => Math.hypot(
    Math.max(box.left - px, 0, px - box.right), Math.max(box.top - py, 0, py - box.bottom)
  );
  const nearest = [around, ...boxes.map((box) => ({ box, distance: distanceTo(box) }))
    .sort((a, b) => a.distance - b.distance).slice(0, NEAREST_BOXES).map(({ box }) => box)];
  /** @param {number} from @param {number} to @param {readonly number[]} edges */
  const stops = (from, to, edges) => {
    const end = Math.max(from, to);
    const values = new Set(edges.map((value) => clamp(value, from, end)));
    for (let step = 0; step < GRID_STOPS; step += 1) values.add(from + ((end - from) * step) / (GRID_STOPS - 1));
    return [...values];
  };
  const xs = stops(left, right - w, [px - w / 2, besidePoint.x, ...nearest.flatMap((box) => [box.right, box.left - w])]);
  const ys = stops(top, bottom - h, [py - h / 2, besidePoint.y, ...nearest.flatMap((box) => [box.bottom, box.top - h])]);
  return xs.flatMap((x) => ys.map((y) => ({ x, y })));
};

/**
 * `placePopupNearPoint` with obstacles: `besidePoint` when it is clear of them and of the point, else the clear place
 * nearest to the point, else the place that covers the least of them (the nearer one on a tie). A place is clear when
 * it keeps `gap` from every box. The places tried are `popupPlaces`.
 * @param {Parameters<typeof placePopupNearPoint>[0]} request
 * @param {PopupPlacement} besidePoint
 * @param {readonly PopupBounds[]} obstacles
 * @returns {PopupPlacement}
 */
const placeClearOf = (request, besidePoint, obstacles) => {
  if (!request.point) return besidePoint;
  const { top, bottom, left, right, w, h, px, py, keep } = placementFrame(request, obstacles);
  const cellWidth = Math.max(1, (right - left) / INDEX_CELLS);
  const cellHeight = Math.max(1, (bottom - top) / INDEX_CELLS);
  const column = (/** @type {number} */ x) => clamp(Math.floor((x - left) / cellWidth), 0, INDEX_CELLS - 1);
  const row = (/** @type {number} */ y) => clamp(Math.floor((y - top) / cellHeight), 0, INDEX_CELLS - 1);
  /** @type {number[][]} The boxes that reach each cell, by index into `keep`. */
  const cells = Array.from({ length: INDEX_CELLS * INDEX_CELLS }, () => []);
  keep.forEach((box, index) => {
    if (box.right <= left || box.left >= right || box.bottom <= top || box.top >= bottom) return;
    for (let r = row(box.top); r <= row(box.bottom); r += 1) {
      for (let c = column(box.left); c <= column(box.right); c += 1) cells[r * INDEX_CELLS + c].push(index);
    }
  });
  const seen = new Int32Array(keep.length).fill(-1);
  let measured = 0;
  /** @param {number} x @param {number} y @param {number} limit @returns {number} The area covered, counted up to `limit`. */
  const covered = (x, y, limit) => {
    const rect = { left: x, right: x + w, top: y, bottom: y + h };
    const stamp = measured;
    measured += 1;
    let area = 0;
    for (let r = row(rect.top); r <= row(rect.bottom); r += 1) {
      for (let c = column(rect.left); c <= column(rect.right); c += 1) {
        for (const index of cells[r * INDEX_CELLS + c]) {
          if (seen[index] === stamp) continue;
          seen[index] = stamp;
          area += overlapArea(rect, keep[index]);
          if (area > limit) return area;
        }
      }
    }
    return area;
  };
  if (covered(besidePoint.x, besidePoint.y, 0) === 0) return besidePoint;
  /** @param {number} x @param {number} y */
  const distance = (x, y) => Math.hypot(Math.max(x - px, 0, px - x - w), Math.max(y - py, 0, py - y - h));
  const places = popupPlaces(request, besidePoint, obstacles).map(({ x, y }) => ({ x, y, distance: distance(x, y) }))
    .sort((a, b) => a.distance - b.distance);
  let best = places[0];
  let bestArea = Infinity;
  for (const place of places) {
    const area = covered(place.x, place.y, bestArea);
    if (area < bestArea) {
      best = place;
      bestArea = area;
      if (area === 0) break;
    }
  }
  return { x: best.x, y: best.y, maxHeight: Math.max(0, bottom - best.y) };
};
