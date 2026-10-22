import assert from 'node:assert/strict';
import test from 'node:test';

import { placePopupNearPoint, placePopupWithin, popupViewportBounds } from '../../gbdraw/web/js/utils/popup-bounds.js';

// OV-371 and OV-384: the first placement of a popup opened at a clicked point.
const desktop = { top: 56, bottom: 968, left: 0, right: 1600 };
const inside = (placed, width, height, bounds) => (
  placed.x >= bounds.left + 12 && placed.x + width <= bounds.right - 12
  && placed.y >= bounds.top + 12 && placed.y + Math.min(height, placed.maxHeight) <= bounds.bottom - 12
);
const covers = (placed, width, height, point) => (
  point.x >= placed.x && point.x <= placed.x + width
  && point.y >= placed.y && point.y <= placed.y + Math.min(height, placed.maxHeight)
);

test('a popup opens beside the clicked point, right first, and inside the bounds', () => {
  for (const point of [{ x: 300, y: 500 }, { x: 900, y: 430 }, { x: 1500, y: 900 }, { x: 800, y: 60 }]) {
    const placed = placePopupNearPoint({ point, width: 720, height: 600, bounds: desktop });
    assert.ok(inside(placed, 720, 600, desktop), JSON.stringify({ point, placed }));
    assert.equal(covers(placed, 720, 600, point), false, JSON.stringify({ point, placed }));
  }
  assert.equal(placePopupNearPoint({ point: { x: 300, y: 500 }, width: 720, height: 600, bounds: desktop }).x, 312);
  assert.equal(placePopupNearPoint({ point: { x: 1200, y: 500 }, width: 720, height: 600, bounds: desktop }).x, 468);
});

test('a popup taller than the bounds gets their full height and never reaches the header', () => {
  const placed = placePopupNearPoint({ point: { x: 900, y: 430 }, width: 720, height: 2000, bounds: desktop });
  assert.deepEqual(placed, { x: 168, y: 68, maxHeight: 888 });
});

test('a point outside the window keeps the popup inside the bounds (below the header)', () => {
  const placed = placePopupNearPoint({ point: { x: -100, y: -100 }, width: 720, height: 520, bounds: desktop });
  assert.deepEqual(placed, { x: 12, y: 68, maxHeight: 888 });
});

test('without room beside the point the popup opens below or above it', () => {
  const narrow = { top: 56, bottom: 1400, left: 0, right: 900 };
  const below = placePopupNearPoint({ point: { x: 450, y: 200 }, width: 720, height: 600, bounds: narrow });
  assert.deepEqual(below, { x: 90, y: 212, maxHeight: 1176 });
  const above = placePopupNearPoint({ point: { x: 450, y: 1300 }, width: 720, height: 600, bounds: narrow });
  assert.deepEqual(above, { x: 90, y: 688, maxHeight: 600 });
  assert.equal(covers(above, 720, 2000, { x: 450, y: 1300 }), false);
});

test('on a phone the popup covers the point at full height when no side holds most of it', () => {
  const phone = { top: 135, bottom: 685, left: 0, right: 390 };
  const placed = placePopupNearPoint({ point: { x: 174, y: 330 }, width: 720, height: 1100, bounds: phone });
  assert.deepEqual(placed, { x: 12, y: 147, maxHeight: 526 });
});

test('without a point the popup is centred in the bounds', () => {
  assert.deepEqual(placePopupNearPoint({ point: null, width: 460, height: 400, bounds: desktop }), { x: 570, y: 312, maxHeight: 644 });
});

// A record definition popup leaves the whole clicked definition uncovered: a
// click on the left of a Circular center text still opens it right of the text.
test('a popup leaves the clicked item rectangle uncovered', () => {
  const avoid = { left: 700, right: 900, top: 450, bottom: 560 };
  const placed = placePopupNearPoint({ point: { x: 720, y: 500, avoid }, width: 400, height: 500, bounds: desktop });
  assert.equal(placed.x, 912);
  const narrow = { ...desktop, right: 1200 };
  const left = placePopupNearPoint({ point: { x: 880, y: 500, avoid }, width: 400, height: 500, bounds: narrow });
  assert.equal(left.x, 700 - 12 - 400);
});

// OV-503: a control inside the Editor drawer opens the popup beside the whole
// drawer. Where the bounds have no room beside it, the popup opens as it does
// without a panel.
test('a popup opens beside the panel that holds the clicked control when the bounds have room', () => {
  const request = { width: 720, height: 520, bounds: desktop };
  const panel = { left: 1223, right: 1583, top: 176, bottom: 955 };
  assert.deepEqual(placePopupNearPoint({ ...request, point: { x: 1562, y: 393, panel } }), { x: 491, y: 133, maxHeight: 823 });
  // A click made with the keyboard has no point of its own.
  assert.equal(placePopupNearPoint({ ...request, point: { x: 0, y: 0, panel } }).x, 491);
  // A panel on the left edge: the popup lies right of it.
  assert.equal(placePopupNearPoint({ ...request, point: { x: 40, y: 393, panel: { ...panel, left: 0, right: 360 } } }).x, 372);

  const cases = [
    // A window narrower than the popup and the drawer together.
    [{ ...desktop, right: 1100 }, { x: 1062, y: 427 }, { left: 723, right: 1083, top: 210, bottom: 755 }],
    // The docked Editor of a one-column page, as wide as the window.
    [{ top: 135, bottom: 844, left: 0, right: 390 }, { x: 360, y: 420 }, { left: 9, right: 381, top: 205, bottom: 601 }]
  ];
  for (const [bounds, point, docked] of cases) {
    assert.deepEqual(
      placePopupNearPoint({ ...request, bounds, point: { ...point, panel: docked } }),
      placePopupNearPoint({ ...request, bounds, point })
    );
  }
});

// A drag, and a window that changed size, keep a popup inside the bounds.
test('a popup is placed the least distance from the wanted place that puts it inside the bounds', () => {
  const size = { width: 320, height: 400 };
  assert.deepEqual(placePopupWithin({ x: 500, y: 300, ...size, bounds: desktop }), { x: 500, y: 300, maxHeight: 0 });
  assert.deepEqual(placePopupWithin({ x: 1500, y: 800, ...size, bounds: desktop }), { x: 1268, y: 556, maxHeight: 0 });
  assert.deepEqual(placePopupWithin({ x: -40, y: 10, ...size, bounds: desktop }), { x: 12, y: 68, maxHeight: 0 });
  // A popup larger than the bounds starts at their top left.
  assert.deepEqual(placePopupWithin({ x: 200, y: 300, width: 2000, height: 2000, bounds: desktop }), { x: 12, y: 68, maxHeight: 0 });
});

// OV-496: the notice banner is a top limit only where it is a banner above the
// popups. On a one-column page it lies below the sidebar: in the lower half of
// the window, or far under it while the page is scrolled above the Result.
test('the notice limits the popup top only while it ends in the upper half under the header', () => {
  const boundsWith = (noticeBottom) => {
    const bottoms = { 'header.app-header': 135, '[data-app-notice]': noticeBottom };
    Object.assign(globalThis, {
      window: { innerWidth: 390, innerHeight: 844 },
      document: {
        querySelector: (selector) => (
          typeof bottoms[selector] === 'number' ? { getBoundingClientRect: () => ({ top: 0, bottom: bottoms[selector] }) } : null
        )
      }
    });
    try {
      return popupViewportBounds();
    } finally {
      Reflect.deleteProperty(globalThis, 'window');
      Reflect.deleteProperty(globalThis, 'document');
    }
  };
  assert.deepEqual(boundsWith(193), { top: 193, bottom: 844, left: 0, right: 390 });
  assert.equal(boundsWith(489).top, 489);
  for (const noticeBottom of [undefined, -2, 100, 490, 700, 843, 4924]) {
    assert.equal(boundsWith(noticeBottom).top, 135, `notice bottom ${noticeBottom}`);
  }
});
