import assert from 'node:assert/strict';
import test from 'node:test';

import {
  MAX_POPUP_PLACES, placePopupNearPoint, placePopupWithin, popupPlaces, popupViewportBounds, primaryPressAfter
} from '../../gbdraw/web/js/utils/popup-bounds.js';

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

// Select region (D-28): the Region popup keeps clear of the band (`avoid`), the
// Legend and the label boxes (`obstacles`); when no place inside the bounds is
// clear of them all, it takes the place that covers the least of them.
const overlap = (placed, width, height, box) => (
  Math.max(0, Math.min(placed.x + width, box.right) - Math.max(placed.x, box.left))
  * Math.max(0, Math.min(placed.y + height, box.bottom) - Math.max(placed.y, box.top))
);

test('a popup keeps clear of the obstacles it is given, and opens where it would without them when that place is clear', () => {
  const avoid = { left: 600, right: 700, top: 400, bottom: 500 };
  const point = { x: 650, y: 450, avoid };
  const request = { width: 400, height: 300, bounds: desktop };
  const free = placePopupNearPoint({ ...request, point });
  assert.equal(free.x, 712);
  // An obstacle elsewhere changes nothing.
  assert.deepEqual(placePopupNearPoint({ ...request, point: { ...point, obstacles: [{ left: 0, right: 100, top: 60, bottom: 160 }] } }), free);
  // The Legend right of the band: the popup goes to the nearest clear place.
  const legend = { left: 720, right: 1300, top: 200, bottom: 700 };
  const placed = placePopupNearPoint({ ...request, point: { ...point, obstacles: [legend] } });
  for (const box of [legend, avoid]) assert.equal(overlap(placed, 400, 300, box), 0, JSON.stringify({ placed, box }));
  assert.ok(placed.x >= desktop.left + 12 && placed.x + 400 <= desktop.right - 12);
  assert.ok(placed.y >= desktop.top + 12 && placed.y + 300 <= desktop.bottom - 12);
  assert.ok(placed.maxHeight >= 300);
  assert.equal(placed.x, 600 - 12 - 400, 'left of the band, the nearest clear side');
});

test('with no clear place the popup covers the least of the obstacles', () => {
  const bounds = { top: 0, bottom: 600, left: 0, right: 1000 };
  // The left half is covered, and the right half is clear only 250 px down from 350.
  const wide = { left: 0, right: 500, top: 0, bottom: 600 };
  const upper = { left: 500, right: 1000, top: 0, bottom: 350 };
  const point = { x: 250, y: 300, obstacles: [wide, upper] };
  const placed = placePopupNearPoint({ point, width: 400, height: 300, bounds });
  assert.deepEqual(placed, { x: 512, y: 288, maxHeight: 300 });
  assert.equal(overlap(placed, 400, 300, wide), 0);
  assert.equal(overlap(placed, 400, 300, upper), (350 - 288) * 400);
});

// D-30 (e): a Result can draw thousands of labels; the places tried stay a fixed number, and each place is measured
// only against the boxes near it, so a release opens the popup at once.
test('the places a popup tries are a fixed number, however many labels it keeps clear of', () => {
  const bounds = { top: 56, bottom: 968, left: 0, right: 1600 };
  const ring = (count) => Array.from({ length: count }, (_, index) => {
    const angle = (index / count) * 2 * Math.PI;
    const x = 800 + 420 * Math.cos(angle);
    const y = 500 + 420 * Math.sin(angle);
    return { left: x - 30, right: x + 30, top: y - 5, bottom: y + 5 };
  });
  const legend = { left: 1300, right: 1500, top: 100, bottom: 300 };
  const point = { x: 900, y: 300, avoid: { left: 850, right: 1000, top: 200, bottom: 320 } };
  const request = { width: 360, height: 300, bounds };
  const beside = { x: 1012, y: 150, maxHeight: 800 };
  const few = popupPlaces({ ...request, point: { ...point, obstacles: ring(4) } }, beside, ring(4));
  const many = popupPlaces({ ...request, point: { ...point, obstacles: ring(1_500) } }, beside, ring(1_500));
  assert.ok(many.length <= MAX_POPUP_PLACES, `${many.length} places`);
  assert.ok(few.length <= MAX_POPUP_PLACES && MAX_POPUP_PLACES <= 21 * 21);
  const obstacles = [legend, ...ring(1_500)];
  const placed = placePopupNearPoint({ ...request, point: { ...point, obstacles } });
  // Inside the ring of labels there is room: the popup is clear of every label, the Legend and the band.
  for (const box of [...obstacles, point.avoid]) assert.equal(overlap(placed, 360, 300, box), 0, JSON.stringify({ placed, box }));
});

// D-30 (c): a selected region's popup stays above the preview toolbar; `within` cuts the bounds.
test('a popup stays inside the rectangle its anchor gives besides the bounds', () => {
  const toolbarTop = 880;
  const within = { top: -Infinity, left: -Infinity, right: Infinity, bottom: toolbarTop };
  for (const y of [300, 700, 870, 950]) {
    const placed = placePopupNearPoint({ point: { x: 700, y, within }, width: 400, height: 500, bounds: desktop });
    assert.ok(placed.y + Math.min(500, placed.maxHeight) <= toolbarTop - 12, JSON.stringify({ y, placed }));
    assert.ok(placed.y >= desktop.top + 12);
  }
  const obstacles = [{ left: 0, right: 1600, top: 56, bottom: 400 }];
  const clear = placePopupNearPoint({ point: { x: 700, y: 600, within, obstacles }, width: 400, height: 300, bounds: desktop });
  assert.ok(clear.y >= 412 && clear.y + 300 <= toolbarTop - 12, JSON.stringify(clear));
});

// OV-2704: a click is measured from the press of the primary button, which a
// chorded press (another button already held) sends as a pointermove.
test('the press a click is measured from is that of the primary button', () => {
  const kept = { clientX: 5, clientY: 5 };
  const event = (type, button, buttons) => ({ type, button, buttons, clientX: 100, clientY: 80 });
  const at = { clientX: 100, clientY: 80 };
  assert.deepEqual(primaryPressAfter(kept, event('pointerdown', 0, 1)), at);
  assert.deepEqual(primaryPressAfter(null, event('pointerdown', 0, 0)), at, 'a synthetic press sets no buttons');
  assert.equal(primaryPressAfter(kept, event('pointerdown', 2, 2)), kept, 'another button does not click');
  assert.deepEqual(primaryPressAfter(kept, event('pointermove', 0, 3)), at, 'the primary button added to a held one');
  assert.equal(primaryPressAfter(kept, event('pointermove', 0, 2)), kept, 'the primary button released, another held');
  assert.equal(primaryPressAfter(kept, event('pointermove', -1, 1)), kept, 'a move with no button change');
});
