import assert from 'node:assert/strict';
import test from 'node:test';

import { placePopupNearPoint } from '../../gbdraw/web/js/utils/popup-bounds.js';

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
