// Which drawings Save writes and which one Load opens (R16-8, DESIGN §4.4).
import test from 'node:test';
import assert from 'node:assert/strict';

import { savedDrawingsOf } from '../../gbdraw/web/js/services/session-drawings.js';

const tabs = (...ids) => ids.map((id) => ({ id }));
/** @param {...string} used */
const usedOnly = (...used) => (/** @type {{ id: string }} */ drawing) => used.includes(drawing.id);
/** @param {{ id: string }[]} list @param {string} id */
const shownOf = (list, id) => /** @type {const} */ ([list, /** @type {{ id: string }} */ (list.find((tab) => tab.id === id))]);
const idsOf = (/** @type {{ drawings: readonly { id: string }[] }} */ saved) => saved.drawings.map(({ id }) => id);

test('a used shown drawing is saved and opens', () => {
  const saved = savedDrawingsOf(...shownOf(tabs('circular', 'linear', 'linear-2'), 'linear'), usedOnly('circular', 'linear'));
  assert.deepEqual(idsOf(saved), ['circular', 'linear']);
  assert.equal(saved.activeDrawingId, 'linear');
});

test('an unused shown drawing is not saved; the next saved tab opens', () => {
  const saved = savedDrawingsOf(...shownOf(tabs('circular', 'linear', 'linear-2'), 'linear'), usedOnly('circular', 'linear-2'));
  assert.deepEqual(idsOf(saved), ['circular', 'linear-2']);
  assert.equal(saved.activeDrawingId, 'linear-2');
});

test('an unused shown last tab opens the previous saved tab', () => {
  const saved = savedDrawingsOf(...shownOf(tabs('circular', 'linear', 'linear-2'), 'linear-2'), usedOnly('circular'));
  assert.deepEqual(idsOf(saved), ['circular']);
  assert.equal(saved.activeDrawingId, 'circular');
});

test('when no drawing is used, the shown one is saved', () => {
  const saved = savedDrawingsOf(...shownOf(tabs('circular', 'linear'), 'linear'), usedOnly());
  assert.deepEqual(idsOf(saved), ['linear']);
  assert.equal(saved.activeDrawingId, 'linear');
});

test('a shown drawing is found in the list by its ID', () => {
  const saved = savedDrawingsOf(tabs('circular', 'linear'), { id: 'linear' }, usedOnly('circular'));
  assert.deepEqual(idsOf(saved), ['circular']);
  assert.equal(saved.activeDrawingId, 'circular');
});
