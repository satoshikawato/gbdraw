import assert from 'node:assert/strict';
import test from 'node:test';

import { SOURCE_LOCALE, setUiLocale } from '../../gbdraw/web/js/utils/i18n.js';
import { limitedList, moreLine } from '../../gbdraw/web/js/utils/limited-list.js';

// OV-3495: a list longer than its line limit shows one line fewer and one moreLine for the rest, so it never takes
// more lines than the limit (the Load dialog, the preview's input-summary notice).
test('a list over its line limit keeps the first lines and sets the rest apart', () => {
  const items = ['A', 'B', 'C', 'D', 'E', 'F'];
  assert.deepEqual(limitedList(items, 6), { shown: items, rest: [] });
  assert.deepEqual(limitedList(items, 7), { shown: items, rest: [] });
  assert.deepEqual(limitedList(items, 4), { shown: ['A', 'B', 'C'], rest: ['D', 'E', 'F'] });
  assert.deepEqual(limitedList(items, 1), { shown: [], rest: items });
  assert.deepEqual(limitedList([], 4), { shown: [], rest: [] });
});

test('a limited list and its parts are frozen, and the caller\'s list is not', () => {
  const items = ['A', 'B'];
  for (const list of [limitedList(items, 4), limitedList(items, 1)]) {
    assert.ok(Object.isFrozen(list) && Object.isFrozen(list.shown) && Object.isFrozen(list.rest));
  }
  assert.ok(!Object.isFrozen(items));
});

test('one line stands for the rest, in the shown UI language', async () => {
  assert.equal(moreLine(1), 'and one more');
  assert.equal(moreLine(53), 'and 53 more');
  try {
    await setUiLocale('ja');
    assert.notEqual(moreLine(53), 'and 53 more');
    assert.match(moreLine(53), /53/);
  } finally {
    await setUiLocale(SOURCE_LOCALE);
  }
});
