// OV-2792: a list typed into one text field (ring names, the Label filtering
// Blacklist) is cut at the comma, at the Arabic comma ، (U+060C) that fa and ar
// keyboards type and their catalogs write for "comma", and at a line break.
// Text without U+060C reads as it did before.
import assert from 'node:assert/strict';
import test from 'node:test';

const {
  endsWithTypedListSeparator,
  splitTypedList,
  splitTypedListKeepingSeparators,
  trimTypedListItem
} = await import('../../gbdraw/web/js/utils/typed-list.js');

test('a typed list is cut at the comma, the Arabic comma and a line break', () => {
  assert.deepEqual(splitTypedList('a،b'), ['a', 'b']);
  assert.deepEqual(splitTypedList(' حلقه الف ، حلقه ب '), ['حلقه الف', 'حلقه ب']);
  assert.deepEqual(splitTypedList('a, b،c\nd\r\ne\rf'), ['a', 'b', 'c', 'd', 'e', 'f']);
  // Runs of separators and blank items give no item.
  assert.deepEqual(splitTypedList(',،\n a ,, \n،b,'), ['a', 'b']);
  assert.deepEqual(splitTypedList(''), []);
  assert.deepEqual(splitTypedList(undefined), []);
  assert.deepEqual(splitTypedList(null), []);
});

test('the semicolon, the Arabic semicolon and spaces are no separators', () => {
  assert.deepEqual(splitTypedList('a;b؛c d'), ['a;b؛c d']);
});

test('text without the Arabic comma reads as it did before', () => {
  const before = (text) => String(text || '').split(/[\r\n,]+/).map((item) => item.trim()).filter(Boolean);
  for (const text of ['hypothetical, uncharacterized, putative, unknown', 'Ring A\nRing B', 'a,,b', ' , ', 'x\r\ny', '']) {
    assert.deepEqual(splitTypedList(text), before(text), JSON.stringify(text));
  }
});

test('the text keeps its separators: items at even places, one separator at each odd place', () => {
  assert.deepEqual(splitTypedListKeepingSeparators('a, b،c\nd\r\ne'), ['a', ',', ' b', '،', 'c', '\n', 'd', '\r\n', 'e']);
  assert.deepEqual(splitTypedListKeepingSeparators('a'), ['a']);
  assert.equal(splitTypedListKeepingSeparators('x,،\ny\r\n').join(''), 'x,،\ny\r\n');
});

test('a text that ends with a separator', () => {
  for (const text of ['a,', 'a،', 'a\n', 'a\r']) assert.equal(endsWithTypedListSeparator(text), true, JSON.stringify(text));
  for (const text of ['', 'a', 'a, b', 'a;']) assert.equal(endsWithTypedListSeparator(text), false, JSON.stringify(text));
});

test('an item loses the characters Python trims too, and keeps the rest', () => {
  // A byte order mark (a Notepad file), U+0085 and U+001C, which Python strips and trim() kept.
  assert.deepEqual(splitTypedList('\uFEFFa,\u0085b\u001C'), ['a', 'b']);
  assert.equal(trimTypedListItem(' \u3000a b\u00A0\t'), 'a b');
  assert.equal(trimTypedListItem('a\u200Bb'), 'a\u200Bb');
});
