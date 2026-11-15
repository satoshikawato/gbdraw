// OV-2792: a list typed into one text field (ring names, the Label filtering
// Blacklist) is cut at the comma, at the Arabic comma ، (U+060C) that fa and ar
// keyboards type and their catalogs write for "comma", and at a line break.
// Each old Web reader gives the items it gave before, except for the characters listed below.
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

// The three Web readers before OV-2792, each with its own cut and JavaScript's trim().
const OLD_READERS = {
  'ring names (parseConservationLabelText)': (text) => String(text || '').split(/[\n,]+/)
    .map((item) => item.trim()).filter(Boolean),
  'Blacklist text in the request (labelBlacklistOverride)': (text) => String(text || '').split(/[,\n]/)
    .map((item) => item.trim()).filter(Boolean),
  'Blacklist file words (parseBlacklistWords)': (text) => String(text || '').split(/[\r\n,]+/)
    .map((item) => item.trim()).filter(Boolean)
};
const CR = '\r';
const ARABIC_COMMA = '\u060C';
// The characters an item now loses at its ends that trim() kept (Python's strip takes them).
const NEWLY_TRIMMED = ['\u001C', '\u001D', '\u001E', '\u001F', '\u0085'];
// Every character that reads differently from before, for each old reader: the Arabic comma and
// the newly trimmed characters for all three, and a lone CR, which cuts two items now, for the two
// readers that cut only at LF. A textarea holds no lone CR (browsers write LF), so only a
// hand-edited Session reaches that one.
const CHANGED = {
  'ring names (parseConservationLabelText)': [CR, ...NEWLY_TRIMMED, ARABIC_COMMA],
  'Blacklist text in the request (labelBlacklistOverride)': [CR, ...NEWLY_TRIMMED, ARABIC_COMMA],
  'Blacklist file words (parseBlacklistWords)': [...NEWLY_TRIMMED, ARABIC_COMMA]
};
const BMP = Array.from({ length: 0x10000 }, (_, code) => code)
  .filter((code) => code < 0xD800 || code > 0xDFFF).map((code) => String.fromCharCode(code));

test('each old Web reader gives the items it gave before, except for the characters listed as changed', () => {
  for (const [reader, before] of Object.entries(OLD_READERS)) {
    for (const text of ['hypothetical, uncharacterized, putative, unknown', 'Ring A\nRing B', 'a,,b', ' , ', 'x\r\ny', '']) {
      assert.deepEqual(splitTypedList(text), before(text), `${reader}: ${JSON.stringify(text)}`);
    }
    // Every BMP character between two items and around them: only the listed ones read otherwise.
    const differing = BMP.filter((c) => [`a${c}b`, `${c}a${c},${c}b${c}`]
      .some((text) => JSON.stringify(splitTypedList(text)) !== JSON.stringify(before(text))));
    assert.deepEqual(differing, CHANGED[reader], reader);
  }
});

test('the intended differences from the old Web readers', () => {
  const lfOnly = ['ring names (parseConservationLabelText)', 'Blacklist text in the request (labelBlacklistOverride)'];
  /** @type {{ text: string, readers: string[], before: string[], now: string[] }[]} */
  const cases = [
    // A lone CR (a hand-edited Session) was inside an item for the two readers that cut at LF only.
    { text: 'a\rb', readers: lfOnly, before: ['a\rb'], now: ['a', 'b'] },
    // U+0085 and U+001C-U+001F at an item's end stayed; Python's strip takes them, and so does the Web now.
    ...NEWLY_TRIMMED.map((c) => ({ text: `a${c}, ${c}b`, readers: Object.keys(OLD_READERS), before: [`a${c}`, `${c}b`], now: ['a', 'b'] })),
    // The Arabic comma was inside an item for all three.
    { text: 'a\u060Cb', readers: Object.keys(OLD_READERS), before: ['a\u060Cb'], now: ['a', 'b'] }
  ];
  for (const { text, readers, before, now } of cases) {
    for (const reader of readers) assert.deepEqual(OLD_READERS[reader](text), before, `${reader}: ${JSON.stringify(text)}`);
    assert.deepEqual(splitTypedList(text), now, JSON.stringify(text));
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
