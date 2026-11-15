import assert from 'node:assert/strict';
import test from 'node:test';
import { asciiNumberText } from '../../gbdraw/web/js/utils/number-field-digits.js';

test('each digit of the four digit sets reads as its 0-9 digit', () => {
  for (const zero of [0x0660, 0x06F0, 0x0966, 0x09E6]) {
    const digits = Array.from({ length: 10 }, (_, value) => String.fromCodePoint(zero + value)).join('');
    assert.equal(asciiNumberText(digits), '0123456789', `digits from U+${zero.toString(16).toUpperCase()}`);
  }
});

test('the Arabic decimal separator reads as a point; other characters stay', () => {
  assert.equal(asciiNumberText('۱٫۵'), '1.5');
  assert.equal(asciiNumberText('-12.5e3'), '-12.5e3');
  assert.equal(asciiNumberText('۱٬۲۳۴'), '1٬234', 'the Arabic thousands separator is not a decimal point');
  assert.equal(asciiNumberText('abc'), 'abc');
});
