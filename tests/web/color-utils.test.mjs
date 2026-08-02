import assert from 'node:assert/strict';
import {
  colorValueForMode,
  colorValueMode,
  resolveColorToHex,
  toNativeColorInputValue
} from '../../gbdraw/web/js/app/color-utils.js';

assert.equal(colorValueMode(null), 'auto');
assert.equal(colorValueMode(''), 'auto');
assert.equal(colorValueMode('none'), 'none');
assert.equal(colorValueMode('red'), 'color');

assert.equal(toNativeColorInputValue(null, '#aabbcc'), '#aabbcc');
assert.equal(toNativeColorInputValue('none', '#aabbcc'), '#aabbcc');
assert.equal(toNativeColorInputValue('red'), '#ff0000');
assert.equal(toNativeColorInputValue('#AbC'), '#aabbcc');
assert.equal(toNativeColorInputValue('invalid', 'invalid'), '#000000');
assert.equal(resolveColorToHex('rebeccapurple'), '#663399');
assert.equal(resolveColorToHex('seashell'), '#FFF5EE');

assert.equal(colorValueForMode('auto', '#123456'), null);
assert.equal(colorValueForMode('none', '#123456'), 'none');
assert.equal(colorValueForMode('color', 'navy'), '#000080');
