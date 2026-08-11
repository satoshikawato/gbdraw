import assert from 'node:assert/strict';
import {
  buildDefaultColorOverrideTsv,
  buildPaletteColorOverrideRows,
  colorValueForMode,
  colorValueMode,
  interpolateColor,
  normalizePaletteColors,
  resolveCollinearMatchColor,
  resolvePairwiseLegendGradientColorKeys,
  toNativeColorInputValue
} from '../../gbdraw/web/js/app/color-utils.js';

assert.equal(colorValueMode(null), 'auto');
assert.equal(colorValueMode(''), 'auto');
assert.equal(colorValueMode('none'), 'none');
assert.equal(colorValueMode('red'), 'color');

assert.equal(toNativeColorInputValue(null, '#aabbcc'), '#aabbcc');
assert.equal(toNativeColorInputValue('none', '#aabbcc'), '#aabbcc');
assert.equal(toNativeColorInputValue('#AbC'), '#aabbcc');
assert.equal(toNativeColorInputValue('invalid', 'invalid'), '#000000');

assert.equal(colorValueForMode('auto', '#123456'), null);
assert.equal(colorValueForMode('none', '#123456'), 'none');
assert.equal(colorValueForMode('color', '#000080'), '#000080');

const colors = {
  pairwise_match_min: '#ffffff',
  pairwise_match_max: '#000000',
  collinear_block_plus_min: '#eeeeee',
  collinear_block_plus: '#808080',
  collinear_block_minus_min: '#ffeeee',
  collinear_block_minus: '#ff0000'
};
assert.equal(resolveCollinearMatchColor({
  blockId: 'block_1',
  colorMode: 'orientation_identity',
  orientation: 'minus',
  identityFactor: 0.5,
  colors
}), interpolateColor('#ffeeee', '#ff0000', 0.5));
assert.equal(resolveCollinearMatchColor({
  blockId: 'block_1',
  colorMode: 'orientation',
  orientation: 'minus',
  identityFactor: 0.5,
  colors
}), '#ff0000');
assert.equal(resolveCollinearMatchColor({
  blockId: 'block_1',
  colorMode: 'average_identity',
  orientation: 'minus',
  identityFactor: 0.5,
  colors
}), null);
assert.deepEqual(resolvePairwiseLegendGradientColorKeys('Collinear'), {
  minKey: 'collinear_block_plus_min',
  maxKey: 'collinear_block_plus'
});
assert.deepEqual(resolvePairwiseLegendGradientColorKeys('Inverted'), {
  minKey: 'collinear_block_minus_min',
  maxKey: 'collinear_block_minus'
});

assert.equal(
  normalizePaletteColors({ collinear_block_plus_max: '#123456' })
    .collinear_block_plus,
  '#123456'
);
const paletteColors = normalizePaletteColors({
  CDS: '#54bcf8',
  rRNA: '#71ee7d',
  pairwise_match_min: '#FFE7E7',
  pairwise_match_max: '#FF7272'
});
assert.deepEqual(buildPaletteColorOverrideRows({
  colors: {
    CDS: '#54bcf8',
    rRNA: '#71EE7D',
    pairwise_match_min: '#ffe7e7',
    pairwise_match_max: '#ff7272'
  },
  paletteColors
}), []);
assert.equal(buildDefaultColorOverrideTsv({
  colors: { CDS: '#000000', rRNA: '#71ee7d', custom_feature: '#abcdef' },
  paletteColors
}), 'CDS\t#000000\ncustom_feature\t#abcdef');
