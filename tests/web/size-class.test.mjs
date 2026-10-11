import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import test from 'node:test';

import { createAutoValueDisplay } from '../../gbdraw/web/js/app/auto-value-display.js';
import { sizeClassForRecordLengths } from '../../gbdraw/web/js/services/track-slot-display.js';
import { withDrawings } from './helpers/drawing-state.mjs';

const packagedThresholds = () => {
  const config = readFileSync(new URL('../../gbdraw/data/config.toml', import.meta.url), 'utf8');
  const table = config.split('[labels.length_threshold]')[1].split('[')[0];
  return Object.fromEntries(
    [...table.matchAll(/^\s*(circular|linear)\s*=\s*(\d+)\s*$/gm)].map(([, mode, bp]) => [mode, Number(bp)])
  );
};

test('the Web size class uses the packaged short/long boundary', () => {
  const { circular } = packagedThresholds();
  assert.equal(circular, 100000);
  assert.equal(sizeClassForRecordLengths([3000, circular - 1]), 'short');
  assert.equal(sizeClassForRecordLengths([3000, circular]), 'long');
  assert.equal(sizeClassForRecordLengths([]), null);
});

test('a Circular record under 100 kb shows the short Auto sizes', () => {
  const autoTextFor = (recordLength) => createAutoValueDisplay(withDrawings({
    mode: { value: 'circular' },
    circularRecordList: { value: [{ record_length: recordLength }] },
    form: {},
    adv: {}
  })).autoValueText('labelFontSize');
  assert.equal(autoTextFor(99_999), '14 px (auto)');
  assert.equal(autoTextFor(100_000), '8 px (auto)');
});
