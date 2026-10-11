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

test('the Web size class uses the packaged short/long boundary of each mode', () => {
  const thresholds = packagedThresholds();
  assert.deepEqual(thresholds, { circular: 75000, linear: 100000 });
  for (const [mode, boundary] of Object.entries(thresholds)) {
    assert.equal(sizeClassForRecordLengths(mode, [3000, boundary - 1]), 'short', mode);
    assert.equal(sizeClassForRecordLengths(mode, [3000, boundary]), 'long', mode);
    assert.equal(sizeClassForRecordLengths(mode, []), null, mode);
  }
});

test('a Circular record under 75 kb shows the short Auto sizes', () => {
  const autoTextFor = (recordLength) => createAutoValueDisplay(withDrawings({
    mode: { value: 'circular' },
    circularRecordList: { value: [{ record_length: recordLength }] },
    form: {},
    adv: {}
  })).autoValueText('labelFontSize');
  assert.equal(autoTextFor(74_999), '14 px (auto)');
  assert.equal(autoTextFor(75_000), '8 px (auto)');
});
