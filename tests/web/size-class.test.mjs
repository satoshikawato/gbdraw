import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import test from 'node:test';

import { createAutoValueDisplay } from '../../gbdraw/web/js/app/auto-value-display.js';
import { sizeClassForRecordLengths } from '../../gbdraw/web/js/services/track-slot-display.js';
import { withDrawings } from './helpers/drawing-state.mjs';

test('the Web size class uses the packaged short/long boundary of each mode', () => {
  // The boundaries reach the Web through mode-profiles.generated.js (labels.length_threshold).
  for (const [mode, boundary] of /** @type {const} */ ([['circular', 75_000], ['linear', 100_000]])) {
    assert.equal(sizeClassForRecordLengths(mode, [3000, boundary - 1]), 'short', mode);
    assert.equal(sizeClassForRecordLengths(mode, [3000, boundary]), 'long', mode);
    assert.equal(sizeClassForRecordLengths(mode, []), null, mode);
  }
});

test('the Circular label font note follows the shared Auto size vectors, mixed size classes included', () => {
  const vectors = JSON.parse(readFileSync(new URL('../fixtures/auto_size_vectors.json', import.meta.url), 'utf8'));
  const cases = vectors.autoValues.filter((entry) => entry.mode === 'circular' && 'labels.font_size' in entry.expected);
  assert.ok(cases.some((entry) => entry.recordLengths.length > 1), 'a vector with records of both size classes');
  for (const entry of cases) {
    const display = createAutoValueDisplay(withDrawings({
      mode: { value: 'circular' },
      circularRecordList: { value: entry.recordLengths.map((record_length) => ({ record_length })) },
      form: {},
      adv: {}
    }));
    assert.equal(display.autoValueText('labelFontSize'), `${entry.expected['labels.font_size']} px (auto)`, entry.name);
  }
});
