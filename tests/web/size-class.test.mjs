import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import test from 'node:test';

import { createAutoValueDisplay } from '../../gbdraw/web/js/app/auto-value-display.js';
import { sizeClassForRecordLengths } from '../../gbdraw/web/js/services/track-slot-display.js';
import { withDrawings as withFixtureDrawings } from './helpers/drawing-state.mjs';

// Every drawing holds the settings its Session preserved (state.js); a case
// that keeps a boundary passes its own.
const withDrawings = (/** @type {Record<string, unknown>} */ fixture) => (
  withFixtureDrawings({ unmanagedConfigOverrides: {}, ...fixture })
);

test('the Web size class uses the packaged short/long boundary of each mode', () => {
  // The boundaries reach the Web through mode-profiles.generated.js (labels.length_threshold).
  for (const [mode, boundary] of /** @type {const} */ ([['circular', 75_000], ['linear', 100_000]])) {
    const drawing = { unmanagedConfigOverrides: {} };
    assert.equal(sizeClassForRecordLengths(mode, [3000, boundary - 1], drawing), 'short', mode);
    assert.equal(sizeClassForRecordLengths(mode, [3000, boundary], drawing), 'long', mode);
    assert.equal(sizeClassForRecordLengths(mode, [], drawing), null, mode);
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

// A Session that stores its whole configuration keeps the boundary it was saved
// with (labels.length_threshold.<mode> among the drawing's preserved settings),
// and Generate draws with it. The Auto notes and the track preview read the same
// boundary: 60,000 bp is short by both packaged boundaries and long by a kept 50 kb.
test('a drawing that keeps a stored boundary shows the Auto values Generate draws', () => {
  const lengths = [60_000];
  const drawing = (/** @type {Record<string, unknown>} */ unmanagedConfigOverrides) => ({
    form: {}, adv: {}, unmanagedConfigOverrides
  });
  const display = (/** @type {'circular' | 'linear'} */ mode, /** @type {Record<string, unknown>} */ overrides) => (
    createAutoValueDisplay(withDrawings({
      mode: { value: mode },
      circularRecordList: { value: lengths.map((record_length) => ({ record_length })) },
      ...drawing(overrides)
    }), { linearRecordLengths: () => lengths })
  );

  for (const mode of /** @type {const} */ (['circular', 'linear'])) {
    const kept = { [`labels.length_threshold.${mode}`]: 50_000 };
    assert.equal(sizeClassForRecordLengths(mode, lengths, drawing({})), 'short', mode);
    assert.equal(sizeClassForRecordLengths(mode, lengths, drawing(kept)), 'long', mode);
    // The other mode's boundary is not this drawing's.
    const other = mode === 'circular' ? 'linear' : 'circular';
    assert.equal(
      sizeClassForRecordLengths(mode, lengths, drawing({ [`labels.length_threshold.${other}`]: 50_000 })),
      'short',
      mode
    );
  }

  const circular = { packaged: display('circular', {}), kept: display('circular', { 'labels.length_threshold.circular': 50_000 }) };
  assert.equal(circular.packaged.autoValueText('labelFontSize'), '14 px (auto)');
  assert.equal(circular.kept.autoValueText('labelFontSize'), '8 px (auto)');
  assert.equal(circular.kept.autoValueText('axisStrokeWidth'), '1 px (auto)');
  assert.equal(circular.kept.autoValueText('circularFeatureWidth'), '18.5 px (auto)');
  assert.equal(circular.packaged.autoValueText('circularFeatureWidth'), '37 px (auto)');

  const linear = { packaged: display('linear', {}), kept: display('linear', { 'labels.length_threshold.linear': 50_000 }) };
  assert.equal(linear.packaged.autoValueText('labelFontSize'), '24 px (auto)');
  assert.equal(linear.kept.autoValueText('labelFontSize'), '5 px (auto)');
  assert.equal(linear.kept.autoValueText('linearFeatureHeight'), '20 px (auto)');
  assert.equal(linear.kept.autoValueText('legendFontSize'), '16 px (auto)');
});
