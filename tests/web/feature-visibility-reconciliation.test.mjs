import assert from 'node:assert/strict';
import test from 'node:test';
import { reconcileFeatureVisibilityOverrides } from '../../gbdraw/web/js/app/feature-visibility.js';

test('source replacement discards absent feature intent and cannot resurrect it on return', () => {
  const overrides = { fa: 'off', fb: 'exclude_matching', shared: 'on' };
  const sourceB = [{ svg_id: 'fb' }, { svg_id: 'shared' }];
  assert.equal(reconcileFeatureVisibilityOverrides(overrides, sourceB), 1);
  assert.deepEqual(overrides, { fb: 'exclude_matching', shared: 'on' });
  const restored = JSON.parse(JSON.stringify(overrides));
  reconcileFeatureVisibilityOverrides(restored, [{ svg_id: 'fa' }, { svg_id: 'shared' }]);
  assert.deepEqual(restored, { shared: 'on' });
});

test('hidden, cropped, and unchanged semantic targets survive without a rendered path', () => {
  const overrides = { hidden_record_2: 'off', croppedDisplay: 'on' };
  const biological = [{ svg_id: 'hidden' }, { svg_id: 'original' }];
  const previous = [{ svg_id: 'croppedDisplay', stable_svg_id: 'original' }];
  assert.equal(reconcileFeatureVisibilityOverrides(overrides, biological, previous), 0);
  assert.deepEqual(overrides, { hidden_record_2: 'off', croppedDisplay: 'on' });
});

test('a different biological target cannot inherit an old rendered ID', () => {
  const overrides = { sameDisplay: 'off' };
  reconcileFeatureVisibilityOverrides(overrides, [{ svg_id: 'newBiology' }], [
    { svg_id: 'sameDisplay', stable_svg_id: 'oldBiology' }
  ]);
  assert.deepEqual(overrides, {});
});

test('hidden rendered instances retain the stable selector owned by the visibility cache', () => {
  const overrides = { opaqueInstance: 'off' };
  const cache = { opaqueInstance: { qualifier: 'hash', value: 'shared' } };
  reconcileFeatureVisibilityOverrides(overrides, [{ svg_id: 'shared' }], [], cache);
  assert.deepEqual(overrides, { opaqueInstance: 'off' });
  reconcileFeatureVisibilityOverrides(overrides, [{ svg_id: 'different' }], [], cache);
  assert.deepEqual(overrides, {});
});

test('no individual intent requires no source traversal', () => {
  assert.equal(reconcileFeatureVisibilityOverrides({}, new Proxy([], {
    get() { throw new Error('unnecessary catalog traversal'); }
  })), 0);
});
