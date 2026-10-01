import assert from 'node:assert/strict';
import { createRequire } from 'node:module';
import test from 'node:test';

// G-C (Web GUI audit 2026-09-30): the browser spec
// non-edit-state-preservation.playwright.spec.js relies on this diff to report
// every silent change to user-owned state, so the diff itself must not hide one.
const require = createRequire(import.meta.url);
const { diffUserOwnedState, expectNoSilentStateChange } = require('./helpers/app-lifecycle.cjs');

test('the user-owned state diff reports changed, added, removed and resized values', () => {
  const before = {
    config: { adv: { multi_record_positions: [{ selector: '#2' }, { selector: '#1' }] } },
    features: { labelTextFeatureOverrides: { f1: 'EDITED' } },
    ui: { canvasPadding: { right: 150 } }
  };
  const after = {
    config: { adv: { multi_record_positions: [{ selector: '#1' }, { selector: '#2' }] } },
    features: { labelTextFeatureOverrides: {}, labelVisibilityOverrides: { f2: 'off' } },
    ui: { canvasPadding: { right: 0 }, trackSlots: [1, 2, 3] }
  };
  assert.deepEqual(diffUserOwnedState(before, structuredClone(before)), []);
  assert.deepEqual(diffUserOwnedState(before, after).map(({ path }) => path), [
    'config.adv.multi_record_positions[0].selector',
    'config.adv.multi_record_positions[1].selector',
    'features.labelTextFeatureOverrides.f1',
    'features.labelVisibilityOverrides',
    'ui.canvasPadding.right',
    'ui.trackSlots'
  ]);
  assert.deepEqual(diffUserOwnedState({ rows: [1] }, { rows: [1, 2] }), [
    { path: 'rows', before: [1], after: [1, 2] }
  ]);
});

test('only named owned paths or a recorded History step excuse a change', () => {
  const observation = {
    changes: [{ path: 'ui.canvasPadding.right', before: 150, after: 0 }],
    undoDelta: 0
  };
  assert.deepEqual(expectNoSilentStateChange(observation, {
    label: 'owned', allowedPaths: ['ui.canvasPadding']
  }), []);
  assert.throws(
    () => expectNoSilentStateChange(observation, { label: 'silent', allowedPaths: ['ui.canvasPad'] }),
    /silent silently changed user-owned state/
  );
  assert.throws(
    () => expectNoSilentStateChange(observation, { label: 'unrecorded', allowRecorded: true }),
    /unrecorded silently changed user-owned state/
  );
  assert.equal(
    expectNoSilentStateChange({ ...observation, undoDelta: 1 }, { label: 'recorded', allowRecorded: true }).length,
    1
  );
});
