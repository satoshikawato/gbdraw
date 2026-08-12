import assert from 'node:assert/strict';
import test from 'node:test';

import {
  captureRightDrawerState,
  createRightDrawerController,
  restoreRightDrawerState
} from '../../gbdraw/web/js/app/right-drawer.js';

const ref = (value) => ({ value });

const createWatchHarness = () => {
  const registrations = [];
  const watch = (source, callback, options = {}) => {
    const registration = { source, callback, options };
    registrations.push(registration);
    if (options.immediate) callback(source(), undefined);
    return () => {};
  };
  const flush = () => registrations.forEach(({ source, callback }) => callback(source()));
  return { watch, registrations, flush };
};

const createState = ({ open = false, tab = 'features', groups = [] } = {}) => ({
  showRightDrawer: ref(open),
  rightDrawerTab: ref(tab),
  orthogroups: ref(groups)
});

test('unavailable and unknown tabs resolve to Features instead of becoming no-ops', () => {
  const state = createState();
  const harness = createWatchHarness();
  const drawer = createRightDrawerController({
    state,
    watch: harness.watch
  });

  assert.equal(harness.registrations[0].options.flush, 'sync');
  assert.equal(drawer.openRightDrawerTab('orthogroups'), 'features');
  assert.equal(state.showRightDrawer.value, true);
  assert.equal(state.rightDrawerTab.value, 'features');

  drawer.closeRightDrawer();
  assert.equal(drawer.openRightDrawerTab('not-a-tab'), 'features');
  assert.equal(state.showRightDrawer.value, true);
  assert.equal(state.rightDrawerTab.value, 'features');
});

test('capability loss normalizes the selected tab synchronously while open or closed', () => {
  const state = createState({
    tab: 'orthogroups',
    groups: [{ id: 'group-a', members: [] }]
  });
  const harness = createWatchHarness();
  const drawer = createRightDrawerController({
    state,
    watch: harness.watch
  });

  drawer.openRightDrawerTab('orthogroups');
  state.orthogroups.value = [];
  harness.flush();
  assert.equal(state.showRightDrawer.value, true);
  assert.equal(state.rightDrawerTab.value, 'features');

  state.orthogroups.value = [{ id: 'group-a', members: [] }];
  harness.flush();
  assert.equal(state.rightDrawerTab.value, 'features');

  drawer.openRightDrawerTab('orthogroups');
  drawer.closeRightDrawer();
  state.orthogroups.value = [];
  harness.flush();
  assert.equal(state.showRightDrawer.value, false);
  assert.equal(state.rightDrawerTab.value, 'features');

  drawer.toggleRightDrawer();
  assert.equal(state.showRightDrawer.value, true);
  assert.equal(state.rightDrawerTab.value, 'features');
});

test('close preserves a valid selection, reset clears it, and rollback restores then validates it', () => {
  const state = createState({
    open: true,
    tab: 'orthogroups',
    groups: [{ id: 'group-a', members: [] }]
  });
  const harness = createWatchHarness();
  const drawer = createRightDrawerController({
    state,
    watch: harness.watch
  });
  const snapshot = captureRightDrawerState(state);

  drawer.closeRightDrawer();
  assert.equal(state.showRightDrawer.value, false);
  assert.equal(state.rightDrawerTab.value, 'orthogroups');

  drawer.resetRightDrawer();
  assert.equal(state.showRightDrawer.value, false);
  assert.equal(state.rightDrawerTab.value, 'features');

  restoreRightDrawerState(state, snapshot);
  assert.deepEqual(captureRightDrawerState(state), snapshot);

  state.orthogroups.value = [];
  restoreRightDrawerState(state, snapshot);
  assert.deepEqual(captureRightDrawerState(state), {
    showRightDrawer: true,
    rightDrawerTab: 'features'
  });
});
