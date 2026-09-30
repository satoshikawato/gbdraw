import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import test from 'node:test';

import {
  captureRightDrawerState,
  createRightDrawerController,
  restoreRightDrawerState
} from '../../gbdraw/web/js/app/right-drawer.js';
import { assertKnownDefect } from './helpers/known-defect.mjs';

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

const createState = ({ open = false, tab = 'features', groups = [], plan = null } = {}) => ({
  showRightDrawer: ref(open),
  rightDrawerTab: ref(tab),
  orthogroups: ref(groups),
  similarityAlignmentPlan: ref(plan)
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

test('an active alignment keeps its inspector and Reset tab reachable without group rows', () => {
  const state = createState({ tab: 'orthogroups', plan: { schema: 1 } });
  const harness = createWatchHarness();
  const drawer = createRightDrawerController({ state, watch: harness.watch });

  assert.equal(drawer.isRightDrawerTabAvailable('orthogroups'), true);
  assert.equal(drawer.openRightDrawerTab('orthogroups'), 'orthogroups');
  state.orthogroups.value = [];
  harness.flush();
  assert.equal(state.rightDrawerTab.value, 'orthogroups');

  state.similarityAlignmentPlan.value = null;
  harness.flush();
  assert.equal(drawer.isRightDrawerTabAvailable('orthogroups'), false);
  assert.equal(state.rightDrawerTab.value, 'features');
});

// PV-12 (OIC-024): the Legend tab help names only edits that the tab offers.
// The entry swatch displays the fill color; no control in the tab changes it.
test('the Legend tab help does not claim a fill color edit it lacks (PV-12 known defect)', async () => {
  const html = readFileSync(new URL('../../gbdraw/web/index.html', import.meta.url), 'utf8');
  const start = html.indexOf(`v-show="rightDrawerTab === 'legend'"`);
  const end = html.indexOf('v-if="isFeatureDrawerMounted"', start);
  assert.ok(start >= 0 && end > start, 'the Legend tab markup is found');
  const legendTab = html.slice(start, end);
  const help = legendTab.match(/<p\b[^>]*>\s*(Live edit:[^<]*)<\/p>/)?.[1] || '';
  assert.match(help, /^Live edit:/);
  assert.match(legendTab, /renameLegendEntry\(|deleteLegendEntry\(|setLegendEntryStrokeColorValue\(/);
  await assertKnownDefect('PV-12', () => {
    const claimsFillColor = /\bcolou?rs?\b/i.test(help.replace(/\bstroke colou?rs?\b/gi, ''));
    assert.ok(!claimsFillColor || /updateLegendEntryColor\(/.test(legendTab), help);
  });
});
