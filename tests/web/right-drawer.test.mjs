import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import test from 'node:test';

import { createRightDrawerController } from '../../gbdraw/web/js/app/right-drawer.js';
import {
  captureRightDrawerState,
  restoreRightDrawerState
} from '../../gbdraw/web/js/services/right-drawer-state.js';
import { settingApplicationNote } from '../../gbdraw/web/js/services/setting-application.js';
import { loadPageVue } from './helpers/vue-runtime.mjs';

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
  // No cache: each read evaluates, as the scans did before OV-921.
  const computed = (getter) => ({ get value() { return getter(); } });
  return { watch, computed, registrations, flush };
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
    watch: harness.watch, computed: harness.computed
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
    watch: harness.watch, computed: harness.computed
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
    watch: harness.watch, computed: harness.computed
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
  const drawer = createRightDrawerController({ state, watch: harness.watch, computed: harness.computed });

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
test('the Legend tab help does not claim a fill color edit it lacks (PV-12)', () => {
  const html = readFileSync(new URL('../../gbdraw/web/index.html', import.meta.url), 'utf8');
  const start = html.indexOf(`v-show="rightDrawerTab === 'legend'"`);
  const end = html.indexOf('v-if="isFeatureDrawerMounted"', start);
  assert.ok(start >= 0 && end > start, 'the Legend tab markup is found');
  const legendTab = html.slice(start, end);
  // The tab's help binds the SETTING_APPLICATION owner's note for the legend editor.
  assert.match(legendTab, /<help-tip\b[^>]*:text="settingApplicationNote\('legend-editor'\)\.text/);
  const help = settingApplicationNote('legend-editor').text;
  assert.match(help, /^Live edit:/);
  assert.match(legendTab, /renameLegendEntry\(|deleteLegendEntry\(|setLegendEntryStrokeColorValue\(/);
  const claimsFillColor = /\bcolou?rs?\b/i.test(help.replace(/\bstroke colou?rs?\b/gi, ''));
  assert.ok(!claimsFillColor || /updateLegendEntryColor\(/.test(legendTab), help);
});

// PV-11 (PD-OI-038): Close and Escape return keyboard focus from inside the
// drawer to the Editor toggle; focus outside the drawer stays where it is.
test('closing the drawer returns focus from inside it to the Editor toggle', () => {
  for (const focusInDrawer of [true, false]) {
    const state = createState({ open: true });
    const harness = createWatchHarness();
    const events = [];
    const drawer = createRightDrawerController({
      state,
      watch: harness.watch, computed: harness.computed,
      focusReturn: {
        isFocusInDrawer: () => focusInDrawer,
        focusToggle: () => events.push(`focus:${state.showRightDrawer.value}`)
      }
    });
    drawer.closeRightDrawer();
    assert.equal(state.showRightDrawer.value, false);
    assert.deepEqual(events, focusInDrawer ? ['focus:false'] : []);
  }
});

// OV-921: the template asks for the Similarity groups tab's availability three
// times per render. The group list is scanned once per change of the groups,
// not per question.
test('the tab availability scans the groups once per change, not per read (OV-921)', () => {
  const { ref: vueRef, computed, watch } = loadPageVue();
  let idReads = 0;
  const counted = (id) => new Proxy({ id, members: [] }, {
    get(target, key, receiver) {
      if (key === 'id') idReads += 1;
      return Reflect.get(target, key, receiver);
    }
  });
  const state = {
    showRightDrawer: vueRef(true),
    rightDrawerTab: vueRef('features'),
    orthogroups: vueRef([counted('og_1'), counted('og_2')]),
    similarityAlignmentPlan: vueRef(null)
  };
  const drawer = createRightDrawerController({ state, watch, computed });
  const afterSetup = idReads;
  for (let render = 0; render < 5; render += 1) {
    for (let read = 0; read < 3; read += 1) assert.equal(drawer.isRightDrawerTabAvailable('orthogroups'), true);
  }
  assert.equal(idReads, afterSetup, 'reads after setup reuse one scan');
  state.orthogroups.value = [];
  assert.equal(drawer.isRightDrawerTabAvailable('orthogroups'), false);
  assert.equal(state.rightDrawerTab.value, 'features');
});
