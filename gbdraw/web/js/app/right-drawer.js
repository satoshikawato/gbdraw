import { uniqueOrthogroupEntries } from '../services/feature-identity.js';

const DEFAULT_RIGHT_DRAWER_TAB = 'features';
const ALWAYS_AVAILABLE_TABS = new Set(['legend', DEFAULT_RIGHT_DRAWER_TAB]);

const normalizedOrthogroupCount = (value) => {
  const count = Number(value);
  return Number.isFinite(count) && count > 0 ? count : 0;
};

const orthogroupCountFromState = (state) => uniqueOrthogroupEntries(
  state.orthogroups?.value
).length;

const isRightDrawerTabAvailable = (tab, orthogroupCount = 0) => {
  if (ALWAYS_AVAILABLE_TABS.has(tab)) return true;
  return tab === 'orthogroups' && normalizedOrthogroupCount(orthogroupCount) > 0;
};

const resolveRightDrawerTab = (tab, orthogroupCount = 0) => (
  isRightDrawerTabAvailable(tab, orthogroupCount)
    ? tab
    : DEFAULT_RIGHT_DRAWER_TAB
);

export const captureRightDrawerState = (state) => ({
  showRightDrawer: Boolean(state.showRightDrawer.value),
  rightDrawerTab: state.rightDrawerTab.value
});

const reconcileRightDrawerState = (
  state,
  orthogroupCount = orthogroupCountFromState(state)
) => {
  const resolvedTab = resolveRightDrawerTab(
    state.rightDrawerTab.value,
    orthogroupCount
  );
  if (state.rightDrawerTab.value !== resolvedTab) {
    state.rightDrawerTab.value = resolvedTab;
  }
  return resolvedTab;
};

const openRightDrawerState = (
  state,
  tab = state.rightDrawerTab.value,
  orthogroupCount = orthogroupCountFromState(state)
) => {
  const resolvedTab = resolveRightDrawerTab(tab, orthogroupCount);
  state.rightDrawerTab.value = resolvedTab;
  state.showRightDrawer.value = true;
  return resolvedTab;
};

const closeRightDrawerState = (state) => {
  state.showRightDrawer.value = false;
};

export const resetRightDrawerState = (state) => {
  state.showRightDrawer.value = false;
  state.rightDrawerTab.value = DEFAULT_RIGHT_DRAWER_TAB;
};

export const restoreRightDrawerState = (
  state,
  snapshot,
  orthogroupCount = orthogroupCountFromState(state)
) => {
  state.rightDrawerTab.value = resolveRightDrawerTab(
    snapshot?.rightDrawerTab,
    orthogroupCount
  );
  state.showRightDrawer.value = Boolean(snapshot?.showRightDrawer);
};

export const createRightDrawerController = ({ state, watch }) => {
  const currentOrthogroupCount = () => orthogroupCountFromState(state);
  const isTabAvailable = (tab) => isRightDrawerTabAvailable(
    tab,
    currentOrthogroupCount()
  );
  const reconcile = () => reconcileRightDrawerState(
    state,
    currentOrthogroupCount()
  );
  const openRightDrawerTab = (tab = state.rightDrawerTab.value) => openRightDrawerState(
    state,
    tab,
    currentOrthogroupCount()
  );
  const closeRightDrawer = () => closeRightDrawerState(state);
  const resetRightDrawer = () => resetRightDrawerState(state);
  const toggleRightDrawer = () => {
    if (state.showRightDrawer.value) {
      closeRightDrawer();
      return;
    }
    openRightDrawerTab(state.rightDrawerTab.value);
  };

  watch(
    currentOrthogroupCount,
    reconcile,
    { flush: 'sync', immediate: true }
  );

  return {
    isRightDrawerTabAvailable: isTabAvailable,
    openRightDrawerTab,
    toggleRightDrawer,
    closeRightDrawer,
    resetRightDrawer
  };
};
