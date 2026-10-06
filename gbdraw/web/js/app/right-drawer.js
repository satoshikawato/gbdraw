// @ts-check
import { uniqueOrthogroupEntries } from '../services/feature-identity.js';

const DEFAULT_RIGHT_DRAWER_TAB = 'features';
const ALWAYS_AVAILABLE_TABS = new Set(['legend', DEFAULT_RIGHT_DRAWER_TAB]);

const normalizedOrthogroupCount = (value) => {
  const count = Number(value);
  return Number.isFinite(count) && count > 0 ? count : 0;
};

const orthogroupTabContentCountFromState = (state) => (
  uniqueOrthogroupEntries(state.orthogroups?.value).length
  + Number(Boolean(state.similarityAlignmentPlan?.value))
);

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
  orthogroupCount = orthogroupTabContentCountFromState(state)
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
  orthogroupCount = orthogroupTabContentCountFromState(state)
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
  orthogroupCount = orthogroupTabContentCountFromState(state)
) => {
  state.rightDrawerTab.value = resolveRightDrawerTab(
    snapshot?.rightDrawerTab,
    orthogroupCount
  );
  state.showRightDrawer.value = Boolean(snapshot?.showRightDrawer);
};

/**
 * @typedef {Object} RightDrawerFocusReturn
 * @property {() => boolean} isFocusInDrawer
 * @property {() => void} focusToggle
 */

/**
 * @typedef {Object} RightDrawerControllerOptions
 * @property {Record<string, any>} state Shape owned by state.js.
 * @property {(source: () => unknown, callback: () => void, options?: Record<string, any>) => unknown} watch
 * @property {() => string} [getOpenDisabledReason]
 * @property {() => void} [onClose]
 * @property {RightDrawerFocusReturn | null} [focusReturn]
 */

/** @param {RightDrawerControllerOptions} options */
export const createRightDrawerController = ({
  state,
  watch,
  getOpenDisabledReason = () => '',
  onClose = () => {},
  focusReturn = null
}) => {
  const currentOrthogroupCount = () => orthogroupTabContentCountFromState(state);
  const isTabAvailable = (tab) => ALWAYS_AVAILABLE_TABS.has(tab) || isRightDrawerTabAvailable(
    tab,
    currentOrthogroupCount()
  );
  const reconcile = () => reconcileRightDrawerState(
    state,
    currentOrthogroupCount()
  );
  const openRightDrawerTab = (tab = state.rightDrawerTab.value) => {
    if (getOpenDisabledReason()) return false;
    return openRightDrawerState(state, tab, currentOrthogroupCount());
  };
  // Closing returns keyboard focus from inside the drawer to the Editor toggle
  // (PD-OI-038); focus elsewhere stays where it is.
  const closeRightDrawer = () => {
    onClose();
    const returnFocus = Boolean(focusReturn?.isFocusInDrawer?.());
    closeRightDrawerState(state);
    if (returnFocus) focusReturn.focusToggle();
  };
  const resetRightDrawer = () => {
    onClose();
    resetRightDrawerState(state);
  };
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
