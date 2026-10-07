// @ts-check
import {
  ALWAYS_AVAILABLE_TABS,
  closeRightDrawerState,
  isRightDrawerTabAvailable,
  openRightDrawerState,
  orthogroupTabContentCountFromState,
  reconcileRightDrawerState,
  resetRightDrawerState
} from '../services/right-drawer-state.js';

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
    if (returnFocus) focusReturn?.focusToggle();
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
