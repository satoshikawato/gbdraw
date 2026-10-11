// @ts-check
import { uniqueOrthogroupEntries } from './feature-identity.js';

const DEFAULT_RIGHT_DRAWER_TAB = 'features';
export const ALWAYS_AVAILABLE_TABS = new Set(['legend', DEFAULT_RIGHT_DRAWER_TAB]);

/**
 * @param {unknown} value
 * @returns {number}
 */
const normalizedOrthogroupCount = (value) => {
  const count = Number(value);
  return Number.isFinite(count) && count > 0 ? count : 0;
};

/**
 * @param {Record<string, any>} state Shape owned by `state.js`.
 * @returns {number}
 */
export const orthogroupTabContentCountFromState = (state) => (
  uniqueOrthogroupEntries(state.orthogroups?.value).length
  + Number(Boolean(state.similarityAlignmentPlan?.value))
);

/**
 * @param {string} tab
 * @param {number} [orthogroupCount]
 * @returns {boolean}
 */
export const isRightDrawerTabAvailable = (tab, orthogroupCount = 0) => {
  if (ALWAYS_AVAILABLE_TABS.has(tab)) return true;
  return tab === 'orthogroups' && normalizedOrthogroupCount(orthogroupCount) > 0;
};

/**
 * @param {string} tab
 * @param {number} [orthogroupCount]
 * @returns {string}
 */
const resolveRightDrawerTab = (tab, orthogroupCount = 0) => (
  isRightDrawerTabAvailable(tab, orthogroupCount)
    ? tab
    : DEFAULT_RIGHT_DRAWER_TAB
);

/**
 * @param {Record<string, any>} state Shape owned by `state.js`.
 */
export const captureRightDrawerState = (state) => ({
  showRightDrawer: Boolean(state.showRightDrawer.value),
  rightDrawerTab: state.rightDrawerTab.value
});

/**
 * @param {Record<string, any>} state Shape owned by `state.js`.
 * @param {number} [orthogroupCount]
 * @returns {string}
 */
export const reconcileRightDrawerState = (
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

/**
 * @param {Record<string, any>} state Shape owned by `state.js`.
 * @param {string} [tab]
 * @param {number} [orthogroupCount]
 * @returns {string}
 */
export const openRightDrawerState = (
  state,
  tab = state.rightDrawerTab.value,
  orthogroupCount = orthogroupTabContentCountFromState(state)
) => {
  const resolvedTab = resolveRightDrawerTab(tab, orthogroupCount);
  state.rightDrawerTab.value = resolvedTab;
  state.showRightDrawer.value = true;
  return resolvedTab;
};

/**
 * @param {Record<string, any>} state Shape owned by `state.js`.
 */
export const closeRightDrawerState = (state) => {
  state.showRightDrawer.value = false;
};

/**
 * @param {Record<string, any>} state Shape owned by `state.js`.
 */
export const resetRightDrawerState = (state) => {
  state.showRightDrawer.value = false;
  state.rightDrawerTab.value = DEFAULT_RIGHT_DRAWER_TAB;
};

/**
 * @param {Record<string, any>} state Shape owned by `state.js`.
 * @param {Record<string, any> | null | undefined} snapshot From `captureRightDrawerState`.
 * @param {number} [orthogroupCount]
 */
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
