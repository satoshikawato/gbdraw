// @ts-check
import { uniqueOrthogroupEntries } from './feature-identity.js';
/** @import { WebState } from '../state.js' */

/** @typedef {Pick<WebState, 'showRightDrawer' | 'rightDrawerTab'>} RightDrawerMembers Whether the drawer is open, and its tab. */
/** @typedef {Pick<WebState, 'orthogroups' | 'similarityAlignmentPlan'>} OrthogroupTabContent What the Orthogroups tab shows. */
/** @typedef {RightDrawerMembers & OrthogroupTabContent} RightDrawerState */
/** @typedef {'legend' | 'features' | 'orthogroups'} RightDrawerTab */

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
 * @param {OrthogroupTabContent} state
 * @returns {number}
 */
export const orthogroupTabContentCountFromState = (state) => (
  uniqueOrthogroupEntries(state.orthogroups.value).length
  + Number(Boolean(state.similarityAlignmentPlan.value))
);

/**
 * @param {string} tab
 * @param {number} [orthogroupCount]
 * @returns {tab is RightDrawerTab}
 */
export const isRightDrawerTabAvailable = (tab, orthogroupCount = 0) => {
  if (ALWAYS_AVAILABLE_TABS.has(tab)) return true;
  return tab === 'orthogroups' && normalizedOrthogroupCount(orthogroupCount) > 0;
};

/**
 * @param {string} tab
 * @param {number} [orthogroupCount]
 * @returns {RightDrawerTab}
 */
const resolveRightDrawerTab = (tab, orthogroupCount = 0) => (
  isRightDrawerTabAvailable(tab, orthogroupCount)
    ? tab
    : DEFAULT_RIGHT_DRAWER_TAB
);

/** @param {RightDrawerMembers} state */
export const captureRightDrawerState = (state) => ({
  showRightDrawer: state.showRightDrawer.value,
  rightDrawerTab: state.rightDrawerTab.value
});

/**
 * @param {RightDrawerState} state
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
 * @param {RightDrawerState} state
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

/** @param {RightDrawerMembers} state */
export const closeRightDrawerState = (state) => {
  state.showRightDrawer.value = false;
};

/** @param {RightDrawerMembers} state */
export const resetRightDrawerState = (state) => {
  state.showRightDrawer.value = false;
  state.rightDrawerTab.value = DEFAULT_RIGHT_DRAWER_TAB;
};

/**
 * @param {RightDrawerState} state
 * @param {ReturnType<typeof captureRightDrawerState>} snapshot
 * @param {number} [orthogroupCount]
 */
export const restoreRightDrawerState = (
  state,
  snapshot,
  orthogroupCount = orthogroupTabContentCountFromState(state)
) => {
  state.rightDrawerTab.value = resolveRightDrawerTab(
    snapshot.rightDrawerTab,
    orthogroupCount
  );
  state.showRightDrawer.value = snapshot.showRightDrawer;
};
