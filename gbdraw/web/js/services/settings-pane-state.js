// @ts-check
// The display of the settings pane and its transitions. The pane shows the
// stacked settings, or an icon rail with at most one docked entry: the
// sections of that entry beside the rail. Nothing is docked in the stacked
// display, and a docked entry is always one the mode has. One section may
// leave the pane for a popup beside the Result; a section is in one place, so
// docking that section closes the popup, and so does a mode change. The popup's
// section and its place are written here alone.
import { SETTINGS_SECTIONS, settingsEntriesOf, settingsEntryIn } from './settings-sections.js';

/** @import { PopupPlacement } from '../utils/popup-bounds.js' */
/** @import { SettingsSectionKey } from './settings-sections.js' */

/** @typedef {'stacked' | 'rail'} SettingsPaneDisplay */
/** @typedef {'circular' | 'linear'} DiagramMode */
/** @typedef {PopupPlacement & { section: SettingsSectionKey }} SettingsPopup The section in the popup, and where the popup is. */
/**
 * The three refs of `state.js` that hold the display.
 * @typedef {object} SettingsPaneRefs
 * @property {{ value: SettingsPaneDisplay }} settingsPaneDisplay
 * @property {{ value: SettingsSectionKey | null }} settingsPaneDocked
 * @property {{ value: SettingsPopup | null }} settingsPopup
 */

// Closes the popup: its section is in the pane again.
/** @param {SettingsPaneRefs} state */
export const returnSettingsSection = (state) => {
  state.settingsPopup.value = null;
};

/**
 * Shows a section of the shown mode in the popup, in place of the one there.
 * @param {SettingsPaneRefs} state
 * @param {SettingsSectionKey} section
 * @param {PopupPlacement} placement
 */
export const floatSettingsSection = (state, section, placement) => {
  state.settingsPopup.value = { ...placement, section };
};

/**
 * Moves the popup (a drag, its rendered size, a window that changed size); its
 * section stays. Nothing moves while no section is in it.
 * @param {SettingsPaneRefs} state
 * @param {PopupPlacement} placement
 */
export const moveSettingsPopup = (state, placement) => {
  if (state.settingsPopup.value) Object.assign(state.settingsPopup.value, placement);
};

/** @param {SettingsPaneRefs} state @param {SettingsSectionKey | null} docked The entry docked now. */
const dock = (state, docked) => {
  state.settingsPaneDocked.value = docked;
  if (docked && state.settingsPopup.value?.section === docked) returnSettingsSection(state);
  return docked;
};

/** @type {readonly SettingsPaneDisplay[]} */
export const SETTINGS_PANE_DISPLAYS = Object.freeze(['stacked', 'rail']);
export const DEFAULT_SETTINGS_PANE_DISPLAY = SETTINGS_PANE_DISPLAYS[0];

// The rail entry that shows a section in a mode, or the first entry (Input)
// when the mode has no entry for it.
/** @param {DiagramMode} mode @param {SettingsSectionKey} section */
const entryIn = (mode, section) => settingsEntryIn(mode, section) ?? settingsEntriesOf(mode, 'rail')[0].key;

/**
 * Switches to the rail with the entry of `section` docked.
 * @param {SettingsPaneRefs} state
 * @param {DiagramMode} mode
 * @param {SettingsSectionKey} [section]
 */
export const showSettingsRail = (state, mode, section = SETTINGS_SECTIONS[0].key) => {
  state.settingsPaneDisplay.value = 'rail';
  return dock(state, entryIn(mode, section));
};

/** @param {SettingsPaneRefs} state */
export const showSettingsStacked = (state) => {
  state.settingsPaneDisplay.value = 'stacked';
  state.settingsPaneDocked.value = null;
};

/**
 * Docks the rail entry that shows `section`. The stacked display docks nothing.
 * @param {SettingsPaneRefs} state
 * @param {DiagramMode} mode
 * @param {SettingsSectionKey} section
 */
export const dockSettingsSection = (state, mode, section) => {
  return state.settingsPaneDisplay.value === 'rail' ? dock(state, entryIn(mode, section)) : state.settingsPaneDocked.value;
};

/** @param {SettingsPaneRefs} state */
export const closeDockedSettingsEntry = (state) => {
  state.settingsPaneDocked.value = null;
};

/**
 * A mode change: the popup closes, and a docked entry the mode lacks gives way
 * to the mode's entry of the same title, or to Input when the mode has none.
 * @param {SettingsPaneRefs} state
 * @param {DiagramMode} mode
 */
export const reconcileSettingsPaneState = (state, mode) => {
  returnSettingsSection(state);
  const docked = state.settingsPaneDocked.value;
  if (docked && entryIn(mode, docked) !== docked) state.settingsPaneDocked.value = entryIn(mode, docked);
  return state.settingsPaneDocked.value;
};
