// @ts-check
// The display of the settings pane and its transitions. The pane shows the
// stacked settings, or an icon rail with at most one docked entry: the
// sections of that entry beside the rail. Nothing is docked in the stacked
// display, and a docked entry is always one the mode has.
import { SETTINGS_SECTIONS, settingsEntriesOf, settingsEntryIn } from './settings-sections.js';

/** @import { SettingsSectionKey } from './settings-sections.js' */

/** @typedef {'stacked' | 'rail'} SettingsPaneDisplay */
/** @typedef {'circular' | 'linear'} DiagramMode */
/**
 * The two refs of `state.js` that hold the display.
 * @typedef {object} SettingsPaneRefs
 * @property {{ value: SettingsPaneDisplay }} settingsPaneDisplay
 * @property {{ value: SettingsSectionKey | null }} settingsPaneDocked
 */

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
  state.settingsPaneDocked.value = entryIn(mode, section);
  return state.settingsPaneDocked.value;
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
  if (state.settingsPaneDisplay.value === 'rail') state.settingsPaneDocked.value = entryIn(mode, section);
  return state.settingsPaneDocked.value;
};

/** @param {SettingsPaneRefs} state */
export const closeDockedSettingsEntry = (state) => {
  state.settingsPaneDocked.value = null;
};

/**
 * A mode change: a docked entry the mode lacks gives way to the mode's entry
 * of the same title, or to Input when the mode has none.
 * @param {SettingsPaneRefs} state
 * @param {DiagramMode} mode
 */
export const reconcileSettingsPaneState = (state, mode) => {
  const docked = state.settingsPaneDocked.value;
  if (docked && entryIn(mode, docked) !== docked) state.settingsPaneDocked.value = entryIn(mode, docked);
  return state.settingsPaneDocked.value;
};
