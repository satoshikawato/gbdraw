// @ts-check
// What the app keeps in the browser between visits: the display of the
// settings pane, and nothing else. This module is the only one that touches
// browser storage (tests/web/browser-preferences.test.mjs). Storage that is
// absent, blocked, full, or holds another value never reaches the user: a
// read gives the default display and a write is dropped.
import { DEFAULT_SETTINGS_PANE_DISPLAY, SETTINGS_PANE_DISPLAYS } from './settings-pane-state.js';

/** @import { SettingsPaneDisplay } from './settings-pane-state.js' */

const SETTINGS_PANE_DISPLAY_KEY = 'gbdraw.settings-pane-display';

// Reading the property itself can throw when site data is blocked.
/** @returns {Pick<Storage, 'getItem' | 'setItem'> | undefined} */
const browserStorage = () => window.localStorage;

/** @param {typeof browserStorage} [storage] */
export const readSettingsPaneDisplay = (storage = browserStorage) => {
  try {
    const stored = storage()?.getItem(SETTINGS_PANE_DISPLAY_KEY);
    return SETTINGS_PANE_DISPLAYS.find((display) => display === stored) ?? DEFAULT_SETTINGS_PANE_DISPLAY;
  } catch {
    return DEFAULT_SETTINGS_PANE_DISPLAY;
  }
};

/** @param {SettingsPaneDisplay} display @param {typeof browserStorage} [storage] */
export const writeSettingsPaneDisplay = (display, storage = browserStorage) => {
  try {
    storage()?.setItem(SETTINGS_PANE_DISPLAY_KEY, display);
  } catch {
    // The display still switches for this page.
  }
};
