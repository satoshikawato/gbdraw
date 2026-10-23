// @ts-check
// What the app keeps in the browser between visits: the display of the
// settings pane (PD-OI-093) and the UI language (PD-OI-106), and nothing else.
// This module is the only one that touches browser storage
// (tests/web/browser-preferences.test.mjs). Storage that is absent, blocked,
// full, or holds another value never reaches the user: a read gives the
// default and a write is dropped.
import { DEFAULT_SETTINGS_PANE_DISPLAY, SETTINGS_PANE_DISPLAYS } from './settings-pane-state.js';
import { SOURCE_LOCALE, uiLocaleOf } from '../utils/i18n.js';

/** @import { SettingsPaneDisplay } from './settings-pane-state.js' */
/** @import { UiLocale } from '../utils/i18n.js' */

const SETTINGS_PANE_DISPLAY_KEY = 'gbdraw.settings-pane-display';
const UI_LANGUAGE_KEY = 'gbdraw.ui-language';

// Reading the property itself can throw when site data is blocked.
/** @returns {Pick<Storage, 'getItem' | 'setItem'> | undefined} */
const browserStorage = () => window.localStorage;

/** @param {string} key @param {typeof browserStorage} storage */
const readStored = (key, storage) => {
  try {
    return storage()?.getItem(key) ?? null;
  } catch {
    return null;
  }
};

/** @param {string} key @param {string} value @param {typeof browserStorage} storage */
const writeStored = (key, value, storage) => {
  try {
    storage()?.setItem(key, value);
  } catch {
    // The choice still applies to this page.
  }
};

/** @param {typeof browserStorage} [storage] @returns {SettingsPaneDisplay} */
export const readSettingsPaneDisplay = (storage = browserStorage) => {
  const stored = readStored(SETTINGS_PANE_DISPLAY_KEY, storage);
  return SETTINGS_PANE_DISPLAYS.find((display) => display === stored) ?? DEFAULT_SETTINGS_PANE_DISPLAY;
};

/** @param {SettingsPaneDisplay} display @param {typeof browserStorage} [storage] */
export const writeSettingsPaneDisplay = (display, storage = browserStorage) => {
  writeStored(SETTINGS_PANE_DISPLAY_KEY, display, storage);
};

/**
 * The UI language this browser chose. A first visit, or anything else
 * stored, gives English: the browser's own language never chooses.
 * @param {typeof browserStorage} [storage]
 * @returns {UiLocale}
 */
export const readUiLanguage = (storage = browserStorage) => (
  uiLocaleOf(readStored(UI_LANGUAGE_KEY, storage)) ?? SOURCE_LOCALE
);

/** @param {UiLocale} language @param {typeof browserStorage} [storage] */
export const writeUiLanguage = (language, storage = browserStorage) => {
  writeStored(UI_LANGUAGE_KEY, language, storage);
};
