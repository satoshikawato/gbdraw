// @ts-check
// The header Language menu (PD-OI-106): the language of the Web UI text. A
// choice applies at once, without a reload, and is remembered per browser
// (services/browser-preferences.js). It changes no Result, setting, or
// History: the menu is `data-history-ignore`, and nothing here writes a
// drawing. The Session, the figure, and every export stay as they are.
// Beside the menu, a link opens the GitHub translation feedback form with the
// shown language and the gbdraw version filled in, and nothing else.
import { GBDRAW_VERSION } from '../config.js';
import { readUiLanguage, writeUiLanguage } from '../services/browser-preferences.js';
import { setUiLocale, uiLocale, uiLocaleOf } from '../utils/i18n.js';

/** @import { UiLocale } from '../utils/i18n.js' */

/**
 * Shows the UI in `language`. A catalog that fails to load leaves the shown
 * language as it was.
 * @param {UiLocale} language
 * @returns {Promise<boolean>} whether the UI now shows `language`
 */
const showUiLanguage = async (language) => {
  try {
    return await setUiLocale(language);
  } catch (error) {
    console.warn(`The ${language} UI text could not be loaded; the UI stays in ${uiLocale.value}.`, error);
    return false;
  }
};

/**
 * Shows the language this browser chose, before the app mounts. English (a
 * first visit) needs no catalog.
 * @returns {Promise<void>}
 */
export const showStoredUiLanguage = async () => {
  const stored = readUiLanguage();
  if (stored !== uiLocale.value) await showUiLanguage(stored);
};

/**
 * The menu's change: shows and remembers the chosen language. When it cannot
 * be shown, the menu shows the language the UI still has.
 * @param {Event} event the `change` of the Language `<select>`
 * @returns {Promise<void>}
 */
export const chooseUiLanguage = async (event) => {
  const menu = /** @type {HTMLSelectElement} */ (event.target);
  const language = uiLocaleOf(menu.value);
  if (language && await showUiLanguage(language)) {
    writeUiLanguage(language);
    return;
  }
  menu.value = uiLocale.value;
};

const TRANSLATION_FEEDBACK_FORM = 'https://github.com/satoshikawato/gbdraw/issues/new?template=translation-feedback.yml';

/**
 * The link of the translation feedback form (`.github/ISSUE_TEMPLATE/translation-feedback.yml`)
 * with its `language` and `version` fields filled in and the language tag in
 * the issue title. No file, setting,
 * Result or Session content goes into it (PD-OI-106).
 * @param {UiLocale} language the shown UI language
 * @returns {string}
 */
export const translationFeedbackUrl = (language) => {
  const url = new URL(TRANSLATION_FEEDBACK_FORM);
  url.searchParams.set('language', language);
  url.searchParams.set('version', GBDRAW_VERSION);
  url.searchParams.set('title', `[Translation] ${language}: `);
  return url.href;
};
