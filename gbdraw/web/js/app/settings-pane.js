// @ts-check
// The settings pane: which display it shows (the stacked settings or the icon
// rail), the entry docked beside the rail, the one way to bring a setting
// into view, and the display the browser remembers. The transitions are in
// services/settings-pane-state.js; templates and other owners call the
// functions returned here and never assign the two refs.
import { readSettingsPaneDisplay, writeSettingsPaneDisplay } from '../services/browser-preferences.js';
import {
  closeDockedSettingsEntry, dockSettingsSection, reconcileSettingsPaneState, showSettingsRail, showSettingsStacked
} from '../services/settings-pane-state.js';
import { settingsScrollOwner, settingsView } from '../services/settings-sections.js';

/** @import { DiagramMode, SettingsPaneDisplay, SettingsPaneRefs } from '../services/settings-pane-state.js' */
/** @import { SettingsSectionKey } from '../services/settings-sections.js' */

const { computed, nextTick, onMounted, onUnmounted, ref, watch } = window.Vue;

/**
 * What the pane shows now. In the phone layout it is the stacked settings
 * without the switch, whatever display is remembered.
 * @typedef {object} SettingsPaneView
 * @property {SettingsPaneDisplay} display
 * @property {SettingsSectionKey | null} docked The docked rail entry.
 * @property {boolean} switchable Whether the display can be switched here.
 */

/**
 * Brings a setting into view and focuses it: opens the disclosures around
 * `focus`, docks the entry of the section when the rail shows, and scrolls
 * the settings alone. `align: 'start'` puts `target` at the top edge;
 * otherwise a target that is not fully in view is centred.
 * @typedef {(target: HTMLElement, options?: { align?: 'start', focus?: HTMLElement }) => Promise<void>} RevealSetting
 */

/** @param {{ state: SettingsPaneRefs & { mode: { value: DiagramMode } } }} options */
export const createSettingsPaneController = ({ state }) => {
  // The phone layout is read from the layout itself (the settings do not
  // scroll there), so its breakpoint stays in the stylesheet.
  /** @type {{ value: boolean }} */
  const phoneLayout = ref(false);
  const readLayout = () => {
    const area = document.querySelector('.settings-scroll');
    phoneLayout.value = Boolean(area && settingsScrollOwner(area) !== area);
  };

  // A page that opens in the remembered rail docks the first entry.
  if (readSettingsPaneDisplay() === 'rail') showSettingsRail(state, state.mode.value);

  /** @type {{ readonly value: SettingsPaneView }} */
  const settingsPane = computed(() => {
    const rail = state.settingsPaneDisplay.value === 'rail' && !phoneLayout.value;
    return {
      display: rail ? 'rail' : 'stacked',
      docked: rail ? state.settingsPaneDocked.value : null,
      switchable: !phoneLayout.value
    };
  });

  /** @param {SettingsSectionKey} [section] The section whose entry is docked. */
  const showSettingsPaneRail = (section) => {
    showSettingsRail(state, state.mode.value, section);
    writeSettingsPaneDisplay('rail');
  };
  const showSettingsPaneStacked = () => {
    showSettingsStacked(state);
    writeSettingsPaneDisplay('stacked');
  };
  const closeDockedSettings = () => closeDockedSettingsEntry(state);

  /** @type {RevealSetting} */
  const revealSetting = async (target, { align, focus = target } = {}) => {
    const section = /** @type {HTMLElement | null} */ (target.closest('[data-settings-section]'));
    const area = section?.closest('.settings-scroll');
    if (!section || !area) return;
    for (let details = focus.closest('details'); details && area.contains(details);
      details = details.parentElement?.closest('details') ?? null) {
      details.open = true;
    }
    if (settingsPane.value.display === 'rail') {
      dockSettingsSection(state, state.mode.value, /** @type {SettingsSectionKey} */ (section.dataset.settingsSection));
    }
    await nextTick();
    const { owner, top, bottom } = settingsView(area);
    const box = target.getBoundingClientRect();
    const inView = box.top >= top && box.bottom <= bottom;
    const distance = align === 'start' ? box.top - top : box.top + box.height / 2 - (top + bottom) / 2;
    if (align === 'start' || !inView) owner?.scrollBy({ top: distance, behavior: 'instant' });
    focus.focus({ preventScroll: true });
  };

  // R10: a mode written by anything (the mode buttons, Undo, a loaded Session) reconciles the docked entry.
  watch(() => state.mode.value, (/** @type {DiagramMode} */ mode) => reconcileSettingsPaneState(state, mode), { flush: 'sync' });
  onMounted(() => {
    readLayout();
    window.addEventListener('resize', readLayout);
  });
  onUnmounted(() => window.removeEventListener('resize', readLayout));

  return { settingsPane, showSettingsPaneRail, showSettingsPaneStacked, closeDockedSettings, revealSetting };
};
