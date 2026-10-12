// @ts-check
// The settings pane: which display it shows (the stacked settings or the icon
// rail), the entry docked beside the rail, the one way to bring a setting
// into view, the one way to open a section by its key, and the display the
// browser remembers. The transitions are in
// services/settings-pane-state.js; templates and other owners call the
// functions returned here and never assign the two refs.
import { readSettingsPaneDisplay, writeSettingsPaneDisplay } from '../services/browser-preferences.js';
import {
  closeDockedSettingsEntry, dockSettingsSection, reconcileSettingsPaneState, showSettingsRail, showSettingsStacked
} from '../services/settings-pane-state.js';
import { settingsEntryIn, settingsScrollOwner, settingsSectionElement, settingsView } from '../services/settings-sections.js';

/** @import { DiagramMode, SettingsPaneDisplay, SettingsPaneRefs } from '../services/settings-pane-state.js' */
/** @import { SettingsSectionKey } from '../services/settings-sections.js' */

const { computed, nextTick, onMounted, onUnmounted, reactive, ref, watch } = window.Vue;

/**
 * What the pane shows now. In the phone layout it is the stacked settings
 * without the switch, whatever display is remembered.
 * @typedef {object} SettingsPaneView
 * @property {SettingsPaneDisplay} display
 * @property {SettingsSectionKey | null} docked The docked rail entry.
 * @property {boolean} switchable Whether the display can be switched here.
 */

/**
 * The section the pane last opened by its key, and where that left its top
 * edge. The section index reads it to hold its mark. It is an object of its
 * own, not a part of the view the root template binds, so opening a section
 * re-renders the index alone.
 * @typedef {{ last: { key: SettingsSectionKey, top: number } | null }} SettingsReveal
 */

/**
 * Brings a setting into view and focuses it: opens the disclosures around
 * `focus`, docks the entry of the section when the rail shows, and scrolls
 * the settings alone. `align: 'start'` puts `target` at the top edge;
 * otherwise a target that is not fully in view is centred. A target outside
 * the settings takes the focus where it is.
 * @typedef {(target: HTMLElement, options?: { align?: 'start', focus?: HTMLElement }) => Promise<void>} RevealSetting
 */

/**
 * Opens the section of `key` in the shown mode (for a key of the other mode,
 * the shown mode's section of the same title), brings its top to the top edge
 * of the settings, and focuses its heading, or `focus` when one is given. The
 * section index and a click on a drawn element of the Result both open a
 * section this way. A key without a section in the shown mode opens nothing.
 * @typedef {(key: SettingsSectionKey, focus?: HTMLElement | null) => Promise<void>} RevealSettingsSection
 */

/** @param {{ state: SettingsPaneRefs & { mode: { value: DiagramMode } } }} options */
export const createSettingsPaneController = ({ state }) => {
  // The phone layout is read from the layout itself (the settings do not
  // scroll there), so its breakpoint stays in the stylesheet.
  /** @type {{ value: boolean }} */
  const phoneLayout = ref(false);
  /** @type {SettingsReveal} */
  const settingsReveal = reactive({ last: null });
  const settingsArea = () => document.querySelector('.settings-scroll');
  const readLayout = () => {
    const area = settingsArea();
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
    if (!section || !area) return focus.focus();
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

  // The heading of a section is the summary of a card that is one disclosure,
  // or the heading element that takes the focus (Input Genomes, Comparison).
  /** @type {RevealSettingsSection} */
  const revealSettingsSection = async (key, focus) => {
    const area = settingsArea();
    const entry = settingsEntryIn(state.mode.value, key);
    const section = area && entry ? settingsSectionElement(area, entry) : null;
    const heading = /** @type {HTMLElement | null | undefined} */ (
      section?.querySelector(':scope > details:first-child > summary, :scope > [tabindex="-1"]'));
    if (!entry || !section || !heading) return;
    if (heading.parentElement instanceof HTMLDetailsElement) heading.parentElement.open = true;
    await revealSetting(section, { align: 'start', focus: focus ?? heading });
    settingsReveal.last = { key: entry, top: section.getBoundingClientRect().top };
  };

  // R10: a mode written by anything (the mode buttons, Undo, a loaded Session) reconciles the docked entry.
  watch(() => state.mode.value, (/** @type {DiagramMode} */ mode) => reconcileSettingsPaneState(state, mode), { flush: 'sync' });
  onMounted(() => {
    readLayout();
    window.addEventListener('resize', readLayout);
  });
  onUnmounted(() => window.removeEventListener('resize', readLayout));

  return { settingsPane, showSettingsPaneRail, showSettingsPaneStacked, closeDockedSettings, revealSetting, revealSettingsSection, settingsReveal };
};
