// @ts-check
// The settings pane: which display it shows (the stacked settings or the icon
// rail), the entry docked beside the rail, the section shown in a popup beside
// the Result, the one way to bring a setting into view, the one way to open a
// section by its key, and the display the browser remembers. Where a section
// is shown (the pane, the dock, the popup) is decided here, and so is how long
// the popup stays: it belongs to the shown Result. The transitions are in
// services/settings-pane-state.js; templates and other owners call the
// functions returned here and never assign the three refs.
import { readSettingsPaneDisplay, writeSettingsPaneDisplay } from '../services/browser-preferences.js';
import {
  closeDockedSettingsEntry, dockSettingsSection, floatSettingsSection, moveSettingsPopup, reconcileSettingsPaneState,
  returnSettingsSection, showSettingsRail, showSettingsStacked
} from '../services/settings-pane-state.js';
import { settingsEntryIn, settingsScrollOwner, settingsSectionElement, settingsView } from '../services/settings-sections.js';
import { placePopupNearPoint } from '../utils/popup-bounds.js';

/** @import { DiagramMode, SettingsPaneDisplay, SettingsPaneRefs } from '../services/settings-pane-state.js' */
/** @import { DrawnElementSectionKey, SettingsSectionKey } from '../services/settings-sections.js' */
/** @import { PopupAnchor, PopupBounds, PopupPlacement } from '../utils/popup-bounds.js' */

const { computed, nextTick, onMounted, onUnmounted, ref, watch } = window.Vue;

// The popup is as wide as the docked panel of the rail (the stylesheet narrows
// it in a narrow window); its height is its content's, read once it is shown.
const POPUP_WIDTH = 320;
const POPUP_EXPECTED_HEIGHT = 360;

/**
 * What the pane shows now. In the phone layout it is the stacked settings
 * without the switch, whatever display is remembered.
 * @typedef {object} SettingsPaneView
 * @property {SettingsPaneDisplay} display
 * @property {SettingsSectionKey | null} docked The docked rail entry.
 * @property {SettingsSectionKey | null} floating The section in the popup beside the Result; the pane does not show it.
 * @property {boolean} switchable Whether the display can be switched here.
 * @property {PhonePane} pane The pane the phone layout shows; wider layouts show both and ignore it.
 */

/**
 * The pane the phone layout (<= 768 px, decided by the stylesheet) shows: the
 * settings or the Preview. A view choice only: not in the Session, not a
 * History step, the same for every drawing, and not remembered by the browser.
 * @typedef {'settings' | 'preview'} PhonePane
 */

/** @typedef {(pane: PhonePane) => void} ShowPhonePane */

/**
 * Brings a setting into view and focuses it: opens the disclosures around
 * `focus`, docks the entry of the section when the rail shows, and scrolls
 * the settings alone. `align: 'start'` puts `target` at the top edge;
 * otherwise a target that is not fully in view is centred. A section in the
 * popup returns to the pane first. A target outside the settings takes the
 * focus where it is.
 * @typedef {(target: HTMLElement, options?: { align?: 'start', focus?: HTMLElement }) => Promise<void>} RevealSetting
 */

/**
 * Opens the section of `key` in the shown mode (for a key of the other mode,
 * the shown mode's section of the same title), brings its top to the top edge
 * of the settings, and focuses its heading, or `focus` when one is given. The
 * section index opens a section this way, and so does a click on the drawn
 * element of the docked section. A section in the popup returns to the pane
 * first. A key without a section in the shown mode opens nothing.
 * @typedef {(key: SettingsSectionKey, focus?: HTMLElement | null) => Promise<void>} RevealSettingsSection
 */

/**
 * Shows the section of `key` (as `RevealSettingsSection` reads a key) beside
 * the point of the Result that was clicked: in the popup, placed inside
 * `bounds` beside `anchor` and clear of the clicked item (`anchor.avoid`:
 * beside it, else below or above it) when there is room, with the focus on
 * its heading. The section docked beside the rail stays there and takes the
 * focus. Returns whether the shown mode has the section.
 * @typedef {(key: DrawnElementSectionKey, anchor: PopupAnchor, bounds: PopupBounds) => boolean} ShowSettingsBesideResult
 */

/**
 * @typedef {object} SettingsPaneState The refs of `state.js` the pane reads and, through its transitions, writes.
 * @property {{ value: DiagramMode }} mode
 * @property {{ readonly value: unknown }} svgContent The shown Result; none before the first Generate.
 * @property {{ value: string }} resultPanelTab The tab of the Result panel; `preview` shows the Result.
 * @property {{ readonly value: boolean }} layoutRepositionMode Layout edit: a press on the figure drags its items.
 */

/** @param {{ state: SettingsPaneRefs & SettingsPaneState }} options */
export const createSettingsPaneController = ({ state }) => {
  // The phone layout is read from the layout itself (the settings do not
  // scroll there), so its breakpoint stays in the stylesheet.
  /** @type {{ value: boolean }} */
  const phoneLayout = ref(false);
  /** @type {{ value: PhonePane }} */
  const phonePane = ref('settings');
  /** @type {ShowPhonePane} */
  const showPhonePane = (pane) => {
    phonePane.value = pane;
  };
  const settingsArea = () => document.querySelector('.settings-scroll');
  const popupElement = () => document.getElementById('settings-popup');
  // What had the focus when the popup opened: it takes the focus back.
  /** @type {Element | null} */
  let popupOpener = null;
  const readLayout = () => {
    const area = settingsArea();
    phoneLayout.value = Boolean(area && settingsScrollOwner(area) !== area);
  };

  // A page that opens in the remembered rail docks the first entry.
  if (readSettingsPaneDisplay() === 'rail') showSettingsRail(state, state.mode.value);

  /** @type {{ readonly value: SettingsPaneView }} */
  // Where a settings section can float beside the Result: a Result shown in Preview, Layout edit off.
  const settingsPopupAvailable = computed(() => Boolean(state.svgContent.value)
    && state.resultPanelTab.value === 'preview' && !state.layoutRepositionMode.value);
  const settingsPane = computed(() => {
    const rail = state.settingsPaneDisplay.value === 'rail' && !phoneLayout.value;
    return {
      display: rail ? 'rail' : 'stacked',
      docked: rail ? state.settingsPaneDocked.value : null,
      floating: state.settingsPopup.value?.section ?? null,
      switchable: !phoneLayout.value,
      pane: phonePane.value
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

  // Closes the popup. The focus it held goes back to what had it before, as
  // the focus of the Result's other popups does. Returns whether it was open.
  const closeSettingsPopup = () => {
    if (!state.settingsPopup.value) return false;
    const held = popupElement()?.contains(document.activeElement);
    returnSettingsSection(state);
    if (held) {
      void nextTick(() => {
        // Nothing else took the focus: it is nowhere (its section left the popup), or still on the popup's hidden close button.
        const active = document.activeElement;
        if ((active === document.body || popupElement()?.contains(active)) && popupOpener instanceof HTMLElement && popupOpener.isConnected) {
          popupOpener.focus({ preventScroll: true });
        }
      });
    }
    return true;
  };
  // A pane action that needs the section in the popup has it back first.
  /** @param {string | undefined} key */
  const returnToPane = async (key) => {
    if (!key || state.settingsPopup.value?.section !== key) return;
    returnSettingsSection(state);
    await nextTick();
  };

  /** @type {RevealSetting} */
  const revealSetting = async (target, { align, focus = target } = {}) => {
    const section = /** @type {HTMLElement | null} */ (target.closest('[data-settings-section]'));
    await returnToPane(section?.dataset.settingsSection);
    const area = section?.closest('.settings-scroll');
    if (!section || !area) return focus.focus();
    // A setting is shown in the settings pane, so the phone layout shows that pane before it is measured.
    showPhonePane('settings');
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
  /** @param {Element | null | undefined} section */
  const headingOf = (section) => /** @type {HTMLElement | null | undefined} */ (
    section?.querySelector(':scope > details:first-child > summary, :scope > [tabindex="-1"]'));
  /** @type {RevealSettingsSection} */
  const revealSettingsSection = async (key, focus) => {
    const area = settingsArea();
    const entry = settingsEntryIn(state.mode.value, key);
    await returnToPane(entry);
    const section = area && entry ? settingsSectionElement(area, entry) : null;
    const heading = headingOf(section);
    if (!section || !heading) return;
    if (heading.parentElement instanceof HTMLDetailsElement) heading.parentElement.open = true;
    await revealSetting(section, { align: 'start', focus: focus ?? heading });
  };

  /** @type {ShowSettingsBesideResult} */
  const showSettingsBesideResult = (key, anchor, bounds) => {
    const entry = settingsEntryIn(state.mode.value, key);
    if (!entry || !settingsPopupAvailable.value) return false;
    if (settingsPane.value.docked === entry) {
      closeSettingsPopup();
      void revealSettingsSection(entry);
      return true;
    }
    // Clear of the clicked item when the popup fits beside, below, or above it; a placement
    // that lies over the item (a ring fills the Result) gives way to one beside the point.
    /** @param {number} width @param {number} height */
    const place = (width, height) => {
      /** @param {PopupAnchor} point */
      const near = (point) => placePopupNearPoint({ point, width, height, bounds });
      const placed = near(anchor);
      const over = anchor.avoid && placed.x < anchor.avoid.right && placed.x + width > anchor.avoid.left
        && placed.y < anchor.avoid.bottom && placed.y + Math.min(height, placed.maxHeight) > anchor.avoid.top;
      return over ? near({ x: anchor.x, y: anchor.y }) : placed;
    };
    const active = document.activeElement;
    if (!state.settingsPopup.value) popupOpener = active === document.body ? null : active;
    floatSettingsSection(state, entry, place(POPUP_WIDTH, POPUP_EXPECTED_HEIGHT));
    void nextTick(() => {
      const popup = popupElement();
      const section = popup && settingsSectionElement(popup, entry);
      if (!popup || !section || state.settingsPopup.value?.section !== entry) return;
      // The content's own height: what the section holds, not what the popup shows of it.
      moveSettingsPopup(state, place(popup.offsetWidth, popup.offsetHeight + section.scrollHeight - section.clientHeight));
      headingOf(section)?.focus({ preventScroll: true });
    });
    return true;
  };

  // R10: a mode written by anything (the mode buttons, Undo, a loaded Session) reconciles the docked entry and closes the popup.
  watch(() => state.mode.value, (/** @type {DiagramMode} */ mode) => reconcileSettingsPaneState(state, mode), { flush: 'sync' });
  // The popup belongs to the shown Result while Layout edit is off: whatever ends that (another tab
  // of the Result panel, an Undo to before the first Generate, a load, Layout edit turned on) closes
  // it, and nothing opens it meanwhile. The same predicate decides both.
  watch(() => settingsPopupAvailable.value, (/** @type {boolean} */ available) => {
    if (!available) closeSettingsPopup();
  }, { flush: 'sync' });
  onMounted(() => {
    readLayout();
    window.addEventListener('resize', readLayout);
  });
  onUnmounted(() => window.removeEventListener('resize', readLayout));

  return {
    settingsPane, showPhonePane, showSettingsPaneRail, showSettingsPaneStacked, closeDockedSettings, revealSetting, revealSettingsSection,
    showSettingsBesideResult, closeSettingsPopup,
    // A drag of the popup by its heading.
    moveSettingsPopup: (/** @type {PopupPlacement} */ placement) => moveSettingsPopup(state, placement)
  };
};
