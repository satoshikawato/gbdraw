// @ts-check
import { state } from './state.js';
import {
  colorValueForMode,
  colorValueMode,
  toNativeColorInputValue
} from './utils/color-utils.js';

import { normalizeUserFacingError, operationErrorTitle, generationRecoveryGuidance } from './utils/error-normalization.js';
import { italicSegments } from './utils/markup-preview.js';
import { parseMarkupInertly } from './utils/inert-markup.js';
import { getLinearSourceDefaultDefinition, getLinearSourceDefaultSubtitle, resolveLinearRecordEffectiveDefinition } from './services/linear-sources.js';
import { SETTINGS_SECTIONS, settingsEntriesOf, settingsEntryOf, settingsSectionElement, settingsView } from './services/settings-sections.js';
import { SETTING_NOTE_APPLIES, settingApplicationNote } from './services/setting-application.js';
/** @import { UserFacingError } from './utils/error-normalization.js' */
/** @import { SettingApplicationId } from './services/setting-application.js' */

/** @import { SettingsEntry, SettingsSectionKey } from './services/settings-sections.js' */
/** @import { SettingsPaneView } from './app/settings-pane.js' */

const { ref, reactive, computed, inject, nextTick, watch, useId, onMounted, onBeforeUnmount } = window.Vue;

/** @typedef {{ element: HTMLElement, returnFocus: HTMLElement | SVGElement | null }} FocusScope */
/** @type {FocusScope[]} */
const focusScopes = [];

// UI-02, UI-03: a dialog takes focus when it mounts (its element marked
// data-dialog-initial-focus, else its first enabled button, else the dialog)
// and, when it unmounts with focus inside it or lost, returns focus to the
// element that had it before. A dialog that closes while it holds the return
// target of a later one hands that one its own return target.
export const dialogFocus = {
  /** @param {HTMLElement} element */
  beforeMount(element) {
    const active = document.activeElement;
    // An SVG opener (a preview Legend row or record definition) takes focus back too.
    const opener = active instanceof HTMLElement || active instanceof SVGElement ? active : null;
    focusScopes.push({ element, returnFocus: opener !== document.body ? opener : null });
  },
  /** @param {HTMLElement} element */
  mounted(element) {
    const dialog = element.matches('[role="dialog"]') ? element : element.querySelector('[role="dialog"]');
    const target = element.querySelector('[data-dialog-initial-focus]')
      || dialog?.querySelector('button:not(:disabled)') || dialog;
    if (target instanceof HTMLElement) target.focus({ preventScroll: true });
  },
  /** @param {HTMLElement} element */
  beforeUnmount(element) {
    const index = focusScopes.findIndex((scope) => scope.element === element);
    if (index < 0) return;
    const [{ returnFocus }] = focusScopes.splice(index, 1);
    for (const scope of focusScopes) {
      if (scope.returnFocus && element.contains(scope.returnFocus)) scope.returnFocus = returnFocus;
    }
    const active = document.activeElement;
    if (active && active !== document.body && !element.contains(active)) return;
    if (returnFocus?.isConnected) returnFocus.focus({ preventScroll: true });
  }
};

// UI-02: a modal choice. The panel is the dialog, named by the heading and
// described by the text that `labelledby` and `describedby` name. Escape and a
// backdrop click emit `cancel`, which a use binds to its Cancel handler; a use
// binds Tab to `trapDialogFocus`. While `busy` (a choice still applying), its
// controls are disabled and it says so (D-12), or says `busyStatus` when a use
// is busy with something else (the label fit dialog finishing a row, LP-D60). Focus stays in it: the panel
// takes focus before the controls are disabled, so Escape still reaches the
// dialog, and the first enabled control gets it back if the dialog stays open.
export const ChoiceDialog = {
  template: '#choice-dialog-template',
  props: {
    labelledby: { type: String, required: true },
    describedby: { type: String, default: null },
    busy: { type: Boolean, default: false },
    busyStatus: { type: String, default: null }
  },
  emits: ['cancel'],
  /** @param {{ busy: boolean }} props */
  setup(props) {
    /** @type {{ value: HTMLElement | null }} */
    const panel = ref(null);
    const focusLost = () => !document.activeElement || document.activeElement === document.body;
    watch(() => props.busy, async (/** @type {boolean} */ busy) => {
      const element = panel.value;
      if (!element) return;
      if (busy) {
        if (focusLost() || element.contains(document.activeElement)) element.focus({ preventScroll: true });
        return;
      }
      await nextTick();
      if (element.isConnected && (focusLost() || document.activeElement === element)) {
        /** @type {HTMLElement | null} */ (element.querySelector('button:not(:disabled)'))?.focus({ preventScroll: true });
      }
    });
    return { panel };
  }
};

/** @param {Node} node @returns {string} */
const visibleText = (node) => {
  if (node.nodeType === Node.TEXT_NODE) return node.textContent || '';
  if (!(node instanceof Element) || node.matches('button, input, select, textarea, [hidden], [aria-hidden="true"]')) return '';
  return Array.from(node.childNodes, visibleText).join(' ');
};

// R12 (UI-04): a tip is named after what it explains, "Help: <label>": the
// author name or visible text of the element that holds it, after the name of
// the row group it is in. `label` names a tip whose holder has no such text;
// the layout warning beside Generate is named "Layout warning" (not help).
/** @param {Element | null} tip */
const helpTipName = (tip) => {
  const holder = tip?.parentElement;
  const group = holder?.closest('[role="group"][aria-label]')?.getAttribute('aria-label') || '';
  const label = holder?.getAttribute('aria-label') || (holder ? visibleText(holder) : '');
  const name = `${group} ${label}`.replace(/\s+/g, ' ').trim();
  return name ? `Help: ${name}` : 'Help';
};

// Loading indicator (Owner decision 2026-10-11): the logo drawn in one stroke. The markup and the
// motion live in index.html (#gbdraw-loader-template and the .gbdraw-loader
// styles); each instance gets its own mask ids. The boot splash holds a static
// copy of the same markup because it renders before Vue mounts.
export const GbdrawLoader = {
  template: '#gbdraw-loader-template',
  setup() {
    const id = `gbdraw-loader-${useId()}`;
    return { ids: { a: `${id}-a`, b: `${id}-b` } };
  }
};

// Every tip is a disclosure button (PD-OI-057). One text source feeds the
// visible tooltip and the accessible description, identified by `id` or an
// automatic id; the control that a tip explains references that id through
// aria-describedby. The description is a hidden element teleported to <body>,
// so the text around the tip stays its visible label and screen readers do
// not read the tip text twice. Mouse hover and keyboard focus show the
// tooltip, a click or tap toggles it, and Escape or blur closes it. Place
// tips outside <label> so they never change the accessible name of a label or
// its control. The default slot replaces the question mark (the layout warning
// beside the compact Generate button shows a warning mark).
export const HelpTip = {
  template: '#help-tip-template',
  props: ['text', 'id', 'label'],
  /** @param {{ text?: string, id?: string | null, label?: string }} props */
  setup(props) {
    const automaticId = `help-tip-${useId()}`;
    const descriptionId = computed(() => props.id || automaticId);
    const name = ref('Help');
    const hovered = ref(false);
    const keyboardFocus = ref(false);
    const pinned = ref(false);
    const visible = computed(() => hovered.value || keyboardFocus.value || pinned.value);
    const style = reactive({ top: '0px', left: '0px' });
    const trigger = ref(null);
    const position = () => {
      if (!trigger.value) return;
      const rect = trigger.value.getBoundingClientRect();

      // Get viewport width and tooltip max width
      const viewportWidth = window.innerWidth;
      const tooltipMaxWidth = 260; // max-w-250px + extra margin
      const halfWidth = tooltipMaxWidth / 2;
      const gap = 12; // gap between icon and tooltip

      // Centered on the icon
      let left = rect.left + rect.width / 2;
      let top = rect.top - gap;
      let transform = 'translate(-50%, -100%)'; // above the icon

      // Prevent overflow on left/right (keep at least 10px from viewport edges)
      if (left < halfWidth + 10) {
        left = halfWidth + 10;
      } else if (left > viewportWidth - halfWidth - 10) {
        left = viewportWidth - halfWidth - 10;
      }

      // Prevent overflow at the top (if too close to top edge, show below)
      // Considering header and browser frame, if y < 60px, show below
      if (rect.top < 60) {
        top = rect.bottom + gap;
        transform = 'translate(-50%, 0)'; // transform for below display
      }

      style.top = `${top}px`;
      style.left = `${left}px`;
      style.transform = transform;
    };
    watch(visible, (/** @type {boolean} */ open) => { if (open) position(); });
    // A touch tap also reports pointerenter; only a hovering pointer shows it.
    const onEnter = (/** @type {PointerEvent} */ event) => { hovered.value = event.pointerType !== 'touch'; };
    const onLeave = () => { hovered.value = false; };
    const onFocus = () => { keyboardFocus.value = Boolean(trigger.value?.matches?.(':focus-visible')); };
    const close = () => {
      hovered.value = false;
      keyboardFocus.value = false;
      pinned.value = false;
    };
    const toggle = () => { pinned.value = !pinned.value; };
    onMounted(() => { name.value = helpTipName(trigger.value?.closest('.help-tip') || null); });
    return { descriptionId, name, visible, style, trigger, onEnter, onLeave, onFocus, close, toggle };
  }
};

// OV-912: the preview zoom and pan change on every wheel tick and pan. Only
// the zoom label and the preview surface read them, so a change renders the
// label and writes the surface transform instead of rendering the whole app.
// The values stay in `state` for History, the Session, and the pan and zoom
// owner (app/ui.js), which also writes the transform during a gesture.
export const PreviewZoomPercent = {
  setup: () => ({ zoom: state.zoom }),
  template: '{{ Math.round(zoom * 100) }}%'
};

/** @type {WeakMap<HTMLElement, () => void>} */
const previewTransformStops = new WeakMap();
// `v-preview-transform` on the preview surface: its transform follows the pan
// and zoom from its creation (before it is inserted) until it unmounts.
export const previewTransform = {
  /** @param {HTMLElement} element */
  beforeMount: (element) => {
    previewTransformStops.set(element, watch(
      () => /** @type {const} */ ([state.canvasPan.x, state.canvasPan.y, state.zoom.value, state.isPanning.value]),
      (/** @type {readonly [number, number, number, boolean]} */ [x, y, zoom, panning]) => {
        element.style.transform = `translate(${x}px, ${y}px) scale(${zoom})`;
        element.style.transformOrigin = 'top center';
        element.style.transition = panning ? 'none' : 'transform 0.2s';
      },
      { immediate: true }
    ));
  },
  /** @param {HTMLElement} element */
  unmounted: (element) => {
    previewTransformStops.get(element)?.();
    previewTransformStops.delete(element);
  }
};

// `v-preview-markup="svgContent"` on the preview surface: its children are the
// Result's markup, as `v-html` would set them, at the same points of the
// patch. The markup is parsed outside the page document (OV-1030).
/** @param {HTMLElement} element @param {string | null | undefined} markup */
const setPreviewMarkup = (element, markup) => {
  element.replaceChildren(parseMarkupInertly(markup));
};
export const previewMarkup = {
  /** @param {HTMLElement} element @param {{ value: string | null | undefined }} binding */
  beforeMount: (element, { value }) => setPreviewMarkup(element, value),
  /** @param {HTMLElement} element @param {{ value: string | null | undefined, oldValue: string | null | undefined }} binding */
  beforeUpdate: (element, { value, oldValue }) => {
    if (value !== oldValue) setPreviewMarkup(element, value);
  }
};

// The heading of a settings section: its icon and title from the catalog
// (services/settings-sections.js), which the section index shows too.
export const SettingsSectionTitle = {
  props: { section: { type: String, required: true } },
  /** @param {{ section: SettingsSectionKey }} props */
  setup: (props) => ({ entry: SETTINGS_SECTIONS.find((section) => section.key === props.section) }),
  template: '<i class="ph" :class="entry.icon" aria-hidden="true"></i> {{ entry.title }}'
};

// The timing note of a `SETTING_APPLICATION` row (services/setting-application.js),
// placed inside that row's `data-setting` scope: "Applies on Generate" with
// `data-applies-on-generate="<key>"`, "Live edit" for a live block inside a
// Generate scope, or the row's own qualified note as a hint line. A live
// setting row the shown Result cannot draw says "Applies on Generate"
// (the injected `SETTING_NOTE_APPLIES` predicate of the committed request).
export const SettingNote = {
  props: { setting: { type: String, required: true } },
  /** @param {{ setting: SettingApplicationId }} props */
  setup: (props) => {
    const applies = /** @type {(requestPath: string) => boolean} */ (inject(SETTING_NOTE_APPLIES, () => true));
    return { note: computed(() => settingApplicationNote(props.setting, applies)) };
  },
  template: `<span v-if="note.text" :class="note.qualified ? 'ui-hint block' : 'apply-note'"
    :data-applies-on-generate="note.application === 'generate' ? setting : null"><i v-if="note.application === 'generate'"
    class="ph ph-clock" aria-hidden="true"></i> {{ note.text }}</span>`
};

// `v-settings-section="'<key>'"` on the element of a settings section: it
// names the section (`data-settings-section`) and tells the stylesheet whether
// the docked rail entry shows it (`data-settings-in-dock`). The pane says which
// display shows, so the stacked display ignores the second attribute.
/** @param {HTMLElement} element @param {{ value: SettingsSectionKey }} binding */
export const settingsSection = (element, { value }) => {
  element.dataset.settingsSection = value;
  element.dataset.settingsInDock = String(settingsEntryOf(value) === state.settingsPaneDocked.value);
};

// The section index of the settings: one button per section of the shown
// mode, in card order, in the form of the pane's display (`view`, from
// app/settings-pane.js). `controls` is the id of the element that holds the
// sections (`data-settings-section`).
// - Stacked: a row of icons. An entry opens its section, brings it to the top
//   of the visible settings at once, and moves the focus to its heading; the
//   marked entry (`aria-current`) is the last section whose top has reached
//   the top edge of the visible settings, and a section the pane opened (an
//   activated entry, or a click on its drawn element in the Result) stays
//   marked while it stays where that left it, so a short last section keeps
//   its mark. The switch to the rail ends the row.
// - Rail: a column of icons over short names, the switch back first. An entry
//   docks its section beside the rail, opened and focused in the same way; the
//   marked entry is the docked one, and the close button leaves the rail alone.
//   A column too low for its entries scrolls: the docked entry and the focused
//   one are brought into view.
// An entry never closes a section and writes no setting. Both forms open a
// section through `reveal`, the pane's one way to open a section by its key.
// The text shown for a hovered or keyboard-focused button (the section title,
// or the words of a switch or of the close button) is the index's own: each
// button is named by its aria-label ("Go to <title>" for an entry; the
// section's own disclosure keeps the bare title), so the text is hidden from
// assistive technology and no `title` attribute repeats it. When the rail
// column scrolls, the text of a hovered entry goes and the text of the focused
// entry follows it.
/**
 * The pane's way to open a section: open, its top at the top edge of the
 * settings, the focus on its heading, or on `focus` when one is given.
 * @typedef {(key: SettingsSectionKey, focus?: HTMLElement | null) => unknown} RevealSection
 */
export const SettingsSectionIndex = {
  template: '#settings-section-index-template',
  props: {
    mode: { type: String, required: true },
    controls: { type: String, required: true },
    view: { type: Object, required: true },
    reveal: { type: Function, required: true }
  },
  emits: ['rail', 'stacked', 'close'],
  /**
   * @param {{ mode: 'circular' | 'linear', controls: string, view: SettingsPaneView, reveal: RevealSection }} props
   * @param {{ emit: (event: 'rail' | 'stacked' | 'close', section?: SettingsSectionKey) => void }} context
   */
  setup(props, { emit }) {
    /** @type {{ value: boolean }} */
    const rail = computed(() => props.view.display === 'rail');
    /** @type {{ value: readonly SettingsEntry[] }} */
    const entries = computed(() => settingsEntriesOf(props.mode, props.view.display));
    // The section at the top of the stacked settings.
    /** @type {{ value: SettingsSectionKey }} */
    const current = ref(SETTINGS_SECTIONS[0].key);
    /** @type {{ value: SettingsSectionKey | null }} */
    const marked = computed(() => (rail.value ? props.view.docked : current.value));
    /** @type {{ value: { text: string, side: 'below' | 'right', left: number, top: number } | null }} */
    const tip = ref(null);
    /** @type {{ value: HTMLElement | null }} */
    const tipElement = ref(null);
    // The button whose text is shown.
    /** @type {HTMLElement | null} */
    let tipButton = null;
    /** @type {{ value: HTMLElement | null }} */
    const root = ref(null);
    /** @type {{ key: SettingsSectionKey, top: number } | null} */
    let held = null;
    // Where the stacked settings were scrolled when the rail took their place.
    let stackedTop = 0;
    // The space right of the cards inside the scroll area (its padding and a
    // classic scrollbar, when it shows one): the row ends where the cards end.
    /** @type {{ value: number }} */
    const inset = ref(0);
    const cardEdge = new ResizeObserver(() => {
      inset.value = area().offsetWidth - area().clientWidth + parseFloat(getComputedStyle(area()).paddingRight);
    });

    const area = () => /** @type {HTMLElement} */ (document.getElementById(props.controls));
    /** @param {SettingsSectionKey} key */
    const sectionOf = (key) => settingsSectionElement(area(), key);
    /** @param {string} selector A button of the index. */
    const focusButton = (selector) => /** @type {HTMLElement | null | undefined} */ (root.value?.querySelector(selector))
      ?.focus({ preventScroll: true });
    const mark = () => {
      if (rail.value || (held && sectionOf(held.key)?.getBoundingClientRect().top === held.top)) return;
      held = null;
      const edge = settingsView(area()).top;
      const reached = entries.value.filter((section) => (sectionOf(section.key)?.getBoundingClientRect().top ?? Infinity) <= edge + 1);
      current.value = (reached.at(-1) ?? entries.value[0]).key;
    };
    // Opens a section and brings it to the top. The focus goes to its heading
    // unless `focus` names a button of the index.
    /** @param {SettingsSectionKey} key @param {string} [focus] */
    const bringToTop = (key, focus) =>
      props.reveal(key, focus ? /** @type {HTMLElement | null | undefined} */ (root.value?.querySelector(focus)) : null);
    // A section the pane opened holds the mark. A section without an entry in
    // this form (About & Citation in the row) leaves the mark to the geometry.
    watch(() => props.view.revealed, (/** @type {SettingsPaneView['revealed']} */ revealed) => {
      held = revealed && entries.value.some((entry) => entry.key === revealed.key) ? revealed : null;
      if (held) current.value = held.key; else mark();
    }, { flush: 'sync' });
    /** @param {SettingsSectionKey} key */
    const activate = (key) => {
      hideTip();
      return bringToTop(key);
    };
    // Switching keeps the place: the rail docks the marked section, and the
    // stacked settings show the docked section at the top or, after the rail
    // alone, return to where they were. The focus goes to the switch of the
    // new display.
    const showRail = async () => {
      hideTip();
      stackedTop = area().scrollTop;
      emit('rail', current.value);
      await nextTick();
      if (props.view.docked) await bringToTop(props.view.docked, '.settings-index-switch');
    };
    const showStacked = async () => {
      hideTip();
      const docked = props.view.docked;
      emit('stacked');
      await nextTick();
      if (docked) return bringToTop(docked, '.settings-index-switch');
      held = null;
      area().scrollTop = stackedTop;
      mark();
      return focusButton('.settings-index-switch');
    };
    const closeDocked = async () => {
      hideTip();
      const docked = props.view.docked;
      emit('close');
      await nextTick();
      // Without `preventScroll`: a rail that scrolls brings the entry into view.
      /** @type {HTMLElement | null | undefined} */ (root.value?.querySelector(`[data-settings-entry="${docked}"]`))?.focus();
    };

    // Beside a button of the rail; else below the button and centred on it,
    // then kept 8 px inside the viewport from the rendered width.
    /** @param {string} text @param {HTMLElement} button */
    const showTip = async (text, button) => {
      tipButton = button;
      const box = button.getBoundingClientRect();
      tip.value = rail.value && root.value?.contains(button)
        ? { text, side: 'right', left: box.right + 6, top: box.top + box.height / 2 }
        : { text, side: 'below', left: box.left + box.width / 2, top: box.bottom + 6 };
      await nextTick();
      const half = (tipElement.value?.offsetWidth ?? 0) / 2;
      if (tip.value?.side === 'below') tip.value.left = Math.max(8 + half, Math.min(tip.value.left, window.innerWidth - 8 - half));
    };
    const hideTip = () => {
      tipButton = null;
      tip.value = null;
    };
    // Escape hides a shown text and stops there; with none shown it goes on to the page's own handler.
    /** @param {KeyboardEvent} event */
    const onEscape = (event) => {
      if (!tip.value) return;
      event.stopPropagation();
      hideTip();
    };
    // The listeners of a button that shows `text`. A touch tap also reports
    // pointerenter, so only a hovering pointer shows it, and focus shows it
    // only for the keyboard.
    /** @param {string} text */
    const tipOn = (text) => ({
      pointerenter: (/** @type {PointerEvent} */ event) => {
        if (event.pointerType !== 'touch') showTip(text, /** @type {HTMLElement} */ (event.currentTarget));
      },
      pointerleave: hideTip,
      focus: (/** @type {FocusEvent} */ event) => {
        const button = /** @type {HTMLElement} */ (event.currentTarget);
        if (button.matches(':focus-visible')) showTip(text, button);
      },
      blur: hideTip
    });

    // Scroll events do not bubble: one capturing listener hears the scroll area, the page, and the rail column.
    /** @param {Event} event */
    const onScroll = (event) => {
      if (event.target === area() || event.target === document) return mark();
      if (!tip.value || !tipButton || event.target !== tipButton.parentElement) return;
      // The column moved under the pointer; the focused entry it brought into view keeps its text.
      if (tipButton === document.activeElement) showTip(tip.value.text, tipButton); else hideTip();
    };
    onMounted(() => {
      cardEdge.observe(area());
      window.addEventListener('scroll', onScroll, { capture: true, passive: true });
    });
    onBeforeUnmount(() => {
      cardEdge.disconnect();
      window.removeEventListener('scroll', onScroll, { capture: true });
    });
    // The docked entry stays in view in a column that scrolls.
    watch(marked, () => nextTick(() => {
      if (rail.value) root.value?.querySelector('[aria-current="true"]')?.scrollIntoView({ block: 'nearest' });
    }));
    // A mode change may remove the entry whose text is shown, and the marked
    // section; so may a display that changes with the layout (the phone width).
    // Two sources: one getter that returns both would run on every change of the view.
    watch([() => props.mode, () => props.view.display], () => {
      hideTip();
      held = null;
      nextTick(mark);
    });

    return { rail, entries, marked, tip, tipElement, root, inset, activate, showRail, showStacked, closeDocked, onEscape, tipOn };
  }
};

export const AutoValueField = {
  props: ['visible', 'text'],
  template: `
    <div class="auto-value-field">
      <span v-if="visible" class="auto-value-placeholder">{{ text }}</span>
      <slot></slot>
    </div>
  `
};

/**
 * One item of an action menu. `live` announces a label that changes after the
 * action (Copied!, Copy failed).
 * @typedef {{ label: string, title?: string, icon?: string, disabled?: boolean, live?: boolean, describedby?: string, run: () => unknown }} ActionMenuItem
 * @typedef {{ section: string, help?: string, helpLabel?: string, helpId?: string }} ActionMenuSection
 *   A heading that groups the items after it, with an optional help tip.
 */

// A disclosure that holds a group of secondary actions, such as the copy and
// FASTA actions of the similarity popups (UX-B01..B03). The items open in
// place, inside the popup's scroll area, so no item is clipped.
export const ActionMenu = {
  components: { HelpTip },
  props: {
    label: { type: String, required: true },
    icon: { type: String, default: '' },
    title: { type: String, default: '' },
    /** @type {{ type: ArrayConstructor, default: () => (ActionMenuItem | ActionMenuSection)[] }} */
    items: { type: Array, default: () => [] }
  },
  template: `
    <details class="action-menu" data-nodrag="true">
      <summary class="action-menu-toggle" :title="title || label">
        <i v-if="icon" :class="['ph', icon]" aria-hidden="true"></i><span>{{ label }}</span>
        <i class="ph ph-caret-down action-menu-caret" aria-hidden="true"></i>
      </summary>
      <div class="action-menu-items">
        <template v-for="(item, index) in items" :key="index + '-' + (item.section || item.label)">
          <div v-if="item.section" class="action-menu-section">
            <span>{{ item.section }}</span>
            <help-tip v-if="item.help" :id="item.helpId" :label="item.helpLabel" :text="item.help"></help-tip>
          </div>
          <button v-else type="button" class="action-menu-item"
                  :disabled="item.disabled" :title="item.title || item.label" :aria-describedby="item.describedby || null"
                  :aria-live="item.live ? 'polite' : null" :aria-atomic="item.live ? 'true' : null"
                  @click="item.run()">
            <i v-if="item.icon" :class="['ph', item.icon]" aria-hidden="true"></i><span>{{ item.label }}</span>
          </button>
        </template>
      </div>
    </details>
  `
};

export const ColorValueControl = {
  props: {
    modelValue: { default: null },
    fallback: { type: String, default: '#000000' },
    allowNone: { type: Boolean, default: true },
    ariaLabel: { type: String, default: 'Color value' },
    ariaDescribedby: { type: String, default: null }
  },
  emits: ['update:modelValue'],
  /**
   * @param {{ modelValue: unknown, fallback: string }} props
   * @param {{ emit: (event: 'update:modelValue', value: string | null) => void }} context
   */
  setup(props, { emit }) {
    const controlAvailable = computed(() => !state.sessionOperationAvailability());
    const mode = computed(() => colorValueMode(props.modelValue));
    const pickerValue = computed(() => (
      toNativeColorInputValue(props.modelValue, props.fallback)
    ));
    const updateMode = (/** @type {{ target: HTMLSelectElement }} */ event) => {
      const busy = state.sessionOperationAvailability();
      if (busy) return busy;
      const nextMode = String(event?.target?.value || 'auto');
      emit(
        'update:modelValue',
        colorValueForMode(nextMode, props.modelValue, props.fallback)
      );
    };
    const updateColor = (/** @type {{ target: HTMLInputElement }} */ event) => {
      const busy = state.sessionOperationAvailability();
      if (busy) return busy;
      emit(
        'update:modelValue',
        toNativeColorInputValue(event?.target?.value, props.fallback)
      );
    };
    return { mode, pickerValue, updateMode, updateColor, controlAvailable };
  },
  template: `
    <div class="grid grid-cols-[minmax(0,1fr)_2.25rem] gap-1 items-center">
      <select
        :value="mode"
        :disabled="!controlAvailable"
        @change="updateMode"
        class="form-input form-input-compact min-w-0"
        :aria-label="\`\${ariaLabel} mode\`"
        :aria-describedby="ariaDescribedby"
      >
        <option value="auto">Auto</option>
        <option v-if="allowNone" value="none">None</option>
        <option value="color">Color</option>
      </select>
      <input
        type="color"
        :value="pickerValue"
        @input="updateColor"
        :disabled="!controlAvailable || mode !== 'color'"
        class="h-8 w-full p-0 border rounded disabled:opacity-40"
        :aria-label="ariaLabel"
        :aria-describedby="ariaDescribedby"
      >
    </div>
  `
};

/**
 * What FileUploader reads of its props. A static `request-clear="true"`
 * attribute arrives as the string `'true'`.
 * @typedef {object} FileUploaderProps
 * @property {{ name?: string } | { name?: string }[] | null} [modelValue] the selected file, or the files of a `multiple` slot
 * @property {boolean | string} [multiple]
 * @property {boolean | string} [requestClear]
 * @property {boolean} [artifactHistory]
 * @property {(file: File | undefined) => unknown} [afterChange] runs after the selection is stored; a result of `false` is not committed
 */

// `inspectionFailed` is the slot's existing inspection status, passed by the
// parent that shows its message (for example "No records were found"): a
// selected file whose inspection failed is not shown as ready (UI-08).
export const FileUploader = {
  template: '#file-uploader-template',
  props: ['label', 'accept', 'modelValue', 'small', 'multiple', 'testId', 'afterChange', 'requestClear', 'artifactHistory', 'inspectionFailed'],
  emits: ['update:modelValue', 'clearRequest'],
  /**
   * @param {FileUploaderProps} props
   * @param {{ emit: (event: 'update:modelValue' | 'clearRequest', value: File | File[] | never[] | EventTarget | null) => void }} context
   */
  setup(props, { emit }) {
    const input = ref(null);
    const uploadAvailable = computed(() => !state.sessionOperationAvailability());
    const selectedFiles = computed(() => {
      if (Array.isArray(props.modelValue)) return props.modelValue.filter(Boolean);
      return props.modelValue ? [props.modelValue] : [];
    });
    const hasSelection = computed(() => selectedFiles.value.length > 0);
    const selectedLabel = computed(() => {
      const items = selectedFiles.value;
      if (items.length === 0) return '';
      if (items.length === 1) return items[0]?.name || 'Selected file';
      const firstNames = items.slice(0, 2).map((/** @type {{ name?: string }} */ file) => file?.name || 'file').join(', ');
      const suffix = items.length > 2 ? ` +${items.length - 2}` : '';
      return `${items.length} files: ${firstNames}${suffix}`;
    });
    const handleFile = (/** @type {{ target: HTMLInputElement }} */ e) => {
      const busy = state.sessionOperationAvailability();
      if (busy) { e.target.value = ''; return busy; }
      const nextFiles = Array.from(e.target.files || []);
      const update = () => {
        if (props.multiple) {
          emit('update:modelValue', nextFiles);
        } else if (nextFiles[0]) {
          emit('update:modelValue', nextFiles[0]);
        }
      };
      const history = window.__GBDRAW_HISTORY__;
      const apply = async () => {
        update();
        await nextTick();
        return await props.afterChange?.(nextFiles[0]);
      };
      if (props.artifactHistory) {
        void history.runUndoableCheckpoint('Change uploaded file', apply, { shouldCommit: (/** @type {unknown} */ result) => result !== false });
      } else if (history?.runUndoable) {
        void history.runUndoable('Change uploaded file', apply);
      } else {
        update();
      }
      e.target.value = '';
    };
    const clearFile = (/** @type {Event} */ event) => {
      const busy = state.sessionOperationAvailability();
      if (busy) return busy;
      if (props.requestClear) {
        emit('clearRequest', event?.currentTarget || null);
        return;
      }
      const history = window.__GBDRAW_HISTORY__;
      if (history?.runUndoable) {
        void history.runUndoable('Change uploaded file', async () => {
          emit('update:modelValue', props.multiple ? [] : null);
          await nextTick();
        });
      } else {
        emit('update:modelValue', props.multiple ? [] : null);
      }
    };
    return { input, handleFile, clearFile, hasSelection, selectedLabel, uploadAvailable };
  }
};

// A record's name in a narrow row: the name shrinks with an ellipsis and its
// detail (the accession) never does, so records of one File stay distinct
// (OV-436). The tooltip and the text hold the whole name.
export const RecordName = {
  props: { prefix: { type: String, default: '' }, name: { type: String, default: '' }, detail: { type: String, default: '' } },
  /** @param {{ prefix: string, name: string, detail: string }} props */
  setup(props) {
    const title = computed(() => [props.prefix, [props.name, props.detail].filter(Boolean).join(' · ')].filter(Boolean).join(' '));
    return { title };
  },
  template: `<span class="record-name" :title="title"><span v-if="prefix" class="record-name-part">{{ prefix + ' ' }}</span><span class="record-name-main">{{ name }}</span><span v-if="detail" class="record-name-part" data-record-name-detail>{{ ' · ' + detail }}</span></span>`
};

/**
 * What RecordDisplayControl calls on `createRecordDisplayControls`' result,
 * which the template passes as `controller`; the row is passed back untouched.
 * @typedef {object} RecordDisplayControlPort
 * @property {(row: unknown, value: string | number) => Promise<unknown>} setStart
 * @property {(row: unknown) => { displayStart?: string | number | null }} surfaceFor
 * @property {(row: unknown) => unknown} draftFor
 */

export const RecordDisplayControl = {
  template: '#record-display-control-template',
  props: ['row', 'controller'],
  /** @param {{ row: unknown, controller: RecordDisplayControlPort }} props */
  setup(props) {
    const error = ref('');
    // The Display start `<input>`: assigning a number sets its decimal text.
    const setStart = async (/** @type {{ target: { value: string | number, checkValidity: () => boolean, validationMessage: string } }} */ event) => {
      const busy = state.sessionOperationAvailability();
      if (busy) return busy;
      try {
        if (!event.target.checkValidity()) throw new Error(event.target.validationMessage);
        await props.controller.setStart(props.row, event.target.value);
        error.value = '';
      } catch (failure) {
        error.value = failure.message;
      }
      // The field shows the model's start after every change: a rejected value
      // goes back, and a cleared field shows the current start (OV-421).
      event.target.value = props.controller.surfaceFor(props.row).displayStart ?? '';
    };
    return { error, setStart, controlAvailable: computed(() => !state.sessionOperationAvailability()), draft: computed(() => props.controller.draftFor(props.row)),
      surface: computed(() => props.controller.surfaceFor(props.row)) };
  }
};

// `compact` shows the error beside its field: the summary and a small Details
// disclosure with the same diagnostics and Copy (PD-OI-046).
export const OperationError = {
  template: '#operation-error-template',
  props: { error: { type: null, default: null }, recovery: { type: null, default: null }, compact: Boolean },
  /** @param {{ error: string | Record<string, unknown> | null | undefined, recovery: string | null | undefined }} props */
  setup(props) {
    const model = computed(() => normalizeUserFacingError(props.error));
    const title = computed(() => operationErrorTitle(model.value?.operation));
    const recoveryText = computed(() => generationRecoveryGuidance(props.recovery));
    const diagnostics = computed(() => (model.value?.details || []).map((/** @type {UserFacingError['details'][number]} */ section) =>
      `${section.label}\n${section.text}`).join('\n\n'));
    const details = ref(null);
    const diagnosticText = ref(null);
    const copyStatus = ref('');
    watch(() => props.error, () => {
      if (details.value) details.value.open = false;
      copyStatus.value = '';
    }, { flush: 'sync' });
    const selectDiagnostics = () => {
      diagnosticText.value?.focus();
      diagnosticText.value?.select();
    };
    const copyDiagnostics = async () => {
      const displayed = props.error;
      const text = diagnostics.value;
      if (!details.value?.open || !text) return;
      try {
        if (!navigator.clipboard?.writeText) throw new Error('Clipboard unavailable');
        await navigator.clipboard.writeText(text);
        if (props.error === displayed) copyStatus.value = 'Diagnostics copied.';
      } catch (_) {
        if (props.error !== displayed) return;
        copyStatus.value = 'Clipboard unavailable. Select diagnostics and copy manually.';
      }
    };
    return { model, title, recoveryText, diagnostics, details, diagnosticText, copyStatus,
      copyDiagnostics, selectDiagnostics };
  }
};

// The Specific color rules rows: a child of the root, so a root update that
// leaves the rules unchanged renders no row (OV-251). The root keeps the
// actions; a row emits them with the arguments it passed before.
export const SpecificRuleRows = {
  template: '#specific-rule-rows-template',
  components: { OperationError, HelpTip },
  props: ['rules', 'featureKeys', 'available', 'patternText', 'patternDraft', 'fieldId'],
  emits: ['set-field', 'move-up', 'move-down', 'remove', 'edit-pattern', 'retry-pattern', 'revert-pattern']
};

// The preview of record or title text with italic markup, shown under its
// field; nothing when the text has no markup.
export const MarkupPreview = {
  props: ['text'],
  template: `<p v-if="segments.length" class="markup-preview" :title="segments.map((part) => part.text).join('')" data-markup-preview><span class="sr-only">Preview: </span><span v-for="(part, index) in segments" :key="index" :class="{ italic: part.italic }">{{ part.text }}</span></p>`,
  /** @param {{ text?: string | null }} props */
  setup(props) {
    return { segments: computed(() => italicSegments(props.text)) };
  }
};

// A Linear record's own name fields. The File card of a one-record File and
// the Record options of every other record render this one component on the
// record's state; `linearFileCardRecord` decides which place shows it.
export const LinearRecordNames = {
  template: '#linear-record-names-template',
  components: { MarkupPreview },
  props: {
    seq: { type: Object, required: true },
    source: { type: Object, required: true },
    idx: { type: Number, required: true },
    disabled: Boolean
  },
  emits: ['reset-definition', 'reset-subtitle'],
  setup() {
    return { getLinearSourceDefaultDefinition, getLinearSourceDefaultSubtitle, resolveLinearRecordEffectiveDefinition };
  }
};
