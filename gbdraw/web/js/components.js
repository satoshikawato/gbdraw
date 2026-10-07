// @ts-check
import { state } from './state.js';
import {
  colorValueForMode,
  colorValueMode,
  toNativeColorInputValue
} from './utils/color-utils.js';

import { normalizeUserFacingError, operationErrorTitle, generationRecoveryGuidance } from './utils/error-normalization.js';

const { ref, reactive, computed, nextTick, watch, useId, onMounted } = window.Vue;

/** @typedef {{ element: HTMLElement, returnFocus: HTMLElement | null }} FocusScope */
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
    focusScopes.push({ element, returnFocus: active instanceof HTMLElement && active !== document.body ? active : null });
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
// binds Tab to `trapDialogFocus`.
export const ChoiceDialog = {
  template: '#choice-dialog-template',
  props: {
    labelledby: { type: String, required: true },
    describedby: { type: String, default: null }
  },
  emits: ['cancel']
};

/** @param {Node} node @returns {string} */
const visibleText = (node) => {
  if (node.nodeType === Node.TEXT_NODE) return node.textContent || '';
  if (!(node instanceof Element) || node.matches('button, input, select, textarea, [hidden], [aria-hidden="true"]')) return '';
  return Array.from(node.childNodes, visibleText).join(' ');
};

// R12 (UI-04): a tip is named after what it explains, "Help: <label>": the
// author name or visible text of the element that holds it, after the name of
// the row group it is in. `label` names a tip whose holder has no such text.
/** @param {Element | null} tip */
const helpTipName = (tip) => {
  const holder = tip?.parentElement;
  const group = holder?.closest('[role="group"][aria-label]')?.getAttribute('aria-label') || '';
  const label = holder?.getAttribute('aria-label') || (holder ? visibleText(holder) : '');
  const name = `${group} ${label}`.replace(/\s+/g, ' ').trim();
  return name ? `Help: ${name}` : 'Help';
};

// Every tip is a disclosure button (PD-OI-057). One text source feeds the
// visible tooltip and the accessible description, identified by `id` or an
// automatic id; the control that a tip explains references that id through
// aria-describedby. The description is a hidden element teleported to <body>,
// so the text around the tip stays its visible label and screen readers do
// not read the tip text twice. Mouse hover and keyboard focus show the
// tooltip, a click or tap toggles it, and Escape or blur closes it. Place
// tips outside <label> so they never change the accessible name of a label or
// its control.
export const HelpTip = {
  template: '#help-tip-template',
  props: ['text', 'id', 'label'],
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
    watch(visible, (open) => { if (open) position(); });
    // A touch tap also reports pointerenter; only a hovering pointer shows it.
    const onEnter = (event) => { hovered.value = event.pointerType !== 'touch'; };
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

export const AutoValueField = {
  props: ['visible', 'text'],
  template: `
    <div class="auto-value-field">
      <span v-if="visible" class="auto-value-placeholder">{{ text }}</span>
      <slot></slot>
    </div>
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
  setup(props, { emit }) {
    const controlAvailable = computed(() => !state.sessionOperationAvailability());
    const mode = computed(() => colorValueMode(props.modelValue));
    const pickerValue = computed(() => (
      toNativeColorInputValue(props.modelValue, props.fallback)
    ));
    const updateMode = (event) => {
      const busy = state.sessionOperationAvailability();
      if (busy) return busy;
      const nextMode = String(event?.target?.value || 'auto');
      emit(
        'update:modelValue',
        colorValueForMode(nextMode, props.modelValue, props.fallback)
      );
    };
    const updateColor = (event) => {
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

export const FileUploader = {
  template: '#file-uploader-template',
  props: ['label', 'accept', 'modelValue', 'small', 'multiple', 'testId', 'afterChange', 'requestClear', 'artifactHistory'],
  emits: ['update:modelValue', 'clearRequest'],
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
      const firstNames = items.slice(0, 2).map((file) => file?.name || 'file').join(', ');
      const suffix = items.length > 2 ? ` +${items.length - 2}` : '';
      return `${items.length} files: ${firstNames}${suffix}`;
    });
    const handleFile = (e) => {
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
        void history.runUndoableCheckpoint('Change uploaded file', apply, { shouldCommit: result => result !== false });
      } else if (history?.runUndoable) {
        void history.runUndoable('Change uploaded file', apply);
      } else {
        update();
      }
      e.target.value = '';
    };
    const clearFile = (event) => {
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

export const RecordDisplayControl = {
  template: '#record-display-control-template',
  props: ['row', 'controller'],
  setup(props) {
    const error = ref('');
    const setStart = async (event) => {
      const busy = state.sessionOperationAvailability();
      if (busy) return busy;
      try {
        if (!event.target.checkValidity()) throw new Error(event.target.validationMessage);
        await props.controller.setStart(props.row, event.target.value);
        error.value = '';
      } catch (failure) {
        error.value = failure.message;
        event.target.value = props.controller.draftFor(props.row).startCoordinate ?? '';
      }
    };
    return { error, setStart, controlAvailable: computed(() => !state.sessionOperationAvailability()), draft: computed(() => props.controller.draftFor(props.row)),
      surface: computed(() => props.controller.surfaceFor(props.row)) };
  }
};

export const OperationError = {
  template: '#operation-error-template',
  props: ['error', 'recovery'],
  setup(props) {
    const model = computed(() => normalizeUserFacingError(props.error));
    const title = computed(() => operationErrorTitle(model.value?.operation));
    const recoveryText = computed(() => generationRecoveryGuidance(props.recovery));
    const diagnostics = computed(() => (model.value?.details || []).map(section =>
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
