// @ts-check
// PD-OI-112: the header's History controls. Undo and Redo walk one step each.
// The History list shows the steps they walk, the shown drawing's and the
// project's, newest first, with the current state marked; a pick walks Undo or
// Redo to the state right after the picked step (`travelTo`), and the pick is
// no step. A step reads as its detail, or as its label translated where it is
// shown (`stepText`), in the list and in the Undo and Redo titles. The list
// stays open until a press outside it or Escape. While a pick walks, History
// reports Undo and Redo unavailable, so the rows wait as the buttons do; a
// row that cannot be picked stays focusable (`aria-disabled`).
/** @import { HistoryTimeline, HistoryTimelineStep } from '../services/history.js' */
import { controlValueText } from '../utils/control-description.js';
import { t, tm, tn } from '../utils/i18n.js';

/**
 * @typedef {object} HistoryControlsOptions
 *   The History functions the header's History controls call (R11), and Vue's.
 * @property {() => number} revision A reactive read that changes with every step, Undo, Redo, and baseline.
 * @property {() => boolean} canUndo
 * @property {() => boolean} canRedo
 * @property {() => Promise<unknown>} undo
 * @property {() => Promise<unknown>} redo
 * @property {() => HistoryTimeline} timeline
 * @property {(target: number | null) => Promise<unknown>} travelTo
 * @property {<T>(getter: () => T) => { readonly value: T }} computed Vue `computed`
 * @property {<T>(value: T) => { value: T }} ref Vue `ref`
 * @property {() => Promise<void>} nextTick Vue `nextTick`
 *
 * @typedef {{ name: string, value: string, color: string }} HistoryStepText
 *   A step as shown: the setting's name (or the step's label), and the value it set, with its color for a color.
 * @typedef {{
 *   key: string, target: number | null, text: HistoryStepText,
 *   place: 'undone' | 'current' | 'done' | 'start', direction: 'undo' | 'redo', steps: number
 * }} HistoryListRow
 *   One row of the list: `target` is the step's sequence (`null`: the state before every kept step),
 *   and a pick takes `steps` Undo or Redo steps (`direction`).
 */

const LIST_ID = 'history-step-list';
const LIST_GAP_PX = 4;
const VIEWPORT_MARGIN_PX = 8;

/**
 * Where a list of `width` starts under its toggle: at the toggle's left edge,
 * or so that it ends at the toggle's right edge in a right-to-left page.
 * @param {HTMLElement} toggle @param {DOMRect} box The toggle's box. @param {number} width
 */
const startEdge = (toggle, box, width) => (getComputedStyle(toggle).direction === 'rtl' ? box.right - width : box.left);

/**
 * The one text of a step, for the list rows and the Undo and Redo titles.
 * @param {HistoryTimelineStep} step
 * @returns {HistoryStepText}
 */
const stepText = ({ label, detail }) => ({
  name: detail?.name || tm(label),
  value: detail ? controlValueText(detail) : '',
  color: detail?.color ?? ''
});

/** @param {HistoryControlsOptions} options */
export const createHistoryControls = ({ revision, canUndo, canRedo, undo, redo, timeline: readTimeline, travelTo, computed, ref, nextTick }) => {
  const canUndoHistory = computed(() => {
    void revision();
    return canUndo();
  });
  const canRedoHistory = computed(() => {
    void revision();
    return canRedo();
  });
  const timeline = computed(() => {
    void revision();
    return readTimeline();
  });
  const undoHistoryTitle = computed(() => {
    const step = timeline.value.undo[0];
    return canUndoHistory.value && step ? t('Undo {action}', { action: stepText(step).name }) : t('Undo');
  });
  const redoHistoryTitle = computed(() => {
    const step = timeline.value.redo[0];
    return canRedoHistory.value && step ? t('Redo {action}', { action: stepText(step).name }) : t('Redo');
  });

  /**
   * @param {HistoryTimelineStep} step @param {HistoryListRow['place']} place
   * @param {HistoryListRow['direction']} direction @param {number} steps
   * @returns {HistoryListRow}
   */
  const stepRow = (step, place, direction, steps) => ({
    key: String(step.sequence), target: step.sequence, text: stepText(step), place, direction, steps
  });
  /** @type {{ readonly value: HistoryListRow[] }} */
  const rows = computed(() => {
    const { undo, redo, start } = timeline.value;
    return [
      ...redo.map((step, index) => stepRow(step, 'undone', 'redo', index + 1)).reverse(),
      ...undo.map((step, index) => stepRow(step, index === 0 ? 'current' : 'done', 'undo', index)),
      {
        key: 'start', target: null, direction: 'undo', steps: undo.length, place: undo.length === 0 ? 'current' : 'start',
        text: { name: start.trimmed ? t('Oldest kept state') : tm(start.label), value: '', color: '' }
      }
    ];
  });
  const available = computed(() => canUndoHistory.value || canRedoHistory.value);
  /** @param {HistoryListRow} row */
  const rowDisabled = (row) => row.steps === 0 || !(row.direction === 'undo' ? canUndoHistory.value : canRedoHistory.value);

  const open = ref(false);
  // The row under the pointer or focus, by its key.
  const pointed = ref(/** @type {string | null} */ (null));
  const hint = computed(() => {
    const row = rows.value.find((candidate) => candidate.key === pointed.value);
    if (!row || row.steps === 0) return t('Pick a row to go back or forward to it');
    return row.direction === 'undo'
      ? tn(row.steps, 'Undo {count} step', 'Undo {count} steps')
      : tn(row.steps, 'Redo {count} step', 'Redo {count} steps');
  });
  const place = ref({ top: 0, left: 0 });
  /** @type {HTMLElement | null} */
  let anchor = null;
  const listElement = () => document.getElementById(LIST_ID);

  // Below its toggle and starting where it starts (its right edge in a
  // right-to-left page), inside the window.
  const placeList = () => {
    const list = listElement();
    if (!anchor || !list) return;
    if (!anchor.getClientRects().length) { closeList(); return; }
    const box = anchor.getBoundingClientRect();
    const right = window.innerWidth - list.offsetWidth - VIEWPORT_MARGIN_PX;
    place.value = { top: box.bottom + LIST_GAP_PX, left: Math.max(VIEWPORT_MARGIN_PX, Math.min(startEdge(anchor, box, list.offsetWidth), right)) };
  };
  /** @param {PointerEvent} event */
  const onOutsidePress = (event) => {
    const target = /** @type {Node | null} */ (event.target);
    if (target && (listElement()?.contains(target) || anchor?.contains(target))) return;
    closeList();
  };
  // Escape closes the list alone, not the layer under it as well.
  /** @param {KeyboardEvent} event */
  const onEscape = (event) => {
    if (event.key !== 'Escape') return;
    event.stopPropagation();
    const toggle = anchor;
    closeList();
    toggle?.focus({ preventScroll: true });
  };
  const closeList = () => {
    if (!open.value) return;
    open.value = false;
    pointed.value = null;
    document.removeEventListener('pointerdown', onOutsidePress, true);
    document.removeEventListener('keydown', onEscape, true);
    window.removeEventListener('resize', placeList);
  };
  /** @param {MouseEvent} event */
  const toggleList = async (event) => {
    if (open.value) { closeList(); return; }
    anchor = /** @type {HTMLElement} */ (event.currentTarget);
    open.value = true;
    document.addEventListener('pointerdown', onOutsidePress, true);
    document.addEventListener('keydown', onEscape, true);
    window.addEventListener('resize', placeList);
    await nextTick();
    placeList();
    /** @type {HTMLElement | null | undefined} */ (listElement()?.querySelector('[aria-current="step"]'))?.focus();
  };

  /** @param {HistoryListRow} row */
  const pick = async (row) => {
    if (rowDisabled(row)) return;
    await travelTo(row.target);
  };
  // Arrow keys move through the rows.
  /** @param {KeyboardEvent} event */
  const onListKeydown = (event) => {
    const step = event.key === 'ArrowDown' ? 1 : event.key === 'ArrowUp' ? -1 : 0;
    if (!step) return;
    event.preventDefault();
    const buttons = /** @type {HTMLButtonElement[]} */ ([...(listElement()?.querySelectorAll('.history-list-row') ?? [])]);
    const at = buttons.indexOf(/** @type {HTMLButtonElement} */ (document.activeElement));
    buttons[Math.min(buttons.length - 1, Math.max(0, at + step))]?.focus();
  };

  return {
    canUndoHistory,
    canRedoHistory,
    undoHistoryTitle,
    redoHistoryTitle,
    undoHistory: undo,
    redoHistory: redo,
    historyList: {
      id: LIST_ID,
      available,
      open,
      rows,
      hint,
      place,
      pointed,
      rowDisabled,
      toggleList,
      pick,
      onListKeydown
    }
  };
};
