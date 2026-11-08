// @ts-check
// The Circular drawing's source file and record controls: the Multi-Record
// Canvas checkbox, the source file upload, replace and Remove with the record
// presentation they retire, the record presentation selector, the links that
// reveal the Multi-Record Canvas setting and the Ticks row, and the scroll
// anchoring when the record presentation panel appears. The composition root
// passes the state, the shown-drawing readers and admission, the record
// display's source change, the record discovery, the rebind of annotations,
// the settings pane's reveal, the panels' refs, and the record form defaults
// as ports.
import { bindDrawingInputFile } from '../services/drawing/inputs.js';
import { circularMultiRecordDefaults, updateActiveLayoutPreference } from '../services/layout-preferences.js';
import { settingsScrollOwner } from '../services/settings-sections.js';
/** @import { DiagramMode, DrawingState, SessionBusy, WebState } from '../state.js' */
/** @import { InputFileCandidate, InputFileRef } from '../services/history-files.js' */
/** @import { CircularTrackSlot } from '../services/circular-track-slot-model.js' */
/** @import { LinearSeqWindow } from '../services/drawing/inputs.js' */
/** @import { RevealSetting } from './settings-pane.js' */

const { computed, nextTick, watch } = window.Vue;

/** @typedef {(drawingMode: DiagramMode) => DrawingState | null} ShownDrawingOfMode The shown drawing when it is of that mode. */

/**
 * @typedef {object} CircularMultiRecordCanvasPorts
 * @property {WebState} state
 * @property {() => void} rebindAnnotationRecords Rebinds the shown drawing's annotations to its record catalog (R10).
 */

/**
 * Registers nothing; the root calls it after its `rebindAnnotationRecords`.
 * @param {CircularMultiRecordCanvasPorts} ports
 */
export const createCircularMultiRecordCanvas = ({ state, rebindAnnotationRecords }) => {
  // A Circular drawing whose Multi-Record Canvas is on and has no multi-record
  // sides stored gets the toggle's smart defaults. Load, its rollback, and a
  // History restore install the sides they read or captured (R11, D-15-48).
  /** @param {DrawingState} drawing */
  const storeCircularMultiRecordDefaults = (drawing) => {
    const defaults = drawing.mode === 'circular' && drawing.form.multi_record_canvas
      ? circularMultiRecordDefaults(drawing.layoutPreferences) : null;
    if (defaults) updateActiveLayoutPreference(drawing.layoutPreferences, 'circular', true, defaults);
  };
  // The Multi-Record Canvas checkbox, after its `v-model` wrote the shown
  // drawing: in the control's History step.
  const changeCircularMultiRecordCanvas = () => {
    storeCircularMultiRecordDefaults(state.activeDrawing());
    rebindAnnotationRecords();
  };
  return { changeCircularMultiRecordCanvas };
};

/**
 * A mode panel's mutating handler: the session's busy answer, nothing when the
 * shown drawing is of another mode, else `action` on the shown drawing.
 * @typedef {<T>(drawingMode: DiagramMode, action: (drawing: DrawingState) => T) => T | SessionBusy | undefined} ShownAction
 */
/**
 * The record form fields of a Circular drawing as a new drawing has them.
 * @typedef {Readonly<{ circular_record_selector: string, circular_region_start: number | null,
 *   circular_region_end: number | null, circular_reverse: boolean, circular_windows: readonly LinearSeqWindow[],
 *   circular_record_label: string, circular_record_subtitle: string }>} CircularRecordFormDefaults
 */

/**
 * @typedef {object} CircularSourceEditsPorts
 * @property {WebState} state
 * @property {CircularRecordFormDefaults} CIRCULAR_RECORD_FORM_DEFAULTS
 * @property {<T>(drawing: DrawingState, mutate: () => T) => T} changeRecordSources
 *   Runs a write of the drawing's record sources; the replaced source's display drafts go.
 * @property {(drawing: DrawingState) => { circularUpload: InputFileRef | null }} pendingDiscoveryOf
 *   What the drawing's next record discovery does with a Circular upload.
 * @property {(drawing: DrawingState) => unknown} discoverDrawingRecords Reads the records of the drawing's inputs.
 * @property {ShownAction} shownAction
 */

/**
 * Registers nothing.
 * @param {CircularSourceEditsPorts} ports
 */
export const createCircularSourceEdits = ({
  state, CIRCULAR_RECORD_FORM_DEFAULTS, changeRecordSources, pendingDiscoveryOf, discoverDrawingRecords, shownAction
}) => {
  // CI-03: the per-record presentation a source file carries in each mode.
  // Replacing or removing the file resets it in the replacement's History step.
  // The OFF records name records of the old file too (record selection).
  /** @param {DrawingState} drawing */
  const retireCircularRecordSelector = ({ form, adv, recordsOff }) => {
    form.circular_record_selector = CIRCULAR_RECORD_FORM_DEFAULTS.circular_record_selector;
    if (adv.circular_grouping_intent === 'single') adv.circular_grouping_intent = 'auto';
    if (recordsOff.length) recordsOff.splice(0);
  };
  // The whole record presentation: the record form fields return to their
  // defaults (the one list a style copy leaves out too). Each field is written
  // by name, so the track-layout writer guard reads no whole-form write; a
  // node case pins that every field of the list is written.
  /** @param {DrawingState} drawing */
  const retireCircularRecordPresentation = (drawing) => {
    const { form } = drawing;
    retireCircularRecordSelector(drawing);
    form.circular_region_start = CIRCULAR_RECORD_FORM_DEFAULTS.circular_region_start;
    form.circular_region_end = CIRCULAR_RECORD_FORM_DEFAULTS.circular_region_end;
    form.circular_reverse = CIRCULAR_RECORD_FORM_DEFAULTS.circular_reverse;
    form.circular_windows = [...CIRCULAR_RECORD_FORM_DEFAULTS.circular_windows];
    form.circular_record_label = CIRCULAR_RECORD_FORM_DEFAULTS.circular_record_label;
    form.circular_record_subtitle = CIRCULAR_RECORD_FORM_DEFAULTS.circular_record_subtitle;
  };
  // D-32: one record replacing one record (a new version of the same genome)
  // keeps its crop and titles; a removed file or a replaced multi-record file
  // retires them. The record selector names a record of the old file, so any
  // replacement retires it.
  /** @param {{ removed: boolean, previousRecordCount: number }} replacement */
  const sourceReplacementRetiresPresentation = ({ removed, previousRecordCount }) => removed || previousRecordCount > 1;
  /**
   * A Circular drawing's source file (upload, replace, Remove), with the
   * drawing's reactions to its new source: the shown drawing's controls and
   * the replacement for every drawing that uses a file write through it.
   * @param {DrawingState} drawing
   * @param {'c_gb' | 'c_gff' | 'c_fasta'} field
   * @param {InputFileCandidate | InputFileRef} value
   */
  const replaceCircularSource = (drawing, field, value) => {
    const nextValue = value ?? null;
    const previous = drawing.files[field];
    const nextRef = state.inputFiles.registerFile(nextValue);
    if (previous === nextRef) return;
    const previousRecordCount = drawing.circularRecordList.value.length;
    changeRecordSources(drawing, () => { bindDrawingInputFile(state.inputFiles, drawing.slices.inputs, field, nextValue); });
    if (nextRef) pendingDiscoveryOf(drawing).circularUpload = nextRef;
    void discoverDrawingRecords(drawing);
    // An upload into an empty slot starts reading another input (after an
    // Input switch), whose records the OFF list does not name.
    if (!previous) {
      if (drawing.recordsOff.length) drawing.recordsOff.splice(0);
      return;
    }
    if (sourceReplacementRetiresPresentation({
      removed: !nextValue, previousRecordCount
    })) retireCircularRecordPresentation(drawing);
    else retireCircularRecordSelector(drawing);
  };
  /**
   * The Circular source file controls of the shown drawing.
   * @param {'c_gb' | 'c_gff' | 'c_fasta'} field
   * @param {File | null} value
   */
  const setCircularSourceFile = (field, value) => shownAction('circular', (drawing) => replaceCircularSource(drawing, field, value));
  return { sourceReplacementRetiresPresentation, replaceCircularSource, setCircularSourceFile };
};

/**
 * @typedef {object} CircularRecordDisplayPorts
 * @property {WebState} state
 * @property {{ readonly value: boolean }} semanticMutationAvailable No session operation holds the inputs.
 * @property {{ readonly value: { status: string, primaryFile: InputFileRef | null, pairedFile: InputFileRef | null } }} circularRecordDiscoveryState
 *   The shown drawing's Circular record discovery.
 * @property {{ readonly value: boolean }} circularSingleRecordPresentationEnabled
 *   Whether the shown Circular drawing presents one record.
 * @property {() => SessionBusy | null} sessionOperationAvailability
 * @property {ShownDrawingOfMode} shownDrawingOfMode
 * @property {() => void} rebindAnnotationRecords Rebinds the shown drawing's annotations to its record catalog (R10).
 * @property {RevealSetting} revealSetting Brings a setting into view.
 * @property {{ value: boolean }} circularTrackSlotsPanelOpen Whether the Custom Track Slots panel is open.
 * @property {{ readonly value: HTMLDetailsElement | null }} circularRecordPresentationPanel The record presentation panel.
 * @property {() => Promise<unknown>} afterPaint Settles after the next paint.
 */

/**
 * Registers the `circularRecordSelectionEnabled` computed and the scroll
 * anchoring watch; the root calls it where they stood.
 * @param {CircularRecordDisplayPorts} ports
 */
export const createCircularRecordDisplay = ({
  state, semanticMutationAvailable, circularRecordDiscoveryState, circularSingleRecordPresentationEnabled,
  sessionOperationAvailability, shownDrawingOfMode, rebindAnnotationRecords, revealSetting, circularTrackSlotsPanelOpen,
  circularRecordPresentationPanel, afterPaint
}) => {
  const { mode } = state;
  const circularRecordSelectionEnabled = computed(() => semanticMutationAvailable.value
    && mode.value === 'circular'
    && !state.activeDrawing().form.multi_record_canvas && circularRecordDiscoveryState.value.status === 'ready');
  /** @param {unknown} value */
  const setCircularRecordPresentationSelector = (value) => {
    const busy = sessionOperationAvailability();
    if (busy) return busy;
    const drawing = shownDrawingOfMode('circular');
    if (!drawing) return { status: 'unavailable' };
    if (!circularRecordSelectionEnabled.value) return { status: 'unavailable' };
    const normalized = String(value || '').trim();
    drawing.form.circular_record_selector = normalized;
    drawing.adv.circular_grouping_intent = normalized ? 'single' : 'auto';
    rebindAnnotationRecords();
    return { status: 'ok' };
  };
  const showCircularCanvasSetting = async () => {
    const drawing = shownDrawingOfMode('circular');
    const control = /** @type {HTMLElement | null} */ (document.querySelector('[data-circular-canvas-setting]'));
    if (!control || !drawing || !drawing.form.multi_record_canvas) {
      return { status: 'unavailable' };
    }
    await revealSetting(control);
    return { status: 'ok' };
  };
  // With Custom Track Slots on, a Ticks row owns the coordinate scale: the link beside Show Coordinate Scale
  // reveals that row (the first enabled one, else the first), or the Custom Track Slots panel when the stack
  // has none. It changes no setting.
  const showCircularTicksSlot = async () => {
    const drawing = state.activeDrawing();
    if (mode.value !== 'circular' || !drawing.adv.circular_track_slots_enabled) return { status: 'unavailable' };
    circularTrackSlotsPanelOpen.value = true;
    await nextTick();
    /** @type {CircularTrackSlot[]} */
    const slots = Array.isArray(drawing.adv.circular_track_slots) ? drawing.adv.circular_track_slots : [];
    const ticks = slots.filter((slot) => slot.renderer === 'ticks');
    const slot = ticks.find((candidate) => candidate.enabled) || ticks[0];
    const row = slot ? document.querySelector(`[data-capture="circular-track-slot-${CSS.escape(slot.id)}"]`) : null;
    const panel = document.getElementById('circular-custom-track-slots-panel');
    const toggle = /** @type {HTMLElement | null} */ (document.querySelector('[data-track-slots-toggle]'));
    const control = /** @type {HTMLElement | null} */ (row?.querySelector('input, select, button') ?? null);
    if (row instanceof HTMLElement && control) {
      await revealSetting(row, { focus: control });
    } else if (panel && toggle) {
      await revealSetting(panel, { align: 'start', focus: toggle });
    } else {
      return { status: 'unavailable' };
    }
    return { status: 'ok' };
  };
  watch(() => [mode.value, circularSingleRecordPresentationEnabled.value,
    circularRecordDiscoveryState.value.primaryFile, circularRecordDiscoveryState.value.pairedFile,
    state.activeDrawing().form.circular_record_selector], async (/** @type {unknown[]} */ current, /** @type {unknown[]} */ previous = []) => {
    if (!current[1] || current.every((value, index) => Object.is(value, previous[index]))) return;
    const origin = document.activeElement;
    const pane = document.querySelector('.settings-scroll');
    const scrollOwner = pane ? settingsScrollOwner(pane) : document.scrollingElement;
    const anchored = origin && pane?.contains(origin) && origin.getClientRects().length;
    const top = anchored ? origin.getBoundingClientRect().top : null;
    const scrollTop = scrollOwner?.scrollTop;
    await nextTick();
    if (!circularSingleRecordPresentationEnabled.value
      || circularRecordDiscoveryState.value.primaryFile !== current[2]
      || circularRecordDiscoveryState.value.pairedFile !== current[3]) return;
    if (circularRecordPresentationPanel.value) circularRecordPresentationPanel.value.open = true;
    await afterPaint();
    if (scrollOwner && document.activeElement === origin) {
      // `top` is a number whenever `anchored` is truthy; `scrollTop` was read from
      // the same non-null `scrollOwner`.
      scrollOwner.scrollTop = anchored && origin.isConnected
        ? scrollOwner.scrollTop + origin.getBoundingClientRect().top - /** @type {number} */ (top)
        : /** @type {number} */ (scrollTop);
    }
  }, { flush: 'sync' });
  return {
    circularRecordSelectionEnabled, setCircularRecordPresentationSelector, showCircularCanvasSetting,
    showCircularTicksSlot
  };
};
