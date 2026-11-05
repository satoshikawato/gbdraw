// @ts-check
// Select region (PD-OI-115): the preview toolbar toggle, the drag along a
// record's axis, the band over the span, and the actions of the Region popup.
// The span comes only from the axis Python wrote for the shown Result
// (`services/region-axis.js`); the selection lives in the open Region popup
// (the one Result popup, PD-OI-105), so it is no History step, is not saved,
// and is not exported.
/**
 * @import { DisplaySpan, RegionAxis, RegionAxisRecord, RegionFeatureRow, RegionSelection } from '../services/region-axis.js'
 * @import { PopupBounds } from '../utils/popup-bounds.js'
 */
import { FEATURE_LABEL_SELECTOR, LEGEND_GROUP_SELECTOR, RECORD_DEFINITION_SELECTOR } from '../services/drawn-groups.js';
import {
  axisPosition, circularDragSpan, featureParts, linearDragSpan, linearRowPosition, partsInsideSpan, regionAxisFrame,
  regionBandPath, regionSpanText, sourceSpanOf
} from '../services/region-axis.js';
import { copyPendingTextToClipboard } from '../utils/clipboard.js';
import { normalizeCaughtError } from '../utils/error-normalization.js';
import { errorText } from '../utils/error-text.js';
import { capturePointer, holdTextSelection } from '../utils/pointer-gesture.js';
import { isClickAtPress } from '../utils/popup-bounds.js';
import { computed, reactive, ref, shallowRef } from '../utils/reactivity.js';
import { t, tn } from '../utils/i18n.js';

/** How far beside the feature tracks, in screen pixels, a press still starts a drag. */
const BAND_SLACK_PX = 6;

/**
 * The region an action works on: the record by its request key, in source coordinates (1-based, both ends included;
 * `start > end` across the origin of a circular record).
 * @typedef {object} RegionTarget
 * @property {string} recordKey
 * @property {number} start
 * @property {number} end
 */

/**
 * @typedef {'copy' | 'add'} RegionAction An action of the Region popup that waits on Python.
 * @typedef {object} RegionSelectPorts
 * @property {() => RegionAxisRecord[]} readAxisRecords The records of the shown Result whose axis Python wrote.
 * @property {() => SVGSVGElement | null} shownRoot The mounted root of the shown Result.
 * @property {{ value: boolean }} layoutEdit Layout edit; this owner is the one that turns it on and off.
 * @property {() => RegionSelection | null} openSelection The selection the open Region popup holds, else null.
 * @property {(selection: RegionSelection, element: Element, anchor: RegionPopupAnchor) => boolean} openRegionPopup
 *   Opens the Region popup (the one Result popup) beside the band, clear of the anchor's obstacles; false when popups
 *   cannot open now.
 * @property {() => void} closeRegionPopup
 * @property {() => readonly RegionFeatureRow[]} drawnFeatures The features the shown Result draws.
 * @property {(featureIds: readonly string[]) => void} selectFeatures The feature selection owner's selection by id.
 * @property {(region: RegionTarget) => Promise<unknown>} addRegionAnnotation The annotation owner's "Region N" row,
 *   bound by the record catalog's key; it reads the drawing's records first when a Load left them unread, and adds
 *   nothing (it rejects) when they cannot be read.
 * @property {(region: RegionTarget) => Promise<unknown>} newDrawingFromRegion The region drawing owner's derivation.
 * @property {(region: RegionTarget) => Promise<string>} readRegionFasta The FASTA text of the region, which Python
 *   reads from the request the shown Result was made from (a Worker helper, started by Copy FASTA only).
 * @property {PanPort} panStart The preview's pan starts at a press (the pan owner's, which pans on mouse events).
 * @property {PanPort} panMove The pan follows the pointer.
 * @property {PanPort} panEnd The pan ends.
 */

/**
 * The preview pan owner's handler of one pan event; a pointer event is a mouse event with the same position.
 * @typedef {(event: MouseEvent & { target: Element }) => void} PanPort
 */

/**
 * Where the Region popup opens: the release point, the boxes it keeps clear of besides the band (D-28), and the
 * rectangle it stays inside (above the preview toolbar, D-30).
 * @typedef {{ clientX: number, clientY: number, obstacles: PopupBounds[], within: PopupBounds | null }} RegionPopupAnchor
 */

// The preview's own controls (the page's markup): the export bar and the search bar are obstacles like the Legend;
// the toolbar of the preview, where Select region is turned off, is a limit the popup stays above (D-30).
const PREVIEW_BARS_SELECTOR = '.preview-toolbar, .preview-feature-search';
const PREVIEW_CONTROLS_SELECTOR = '.preview-controls';

/** @param {Element} element @returns {PopupBounds | null} The element's screen box; null when it is not shown. */
const shownBox = (element) => {
  const { left, right, top, bottom, width, height } = element.getBoundingClientRect();
  return width > 0 && height > 0 ? { left, right, top, bottom } : null;
};

/**
 * Where the Region popup may go: clear of the shown Result's Legend, feature labels and record definitions and of
 * the export and search bars, and above the preview toolbar (an obstacle instead when the toolbar is in the upper
 * half of the window).
 * @param {SVGSVGElement | null} root
 * @returns {{ obstacles: PopupBounds[], within: PopupBounds | null }}
 */
const popupLimits = (root) => {
  const drawn = root ? Array.from(root.querySelectorAll(`${LEGEND_GROUP_SELECTOR}, ${FEATURE_LABEL_SELECTOR}, ${RECORD_DEFINITION_SELECTOR}`)) : [];
  const obstacles = [...drawn, ...document.querySelectorAll(PREVIEW_BARS_SELECTOR)].flatMap((element) => shownBox(element) ?? []);
  const controls = document.querySelector(PREVIEW_CONTROLS_SELECTOR);
  const toolbar = controls ? shownBox(controls) : null;
  if (toolbar && toolbar.top <= window.innerHeight / 2) obstacles.push(toolbar);
  const within = toolbar && toolbar.top > window.innerHeight / 2
    ? { top: -Infinity, left: -Infinity, right: Infinity, bottom: toolbar.top }
    : null;
  return { obstacles, within };
};

/** @param {RegionSelection} selection @returns {RegionTarget} */
const targetOf = (selection) => ({ recordKey: selection.recordKey, start: selection.source.start, end: selection.source.end });

/** @param {RegionAxis} axis @param {number} start @param {number} current @param {number} sweep */
const dragSpan = (axis, start, current, sweep) => (axis.kind === 'circular'
  ? circularDragSpan(axis, start, sweep)
  : linearDragSpan(start, current));

/**
 * @param {RegionSelectPorts} ports
 */
export const createRegionSelect = ({
  readAxisRecords, shownRoot, layoutEdit, openSelection, openRegionPopup, closeRegionPopup, drawnFeatures,
  selectFeatures, addRegionAnnotation, newDrawingFromRegion, readRegionFasta, panStart, panMove, panEnd
}) => {
  const on = ref(false);
  /** Whether the shown Result has an axis to select on (one predicate for the toggle and the gesture). */
  const available = computed(() => readAxisRecords().length > 0);
  const active = computed(() => on.value && available.value);
  /** Why the toggle is unavailable: a shown Result without the axis (saved before Python wrote it, Q-A1). */
  const unavailableReason = computed(() => (available.value ? '' : t('Generate to enable region selection')));

  /** Turning Select region on turns Layout edit off; turning it off cancels a drag and closes the popup. */
  const toggle = () => {
    if (!on.value && !available.value) return;
    on.value = !on.value;
    if (on.value) layoutEdit.value = false;
    else turnOff();
  };
  /** Turning Layout edit on turns Select region off. */
  const toggleLayoutEdit = () => {
    layoutEdit.value = !layoutEdit.value;
    if (layoutEdit.value && on.value) {
      on.value = false;
      turnOff();
    }
  };
  const turnOff = () => {
    cancelDrag();
    closeRegionPopup();
  };

  /**
   * The drag in progress: its record, where it started along the axis, and the span so far.
   * @type {{ pointerId: number, startX: number, startY: number, moved: boolean, record: RegionAxisRecord,
   *   frame: SVGGraphicsElement, start: number, previous: number, sweep: number, row: number, span: DisplaySpan,
   *   clientX: number, clientY: number, owner: Element } | null}
   */
  let pending = null;
  const drag = reactive({ active: false, recordIndex: -1, span: /** @type {DisplaySpan} */ ({ start: 0, length: 0 }), x: 0, y: 0 });

  /** @param {SVGGraphicsElement} frame @param {number} clientX @param {number} clientY */
  const framePoint = (frame, clientX, clientY) => {
    const matrix = frame.getScreenCTM();
    if (!matrix) return null;
    const point = new DOMPoint(clientX, clientY).matrixTransform(matrix.inverse());
    return { point: { x: point.x, y: point.y }, scale: Math.hypot(matrix.a, matrix.b) || 1 };
  };

  /**
   * The record whose feature tracks are under a screen point.
   * @param {number} clientX
   * @param {number} clientY
   */
  const recordAt = (clientX, clientY) => {
    const root = shownRoot();
    if (!root) return null;
    for (const record of readAxisRecords()) {
      const frame = regionAxisFrame(root, record.recordIndex);
      const mapped = frame && framePoint(frame, clientX, clientY);
      if (!frame || !mapped) continue;
      const position = axisPosition(record.axis, mapped.point, BAND_SLACK_PX / mapped.scale);
      if (position.inBand) return { record, frame, position };
    }
    return null;
  };

  /** Ends the hold a drag takes on the pointer: on text selection from the press, and its capture once it moves. */
  let releaseGesture = () => {};

  /** The touch that pans the preview while the tool is on, by its pointer id; -1 when none does. */
  let panningTouch = -1;

  /**
   * A press with Select region on: on a record's feature tracks it starts a drag (and no pan); anywhere else, or
   * with Shift, Ctrl, or Meta held, it keeps its meaning. A touch drags the same way: the canvas takes no native touch
   * scroll while the tool is on (`.region-select-on`), so a touch elsewhere pans through the pan owner.
   * @param {PointerEvent} event
   */
  const pointerDown = (event) => {
    if (!active.value || event.button !== 0 || event.isPrimary === false || !(event.currentTarget instanceof Element)) return;
    if (event.shiftKey || event.ctrlKey || event.metaKey || event.altKey) return;
    const hit = recordAt(event.clientX, event.clientY);
    if (!hit) {
      if (event.pointerType === 'touch' && event.target instanceof Element) {
        panningTouch = event.pointerId;
        panStart(/** @type {PointerEvent & { target: Element }} */ (event));
      }
      return;
    }
    // A press on the tracks is the drag's: no pan starts and no text is selected.
    event.preventDefault();
    pending = {
      pointerId: event.pointerId, startX: event.clientX, startY: event.clientY, moved: false,
      record: hit.record, frame: hit.frame, start: hit.position.position, previous: hit.position.position, sweep: 0,
      row: hit.position.row, span: { start: Math.floor(hit.position.position), length: 1 },
      clientX: event.clientX, clientY: event.clientY, owner: event.currentTarget
    };
    releaseGesture = holdTextSelection();
    // Esc and leaving the window cancel the drag, as they cancel the marquee; held only while a press is held.
    window.addEventListener('keydown', cancelOnEscape, true);
    window.addEventListener('blur', cancelDrag);
  };

  /**
   * Whether the drag can still select: the tool is on and the record it started on is still drawn (a drawing switch,
   * a Result switch, a Generate, or a Load mounts another Result, and the frame leaves the page).
   * @param {{ frame: Element }} drag
   */
  const dragLives = (drag) => active.value && drag.frame.isConnected;

  /** @param {PointerEvent} event */
  const pointerMove = (event) => {
    if (event.pointerId === panningTouch) panMove(/** @type {PointerEvent & { target: Element }} */ (event));
    if (!pending || event.pointerId !== pending.pointerId) return;
    if (!dragLives(pending)) {
      cancelDrag();
      return;
    }
    // A press that stays within the click's slop is the click's (a feature popup), not a drag.
    if (!pending.moved && isClickAtPress({ clientX: pending.startX, clientY: pending.startY }, event)) return;
    const mapped = framePoint(pending.frame, event.clientX, event.clientY);
    if (!mapped) return;
    const { axis } = pending.record;
    let current;
    if (axis.kind === 'circular') {
      current = axisPosition(axis, mapped.point).position;
      // The sweep follows the pointer around the circle, so a drag passes the origin.
      const turn = axis.source.drawnLength;
      const step = ((current - pending.previous + turn * 1.5) % turn) - turn / 2;
      pending.sweep += step;
      pending.previous = current;
    } else {
      current = linearRowPosition(axis, pending.row, mapped.point.x).position;
    }
    if (!pending.moved) {
      // From here the press is a drag: the canvas captures the pointer, so a release over the toolbar or outside the
      // window still ends it. A press that stays put keeps no capture, so its click reaches what it was on.
      pending.moved = true;
      const releaseCapture = capturePointer(pending.owner, pending.pointerId);
      const releaseText = releaseGesture;
      releaseGesture = () => {
        releaseCapture();
        releaseText();
      };
      closeRegionPopup();
    }
    pending.span = dragSpan(axis, pending.start, current, pending.sweep);
    pending.clientX = event.clientX;
    pending.clientY = event.clientY;
    Object.assign(drag, { active: true, recordIndex: pending.record.recordIndex, span: pending.span, x: event.clientX, y: event.clientY });
    event.preventDefault();
    event.stopPropagation();
  };

  const endDrag = () => {
    pending = null;
    drag.active = false;
    window.removeEventListener('keydown', cancelOnEscape, true);
    window.removeEventListener('blur', cancelDrag);
    const release = releaseGesture;
    releaseGesture = () => {};
    release();
  };
  /** Ends a drag in progress with no selection and no popup. */
  const cancelDrag = () => {
    if (pending) endDrag();
  };
  /**
   * Esc during a drag cancels it, and is the drag's alone (the marquee's Esc does not also clear the selection).
   * @param {KeyboardEvent} event
   */
  const cancelOnEscape = (event) => {
    if (event.key !== 'Escape' || !pending) return;
    cancelDrag();
    event.preventDefault();
    event.stopPropagation();
  };

  /** @param {PointerEvent} event The release or cancel of the touch that pans ends its pan. */
  const endTouchPan = (event) => {
    if (event.pointerId !== panningTouch) return;
    panningTouch = -1;
    panEnd(/** @type {PointerEvent & { target: Element }} */ (event));
  };

  /** The band element the popup opens beside; bound by the template. @type {Element | null} */
  let bandElement = null;
  /** @param {Element | null} element */
  const bindBand = (element) => { bandElement = element; };

  /**
   * The release: a drag selects its span and opens the Region popup beside the band; a press that did not move is
   * left to the click.
   * @param {PointerEvent} event
   */
  const pointerUp = (event) => {
    endTouchPan(event);
    if (!pending || event.pointerId !== pending.pointerId) return;
    const finished = pending;
    endDrag();
    if (!finished.moved || !dragLives(finished)) return;
    event.preventDefault();
    event.stopPropagation();
    const { record } = finished;
    /** @type {RegionSelection} */
    const selection = {
      kind: 'region', resultIndex: record.resultIndex, recordIndex: record.recordIndex, recordKey: record.recordKey,
      recordId: record.recordId, span: finished.span, source: sourceSpanOf(record.axis.source, finished.span)
    };
    shownSelection.value = selection;
    window.requestAnimationFrame(() => {
      if (shownSelection.value !== selection) return;
      // Once open, the popup holds the selection: it goes when the popup closes.
      if (bandElement) {
        openRegionPopup(selection, bandElement, { clientX: finished.clientX, clientY: finished.clientY, ...popupLimits(shownRoot()) });
      }
      shownSelection.value = null;
    });
  };

  /**
   * A cancelled pointer (a touch the browser took), or a capture the drag lost, ends the drag with no selection. Only
   * the canvas's own capture counts: a touch press is captured by the element under it, and that capture moves to
   * the canvas when the drag starts (OV-2483).
   * @param {PointerEvent} event
   */
  const pointerCancel = (event) => {
    // The canvas's own cancel and capture-loss handlers end a pan (the pan owner's `endPan`).
    if (event.pointerId === panningTouch) panningTouch = -1;
    if (!pending || event.pointerId !== pending.pointerId) return;
    if (event.type === 'lostpointercapture' && event.target !== pending.owner) return;
    endDrag();
  };

  /**
   * The selection whose band shows between the release and the popup's open (shallow: the release compares it).
   * @type {{ value: RegionSelection | null }}
   */
  const shownSelection = shallowRef(null);
  /** The selection the band shows: the open popup's, or the one being opened. */
  const selection = computed(() => openSelection() ?? shownSelection.value);

  /** The record of a selection in the shown Result, found by its key; null once the Result no longer draws it. */
  /** @param {RegionSelection} chosen */
  const recordOf = (chosen) => readAxisRecords().find((record) => (
    record.resultIndex === chosen.resultIndex && record.recordKey === chosen.recordKey
  )) ?? null;

  // The band is drawn over the Result in the overlay's own SVG (never in the Result), placed by the record group's
  // screen matrix, which a pan, a zoom, or a Layout edit move changes: it is read again on each frame while shown.
  const overlayMatrix = ref('');
  /** @type {Element | null} */
  let overlayElement = null;
  let frameRequest = 0;
  /** @param {Element | null} element */
  const bindOverlay = (element) => {
    overlayElement = element;
    if (element && !frameRequest) frameRequest = window.requestAnimationFrame(followFrame);
  };
  const followFrame = () => {
    frameRequest = 0;
    const root = shownRoot();
    const recordIndex = drag.active ? drag.recordIndex : (selection.value?.recordIndex ?? -1);
    const frame = root && recordIndex >= 0 ? regionAxisFrame(root, recordIndex) : null;
    const matrix = frame?.getScreenCTM();
    if (!overlayElement || !matrix) {
      overlayMatrix.value = '';
      return;
    }
    const box = overlayElement.getBoundingClientRect();
    const next = `matrix(${matrix.a} ${matrix.b} ${matrix.c} ${matrix.d} ${matrix.e - box.left} ${matrix.f - box.top})`;
    if (next !== overlayMatrix.value) overlayMatrix.value = next;
    frameRequest = window.requestAnimationFrame(followFrame);
  };

  /** The band and, while dragging, the readout "start-end · N bp". */
  const overlay = computed(() => {
    if (drag.active && pending) {
      const source = sourceSpanOf(pending.record.axis.source, drag.span);
      return {
        path: regionBandPath(pending.record.axis, drag.span),
        readout: t('{start}-{end} · {length} bp', regionSpanText(source)),
        readoutX: drag.x,
        readoutY: drag.y
      };
    }
    const chosen = selection.value;
    const record = chosen && recordOf(chosen);
    if (!chosen || !record) return null;
    return { path: regionBandPath(record.axis, chosen.span), readout: '', readoutX: 0, readoutY: 0 };
  });

  /**
   * The drawn features fully inside the selection, by feature id.
   * @param {RegionSelection} chosen
   * @returns {string[]}
   */
  const featuresInside = (chosen) => {
    const record = recordOf(chosen);
    if (!record) return [];
    return drawnFeatures().flatMap((feature) => {
      const id = typeof feature.svg_id === 'string' ? feature.svg_id : '';
      if (!id || feature.record_key !== chosen.recordKey) return [];
      return partsInsideSpan(record.axis.source, chosen.span, featureParts(feature)) ? [id] : [];
    });
  };

  /**
   * The open selection's action that waits on Python (Copy FASTA, Add to annotation set): busy while it runs, then its
   * outcome text or the failure's.
   * @type {{ selection: RegionSelection | null, name: RegionAction | '', busy: boolean, status: string, failed: boolean }}
   */
  const action = reactive({ selection: null, name: '', busy: false, status: '', failed: false });
  /** @param {RegionSelection} chosen @param {RegionAction} name The state of one action of this selection. */
  const actionView = (chosen, name) => {
    const own = action.selection === chosen && action.name === name;
    return { busy: own && action.busy, status: own ? action.status : '', failed: own && action.failed };
  };

  /** The Region popup's body: the feature count and the actions. */
  const popup = computed(() => {
    const chosen = openSelection();
    if (!chosen) return null;
    const featureIds = featuresInside(chosen);
    const featureText = tn(featureIds.length, '{count} feature', '{count} features');
    return {
      featureCount: featureIds.length,
      featureText: chosen.source.wrapsOrigin ? t('{features} · crosses the origin', { features: featureText }) : featureText,
      selectFeaturesText: tn(featureIds.length, 'Select {count} feature', 'Select {count} features'),
      newDrawingReason: chosen.source.wrapsOrigin ? t('A region that crosses the origin cannot be drawn on its own.') : '',
      // The state of this selection's actions only.
      copy: actionView(chosen, 'copy'),
      add: actionView(chosen, 'add')
    };
  });

  /** Select N features: the feature selection owner selects them by id. */
  const selectFeaturesInside = () => {
    const chosen = openSelection();
    if (chosen) selectFeatures(featuresInside(chosen));
  };
  /**
   * Runs one action of the open selection; another action of it waits until this one ends.
   * @param {RegionAction} name
   * @param {(chosen: RegionSelection) => Promise<unknown>} run
   * @param {string} doneText
   * @param {string} operation The operation a failure names.
   */
  const runAction = async (name, run, doneText, operation) => {
    const chosen = openSelection();
    if (!chosen || (action.selection === chosen && action.busy)) return false;
    Object.assign(action, { selection: chosen, name, busy: true, status: '', failed: false });
    try {
      await run(chosen);
      if (action.selection === chosen) Object.assign(action, { busy: false, status: doneText, failed: false });
      return true;
    } catch (error) {
      if (action.selection === chosen) {
        const failure = normalizeCaughtError(error, { operation, stage: 'helper' });
        Object.assign(action, { busy: false, status: errorText(failure.summaryParts), failed: true });
      }
      return false;
    }
  };
  /** Add to annotation set: the annotation owner adds one "Region N" row (after reading the records, when needed). */
  const addToAnnotationSet = () => runAction('add', (chosen) => addRegionAnnotation(targetOf(chosen)), '', 'addRegionAnnotation');
  /** New drawing from region; refused by Python, and disabled here, for a span across the origin. */
  const newDrawing = async () => {
    const chosen = openSelection();
    if (!chosen || chosen.source.wrapsOrigin) return null;
    const outcome = await newDrawingFromRegion(targetOf(chosen));
    closeRegionPopup();
    return outcome;
  };
  /**
   * Copy FASTA: Python reads the region's bases on the source strand; the clipboard write starts with the click, so it
   * is allowed when the bases arrive.
   */
  const copyFasta = () => runAction(
    'copy', (chosen) => copyPendingTextToClipboard(readRegionFasta(targetOf(chosen))), t('Copied FASTA'), 'readRecordSpanSequences'
  );

  return {
    on, available, active, unavailableReason, toggle, toggleLayoutEdit,
    pointerDown, pointerMove, pointerUp, pointerCancel,
    overlay, overlayMatrix, bindOverlay, bindBand, popup,
    selectFeaturesInside, addToAnnotationSet, newDrawing, copyFasta
  };
};
