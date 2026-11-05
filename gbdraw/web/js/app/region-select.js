// @ts-check
// Select region (PD-OI-115): the preview toolbar toggle, the drag along a
// record's axis, the band over the span, and the actions of the Region popup.
// The span comes only from the axis Python wrote for the shown Result
// (`services/region-axis.js`); the selection lives in the open Region popup
// (the one Result popup, PD-OI-105), so it is no History step, is not saved,
// and is not exported.
/**
 * @import { DisplaySpan, RegionAxis, RegionAxisRecord, RegionSelection, RecordPart } from '../services/region-axis.js'
 */
import {
  axisPosition, circularDragSpan, linearDragSpan, linearRowPosition, partsInsideSpan, regionAxisFrame,
  regionBandPath, regionSpanText, sourceSpanOf
} from '../services/region-axis.js';
import { holdTextSelection } from '../utils/gesture-text-selection.js';
import { isClickAtPress } from '../utils/popup-bounds.js';
import { computed, reactive, ref, shallowRef } from '../utils/reactivity.js';
import { t, tn } from '../utils/i18n.js';

/** How far beside the feature tracks, in screen pixels, a press still starts a drag. */
const BAND_SLACK_PX = 6;

/**
 * A drawn feature of the shown Result as Select region reads it (a catalog row of `extractedFeatures`).
 * @typedef {object} RegionFeatureRow
 * @property {unknown} [svg_id]
 * @property {unknown} [record_key]
 * @property {unknown} [start]
 * @property {unknown} [end]
 * @property {unknown} [location_parts]
 */

/**
 * The region an action works on, in source coordinates (1-based, both ends included).
 * @typedef {object} RegionTarget
 * @property {string} recordKey
 * @property {string} recordId
 * @property {number} recordIndex
 * @property {number} start
 * @property {number} end
 */

/**
 * @typedef {object} RegionSelectPorts
 * @property {() => RegionAxisRecord[]} readAxisRecords The records of the shown Result whose axis Python wrote.
 * @property {() => SVGSVGElement | null} shownRoot The mounted root of the shown Result.
 * @property {{ value: boolean }} layoutEdit Layout edit; this owner is the one that turns it on and off.
 * @property {() => RegionSelection | null} openSelection The selection the open Region popup holds, else null.
 * @property {(selection: RegionSelection, element: Element, anchor: { clientX: number, clientY: number }) => boolean} openRegionPopup
 *   Opens the Region popup (the one Result popup); false when popups cannot open now.
 * @property {() => void} closeRegionPopup
 * @property {() => readonly RegionFeatureRow[]} drawnFeatures The features the shown Result draws.
 * @property {(featureIds: readonly string[]) => void} selectFeatures The feature selection owner's selection by id.
 * @property {(region: RegionTarget) => unknown} addRegionAnnotation The annotation owner's "Region N" row.
 * @property {(region: RegionTarget) => Promise<unknown>} newDrawingFromRegion The region drawing owner's derivation.
 */

/**
 * The parts of a feature in its record's own bases.
 * @param {RegionFeatureRow} feature
 * @returns {RecordPart[]}
 */
const featureParts = (feature) => {
  const listed = Array.isArray(feature.location_parts) && feature.location_parts.length ? feature.location_parts : [feature];
  return listed.flatMap((/** @type {{ start?: unknown, end?: unknown }} */ part) => {
    const start = Number(part?.start);
    const end = Number(part?.end);
    return Number.isSafeInteger(start) && Number.isSafeInteger(end) && end > start ? [{ start, end }] : [];
  });
};

/** @param {RegionSelection} selection @returns {RegionTarget} */
const targetOf = (selection) => ({
  recordKey: selection.recordKey, recordId: selection.recordId, recordIndex: selection.recordIndex,
  start: selection.source.start, end: selection.source.end
});

/** @param {RegionAxis} axis @param {number} start @param {number} current @param {number} sweep */
const dragSpan = (axis, start, current, sweep) => (axis.kind === 'circular'
  ? circularDragSpan(axis, start, sweep)
  : linearDragSpan(start, current));

/**
 * @param {RegionSelectPorts} ports
 */
export const createRegionSelect = ({
  readAxisRecords, shownRoot, layoutEdit, openSelection, openRegionPopup, closeRegionPopup, drawnFeatures,
  selectFeatures, addRegionAnnotation, newDrawingFromRegion
}) => {
  const on = ref(false);
  /** Whether the shown Result has an axis to select on (one predicate for the toggle and the gesture). */
  const available = computed(() => readAxisRecords().length > 0);
  const active = computed(() => on.value && available.value);

  /** Turning Select region on turns Layout edit off. */
  const toggle = () => {
    if (!on.value && !available.value) return;
    on.value = !on.value;
    if (on.value) layoutEdit.value = false;
    else closeRegionPopup();
  };
  /** Turning Layout edit on turns Select region off. */
  const toggleLayoutEdit = () => {
    layoutEdit.value = !layoutEdit.value;
    if (layoutEdit.value && on.value) {
      on.value = false;
      closeRegionPopup();
    }
  };

  /**
   * The drag in progress: its record, where it started along the axis, and the span so far.
   * @type {{ pointerId: number, startX: number, startY: number, moved: boolean, record: RegionAxisRecord,
   *   frame: SVGGraphicsElement, start: number, previous: number, sweep: number, row: number, span: DisplaySpan,
   *   clientX: number, clientY: number } | null}
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

  /** Ends the hold on text selection that a drag takes. */
  let releaseTextSelection = () => {};

  /**
   * A press with Select region on: on a record's feature tracks it starts a drag (and no pan); anywhere else, or
   * with Shift, Ctrl, or Meta held, it keeps its meaning.
   * @param {PointerEvent} event
   */
  const pointerDown = (event) => {
    if (!active.value || event.button !== 0 || event.isPrimary === false) return;
    if (event.shiftKey || event.ctrlKey || event.metaKey || event.altKey) return;
    const hit = recordAt(event.clientX, event.clientY);
    if (!hit) return;
    // A press on the tracks is the drag's: no pan starts and no text is selected.
    event.preventDefault();
    pending = {
      pointerId: event.pointerId, startX: event.clientX, startY: event.clientY, moved: false,
      record: hit.record, frame: hit.frame, start: hit.position.position, previous: hit.position.position, sweep: 0,
      row: hit.position.row, span: { start: Math.floor(hit.position.position), length: 1 },
      clientX: event.clientX, clientY: event.clientY
    };
    releaseTextSelection = holdTextSelection();
  };

  /** @param {PointerEvent} event */
  const pointerMove = (event) => {
    if (!pending || event.pointerId !== pending.pointerId) return;
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
      pending.moved = true;
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
    releaseTextSelection();
    releaseTextSelection = () => {};
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
    if (!pending || event.pointerId !== pending.pointerId) return;
    const finished = pending;
    endDrag();
    if (!finished.moved) return;
    event.preventDefault();
    event.stopPropagation();
    const { record } = finished;
    /** @type {RegionSelection} */
    const selection = {
      kind: 'region', resultIndex: record.resultIndex, recordIndex: record.recordIndex, recordKey: record.recordKey,
      recordId: record.recordId, span: finished.span, source: sourceSpanOf(record.axis.source, finished.span)
    };
    shownSelection.value = selection;
    const anchor = { clientX: finished.clientX, clientY: finished.clientY };
    window.requestAnimationFrame(() => {
      if (shownSelection.value !== selection) return;
      // Once open, the popup holds the selection: it goes when the popup closes.
      if (bandElement) openRegionPopup(selection, bandElement, anchor);
      shownSelection.value = null;
    });
  };

  const pointerCancel = () => { if (pending) endDrag(); };

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

  /** The Region popup's body: the feature count and the actions. */
  const popup = computed(() => {
    const chosen = openSelection();
    if (!chosen) return null;
    const featureIds = featuresInside(chosen);
    return {
      featureCount: featureIds.length,
      featureText: tn(featureIds.length, '{count} feature', '{count} features'),
      selectFeaturesText: tn(featureIds.length, 'Select {count} feature', 'Select {count} features'),
      newDrawingReason: chosen.source.wrapsOrigin ? t('A region that crosses the origin cannot be drawn on its own.') : ''
    };
  });

  /** Select N features: the feature selection owner selects them by id. */
  const selectFeaturesInside = () => {
    const chosen = openSelection();
    if (chosen) selectFeatures(featuresInside(chosen));
  };
  /** Add to annotation set: the annotation owner adds one "Region N" row. */
  const addToAnnotationSet = () => {
    const chosen = openSelection();
    return chosen ? addRegionAnnotation(targetOf(chosen)) : null;
  };
  /** New drawing from region; refused by Python, and disabled here, for a span across the origin. */
  const newDrawing = async () => {
    const chosen = openSelection();
    if (!chosen || chosen.source.wrapsOrigin) return null;
    const outcome = await newDrawingFromRegion(targetOf(chosen));
    closeRegionPopup();
    return outcome;
  };

  return {
    on, available, active, toggle, toggleLayoutEdit,
    pointerDown, pointerMove, pointerUp, pointerCancel,
    overlay, overlayMatrix, bindOverlay, bindBand, popup,
    selectFeaturesInside, addToAnnotationSet, newDrawing
  };
};
