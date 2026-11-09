// @ts-check
/** @import { PersistedPreviewViewport, StoredFigurePreview } from '../services/config.js' */
/** @import { PressPoint } from '../utils/popup-bounds.js' */
/** @import { WebState } from '../state.js' */
import { recordStructuralMetric } from '../services/runtime-test-hooks.js';
import { isClickAtPress, primaryPressAfter } from '../utils/popup-bounds.js';
import { watch } from '../utils/reactivity.js';

// Keep one burst alive through the existing 0.2 s transform transition. This
// also absorbs browser delivery delay between nominally 80 ms wheel inputs.
const WHEEL_BURST_QUIET_MS = 220;
const WHEEL_TRANSITION_FALLBACK_MS = 260;
// The preview zoom range; the zoom-in and zoom-out buttons repeat it in index.html.
const MIN_ZOOM = 0.1;
const MAX_ZOOM = 5;
// Fit leaves the canvas padding (p-2) free around the diagram.
const FIT_MARGIN_PX = 8;

/** @param {number} value The zoom range and the wheel's 0.1 steps. */
const clampZoom = (value) => Math.round(Math.max(MIN_ZOOM, Math.min(MAX_ZOOM, value)) * 10) / 10;

/**
 * The width the open Editor drawer covers at the right of the canvas: the
 * stylesheet's `--preview-editor-reserve`, which the search and controls keep
 * free too (0 when the drawer is closed or docked below the canvas).
 * @param {HTMLElement} container The preview canvas.
 */
const editorReserveWidth = (container) => {
  const probe = container.ownerDocument?.createElement?.('div');
  if (!probe) return 0;
  probe.style.cssText = 'position:absolute;visibility:hidden;height:0;width:var(--preview-editor-reserve,0px)';
  container.appendChild(probe);
  const { width } = probe.getBoundingClientRect();
  probe.remove();
  return width;
};

// The room of the figure is the canvas box, which a zoomed figure's scrollbars do not change.
/** @param {HTMLElement} container The preview canvas. */
const previewRoom = (container) => `${container.offsetWidth}x${container.offsetHeight}`;

// A press on these starts no pan: a form control, or a figure element that its own editing owns.
/** @param {Element | null | undefined} target */
const isFormControlTarget = (target) =>
  Boolean(target?.closest?.('button, input, textarea, select, a, [role="button"]'));

/** @param {Element | null | undefined} target */
const isSvgEditingTarget = (target) =>
  Boolean(
    target?.closest?.(
      [
        'text[data-label-editable="true"]',
        '[data-gbdraw-feature-id]',
        'path[id^="f"]',
        'polygon[id^="f"]',
        'rect[id^="f"]',
        '[data-gbdraw-match-id]',
        '[data-match-kind]',
        '[data-pairwise-match-style]',
        '[data-collinearity-block-id]',
        '[data-collinear-group-scope]'
      ].join(', ')
    )
  );

/**
 * What a preview transform interaction (wheel zoom or pan) tells its subscribers when it starts or ends.
 * @typedef {object} PreviewTransformInteractionChange
 * @property {boolean} active
 * @property {string} kind The interaction: `wheel`, `pan`, or `cancel`.
 * @property {Event | null} event The event that started it; null when it ends.
 * @property {boolean} reconcile Whether the subscriber reconciles what the interaction deferred.
 */

/** @typedef {(change: PreviewTransformInteractionChange) => void} PreviewTransformInteractionListener */

/**
 * The figure a file Session Load mounts stays hidden until the Load ends, and is then shown, fitted at once
 * to the Preview the Load's notice and phone pane switch sized, or as it is (OV-2702). The stored figure a
 * Session 27-44 Load shows while Python converts the file is hidden and shown fitted in the same way, and a
 * Load that withdraws it shows the view from before it (OV-2703). A figure to show fitted while its canvas
 * has no room (a phone shows the Settings pane) stays hidden until the canvas gets one, and is then fitted
 * and shown; hiding a figure again, or showing it without the Fit, drops that wait.
 * @param {{ hidden: { value: boolean }, canvas: { value: HTMLElement | null }, fitAtOnce: () => void,
 *   setPreviewViewport: (pan: { x: number, y: number }, zoomLevel: number) => void }} owner The pan and zoom
 *   owner's hold, canvas, Fit shown at once, and whole-view setter.
 */
const createHiddenPreview = ({ hidden, canvas, fitAtOnce, setPreviewViewport }) => {
  /** @param {HTMLElement | null} container */
  const hasRoom = (container) => Boolean(container?.offsetWidth && container.offsetHeight);
  const roomObserver = typeof ResizeObserver === 'function'
    ? new ResizeObserver(() => { if (hasRoom(canvas.value)) showHiddenPreview({ fit: true }); }) : null;
  const hidePreviewUntilShown = () => {
    roomObserver?.disconnect();
    hidden.value = true;
  };
  /** @param {{ fit: boolean }} options `fit` shows the figure fitted, with no transition to the fitted view. */
  const showHiddenPreview = ({ fit }) => {
    roomObserver?.disconnect();
    if (fit && roomObserver && canvas.value && !hasRoom(canvas.value)) {
      roomObserver.observe(canvas.value);
      return;
    }
    if (fit) fitAtOnce();
    hidden.value = false;
  };
  /** @type {StoredFigurePreview} */
  const storedFigurePreview = Object.freeze({
    hide: hidePreviewUntilShown, show: showHiddenPreview,
    restore: (/** @type {PersistedPreviewViewport} */ view) => setPreviewViewport(view.canvasPan, view.zoom)
  });
  return { hidePreviewUntilShown, showHiddenPreview, storedFigurePreview, disposeHiddenPreview: () => roomObserver?.disconnect() };
};

/**
 * The view the Preview surface shows.
 * @typedef {Pick<WebState, 'canvasPan' | 'zoom' | 'isPanning' | 'previewTransitionHeld' | 'previewHiddenUntilShown'>}
 *   PreviewSurfaceView
 */

/**
 * The Preview surface's style, projected from the view: the one writer of the surface's transform, transform
 * origin, transition, will-change, and visibility, from `v-preview-transform`'s creation of the surface to its
 * unmount. Outside a gesture or a held transition the transform moves through a 0.2 s transition; during one
 * it follows at once, as a promoted layer. A change of the view shows at Vue's next flush, before the frame is
 * drawn, once however many members changed (a pan frame changes both coordinates); a held transition shows at
 * once, so the Fit measures, and the Fit shown at once reads, the style its write set. No app template reads
 * the view (OV-912).
 * @param {HTMLElement} element The Preview surface.
 * @param {PreviewSurfaceView} view
 * @returns {() => void} Stops the projection.
 */
export const projectPreviewSurface = (element, view) => {
  const show = () => {
    const immediate = view.isPanning.value || view.previewTransitionHeld.value;
    element.style.transform = `translate(${view.canvasPan.x}px, ${view.canvasPan.y}px) scale(${view.zoom.value})`;
    element.style.transformOrigin = 'top center';
    element.style.transition = immediate ? 'none' : 'transform 0.2s';
    element.style.willChange = immediate ? 'transform' : '';
    element.style.visibility = view.previewHiddenUntilShown.value ? 'hidden' : '';
  };
  const stops = [
    watch(() => [view.canvasPan.x, view.canvasPan.y, view.zoom.value, view.isPanning.value,
      view.previewHiddenUntilShown.value], show, { immediate: true }),
    watch(() => view.previewTransitionHeld.value, show, { flush: 'sync' })
  ];
  return () => stops.forEach((stop) => stop());
};

/**
 * What the preview transform interaction shows its readers: whether a wheel zoom or a pan is active, and its
 * start and end.
 * @typedef {Readonly<{ isActive: () => boolean,
 *   subscribe: (listener: PreviewTransformInteractionListener) => () => void }>} PreviewTransformInteraction
 */

/**
 * The preview transform interaction: the wheel zoom and the pan that hold it, as sources, and its readers. A
 * wheel zoom ends once its burst is quiet and the surface's transform transition has ended, or past that
 * transition when a clamped zoom runs none. It writes no view state; a pan the caller stops itself.
 * @param {WebState['svgContainer']} svgContainer The Preview surface, whose transform transition ends a wheel zoom.
 * @returns {{ previewTransformInteraction: PreviewTransformInteraction,
 *   beginPreviewTransformInteraction: (kind: string, event: Event) => boolean,
 *   endPreviewTransformInteraction: (kind: string) => boolean, beginWheelInteraction: (event: WheelEvent) => void,
 *   endEveryInteraction: (reconcile: boolean) => boolean, clearInteractionListeners: () => void }}
 */
const createPreviewTransformInteraction = (svgContainer) => {
  /** @type {Set<string>} */
  const previewTransformInteractionSources = new Set();
  /** @type {number | null} */
  let wheelBurstTimerId = null;
  /** @type {number | null} */
  let wheelFallbackTimerId = null;
  /** @type {HTMLElement | null} */
  let wheelTransitionTarget = null;
  let wheelBurstComplete = false;
  let wheelTransitionComplete = false;
  /** @type {Set<PreviewTransformInteractionListener>} */
  const previewTransformInteractionListeners = new Set();

  /** @type {PreviewTransformInteraction} */
  const previewTransformInteraction = Object.freeze({
    isActive: () => previewTransformInteractionSources.size > 0,
    /** @param {PreviewTransformInteractionListener} listener */
    subscribe(listener) {
      if (typeof listener !== 'function') return () => {};
      previewTransformInteractionListeners.add(listener);
      return () => previewTransformInteractionListeners.delete(listener);
    }
  });

  /**
   * @param {boolean} active
   * @param {string} kind
   * @param {Event | null} event
   * @param {boolean} reconcile
   */
  const notifyPreviewTransformInteraction = (active, kind, event, reconcile) => {
    previewTransformInteractionListeners.forEach((listener) => {
      listener({ active, kind, event, reconcile });
    });
  };

  /**
   * @param {string} kind
   * @param {Event} event
   */
  const beginPreviewTransformInteraction = (kind, event) => {
    const wasActive = previewTransformInteractionSources.size > 0;
    previewTransformInteractionSources.add(kind);
    if (wasActive) return false;
    recordStructuralMetric('previewTransformInteractionStartCount', 1, { kind });
    notifyPreviewTransformInteraction(true, kind, event, false);
    return true;
  };

  /** @param {string} kind The source that ends; the interaction ends, reconciled, with the last one. */
  const endPreviewTransformInteraction = (kind) => {
    if (!previewTransformInteractionSources.delete(kind)) return false;
    if (previewTransformInteractionSources.size > 0) return false;
    recordStructuralMetric('previewTransformInteractionEndCount', 1, { kind });
    notifyPreviewTransformInteraction(false, kind, null, true);
    return true;
  };

  const clearWheelTimers = () => {
    if (wheelBurstTimerId !== null) {
      window.clearTimeout(wheelBurstTimerId);
      wheelBurstTimerId = null;
    }
    if (wheelFallbackTimerId !== null) {
      window.clearTimeout(wheelFallbackTimerId);
      wheelFallbackTimerId = null;
    }
  };

  // A wheel zoom's timers, transition listener, and end flags, dropped when it ends or is cancelled.
  const dropWheelEnd = () => {
    clearWheelTimers();
    wheelTransitionTarget?.removeEventListener?.('transitionend', handleWheelTransitionEnd);
    wheelTransitionTarget = null;
    wheelBurstComplete = false;
    wheelTransitionComplete = false;
  };

  const finishWheelInteraction = () => {
    if (!previewTransformInteractionSources.has('wheel')) return false;
    dropWheelEnd();
    return endPreviewTransformInteraction('wheel');
  };

  /** @param {TransitionEvent} event */
  function handleWheelTransitionEnd(event) {
    if (
      !previewTransformInteractionSources.has('wheel')
      || event?.target !== wheelTransitionTarget
      || event?.propertyName !== 'transform'
    ) return;
    wheelTransitionComplete = true;
    if (wheelBurstComplete) finishWheelInteraction();
  }

  const scheduleWheelInteractionEnd = () => {
    clearWheelTimers();
    wheelBurstComplete = false;
    wheelTransitionComplete = false;
    wheelBurstTimerId = window.setTimeout(() => {
      wheelBurstTimerId = null;
      wheelBurstComplete = true;
      if (wheelTransitionComplete) finishWheelInteraction();
    }, WHEEL_BURST_QUIET_MS);
    // A clamped zoom may not produce transitionend. Bound that path just past
    // the existing 0.2 s transform transition without changing the transition.
    wheelFallbackTimerId = window.setTimeout(() => {
      wheelFallbackTimerId = null;
      finishWheelInteraction();
    }, WHEEL_TRANSITION_FALLBACK_MS);
  };

  /** @param {WheelEvent} event */
  const beginWheelInteraction = (event) => {
    const isNewWheelInteraction = !previewTransformInteractionSources.has('wheel');
    beginPreviewTransformInteraction('wheel', event);
    if (isNewWheelInteraction) {
      wheelTransitionTarget = svgContainer.value;
      wheelTransitionTarget?.addEventListener?.('transitionend', handleWheelTransitionEnd);
    }
    scheduleWheelInteractionEnd();
  };

  /** Ends every source at once, as one `cancel` end. @param {boolean} reconcile */
  const endEveryInteraction = (reconcile) => {
    const wasActive = previewTransformInteractionSources.size > 0;
    dropWheelEnd();
    previewTransformInteractionSources.clear();
    if (!wasActive) return false;
    recordStructuralMetric('previewTransformInteractionEndCount', 1, { kind: 'cancel' });
    notifyPreviewTransformInteraction(false, 'cancel', null, reconcile);
    return true;
  };

  return {
    previewTransformInteraction, beginPreviewTransformInteraction, endPreviewTransformInteraction,
    beginWheelInteraction, endEveryInteraction, clearInteractionListeners: () => previewTransformInteractionListeners.clear()
  };
};

/**
 * The pan gesture: a press on the canvas that no control or figure editing owns starts it, the pointer moves
 * the view once a frame, and its release or cancel ends it. It writes the pan and the canvas cursor.
 * @param {Pick<WebState, 'layoutRepositionMode' | 'isPanning' | 'panStart' | 'canvasPan' | 'canvasContainerRef'
 *   | 'svgContainer'>} state
 * @param {{ beginInteraction: (kind: string, event: Event) => boolean, endInteraction: (kind: string) => boolean }}
 *   interaction The transform interaction's start and end of the `pan` source.
 * @returns {{ startPan: (event: MouseEvent & { target: Element }) => void, doPan: (event: MouseEvent) => void,
 *   endPan: (event: MouseEvent) => void, stopPan: () => void }}
 */
const createPanGesture = (state, { beginInteraction, endInteraction }) => {
  const { layoutRepositionMode, isPanning, panStart, canvasPan, canvasContainerRef, svgContainer } = state;
  /** @type {number | null} */
  let panFrameId = null;
  /** @type {{ x: number, y: number } | null} */
  let pendingPanPointer = null;

  const cancelPanFrame = () => {
    if (panFrameId !== null) {
      cancelAnimationFrame(panFrameId);
      panFrameId = null;
    }
  };

  /** @param {number} clientX @param {number} clientY */
  const getPanPosition = (clientX, clientY) => {
    const dx = clientX - panStart.x;
    const dy = clientY - panStart.y;
    return {
      x: panStart.panX + dx,
      y: panStart.panY + dy
    };
  };

  // A pan writes its live position to the view, so every reader (the surface, a wheel zoom during the
  // pan) sees where the figure is (OV-2722).
  /** @param {number} clientX @param {number} clientY */
  const flushPanUpdate = (clientX, clientY) => {
    const nextPan = getPanPosition(clientX, clientY);
    canvasPan.x = nextPan.x;
    canvasPan.y = nextPan.y;
  };

  // The pan stops: no frame, no pointer, and the canvas cursor shows it can pan again.
  const stopPan = () => {
    cancelPanFrame();
    pendingPanPointer = null;
    isPanning.value = false;
    if (canvasContainerRef.value) canvasContainerRef.value.style.cursor = 'grab';
  };

  /** @param {MouseEvent & { target: Element }} event A press on the canvas; its target is an element. */
  const startPan = (event) => {
    if (event.button !== 0) return;
    if (event.shiftKey && svgContainer.value?.querySelector?.('svg')) return;
    const container = canvasContainerRef.value;
    if (!container) return;

    const target = event.target;
    if (isFormControlTarget(target) || isSvgEditingTarget(target)) return;

    const closestGroup = target.closest?.('g[id]');
    if (closestGroup) {
      const groupId = closestGroup.id;
      if (groupId.startsWith('f') && !target.closest('.gbdraw-preview-layout-target')) {
        return;
      }
      if (layoutRepositionMode.value && target.closest('svg')) return;
    }
    if (layoutRepositionMode.value && target.tagName === 'path' && target.closest('svg')) {
      return;
    }

    // Own the accepted pan instead of starting native SVG text selection/drag.
    event.preventDefault?.();
    cancelPanFrame();
    pendingPanPointer = null;
    beginInteraction('pan', event);
    isPanning.value = true;
    panStart.x = event.clientX;
    panStart.y = event.clientY;
    panStart.panX = canvasPan.x;
    panStart.panY = canvasPan.y;
    container.style.cursor = 'grabbing';
  };

  /** @param {MouseEvent} event */
  const doPan = (event) => {
    if (!isPanning.value) return;
    pendingPanPointer = { x: event.clientX, y: event.clientY };
    if (panFrameId !== null) return;
    panFrameId = requestAnimationFrame(() => {
      panFrameId = null;
      if (!isPanning.value || !pendingPanPointer) return;
      flushPanUpdate(pendingPanPointer.x, pendingPanPointer.y);
    });
  };

  /** @param {MouseEvent} event A mouse event or a pointer cancel / capture loss. */
  const endPan = (event) => {
    const wasPanning = isPanning.value;
    // Cancellation and capture loss may carry (0, 0), not a pointer position.
    const finalPointer =
      event?.type !== 'pointercancel' && event?.type !== 'lostpointercapture'
      && Number.isFinite(event?.clientX) && Number.isFinite(event?.clientY)
        ? { x: event.clientX, y: event.clientY }
        : pendingPanPointer;
    if (wasPanning && finalPointer) flushPanUpdate(finalPointer.x, finalPointer.y);
    stopPan();
    if (wasPanning) endInteraction('pan');
  };

  return { startPan, doPan, endPan, stopPan };
};

/**
 * Fit (UI-11): the largest whole-percent zoom that shows the whole diagram in the visible preview frame (the
 * canvas left of an open Editor drawer), and the pan that centres it under the top-center origin. It keeps
 * the room it measured, which the held view reads.
 * @param {Pick<WebState, 'zoom' | 'canvasPan' | 'canvasContainerRef' | 'svgContainer' | 'previewTransitionHeld'>
 *   & { setPreviewViewport: (pan: { x: number, y: number }, zoomLevel: number) => void }} owner The view, the
 *   canvas and surface, and the whole-view setter.
 * @returns {{ fitPreviewToViewport: () => void, inFittedRoom: (container: HTMLElement) => boolean,
 *   forgetFittedRoom: () => void }}
 */
const createPreviewFit = ({ zoom, canvasPan, canvasContainerRef, svgContainer, previewTransitionHeld, setPreviewViewport }) => {
  // The room the last Fit measured.
  let fittedRoom = '';

  const fitPreviewToViewport = () => {
    const container = canvasContainerRef.value;
    const surface = svgContainer.value;
    const svg = surface?.querySelector?.('svg');
    if (!container || !surface || !svg) return;
    // Measure the current transform itself, not a frame of its transition.
    previewTransitionHeld.value = true;
    const surfaceBox = surface.getBoundingClientRect();
    const svgBox = svg.getBoundingClientRect();
    const frameBox = container.getBoundingClientRect();
    const frameWidth = container.clientWidth - editorReserveWidth(container);
    const width = svgBox.width / zoom.value;
    const height = svgBox.height / zoom.value;
    if (!(width > 0 && height > 0)) {
      previewTransitionHeld.value = false;
      return;
    }
    fittedRoom = previewRoom(container);
    const largest = Math.min(
      (frameWidth - 2 * FIT_MARGIN_PX) / width,
      (container.clientHeight - 2 * FIT_MARGIN_PX) / height
    );
    // Round down to a whole percent so the diagram fills the frame without
    // overflowing it; the wheel and buttons step by 0.1 from there.
    const nextZoom = Math.max(MIN_ZOOM, Math.min(MAX_ZOOM, Math.floor(largest * 100 + 1e-9) / 100));
    // The transform origin without the current pan, and the diagram centre's
    // unscaled offset from it.
    const originX = surfaceBox.left + surfaceBox.width / 2 - canvasPan.x;
    const originY = surfaceBox.top - canvasPan.y;
    const centerX = (svgBox.left + svgBox.width / 2 - surfaceBox.left - surfaceBox.width / 2) / zoom.value;
    const centerY = (svgBox.top + svgBox.height / 2 - surfaceBox.top) / zoom.value;
    setPreviewViewport({
      x: frameBox.left + container.clientLeft + frameWidth / 2 - originX - nextZoom * centerX,
      y: frameBox.top + container.clientTop + container.clientHeight / 2 - originY - nextZoom * centerY
    }, nextZoom);
    previewTransitionHeld.value = false;
  };

  return {
    fitPreviewToViewport,
    inFittedRoom: (container) => previewRoom(container) === fittedRoom,
    forgetFittedRoom: () => { fittedRoom = ''; }
  };
};

/**
 * A view shown only for a while (the fitted figure of Review alignment) is not what a Session or an artifact
 * checkpoint keeps. The view from before it is held here, and the Session writer and the UI capture read it
 * through `heldPreviewViewport`. While a view is held and the figure fitted in its place is shown, that figure
 * follows its room: it is fitted again each time the canvas changes size, whatever changed it (the window, the
 * settings pane, the docked panel changing sides). A zoom or a pan made in between lasts until the next such
 * change. A hidden canvas (the Run info tab) has no room, and the room of the last Fit is no change. Another
 * figure shown under the hold is not fitted.
 * @param {Pick<WebState, 'zoom' | 'canvasPan' | 'canvasContainerRef'> & { fit: () => void,
 *   inFittedRoom: (container: HTMLElement) => boolean, forgetFittedRoom: () => void,
 *   setPreviewViewport: (pan: { x: number, y: number }, zoomLevel: number) => void }} owner The view, the
 *   canvas, the Fit and the room it measured, and the whole-view setter.
 * @returns {{ holdPreviewViewport: (figureShown: () => boolean) => void,
 *   releasePreviewViewport: (options: { restore: boolean }) => void,
 *   heldPreviewViewport: () => PersistedPreviewViewport | null, disposeHeldPreviewView: () => void }}
 */
const createHeldPreviewView = ({ zoom, canvasPan, canvasContainerRef, fit, inFittedRoom, forgetFittedRoom, setPreviewViewport }) => {
  /** @type {PersistedPreviewViewport | null} */
  let held = null;
  /** @type {() => boolean} */
  let heldFigureShown = () => false;
  const heldRoomObserver = typeof ResizeObserver === 'function' ? new ResizeObserver(() => {
    const container = canvasContainerRef.value;
    if (!container?.offsetWidth || !container.offsetHeight || inFittedRoom(container)
      || !heldFigureShown()) return;
    fit();
  }) : null;
  /** @param {() => boolean} figureShown Whether the figure fitted in place of the held view is shown. */
  const holdPreviewViewport = (figureShown) => {
    if (held) return;
    held = { zoom: zoom.value, canvasPan: { x: canvasPan.x, y: canvasPan.y } };
    heldFigureShown = figureShown;
    forgetFittedRoom();
    if (canvasContainerRef.value) heldRoomObserver?.observe(canvasContainerRef.value);
  };
  /** @param {{ restore: boolean }} options `restore` shows the held view again; without it the hold is dropped. */
  const releasePreviewViewport = ({ restore }) => {
    heldRoomObserver?.disconnect();
    const view = held;
    held = null;
    if (view && restore) setPreviewViewport(view.canvasPan, view.zoom);
  };
  return {
    holdPreviewViewport, releasePreviewViewport, heldPreviewViewport: () => held,
    disposeHeldPreviewView: () => heldRoomObserver?.disconnect()
  };
};

/**
 * The pan and zoom owner. It writes the view state only; `projectPreviewSurface` shows it. Its steps (the
 * transform interaction, the pan gesture, the Fit, and the held view) each own their state; this owner joins
 * them and owns the whole-view setter.
 * @param {Pick<WebState, 'zoom' | 'layoutRepositionMode' | 'isPanning' | 'panStart' | 'canvasPan' | 'canvasContainerRef'
 *   | 'svgContainer' | 'previewHiddenUntilShown' | 'previewTransitionHeld'>} state
 */
export const createPanZoom = (state) => {
  const { zoom, panStart, canvasPan, canvasContainerRef, svgContainer, previewHiddenUntilShown, previewTransitionHeld } = state;
  const {
    previewTransformInteraction, beginPreviewTransformInteraction, endPreviewTransformInteraction, beginWheelInteraction,
    endEveryInteraction, clearInteractionListeners
  } = createPreviewTransformInteraction(svgContainer);
  const { startPan, doPan, endPan, stopPan } = createPanGesture(state, {
    beginInteraction: beginPreviewTransformInteraction, endInteraction: endPreviewTransformInteraction
  });

  // A cancel stops a pan, then ends every source of the interaction.
  const cancelPreviewTransformInteraction = ({ reconcile = false } = {}) => {
    stopPan();
    return endEveryInteraction(reconcile);
  };

  /**
   * Replace the whole preview viewport; Reset and Fit both end here.
   * @param {{ x: number, y: number }} pan
   * @param {number} zoomLevel
   */
  const setPreviewViewport = (pan, zoomLevel) => {
    cancelPreviewTransformInteraction({ reconcile: false });
    Object.assign(panStart, { x: 0, y: 0, panX: 0, panY: 0 });
    canvasPan.x = pan.x;
    canvasPan.y = pan.y;
    zoom.value = zoomLevel;
  };

  /** @param {{ resetZoom?: boolean, pan?: { x?: number, y?: number } | null }} [options] */
  const resetPreviewViewport = ({ resetZoom = false, pan = null } = {}) => setPreviewViewport(
    { x: Number(pan?.x) || 0, y: Number(pan?.y) || 0 },
    resetZoom ? 1.0 : zoom.value
  );

  /** @param {WheelEvent} event */
  const handleWheel = (event) => {
    beginWheelInteraction(event);
    zoom.value = clampZoom(zoom.value + (event.deltaY > 0 ? -0.1 : 0.1));
  };

  const { fitPreviewToViewport, inFittedRoom, forgetFittedRoom } = createPreviewFit({
    zoom, canvasPan, canvasContainerRef, svgContainer, previewTransitionHeld, setPreviewViewport
  });
  const { holdPreviewViewport, releasePreviewViewport, heldPreviewViewport, disposeHeldPreviewView } = createHeldPreviewView({
    zoom, canvasPan, canvasContainerRef, fit: fitPreviewToViewport, inFittedRoom, forgetFittedRoom, setPreviewViewport
  });

  // The Fit shown at once: no transition runs to the fitted view, and the browser takes that view
  // as shown before the transition comes back.
  const fitAtOnce = () => {
    if (!svgContainer.value) return;
    fitPreviewToViewport();
    previewTransitionHeld.value = true;
    void getComputedStyle(svgContainer.value).transform;
    previewTransitionHeld.value = false;
  };
  const { hidePreviewUntilShown, showHiddenPreview, storedFigurePreview, disposeHiddenPreview } = createHiddenPreview({
    hidden: previewHiddenUntilShown, canvas: canvasContainerRef, fitAtOnce, setPreviewViewport
  });

  const disposePanZoom = () => {
    disposeHeldPreviewView();
    disposeHiddenPreview();
    cancelPreviewTransformInteraction({ reconcile: false });
    clearInteractionListeners();
  };

  return {
    handleWheel,
    startPan,
    doPan,
    endPan,
    resetPreviewViewport,
    fitPreviewToViewport,
    holdPreviewViewport,
    releasePreviewViewport,
    heldPreviewViewport,
    hidePreviewUntilShown,
    showHiddenPreview,
    storedFigurePreview,
    cancelPreviewTransformInteraction,
    disposePanZoom,
    previewTransformInteraction
  };
};

/** @param {Pick<WebState, 'sidebarWidth' | 'isResizing'>} state */
export const createSidebarResize = (state) => {
  const { sidebarWidth, isResizing } = state;
  // Where the pointer pressed the handle, as its distance from the pane's width: the pane
  // follows the pointer from the first move (OV-620).
  let grabOffset = 0;

  /** @param {MouseEvent} event */
  const doResize = (event) => {
    if (!isResizing.value) return;
    sidebarWidth.value = Math.max(240, Math.min(500, event.clientX - grabOffset));
  };

  const stopResizing = () => {
    isResizing.value = false;
    document.removeEventListener('mousemove', doResize);
    document.removeEventListener('mouseup', stopResizing);
  };

  /** @param {MouseEvent} event The press on the handle. */
  const startResizing = (event) => {
    // The press is default-prevented, so the field is blurred here.
    if (document.activeElement instanceof HTMLElement) document.activeElement.blur();
    grabOffset = event.clientX - sidebarWidth.value;
    isResizing.value = true;
    document.addEventListener('mousemove', doResize);
    document.addEventListener('mouseup', stopResizing);
  };

  return { startResizing };
};

/**
 * @typedef {Object} GlobalUiEventsOptions
 * @property {Pick<WebState, 'clickedFeature' | 'clickedPairwiseMatch' | 'showCanvasControls'>} state
 *   The popups that Escape and a click outside close.
 * @property {(callback: () => void) => void} onMounted
 * @property {(callback: () => void) => void} onUnmounted
 * @property {() => void} closeRightDrawer
 * @property {() => boolean} previewPopupsHeld While it holds, a click outside leaves the feature and match popups as they are.
 * @property {() => boolean} closeSettingsPopup The pane's closing of the popup that holds a settings section; whether it was open.
 */

/** @param {GlobalUiEventsOptions} options */
export const setupGlobalUiEvents = ({
  state,
  onMounted,
  onUnmounted,
  closeRightDrawer,
  previewPopupsHeld,
  closeSettingsPopup
}) => {
  const {
    clickedFeature,
    clickedPairwiseMatch,
    showCanvasControls
  } = state;

  // UI-02: a modal dialog answers Escape and its own clicks (its Cancel), so
  // nothing behind it closes. The target is read through closest() because
  // a choice unmounts the dialog before this document listener runs.
  const modalDialogOpen = () => Boolean(document.querySelector('[role="dialog"][aria-modal="true"]'));

  // A click whose press began in a popup (a drag released past the popup's
  // clamp, a text selection) is not a click outside it.
  let pressStartedInPopup = false;
  /** @type {PressPoint | null} */
  let pressPoint = null;
  // Only the press of the primary button, which clicks, counts (`primaryPressAfter`): from Pointer
  // Events, so a press whose pointerdown a tool cancelled (no mousedown follows) counts too (OV-2705).
  /** @param {PointerEvent} e A pointerdown or a pointermove. */
  const rememberPressTarget = (e) => {
    const press = primaryPressAfter(pressPoint, e);
    if (press === pressPoint) return;
    pressPoint = press;
    pressStartedInPopup = e.target instanceof Element
      && Boolean(e.target.closest('.feature-popup, .pairwise-match-popup'));
  };

  /** @param {MouseEvent & { target: Element }} e A click; the popups and dialogs are elements. */
  const closeFeaturePopup = (e) => {
    const startedInPopup = pressStartedInPopup;
    pressStartedInPopup = false;
    if (startedInPopup) return;
    if (
      !e.target.closest('[data-modal-overlay], [aria-modal="true"]')
      && !e.target.closest('.feature-popup')
      && !e.target.closest('.pairwise-match-popup')
      && !e.target.closest('[data-similarity-alignment-overlay]')
      // The notice banner's Dismiss is not a click on the diagram.
      && !e.target.closest('[data-app-notice]')
    ) {
      if (previewPopupsHeld()) return;
      if (clickedFeature.value) clickedFeature.value = null;
      if (clickedPairwiseMatch.value) clickedPairwiseMatch.value = null;
      // The popup of a settings section stays through clicks elsewhere (Generate, Undo, the
      // pane) and through a pan: a click on the Result's canvas that nothing took closes it.
      if (e.target.closest('.preview-canvas') && isClickAtPress(pressPoint, e)) closeSettingsPopup();
    }
  };

  // Escape closes the top layer only: an open popup first, so the drawer that
  // opened it stays and focus returns to its opener there; then the canvas
  // controls and the drawer.
  /** @param {KeyboardEvent} e */
  const handleEscapeKey = (e) => {
    if (e.key !== 'Escape' || modalDialogOpen()) return;
    if (clickedFeature.value || clickedPairwiseMatch.value) {
      if (clickedFeature.value) clickedFeature.value = null;
      if (clickedPairwiseMatch.value) clickedPairwiseMatch.value = null;
      return;
    }
    if (closeSettingsPopup()) return;
    if (showCanvasControls.value) showCanvasControls.value = false;
    closeRightDrawer();
  };

  // A click that reaches the document comes from an element, as `closeFeaturePopup` reads it:
  // the cast to EventListener widens its parameter back to Event.
  onMounted(() => {
    document.addEventListener('pointerdown', rememberPressTarget, true);
    document.addEventListener('pointermove', rememberPressTarget, true);
    document.addEventListener('click', /** @type {EventListener} */ (closeFeaturePopup));
    document.addEventListener('keydown', handleEscapeKey);
  });

  onUnmounted(() => {
    document.removeEventListener('pointerdown', rememberPressTarget, true);
    document.removeEventListener('pointermove', rememberPressTarget, true);
    document.removeEventListener('click', /** @type {EventListener} */ (closeFeaturePopup));
    document.removeEventListener('keydown', handleEscapeKey);
  });
};
