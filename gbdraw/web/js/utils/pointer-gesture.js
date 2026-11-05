// @ts-check
// A drag on the preview (the Shift+drag feature marquee, a Select region drag)
// owns the pointer: while it runs, the browser starts no text selection over
// the figure's labels (OV-2478), and the element the drag belongs to captures
// the pointer, so the move and the release reach the drag also over the
// toolbar or outside the window.

/** @param {Event} event */
const refuseSelection = (event) => event.preventDefault();

/**
 * Holds off native text selection until the returned release is called.
 * @returns {() => void}
 */
export const holdTextSelection = () => {
  if (typeof document === 'undefined') return () => {};
  document.addEventListener('selectstart', refuseSelection, true);
  return () => document.removeEventListener('selectstart', refuseSelection, true);
};

/**
 * Captures the pointer on `element` until the returned release is called. Capture is best effort (an SVG element
 * may refuse it); the release then does nothing.
 * @param {Element} element
 * @param {number} pointerId
 * @returns {() => void}
 */
export const capturePointer = (element, pointerId) => {
  if (typeof element.setPointerCapture !== 'function' || !Number.isFinite(pointerId)) return () => {};
  try {
    element.setPointerCapture(pointerId);
  } catch (_err) {
    return () => {};
  }
  return () => {
    try {
      element.releasePointerCapture(pointerId);
    } catch (_err) {
      // The browser released it already (the pointer went up or was cancelled).
    }
  };
};
