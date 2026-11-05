// @ts-check
// A drag on the preview (the Shift+drag feature marquee, a Select region drag)
// owns the pointer: while it runs, the browser starts no text selection over
// the figure's labels (OV-2478).

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
