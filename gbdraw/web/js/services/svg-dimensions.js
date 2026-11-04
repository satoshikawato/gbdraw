// @ts-check
// The size of an SVG root in px, the one reader for export (services/export.js) and the Figure size line
// (app/figure-size-field.js). A root `width`/`height` in px (a bare number or `px`) reads as it is. The page's own
// Results are px (gbdraw/web_support/request_render.py); a Result a CLI `--save_session` Session holds is the written
// file, whose root states a requested size in its unit (`180mm`, gbdraw/render/document.py), and then the px
// `viewBox` is the size.

// A length in px as an SVG root writes it: a number with no unit or with `px`.
const PX_LENGTH = /^\s*([0-9.eE+-]+)\s*(?:px)?\s*$/;

/** @param {string | null} text */
const pxLength = (text) => {
  const match = PX_LENGTH.exec(text ?? '');
  return match ? parseFloat(match[1]) : NaN;
};

/**
 * @param {Element | null} svg
 * @returns {{ width: number, height: number } | null}
 */
export const getSvgDimensions = (svg) => {
  if (!svg) return null;
  let w = pxLength(svg.getAttribute('width'));
  let h = pxLength(svg.getAttribute('height'));
  if (!w || !h) {
    const viewBox = svg.getAttribute('viewBox');
    if (!viewBox) return null;
    const parts = viewBox.trim().split(/[\s,]+/);
    if (parts.length < 4) return null;
    w = parseFloat(parts[2]);
    h = parseFloat(parts[3]);
  }
  if (!w || !h) return null;
  return { width: w, height: h };
};
