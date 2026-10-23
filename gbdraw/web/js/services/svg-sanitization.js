// @ts-check
export const SVG_SANITIZE_OPTIONS = Object.freeze({
  USE_PROFILES: { svg: true },
  ADD_TAGS: [
    'use',
    'g',
    'defs',
    'linearGradient',
    'radialGradient',
    'stop',
    'path',
    'rect',
    'circle',
    'line',
    'polyline',
    'polygon',
    'text',
    'tspan'
  ],
  ADD_ATTR: [
    'xlink:href',
    'href',
    'id',
    'class',
    'data-definition-line-kind',
    'data-gbdraw-definition-part',
    'data-gbdraw-feature-id',
    'data-gbdraw-rendered-feature-id',
    'data-gbdraw-stable-feature-id',
    'data-gbdraw-feature-part',
    'data-gbdraw-feature-offset-y',
    'data-gbdraw-auto-feature-underlay',
    'data-gbdraw-slot-id',
    'data-gbdraw-slot-renderer',
    'data-gbdraw-annotation-clipped',
    'data-legend-key',
    'data-legend-name',
    'data-legend-owner',
    'data-legend-merged',
    'data-label-key',
    'data-label-feature-id',
    'data-gbdraw-label-binding-schema',
    'data-label-source-text',
    'data-label-editable',
    'data-gbdraw-label-visibility-preview',
    'data-collinearity-block-id',
    'data-collinearity-block-kind',
    'data-collinearity-orientation',
    'data-collinearity-block-evalue',
    'data-collinearity-color-mode',
    'data-edge-kind',
    'data-group-kind',
    'data-group-scope',
    'data-collinear-group-scope',
    'data-ortholog-path-id',
    'data-orthogroup-id',
    'data-rbh-orthogroup-id',
    'data-render-role',
    'data-query-protein-id',
    'data-subject-protein-id',
    'data-query-feature-svg-id',
    'data-subject-feature-svg-id',
    'data-query-stable-feature-svg-id',
    'data-subject-stable-feature-svg-id',
    'data-query-unit-id',
    'data-subject-unit-id',
    'data-query-orthogroup-assignment-reason',
    'data-subject-orthogroup-assignment-reason',
    'data-query-orthogroup-confidence',
    'data-subject-orthogroup-confidence',
    'data-query-orthogroup-member-count',
    'data-subject-orthogroup-member-count',
    'data-query-orthogroup-representative',
    'data-subject-orthogroup-representative',
    'data-query-orthogroup-role',
    'data-subject-orthogroup-role',
    'data-gbdraw-match-id',
    'data-gbdraw-pairwise-match-id',
    'data-match-kind',
    'data-query-record-index',
    'data-subject-record-index',
    'data-query-row',
    'data-subject-row',
    'data-query-record-id',
    'data-subject-record-id',
    'data-qstart',
    'data-qend',
    'data-sstart',
    'data-send',
    'data-alignment-length',
    'data-mismatches',
    'data-gap-opens',
    'data-collinearity-block-score',
    'data-collinearity-anchor-index',
    'data-collinearity-anchor-count',
    'data-query-locus-id',
    'data-subject-locus-id',
    'data-query-display-name',
    'data-subject-display-name',
    'data-pairwise-match-style',
    'data-identity-factor',
    'data-source-index',
    'data-track-index',
    'data-track-label',
    'data-track-color',
    'data-reference-side',
    'data-identity',
    'data-query',
    'data-subject',
    'data-evalue',
    'data-bitscore',
    'data-orientation',
    'data-reference-record-id',
    'data-gbdraw-annotation-id',
    'data-gbdraw-annotation-set-id',
    'data-gbdraw-annotation-track-id',
    'data-gbdraw-record-id',
    'data-gbdraw-record-index',
    'data-record-column',
    'data-record-index',
    'data-record-key',
    'data-record-row',
    'data-record-translation-x',
    'data-record-translation-y',
    'data-gbdraw-annotation-mark',
    'data-gbdraw-annotation-label',
    'data-gbdraw-region-label-for',
    'data-gbdraw-role',
    'data-gbdraw-composition-schema',
    'data-gbdraw-composition',
    'data-gbdraw-composition-role',
    'data-gbdraw-orientation',
    'data-annotation-set-id',
    'data-annotation-track-id',
    'data-annotation-record-id',
    'data-annotation-record-index',
    'data-annotation-mark',
    'data-annotation-label',
    'fill',
    'fill-opacity',
    'stroke',
    'stroke-width',
    'stroke-opacity',
    'stroke-dasharray',
    'stroke-linecap',
    'stroke-linejoin',
    'd',
    'x',
    'y',
    'x1',
    'y1',
    'x2',
    'y2',
    'cx',
    'cy',
    'r',
    'rx',
    'ry',
    'width',
    'height',
    'transform',
    'viewBox',
    'preserveAspectRatio',
    'font-family',
    'font-size',
    'font-weight',
    'text-anchor',
    'dominant-baseline',
    'writing-mode',
    'letter-spacing',
    'display',
    'pointer-events'
  ],
  FORBID_TAGS: [
    'style',
    'script',
    'foreignObject',
    'iframe',
    'embed',
    'object',
    'animate',
    'set',
    'animateTransform',
    'image'
  ],
  FORBID_ATTR: ['name', 'onload', 'onclick', 'onmouseover', 'onfocus', 'onerror']
});

/**
 * The part of DOMPurify that this module calls.
 * @typedef {{ sanitize?: (source: string, options: Readonly<Record<string, unknown>>) => unknown }} SvgSanitizer
 */

/**
 * @param {unknown} rawSvg The SVG text; anything else is rejected.
 * @param {SvgSanitizer} [sanitizer] A DOMPurify-compatible sanitizer.
 */
export const sanitizeSvgContent = (
  rawSvg,
  sanitizer = globalThis.DOMPurify
) => {
  if (!sanitizer || typeof sanitizer.sanitize !== 'function') {
    throw new Error('SVG sanitization is unavailable.');
  }
  if (typeof rawSvg !== 'string' || !rawSvg.trim()) {
    throw new Error('The diagram engine returned empty SVG content.');
  }
  const sanitized = sanitizer.sanitize(rawSvg, SVG_SANITIZE_OPTIONS);
  const content = typeof sanitized === 'string'
    ? sanitized
    : String(sanitized ?? '');
  if (!content.trim()) {
    throw new Error('The diagram engine returned SVG content that could not be sanitized.');
  }
  return content;
};

const SVG_FRAGMENT_OPEN = '<svg xmlns="http://www.w3.org/2000/svg">';
const SVG_END = '</svg>';

/**
 * A fragment of SVG elements inside an SVG root, as the sanitizer reads a
 * whole document. Pure: calls nothing.
 * @param {string} fragment
 */
export const wrapSvgFragment = (fragment) => `${SVG_FRAGMENT_OPEN}${fragment}${SVG_END}`;

/**
 * The fragment of a sanitized wrapped fragment: the text between `prefix`
 * (the sanitized form of the wrapper's start tag) and the root end tag, or
 * null when either is missing. Pure: calls nothing.
 * @param {string} sanitized
 * @param {string} prefix
 * @returns {string | null}
 */
export const unwrapSvgFragment = (sanitized, prefix) => (
  sanitized.length >= prefix.length + SVG_END.length
  && sanitized.startsWith(prefix)
  && sanitized.endsWith(SVG_END)
    ? sanitized.slice(prefix.length, sanitized.length - SVG_END.length)
    : null
);
