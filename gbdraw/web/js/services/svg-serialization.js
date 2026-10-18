// @ts-check
/** @import { FeatureCatalog } from './feature-catalog.js' */
/** @import { FeatureOverrideDraft } from './feature-placement.js' */

/**
 * A Vue ref, or the `{ value }` cell that stands in for one (a Gallery
 * Session's request state), as the export capture and the request writer
 * unwrap it.
 * @template T
 * @typedef {{ value: T }} ExportRef
 */

/**
 * One mode's drawing, as far as the interactive export reads its edits.
 * @typedef {object} SvgExportDrawing
 * @property {FeatureOverrideDraft} [featureOverrides]
 * @property {Record<string, string>} labelTextBulkOverrides
 * @property {Record<string, string>} orthogroupNameOverrides
 * @property {Record<string, string>} orthogroupDescriptionOverrides
 */

/**
 * What the export capture reads of the app state: the refs it unwraps.
 * @typedef {object} SvgExportState
 * @property {ExportRef<number | string>} selectedResultIndex
 * @property {ExportRef<ReadonlyArray<{ name?: string } | undefined> | null>} results
 * @property {ExportRef<Element | null>} svgContainer
 * @property {ExportRef<string | null>} svgContent
 * @property {ExportRef<number | string>} downloadDpi
 * @property {ExportRef<boolean>} [richFeaturePopup]
 * @property {ExportRef<FeatureCatalog | null>} [featureCatalog]
 * @property {ExportRef<ReadonlyArray<{ featureId: string, text: string, sourceText: string }>>} [editableLabels]
 */

const TRANSIENT_PREVIEW_CLASSES = Object.freeze([
  'gbdraw-preview-layout-target',
  'gbdraw-preview-feature-search-match',
  'gbdraw-preview-feature-search-active-match',
  'gbdraw-preview-feature-search-dimmed',
  'gbdraw-preview-feature-search-results-active',
  'gbdraw-preview-feature-search-updating',
  'gbdraw-feature-selected',
  'gbdraw-feature-selection-anchor',
  'gbdraw-feature-selection-candidate',
  'gbdraw-match-pending',
  'gbdraw-match-selected',
  'feature-selection-marquee',
  'feature-selection-status'
]);

/**
 * @param {Element | null | undefined} element
 * @param {string} token
 * @param {unknown} enabled
 */
export const setClassToken = (element, token, enabled) => {
  if (!element) return;
  if (element.classList?.toggle) {
    element.classList.toggle(token, Boolean(enabled));
    return;
  }
  const existing = String(element.getAttribute('class') || '').trim();
  const tokens = existing ? existing.split(/\s+/).filter(Boolean) : [];
  const nextTokens = tokens.filter((entry) => entry !== token);
  if (enabled) nextTokens.push(token);
  if (nextTokens.length) {
    element.setAttribute('class', nextTokens.join(' '));
  } else {
    element.removeAttribute('class');
  }
};

/**
 * @param {Element | null | undefined} element
 * @param {string} token
 */
const removeClassToken = (element, token) => {
  if (!element) return;
  if (element.classList?.remove) {
    element.classList.remove(token);
    if (element.classList.length === 0) element.removeAttribute('class');
    return;
  }
  const tokens = String(element.getAttribute('class') || '').split(/\s+/).filter((entry) => entry && entry !== token);
  if (tokens.length) {
    element.setAttribute('class', tokens.join(' '));
  } else {
    element.removeAttribute('class');
  }
};

// Also drops an inline style that Preview feedback (a cursor, a hover, a drag)
// left empty when it was removed.
/** @param {Element | null | undefined} svg */
const stripEditorOnlyCursorStyles = (svg) => {
  if (!svg) return;
  svg.querySelectorAll('[style]').forEach((element) => {
    if (/\bcursor\s*:/i.test(element.getAttribute('style') || '')) {
      // `[style]` selects HTML and SVG elements, which carry inline styles.
      /** @type {HTMLElement | SVGElement} */ (element).style.removeProperty('cursor');
    }
    if (!element.getAttribute('style')?.trim()) {
      element.removeAttribute('style');
    }
  });
};

/**
 * @param {Element | null | undefined} svg
 * @param {{ stripCursor?: boolean }} [options]
 */
export const stripTransientPreviewState = (svg, { stripCursor = true } = {}) => {
  if (!svg) return;
  /** @param {Element} element */
  const stripClasses = (element) => {
    TRANSIENT_PREVIEW_CLASSES.forEach((className) => removeClassToken(element, className));
  };
  stripClasses(svg);
  svg.querySelectorAll(TRANSIENT_PREVIEW_CLASSES.map((className) => `.${className}`).join(','))
    .forEach(stripClasses);
  if (stripCursor) {
    stripEditorOnlyCursorStyles(svg);
    svg.querySelectorAll('[data-gbdraw-pairwise-match-id][role="button"][tabindex="0"]')
      .forEach((element) => {
        if (!/^Pairwise match \d+$/.test(element.getAttribute('aria-label') || '')) return;
        element.removeAttribute('role');
        element.removeAttribute('tabindex');
        element.removeAttribute('aria-label');
      });
    svg.querySelectorAll('[data-gbdraw-preview-item]').forEach((element) => {
      ['data-gbdraw-preview-item', 'role', 'tabindex', 'aria-haspopup'].forEach((name) => element.removeAttribute(name));
    });
  }
};

// The preview feature-search classes (app/feature-search/preview-svg.js) that
// an export strips from its clone.
export const PREVIEW_FEATURE_SEARCH_MATCH_CLASS = 'gbdraw-preview-feature-search-match';
export const PREVIEW_FEATURE_SEARCH_ACTIVE_CLASS = 'gbdraw-preview-feature-search-active-match';
export const PREVIEW_FEATURE_SEARCH_DIMMED_CLASS = 'gbdraw-preview-feature-search-dimmed';
export const PREVIEW_FEATURE_SEARCH_ROOT_ACTIVE_CLASS = 'gbdraw-preview-feature-search-results-active';
export const PREVIEW_FEATURE_SEARCH_ROOT_UPDATING_CLASS = 'gbdraw-preview-feature-search-updating';

export const PREVIEW_FEATURE_SEARCH_CLASSES = Object.freeze([
  PREVIEW_FEATURE_SEARCH_MATCH_CLASS,
  PREVIEW_FEATURE_SEARCH_ACTIVE_CLASS,
  PREVIEW_FEATURE_SEARCH_DIMMED_CLASS
]);

/** @param {Element | null | undefined} svg */
export const stripPreviewFeatureSearchClasses = (svg) => {
  if (!svg) return;
  setClassToken(svg, PREVIEW_FEATURE_SEARCH_ROOT_ACTIVE_CLASS, false);
  setClassToken(svg, PREVIEW_FEATURE_SEARCH_ROOT_UPDATING_CLASS, false);
  PREVIEW_FEATURE_SEARCH_CLASSES.forEach((className) => {
    svg.querySelectorAll(`.${className}`).forEach((element) => {
      setClassToken(element, className, false);
    });
  });
};

/**
 * @param {Element | null | undefined} svg
 * @param {{ stripCursor?: boolean }} [options]
 */
export const serializeCleanSvg = (svg, options = {}) => {
  if (!svg) return '';
  // A deep clone of an element is an element.
  const clone = /** @type {Element} */ (svg.cloneNode(true));
  stripTransientPreviewState(clone, options);
  if (!clone.getAttribute('xmlns')) {
    clone.setAttribute('xmlns', 'http://www.w3.org/2000/svg');
  }
  if (!clone.getAttribute('xmlns:xlink')) {
    clone.setAttribute('xmlns:xlink', 'http://www.w3.org/1999/xlink');
  }
  return new XMLSerializer().serializeToString(clone);
};

// The standalone popup reads label text by rendered ID; each per-feature label
// edit of the drawing of the Result's mode (identity-keyed, design Q4 6.4, R2)
// is projected onto the rendered IDs of the exported Result's catalog item.
/**
 * @param {SvgExportState} state
 * @param {SvgExportDrawing} drawing
 * @param {number} resultIndex
 */
const renderedLabelTextOverrides = (state, drawing, resultIndex) => {
  /** @type {Record<string, string>} */
  const overrides = {};
  const textByFeature = new Map(Object.values(drawing.featureOverrides || {})
    .filter((row) => typeof row?.labelText === 'string' && row.labelText)
    .map((row) => [`${row.recordKey}\u0000${row.biologicalFeatureId}`, row.labelText]));
  const item = state.featureCatalog?.value?.items?.[resultIndex];
  (Array.isArray(item?.features) ? item.features : []).forEach((feature) => {
    const text = textByFeature.get(`${feature?.recordKey}\u0000${feature?.biologicalFeatureId}`);
    if (text) overrides[feature.svgId] = text;
  });
  return overrides;
};

// Capture the selected Result before any lazy export modules or libraries load.
// Its popup settings and edits are those of the drawing the Result belongs to,
// which the caller resolves (`displayedArtifactDrawing`).
/**
 * @param {SvgExportState} state
 * @param {SvgExportDrawing} drawing
 * @param {{ interactive?: boolean }} [options]
 */
export const captureSvgExport = (state, drawing, { interactive = false } = {}) => {
  const resultIndex = Number(state.selectedResultIndex.value);
  const name = state.results.value?.[resultIndex]?.name || 'gbdraw.svg';
  return {
    // A deep clone of the preview's <svg> element is an <svg> element.
    svg: /** @type {SVGSVGElement | null} */ (state.svgContainer.value?.querySelector('svg')?.cloneNode(true) || null),
    svgContent: state.svgContent.value,
    name,
    dpi: state.downloadDpi.value,
    ...(interactive ? { interactivity: {
      popupMode: state.richFeaturePopup?.value === false ? 'simple' : 'rich',
      featureCatalog: state.featureCatalog?.value,
      catalogResultIndex: resultIndex,
      catalogResultName: name,
      requireFeatureCatalog: true,
      editableLabels: (state.editableLabels?.value || []).map(({ featureId, text, sourceText }) => ({
        featureId, text, sourceText
      })),
      labelTextFeatureOverrides: renderedLabelTextOverrides(state, drawing, resultIndex),
      labelTextBulkOverrides: { ...drawing.labelTextBulkOverrides },
      orthogroupNameOverrides: { ...drawing.orthogroupNameOverrides },
      orthogroupDescriptionOverrides: { ...drawing.orthogroupDescriptionOverrides }
    } } : {})
  };
};

/**
 * @param {SVGSVGElement} svg
 * @returns {SVGDefsElement}
 */
export const ensureSvgDefs = (svg) => {
  let defs = svg.querySelector('defs');
  if (!defs) {
    defs = document.createElementNS('http://www.w3.org/2000/svg', 'defs');
    svg.insertBefore(defs, svg.firstChild);
  }
  return defs;
};
