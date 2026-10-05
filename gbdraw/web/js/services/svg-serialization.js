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

const stripEditorOnlyCursorStyles = (svg) => {
  if (!svg) return;
  svg.querySelectorAll('[style]').forEach((element) => {
    const style = element.getAttribute('style');
    if (!style || !/\bcursor\s*:/i.test(style)) return;
    element.style.removeProperty('cursor');
    if (!element.getAttribute('style')?.trim()) {
      element.removeAttribute('style');
    }
  });
};

export const stripTransientPreviewState = (svg, { stripCursor = true } = {}) => {
  if (!svg) return;
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
  }
};

export const serializeCleanSvg = (svg, options = {}) => {
  if (!svg) return '';
  const clone = svg.cloneNode(true);
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
// edit of the Result's mode (identity-keyed, design Q4 6.4, R2) is projected
// onto the rendered IDs of the exported Result's catalog item.
const renderedLabelTextOverrides = (state, resultIndex) => {
  const overrides = {};
  const mode = state.generatedMode?.value;
  const textByFeature = new Map(Object.values(state.featureOverrides || {})
    .filter((row) => row?.scope === mode && typeof row.labelText === 'string' && row.labelText)
    .map((row) => [`${row.recordKey}\u0000${row.biologicalFeatureId}`, row.labelText]));
  const item = state.featureCatalog?.value?.items?.[resultIndex];
  (Array.isArray(item?.features) ? item.features : []).forEach((feature) => {
    const text = textByFeature.get(`${feature?.recordKey}\u0000${feature?.biologicalFeatureId}`);
    if (text) overrides[feature.svgId] = text;
  });
  return overrides;
};

// Capture the selected Result before any lazy export modules or libraries load.
export const captureSvgExport = (state, { interactive = false } = {}) => {
  const resultIndex = Number(state.selectedResultIndex.value);
  const name = state.results.value?.[resultIndex]?.name || 'gbdraw.svg';
  return {
    svg: state.svgContainer.value?.querySelector('svg')?.cloneNode(true) || null,
    svgContent: state.svgContent.value,
    name,
    dpi: state.downloadDpi.value,
    ...(interactive ? { interactivity: {
      popupMode: state.adv.rich_feature_popup === false ? 'simple' : 'rich',
      featureCatalog: state.featureCatalog?.value,
      catalogResultIndex: resultIndex,
      catalogResultName: name,
      requireFeatureCatalog: true,
      editableLabels: (state.editableLabels?.value || []).map(({ featureId, text, sourceText }) => ({
        featureId, text, sourceText
      })),
      labelTextFeatureOverrides: renderedLabelTextOverrides(state, resultIndex),
      labelTextBulkOverrides: { ...state.labelTextBulkOverrides },
      orthogroupNameOverrides: { ...state.orthogroupNameOverrides },
      orthogroupDescriptionOverrides: { ...state.orthogroupDescriptionOverrides }
    } } : {})
  };
};

export const ensureSvgDefs = (svg) => {
  let defs = svg.querySelector('defs');
  if (!defs) {
    defs = document.createElementNS('http://www.w3.org/2000/svg', 'defs');
    svg.insertBefore(defs, svg.firstChild);
  }
  return defs;
};
