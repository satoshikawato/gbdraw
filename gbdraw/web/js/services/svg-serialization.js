const TRANSIENT_PREVIEW_CLASSES = Object.freeze([
  'gbdraw-preview-feature-search-match',
  'gbdraw-preview-feature-search-active-match',
  'gbdraw-preview-feature-search-dimmed',
  'gbdraw-feature-selected',
  'gbdraw-feature-selection-anchor',
  'gbdraw-feature-selection-candidate',
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
  TRANSIENT_PREVIEW_CLASSES.forEach((className) => {
    svg.querySelectorAll(`.${className}`).forEach((element) => removeClassToken(element, className));
  });
  if (stripCursor) stripEditorOnlyCursorStyles(svg);
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

export const ensureSvgDefs = (svg) => {
  let defs = svg.querySelector('defs');
  if (!defs) {
    defs = document.createElementNS('http://www.w3.org/2000/svg', 'defs');
    svg.insertBefore(defs, svg.firstChild);
  }
  return defs;
};
