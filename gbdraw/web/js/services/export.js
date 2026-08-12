import { state } from '../state.js';
import { setDpiInPng } from '../utils/png.js';
import { stripPreviewFeatureSearchClasses } from '../app/feature-search/preview-svg.js';
import { ensureSvgDefs, stripTransientPreviewState } from './svg-serialization.js';
import { downloadBlob } from './text-download.js';

const SVG_NS = 'http://www.w3.org/2000/svg';
const JSPDF_SCRIPT_URL = new URL('../../vendor/jspdf/jspdf.umd.min.js', import.meta.url);
const SVG2PDF_SCRIPT_URL = new URL(
  '../../vendor/svg2pdf.js/svg2pdf.umd.min.js',
  import.meta.url
);

let standaloneInteractivityPromise = null;
let pdfLibrariesPromise = null;

const loadStandaloneInteractivity = () => {
  standaloneInteractivityPromise ??= import('./standalone-interactivity.js');
  return standaloneInteractivityPromise;
};

const loadSameOriginScript = (url, label) => new Promise((resolve, reject) => {
  if (url.origin !== window.location.origin) {
    reject(new Error(`${label} must be loaded from the application origin.`));
    return;
  }

  const script = document.createElement('script');
  script.src = url.href;
  script.async = false;
  script.dataset.gbdrawExportLibrary = label;
  script.addEventListener('load', resolve, { once: true });
  script.addEventListener('error', () => {
    reject(new Error(`Failed to load the vendored ${label} library.`));
  }, { once: true });
  document.head.appendChild(script);
});

const loadPdfLibraries = () => {
  pdfLibrariesPromise ??= (async () => {
    await loadSameOriginScript(JSPDF_SCRIPT_URL, 'jsPDF');
    if (typeof window.jspdf?.jsPDF !== 'function') {
      throw new Error('The vendored jsPDF library did not initialize.');
    }

    await loadSameOriginScript(SVG2PDF_SCRIPT_URL, 'svg2pdf');
    if (typeof window.jspdf.jsPDF.API?.svg !== 'function') {
      throw new Error('The vendored svg2pdf library did not initialize.');
    }
  })();
  return pdfLibrariesPromise;
};

const getDownloadName = (extension) => {
  const baseName =
    state.results.value?.[state.selectedResultIndex.value]?.name ||
    (extension ? `gbdraw.${extension}` : 'gbdraw.svg');
  if (!extension) return baseName;
  const normalized = baseName.replace(/\.svg$/i, `.${extension}`);
  if (normalized === baseName && !baseName.toLowerCase().endsWith(`.${extension}`)) {
    return `${baseName}.${extension}`;
  }
  return normalized;
};

const cloneCurrentSvg = () => {
  const liveSvg = state.svgContainer.value?.querySelector('svg');
  if (liveSvg) {
    const clone = liveSvg.cloneNode(true);
    if (!clone.getAttribute('xmlns')) {
      clone.setAttribute('xmlns', 'http://www.w3.org/2000/svg');
    }
    if (!clone.getAttribute('xmlns:xlink')) {
      clone.setAttribute('xmlns:xlink', 'http://www.w3.org/1999/xlink');
    }
    return clone;
  }
  return getSvgFromString(state.svgContent.value);
};

const getCurrentSvgClone = ({ interactive = false } = {}) => {
  const clone = cloneCurrentSvg();
  if (!clone) {
    throw new Error('No SVG result is available for export.');
  }
  stripTransientPreviewState(clone, { stripCursor: !interactive });
  stripPreviewFeatureSearchClasses(clone);
  return clone;
};

const getCurrentSvgString = () => (
  new XMLSerializer().serializeToString(getCurrentSvgClone())
);

const getInteractiveSvgString = async () => {
  const { enrichSvgWithStandaloneInteractivity } = await loadStandaloneInteractivity();
  const clone = getCurrentSvgClone({ interactive: true });
  const resultIndex = Number(state.selectedResultIndex.value);
  const result = state.results.value?.[resultIndex] || null;
  const enriched = enrichSvgWithStandaloneInteractivity(clone, {
    popupMode: state.adv.rich_feature_popup === false ? 'simple' : 'rich',
    featureCatalog: state.featureCatalog?.value,
    catalogResultIndex: resultIndex,
    catalogResultName: result?.name,
    requireFeatureCatalog: true,
    editableLabels: state.editableLabels?.value,
    labelTextFeatureOverrides: state.labelTextFeatureOverrides,
    labelTextBulkOverrides: state.labelTextBulkOverrides,
    orthogroupNameOverrides: state.orthogroupNameOverrides,
    orthogroupDescriptionOverrides: state.orthogroupDescriptionOverrides
  });
  if (!enriched) {
    throw new Error('Interactive SVG export requires the committed feature catalog.');
  }
  return new XMLSerializer().serializeToString(clone);
};

const getSvgFromString = (svgString) => {
  if (!svgString) return null;
  const svgEl = document.createElement('div');
  svgEl.innerHTML = svgString;
  return svgEl.querySelector('svg');
};

const getSvgDimensions = (svg) => {
  if (!svg) return null;
  let w = parseFloat(svg.getAttribute('width'));
  let h = parseFloat(svg.getAttribute('height'));
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

const moveGradientsToDefs = (svg) => {
  const defs = ensureSvgDefs(svg);
  svg.querySelectorAll('linearGradient, radialGradient').forEach((gradient) => {
    if (gradient.parentNode !== defs) {
      defs.appendChild(gradient);
    }
  });
};

const copyAttributes = (target, source, attributes) => {
  attributes.forEach((attr) => {
    const value = source.getAttribute(attr);
    if (value !== null) {
      target.setAttribute(attr, value);
    }
  });
};

const applyFontAndPaintStyles = (target, source) => {
  if (!source) return;
  const style = window.getComputedStyle(source);
  [
    'font-family',
    'font-size',
    'font-weight',
    'font-style',
    'font-variant',
    'fill',
    'fill-opacity',
    'stroke',
    'stroke-opacity',
    'stroke-width'
  ].forEach((prop) => {
    const value = style.getPropertyValue(prop);
    if (value) {
      target.setAttribute(prop, value.trim());
    }
  });
};

const flattenTextPathsForPdf = (svg) => {
  const textElements = Array.from(svg.querySelectorAll('text')).filter((el) =>
    el.querySelector('textPath')
  );
  if (!textElements.length) return;

  textElements.forEach((textEl) => {
    const textPath = textEl.querySelector('textPath');
    if (!textPath) return;
    const textContent = textPath.textContent || textEl.textContent || '';
    const charCount =
      typeof textEl.getNumberOfChars === 'function' ? textEl.getNumberOfChars() : 0;
    if (!textContent || !charCount) return;

    const chars = Array.from(textContent);
    const group = document.createElementNS(SVG_NS, 'g');
    copyAttributes(group, textEl, ['transform', 'opacity', 'display', 'visibility']);
    applyFontAndPaintStyles(group, textPath);

    for (let i = 0; i < charCount && i < chars.length; i += 1) {
      const char = chars[i];
      if (char === '\n' || char === '\r') continue;
      let pos;
      try {
        pos = textEl.getStartPositionOfChar(i);
      } catch (error) {
        pos = null;
      }
      if (!pos) continue;
      const rotation = textEl.getRotationOfChar(i) || 0;
      const charText = document.createElementNS(SVG_NS, 'text');
      charText.textContent = char;
      charText.setAttribute('x', pos.x);
      charText.setAttribute('y', pos.y);
      charText.setAttribute('text-anchor', 'start');
      charText.setAttribute('dominant-baseline', 'alphabetic');
      if (rotation) {
        charText.setAttribute('transform', `rotate(${rotation} ${pos.x} ${pos.y})`);
      }
      group.appendChild(charText);
    }

    if (textEl.parentNode) {
      textEl.parentNode.replaceChild(group, textEl);
    }
  });
};

const normalizeTextBaselinesForPdf = (svg) => {
  svg.querySelectorAll('text').forEach((textEl) => {
    const baseline = textEl.getAttribute('dominant-baseline');
    if (!baseline || baseline === 'alphabetic' || baseline === 'auto') return;

    let before;
    try {
      before = textEl.getBBox();
    } catch (error) {
      return;
    }

    textEl.setAttribute('dominant-baseline', 'alphabetic');

    let after;
    try {
      after = textEl.getBBox();
    } catch (error) {
      textEl.setAttribute('dominant-baseline', baseline);
      return;
    }

    const dy = before.y - after.y;
    if (!Number.isFinite(dy) || Math.abs(dy) < 0.01) {
      textEl.removeAttribute('dominant-baseline');
      return;
    }

    const yAttr = textEl.getAttribute('y');
    if (yAttr) {
      const yValues = yAttr.split(/[\s,]+/).map((value) => parseFloat(value));
      const adjusted = yValues.map((value) =>
        Number.isFinite(value) ? value + dy : value
      );
      textEl.setAttribute('y', adjusted.join(' '));
    } else {
      const dyAttr = textEl.getAttribute('dy');
      const dyValue = dyAttr ? parseFloat(dyAttr) : 0;
      const nextDy = Number.isFinite(dyValue) ? dyValue + dy : dy;
      textEl.setAttribute('dy', nextDy);
    }

    textEl.removeAttribute('dominant-baseline');
  });
};

const prepareSvgForPdf = (svg) => {
  if (!svg) {
    throw new Error('The current SVG could not be prepared for PDF export.');
  }
  const clone = svg.cloneNode(true);
  if (!clone.getAttribute('xmlns')) {
    clone.setAttribute('xmlns', SVG_NS);
  }
  if (!clone.getAttribute('xmlns:xlink')) {
    clone.setAttribute('xmlns:xlink', 'http://www.w3.org/1999/xlink');
  }

  const staging = document.createElement('div');
  staging.style.position = 'absolute';
  staging.style.left = '-10000px';
  staging.style.top = '-10000px';
  staging.style.visibility = 'hidden';
  staging.appendChild(clone);
  document.body.appendChild(staging);

  try {
    moveGradientsToDefs(clone);
    flattenTextPathsForPdf(clone);
    normalizeTextBaselinesForPdf(clone);
  } catch (error) {
    staging.remove();
    throw error;
  }

  return {
    svg: clone,
    cleanup: () => staging.remove()
  };
};

const downloadSvgString = (svgString, filename) => {
  if (!svgString) {
    throw new Error('The SVG export produced no content.');
  }
  downloadBlob(new Blob([svgString], { type: 'image/svg+xml' }), filename);
};

export const downloadSVG = () => {
  const svgString = getCurrentSvgString();
  downloadSvgString(svgString, getDownloadName('svg'));
};

export const downloadInteractiveSVG = async () => {
  const svgString = await getInteractiveSvgString();
  downloadSvgString(svgString, getDownloadName('interactive.svg'));
};

export const downloadPNG = async () => {
  const svgString = getCurrentSvgString();
  const svg = getSvgFromString(svgString);
  if (!svg) {
    throw new Error('The current SVG could not be parsed for PNG export.');
  }
  const dims = getSvgDimensions(svg);
  if (!dims) {
    throw new Error('The current SVG has no usable dimensions for PNG export.');
  }
  const canvas = document.createElement('canvas');
  const dpi = parseInt(state.downloadDpi.value, 10);
  if (!Number.isFinite(dpi) || dpi <= 0) {
    throw new Error('The selected PNG DPI is invalid.');
  }
  const scale = dpi / 96;
  canvas.width = dims.width * scale;
  canvas.height = dims.height * scale;
  const ctx = canvas.getContext('2d');
  if (!ctx) {
    throw new Error('The browser could not initialize PNG conversion.');
  }
  const img = new Image();
  const blob = new Blob([svgString], { type: 'image/svg+xml;charset=utf-8' });
  const url = URL.createObjectURL(blob);
  try {
    await new Promise((resolve, reject) => {
      img.onload = resolve;
      img.onerror = () => reject(new Error('The browser could not load the SVG for PNG export.'));
      img.src = url;
    });
    ctx.fillStyle = 'white';
    ctx.fillRect(0, 0, canvas.width, canvas.height);
    ctx.drawImage(img, 0, 0, canvas.width, canvas.height);
    const pngBlob = await new Promise((resolve, reject) => {
      canvas.toBlob((value) => {
        if (value) {
          resolve(value);
        } else {
          reject(new Error('The browser produced no PNG export data.'));
        }
      }, 'image/png');
    });
    const fixedBlob = await setDpiInPng(pngBlob, dpi);
    downloadBlob(fixedBlob, getDownloadName('png'));
  } finally {
    URL.revokeObjectURL(url);
  }
};

export const downloadPDF = async () => {
  const svgString = getCurrentSvgString();

  // 1. Create temporary SVG element to read dimensions
  const svg = getSvgFromString(svgString);
  if (!svg) {
    throw new Error('The current SVG could not be parsed for PDF export.');
  }
  const prepared = prepareSvgForPdf(svg);
  const { svg: pdfSvg, cleanup } = prepared;

  try {
    // 2. Determine Width/Height (Use attributes or viewBox)
    const dims = getSvgDimensions(pdfSvg);
    if (!dims) {
      throw new Error('The current SVG has no usable dimensions for PDF export.');
    }

    // 3. Load and initialize jsPDF/svg2pdf only for this export format.
    await loadPdfLibraries();
    const { jsPDF } = window.jspdf;
    const doc = new jsPDF({
      orientation: dims.width > dims.height ? 'l' : 'p',
      unit: 'pt',
      format: [dims.width, dims.height]
    });

    // 4. Convert SVG to PDF
    await doc.svg(pdfSvg, {
      x: 0,
      y: 0,
      width: dims.width,
      height: dims.height
    });
    // 5. Save File
    doc.save(getDownloadName('pdf'));
  } finally {
    cleanup();
  }
};
