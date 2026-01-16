import { state } from '../state.js';
import { setDpiInPng } from '../utils/png.js';

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

const getCurrentSvgString = () => {
  const liveSvg = state.svgContainer.value?.querySelector('svg');
  if (liveSvg) {
    const clone = liveSvg.cloneNode(true);
    if (!clone.getAttribute('xmlns')) {
      clone.setAttribute('xmlns', 'http://www.w3.org/2000/svg');
    }
    if (!clone.getAttribute('xmlns:xlink')) {
      clone.setAttribute('xmlns:xlink', 'http://www.w3.org/1999/xlink');
    }
    return new XMLSerializer().serializeToString(clone);
  }
  return state.svgContent.value;
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

const SVG_NS = 'http://www.w3.org/2000/svg';

const ensureSvgDefs = (svg) => {
  let defs = svg.querySelector('defs');
  if (!defs) {
    defs = document.createElementNS(SVG_NS, 'defs');
    svg.insertBefore(defs, svg.firstChild);
  }
  return defs;
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
  if (!svg) return null;
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

  moveGradientsToDefs(clone);
  flattenTextPathsForPdf(clone);
  normalizeTextBaselinesForPdf(clone);

  return {
    svg: clone,
    cleanup: () => staging.remove()
  };
};

export const downloadSVG = () => {
  const svgString = getCurrentSvgString();
  if (!svgString) return;
  const blob = new Blob([svgString], { type: 'image/svg+xml' });
  const url = URL.createObjectURL(blob);
  const a = document.createElement('a');
  a.href = url;
  a.download = getDownloadName('svg');
  a.click();
  URL.revokeObjectURL(url);
};

export const downloadPNG = () => {
  const svgString = getCurrentSvgString();
  if (!svgString) return;
  const svg = getSvgFromString(svgString);
  if (!svg) return;
  const dims = getSvgDimensions(svg);
  if (!dims) return;
  const canvas = document.createElement('canvas');
  const dpi = parseInt(state.downloadDpi.value, 10);
  const scale = dpi / 96;
  canvas.width = dims.width * scale;
  canvas.height = dims.height * scale;
  const ctx = canvas.getContext('2d');
  const img = new Image();
  const blob = new Blob([svgString], { type: 'image/svg+xml;charset=utf-8' });
  const url = URL.createObjectURL(blob);
  img.onload = () => {
    ctx.fillStyle = 'white';
    ctx.fillRect(0, 0, canvas.width, canvas.height);
    ctx.drawImage(img, 0, 0, canvas.width, canvas.height);
    canvas.toBlob(async (pngBlob) => {
      const fixedBlob = await setDpiInPng(pngBlob, dpi);
      const downloadUrl = URL.createObjectURL(fixedBlob);
      const link = document.createElement('a');
      link.download = getDownloadName('png');
      link.href = downloadUrl;
      link.click();
      URL.revokeObjectURL(downloadUrl);
      URL.revokeObjectURL(url);
    }, 'image/png');
  };
  img.onerror = () => {
    console.error('Failed to load SVG for PNG conversion');
    URL.revokeObjectURL(url);
  };
  img.src = url;
};

export const downloadPDF = async () => {
  const svgString = getCurrentSvgString();
  if (!svgString) return;

  // 1. Create temporary SVG element to read dimensions
  const svg = getSvgFromString(svgString);
  if (!svg) return;
  const prepared = prepareSvgForPdf(svg);
  if (!prepared) return;
  const { svg: pdfSvg, cleanup } = prepared;

  // 2. Determine Width/Height (Use attributes or viewBox)
  const dims = getSvgDimensions(pdfSvg);
  if (!dims) {
    cleanup();
    return;
  }

  // 3. Initialize jsPDF (Use 'pt' as unit for vector fidelity)
  const { jsPDF } = window.jspdf;
  const doc = new jsPDF({
    orientation: dims.width > dims.height ? 'l' : 'p',
    unit: 'pt',
    format: [dims.width, dims.height]
  });

  // 4. Convert SVG to PDF
  try {
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
