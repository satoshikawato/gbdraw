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

  // 2. Determine Width/Height (Use attributes or viewBox)
  const dims = getSvgDimensions(svg);
  if (!dims) return;

  // 3. Initialize jsPDF (Use 'pt' as unit for vector fidelity)
  const { jsPDF } = window.jspdf;
  const doc = new jsPDF({
    orientation: dims.width > dims.height ? 'l' : 'p',
    unit: 'pt',
    format: [dims.width, dims.height]
  });

  // 4. Convert SVG to PDF
  await doc.svg(svg, {
    x: 0,
    y: 0,
    width: dims.width,
    height: dims.height
  });

  // 5. Save File
  doc.save(getDownloadName('pdf'));
};
