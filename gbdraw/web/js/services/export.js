import { state } from '../state.js';
import { setDpiInPng } from '../utils/png.js';

export const downloadSVG = () => {
  if (!state.svgContent.value) return;
  const blob = new Blob([state.svgContent.value], { type: 'image/svg+xml' });
  const url = URL.createObjectURL(blob);
  const a = document.createElement('a');
  a.href = url;
  a.download = state.results.value[state.selectedResultIndex.value].name;
  a.click();
  URL.revokeObjectURL(url);
};

export const downloadPNG = () => {
  if (!state.svgContent.value) return;
  const svgEl = document.createElement('div');
  svgEl.innerHTML = state.svgContent.value;
  const svg = svgEl.querySelector('svg');
  if (!svg) return;
  let w = parseFloat(svg.getAttribute('width'));
  let h = parseFloat(svg.getAttribute('height'));
  if (!w || !h) {
    const parts = svg.getAttribute('viewBox').split(' ');
    w = parseFloat(parts[2]);
    h = parseFloat(parts[3]);
  }
  const canvas = document.createElement('canvas');
  const dpi = parseInt(state.downloadDpi.value, 10);
  const scale = dpi / 96;
  canvas.width = w * scale;
  canvas.height = h * scale;
  const ctx = canvas.getContext('2d');
  const img = new Image();
  const svgData = new XMLSerializer().serializeToString(svg);
  const blob = new Blob([svgData], { type: 'image/svg+xml;charset=utf-8' });
  const url = URL.createObjectURL(blob);
  img.onload = () => {
    ctx.fillStyle = 'white';
    ctx.fillRect(0, 0, canvas.width, canvas.height);
    ctx.drawImage(img, 0, 0, canvas.width, canvas.height);
    canvas.toBlob(async (pngBlob) => {
      const fixedBlob = await setDpiInPng(pngBlob, dpi);
      const downloadUrl = URL.createObjectURL(fixedBlob);
      const link = document.createElement('a');
      let fname = state.results.value[state.selectedResultIndex.value].name;
      fname = fname.replace(/\.svg$/i, '.png');
      link.download = fname;
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
  if (!state.svgContent.value) return;

  // 1. Create temporary SVG element to read dimensions
  const svgEl = document.createElement('div');
  svgEl.innerHTML = state.svgContent.value;
  const svg = svgEl.querySelector('svg');
  if (!svg) return;

  // 2. Determine Width/Height (Use attributes or viewBox)
  let w = parseFloat(svg.getAttribute('width'));
  let h = parseFloat(svg.getAttribute('height'));
  if (!w || !h) {
    const parts = svg.getAttribute('viewBox').split(' ');
    w = parseFloat(parts[2]);
    h = parseFloat(parts[3]);
  }

  // 3. Initialize jsPDF (Use 'pt' as unit for vector fidelity)
  const { jsPDF } = window.jspdf;
  const doc = new jsPDF({
    orientation: w > h ? 'l' : 'p',
    unit: 'pt',
    format: [w, h]
  });

  // 4. Convert SVG to PDF
  await doc.svg(svg, {
    x: 0,
    y: 0,
    width: w,
    height: h
  });

  // 5. Save File
  let fname = state.results.value[state.selectedResultIndex.value].name;
  fname = fname.replace(/\.svg$/i, '.pdf');
  doc.save(fname);
};
