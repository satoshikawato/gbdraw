// @ts-check
import { t } from '../utils/i18n.js';
import { setDpiInPng } from '../utils/png.js';
import {
  ensureSvgDefs,
  stripPreviewFeatureSearchClasses,
  stripTransientPreviewState
} from './svg-serialization.js';
import { stripResultBaseAttributes } from './result-paint-bases.js';
import { downloadBlob } from './text-download.js';
import { preparePdfFonts } from './pdf-fonts.js';
import { diagnosticError } from '../utils/error-normalization.js';
import { parseMarkupInertly } from '../utils/inert-markup.js';
import { getSvgDimensions } from './svg-dimensions.js';
/** @import { ResolvePdfFonts } from './pdf-fonts.js' */

const SVG_NS = 'http://www.w3.org/2000/svg';
const XML_NS = 'http://www.w3.org/XML/1998/namespace';
// SVG lengths are CSS px: 96 px = 72 pt = 1 in (PNG DPI and PDF points share it).
const CSS_PX_PER_INCH = 96;
const PT_PER_CSS_PX = 72 / CSS_PX_PER_INCH;
const JSPDF_SCRIPT_URL = new URL('../../vendor/jspdf/jspdf.umd.min.js', import.meta.url);
const SVG2PDF_SCRIPT_URL = new URL(
  '../../vendor/svg2pdf.js/svg2pdf.umd.min.js',
  import.meta.url
);

/** @type {Promise<typeof import('./standalone-interactivity.js')> | null} */
let standaloneInteractivityPromise = null;
/** @type {Promise<void> | null} */
let pdfLibrariesPromise = null;

const loadStandaloneInteractivity = () => {
  standaloneInteractivityPromise ??= import('./standalone-interactivity.js');
  return standaloneInteractivityPromise;
};

/**
 * @param {URL} url
 * @param {string} label
 * @returns {Promise<Event>} Resolves with the script's load event.
 */
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
    script.remove();
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
  })().catch((error) => {
    pdfLibrariesPromise = null;
    throw error;
  });
  return pdfLibrariesPromise;
};

/** @import { ExportInteractivity } from './svg-serialization.js' */

/**
 * The fonts an exported SVG carries (custom fonts D-05, D-13): Python's
 * `svg_font_embedding` in the Worker. `block` is the
 * `<defs id="gbdraw-embedded-fonts">` block Python builds (`embedded_fonts_block`,
 * "" for none); `refusedFamilies` names the fonts left out because their
 * license forbids embedding. `kind`: `svg` for SVG and Interactive SVG files,
 * `raster` for the SVG the PNG export draws (every custom font, whatever its
 * license: rasterizing is not embedding).
 * @typedef {{ block: string, refusedFamilies: string[] }} FontEmbedding
 * @typedef {'svg' | 'raster'} FontEmbeddingKind
 * @typedef {(svgText: string, kind: FontEmbeddingKind) => Promise<FontEmbedding>} EmbedFonts
 */

/**
 * A PNG or PDF size, in SVG pixels, for the Worker's limit check
 * (gbdraw/render/export_limits.py); `dpi` is the PNG's, null for a PDF.
 * @typedef {{ format: 'png' | 'pdf', width: number, height: number, dpi: number | null }} ExportLimitRequest
 * A PNG's width and height in pixels: Python's `png_pixel_size`, the size the
 * command line's PNG of the same figure and DPI has.
 * @typedef {{ width: number, height: number }} PngPixelSize
 * The check's answer for a figure that fits (`ExportLimitReply` in
 * gbdraw/web_support/export_limit.py): `pixels` is the PNG's size, null for a PDF.
 * @typedef {{ ok: true, pixels: PngPixelSize | null }} ExportLimitReply
 * @typedef {(request: ExportLimitRequest) => Promise<ExportLimitReply>} CheckExportLimit
 *   Rejects with the EXPORT_LIMIT error when the figure is over the format's limit.
 */

/**
 * The root `width` and `height` an exported SVG of a Result drawn to a
 * requested size states, in the requested unit (`180mm`; the viewBox stays px):
 * Python's `export_root_size` (gbdraw/render/document.py), which the CLI
 * writes, asked through the export check with the root's px size.
 * @typedef {{ width: string, height: string }} ExportRootSize
 * The check's answer for that SVG (`SvgRootReply` in gbdraw/web_support/export_limit.py).
 * @typedef {{ ok: true, root: ExportRootSize }} SvgRootReply
 * @typedef {(size: { width: number, height: number }) => Promise<ExportRootSize>} SizedSvgRoot
 */

/**
 * The Worker's helpers an export uses. `embedFonts` is null when no font
 * is loaded in this page. `embedInSvg` is Embed font in exported SVG: whether
 * SVG and Interactive SVG carry the fonts; a PNG always draws with them.
 * `checkExportLimit` is the size check of a PNG or PDF, which knows what the
 * exported Result drew; SVG and Interactive SVG have no limit. `sizedSvgRoot`
 * is null for a Result drawn at Auto, whose SVG root stays as it is shown.
 * @typedef {object} ExportPorts
 * @property {EmbedFonts | null} [embedFonts]
 * @property {boolean} [embedInSvg]
 * @property {ResolvePdfFonts} [resolvePdfFonts]
 * @property {CheckExportLimit} [checkExportLimit]
 * @property {SizedSvgRoot | null} [sizedSvgRoot]
 */

/**
 * The size check, before the browser draws anything: a refusal throws, and
 * nothing is downloaded. A missing port fails here, as a missing PDF font port does.
 * @param {CheckExportLimit | undefined} checkExportLimit
 * @param {ExportLimitRequest} request
 * @returns {Promise<ExportLimitReply>}
 */
const requireExportWithinLimit = async (checkExportLimit, request) => {
  if (typeof checkExportLimit !== 'function') throw new Error('The export size check is not available.');
  return checkExportLimit(request);
};

/**
 * What a finished export reports: the file name, and one warning per font the
 * file does not carry as chosen (shown in the export notice).
 * @typedef {{ filename: string, fontWarnings: string[] }} ExportOutcome
 */

/** @type {FontEmbedding} */
const NO_FONT_EMBEDDING = Object.freeze({ block: '', refusedFamilies: [] });

/**
 * Insert the block Python built (`embedded_fonts_block`) as the first child of
 * the root; the export adds nothing of its own.
 * @param {SVGSVGElement} svg
 * @param {string} block
 */
const insertEmbeddedFontBlock = (svg, block) => {
  if (!block) return;
  const parsed = new DOMParser().parseFromString(`<svg xmlns="${SVG_NS}">${block}</svg>`, 'image/svg+xml');
  const fragment = parsed.documentElement.firstElementChild;
  if (!fragment || fragment.namespaceURI !== SVG_NS || fragment.localName !== 'defs') {
    throw diagnosticError('RESULT_INVALID', { field: 'font_family' }, { stage: 'result-admission', operation: 'embedFonts' });
  }
  svg.insertBefore(svg.ownerDocument.importNode(fragment, true), svg.firstChild);
};

/**
 * The fonts an export clone carries, computed from its text before any other
 * change (Interactive SVG adds its script and styles afterwards, as Python
 * computes the embedding of the plain SVG).
 * @param {SVGSVGElement} clone
 * @param {EmbedFonts | null | undefined} embedFonts
 * @param {FontEmbeddingKind} kind
 * @returns {Promise<FontEmbedding>}
 */
const fontEmbeddingOf = async (clone, embedFonts, kind) => (
  embedFonts ? embedFonts(new XMLSerializer().serializeToString(clone), kind) : NO_FONT_EMBEDDING
);

/**
 * @param {ExportPorts} ports
 * @returns {EmbedFonts | null}
 */
const svgFontPort = ({ embedFonts, embedInSvg }) => (embedInSvg && embedFonts ? embedFonts : null);

/**
 * @param {readonly string[]} families
 * @returns {string[]}
 */
const refusedFontWarnings = (families) => families.map((family) => (
  t('Font {family} is not embedded: its license (OS/2 fsType) does not allow embedding. The SVG shows it only where {family} is installed.', { family })
));

/**
 * @typedef {object} ExportSnapshot The captured preview an export action reads.
 * @property {string} name
 * @property {SVGSVGElement | null} [svg]
 * @property {string | null} [svgContent] `getSvgFromString` reads a null content as no SVG.
 * @property {ExportInteractivity} [interactivity] An Interactive SVG's popup mode and edits, and the catalog item it
 *   embeds with its sequences (`embeddedSequences`, which the export reads from Python before it builds the file).
 * @property {number | string} [dpi]
 * @property {boolean} [whiteBackground] the PNG paints the drawing's transparent background white
 */

/**
 * @param {ExportSnapshot} snapshot
 * @param {string} extension
 */
const getDownloadName = (snapshot, extension) => {
  const baseName = snapshot.name;
  if (!extension) return baseName;
  const normalized = baseName.replace(/\.svg$/i, `.${extension}`);
  if (normalized === baseName && !baseName.toLowerCase().endsWith(`.${extension}`)) {
    return `${baseName}.${extension}`;
  }
  return normalized;
};

/**
 * @param {ExportSnapshot} snapshot
 * @returns {SVGSVGElement}
 */
const getCurrentSvgClone = (snapshot) => {
  const clone = snapshot.svg || getSvgFromString(snapshot.svgContent);
  if (clone) {
    if (!clone.getAttribute('xmlns')) {
      clone.setAttribute('xmlns', 'http://www.w3.org/2000/svg');
    }
    if (!clone.getAttribute('xmlns:xlink')) {
      clone.setAttribute('xmlns:xlink', 'http://www.w3.org/1999/xlink');
    }
  }
  if (!clone) {
    throw new Error('No SVG result is available for export.');
  }
  stripTransientPreviewState(clone);
  stripPreviewFeatureSearchClasses(clone);
  stripResultBaseAttributes(clone);
  return clone;
};

/** @param {ExportSnapshot} snapshot */
const getCurrentSvgString = (snapshot) => (
  new XMLSerializer().serializeToString(getCurrentSvgClone(snapshot))
);

/**
 * State the requested size on an export clone's root (`sizedSvgRoot`); the
 * page's root stays in px.
 * @param {SVGSVGElement} clone
 * @param {SizedSvgRoot} sizedSvgRoot
 */
const stateRequestedSize = async (clone, sizedSvgRoot) => {
  const px = getSvgDimensions(clone);
  if (!px) throw new Error('The SVG export has no root size.');
  const root = await sizedSvgRoot(px);
  if (typeof root?.width !== 'string' || typeof root.height !== 'string') {
    throw diagnosticError('RESULT_INVALID', { field: 'figure_size' }, { stage: 'result-admission', operation: 'checkExportLimit' });
  }
  clone.setAttribute('width', root.width);
  clone.setAttribute('height', root.height);
};

/**
 * The exported SVG text with its embedded fonts; a Result drawn to a
 * requested size states it (`sizedSvgRoot`, plain SVG only).
 * @param {ExportSnapshot} snapshot
 * @param {EmbedFonts | null | undefined} embedFonts
 * @param {FontEmbeddingKind} kind
 * @param {SizedSvgRoot | null} [sizedSvgRoot]
 * @returns {Promise<{ svgString: string, refusedFamilies: string[] }>}
 */
const getEmbeddedSvgString = async (snapshot, embedFonts, kind, sizedSvgRoot = null) => {
  const clone = getCurrentSvgClone(snapshot);
  const { block, refusedFamilies } = await fontEmbeddingOf(clone, embedFonts, kind);
  insertEmbeddedFontBlock(clone, block);
  if (sizedSvgRoot) await stateRequestedSize(clone, sizedSvgRoot);
  return { svgString: new XMLSerializer().serializeToString(clone), refusedFamilies };
};

/**
 * @param {ExportSnapshot} snapshot
 * @param {EmbedFonts | null | undefined} embedFonts
 * @returns {Promise<{ svgString: string, refusedFamilies: string[] }>}
 */
const getInteractiveSvgString = async (snapshot, embedFonts) => {
  const { enrichSvgWithStandaloneInteractivity } = await loadStandaloneInteractivity();
  const clone = getCurrentSvgClone(snapshot);
  const { block, refusedFamilies } = await fontEmbeddingOf(clone, embedFonts, 'svg');
  const enriched = enrichSvgWithStandaloneInteractivity(clone, snapshot.interactivity);
  if (!enriched) {
    throw new Error('Interactive SVG export requires the committed feature catalog.');
  }
  insertEmbeddedFontBlock(clone, block);
  return { svgString: new XMLSerializer().serializeToString(clone), refusedFamilies };
};

/**
 * @param {string | null | undefined} svgString
 * @returns {SVGSVGElement | null}
 */
const getSvgFromString = (svgString) => {
  if (!svgString) return null;
  return parseMarkupInertly(svgString).querySelector('svg');
};

/** @param {SVGSVGElement} svg */
const moveGradientsToDefs = (svg) => {
  const defs = ensureSvgDefs(svg);
  svg.querySelectorAll('linearGradient, radialGradient').forEach((gradient) => {
    if (gradient.parentNode !== defs) {
      defs.appendChild(gradient);
    }
  });
};

/**
 * @param {Element} target
 * @param {Element} source
 * @param {string[]} attributes
 */
const copyAttributes = (target, source, attributes) => {
  attributes.forEach((attr) => {
    const value = source.getAttribute(attr);
    if (value !== null) {
      target.setAttribute(attr, value);
    }
  });
};

/**
 * @param {Element} target
 * @param {Element | null | undefined} source
 */
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

/**
 * The element that styles each UTF-16 unit of `textPath`'s text, in `textContent` order: the textPath itself or the
 * tspan run around the unit (an italic run of a region label, OV-2808).
 * @param {Element} textPath
 * @returns {Element[]}
 */
const textPathUnitStyles = (textPath) => {
  /** @type {Element[]} */
  const sources = [];
  const walker = document.createTreeWalker(textPath, NodeFilter.SHOW_TEXT);
  for (let node = walker.nextNode(); node; node = walker.nextNode()) {
    const source = node.parentElement || textPath;
    for (let unit = 0; unit < (node.nodeValue || '').length; unit += 1) sources.push(source);
  }
  return sources;
};

/** @param {SVGSVGElement} svg */
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

    const group = document.createElementNS(SVG_NS, 'g');
    copyAttributes(group, textEl, ['transform', 'opacity', 'display', 'visibility']);
    applyFontAndPaintStyles(group, textPath);
    const unitStyles = textPathUnitStyles(textPath);

    // SVG character positions are UTF-16 indexes; keep surrogate pairs together.
    for (let i = 0; i < charCount && i < textContent.length; i += 1) {
      const code = textContent.charCodeAt(i);
      const char = code >= 0xd800 && code <= 0xdbff ? textContent.slice(i, i + 2) : textContent[i];
      const index = i;
      i += char.length - 1;
      if (char === '\n' || char === '\r') continue;
      let pos;
      try {
        pos = textEl.getStartPositionOfChar(index);
      } catch (error) {
        pos = null;
      }
      if (!pos) continue;
      const rotation = textEl.getRotationOfChar(index) || 0;
      const charText = document.createElementNS(SVG_NS, 'text');
      charText.textContent = char;
      // Without preserve, svg2pdf drops a whitespace-only text run from the text layer.
      if (!char.trim()) charText.setAttributeNS(XML_NS, 'xml:space', 'preserve');
      charText.setAttribute('x', String(pos.x));
      charText.setAttribute('y', String(pos.y));
      charText.setAttribute('text-anchor', 'start');
      charText.setAttribute('dominant-baseline', 'alphabetic');
      const unitStyle = unitStyles[index];
      if (unitStyle && unitStyle !== textPath) applyFontAndPaintStyles(charText, unitStyle);
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

// Reads and writes run in separate phases: interleaving them forces a layout per text.
/** @param {ParentNode} svg */
export const normalizeTextBaselinesForPdf = (svg) => {
  /** @type {{ textEl: SVGTextElement, baseline: string, before: DOMRect }[]} */
  const pending = [];
  svg.querySelectorAll('text').forEach((textEl) => {
    const baseline = textEl.getAttribute('dominant-baseline');
    if (!baseline || baseline === 'alphabetic' || baseline === 'auto') return;
    try {
      pending.push({ textEl, baseline, before: textEl.getBBox() });
    } catch (error) {
      // Text without a measurable box keeps its baseline.
    }
  });

  pending.forEach(({ textEl }) => textEl.setAttribute('dominant-baseline', 'alphabetic'));

  /** @type {{ textEl: SVGTextElement, baseline: string, before: DOMRect, after: DOMRect }[]} */
  const measured = [];
  /** @type {typeof pending} */
  const unmeasurable = [];
  pending.forEach((entry) => {
    try {
      measured.push({ ...entry, after: entry.textEl.getBBox() });
    } catch (error) {
      unmeasurable.push(entry);
    }
  });
  unmeasurable.forEach(({ textEl, baseline }) => textEl.setAttribute('dominant-baseline', baseline));

  measured.forEach(({ textEl, before, after }) => {
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
      textEl.setAttribute('dy', String(nextDy));
    }

    textEl.removeAttribute('dominant-baseline');
  });
};

// svg2pdf draws a `<use>` of a text as a filled copy of its letters, whatever
// the `<use>` paints: a label's halo (a white stroke of its text) would print
// each label twice, in black. The PDF draws labels without halos (LP-D102).
/** @param {ParentNode} svg */
const dropTextUsesForPdf = (svg) => {
  svg.querySelectorAll('use').forEach((use) => {
    const href = use.getAttribute('href') || use.getAttributeNS('http://www.w3.org/1999/xlink', 'href') || '';
    const target = href.startsWith('#') ? svg.querySelector(`[id="${CSS.escape(href.slice(1))}"]`) : null;
    if (target?.tagName.toLowerCase() === 'text') use.remove();
  });
};

/** @param {SVGSVGElement | null} svg */
const prepareSvgForPdf = (svg) => {
  if (!svg) {
    throw new Error('The current SVG could not be prepared for PDF export.');
  }
  // The deep clone of an svg root element is an svg root element.
  const clone = /** @type {SVGSVGElement} */ (svg.cloneNode(true));
  if (!clone.getAttribute('xmlns')) {
    clone.setAttribute('xmlns', SVG_NS);
  }
  if (!clone.getAttribute('xmlns:xlink')) {
    clone.setAttribute('xmlns:xlink', 'http://www.w3.org/1999/xlink');
  }

  const staging = document.createElement('div');
  // The figure's text is English in every UI language and reads left to
  // right: the browser measures it, and picks fallback fonts, as English (as
  // the preview does).
  staging.lang = 'en';
  staging.dir = 'ltr';
  staging.style.position = 'absolute';
  staging.style.left = '-10000px';
  staging.style.top = '-10000px';
  staging.style.visibility = 'hidden';
  staging.appendChild(clone);
  document.body.appendChild(staging);

  try {
    moveGradientsToDefs(clone);
    dropTextUsesForPdf(clone);
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

/**
 * @param {string} svgString
 * @param {string} filename
 */
const downloadSvgString = (svgString, filename) => {
  if (!svgString) {
    throw new Error('The SVG export produced no content.');
  }
  downloadBlob(new Blob([svgString], { type: 'image/svg+xml' }), filename);
};

/**
 * @param {ExportSnapshot} snapshot
 * @param {ExportPorts} [ports]
 * @returns {Promise<ExportOutcome>}
 */
export const downloadSVG = async (snapshot, ports = {}) => {
  const filename = getDownloadName(snapshot, 'svg');
  const { svgString, refusedFamilies } = await getEmbeddedSvgString(snapshot, svgFontPort(ports), 'svg', ports.sizedSvgRoot ?? null);
  downloadSvgString(svgString, filename);
  return { filename, fontWarnings: refusedFontWarnings(refusedFamilies) };
};

/**
 * @param {ExportSnapshot} snapshot
 * @param {ExportPorts} [ports]
 * @returns {Promise<ExportOutcome>}
 */
export const downloadInteractiveSVG = async (snapshot, ports = {}) => {
  const filename = getDownloadName(snapshot, 'interactive.svg');
  const { svgString, refusedFamilies } = await getInteractiveSvgString(snapshot, svgFontPort(ports));
  downloadSvgString(svgString, filename);
  return { filename, fontWarnings: refusedFontWarnings(refusedFamilies) };
};

/**
 * The browser rasterizes the exported SVG as an image, which sees the
 * embedded fonts and the installed fonts, not the fonts of this page; so the
 * SVG carries every custom font it draws with (D-13: whatever Embed font in
 * exported SVG and the font's license; the bytes never leave the page).
 * @param {ExportSnapshot} snapshot
 * @param {ExportPorts} [ports]
 * @returns {Promise<ExportOutcome>}
 */
export const downloadPNG = async (snapshot, { embedFonts, checkExportLimit } = {}) => {
  const filename = getDownloadName(snapshot, 'png');
  const { svgString } = await getEmbeddedSvgString(snapshot, embedFonts, 'raster');
  const svg = getSvgFromString(svgString);
  if (!svg) {
    throw new Error('The current SVG could not be parsed for PNG export.');
  }
  const dims = getSvgDimensions(svg);
  if (!dims) {
    throw new Error('The current SVG has no usable dimensions for PNG export.');
  }
  // The Web sets a number (default 300) or a string; parseInt accepts both.
  const dpi = parseInt(/** @type {string} */ (snapshot.dpi), 10);
  if (!Number.isFinite(dpi) || dpi <= 0) {
    throw new Error('The selected PNG DPI is invalid.');
  }
  const { pixels } = await requireExportWithinLimit(checkExportLimit, { format: 'png', width: dims.width, height: dims.height, dpi });
  if (!pixels) throw new Error('The export size check named no PNG size.');
  // The size the command line's PNG has: Python rounds it as CairoSVG does, where a canvas would truncate.
  const canvas = document.createElement('canvas');
  canvas.width = pixels.width;
  canvas.height = pixels.height;
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
    // The drawing has no background; the PNG keeps it transparent unless asked to paint it white.
    if (snapshot.whiteBackground) {
      ctx.fillStyle = 'white';
      ctx.fillRect(0, 0, canvas.width, canvas.height);
    }
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
    downloadBlob(fixedBlob, filename);
  } finally {
    URL.revokeObjectURL(url);
  }
  return { filename, fontWarnings: [] };
};

/**
 * @param {ExportSnapshot} snapshot
 * @param {ExportPorts} [ports]
 * @returns {Promise<ExportOutcome>}
 */
export const downloadPDF = async (snapshot, { resolvePdfFonts, checkExportLimit } = {}) => {
  const filename = getDownloadName(snapshot, 'pdf');
  const svgString = getCurrentSvgString(snapshot);

  // 1. Create temporary SVG element to read dimensions
  const svg = getSvgFromString(svgString);
  if (!svg) {
    throw new Error('The current SVG could not be parsed for PDF export.');
  }
  // 2. Determine Width/Height (Use attributes or viewBox); the PDF copy keeps them.
  const dims = getSvgDimensions(svg);
  if (!dims) {
    throw new Error('The current SVG has no usable dimensions for PDF export.');
  }
  await requireExportWithinLimit(checkExportLimit, { format: 'pdf', width: dims.width, height: dims.height, dpi: null });
  const prepared = prepareSvgForPdf(svg);
  const { svg: pdfSvg, cleanup } = prepared;

  try {
    if (!pdfSvg.getAttribute('viewBox')) {
      pdfSvg.setAttribute('viewBox', `0 0 ${dims.width} ${dims.height}`);
    }
    const page = { width: dims.width * PT_PER_CSS_PX, height: dims.height * PT_PER_CSS_PX };

    // 3. Load and initialize jsPDF/svg2pdf only for this export format.
    await loadPdfLibraries();
    const { jsPDF } = window.jspdf;
    const doc = new jsPDF({
      orientation: page.width > page.height ? 'l' : 'p',
      unit: 'pt',
      format: [page.width, page.height]
    });

    const fontWarnings = await preparePdfFonts(doc, pdfSvg, resolvePdfFonts);

    // 4. Convert SVG to PDF
    await doc.svg(pdfSvg, {
      x: 0,
      y: 0,
      width: page.width,
      height: page.height
    });
    // 5. Save File
    doc.save(filename);
    return { filename, fontWarnings };
  } finally {
    cleanup();
  }
};
