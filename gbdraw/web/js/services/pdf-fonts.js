// @ts-check
import { t } from '../utils/i18n.js';
/** @typedef {'normal' | 'bold' | 'italic' | 'bolditalic'} PdfFontStyle */

/**
 * The metrics of a registered font, as far as the glyph check reads them.
 * @typedef {{ metadata: { characterToGlyph: (codePoint: number) => unknown } }} PdfFont
 */

/**
 * The part of a jsPDF document that registers fonts.
 * @typedef {object} PdfFontDocument
 * @property {(filename: string, base64: string) => void} addFileToVFS
 * @property {(filename: string, family: string, style: PdfFontStyle) => void} addFont
 * @property {(family: string, style: PdfFontStyle) => PdfFont} getFont
 * @property {(text: string, options: { font: PdfFont, fontSize: number }) => number} getStringUnitWidth
 */

/**
 * A text run of the PDF: the CSS font-family list its text was drawn with,
 * its jsPDF style slot, and the characters it draws (a custom face is sent
 * as a subset of them, S3).
 * @typedef {{ fontFamily: string, bold: boolean, italic: boolean, text: string }} PdfFontRun
 */

/**
 * A custom font the PDF draws with a Liberation face instead: its license
 * forbids embedding, or its file cannot become TrueType outlines.
 * @typedef {{ family: string, substitute: string, reason: 'license' | 'font-file' }} PdfFontSubstitution
 */

/**
 * Python's `pdf_fonts` (gbdraw/render/font_embedding.py) in the Worker: per
 * run, the jsPDF family and the TrueType file of the face its text was
 * measured with (a custom face: a subset of the characters drawn with it, and
 * a file name that names them); the base64 bytes of each file not in `known`;
 * and the custom fonts drawn with a substitute.
 * @typedef {object} PdfFontsReply
 * @property {{ family: string, file: string }[]} faces
 * @property {Record<string, string>} fonts
 * @property {PdfFontSubstitution[]} substitutions
 * @typedef {(runs: PdfFontRun[], known: string[]) => Promise<PdfFontsReply>} ResolvePdfFonts
 */

/**
 * The base64 TrueType bytes of each font file the PDF export received, by file
 * name. A name is a face's measurement name, which identifies its bytes
 * (a custom face's name holds the hash of its file).
 * @type {Map<string, string>}
 */
const payloads = new Map();

/**
 * @param {CSSStyleDeclaration} style
 * @returns {PdfFontStyle}
 */
const fontStyle = (style) => {
  const bold = style.fontWeight === 'bold' || Number(style.fontWeight) >= 600;
  const italic = ['italic', 'oblique'].includes(style.fontStyle);
  // The two parts spell only the style names of `PdfFontStyle`.
  return /** @type {PdfFontStyle} */ (`${bold ? 'bold' : ''}${italic ? 'italic' : ''}` || 'normal');
};

// The attributes that place part of a text on its own: such a text keeps its
// anchor for svg2pdf to resolve.
const RUN_LAYOUT_ATTRIBUTES = ['x', 'dx', 'rotate', 'textLength', 'letter-spacing', 'word-spacing'];
/** @type {Readonly<Record<string, number | undefined>>} */
const ANCHOR_SHIFT = Object.freeze({ middle: 0.5, end: 1 });

/**
 * Draw each middle- or end-anchored `<text>` from its start, at the x its
 * anchor gives with the advances of the faces the PDF draws it with.
 * svg2pdf would measure the anchor with the browser, which knows the jsPDF
 * family under no name (or with other advances), and shift the text (OV-1303).
 * A text that places a part on its own keeps its anchor.
 * @param {PdfFontDocument} doc
 * @param {Element} svg
 * @param {Map<Node, { font: PdfFont, size: number }>} drawnWith The face and font size (px) of each text node.
 */
const anchorTextsAtStart = (doc, svg, drawnWith) => {
  svg.querySelectorAll('text').forEach((text) => {
    const shift = ANCHOR_SHIFT[getComputedStyle(text).textAnchor];
    const x = Number(text.getAttribute('x') ?? 0);
    const parts = [...text.querySelectorAll('*')];
    if (shift === undefined || !Number.isFinite(x) || parts.some((part) => RUN_LAYOUT_ATTRIBUTES.some((name) => part.hasAttribute(name)))) return;
    let width = 0;
    let afterSpace = true;
    /** @type {{ text: string, font: PdfFont, size: number }[]} */
    const runs = [];
    const walker = document.createTreeWalker(text, NodeFilter.SHOW_TEXT);
    while (walker.nextNode()) {
      const drawn = drawnWith.get(walker.currentNode);
      // Spaces collapse across runs, as in the SVG.
      let characters = (walker.currentNode.textContent ?? '').replace(/\s+/g, ' ');
      if (afterSpace) characters = characters.replace(/^ /, '');
      if (!drawn || !characters) continue;
      afterSpace = characters.endsWith(' ');
      runs.push({ text: characters, ...drawn });
    }
    if (runs.length === 0) return;
    // A trailing space does not count toward the anchored width.
    runs[runs.length - 1].text = runs[runs.length - 1].text.replace(/ $/, '');
    for (const run of runs) width += doc.getStringUnitWidth(run.text, { font: run.font, fontSize: run.size }) * run.size;
    text.setAttribute('x', String(x - shift * width));
    text.setAttribute('text-anchor', 'start');
  });
};

/**
 * @param {PdfFontSubstitution} substitution
 * @returns {string}
 */
const substitutionWarning = ({ family, substitute, reason }) => (reason === 'license'
  ? t('Font {family} is not embedded in the PDF: its license (OS/2 fsType) does not allow embedding. Its text is drawn in {substitute}.', { family, substitute })
  : t('Font {family} could not be prepared for the PDF. Its text is drawn in {substitute}.', { family, substitute }));

// jsPDF's built-in PDF fonts use WinAnsi. Register the Unicode TrueType face
// each text was measured with (a Liberation font, or a custom font loaded in
// this page) before svg2pdf resolves any text metrics.
/**
 * @param {PdfFontDocument} doc
 * @param {Element} svg
 * @param {ResolvePdfFonts | undefined} resolveFonts The Worker's `pdfFonts` helper; PDF export passes it, and a missing one fails here.
 * @returns {Promise<string[]>} One warning per custom font drawn with a substitute.
 */
export const preparePdfFonts = async (doc, svg, resolveFonts) => {
  /** @type {Array<Node & { parentElement: HTMLElement, textContent: string }>} */
  const nodes = [];
  const walker = document.createTreeWalker(svg, NodeFilter.SHOW_TEXT);
  while (walker.nextNode()) {
    const node = walker.currentNode;
    // The condition checks the parent element and the text.
    if (node.parentElement?.closest('text') && node.textContent) nodes.push(/** @type {(typeof nodes)[number]} */ (node));
  }
  // Capture computed styles before replacing inherited font families.
  const runs = nodes.map((node) => {
    const style = getComputedStyle(node.parentElement);
    return { node, fontFamily: style.fontFamily, style: fontStyle(style), size: parseFloat(style.fontSize) };
  });
  /** @type {Map<string, number>} */
  const runIndex = new Map();
  /** @type {PdfFontRun[]} */
  const requested = [];
  /** @type {Set<string>[]} */
  const characters = [];
  for (const run of runs) {
    const key = `${run.style}\n${run.fontFamily}`;
    if (!runIndex.has(key)) {
      runIndex.set(key, requested.length);
      requested.push({ fontFamily: run.fontFamily, bold: run.style.startsWith('bold'), italic: run.style.endsWith('italic'), text: '' });
      characters.push(new Set());
    }
    const drawn = characters[/** @type {number} */ (runIndex.get(key))];
    for (const character of run.node.textContent) drawn.add(character);
  }
  requested.forEach((run, index) => { run.text = [...characters[index]].sort().join(''); });
  if (requested.length === 0) return [];
  // PDF export always passes the helper; a missing one rejects here as a TypeError.
  const reply = await /** @type {ResolvePdfFonts} */ (resolveFonts)(requested, [...payloads.keys()]);
  Object.entries(reply.fonts).forEach(([file, base64]) => payloads.set(file, base64));

  /** @type {Set<string>} */
  const files = new Set();
  /** @type {Map<string, PdfFont>} */
  const fonts = new Map();
  /** @type {Map<Node, { font: PdfFont, size: number }>} */
  const drawnWith = new Map();
  /**
   * @param {{ family: string, file: string }} face
   * @param {PdfFontStyle} style
   */
  const register = (face, style) => {
    const key = `${face.family}:${style}`;
    if (!fonts.has(key)) {
      if (!files.has(face.file)) {
        // Python sent the bytes of every file this page did not hold.
        doc.addFileToVFS(face.file, /** @type {string} */ (payloads.get(face.file)));
        files.add(face.file);
      }
      doc.addFont(face.file, face.family, style);
      fonts.set(key, doc.getFont(face.family, style));
    }
    // The key was set above when it was missing.
    return /** @type {PdfFont} */ (fonts.get(key));
  };
  for (const { node, fontFamily, style, size } of runs) {
    const face = reply.faces[/** @type {number} */ (runIndex.get(`${style}\n${fontFamily}`))];
    const primary = register(face, style);
    drawnWith.set(node, { font: primary, size });
    for (const character of node.textContent) {
      // Iterating a string yields characters that each hold a code point.
      if (!/\s/.test(character) && !primary.metadata.characterToGlyph(/** @type {number} */ (character.codePointAt(0)))) {
        throw new Error(`PDF fonts do not contain U+${/** @type {number} */ (character.codePointAt(0)).toString(16).toUpperCase()}. Use SVG to retain this text.`);
      }
    }
    // The parent controls anchor and metrics even when its content has tspans.
    node.parentElement.setAttribute('font-family', face.family);
  }
  anchorTextsAtStart(doc, svg, drawnWith);
  return reply.substitutions.map(substitutionWarning);
};
