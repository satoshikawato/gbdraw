// @ts-check
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
 */

/**
 * A text run of the PDF: the CSS font-family list its text was drawn with,
 * and its jsPDF style slot.
 * @typedef {{ fontFamily: string, bold: boolean, italic: boolean }} PdfFontRun
 */

/**
 * A custom font the PDF draws with a Liberation face instead: its license
 * forbids embedding, or its file cannot become TrueType outlines.
 * @typedef {{ family: string, substitute: string, reason: 'license' | 'font-file' }} PdfFontSubstitution
 */

/**
 * Python's `pdf_fonts` (gbdraw/render/font_embedding.py) in the Worker: per
 * run, the jsPDF family and the TrueType file of the face its text was
 * measured with; the base64 bytes of each file not in `known`; and the custom
 * fonts drawn with a substitute.
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

/**
 * @param {PdfFontSubstitution} substitution
 * @returns {string}
 */
const substitutionWarning = ({ family, substitute, reason }) => (reason === 'license'
  ? `Font ${family} is not embedded in the PDF: its license (OS/2 fsType) does not allow embedding. Its text is drawn in ${substitute}.`
  : `Font ${family} could not be prepared for the PDF. Its text is drawn in ${substitute}.`);

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
    return { node, fontFamily: style.fontFamily, style: fontStyle(style) };
  });
  /** @type {Map<string, number>} */
  const runIndex = new Map();
  /** @type {PdfFontRun[]} */
  const requested = [];
  for (const run of runs) {
    const key = `${run.style}\n${run.fontFamily}`;
    if (!runIndex.has(key)) {
      runIndex.set(key, requested.length);
      requested.push({ fontFamily: run.fontFamily, bold: run.style.startsWith('bold'), italic: run.style.endsWith('italic') });
    }
  }
  if (requested.length === 0) return [];
  // PDF export always passes the helper; a missing one rejects here as a TypeError.
  const reply = await /** @type {ResolvePdfFonts} */ (resolveFonts)(requested, [...payloads.keys()]);
  Object.entries(reply.fonts).forEach(([file, base64]) => payloads.set(file, base64));

  /** @type {Set<string>} */
  const files = new Set();
  /** @type {Map<string, PdfFont>} */
  const fonts = new Map();
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
  for (const { node, fontFamily, style } of runs) {
    const face = reply.faces[/** @type {number} */ (runIndex.get(`${style}\n${fontFamily}`))];
    const primary = register(face, style);
    for (const character of node.textContent) {
      // Iterating a string yields characters that each hold a code point.
      if (!/\s/.test(character) && !primary.metadata.characterToGlyph(/** @type {number} */ (character.codePointAt(0)))) {
        throw new Error(`PDF fonts do not contain U+${/** @type {number} */ (character.codePointAt(0)).toString(16).toUpperCase()}. Use SVG to retain this text.`);
      }
    }
    // The parent controls anchor and metrics even when its content has tspans.
    node.parentElement.setAttribute('font-family', face.family);
  }
  return reply.substitutions.map(substitutionWarning);
};
