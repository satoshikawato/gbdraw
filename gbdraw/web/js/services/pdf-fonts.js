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
 * Reads the base64 bytes of one packaged font file.
 * @typedef {(filename: string) => Promise<string> | string} PdfFontLoader
 */

/** @type {Record<PdfFontStyle, string>} */
const styles = { normal: 'Regular', bold: 'Bold', italic: 'Italic', bolditalic: 'BoldItalic' };
/** @type {Map<string, Promise<string>>} */
const payloads = new Map();

/**
 * @param {string} filename
 * @param {PdfFontLoader | undefined} loadFont
 * @returns {Promise<string>}
 */
const fontPayload = (filename, loadFont) => {
  if (!payloads.has(filename)) {
    // PDF export always passes the loader; a missing one rejects here as a TypeError.
    const pending = Promise.resolve().then(() => /** @type {PdfFontLoader} */ (loadFont)(filename)).catch((error) => {
      payloads.delete(filename);
      throw error;
    });
    payloads.set(filename, pending);
  }
  // The entry was set above when it was missing.
  return /** @type {Promise<string>} */ (payloads.get(filename));
};

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

/** @param {string} family */
const fontFamily = (family) => {
  const first = family.split(',')[0].replaceAll(/['"]/g, '').trim().toLowerCase();
  if (/mono|courier/.test(first)) return 'LiberationMono';
  if (/serif|times/.test(first) && !first.includes('sans')) return 'LiberationSerif';
  return 'LiberationSans';
};

// jsPDF's built-in PDF fonts use WinAnsi. Register Unicode fonts explicitly,
// from the existing Python package before svg2pdf resolves any text metrics.
/**
 * @param {PdfFontDocument} doc
 * @param {Element} svg
 * @param {PdfFontLoader | undefined} loadFont The loader; PDF export passes it, and a missing one fails the first font read.
 */
export const preparePdfFonts = async (doc, svg, loadFont) => {
  /** @type {Map<string, PdfFont>} */
  const fonts = new Map();
  /**
   * @param {string} family
   * @param {PdfFontStyle} style
   */
  const register = async (family, style) => {
    const key = `${family}:${style}`;
    if (!fonts.has(key)) {
      const filename = `${family}-${styles[style]}.ttf`;
      doc.addFileToVFS(filename, await fontPayload(filename, loadFont));
      doc.addFont(filename, family, style);
      fonts.set(key, doc.getFont(family, style));
    }
    // The key was set above when it was missing.
    return /** @type {PdfFont} */ (fonts.get(key));
  };
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
    return { node, family: fontFamily(style.fontFamily), style: fontStyle(style) };
  });
  for (const { node, family, style } of runs) {
    const primary = await register(family, style);
    for (const character of node.textContent) {
      // Iterating a string yields characters that each hold a code point.
      if (!/\s/.test(character) && !primary.metadata.characterToGlyph(/** @type {number} */ (character.codePointAt(0)))) {
        throw new Error(`PDF fonts do not contain U+${/** @type {number} */ (character.codePointAt(0)).toString(16).toUpperCase()}. Use SVG to retain this text.`);
      }
    }
    // The parent controls anchor and metrics even when its content has tspans.
    node.parentElement.setAttribute('font-family', family);
  }
};
