const styles = { normal: 'Regular', bold: 'Bold', italic: 'Italic', bolditalic: 'BoldItalic' };
const payloads = new Map();

const fontPayload = (filename, loadFont) => {
  if (!payloads.has(filename)) {
    const pending = Promise.resolve().then(() => loadFont(filename)).catch((error) => {
      payloads.delete(filename);
      throw error;
    });
    payloads.set(filename, pending);
  }
  return payloads.get(filename);
};

const fontStyle = (style) => {
  const bold = style.fontWeight === 'bold' || Number(style.fontWeight) >= 600;
  const italic = ['italic', 'oblique'].includes(style.fontStyle);
  return `${bold ? 'bold' : ''}${italic ? 'italic' : ''}` || 'normal';
};

const fontFamily = (family) => {
  const first = family.split(',')[0].replaceAll(/['"]/g, '').trim().toLowerCase();
  if (/mono|courier/.test(first)) return 'LiberationMono';
  if (/serif|times/.test(first) && !first.includes('sans')) return 'LiberationSerif';
  return 'LiberationSans';
};

// jsPDF's built-in PDF fonts use WinAnsi. Register Unicode fonts explicitly,
// from the existing Python package before svg2pdf resolves any text metrics.
export const preparePdfFonts = async (doc, svg, loadFont) => {
  const fonts = new Map();
  const register = async (family, style) => {
    const key = `${family}:${style}`;
    if (!fonts.has(key)) {
      const filename = `${family}-${styles[style]}.ttf`;
      doc.addFileToVFS(filename, await fontPayload(filename, loadFont));
      doc.addFont(filename, family, style);
      fonts.set(key, doc.getFont(family, style));
    }
    return fonts.get(key);
  };
  const nodes = [];
  const walker = document.createTreeWalker(svg, NodeFilter.SHOW_TEXT);
  while (walker.nextNode()) {
    const node = walker.currentNode;
    if (node.parentElement?.closest('text') && node.textContent) nodes.push(node);
  }
  // Capture computed styles before replacing inherited font families.
  const runs = nodes.map((node) => {
    const style = getComputedStyle(node.parentElement);
    return { node, family: fontFamily(style.fontFamily), style: fontStyle(style) };
  });
  for (const { node, family, style } of runs) {
    const primary = await register(family, style);
    for (const character of node.textContent) {
      if (!/\s/.test(character) && !primary.metadata.characterToGlyph(character.codePointAt(0))) {
        throw new Error(`PDF fonts do not contain U+${character.codePointAt(0).toString(16).toUpperCase()}. Use SVG to retain this text.`);
      }
    }
    // The parent controls anchor and metrics even when its content has tspans.
    node.parentElement.setAttribute('font-family', family);
  }
};
