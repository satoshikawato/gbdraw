import { bytesToBase64 } from './byte-utils.js';

const FONT_ROOT = new URL('../../vendor/fonts/', import.meta.url);
const payloads = new Map();
const styles = { normal: 'Regular', bold: 'Bold', italic: 'Italic', bolditalic: 'BoldItalic' };

const fontPayload = (filename) => {
  if (!payloads.has(filename)) {
    const pending = fetch(new URL(filename, FONT_ROOT)).then(async (response) => {
      if (!response.ok) throw new Error(`Failed to load PDF font (${response.status}).`);
      return bytesToBase64(new Uint8Array(await response.arrayBuffer()));
    }).catch((error) => {
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
// including glyph fallback, before svg2pdf resolves any text metrics.
export const preparePdfFonts = async (doc, svg) => {
  const fonts = new Map();
  const register = async (family, style) => {
    const key = `${family}:${style}`;
    if (!fonts.has(key)) {
      const filename = family === 'NotoSansJP'
        ? `noto-sans-jp/noto-sans-jp-japanese-${style.includes('bold') ? 700 : 400}-normal.ttf`
        : `${family}-${styles[style]}.ttf`;
      doc.addFileToVFS(filename, await fontPayload(filename));
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
    const fragments = [];
    for (const character of node.textContent) {
      let resolvedFamily = family;
      if (!/\s/.test(character) && !primary.metadata.characterToGlyph(character.codePointAt(0))) {
        const fallback = await register('NotoSansJP', style);
        if (!fallback.metadata.characterToGlyph(character.codePointAt(0))) {
          throw new Error(`PDF fonts do not contain U+${character.codePointAt(0).toString(16).toUpperCase()}. Use SVG to retain this text.`);
        }
        resolvedFamily = 'NotoSansJP';
      }
      const previous = fragments.at(-1);
      if (previous?.family === resolvedFamily) previous.text += character;
      else fragments.push({ family: resolvedFamily, text: character });
    }
    const spans = fragments.map((fragment) => {
      const span = document.createElementNS('http://www.w3.org/2000/svg', 'tspan');
      span.setAttribute('font-family', fragment.family);
      span.textContent = fragment.text;
      return span;
    });
    // The parent controls anchor and metrics even when its content has tspans.
    node.parentElement.setAttribute('font-family', family);
    node.replaceWith(...spans);
  }
};
