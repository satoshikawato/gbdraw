// Read text operators and ToUnicode maps from the uncompressed PDFs under test.
// This verifies the downloaded bytes, independently of the browser's SVG input.
const readPdfText = (bytes) => {
  const source = bytes.toString('latin1');
  const objects = new Map([...source.matchAll(/(\d+) 0 obj\s*([\s\S]*?)endobj/g)].map((match) => [match[1], match[2]]));
  const fonts = new Map();
  for (const [, font, id] of source.matchAll(/\/(F\d+)\s+(\d+) 0 R/g)) {
    const cmapId = objects.get(id)?.match(/\/ToUnicode (\d+) 0 R/)?.[1];
    const cmap = objects.get(cmapId)?.match(/beginbfchar([\s\S]*?)endbfchar/)?.[1] || '';
    fonts.set(font, new Map([...cmap.matchAll(/<([0-9a-f]+)>\s*<([0-9a-f]+)>/gi)]
      .map(([, glyph, unicode]) => [glyph.toLowerCase(), String.fromCharCode(...unicode.match(/.{4}/g).map((code) => parseInt(code, 16)))])));
  }
  let text = '';
  for (const body of objects.values()) {
    if (!body.includes('stream') || !body.includes('BT')) continue;
    let font = '';
    for (const [, nextFont, hex, literal] of body.matchAll(/\/(F\d+)\s+[\d.]+\s+Tf|<([0-9a-f]+)>\s*Tj|\(([^()]*)\)\s*Tj/gi)) {
      if (nextFont) font = nextFont;
      else if (hex) text += hex.match(/.{4}/g).map((code) => fonts.get(font)?.get(code.toLowerCase()) || '\uFFFD').join('');
      else text += literal;
    }
  }
  return text;
};
module.exports = { readPdfText };
