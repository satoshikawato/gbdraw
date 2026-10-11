// @ts-check
// Record and title text may carry <i>…</i> markup, which the diagram draws in
// italics (UX-C06, UX-C15). A field shows a preview of that text so the raw
// tags are not the only reading; any other markup stays literal text.

/** @typedef {{ text: string, italic: boolean }} MarkupSegment */

const ITALIC_TAG = /(<\/?i>)/i;

/**
 * The display segments of text with italic markup, or [] when it has none.
 * @param {unknown} value
 * @returns {MarkupSegment[]}
 */
export const italicSegments = (value) => {
  const text = String(value ?? '');
  if (!ITALIC_TAG.test(text)) return [];
  /** @type {MarkupSegment[]} */
  const segments = [];
  let italic = false;
  for (const part of text.split(ITALIC_TAG)) {
    if (/^<i>$/i.test(part)) italic = true;
    else if (/^<\/i>$/i.test(part)) italic = false;
    else if (part) segments.push({ text: part, italic });
  }
  return segments;
};
