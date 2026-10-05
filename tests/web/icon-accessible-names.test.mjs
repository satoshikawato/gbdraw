// R12 / OV-30: a Phosphor icon is decorative (its glyph is a private-use
// character from a CSS ::before), and a button whose only content is an icon
// carries an author-provided aria-label instead of taking the glyph as its name.
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { join } from 'node:path';
import test from 'node:test';

const html = readFileSync(join(process.cwd(), 'gbdraw', 'web', 'index.html'), 'utf8');
const lineOf = (offset) => html.slice(0, offset).split('\n').length;

// The end of the opening tag, skipping '>' inside quoted attribute values.
const openingTagEnd = (from) => {
  let quote = null;
  for (let index = from; index < html.length; index += 1) {
    const char = html[index];
    if (quote) {
      if (char === quote) quote = null;
    } else if (char === '"' || char === "'") {
      quote = char;
    } else if (char === '>') {
      return index;
    }
  }
  throw new Error(`Unterminated tag at line ${lineOf(from)}`);
};

test('every Phosphor icon in index.html is aria-hidden', () => {
  const missing = [];
  for (const match of html.matchAll(/<i\b/g)) {
    const tag = html.slice(match.index, openingTagEnd(match.index) + 1);
    if (!/\bph\b|\bph-[a-z]/.test(tag)) continue;
    if (!/\baria-hidden="true"/.test(tag)) missing.push(`line ${lineOf(match.index)}: ${tag.slice(0, 100)}`);
  }
  assert.deepEqual(missing, []);
});

// The characters of an HTML fragment outside its tags (a scan, not a regex
// strip, so no tag remnant can survive).
const textOutsideTags = (fragment) => {
  let text = '';
  let inTag = false;
  for (const ch of fragment) {
    if (ch === '<') inTag = true;
    else if (ch === '>') inTag = false;
    else if (!inTag) text += ch;
  }
  return text;
};

test('every icon-only button in index.html has an aria-label', () => {
  const unnamed = [];
  for (const match of html.matchAll(/<button\b/g)) {
    const tagEnd = openingTagEnd(match.index);
    const tag = html.slice(match.index, tagEnd + 1);
    const inner = html.slice(tagEnd + 1, html.indexOf('</button>', tagEnd));
    if (!/<i\b/.test(inner) || /\sv-text=/.test(inner)) continue;
    const text = textOutsideTags(inner).replace(/\{\{[\s\S]*?\}\}/g, 'x').trim();
    if (text) continue;
    if (!/\saria-label(ledby)?=|\s:aria-label=/.test(tag)) unnamed.push(`line ${lineOf(match.index)}: ${tag.slice(0, 100)}`);
  }
  assert.deepEqual(unnamed, []);
});
