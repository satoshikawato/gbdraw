// The settings-section catalog (services/settings-sections.js) is the one
// owner of a section's key, title, icon, modes, and rail entry: index.html
// names a section only by its key, once on the element and once in its
// heading, in catalog order, and every top-level child of the settings scroll
// area is a section that one rail entry docks.
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { join } from 'node:path';
import test from 'node:test';

import {
  SETTINGS_SECTIONS, settingsEntriesOf, settingsEntryIn, settingsEntryOf
} from '../../gbdraw/web/js/services/settings-sections.js';

const html = readFileSync(join(process.cwd(), 'gbdraw', 'web', 'index.html'), 'utf8');
const keys = SETTINGS_SECTIONS.map((section) => section.key);
const MODES = ['circular', 'linear'];

test('the stacked row lists the sections of each mode in card order', () => {
  assert.deepEqual(settingsEntriesOf('circular', 'stacked').map((section) => section.title), [
    'Input Genomes', 'Region Annotations', 'Pairwise Comparisons', 'Layout', 'Font', 'Titles & Record Labels', 'Legend',
    'Colors', 'Features', 'Labels', 'Axis & Scale', 'Dinucleotide content/skew'
  ]);
  assert.deepEqual(settingsEntriesOf('linear', 'stacked').map((section) => section.title), [
    'Input Genomes', 'Comparison', 'Advanced comparison and layout', 'Region Annotations', 'Layout',
    'Font', 'Titles & Record Labels', 'Legend', 'Colors', 'Features', 'Labels', 'Axis & Scale', 'Dinucleotide content/skew'
  ]);
  assert.equal(new Set(keys).size, keys.length, 'section keys are unique');
  for (const mode of MODES) {
    for (const section of settingsEntriesOf(mode, 'rail')) {
      assert.match(section.icon, /^ph-[a-z-]+$/, section.key);
      assert.ok(section.title && section.short, section.key);
    }
  }
  for (const section of SETTINGS_SECTIONS) assert.ok(section.modes.length > 0, section.key);
});

test('the rail lists the same entries under short names, with About & Citation last', () => {
  assert.deepEqual(settingsEntriesOf('circular', 'rail').map((section) => section.short), [
    'Input', 'Regions', 'Compare', 'Layout', 'Font', 'Titles', 'Legend', 'Colors', 'Features', 'Labels', 'Axis',
    'GC / skew', 'About'
  ]);
  assert.deepEqual(settingsEntriesOf('linear', 'rail').map((section) => section.short), [
    'Input', 'Compare', 'Advanced', 'Regions', 'Layout', 'Font', 'Titles', 'Legend', 'Colors', 'Features', 'Labels',
    'Axis', 'GC / skew', 'About'
  ]);
  for (const mode of MODES) {
    const rail = settingsEntriesOf(mode, 'rail');
    // The row is the rail without the entries the rail alone has, in the same order.
    assert.deepEqual(settingsEntriesOf(mode, 'stacked'), rail.filter((section) => !section.railOnly), mode);
    assert.deepEqual(rail.filter((section) => section.railOnly).map((section) => section.title), ['About & Citation'], mode);
    assert.equal(rail.at(-1).key, 'about', mode);
  }
});

test('every section is docked by exactly one rail entry of each mode that shows it', () => {
  for (const mode of MODES) {
    const entries = settingsEntriesOf(mode, 'rail').map((section) => section.key);
    for (const section of SETTINGS_SECTIONS.filter((candidate) => candidate.modes.includes(mode))) {
      assert.ok(entries.includes(settingsEntryOf(section.key)), `${mode} ${section.key} -> ${settingsEntryOf(section.key)}`);
    }
  }
  // An entry docks its own section; Basic and the preserved-settings notice go with Input.
  for (const section of settingsEntriesOf('linear', 'rail')) assert.equal(settingsEntryOf(section.key), section.key);
  assert.equal(settingsEntryOf('basic'), 'input-genomes');
  assert.equal(settingsEntryOf('preserved-settings'), 'input-genomes');
});

// A docked entry follows a mode change to the entry of the same section: its
// own, or the one with the same title (Axis & Scale is one card per mode).
// Comparison and Pairwise Comparisons share a short name, not a title.
test('the entry of a section in a mode is its own or the one of the same title', () => {
  for (const mode of MODES) {
    const entries = settingsEntriesOf(mode, 'rail');
    for (const section of SETTINGS_SECTIONS) {
      const entry = settingsEntryOf(section.key);
      const title = SETTINGS_SECTIONS.find((candidate) => candidate.key === entry).title;
      const expected = entries.find((candidate) => candidate.key === entry) ?? entries.find((candidate) => candidate.title === title);
      assert.equal(settingsEntryIn(mode, section.key), expected?.key, `${mode} ${section.key}`);
    }
  }
  assert.equal(settingsEntryIn('linear', 'circular-axis'), 'linear-axis');
  assert.equal(settingsEntryIn('circular', 'linear-axis'), 'circular-axis');
  assert.equal(settingsEntryIn('circular', 'comparison'), undefined);
  assert.equal(settingsEntryIn('linear', 'pairwise-comparisons'), undefined);
  assert.equal(settingsEntryIn('circular', 'advanced-comparison'), undefined);
  assert.equal(settingsEntryIn('linear', 'basic'), 'input-genomes');
});

// The elements of a template slice with their nesting depth, and the depth
// at its end. A quoted attribute value may hold `>`.
const VOID = new Set(['area', 'base', 'br', 'col', 'embed', 'hr', 'img', 'input', 'link', 'meta', 'source', 'track', 'wbr']);
const elementsOf = (markup) => {
  const found = [];
  let depth = 0;
  const tag = /<!--[\s\S]*?-->|<(\/?)([a-zA-Z][\w-]*)((?:\s+(?:[^"'>]|"[^"]*"|'[^']*')*)?)>/g;
  for (const match of markup.matchAll(tag)) {
    const [, closing, name, attributes = ''] = match;
    if (!name) continue;
    if (closing) { depth -= 1; continue; }
    found.push({ name, attributes, depth });
    if (!VOID.has(name) && !attributes.trimEnd().endsWith('/')) depth += 1;
  }
  return { elements: found, depth };
};

test('index.html names each section by its key, once on the element and once in the heading, in catalog order', () => {
  const marked = Array.from(html.matchAll(/\sv-settings-section="'([^']*)'"/g), (match) => match[1]);
  const headings = Array.from(html.matchAll(/<settings-section-title section="([^"]*)">/g), (match) => match[1]);
  assert.deepEqual(marked, keys, 'v-settings-section values');
  assert.equal(html.includes('data-settings-section='), false, 'the directive sets data-settings-section');
  // A section with a title shows it through the catalog.
  const titled = SETTINGS_SECTIONS.filter((section) => 'title' in section).map((section) => section.key);
  assert.deepEqual(headings, titled, 'settings-section-title headings');
  // Each heading sits inside its own section: the key comes first, and the
  // next titled section starts after the heading.
  const offsets = titled.map((key) => [
    html.indexOf(`v-settings-section="'${key}'"`), html.indexOf(`<settings-section-title section="${key}">`)
  ]);
  for (const [order, [card, heading]] of offsets.entries()) {
    assert.ok(card < heading, `${titled[order]}: heading after its key`);
    if (order > 0) assert.ok(offsets[order - 1][1] < card, `${titled[order]}: after the previous heading`);
  }
});

// A section that the popup beside the Result can hold stands inside a
// `<teleport>`, which renders it in place while it is in the pane.
const paneMarkup = () => {
  const start = html.indexOf('<div id="settings-sections"');
  return html.slice(start, html.indexOf('<div class="generate-bar', start)).replace(/<\/?teleport\b[^>]*>/g, '');
};

test('every top-level child of the settings scroll area is a section, so a rail entry docks it', () => {
  const { elements: area, depth } = elementsOf(paneMarkup());
  assert.equal(area[0].depth, 0);
  // The scroll area closes before the Generate bar: the bar shows in every display.
  assert.equal(depth, 0, 'the Generate bar is outside the scroll area');
  const children = area.filter((element) => element.depth === 1);
  assert.ok(children.length >= 15, `${children.length} children found`);
  const unkeyed = children.filter((element) => !/\sv-settings-section="'[^']*'"/.test(element.attributes))
    .map((element) => `<${element.name}${element.attributes.replace(/\s+/g, ' ').slice(0, 80)}>`);
  assert.deepEqual(unkeyed, []);
  const nested = area.filter((element) => element.depth > 1 && /\sv-settings-section=/.test(element.attributes))
    .map((element) => element.attributes.match(/v-settings-section="'([^']*)'"/)[1]);
  assert.deepEqual(nested, ['comparison'], 'Comparison is the one section inside another');
});

// With a nested section docked, the rail folds away what stands beside it:
// the other children of the card, and the blocks beside the section inside the
// card's block. The style rule names their tags, so a child with another tag
// would stay on screen.
test('what stands beside a nested section has a tag the rail folds away', () => {
  const rule = html.match(/\[data-settings-in-dock="false"\]:not\(:has[^{]*\{ display: none; \}/)[0];
  const tagsAfter = (path) => new Set(rule.split(',').map((selector) => selector.trim().split(`[data-settings-in-dock="false"] > ${path}`)[1])
    .filter((rest) => rest !== undefined && /^[a-z]+(?![a-z >])/.test(rest)).map((rest) => rest.match(/^[a-z]+/)[0]));
  const { elements: area } = elementsOf(paneMarkup());
  const nested = area.findIndex((element) => element.depth > 1 && /\sv-settings-section=/.test(element.attributes));
  assert.equal(area[nested].depth, 3, 'the nested section is a grandchild of its card');
  const parentOf = (index) => area.findLastIndex((element, order) => order < index && element.depth === area[index].depth - 1);
  const childrenOf = (index) => {
    const end = area.findIndex((element, order) => order > index && element.depth <= area[index].depth);
    return area.slice(index + 1, end < 0 ? undefined : end).filter((element) => element.depth === area[index].depth + 1);
  };
  const block = parentOf(nested);
  const card = parentOf(block);
  assert.equal(area[block].name, 'div');
  assert.match(area[card].attributes, /\sv-settings-section=/);
  for (const [holder, path] of [[card, ''], [block, 'div > ']]) {
    const folded = tagsAfter(path);
    assert.ok(folded.size > 0, `tags of "${path}"`);
    assert.deepEqual([...new Set(childrenOf(holder).map((element) => element.name))].filter((name) => !folded.has(name)), [], `children of "${path || 'the card'}"`);
  }
});

test('the settings pane keeps no second copy of a section title beside its icon', () => {
  const pane = html.slice(html.indexOf('<div class="settings-pane'), html.indexOf('<!-- Resize Handle -->'));
  assert.ok(pane.includes('v-settings-section='), 'the settings pane was found');
  const copies = [];
  for (const { title, icon } of settingsEntriesOf('circular', 'rail').concat(settingsEntriesOf('linear', 'rail'))) {
    const escaped = title.replace(/&/g, '&amp;');
    for (const text of new Set([title, escaped])) {
      if (new RegExp(`<i class="ph ${icon}[^>]*></i>\\s*(?:<span[^>]*>)?${text.replace(/[.*+?^${}()|[\]\\]/g, '\\$&')}\\b`).test(pane)) {
        copies.push(`${icon} ${text}`);
      }
    }
  }
  assert.deepEqual(copies, []);
});

test('an index entry is named "Go to <title>" in one place, so it shares no name with the section it opens', () => {
  const template = html.slice(html.indexOf('id="settings-section-index-template"'));
  const index = template.slice(0, template.indexOf('</script>'));
  assert.equal(index.split(':aria-label="t(\'Go to {section}\', { section: tm(section.title) })"').length - 1, 1, 'the entry name');
  // The hover text of an entry stays the bare title; the rail shows the short name.
  assert.ok(index.includes('tipOn(tm(section.title))'));
  assert.ok(index.includes('{{ tm(section.short) }}'));
});
