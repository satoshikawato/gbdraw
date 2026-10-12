// The settings-section catalog (services/settings-sections.js) is the one
// owner of a section's key, title, icon, and modes: index.html names a section
// only by its key, once on the card and once in its heading, in catalog order.
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { join } from 'node:path';
import test from 'node:test';

import { SETTINGS_SECTIONS, settingsSectionsOf } from '../../gbdraw/web/js/services/settings-sections.js';

const html = readFileSync(join(process.cwd(), 'gbdraw', 'web', 'index.html'), 'utf8');
const keys = SETTINGS_SECTIONS.map((section) => section.key);

test('each mode lists its sections in card order', () => {
  assert.deepEqual(settingsSectionsOf('circular').map((section) => section.title), [
    'Input Genomes', 'Region Annotations', 'Pairwise Comparisons', 'Layout', 'Titles & Record Labels', 'Legend',
    'Colors', 'Features', 'Labels', 'Axis & Scale', 'Dinucleotide content/skew'
  ]);
  assert.deepEqual(settingsSectionsOf('linear').map((section) => section.title), [
    'Input Genomes', 'Comparison', 'Advanced comparison and layout', 'Region Annotations', 'Layout',
    'Titles & Record Labels', 'Legend', 'Colors', 'Features', 'Labels', 'Axis & Scale', 'Dinucleotide content/skew'
  ]);
  assert.equal(new Set(keys).size, keys.length, 'section keys are unique');
  for (const section of SETTINGS_SECTIONS) {
    assert.match(section.icon, /^ph-[a-z-]+$/, section.key);
    assert.ok(section.modes.length > 0, section.key);
  }
});

test('index.html names each section by its key, once on the card and once in the heading, in catalog order', () => {
  const marked = Array.from(html.matchAll(/\sdata-settings-section="([^"]*)"/g), (match) => match[1]);
  const headings = Array.from(html.matchAll(/<settings-section-title section="([^"]*)">/g), (match) => match[1]);
  assert.deepEqual(marked, keys, 'data-settings-section values');
  assert.deepEqual(headings, keys, 'settings-section-title headings');
  // Each heading sits inside its own card: the key attribute comes first, and
  // the next card starts after the heading.
  const offsets = keys.map((key) => [
    html.indexOf(`data-settings-section="${key}"`), html.indexOf(`<settings-section-title section="${key}">`)
  ]);
  for (const [order, [card, heading]] of offsets.entries()) {
    assert.ok(card < heading, `${keys[order]}: heading after its key`);
    if (order > 0) assert.ok(offsets[order - 1][1] < card, `${keys[order]}: after the previous heading`);
  }
});

test('the settings pane keeps no second copy of a section title beside its icon', () => {
  const pane = html.slice(html.indexOf('<div class="settings-pane'), html.indexOf('<!-- Resize Handle -->'));
  assert.ok(pane.includes('data-settings-section='), 'the settings pane was found');
  const copies = [];
  for (const { title, icon } of SETTINGS_SECTIONS) {
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
  assert.equal(index.split(':aria-label="\'Go to \' + section.title"').length - 1, 1, 'the entry name');
  // The hover text stays the bare title.
  assert.ok(index.includes('{{ tip.title }}'));
});
