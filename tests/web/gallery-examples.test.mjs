// UJ-06 (Owner 2026-10-08): the Gallery catalog this origin serves is the one
// source of Load an example: the button shows when it gives an example, and
// the chooser lists its usable entries in catalog order.
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import test from 'node:test';

import {
  galleryExamplesFromCatalog,
  galleryTitleParts
} from '../../gbdraw/web/js/services/gallery-examples.js';

const CATALOG_URL = new URL('https://gbdraw.app/gallery/examples.json');
const catalog = JSON.parse(readFileSync(new URL('../../gbdraw/web/gallery/examples.json', import.meta.url), 'utf8'));

test('every shipped catalog entry is an example, in catalog order, with its Session and mode', () => {
  const examples = galleryExamplesFromCatalog(catalog, CATALOG_URL);
  assert.deepEqual(examples.map((example) => example.id), catalog.map((entry) => entry.id));
  for (const [index, example] of examples.entries()) {
    const entry = catalog[index];
    assert.equal(example.sessionUrl, new URL(entry.session, CATALOG_URL).href);
    assert.equal(example.sessionName, entry.session.split('/').pop());
    assert.equal(example.thumbnailUrl, new URL(entry.thumbnail, CATALOG_URL).href);
    assert.equal(example.mode, entry.tags.includes('Linear') ? 'Linear' : 'Circular');
  }
  assert.ok(examples.some((example) => example.sessionName.endsWith('.json.gz')));
});

test('a title shows its <i> runs as italics and every other character as text', () => {
  assert.deepEqual(galleryTitleParts('<i>Vibrio parahaemolyticus</i> and <i>V. alginolyticus</i> collinearity'), [
    { text: 'Vibrio parahaemolyticus', italic: true },
    { text: ' and ', italic: false },
    { text: 'V. alginolyticus', italic: true },
    { text: ' collinearity', italic: false }
  ]);
  assert.deepEqual(galleryTitleParts('<b>x</b> & <img src=x>'), [{ text: '<b>x</b> & <img src=x>', italic: false }]);
});

test('an unusable catalog or entry gives no example', () => {
  for (const value of [null, {}, 'x', [{}]]) assert.deepEqual(galleryExamplesFromCatalog(value, CATALOG_URL), []);
  const valid = { id: 'a', title: 'A', session: './sessions/a.gbdraw-session.json', tags: [] };
  const examples = galleryExamplesFromCatalog([
    valid,
    { ...valid, id: '' },
    { ...valid, title: ' ' },
    { ...valid, session: 'https://example.org/a.gbdraw-session.json' },
    { ...valid, session: './sessions/a.svg' },
    { ...valid, id: 'b', mode: 'linear', thumbnail: 'https://example.org/b.webp' }
  ], CATALOG_URL);
  assert.deepEqual(examples.map(({ id, mode, thumbnailUrl }) => ({ id, mode, thumbnailUrl })), [
    { id: 'a', mode: '', thumbnailUrl: '' },
    { id: 'b', mode: 'Linear', thumbnailUrl: '' }
  ]);
});
