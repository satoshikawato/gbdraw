// UJ-06 (Owner 2026-10-08): the Gallery catalog this origin serves is the one
// source of Open Gallery example: the button shows when it gives an example, and
// the chooser lists its usable entries in catalog order.
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import test from 'node:test';
import { gunzipSync } from 'node:zlib';

import {
  galleryExamplesFromCatalog,
  galleryTitleParts
} from '../../gbdraw/web/js/services/gallery-examples.js';
import { normalizeConservationSeriesColor } from '../../gbdraw/web/js/services/conservation-series.js';

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
    // #1383 (PD-OI-108): every shipped example links to the Tutorial tab of its Gallery page.
    assert.equal(entry.tutorialStatus, 'ready');
    assert.equal(example.tutorialUrl, `https://gbdraw.app/gallery/${encodeURIComponent(entry.id)}#Tutorial`);
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

// #1383 (PD-OI-108): only a tutorial the Gallery publishes (named and ready) is
// linked, under the catalog's own Gallery path.
test('an example links to its tutorial only when the Gallery publishes it', () => {
  const valid = { id: 'a b', title: 'A', session: './sessions/a.gbdraw-session.json', tags: [] };
  const tutorial = './tutorials/a.json';
  const examples = galleryExamplesFromCatalog([
    { ...valid, tutorial, tutorialStatus: 'ready' },
    { ...valid, tutorial },
    { ...valid, tutorial, tutorialStatus: 'draft' },
    { ...valid, tutorial, tutorialStatus: 'planned' },
    { ...valid, tutorialStatus: 'ready' },
    valid
  ], new URL('http://127.0.0.1:4173/gbdraw/web/gallery/examples.json'));
  assert.deepEqual(examples.map((example) => example.tutorialUrl), [
    'http://127.0.0.1:4173/gbdraw/web/gallery/a%20b#Tutorial',
    'http://127.0.0.1:4173/gbdraw/web/gallery/a%20b#Tutorial',
    '', '', '', ''
  ]);
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

// OV-3073: a tutorial replays its steps from Reset Settings, so a ring the
// reader adds takes the Web's automatic colour of its row. The published
// Session must draw the same colours unless the tutorial sets one, or the
// tutorial's previews and the Gallery figure disagree (as after #1803).
test('a Gallery Session whose tutorial sets no ring colour draws the automatic ring colours', () => {
  const galleryFile = (path) => new URL(`../../gbdraw/web/gallery/${path.replace(/^\.\//, '')}`, import.meta.url);
  let checked = 0;
  for (const entry of catalog) {
    const tutorial = JSON.parse(readFileSync(galleryFile(entry.tutorial), 'utf8'));
    const setsColour = (tutorial.manualSteps ?? []).some((step) => (step.replay ?? [])
      .some((action) => /ring colou?r/i.test(String(action.selector ?? ''))));
    if (setsColour) continue;
    const bytes = readFileSync(galleryFile(entry.session));
    const session = JSON.parse((entry.session.endsWith('.gz') ? gunzipSync(bytes) : bytes).toString('utf8'));
    for (const drawing of session.drawings ?? []) {
      const series = drawing.config?.circularConservation?.series ?? [];
      series.forEach((row, index) => {
        assert.equal(normalizeConservationSeriesColor(row.color, index), normalizeConservationSeriesColor('', index),
          `${entry.id}: comparison ring ${index + 1}`);
        checked += 1;
      });
    }
  }
  assert.ok(checked > 0, 'no Gallery Session has comparison rings');
});
