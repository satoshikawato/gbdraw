// UJ-06 (Owner 2026-10-08): the Gallery catalog this origin serves is the one
// source of Open Gallery example: the button shows when it gives an example, and
// the chooser lists its usable entries in catalog order.
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import test from 'node:test';

import {
  galleryExamplesFromCatalog,
  galleryTitleParts
} from '../../gbdraw/web/js/services/gallery-examples.js';
import { isReusableLosatEntry } from '../../gbdraw/web/js/services/losat-cache.js';
import { readSessionContainer } from '../../gbdraw/web/js/services/session-container.js';
import { normalizeOptionalHexColor } from '../../gbdraw/web/js/utils/color-utils.js';

const CATALOG_URL = new URL('https://gbdraw.app/gallery/examples.json');
const catalog = JSON.parse(readFileSync(new URL('../../gbdraw/web/gallery/examples.json', import.meta.url), 'utf8'));
const galleryFile = (path) => new URL(`../../gbdraw/web/gallery/${path.replace(/^\.\//, '')}`, import.meta.url);
const sessions = new Map();
/**
 * The document of a catalog entry's Session file (a `.gbdraw` container), read once for every test of this file; its
 * inputs stay unread.
 */
const gallerySession = (entry) => {
  if (!sessions.has(entry.session)) {
    const file = new Blob([readFileSync(galleryFile(entry.session))]);
    sessions.set(entry.session, readSessionContainer(file).then(({ document }) => document));
  }
  return sessions.get(entry.session);
};

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
  // Load reads a Gallery Session by its content; its name is one the Load Session picker takes.
  assert.ok(examples.every((example) => /\.(?:gbdraw|zip|json|gz)$/.test(example.sessionName)));
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
    { ...valid, id: 'b', mode: 'linear', thumbnail: 'https://example.org/b.webp' },
    { ...valid, id: 'c', session: './sessions/c.gbdraw' }
  ], CATALOG_URL);
  assert.deepEqual(examples.map(({ id, mode, thumbnailUrl }) => ({ id, mode, thumbnailUrl })), [
    { id: 'a', mode: '', thumbnailUrl: '' },
    { id: 'b', mode: 'Linear', thumbnailUrl: '' },
    { id: 'c', mode: '', thumbnailUrl: '' }
  ]);
});

// OV-3073: a tutorial replays its steps from Reset Settings, so a ring the
// reader adds takes the Web's automatic colour of its row. The published
// Session must draw the same colours unless the tutorial sets one, or the
// tutorial's previews and the Gallery figure disagree (as after #1803).
// PD-OI-110: a ring's colour in a Gallery Session is its own exactly when the
// tutorial sets it (D-15 keeps the E. coli tutorial's three colours); a ring
// the tutorial leaves alone takes the Ring colors scheme (null).
test('a Gallery Session ring has its own colour exactly when its tutorial sets one', async () => {
  let checked = 0;
  for (const entry of catalog) {
    const tutorial = JSON.parse(readFileSync(galleryFile(entry.tutorial), 'utf8'));
    /** @type {Map<string, string>} The colour the tutorial sets, by the ring's file name. */
    const setColours = new Map((tutorial.manualSteps ?? []).flatMap((step) => (step.replay ?? []))
      .flatMap((action) => {
        const match = String(action.selector ?? '').match(/Comparison ring color for (.+?)'/);
        return match ? [[match[1], String(action.value).toLowerCase()]] : [];
      }));
    for (const drawing of (await gallerySession(entry)).drawings ?? []) {
      const series = drawing.config?.circularConservation?.series ?? [];
      series.forEach((row, index) => {
        const name = [...setColours.keys()].find((fileName) => String(row.fileName).startsWith(fileName.replace(/\.[^.]+$/, '')));
        assert.equal(normalizeOptionalHexColor(row.color), name === undefined ? null : setColours.get(name),
          `${entry.id}: comparison ring ${index + 1}`);
        checked += 1;
      });
    }
  }
  assert.ok(checked > 0, 'no Gallery Session has comparison rings');
});

// A Gallery example draws on first Generate from the LOSAT results its Session
// stores, with no search. Generate reuses a stored result only when this Web
// build would (isReusableLosatEntry): rows of another LOSAT Wasm build are
// searched again, so a LOSAT update re-searches them in the Session with
// tools/refresh_gallery_sessions.py.
test('every LOSAT result a Gallery Session stores is one that Generate reuses', async () => {
  let checked = 0;
  for (const entry of catalog) {
    for (const drawing of (await gallerySession(entry)).drawings ?? []) {
      for (const result of drawing.losatCache?.entries ?? []) {
        assert.ok(isReusableLosatEntry(result), `${entry.id}: ${result.filename ?? result.key} (${JSON.stringify(result.runtime)}); `
          + `search it again: python tools/refresh_gallery_sessions.py --session ${entry.id}`);
        checked += 1;
      }
    }
  }
  assert.ok(checked > 0, 'no Gallery Session stores a LOSAT result');
});
