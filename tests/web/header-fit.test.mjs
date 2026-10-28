import assert from 'node:assert/strict';
import { test } from 'node:test';

import { HEADER_FITS, SESSION_TITLE_MIN_WIDTH, chooseHeaderFit } from '../../gbdraw/web/js/app/header-fit.js';

// #1383 (PD-OI-108): the header keeps the first state that leaves the tagline
// room a pixel or more: all words, then the Gallery cap alone, then the file
// buttons' icons with the Gallery words, then icons alone, then icons alone
// and the Session title without its caption, then also the Language menu as
// wide as the shown name. OV-2177, OV-2178: where every one-row state narrows
// the title, the last keeps one row while the title is at least
// SESSION_TITLE_MIN_WIDTH wide; narrower, the header wraps.
test('the header keeps the first state, in its order, that narrows nothing', () => {
  const fileWords = 400;
  const galleryWords = 150;
  const caption = 60;
  const wideLanguage = 50;
  const title = 128;
  // What the header leaves without any words; a state's words take from the
  // room, then from the title.
  const measureIn = (free) => (fit) => {
    const left = free - (fit.fileWords ? fileWords : 0) - (fit.galleryWords ? galleryWords : 0)
      - (fit.sessionCaption ? caption : 0) - (fit.languageWide ? wideLanguage : 0);
    return { room: Math.max(0, left), title: Math.max(10, Math.min(title, title + left)) };
  };
  const choose = (free) => HEADER_FITS.indexOf(chooseHeaderFit(measureIn(free)));
  assert.equal(choose(661), 0);
  assert.equal(choose(660), 1);
  assert.equal(choose(511), 1);
  assert.equal(choose(510.5), 2);
  assert.equal(choose(261), 2);
  assert.equal(choose(260), 3);
  assert.equal(choose(111), 3);
  // The icons alone keep the caption while they leave the room any width.
  assert.equal(choose(110.5), 3);
  assert.equal(choose(110), 4);
  assert.equal(choose(51), 4);
  // Without the caption the title narrows: the Language menu narrows first.
  assert.equal(choose(50), 5);
  assert.equal(choose(1), 5);
  // The last one-row state narrows the title, down to its minimum.
  assert.equal(choose(SESSION_TITLE_MIN_WIDTH - title), 5);
  assert.equal(choose(SESSION_TITLE_MIN_WIDTH - title - 0.5), 6);
  assert.equal(choose(-200), 6);
  assert.deepEqual(HEADER_FITS.map(({ fileWords: file, galleryWords: gallery, sessionCaption, languageWide, oneRow }) =>
    [file, gallery, sessionCaption, languageWide, oneRow]), [
    [true, true, true, true, true], [true, false, true, true, true], [false, true, true, true, true],
    [false, false, true, true, true], [false, false, false, true, true], [false, false, false, false, true],
    [false, false, false, false, false]
  ]);
});
