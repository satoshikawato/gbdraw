import assert from 'node:assert/strict';
import { test } from 'node:test';

import { HEADER_FITS, chooseHeaderFit } from '../../gbdraw/web/js/app/header-fit.js';

// #1383 (PD-OI-108): the header keeps the first state that leaves the tagline
// room a pixel or more: all words, then the Gallery cap alone, then the file
// buttons' icons with the Gallery words, then icons alone, then icons alone
// and the Session title without its caption.
test('the header keeps the first state, in its order, that narrows nothing', () => {
  const fileWords = 400;
  const galleryWords = 150;
  const caption = 60;
  // The room the header leaves without any words; a state's words take from it.
  const roomIn = (free) => (fit) => Math.max(0, free - (fit.fileWords ? fileWords : 0)
    - (fit.galleryWords ? galleryWords : 0) - (fit.sessionCaption ? caption : 0));
  const choose = (free) => HEADER_FITS.indexOf(chooseHeaderFit(roomIn(free)));
  assert.equal(choose(611), 0);
  assert.equal(choose(610), 1);
  assert.equal(choose(461), 1);
  assert.equal(choose(460.5), 2);
  assert.equal(choose(211), 2);
  assert.equal(choose(210), 3);
  assert.equal(choose(61), 3);
  // The icons alone keep the caption while they leave the room any width.
  assert.equal(choose(60.5), 3);
  assert.equal(choose(60), 4);
  assert.equal(choose(0), 4);
  assert.deepEqual(HEADER_FITS.map(({ fileWords: file, galleryWords: gallery, sessionCaption }) => [file, gallery, sessionCaption]),
    [[true, true, true], [true, false, true], [false, true, true], [false, false, true], [false, false, false]]);
});
