import assert from 'node:assert/strict';
import { test } from 'node:test';

import { HEADER_FITS, chooseHeaderFit } from '../../gbdraw/web/js/app/header-fit.js';

// #1383 (PD-OI-108): the header keeps the first state that leaves the tagline
// room a pixel or more: all words, then the Gallery cap alone, then the file
// buttons' icons with the Gallery words, then icons alone.
test('the header keeps the first state, in its order, that narrows nothing', () => {
  const fileWords = 400;
  const galleryWords = 150;
  // The room the header leaves without any words; a state's words take from it.
  const roomIn = (free) => (fit) => Math.max(0, free - (fit.fileWords ? fileWords : 0) - (fit.galleryWords ? galleryWords : 0));
  const choose = (free) => HEADER_FITS.indexOf(chooseHeaderFit(roomIn(free)));
  assert.equal(choose(551), 0);
  assert.equal(choose(550), 1);
  assert.equal(choose(401), 1);
  assert.equal(choose(400.5), 2);
  assert.equal(choose(151), 2);
  assert.equal(choose(150), 3);
  assert.equal(choose(0), 3);
  assert.deepEqual(HEADER_FITS.map(({ fileWords: file, galleryWords: gallery }) => [file, gallery]),
    [[true, true], [true, false], [false, true], [false, false]]);
});
