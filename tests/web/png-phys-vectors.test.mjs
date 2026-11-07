// The pHYs chunk the web app's PNG export writes (setDpiInPng, utils/png.js) is the chunk Python writes
// (png_with_resolution in gbdraw/render/export.py): both read tests/fixtures/png_phys_vectors.json (tools/twins.json).
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import test from 'node:test';

import { setDpiInPng } from '../../gbdraw/web/js/utils/png.js';

/** @type {{ pixelsPerMetre: { dpi: number, pixelsPerMetre: number }[], chunks: { case: string, dpi: number, png: string, expected: string }[] }} */
const VECTORS = JSON.parse(readFileSync(new URL('../fixtures/png_phys_vectors.json', import.meta.url), 'utf8'));

/** @param {string} hex */
const pngBlob = (hex) => new Blob([Buffer.from(hex, 'hex')], { type: 'image/png' });
/** @param {Blob} blob */
const hexOf = async (blob) => Buffer.from(await blob.arrayBuffer()).toString('hex');

test('setDpiInPng writes the chunk Python writes (D-84: half-way pixels per metre round up)', async () => {
  for (const { case: name, dpi, png, expected } of VECTORS.chunks) {
    assert.equal(await hexOf(await setDpiInPng(pngBlob(png), dpi)), expected, name);
  }
  const bare = VECTORS.chunks[0].png;
  for (const { dpi, pixelsPerMetre } of VECTORS.pixelsPerMetre) {
    const bytes = Buffer.from(await (await setDpiInPng(pngBlob(bare), dpi)).arrayBuffer());
    const at = bytes.indexOf('pHYs');
    assert.deepEqual([bytes.readUInt32BE(at + 4), bytes.readUInt32BE(at + 8), bytes[at + 12]], [pixelsPerMetre, pixelsPerMetre, 1], String(dpi));
  }
});
