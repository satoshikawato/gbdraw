// The Web port of Python's font-family rule (services/font-families.js) on the
// vectors Python generated (tools/generate_legend_font_metrics.py), the font
// store (services/font-store.js), and the Legend face of a registered font
// (services/legend-layout.js, OV-683).
import assert from 'node:assert/strict';
import { readFile } from 'node:fs/promises';
import test from 'node:test';

import { fontFamilyListFor, resolveWebFontFamily } from '../../gbdraw/web/js/services/font-families.js';
import {
  fontIdentityOf,
  mergeWebFonts,
  registerWebFonts,
  registeredLegendFaces,
  webFontFamilies,
  webFontIdentitiesFor
} from '../../gbdraw/web/js/services/font-store.js';
import {
  createLegendFontMetrics,
  legendMeasurementFace,
  measureTextBox,
  resolveLegendFontFace,
  withRegisteredFontFaces
} from '../../gbdraw/web/js/services/legend-layout.js';
import { LEGEND_FONT_METRICS } from '../../gbdraw/web/js/utils/legend-font-metrics.generated.js';

const vectors = JSON.parse(await readFile(new URL('../fixtures/font_family_resolution_vectors.json', import.meta.url), 'utf8'));

test('the port resolves every Python vector exactly', () => {
  assert.equal(vectors.schema, 1);
  assert.ok(vectors.vectors.length >= 20);
  vectors.vectors.forEach((vector, index) => {
    const resolution = resolveWebFontFamily(vector.requested, vector.registered);
    assert.equal(resolution.effective, vector.effective, `vector ${index}: ${vector.requested}`);
    assert.deepEqual([...resolution.missing], vector.missing, `vector ${index} missing`);
  });
});

test('a stored list quotes the family as Python does', () => {
  assert.equal(fontFamilyListFor('DejaVu Serif', 'serif'), "'DejaVu Serif', serif");
  assert.equal(fontFamilyListFor("O'Brien Sans", 'sans-serif'), "\"O'Brien Sans\", sans-serif");
});

const bytesOf = (text) => new TextEncoder().encode(text).buffer;
const SERIF_TABLE = LEGEND_FONT_METRICS.faces['LiberationSerif-Regular'];

const fakeWorker = (weightOf = () => 400) => async (fonts) => fonts.map((font, index) => ({
  identity: font.identity,
  family: 'Gbdraw Test Serif',
  weight: weightOf(index),
  italic: false,
  generic: 'serif',
  measurementName: `user-${font.identity.slice(0, 16)}-Regular`,
  metrics: SERIF_TABLE
}));

test('the store keys fonts by the SHA-256 of their bytes and registers identical bytes once', async () => {
  const bytes = bytesOf('font bytes');
  const identity = await fontIdentityOf(bytes);
  assert.match(identity, /^[0-9a-f]{64}$/);
  let requested = 0;
  const entries = await registerWebFonts(
    [{ bytes, source: 'file', label: 'a.ttf' }, { bytes: bytesOf('font bytes'), source: 'file', label: 'b.ttf' }],
    async (fonts) => { requested = fonts.length; return fakeWorker()(fonts); }
  );
  assert.equal(requested, 1);
  assert.equal(entries.length, 1);
  assert.equal(entries[0].identity, identity);
  assert.equal(entries[0].fontFamily, "'Gbdraw Test Serif', serif");
  assert.equal(entries[0].label, 'a.ttf');
  assert.deepEqual(webFontFamilies(entries).map((entry) => entry.family), ['Gbdraw Test Serif']);
  assert.deepEqual(webFontIdentitiesFor(entries, "'gbdraw test serif', serif"), [identity]);
  assert.deepEqual(webFontIdentitiesFor(entries, "'Liberation Sans', sans-serif"), []);
  assert.deepEqual(mergeWebFonts(entries, entries).map((entry) => entry.identity), [identity]);
});

test('a malformed Worker reply is refused', async () => {
  await assert.rejects(
    registerWebFonts([{ bytes: bytesOf('x'), source: 'file', label: 'x.ttf' }], async () => [{ identity: 'other' }]),
    (error) => error.code === 'RESULT_INVALID'
  );
});

test('the Legend measures a registered family with its face, and an unknown face by family (OV-683)', async () => {
  const entries = await registerWebFonts(
    [{ bytes: bytesOf('light'), source: 'local', label: 'Light' }, { bytes: bytesOf('regular'), source: 'local', label: 'Regular' }],
    fakeWorker((index) => (index === 0 ? 300 : 400))
  );
  const metrics = withRegisteredFontFaces(createLegendFontMetrics(LEGEND_FONT_METRICS), registeredLegendFaces(entries));
  const face = resolveLegendFontFace(metrics, "'Gbdraw Test Serif', serif");
  assert.equal(face, entries[1].measurementName, 'CSS matching takes weight 400 over 300');
  assert.deepEqual(
    measureTextBox(metrics, { text: 'Hypothetical protein', fontFile: face, fontSize: 14, dpi: 72 }),
    measureTextBox(metrics, { text: 'Hypothetical protein', fontFile: 'LiberationSerif-Regular', fontSize: 14, dpi: 72 })
  );
  const replaced = withRegisteredFontFaces(createLegendFontMetrics(LEGEND_FONT_METRICS), [
    ...registeredLegendFaces(entries), { ...registeredLegendFaces(entries)[1], measurementName: 'user-later-Regular' }
  ]);
  assert.equal(resolveLegendFontFace(replaced, "'Gbdraw Test Serif', serif"), 'user-later-Regular', 'a later face of the same weight and style replaces');
  const bundled = createLegendFontMetrics(LEGEND_FONT_METRICS);
  assert.equal(resolveLegendFontFace(bundled, "'Gbdraw Test Serif', serif"), 'LiberationSerif-Regular');
  // A face recorded on another machine (OV-683) falls back to the registered family's face, not a bundled one.
  assert.equal(legendMeasurementFace(metrics, 'GbdrawTestSerif-Regular', "'Gbdraw Test Serif', serif"), entries[1].measurementName);
  assert.equal(legendMeasurementFace(metrics, entries[0].measurementName, "'Gbdraw Test Serif', serif"), entries[0].measurementName);
});
