// The Web port of Python's font-family rule (services/font-families.js) on the
// vectors Python generated (tools/generate_legend_font_metrics.py), the font
// store (services/font-store.js), and the Legend face of a registered font
// (services/legend-layout.js, OV-683).
import assert from 'node:assert/strict';
import { readFile } from 'node:fs/promises';
import test from 'node:test';

import { comparableFamilyName, customFontFamilies, resolveWebFontFamily } from '../../gbdraw/web/js/services/font-families.js';
import {
  fontIdentityOf,
  identifyWebFontSources,
  loadWebFonts,
  mergeWebFonts,
  registerWebFonts,
  registeredLegendFaces,
  standardInstalledFaces,
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
    assert.deepEqual([...resolution.replaced], vector.replaced, `vector ${index} replaced (OV-743)`);
  });
});

// R3 (b), D-10: one rule says which fonts a Session cannot carry, Python's
// `custom_font_families_among`, for the note under Font and the Save notice;
// a loaded alias name counts.
test('the port finds the custom families of each registration state as Python does', () => {
  assert.ok(vectors.vectors.some((vector) => vector.registered.length > 0 && vector.custom.length > 0));
  vectors.vectors.forEach((vector, index) => {
    assert.ok(Array.isArray(vector.custom), `vector ${index} has Python's custom families`);
    assert.deepEqual([...customFontFamilies(vector.requested, vector.registered)], vector.custom, `vector ${index}: ${vector.requested}`);
  });
});

const bytesOf = (text) => new TextEncoder().encode(text).buffer;
const SERIF_TABLE = LEGEND_FONT_METRICS.faces['LiberationSerif-Regular'];

const fakeWorker = (weightOf = () => 400) => async (fonts) => fonts.map((font, index) => ({
  identity: font.identity,
  family: 'Gbdraw Test Serif',
  fontFamily: "'Gbdraw Test Serif', serif",
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
    await identifyWebFontSources([{ bytes, source: 'file', label: 'a.ttf' }, { bytes: bytesOf('font bytes'), source: 'file', label: 'b.ttf' }]),
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
    registerWebFonts(await identifyWebFontSources([{ bytes: bytesOf('x'), source: 'file', label: 'x.ttf' }]), async () => [{ identity: 'other' }]),
    (error) => error.code === 'RESULT_INVALID'
  );
});

// R3 (a): the stored list of a family is the one Python wrote
// (`font_family_list_for`), not a JS copy of that rule.
test('an entry keeps the font-family list of the Python reply, and a reply without one is refused', async () => {
  const sources = await identifyWebFontSources([{ bytes: bytesOf('python list'), source: 'file', label: 'p.ttf' }]);
  const [entry] = await registerWebFonts(sources, async (fonts) => (await fakeWorker()(fonts))
    .map((fact) => ({ ...fact, fontFamily: '"Gbdraw Test Serif", serif' })));
  assert.equal(entry.fontFamily, '"Gbdraw Test Serif", serif');
  await assert.rejects(
    registerWebFonts(sources, async (fonts) => (await fakeWorker()(fonts)).map(({ fontFamily: _list, ...fact }) => fact)),
    (error) => error.code === 'RESULT_INVALID'
  );
});

// A FontFace stand-in: `load()` rejects bytes that start with "bad".
class FakeFontFace {
  constructor(family, bytes, descriptors = {}) {
    this.family = family;
    this.bytes = bytes;
    this.weight = descriptors.weight || 'normal';
    this.style = descriptors.style || 'normal';
  }

  async load() {
    if (new TextDecoder().decode(this.bytes).startsWith('bad')) throw new Error('rejected by the browser');
    return this;
  }
}
const fakeFontSet = () => ({ faces: [], add(face) { this.faces.push(face); } });

// R4: a font already in the store is not read, sent, or added again; the
// Worker gets the store's bytes (cloned once by postMessage, not sliced first).
test('loading a font already in the store registers nothing and selects the stored entry', async () => {
  const fontSet = fakeFontSet();
  /** @type {unknown[][]} */
  const calls = [];
  const registerInWorker = async (fonts) => { calls.push(fonts); return fakeWorker()(fonts); };
  const bytes = bytesOf('one font');
  const first = await loadWebFonts([{ bytes, source: 'file', label: 'one.ttf' }], {
    known: [], registerInWorker, fontSet, FontFaceClass: FakeFontFace
  });
  assert.equal(calls.length, 1);
  assert.equal(calls[0][0].bytes, bytes, 'the Worker request carries the store bytes, not a copy');
  assert.equal(first.chosen[0].bytes, bytes);
  assert.deepEqual(fontSet.faces.map((face) => [face.family, face.weight, face.style]), [['Gbdraw Test Serif', '400', 'normal']]);
  assert.deepEqual(Object.keys(first), ['added', 'chosen'], 'the caller merges `added` with mergeWebFonts');
  const again = await loadWebFonts([{ bytes: bytesOf('one font'), source: 'file', label: 'copy.ttf' }], {
    known: first.added, registerInWorker, fontSet, FontFaceClass: FakeFontFace
  });
  assert.equal(calls.length, 1, 'no second registration');
  assert.equal(fontSet.faces.length, 1, 'one FontFace per identity');
  assert.deepEqual(again.added, []);
  assert.equal(again.chosen[0], first.chosen[0]);
});

// M1: the one merge of a load into the page's fonts. A reloaded font (after
// the Worker lost it) replaces its entry; the other fonts keep their order.
test('mergeWebFonts replaces a reloaded entry and keeps the order of the rest', () => {
  const entry = (identity, label) => ({ identity, label });
  const merged = mergeWebFonts(
    [entry('a', 'a.ttf'), entry('b', 'b.ttf'), entry('c', 'c.ttf')],
    [entry('b', 'b-again.ttf'), entry('d', 'd.ttf')]
  );
  assert.deepEqual(merged.map((item) => item.label), ['a.ttf', 'c.ttf', 'b-again.ttf', 'd.ttf']);
  assert.deepEqual(mergeWebFonts([entry('a', 'a.ttf')], []).map((item) => item.label), ['a.ttf']);
});

// R7 (OV-687): the browser checks the font before Python registers it, so the
// store, Python, and the preview never disagree about a rejected font.
test('a font the browser rejects is not registered in Python', async () => {
  const fontSet = fakeFontSet();
  let registered = 0;
  await assert.rejects(
    loadWebFonts([{ bytes: bytesOf('bad font'), source: 'file', label: 'bad.ttf' }], {
      known: [], registerInWorker: async (fonts) => { registered += 1; return fakeWorker()(fonts); }, fontSet, FontFaceClass: FakeFontFace
    }),
    (error) => error.code === 'INPUT_INVALID' && error.context?.reason === 'FONT_FILE' && error.stage === 'font-validation'
  );
  assert.equal(registered, 0);
  assert.equal(fontSet.faces.length, 0);
});

// OV-685: a family the user loads wins over the bundled font of the same alias
// name (Times New Roman, Arial, Helvetica, Courier): it is measured with its own face.
test('OV-685: a loaded family that is a bundled alias is measured with the loaded face', () => {
  const metrics = withRegisteredFontFaces(createLegendFontMetrics(LEGEND_FONT_METRICS), [{
    family: 'Times New Roman', weight: 400, italic: false, measurementName: 'user-0123456789abcdef-Regular', table: SERIF_TABLE
  }]);
  assert.equal(resolveLegendFontFace(metrics, "'Times New Roman', serif"), 'user-0123456789abcdef-Regular');
  assert.equal(resolveLegendFontFace(createLegendFontMetrics(LEGEND_FONT_METRICS), "'Times New Roman', serif"), 'LiberationSerif-Regular');
});

test('the Legend measures a registered family with its face, and an unknown face by family (OV-683)', async () => {
  const entries = await registerWebFonts(
    await identifyWebFontSources([{ bytes: bytesOf('light'), source: 'local', label: 'Light' }, { bytes: bytesOf('regular'), source: 'local', label: 'Regular' }]),
    fakeWorker((index) => (index === 0 ? 300 : 400))
  );
  const metrics = withRegisteredFontFaces(createLegendFontMetrics(LEGEND_FONT_METRICS), registeredLegendFaces(entries));
  const face = resolveLegendFontFace(metrics, "'Gbdraw Test Serif', serif");
  assert.equal(face, entries[1].measurementName, 'CSS matching takes weight 400 over 300');
  assert.deepEqual(
    measureTextBox(metrics, { text: 'Hypothetical protein', fontFile: face, fontSize: 14 }),
    measureTextBox(metrics, { text: 'Hypothetical protein', fontFile: 'LiberationSerif-Regular', fontSize: 14 })
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

test('the Legend picks the face of a registered family Python picks for the normal weight and style', () => {
  assert.ok(vectors.faceMatch.length >= 10);
  vectors.faceMatch.forEach((vector, index) => {
    const faces = vector.faces.map((face, position) => ({
      family: 'Vector Family', weight: face.weight, italic: face.italic, measurementName: `face${position}`, table: SERIF_TABLE
    }));
    const metrics = withRegisteredFontFaces(createLegendFontMetrics(LEGEND_FONT_METRICS), faces);
    assert.equal(resolveLegendFontFace(metrics, "'Vector Family', serif"), `face${vector.normal}`, `faceMatch[${index}]`);
  });
});

// Review R2 / OV-685: the face Python measures a list with for each
// registration state; a registered family wins over the bundled alias of its
// name, generic keywords keep their Liberation face.
test('the Legend face of every vector is the face Python measures with', () => {
  vectors.vectors.forEach((vector, index) => {
    const faces = vector.registered.map((family) => ({
      family, weight: 400, italic: false, measurementName: `registered:${comparableFamilyName(family)}`, table: SERIF_TABLE
    }));
    const metrics = withRegisteredFontFaces(createLegendFontMetrics(LEGEND_FONT_METRICS), faces);
    const expected = vector.face.startsWith('registered:')
      ? `registered:${comparableFamilyName(vector.face.slice('registered:'.length))}`
      : vector.face;
    assert.equal(resolveLegendFontFace(metrics, vector.requested), expected, `vector ${index}: ${vector.requested}`);
  });
});

test('an installed family registers its Regular, Bold, Italic, and Bold Italic faces, upright first', () => {
  const face = (style) => ({ style });
  assert.deepEqual(
    standardInstalledFaces(['Light', 'Bold Italic', 'Black', 'Italic', 'Bold', 'Regular', 'Medium'].map(face)).map((entry) => entry.style),
    ['Regular', 'Bold', 'Italic', 'Bold Italic']
  );
  assert.deepEqual(standardInstalledFaces(['Book', 'Oblique', 'BoldOblique'].map(face)).map((entry) => entry.style), ['Book', 'Oblique', 'BoldOblique']);
  assert.deepEqual(standardInstalledFaces(['Light', 'Thin'].map(face)).map((entry) => entry.style), ['Light'], 'a family without them registers its first face');
  assert.deepEqual(standardInstalledFaces([]), []);
});
