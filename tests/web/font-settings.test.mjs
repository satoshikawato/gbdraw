// The Font section's owner (app/font-settings.js): the select's entries, the
// custom-font note, loading a font (review R6), and the Session font question
// before a Session load (receipt diagram-generation.diagram-font rev 1,
// 〔読み込むとき〕).
import assert from 'node:assert/strict';
import test from 'node:test';

globalThis.window = {
  Vue: {
    ref: (value) => ({ value }),
    reactive: (value) => value,
    computed: (definition) => (typeof definition === 'function'
      ? { get value() { return definition(); } }
      : { get value() { return definition.get(); }, set value(next) { definition.set(next); } })
  }
};

// Load font file… opens a file input; the test answers with `nextFiles`.
/** @type {{ name: string, text: string, unreadable?: boolean }[]} */
let nextFiles = [];
globalThis.document = {
  body: { append: () => {} },
  createElement: () => {
    const listeners = new Map();
    const input = {
      dataset: {},
      files: null,
      remove: () => {},
      addEventListener: (type, listener) => { listeners.set(type, listener); },
      click: () => {
        input.files = nextFiles.map(({ name, text, unreadable = false }) => ({
          name,
          arrayBuffer: async () => {
            if (unreadable) throw new Error('NotReadableError');
            return new TextEncoder().encode(text).buffer;
          }
        }));
        queueMicrotask(() => listeners.get(input.files.length ? 'change' : 'cancel')?.());
      }
    };
    return input;
  }
};

const { createFontSettings } = await import('../../gbdraw/web/js/app/font-settings.js');
const { LEGEND_FONT_METRICS } = await import('../../gbdraw/web/js/utils/legend-font-metrics.generated.js');

const SERIF_TABLE = LEGEND_FONT_METRICS.faces['LiberationSerif-Regular'];
const quoted = (family) => `'${family}', serif`;

// The Worker stand-in names each font by its bytes: "family:<name>".
const registerFontsInWorker = async (fonts) => fonts.map((font) => {
  const family = new TextDecoder().decode(font.bytes).replace(/^family:/, '');
  return {
    identity: font.identity, family, fontFamily: quoted(family), weight: 400, italic: false, generic: 'serif',
    measurementName: `user-${font.identity.slice(0, 16)}-Regular`, metrics: SERIF_TABLE
  };
});

const setup = ({ fontFamily = null, loaded = [], answers = [] } = {}) => {
  const adv = { font_family: fontFamily };
  const state = { loadedFonts: { value: loaded }, activeDrawing: () => ({ adv }) };
  /** @type {object[]} */
  const questions = [];
  const steps = [];
  const settings = createFontSettings({
    state,
    runUndoable: async (label, fn) => { steps.push(label); return fn(); },
    registerFontsInWorker,
    askChoice: async (request) => { questions.push(request); return answers.shift() || 'cancel'; }
  });
  return { settings, state, adv, questions, steps };
};

const loadedEntry = (family) => ({
  identity: 'f'.repeat(64), family, weight: 400, italic: false, generic: 'serif', fontFamily: quoted(family),
  measurementName: 'user-ffffffffffffffff-Regular', source: 'local', label: `${family} Regular`, bytes: new ArrayBuffer(1), metrics: SERIF_TABLE
});

test('OV-685: a loaded font named like a bundled alias is the selected entry, not Liberation', () => {
  const { settings } = setup({ fontFamily: quoted('Times New Roman'), loaded: [loadedEntry('Times New Roman')] });
  assert.equal(settings.selectedChoice.value, quoted('Times New Roman'));
  assert.deepEqual(settings.choices.value.filter((choice) => choice.kind === 'loaded').map((choice) => choice.label), ['Times New Roman']);
  const unloaded = setup({ fontFamily: quoted('Times New Roman') });
  assert.equal(unloaded.settings.selectedChoice.value, "'Liberation Serif', serif", 'without the font, the alias is the bundled font');
});

test('a font this page cannot draw shows its name and the font it is drawn with', () => {
  const { settings } = setup({ fontFamily: quoted('Garamond Premier Pro') });
  const current = settings.choices.value.find((choice) => choice.kind === 'unavailable');
  assert.equal(current.label, 'Garamond Premier Pro — not available, drawn in Liberation Serif');
  assert.equal(settings.selectedChoice.value, quoted('Garamond Premier Pro'));
});

test('the custom-font note names any family gbdraw does not bundle, not only the first', () => {
  assert.equal(setup({ fontFamily: "'Liberation Serif', 'DejaVu Serif', serif" }).settings.customFontSelected.value, true);
  assert.equal(setup({ fontFamily: "'Arial', sans-serif" }).settings.customFontSelected.value, false);
  assert.equal(setup().settings.customFontSelected.value, false);
});

test('Load font file… reads, loads, and selects the font as one History step, loading from the read on (R6)', async () => {
  const { settings, adv, steps, state } = setup();
  nextFiles = [{ name: 'serif.ttf', text: 'family:Gbdraw Test Serif' }];
  const pending = settings.loadFontFile();
  assert.equal(settings.loading.value, true, 'loading while the file is read');
  await pending;
  assert.equal(settings.loading.value, false);
  assert.equal(settings.loadError.value, '');
  assert.deepEqual(steps, ['Change font']);
  assert.equal(adv.font_family, quoted('Gbdraw Test Serif'));
  assert.equal(state.loadedFonts.value.length, 1);
});

test('a file the browser cannot read is a font-file error, not an uncaught failure (R6)', async () => {
  const { settings, steps } = setup();
  nextFiles = [{ name: 'gone.ttf', text: '', unreadable: true }];
  await settings.loadFontFile();
  assert.match(settings.loadError.value, /TrueType or OpenType font file/);
  assert.equal(settings.loading.value, false);
  assert.deepEqual(steps, []);
});

test('Load with <nearest> loads the Session and keeps the names; Cancel loads nothing', async () => {
  const substitute = setup({ answers: ['substitute'] });
  assert.equal(await substitute.settings.confirmSessionFonts([quoted('Garamond Premier Pro'), null]), true);
  assert.equal(substitute.questions.length, 1);
  assert.equal(substitute.questions[0].title, 'This Session uses a font that is not available');
  const cancel = setup({ answers: ['cancel'] });
  assert.equal(await cancel.settings.confirmSessionFonts([quoted('Garamond Premier Pro')]), false);
  const none = setup();
  assert.equal(await none.settings.confirmSessionFonts([null, "'Liberation Mono', monospace"]), true);
  assert.equal(none.questions.length, 0, 'no question when every font resolves');
});

test('Load font file… in the question asks again, naming the family loaded, until the Session font is loaded', async () => {
  const { settings, questions, steps, adv } = setup({ answers: ['file', 'file'] });
  const files = [[{ name: 'other.ttf', text: 'family:DejaVu Sans' }], [{ name: 'garamond.otf', text: 'family:Garamond Premier Pro' }]];
  const original = globalThis.document.createElement;
  globalThis.document.createElement = () => { nextFiles = files.shift(); return original(); };
  try {
    assert.equal(await settings.confirmSessionFonts([quoted('Garamond Premier Pro')]), true);
  } finally {
    globalThis.document.createElement = original;
  }
  assert.equal(questions.length, 2);
  assert.deepEqual(questions[1].items, ['The loaded font is DejaVu Sans, not Garamond Premier Pro.']);
  assert.deepEqual(steps, [], 'loading fonts for a Session is no History step');
  assert.equal(adv.font_family, null, 'the current drawing keeps its font');
});
