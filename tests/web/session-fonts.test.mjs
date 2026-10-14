// The Session's custom-font notices (custom fonts D-02, D-03; receipt
// diagram-generation.diagram-font rev 1, 〔警告〕 and 〔読み込むとき〕): the Save
// Session warning, and the question Load asks before a Session whose font is
// not available replaces the current one. The texts are the approved ones
// (pr2/images/mock-q2q3/captures.json of the custom-fonts campaign).
import assert from 'node:assert/strict';
import { readFile, readdir } from 'node:fs/promises';
import { gunzipSync } from 'node:zlib';
import test from 'node:test';

import {
  SESSION_FONT_ACTIONS,
  customFontSaveWarnings,
  missingSessionFonts,
  sessionFontQuestion
} from '../../gbdraw/web/js/services/session-fonts.js';

test('Save warns once per font a Session cannot carry, with the approved text', () => {
  assert.deepEqual(customFontSaveWarnings([null, "'Liberation Serif', serif"], []), []);
  assert.deepEqual(customFontSaveWarnings(["'Arial', sans-serif", undefined], []), [], 'a bundled alias is carried by name');
  assert.equal(customFontSaveWarnings(["'Arial', sans-serif"], ['arial']).length, 1,
    'D-10: a loaded alias name is custom (here its file draws it; elsewhere Liberation)');
  assert.deepEqual(customFontSaveWarnings(["'DejaVu Serif', serif", "'DejaVu Serif', serif"], []), [
    'The Session keeps the font name DejaVu Serif, not the font file. To open it elsewhere, the same font must be available there, or it is replaced with the nearest available font.'
  ]);
  assert.equal(customFontSaveWarnings(["'Liberation Serif', 'DejaVu Serif', serif", "'Fira Code', monospace"], []).length, 2,
    'every family of every drawing that gbdraw does not bundle');
});

test('Load asks only when a named family is neither bundled nor loaded', () => {
  assert.equal(missingSessionFonts([null, "'Liberation Sans', 'Arial', 'Helvetica', 'Nimbus Sans L', sans-serif"], []), null);
  assert.equal(missingSessionFonts(["'DejaVu Serif', serif"], ['DejaVu Serif']), null);
  assert.deepEqual(missingSessionFonts(["'Garamond Premier Pro', serif", "'Garamond Premier Pro', serif"], []), {
    missing: ['Garamond Premier Pro'], nearest: ['Liberation Serif']
  });
  assert.deepEqual(missingSessionFonts(["'Garamond Premier Pro', serif", "'Fira Code', monospace"], []), {
    missing: ['Garamond Premier Pro', 'Fira Code'], nearest: ['Liberation Serif', 'Liberation Mono']
  }, 'one question for both drawings');
});

test('OV-743: Load asks only when the drawn family is not the requested first family', () => {
  assert.equal(missingSessionFonts(["'Liberation Serif', 'Foo', serif"], []), null, 'Liberation Serif is drawn; nothing is replaced');
  assert.deepEqual(missingSessionFonts(["'Foo', 'DejaVu Serif', 'Bar', serif"], ['DejaVu Serif']), {
    missing: ['Foo'], nearest: ['DejaVu Serif']
  }, 'names only the families before the one drawn');
});

test('the question has the approved title, texts, and buttons', () => {
  const missing = missingSessionFonts(["'Garamond Premier Pro', serif"], []);
  assert.deepEqual(sessionFontQuestion(missing, { installedFontsSupported: true, loadedFamily: null }), {
    title: 'This Session uses a font that is not available',
    description: 'Its text is set in Garamond Premier Pro, which is not part of gbdraw and is not loaded in this page.',
    items: [],
    note: 'Load with Liberation Serif: the Session opens and its text is drawn in Liberation Serif, the nearest available font. The Font setting keeps the name Garamond Premier Pro, so the font is used again once it is loaded. The saved preview keeps its look until you Generate.',
    actions: [
      { id: SESSION_FONT_ACTIONS.FILE, label: 'Load font file…', tone: 'primary' },
      { id: SESSION_FONT_ACTIONS.INSTALLED, label: 'Use installed font…', tone: 'secondary' },
      { id: SESSION_FONT_ACTIONS.SUBSTITUTE, label: 'Load with Liberation Serif', tone: 'secondary' }
    ],
    cancelLabel: 'Cancel'
  });
  const firefox = sessionFontQuestion(missing, { installedFontsSupported: false, loadedFamily: 'DejaVu Sans' });
  assert.deepEqual(firefox.actions.map((action) => action.label), ['Load font file…', 'Load with Liberation Serif'],
    'Use installed font… only with Local Font Access');
  assert.deepEqual(firefox.items, ['The loaded font is DejaVu Sans, not Garamond Premier Pro.']);
});

test('the question names every missing family and nearest font', () => {
  const question = sessionFontQuestion(
    missingSessionFonts(["'Garamond Premier Pro', serif", "'Fira Code', monospace"], []),
    { installedFontsSupported: false, loadedFamily: null }
  );
  assert.equal(question.description,
    'Its text is set in Garamond Premier Pro and Fira Code, which are not part of gbdraw and are not loaded in this page.');
  assert.equal(question.note, 'Load with Liberation Serif and Liberation Mono: the Session opens and its text is drawn in Liberation Serif and Liberation Mono, the nearest available fonts. The Font setting keeps the names Garamond Premier Pro and Fira Code, so the fonts are used again once they are loaded. The saved preview keeps its look until you Generate.');
  assert.equal(question.actions.at(-1).label, 'Load with Liberation Serif and Liberation Mono');
});

// Gallery Sessions never ask: every font list they hold resolves with no font loaded.
test('no Gallery Session names a font that is not bundled', async () => {
  const directory = new URL('../../gbdraw/web/gallery/sessions/', import.meta.url);
  const names = (await readdir(directory)).filter((name) => /\.gbdraw-session\.json(?:\.gz)?$/.test(name));
  assert.ok(names.length >= 5);
  for (const name of names) {
    const bytes = await readFile(new URL(name, directory));
    const text = (name.endsWith('.gz') ? gunzipSync(bytes) : bytes).toString('utf8');
    const lists = [...text.matchAll(/"font_family"\s*:\s*("(?:[^"\\]|\\.)*")/g)].map((match) => JSON.parse(match[1]));
    assert.equal(missingSessionFonts(lists, []), null, name);
  }
});
