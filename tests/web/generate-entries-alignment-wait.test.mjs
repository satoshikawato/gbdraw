// D-44 (S-2), D-57 (S-b): every entry that starts a Generate waits for the Similarity alignment work in flight, as the
// main Generate button does: Retry Generate, the error card's Generate, the label-fit dialog's "Apply and Generate",
// and the missing-BLAST-pairs dialog's "Skip these comparisons and Generate" and "Generate". Each control is disabled
// with the work's reason, read from the one work-in-flight owner (`similarityAlignmentBusy`,
// `similarityAlignmentWorkReason`); an action reached while the work runs changes nothing and returns the work's reason.
// The template's controls are found by their handlers, and the code that starts a Generate by its references to
// `runAnalysis` and `runGeneratedDiagramAnalysis`, so a new entry without the reason fails here.
import assert from 'node:assert/strict';
import { readdirSync, readFileSync } from 'node:fs';
import { join } from 'node:path';
import { fileURLToPath } from 'node:url';
import test from 'node:test';

globalThis.window = {
  Vue: { ref: (value) => ({ value }), reactive: (value) => value, computed: (getter) => ({ get value() { return getter(); } }) }
};
const { createMissingBlastTsvChoice } = await import('../../gbdraw/web/js/app/result-actions.js');

const INDEX = readFileSync(new URL('../../gbdraw/web/index.html', import.meta.url), 'utf8');
/** The template's opening tags with the handlers of their events (a quoted value may hold `>`). */
const CONTROLS = [...INDEX.matchAll(/<[a-z][\w-]*\b(?:[^>"]|"[^"]*")*>/gs)].map(([tag]) => ({
  tag, handlers: [...tag.matchAll(/\s(?:@|v-on:)[\w.:-]+="([^"]*)"/g)].map(([, handler]) => handler.trim())
}));
/** @param {RegExp} handler */
const controlsCalling = (handler) => CONTROLS.filter(({ handlers }) => handlers.some((each) => handler.test(each)))
  .map(({ tag }) => tag);
const assertWaits = (/** @type {string} */ tag) => {
  assert.match(tag, /:disabled="[^"]*\bsimilarityAlignmentBusy\b[^"]*"/, tag);
  assert.match(tag, /:title="similarityAlignmentBusy \? similarityAlignmentWorkReason : null"/, tag);
};

// The template's handlers that start a Generate.
const GENERATE_ENTRIES = [
  ['Retry Generate and the main Generate button', /^runAnalysis(?:\(\))?$/],
  ['the error card\'s Generate', /^generateFromError(?:\(\))?$/],
  ['the label-fit dialog\'s Apply and Generate', /^resolveLabelFitChoice\('apply'\)$/],
  ['the missing-BLAST-pairs dialog\'s Skip and Generate', /^resolveMissingBlastTsv\('none'\)$/],
  ['the missing-BLAST-pairs dialog\'s Generate', /^resolveMissingBlastTsv\('generate'\)$/]
];
for (const [entry, handler] of GENERATE_ENTRIES) {
  test(`${entry} is disabled with the alignment work's reason`, () => {
    const controls = controlsCalling(handler);
    assert.ok(controls.length > 0, `a control calls ${handler}`);
    controls.forEach(assertWaits);
  });
}

test('every template control whose handler starts a Generate is disabled with the alignment work\'s reason', () => {
  const starts = /\b(?:runAnalysis|generateFromError)\b|\bresolveLabelFitChoice\('apply'\)|\bresolveMissingBlastTsv\('(?:none|generate)'\)/;
  const controls = CONTROLS.filter(({ handlers }) => handlers.some((handler) => starts.test(handler)));
  assert.equal(controls.length, GENERATE_ENTRIES.reduce((count, [, handler]) => count + controlsCalling(handler).length, 0),
    'a control starts a Generate through a handler this test does not list');
  controls.forEach(({ tag }) => assertWaits(tag));
});

// The code that reaches the Generate the template's entries start: every reference, outside comment lines, to
// `runAnalysis` (the root's Generate, and the engine's, which the root binds as `runGeneratedDiagramAnalysis`) and to
// `runGeneratedDiagramAnalysis`, by file and by form: a call, a value (passed or returned, as `{ generate: runAnalysis }`
// would be), a declaration, or an object key. Today they are the root's Generate and its error card and label-fit
// port, the label-fit dialog's Apply, the missing-BLAST-pairs choices, and the engine's own. A new reference is a new
// way to start a Generate: give its template control the work's reason and list it.
const GENERATE_NAMES = ['runAnalysis', 'runGeneratedDiagramAnalysis'];
/** @param {string} source @returns {Record<string, Record<string, number>>} */
const generateReferences = (source) => {
  const code = source.split('\n').filter((line) => !/^\s*(?:\/\/|\/?\*)/.test(line)).join('\n');
  /** @type {Record<string, Record<string, number>>} */
  const forms = {};
  for (const name of GENERATE_NAMES) {
    for (const match of code.matchAll(new RegExp(`(?<![\\w$])${name}(?![\\w$])`, 'g'))) {
      const before = code.slice(Math.max(0, match.index - 20), match.index);
      const after = code.slice(match.index + name.length);
      const form = /\b(?:const|let|var|function)\s+$/.test(before) ? 'declaration'
        : /^\s*\(/.test(after) ? 'call' : /^\s*:/.test(after) ? 'key' : 'value';
      forms[name] = { ...forms[name], [form]: (forms[name]?.[form] ?? 0) + 1 };
    }
  }
  return forms;
};
test('the code that reaches a Generate is the code of the listed entries', () => {
  const js = fileURLToPath(new URL('../../gbdraw/web/js/', import.meta.url));
  const references = Object.fromEntries(readdirSync(js, { recursive: true }).map(String).filter((file) => file.endsWith('.js'))
    .sort().map((file) => [file, generateReferences(readFileSync(join(js, file), 'utf8'))])
    .filter(([, forms]) => Object.keys(forms).length > 0));
  assert.deepEqual(references, {
    // The root's Generate (declared, called by generateFromError and the label-fit port, given to the
    // missing-BLAST-pairs choice and the template), and the engine's, bound and called once.
    'app/app-setup.js': {
      runAnalysis: { declaration: 1, call: 2, key: 2, value: 2 },
      runGeneratedDiagramAnalysis: { value: 1, call: 1 }
    },
    // The label-fit dialog's Apply and the missing-BLAST-pairs choices: the port each takes and its calls.
    'app/label-fit-dialog.js': { runAnalysis: { value: 1, call: 1 } },
    'app/result-actions.js': { runAnalysis: { value: 1, call: 2 } },
    // The engine's Generate: declared and returned.
    'app/run-analysis.js': { runAnalysis: { declaration: 1, value: 1 } }
  });
});

const BUSY = Object.freeze({ status: 'busy', reason: 'Resolving the current alignment.' });
const missingBlastHarness = (busy) => {
  const steps = [];
  const omitted = [];
  const generates = [];
  const drawing = { linearComparisonResolution: { value: { errors: [{ code: 'missing-upload', edgeKey: 'a|b' }] } } };
  const choice = createMissingBlastTsvChoice({
    state: { activeDrawingId: { value: 'linear' }, drawingOf: () => drawing },
    withHistory: (label, handler) => async () => { steps.push(label()); return handler(); },
    closeAfterChoice: (close) => close(),
    choicePending: { value: false },
    omitLinearComparisonPairs: (_drawing, edgeKeys) => omitted.push(...edgeKeys),
    linearComparisonPairForEdgeKey: () => null,
    runAnalysis: async () => { generates.push(true); return { status: 'ok' }; },
    alignment: { busy: () => busy }
  });
  choice.missingBlastTsvDialog.show = true;
  return { choice, steps, omitted, generates };
};

test('the missing-BLAST-pairs choices change nothing and return the reason while alignment work runs', async () => {
  for (const option of ['generate', 'none']) {
    const { choice, steps, omitted, generates } = missingBlastHarness(BUSY);
    assert.deepEqual(await choice.resolveMissingBlastTsv(option), BUSY, option);
    assert.equal(choice.missingBlastTsvDialog.show, true, `${option}: the dialog stays open`);
    assert.deepEqual({ steps, omitted, generates }, { steps: [], omitted: [], generates: [] }, option);
  }
  const { choice, generates } = missingBlastHarness(null);
  assert.deepEqual(await choice.resolveMissingBlastTsv('generate'), { status: 'ok' }, 'with no work, it Generates');
  assert.equal(generates.length, 1);
});
