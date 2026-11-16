// D-44 (S-2), D-57 (S-b): every entry that starts a Generate waits for the Similarity alignment work in flight, as the
// main Generate button does: Retry Generate, the error card's Generate, the label-fit dialog's "Apply and Generate",
// and the missing-BLAST-pairs dialog's "Skip these comparisons and Generate" and "Generate". Each control is disabled
// with the work's reason, read from the one work-in-flight owner (`similarityAlignmentBusy`,
// `similarityAlignmentWorkReason`); an action reached while the work runs changes nothing and returns the work's reason.
// The template's controls are found by their handlers, and the code that starts a Generate by its calls of
// `runAnalysis(`, so a new entry without the reason fails here.
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

// The code that calls `runAnalysis(` (the Generate the template's entries start), by file: the error card's
// `generateFromError` and the label-fit port in the root, the label-fit dialog's Apply, and the missing-BLAST-pairs
// choices. A new caller is a new way to start a Generate: give its template control the work's reason and list it.
test('the code that starts a Generate is the code of the listed entries', () => {
  const js = fileURLToPath(new URL('../../gbdraw/web/js/', import.meta.url));
  const callers = Object.fromEntries(readdirSync(js, { recursive: true }).map(String).filter((file) => file.endsWith('.js'))
    .map((file) => [file, (readFileSync(join(js, file), 'utf8').match(/(?<![\w.])runAnalysis\(/g) ?? []).length])
    .filter(([, count]) => count > 0));
  assert.deepEqual(callers, { 'app/app-setup.js': 2, 'app/label-fit-dialog.js': 1, 'app/result-actions.js': 2 });
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
