// D-44 (S-2): every entry that starts a Generate waits for the Similarity alignment work in flight, as the main
// Generate button does: Retry Generate, the label-fit dialog's "Apply and Generate", and the missing-BLAST-pairs
// dialog's "Skip these comparisons and Generate" and "Generate". Each button is disabled with the work's reason, read
// from the one work-in-flight owner (`similarityAlignmentBusy`, `similarityAlignmentWorkReason`); an action reached
// while the work runs changes nothing and returns the work's reason.
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import test from 'node:test';

globalThis.window = {
  Vue: { ref: (value) => ({ value }), reactive: (value) => value, computed: (getter) => ({ get value() { return getter(); } }) }
};
const { createMissingBlastTsvChoice } = await import('../../gbdraw/web/js/app/result-actions.js');

const INDEX = readFileSync(new URL('../../gbdraw/web/index.html', import.meta.url), 'utf8');
/** The opening tag of each button whose click runs `handler` (a quoted value may hold `>`). @param {RegExp} handler */
const buttonsCalling = (handler) => [...INDEX.matchAll(/<button\b(?:[^>"]|"[^"]*")*>/gs)].map(([tag]) => tag)
  .filter((tag) => handler.test(tag.match(/@click="([^"]*)"/)?.[1] ?? ''));

for (const [entry, handler] of [
  ['Retry Generate and the main Generate button', /^runAnalysis$/],
  ['the label-fit dialog\'s Apply and Generate', /^resolveLabelFitChoice\('apply'\)$/],
  ['the missing-BLAST-pairs dialog\'s Skip and Generate', /^resolveMissingBlastTsv\('none'\)$/],
  ['the missing-BLAST-pairs dialog\'s Generate', /^resolveMissingBlastTsv\('generate'\)$/]
]) {
  test(`${entry} is disabled with the alignment work's reason`, () => {
    const buttons = buttonsCalling(handler);
    assert.ok(buttons.length > 0, `a button calls ${handler}`);
    for (const tag of buttons) {
      assert.match(tag, /:disabled="[^"]*\bsimilarityAlignmentBusy\b[^"]*"/, tag);
      assert.match(tag, /:title="similarityAlignmentBusy \? similarityAlignmentWorkReason : null"/, tag);
    }
  });
}

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
