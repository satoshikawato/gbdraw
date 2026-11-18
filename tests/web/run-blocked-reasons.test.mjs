// Carried (b): a control that waits for a run (a Generate, a Rotate or Align run, a Session Save or Load) is disabled
// with the reason of the one owner of that availability, `sessionOperationAvailability('not-during-run')`. A control
// whose `:disabled` reads it names, in its aria-describedby, a `<run-blocked-reason>` (components.js), which shows that
// owner's reason; no control copies the reason or reads the availability without showing it. The Similarity
// alignment owner decides its own Align, Review and Apply (`actionBlock` in app/similarity-alignment.js) and words the
// run's reason on their reason lines, so those controls read only the owner.
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import test from 'node:test';

const INDEX = readFileSync(new URL('../../gbdraw/web/index.html', import.meta.url), 'utf8');
/** The template's opening tags (a quoted value may hold `>`). */
const TAGS = [...INDEX.matchAll(/<[a-z][\w-]*\b(?:[^>"]|"[^"]*")*>/gs)].map(([tag]) => tag);
const RUN_AVAILABILITY = "sessionOperationAvailability('not-during-run')";
const REASON_IDS = new Set(TAGS.filter((tag) => tag.startsWith('<run-blocked-reason'))
  .map((tag) => /\sid="([^"]+)"/.exec(tag)?.[1]).filter(Boolean));

test('every control disabled while a run holds it names the run\'s reason in its description', () => {
  const controls = TAGS.filter((tag) => (/\s:disabled="([^"]*)"/.exec(tag)?.[1] || '').includes(RUN_AVAILABILITY));
  assert.ok(controls.length >= 6, `the run-blocked controls are found (${controls.length})`);
  for (const tag of controls) {
    const described = (/\saria-describedby="([^"]*)"/.exec(tag)?.[1] || '').split(/\s+/).filter(Boolean);
    assert.ok(described.some((id) => REASON_IDS.has(id)),
      `a run-blocked control names a <run-blocked-reason> in aria-describedby: ${tag.replace(/\s+/g, ' ')}`);
  }
});

// An aria-describedby that names a missing element describes nothing (accessibility.playwright.spec.js, 'aria-describedby
// without a target'): the element with a reason's id is there whenever a control that names it is.
test('a run-blocked reason element is there whenever a control that names it is', () => {
  const components = readFileSync(new URL('../../gbdraw/web/js/components.js', import.meta.url), 'utf8');
  const component = /export const RunBlockedReason = \{[\s\S]*?\n\};/.exec(components)?.[0] || '';
  const root = /template: `(<[a-z][\w-]*\b[^>]*>)/.exec(component)?.[1] || '';
  assert.match(root, /\s:id="id"/, 'the root element carries the id');
  assert.doesNotMatch(root, /\sv-(if|else|else-if|for)=/, `the root element with the id is always rendered: ${root}`);
  const conditionOf = (tag) => /\sv-if="([^"]*)"/.exec(tag)?.[1] ?? null;
  for (const tag of TAGS.filter((entry) => entry.startsWith('<run-blocked-reason'))) {
    const id = /\sid="([^"]+)"/.exec(tag)?.[1];
    const controls = TAGS.filter((entry) => (/\saria-describedby="([^"]*)"/.exec(entry)?.[1] || '').split(/\s+/).includes(id));
    assert.ok(controls.length > 0, `a control names ${id}`);
    for (const control of controls) {
      assert.ok(conditionOf(tag) === null || conditionOf(tag) === conditionOf(control),
        `<run-blocked-reason id="${id}"> is rendered whenever its control is: ${tag} / ${control.replace(/\s+/g, ' ')}`);
    }
  }
});

test('the alignment owner\'s Align, Review and Apply read only the owner, which words the run\'s reason', () => {
  const drawer = TAGS.filter((tag) => /@click="startSimilarityAlignmentFromDrawer\(/.test(tag));
  assert.equal(drawer.length, 2);
  for (const tag of drawer) {
    assert.match(tag, /\s:disabled="Boolean\(similarityAlignmentDrawer\.disabledReason\)"/, tag);
    assert.match(tag, /\saria-describedby="similarity-alignment-drawer-reason"/, tag);
  }
  const apply = TAGS.filter((tag) => /@click="applySimilarityAlignmentDraft"/.test(tag));
  assert.equal(apply.length, 1);
  assert.match(apply[0], /\s:disabled="!canApplySimilarityAlignment"/);
  assert.match(apply[0], /\saria-describedby="similarity-alignment-apply-reason"/);
  const owner = readFileSync(new URL('../../gbdraw/web/js/app/similarity-alignment.js', import.meta.url), 'utf8');
  assert.equal(owner.split(RUN_AVAILABILITY).length - 1, 1, 'the owner reads the run availability once (actionBlock)');
});
