import assert from 'node:assert/strict';
import { test } from 'node:test';
import { describeGenerationApplication } from '../../gbdraw/web/js/app/generation-status.js';

const facts = (status, extra = {}) => ({ status, differences: [], unknown: [], operations: {
  generating: false, cancelRequested: false, liveApplying: false, liveError: null
}, ...extra });

test('application feedback keeps uncertain/invalid/unrendered settings distinct from Applied', () => {
  for (const status of ['pending', 'unknown', 'invalid', 'ungenerated']) {
    const feedback = describeGenerationApplication(facts(status));
    assert.notEqual(feedback.label, 'Applied');
    assert.ok(feedback.message);
  }
  assert.equal(describeGenerationApplication(facts('clean')).label, 'Applied');
});

test('Pending, uncertainty, live failure and generation operations remain independent', () => {
  const input = facts('pending', { differences: ['$.scale'], unknown: ['$live-render'],
    operations: { generating: true, cancelRequested: true, liveApplying: false, liveError: 'failure' } });
  const original = structuredClone(input);
  const feedback = describeGenerationApplication(input);
  assert.equal(feedback.label, 'Pending');
  assert.match(feedback.unknownMessage, /unknown/);
  assert.match(feedback.liveMessage, /failed.*direct edits already applied are kept/i);
  assert.match(feedback.generationMessage, /Canceling/);
  assert.deepEqual(input, original);
  input.operations.liveApplying = true;
  assert.match(describeGenerationApplication(input).liveMessage, /applying/);
});

test('screen-reader feedback ignores path/count and validation-detail keystrokes', () => {
  const initial = facts('pending', { differences: ['$.scale'], unknown: ['$.table'] });
  const changed = { ...initial, differences: ['$.scale', '$.crop'], unknown: ['$.source'] };
  assert.equal(describeGenerationApplication(initial).announcement,
    describeGenerationApplication(changed).announcement);
  const invalid = describeGenerationApplication(facts('invalid', { error: 'Input A is invalid' }));
  const correction = describeGenerationApplication(facts('invalid', { error: 'Input B is invalid' }));
  assert.notEqual(invalid.error, correction.error);
  assert.equal(invalid.announcement, correction.announcement);
  assert.notEqual(describeGenerationApplication(facts('clean')).announcement,
    describeGenerationApplication(initial).announcement);
});
