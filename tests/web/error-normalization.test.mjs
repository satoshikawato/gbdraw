import assert from 'node:assert/strict';
import {
  normalizeUserFacingError,
  safeErrorText
} from '../../gbdraw/web/js/services/error-normalization.js';
import {
  deserializeWorkerError
} from '../../gbdraw/web/js/services/diagram-generation.js';

const structured = normalizeUserFacingError({
  type: 'ValidationError',
  message: { summary: 'Invalid annotation', field: 'set_id' },
  stderr: 'Traceback (most recent call last):\n  File "request_render.py", line 1\nValueError: bad row',
  stdout: 'validation stopped',
  traceback: 'Traceback: private diagnostic',
  notes: [{ message: 'temporary workspace cleanup failed' }]
});
assert.equal(structured.summary, 'ValidationError: Invalid annotation');
assert.deepEqual(structured.details, [
  { label: 'STDOUT', text: 'validation stopped' },
  { label: 'Cleanup note', text: 'temporary workspace cleanup failed' }
]);

const arbitrary = normalizeUserFacingError({
  code: 'bad_input',
  sequence: 'ACGT'.repeat(5000),
  content: 'private file contents',
  field: 'track_slots'
});
assert.equal(arbitrary.summary, 'code: bad_input\nfield: track_slots');
assert.doesNotMatch(JSON.stringify(arbitrary), /ACGT|private file contents|\[object Object\]|Traceback/);

const workerError = deserializeWorkerError({
  name: 'ValidationError',
  message: 'Annotation target is missing.',
  details: [{ label: 'Annotation row', text: 'Row 2 references an unknown set.' }]
});
assert.deepEqual(workerError.details, [
  { label: 'Annotation row', text: 'Row 2 references an unknown set.' }
]);
assert.doesNotMatch(JSON.stringify(normalizeUserFacingError(workerError)), /\[object Object\]/);

const nestedJson = normalizeUserFacingError(JSON.stringify({
  summary: 'Could not render diagram',
  details: [
    { label: 'Details', text: 'Annotation set "missing" does not exist.' },
    { label: 'Traceback', text: 'Traceback: hidden' }
  ]
}));
assert.deepEqual(nestedJson, {
  summary: 'Could not render diagram',
  details: [{ label: 'Details', text: 'Annotation set "missing" does not exist.' }]
});

const bounded = safeErrorText('x'.repeat(5000), { limit: 80 });
assert.ok(bounded.length <= 80);
assert.match(bounded, /details truncated/);
