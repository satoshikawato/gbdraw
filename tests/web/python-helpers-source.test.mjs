import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import test from 'node:test';

const HELPERS_URL = new URL('../../gbdraw/web/js/app/python-helpers.js', import.meta.url);

// The diagram Worker runs PYTHON_HELPERS, the value of a JavaScript template literal, while the Python tests read the
// text between its backticks. A backtick, `${` or backslash escape in the Python changes or breaks that value (a
// docstring with ``name`` made the module throw and no Worker start, OV-3493), so the value is the text as written.
test('the Worker runs the Python helper source as written', async () => {
  const { PYTHON_HELPERS } = await import(HELPERS_URL.href);
  const source = readFileSync(HELPERS_URL, 'utf8');
  const start = 'export const PYTHON_HELPERS = `';
  const end = '`;\n';
  assert.equal(source.indexOf(start), source.lastIndexOf(start));
  const written = source.slice(source.indexOf(start) + start.length, source.lastIndexOf(end));
  assert.equal(typeof PYTHON_HELPERS, 'string');
  assert.ok(PYTHON_HELPERS === written, 'PYTHON_HELPERS differs from the Python text written in python-helpers.js');
});
