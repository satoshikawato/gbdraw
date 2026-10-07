// OV-160: the Web's CSS named-color table is Python's
// (`gbdraw.io.colors._COLOR_NAME_MAP`, pinned to CSS by tests/test_named_colors.py),
// so Load and the Session split resolve a color name to the hex Python writes.
import assert from 'node:assert/strict';
import { execFileSync } from 'node:child_process';
import { fileURLToPath } from 'node:url';
import { test } from 'node:test';

const { CSS_NAMED_COLORS, namedColorHex } = await import('../../gbdraw/web/js/utils/named-colors.js');
const { resolveColorToHex } = await import('../../gbdraw/web/js/utils/color-utils.js');

const REPO = fileURLToPath(new URL('../../', import.meta.url));

test('the Web named-color table equals the Python table', () => {
  const python = JSON.parse(execFileSync('python', ['-c', `
import json
from gbdraw.io.colors import _COLOR_NAME_MAP
print(json.dumps(_COLOR_NAME_MAP))
`], { encoding: 'utf8', cwd: REPO, env: { ...process.env, PYTHONPATH: REPO } }));
  assert.deepEqual({ ...CSS_NAMED_COLORS }, python);
  assert.equal(Object.keys(CSS_NAMED_COLORS).length, 148);
});

test('a color name resolves without a browser canvas, in any case; other values pass through', () => {
  assert.equal(globalThis.document, undefined);
  assert.equal(namedColorHex('Gray'), '#808080');
  assert.equal(resolveColorToHex(' gray '), '#808080');
  assert.equal(resolveColorToHex('REBECCAPURPLE'), '#663399');
  assert.equal(resolveColorToHex('#abc'), '#abc');
  assert.equal(namedColorHex('notacolor'), null);
  assert.equal(resolveColorToHex('notacolor'), 'notacolor');
  assert.equal(namedColorHex('constructor'), null);
});
