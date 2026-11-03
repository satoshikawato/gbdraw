import assert from 'node:assert/strict';
import { spawnSync } from 'node:child_process';
import test from 'node:test';

// Python decides whether a Session 46-48 used the mode it did not show from
// each mode's default settings (gbdraw/data/mode-slice-defaults.json), which
// the Web owns: a new drawing of `createDrawingState`. A change to a default
// or to the settings registry must regenerate the file.
test('the packaged mode slice defaults are the Web defaults', () => {
  const result = spawnSync(
    process.execPath,
    ['tools/generate_mode_slice_defaults.mjs', '--check'],
    { cwd: process.cwd(), encoding: 'utf8' }
  );
  assert.equal(result.status, 0, `${result.stdout}${result.stderr}`);
});
