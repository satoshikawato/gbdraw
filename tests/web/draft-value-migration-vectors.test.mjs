import assert from 'node:assert/strict';
import { spawnSync } from 'node:child_process';
import test from 'node:test';

// The Python twins of the Web Load value migrations of a Session 27-44 draft
// (gbdraw/session_io.py) read vectors that the Web functions of
// services/config.js wrote. A change to those functions must regenerate them.
test('the draft value migration vectors are the Web functions output', () => {
  const result = spawnSync(
    process.execPath,
    ['tools/generate_draft_value_migration_vectors.mjs', '--check'],
    { cwd: process.cwd(), encoding: 'utf8' }
  );
  assert.equal(result.status, 0, `${result.stdout}${result.stderr}`);
});
