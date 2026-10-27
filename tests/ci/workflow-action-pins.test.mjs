import assert from 'node:assert/strict';
import { readdirSync, readFileSync } from 'node:fs';
import test from 'node:test';

const dir = new URL('../../.github/workflows/', import.meta.url);

// Every third-party action is pinned to a commit and carries its full version,
// so a moved tag cannot change what runs. Local actions (./path) are exempt.
const PINNED = /^\s*(?:- )?uses: [\w.-]+\/[\w./-]+@[0-9a-f]{40} # v\d+\.\d+\.\d+$/;
const LOCAL = /^\s*(?:- )?uses: \.\//;

test('every workflow pins each action to a commit with a full version comment', () => {
  const files = readdirSync(dir).filter(name => /\.ya?ml$/.test(name));
  assert.ok(files.length > 0);
  for (const name of files) {
    const lines = readFileSync(new URL(name, dir), 'utf8').split('\n');
    lines.forEach((line, index) => {
      if (!/^\s*(?:- )?uses:/.test(line) || LOCAL.test(line)) return;
      assert.match(line, PINNED, `${name}:${index + 1}`);
    });
  }
});
