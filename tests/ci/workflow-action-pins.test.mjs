import assert from 'node:assert/strict';
import { readdirSync, readFileSync } from 'node:fs';
import test from 'node:test';

const dir = new URL('../../.github/workflows/', import.meta.url);
const actions = new URL('../../.github/actions/', import.meta.url);

// Every third-party action is pinned to a commit and carries its full version,
// so a moved tag cannot change what runs. Local actions (./path) are exempt.
const PINNED = /^\s*(?:- )?uses: [\w.-]+\/[\w./-]+@[0-9a-f]{40} # v\d+\.\d+\.\d+$/;
const LOCAL = /^\s*(?:- )?uses: \.\//;

// Workflows and the local composite actions they run (.github/actions/<name>/action.yml).
const sources = () => [
  ...readdirSync(dir).filter(name => /\.ya?ml$/.test(name)).map(name => new URL(name, dir)),
  ...readdirSync(actions, { withFileTypes: true }).filter(entry => entry.isDirectory())
    .map(entry => new URL(`${entry.name}/action.yml`, actions))
];

test('every workflow and local action pins each action to a commit with a full version comment', () => {
  const files = sources();
  assert.ok(files.some(file => file.pathname.endsWith('/managed-losat/action.yml')));
  for (const file of files) {
    const name = file.pathname.split('/.github/')[1];
    const lines = readFileSync(file, 'utf8').split('\n');
    lines.forEach((line, index) => {
      if (!/^\s*(?:- )?uses:/.test(line) || LOCAL.test(line)) return;
      assert.match(line, PINNED, `${name}:${index + 1}`);
    });
  }
});
