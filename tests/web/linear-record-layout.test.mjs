import assert from 'node:assert/strict';
import { cp, mkdtemp, writeFile } from 'node:fs/promises';
import { tmpdir } from 'node:os';
import { join } from 'node:path';
import { pathToFileURL } from 'node:url';

const repoRoot = process.cwd();
const tempRoot = await mkdtemp(join(tmpdir(), 'gbdraw-linear-layout-'));
await cp(
  join(repoRoot, 'gbdraw', 'web', 'js', 'app', 'linear-record-layout.js'),
  join(tempRoot, 'linear-record-layout.js')
);
await writeFile(join(tempRoot, 'package.json'), '{"type":"module"}', 'utf8');

const {
  reconcileLinearRecordLayout,
  linearRecordPositionTokens,
  moveLinearRecordInRow,
  setLinearRecordRow
} = await import(pathToFileURL(join(tempRoot, 'linear-record-layout.js')));

const sequences = [{ uid: 'a' }, { uid: 'b' }, { uid: 'c' }];
const layout = reconcileLinearRecordLayout(sequences, [{ uid: 'a', row: 1 }, { uid: 'b', row: 1 }]);
assert.deepEqual(layout, [{ uid: 'a', row: 1 }, { uid: 'b', row: 1 }, { uid: 'c', row: 3 }]);
setLinearRecordRow(layout, 'c', 2);
assert.deepEqual(linearRecordPositionTokens(sequences, layout), ['#1@1', '#2@1', '#3@2']);

const moved = moveLinearRecordInRow(sequences, layout, 'b', -1);
assert.deepEqual(sequences.map((sequence) => sequence.uid), ['b', 'a', 'c']);
assert.deepEqual(moved.map((entry) => entry.uid), ['b', 'a', 'c']);
