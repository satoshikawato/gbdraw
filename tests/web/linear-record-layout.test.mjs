import assert from 'node:assert/strict';
import { cp, mkdtemp, writeFile } from 'node:fs/promises';
import { tmpdir } from 'node:os';
import { join } from 'node:path';
import test from 'node:test';
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

const freshSequences = () => [{ uid: 'a' }, { uid: 'b' }, { uid: 'c' }];

test('shared rows remain UID-keyed while records are added and removed', () => {
  const sequences = freshSequences();
  const originalEntries = [{ uid: 'a', row: 1 }, { uid: 'b', row: 1 }];
  const layout = reconcileLinearRecordLayout(sequences, originalEntries);
  assert.deepEqual(layout, [
    { uid: 'a', row: 1 },
    { uid: 'b', row: 1 },
    { uid: 'c', row: 3 }
  ]);
  assert.deepEqual(originalEntries, [{ uid: 'a', row: 1 }, { uid: 'b', row: 1 }]);

  setLinearRecordRow(layout, 'c', 1);
  assert.deepEqual(linearRecordPositionTokens(sequences, layout), [
    '#1@1', '#2@1', '#3@1'
  ]);
  const withNewRecord = reconcileLinearRecordLayout(
    [...sequences, { uid: 'd' }],
    layout
  );
  assert.deepEqual(withNewRecord, [
    { uid: 'a', row: 1 },
    { uid: 'b', row: 1 },
    { uid: 'c', row: 1 },
    { uid: 'd', row: 4 }
  ]);
  assert.deepEqual(
    reconcileLinearRecordLayout([sequences[0], sequences[2]], withNewRecord),
    [{ uid: 'a', row: 1 }, { uid: 'c', row: 1 }]
  );
});

test('an invalid row preserves only the edited UID row', () => {
  const layout = [
    { uid: 'a', row: 1 },
    { uid: 'b', row: 2 },
    { uid: 'c', row: 1 }
  ];
  setLinearRecordRow(layout, 'b', 0);
  assert.deepEqual(layout, [
    { uid: 'a', row: 1 },
    { uid: 'b', row: 2 },
    { uid: 'c', row: 1 }
  ]);
});

test('moving within a row changes card order without changing membership', () => {
  const sequences = freshSequences();
  const layout = sequences.map(({ uid }) => ({ uid, row: 1 }));
  const moved = moveLinearRecordInRow(sequences, layout, 'b', -1);
  assert.deepEqual(sequences.map((sequence) => sequence.uid), ['b', 'a', 'c']);
  assert.deepEqual(moved, [
    { uid: 'b', row: 1 },
    { uid: 'a', row: 1 },
    { uid: 'c', row: 1 }
  ]);
});
