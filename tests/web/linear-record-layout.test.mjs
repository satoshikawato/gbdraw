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
  linearRecordLayoutHasSharedRow,
  reconcileLinearRecordLayout,
  resolveEffectiveLinearRecordRows,
  linearRecordPositionTokens,
  moveLinearRecordInRow,
  planLinearSourceRowMove,
  setLinearRecordRow
} = await import(pathToFileURL(join(tempRoot, 'linear-record-layout.js')));

const sequences = [{ uid: 'a' }, { uid: 'b' }, { uid: 'c' }];
const layout = reconcileLinearRecordLayout(sequences, [{ uid: 'a', row: 1 }, { uid: 'b', row: 1 }]);
assert.deepEqual(layout, [{ uid: 'a', row: 1 }, { uid: 'b', row: 1 }, { uid: 'c', row: 3 }]);
assert.equal(linearRecordLayoutHasSharedRow(sequences, layout), true);
assert.equal(linearRecordLayoutHasSharedRow(sequences, layout, { enabled: false }), false);
assert.deepEqual(resolveEffectiveLinearRecordRows(sequences, layout, { enabled: false }), [
  { uid: 'a', row: 1 }, { uid: 'b', row: 2 }, { uid: 'c', row: 3 }
]);
const fourSequences = ['a', 'b', 'c', 'd'].map((uid) => ({ uid }));
assert.equal(linearRecordLayoutHasSharedRow(fourSequences, [
  { uid: 'a', row: 1 }, { uid: 'b', row: 2 },
  { uid: 'c', row: 3 }, { uid: 'd', row: 4 }
]), false, 'four records on four rendered rows are not shared');
assert.equal(linearRecordLayoutHasSharedRow(fourSequences, [
  { uid: 'a', row: 1 }, { uid: 'b', row: 1 },
  { uid: 'c', row: 2 }, { uid: 'd', row: 3 }
]), true, 'rows 1,1,2,3 contain a shared rendered row');
assert.equal(linearRecordLayoutHasSharedRow(fourSequences, [
  { uid: 'a', row: 1 }, { uid: 'b', row: 1 },
  { uid: 'c', row: 2 }, { uid: 'd', row: 3 }
], { enabled: false }), false, 'disabled layout ignores dormant shared rows');
setLinearRecordRow(layout, 'c', 2);
assert.deepEqual(linearRecordPositionTokens(sequences, layout), ['#1@1', '#2@1', '#3@2']);

const moved = moveLinearRecordInRow(sequences, layout, 'b', -1);
assert.deepEqual(sequences.map((sequence) => sequence.uid), ['b', 'a', 'c']);
assert.deepEqual(moved.map((entry) => entry.uid), ['b', 'a', 'c']);

const sourceGroups = [
  { records: [{ sequence: { uid: 'a-1' } }, { sequence: { uid: 'a-2' } }] },
  { records: [{ sequence: { uid: 'b-1' } }, { sequence: { uid: 'b-2' } }, { sequence: { uid: 'b-3' } }] }
];
const sourceRows = [
  { uid: 'a-1', row: 1 }, { uid: 'a-2', row: 1 },
  { uid: 'b-1', row: 3 }, { uid: 'b-2', row: 3 }, { uid: 'b-3', row: 3 }
];
const rowMove = planLinearSourceRowMove({
  sourceGroups,
  entries: sourceRows,
  sourceIndex: 1,
  direction: -1
});
assert.deepEqual(rowMove, {
  allowed: true,
  reason: '',
  rows: [
    { uid: 'b-1', row: 1 }, { uid: 'b-2', row: 1 }, { uid: 'b-3', row: 1 },
    { uid: 'a-1', row: 3 }, { uid: 'a-2', row: 3 }
  ]
});
assert.deepEqual(sourceRows, [
  { uid: 'a-1', row: 1 }, { uid: 'a-2', row: 1 },
  { uid: 'b-1', row: 3 }, { uid: 'b-2', row: 3 }, { uid: 'b-3', row: 3 }
], 'planning a File move does not mutate the current layout');

for (const customRows of [
  [
    { uid: 'a-1', row: 1 }, { uid: 'a-2', row: 2 },
    { uid: 'b-1', row: 3 }, { uid: 'b-2', row: 3 }, { uid: 'b-3', row: 3 }
  ],
  [
    { uid: 'a-1', row: 1 }, { uid: 'a-2', row: 1 },
    { uid: 'b-1', row: 1 }, { uid: 'b-2', row: 1 }, { uid: 'b-3', row: 1 }
  ]
]) {
  const blocked = planLinearSourceRowMove({
    sourceGroups,
    entries: customRows,
    sourceIndex: 1,
    direction: -1
  });
  assert.equal(blocked.allowed, false);
  assert.equal(blocked.reason, 'custom-layout');
  assert.deepEqual(blocked.rows, customRows);
}

assert.deepEqual(planLinearSourceRowMove({
  sourceGroups,
  entries: sourceRows,
  sourceIndex: 0,
  direction: -1
}), {
  allowed: false,
  reason: 'boundary',
  rows: sourceRows
}, 'a boundary move does not change the current layout');
