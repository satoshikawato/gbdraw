import assert from 'node:assert/strict';
import { cp, mkdtemp, writeFile } from 'node:fs/promises';
import { tmpdir } from 'node:os';
import { join } from 'node:path';
import { pathToFileURL } from 'node:url';

const repoRoot = process.cwd();
const tempRoot = await mkdtemp(join(tmpdir(), 'gbdraw-linear-layout-'));
await cp(
  join(repoRoot, 'gbdraw', 'web', 'js', 'services', 'linear-record-layout.js'),
  join(tempRoot, 'linear-record-layout.js')
);
await writeFile(join(tempRoot, 'package.json'), '{"type":"module"}', 'utf8');

const {
  linearRecordLayoutHasSharedRow,
  recordTranslationsFor,
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
    { uid: 'a-1', row: 1 }, { uid: 'a-2', row: 3 },
    { uid: 'b-1', row: 2 }, { uid: 'b-2', row: 2 }, { uid: 'b-3', row: 2 }
  ],
  [
    { uid: 'a-1', row: 1 }, { uid: 'a-2', row: 2 },
    { uid: 'b-1', row: 2 }, { uid: 'b-2', row: 3 }, { uid: 'b-3', row: 3 }
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

// A File may span several consecutive rows. A move exchanges whole row
// blocks and keeps each record's row relative to its own File.
const blockGroups = [
  { records: [{ sequence: { uid: 'one' } }] },
  { records: [{ sequence: { uid: 'two-b' } }, { sequence: { uid: 'two-a' } }, { sequence: { uid: 'two-c' } }] },
  { records: [{ sequence: { uid: 'three' } }] }
];
const blockRows = [
  { uid: 'one', row: 1 },
  { uid: 'two-b', row: 3 }, { uid: 'two-a', row: 2 }, { uid: 'two-c', row: 3 },
  { uid: 'three', row: 4 }
];
assert.deepEqual(planLinearSourceRowMove({
  sourceGroups: blockGroups, entries: blockRows, sourceIndex: 1, direction: -1
}), {
  allowed: true,
  reason: '',
  rows: [
    { uid: 'two-b', row: 2 }, { uid: 'two-a', row: 1 }, { uid: 'two-c', row: 2 },
    { uid: 'one', row: 3 },
    { uid: 'three', row: 4 }
  ]
}, 'a two-row File moves up as one block');
assert.deepEqual(planLinearSourceRowMove({
  sourceGroups: blockGroups, entries: blockRows, sourceIndex: 1, direction: 1
}), {
  allowed: true,
  reason: '',
  rows: [
    { uid: 'one', row: 1 },
    { uid: 'three', row: 2 },
    { uid: 'two-b', row: 4 }, { uid: 'two-a', row: 3 }, { uid: 'two-c', row: 4 }
  ]
}, 'a two-row File moves down as one block');

// A CLI Session draws each record on its own row: File 1 (one record) on
// row 1 and File 2 (two records) on rows 2 and 3 form consecutive blocks.
const cliGroups = [
  { records: [{ sequence: { uid: 'single' } }] },
  { records: [{ sequence: { uid: 'pair-1' } }, { sequence: { uid: 'pair-2' } }] }
];
assert.deepEqual(planLinearSourceRowMove({
  sourceGroups: cliGroups, entries: [], sourceIndex: 1, direction: -1
}), {
  allowed: true,
  reason: '',
  rows: [{ uid: 'pair-1', row: 1 }, { uid: 'pair-2', row: 2 }, { uid: 'single', row: 3 }]
}, 'per-record CLI rows move File blocks');

// OV-2291: the translations of a record set, by record key.
{
  const stored = [
    { recordKey: 'record-1', x: 10, y: 0 },
    { recordKey: 'record-2', x: 0, y: -5 },
    { recordKey: 'record-3', x: 7, y: 7 }
  ];
  // A removed record's translation is left out; the others keep their values, in the drawn order.
  assert.deepEqual(recordTranslationsFor(stored, ['record-3', 'record-1']), [
    { recordKey: 'record-3', x: 7, y: 7 },
    { recordKey: 'record-1', x: 10, y: 0 }
  ]);
  // An added record starts at (0, 0).
  assert.deepEqual(recordTranslationsFor(stored, ['record-1', 'record-2', 'record-3', 'record-4']), [
    { recordKey: 'record-1', x: 10, y: 0 },
    { recordKey: 'record-2', x: 0, y: -5 },
    { recordKey: 'record-3', x: 7, y: 7 },
    { recordKey: 'record-4', x: 0, y: 0 }
  ]);
  // Nothing stored: every record at (0, 0); no records: nothing.
  assert.deepEqual(recordTranslationsFor([], ['a', 'b']), [{ recordKey: 'a', x: 0, y: 0 }, { recordKey: 'b', x: 0, y: 0 }]);
  assert.deepEqual(recordTranslationsFor(stored, []), []);
  // A stored entry of a drawn record passes through as it is, also a duplicate or a
  // value that is not a number: the request's validation refuses those, not this function.
  assert.deepEqual(
    recordTranslationsFor([{ recordKey: 'a', x: '3', y: Infinity }, { recordKey: 'a', x: 1, y: 1 }, null], ['a']),
    [{ recordKey: 'a', x: '3', y: Infinity }, { recordKey: 'a', x: 1, y: 1 }]
  );
}
