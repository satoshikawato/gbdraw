import assert from 'node:assert/strict';
import { cp, mkdtemp, writeFile } from 'node:fs/promises';
import { tmpdir } from 'node:os';
import { join } from 'node:path';
import { pathToFileURL } from 'node:url';

const repoRoot = process.cwd();
const tempRoot = await mkdtemp(join(tmpdir(), 'gbdraw-linear-comparisons-'));
await cp(
  join(repoRoot, 'gbdraw', 'web', 'js', 'app', 'linear-comparisons.js'),
  join(tempRoot, 'linear-comparisons.js')
);
await writeFile(join(tempRoot, 'package.json'), '{"type":"module"}', 'utf8');

const {
  addLinearComparison,
  adjacentRowPairs,
  reconcileLinearComparisons,
  validateLinearComparisons
} = await import(pathToFileURL(join(tempRoot, 'linear-comparisons.js')));

const sequences = [{ uid: 'a' }, { uid: 'b' }, { uid: 'c' }, { uid: 'd' }];
const layout = [
  { uid: 'a', row: 1 }, { uid: 'b', row: 1 },
  { uid: 'c', row: 2 }, { uid: 'd', row: 2 }
];
assert.deepEqual(adjacentRowPairs(sequences, layout), [['a', 'c'], ['b', 'd']]);
assert.deepEqual(adjacentRowPairs(sequences, layout, true), [
  ['a', 'c'], ['a', 'd'], ['b', 'c'], ['b', 'd']
]);

const comparisons = [];
assert.equal(addLinearComparison(comparisons, 'a', 'c'), true);
assert.equal(addLinearComparison(comparisons, 'a', 'c'), false);
comparisons[0].file = { name: 'a-c.tsv' };
assert.equal(validateLinearComparisons(sequences, layout, comparisons), '');
comparisons[0].source = 'losat';
comparisons[0].file = null;
assert.equal(validateLinearComparisons(sequences, layout, comparisons), '');
comparisons[0].subjectUid = 'b';
assert.match(validateLinearComparisons(sequences, layout, comparisons), /adjacent rows/);
assert.deepEqual(reconcileLinearComparisons(sequences.slice(1), comparisons), []);
