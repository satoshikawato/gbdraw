import assert from 'node:assert/strict';
import { mkdtemp, readFile, writeFile } from 'node:fs/promises';
import { tmpdir } from 'node:os';
import { pathToFileURL } from 'node:url';
import { join } from 'node:path';

const repoRoot = process.cwd();
const tempRoot = await mkdtemp(join(tmpdir(), 'gbdraw-optional-positive-number-'));
await writeFile(join(tempRoot, 'package.json'), '{"type":"module"}\n', 'utf8');
await writeFile(
  join(tempRoot, 'optional-positive-number.js'),
  await readFile(join(repoRoot, 'gbdraw', 'web', 'js', 'utils', 'optional-positive-number.js'), 'utf8'),
  'utf8'
);
const { classifyOptionalPositiveNumber } = await import(pathToFileURL(
  join(tempRoot, 'optional-positive-number.js')
));
const cases = JSON.parse(await readFile(
  join(repoRoot, 'tests', 'fixtures', 'optional-positive-number.json'),
  'utf8'
));

for (const testCase of cases) {
  const actual = classifyOptionalPositiveNumber(testCase.value);
  assert.equal(actual.status, testCase.status, JSON.stringify(testCase.value));
  if (testCase.status !== 'invalid') {
    assert.equal(actual.value, testCase.normalized, JSON.stringify(testCase.value));
  }
}

for (const value of [undefined, Number.NaN, Number.POSITIVE_INFINITY, Number.NEGATIVE_INFINITY]) {
  const expected = value === undefined ? 'auto' : 'invalid';
  assert.equal(classifyOptionalPositiveNumber(value).status, expected);
}
