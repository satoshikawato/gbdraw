import assert from 'node:assert/strict';
import { readFile } from 'node:fs/promises';
import { join } from 'node:path';
import {
  classifyOptionalNumber,
  classifyOptionalPositiveNumber,
  projectOptionalNumber
} from '../../gbdraw/web/js/utils/optional-positive-number.js';
import { normalizeUserFacingError } from '../../gbdraw/web/js/services/error-normalization.js';

const repoRoot = process.cwd();
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

// R7: the projection passes every finite number literally (Python owns the
// range) and rejects only what JSON cannot carry as a number.
for (const [value, expected] of [
  [null, null], [undefined, null], ['', null], ['  ', null], ['Auto', null], ['none', null], ['NULL', null],
  [0, 0], [-5, -5], [12.5, 12.5], ['-100', -100], [' 500.5 ', 500.5], ['1e3', 1000]
]) {
  assert.equal(classifyOptionalNumber(value).status === 'invalid', false, JSON.stringify(value));
  assert.equal(projectOptionalNumber(value, { field: 'window' }), expected, JSON.stringify(value));
}
for (const value of [Number.NaN, Infinity, 'NaN', '1e-50x', 'abc', true, {}, []]) {
  assert.throws(() => projectOptionalNumber(value, { field: 'window' }), (error) => {
    const model = normalizeUserFacingError(error);
    assert.equal(model.code, 'INPUT_INVALID');
    assert.deepEqual(model.context, { field: 'window', reason: 'FINITE' });
    return true;
  }, JSON.stringify(value));
}
assert.throws(() => projectOptionalNumber('x', { configPath: 'objects.scale.interval' }),
  { code: 'INPUT_INVALID', context: { configPath: 'objects.scale.interval', reason: 'FINITE' } });
