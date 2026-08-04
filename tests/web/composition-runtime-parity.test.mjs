import assert from 'node:assert/strict';
import { spawnSync } from 'node:child_process';
import { fileURLToPath } from 'node:url';

import {
  COMPOSITION_METADATA_ATTRIBUTE,
  COMPOSITION_SCHEMA_ATTRIBUTE,
  parseCompositionMetadata,
  replanCompositionMetadata
} from '../../gbdraw/web/js/app/legend-layout/composition-actions.js';

const repoRoot = fileURLToPath(new URL('../../', import.meta.url));
const oracle = spawnSync(
  process.env.PYTHON || 'python',
  ['tests/web/composition-python-oracle.py'],
  { cwd: repoRoot, encoding: 'utf8' }
);
assert.equal(
  oracle.status,
  0,
  `Python composition oracle failed:\n${oracle.stderr || oracle.stdout}`
);
const cases = JSON.parse(oracle.stdout);

const assertEquivalent = (actual, expected, path = 'plan') => {
  if (typeof expected === 'number') {
    assert.equal(typeof actual, 'number', `${path} is not numeric`);
    assert.ok(
      Math.abs(actual - expected) <= 1e-9,
      `${path}: expected ${expected}, received ${actual}`
    );
    return;
  }
  if (Array.isArray(expected)) {
    assert.ok(Array.isArray(actual), `${path} is not an array`);
    assert.equal(actual.length, expected.length, `${path} length differs`);
    expected.forEach((value, index) => assertEquivalent(actual[index], value, `${path}[${index}]`));
    return;
  }
  if (expected && typeof expected === 'object') {
    assert.deepEqual(Object.keys(actual).sort(), Object.keys(expected).sort(), `${path} keys differ`);
    Object.entries(expected).forEach(([key, value]) => {
      assertEquivalent(actual[key], value, `${path}.${key}`);
    });
    return;
  }
  assert.equal(actual, expected, `${path} differs`);
};

for (const parityCase of cases) {
  const attributes = new Map([
    [COMPOSITION_SCHEMA_ATTRIBUTE, parityCase.schema],
    [COMPOSITION_METADATA_ATTRIBUTE, JSON.stringify(parityCase.metadata)]
  ]);
  const svg = {
    getAttribute: (name) => attributes.has(name) ? attributes.get(name) : null
  };
  const metadata = parseCompositionMetadata(svg);
  const replayed = replanCompositionMetadata(metadata);
  assertEquivalent(replayed, parityCase.expected, parityCase.name);
}

assert.deepEqual(
  cases.map((parityCase) => parityCase.name),
  ['dock', 'overlay', 'title_stack', 'no_legend', 'canvas_growth']
);

console.log('Python and JavaScript composition runtime parity tests passed');
