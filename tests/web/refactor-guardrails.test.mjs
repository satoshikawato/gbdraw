import assert from 'node:assert/strict';
import test from 'node:test';
import {
  existsSync,
  readFileSync,
  readdirSync,
  statSync
} from 'node:fs';
import { join, relative, resolve } from 'node:path';

import { analyzeWebRefactorSource } from './helpers/refactor-guard-ast.mjs';

const root = resolve(process.env.GBDRAW_REPO || process.cwd());
const webJsRoot = join(root, 'gbdraw', 'web', 'js');
const fixtureRoot = join(root, 'tests', 'web', 'fixtures', 'refactor-guardrails');
const read = (path) => readFileSync(path, 'utf8');

const walk = (directory) => {
  if (!existsSync(directory)) return [];
  return readdirSync(directory)
    .flatMap((name) => {
      const path = join(directory, name);
      return statSync(path).isDirectory() ? walk(path) : [path];
    });
};

const productionJs = walk(webJsRoot).filter((path) => path.endsWith('.js'));
const sourceByRelativePath = new Map(
  productionJs.map((path) => [relative(webJsRoot, path).replaceAll('\\', '/'), read(path)])
);
const projectionPath = join(webJsRoot, 'app', 'editor-svg-projection.js');
const resultsPath = join(webJsRoot, 'app', 'results.js');

const assertNoViolations = (owner, source, options) => {
  const violations = analyzeWebRefactorSource(source, options);
  assert.deepEqual(
    violations,
    [],
    `${owner}: ${violations.map(({ code, message }) => `${code}: ${message}`).join('; ')}`
  );
};

test('fresh-SVG projection borrows large owners and emits no synthetic zero evidence', () => {
  assertNoViolations('app/editor-svg-projection.js', read(projectionPath), {
    largeOwnerCopies: true,
    syntheticMetrics: true
  });
});

test('production live-edit commits have no placeholder mutation callbacks', () => {
  for (const [owner, source] of sourceByRelativePath) {
    assertNoViolations(owner, source, { placeholderCommits: true });
  }
});

test('direct editor modules do not serialize or replace Result content themselves', () => {
  const directOwners = [...sourceByRelativePath.entries()].filter(([owner]) => (
    owner.startsWith('app/feature-editor/')
    || owner.startsWith('app/legend/')
    || owner.startsWith('app/legend-layout/')
  ));
  for (const [owner, source] of directOwners) {
    assertNoViolations(owner, source, { directResultOwnership: true });
  }
});

test('Label and SVG style commit paths use the supplied mutation journal', () => {
  for (const owner of ['app/feature-editor/label-actions.js', 'app/svg-styles.js']) {
    assertNoViolations(owner, sourceByRelativePath.get(owner), { journalDomWrites: true });
  }
});

test('definition helper settlements use a pre-await token and revalidate every await', () => {
  assertNoViolations('app/results.js', read(resultsPath), { asyncTokenProtocol: true });
});

const badFixtures = [
  {
    file: 'alias-clone.js',
    options: { largeOwnerCopies: true },
    violation: 'large-owner-clone'
  },
  {
    file: 'structured-clone.js',
    options: { largeOwnerCopies: true },
    violation: 'large-owner-structured-clone'
  },
  {
    file: 'aliased-serialization.js',
    options: { largeOwnerCopies: true },
    violation: 'large-owner-serialization'
  },
  {
    file: 'placeholder-mutation.js',
    options: { placeholderCommits: true },
    violation: 'placeholder-mutation'
  },
  {
    file: 'unjournaled-dom-mutation.js',
    options: { journalDomWrites: true },
    violation: 'unjournaled-dom-mutation'
  },
  {
    file: 'unjournaled-dom-property.js',
    options: { journalDomWrites: true },
    violation: 'unjournaled-dom-mutation'
  },
  {
    file: 'direct-result-write.js',
    options: { directResultOwnership: true },
    violation: 'direct-result-write'
  },
  {
    file: 'direct-result-serialization.js',
    options: { directResultOwnership: true },
    violation: 'direct-result-serialization'
  },
  {
    file: 'async-commit-without-token.js',
    options: { asyncTokenProtocol: true },
    violation: 'async-target-without-token'
  },
  {
    file: 'async-commit-without-revalidation.js',
    options: { asyncTokenProtocol: true },
    violation: 'async-target-without-revalidation'
  },
  {
    file: 'async-commit-without-target-token.js',
    options: { asyncTokenProtocol: true },
    violation: 'async-commit-without-token'
  },
  {
    file: 'synthetic-zero-metric.js',
    options: { syntheticMetrics: true },
    violation: 'synthetic-zero-metric'
  }
];

for (const badFixture of badFixtures) {
  test(`AST guard rejects ${badFixture.file}`, () => {
    const violations = analyzeWebRefactorSource(
      read(join(fixtureRoot, badFixture.file)),
      badFixture.options
    );
    assert.ok(
      violations.some(({ code }) => code === badFixture.violation),
      `${badFixture.file} did not trigger ${badFixture.violation}: ${JSON.stringify(violations)}`
    );
  });
}

test('empty editor projection preserves renderer-owned label visibility', async () => {
  const helperPath = join(root, 'tests', 'web', 'helpers', 'editor-svg-fixture.mjs');
  const {
    createEditorSvgProjection
  } = await import(`${projectionPath}?guard=${Date.now()}`);
  const {
    FakeSvgElement,
    appendEditableLabel
  } = await import(`${helperPath}?guard=${Date.now()}`);

  const svg = new FakeSvgElement('svg');
  const label = appendEditableLabel(svg, 'rendered-a', 'renderer label');
  label.setAttribute('display', 'none');

  createEditorSvgProjection({}).project(svg, { resetLabelState: true });

  assert.equal(
    label.getAttribute('display'),
    'none',
    'clearing editor state must preserve a renderer-hidden baseline'
  );
});
