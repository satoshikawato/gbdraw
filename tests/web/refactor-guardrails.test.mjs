import assert from 'node:assert/strict';
import test from 'node:test';
import {
  existsSync,
  readFileSync,
  readdirSync,
  statSync
} from 'node:fs';
import { join, relative, resolve } from 'node:path';

const root = resolve(process.env.GBDRAW_REPO || process.cwd());
const webJsRoot = join(root, 'gbdraw', 'web', 'js');

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
const previewRuntimePath = join(webJsRoot, 'app', 'preview-runtime.js');
const agentsPath = join(root, 'AGENTS.md');
const skillPath = join(
  root,
  '.agents',
  'skills',
  'refactor-gbdraw-web-safely',
  'SKILL.md'
);

test('the behavior-preserving Web refactor Skill is discoverable from AGENTS.md', () => {
  assert.equal(existsSync(skillPath), true, 'missing refactor-gbdraw-web-safely Skill');
  const agents = read(agentsPath);
  assert.match(agents, /\.agents\/skills\/refactor-gbdraw-web-safely\/SKILL\.md/);
  assert.match(agents, /characterization tests/i);
  assert.match(read(skillPath), /The new implementation is not an independent oracle/i);
});

test('fresh-SVG projection does not defensively JSON-clone large owners', () => {
  assert.equal(existsSync(projectionPath), true, 'missing editor SVG projection owner');
  const source = read(projectionPath);
  const largeOwners = [
    'features',
    'featureCatalog',
    'extractedFeatures',
    'biologicalFeatures',
    'orthogroups',
    'results',
    'svg'
  ];
  for (const owner of largeOwners) {
    assert.doesNotMatch(
      source,
      new RegExp(`cloneJson(?:Value|Data)\\s*\\(\\s*input\\.${owner}\\b`),
      `${owner} must be borrowed read-only or reduced to a compact index`
    );
    assert.doesNotMatch(
      source,
      new RegExp(`JSON\\.stringify\\s*\\(\\s*input\\.${owner}\\b`),
      `${owner} must not be serialized by the projection boundary`
    );
  }
});

test('production live-edit code has no placeholder mutate callback', () => {
  for (const [owner, source] of sourceByRelativePath) {
    assert.doesNotMatch(
      source,
      /(?:mutate\s*:|commit[A-Za-z]*Mutation\s*\([^,\n]+,)\s*(?:\(\s*\)|\([^)]*\))\s*=>\s*true\b/,
      `${owner} mutates outside the commit owner and reports a placeholder change`
    );
  }
});

test('direct editor modules do not serialize or replace Result content themselves', () => {
  const directOwners = [...sourceByRelativePath.entries()].filter(([owner]) => (
    owner.startsWith('app/feature-editor/')
    || owner.startsWith('app/legend/')
    || owner.startsWith('app/legend-layout/')
  ));
  for (const [owner, source] of directOwners) {
    assert.doesNotMatch(source, /\bserializeCleanSvg\b/, owner);
    assert.doesNotMatch(
      source,
      /(?:state\.)?results\.value\[[^\]]+\]\s*=/,
      owner
    );
  }
});

test('manual visibility rules participate in effective fresh-SVG replay', () => {
  const source = read(projectionPath);
  assert.match(
    source,
    /resolveEffectiveFeatureVisibility|effectiveFeatureVisibility|featureVisibilityManualRules/,
    'fresh-SVG projection must not replay only direct visibility overrides'
  );
});

test('mounted Result ownership does not cross await without lease revalidation', () => {
  const source = read(previewRuntimePath);
  if (!/await\s+action\s*\(/.test(source)) return;
  assert.match(
    source,
    /lease|targetToken|revalidate|validateCommitTarget|assertCurrentTarget/i,
    'an async mounted-edit batch must revalidate an explicit Result/SVG lease'
  );
  assert.doesNotMatch(
    source,
    /const\s+target\s*=\s*resolveActiveCommitTarget\(\)[\s\S]{0,1600}await\s+action\s*\(\)[\s\S]{0,1200}flushCommitTarget\s*\(\s*target\s*\)/,
    'a target captured before await must not be flushed without revalidation'
  );
});

test('test evidence does not use disconnected zero-valued execution counters', () => {
  const candidate = join(root, 'tests', 'web', 'live-svg-edit-commit-matrix.test.mjs');
  if (!existsSync(candidate)) return;
  const source = read(candidate);
  assert.doesNotMatch(
    source,
    /const\s+(?:execution|metrics|counters)\s*=\s*\{[\s\S]{0,600}workerConstruction\s*:\s*0[\s\S]{0,600}pythonCalls\s*:\s*0/,
    'Worker/Python evidence must be connected to production instrumentation'
  );
});

test('empty editor projection preserves renderer-owned label visibility', async () => {
  const helperPath = join(root, 'tests', 'web', 'helpers', 'editor-svg-fixture.mjs');
  if (!existsSync(helperPath)) return;

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
