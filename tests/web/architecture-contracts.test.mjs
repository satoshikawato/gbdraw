import assert from 'node:assert/strict';
import {
  existsSync,
  readFileSync,
  readdirSync
} from 'node:fs';
import {
  dirname,
  extname,
  join,
  relative,
  resolve,
  sep
} from 'node:path';
import test from 'node:test';
import { fileURLToPath } from 'node:url';

const REPOSITORY_ROOT = resolve(dirname(fileURLToPath(import.meta.url)), '../..');
const JAVASCRIPT_ROOT = join(REPOSITORY_ROOT, 'gbdraw/web/js');
const APP_ENTRY = join(JAVASCRIPT_ROOT, 'app.js');

const normalizePath = (path) => path.split(sep).join('/');
const relativeModulePath = (path) => normalizePath(relative(JAVASCRIPT_ROOT, path));

const walkJavaScriptFiles = (directory) => readdirSync(directory, { withFileTypes: true })
  .flatMap((entry) => {
    const path = join(directory, entry.name);
    if (entry.isDirectory()) return walkJavaScriptFiles(path);
    return entry.isFile() && entry.name.endsWith('.js') ? [path] : [];
  });

const productionSources = new Map(
  walkJavaScriptFiles(JAVASCRIPT_ROOT).map((path) => [
    relativeModulePath(path),
    readFileSync(path, 'utf8')
  ])
);

const literalImportSpecifiers = (source) => {
  const specifiers = [];
  const patterns = [
    /(?:^|\n)\s*(?:import|export)\s+[^;]*?\s+from\s+(['"])([^'"]+)\1\s*;?/g,
    /(?:^|\n)\s*import\s*(['"])([^'"]+)\1\s*;?/g,
    /\bimport\s*\(\s*(['"])([^'"]+)\1\s*\)/g
  ];
  for (const pattern of patterns) {
    for (const match of source.matchAll(pattern)) specifiers.push(match[2]);
  }
  return specifiers;
};

const resolveRelativeImport = (ownerPath, specifier) => {
  if (!specifier.startsWith('.')) return null;
  const unresolved = resolve(dirname(ownerPath), specifier);
  const candidates = extname(unresolved)
    ? [unresolved]
    : [`${unresolved}.js`, join(unresolved, 'index.js')];
  const resolved = candidates.find((candidate) => existsSync(candidate));
  assert.ok(resolved, `Cannot resolve ${specifier} imported by ${relativeModulePath(ownerPath)}`);
  assert.ok(
    resolved.startsWith(`${JAVASCRIPT_ROOT}${sep}`),
    `First-party Web import escapes js/: ${specifier}`
  );
  return resolved;
};

const directImports = new Map(
  [...productionSources].map(([owner, source]) => {
    const ownerPath = join(JAVASCRIPT_ROOT, owner);
    const imports = literalImportSpecifiers(source)
      .map((specifier) => resolveRelativeImport(ownerPath, specifier))
      .filter(Boolean)
      .map(relativeModulePath);
    return [owner, new Set(imports)];
  })
);

const reachableModules = (entryPath) => {
  const visited = new Set();
  const pending = [relativeModulePath(entryPath)];
  while (pending.length) {
    const current = pending.pop();
    if (visited.has(current)) continue;
    visited.add(current);
    for (const dependency of directImports.get(current) || []) pending.push(dependency);
  }
  return visited;
};

const occurrenceOwners = (pattern) => new Map(
  [...productionSources]
    .map(([owner, source]) => [owner, [...source.matchAll(pattern)].length])
    .filter(([, count]) => count > 0)
);

const importersOf = (target) => [...directImports]
  .filter(([, imports]) => imports.has(target))
  .map(([owner]) => owner)
  .sort();

test('the main application import graph excludes Worker-only modules', () => {
  const reachable = reachableModules(APP_ENTRY);
  assert.ok(reachable.has('app/run-analysis.js'), 'render owner must remain reachable');
  assert.deepEqual(
    [...reachable].filter((path) => path.startsWith('workers/')),
    []
  );
});

test('Worker construction and the diagram-generation client have explicit owners', () => {
  assert.deepEqual(
    occurrenceOwners(/\bnew\s+Worker\s*\(/g),
    new Map([
      ['services/diagram-generation.js', 1],
      ['services/losat.js', 2],
      ['workers/losat-threaded-worker.js', 2]
    ])
  );
  assert.deepEqual(
    occurrenceOwners(/diagram-generation-worker\.js/g),
    new Map([['services/diagram-generation.js', 1]])
  );
  assert.deepEqual(
    occurrenceOwners(/\brunDiagramGeneration\b/g),
    new Map([
      ['app/run-analysis.js', 2],
      ['services/diagram-generation.js', 1]
    ])
  );
  assert.deepEqual(
    occurrenceOwners(/\brunFeatureExtraction\b/g),
    new Map([
      ['app/feature-metadata-extraction.js', 2],
      ['services/diagram-generation.js', 1],
      ['workers/diagram-generation-worker.js', 2]
    ])
  );
  assert.deepEqual(importersOf('services/diagram-generation.js'), [
    'app/app-setup.js',
    'app/feature-metadata-extraction.js',
    'app/run-analysis.js'
  ]);
});

test('Pyodide initialization remains on its explicit main and Worker owners', () => {
  assert.deepEqual(
    occurrenceOwners(/\bloadPyodide\b/g),
    new Map([
      ['app/pyodide.js', 1],
      ['workers/diagram-generation-worker.js', 2]
    ])
  );
  assert.deepEqual(
    occurrenceOwners(/\bloadPyodide\s*\(/g),
    new Map([
      ['app/pyodide.js', 1],
      ['workers/diagram-generation-worker.js', 1]
    ])
  );
  const reachable = reachableModules(APP_ENTRY);
  assert.deepEqual(
    [...occurrenceOwners(/\bloadPyodide\s*\(/g).keys()].filter((path) => reachable.has(path)),
    ['app/pyodide.js']
  );
});

test('production rendering crosses the canonical request and Worker boundary', () => {
  assert.deepEqual(
    occurrenceOwners(/\brunDiagramGeneration\s*\(/g),
    new Map([['app/run-analysis.js', 1]])
  );
  assert.deepEqual(
    occurrenceOwners(/\bbuildCanonicalRenderRequest\s*\(/g),
    new Map([
      ['app/run-analysis.js', 1],
      ['services/config.js', 1],
      ['services/gallery-session-migration.js', 1]
    ])
  );
  assert.deepEqual(directImports.get('app/run-analysis.js')?.has('services/session-request.js'), true);
  assert.deepEqual(directImports.get('app/run-analysis.js')?.has('services/diagram-generation.js'), true);
  assert.match(
    productionSources.get('app/run-analysis.js'),
    /runDiagramGeneration\s*\(\s*\{\s*request:\s*canonical\.renderRequest,\s*resources:\s*canonical\.resources\s*\}\s*\)/s
  );
});

test('the embedded Python render bridge has no alternate production caller', () => {
  assert.deepEqual(
    occurrenceOwners(/globals\.get\(\s*['"]run_canonical_request_wrapper['"]\s*\)/g),
    new Map([['workers/diagram-generation-worker.js', 1]])
  );
  assert.deepEqual(
    occurrenceOwners(/\brender_embedded_canonical_web_request\b/g),
    new Map([['app/python-helpers.js', 2]])
  );
  assert.deepEqual(
    occurrenceOwners(/\bdef\s+run_canonical_request_wrapper\s*\(/g),
    new Map([['app/python-helpers.js', 1]])
  );
});
