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
const INDEX_HTML = readFileSync(
  join(REPOSITORY_ROOT, 'gbdraw/web/index.html'),
  'utf8'
);

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

const staticImportSpecifiers = (source) => {
  const specifiers = [];
  const patterns = [
    /(?:^|\n)\s*(?:import|export)\s+[^;]*?\s+from\s+(['"])([^'"]+)\1\s*;?/g,
    /(?:^|\n)\s*import\s*(['"])([^'"]+)\1\s*;?/g
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

const staticDirectImports = new Map(
  [...productionSources].map(([owner, source]) => {
    const ownerPath = join(JAVASCRIPT_ROOT, owner);
    const imports = staticImportSpecifiers(source)
      .map((specifier) => resolveRelativeImport(ownerPath, specifier))
      .filter(Boolean)
      .map(relativeModulePath);
    return [owner, new Set(imports)];
  })
);

const reachableModules = (entryPath, importGraph = directImports) => {
  const visited = new Set();
  const pending = [relativeModulePath(entryPath)];
  while (pending.length) {
    const current = pending.pop();
    if (visited.has(current)) continue;
    visited.add(current);
    for (const dependency of importGraph.get(current) || []) pending.push(dependency);
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
    'app/record-discovery.js',
    'app/results.js',
    'app/run-analysis.js'
  ]);
});

test('Pyodide initialization and helper execution are Worker-only', () => {
  assert.deepEqual(
    occurrenceOwners(/\bloadPyodide\b/g),
    new Map([['workers/diagram-generation-worker.js', 2]])
  );
  assert.deepEqual(
    occurrenceOwners(/\bloadPyodide\s*\(/g),
    new Map([['workers/diagram-generation-worker.js', 1]])
  );
  const reachable = reachableModules(APP_ENTRY);
  assert.deepEqual(
    [...occurrenceOwners(/\bloadPyodide\s*\(/g).keys()].filter((path) => reachable.has(path)),
    []
  );
  assert.equal(productionSources.has('app/pyodide.js'), false);
  assert.doesNotMatch(INDEX_HTML, /<script\b[^>]*\bsrc\s*=\s*['"][^'"]*pyodide\.js/i);
  const mainThreadSource = [...productionSources]
    .filter(([path]) => reachable.has(path))
    .map(([, source]) => source)
    .join('\n');
  assert.doesNotMatch(mainThreadSource, /\bcreatePyodideManager\b|\bgetPyodide\b|\bensurePyodide\b/);
  assert.doesNotMatch(
    mainThreadSource,
    /\bloadPyodide\s*\(|\bPYTHON_HELPERS\b|\brunPython(?:Async)?\s*\(|\.globals\.get\s*\(/
  );
  assert.match(
    productionSources.get('services/diagram-worker-protocol.js'),
    /DIAGRAM_HELPER_OPERATIONS\s*=\s*Object\.freeze/
  );
  assert.match(
    productionSources.get('workers/diagram-generation-worker.js'),
    /HELPER_OPERATION_SPECS[\s\S]+assertAllowedPayloadKeys/
  );
});

test('export payloads stay outside the eager application graph', () => {
  const eagerModules = reachableModules(APP_ENTRY, staticDirectImports);
  const appSetupSource = productionSources.get('app/app-setup.js');
  const exportSource = productionSources.get('services/export.js');

  assert.ok(eagerModules.has('app/app-setup.js'));
  assert.equal(eagerModules.has('services/export.js'), false);
  assert.equal(eagerModules.has('services/standalone-interactivity.js'), false);
  assert.equal(eagerModules.has('services/standalone-interactivity-assets.js'), false);
  assert.match(
    appSetupSource,
    /exportServicePromise\s*\?\?=\s*import\(\s*['"]\.\.\/services\/export\.js['"]\s*\)/
  );
  assert.match(
    exportSource,
    /standaloneInteractivityPromise\s*\?\?=\s*import\(\s*['"]\.\/standalone-interactivity\.js['"]\s*\)/
  );
  assert.match(exportSource, /pdfLibrariesPromise\s*\?\?=/);
  assert.doesNotMatch(
    INDEX_HTML,
    /<script\b[^>]*\bsrc\s*=\s*['"][^'"]*(?:jspdf|svg2pdf)[^'"]*['"][^>]*>/i
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
    occurrenceOwners(/\brender_staged_canonical_web_request\b/g),
    new Map([['app/python-helpers.js', 2]])
  );
  assert.doesNotMatch(
    productionSources.get('workers/diagram-generation-worker.js'),
    /JSON\.stringify\(resources\)/
  );
  assert.deepEqual(
    occurrenceOwners(/\bdef\s+run_canonical_request_wrapper\s*\(/g),
    new Map([['app/python-helpers.js', 1]])
  );
});

test('History intent and SVG admission have one production ownership path', () => {
  assert.deepEqual(
    importersOf('services/svg-sanitization.js'),
    ['services/svg-result-ingestion.js']
  );
  assert.deepEqual(
    occurrenceOwners(/\bsanitizeSvgContent\s*\(/g),
    new Map([['services/svg-result-ingestion.js', 1]])
  );
  assert.deepEqual(importersOf('services/svg-result-ingestion.js'), [
    'app/candidate-render.js',
    'app/preview-runtime.js',
    'app/session-feature-metadata.js',
    'app/watchers.js',
    'services/config.js',
    'state.js'
  ]);
  assert.deepEqual(importersOf('services/session-feature-metadata.js'), [
    'app/session-feature-metadata.js',
    'services/svg-result-ingestion.js'
  ]);
  assert.deepEqual(importersOf('services/svg-result-normalization.js'), [
    'app/svg-styles.js',
    'services/config.js',
    'services/svg-result-ingestion.js'
  ]);
  assert.doesNotMatch(
    productionSources.get('app/session-feature-metadata.js'),
    /DOMParser|parseFromString|result\?\.content/
  );
  assert.match(
    productionSources.get('app/feature-search/preview-svg.js'),
    /getPreviewFeatureElementIndex\s*=\s*\(svg\)\s*=>\s*getFeatureElementIndex\(svg\)/
  );
  assert.doesNotMatch(
    productionSources.get('app/feature-search/preview-svg.js'),
    /rebuild\s*:\s*true|buildFeatureElementIndex/
  );
  assert.doesNotMatch(
    productionSources.get('app/watchers.js'),
    /ensureUniqueSkewClipPathIds|ensureUniquePairwiseGradientIds|reapplyStrokeOverrides/
  );
  assert.doesNotMatch(
    productionSources.get('app/legend/entry-actions.js'),
    /normalizeLegacyLegendEntryGroups/
  );
  assert.doesNotMatch(
    productionSources.get('app/app-setup.js'),
    /normalizeLegacySvg|normalizeLegacyComposition/
  );
  assert.match(
    productionSources.get('services/config.js'),
    /normalizeLegacyLegendEntryGroups\(svg\)[\s\S]+normalizeLegacyComposition\(svg,/
  );
  assert.deepEqual(occurrenceOwners(/\bbuildHistorySnapshot\b/g), new Map());
  assert.deepEqual(occurrenceOwners(/\bapplyHistorySnapshot\b/g), new Map());
});

test('live SVG edits and fresh-SVG replay have explicit single owners', () => {
  const directEditorModules = [
    'app/feature-editor/color-actions.js',
    'app/feature-editor/label-actions.js',
    'app/feature-editor/svg-actions.js',
    'app/feature-editor/visibility-actions.js',
    'app/legend/entry-actions.js',
    'app/legend/stroke-actions.js',
    'app/legend/drag-actions.js',
    'app/legend/sort-actions.js',
    'app/legend-layout/canvas-actions.js',
    'app/legend-layout/diagram-drag.js',
    'app/legend-layout/reposition-actions.js',
    'app/results.js',
    'app/svg-styles.js'
  ];
  const ordinaryDirectActionModules = directEditorModules.filter(
    (owner) => owner !== 'app/results.js'
  );
  directEditorModules.forEach((owner) => {
    const source = productionSources.get(owner) || '';
    assert.doesNotMatch(source, /\bserializeCleanSvg\b/, `${owner} must use the live edit owner`);
    assert.doesNotMatch(
      source,
      /(?:state\.)?results\.value\[[^\]]+\]\s*=/,
      `${owner} must not replace Result SVG content`
    );
  });
  ordinaryDirectActionModules.forEach((owner) => {
    assert.doesNotMatch(
      productionSources.get(owner) || '',
      /diagram-generation|run-analysis/,
      `${owner} must not own rendering`
    );
  });

  assert.deepEqual(
    occurrenceOwners(/state\.results\.value\[[^\]]+\]\s*=/g),
    new Map([['app/preview-runtime.js', 1]])
  );
  assert.deepEqual(
    occurrenceOwners(/\.flushActiveResult\s*\(/g),
    new Map([
      ['app/app-setup.js', 1],
      ['services/config.js', 1]
    ])
  );

  const projectionSource = productionSources.get('app/editor-svg-projection.js') || '';
  [
    'featureColorOverrides',
    'featureStrokeOverrides',
    'featureVisibilityOverrides',
    'labelTextFeatureOverrides',
    'labelVisibilityOverrides',
    'legendColorOverrides',
    'legendStrokeOverrides',
    'manualSpecificRules',
    'suppressPairwiseIdentityLegend'
  ].forEach((ownerToken) => assert.match(projectionSource, new RegExp(`\\b${ownerToken}\\b`)));
  assert.ok(
    directImports.get('app/candidate-render.js')?.has('app/editor-svg-projection.js'),
    'Generate and reflow admission must use the projection owner'
  );
  assert.ok(
    directImports.get('services/config.js')?.has('app/editor-svg-projection.js'),
    'session replay must use the projection owner when replay is required'
  );

  ordinaryDirectActionModules.forEach((owner) => {
    const reachable = reachableModules(join(JAVASCRIPT_ROOT, owner), staticDirectImports);
    assert.equal(reachable.has('services/diagram-generation.js'), false, owner);
    assert.equal(reachable.has('app/run-analysis.js'), false, owner);
    assert.deepEqual([...reachable].filter((path) => path.startsWith('workers/')), [], owner);
  });
});

test('right drawer availability and transitions have one production owner', () => {
  const transitionOwners = [
    ...occurrenceOwners(/\b(?:showRightDrawer|rightDrawerTab)\.value\s*=(?!=)/g).keys()
  ];

  assert.deepEqual(transitionOwners, ['app/right-drawer.js']);
  assert.deepEqual(
    occurrenceOwners(/\b(?:showFeaturePanel|showLegendPanel)\b/g),
    new Map()
  );
  assert.match(INDEX_HTML, /@click="toggleRightDrawer"/);
  assert.match(
    INDEX_HTML,
    /:disabled="!isRightDrawerTabAvailable\('orthogroups'\)"/
  );
  assert.doesNotMatch(INDEX_HTML, /rightDrawerTab\s*\|\|/);
});
