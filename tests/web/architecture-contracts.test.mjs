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
    occurrenceOwners(/\brender_embedded_canonical_web_request\b/g),
    new Map([['app/python-helpers.js', 2]])
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
    'app/feature-editor/bulk-style-actions.js',
    'app/feature-editor/fill-actions.js',
    'app/feature-editor/fill-result-projection.js',
    'app/feature-editor/label-result-projection.js',
    'app/feature-editor/label-style-actions.js',
    'app/feature-editor/stroke-actions.js',
    'app/feature-editor/stroke-result-projection.js',
    'app/legend/manual-intent-command.js',
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

test('Feature fill UI exposes only planned request operations', () => {
  const editorSource = productionSources.get('app/feature-editor.js');
  const setupSource = productionSources.get('app/app-setup.js');
  const fillSource = productionSources.get('app/feature-editor/fill-actions.js');

  assert.doesNotMatch(
    INDEX_HTML,
    /\b(?:setFeatureColor|setFeatureColorValue|handleColorScopeChoice|resetClickedFeatureColor)\b/
  );
  assert.doesNotMatch(
    `${editorSource}\n${setupSource}`,
    /\b(?:setFeatureColor|setFeatureColorValue|handleColorScopeChoice|resetClickedFeatureColor)\s*[:,]/
  );
  assert.match(INDEX_HTML, /<feature-fill-color-control[\s\S]+@request="requestFeatureFillChange/);
  assert.match(INDEX_HTML, /v-for="candidate in pendingFeatureFillPlan\.candidates"/);
  assert.match(fillSource, /state\.manualLegendEntries\?\.value/);
  assert.doesNotMatch(fillSource, /featureIds\.length\s*===\s*0/);
});

test('atomic Feature fill preparation is browser-owned and has no Pyodide path', () => {
  for (const owner of [
    'app/feature-editor/fill-actions.js',
    'app/feature-editor/fill-command.js',
    'app/feature-editor/fill-result-projection.js',
    'app/legend/feature-fill-projection.js'
  ]) {
    assert.doesNotMatch(
      productionSources.get(owner),
      /\b(?:getPyodide|ensurePyodide|loadPyodide|runPython)\b/,
      owner
    );
  }
});

test('Feature style UI exposes planned facades and no retired raw setters', () => {
  const editorSource = productionSources.get('app/feature-editor.js');
  const setupSource = productionSources.get('app/app-setup.js');
  const publicSources = `${INDEX_HTML}\n${editorSource}\n${setupSource}`;
  assert.doesNotMatch(
    publicSources,
    /\b(?:setFeatureColor|setFeatureColorValue|handleColorScopeChoice|resetClickedFeatureColor|setClickedFeatureStrokeColorValue|updateClickedFeatureStroke|resetClickedFeatureStroke|applyStrokeToAllSiblings|applyStrokeToSelectedFeatures|handleLegendRenameChoice)\b/
  );
  assert.match(INDEX_HTML, /<feature-stroke-color-control[\s\S]+@request="requestFeatureStrokeChange/);
  assert.match(INDEX_HTML, /v-for="candidate in pendingFeatureStrokePlan\.candidates"/);
  assert.match(setupSource, /setFeatureLegendIntentHandler\(\(intent\) =>/);
  assert.match(setupSource, /requestManualLegendIntent\(/);
  assert.equal(productionSources.has('app/feature-editor/color-actions.js'), false);
});

test('Feature style command owners do not initialize Pyodide', () => {
  for (const owner of [
    'app/feature-editor/bulk-style-actions.js',
    'app/feature-editor/bulk-style-command.js',
    'app/feature-editor/stroke-actions.js',
    'app/feature-editor/stroke-command.js',
    'app/feature-editor/stroke-result-projection.js',
    'app/legend/manual-intent-command.js'
  ]) {
    assert.doesNotMatch(
      productionSources.get(owner),
      /\b(?:getPyodide|ensurePyodide|loadPyodide|runPython)\b/,
      owner
    );
  }
});

test('bulk rule and palette writers use the atomic snapshot facade', () => {
  const rules = productionSources.get('app/feature-editor/rule-actions.js');
  const watchers = productionSources.get('app/watchers.js');
  const results = productionSources.get('app/results.js');
  assert.match(rules, /requestFeatureBulkStyleChange/);
  assert.match(watchers, /writerKind: 'specific-color-tsv-import'/);
  assert.match(watchers, /writerKind: 'default-color-tsv-import'/);
  assert.match(results, /writerKind: 'palette-default-acceptance'/);
  assert.doesNotMatch(
    watchers,
    /watch\(\(\) => \[\.\.\.manualSpecificRules\][\s\S]{0,500}applySpecificRulesToSvg/
  );
});

test('managed exact and semantic Feature rules are read-only in the generic editor', () => {
  const rules = productionSources.get('app/feature-editor/rule-actions.js');
  const editor = productionSources.get('app/feature-editor.js');
  const setup = productionSources.get('app/app-setup.js');

  assert.match(rules, /isExactFeatureRule\(rule\) \|\| isFeatureSemanticScopeRule\(rule\)/);
  assert.match(rules, /if \(isOpaqueSpecificRule\(current\)\) return false/);
  assert.match(rules, /manualSpecificRules\.filter\(isOpaqueSpecificRule\)/);
  assert.match(`${editor}\n${setup}`, /canMoveSpecificRule/);
  assert.match(`${editor}\n${setup}`, /isOpaqueSpecificRule/);
  assert.match(INDEX_HTML, /data-managed-specific-rule/);
  assert.match(INDEX_HTML, /:disabled="isOpaqueSpecificRule\(rule\)"/);
  assert.match(INDEX_HTML, /:disabled="!canMoveSpecificRule\(i, -1\)"/);
  assert.match(INDEX_HTML, /:disabled="!canMoveSpecificRule\(i, 1\)"/);
});

test('legend order is document-global intent projected only by atomic owners', () => {
  const sortActions = productionSources.get('app/legend/sort-actions.js');
  const command = productionSources.get('app/legend/manual-intent-command.js');
  const candidate = productionSources.get('app/candidate-render.js');
  const config = productionSources.get('services/config.js');

  assert.match(sortActions, /requestLegendOrderIntent\(\{ kind: 'order', order: normalized \}/);
  assert.doesNotMatch(sortActions, /setAttribute|serializeCleanSvg|results\.value/);
  assert.match(command, /prepareLegendOrderResultProjection/);
  assert.match(command, /writeRef\(state\.results, prepared\.projection\.nextResults/);
  assert.match(command, /owner: isOrderIntent \? ORDER_OWNER : MANUAL_OWNER/);
  assert.doesNotMatch(command, /\b(?:getPyodide|ensurePyodide|loadPyodide|runPython)\b/);
  assert.match(candidate, /applyLegendOrderToSvg\(svg, legendOrderIntent\)/);
  assert.match(config, /orderIntent: cloneJsonArray\(state\.legendOrderIntent\?\.value\)/);
  assert.match(config, /orderIntent: normalizeStringArray\(legend\.orderIntent\)/);
});
