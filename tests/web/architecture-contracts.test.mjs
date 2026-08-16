import assert from 'node:assert/strict';
import { spawnSync } from 'node:child_process';
import {
  existsSync,
  mkdirSync,
  mkdtempSync,
  readFileSync,
  readdirSync,
  rmSync,
  writeFileSync
} from 'node:fs';
import { tmpdir } from 'node:os';
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
    'app/legend/entry-actions.js',
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

test('privileged capabilities keep their current import boundaries', () => {
  const expectedImporters = new Map([
    ['services/session-request.js', [
      'app/run-analysis.js',
      'services/config.js',
      'services/gallery-session-migration.js'
    ]],
    ['app/python-helpers.js', ['workers/diagram-generation-worker.js']],
    ['services/diagram-worker-protocol.js', [
      'services/diagram-generation.js',
      'workers/diagram-generation-worker.js'
    ]],
    ['services/diagram-resource-staging.js', ['services/diagram-generation.js']],
    ['services/resource-payload-owner.js', [
      'services/config.js',
      'services/diagram-resource-staging.js',
      'services/session-request.js'
    ]],
    ['services/svg-serialization.js', [
      'app/app-setup.js',
      'app/feature-editor/label-actions.js',
      'app/feature-editor/svg-actions.js',
      'app/feature-search/preview-svg.js',
      'app/feature-selection.js',
      'app/legend-layout/canvas-actions.js',
      'app/legend-layout/diagram-drag.js',
      'app/legend-layout/reposition-actions.js',
      'app/legend/drag-actions.js',
      'app/legend/entry-actions.js',
      'app/legend/sort-actions.js',
      'app/legend/stroke-actions.js',
      'app/results.js',
      'app/svg-styles.js',
      'services/config.js',
      'services/export.js',
      'services/history-snapshot.js',
      'services/standalone-interactivity.js',
      'services/svg-result-ingestion.js'
    ]],
    ['services/history.js', ['app/app-setup.js']],
    ['services/history-snapshot.js', ['app/app-setup.js']],
    ['services/config.js', ['app/app-setup.js']],
    ['services/gallery-session-migration.js', ['services/config.js']],
    ['services/session-authority.js', [
      'services/config.js',
      'services/session-resources.js'
    ]],
    ['services/session-file.js', ['services/config.js']],
    ['app/feature-editor.js', ['app/app-setup.js']],
    ['app/legend.js', ['app/app-setup.js']],
    ['app/preview-runtime.js', ['app/app-setup.js']],
    ['app/right-drawer.js', ['app/app-setup.js', 'services/config.js']]
  ]);

  expectedImporters.forEach((expected, target) => {
    assert.deepEqual(importersOf(target), expected, target);
  });

  assert.deepEqual(
    [...occurrenceOwners(
      /\b(?:results|state\.results)\.value(?:\[[^\]]+\])?\s*=(?!=)/g
    ).keys()].sort(),
    [
      'app/feature-editor/label-actions.js',
      'app/feature-editor/svg-actions.js',
      'app/legend-layout/canvas-actions.js',
      'app/legend-layout/diagram-drag.js',
      'app/legend-layout/reposition-actions.js',
      'app/legend/drag-actions.js',
      'app/legend/entry-actions.js',
      'app/legend/sort-actions.js',
      'app/legend/stroke-actions.js',
      'app/preview-runtime.js',
      'app/results.js',
      'app/run-analysis.js',
      'app/svg-styles.js',
      'services/config.js'
    ]
  );
  assert.deepEqual(
    [...occurrenceOwners(
      /\b(?:createHistoryManager|beginCheckpoint|commitCheckpoint|buildArtifactCheckpoint)\b/g
    ).keys()].sort(),
    ['app/app-setup.js', 'services/history-snapshot.js', 'services/history.js']
  );
});

const CHANGE_BUDGET_CHECKER = join(REPOSITORY_ROOT, 'tools/check-web-change-budget.mjs');
const BUDGET_FIXTURE = Object.freeze({
  '.github/workflows/test.yml': 'name: baseline\n',
  'package.json': '{"private":true}\n',
  'tests/web/architecture-contracts.test.mjs': '// baseline contract\n',
  'tools/check-web-change-budget.mjs': '// baseline checker\n',
  'gbdraw/web/js/app/editor.js': 'export const editExistingOwner = () => 1;\n',
  'gbdraw/web/js/services/diagram-generation.js': (
    'export const runDiagramGeneration = () => 1;\n'
  )
});

const writeFixtureFile = (root, path, content) => {
  const target = join(root, path);
  mkdirSync(dirname(target), { recursive: true });
  writeFileSync(target, content, 'utf8');
};

const runChangeBudgetCase = (mutate, extraEnvironment = {}) => {
  const root = mkdtempSync(join(tmpdir(), 'gbdraw-web-budget-'));
  try {
    Object.entries(BUDGET_FIXTURE).forEach(([path, content]) => {
      writeFixtureFile(root, path, content);
    });
    const git = (args) => spawnSync('git', args, {
      cwd: root,
      encoding: 'utf8',
      stdio: ['ignore', 'pipe', 'pipe']
    });
    assert.equal(git(['init', '--quiet']).status, 0);
    assert.equal(git(['config', 'user.email', 'web-budget@example.invalid']).status, 0);
    assert.equal(git(['config', 'user.name', 'Web Budget Test']).status, 0);
    assert.equal(git(['add', '.']).status, 0);
    assert.equal(git(['commit', '--quiet', '-m', 'baseline']).status, 0);
    mutate((path, content) => writeFixtureFile(root, path, content));
    const result = spawnSync(process.execPath, [CHANGE_BUDGET_CHECKER], {
      cwd: root,
      encoding: 'utf8',
      env: { ...process.env, ...extraEnvironment },
      stdio: ['ignore', 'pipe', 'pipe']
    });
    return {
      status: result.status,
      output: `${result.stdout || ''}${result.stderr || ''}`
    };
  } finally {
    rmSync(root, { recursive: true, force: true });
  }
};

test('the Web change-budget checker rejects ordinary surface increases', () => {
  const cases = [
    {
      name: 'dummy production module',
      mutate: (write) => write('gbdraw/web/js/app/dummy.js', 'export const dummy = true;\n'),
      expected: /new production JavaScript modules are not allowed/
    },
    {
      name: 'diagram-generation importer from an editor module',
      mutate: (write) => write(
        'gbdraw/web/js/app/editor.js',
        "import { runDiagramGeneration } from '../services/diagram-generation.js';\n"
          + 'export const editExistingOwner = () => runDiagramGeneration();\n'
      ),
      expected: /Diagram Worker: importer .*app\/editor\.js/
    },
    {
      name: 'production dependency',
      mutate: (write) => write(
        'package.json',
        '{"private":true,"dependencies":{"left-pad":"1.3.0"}}\n'
      ),
      expected: /new production dependencies are not allowed/
    },
    {
      name: 'production and guard co-change',
      mutate: (write) => {
        write('gbdraw/web/js/app/editor.js', 'export const editExistingOwner = () => 2;\n');
        write('.github/workflows/test.yml', 'name: changed guard\n');
      },
      expected: /production runtime files and Web guard\/CI files changed together/
    }
  ];

  cases.forEach(({ name, mutate, expected }) => {
    const result = runChangeBudgetCase(mutate);
    assert.equal(result.status, 1, `${name}\n${result.output}`);
    assert.match(result.output, expected, name);
  });
});

test('the Web change-budget checker allows an owner-internal edit', () => {
  const result = runChangeBudgetCase((write) => {
    write(
      'gbdraw/web/js/services/diagram-generation.js',
      'export const runDiagramGeneration = () => 2;\n'
    );
  });
  assert.equal(result.status, 0, result.output);
  assert.match(result.output, /Result: \*\*PASS\*\*/);
});

test('architecture-change metadata waives budgets but not guard integrity', () => {
  const exception = { WEB_ARCHITECTURE_CHANGE: 'true' };
  const waived = runChangeBudgetCase((write) => {
    write('gbdraw/web/js/app/dummy.js', 'export const dummy = true;\n');
  }, exception);
  assert.equal(waived.status, 0, waived.output);
  assert.match(waived.output, /Rules waived by architecture-change/);

  const mixed = runChangeBudgetCase((write) => {
    write('gbdraw/web/js/app/editor.js', 'export const editExistingOwner = () => 2;\n');
    write('tools/check-web-change-budget.mjs', '// changed checker\n');
  }, exception);
  assert.equal(mixed.status, 1, mixed.output);
  assert.match(mixed.output, /production runtime files and Web guard\/CI files changed together/);
});
