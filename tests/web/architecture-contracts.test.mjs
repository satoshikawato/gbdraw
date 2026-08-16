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
import {
  literalImportSpecifiers,
  maskJavaScript
} from '../../tools/web-change-source.mjs';

const REPOSITORY_ROOT = resolve(dirname(fileURLToPath(import.meta.url)), '../..');
const JAVASCRIPT_ROOT = join(REPOSITORY_ROOT, 'gbdraw/web/js');
const APP_ENTRY = join(JAVASCRIPT_ROOT, 'app.js');
const INDEX_HTML = readFileSync(
  join(REPOSITORY_ROOT, 'gbdraw/web/index.html'),
  'utf8'
);
const WEB_CHANGE_POLICY = JSON.parse(readFileSync(
  join(REPOSITORY_ROOT, 'tools/web-change-policy.json'),
  'utf8'
));
const TEST_WORKFLOW = readFileSync(
  join(REPOSITORY_ROOT, '.github/workflows/test.yml'),
  'utf8'
);
const BASE_POLICY_WORKFLOW = readFileSync(
  join(REPOSITORY_ROOT, '.github/workflows/web-base-policy.yml'),
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

const staticImportSpecifiers = (source) => literalImportSpecifiers(source, { dynamic: false });

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

const executableOccurrenceOwners = (pattern) => new Map(
  [...productionSources]
    .map(([owner, source]) => [owner, [...maskJavaScript(source).matchAll(pattern)].length])
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

test('privileged capability owners and importers stay within their allowlists', () => {
  const ownerPatterns = new Map([
    ['Render request', /\bbuildCanonicalRenderRequest\s*\(/g],
    ['Diagram Worker', /\b(?:new\s+Worker|runDiagramGeneration|loadPyodide)\s*\(/g],
    ['Python helper', /\b(?:runPythonAsync|runPython|PYTHON_HELPERS|globals\.get)\b/g],
    [
      'Resource staging',
      /\b(?:createDiagramResourceTransport|stageRenderResources|setResourcePayloadOwner)\b/g
    ],
    [
      'SVG/Result admission',
      /\b(?:sanitizeSvgContent|ingestSvgResult|markCommitted)\b/g
    ],
    [
      'Mounted SVG/Result replacement',
      /\b(?:serializeCleanSvg|flushActiveResult)\s*\(|\b(?:results|state\.results)\.value(?:\[[^\]]+\])?\s*=(?!=)/g
    ],
    [
      'History',
      /\b(?:createHistoryManager|beginCheckpoint|commitCheckpoint|buildArtifactCheckpoint)\b/g
    ],
    [
      'Session',
      /\b(?:exportSession|compressSessionData|migrateSessionDataToCurrent|migratePersistedGalleryConfig)\b/g
    ],
    [
      'Canonical editor state',
      /\b(?:createFeatureEditor|createLegendManager|createRightDrawerController|createPreviewRuntime)\b/g
    ]
  ]);
  const assertSubset = (actual, allowed, label) => {
    const allowedSet = new Set(allowed);
    assert.deepEqual(actual.filter((path) => !allowedSet.has(path)), [], label);
  };

  Object.entries(WEB_CHANGE_POLICY.allowedPrivilegedImporters).forEach(([target, allowed]) => {
    assertSubset(importersOf(target), allowed, `${target} importers`);
  });
  assert.deepEqual(
    Object.keys(WEB_CHANGE_POLICY.allowedPrivilegedOwners).sort(),
    [...ownerPatterns.keys()].sort()
  );
  ownerPatterns.forEach((pattern, capability) => {
    const actual = [...executableOccurrenceOwners(pattern).keys()].sort();
    assertSubset(actual, WEB_CHANGE_POLICY.allowedPrivilegedOwners[capability], capability);
  });
});

test('PR workflows separate normal tests from trusted base-policy execution', () => {
  const activityTypes = /types: \[opened, synchronize, reopened, labeled, unlabeled\]/;
  assert.match(TEST_WORKFLOW, /\n  pull_request:\n/);
  assert.match(TEST_WORKFLOW, activityTypes);
  assert.match(TEST_WORKFLOW, /name: Run core tests/);

  assert.match(BASE_POLICY_WORKFLOW, /\n  pull_request_target:\n/);
  assert.match(BASE_POLICY_WORKFLOW, activityTypes);
  assert.match(BASE_POLICY_WORKFLOW, /permissions:\n  contents: read/);
  assert.match(BASE_POLICY_WORKFLOW, /name: Web base policy \(trusted base\)/);
  assert.match(BASE_POLICY_WORKFLOW, /ref: \$\{\{ github\.event\.pull_request\.base\.sha \}\}/);
  assert.match(BASE_POLICY_WORKFLOW, /persist-credentials: false/);
  assert.equal([...BASE_POLICY_WORKFLOW.matchAll(/uses: actions\/checkout@/g)].length, 1);
  assert.match(BASE_POLICY_WORKFLOW, /git fetch --no-tags --depth=1 origin "\$HEAD_SHA"/);
  assert.match(
    BASE_POLICY_WORKFLOW,
    /node tools\/check-web-change-budget\.mjs --base "\$BASE_SHA" --head "\$HEAD_SHA"/
  );
  assert.doesNotMatch(
    BASE_POLICY_WORKFLOW,
    /(?:checkout|run):[^\n]*pull_request\.head|ref:.*pull_request\.head|git (?:checkout|switch) /
  );
});

const CHANGE_BUDGET_CHECKER = join(REPOSITORY_ROOT, 'tools/check-web-change-budget.mjs');
const BUDGET_POLICY = Object.freeze({
  allowedPrivilegedImporters: {
    'services/diagram-generation.js': ['app/editor.js'],
    'workers/diagram-generation-worker.js': []
  },
  allowedPrivilegedOwners: {
    'Diagram Worker': [
      'app/editor.js',
      'app/future-owner.js',
      'services/diagram-generation.js'
    ]
  }
});
const BUDGET_FIXTURE = Object.freeze({
  '.github/workflows/test.yml': 'name: baseline\n',
  '.github/workflows/web-base-policy.yml': 'name: baseline policy\n',
  'package.json': '{"private":true}\n',
  'tests/web/architecture-contracts.test.mjs': '// baseline contract\n',
  'tools/check-web-change-budget.mjs': readFileSync(CHANGE_BUDGET_CHECKER, 'utf8'),
  'tools/web-change-source.mjs': readFileSync(
    join(REPOSITORY_ROOT, 'tools/web-change-source.mjs'),
    'utf8'
  ),
  'tools/web-change-policy.json': `${JSON.stringify(BUDGET_POLICY, null, 2)}\n`,
  'gbdraw/web/index.html': '<main>baseline</main>\n',
  'gbdraw/web/js/app/editor.js': (
    "import { runDiagramGeneration } from '../services/diagram-generation.js';\n"
    + 'export const editExistingOwner = () => runDiagramGeneration();\n'
  ),
  'gbdraw/web/js/app/future-owner.js': 'export const futureOwnerPlaceholder = () => 1;\n',
  'gbdraw/web/js/app/secondary.js': 'export const editSecondaryOwner = () => 1;\n',
  'gbdraw/web/js/services/diagram-generation.js': (
    'export const runDiagramGeneration = () => 1;\n'
  ),
  'gbdraw/web/js/services/session-file.js': (
    'export const readSession = () => ({ version: 1 });\n'
  ),
  'gbdraw/web/vendor/library.js': 'globalThis.vendorLibrary = true;\n'
});

const writeFixtureFile = (root, path, content) => {
  const target = join(root, path);
  mkdirSync(dirname(target), { recursive: true });
  if (Buffer.isBuffer(content)) writeFileSync(target, content);
  else writeFileSync(target, content, 'utf8');
};

const withChangeBudgetRepository = (runCase) => {
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
    const commit = (message) => {
      assert.equal(git(['add', '.']).status, 0);
      assert.equal(git(['commit', '--quiet', '-m', message]).status, 0);
      const revision = git(['rev-parse', 'HEAD']);
      assert.equal(revision.status, 0);
      return revision.stdout.trim();
    };
    const execute = ({ base = '', head = '', environment = {} } = {}) => {
      const arguments_ = [CHANGE_BUDGET_CHECKER];
      if (base) arguments_.push('--base', base);
      if (head) arguments_.push('--head', head);
      const result = spawnSync(process.execPath, arguments_, {
        cwd: root,
        encoding: 'utf8',
        env: { ...process.env, ...environment },
        stdio: ['ignore', 'pipe', 'pipe']
      });
      return {
        status: result.status,
        output: `${result.stdout || ''}${result.stderr || ''}`
      };
    };
    return runCase({
      commit,
      execute,
      git,
      root,
      write: (path, content) => writeFixtureFile(root, path, content)
    });
  } finally {
    rmSync(root, { recursive: true, force: true });
  }
};

const runChangeBudgetCase = (mutate, environment = {}) => withChangeBudgetRepository(
  ({ execute, write }) => {
    mutate(write);
    return execute({ environment });
  }
);

const assertNonWaivableRevisionFailure = (execute, base, head, expectedPatterns) => {
  for (const environment of [{}, { WEB_ARCHITECTURE_CHANGE: 'true' }]) {
    const result = execute({ base, head, environment });
    assert.equal(result.status, 1, result.output);
    expectedPatterns.forEach((pattern) => assert.match(result.output, pattern));
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
      name: 'production dependency',
      mutate: (write) => write(
        'package.json',
        '{"private":true,"dependencies":{"left-pad":"1.3.0"}}\n'
      ),
      expected: /new production dependencies are not allowed/
    },
    {
      name: 'binary runtime file',
      mutate: (write) => write('gbdraw/web/js/runtime.bin', Buffer.from([0, 1, 2, 3])),
      expected: /added binary runtime files are not allowed/
    },
    {
      name: 'vendored runtime change',
      mutate: (write) => write(
        'gbdraw/web/vendor/library.js',
        'globalThis.vendorLibrary = false;\n'
      ),
      expected: /changes under gbdraw\/web\/vendor\/ are not allowed/
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

test('privileged ownership and import contractions do not require a policy edit', () => {
  const result = runChangeBudgetCase((write) => {
    write(
      'gbdraw/web/js/app/editor.js',
      'export const editExistingOwner = () => 1;\n'
    );
  });
  assert.equal(result.status, 0, result.output);
  assert.match(result.output, /Result: \*\*PASS\*\*/);
});

test('unapproved privileged expansion fails against the base allowlist', () => {
  const result = runChangeBudgetCase((write) => {
    write(
      'gbdraw/web/js/app/secondary.js',
      "import { runDiagramGeneration } from '../services/diagram-generation.js';\n"
        + 'export const editSecondaryOwner = () => runDiagramGeneration();\n'
    );
  });
  assert.equal(result.status, 1, result.output);
  assert.match(
    result.output,
    /services\/diagram-generation\.js: importer app\/secondary\.js/
  );
});

test('an unused privileged owner allowlist entry can contract', () => {
  withChangeBudgetRepository(({ commit, execute, git, write }) => {
    const base = git(['rev-parse', 'HEAD']).stdout.trim();
    const policy = JSON.parse(JSON.stringify(BUDGET_POLICY));
    policy.allowedPrivilegedOwners['Diagram Worker'] = policy.allowedPrivilegedOwners[
      'Diagram Worker'
    ].filter((path) => path !== 'app/future-owner.js');
    write('tools/web-change-policy.json', `${JSON.stringify(policy, null, 2)}\n`);

    const head = commit('remove inactive owner allowlist entry');
    const result = execute({ base, head });
    assert.equal(result.status, 0, result.output);
    assert.match(result.output, /Result: \*\*PASS\*\*/);
    assert.match(
      result.output,
      /allowedPrivilegedOwners\.Diagram Worker: app\/future-owner\.js/
    );
    assert.doesNotMatch(
      result.output,
      /proposed privileged capability allowlist excludes active owners or importers/
    );
  });
});

test('an active privileged owner allowlist entry cannot contract', () => {
  withChangeBudgetRepository(({ commit, execute, git, write }) => {
    const base = git(['rev-parse', 'HEAD']).stdout.trim();
    const policy = JSON.parse(JSON.stringify(BUDGET_POLICY));
    policy.allowedPrivilegedOwners['Diagram Worker'] = policy.allowedPrivilegedOwners[
      'Diagram Worker'
    ].filter((path) => path !== 'app/editor.js');
    write('tools/web-change-policy.json', `${JSON.stringify(policy, null, 2)}\n`);

    const head = commit('remove active owner allowlist entry');
    assertNonWaivableRevisionFailure(execute, base, head, [
      /Diagram Worker: owner app\/editor\.js/,
      /proposed privileged capability allowlist excludes active owners or importers/
    ]);
  });
});

test('an active privileged importer allowlist entry cannot contract', () => {
  withChangeBudgetRepository(({ commit, execute, git, write }) => {
    const base = git(['rev-parse', 'HEAD']).stdout.trim();
    const policy = JSON.parse(JSON.stringify(BUDGET_POLICY));
    policy.allowedPrivilegedImporters['services/diagram-generation.js'] = [];
    write('tools/web-change-policy.json', `${JSON.stringify(policy, null, 2)}\n`);

    const head = commit('remove active importer allowlist entry');
    assertNonWaivableRevisionFailure(execute, base, head, [
      /services\/diagram-generation\.js: importer app\/editor\.js/,
      /proposed privileged capability allowlist excludes active owners or importers/
    ]);
  });
});

test('a policy change cannot self-authorize privileged runtime expansion', () => {
  withChangeBudgetRepository(({ commit, execute, git, write }) => {
    const base = git(['rev-parse', 'HEAD']).stdout.trim();
    const policy = JSON.parse(JSON.stringify(BUDGET_POLICY));
    policy.allowedPrivilegedImporters['services/diagram-generation.js'].push('app/secondary.js');
    policy.allowedPrivilegedOwners['Diagram Worker'].push('app/secondary.js');
    write('tools/web-change-policy.json', `${JSON.stringify(policy, null, 2)}\n`);
    write(
      'gbdraw/web/js/app/secondary.js',
      "import { runDiagramGeneration } from '../services/diagram-generation.js';\n"
        + 'export const editSecondaryOwner = () => runDiagramGeneration();\n'
    );

    const head = commit('expand runtime and proposed allowlist together');
    assertNonWaivableRevisionFailure(execute, base, head, [
      /production runtime files and Web guard\/CI files changed together/,
      /privileged capability owners or importers exceed the base allowlist/,
      /services\/diagram-generation\.js: importer app\/secondary\.js/,
      /Diagram Worker: owner app\/secondary\.js/
    ]);
  });
});

test('base privileged allowlist keys cannot disappear', () => {
  withChangeBudgetRepository(({ commit, execute, git, write }) => {
    const base = git(['rev-parse', 'HEAD']).stdout.trim();
    const policy = JSON.parse(JSON.stringify(BUDGET_POLICY));
    delete policy.allowedPrivilegedImporters['workers/diagram-generation-worker.js'];
    write('tools/web-change-policy.json', `${JSON.stringify(policy, null, 2)}\n`);

    const head = commit('delete empty importer allowlist key');
    assertNonWaivableRevisionFailure(execute, base, head, [
      /allowedPrivilegedImporters\.workers\/diagram-generation-worker\.js/,
      /proposed privileged capability policy is missing base allowlist keys/
    ]);
  });
});

test('privileged expansion requires policy preauthorization on the base revision', () => {
  withChangeBudgetRepository(({ commit, execute, git, write }) => {
    const policy = JSON.parse(JSON.stringify(BUDGET_POLICY));
    policy.allowedPrivilegedImporters['services/diagram-generation.js'].push('app/secondary.js');
    policy.allowedPrivilegedOwners['Diagram Worker'].push('app/secondary.js');
    write('tools/web-change-policy.json', `${JSON.stringify(policy, null, 2)}\n`);

    const preauthorization = execute();
    assert.equal(preauthorization.status, 0, preauthorization.output);
    const preauthorizedBase = commit('preauthorize secondary owner');

    write(
      'gbdraw/web/js/app/secondary.js',
      "import { runDiagramGeneration } from '../services/diagram-generation.js';\n"
        + 'export const editSecondaryOwner = () => runDiagramGeneration();\n'
    );
    const implementationHead = commit('use preauthorized secondary owner');
    const implementation = execute({ base: preauthorizedBase, head: implementationHead });
    assert.equal(implementation.status, 0, implementation.output);
    assert.match(implementation.output, /Result: \*\*PASS\*\*/);

    assert.equal(git(['diff', '--quiet', `${preauthorizedBase}^`, preauthorizedBase, '--', 'gbdraw/web']).status, 0);
  });
});

test('comments, strings, and local session object keys are not hard failures', () => {
  const result = runChangeBudgetCase((write) => {
    write(
      'gbdraw/web/js/services/session-file.js',
      'export const readSession = () => {\n'
        + '  // export const createCacheManager = () => watch(runDiagramGeneration);\n'
        + '  const note = "import runPython from \'./app/python-helpers.js\'";\n'
        + '  return { version: 1, journalToken: note };\n'
        + '};\n'
    );
  });
  assert.equal(result.status, 0, result.output);
  assert.match(result.output, /Report-only session object keys and compatibility names/);
});

test('index.html growth counts against the runtime line budget', () => {
  const result = runChangeBudgetCase((write) => {
    const additions = Array.from({ length: 110 }, (_, index) => `<div>${index}</div>`);
    write('gbdraw/web/index.html', `<main>baseline</main>\n${additions.join('\n')}\n`);
  });
  assert.equal(result.status, 1, result.output);
  assert.match(result.output, /M gbdraw\/web\/index\.html/);
  assert.match(result.output, /production additions exceed deletions by more than 100 lines/);
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
    write('tools/web-change-policy.json', `${JSON.stringify(BUDGET_POLICY)}\n`);
  }, exception);
  assert.equal(mixed.status, 1, mixed.output);
  assert.match(mixed.output, /production runtime files and Web guard\/CI files changed together/);
});

test('base-policy execution ignores a PR-modified checker and detects its runtime co-change', () => {
  withChangeBudgetRepository(({ commit, execute, git, write }) => {
    const base = git(['rev-parse', 'HEAD']).stdout.trim();
    write('gbdraw/web/js/app/editor.js', 'export const editExistingOwner = () => 2;\n');
    write('tools/check-web-change-budget.mjs', 'process.exit(0);\n');
    const head = commit('modify runtime and checker');

    const result = execute({
      base,
      head,
      environment: { WEB_ARCHITECTURE_CHANGE: 'true' }
    });
    assert.equal(result.status, 1, result.output);
    assert.match(result.output, /production runtime files and Web guard\/CI files changed together/);
  });
});
