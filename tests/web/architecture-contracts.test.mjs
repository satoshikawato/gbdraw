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
  detectPrivilegedWebCapabilities,
  detectReportOnlySourceFacts,
  WEB_ARCHITECTURE_DETECTORS,
  WEB_PRIVILEGED_CAPABILITY_KEYS,
  WEB_PRIVILEGED_IMPORT_TARGETS
} from '../../tools/web-architecture-detectors.mjs';
import {
  classifyArchitectureRuleObservation,
  evaluateArchitectureRuleResult,
  validateArchitectureRuleRegistry
} from '../../tools/web-architecture-evaluation.mjs';
import {
  validateProductDecisionAuthority,
  validateProductImpactMap
} from '../../tools/web-product-impact-evaluation.mjs';
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
const GALLERY_PUBLICATION_WORKFLOW = readFileSync(
  join(REPOSITORY_ROOT, '.github/workflows/gallery-publication.yml'),
  'utf8'
);
const DEPLOY_WORKFLOW = readFileSync(
  join(REPOSITORY_ROOT, '.github/workflows/deploy_web.yml'),
  'utf8'
);
const PACKAGE_SCRIPTS = JSON.parse(readFileSync(
  join(REPOSITORY_ROOT, 'package.json'),
  'utf8'
)).scripts;
const BASE_PLAYWRIGHT_CONFIG = readFileSync(
  join(REPOSITORY_ROOT, 'playwright.config.js'),
  'utf8'
);
const FUNCTIONAL_PLAYWRIGHT_CONFIG = readFileSync(
  join(REPOSITORY_ROOT, 'playwright.functional.config.js'),
  'utf8'
);
const PR_SMOKE_PLAYWRIGHT_CONFIG = readFileSync(
  join(REPOSITORY_ROOT, 'playwright.pr-smoke.config.js'),
  'utf8'
);
const WORKFLOW_NAMES = readdirSync(
  join(REPOSITORY_ROOT, '.github/workflows'),
  { withFileTypes: true }
).filter((entry) => entry.isFile() && /\.ya?ml$/.test(entry.name))
  .map((entry) => entry.name)
  .sort();
const WORKFLOW_SOURCES = WORKFLOW_NAMES
  .map((name) => readFileSync(
    join(REPOSITORY_ROOT, '.github/workflows', name),
    'utf8'
  ));
const PULL_REQUEST_TEMPLATE = readFileSync(
  join(REPOSITORY_ROOT, '.github/pull_request_template.md'),
  'utf8'
);
const WEB_ARCHITECTURE_EVALUATION_SOURCE = readFileSync(
  join(REPOSITORY_ROOT, 'tools/web-architecture-evaluation.mjs'),
  'utf8'
);
const WEB_ARCHITECTURE_DETECTOR_SOURCE = readFileSync(
  join(REPOSITORY_ROOT, 'tools/web-architecture-detectors.mjs'),
  'utf8'
);
const WEB_ARCHITECTURE_RULES_SOURCE = readFileSync(
  join(REPOSITORY_ROOT, 'tools/web-architecture-rules.json'),
  'utf8'
);
const WEB_ARCHITECTURE_RULES = JSON.parse(WEB_ARCHITECTURE_RULES_SOURCE);
const WEB_CHANGE_BUDGET_SOURCE = readFileSync(
  join(REPOSITORY_ROOT, 'tools/check-web-change-budget.mjs'),
  'utf8'
);
const PRODUCT_IMPACT_EVALUATION_SOURCE = readFileSync(
  join(REPOSITORY_ROOT, 'tools/web-product-impact-evaluation.mjs'),
  'utf8'
);
const PRODUCT_IMPACT_DECISION_SOURCE = readFileSync(
  join(REPOSITORY_ROOT, 'tools/web-product-impact-decision-source.mjs'),
  'utf8'
);
const PRODUCT_IMPACT_MAP_PATH = join(
  REPOSITORY_ROOT,
  'tools/web-product-impact-map.json'
);
const PRODUCT_DECISIONS_PATH = join(
  REPOSITORY_ROOT,
  'tools/web-product-decisions.json'
);
const PRODUCT_IMPACT_MAP_SOURCE = readFileSync(PRODUCT_IMPACT_MAP_PATH, 'utf8');
const PRODUCT_DECISIONS_SOURCE = readFileSync(PRODUCT_DECISIONS_PATH, 'utf8');

const normalizePath = (path) => path.split(sep).join('/');
const relativeModulePath = (path) => normalizePath(relative(JAVASCRIPT_ROOT, path));

const walkJavaScriptFiles = (directory) => readdirSync(directory, { withFileTypes: true })
  .flatMap((entry) => {
    const path = join(directory, entry.name);
    if (entry.isDirectory()) return walkJavaScriptFiles(path);
    return entry.isFile() && entry.name.endsWith('.js') ? [path] : [];
  });

const PLAYWRIGHT_SPEC_SOURCES = walkJavaScriptFiles(
  join(REPOSITORY_ROOT, 'tests/web')
).filter((path) => path.endsWith('.playwright.spec.js'))
  .map((path) => readFileSync(path, 'utf8'));

const workflowJob = (jobId, workflow = TEST_WORKFLOW) => {
  const job = workflow.match(
    new RegExp(`\\n  ${jobId}:\\n[\\s\\S]*?(?=\\n  [a-z0-9-]+:\\n|$)`)
  )?.[0];
  assert.ok(job, `${jobId} job must exist`);
  return job;
};

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
  const detected = detectPrivilegedWebCapabilities(productionSources);
  const assertSubset = (actual, allowed, label) => {
    const allowedSet = new Set(allowed);
    assert.deepEqual(actual.filter((path) => !allowedSet.has(path)), [], label);
  };

  Object.entries(WEB_CHANGE_POLICY.allowedPrivilegedImporters).forEach(([target, allowed]) => {
    assertSubset(detected.importersByTarget[target], allowed, `${target} importers`);
  });
  assert.deepEqual(
    Object.keys(WEB_CHANGE_POLICY.allowedPrivilegedOwners).sort(),
    WEB_PRIVILEGED_CAPABILITY_KEYS
  );
  assert.deepEqual(
    Object.keys(WEB_CHANGE_POLICY.allowedPrivilegedImporters).sort(),
    WEB_PRIVILEGED_IMPORT_TARGETS
  );
  Object.entries(detected.operatorMatchesByCapability).forEach(([capability, matches_]) => {
    const actual = matches_.map(({ path }) => path);
    assertSubset(actual, WEB_CHANGE_POLICY.allowedPrivilegedOwners[capability], capability);
  });
});

test('shared privileged detectors preserve the characterized current-source facts', () => {
  assert.deepEqual(WEB_PRIVILEGED_CAPABILITY_KEYS, [
    'Canonical editor state',
    'Diagram Worker',
    'History',
    'Mounted SVG/Result replacement',
    'Python helper',
    'Render request',
    'Resource staging',
    'SVG/Result admission',
    'Session'
  ]);
  assert.deepEqual(WEB_PRIVILEGED_IMPORT_TARGETS, [
    'app/feature-editor.js',
    'app/legend.js',
    'app/preview-runtime.js',
    'app/python-helpers.js',
    'app/right-drawer.js',
    'services/config.js',
    'services/diagram-generation.js',
    'services/diagram-resource-staging.js',
    'services/diagram-worker-protocol.js',
    'services/gallery-session-migration.js',
    'services/history-snapshot.js',
    'services/history.js',
    'services/resource-payload-owner.js',
    'services/session-authority.js',
    'services/session-file.js',
    'services/session-request.js',
    'services/svg-result-ingestion.js',
    'services/svg-sanitization.js',
    'services/svg-serialization.js',
    'workers/diagram-generation-worker.js'
  ]);

  const detected = detectPrivilegedWebCapabilities(productionSources);
  assert.deepEqual(
    detected.importersByTarget,
    Object.fromEntries(
      WEB_PRIVILEGED_IMPORT_TARGETS.map((target) => [target, importersOf(target)])
    )
  );
  assert.deepEqual(detected.operatorMatchesByCapability, {
    'Canonical editor state': [
      { path: 'app/app-setup.js', count: 8 },
      { path: 'app/feature-editor.js', count: 1 },
      { path: 'app/legend.js', count: 1 },
      { path: 'app/preview-runtime.js', count: 1 },
      { path: 'app/right-drawer.js', count: 1 }
    ],
    'Diagram Worker': [
      { path: 'app/run-analysis.js', count: 1 },
      { path: 'services/diagram-generation.js', count: 1 },
      { path: 'services/losat.js', count: 2 },
      { path: 'workers/diagram-generation-worker.js', count: 1 },
      { path: 'workers/losat-threaded-worker.js', count: 2 }
    ],
    History: [
      { path: 'app/app-setup.js', count: 3 },
      { path: 'services/history-snapshot.js', count: 2 },
      { path: 'services/history.js', count: 7 }
    ],
    'Mounted SVG/Result replacement': [
      { path: 'app/app-setup.js', count: 1 },
      { path: 'app/feature-editor/color-actions.js', count: 1 },
      { path: 'app/feature-editor/label-actions.js', count: 2 },
      { path: 'app/feature-editor/svg-actions.js', count: 4 },
      { path: 'app/legend-layout/canvas-actions.js', count: 2 },
      { path: 'app/legend-layout/diagram-drag.js', count: 2 },
      { path: 'app/legend-layout/reposition-actions.js', count: 2 },
      { path: 'app/legend/drag-actions.js', count: 4 },
      { path: 'app/legend/entry-actions.js', count: 5 },
      { path: 'app/legend/sort-actions.js', count: 2 },
      { path: 'app/legend/stroke-actions.js', count: 2 },
      { path: 'app/preview-runtime.js', count: 4 },
      { path: 'app/results.js', count: 4 },
      { path: 'app/run-analysis.js', count: 2 },
      { path: 'app/svg-styles.js', count: 2 },
      { path: 'services/config.js', count: 5 },
      { path: 'services/history-snapshot.js', count: 1 },
      { path: 'services/svg-result-ingestion.js', count: 1 }
    ],
    'Python helper': [
      { path: 'app/python-helpers.js', count: 1 },
      { path: 'workers/diagram-generation-worker.js', count: 6 }
    ],
    'Render request': [
      { path: 'app/run-analysis.js', count: 1 },
      { path: 'services/config.js', count: 1 },
      { path: 'services/gallery-session-migration.js', count: 1 }
    ],
    'Resource staging': [
      { path: 'services/config.js', count: 3 },
      { path: 'services/diagram-generation.js', count: 2 },
      { path: 'services/diagram-resource-staging.js', count: 1 },
      { path: 'services/resource-payload-owner.js', count: 1 },
      { path: 'services/session-request.js', count: 2 },
      { path: 'workers/diagram-generation-worker.js', count: 2 }
    ],
    'SVG/Result admission': [
      { path: 'services/svg-result-ingestion.js', count: 7 },
      { path: 'services/svg-sanitization.js', count: 1 }
    ],
    Session: [
      { path: 'app/app-setup.js', count: 2 },
      { path: 'services/config.js', count: 5 },
      { path: 'services/gallery-session-migration.js', count: 3 },
      { path: 'services/session-file.js', count: 1 }
    ]
  });
});

test('versioned architecture detectors expose stable normalized subjects', () => {
  const fixture = new Map([
    [
      'services/session-request.js',
      'export const buildCanonicalRenderRequest = () => ({});\n'
        + '// export const buildCanonicalRenderRequest = () => fake();\n'
    ],
    [
      'app/run-analysis.js',
      "import { buildCanonicalRenderRequest } from '../services/session-request.js';\n"
        + 'const note = "buildCanonicalRenderRequest(); runDiagramGeneration();";\n'
        + 'export const generate = () => {\n'
        + '  const canonical = buildCanonicalRenderRequest();\n'
        + '  return runDiagramGeneration(canonical);\n'
        + '};\n'
    ],
    [
      'services/config.js',
      "import { buildCanonicalRenderRequest } from './session-request.js';\n"
        + 'export const save = () => buildCanonicalRenderRequest();\n'
    ],
    [
      'app/ignored.js',
      '// buildCanonicalRenderRequest(); runDiagramGeneration();\n'
        + 'export const note = "buildCanonicalRenderRequest";\n'
    ]
  ]);
  const semantic = WEB_ARCHITECTURE_DETECTORS['semantic-owner.render-request.v1'];
  const canonical = WEB_ARCHITECTURE_DETECTORS['canonical-path.render-request.v1'];

  assert.deepEqual(Object.keys(WEB_ARCHITECTURE_DETECTORS).sort(), [
    'canonical-path.render-request.v1',
    'semantic-owner.render-request.v1'
  ]);
  assert.equal(semantic.subjectCategory, 'definition-path');
  assert.deepEqual(semantic.detect(fixture), {
    definitionCount: 1,
    observedDefinitions: [{
      path: 'services/session-request.js',
      count: 1,
      subject: 'services/session-request.js'
    }],
    subjects: ['services/session-request.js']
  });
  assert.equal(
    semantic.encodeSubject({ path: 'gbdraw/web/js/services/session-request.js' }),
    'services/session-request.js'
  );

  assert.equal(canonical.subjectCategory, 'canonical-entry-edge');
  assert.deepEqual(canonical.detect(fixture), {
    observedEdges: [{
      from: 'app/run-analysis.js',
      to: 'services/session-request.js',
      subject: 'app/run-analysis.js -> services/session-request.js'
    }],
    subjects: ['app/run-analysis.js -> services/session-request.js']
  });
  assert.equal(
    canonical.encodeSubject({
      from: 'gbdraw/web/js/app/run-analysis.js',
      to: 'gbdraw/web/js/services/session-request.js'
    }),
    'app/run-analysis.js -> services/session-request.js'
  );
});

test('declared architecture rules actively govern the current source tree', () => {
  assert.deepEqual(validateArchitectureRuleRegistry(
    WEB_ARCHITECTURE_RULES,
    WEB_ARCHITECTURE_DETECTORS,
    { availableEnforcements: ['hard', 'report-only'] }
  ), { valid: true, errors: [] });

  const results = WEB_ARCHITECTURE_RULES.rules.map((rule) => {
    const detector = WEB_ARCHITECTURE_DETECTORS[rule.detector];
    const observed = classifyArchitectureRuleObservation(rule, detector.detect(productionSources));
    return {
      key: rule.key,
      ...evaluateArchitectureRuleResult({
        observation: observed.observation,
        mode: rule.enforcement.replace('-', '_').toUpperCase(),
        baselineRelation: 'NOT_APPLICABLE',
        authorityResolution: 'NOT_APPLICABLE'
      })
    };
  });
  assert.deepEqual(results.map(({ key, decision }) => ({ key, decision })), [
    { key: 'canonical-path.render-request', decision: 'PASS' },
    { key: 'semantic-owner.render-request', decision: 'PASS' }
  ]);
});

test('pure architecture evaluation owns no I/O, source detection, or authority data', () => {
  assert.doesNotMatch(WEB_ARCHITECTURE_EVALUATION_SOURCE, /^\s*import\s/m);
  assert.doesNotMatch(
    WEB_ARCHITECTURE_EVALUATION_SOURCE,
    /node:(?:fs|path|child_process|process)|process\.|console\.|(?:read|write|append)File|execFile|spawn/
  );
  assert.doesNotMatch(
    WEB_ARCHITECTURE_EVALUATION_SOURCE,
    /web-architecture-(?:detectors|rules|violations)|web-change-policy/
  );
  assert.doesNotMatch(WEB_ARCHITECTURE_EVALUATION_SOURCE, /^#!|process\.argv/);
});

test('pure Product Impact mechanics own no I/O, Git, environment, reporting, or CLI', () => {
  for (const source of [PRODUCT_IMPACT_EVALUATION_SOURCE, PRODUCT_IMPACT_DECISION_SOURCE]) {
    assert.doesNotMatch(source, /^\s*import\s/m);
    assert.doesNotMatch(
      source,
      /node:(?:fs|path|child_process|process)|process\.|console\.|(?:read|write|append)File|execFile|spawn|fetch\s*\(/
    );
    assert.doesNotMatch(source, /^#!|process\.argv/);
  }
  assert.doesNotMatch(
    PRODUCT_IMPACT_EVALUATION_SOURCE,
    /web-product-impact-(?:map|decisions)\.json|web-architecture-rules\.json/
  );
});

test('Product Impact authority bootstrap accepts zero or two inert JSON files', () => {
  const mapExists = existsSync(PRODUCT_IMPACT_MAP_PATH);
  const decisionsExist = existsSync(PRODUCT_DECISIONS_PATH);
  assert.equal(
    mapExists,
    decisionsExist,
    'Product Impact map and durable decision authority must be introduced together'
  );
  if (!mapExists) return;

  const map = JSON.parse(readFileSync(PRODUCT_IMPACT_MAP_PATH, 'utf8'));
  const decisions = JSON.parse(readFileSync(PRODUCT_DECISIONS_PATH, 'utf8'));
  assert.deepEqual(
    validateProductImpactMap(map, WEB_ARCHITECTURE_RULES, WEB_ARCHITECTURE_DETECTORS),
    { valid: true, errors: [] }
  );
  assert.deepEqual(validateProductDecisionAuthority(decisions, map), {
    valid: true,
    errors: []
  });
  const pilotEnforcement = map.ruleCoverage
    .filter(({ ruleKey }) => [
      'canonical-path.render-request',
      'semantic-owner.render-request'
    ].includes(ruleKey))
    .map(({ enforcement }) => enforcement);
  assert.equal(pilotEnforcement.length, 2);
  assert.equal(
    new Set(pilotEnforcement).size,
    1,
    'the two render-request pilot rules must activate in the same stage'
  );
  assert.deepEqual(decisions.decisions, []);
  map.concerns.flatMap(({ contracts }) => contracts).forEach(({ ref }) => {
    assert.equal(existsSync(join(REPOSITORY_ROOT, ref.split('::', 1)[0])), true, ref);
  });
});

test('registered architecture authority and executable matching have separate owners', () => {
  assert.deepEqual(Object.keys(WEB_ARCHITECTURE_RULES), ['schemaVersion', 'rules']);
  assert.doesNotMatch(
    WEB_CHANGE_BUDGET_SOURCE,
    /app\/run-analysis\.js|services\/session-request\.js/
  );
  assert.doesNotMatch(
    WEB_ARCHITECTURE_EVALUATION_SOURCE,
    /app\/run-analysis\.js|services\/session-request\.js/
  );
  assert.doesNotMatch(
    WEB_ARCHITECTURE_DETECTOR_SOURCE,
    /\ballowedEdges\b|\ballowedDefinitionPaths\b|\bexactActiveEdgeCount\b|\bexactDefinitionCount\b|\bbaselineEligible\b/
  );
  assert.match(WEB_ARCHITECTURE_DETECTOR_SOURCE, /RENDER_REQUEST_DEFINITION_PATTERN/);
  assert.match(WEB_ARCHITECTURE_DETECTOR_SOURCE, /RENDER_REQUEST_CALL_PATTERN/);
});

test('the Web checker remains the sole evaluator orchestrator and CLI entry point', () => {
  assert.match(
    WEB_CHANGE_BUDGET_SOURCE,
    /from '\.\/web-architecture-evaluation\.mjs'/
  );
  [
    'validateArchitectureRuleRegistry',
    'classifyArchitectureAuthorityDelta',
    'classifyArchitectureRuleObservation',
    'evaluateArchitectureRuleResult'
  ].forEach((name) => {
    assert.match(WEB_CHANGE_BUDGET_SOURCE, new RegExp(`\\b${name}\\s*\\(`));
  });
  assert.match(WEB_CHANGE_BUDGET_SOURCE, /^#!\/usr\/bin\/env node/);
  assert.match(WEB_CHANGE_BUDGET_SOURCE, /WEB_ARCHITECTURE_DETECTORS\[rule\.detector\]/);
  assert.doesNotMatch(WEB_CHANGE_BUDGET_SOURCE, /import\s*\(.*web-architecture-detectors/);
  assert.match(
    WEB_CHANGE_BUDGET_SOURCE,
    /from '\.\/web-product-impact-(?:evaluation|decision-source)\.mjs'/
  );
  [
    'createArchitectureSubjectDelta',
    'evaluateProductImpact',
    'parseCurrentProductImpactDecisions',
    'validateProductDecisionAuthority',
    'validateProductImpactMap'
  ].forEach((name) => {
    assert.match(WEB_CHANGE_BUDGET_SOURCE, new RegExp(`\\b${name}\\s*\\(`));
  });
  assert.match(WEB_CHANGE_BUDGET_SOURCE, /## Product impact/);
  assert.match(WEB_CHANGE_BUDGET_SOURCE, /NOT_AUTHORITATIVE_CANDIDATE/);
  assert.match(WEB_CHANGE_BUDGET_SOURCE, /candidate data does not alter this head runtime admission/);
});

test('report-only source facts preserve the characterized masking behavior', () => {
  const facts = detectReportOnlySourceFacts(
    "import helper from './helper.js';\n"
      + '// export const createCommentManager = () => watch(commentCache);\n'
      + 'const note = "migrateLegacyString watch(stringCache)";\n'
      + 'export const createFixtureManager = () => {\n'
      + '  const activeCache = ref(0);\n'
      + '  watch(activeCache, () => migrateLegacyFixture());\n'
      + '  return { activeCache: activeCache };\n'
      + '};\n'
  );

  assert.deepEqual(facts, {
    exportedNames: ['createFixtureManager'],
    declaredNames: ['note', 'createFixtureManager', 'activeCache'],
    reactiveDeclarations: ['activeCache (ref)'],
    watcherCount: 1,
    objectKeys: ['activeCache'],
    compatibilityNames: ['migrateLegacyFixture'],
    namedResourceNames: ['createFixtureManager', 'activeCache'],
    importSpecifiers: ['./helper.js']
  });
});

test('workflow triggers separate dev admission, dev staging, promotion, and deployment', () => {
  const candidateActivityTypes = /types: \[opened, synchronize, reopened, labeled, unlabeled\]/;
  const trustedActivityTypes = /types: \[opened, synchronize, reopened, edited, labeled, unlabeled\]/;
  const triggerBranches = (workflow, trigger) => {
    const match = workflow.match(new RegExp(
      `\\r?\\n  ${trigger}:\\r?\\n    branches: (\\[[^\\r\\n]+\\])`
    ));
    assert.ok(match, `${trigger} branch filter must exist`);
    return JSON.parse(match[1]);
  };

  assert.match(TEST_WORKFLOW, /\n  pull_request:\n/);
  assert.deepEqual(triggerBranches(TEST_WORKFLOW, 'pull_request'), ['dev']);
  assert.deepEqual(triggerBranches(TEST_WORKFLOW, 'push'), ['dev']);
  assert.match(TEST_WORKFLOW, candidateActivityTypes);
  assert.doesNotMatch(TEST_WORKFLOW, /types: \[[^\]\r\n]*\bedited\b/);
  assert.match(TEST_WORKFLOW, /name: Run core tests/);
  assert.match(
    TEST_WORKFLOW,
    /group: tests-\$\{\{ github\.event_name \}\}-\$\{\{ github\.event\.pull_request\.number \|\| github\.ref \}\}/
  );

  assert.match(BASE_POLICY_WORKFLOW, /\n  pull_request_target:\n/);
  assert.deepEqual(
    triggerBranches(BASE_POLICY_WORKFLOW, 'pull_request_target'),
    ['dev', 'main']
  );
  assert.match(BASE_POLICY_WORKFLOW, trustedActivityTypes);
  assert.match(
    BASE_POLICY_WORKFLOW,
    /permissions:\n  contents: read\n  pull-requests: read\n  actions: read/
  );
  assert.match(
    BASE_POLICY_WORKFLOW,
    /group: web-base-policy-\$\{\{ github\.event\.pull_request\.number \}\}/
  );
  assert.deepEqual(
    WORKFLOW_NAMES.filter((name, index) => (
      /types: \[[^\]\r\n]*\bedited\b/.test(WORKFLOW_SOURCES[index])
    )),
    ['web-base-policy.yml']
  );

  assert.deepEqual(triggerBranches(GALLERY_PUBLICATION_WORKFLOW, 'push'), ['dev']);
  assert.doesNotMatch(GALLERY_PUBLICATION_WORKFLOW, /\n  pull_request:\n/);
  assert.match(GALLERY_PUBLICATION_WORKFLOW, /\n  workflow_dispatch:\n/);

  assert.deepEqual(triggerBranches(DEPLOY_WORKFLOW, 'push'), ['main']);
  assert.doesNotMatch(TEST_WORKFLOW, /"main"|"refactoring"/);
  assert.doesNotMatch(GALLERY_PUBLICATION_WORKFLOW, /"main"|"refactoring"/);
});

test('PR-to-dev jobs and aggregate use the trusted selective plan', () => {
  const planner = workflowJob('ci-impact');
  assert.match(planner, /Checkout complete history[\s\S]*fetch-depth: 0/);
  assert.match(planner, /node --test tests\/ci\/\*\.test\.mjs/);
  assert.match(
    planner,
    /ref: \$\{\{ github\.event\.pull_request\.base\.sha \}\}[\s\S]*path: \.ci-trusted-base[\s\S]*persist-credentials: false/
  );
  assert.match(planner, /CI_IMPACT_REPOSITORY_ROOT: \$\{\{ github\.workspace \}\}/);
  assert.match(planner, /node \.ci-trusted-base\/tools\/ci-impact\.mjs plan/);
  assert.match(planner, /Build shadow dev CI impact plan[\s\S]*node tools\/ci-impact\.mjs plan/);
  assert.doesNotMatch(TEST_WORKFLOW, /\n    paths(?:-ignore)?:/);

  const corePr = workflowJob('core-pr');
  assert.match(corePr, /\n    name: Core PR \(Python 3\.11\)\n/);
  assert.match(corePr, /needs: ci-impact/);
  assert.match(corePr, /needs\.ci-impact\.result == 'success'/);
  assert.match(corePr, /requiredJobs, 'core-pr'/);
  assert.match(corePr, /-m "not slow and not \(recipe or gallery or browser\)"/);
  assert.match(corePr, /git diff --exit-code -- tests\/reference_outputs\//);
  assert.doesNotMatch(corePr, /matrix:/);

  const webPrSmoke = workflowJob('web-pr-smoke');
  assert.match(webPrSmoke, /\n    name: Web PR smoke\n/);
  assert.match(webPrSmoke, /needs: ci-impact/);
  assert.match(webPrSmoke, /needs\.ci-impact\.result == 'success'/);
  assert.match(webPrSmoke, /requiredJobs, 'web-pr-smoke'/);
  assert.match(webPrSmoke, /timeout-minutes: 5/);
  assert.equal([...webPrSmoke.matchAll(/uses: actions\/checkout@/g)].length, 1);
  assert.equal([...webPrSmoke.matchAll(/npm ci/g)].length, 1);
  assert.equal([...webPrSmoke.matchAll(/playwright install --with-deps chromium/g)].length, 1);
  assert.equal([...webPrSmoke.matchAll(/python tools\/prepare_browser_wheel\.py/g)].length, 1);
  assert.match(webPrSmoke, /Run fast Web JavaScript contracts/);
  assert.match(webPrSmoke, /! -name 'gallery-session-publication\.test\.mjs'/);
  assert.match(webPrSmoke, /-m "browser and not slow"/);
  assert.match(webPrSmoke, /npm run test:web:pr-smoke/);
  assert.match(webPrSmoke, /if: failure\(\)[\s\S]+path: test-results\//);
  assert.doesNotMatch(
    webPrSmoke,
    /functional-full|perf-smoke|losat-cache|gallery-publication|Vibrio/
  );

  const gate = workflowJob('pr-gate');
  assert.match(gate, /\n    name: PR \/ gate\n/);
  const dependencies = gate.match(
    /\n    needs:\n((?:      - [a-z0-9-]+\n)+)/
  )?.[1].trim().split('\n').map((line) => line.replace('- ', '').trim());
  assert.deepEqual(dependencies, [
    'ci-impact',
    'web-change-budget',
    'core-pr',
    'recipes-standard',
    'gallery',
    'lint',
    'web-pr-smoke'
  ]);
  assert.match(
    gate,
    /if: >-\n      always\(\) &&\n      github\.event_name == 'pull_request' &&\n      github\.base_ref == 'dev'/
  );
  assert.match(
    gate,
    /ref: \$\{\{ github\.event\.pull_request\.base\.sha \}\}[\s\S]*path: \.ci-trusted-base[\s\S]*persist-credentials: false/
  );
  assert.match(gate, /CI_IMPACT_PLAN_JSON: \$\{\{ needs\.ci-impact\.outputs\.plan \}\}/);
  assert.match(gate, /CI_IMPACT_NEEDS_JSON: \$\{\{ toJSON\(needs\) \}\}/);
  assert.match(gate, /CI_IMPACT_EXPECTED_PROFILE: pr/);
  assert.match(gate, /CI_IMPACT_EXPECTED_WORKFLOW_SHA: \$\{\{ github\.sha \}\}/);
  assert.match(gate, /node \.ci-trusted-base\/tools\/ci-impact\.mjs gate/);
  assert.doesNotMatch(gate, /test "\$\{\{ needs\./);
  assert.doesNotMatch(gate, /gh api|curl|check-web-change-budget|pull_request_target/);

  for (const jobId of [
    'web-change-budget',
    'core-pr',
    'recipes-standard',
    'gallery',
    'lint',
    'web-pr-smoke'
  ]) {
    const job = workflowJob(jobId);
    assert.match(job, /needs: ci-impact/, `${jobId} must wait for the plan`);
    assert.match(
      job,
      new RegExp(`requiredJobs, '${jobId}'`),
      `${jobId} must read its selection from the plan`
    );
    if (['web-change-budget', 'recipes-standard', 'gallery', 'lint'].includes(jobId)) {
      assert.match(job, /!cancelled\(\)/, `${jobId} must remain runnable on dev`);
      assert.match(job, /github\.event_name == 'push' && github\.ref == 'refs\/heads\/dev'/);
      assert.match(
        job,
        /github\.event_name == 'workflow_dispatch' && github\.ref == 'refs\/heads\/dev'/
      );
    }
  }
  for (const jobId of [
    'core',
    'browser',
    'playwright-functional',
    'playwright-performance',
    'acceptance-supported-main',
    'slow-main',
    'losat-cache-browser-acceptance'
  ]) {
    const job = workflowJob(jobId);
    assert.doesNotMatch(
      job,
      /github\.event_name == 'pull_request' && github\.base_ref == 'dev'/,
      `${jobId} must not run on a dev pull request`
    );
  }
});

test('PR smoke selection is explicit while the full functional inventory stays wide', () => {
  assert.equal(
    PACKAGE_SCRIPTS['test:web:functional-full'],
    'playwright test --config=playwright.functional.config.js'
  );
  assert.equal(
    PACKAGE_SCRIPTS['test:web:pr-smoke'],
    'playwright test --config=playwright.pr-smoke.config.js'
  );
  assert.equal(
    PACKAGE_SCRIPTS['test:web:functional-smoke'],
    'npm run test:web:functional-full'
  );
  assert.doesNotMatch(FUNCTIONAL_PLAYWRIGHT_CONFIG, /\bgrep\s*:/);
  assert.match(FUNCTIONAL_PLAYWRIGHT_CONFIG, /performance\\\.playwright\\\.spec\\\.js/);
  assert.match(FUNCTIONAL_PLAYWRIGHT_CONFIG, /losat-cache-migration/);

  assert.match(PR_SMOKE_PLAYWRIGHT_CONFIG, /grep: \/@pr-smoke\//);
  assert.match(PR_SMOKE_PLAYWRIGHT_CONFIG, /retries: 0/);
  assert.match(PR_SMOKE_PLAYWRIGHT_CONFIG, /workers: 1/);
  assert.match(PR_SMOKE_PLAYWRIGHT_CONFIG, /reporter: 'list'/);
  assert.match(PR_SMOKE_PLAYWRIGHT_CONFIG, /trace: 'retain-on-failure'/);
  assert.match(PR_SMOKE_PLAYWRIGHT_CONFIG, /name === 'chromium'/);
  assert.match(PR_SMOKE_PLAYWRIGHT_CONFIG, /performance\\\.playwright\\\.spec\\\.js/);
  assert.match(PR_SMOKE_PLAYWRIGHT_CONFIG, /losat-cache-migration/);
  assert.doesNotMatch(
    `${PACKAGE_SCRIPTS['test:web:pr-smoke']}\n${PR_SMOKE_PLAYWRIGHT_CONFIG}`,
    /pass-with-no-tests/
  );

  const selectedCount = PLAYWRIGHT_SPEC_SOURCES.reduce(
    (count, source) => count + [...source.matchAll(/tag: ['"]@pr-smoke['"]/g)].length,
    0
  );
  assert.ok(selectedCount >= 6 && selectedCount <= 10, `selected ${selectedCount} smoke tests`);
});

test('exact dev staging runs the full inventory with four mandatory Playwright shards', () => {
  const stagingOnly = [
    'core',
    'browser',
    'playwright-functional',
    'playwright-performance',
    'acceptance-supported-main',
    'slow-main',
    'losat-cache-browser-acceptance'
  ];
  for (const jobId of stagingOnly) {
    const job = workflowJob(jobId);
    assert.match(job, /github\.event_name == 'push' && github\.ref == 'refs\/heads\/dev'/);
    assert.match(
      job,
      /github\.event_name == 'workflow_dispatch' && github\.ref == 'refs\/heads\/dev'/
    );
  }

  const fullPlaywright = workflowJob('playwright-functional');
  assert.match(fullPlaywright, /\n    name: Playwright functional \(shard \$\{\{ matrix\.shard \}\}\/4\)\n/);
  assert.match(fullPlaywright, /shard: \[1, 2, 3, 4\]/);
  assert.match(
    fullPlaywright,
    /npm run test:web:functional-full -- --shard=\$\{\{ matrix\.shard \}\}\/4/
  );
  assert.match(fullPlaywright, /playwright-functional-shard-\$\{\{ matrix\.shard \}\}-traces-/);
  assert.match(BASE_PLAYWRIGHT_CONFIG, /retries: process\.env\.CI \? 2 : 0/);
  assert.match(BASE_PLAYWRIGHT_CONFIG, /workers: process\.env\.CI \? 1 : undefined/);

  const gate = workflowJob('dev-staging-gate');
  const dependencies = gate.match(
    /\n    needs:\n((?:      - [a-z0-9-]+\n)+)/
  )?.[1].trim().split('\n').map((line) => line.replace('- ', '').trim());
  assert.deepEqual(dependencies, [
    'ci-impact',
    'web-change-budget',
    'core',
    'recipes-standard',
    'gallery',
    'browser',
    'playwright-functional',
    'playwright-performance',
    'acceptance-supported-main',
    'slow-main',
    'lint',
    'losat-cache-browser-acceptance'
  ]);
  assert.match(gate, /\n    name: Dev staging \/ gate\n/);
  assert.match(gate, /always\(\)/);
  assert.match(gate, /github\.event_name == 'push' && github\.ref == 'refs\/heads\/dev'/);
  assert.match(
    gate,
    /github\.event_name == 'workflow_dispatch' && github\.ref == 'refs\/heads\/dev'/
  );
  dependencies.forEach((jobId) => {
    assert.ok(
      gate.includes(`test "\${{ needs.${jobId}.result }}" = "success"`),
      `${jobId} must fail the staging aggregate unless it succeeds`
    );
  });
  assert.doesNotMatch(gate, /gh api|curl|check-promotion-readiness/);
});

test('Gallery readiness aggregates common, Vibrio, and projection evidence on dev', () => {
  const browser = workflowJob('browser', GALLERY_PUBLICATION_WORKFLOW);
  assert.match(browser, /example: common 9\n            command: test:web:gallery-publication/);
  assert.match(browser, /example: Vibrio\n            command: test:web:vibrio-generate/);
  assert.match(browser, /ref: \$\{\{ github\.sha \}\}/);

  const performance = workflowJob('performance', GALLERY_PUBLICATION_WORKFLOW);
  assert.match(performance, /\n    name: Gallery publication performance \(projection\)\n/);
  assert.match(performance, /measure_gallery_publication_performance\.py projection/);
  assert.match(performance, /github\.event_name == 'workflow_dispatch' && inputs\.complete_refresh/);
  assert.match(performance, /measure_gallery_publication_performance\.py refresh/);

  const gate = workflowJob('readiness-gate', GALLERY_PUBLICATION_WORKFLOW);
  assert.match(gate, /\n    name: Gallery readiness \/ gate\n/);
  assert.match(gate, /\n    needs:\n      - browser\n      - performance\n/);
  assert.match(gate, /\n    if: always\(\)\n/);
  for (const jobId of ['browser', 'performance']) {
    assert.ok(
      gate.includes(`test "\${{ needs.${jobId}.result }}" = "success"`),
      `${jobId} must fail the Gallery aggregate unless all matrix jobs succeed`
    );
  }
  assert.doesNotMatch(gate, /gh api|curl|check-promotion-readiness/);
});

test('trusted promotion uses base code for topology, exact-SHA evidence, and tree proof', () => {
  const devPolicy = workflowJob('trusted-base-policy', BASE_POLICY_WORKFLOW);
  assert.match(devPolicy, /\n    name: Web base policy \(trusted base\)\n/);
  assert.match(devPolicy, /if: github\.event\.pull_request\.base\.ref == 'dev'/);
  assert.equal([...devPolicy.matchAll(/uses: actions\/checkout@/g)].length, 1);
  assert.match(devPolicy, /ref: \$\{\{ github\.event\.pull_request\.base\.sha \}\}/);
  assert.match(devPolicy, /fetch-depth: 0/);
  assert.match(devPolicy, /persist-credentials: false/);
  assert.match(devPolicy, /git fetch --no-tags origin "\$HEAD_SHA"/);
  assert.doesNotMatch(
    devPolicy,
    /ref:.*pull_request\.head|git (?:checkout|switch) |uses: \.\//
  );

  const promotion = workflowJob('promotion-gate', BASE_POLICY_WORKFLOW);
  assert.match(promotion, /\n    name: Promotion \/ gate\n/);
  assert.match(promotion, /if: github\.event\.pull_request\.base\.ref == 'main'/);
  assert.equal([...promotion.matchAll(/uses: actions\/checkout@/g)].length, 1);
  assert.match(promotion, /ref: \$\{\{ github\.event\.pull_request\.base\.sha \}\}/);
  assert.match(promotion, /fetch-depth: 0/);
  assert.match(promotion, /persist-credentials: false/);
  assert.match(promotion, /EXPECTED_REPOSITORY: satoshikawato\/gbdraw/);
  assert.match(promotion, /test "\$BASE_REF" = "main"/);
  assert.match(promotion, /test "\$HEAD_REF" = "dev"/);
  assert.match(promotion, /test "\$HEAD_REPOSITORY" = "\$EXPECTED_REPOSITORY"/);
  assert.match(
    promotion,
    /git fetch --no-tags origin "refs\/heads\/dev:refs\/remotes\/origin\/dev"/
  );
  assert.match(promotion, /test "\$HEAD_SHA" = "\$REMOTE_DEV_SHA"/);
  assert.match(
    promotion,
    /node tools\/check-web-change-budget\.mjs --base "\$BASE_SHA" --head "\$HEAD_SHA"/
  );
  assert.equal([...promotion.matchAll(/workflow-evidence/g)].length, 2);
  assert.match(
    promotion,
    /--workflow-path \.github\/workflows\/test\.yml[\s\S]*--aggregate-name "Dev staging \/ gate"/
  );
  assert.match(
    promotion,
    /--workflow-path \.github\/workflows\/gallery-publication\.yml[\s\S]*--aggregate-name "Gallery readiness \/ gate"/
  );
  assert.match(promotion, /check-promotion-readiness\.mjs merge-tree/);
  assert.match(promotion, /--base-sha "\$BASE_SHA"[\s\S]*--head-sha "\$HEAD_SHA"/);
  assert.doesNotMatch(
    promotion,
    /ref:.*pull_request\.head|git (?:checkout|switch) |npm ci|pip install|functional-full|perf-smoke|pytest/
  );
  assert.doesNotMatch(promotion, /uses: \.\//);
});

test('steady-state topology removes legacy main producers without changing final owners', () => {
  assert.doesNotMatch(TEST_WORKFLOW, /\n  playwright-functional-main-transition:\n/);
  assert.doesNotMatch(TEST_WORKFLOW, /main-transition|PR 7 removes/);
  assert.doesNotMatch(
    TEST_WORKFLOW,
    /github\.base_ref == 'main'|github\.head_ref == 'dev'|head\.repo\.full_name == github\.repository/
  );
  assert.doesNotMatch(BASE_POLICY_WORKFLOW, /trusted-base-policy-main-transition|PR 7 removes/);
  assert.equal(
    [...BASE_POLICY_WORKFLOW.matchAll(/name: Web base policy \(trusted base\)/g)].length,
    1
  );
  assert.equal([...BASE_POLICY_WORKFLOW.matchAll(/name: Promotion \/ gate/g)].length, 1);
  assert.deepEqual(WORKFLOW_NAMES, [
    'deploy_web.yml',
    'gallery-publication.yml',
    'test.yml',
    'web-base-policy.yml'
  ]);
  assert.doesNotMatch(TEST_WORKFLOW, /needs: web-change-budget/);
  assert.equal(
    existsSync(join(REPOSITORY_ROOT, '.github/workflows/dev-staging.yml')),
    false
  );
  assert.equal(
    existsSync(join(REPOSITORY_ROOT, '.github/workflows/gallery-readiness.yml')),
    false
  );
  assert.doesNotMatch(TEST_WORKFLOW, /dev-staging\.yml|gallery-readiness\.yml/);
  assert.doesNotMatch(
    WORKFLOW_SOURCES.join('\n'),
    /branches\/[^"]+\/protection|required_status_checks|\/rulesets(?:\b|\/)/
  );
});

test('the PR template retains architecture evidence anchors', () => {
  [
    'Architecture impact',
    'This is not architecture-bearing',
    'This is architecture-bearing',
    'Semantic owners before',
    'Semantic owners after',
    'Canonical production paths before',
    'Canonical production paths after',
    'Maintainer architecture decision comment permalink'
  ].forEach((anchor) => {
    assert.ok(PULL_REQUEST_TEMPLATE.includes(anchor), `Missing PR template anchor: ${anchor}`);
  });
});

test('the PR template contains one optional bounded Product Impact decision block', () => {
  const startMarkers = PULL_REQUEST_TEMPLATE.match(
    /<!-- gbdraw-product-impact-decision:start -->/g
  ) || [];
  const endMarkers = PULL_REQUEST_TEMPLATE.match(
    /<!-- gbdraw-product-impact-decision:end -->/g
  ) || [];
  assert.equal(startMarkers.length, 1);
  assert.equal(endMarkers.length, 1);
  assert.ok(
    PULL_REQUEST_TEMPLATE.indexOf(startMarkers[0])
      < PULL_REQUEST_TEMPLATE.indexOf(endMarkers[0])
  );
  assert.match(
    PULL_REQUEST_TEMPLATE,
    /\{"schemaVersion":1,"headSha":"","decisions":\[\]\}/
  );
});

test('the PR template requires exact changed-scope architecture debt arithmetic', () => {
  [
    'Owner-excess rows (repeat for every changed capability)',
    'O before',
    'T before',
    'OE before',
    'O after',
    'T after',
    'OE after',
    'Path-excess rows (repeat for every changed behavior)',
    'P before',
    'PE before',
    'P after',
    'PE after',
    'Compatibility-burden rows (repeat for every changed compatibility namespace)',
    'Compatibility path stable IDs before',
    'CB before',
    'Compatibility path stable IDs after',
    'CB after',
    'Changed-scope totals',
    'OE before -> OE after; delta(OE) = <integer>',
    'PE before -> PE after; delta(PE) = <integer>',
    'CB before -> CB after; delta(CB) = <integer>',
    'Superseded semantic owners',
    'Superseded canonical production paths',
    'Superseded compatibility paths, with stable IDs',
    '`<= 0`, `non-positive`, or `yes` alone is insufficient evidence.'
  ].forEach((anchor) => {
    assert.ok(PULL_REQUEST_TEMPLATE.includes(anchor), `Missing exact-debt anchor: ${anchor}`);
  });
});

test('normal CI uses only dev as the staging push branch', () => {
  assert.match(
    TEST_WORKFLOW,
    /push:\n    branches: \["dev"\]/
  );
  assert.match(TEST_WORKFLOW, /\n  workflow_dispatch:\n/);
});

test('dev staging Web checks scope only the newly integrated change', () => {
  const step = TEST_WORKFLOW.match(
    /      - name: Check Web change budget\n[\s\S]*?(?=\n      - name: Check Web architecture contracts)/
  )?.[0];
  assert.ok(step, 'Web change-budget step must exist');

  assert.match(
    step,
    /DEV_INTEGRATED_CHANGE: \$\{\{ github\.ref == 'refs\/heads\/dev' && \(github\.event_name == 'push' \|\| github\.event_name == 'workflow_dispatch'\) \}\}/
  );
  assert.match(
    step,
    /DEV_PUSH_BASE: \$\{\{ github\.event_name == 'push' && github\.event\.before \|\| '' \}\}/
  );
  assert.match(
    step,
    /DEV_PUSH_HEAD: \$\{\{ github\.event_name == 'push' && github\.sha \|\| '' \}\}/
  );
  assert.match(step, /BASE_SHA="\$DEV_PUSH_BASE"/);
  assert.match(step, /HEAD_SHA="\$DEV_PUSH_HEAD"/);
  assert.match(step, /BASE_SHA="\$\(git rev-parse HEAD\^1\)"/);
  assert.match(step, /HEAD_SHA="\$\(git rev-parse HEAD\)"/);
  assert.doesNotMatch(step, /git fetch origin main|git merge-base origin\/main HEAD/);
  assert.match(
    step,
    /node tools\/check-web-change-budget\.mjs --base "\$BASE_SHA" --head "\$HEAD_SHA"/
  );
  assert.match(
    step,
    /WEB_CHANGE_BASE: \$\{\{ github\.event_name == 'pull_request' && github\.event\.pull_request\.base\.sha \|\| '' \}\}/
  );
  assert.match(
    step,
    /WEB_CHANGE_HEAD: \$\{\{ github\.event_name == 'pull_request' && github\.event\.pull_request\.head\.sha \|\| '' \}\}/
  );
  assert.match(
    step,
    /else\n            node tools\/check-web-change-budget\.mjs\n          fi/
  );
});

test('supported-version and slow matrices cover exact dev staging', () => {
  const stagingCondition = [
    "    if: >-",
    "      (github.event_name == 'push' && github.ref == 'refs/heads/dev') ||",
    "      (github.event_name == 'workflow_dispatch' && github.ref == 'refs/heads/dev')"
  ].join('\n');

  for (const jobName of ['acceptance-supported-main', 'slow-main']) {
    const job = TEST_WORKFLOW.match(
      new RegExp(`\\n  ${jobName}:\\n[\\s\\S]*?(?=\\n  [a-z0-9-]+:\\n|$)`)
    )?.[0];
    assert.ok(job, `${jobName} job must exist`);
    assert.ok(job.includes(stagingCondition), `${jobName} must use the exact dev staging condition`);
  }
});

const CHANGE_BUDGET_CHECKER = join(REPOSITORY_ROOT, 'tools/check-web-change-budget.mjs');
const WEB_ARCHITECTURE_DETECTOR_MODULE = join(
  REPOSITORY_ROOT,
  'tools/web-architecture-detectors.mjs'
);
const WEB_ARCHITECTURE_EVALUATION_MODULE = join(
  REPOSITORY_ROOT,
  'tools/web-architecture-evaluation.mjs'
);
const PRODUCT_IMPACT_EVALUATION_MODULE = join(
  REPOSITORY_ROOT,
  'tools/web-product-impact-evaluation.mjs'
);
const PRODUCT_IMPACT_DECISION_MODULE = join(
  REPOSITORY_ROOT,
  'tools/web-product-impact-decision-source.mjs'
);
const WEB_ARCHITECTURE_RATCHET_FIXTURES = join(
  REPOSITORY_ROOT,
  'tests/web/architecture-ratchet-fixtures.test.mjs'
);
const PRODUCT_IMPACT_RATCHET_FIXTURES = join(
  REPOSITORY_ROOT,
  'tests/web/product-impact-ratchet-fixtures.test.mjs'
);
const WEB_PROMOTION_CONTEXT_MODULE = join(
  REPOSITORY_ROOT,
  'tools/web-promotion-context.mjs'
);
const PROMOTION_READINESS_HELPER = join(
  REPOSITORY_ROOT,
  'tools/check-promotion-readiness.mjs'
);
const PROMOTION_READINESS_TEST = join(
  REPOSITORY_ROOT,
  'tests/web/promotion-readiness.test.mjs'
);
const BUDGET_POLICY = Object.freeze({
  allowedPrivilegedImporters: {
    'services/diagram-generation.js': ['app/editor.js', 'app/run-analysis.js'],
    'services/session-request.js': ['app/run-analysis.js'],
    'workers/diagram-generation-worker.js': []
  },
  allowedPrivilegedOwners: {
    'Diagram Worker': [
      'app/editor.js',
      'app/future-owner.js',
      'app/run-analysis.js',
      'services/diagram-generation.js'
    ],
    'Render request': ['app/run-analysis.js']
  }
});
const BUDGET_PROFILES = Object.freeze([
  Object.freeze({
    name: 'ordinary',
    environment: Object.freeze({}),
    productionFiles: 8,
    grossChurn: 800,
    netAdditions: 100
  }),
  Object.freeze({
    name: 'architecture',
    environment: Object.freeze({ WEB_ARCHITECTURE_CHANGE: 'true' }),
    productionFiles: 12,
    grossChurn: 1500,
    netAdditions: 400
  })
]);
const BUDGET_LINE_CAPACITY = 800;
const budgetLineSource = (changedLines = 0, appendedLines = 0) => [
  ...Array.from(
    { length: BUDGET_LINE_CAPACITY },
    (_, index) => `const budgetLine${index} = ${index < changedLines ? 1 : 0};`
  ),
  ...Array.from(
    { length: appendedLines },
    (_, index) => `const appendedBudgetLine${index} = 1;`
  )
].join('\n') + '\n';
const BUDGET_FIXTURE = Object.freeze({
  '.github/workflows/deploy_web.yml': 'name: baseline deploy\n',
  '.github/workflows/gallery-publication.yml': 'name: baseline gallery\n',
  '.github/workflows/test.yml': 'name: baseline\n',
  '.github/workflows/web-base-policy.yml': 'name: baseline policy\n',
  'docs/internal/WEB_CHANGE_POLICY.md': '# Baseline Web change policy\n',
  'package.json': '{"private":true}\n',
  'tests/web/architecture-contracts.test.mjs': '// baseline contract\n',
  'tests/web/contracts/session-regenerate-intent.playwright.spec.js': (
    "test('no-draft session preserves preview, active config, canonical request, and catalog', () => {});\n"
  ),
  'tests/web/session-request.test.mjs': (
    "assert.deepEqual(roundTripCanonical.renderRequest, canonical.renderRequest,\n"
      + "  'canonical Session projection rebuilds the same render request');\n"
  ),
  'tests/web/promotion-readiness.test.mjs': readFileSync(
    PROMOTION_READINESS_TEST,
    'utf8'
  ),
  'tools/check-promotion-readiness.mjs': readFileSync(
    PROMOTION_READINESS_HELPER,
    'utf8'
  ),
  'tools/check-web-change-budget.mjs': readFileSync(CHANGE_BUDGET_CHECKER, 'utf8'),
  'tools/web-architecture-detectors.mjs': readFileSync(
    WEB_ARCHITECTURE_DETECTOR_MODULE,
    'utf8'
  ),
  'tools/web-architecture-evaluation.mjs': readFileSync(
    WEB_ARCHITECTURE_EVALUATION_MODULE,
    'utf8'
  ),
  'tools/web-product-impact-evaluation.mjs': readFileSync(
    PRODUCT_IMPACT_EVALUATION_MODULE,
    'utf8'
  ),
  'tools/web-product-impact-decision-source.mjs': readFileSync(
    PRODUCT_IMPACT_DECISION_MODULE,
    'utf8'
  ),
  'tools/web-product-impact-map.json': PRODUCT_IMPACT_MAP_SOURCE,
  'tools/web-product-decisions.json': PRODUCT_DECISIONS_SOURCE,
  'tools/web-architecture-rules.json': WEB_ARCHITECTURE_RULES_SOURCE,
  'tools/web-change-source.mjs': readFileSync(
    join(REPOSITORY_ROOT, 'tools/web-change-source.mjs'),
    'utf8'
  ),
  'tools/web-promotion-context.mjs': readFileSync(
    WEB_PROMOTION_CONTEXT_MODULE,
    'utf8'
  ),
  'tools/web-change-policy.json': `${JSON.stringify(BUDGET_POLICY, null, 2)}\n`,
  'tests/web/architecture-ratchet-fixtures.test.mjs': readFileSync(
    WEB_ARCHITECTURE_RATCHET_FIXTURES,
    'utf8'
  ),
  'tests/web/product-impact-ratchet-fixtures.test.mjs': readFileSync(
    PRODUCT_IMPACT_RATCHET_FIXTURES,
    'utf8'
  ),
  'gbdraw/web/index.html': '<main>baseline</main>\n',
  'gbdraw/web/js/app/editor.js': (
    "import { runDiagramGeneration } from '../services/diagram-generation.js';\n"
    + 'export const editExistingOwner = () => runDiagramGeneration();\n'
  ),
  'gbdraw/web/js/app/future-owner.js': 'export const futureOwnerPlaceholder = () => 1;\n',
  'gbdraw/web/js/app/run-analysis.js': (
    "import { runDiagramGeneration } from '../services/diagram-generation.js';\n"
    + "import { buildCanonicalRenderRequest } from '../services/session-request.js';\n"
    + 'export const run = () => {\n'
    + '  const request = buildCanonicalRenderRequest();\n'
    + '  return runDiagramGeneration(request);\n'
    + '};\n'
  ),
  'gbdraw/web/js/app/secondary.js': 'export const editSecondaryOwner = () => 1;\n',
  'gbdraw/web/js/app/budget-lines.js': budgetLineSource(),
  'gbdraw/web/js/services/diagram-generation.js': (
    'export const runDiagramGeneration = () => 1;\n'
  ),
  'gbdraw/web/js/services/session-file.js': (
    'export const readSession = () => ({ version: 1 });\n'
  ),
  'gbdraw/web/js/services/session-request.js': (
    'export const buildCanonicalRenderRequest = () => ({});\n'
  ),
  'gbdraw/web/vendor/library.js': 'globalThis.vendorLibrary = true;\n'
});
const FUTURE_GUARD_PATHS = Object.freeze([
  'docs/internal/ARCHITECTURE_FITNESS_FUNCTION_RATCHET.md',
  '.github/pull_request_template.md',
  'tools/web-architecture-violations.json'
]);
const PROTECTED_ARCHITECTURE_GUARD_PATHS = Object.freeze([
  ...FUTURE_GUARD_PATHS,
  'tools/web-architecture-rules.json'
]);
const FUTURE_AUTHORITY_PATHS = Object.freeze([
  'docs/internal/ARCHITECTURE_FITNESS_FUNCTION_RATCHET.md',
  'tools/web-architecture-rules.json',
  'tools/web-architecture-violations.json'
]);
const PRODUCT_IMPACT_GUARD_PATHS = Object.freeze([
  'docs/internal/PRODUCT_IMPACT_RATCHET.md',
  'tools/web-product-impact-map.json',
  'tools/web-product-decisions.json',
  'tools/web-product-impact-evaluation.mjs',
  'tools/web-product-impact-decision-source.mjs',
  'tests/web/product-impact-ratchet-fixtures.test.mjs'
]);
const PRODUCT_IMPACT_CHECKER_IMPLEMENTATION_PATHS = Object.freeze([
  'tools/web-product-impact-evaluation.mjs',
  'tools/web-product-impact-decision-source.mjs'
]);
const PRODUCT_IMPACT_AUTHORITY_PATHS = Object.freeze([
  'docs/internal/PRODUCT_IMPACT_RATCHET.md',
  'tools/web-product-impact-map.json',
  'tools/web-product-decisions.json'
]);
const NARROW_PRODUCT_IMPACT_AUTHORITY_BUNDLE = Object.freeze([
  'tools/web-architecture-rules.json',
  'tools/web-product-impact-map.json',
  'tools/web-product-decisions.json'
]);
const reservedPathContent = (path) => {
  if (path === 'tools/web-architecture-rules.json') {
    return '{"schemaVersion":1,"rules":[]}\n';
  }
  return path.endsWith('.json') ? '{}\n' : `// reserved fixture for ${path}\n`;
};

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
        env: {
          ...process.env,
          GITHUB_ACTIONS: '',
          GITHUB_EVENT_NAME: '',
          GITHUB_EVENT_PATH: '',
          GITHUB_REPOSITORY: '',
          GITHUB_STEP_SUMMARY: '',
          ...environment
        },
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

const promotionReportSection = (output, heading) => {
  const escapedHeading = heading.replace(/[.*+?^${}()|[\]\\]/g, '\\$&');
  return output.match(new RegExp(`## ${escapedHeading}\\n\\n([\\s\\S]*?)(?=\\n## |$)`))?.[1] || '';
};

const promotionEventSource = ({
  base,
  head,
  baseRef = 'main',
  headRef = 'dev',
  baseRepository = 'satoshikawato/gbdraw',
  headRepository = 'satoshikawato/gbdraw'
}) => JSON.stringify({
  pull_request: {
    base: {
      ref: baseRef,
      sha: base,
      repo: { full_name: baseRepository }
    },
    head: {
      ref: headRef,
      sha: head,
      repo: { full_name: headRepository }
    }
  }
});

const executePromotion = ({
  execute,
  root,
  base,
  head,
  event = {},
  environment = {}
}) => {
  const eventPath = join(root, '.git', 'promotion-event.json');
  writeFileSync(eventPath, promotionEventSource({ base, head, ...event }), 'utf8');
  return execute({
    base,
    head,
    environment: {
      GITHUB_EVENT_NAME: 'pull_request',
      GITHUB_EVENT_PATH: eventPath,
      GITHUB_REPOSITORY: 'satoshikawato/gbdraw',
      ...environment
    }
  });
};

const writeBudgetFiles = (write, count) => {
  for (let index = 0; index < count; index += 1) {
    write(
      `gbdraw/web/js/app/budget-file-${index}.js`,
      `const budgetFile${index} = ${index};\n`
    );
  }
};

const cloneBudgetPolicy = () => JSON.parse(JSON.stringify(BUDGET_POLICY));
const writeBudgetPolicy = (write, policy, space = 2) => {
  write('tools/web-change-policy.json', `${JSON.stringify(policy, null, space)}\n`);
};
const removePolicyPath = (policy, section, key, path) => {
  policy[section][key] = policy[section][key].filter((candidate) => candidate !== path);
};
const removeEditorOwnerUse = (write) => write(
  'gbdraw/web/js/app/editor.js',
  "import { runDiagramGeneration } from '../services/diagram-generation.js';\n"
    + 'export const editExistingOwner = () => 1;\n'
);
const removeEditorOwnerAndImporterUse = (write) => write(
  'gbdraw/web/js/app/editor.js',
  'export const editExistingOwner = () => 1;\n'
);
const applySingleOwnerContraction = (write) => {
  const policy = cloneBudgetPolicy();
  removePolicyPath(
    policy,
    'allowedPrivilegedOwners',
    'Diagram Worker',
    'app/editor.js'
  );
  writeBudgetPolicy(write, policy);
  removeEditorOwnerUse(write);
};
const runRevisionChangeBudgetCase = (mutate, environment = {}) => withChangeBudgetRepository(
  ({ commit, execute, git, write }) => {
    const base = git(['rev-parse', 'HEAD']).stdout.trim();
    mutate(write);
    const head = commit('runtime and policy change');
    return execute({ base, head, environment });
  }
);

const reportSection = (output, heading) => {
  const marker = `## ${heading}\n\n`;
  const start = output.indexOf(marker);
  assert.notEqual(start, -1, `missing report section: ${heading}\n${output}`);
  const contentStart = start + marker.length;
  const nextHeading = output.indexOf('\n## ', contentStart);
  return output.slice(contentStart, nextHeading === -1 ? output.length : nextHeading);
};

const workflowWarningCount = (output) => (
  output.match(/::warning title=Web policy review required::/g) || []
).length;

const assertNonWaivableWorkingTreeFailure = (mutate, expectedPatterns, context = '') => {
  BUDGET_PROFILES.forEach(({ environment, name }) => {
    const result = runChangeBudgetCase(mutate, environment);
    assert.equal(result.status, 1, `${context} ${name}\n${result.output}`);
    expectedPatterns.forEach((pattern) => assert.match(result.output, pattern));
  });
};

const assertNonWaivableRevisionFailure = (execute, base, head, expectedPatterns) => {
  for (const environment of [{}, { WEB_ARCHITECTURE_CHANGE: 'true' }]) {
    const result = execute({ base, head, environment });
    assert.equal(result.status, 1, result.output);
    expectedPatterns.forEach((pattern) => assert.match(result.output, pattern));
  }
};

const assertRevisionChangeBudgetFailure = (mutate, expectedPatterns) => {
  withChangeBudgetRepository(({ commit, execute, git, write }) => {
    const base = git(['rev-parse', 'HEAD']).stdout.trim();
    mutate(write);
    const head = commit('invalid runtime and policy change');
    assertNonWaivableRevisionFailure(execute, base, head, expectedPatterns);
  });
};

const TRUSTED_ARCHITECTURE_FIXTURE_FILES = Object.freeze({
  'package.json': '{"private":true}\n',
  'tests/web/architecture-contracts.test.mjs': (
    "test('production rendering crosses the canonical request and Worker boundary', () => {});\n"
  ),
  'tests/web/contracts/session-regenerate-intent.playwright.spec.js': (
    "test('no-draft session preserves preview, active config, canonical request, and catalog', () => {});\n"
  ),
  'tests/web/session-request.test.mjs': (
    "assert.deepEqual(roundTripCanonical.renderRequest, canonical.renderRequest,\n"
      + "  'canonical Session projection rebuilds the same render request');\n"
  ),
  'tools/check-web-change-budget.mjs': readFileSync(CHANGE_BUDGET_CHECKER, 'utf8'),
  'tools/web-architecture-detectors.mjs': readFileSync(
    WEB_ARCHITECTURE_DETECTOR_MODULE,
    'utf8'
  ),
  'tools/web-architecture-evaluation.mjs': readFileSync(
    WEB_ARCHITECTURE_EVALUATION_MODULE,
    'utf8'
  ),
  'tools/web-architecture-rules.json': WEB_ARCHITECTURE_RULES_SOURCE,
  'tools/web-product-impact-evaluation.mjs': readFileSync(
    PRODUCT_IMPACT_EVALUATION_MODULE,
    'utf8'
  ),
  'tools/web-product-impact-decision-source.mjs': readFileSync(
    PRODUCT_IMPACT_DECISION_MODULE,
    'utf8'
  ),
  'tools/web-product-impact-map.json': PRODUCT_IMPACT_MAP_SOURCE,
  'tools/web-product-decisions.json': PRODUCT_DECISIONS_SOURCE,
  'tools/web-change-policy.json': readFileSync(
    join(REPOSITORY_ROOT, 'tools/web-change-policy.json'),
    'utf8'
  ),
  'tools/web-change-source.mjs': readFileSync(
    join(REPOSITORY_ROOT, 'tools/web-change-source.mjs'),
    'utf8'
  ),
  'tools/web-promotion-context.mjs': readFileSync(
    WEB_PROMOTION_CONTEXT_MODULE,
    'utf8'
  ),
  'gbdraw/web/js/app/run-analysis.js': (
    "import { runDiagramGeneration } from '../services/diagram-generation.js';\n"
    + "import { buildCanonicalRenderRequest } from '../services/session-request.js';\n"
    + 'export const run = () => {\n'
    + '  const request = buildCanonicalRenderRequest();\n'
    + '  return runDiagramGeneration(request);\n'
    + '};\n'
  ),
  'gbdraw/web/js/services/diagram-generation.js': (
    'export const runDiagramGeneration = () => 1;\n'
  ),
  'gbdraw/web/js/services/session-request.js': (
    'export const buildCanonicalRenderRequest = () => ({});\n'
  )
});

const canonicalCallerSource = (
  "import { runDiagramGeneration } from '../services/diagram-generation.js';\n"
  + "import { buildCanonicalRenderRequest } from '../services/session-request.js';\n"
  + 'export const run = () => {\n'
  + '  const request = buildCanonicalRenderRequest();\n'
  + '  return runDiagramGeneration(request);\n'
  + '};\n'
);
const rulesSource = (mutate) => {
  const rules = JSON.parse(WEB_ARCHITECTURE_RULES_SOURCE);
  mutate(rules);
  return `${JSON.stringify(rules, null, 2)}\n`;
};
const canonicalRule = (rules) => rules.rules.find(
  ({ key }) => key === 'canonical-path.render-request'
);
const preauthorizedCanonicalRulesSource = rulesSource((rules) => {
  canonicalRule(rules).allowedEdges = [
    { from: 'app/alternate.js', to: 'services/session-request.js' },
    { from: 'app/run-analysis.js', to: 'services/session-request.js' }
  ];
});
const reportOnlyCanonicalRulesSource = rulesSource((rules) => {
  canonicalRule(rules).enforcement = 'report-only';
});
const reportOnlyProductImpactMapSource = (() => {
  const map = JSON.parse(PRODUCT_IMPACT_MAP_SOURCE);
  map.ruleCoverage.forEach((coverage) => { coverage.enforcement = 'report-only'; });
  return `${JSON.stringify(map, null, 2)}\n`;
})();
const productImpactMapSourceForRules = (architectureSource) => {
  const map = JSON.parse(PRODUCT_IMPACT_MAP_SOURCE);
  const architecture = JSON.parse(architectureSource);
  const rules = new Map(architecture.rules.map((rule) => [rule.key, rule]));
  const concern = map.concerns[0];
  const requirements = new Map(
    concern.options[0].requirements.map((requirement) => [requirement.id, requirement])
  );
  requirements.get('REQ-RENDER-REQUEST-ENTRY').anyOf = rules
    .get('canonical-path.render-request').allowedEdges.map(({ from, to }) => ({
      ruleKey: 'canonical-path.render-request',
      subjectCategory: 'canonical-entry-edge',
      subject: `${from} -> ${to}`
    }));
  requirements.get('REQ-RENDER-REQUEST-OWNER').anyOf = rules
    .get('semantic-owner.render-request').allowedDefinitionPaths.map((subject) => ({
      ruleKey: 'semantic-owner.render-request',
      subjectCategory: 'definition-path',
      subject
    }));
  return `${JSON.stringify(map, null, 2)}\n`;
};
const canonicalTransitionProductImpactMapSource = ({
  enforcement = 'report-only',
  resolution = { kind: 'decision-required' }
} = {}) => {
  const map = JSON.parse(PRODUCT_IMPACT_MAP_SOURCE);
  const concern = map.concerns[0];
  concern.scenarioRevision = 2;
  concern.options = [
    {
      id: 'after-entry-option',
      summary: 'Generation uses the preauthorized alternate entry and preserves the supported Result and Session continuation.',
      effectRefs: [
        'EFF-RENDER-REQUEST-CURRENT-STATE',
        'EFF-RENDER-REQUEST-ROUNDTRIP'
      ],
      requirements: [
        {
          id: 'REQ-AFTER-ENTRY',
          category: 'SEMANTIC',
          checkpointRefs: ['JRN-RENDER-001:generate-submit'],
          anyOf: [{
            ruleKey: 'canonical-path.render-request',
            subjectCategory: 'canonical-entry-edge',
            subject: 'app/alternate.js -> services/session-request.js'
          }]
        },
        {
          id: 'REQ-AFTER-OWNER',
          category: 'SEMANTIC',
          checkpointRefs: [
            'JRN-RENDER-001:generate-submit',
            'JRN-RENDER-001:roundtrip-regeneration'
          ],
          anyOf: [{
            ruleKey: 'semantic-owner.render-request',
            subjectCategory: 'definition-path',
            subject: 'services/session-request.js'
          }]
        }
      ]
    },
    {
      id: 'before-entry-option',
      summary: 'Generation uses the current entry and preserves the supported Result and Session continuation.',
      effectRefs: [
        'EFF-RENDER-REQUEST-CURRENT-STATE',
        'EFF-RENDER-REQUEST-ROUNDTRIP'
      ],
      requirements: [
        {
          id: 'REQ-BEFORE-ENTRY',
          category: 'SEMANTIC',
          checkpointRefs: ['JRN-RENDER-001:generate-submit'],
          anyOf: [{
            ruleKey: 'canonical-path.render-request',
            subjectCategory: 'canonical-entry-edge',
            subject: 'app/run-analysis.js -> services/session-request.js'
          }]
        },
        {
          id: 'REQ-BEFORE-OWNER',
          category: 'SEMANTIC',
          checkpointRefs: [
            'JRN-RENDER-001:generate-submit',
            'JRN-RENDER-001:roundtrip-regeneration'
          ],
          anyOf: [{
            ruleKey: 'semantic-owner.render-request',
            subjectCategory: 'definition-path',
            subject: 'services/session-request.js'
          }]
        }
      ]
    }
  ];
  concern.resolution = resolution;
  map.ruleCoverage.forEach((coverage) => { coverage.enforcement = enforcement; });
  if (enforcement === 'hard') {
    concern.contracts.forEach((contract) => {
      contract.execution = 'PR_GATE';
    });
  }
  return `${JSON.stringify(map, null, 2)}\n`;
};
const preauthorizedProductPolicy = JSON.parse(JSON.stringify(WEB_CHANGE_POLICY));
for (const target of ['services/diagram-generation.js', 'services/session-request.js']) {
  preauthorizedProductPolicy.allowedPrivilegedImporters[target].push('app/alternate.js');
  preauthorizedProductPolicy.allowedPrivilegedImporters[target].sort();
}
for (const capability of ['Diagram Worker', 'Render request']) {
  preauthorizedProductPolicy.allowedPrivilegedOwners[capability].push('app/alternate.js');
  preauthorizedProductPolicy.allowedPrivilegedOwners[capability].sort();
}
const preauthorizedProductBaseFiles = (mapSource) => ({
  'tools/web-architecture-rules.json': preauthorizedCanonicalRulesSource,
  'tools/web-change-policy.json': `${JSON.stringify(preauthorizedProductPolicy, null, 2)}\n`,
  'tools/web-product-impact-map.json': mapSource
});
const moveCanonicalEntry = (write) => {
  write('gbdraw/web/js/app/run-analysis.js', 'export const run = () => 1;\n');
  write('gbdraw/web/js/app/alternate.js', canonicalCallerSource);
};
const currentDecisionBody = (head, overrides = {}) => {
  const decision = {
    concernKey: 'product.canonical-render-request-boundary',
    scenarioRevision: 2,
    selectedOptionId: 'after-entry-option',
    acceptedImpactClass: 'AFFORDANCE_PRESERVED',
    rationale: 'The entry changes while the supported Result and Session continuation remains intact.',
    acknowledgedChangedRequirementRefs: ['REQ-AFTER-ENTRY', 'REQ-BEFORE-ENTRY'],
    evidenceRefs: [
      'tests/web/architecture-contracts.test.mjs::production rendering crosses the canonical request and Worker boundary'
    ],
    residualRisk: 'Alternate modes remain covered by the mapped staging contract.',
    ...overrides
  };
  return [
    '<!-- gbdraw-product-impact-decision:start -->',
    JSON.stringify({ schemaVersion: 1, headSha: head, decisions: [decision] }),
    '<!-- gbdraw-product-impact-decision:end -->'
  ].join('\n');
};
const productImpactPullRequestEnvironment = ({
  author = 'satoshikawato',
  baseRef = 'dev',
  body = '',
  eventName = 'pull_request_target',
  summary = false
} = {}) => ({ base, head, root }) => {
  const eventPath = join(root, '.git', 'product-impact-event.json');
  const resolvedBody = typeof body === 'function' ? body({ base, head }) : body;
  writeFileSync(eventPath, JSON.stringify({
    pull_request: {
      base: {
        ref: baseRef,
        sha: base,
        repo: { full_name: 'satoshikawato/gbdraw' }
      },
      head: {
        ref: 'feature/product-impact-fixture',
        sha: head,
        repo: { full_name: 'satoshikawato/gbdraw' }
      },
      user: { login: author },
      body: resolvedBody
    }
  }), 'utf8');
  return {
    GITHUB_EVENT_NAME: eventName,
    GITHUB_EVENT_PATH: eventPath,
    GITHUB_REPOSITORY: 'satoshikawato/gbdraw',
    ...(summary ? { GITHUB_STEP_SUMMARY: join(root, 'product-impact-summary.md') } : {})
  };
};

const withTrustedArchitectureRepository = (
  mutateHead,
  runCase,
  baseFiles = {},
  executionEnvironment = {}
) => {
  const root = mkdtempSync(join(tmpdir(), 'gbdraw-architecture-ratchet-'));
  try {
    const fixtureFiles = { ...TRUSTED_ARCHITECTURE_FIXTURE_FILES, ...baseFiles };
    if (!Object.hasOwn(baseFiles, 'tools/web-product-impact-map.json')) {
      fixtureFiles['tools/web-product-impact-map.json'] = productImpactMapSourceForRules(
        fixtureFiles['tools/web-architecture-rules.json']
      );
    }
    Object.entries(fixtureFiles)
      .forEach(([path, content]) => writeFixtureFile(root, path, content));
    const git = (args) => spawnSync('git', args, {
      cwd: root,
      encoding: 'utf8',
      stdio: ['ignore', 'pipe', 'pipe']
    });
    assert.equal(git(['init', '--quiet']).status, 0);
    assert.equal(git(['config', 'user.email', 'ratchet@example.invalid']).status, 0);
    assert.equal(git(['config', 'user.name', 'Ratchet Fixture']).status, 0);
    assert.equal(git(['add', '.']).status, 0);
    assert.equal(git(['commit', '--quiet', '-m', 'trusted base']).status, 0);
    const base = git(['rev-parse', 'HEAD']).stdout.trim();
    mutateHead((path, content) => writeFixtureFile(root, path, content));
    assert.equal(git(['add', '-A']).status, 0);
    assert.equal(git(['commit', '--quiet', '-m', 'candidate head']).status, 0);
    const head = git(['rev-parse', 'HEAD']).stdout.trim();
    assert.equal(git(['checkout', '--quiet', '--detach', base]).status, 0);
    const extraEnvironment = typeof executionEnvironment === 'function'
      ? executionEnvironment({ base, head, root })
      : executionEnvironment;
    const result = spawnSync(process.execPath, [
      'tools/check-web-change-budget.mjs',
      '--base', base,
      '--head', head
    ], {
      cwd: root,
      encoding: 'utf8',
      env: {
        ...process.env,
        GITHUB_ACTIONS: '',
        GITHUB_EVENT_NAME: '',
        GITHUB_EVENT_PATH: '',
        GITHUB_REPOSITORY: '',
        GITHUB_STEP_SUMMARY: '',
        ...extraEnvironment
      },
      stdio: ['ignore', 'pipe', 'pipe']
    });
    runCase({
      base,
      head,
      root,
      status: result.status,
      output: `${result.stdout || ''}${result.stderr || ''}`
    });
  } finally {
    rmSync(root, { recursive: true, force: true });
  }
};

test('active base rules accept one canonical entry edge', () => {
  withTrustedArchitectureRepository(
    (write) => write('fixture-note.txt', 'unchanged architecture\n'),
    ({ status, output }) => {
      assert.equal(status, 0, output);
      assert.match(
        output,
        /canonical-path\.render-request .* observation=CONFORMING .* decision=PASS/
      );
      assert.match(
        output,
        /semantic-owner\.render-request .* registeredDefinitionLocations=1/
      );
      assert.match(
        output,
        /First-party static import graph[\s\S]*Before: 3 modules, 2 edges, 0 cycles[\s\S]*After: 3 modules, 2 edges, 0 cycles/
      );
      assert.match(
        output,
        /Registered Authority Location Count \| 1 \| 0 \| 0 \| 1 \| 0 \| registered rule/
      );
    }
  );
});

test('an extra registered authority location fails and remains visible in the delta', () => {
  withTrustedArchitectureRepository(
    (write) => write(
      'gbdraw/web/js/app/extra-authority.js',
      'export const buildCanonicalRenderRequest = () => ({});\n'
    ),
    ({ status, output }) => {
      assert.equal(status, 1, output);
      assert.match(
        output,
        /semantic-owner\.render-request .* observation=DIVERGENT .* decision=FAIL/
      );
      assert.match(
        output,
        /Registered Authority Location Count \| 1 \| 1 \| 0 \| 2 \| \+1 \| registered rule/
      );
      assert.match(output, /semantic-owner\.render-request: app\/extra-authority\.js/);
    }
  );
});

test('an unpreauthorized registered authority move fails', () => {
  withTrustedArchitectureRepository(
    (write) => {
      write(
        'gbdraw/web/js/services/session-request.js',
        'export const retiredRequestHelper = () => ({});\n'
      );
      write(
        'gbdraw/web/js/app/moved-authority.js',
        'export const buildCanonicalRenderRequest = () => ({});\n'
      );
    },
    ({ status, output }) => {
      assert.equal(status, 1, output);
      assert.match(
        output,
        /semantic-owner\.render-request \| subject=app\/moved-authority\.js .* observation=DIVERGENT .* decision=FAIL/
      );
      assert.match(output, /Removed \(1\):[\s\S]*services\/session-request\.js/);
    }
  );
});

test('a first-party static self-import fails as a one-node cycle', () => {
  withTrustedArchitectureRepository(
    (write) => write(
      'gbdraw/web/js/app/self-cycle.js',
      "import './self-cycle.js';\nexport const selfCycle = true;\n"
    ),
    ({ status, output }) => {
      assert.equal(status, 1, output);
      assert.match(output, /first-party static import cycles are not allowed/);
      assert.match(
        output,
        /nodes=\[app\/self-cycle\.js\]; edges=\[app\/self-cycle\.js -> app\/self-cycle\.js\]/
      );
    }
  );
});

test('an unresolved first-party static import fails closed', () => {
  withTrustedArchitectureRepository(
    (write) => write(
      'gbdraw/web/js/app/unresolved.js',
      "import './missing-module.js';\nexport const unresolved = true;\n"
    ),
    ({ status, output }) => {
      assert.equal(status, 1, output);
      assert.match(output, /Gate: \*\*FAIL\*\*/);
      assert.match(output, /head first-party static import graph is incomplete/);
      assert.match(output, /cannot resolve first-party static import \.\/missing-module\.js/);
    }
  );
});

test('two-node and three-node first-party static cycles fail', () => {
  const cases = [
    {
      name: 'two-node',
      files: {
        'gbdraw/web/js/app/two-a.js': "import './two-b.js';\n",
        'gbdraw/web/js/app/two-b.js': "import './two-a.js';\n"
      },
      expected: /nodes=\[app\/two-a\.js, app\/two-b\.js\]; edges=\[app\/two-a\.js -> app\/two-b\.js, app\/two-b\.js -> app\/two-a\.js\]/
    },
    {
      name: 'three-node',
      files: {
        'gbdraw/web/js/app/three-a.js': "import './three-b.js';\n",
        'gbdraw/web/js/app/three-b.js': "import './three-c.js';\n",
        'gbdraw/web/js/app/three-c.js': "import './three-a.js';\n"
      },
      expected: /nodes=\[app\/three-a\.js, app\/three-b\.js, app\/three-c\.js\]; edges=\[app\/three-a\.js -> app\/three-b\.js, app\/three-b\.js -> app\/three-c\.js, app\/three-c\.js -> app\/three-a\.js\]/
    }
  ];

  cases.forEach(({ name, files, expected }) => {
    withTrustedArchitectureRepository(
      (write) => Object.entries(files).forEach(([path, source]) => write(path, source)),
      ({ status, output }) => {
        assert.equal(status, 1, `${name}\n${output}`);
        assert.match(output, expected);
      }
    );
  });
});

test('an acyclic graph resolves extensionless modules and directory indexes', () => {
  withTrustedArchitectureRepository(
    (write) => {
      write(
        'gbdraw/web/js/app/acyclic-root.js',
        "import './acyclic-leaf';\nimport './acyclic-directory';\n"
      );
      write('gbdraw/web/js/app/acyclic-leaf.js', 'export const leaf = true;\n');
      write('gbdraw/web/js/app/acyclic-directory/index.js', 'export const index = true;\n');
    },
    ({ status, output }) => {
      assert.equal(status, 0, output);
      assert.match(output, /After: 6 modules, 4 edges, 0 cycles/);
      assert.doesNotMatch(output, /first-party static import graph is incomplete/);
    }
  );
});

test('multiple cycles report stable sorted nodes and internal edges', () => {
  withTrustedArchitectureRepository(
    (write) => {
      write('gbdraw/web/js/app/zeta-a.js', "import './zeta-b.js';\n");
      write('gbdraw/web/js/app/zeta-b.js', "import './zeta-a.js';\n");
      write('gbdraw/web/js/app/alpha-a.js', "import './alpha-b.js';\n");
      write('gbdraw/web/js/app/alpha-b.js', "import './alpha-a.js';\n");
    },
    ({ status, output }) => {
      assert.equal(status, 1, output);
      const alpha = output.indexOf(
        'nodes=[app/alpha-a.js, app/alpha-b.js]; edges=[app/alpha-a.js -> app/alpha-b.js, app/alpha-b.js -> app/alpha-a.js]'
      );
      const zeta = output.indexOf(
        'nodes=[app/zeta-a.js, app/zeta-b.js]; edges=[app/zeta-a.js -> app/zeta-b.js, app/zeta-b.js -> app/zeta-a.js]'
      );
      assert.ok(alpha >= 0 && alpha < zeta, output);
      assert.match(
        output,
        /First-party static import cycles \| 0 \| 2 \| 0 \| 2 \| \+2 \| hard/
      );
    }
  );
});

test('compatibility-like names report and a private helper adds no authority', () => {
  withTrustedArchitectureRepository(
    (write) => {
      write(
        'gbdraw/web/js/app/private-helper.js',
        'const migrateLegacyFixture = () => 1;\nexport const usePrivateHelper = () => migrateLegacyFixture();\n'
      );
    },
    ({ status, output }) => {
      assert.equal(status, 0, output);
      assert.match(
        output,
        /Compatibility-like names \| 0 \| 1 \| 0 \| 1 \| \+1 \| report-only/
      );
      assert.match(output, /gbdraw\/web\/js\/app\/private-helper\.js: migrateLegacyFixture/);
      assert.match(
        output,
        /Registered Authority Location Count \| 1 \| 0 \| 0 \| 1 \| 0 \| registered rule/
      );
      assert.match(output, /Gate: \*\*PASS\*\*/);
      assert.match(output, /Review: \*\*REQUIRED\*\*/);
    }
  );
});

test('active hard rule fails when the canonical entry edge is absent', () => {
  withTrustedArchitectureRepository(
    (write) => write(
      'gbdraw/web/js/app/run-analysis.js',
      'export const run = () => 1;\n'
    ),
    ({ status, output }) => {
      assert.equal(status, 1, output);
      assert.match(
        output,
        /canonical-path\.render-request .* observation=ABSENT_REQUIRED .* decision=FAIL/
      );
    }
  );
});

test('active hard rule rejects an unapproved canonical entry edge', () => {
  withTrustedArchitectureRepository(
    (write) => {
      write('gbdraw/web/js/app/run-analysis.js', 'export const run = () => 1;\n');
      write('gbdraw/web/js/app/alternate.js', canonicalCallerSource);
    },
    ({ status, output }) => {
      assert.equal(status, 1, output);
      assert.match(
        output,
        /canonical-path\.render-request \| subject=app\/alternate\.js -> services\/session-request\.js \| observation=DIVERGENT .* decision=FAIL/
      );
    }
  );
});

test('two simultaneously active preauthorized edges fail in stable subject order', () => {
  withTrustedArchitectureRepository(
    (write) => write('gbdraw/web/js/app/alternate.js', canonicalCallerSource),
    ({ status, output }) => {
      assert.equal(status, 1, output);
      const alternate = output.indexOf(
        'canonical-path.render-request | subject=app/alternate.js -> services/session-request.js'
      );
      const current = output.indexOf(
        'canonical-path.render-request | subject=app/run-analysis.js -> services/session-request.js'
      );
      const semantic = output.indexOf('semantic-owner.render-request | subject=');
      assert.ok(alternate >= 0 && alternate < current && current < semantic, output);
      assert.match(output.slice(alternate, semantic), /observation=DIVERGENT[\s\S]+decision=FAIL/);
    },
    { 'tools/web-architecture-rules.json': preauthorizedCanonicalRulesSource }
  );
});

test('an inactive preauthorized edge does not fail', () => {
  withTrustedArchitectureRepository(
    (write) => write('fixture-note.txt', 'preauthorization remains inactive\n'),
    ({ status, output }) => {
      assert.equal(status, 0, output);
      assert.match(
        output,
        /canonical-path\.render-request .* observation=CONFORMING .* decision=PASS/
      );
    },
    { 'tools/web-architecture-rules.json': preauthorizedCanonicalRulesSource }
  );
});

test('a preauthorized one-for-one canonical path move passes and requires review', () => {
  const preauthorizedPolicy = JSON.parse(JSON.stringify(WEB_CHANGE_POLICY));
  for (const target of ['services/diagram-generation.js', 'services/session-request.js']) {
    preauthorizedPolicy.allowedPrivilegedImporters[target].push('app/alternate.js');
    preauthorizedPolicy.allowedPrivilegedImporters[target].sort();
  }
  for (const capability of ['Diagram Worker', 'Render request']) {
    preauthorizedPolicy.allowedPrivilegedOwners[capability].push('app/alternate.js');
    preauthorizedPolicy.allowedPrivilegedOwners[capability].sort();
  }
  withTrustedArchitectureRepository(
    (write) => {
      write('gbdraw/web/js/app/run-analysis.js', 'export const run = () => 1;\n');
      write('gbdraw/web/js/app/alternate.js', canonicalCallerSource);
    },
    ({ status, output }) => {
      assert.equal(status, 0, output);
      assert.match(output, /Gate: \*\*PASS\*\*/);
      assert.match(output, /Review: \*\*REQUIRED\*\*/);
      assert.match(
        reportSection(output, 'Review reasons'),
        /Architecture-bearing signals:.*registered canonical contracts/
      );
      assert.match(
        output,
        /canonical-path\.render-request \| subject=app\/alternate\.js -> services\/session-request\.js .* decision=PASS/
      );
    },
    {
      'tools/web-architecture-rules.json': preauthorizedCanonicalRulesSource,
      'tools/web-change-policy.json': `${JSON.stringify(preauthorizedPolicy, null, 2)}\n`
    }
  );
});

test('trusted Product Impact stays absent for unchanged source and non-authoritative pull requests', () => {
  withTrustedArchitectureRepository(
    (write) => write('fixture-note.txt', 'no registered Product Impact delta\n'),
    ({ status, output }) => {
      assert.equal(status, 0, output);
      assert.doesNotMatch(output, /## Product impact/);
    },
    {},
    productImpactPullRequestEnvironment()
  );

  withTrustedArchitectureRepository(
    moveCanonicalEntry,
    ({ status, output }) => {
      assert.equal(status, 0, output);
      assert.doesNotMatch(output, /## Product impact/);
      assert.doesNotMatch(output, /PRIVATE-CANDIDATE-BODY/);
    },
    preauthorizedProductBaseFiles(
      productImpactMapSourceForRules(preauthorizedCanonicalRulesSource)
    ),
    productImpactPullRequestEnvironment({
      eventName: 'pull_request',
      body: '<!-- gbdraw-product-impact-decision:start -->PRIVATE-CANDIDATE-BODY'
    })
  );

  withTrustedArchitectureRepository(
    moveCanonicalEntry,
    ({ status, output }) => {
      assert.equal(status, 0, output);
      assert.doesNotMatch(output, /## Product impact/);
      assert.doesNotMatch(output, /PRIVATE-NON-DEV-BODY/);
    },
    preauthorizedProductBaseFiles(
      productImpactMapSourceForRules(preauthorizedCanonicalRulesSource)
    ),
    productImpactPullRequestEnvironment({
      baseRef: 'feature-base',
      body: '<!-- gbdraw-product-impact-decision:start -->PRIVATE-NON-DEV-BODY'
    })
  );
});

test('trusted Product Impact reports a mapped provider substitution as one conforming concern', () => {
  withTrustedArchitectureRepository(
    moveCanonicalEntry,
    ({ status, output }) => {
      assert.equal(status, 0, output);
      assert.match(output, /## Product impact/);
      assert.match(output, /Runtime context: TRUSTED_BASE_PULL_REQUEST/);
      assert.match(output, /Impact: NO_USER_VISIBLE_DIFFERENCE/);
      assert.match(output, /Observation: CONFORMING/);
      assert.match(output, /Gate contribution: None \(report-only/);
      assert.equal(
        [...output.matchAll(/### product\.canonical-render-request-boundary/g)].length,
        1
      );
    },
    preauthorizedProductBaseFiles(
      productImpactMapSourceForRules(preauthorizedCanonicalRulesSource)
    ),
    productImpactPullRequestEnvironment()
  );
});

test('two changed architecture rules still produce one trusted Product Impact concern packet', () => {
  const rules = JSON.parse(preauthorizedCanonicalRulesSource);
  canonicalRule(rules).allowedEdges = [
    { from: 'app/alternate.js', to: 'services/future-session-request.js' },
    { from: 'app/run-analysis.js', to: 'services/session-request.js' }
  ];
  rules.rules.find(({ key }) => key === 'semantic-owner.render-request')
    .allowedDefinitionPaths = [
      'services/future-session-request.js',
      'services/session-request.js'
    ];
  const twoRuleSource = `${JSON.stringify(rules, null, 2)}\n`;
  withTrustedArchitectureRepository(
    (write) => {
      write('gbdraw/web/js/app/run-analysis.js', 'export const run = () => 1;\n');
      write(
        'gbdraw/web/js/app/alternate.js',
        canonicalCallerSource.replaceAll('session-request.js', 'future-session-request.js')
      );
      write(
        'gbdraw/web/js/services/session-request.js',
        'export const retiredRequestHelper = () => ({});\n'
      );
      write(
        'gbdraw/web/js/services/future-session-request.js',
        'export const buildCanonicalRenderRequest = () => ({});\n'
      );
    },
    ({ status, output }) => {
      assert.equal(status, 0, output);
      assert.match(output, /Rule: `canonical-path\.render-request`/);
      assert.match(output, /Rule: `semantic-owner\.render-request`/);
      assert.equal(
        [...output.matchAll(/### product\.canonical-render-request-boundary/g)].length,
        1
      );
      assert.match(output, /Observation: CONFORMING/);
    },
    {
      ...preauthorizedProductBaseFiles(productImpactMapSourceForRules(twoRuleSource)),
      'tools/web-architecture-rules.json': twoRuleSource
    },
    productImpactPullRequestEnvironment()
  );
});

test('report-only Product Impact regression requires review without changing a passing Gate', () => {
  withTrustedArchitectureRepository(
    (write) => write('gbdraw/web/js/app/run-analysis.js', 'export const run = () => 1;\n'),
    ({ status, output }) => {
      assert.equal(status, 0, output);
      assert.match(output, /Gate: \*\*PASS\*\*/);
      assert.match(output, /Review: \*\*REQUIRED\*\*/);
      assert.match(output, /Impact: RETIREMENT/);
      assert.match(output, /Observation: ORDINARY_REGRESSION/);
      assert.match(output, /Gate contribution: None \(report-only/);
    },
    {
      'tools/web-architecture-rules.json': reportOnlyCanonicalRulesSource,
      'tools/web-product-impact-map.json': reportOnlyProductImpactMapSource
    },
    productImpactPullRequestEnvironment()
  );
});

test('Decision Pack routing is outcome-oriented and accepts only an eligible exact-head decision', () => {
  const transitionMap = canonicalTransitionProductImpactMapSource();
  const baseFiles = preauthorizedProductBaseFiles(transitionMap);
  withTrustedArchitectureRepository(
    moveCanonicalEntry,
    ({ status, output }) => {
      assert.equal(status, 0, output);
      assert.match(output, /Observation: UNRESOLVED_DECISION/);
      assert.match(output, /#### Decision Pack/);
      for (const code of ['A', 'B', 'C', 'D', 'E', 'F']) {
        assert.match(output, new RegExp(`  - ${code}\\.`));
      }
      assert.match(output, /B\..*after-entry-option.*PR_LOCAL_ALLOWED/);
      assert.match(output, /Action: Configure the diagram inputs and layout -> Select Generate/);
      assert.match(output, /Effects: `EFF-RENDER-REQUEST-CURRENT-STATE`:/);
      assert.match(output, /Concern: product\.canonical-render-request-boundary/);
      assert.match(output, /Scenario revision: 2/);
      assert.match(output, /The human does not edit JSON, SHA, requirement refs, or evidence refs/);
    },
    baseFiles,
    productImpactPullRequestEnvironment()
  );

  withTrustedArchitectureRepository(
    moveCanonicalEntry,
    ({ status, output }) => {
      assert.equal(status, 0, output);
      assert.match(output, /Resolution: CURRENT_MAINTAINER_DECISION/);
      assert.match(output, /Observation: CONFORMING/);
      assert.match(
        output,
        /entry changes while the supported Result and Session continuation remains intact/
      );
      assert.match(output, /\\<script\\>not markup\\<\/script\\>/);
      assert.equal(output.includes('<script>'), false);
      assert.doesNotMatch(output, /#### Decision Pack/);
    },
    baseFiles,
    productImpactPullRequestEnvironment({
      body: ({ head }) => currentDecisionBody(head, {
        rationale: (
          'The entry changes while the supported Result and Session continuation remains intact; '
          + '<script>not markup</script>.'
        )
      })
    })
  );
});

test('stale and non-maintainer Product Impact decisions remain unresolved and routed', () => {
  const baseFiles = preauthorizedProductBaseFiles(
    canonicalTransitionProductImpactMapSource()
  );
  withTrustedArchitectureRepository(
    moveCanonicalEntry,
    ({ status, output }) => {
      assert.equal(status, 0, output);
      assert.match(output, /Current decision declaration/);
      assert.match(output, /Observation: INSUFFICIENT_EVIDENCE/);
      assert.match(output, /Observation: UNRESOLVED_DECISION/);
      assert.match(output, /review the exact current head/);
    },
    baseFiles,
    productImpactPullRequestEnvironment({
      body: () => currentDecisionBody('a'.repeat(40))
    })
  );

  withTrustedArchitectureRepository(
    moveCanonicalEntry,
    ({ status, output }) => {
      assert.equal(status, 0, output);
      assert.match(output, /not an eligible Product Decision Owner/);
      assert.match(output, /B\..*after-entry-option.*DURABLE_AUTHORITY_REQUIRED/);
    },
    baseFiles,
    productImpactPullRequestEnvironment({
      author: 'contributor',
      body: ({ head }) => currentDecisionBody(head)
    })
  );
});

test('invalid current Product Impact authority fails closed without exposing body content', () => {
  const secret = 'PRIVATE-SEQUENCE-DO-NOT-ECHO';
  const baseFiles = preauthorizedProductBaseFiles(
    canonicalTransitionProductImpactMapSource()
  );
  withTrustedArchitectureRepository(
    moveCanonicalEntry,
    ({ status, output }) => {
      assert.equal(status, 1, output);
      assert.match(output, /invalid-current-impact-class/);
      assert.doesNotMatch(output, new RegExp(secret));
    },
    baseFiles,
    productImpactPullRequestEnvironment({
      body: ({ head }) => currentDecisionBody(head, {
        acceptedImpactClass: 'PRODUCT_CHANGE',
        rationale: secret
      })
    })
  );
});

test('conflicting static and current Product Impact authority is report-only but explicit', () => {
  const map = canonicalTransitionProductImpactMapSource({
    resolution: {
      kind: 'authority-covered',
      selectedOptionId: 'before-entry-option',
      authorityRefs: ['gbdraw/web/CLAUDE.md#runtime-data-flow']
    }
  });
  withTrustedArchitectureRepository(
    moveCanonicalEntry,
    ({ status, output }) => {
      assert.equal(status, 0, output);
      assert.match(output, /Resolution: CONFLICT/);
      assert.match(output, /Observation: AUTHORITY_CONFLICT/);
      assert.match(output, /`STATIC_AUTHORITY` -> `before-entry-option`/);
      assert.match(output, /`CURRENT_MAINTAINER_DECISION` -> `after-entry-option`/);
      assert.match(output, /authority-only change/);
    },
    preauthorizedProductBaseFiles(map),
    productImpactPullRequestEnvironment({
      body: ({ head }) => currentDecisionBody(head)
    })
  );
});

test('candidate Product Impact authority is validation-only and cannot authorize candidate runtime', () => {
  withTrustedArchitectureRepository(
    (write) => {
      moveCanonicalEntry(write);
      write('tools/web-architecture-rules.json', preauthorizedCanonicalRulesSource);
      write(
        'tools/web-product-impact-map.json',
        productImpactMapSourceForRules(preauthorizedCanonicalRulesSource)
      );
    },
    ({ status, output }) => {
      assert.equal(status, 1, output);
      assert.match(output, /Candidate authority validation: VALID \(inert data only\)/);
      assert.match(output, /candidate data does not alter this head runtime admission/);
      assert.match(output, /Observation: INSUFFICIENT_EVIDENCE/);
      assert.doesNotMatch(output, /Observation: CONFORMING/);
    },
    {},
    productImpactPullRequestEnvironment()
  );
});

test('malformed candidate Product Impact authority fails closed', () => {
  withTrustedArchitectureRepository(
    (write) => write('tools/web-product-decisions.json', '{'),
    ({ status, output }) => {
      assert.equal(status, 1, output);
      assert.match(output, /candidate Product Impact authority: Cannot parse/);
      assert.match(output, /Gate: \*\*FAIL\*\*/);
    },
    {},
    productImpactPullRequestEnvironment()
  );
});

test('mapped contract changes are surfaced while unrelated tests add no Product Impact burden', () => {
  const rules = JSON.parse(WEB_ARCHITECTURE_RULES_SOURCE);
  rules.rules.forEach((rule) => { rule.enforcement = 'report-only'; });
  const reportOnlyRules = `${JSON.stringify(rules, null, 2)}\n`;
  withTrustedArchitectureRepository(
    (write) => {
      write(
        'gbdraw/web/js/services/session-request.js',
        'export const retiredRequestHelper = () => ({});\n'
      );
      write(
        'tests/web/contracts/session-regenerate-intent.playwright.spec.js',
        "test('no-draft session preserves preview, active config, canonical request, and catalog', () => { /* changed */ });\n"
      );
    },
    ({ status, output }) => {
      assert.equal(status, 0, output);
      assert.match(output, /integrity=CANDIDATE_MODIFIED/);
      assert.match(output, /mapped contract changed in the candidate/);
      assert.match(output, /Gate contribution: None \(report-only/);
    },
    {
      'tools/web-architecture-rules.json': reportOnlyRules,
      'tools/web-product-impact-map.json': reportOnlyProductImpactMapSource
    },
    productImpactPullRequestEnvironment()
  );

  withTrustedArchitectureRepository(
    (write) => write('tests/web/unrelated.test.mjs', 'test("unrelated", () => {});\n'),
    ({ status, output }) => {
      assert.equal(status, 0, output);
      assert.doesNotMatch(output, /## Product impact/);
    },
    {},
    productImpactPullRequestEnvironment()
  );
});

test('Product Impact stdout and GitHub summary use the same deterministic report', () => {
  withTrustedArchitectureRepository(
    (write) => write('gbdraw/web/js/app/run-analysis.js', 'export const run = () => 1;\n'),
    ({ root, status, output }) => {
      assert.equal(status, 0, output);
      const summary = readFileSync(join(root, 'product-impact-summary.md'), 'utf8');
      assert.equal(output, summary);
      assert.match(summary, /Observation: ORDINARY_REGRESSION/);
    },
    {
      'tools/web-architecture-rules.json': reportOnlyCanonicalRulesSource,
      'tools/web-product-impact-map.json': reportOnlyProductImpactMapSource
    },
    productImpactPullRequestEnvironment({ summary: true })
  );
});

test('hard Product Impact fixtures block every blocking observation', () => {
  const reportOnlyRules = rulesSource((rules) => {
    rules.rules.forEach((rule) => { rule.enforcement = 'report-only'; });
  });
  const hardCurrentMap = (() => {
    const map = JSON.parse(PRODUCT_IMPACT_MAP_SOURCE);
    map.ruleCoverage.forEach((coverage) => { coverage.enforcement = 'hard'; });
    map.concerns[0].contracts.forEach((contract) => { contract.execution = 'PR_GATE'; });
    return `${JSON.stringify(map, null, 2)}\n`;
  })();
  withTrustedArchitectureRepository(
    (write) => write('gbdraw/web/js/app/run-analysis.js', 'export const run = () => 1;\n'),
    ({ status, output }) => {
      assert.equal(status, 1, output);
      assert.match(output, /Product Impact hard enforcement:.*ORDINARY_REGRESSION/);
    },
    {
      'tools/web-architecture-rules.json': reportOnlyRules,
      'tools/web-product-impact-map.json': hardCurrentMap
    },
    productImpactPullRequestEnvironment()
  );

  const hardTransition = canonicalTransitionProductImpactMapSource({ enforcement: 'hard' });
  withTrustedArchitectureRepository(
    moveCanonicalEntry,
    ({ status, output }) => {
      assert.equal(status, 1, output);
      assert.match(output, /Product Impact hard enforcement:.*UNRESOLVED_DECISION/);
    },
    preauthorizedProductBaseFiles(hardTransition),
    productImpactPullRequestEnvironment()
  );

  const hardConflict = canonicalTransitionProductImpactMapSource({
    enforcement: 'hard',
    resolution: {
      kind: 'authority-covered',
      selectedOptionId: 'before-entry-option',
      authorityRefs: ['gbdraw/web/CLAUDE.md#runtime-data-flow']
    }
  });
  withTrustedArchitectureRepository(
    moveCanonicalEntry,
    ({ status, output }) => {
      assert.equal(status, 1, output);
      assert.match(output, /Product Impact hard enforcement:.*AUTHORITY_CONFLICT/);
    },
    preauthorizedProductBaseFiles(hardConflict),
    productImpactPullRequestEnvironment({ body: ({ head }) => currentDecisionBody(head) })
  );

  const hardEvidenceMap = JSON.parse(hardTransition);
  hardEvidenceMap.concerns[0].contracts[0].ref = (
    'tests/web/contracts/product-impact-entry.test.mjs::canonical alternate entry'
  );
  const hardEvidenceSource = `${JSON.stringify(hardEvidenceMap, null, 2)}\n`;
  withTrustedArchitectureRepository(
    (write) => {
      moveCanonicalEntry(write);
      write(
        'tests/web/contracts/product-impact-entry.test.mjs',
        "test('canonical alternate entry', () => { /* candidate modified */ });\n"
      );
    },
    ({ status, output }) => {
      assert.equal(status, 1, output);
      assert.match(output, /Product Impact hard enforcement:.*INSUFFICIENT_EVIDENCE/);
      assert.match(output, /integrity=CANDIDATE_MODIFIED/);
    },
    {
      ...preauthorizedProductBaseFiles(hardEvidenceSource),
      'tests/web/contracts/product-impact-entry.test.mjs': (
        "test('canonical alternate entry', () => {});\n"
      )
    },
    productImpactPullRequestEnvironment()
  );
});

test('report-only rules report absence without consulting an accepted store or failing', () => {
  withTrustedArchitectureRepository(
    (write) => write(
      'gbdraw/web/js/app/run-analysis.js',
      'export const run = () => 1;\n'
    ),
    ({ status, output }) => {
      assert.equal(status, 0, output);
      assert.match(
        output,
        /canonical-path\.render-request .* observation=ABSENT_REQUIRED .* mode=REPORT_ONLY .* decision=REPORT/
      );
      assert.match(output, /Gate: \*\*PASS\*\*/);
      assert.match(output, /Review: \*\*REQUIRED\*\*/);
    },
    {
      'tools/web-architecture-rules.json': reportOnlyCanonicalRulesSource,
      'tools/web-product-impact-map.json': reportOnlyProductImpactMapSource,
      'tools/web-architecture-violations.json': 'not valid JSON and must not be loaded\n'
    }
  );
});

test('hard rules do not consult an accepted-violation store', () => {
  withTrustedArchitectureRepository(
    (write) => write('fixture-note.txt', 'hard rule remains conforming\n'),
    ({ status, output }) => {
      assert.equal(status, 0, output);
      assert.doesNotMatch(output, /Cannot parse tools\/web-architecture-violations\.json/);
      assert.match(output, /Gate: \*\*PASS\*\*/);
      assert.match(output, /Review: \*\*CLEAR\*\*/);
    },
    { 'tools/web-architecture-violations.json': 'not valid JSON and must not be loaded\n' }
  );
});

test('trusted checker rejects malformed proposed rule JSON as inert data', () => {
  withTrustedArchitectureRepository(
    (write) => write('tools/web-architecture-rules.json', '{"schemaVersion":1,"rules":['),
    ({ status, output }) => {
      assert.equal(status, 1, output);
      assert.match(output, /Cannot parse tools\/web-architecture-rules\.json/);
      assert.match(output, /Gate: \*\*FAIL\*\*/);
      assert.match(output, /Review: \*\*REQUIRED\*\*/);
    }
  );
});

test('a proposed detector that throws is never loaded or executed', () => {
  withTrustedArchitectureRepository(
    (write) => write(
      'tools/web-architecture-detectors.mjs',
      'throw new Error("MALICIOUS HEAD DETECTOR EXECUTED");\n'
    ),
    ({ status, output }) => {
      assert.equal(status, 0, output);
      assert.doesNotMatch(output, /MALICIOUS HEAD DETECTOR EXECUTED/);
      assert.match(output, /canonical-path\.render-request .* decision=PASS/);
    }
  );
});

test('trusted fixture execution imports detector and evaluator modules from its base revision', () => {
  const localDetector = `${WEB_ARCHITECTURE_DETECTOR_SOURCE}\nprocess.stderr.write('FIXTURE_LOCAL_DETECTOR\\n');\n`;
  const localEvaluator = `${WEB_ARCHITECTURE_EVALUATION_SOURCE}\nprocess.stderr.write('FIXTURE_LOCAL_EVALUATOR\\n');\n`;
  withTrustedArchitectureRepository(
    (write) => write('fixture-note.txt', 'exercise fixture-local modules\n'),
    ({ status, output }) => {
      assert.equal(status, 0, output);
      assert.match(output, /FIXTURE_LOCAL_DETECTOR/);
      assert.match(output, /FIXTURE_LOCAL_EVALUATOR/);
    },
    {
      'tools/web-architecture-detectors.mjs': localDetector,
      'tools/web-architecture-evaluation.mjs': localEvaluator
    }
  );
});

test('proposed detector ID and implementation cannot alter trusted source observation', () => {
  const proposedRules = rulesSource((rules) => {
    canonicalRule(rules).detector = 'canonical-path.render-request.v999';
  });
  withTrustedArchitectureRepository(
    (write) => {
      write('gbdraw/web/js/app/run-analysis.js', 'export const run = () => 1;\n');
      write('gbdraw/web/js/app/alternate.js', canonicalCallerSource);
      write('tools/web-architecture-rules.json', proposedRules);
      write(
        'tools/web-architecture-detectors.mjs',
        'throw new Error("MALICIOUS HEAD DETECTOR EXECUTED");\n'
      );
    },
    ({ status, output }) => {
      assert.equal(status, 1, output);
      assert.doesNotMatch(output, /MALICIOUS HEAD DETECTOR EXECUTED/);
      assert.match(output, /unknown-detector/);
      assert.match(
        output,
        /canonical-path\.render-request .* observation=DIVERGENT .* decision=FAIL/
      );
    }
  );
});

test('cross-validated inert architecture and Product Impact preauthorization passes and requires review', () => {
  withTrustedArchitectureRepository(
    (write) => {
      write('tools/web-architecture-rules.json', preauthorizedCanonicalRulesSource);
      write(
        'tools/web-product-impact-map.json',
        productImpactMapSourceForRules(preauthorizedCanonicalRulesSource)
      );
    },
    ({ status, output }) => {
      assert.equal(status, 0, output);
      assert.match(output, /Gate: \*\*PASS\*\*/);
      assert.match(output, /Review: \*\*REQUIRED\*\*/);
      assert.match(
        reportSection(output, 'Review reasons'),
        /registered architecture-rule authority changed/
      );
      assert.match(output, /Classification: EXPANSION/);
    }
  );
});

test('mutually conforming proposed rules and source cannot replace base authority', () => {
  const proposedRules = rulesSource((rules) => {
    canonicalRule(rules).allowedEdges = [{
      from: 'app/alternate.js',
      to: 'services/session-request.js'
    }];
  });
  withTrustedArchitectureRepository(
    (write) => {
      write('gbdraw/web/js/app/run-analysis.js', 'export const run = () => 1;\n');
      write('gbdraw/web/js/app/alternate.js', canonicalCallerSource);
      write('tools/web-architecture-rules.json', proposedRules);
    },
    ({ status, output }) => {
      assert.equal(status, 1, output);
      assert.match(
        output,
        /canonical-path\.render-request \| subject=app\/alternate\.js -> services\/session-request\.js \| observation=DIVERGENT .* decision=FAIL/
      );
      assert.match(output, /architecture rule authority changes must be isolated/);
      assert.match(output, /Gate: \*\*FAIL\*\*/);
      assert.match(output, /Review: \*\*REQUIRED\*\*/);
    }
  );
});

test('fixed profiles classify exact and excess production-file review thresholds', () => {
  BUDGET_PROFILES.forEach((profile) => {
    const exact = runChangeBudgetCase(
      (write) => writeBudgetFiles(write, profile.productionFiles),
      profile.environment
    );
    assert.equal(exact.status, 0, `${profile.name} exact\n${exact.output}`);
    assert.match(exact.output, new RegExp(`Selected profile: ${profile.name}`));
    assert.match(
      exact.output,
      new RegExp(`Size-review threshold for production files: ${profile.productionFiles}`)
    );
    assert.match(exact.output, /Gate: \*\*PASS\*\*/);
    assert.match(exact.output, /Review: \*\*REQUIRED\*\*/);
    assert.doesNotMatch(reportSection(exact.output, 'Review reasons'), /Size and scope:/);

    const excess = runChangeBudgetCase(
      (write) => writeBudgetFiles(write, profile.productionFiles + 1),
      profile.environment
    );
    assert.equal(excess.status, 0, `${profile.name} excess\n${excess.output}`);
    assert.match(excess.output, /Gate: \*\*PASS\*\*/);
    assert.match(excess.output, /Review: \*\*REQUIRED\*\*/);
    assert.match(
      reportSection(excess.output, 'Review reasons'),
      new RegExp(
        `production files changed exceed ${profile.name} size-review threshold `
        + `\\(${profile.productionFiles + 1} > ${profile.productionFiles}\\)`
      )
    );
    assert.doesNotMatch(
      reportSection(excess.output, 'Blocking violations'),
      /production files changed exceed/
    );
  });
});

test('fixed profiles classify exact and excess gross-churn review thresholds', () => {
  BUDGET_PROFILES.forEach((profile) => {
    const replacementLines = profile.grossChurn / 2;
    const exact = runChangeBudgetCase((write) => {
      write('gbdraw/web/js/app/budget-lines.js', budgetLineSource(replacementLines));
    }, profile.environment);
    assert.equal(exact.status, 0, `${profile.name} exact\n${exact.output}`);
    assert.match(exact.output, new RegExp(`Gross churn: ${profile.grossChurn}`));
    assert.match(exact.output, /Gate: \*\*PASS\*\*/);
    assert.match(exact.output, /Review: \*\*CLEAR\*\*/);
    assert.equal(reportSection(exact.output, 'Review reasons').trim(), '- None');

    const excess = runChangeBudgetCase((write) => {
      write('gbdraw/web/js/app/budget-lines.js', budgetLineSource(replacementLines, 1));
    }, profile.environment);
    assert.equal(excess.status, 0, `${profile.name} excess\n${excess.output}`);
    assert.match(excess.output, /Gate: \*\*PASS\*\*/);
    assert.match(excess.output, /Review: \*\*REQUIRED\*\*/);
    assert.match(
      reportSection(excess.output, 'Review reasons'),
      new RegExp(
        `production gross churn exceeds ${profile.name} size-review threshold `
        + `\\(${profile.grossChurn + 1} > ${profile.grossChurn}\\)`
      )
    );
    assert.doesNotMatch(
      reportSection(excess.output, 'Blocking violations'),
      /production gross churn exceeds/
    );
  });
});

test('fixed profiles classify exact and excess net-addition review thresholds', () => {
  BUDGET_PROFILES.forEach((profile) => {
    const exact = runChangeBudgetCase((write) => {
      write('gbdraw/web/js/app/budget-lines.js', budgetLineSource(0, profile.netAdditions));
    }, profile.environment);
    assert.equal(exact.status, 0, `${profile.name} exact\n${exact.output}`);
    assert.match(exact.output, new RegExp(`Net additions: ${profile.netAdditions}`));
    assert.match(exact.output, /Gate: \*\*PASS\*\*/);
    assert.match(exact.output, /Review: \*\*CLEAR\*\*/);
    assert.equal(reportSection(exact.output, 'Review reasons').trim(), '- None');

    const excess = runChangeBudgetCase((write) => {
      write('gbdraw/web/js/app/budget-lines.js', budgetLineSource(0, profile.netAdditions + 1));
    }, profile.environment);
    assert.equal(excess.status, 0, `${profile.name} excess\n${excess.output}`);
    assert.match(excess.output, /Gate: \*\*PASS\*\*/);
    assert.match(excess.output, /Review: \*\*REQUIRED\*\*/);
    assert.match(
      reportSection(excess.output, 'Review reasons'),
      new RegExp(
        `production net additions exceed ${profile.name} size-review threshold `
        + `\\(${profile.netAdditions + 1} > ${profile.netAdditions}\\)`
      )
    );
    assert.doesNotMatch(
      reportSection(excess.output, 'Blocking violations'),
      /production net additions exceed/
    );
  });
});

test('the architecture label selects a larger advisory profile without waiving blockers', () => {
  const mutate = (write) => writeBudgetFiles(write, 9);
  const ordinary = runChangeBudgetCase(mutate);
  assert.equal(ordinary.status, 0, ordinary.output);
  assert.match(ordinary.output, /Gate: \*\*PASS\*\*/);
  assert.match(ordinary.output, /Review: \*\*REQUIRED\*\*/);
  assert.match(
    ordinary.output,
    /production files changed exceed ordinary size-review threshold \(9 > 8\)/
  );
  assert.doesNotMatch(ordinary.output, /waived/i);

  const architecture = runChangeBudgetCase(mutate, { WEB_ARCHITECTURE_CHANGE: 'true' });
  assert.equal(architecture.status, 0, architecture.output);
  assert.match(architecture.output, /Selected profile: architecture/);
  assert.match(architecture.output, /Gate: \*\*PASS\*\*/);
  assert.match(architecture.output, /Review: \*\*REQUIRED\*\*/);
  assert.doesNotMatch(
    reportSection(architecture.output, 'Review reasons'),
    /Size and scope:/
  );
  assert.doesNotMatch(architecture.output, /waived/i);
});

test('net-zero rewrites above gross thresholds require review in both profiles', () => {
  BUDGET_PROFILES.forEach((profile) => {
    const result = runChangeBudgetCase((write) => {
      write(
        'gbdraw/web/js/app/budget-lines.js',
        budgetLineSource((profile.grossChurn / 2) + 1)
      );
    }, profile.environment);
    assert.equal(result.status, 0, `${profile.name}\n${result.output}`);
    assert.match(result.output, /Gate: \*\*PASS\*\*/);
    assert.match(result.output, /Review: \*\*REQUIRED\*\*/);
    assert.match(result.output, /Net additions: 0/);
    assert.match(result.output, /production gross churn exceeds .* size-review threshold/);
  });
});

test('new module, export, create, reactive, and watcher signals are report-only', () => {
  BUDGET_PROFILES.forEach(({ environment, name }) => {
    const result = runChangeBudgetCase((write) => {
      write(
        'gbdraw/web/js/app/report-signals.js',
        'export const createBudgetOwner = () => 1;\n'
          + 'export const budgetState = ref(0);\n'
          + 'watch(budgetState, () => {});\n'
      );
    }, environment);
    assert.equal(result.status, 0, `${name}\n${result.output}`);
    assert.match(result.output, /Report-only new production modules/);
    assert.match(result.output, /gbdraw\/web\/js\/app\/report-signals\.js/);
    assert.match(result.output, /createBudgetOwner/);
    assert.match(result.output, /budgetState \(ref\)/);
    assert.match(result.output, /\+1 watcher call/);
    assert.match(result.output, /Gate: \*\*PASS\*\*/);
    assert.match(result.output, /Review: \*\*REQUIRED\*\*/);
  });
});

test('registered performance, reference-output, and session paths require review', () => {
  const cases = [
    [
      'performance',
      (write) => write('tests/performance/gallery-publication-baseline.json', '{}\n'),
      /registered performance evidence paths changed/
    ],
    [
      'reference output',
      (write) => write('tests/reference_outputs/fixture.svg', '<svg/>\n'),
      /reference output paths changed/
    ],
    [
      'session contract',
      (write) => write(
        'gbdraw/web/js/services/session-file.js',
        'export const readSession = () => ({ version: 2 });\n'
      ),
      /registered session or compatibility paths changed/
    ]
  ];
  cases.forEach(([name, mutate, reason]) => {
    const result = runChangeBudgetCase(mutate);
    assert.equal(result.status, 0, `${name}\n${result.output}`);
    assert.match(result.output, /Gate: \*\*PASS\*\*/);
    assert.match(result.output, /Review: \*\*REQUIRED\*\*/);
    assert.match(reportSection(result.output, 'Review reasons'), reason);
  });
});

test('review reasons are deduplicated and deterministically sorted', () => {
  const result = runChangeBudgetCase((write) => {
    write('tests/performance/gallery-publication-baseline.json', '{}\n');
    write('tests/reference_outputs/fixture.svg', '<svg/>\n');
    write('gbdraw/web/js/app/review-signal.js', 'export const createCacheManager = () => 1;\n');
  });
  assert.equal(result.status, 0, result.output);
  const reasons = reportSection(result.output, 'Review reasons')
    .trim()
    .split('\n')
    .map((line) => line.replace(/^- /, ''));
  assert.deepEqual(reasons, [...new Set(reasons)].sort((left, right) => left.localeCompare(right)));
});

test('report-only inventory contractions display removals without failing', () => {
  withTrustedArchitectureRepository(
    (write) => write(
      'gbdraw/web/js/app/report-base.js',
      'export const plainReportValue = 1;\n'
    ),
    ({ status, output }) => {
      assert.equal(status, 0, output);
      assert.match(
        output,
        /create\* declarations \| 1 \| 0 \| 1 \| 0 \| -1 \| report-only/
      );
      assert.match(output, /Watcher calls \| 1 \| 0 \| 1 \| 0 \| -1 \| report-only/);
      assert.match(
        output,
        /Compatibility-like names \| 1 \| 0 \| 1 \| 0 \| -1 \| report-only/
      );
      assert.match(output, /Removed \(1\):[\s\S]*createCacheManager/);
    },
    {
      'gbdraw/web/js/app/report-base.js': (
        'export const createCacheManager = () => 1;\n'
        + 'const reactiveReport = ref(0);\n'
        + 'watch(reactiveReport, () => migrateLegacyReport());\n'
      )
    }
  );
});

test('the GitHub step summary leads with independent Gate and Review results', () => {
  BUDGET_PROFILES.forEach((profile) => {
    withChangeBudgetRepository(({ execute, root, write }) => {
      const summaryPath = join(root, 'step-summary.md');
      writeBudgetFiles(write, profile.productionFiles + 1);
      const result = execute({
        environment: {
          ...profile.environment,
          GITHUB_STEP_SUMMARY: summaryPath
        }
      });
      assert.equal(result.status, 0, `${profile.name}\n${result.output}`);
      const summary = readFileSync(summaryPath, 'utf8');
      assert.match(summary, /Gate: \*\*PASS\*\*/);
      assert.match(summary, /Review: \*\*REQUIRED\*\*/);
      assert.doesNotMatch(summary, /- Result:|- Size review:/);
      assert.match(
        summary,
        new RegExp(`Size-review threshold for production files: ${profile.productionFiles}`)
      );
      assert.match(summary, /## Key architecture differential summary/);
      assert.match(
        summary,
        /\| Inventory \| Before \| Added \| Removed \| After \| Delta \| Classification \|/
      );
      assert.match(summary, /## First-party static import graph/);
      assert.match(summary, /## Production files touched/);
      assert.match(summary, /## Production additions\/deletions/);
      assert.match(
        reportSection(summary, 'Review reasons'),
        new RegExp(
          `production files changed exceed ${profile.name} size-review threshold `
          + `\\(${profile.productionFiles + 1} > ${profile.productionFiles}\\)`
        )
      );
      assert.equal(reportSection(summary, 'Blocking violations').trim(), '- None');
    });
  });
});

test('review reasons and blocking violations produce the four independent decision states', () => {
  BUDGET_PROFILES.forEach((profile) => {
    const clear = runChangeBudgetCase((write) => {
      write(
        'gbdraw/web/js/services/diagram-generation.js',
        'export const runDiagramGeneration = () => 2;\n'
      );
    }, profile.environment);
    assert.equal(clear.status, 0, `${profile.name} clear\n${clear.output}`);
    assert.match(clear.output, /Gate: \*\*PASS\*\*/);
    assert.match(clear.output, /Review: \*\*CLEAR\*\*/);

    const reviewOnly = runChangeBudgetCase((write) => {
      write(
        'gbdraw/web/js/app/budget-lines.js',
        budgetLineSource(profile.grossChurn / 2, 1)
      );
    }, profile.environment);
    assert.equal(reviewOnly.status, 0, `${profile.name} review only\n${reviewOnly.output}`);
    assert.match(reviewOnly.output, /Gate: \*\*PASS\*\*/);
    assert.match(reviewOnly.output, /Review: \*\*REQUIRED\*\*/);

    const mixed = runChangeBudgetCase((write) => {
      writeBudgetFiles(write, profile.productionFiles + 1);
      write(
        'package.json',
        '{"private":true,"dependencies":{"left-pad":"1.3.0"}}\n'
      );
    }, profile.environment);
    assert.equal(mixed.status, 1, `${profile.name} mixed\n${mixed.output}`);
    assert.match(mixed.output, /Gate: \*\*FAIL\*\*/);
    assert.match(mixed.output, /Review: \*\*REQUIRED\*\*/);
    assert.match(
      reportSection(mixed.output, 'Review reasons'),
      /production files changed exceed/
    );
    assert.match(
      reportSection(mixed.output, 'Blocking violations'),
      /new production dependencies are not allowed/
    );
    assert.doesNotMatch(
      reportSection(mixed.output, 'Blocking violations'),
      /production files changed exceed/
    );

    const blockerOnly = runChangeBudgetCase((write) => {
      write(
        'package.json',
        '{"private":true,"dependencies":{"left-pad":"1.3.0"}}\n'
      );
    }, profile.environment);
    assert.equal(blockerOnly.status, 1, `${profile.name} blocker only\n${blockerOnly.output}`);
    assert.match(blockerOnly.output, /Gate: \*\*FAIL\*\*/);
    assert.match(blockerOnly.output, /Review: \*\*CLEAR\*\*/);
    assert.equal(reportSection(blockerOnly.output, 'Review reasons').trim(), '- None');
    assert.match(
      reportSection(blockerOnly.output, 'Blocking violations'),
      /new production dependencies are not allowed/
    );
  });
});

test('GitHub Actions emits one bounded warning only when Review is required', () => {
  BUDGET_PROFILES.forEach((profile) => {
    const required = runChangeBudgetCase((write) => {
      write(
        'gbdraw/web/js/app/budget-lines.js',
        budgetLineSource(profile.grossChurn / 2, profile.netAdditions + 1)
      );
      writeBudgetFiles(write, profile.productionFiles);
    }, { ...profile.environment, GITHUB_ACTIONS: 'true' });
    assert.equal(required.status, 0, `${profile.name} required\n${required.output}`);
    assert.match(required.output, /Review: \*\*REQUIRED\*\*/);
    assert.equal(workflowWarningCount(required.output), 1, required.output);
    const warning = required.output.split('\n').find((line) => line.startsWith('::warning'));
    assert.match(
      warning,
      new RegExp(
        '^::warning title=Web policy review required::Web policy review required: '
        + 'reasons=4; categories=Architecture-bearing signals=1, Size and scope=3; '
        + 'see the step summary\\.$'
      )
    );

    const clear = runChangeBudgetCase(
      (write) => write(
        'gbdraw/web/js/services/diagram-generation.js',
        'export const runDiagramGeneration = () => 2;\n'
      ),
      { ...profile.environment, GITHUB_ACTIONS: 'true' }
    );
    assert.equal(clear.status, 0, `${profile.name} clear\n${clear.output}`);
    assert.match(clear.output, /Review: \*\*CLEAR\*\*/);
    assert.equal(workflowWarningCount(clear.output), 0, clear.output);
  });
});

test('dependency, vendor, binary, and runtime/guard rules are non-waivable', () => {
  const cases = [
    {
      name: 'manifest production dependency',
      mutate: (write) => write(
        'package.json',
        '{"private":true,"dependencies":{"left-pad":"1.3.0"}}\n'
      ),
      expected: /new production dependencies are not allowed/
    },
    {
      name: 'bare production import',
      mutate: (write) => write(
        'gbdraw/web/js/app/secondary.js',
        "import 'left-pad';\nexport const editSecondaryOwner = () => 1;\n"
      ),
      expected: /bare production import: gbdraw\/web\/js\/app\/secondary\.js: left-pad/
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
    },
    {
      name: 'production and detector co-change',
      mutate: (write) => {
        write('gbdraw/web/js/app/editor.js', 'export const editExistingOwner = () => 2;\n');
        write(
          'tools/web-architecture-detectors.mjs',
          `${BUDGET_FIXTURE['tools/web-architecture-detectors.mjs']}\n// changed\n`
        );
      },
      expected: /production runtime files and Web guard\/CI files changed together/
    }
  ];

  cases.forEach(({ name, mutate, expected }) => {
    assertNonWaivableWorkingTreeFailure(mutate, [expected], name);
  });
});

test('malformed authority and unsupported accepted-violation authority fail canonically', () => {
  const malformedPolicy = runChangeBudgetCase((write) => {
    write('tools/web-change-policy.json', '{\n');
  });
  assert.equal(malformedPolicy.status, 1, malformedPolicy.output);
  assert.match(malformedPolicy.output, /Gate: \*\*FAIL\*\*/);
  assert.match(malformedPolicy.output, /Review: \*\*REQUIRED\*\*/);
  assert.match(malformedPolicy.output, /Web policy evaluation failed closed: Cannot parse/);

  const acceptedViolation = runChangeBudgetCase((write) => {
    write('tools/web-architecture-violations.json', '{"schemaVersion":1,"violations":[]}\n');
  });
  assert.equal(acceptedViolation.status, 1, acceptedViolation.output);
  assert.match(acceptedViolation.output, /Gate: \*\*FAIL\*\*/);
  assert.match(acceptedViolation.output, /Review: \*\*REQUIRED\*\*/);
  assert.match(acceptedViolation.output, /accepted-violation authority cannot be introduced/);
});

test('a trusted workflow cannot check out candidate head code', () => {
  const result = runChangeBudgetCase((write) => {
    write(
      '.github/workflows/web-base-policy.yml',
      'name: unsafe trusted workflow\n'
        + 'on: pull_request_target\n'
        + 'jobs:\n'
        + '  policy:\n'
        + '    steps:\n'
        + '      - uses: actions/checkout@v4\n'
        + '        with:\n'
        + '          ref: ${{ github.event.pull_request.head.sha }}\n'
    );
  });
  assert.equal(result.status, 1, result.output);
  assert.match(result.output, /Gate: \*\*FAIL\*\*/);
  assert.match(result.output, /Review: \*\*REQUIRED\*\*/);
  assert.match(result.output, /candidate trusted workflow would check out/);
});

test('reserved future guard paths need not exist before their rollout', () => {
  assert.deepEqual(
    FUTURE_GUARD_PATHS.filter((path) => Object.hasOwn(BUDGET_FIXTURE, path)),
    []
  );
  const result = runChangeBudgetCase(() => {});
  assert.equal(result.status, 0, result.output);
  assert.match(result.output, /Gate: \*\*PASS\*\*/);
  assert.match(result.output, /Review: \*\*CLEAR\*\*/);
  assert.doesNotMatch(result.output, /- Result:|- Size review:/);
});

test('every Product Impact policy, checker, authority, parser, and fixture path is guarded', () => {
  PRODUCT_IMPACT_GUARD_PATHS.forEach((guardPath) => {
    assertNonWaivableWorkingTreeFailure((write) => {
      write(
        'gbdraw/web/js/services/session-file.js',
        'export const readSession = () => ({ version: 2 });\n'
      );
      write(
        guardPath,
        BUDGET_FIXTURE[guardPath]
          ? `${BUDGET_FIXTURE[guardPath]}\n// changed\n`
          : reservedPathContent(guardPath)
      );
    }, [/production runtime files and Web guard\/CI files changed together/], guardPath);
  });
});

test('Product Impact checker implementation is separate from Product Impact authority', () => {
  PRODUCT_IMPACT_CHECKER_IMPLEMENTATION_PATHS.forEach((implementationPath) => {
    PRODUCT_IMPACT_AUTHORITY_PATHS.forEach((authorityPath) => {
      assertNonWaivableWorkingTreeFailure((write) => {
        write(
          implementationPath,
          `${BUDGET_FIXTURE[implementationPath]}\n// changed\n`
        );
        write(authorityPath, reservedPathContent(authorityPath));
      }, [/Web checker\/source parser and authority policy\/workflow files changed together/]);
    });
  });
});

test('the narrow inert authority bundle contains only architecture rules, map, and decisions', () => {
  assert.deepEqual(NARROW_PRODUCT_IMPACT_AUTHORITY_BUNDLE, [
    'tools/web-architecture-rules.json',
    'tools/web-product-impact-map.json',
    'tools/web-product-decisions.json'
  ]);

  const productAuthorityOnly = runChangeBudgetCase((write) => {
    write('tools/web-product-impact-map.json', `${PRODUCT_IMPACT_MAP_SOURCE}\n`);
    write('tools/web-product-decisions.json', `${PRODUCT_DECISIONS_SOURCE}\n`);
  });
  assert.equal(productAuthorityOnly.status, 0, productAuthorityOnly.output);
  assert.match(productAuthorityOnly.output, /Gate: \*\*PASS\*\*/);
  assert.match(productAuthorityOnly.output, /Review: \*\*REQUIRED\*\*/);

  const completeBundle = runChangeBudgetCase((write) => {
    write('tools/web-architecture-rules.json', preauthorizedCanonicalRulesSource);
    write(
      'tools/web-product-impact-map.json',
      productImpactMapSourceForRules(preauthorizedCanonicalRulesSource)
    );
    write('tools/web-product-decisions.json', `${PRODUCT_DECISIONS_SOURCE}\n`);
  });
  assert.equal(completeBundle.status, 0, completeBundle.output);
  assert.match(completeBundle.output, /Classification: EXPANSION/);

  const expandedBundle = runChangeBudgetCase((write) => {
    write('tools/web-architecture-rules.json', preauthorizedCanonicalRulesSource);
    write(
      'tools/web-product-impact-map.json',
      productImpactMapSourceForRules(preauthorizedCanonicalRulesSource)
    );
    write('tools/web-product-decisions.json', `${PRODUCT_DECISIONS_SOURCE}\n`);
    write('docs/internal/PRODUCT_IMPACT_RATCHET.md', '# changed policy\n');
  });
  assert.equal(expandedBundle.status, 1, expandedBundle.output);
  assert.match(expandedBundle.output, /architecture rule authority changes must be isolated/);
});

test('runtime cannot co-change with any protected architecture guard path', () => {
  PROTECTED_ARCHITECTURE_GUARD_PATHS.forEach((guardPath) => {
    assertNonWaivableWorkingTreeFailure((write) => {
      write(
        'gbdraw/web/js/services/session-file.js',
        'export const readSession = () => ({ version: 2 });\n'
      );
      write(guardPath, reservedPathContent(guardPath));
    }, [/production runtime files and Web guard\/CI files changed together/], guardPath);
  });
});

test('checker implementation and authority files cannot change together', () => {
  const implementationPaths = [
    'tools/check-web-change-budget.mjs',
    'tools/web-change-source.mjs',
    'tools/web-architecture-detectors.mjs',
    'tools/web-architecture-evaluation.mjs',
    'tools/web-product-impact-evaluation.mjs',
    'tools/web-product-impact-decision-source.mjs',
    'tools/web-promotion-context.mjs',
    'tools/check-promotion-readiness.mjs'
  ];
  const authorityPaths = [
    'docs/internal/ARCHITECTURE_FITNESS_FUNCTION_RATCHET.md',
    'docs/internal/PRODUCT_IMPACT_RATCHET.md',
    'docs/internal/WEB_CHANGE_POLICY.md',
    'tools/web-change-policy.json',
    'tools/web-architecture-rules.json',
    'tools/web-architecture-violations.json',
    'tools/web-product-impact-map.json',
    'tools/web-product-decisions.json',
    '.github/workflows/gallery-publication.yml',
    '.github/workflows/deploy_web.yml',
    '.github/workflows/test.yml',
    '.github/workflows/web-base-policy.yml'
  ];

  implementationPaths.forEach((implementationPath) => {
    authorityPaths.forEach((authorityPath) => {
      assertNonWaivableWorkingTreeFailure((write) => {
        write(implementationPath, `${BUDGET_FIXTURE[implementationPath]}\n// changed\n`);
        write(
          authorityPath,
          authorityPath === 'tools/web-change-policy.json'
            ? `${JSON.stringify(BUDGET_POLICY)}\n`
            : `${BUDGET_FIXTURE[authorityPath] || reservedPathContent(authorityPath)}# changed\n`
        );
      }, [/Web checker\/source parser and authority policy\/workflow files changed together/]);
    });
  });
});

test('all checker implementations are separated from reserved future authority', () => {
  const implementationPaths = [
    'tools/check-web-change-budget.mjs',
    'tools/web-change-source.mjs',
    'tools/web-architecture-detectors.mjs',
    'tools/web-architecture-evaluation.mjs',
    'tools/web-product-impact-evaluation.mjs',
    'tools/web-product-impact-decision-source.mjs',
    'tools/web-promotion-context.mjs',
    'tools/check-promotion-readiness.mjs'
  ];

  implementationPaths.forEach((implementationPath) => {
    FUTURE_AUTHORITY_PATHS.forEach((authorityPath) => {
      assertNonWaivableWorkingTreeFailure((write) => {
        write(
          implementationPath,
          BUDGET_FIXTURE[implementationPath]
            ? `${BUDGET_FIXTURE[implementationPath]}\n// changed\n`
            : reservedPathContent(implementationPath)
        );
        write(authorityPath, reservedPathContent(authorityPath));
      }, [/Web checker\/source parser and authority policy\/workflow files changed together/]);
    });
  });
});

test('promotion helper stays checker-only while evidence workflows stay authority-only', () => {
  const helperOnly = runChangeBudgetCase((write) => {
    write(
      'tools/check-promotion-readiness.mjs',
      `${BUDGET_FIXTURE['tools/check-promotion-readiness.mjs']}\n// helper changed\n`
    );
  });
  assert.equal(helperOnly.status, 0, helperOnly.output);

  for (const workflowPath of [
    '.github/workflows/gallery-publication.yml',
    '.github/workflows/deploy_web.yml'
  ]) {
    const workflowOnly = runChangeBudgetCase((write) => {
      write(workflowPath, `${BUDGET_FIXTURE[workflowPath]}# workflow changed\n`);
    });
    assert.equal(workflowOnly.status, 0, workflowOnly.output);
    assert.match(workflowOnly.output, new RegExp(workflowPath.replaceAll('.', '\\.')));
  }
});

test('evaluator fixtures stay implementation-only while active rule authority stays isolated', () => {
  const evaluatorAndFixture = runChangeBudgetCase((write) => {
    write(
      'tools/web-architecture-evaluation.mjs',
      reservedPathContent('tools/web-architecture-evaluation.mjs')
    );
    write(
      'tests/web/architecture-ratchet-fixtures.test.mjs',
      reservedPathContent('tests/web/architecture-ratchet-fixtures.test.mjs')
    );
  });
  assert.equal(evaluatorAndFixture.status, 0, evaluatorAndFixture.output);

  const fixtureAndAuthority = runChangeBudgetCase((write) => {
    write(
      'tests/web/architecture-ratchet-fixtures.test.mjs',
      reservedPathContent('tests/web/architecture-ratchet-fixtures.test.mjs')
    );
    write(
      'tools/web-architecture-rules.json',
      reservedPathContent('tools/web-architecture-rules.json')
    );
  });
  assert.equal(fixtureAndAuthority.status, 1, fixtureAndAuthority.output);
  assert.match(
    fixtureAndAuthority.output,
    /architecture rule authority changes must be isolated/
  );
});

test('the reserved PR template is guard-only', () => {
  const templateAndChecker = runChangeBudgetCase((write) => {
    write(
      '.github/pull_request_template.md',
      reservedPathContent('.github/pull_request_template.md')
    );
    write(
      'tools/web-architecture-evaluation.mjs',
      reservedPathContent('tools/web-architecture-evaluation.mjs')
    );
  });
  assert.equal(templateAndChecker.status, 0, templateAndChecker.output);

  const templateAndAuthority = runChangeBudgetCase((write) => {
    write(
      '.github/pull_request_template.md',
      reservedPathContent('.github/pull_request_template.md')
    );
    write(
      'tools/web-architecture-rules.json',
      reservedPathContent('tools/web-architecture-rules.json')
    );
  });
  assert.equal(templateAndAuthority.status, 1, templateAndAuthority.output);
  assert.match(
    templateAndAuthority.output,
    /architecture rule authority changes must be isolated/
  );
});

test('an unregistered path acquires no guard or authority classification', () => {
  const unregisteredPath = 'docs/internal/UNREGISTERED_ARCHITECTURE_NOTE.md';
  const runtimeChange = runChangeBudgetCase((write) => {
    write(
      'gbdraw/web/js/services/session-file.js',
      'export const readSession = () => ({ version: 2 });\n'
    );
    write(unregisteredPath, '# Unregistered fixture\n');
  });
  assert.equal(runtimeChange.status, 0, runtimeChange.output);

  const checkerChange = runChangeBudgetCase((write) => {
    write(
      'tools/web-architecture-evaluation.mjs',
      reservedPathContent('tools/web-architecture-evaluation.mjs')
    );
    write(unregisteredPath, '# Unregistered fixture\n');
  });
  assert.equal(checkerChange.status, 0, checkerChange.output);
});

test('the Web change-budget checker allows an owner-internal edit', () => {
  const result = runChangeBudgetCase((write) => {
    write(
      'gbdraw/web/js/services/diagram-generation.js',
      'export const runDiagramGeneration = () => 2;\n'
    );
  });
  assert.equal(result.status, 0, result.output);
  assert.match(result.output, /Gate: \*\*PASS\*\*/);
  assert.match(result.output, /Review: \*\*CLEAR\*\*/);
});

test('privileged ownership and import contractions do not require a policy edit', () => {
  const result = runChangeBudgetCase((write) => {
    write(
      'gbdraw/web/js/app/editor.js',
      'export const editExistingOwner = () => 1;\n'
    );
  });
  assert.equal(result.status, 0, result.output);
  assert.match(result.output, /Gate: \*\*PASS\*\*/);
  assert.match(result.output, /Review: \*\*REQUIRED\*\*/);
});

test('runtime removal may contract its now-inactive owner authorization in one revision', () => {
  const result = runRevisionChangeBudgetCase(applySingleOwnerContraction);
  assert.equal(result.status, 0, result.output);
  assert.match(result.output, /Gate: \*\*PASS\*\*/);
  assert.match(result.output, /Review: \*\*REQUIRED\*\*/);
  assert.match(
    result.output,
    /allowedPrivilegedOwners\.Diagram Worker: app\/editor\.js/
  );
});

test('runtime removal may contract exactly its inactive owner and importer authorizations', () => {
  const result = runRevisionChangeBudgetCase((write) => {
    const policy = cloneBudgetPolicy();
    removePolicyPath(
      policy,
      'allowedPrivilegedOwners',
      'Diagram Worker',
      'app/editor.js'
    );
    removePolicyPath(
      policy,
      'allowedPrivilegedImporters',
      'services/diagram-generation.js',
      'app/editor.js'
    );
    writeBudgetPolicy(write, policy);
    removeEditorOwnerAndImporterUse(write);
  });
  assert.equal(result.status, 0, result.output);
  assert.match(result.output, /Gate: \*\*PASS\*\*/);
  assert.match(result.output, /Review: \*\*REQUIRED\*\*/);
  assert.match(
    result.output,
    /allowedPrivilegedImporters\.services\/diagram-generation\.js: app\/editor\.js/
  );
});

test('safe contraction passes but requires review at exact size thresholds', () => {
  BUDGET_PROFILES.forEach((profile) => {
    const result = runRevisionChangeBudgetCase((write) => {
      applySingleOwnerContraction(write);
      const additionalFiles = profile.productionFiles - 2;
      const appendedLines = profile.netAdditions - additionalFiles;
      const replacementLines = (
        profile.grossChurn - profile.netAdditions - 2
      ) / 2;
      write(
        'gbdraw/web/js/app/budget-lines.js',
        budgetLineSource(replacementLines, appendedLines)
      );
      writeBudgetFiles(write, additionalFiles);
    }, profile.environment);
    assert.equal(result.status, 0, `${profile.name}\n${result.output}`);
    assert.match(result.output, /Gate: \*\*PASS\*\*/);
    assert.match(result.output, /Review: \*\*REQUIRED\*\*/);
    assert.doesNotMatch(reportSection(result.output, 'Review reasons'), /Size and scope:/);
    assert.match(
      result.output,
      new RegExp(`Size-review threshold for production files: ${profile.productionFiles}`)
    );
    assert.match(result.output, new RegExp(`Gross churn: ${profile.grossChurn}`));
    assert.match(result.output, new RegExp(`Net additions: ${profile.netAdditions}`));
  });
});

test('runtime and policy cannot self-authorize a new privileged owner', () => {
  assertRevisionChangeBudgetFailure((write) => {
    const policy = cloneBudgetPolicy();
    policy.allowedPrivilegedOwners['Diagram Worker'].push('app/secondary.js');
    writeBudgetPolicy(write, policy);
    write(
      'gbdraw/web/js/app/secondary.js',
      'export const editSecondaryOwner = () => runDiagramGeneration();\n'
    );
  }, [
    /production runtime files and Web guard\/CI files changed together/,
    /Diagram Worker: owner app\/secondary\.js/,
    /privileged capability owners or importers exceed the base allowlist/,
    /allowedPrivilegedOwners\.Diagram Worker: app\/secondary\.js/
  ]);
});

test('runtime and policy cannot self-authorize a new privileged importer', () => {
  assertRevisionChangeBudgetFailure((write) => {
    const policy = cloneBudgetPolicy();
    policy.allowedPrivilegedImporters['services/diagram-generation.js'].push(
      'app/secondary.js'
    );
    writeBudgetPolicy(write, policy);
    write(
      'gbdraw/web/js/app/secondary.js',
      "import '../services/diagram-generation.js';\n"
        + 'export const editSecondaryOwner = () => 1;\n'
    );
  }, [
    /production runtime files and Web guard\/CI files changed together/,
    /services\/diagram-generation\.js: importer app\/secondary\.js/,
    /privileged capability owners or importers exceed the base allowlist/,
    /allowedPrivilegedImporters\.services\/diagram-generation\.js: app\/secondary\.js/
  ]);
});

test('runtime contraction cannot mix removal with a new authorization', () => {
  assertRevisionChangeBudgetFailure((write) => {
    const policy = cloneBudgetPolicy();
    removePolicyPath(
      policy,
      'allowedPrivilegedOwners',
      'Diagram Worker',
      'app/editor.js'
    );
    policy.allowedPrivilegedOwners['Diagram Worker'].push('app/secondary.js');
    writeBudgetPolicy(write, policy);
    removeEditorOwnerUse(write);
  }, [
    /production runtime files and Web guard\/CI files changed together/,
    /allowedPrivilegedOwners\.Diagram Worker: app\/editor\.js/,
    /allowedPrivilegedOwners\.Diagram Worker: app\/secondary\.js/
  ]);
});

test('runtime contraction cannot remove an authorization still active at head', () => {
  assertRevisionChangeBudgetFailure((write) => {
    const policy = cloneBudgetPolicy();
    removePolicyPath(
      policy,
      'allowedPrivilegedOwners',
      'Diagram Worker',
      'app/editor.js'
    );
    writeBudgetPolicy(write, policy);
    write(
      'gbdraw/web/js/app/editor.js',
      "import { runDiagramGeneration } from '../services/diagram-generation.js';\n"
        + 'export const editExistingOwner = () => runDiagramGeneration(1);\n'
    );
  }, [
    /production runtime files and Web guard\/CI files changed together/,
    /Diagram Worker: owner app\/editor\.js/,
    /proposed privileged capability allowlist excludes active owners or importers/
  ]);
});

test('runtime contraction cannot delete a capability key', () => {
  assertRevisionChangeBudgetFailure((write) => {
    const policy = cloneBudgetPolicy();
    delete policy.allowedPrivilegedOwners['Diagram Worker'];
    writeBudgetPolicy(write, policy);
    removeEditorOwnerUse(write);
  }, [
    /production runtime files and Web guard\/CI files changed together/,
    /allowedPrivilegedOwners\.Diagram Worker/,
    /proposed privileged capability policy is missing base allowlist keys/
  ]);
});

test('runtime contraction cannot delete a privileged import-target key', () => {
  assertRevisionChangeBudgetFailure((write) => {
    const policy = cloneBudgetPolicy();
    delete policy.allowedPrivilegedImporters['workers/diagram-generation-worker.js'];
    writeBudgetPolicy(write, policy);
    removeEditorOwnerUse(write);
  }, [
    /production runtime files and Web guard\/CI files changed together/,
    /allowedPrivilegedImporters\.workers\/diagram-generation-worker\.js/,
    /proposed privileged capability policy is missing base allowlist keys/
  ]);
});

test('runtime contraction cannot add privileged policy keys', () => {
  const keyCases = [
    {
      key: 'allowedPrivilegedOwners',
      name: 'Future capability',
      expected: /allowedPrivilegedOwners\.Future capability/
    },
    {
      key: 'allowedPrivilegedImporters',
      name: 'services/future-privileged-target.js',
      expected: /allowedPrivilegedImporters\.services\/future-privileged-target\.js/
    }
  ];
  keyCases.forEach(({ key, name, expected }) => {
    assertRevisionChangeBudgetFailure((write) => {
      const policy = cloneBudgetPolicy();
      removePolicyPath(
        policy,
        'allowedPrivilegedOwners',
        'Diagram Worker',
        'app/editor.js'
      );
      policy[key][name] = [];
      writeBudgetPolicy(write, policy);
      removeEditorOwnerUse(write);
    }, [
      /production runtime files and Web guard\/CI files changed together/,
      expected
    ]);
  });
});

test('runtime plus a semantically unchanged policy remains separated', () => {
  assertRevisionChangeBudgetFailure((write) => {
    const reorderedPolicy = {
      allowedPrivilegedOwners: BUDGET_POLICY.allowedPrivilegedOwners,
      allowedPrivilegedImporters: BUDGET_POLICY.allowedPrivilegedImporters
    };
    writeBudgetPolicy(write, reorderedPolicy, 0);
    write(
      'gbdraw/web/js/services/diagram-generation.js',
      'export const runDiagramGeneration = () => 2;\n'
    );
  }, [/production runtime files and Web guard\/CI files changed together/]);
});

test('runtime plus policy contraction rejects every additional guard change', () => {
  const guardCases = [
    ['checker', 'tools/check-web-change-budget.mjs'],
    ['promotion helper', 'tools/check-promotion-readiness.mjs'],
    ['architecture detectors', 'tools/web-architecture-detectors.mjs'],
    ['source parser', 'tools/web-change-source.mjs'],
    ['architecture contracts', 'tests/web/architecture-contracts.test.mjs'],
    ['promotion readiness tests', 'tests/web/promotion-readiness.test.mjs'],
    ['Gallery workflow', '.github/workflows/gallery-publication.yml'],
    ['deploy workflow', '.github/workflows/deploy_web.yml'],
    ['normal workflow', '.github/workflows/test.yml'],
    ['trusted-base workflow', '.github/workflows/web-base-policy.yml'],
    ['policy documentation', 'docs/internal/WEB_CHANGE_POLICY.md']
  ];
  guardCases.forEach(([name, path]) => {
    assertRevisionChangeBudgetFailure((write) => {
      applySingleOwnerContraction(write);
      const comment = path.endsWith('.mjs') ? '//' : '#';
      write(path, `${BUDGET_FIXTURE[path]}\n${comment} ${name} changed\n`);
    }, [/production runtime files and Web guard\/CI files changed together/]);
  });
});

test('safe contraction requires review above the selected production-file threshold', () => {
  BUDGET_PROFILES.forEach((profile) => {
    const result = runRevisionChangeBudgetCase((write) => {
      applySingleOwnerContraction(write);
      writeBudgetFiles(write, profile.productionFiles);
    }, profile.environment);
    assert.equal(result.status, 0, `${profile.name}\n${result.output}`);
    assert.match(result.output, /Gate: \*\*PASS\*\*/);
    assert.match(result.output, /Review: \*\*REQUIRED\*\*/);
    assert.match(
      result.output,
      new RegExp(
        `production files changed exceed ${profile.name} size-review threshold `
        + `\\(${profile.productionFiles + 1} > ${profile.productionFiles}\\)`
      )
    );
  });
});

test('safe contraction does not waive dependency, vendor, or binary integrity', () => {
  const cases = [
    {
      mutate: (write) => write(
        'package.json',
        '{"private":true,"dependencies":{"left-pad":"1.3.0"}}\n'
      ),
      expected: /new production dependencies are not allowed/
    },
    {
      mutate: (write) => write(
        'gbdraw/web/js/app/secondary.js',
        "import 'left-pad';\nexport const editSecondaryOwner = () => 1;\n"
      ),
      expected: /bare production import: gbdraw\/web\/js\/app\/secondary\.js: left-pad/
    },
    {
      mutate: (write) => write(
        'gbdraw/web/vendor/library.js',
        'globalThis.vendorLibrary = false;\n'
      ),
      expected: /changes under gbdraw\/web\/vendor\/ are not allowed/
    },
    {
      mutate: (write) => write(
        'gbdraw/web/js/runtime.bin',
        Buffer.from([0, 1, 2, 3])
      ),
      expected: /added binary runtime files are not allowed/
    }
  ];
  cases.forEach(({ mutate, expected }) => {
    assertRevisionChangeBudgetFailure((write) => {
      applySingleOwnerContraction(write);
      mutate(write);
    }, [expected]);
  });
});

test('unapproved privileged expansion fails against the base allowlist', () => {
  const mutate = (write) => {
    write(
      'gbdraw/web/js/app/secondary.js',
      "import { runDiagramGeneration } from '../services/diagram-generation.js';\n"
        + 'export const editSecondaryOwner = () => runDiagramGeneration();\n'
    );
  };
  assertNonWaivableWorkingTreeFailure(mutate, [
    /services\/diagram-generation\.js: importer app\/secondary\.js/,
    /privileged capability owners or importers exceed the base allowlist/
  ]);
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
    assert.match(result.output, /Gate: \*\*PASS\*\*/);
    assert.match(result.output, /Review: \*\*REQUIRED\*\*/);
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
    assert.match(preauthorization.output, /Gate: \*\*PASS\*\*/);
    assert.match(preauthorization.output, /Review: \*\*REQUIRED\*\*/);
    const preauthorizedBase = commit('preauthorize secondary owner');

    write(
      'gbdraw/web/js/app/secondary.js',
      "import { runDiagramGeneration } from '../services/diagram-generation.js';\n"
        + 'export const editSecondaryOwner = () => runDiagramGeneration();\n'
    );
    const implementationHead = commit('use preauthorized secondary owner');
    const implementation = execute({ base: preauthorizedBase, head: implementationHead });
    assert.equal(implementation.status, 0, implementation.output);
    assert.match(implementation.output, /Gate: \*\*PASS\*\*/);
    assert.match(implementation.output, /Review: \*\*REQUIRED\*\*/);

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
        + '  const template = `\nimport \'./session-file.js\';\n`;\n'
        + '  return { version: 1, journalToken: note + template };\n'
        + '};\n'
    );
  });
  assert.equal(result.status, 0, result.output);
  assert.match(result.output, /Report-only session object keys and compatibility names/);
  assert.match(result.output, /Review: \*\*REQUIRED\*\*/);
});

test('index.html growth counts toward the net-addition review threshold', () => {
  const result = runChangeBudgetCase((write) => {
    const additions = Array.from({ length: 110 }, (_, index) => `<div>${index}</div>`);
    write('gbdraw/web/index.html', `<main>baseline</main>\n${additions.join('\n')}\n`);
  });
  assert.equal(result.status, 0, result.output);
  assert.match(result.output, /Gate: \*\*PASS\*\*/);
  assert.match(result.output, /Review: \*\*REQUIRED\*\*/);
  assert.match(result.output, /M gbdraw\/web\/index\.html/);
  assert.match(
    result.output,
    /production net additions exceed ordinary size-review threshold \(110 > 100\)/
  );
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

test('trusted-base execution rejects self-authorization despite a successful head checker', () => {
  withChangeBudgetRepository(({ commit, execute, git, root, write }) => {
    const base = git(['rev-parse', 'HEAD']).stdout.trim();
    const policy = cloneBudgetPolicy();
    policy.allowedPrivilegedOwners['Diagram Worker'].push('app/secondary.js');
    writeBudgetPolicy(write, policy);
    write(
      'gbdraw/web/js/app/secondary.js',
      'export const editSecondaryOwner = () => runDiagramGeneration();\n'
    );
    write('tools/check-web-change-budget.mjs', 'process.exit(0);\n');
    const head = commit('malicious head checker self-authorizes expansion');

    const headChecker = spawnSync(
      process.execPath,
      [join(root, 'tools/check-web-change-budget.mjs'), '--base', base, '--head', head],
      { cwd: root, encoding: 'utf8', stdio: ['ignore', 'pipe', 'pipe'] }
    );
    assert.equal(headChecker.status, 0, `${headChecker.stdout}${headChecker.stderr}`);

    const trustedBase = execute({ base, head });
    assert.equal(trustedBase.status, 1, trustedBase.output);
    assert.match(
      trustedBase.output,
      /production runtime files and Web guard\/CI files changed together/
    );
    assert.match(
      trustedBase.output,
      /privileged capability owners or importers exceed the base allowlist/
    );
    assert.match(trustedBase.output, /Diagram Worker: owner app\/secondary\.js/);
  });
});

test('promotion omits runtime and guard inventories without weakening ordinary PRs', () => {
  withChangeBudgetRepository(({ commit, execute, git, root, write }) => {
    const base = git(['rev-parse', 'HEAD']).stdout.trim();
    write('gbdraw/web/js/app/editor.js', 'export const editExistingOwner = () => 2;\n');
    write('.github/workflows/test.yml', 'name: aggregated workflow\n');
    const head = commit('aggregate runtime and guard changes');

    const ordinary = execute({
      base,
      head,
      environment: {
        GITHUB_EVENT_NAME: '',
        GITHUB_EVENT_PATH: '',
        GITHUB_REPOSITORY: '',
        WEB_PROMOTION: 'true'
      }
    });
    assert.equal(ordinary.status, 1, ordinary.output);
    assert.match(ordinary.output, /Context: ORDINARY/);
    assert.match(
      ordinary.output,
      /production runtime files and Web guard\/CI files changed together/
    );

    const promotion = executePromotion({ execute, root, base, head });
    assert.equal(promotion.status, 0, promotion.output);
    assert.match(promotion.output, /Context: PROMOTION/);
    assert.match(promotion.output, /Gate: \*\*PASS\*\*/);
    assert.match(promotion.output, /Review: \*\*CLEAR\*\*/);
    assert.match(promotionReportSection(promotion.output, 'Promotion source coverage'), /Status: PASS/);
    assert.match(promotionReportSection(promotion.output, 'Promotion source coverage'), /Basis: DIRECT_ANCESTRY/);
    assert.doesNotMatch(promotion.output, /Production files touched|Architecture differential/);
    assert.equal(
      promotionReportSection(promotion.output, 'Blocking violations').trim(),
      '- None'
    );
  });
});

test('promotion omits checker and authority aggregation without authorizing ordinary PRs', () => {
  withChangeBudgetRepository(({ commit, execute, git, root, write }) => {
    const base = git(['rev-parse', 'HEAD']).stdout.trim();
    write(
      'tools/web-change-source.mjs',
      `${BUDGET_FIXTURE['tools/web-change-source.mjs']}\n// aggregate checker change\n`
    );
    write('.github/workflows/test.yml', 'name: aggregated authority\n');
    const head = commit('aggregate checker and authority changes');

    const ordinary = execute({ base, head });
    assert.equal(ordinary.status, 1, ordinary.output);
    assert.match(
      ordinary.output,
      /Web checker\/source parser and authority policy\/workflow files changed together/
    );

    const promotion = executePromotion({ execute, root, base, head });
    assert.equal(promotion.status, 0, promotion.output);
    assert.match(promotion.output, /Gate: \*\*PASS\*\*/);
    assert.match(promotion.output, /Review: \*\*CLEAR\*\*/);
    assert.doesNotMatch(promotion.output, /tools\/web-change-source\.mjs|\.github\/workflows\/test\.yml/);
  });
});

test('promotion omits ordinary size profiling while ordinary review remains required', () => {
  withChangeBudgetRepository(({ commit, execute, git, root, write }) => {
    const base = git(['rev-parse', 'HEAD']).stdout.trim();
    writeBudgetFiles(write, 9);
    write('gbdraw/web/js/app/budget-lines.js', budgetLineSource(401, 101));
    const head = commit('aggregate oversized production change');

    const ordinary = execute({ base, head });
    assert.equal(ordinary.status, 0, ordinary.output);
    assert.match(ordinary.output, /Review: \*\*REQUIRED\*\*/);
    const ordinarySizeReasons = promotionReportSection(ordinary.output, 'Review reasons');
    assert.match(
      ordinarySizeReasons,
      /production files changed exceed ordinary size-review threshold \(10 > 8\)/
    );
    assert.match(
      ordinarySizeReasons,
      /production gross churn exceeds ordinary size-review threshold \(912 > 800\)/
    );
    assert.match(
      ordinarySizeReasons,
      /production net additions exceed ordinary size-review threshold \(110 > 100\)/
    );
    assert.equal(
      promotionReportSection(ordinary.output, 'Blocking violations').trim(),
      '- None'
    );

    const promotion = executePromotion({ execute, root, base, head });
    assert.equal(promotion.status, 0, promotion.output);
    assert.match(promotion.output, /Review: \*\*CLEAR\*\*/);
    assert.doesNotMatch(promotion.output, /size-review threshold|Production files touched/);
    assert.equal(
      promotionReportSection(promotion.output, 'Blocking violations').trim(),
      '- None'
    );
  });
});

test('promotion relies on exact staging instead of rerunning ordinary import-cycle checks', () => {
  withChangeBudgetRepository(({ commit, execute, git, root, write }) => {
    const base = git(['rev-parse', 'HEAD']).stdout.trim();
    write(
      'gbdraw/web/js/app/editor.js',
      "import './secondary.js';\nexport const editExistingOwner = () => 1;\n"
    );
    write(
      'gbdraw/web/js/app/secondary.js',
      "import './editor.js';\nexport const editSecondaryOwner = () => 1;\n"
    );
    const head = commit('introduce import cycle');

    const ordinary = execute({ base, head });
    assert.equal(ordinary.status, 1, ordinary.output);
    assert.match(ordinary.output, /first-party static import cycles are not allowed/);

    const promotion = executePromotion({ execute, root, base, head });
    assert.equal(promotion.status, 0, promotion.output);
    assert.doesNotMatch(promotion.output, /first-party static import/);
  });
});

test('promotion relies on exact staging instead of rerunning ordinary privileged checks', () => {
  withChangeBudgetRepository(({ commit, execute, git, root, write }) => {
    const base = git(['rev-parse', 'HEAD']).stdout.trim();
    write(
      'gbdraw/web/js/app/secondary.js',
      "import { runDiagramGeneration } from '../services/diagram-generation.js';\n"
        + 'export const editSecondaryOwner = () => runDiagramGeneration();\n'
    );
    const head = commit('add unapproved privileged caller');

    const ordinary = execute({ base, head });
    assert.equal(ordinary.status, 1, ordinary.output);
    assert.match(
      ordinary.output,
      /privileged capability owners or importers exceed the base allowlist/
    );
    assert.match(ordinary.output, /Diagram Worker: owner app\/secondary\.js/);

    const promotion = executePromotion({ execute, root, base, head });
    assert.equal(promotion.status, 0, promotion.output);
    assert.doesNotMatch(promotion.output, /privileged capability/);
  });
});

test('promotion relies on exact staging instead of parsing ordinary architecture authority', () => {
  withChangeBudgetRepository(({ commit, execute, git, root, write }) => {
    const base = git(['rev-parse', 'HEAD']).stdout.trim();
    write('tools/web-architecture-rules.json', '{\n');
    const head = commit('malform architecture registry');

    const ordinary = execute({ base, head });
    assert.equal(ordinary.status, 1, ordinary.output);
    assert.match(ordinary.output, /candidate architecture rules: Cannot parse/);

    const promotion = executePromotion({ execute, root, base, head });
    assert.equal(promotion.status, 0, promotion.output);
    assert.doesNotMatch(promotion.output, /architecture rules/);
  });
});

test('arbitrary-head and fork pull requests to main fail the bounded promotion gate', () => {
  withChangeBudgetRepository(({ commit, execute, git, root, write }) => {
    const base = git(['rev-parse', 'HEAD']).stdout.trim();
    write('fixture-note.txt', 'safe source change\n');
    const head = commit('candidate source');
    const cases = [
      [
        'arbitrary head',
        { headRef: 'feature/example' },
        /Promotion pull requests must use dev as the head branch/
      ],
      [
        'fork head',
        { headRepository: 'contributor/gbdraw' },
        /Promotion pull requests must use dev from the current repository/
      ]
    ];
    cases.forEach(([name, event, reason]) => {
      const result = executePromotion({ execute, root, base, head, event });
      assert.equal(result.status, 1, `${name}\n${result.output}`);
      assert.match(result.output, /Context: PROMOTION/);
      assert.match(result.output, /Gate: \*\*FAIL\*\*/);
      assert.match(result.output, /Review: \*\*CLEAR\*\*/);
      assert.match(result.output, reason);
      assert.match(result.output, /Status: NOT_EVALUATED/);
    });
  });
});

test('malformed pull-request metadata fails before ordinary inventory evaluation', () => {
  withChangeBudgetRepository(({ commit, execute, git, root, write }) => {
    const base = git(['rev-parse', 'HEAD']).stdout.trim();
    write('fixture-note.txt', 'safe source change\n');
    const head = commit('candidate source');
    const eventPath = join(root, '.git', 'malformed-event.json');
    writeFileSync(eventPath, '{', 'utf8');
    const result = execute({
      base,
      head,
      environment: {
        GITHUB_EVENT_NAME: 'pull_request_target',
        GITHUB_EVENT_PATH: eventPath,
        GITHUB_REPOSITORY: 'satoshikawato/gbdraw'
      }
    });
    assert.equal(result.status, 1, result.output);
    assert.match(result.output, /Gate: \*\*FAIL\*\*/);
    assert.match(result.output, /GitHub pull request event payload is missing or malformed/);
    assert.doesNotMatch(result.output, /Architecture differential|Production files touched/);
  });
});

test('promotion accepts a content-neutral main merge commit without a sync PR', () => {
  withChangeBudgetRepository(({ commit, execute, git, root, write }) => {
    const common = git(['rev-parse', 'HEAD']).stdout.trim();
    write('gbdraw/web/js/app/editor.js', 'export const editExistingOwner = () => 2;\n');
    const promotedDevHead = commit('prior promoted dev head');

    assert.equal(git(['checkout', '--quiet', '--detach', common]).status, 0);
    assert.equal(
      git(['merge', '--quiet', '--no-ff', '-m', 'merge prior dev promotion', promotedDevHead]).status,
      0
    );
    const base = git(['rev-parse', 'HEAD']).stdout.trim();
    assert.equal(git(['merge-base', '--is-ancestor', base, promotedDevHead]).status, 1);
    assert.equal(
      git(['rev-parse', `${base}^{tree}`]).stdout.trim(),
      git(['rev-parse', `${promotedDevHead}^{tree}`]).stdout.trim()
    );

    assert.equal(git(['checkout', '--quiet', '--detach', promotedDevHead]).status, 0);
    write('gbdraw/web/js/app/secondary.js', 'export const editSecondaryOwner = () => 2;\n');
    const head = commit('continue dev after promotion');

    const promotion = executePromotion({ execute, root, base, head });
    assert.equal(promotion.status, 0, promotion.output);
    assert.match(
      promotionReportSection(promotion.output, 'Promotion source coverage'),
      /Status: PASS[\s\S]*Basis: MERGE_PARENT_TREE/
    );
    assert.equal(
      promotionReportSection(promotion.output, 'Blocking violations').trim(),
      '- None'
    );
  });
});

test('promotion rejects a main merge commit with main-only tree content', () => {
  withChangeBudgetRepository(({ commit, execute, git, root, write }) => {
    const common = git(['rev-parse', 'HEAD']).stdout.trim();
    write('gbdraw/web/js/app/editor.js', 'export const editExistingOwner = () => 2;\n');
    const promotedDevHead = commit('prior promoted dev head');

    assert.equal(git(['checkout', '--quiet', '--detach', common]).status, 0);
    assert.equal(
      git(['merge', '--quiet', '--no-ff', '-m', 'merge prior dev promotion', promotedDevHead]).status,
      0
    );
    write('.github/workflows/test.yml', 'name: main-only resolution\n');
    assert.equal(git(['add', '.github/workflows/test.yml']).status, 0);
    assert.equal(git(['commit', '--quiet', '--amend', '--no-edit']).status, 0);
    const base = git(['rev-parse', 'HEAD']).stdout.trim();
    assert.notEqual(
      git(['rev-parse', `${base}^{tree}`]).stdout.trim(),
      git(['rev-parse', `${promotedDevHead}^{tree}`]).stdout.trim()
    );

    assert.equal(git(['checkout', '--quiet', '--detach', promotedDevHead]).status, 0);
    write('gbdraw/web/js/app/secondary.js', 'export const editSecondaryOwner = () => 2;\n');
    const head = commit('continue dev without main-only content');

    const promotion = executePromotion({ execute, root, base, head });
    assert.equal(promotion.status, 1, promotion.output);
    assert.match(
      promotionReportSection(promotion.output, 'Promotion source coverage'),
      /Status: FAIL[\s\S]*Basis: MAIN_CONTENT_MISSING/
    );
    assert.match(
      promotion.output,
      /The promotion source does not contain the current main content\. Merge or rebase main into dev, then rerun the promotion\./
    );
  });
});

test('promotion fails when dev does not contain current main content', () => {
  withChangeBudgetRepository(({ commit, execute, git, root, write }) => {
    const common = git(['rev-parse', 'HEAD']).stdout.trim();
    write('gbdraw/web/js/app/editor.js', 'export const editExistingOwner = () => 2;\n');
    const head = commit('promotion source change');

    assert.equal(git(['checkout', '--quiet', '--detach', common]).status, 0);
    write('.github/workflows/test.yml', 'name: current main\n');
    const base = commit('advance main separately');

    const promotion = executePromotion({ execute, root, base, head });
    assert.equal(promotion.status, 1, promotion.output);
    assert.match(
      promotionReportSection(promotion.output, 'Promotion source coverage'),
      /Status: FAIL[\s\S]*Basis: MAIN_CONTENT_MISSING/
    );
    assert.match(
      promotion.output,
      /The promotion source does not contain the current main content\. Merge or rebase main into dev, then rerun the promotion\./
    );
  });
});

test('pull request SHA mismatches are blocking metadata errors', () => {
  withChangeBudgetRepository(({ commit, execute, git, root, write }) => {
    const base = git(['rev-parse', 'HEAD']).stdout.trim();
    write('gbdraw/web/js/app/editor.js', 'export const editExistingOwner = () => 2;\n');
    const head = commit('candidate change');

    const result = executePromotion({
      execute,
      root,
      base,
      head,
      event: { base: 'f'.repeat(40) }
    });
    assert.equal(result.status, 1, result.output);
    assert.match(result.output, /Context: PROMOTION/);
    assert.match(result.output, /Checker base SHA does not match the GitHub event payload/);
  });
});
