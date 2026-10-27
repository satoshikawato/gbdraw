import assert from 'node:assert/strict';
import { spawnSync } from 'node:child_process';
import { createHash } from 'node:crypto';
import {
  existsSync,
  mkdirSync,
  mkdtempSync,
  readFileSync,
  readdirSync,
  rmSync,
  symlinkSync,
  writeFileSync
} from 'node:fs';
import { tmpdir } from 'node:os';
import {
  dirname,
  extname,
  isAbsolute,
  join,
  posix,
  relative,
  resolve,
  sep
} from 'node:path';
import test from 'node:test';
import { fileURLToPath } from 'node:url';
import { knownJobsFor } from '../../tools/ci-impact-policy.mjs';
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
  join(REPOSITORY_ROOT, 'tests/web/playwright/functional.config.js'),
  'utf8'
);
const PR_SMOKE_PLAYWRIGHT_CONFIG = readFileSync(
  join(REPOSITORY_ROOT, 'tests/web/playwright/pr-smoke.config.js'),
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

const specifierFrom = (owner, target) => {
  const specifier = posix.relative(posix.dirname(owner), target);
  return specifier.startsWith('.') ? specifier : `./${specifier}`;
};

// These two request/Session callers may use the pure warning validator, not
// rendered-identity collection. The metadata admission callers remain exact.
const assertMetadataImportOwners = (sources) => {
  const target = 'services/session-feature-metadata.js';
  const validationCallers = new Set(['app/run-analysis.js', 'services/session-authority.js']);
  const identityCallers = [];
  for (const [owner, source] of sources) {
    const imports = literalImportSpecifiers(source).filter(specifier => {
      const resolved = resolveRelativeImport(join(JAVASCRIPT_ROOT, owner), specifier);
      return resolved && relativeModulePath(resolved) === target;
    });
    if (!imports.length) continue;
    if (!validationCallers.has(owner)) { identityCallers.push(owner); continue; }
    assert.equal(imports.length, 1, `${owner}: one warning validation import`);
    const code = maskJavaScript(source);
    const declarations = [...maskJavaScript(source, { strings: false }).matchAll(
      /(?:^|\n)\s*import\s*\{\s*validateAnnotationWarnings\s*\}\s*from\s*(['"])([^'"]+)\1/g
    )];
    assert.ok(declarations.some(match => match[2] === imports[0]
      && code.slice(match.index + match[0].indexOf('import')).startsWith('import')),
    `${owner}: only the named pure warning validator is permitted`);
  }
  assert.deepEqual(identityCallers.sort(), [
    'services/session-feature-recovery.js', 'services/svg-result-ingestion.js'
  ]);
};

// The Result replacement ceilings below count Result writes and SVG serialization
// together, as the single Mounted SVG/Result replacement capability did before G-J(1).
const RESULT_REPLACEMENT_CAPABILITIES = Object.freeze([
  'Result content commit',
  'SVG serialization'
]);
const operatorCount = (detected, path, capability) => (
  detected.operatorMatchesByCapability[capability]
    .find((entry) => entry.path === path)?.count || 0
);
const operatorEntries = (detected, path, capability) => {
  const count = operatorCount(detected, path, capability);
  return count ? [{ path, count }] : [];
};
const totalResultReplacementCount = (detected, path) => RESULT_REPLACEMENT_CAPABILITIES
  .reduce((total, capability) => total + operatorCount(detected, path, capability), 0);

// IN-01 and GE-02 (Web GUI audit 2026-09-30, D-04 and D-05) retire the live
// definition and stroke rewrites, and A1 of the post-audit follow-ups (W8 G-J)
// moves the editor Result commits into app/preview-runtime.js. These operators
// and imports may shrink or disappear without editing this guard; they may not
// grow, and no owner is added. Each ceiling counts both capabilities together.
const RETIRING_REPLACEMENT_CEILINGS = Object.freeze({
  'app/app-setup.js': 1,
  'app/feature-editor/color-actions.js': 1,
  'app/feature-editor/label-actions.js': 2,
  'app/feature-editor/svg-actions.js': 4,
  'app/legend-layout/canvas-actions.js': 2,
  'app/legend-layout/diagram-drag.js': 2,
  'app/legend-layout/reposition-actions.js': 2,
  'app/legend/drag-actions.js': 4,
  'app/legend/entry-actions.js': 5,
  'app/legend/sort-actions.js': 2,
  'app/legend/stroke-actions.js': 2,
  'app/results.js': 4,
  'app/svg-styles.js': 2
});
const retiringReplacementOperators = (detected, path, capability) => {
  const ceiling = RETIRING_REPLACEMENT_CEILINGS[path];
  assert.ok(totalResultReplacementCount(detected, path) <= ceiling,
    `${path}: retiring Result replacement may only contract (ceiling ${ceiling})`);
  return operatorEntries(detected, path, capability);
};
// The retiring owners that sort before `boundary`, in detector (path) order.
const retiringOperatorsBefore = (detected, boundary, capability) => Object.keys(
  RETIRING_REPLACEMENT_CEILINGS
).filter((path) => path < boundary)
  .flatMap((path) => retiringReplacementOperators(detected, path, capability));
// After A1, the Result switch and Save flushes of the remaining owners never
// write (T6 of the post-audit follow-ups). Their Result content commit counts
// may shrink without editing this guard; they may not grow. Generate's second
// commit (the label reflow) may leave app/run-analysis.js for the artifact
// owner of its Result's drawing (OV-111, D-15-49).
const RESULT_COMMIT_KEEPER_CEILINGS = Object.freeze({
  'app/preview-runtime.js': 4,
  'app/run-analysis.js': 2,
  'services/config.js': 4
});
// D6 (perf audit PR-9): the Save and History copies of the mounted root's
// serialization retire into the preview runtime's write on demand. Their
// counts may shrink or disappear without editing this guard; they may not grow.
const SERIALIZATION_RETIRING_CEILINGS = Object.freeze({
  'services/config.js': 1,
  'services/history-snapshot.js': 1
});
const retiringSerializationOperators = (detected, path) => {
  const ceiling = SERIALIZATION_RETIRING_CEILINGS[path];
  assert.ok(operatorCount(detected, path, 'SVG serialization') <= ceiling,
    `${path}: SVG serialization may only contract (ceiling ${ceiling})`);
  return operatorEntries(detected, path, 'SVG serialization');
};
const keeperCommitOperators = (detected, path) => {
  const ceiling = RESULT_COMMIT_KEEPER_CEILINGS[path];
  assert.ok(operatorCount(detected, path, 'Result content commit') <= ceiling,
    `${path}: Result content commit may only contract (ceiling ${ceiling})`);
  return operatorEntries(detected, path, 'Result content commit');
};
// DRAWINGS-A moves the History transactions out of the composition root into
// History's own owner (D-15-55): app/app-setup.js may only shed its History
// operators, and services/history.js may hold only what app/app-setup.js sheds.
const HISTORY_SHED_CEILING = 3;
const HISTORY_ROUTER_CEILING = 10;
const historyOperators = (detected, path) => {
  const shed = operatorCount(detected, 'app/app-setup.js', 'History');
  assert.ok(shed <= HISTORY_SHED_CEILING,
    `app/app-setup.js: History operators may only contract (ceiling ${HISTORY_SHED_CEILING})`);
  assert.ok(shed + operatorCount(detected, 'services/history.js', 'History') <= HISTORY_ROUTER_CEILING,
    `services/history.js: History operators may only take what app/app-setup.js sheds (ceiling ${HISTORY_ROUTER_CEILING} together)`);
  return operatorEntries(detected, path, 'History');
};
const importersWithRetiringOwners = (target, required, retiring) => {
  const actual = importersOf(target);
  const allowed = new Set([...required, ...retiring]);
  assert.deepEqual(actual.filter((path) => !allowed.has(path)), [],
    `${target}: no importer may be added`);
  return [...required, ...retiring.filter((path) => actual.includes(path))].sort();
};

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
      ...(productionSources.has('services/session-import-client.js')
        ? [['services/session-import-client.js', 1]] : []),
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
      ['services/feature-metadata-extraction.js', 2],
      ['services/diagram-generation.js', 1],
      ['workers/diagram-generation-worker.js', 2]
    ])
  );
  assert.deepEqual(importersOf('services/diagram-generation.js'), [
    'app/app-setup.js',
    'app/run-analysis.js',
    'services/feature-metadata-extraction.js',
    'services/record-discovery.js'
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
  const canonicalPayloadCall = /runDiagramGeneration\s*\(\s*\{\s*request:\s*canonical\.renderRequest,\s*resources:\s*canonical\.resources\s*\}\s*(?=\)|,\s*\{\s*onProgress\s*:)/s;
  assert.match(productionSources.get('app/run-analysis.js'), canonicalPayloadCall);
  for (const ending of [')', ', { onProgress: observer })']) {
    assert.match(
      `runDiagramGeneration({ request: canonical.renderRequest, resources: canonical.resources }${ending}`,
      canonicalPayloadCall
    );
  }
  for (const call of [
    'runDiagramGeneration({ request: alternateRequest, resources: canonical.resources })',
    'runDiagramGeneration({ request: canonical.renderRequest, resources: alternateResources })',
    'runDiagramGeneration({ request: canonical.renderRequest, resources: canonical.resources, argv: [] })',
    'runDiagramGeneration({ request: canonical.renderRequest, resources: canonical.resources }, alternatePayload)'
  ]) {
    assert.doesNotMatch(call, canonicalPayloadCall);
  }
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
    'app/watchers.js',
    'services/config.js',
    'services/session-feature-recovery.js',
    'state.js'
  ]);
  assertMetadataImportOwners(productionSources);
  assert.deepEqual(importersOf('services/svg-result-normalization.js'), importersWithRetiringOwners(
    'services/svg-result-normalization.js',
    ['services/config.js', 'services/svg-result-ingestion.js'],
    ['app/svg-styles.js']
  ));
  assert.doesNotMatch(
    productionSources.get('services/session-feature-recovery.js'),
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
    /normalizeLegacySvg/
  );
  // An older Session's Result gets the legacy composition once, in the
  // admission transform of services/config.js after its legend entry groups
  // (and, in the stroke form, before its strokes): services/config.js calls the
  // composition root's transformLegacyResultSvg port (R13 layering plan D4 B1),
  // and app/app-setup.js calls the owners, which stay in app/. In the
  // stroke-free form (edit-unify U2a) the saved bytes show the strokes, the
  // executor is the only writer of a live stroke edit, and no stroke
  // projection or reapply is left anywhere.
  const configSource = productionSources.get('services/config.js');
  const setupSource = productionSources.get('app/app-setup.js');
  const configImports = directImports.get('services/config.js');
  const strokeFree = occurrenceOwners(
    /\b(?:applyStrokeOverridesToSvg|reapplyStrokeOverrides|reconcileStrokeOverrides)\b/g
  ).size === 0;
  const legacyOwnerModules = [
    'app/legend-layout/composition-actions.js',
    ...(strokeFree ? [] : ['app/legend/stroke-actions.js'])
  ];
  const legacyCompositionCalls = occurrenceOwners(/\bnormalizeLegacyComposition\(svg,/g);
  const legacyStrokeCalls = occurrenceOwners(/\bapplyStrokeOverridesToSvg\(\{/g);
  const strokeCall = strokeFree ? '' : String.raw`[\s\S]+applyStrokeOverridesToSvg\(\{`;
  assert.deepEqual(
    occurrenceOwners(/\bexport const (?:normalizeLegacyComposition|applyStrokeOverridesToSvg)\b/g),
    new Map(legacyOwnerModules.map((module) => [module, 1]))
  );
  assert.match(configSource, /normalizeLegacyLegendEntryGroups\(svg\)[\s\S]+transformLegacyResultSvg\(svg,/);
  assert.deepEqual(legacyOwnerModules.filter((module) => configImports.has(module)), []);
  assert.doesNotMatch(configSource, /\bnormalizeLegacyComposition\b|\bapplyStrokeOverridesToSvg\b/);
  assert.match(
    setupSource,
    new RegExp(String.raw`export const transformLegacyResultSvg\b[\s\S]+COMPOSITION_SCHEMA_ATTRIBUTE\) === null[\s\S]+COMPOSITION_METADATA_ATTRIBUTE\) === null[\s\S]+normalizeLegacyComposition\(svg,` + strokeCall)
  );
  assert.deepEqual(legacyOwnerModules.filter((module) => !directImports.get('app/app-setup.js').has(module)), []);
  assert.deepEqual(legacyCompositionCalls, new Map([['app/app-setup.js', 1]]));
  assert.equal(legacyStrokeCalls.get('app/app-setup.js'), strokeFree ? undefined : 1);
  assert.equal(legacyStrokeCalls.has('services/config.js'), false);
  assert.deepEqual(occurrenceOwners(/\btransformLegacyResultSvg\(svg,/g), new Map([['services/config.js', 1]]));
  assert.deepEqual(occurrenceOwners(/\bbuildHistorySnapshot\b/g), new Map());
  assert.deepEqual(occurrenceOwners(/\bapplyHistorySnapshot\b/g), new Map());
});

test('right drawer availability and transitions have one production owner', () => {
  const transitionOwners = [
    ...occurrenceOwners(/\b(?:showRightDrawer|rightDrawerTab)\.value\s*=(?!=)/g).keys()
  ];

  // Exactly one module writes the drawer state: services/right-drawer-state.js
  // (R13 layering slice S8). The controller stays in app/right-drawer.js.
  assert.deepEqual(transitionOwners, ['services/right-drawer-state.js']);
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

// Capability keys that the policy holds while no detector reports them, during a
// split ("Splitting a capability" in docs/internal/WEB_CHANGE_POLICY.md).
// None is pending: the G-J(1) split of Mounted SVG/Result replacement is complete.
const UNDETECTED_POLICY_CAPABILITY_KEYS = Object.freeze([]);
const policyCapabilityKeysMatchDetectors = (policy) => assert.deepEqual(
  Object.keys(policy.allowedPrivilegedOwners).sort(),
  [...WEB_PRIVILEGED_CAPABILITY_KEYS, ...UNDETECTED_POLICY_CAPABILITY_KEYS].sort()
);

test('privileged capability owners and importers stay within their allowlists', () => {
  const detected = detectPrivilegedWebCapabilities(productionSources);
  const assertSubset = (actual, allowed, label) => {
    const allowedSet = new Set(allowed);
    assert.deepEqual(actual.filter((path) => !allowedSet.has(path)), [], label);
  };

  Object.entries(WEB_CHANGE_POLICY.allowedPrivilegedImporters).forEach(([target, allowed]) => {
    assertSubset(detected.importersByTarget[target], allowed, `${target} importers`);
  });
  policyCapabilityKeysMatchDetectors(WEB_CHANGE_POLICY);
  assert.deepEqual(
    Object.keys(WEB_CHANGE_POLICY.allowedPrivilegedImporters).sort(),
    WEB_PRIVILEGED_IMPORT_TARGETS
  );
  Object.entries(detected.operatorMatchesByCapability).forEach(([capability, matches_]) => {
    const actual = matches_.map(({ path }) => path);
    assertSubset(actual, WEB_CHANGE_POLICY.allowedPrivilegedOwners[capability], capability);
  });
});

test('the policy capability key set is exact', () => {
  const withPolicyKeys = (keys) => ({
    allowedPrivilegedOwners: Object.fromEntries(keys.map((key) => [key, []]))
  });
  const expectedKeys = [...WEB_PRIVILEGED_CAPABILITY_KEYS, ...UNDETECTED_POLICY_CAPABILITY_KEYS];
  policyCapabilityKeysMatchDetectors(withPolicyKeys(expectedKeys));
  assert.throws(() => policyCapabilityKeysMatchDetectors(
    withPolicyKeys([...expectedKeys, 'Unregistered capability'])
  ));
  assert.throws(() => policyCapabilityKeysMatchDetectors(
    withPolicyKeys(expectedKeys.filter((key) => key !== 'Result content commit'))
  ));
  assert.throws(() => policyCapabilityKeysMatchDetectors(
    withPolicyKeys([...expectedKeys, 'Mounted SVG/Result replacement'])
  ));
});

const privilegedOwnerViolations = (sources) => {
  const detected = detectPrivilegedWebCapabilities(sources);
  return Object.entries(detected.operatorMatchesByCapability).flatMap(([capability, matches_]) => {
    const allowed = new Set(WEB_CHANGE_POLICY.allowedPrivilegedOwners[capability]);
    return matches_.filter(({ path }) => !allowed.has(path))
      .map(({ path }) => `${capability}: owner ${path}`);
  });
};
const withSeededLine = (path, line) => new Map([
  ...productionSources,
  [path, `${productionSources.get(path) ?? ''}\n${line}\n`]
]);

test('Result content commit and SVG serialization detect different operators', () => {
  const detected = detectPrivilegedWebCapabilities(new Map([
    ['app/seeded-result-writer.js', [
      'results.value = next;',
      'state.results.value[index] = next;',
      'previewRuntime.flushActiveResult();',
      'if (results.value === next) previous = results.value;',
      '// results.value = comment;',
      'const text = "serializeCleanSvg(svg)";'
    ].join('\n')],
    ['app/seeded-serializer.js', 'const content = serializeCleanSvg(svg);\n']
  ]));
  assert.deepEqual(detected.operatorMatchesByCapability['Result content commit'],
    [{ path: 'app/seeded-result-writer.js', count: 3 }]);
  assert.deepEqual(detected.operatorMatchesByCapability['SVG serialization'],
    [{ path: 'app/seeded-serializer.js', count: 1 }]);
});

test('a Result writer or SVG serializer outside its allowlist fails the policy', () => {
  assert.deepEqual(privilegedOwnerViolations(productionSources), []);
  const unlisted = 'app/feature-search/preview-svg.js';
  assert.ok(productionSources.has(unlisted), unlisted);
  for (const [line, capability] of [
    ['results.value = nextResults;', 'Result content commit'],
    ['state.results.value[0] = nextResult;', 'Result content commit'],
    ['previewRuntime.flushActiveResult();', 'Result content commit'],
    ['const content = serializeCleanSvg(svg);', 'SVG serialization']
  ]) {
    assert.deepEqual(privilegedOwnerViolations(withSeededLine(unlisted, line)),
      [`${capability}: owner ${unlisted}`], line);
  }
  assert.deepEqual(privilegedOwnerViolations(
    withSeededLine('app/preview-runtime.js', 'state.results.value = nextResults;')
  ), []);
  assert.deepEqual(privilegedOwnerViolations(
    withSeededLine('services/svg-result-ingestion.js', 'const content = serializeCleanSvg(svg);')
  ), []);
});

// R1 (gbdraw/web/CLAUDE.md): the Result content commit list only contracts.
// Removing an owner needs no edit here; adding one fails.
const RESULT_CONTENT_COMMIT_OWNER_CEILING = Object.freeze([
  'app/app-setup.js',
  'app/feature-editor/color-actions.js',
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
  'app/run-analysis.js',
  'app/svg-styles.js',
  'services/config.js'
]);
const assertResultCommitOwnersOnlyContract = (owners) => assert.deepEqual(
  owners.filter((path) => !RESULT_CONTENT_COMMIT_OWNER_CEILING.includes(path)),
  [],
  'Result content commit owners may only be removed (R1)'
);

test('the Result content commit allowlist can only contract', () => {
  assertResultCommitOwnersOnlyContract(
    WEB_CHANGE_POLICY.allowedPrivilegedOwners['Result content commit']
  );
  assertResultCommitOwnersOnlyContract(RESULT_CONTENT_COMMIT_OWNER_CEILING.slice(1));
  assertResultCommitOwnersOnlyContract([]);
  for (const added of ['services/history-snapshot.js', 'app/feature-search/preview-svg.js']) {
    assert.throws(() => assertResultCommitOwnersOnlyContract(
      [...RESULT_CONTENT_COMMIT_OWNER_CEILING, added]
    ), /only be removed/, added);
  }
});

test('each Result write allowlist admits only its own kind of operator', () => {
  const owners = WEB_CHANGE_POLICY.allowedPrivilegedOwners;
  const commitOwners = owners['Result content commit'];
  const serializationOwners = owners['SVG serialization'];
  const seededViolations = (path, line) => privilegedOwnerViolations(
    new Map([[path, `${productionSources.get(path) ?? ''}\n${line}\n`]])
  );
  const commitLine = 'state.results.value = nextResults;';
  const serializationLine = 'const content = serializeCleanSvg(svg);';

  // The shared Mounted SVG/Result replacement list let all of these owners do both.
  const crossKindCases = [
    ...serializationOwners.filter((path) => !commitOwners.includes(path))
      .map((path) => [path, commitLine, 'Result content commit']),
    ...commitOwners.filter((path) => !serializationOwners.includes(path))
      .map((path) => [path, serializationLine, 'SVG serialization'])
  ];
  assert.ok(crossKindCases.length > 0);
  for (const [path, line, capability] of crossKindCases) {
    assert.deepEqual(seededViolations(path, line), [`${capability}: owner ${path}`],
      `${path}: ${line}`);
  }
  commitOwners.forEach((path) => assert.deepEqual(seededViolations(path, commitLine), [], path));
  serializationOwners.forEach((path) => (
    assert.deepEqual(seededViolations(path, serializationLine), [], path)
  ));
});

test('shared privileged detectors preserve the characterized current-source facts', () => {
  assert.deepEqual(WEB_PRIVILEGED_CAPABILITY_KEYS, [
    'Canonical editor state',
    'Diagram Worker',
    'History',
    'Python helper',
    'Render request',
    'Resource staging',
    'Result content commit',
    'SVG serialization',
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
      ...(productionSources.has('services/session-import-client.js')
        ? [{ path: 'services/session-import-client.js', count: 1 }] : []),
      { path: 'workers/diagram-generation-worker.js', count: 1 },
      { path: 'workers/losat-threaded-worker.js', count: 2 }
    ],
    History: [
      ...historyOperators(detected, 'app/app-setup.js'),
      { path: 'services/history-snapshot.js', count: 2 },
      ...historyOperators(detected, 'services/history.js')
    ],
    'Result content commit': [
      ...retiringOperatorsBefore(detected, 'app/preview-runtime.js', 'Result content commit'),
      ...keeperCommitOperators(detected, 'app/preview-runtime.js'),
      ...retiringReplacementOperators(detected, 'app/results.js', 'Result content commit'),
      ...keeperCommitOperators(detected, 'app/run-analysis.js'),
      ...retiringReplacementOperators(detected, 'app/svg-styles.js', 'Result content commit'),
      ...keeperCommitOperators(detected, 'services/config.js')
    ],
    'SVG serialization': [
      ...retiringOperatorsBefore(detected, 'services/config.js', 'SVG serialization'),
      ...retiringSerializationOperators(detected, 'services/config.js'),
      ...retiringSerializationOperators(detected, 'services/history-snapshot.js'),
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

test('metadata validation callers cannot acquire identity or namespace access', () => {
  const recovery = 'services/session-feature-recovery.js';
  const admission = new Map([
    [recovery, `import { normalizeRenderedFeatureId } from '${specifierFrom(recovery, 'services/session-feature-metadata.js')}';`],
    ['services/svg-result-ingestion.js', "import { collectRenderedFeatureIdentitiesFromSvgRoot } from './session-feature-metadata.js';"]
  ]);
  assertMetadataImportOwners(admission);
  for (const owner of ['app/run-analysis.js', 'services/session-authority.js']) {
    const specifier = owner.startsWith('app/') ? '../services/session-feature-metadata.js' : './session-feature-metadata.js';
    const valid = `import { validateAnnotationWarnings } from '${specifier}';`;
    assertMetadataImportOwners(new Map([...admission, [owner, valid]]));
    for (const invalid of [
      `import { collectRenderedFeatureIdentitiesFromSvgRoot } from '${specifier}';`,
      `import { validateAnnotationWarnings, normalizeRenderedFeatureId } from '${specifier}';`,
      `import * as metadata from '${specifier}';`,
      `import metadata from '${specifier}';`,
      `const metadata = await import('${specifier}');`,
      `${valid}\nimport { normalizeRenderedFeatureId } from '${specifier}';`,
      `const fake = \`${valid}\`;\nimport metadata from '${specifier}';`
    ]) assert.throws(() => assertMetadataImportOwners(new Map([...admission, [owner, invalid]])));
  }
  assert.throws(() => assertMetadataImportOwners(new Map([...admission,
    ['app/unapproved.js', "import { validateAnnotationWarnings } from '../services/session-feature-metadata.js';"]
  ])));
  assert.throws(() => assertMetadataImportOwners(new Map([...admission].slice(1))));
});

test('legend replacement characterization permits contraction but rejects growth', () => {
  const path = 'app/legend/entry-actions.js';
  for (const [commits, serializations] of [[4, 0], [5, 0], [2, 2], [3, 2], [6, 0], [3, 3], [0, 0]]) {
    const code = 'results.value = candidate;\n'.repeat(commits)
      + 'serializeCleanSvg(svg);\n'.repeat(serializations)
      + '// results.value = comment;\nconst text = "results.value = string;";';
    const detected = detectPrivilegedWebCapabilities(new Map([[path, code]]));
    assert.deepEqual(detected.operatorMatchesByCapability['Result content commit'],
      commits ? [{ path, count: commits }] : []);
    assert.deepEqual(detected.operatorMatchesByCapability['SVG serialization'],
      serializations ? [{ path, count: serializations }] : []);
    const total = commits + serializations;
    if (total <= 5) {
      assert.deepEqual(retiringReplacementOperators(detected, path, 'Result content commit'),
        commits ? [{ path, count: commits }] : []);
    } else {
      assert.throws(() => retiringReplacementOperators(detected, path, 'Result content commit'),
        /may only contract/);
    }
  }
});

test('retiring Result replacement owners may contract or disappear but never grow', () => {
  for (const [path, ceiling] of Object.entries(RETIRING_REPLACEMENT_CEILINGS)) {
    for (const count of [0, 1, ceiling, ceiling + 1]) {
      const code = 'results.value = candidate;\n'.repeat(count);
      const detected = detectPrivilegedWebCapabilities(new Map([[path, code]]));
      if (count <= ceiling) {
        assert.deepEqual(retiringReplacementOperators(detected, path, 'Result content commit'),
          count ? [{ path, count }] : []);
      } else {
        assert.throws(
          () => retiringReplacementOperators(detected, path, 'Result content commit'),
          /may only contract/
        );
      }
    }
    const split = detectPrivilegedWebCapabilities(new Map([[
      path,
      'results.value = candidate;\n'.repeat(ceiling) + 'serializeCleanSvg(svg);\n'
    ]]));
    assert.throws(() => retiringReplacementOperators(split, path, 'SVG serialization'),
      /may only contract/);
  }
  const required = importersOf('services/diagram-generation.js').filter((path) => path !== 'app/results.js');
  assert.throws(
    () => importersWithRetiringOwners('services/diagram-generation.js', required.slice(1), ['app/results.js']),
    /no importer may be added/
  );
});

test('kept Result content commit owners may contract but never grow', () => {
  for (const [path, ceiling] of Object.entries(RESULT_COMMIT_KEEPER_CEILINGS)) {
    for (const count of [0, 1, ceiling, ceiling + 1]) {
      const code = 'results.value = candidate;\n'.repeat(count) + 'serializeCleanSvg(svg);\n';
      const detected = detectPrivilegedWebCapabilities(new Map([[path, code]]));
      if (count <= ceiling) {
        assert.deepEqual(keeperCommitOperators(detected, path), count ? [{ path, count }] : []);
      } else {
        assert.throws(() => keeperCommitOperators(detected, path), /may only contract/);
      }
    }
  }
});

test('History operators may move from the composition root into History but never grow', () => {
  const historyCode = (count) => 'beginCheckpoint(label);\n'.repeat(count);
  for (const [shed, owner] of [[3, 7], [0, 10], [2, 8], [0, 0], [3, 0], [4, 6], [3, 8], [0, 11]]) {
    const detected = detectPrivilegedWebCapabilities(new Map([
      ['app/app-setup.js', historyCode(shed)], ['services/history.js', historyCode(owner)]
    ]));
    assert.equal(operatorCount(detected, 'app/app-setup.js', 'History'), shed);
    assert.equal(operatorCount(detected, 'services/history.js', 'History'), owner);
    if (shed <= HISTORY_SHED_CEILING && shed + owner <= HISTORY_ROUTER_CEILING) {
      assert.deepEqual(historyOperators(detected, 'services/history.js'),
        owner ? [{ path: 'services/history.js', count: owner }] : []);
    } else {
      assert.throws(() => historyOperators(detected, 'services/history.js'), /may only (contract|take)/);
    }
  }
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
    'canonical-path.current-result-admission.v1',
    'canonical-path.render-request.v1',
    'semantic-owner.current-result-admission.v1',
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

test('current Result admission detectors characterize the merged runtime source facts', () => {
  const semantic = WEB_ARCHITECTURE_DETECTORS[
    'semantic-owner.current-result-admission.v1'
  ];
  const canonical = WEB_ARCHITECTURE_DETECTORS[
    'canonical-path.current-result-admission.v1'
  ];

  assert.deepEqual(semantic.detect(productionSources), {
    definitionCount: 1,
    observedDefinitions: [{
      path: 'services/svg-result-ingestion.js',
      count: 1,
      subject: 'services/svg-result-ingestion.js'
    }],
    subjects: ['services/svg-result-ingestion.js']
  });
  assert.deepEqual(canonical.detect(productionSources), {
    observedEdges: [{
      from: 'app/candidate-render.js',
      to: 'services/svg-result-ingestion.js',
      subject: 'app/candidate-render.js -> services/svg-result-ingestion.js'
    }],
    subjects: ['app/candidate-render.js -> services/svg-result-ingestion.js']
  });
});

test('declared architecture rules actively govern the current source tree', () => {
  assert.deepEqual(validateArchitectureRuleRegistry(
    WEB_ARCHITECTURE_RULES,
    WEB_ARCHITECTURE_DETECTORS,
    { availableEnforcements: ['hard', 'report-only'] }
  ), { valid: true, errors: [] });

  const ruleKeys = WEB_ARCHITECTURE_RULES.rules.map(({ key }) => key);
  assert.deepEqual(ruleKeys, [...ruleKeys].sort());

  const evaluateRules = (sources) => WEB_ARCHITECTURE_RULES.rules.map((rule) => {
    const detector = WEB_ARCHITECTURE_DETECTORS[rule.detector];
    assert.ok(detector, `Active rule ${rule.key} must have a trusted detector`);
    const observed = classifyArchitectureRuleObservation(rule, detector.detect(sources));
    return {
      key: rule.key,
      enforcement: rule.enforcement,
      subjects: observed.subjects,
      ...evaluateArchitectureRuleResult({
        observation: observed.observation,
        mode: rule.enforcement.replace('-', '_').toUpperCase(),
        baselineRelation: 'NOT_APPLICABLE',
        authorityResolution: 'NOT_APPLICABLE'
      })
    };
  });
  const results = evaluateRules(productionSources);

  assert.deepEqual(results.map(({ key }) => key), ruleKeys);
  assert.deepEqual(
    results
      .filter(({ enforcement }) => enforcement === 'hard')
      .map(({ key, observation, decision }) => ({ key, observation, decision })),
    WEB_ARCHITECTURE_RULES.rules
      .filter(({ enforcement }) => enforcement === 'hard')
      .map(({ key }) => ({ key, observation: 'CONFORMING', decision: 'PASS' }))
  );
  assert.deepEqual(
    evaluateRules(new Map([...productionSources].reverse())),
    results
  );
});

test('four-rule authority activation validates and passes on the current source tree', () => {
  const activationReadinessRegistry = {
    schemaVersion: 1,
    rules: [
      {
        key: 'canonical-path.current-result-admission',
        kind: 'single-canonical-entry-edge',
        detector: 'canonical-path.current-result-admission.v1',
        allowedEdges: [{
          from: 'app/candidate-render.js',
          to: 'services/svg-result-ingestion.js'
        }],
        exactActiveEdgeCount: 1,
        enforcement: 'hard',
        baselineEligible: false
      },
      {
        key: 'canonical-path.render-request',
        kind: 'single-canonical-entry-edge',
        detector: 'canonical-path.render-request.v1',
        allowedEdges: [{
          from: 'app/run-analysis.js',
          to: 'services/session-request.js'
        }],
        exactActiveEdgeCount: 1,
        enforcement: 'hard',
        baselineEligible: false
      },
      {
        key: 'semantic-owner.current-result-admission',
        kind: 'single-semantic-owner',
        detector: 'semantic-owner.current-result-admission.v1',
        allowedDefinitionPaths: ['services/svg-result-ingestion.js'],
        exactDefinitionCount: 1,
        enforcement: 'hard',
        baselineEligible: false
      },
      {
        key: 'semantic-owner.render-request',
        kind: 'single-semantic-owner',
        detector: 'semantic-owner.render-request.v1',
        allowedDefinitionPaths: ['services/session-request.js'],
        exactDefinitionCount: 1,
        enforcement: 'hard',
        baselineEligible: false
      }
    ]
  };
  assert.deepEqual(validateArchitectureRuleRegistry(
    activationReadinessRegistry,
    WEB_ARCHITECTURE_DETECTORS,
    { availableEnforcements: ['hard', 'report-only'] }
  ), { valid: true, errors: [] });

  const evaluateRegistry = (sources) => activationReadinessRegistry.rules.map((rule) => {
    const detector = WEB_ARCHITECTURE_DETECTORS[rule.detector];
    assert.ok(detector, `Activation rule ${rule.key} must have a trusted detector`);
    const observed = classifyArchitectureRuleObservation(rule, detector.detect(sources));
    const result = evaluateArchitectureRuleResult({
      observation: observed.observation,
      mode: 'HARD',
      baselineRelation: 'NOT_APPLICABLE',
      authorityResolution: 'NOT_APPLICABLE'
    });
    return {
      key: rule.key,
      detector: rule.detector,
      subjects: observed.subjects,
      observation: result.observation,
      decision: result.decision
    };
  });
  const results = evaluateRegistry(productionSources);

  assert.deepEqual(results, [
    {
      key: 'canonical-path.current-result-admission',
      detector: 'canonical-path.current-result-admission.v1',
      subjects: ['app/candidate-render.js -> services/svg-result-ingestion.js'],
      observation: 'CONFORMING',
      decision: 'PASS'
    },
    {
      key: 'canonical-path.render-request',
      detector: 'canonical-path.render-request.v1',
      subjects: ['app/run-analysis.js -> services/session-request.js'],
      observation: 'CONFORMING',
      decision: 'PASS'
    },
    {
      key: 'semantic-owner.current-result-admission',
      detector: 'semantic-owner.current-result-admission.v1',
      subjects: ['services/svg-result-ingestion.js'],
      observation: 'CONFORMING',
      decision: 'PASS'
    },
    {
      key: 'semantic-owner.render-request',
      detector: 'semantic-owner.render-request.v1',
      subjects: ['services/session-request.js'],
      observation: 'CONFORMING',
      decision: 'PASS'
    }
  ]);
  assert.deepEqual(
    evaluateRegistry(new Map([...productionSources].reverse())),
    results
  );
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

test('source masking follows template expressions as code', () => {
  const lines = (source) => source.split('\n');
  const declarationsOf = (masked) => [...masked.matchAll(/^const\s+(\w+)/gm)].map(match => match[1]);

  // A nested template inside `${ ... }` must not end the outer template early.
  const nested = 'const a = `x ${ `y ${ 1 }` } z`;\nconst b = 2;\n';
  assert.equal(maskJavaScript(nested).split('\n')[0].replace(/\s/g, ''), 'consta=1;');
  assert.deepEqual(declarationsOf(maskJavaScript(nested)), ['a', 'b']);

  // A quote inside a regex literal inside an expression does not open a string.
  const regex = "const a = `${ s.replace(/'/g, '') }`;\nconst b = 2;\n";
  assert.deepEqual(declarationsOf(maskJavaScript(regex)), ['a', 'b']);
  assert.match(maskJavaScript(regex), /s\.replace\(/);

  // Braces of object literals inside an expression do not end it early.
  const braces = 'const a = `${ f({ k: { v: 1 } }) } tail`;\nconst b = 2;\n';
  const bracesMasked = maskJavaScript(braces);
  assert.deepEqual(declarationsOf(bracesMasked), ['a', 'b']);
  assert.ok(!bracesMasked.includes('tail'));
  assert.ok(bracesMasked.includes('f({ k: { v: 1 } })'));

  // Comments inside an expression are comments, including a `}` or backtick.
  const comment = 'const a = `${ f( // ` }\n  1) } tail`;\nconst b = `${ g(/* } ` */ 2) }`;\nconst c = 3;\n';
  const commentMasked = maskJavaScript(comment, { strings: false });
  assert.deepEqual(declarationsOf(commentMasked), ['a', 'b', 'c']);
  assert.ok(!commentMasked.includes('// `'));
  assert.ok(!commentMasked.includes('/* }'));
  assert.deepEqual(declarationsOf(maskJavaScript(comment)), ['a', 'b', 'c']);

  // Template text stays hidden under strings: true and visible otherwise.
  assert.ok(!maskJavaScript('const a = `watch(cache) ${ 1 }`;').includes('watch'));
  assert.ok(maskJavaScript('const a = `watch(cache) ${ 1 }`;', { strings: false }).includes('watch'));

  // The run-info.js shell quoting line must not desynchronize what follows.
  const quoting = [
    'export const quoteShellArg = (value) => {',
    "  return `'${token.replace(/'/g, `'\\\\''`)}'`;",
    '};',
    '',
    '// comment after the template',
    'export const buildShellCommand = (tokens) => tokens.length;',
    ''
  ].join('\n');
  const quotingCode = lines(maskJavaScript(quoting));
  const quotingComments = lines(maskJavaScript(quoting, { strings: false }));
  assert.equal(quotingCode[5], quoting.split('\n')[5]);
  assert.equal(quotingComments[4].trim(), '');
  assert.match(quotingCode[1], /token\.replace\(/);
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
  // Pull request and dispatch runs replace older runs; dev push staging runs finish.
  assert.match(TEST_WORKFLOW, /\n  cancel-in-progress: \$\{\{ github\.event_name != 'push' \}\}\n/);

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
  assert.match(
    GALLERY_PUBLICATION_WORKFLOW,
    /group: gallery-publication-\$\{\{ github\.event_name \}\}-\$\{\{ github\.ref \}\}/
  );
  // Like Tests, a started Gallery push run finishes so its SHA keeps exact evidence.
  assert.match(GALLERY_PUBLICATION_WORKFLOW, /\n  cancel-in-progress: \$\{\{ github\.event_name != 'push' \}\}\n/);
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
  assert.match(planner, /Build dev CI impact plan[\s\S]*node tools\/ci-impact\.mjs plan/);
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
  assert.match(webPrSmoke, /\n    timeout-minutes: 20\n/);
  assert.equal([...webPrSmoke.matchAll(/uses: actions\/checkout@/g)].length, 1);
  assert.equal([...webPrSmoke.matchAll(/npm ci/g)].length, 1);
  assert.equal([...webPrSmoke.matchAll(/playwright install --with-deps chromium/g)].length, 1);
  assert.equal([...webPrSmoke.matchAll(/python tools\/prepare_browser_wheel\.py/g)].length, 1);
  const webContracts = workflowJob('web-contracts-pr');
  assert.match(webContracts, /Run fast Web JavaScript contracts/);
  assert.match(webContracts, /-m "browser and not slow"/);
  assert.doesNotMatch(webPrSmoke, /Run fast Web JavaScript contracts|python -m pytest/);
  assert.match(webPrSmoke, /npm run test:web:pr-smoke/);
  assert.match(
    webPrSmoke,
    /Prepare browser wheel[\s\S]+Run Playwright PR smoke[\s\S]+npm run test:web:pr-smoke/
  );
  assert.doesNotMatch(webPrSmoke, /test:web:gallery-publication/);
  assert.match(
    workflowJob('gallery'),
    /Prepare Gallery browser wheel[\s\S]+Verify Gallery first-Generate parity[\s\S]+npm run test:web:gallery-publication/
  );
  assert.match(webPrSmoke, /if: failure\(\)[\s\S]+path: test-results\//);
  assert.doesNotMatch(
    webPrSmoke,
    /functional-full|perf-smoke|losat-cache|Vibrio/
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
    'web-contracts-pr',
    'web-pr-smoke',
    'playwright-functional'
  ]);
  assert.match(
    gate,
    /if: >-\n      always\(\) &&\n      github\.event_name == 'pull_request' &&\n      github\.base_ref == 'dev'/
  );
  assert.match(
    gate,
    /ref: \$\{\{ github\.event\.pull_request\.base\.sha \}\}[\s\S]*path: \.ci-trusted-base[\s\S]*persist-credentials: false/
  );
  // Node comes from the base branch tip: github.sha is the PR merge commit and base.sha can predate .node-version.
  assert.equal([...gate.matchAll(/uses: actions\/checkout@/g)].length, 2);
  assert.match(
    gate,
    /ref: refs\/heads\/\$\{\{ github\.base_ref \}\}\n {10}path: \.ci-node-version\n {10}sparse-checkout: \.node-version\n {10}sparse-checkout-cone-mode: false\n {10}persist-credentials: false\n[\s\S]*node-version-file: \.ci-node-version\/\.node-version\n[\s\S]*Validate required PR evidence/
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
    'web-contracts-pr',
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
  // The trusted plan selects functional Playwright on pull requests; the job is shared with dev.
  assert.match(
    workflowJob('playwright-functional'),
    /\(\(github\.event_name == 'pull_request' && github\.base_ref == 'dev'\) \|\|/
  );
  for (const jobId of [
    'core',
    'browser',
    'playwright-performance',
    'acceptance-supported-main',
    'slow-main',
    'vibrio-generate-release',
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
    'playwright test --config=tests/web/playwright/functional.config.js'
  );
  assert.equal(
    PACKAGE_SCRIPTS['test:web:pr-smoke'],
    'playwright test --config=tests/web/playwright/pr-smoke.config.js'
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

  // Expanded Playwright discovery (including nested/parameterized tests) is checked
  // by tests/ci/playwright-inventory.test.mjs, not a tag-source regex.
});

test('exact dev staging routes every job through the protected-branch plan', () => {
  const stagingJobs = [
    'web-change-budget',
    'core',
    'recipes-standard',
    'gallery',
    'browser',
    'playwright-functional',
    'playwright-performance',
    'lint',
    'losat-cache-browser-acceptance'
  ];
  for (const jobId of stagingJobs) {
    const job = workflowJob(jobId);
    assert.match(job, /needs: ci-impact/, jobId);
    assert.match(job, /needs\.ci-impact\.result == 'success'/, jobId);
    assert.match(job, new RegExp(`requiredJobs, '${jobId}'`), jobId);
    assert.match(job, /github\.event_name == 'push' && github\.ref == 'refs\/heads\/dev'/);
    assert.match(
      job,
      /github\.event_name == 'workflow_dispatch' && github\.ref == 'refs\/heads\/dev'/
    );
  }

  const fullPlaywright = workflowJob('playwright-functional');
  assert.match(fullPlaywright, /\n    name: Playwright functional \(shard \$\{\{ matrix\.shard \}\}\/8\)\n/);
  assert.match(fullPlaywright, /shard: \[1, 2, 3, 4, 5, 6, 7, 8\]/);
  assert.match(
    fullPlaywright,
    /npm run test:web:functional-full -- --reporter=line,github,json \$files\n/
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
    'vibrio-generate-release',
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
  assert.match(gate, /ref: \$\{\{ github\.sha \}\}/);
  assert.match(gate, /CI_IMPACT_PLAN_JSON: \$\{\{ needs\.ci-impact\.outputs\.plan \}\}/);
  assert.match(gate, /CI_IMPACT_NEEDS_JSON: \$\{\{ toJSON\(needs\) \}\}/);
  assert.match(gate, /CI_IMPACT_EXPECTED_PROFILE: dev/);
  assert.match(gate, /CI_IMPACT_EXPECTED_WORKFLOW_SHA: \$\{\{ github\.sha \}\}/);
  assert.match(gate, /node tools\/ci-impact\.mjs gate/);
  assert.doesNotMatch(gate, /test "\$\{\{ needs\.|gh api|curl|check-promotion-readiness/);
});

test('Gallery readiness routes jobs from direct-parent evidence and aggregates the current SHA', () => {
  assert.match(
    GALLERY_PUBLICATION_WORKFLOW,
    /permissions:\n  contents: read\n  actions: read/
  );
  assert.doesNotMatch(GALLERY_PUBLICATION_WORKFLOW, /\n    paths(?:-ignore)?:/);

  const planner = workflowJob('ci-impact', GALLERY_PUBLICATION_WORKFLOW);
  assert.match(planner, /\n    name: CI impact plan\n/);
  assert.match(planner, /timeout-minutes: 5/);
  assert.match(planner, /plan: \$\{\{ steps\.plan\.outputs\.plan \}\}/);
  assert.match(
    planner,
    /Checkout exact workflow SHA[\s\S]*ref: \$\{\{ github\.sha \}\}[\s\S]*fetch-depth: 0[\s\S]*persist-credentials: false/
  );
  assert.match(planner, /node-version-file: \.node-version/);
  assert.match(planner, /CI_IMPACT_PROFILE: gallery/);
  assert.match(planner, /CI_IMPACT_EVENT_NAME: \$\{\{ github\.event_name \}\}/);
  assert.match(
    planner,
    /CI_IMPACT_CHANGE_BASE_SHA: \$\{\{ github\.event_name == 'push' && github\.event\.before \|\| github\.sha \}\}/
  );
  assert.match(planner, /CI_IMPACT_CHANGE_HEAD_SHA: \$\{\{ github\.sha \}\}/);
  assert.match(planner, /CI_IMPACT_WORKFLOW_SHA: \$\{\{ github\.sha \}\}/);
  assert.match(planner, /GITHUB_TOKEN: \$\{\{ github\.token \}\}/);
  assert.match(planner, /node tools\/ci-impact\.mjs plan/);
  assert.doesNotMatch(planner, /node --test/);

  const browser = workflowJob('browser', GALLERY_PUBLICATION_WORKFLOW);
  assert.match(browser, /\n    name: Gallery browser \(common 9\)\n/);
  assert.match(browser, /needs: ci-impact/);
  assert.match(browser, /needs\.ci-impact\.result == 'success'/);
  assert.match(browser, /requiredJobs, 'browser'/);
  assert.match(browser, /ref: \$\{\{ github\.sha \}\}/);
  assert.match(browser, /npm run test:web:gallery-publication/);
  assert.equal(
    [...TEST_WORKFLOW.matchAll(/npm run test:web:gallery-publication/g)].length,
    1
  );
  assert.doesNotMatch(browser, /matrix|Vibrio|test:web:vibrio-generate/);
  assert.equal(
    PACKAGE_SCRIPTS['test:web:vibrio-generate'],
    'playwright test --config=tests/web/playwright/vibrio.config.js'
  );
  assert.match(DEPLOY_WORKFLOW, /npm run test:web:vibrio-generate/);

  const performance = workflowJob('performance', GALLERY_PUBLICATION_WORKFLOW);
  assert.match(performance, /\n    name: Gallery publication performance \(projection\)\n/);
  assert.match(performance, /needs: ci-impact/);
  assert.match(performance, /needs\.ci-impact\.result == 'success'/);
  assert.match(performance, /requiredJobs, 'performance'/);
  assert.match(performance, /measure_gallery_publication_performance\.py projection/);
  assert.match(performance, /github\.event_name == 'workflow_dispatch' && inputs\.complete_refresh/);
  assert.match(performance, /measure_gallery_publication_performance\.py refresh/);
  assert.equal(
    [...performance.matchAll(/--baseline-revision 574b33b83962949397839e2aaa862a8b96667625/g)].length,
    2
  );

  const gate = workflowJob('readiness-gate', GALLERY_PUBLICATION_WORKFLOW);
  assert.match(gate, /\n    name: Gallery readiness \/ gate\n/);
  assert.match(gate, /\n    needs:\n      - ci-impact\n      - browser\n      - performance\n/);
  assert.match(gate, /\n    if: always\(\)\n/);
  assert.match(gate, /ref: \$\{\{ github\.sha \}\}/);
  assert.match(gate, /CI_IMPACT_PLAN_JSON: \$\{\{ needs\.ci-impact\.outputs\.plan \}\}/);
  assert.match(gate, /CI_IMPACT_NEEDS_JSON: \$\{\{ toJSON\(needs\) \}\}/);
  assert.match(gate, /CI_IMPACT_EXPECTED_PROFILE: gallery/);
  assert.match(gate, /CI_IMPACT_EXPECTED_WORKFLOW_SHA: \$\{\{ github\.sha \}\}/);
  assert.match(gate, /node tools\/ci-impact\.mjs gate/);
  assert.doesNotMatch(
    gate,
    /test "\$\{\{ needs\.|gh api|curl|check-promotion-readiness/
  );
});

test('trusted promotion uses base code for topology, exact-SHA evidence, and tree proof', () => {
  const devPolicy = workflowJob('trusted-base-policy', BASE_POLICY_WORKFLOW);
  assert.match(devPolicy, /\n    name: Web base policy \(trusted base\)\n/);
  assert.match(devPolicy, /if: github\.event\.pull_request\.base\.ref == 'dev'/);
  // The trusted base checkout plus a sparse checkout of .node-version from the base branch tip.
  assert.equal([...devPolicy.matchAll(/uses: actions\/checkout@/g)].length, 2);
  assert.match(
    devPolicy,
    /ref: \$\{\{ github\.sha \}\}\n {10}path: \.ci-node-version\n {10}sparse-checkout: \.node-version\n {10}sparse-checkout-cone-mode: false\n {10}persist-credentials: false\n/
  );
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
  // The trusted base checkout plus a sparse checkout of .node-version from the base branch tip.
  assert.equal([...promotion.matchAll(/uses: actions\/checkout@/g)].length, 2);
  assert.match(
    promotion,
    /ref: \$\{\{ github\.sha \}\}\n {10}path: \.ci-node-version\n {10}sparse-checkout: \.node-version\n {10}sparse-checkout-cone-mode: false\n {10}persist-credentials: false\n/
  );
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
    'ci-inventory.yml',
    'codeql.yml',
    'deploy_web.yml',
    'gallery-publication.yml',
    'losat-distribution.yml',
    'release.yml',
    'test.yml',
    'translation-review.yml',
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

test('every workflow reads its Node version from the root .node-version', () => {
  assert.match(readFileSync(join(REPOSITORY_ROOT, '.node-version'), 'utf8'), /^\d+\n$/);
  WORKFLOW_NAMES.forEach((name, index) => {
    const source = WORKFLOW_SOURCES[index];
    assert.doesNotMatch(source, /^\s*node-version:/m, name);
    assert.equal(
      (source.match(/uses: actions\/setup-node@(?:v\d+|[0-9a-f]{40} # v\d+\.\d+\.\d+)\n\s+with:\n\s+node-version-file: (?:\.ci-node-version\/)?\.node-version\n/g) || []).length,
      (source.match(/uses: actions\/setup-node@/g) || []).length,
      name
    );
    // pull_request.base.sha can predate the .node-version file, so no job reads it from that checkout.
    for (const job of source.split(/\n(?=  [a-z0-9-]+:\n)/)) {
      for (const checkout of job.match(/uses: actions\/checkout@[^\n]*\n {8}with:\n(?: {10}[^\n]*\n)*/g) || []) {
        if (!/ref: \$\{\{ github\.event\.pull_request\.base\.sha \}\}/.test(checkout)) continue;
        const path = checkout.match(/\n {10}path: (\S+)/)?.[1];
        const file = path ? `${path}/.node-version` : '.node-version';
        assert.ok(!job.includes(`node-version-file: ${file}\n`), `${name}: ${file} comes from base.sha`);
      }
    }
  });
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
    'Maintainer approval (review or comment link)'
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
  const job = TEST_WORKFLOW.match(/\n  web-change-budget:\n[\s\S]*?(?=\n  [a-z0-9-]+:\n|$)/)?.[0];
  assert.ok(job, 'Web change-budget job must exist');
  assert.match(job, /\n    timeout-minutes: 10\n/);

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

test('supported-version and slow matrices are mandatory in explicit release acceptance', () => {
  for (const jobName of ['acceptance-supported-main', 'slow-main', 'vibrio-generate-release']) {
    const job = TEST_WORKFLOW.match(
      new RegExp(`\\n  ${jobName}:\\n[\\s\\S]*?(?=\\n  [a-z0-9-]+:\\n|$)`)
    )?.[0];
    assert.ok(job, `${jobName} job must exist`);
    assert.match(job, /needs: ci-impact/);
    assert.match(job, /needs\.ci-impact\.result == 'success'/);
    assert.match(job, new RegExp(`requiredJobs, '${jobName}'`));
    assert.match(job, /github\.event_name == 'push' && github\.ref == 'refs\/heads\/dev'/);
    assert.match(
      job,
      /github\.event_name == 'workflow_dispatch' && github\.ref == 'refs\/heads\/dev'/
    );
  }
  const gate = workflowJob('release-gate');
  assert.match(gate, /name: Release \/ gate/);
  assert.match(gate, /inputs\.tier == 'release'/);
  assert.match(gate, /CI_IMPACT_EXPECTED_PROFILE: release/);
  assert.match(gate, /Require exact candidate Gallery readiness/);
  assert.match(workflowJob('dev-staging-gate'), /inputs\.tier != 'release'/);
});

test('release acceptance runs the Vibrio full-generation suite that main deploy also keeps', () => {
  const job = workflowJob('vibrio-generate-release');
  assert.match(job, /\n    name: Vibrio full generation\n/);
  assert.match(job, /\n {8}run: npm run test:web:vibrio-generate\n/);
  assert.match(job, /python tools\/prepare_browser_wheel\.py/);
  // The release profile requires it and the dev profile does not.
  assert.deepEqual(knownJobsFor('release').filter((id) => !knownJobsFor('dev').includes(id)), [
    'acceptance-supported-main',
    'slow-main',
    'vibrio-generate-release'
  ]);
  // Both aggregate gates list it, so a skipped or failed run blocks release evidence.
  for (const gateId of ['release-gate', 'dev-staging-gate']) {
    assert.match(workflowJob(gateId), /\n {6}- vibrio-generate-release\n/, gateId);
  }
  // The hosted-Web verification on main keeps its own run.
  assert.match(DEPLOY_WORKFLOW, /npm run test:web:vibrio-generate/);
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
const WEB_OWNER_GRAPH_DETECTOR_MODULE = join(
  REPOSITORY_ROOT,
  'tools/web-owner-graph-detectors.mjs'
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
    'services/svg-result-ingestion.js': (
      WEB_CHANGE_POLICY.allowedPrivilegedImporters['services/svg-result-ingestion.js']
    ),
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
const PRODUCT_CONTRACT_AUTHORITY_PATH = 'docs/internal/OPTION_INTEGRITY_PRODUCT_CONTRACT.md';
const CURRENT_RESULT_ADMISSION_FIXTURE_SOURCES = Object.freeze({
  'gbdraw/web/js/services/svg-result-ingestion.js': (
    '// Fixture-only detector input; independent of architecture authority.\n'
      + 'export const admitCurrentGeneratedResults = () => ({ admitted: true });\n'
  ),
  'gbdraw/web/js/app/candidate-render.js': (
    '// Fixture-only detector input; independent of architecture authority.\n'
      + "import { admitCurrentGeneratedResults } from '../services/svg-result-ingestion.js';\n"
      + 'export const admitCandidateResult = () => admitCurrentGeneratedResults();\n'
  )
});

const normalizeMappedEvidenceFixturePath = (path) => {
  const portablePath = String(path || '').replaceAll('\\', '/');
  if (!portablePath
      || isAbsolute(portablePath)
      || /^[A-Za-z]:\//.test(portablePath)) {
    throw new Error(`Mapped evidence fixture path must be repository-relative: ${path}`);
  }
  const segments = portablePath.split('/');
  if (segments.includes('..')) {
    throw new Error(`Mapped evidence fixture path must not traverse parents: ${path}`);
  }
  const normalized = segments.filter((segment) => segment && segment !== '.').join('/');
  if (!normalized) {
    throw new Error(`Mapped evidence fixture path must name a file: ${path}`);
  }
  return normalized;
};

const mappedEvidenceFixtureFiles = (productImpactMapSource) => {
  const productImpactMap = JSON.parse(productImpactMapSource);
  const paths = new Set(
    productImpactMap.concerns.flatMap(({ contracts }) => contracts.map(({ ref }) => (
      normalizeMappedEvidenceFixturePath(ref.split('::', 1)[0])
    )))
  );
  return Object.fromEntries([...paths].sort().map((path) => [
    path,
    '// Generic mapped-evidence dependency for a synthetic repository.\n'
  ]));
};

const syntheticArchitectureFixtureDependencies = (productImpactMapSource) => Object.freeze({
  ...mappedEvidenceFixtureFiles(productImpactMapSource),
  ...CURRENT_RESULT_ADMISSION_FIXTURE_SOURCES
});
const SYNTHETIC_ARCHITECTURE_FIXTURE_DEPENDENCIES = (
  syntheticArchitectureFixtureDependencies(PRODUCT_IMPACT_MAP_SOURCE)
);
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
  ...SYNTHETIC_ARCHITECTURE_FIXTURE_DEPENDENCIES,
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
  'tools/web-owner-graph-detectors.mjs': readFileSync(
    WEB_OWNER_GRAPH_DETECTOR_MODULE,
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
const WEB_DESIGN_RULES_PATH = 'gbdraw/web/CLAUDE.md';
const DESIGN_RULE_GUARD_REGISTRY_PATH = 'tools/web-design-rule-guards.json';
const OWNER_GRAPH_REGISTRY_PATH = 'tools/web-owner-graph.json';
const FUTURE_GUARD_PATHS = Object.freeze([
  'docs/internal/ARCHITECTURE_FITNESS_FUNCTION_RATCHET.md',
  '.github/pull_request_template.md',
  'tools/web-architecture-violations.json',
  WEB_DESIGN_RULES_PATH,
  DESIGN_RULE_GUARD_REGISTRY_PATH,
  OWNER_GRAPH_REGISTRY_PATH
]);
const PROTECTED_ARCHITECTURE_GUARD_PATHS = Object.freeze([
  ...FUTURE_GUARD_PATHS,
  'tools/web-architecture-rules.json'
]);
const FUTURE_AUTHORITY_PATHS = Object.freeze([
  'docs/internal/ARCHITECTURE_FITNESS_FUNCTION_RATCHET.md',
  'tools/web-architecture-rules.json',
  'tools/web-architecture-violations.json',
  WEB_DESIGN_RULES_PATH,
  DESIGN_RULE_GUARD_REGISTRY_PATH,
  OWNER_GRAPH_REGISTRY_PATH
]);
const PRODUCT_IMPACT_GUARD_PATHS = Object.freeze([
  'docs/internal/PRODUCT_IMPACT_RATCHET.md',
  PRODUCT_CONTRACT_AUTHORITY_PATH,
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
  PRODUCT_CONTRACT_AUTHORITY_PATH,
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
  ...SYNTHETIC_ARCHITECTURE_FIXTURE_DEPENDENCIES,
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
  'tools/web-owner-graph-detectors.mjs': readFileSync(
    WEB_OWNER_GRAPH_DETECTOR_MODULE,
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
  const replacementsByRule = new Map();
  replacementsByRule.set(
    'canonical-path.render-request',
    rules.get('canonical-path.render-request').allowedEdges.map(({ from, to }) => ({
      ruleKey: 'canonical-path.render-request',
      subjectCategory: 'canonical-entry-edge',
      subject: `${from} -> ${to}`
    }))
  );
  replacementsByRule.set(
    'semantic-owner.render-request',
    rules.get('semantic-owner.render-request').allowedDefinitionPaths.map((subject) => ({
      ruleKey: 'semantic-owner.render-request',
      subjectCategory: 'definition-path',
      subject
    }))
  );
  map.concerns.flatMap(({ options }) => options)
    .flatMap(({ requirements }) => requirements)
    .forEach((requirement) => {
      const replacedRules = new Set();
      requirement.anyOf = requirement.anyOf.flatMap((candidate) => {
        const replacements = replacementsByRule.get(candidate.ruleKey);
        if (!replacements) return [candidate];
        if (replacedRules.has(candidate.ruleKey)) return [];
        replacedRules.add(candidate.ruleKey);
        return replacements;
      });
    });
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

test('synthetic architecture fixtures grow through evidence dependencies, not alternate authority', () => {
  const assemblySource = [
    normalizeMappedEvidenceFixturePath,
    mappedEvidenceFixtureFiles,
    syntheticArchitectureFixtureDependencies,
    ...Object.values(CURRENT_RESULT_ADMISSION_FIXTURE_SOURCES)
  ].join('\n');
  for (const authorityField of [
    'allowedEdges',
    'allowedDefinitionPaths',
    'exactActiveEdgeCount',
    'exactDefinitionCount',
    'baselineEligible'
  ]) {
    assert.equal(assemblySource.includes(authorityField), false, authorityField);
  }

  const mappedPaths = Object.keys(mappedEvidenceFixtureFiles(PRODUCT_IMPACT_MAP_SOURCE));
  for (const fixture of [BUDGET_FIXTURE, TRUSTED_ARCHITECTURE_FIXTURE_FILES]) {
    mappedPaths.forEach((path) => assert.equal(Object.hasOwn(fixture, path), true, path));
  }
  assert.notEqual(
    BUDGET_FIXTURE['tests/web/architecture-contracts.test.mjs'],
    SYNTHETIC_ARCHITECTURE_FIXTURE_DEPENDENCIES[
      'tests/web/architecture-contracts.test.mjs'
    ]
  );
  assert.notEqual(
    TRUSTED_ARCHITECTURE_FIXTURE_FILES['tests/web/session-request.test.mjs'],
    SYNTHETIC_ARCHITECTURE_FIXTURE_DEPENDENCIES['tests/web/session-request.test.mjs']
  );

  const syntheticSources = new Map(Object.entries(
    SYNTHETIC_ARCHITECTURE_FIXTURE_DEPENDENCIES
  ).filter(([path]) => path.startsWith('gbdraw/web/js/'))
    .map(([path, source]) => [path.replace('gbdraw/web/js/', ''), source]));
  assert.deepEqual(
    WEB_ARCHITECTURE_DETECTORS['semantic-owner.current-result-admission.v1']
      .detect(syntheticSources).subjects,
    ['services/svg-result-ingestion.js']
  );
  assert.deepEqual(
    WEB_ARCHITECTURE_DETECTORS['canonical-path.current-result-admission.v1']
      .detect(syntheticSources).subjects,
    ['app/candidate-render.js -> services/svg-result-ingestion.js']
  );

  const duplicateMap = JSON.stringify({
    concerns: [{
      contracts: [
        { ref: './tests/web/evidence.test.mjs::first' },
        { ref: 'tests/web/evidence.test.mjs::second' }
      ]
    }]
  });
  assert.deepEqual(Object.keys(mappedEvidenceFixtureFiles(duplicateMap)), [
    'tests/web/evidence.test.mjs'
  ]);
  for (const ref of [
    '/tmp/outside.test.mjs::case',
    'C:/outside.test.mjs::case',
    '../outside.test.mjs::case',
    'tests/../../outside.test.mjs::case'
  ]) {
    assert.throws(
      () => mappedEvidenceFixtureFiles(JSON.stringify({
        concerns: [{ contracts: [{ ref }] }]
      })),
      /Mapped evidence fixture path/
    );
  }
});

test('otherwise unmodified synthetic repositories pass current inherited authority', () => {
  const budget = runChangeBudgetCase(() => {});
  assert.equal(budget.status, 0, budget.output);

  withTrustedArchitectureRepository(
    (write) => write('fixture-note.txt', 'current inherited authority\n'),
    ({ status, output }) => {
      assert.equal(status, 0, output);
      WEB_ARCHITECTURE_RULES.rules.forEach(({ key }) => {
        assert.match(
          output,
          new RegExp(`${key.replaceAll('.', '\\.')} .* observation=CONFORMING .* decision=PASS`)
        );
      });
    }
  );
});

test('active base rules accept one canonical entry edge', () => {
  const registeredOwnerCount = WEB_ARCHITECTURE_RULES.rules.filter(
    ({ kind }) => kind === 'single-semantic-owner'
  ).length;
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
        /First-party static import graph[\s\S]*Before: 5 modules, 3 edges, 0 cycles[\s\S]*After: 5 modules, 3 edges, 0 cycles/
      );
      assert.match(
        output,
        new RegExp(
          `Registered Authority Location Count \\| ${registeredOwnerCount} \\| 0 \\| 0 \\| `
            + `${registeredOwnerCount} \\| 0 \\| registered rule`
        )
      );
    }
  );
});

test('an extra registered authority location fails and remains visible in the delta', () => {
  const registeredOwnerCount = WEB_ARCHITECTURE_RULES.rules.filter(
    ({ kind }) => kind === 'single-semantic-owner'
  ).length;
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
        new RegExp(
          `Registered Authority Location Count \\| ${registeredOwnerCount} \\| 1 \\| 0 \\| `
            + `${registeredOwnerCount + 1} \\| \\+1 \\| registered rule`
        )
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
      assert.match(output, /After: 8 modules, 5 edges, 0 cycles/);
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
  const registeredOwnerCount = WEB_ARCHITECTURE_RULES.rules.filter(
    ({ kind }) => kind === 'single-semantic-owner'
  ).length;
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
        new RegExp(
          `Registered Authority Location Count \\| ${registeredOwnerCount} \\| 0 \\| 0 \\| `
            + `${registeredOwnerCount} \\| 0 \\| registered rule`
        )
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
      assert.match(output, /Candidate authority validation: VALID \(inert map and decision data only\)/);
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

const VENDOR_LIBRARY_PATH = 'gbdraw/web/vendor/library.js';
const VENDOR_LIBRARY_UPDATE = 'globalThis.vendorLibrary = false;\n';
const sha256Hex = (content) => createHash('sha256').update(content).digest('hex');
const vendorChangeEntry = (overrides = {}) => ({
  path: VENDOR_LIBRARY_PATH,
  change: 'replace',
  sha256: sha256Hex(VENDOR_LIBRARY_UPDATE),
  package: 'library',
  version: '1.0.1',
  source: 'https://example.invalid/library@1.0.1/library.js',
  ...overrides
});
const writeVendorChangeListing = (write, entries) => {
  writeBudgetPolicy(write, { ...cloneBudgetPolicy(), allowedVendorChanges: entries });
};
// Merges the listing as an earlier authority-only pull request (and, when given, a later
// merged change), then evaluates the next pull request against that trusted base.
const runListedVendorChangeCase = (entries, mutate, merged = null) => withChangeBudgetRepository(
  ({ commit, execute, git, root, write }) => {
    const writeOrRemove = (path, content) => (
      content === null ? rmSync(join(root, path)) : write(path, content)
    );
    if (entries) writeVendorChangeListing(write, entries);
    if (merged) merged(writeOrRemove);
    const base = entries
      ? commit('authority-only vendored runtime listing')
      : git(['rev-parse', 'HEAD']).stdout.trim();
    mutate(writeOrRemove);
    const head = commit('vendored runtime change');
    return execute({ base, head });
  }
);
const assertVendorGateFailure = (result, pattern, context) => {
  assert.equal(result.status, 1, `${context}\n${result.output}`);
  assert.match(result.output, /Gate: \*\*FAIL\*\*/, context);
  assert.match(result.output, pattern, context);
};
const UNLISTED_VENDOR_CHANGE = /changes under gbdraw\/web\/vendor\/ are not allowed unless the trusted base lists them/;

test('a vendored runtime change passes only when the trusted base lists its exact bytes', () => {
  const unlisted = runListedVendorChangeCase(
    [vendorChangeEntry({ path: 'gbdraw/web/vendor/other.js', change: 'add' })],
    (write) => write(VENDOR_LIBRARY_PATH, VENDOR_LIBRARY_UPDATE)
  );
  assertVendorGateFailure(unlisted, UNLISTED_VENDOR_CHANGE, 'unlisted path');
  assert.match(unlisted.output, /gbdraw\/web\/vendor\/library\.js: not listed/);

  const listedCases = [
    {
      name: 'replace',
      entry: vendorChangeEntry(),
      mutate: (write) => write(VENDOR_LIBRARY_PATH, VENDOR_LIBRARY_UPDATE)
    },
    {
      name: 'add',
      entry: vendorChangeEntry({
        path: 'gbdraw/web/vendor/library.prod.js',
        change: 'add'
      }),
      mutate: (write) => write('gbdraw/web/vendor/library.prod.js', VENDOR_LIBRARY_UPDATE)
    },
    {
      name: 'remove',
      entry: vendorChangeEntry({
        change: 'remove',
        sha256: sha256Hex(BUDGET_FIXTURE[VENDOR_LIBRARY_PATH])
      }),
      mutate: (write) => write(VENDOR_LIBRARY_PATH, null)
    }
  ];
  listedCases.forEach(({ name, entry, mutate }) => {
    const result = runListedVendorChangeCase([entry], mutate);
    assert.equal(result.status, 0, `${name}\n${result.output}`);
    assert.match(result.output, /Gate: \*\*PASS\*\*/, name);
    assert.match(result.output, /Review: \*\*REQUIRED\*\*/, name);
    assert.ok(
      result.output.includes(`Vendored runtime: listed ${name} of ${entry.path} (library 1.0.1)`),
      `${name}\n${result.output}`
    );
  });
});

test('a listed vendored runtime change fails on a hash or kind mismatch', () => {
  const mismatch = runListedVendorChangeCase(
    [vendorChangeEntry({ sha256: sha256Hex('globalThis.vendorLibrary = 0;\n') })],
    (write) => write(VENDOR_LIBRARY_PATH, VENDOR_LIBRARY_UPDATE)
  );
  assertVendorGateFailure(mismatch, UNLISTED_VENDOR_CHANGE, 'hash mismatch');
  assert.match(
    mismatch.output,
    new RegExp(`library\\.js: sha256 ${sha256Hex(VENDOR_LIBRARY_UPDATE)} does not match the listed`)
  );

  const kind = runListedVendorChangeCase(
    [vendorChangeEntry()],
    (write) => write(VENDOR_LIBRARY_PATH, null)
  );
  assertVendorGateFailure(kind, /library\.js: listed as replace, but the change is remove/, 'kind');
});

const IMMUTABLE_ASSETS_MANIFEST_PATH = 'tools/immutable-web-assets.json';
const immutableManifestIn = (files) => (write) => write(
  IMMUTABLE_ASSETS_MANIFEST_PATH,
  typeof files === 'string' ? files : `${JSON.stringify({ prefixes: ['vendor/'], files }, null, 2)}\n`
);

test('a listed replace of a file the base pins as immutable fails, other listed changes pass', () => {
  const pinned = immutableManifestIn({ 'vendor/library.js': sha256Hex(BUDGET_FIXTURE[VENDOR_LIBRARY_PATH]) });
  const replaced = runListedVendorChangeCase(
    [vendorChangeEntry()],
    (write) => write(VENDOR_LIBRARY_PATH, VENDOR_LIBRARY_UPDATE),
    pinned
  );
  assertVendorGateFailure(replaced, UNLISTED_VENDOR_CHANGE, 'pinned replace');
  assert.match(
    replaced.output,
    /library\.js: pinned in tools\/immutable-web-assets\.json and served as immutable; publish the new bytes under a new path/
  );

  const unpinned = runListedVendorChangeCase(
    [vendorChangeEntry()],
    (write) => write(VENDOR_LIBRARY_PATH, VENDOR_LIBRARY_UPDATE),
    immutableManifestIn({ 'vendor/other.js': sha256Hex('other') })
  );
  assert.equal(unpinned.status, 0, unpinned.output);
  assert.match(unpinned.output, /Gate: \*\*PASS\*\*/);

  const added = runListedVendorChangeCase(
    [vendorChangeEntry({ path: 'gbdraw/web/vendor/library.v2.js', change: 'add' })],
    (write) => write('gbdraw/web/vendor/library.v2.js', VENDOR_LIBRARY_UPDATE),
    pinned
  );
  assert.equal(added.status, 0, added.output);
  assert.match(added.output, /Gate: \*\*PASS\*\*/);
});

test('a malformed immutable-asset manifest in the base fails a vendored replace closed', () => {
  ['{not json', '{"files": []}'].forEach((manifest) => {
    const result = runListedVendorChangeCase(
      [vendorChangeEntry()],
      (write) => write(VENDOR_LIBRARY_PATH, VENDOR_LIBRARY_UPDATE),
      immutableManifestIn(manifest)
    );
    assertVendorGateFailure(
      result,
      /tools\/immutable-web-assets\.json in the base is malformed/,
      manifest
    );
  });
});

test('an added binary vendored file passes only when the trusted base lists it as add with its bytes', () => {
  const fontPath = 'gbdraw/web/vendor/fonts/latin.woff2';
  const fontBytes = Buffer.from([0x77, 0x4f, 0x46, 0x32, 0, 1, 2, 3]);
  const entry = vendorChangeEntry({ path: fontPath, change: 'add', sha256: sha256Hex(fontBytes) });
  const addedBinarySection = (output) => output.split('## Added binary runtime files')[1].split('\n## ')[0];

  const listed = runListedVendorChangeCase([entry], (write) => write(fontPath, fontBytes));
  assert.equal(listed.status, 0, `listed\n${listed.output}`);
  assert.match(listed.output, /Gate: \*\*PASS\*\*/);
  assert.match(listed.output, /Review: \*\*REQUIRED\*\*/);
  assert.doesNotMatch(listed.output, /added binary runtime files are not allowed/);
  assert.ok(
    listed.output.includes(`Vendored runtime: listed add of ${fontPath} (library 1.0.1)`),
    listed.output
  );
  assert.equal(addedBinarySection(listed.output).trim(), '- None', listed.output);

  const otherBytes = runListedVendorChangeCase([entry], (write) => write(fontPath, Buffer.from([0, 9, 9, 9])));
  assertVendorGateFailure(otherBytes, /added binary runtime files are not allowed/, 'other bytes');
  assert.match(otherBytes.output, UNLISTED_VENDOR_CHANGE);
  assert.match(otherBytes.output, /latin\.woff2: sha256 [0-9a-f]{64} does not match the listed/);
  assert.ok(addedBinarySection(otherBytes.output).includes(`- ${fontPath}`), otherBytes.output);

  const unlisted = runListedVendorChangeCase(
    [vendorChangeEntry({ path: 'gbdraw/web/vendor/other.js', change: 'add' })],
    (write) => write(fontPath, fontBytes)
  );
  assertVendorGateFailure(unlisted, /added binary runtime files are not allowed/, 'unlisted');
  assert.match(unlisted.output, /latin\.woff2: not listed/);
});

test('a vendored runtime listing in the pull request head authorizes nothing', () => {
  // The candidate lists its own change: the trusted base has no entry, so the vendored
  // file stays unlisted, and the listing beside runtime is a guard co-change as well.
  const result = runListedVendorChangeCase(null, (write) => {
    writeVendorChangeListing(write, [vendorChangeEntry()]);
    write(VENDOR_LIBRARY_PATH, VENDOR_LIBRARY_UPDATE);
  });
  assertVendorGateFailure(result, UNLISTED_VENDOR_CHANGE, 'head-only listing');
  assert.match(result.output, /gbdraw\/web\/vendor\/library\.js: not listed/);
  assert.match(result.output, /production runtime files and Web guard\/CI files changed together/);
});

test('a vendored runtime change cannot edit the listing in the same pull request', () => {
  // The trusted base lists the change, but the pull request also edits the listing.
  const consumed = runListedVendorChangeCase([vendorChangeEntry()], (write) => {
    writeVendorChangeListing(write, []);
    write(VENDOR_LIBRARY_PATH, VENDOR_LIBRARY_UPDATE);
  });
  assertVendorGateFailure(
    consumed,
    /production runtime files and Web guard\/CI files changed together/,
    'listing removed beside the vendored file'
  );

  // A safe privileged contraction may change the policy beside runtime, but not the listing.
  const otherEntry = vendorChangeEntry({ path: 'gbdraw/web/vendor/other.js', change: 'add' });
  const contraction = runListedVendorChangeCase([otherEntry], (write) => {
    const policy = cloneBudgetPolicy();
    removePolicyPath(policy, 'allowedPrivilegedOwners', 'Diagram Worker', 'app/editor.js');
    writeBudgetPolicy(write, {
      ...policy,
      allowedVendorChanges: [otherEntry, vendorChangeEntry()]
    });
    removeEditorOwnerUse(write);
  });
  assertVendorGateFailure(
    contraction,
    /production runtime files and Web guard\/CI files changed together/,
    'listing added beside a safe contraction'
  );
});

test('a malformed vendored runtime listing fails closed in the base and the candidate', () => {
  const malformed = [
    [vendorChangeEntry({ sha256: 'ABC' }), /sha256 must be 64 lowercase hexadecimal digits/],
    [vendorChangeEntry({ path: 'gbdraw/web/js/app/editor.js' }), /path must be a normalized path under gbdraw\/web\/vendor\//],
    [vendorChangeEntry({ path: 'gbdraw/web/vendor/../js/app/editor.js' }), /path must be a normalized path under gbdraw\/web\/vendor\//],
    [vendorChangeEntry({ change: 'move' }), /change must be one of add, replace, remove/],
    [vendorChangeEntry({ source: 'http://example.invalid/library.js' }), /source must be an https URL/],
    [vendorChangeEntry({ version: '' }), /version must be a non-empty string/],
    [{ ...vendorChangeEntry(), note: 'extra' }, /must have exactly the keys change, package, path, sha256, source, version/]
  ];
  malformed.forEach(([entry, pattern]) => {
    const candidate = runChangeBudgetCase((write) => writeVendorChangeListing(write, [entry]));
    assertVendorGateFailure(candidate, /candidate vendored runtime listing: /, `candidate ${pattern}`);
    assert.match(candidate.output, pattern);
  });

  const duplicate = runChangeBudgetCase((write) => writeVendorChangeListing(
    write,
    [vendorChangeEntry(), vendorChangeEntry({ change: 'remove' })]
  ));
  assertVendorGateFailure(duplicate, /lists gbdraw\/web\/vendor\/library\.js twice/, 'duplicate');

  const trusted = runListedVendorChangeCase(
    [vendorChangeEntry({ sha256: 'ABC' })],
    (write) => write('gbdraw/web/js/app/secondary.js', 'export const editSecondaryOwner = () => 2;\n')
  );
  assertVendorGateFailure(trusted, /trusted vendored runtime listing: /, 'trusted base');
});

const describedVendorEntry = (entry) => (
  `- ${entry.path}: ${entry.change} sha256=${entry.sha256} package=${entry.package} `
  + `version=${entry.version} source=${entry.source}`
);

test('an authority-only vendored runtime listing passes with review and reports each entry', () => {
  const entry = vendorChangeEntry();
  const added = runChangeBudgetCase((write) => writeVendorChangeListing(write, [entry]));
  assert.equal(added.status, 0, added.output);
  assert.match(added.output, /Review: \*\*REQUIRED\*\*/);
  assert.match(added.output, /vendored runtime listing changed \(added=1, removed=0;/);
  assert.ok(
    reportSection(added.output, 'Added vendored runtime listing entries').includes(describedVendorEntry(entry)),
    added.output
  );

  const removed = runListedVendorChangeCase([entry], (write) => writeVendorChangeListing(write, []));
  assert.equal(removed.status, 0, removed.output);
  assert.match(removed.output, /vendored runtime listing changed \(added=0, removed=1;/);
  assert.ok(
    reportSection(removed.output, 'Removed vendored runtime listing entries').includes(describedVendorEntry(entry)),
    removed.output
  );
});

// Persisted contract co-change (WEB_CHANGE_POLICY.md). The base carries the real schema module,
// small description modules in the generator's format, the two documents, and an old fixture.
// Every head edit derives from the values the real module holds today.
const CONTRACT_SCHEMA_PATH = 'gbdraw/web/js/contracts/schema.generated.js';
const CONTRACT_REQUEST_SHAPE_PATH = 'gbdraw/web/js/contracts/request-shape.generated.js';
const CONTRACT_SESSION_SHAPE_PATH = 'gbdraw/web/js/contracts/session-shape.generated.js';
const SESSION_HISTORY_PATH = 'docs/internal/SESSION_FORMAT_HISTORY.md';
const COMPATIBILITY_REFERENCE_PATH = 'docs/REFERENCE/session-and-request-compatibility.md';
const CONTRACT_SCHEMA_SOURCE = readFileSync(join(REPOSITORY_ROOT, CONTRACT_SCHEMA_PATH), 'utf8');
const PERSISTED_CONTRACT_KEY = Object.freeze({
  exemptGeneratedExports: ['EXPECTED_WEB_RUNTIME_CAPABILITIES']
});
const contractExportPattern = (name) => new RegExp(`export const ${name} = ([^;]*);`);
const contractExportLiteral = (name) => {
  const literal = contractExportPattern(name).exec(CONTRACT_SCHEMA_SOURCE)?.[1];
  assert.ok(literal, name);
  return literal;
};
const withContractExport = (source, name, value) => {
  const changed = source.replace(contractExportPattern(name), `export const ${name} = ${value};`);
  assert.notEqual(changed, source, `${name} = ${value} leaves the module unchanged`);
  return changed;
};
const CURRENT_SESSION_VERSION = Number(contractExportLiteral('SESSION_VERSION'));
const NEXT_SESSION_VERSION = CURRENT_SESSION_VERSION + 1;
const CURRENT_SUPPORTED_SESSION_VERSIONS = JSON.parse(
  contractExportLiteral('SUPPORTED_SESSION_VERSIONS').replace(/^Object\.freeze\(([\s\S]*)\)$/, '$1')
);
const NEXT_LOSAT_CACHE_SCHEMA = Number(contractExportLiteral('PROTEIN_LOSAT_CACHE_SCHEMA')) + 1;
const requestShapeModule = ({ modeType = 'Mode', modes = '"circular", "linear"', reflow = false } = {}) => {
  const fields = `Object.freeze({mode: Object.freeze({kind: "ref", name: "${modeType}"}), "x-y": Object.freeze({kind: "string"})})`;
  const separator = reflow ? ',\n      ' : ', ';
  return [
    '// @ts-check',
    '// Generated by tools/generate_web_contracts.py from describe_canonical_request(); do not edit manually.',
    '',
    '/** @type {ContractDescription} */',
    'export const CANONICAL_REQUEST_DESCRIPTION = Object.freeze({',
    '  root: "Request",',
    '  types: Object.freeze({',
    `    ${modeType}: Object.freeze({kind: "literal", values: Object.freeze([${modes}])}),`,
    `    Request: Object.freeze({kind: "object"${separator}fields: ${fields}${separator}optional: Object.freeze([])${separator}open: false})`,
    '  })',
    '});',
    ''
  ].join('\n');
};
const sessionShapeModule = ({
  versions = [CURRENT_SESSION_VERSION - 1, CURRENT_SESSION_VERSION],
  extraField = ''
} = {}) => [
  '// @ts-check',
  '// Generated by tools/generate_web_contracts.py from describe_session_admission(); do not edit manually.',
  '',
  '/** @type {BandedContractDescription} */',
  'export const SESSION_ADMISSION_DESCRIPTION = Object.freeze({',
  `  bands: Object.freeze([Object.freeze({minVersion: ${versions[0]}, versions: Object.freeze([${versions.join(', ')}]), root: "SessionCurrent"})]),`,
  '  types: Object.freeze({',
  `    SessionCurrent: Object.freeze({kind: "object", fields: Object.freeze({tree: Object.freeze({kind: "ref", name: "Node"})${extraField}}), optional: Object.freeze([]), open: true}),`,
  '    Node: Object.freeze({kind: "array", items: Object.freeze({kind: "ref", name: "Node"})})',
  '  })',
  '});',
  ''
].join('\n');
const PERSISTED_CONTRACT_BASE_FILES = Object.freeze({
  [CONTRACT_SCHEMA_PATH]: CONTRACT_SCHEMA_SOURCE,
  [CONTRACT_REQUEST_SHAPE_PATH]: requestShapeModule(),
  [CONTRACT_SESSION_SHAPE_PATH]: sessionShapeModule(),
  [SESSION_HISTORY_PATH]: `# Session format history\n\n## Unreleased: Session ${CURRENT_SESSION_VERSION} keys the Legend edits\n`,
  [COMPATIBILITY_REFERENCE_PATH]: `# Compatibility\n\nSession ${CURRENT_SESSION_VERSION}.\n`,
  [`tests/fixtures/sessions/x.v${CURRENT_SESSION_VERSION - 1}.gbdraw-session.json`]: '{}\n',
  'gbdraw/session_io.py': `CURRENT_SESSION_VERSION = ${CURRENT_SESSION_VERSION}\n`
});
const contractPolicyFile = (key) => ({
  'tools/web-change-policy.json': `${JSON.stringify(
    key === null ? cloneBudgetPolicy() : { ...cloneBudgetPolicy(), persistedContractCoChange: key },
    null,
    2
  )}\n`
});
// Commits the base (with the rule data when `policyKey` is not null), then the head files. With
// `dev`, the base moves on after the head's branch point, so the head is behind its base. With
// `main`, refs/remotes/origin/main points at the base with those files applied (`{}`: the base).
const runPersistedContractCase = ({
  policyKey = PERSISTED_CONTRACT_KEY,
  main = null,
  base = {},
  dev = null,
  head
}) => (
  withChangeBudgetRepository(({ commit, execute, git, root, write }) => {
    const writeAll = (files) => Object.entries(files).forEach(([path, content]) => (
      content === null ? rmSync(join(root, path)) : write(path, content)
    ));
    const baseFiles = { ...PERSISTED_CONTRACT_BASE_FILES, ...contractPolicyFile(policyKey), ...base };
    let mainRevision = '';
    if (main && Object.keys(main).length) {
      writeAll({ ...baseFiles, ...main });
      mainRevision = commit('main');
    }
    writeAll(baseFiles);
    const branchPoint = commit('persisted contract base');
    if (main) {
      assert.equal(git(['update-ref', 'refs/remotes/origin/main', mainRevision || branchPoint]).status, 0);
    }
    let baseRevision = branchPoint;
    if (dev) {
      writeAll(dev);
      baseRevision = commit('dev moves on');
      assert.equal(git(['checkout', '--quiet', branchPoint]).status, 0);
    }
    writeAll(head);
    return execute({ base: baseRevision, head: commit('persisted contract change') });
  })
);
const PERSISTED_CONTRACT_FAILURE = /persisted-contract\.co-change: /;
const PERSISTED_CONTRACT_REVIEW = /Material behavior\/output risk: persisted contract changed/;
const historyWith = (text) => `${PERSISTED_CONTRACT_BASE_FILES[SESSION_HISTORY_PATH]}\n${text}\n`;
const LOSAT_SCHEMA_CHANGE = Object.freeze({
  [CONTRACT_SCHEMA_PATH]: withContractExport(
    CONTRACT_SCHEMA_SOURCE,
    'PROTEIN_LOSAT_CACHE_SCHEMA',
    String(NEXT_LOSAT_CACHE_SCHEMA)
  ),
  [COMPATIBILITY_REFERENCE_PATH]: `# Compatibility\n\nLOSAT cache schema ${NEXT_LOSAT_CACHE_SCHEMA}.\n`,
  'tests/fixtures/losat-cache.json': `{"schema": ${NEXT_LOSAT_CACHE_SCHEMA}}\n`
});

test('persisted contract: a changed schema constant without the history fails the Gate', () => {
  const result = runPersistedContractCase({ head: LOSAT_SCHEMA_CHANGE });
  assertVendorGateFailure(result, PERSISTED_CONTRACT_FAILURE, 'missing history');
  const failure = result.output.split('\n').find((line) => PERSISTED_CONTRACT_FAILURE.test(line));
  assert.match(failure, /schema-constant: PROTEIN_LOSAT_CACHE_SCHEMA/);
  assert.ok(failure.includes(SESSION_HISTORY_PATH), failure);
  assert.ok(!failure.includes(COMPATIBILITY_REFERENCE_PATH), failure);
  assert.ok(!failure.includes('tests/fixtures/'), failure);
  assert.ok(!failure.includes('behind dev'), failure);

  const complete = runPersistedContractCase({
    head: { ...LOSAT_SCHEMA_CHANGE, [SESSION_HISTORY_PATH]: historyWith('LOSAT cache schema changed.') }
  });
  assert.equal(complete.status, 0, complete.output);
  assert.match(complete.output, /Review: \*\*REQUIRED\*\*/);
  assert.match(complete.output, /persisted contract changed \(schema-constant: PROTEIN_LOSAT_CACHE_SCHEMA\)/);
});

test('persisted contract: a refactor that keeps every description passes without the rule', () => {
  const reflowedSupported = withContractExport(
    CONTRACT_SCHEMA_SOURCE,
    'SUPPORTED_SESSION_VERSIONS',
    `Object.freeze([\n  ${CURRENT_SUPPORTED_SESSION_VERSIONS.join(',\n  ')}\n])`
  ).replace(
    '// gbdraw/session_io.py: CURRENT_SESSION_VERSION',
    '// gbdraw/session_version.py: CURRENT_SESSION_VERSION'
  );
  const result = runPersistedContractCase({
    head: {
      'gbdraw/session_io.py': 'from gbdraw.session_version import CURRENT_SESSION_VERSION\n',
      [CONTRACT_SCHEMA_PATH]: reflowedSupported,
      [CONTRACT_REQUEST_SHAPE_PATH]: requestShapeModule({ modeType: 'RenderMode', reflow: true })
    }
  });
  assert.equal(result.status, 0, result.output);
  assert.doesNotMatch(result.output, PERSISTED_CONTRACT_FAILURE);
  assert.doesNotMatch(result.output, /persisted contract changed|new generated contract export|renamed \(/);
});

test('persisted contract: the rule is off while the base policy has no rule data', () => {
  const result = runPersistedContractCase({ policyKey: null, head: LOSAT_SCHEMA_CHANGE });
  assert.equal(result.status, 0, result.output);
  assert.doesNotMatch(result.output, PERSISTED_CONTRACT_FAILURE);
  assert.doesNotMatch(result.output, PERSISTED_CONTRACT_REVIEW);
});

test('persisted contract: a new Session version needs its heading and an old-version fixture', () => {
  const sessionBump = {
    [CONTRACT_SCHEMA_PATH]: withContractExport(
      withContractExport(CONTRACT_SCHEMA_SOURCE, 'SESSION_VERSION', String(NEXT_SESSION_VERSION)),
      'SUPPORTED_SESSION_VERSIONS',
      `Object.freeze([${[...CURRENT_SUPPORTED_SESSION_VERSIONS, NEXT_SESSION_VERSION].join(', ')}])`
    ),
    [CONTRACT_SESSION_SHAPE_PATH]: sessionShapeModule({
      versions: [CURRENT_SESSION_VERSION - 1, CURRENT_SESSION_VERSION, NEXT_SESSION_VERSION]
    }),
    [COMPATIBILITY_REFERENCE_PATH]: `# Compatibility\n\nSession ${NEXT_SESSION_VERSION}.\n`,
    [`tests/fixtures/sessions/x.v${NEXT_SESSION_VERSION}.gbdraw-session.json`]: '{}\n'
  };

  const noHeading = runPersistedContractCase({
    head: { ...sessionBump, [SESSION_HISTORY_PATH]: historyWith(`Session ${NEXT_SESSION_VERSION} in prose only.`) }
  });
  assertVendorGateFailure(noHeading, PERSISTED_CONTRACT_FAILURE, 'no heading');
  assert.match(noHeading.output, /session-version: SESSION_VERSION, SUPPORTED_SESSION_VERSIONS/);
  assert.ok(noHeading.output.includes(
    `a "## " heading in ${SESSION_HISTORY_PATH} naming Session ${NEXT_SESSION_VERSION}`
  ), noHeading.output);

  const withHeading = {
    ...sessionBump,
    [SESSION_HISTORY_PATH]: historyWith(`## Unreleased: Session ${NEXT_SESSION_VERSION} adds a key`)
  };
  const noOldFixture = runPersistedContractCase({ head: withHeading });
  assertVendorGateFailure(noOldFixture, PERSISTED_CONTRACT_FAILURE, 'no old fixture');
  assert.ok(noOldFixture.output.includes(
    `a fixture of Session ${CURRENT_SESSION_VERSION} under tests/fixtures/sessions/`
  ), noOldFixture.output);
  assert.doesNotMatch(noOldFixture.output, /naming Session /);

  const complete = runPersistedContractCase({
    head: {
      ...withHeading,
      [`tests/fixtures/sessions/x.v${CURRENT_SESSION_VERSION}-schema9.gbdraw-session.json.gz`]: 'old\n'
    }
  });
  assert.equal(complete.status, 0, complete.output);
  assert.match(complete.output, /Review: \*\*REQUIRED\*\*/);
  assert.match(
    complete.output,
    /persisted contract changed \(session-version: SESSION_VERSION, SUPPORTED_SESSION_VERSIONS; session-admission: /
  );
});

test('persisted contract: an exempt export and a new export need review only', () => {
  const capabilities = contractExportLiteral('EXPECTED_WEB_RUNTIME_CAPABILITIES');
  const renderProtocol = Number(/renderProtocol: (\d+)/.exec(capabilities)?.[1]);
  assert.ok(Number.isInteger(renderProtocol), capabilities);
  const schema = `${withContractExport(
    CONTRACT_SCHEMA_SOURCE,
    'EXPECTED_WEB_RUNTIME_CAPABILITIES',
    capabilities.replace(`renderProtocol: ${renderProtocol}`, `renderProtocol: ${renderProtocol + 1}`)
  )}\n// gbdraw/session_io.py: NEW_SCHEMA\nexport const NEW_SCHEMA = "new";\n`;
  const result = runPersistedContractCase({ head: { [CONTRACT_SCHEMA_PATH]: schema } });
  assert.equal(result.status, 0, result.output);
  assert.doesNotMatch(result.output, PERSISTED_CONTRACT_REVIEW);
  assert.match(result.output, /new generated contract export\(s\) \(NEW_SCHEMA\)/);
});

test('persisted contract: malformed rule data fails closed in the base and the candidate', () => {
  const trusted = runPersistedContractCase({
    policyKey: { exemptGeneratedExports: 'EXPECTED_WEB_RUNTIME_CAPABILITIES' },
    head: { 'gbdraw/web/js/app/secondary.js': 'export const editSecondaryOwner = () => 2;\n' }
  });
  assertVendorGateFailure(trusted, /trusted persisted-contract policy: /, 'trusted base');
  const trustedUnit = runPersistedContractCase({
    policyKey: { ...PERSISTED_CONTRACT_KEY, sessionVersionUnit: 'dev' },
    head: { 'gbdraw/web/js/app/secondary.js': 'export const editSecondaryOwner = () => 2;\n' }
  });
  assertVendorGateFailure(
    trustedUnit,
    /trusted persisted-contract policy: .*sessionVersionUnit must be "main"/,
    'trusted unit'
  );

  const candidates = [
    [{ exemptGeneratedExports: ['SESSION_VERSION'] }, /exempts the version export SESSION_VERSION/],
    [{ exemptGeneratedExports: ['NOT_AN_EXPORT'] }, /exempts NOT_AN_EXPORT, which .* does not export/],
    [{ exemptGeneratedExports: [], note: 'extra' }, /must be an object with the key exemptGeneratedExports and/],
    [{ sessionVersionUnit: 'main' }, /must be an object with the key exemptGeneratedExports and/],
    [{ exemptGeneratedExports: [], sessionVersionUnit: 'dev' }, /sessionVersionUnit must be "main"/],
    [{ exemptGeneratedExports: ['lower', 'lower'] }, /must be unique export names/]
  ];
  candidates.forEach(([key, pattern]) => {
    const candidate = runPersistedContractCase({ policyKey: null, head: contractPolicyFile(key) });
    assertVendorGateFailure(candidate, /candidate persisted-contract policy: /, `candidate ${pattern}`);
    assert.match(candidate.output, pattern);
  });

  const activation = runPersistedContractCase({ policyKey: null, head: contractPolicyFile(PERSISTED_CONTRACT_KEY) });
  assert.equal(activation.status, 0, activation.output);
  assert.match(activation.output, /persisted-contract rule turned on/);
  assert.match(activation.output, /persisted-contract exemptions changed \(added=1, removed=0\)/);
});

test('persisted contract: removing the rule data needs governance review', () => {
  [PERSISTED_CONTRACT_KEY, { exemptGeneratedExports: [] }].forEach((policyKey) => {
    const result = runPersistedContractCase({ policyKey, head: contractPolicyFile(null) });
    assert.equal(result.status, 0, result.output);
    assert.match(result.output, /Governance and authority: persisted-contract rule turned off/);
  });
});

test('persisted contract: a changed Session admission description needs the compatibility reference', () => {
  const result = runPersistedContractCase({
    head: {
      [CONTRACT_SESSION_SHAPE_PATH]: sessionShapeModule({ extraField: ', note: Object.freeze({kind: "string"})' }),
      [SESSION_HISTORY_PATH]: historyWith('The current Session admits a note.'),
      'tests/fixtures/session_admission_vectors.json': '{}\n'
    }
  });
  assertVendorGateFailure(result, PERSISTED_CONTRACT_FAILURE, 'missing reference');
  assert.match(result.output, /session-admission: gbdraw\/web\/js\/contracts\/session-shape\.generated\.js/);
  assert.ok(result.output.includes(`${COMPATIBILITY_REFERENCE_PATH} in the diff`), result.output);
});

test('persisted contract: a branch behind its base is measured from the merge base', () => {
  // The base gained an export and a history entry after the branch point; the branch did not.
  const devAddsExport = {
    [CONTRACT_SCHEMA_PATH]: `${CONTRACT_SCHEMA_SOURCE}\n// gbdraw/session_io.py: NEW_SCHEMA\nexport const NEW_SCHEMA = "new";\n`,
    [SESSION_HISTORY_PATH]: historyWith('NEW_SCHEMA added.'),
    [COMPATIBILITY_REFERENCE_PATH]: '# Compatibility\n\nNEW_SCHEMA.\n',
    'tests/fixtures/new-schema.json': '{}\n'
  };
  const unrelated = runPersistedContractCase({
    dev: devAddsExport,
    head: { 'gbdraw/web/js/app/secondary.js': 'export const editSecondaryOwner = () => 2;\n' }
  });
  assert.equal(unrelated.status, 0, unrelated.output);
  assert.doesNotMatch(unrelated.output, PERSISTED_CONTRACT_FAILURE);
  assert.doesNotMatch(unrelated.output, /persisted contract changed|new generated contract export/);

  // The branch's own change still needs its own history entry; the base's entry does not count.
  const ownChange = runPersistedContractCase({ dev: devAddsExport, head: LOSAT_SCHEMA_CHANGE });
  assertVendorGateFailure(ownChange, PERSISTED_CONTRACT_FAILURE, 'behind with its own change');
  assert.match(ownChange.output, /schema-constant: PROTEIN_LOSAT_CACHE_SCHEMA\) without: docs\/internal\/SESSION_FORMAT_HISTORY\.md in the diff$/m);
});

test('persisted contract: a renamed export with the same value needs review only', () => {
  const rename = (from, to) => {
    const source = CONTRACT_SCHEMA_SOURCE.replace(
      new RegExp(`^export const ${from} = `, 'm'),
      `export const ${to} = `
    );
    assert.notEqual(source, CONTRACT_SCHEMA_SOURCE, from);
    return source;
  };
  const renamed = runPersistedContractCase({
    head: { [CONTRACT_SCHEMA_PATH]: rename('PROTEIN_LOSAT_CACHE_SCHEMA', 'PROTEIN_LOSAT_CACHE_VERSION') }
  });
  assert.equal(renamed.status, 0, renamed.output);
  assert.match(
    renamed.output,
    /generated contract export\(s\) renamed \(PROTEIN_LOSAT_CACHE_SCHEMA -> PROTEIN_LOSAT_CACHE_VERSION\)/
  );
  assert.doesNotMatch(renamed.output, PERSISTED_CONTRACT_REVIEW);
  assert.doesNotMatch(renamed.output, /new generated contract export/);

  // The checker reads version exports by name, so renaming any of them stays a failure, even
  // with every document and fixture present.
  const documented = {
    [SESSION_HISTORY_PATH]: historyWith('A version export renamed.'),
    [COMPATIBILITY_REFERENCE_PATH]: '# Compatibility\n\nRenamed.\n',
    'tests/fixtures/renamed.json': '{}\n'
  };
  const assertVersionRenameFails = (name, source, extra = {}) => {
    const result = runPersistedContractCase({ head: { ...documented, [CONTRACT_SCHEMA_PATH]: source, ...extra } });
    assertVendorGateFailure(result, PERSISTED_CONTRACT_FAILURE, `renamed ${name}`);
    assert.ok(
      result.output.includes(`a renamed or removed version export (${name}) needs a checker change`),
      result.output
    );
  };
  assertVersionRenameFails('SESSION_VERSION', rename('SESSION_VERSION', 'CURRENT_SESSION_VERSION'));
  assertVersionRenameFails(
    'SUPPORTED_CANONICAL_REQUEST_SCHEMAS',
    rename('SUPPORTED_CANONICAL_REQUEST_SCHEMAS', 'SUPPORTED_REQUEST_SCHEMAS')
  );
  // A version export that keeps its name but changes shape is unreadable to the checker too.
  assertVersionRenameFails(
    'SESSION_VERSION',
    withContractExport(CONTRACT_SCHEMA_SOURCE, 'SESSION_VERSION', `Object.freeze([${NEXT_SESSION_VERSION}])`)
  );
  assertVersionRenameFails(
    'SUPPORTED_SESSION_VERSIONS',
    withContractExport(CONTRACT_SCHEMA_SOURCE, 'SUPPORTED_SESSION_VERSIONS', String(NEXT_SESSION_VERSION))
  );
  // A new Session version whose supported list is renamed cannot skip the old-version fixture.
  assertVersionRenameFails(
    'SUPPORTED_SESSION_VERSIONS',
    withContractExport(
      withContractExport(
        rename('SUPPORTED_SESSION_VERSIONS', 'SESSION_VERSIONS_SUPPORTED'),
        'SESSION_VERSIONS_SUPPORTED',
        `Object.freeze([${[...CURRENT_SUPPORTED_SESSION_VERSIONS, NEXT_SESSION_VERSION].join(', ')}])`
      ),
      'SESSION_VERSION',
      String(NEXT_SESSION_VERSION)
    ),
    {
      [SESSION_HISTORY_PATH]: historyWith(`## Unreleased: Session ${NEXT_SESSION_VERSION} adds a key`),
      [`tests/fixtures/sessions/x.v${NEXT_SESSION_VERSION}.gbdraw-session.json`]: '{}\n'
    }
  );
});

// Session version unit: main's SESSION_VERSION (refs/remotes/origin/main) is the unit.
const SESSION_UNIT_KEY = Object.freeze({ ...PERSISTED_CONTRACT_KEY, sessionVersionUnit: 'main' });
const SESSION_UNIT_FAILURE = /persisted-contract\.fold: /;
const SESSION_UNIT_REVIEW = /the head and main both write Session /;
const sessionVersionSchema = (version) => withContractExport(
  withContractExport(CONTRACT_SCHEMA_SOURCE, 'SESSION_VERSION', String(version)),
  'SUPPORTED_SESSION_VERSIONS',
  `Object.freeze([${[...CURRENT_SUPPORTED_SESSION_VERSIONS, NEXT_SESSION_VERSION, version]
    .filter((item, index, items) => item <= version && items.indexOf(item) === index).join(', ')}])`
);
// Raises the Session version from `from` to `to` with every co-change the Gate asks for.
const sessionRaise = (from, to) => ({
  [CONTRACT_SCHEMA_PATH]: sessionVersionSchema(to),
  [SESSION_HISTORY_PATH]: historyWith(`## Unreleased: Session ${to} adds a key`),
  [COMPATIBILITY_REFERENCE_PATH]: `# Compatibility\n\nSession ${to}.\n`,
  [`tests/fixtures/sessions/x.v${from}.gbdraw-session.json`]: '{}\n'
});
const SESSION_ADMISSION_CHANGE = Object.freeze({
  [CONTRACT_SESSION_SHAPE_PATH]: sessionShapeModule({ extraField: ', note: Object.freeze({kind: "string"})' }),
  [SESSION_HISTORY_PATH]: historyWith('The current Session admits a note.'),
  [COMPATIBILITY_REFERENCE_PATH]: '# Compatibility\n\nA note.\n',
  'tests/fixtures/session_admission_vectors.json': '{}\n'
});

test('persisted contract: a Session version above main takes later breaking changes', () => {
  // dev already writes NEXT, which main (CURRENT) has not reached; raising it again fails.
  const raiseAgain = {
    main: { [CONTRACT_SCHEMA_PATH]: CONTRACT_SCHEMA_SOURCE },
    base: { [CONTRACT_SCHEMA_PATH]: sessionVersionSchema(NEXT_SESSION_VERSION) },
    head: sessionRaise(NEXT_SESSION_VERSION, NEXT_SESSION_VERSION + 1)
  };
  const folded = runPersistedContractCase({ ...raiseAgain, policyKey: SESSION_UNIT_KEY });
  assertVendorGateFailure(folded, SESSION_UNIT_FAILURE, 'raised again');
  assert.ok(folded.output.includes(
    `SESSION_VERSION changes ${NEXT_SESSION_VERSION} -> ${NEXT_SESSION_VERSION + 1}, `
    + `but main (refs/remotes/origin/main) writes Session ${CURRENT_SESSION_VERSION}, so the head may write `
    + `only Session ${CURRENT_SESSION_VERSION} or ${NEXT_SESSION_VERSION}`
  ), folded.output);
  assert.doesNotMatch(folded.output, PERSISTED_CONTRACT_FAILURE);
  assert.match(folded.output, new RegExp(`Session version: main ${CURRENT_SESSION_VERSION}, head `));

  // A jump of two from main fails the same way.
  const jump = runPersistedContractCase({
    policyKey: SESSION_UNIT_KEY,
    main: {},
    head: sessionRaise(CURRENT_SESSION_VERSION, NEXT_SESSION_VERSION + 1)
  });
  assertVendorGateFailure(jump, SESSION_UNIT_FAILURE, 'jump of two');

  // Lowering it below main fails too.
  const lowered = runPersistedContractCase({
    policyKey: SESSION_UNIT_KEY,
    main: {},
    base: { [CONTRACT_SCHEMA_PATH]: sessionVersionSchema(NEXT_SESSION_VERSION) },
    head: {
      [CONTRACT_SCHEMA_PATH]: CONTRACT_SCHEMA_SOURCE,
      [SESSION_HISTORY_PATH]: historyWith('Session version lowered.'),
      [COMPATIBILITY_REFERENCE_PATH]: '# Compatibility\n\nLowered.\n',
      'tests/fixtures/sessions/lowered.json': '{}\n'
    }
  });
  assertVendorGateFailure(lowered, SESSION_UNIT_FAILURE, 'lowered below main');
  assert.doesNotMatch(lowered.output, PERSISTED_CONTRACT_FAILURE);

  // A version already more than one above main that the pull request leaves alone passes.
  const inherited = runPersistedContractCase({
    policyKey: SESSION_UNIT_KEY,
    main: { [CONTRACT_SCHEMA_PATH]: CONTRACT_SCHEMA_SOURCE },
    base: { [CONTRACT_SCHEMA_PATH]: sessionVersionSchema(NEXT_SESSION_VERSION + 1) },
    head: SESSION_ADMISSION_CHANGE
  });
  assert.equal(inherited.status, 0, inherited.output);
  assert.doesNotMatch(inherited.output, /persisted-contract\.fold|both write Session/);

  // Without sessionVersionUnit in the base policy the rule is off.
  const inert = runPersistedContractCase(raiseAgain);
  assert.equal(inert.status, 0, inert.output);
  assert.doesNotMatch(inert.output, /persisted-contract\.fold|both write Session|Session version: main/);

  // Folding a change into the version above main passes without the main question.
  const fold = runPersistedContractCase({
    policyKey: SESSION_UNIT_KEY,
    main: raiseAgain.main,
    base: raiseAgain.base,
    head: SESSION_ADMISSION_CHANGE
  });
  assert.equal(fold.status, 0, fold.output);
  assert.match(fold.output, PERSISTED_CONTRACT_REVIEW);
  assert.doesNotMatch(fold.output, SESSION_UNIT_REVIEW);
});

test('persisted contract: raising the Session version above main passes', () => {
  const result = runPersistedContractCase({
    policyKey: SESSION_UNIT_KEY,
    main: {},
    head: sessionRaise(CURRENT_SESSION_VERSION, NEXT_SESSION_VERSION)
  });
  assert.equal(result.status, 0, result.output);
  assert.match(result.output, /persisted contract changed \(session-version: /);
  assert.doesNotMatch(result.output, /persisted-contract\.fold|both write Session/);
});

test('persisted contract: a change while the head writes main\'s Session version asks whether it breaks main', () => {
  const result = runPersistedContractCase({ policyKey: SESSION_UNIT_KEY, main: {}, head: SESSION_ADMISSION_CHANGE });
  assert.equal(result.status, 0, result.output);
  assert.ok(result.output.includes(
    `the head and main both write Session ${CURRENT_SESSION_VERSION}; confirm that main reads what the head `
    + `writes the same way, or raise SESSION_VERSION to ${NEXT_SESSION_VERSION}`
  ), result.output);

  // A contract change that leaves the schema module alone reads the head's version too.
  const losat = runPersistedContractCase({
    policyKey: SESSION_UNIT_KEY,
    main: {},
    head: { ...LOSAT_SCHEMA_CHANGE, [SESSION_HISTORY_PATH]: historyWith('LOSAT cache schema raised.') }
  });
  assert.equal(losat.status, 0, losat.output);
  assert.match(losat.output, SESSION_UNIT_REVIEW);

  // A new export alone is not a persisted contract change and asks nothing.
  const added = runPersistedContractCase({
    policyKey: SESSION_UNIT_KEY,
    head: {
      [CONTRACT_SCHEMA_PATH]: `${CONTRACT_SCHEMA_SOURCE}\n// gbdraw/session_io.py: NEW_SCHEMA\n`
        + 'export const NEW_SCHEMA = "new";\n'
    }
  });
  assert.equal(added.status, 0, added.output);
  assert.match(added.output, /new generated contract export\(s\) \(NEW_SCHEMA\)/);
  assert.doesNotMatch(added.output, /persisted-contract\.fold|both write Session/);
});

test('persisted contract: main\'s Session version that cannot be read fails closed', () => {
  const result = runPersistedContractCase({ policyKey: SESSION_UNIT_KEY, head: SESSION_ADMISSION_CHANGE });
  assertVendorGateFailure(result, SESSION_UNIT_FAILURE, 'no origin/main');
  assert.match(result.output, /SESSION_VERSION cannot be read at main \(refs\/remotes\/origin\/main; fetch main/);
  assert.match(result.output, /Session version: main unreadable/);

  // An unreadable head names the head, not main.
  const head = runPersistedContractCase({
    policyKey: SESSION_UNIT_KEY,
    main: {},
    head: {
      ...SESSION_ADMISSION_CHANGE,
      [CONTRACT_SCHEMA_PATH]: withContractExport(
        CONTRACT_SCHEMA_SOURCE,
        'SESSION_VERSION',
        `Object.freeze([${NEXT_SESSION_VERSION}])`
      )
    }
  });
  assertVendorGateFailure(head, SESSION_UNIT_FAILURE, 'unreadable head');
  assert.match(head.output, /persisted-contract\.fold: SESSION_VERSION cannot be read at the head/);
});

test('persisted contract: setting or removing the Session version unit needs governance review', () => {
  const set = runPersistedContractCase({ head: contractPolicyFile(SESSION_UNIT_KEY) });
  assert.equal(set.status, 0, set.output);
  assert.match(set.output, /Governance and authority: Session version unit changed \(none -> main\)/);
  assert.match(set.output, /Candidate rule data: valid; exempt: [A-Z_]+; Session version unit: main/);

  const removed = runPersistedContractCase({
    policyKey: SESSION_UNIT_KEY,
    head: contractPolicyFile(PERSISTED_CONTRACT_KEY)
  });
  assert.equal(removed.status, 0, removed.output);
  assert.match(removed.output, /Governance and authority: Session version unit changed \(main -> none\)/);
});

test('a repeated key in the change policy fails closed', () => {
  const policyBody = JSON.stringify(cloneBudgetPolicy(), null, 2).slice(0, -2);
  const entry = JSON.stringify(vendorChangeEntry());
  [
    `${policyBody},\n  "allowedVendorChanges": [],\n  "allowedVendorChanges": [${entry}]\n}\n`,
    `${policyBody},\n  "allowedVendorChanges": [],\n  "allowedVendor\\u0043hanges": [${entry}]\n}\n`
  ].forEach((source, index) => {
    const result = runChangeBudgetCase((write) => write('tools/web-change-policy.json', source));
    assertVendorGateFailure(result, /repeats the key\(s\) allowedVendorChanges/, `variant ${index}`);
  });
});

test('a vendored runtime listing entry that cannot apply to the base fails where it is added', () => {
  const cases = [
    [vendorChangeEntry({ path: 'gbdraw/web/vendor/missing.js' }), /missing\.js: replace for a path the base lacks/],
    [
      vendorChangeEntry({ path: 'gbdraw/web/vendor/missing.js', change: 'remove' }),
      /missing\.js: remove for a path the base lacks/
    ],
    [vendorChangeEntry({ change: 'add' }), /library\.js: add for a path the base already has/],
    [vendorChangeEntry({ change: 'remove' }), /library\.js: remove lists sha256 [0-9a-f]{64}, but the base bytes hash to/],
    [
      vendorChangeEntry({ sha256: sha256Hex(BUDGET_FIXTURE[VENDOR_LIBRARY_PATH]) }),
      /library\.js: replace lists the bytes the base already has/
    ]
  ];
  cases.forEach(([entry, pattern]) => {
    const result = runChangeBudgetCase((write) => writeVendorChangeListing(write, [entry]));
    assertVendorGateFailure(result, /candidate vendored runtime listing: /, String(pattern));
    assert.match(result.output, pattern);
  });
});

test('vendored runtime listing edge cases fail closed', () => {
  // Bootstrap: the base has no policy file, so the head policy supplies no vendored entries.
  const bootstrap = withChangeBudgetRepository(({ commit, execute, root, write }) => {
    rmSync(join(root, 'tools/web-change-policy.json'));
    const base = commit('base without a change policy');
    writeVendorChangeListing(write, [vendorChangeEntry()]);
    write(VENDOR_LIBRARY_PATH, VENDOR_LIBRARY_UPDATE);
    const head = commit('policy and vendored change');
    return execute({ base, head });
  });
  assertVendorGateFailure(bootstrap, /gbdraw\/web\/vendor\/library\.js: not listed/, 'bootstrap');

  // A rename is a remove and an add; each side needs its own entry.
  const renamedPath = 'gbdraw/web/vendor/library-renamed.js';
  const original = BUDGET_FIXTURE[VENDOR_LIBRARY_PATH];
  const addRenamed = vendorChangeEntry({ path: renamedPath, change: 'add', sha256: sha256Hex(original) });
  const rename = (write) => {
    write(VENDOR_LIBRARY_PATH, null);
    write(renamedPath, original);
  };
  const renameHalf = runListedVendorChangeCase([addRenamed], rename);
  assertVendorGateFailure(renameHalf, /library\.js: not listed/, 'rename with one entry');
  const renameBoth = runListedVendorChangeCase(
    [vendorChangeEntry({ change: 'remove', sha256: sha256Hex(original) }), addRenamed],
    rename
  );
  assert.equal(renameBoth.status, 0, renameBoth.output);

  const replaceForAdd = runListedVendorChangeCase(
    [vendorChangeEntry({ path: 'gbdraw/web/vendor/new.js' })],
    (write) => write('gbdraw/web/vendor/new.js', VENDOR_LIBRARY_UPDATE)
  );
  assertVendorGateFailure(replaceForAdd, /new\.js: listed as replace, but the change is add/, 'replace for add');

  const typeChange = withChangeBudgetRepository(({ commit, execute, root, write }) => {
    writeVendorChangeListing(write, [vendorChangeEntry()]);
    const base = commit('listing');
    rmSync(join(root, VENDOR_LIBRARY_PATH));
    symlinkSync('library-target.js', join(root, VENDOR_LIBRARY_PATH));
    const head = commit('vendored file becomes a symlink');
    return execute({ base, head });
  });
  assertVendorGateFailure(
    typeChange,
    /library\.js: listed as replace, but the change is unsupported status T/,
    'type change'
  );

  const malformedTrusted = runListedVendorChangeCase(
    [vendorChangeEntry({ sha256: 'ABC' })],
    (write) => write(VENDOR_LIBRARY_PATH, VENDOR_LIBRARY_UPDATE)
  );
  assertVendorGateFailure(malformedTrusted, /trusted vendored runtime listing: /, 'malformed trusted');
  assert.match(malformedTrusted.output, /library\.js: not listed/);

  [
    [{}, /allowedVendorChanges must be an array/],
    [[1], /allowedVendorChanges\[0\] must be an object/]
  ].forEach(([listing, pattern]) => {
    const result = runChangeBudgetCase((write) => writeVendorChangeListing(write, listing));
    assertVendorGateFailure(result, pattern, String(pattern));
  });
});

// The listing merges first (fork), the listed change merges after it (base); the head either
// forks from the listing commit (behind) or builds on the base (up to date).
const runForkedVendorChangeCase = ({ entry, merged, branch, behind }) => withChangeBudgetRepository(
  ({ commit, execute, git, root, write }) => {
    const writeOrRemove = (path, content) => (
      content === null ? rmSync(join(root, path)) : write(path, content)
    );
    writeVendorChangeListing(write, [entry]);
    const fork = commit('authority-only vendored runtime listing');
    let head = '';
    if (behind) {
      branch(writeOrRemove);
      head = commit('branch forked before the vendored change merged');
      assert.equal(git(['checkout', '--quiet', fork]).status, 0);
    }
    merged(writeOrRemove);
    const base = commit('listed vendored change merged');
    if (!behind) {
      branch(writeOrRemove);
      head = commit('branch on the current base');
    }
    return execute({ base, head });
  }
);
const BEHIND_BASE = /the branch is behind its base: update it from dev/;

test('only a branch that is really behind its base is told to update', () => {
  const prodPath = 'gbdraw/web/vendor/library.prod.js';
  const unrelated = (write) => write(
    'gbdraw/web/js/app/secondary.js',
    'export const editSecondaryOwner = () => 2;\n'
  );
  [
    {
      name: 'add',
      entry: vendorChangeEntry({ path: prodPath, change: 'add' }),
      merged: (write) => write(prodPath, VENDOR_LIBRARY_UPDATE)
    },
    {
      name: 'replace',
      entry: vendorChangeEntry(),
      merged: (write) => write(VENDOR_LIBRARY_PATH, VENDOR_LIBRARY_UPDATE)
    },
    {
      name: 'remove',
      entry: vendorChangeEntry({
        change: 'remove',
        sha256: sha256Hex(BUDGET_FIXTURE[VENDOR_LIBRARY_PATH])
      }),
      merged: (write) => write(VENDOR_LIBRARY_PATH, null)
    }
  ].forEach(({ name, entry, merged }) => {
    const result = runForkedVendorChangeCase({ entry, merged, branch: unrelated, behind: true });
    assertVendorGateFailure(result, BEHIND_BASE, `behind ${name}`);
  });

  // Up to date: the base merged the listed replace, and the head changes the file again.
  const tampered = runForkedVendorChangeCase({
    entry: vendorChangeEntry(),
    merged: (write) => write(VENDOR_LIBRARY_PATH, VENDOR_LIBRARY_UPDATE),
    branch: (write) => write(VENDOR_LIBRARY_PATH, 'EVIL\n'),
    behind: false
  });
  assertVendorGateFailure(
    tampered,
    new RegExp(`library\\.js: sha256 ${sha256Hex('EVIL\n')} does not match the listed`),
    'up-to-date tampering'
  );
  assert.doesNotMatch(tampered.output, BEHIND_BASE);

  // Up to date: the head deletes a file that a spent add entry still lists.
  const spentAdd = runForkedVendorChangeCase({
    entry: vendorChangeEntry({ path: prodPath, change: 'add' }),
    merged: (write) => write(prodPath, VENDOR_LIBRARY_UPDATE),
    branch: (write) => write(prodPath, null),
    behind: false
  });
  assertVendorGateFailure(
    spentAdd,
    /library\.prod\.js: listed as add, but the change is remove/,
    'spent add deleted'
  );
  assert.doesNotMatch(spentAdd.output, BEHIND_BASE);
});

test('an add entry for a path that is a directory in the base fails', () => {
  const result = withChangeBudgetRepository(({ commit, execute, write }) => {
    write('gbdraw/web/vendor/bundle/index.js', 'globalThis.bundle = true;\n');
    const base = commit('vendored directory');
    writeVendorChangeListing(write, [vendorChangeEntry({ path: 'gbdraw/web/vendor/bundle', change: 'add' })]);
    const head = commit('add entry for the directory path');
    return execute({ base, head });
  });
  assertVendorGateFailure(
    result,
    /bundle: add for a path that is not a file in the base/,
    'directory add'
  );
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

test('the Product Contract authority path is exact and isolated', () => {
  const authorityOnly = runChangeBudgetCase((write) => {
    write(PRODUCT_CONTRACT_AUTHORITY_PATH, '# Option Integrity Product Contract\n');
  });
  assert.equal(authorityOnly.status, 0, authorityOnly.output);
  assert.match(authorityOnly.output, /Gate: \*\*PASS\*\*/);
  assert.match(authorityOnly.output, /Review: \*\*REQUIRED\*\*/);
  assert.match(
    authorityOnly.output,
    /Changed Product Impact authority paths: `docs\/internal\/OPTION_INTEGRITY_PRODUCT_CONTRACT\.md`/
  );
  assert.match(
    reportSection(authorityOnly.output, 'Guard files touched'),
    /docs\/internal\/OPTION_INTEGRITY_PRODUCT_CONTRACT\.md/
  );

  const similarPath = 'docs/internal/OPTION_INTEGRITY_PRODUCT_CONTRACT_COPY.md';
  const lookalikeOnly = runChangeBudgetCase((write) => {
    write(similarPath, '# Similar name, unrelated path\n');
  });
  assert.equal(lookalikeOnly.status, 0, lookalikeOnly.output);
  assert.match(lookalikeOnly.output, /Review: \*\*CLEAR\*\*/);
  assert.equal(reportSection(lookalikeOnly.output, 'Guard files touched').trim(), '- None');

  const ambiguousBundle = runChangeBudgetCase((write) => {
    write(PRODUCT_CONTRACT_AUTHORITY_PATH, '# Option Integrity Product Contract\n');
    write(similarPath, '# Similar name, unrelated path\n');
  });
  assert.equal(ambiguousBundle.status, 1, ambiguousBundle.output);
  assert.match(ambiguousBundle.output, /Gate: \*\*FAIL\*\*/);
  assert.match(
    ambiguousBundle.output,
    /Product Contract authority changes must be isolated from other changed paths: docs\/internal\/OPTION_INTEGRITY_PRODUCT_CONTRACT_COPY\.md/
  );
});

const PRODUCT_CONTRACT_BASE_SOURCE = '# Option Integrity Product Contract\n\n### PD-OI-900: Fixture record\n';
const PRODUCT_CONTRACT_CHANGED_SOURCE = `${PRODUCT_CONTRACT_BASE_SOURCE}\nRevised outcome with its receipt.\n`;
const PRODUCT_CONTRACT_CO_CHANGE_REVIEW =
  /static Product Contract co-change requires human review that each changed record serializes the explicit Product Decision Owner receipt/;
const writeProductContractCoChangeRuntime = (write) => {
  write('gbdraw/web/js/services/session-file.js', 'export const readSession = () => ({ version: 2 });\n');
  write('tests/web/session-file.test.mjs', "test('reads the revised Session', () => {});\n");
  write('docs/SESSION_COMPATIBILITY.md', '# Session compatibility\n\nRevised outcome.\n');
};
// The Contract exists at the base, as it does in the repository, so the
// candidate modifies it in place.
const runProductContractRevisionCase = (mutate) => withChangeBudgetRepository(
  ({ commit, execute, write }) => {
    write(PRODUCT_CONTRACT_AUTHORITY_PATH, PRODUCT_CONTRACT_BASE_SOURCE);
    const base = commit('Product Contract baseline');
    mutate(write);
    const head = commit('candidate');
    return execute({ base, head });
  }
);

test('a static Product Contract co-change with runtime, tests, and docs requires review only', () => {
  const coChange = runProductContractRevisionCase((write) => {
    write(PRODUCT_CONTRACT_AUTHORITY_PATH, PRODUCT_CONTRACT_CHANGED_SOURCE);
    writeProductContractCoChangeRuntime(write);
  });
  assert.equal(coChange.status, 0, coChange.output);
  assert.match(coChange.output, /Gate: \*\*PASS\*\*/);
  assert.match(coChange.output, /Review: \*\*REQUIRED\*\*/);
  assert.match(reportSection(coChange.output, 'Review reasons'), PRODUCT_CONTRACT_CO_CHANGE_REVIEW);
  assert.match(coChange.output, /Static Product Contract: co-change with its implementation; human Review verifies the receipt/);
  assert.doesNotMatch(coChange.output, /production runtime files and Web guard\/CI files changed together/);
  assert.doesNotMatch(coChange.output, /Product Contract authority changes must be isolated/);

  const contractOnly = runProductContractRevisionCase((write) => {
    write(PRODUCT_CONTRACT_AUTHORITY_PATH, PRODUCT_CONTRACT_CHANGED_SOURCE);
  });
  assert.equal(contractOnly.status, 0, contractOnly.output);
  assert.match(contractOnly.output, /Review: \*\*REQUIRED\*\*/);
  assert.doesNotMatch(contractOnly.output, PRODUCT_CONTRACT_CO_CHANGE_REVIEW);
  assert.match(contractOnly.output, /Static Product Contract: isolated change for human Review/);

  const runtimeOnly = runProductContractRevisionCase(writeProductContractCoChangeRuntime);
  assert.equal(runtimeOnly.status, 0, runtimeOnly.output);
  assert.doesNotMatch(runtimeOnly.output, PRODUCT_CONTRACT_CO_CHANGE_REVIEW);
  assert.equal(reportSection(runtimeOnly.output, 'Guard files touched').trim(), '- None');

  const lookalikeCompanion = runProductContractRevisionCase((write) => {
    write(PRODUCT_CONTRACT_AUTHORITY_PATH, PRODUCT_CONTRACT_CHANGED_SOURCE);
    write('docs/internal/OPTION_INTEGRITY_PRODUCT_CONTRACT_COPY.md', '# Similar name, ordinary document\n');
  });
  assert.equal(lookalikeCompanion.status, 0, lookalikeCompanion.output);
  assert.match(reportSection(lookalikeCompanion.output, 'Review reasons'), PRODUCT_CONTRACT_CO_CHANGE_REVIEW);
  assert.equal(
    reportSection(lookalikeCompanion.output, 'Guard files touched').trim(),
    `- ${PRODUCT_CONTRACT_AUTHORITY_PATH}`
  );
});

test('a static Product Contract co-change cannot carry any other guard or authority change', () => {
  [
    'tools/check-web-change-budget.mjs',
    'tools/web-product-impact-evaluation.mjs',
    'tools/web-architecture-detectors.mjs',
    '.github/workflows/test.yml',
    '.github/workflows/web-base-policy.yml',
    'tools/web-product-impact-map.json',
    'tools/web-product-decisions.json',
    'tools/web-architecture-rules.json',
    'tools/web-change-policy.json',
    'docs/internal/PRODUCT_IMPACT_RATCHET.md',
    'docs/internal/WEB_CHANGE_POLICY.md',
    '.github/pull_request_template.md',
    'tests/web/architecture-contracts.test.mjs'
  ].forEach((extraPath) => {
    const result = runProductContractRevisionCase((write) => {
      write(PRODUCT_CONTRACT_AUTHORITY_PATH, PRODUCT_CONTRACT_CHANGED_SOURCE);
      writeProductContractCoChangeRuntime(write);
      const baseline = BUDGET_FIXTURE[extraPath];
      write(
        extraPath,
        !baseline ? reservedPathContent(extraPath)
          : extraPath.endsWith('.json') ? `${baseline}\n` : `${baseline}\n// changed\n`
      );
    });
    assert.equal(result.status, 1, `${extraPath}\n${result.output}`);
    assert.match(result.output, /Gate: \*\*FAIL\*\*/, extraPath);
    assert.match(result.output, /production runtime files and Web guard\/CI files changed together/, extraPath);
    assert.match(result.output, /Product Contract authority changes must be isolated from other changed paths/, extraPath);
    assert.doesNotMatch(result.output, PRODUCT_CONTRACT_CO_CHANGE_REVIEW, extraPath);
  });
});

test('deleting or moving the static Product Contract never uses the co-change route', () => {
  [
    ['deleted', (write, root) => rmSync(join(root, PRODUCT_CONTRACT_AUTHORITY_PATH))],
    ['moved', (write, root) => {
      rmSync(join(root, PRODUCT_CONTRACT_AUTHORITY_PATH));
      write('docs/internal/OPTION_INTEGRITY_PRODUCT_CONTRACT_COPY.md', PRODUCT_CONTRACT_CHANGED_SOURCE);
    }]
  ].forEach(([label, mutate]) => {
    const result = withChangeBudgetRepository(({ commit, execute, root, write }) => {
      write(PRODUCT_CONTRACT_AUTHORITY_PATH, PRODUCT_CONTRACT_BASE_SOURCE);
      const base = commit('Product Contract baseline');
      mutate(write, root);
      writeProductContractCoChangeRuntime(write);
      const head = commit(`candidate ${label}`);
      return execute({ base, head });
    });
    assert.equal(result.status, 1, `${label}\n${result.output}`);
    assert.match(result.output, /production runtime files and Web guard\/CI files changed together/, label);
    assert.match(result.output, /Product Contract authority changes must be isolated from other changed paths/, label);
  });
});

test('a static Product Contract co-change cannot replace mapped hard evidence', () => {
  const hardEvidenceMap = JSON.parse(canonicalTransitionProductImpactMapSource({ enforcement: 'hard' }));
  hardEvidenceMap.concerns[0].contracts[0].ref = (
    'tests/web/contracts/product-impact-entry.test.mjs::canonical alternate entry'
  );
  withTrustedArchitectureRepository(
    (write) => {
      write(PRODUCT_CONTRACT_AUTHORITY_PATH, PRODUCT_CONTRACT_CHANGED_SOURCE);
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
      assert.match(reportSection(output, 'Review reasons'), PRODUCT_CONTRACT_CO_CHANGE_REVIEW);
      assert.doesNotMatch(output, /Product Contract authority changes must be isolated/);
    },
    {
      ...preauthorizedProductBaseFiles(`${JSON.stringify(hardEvidenceMap, null, 2)}\n`),
      [PRODUCT_CONTRACT_AUTHORITY_PATH]: PRODUCT_CONTRACT_BASE_SOURCE,
      'tests/web/contracts/product-impact-entry.test.mjs': "test('canonical alternate entry', () => {});\n"
    },
    productImpactPullRequestEnvironment()
  );
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
    'tools/web-owner-graph-detectors.mjs',
    'tools/web-architecture-evaluation.mjs',
    'tools/web-product-impact-evaluation.mjs',
    'tools/web-product-impact-decision-source.mjs',
    'tools/web-promotion-context.mjs',
    'tools/check-promotion-readiness.mjs'
  ];
  const authorityPaths = [
    'docs/internal/ARCHITECTURE_FITNESS_FUNCTION_RATCHET.md',
    PRODUCT_CONTRACT_AUTHORITY_PATH,
    'docs/internal/PRODUCT_IMPACT_RATCHET.md',
    'docs/internal/WEB_CHANGE_POLICY.md',
    'tools/web-change-policy.json',
    'tools/web-architecture-rules.json',
    'tools/web-architecture-violations.json',
    'tools/web-product-impact-map.json',
    'tools/web-product-decisions.json',
    WEB_DESIGN_RULES_PATH,
    DESIGN_RULE_GUARD_REGISTRY_PATH,
    OWNER_GRAPH_REGISTRY_PATH,
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
    'tools/web-owner-graph-detectors.mjs',
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

test('every file under .github/actions/ is a guard and authority path like a workflow', () => {
  const action = '.github/actions/example/action.yml';
  const actionYaml = 'name: Example\nruns:\n  using: composite\n  steps: []\n';
  const actionOnly = runChangeBudgetCase((write) => write(action, actionYaml));
  assert.equal(actionOnly.status, 0, actionOnly.output);
  assert.match(
    reportSection(actionOnly.output, 'Review reasons'),
    /Governance and authority: registered governance or authority files changed \(\.github\/actions\/example\/action\.yml\)/
  );
  assert.equal(reportSection(actionOnly.output, 'Guard files touched').trim(), `- ${action}`);

  // An action, like the workflow that runs it, does not change with Web runtime or the checker.
  const withRuntime = runChangeBudgetCase((write) => {
    write(action, actionYaml);
    write('gbdraw/web/js/services/session-file.js', 'export const readSession = () => ({ version: 2 });\n');
  });
  assert.equal(withRuntime.status, 1, withRuntime.output);
  assert.match(withRuntime.output, /production runtime files and Web guard\/CI files changed together/);
  const withChecker = runChangeBudgetCase((write) => {
    write(action, actionYaml);
    write('tools/web-architecture-evaluation.mjs', reservedPathContent('tools/web-architecture-evaluation.mjs'));
  });
  assert.equal(withChecker.status, 1, withChecker.output);
  assert.match(withChecker.output, /Web checker\/source parser and authority policy\/workflow files changed together/);

  // A look-alike outside the prefix stays unregistered.
  const lookalike = runChangeBudgetCase((write) => write('.github/actions-notes.md', '# Notes\n'));
  assert.doesNotMatch(lookalike.output, /registered governance or authority files changed/);
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

const withRetiredCapabilityBase = (runCase) => withChangeBudgetRepository((repository) => {
  assert.ok(!WEB_PRIVILEGED_CAPABILITY_KEYS.includes('Retired capability'));
  const policy = cloneBudgetPolicy();
  policy.allowedPrivilegedOwners['Retired capability'] = ['app/editor.js', 'app/run-analysis.js'];
  writeBudgetPolicy(repository.write, policy);
  const base = repository.commit('hold a capability key that no detector defines');
  delete policy.allowedPrivilegedOwners['Retired capability'];
  return runCase({ ...repository, base, policy });
});

test('an authority-only change removes a capability key that no detector defines', () => {
  withRetiredCapabilityBase(({ base, commit, execute, policy, write }) => {
    writeBudgetPolicy(write, policy);
    const head = commit('remove retired capability key');
    for (const environment of [{}, { WEB_ARCHITECTURE_CHANGE: 'true' }]) {
      const result = execute({ base, head, environment });
      assert.equal(result.status, 0, result.output);
      assert.match(result.output, /Gate: \*\*PASS\*\*/);
      assert.match(result.output, /Review: \*\*REQUIRED\*\*/);
      assert.match(result.output,
        /## Removed retired privileged allowlist keys\n\n- allowedPrivilegedOwners\.Retired capability\n/);
      assert.match(result.output, /allowedPrivilegedOwners\.Retired capability: app\/editor\.js/);
      assert.doesNotMatch(result.output,
        /proposed privileged capability policy is missing base allowlist keys/);
    }
  });
});

test('a capability key that a detector defines still cannot be removed', () => {
  withRetiredCapabilityBase(({ base, commit, execute, policy, write }) => {
    delete policy.allowedPrivilegedOwners['Render request'];
    writeBudgetPolicy(write, policy);
    const head = commit('remove retired and detected capability keys');
    assertNonWaivableRevisionFailure(execute, base, head, [
      /## Missing base privileged allowlist keys\n\n- allowedPrivilegedOwners\.Render request\n/,
      /proposed privileged capability policy is missing base allowlist keys/
    ]);
  });
});

test('a runtime contraction cannot also remove a retired capability key', () => {
  withRetiredCapabilityBase(({ base, commit, execute, policy, write }) => {
    removePolicyPath(policy, 'allowedPrivilegedOwners', 'Diagram Worker', 'app/editor.js');
    writeBudgetPolicy(write, policy);
    removeEditorOwnerUse(write);
    const head = commit('remove runtime owner, its permission, and a retired key');
    assertNonWaivableRevisionFailure(execute, base, head, [
      /production runtime files and Web guard\/CI files changed together/,
      /## Removed retired privileged allowlist keys\n\n- allowedPrivilegedOwners\.Retired capability\n/
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

test('in-memory Session-import snapshot keys are not Session object keys', () => {
  const configSource = (snapshotKeys, sessionKeys) => (
    'const captureSessionImportSnapshot = () => ({\n'
      + snapshotKeys.map((key) => `  ${key}: state.${key}.value,\n`).join('')
      + '  transients: captureTransientState()\n'
      + '});\n'
      + 'const captureTransientState = () => ({ pan: state.pan.value });\n'
      + 'export const buildRunStateData = () => ({\n'
      + sessionKeys.map((key) => `  ${key}: 1,\n`).join('')
      + '  version: 1\n'
      + '});\n'
  );
  const signalLine = /Architecture-bearing signals: .*session schema fields/;
  withChangeBudgetRepository(({ commit, execute, write }) => {
    write('gbdraw/web/js/services/config.js', configSource(['errorLog'], ['runId']));
    const base = commit('config baseline');

    write('gbdraw/web/js/services/config.js', configSource(['errorLog', 'newSnapshotField'], ['runId']));
    const snapshotOnly = execute({ base, head: commit('snapshot key only') });
    assert.equal(snapshotOnly.status, 0, snapshotOnly.output);
    assert.doesNotMatch(snapshotOnly.output, signalLine);
    assert.doesNotMatch(snapshotOnly.output, /newSnapshotField/);

    write('gbdraw/web/js/services/config.js', configSource(['errorLog'], ['runId', 'newSessionField']));
    const realKey = execute({ base, head: commit('real Session key') });
    assert.equal(realKey.status, 0, realKey.output);
    assert.match(realKey.output, signalLine);
    assert.match(realKey.output, /services\/config\.js: newSessionField/);
  });
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

test('trusted-base execution rejects runtime plus a candidate Product Contract and checker', () => {
  withChangeBudgetRepository(({ commit, execute, git, root, write }) => {
    const base = git(['rev-parse', 'HEAD']).stdout.trim();
    write(
      'gbdraw/web/js/services/session-file.js',
      'export const readSession = () => ({ version: 2 });\n'
    );
    write(PRODUCT_CONTRACT_AUTHORITY_PATH, '# Candidate Product Contract\n');
    write('tools/check-web-change-budget.mjs', 'process.exit(0);\n');
    const head = commit('candidate Product Contract self-authorizes runtime');

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
      /Web checker\/source parser and authority policy\/workflow files changed together/
    );
    assert.match(
      trustedBase.output,
      /Product Contract authority changes must be isolated from other changed paths/
    );
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

// Design-rule guard registry (WEB_CHANGE_POLICY.md "Design-rule co-change"):
// the registry names the allowlist and baseline literals inside the guard
// tests of the numbered rules in gbdraw/web/CLAUDE.md. A runtime diff may contract them and
// may not expand them; an expansion is an authority-only change.
const DESIGN_RULE_GUARD_TEST_PATH = 'tests/web/design-rule-guard-fixture.test.mjs';
const designRuleGuardRegistrySource = (rules, schemaVersion = 1) => `${JSON.stringify({
  schemaVersion,
  rulesPath: WEB_DESIGN_RULES_PATH,
  rules
}, null, 2)}\n`;
const DESIGN_RULE_REGISTRY = designRuleGuardRegistrySource([{
  id: 'R10',
  heading: 'R10: A watcher does not repair state',
  guards: [DESIGN_RULE_GUARD_TEST_PATH],
  allowlists: [
    { path: DESIGN_RULE_GUARD_TEST_PATH, symbol: 'RECONCILE', kind: 'set' },
    { path: DESIGN_RULE_GUARD_TEST_PATH, symbol: 'WRITERS', kind: 'writer-map' },
    { path: DESIGN_RULE_GUARD_TEST_PATH, symbol: 'THROW_BASELINE', kind: 'count-map' },
    { path: DESIGN_RULE_GUARD_TEST_PATH, symbol: 'BASELINE', kind: 'count' }
  ]
}]);
const designRuleGuardTestSource = ({
  reconcile = ['normalizeCircularTrackSlots', 'changeCircularDepthSources'],
  writers = "{\n  'app/app-setup.js': { removeCircularDepthTrack: 'reconcile' }, // keeps every lane\n  'services/reset.js': { resetSettings: 'restore' }\n}",
  throwBaseline = "{ 'app/a.js': 3, 'app/b.js': 0 }",
  baseline = 46
} = {}) => [
  "import test from 'node:test';",
  '',
  '// Exports that may write outside the transition (R10 reconciles).',
  `const RECONCILE = new Set([${reconcile.map((name) => `'${name}'`).join(', ')}]);`,
  `const WRITERS = ${writers};`,
  `const THROW_BASELINE = ${throwBaseline};`,
  '',
  "test('baseline', () => {",
  `  const BASELINE = ${baseline};`,
  '  void RECONCILE; void WRITERS; void THROW_BASELINE; void BASELINE;',
  '});',
  ''
].join('\n');

const withDesignRuleGuardRepository = (runCase, {
  registry = DESIGN_RULE_REGISTRY,
  sources = { [DESIGN_RULE_GUARD_TEST_PATH]: designRuleGuardTestSource() }
} = {}) => (
  withChangeBudgetRepository(({ commit, execute, git, root, write }) => {
    write(DESIGN_RULE_GUARD_REGISTRY_PATH, registry);
    Object.entries(sources).forEach(([path, content]) => write(path, content));
    const base = commit('register design-rule guards');
    return runCase({ base, commit, execute, git, root, write });
  })
);
const runtimeChange = (write) => {
  write('gbdraw/web/js/services/session-file.js', 'export const readSession = () => ({ version: 2 });\n');
};

test('design-rule allowlists may contract alongside runtime files', () => {
  withDesignRuleGuardRepository(({ execute, write }) => {
    runtimeChange(write);
    write(DESIGN_RULE_GUARD_TEST_PATH, designRuleGuardTestSource({
      reconcile: ['normalizeCircularTrackSlots'],
      writers: "{\n  'services/reset.js': { resetSettings: 'restore' }\n}",
      throwBaseline: "{ 'app/a.js': 2, 'app/b.js': 0 }",
      baseline: 45
    }));
    const result = execute();
    assert.equal(result.status, 0, result.output);
    assert.match(result.output, /Gate: \*\*PASS\*\*/);
    assert.doesNotMatch(result.output, /design-rule\.co-change/);
    assert.match(result.output, /- contraction: R10 tests\/web\/design-rule-guard-fixture\.test\.mjs#RECONCILE/);
    assert.match(result.output, /- contraction: R10 tests\/web\/design-rule-guard-fixture\.test\.mjs#WRITERS/);
    assert.match(result.output, /- contraction: R10 tests\/web\/design-rule-guard-fixture\.test\.mjs#THROW_BASELINE/);
    assert.match(result.output, /- contraction: R10 tests\/web\/design-rule-guard-fixture\.test\.mjs#BASELINE/);
    assert.match(result.output, /Registry: 1 rule\(s\) from the trusted base/);
  });
});

test('design-rule allowlists cannot expand alongside runtime files', () => {
  const expansions = [
    ['RECONCILE', { reconcile: ['normalizeCircularTrackSlots', 'changeCircularDepthSources', 'applyTrackLayout'] }],
    ['WRITERS', { writers: "{\n  'app/app-setup.js': { removeCircularDepthTrack: 'reconcile' },\n  'services/reset.js': { resetSettings: 'restore' },\n  'app/feature-editor/placement-actions.js': { changeLayoutSetting: 'transition' }\n}" }],
    ['WRITERS', { writers: "{\n  'app/app-setup.js': { removeCircularDepthTrack: 'transition' },\n  'services/reset.js': { resetSettings: 'restore' }\n}" }],
    ['THROW_BASELINE', { throwBaseline: "{ 'app/a.js': 4, 'app/b.js': 0 }" }],
    ['THROW_BASELINE', { throwBaseline: "{ 'app/a.js': 3, 'app/b.js': 0, 'app/c.js': 1 }" }],
    // A new key is new even when Object.prototype has a member of that name.
    ['THROW_BASELINE', { throwBaseline: "{ 'app/a.js': 3, 'app/b.js': 0, constructor: 1 }" }],
    ['THROW_BASELINE', { throwBaseline: "{ 'app/a.js': 3, 'app/b.js': 0, '__proto__': 1 }" }],
    ['WRITERS', { writers: "{\n  'app/app-setup.js': { removeCircularDepthTrack: 'reconcile' },\n  'services/reset.js': { resetSettings: 'restore' },\n  constructor: {}\n}" }],
    ['BASELINE', { baseline: 47 }]
  ];
  expansions.forEach(([symbol, overrides]) => {
    withDesignRuleGuardRepository(({ execute, write }) => {
      runtimeChange(write);
      write(DESIGN_RULE_GUARD_TEST_PATH, designRuleGuardTestSource(overrides));
      for (const environment of [{}, { WEB_ARCHITECTURE_CHANGE: 'true' }]) {
        const result = execute({ environment });
        assert.equal(result.status, 1, `${symbol}\n${result.output}`);
        assert.match(result.output, /Gate: \*\*FAIL\*\*/);
        assert.match(
          result.output,
          new RegExp(`design-rule\\.co-change: a registered design-rule allowlist expanded alongside production runtime files \\([^)]*R10 tests/web/design-rule-guard-fixture\\.test\\.mjs#${symbol} \\(expansion\\)`)
        );
      }
    });
  });
});

test('a design-rule allowlist expansion without runtime files is an authority change', () => {
  withDesignRuleGuardRepository(({ execute, write }) => {
    write(DESIGN_RULE_GUARD_TEST_PATH, designRuleGuardTestSource({ baseline: 47 }));
    const result = execute();
    assert.equal(result.status, 0, result.output);
    assert.match(result.output, /Gate: \*\*PASS\*\*/);
    assert.match(result.output, /Review: \*\*REQUIRED\*\*/);
    assert.match(result.output, /Governance and authority: registered design-rule allowlists expanded \(R10 tests\/web\/design-rule-guard-fixture\.test\.mjs#BASELINE \(expansion\)\)/);
  });
});

test('an unreadable design-rule allowlist counts as an expansion', () => {
  withDesignRuleGuardRepository(({ execute, write }) => {
    runtimeChange(write);
    write(DESIGN_RULE_GUARD_TEST_PATH, designRuleGuardTestSource({ baseline: 'computeBaseline()' }));
    const result = execute();
    assert.equal(result.status, 1, result.output);
    assert.match(result.output, /design-rule\.co-change: .*#BASELINE \(unreadable\)/);
  });
});

test('the design-rule guard registry is read from the trusted base and fails closed when malformed', () => {
  withDesignRuleGuardRepository(({ execute, write }) => {
    // A head-side registry rewrite cannot relax the base registry.
    runtimeChange(write);
    write(DESIGN_RULE_GUARD_REGISTRY_PATH, designRuleGuardRegistrySource([]));
    write(DESIGN_RULE_GUARD_TEST_PATH, designRuleGuardTestSource({ baseline: 47 }));
    const result = execute();
    assert.equal(result.status, 1, result.output);
    assert.match(result.output, /design-rule\.co-change: .*#BASELINE \(expansion\)/);
    assert.match(result.output, /production runtime files and Web guard\/CI files changed together/);
  });
  withDesignRuleGuardRepository(({ execute }) => {
    const result = execute();
    assert.equal(result.status, 1, result.output);
    assert.match(result.output, /trusted design-rule guard registry: tools\/web-design-rule-guards\.json: rules\[0\]\.allowlists\[0\]\.kind: expected one of/);
  }, {
    registry: designRuleGuardRegistrySource([{
      id: 'R10',
      heading: 'R10',
      guards: [DESIGN_RULE_GUARD_TEST_PATH],
      allowlists: [{ path: DESIGN_RULE_GUARD_TEST_PATH, symbol: 'RECONCILE', kind: 'list' }]
    }])
  });
  withDesignRuleGuardRepository(({ execute }) => {
    const result = execute();
    assert.equal(result.status, 1, result.output);
    assert.match(result.output, /trusted design-rule guard registry: tools\/web-design-rule-guards\.json: invalid JSON/);
    assert.match(result.output, /Review: \*\*REQUIRED\*\*/);
  }, { registry: '{ not json\n' });
});

// A literal is read only when `;` closes it: `[...].concat(MORE)` is not the
// list between the brackets, on the same line or on the next one.
test('a design-rule allowlist literal that continues as an expression is unreadable', () => {
  const continued = [
    ['RECONCILE', (source) => source.replace("'changeCircularDepthSources']);", "'changeCircularDepthSources'].concat(MORE));")],
    ['RECONCILE', (source) => source.replace("'changeCircularDepthSources']);", "'changeCircularDepthSources']).add(MORE);")],
    ['THROW_BASELINE', (source) => source.replace("'app/b.js': 0 };", "'app/b.js': 0 }.merged;")],
    ['RECONCILE', (source) => source.replace("'changeCircularDepthSources']);", "'changeCircularDepthSources'])\n  .union(MORE);")],
    ['RECONCILE', (source) => source.replace("'changeCircularDepthSources']);", "'changeCircularDepthSources']\n  .concat(MORE));")],
    ['THROW_BASELINE', (source) => source.replace("'app/b.js': 0 };", "'app/b.js': 0 }\n  // merged below\n  .merged;")],
    ['THROW_BASELINE', (source) => source.replace("'app/b.js': 0 };", "'app/b.js': 0 }")]
  ];
  continued.forEach(([symbol, rewrite]) => {
    withDesignRuleGuardRepository(({ execute, write }) => {
      runtimeChange(write);
      const source = rewrite(designRuleGuardTestSource());
      assert.notEqual(source, designRuleGuardTestSource());
      write(DESIGN_RULE_GUARD_TEST_PATH, source);
      const result = execute();
      assert.equal(result.status, 1, `${symbol}\n${result.output}`);
      assert.match(result.output, new RegExp(`design-rule\\.co-change: .*#${symbol} \\(unreadable\\)`));
    });
  });
});

// Registry schema 2: a guard records its root cause, its CI job and when it
// can retire; a rule outside gbdraw/web/CLAUDE.md names its own rulesPath; and
// an allowlist may be a frozen JS literal, a path -> names map, the entries of
// a JSON file, or a module-level Python literal.
const DESIGN_RULE_PYTHON_GUARD_PATH = 'tests/test_design_rule_guard_fixture.py';
const DESIGN_RULE_JSON_GUARD_PATH = 'tests/web/contracts/design-rule-guard-fixture.json';
const designRuleGuard = (path, overrides = {}) => ({
  path, rootCause: 'P3', job: 'web-contracts-pr', retireWhen: 'never', ...overrides
});
const designRuleRules2 = () => [{
  id: 'R10',
  heading: 'R10: A watcher does not repair state',
  guards: [
    designRuleGuard(DESIGN_RULE_GUARD_TEST_PATH),
    designRuleGuard(DESIGN_RULE_JSON_GUARD_PATH, {
      rootCause: 'new', job: 'promotion', retireWhen: 'when live edits and Generate share one projection'
    })
  ],
  allowlists: [
    { path: DESIGN_RULE_GUARD_TEST_PATH, symbol: 'RULE_WRITERS', kind: 'set-map' },
    { path: DESIGN_RULE_GUARD_TEST_PATH, symbol: 'FROZEN', kind: 'set' },
    { path: DESIGN_RULE_JSON_GUARD_PATH, symbol: 'entries', kind: 'json-entries' }
  ]
}, {
  id: 'TYPE-DEBT',
  heading: 'Type checking',
  rulesPath: 'CLAUDE.md',
  guards: [designRuleGuard(DESIGN_RULE_PYTHON_GUARD_PATH, { rootCause: 'P10', job: 'core-pr' })],
  allowlists: [
    { path: DESIGN_RULE_PYTHON_GUARD_PATH, symbol: 'TYPE_DEBT_BASELINE', kind: 'py-count-map' },
    { path: DESIGN_RULE_PYTHON_GUARD_PATH, symbol: 'UNPROVEN_BASELINE', kind: 'py-count-map' },
    { path: DESIGN_RULE_PYTHON_GUARD_PATH, symbol: 'ALLOWED_FILES', kind: 'py-set' },
    { path: DESIGN_RULE_PYTHON_GUARD_PATH, symbol: 'CAP', kind: 'py-count' }
  ]
}];
const DESIGN_RULE_SUBJECTS_2 = Object.freeze({
  RULE_WRITERS: `R10 ${DESIGN_RULE_GUARD_TEST_PATH}#RULE_WRITERS`,
  FROZEN: `R10 ${DESIGN_RULE_GUARD_TEST_PATH}#FROZEN`,
  entries: `R10 ${DESIGN_RULE_JSON_GUARD_PATH}#entries`,
  TYPE_DEBT_BASELINE: `TYPE-DEBT ${DESIGN_RULE_PYTHON_GUARD_PATH}#TYPE_DEBT_BASELINE`,
  UNPROVEN_BASELINE: `TYPE-DEBT ${DESIGN_RULE_PYTHON_GUARD_PATH}#UNPROVEN_BASELINE`,
  ALLOWED_FILES: `TYPE-DEBT ${DESIGN_RULE_PYTHON_GUARD_PATH}#ALLOWED_FILES`,
  CAP: `TYPE-DEBT ${DESIGN_RULE_PYTHON_GUARD_PATH}#CAP`
});
const designRuleGuardSources2 = ({
  ruleWriters = "Object.freeze({\n  // The edit transition and the popup's rule scopes.\n  'js/app/visibility-actions.js': ['editRules'],\n  'js/services/config.js': ['resetBaseline', 'applyDrawingData']\n})",
  frozen = "Object.freeze(new Set(['normalizeSlots', 'changeSources']))",
  entries = [
    { id: 'label-placement-awaits-reflow', kind: 'leader', condition: { autoReflow: false } },
    { id: 'legend-order-awaits-generate', kind: 'legend', condition: { autoReflow: true } }
  ],
  typeDebt = '{"gbdraw/api/options.py": 2, "gbdraw/render/export.py": 1}',
  unproven = '{\n    "gbdraw/analysis/collinearity.py": (10, _COMPARISON),  # (sites, what they check)\n    "gbdraw/api/options.py": (3, _TYPED_API),\n}',
  allowedFiles = 'frozenset({\n    "gbdraw/web/js/a.js",\n    "gbdraw/web/js/b.js",\n})',
  cap = '12'
} = {}) => ({
  [DESIGN_RULE_GUARD_TEST_PATH]: [
    "import test from 'node:test';",
    '',
    `export const RULE_WRITERS = ${ruleWriters};`,
    `const FROZEN = ${frozen};`,
    '',
    "test('writers', () => { void RULE_WRITERS; void FROZEN; });",
    ''
  ].join('\n'),
  [DESIGN_RULE_JSON_GUARD_PATH]: typeof entries === 'string'
    ? entries
    : `${JSON.stringify({ description: 'Differences a decision allows.', entries }, null, 2)}\n`,
  [DESIGN_RULE_PYTHON_GUARD_PATH]: [
    '"""Shrink-only baselines.',
    '',
    'CAP = 99 in this docstring is not the assignment.',
    '"""',
    'from __future__ import annotations',
    '',
    '_COMPARISON = "comparison checks"',
    '_TYPED_API = "typed API checks"',
    '# TYPE_DEBT_BASELINE = {} in a comment is not the assignment either.',
    `TYPE_DEBT_BASELINE: dict[str, int] = ${typeDebt}`,
    '',
    '# Shrink-only: path -> (unproven sites, what they check).',
    `UNPROVEN_BASELINE: dict[str, tuple[int, str]] = ${unproven}`,
    `ALLOWED_FILES = ${allowedFiles}`,
    `CAP = ${cap}  # the cap`,
    '',
    '',
    'def test_baselines() -> None:',
    '    assert CAP and TYPE_DEBT_BASELINE is not None and UNPROVEN_BASELINE and ALLOWED_FILES',
    '    for path in UNPROVEN_BASELINE:',
    '        assert UNPROVEN_BASELINE.get(path, (0, ""))[0] <= len(UNPROVEN_BASELINE) * CAP',
    '    counts = {path: count for path, (count, _note) in UNPROVEN_BASELINE.items()}',
    '    assert counts == {path: UNPROVEN_BASELINE[path][0] for path in sorted(UNPROVEN_BASELINE)}, (',
    '        f"lower UNPROVEN_BASELINE = {counts}"',
    '    )',
    '    assert check(UNPROVEN_BASELINE, strict=True) and check(baseline=UNPROVEN_BASELINE)',
    '    assert "gbdraw/web/js/c.js" not in ALLOWED_FILES and ALLOWED_FILES != set()',
    ''
  ].join('\n')
});
const withDesignRuleGuardRepository2 = (runCase, { rules = designRuleRules2(), schemaVersion = 2 } = {}) => (
  withDesignRuleGuardRepository(({ execute, root, write }) => runCase({
    execute,
    root,
    write,
    writeGuards: (overrides) => Object.entries(designRuleGuardSources2(overrides))
      .forEach(([path, source]) => write(path, source))
  }), { registry: designRuleGuardRegistrySource(rules, schemaVersion), sources: designRuleGuardSources2() })
);
const pythonProductionChange = (write) => {
  write('gbdraw/api/options.py', 'LIMIT = 2\n');
};
const DESIGN_RULE_CONTRACTIONS_2 = Object.freeze({
  ruleWriters: "Object.freeze({\n  'js/services/config.js': ['applyDrawingData']\n})",
  frozen: "Object.freeze(new Set(['normalizeSlots']))",
  entries: [{ id: 'label-placement-awaits-reflow', kind: 'leader', condition: { autoReflow: false } }],
  typeDebt: '{"gbdraw/api/options.py": 1}',
  unproven: '{"gbdraw/analysis/collinearity.py": (9, _COMPARISON), "gbdraw/api/options.py": (3, _TYPED_API)}',
  allowedFiles: 'frozenset({"gbdraw/web/js/a.js"})',
  cap: '11'
});
// [symbol, source overrides, direction when it is not a readable expansion]
const DESIGN_RULE_EXPANSIONS_2 = Object.freeze([
  ['RULE_WRITERS', { ruleWriters: "Object.freeze({\n  'js/app/visibility-actions.js': ['editRules'],\n  'js/services/config.js': ['resetBaseline', 'applyDrawingData'],\n  'js/services/reset.js': ['resetDraft']\n})" }],
  ['RULE_WRITERS', { ruleWriters: "Object.freeze({\n  'js/app/visibility-actions.js': ['editRules', 'replaceRules'],\n  'js/services/config.js': ['resetBaseline', 'applyDrawingData']\n})" }],
  ['RULE_WRITERS', { ruleWriters: "Object.freeze({ ...OTHER_WRITERS, 'js/services/config.js': ['resetBaseline'] })" }, 'unreadable'],
  ['RULE_WRITERS', { ruleWriters: "Object.freeze({ 'js/services/config.js': 'resetBaseline' })" }, 'unreadable'],
  ['FROZEN', { frozen: "Object.freeze(new Set(['normalizeSlots', 'changeSources', 'applyLayout']))" }],
  ['FROZEN', { frozen: "Object.freeze(new Set(['normalizeSlots']).add(MORE))" }, 'unreadable'],
  ['FROZEN', { frozen: "Object.freeze(['normalizeSlots', 'changeSources'])\n  .concat(MORE)" }, 'unreadable'],
  ['entries', { entries: [
    { id: 'label-placement-awaits-reflow', kind: 'leader', condition: { autoReflow: false } },
    { id: 'legend-order-awaits-generate', kind: 'legend', condition: { autoReflow: true } },
    { id: 'track-order-awaits-generate', kind: 'track', condition: {} }
  ] }],
  ['entries', { entries: [
    { id: 'label-placement-awaits-reflow', kind: 'leader', condition: {} },
    { id: 'legend-order-awaits-generate', kind: 'legend', condition: { autoReflow: true } }
  ] }],
  ['entries', { entries: [
    { id: 'legend-order-awaits-generate', kind: 'legend', condition: { autoReflow: true } },
    { id: 'legend-order-awaits-generate', kind: 'leader', condition: {} }
  ] }, 'unreadable'],
  ['entries', { entries: [{ kind: 'leader', condition: {} }] }, 'unreadable'],
  ['entries', { entries: '{ "entries": [\n' }, 'unreadable'],
  ['TYPE_DEBT_BASELINE', { typeDebt: '{"gbdraw/api/options.py": 3, "gbdraw/render/export.py": 1}' }],
  ['TYPE_DEBT_BASELINE', { typeDebt: '{"gbdraw/api/options.py": 2, "gbdraw/render/export.py": 1, "gbdraw/cli.py": 1}' }],
  ['TYPE_DEBT_BASELINE', { typeDebt: '{**_OTHER, "gbdraw/api/options.py": 2}' }, 'unreadable'],
  ['TYPE_DEBT_BASELINE', { typeDebt: 'dict(_OTHER)' }, 'unreadable'],
  ['TYPE_DEBT_BASELINE', { typeDebt: '{"gbdraw/api/options.py": 2, "gbdraw/render/export.py": 1} | _OTHER' }, 'unreadable'],
  ['TYPE_DEBT_BASELINE', { typeDebt: '{"gbdraw/api/options.py": 2, "gbdraw/render/export.py": 1}\nTYPE_DEBT_BASELINE["gbdraw/cli.py"] = 4' }, 'unreadable'],
  ['TYPE_DEBT_BASELINE', { typeDebt: '{"gbdraw/api/options.py": -1}' }, 'unreadable'],
  ['UNPROVEN_BASELINE', { unproven: '{"gbdraw/analysis/collinearity.py": (11, _COMPARISON), "gbdraw/api/options.py": (3, _TYPED_API)}' }],
  ['UNPROVEN_BASELINE', { unproven: '{"gbdraw/analysis/collinearity.py": (_COMPARISON, 10), "gbdraw/api/options.py": (3, _TYPED_API)}' }, 'unreadable'],
  ['ALLOWED_FILES', { allowedFiles: '{"gbdraw/web/js/a.js", "gbdraw/web/js/b.js", "gbdraw/web/js/c.js"}' }],
  ['ALLOWED_FILES', { allowedFiles: 'frozenset(_GENERATED | {"gbdraw/web/js/a.js"})' }, 'unreadable'],
  ['ALLOWED_FILES', { allowedFiles: '("gbdraw/web/js/a.js", "gbdraw/web/js/b.js" "gbdraw/web/js/c.js")' }, 'unreadable'],
  ['CAP', { cap: '13' }],
  ['CAP', { cap: 'len(_ALL)' }, 'unreadable'],
  ['CAP', { cap: '12 + _MORE' }, 'unreadable']
]);
const designRuleSubjectPattern = (symbol, direction = 'expansion') => (
  `${DESIGN_RULE_SUBJECTS_2[symbol].replace(/[.*+?^${}()|[\]\\]/g, '\\$&')} \\(${direction}\\)`
);

test('schema 2 design-rule allowlists of every kind may contract alongside runtime files', () => {
  withDesignRuleGuardRepository2(({ execute, write, writeGuards }) => {
    runtimeChange(write);
    pythonProductionChange(write);
    writeGuards(DESIGN_RULE_CONTRACTIONS_2);
    const result = execute();
    assert.equal(result.status, 0, result.output);
    assert.match(result.output, /Gate: \*\*PASS\*\*/);
    assert.doesNotMatch(result.output, /design-rule\.co-change|- expansion: /);
    Object.values(DESIGN_RULE_SUBJECTS_2).forEach((subject) => {
      assert.ok(result.output.includes(`- contraction: ${subject}\n`), `${subject}\n${result.output}`);
    });
    assert.match(result.output, /Registry: 2 rule\(s\) from the trusted base/);
  });
  // The same values written another way are neither a contraction nor an expansion.
  withDesignRuleGuardRepository2(({ execute, write, writeGuards }) => {
    runtimeChange(write);
    pythonProductionChange(write);
    writeGuards({
      ruleWriters: "{ 'js/services/config.js': ['applyDrawingData', 'resetBaseline'], 'js/app/visibility-actions.js': ['editRules'], }",
      frozen: "new Set(['changeSources', 'normalizeSlots'])",
      typeDebt: "{\n    'gbdraw/render/export.py': 1,  # export\n    'gbdraw/api/options.py': 2,\n}",
      unproven: '{"gbdraw/api/options.py": (3, "typed API checks"), "gbdraw/analysis/collinearity.py": (10, _COMPARISON)}',
      allowedFiles: '["gbdraw/web/js/b.js", "gbdraw/web/js/a.js"]',
      cap: '12'
    });
    const result = execute();
    assert.equal(result.status, 0, result.output);
    assert.doesNotMatch(result.output, /design-rule\.co-change|- expansion: |- contraction: /);
  });
});

test('schema 2 design-rule allowlists cannot expand alongside runtime files', () => {
  DESIGN_RULE_EXPANSIONS_2.forEach(([symbol, overrides, direction]) => {
    withDesignRuleGuardRepository2(({ execute, write, writeGuards }) => {
      runtimeChange(write);
      writeGuards(overrides);
      const result = execute();
      const label = `${symbol} ${JSON.stringify(overrides)}\n${result.output}`;
      assert.equal(result.status, 1, label);
      assert.match(result.output, new RegExp(
        'design-rule\\.co-change: a registered design-rule allowlist expanded alongside production runtime files '
        + `\\(${designRuleSubjectPattern(symbol, direction)}\\)`
      ), label);
    });
  });
});

test('a schema 2 design-rule allowlist expansion without runtime files is an authority change', () => {
  DESIGN_RULE_EXPANSIONS_2.filter(([, , direction]) => !direction).forEach(([symbol, overrides]) => {
    withDesignRuleGuardRepository2(({ execute, writeGuards }) => {
      writeGuards(overrides);
      const result = execute();
      const label = `${symbol} ${JSON.stringify(overrides)}\n${result.output}`;
      assert.equal(result.status, 0, label);
      assert.match(result.output, /Gate: \*\*PASS\*\*/, label);
      assert.match(result.output, /Review: \*\*REQUIRED\*\*/, label);
      assert.match(result.output, new RegExp(
        `Governance and authority: registered design-rule allowlists expanded \\(${designRuleSubjectPattern(symbol)}\\)`
      ), label);
    });
  });
});

// Counts compare per key with <=: lowering a count or removing a key needs
// nothing else, with or without a runtime file in the diff.
test('a design-rule count contraction without runtime files passes and is not an authority change', () => {
  withDesignRuleGuardRepository2(({ execute, writeGuards }) => {
    const { typeDebt, unproven, cap } = DESIGN_RULE_CONTRACTIONS_2;
    writeGuards({ typeDebt, unproven, cap });
    const result = execute();
    assert.equal(result.status, 0, result.output);
    assert.match(result.output, /Gate: \*\*PASS\*\*/);
    assert.doesNotMatch(result.output, /design-rule\.co-change|allowlists expanded|- expansion: /);
    ['TYPE_DEBT_BASELINE', 'UNPROVEN_BASELINE', 'CAP'].forEach((symbol) => {
      assert.ok(result.output.includes(`- contraction: ${DESIGN_RULE_SUBJECTS_2[symbol]}\n`), result.output);
    });
  });
  withDesignRuleGuardRepository(({ execute, write }) => {
    write(DESIGN_RULE_GUARD_TEST_PATH, designRuleGuardTestSource({ throwBaseline: "{ 'app/a.js': 2 }", baseline: 45 }));
    const result = execute();
    assert.equal(result.status, 0, result.output);
    assert.doesNotMatch(result.output, /design-rule\.co-change|allowlists expanded|- expansion: /);
    assert.match(result.output, /- contraction: R10 tests\/web\/design-rule-guard-fixture\.test\.mjs#THROW_BASELINE/);
    assert.match(result.output, /- contraction: R10 tests\/web\/design-rule-guard-fixture\.test\.mjs#BASELINE/);
  });
});

// A Python baseline guards the Python package: for a py-* allowlist a changed
// gbdraw/**/*.py outside gbdraw/web/ is production too.
test('a Python design-rule baseline cannot expand alongside Python package files', () => {
  const pythonExpansion = { typeDebt: '{"gbdraw/api/options.py": 3, "gbdraw/render/export.py": 1}' };
  const webExpansion = DESIGN_RULE_EXPANSIONS_2[0][1];
  const blocked = /design-rule\.co-change: a registered design-rule allowlist expanded alongside production runtime files /;
  withDesignRuleGuardRepository2(({ execute, write, writeGuards }) => {
    pythonProductionChange(write);
    writeGuards({ ...pythonExpansion, ...webExpansion });
    const result = execute();
    assert.equal(result.status, 1, result.output);
    assert.match(result.output, /Gate: \*\*FAIL\*\*/);
    // Only the Python baseline blocks; the Web allowlist stays an authority change.
    assert.match(result.output, new RegExp(`${blocked.source}\\(${designRuleSubjectPattern('TYPE_DEBT_BASELINE')}\\)`));
    assert.match(result.output, new RegExp(`- expansion: ${designRuleSubjectPattern('RULE_WRITERS')}`));
  });
  // Not Python production: the Web tree, tests, tools, and a non-.py package file.
  ['gbdraw/web/server.py', 'tests/test_other.py', 'tools/prepare.py', 'gbdraw/data/config.toml'].forEach((path) => {
    withDesignRuleGuardRepository2(({ execute, write, writeGuards }) => {
      write(path, 'LIMIT = 2\n');
      writeGuards(pythonExpansion);
      const result = execute();
      assert.equal(result.status, 0, `${path}\n${result.output}`);
      assert.doesNotMatch(result.output, blocked, path);
      assert.match(result.output, new RegExp(
        `registered design-rule allowlists expanded \\(${designRuleSubjectPattern('TYPE_DEBT_BASELINE')}\\)`
      ), path);
    });
  });
  // Python package files do not make a Web allowlist expansion blocking, in either schema.
  withDesignRuleGuardRepository2(({ execute, write, writeGuards }) => {
    pythonProductionChange(write);
    writeGuards(webExpansion);
    const result = execute();
    assert.equal(result.status, 0, result.output);
    assert.doesNotMatch(result.output, blocked);
    assert.match(result.output, /Review: \*\*REQUIRED\*\*/);
  });
  withDesignRuleGuardRepository(({ execute, write }) => {
    pythonProductionChange(write);
    write(DESIGN_RULE_GUARD_TEST_PATH, designRuleGuardTestSource({ baseline: 47 }));
    const result = execute();
    assert.equal(result.status, 0, result.output);
    assert.doesNotMatch(result.output, blocked);
  });
  // Python has no `const`: a baseline is read only when one module-level
  // assignment is the single place that binds or changes the name. A removed
  // baseline contracts, but alongside production code its test has lost a
  // registered literal (design-rule.count-scope).
  const pythonSource = designRuleGuardSources2()[DESIGN_RULE_PYTHON_GUARD_PATH];
  const declaration = 'TYPE_DEBT_BASELINE: dict[str, int] = ';
  const unreadable = `${blocked.source}\\(${designRuleSubjectPattern('TYPE_DEBT_BASELINE', 'unreadable')}\\)`;
  const withStatement = (statement) => pythonSource.replace('\n\ndef test_baselines', `\n${statement}\n\n\ndef test_baselines`);
  [
    [pythonSource.replace(declaration, '_, TYPE_DEBT_BASELINE = 0, ')],
    [pythonSource.replace(declaration, `if True:\n    ${declaration}`)],
    // The module-level literal stays, and the name is bound or changed elsewhere.
    [withStatement('if sys.version_info < (3, 11):\n    TYPE_DEBT_BASELINE["gbdraw/cli.py"] = 3')],
    [withStatement('try:\n    import tomllib\nexcept ImportError:\n    TYPE_DEBT_BASELINE = {"gbdraw/cli.py": 3}')],
    [withStatement('for _path in ("gbdraw/cli.py",):\n    TYPE_DEBT_BASELINE.setdefault(_path, 3)')],
    [withStatement('if True: TYPE_DEBT_BASELINE.update({"gbdraw/cli.py": 3})')],
    [withStatement('_X = 1; TYPE_DEBT_BASELINE = {"gbdraw/cli.py": 3}')],
    [withStatement('_X, TYPE_DEBT_BASELINE = 1, {"gbdraw/cli.py": 3}')],
    [withStatement('_A, TYPE_DEBT_BASELINE, _B = 1, {"gbdraw/cli.py": 3}, 2')],
    [withStatement('_X = TYPE_DEBT_BASELINE = {"gbdraw/cli.py": 3}')],
    [withStatement('_ALL |= TYPE_DEBT_BASELINE; TYPE_DEBT_BASELINE |= {"gbdraw/cli.py": 3}')],
    [withStatement('from tests.other_baselines import TYPE_DEBT_BASELINE')],
    [withStatement('import tests.other_baselines as TYPE_DEBT_BASELINE')],
    [withStatement('def _grow():\n    global TYPE_DEBT_BASELINE\n    TYPE_DEBT_BASELINE = {"gbdraw/cli.py": 3}')],
    [pythonSource.replace(/^TYPE_DEBT_BASELINE.*\n/m, '').replace('TYPE_DEBT_BASELINE is not None and ', ''), 1,
      '^(?=[\\s\\S]*- contraction: TYPE-DEBT [^\\n]*#TYPE_DEBT_BASELINE\\n)(?=[\\s\\S]*design-rule\\.count-scope: [^\\n]*'
      + 'tests/test_design_rule_guard_fixture\\.py \\(registered literal missing or unreadable\\))']
  ].forEach(([source, status = 1, expected = unreadable]) => {
    withDesignRuleGuardRepository2(({ execute, write }) => {
      assert.notEqual(source, pythonSource);
      pythonProductionChange(write);
      write(DESIGN_RULE_PYTHON_GUARD_PATH, source);
      const result = execute();
      assert.equal(result.status, status, `${source}\n${result.output}`);
      assert.match(result.output, new RegExp(expected), source);
    });
  });
});

test('the design-rule guard registry parser validates each schema 2 field and keeps schema 1 strict', () => {
  const malformed = (mutate, expected, schemaVersion = 2) => {
    const rules = designRuleRules2();
    mutate(rules);
    withDesignRuleGuardRepository2(({ execute }) => {
      const result = execute();
      const label = `${expected}\n${result.output}`;
      assert.equal(result.status, 1, label);
      assert.ok(result.output.includes(
        `trusted design-rule guard registry: tools/web-design-rule-guards.json: ${expected}`
      ), label);
    }, { rules, schemaVersion });
  };
  malformed(() => {}, 'schemaVersion must be 1 or 2', 3);
  malformed((rules) => { rules[0].guards[0] = DESIGN_RULE_GUARD_TEST_PATH; }, 'rules[0].guards[0]: expected an object');
  malformed((rules) => { delete rules[0].guards[0].path; }, 'rules[0].guards[0].path: expected a path');
  malformed((rules) => { rules[0].guards[0].rootCause = 'P11'; }, 'rules[0].guards[0].rootCause: expected P1..P10 or new');
  malformed((rules) => { delete rules[0].guards[1].rootCause; }, 'rules[0].guards[1].rootCause: expected P1..P10 or new');
  malformed((rules) => { rules[0].guards[0].job = 'Web contracts'; }, 'rules[0].guards[0].job: expected a CI job id');
  malformed((rules) => { delete rules[1].guards[0].job; }, 'rules[1].guards[0].job: expected a CI job id');
  malformed((rules) => { rules[0].guards[0].retireWhen = ' '; }, 'rules[0].guards[0].retireWhen: expected text');
  malformed((rules) => { delete rules[0].guards[0].retireWhen; }, 'rules[0].guards[0].retireWhen: expected text');
  malformed((rules) => { rules[0].guards = DESIGN_RULE_GUARD_TEST_PATH; }, 'rules[0].guards: expected an array');
  malformed((rules) => { delete rules[1].rulesPath; }, 'rules[1].id: expected R<n>, or an upper-case id with its own rulesPath');
  malformed((rules) => { rules[1].id = 'type-debt'; }, 'rules[1].id: expected R<n>, or an upper-case id with its own rulesPath');
  malformed((rules) => { rules[1].rulesPath = ''; }, 'rules[1].rulesPath: expected a path');
  malformed((rules) => { rules[1].id = 'R10'; }, 'rules[1].id: duplicate R10');
  malformed((rules) => { rules[1].allowlists[0].kind = 'py-map'; }, 'rules[1].allowlists[0].kind: expected one of count, count-map, set, writer-map, set-map, json-entries, json-ids, py-count, py-count-map, py-set');
  // Schema 1 reads as before: bare path guards, R<n> ids, and its four kinds.
  malformed((rules) => { rules.pop(); }, 'rules[0].guards: expected an array of paths', 1);
  malformed((rules) => {
    rules[0].guards = [DESIGN_RULE_GUARD_TEST_PATH];
    rules[0].allowlists = [];
    rules[1].guards = [DESIGN_RULE_PYTHON_GUARD_PATH];
    rules[1].allowlists = [];
  }, 'rules[1].id: expected R<n>', 1);
  malformed((rules) => {
    rules.pop();
    rules[0].guards = [DESIGN_RULE_GUARD_TEST_PATH];
  }, 'rules[0].allowlists[0].kind: expected one of count, count-map, set, writer-map', 1);
});

// Every rule a registry names exists in its rules file and every guard exists.
const assertDesignRuleRegistryResolves = (registry, root) => {
  assert.ok([1, 2].includes(registry.schemaVersion), 'schemaVersion');
  registry.rules.forEach((rule) => {
    // The numbered Web rules are `###` headings; a rule with its own rulesPath may sit at any level.
    const level = rule.rulesPath ? /^#{1,6} / : /^### /;
    const headings = readFileSync(join(root, rule.rulesPath || registry.rulesPath), 'utf8').split(/\r?\n/)
      .filter((line) => level.test(line)).map((line) => line.replace(/^#+ /, ''));
    assert.ok(headings.includes(rule.heading), `${rule.id}: ${rule.heading}`);
    rule.guards.forEach((guard) => {
      const path = typeof guard === 'string' ? guard : guard.path;
      assert.ok(existsSync(join(root, path)), `${rule.id}: ${path}`);
    });
  });
};
// The Gate reads the registry from the trusted base and fails closed: a
// malformed registry fails every later pull request, and a literal the Gate
// cannot read fails every later pull request that changes its file together
// with production code. So a candidate registry goes through the Gate itself
// before it merges. The registry and every registered literal's file are the
// base, and each file then gains one blank line, which makes the Gate read
// every literal with its registered kind on both sides.
const assertGateReadsDesignRuleRegistry = (registrySource, root) => {
  const paths = [...new Set(
    JSON.parse(registrySource).rules.flatMap((rule) => rule.allowlists.map(({ path }) => path))
  )];
  const sources = Object.fromEntries(paths.map((path) => [path, readFileSync(join(root, path), 'utf8')]));
  withDesignRuleGuardRepository(({ execute, write }) => {
    paths.forEach((path) => write(path, `${sources[path]}\n`));
    const { output } = execute();
    assert.doesNotMatch(output, /trusted design-rule guard registry:/);
    assert.match(output, /Registry: \d+ rule\(s\) from the trusted base/);
    // A literal that is missing or unreadable shows as an expansion; a readable one is unchanged.
    assert.deepEqual(output.split('\n').filter((line) => /^- (?:expansion|contraction): /.test(line)), []);
  }, { registry: registrySource, sources });
};

test('the head-registry check resolves schema 2 guards, a rule\'s own rulesPath, and Python and JSON literals', () => {
  withDesignRuleGuardRepository2(({ root, write, writeGuards }) => {
    write(WEB_DESIGN_RULES_PATH, '## Design rules\n\n### R10: A watcher does not repair state\n');
    write('CLAUDE.md', '# CLAUDE.md\n\n## Coding Conventions\n\n#### Type checking\n');
    const registry = JSON.parse(designRuleGuardRegistrySource(designRuleRules2(), 2));
    assertDesignRuleRegistryResolves(registry, root);
    assertGateReadsDesignRuleRegistry(JSON.stringify(registry), root);
    const broken = (mutate, expected, check = assertDesignRuleRegistryResolves) => {
      const candidate = structuredClone(registry);
      mutate(candidate);
      assert.throws(() => check(check === assertDesignRuleRegistryResolves ? candidate : JSON.stringify(candidate), root), expected);
    };
    broken((candidate) => { candidate.schemaVersion = 3; }, /schemaVersion/);
    broken((candidate) => { candidate.rules[1].heading = 'Type debt'; }, /TYPE-DEBT: Type debt/);
    broken((candidate) => { delete candidate.rules[1].rulesPath; }, /TYPE-DEBT: Type checking/);
    broken((candidate) => { candidate.rules[1].guards[0].path = 'tests/test_missing.py'; }, /tests\/test_missing\.py/);
    // A registered literal the Gate does not find, or finds and cannot read.
    const unreadable = (symbol, status) => new RegExp(`#${symbol} \\(unreadable in the trusted base \\(${status}\\)\\)`);
    const gate = assertGateReadsDesignRuleRegistry;
    broken((candidate) => { candidate.schemaVersion = 3; }, /trusted design-rule guard registry:/, gate);
    broken((candidate) => { candidate.rules[1].allowlists[0].symbol = 'DEBT_BASELINE'; }, unreadable('DEBT_BASELINE', 'missing'), gate);
    broken((candidate) => { candidate.rules[0].allowlists[0].symbol = 'WRITERS'; }, unreadable('WRITERS', 'missing'), gate);
    broken((candidate) => { candidate.rules[0].allowlists[2].symbol = 'allowed'; }, unreadable('allowed', 'missing'), gate);
    // A literal that is read but is not of its registered kind.
    broken((candidate) => { candidate.rules[1].allowlists[3].kind = 'py-set'; }, /#CAP \(unreadable\)/, gate);
    broken((candidate) => { candidate.rules[1].allowlists[0].path = 'tests/test_missing.py'; }, /ENOENT/, gate);
    writeGuards({ typeDebt: 'dict(_OTHER)', frozen: "Object.freeze(['normalizeSlots'])\n  .concat(MORE)" });
    assert.throws(() => gate(JSON.stringify(registry), root), unreadable('TYPE_DEBT_BASELINE', 'unreadable'));
    assert.throws(() => gate(JSON.stringify(registry), root), unreadable('FROZEN', 'unreadable'));
  });
});

// design-rule.count-scope (BA-OD40 (d), form B): a test that holds a count
// baseline counts with its own code, so with production code it may change
// only inside its registered literals.
const countScope = (detail) => new RegExp(
  'design-rule\\.count-scope: the test or counting code of a registered count baseline changed '
  + `alongside production runtime files \\([^\\n]*${detail}`
);
const escaped = (path) => path.replace(/[.*+?^${}()|[\]\\]/g, '\\$&');

test('a count-baseline test changes only inside its registered literals alongside runtime files', () => {
  // The counting code narrows (here: one more line of test code).
  const narrowed = `${designRuleGuardTestSource()}const COUNTED = (path) => path.startsWith('app/');\nvoid COUNTED;\n`;
  withDesignRuleGuardRepository(({ execute, write }) => {
    runtimeChange(write);
    write(DESIGN_RULE_GUARD_TEST_PATH, narrowed);
    const result = execute();
    assert.equal(result.status, 1, result.output);
    assert.match(result.output, countScope(escaped(DESIGN_RULE_GUARD_TEST_PATH)));
  });
  withDesignRuleGuardRepository(({ execute, root, write }) => {
    runtimeChange(write);
    rmSync(join(root, DESIGN_RULE_GUARD_TEST_PATH));
    const result = execute();
    assert.equal(result.status, 1, result.output);
    assert.match(result.output, countScope(`${escaped(DESIGN_RULE_GUARD_TEST_PATH)} \\(deleted\\)`));
  });
  // Without production code the same change is a test change; inside the
  // literals (a lowered count) it passes with production code.
  withDesignRuleGuardRepository(({ execute, write }) => {
    write(DESIGN_RULE_GUARD_TEST_PATH, narrowed);
    const result = execute();
    assert.equal(result.status, 0, result.output);
    assert.doesNotMatch(result.output, /design-rule\.count-scope/);
  });
  withDesignRuleGuardRepository(({ execute, write }) => {
    runtimeChange(write);
    write(DESIGN_RULE_GUARD_TEST_PATH, designRuleGuardTestSource({ baseline: 45, throwBaseline: "{ 'app/a.js': 2 }" }));
    const result = execute();
    assert.equal(result.status, 0, result.output);
    assert.doesNotMatch(result.output, /design-rule\.count-scope/);
  });
  // A Python count baseline: a package change counts as production; a test
  // with set literals only is not a count-baseline test.
  withDesignRuleGuardRepository2(({ execute, write }) => {
    pythonProductionChange(write);
    const source = designRuleGuardSources2();
    write(DESIGN_RULE_PYTHON_GUARD_PATH, source[DESIGN_RULE_PYTHON_GUARD_PATH].replace(
      '    for path in UNPROVEN_BASELINE:',
      '    for path in [p for p in UNPROVEN_BASELINE if p.startswith("gbdraw/api/")]:'
    ));
    write(DESIGN_RULE_GUARD_TEST_PATH, `${source[DESIGN_RULE_GUARD_TEST_PATH]}void 0;\n`);
    const result = execute();
    assert.equal(result.status, 1, result.output);
    assert.match(result.output, countScope(escaped(DESIGN_RULE_PYTHON_GUARD_PATH)));
    assert.doesNotMatch(result.output, /count-scope[^\n]*design-rule-guard-fixture\.test\.mjs/);
  });
});

// F1: a registered literal that the head no longer holds as a readable
// literal (a `let`, an indirection, a rename everywhere, a move into another
// file) leaves the counting code free, so it is a count-scope change.
test('a count-baseline test keeps its registered literals readable alongside runtime files', () => {
  const missing = (path) => countScope(`${escaped(path)} \\(registered literal missing or unreadable\\)`);
  const literal = '  const BASELINE = 46;';
  const cases = {
    let: (source) => source.replace(literal, '  let BASELINE = 99;'),
    indirection: (source) => source.replace(literal, '  const BASELINE_NEXT = 99;\n  const BASELINE = BASELINE_NEXT;'),
    renamed: (source) => source.replaceAll('BASELINE', 'CEILING').replace('const CEILING = 46;', 'const CEILING = 99;'),
    moved: (source) => `import { BASELINE } from './count-literal.mjs';\n${source.replace(`${literal}\n`, '')}`
  };
  Object.entries(cases).forEach(([name, edit]) => withDesignRuleGuardRepository(({ execute, write }) => {
    runtimeChange(write);
    write(DESIGN_RULE_GUARD_TEST_PATH, edit(designRuleGuardTestSource()));
    if (name === 'moved') write('tests/web/count-literal.mjs', 'export const BASELINE = 99;\n');
    const result = execute();
    assert.equal(result.status, 1, `${name}: ${result.output}`);
    assert.match(result.output, missing(DESIGN_RULE_GUARD_TEST_PATH), name);
  }));
  // Python: TYPE_DEBT_BASELINE renamed everywhere in the file.
  withDesignRuleGuardRepository2(({ execute, write }) => {
    pythonProductionChange(write);
    const source = designRuleGuardSources2()[DESIGN_RULE_PYTHON_GUARD_PATH];
    write(DESIGN_RULE_PYTHON_GUARD_PATH, source.replaceAll('TYPE_DEBT_BASELINE', 'TYPE_DEBT_CEILING'));
    const result = execute();
    assert.equal(result.status, 1, result.output);
    assert.match(result.output, missing(DESIGN_RULE_PYTHON_GUARD_PATH));
  });
  // Without production code the same edits are test changes.
  withDesignRuleGuardRepository(({ execute, write }) => {
    write(DESIGN_RULE_GUARD_TEST_PATH, cases.let(designRuleGuardTestSource()));
    assert.doesNotMatch(execute().output, /design-rule\.count-scope/);
  });
});

// F2: the code that counts for a registered count baseline outside its test
// file (the shared comparison helper, the owner-graph report, the tsconfig the
// type ratchets count with) does not change alongside production code either.
const COUNTING_CODE_PATHS = Object.freeze([
  'tests/web/helpers/count-baseline.mjs',
  'tests/web/types/tsconfig.json',
  'tests/web/types/web-globals.d.ts',
  'tools/report-web-owner-graph.mjs'
]);
test('the counting code of the count baselines does not change alongside runtime files', () => {
  const source = (path) => (path.endsWith('.json')
    ? '{ "include": ["../../../gbdraw/web/js/**/*.js"] }\n'
    : 'export const compareCountBaseline = (row) => row.observed > row.recorded;\n');
  const weakened = (path) => source(path)
    .replace('"../../../gbdraw/web/js/**/*.js"', '"../../../gbdraw/web/js/app/*.js"')
    .replace('row.observed > row.recorded', 'false');
  COUNTING_CODE_PATHS.forEach((path) => {
    withDesignRuleGuardRepository(({ execute, root, write }) => {
      runtimeChange(write);
      write(path, weakened(path));
      const result = execute();
      assert.equal(result.status, 1, `${path}: ${result.output}`);
      assert.match(result.output, countScope(escaped(path)), path);
      write(path, source(path));
      rmSync(join(root, path));
      const deleted = execute();
      assert.equal(deleted.status, 1, `${path}: ${deleted.output}`);
      assert.match(deleted.output, countScope(`${escaped(path)} \\(deleted\\)`), path);
    }, {
      sources: { [DESIGN_RULE_GUARD_TEST_PATH]: designRuleGuardTestSource(), [path]: source(path) }
    });
    withDesignRuleGuardRepository(({ execute, write }) => {
      write(path, weakened(path));
      assert.doesNotMatch(execute().output, /design-rule\.count-scope/, path);
    }, {
      sources: { [DESIGN_RULE_GUARD_TEST_PATH]: designRuleGuardTestSource(), [path]: source(path) }
    });
  });
});

test('every count-baseline test counts only with counting code or guard paths', () => {
  const checker = readFileSync(CHANGE_BUDGET_CHECKER, 'utf8');
  const pathSet = (name) => (checker.match(new RegExp(`const ${name} = new Set\\(\\[([\\s\\S]*?)\\]\\);`))?.[1] || '')
    .split(',').map((item) => item.trim()).filter(Boolean)
    .map((item) => item.match(/^'([^']+)'$/)?.[1] || checker.match(new RegExp(`const ${item} = '([^']+)';`))?.[1]);
  assert.deepEqual(pathSet('COUNTING_CODE_PATHS').sort(), [...COUNTING_CODE_PATHS]);
  const guards = new Set(pathSet('guardPaths'));
  const registry = JSON.parse(readFileSync(join(REPOSITORY_ROOT, DESIGN_RULE_GUARD_REGISTRY_PATH), 'utf8'));
  const countFiles = new Set(registry.rules.flatMap((rule) => rule.allowlists)
    .filter(({ kind }) => /^(py-)?count(-map)?$/.test(kind))
    .map(({ path }) => path));
  // Every module a count-baseline test imports outside the package, directly
  // or through such a module, can change the count: it is counting code or a
  // guard path.
  const pending = [...countFiles].filter((path) => path.endsWith('.mjs'));
  const seen = new Set(pending);
  while (pending.length) {
    const path = pending.pop();
    const source = readFileSync(join(REPOSITORY_ROOT, path), 'utf8');
    for (const [, specifier] of source.matchAll(/^import [^;]*? from '(\.{1,2}\/[^']+)';$/gm)) {
      const target = posix.join(posix.dirname(path), specifier);
      if (target.startsWith('gbdraw/') || seen.has(target)) continue;
      seen.add(target);
      assert.ok(
        COUNTING_CODE_PATHS.includes(target) || guards.has(target),
        `${target} (imported by ${path}) counts for a registered count baseline but is neither counting code nor a guard path`
      );
      if (COUNTING_CODE_PATHS.includes(target)) pending.push(target);
    }
  }
  assert.ok(seen.has('tests/web/helpers/count-baseline.mjs') && seen.has('tools/report-web-owner-graph.mjs'));
  // The type ratchet counts with tsc (explicit-any-ratchet's TSCONFIG and
  // what it extends): every file its config includes outside the package (a
  // declaration file can type a read as `any`) is counting code, and no glob
  // there reaches outside it.
  const typeConfigs = ['tests/web/types/tsconfig.json'];
  const included = [];
  for (let config = typeConfigs.pop(); config; config = typeConfigs.pop()) {
    const parsed = JSON.parse(readFileSync(join(REPOSITORY_ROOT, config), 'utf8'));
    const resolveEntry = (entry) => posix.join(posix.dirname(config), entry);
    if (parsed.extends) typeConfigs.push(resolveEntry(parsed.extends));
    included.push(...[...(parsed.include || []), ...(parsed.files || [])].map(resolveEntry));
  }
  const outsidePackage = included.filter((path) => !path.startsWith('gbdraw/'));
  assert.ok(outsidePackage.length, included.join(', '));
  outsidePackage.forEach((path) => {
    assert.doesNotMatch(path, /[*?]/, `${path}: a type-ratchet glob reaches outside the package`);
    assert.ok(COUNTING_CODE_PATHS.includes(path), `${path} is included by a type-ratchet config but is not counting code`);
  });
});

// F3 (remedy b): a runtime move that a count-baseline test follows by path
// (#927 moved app/linear-comparisons.js to services/) cannot be split, so a
// path string that names the old path of a rename in the same diff may name
// the new path. Any other edit is still a count-scope change.
test('a count-baseline test follows a runtime move by path alongside runtime files', () => {
  const from = 'gbdraw/web/js/app/linear-comparisons.js';
  const to = 'gbdraw/web/js/services/linear-comparisons.js';
  const module = 'export const createDefaultLinearComparisonPlan = () => ({\n  comparisons: [],\n  thresholds: { identity: 0, length: 0 }\n});\n';
  const withPaths = (specifier, source) => [
    `import { createDefaultLinearComparisonPlan } from '${specifier}';`,
    `const SOURCE = '${source}'; // read by path`,
    'void createDefaultLinearComparisonPlan; void SOURCE;',
    designRuleGuardTestSource()
  ].join('\n');
  const before = withPaths('../../gbdraw/web/js/app/linear-comparisons.js', 'app/linear-comparisons.js');
  const after = withPaths('../../gbdraw/web/js/services/linear-comparisons.js', 'services/linear-comparisons.js');
  const run = (edit, { rename = true } = {}) => withDesignRuleGuardRepository(({ base, commit, execute, root, write }) => {
    runtimeChange(write);
    if (rename) {
      rmSync(join(root, from));
      write(to, module);
    } else {
      write(from, 'export const createDefaultLinearComparisonPlan = () => null;\n');
      write(to, 'export const other = 1;\n');
    }
    write(DESIGN_RULE_GUARD_TEST_PATH, edit);
    return execute({ base, head: commit('move') });
  }, { sources: { [DESIGN_RULE_GUARD_TEST_PATH]: before, [from]: module } });
  const moved = run(after);
  assert.doesNotMatch(moved.output, /design-rule\.count-scope/, moved.output);
  // A lowered literal is still allowed with the move.
  assert.doesNotMatch(run(after.replace('const BASELINE = 46;', 'const BASELINE = 45;')).output, /design-rule\.count-scope/);
  const refused = {
    'another path': after.replace("'services/linear-comparisons.js'", "'services/other.js'"),
    'a path without a rename': after,
    'another edit': after.replace('// read by path', '').replace('void SOURCE;', 'void SOURCE; void 0;'),
    'a narrowed count': after.replace('test(\'baseline\'', 'test.skip(\'baseline\'')
  };
  Object.entries(refused).forEach(([name, edit]) => {
    const result = run(edit, { rename: name !== 'a path without a rename' });
    assert.equal(result.status, 1, `${name}: ${result.output}`);
    assert.match(result.output, countScope(escaped(DESIGN_RULE_GUARD_TEST_PATH)), name);
  });
});

// BA-D10: a rule's text is authority. The Gate checks the candidate registry
// whenever it changes; an unchanged trusted registry is not re-checked.
// json-ids (schema 2): the ids of a JSON list are a shrink-only set, and so is
// the length of each row's lists (R9: an owner added to a twin row is
// authority-only); a renamed item and the row's other fields change freely.
test('a json-ids allowlist compares the entry ids and the length of each row list', () => {
  const TWINS_PATH = 'tools/twins.json';
  const row = (id, extra = {}) => ({ id, pythonOwner: ['gbdraw/a.py:f'], jsOwner: ['gbdraw/web/js/services/a.js:f'], reason: 'r', retireWhen: 'w', ...extra });
  const twins = (rows) => `${JSON.stringify({ schemaVersion: 1, twins: rows }, null, 2)}\n`;
  const rules = designRuleRules2();
  rules[0].allowlists.push({ path: TWINS_PATH, symbol: 'twins', kind: 'json-ids' });
  const subject = `R10 ${TWINS_PATH}#twins`.replace(/[.*+?^${}()|[\]\\]/g, '\\$&');
  const rowA = (extra = {}) => row('function:a', { pythonOwner: ['gbdraw/a.py:f', 'gbdraw/a.py:g'], ...extra });
  const base = [rowA(), row('function:b', { pending: true })];
  const run = (rows) => withDesignRuleGuardRepository(({ execute, write }) => {
    runtimeChange(write);
    write(TWINS_PATH, twins(rows));
    return execute();
  }, {
    registry: designRuleGuardRegistrySource(rules, 2),
    sources: { ...designRuleGuardSources2(), [TWINS_PATH]: twins(base) }
  });
  const added = run([...base, row('function:c')]);
  assert.equal(added.status, 1, added.output);
  assert.match(added.output, new RegExp(`design-rule\\.co-change: [^\\n]*${subject} \\(expansion\\)`));
  const removed = run([base[0]]);
  assert.equal(removed.status, 0, removed.output);
  assert.match(removed.output, new RegExp(`- contraction: ${subject}\\n`));
  const edited = run([rowA({ jsOwner: ['gbdraw/web/js/services/b.js:f'], reason: 'moved', retireWhen: 'later' }), row('function:b')]);
  assert.equal(edited.status, 0, edited.output);
  assert.doesNotMatch(edited.output, new RegExp(subject));
  const ownerAdded = run([rowA({ jsOwner: ['gbdraw/web/js/services/a.js:f', 'gbdraw/web/js/services/b.js:f'] }), base[1]]);
  assert.equal(ownerAdded.status, 1, ownerAdded.output);
  assert.match(ownerAdded.output, new RegExp(`design-rule\\.co-change: [^\\n]*${subject} \\(expansion\\)`));
  const ownerMoved = run([rowA({ pythonOwner: ['gbdraw/a.py:f'], jsOwner: ['gbdraw/web/js/services/a.js:f', 'gbdraw/web/js/services/b.js:f'] }), base[1]]);
  assert.equal(ownerMoved.status, 1, ownerMoved.output);
  assert.match(ownerMoved.output, new RegExp(`design-rule\\.co-change: [^\\n]*${subject} \\(expansion\\)`));
  const ownerRemoved = run([rowA({ pythonOwner: ['gbdraw/a.py:f'] }), base[1]]);
  assert.equal(ownerRemoved.status, 0, ownerRemoved.output);
  assert.match(ownerRemoved.output, new RegExp(`- contraction: ${subject}\\n`));
});

test('a candidate guard registry names only authority paths as rule files', () => {
  const rulesPathError = /candidate design-rule guard registry: rules\[1\] TYPE-DEBT: rulesPath CLAUDE\.md is not an authority path/;
  withDesignRuleGuardRepository2(({ execute, write }) => {
    const unchanged = execute();
    assert.doesNotMatch(unchanged.output, /candidate design-rule guard registry/);
    write(DESIGN_RULE_GUARD_REGISTRY_PATH, `${designRuleGuardRegistrySource(designRuleRules2(), 2)}\n`);
    const reformatted = execute();
    assert.equal(reformatted.status, 1, reformatted.output);
    assert.match(reformatted.output, rulesPathError);
    const rules = designRuleRules2();
    rules[1].rulesPath = 'docs/internal/WEB_CHANGE_POLICY.md';
    write(DESIGN_RULE_GUARD_REGISTRY_PATH, designRuleGuardRegistrySource(rules, 2));
    const moved = execute();
    assert.doesNotMatch(moved.output, /candidate design-rule guard registry/);
  });
});

test('the twin detector is pre-registered checker implementation and a guard path', () => {
  const checker = readFileSync(CHANGE_BUDGET_CHECKER, 'utf8');
  assert.match(checker, /const twinDetectorsPath = 'tools\/web-twin-detectors\.mjs';/);
  const guardBlock = checker.match(/const guardPaths = new Set\(\[([\s\S]*?)\]\);/)[1];
  const implementationBlock = checker.match(/const checkerImplementationPaths = new Set\(\[([\s\S]*?)\]\);/)[1];
  assert.match(guardBlock, /twinDetectorsPath/);
  assert.match(implementationBlock, /twinDetectorsPath/);
  const result = runChangeBudgetCase((write) => {
    write('gbdraw/web/js/services/session-file.js', 'export const readSession = () => ({ version: 2 });\n');
    write('tools/web-twin-detectors.mjs', 'export const TWIN_DETECTOR_ID = "twin.unregistered.v1";\n');
  });
  assert.equal(result.status, 1, result.output);
  assert.match(result.output, /production runtime files and Web guard\/CI files changed together/);
});

test('the Web design rules and the guard registry are pre-registered guard and authority paths', () => {
  const checker = readFileSync(CHANGE_BUDGET_CHECKER, 'utf8');
  assert.match(checker, /const webDesignRulesPath = 'gbdraw\/web\/CLAUDE\.md';/);
  assert.match(checker, /const designRuleGuardRegistryPath = 'tools\/web-design-rule-guards\.json';/);
  const guardBlock = checker.match(/const guardPaths = new Set\(\[([\s\S]*?)\]\);/)[1];
  const authorityBlock = checker.match(/const authorityPaths = new Set\(\[([\s\S]*?)\]\);/)[1];
  for (const block of [guardBlock, authorityBlock]) {
    assert.match(block, /webDesignRulesPath/);
    assert.match(block, /designRuleGuardRegistryPath/);
  }
  if (existsSync(join(REPOSITORY_ROOT, DESIGN_RULE_GUARD_REGISTRY_PATH))) {
    // Once the registry is activated, every rule it names exists in its rules
    // file, every guard exists, and the Gate parses the registry and reads
    // every allowlist literal with its registered kind.
    const source = readFileSync(join(REPOSITORY_ROOT, DESIGN_RULE_GUARD_REGISTRY_PATH), 'utf8');
    assertDesignRuleRegistryResolves(JSON.parse(source), REPOSITORY_ROOT);
    assertGateReadsDesignRuleRegistry(source, REPOSITORY_ROOT);
  }
});

// The guard registry in both directions (boundary audit M5a). (i) Every
// `tests/` path the Web rules name is registered, and every registered path is
// named by its rule. (ii) Every guard is selected by the CI job it declares.
// (iii) Every node test under tests/web is selected by a pull-request job.
// Selection is read from the workflow's commands, the functional shard file,
// and the pytest markers of the file, never from a table kept here; a command
// the resolver cannot read fails by name.
const expandBraces = (path) => {
  const braces = path.match(/^(.*?)\{([^{}]*)\}(.*)$/);
  return braces ? braces[2].split(',').flatMap((part) => expandBraces(`${braces[1]}${part}${braces[3]}`)) : [path];
};
const namedTestPaths = (text) => [...new Set(
  (text.match(/\btests\/[\w./{},-]*\.(?:mjs|cjs|js|json|py|ts)\b/g) || []).flatMap(expandBraces)
)];
// A heading's text up to the next heading of its level or above, without code blocks.
const ruleSection = (text, heading) => {
  const lines = text.replace(/^```[\s\S]*?^```/gm, '').split(/\r?\n/);
  const start = lines.findIndex((line) => /^#{1,6} /.test(line) && line.replace(/^#+ /, '') === heading);
  if (start < 0) return '';
  const closes = new RegExp(`^#{1,${lines[start].indexOf(' ')}} `);
  const length = lines.slice(start + 1).findIndex((line) => closes.test(line));
  return lines.slice(start, length < 0 ? lines.length : start + 1 + length).join('\n');
};
// The shell command lines of a job's `run:` steps: a folded scalar is one
// line, and a literal one is its lines with `\` continuations joined.
const workflowCommands = (job) => [...`${job}\n`.matchAll(/\n( +)run: *([>|]-?)? *(.*)\n((?:\1 +.*\n|[ \t]*\n)*)/g)]
  .flatMap(([, , style, inline, block]) => (!style ? [inline]
    : (style[0] === '|' ? block.replace(/\\\n/g, ' ') : block.replace(/\n/g, ' ')).split('\n')))
  .map((line) => line.replace(/\s+/g, ' ').trim()).filter(Boolean);
const globMatcher = (glob) => {
  const pattern = new RegExp(`^${glob.split('*').map((part) => part.replace(/[.+?^${}()|[\]\\]/g, '\\$&')).join('[^/]*')}$`);
  return (text) => pattern.test(text);
};
// What a job's `node --test` commands run: `find <root> [-maxdepth N] -type f
// [!] -name '<glob>' ... | xargs -0 node --test`, or the files named after it.
const nodeTestSelectors = (commands) => commands.filter((command) => /\bnode --test\b/.test(command)).map((command) => {
  const found = command.match(/^find (\S+) (.*)-print0 \| xargs -0 node --test$/);
  if (found) {
    const [, root, predicates] = found;
    const tests = [];
    const unread = predicates.replace(/(! )?-(maxdepth|type|name) (?:'([^']*)'|(\S+)) /g, (all, negated, key, quoted, bare) => {
      const value = quoted ?? bare;
      if (key === 'name') {
        const matches = globMatcher(value);
        tests.push((path) => matches(posix.basename(path)) !== Boolean(negated));
      } else if (key === 'maxdepth' && !negated) tests.push((path) => path.split('/').length <= Number(value));
      else if (key !== 'type' || value !== 'f' || negated) return all;
      return '';
    });
    assert.equal(unread.trim(), '', `cannot resolve the find predicates of: ${command}`);
    return (path) => path.startsWith(`${root}/`) && tests.every((matches) => matches(path.slice(root.length + 1)));
  }
  const named = command.match(/^node --test((?: \S+)*)$/);
  assert.ok(named, `cannot resolve the node test command: ${command}`);
  const files = named[1].split(' ').filter((argument) => argument && !argument.startsWith('--')).map(globMatcher);
  return (path) => files.some((matches) => matches(path));
});
// The marker sets of a pytest file's tests: the file's `pytestmark`, the
// decorators of a test's class, and the test's own decorators.
const pytestMarkerSets = (source) => {
  const markersIn = (text) => [...text.matchAll(/pytest\.mark\.(\w+)/g)].map(([, name]) => name);
  const fileLevel = markersIn(source.match(/^pytestmark\s*=.*(?:\n[ \t)\]].*)*/m)?.[0] ?? '');
  const sets = [];
  let pending = [];
  let classMarkers = [];
  source.split(/\r?\n/).forEach((line) => {
    if (/^\s*@/.test(line)) pending.push(...markersIn(line));
    const [, indent, keyword, name] = line.match(/^(\s*)(?:async\s+)?(class|def)\s+(\w+)/) || [];
    if (!keyword) return;
    if (!indent) classMarkers = keyword === 'class' ? pending : [];
    if (keyword === 'def' && name.startsWith('test')) sets.push(new Set([...fileLevel, ...(indent ? classMarkers : []), ...pending]));
    pending = [];
  });
  return sets;
};
// What a job's `python -m pytest <target> [-m "<expression>"]` commands run:
// a file under the target with a test whose markers satisfy the expression.
const pytestSelectors = (commands) => commands.filter((command) => /\bpython -m pytest\b/.test(command)).map((command) => {
  const [, target, options = ''] = command.match(/^python -m pytest (\S+)(.*)$/) || [];
  const expression = options.match(/ -m "([\w\s()]*)"/);
  assert.ok(
    target && !/ (?:-m|-k|--deselect|--ignore)\b/.test(options.replace(expression?.[0] ?? '', '')),
    `cannot resolve the pytest command: ${command}`
  );
  const source = (expression?.[1] ?? '').replace(/\w+/g, (word) => (
    { and: '&&', or: '||', not: '!' }[word] ?? `markers.has('${word}')`
  ));
  const selects = source.trim() ? new Function('markers', `return ${source};`) : () => true;
  return (path, markerSets) => (target.endsWith('/') ? path.startsWith(target) : path === target)
    && markerSets.some((markers) => selects(markers));
});
const guardRegistryProblems = (registry, {
  read = (path) => readFileSync(join(REPOSITORY_ROOT, path), 'utf8'),
  workflow = TEST_WORKFLOW,
  nodeTests = readdirSync(join(REPOSITORY_ROOT, 'tests/web'), { recursive: true })
    .map((path) => `tests/web/${normalizePath(path)}`).filter((path) => path.endsWith('.test.mjs')).sort()
} = {}) => {
  const problems = { unregistered: [], unnamed: [], unrun: [] };
  const registeredPaths = (rule) => [...rule.guards, ...rule.allowlists].map(({ path }) => path);
  const registered = new Set(registry.rules.flatMap(registeredPaths));
  namedTestPaths(read(registry.rulesPath)).filter((path) => !registered.has(path)).forEach((path) => {
    problems.unregistered.push(`${registry.rulesPath} names ${path}, which no rule registers`);
  });
  const jobIds = [...workflow.slice(workflow.indexOf('\njobs:\n')).matchAll(/\n {2}([a-z][a-z0-9-]*):\n/g)].map(([, id]) => id);
  const jobCommands = new Map();
  const commandsOf = (job) => jobCommands.get(job)
    ?? jobCommands.set(job, workflowCommands(workflowJob(job, workflow))).get(job);
  const shards = new Set(JSON.parse(read('tests/ci/functional-shards.json')).shards.flat());
  // true or false when the job's commands decide it; undefined when the path is not a test of a form they run.
  const selectedBy = (job, path) => {
    const commands = commandsOf(job);
    if (path.endsWith('.test.mjs')) return nodeTestSelectors(commands).some((selects) => selects(path));
    if (/(?:^|\/)test_[^/]*\.py$/.test(path)) {
      const markerSets = pytestMarkerSets(read(path));
      return pytestSelectors(commands).some((selects) => selects(path, markerSets));
    }
    if (path.endsWith('.playwright.spec.js') && commands.some((command) => /\bci-impact\.mjs shard-files\b/.test(command))) {
      return shards.has(path);
    }
    return undefined;
  };
  registry.rules.forEach((rule) => {
    // A rule may also register a JSON allowlist under tools/ that it names.
    const section = ruleSection(read(rule.rulesPath || registry.rulesPath), rule.heading);
    const named = new Set([...namedTestPaths(section), ...(section.match(/\btools\/[\w.-]+\.json\b/g) || [])]);
    [...new Set(registeredPaths(rule))].filter((path) => !named.has(path)).forEach((path) => {
      problems.unnamed.push(`${rule.id} registers ${path}, which "${rule.heading}" does not name`);
    });
    const unrunReason = ({ path, job }, seen = []) => {
      if (!jobIds.includes(job)) return 'is not a job of .github/workflows/test.yml';
      const selected = selectedBy(job, path);
      if (selected !== undefined) return selected ? '' : 'does not select it';
      if (/\.spec\.js$/.test(path)) return 'has no selection this check can read for it';
      // A fixture, helper or config runs with a guard of its rule and job that names it.
      return rule.guards.some((other) => other.job === job && other.path !== path && !seen.includes(other.path)
        && read(other.path).includes(posix.basename(path)) && !unrunReason(other, [...seen, path]))
        ? '' : 'runs no guard of this rule that names it';
    };
    rule.guards.forEach((guard) => {
      const reason = unrunReason(guard);
      if (reason) problems.unrun.push(`${rule.id} ${guard.path}: ${guard.job} ${reason}`);
    });
  });
  nodeTests.filter((path) => !knownJobsFor('pr').some((job) => selectedBy(job, path))).forEach((path) => {
    problems.unrun.push(`${path}: no pull-request job runs this node test`);
  });
  return problems;
};

test('every test path the Web rules name is registered, and every guard and node test runs in its CI job', (t) => {
  const registry = JSON.parse(readFileSync(join(REPOSITORY_ROOT, DESIGN_RULE_GUARD_REGISTRY_PATH), 'utf8'));
  const problems = guardRegistryProblems(registry);
  t.diagnostic(`guard registry: ${registry.rules.length} rules, ${registry.rules.flatMap((rule) => rule.guards).length} guards; `
    + `${problems.unregistered.length} unregistered named paths, ${problems.unnamed.length} registered paths a rule does not name, `
    + `${problems.unrun.length} unrun guards or node tests`);
  assert.deepEqual(problems, { unregistered: [], unnamed: [], unrun: [] });
});

test('the guard registry check reports an unregistered path, an unnamed guard, an unrun guard, and an unreadable command', () => {
  const registry = JSON.parse(readFileSync(join(REPOSITORY_ROOT, DESIGN_RULE_GUARD_REGISTRY_PATH), 'utf8'));
  const read = (path) => readFileSync(join(REPOSITORY_ROOT, path), 'utf8');
  const rule = (id) => registry.rules.find((candidate) => candidate.id === id);
  assert.deepEqual(guardRegistryProblems(registry, {
    read: (path) => read(path) + (path === registry.rulesPath ? '\nSee `tests/web/{unregistered,run-info}.test.mjs`.\n' : '')
  }).unregistered, ['gbdraw/web/CLAUDE.md names tests/web/unregistered.test.mjs, which no rule registers']);

  const moved = structuredClone(registry);
  moved.rules.find(({ id }) => id === 'R8').guards.push(rule('R12').guards.find(({ path }) => path.endsWith('run-info.test.mjs')));
  assert.deepEqual(guardRegistryProblems(moved).unnamed, [
    'R8 registers tests/web/run-info.test.mjs, which "R8: The Multi-Record Canvas owns placement only" does not name'
  ]);
  // A rule may register a tools/*.json allowlist that its section names; an
  // unnamed one is still reported.
  const withTwins = structuredClone(registry);
  withTwins.rules.find(({ id }) => id === 'R13').allowlists.push({ path: 'tools/twins.json', symbol: 'twins', kind: 'json-ids' });
  const namingTwins = (path) => read(path).replace('### R13: Owner layering and ports\n', '$&\nTwins: `tools/twins.json`.\n');
  assert.deepEqual(guardRegistryProblems(withTwins, { read: namingTwins }).unnamed, []);
  assert.deepEqual(guardRegistryProblems(withTwins).unnamed, [
    'R13 registers tools/twins.json, which "R13: Owner layering and ports" does not name'
  ]);

  // A node test no job selects, and the fixture that only that test reads.
  const excluding = (name) => TEST_WORKFLOW.replaceAll("-name '*.test.mjs'", `$& ! -name '${name}'`);
  assert.deepEqual(guardRegistryProblems(registry, { workflow: excluding('record-metadata-inference.test.mjs') }).unrun, [
    'R4 tests/fixtures/record_metadata_inference_cases.json: web-contracts-pr runs no guard of this rule that names it',
    'R4 tests/web/record-metadata-inference.test.mjs: web-contracts-pr does not select it',
    'tests/web/record-metadata-inference.test.mjs: no pull-request job runs this node test'
  ]);
  const otherJob = structuredClone(registry);
  otherJob.rules.find(({ id }) => id === 'R8').guards[0].job = 'web-contracts-pr';
  otherJob.rules.find(({ id }) => id === 'R12').guards[0].job = 'promotion';
  assert.deepEqual(guardRegistryProblems(otherJob).unrun, [
    'R8 tests/test_circular_multi_record_parity.py: web-contracts-pr does not select it',
    `R12 ${rule('R12').guards[0].path}: promotion is not a job of .github/workflows/test.yml`
  ]);
  // A pytest file runs in a job when one of its tests does: one `browser`
  // test leaves the others in core-pr, and a file of only such tests leaves it.
  const withBrowserTests = (flags) => (path) => (path === 'tests/test_mode_scoped_settings.py'
    ? read(path).replace(new RegExp('^def test_', flags), '@pytest.mark.browser\n$&') : read(path));
  assert.deepEqual(guardRegistryProblems(registry, { read: withBrowserTests('m') }).unrun, []);
  assert.deepEqual(guardRegistryProblems(registry, { read: withBrowserTests('gm') }).unrun, [
    'R15 tests/test_mode_scoped_settings.py: core-pr does not select it'
  ]);
  assert.throws(
    () => guardRegistryProblems(registry, { workflow: TEST_WORKFLOW.replaceAll("-name '*.test.mjs'", "-newer package.json $&") }),
    /cannot resolve the find predicates of: find tests\/web /
  );
});

// Owner-graph report (implementation plan Phase B1): the checker prints the
// detector summaries for the base and the head and the NEW and REMOVED
// subjects, and nothing in that section affects the Gate until Phase B3.
test('the owner-graph report lists base and head counts and new subjects without gating', () => {
  const registry = JSON.stringify({
    schemaVersion: 1,
    compositionRoots: ['app/wiring.js'],
    ownerObjectNames: ['store', 'legendOwner'],
    stateModule: 'state.js',
    projectionDomains: [],
    heavyProducers: []
  });
  const root = [
    "import { createLegendOwner } from './legend-owner.js';",
    "import { createStore } from '../services/store.js';",
    'export const createWiring = () => {',
    '  const store = createStore({ state });',
    '  const legendOwner = createLegendOwner({ state, store });',
    '  return { store, legendOwner };',
    '};',
    ''
  ].join('\n');
  const withClosure = root.replace(
    '  const store = createStore({ state });',
    '  const store = createStore({ state, buildLegendEntryOwners: () => legendOwner.captureLegendEntryOwners() });'
  );
  withChangeBudgetRepository(({ commit, execute, write }) => {
    write(OWNER_GRAPH_REGISTRY_PATH, `${registry}\n`);
    write('gbdraw/web/js/app/wiring.js', root);
    write('gbdraw/web/js/app/legend-owner.js', 'export const createLegendOwner = ({ state, store }) => ({ state, store });\n');
    write('gbdraw/web/js/services/store.js', 'export const createStore = ({ state }) => ({ state });\n');
    commit('baseline composition');
    const unchanged = execute();
    assert.equal(unchanged.status, 0, unchanged.output);
    assert.match(unchanged.output, /## Owner graph \(report\)\n\n- Registry: tools\/web-owner-graph\.json from the trusted base\n- injectionEdges: 1 -> 1\n- forwardClosures: 0 -> 0\n- wholeObjectPorts: 1 -> 1/);
    write('gbdraw/web/js/app/wiring.js', withClosure);
    const result = execute();
    assert.equal(result.status, 0, result.output);
    assert.match(result.output, /Gate: \*\*PASS\*\*/);
    assert.match(result.output, /- forwardClosures: 0 -> 1\n/);
    assert.match(result.output, /- owner-graph\.forward-closure\.v1: NEW app\/wiring\.js\|store\.buildLegendEntryOwners->legendOwner\.captureLegendEntryOwners/);
    assert.doesNotMatch(result.output, /Blocking violations\n\n- .*owner-graph/);
  });
});

test('the owner-graph registry is pre-registered authority and is read from the trusted base', () => {
  const checker = readFileSync(CHANGE_BUDGET_CHECKER, 'utf8');
  assert.match(checker, /const ownerGraphRegistryPath = 'tools\/web-owner-graph\.json';/);
  assert.match(checker, /const ownerGraphDetectorsPath = 'tools\/web-owner-graph-detectors\.mjs';/);
  const guardBlock = checker.match(/const guardPaths = new Set\(\[([\s\S]*?)\]\);/)[1];
  const authorityBlock = checker.match(/const authorityPaths = new Set\(\[([\s\S]*?)\]\);/)[1];
  const implementationBlock = checker.match(/const checkerImplementationPaths = new Set\(\[([\s\S]*?)\]\);/)[1];
  assert.match(guardBlock, /ownerGraphRegistryPath/);
  assert.match(guardBlock, /ownerGraphDetectorsPath/);
  assert.match(authorityBlock, /ownerGraphRegistryPath/);
  assert.match(implementationBlock, /ownerGraphDetectorsPath/);
  assert.match(checker, /readRevisionFile\(base, ownerGraphRegistryPath\)/);
  withChangeBudgetRepository(({ commit, execute, write }) => {
    write(OWNER_GRAPH_REGISTRY_PATH, '{"schemaVersion":1,"compositionRoots":["app/app-setup.js"]}\n');
    commit('register owner graph');
    write(OWNER_GRAPH_REGISTRY_PATH, '{ malformed\n');
    const result = execute();
    assert.match(result.output, /- Registry: tools\/web-owner-graph\.json from the trusted base/);
  });
  const malformed = runChangeBudgetCase((write) => {
    write(OWNER_GRAPH_REGISTRY_PATH, '{ malformed\n');
  });
  // An uncommitted registry is read from the base (absent) and reported as defaults.
  assert.match(malformed.output, /- Registry: defaults \(registry absent\)/);
});
