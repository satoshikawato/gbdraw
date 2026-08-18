import assert from 'node:assert/strict';
import { spawnSync } from 'node:child_process';
import {
  mkdirSync,
  mkdtempSync,
  readFileSync,
  rmSync,
  writeFileSync
} from 'node:fs';
import { tmpdir } from 'node:os';
import { dirname, join, resolve } from 'node:path';
import test from 'node:test';
import { fileURLToPath } from 'node:url';

import { WEB_ARCHITECTURE_DETECTORS } from '../../tools/web-architecture-detectors.mjs';
import {
  classifyArchitectureAuthorityDelta,
  classifyArchitectureRuleObservation,
  validateArchitectureRuleRegistry,
  WEB_ARCHITECTURE_RULE_SCHEMA
} from '../../tools/web-architecture-evaluation.mjs';

const REPOSITORY_ROOT = resolve(dirname(fileURLToPath(import.meta.url)), '../..');
const AVAILABLE_ENFORCEMENTS = Object.freeze({
  availableEnforcements: Object.freeze(['hard', 'report-only'])
});

const canonicalRule = (overrides = {}) => ({
  key: 'canonical-path.render-request',
  kind: 'single-canonical-entry-edge',
  detector: 'canonical-path.render-request.v1',
  allowedEdges: [{
    from: 'app/run-analysis.js',
    to: 'services/session-request.js'
  }],
  exactActiveEdgeCount: 1,
  enforcement: 'hard',
  baselineEligible: false,
  ...overrides
});

const semanticRule = (overrides = {}) => ({
  key: 'semantic-owner.render-request',
  kind: 'single-semantic-owner',
  detector: 'semantic-owner.render-request.v1',
  allowedDefinitionPaths: ['services/session-request.js'],
  exactDefinitionCount: 1,
  enforcement: 'hard',
  baselineEligible: false,
  ...overrides
});

const registry = (rules = [canonicalRule(), semanticRule()]) => ({
  schemaVersion: 1,
  rules
});
const clone = (value) => JSON.parse(JSON.stringify(value));
const validationCodes = (value, catalog = WEB_ARCHITECTURE_DETECTORS, options = {}) => (
  validateArchitectureRuleRegistry(value, catalog, options).errors.map(({ code }) => code)
);

test('schema version 1 accepts only the two discriminated initial rule kinds', () => {
  assert.deepEqual(WEB_ARCHITECTURE_RULE_SCHEMA, {
    schemaVersion: 1,
    maximumRuleCount: 3,
    kinds: ['single-canonical-entry-edge', 'single-semantic-owner'],
    enforcementModes: ['frozen', 'hard', 'report-only']
  });
  assert.deepEqual(
    validateArchitectureRuleRegistry(
      registry(),
      WEB_ARCHITECTURE_DETECTORS,
      AVAILABLE_ENFORCEMENTS
    ),
    { valid: true, errors: [] }
  );
});

test('strict schema rejects unknown, executable-looking, duplicate, and unsorted data', () => {
  const cases = [
    {
      name: 'unknown top-level field',
      value: { ...registry(), predicate: 'process.exit(0)' },
      code: 'unknown-field'
    },
    {
      name: 'wrong kind-specific field',
      value: registry([{
        ...semanticRule(),
        allowedEdges: canonicalRule().allowedEdges
      }]),
      code: 'unknown-field'
    },
    {
      name: 'wildcard path',
      value: registry([semanticRule({ allowedDefinitionPaths: ['services/*.js'] })]),
      code: 'invalid-web-module-path'
    },
    {
      name: 'JavaScript-looking path',
      value: registry([semanticRule({
        allowedDefinitionPaths: ['services/session-request.js;process.exit.js']
      })]),
      code: 'invalid-web-module-path'
    },
    {
      name: 'repository-prefixed path',
      value: registry([semanticRule({
        allowedDefinitionPaths: ['gbdraw/web/js/services/session-request.js']
      })]),
      code: 'invalid-web-module-path'
    },
    {
      name: 'duplicate path',
      value: registry([semanticRule({
        allowedDefinitionPaths: [
          'services/session-request.js',
          'services/session-request.js'
        ]
      })]),
      code: 'duplicate-value'
    },
    {
      name: 'duplicate rule key',
      value: registry([semanticRule(), semanticRule()]),
      code: 'duplicate-rule-key'
    },
    {
      name: 'unsupported enforcement',
      value: registry([semanticRule({ enforcement: 'disabled' })]),
      code: 'unsupported-enforcement'
    },
    {
      name: 'unsupported definition count',
      value: registry([semanticRule({ exactDefinitionCount: 2 })]),
      code: 'invalid-definition-count'
    },
    {
      name: 'unsorted rules',
      value: registry([semanticRule(), canonicalRule()]),
      code: 'unsorted-rules'
    },
    {
      name: 'more than three rules',
      value: registry([
        canonicalRule(),
        semanticRule({ key: 'semantic-owner.alpha' }),
        semanticRule({ key: 'semantic-owner.beta' }),
        semanticRule({ key: 'semantic-owner.gamma' })
      ]),
      code: 'too-many-rules'
    }
  ];

  cases.forEach(({ name, value, code }) => {
    assert.ok(validationCodes(value).includes(code), name);
  });
});

test('canonical edge schema requires nonempty unique sorted edges and one active edge', () => {
  const edgeA = { from: 'app/alternate.js', to: 'services/session-request.js' };
  const edgeB = { from: 'app/run-analysis.js', to: 'services/session-request.js' };
  const cases = [
    [canonicalRule({ allowedEdges: [] }), 'empty-array'],
    [canonicalRule({ allowedEdges: [edgeA, edgeA] }), 'duplicate-edge'],
    [canonicalRule({ allowedEdges: [edgeB, edgeA] }), 'unsorted-edges'],
    [canonicalRule({ exactActiveEdgeCount: 0 }), 'invalid-active-edge-count'],
    [canonicalRule({ exactActiveEdgeCount: 2 }), 'invalid-active-edge-count']
  ];

  cases.forEach(([rule, code]) => {
    assert.ok(validationCodes(registry([rule])).includes(code));
  });
});

test('detector compatibility and staged enforcement fail closed', () => {
  const wrongCategoryCatalog = {
    ...WEB_ARCHITECTURE_DETECTORS,
    'semantic-owner.render-request.v1': {
      ...WEB_ARCHITECTURE_DETECTORS['semantic-owner.render-request.v1'],
      subjectCategory: 'canonical-entry-edge'
    }
  };
  const missingEncoderCatalog = {
    ...WEB_ARCHITECTURE_DETECTORS,
    'semantic-owner.render-request.v1': {
      subjectCategory: 'definition-path'
    }
  };
  assert.ok(validationCodes(
    registry([semanticRule({ detector: 'semantic-owner.unknown.v1' })])
  ).includes('unknown-detector'));
  assert.ok(validationCodes(
    registry([semanticRule()]),
    wrongCategoryCatalog
  ).includes('detector-kind-mismatch'));
  assert.ok(validationCodes(
    registry([semanticRule()]),
    missingEncoderCatalog
  ).includes('missing-subject-encoder'));
  assert.ok(validationCodes(
    registry([semanticRule({ enforcement: 'frozen', baselineEligible: true })]),
    WEB_ARCHITECTURE_DETECTORS,
    AVAILABLE_ENFORCEMENTS
  ).includes('unavailable-enforcement'));
  assert.ok(validationCodes(
    registry([semanticRule({ baselineEligible: true })])
  ).includes('invalid-hard-eligibility'));
});

test('single-canonical-entry observation separates preauthorization from active edges', () => {
  const rule = canonicalRule({
    allowedEdges: [
      { from: 'app/alternate.js', to: 'services/session-request.js' },
      { from: 'app/run-analysis.js', to: 'services/session-request.js' }
    ]
  });
  const edge = (from, to = 'services/session-request.js') => ({
    from,
    to,
    subject: `${from} -> ${to}`
  });

  assert.equal(classifyArchitectureRuleObservation(rule, {
    observedEdges: [edge('app/run-analysis.js')]
  }).observation, 'CONFORMING');
  assert.equal(classifyArchitectureRuleObservation(rule, {
    observedEdges: []
  }).observation, 'ABSENT_REQUIRED');
  assert.equal(classifyArchitectureRuleObservation(rule, {
    observedEdges: [edge('app/unapproved.js')]
  }).observation, 'DIVERGENT');
  assert.equal(classifyArchitectureRuleObservation(rule, {
    observedEdges: [edge('app/alternate.js'), edge('app/run-analysis.js')]
  }).observation, 'DIVERGENT');
});

test('single-owner observation requires one allowed active definition', () => {
  const rule = semanticRule({
    allowedDefinitionPaths: [
      'services/future-session-request.js',
      'services/session-request.js'
    ]
  });
  const definition = (path, count = 1) => ({ path, count, subject: path });

  assert.equal(classifyArchitectureRuleObservation(rule, {
    definitionCount: 1,
    observedDefinitions: [definition('services/session-request.js')]
  }).observation, 'CONFORMING');
  assert.equal(classifyArchitectureRuleObservation(rule, {
    definitionCount: 0,
    observedDefinitions: []
  }).observation, 'ABSENT_REQUIRED');
  assert.equal(classifyArchitectureRuleObservation(rule, {
    definitionCount: 1,
    observedDefinitions: [definition('services/unapproved.js')]
  }).observation, 'DIVERGENT');
  assert.equal(classifyArchitectureRuleObservation(rule, {
    definitionCount: 2,
    observedDefinitions: [definition('services/session-request.js', 2)]
  }).observation, 'DIVERGENT');
});

test('authority delta reports every planned direction without treating removal as contraction', () => {
  const base = registry([semanticRule({
    allowedDefinitionPaths: [
      'services/future-session-request.js',
      'services/session-request.js'
    ]
  })]);
  const contraction = registry([semanticRule()]);
  assert.deepEqual(classifyArchitectureAuthorityDelta(base, contraction).directions, [
    'CONTRACTION'
  ]);
  assert.deepEqual(classifyArchitectureAuthorityDelta(contraction, base).directions, [
    'EXPANSION'
  ]);
  assert.deepEqual(classifyArchitectureAuthorityDelta(registry([]), contraction).directions, [
    'EXPANSION'
  ]);
  assert.deepEqual(classifyArchitectureAuthorityDelta(contraction, registry([])).directions, [
    'RELAXATION'
  ]);

  const reportOnly = registry([semanticRule({
    enforcement: 'report-only',
    baselineEligible: true
  })]);
  assert.deepEqual(classifyArchitectureAuthorityDelta(contraction, reportOnly).directions, [
    'RELAXATION'
  ]);
  assert.deepEqual(classifyArchitectureAuthorityDelta(reportOnly, contraction).directions, [
    'TIGHTENING'
  ]);

  const reinterpreted = registry([semanticRule({ detector: 'semantic-owner.render-request.v2' })]);
  assert.deepEqual(classifyArchitectureAuthorityDelta(contraction, reinterpreted).directions, [
    'REINTERPRETATION'
  ]);
});

test('evaluation leaves frozen inputs unchanged and returns stable plain data', () => {
  const deepFreeze = (value) => {
    if (value && typeof value === 'object') {
      Object.values(value).forEach(deepFreeze);
      Object.freeze(value);
    }
    return value;
  };
  const value = deepFreeze(registry());
  const before = JSON.stringify(value);
  const first = validateArchitectureRuleRegistry(
    value,
    WEB_ARCHITECTURE_DETECTORS,
    AVAILABLE_ENFORCEMENTS
  );
  const second = validateArchitectureRuleRegistry(
    value,
    WEB_ARCHITECTURE_DETECTORS,
    AVAILABLE_ENFORCEMENTS
  );
  assert.deepEqual(first, second);
  assert.equal(JSON.stringify(value), before);
  assert.equal(Object.getPrototypeOf(first), Object.prototype);
  assert.equal(Object.getPrototypeOf(first.errors), Array.prototype);
});

const TRUSTED_FIXTURE_FILES = Object.freeze({
  'package.json': '{"private":true}\n',
  'tools/check-web-change-budget.mjs': readFileSync(
    join(REPOSITORY_ROOT, 'tools/check-web-change-budget.mjs'),
    'utf8'
  ),
  'tools/web-architecture-detectors.mjs': readFileSync(
    join(REPOSITORY_ROOT, 'tools/web-architecture-detectors.mjs'),
    'utf8'
  ),
  'tools/web-architecture-evaluation.mjs': readFileSync(
    join(REPOSITORY_ROOT, 'tools/web-architecture-evaluation.mjs'),
    'utf8'
  ),
  'tools/web-change-policy.json': readFileSync(
    join(REPOSITORY_ROOT, 'tools/web-change-policy.json'),
    'utf8'
  ),
  'tools/web-change-source.mjs': readFileSync(
    join(REPOSITORY_ROOT, 'tools/web-change-source.mjs'),
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

const writeFixtureFile = (root, path, content) => {
  const target = join(root, path);
  mkdirSync(dirname(target), { recursive: true });
  writeFileSync(target, content, 'utf8');
};

const withTrustedCheckerRepository = (mutateHead, runCase, baseFiles = {}) => {
  const root = mkdtempSync(join(tmpdir(), 'gbdraw-architecture-ratchet-'));
  try {
    Object.entries({ ...TRUSTED_FIXTURE_FILES, ...baseFiles }).forEach(([path, content]) => {
      writeFixtureFile(root, path, content);
    });
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
    const result = spawnSync(process.execPath, [
      'tools/check-web-change-budget.mjs',
      '--base', base,
      '--head', head
    ], {
      cwd: root,
      encoding: 'utf8',
      stdio: ['ignore', 'pipe', 'pipe']
    });
    runCase({
      base,
      head,
      status: result.status,
      output: `${result.stdout || ''}${result.stderr || ''}`
    });
  } finally {
    rmSync(root, { recursive: true, force: true });
  }
};

const writeRules = (write, value) => write(
  'tools/web-architecture-rules.json',
  `${JSON.stringify(value, null, 2)}\n`
);

test('trusted checker admits an isolated clean-base candidate registry as inert data', () => {
  withTrustedCheckerRepository(
    (write) => writeRules(write, registry()),
    ({ status, output }) => {
      assert.equal(status, 0, output);
      assert.match(output, /Classification: EXPANSION/);
      assert.match(output, /canonical-path\.render-request against untouched base: CONFORMING/);
      assert.match(output, /semantic-owner\.render-request against untouched base: CONFORMING/);
      assert.match(output, /Result: \*\*PASS\*\*/);
    }
  );
});

test('trusted checker rejects malformed candidate JSON', () => {
  withTrustedCheckerRepository(
    (write) => write('tools/web-architecture-rules.json', '{"schemaVersion":1,"rules":['),
    ({ status, output }) => {
      assert.equal(status, 1, output);
      assert.match(output, /Cannot parse tools\/web-architecture-rules\.json/);
      assert.match(output, /Result: \*\*FAIL\*\*/);
    }
  );
});

test('trusted checker rejects an authority contraction that excludes active head source', () => {
  const baseRules = registry([semanticRule({
    allowedDefinitionPaths: [
      'services/future-session-request.js',
      'services/session-request.js'
    ]
  })]);
  const candidateRules = registry([semanticRule({
    allowedDefinitionPaths: ['services/future-session-request.js']
  })]);
  withTrustedCheckerRepository(
    (write) => writeRules(write, candidateRules),
    ({ status, output }) => {
      assert.equal(status, 1, output);
      assert.match(output, /Classification: CONTRACTION/);
      assert.match(output, /authority contraction or tightening is DIVERGENT on head source data/);
    },
    { 'tools/web-architecture-rules.json': `${JSON.stringify(baseRules, null, 2)}\n` }
  );
});

test('trusted checker never executes a proposed head detector during self-authorization', () => {
  withTrustedCheckerRepository(
    (write) => {
      writeRules(write, registry([semanticRule({
        allowedDefinitionPaths: ['app/secondary.js']
      })]));
      write(
        'gbdraw/web/js/app/secondary.js',
        'export const buildCanonicalRenderRequest = () => ({});\n'
      );
      write(
        'tools/web-architecture-detectors.mjs',
        'throw new Error("MALICIOUS HEAD DETECTOR EXECUTED");\n'
      );
    },
    ({ status, output }) => {
      assert.equal(status, 1, output);
      assert.doesNotMatch(output, /MALICIOUS HEAD DETECTOR EXECUTED/);
      assert.match(output, /claims a clean base but is DIVERGENT/);
      assert.match(output, /architecture rule authority changes must be isolated/);
      assert.match(output, /production runtime files and Web guard\/CI files changed together/);
    }
  );
});
