import assert from 'node:assert/strict';
import test from 'node:test';

import { WEB_ARCHITECTURE_DETECTORS } from '../../tools/web-architecture-detectors.mjs';
import {
  classifyArchitectureAuthorityDelta,
  classifyArchitectureRuleObservation,
  evaluateArchitectureRuleResult,
  validateArchitectureRuleRegistry,
  WEB_ARCHITECTURE_RULE_SCHEMA
} from '../../tools/web-architecture-evaluation.mjs';

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

test('hard and report-only decisions use the five-field result envelope', () => {
  const cases = [
    ['HARD', 'CONFORMING', 'PASS'],
    ['HARD', 'DIVERGENT', 'FAIL'],
    ['HARD', 'ABSENT_REQUIRED', 'FAIL'],
    ['REPORT_ONLY', 'CONFORMING', 'PASS'],
    ['REPORT_ONLY', 'DIVERGENT', 'REPORT'],
    ['REPORT_ONLY', 'ABSENT_REQUIRED', 'REPORT']
  ];

  cases.forEach(([mode, observation, decision]) => {
    const semanticInputs = Object.freeze({
      observation,
      mode,
      baselineRelation: 'NOT_APPLICABLE',
      authorityResolution: 'NOT_APPLICABLE'
    });
    assert.deepEqual(evaluateArchitectureRuleResult(semanticInputs), {
      ...semanticInputs,
      decision
    });
  });
});

test('hard and report-only evaluation rejects accepted-store relations', () => {
  for (const mode of ['HARD', 'REPORT_ONLY']) {
    for (const baselineRelation of ['ACCEPTED', 'NEW', 'FIXED']) {
      assert.throws(() => evaluateArchitectureRuleResult({
        observation: 'DIVERGENT',
        mode,
        baselineRelation,
        authorityResolution: 'NOT_APPLICABLE'
      }), /requires baselineRelation NOT_APPLICABLE/);
    }
    for (const authorityResolution of ['RETAINED', 'EXACT_CONTRACTION', 'INVALID_CHANGE']) {
      assert.throws(() => evaluateArchitectureRuleResult({
        observation: 'DIVERGENT',
        mode,
        baselineRelation: 'NOT_APPLICABLE',
        authorityResolution
      }), /requires authorityResolution NOT_APPLICABLE/);
    }
  }
});

test('malformed four-input evaluation combinations fail closed', () => {
  const valid = {
    observation: 'CONFORMING',
    mode: 'HARD',
    baselineRelation: 'NOT_APPLICABLE',
    authorityResolution: 'NOT_APPLICABLE'
  };
  const unsupported = [
    ['observation', 'UNKNOWN_OBSERVATION'],
    ['mode', 'UNKNOWN_MODE'],
    ['baselineRelation', 'UNKNOWN_BASELINE'],
    ['authorityResolution', 'UNKNOWN_AUTHORITY']
  ];

  unsupported.forEach(([field, value]) => {
    assert.throws(
      () => evaluateArchitectureRuleResult({ ...valid, [field]: value }),
      new RegExp(`Unsupported architecture rule ${field}`)
    );
  });
  assert.throws(
    () => evaluateArchitectureRuleResult({ ...valid, mode: 'FROZEN' }),
    /mode FROZEN is unavailable/
  );
  assert.throws(
    () => evaluateArchitectureRuleResult({ ...valid, decision: 'PASS' }),
    /requires exactly/
  );
  const { authorityResolution: _omitted, ...missingField } = valid;
  assert.throws(() => evaluateArchitectureRuleResult(missingField), /requires exactly/);
  assert.throws(() => evaluateArchitectureRuleResult(null), /semantic input object/);
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
