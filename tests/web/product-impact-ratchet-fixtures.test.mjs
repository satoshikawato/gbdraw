import assert from 'node:assert/strict';
import test from 'node:test';

import { WEB_ARCHITECTURE_DETECTORS } from '../../tools/web-architecture-detectors.mjs';
import { createArchitectureSubjectDelta } from '../../tools/web-architecture-evaluation.mjs';
import { parseCurrentProductImpactDecisions } from '../../tools/web-product-impact-decision-source.mjs';
import {
  evaluateProductImpact,
  validateProductDecisionAuthority,
  validateProductImpactMap
} from '../../tools/web-product-impact-evaluation.mjs';

const HEAD_SHA = 'a'.repeat(40);
const START_MARKER = '<!-- gbdraw-product-impact-decision:start -->';
const END_MARKER = '<!-- gbdraw-product-impact-decision:end -->';
const OWNER_CURRENT = 'services/session-request.js';
const OWNER_FUTURE = 'services/future-session-request.js';
const ENTRY_CURRENT = 'app/run-analysis.js -> services/session-request.js';
const ENTRY_FUTURE = 'app/alternate.js -> services/session-request.js';

const architectureRegistry = () => ({
  schemaVersion: 1,
  rules: [
    {
      key: 'canonical-path.render-request',
      kind: 'single-canonical-entry-edge',
      detector: 'canonical-path.render-request.v1',
      allowedEdges: [
        { from: 'app/alternate.js', to: 'services/session-request.js' },
        { from: 'app/run-analysis.js', to: 'services/session-request.js' }
      ],
      exactActiveEdgeCount: 1,
      enforcement: 'hard',
      baselineEligible: false
    },
    {
      key: 'semantic-owner.render-request',
      kind: 'single-semantic-owner',
      detector: 'semantic-owner.render-request.v1',
      allowedDefinitionPaths: [OWNER_FUTURE, OWNER_CURRENT],
      exactDefinitionCount: 1,
      enforcement: 'hard',
      baselineEligible: false
    }
  ]
});

const candidate = (ruleKey, subjectCategory, subject) => ({
  ruleKey,
  subjectCategory,
  subject
});
const entryCandidates = () => [
  candidate('canonical-path.render-request', 'canonical-entry-edge', ENTRY_FUTURE),
  candidate('canonical-path.render-request', 'canonical-entry-edge', ENTRY_CURRENT)
];
const ownerCandidates = () => [
  candidate('semantic-owner.render-request', 'definition-path', OWNER_FUTURE),
  candidate('semantic-owner.render-request', 'definition-path', OWNER_CURRENT)
];
const entryRequirement = (overrides = {}) => ({
  id: 'REQ-RENDER-ENTRY',
  category: 'SEMANTIC',
  checkpointRefs: ['JRN-RENDER-001:generate-submit'],
  anyOf: entryCandidates(),
  ...overrides
});
const ownerRequirement = (overrides = {}) => ({
  id: 'REQ-RENDER-OWNER',
  category: 'SEMANTIC',
  checkpointRefs: [
    'JRN-RENDER-001:generate-submit',
    'JRN-RENDER-001:roundtrip-regeneration'
  ],
  anyOf: ownerCandidates(),
  ...overrides
});
const contracts = () => [
  {
    ref: 'tests/web/request.test.mjs::canonical request projection',
    checkpointRefs: ['JRN-RENDER-001:generate-submit'],
    sensitivity: 'DIRECT_ASSERTION',
    execution: 'PR_GATE',
    residualRisk: 'Browser timing is reviewed separately.'
  },
  {
    ref: 'tests/web/session.test.mjs::request round trip',
    checkpointRefs: ['JRN-RENDER-001:roundtrip-regeneration'],
    sensitivity: 'MUTATION_FAIL',
    execution: 'PR_GATE',
    residualRisk: 'Future schema migrations require separate coverage.'
  }
];
const baseConcern = (overrides = {}) => ({
  key: 'product.render-request',
  scenarioRevision: 1,
  title: 'Canonical render request continuity',
  layer: 'STATE',
  journeys: [{
    id: 'JRN-RENDER-001',
    actor: 'Web user',
    context: 'A configured diagram is ready to generate.',
    goal: 'Generate and reproduce a diagram from canonical state.',
    steps: ['Configure a diagram', 'Select Generate', 'Reload and regenerate'],
    checkpoints: [
      { id: 'generate-submit', label: 'Generate request construction' },
      { id: 'roundtrip-regeneration', label: 'Regeneration after Session round trip' }
    ]
  }],
  effects: [
    {
      id: 'EFF-RENDER-CURRENT',
      checkpointRef: 'JRN-RENDER-001:generate-submit',
      statement: 'Generation uses the current canonical configuration.'
    },
    {
      id: 'EFF-RENDER-ROUNDTRIP',
      checkpointRef: 'JRN-RENDER-001:roundtrip-regeneration',
      statement: 'Reloaded state regenerates with the same request semantics.'
    }
  ],
  options: [{
    id: 'canonical-boundary',
    summary: 'Generation and replay share one canonical request boundary.',
    effectRefs: ['EFF-RENDER-CURRENT', 'EFF-RENDER-ROUNDTRIP'],
    requirements: [entryRequirement(), ownerRequirement()]
  }],
  resolution: {
    kind: 'authority-covered',
    selectedOptionId: 'canonical-boundary',
    authorityRefs: [
      'gbdraw/web/CLAUDE.md#request-and-session-boundary',
      'gbdraw/web/CLAUDE.md#runtime-data-flow'
    ]
  },
  contracts: contracts(),
  ...overrides
});
const productMap = (overrides = {}) => ({
  schemaVersion: 1,
  ruleCoverage: [
    {
      ruleKey: 'canonical-path.render-request',
      coverage: 'complete',
      enforcement: 'report-only'
    },
    {
      ruleKey: 'semantic-owner.render-request',
      coverage: 'complete',
      enforcement: 'report-only'
    }
  ],
  concerns: [baseConcern()],
  ...overrides
});
const decisionAuthority = (overrides = {}) => ({
  schemaVersion: 1,
  maintainerLogins: ['maintainer'],
  decisions: [],
  ...overrides
});
const validationCodes = (map, registry = architectureRegistry()) => (
  validateProductImpactMap(map, registry, WEB_ARCHITECTURE_DETECTORS)
    .errors.map(({ code }) => code)
);
const ownerDelta = (beforeSubjects, afterSubjects) => createArchitectureSubjectDelta({
  ruleKey: 'semantic-owner.render-request',
  kind: 'single-semantic-owner',
  detector: 'semantic-owner.render-request.v1',
  subjectCategory: 'definition-path',
  beforeSubjects,
  afterSubjects
});
const entryDelta = (beforeSubjects, afterSubjects) => createArchitectureSubjectDelta({
  ruleKey: 'canonical-path.render-request',
  kind: 'single-canonical-entry-edge',
  detector: 'canonical-path.render-request.v1',
  subjectCategory: 'canonical-entry-edge',
  beforeSubjects,
  afterSubjects
});
const ruleFacts = ({
  ownerBefore = [OWNER_CURRENT],
  ownerAfter = [OWNER_FUTURE],
  entryBefore = [ENTRY_CURRENT],
  entryAfter = [ENTRY_CURRENT]
} = {}) => [
  {
    ruleKey: 'canonical-path.render-request',
    subjectCategory: 'canonical-entry-edge',
    beforeSubjects: entryBefore,
    afterSubjects: entryAfter
  },
  {
    ruleKey: 'semantic-owner.render-request',
    subjectCategory: 'definition-path',
    beforeSubjects: ownerBefore,
    afterSubjects: ownerAfter
  }
];
const emptyCurrentDecisions = () => parseCurrentProductImpactDecisions({
  body: '',
  eventHeadSha: HEAD_SHA,
  eventAuthorLogin: 'maintainer',
  maxBodyBytes: 65536,
  maxDecisionCount: 8
});
const decisionBlock = (decisions, { headSha = HEAD_SHA } = {}) => (
  `${START_MARKER}\n${JSON.stringify({ schemaVersion: 1, headSha, decisions })}\n${END_MARKER}`
);
const currentDecision = (overrides = {}) => ({
  concernKey: 'product.render-request',
  scenarioRevision: 2,
  selectedOptionId: 'after-option',
  acceptedImpactClass: 'AFFORDANCE_PRESERVED',
  rationale: 'The internal request owner changes while the supported workflow remains intact.',
  acknowledgedChangedRequirementRefs: ['REQ-AFTER-OWNER', 'REQ-BEFORE-OWNER'],
  evidenceRefs: ['tests/web/request.test.mjs::owner transition'],
  residualRisk: 'No additional residual risk is accepted.',
  ...overrides
});
const parseDecisions = (decisions, overrides = {}) => parseCurrentProductImpactDecisions({
  body: decisionBlock(decisions, overrides),
  eventHeadSha: overrides.eventHeadSha || HEAD_SHA,
  eventAuthorLogin: overrides.eventAuthorLogin || 'maintainer',
  maxBodyBytes: overrides.maxBodyBytes || 65536,
  maxDecisionCount: overrides.maxDecisionCount || 8
});

const transitionConcern = ({ afterEffects = ['EFF-SHARED'], resolution = null } = {}) => ({
  key: 'product.render-request',
  scenarioRevision: 2,
  title: 'Equivalent render request realization',
  layer: 'STATE',
  journeys: [{
    id: 'JRN-RENDER-001',
    actor: 'Web user',
    context: 'Generation is requested.',
    goal: 'Keep the supported generation outcome.',
    steps: ['Select Generate', 'Inspect the Result'],
    checkpoints: [{ id: 'generate-submit', label: 'Generate submission' }]
  }],
  effects: [
    {
      id: 'EFF-CHANGED',
      checkpointRef: 'JRN-RENDER-001:generate-submit',
      statement: 'Generation exposes a materially different continuation.'
    },
    {
      id: 'EFF-SHARED',
      checkpointRef: 'JRN-RENDER-001:generate-submit',
      statement: 'Generation preserves the supported Result continuation.'
    }
  ],
  options: [
    {
      id: 'after-option',
      summary: 'The preauthorized future owner provides the request boundary.',
      effectRefs: afterEffects,
      requirements: [
        entryRequirement({ id: 'REQ-AFTER-ENTRY' }),
        ownerRequirement({
          id: 'REQ-AFTER-OWNER',
          checkpointRefs: ['JRN-RENDER-001:generate-submit'],
          anyOf: [candidate(
            'semantic-owner.render-request',
            'definition-path',
            OWNER_FUTURE
          )]
        })
      ]
    },
    {
      id: 'before-option',
      summary: 'The current owner provides the request boundary.',
      effectRefs: ['EFF-SHARED'],
      requirements: [
        entryRequirement({ id: 'REQ-BEFORE-ENTRY' }),
        ownerRequirement({
          id: 'REQ-BEFORE-OWNER',
          checkpointRefs: ['JRN-RENDER-001:generate-submit'],
          anyOf: [candidate(
            'semantic-owner.render-request',
            'definition-path',
            OWNER_CURRENT
          )]
        })
      ]
    }
  ],
  resolution: resolution || { kind: 'decision-required' },
  contracts: [{
    ref: 'tests/web/request.test.mjs::owner transition',
    checkpointRefs: ['JRN-RENDER-001:generate-submit'],
    sensitivity: 'DIRECT_ASSERTION',
    execution: 'PR_GATE',
    residualRisk: 'Timing remains outside this direct contract.'
  }]
});
const transitionMap = (options = {}) => productMap({
  concerns: [transitionConcern(options)]
});
const evaluate = ({
  map = productMap(),
  authority = decisionAuthority(),
  deltas = [ownerDelta([OWNER_CURRENT], [OWNER_FUTURE])],
  facts = ruleFacts(),
  current = emptyCurrentDecisions(),
  changedPaths = []
} = {}) => evaluateProductImpact({
  subjectDeltas: deltas,
  ruleFacts: facts,
  map,
  decisionAuthority: authority,
  currentDecisions: current,
  changedPaths
});

test('structured subject deltas normalize unchanged, addition, removal, and replacement', () => {
  assert.deepEqual(ownerDelta([OWNER_CURRENT, OWNER_CURRENT], [OWNER_CURRENT]), {
    ruleKey: 'semantic-owner.render-request',
    kind: 'single-semantic-owner',
    detector: 'semantic-owner.render-request.v1',
    subjectCategory: 'definition-path',
    beforeSubjects: [OWNER_CURRENT],
    addedSubjects: [],
    removedSubjects: [],
    afterSubjects: [OWNER_CURRENT],
    changed: false
  });
  assert.deepEqual(ownerDelta([], [OWNER_CURRENT]).addedSubjects, [OWNER_CURRENT]);
  assert.deepEqual(ownerDelta([OWNER_CURRENT], []).removedSubjects, [OWNER_CURRENT]);
  assert.deepEqual(
    ownerDelta([OWNER_CURRENT], [OWNER_FUTURE]),
    {
      ruleKey: 'semantic-owner.render-request',
      kind: 'single-semantic-owner',
      detector: 'semantic-owner.render-request.v1',
      subjectCategory: 'definition-path',
      beforeSubjects: [OWNER_CURRENT],
      addedSubjects: [OWNER_FUTURE],
      removedSubjects: [OWNER_CURRENT],
      afterSubjects: [OWNER_FUTURE],
      changed: true
    }
  );
});

test('structured subject deltas reject unknown fields and invalid canonical identifiers', () => {
  assert.throws(() => createArchitectureSubjectDelta({
    ruleKey: 'semantic-owner.render-request',
    kind: 'single-semantic-owner',
    detector: 'semantic-owner.render-request.v1',
    subjectCategory: 'definition-path',
    beforeSubjects: [],
    afterSubjects: [],
    extra: true
  }), /requires exactly/);
  assert.throws(() => ownerDelta([' padded '], []), /canonical nonempty subject strings/);
  assert.throws(() => createArchitectureSubjectDelta({
    ruleKey: 'INVALID',
    kind: 'single-semantic-owner',
    detector: 'semantic-owner.render-request.v1',
    subjectCategory: 'definition-path',
    beforeSubjects: [],
    afterSubjects: []
  }), /canonical rule key/);
});

test('structured subject delta inputs remain unchanged and outputs are deeply frozen', () => {
  const input = {
    ruleKey: 'semantic-owner.render-request',
    kind: 'single-semantic-owner',
    detector: 'semantic-owner.render-request.v1',
    subjectCategory: 'definition-path',
    beforeSubjects: [OWNER_CURRENT, OWNER_FUTURE],
    afterSubjects: [OWNER_CURRENT]
  };
  const before = structuredClone(input);
  const delta = createArchitectureSubjectDelta(input);
  assert.deepEqual(input, before);
  assert.equal(Object.isFrozen(delta), true);
  assert.equal(Object.isFrozen(delta.beforeSubjects), true);
  assert.throws(() => delta.beforeSubjects.push('x'));
});

test('the canonical Product Impact fixture validates', () => {
  assert.deepEqual(
    validateProductImpactMap(productMap(), architectureRegistry(), WEB_ARCHITECTURE_DETECTORS),
    { valid: true, errors: [] }
  );
  assert.deepEqual(validateProductDecisionAuthority(decisionAuthority(), productMap()), {
    valid: true,
    errors: []
  });
});

test('strict map validation rejects unknown fields at every schema level', () => {
  const mutations = [
    (map) => { map.extra = true; },
    (map) => { map.ruleCoverage[0].extra = true; },
    (map) => { map.concerns[0].extra = true; },
    (map) => { map.concerns[0].journeys[0].extra = true; },
    (map) => { map.concerns[0].journeys[0].checkpoints[0].extra = true; },
    (map) => { map.concerns[0].effects[0].extra = true; },
    (map) => { map.concerns[0].options[0].extra = true; },
    (map) => { map.concerns[0].options[0].requirements[0].extra = true; },
    (map) => { map.concerns[0].options[0].requirements[0].anyOf[0].extra = true; },
    (map) => { map.concerns[0].resolution.extra = true; },
    (map) => { map.concerns[0].contracts[0].extra = true; }
  ];
  mutations.forEach((mutate) => {
    const map = productMap();
    mutate(map);
    assert.ok(validationCodes(map).includes('unknown-field'));
  });
});

test('canonical IDs and references must be sorted and unique', () => {
  const duplicateEffects = productMap();
  duplicateEffects.concerns[0].options[0].effectRefs = [
    'EFF-RENDER-ROUNDTRIP',
    'EFF-RENDER-ROUNDTRIP'
  ];
  assert.ok(validationCodes(duplicateEffects).includes('duplicate-value'));

  const unsortedCoverage = productMap();
  unsortedCoverage.ruleCoverage.reverse();
  assert.ok(validationCodes(unsortedCoverage).includes('unsorted-array'));

  const unsortedRequirements = productMap();
  unsortedRequirements.concerns[0].options[0].requirements.reverse();
  assert.ok(validationCodes(unsortedRequirements).includes('unsorted-array'));
});

test('invalid architecture rule, subject category, and subject references fail closed', () => {
  const cases = [
    ['unknown-architecture-rule', (candidate_) => { candidate_.ruleKey = 'unknown.rule'; }],
    ['subject-category-mismatch', (candidate_) => { candidate_.subjectCategory = 'definition-path'; }],
    ['unapproved-architecture-subject', (candidate_) => { candidate_.subject = 'app/unknown.js'; }]
  ];
  cases.forEach(([code, mutate]) => {
    const map = productMap();
    mutate(map.concerns[0].options[0].requirements[0].anyOf[0]);
    assert.ok(validationCodes(map).includes(code), code);
  });
});

test('hard enforcement rejects partial coverage and manual-only checkpoint evidence', () => {
  const hardPartial = productMap();
  hardPartial.ruleCoverage[0] = {
    ...hardPartial.ruleCoverage[0],
    coverage: 'partial',
    enforcement: 'hard'
  };
  assert.ok(validationCodes(hardPartial).includes('hard-requires-complete-coverage'));

  const manualOnly = productMap();
  manualOnly.ruleCoverage.forEach((coverage) => { coverage.enforcement = 'hard'; });
  manualOnly.concerns[0].contracts.forEach((contract) => {
    contract.sensitivity = 'MANUAL_ACCEPTANCE';
    contract.execution = 'MANUAL';
  });
  const codes = validationCodes(manualOnly);
  assert.ok(codes.includes('missing-hard-contract-coverage'));
  assert.ok(codes.includes('manual-only-hard-concern'));
});

test('complete coverage rejects an unmapped allowed subject', () => {
  const map = productMap();
  map.concerns[0].options[0].requirements[0].anyOf = [
    candidate('canonical-path.render-request', 'canonical-entry-edge', ENTRY_CURRENT)
  ];
  assert.ok(validationCodes(map).includes('unmapped-allowed-subject'));
});

test('journey, checkpoint, effect, option, and contract references are cross-validated', () => {
  const cases = [
    ['unknown-checkpoint-ref', (map) => {
      map.concerns[0].effects[0].checkpointRef = 'JRN-RENDER-001:unknown';
    }],
    ['unknown-effect-ref', (map) => {
      map.concerns[0].options[0].effectRefs[0] = 'EFF-UNKNOWN';
    }],
    ['unknown-selected-option', (map) => {
      map.concerns[0].resolution.selectedOptionId = 'unknown-option';
    }],
    ['unknown-checkpoint-ref', (map) => {
      map.concerns[0].contracts[0].checkpointRefs = ['JRN-RENDER-001:unknown'];
    }]
  ];
  cases.forEach(([code, mutate]) => {
    const map = productMap();
    mutate(map);
    assert.ok(validationCodes(map).includes(code), code);
  });
});

test('AND-of-OR alternatives preserve one option when its active provider changes', () => {
  const results = evaluate();
  assert.equal(results.length, 1);
  assert.equal(results[0].impactClass, 'NO_USER_VISIBLE_DIFFERENCE');
  assert.equal(results[0].observation, 'CONFORMING');
  assert.deepEqual(results[0].beforeOptions, ['canonical-boundary']);
  assert.deepEqual(results[0].afterOptions, ['canonical-boundary']);
  assert.deepEqual(
    results[0].requirementResults.find(({ requirementRef }) => (
      requirementRef === 'REQ-RENDER-OWNER'
    )),
    {
      requirementRef: 'REQ-RENDER-OWNER',
      category: 'SEMANTIC',
      beforeSatisfied: true,
      afterSatisfied: true,
      beforeActiveSubjects: [OWNER_CURRENT],
      afterActiveSubjects: [OWNER_FUTURE],
      changed: true
    }
  );
});

test('separate jointly required contributions cannot be hidden by a shared option ID', () => {
  const results = evaluate({
    deltas: [entryDelta([ENTRY_CURRENT], [])],
    facts: ruleFacts({
      ownerBefore: [OWNER_CURRENT],
      ownerAfter: [OWNER_CURRENT],
      entryBefore: [ENTRY_CURRENT],
      entryAfter: []
    })
  });
  assert.equal(results[0].impactClass, 'RETIREMENT');
  assert.equal(results[0].observation, 'ORDINARY_REGRESSION');
  assert.deepEqual(results[0].beforeOptions, ['canonical-boundary']);
  assert.deepEqual(results[0].afterOptions, []);
  assert.equal(
    results[0].requirementResults.find(({ requirementRef }) => (
      requirementRef === 'REQ-RENDER-OWNER'
    )).afterSatisfied,
    true
  );
});

test('effect-equivalent option transitions require authority and retain stable effects', () => {
  const map = transitionMap();
  assert.equal(validationCodes(map).length, 0);
  const result = evaluate({ map })[0];
  assert.equal(result.impactClass, 'AFFORDANCE_PRESERVED');
  assert.equal(result.observation, 'UNRESOLVED_DECISION');
  assert.deepEqual(result.beforeEffects, ['EFF-SHARED']);
  assert.deepEqual(result.afterEffects, ['EFF-SHARED']);
});

test('a stable effect change or retirement is not eligible for a current decision', () => {
  const map = transitionMap({ afterEffects: ['EFF-CHANGED'] });
  const current = parseDecisions([currentDecision()]);
  assert.equal(current.valid, true);
  const result = evaluate({ map, current })[0];
  assert.equal(result.impactClass, 'RETIREMENT');
  assert.equal(result.observation, 'UNRESOLVED_DECISION');
  assert.match(result.currentDecisionIssues.join(' '), /not affordance-preserving/);
});

test('static and domain authority are represented separately', () => {
  const staticResult = evaluate()[0];
  assert.equal(staticResult.authorityResolution, 'STATIC_AUTHORITY');

  const domainMap = productMap();
  domainMap.concerns[0].resolution = {
    kind: 'domain-derived',
    selectedOptionId: 'canonical-boundary',
    authorityRefs: ['docs/domain.md#canonical-request-invariant']
  };
  assert.equal(validationCodes(domainMap).length, 0);
  assert.equal(evaluate({ map: domainMap })[0].authorityResolution, 'DOMAIN_AUTHORITY');
});

test('durable decisions validate active uniqueness, concern revision, and contracts', () => {
  const map = transitionMap();
  const durable = {
    id: 'BD-001',
    concernKey: 'product.render-request',
    scenarioRevision: 2,
    selectedOptionId: 'after-option',
    acceptedImpactClass: 'AFFORDANCE_PRESERVED',
    acknowledgedRetiredRequirementRefs: [],
    acceptanceCheckpointRefs: ['JRN-RENDER-001:generate-submit'],
    rationale: 'The future owner preserves the supported generation outcome.'
  };
  const authority = decisionAuthority({ decisions: [durable] });
  assert.deepEqual(validateProductDecisionAuthority(authority, map), {
    valid: true,
    errors: []
  });
  assert.equal(evaluate({ map, authority })[0].authorityResolution, 'DURABLE_DECISION');

  const stale = decisionAuthority({
    decisions: [{ ...durable, scenarioRevision: 1 }]
  });
  assert.ok(
    validateProductDecisionAuthority(stale, map).errors
      .some(({ code }) => code === 'stale-scenario-revision')
  );

  const duplicate = decisionAuthority({
    decisions: [durable, { ...durable, id: 'BD-002' }]
  });
  assert.ok(
    validateProductDecisionAuthority(duplicate, map).errors
      .some(({ code }) => code === 'duplicate-active-decision')
  );
});

test('durable decision authority rejects unknown fields and invalid retirement scope', () => {
  const map = transitionMap();
  const durable = {
    id: 'BD-001',
    concernKey: 'product.render-request',
    scenarioRevision: 2,
    selectedOptionId: 'after-option',
    acceptedImpactClass: 'AFFORDANCE_PRESERVED',
    acknowledgedRetiredRequirementRefs: [],
    acceptanceCheckpointRefs: ['JRN-RENDER-001:generate-submit'],
    rationale: 'The selected outcome preserves the supported generation continuation.'
  };
  const unknownTop = decisionAuthority({ extra: true });
  assert.ok(
    validateProductDecisionAuthority(unknownTop, map).errors
      .some(({ code }) => code === 'unknown-field')
  );
  const unknownDecision = decisionAuthority({ decisions: [{ ...durable, extra: true }] });
  assert.ok(
    validateProductDecisionAuthority(unknownDecision, map).errors
      .some(({ code }) => code === 'unknown-field')
  );
  const invalidRetirement = decisionAuthority({
    decisions: [{
      ...durable,
      acknowledgedRetiredRequirementRefs: ['REQ-BEFORE-OWNER']
    }]
  });
  assert.ok(
    validateProductDecisionAuthority(invalidRetirement, map).errors
      .some(({ code }) => code === 'retirement-impact-mismatch')
  );
  assert.ok(
    validateProductDecisionAuthority(
      decisionAuthority({ maintainerLogins: [] }),
      map
    ).errors.some(({ code }) => code === 'empty-maintainer-logins')
  );
});

test('durable decisions are forbidden for static authority and automatic no-difference', () => {
  const decision = {
    id: 'BD-001',
    concernKey: 'product.render-request',
    scenarioRevision: 1,
    selectedOptionId: 'canonical-boundary',
    acceptedImpactClass: 'PRODUCT_CHANGE',
    acknowledgedRetiredRequirementRefs: [],
    acceptanceCheckpointRefs: ['JRN-RENDER-001:generate-submit'],
    rationale: 'A durable choice is deliberately attempted against static authority.'
  };
  const validation = validateProductDecisionAuthority(
    decisionAuthority({ decisions: [decision] }),
    productMap()
  );
  assert.ok(validation.errors.some(({ code }) => code === 'decision-conflicts-with-resolution'));
  assert.ok(
    validateProductDecisionAuthority(
      decisionAuthority({ decisions: [{ ...decision, acceptedImpactClass: 'NO_USER_VISIBLE_DIFFERENCE' }] }),
      productMap()
    ).errors.some(({ code }) => code === 'invalid-accepted-impact-class')
  );
});

test('bounded parser accepts empty declarations without a head SHA', () => {
  const parsed = parseCurrentProductImpactDecisions({
    body: decisionBlock([], { headSha: '' }),
    eventHeadSha: '',
    eventAuthorLogin: 'Maintainer',
    maxBodyBytes: 4096,
    maxDecisionCount: 2
  });
  assert.deepEqual(parsed, {
    valid: true,
    present: true,
    metadata: { eventAuthorLogin: 'maintainer', eventHeadSha: '' },
    declaration: { schemaVersion: 1, headSha: '', decisions: [] },
    errors: []
  });
});

test('bounded parser normalizes exact-head decisions and product-level rationale', () => {
  const parsed = parseCurrentProductImpactDecisions({
    body: decisionBlock([currentDecision({ rationale: '  Preserve   the supported outcome.  ' })], {
      headSha: HEAD_SHA.toUpperCase()
    }),
    eventHeadSha: HEAD_SHA,
    eventAuthorLogin: 'Maintainer',
    maxBodyBytes: 65536,
    maxDecisionCount: 8
  });
  assert.equal(parsed.valid, true);
  assert.equal(parsed.declaration.headSha, HEAD_SHA);
  assert.equal(parsed.declaration.decisions[0].rationale, 'Preserve the supported outcome.');
  assert.equal(parsed.metadata.eventAuthorLogin, 'maintainer');
});

test('bounded parser rejects stale heads, wrong impact class, and incomplete rationale', () => {
  const stale = parseDecisions([currentDecision()], { eventHeadSha: 'b'.repeat(40) });
  assert.ok(stale.errors.some(({ code }) => code === 'stale-head-sha'));

  const wrongImpact = parseDecisions([currentDecision({ acceptedImpactClass: 'PRODUCT_CHANGE' })]);
  assert.ok(wrongImpact.errors.some(({ code }) => code === 'invalid-current-impact-class'));

  const emptyRationale = parseDecisions([currentDecision({ rationale: '   ' })]);
  assert.ok(emptyRationale.errors.some(({ code }) => code === 'invalid-string'));
});

test('bounded parser rejects unknown fields and bounded collection overflows', () => {
  const unknownTopBody = `${START_MARKER}\n${JSON.stringify({
    schemaVersion: 1,
    headSha: '',
    decisions: [],
    extra: true
  })}\n${END_MARKER}`;
  const unknownTop = parseCurrentProductImpactDecisions({
    body: unknownTopBody,
    eventHeadSha: HEAD_SHA,
    eventAuthorLogin: 'maintainer',
    maxBodyBytes: 65536,
    maxDecisionCount: 8
  });
  assert.ok(unknownTop.errors.some(({ code }) => code === 'unknown-field'));

  const unknownDecision = parseDecisions([currentDecision({ extra: true })]);
  assert.ok(unknownDecision.errors.some(({ code }) => code === 'unknown-field'));

  const tooManyDecisions = parseDecisions(
    [currentDecision(), currentDecision({ concernKey: 'product.second' })],
    { maxDecisionCount: 1 }
  );
  assert.ok(tooManyDecisions.errors.some(({ code }) => code === 'too-many-decisions'));

  const oversizedRefs = parseDecisions([currentDecision({
    evidenceRefs: Array.from({ length: 33 }, (_, index) => `evidence-${String(index).padStart(2, '0')}`)
  })]);
  assert.ok(oversizedRefs.errors.some(({ code }) => code === 'oversized-array'));
});

test('bounded parser rejects unmatched, reversed, duplicate, malformed, and oversized blocks', () => {
  const inputs = [
    START_MARKER,
    `${END_MARKER}\n{}\n${START_MARKER}`,
    `${decisionBlock([])}\n${decisionBlock([])}`,
    `${START_MARKER}\nnot-json\n${END_MARKER}`
  ];
  inputs.forEach((body) => {
    const parsed = parseCurrentProductImpactDecisions({
      body,
      eventHeadSha: HEAD_SHA,
      eventAuthorLogin: 'maintainer',
      maxBodyBytes: 65536,
      maxDecisionCount: 8
    });
    assert.equal(parsed.valid, false, body.slice(0, 40));
  });
  const oversized = parseCurrentProductImpactDecisions({
    body: 'x'.repeat(101),
    eventHeadSha: HEAD_SHA,
    eventAuthorLogin: 'maintainer',
    maxBodyBytes: 100,
    maxDecisionCount: 8
  });
  assert.ok(oversized.errors.some(({ code }) => code === 'oversized-body'));
});

test('parser errors are sanitized and never expose the untrusted body', () => {
  const secret = 'PRIVATE-SEQUENCE-DO-NOT-ECHO';
  const parsed = parseCurrentProductImpactDecisions({
    body: `${START_MARKER}\n${secret}\n${END_MARKER}`,
    eventHeadSha: HEAD_SHA,
    eventAuthorLogin: 'maintainer',
    maxBodyBytes: 65536,
    maxDecisionCount: 8
  });
  assert.equal(parsed.valid, false);
  assert.doesNotMatch(JSON.stringify(parsed.errors), new RegExp(secret));
});

test('free-form PRODUCT_DECISION text is never parsed as checker authority', () => {
  const human = [
    'PRODUCT_DECISION',
    'Concern: product.render-request',
    'Choice: after-option',
    'Rationale: preserve behavior'
  ].join('\n');
  const outsideMarkers = parseCurrentProductImpactDecisions({
    body: human,
    eventHeadSha: HEAD_SHA,
    eventAuthorLogin: 'maintainer',
    maxBodyBytes: 65536,
    maxDecisionCount: 8
  });
  assert.equal(outsideMarkers.valid, true);
  assert.equal(outsideMarkers.present, false);
  assert.deepEqual(outsideMarkers.declaration.decisions, []);

  const insideMarkers = parseCurrentProductImpactDecisions({
    body: `${START_MARKER}\n${human}\n${END_MARKER}`,
    eventHeadSha: HEAD_SHA,
    eventAuthorLogin: 'maintainer',
    maxBodyBytes: 65536,
    maxDecisionCount: 8
  });
  assert.equal(insideMarkers.valid, false);
  assert.ok(insideMarkers.errors.some(({ code }) => code === 'invalid-json'));
});

test('an exact-head current maintainer decision resolves an effect-equivalent transition', () => {
  const map = transitionMap();
  const current = parseDecisions([currentDecision()]);
  assert.equal(current.valid, true);
  const result = evaluate({ map, current })[0];
  assert.equal(result.observation, 'CONFORMING');
  assert.equal(result.authorityResolution, 'CURRENT_MAINTAINER_DECISION');
  assert.equal(result.selectedOptionId, 'after-option');
  assert.equal(result.decisionRationale, currentDecision().rationale);
});

test('non-maintainers and incomplete changed-requirement acknowledgements are ineligible', () => {
  const map = transitionMap();
  const nonMaintainer = parseDecisions([currentDecision()], {
    eventAuthorLogin: 'contributor'
  });
  const first = evaluate({ map, current: nonMaintainer })[0];
  assert.equal(first.observation, 'UNRESOLVED_DECISION');
  assert.match(first.currentDecisionIssues.join(' '), /not an eligible Product Decision Owner/);

  const incomplete = parseDecisions([currentDecision({
    acknowledgedChangedRequirementRefs: ['REQ-AFTER-OWNER']
  })]);
  const second = evaluate({ map, current: incomplete })[0];
  assert.equal(second.observation, 'UNRESOLVED_DECISION');
  assert.match(second.currentDecisionIssues.join(' '), /acknowledgement is incomplete/);
});

test('incompatible durable and current selections produce authority conflict', () => {
  const map = transitionMap();
  const durable = decisionAuthority({
    decisions: [{
      id: 'BD-001',
      concernKey: 'product.render-request',
      scenarioRevision: 2,
      selectedOptionId: 'after-option',
      acceptedImpactClass: 'AFFORDANCE_PRESERVED',
      acknowledgedRetiredRequirementRefs: [],
      acceptanceCheckpointRefs: ['JRN-RENDER-001:generate-submit'],
      rationale: 'The future owner is the durable selected outcome.'
    }]
  });
  const current = parseDecisions([currentDecision({ selectedOptionId: 'before-option' })]);
  const result = evaluate({ map, authority: durable, current })[0];
  assert.equal(result.observation, 'AUTHORITY_CONFLICT');
  assert.equal(result.authorityResolution, 'CONFLICT');
  assert.equal(result.selectedOptionId, null);
});

test('a current decision incompatible with static authority also produces conflict', () => {
  const map = transitionMap({
    resolution: {
      kind: 'authority-covered',
      selectedOptionId: 'before-option',
      authorityRefs: ['gbdraw/web/CLAUDE.md#runtime-data-flow']
    }
  });
  const current = parseDecisions([currentDecision()]);
  const result = evaluate({ map, current })[0];
  assert.equal(result.observation, 'AUTHORITY_CONFLICT');
  assert.equal(result.authorityResolution, 'CONFLICT');
});

test('a selected static option missing after the change is an ordinary regression', () => {
  const result = evaluate({
    deltas: [entryDelta([ENTRY_CURRENT], [])],
    facts: ruleFacts({
      ownerBefore: [OWNER_CURRENT],
      ownerAfter: [OWNER_CURRENT],
      entryBefore: [ENTRY_CURRENT],
      entryAfter: []
    })
  })[0];
  assert.equal(result.observation, 'ORDINARY_REGRESSION');
  assert.equal(result.selectedOptionId, 'canonical-boundary');
});

test('missing complete inventories and partial mappings yield insufficient evidence', () => {
  const partial = productMap();
  partial.ruleCoverage[1].coverage = 'partial';
  const partialResult = evaluate({ map: partial })[0];
  assert.equal(partialResult.observation, 'INSUFFICIENT_EVIDENCE');

  const missingFactResult = evaluate({ facts: ruleFacts().slice(0, 1) })[0];
  assert.equal(missingFactResult.observation, 'INSUFFICIENT_EVIDENCE');
  assert.match(missingFactResult.evidenceGaps.join(' '), /Missing complete/);
});

test('candidate-modified mapped contracts cannot be sole hard evidence', () => {
  const map = productMap();
  map.ruleCoverage.forEach((coverage) => { coverage.enforcement = 'hard'; });
  assert.equal(validationCodes(map).length, 0);
  const result = evaluate({
    map,
    changedPaths: ['tests/web/request.test.mjs']
  })[0];
  assert.equal(result.observation, 'INSUFFICIENT_EVIDENCE');
  assert.equal(
    result.contractResults.find(({ ref }) => ref.startsWith('tests/web/request.test.mjs'))
      .candidateModified,
    true
  );
});

test('one concern triggered by multiple rules is evaluated exactly once', () => {
  const results = evaluate({
    deltas: [
      entryDelta([ENTRY_CURRENT], [ENTRY_FUTURE]),
      ownerDelta([OWNER_CURRENT], [OWNER_FUTURE])
    ],
    facts: ruleFacts({ entryAfter: [ENTRY_FUTURE] })
  });
  assert.equal(results.length, 1);
  assert.deepEqual(results[0].triggeringRuleKeys, [
    'canonical-path.render-request',
    'semantic-owner.render-request'
  ]);
  assert.equal(results[0].impactClass, 'NO_USER_VISIBLE_DIFFERENCE');
});

test('unchanged registered subjects with no current decision are not applicable', () => {
  const results = evaluate({
    deltas: [ownerDelta([OWNER_CURRENT], [OWNER_CURRENT])],
    facts: ruleFacts({
      ownerBefore: [OWNER_CURRENT],
      ownerAfter: [OWNER_CURRENT]
    })
  });
  assert.deepEqual(results, []);
});

test('evaluation ordering is deterministic, inputs remain unchanged, and outputs are deeply frozen', () => {
  const input = {
    subjectDeltas: [ownerDelta([OWNER_CURRENT], [OWNER_FUTURE])],
    ruleFacts: ruleFacts(),
    map: productMap(),
    decisionAuthority: decisionAuthority(),
    currentDecisions: emptyCurrentDecisions(),
    changedPaths: []
  };
  const before = structuredClone(input);
  const first = evaluateProductImpact(input);
  const second = evaluateProductImpact(input);
  assert.deepEqual(first, second);
  assert.deepEqual(input, before);
  assert.equal(Object.isFrozen(first), true);
  assert.equal(Object.isFrozen(first[0]), true);
  assert.equal(Object.isFrozen(first[0].requirementResults), true);
  assert.throws(() => first[0].requirementResults.push({}));
});

test('validation results are deterministic deeply frozen plain data', () => {
  const map = productMap();
  const first = validateProductImpactMap(map, architectureRegistry(), WEB_ARCHITECTURE_DETECTORS);
  const second = validateProductImpactMap(map, architectureRegistry(), WEB_ARCHITECTURE_DETECTORS);
  assert.deepEqual(first, second);
  assert.equal(Object.getPrototypeOf(first), Object.prototype);
  assert.equal(Object.isFrozen(first), true);
  assert.equal(Object.isFrozen(first.errors), true);
});
