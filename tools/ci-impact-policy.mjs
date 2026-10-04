const freeze = (value) => Object.freeze(value);

export const IMPACT_PLAN_SCHEMA_VERSION = 3;

// Ordered only for the summary label. Routing unions every affected capability.
// `none` marks an identical tree; no path classifies as `none`.
export const IMPACT_CLASSES = freeze([
  'none', 'metadata', 'documentation', 'policy-documentation', 'tests-only', 'python-core', 'renderer',
  'web-runtime', 'session-persistence', 'gallery', 'losat-integration',
  'packaging', 'ci-only', 'full'
]);
export const IMPACT_DECISIONS = freeze(['selective', 'full']);
export const IMPACT_PLAN_BASES = freeze([
  'FULL_CHANGE', 'UNKNOWN_OR_INVALID_CHANGE', 'MANUAL_FULL_RUN',
  'ARCHITECTURE_CHANGE', 'DOCUMENTATION_ONLY_PR', 'LIGHT_CHANGE_WITH_DIRECT_BASE_EVIDENCE',
  'LIGHT_CHANGE_WITH_DIRECT_PARENT_EVIDENCE', 'INHERITED_EVIDENCE_UNAVAILABLE',
  'IDENTICAL_TREE_WITH_DIRECT_PARENT_EVIDENCE', 'LEAF_TEST_CHANGE_WITH_DIRECT_PARENT_EVIDENCE'
]);
const IDENTICAL_TREE_BASIS = 'IDENTICAL_TREE_WITH_DIRECT_PARENT_EVIDENCE';
const LEAF_TEST_BASIS = 'LEAF_TEST_CHANGE_WITH_DIRECT_PARENT_EVIDENCE';
const PARENT_EVIDENCE_BASES = freeze([
  'LIGHT_CHANGE_WITH_DIRECT_PARENT_EVIDENCE', IDENTICAL_TREE_BASIS, LEAF_TEST_BASIS
]);

// The full PR tier. Functional Playwright joins it only through the capabilities below.
const PR_FULL_TIER_JOBS = freeze([
  'web-change-budget', 'core-pr', 'recipes-standard', 'gallery', 'lint',
  'web-contracts-pr', 'web-pr-smoke'
]);
const PR_JOBS = freeze([...PR_FULL_TIER_JOBS, 'playwright-functional']);
const DEV_JOBS = freeze([
  'web-change-budget', 'core', 'recipes-standard', 'gallery', 'browser',
  'playwright-functional', 'playwright-performance', 'lint',
  'losat-cache-browser-acceptance'
]);
const PROFILE_REQUIRED_JOBS = freeze({
  pr: PR_JOBS,
  dev: DEV_JOBS,
  release: freeze([...DEV_JOBS, 'acceptance-supported-main', 'slow-main', 'vibrio-generate-release']),
  gallery: freeze(['browser', 'performance'])
});
const PR_CAPABILITY_JOBS = freeze({
  none: freeze([]),
  metadata: freeze([]),
  documentation: freeze(['recipes-standard']),
  'policy-documentation': freeze(['web-change-budget']),
  'tests-only': PR_JOBS,
  'python-core': freeze(['web-change-budget', 'core-pr', 'gallery', 'lint', 'web-contracts-pr', 'web-pr-smoke']),
  renderer: freeze(['web-change-budget', 'core-pr', 'gallery', 'lint', 'web-contracts-pr', 'web-pr-smoke']),
  'web-runtime': freeze(['web-change-budget', 'gallery', 'web-contracts-pr', 'web-pr-smoke', 'playwright-functional']),
  'session-persistence': freeze(['web-change-budget', 'core-pr', 'recipes-standard', 'gallery', 'lint', 'web-contracts-pr', 'web-pr-smoke', 'playwright-functional']),
  gallery: freeze(['web-change-budget', 'gallery', 'web-contracts-pr', 'web-pr-smoke', 'playwright-functional']),
  'losat-integration': freeze(['web-change-budget', 'core-pr', 'gallery', 'lint', 'web-contracts-pr', 'web-pr-smoke', 'playwright-functional']),
  packaging: PR_FULL_TIER_JOBS,
  'ci-only': PR_FULL_TIER_JOBS,
  full: PR_JOBS
});

// Runtime changes always receive comprehensive integrated-dev/Gallery validation.
// Control-plane, dependency, and unknown changes cannot inherit a narrower route.
const LIGHT_CLASSES = freeze(['metadata', 'documentation', 'policy-documentation']);
// Leaf tests carry the capability of the subject they exercise (or tests-only).
const LEAF_SUBJECT_CLASSES = freeze(['tests-only', 'web-runtime', 'session-persistence', 'gallery', 'losat-integration']);
export const requiresFullCoverage = (profile, capabilities) => profile === 'release'
  || capabilities.some((capability) => ['full', 'ci-only', 'packaging', 'tests-only'].includes(capability))
  || (profile !== 'pr' && capabilities.some((capability) => !LIGHT_CLASSES.includes(capability)));
export const isDocumentationOnly = (capabilities) => Array.isArray(capabilities)
  && capabilities.some((capability) => ['documentation', 'policy-documentation'].includes(capability))
  && capabilities.every((capability) => LIGHT_CLASSES.includes(capability));
const orderedCapabilities = (capabilities) => IMPACT_CLASSES.filter((capability) => capabilities.includes(capability));
const primaryImpact = (capabilities) => capabilities.at(-1);
const FULL_OBJECT_ID = /^(?:[0-9a-f]{40}|[0-9a-f]{64})$/i;
const ROOT_MARKDOWN = /^[^/]+\.md$/;
const POLICY_DOCUMENTS = new Set([
  'docs/internal/ARCHITECTURE_FITNESS_FUNCTION_RATCHET.md',
  'docs/internal/OPTION_INTEGRITY_PRODUCT_CONTRACT.md',
  'docs/internal/PRODUCT_IMPACT_RATCHET.md',
  'docs/internal/SELECTIVE_CI.md',
  'docs/internal/WEB_CHANGE_POLICY.md'
]);
const METADATA_DIRECTORIES = freeze(['.agents', '.claude', '.codex', '.cursor']);
const METADATA_FILES = new Set([
  '.github/pull_request_template.md',
  '.dockerignore',
  '.gitattributes',
  '.gitignore',
  'CITATION.cff',
  'LICENSE.txt',
  'LICENSE_LIBERATION_FONTS.txt'
]);
const PLAN_KEYS = freeze([
  'schemaVersion',
  'profile',
  'impact',
  'capabilities',
  'decision',
  'basis',
  'changeBaseSha',
  'changeHeadSha',
  'workflowSha',
  'requiredJobs',
  'changedPathCount',
  'inheritedEvidence',
  'leafTests'
]);
const LEAF_TEST_KEYS = freeze(['path', 'kind', 'jobs']);
const EVIDENCE_KEYS = freeze([
  'workflowPath',
  'aggregateName',
  'headSha',
  'runId',
  'aggregateJobId',
  'runUrl',
  'aggregateJobUrl'
]);

const isPlainObject = (value) => value !== null
  && typeof value === 'object'
  && !Array.isArray(value)
  && (Object.getPrototypeOf(value) === Object.prototype
    || Object.getPrototypeOf(value) === null);

const sameKeys = (value, expectedKeys) => isPlainObject(value)
  && Object.keys(value).length === expectedKeys.length
  && expectedKeys.every((key) => Object.hasOwn(value, key));

const fail = (code, message, details = {}) => {
  throw Object.assign(new Error(message), {
    name: 'CiImpactPolicyError',
    code,
    details
  });
};

const assertProfile = (profile) => {
  if (!Object.hasOwn(PROFILE_REQUIRED_JOBS, profile)) {
    fail('UNKNOWN_PROFILE', 'CI impact profile is not supported.', { profile });
  }
};

const assertImpact = (impact) => {
  if (!IMPACT_CLASSES.includes(impact)) {
    fail('UNKNOWN_IMPACT', 'CI impact class is not supported.', { impact });
  }
};

export const isFullObjectId = (value) => typeof value === 'string'
  && FULL_OBJECT_ID.test(value);

const isValidRepositoryPath = (path) => typeof path === 'string'
  && path.length > 0
  && !path.startsWith('/')
  && !path.endsWith('/')
  && !path.includes('\\')
  && !path.includes('\0')
  && path.split('/').every((part) => part && part !== '.' && part !== '..');

export const classifyPath = (path) => {
  const classified = (impact, reason) => freeze({ path, impact, reason });
  if (!isValidRepositoryPath(path)) return classified('full', 'INVALID_REPOSITORY_PATH');
  if (POLICY_DOCUMENTS.has(path)) return classified('policy-documentation', 'POLICY_DOCUMENT');
  if (path.startsWith('.github/workflows/') || path.startsWith('tests/ci/')
      || /^tools\/(?:ci-impact|check-web|web-(?:architecture|product|change)|check-promotion)/.test(path)
      || /^tests\/web\/(?:architecture|product-impact|promotion-readiness).*\.test\.mjs$/.test(path)
      || /^playwright.*\.config\.js$/.test(path)) {
    return classified('ci-only', 'CI_CONTROL_PLANE');
  }
  if (METADATA_FILES.has(path)) return classified('metadata', 'METADATA_FILE_ALLOWLIST');
  if (METADATA_DIRECTORIES.some((directory) => path.startsWith(`${directory}/`))) {
    return classified('metadata', 'METADATA_DIRECTORY_ALLOWLIST');
  }
  if (ROOT_MARKDOWN.test(path)) return classified('documentation', 'ROOT_MARKDOWN');
  if (path.startsWith('docs/')) return classified('documentation', 'DOCUMENTATION_TREE');
  if (['pyproject.toml', 'setup.py', 'MANIFEST.in', 'package.json', 'package-lock.json',
    'wrangler.toml', 'gbdraw/_build_support.py'].includes(path)
      || path.startsWith('recipe/') || path.startsWith('gbdraw/web/vendor/')
      || /^tools\/(?:prepare_browser_wheel|prepare_cloudflare_pages|verify_gui_offline|check_release_source)\.py$/.test(path)) {
    return classified('packaging', 'PACKAGE_OR_DEPENDENCY_INPUT');
  }
  if (path.startsWith('gbdraw/web/wasm/') || /^gbdraw\/web\/js\/.*losat[^/]*\.js$/.test(path)
      || /^gbdraw\/(?:comparisons|losat)(?:\/|\.)/.test(path)) {
    return classified('losat-integration', 'LOSAT_OWNER');
  }
  if (path.startsWith('gbdraw/web/gallery/')
      || /^tools\/(?:build_web_gallery|prepare_interactive_gallery_assets|refresh_gallery_sessions)\.py$/.test(path)) {
    return classified('gallery', 'GALLERY_INPUT');
  }
  if (/^gbdraw\/(?:session[^/]*|api\/session)\.py$/.test(path)
      || /^gbdraw\/web\/js\/(?:services\/(?:session[^/]*|config)|app\/(?:session[^/]*|history[^/]*))\.js$/.test(path)) {
    return classified('session-persistence', 'SESSION_OWNER');
  }
  if (path === 'gbdraw/web/index.html' || /^gbdraw\/web\/js\/.*\.js$/.test(path)) {
    return classified('web-runtime', 'WEB_SOURCE');
  }
  if (/^gbdraw\/(?:render|svg|diagrams|canvas|features|labels|legend|layout|tracks)\/.*\.py$/.test(path)) {
    return classified('renderer', 'RENDERER_SOURCE');
  }
  if (/^gbdraw\/(?:api|config|configurators|core|io)\/.*\.py$/.test(path)
      || /^gbdraw\/(?:__init__|cli|circular|linear)\.py$/.test(path)
      || path.startsWith('gbdraw/data/')) {
    return classified('python-core', 'PYTHON_SOURCE');
  }
  if (/^tests\/web\/.*(?:\.test\.mjs|\.playwright\.spec\.js|\.cjs)$/.test(path)) {
    if (/\/(?:session|current-session|history)[^/]*\./.test(path)) return classified('session-persistence', 'SESSION_TEST');
    if (/\/losat[^/]*\./.test(path)) return classified('losat-integration', 'LOSAT_TEST');
    if (/\/gallery[^/]*\./.test(path)) return classified('gallery', 'GALLERY_TEST');
    return classified('web-runtime', 'WEB_TEST');
  }
  if (path.startsWith('tests/')) return classified('tests-only', 'SHARED_OR_UNCLASSIFIED_TEST');
  return classified('full', 'FULL_BY_DEFAULT');
};

// Specs that only a dedicated Playwright configuration runs (SELECTIVE_CI.md "Leaf tests").
export const GALLERY_PARITY_SPEC = 'tests/web/contracts/gallery-publication-parity.serial.spec.js';
export const VIBRIO_FULL_GENERATION_SPEC = 'tests/web/contracts/vibrio-full-generation.serial.spec.js';
export const LEAF_TEST_KINDS = freeze(['python', 'node', 'functional', 'performance', 'gallery-parity', 'release-only']);
const PYTEST_JOBS = freeze(['core', 'recipes-standard', 'gallery', 'browser']);
const PYTEST_MARKER_JOBS = freeze({ recipe: 'recipes-standard', gallery: 'gallery', browser: 'browser' });
const LEAF_KIND_JOBS = freeze({
  dev: freeze({
    node: freeze(['browser']),
    functional: freeze(['playwright-functional']),
    performance: freeze(['playwright-performance']),
    'gallery-parity': freeze([]),
    'release-only': freeze([])
  }),
  gallery: freeze({
    node: freeze([]),
    functional: freeze([]),
    performance: freeze([]),
    'gallery-parity': freeze(['browser']),
    'release-only': freeze([])
  })
});

// The runner kind of a path that may be a leaf test, from its path alone. A functional spec
// must also be listed in tests/ci/functional-shards.json, and no other file may name it.
export const leafTestKind = (path) => {
  if (!isValidRepositoryPath(path) || classifyPath(path).impact === 'ci-only') return null;
  if (/^tests\/test_[^/]+\.py$/.test(path)) return 'python';
  if (path === GALLERY_PARITY_SPEC) return 'gallery-parity';
  if (path === VIBRIO_FULL_GENERATION_SPEC) return 'release-only';
  if (/^tests\/web\/(?:[^/]+\/)*[^/]*performance\.playwright\.spec\.js$/.test(path)) return 'performance';
  if (/^tests\/web\/(?:[^/]+\/)*[^/]+\.playwright\.spec\.js$/.test(path)) return 'functional';
  if (/^tests\/web\/(?:[^/]+\/)*[^/]+\.test\.mjs$/.test(path)) return 'node';
  return null;
};

// `markers` lists the selection markers named in a Python test file, or is null when they
// cannot be read; null selects every pytest job of the dev tier.
export const leafJobsFor = ({ profile, kind, markers = null }) => {
  if (!['dev', 'gallery'].includes(profile) || !LEAF_TEST_KINDS.includes(kind)) {
    fail('INVALID_LEAF_TESTS', 'Leaf tests route only dev and Gallery pushes of known kinds.', { profile, kind });
  }
  if (kind !== 'python') return LEAF_KIND_JOBS[profile][kind];
  if (profile === 'gallery') return freeze([]);
  if (markers === null) return PYTEST_JOBS;
  const selected = ['core', ...markers.map((marker) => PYTEST_MARKER_JOBS[marker]).filter(Boolean)];
  return freeze(PYTEST_JOBS.filter((job) => selected.includes(job)));
};

const validScoredStatus = (status) => {
  const match = status.match(/^([RC])(\d{1,3})$/);
  return Boolean(match) && Number(match[2]) <= 100;
};

const pathsForChange = (change) => {
  if (!isPlainObject(change) || typeof change.status !== 'string'
      || !Array.isArray(change.paths)) {
    return null;
  }
  if (['A', 'M', 'D'].includes(change.status) && change.paths.length === 1) {
    return change.paths;
  }
  if (validScoredStatus(change.status) && change.paths.length === 2) {
    return change.paths;
  }
  return null;
};

export const classifyChanges = (changes) => {
  if (!Array.isArray(changes) || changes.length === 0) {
    return freeze({
      impact: 'full',
      capabilities: freeze(['full']),
      valid: false,
      changedPathCount: 0,
      paths: freeze([]),
      reason: 'EMPTY_OR_INVALID_DIFF'
    });
  }

  const classifiedPaths = [];
  for (const change of changes) {
    const paths = pathsForChange(change);
    if (paths === null) {
      return freeze({
        impact: 'full',
        capabilities: freeze(['full']),
        valid: false,
        changedPathCount: classifiedPaths.length,
        paths: freeze(classifiedPaths),
        reason: 'UNKNOWN_OR_INVALID_GIT_STATUS'
      });
    }
    for (const path of paths) classifiedPaths.push(classifyPath(path));
  }

  const invalidPath = classifiedPaths.some(({ reason }) => reason === 'INVALID_REPOSITORY_PATH');
  const capabilities = invalidPath ? ['full'] : orderedCapabilities(classifiedPaths.map(({ impact }) => impact));
  return freeze({
    impact: primaryImpact(capabilities),
    capabilities: freeze(capabilities),
    valid: !invalidPath,
    changedPathCount: classifiedPaths.length,
    paths: freeze(classifiedPaths),
    reason: invalidPath ? 'INVALID_REPOSITORY_PATH' : 'CLASSIFIED'
  });
};

const validateLeafTests = ({ profile, capabilities, leafTests }) => {
  const invalid = (message, details = {}) => fail('INVALID_LEAF_TESTS', message, details);
  if (!Array.isArray(leafTests) || leafTests.length === 0) invalid('A leaf-test plan lists at least one leaf test.');
  leafTests.forEach((entry, index) => {
    if (!sameKeys(entry, LEAF_TEST_KEYS) || leafTestKind(entry.path) === null || entry.kind !== leafTestKind(entry.path)
        || !Array.isArray(entry.jobs)) {
      invalid('Leaf test entry has an invalid path, kind, or schema.', { index });
    }
    if (index > 0 && leafTests[index - 1].path >= entry.path) invalid('Leaf tests must be unique and sorted.', { index });
    const expected = entry.kind === 'python' && profile === 'dev'
      ? PYTEST_JOBS.filter((job) => entry.jobs.includes(job))
      : leafJobsFor({ profile, kind: entry.kind });
    if (JSON.stringify(entry.jobs) !== JSON.stringify(expected)
        || (entry.kind === 'python' && profile === 'dev' && entry.jobs[0] !== 'core')) {
      invalid('Leaf test jobs do not match its kind.', { path: entry.path, jobs: entry.jobs });
    }
  });
  const leafClasses = leafTests.map(({ path }) => classifyPath(path).impact);
  if (capabilities.some((capability) => !LIGHT_CLASSES.includes(capability) && !leafClasses.includes(capability))
      || leafClasses.some((capability) => !capabilities.includes(capability))
      || capabilities.some((capability) => ![...LIGHT_CLASSES, ...LEAF_SUBJECT_CLASSES].includes(capability))) {
    invalid('Leaf-test plan capabilities must come from documentation and the listed leaf tests.', { capabilities });
  }
};

export const requiredJobsFor = ({
  profile,
  impact,
  decision,
  capabilities = [impact],
  basis,
  leafTests = null
}) => {
  assertProfile(profile);
  assertImpact(impact);
  if (!IMPACT_DECISIONS.includes(decision)) {
    fail('UNKNOWN_DECISION', 'CI impact decision is not supported.', { decision });
  }
  if (!Array.isArray(capabilities) || !capabilities.length
      || capabilities.some((capability) => !IMPACT_CLASSES.includes(capability))
      || JSON.stringify(orderedCapabilities(capabilities)) !== JSON.stringify(capabilities)
      || primaryImpact(capabilities) !== impact) {
    fail('INVALID_CAPABILITIES', 'Capabilities must be known, ordered, unique, and match the impact.');
  }
  const all = PROFILE_REQUIRED_JOBS[profile];
  if (decision === 'selective' && basis === IDENTICAL_TREE_BASIS) return freeze([]);
  if (decision === 'selective' && basis === LEAF_TEST_BASIS) {
    validateLeafTests({ profile, capabilities, leafTests });
    const selected = [
      ...(profile === 'dev'
        ? capabilities.filter((capability) => LIGHT_CLASSES.includes(capability))
          .flatMap((capability) => PR_CAPABILITY_JOBS[capability])
        : []),
      ...leafTests.flatMap(({ jobs }) => jobs)
    ];
    return freeze(all.filter((job) => selected.includes(job)));
  }
  if (decision === 'selective' && requiresFullCoverage(profile, capabilities)) {
    fail('INVALID_SELECTIVE_PLAN', 'This impact requires full coverage.');
  }
  if (decision === 'full' && profile !== 'pr') return freeze([...all]);
  const selected = profile === 'pr'
    ? [
      ...(decision === 'full' ? PR_FULL_TIER_JOBS : []),
      ...capabilities.flatMap((capability) => PR_CAPABILITY_JOBS[capability])
    ]
    : profile === 'dev'
      ? capabilities.flatMap((capability) => PR_CAPABILITY_JOBS[capability] || [])
      : [];
  return freeze(all.filter((job) => selected.includes(job)));
};

export const knownJobsFor = (profile) => {
  assertProfile(profile);
  return freeze([...PROFILE_REQUIRED_JOBS[profile]]);
};

const validateInheritedEvidence = (evidence, expectedHeadSha) => {
  if (!sameKeys(evidence, EVIDENCE_KEYS)) {
    fail('INVALID_EVIDENCE_SCHEMA', 'Inherited evidence has an invalid schema.');
  }
  if (typeof evidence.workflowPath !== 'string' || !evidence.workflowPath
      || typeof evidence.aggregateName !== 'string' || !evidence.aggregateName
      || evidence.headSha !== expectedHeadSha
      || !Number.isSafeInteger(evidence.runId) || evidence.runId <= 0
      || !Number.isSafeInteger(evidence.aggregateJobId) || evidence.aggregateJobId <= 0
      || typeof evidence.runUrl !== 'string' || !evidence.runUrl
      || typeof evidence.aggregateJobUrl !== 'string' || !evidence.aggregateJobUrl) {
    fail('INVALID_EVIDENCE', 'Inherited evidence identity is invalid.');
  }
};

const validateBasis = (plan) => {
  if (!IMPACT_PLAN_BASES.includes(plan.basis)) {
    fail('UNKNOWN_BASIS', 'CI impact basis is not supported.', { basis: plan.basis });
  }
  const selectiveBases = plan.profile === 'pr'
    ? ['LIGHT_CHANGE_WITH_DIRECT_BASE_EVIDENCE']
    : plan.profile === 'release' ? [] : PARENT_EVIDENCE_BASES;
  const documentationOnlyPr = plan.basis === 'DOCUMENTATION_ONLY_PR';
  if (documentationOnlyPr && (plan.profile !== 'pr' || plan.decision !== 'selective'
      || !isDocumentationOnly(plan.capabilities))) {
    fail('BASIS_DECISION_MISMATCH', 'Documentation-only PR basis requires only documentation and metadata.');
  }
  if (plan.decision === 'selective' && !documentationOnlyPr && !selectiveBases.includes(plan.basis)) {
    fail('BASIS_DECISION_MISMATCH', 'Selective decision does not match its evidence basis.');
  }
  if (plan.decision === 'full' && [
    'LIGHT_CHANGE_WITH_DIRECT_BASE_EVIDENCE',
    ...PARENT_EVIDENCE_BASES
  ].includes(plan.basis)) {
    fail('BASIS_DECISION_MISMATCH', 'Full decision cannot claim inherited evidence.');
  }
  const identicalTree = Array.isArray(plan.capabilities) && plan.capabilities.includes('none');
  if ((plan.basis === IDENTICAL_TREE_BASIS || identicalTree)
      && (!identicalTree || plan.capabilities.length !== 1 || plan.changedPathCount !== 0
        || !['dev', 'gallery'].includes(plan.profile)
        || ![IDENTICAL_TREE_BASIS, 'INHERITED_EVIDENCE_UNAVAILABLE'].includes(plan.basis))) {
    fail('BASIS_IMPACT_MISMATCH', 'Only a proven identical tree has no changed path.');
  }
  if (plan.basis !== LEAF_TEST_BASIS && plan.leafTests !== null) {
    fail('UNEXPECTED_LEAF_TESTS', 'Only a leaf-test plan lists leaf tests.');
  }
  if (['UNKNOWN_OR_INVALID_CHANGE', 'MANUAL_FULL_RUN'].includes(plan.basis)
      && plan.impact !== 'full') {
    fail('BASIS_IMPACT_MISMATCH', 'Full or invalid change basis requires full impact.');
  }
  if (plan.basis === 'FULL_CHANGE' && !requiresFullCoverage(plan.profile, plan.capabilities)) {
    fail('BASIS_IMPACT_MISMATCH', 'Full change basis requires a full-coverage impact.');
  }
  if (plan.basis === 'INHERITED_EVIDENCE_UNAVAILABLE' && plan.impact === 'full') {
    fail('BASIS_IMPACT_MISMATCH', 'Evidence fallback requires a light impact candidate.');
  }
};

export const validateImpactPlan = (plan, expected = {}) => {
  if (!sameKeys(plan, PLAN_KEYS)) {
    fail('INVALID_PLAN_SCHEMA', 'Impact plan has an invalid schema.');
  }
  if (plan.schemaVersion !== IMPACT_PLAN_SCHEMA_VERSION) {
    fail('UNKNOWN_SCHEMA_VERSION', 'Impact plan schema version is not supported.', {
      schemaVersion: plan.schemaVersion
    });
  }
  assertProfile(plan.profile);
  assertImpact(plan.impact);
  if (!IMPACT_DECISIONS.includes(plan.decision)) {
    fail('UNKNOWN_DECISION', 'CI impact decision is not supported.', {
      decision: plan.decision
    });
  }
  validateBasis(plan);
  for (const field of ['changeBaseSha', 'changeHeadSha', 'workflowSha']) {
    if (!isFullObjectId(plan[field])) {
      fail('INVALID_OBJECT_ID', `${field} must be a full Git object ID.`, { field });
    }
  }
  if (!Number.isSafeInteger(plan.changedPathCount) || plan.changedPathCount < 0) {
    fail('INVALID_PATH_COUNT', 'changedPathCount must be a nonnegative integer.');
  }
  const expectedJobs = requiredJobsFor(plan);
  if (!Array.isArray(plan.requiredJobs)
      || plan.requiredJobs.length !== expectedJobs.length
      || plan.requiredJobs.some((job, index) => job !== expectedJobs[index])) {
    fail('REQUIRED_JOBS_MISMATCH', 'Impact plan required jobs do not match policy.', {
      expected: expectedJobs,
      observed: plan.requiredJobs
    });
  }
  if (plan.decision === 'selective' && plan.basis !== 'DOCUMENTATION_ONLY_PR') {
    validateInheritedEvidence(plan.inheritedEvidence, plan.changeBaseSha);
  } else if (plan.inheritedEvidence !== null) {
    fail('UNEXPECTED_EVIDENCE', 'Plans without an evidence basis cannot contain inherited evidence.');
  }
  if (expected.profile !== undefined && plan.profile !== expected.profile) {
    fail('PROFILE_MISMATCH', 'Impact plan profile does not match the gate.', {
      expected: expected.profile,
      observed: plan.profile
    });
  }
  if (expected.workflowSha !== undefined
      && plan.workflowSha.toLowerCase() !== String(expected.workflowSha).toLowerCase()) {
    fail('WORKFLOW_SHA_MISMATCH', 'Impact plan workflow SHA does not match the gate.', {
      expected: expected.workflowSha,
      observed: plan.workflowSha
    });
  }
  return true;
};

export const createImpactPlan = (fields) => {
  validateBasis({ ...fields, capabilities: fields.capabilities ?? [fields.impact], leafTests: fields.leafTests ?? null });
  const plan = {
    schemaVersion: IMPACT_PLAN_SCHEMA_VERSION,
    profile: fields.profile,
    impact: fields.impact,
    capabilities: freeze([...(fields.capabilities ?? [fields.impact])]),
    decision: fields.decision,
    basis: fields.basis,
    changeBaseSha: fields.changeBaseSha,
    changeHeadSha: fields.changeHeadSha,
    workflowSha: fields.workflowSha,
    requiredJobs: requiredJobsFor(fields),
    changedPathCount: fields.changedPathCount,
    inheritedEvidence: fields.inheritedEvidence,
    leafTests: fields.leafTests ?? null
  };
  validateImpactPlan(plan);
  if (plan.inheritedEvidence !== null) freeze(plan.inheritedEvidence);
  if (plan.leafTests !== null) {
    plan.leafTests = freeze(plan.leafTests.map((entry) => freeze({ ...entry, jobs: freeze([...entry.jobs]) })));
  }
  freeze(plan.requiredJobs);
  return freeze(plan);
};

const jobResult = (needs, jobId) => {
  const entry = needs[jobId];
  if (!isPlainObject(entry) || typeof entry.result !== 'string') {
    fail('MISSING_OR_INVALID_JOB', 'Expected CI job result is missing or invalid.', { jobId });
  }
  return entry.result;
};

export const validateGateResults = ({
  plan,
  needs,
  knownJobs = knownJobsFor(plan?.profile),
  expected = {}
}) => {
  validateImpactPlan(plan, expected);
  if (!isPlainObject(needs)) {
    fail('INVALID_NEEDS', 'Gate needs payload must be an object.');
  }
  if (!Array.isArray(knownJobs) || knownJobs.some((job) => typeof job !== 'string' || !job)
      || new Set(knownJobs).size !== knownJobs.length) {
    fail('INVALID_KNOWN_JOBS', 'Known CI job IDs must be a unique string array.');
  }

  if (jobResult(needs, 'ci-impact') !== 'success') {
    fail('PLANNER_NOT_SUCCESSFUL', 'CI impact planner did not succeed.');
  }

  const required = new Set(plan.requiredJobs);
  const observed = [];
  for (const jobId of knownJobs) {
    const result = jobResult(needs, jobId);
    if (required.has(jobId) && result !== 'success') {
      fail('REQUIRED_JOB_NOT_SUCCESSFUL', 'A required CI job did not succeed.', {
        jobId,
        result
      });
    }
    if (!required.has(jobId) && !['success', 'skipped'].includes(result)) {
      fail('UNREQUIRED_JOB_FAILED', 'An unrequired CI job failed unexpectedly.', {
        jobId,
        result
      });
    }
    observed.push(freeze({ jobId, required: required.has(jobId), result }));
  }

  for (const [jobId, entry] of Object.entries(needs)) {
    if (jobId === 'ci-impact' || knownJobs.includes(jobId)) continue;
    const result = isPlainObject(entry) ? entry.result : undefined;
    if (!['success', 'skipped'].includes(result)) {
      fail('UNKNOWN_JOB_FAILED', 'An unknown CI job failed unexpectedly.', { jobId, result });
    }
    observed.push(freeze({ jobId, required: false, result }));
  }

  return freeze({
    ok: true,
    profile: plan.profile,
    workflowSha: plan.workflowSha,
    inheritedEvidence: plan.inheritedEvidence !== null,
    jobs: freeze(observed)
  });
};

// SELECTIVE_CI.md "Carrying evidence to a later commit". `entries` hold each changed path with
// its leaf-test facts ({ kind, markers }) or null.
const DOC_PREFIXES = freeze({
  release: freeze(['docs/recipes/', 'docs/capture/']),
  generated: freeze(['docs/images/', 'docs/capture/', 'docs/recipes/']),
  local: freeze(['docs/recipes/'])
});
const DEV_TIER_LEAF_KINDS = freeze(['node', 'functional', 'performance', 'gallery-parity']);

export const carryForwardVerdicts = ({ ancestor, valid, identicalTree, entries }) => {
  const verdict = (accepts) => Boolean(ancestor) && Boolean(valid)
    && (Boolean(identicalTree) || (Array.isArray(entries) && entries.length > 0 && entries.every(accepts)));
  const documentary = ({ path }, prefixes) => LIGHT_CLASSES.includes(classifyPath(path).impact)
    && !prefixes.some((prefix) => path.startsWith(prefix));
  const markersKnown = (leaf) => Array.isArray(leaf.markers);
  return freeze({
    releaseEvidenceCarries: verdict((entry) => documentary(entry, DOC_PREFIXES.release)
      || (entry.leaf !== null && (DEV_TIER_LEAF_KINDS.includes(entry.leaf.kind)
        || (entry.leaf.kind === 'python' && markersKnown(entry.leaf) && entry.leaf.markers.length === 0)))),
    generatedArtifactChecksCarry: verdict((entry) => documentary(entry, DOC_PREFIXES.generated)
      || classifyPath(entry.path).impact === 'ci-only' || entry.leaf !== null),
    localTestEvidenceCarries: verdict((entry) => documentary(entry, DOC_PREFIXES.local)
      || (entry.leaf !== null && (entry.leaf.kind !== 'python'
        || (markersKnown(entry.leaf) && !entry.leaf.markers.includes('slow')))))
  });
};
