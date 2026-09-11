const freeze = (value) => Object.freeze(value);

export const IMPACT_PLAN_SCHEMA_VERSION = 2;

// Ordered only for the summary label. Routing unions every affected capability.
export const IMPACT_CLASSES = freeze([
  'metadata', 'documentation', 'tests-only', 'python-core', 'renderer',
  'web-runtime', 'session-persistence', 'gallery', 'losat-integration',
  'packaging', 'ci-only', 'full'
]);
export const IMPACT_DECISIONS = freeze(['selective', 'full']);
export const IMPACT_PLAN_BASES = freeze([
  'FULL_CHANGE', 'UNKNOWN_OR_INVALID_CHANGE', 'MANUAL_FULL_RUN',
  'ARCHITECTURE_CHANGE', 'LIGHT_CHANGE_WITH_DIRECT_BASE_EVIDENCE',
  'LIGHT_CHANGE_WITH_DIRECT_PARENT_EVIDENCE', 'INHERITED_EVIDENCE_UNAVAILABLE'
]);

const PR_JOBS = freeze([
  'web-change-budget', 'core-pr', 'recipes-standard', 'gallery', 'lint',
  'web-contracts-pr', 'web-pr-smoke'
]);
const DEV_JOBS = freeze([
  'web-change-budget', 'core', 'recipes-standard', 'gallery', 'browser',
  'playwright-functional', 'playwright-performance', 'lint',
  'losat-cache-browser-acceptance'
]);
const PROFILE_REQUIRED_JOBS = freeze({
  pr: PR_JOBS,
  dev: DEV_JOBS,
  release: freeze([...DEV_JOBS, 'acceptance-supported-main', 'slow-main']),
  gallery: freeze(['browser', 'performance'])
});
const PR_CAPABILITY_JOBS = freeze({
  metadata: freeze([]),
  documentation: freeze(['recipes-standard']),
  'tests-only': PR_JOBS,
  'python-core': freeze(['web-change-budget', 'core-pr', 'lint', 'web-contracts-pr', 'web-pr-smoke']),
  renderer: freeze(['web-change-budget', 'core-pr', 'lint', 'web-contracts-pr', 'web-pr-smoke']),
  'web-runtime': freeze(['web-change-budget', 'web-contracts-pr', 'web-pr-smoke']),
  'session-persistence': freeze(['web-change-budget', 'core-pr', 'recipes-standard', 'lint', 'web-contracts-pr', 'web-pr-smoke']),
  gallery: freeze(['web-change-budget', 'gallery', 'web-contracts-pr', 'web-pr-smoke']),
  'losat-integration': freeze(['web-change-budget', 'core-pr', 'lint', 'web-contracts-pr', 'web-pr-smoke']),
  packaging: PR_JOBS,
  'ci-only': PR_JOBS,
  full: PR_JOBS
});

// Runtime changes always receive comprehensive integrated-dev/Gallery validation.
// Control-plane, dependency, and unknown changes cannot inherit a narrower route.
export const requiresFullCoverage = (profile, capabilities) => profile === 'release'
  || capabilities.some((capability) => ['full', 'ci-only', 'packaging', 'tests-only'].includes(capability))
  || (profile !== 'pr' && capabilities.some((capability) => !['metadata', 'documentation'].includes(capability)));
const orderedCapabilities = (capabilities) => IMPACT_CLASSES.filter((capability) => capabilities.includes(capability));
const primaryImpact = (capabilities) => capabilities.at(-1);
const FULL_OBJECT_ID = /^(?:[0-9a-f]{40}|[0-9a-f]{64})$/i;
const ROOT_MARKDOWN = /^[^/]+\.md$/;
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
  'inheritedEvidence'
]);
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
  // Policy documents are executable authority even though they are Markdown.
  if (path.startsWith('.github/workflows/') || path.startsWith('tests/ci/')
      || /^tools\/(?:ci-impact|check-web|web-(?:architecture|product|change)|check-promotion)/.test(path)
      || /^tests\/web\/(?:architecture|product-impact|promotion-readiness).*\.test\.mjs$/.test(path)
      || /^playwright.*\.config\.js$/.test(path)
      || /^docs\/internal\/(?:WEB_CHANGE_POLICY|ARCHITECTURE_FITNESS_FUNCTION_RATCHET|PRODUCT_|OPTION_INTEGRITY_PRODUCT_CONTRACT|SELECTIVE_CI)/.test(path)) {
    return classified('ci-only', 'CI_OR_POLICY_AUTHORITY');
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

export const requiredJobsFor = ({ profile, impact, decision, capabilities = [impact] }) => {
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
  if (decision === 'selective' && requiresFullCoverage(profile, capabilities)) {
    fail('INVALID_SELECTIVE_PLAN', 'This impact requires full coverage.');
  }
  const all = PROFILE_REQUIRED_JOBS[profile];
  if (decision === 'full') return freeze([...all]);
  const selected = profile === 'pr'
    ? capabilities.flatMap((capability) => PR_CAPABILITY_JOBS[capability])
    : profile === 'dev' && capabilities.includes('documentation') ? ['recipes-standard'] : [];
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
  const selectiveBasis = plan.profile === 'pr'
    ? 'LIGHT_CHANGE_WITH_DIRECT_BASE_EVIDENCE'
    : 'LIGHT_CHANGE_WITH_DIRECT_PARENT_EVIDENCE';
  if (plan.decision === 'selective' && plan.basis !== selectiveBasis) {
    fail('BASIS_DECISION_MISMATCH', 'Selective decision does not match its evidence basis.');
  }
  if (plan.decision === 'full' && [
    'LIGHT_CHANGE_WITH_DIRECT_BASE_EVIDENCE',
    'LIGHT_CHANGE_WITH_DIRECT_PARENT_EVIDENCE'
  ].includes(plan.basis)) {
    fail('BASIS_DECISION_MISMATCH', 'Full decision cannot claim inherited evidence.');
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
  if (plan.decision === 'selective') {
    validateInheritedEvidence(plan.inheritedEvidence, plan.changeBaseSha);
  } else if (plan.inheritedEvidence !== null) {
    fail('UNEXPECTED_EVIDENCE', 'Full plans cannot contain inherited evidence.');
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
    inheritedEvidence: fields.inheritedEvidence
  };
  validateImpactPlan(plan);
  if (plan.inheritedEvidence !== null) freeze(plan.inheritedEvidence);
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
