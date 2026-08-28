const freeze = (value) => Object.freeze(value);

export const IMPACT_PLAN_SCHEMA_VERSION = 1;

export const IMPACT_CLASSES = freeze([
  'metadata',
  'documentation',
  'full'
]);

export const IMPACT_DECISIONS = freeze(['selective', 'full']);

export const IMPACT_PLAN_BASES = freeze([
  'FULL_CHANGE',
  'UNKNOWN_OR_INVALID_CHANGE',
  'MANUAL_FULL_RUN',
  'ARCHITECTURE_CHANGE',
  'LIGHT_CHANGE_WITH_DIRECT_BASE_EVIDENCE',
  'LIGHT_CHANGE_WITH_DIRECT_PARENT_EVIDENCE',
  'INHERITED_EVIDENCE_UNAVAILABLE'
]);

const PROFILE_REQUIRED_JOBS = freeze({
  pr: freeze({
    selective: freeze({
      metadata: freeze([]),
      documentation: freeze(['recipes-standard'])
    }),
    full: freeze([
      'web-change-budget',
      'core-pr',
      'recipes-standard',
      'gallery',
      'lint',
      'web-pr-smoke'
    ])
  }),
  dev: freeze({
    selective: freeze({
      metadata: freeze([]),
      documentation: freeze(['recipes-standard'])
    }),
    full: freeze([
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
    ])
  }),
  gallery: freeze({
    selective: freeze({
      metadata: freeze([]),
      documentation: freeze([])
    }),
    full: freeze(['browser', 'performance'])
  })
});

const IMPACT_WEIGHT = freeze({ metadata: 0, documentation: 1, full: 2 });
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
  if (!isValidRepositoryPath(path)) {
    return freeze({ path, impact: 'full', reason: 'INVALID_REPOSITORY_PATH' });
  }
  if (METADATA_FILES.has(path)) {
    return freeze({ path, impact: 'metadata', reason: 'METADATA_FILE_ALLOWLIST' });
  }
  if (METADATA_DIRECTORIES.some((directory) => path.startsWith(`${directory}/`))) {
    return freeze({ path, impact: 'metadata', reason: 'METADATA_DIRECTORY_ALLOWLIST' });
  }
  if (ROOT_MARKDOWN.test(path)) {
    return freeze({ path, impact: 'documentation', reason: 'ROOT_MARKDOWN' });
  }
  if (path.startsWith('docs/')) {
    return freeze({ path, impact: 'documentation', reason: 'DOCUMENTATION_TREE' });
  }
  return freeze({ path, impact: 'full', reason: 'FULL_BY_DEFAULT' });
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
        valid: false,
        changedPathCount: classifiedPaths.length,
        paths: freeze(classifiedPaths),
        reason: 'UNKNOWN_OR_INVALID_GIT_STATUS'
      });
    }
    for (const path of paths) classifiedPaths.push(classifyPath(path));
  }

  const invalidPath = classifiedPaths.some(({ reason }) => reason === 'INVALID_REPOSITORY_PATH');
  const impact = classifiedPaths.reduce(
    (current, entry) => IMPACT_WEIGHT[entry.impact] > IMPACT_WEIGHT[current]
      ? entry.impact
      : current,
    'metadata'
  );
  return freeze({
    impact: invalidPath ? 'full' : impact,
    valid: !invalidPath,
    changedPathCount: classifiedPaths.length,
    paths: freeze(classifiedPaths),
    reason: invalidPath ? 'INVALID_REPOSITORY_PATH' : 'CLASSIFIED'
  });
};

export const requiredJobsFor = ({ profile, impact, decision }) => {
  assertProfile(profile);
  assertImpact(impact);
  if (!IMPACT_DECISIONS.includes(decision)) {
    fail('UNKNOWN_DECISION', 'CI impact decision is not supported.', { decision });
  }
  if (decision === 'selective' && impact === 'full') {
    fail('INVALID_SELECTIVE_PLAN', 'A full impact cannot use selective execution.');
  }
  const jobs = decision === 'full'
    ? PROFILE_REQUIRED_JOBS[profile].full
    : PROFILE_REQUIRED_JOBS[profile].selective[impact];
  return freeze([...jobs]);
};

export const knownJobsFor = (profile) => {
  assertProfile(profile);
  return freeze([...PROFILE_REQUIRED_JOBS[profile].full]);
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
  if (['FULL_CHANGE', 'UNKNOWN_OR_INVALID_CHANGE', 'MANUAL_FULL_RUN'].includes(plan.basis)
      && plan.impact !== 'full') {
    fail('BASIS_IMPACT_MISMATCH', 'Full or invalid change basis requires full impact.');
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
