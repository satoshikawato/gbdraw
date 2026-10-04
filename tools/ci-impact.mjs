#!/usr/bin/env node

import { spawnSync } from 'node:child_process';
import { appendFileSync, readFileSync } from 'node:fs';
import { resolve } from 'node:path';
import { TextDecoder } from 'node:util';
import { pathToFileURL } from 'node:url';
import {
  PromotionReadinessError,
  verifyWorkflowEvidence
} from './check-promotion-readiness.mjs';
import {
  IMPACT_PLAN_SCHEMA_VERSION,
  carryForwardVerdicts,
  classifyChanges,
  classifyPath,
  createImpactPlan,
  isFullObjectId,
  isDocumentationOnly,
  knownJobsFor,
  leafJobsFor,
  leafTestKind,
  requiresFullCoverage,
  validateGateResults,
  validateImpactPlan
} from './ci-impact-policy.mjs';

const SUPPORTED_PROFILES = new Set(['pr', 'dev', 'gallery', 'release']);
const SUPPORTED_EVENTS = new Set(['pull_request', 'push', 'workflow_dispatch']);
const REPOSITORY_NAME = /^[A-Za-z0-9](?:[A-Za-z0-9._-]{0,99})\/[A-Za-z0-9](?:[A-Za-z0-9._-]{0,99})$/;
const SUMMARY_PATH_LIMIT = 20;
const utf8Decoder = new TextDecoder('utf-8', { fatal: true });

const EVIDENCE_CONTRACTS = Object.freeze({
  pr: Object.freeze({
    workflowPath: '.github/workflows/test.yml',
    aggregateName: 'Dev staging / gate'
  }),
  dev: Object.freeze({
    workflowPath: '.github/workflows/test.yml',
    aggregateName: 'Dev staging / gate'
  }),
  gallery: Object.freeze({
    workflowPath: '.github/workflows/gallery-publication.yml',
    aggregateName: 'Gallery readiness / gate'
  })
});

const isPlainObject = (value) => value !== null
  && typeof value === 'object'
  && !Array.isArray(value);

const fail = (code, message, details = {}) => {
  throw Object.assign(new Error(message), {
    name: 'CiImpactCliError',
    code,
    details
  });
};

const boundedText = (value, limit = 300) => {
  const source = typeof value === 'string' ? value : String(value ?? '');
  const normalized = source.replace(/[\r\n\t]+/g, ' ').trim();
  return normalized.length <= limit ? normalized : `${normalized.slice(0, limit)}...`;
};

const redact = (text, env) => {
  const token = typeof env.GITHUB_TOKEN === 'string' ? env.GITHUB_TOKEN : '';
  return token ? text.replaceAll(token, '[REDACTED]') : text;
};

const requiredEnvironment = (env, name) => {
  const value = env[name];
  if (typeof value !== 'string' || !value) {
    fail('MISSING_ENVIRONMENT', `Required environment variable ${name} is missing.`, { name });
  }
  return value;
};

const booleanEnvironment = (env, name) => {
  const value = requiredEnvironment(env, name);
  if (!['true', 'false'].includes(value)) {
    fail('INVALID_BOOLEAN', `${name} must be true or false.`, { name, value: boundedText(value) });
  }
  return value === 'true';
};

const objectIdEnvironment = (env, name) => {
  const value = requiredEnvironment(env, name).toLowerCase();
  if (!isFullObjectId(value)) {
    fail('INVALID_OBJECT_ID', `${name} must be a full Git object ID.`, { name });
  }
  return value;
};

export const readPlanConfiguration = (env, cwd = process.cwd()) => {
  const profile = requiredEnvironment(env, 'CI_IMPACT_PROFILE');
  const eventName = requiredEnvironment(env, 'CI_IMPACT_EVENT_NAME');
  const repository = requiredEnvironment(env, 'CI_IMPACT_REPOSITORY');
  if (!SUPPORTED_PROFILES.has(profile)) {
    fail('INVALID_PROFILE', 'CI_IMPACT_PROFILE is not supported.', { profile });
  }
  if (!SUPPORTED_EVENTS.has(eventName)) {
    fail('INVALID_EVENT', 'CI_IMPACT_EVENT_NAME is not supported.', { eventName });
  }
  if ((profile === 'pr') !== (eventName === 'pull_request')) {
    fail('PROFILE_EVENT_MISMATCH', 'PR profile must match a pull_request event.');
  }
  if (profile === 'release' && eventName !== 'workflow_dispatch') {
    fail('PROFILE_EVENT_MISMATCH', 'Release acceptance requires an explicit dispatch.');
  }
  if (!REPOSITORY_NAME.test(repository)) {
    fail('INVALID_REPOSITORY', 'CI_IMPACT_REPOSITORY must be OWNER/REPOSITORY.');
  }
  return Object.freeze({
    profile,
    eventName,
    repository,
    repositoryRoot: resolve(env.CI_IMPACT_REPOSITORY_ROOT || cwd),
    changeBaseSha: objectIdEnvironment(env, 'CI_IMPACT_CHANGE_BASE_SHA'),
    changeHeadSha: objectIdEnvironment(env, 'CI_IMPACT_CHANGE_HEAD_SHA'),
    workflowSha: objectIdEnvironment(env, 'CI_IMPACT_WORKFLOW_SHA'),
    architectureChange: booleanEnvironment(env, 'CI_IMPACT_ARCHITECTURE_CHANGE')
  });
};

const splitNulTokens = (source) => {
  const bytes = Buffer.isBuffer(source) ? source : Buffer.from(source);
  if (bytes.length === 0 || bytes[bytes.length - 1] !== 0) {
    fail('MALFORMED_GIT_DIFF', 'Git name-status output is empty or not NUL-terminated.');
  }
  const tokens = [];
  let start = 0;
  for (let index = 0; index < bytes.length; index += 1) {
    if (bytes[index] !== 0) continue;
    const token = bytes.subarray(start, index);
    try {
      tokens.push(utf8Decoder.decode(token));
    } catch (_error) {
      fail('INVALID_PATH_ENCODING', 'Git path is not valid UTF-8.');
    }
    start = index + 1;
  }
  if (tokens.at(-1) === '') tokens.pop();
  return tokens;
};

export const parseNameStatusZ = (source) => {
  const tokens = splitNulTokens(source);
  if (tokens.length === 0) {
    fail('MALFORMED_GIT_DIFF', 'Git name-status output contains no changes.');
  }
  const changes = [];
  for (let index = 0; index < tokens.length;) {
    const status = tokens[index];
    index += 1;
    const pathCount = /^[RC]/.test(status) ? 2 : 1;
    if (index + pathCount > tokens.length) {
      fail('MALFORMED_GIT_DIFF', 'Git name-status output ended before its paths.');
    }
    changes.push(Object.freeze({
      status,
      paths: Object.freeze(tokens.slice(index, index + pathCount))
    }));
    index += pathCount;
  }
  return Object.freeze(changes);
};

export const runGit = (repositoryRoot, args) => spawnSync('git', args, {
  cwd: repositoryRoot,
  encoding: null,
  stdio: ['ignore', 'pipe', 'pipe']
});

const gitText = (runGitImpl, repositoryRoot, args) => {
  try {
    const result = runGitImpl(repositoryRoot, args);
    if (!isPlainObject(result) || !Number.isInteger(result.status)
        || !(Buffer.isBuffer(result.stdout) || result.stdout instanceof Uint8Array)) {
      return null;
    }
    return { status: result.status, text: Buffer.from(result.stdout).toString('utf8') };
  } catch (_error) {
    return null;
  }
};

// An empty diff counts only when both tree object IDs are proven equal.
const sameTrees = ({ runGitImpl, repositoryRoot, baseSha, headSha }) => {
  const result = gitText(runGitImpl, repositoryRoot, ['rev-parse', `${baseSha}^{tree}`, `${headSha}^{tree}`]);
  if (result === null || result.status !== 0) return false;
  const trees = result.text.trim().split('\n');
  return trees.length === 2 && isFullObjectId(trees[0]) && trees[0] === trees[1];
};

const readGitChanges = ({ configuration, runGitImpl }) => {
  const range = configuration.profile === 'pr'
    ? [`${configuration.changeBaseSha}...${configuration.changeHeadSha}`]
    : [configuration.changeBaseSha, configuration.changeHeadSha];
  let result;
  try {
    result = runGitImpl(configuration.repositoryRoot, [
      'diff',
      '--name-status',
      '-z',
      '--find-renames',
      ...range,
      '--'
    ]);
  } catch (error) {
    return Object.freeze({
      valid: false,
      reason: 'GIT_DIFF_FAILED',
      diagnostic: boundedText(error?.message || error)
    });
  }
  if (!isPlainObject(result) || !Number.isInteger(result.status)
      || !(Buffer.isBuffer(result.stdout) || result.stdout instanceof Uint8Array)) {
    return Object.freeze({
      valid: false,
      reason: 'MALFORMED_GIT_RESULT',
      diagnostic: 'Git runner returned an invalid result.'
    });
  }
  if (result.status !== 0) {
    return Object.freeze({
      valid: false,
      reason: 'GIT_DIFF_FAILED',
      diagnostic: boundedText(Buffer.from(result.stderr || '').toString('utf8'))
    });
  }
  if (result.stdout.length === 0 && configuration.profile !== 'pr' && sameTrees({
    runGitImpl,
    repositoryRoot: configuration.repositoryRoot,
    baseSha: configuration.changeBaseSha,
    headSha: configuration.changeHeadSha
  })) {
    return Object.freeze({ valid: true, identicalTree: true, changes: Object.freeze([]) });
  }
  try {
    return Object.freeze({ valid: true, identicalTree: false, changes: parseNameStatusZ(result.stdout) });
  } catch (error) {
    return Object.freeze({
      valid: false,
      reason: error?.code || 'MALFORMED_GIT_DIFF',
      diagnostic: boundedText(error?.message || error)
    });
  }
};

const LIGHT_CLASSES = new Set(['metadata', 'documentation', 'policy-documentation']);
const isLightPath = (path) => LIGHT_CLASSES.has(classifyPath(path).impact);
const IDENTICAL_TREE_CLASSIFICATION = Object.freeze({
  impact: 'none',
  capabilities: Object.freeze(['none']),
  valid: true,
  changedPathCount: 0,
  paths: Object.freeze([]),
  reason: 'IDENTICAL_TREE'
});
const FUNCTIONAL_SHARDS_PATH = 'tests/ci/functional-shards.json';
// The reference check ignores documentation, metadata, and the runner lists.
const REFERENCE_SCAN_EXCLUSIONS = Object.freeze([
  'docs/', '.agents/', '.claude/', '.codex/', '.cursor/', 'tests/ci/', '.github/workflows/'
].map((path) => `:(exclude)${path}`));
const isRunnerList = (path) => path.startsWith('tests/ci/') || path.startsWith('.github/workflows/')
  || /^playwright[^/]*\.config\.js$/.test(path) || path === 'package.json';
const PYTEST_SELECTION_MARKER = /\bmark\.(recipe|gallery|browser|slow)\b/g;
const UNREADABLE_PYTEST_MARKER = /getattr\(\s*(?:pytest\.)?mark\b|\bmark\s*\[|add_marker|MarkDecorator/;
const escapeRegExp = (value) => value.replace(/[.*+?^${}()|[\]\\]/g, '\\$&');

const changeEntries = (changes) => changes.flatMap(({ status, paths }) => {
  if (status === 'D') return [{ path: paths[0], status, present: false }];
  if (status.startsWith('R')) {
    return [{ path: paths[0], status, present: false }, { path: paths[1], status, present: true }];
  }
  return paths.map((path) => ({ path, status, present: true }));
});

const readBlob = ({ runGitImpl, repositoryRoot, sha, path }) => {
  const result = gitText(runGitImpl, repositoryRoot, ['cat-file', 'blob', `${sha}:${path}`]);
  return result !== null && result.status === 0 ? result.text : null;
};

const pythonMarkers = (source) => {
  if (source === null || UNREADABLE_PYTEST_MARKER.test(source)) return null;
  return Object.freeze([...new Set([...source.matchAll(PYTEST_SELECTION_MARKER)].map((match) => match[1]))].sort());
};

const referencePatterns = (path) => {
  const name = path.split('/').pop();
  const patterns = [new RegExp(`(?:^|[^A-Za-z0-9_.-])${escapeRegExp(name)}(?![A-Za-z0-9_-])`)];
  if (path.endsWith('.py')) {
    const module = escapeRegExp(name.slice(0, -3));
    patterns.push(
      new RegExp(`(?:^|[^A-Za-z0-9_])tests\\.${module}(?![A-Za-z0-9_])`),
      new RegExp(`\\bfrom\\s+\\.?${module}\\s+import\\b`),
      new RegExp(`\\bimport\\s+${module}(?![A-Za-z0-9_])`),
      new RegExp(`\\bfrom\\s+tests\\s+import\\b[^\\n]*\\b${module}(?![A-Za-z0-9_])`)
    );
  }
  return patterns;
};

// Returns the candidate paths that another tracked file names, or null when the scan fails.
const referencedCandidates = ({ runGitImpl, repositoryRoot, headSha, paths }) => {
  const needles = [...new Set(paths.flatMap((path) => {
    const name = path.split('/').pop();
    return path.endsWith('.py') ? [name, name.slice(0, -3)] : [name];
  }))];
  const result = gitText(runGitImpl, repositoryRoot, [
    'grep', '-I', '-z', '-F', ...needles.flatMap((needle) => ['-e', needle]),
    headSha, '--', '.', ...REFERENCE_SCAN_EXCLUSIONS
  ]);
  if (result === null || ![0, 1].includes(result.status)) return null;
  const prefix = `${headSha}:`;
  const records = result.status === 0 ? result.text.split('\n').filter(Boolean) : [];
  const patterns = new Map(paths.map((path) => [path, referencePatterns(path)]));
  const referenced = new Set();
  for (const record of records) {
    const separator = record.indexOf('\0');
    if (separator < 0) continue;
    const name = record.slice(0, separator);
    const file = name.startsWith(prefix) ? name.slice(prefix.length) : name;
    const line = record.slice(separator + 1);
    if (isLightPath(file) || isRunnerList(file)) continue;
    for (const [path, expressions] of patterns) {
      if (file !== path && expressions.some((expression) => expression.test(line))) referenced.add(path);
    }
  }
  return referenced;
};

// Leaf-test facts for every changed path whose pattern can be a leaf test (SELECTIVE_CI.md
// "Leaf tests"): its kind, its pytest markers, and whether it is a leaf test.
export const inspectLeafTests = ({ runGitImpl, repositoryRoot, headSha, changes }) => {
  const candidates = new Map();
  for (const entry of changeEntries(changes)) {
    if (leafTestKind(entry.path) === null) continue;
    candidates.set(entry.path, Boolean(candidates.get(entry.path)) || entry.present);
  }
  const facts = new Map();
  if (candidates.size === 0) return facts;
  const referenced = referencedCandidates({ runGitImpl, repositoryRoot, headSha, paths: [...candidates.keys()] });
  let shardSpecs = null;
  if ([...candidates.keys()].some((path) => leafTestKind(path) === 'functional')) {
    try {
      const shards = JSON.parse(readBlob({ runGitImpl, repositoryRoot, sha: headSha, path: FUNCTIONAL_SHARDS_PATH }));
      shardSpecs = new Set(shards.shards.flat());
    } catch (_error) {
      shardSpecs = null;
    }
  }
  for (const [path, present] of candidates) {
    const kind = leafTestKind(path);
    const markers = kind === 'python' && present
      ? pythonMarkers(readBlob({ runGitImpl, repositoryRoot, sha: headSha, path }))
      : null;
    const leaf = referenced !== null && !referenced.has(path)
      && (kind !== 'functional' || Boolean(shardSpecs?.has(path)));
    facts.set(path, Object.freeze({ kind, markers, leaf }));
  }
  return facts;
};

// The leaf tests of a dev or Gallery push that changes only documentation, metadata, policy
// documentation, and leaf tests; null when any changed path keeps the complete tier.
const leafTestRoute = ({ configuration, changes, runGitImpl }) => {
  const entries = changeEntries(changes);
  if (!entries.some(({ path }) => leafTestKind(path) !== null)
      || !entries.every(({ path }) => isLightPath(path) || leafTestKind(path) !== null)) {
    return null;
  }
  const facts = inspectLeafTests({
    runGitImpl,
    repositoryRoot: configuration.repositoryRoot,
    headSha: configuration.changeHeadSha,
    changes
  });
  if (entries.some(({ path }) => !isLightPath(path) && !facts.get(path)?.leaf)) return null;
  return Object.freeze([...facts.entries()]
    .sort(([left], [right]) => (left < right ? -1 : 1))
    .map(([path, { kind, markers }]) => Object.freeze({
      path,
      kind,
      jobs: leafJobsFor({ profile: configuration.profile, kind, markers })
    })));
};

const planFields = ({ configuration, classification, decision, basis, inheritedEvidence, leafTests = null }) => ({
  profile: configuration.profile,
  impact: classification.impact,
  capabilities: classification.capabilities,
  decision,
  basis,
  changeBaseSha: configuration.changeBaseSha,
  changeHeadSha: configuration.changeHeadSha,
  workflowSha: configuration.workflowSha,
  changedPathCount: classification.changedPathCount,
  inheritedEvidence,
  leafTests
});

const fullClassification = (reason) => Object.freeze({
  impact: 'full',
  capabilities: Object.freeze(['full']),
  valid: false,
  changedPathCount: 0,
  paths: Object.freeze([]),
  reason
});

const normalizeEvidence = (evidence, contract, expectedHeadSha) => {
  if (!isPlainObject(evidence) || !isPlainObject(evidence.run)
      || !isPlainObject(evidence.aggregateJob)
      || evidence.workflow?.path !== contract.workflowPath
      || evidence.run.headSha !== expectedHeadSha
      || !Number.isSafeInteger(evidence.run.id) || evidence.run.id <= 0
      || !Number.isSafeInteger(evidence.aggregateJob.id) || evidence.aggregateJob.id <= 0
      || evidence.aggregateJob.name !== contract.aggregateName
      || typeof evidence.run.url !== 'string' || !evidence.run.url
      || typeof evidence.aggregateJob.url !== 'string' || !evidence.aggregateJob.url) {
    fail('MALFORMED_EVIDENCE', 'Evidence verifier returned an invalid success result.');
  }
  return {
    workflowPath: contract.workflowPath,
    aggregateName: contract.aggregateName,
    headSha: expectedHeadSha,
    runId: evidence.run.id,
    aggregateJobId: evidence.aggregateJob.id,
    runUrl: evidence.run.url,
    aggregateJobUrl: evidence.aggregateJob.url
  };
};

export const buildImpactPlan = async ({
  configuration,
  token,
  runGitImpl = runGit,
  verifyWorkflowEvidenceImpl = verifyWorkflowEvidence
}) => {
  if (configuration.eventName === 'workflow_dispatch') {
    const classification = fullClassification('MANUAL_FULL_RUN');
    return Object.freeze({
      plan: createImpactPlan(planFields({
        configuration,
        classification,
        decision: 'full',
        basis: 'MANUAL_FULL_RUN',
        inheritedEvidence: null
      })),
      classification
    });
  }

  const gitChanges = readGitChanges({ configuration, runGitImpl });
  const classification = !gitChanges.valid
    ? fullClassification(gitChanges.reason)
    : gitChanges.identicalTree ? IDENTICAL_TREE_CLASSIFICATION : classifyChanges(gitChanges.changes);

  if (configuration.architectureChange) {
    return Object.freeze({
      plan: createImpactPlan(planFields({
        configuration,
        classification,
        decision: 'full',
        basis: 'ARCHITECTURE_CHANGE',
        inheritedEvidence: null
      })),
      classification
    });
  }

  const parentEvidenceProfile = ['dev', 'gallery'].includes(configuration.profile);
  const identicalTree = parentEvidenceProfile && classification.reason === 'IDENTICAL_TREE';
  const leafTests = classification.valid && parentEvidenceProfile && !identicalTree
      && requiresFullCoverage(configuration.profile, classification.capabilities)
    ? leafTestRoute({ configuration, changes: gitChanges.changes, runGitImpl })
    : null;
  if (!classification.valid || (!identicalTree && leafTests === null
      && requiresFullCoverage(configuration.profile, classification.capabilities))) {
    return Object.freeze({
      plan: createImpactPlan(planFields({
        configuration,
        classification,
        decision: 'full',
        basis: classification.valid ? 'FULL_CHANGE' : 'UNKNOWN_OR_INVALID_CHANGE',
        inheritedEvidence: null
      })),
      classification
    });
  }

  if (configuration.profile === 'pr' && isDocumentationOnly(classification.capabilities)) {
    return Object.freeze({
      plan: createImpactPlan(planFields({
        configuration,
        classification,
        decision: 'selective',
        basis: 'DOCUMENTATION_ONLY_PR',
        inheritedEvidence: null
      })),
      classification
    });
  }

  const contract = EVIDENCE_CONTRACTS[configuration.profile];
  let evidence;
  try {
    evidence = await verifyWorkflowEvidenceImpl({
      repository: configuration.repository,
      workflowPath: contract.workflowPath,
      expectedHeadSha: configuration.changeBaseSha,
      expectedAggregateName: contract.aggregateName,
      token
    });
  } catch (error) {
    if (!(error instanceof PromotionReadinessError)) throw error;
    // A dev push falls back to the full dev tier; Gallery publication keeps failing closed.
    if (configuration.profile === 'gallery'
        && isDocumentationOnly(classification.capabilities)
        && !classification.capabilities.includes('metadata')) {
      fail('DOCUMENTATION_BASE_EVIDENCE_UNAVAILABLE',
        'Documentation-only changes cannot inherit the required baseline CI evidence.',
        { evidenceCode: error.code });
    }
    return Object.freeze({
      plan: createImpactPlan(planFields({
        configuration,
        classification,
        decision: 'full',
        basis: 'INHERITED_EVIDENCE_UNAVAILABLE',
        inheritedEvidence: null
      })),
      classification,
      evidenceFailure: Object.freeze({
        code: error.code,
        reason: token
          ? boundedText(error.message).replaceAll(token, '[REDACTED]')
          : boundedText(error.message)
      })
    });
  }

  const inheritedEvidence = normalizeEvidence(
    evidence,
    contract,
    configuration.changeBaseSha
  );
  const basis = configuration.profile === 'pr'
    ? 'LIGHT_CHANGE_WITH_DIRECT_BASE_EVIDENCE'
    : identicalTree ? 'IDENTICAL_TREE_WITH_DIRECT_PARENT_EVIDENCE'
      : leafTests !== null ? 'LEAF_TEST_CHANGE_WITH_DIRECT_PARENT_EVIDENCE'
        : 'LIGHT_CHANGE_WITH_DIRECT_PARENT_EVIDENCE';
  return Object.freeze({
    plan: createImpactPlan(planFields({
      configuration,
      classification,
      decision: 'selective',
      basis,
      inheritedEvidence,
      leafTests
    })),
    classification
  });
};

const validShards = (shards) => Array.isArray(shards) && shards.length > 0
  && shards.every((paths) => Array.isArray(paths) && paths.every((path) => typeof path === 'string' && /^\S+$/.test(path)));

// SELECTIVE_CI.md "Jobs for changed leaf tests": a leaf-test plan runs only the changed specs
// assigned to the shard; every other plan runs the whole shard. A pull request runs the
// trusted base's plan, which may use an older schema; such a plan runs the whole shard.
export const functionalShardFiles = ({ plan, shards, shard }) => {
  const currentSchema = isPlainObject(plan) && plan.schemaVersion === IMPACT_PLAN_SCHEMA_VERSION;
  if (currentSchema) validateImpactPlan(plan);
  if (!isPlainObject(plan) || !Array.isArray(plan.requiredJobs) || !plan.requiredJobs.includes('playwright-functional')) {
    fail('FUNCTIONAL_NOT_REQUIRED', 'The plan does not require functional Playwright.');
  }
  if (!validShards(shards)) fail('INVALID_SHARD_MAP', 'The functional shard map is invalid.');
  if (!Number.isInteger(shard) || shard < 1 || shard > shards.length) {
    fail('INVALID_SHARD', 'The shard number is outside the shard map.', { shard });
  }
  const assigned = shards[shard - 1];
  if (!currentSchema || plan.leafTests === null) return Object.freeze([...assigned]);
  const changed = new Set(plan.leafTests.filter(({ kind }) => kind === 'functional').map(({ path }) => path));
  return Object.freeze(assigned.filter((path) => changed.has(path)));
};

const resolveCommit = ({ runGitImpl, repositoryRoot, revision }) => {
  const result = gitText(runGitImpl, repositoryRoot, ['rev-parse', '--verify', '--end-of-options', `${revision}^{commit}`]);
  const sha = result?.status === 0 ? result.text.trim().toLowerCase() : '';
  if (!isFullObjectId(sha)) fail('UNKNOWN_REVISION', 'Revision does not name a commit.', { revision: boundedText(revision) });
  return sha;
};

// `classify --base E --head H`: local Git data only; no workflow evidence.
export const classifyRange = ({ repositoryRoot, base, head, runGitImpl = runGit }) => {
  const baseSha = resolveCommit({ runGitImpl, repositoryRoot, revision: base });
  const headSha = resolveCommit({ runGitImpl, repositoryRoot, revision: head });
  const ancestry = gitText(runGitImpl, repositoryRoot, ['merge-base', '--is-ancestor', baseSha, headSha]);
  if (ancestry === null || ![0, 1].includes(ancestry.status)) {
    fail('ANCESTRY_CHECK_FAILED', 'Git could not decide whether the base is an ancestor of the head.');
  }
  const gitChanges = readGitChanges({
    configuration: { profile: 'dev', repositoryRoot, changeBaseSha: baseSha, changeHeadSha: headSha },
    runGitImpl
  });
  const classification = !gitChanges.valid
    ? fullClassification(gitChanges.reason)
    : gitChanges.identicalTree ? IDENTICAL_TREE_CLASSIFICATION : classifyChanges(gitChanges.changes);
  const facts = gitChanges.valid
    ? inspectLeafTests({ runGitImpl, repositoryRoot, headSha, changes: gitChanges.changes })
    : new Map();
  const statuses = new Map(gitChanges.valid
    ? changeEntries(gitChanges.changes).map(({ path, status }) => [path, status])
    : []);
  const paths = classification.paths.map(({ path, impact, reason }) => {
    const fact = facts.get(path);
    return Object.freeze({
      path,
      status: statuses.get(path),
      impact,
      reason,
      leaf: fact?.leaf
        ? Object.freeze({
          kind: fact.kind,
          markers: fact.markers,
          jobs: leafJobsFor({ profile: 'dev', kind: fact.kind, markers: fact.markers })
        })
        : null
    });
  });
  return Object.freeze({
    base: baseSha,
    head: headSha,
    ancestor: ancestry.status === 0,
    identicalTree: classification.reason === 'IDENTICAL_TREE',
    valid: classification.valid,
    reason: classification.reason,
    capabilities: classification.capabilities,
    paths: Object.freeze(paths),
    verdicts: carryForwardVerdicts({
      ancestor: ancestry.status === 0,
      valid: classification.valid,
      identicalTree: classification.reason === 'IDENTICAL_TREE',
      entries: paths
    })
  });
};

const visiblePath = (path) => String(path)
  .replaceAll('\\', '\\\\')
  .replaceAll('\r', '\\r')
  .replaceAll('\n', '\\n')
  .replaceAll('\t', '\\t');

const html = (source) => visiblePath(source)
  .replaceAll('&', '&amp;')
  .replaceAll('<', '&lt;')
  .replaceAll('>', '&gt;');

export const formatPlanSummary = ({ plan, classification, evidenceFailure }) => {
  const requiredJobs = plan.requiredJobs.length ? plan.requiredJobs.map((job) => `\`${job}\``).join(', ') : 'none';
  const routing = plan.profile === 'pr'
    ? 'active; pull-request jobs use the trusted-base plan'
    : ['dev', 'release'].includes(plan.profile)
      ? 'active; dev staging jobs use the protected-branch plan'
      : 'active; Gallery readiness jobs use the protected-branch plan';
  const lines = [
    '## CI impact plan',
    '',
    `- Profile: \`${plan.profile}\``,
    `- Impact: \`${plan.impact}\``,
    `- Capabilities: ${plan.capabilities.join(', ')}`,
    `- Decision: \`${plan.decision}\``,
    `- Basis: \`${plan.basis}\``,
    `- Change base/head: \`${plan.changeBaseSha}\` / \`${plan.changeHeadSha}\``,
    `- Workflow SHA: \`${plan.workflowSha}\``,
    `- Changed paths: ${plan.changedPathCount}`,
    `- Planned required job IDs: ${requiredJobs}`,
    `- Routing: ${routing}`,
    ''
  ];
  if (plan.inheritedEvidence) {
    lines.push(
      `- Inherited run: [${plan.inheritedEvidence.runId}](${plan.inheritedEvidence.runUrl})`,
      `- Inherited aggregate: [${html(plan.inheritedEvidence.aggregateName)}](${plan.inheritedEvidence.aggregateJobUrl})`,
      `- Evidence SHA: \`${plan.inheritedEvidence.headSha}\``,
      ''
    );
  }
  if (plan.leafTests) {
    lines.push('### Leaf tests', '');
    plan.leafTests.forEach(({ path, kind, jobs }) => {
      lines.push(`- <code>${html(path)}</code> — \`${kind}\` → ${jobs.length ? jobs.map((job) => `\`${job}\``).join(', ') : 'no job'}`);
    });
    lines.push('');
  }
  if (evidenceFailure) {
    lines.push(
      `- Full fallback: \`${html(evidenceFailure.code)}\` — ${html(evidenceFailure.reason)}`,
      ''
    );
  }
  if (classification.paths.length) {
    lines.push('### Classified paths', '');
    classification.paths.slice(0, SUMMARY_PATH_LIMIT).forEach((entry) => {
      lines.push(`- <code>${html(entry.path)}</code> — \`${entry.impact}\` / \`${entry.reason}\``);
    });
    if (classification.paths.length > SUMMARY_PATH_LIMIT) {
      lines.push(`- … ${classification.paths.length - SUMMARY_PATH_LIMIT} more path(s)`);
    }
    lines.push('');
  } else if (!classification.valid) {
    lines.push(`- Diff fallback reason: \`${html(classification.reason)}\``, '');
  }
  return lines.join('\n');
};

export const formatGateSummary = ({ plan, result }) => [
  '## CI impact aggregate validation',
  '',
  `- Profile: \`${plan.profile}\``,
  `- Workflow SHA: \`${plan.workflowSha}\``,
  `- Inherited evidence: ${result.inheritedEvidence ? 'yes' : 'no'}`,
  ...result.jobs.map(({ jobId, required, result: jobStatus }) => (
    `- \`${jobId}\`: \`${jobStatus}\` (${required ? 'required' : 'not required'})`
  )),
  '- Gate result: pass',
  ''
].join('\n');

const appendOutput = (path, content, appendFileImpl) => {
  try {
    appendFileImpl(resolve(path), content, 'utf8');
  } catch (_error) {
    fail('OUTPUT_WRITE_FAILED', 'CI output file could not be written.');
  }
};

const parseJsonEnvironment = (env, name) => {
  const source = requiredEnvironment(env, name);
  try {
    return JSON.parse(source);
  } catch (_error) {
    fail('MALFORMED_JSON', `${name} must contain valid JSON.`, { name });
  }
};

const runPlanCommand = async ({
  env,
  cwd,
  stdout,
  appendFileImpl,
  runGitImpl,
  verifyWorkflowEvidenceImpl
}) => {
  const configuration = readPlanConfiguration(env, cwd);
  const outcome = await buildImpactPlan({
    configuration,
    token: env.GITHUB_TOKEN,
    runGitImpl,
    verifyWorkflowEvidenceImpl
  });
  const compact = JSON.stringify(outcome.plan);
  stdout.write(`${compact}\n`);
  if (env.GITHUB_OUTPUT) appendOutput(env.GITHUB_OUTPUT, `plan=${compact}\n`, appendFileImpl);
  if (env.GITHUB_STEP_SUMMARY) {
    appendOutput(
      env.GITHUB_STEP_SUMMARY,
      formatPlanSummary(outcome),
      appendFileImpl
    );
  }
};

const runGateCommand = ({ env, stdout, appendFileImpl }) => {
  const plan = parseJsonEnvironment(env, 'CI_IMPACT_PLAN_JSON');
  const needs = parseJsonEnvironment(env, 'CI_IMPACT_NEEDS_JSON');
  const expectedProfile = requiredEnvironment(env, 'CI_IMPACT_EXPECTED_PROFILE');
  const expectedWorkflowSha = objectIdEnvironment(env, 'CI_IMPACT_EXPECTED_WORKFLOW_SHA');
  const result = validateGateResults({
    plan,
    needs,
    knownJobs: knownJobsFor(expectedProfile),
    expected: { profile: expectedProfile, workflowSha: expectedWorkflowSha }
  });
  stdout.write(`${JSON.stringify(result)}\n`);
  if (env.GITHUB_STEP_SUMMARY) {
    appendOutput(env.GITHUB_STEP_SUMMARY, formatGateSummary({ plan, result }), appendFileImpl);
  }
};

const optionValues = (argv, names) => {
  const values = {};
  for (let index = 1; index < argv.length; index += 2) {
    const name = argv[index]?.replace(/^--/, '');
    if (!argv[index]?.startsWith('--') || !names.includes(name) || Object.hasOwn(values, name)
        || index + 1 >= argv.length || !argv[index + 1]) {
      return null;
    }
    values[name] = argv[index + 1];
  }
  return names.every((name) => Object.hasOwn(values, name)) ? values : null;
};

const runShardFilesCommand = ({ argv, env, cwd, stdout, appendFileImpl }) => {
  if (argv.length !== 2 || !/^[1-9][0-9]*$/.test(argv[1])) {
    fail('INVALID_ARGUMENTS', 'Usage: shard-files <shard number>.');
  }
  const plan = parseJsonEnvironment(env, 'CI_IMPACT_PLAN_JSON');
  const repositoryRoot = resolve(env.CI_IMPACT_REPOSITORY_ROOT || cwd);
  let shards;
  try {
    shards = JSON.parse(readFileSync(resolve(repositoryRoot, FUNCTIONAL_SHARDS_PATH), 'utf8')).shards;
  } catch (_error) {
    fail('INVALID_SHARD_MAP', 'The functional shard map could not be read.');
  }
  const files = functionalShardFiles({ plan, shards, shard: Number(argv[1]) }).join(' ');
  stdout.write(`${files}\n`);
  if (env.GITHUB_OUTPUT) appendOutput(env.GITHUB_OUTPUT, `files=${files}\n`, appendFileImpl);
};

const runClassifyCommand = ({ argv, env, cwd, stdout, runGitImpl }) => {
  const options = optionValues(argv, ['base', 'head']);
  if (!options) fail('INVALID_ARGUMENTS', 'Usage: classify --base <E> --head <H>.');
  const result = classifyRange({
    repositoryRoot: resolve(env.CI_IMPACT_REPOSITORY_ROOT || cwd),
    base: options.base,
    head: options.head,
    runGitImpl
  });
  stdout.write(`${JSON.stringify(result, null, 2)}\n`);
};

export const runCiImpactCli = async ({
  argv = process.argv.slice(2),
  env = process.env,
  cwd = process.cwd(),
  stdout = process.stdout,
  stderr = process.stderr,
  appendFileImpl = appendFileSync,
  runGitImpl = runGit,
  verifyWorkflowEvidenceImpl = verifyWorkflowEvidence
} = {}) => {
  try {
    const command = argv[0];
    if (!['plan', 'gate', 'shard-files', 'classify'].includes(command)
        || (['plan', 'gate'].includes(command) && argv.length !== 1)) {
      fail('INVALID_ARGUMENTS', 'Command must be plan, gate, shard-files, or classify.', {
        arguments: argv.map((argument) => boundedText(argument))
      });
    }
    if (command === 'shard-files') {
      runShardFilesCommand({ argv, env, cwd, stdout, appendFileImpl });
    } else if (command === 'classify') {
      runClassifyCommand({ argv, env, cwd, stdout, runGitImpl });
    } else if (command === 'plan') {
      await runPlanCommand({
        env,
        cwd,
        stdout,
        appendFileImpl,
        runGitImpl,
        verifyWorkflowEvidenceImpl
      });
    } else {
      runGateCommand({ env, stdout, appendFileImpl });
    }
    return 0;
  } catch (error) {
    const code = typeof error?.code === 'string' ? error.code : 'UNEXPECTED_ERROR';
    const message = error?.name === 'CiImpactPolicyError'
      || error?.name === 'CiImpactCliError'
      || error instanceof PromotionReadinessError
      ? error.message
      : 'Unexpected CI impact planner failure.';
    const details = isPlainObject(error?.details) ? error.details : {};
    const diagnostic = [
      `CI impact command failed [${code}]: ${message}`,
      `Details: ${JSON.stringify(details)}`,
      ''
    ].join('\n');
    stderr.write(redact(diagnostic, env));
    return 1;
  }
};

const isDirectExecution = process.argv[1]
  && import.meta.url === pathToFileURL(resolve(process.argv[1])).href;
if (isDirectExecution) process.exitCode = await runCiImpactCli();
