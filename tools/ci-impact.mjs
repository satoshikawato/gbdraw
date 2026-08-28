#!/usr/bin/env node

import { spawnSync } from 'node:child_process';
import { appendFileSync } from 'node:fs';
import { resolve } from 'node:path';
import { TextDecoder } from 'node:util';
import { pathToFileURL } from 'node:url';
import {
  PromotionReadinessError,
  verifyWorkflowEvidence
} from './check-promotion-readiness.mjs';
import {
  classifyChanges,
  createImpactPlan,
  isFullObjectId,
  knownJobsFor,
  validateGateResults
} from './ci-impact-policy.mjs';

const SUPPORTED_PROFILES = new Set(['pr', 'dev', 'gallery']);
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
  try {
    return Object.freeze({ valid: true, changes: parseNameStatusZ(result.stdout) });
  } catch (error) {
    return Object.freeze({
      valid: false,
      reason: error?.code || 'MALFORMED_GIT_DIFF',
      diagnostic: boundedText(error?.message || error)
    });
  }
};

const planFields = ({ configuration, classification, decision, basis, inheritedEvidence }) => ({
  profile: configuration.profile,
  impact: classification.impact,
  decision,
  basis,
  changeBaseSha: configuration.changeBaseSha,
  changeHeadSha: configuration.changeHeadSha,
  workflowSha: configuration.workflowSha,
  changedPathCount: classification.changedPathCount,
  inheritedEvidence
});

const fullClassification = (reason) => Object.freeze({
  impact: 'full',
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
  const classification = gitChanges.valid
    ? classifyChanges(gitChanges.changes)
    : fullClassification(gitChanges.reason);

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

  if (!classification.valid || classification.impact === 'full') {
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
  return Object.freeze({
    plan: createImpactPlan(planFields({
      configuration,
      classification,
      decision: 'selective',
      basis: configuration.profile === 'pr'
        ? 'LIGHT_CHANGE_WITH_DIRECT_BASE_EVIDENCE'
        : 'LIGHT_CHANGE_WITH_DIRECT_PARENT_EVIDENCE',
      inheritedEvidence
    })),
    classification
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
    : plan.profile === 'dev'
      ? 'active; dev staging jobs use the protected-branch plan'
      : 'observation only; existing test job conditions are unchanged';
  const lines = [
    `## CI impact plan${plan.profile === 'gallery' ? ' (shadow mode)' : ''}`,
    '',
    `- Profile: \`${plan.profile}\``,
    `- Impact: \`${plan.impact}\``,
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
    if (argv.length !== 1 || !['plan', 'gate'].includes(argv[0])) {
      fail('INVALID_ARGUMENTS', 'Command must be exactly plan or gate.', {
        arguments: argv.map((argument) => boundedText(argument))
      });
    }
    if (argv[0] === 'plan') {
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
