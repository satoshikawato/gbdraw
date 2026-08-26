#!/usr/bin/env node

import { spawnSync } from 'node:child_process';
import { appendFileSync } from 'node:fs';
import { resolve } from 'node:path';
import { pathToFileURL } from 'node:url';

const API_ORIGIN = 'https://api.github.com';
const API_VERSION = '2026-03-10';
const PER_PAGE = 100;
const MAX_PAGES = 10;
const MAX_SEARCH_RESULTS = PER_PAGE * MAX_PAGES;
const FULL_OBJECT_ID = /^(?:[0-9a-f]{40}|[0-9a-f]{64})$/i;
const REPOSITORY_NAME = /^[A-Za-z0-9](?:[A-Za-z0-9._-]{0,99})\/[A-Za-z0-9](?:[A-Za-z0-9._-]{0,99})$/;

export const PROMOTION_WORKFLOW_CONTRACTS = Object.freeze({
  '.github/workflows/gallery-publication.yml': 'Gallery readiness / gate',
  '.github/workflows/test.yml': 'Dev staging / gate'
});

const isObject = (value) => value !== null
  && typeof value === 'object'
  && !Array.isArray(value);

const boundedText = (value, limit = 240) => {
  if (typeof value !== 'string') return '';
  const normalized = value.replace(/\s+/g, ' ').trim();
  return normalized.length <= limit ? normalized : `${normalized.slice(0, limit)}...`;
};

export class PromotionReadinessError extends Error {
  constructor(code, reason, { expected = {}, observed = {} } = {}) {
    super(reason);
    this.name = 'PromotionReadinessError';
    this.code = code;
    this.expected = expected;
    this.observed = observed;
  }
}

const fail = (code, reason, details) => {
  throw new PromotionReadinessError(code, reason, details);
};

const positiveInteger = (value) => Number.isSafeInteger(value) && value > 0;

const normalizeObjectId = (value, field, expected = {}) => {
  if (typeof value !== 'string' || !FULL_OBJECT_ID.test(value)) {
    fail('INVALID_OBJECT_ID', `${field} must be a full 40- or 64-character hexadecimal object ID.`, {
      expected,
      observed: { [field]: typeof value === 'string' ? boundedText(value) : typeof value }
    });
  }
  return value.toLowerCase();
};

const repositoryParts = (repository) => {
  if (typeof repository !== 'string' || !REPOSITORY_NAME.test(repository)) {
    fail('INVALID_REPOSITORY', 'repository must be an exact OWNER/REPOSITORY full name.', {
      expected: { format: 'OWNER/REPOSITORY' },
      observed: { repository: boundedText(repository) }
    });
  }
  return repository.split('/').map(encodeURIComponent);
};

const workflowContract = (workflowPath, expectedAggregateName) => {
  if (!Object.hasOwn(PROMOTION_WORKFLOW_CONTRACTS, workflowPath)) {
    fail('UNSUPPORTED_WORKFLOW', 'workflow path is not a promotion evidence producer.', {
      expected: { workflowPaths: Object.keys(PROMOTION_WORKFLOW_CONTRACTS) },
      observed: { workflowPath: boundedText(workflowPath) }
    });
  }
  const aggregateName = PROMOTION_WORKFLOW_CONTRACTS[workflowPath];
  if (expectedAggregateName !== aggregateName) {
    fail('AGGREGATE_CONTRACT_MISMATCH', 'aggregate job name does not match the workflow contract.', {
      expected: { workflowPath, aggregateName },
      observed: { expectedAggregateName: boundedText(expectedAggregateName) }
    });
  }
  return aggregateName;
};

const responseHeader = (response, name) => {
  if (response.headers && typeof response.headers.get === 'function') {
    return response.headers.get(name);
  }
  if (isObject(response.headers)) {
    const entry = Object.entries(response.headers)
      .find(([key]) => key.toLowerCase() === name.toLowerCase());
    return entry?.[1] ?? null;
  }
  return null;
};

const requestJson = async ({ url, token, fetchImpl, expected, notFoundCode = '' }) => {
  let response;
  try {
    response = await fetchImpl(url, {
      method: 'GET',
      headers: {
        Accept: 'application/vnd.github+json',
        Authorization: `Bearer ${token.trim()}`,
        'X-GitHub-Api-Version': API_VERSION
      }
    });
  } catch (error) {
    fail('API_REQUEST_FAILED', 'GitHub Actions API request failed.', {
      expected,
      observed: { errorType: boundedText(error?.name || typeof error) }
    });
  }
  if (!isObject(response) || typeof response.ok !== 'boolean'
      || !Number.isInteger(response.status) || typeof response.json !== 'function') {
    fail('MALFORMED_API_RESPONSE', 'GitHub Actions API returned a malformed response envelope.', {
      expected,
      observed: { responseType: typeof response }
    });
  }

  let body;
  try {
    body = await response.json();
  } catch (_error) {
    fail('MALFORMED_API_RESPONSE', 'GitHub Actions API response was not valid JSON.', {
      expected,
      observed: { status: response.status }
    });
  }
  if (!response.ok) {
    const code = response.status === 404 && notFoundCode ? notFoundCode : 'API_REQUEST_FAILED';
    fail(code, response.status === 404 && notFoundCode
      ? 'Expected GitHub Actions resource was not found.'
      : 'GitHub Actions API request did not succeed.', {
      expected,
      observed: {
        status: response.status,
        message: boundedText(isObject(body) ? body.message : '')
      }
    });
  }
  return { body, link: responseHeader(response, 'link') };
};

export const parseLinkHeader = (source) => {
  if (source === null || source === undefined) return new Map();
  if (typeof source !== 'string' || !source.trim()) {
    fail('MALFORMED_PAGINATION', 'GitHub pagination Link header is malformed.', {
      observed: { link: boundedText(source) }
    });
  }
  const relations = new Map();
  for (const entry of source.split(/,\s*(?=<)/)) {
    const match = entry.trim().match(/^<([^<>]+)>;\s*rel="(next|prev|first|last)"$/);
    if (!match || relations.has(match?.[2])) {
      fail('MALFORMED_PAGINATION', 'GitHub pagination Link header is malformed.', {
        observed: { link: boundedText(source) }
      });
    }
    relations.set(match[2], match[1]);
  }
  return relations;
};

const validatePaginationUrl = (source, { pathname, immutableParameters, expected }) => {
  let url;
  try {
    url = new URL(source);
  } catch (_error) {
    fail('MALFORMED_PAGINATION', 'GitHub pagination next URL is invalid.', {
      expected,
      observed: { nextUrl: boundedText(source) }
    });
  }
  const allowedParameters = new Set([...immutableParameters.keys(), 'page']);
  const allParameterNames = [...url.searchParams.keys()];
  const parameterNames = [...new Set(allParameterNames)];
  const page = url.searchParams.get('page');
  if (url.origin !== API_ORIGIN || url.pathname !== pathname || url.username || url.password
      || url.hash || parameterNames.some((name) => !allowedParameters.has(name))
      || parameterNames.length !== allParameterNames.length
      || [...immutableParameters].some(([name, value]) => url.searchParams.get(name) !== value)
      || (page !== null && (!/^\d+$/.test(page) || !Number.isSafeInteger(Number(page))
        || Number(page) < 1))) {
    fail('MALFORMED_PAGINATION', 'GitHub pagination next URL changed the evidence query.', {
      expected,
      observed: { nextUrl: boundedText(source) }
    });
  }
  return url.toString();
};

const collectPaginated = async ({
  initialUrl,
  pathname,
  immutableParameters,
  arrayField,
  token,
  fetchImpl,
  expected
}) => {
  let nextUrl = initialUrl;
  let totalCount = null;
  const pages = new Set();
  const itemIds = new Set();
  const items = [];

  while (nextUrl) {
    if (pages.size >= MAX_PAGES) {
      fail('INCOMPLETE_PAGINATION', 'GitHub evidence exceeded the bounded pagination limit.', {
        expected,
        observed: { pages: pages.size, limit: MAX_PAGES }
      });
    }
    const validatedUrl = validatePaginationUrl(nextUrl, {
      pathname,
      immutableParameters,
      expected
    });
    if (pages.has(validatedUrl)) {
      fail('PAGINATION_LOOP', 'GitHub evidence pagination repeated a page.', {
        expected,
        observed: { nextUrl: boundedText(validatedUrl), pages: pages.size }
      });
    }
    const page = Number(new URL(validatedUrl).searchParams.get('page') || 1);
    if (page !== pages.size + 1) {
      fail('INCOMPLETE_PAGINATION', 'GitHub evidence pagination skipped or reordered a page.', {
        expected,
        observed: { expectedPage: pages.size + 1, page }
      });
    }
    pages.add(validatedUrl);

    const { body, link } = await requestJson({
      url: validatedUrl,
      token,
      fetchImpl,
      expected
    });
    if (!isObject(body) || !Number.isSafeInteger(body.total_count) || body.total_count < 0
        || !Array.isArray(body[arrayField])) {
      fail('MALFORMED_API_RESPONSE', `GitHub ${arrayField} response shape is malformed.`, {
        expected,
        observed: {
          totalCount: body?.total_count,
          collectionType: Array.isArray(body?.[arrayField]) ? 'array' : typeof body?.[arrayField]
        }
      });
    }
    if (body.total_count > MAX_SEARCH_RESULTS) {
      fail('INCOMPLETE_PAGINATION', 'GitHub evidence result count exceeds the bounded search limit.', {
        expected,
        observed: { totalCount: body.total_count, limit: MAX_SEARCH_RESULTS }
      });
    }
    if (totalCount === null) totalCount = body.total_count;
    else if (totalCount !== body.total_count) {
      fail('INCOMPLETE_PAGINATION', 'GitHub evidence total changed during pagination.', {
        expected,
        observed: { initialTotalCount: totalCount, currentTotalCount: body.total_count }
      });
    }

    for (const item of body[arrayField]) {
      if (!isObject(item) || !positiveInteger(item.id) || itemIds.has(item.id)) {
        fail('MALFORMED_API_RESPONSE', `GitHub ${arrayField} contains malformed or duplicate items.`, {
          expected,
          observed: { itemId: item?.id }
        });
      }
      itemIds.add(item.id);
      items.push(item);
    }
    if (items.length > totalCount) {
      fail('INCOMPLETE_PAGINATION', 'GitHub evidence returned more items than total_count.', {
        expected,
        observed: { totalCount, received: items.length }
      });
    }

    const relations = parseLinkHeader(link);
    nextUrl = relations.get('next') || '';
    if (nextUrl && items.length >= totalCount) {
      fail('INCOMPLETE_PAGINATION', 'GitHub evidence advertised another page after total_count.', {
        expected,
        observed: { totalCount, received: items.length }
      });
    }
  }

  if (items.length !== totalCount) {
    fail('INCOMPLETE_PAGINATION', 'GitHub evidence pagination ended before total_count.', {
      expected,
      observed: { totalCount, received: items.length }
    });
  }
  return items;
};

const validTimestamp = (value) => typeof value === 'string'
  && value.length > 0
  && Number.isFinite(Date.parse(value));

const runIdentity = (run) => ({
  id: run?.id,
  runNumber: run?.run_number,
  runAttempt: run?.run_attempt,
  workflowId: run?.workflow_id,
  workflowPath: run?.path,
  repository: run?.repository?.full_name,
  event: run?.event,
  branch: run?.head_branch,
  headSha: run?.head_sha,
  status: run?.status,
  conclusion: run?.conclusion
});

const validateRun = (run, expected) => {
  const observed = runIdentity(run);
  if (!isObject(run) || !positiveInteger(run.id) || !positiveInteger(run.run_number)
      || !positiveInteger(run.run_attempt) || !positiveInteger(run.workflow_id)
      || typeof run.path !== 'string' || typeof run.event !== 'string'
      || typeof run.head_branch !== 'string' || typeof run.head_sha !== 'string'
      || !isObject(run.repository) || typeof run.repository.full_name !== 'string'
      || typeof run.status !== 'string'
      || !(run.conclusion === null || typeof run.conclusion === 'string')
      || typeof run.url !== 'string' || typeof run.html_url !== 'string'
      || !validTimestamp(run.created_at) || !validTimestamp(run.updated_at)) {
    fail('MALFORMED_RUN', 'GitHub workflow run metadata is malformed.', { expected, observed });
  }
  return run;
};

const runMatches = (run, expected) => run.workflow_id === expected.workflowId
  && run.path === expected.workflowPath
  && run.repository.full_name === expected.repository
  && run.event === expected.event
  && run.head_branch === expected.branch
  && run.head_sha.toLowerCase() === expected.headSha;

const validateRunUrls = (run, repositoryApiPath, expected) => {
  const expectedApiUrl = `${API_ORIGIN}${repositoryApiPath}/actions/runs/${run.id}`;
  const expectedHtmlUrl = `https://github.com/${expected.repository}/actions/runs/${run.id}`;
  if (run.url !== expectedApiUrl || run.html_url !== expectedHtmlUrl) {
    fail('RUN_IDENTITY_MISMATCH', 'GitHub workflow run URLs do not match its repository and run ID.', {
      expected: { ...expected, runApiUrl: expectedApiUrl, runUrl: expectedHtmlUrl },
      observed: { ...runIdentity(run), runApiUrl: run.url, runUrl: run.html_url }
    });
  }
};

const requireSuccessfulRun = (run, expected) => {
  if (run.status !== 'completed' || run.conclusion !== 'success') {
    fail('RUN_NOT_SUCCESSFUL', 'Newest exact-identity workflow run is not completed successfully.', {
      expected: { ...expected, status: 'completed', conclusion: 'success' },
      observed: runIdentity(run)
    });
  }
};

const validateJob = (job, run, expected, repositoryApiPath) => {
  const observed = {
    id: job?.id,
    runId: job?.run_id,
    name: job?.name,
    headSha: job?.head_sha,
    status: job?.status,
    conclusion: job?.conclusion
  };
  const expectedApiUrl = `${API_ORIGIN}${repositoryApiPath}/actions/jobs/${job?.id}`;
  const expectedHtmlUrl = `https://github.com/${expected.repository}/actions/runs/${run.id}/job/${job?.id}`;
  if (!isObject(job) || !positiveInteger(job.id) || job.run_id !== run.id
      || typeof job.name !== 'string' || typeof job.head_sha !== 'string'
      || job.head_sha.toLowerCase() !== expected.headSha || typeof job.status !== 'string'
      || !(job.conclusion === null || typeof job.conclusion === 'string')
      || job.url !== expectedApiUrl || typeof job.html_url !== 'string'
      || job.html_url !== expectedHtmlUrl) {
    fail('MALFORMED_JOB', 'GitHub workflow job metadata is malformed or identity-mismatched.', {
      expected: { ...expected, runId: run.id },
      observed
    });
  }
  return job;
};

const exactRunSnapshot = ({ run, selectedRun, expected, repositoryApiPath }) => {
  validateRun(run, expected);
  validateRunUrls(run, repositoryApiPath, expected);
  if (!runMatches(run, expected) || run.id !== selectedRun.id
      || run.run_number !== selectedRun.run_number
      || run.run_attempt < selectedRun.run_attempt) {
    fail('RUN_IDENTITY_MISMATCH', 'Current workflow run no longer matches the selected evidence identity.', {
      expected: { ...expected, runId: selectedRun.id, minimumAttempt: selectedRun.run_attempt },
      observed: runIdentity(run)
    });
  }
  return run;
};

export const verifyWorkflowEvidence = async ({
  repository,
  workflowPath,
  expectedEvent = 'push',
  expectedBranch = 'dev',
  expectedHeadSha,
  expectedAggregateName,
  token,
  fetchImpl = globalThis.fetch
}) => {
  const [owner, repo] = repositoryParts(repository);
  workflowContract(workflowPath, expectedAggregateName);
  const headSha = normalizeObjectId(expectedHeadSha, 'expectedHeadSha');
  const expectedIdentity = {
    repository,
    workflowPath,
    event: 'push',
    branch: 'dev',
    headSha,
    aggregateName: expectedAggregateName
  };
  if (expectedEvent !== 'push' || expectedBranch !== 'dev') {
    fail('UNSUPPORTED_STAGING_IDENTITY', 'Promotion evidence is limited to push runs on dev.', {
      expected: expectedIdentity,
      observed: { event: expectedEvent, branch: expectedBranch }
    });
  }
  if (typeof token !== 'string' || !token.trim()) {
    fail('MISSING_TOKEN', 'A GitHub token with Actions read access is required.', {
      expected: expectedIdentity,
      observed: { token: 'missing' }
    });
  }
  if (typeof fetchImpl !== 'function') {
    fail('MISSING_FETCH', 'A fetch implementation is required.', { expected: expectedIdentity });
  }

  const repositoryApiPath = `/repos/${owner}/${repo}`;
  const workflowResolutionUrl = `${API_ORIGIN}${repositoryApiPath}/actions/workflows/${encodeURIComponent(workflowPath)}`;
  const { body: workflow } = await requestJson({
    url: workflowResolutionUrl,
    token,
    fetchImpl,
    expected: expectedIdentity,
    notFoundCode: 'WORKFLOW_NOT_FOUND'
  });
  const workflowObserved = {
    id: workflow?.id,
    path: workflow?.path,
    state: workflow?.state,
    url: workflow?.url
  };
  if (!isObject(workflow) || !positiveInteger(workflow.id)
      || workflow.path !== workflowPath || workflow.state !== 'active'
      || workflow.url !== `${API_ORIGIN}${repositoryApiPath}/actions/workflows/${workflow.id}`) {
    fail('WORKFLOW_IDENTITY_MISMATCH', 'Resolved workflow is disabled or does not match the exact workflow path.', {
      expected: expectedIdentity,
      observed: workflowObserved
    });
  }

  const runExpected = { ...expectedIdentity, workflowId: workflow.id };
  const runsPath = `${repositoryApiPath}/actions/workflows/${workflow.id}/runs`;
  const runParameters = new Map([
    ['branch', 'dev'],
    ['event', 'push'],
    ['exclude_pull_requests', 'true'],
    ['head_sha', headSha],
    ['per_page', String(PER_PAGE)]
  ]);
  const runsUrl = new URL(`${API_ORIGIN}${runsPath}`);
  runParameters.forEach((value, name) => runsUrl.searchParams.set(name, value));
  const runs = await collectPaginated({
    initialUrl: runsUrl.toString(),
    pathname: runsPath,
    immutableParameters: runParameters,
    arrayField: 'workflow_runs',
    token,
    fetchImpl,
    expected: runExpected
  });
  runs.forEach((run) => {
    validateRun(run, runExpected);
    validateRunUrls(run, repositoryApiPath, runExpected);
  });
  const matchingRuns = runs.filter((run) => runMatches(run, runExpected)).sort((left, right) => (
    Date.parse(right.created_at) - Date.parse(left.created_at)
    || right.run_number - left.run_number
    || right.id - left.id
  ));
  if (!matchingRuns.length) {
    fail('NO_MATCHING_RUN', 'No workflow run matches the exact promotion evidence identity.', {
      expected: runExpected,
      observed: { runs: runs.slice(0, 5).map(runIdentity), observedCount: runs.length }
    });
  }
  const selectedRun = matchingRuns[0];
  const currentRunUrl = `${API_ORIGIN}${repositoryApiPath}/actions/runs/${selectedRun.id}`;
  const { body: currentRunBody } = await requestJson({
    url: currentRunUrl,
    token,
    fetchImpl,
    expected: { ...runExpected, runId: selectedRun.id }
  });
  const currentRun = exactRunSnapshot({
    run: currentRunBody,
    selectedRun,
    expected: runExpected,
    repositoryApiPath
  });
  requireSuccessfulRun(currentRun, runExpected);

  const jobsPath = `${repositoryApiPath}/actions/runs/${currentRun.id}/attempts/${currentRun.run_attempt}/jobs`;
  const jobParameters = new Map([['per_page', String(PER_PAGE)]]);
  const jobsUrl = new URL(`${API_ORIGIN}${jobsPath}`);
  jobParameters.forEach((value, name) => jobsUrl.searchParams.set(name, value));
  const jobs = await collectPaginated({
    initialUrl: jobsUrl.toString(),
    pathname: jobsPath,
    immutableParameters: jobParameters,
    arrayField: 'jobs',
    token,
    fetchImpl,
    expected: { ...runExpected, runId: currentRun.id, runAttempt: currentRun.run_attempt }
  });
  jobs.forEach((job) => validateJob(job, currentRun, runExpected, repositoryApiPath));
  const aggregateJobs = jobs.filter((job) => job.name === expectedAggregateName);
  if (aggregateJobs.length !== 1) {
    fail('AGGREGATE_JOB_CARDINALITY', 'Expected exactly one aggregate job on the current run attempt.', {
      expected: {
        ...runExpected,
        runId: currentRun.id,
        runAttempt: currentRun.run_attempt,
        aggregateMatches: 1
      },
      observed: {
        aggregateMatches: aggregateJobs.length,
        jobs: jobs.slice(0, 5).map(({ id, name, status, conclusion }) => ({
          id, name, status, conclusion
        }))
      }
    });
  }
  const aggregateJob = aggregateJobs[0];
  if (aggregateJob.status !== 'completed' || aggregateJob.conclusion !== 'success') {
    fail('AGGREGATE_JOB_NOT_SUCCESSFUL', 'Aggregate job is not completed successfully.', {
      expected: { ...runExpected, status: 'completed', conclusion: 'success' },
      observed: {
        id: aggregateJob.id,
        name: aggregateJob.name,
        status: aggregateJob.status,
        conclusion: aggregateJob.conclusion
      }
    });
  }

  const { body: finalRunBody } = await requestJson({
    url: currentRunUrl,
    token,
    fetchImpl,
    expected: { ...runExpected, runId: currentRun.id, runAttempt: currentRun.run_attempt }
  });
  const finalRun = exactRunSnapshot({
    run: finalRunBody,
    selectedRun: currentRun,
    expected: runExpected,
    repositoryApiPath
  });
  if (finalRun.run_attempt !== currentRun.run_attempt) {
    fail('RUN_ATTEMPT_CHANGED', 'Workflow run attempt changed during evidence verification.', {
      expected: { ...runExpected, runId: currentRun.id, runAttempt: currentRun.run_attempt },
      observed: runIdentity(finalRun)
    });
  }
  requireSuccessfulRun(finalRun, runExpected);

  return {
    kind: 'promotion-workflow-evidence',
    repository,
    workflow: {
      path: workflowPath,
      id: workflow.id,
      state: workflow.state
    },
    run: {
      id: finalRun.id,
      number: finalRun.run_number,
      attempt: finalRun.run_attempt,
      url: finalRun.html_url,
      apiUrl: finalRun.url,
      createdAt: finalRun.created_at,
      updatedAt: finalRun.updated_at,
      event: finalRun.event,
      branch: finalRun.head_branch,
      headSha,
      status: finalRun.status,
      conclusion: finalRun.conclusion
    },
    aggregateJob: {
      id: aggregateJob.id,
      name: aggregateJob.name,
      url: aggregateJob.html_url,
      apiUrl: aggregateJob.url,
      status: aggregateJob.status,
      conclusion: aggregateJob.conclusion
    },
    selection: {
      basis: 'newest exact-identity run by created_at, then run_number, then run ID; current attempt revalidated after job inspection',
      matchingRunIds: matchingRuns.map(({ id }) => id),
      ignoredOlderRunIds: matchingRuns.slice(1).map(({ id }) => id)
    }
  };
};

export const runGit = (repositoryRoot, args) => spawnSync('git', args, {
  cwd: repositoryRoot,
  encoding: 'utf8',
  stdio: ['ignore', 'pipe', 'pipe']
});

const gitProbe = (runGitImpl, repositoryRoot, args, expected) => {
  const result = runGitImpl(repositoryRoot, args);
  if (!isObject(result) || !Number.isInteger(result.status)
      || typeof result.stdout !== 'string' || typeof result.stderr !== 'string') {
    fail('MALFORMED_GIT_RESULT', 'Injected Git runner returned a malformed result.', {
      expected,
      observed: { args, resultType: typeof result }
    });
  }
  return result;
};

export const verifyPromotionResultTree = ({
  repositoryRoot,
  baseSha,
  headSha,
  runGit: runGitImpl = runGit
}) => {
  const expected = {
    baseSha: normalizeObjectId(baseSha, 'baseSha'),
    headSha: normalizeObjectId(headSha, 'headSha')
  };
  if (typeof repositoryRoot !== 'string' || !repositoryRoot.trim()) {
    fail('INVALID_REPOSITORY_ROOT', 'repositoryRoot must be a nonempty path.', { expected });
  }
  if (typeof runGitImpl !== 'function') {
    fail('MISSING_GIT_RUNNER', 'A Git runner is required.', { expected });
  }

  const shallow = gitProbe(
    runGitImpl,
    repositoryRoot,
    ['rev-parse', '--is-shallow-repository'],
    expected
  );
  if (shallow.status !== 0 || shallow.stdout.trim() !== 'false') {
    fail('INCOMPLETE_GIT_HISTORY', 'Complete Git history is required for promotion tree proof.', {
      expected,
      observed: {
        shallow: shallow.stdout.trim() || 'unknown',
        error: boundedText(shallow.stderr)
      }
    });
  }
  for (const [name, sha] of [['baseSha', expected.baseSha], ['headSha', expected.headSha]]) {
    const object = gitProbe(
      runGitImpl,
      repositoryRoot,
      ['cat-file', '-e', `${sha}^{commit}`],
      expected
    );
    if (object.status !== 0) {
      fail('MISSING_GIT_OBJECT', `Complete Git object for ${name} is unavailable.`, {
        expected,
        observed: { missing: sha, error: boundedText(object.stderr) }
      });
    }
  }

  const headTreeResult = gitProbe(
    runGitImpl,
    repositoryRoot,
    ['rev-parse', '--verify', `${expected.headSha}^{tree}`],
    expected
  );
  const headTreeSha = headTreeResult.stdout.trim().toLowerCase();
  if (headTreeResult.status !== 0 || !FULL_OBJECT_ID.test(headTreeSha)) {
    fail('MISSING_GIT_OBJECT', 'Head tree object could not be resolved.', {
      expected,
      observed: { error: boundedText(headTreeResult.stderr) }
    });
  }

  const mergeResult = gitProbe(
    runGitImpl,
    repositoryRoot,
    ['merge-tree', '--write-tree', '--no-messages', expected.baseSha, expected.headSha],
    expected
  );
  if (mergeResult.status !== 0) {
    const error = boundedText(mergeResult.stderr || mergeResult.stdout);
    if (/missing|not a valid object|could not read|unable to read (?:tree|object)|promisor/i.test(error)) {
      fail('MISSING_GIT_OBJECT', 'Synthetic merge could not read complete Git history.', {
        expected,
        observed: { error }
      });
    }
    if (/read-only file system|unable to create temporary file|permission denied/i.test(error)) {
      fail('GIT_MERGE_TREE_UNAVAILABLE', 'Git could not create the temporary synthetic tree object.', {
        expected,
        observed: { error }
      });
    }
    fail('MERGE_CONFLICT', 'Synthetic promotion merge is not conflict-free.', {
      expected,
      observed: { error }
    });
  }
  const mergeLines = mergeResult.stdout.trim().split(/\r?\n/).filter(Boolean);
  const syntheticMergeTreeSha = mergeLines.length === 1
    ? mergeLines[0].toLowerCase()
    : '';
  if (!FULL_OBJECT_ID.test(syntheticMergeTreeSha)) {
    fail('MALFORMED_MERGE_TREE', 'Synthetic merge did not return one full tree object ID.', {
      expected,
      observed: { output: boundedText(mergeResult.stdout) }
    });
  }
  if (syntheticMergeTreeSha !== headTreeSha) {
    fail('RESULT_TREE_MISMATCH', 'Synthetic promotion merge tree differs from the dev head tree.', {
      expected: { ...expected, headTreeSha },
      observed: { syntheticMergeTreeSha }
    });
  }

  return {
    kind: 'promotion-result-tree-evidence',
    repositoryRoot: resolve(repositoryRoot),
    baseSha: expected.baseSha,
    headSha: expected.headSha,
    syntheticMergeTreeSha,
    headTreeSha,
    conflictFree: true,
    treeMatchesHead: true,
    method: 'git merge-tree --write-tree --no-messages'
  };
};

const parseArguments = (argv) => {
  const [command, ...rest] = argv;
  if (!['workflow-evidence', 'merge-tree'].includes(command)) {
    fail('INVALID_COMMAND', 'Command must be workflow-evidence or merge-tree.', {
      expected: { commands: ['workflow-evidence', 'merge-tree'] },
      observed: { command: boundedText(command) }
    });
  }
  const options = new Map();
  for (let index = 0; index < rest.length; index += 2) {
    const name = rest[index];
    const value = rest[index + 1];
    if (!name?.startsWith('--') || !value || value.startsWith('--') || options.has(name)) {
      fail('INVALID_ARGUMENTS', 'CLI arguments must be unique --name value pairs.', {
        observed: { argument: boundedText(name) }
      });
    }
    options.set(name, value);
  }
  return { command, options };
};

const requiredOption = (options, name) => {
  const value = options.get(name);
  if (!value) {
    fail('MISSING_ARGUMENT', `Missing required CLI argument ${name}.`, {
      expected: { argument: name }
    });
  }
  return value;
};

export const formatEvidenceSummary = (evidence) => {
  if (evidence.kind === 'promotion-result-tree-evidence') {
    return [
      '## Promotion result tree evidence',
      '',
      `- Base SHA: \`${evidence.baseSha}\``,
      `- Head SHA: \`${evidence.headSha}\``,
      `- Synthetic/head tree: \`${evidence.headTreeSha}\``,
      '- Result: conflict-free and tree-identical',
      ''
    ].join('\n');
  }
  return [
    '## Promotion workflow evidence',
    '',
    `- Repository: \`${evidence.repository}\``,
    `- Workflow: \`${evidence.workflow.path}\` (ID ${evidence.workflow.id})`,
    `- Run: [${evidence.run.id}](${evidence.run.url}), attempt ${evidence.run.attempt}`,
    `- Identity: \`${evidence.run.event}\` / \`${evidence.run.branch}\` / \`${evidence.run.headSha}\``,
    `- Aggregate: [${evidence.aggregateJob.name}](${evidence.aggregateJob.url})`,
    '- Result: run and aggregate completed successfully',
    ''
  ].join('\n');
};

export const runPromotionReadinessCli = async ({
  argv = process.argv.slice(2),
  env = process.env,
  cwd = process.cwd(),
  fetchImpl = globalThis.fetch,
  runGit: runGitImpl = runGit,
  stdout = process.stdout,
  stderr = process.stderr
} = {}) => {
  try {
    const { command, options } = parseArguments(argv);
    let evidence;
    if (command === 'workflow-evidence') {
      const allowed = new Set([
        '--aggregate-name', '--head-sha', '--repository', '--workflow-path'
      ]);
      if ([...options.keys()].some((name) => !allowed.has(name))) {
        fail('INVALID_ARGUMENTS', 'Unknown workflow-evidence CLI argument.', {
          observed: { arguments: [...options.keys()] }
        });
      }
      evidence = await verifyWorkflowEvidence({
        repository: requiredOption(options, '--repository'),
        workflowPath: requiredOption(options, '--workflow-path'),
        expectedHeadSha: requiredOption(options, '--head-sha'),
        expectedAggregateName: requiredOption(options, '--aggregate-name'),
        token: env.GITHUB_TOKEN,
        fetchImpl
      });
    } else {
      const allowed = new Set(['--base-sha', '--head-sha', '--repository-root']);
      if ([...options.keys()].some((name) => !allowed.has(name))) {
        fail('INVALID_ARGUMENTS', 'Unknown merge-tree CLI argument.', {
          observed: { arguments: [...options.keys()] }
        });
      }
      evidence = verifyPromotionResultTree({
        repositoryRoot: options.get('--repository-root') || cwd,
        baseSha: requiredOption(options, '--base-sha'),
        headSha: requiredOption(options, '--head-sha'),
        runGit: runGitImpl
      });
    }
    stdout.write(`${JSON.stringify(evidence, null, 2)}\n`);
    if (env.GITHUB_STEP_SUMMARY) {
      appendFileSync(resolve(env.GITHUB_STEP_SUMMARY), formatEvidenceSummary(evidence), 'utf8');
    }
    return 0;
  } catch (error) {
    const failure = error instanceof PromotionReadinessError
      ? error
      : new PromotionReadinessError('UNEXPECTED_ERROR', 'Unexpected promotion verification failure.', {
        observed: { error: boundedText(error?.message || String(error)) }
      });
    const diagnostic = [
      `Promotion readiness verification failed [${failure.code}]: ${failure.message}`,
      `Expected: ${JSON.stringify(failure.expected)}`,
      `Observed: ${JSON.stringify(failure.observed)}`,
      ''
    ].join('\n');
    const token = typeof env.GITHUB_TOKEN === 'string' ? env.GITHUB_TOKEN : '';
    stderr.write(token ? diagnostic.replaceAll(token, '[REDACTED]') : diagnostic);
    return 1;
  }
};

const isDirectExecution = process.argv[1]
  && import.meta.url === pathToFileURL(resolve(process.argv[1])).href;
if (isDirectExecution) process.exitCode = await runPromotionReadinessCli();
