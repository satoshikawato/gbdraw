import assert from 'node:assert/strict';
import test from 'node:test';

import {
  PromotionReadinessError,
  runPromotionReadinessCli,
  verifyPromotionResultTree,
  verifyWorkflowEvidence
} from '../../tools/check-promotion-readiness.mjs';

const API_ORIGIN = 'https://api.github.com';
const REPOSITORY = 'satoshikawato/gbdraw';
const REPOSITORY_API_PATH = '/repos/satoshikawato/gbdraw';
const HEAD_SHA = '2'.repeat(40);
const BASE_SHA = '1'.repeat(40);
const HEAD_TREE_SHA = '3'.repeat(40);
const OTHER_TREE_SHA = '4'.repeat(40);
const TEST_WORKFLOW_PATH = '.github/workflows/test.yml';
const TEST_AGGREGATE = 'Dev staging / gate';
const GALLERY_WORKFLOW_PATH = '.github/workflows/gallery-publication.yml';
const GALLERY_AGGREGATE = 'Gallery readiness / gate';

const workflowFixture = ({ path = TEST_WORKFLOW_PATH, id = 501, state = 'active' } = {}) => ({
  id,
  path,
  state,
  url: `${API_ORIGIN}${REPOSITORY_API_PATH}/actions/workflows/${id}`
});

const runFixture = (overrides = {}) => {
  const id = overrides.id ?? 701;
  return {
    id,
    run_number: 41,
    run_attempt: 1,
    workflow_id: 501,
    path: TEST_WORKFLOW_PATH,
    event: 'push',
    head_branch: 'dev',
    head_sha: HEAD_SHA,
    repository: { full_name: REPOSITORY },
    status: 'completed',
    conclusion: 'success',
    url: `${API_ORIGIN}${REPOSITORY_API_PATH}/actions/runs/${id}`,
    html_url: `https://github.com/${REPOSITORY}/actions/runs/${id}`,
    created_at: '2026-08-27T10:00:00Z',
    updated_at: '2026-08-27T10:10:00Z',
    ...overrides
  };
};

const jobFixture = (overrides = {}) => {
  const id = overrides.id ?? 801;
  const runId = overrides.run_id ?? 701;
  return {
    id,
    run_id: runId,
    head_sha: HEAD_SHA,
    name: TEST_AGGREGATE,
    status: 'completed',
    conclusion: 'success',
    url: `${API_ORIGIN}${REPOSITORY_API_PATH}/actions/jobs/${id}`,
    html_url: `https://github.com/${REPOSITORY}/actions/runs/${runId}/job/${id}`,
    ...overrides
  };
};

const response = (body, { status = 200, link = null } = {}) => ({
  ok: status >= 200 && status < 300,
  status,
  headers: { get: (name) => name.toLowerCase() === 'link' ? link : null },
  json: async () => body
});

const nextLink = (source, page) => {
  const url = new URL(source);
  url.searchParams.set('page', String(page));
  return `<${url}>; rel="next"`;
};

const createFetch = ({
  workflow = workflowFixture(),
  workflowStatus = 200,
  runPages = [[runFixture()]],
  runTotalCount = runPages.reduce((count, page) => count + page.length, 0),
  runLinks = [],
  currentRuns = [runFixture()],
  jobPages = [[jobFixture()]],
  jobTotalCount = jobPages.reduce((count, page) => count + page.length, 0),
  jobLinks = []
} = {}) => {
  const calls = [];
  let currentRunIndex = 0;
  const fetchImpl = async (source, options) => {
    const url = new URL(source);
    calls.push({ url: url.toString(), options });
    assert.equal(options.headers.Authorization, 'Bearer fixture-token');

    if (/\/actions\/workflows\/[^/]+$/.test(url.pathname)) {
      return response(
        workflowStatus === 200 ? workflow : { message: 'Not Found' },
        { status: workflowStatus }
      );
    }
    if (/\/actions\/workflows\/\d+\/runs$/.test(url.pathname)) {
      const page = Number(url.searchParams.get('page') || 1);
      const link = Object.hasOwn(runLinks, page - 1)
        ? runLinks[page - 1]
        : page < runPages.length ? nextLink(url, page + 1) : null;
      return response({
        total_count: runTotalCount,
        workflow_runs: runPages[page - 1] ?? []
      }, { link });
    }
    if (/\/actions\/runs\/\d+\/attempts\/\d+\/jobs$/.test(url.pathname)) {
      const page = Number(url.searchParams.get('page') || 1);
      const link = Object.hasOwn(jobLinks, page - 1)
        ? jobLinks[page - 1]
        : page < jobPages.length ? nextLink(url, page + 1) : null;
      return response({
        total_count: jobTotalCount,
        jobs: jobPages[page - 1] ?? []
      }, { link });
    }
    if (/\/actions\/runs\/\d+$/.test(url.pathname)) {
      const snapshot = currentRuns[Math.min(currentRunIndex, currentRuns.length - 1)];
      currentRunIndex += 1;
      return response(snapshot);
    }
    assert.fail(`Unexpected GitHub API request: ${url}`);
  };
  fetchImpl.calls = calls;
  return fetchImpl;
};

const verify = (overrides = {}) => verifyWorkflowEvidence({
  repository: REPOSITORY,
  workflowPath: TEST_WORKFLOW_PATH,
  expectedEvent: 'push',
  expectedBranch: 'dev',
  expectedHeadSha: HEAD_SHA,
  expectedAggregateName: TEST_AGGREGATE,
  token: 'fixture-token',
  fetchImpl: createFetch(),
  ...overrides
});

const rejectsWith = async (promise, code, reason) => assert.rejects(promise, (error) => {
  assert.ok(error instanceof PromotionReadinessError);
  assert.equal(error.code, code);
  assert.match(error.message, reason);
  assert.ok(Number.isFinite(JSON.stringify(error.observed).length));
  return true;
});

test('exact Tests and Gallery workflow evidence succeeds with auditable output', async () => {
  const cases = [
    [TEST_WORKFLOW_PATH, TEST_AGGREGATE, 501, 701, 801],
    [GALLERY_WORKFLOW_PATH, GALLERY_AGGREGATE, 502, 702, 802]
  ];
  for (const [workflowPath, aggregateName, workflowId, runId, jobId] of cases) {
    const workflow = workflowFixture({ path: workflowPath, id: workflowId });
    const run = runFixture({
      id: runId,
      workflow_id: workflowId,
      path: workflowPath,
      url: `${API_ORIGIN}${REPOSITORY_API_PATH}/actions/runs/${runId}`,
      html_url: `https://github.com/${REPOSITORY}/actions/runs/${runId}`
    });
    const job = jobFixture({
      id: jobId,
      run_id: runId,
      name: aggregateName,
      url: `${API_ORIGIN}${REPOSITORY_API_PATH}/actions/jobs/${jobId}`,
      html_url: `https://github.com/${REPOSITORY}/actions/runs/${runId}/job/${jobId}`
    });
    const fetchImpl = createFetch({
      workflow,
      runPages: [[run]],
      currentRuns: [run],
      jobPages: [[job]]
    });
    const evidence = await verify({
      workflowPath,
      expectedAggregateName: aggregateName,
      fetchImpl
    });
    assert.equal(evidence.repository, REPOSITORY);
    assert.deepEqual(evidence.workflow, { path: workflowPath, id: workflowId, state: 'active' });
    assert.equal(evidence.run.id, runId);
    assert.equal(evidence.run.headSha, HEAD_SHA);
    assert.equal(evidence.aggregateJob.id, jobId);
    assert.equal(evidence.aggregateJob.name, aggregateName);
    assert.match(evidence.selection.basis, /newest exact-identity run/);
    assert.equal(fetchImpl.calls.length, 5);
  }
});

test('latest successful rerun attempt is inspected through the attempt-specific endpoint', async () => {
  const rerun = runFixture({ run_attempt: 2 });
  const fetchImpl = createFetch({ runPages: [[rerun]], currentRuns: [rerun] });
  const evidence = await verify({ fetchImpl });
  assert.equal(evidence.run.attempt, 2);
  assert.ok(fetchImpl.calls.some(({ url }) => (
    url.includes('/actions/runs/701/attempts/2/jobs')
  )));
  assert.ok(!fetchImpl.calls.some(({ url }) => url.includes('filter=latest')));
});

test('multiple successful exact-identity runs report every ignored older run ID', async () => {
  const older = runFixture({
    id: 700,
    run_number: 40,
    created_at: '2026-08-27T09:00:00Z',
    url: `${API_ORIGIN}${REPOSITORY_API_PATH}/actions/runs/700`,
    html_url: `https://github.com/${REPOSITORY}/actions/runs/700`
  });
  const evidence = await verify({
    fetchImpl: createFetch({ runPages: [[older, runFixture()]] })
  });
  assert.deepEqual(evidence.selection.matchingRunIds, [701, 700]);
  assert.deepEqual(evidence.selection.ignoredOlderRunIds, [700]);
});

test('run and aggregate evidence on later pages is not missed', async () => {
  const wrongRun = runFixture({
    id: 700,
    run_number: 40,
    head_sha: '9'.repeat(40),
    url: `${API_ORIGIN}${REPOSITORY_API_PATH}/actions/runs/700`,
    html_url: `https://github.com/${REPOSITORY}/actions/runs/700`
  });
  const leafJob = jobFixture({ id: 800, name: 'leaf' });
  const fetchImpl = createFetch({
    runPages: [[wrongRun], [runFixture()]],
    jobPages: [[leafJob], [jobFixture()]]
  });
  const evidence = await verify({ fetchImpl });
  assert.equal(evidence.run.id, 701);
  assert.equal(evidence.aggregateJob.id, 801);
  assert.equal(fetchImpl.calls.filter(({ url }) => url.includes('/workflows/501/runs')).length, 2);
  assert.equal(fetchImpl.calls.filter(({ url }) => url.includes('/attempts/1/jobs')).length, 2);
});

test('workflow resolution rejects missing, disabled, and path-mismatched workflows', async () => {
  await rejectsWith(
    verify({ fetchImpl: createFetch({ workflowStatus: 404 }) }),
    'WORKFLOW_NOT_FOUND',
    /not found/
  );
  await rejectsWith(
    verify({ fetchImpl: createFetch({ workflow: workflowFixture({ state: 'disabled_manually' }) }) }),
    'WORKFLOW_IDENTITY_MISMATCH',
    /disabled or does not match/
  );
  await rejectsWith(
    verify({
      fetchImpl: createFetch({
        workflow: workflowFixture({ path: GALLERY_WORKFLOW_PATH })
      })
    }),
    'WORKFLOW_IDENTITY_MISMATCH',
    /exact workflow path/
  );
});

test('unsupported workflow, manual evidence, and wrong aggregate contracts fail before fetch', async () => {
  await rejectsWith(
    verify({ workflowPath: '.github/workflows/deploy_web.yml' }),
    'UNSUPPORTED_WORKFLOW',
    /not a promotion evidence producer/
  );
  await rejectsWith(
    verify({ expectedEvent: 'workflow_dispatch' }),
    'UNSUPPORTED_STAGING_IDENTITY',
    /push runs on dev/
  );
  await rejectsWith(
    verify({ expectedBranch: 'main' }),
    'UNSUPPORTED_STAGING_IDENTITY',
    /push runs on dev/
  );
  await rejectsWith(
    verify({ expectedAggregateName: 'Tests' }),
    'AGGREGATE_CONTRACT_MISMATCH',
    /does not match the workflow contract/
  );
});

test('wrong repository, event, branch, SHA, workflow ID, and absent runs never match', async () => {
  const cases = [
    ['wrong repository', { repository: { full_name: 'attacker/gbdraw' } }],
    ['manual dispatch', { event: 'workflow_dispatch' }],
    ['wrong event', { event: 'pull_request' }],
    ['wrong branch', { head_branch: 'main' }],
    ['wrong SHA', { head_sha: '8'.repeat(40) }],
    ['wrong workflow ID', { workflow_id: 999 }]
  ];
  for (const [, changes] of cases) {
    await rejectsWith(
      verify({ fetchImpl: createFetch({ runPages: [[runFixture(changes)]] }) }),
      'NO_MATCHING_RUN',
      /No workflow run matches/
    );
  }
  await rejectsWith(
    verify({ fetchImpl: createFetch({ runPages: [], runTotalCount: 0 }) }),
    'NO_MATCHING_RUN',
    /No workflow run matches/
  );
});

test('a newer queued or failed run is selected instead of an older success', async () => {
  const older = runFixture({
    id: 700,
    run_number: 40,
    created_at: '2026-08-27T09:00:00Z',
    url: `${API_ORIGIN}${REPOSITORY_API_PATH}/actions/runs/700`,
    html_url: `https://github.com/${REPOSITORY}/actions/runs/700`
  });
  for (const [status, conclusion] of [['queued', null], ['completed', 'failure']]) {
    const newer = runFixture({ status, conclusion });
    await rejectsWith(
      verify({
        fetchImpl: createFetch({
          runPages: [[older, newer]],
          currentRuns: [newer]
        })
      }),
      'RUN_NOT_SUCCESSFUL',
      /Newest exact-identity workflow run/
    );
  }
});

test('cancelled, timed-out, skipped, and conclusion-less current runs fail closed', async () => {
  for (const conclusion of ['cancelled', 'timed_out', 'skipped', null]) {
    const run = runFixture({ conclusion });
    await rejectsWith(
      verify({ fetchImpl: createFetch({ runPages: [[run]], currentRuns: [run] }) }),
      'RUN_NOT_SUCCESSFUL',
      /not completed successfully/
    );
  }
});

test('a newly started current attempt cannot reuse an older successful attempt', async () => {
  const olderAttempt = runFixture({ run_attempt: 1 });
  const currentAttempt = runFixture({ run_attempt: 2, status: 'queued', conclusion: null });
  await rejectsWith(
    verify({
      fetchImpl: createFetch({
        runPages: [[olderAttempt]],
        currentRuns: [currentAttempt]
      })
    }),
    'RUN_NOT_SUCCESSFUL',
    /not completed successfully/
  );
});

test('an attempt change during verification invalidates otherwise successful evidence', async () => {
  const attemptOne = runFixture({ run_attempt: 1 });
  const attemptTwo = runFixture({ run_attempt: 2 });
  await rejectsWith(
    verify({
      fetchImpl: createFetch({
        currentRuns: [attemptOne, attemptTwo]
      })
    }),
    'RUN_ATTEMPT_CHANGED',
    /changed during evidence verification/
  );
});

test('aggregate job must be present exactly once and successful', async () => {
  await rejectsWith(
    verify({ fetchImpl: createFetch({ jobPages: [[jobFixture({ name: 'leaf' })]] }) }),
    'AGGREGATE_JOB_CARDINALITY',
    /exactly one aggregate job/
  );
  await rejectsWith(
    verify({
      fetchImpl: createFetch({
        jobPages: [[jobFixture(), jobFixture({ id: 802 })]]
      })
    }),
    'AGGREGATE_JOB_CARDINALITY',
    /exactly one aggregate job/
  );
  for (const conclusion of ['skipped', 'cancelled', 'failure']) {
    await rejectsWith(
      verify({
        fetchImpl: createFetch({
          jobPages: [[jobFixture({ conclusion })]]
        })
      }),
      'AGGREGATE_JOB_NOT_SUCCESSFUL',
      /not completed successfully/
    );
  }
});

test('malformed run and job metadata fails closed with a bounded reason', async () => {
  const malformedRunFetch = createFetch();
  const originalRunResponse = malformedRunFetch;
  await rejectsWith(
    verify({
      fetchImpl: async (source, options) => {
        const result = await originalRunResponse(source, options);
        if (String(source).includes('/workflows/501/runs')) {
          return response({ total_count: '1', workflow_runs: [] });
        }
        return result;
      }
    }),
    'MALFORMED_API_RESPONSE',
    /response shape is malformed/
  );
  await rejectsWith(
    verify({
      fetchImpl: createFetch({
        jobPages: [[{ ...jobFixture(), name: null }]]
      })
    }),
    'MALFORMED_JOB',
    /job metadata is malformed/
  );
});

test('incomplete, malformed, and looping pagination fails closed', async () => {
  await rejectsWith(
    verify({ fetchImpl: createFetch({ runTotalCount: 2 }) }),
    'INCOMPLETE_PAGINATION',
    /ended before total_count/
  );
  await rejectsWith(
    verify({ fetchImpl: createFetch({ runLinks: ['not-a-link'], runTotalCount: 2 }) }),
    'MALFORMED_PAGINATION',
    /Link header is malformed/
  );
  const skippedPageUrl = new URL(`${API_ORIGIN}${REPOSITORY_API_PATH}/actions/workflows/501/runs`);
  skippedPageUrl.searchParams.set('branch', 'dev');
  skippedPageUrl.searchParams.set('event', 'push');
  skippedPageUrl.searchParams.set('exclude_pull_requests', 'true');
  skippedPageUrl.searchParams.set('head_sha', HEAD_SHA);
  skippedPageUrl.searchParams.set('per_page', '100');
  skippedPageUrl.searchParams.set('page', '3');
  await rejectsWith(
    verify({
      fetchImpl: createFetch({
        runPages: [[runFixture()], [], []],
        runTotalCount: 2,
        runLinks: [`<${skippedPageUrl}>; rel="next"`]
      })
    }),
    'INCOMPLETE_PAGINATION',
    /skipped or reordered a page/
  );
  const initialRunsUrl = new URL(`${API_ORIGIN}${REPOSITORY_API_PATH}/actions/workflows/501/runs`);
  initialRunsUrl.searchParams.set('branch', 'dev');
  initialRunsUrl.searchParams.set('event', 'push');
  initialRunsUrl.searchParams.set('exclude_pull_requests', 'true');
  initialRunsUrl.searchParams.set('head_sha', HEAD_SHA);
  initialRunsUrl.searchParams.set('per_page', '100');
  await rejectsWith(
    verify({
      fetchImpl: createFetch({
        runPages: [[runFixture()], []],
        runTotalCount: 2,
        runLinks: [`<${initialRunsUrl}>; rel="next"`]
      })
    }),
    'PAGINATION_LOOP',
    /repeated a page/
  );
});

const gitResult = (status = 0, stdout = '', stderr = '') => ({ status, stdout, stderr });

const createGitRunner = ({
  shallow = false,
  missing = '',
  mergeStatus = 0,
  mergeTree = HEAD_TREE_SHA
} = {}) => (_repositoryRoot, args) => {
  if (args.join(' ') === 'rev-parse --is-shallow-repository') {
    return gitResult(0, `${shallow}\n`);
  }
  if (args[0] === 'cat-file') {
    return missing && args.at(-1).startsWith(missing)
      ? gitResult(128, '', 'missing object')
      : gitResult();
  }
  if (args[0] === 'rev-parse' && args.includes('--verify')) {
    return gitResult(0, `${HEAD_TREE_SHA}\n`);
  }
  if (args[0] === 'merge-tree') {
    return gitResult(mergeStatus, mergeStatus ? '' : `${mergeTree}\n`, mergeStatus ? 'CONFLICT' : '');
  }
  assert.fail(`Unexpected Git arguments: ${args.join(' ')}`);
};

test('conflict-free synthetic merge tree equal to the head tree succeeds without mutation commands', () => {
  const calls = [];
  const inner = createGitRunner();
  const evidence = verifyPromotionResultTree({
    repositoryRoot: '/fixture/repository',
    baseSha: BASE_SHA,
    headSha: HEAD_SHA,
    runGit: (root, args) => {
      calls.push({ root, args });
      return inner(root, args);
    }
  });
  assert.equal(evidence.conflictFree, true);
  assert.equal(evidence.treeMatchesHead, true);
  assert.equal(evidence.syntheticMergeTreeSha, HEAD_TREE_SHA);
  assert.deepEqual(calls.at(-1).args, [
    'merge-tree', '--write-tree', '--no-messages', BASE_SHA, HEAD_SHA
  ]);
  assert.ok(calls.every(({ args }) => !['checkout', 'merge', 'switch'].includes(args[0])));
});

test('shallow history, missing objects, conflicts, and different trees fail closed', () => {
  const cases = [
    ['INCOMPLETE_GIT_HISTORY', /Complete Git history/, createGitRunner({ shallow: true })],
    ['MISSING_GIT_OBJECT', /Complete Git object/, createGitRunner({ missing: BASE_SHA })],
    ['MERGE_CONFLICT', /not conflict-free/, createGitRunner({ mergeStatus: 1 })],
    ['RESULT_TREE_MISMATCH', /differs from the dev head tree/, createGitRunner({ mergeTree: OTHER_TREE_SHA })]
  ];
  for (const [code, reason, runGit] of cases) {
    assert.throws(
      () => verifyPromotionResultTree({
        repositoryRoot: '/fixture/repository',
        baseSha: BASE_SHA,
        headSha: HEAD_SHA,
        runGit
      }),
      (error) => error instanceof PromotionReadinessError
        && error.code === code
        && reason.test(error.message)
    );
  }
});

test('invalid SHA arguments are rejected before Git executes', () => {
  let executed = false;
  assert.throws(
    () => verifyPromotionResultTree({
      repositoryRoot: '/fixture/repository',
      baseSha: 'HEAD; touch /tmp/injected',
      headSha: HEAD_SHA,
      runGit: () => {
        executed = true;
        return gitResult();
      }
    }),
    (error) => error instanceof PromotionReadinessError
      && error.code === 'INVALID_OBJECT_ID'
  );
  assert.equal(executed, false);
});

test('CLI failure returns nonzero with expected identity and bounded diagnostics', async () => {
  let stdout = '';
  let stderr = '';
  const status = await runPromotionReadinessCli({
    argv: ['workflow-evidence', '--repository', REPOSITORY],
    env: { GITHUB_TOKEN: 'fixture-token' },
    stdout: { write: (value) => { stdout += value; } },
    stderr: { write: (value) => { stderr += value; } }
  });
  assert.equal(status, 1);
  assert.equal(stdout, '');
  assert.match(stderr, /MISSING_ARGUMENT/);
  assert.match(stderr, /Expected:/);
  assert.match(stderr, /Observed:/);
  assert.ok(!stderr.includes('fixture-token'));
});
