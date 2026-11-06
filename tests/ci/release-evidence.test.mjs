import assert from 'node:assert/strict';
import { execFileSync } from 'node:child_process';
import { mkdirSync, mkdtempSync, readFileSync, rmSync, writeFileSync } from 'node:fs';
import { tmpdir } from 'node:os';
import { dirname, join } from 'node:path';
import test from 'node:test';

import { RELEASE_GATE_JOB, runReleaseEvidenceCli } from '../../tools/check-release-evidence.mjs';

const REPOSITORY = 'satoshikawato/gbdraw';

// History: root A on main and dev; dev adds E (docs), R (runtime), H (docs); main merges H as
// the tagged promotion commit T; dev then moves to D (runtime).
const createHistory = () => {
  const root = mkdtempSync(join(tmpdir(), 'release-evidence-'));
  const git = (...args) => execFileSync('git', ['-c', 'user.name=t', '-c', 'user.email=t@example.com', ...args], {
    cwd: root, encoding: 'utf8', stdio: ['ignore', 'pipe', 'pipe']
  }).trim();
  const commit = (path, text, message) => {
    mkdirSync(dirname(join(root, path)), { recursive: true });
    writeFileSync(join(root, path), text);
    git('add', path);
    git('commit', '-q', '-m', message);
    return git('rev-parse', 'HEAD');
  };
  git('init', '-q', '-b', 'main');
  const a = commit('gbdraw/draw.py', 'x = 1\n', 'A');
  git('switch', '-q', '-c', 'dev');
  const e = commit('docs/notes.md', 'one\n', 'E');
  const r = commit('gbdraw/draw.py', 'x = 2\n', 'R');
  const h = commit('docs/notes.md', 'two\n', 'H');
  git('switch', '-q', 'main');
  git('merge', '-q', '--no-ff', '-m', 'Promote', h);
  const t = git('rev-parse', 'HEAD');
  git('switch', '-q', 'dev');
  const d = commit('gbdraw/draw.py', 'x = 3\n', 'D');
  git('update-ref', 'refs/remotes/origin/dev', d);
  return { root, a, e, r, h, t, d, cleanup: () => rmSync(root, { recursive: true, force: true }) };
};

let minute = 0;
const dispatch = (id, headSha, overrides = {}) => ({
  databaseId: id,
  attempt: 1,
  headSha,
  url: `https://github.com/${REPOSITORY}/actions/runs/${id}`,
  createdAt: new Date(Date.UTC(2026, 10, 1, 0, minute++)).toISOString(),
  event: 'workflow_dispatch',
  headBranch: 'dev',
  status: 'completed',
  workflowName: 'Tests',
  ...overrides
});
const gate = (conclusion) => [{ name: RELEASE_GATE_JOB, status: 'completed', conclusion, url: 'https://example.invalid/job' }];

// A fake `gh` that answers `run list` and `run view <id>`; `failList` makes `run list` fail.
// `listOutput` replaces the `run list` output; `views` maps a run ID to its `run view` payload.
const fakeGh = ({ runs, jobs = {}, views = {}, failList = false, listOutput = null }) => {
  const calls = [];
  const run = (command, args) => {
    calls.push(args.join(' '));
    assert.equal(command, 'gh');
    if (args[0] === 'run' && args[1] === 'list') {
      if (failList) return { status: 1, stdout: '', stderr: 'HTTP 401: Bad credentials' };
      return { status: 0, stdout: listOutput ?? JSON.stringify(runs), stderr: '' };
    }
    if (args[0] === 'run' && args[1] === 'view') {
      const view = Object.hasOwn(views, args[2]) ? views[args[2]] : { jobs: jobs[args[2]] ?? [] };
      return { status: 0, stdout: JSON.stringify(view), stderr: '' };
    }
    return assert.fail(`unexpected gh call: ${args.join(' ')}`);
  };
  run.calls = calls;
  return run;
};

const check = (history, gh) => {
  let out = '';
  let err = '';
  const summaries = [];
  const status = runReleaseEvidenceCli({
    argv: ['--tag-commit', history.t],
    env: { GITHUB_REPOSITORY: REPOSITORY, GITHUB_STEP_SUMMARY: 'summary.md' },
    cwd: history.root,
    stdout: { write: (text) => { out += text; } },
    stderr: { write: (text) => { err += text; } },
    run: gh,
    appendFileImpl: (path, text) => summaries.push([path, text])
  });
  return { status, out, err, summaries };
};

test('a release-tier run on the promotion commit carries to the tag, whatever the run result', (t) => {
  const history = createHistory();
  t.after(history.cleanup);
  // A report-only promotion audit job failed the run; its `Release / gate` passed.
  const gh = fakeGh({
    runs: [dispatch(11, history.h, { attempt: 2, conclusion: 'failure' })],
    jobs: { 11: gate('success') }
  });
  const result = check(history, gh);
  assert.equal(result.status, 0, result.err);
  assert.deepEqual(gh.calls, [
    `run list --repo ${REPOSITORY} --workflow test.yml --branch dev --event workflow_dispatch --status completed `
      + '--limit 100 --json databaseId,attempt,headSha,url,createdAt,event,headBranch,status,workflowName',
    `run view 11 --repo ${REPOSITORY} --attempt 2 --json jobs`
  ]);
  assert.match(result.out, /^### Release evidence: carries to the tag/);
  assert.match(result.out, new RegExp(`run 11\\]\\(https://github.com/${REPOSITORY}/actions/runs/11\\)`));
  assert.match(result.out, /`releaseEvidenceCarries` is true \(no changed path\)/);
  assert.deepEqual(result.summaries, [['summary.md', result.out]]);
});

test('an older release-tier run counts when a newer dispatch is dev tier and only documentation changed since', (t) => {
  const history = createHistory();
  t.after(history.cleanup);
  const gh = fakeGh({
    runs: [dispatch(21, history.r), dispatch(22, history.h)],
    jobs: { 21: gate('success'), 22: gate('skipped') }
  });
  const result = check(history, gh);
  assert.equal(result.status, 0, result.err);
  assert.match(result.out, /run 21\]/);
  assert.match(result.out, /\(1 changed path\)/);
  assert.match(result.out, /run 22\].*: `Release \/ gate` was skipped \(a dev-tier dispatch\)/);
  assert.deepEqual(gh.calls.filter((call) => call.startsWith('run view')).map((call) => call.split(' ')[2]), ['22', '21']);
});

test('a release without carrying evidence stops before publication', (t) => {
  const history = createHistory();
  t.after(history.cleanup);
  const cases = [
    // Evidence on E does not carry over the runtime change R.
    [{ runs: [dispatch(31, history.e)], jobs: { 31: gate('success') } }, /does not carry over the 2 changed paths/],
    // A run on dev after the tag tests another tree.
    [{ runs: [dispatch(32, history.d)], jobs: { 32: gate('success') } }, /not an ancestor of the tag commit/],
    // The release tier ran on H but its gate failed.
    [{ runs: [dispatch(33, history.h)], jobs: { 33: gate('failure') } }, /`Release \/ gate` ended failure/],
    // The newest gate result decides: a newer failure on R stops the check although H passed earlier.
    [{ runs: [dispatch(34, history.h), dispatch(35, history.r)], jobs: { 34: gate('success'), 35: gate('timed_out') } },
      /run 35\].*ended timed_out; the newest release-tier result decides/],
    // A cancelled gate, two gates, and a run without a job list decide nothing.
    [{
      runs: [dispatch(36, history.h), dispatch(37, history.h), dispatch(38, history.h)],
      jobs: { 36: gate('cancelled'), 37: [...gate('success'), ...gate('success')] },
      views: { 38: { jobs: [] } }
    }, /run 38\].*no single[\s\S]*run 37\].*no single[\s\S]*run 36\].*ended cancelled/],
    // No dispatch run at all.
    [{ runs: [] }, /read: 0/]
  ];
  for (const [fixture, reason] of cases) {
    const result = check(history, fakeGh(fixture));
    assert.equal(result.status, 1, reason.source);
    assert.match(result.out, /^### Release evidence: missing/);
    assert.match(result.out, reason);
    assert.match(result.out, new RegExp(`second parent\\): \`${history.h}\``));
    assert.match(result.out, /dev now: `[0-9a-f]{12}`, not the promotion's dev commit/);
    assert.match(result.out, /"Release evidence check"/);
  }
});

test('the check fails closed when gh fails or lists another kind of run', (t) => {
  const history = createHistory();
  t.after(history.cleanup);
  const failed = check(history, fakeGh({ runs: [], failList: true }));
  assert.equal(failed.status, 1);
  assert.match(failed.err, /\[GH_FAILED\].*Bad credentials/);
  const pushRun = check(history, fakeGh({ runs: [dispatch(41, history.h, { event: 'push' })] }));
  assert.equal(pushRun.status, 1);
  assert.match(pushRun.err, /\[UNEXPECTED_RUN\]/);
  const malformed = check(history, fakeGh({ runs: [], listOutput: 'not json' }));
  assert.match(malformed.err, /\[GH_MALFORMED\]/);
  const noJobs = check(history, fakeGh({ runs: [dispatch(42, history.h)], views: { 42: { jobs: null } } }));
  assert.match(noJobs.err, /\[GH_MALFORMED\].*42/);
  assert.deepEqual([malformed.status, noJobs.status], [1, 1]);
  assert.equal(failed.out + pushRun.out + malformed.out + noJobs.out, '');
});

test('the gate job the check reads is the release-gate job of the Tests workflow', () => {
  const workflow = readFileSync(new URL('../../.github/workflows/test.yml', import.meta.url), 'utf8');
  assert.match(workflow, /^name: Tests$/m);
  assert.match(workflow, new RegExp(`\\n  release-gate:\\n    name: ${RELEASE_GATE_JOB.replace('/', '\\/')}\\n`));
});
