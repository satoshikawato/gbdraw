#!/usr/bin/env node
// Release evidence check (docs/internal/WEB_PERIODIC_AUDIT.md, "Release evidence check").
//
//   node tools/check-release-evidence.mjs --tag-commit <commit> [--repository OWNER/REPO]
//
// `pypi-publish` in .github/workflows/release.yml waits for this check. It passes when a
// release-tier `Tests` run carries its evidence to the tagged commit T: a completed
// `workflow_dispatch` run on dev, among the newest RUN_LIMIT, whose `Release / gate` job
// succeeded, on a commit E for which `classifyRange` of tools/ci-impact.mjs (the owner of the
// carry rule, `classify --base E --head T`) reports `releaseEvidenceCarries`, and no newer run on
// such a commit failed that job. The run's own conclusion does not count: the report-only
// promotion audit jobs can fail it.
//
// Reads runs with the `gh` CLI (GH_TOKEN or GITHUB_TOKEN) and history from the Git clone in
// the working directory. Writes only stdout, stderr, and $GITHUB_STEP_SUMMARY.
import { spawnSync } from 'node:child_process';
import { appendFileSync } from 'node:fs';
import { resolve } from 'node:path';
import { pathToFileURL } from 'node:url';
import { classifyRange, runGit } from './ci-impact.mjs';

// The name of the `release-gate` job in .github/workflows/test.yml.
export const RELEASE_GATE_JOB = 'Release / gate';
export const RUN_LIMIT = 100;
const SUMMARY_RUN_LIMIT = 20;
const SHA = /^[0-9a-f]{40}$/;
const REPOSITORY = /^[A-Za-z0-9][A-Za-z0-9._-]*\/[A-Za-z0-9][A-Za-z0-9._-]*$/;
const RUN_FIELDS = 'databaseId,attempt,headSha,url,createdAt,event,headBranch,status,workflowName';

class ReleaseEvidenceError extends Error {
  constructor(code, message) {
    super(message);
    this.code = code;
  }
}
const fail = (code, message) => {
  throw new ReleaseEvidenceError(code, message);
};
const bounded = (text) => String(text ?? '').replace(/\s+/g, ' ').trim().slice(0, 300);
const short = (sha) => sha.slice(0, 12);
const paths = (count) => `${count} changed path${count === 1 ? '' : 's'}`;

export const runCommand = (command, args, cwd) => spawnSync(command, args, {
  cwd,
  encoding: 'utf8',
  maxBuffer: 64 * 1024 * 1024,
  stdio: ['ignore', 'pipe', 'pipe']
});

const readGh = (run, args, what) => {
  const result = run('gh', args);
  if (result.error || result.status !== 0) {
    fail('GH_FAILED', `gh could not read ${what}: ${bounded(result.stderr || result.error?.message)}`);
  }
  try {
    return JSON.parse(result.stdout);
  } catch {
    fail('GH_MALFORMED', `gh returned malformed JSON for ${what}.`);
  }
};

// The newest completed `Tests` dispatch runs on dev. A run of another identity fails the check.
export const listDispatchRuns = ({ repository, run }) => {
  const runs = readGh(run, [
    'run', 'list', '--repo', repository, '--workflow', 'test.yml', '--branch', 'dev',
    '--event', 'workflow_dispatch', '--status', 'completed', '--limit', String(RUN_LIMIT),
    '--json', RUN_FIELDS
  ], 'the Tests dispatch runs on dev');
  if (!Array.isArray(runs)) fail('GH_MALFORMED', 'gh run list did not return a list.');
  for (const entry of runs) {
    if (entry?.workflowName !== 'Tests' || entry.event !== 'workflow_dispatch' || entry.headBranch !== 'dev'
        || entry.status !== 'completed' || !SHA.test(entry.headSha) || !Number.isSafeInteger(entry.databaseId)
        || !Number.isSafeInteger(entry.attempt) || entry.attempt < 1 || !Number.isFinite(Date.parse(entry.createdAt))
        || typeof entry.url !== 'string') {
      fail('UNEXPECTED_RUN', `gh run list returned a run that is not a completed Tests dispatch on dev: ${bounded(JSON.stringify(entry))}`);
    }
  }
  return runs;
};

// The `Release / gate` job of the listed attempt, or null when the attempt has none or several.
export const readReleaseGate = ({ repository, run }) => (entry) => {
  const view = readGh(run, [
    'run', 'view', String(entry.databaseId), '--repo', repository, '--attempt', String(entry.attempt),
    '--json', 'jobs'
  ], `the jobs of run ${entry.databaseId}`);
  if (!Array.isArray(view?.jobs)) fail('GH_MALFORMED', `gh run view ${entry.databaseId} returned no job list.`);
  const gates = view.jobs.filter((job) => job?.name === RELEASE_GATE_JOB);
  return gates.length === 1 ? gates[0] : null;
};

// Whether the evidence of a run on `base` carries to the tag commit, by the one carry rule.
export const carryVerdict = ({ repositoryRoot, tagCommit, runGitImpl = runGit }) => (base) => {
  // classify reports a non-ancestor as not carrying only after diffing it; skip that diff.
  if (runGitImpl(repositoryRoot, ['merge-base', '--is-ancestor', base, tagCommit])?.status === 1) {
    return { carries: false, reason: 'its commit is not an ancestor of the tag commit' };
  }
  let result;
  try {
    result = classifyRange({ repositoryRoot, base, head: tagCommit, runGitImpl });
  } catch (error) {
    if (error?.code === 'UNKNOWN_REVISION') return { carries: false, reason: 'its commit is not in this clone' };
    throw error;
  }
  if (!result.ancestor) return { carries: false, reason: 'its commit is not an ancestor of the tag commit' };
  if (!result.verdicts.releaseEvidenceCarries) {
    return {
      carries: false,
      reason: `release evidence does not carry over the ${paths(result.paths.length)} since its commit`
    };
  }
  return { carries: true, changedPaths: result.paths.length };
};

// Gate results that stop the check. Any other gate (skipped by a dev-tier dispatch, cancelled,
// missing, or not finished) has no result and is passed over.
const FAILED_GATE = new Set(['failure', 'timed_out', 'startup_failure']);

const gateReason = (gate) => {
  if (!gate) return `it has no single \`${RELEASE_GATE_JOB}\` job`;
  if (gate.conclusion === 'skipped') return `\`${RELEASE_GATE_JOB}\` was skipped (a dev-tier dispatch)`;
  const result = bounded(gate.conclusion || gate.status) || 'without a result';
  return FAILED_GATE.has(gate.conclusion)
    ? `\`${RELEASE_GATE_JOB}\` ended ${result}; the newest release-tier result decides, so older runs do not count`
    : `\`${RELEASE_GATE_JOB}\` ended ${result}`;
};

// Newest run first. Among runs whose commit carries to the tag, the newest whose
// `Release / gate` has a result decides: success is the evidence, a failure stops the check.
export const selectReleaseEvidence = ({ runs, carries, releaseGate }) => {
  const examined = [];
  const ordered = [...runs].sort((left, right) => (
    Date.parse(right.createdAt) - Date.parse(left.createdAt) || right.databaseId - left.databaseId
  ));
  for (const entry of ordered) {
    const verdict = carries(entry.headSha);
    if (!verdict.carries) {
      examined.push({ run: entry, reason: verdict.reason });
      continue;
    }
    const gate = releaseGate(entry);
    if (gate?.status === 'completed' && gate.conclusion === 'success') {
      return { evidence: { run: entry, gate, verdict }, examined };
    }
    examined.push({ run: entry, reason: gateReason(gate) });
    if (gate?.status === 'completed' && FAILED_GATE.has(gate.conclusion)) break;
  }
  return { evidence: null, examined };
};

const runLink = (entry) => `[run ${entry.databaseId}](${entry.url}) on \`${short(entry.headSha)}\``;

export const renderSummary = ({ tagCommit, promotionHead, devHead, runCount, evidence, examined }) => {
  const passedOver = examined.slice(0, SUMMARY_RUN_LIMIT).map(({ run, reason }) => `  - ${runLink(run)}: ${reason}`);
  const more = examined.length > SUMMARY_RUN_LIMIT ? [`  - ${examined.length - SUMMARY_RUN_LIMIT} more`] : [];
  if (evidence) {
    const { run, gate, verdict } = evidence;
    return [
      '### Release evidence: carries to the tag',
      '',
      `- Tag commit: \`${tagCommit}\``,
      `- Release-tier \`Tests\` ${runLink(run)}, attempt ${run.attempt}: \`${RELEASE_GATE_JOB}\` succeeded${gate.url ? ` ([job](${gate.url}))` : ''}`,
      `- \`node tools/ci-impact.mjs classify --base ${run.headSha} --head ${tagCommit}\`: \`releaseEvidenceCarries\` is true`
        + ` (${verdict.changedPaths ? paths(verdict.changedPaths) : 'no changed path'})`,
      ...(examined.length ? ['- Newer dispatch runs passed over:', ...passedOver, ...more] : []),
      ''
    ].join('\n');
  }
  const devState = !devHead ? 'not in this clone'
    : devHead === promotionHead ? `\`${short(devHead)}\`, still the promotion's dev commit`
      : `\`${short(devHead)}\`, not the promotion's dev commit`;
  return [
    '### Release evidence: missing',
    '',
    'PyPI publication waits: no release-tier `Tests` run carries its evidence to the tag commit.',
    '',
    `- Tag commit: \`${tagCommit}\``,
    `- Promotion's dev commit (the tag commit's second parent): ${promotionHead ? `\`${promotionHead}\`` : 'none'}`,
    `- dev now: ${devState}`,
    `- Completed \`Tests\` dispatch runs on dev read: ${runCount} (newest ${RUN_LIMIT} at most)`,
    ...(examined.length ? ['- Why each run does not count, newest first:', ...passedOver, ...more] : []),
    '- What to do: `docs/internal/WEB_PERIODIC_AUDIT.md`, "Release evidence check".',
    ''
  ].join('\n');
};

const parseArgs = (argv, env) => {
  const options = { repository: env.GITHUB_REPOSITORY || '' };
  for (let index = 0; index < argv.length; index += 2) {
    const name = { '--tag-commit': 'tagCommit', '--repository': 'repository' }[argv[index]];
    if (!name || argv[index + 1] === undefined) fail('INVALID_ARGUMENTS', `Unknown or incomplete argument: ${bounded(argv[index])}`);
    options[name] = argv[index + 1];
  }
  if (!options.tagCommit) fail('INVALID_ARGUMENTS', 'Usage: check-release-evidence.mjs --tag-commit <commit> [--repository OWNER/REPO]');
  if (!REPOSITORY.test(options.repository)) fail('INVALID_ARGUMENTS', 'Set --repository or GITHUB_REPOSITORY to OWNER/REPO.');
  return options;
};

export const runReleaseEvidenceCli = ({
  argv = process.argv.slice(2),
  env = process.env,
  cwd = process.cwd(),
  stdout = process.stdout,
  stderr = process.stderr,
  run = runCommand,
  runGitImpl = runGit,
  appendFileImpl = appendFileSync
} = {}) => {
  try {
    const { tagCommit: tagRevision, repository } = parseArgs(argv, env);
    const repositoryRoot = resolve(cwd);
    const revision = (name) => {
      const result = runGitImpl(repositoryRoot, ['rev-parse', '--verify', '--quiet', '--end-of-options', `${name}^{commit}`]);
      const sha = result?.status === 0 ? Buffer.from(result.stdout).toString('utf8').trim() : '';
      return SHA.test(sha) ? sha : null;
    };
    const tagCommit = revision(tagRevision);
    if (!tagCommit) fail('UNKNOWN_TAG_COMMIT', `The tag commit ${bounded(tagRevision)} is not in this clone.`);
    const runs = listDispatchRuns({ repository, run: (command, args) => run(command, args, repositoryRoot) });
    const { evidence, examined } = selectReleaseEvidence({
      runs,
      carries: carryVerdict({ repositoryRoot, tagCommit, runGitImpl }),
      releaseGate: readReleaseGate({ repository, run: (command, args) => run(command, args, repositoryRoot) })
    });
    const summary = renderSummary({
      tagCommit,
      promotionHead: revision(`${tagCommit}^2`),
      devHead: revision('refs/remotes/origin/dev'),
      runCount: runs.length,
      evidence,
      examined
    });
    stdout.write(summary);
    if (env.GITHUB_STEP_SUMMARY) appendFileImpl(env.GITHUB_STEP_SUMMARY, summary);
    return evidence ? 0 : 1;
  } catch (error) {
    const code = typeof error?.code === 'string' ? error.code : 'UNEXPECTED_ERROR';
    stderr.write(`Release evidence check failed [${code}]: ${bounded(error?.message)}\n`);
    return 1;
  }
};

if (process.argv[1] && import.meta.url === pathToFileURL(resolve(process.argv[1])).href) {
  process.exitCode = runReleaseEvidenceCli();
}
