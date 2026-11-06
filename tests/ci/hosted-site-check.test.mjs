import assert from 'node:assert/strict';
import { execFileSync } from 'node:child_process';
import { readFileSync } from 'node:fs';
import { fileURLToPath } from 'node:url';
import test from 'node:test';

import {
  CHECK_APP, CHECK_NAME, latestBuildRun, parseArgs, parseBuildStamp, parseVersionId, renderSummary,
  responseProblem, servedState
} from '../../tools/check-hosted-site.mjs';

const SHA = '033a3aa899a6966925df6f413bd29e4a2d342777';
const LATER = 'b22e4b6cc0000000000000000000000000000000';
const workflow = readFileSync(new URL('../../.github/workflows/deploy_web.yml', import.meta.url), 'utf8');
// The version label as its owner, tools/stamp_web_build.py, renders it.
const renderedLabel = (commit) => execFileSync('python3', ['-c', [
  'import runpy, sys',
  "render = runpy.run_path(sys.argv[1])['render_build_label']",
  "print(render(project_version='0.14.5', commit_sha=sys.argv[2]))"
].join('\n'), fileURLToPath(new URL('../../tools/stamp_web_build.py', import.meta.url)), commit], { encoding: 'utf8' });
const job = (id) => workflow.match(new RegExp(`\\n  ${id}:\\n[\\s\\S]*?(?=\\n  [a-z0-9-]+:\\n|$)`))?.[0] || '';

test('the build stamp names the full commit that tools/stamp_web_build.py writes', () => {
  assert.equal(parseBuildStamp(`<footer>${renderedLabel(SHA)}</footer>`), SHA);
  assert.equal(parseBuildStamp(renderedLabel('unknown')), null);
  assert.equal(parseBuildStamp('<html></html>'), null);
});

test('the Version ID comes from its own line of the Workers Builds summary', () => {
  const summary = '\nBuild ID: [2b68](https://dash.cloudflare.com/x)\nScript: [gbdraw](https://x)\n'
    + 'Version ID: bbe3fa12-e929-4348-84a5-4b7c6ebef4ef\n';
  assert.equal(parseVersionId(summary), 'bbe3fa12-e929-4348-84a5-4b7c6ebef4ef');
  assert.equal(parseVersionId('Build ID: x\n'), null);
  assert.equal(parseVersionId(null), null);
});

test('the newest Workers Builds check run of the Cloudflare app is used', () => {
  const run = (id, overrides = {}) => ({ id, name: CHECK_NAME, app: { slug: CHECK_APP }, ...overrides });
  assert.equal(latestBuildRun([run(1), run(3), run(2)]).id, 3);
  assert.equal(latestBuildRun([run(1), run(5, { name: 'Other' }), run(6, { app: { slug: 'other' } })]).id, 1);
  assert.equal(latestBuildRun([]), null);
});

test('a served later main commit is superseded; an older or missing one is still pending', () => {
  assert.equal(servedState({ expected: SHA, served: SHA, compareStatus: 'identical' }), 'match');
  assert.equal(servedState({ expected: SHA, served: LATER, compareStatus: 'ahead' }), 'superseded');
  assert.equal(servedState({ expected: SHA, served: LATER, compareStatus: 'behind' }), 'pending');
  assert.equal(servedState({ expected: SHA, served: LATER, compareStatus: 'diverged' }), 'pending');
  assert.equal(servedState({ expected: SHA, served: null, compareStatus: undefined }), 'pending');
});

test('a client error such as a Cloudflare challenge fails at once; a server error is retried', () => {
  assert.deepEqual(responseProblem(403, 'challenge'), { retry: false, text: 'HTTP 403 (cf-mitigated: challenge)' });
  assert.deepEqual(responseProblem(404, null), { retry: false, text: 'HTTP 404' });
  assert.deepEqual(responseProblem(503, null), { retry: true, text: 'HTTP 503' });
});

test('arguments take a full SHA and second-valued timeouts', () => {
  const options = parseArgs(['--sha', SHA, '--serve-timeout', '60', '--interval', '5']);
  assert.equal(options.sha, SHA);
  assert.equal(options.url, 'https://gbdraw.app/');
  assert.equal(options.serveTimeout, 60);
  assert.equal(options.interval, 5);
  assert.throws(() => parseArgs(['--sha', SHA.slice(0, 7)]), /40-character/);
  assert.throws(() => parseArgs(['--sha', SHA, '--interval', 'soon']), /seconds/);
  assert.throws(() => parseArgs(['--sha', SHA, '--branch', 'main']), /unknown/);
});

test('the summary names the commit, Version ID, served state, startup, and any failure', () => {
  const passed = renderSummary({
    sha: SHA, url: 'https://gbdraw.app/', versionId: 'v-1', served: 'superseded', servedSha: LATER, startupSeconds: 12.34
  });
  assert.match(passed, /^### Hosted site check: passed/);
  assert.match(passed, /Version ID `v-1`/);
  assert.match(passed, new RegExp(`serves \`${LATER}\`, a later main commit \\(superseded\\)`));
  assert.match(passed, /started in 12\.3 s/);
  const failed = renderSummary({ sha: SHA, url: 'https://gbdraw.app/', failure: 'still serves x' });
  assert.match(failed, /^### Hosted site check: failed/);
  assert.match(failed, /- Failure: still serves x/);
  assert.doesNotMatch(failed, /Version ID/);
});

test('the hosted job runs on main pushes, daily, and on a dispatch that names site_sha, without gbdraw', () => {
  const hosted = job('hosted-site');
  const gallery = job('gallery-publication-browser');
  assert.match(workflow, /\n  workflow_dispatch:\n    inputs:\n      site_sha:\n/);
  assert.match(workflow, /\n  schedule:\n    - cron: "17 21 \* \* \*"\n/);
  assert.match(hosted, /\n    if: github\.event_name == 'push' \|\| github\.event_name == 'schedule' \|\| inputs\.site_sha != ''\n/);
  // The schedule has no inputs, so the Gallery job must name the dispatch event.
  assert.match(gallery, /\n    if: github\.event_name == 'push' \|\| \(github\.event_name == 'workflow_dispatch' && inputs\.site_sha == ''\)\n/);
  assert.match(hosted, /\n    permissions:\n      contents: read\n      checks: read\n/);
  assert.match(hosted, /SITE_SHA: \$\{\{ inputs\.site_sha \|\| github\.sha \}\}/);
  assert.match(hosted, /node tools\/check-hosted-site\.mjs --sha "\$SITE_SHA" --url https:\/\/gbdraw\.genome-kawato\.workers\.dev\//);
  assert.doesNotMatch(hosted, /pip install|setup-python|managed-losat/);
});
