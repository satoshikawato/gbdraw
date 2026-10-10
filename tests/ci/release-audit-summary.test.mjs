import assert from 'node:assert/strict';
import { mkdirSync, mkdtempSync, readFileSync, writeFileSync } from 'node:fs';
import { tmpdir } from 'node:os';
import { join } from 'node:path';
import test from 'node:test';
import { readArtifact, renderSummary } from '../../tools/audit/release-audit-summary.mjs';

const workflow = readFileSync(new URL('../../.github/workflows/test.yml', import.meta.url), 'utf8');
const job = (id) => workflow.match(new RegExp(`\\n  ${id}:\\n[\\s\\S]*?(?=\\n  [a-z0-9-]+:\\n|$)`))?.[0] || '';
const AUDIT_JOBS = ['promotion-audit-recipes', 'promotion-audit-sweeps', 'promotion-audit-journeys', 'promotion-audit-random-walk'];

const artifact = (files) => {
  const dir = mkdtempSync(join(tmpdir(), 'audit-artifact-'));
  for (const [path, content] of Object.entries(files)) {
    mkdirSync(join(dir, path, '..'), { recursive: true });
    writeFileSync(join(dir, path), typeof content === 'string' ? content : JSON.stringify(content));
  }
  return dir;
};

test('an artifact counts failed tests, failed recipes, and parity replay differences', () => {
  const dir = artifact({
    'report.json': { stats: { expected: 3, unexpected: 2 } },
    'recipes/recipes-status.tsv': 'cli\tT-CLI-01\tpass\ncli\tH-CLI-02\tfail\npython\tT-PY-01\tpass\n'
      + 'cli\tH-CLI-06\tskipped-run-locally\npython\tT-PY-05\tskipped-run-locally\n',
    'parity/linear-MJNV/replay-report.json': [{ outcome: 'DIFF' }, { outcome: 'MATCH' }, { outcome: 'NO_EFFECT' }],
    'losat.txt': 'unavailable: no bundled LOSAT\n',
    'seed.txt': '1234\n'
  });
  const read = readArtifact(dir);
  assert.equal(read.findings, 2 + 1 + 1);
  assert.ok(read.notes.includes('2 failed test(s)'));
  assert.ok(read.notes.includes('2 of 3 recipes pass; failed: H-CLI-02; run locally (no LOSAT runtime): H-CLI-06, T-PY-05'));
  assert.ok(read.notes.includes('linear-MJNV: 1 DIFF, 1 NO_EFFECT of 3'));
  assert.ok(read.notes.some((note) => note.startsWith('LOSAT runtime unavailable')));
  assert.ok(read.notes.includes('seed 1234'));
  assert.deepEqual(readArtifact(join(dir, 'missing')), { findings: null, notes: ['no artifact'], coverage: null });
});

test('the summary lists each audit job with its artifact link and the uncovered capabilities', () => {
  const root = mkdtempSync(join(tmpdir(), 'audit-artifacts-'));
  const journeys = join(root, 'promotion-audit-journeys-J1');
  mkdirSync(join(journeys, 'journeys'), { recursive: true });
  writeFileSync(join(journeys, 'report.json'), JSON.stringify({ stats: { unexpected: 0 } }));
  writeFileSync(join(journeys, 'journeys', 'coverage.json'), JSON.stringify({
    changed: ['web-runtime', 'packaging', 'full'], uncovered: ['packaging', 'full'],
    uncoveredPaths: { packaging: ['pyproject.toml'], full: ['odd\\|name'] }
  }));
  const markdown = renderSummary({
    jobs: [
      { name: 'Core (Python 3.11)', conclusion: 'success' },
      { name: 'Promotion audit / summary', status: 'in_progress', conclusion: null },
      { name: 'Promotion audit / journeys-J1', conclusion: 'success' },
      { name: 'Promotion audit / recipes', conclusion: 'failure' }
    ],
    artifacts: [{ name: 'promotion-audit-journeys-J1', id: 7 }],
    artifactsDir: root,
    runUrl: 'https://github.com/o/r/actions/runs/9'
  });
  assert.match(markdown, /\| journeys-J1 \| success \| \[promotion-audit-journeys-J1\]\(https:\/\/github\.com\/o\/r\/actions\/runs\/9\/artifacts\/7\) \| 0 \| {2}\|/);
  assert.match(markdown, /\| recipes \| failure \| none \| unknown \| no artifact \|/);
  assert.doesNotMatch(markdown, /Core|\| summary \|/);
  assert.match(markdown, /Uncovered by the user journeys \(hand look with a time limit, or an Owner waiver\): packaging/);
  assert.match(markdown, /\| packaging \| `pyproject\.toml` \|/);
  assert.ok(markdown.includes('| full | `odd\\\\\\|name` |'), 'backslash and pipe are both escaped');
});

test('promotion audit jobs run only on the release-tier dispatch and stay outside every gate', () => {
  const dispatch = workflow.split('\n  workflow_dispatch:\n')[1].split('\nconcurrency:')[0];
  assert.deepEqual([...dispatch.matchAll(/\n {6}([a-z_]+):\n/g)].map(([, name]) => name), ['tier']);
  assert.match(dispatch, /options: \[dev, release\]/);
  for (const id of [...AUDIT_JOBS, 'promotion-audit-summary']) {
    const source = job(id);
    assert.ok(source, id);
    assert.match(source, /github\.event_name == 'workflow_dispatch' &&\n\s+github\.ref == 'refs\/heads\/dev' && inputs\.tier == 'release'/, id);
    assert.doesNotMatch(source, /needs: ci-impact|requiredJobs|continue-on-error/, id);
    for (const gate of ['pr-gate', 'dev-staging-gate', 'release-gate']) {
      assert.ok(!job(gate).includes(`- ${id}\n`), `${gate} must not need ${id}`);
    }
  }
  for (const id of AUDIT_JOBS) {
    assert.ok(job('promotion-audit-summary').includes(`      - ${id}\n`), `the summary needs ${id}`);
    assert.match(job(id), /name: promotion-audit-/, `${id} uploads a promotion-audit- artifact`);
  }
  assert.match(job('promotion-audit-summary'), /always\(\) && github\.event_name == 'workflow_dispatch'/);
  assert.match(job('promotion-audit-journeys'), /journey: \[J1, J2, J3, J4, J5, J6, J7\]/);
  assert.match(job('promotion-audit-random-walk'), /- name: Print the random walk seed\n[\s\S]*GBDRAW_RANDOM_WALK_SEED=\$GITHUB_RUN_ID/);
  assert.ok(job('promotion-audit-random-walk').indexOf('Print the random walk seed') < job('promotion-audit-random-walk').indexOf('- name: Checkout'));
  // A failed LOSAT install leaves the job running; the probe only resolves the runtime.
  assert.match(job('promotion-audit-recipes'), /run: gbdraw setup-losat \|\| echo "::warning::/);
  assert.doesNotMatch(job('promotion-audit-recipes'), /--version/);
  assert.match(job('promotion-audit-recipes'), /skipped-run-locally/);
  // The summary keeps a failed table write visible and still publishes the contact sheet.
  assert.match(job('promotion-audit-summary'), /- name: Write the audit table to the run summary\n\s+shell: bash\n/);
  for (const name of ['Merge the journey contact sheets', 'Upload the journey contact sheet']) {
    assert.match(job('promotion-audit-summary'), new RegExp(`- name: ${name}\\n\\s+if: always\\(\\)\\n`), name);
  }
});
