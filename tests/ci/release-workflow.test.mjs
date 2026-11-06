import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import test from 'node:test';

const source = readFileSync(new URL('../../.github/workflows/release.yml', import.meta.url), 'utf8');
const build = source.split('  build-release-distributions:\n')[1].split('  release-evidence:\n')[0];
const evidence = source.split('  release-evidence:\n')[1].split('  pypi-publish:\n')[0];
const publish = source.split('  pypi-publish:\n')[1];

test('publishing starts only when a GitHub Release is published', () => {
  assert.match(source, /\non:\n  release:\n    types: \[published\]\n\n/);
  assert.doesNotMatch(source, /workflow_dispatch|workflow_call|pull_request|\n  push:|branches:/);
  assert.match(source, /\npermissions: \{\}/);
  assert.equal((source.match(/id-token: write/g) || []).length, 1);
  assert.doesNotMatch(source, /secrets\.|password:|username:|user:|PYPI_TOKEN|TWINE_|repository-url:|continue-on-error|retry|always\(\)/);
  assert.equal((source.match(/^  [\w-]+:\n    (?:needs|runs-on|name):/gm) || []).length, 3);
});

test('publication waits for release-tier evidence that carries to the tagged commit', () => {
  assert.match(evidence, /permissions:\n      contents: read\n      actions: read\n    steps:/);
  assert.match(evidence, /ref: \$\{\{ github.sha \}\}/);
  assert.match(evidence, /fetch-depth: 0/);
  assert.match(evidence, /persist-credentials: false/);
  assert.match(evidence, /GH_TOKEN: \$\{\{ github.token \}\}/);
  assert.match(evidence, /run: node tools\/check-release-evidence.mjs --tag-commit "\$GITHUB_SHA" --repository "\$GITHUB_REPOSITORY"\n/);
  assert.doesNotMatch(evidence, /\bif:|environment:|id-token/);
});

test('build checks exact source before one standard build and package verification', () => {
  assert.match(build, /permissions:\n      contents: read/);
  assert.match(build, /ref: \$\{\{ github.sha \}\}/);
  assert.match(build, /fetch-depth: 0/);
  assert.match(build, /persist-credentials: false/);
  const commands = [
    'python tools/check_release_source.py',
    'python tools/prepare_browser_wheel.py',
    'python -m build --outdir "$RUNNER_TEMP/release-dist"',
    'python -m twine check --strict',
    'tools/verify_gui_offline.py inspect-distributions "$RUNNER_TEMP/release-dist"',
    '-I installed_package_smoke.py',
    'uses: actions/upload-artifact@'
  ];
  let previous = -1;
  for (const command of commands) {
    const index = build.indexOf(command);
    assert.ok(index > previous, command);
    previous = index;
  }
  assert.equal((build.match(/python -m build /g) || []).length, 1);
  assert.match(build, /artifact-id: \$\{\{ steps.distributions.outputs.artifact-id \}\}/);
  assert.match(build, /path: \$\{\{ runner.temp \}\}\/release-dist\//);
  assert.match(build, /if-no-files-found: error/);
});

test('credential job only downloads the same verified artifact and publishes with OIDC', () => {
  assert.match(publish, /needs: \[build-release-distributions, release-evidence\]/);
  assert.match(publish, /environment: pypi/);
  assert.match(publish, /permissions:\n      id-token: write\n    steps:/);
  assert.equal((publish.match(/- name:/g) || []).length, 2);
  assert.doesNotMatch(publish, /\brun:|checkout|setup-python|build --|pytest|git |contents:|github-token:|run-id:|repository:|\bif:/);
  assert.match(publish, /uses: actions\/download-artifact@/);
  assert.match(publish, /artifact-ids: \$\{\{ needs.build-release-distributions.outputs.artifact-id \}\}/);
  assert.match(publish, /path: dist\//);
  assert.match(publish, /digest-mismatch: error/);
  assert.match(publish, /uses: pypa\/gh-action-pypi-publish@/);
  assert.match(publish, /packages-dir: dist\//);
  assert.match(publish, /skip-existing: false/);
});
