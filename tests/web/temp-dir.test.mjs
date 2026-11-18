// A test's temporary folders are made and removed by tests/web/helpers/temp-dir.cjs
// (OV-3858: tests that removed their own folders, or did not, left 124 GiB in
// one machine's TMPDIR). This file checks that the helper removes them on every
// way a test process ends.
import assert from 'node:assert/strict';
import { spawn } from 'node:child_process';
import { existsSync, mkdirSync, readFileSync, readdirSync, writeFileSync } from 'node:fs';
import { join } from 'node:path';
import { pathToFileURL } from 'node:url';
import test from 'node:test';
import { makeTempDirSync } from './helpers/temp-dir.cjs';

const OWNER = 'tests/web/helpers/temp-dir.cjs';
const HELPER_URL = pathToFileURL(join(process.cwd(), OWNER)).href;

// One test process: two folders made and reported, then the end that
// `TEMP_DIR_CASE` names. 'handled' stands for a Playwright worker, which
// handles SIGTERM itself and exits normally.
const FIXTURE = `
import { writeFileSync } from 'node:fs';
import test from 'node:test';
import { makeTempDir, makeTempDirSync } from ${JSON.stringify(HELPER_URL)};
const end = process.env.TEMP_DIR_CASE;
if (end === 'handled') process.on('SIGTERM', () => process.exit(0));
const made = [makeTempDirSync('case-sync-'), await makeTempDir('case-async-')];
writeFileSync(process.env.TEMP_DIR_REPORT, JSON.stringify({ pid: process.pid, made }));
test('case', async () => {
  if (end === 'fail') throw new Error('failed on purpose');
  if (end === 'throw') setTimeout(() => { throw new Error('uncaught on purpose'); }, 0);
  if (end === 'exit') process.exit(3);
  if (['hang', 'handled'].includes(end)) await new Promise((resolve) => setTimeout(resolve, 60_000));
});
`;
const fixtureDirectory = makeTempDirSync('gbdraw-temp-dir-fixture-');
const fixture = join(fixtureDirectory, 'case.test.mjs');
writeFileSync(fixture, FIXTURE);

/** @param {() => boolean} done @param {string} message */
const waitFor = async (done, message) => {
  for (const started = Date.now(); !done();) {
    assert.ok(Date.now() - started < 10_000, message);
    await new Promise((resolve) => setTimeout(resolve, 20));
  }
};

/** @param {number} pid */
const running = (pid) => {
  try {
    process.kill(pid, 0);
    return true;
  } catch {
    return false;
  }
};

/** @param {string} report @returns {Promise<{ pid: number, made: string[] }>} */
const reported = async (report) => {
  await waitFor(() => existsSync(report), 'the fixture never made its folders');
  return JSON.parse(readFileSync(report, 'utf8'));
};

/**
 * Runs the fixture in its own process group with its own TMPDIR, outside the
 * test runner that runs this file, and stops it the way `stop` says once both
 * folders exist. Returns the folders it made and what its TMPDIR still holds.
 * @param {string} end
 * @param {{ runner?: boolean, stop?: (pid: number) => void }} options
 */
const runFixture = async (end, { runner = true, stop } = {}) => {
  const root = makeTempDirSync(`gbdraw-temp-dir-${end}-`);
  const directory = join(root, 'tmp');
  const report = join(root, 'made.json');
  mkdirSync(directory);
  const env = Object.fromEntries(Object.entries(process.env).filter(([name]) => !name.startsWith('NODE_TEST')));
  const child = spawn(process.execPath, runner ? ['--test', fixture] : [fixture], {
    env: { ...env, TMPDIR: directory, TEMP_DIR_CASE: end, TEMP_DIR_REPORT: report },
    stdio: 'ignore',
    detached: true
  });
  const exited = new Promise((resolve) => child.once('exit', resolve));
  if (stop) {
    await reported(report);
    stop(/** @type {number} */ (child.pid));
  }
  await exited;
  // `node --test` may exit before the test-file process that made the folders.
  const { pid, made } = await reported(report);
  await waitFor(() => !running(pid), 'the test-file process never ended');
  return { made, left: readdirSync(directory), existing: made.filter((folder) => existsSync(folder)) };
};

/** @param {{ made: string[], left: string[], existing: string[] }} run */
const assertRemoved = ({ made, left, existing }) => {
  assert.equal(made.length, 2);
  assert.deepEqual({ left, existing }, { left: [], existing: [] });
};

for (const end of ['pass', 'fail', 'throw', 'exit']) {
  test(`a test file that ends by ${end} leaves no temporary folder`, async () => {
    assertRemoved(await runFixture(end));
  });
}

test('a test file the runner stops with SIGTERM leaves no temporary folder', async () => {
  assertRemoved(await runFixture('hang', { stop: (pid) => process.kill(pid, 'SIGTERM') }));
});

test('a test run stopped with Ctrl-C (SIGINT to its process group) leaves no temporary folder', async () => {
  assertRemoved(await runFixture('hang', { stop: (pid) => process.kill(-pid, 'SIGINT') }));
});

test('a process that handles SIGTERM itself and exits leaves no temporary folder', async () => {
  assertRemoved(await runFixture('handled', { runner: false, stop: (pid) => process.kill(pid, 'SIGTERM') }));
});
