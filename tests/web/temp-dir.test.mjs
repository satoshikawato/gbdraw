// A test's temporary folders are made and removed by tests/web/helpers/temp-dir.cjs
// (OV-3858: tests that removed their own folders, or did not, left 124 GiB in
// one machine's TMPDIR). This file checks that the helper removes them on every
// way a test process ends: at its exit, or, when a signal or a crash ends it,
// in the next process that makes a folder.
import assert from 'node:assert/strict';
import { spawn } from 'node:child_process';
import { existsSync, mkdirSync, readFileSync, readdirSync, writeFileSync } from 'node:fs';
import { join } from 'node:path';
import { pathToFileURL } from 'node:url';
import test from 'node:test';
import { makeTempDirSync } from './helpers/temp-dir.cjs';

const OWNER = 'tests/web/helpers/temp-dir.cjs';
const HELPER_URL = pathToFileURL(join(process.cwd(), OWNER)).href;

// One test process: two folders made and reported with the number of signal
// listeners the helper added, then the end that `TEMP_DIR_CASE` names.
// 'handled' stands for a Playwright worker, which handles SIGTERM itself and
// exits normally.
const FIXTURE = `
import { renameSync, writeFileSync } from 'node:fs';
import test from 'node:test';
import { makeTempDir, makeTempDirSync } from ${JSON.stringify(HELPER_URL)};
const end = process.env.TEMP_DIR_CASE;
const made = [makeTempDirSync('case-sync-'), await makeTempDir('case-async-')];
const listeners = ['SIGINT', 'SIGTERM', 'SIGHUP'].reduce((sum, signal) => sum + process.listenerCount(signal), 0);
if (end === 'handled') process.on('SIGTERM', () => process.exit(0));
const report = process.env.TEMP_DIR_REPORT;
writeFileSync(report + '.part', JSON.stringify({ pid: process.pid, made, listeners }));
renameSync(report + '.part', report);
test('case', async () => {
  if (end === 'fail') throw new Error('failed on purpose');
  if (end === 'throw') setTimeout(() => { throw new Error('uncaught on purpose'); }, 0);
  if (end === 'exit') process.exit(3);
  if (['hang', 'handled'].includes(end)) await new Promise((resolve) => setTimeout(resolve, 60_000));
});
`;
const fixture = join(makeTempDirSync('gbdraw-temp-dir-fixture-'), 'case.test.mjs');
writeFileSync(fixture, FIXTURE);

/** @param {() => boolean} done @param {string} message */
const waitFor = async (done, message) => {
  for (const started = Date.now(); !done();) {
    assert.ok(Date.now() - started < 30_000, message);
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

/** A folder whose `tmp` folder serves as the TMPDIR of the fixtures of one case. */
const makeRoot = () => {
  const root = makeTempDirSync('gbdraw-temp-dir-case-');
  mkdirSync(join(root, 'tmp'));
  return root;
};

/**
 * @typedef {object} FixtureRun
 * @property {number} pid The test-file process, which made the folders.
 * @property {number} group The process group the fixture runs in.
 * @property {string[]} made
 * @property {number} listeners
 * @property {() => Promise<void>} ended Resolves when the test-file process has ended.
 */

let reports = 0;

/**
 * Starts the fixture in its own process group with `<root>/tmp` as its TMPDIR,
 * outside the test runner that runs this file. Resolves once its folders exist.
 * @param {string} root
 * @param {string} end
 * @param {{ runner?: boolean }} [options]
 * @returns {Promise<FixtureRun>}
 */
const startFixture = async (root, end, { runner = true } = {}) => {
  const report = join(root, `made-${reports += 1}.json`);
  const env = Object.fromEntries(Object.entries(process.env).filter(([name]) => !name.startsWith('NODE_TEST')));
  const child = spawn(process.execPath, runner ? ['--test', fixture] : [fixture], {
    env: { ...env, TMPDIR: join(root, 'tmp'), TEMP_DIR_CASE: end, TEMP_DIR_REPORT: report },
    stdio: 'ignore',
    detached: true
  });
  const exited = new Promise((resolve) => child.once('exit', resolve));
  await waitFor(() => existsSync(report), 'the fixture never made its folders');
  /** @type {{ pid: number, made: string[], listeners: number }} */
  const made = JSON.parse(readFileSync(report, 'utf8'));
  return {
    ...made,
    group: /** @type {number} */ (child.pid),
    ended: async () => {
      await exited;
      // `node --test` may exit before the test-file process it ran.
      await waitFor(() => !running(made.pid), 'the test-file process never ended');
    }
  };
};

/** @param {string} root @param {FixtureRun[]} runs */
const assertRemoved = (root, ...runs) => {
  for (const run of runs) assert.equal(run.made.length, 2);
  assert.deepEqual({
    existing: runs.flatMap((run) => run.made).filter((folder) => existsSync(folder)),
    left: readdirSync(join(root, 'tmp'))
  }, { existing: [], left: [] });
};

for (const end of ['pass', 'fail', 'throw', 'exit']) {
  test(`a test file that ends by ${end} removes its temporary folders`, async () => {
    const root = makeRoot();
    const run = await startFixture(root, end);
    await run.ended();
    assertRemoved(root, run);
  });
}

test('a process that handles SIGTERM itself and exits removes its temporary folders', async () => {
  const root = makeRoot();
  const run = await startFixture(root, 'handled', { runner: false });
  process.kill(run.pid, 'SIGTERM');
  await run.ended();
  assertRemoved(root, run);
});

test('the helper handles no signal, so a signal ends a test process as before', async () => {
  const root = makeRoot();
  const run = await startFixture(root, 'pass');
  await run.ended();
  assert.equal(run.listeners, 0);
});

for (const [how, stop] of /** @type {Array<[string, (run: FixtureRun) => void]>} */ ([
  ['SIGKILL', (run) => process.kill(run.pid, 'SIGKILL')],
  ['Ctrl-C (SIGINT to its process group)', (run) => process.kill(-run.group, 'SIGINT')]
])) {
  test(`the folders of a test file stopped by ${how} are removed by the next test process`, async () => {
    const root = makeRoot();
    const stopped = await startFixture(root, 'hang');
    stop(stopped);
    await stopped.ended();
    assert.deepEqual(stopped.made.filter((folder) => existsSync(folder)), stopped.made);
    const next = await startFixture(root, 'pass', { runner: false });
    await next.ended();
    assertRemoved(root, stopped, next);
  });
}

test('the folders of a running test process stay when another test process starts', async () => {
  const root = makeRoot();
  const alive = await startFixture(root, 'hang');
  const next = await startFixture(root, 'pass', { runner: false });
  await next.ended();
  assert.deepEqual(alive.made.filter((folder) => existsSync(folder)), alive.made);
  process.kill(-alive.group, 'SIGKILL');
  await alive.ended();
});
