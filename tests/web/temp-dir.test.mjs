// A test's temporary folders are made and removed by tests/web/helpers/temp-dir.cjs
// (OV-3858: tests that removed their own folders, or did not, left 124 GiB in
// one machine's TMPDIR). This file checks that the helper removes them on every
// way a test process ends: at its exit, or, when a signal or a crash ends it,
// in the next process that makes a folder. It also checks that no other test
// file makes one.
import assert from 'node:assert/strict';
import { execFileSync, spawn } from 'node:child_process';
import { existsSync, mkdirSync, readFileSync, readdirSync, writeFileSync } from 'node:fs';
import { basename, dirname, join } from 'node:path';
import { pathToFileURL } from 'node:url';
import test from 'node:test';
import { makeTempDirSync, removeTempDir } from './helpers/temp-dir.cjs';

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

/** @param {number} group */
const killGroup = (group) => {
  try {
    process.kill(-group, 'SIGKILL');
  } catch {
    // The group has ended.
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
 * A fixture that does not get there, or does not end when awaited, is killed so
 * that a failing case does not leave it running for a minute.
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
  const group = /** @type {number} */ (child.pid);
  /** @param {() => boolean} done @param {string} message */
  const waitOrKill = async (done, message) => {
    try {
      await waitFor(done, message);
    } catch (error) {
      killGroup(group);
      throw error;
    }
  };
  await waitOrKill(() => existsSync(report), 'the fixture never made its folders');
  /** @type {{ pid: number, made: string[], listeners: number }} */
  const made = JSON.parse(readFileSync(report, 'utf8'));
  return {
    ...made,
    group,
    // `node --test` may exit before the test-file process it ran.
    ended: () => waitOrKill(
      () => (child.exitCode !== null || child.signalCode !== null) && !running(made.pid),
      'the fixture never ended'
    )
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
  try {
    const next = await startFixture(root, 'pass', { runner: false });
    await next.ended();
    assert.deepEqual(alive.made.filter((folder) => existsSync(folder)), alive.made);
  } finally {
    killGroup(alive.group);
    await alive.ended();
  }
});

test('a test process keeps its folders in a root of its user, in a folder named by its pid and start time', async () => {
  const root = makeRoot();
  const before = Date.now();
  const run = await startFixture(root, 'pass', { runner: false });
  const after = Date.now();
  await run.ended();
  const own = dirname(run.made[0]);
  assert.equal(dirname(own), join(root, 'tmp', `gbdraw-tests-${process.getuid ? process.getuid() : 'user'}`));
  const [, pid, started] = (/-(\d+)-(\d+)$/.exec(basename(own)) ?? []).map(Number);
  assert.equal(pid, run.pid);
  assert.ok(started >= before - 1 && started <= after, `start time ${started} not in [${before}, ${after}]`);
});

test('the next test process removes the folders of ended processes on its host and in its PID namespace only', async () => {
  const root = makeRoot();
  const ended = await startFixture(root, 'pass', { runner: false });
  await ended.ended();
  // This file's processes and the fixtures they start share a host and a PID namespace.
  const own = dirname(makeTempDirSync('gbdraw-temp-dir-scope-'));
  const scope = basename(own).replace(/-\d+-\d+$/, '');
  const userRoot = join(root, 'tmp', basename(dirname(own)));
  const sameScope = join(userRoot, `${scope}-${ended.pid}-1`);
  const otherScope = join(userRoot, `other-${scope}-${ended.pid}-1`);
  for (const folder of [sameScope, otherScope]) mkdirSync(join(folder, 'left'), { recursive: true });
  const next = await startFixture(root, 'pass', { runner: false });
  await next.ended();
  assert.deepEqual({ sameScope: existsSync(sameScope), otherScope: existsSync(otherScope) }, { sameScope: false, otherScope: true });
});

test('a folder removed before its process ends is gone, and only the helper\'s folders can be removed', () => {
  const folder = makeTempDirSync('gbdraw-temp-dir-early-');
  removeTempDir(folder);
  assert.equal(existsSync(folder), false);
  assert.throws(() => removeTempDir(folder), /is not a temporary folder of this process/);
  const other = join(makeRoot(), 'not-made-by-the-helper');
  mkdirSync(other);
  assert.throws(() => removeTempDir(other), /is not a temporary folder of this process/);
  assert.equal(existsSync(other), true);
});

// Node's folder makers in every form (the promise and sync ones, the
// disposable ones of Node 24, imported or as a member) and a path built on the
// system temporary directory itself. This file is scanned too. Python's
// tempfile.mkdtemp in a snippet runs under the evaluator's TMPDIR, a folder
// the helper made (tests/web/helpers/python-evaluator.mjs).
const OUTSIDE_THE_HELPER = /(?<!tempfile\.)\bmkdtemp\w*|\btmpdir\(/;

test('only the helper makes temporary folders for tests', () => {
  const files = execFileSync('git', ['ls-files', '-z', '--', 'tests'], { encoding: 'utf8' })
    .split('\0')
    .filter((path) => /\.(?:c|m)?js$/.test(path) && path !== OWNER);
  const offenders = files.filter((path) => OUTSIDE_THE_HELPER.test(readFileSync(path, 'utf8')));
  assert.deepEqual(offenders, [],
    `Make a test's temporary folder with makeTempDir or makeTempDirSync from ${OWNER}:\n${offenders.join('\n')}`);
});
