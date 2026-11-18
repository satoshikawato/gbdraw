// The one owner of the temporary folders that tests make. A process keeps its
// folders in a folder of its own,
// <os.tmpdir()>/gbdraw-tests-<uid>/<scope>-<pid>-<start time>, and removes it
// (and gbdraw-tests-<uid> once that is empty) when it exits: after the last
// test, after a failed test, an uncaught error or process.exit(). A process
// that a signal or a crash ends cannot remove it, so each process, when it
// makes its first folder, removes the folders of the processes in its scope
// that have ended. `node --test` runs each file in its own process and
// Playwright runs each worker as one (tests/web/temp-dir.test.mjs).
//
// - The root belongs to one user, so a folder another user left never blocks
//   this user's tests.
// - The scope is the host name and, on Linux, the PID namespace: the processes
//   whose pids this process can check.
// - The start time gives a reused pid a new folder, so a removal of an ended
//   process's folders never takes the folders of a later process with its pid.
//
// The helper handles no signal: a signal ends a test process as it did before.
// Node 26 can abort a `node --test` file that handles one while its run is
// stopped (OV-3920).
const { mkdirSync, mkdtempSync, readdirSync, readlinkSync, rmSync, rmdirSync } = require('node:fs');
const { mkdtemp } = require('node:fs/promises');
const { hostname, tmpdir } = require('node:os');
const { join } = require('node:path');

/** @param {string} folder */
const removeQuietly = (folder) => {
  try {
    rmSync(folder, { recursive: true, force: true, maxRetries: 2 });
  } catch {
    // A child process may still write into it: the next process removes it.
  }
};

/** @param {number} pid */
const hasEnded = (pid) => {
  try {
    process.kill(pid, 0);
    return false;
  } catch (error) {
    return /** @type {NodeJS.ErrnoException} */ (error).code === 'ESRCH';
  }
};

/** @returns {string} The PID namespace's number on Linux, '' elsewhere. */
const pidNamespace = () => {
  try {
    return readlinkSync('/proc/self/ns/pid').replace(/\D/g, '');
  } catch {
    return '';
  }
};

/** @type {string | null} */
let processFolder = null;

/** @returns {string} */
const ownFolder = () => {
  if (processFolder !== null) return processFolder;
  // Windows has no uid and gives each user a temporary folder of its own.
  const root = join(tmpdir(), `gbdraw-tests-${process.getuid ? process.getuid() : 'user'}`);
  const scope = `${hostname()}.${pidNamespace()}`;
  const folder = join(root, `${scope}-${process.pid}-${Math.round(performance.timeOrigin)}`);
  mkdirSync(folder, { recursive: true });
  for (const name of readdirSync(root)) {
    const owner = /^(.+)-([1-9]\d*)-\d+$/.exec(name);
    if (owner && owner[1] === scope && hasEnded(Number(owner[2]))) removeQuietly(join(root, name));
  }
  process.once('exit', () => {
    removeQuietly(folder);
    try {
      rmdirSync(root);
    } catch {
      // Another process still keeps its folders there.
    }
  });
  processFolder = folder;
  return folder;
};

/** @param {string} prefix @returns {string} */
const makeTempDirSync = (prefix) => mkdtempSync(join(ownFolder(), prefix));

/** @param {string} prefix @returns {Promise<string>} */
const makeTempDir = (prefix) => mkdtemp(join(ownFolder(), prefix));

module.exports = { makeTempDir, makeTempDirSync };
