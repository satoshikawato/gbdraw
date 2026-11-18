// The one owner of the temporary folders that tests make. A process keeps its
// folders in a folder of its own, <os.tmpdir()>/gbdraw-tests/<host>-<pid>, and
// removes it (and gbdraw-tests once that is empty) when it exits: after the
// last test, after a failed test, an uncaught error or process.exit(). A
// process that a signal or a crash ends cannot remove it, so each process, when
// it makes its first folder, removes the folders of the processes on this host
// that have ended. `node --test` runs each file in its own process and
// Playwright runs each worker as one (tests/web/temp-dir.test.mjs).
//
// The helper handles no signal: a signal ends a test process as it did before.
// Node 26 can abort a `node --test` file that handles one while its run is
// stopped (OV-3920).
const { mkdirSync, mkdtempSync, readdirSync, rmSync, rmdirSync } = require('node:fs');
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

/** @type {string | null} */
let processFolder = null;

/** @returns {string} */
const ownFolder = () => {
  if (processFolder !== null) return processFolder;
  const root = join(tmpdir(), 'gbdraw-tests');
  const host = hostname();
  const folder = join(root, `${host}-${process.pid}`);
  mkdirSync(folder, { recursive: true });
  for (const name of readdirSync(root)) {
    const owner = /^(.+)-([1-9]\d*)$/.exec(name);
    if (owner && owner[1] === host && Number(owner[2]) !== process.pid && hasEnded(Number(owner[2]))) {
      removeQuietly(join(root, name));
    }
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
