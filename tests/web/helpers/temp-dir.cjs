// The one owner of the temporary folders that tests make. Each folder is made
// under os.tmpdir() and removed when the process that made it exits: after the
// last test, after a failed test or an uncaught error, and when a signal stops
// the process. `node --test` runs each file in its own process and Playwright
// runs each worker as one, so no folder outlives its test file or worker
// (tests/web/temp-dir.test.mjs).
const { mkdtempSync, rmSync } = require('node:fs');
const { mkdtemp } = require('node:fs/promises');
const { tmpdir } = require('node:os');
const { join } = require('node:path');

/** @type {Set<string>} */
const folders = new Set();
let watching = false;

// Each folder is removed on its own: a folder that cannot be removed (a child
// process still writing into it) must not change the exit status, keep the
// other folders, or skip the other 'exit' listeners.
const removeTempDirs = () => {
  for (const folder of folders) {
    try {
      rmSync(folder, { recursive: true, force: true, maxRetries: 2 });
    } catch {
      // Left for the next cleanup of TMPDIR.
    }
  }
  folders.clear();
};

// A signal nobody handles ends the process without its 'exit' event, so remove
// the folders first and then let the signal end the process as before. A
// Playwright worker handles SIGINT and SIGTERM itself and exits normally.
/** @param {NodeJS.Signals} signal */
const removeAndResignal = (signal) => {
  removeTempDirs();
  process.kill(process.pid, signal);
};

/** @param {string} folder @returns {string} */
const own = (folder) => {
  if (!watching) {
    watching = true;
    process.once('exit', removeTempDirs);
    for (const signal of /** @type {NodeJS.Signals[]} */ (['SIGINT', 'SIGTERM', 'SIGHUP'])) {
      if (process.listenerCount(signal) === 0) process.once(signal, removeAndResignal);
    }
  }
  folders.add(folder);
  return folder;
};

/** @param {string} prefix @returns {string} */
const makeTempDirSync = (prefix) => own(mkdtempSync(join(tmpdir(), prefix)));

/** @param {string} prefix @returns {Promise<string>} */
const makeTempDir = async (prefix) => own(await mkdtemp(join(tmpdir(), prefix)));

module.exports = { makeTempDir, makeTempDirSync };
