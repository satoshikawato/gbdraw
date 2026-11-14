// SELECTIVE_CI.md "Fresh guards at the gate": the fast structure guards that `PR / gate` runs
// again on the pull request merged with dev as dev is when the gate runs. The other jobs of a run
// test the run's own merge commit, and a rerun of the gate alone, as after dev turns green, keeps
// their results; without this, a rule that dev gained after the run's base never reaches the pull
// request before it merges.
//
//   node tools/ci-impact-fresh-guards.mjs     (in the merged working tree; exit 1 when one fails)
//
// Each guard reads files of the working tree only, needs no installed gbdraw, and takes seconds.
import { spawn } from 'node:child_process';
import { appendFileSync } from 'node:fs';
import { pathToFileURL } from 'node:url';

export const FRESH_GUARD_FILES = Object.freeze([
  'tests/test_python_structure_ratchet.py',
  'tests/web/function-size-ratchet.test.mjs',
  'tests/web/explicit-any-ratchet.test.mjs',
  'tests/web/owner-graph-baseline.test.mjs'
]);

// The gate installs no pytest: a Python guard file runs each of its module-level test_* functions,
// which take no argument, with the standard library alone.
const PYTHON_TEST_FUNCTIONS = [
  'import runpy, sys, traceback',
  'failed = []',
  'for name, test in runpy.run_path(sys.argv[1]).items():',
  '    if name.startswith("test_") and callable(test):',
  '        try:',
  '            test()',
  '        except Exception:',
  '            failed.append(name)',
  '            traceback.print_exc()',
  'print(f"failed: {failed}" if failed else "all test functions passed")',
  'sys.exit(1 if failed else 0)'
].join('\n');

/** @returns {[string, string[]]} The command and its arguments that run one guard file. */
export const guardCommand = (file) => (file.endsWith('.py')
  ? ['python3', ['-I', '-c', PYTHON_TEST_FUNCTIONS, file]]
  : [process.execPath, ['--test', file]]);

/**
 * @typedef {object} GuardResult
 * @property {string} file
 * @property {boolean} passed
 * @property {number} seconds
 * @property {string} output Its standard output and error, interleaved.
 */

/** @returns {Promise<GuardResult>} */
const runGuard = (file, cwd) => new Promise((resolve) => {
  const started = Date.now();
  const [command, args] = guardCommand(file);
  const child = spawn(command, args, { cwd, stdio: ['ignore', 'pipe', 'pipe'] });
  let output = '';
  child.stdout.on('data', (chunk) => { output += chunk; });
  child.stderr.on('data', (chunk) => { output += chunk; });
  const finish = (passed, extra = '') => resolve({
    file, passed, seconds: (Date.now() - started) / 1000, output: output + extra
  });
  child.on('error', (error) => finish(false, `\n${error.message}\n`));
  child.on('close', (code) => finish(code === 0));
});

/** @returns {Promise<GuardResult[]>} Every guard of FRESH_GUARD_FILES, run at once in cwd. */
export const runFreshGuards = (cwd = process.cwd()) => Promise.all(FRESH_GUARD_FILES.map((file) => runGuard(file, cwd)));

const main = async () => {
  const results = await runFreshGuards();
  const rows = ['| Guard | Result | Seconds |', '| --- | --- | --- |'];
  for (const { file, passed, seconds, output } of results) {
    console.log(`::group::${file}: ${passed ? 'passed' : 'FAILED'} in ${seconds.toFixed(1)} s`);
    console.log(output.trimEnd());
    console.log('::endgroup::');
    rows.push(`| \`${file}\` | ${passed ? 'passed' : 'failed'} | ${seconds.toFixed(1)} |`);
  }
  const failed = results.filter(({ passed }) => !passed);
  for (const { file } of failed) {
    console.log(`::error title=Fresh guard failed::${file} fails on this pull request merged with the current dev. `
      + 'Update the branch from dev and fix it (docs/internal/SELECTIVE_CI.md "Fresh guards at the gate").');
  }
  if (process.env.GITHUB_STEP_SUMMARY) appendFileSync(process.env.GITHUB_STEP_SUMMARY, `${rows.join('\n')}\n\n`);
  process.exitCode = failed.length ? 1 : 0;
};

if (process.argv[1] && import.meta.url === pathToFileURL(process.argv[1]).href) await main();
