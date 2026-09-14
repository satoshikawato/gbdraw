import assert from 'node:assert/strict';
import { readFile, mkdtemp, writeFile, rm } from 'node:fs/promises';
import { tmpdir } from 'node:os';
import { join } from 'node:path';
import { spawnSync } from 'node:child_process';
import vm from 'node:vm';
import { runLosatPairWasi } from '../../gbdraw/web/js/services/losat-runtime.js';

// NCBI c++/src/algo/blast/blastinput/cmdline_flags.cpp:46-75: kArgQuery("query"),
// kArgSubject("subject"), kArgNumThreads("num_threads").
// Exercise the actual two command argv builders. The shim captures serial argv;
// the threaded worker's pure builder is evaluated without starting a worker.
const source = await readFile(new URL('../../gbdraw/web/js/workers/losat-threaded-worker.js', import.meta.url), 'utf8');
const builder = vm.runInNewContext(`${source.slice(source.indexOf('const hasNumThreadsArg'), source.indexOf('const getChildWorkerCount'))}\nbuildLosatArgs`);
let captured;
class StubFile { constructor() {} }
class StubWasi { constructor(args) { captured = args; } start() { return 0; } }
const shim = { WASI: StubWasi, File: StubFile, OpenFile: StubFile, PreopenDirectory: StubFile, ConsoleStdout: StubFile };
const originalInstantiate = WebAssembly.instantiate;
const directory = await mkdtemp(join(tmpdir(), 'gbdraw-cli-v2-'));
try {
  await writeFile(join(directory, 'query.fa'), '>q\nACGTACGTGGTACCGTACGTAACCGGTTAACCGGTTAACCGGTT\n');
  await writeFile(join(directory, 'subject.fa'), '>s\nACGTACGTGGTACCGTACGTAACCGGTTAACCGGTTAACCGGTT\n');
  WebAssembly.instantiate = async () => ({ exports: {} });
  for (const program of ['blastn', 'blastp', 'tblastx']) {
    const extraArgs = program === 'blastn' ? ['-task', 'blastn'] : program === 'tblastx' ? ['-query_gencode', '1', '-db_gencode', '1'] : ['-max_hsps', '1', '-max_target_seqs', '5'];
    await runLosatPairWasi({ program, queryFasta: '>q\nACGT\n', subjectFasta: '>s\nACGT\n', extraArgs, wasiShim: shim, wasmModule: {} });
    const serial = [...captured];
    const threaded = [...builder({ program, outfmt: '6', extraArgs, threadsPerJob: 2 })];
    for (const [args, threads] of [[serial, '1'], [threaded, '2']]) {
      assert.equal(args[args.indexOf('-num_threads') + 1], threads);
      assert.equal(args[args.indexOf('-outfmt') + 1], '6');
      assert.ok(args.every((arg) => !arg.startsWith('--')));
      if (process.env.LOSAT_CLI_V2_BINARY) {
        const run = spawnSync(process.env.LOSAT_CLI_V2_BINARY, args.slice(1), { cwd: directory, encoding: 'utf8' });
        assert.equal(run.status, 0, `${args.join(' ')}\n${run.stderr}`);
      }
    }
    assert.equal(builder({ program, outfmt: '6', extraArgs: ['-num_threads', '3'], threadsPerJob: 2 }).filter((arg) => arg === '-num_threads').length, 1);
  }
} finally {
  WebAssembly.instantiate = originalInstantiate;
  await rm(directory, { recursive: true, force: true });
}
console.log('CLI v2 serial/threaded argv passed for all three programs');
