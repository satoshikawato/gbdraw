import assert from 'node:assert/strict';
import { readFile, writeFile, mkdtemp } from 'node:fs/promises';
import { tmpdir } from 'node:os';
import { join } from 'node:path';
import { pathToFileURL } from 'node:url';

const repoRoot = process.cwd();
const sourcePath = join(repoRoot, 'gbdraw', 'web', 'js', 'app', 'run-info.js');
const tempDir = await mkdtemp(join(tmpdir(), 'gbdraw-run-info-'));
const modulePath = join(tempDir, 'run-info.mjs');
await writeFile(modulePath, await readFile(sourcePath, 'utf8'), 'utf8');

const {
  buildRunInfo,
  isCliInvocationSessionExportable,
  quoteShellArg
} = await import(pathToFileURL(modulePath));

assert.equal(quoteShellArg('simple.tsv'), 'simple.tsv');
assert.equal(quoteShellArg('two words.tsv'), "'two words.tsv'");
assert.equal(quoteShellArg("Bob's.gbk"), "'Bob'\\''s.gbk'");
assert.equal(quoteShellArg(''), "''");

{
  const info = buildRunInfo({
    mode: 'circular',
    args: ['-o', 'my diagram', '-d', '/combined_d.tsv', '--track_type', 'middle', '--gbk', '/input.gb'],
    fileMetadata: new Map([
      ['/input.gb', { name: 'input file.gbk', slot: 'files.c_gb', kind: 'uploaded' }],
      ['/combined_d.tsv', { name: 'combined_d.tsv', slot: 'generatedFiles.combined_d', kind: 'generated' }]
    ]),
    elapsedMs: 1234.5,
    resultCount: 1,
    startedAtIso: '2026-06-23T00:00:00.000Z'
  });

  assert.equal(
    info.command,
    "gbdraw circular -o 'my diagram' -d combined_d.tsv --track_type middle --gbk 'input file.gbk' -f svg"
  );
  assert.deepEqual(info.invocation.args.slice(-2), ['-f', 'svg']);
  assert.deepEqual(info.helperFiles, [
    { path: '/combined_d.tsv', name: 'combined_d.tsv', slot: 'generatedFiles.combined_d' }
  ]);
  assert.equal(info.reproducibility.level, 'requires-helper-files');
  assert.equal(isCliInvocationSessionExportable(info.invocation), false);
}

{
  const info = buildRunInfo({
    mode: 'linear',
    args: ['--gbk', '/seq_0.gb', '/seq_1.gb', '-b', '/blast_0.txt'],
    fileMetadata: {
      '/seq_0.gb': { name: 'alpha.gb', slot: 'files.linearSeqs[0].gb', kind: 'uploaded' },
      '/seq_1.gb': { name: 'beta.gb', slot: 'files.linearSeqs[1].gb', kind: 'uploaded' },
      '/blast_0.txt': { name: 'alpha_vs_beta.tsv', slot: 'files.linearSeqs[0].blast', kind: 'uploaded' }
    },
    elapsedMs: 900,
    resultCount: 2
  });

  assert.equal(info.command, 'gbdraw linear --gbk alpha.gb beta.gb -b alpha_vs_beta.tsv -f svg');
  assert.deepEqual(info.invocation.fileBindings, [
    { argIndex: 1, slot: 'files.linearSeqs[0].gb', name: 'alpha.gb' },
    { argIndex: 2, slot: 'files.linearSeqs[1].gb', name: 'beta.gb' },
    { argIndex: 4, slot: 'files.linearSeqs[0].blast', name: 'alpha_vs_beta.tsv' }
  ]);
  assert.equal(info.reproducibility.level, 'exact-uploaded-files');
  assert.equal(isCliInvocationSessionExportable(info.invocation), true);
}

