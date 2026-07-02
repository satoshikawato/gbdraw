import assert from 'node:assert/strict';
import { mkdir, mkdtemp, readFile, writeFile } from 'node:fs/promises';
import { tmpdir } from 'node:os';
import { join } from 'node:path';
import { pathToFileURL } from 'node:url';

const repoRoot = process.cwd();
const sourcePath = join(repoRoot, 'gbdraw', 'web', 'js', 'app', 'history-inputs.js');
const tempDir = await mkdtemp(join(tmpdir(), 'gbdraw-history-inputs-'));
await writeFile(join(tempDir, 'package.json'), '{"type":"module"}\n', 'utf8');
await mkdir(join(tempDir, 'app'), { recursive: true });
await writeFile(join(tempDir, 'app', 'history-inputs.js'), await readFile(sourcePath, 'utf8'), 'utf8');

const { isIgnoredTarget } = await import(pathToFileURL(join(tempDir, 'app', 'history-inputs.js')));

const targetMatching = (attribute) => ({
  closest: (selector) => (String(selector).includes(attribute) ? {} : null)
});

assert.equal(isIgnoredTarget(null), false);
assert.equal(isIgnoredTarget({ closest: () => null }), false);
assert.equal(isIgnoredTarget(targetMatching('data-history-ignore')), true);
assert.equal(isIgnoredTarget(targetMatching('data-history-managed')), true);
assert.equal(isIgnoredTarget(targetMatching('data-history-scope="transient"')), true);

console.log('history input tests passed');
