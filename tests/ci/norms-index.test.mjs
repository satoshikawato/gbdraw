import assert from 'node:assert/strict';
import { spawnSync } from 'node:child_process';
import { mkdtempSync, readFileSync, rmSync, writeFileSync } from 'node:fs';
import { tmpdir } from 'node:os';
import { dirname, join } from 'node:path';
import test from 'node:test';
import { fileURLToPath } from 'node:url';
import { buildNormsIndex, citedAnchors, headingsOf } from '../../tools/generate_norms_index.mjs';

// docs/internal/NORMS_INDEX.md is generated from the guard registry, the rule
// headings, the CI routing policy, the Product Contract, and the tests that cite
// it. This file runs in the CI impact job, which every pull request and dev push
// runs, so a change to any of those inputs regenerates the index in the same
// pull request.
// Record ids and anchors below are built at run time: the generator reads this
// file, and a literal one would enter the index.
const ROOT = join(dirname(fileURLToPath(import.meta.url)), '..', '..');
const read = (file) => readFileSync(join(ROOT, file), 'utf8');
const registry = () => JSON.parse(read('tools/web-design-rule-guards.json'));
const RECORD = (number) => `PD-OI-${number}`;
const ANCHOR = (document, slug) => `${document}${'#'}${slug}`;

test('docs/internal/NORMS_INDEX.md is the generator output', () => {
  const result = spawnSync(process.execPath, ['tools/generate_norms_index.mjs', '--check'], { cwd: ROOT, encoding: 'utf8' });
  assert.equal(result.status, 0, `${result.stdout}${result.stderr}`);
});

test('the index refuses a declared job that a pull request changing the guard does not require', () => {
  const changed = registry();
  const guard = changed.rules.find(({ id }) => id === 'R1').guards.find(({ path }) => path.endsWith('.test.mjs'));
  const declared = guard.job;
  guard.job = 'core-pr';
  const { problems } = buildNormsIndex({ registry: changed });
  assert.equal(problems.length, 1, problems.join('\n'));
  const [, required] = problems[0].match(new RegExp(
    `^R1 ${guard.path}: the registry declares core-pr, but a pull request that changes only this file requires (.+)$`
  ));
  assert.ok(required.split(', ').includes(declared));
  assert.ok(!required.split(', ').includes('core-pr'));
});

test('the index refuses a rule heading that its document does not have', () => {
  const changed = registry();
  changed.rules.find(({ id }) => id === 'TYPE-DEBT').heading = 'Type debt';
  assert.deepEqual(buildNormsIndex({ registry: changed }).problems, [
    'TYPE-DEBT: docs/internal/ARCHITECTURE_FITNESS_FUNCTION_RATCHET.md has no heading "Type debt"'
  ]);
});

test('the index refuses an anchor that is not a heading of an existing document and lists a missing document', () => {
  const tool = 'tools/web-product-impact-map.json';
  const citing = (text) => (file) => read(file) + (file === tool ? text : '');
  const { problems } = buildNormsIndex({ read: citing(` ${ANCHOR('gbdraw/web/CLAUDE.md', 'no-such-heading')}`) });
  assert.deepEqual(problems, [
    `${tool} cites ${ANCHOR('gbdraw/web/CLAUDE.md', 'no-such-heading')}: gbdraw/web/CLAUDE.md has no heading with that anchor`
  ]);
  const absent = buildNormsIndex({ read: citing(` ${ANCHOR('docs/absent.md', 'anything')}`) });
  assert.deepEqual(absent.problems, []);
  assert.match(absent.text, new RegExp(`^- \`${ANCHOR('docs/absent.md', 'anything')}\`: .*\`${tool}\``, 'm'));
});

test('the index lists the tests that cite an id without a record, and the records that mention it', () => {
  const spec = 'tests/web/drawing-list.test.mjs';
  const { text, problems } = buildNormsIndex({
    read: (file) => read(file) + (file === spec ? `\n// ${RECORD(999)}\n` : '')
  });
  assert.deepEqual(problems, []);
  const missing = text.slice(text.indexOf('### Cited ids without a record'));
  assert.match(missing, new RegExp(`^- \`${RECORD(999)}\`: no record\\. No record mentions it\\.\\n  - \`${spec}\`$`, 'm'));
});

test('pull requests that make tests cite two neighbouring records change the index without a merge conflict', () => {
  const spec = 'tests/web/drawing-list.test.mjs';
  const base = buildNormsIndex();
  const number = (id) => Number(id.slice(RECORD('').length));
  const ids = headingsOf(read('docs/internal/OPTION_INTEGRITY_PRODUCT_CONTRACT.md'))
    .filter(({ level }) => level === 3)
    .map(({ text }) => text.match(new RegExp(`^${RECORD('\\d+')}(?=: )`))?.[0])
    .filter(Boolean).sort((left, right) => number(left) - number(right));
  const cited = new Set(read(spec).match(new RegExp(RECORD('\\d+'), 'g')) ?? []);
  const index = ids.findIndex((id, at) => at > 0 && !cited.has(id) && !cited.has(ids[at - 1]));
  const citing = (id) => buildNormsIndex({ read: (file) => read(file) + (file === spec ? `\n// ${id}\n` : '') }).text;
  const directory = mkdtempSync(join(tmpdir(), 'norms-index-'));
  try {
    const [ours, original, theirs] = ['ours', 'base', 'theirs'].map((name) => join(directory, name));
    writeFileSync(ours, citing(ids[index - 1]));
    writeFileSync(original, base.text);
    writeFileSync(theirs, citing(ids[index]));
    const merged = spawnSync('git', ['merge-file', '-p', ours, original, theirs], { encoding: 'utf8' });
    assert.equal(merged.status, 0, `conflicting hunks: ${merged.status}`);
    const listings = (text) => text.split(`\`${spec}\``).length - 1;
    assert.equal(listings(merged.stdout), listings(base.text) + 2);
  } finally {
    rmSync(directory, { recursive: true, force: true });
  }
});

test('heading anchors follow GitHub: rendered text, no punctuation, numbered repeats, no fenced code', () => {
  const slugs = (source) => headingsOf(source).map(({ slug }) => slug);
  assert.deepEqual(slugs('# Leaf tests\n## Leaf tests\n### Leaf tests\n'), ['leaf-tests', 'leaf-tests-1', 'leaf-tests-2']);
  assert.deepEqual(slugs('## `grid_column` and [the link](page.md) 日本語 _em_ — R14: Typed boundaries\n'), [
    'grid_column-and-the-link-日本語-em--r14-typed-boundaries'
  ]);
  assert.deepEqual(slugs('```js\n# inside\n```\n# after\n'), ['after']);
  assert.deepEqual(slugs('````\n```\n# inside\n````\n# after\n'), ['after']);
  assert.deepEqual(slugs('~~~\n```\n# inside\n~~~\n# after\n'), ['after']);
});

test('a cited anchor is read up to its last anchor character, and a URL is not a citation', () => {
  const page = ANCHOR('docs/internal/page.md', '日本語-見出し');
  assert.deepEqual(citedAnchors(`see ${page}, and https://example.org/${ANCHOR('a.md', 'b')}`), [
    { document: 'docs/internal/page.md', slug: '日本語-見出し' }
  ]);
});
