import assert from 'node:assert/strict';
import { cp, mkdtemp, writeFile } from 'node:fs/promises';
import { tmpdir } from 'node:os';
import { join } from 'node:path';
import { pathToFileURL } from 'node:url';

const repoRoot = process.cwd();
const tempRoot = await mkdtemp(join(tmpdir(), 'gbdraw-legend-sync-'));
await cp(join(repoRoot, 'gbdraw', 'web', 'js', 'app'), join(tempRoot, 'app'), { recursive: true });
await writeFile(join(tempRoot, 'package.json'), '{"type":"module"}\n', 'utf8');

const {
  SPECIFIC_COLOR_FILE_OWNER,
  buildLegendIntents,
  diffLegendIntents
} = await import(pathToFileURL(join(tempRoot, 'app', 'specific-color-rules.js')));
const { parseTransformXY } = await import(
  pathToFileURL(join(tempRoot, 'app', 'legend', 'utils.js'))
);

assert.equal(SPECIFIC_COLOR_FILE_OWNER, 'specific-color-file');
assert.deepEqual(parseTransformXY('translate(12.5,-3.25)'), { x: 12.5, y: -3.25 });
assert.deepEqual(parseTransformXY('translate(.5 2e1)'), { x: 0.5, y: 20 });

const rules = [
  { feat: 'CDS', qual: 'gene', val: 'a', color: '#112233', cap: 'Shared' },
  { feat: 'CDS', qual: 'product', val: 'b', color: '#112233', cap: 'Shared' }
];
const desired = buildLegendIntents(rules).intents;
assert.deepEqual(desired, [{ caption: 'Shared', color: '#112233' }]);

const first = diffLegendIntents([], desired);
assert.deepEqual(first, {
  add: [{ caption: 'Shared', color: '#112233' }],
  update: [],
  remove: [],
  unchanged: []
});
const second = diffLegendIntents(first.add, desired);
assert.deepEqual(second, {
  add: [],
  update: [],
  remove: [],
  unchanged: [{ caption: 'Shared', color: '#112233' }]
});

assert.deepEqual(
  buildLegendIntents([
    { feat: 'CDS', qual: 'gene', val: 'a', color: '#112233', cap: 'Historical' },
    { feat: 'CDS', qual: 'gene', val: 'b', color: '#445566', cap: 'Historical' }
  ], { conflictPolicy: 'last-wins' }),
  {
    intents: [{ caption: 'Historical', color: '#445566' }],
    conflicts: [{
      caption: 'Historical',
      previousColor: '#112233',
      nextColor: '#445566',
      ruleIndex: 1
    }]
  }
);
