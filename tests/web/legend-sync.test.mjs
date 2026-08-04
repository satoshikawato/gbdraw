import assert from 'node:assert/strict';
import { cp, mkdtemp, readFile, writeFile } from 'node:fs/promises';
import { tmpdir } from 'node:os';
import { join } from 'node:path';
import { pathToFileURL } from 'node:url';

const repoRoot = process.cwd();
const tempRoot = await mkdtemp(join(tmpdir(), 'gbdraw-legend-sync-'));
await cp(join(repoRoot, 'gbdraw', 'web', 'js', 'app'), join(tempRoot, 'app'), { recursive: true });
await cp(join(repoRoot, 'gbdraw', 'web', 'js', 'services'), join(tempRoot, 'services'), { recursive: true });
await writeFile(join(tempRoot, 'package.json'), '{"type":"module"}\n', 'utf8');

const {
  SPECIFIC_COLOR_FILE_OWNER,
  buildLegendIntents,
  diffLegendIntents
} = await import(pathToFileURL(join(tempRoot, 'app', 'specific-color-rules.js')));
const {
  COMPARISON_LEGEND_SELECTOR,
  PAIRWISE_LEGEND_SELECTOR,
  getComparisonLegendGroup,
  getLegendChildById,
  parseTransformXY
} = await import(
  pathToFileURL(join(tempRoot, 'app', 'legend', 'utils.js'))
);
assert.equal(SPECIFIC_COLOR_FILE_OWNER, 'specific-color-file');
assert.deepEqual(parseTransformXY('translate(12.5,-3.25)'), { x: 12.5, y: -3.25 });
assert.deepEqual(parseTransformXY('translate(.5 2e1)'), { x: 0.5, y: 20 });
assert.match(COMPARISON_LEGEND_SELECTOR, /^\[data-gbdraw-role="comparison-legend"\]/);
assert.doesNotMatch(PAIRWISE_LEGEND_SELECTOR, /conservation_identity_legend/);

const comparisonLegend = {
  id: 'pairwise_legend_h',
  getAttribute: (name) => name === 'data-gbdraw-role' ? 'comparison-legend' : null
};
const parent = { children: [comparisonLegend] };
assert.equal(getLegendChildById(parent, 'pairwise_legend'), comparisonLegend);
assert.equal(getComparisonLegendGroup(parent), comparisonLegend);

const svgStylesSource = await readFile(join(tempRoot, 'app', 'svg-styles.js'), 'utf8');
const repositionSource = await readFile(
  join(tempRoot, 'app', 'legend-layout', 'reposition-actions.js'),
  'utf8'
);
const legendLayoutSource = await readFile(join(tempRoot, 'app', 'legend-layout.js'), 'utf8');
const entryActionsSource = await readFile(
  join(tempRoot, 'app', 'legend', 'entry-actions.js'),
  'utf8'
);
const appSetupSource = await readFile(join(tempRoot, 'app', 'app-setup.js'), 'utf8');
const watchersSource = await readFile(join(tempRoot, 'app', 'watchers.js'), 'utf8');
const configSource = await readFile(join(tempRoot, 'services', 'config.js'), 'utf8');
assert.match(svgStylesSource, /querySelectorAll\(PAIRWISE_LEGEND_SELECTOR\)/);
assert.match(repositionSource, /applyCompositionEdit/);
assert.match(repositionSource, /bindCompositionMetadata/);
assert.doesNotMatch(repositionSource, /data-horizontal-viewbox|data-vertical-viewbox/);
assert.doesNotMatch(repositionSource, /0\.025|0\.85|0\.875|0\.75/);
assert.match(
  legendLayoutSource,
  /resetAllPositions[\s\S]+resetCompositionUserDeltas[\s\S]+persistCurrentSvg\(svg\)/
);
assert.match(entryActionsSource, /setLegendGeometryChangedHandler/);
assert.ok((entryActionsSource.match(/onLegendGeometryChanged\(\);/g) || []).length >= 4);
assert.match(
  appSetupSource,
  /setLegendGeometryChangedHandler\(legendLayout\.refreshLegendGeometry\)/
);
assert.match(
  watchersSource,
  /bind composition metadata[\s\S]+captureBaseConfig/
);
assert.match(configSource, /skipCaptureBaseConfig\.value = true;\s+state\.skipPositionReapply\.value = true;\s+applyResultsData/);

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
