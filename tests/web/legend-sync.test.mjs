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
const { createLegendRepositionActions } = await import(
  pathToFileURL(join(tempRoot, 'app', 'legend-layout', 'reposition-actions.js'))
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
const entryActionsSource = await readFile(
  join(tempRoot, 'app', 'legend', 'entry-actions.js'),
  'utf8'
);
const appSetupSource = await readFile(join(tempRoot, 'app', 'app-setup.js'), 'utf8');
const watchersSource = await readFile(join(tempRoot, 'app', 'watchers.js'), 'utf8');
const configSource = await readFile(join(tempRoot, 'services', 'config.js'), 'utf8');
assert.match(svgStylesSource, /querySelectorAll\(PAIRWISE_LEGEND_SELECTOR\)/);
assert.match(repositionSource, /querySelector\(PAIRWISE_LEGEND_SELECTOR\)/);
assert.match(repositionSource, /generatedLegendPosition\.value \|\| form\.legend/);
assert.match(entryActionsSource, /setLegendGeometryChangedHandler/);
assert.ok((entryActionsSource.match(/onLegendGeometryChanged\(\);/g) || []).length >= 4);
assert.match(
  appSetupSource,
  /setLegendGeometryChangedHandler\(legendLayout\.refreshLegendGeometry\)/
);
assert.match(
  watchersSource,
  /captureBaseConfig[\s\S]+mode\.value === 'circular' && !shouldSkipPositionReapply[\s\S]+refreshLegendGeometry/
);
assert.match(configSource, /skipPositionReapply\.value = true;\s+applyResultsData/);

const ref = (value) => ({ value });
const attributes = new Map([['viewBox', '0 0 1000 600']]);
const legendAttributes = new Map([['transform', 'translate(1015, 240)']]);
const legendGroup = {
  style: {},
  getAttribute: (name) => legendAttributes.get(name) || null,
  setAttribute: (name, value) => legendAttributes.set(name, value),
  getBBox: () => ({ x: 0, y: 0, width: 300, height: 120 })
};
let hasLegend = true;
const diagramAttributes = new Map([['transform', 'translate(10, 20)']]);
const diagramElement = {
  id: 'record_1',
  getAttribute: (name) => diagramAttributes.get(name) || null,
  setAttribute: (name, value) => diagramAttributes.set(name, value)
};
const svg = {
  getAttribute: (name) => attributes.get(name) || null,
  setAttribute: (name, value) => attributes.set(name, value),
  getElementById: (id) => id === 'legend' && hasLegend ? legendGroup : null
};
const state = {
  svgContent: ref('<svg/>'),
  svgContainer: ref({ querySelector: (selector) => selector === 'svg' ? svg : null }),
  mode: ref('circular'),
  generatedLegendPosition: ref('right'),
  circularBaseConfig: ref({
    viewBoxWidth: 1000,
    viewBoxHeight: 600,
    generatedViewBoxWidth: 1000,
    generatedViewBoxHeight: 600,
    legendWidth: 100,
    legendHeight: 120
  }),
  linearBaseConfig: ref({}),
  diagramElements: ref([diagramElement]),
  diagramElementOriginalTransforms: ref(new Map()),
  diagramElementBaseTransforms: ref(new Map([[diagramElement, { x: 10, y: 20 }]])),
  diagramOffset: { x: 4, y: -3 },
  legendInitialTransform: ref({ x: 1015, y: 240 }),
  legendCurrentOffset: { x: 7, y: 5 },
  plotTitleAutoTransform: ref({ x: 0, y: 0 }),
  pairwiseMatchFactors: ref({}),
  appliedPaletteColors: ref({}),
  legendColorOverrides: {},
  selectedResultIndex: ref(-1),
  results: ref([]),
  skipCaptureBaseConfig: ref(false),
  skipPositionReapply: ref(false),
  form: { legend: 'left' },
  adv: { plot_title_position: 'none' }
};
const repositionActions = createLegendRepositionActions({
  state,
  debugLog: () => {},
  legendActions: {
    getAllFeatureLegendGroups: () => [],
    getLegendLayoutLocalBounds: () => ({ x: 0, y: 0, width: 300, height: 120 }),
    recenterCurrentLegendRoot: () => {},
    reflowDualLegendLayout: () => {},
    reflowSingleLegendLayout: () => ({ legendWidth: 300, legendHeight: 120 })
  },
  svgActions: {
    ensureUniquePairwiseGradientIds: () => {},
    ensureUniqueSkewClipPathIds: () => {}
  },
  diagramActions: {
    applyDiagramShift: () => {},
    clearPlotTitleState: () => {},
    setPlotTitleAutoTransform: () => {}
  }
});

assert.equal(repositionActions.refreshLegendGeometry(), true);
const [, , viewBoxWidth, viewBoxHeight] = attributes.get('viewBox').split(/\s+/).map(Number);
const [, legendX, legendY] = legendAttributes.get('transform').match(
  /translate\(\s*([\d.-]+)\s*,\s*([\d.-]+)\s*\)/
).map(Number);
assert.equal(viewBoxHeight, 600);
assert.ok(viewBoxWidth > 1300);
assert.equal(state.generatedLegendPosition.value, 'right');
assert.deepEqual(state.legendCurrentOffset, { x: 7, y: 5 });
assert.deepEqual(state.diagramOffset, { x: 4, y: -3 });
assert.equal(legendX, state.legendInitialTransform.value.x + 7);
assert.equal(legendY, state.legendInitialTransform.value.y + 5);
assert.ok(legendX >= 0 && legendX + 300 <= viewBoxWidth);
assert.equal(diagramAttributes.get('transform'), 'translate(14, 17)');
assert.deepEqual(
  state.diagramElementOriginalTransforms.value.get(diagramElement),
  { x: 10, y: 20 }
);

const refreshedViewBox = attributes.get('viewBox');
state.generatedLegendPosition.value = 'none';
assert.equal(repositionActions.refreshLegendGeometry(), false);
assert.equal(attributes.get('viewBox'), refreshedViewBox);
state.generatedLegendPosition.value = 'right';
hasLegend = false;
assert.equal(repositionActions.refreshLegendGeometry(), false);
assert.equal(attributes.get('viewBox'), refreshedViewBox);

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
