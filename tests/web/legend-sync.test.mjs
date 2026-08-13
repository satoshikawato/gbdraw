import assert from 'node:assert/strict';
import { cp, mkdtemp, readFile, writeFile } from 'node:fs/promises';
import { tmpdir } from 'node:os';
import { join } from 'node:path';
import { pathToFileURL } from 'node:url';

const repoRoot = process.cwd();
globalThis.CSS = { escape: (value) => String(value) };
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
const { createLegendEntryActions } = await import(
  pathToFileURL(join(tempRoot, 'app', 'legend', 'entry-actions.js'))
);
const { createLegendStrokeActions } = await import(
  pathToFileURL(join(tempRoot, 'app', 'legend', 'stroke-actions.js'))
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
const sessionLegendSyncSource = appSetupSource.match(
  /const synchronizeLoadedSessionLegendEntries[\s\S]*?\n  const importSession/
)?.[0] || '';
assert.match(sessionLegendSyncSource, /extractLegendEntries\(\)/);
assert.doesNotMatch(sessionLegendSyncSource, /initPyodide|addLegendEntry|removeLegendEntry/);
assert.doesNotMatch(appSetupSource, /restoreLoadedSessionLegendEntries/);
assert.match(configSource, /entries: normalizeSessionLegendEntries\(legend\.entries\)/);
assert.match(configSource, /deletedEntries: normalizeSessionLegendEntries\(legend\.deletedEntries\)/);

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

class MockElement {
  constructor(tagName, attributes = {}, textContent = '') {
    this.tagName = tagName;
    this.attributes = new Map(Object.entries(attributes));
    this.children = [];
    this.parentElement = null;
    this.textContent = textContent;
  }

  get id() { return this.getAttribute('id') || ''; }
  getAttribute(name) { return this.attributes.has(name) ? this.attributes.get(name) : null; }
  hasAttribute(name) { return this.attributes.has(name); }
  setAttribute(name, value) { this.attributes.set(name, String(value)); }
  removeAttribute(name) { this.attributes.delete(name); }
  appendChild(child) {
    if (child.parentElement) {
      child.parentElement.children = child.parentElement.children.filter((entry) => entry !== child);
    }
    child.parentElement = this;
    this.children.push(child);
    return child;
  }
  remove() {
    if (!this.parentElement) return;
    this.parentElement.children = this.parentElement.children.filter((entry) => entry !== this);
    this.parentElement = null;
  }
  cloneNode(deep = false) {
    const clone = new MockElement(this.tagName, Object.fromEntries(this.attributes), this.textContent);
    if (deep) this.children.forEach((child) => clone.appendChild(child.cloneNode(true)));
    return clone;
  }
  matchesSelector(selector) {
    if (selector === 'text' || selector === 'path') return this.tagName === selector;
    if (selector === '[transform]') return this.hasAttribute('transform');
    if (selector === 'g[data-legend-key]') {
      return this.tagName === 'g' && this.hasAttribute('data-legend-key');
    }
    if (selector.startsWith('#')) return this.id === selector.slice(1);
    return false;
  }
  querySelectorAll(selector) {
    const matches = [];
    const visit = (node) => {
      node.children.forEach((child) => {
        if (child.matchesSelector(selector)) matches.push(child);
        visit(child);
      });
    };
    visit(this);
    return matches;
  }
  querySelector(selector) { return this.querySelectorAll(selector)[0] || null; }
  getElementById(id) { return this.id === id ? this : this.querySelector(`#${id}`); }
}

const mockLegendEntry = (caption, color, x) => {
  const group = new MockElement('g', { 'data-legend-key': caption });
  group.appendChild(new MockElement('path', { fill: color, transform: `translate(${x},7)` }));
  group.appendChild(new MockElement('text', { transform: `translate(${x + 22},7)` }, caption));
  return group;
};

{
  const ref = (value) => ({ value });
  const svg = new MockElement('svg');
  const legend = new MockElement('g', { id: 'legend' });
  const featureLegend = new MockElement('g', { id: 'feature_legend' });
  featureLegend.appendChild(mockLegendEntry('Alpha', '#112233', 0));
  featureLegend.appendChild(mockLegendEntry('Beta', '#445566', 70));
  legend.appendChild(featureLegend);
  svg.appendChild(legend);

  let pyodideInitializations = 0;
  let dirtyMarks = 0;
  let layoutRefreshes = 0;
  const state = {
    pyodideReady: ref(false),
    results: ref([{ name: 'diagram.svg', content: 'unchanged' }]),
    selectedResultIndex: ref(0),
    svgContainer: ref({ querySelector: () => svg }),
    adv: {},
    legendEntries: ref([
      { caption: 'Alpha', color: '#112233' },
      { caption: 'Beta', color: '#445566' }
    ]),
    deletedLegendEntries: ref([]),
    originalLegendOrder: ref(['Alpha', 'Beta']),
    originalLegendColors: ref({ Alpha: '#112233', Beta: '#445566' }),
    newLegendCaption: ref(''),
    newLegendColor: ref('#808080'),
    legendStrokeOverrides: {},
    legendColorOverrides: {},
    manualSpecificRules: [],
    skipCaptureBaseConfig: ref(false)
  };
  const actions = createLegendEntryActions({
    state,
    getPyodide: () => null,
    ensurePyodide: async () => { pyodideInitializations += 1; },
    layoutActions: {
      compactLegendEntries: () => {},
      reflowDualLegendLayout: () => { layoutRefreshes += 1; },
      updatePairwiseLegendPositions: () => { layoutRefreshes += 1; }
    },
    previewRuntime: {
      applyLegendChanges: () => {
        dirtyMarks += 1;
        return true;
      }
    }
  });

  assert.equal(actions.reconcileLegendEntries(), false);
  assert.equal(pyodideInitializations, 0);
  assert.equal(dirtyMarks, 0);

  state.legendEntries.value = [
    { caption: 'Beta', color: '#abcdef' },
    { caption: 'Gamma', color: '#778899' }
  ];
  assert.equal(actions.reconcileLegendEntries(), true);
  assert.deepEqual(
    featureLegend.children.map((entry) => entry.getAttribute('data-legend-key')),
    ['Beta', 'Gamma']
  );
  assert.equal(featureLegend.children[0].querySelector('path').getAttribute('fill'), '#abcdef');
  assert.equal(dirtyMarks, 1);
  assert.equal(layoutRefreshes, 1);
  assert.equal(pyodideInitializations, 0);
  assert.equal(state.results.value[0].content, 'unchanged');

  state.legendEntries.value = [{ caption: 'Beta', color: '#abcdef' }];
  assert.equal(actions.reconcileLegendEntries(), true);
  state.legendEntries.value.push({ caption: 'Gamma', color: '#778899' });
  assert.equal(actions.reconcileLegendEntries(), true);
  assert.deepEqual(
    featureLegend.children.map((entry) => entry.getAttribute('data-legend-key')),
    ['Beta', 'Gamma']
  );
  assert.equal(dirtyMarks, 3);
  assert.equal(pyodideInitializations, 0);

  state.legendEntries.value = [
    {
      caption: 'Beta',
      color: '#abcdef',
      showStroke: true,
      featureIds: ['feature-safe']
    },
    {
      caption: 'Gamma',
      color: 'url(javascript:unsafe)',
      showStroke: true,
      featureIds: ['feature-unsafe']
    }
  ];
  actions.extractLegendEntries();
  assert.deepEqual(
    state.legendEntries.value.map((entry) => ({
      caption: entry.caption,
      color: entry.color,
      showStroke: entry.showStroke,
      featureIds: entry.featureIds
    })),
    [
      {
        caption: 'Beta',
        color: '#abcdef',
        showStroke: true,
        featureIds: ['feature-safe']
      },
      {
        caption: 'Gamma',
        color: '#778899',
        showStroke: false,
        featureIds: []
      }
    ],
    'the sanitized mounted legend remains visual authority and mismatched metadata is ignored'
  );
  assert.equal(pyodideInitializations, 0);

  const noOpDirtyMarks = dirtyMarks;
  assert.equal(actions.updateLegendEntryColor(0, '#abcdef'), false);
  assert.equal(actions.updateLegendEntryCaption(0, 'Beta'), false);
  assert.equal(dirtyMarks, noOpDirtyMarks);

  state.extractedFeatures = ref([]);
  state.featureStrokeOverrides = {};
  state.originalSvgStroke = ref({ color: null, width: null });
  const strokeActions = createLegendStrokeActions({
    state,
    previewRuntime: {
      markActiveResultDirty: () => {
        dirtyMarks += 1;
        return true;
      }
    }
  });
  assert.equal(strokeActions.updateLegendEntryStrokeColor(0, '#222222'), true);
  const strokeColorDirtyMarks = dirtyMarks;
  assert.equal(strokeActions.updateLegendEntryStrokeColor(0, '#222222'), false);
  assert.equal(dirtyMarks, strokeColorDirtyMarks);
  assert.equal(strokeActions.updateLegendEntryStrokeWidth(0, 2), true);
  const strokeWidthDirtyMarks = dirtyMarks;
  assert.equal(strokeActions.updateLegendEntryStrokeWidth(0, 2), false);
  assert.equal(dirtyMarks, strokeWidthDirtyMarks);

  const betaSwatch = featureLegend.children[0].querySelector('path');
  betaSwatch.setAttribute('stroke', '#222222');
  betaSwatch.setAttribute('stroke-width', '2');
  assert.equal(strokeActions.resetLegendEntryStroke(0), true);
  assert.equal(betaSwatch.getAttribute('stroke'), null);
  assert.equal(betaSwatch.getAttribute('stroke-width'), null);
  assert.equal(dirtyMarks, strokeWidthDirtyMarks + 1);
  const resetDirtyMarks = dirtyMarks;
  assert.equal(strokeActions.resetLegendEntryStroke(0), false);
  assert.equal(dirtyMarks, resetDirtyMarks);

  state.legendStrokeOverrides.Beta = { strokeColor: '#222222', strokeWidth: 2 };
  assert.equal(strokeActions.resetLegendEntryStroke(0), true);
  assert.deepEqual(state.legendStrokeOverrides, {});
  assert.equal(
    dirtyMarks,
    resetDirtyMarks,
    'removing a semantic override from an already-default SVG does not dirty the artifact'
  );
  assert.equal(strokeActions.resetAllStrokes(), false);
  assert.equal(dirtyMarks, resetDirtyMarks);
  state.featureStrokeOverrides['record:0:feature:1'] = { strokeColor: '#222222' };
  assert.equal(strokeActions.resetAllStrokes(), true);
  assert.deepEqual(state.featureStrokeOverrides, {});
  assert.equal(dirtyMarks, resetDirtyMarks);

  betaSwatch.setAttribute('stroke', 'gray');
  betaSwatch.setAttribute('stroke-width', '3');
  const appliedOverride = {
    originalStrokeColor: 'gray',
    originalStrokeWidth: 2,
    strokeWidth: 3
  };
  const historyReconcileDirtyMarks = dirtyMarks;
  assert.equal(strokeActions.reconcileStrokeOverrides({
    changes: [{
      path: ['editorState', 'legend', 'strokeOverrides', 'Beta'],
      before: undefined,
      after: appliedOverride
    }]
  }), true);
  assert.equal(betaSwatch.getAttribute('stroke'), 'gray');
  assert.equal(betaSwatch.getAttribute('stroke-width'), '2');
  assert.equal(dirtyMarks, historyReconcileDirtyMarks + 1);
  assert.equal(strokeActions.reconcileStrokeOverrides({
    changes: [{ path: ['editorState', 'legend', 'entries', '0', 'caption'] }]
  }), false);
}
