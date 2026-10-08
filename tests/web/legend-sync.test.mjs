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
await cp(join(repoRoot, 'gbdraw', 'web', 'js', 'utils'), join(tempRoot, 'utils'), { recursive: true });
await cp(join(repoRoot, 'gbdraw', 'web', 'js', 'config.js'), join(tempRoot, 'config.js'));
await writeFile(join(tempRoot, 'package.json'), '{"type":"module"}\n', 'utf8');

const {
  SPECIFIC_COLOR_FILE_OWNER,
  buildLegendIntents,
  diffLegendIntents
} = await import(pathToFileURL(join(tempRoot, 'services', 'specific-color-rules.js')));
const {
  COMPARISON_LEGEND_SELECTOR,
  PAIRWISE_LEGEND_SELECTOR,
  getComparisonLegendGroup,
  getLegendChildById,
  parseTransformXY
} = await import(
  pathToFileURL(join(tempRoot, 'services', 'legend-svg.js'))
);
const { createLegendEntryActions } = await import(
  pathToFileURL(join(tempRoot, 'app', 'legend', 'entry-actions.js'))
);
const { createLegendStrokeActions } = await import(
  pathToFileURL(join(tempRoot, 'app', 'legend', 'stroke-actions.js'))
);
const { createLegendSortActions } = await import(
  pathToFileURL(join(tempRoot, 'app', 'legend', 'sort-actions.js'))
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
  /resetAllPositions[\s\S]+resetCompositionUserDeltas[\s\S]+commitActiveResultEdit\('layout-position-reset'\)/
);
assert.match(entryActionsSource, /setLegendGeometryChangedHandler/);
assert.ok((entryActionsSource.match(/onLegendGeometryChanged\(\);/g) || []).length >= 4);
assert.match(
  appSetupSource,
  /setLegendGeometryChangedHandler\(legendLayout\.refreshLegendGeometry\)/
);
assert.match(
  appSetupSource,
  /bindComposition\(context\)[\s\S]+captureBaseConfig\(\)/
);
assert.doesNotMatch(watchersSource, /captureBaseConfig/);
assert.match(configSource, /skipCaptureBaseConfig\.value = true;\s+applyResultsData/);
const sessionLegendSyncSource = appSetupSource.match(
  /adoptLegend\(context\)[\s\S]*?\n    bindComposition/
)?.[0] || '';
assert.match(
  sessionLegendSyncSource,
  /const drawn = !context\.bindingOptions\.isIncrementalEdit\s*\|\| Boolean\(context\.bindingOptions\.replaceGeneratedLegend\);[\s\S]+extractLegendEntries\(\{\s*replaceGeneratedInventory: !selecting && drawn,/
);
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

assert.deepEqual(buildLegendIntents([
  { feat: 'CDS', qual: 'gene', val: 'a', color: '#112233', cap: 'Historical [#112233]' },
  { feat: 'CDS', qual: 'gene', val: 'b', color: '#445566', cap: 'Historical [#445566]' }
]).intents, [
  { caption: 'Historical [#112233]', color: '#112233' },
  { caption: 'Historical [#445566]', color: '#445566' }
]);

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

  let dirtyMarks = 0;
  let layoutRefreshes = 0;
  const state = {
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
    commitActiveResultEdit: () => {
      dirtyMarks += 1;
      return true;
    }
  });
  // The layout owner lays the Legend out after a restore (zero shift).
  actions.setLegendGeometryChangedHandler(() => { layoutRefreshes += 1; });

  assert.equal(actions.reconcileLegendEntries(), false);
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
  // R3: History restore takes the Legend slots through orderLegendEntries, as
  // live Sort and Move do.
  assert.deepEqual(
    featureLegend.children.map((entry) => entry.querySelector('text').getAttribute('transform')),
    ['translate(22, 7)', 'translate(92, 7)']
  );
  assert.equal(featureLegend.children[0].querySelector('path').getAttribute('fill'), '#abcdef');
  assert.equal(dirtyMarks, 1);
  assert.equal(layoutRefreshes, 1);
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
      featureIds: entry.featureIds
    })),
    [
      {
        caption: 'Beta',
        color: '#abcdef',
        featureIds: ['feature-safe']
      },
      {
        caption: 'Gamma',
        color: '#778899',
        featureIds: []
      }
    ],
    'the sanitized mounted legend remains visual authority and mismatched metadata is ignored'
  );
  assert.equal(state.legendEntries.value.some((entry) => Object.hasOwn(entry, 'showStroke')), false,
    'the Stroke options disclosure is view state, not a Legend entry field (OV-157)');

  const capturedOwners = actions.captureLegendEntryOwners();
  featureLegend.children[0].setAttribute('data-legend-owner', 'specific-color-file');
  featureLegend.children[0].querySelector('path').setAttribute('fill', '#112233');
  assert.equal(actions.reconcileLegendEntries({ entryOwners: capturedOwners }), true);
  assert.equal(featureLegend.children[0].getAttribute('data-legend-owner'), null);
  assert.equal(featureLegend.children[0].querySelector('path').getAttribute('fill'), '#abcdef',
    'History restores the captured swatch, even when the original palette differs');
  assert.equal(state.legendEntries.value[0].color, '#abcdef');
  assert.equal(Object.hasOwn(state.legendEntries.value[0], 'owner'), false);

  const noOpDirtyMarks = dirtyMarks;
  assert.equal(actions.updateLegendEntryColor(0, '#abcdef'), false);
  assert.equal(actions.updateLegendEntryCaption(0, 'Beta'), false);
  assert.equal(dirtyMarks, noOpDirtyMarks);

  state.extractedFeatures = ref([]);
  state.featureStrokeOverrides = {};
  state.originalSvgStroke = ref({ color: null, width: null });
  const strokeActions = createLegendStrokeActions({
    state,
    commitActiveResultEdit: () => {
      dirtyMarks += 1;
      return true;
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
  betaSwatch.setAttribute('fill', '#aabbcc');
  assert.equal(actions.updateLegendEntryColorByCaption('Beta', '#abc', {commit:false}), false);
  assert.equal(betaSwatch.getAttribute('fill'), '#aabbcc');
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

  // The generated inventory follows the current diagram, independently of
  // explicitly owned rows and the default order of surviving categories.
  state.legendEntries.value = [];
  state.originalLegendOrder.value = ['Alpha', 'Beta'];
  state.deletedLegendEntries.value = [{ caption: 'Deleted', originalCaption: 'Deleted' }];
  state.originalLegendOrder.value.push('Deleted');
  featureLegend.children.forEach(entry => { entry.parentElement = null; });
  featureLegend.children = [];
  featureLegend.appendChild(mockLegendEntry('Gamma', '#334455', 0));
  featureLegend.appendChild(mockLegendEntry('Beta', '#445566', 70));
  const manual = mockLegendEntry('Manual', '#884422', 140);
  manual.setAttribute('data-legend-owner', 'direct-editor');
  featureLegend.appendChild(manual);
  actions.extractLegendEntries();
  assert.deepEqual(state.originalLegendOrder.value, ['Alpha', 'Beta', 'Deleted']);
  actions.extractLegendEntries({ replaceGeneratedInventory: true });
  assert.deepEqual(state.originalLegendOrder.value, ['Beta', 'Deleted', 'Gamma']);
  assert.deepEqual(state.legendEntries.value.map(e => e.caption), ['Gamma', 'Beta', 'Manual']);
  actions.extractLegendEntries();
  assert.deepEqual(state.originalLegendOrder.value, ['Beta', 'Deleted', 'Gamma']);
  featureLegend.children[0].remove();
  actions.extractLegendEntries();
  assert.deepEqual(state.originalLegendOrder.value, ['Beta', 'Deleted', 'Gamma']);
  actions.extractLegendEntries({ replaceGeneratedInventory: true });
  assert.deepEqual(state.originalLegendOrder.value, ['Beta', 'Deleted']);
  assert.deepEqual(state.legendEntries.value.map(e => e.caption), ['Beta', 'Manual']);
}

{
  // D-08 (PD-OI-063): without an order edit, Generate shows the renderer's
  // order, including a category it adds between existing ones, and the next
  // Generate replays no order. An edited order keeps being replayed.
  const { compileDirectEditorMutationPlan } = await import(
    pathToFileURL(join(tempRoot, 'app', 'candidate-render.js'))
  );
  const ref = (value) => ({ value });
  const svg = new MockElement('svg');
  const legend = new MockElement('g', { id: 'legend' });
  const featureLegend = new MockElement('g', { id: 'feature_legend' });
  legend.appendChild(featureLegend);
  svg.appendChild(legend);
  const render = (...captions) => {
    featureLegend.children.forEach(entry => { entry.parentElement = null; });
    featureLegend.children = [];
    captions.forEach((caption, index) => featureLegend.appendChild(mockLegendEntry(caption, '#112233', index * 70)));
  };
  const state = {
    results: ref([{ name: 'diagram.svg', content: 'unchanged' }]),
    selectedResultIndex: ref(0),
    svgContainer: ref({ querySelector: () => svg }),
    adv: {},
    legendEntries: ref([]),
    deletedLegendEntries: ref([]),
    originalLegendOrder: ref([]),
    originalLegendColors: ref({}),
    newLegendCaption: ref(''),
    newLegendColor: ref('#808080'),
    legendStrokeOverrides: {},
    legendColorOverrides: {},
    manualSpecificRules: [],
    skipCaptureBaseConfig: ref(false)
  };
  const actions = createLegendEntryActions({
    state,
    commitActiveResultEdit: () => true
  });
  const generate = (...rendered) => {
    const replayed = compileDirectEditorMutationPlan({
      catalogAdmission: { resultNames: ['diagram.svg'], renderedTargetsByOverrideKey: new Map(), resultIndexesByRenderedId: new Map() },
      legendEntries: state.legendEntries.value,
      originalLegendOrder: state.originalLegendOrder.value
    }).operationsByResult[0].legendOrder.map(({ captions }) => [...captions]);
    render(...rendered);
    actions.extractLegendEntries({ replaceGeneratedInventory: true });
    return replayed;
  };

  assert.deepEqual(generate('Core', 'Other'), []);
  assert.deepEqual(state.originalLegendOrder.value, ['Core', 'Other']);
  assert.deepEqual(generate('Core', 'Added', 'Other'), []);
  assert.deepEqual(state.originalLegendOrder.value, ['Core', 'Added', 'Other']);
  assert.deepEqual(generate('Core', 'Added', 'Other'), [], 'an unedited Legend order replays no order');

  // Sort Z-A, then Generate replays it; a new category follows the edited order.
  render('Other', 'Core', 'Added');
  actions.extractLegendEntries();
  assert.deepEqual(generate('Other', 'Core', 'Added', 'Late'), [['Other', 'Core', 'Added']]);
  assert.deepEqual(state.originalLegendOrder.value, ['Core', 'Added', 'Other', 'Late']);
  assert.deepEqual(generate('Other', 'Core', 'Added', 'Late'), [['Other', 'Core', 'Added', 'Late']]);

  // B20 (D-07): a displayed batch Result that may still show an earlier edited
  // order receives the default order; Generate still replays none.
  render('Core', 'Added', 'Other', 'Late');
  actions.extractLegendEntries();
  const compile = (options) => compileDirectEditorMutationPlan({
    catalogAdmission: { resultNames: ['diagram.svg'], renderedTargetsByOverrideKey: new Map(), resultIndexesByRenderedId: new Map() },
    legendEntries: state.legendEntries.value,
    originalLegendOrder: ['Core', 'Added', 'Other', 'Late'],
    ...options
  }).operationsByResult[0].legendOrder.map(({ captions }) => [...captions]);
  assert.deepEqual(compile({}), []);
  assert.deepEqual(compile({ replayDefaultLegendOrder: ['Core', 'Added', 'Other', 'Late'] }), [['Core', 'Added', 'Other', 'Late']]);
}

{
  // B19 (D-07, D-08): a History step made on another batch Result projects its
  // shared Legend intent onto the displayed Result, which never gains an entry
  // only the other Result draws and keeps its own entries.
  const ref = (value) => ({ value });
  const svg = new MockElement('svg');
  const legend = new MockElement('g', { id: 'legend' });
  const featureLegend = new MockElement('g', { id: 'feature_legend' });
  legend.appendChild(featureLegend);
  svg.appendChild(legend);
  const render = (...captions) => {
    featureLegend.children.forEach(entry => { entry.parentElement = null; });
    featureLegend.children = [];
    captions.forEach((caption, index) => featureLegend.appendChild(mockLegendEntry(caption, '#112233', index * 70)));
  };
  const drawn = () => featureLegend.children.map((entry) => entry.getAttribute('data-legend-key'));
  const listed = () => state.legendEntries.value.map((entry) => entry.caption);
  const list = (...captions) => captions.map((caption) => ({ caption, originalCaption: caption, color: '#112233' }));
  const state = {
    results: ref([{ name: 'r1.svg', content: 'unchanged' }, { name: 'r2.svg', content: 'unchanged' }]),
    selectedResultIndex: ref(0),
    svgContainer: ref({ querySelector: () => svg }),
    adv: {},
    legendEntries: ref([]),
    deletedLegendEntries: ref([]),
    originalLegendOrder: ref(['Alpha', 'Beta']),
    originalLegendColors: ref({}),
    newLegendCaption: ref(''),
    newLegendColor: ref('#808080'),
    legendStrokeOverrides: {},
    legendColorOverrides: {},
    manualSpecificRules: [],
    skipCaptureBaseConfig: ref(false)
  };
  let commits = 0;
  const actions = createLegendEntryActions({
    state,
    commitActiveResultEdit: () => { commits += 1; return true; },
    readActiveResultIdentity: () => 'result-1'
  });
  // Result 1 draws Own1; Result 2 drew Only2 and was sorted Z-A, which Result 1
  // shows. Undo restores Result 2's list from before the sort.
  render('Beta', 'Alpha', 'Own1');
  state.legendEntries.value = list('Alpha', 'Beta', 'Only2');
  assert.equal(actions.reconcileLegendEntries({ from: list('Only2', 'Beta', 'Alpha') }), true);
  assert.deepEqual(drawn(), ['Alpha', 'Beta', 'Own1'], 'Undo orders Result 1 without the Result 2 entry');
  assert.deepEqual(listed(), ['Alpha', 'Beta', 'Own1'], 'the Legend panel lists Result 1');
  assert.equal(commits, 1);

  // Redo of a deletion made on Result 2 removes the entry here too, and its
  // Undo returns this Result's entry; Only2 never appears.
  state.legendEntries.value = list('Alpha', 'Only2');
  actions.reconcileLegendEntries({ from: list('Alpha', 'Beta', 'Only2') });
  assert.deepEqual(drawn(), ['Alpha', 'Own1']);
  state.legendEntries.value = list('Alpha', 'Beta', 'Only2');
  actions.reconcileLegendEntries({ from: list('Alpha', 'Only2') });
  assert.deepEqual(drawn(), ['Alpha', 'Beta', 'Own1'], 'Undo returns the deleted entry to Result 1');
  assert.deepEqual(listed(), ['Alpha', 'Beta', 'Own1']);

  // A list that describes the displayed Result installs as is.
  state.legendEntries.value = list('Own1', 'Beta', 'Alpha');
  actions.reconcileLegendEntries({ from: list('Alpha', 'Beta', 'Own1') });
  assert.deepEqual(drawn(), ['Own1', 'Beta', 'Alpha']);

  // Live Sort orders the mounted Legend through the entry owner's port, the
  // ordering History restore uses (R3).
  let sortCommits = 0;
  const sortActions = createLegendSortActions({
    state,
    extractLegendEntries: actions.extractLegendEntries,
    orderMountedLegend: actions.orderMountedLegend,
    commitActiveResultEdit: (reason) => { if (reason === 'legend-order') sortCommits += 1; }
  });
  actions.extractLegendEntries();
  sortActions.sortLegendEntries('asc');
  assert.deepEqual(drawn(), ['Alpha', 'Beta', 'Own1']);
  assert.equal(sortCommits, 1);
  sortActions.sortLegendEntries('asc');
  assert.equal(sortCommits, 1, 'an order already shown records nothing');

  // Outside a batch the restored list describes the displayed Result, so it
  // installs as is whatever the step's other side lists (B19).
  state.results.value = [{ name: 'r1.svg', content: 'unchanged' }];
  state.legendEntries.value = list('Beta', 'Alpha', 'Only2');
  actions.reconcileLegendEntries({ from: list('Only2', 'Alpha', 'Beta') });
  assert.deepEqual(drawn(), ['Beta', 'Alpha', 'Only2']);

  // Without a mounted Legend there is nothing to order.
  const mounted = state.svgContainer.value;
  state.svgContainer.value = null;
  assert.equal(actions.orderMountedLegend(['Alpha']), null);
  sortActions.sortLegendEntries('desc');
  assert.equal(sortCommits, 1);
  state.svgContainer.value = mounted;
}

{
  // OV-47: each batch Result has its own generated Legend order. A draw made
  // while Result 2 is displayed does not replace the order of Result 1, which
  // is read from its Legend when it is displayed first, and each Result gets
  // its own order back when it is displayed again.
  const ref = (value) => ({ value });
  const svg = new MockElement('svg');
  const legend = new MockElement('g', { id: 'legend' });
  const featureLegend = new MockElement('g', { id: 'feature_legend' });
  legend.appendChild(featureLegend);
  svg.appendChild(legend);
  const render = (...captions) => {
    featureLegend.children.forEach(entry => { entry.parentElement = null; });
    featureLegend.children = [];
    captions.forEach((caption, index) => featureLegend.appendChild(mockLegendEntry(caption, '#112233', index * 70)));
  };
  const state = {
    results: ref([{ name: 'r1.svg', content: 'unchanged' }, { name: 'r2.svg', content: 'unchanged' }]),
    selectedResultIndex: ref(0),
    svgContainer: ref({ querySelector: () => svg }),
    adv: {},
    legendEntries: ref([]),
    deletedLegendEntries: ref([]),
    originalLegendOrder: ref([]),
    originalLegendColors: ref({}),
    newLegendCaption: ref(''),
    newLegendColor: ref('#808080'),
    legendStrokeOverrides: {},
    legendColorOverrides: {},
    manualSpecificRules: [],
    skipCaptureBaseConfig: ref(false)
  };
  let identity = 'result-2';
  const live = ['result-1', 'result-2'];
  const actions = createLegendEntryActions({
    state,
    commitActiveResultEdit: () => true,
    readActiveResultIdentity: () => identity
  });
  // Generate draws while Result 2 is displayed: it keeps Result 2's order.
  render('Rule', 'Other', 'tRNA', 'GC');
  actions.extractLegendEntries({ replaceGeneratedInventory: true, liveResultIdentities: live });
  assert.deepEqual(state.originalLegendOrder.value, ['Rule', 'Other', 'tRNA', 'GC']);
  // Result 1 is displayed first: its order is read before any projection, and
  // the inventory of Result 2 is not its own.
  render('CDS', 'tRNA', 'GC');
  identity = 'result-1';
  assert.deepEqual(actions.captureResultInventory(svg, { resultIdentity: 'result-1', liveResultIdentities: live }), ['CDS', 'tRNA', 'GC']);
  actions.adoptResultInventory('result-1');
  assert.deepEqual(state.originalLegendOrder.value, ['CDS', 'tRNA', 'GC']);
  // Result 2 comes back with its own inventory, not Result 1's.
  identity = 'result-2';
  assert.deepEqual(actions.captureResultInventory(svg, { resultIdentity: 'result-2', liveResultIdentities: live }), ['Rule', 'Other', 'tRNA', 'GC']);
  actions.adoptResultInventory('result-2');
  assert.deepEqual(state.originalLegendOrder.value, ['Rule', 'Other', 'tRNA', 'GC']);
  // A Result that no longer exists is forgotten.
  actions.captureResultInventory(svg, { resultIdentity: 'result-2', liveResultIdentities: ['result-2'] });
  identity = 'result-1';
  render('Z', 'CDS');
  assert.deepEqual(actions.captureResultInventory(svg, { resultIdentity: 'result-1', liveResultIdentities: ['result-1', 'result-2'] }), ['Z', 'CDS']);
  // A draw that replays an edited order shows it on the other Results, whose
  // order is then the displayed one with their own entries after it.
  state.originalLegendOrder.value = ['A', 'B'];
  render('B', 'A');
  actions.extractLegendEntries();
  identity = 'result-2';
  actions.extractLegendEntries({ replaceGeneratedInventory: true, liveResultIdentities: ['result-1', 'result-2', 'result-3'] });
  render('B', 'A', 'Own');
  assert.deepEqual(
    actions.captureResultInventory(svg, { resultIdentity: 'result-3', liveResultIdentities: ['result-1', 'result-2', 'result-3'] }),
    ['A', 'B', 'Own'],
    'a Result drawn with a replayed edited order keeps the displayed default order'
  );
}
