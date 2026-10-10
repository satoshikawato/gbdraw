import assert from 'node:assert/strict';
import { cp, mkdtemp, readFile, writeFile } from 'node:fs/promises';
import { tmpdir } from 'node:os';
import { join } from 'node:path';
import { pathToFileURL } from 'node:url';
import { withDrawings } from './helpers/drawing-state.mjs';

const repoRoot = process.cwd();
globalThis.CSS = { escape: (value) => String(value) };
const tempRoot = await mkdtemp(join(tmpdir(), 'gbdraw-legend-sync-'));
await cp(join(repoRoot, 'gbdraw', 'web', 'js', 'app'), join(tempRoot, 'app'), { recursive: true });
await cp(join(repoRoot, 'gbdraw', 'web', 'js', 'services'), join(tempRoot, 'services'), { recursive: true });
await cp(join(repoRoot, 'gbdraw', 'web', 'js', 'utils'), join(tempRoot, 'utils'), { recursive: true });
await cp(join(repoRoot, 'gbdraw', 'web', 'js', 'config.js'), join(tempRoot, 'config.js'));
await writeFile(join(tempRoot, 'package.json'), '{"type":"module"}\n', 'utf8');

const {
  buildLegendIntents,
  diffLegendIntents
} = await import(pathToFileURL(join(tempRoot, 'services', 'specific-color-rules.js')));
const {
  SPECIFIC_COLOR_FILE_OWNER,
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
const sortActionsSource = await readFile(join(tempRoot, 'app', 'legend', 'sort-actions.js'), 'utf8');
const colorActionsSource = await readFile(join(tempRoot, 'app', 'feature-editor', 'color-actions.js'), 'utf8');
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
// U3a A2a (R1): the Legend editor's writers edit the drawing's intent; only
// the Result executor writes Legend rows, shown through the root's port.
[['entry-actions', entryActionsSource], ['sort-actions', sortActionsSource], ['color-actions', colorActionsSource]]
  .forEach(([name, source]) => assert.doesNotMatch(
    source,
    /\.setAttribute\('data-legend-|\.appendChild\(|\.replaceWith\(|\borderLegendEntries\(|textContent = |\.remove\(\)/,
    `${name} writes no Legend row`
  ));
assert.match(appSetupSource, /deleteLegendEntry: editEditorIntent\('Delete legend item', deleteLegendEntry, \{ domains: \[\.\.\.LEGEND_STRUCTURE_DOMAINS, \.\.\.STROKE_DOMAINS\] \}\)/);
assert.match(appSetupSource, /const restoreLegendItems = \(label, restore\) => history\.runUndoable\(label,/);
assert.doesNotMatch(appSetupSource, /reconcileLegendEntries|prepareDisplayedResultLegend|hasRetiredResultLegend|legendChanged/);
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
assert.match(configSource, /const entries = normalizeSessionLegendEntries\(legend\.entries, 'entries', droppedLegendColors\)/);
assert.match(configSource, /entries: entries\.filter\(\(entry\) => entry\.dormant !== true\)/);
assert.match(configSource, /deletedEntries: normalizeSessionLegendEntries\(legend\.deletedEntries, 'deletedEntries', droppedLegendColors\)/);

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
    dormantLegendEntries: ref([]),
    originalLegendOrder: ref(['Alpha', 'Beta']),
    originalLegendColors: ref({ Alpha: '#112233', Beta: '#445566' }),
    legendStrokeOverrides: {},
    legendColorOverrides: {},
    manualSpecificRules: [],
    skipCaptureBaseConfig: ref(false)
  };
  const actions = createLegendEntryActions({
    state: withDrawings(state),
    commitActiveResultEdit: () => {
      dirtyMarks += 1;
      return true;
    },
    // The root's reader without a displayed Result: the row's recorded color.
    readShownLegendColor: (entry) => entry?.color
  });
  // The layout owner lays the Legend out after a restore (zero shift).
  actions.setLegendGeometryChangedHandler(() => { layoutRefreshes += 1; });

  // U3a A2a: Delete and Restore write the drawing's intent only; the
  // composition root shows it on the Result through the port in one compile.
  // Gamma and Gamma (1) are editor rows a Session keeps from the retired Add
  // legend item (R15-2): no generated row has their captions.
  const shownRows = () => JSON.stringify(featureLegend.children.map((entry) => [
    [...entry.attributes], entry.children.map((child) => [[...child.attributes], child.textContent])
  ]));
  const drawnRows = shownRows();
  const listedRows = () => state.legendEntries.value.map((entry) => entry.caption);
  state.legendEntries.value = [
    ...state.legendEntries.value,
    { caption: 'Gamma', originalCaption: 'Gamma', color: '#778899', featureIds: [] },
    { caption: 'Gamma (1)', originalCaption: 'Gamma (1)', color: '#000000', featureIds: [] }
  ];
  assert.deepEqual(listedRows(), ['Alpha', 'Beta', 'Gamma', 'Gamma (1)']);
  assert.equal(actions.deleteLegendEntry(3), true);
  assert.equal(actions.deleteLegendEntry(0), true);
  assert.deepEqual(listedRows(), ['Beta', 'Gamma']);
  assert.deepEqual(state.deletedLegendEntries.value.map((entry) => entry.caption), ['Gamma (1)', 'Alpha']);
  // Restore returns a row at its place in the default order and names the
  // rows of Python's it returns, which the root asks Python for when the
  // Result lacks them (O-2); an editor row is the editor's own.
  assert.deepEqual(await actions.restoreDeletedLegendEntries([1]), ['Alpha']);
  assert.deepEqual(listedRows(), ['Alpha', 'Beta', 'Gamma']);
  assert.deepEqual(await actions.restoreDeletedLegendEntries(), []);
  assert.deepEqual(listedRows(), ['Alpha', 'Beta', 'Gamma', 'Gamma (1)']);
  assert.deepEqual(state.deletedLegendEntries.value, []);
  assert.equal(await actions.restoreDeletedLegendEntries(), false);
  assert.equal(shownRows(), drawnRows, 'the writers leave the Result to the executor');
  assert.equal(dirtyMarks, 0);
  assert.equal(layoutRefreshes, 0);

  featureLegend.children.forEach((entry) => { entry.parentElement = null; });
  featureLegend.children = [];
  featureLegend.appendChild(mockLegendEntry('Beta', '#abcdef', 0));
  featureLegend.appendChild(mockLegendEntry('Gamma', '#778899', 70));
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
        featureIds: ['feature-unsafe']
      }
    ],
    'the paint is the sanitized Result\'s, and each row keeps the feature ids of its key (U3b)'
  );
  assert.equal(state.legendEntries.value.some((entry) => Object.hasOwn(entry, 'showStroke')), false,
    'the Stroke options disclosure is view state, not a Legend entry field (OV-157)');

  const noOpDirtyMarks = dirtyMarks;
  assert.equal(actions.updateLegendEntryColor(0, '#abcdef'), false);
  assert.equal(dirtyMarks, noOpDirtyMarks);

  // A Legend row stroke edit writes the intent only, from the stroke Python
  // drew on the swatch (its base record when the Result has one); the
  // composition root shows it through the executor (EU U2a).
  state.extractedFeatures = ref([]);
  state.featureStrokeOverrides = {};
  const strokeActions = createLegendStrokeActions({ state: withDrawings(state) });
  const betaSwatch = featureLegend.children[0].querySelector('path');
  betaSwatch.setAttribute('stroke', '#999999');
  betaSwatch.setAttribute('data-gbdraw-base-stroke', 'gray');
  betaSwatch.setAttribute('stroke-width', '2');
  const strokeDirtyMarks = dirtyMarks;
  assert.equal(strokeActions.updateLegendEntryStrokeColor(0, '#222222'), true);
  assert.equal(strokeActions.updateLegendEntryStrokeColor(0, '#222222'), false);
  assert.equal(strokeActions.updateLegendEntryStrokeWidth(0, 2), true);
  assert.equal(strokeActions.updateLegendEntryStrokeWidth(0, 2), false);
  assert.deepEqual(state.legendStrokeOverrides.Beta, {
    originalStrokeColor: 'gray', originalStrokeWidth: 2, strokeColor: '#222222', strokeWidth: 2
  });
  assert.equal(betaSwatch.getAttribute('stroke'), '#999999', 'the action leaves the Result to the executor');
  assert.equal(dirtyMarks, strokeDirtyMarks);
  assert.equal(strokeActions.resetLegendEntryStroke(0), true);
  assert.deepEqual(state.legendStrokeOverrides, {});
  assert.equal(strokeActions.resetLegendEntryStroke(0), false);
  assert.equal(strokeActions.resetAllStrokes(), false);
  state.featureStrokeOverrides['record:0:feature:1'] = { strokeColor: '#222222' };
  assert.equal(strokeActions.resetAllStrokes(), true);
  assert.deepEqual(state.featureStrokeOverrides, {});
  assert.equal(dirtyMarks, strokeDirtyMarks);
  betaSwatch.removeAttribute('data-gbdraw-base-stroke');

  // The generated inventory follows the current diagram, independently of
  // explicitly owned rows and the default order of surviving categories.
  // The intent lists the editor row.
  state.legendEntries.value = [{ caption: 'Manual', originalCaption: 'Manual', color: '#884422', featureIds: [] }];
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
    dormantLegendEntries: ref([]),
    originalLegendOrder: ref([]),
    originalLegendColors: ref({}),
    legendStrokeOverrides: {},
    legendColorOverrides: {},
    manualSpecificRules: [],
    skipCaptureBaseConfig: ref(false)
  };
  const actions = createLegendEntryActions({
    state: withDrawings(state),
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
  // B19 (D-07, D-08): a History step made on another batch Result gives the
  // displayed Result the step's shared Legend intent; it never gains an entry
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
    dormantLegendEntries: ref([]),
    originalLegendOrder: ref(['Alpha', 'Beta']),
    originalLegendColors: ref({}),
    legendStrokeOverrides: {},
    legendColorOverrides: {},
    manualSpecificRules: [],
    skipCaptureBaseConfig: ref(false)
  };
  let commits = 0;
  const actions = createLegendEntryActions({
    state: withDrawings(state),
    commitActiveResultEdit: () => { commits += 1; return true; },
    readActiveResultIdentity: () => 'result-1'
  });
  // U3a A2a: the Legend owner writes the list the displayed Result shows into
  // the intent, and the port compiles it; it writes no row.
  const shownRows = () => JSON.stringify(featureLegend.children.map((entry) => [...entry.attributes]));
  // Result 1 draws Own1; Result 2 drew Only2 and was sorted Z-A, which Result 1
  // shows. Undo restores Result 2's list from before the sort.
  state.originalLegendOrder.value = ['Alpha', 'Beta', 'Own1'];
  actions.rememberResultInventory('result-2', ['Alpha', 'Beta', 'Only2']);
  render('Beta', 'Alpha', 'Own1');
  let drawnRows = shownRows();
  state.legendEntries.value = list('Alpha', 'Beta', 'Only2');
  assert.equal(actions.adoptRestoredLegend({ from: list('Only2', 'Beta', 'Alpha') }), true);
  assert.deepEqual(listed(), ['Alpha', 'Beta', 'Own1'], 'Undo lists Result 1 in its order, without the Result 2 entry');
  assert.equal(shownRows(), drawnRows);

  // Redo of a deletion made on Result 2 removes the entry here too, and its
  // Undo returns this Result's entry; Only2 never appears.
  render('Alpha', 'Beta', 'Own1');
  state.legendEntries.value = list('Alpha', 'Only2');
  actions.adoptRestoredLegend({ from: list('Alpha', 'Beta', 'Only2') });
  assert.deepEqual(listed(), ['Alpha', 'Own1']);
  featureLegend.children[1].setAttribute('display', 'none');
  state.legendEntries.value = list('Alpha', 'Beta', 'Only2');
  actions.adoptRestoredLegend({ from: list('Alpha', 'Only2') });
  assert.deepEqual(listed(), ['Alpha', 'Beta', 'Own1'], 'Undo returns the deleted entry to Result 1');
  featureLegend.children[1].removeAttribute('display');

  // A step that renames and recolors a shared entry renames and recolors it
  // here; an editor row the step adds is listed, after this Result's own.
  state.legendEntries.value = [
    { caption: 'A2', originalCaption: 'Alpha', color: '#445566' }, ...list('Beta', 'Only2', 'Mine')
  ];
  actions.adoptRestoredLegend({ from: list('Alpha', 'Beta', 'Only2') });
  assert.deepEqual(state.legendEntries.value.map(({ caption, originalCaption, color }) => [caption, originalCaption, color]), [
    ['A2', 'Alpha', '#445566'], ['Beta', 'Beta', '#112233'], ['Own1', 'Own1', '#112233'], ['Mine', 'Mine', '#112233']
  ]);

  // A list that describes the displayed Result installs as is.
  state.legendEntries.value = list('Own1', 'Beta', 'Alpha');
  assert.equal(actions.adoptRestoredLegend({ from: list('Alpha', 'Beta', 'Own1') }), false);
  assert.deepEqual(listed(), ['Own1', 'Beta', 'Alpha']);

  // Live Sort and Move write the order into the intent and show it through
  // the port once (R3); an order already listed shows nothing.
  let shows = 0;
  const sortActions = createLegendSortActions({
    state: withDrawings(state),
    showLegendStructure: () => { shows += 1; }
  });
  drawnRows = shownRows();
  sortActions.sortLegendEntries('asc');
  assert.deepEqual(listed(), ['Alpha', 'Beta', 'Own1']);
  assert.equal(shows, 1);
  sortActions.sortLegendEntries('asc');
  assert.equal(shows, 1, 'an order already listed records nothing');
  sortActions.moveLegendEntryDown(0);
  assert.deepEqual(listed(), ['Beta', 'Alpha', 'Own1']);
  assert.equal(shows, 2);
  sortActions.sortLegendEntriesByDefault();
  assert.deepEqual(listed(), ['Alpha', 'Beta', 'Own1']);
  assert.equal(shows, 3);
  assert.equal(shownRows(), drawnRows);

  // Outside a batch the restored list describes the displayed Result, so it
  // installs as is whatever the step's other side lists (B19).
  state.results.value = [{ name: 'r1.svg', content: 'unchanged' }];
  state.legendEntries.value = list('Beta', 'Alpha', 'Only2');
  assert.equal(actions.adoptRestoredLegend({ from: list('Only2', 'Alpha', 'Beta') }), false);
  assert.deepEqual(listed(), ['Beta', 'Alpha', 'Only2']);
  assert.equal(commits, 0);
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
    dormantLegendEntries: ref([]),
    originalLegendOrder: ref([]),
    originalLegendColors: ref({}),
    legendStrokeOverrides: {},
    legendColorOverrides: {},
    manualSpecificRules: [],
    skipCaptureBaseConfig: ref(false)
  };
  let identity = 'result-2';
  const live = ['result-1', 'result-2'];
  const actions = createLegendEntryActions({
    state: withDrawings(state),
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
  // OV-348: a Result arriving from the other diagram mode brings the inventory
  // it left with (a live Legend rename changed it in place), which replaces
  // the copy stored when it was last displayed.
  actions.rememberResultInventory('result-3', ['A', 'Coding', 'Own']);
  actions.adoptResultInventory('result-3');
  assert.deepEqual(state.originalLegendOrder.value, ['A', 'Coding', 'Own']);
}

{
  // U3a (L1): a Result whose rows the executor reordered keeps Python's order
  // as a record, which its inventory reads; a row a Legend delete hid is not
  // listed.
  const ref = (value) => ({ value });
  const svg = new MockElement('svg');
  const legend = new MockElement('g', { id: 'legend' });
  const featureLegend = new MockElement('g', { id: 'feature_legend' });
  legend.appendChild(featureLegend);
  svg.appendChild(legend);
  ['tRNA', 'CDS', 'GC'].forEach((caption, index) => featureLegend.appendChild(mockLegendEntry(caption, '#112233', index * 70)));
  featureLegend.setAttribute('data-gbdraw-base-legend-order', JSON.stringify(['CDS', 'tRNA', 'GC']));
  const state = {
    results: ref([{ name: 'r1.svg', content: 'unchanged' }]),
    selectedResultIndex: ref(0),
    svgContainer: ref({ querySelector: () => svg }),
    adv: {},
    legendEntries: ref([]),
    deletedLegendEntries: ref([]),
    dormantLegendEntries: ref([]),
    originalLegendOrder: ref(['CDS', 'tRNA', 'GC']),
    originalLegendColors: ref({}),
    legendStrokeOverrides: {},
    legendColorOverrides: {},
    manualSpecificRules: [],
    skipCaptureBaseConfig: ref(false)
  };
  const actions = createLegendEntryActions({
    state: withDrawings(state), commitActiveResultEdit: () => true, readActiveResultIdentity: () => 'result-1'
  });
  assert.deepEqual(actions.captureResultInventory(svg, { resultIdentity: 'result-1', liveResultIdentities: ['result-1'] }), ['CDS', 'tRNA', 'GC']);
  featureLegend.children[0].setAttribute('display', 'none');
  actions.extractLegendEntries();
  assert.deepEqual(state.legendEntries.value.map((entry) => entry.caption), ['CDS', 'GC']);
}

{
  // U3b B2: the Legend list is the intent's view of the Result's rows (one
  // direction). A row is Python's key (its record), listed unless the intent
  // deletes that key and named by the intent's rename of it; the Result gives
  // the rows' order and the paint the executor showed.
  const ref = (value) => ({ value });
  const svg = new MockElement('svg');
  const legend = new MockElement('g', { id: 'legend' });
  const featureLegend = new MockElement('g', { id: 'feature_legend' });
  legend.appendChild(featureLegend);
  svg.appendChild(legend);
  const row = (caption, color, index, attributes = {}) => {
    const entry = mockLegendEntry(caption, color, index * 70);
    Object.entries(attributes).forEach(([name, value]) => entry.setAttribute(name, value));
    featureLegend.appendChild(entry);
    return entry;
  };
  const state = {
    results: ref([{ name: 'r1.svg', content: 'unchanged' }]),
    selectedResultIndex: ref(0),
    svgContainer: ref({ querySelector: () => svg }),
    adv: {},
    legendEntries: ref([
      { caption: 'CDS', originalCaption: 'CDS', color: '#111111', featureIds: ['cds-1'] },
      { caption: 'Skew+', originalCaption: 'GC skew (+)', color: '#6dded3', featureIds: [] }
    ]),
    deletedLegendEntries: ref([]),
    dormantLegendEntries: ref([]),
    originalLegendOrder: ref(['CDS', 'GC skew (+)']),
    originalLegendColors: ref({}),
    legendStrokeOverrides: {},
    legendColorOverrides: {},
    manualSpecificRules: [],
    skipCaptureBaseConfig: ref(false)
  };
  const actions = createLegendEntryActions({
    state: withDrawings(state), commitActiveResultEdit: () => true, readActiveResultIdentity: () => 'result-1'
  });
  const listed = () => state.legendEntries.value.map(({ caption, originalCaption, color, featureIds }) => (
    [caption, originalCaption, color, featureIds]
  ));
  // OV-243: a palette change gives the renamed GC skew row another color at
  // Generate; the row stays the rename of Python's "GC skew (+)", so the next
  // Generate renames it again.
  row('CDS', '#222222', 0);
  row('Skew+', '#80b1d3', 1, { 'data-gbdraw-base-data-legend-key': 'GC skew (+)' });
  actions.extractLegendEntries({ replaceGeneratedInventory: true });
  assert.deepEqual(listed(), [['CDS', 'CDS', '#222222', ['cds-1']], ['Skew+', 'GC skew (+)', '#80b1d3', []]]);
  assert.deepEqual(state.originalLegendOrder.value, ['CDS', 'GC skew (+)']);

  // A Result that differs from the intent only in its DOM leaves the list to
  // the intent: a row renamed without the intent keeps its name, a row hidden
  // without a delete is listed, a row the intent deletes is not, nor an
  // editor row the intent does not list, nor a row a rule commit retired.
  featureLegend.children.forEach((entry) => { entry.parentElement = null; });
  featureLegend.children = [];
  state.legendEntries.value = [
    { caption: 'Alpha', originalCaption: 'Alpha', color: '#112233', featureIds: [] },
    { caption: 'Beta', originalCaption: 'Beta', color: '#445566', featureIds: [] },
    { caption: 'Gamma', originalCaption: 'Gamma', color: '#778899', featureIds: [] }
  ];
  state.deletedLegendEntries.value = [{ caption: 'Gamma', originalCaption: 'Gamma', color: '#778899' }];
  state.originalLegendOrder.value = ['Alpha', 'Beta', 'Gamma', 'Retired'];
  row('Stale', '#112233', 0, { 'data-gbdraw-base-data-legend-key': 'Alpha' });
  row('Beta', '#445566', 1, { display: 'none', 'data-gbdraw-base-display': '' });
  row('Gamma', '#778899', 2);
  row('Manual', '#884422', 3, { 'data-legend-owner': 'direct-editor' });
  row('Retired', '#aa0000', 4, { display: 'none' });
  state.legendEntries.value = [...state.legendEntries.value];
  actions.extractLegendEntries();
  assert.deepEqual(listed().map(([caption, originalCaption]) => [caption, originalCaption]), [['Alpha', 'Alpha'], ['Beta', 'Beta']]);

  const named = () => listed().map(([caption, originalCaption]) => [caption, originalCaption]);
  const redraw = (entries, { deleted = [], dormant = [], inventory }, ...shown) => {
    featureLegend.children.forEach((entry) => { entry.parentElement = null; });
    featureLegend.children = [];
    state.legendEntries.value = entries;
    state.deletedLegendEntries.value = deleted;
    state.dormantLegendEntries.value = dormant;
    state.originalLegendOrder.value = inventory;
    shown.forEach(([caption, color, attributes], index) => row(caption, color, index, attributes));
  };
  const cds = { caption: 'CDS', originalCaption: 'CDS', color: '#111111', featureIds: [] };
  const skewRenamed = { caption: 'GC skew (-)', originalCaption: 'GC skew (+)', color: '#6dded3', featureIds: [] };
  const skewMinusDeleted = { deleted: [{ caption: 'GC skew (-)', originalCaption: 'GC skew (-)', color: '#aa0000' }], inventory: ['CDS', 'GC skew (+)', 'GC skew (-)'] };
  // Review M1: Remove "GC skew (-)", rename "GC skew (+)" to "GC skew (-)",
  // Generate. Python's hidden "GC skew (-)" row has no key record; another
  // row of the Result is Python's "GC skew (+)", so the hidden row is not
  // that rename and stays deleted: the row is listed once.
  redraw([cds, skewRenamed], skewMinusDeleted,
    ['CDS', '#111111'],
    ['GC skew (-)', '#6dded3', { 'data-gbdraw-base-data-legend-key': 'GC skew (+)' }],
    ['GC skew (-)', '#aa0000', { display: 'none', 'data-gbdraw-base-display': '' }]);
  actions.extractLegendEntries({ replaceGeneratedInventory: true });
  assert.deepEqual(named(), [['CDS', 'CDS'], ['GC skew (-)', 'GC skew (+)']]);
  // A Result saved before U3a shows the rename without a record and lacks
  // the deleted row: the row is the rename of the intent entry that names it.
  redraw([cds, skewRenamed], skewMinusDeleted, ['CDS', '#111111'], ['GC skew (-)', '#6dded3']);
  actions.extractLegendEntries();
  assert.deepEqual(named(), [['CDS', 'CDS'], ['GC skew (-)', 'GC skew (+)']]);
  // Review L1: a listed entry that names "tRNA" unrenamed outranks a dormant
  // rename of it, as the compile drops that dormant entry; the list names the
  // row as Generate draws it.
  redraw([cds, { caption: 'tRNA', originalCaption: 'tRNA', color: '#222222', featureIds: [] }], {
    dormant: [{ caption: 'T', originalCaption: 'tRNA', color: '#222222', featureIds: [] }], inventory: ['CDS', 'tRNA']
  }, ['CDS', '#111111'], ['tRNA', '#222222']);
  actions.extractLegendEntries();
  assert.deepEqual(named(), [['CDS', 'CDS'], ['tRNA', 'tRNA']]);
}

{
  // U3a A2a (gaps 1 and 2): the Legend rows of a rule commit are intent. The
  // diff and the OV-158 placement are computed on the listed rows; `apply`
  // writes them and names the rows the displayed Result shows at once: a
  // rule's new row where the Result lacks it, before the row it replaces, and
  // the retired row of a removed or renamed rule, which Python drew.
  const ref = (value) => ({ value });
  const svg = new MockElement('svg');
  const legend = new MockElement('g', { id: 'legend' });
  const featureLegend = new MockElement('g', { id: 'feature_legend' });
  [['CDS', '#112233'], ['Old', '#aa0000'], ['tRNA', '#334455']]
    .forEach(([caption, color], index) => featureLegend.appendChild(mockLegendEntry(caption, color, index * 70)));
  legend.appendChild(featureLegend);
  svg.appendChild(legend);
  const shownRows = () => JSON.stringify(featureLegend.children.map((entry) => [
    [...entry.attributes], entry.children.map((child) => [[...child.attributes], child.textContent])
  ]));
  const state = {
    results: ref([{ name: 'r1.svg', content: 'unchanged' }]),
    selectedResultIndex: ref(0),
    svgContainer: ref({ querySelector: () => svg }),
    adv: {},
    legendEntries: ref([
      { caption: 'CDS', originalCaption: 'CDS', color: '#112233' },
      { caption: 'Old', originalCaption: 'Old', color: '#aa0000' },
      { caption: 'tRNA', originalCaption: 'tRNA', color: '#334455' }
    ]),
    deletedLegendEntries: ref([]),
    dormantLegendEntries: ref([]),
    originalLegendOrder: ref(['CDS', 'Old', 'tRNA']),
    originalLegendColors: ref({}),
    legendStrokeOverrides: {},
    legendColorOverrides: {},
    manualSpecificRules: [],
    skipCaptureBaseConfig: ref(false)
  };
  let commits = 0;
  const actions = createLegendEntryActions({
    state: withDrawings(state), commitActiveResultEdit: () => { commits += 1; return true; }, readActiveResultIdentity: () => 'result-1'
  });
  const listed = () => state.legendEntries.value.map(({ caption, color }) => [caption, color]);
  const drawnRows = shownRows();
  const renamed = await actions.prepareFileLegendEntries([{ caption: 'New', color: '#aa0000' }], {
    previousFileIntents: [{ caption: 'Old', color: '#aa0000' }],
    placement: { caption: 'New', at: 'Old' }
  });
  assert.deepEqual(renamed.diff, {
    add: [{ caption: 'New', color: '#aa0000' }], update: [], remove: [{ caption: 'Old', color: '#aa0000' }], unchanged: []
  });
  assert.equal(renamed.isCurrent(), true);
  assert.deepEqual(renamed.apply(), { add: [{ caption: 'New', color: '#aa0000', before: 'Old' }], retire: ['Old'] });
  assert.deepEqual(listed(), [['CDS', '#112233'], ['New', '#aa0000'], ['tRNA', '#334455']]);
  // A rule color change draws no row; a listed row the rules do not own keeps
  // its color or stops the commit.
  const recolored = await actions.prepareFileLegendEntries([{ caption: 'New', color: '#00aa00' }], {
    previousFileIntents: [{ caption: 'New', color: '#aa0000' }]
  });
  assert.deepEqual(recolored.diff.update, [{ caption: 'New', color: '#00aa00' }]);
  assert.equal(recolored.apply(), null);
  assert.deepEqual(listed()[1], ['New', '#00aa00']);
  await assert.rejects(
    actions.prepareFileLegendEntries([{ caption: 'tRNA', color: '#000000' }]),
    /Legend entry "tRNA" already exists with a different color/
  );
  const reused = await actions.prepareFileLegendEntries([{ caption: 'tRNA', color: '#334455' }]);
  assert.deepEqual(reused.diff.unchanged, [{ caption: 'tRNA', color: '#334455' }]);
  assert.equal(reused.apply(), null);
  assert.equal(shownRows(), drawnRows, 'a rule commit leaves the Result to the executor');
  assert.equal(commits, 0);
}
