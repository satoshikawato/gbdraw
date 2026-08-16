import assert from 'node:assert/strict';
import { mkdir, mkdtemp, readFile, writeFile } from 'node:fs/promises';
import { tmpdir } from 'node:os';
import { join } from 'node:path';
import { pathToFileURL } from 'node:url';
import { runInNewContext } from 'node:vm';

import { FakeSvgElement } from './helpers/editor-svg-fixture.mjs';

const repoRoot = process.cwd();
const sourcePath = join(repoRoot, 'gbdraw', 'web', 'js', 'app', 'preview-runtime.js');
const featureDomSourcePath = join(repoRoot, 'gbdraw', 'web', 'js', 'app', 'feature-dom.js');
const featureVisibilityDomSourcePath = join(
  repoRoot,
  'gbdraw',
  'web',
  'js',
  'app',
  'feature-visibility-dom.js'
);
const mutationJournalSourcePath = join(repoRoot, 'gbdraw', 'web', 'js', 'app', 'dom-mutation-journal.js');
const serviceNames = [
  'runtime-test-hooks.js',
  'session-feature-metadata.js',
  'svg-result-normalization.js',
  'svg-result-ingestion.js',
  'svg-sanitization.js',
  'svg-serialization.js'
];
const tempDir = await mkdtemp(join(tmpdir(), 'gbdraw-preview-runtime-'));
await writeFile(join(tempDir, 'package.json'), '{"type":"module"}\n', 'utf8');
await mkdir(join(tempDir, 'app'), { recursive: true });
await mkdir(join(tempDir, 'services'), { recursive: true });
await writeFile(join(tempDir, 'app', 'preview-runtime.js'), await readFile(sourcePath, 'utf8'), 'utf8');
await writeFile(join(tempDir, 'app', 'feature-dom.js'), await readFile(featureDomSourcePath, 'utf8'), 'utf8');
await writeFile(
  join(tempDir, 'app', 'feature-visibility-dom.js'),
  await readFile(featureVisibilityDomSourcePath, 'utf8'),
  'utf8'
);
await writeFile(
  join(tempDir, 'app', 'dom-mutation-journal.js'),
  await readFile(mutationJournalSourcePath, 'utf8'),
  'utf8'
);
await Promise.all(serviceNames.map(async (name) => {
  const servicePath = join(repoRoot, 'gbdraw', 'web', 'js', 'services', name);
  await writeFile(join(tempDir, 'services', name), await readFile(servicePath, 'utf8'), 'utf8');
}));

const { createPreviewRuntime } = await import(pathToFileURL(join(tempDir, 'app', 'preview-runtime.js')));
const { createDomMutationJournal } = await import(
  pathToFileURL(join(tempDir, 'app', 'dom-mutation-journal.js'))
);

const ref = (value) => ({ value });

class FakeFeatureElement {
  constructor(id, { featureId = id, part = 'block', fill = null } = {}) {
    this.id = id;
    this.attributes = {
      'data-gbdraw-feature-id': featureId,
      'data-gbdraw-feature-part': part
    };
    if (fill !== null) this.attributes.fill = fill;
  }

  getAttribute(name) {
    return Object.hasOwn(this.attributes, name) ? this.attributes[name] : null;
  }

  setAttribute(name, value) {
    this.attributes[name] = String(value);
  }

  removeAttribute(name) {
    delete this.attributes[name];
  }
}

const makeSvg = (elements) => ({
  elements,
  querySelectorAll: () => elements,
  getElementById: (id) => elements.find((element) => element.id === id) || null
});

const featureA = new FakeFeatureElement('feature-a');
const featureBBlock = new FakeFeatureElement('feature-b__part1', {
  featureId: 'feature-b',
  fill: '#111111'
});
const featureBConnector = new FakeFeatureElement('feature-b__line1', {
  featureId: 'feature-b',
  part: 'connector',
  fill: 'none'
});
const rendererHiddenFeature = new FakeFeatureElement('renderer-hidden-feature');
rendererHiddenFeature.setAttribute('display', 'none');
const svg = makeSvg([featureA, featureBConnector, featureBBlock, rendererHiddenFeature]);
let replacementCount = 0;
const resultItems = new Proxy([
  { name: 'one.svg', content: '<svg id="old-one"></svg>' },
  { name: 'two.svg', content: '<svg id="old-two"></svg>' }
], {
  set(target, property, value) {
    if (/^\d+$/.test(String(property))) replacementCount += 1;
    target[property] = value;
    return true;
  }
});
const state = {
  results: ref(resultItems),
  selectedResultIndex: ref(0),
  skipCaptureBaseConfig: ref(false),
  svgContainer: ref({
    querySelector: (selector) => (selector === 'svg' ? svg : null)
  })
};
let serializeCount = 0;
let failSerialization = false;
const runtime = createPreviewRuntime({
  state,
  serializeSvg: (targetSvg) => {
    if (failSerialization) throw new Error('serialization failed');
    serializeCount += 1;
    return `<svg data-count="${serializeCount}" data-elements="${targetSvg.elements.length}"></svg>`;
  }
});

runtime.mountResultSvg(0, svg);
assert.equal(runtime.applyFeatureVisibilityChanges([{
  featureId: 'renderer-hidden-feature',
  mode: 'exclude_matching'
}]), false);
assert.equal(rendererHiddenFeature.getAttribute('display'), 'none');
assert.equal(rendererHiddenFeature.getAttribute('data-gbdraw-feature-visibility-preview'), null);
const directCommit = runtime.commitDomEdit({
  reason: 'direct-contract',
  invalidateIndexes: ['legend'],
  mutate: ({ svg: targetSvg, resultIndex }) => {
    assert.equal(targetSvg, svg);
    assert.equal(resultIndex, 0);
    featureA.setAttribute('data-direct-edit', 'true');
    return 1;
  }
});
assert.deepEqual(directCommit, {
  changed: true,
  flushed: true,
  resultIndex: 0,
  reason: 'direct-contract'
});
assert.equal(serializeCount, 1);
assert.equal(replacementCount, 1);
assert.equal(state.skipCaptureBaseConfig.value, true);

state.skipCaptureBaseConfig.value = false;
const noOp = runtime.commitDomEdit({
  reason: 'no-op',
  mutate: () => 0
});
assert.deepEqual(noOp, {
  changed: false,
  flushed: false,
  resultIndex: 0,
  reason: 'no-op'
});
assert.equal(serializeCount, 1);
assert.equal(replacementCount, 1);
assert.equal(state.skipCaptureBaseConfig.value, false);

runtime.getActiveRuntime().indexes.legend = { cached: true };
assert.equal(runtime.applyFeatureVisibilityChanges([{ featureId: 'feature-a', mode: 'off' }]), true);
assert.equal(featureA.getAttribute('display'), 'none');
assert.equal(serializeCount, 2);
assert.equal(replacementCount, 2);
assert.equal(runtime.getActiveRuntime().indexes.features, null);
assert.deepEqual(runtime.getActiveRuntime().indexes.legend, { cached: true });
assert.equal(state.skipCaptureBaseConfig.value, true);

state.skipCaptureBaseConfig.value = false;
assert.equal(runtime.applyFeatureVisibilityChanges([{ featureId: 'feature-a', mode: 'off' }]), false);
assert.equal(runtime.applyFeatureFillChanges([{ featureId: 'feature-b', color: '#111111' }]), false);
featureBBlock.setAttribute('stroke', '#222222');
featureBBlock.setAttribute('stroke-width', '3');
featureBConnector.setAttribute('stroke', '#222222');
featureBConnector.setAttribute('stroke-width', '3');
assert.equal(runtime.applyFeatureStrokeChanges([{
  featureId: 'feature-b',
  strokeColor: '#222222',
  strokeWidth: 3
}]), false);
assert.equal(runtime.getActiveRuntime().dirty, false);
assert.equal(runtime.flushActiveResult(), false);
assert.equal(serializeCount, 2);
assert.equal(replacementCount, 2);
assert.equal(state.skipCaptureBaseConfig.value, false);
state.skipCaptureBaseConfig.value = false;
assert.equal(runtime.flushActiveResult(), false);
assert.equal(serializeCount, 2);

await runtime.runDomEdit({
  reason: 'compound-feature-edit',
  action: async () => {
    assert.equal(runtime.applyFeatureFillChanges([{ featureId: 'feature-b', color: '#abcdef' }]), true);
    assert.equal(runtime.applyFeatureVisibilityChanges([{ featureId: 'feature-b', mode: 'off' }]), true);
    await Promise.resolve();
    assert.equal(runtime.applyFeatureStrokeChanges([{
      featureId: 'feature-b',
      strokeColor: '#333333',
      strokeWidth: 4
    }]), true);
  }
});
assert.equal(featureBBlock.getAttribute('fill'), '#abcdef');
assert.equal(featureBConnector.getAttribute('fill'), 'none');
assert.equal(featureBBlock.getAttribute('display'), 'none');
assert.equal(featureBConnector.getAttribute('display'), 'none');
assert.equal(featureBBlock.getAttribute('stroke'), '#333333');
assert.equal(serializeCount, 4, 'the synchronous prefix and post-await continuation commit separately');
assert.equal(replacementCount, 4);

await assert.rejects(
  runtime.runDomEdit({
    reason: 'failed-compound-edit',
    action: () => {
      runtime.commitDomEdit({
        reason: 'partial-edit',
        mutate: ({ mutation }) => mutation.setAttribute(featureA, 'data-partial-edit', 'true')
      });
      throw new Error('action failed');
    }
  }),
  /action failed/
);
assert.equal(serializeCount, 4);
assert.equal(replacementCount, 4);
assert.equal(featureA.getAttribute('data-partial-edit'), null);
assert.equal(runtime.getActiveRuntime().dirty, false);
assert.equal(runtime.flushActiveResult(), false);

const featureAFillBeforeFailedBulk = featureA.getAttribute('fill');
const featureBSetAttribute = featureBBlock.setAttribute.bind(featureBBlock);
featureBBlock.setAttribute = (name, value) => {
  if (name === 'fill' && value === '#badbad') throw new Error('second target failed');
  featureBSetAttribute(name, value);
};
assert.throws(() => runtime.applyFeatureFillChanges([
  { featureId: 'feature-a', color: '#123456' },
  { featureId: 'feature-b', color: '#badbad' }
]), /second target failed/);
assert.equal(featureA.getAttribute('fill'), featureAFillBeforeFailedBulk);
assert.equal(featureBBlock.getAttribute('fill'), '#abcdef');
assert.equal(runtime.getActiveRuntime().dirty, false);
assert.equal(serializeCount, 4);
assert.equal(replacementCount, 4);
featureBBlock.setAttribute = featureBSetAttribute;

failSerialization = true;
runtime.getActiveRuntime().indexes.legend = { retained: true };
const retainedLegendIndex = runtime.getActiveRuntime().indexes.legend;
const canonicalOverrides = {};
let inheritedStateValue = 'before';
const inheritedState = Object.create(Object.defineProperty({}, 'value', {
  get: () => inheritedStateValue,
  set: (value) => {
    inheritedStateValue = value;
  }
}));
const resultBeforeFailedCommit = state.results.value[0];
state.skipCaptureBaseConfig.value = false;
assert.throws(() => runtime.commitDomEdit({
  reason: 'serialization-failure',
  invalidateIndexes: ['legend'],
  mutate: ({ mutation }) => {
    mutation.setAttribute(featureA, 'data-retry-edit', 'true');
    mutation.setProperty(canonicalOverrides, 'feature-a', '#abcdef');
    mutation.setProperty(inheritedState, 'value', 'after');
    return true;
  }
}), /serialization failed/);
assert.equal(featureA.getAttribute('data-retry-edit'), null);
assert.equal(Object.hasOwn(canonicalOverrides, 'feature-a'), false);
assert.equal(inheritedStateValue, 'before');
assert.equal(state.results.value[0], resultBeforeFailedCommit);
assert.equal(runtime.getActiveRuntime().dirty, false);
assert.deepEqual([...runtime.getActiveRuntime().dirtyReasons], []);
assert.equal(runtime.getActiveRuntime().indexes.legend, retainedLegendIndex);
assert.equal(state.skipCaptureBaseConfig.value, false);
assert.equal(replacementCount, 4);
const editorMetadata = {};
await assert.rejects(runtime.runDomEdit({
  reason: 'batched-serialization-failure',
  action: () => {
    runtime.commitDomEdit({
      reason: 'editor-metadata',
      journalChangesAffectResult: false,
      mutate: ({ mutation }) => {
        mutation.setProperty(editorMetadata, 'labelKey', 'label-1');
        return false;
      }
    });
    runtime.commitDomEdit({
      reason: 'persistent-label-edit',
      mutate: ({ mutation }) => mutation.setAttribute(featureA, 'data-label-edit', 'true')
    });
  }
}), /serialization failed/);
assert.deepEqual(editorMetadata, {});
assert.equal(featureA.getAttribute('data-label-edit'), null);
assert.equal(runtime.getActiveRuntime().dirty, false);
assert.equal(runtime.getActiveRuntime().indexes.legend, retainedLegendIndex);
failSerialization = false;
assert.equal(runtime.flushActiveResult(), false);
assert.equal(serializeCount, 4);
assert.equal(replacementCount, 4);

assert.throws(() => runtime.commitDomEdit({
  reason: 'async-mutation',
  mutate: async () => true
}), /must be synchronous/);

state.skipCaptureBaseConfig.value = false;
runtime.selectResult(1);
assert.equal(serializeCount, 4);
assert.equal(state.skipCaptureBaseConfig.value, false);
assert.equal(state.selectedResultIndex.value, 1);

const legacyConnector = new FakeFeatureElement('feature-c__line1', {
  featureId: '',
  part: '',
  fill: 'none'
});
delete legacyConnector.attributes['data-gbdraw-feature-id'];
delete legacyConnector.attributes['data-gbdraw-feature-part'];
const legacyBlock = new FakeFeatureElement('feature-c__part1', {
  featureId: '',
  part: '',
  fill: '#222222'
});
delete legacyBlock.attributes['data-gbdraw-feature-id'];
delete legacyBlock.attributes['data-gbdraw-feature-part'];
const legacySvg = makeSvg([legacyConnector, legacyBlock]);
state.svgContainer.value = { querySelector: (selector) => (selector === 'svg' ? legacySvg : null) };
runtime.mountResultSvg(1, legacySvg);
runtime.applyFeatureFillChanges([{ featureId: 'feature-c', color: '#fedcba' }]);
assert.equal(legacyBlock.getAttribute('fill'), '#fedcba');
assert.equal(legacyConnector.getAttribute('fill'), 'none');
assert.equal(serializeCount, 5);
assert.equal(replacementCount, 5);

const deferred = () => {
  let resolve;
  const promise = new Promise((resolvePromise) => {
    resolve = resolvePromise;
  });
  return { promise, resolve };
};

const createLeaseRaceFixture = () => {
  const oldFeature = new FakeFeatureElement('race-feature', { fill: '#111111' });
  const oldSvg = makeSvg([oldFeature]);
  const initialResult = { name: 'race.svg', content: '<svg data-owner="old"></svg>' };
  const raceState = {
    results: ref([initialResult, { name: 'second.svg', content: '<svg data-owner="second"></svg>' }]),
    selectedResultIndex: ref(0),
    resultGenerationKey: ref('generation-1'),
    skipCaptureBaseConfig: ref(false),
    svgContainer: ref({ querySelector: (selector) => selector === 'svg' ? oldSvg : null })
  };
  const raceRuntime = createPreviewRuntime({
    state: raceState,
    serializeSvg: (targetSvg) => `<svg data-fill="${targetSvg.elements[0]?.getAttribute('fill') || ''}"></svg>`
  });
  raceRuntime.mountResultSvg(0, oldSvg);
  return { initialResult, oldFeature, oldSvg, runtime: raceRuntime, state: raceState };
};

const raceCases = [
  {
    name: 'selected Result switch',
    replace: ({ runtime: raceRuntime, state: raceState }) => {
      const nextSvg = makeSvg([new FakeFeatureElement('race-feature', { fill: '#222222' })]);
      raceState.selectedResultIndex.value = 1;
      raceState.svgContainer.value = { querySelector: (selector) => selector === 'svg' ? nextSvg : null };
      raceRuntime.mountResultSvg(1, nextSvg);
      return raceState.results.value[1];
    }
  },
  {
    name: 'same-index Result replacement',
    replace: ({ state: raceState }) => {
      const replacement = { name: 'race.svg', content: '<svg data-owner="replacement"></svg>' };
      raceState.results.value[0] = replacement;
      return replacement;
    }
  },
  {
    name: 'Generate replacement',
    replace: ({ runtime: raceRuntime, state: raceState }) => {
      const replacement = { name: 'generated.svg', content: '<svg data-owner="generated"></svg>' };
      const nextSvg = makeSvg([new FakeFeatureElement('race-feature', { fill: '#333333' })]);
      raceState.resultGenerationKey.value = 'generation-2';
      raceState.results.value = [replacement];
      raceState.svgContainer.value = { querySelector: (selector) => selector === 'svg' ? nextSvg : null };
      raceRuntime.mountResultSvg(0, nextSvg);
      return replacement;
    }
  },
  {
    name: 'session Load document replacement',
    replace: ({ runtime: raceRuntime, state: raceState }) => {
      const replacement = { name: 'loaded.svg', content: '<svg data-owner="loaded"></svg>' };
      const nextSvg = makeSvg([new FakeFeatureElement('race-feature', { fill: '#444444' })]);
      raceState.results.value = [replacement];
      raceState.svgContainer.value = { querySelector: (selector) => selector === 'svg' ? nextSvg : null };
      raceRuntime.mountResultSvg(0, nextSvg);
      return replacement;
    }
  },
  {
    name: 'History artifact replacement',
    replace: ({ runtime: raceRuntime, state: raceState }) => {
      const replacement = { name: 'history.svg', content: '<svg data-owner="history"></svg>' };
      const nextSvg = makeSvg([new FakeFeatureElement('race-feature', { fill: '#555555' })]);
      raceState.results.value[0] = replacement;
      raceState.svgContainer.value = { querySelector: (selector) => selector === 'svg' ? nextSvg : null };
      raceRuntime.mountResultSvg(0, nextSvg);
      return replacement;
    }
  },
  {
    name: 'SVG remount',
    replace: ({ runtime: raceRuntime, state: raceState }) => {
      const retainedResult = raceState.results.value[0];
      const nextSvg = makeSvg([new FakeFeatureElement('race-feature', { fill: '#666666' })]);
      raceState.svgContainer.value = { querySelector: (selector) => selector === 'svg' ? nextSvg : null };
      raceRuntime.mountResultSvg(0, nextSvg);
      return retainedResult;
    }
  },
  {
    name: 'same SVG element with a new document admission',
    replace: ({ oldSvg, runtime: raceRuntime, state: raceState }) => {
      const retainedResult = raceState.results.value[0];
      raceRuntime.mountResultSvg(0, oldSvg);
      return retainedResult;
    }
  }
];

for (const raceCase of raceCases) {
  const fixture = createLeaseRaceFixture();
  const lease = fixture.runtime.beginDomEditLease({ reason: `race-${raceCase.name}` });
  assert.ok(lease);
  assert.equal(lease.mutate(({ mutation }) => (
    mutation.setAttribute(fixture.oldFeature, 'fill', '#abcdef')
  )), true);
  const gate = deferred();
  const settlement = (async () => {
    await gate.promise;
    return lease.commit();
  })();
  const currentOwner = raceCase.replace(fixture);
  const currentContent = currentOwner.content;
  gate.resolve();
  await assert.rejects(
    settlement,
    (error) => error?.name === 'StaleDomEditError',
    raceCase.name
  );
  assert.equal(currentOwner.content, currentContent, raceCase.name);
  assert.equal(fixture.oldFeature.getAttribute('fill'), '#111111', raceCase.name);
}

{
  const fixture = createLeaseRaceFixture();
  const lease = fixture.runtime.beginDomEditLease({ reason: 'editor-only-lease-metadata' });
  assert.equal(lease.mutate(({ mutation }) => {
    mutation.setAttribute(fixture.oldFeature, 'data-label-editable', 'true');
    return false;
  }, { journalChangesAffectResult: false }), true);
  assert.equal(lease.mutate(() => false), false, 'prior journal entries do not taint later no-ops');
  assert.equal(lease.commit(), false);
  assert.equal(fixture.oldFeature.getAttribute('data-label-editable'), 'true');
  assert.equal(fixture.state.results.value[0], fixture.initialResult);
  assert.equal(fixture.initialResult.content, '<svg data-owner="old"></svg>');
  assert.equal(fixture.runtime.getActiveRuntime().dirty, false);
}

{
  const fixture = createLeaseRaceFixture();
  const lease = fixture.runtime.beginDomEditLease({
    reason: 'exclusive-async-edit',
    invalidateIndexes: ['features']
  });
  assert.equal(lease.mutate(({ mutation }) => (
    mutation.setAttribute(fixture.oldFeature, 'fill', '#abcdef')
  )), true);
  const competingCanonical = { visibility: 'default' };

  assert.throws(() => fixture.runtime.runDomEditSync({
    reason: 'competing-live-edit',
    canonicalState: [{ target: competingCanonical }],
    action: () => {
      competingCanonical.visibility = 'off';
      return fixture.runtime.commitDomEdit({
        reason: 'competing-visibility',
        mutate: ({ mutation }) => mutation.setAttribute(fixture.oldFeature, 'display', 'none')
      });
    }
  }), (error) => error?.name === 'LiveEditBusyError');

  assert.deepEqual(competingCanonical, { visibility: 'default' });
  assert.equal(fixture.oldFeature.getAttribute('fill'), '#abcdef');
  assert.equal(fixture.oldFeature.getAttribute('display'), null);
  assert.equal(fixture.state.results.value[0], fixture.initialResult);
  assert.equal(fixture.initialResult.content, '<svg data-owner="old"></svg>');
  assert.equal(lease.cancel(), true);
  assert.equal(fixture.oldFeature.getAttribute('fill'), '#111111');
  assert.equal(fixture.runtime.getActiveRuntime().dirty, false);
}

{
  const fixture = createLeaseRaceFixture();
  const active = fixture.runtime.getActiveRuntime();
  const retainedLegendIndex = { retained: true };
  active.dirty = true;
  active.dirtyReasons.add('existing-edit');
  active.indexes.legend = retainedLegendIndex;
  const replacement = { name: 'external.svg', content: '<svg data-owner="external"></svg>' };
  const outcome = fixture.runtime.commitDomEdit({
    reason: 'ownership-change-during-mutation',
    invalidateIndexes: ['legend'],
    mutate: ({ mutation }) => {
      mutation.setAttribute(fixture.oldFeature, 'fill', '#abcdef');
      fixture.state.results.value[0] = replacement;
      return true;
    }
  });
  assert.deepEqual(outcome, {
    changed: false,
    flushed: false,
    resultIndex: 0,
    reason: 'ownership-change-during-mutation'
  });
  assert.equal(fixture.state.results.value[0], replacement);
  assert.equal(fixture.oldFeature.getAttribute('fill'), '#111111');
  assert.equal(active.dirty, true);
  assert.deepEqual([...active.dirtyReasons], ['existing-edit']);
  assert.equal(active.indexes.legend, retainedLegendIndex);
}

{
  const earlier = { value: 'before-earlier' };
  let middleValue = 'before-middle';
  const middle = {};
  Object.defineProperty(middle, 'value', {
    configurable: true,
    enumerable: true,
    get: () => middleValue,
    set: (value) => {
      if (value === 'before-middle') throw new Error('middle undo failed');
      middleValue = value;
    }
  });
  const later = { value: 'before-later' };
  const mutation = createDomMutationJournal();
  mutation.setProperty(earlier, 'value', 'after-earlier');
  mutation.setProperty(middle, 'value', 'after-middle');
  mutation.setProperty(later, 'value', 'after-later');

  const rollbackErrors = mutation.rollback();
  assert.equal(earlier.value, 'before-earlier');
  assert.equal(middle.value, 'after-middle');
  assert.equal(later.value, 'before-later');
  assert.deepEqual(rollbackErrors.map((error) => error.message), ['middle undo failed']);
  assert.equal(mutation.rollback(), rollbackErrors, 'double rollback is idempotent');
  assert.throws(() => mutation.commit(), /rolled back/);
}

{
  const sourceParent = new FakeSvgElement('g', { id: 'source-parent' });
  const sourceBefore = sourceParent.appendChild(new FakeSvgElement('path', { id: 'source-before' }));
  const replacement = sourceParent.appendChild(new FakeSvgElement('path', { id: 'replacement' }));
  const sourceAfter = sourceParent.appendChild(new FakeSvgElement('path', { id: 'source-after' }));
  const targetParent = new FakeSvgElement('g', { id: 'target-parent' });
  const targetBefore = targetParent.appendChild(new FakeSvgElement('path', { id: 'target-before' }));
  const current = targetParent.appendChild(new FakeSvgElement('path', { id: 'current' }));
  const targetAfter = targetParent.appendChild(new FakeSvgElement('path', { id: 'target-after' }));
  const mutation = createDomMutationJournal();

  mutation.replaceChild(targetParent, replacement, current);
  assert.deepEqual(targetParent.children, [targetBefore, replacement, targetAfter]);
  assert.deepEqual(sourceParent.children, [sourceBefore, sourceAfter]);
  assert.deepEqual(mutation.rollback(), []);
  assert.deepEqual(targetParent.children, [targetBefore, current, targetAfter]);
  assert.deepEqual(sourceParent.children, [sourceBefore, replacement, sourceAfter]);
}

{
  const originalParent = new FakeSvgElement('g');
  const before = originalParent.appendChild(new FakeSvgElement('path', { id: 'before' }));
  const moved = originalParent.appendChild(new FakeSvgElement('path', { id: 'moved' }));
  const after = originalParent.appendChild(new FakeSvgElement('path', { id: 'after' }));
  const secondParent = new FakeSvgElement('g');
  const thirdParent = new FakeSvgElement('g');
  const thirdAnchor = thirdParent.appendChild(new FakeSvgElement('path', { id: 'third-anchor' }));
  const outer = createDomMutationJournal();
  const inner = createDomMutationJournal();

  outer.appendChild(secondParent, moved);
  inner.insertBefore(thirdParent, moved, thirdAnchor);
  assert.deepEqual(inner.rollback(), []);
  assert.deepEqual(secondParent.children, [moved]);
  assert.deepEqual(outer.rollback(), []);
  assert.deepEqual(originalParent.children, [before, moved, after]);
}

{
  const rollbackMetrics = [];
  const previousHooks = globalThis.__GBDRAW_TEST_HOOKS__;
  globalThis.__GBDRAW_TEST_HOOKS__ = {
    onStructuralMetric: (metric) => rollbackMetrics.push(metric)
  };
  const rollbackFeature = new FakeFeatureElement('rollback-feature', { fill: '#111111' });
  const rollbackSvg = makeSvg([rollbackFeature]);
  const originalResult = { name: 'rollback.svg', content: '<svg data-owner="original"></svg>' };
  const rollbackState = {
    results: ref([originalResult]),
    selectedResultIndex: ref(0),
    resultGenerationKey: ref('rollback-generation'),
    skipCaptureBaseConfig: ref(false),
    svgContainer: ref({ querySelector: (selector) => selector === 'svg' ? rollbackSvg : null })
  };
  const rollbackRuntime = createPreviewRuntime({
    state: rollbackState,
    serializeSvg: () => {
      throw new Error('serialization failed');
    }
  });
  rollbackRuntime.mountResultSvg(0, rollbackSvg);
  const retainedFeatureIndex = { retained: true };
  rollbackRuntime.getActiveRuntime().indexes.features = retainedFeatureIndex;
  let fragileValue = 'before';
  const fragile = {};
  Object.defineProperty(fragile, 'value', {
    configurable: true,
    enumerable: true,
    get: () => fragileValue,
    set: (value) => {
      if (value === 'before') throw new Error('rollback setter failed');
      fragileValue = value;
    }
  });

  let caught = null;
  try {
    rollbackRuntime.commitDomEdit({
      reason: 'rollback-error-restoration',
      invalidateIndexes: ['features'],
      mutate: ({ mutation }) => {
        mutation.setAttribute(rollbackFeature, 'fill', '#abcdef');
        mutation.setProperty(fragile, 'value', 'after');
        return true;
      }
    });
  } catch (error) {
    caught = error;
  }
  assert.equal(caught?.message, 'serialization failed');
  assert.deepEqual(caught?.rollbackErrors?.map((error) => error.message), ['rollback setter failed']);
  assert.equal(rollbackFeature.getAttribute('fill'), '#111111');
  assert.equal(fragile.value, 'after');
  assert.equal(rollbackState.results.value[0], originalResult);
  assert.equal(rollbackRuntime.getActiveRuntime().dirty, false);
  assert.deepEqual([...rollbackRuntime.getActiveRuntime().dirtyReasons], []);
  assert.equal(rollbackRuntime.getActiveRuntime().indexes.features, retainedFeatureIndex);
  assert.equal(rollbackState.skipCaptureBaseConfig.value, false);
  const rollbackMetricNames = new Set(rollbackMetrics.map(({ name }) => name));
  assert.equal(rollbackMetricNames.has('liveEditDomRollbackCount'), true);
  assert.equal(rollbackMetricNames.has('liveEditCanonicalRollbackCount'), true);
  assert.equal(rollbackMetricNames.has('liveEditCanonicalRollbackFailureCount'), true);
  assert.equal(rollbackMetricNames.has('liveEditResultRestoreCount'), true);
  if (previousHooks === undefined) delete globalThis.__GBDRAW_TEST_HOOKS__;
  else globalThis.__GBDRAW_TEST_HOOKS__ = previousHooks;
}

{
  const resultSetterFeature = new FakeFeatureElement('result-setter-feature', { fill: '#111111' });
  const resultSetterSvg = makeSvg([resultSetterFeature]);
  const originalResult = { name: 'setter.svg', content: '<svg data-owner="original"></svg>' };
  const resultItems = [originalResult];
  let replacementAttempts = 0;
  let restoreAttempts = 0;
  const proxiedResults = new Proxy(resultItems, {
    set(target, property, value) {
      if (String(property) !== '0') return Reflect.set(target, property, value);
      if (value !== originalResult) {
        replacementAttempts += 1;
        Reflect.set(target, property, value);
        throw new Error('primary Result setter failure');
      }
      restoreAttempts += 1;
      if (restoreAttempts === 1) throw new Error('Result rollback setter failure');
      return Reflect.set(target, property, value);
    }
  });
  const resultSetterState = {
    results: ref(proxiedResults),
    selectedResultIndex: ref(0),
    resultGenerationKey: ref('setter-generation'),
    skipCaptureBaseConfig: ref(false),
    svgContainer: ref({
      querySelector: (selector) => selector === 'svg' ? resultSetterSvg : null
    })
  };
  const resultSetterRuntime = createPreviewRuntime({
    state: resultSetterState,
    serializeSvg: () => '<svg data-owner="replacement"></svg>'
  });
  resultSetterRuntime.mountResultSvg(0, resultSetterSvg);
  const retainedFeatureIndex = { retained: true };
  resultSetterRuntime.getActiveRuntime().indexes.features = retainedFeatureIndex;

  let caught = null;
  try {
    resultSetterRuntime.commitDomEdit({
      reason: 'partial-result-setter-failure',
      invalidateIndexes: ['features'],
      mutate: ({ mutation }) => mutation.setAttribute(resultSetterFeature, 'fill', '#abcdef')
    });
  } catch (error) {
    caught = error;
  }

  assert.equal(caught?.message, 'primary Result setter failure');
  assert.deepEqual(
    caught?.rollbackErrors?.map((error) => error.message),
    ['Result rollback setter failure']
  );
  assert.equal(replacementAttempts, 1);
  assert.equal(restoreAttempts, 2, 'the outer transaction retries the failed Result restore');
  assert.equal(resultSetterState.results.value[0], originalResult);
  assert.equal(resultSetterFeature.getAttribute('fill'), '#111111');
  assert.equal(resultSetterRuntime.getActiveRuntime().dirty, false);
  assert.deepEqual([...resultSetterRuntime.getActiveRuntime().dirtyReasons], []);
  assert.equal(resultSetterRuntime.getActiveRuntime().indexes.features, retainedFeatureIndex);
  assert.equal(resultSetterState.skipCaptureBaseConfig.value, false);
}

{
  const vueContext = {};
  runInNewContext(
    await readFile(join(repoRoot, 'gbdraw', 'web', 'vendor', 'vue', 'vue.global.js'), 'utf8'),
    vueContext
  );
  const vueResults = vueContext.Vue.ref([
    { name: 'reactive.svg', content: '<svg data-owner="original"></svg>' }
  ]);
  const originalOwner = vueResults.value[0];
  const reactiveFeature = new FakeFeatureElement('reactive-result-feature', { fill: '#111111' });
  const reactiveSvg = makeSvg([reactiveFeature]);
  const reactiveState = {
    results: vueResults,
    selectedResultIndex: ref(0),
    resultGenerationKey: ref('reactive-generation'),
    skipCaptureBaseConfig: ref(false),
    svgContainer: ref({ querySelector: (selector) => selector === 'svg' ? reactiveSvg : null })
  };
  const previousHooks = globalThis.__GBDRAW_TEST_HOOKS__;
  globalThis.__GBDRAW_TEST_HOOKS__ = {
    onStructuralMetric: ({ name }) => {
      if (name === 'domEditResultReplacementCount') throw new Error('metric failed');
    }
  };
  const reactiveRuntime = createPreviewRuntime({
    state: reactiveState,
    serializeSvg: () => '<svg data-owner="replacement"></svg>'
  });
  reactiveRuntime.mountResultSvg(0, reactiveSvg);

  assert.throws(() => reactiveRuntime.commitDomEdit({
    reason: 'reactive-result-owner-failure',
    mutate: ({ mutation }) => mutation.setAttribute(reactiveFeature, 'fill', '#abcdef')
  }), /metric failed/);

  assert.equal(reactiveState.results.value[0], originalOwner);
  assert.equal(reactiveState.results.value[0].content, '<svg data-owner="original"></svg>');
  assert.equal(reactiveRuntime.getActiveRuntime().resultOwner, originalOwner);
  assert.equal(reactiveFeature.getAttribute('fill'), '#111111');
  assert.equal(reactiveRuntime.getActiveRuntime().dirty, false);
  if (previousHooks === undefined) delete globalThis.__GBDRAW_TEST_HOOKS__;
  else globalThis.__GBDRAW_TEST_HOOKS__ = previousHooks;
}

{
  const nested = { value: 'before' };
  const items = [nested, 'retained'];
  const mapped = new Map([['retained', nested]]);
  const selected = new Set(['retained']);
  const stateOwner = {
    scalar: 'before',
    removed: 'retained',
    nested,
    items,
    mapped,
    selected
  };
  const mutation = createDomMutationJournal();
  mutation.captureState(stateOwner);
  stateOwner.scalar = 'after';
  delete stateOwner.removed;
  stateOwner.added = true;
  stateOwner.nested.value = 'after';
  stateOwner.items = ['replacement'];
  mapped.delete('retained');
  mapped.set('added', 'after');
  selected.delete('retained');
  selected.add('added');

  assert.deepEqual(mutation.rollback(), []);
  assert.equal(stateOwner.scalar, 'before');
  assert.equal(stateOwner.removed, 'retained');
  assert.equal(Object.hasOwn(stateOwner, 'added'), false);
  assert.equal(stateOwner.nested, nested);
  assert.equal(nested.value, 'before');
  assert.equal(stateOwner.items, items);
  assert.deepEqual(items, [nested, 'retained']);
  assert.deepEqual([...mapped], [['retained', nested]]);
  assert.deepEqual([...selected], ['retained']);
}

{
  const originalParent = new FakeSvgElement('g');
  const before = originalParent.appendChild(new FakeSvgElement('path', { id: 'before' }));
  const moved = originalParent.appendChild(new FakeSvgElement('path', { id: 'moved' }));
  const after = originalParent.appendChild(new FakeSvgElement('path', { id: 'after' }));
  const middleParent = new FakeSvgElement('g');
  const finalParent = new FakeSvgElement('g');
  const mutation = createDomMutationJournal();
  mutation.appendChild(middleParent, moved);
  mutation.appendChild(finalParent, moved);
  mutation.setAttribute(moved, 'fill', '#111111');
  mutation.setAttribute(moved, 'fill', '#222222');

  assert.deepEqual(mutation.rollback(), []);
  assert.deepEqual(originalParent.children, [before, moved, after]);
  assert.deepEqual(middleParent.children, []);
  assert.deepEqual(finalParent.children, []);
  assert.equal(moved.getAttribute('fill'), null);
}

{
  const sourceParent = new FakeSvgElement('g');
  const sourceBefore = sourceParent.appendChild(new FakeSvgElement('path', { id: 'source-before' }));
  const incoming = sourceParent.appendChild(new FakeSvgElement('path', { id: 'incoming' }));
  const sourceAfter = sourceParent.appendChild(new FakeSvgElement('path', { id: 'source-after' }));
  const targetParent = new FakeSvgElement('g');
  const first = targetParent.appendChild(new FakeSvgElement('path', { id: 'first' }));
  const second = targetParent.appendChild(new FakeSvgElement('path', { id: 'second' }));
  const mutation = createDomMutationJournal();
  mutation.replaceChildren(targetParent, incoming);

  assert.deepEqual(mutation.rollback(), []);
  assert.deepEqual(targetParent.children, [first, second]);
  assert.deepEqual(sourceParent.children, [sourceBefore, incoming, sourceAfter]);
}

{
  const sourceParent = new FakeSvgElement('g');
  const firstSource = sourceParent.appendChild(new FakeSvgElement('path', { id: 'first-source' }));
  const secondSource = sourceParent.appendChild(new FakeSvgElement('path', { id: 'second-source' }));
  const thirdSource = sourceParent.appendChild(new FakeSvgElement('path', { id: 'third-source' }));
  const targetParent = new FakeSvgElement('g');
  const retained = targetParent.appendChild(new FakeSvgElement('path', { id: 'retained' }));
  const mutation = createDomMutationJournal();
  mutation.replaceChildren(targetParent, secondSource, firstSource);

  assert.deepEqual(mutation.rollback(), []);
  assert.deepEqual(targetParent.children, [retained]);
  assert.deepEqual(sourceParent.children, [firstSource, secondSource, thirdSource]);
}

console.log('preview runtime tests passed');
