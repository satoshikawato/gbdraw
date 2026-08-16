import assert from 'node:assert/strict';
import { mkdir, mkdtemp, readFile, writeFile } from 'node:fs/promises';
import { tmpdir } from 'node:os';
import { join } from 'node:path';
import { pathToFileURL } from 'node:url';

const repoRoot = process.cwd();
const sourcePath = join(repoRoot, 'gbdraw', 'web', 'js', 'app', 'preview-runtime.js');
const featureDomSourcePath = join(repoRoot, 'gbdraw', 'web', 'js', 'app', 'feature-dom.js');
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
  join(tempDir, 'app', 'dom-mutation-journal.js'),
  await readFile(mutationJournalSourcePath, 'utf8'),
  'utf8'
);
await Promise.all(serviceNames.map(async (name) => {
  const servicePath = join(repoRoot, 'gbdraw', 'web', 'js', 'services', name);
  await writeFile(join(tempDir, 'services', name), await readFile(servicePath, 'utf8'), 'utf8');
}));

const { createPreviewRuntime } = await import(pathToFileURL(join(tempDir, 'app', 'preview-runtime.js')));

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
const svg = makeSvg([featureA, featureBConnector, featureBBlock]);
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
  assert.equal(await settlement, false, raceCase.name);
  assert.equal(currentOwner.content, currentContent, raceCase.name);
  assert.equal(fixture.oldFeature.getAttribute('fill'), '#111111', raceCase.name);
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

console.log('preview runtime tests passed');
