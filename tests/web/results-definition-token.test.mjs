import assert from 'node:assert/strict';
import { mkdir, mkdtemp, readFile, writeFile } from 'node:fs/promises';
import { tmpdir } from 'node:os';
import { join } from 'node:path';
import { pathToFileURL } from 'node:url';
import test from 'node:test';

const repoRoot = process.cwd();
const tempDir = await mkdtemp(join(tmpdir(), 'gbdraw-results-definition-token-'));
await writeFile(join(tempDir, 'package.json'), '{"type":"module"}\n', 'utf8');
await mkdir(join(tempDir, 'app', 'legend-layout'), { recursive: true });
await mkdir(join(tempDir, 'app', 'legend'), { recursive: true });
await mkdir(join(tempDir, 'services'), { recursive: true });
await writeFile(
  join(tempDir, 'app', 'results.js'),
  await readFile(join(repoRoot, 'gbdraw', 'web', 'js', 'app', 'results.js'), 'utf8'),
  'utf8'
);
await writeFile(
  join(tempDir, 'app', 'legend-layout', 'transform-utils.js'),
  'export const parseTransform = () => ({ x: 0, y: 0 });\n',
  'utf8'
);
await writeFile(
  join(tempDir, 'app', 'legend-layout', 'composition-actions.js'),
  "export const COMPOSITION_ROLE_ATTRIBUTE = 'data-gbdraw-composition-role';\n",
  'utf8'
);
await writeFile(
  join(tempDir, 'app', 'legend', 'utils.js'),
  "export const COMPARISON_LEGEND_SELECTOR = '[data-gbdraw-role=\"comparison-legend\"]';\n",
  'utf8'
);
await writeFile(
  join(tempDir, 'app', 'record-groups.js'),
  'export const isMultiRecordCanvasSvg = () => false;\n',
  'utf8'
);
await writeFile(
  join(tempDir, 'services', 'diagram-generation.js'),
  [
    "export const DIAGRAM_HELPER_OPERATIONS = { REGENERATE_DEFINITION_SVGS: 'regenerateDefinitionSvgs' };",
    'export const runDiagramHelperOperation = (...args) => globalThis.__GBDRAW_RESULTS_HELPER__(...args);',
    ''
  ].join('\n'),
  'utf8'
);
await writeFile(
  join(tempDir, 'services', 'file-content-cache.js'),
  'export const cloneFileBytesForTransfer = (...args) => globalThis.__GBDRAW_CLONE_FILE_BYTES__(...args);\n',
  'utf8'
);

class FakeElement {
  constructor(tagName, attributes = {}) {
    this.tagName = tagName;
    this.attributes = new Map(
      Object.entries(attributes).map(([name, value]) => [name, String(value)])
    );
    this.children = [];
    this.parentNode = null;
    this.ownerDocument = null;
  }

  get id() {
    return this.getAttribute('id') || '';
  }

  getAttribute(name) {
    return this.attributes.has(name) ? this.attributes.get(name) : null;
  }

  setAttribute(name, value) {
    this.attributes.set(name, String(value));
  }

  appendChild(child) {
    child.parentNode = this;
    child.ownerDocument = this.ownerDocument;
    this.children.push(child);
    return child;
  }

  replaceChild(replacement, current) {
    const index = this.children.indexOf(current);
    if (index < 0) throw new Error('missing replacement target');
    current.parentNode = null;
    replacement.parentNode = this;
    replacement.ownerDocument = this.ownerDocument;
    this.children[index] = replacement;
    return current;
  }

  cloneNode(deep = false) {
    const clone = new FakeElement(this.tagName, Object.fromEntries(this.attributes));
    if (deep) this.children.forEach((child) => clone.appendChild(child.cloneNode(true)));
    return clone;
  }

  querySelectorAll(selector) {
    if (selector === 'g[data-gbdraw-role="record-definition"]') {
      return this.children.filter(
        (child) => child.getAttribute('data-gbdraw-role') === 'record-definition'
      );
    }
    return [];
  }

  getElementById(id) {
    return this.children.find((child) => child.id === id) || null;
  }
}

const createSvg = (documentName) => {
  const ownerDocument = {
    name: documentName,
    importNode: (node, deep) => {
      const clone = node.cloneNode(deep);
      clone.ownerDocument = ownerDocument;
      return clone;
    }
  };
  const svg = new FakeElement('svg');
  svg.ownerDocument = ownerDocument;
  const definition = new FakeElement('g', {
    id: 'definition-a',
    'data-gbdraw-role': 'record-definition',
    'data-gbdraw-record-index': '0',
    'data-version': 'old'
  });
  svg.appendChild(definition);
  return { definition, svg };
};

globalThis.DOMParser = class {
  parseFromString() {
    const group = new FakeElement('g', {
      id: 'definition-a',
      'data-gbdraw-role': 'record-definition',
      'data-gbdraw-record-index': '0',
      'data-version': 'settled'
    });
    return { querySelector: (selector) => selector === 'g' ? group : null };
  }
};

const { createResultsManager } = await import(
  pathToFileURL(join(tempDir, 'app', 'results.js'))
);
const { createPreviewRuntime } = await import(
  pathToFileURL(join(repoRoot, 'gbdraw', 'web', 'js', 'app', 'preview-runtime.js'))
);

const ref = (value) => ({ value });
const deferred = () => {
  let resolve;
  const promise = new Promise((resolvePromise) => {
    resolve = resolvePromise;
  });
  return { promise, resolve };
};

const createPreviewProtocol = ({ state, svg }) => {
  let active = {
    documentAdmissionIdentity: {},
    generationKey: String(state.resultGenerationKey.value),
    mountRevision: 1,
    result: state.results.value[0],
    resultIndex: 0,
    resultsOwner: state.results.value,
    svg
  };
  let commitCount = 0;

  const isDomEditTokenCurrent = (token) => Boolean(
    token
    && token.resultsOwner === active.resultsOwner
    && token.result === active.result
    && token.resultIndex === active.resultIndex
    && token.generationKey === active.generationKey
    && token.mountRevision === active.mountRevision
    && token.svg === active.svg
    && token.documentAdmissionIdentity === active.documentAdmissionIdentity
  );

  return {
    admitSameElementDocument() {
      active = { ...active, documentAdmissionIdentity: {} };
    },
    captureDomEditToken() {
      return { ...active };
    },
    commitDomEdit({ mutate, targetToken = null }) {
      if (targetToken && !isDomEditTokenCurrent(targetToken)) {
        return { changed: false, flushed: false, resultIndex: targetToken.resultIndex };
      }
      const mutation = {
        appendChild: (parent, child) => parent.appendChild(child),
        removeNode: () => false,
        replaceChild: (parent, replacement, current) => parent.replaceChild(replacement, current)
      };
      const changed = Boolean(mutate({ svg: active.svg, mutation }));
      if (changed) commitCount += 1;
      return { changed, flushed: changed, resultIndex: active.resultIndex };
    },
    get commitCount() {
      return commitCount;
    },
    isDomEditTokenCurrent,
    mount(nextSvg) {
      active = {
        documentAdmissionIdentity: {},
        generationKey: String(state.resultGenerationKey.value),
        mountRevision: active.mountRevision + 1,
        result: state.results.value[0],
        resultIndex: 0,
        resultsOwner: state.results.value,
        svg: nextSvg
      };
    }
  };
};

const createFixture = ({ activeSvg, containerSvg }) => {
  const state = {
    adv: {
      def_font_size: null,
      keep_full_definition_with_plot_title: false,
      plot_title_font_size: null,
      plot_title_position: 'none'
    },
    appliedPaletteColors: ref({}),
    appliedPaletteName: ref('default'),
    cInputType: ref('gb'),
    currentColors: ref({}),
    files: { c_gb: { name: 'source.gbk' } },
    form: { plot_title: '', species: '', strain: '' },
    linearSeqs: [],
    mode: ref('circular'),
    normalizePaletteColors: (colors) => colors,
    paletteDefinitions: ref({}),
    paletteInstantPreviewEnabled: ref(true),
    pendingPaletteColors: ref({}),
    pendingPaletteName: ref(''),
    resultGenerationKey: ref('generation-1'),
    results: ref([{ name: 'result.svg', content: '<svg></svg>' }]),
    selectedPalette: ref('default'),
    shouldDeferCircularPreviewUpdates: ref(false),
    svgContainer: ref({ querySelector: (selector) => selector === 'svg' ? containerSvg : null }),
    svgContent: ref('<svg></svg>')
  };
  const previewRuntime = createPreviewProtocol({ state, svg: activeSvg });
  const manager = createResultsManager({
    state,
    legendLayout: { refreshCompositionGeometry() {} },
    previewRuntime
  });
  let scheduled = null;
  const originalSetTimeout = globalThis.setTimeout;
  globalThis.setTimeout = (callback) => {
    scheduled = callback;
    return 1;
  };
  try {
    manager.scheduleDefinitionUpdate();
  } finally {
    globalThis.setTimeout = originalSetTimeout;
  }
  assert.equal(typeof scheduled, 'function');
  return { manager, previewRuntime, runScheduled: scheduled, state };
};

const helperResponse = {
  result: {
    definitions: [{
      definition_group_id: 'definition-a',
      record_index: 0,
      svg: '<g id="definition-a"></g>'
    }]
  }
};

test('definition update started before preview remount cannot settle onto the remounted Result', async () => {
  const oldDocument = createSvg('old-document');
  const nextDocument = createSvg('next-document');
  const cloneGate = deferred();
  let cloneCalls = 0;
  let helperCalls = 0;
  globalThis.__GBDRAW_CLONE_FILE_BYTES__ = () => {
    cloneCalls += 1;
    return cloneGate.promise;
  };
  globalThis.__GBDRAW_RESULTS_HELPER__ = async () => {
    helperCalls += 1;
    return helperResponse;
  };
  const fixture = createFixture({
    activeSvg: oldDocument.svg,
    containerSvg: nextDocument.svg
  });

  fixture.runScheduled();
  await Promise.resolve();
  fixture.previewRuntime.mount(nextDocument.svg);
  cloneGate.resolve(new Uint8Array([1, 2, 3]));
  await new Promise(setImmediate);

  assert.equal(cloneCalls, 0, 'a pre-remount target is rejected before the first await');
  assert.equal(helperCalls, 0);
  assert.equal(fixture.previewRuntime.commitCount, 0);
  assert.equal(nextDocument.svg.getElementById('definition-a').getAttribute('data-version'), 'old');
});

test('same SVG element with a new document admission rejects a suspended definition settlement', async () => {
  const document = createSvg('reused-element-document');
  const helperGate = deferred();
  let helperCalls = 0;
  globalThis.__GBDRAW_CLONE_FILE_BYTES__ = async () => new Uint8Array([1, 2, 3]);
  globalThis.__GBDRAW_RESULTS_HELPER__ = () => {
    helperCalls += 1;
    return helperGate.promise;
  };
  const fixture = createFixture({ activeSvg: document.svg, containerSvg: document.svg });

  fixture.runScheduled();
  await Promise.resolve();
  await Promise.resolve();
  assert.equal(helperCalls, 1);
  fixture.previewRuntime.admitSameElementDocument();
  helperGate.resolve(helperResponse);
  await new Promise(setImmediate);

  assert.equal(fixture.previewRuntime.commitCount, 0);
  assert.equal(document.svg.getElementById('definition-a').getAttribute('data-version'), 'old');
});

const createRealPreviewFixture = ({ document, secondDocument = null } = {}) => {
  const firstResult = { name: 'result.svg', content: '<svg data-owner="captured"></svg>' };
  const secondResult = { name: 'second.svg', content: '<svg data-owner="second"></svg>' };
  const state = {
    adv: {
      def_font_size: null,
      keep_full_definition_with_plot_title: false,
      plot_title_font_size: null,
      plot_title_position: 'none'
    },
    appliedPaletteColors: ref({}),
    appliedPaletteName: ref('default'),
    cInputType: ref('gb'),
    currentColors: ref({}),
    files: { c_gb: { name: 'source.gbk' } },
    form: { plot_title: '', species: '', strain: '' },
    linearSeqs: [],
    mode: ref('circular'),
    normalizePaletteColors: (colors) => colors,
    paletteDefinitions: ref({}),
    paletteInstantPreviewEnabled: ref(true),
    pendingPaletteColors: ref({}),
    pendingPaletteName: ref(''),
    resultGenerationKey: ref('generation-1'),
    results: ref([firstResult, secondResult]),
    selectedPalette: ref('default'),
    selectedResultIndex: ref(0),
    shouldDeferCircularPreviewUpdates: ref(false),
    skipCaptureBaseConfig: ref(false),
    svgContainer: ref({
      querySelector: (selector) => selector === 'svg' ? document.svg : null
    }),
    svgContent: ref('<svg></svg>')
  };
  const previewRuntime = createPreviewRuntime({
    state,
    serializeSvg: (svg) => {
      const version = svg.getElementById('definition-a')?.getAttribute('data-version') || '';
      return `<svg data-version="${version}"></svg>`;
    }
  });
  previewRuntime.mountResultSvg(0, document.svg);
  const manager = createResultsManager({
    state,
    legendLayout: { refreshCompositionGeometry() {} },
    previewRuntime
  });
  let scheduled = null;
  const originalSetTimeout = globalThis.setTimeout;
  globalThis.setTimeout = (callback) => {
    scheduled = callback;
    return 1;
  };
  try {
    manager.scheduleDefinitionUpdate();
  } finally {
    globalThis.setTimeout = originalSetTimeout;
  }
  return {
    firstResult,
    manager,
    previewRuntime,
    runScheduled: scheduled,
    secondDocument,
    secondResult,
    state
  };
};

test('same-object Result content replacement makes a suspended definition token stale', async () => {
  const document = createSvg('same-result-owner');
  const helperGate = deferred();
  let helperCalls = 0;
  globalThis.__GBDRAW_CLONE_FILE_BYTES__ = async () => new Uint8Array([1, 2, 3]);
  globalThis.__GBDRAW_RESULTS_HELPER__ = () => {
    helperCalls += 1;
    return helperGate.promise;
  };
  const fixture = createRealPreviewFixture({ document });
  const canonicalState = { untouched: true };
  const beforeDomVersion = document.definition.getAttribute('data-version');

  fixture.runScheduled();
  await Promise.resolve();
  await Promise.resolve();
  assert.equal(helperCalls, 1);
  fixture.firstResult.content = '<svg data-owner="newer-same-object"></svg>';
  helperGate.resolve(helperResponse);
  await new Promise(setImmediate);
  await new Promise(setImmediate);

  assert.equal(fixture.state.results.value[0], fixture.firstResult);
  assert.equal(
    fixture.firstResult.content,
    '<svg data-owner="newer-same-object"></svg>'
  );
  assert.equal(document.definition.getAttribute('data-version'), beforeDomVersion);
  assert.deepEqual(canonicalState, { untouched: true });
});

test('Result content token controls distinguish value equality and owner replacement', () => {
  const firstDocument = createSvg('token-controls-first');
  const secondDocument = createSvg('token-controls-second');
  const fixture = createRealPreviewFixture({
    document: firstDocument,
    secondDocument
  });

  const sameContentToken = fixture.previewRuntime.captureDomEditToken();
  fixture.firstResult.content = String(fixture.firstResult.content);
  assert.equal(fixture.previewRuntime.isDomEditTokenCurrent(sameContentToken), true);

  const replacedObjectToken = fixture.previewRuntime.captureDomEditToken();
  fixture.state.results.value[0] = { ...fixture.firstResult };
  assert.equal(fixture.previewRuntime.isDomEditTokenCurrent(replacedObjectToken), false);
  fixture.state.results.value[0] = fixture.firstResult;

  const replacedArrayToken = fixture.previewRuntime.captureDomEditToken();
  fixture.state.results.value = [...fixture.state.results.value];
  assert.equal(fixture.previewRuntime.isDomEditTokenCurrent(replacedArrayToken), false);

  fixture.state.results.value = [fixture.firstResult, fixture.secondResult];
  fixture.previewRuntime.mountResultSvg(0, firstDocument.svg);
  const awayAndBackToken = fixture.previewRuntime.captureDomEditToken();
  fixture.state.svgContainer.value = {
    querySelector: (selector) => selector === 'svg' ? secondDocument.svg : null
  };
  fixture.previewRuntime.selectResult(1);
  fixture.previewRuntime.mountResultSvg(1, secondDocument.svg);
  fixture.previewRuntime.selectResult(0);
  fixture.state.svgContainer.value = {
    querySelector: (selector) => selector === 'svg' ? firstDocument.svg : null
  };
  fixture.previewRuntime.mountResultSvg(0, firstDocument.svg);
  assert.equal(fixture.previewRuntime.isDomEditTokenCurrent(awayAndBackToken), false);
});

test('a Result replaced before remount cannot authorize edits from the old mounted SVG', () => {
  const document = createSvg('unmounted-replacement');
  const fixture = createRealPreviewFixture({ document });
  const replacement = {
    ...fixture.firstResult,
    content: '<svg data-owner="replacement-before-remount"></svg>'
  };

  fixture.state.results.value[0] = replacement;

  assert.equal(fixture.previewRuntime.captureDomEditToken(), null);
  assert.deepEqual(fixture.previewRuntime.commitDomEdit({
    reason: 'unmounted-replacement-probe',
    mutate: () => {
      document.definition.setAttribute('data-version', 'must-not-run');
      return true;
    }
  }), {
    changed: false,
    flushed: false,
    resultIndex: null,
    reason: 'unmounted-replacement-probe'
  });
  assert.equal(
    replacement.content,
    '<svg data-owner="replacement-before-remount"></svg>'
  );
  assert.notEqual(document.definition.getAttribute('data-version'), 'must-not-run');
});

test('same-object content replaced before capture cannot authorize the old mounted SVG', async () => {
  const document = createSvg('same-owner-before-capture');
  const fixture = createRealPreviewFixture({ document });

  fixture.firstResult.content = '<svg data-owner="same-owner-new-content"></svg>';

  assert.equal(fixture.previewRuntime.captureDomEditToken(), null);
  assert.equal(fixture.previewRuntime.commitDomEdit({
    reason: 'same-owner-before-capture',
    mutate: () => {
      document.definition.setAttribute('data-version', 'must-not-run');
      return true;
    }
  }).changed, false);
  assert.equal(
    fixture.firstResult.content,
    '<svg data-owner="same-owner-new-content"></svg>'
  );
  assert.notEqual(document.definition.getAttribute('data-version'), 'must-not-run');
  const canonical = { mode: 'default' };
  assert.throws(() => fixture.previewRuntime.runDomEditSync({
    reason: 'same-owner-canonical-sync',
    canonicalState: [{ target: canonical }],
    action: () => {
      canonical.mode = 'off';
    }
  }), (error) => error?.name === 'StaleDomEditError');
  await assert.rejects(fixture.previewRuntime.runDomEditTransaction({
    reason: 'same-owner-canonical-async',
    canonicalState: [{ target: canonical }],
    action: async () => {
      canonical.mode = 'on';
    }
  }), (error) => error?.name === 'StaleDomEditError');
  assert.deepEqual(canonical, { mode: 'default' });
});
