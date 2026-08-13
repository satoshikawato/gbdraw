import assert from 'node:assert/strict';
import { readFile } from 'node:fs/promises';
import test from 'node:test';

import {
  deriveUsedFeatureFillGroupsByResult,
  measureLegendTextInBrowser,
  prepareFeatureFillLegendProjection,
  prepareFeatureFillLegendProjections,
  preserveResultLocalNonFeatureLegendEntries
} from '../../gbdraw/web/js/app/legend/feature-fill-projection.js';
import {
  extractResultLocalLegendEntries
} from '../../gbdraw/web/js/app/feature-editor/fill-actions.js';

class Element {
  constructor(localName, attributes = {}, textContent = '') {
    this.localName = localName;
    this.tagName = localName;
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
    if (child.parentElement) child.parentElement.children = child.parentElement.children.filter((entry) => entry !== child);
    child.parentElement = this;
    this.children.push(child);
    return child;
  }
  remove() {
    if (this.parentElement) this.parentElement.children = this.parentElement.children.filter((entry) => entry !== this);
    this.parentElement = null;
  }
  cloneNode(deep = false) {
    const next = new Element(this.localName, Object.fromEntries(this.attributes), this.textContent);
    if (deep) this.children.forEach((child) => next.appendChild(child.cloneNode(true)));
    return next;
  }
  querySelectorAll(selector) {
    const matches = [];
    const match = (node) => {
      if (selector === 'path') return node.localName === 'path';
      if (selector === 'text') return node.localName === 'text';
      if (selector === '[transform]') return node.hasAttribute('transform');
      if (selector === 'g[data-legend-key]') return node.localName === 'g' && node.hasAttribute('data-legend-key');
      return false;
    };
    const visit = (node) => node.children.forEach((child) => {
      if (match(child)) matches.push(child);
      visit(child);
    });
    visit(this);
    return matches;
  }
  querySelector(selector) { return this.querySelectorAll(selector)[0] || null; }
}

const entry = (caption, color, y, { stroke = null, strokeWidth = null, owner = null } = {}) => {
  const group = new Element('g', { 'data-legend-key': caption });
  if (owner) group.setAttribute('data-legend-owner', owner);
  group.appendChild(new Element('path', {
    fill: color,
    transform: `translate(0, ${y})`,
    ...(stroke ? { stroke } : {}),
    ...(strokeWidth !== null ? { 'stroke-width': strokeWidth } : {})
  }));
  group.appendChild(new Element('text', { transform: `translate(20, ${y})` }, caption));
  return group;
};

const sourceLegend = (...entries) => {
  const svg = new Element('svg');
  const target = new Element('g', { id: 'feature_legend' });
  entries.forEach((item) => target.appendChild(item));
  svg.appendChild(target);
  return { svg, target };
};

const metadata = () => ({ legendReflow: { colorRectSize: 12, lineHeight: 16, textXOffset: 18 } });
const targets = (svg) => [svg.children[0]];
const swatchColor = (svg, caption) => (
  svg.querySelectorAll('g[data-legend-key]')
    .find((group) => group.getAttribute('data-legend-key') === caption)
    ?.querySelector('path')?.getAttribute('fill')
);
const swatchStroke = (svg, caption) => {
  const swatch = svg.querySelectorAll('g[data-legend-key]')
    .find((group) => group.getAttribute('data-legend-key') === caption)
    ?.querySelector('path');
  return {
    color: swatch?.getAttribute('stroke'),
    width: swatch?.getAttribute('stroke-width')
  };
};

const feature = ({ recordKey, id, type = 'CDS', product = 'Kinase', instanceHash, svgId }) => ({
  biological: {
    recordKey,
    biologicalFeatureId: id,
    instanceHash,
    type,
    qualifiers: { product: [product] }
  },
  rendered: { recordKey, biologicalFeatureId: id, svgId }
});

test('used groups are Result-local, omit no-fill/unused rules, and merge global manual-only intent', () => {
  const a = feature({ recordKey: 'a', id: 'fa', instanceHash: 'fi1_aaaaaaaaaaaaaaaaaaaaaaaaaa', svgId: 'svg-a' });
  const hidden = feature({ recordKey: 'a', id: 'fh', type: 'tRNA', product: 'Hidden', instanceHash: 'fi1_hhhhhhhhhhhhhhhhhhhhhhhhhh', svgId: 'svg-h' });
  const b = feature({ recordKey: 'b', id: 'fb', instanceHash: 'fi1_bbbbbbbbbbbbbbbbbbbbbbbbbb', svgId: 'svg-b' });
  const catalog = { items: [{
    resultKey: 'result-a', resultIndex: 0, resultName: 'a.svg',
    biologicalFeatures: [a.biological, hidden.biological], features: [a.rendered, hidden.rendered]
  }, {
    resultKey: 'result-b', resultIndex: 1, resultName: 'b.svg',
    biologicalFeatures: [b.biological], features: [b.rendered]
  }] };
  const projection = deriveUsedFeatureFillGroupsByResult({
    catalog,
    rules: [{ feat: 'CDS', qual: 'product', val: '^Kinase$', color: '#123456', cap: 'Core' }, {
      feat: 'tRNA', qual: 'product', val: '^Hidden$', color: 'none', cap: 'Hidden'
    }, {
      feat: 'gene', qual: 'product', val: '^Unused$', color: '#999999', cap: 'Unused'
    }],
    paletteColors: { CDS: '#cccccc', tRNA: '#dddddd' },
    manualLegendEntries: [{ caption: 'Notes', color: '#fedcba', showStroke: true }]
  });

  assert.deepEqual(projection.projections.map((item) => item.entries.map((legend) => ({
    caption: legend.caption, color: legend.color, featureIds: legend.featureIds, owner: legend.owner
  }))), [[{
    caption: 'Core', color: '#123456', featureIds: ['svg-a'], owner: 'feature'
  }, {
    caption: 'Notes', color: '#fedcba', featureIds: [], owner: 'manual'
  }], [{
    caption: 'Core', color: '#123456', featureIds: ['svg-b'], owner: 'feature'
  }, {
    caption: 'Notes', color: '#fedcba', featureIds: [], owner: 'manual'
  }]]);
  assert.equal(projection.counters.renderedFeatures, 3);
  assert.equal(projection.counters.ruleResolutionUpperBound, 9);
  assert.equal(Object.isFrozen(projection.projections[0].entries), true);
});

test('used groups discard stale effective captions when current rules fall back to the palette', () => {
  const target = feature({
    recordKey: 'a', id: 'fa', instanceHash: 'fi1_aaaaaaaaaaaaaaaaaaaaaaaaaa', svgId: 'svg-a'
  });
  target.biological.effectiveLegendCaption = 'Retired exact caption';
  const projection = deriveUsedFeatureFillGroupsByResult({
    catalog: { items: [{
      resultKey: 'result-a',
      biologicalFeatures: [target.biological],
      features: [target.rendered]
    }] },
    rules: [],
    paletteColors: { CDS: '#cccccc' }
  });
  assert.equal(projection.projections[0].entries[0].caption, 'CDS');
});

test('Feature legend replacement preserves Result-local built-ins but retires old Feature rows', () => {
  const prepared = preserveResultLocalNonFeatureLegendEntries({
    beforeProjection: {
      resultKey: 'result-a',
      entries: [{ caption: 'Core', color: '#111111', featureIds: ['svg-a'] }]
    },
    afterProjection: {
      resultKey: 'result-a',
      entries: [{ caption: 'Renamed core', color: '#abcdef', featureIds: ['svg-a'] }]
    },
    existingEntries: [{ caption: 'Core', color: '#111111', featureIds: ['svg-a'] }, {
      caption: 'GC content', color: '#777777', owner: 'built-in-palette', featureIds: []
    }]
  });
  assert.deepEqual(prepared.entries.map((entry) => entry.caption), [
    'GC content', 'Renamed core'
  ]);
  assert.equal(prepared.entries[0].owner, 'built-in-palette');
  assert.equal(Object.isFrozen(prepared.entries), true);
});

test('Result-local presentation bases preserve independent order and stroke through immediate batch staging', async () => {
  const firstSource = sourceLegend(
    entry('tRNA', '#222222', 0, { stroke: '#aa0000', strokeWidth: 3 }),
    entry('Notes', '#333333', 16, {
      stroke: '#cc0000', strokeWidth: 4, owner: 'manual'
    }),
    entry('CDS', '#111111', 32, { stroke: '#bb0000', strokeWidth: 2 })
  );
  const secondSource = sourceLegend(
    entry('CDS', '#111111', 0, { stroke: '#0000bb', strokeWidth: 1.5 }),
    entry('tRNA', '#222222', 16, { stroke: '#0000aa', strokeWidth: 1 }),
    entry('Notes', '#333333', 32, {
      stroke: '#0000cc', strokeWidth: 2.5, owner: 'manual'
    })
  );
  const firstExisting = extractResultLocalLegendEntries(firstSource.svg);
  const secondExisting = extractResultLocalLegendEntries(secondSource.svg);
  assert.deepEqual(firstExisting.map(({ caption }) => caption), ['tRNA', 'Notes', 'CDS']);
  assert.deepEqual(secondExisting.map(({ caption }) => caption), ['CDS', 'tRNA', 'Notes']);
  assert.deepEqual(
    firstExisting.map(({ strokeColor, strokeWidth }) => [strokeColor, strokeWidth]),
    [['#aa0000', 3], ['#cc0000', 4], ['#bb0000', 2]]
  );
  assert.deepEqual(
    secondExisting.map(({ strokeColor, strokeWidth }) => [strokeColor, strokeWidth]),
    [['#0000bb', 1.5], ['#0000aa', 1], ['#0000cc', 2.5]]
  );

  const firstCds = feature({
    recordKey: 'a', id: 'a-cds', type: 'CDS', instanceHash: 'fi1_aaaaaaaaaaaaaaaaaaaaaaaaaa', svgId: 'a-cds'
  });
  const firstTrna = feature({
    recordKey: 'a', id: 'a-trna', type: 'tRNA', instanceHash: 'fi1_bbbbbbbbbbbbbbbbbbbbbbbbbb', svgId: 'a-trna'
  });
  const secondCds = feature({
    recordKey: 'b', id: 'b-cds', type: 'CDS', instanceHash: 'fi1_cccccccccccccccccccccccccc', svgId: 'b-cds'
  });
  const secondTrna = feature({
    recordKey: 'b', id: 'b-trna', type: 'tRNA', instanceHash: 'fi1_dddddddddddddddddddddddddd', svgId: 'b-trna'
  });
  const catalog = { items: [{
    resultKey: 'result-a', resultIndex: 0,
    biologicalFeatures: [firstCds.biological, firstTrna.biological],
    features: [firstCds.rendered, firstTrna.rendered]
  }, {
    resultKey: 'result-b', resultIndex: 1,
    biologicalFeatures: [secondCds.biological, secondTrna.biological],
    features: [secondCds.rendered, secondTrna.rendered]
  }] };
  const derived = deriveUsedFeatureFillGroupsByResult({
    catalog,
    paletteColors: { CDS: '#abcdef', tRNA: '#fedcba' },
    manualLegendEntries: [{
      caption: 'Notes', color: '#333333', strokeColor: '#00aa00', strokeWidth: 5
    }],
    existingEntriesByResult: {
      'result-a': firstExisting,
      'result-b': secondExisting
    }
  });
  assert.deepEqual(
    derived.projections.map(({ entries }) => entries.map(({ caption }) => caption)),
    [['tRNA', 'Notes', 'CDS'], ['CDS', 'tRNA', 'Notes']]
  );
  assert.deepEqual(
    derived.projections.map(({ entries }) => entries.map((item) => [
      item.caption, item.strokeColor, item.strokeWidth
    ])),
    [[['tRNA', '#aa0000', 3], ['Notes', '#00aa00', 5], ['CDS', '#bb0000', 2]],
      [['CDS', '#0000bb', 1.5], ['tRNA', '#0000aa', 1], ['Notes', '#00aa00', 5]]]
  );

  const prepared = await prepareFeatureFillLegendProjections({
    sourcesByResult: new Map([
      ['result-a', firstSource.svg],
      ['result-b', secondSource.svg]
    ]),
    projections: derived.projections,
    parseMetadata: metadata,
    getTargetGroups: targets,
    measureText: () => { throw new Error('existing entries do not require measurement'); }
  });
  assert.deepEqual(
    prepared.staged.get('result-a').svg.querySelectorAll('g[data-legend-key]')
      .map((group) => group.getAttribute('data-legend-key')),
    ['tRNA', 'Notes', 'CDS']
  );
  assert.deepEqual(
    prepared.staged.get('result-b').svg.querySelectorAll('g[data-legend-key]')
      .map((group) => group.getAttribute('data-legend-key')),
    ['CDS', 'tRNA', 'Notes']
  );
  assert.deepEqual(swatchStroke(prepared.staged.get('result-a').svg, 'CDS'), {
    color: '#bb0000', width: '2'
  });
  assert.deepEqual(swatchStroke(prepared.staged.get('result-b').svg, 'CDS'), {
    color: '#0000bb', width: '1.5'
  });
  assert.deepEqual(swatchStroke(prepared.staged.get('result-a').svg, 'Notes'), {
    color: '#00aa00', width: '5'
  });
  assert.deepEqual(swatchStroke(prepared.staged.get('result-b').svg, 'Notes'), {
    color: '#00aa00', width: '5'
  });
  assert.deepEqual(swatchStroke(firstSource.svg, 'CDS'), { color: '#bb0000', width: 2 });
  assert.deepEqual(swatchStroke(secondSource.svg, 'CDS'), { color: '#0000bb', width: 1.5 });
});

test('existing swatch updates without text measurement and leaves source untouched', async () => {
  const { svg } = sourceLegend(entry('Core', '#111111', 0));
  const prepared = await prepareFeatureFillLegendProjection({
    sourceSvg: svg,
    entries: [{ caption: 'Core', color: '#abcdef', featureIds: ['svg-a'] }],
    parseMetadata: metadata,
    getTargetGroups: targets,
    measureText: () => { throw new Error('measurement should not run'); }
  });
  assert.equal(swatchColor(svg, 'Core'), '#111111');
  assert.equal(swatchColor(prepared.svg, 'Core'), '#abcdef');
  assert.equal(prepared.counters.measurements, 0);
  assert.equal(prepared.counters.colorUpdates, 1);
});

test('fill projection permits a diagram with the Feature legend disabled', async () => {
  const source = new Element('svg');
  source.appendChild(new Element('g', { id: 'features' }));
  const prepared = await prepareFeatureFillLegendProjection({
    sourceSvg: source,
    entries: [{ caption: 'CDS', color: '#123456', featureIds: ['svg-a'] }],
    allowAbsentLegend: true,
    parseMetadata: () => {
      throw new Error('absent legends do not require reflow metadata');
    },
    getTargetGroups: () => []
  });

  assert.notEqual(prepared.svg, source);
  assert.deepEqual(prepared.entries, []);
  assert.equal(prepared.counters.targetGroups, 0);
  assert.equal(prepared.svg.querySelectorAll('g[data-legend-key]').length, 0);
});

test('manual-style projection still rejects a missing Feature legend by default', async () => {
  const source = new Element('svg');
  await assert.rejects(
    prepareFeatureFillLegendProjection({
      sourceSvg: source,
      entries: [{ caption: 'Notes', color: '#123456', featureIds: [] }],
      getTargetGroups: () => []
    }),
    /no Feature legend target/
  );
});

test('browser text measurement waits for fonts, connects its host, and removes it', async () => {
  let fontsReady = false;
  const body = new Element('body');
  const document = {
    body,
    documentElement: body,
    fonts: {
      ready: Promise.resolve().then(() => { fontsReady = true; })
    },
    defaultView: {
      getComputedStyle: () => ({ getPropertyValue: () => '' })
    },
    createElementNS: (_namespace, localName) => {
      const node = new Element(localName);
      node.ownerDocument = document;
      if (localName === 'text') {
        node.getComputedTextLength = () => {
          assert.equal(fontsReady, true);
          assert.equal(node.parentElement?.parentElement, body);
          return 37;
        };
        node.getBBox = () => ({ width: 37, height: 11 });
      }
      return node;
    }
  };
  const svg = new Element('svg');
  svg.ownerDocument = document;
  const sourceText = new Element('text', { 'font-size': '12' }, 'Core');
  sourceText.ownerDocument = document;
  const measured = await measureLegendTextInBrowser({
    caption: 'Core', textElement: sourceText, svg
  });
  assert.deepEqual(measured, { width: 37, height: 11 });
  assert.equal(body.children.length, 0);
});

test('single-feature split creates a measured browser-owned entry from a template', async () => {
  const { svg } = sourceLegend(entry('Core', '#111111', 0));
  let measurements = 0;
  const prepared = await prepareFeatureFillLegendProjection({
    sourceSvg: svg,
    entries: [{ caption: 'Core', color: '#111111', featureIds: ['svg-b'] }, {
      caption: 'Core (2)', color: '#abcdef', featureIds: ['svg-a']
    }],
    parseMetadata: metadata,
    getTargetGroups: targets,
    measureText: async ({ caption }) => {
      measurements += 1;
      assert.equal(caption, 'Core (2)');
      return { width: 44, height: 12 };
    }
  });
  assert.equal(measurements, 1);
  assert.equal(prepared.counters.additions, 1);
  assert.deepEqual(
    prepared.svg.querySelectorAll('g[data-legend-key]').map((group) => group.getAttribute('data-legend-key')),
    ['Core', 'Core (2)']
  );
  assert.equal(swatchColor(prepared.svg, 'Core (2)'), '#abcdef');
});

test('horizontal additions measure the existing boundary caption and do not overlap it', async () => {
  const { svg } = sourceLegend(entry('Core', '#111111', 0));
  svg.children[0].setAttribute('id', 'feature_legend_h');
  const measured = [];
  const prepared = await prepareFeatureFillLegendProjection({
    sourceSvg: svg,
    entries: [{ caption: 'Core', color: '#111111' }, { caption: 'Core (2)', color: '#abcdef' }],
    parseMetadata: metadata,
    getTargetGroups: targets,
    measureText: async ({ caption }) => {
      measured.push(caption);
      return { width: caption === 'Core' ? 30 : 44, height: 12 };
    }
  });
  assert.deepEqual(measured, ['Core (2)', 'Core']);
  const added = prepared.svg.querySelectorAll('g[data-legend-key]')
    .find((group) => group.getAttribute('data-legend-key') === 'Core (2)');
  assert.equal(added.querySelector('text').getAttribute('transform'), 'translate(86, 0)');
});

test('feature rejoin removes the old group and preserves the existing destination group', async () => {
  const { svg } = sourceLegend(entry('Core', '#111111', 0), entry('Split', '#abcdef', 16));
  const prepared = await prepareFeatureFillLegendProjection({
    sourceSvg: svg,
    entries: [{ caption: 'Core', color: '#111111', featureIds: ['svg-a', 'svg-b'] }],
    parseMetadata: metadata,
    getTargetGroups: targets,
    measureText: () => { throw new Error('measurement should not run'); }
  });
  assert.deepEqual(
    prepared.svg.querySelectorAll('g[data-legend-key]').map((group) => group.getAttribute('data-legend-key')),
    ['Core']
  );
  assert.equal(prepared.counters.removals, 1);
});

test('batch legend preparation fails closed before publishing any source mutation', async () => {
  const first = sourceLegend(entry('Core', '#111111', 0));
  const second = sourceLegend(entry('Core', '#111111', 0));
  const beforeFirst = swatchColor(first.svg, 'Core');
  await assert.rejects(() => prepareFeatureFillLegendProjections({
    sourcesByResult: new Map([['result-a', first.svg], ['result-b', second.svg]]),
    projections: [{ resultKey: 'result-a', entries: [{ caption: 'Core', color: '#abcdef' }] }, {
      resultKey: 'result-b', entries: [{ caption: 'New caption', color: '#abcdef' }]
    }],
    parseMetadata: metadata,
    getTargetGroups: targets,
    measureText: async () => { throw new Error('font measurement failed'); }
  }), /font measurement failed/);
  assert.equal(swatchColor(first.svg, 'Core'), beforeFirst);
  assert.equal(swatchColor(second.svg, 'Core'), '#111111');
});

test('caption/color conflicts reject and atomic legend source contains zero Pyodide calls', async () => {
  const duplicate = feature({ recordKey: 'a', id: 'fb', type: 'gene', product: 'Other', instanceHash: 'fi1_cccccccccccccccccccccccccc', svgId: 'svg-b' });
  assert.throws(() => deriveUsedFeatureFillGroupsByResult({
    catalog: { items: [{ resultKey: 'result-a', biologicalFeatures: [
      { ...duplicate.biological, type: 'CDS', qualifiers: { product: ['Kinase'] } },
      duplicate.biological
    ], features: [
      { ...duplicate.rendered, biologicalFeatureId: 'fb', svgId: 'svg-a' },
      duplicate.rendered
    ] }] },
    rules: [{ feat: 'CDS', qual: 'product', val: 'Kinase', color: '#111111', cap: 'Same' }, {
      feat: 'gene', qual: 'product', val: 'Other', color: '#222222', cap: 'Same'
    }]
  }), /more than one color/);

  const source = await readFile(
    new URL('../../gbdraw/web/js/app/legend/feature-fill-projection.js', import.meta.url),
    'utf8'
  );
  assert.doesNotMatch(source, /Pyodide|pyodide|ensurePyodide|getPyodide/);
});

test('fill and bulk adapters derive every Result presentation base from its own SVG', async () => {
  const [fillSource, bulkSource] = await Promise.all([
    readFile(
      new URL('../../gbdraw/web/js/app/feature-editor/fill-actions.js', import.meta.url),
      'utf8'
    ),
    readFile(
      new URL('../../gbdraw/web/js/app/feature-editor/bulk-style-actions.js', import.meta.url),
      'utf8'
    )
  ]);
  assert.match(
    fillSource,
    /existingEntriesByResult\[resultKey\] = extractResultLocalLegendEntries\(source\)/
  );
  assert.doesNotMatch(
    fillSource,
    /existingEntriesByResult\[resultKey\]\s*=\s*cloneJson\(state\.legendEntries/
  );
  assert.match(bulkSource, /entriesByResult\[resultKey\] = extractResultLocalLegendEntries\(source\)/);
  assert.match(bulkSource, /existingEntriesByResult,/);
  assert.doesNotMatch(bulkSource, /entriesByResult\[resultKey\][^\n]+state\.legendEntries/);
});
