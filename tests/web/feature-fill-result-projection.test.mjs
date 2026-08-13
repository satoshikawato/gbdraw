import assert from 'node:assert/strict';
import test from 'node:test';

import { prepareFeatureFillResultProjection } from '../../gbdraw/web/js/app/feature-editor/fill-result-projection.js';
import { ingestSvgResults } from '../../gbdraw/web/js/services/svg-result-ingestion.js';

const featureTarget = (fill = '#111111') => ({
  fill,
  getAttribute(name) { return name === 'fill' ? this.fill : null; },
  setAttribute(name, value) { if (name === 'fill') this.fill = String(value); }
});

const root = ({ fill = '#111111', hasTarget = true } = {}) => {
  const target = featureTarget(fill);
  const svg = {
    localName: 'svg',
    target,
    style: { removeProperty() {} },
    cloneNode() { return root({ fill: this.target.fill, hasTarget }); },
    getAttribute(name) { return name === 'xmlns' ? 'http://www.w3.org/2000/svg' : null; },
    setAttribute() {},
    removeAttribute() {},
    querySelectorAll(selector) {
      if (selector.includes('path') || selector.includes('polygon') || selector.includes('rect')) {
        if (!hasTarget) return [];
        return [{
          ...target,
          getAttribute(name) {
            if (name === 'data-gbdraw-feature-id' || name === 'id') return 'svg-a';
            if (name === 'data-gbdraw-feature-part') return 'block';
            if (name === 'fill') return target.fill;
            return null;
          },
          setAttribute(name, value) { if (name === 'fill') target.fill = String(value); }
        }];
      }
      return [];
    }
  };
  return svg;
};

const catalog = { items: [{
  resultKey: 'result-a', resultIndex: 0, resultName: 'a.svg',
  biologicalFeatures: [{ recordKey: 'record-a', biologicalFeatureId: 'feature-a' }],
  features: [{ recordKey: 'record-a', biologicalFeatureId: 'feature-a', svgId: 'svg-a' }]
}] };
const targetKey = JSON.stringify(['result-a', 'record-a', 'feature-a']);

const installSerialization = () => {
  globalThis.XMLSerializer = class {
    serializeToString(svg) { return `<svg><path id="svg-a" fill="${svg.target.fill}"/></svg>`; }
  };
};

test('Result projection composes a legend-staged root with object fill values once', () => {
  installSerialization();
  const admitted = ingestSvgResults([{ name: 'a.svg', content: '<svg/>' }], {
    sanitizer: { sanitize: () => '<svg/>' },
    parser: class {
      parseFromString() { return { documentElement: root(), querySelector: () => null }; }
    }
  });
  const staged = root();
  const projection = prepareFeatureFillResultProjection({
    results: admitted,
    catalog,
    fillsByTargetKey: { [targetKey]: { color: '#abcdef', caption: 'Core' } },
    affectedResultKeys: ['result-a'],
    mounted: { resultIndex: 0, resultKey: 'result-a', svg: root() },
    preparedSvgByResultKey: new Map([['result-a', staged]]),
    targetFeatureKeysByResult: new Map([['result-a', [targetKey]]])
  });
  assert.equal(projection.nextResults[0].content.includes('#abcdef'), true);
  assert.equal(Object.isFrozen(projection.nextResults), false, 'published Results remain state-mutable');
  assert.equal(projection.preparedMountedSvg, staged);
  assert.ok(projection.admissionMetadataByResultKey['result-a'].before);
  assert.ok(projection.admissionMetadataByResultKey['result-a'].after);
  assert.deepEqual(projection.counters, {
    affectedResults: 1, mountedSerializations: 1, detachedPasses: 0, changedResults: 1
  });
});

test('Result projection fails closed when any planned rendered node is missing', () => {
  installSerialization();
  const admitted = ingestSvgResults([{ name: 'a.svg', content: '<svg/>' }], {
    sanitizer: { sanitize: () => '<svg/>' },
    parser: class {
      parseFromString() { return { documentElement: root(), querySelector: () => null }; }
    }
  });
  assert.throws(() => prepareFeatureFillResultProjection({
    results: admitted,
    catalog,
    fillsByTargetKey: { [targetKey]: { color: '#abcdef', caption: 'Core' } },
    affectedResultKeys: ['result-a'],
    preparedSvgByResultKey: new Map([['result-a', root({ hasTarget: false })]]),
    targetFeatureKeysByResult: new Map([['result-a', [targetKey]]])
  }), /SVG target is unavailable/);
  assert.equal(admitted[0].content.includes('#abcdef'), false);
});

test('semantic targets without rendered nodes do not create false DOM coverage failures', () => {
  installSerialization();
  const hiddenCatalog = structuredClone(catalog);
  hiddenCatalog.items[0].biologicalFeatures.push({
    recordKey: 'record-a', biologicalFeatureId: 'hidden-feature', hidden: true
  });
  const hiddenTargetKey = JSON.stringify(['result-a', 'record-a', 'hidden-feature']);
  const admitted = ingestSvgResults([{ name: 'a.svg', content: '<svg/>' }], {
    sanitizer: { sanitize: () => '<svg/>' },
    parser: class {
      parseFromString() { return { documentElement: root(), querySelector: () => null }; }
    }
  });
  const projection = prepareFeatureFillResultProjection({
    results: admitted,
    catalog: hiddenCatalog,
    fillsByTargetKey: {
      [targetKey]: { color: '#abcdef', caption: 'Core' },
      [hiddenTargetKey]: { color: '#abcdef', caption: 'Core' }
    },
    affectedResultKeys: ['result-a'],
    preparedSvgByResultKey: new Map([['result-a', root()]]),
    targetFeatureKeysByResult: new Map([['result-a', [targetKey, hiddenTargetKey]]])
  });
  assert.equal(projection.nextResults[0].content.includes('#abcdef'), true);
});

test('canonical Result target keys keep same biological pairs Result-local', () => {
  installSerialization();
  const duplicateCatalog = { items: [0, 1].map((resultIndex) => ({
    resultKey: `result-${resultIndex}`,
    resultIndex,
    resultName: `${resultIndex}.svg`,
    biologicalFeatures: [{ recordKey: 'shared-record', biologicalFeatureId: 'shared-feature' }],
    features: [{ recordKey: 'shared-record', biologicalFeatureId: 'shared-feature', svgId: 'svg-a' }]
  })) };
  const admitted = ingestSvgResults([0, 1].map((resultIndex) => ({
    name: `${resultIndex}.svg`, content: '<svg/>'
  })), {
    sanitizer: { sanitize: () => '<svg/>' },
    parser: class {
      parseFromString() { return { documentElement: root(), querySelector: () => null }; }
    }
  });
  const firstKey = JSON.stringify(['result-0', 'shared-record', 'shared-feature']);
  const secondKey = JSON.stringify(['result-1', 'shared-record', 'shared-feature']);
  const projection = prepareFeatureFillResultProjection({
    results: admitted,
    catalog: duplicateCatalog,
    fillsByTargetKey: {
      [firstKey]: { color: '#111111' },
      [secondKey]: { color: '#222222' }
    },
    affectedResultKeys: ['result-0', 'result-1'],
    preparedSvgByResultKey: new Map([
      ['result-0', root()],
      ['result-1', root()]
    ]),
    targetFeatureKeysByResult: new Map([
      ['result-0', [firstKey]],
      ['result-1', [secondKey]]
    ])
  });
  assert.equal(projection.nextResults[0].content.includes('#111111'), true);
  assert.equal(projection.nextResults[1].content.includes('#222222'), true);
});
