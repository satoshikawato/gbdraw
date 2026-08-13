import assert from 'node:assert/strict';
import test from 'node:test';

import { applyFeatureLabelProjectionToSvg } from '../../gbdraw/web/js/app/feature-editor/label-result-projection.js';

class FakeText {
  constructor(value, attributes = {}) {
    this.textContent = value;
    this.attributes = new Map(Object.entries(attributes));
  }

  querySelector() { return null; }
  closest() { return null; }
  getAttribute(name) { return this.attributes.has(name) ? this.attributes.get(name) : null; }
  setAttribute(name, value) { this.attributes.set(name, String(value)); }
}

class FakeSvg {
  constructor(elements) { this.elements = elements; }
  querySelectorAll(selector) {
    return selector === 'text[data-label-editable="true"]' ? this.elements : [];
  }
}

test('source group updates every matching label and preserves unrelated labels', () => {
  const a = new FakeText('Kinase');
  const b = new FakeText('Previously renamed', { 'data-label-source-text': 'Kinase' });
  const c = new FakeText('Transporter', { 'data-label-source-text': 'Transporter' });
  const outcome = applyFeatureLabelProjectionToSvg({
    svg: new FakeSvg([a, b, c]),
    semanticScope: 'source-annotation-label',
    labelsByTargetKey: {
      a: { text: 'All kinases', sourceText: 'Kinase', renderedSvgIds: ['svg-a'] },
      b: { text: 'All kinases', sourceText: 'Kinase', renderedSvgIds: ['svg-b'] }
    },
    targetFeatureKeys: ['a', 'b'],
    sourceText: 'Kinase'
  });
  assert.equal(a.textContent, 'All kinases');
  assert.equal(b.textContent, 'All kinases');
  assert.equal(c.textContent, 'Transporter');
  assert.equal(a.getAttribute('data-label-source-text'), 'Kinase');
  assert.deepEqual(outcome, {
    targetFeatures: 2, coveredTargets: 2, matchedLabels: 2, changedLabels: 2
  });
});

test('rendered group uses current presentation while exact scopes require unique binding', () => {
  const renderedA = new FakeText('Shared');
  const renderedB = new FakeText('Shared');
  const other = new FakeText('Other');
  const rendered = applyFeatureLabelProjectionToSvg({
    svg: new FakeSvg([renderedA, renderedB, other]),
    semanticScope: 'rendered-label',
    labelsByTargetKey: {
      a: { text: 'Together', sourceText: 'A', renderedSvgIds: ['svg-a'] },
      b: { text: 'Together', sourceText: 'B', renderedSvgIds: ['svg-b'] }
    },
    targetFeatureKeys: ['a', 'b'],
    matchText: 'Shared'
  });
  assert.equal(renderedA.textContent, 'Together');
  assert.equal(renderedB.textContent, 'Together');
  assert.equal(other.textContent, 'Other');
  assert.equal(rendered.matchedLabels, 2);

  const exact = new FakeText('Old', { 'data-label-feature-id': 'svg-a' });
  const exactOutcome = applyFeatureLabelProjectionToSvg({
    svg: new FakeSvg([exact, new FakeText('Old')]),
    semanticScope: 'single',
    labelsByTargetKey: {
      a: { text: 'Only this', sourceText: 'Kinase', fromText: 'Old', renderedSvgIds: ['svg-a'] }
    },
    targetFeatureKeys: ['a']
  });
  assert.equal(exact.textContent, 'Only this');
  assert.equal(exactOutcome.coveredTargets, 1);

  assert.throws(() => applyFeatureLabelProjectionToSvg({
    svg: new FakeSvg([new FakeText('Same'), new FakeText('Same')]),
    semanticScope: 'single',
    labelsByTargetKey: {
      a: { text: 'Only this', sourceText: 'Same', fromText: 'Same', renderedSvgIds: ['svg-a'] }
    },
    targetFeatureKeys: ['a']
  }), /target is unavailable/);
});
