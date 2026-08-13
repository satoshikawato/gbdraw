import assert from 'node:assert/strict';
import test from 'node:test';

import {
  prepareLegendOrderResultProjection
} from '../../gbdraw/web/js/app/legend/manual-intent-command.js';
import {
  getAllFeatureLegendGroups,
  parseTransformXY
} from '../../gbdraw/web/js/app/legend/utils.js';

class MockElement {
  constructor(tagName, attributes = {}, textContent = '') {
    this.tagName = tagName;
    this.localName = tagName;
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
  cloneNode(deep = false) {
    const clone = new MockElement(
      this.tagName,
      Object.fromEntries(this.attributes),
      this.textContent
    );
    if (deep) this.children.forEach((child) => clone.appendChild(child.cloneNode(true)));
    return clone;
  }
  matchesSelector(selector) {
    if (selector === 'text' || selector === 'path') return this.tagName === selector;
    if (selector === 'g[data-legend-key]') {
      return this.tagName === 'g' && this.hasAttribute('data-legend-key');
    }
    if (selector.startsWith('#')) return this.id === selector.slice(1);
    return false;
  }
  querySelectorAll(selector) {
    const matches = [];
    const visit = (node) => node.children.forEach((child) => {
      if (child.matchesSelector(selector)) matches.push(child);
      visit(child);
    });
    visit(this);
    return matches;
  }
  querySelector(selector) { return this.querySelectorAll(selector)[0] || null; }
  getElementById(id) { return this.id === id ? this : this.querySelector(`#${id}`); }
}

const legendEntry = ({ caption, color, x }) => {
  const group = new MockElement('g', { 'data-legend-key': caption });
  group.appendChild(new MockElement('path', { fill: color, transform: `translate(${x},7)` }));
  group.appendChild(new MockElement('text', { transform: `translate(${x + 22},7)` }, caption));
  return group;
};

const resultRoot = (entries, transform) => {
  const root = new MockElement('svg');
  const outer = new MockElement('g', { id: 'legend', transform });
  const featureLegend = new MockElement('g', { id: 'feature_legend' });
  entries.forEach((entry) => featureLegend.appendChild(legendEntry(entry)));
  outer.appendChild(featureLegend);
  root.appendChild(outer);
  return root;
};

const rootState = (root) => {
  const outer = root.querySelector('#legend');
  const groups = Array.from(getAllFeatureLegendGroups(root)[0].children);
  return {
    outerTransform: outer.getAttribute('transform'),
    captions: groups.map((group) => group.getAttribute('data-legend-key')),
    slots: groups.map((group) => (
      parseTransformXY(group.querySelector('path').getAttribute('transform')).x
    ))
  };
};

const rootContent = (root) => JSON.stringify(rootState(root));

test('order projection stages every heterogeneous Result and retains its local layout', async () => {
  const roots = [
    resultRoot([
      { caption: 'Alpha', color: '#111111', x: 0 },
      { caption: 'Beta', color: '#222222', x: 70 },
      { caption: 'Only A', color: '#333333', x: 140 }
    ], 'translate(500,20)'),
    resultRoot([
      { caption: 'Gamma', color: '#444444', x: 5 },
      { caption: 'Beta', color: '#222222', x: 80 }
    ], 'translate(30,900)')
  ];
  const results = roots.map((root, index) => ({
    name: index ? 'b.svg' : 'a.svg',
    content: rootContent(root),
    root
  }));
  const catalog = { schema: 3, items: [{
    resultKey: 'result-a', resultIndex: 0, resultName: 'a.svg'
  }, {
    resultKey: 'result-b', resultIndex: 1, resultName: 'b.svg'
  }] };
  const projection = await prepareLegendOrderResultProjection({
    results,
    catalog,
    legendOrderIntent: ['Beta', 'Alpha', 'Only A'],
    selectedResultIndex: 0,
    selectedLegendEntries: [
      { caption: 'Alpha', color: '#111111', xPos: 0, yPos: 7 },
      { caption: 'Beta', color: '#222222', xPos: 70, yPos: 7 },
      { caption: 'Only A', color: '#333333', xPos: 140, yPos: 7 }
    ],
    parseResultRoot: (result) => result.root,
    admitResult: (result, root) => ({ ...result, root, content: rootContent(root) }),
    committedContent: (result) => result.content,
    committedMetadata: () => null
  });

  assert.notEqual(projection.nextResults, results);
  assert.deepEqual(projection.affectedResultKeys, ['result-a', 'result-b']);
  assert.deepEqual(projection.sourceOrderByResult, {
    'result-a': ['Alpha', 'Beta', 'Only A'],
    'result-b': ['Gamma', 'Beta']
  });
  assert.deepEqual(rootState(projection.nextResults[0].root), {
    outerTransform: 'translate(500,20)',
    captions: ['Beta', 'Alpha', 'Only A'],
    slots: [0, 70, 140]
  });
  assert.deepEqual(rootState(projection.nextResults[1].root), {
    outerTransform: 'translate(30,900)',
    captions: ['Beta', 'Gamma'],
    slots: [5, 80]
  });
  assert.deepEqual(
    projection.selectedLegendEntries.map((entry) => [entry.caption, entry.xPos, entry.yPos]),
    [['Beta', 0, 7], ['Alpha', 70, 7], ['Only A', 140, 7]]
  );
  assert.equal(projection.counters.affectedResults, 2);
  assert.equal(projection.counters.detachedPasses, 2);
  assert.equal(projection.counters.mountedSerializations, 0);
});
