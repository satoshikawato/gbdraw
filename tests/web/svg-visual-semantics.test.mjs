import assert from 'node:assert/strict';
import { visualSemanticsFromRoot, compareVisualSemantics, nonTargetFeatures } from './helpers/svg-visual-semantics.mjs';

// Minimal DOM-shaped fixture exercises the same projection as browser SVGs.
const element = (localName, attrs = {}, children = [], textContent = '') => {
  const node = { localName, textContent, getAttribute: key => attrs[key] ?? null,
    querySelectorAll: () => children.flatMap(child => [child, ...child.querySelectorAll()])
      .filter(child => child.localName !== 'g') };
  children.forEach(child => { child.parentElement = node; });
  return node;
};
const fixture = ({ fill = '#d3d3d3', depth = null, legend = 'translate(10,20)', scale = 'translate(0,100)', parts = 3, unrelated = '#d3d3d3', width = '1000', ancestorStyle = null } = {}) =>
  visualSemanticsFromRoot(element('svg', { width, height: '800', viewBox: '0 0 1000 800' }, [
    element('g', { transform: 'translate(1,2)', style: ancestorStyle }, [
      ...Array.from({ length: parts }, () => element('path', { 'data-gbdraw-feature-id': 'target', 'data-gbdraw-feature-part': 'block', d: 'M0 0L1 1', fill })),
      element('path', { 'data-gbdraw-feature-id': 'other', fill: unrelated })]),
    element('g', { display: depth }, [element('path', { d: 'M0 0L10 10', fill: 'blue' })]),
    element('g', { transform: legend }, [element('text', {}, [], 'CDS')]),
    element('g', { transform: scale }, [element('line', { x1: '0', x2: '100' }), element('text', {}, [], '5 kbp')])
  ]));
const base = fixture();
assert.deepEqual(compareVisualSemantics(base, fixture()), []);
for (const change of [{ fill: '#ff0000' }, { depth: 'none' }, { legend: 'translate(11,20.5)' },
  { scale: 'translate(0,144.5)' }, { parts: 2 }, { parts: 4 }, { width: '1001' },
  { ancestorStyle: 'visibility: hidden; fill: red' }]) {
  assert.ok(compareVisualSemantics(base, fixture(change)).length, JSON.stringify(change));
}
const wrong = fixture({ unrelated: '#cccccc' });
assert.deepEqual(compareVisualSemantics(wrong, wrong), []);
assert.ok(compareVisualSemantics(nonTargetFeatures(base, ['target']), nonTargetFeatures(wrong, ['target'])).length);
assert.deepEqual(compareVisualSemantics(base, fixture({ legend: 'translate(10.0,2e1)' })), []);
const enriched = structuredClone(base);
enriched.at(-1).label = 'known';
assert.ok(compareVisualSemantics(base, enriched).length);
assert.deepEqual(compareVisualSemantics(base, enriched, { catalogIds: ['known'] }), []);
assert.ok(compareVisualSemantics(enriched, base, { catalogIds: ['known'] }).length);
console.log('visual comparator controls passed');
