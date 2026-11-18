import assert from 'node:assert/strict';
import { test } from 'node:test';
import { installFakeSvgDom } from './fake-svg-dom.mjs';
import { drawnBlockStroke } from '../../gbdraw/web/js/services/legend-svg.js';
import { createEmptySvgMutationPlan, reconcileMountedResult } from '../../gbdraw/web/js/services/svg-result-ingestion.js';

// The first feature of a record as Python draws it, a `join(101..200,301..400)`
// CDS with the default strokes: its connector line (lightgray, 5) comes before
// its two blocks (gray, 2) in both modes (connectors are drawn under blocks).
const FIRST_FEATURE = {
  linear: [
    '<path d="M 40,0 L60,0 z" data-gbdraw-feature-id="f22f5f5cf" data-gbdraw-feature-part="connector" fill="none" id="f22f5f5cf__line1" stroke="lightgray" stroke-width="5.0" />',
    '<path d="M 20,-20 L 40,-20 L 40,20 L 20,20 z" data-gbdraw-feature-id="f22f5f5cf" data-gbdraw-feature-part="block" fill="#54bcf8" id="f22f5f5cf__part1" stroke="gray" stroke-width="2.0" />',
    '<path d="M 60,-20 L 64,-20 L 80,0 L 64,20 L 60,20 z" data-gbdraw-feature-id="f22f5f5cf" data-gbdraw-feature-part="block" fill="#54bcf8" id="f22f5f5cf__part2" stroke="gray" stroke-width="2.0" />'
  ],
  circular: [
    '<path d="M 46.069,-364.677 A 367.575,367.575 0 0 1 68.877,-361.064" data-gbdraw-feature-id="f22f5f5cf" data-gbdraw-feature-part="connector" fill="none" id="f22f5f5cf__line1" stroke="lightgray" stroke-width="5.0" />',
    '<path d="M 21.917,-348.361 A 349.05,349.05 0 0 1 43.748,-346.298 L 48.391,-383.055 A 386.1,386.1 0 0 0 24.243,-385.338 z" data-gbdraw-feature-id="f22f5f5cf" data-gbdraw-feature-part="block" fill="#54bcf8" id="f22f5f5cf__part1" stroke="gray" stroke-width="2.0" />',
    '<path d="M 65.405,-342.867 A 349.05,349.05 0 0 1 73.621,-341.198 L 91.412,-356.027 L 81.435,-377.414 A 386.1,386.1 0 0 0 72.348,-379.261 z" data-gbdraw-feature-id="f22f5f5cf" data-gbdraw-feature-part="block" fill="#54bcf8" id="f22f5f5cf__part2" stroke="gray" stroke-width="2.0" />'
  ]
};

for (const [mode, elements] of Object.entries(FIRST_FEATURE)) {
  test(`the Result's block stroke is a block's, not the first feature's connector line (${mode})`, () => {
    installFakeSvgDom();
    const svg = new DOMParser().parseFromString(`<svg><g id="features">${elements.join('')}</g></svg>`).documentElement;
    assert.deepEqual(drawnBlockStroke(svg), { color: 'gray', width: 2 });
  });
}

// OV-3725: a partial first feature (`<101..200`) as Python draws it: its fill
// path without a stroke, then an `__outline` path, open at the partial end,
// that carries the block stroke.
const PARTIAL_FIRST_FEATURE = {
  linear: [
    '<path d="M 20,-20 L 24,-20 L 40,0 L 24,20 L 20,20 z" data-gbdraw-feature-id="f6f815dab" data-gbdraw-feature-part="block" fill="#54bcf8" id="f6f815dab" stroke="none" stroke-width="2.0" />',
    '<path d="M 20,-20 L 24,-20 L 40,0 L 24,20 L 20,20 " data-gbdraw-feature-id="f6f815dab" data-gbdraw-feature-part="block" fill="none" id="f6f815dab__outline" stroke="gray" stroke-width="2.0" />'
  ],
  circular: [
    '<path d="M 21.917,-348.361 A 349.05,349.05 0 0 1 30.277,-347.734 L 46.069,-364.677 L 33.49,-384.645 A 386.1,386.1 0 0 0 24.243,-385.338 z" data-gbdraw-feature-id="f6f815dab" data-gbdraw-feature-part="block" fill="#54bcf8" id="f6f815dab" stroke="none" stroke-width="2.0" />',
    '<path d="M 21.917,-348.361 A 349.05,349.05 0 0 1 30.277,-347.734 L 46.069,-364.677 L 33.49,-384.645 A 386.1,386.1 0 0 0 24.243,-385.338 " data-gbdraw-feature-id="f6f815dab" data-gbdraw-feature-part="block" fill="none" id="f6f815dab__outline" stroke="gray" stroke-width="2.0" />'
  ]
};

for (const [mode, elements] of Object.entries(PARTIAL_FIRST_FEATURE)) {
  test(`OV-3725: the Result's block stroke of a partial first feature is its outline's (${mode})`, () => {
    installFakeSvgDom();
    const svg = new DOMParser().parseFromString(`<svg><g id="features">${elements.join('')}</g></svg>`).documentElement;
    assert.deepEqual(drawnBlockStroke(svg), { color: 'gray', width: 2 });
  });
}

// OV-3726: a `join(<101..200,301..400)` CDS as Python draws it: its connector
// line, its partial first block (a fill path without a stroke, then the
// outline open at the partial end), and its second block.
const PARTIAL_JOINED_FEATURE = {
  linear: [
    '<path d="M 40,0 L60,0 z" data-gbdraw-feature-id="f22f5f5cf" data-gbdraw-feature-part="connector" fill="none" id="f22f5f5cf__line1" stroke="lightgray" stroke-width="5.0" />',
    '<path d="M 20,-20 L 40,-20 L 40,20 L 20,20 z" data-gbdraw-feature-id="f22f5f5cf" data-gbdraw-feature-part="block" fill="#54bcf8" id="f22f5f5cf__part1" stroke="none" stroke-width="2.0" />',
    '<path d="M 20,-20 L 40,-20 L 40,20 L 20,20 " data-gbdraw-feature-id="f22f5f5cf" data-gbdraw-feature-part="block" fill="none" id="f22f5f5cf__part1__outline" stroke="gray" stroke-width="2.0" />',
    '<path d="M 60,-20 L 64,-20 L 80,0 L 64,20 L 60,20 z" data-gbdraw-feature-id="f22f5f5cf" data-gbdraw-feature-part="block" fill="#54bcf8" id="f22f5f5cf__part2" stroke="gray" stroke-width="2.0" />'
  ],
  circular: [
    '<path d="M 46.069,-364.677 A 367.575,367.575 0 0 1 68.877,-361.064" data-gbdraw-feature-id="f22f5f5cf" data-gbdraw-feature-part="connector" fill="none" id="f22f5f5cf__line1" stroke="lightgray" stroke-width="5.0" />',
    '<path d="M 21.917,-348.361 A 349.05,349.05 0 0 1 43.748,-346.298 L 48.391,-383.055 A 386.1,386.1 0 0 0 24.243,-385.338 z" data-gbdraw-feature-id="f22f5f5cf" data-gbdraw-feature-part="block" fill="#54bcf8" id="f22f5f5cf__part1" stroke="none" stroke-width="2.0" />',
    '<path d="M 21.917,-348.361 A 349.05,349.05 0 0 1 43.748,-346.298 L 48.391,-383.055 A 386.1,386.1 0 0 0 24.243,-385.338 " data-gbdraw-feature-id="f22f5f5cf" data-gbdraw-feature-part="block" fill="none" id="f22f5f5cf__part1__outline" stroke="gray" stroke-width="2.0" />',
    '<path d="M 65.405,-342.867 A 349.05,349.05 0 0 1 73.621,-341.198 L 91.412,-356.027 L 81.435,-377.414 A 386.1,386.1 0 0 0 72.348,-379.261 z" data-gbdraw-feature-id="f22f5f5cf" data-gbdraw-feature-part="block" fill="#54bcf8" id="f22f5f5cf__part2" stroke="gray" stroke-width="2.0" />'
  ]
};

const LEGEND_ROW = '<g id="legend"><g id="feature_legend"><g data-legend-key="CDS"><path fill="#54bcf8" stroke="gray" stroke-width="2.0" /></g></g></g>';

/** The stroke of each feature element by id: color and width. */
const featureStrokes = (svg) => Object.fromEntries(Array.from(svg.querySelectorAll('[data-gbdraw-feature-id]'))
  .map((element) => [element.getAttribute('id'), `${element.getAttribute('stroke')}/${element.getAttribute('stroke-width')}`]));

// The edit strokes what carries the feature's stroke as Python draws it: the
// partial block's outline, open at the partial end, takes the edit and its
// fill path keeps no stroke, so the end stays open. The connector line keeps
// Python's line stroke (a separate control; Owner 2026-11-18).
const EDITED = {
  f22f5f5cf__line1: 'lightgray/5.0',
  f22f5f5cf__part1: 'none/2.0',
  f22f5f5cf__part1__outline: '#222222/3',
  f22f5f5cf__part2: '#222222/3'
};

for (const [mode, elements] of Object.entries(PARTIAL_JOINED_FEATURE)) {
  const mounted = () => {
    installFakeSvgDom();
    return new DOMParser().parseFromString(`<svg><g id="features">${elements.join('')}</g>${LEGEND_ROW}</svg>`).documentElement;
  };

  test(`OV-3726: a feature stroke edit keeps a partial block's end open (${mode})`, () => {
    const svg = mounted();
    reconcileMountedResult(svg, {
      ...createEmptySvgMutationPlan(1).operationsByResult[0],
      featureStrokes: [{ renderedId: 'f22f5f5cf', strokeColor: '#222222', strokeWidth: 3 }]
    });
    assert.deepEqual(featureStrokes(svg), EDITED);
  });

  test(`OV-3726: a Legend row stroke keeps a partial block's end open (${mode})`, () => {
    const svg = mounted();
    reconcileMountedResult(svg, {
      ...createEmptySvgMutationPlan(1).operationsByResult[0],
      legendStrokes: [{
        caption: 'CDS', strokeColor: '#222222', strokeWidth: 3,
        reach: { listedIds: [], namedIds: [], ownStrokeIds: [], draftColor: null }
      }]
    });
    assert.deepEqual(featureStrokes(svg), EDITED);
  });
}

// An automatic feature underlay is drawn behind the glyphs, so before them,
// without a stroke: the Result's block stroke is a glyph block's.
test('the Result\'s block stroke is not an automatic underlay\'s', () => {
  installFakeSvgDom();
  const svg = new DOMParser().parseFromString(`<svg><g id="features">${[
    '<path d="M 0,-30 L 100,-30 L 100,30 L 0,30 z" data-gbdraw-auto-feature-underlay="true" data-gbdraw-feature-id="f0underlay" data-gbdraw-feature-part="block" fill="#808080" id="f0underlay" stroke="none" stroke-width="0" />',
    ...FIRST_FEATURE.linear
  ].join('')}</g></svg>`).documentElement;
  assert.deepEqual(drawnBlockStroke(svg), { color: 'gray', width: 2 });
});
