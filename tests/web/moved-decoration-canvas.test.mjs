// Owner, 2026-10-23: 「キャンバスを適宜拡張などするのが好ましい。」 A moved Legend or
// plot title keeps its offset and the canvas grows to hold it. The canvas
// is the union of the automatic canvas and each moved item's bounds plus the
// edge padding, never smaller than the automatic canvas, computed when the
// item is placed (drop, position change, Generate, Session load).
import assert from 'node:assert/strict';
import test from 'node:test';

import { installFakeSvgDom } from './fake-svg-dom.mjs';
import {
  applyCanvasToSvg,
  applyCompositionEdit,
  applyCompositionUserDeltas,
  COMPOSITION_METADATA_ATTRIBUTE,
  COMPOSITION_SCHEMA_ATTRIBUTE,
  compositionUserDeltas,
  resetCompositionUserDeltas
} from '../../gbdraw/web/js/app/legend-layout/composition-actions.js';
import { captureDecorationContinuity } from '../../gbdraw/web/js/app/legend-layout/decoration-continuity.js';

installFakeSvgDom();

const spacing = { dockGapPx: 24, edgePaddingPx: 16, overlayClearancePx: 8, stackGapPx: 20, titleGapPx: 20 };
const overlayPolicy = {
  candidateScoreOrder: ['totalAnchorDistance', 'xAnchorDistance', 'yAnchorDistance', 'nearEdgeX', 'nearEdgeY'],
  canvasGrowthCandidateOrder: ['horizontal', 'vertical'],
  canvasGrowthScoreOrder: ['addedArea', 'addedExtent', 'candidateOrder'],
  quadrantBoundaryRatio: 0.5
};
const target = (role, automaticTranslation, boundsKey, bounds) => ({
  automaticTranslation, [boundsKey]: bounds, role, selector: `[data-gbdraw-composition-role="${role}"]`
});

// The automatic layout: primary 100 x 80, Legend on the right, title on top,
// 16 px edge padding: the automatic canvas is 186 x 152.
const resultSvg = ({ legendSide = 'right', titleSide = 'top' } = {}) => {
  const svg = new DOMParser().parseFromString([
    '<svg viewBox="0 0 186 152" width="186px" height="152px">',
    '<g data-gbdraw-composition-role="primary" id="record_a" transform="translate(6,36)"></g>',
    '<g data-gbdraw-composition-role="legend" id="legend" transform="translate(142,91)"></g>',
    '<g data-gbdraw-composition-role="title" id="plot_title" transform="translate(41,26) scale(1)"></g>',
    '</svg>'
  ].join(''), 'image/svg+xml').documentElement;
  svg.setAttribute(COMPOSITION_SCHEMA_ATTRIBUTE, '1');
  svg.setAttribute(COMPOSITION_METADATA_ATTRIBUTE, JSON.stringify({
    legend: target('legend', [142, 91], 'localBounds', { x: -2, y: -5, width: 30, height: 20 }),
    legendReflow: { colorRectSize: 14, lineHeight: 24, textXOffset: 22 },
    legendSide,
    overlayObstacles: [],
    overlayPolicy,
    primary: target('primary', [6, 36], 'finalBounds', { x: 16, y: 56, width: 100, height: 80 }),
    spacing,
    title: target('title', [41, 26], 'localBounds', { x: -5, y: -10, width: 60, height: 20 }),
    titleSide
  }));
  return svg;
};
const root = (svg) => ({
  viewBox: svg.getAttribute('viewBox'), width: svg.getAttribute('width'), height: svg.getAttribute('height'),
  base: svg.getAttribute('data-original-view-box')
});
const AUTOMATIC = { viewBox: '0 0 186 152', width: '186px', height: '152px', base: null };
const legendAt = (svg) => svg.querySelector('#legend').getAttribute('transform');

test('a Legend moved past the left edge grows the canvas and keeps its offset; reset returns the automatic canvas', () => {
  const svg = resultSvg();
  applyCompositionUserDeltas(svg, { legend: [-200, 10] });
  assert.equal(applyCanvasToSvg(svg, {}), true);
  // Legend box x -60..-30, y 96..116; with the 16 px edge padding the canvas starts at x -76.
  assert.deepEqual(root(svg), { viewBox: '-76 0 262 152', width: '262px', height: '152px', base: '0 0 186 152' });
  assert.equal(legendAt(svg), 'translate(-58,101)');
  assert.deepEqual(compositionUserDeltas(svg).legend, [-200, 10]);
  assert.equal(applyCanvasToSvg(svg, {}), false, 'the same placement changes nothing');

  resetCompositionUserDeltas(svg);
  assert.equal(applyCanvasToSvg(svg, {}), true);
  assert.deepEqual(root(svg), { ...AUTOMATIC, base: '0 0 186 152' });
});

test('a drag that stays inside the automatic canvas leaves the root byte-identical', () => {
  const svg = resultSvg();
  applyCompositionUserDeltas(svg, { legend: [-20, -30], title: [10, 5] });
  assert.equal(applyCanvasToSvg(svg, {}), false);
  assert.deepEqual(root(svg), AUTOMATIC);
});

test('a title moved below the canvas grows its height, and the canvas padding goes around the grown canvas', () => {
  const svg = resultSvg();
  applyCompositionUserDeltas(svg, { title: [0, 150] });
  // Title box y 166..186; + 16 px edge padding ends at 202.
  applyCanvasToSvg(svg, { top: 0, right: 10, bottom: 20, left: 0 });
  assert.deepEqual(root(svg), { viewBox: '0 0 196 222', width: '196px', height: '222px', base: '0 0 186 152' });
});

test('a hidden Legend does not grow the canvas', () => {
  const svg = resultSvg({ legendSide: 'none' });
  applyCompositionUserDeltas(svg, { legend: [-200, 10] });
  assert.equal(applyCanvasToSvg(svg, {}), false);
  assert.deepEqual(root(svg), AUTOMATIC);
});

test('OV-1631: a Legend Position change keeps the drag offset and the canvas holds the Legend', () => {
  const svg = resultSvg();
  applyCompositionUserDeltas(svg, { legend: [-134, 67] });
  applyCompositionEdit(svg, { legendSide: 'top', canvasPadding: {} });
  const [x, y, width, height] = svg.getAttribute('viewBox').split(' ').map(Number);
  const legend = svg.querySelector('#legend');
  const [lx, ly] = legend.getAttribute('transform').match(/translate\(([^,]+),([^)]+)\)/).slice(1).map(Number);
  assert.deepEqual(compositionUserDeltas(svg).legend, [-134, 67], 'the offset is kept');
  // Legend local box (-2,-5)..(28,15) at its translate, with the 16 px edge padding.
  assert.ok(x <= lx - 2 - 16 && y <= ly - 5 - 16, `${x},${y} holds the Legend at ${lx},${ly}`);
  assert.ok(x + width >= lx + 28 + 16 && y + height >= ly + 15 + 16);
  assert.ok(x < 0, 'the Legend left of the automatic canvas grows it to the left');
  assert.equal(svg.getAttribute('width'), `${width}px`);
  assert.equal(svg.getAttribute('height'), `${height}px`);
});

test('Generate gives the moved Legend the same canvas as the in-place edit (live == Generate)', () => {
  const live = resultSvg();
  applyCompositionUserDeltas(live, { legend: [-200, 10] });
  applyCanvasToSvg(live, {});
  const plan = captureDecorationContinuity({
    canonical: null,
    results: [{ content: '' }],
    catalog: { items: [{ recordKeys: ['a'] }] },
    mountedSvg: live,
    projectRecordIdentity: () => ({ record: 'a' }),
    canvasPadding: null
  });
  const candidate = resultSvg();
  plan(null, { catalog: { items: [{ recordKeys: ['a'] }] } }).forEach((transform) => transform?.(candidate));
  assert.deepEqual(root(candidate), root(live));
  assert.equal(root(candidate).viewBox, '-76 0 262 152', 'Generate grows the canvas without canvas padding');
  assert.equal(legendAt(candidate), legendAt(live));
});
