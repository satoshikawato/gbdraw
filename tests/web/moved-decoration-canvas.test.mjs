// Owner, 2026-10-23: 「キャンバスを適宜拡張などするのが好ましい。」 A moved Legend or
// plot title keeps its offset and the canvas grows to hold it. The canvas
// is the union of the automatic canvas and each moved item's bounds plus the
// edge padding, never smaller than the automatic canvas, computed when an
// item is placed (drop, position change, Generate and rerenders, resets,
// Undo/Redo). Loading a Session or showing a stored Result keeps the Result as
// it was saved (PD-OI-052).
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
  normalizeLegacyComposition,
  resetCompositionUserDeltas
} from '../../gbdraw/web/js/app/legend-layout/composition-actions.js';
import { captureDecorationContinuity } from '../../gbdraw/web/js/app/legend-layout/decoration-continuity.js';
import { createLegendCanvasActions } from '../../gbdraw/web/js/app/legend-layout/canvas-actions.js';
import { createDiagramDragActions } from '../../gbdraw/web/js/app/legend-layout/diagram-drag.js';
import { createLegendDragActions } from '../../gbdraw/web/js/app/legend/drag-actions.js';
import { createLegendLayout } from '../../gbdraw/web/js/app/legend-layout.js';

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
// The root's attributes in document order (the fake DOM keeps them in a Map).
const rootMarkup = (svg) => JSON.stringify([...svg.attributes.entries()]);
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
  assert.deepEqual(root(svg), AUTOMATIC, 'the automatic canvas keeps no record of the growth');
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

test('Generate fits the canvas of a carried Legend offset as the drop does', () => {
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

test('a legacy-normalized Result keeps its canvas when its Legend moves', () => {
  const svg = new DOMParser().parseFromString([
    '<svg viewBox="0 0 200 100" width="200px" height="100px">',
    '<g id="record" transform="translate(10,20)"></g>',
    '<g id="legend" transform="translate(150,5)"></g>',
    '</svg>'
  ].join(''), 'image/svg+xml').documentElement;
  const binding = normalizeLegacyComposition(svg, { legendSide: 'right', titleSide: 'none' });
  assert.equal(binding.metadata.legacyNormalized, true);
  applyCompositionUserDeltas(svg, { legend: [5, 0] });
  assert.equal(applyCanvasToSvg(svg, {}, { binding }), false);
  assert.deepEqual(root(svg), { viewBox: '0 0 200 100', width: '200px', height: '100px', base: null });
});

const ref = (value) => ({ value });
const appState = (svg, canvasPadding = {}) => {
  const state = Object.fromEntries([
    'results', 'selectedResultIndex', 'diagramElements', 'diagramElementIds', 'diagramElementOriginalTransforms',
    'diagramDragging', 'lengthBarElement', 'lengthBarOriginalTransform', 'plotTitleElement', 'plotTitleDragging',
    'plotTitleAutoTransform', 'generatedLegendPosition', 'legendDragging', 'legendOriginalTransform',
    'legendInitialTransform', 'svgContent'
  ].map((key) => [key, ref(null)]));
  for (const key of ['diagramOffset', 'diagramDragStart', 'lengthBarUserOffset', 'plotTitleDragStart',
    'plotTitleUserOffset', 'legendDragStart', 'legendCurrentOffset']) state[key] = { x: 0, y: 0 };
  Object.assign(state, {
    svgContainer: ref({ querySelector: (selector) => (selector === 'svg' ? svg : null) }),
    layoutRepositionMode: ref(true),
    zoom: ref(1),
    results: ref([{ id: 'r1' }]),
    activeDrawing: () => ({ canvasPadding })
  });
  return state;
};

test('showing a stored Result with an off-canvas offset changes nothing and commits nothing', () => {
  // A Result saved before the canvas grew, with this padding and with none.
  for (const padding of [{}, { top: 0, right: 10, bottom: 0, left: 10 }]) {
    const svg = resultSvg();
    applyCanvasToSvg(svg, padding);
    applyCompositionUserDeltas(svg, { legend: [-200, 10] });
    const stored = root(svg);
    const commits = [];
    const actions = createLegendCanvasActions({ state: appState(svg, padding), commitActiveResultEdit: (reason) => Boolean(commits.push(reason)) });
    assert.equal(actions.applyCanvasPadding(), false);
    assert.deepEqual(root(svg), stored);
    assert.deepEqual(commits, []);
  }
  // A Result saved with a record of its automatic canvas (padded once, then
  // unpadded) keeps that record when shown.
  const recorded = resultSvg();
  for (const name of ['data-original-view-box', 'data-original-width', 'data-original-height']) {
    recorded.setAttribute(name, name === 'data-original-view-box' ? '0 0 186 152' : name.endsWith('width') ? '186' : '152');
  }
  const savedMarkup = rootMarkup(recorded);
  const recordedCommits = [];
  const recordedActions = createLegendCanvasActions({ state: appState(recorded), commitActiveResultEdit: (reason) => Boolean(recordedCommits.push(reason)) });
  assert.equal(recordedActions.applyCanvasPadding(), false);
  assert.equal(rootMarkup(recorded), savedMarkup);
  assert.deepEqual(recordedCommits, []);
  // A padding edit on a Result saved before the canvas grew re-pads it and adds no growth.
  const old = resultSvg();
  applyCanvasToSvg(old, { left: 40 });
  applyCompositionUserDeltas(old, { legend: [-200, 10] });
  const oldActions = createLegendCanvasActions({ state: appState(old, { left: 10 }), commitActiveResultEdit: () => true });
  assert.equal(oldActions.applyCanvasPadding(), true);
  assert.equal(root(old).viewBox, '-10 0 196 152');
  // A Result saved with its grown canvas keeps it, and a padding edit keeps the growth.
  const grown = resultSvg();
  applyCompositionUserDeltas(grown, { legend: [-200, 10] });
  applyCanvasToSvg(grown, {});
  const commits = [];
  const state = appState(grown);
  const actions = createLegendCanvasActions({ state, commitActiveResultEdit: (reason) => Boolean(commits.push(reason)) });
  assert.equal(actions.applyCanvasPadding(), false);
  assert.deepEqual(root(grown), { viewBox: '-76 0 262 152', width: '262px', height: '152px', base: '0 0 186 152' });
  state.activeDrawing = () => ({ canvasPadding: { top: 0, right: 0, bottom: 0, left: 10 } });
  actions.applyCanvasPadding();
  assert.equal(root(grown).viewBox, '-86 0 272 152');
  assert.deepEqual(commits, ['canvas-padding']);
});

test('the plot title drop grows the canvas before it commits', async () => {
  const svg = resultSvg();
  svg.addEventListener = () => {};
  svg.removeEventListener = () => {};
  globalThis.document = { addEventListener: () => {}, removeEventListener: () => {} };
  const commits = [];
  const actions = createDiagramDragActions({
    state: appState(svg),
    commitActiveResultEdit: (reason) => Boolean(commits.push(`${reason} ${root(svg).viewBox}`))
  });
  actions.setupDiagramDrag(false);
  const title = svg.querySelector('#plot_title');
  actions.startDiagramDrag({ target: title, clientX: 0, clientY: 0, shiftKey: false, preventDefault: () => {} });
  await actions.endDiagramDrag({ clientX: 0, clientY: 150 });
  // Title box y 166..186; + 16 px edge padding ends at 202.
  assert.deepEqual(commits, ['diagram-drag 0 0 186 202']);
});

test('the Legend position reset returns the automatic canvas', () => {
  const svg = resultSvg();
  applyCompositionUserDeltas(svg, { legend: [-200, 10] });
  applyCanvasToSvg(svg, {});
  const commits = [];
  const actions = createLegendDragActions({ state: appState(svg), commitActiveResultEdit: (reason) => Boolean(commits.push(reason)) });
  actions.resetLegendPosition();
  assert.equal(legendAt(svg), 'translate(142,91)');
  assert.equal(root(svg).viewBox, '0 0 186 152');
  assert.deepEqual(commits, ['legend-position-reset']);
});

test('Undo and Redo of a Legend offset refit the canvas, and Undo restores the root byte for byte', () => {
  const svg = resultSvg();
  // Python writes integral sizes with a decimal point.
  svg.setAttribute('viewBox', '0.0 0.0 186.0 152.0');
  svg.setAttribute('width', '186.0px');
  svg.setAttribute('height', '152.0px');
  const before = rootMarkup(svg);
  const state = appState(svg);
  const layout = createLegendLayout({
    state,
    layOutLegend: () => null,
    previewRuntime: {
      commitActiveResultEdit: () => true,
      getResultIdentity: (result) => result.id,
      getActiveRuntime: () => null,
      commitResultEdit: (_index, edit) => edit(svg, { mounted: false })
    }
  });
  const identities = new Set(['r1']);
  layout.reconcileCompositionUserDeltas({ r1: { primary: [[0, 0]], legend: [-200, 10], title: [0, 0] } }, identities);
  assert.equal(root(svg).viewBox, '-76 0 262 152');
  layout.reconcileCompositionUserDeltas({ r1: { primary: [[0, 0]], legend: [0, 0], title: [0, 0] } }, identities);
  assert.equal(rootMarkup(svg), before, 'Undo writes no canvas record onto the restored Result');
  layout.reconcileCompositionUserDeltas({ r1: { primary: [[0, 0]], legend: [-200, 10], title: [0, 0] } }, identities);
  assert.deepEqual(root(svg), { viewBox: '-76 0 262 152', width: '262px', height: '152px', base: '0.0 0.0 186.0 152.0' }, 'Redo grows it again');
});
