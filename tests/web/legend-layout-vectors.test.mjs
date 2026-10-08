// The JavaScript Legend layout port equals Python on every shared vector,
// exactly (zero shift): text measurement, both modes' layouts, the Legend's
// local bounds, and the composition. tests/fixtures/legend_layout_vectors.json
// holds Python's results (tools/generate_legend_layout_vectors.py;
// tests/test_legend_layout_vectors.py checks Python against them).
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import test from 'node:test';

import {
  buildCircularLegendLayout,
  buildLinearLegendLayout,
  circularLegendLocalBounds,
  createLegendFontMetrics,
  linearLegendLocalBounds,
  loadLegendFontMetrics,
  measureTextBox,
  planLegendComposition,
  pythonFloatSum,
  resolveBundledFontFace
} from '../../gbdraw/web/js/services/legend-layout.js';
import { LEGEND_FONT_METRICS } from '../../gbdraw/web/js/utils/legend-font-metrics.generated.js';

const vectors = JSON.parse(readFileSync(new URL('../fixtures/legend_layout_vectors.json', import.meta.url), 'utf8'));
const metrics = createLegendFontMetrics(LEGEND_FONT_METRICS);

// Exact equality of every number, string, boolean, and null; keys must match.
const assertExact = (actual, expected, path) => {
  if (typeof expected === 'number') {
    assert.equal(typeof actual, 'number', `${path} is not a number`);
    assert.ok(actual === expected, `${path}: Python ${expected}, JavaScript ${actual} (difference ${actual - expected})`);
    return;
  }
  if (Array.isArray(expected)) {
    assert.ok(Array.isArray(actual), `${path} is not an array`);
    assert.equal(actual.length, expected.length, `${path} length`);
    expected.forEach((value, index) => assertExact(actual[index], value, `${path}[${index}]`));
    return;
  }
  if (expected && typeof expected === 'object') {
    assert.ok(actual && typeof actual === 'object', `${path} is not an object`);
    assert.deepEqual(Object.keys(actual).sort(), Object.keys(expected).sort(), `${path} keys`);
    Object.keys(expected).forEach((key) => assertExact(actual[key], expected[key], `${path}.${key}`));
    return;
  }
  assert.equal(actual, expected, path);
};

const box = (values) => ({ minX: values[0], minY: values[1], maxX: values[2], maxY: values[3] });
const boxValues = (value) => [value.minX, value.minY, value.maxX, value.maxY];
const entries = (list) => list.map(({ key, rectX, rectY, textX, textY }) => ({ key, rectX, rectY, textX, textY }));

const linearResult = (layout) => {
  const orientation = (value) => ({
    feature: { entries: entries(value.feature.entries), width: value.feature.width, height: value.feature.height, numLines: value.feature.numLines },
    gradient: value.gradient && {
      compact: value.gradient.compact,
      entries: value.gradient.entries.map(({ key, titleX, titleY, barX, barY }) => ({ key, titleX, titleY, barX, barY })),
      width: value.gradient.width,
      height: value.gradient.height,
      barWidth: value.gradient.barWidth,
      minLabelText: value.gradient.minLabelText,
      minLabelX: value.gradient.minLabelX,
      maxLabelX: value.gradient.maxLabelX,
      scaleLabelY: value.gradient.scaleLabelY
    },
    featureX: value.featureX,
    featureY: value.featureY,
    gradientX: value.gradientX,
    gradientY: value.gradientY,
    width: value.width,
    height: value.height
  });
  return { horizontal: orientation(layout.horizontal), vertical: orientation(layout.vertical), activeOrientation: layout.activeOrientation };
};

const circularResult = (layout) => ({
  horizontal: layout.horizontal,
  width: layout.width,
  height: layout.height,
  featureWidth: layout.featureWidth,
  featureHeight: layout.featureHeight,
  pairwiseLegendWidth: layout.pairwiseLegendWidth,
  lineMargin: layout.lineMargin,
  xMargin: layout.xMargin,
  numLines: layout.numLines,
  numColumns: layout.numColumns,
  numItemsPerLine: layout.numItemsPerLine,
  entries: entries(layout.entries),
  gradient: layout.gradient && {
    compact: layout.gradient.compact,
    width: layout.gradient.width,
    height: layout.gradient.height,
    barWidth: layout.gradient.barWidth,
    barX: layout.gradient.barX,
    minLabelText: layout.gradient.minLabelText,
    scaleY: layout.gradient.scaleY,
    compactEntries: layout.gradient.compactEntries.map(({ key, labelY, barY }) => ({ key, labelY, barY })),
    singleEntries: layout.gradient.singleEntries.map(({ key, titleX, titleY, barX, barY, minLabelX, maxLabelX, scaleLabelY }) => (
      { key, titleX, titleY, barX, barY, minLabelX, maxLabelX, scaleLabelY }
    ))
  },
  gradientX: layout.gradientX,
  gradientY: layout.gradientY
});

const layoutCase = (vector) => {
  const options = { ...vector.options };
  const result = vector.mode === 'linear'
    ? linearResult(buildLinearLegendLayout(vector.rows, options, metrics))
    : circularResult(buildCircularLegendLayout(vector.rows, options, metrics));
  const legend = vector.mode === 'linear'
    ? linearLegendLocalBounds(buildLinearLegendLayout(vector.rows, options, metrics), options.colorRectSize)
    : circularLegendLocalBounds(buildCircularLegendLayout(vector.rows, options, metrics), options.colorRectSize);
  const inputs = vector.composition;
  const plan = planLegendComposition({
    primary: box(inputs.primary),
    legend,
    title: inputs.title ? box(inputs.title) : null,
    legendSide: options.side,
    titleSide: inputs.titleSide,
    overlayObstacles: inputs.overlayObstacles.map(box),
    spacing: inputs.spacing,
    overlayPolicy: vectors.overlayPolicy
  });
  return {
    ...result,
    localBounds: boxValues(legend),
    composition: {
      canvas: boxValues(plan.canvas),
      placements: plan.placements.map(({ role, translation, finalBounds }) => ({ role, translation, finalBounds: boxValues(finalBounds) })),
      overlayObstacles: plan.overlayObstacles.map(boxValues),
      overlayConflictIndices: plan.overlayConflictIndices,
      overlayResolution: plan.overlayResolution
    }
  };
};

test('the vectors cover every Gallery Legend, its edits, and the synthetic branches', () => {
  const gallery = vectors.layout.filter((vector) => !vector.name.startsWith('synthetic'));
  assert.equal(gallery.filter((vector) => vector.edit === 'none').length, 10);
  assert.ok(gallery.filter((vector) => vector.edit !== 'none').length >= 80);
  assert.deepEqual(new Set(vectors.layout.map((vector) => vector.mode)), new Set(['circular', 'linear']));
  const branches = vectors.layout.map((vector) => vector.expected);
  assert.ok(branches.some((expected) => expected.gradient?.compact === false), 'a single Circular gradient');
  assert.ok(branches.some((expected) => expected.gradient?.compact === true), 'a compact Circular gradient');
  assert.ok(branches.some((expected) => expected.horizontal?.gradient?.compact === true), 'a compact Linear gradient');
  assert.ok(branches.some((expected) => expected.composition.overlayResolution === 'canvas_growth'), 'an overlay that grows the canvas');
  assert.ok(vectors.measurement.length >= 200);
});

test('text measurement equals calculate_bbox_dimensions', () => {
  vectors.measurement.forEach((vector, index) => {
    assert.equal(resolveBundledFontFace(metrics, vector.fontFamily), vector.fontFile, `measurement[${index}] face`);
    const measured = measureTextBox(metrics, vector);
    assertExact(measured, vector.expected, `measurement[${index}] ${JSON.stringify(vector.text)} ${vector.fontFile} ${vector.fontSize}`);
  });
});

test('Legend layouts, local bounds, and composition equal Python on every vector', () => {
  vectors.layout.forEach((vector, index) => {
    assertExact(layoutCase(vector), vector.expected, `layout[${index}] ${vector.name} / ${vector.edit}`);
  });
});

test('the lazy loader indexes the same table', async () => {
  const loaded = await loadLegendFontMetrics();
  assert.equal(loaded, await loadLegendFontMetrics(), 'loaded once');
  const vector = vectors.measurement.find((entry) => entry.text === 'AVATAR' && entry.fontFile === 'LiberationSans-Regular');
  assertExact(measureTextBox(loaded, vector), vector.expected, 'loaded measurement');
});

test('sums follow CPython float sum()', () => {
  assert.equal(pythonFloatSum([]), 0);
  assert.equal(pythonFloatSum(Array(10).fill(0.1)), 1.0);
  assert.equal(pythonFloatSum([1e16, 1.0, -1e16]), 1.0);
});
