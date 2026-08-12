import assert from 'node:assert/strict';

import {
  applyCompositionEdit,
  applyCompositionUserDeltas,
  bindCompositionMetadata,
  COMPOSITION_METADATA_ATTRIBUTE,
  COMPOSITION_ROLE_ATTRIBUTE,
  COMPOSITION_SCHEMA_ATTRIBUTE,
  CompositionMetadataError,
  compositionUserDeltas,
  normalizeLegacyComposition,
  parseCompositionMetadata,
  planComposition,
  reconcileCompositionTitle,
  resetCompositionUserDeltas
} from '../../gbdraw/web/js/app/legend-layout/composition-actions.js';
import {
  prependTranslate,
  readLeadingTranslate,
  replaceLeadingTranslate
} from '../../gbdraw/web/js/app/legend-layout/transform-utils.js';

class FakeElement {
  constructor({ tagName = 'g', id = '', attributes = {}, bbox = { x: 0, y: 0, width: 1, height: 1 } } = {}) {
    this.tagName = tagName;
    this.id = id;
    this.attributes = new Map(Object.entries(attributes));
    if (id) this.attributes.set('id', id);
    this.children = [];
    this.bbox = { ...bbox };
    this.style = {};
  }

  appendChild(child) {
    this.children.push(child);
    return child;
  }

  getAttribute(name) {
    return this.attributes.has(name) ? this.attributes.get(name) : null;
  }

  setAttribute(name, value) {
    this.attributes.set(name, String(value));
    if (name === 'id') this.id = String(value);
  }

  removeAttribute(name) {
    this.attributes.delete(name);
  }

  getBBox() {
    return { ...this.bbox };
  }

  querySelectorAll(selector) {
    const roleMatch = selector.match(/^\[data-gbdraw-composition-role="([^"]+)"\]$/);
    const matches = [];
    const visit = (node) => {
      if (roleMatch && node.getAttribute(COMPOSITION_ROLE_ATTRIBUTE) === roleMatch[1]) matches.push(node);
      node.children.forEach(visit);
    };
    this.children.forEach(visit);
    return matches;
  }

  getElementById(id) {
    let match = null;
    const visit = (node) => {
      if (match) return;
      if (node.id === id) match = node;
      node.children.forEach(visit);
    };
    this.children.forEach(visit);
    return match;
  }
}

const targetPayload = (role, automaticTranslation, boundsKey, bounds) => ({
  automaticTranslation,
  [boundsKey]: bounds,
  role,
  selector: `[data-gbdraw-composition-role="${role}"]`
});
const compositionSpacing = {
  dockGapPx: 24,
  edgePaddingPx: 16,
  overlayClearancePx: 8,
  stackGapPx: 20,
  titleGapPx: 20
};
const overlayPolicy = {
  candidateScoreOrder: [
    'totalAnchorDistance',
    'xAnchorDistance',
    'yAnchorDistance',
    'nearEdgeX',
    'nearEdgeY'
  ],
  canvasGrowthCandidateOrder: ['horizontal', 'vertical'],
  canvasGrowthScoreOrder: ['addedArea', 'addedExtent', 'candidateOrder'],
  quadrantBoundaryRatio: 0.5
};
const legendReflow = {
  colorRectSize: 14,
  lineHeight: 24,
  textXOffset: 22
};

const schemaOneSvg = () => {
  const metadata = {
    legend: targetPayload('legend', [142, 91], 'localBounds', { x: -2, y: -5, width: 30, height: 20 }),
    legendReflow,
    legendSide: 'right',
    overlayObstacles: [],
    overlayPolicy,
    primary: targetPayload('primary', [6, 36], 'finalBounds', { x: 16, y: 56, width: 100, height: 80 }),
    spacing: compositionSpacing,
    title: targetPayload('title', [41, 26], 'localBounds', { x: -5, y: -10, width: 60, height: 20 }),
    titleSide: 'top'
  };
  const svg = new FakeElement({
    tagName: 'svg',
    attributes: {
      [COMPOSITION_SCHEMA_ATTRIBUTE]: '1',
      [COMPOSITION_METADATA_ATTRIBUTE]: JSON.stringify(metadata),
      viewBox: '0 0 186 132',
      width: '186px',
      height: '132px'
    }
  });
  const primaryA = new FakeElement({
    id: 'primary-a',
    attributes: {
      [COMPOSITION_ROLE_ATTRIBUTE]: 'primary',
      transform: 'translate(10,33) scale(2) rotate(5) matrix(1 0 0 1 3 4)'
    }
  });
  const primaryB = new FakeElement({
    id: 'primary-b',
    attributes: {
      [COMPOSITION_ROLE_ATTRIBUTE]: 'primary',
      transform: 'translate(6,36) matrix(1 0 0 1 8 9) translate(2e1,-4)'
    }
  });
  const legend = new FakeElement({
    id: 'legend',
    attributes: {
      [COMPOSITION_ROLE_ATTRIBUTE]: 'legend',
      transform: 'translate(149,96) rotate(2)'
    },
    bbox: { x: -2, y: -5, width: 30, height: 20 }
  });
  const title = new FakeElement({
    id: 'plot_title',
    attributes: {
      [COMPOSITION_ROLE_ATTRIBUTE]: 'title',
      transform: 'translate(39,30) scale(.9)'
    },
    bbox: { x: -5, y: -10, width: 60, height: 20 }
  });
  svg.appendChild(primaryA);
  svg.appendChild(primaryB);
  svg.appendChild(legend);
  svg.appendChild(title);
  return { svg, metadata, primaryA, primaryB, legend, title };
};

assert.deepEqual(readLeadingTranslate('translate(1e2,-2.5E+1) scale(2)'), {
  found: true,
  x: 100,
  y: -25,
  prefix: '',
  tail: ' scale(2)'
});
assert.equal(
  replaceLeadingTranslate(
    'translate(1e2,-2.5E+1) scale(2) rotate(5) matrix(1 0 0 1 3 4) translate(.2 4e1)',
    7,
    -8
  ),
  'translate(7,-8) scale(2) rotate(5) matrix(1 0 0 1 3 4) translate(.2 4e1)'
);
assert.equal(prependTranslate('scale(2) translate(3,4)', 0, 0), 'translate(0,0) scale(2) translate(3,4)');

const rightPlan = planComposition({
  primaryBounds: { x: 10, y: 20, width: 100, height: 80 },
  legendBounds: { x: -2, y: -5, width: 30, height: 20 },
  legendSide: 'right',
  spacing: compositionSpacing,
  overlayPolicy
});
assert.deepEqual(rightPlan.placements.primary, {
  automaticTranslation: [6, -4],
  finalBounds: { x: 16, y: 16, width: 100, height: 80 }
});
assert.deepEqual(rightPlan.placements.legend, {
  automaticTranslation: [142, 51],
  finalBounds: { x: 140, y: 46, width: 30, height: 20 }
});
assert.equal(rightPlan.width, 186);
assert.equal(rightPlan.height, 112);

const stackedPlan = planComposition({
  primaryBounds: { x: 10, y: 20, width: 100, height: 80 },
  legendBounds: { x: -2, y: -5, width: 30, height: 20 },
  legendSide: 'top',
  titleBounds: { x: -5, y: -10, width: 60, height: 20 },
  titleSide: 'top',
  spacing: compositionSpacing,
  overlayPolicy
});
assert.deepEqual(stackedPlan.placements.primary.finalBounds, { x: 16, y: 100, width: 100, height: 80 });
assert.deepEqual(stackedPlan.placements.legend.finalBounds, { x: 51, y: 56, width: 30, height: 20 });
assert.deepEqual(stackedPlan.placements.title.finalBounds, { x: 36, y: 16, width: 60, height: 20 });
assert.deepEqual([stackedPlan.width, stackedPlan.height], [132, 196]);

const overlayPlan = planComposition({
  primaryBounds: { x: 0, y: 0, width: 100, height: 100 },
  legendBounds: { x: 0, y: 0, width: 20, height: 20 },
  legendSide: 'upper_left',
  overlayObstacles: [{ x: 0, y: 0, width: 30, height: 30 }],
  spacing: compositionSpacing,
  overlayPolicy
});
assert.equal(overlayPlan.placements.legend.finalBounds.x, 16);
assert.equal(overlayPlan.placements.legend.finalBounds.y, 54);
const narrowQuadrantOverlayPlan = planComposition({
  primaryBounds: { x: 0, y: 0, width: 100, height: 100 },
  legendBounds: { x: 0, y: 0, width: 20, height: 20 },
  legendSide: 'upper_left',
  overlayObstacles: [{ x: 0, y: 0, width: 30, height: 30 }],
  spacing: compositionSpacing,
  overlayPolicy: { ...overlayPolicy, quadrantBoundaryRatio: 0.25 }
});
assert.equal(narrowQuadrantOverlayPlan.width, 160);
assert.deepEqual(narrowQuadrantOverlayPlan.placements.primary.finalBounds, {
  x: 44, y: 16, width: 100, height: 100
});
const candidateOrderArguments = {
  primaryBounds: { x: 0, y: 0, width: 100, height: 100 },
  legendBounds: { x: 0, y: 0, width: 20, height: 20 },
  legendSide: 'upper_left',
  overlayObstacles: [{ x: 0, y: 0, width: 10, height: 20 }],
  spacing: compositionSpacing
};
const pythonOrderedCandidatePlan = planComposition({
  ...candidateOrderArguments,
  overlayPolicy
});
assert.deepEqual(pythonOrderedCandidatePlan.placements.legend.finalBounds, {
  x: 34, y: 16, width: 20, height: 20
});
const xDistanceFirstCandidatePlan = planComposition({
  ...candidateOrderArguments,
  overlayPolicy: {
    ...overlayPolicy,
    candidateScoreOrder: [
      'xAnchorDistance',
      'yAnchorDistance',
      'totalAnchorDistance',
      'nearEdgeX',
      'nearEdgeY'
    ]
  }
});
assert.deepEqual(xDistanceFirstCandidatePlan.placements.legend.finalBounds, {
  x: 16, y: 44, width: 20, height: 20
});

const growthArguments = {
  primaryBounds: { x: 0, y: 0, width: 40, height: 40 },
  legendBounds: { x: 0, y: 0, width: 20, height: 20 },
  legendSide: 'upper_left',
  overlayObstacles: [{ x: 0, y: 0, width: 40, height: 40 }],
  spacing: compositionSpacing
};
const horizontalFirstGrowth = planComposition({
  ...growthArguments,
  overlayPolicy
});
assert.deepEqual(horizontalFirstGrowth.placements.primary.finalBounds, {
  x: 44, y: 16, width: 40, height: 40
});
const verticalFirstGrowth = planComposition({
  ...growthArguments,
  overlayPolicy: {
    ...overlayPolicy,
    canvasGrowthCandidateOrder: ['vertical', 'horizontal']
  }
});
assert.deepEqual(verticalFirstGrowth.placements.primary.finalBounds, {
  x: 16, y: 44, width: 40, height: 40
});
assert.deepEqual(verticalFirstGrowth.placements.legend.finalBounds, {
  x: 16, y: 16, width: 20, height: 20
});
const growthScoreArguments = {
  primaryBounds: { x: 0, y: 0, width: 100, height: 10 },
  legendBounds: { x: 0, y: 0, width: 20, height: 20 },
  legendSide: 'upper_left',
  overlayObstacles: [{ x: 0, y: 0, width: 100, height: 10 }],
  spacing: compositionSpacing
};
const minimumAreaGrowth = planComposition({
  ...growthScoreArguments,
  overlayPolicy
});
assert.deepEqual(minimumAreaGrowth.placements.primary.finalBounds, {
  x: 44, y: 16, width: 100, height: 10
});
const minimumExtentGrowth = planComposition({
  ...growthScoreArguments,
  overlayPolicy: {
    ...overlayPolicy,
    canvasGrowthScoreOrder: ['addedExtent', 'addedArea', 'candidateOrder']
  }
});
assert.deepEqual(minimumExtentGrowth.placements.primary.finalBounds, {
  x: 16, y: 44, width: 100, height: 10
});

{
  const { svg, primaryA, primaryB, legend, title } = schemaOneSvg();
  const before = svg.getAttribute(COMPOSITION_METADATA_ATTRIBUTE);
  const binding = bindCompositionMetadata(svg);
  assert.equal(binding.primary.targets.length, 2);
  assert.equal(svg.getAttribute(COMPOSITION_METADATA_ATTRIBUTE), before);
  assert.deepEqual(compositionUserDeltas(svg), {
    primary: [[4, -3], [0, 0]],
    legend: [7, 5],
    title: [-2, 4]
  });

  const next = applyCompositionEdit(svg, {
    legendSide: 'bottom',
    legendLocalBounds: { x: -2, y: -5, width: 80, height: 25 },
    titleLocalBounds: { x: -5, y: -10, width: 60, height: 20 }
  });
  assert.equal(next.metadata.legendSide, 'bottom');
  assert.deepEqual(compositionUserDeltas(svg), {
    primary: [[4, -3], [0, 0]],
    legend: [7, 5],
    title: [-2, 4]
  });
  assert.match(primaryA.getAttribute('transform'), /^translate\([^)]*\) scale\(2\) rotate\(5\) matrix/);
  assert.match(primaryB.getAttribute('transform'), /matrix\(1 0 0 1 8 9\) translate\(2e1,-4\)$/);
  assert.match(legend.getAttribute('transform'), / rotate\(2\)$/);
  assert.match(title.getAttribute('transform'), / scale\(\.9\)$/);
  assert.equal(svg.getAttribute('viewBox'), `0 0 ${next.metadata.primary.finalBounds.width + 32} ${Number.parseFloat(svg.getAttribute('height'))}`);
  assert.equal(svg.getAttribute('data-horizontal-viewbox'), null);
  assert.equal(svg.getAttribute('data-vertical-viewbox'), null);
}

{
  const { svg } = schemaOneSvg();
  const desired = {
    primary: [[8, 3], [-4, 2]],
    legend: [-6, 7],
    title: [9, -5]
  };
  assert.equal(applyCompositionUserDeltas(svg, desired).changed, true);
  assert.deepEqual(compositionUserDeltas(svg), desired);
  assert.equal(applyCompositionUserDeltas(svg, desired).changed, false);
}

{
  const missing = new FakeElement({ tagName: 'svg' });
  assert.throws(
    () => parseCompositionMetadata(missing),
    (error) => error instanceof CompositionMetadataError && error.code === 'MISSING_COMPOSITION_METADATA'
  );
  missing.setAttribute(COMPOSITION_SCHEMA_ATTRIBUTE, '1');
  missing.setAttribute(COMPOSITION_METADATA_ATTRIBUTE, '{broken');
  assert.throws(() => parseCompositionMetadata(missing), /not valid JSON/);
}

{
  const { svg } = schemaOneSvg();
  const metadata = JSON.parse(svg.getAttribute(COMPOSITION_METADATA_ATTRIBUTE));
  delete metadata.legendReflow;
  svg.setAttribute(COMPOSITION_METADATA_ATTRIBUTE, JSON.stringify(metadata));
  assert.throws(
    () => parseCompositionMetadata(svg),
    /legendReflow is required when a legend is present/
  );
}

{
  const { svg } = schemaOneSvg();
  const metadata = JSON.parse(svg.getAttribute(COMPOSITION_METADATA_ATTRIBUTE));
  metadata.legend = null;
  metadata.legendReflow = legendReflow;
  metadata.legendSide = 'none';
  svg.setAttribute(COMPOSITION_METADATA_ATTRIBUTE, JSON.stringify(metadata));
  assert.throws(
    () => parseCompositionMetadata(svg),
    /legendReflow requires a legend target/
  );
}

{
  const numericFields = [];
  for (const [role, boundsKey] of [
    ['primary', 'finalBounds'],
    ['legend', 'localBounds'],
    ['title', 'localBounds']
  ]) {
    for (const index of [0, 1]) {
      numericFields.push({
        path: `composition.${role}.automaticTranslation[${index}]`,
        set: (metadata, value) => { metadata[role].automaticTranslation[index] = value; }
      });
    }
    for (const field of ['x', 'y', 'width', 'height']) {
      numericFields.push({
        path: `composition.${role}.${boundsKey}.${field}`,
        set: (metadata, value) => { metadata[role][boundsKey][field] = value; }
      });
    }
  }
  for (const field of ['x', 'y', 'width', 'height']) {
    numericFields.push({
      path: `composition.overlayObstacles[0].${field}`,
      set: (metadata, value) => {
        metadata.overlayObstacles = [{ x: 0, y: 0, width: 10, height: 10 }];
        metadata.overlayObstacles[0][field] = value;
      }
    });
  }
  for (const field of [
    'dockGapPx',
    'edgePaddingPx',
    'overlayClearancePx',
    'stackGapPx',
    'titleGapPx'
  ]) {
    numericFields.push({
      path: `composition.spacing.${field}`,
      set: (metadata, value) => { metadata.spacing[field] = value; }
    });
  }
  for (const field of ['colorRectSize', 'lineHeight', 'textXOffset']) {
    numericFields.push({
      path: `composition.legendReflow.${field}`,
      set: (metadata, value) => { metadata.legendReflow[field] = value; }
    });
  }
  numericFields.push({
    path: 'composition.overlayPolicy.quadrantBoundaryRatio',
    set: (metadata, value) => { metadata.overlayPolicy.quadrantBoundaryRatio = value; }
  });

  const invalidWireNumbers = [
    null,
    true,
    false,
    '12',
    '',
    [],
    {},
    undefined,
    Number.NaN,
    Number.POSITIVE_INFINITY,
    Number.NEGATIVE_INFINITY
  ];
  for (const numericField of numericFields) {
    for (const invalid of invalidWireNumbers) {
      const { svg } = schemaOneSvg();
      const metadata = JSON.parse(svg.getAttribute(COMPOSITION_METADATA_ATTRIBUTE));
      numericField.set(metadata, invalid);
      svg.setAttribute(COMPOSITION_METADATA_ATTRIBUTE, JSON.stringify(metadata));
      assert.throws(
        () => parseCompositionMetadata(svg),
        (error) => (
          error instanceof CompositionMetadataError &&
          (
            error.message.includes(numericField.path) ||
            (
              invalid === undefined &&
              numericField.path.startsWith('composition.overlayPolicy.') &&
              error.message.includes('overlayPolicy has an invalid field set')
            )
          )
        ),
        `${numericField.path} accepted ${String(invalid)}`
      );
    }
  }
}

{
  const { svg } = schemaOneSvg();
  const metadata = JSON.parse(svg.getAttribute(COMPOSITION_METADATA_ATTRIBUTE));
  delete metadata.overlayPolicy;
  svg.setAttribute(COMPOSITION_METADATA_ATTRIBUTE, JSON.stringify(metadata));
  assert.throws(() => parseCompositionMetadata(svg), /overlayPolicy must be an object/);

  for (const [field, invalid] of [
    ['candidateScoreOrder', ['nearEdgeX']],
    ['canvasGrowthCandidateOrder', ['horizontal', 'horizontal']],
    ['canvasGrowthScoreOrder', ['addedArea', 'addedExtent', true]],
    ['quadrantBoundaryRatio', '0.5']
  ]) {
    const next = schemaOneSvg();
    const payload = JSON.parse(next.svg.getAttribute(COMPOSITION_METADATA_ATTRIBUTE));
    payload.overlayPolicy[field] = invalid;
    next.svg.setAttribute(COMPOSITION_METADATA_ATTRIBUTE, JSON.stringify(payload));
    assert.throws(
      () => parseCompositionMetadata(next.svg),
      new RegExp(`overlayPolicy\\.${field}`)
    );
  }
}

{
  const { svg } = schemaOneSvg();
  const metadata = JSON.parse(svg.getAttribute(COMPOSITION_METADATA_ATTRIBUTE));
  metadata.title = null;
  metadata.titleSide = 'none';
  svg.setAttribute(COMPOSITION_METADATA_ATTRIBUTE, JSON.stringify(metadata));
  svg.children = svg.children.filter(
    (child) => child.getAttribute(COMPOSITION_ROLE_ATTRIBUTE) !== 'title'
  );

  const title = new FakeElement({
    id: 'plot_title',
    attributes: { transform: 'scale(1.25) rotate(4)' },
    bbox: { x: -8, y: -10, width: 80, height: 20 }
  });
  svg.appendChild(title);
  const added = reconcileCompositionTitle(svg, title, 'top');
  assert.equal(added.metadata.titleSide, 'top');
  assert.equal(title.getAttribute(COMPOSITION_ROLE_ATTRIBUTE), 'title');
  assert.match(title.getAttribute('transform'), /^translate\([^)]*\) scale\(1\.25\) rotate\(4\)$/);
  assert.deepEqual(compositionUserDeltas(svg).title, [0, 0]);

  const removed = reconcileCompositionTitle(svg, null, 'none');
  assert.equal(removed.metadata.title, null);
  assert.equal(removed.metadata.titleSide, 'none');
  assert.equal(title.getAttribute(COMPOSITION_ROLE_ATTRIBUTE), null);
}

{
  const legacy = new FakeElement({ tagName: 'svg', attributes: { viewBox: '0 0 200 100' } });
  const primary = new FakeElement({
    id: 'record',
    attributes: { transform: 'translate(10,20) scale(2)' },
    bbox: { x: 0, y: 0, width: 40, height: 20 }
  });
  const legend = new FakeElement({
    id: 'legend',
    attributes: { transform: 'translate(100,5) rotate(3)' },
    bbox: { x: 0, y: 0, width: 30, height: 15 }
  });
  legacy.appendChild(primary);
  legacy.appendChild(legend);
  const binding = normalizeLegacyComposition(legacy, { legendSide: 'right', titleSide: 'none' });
  assert.equal(binding.metadata.legacyNormalized, true);
  assert.equal(primary.getAttribute('transform'), 'translate(0,0) translate(10,20) scale(2)');
  assert.equal(legend.getAttribute('transform'), 'translate(0,0) translate(100,5) rotate(3)');
  assert.deepEqual(compositionUserDeltas(legacy), {
    primary: [[0, 0]],
    legend: [0, 0],
    title: null
  });
}

{
  const legacy = new FakeElement({ tagName: 'svg', attributes: { viewBox: '0 0 200 100' } });
  const primary = new FakeElement({
    id: 'record',
    attributes: { transform: 'translate(10,20) scale(2)' },
    bbox: { x: 0, y: 0, width: 40, height: 20 }
  });
  const legend = new FakeElement({
    id: 'legend',
    attributes: { transform: 'translate(100,5) rotate(3)' },
    bbox: { x: 0, y: 0, width: 30, height: 15 }
  });
  const lengthBar = new FakeElement({
    id: 'length_bar',
    attributes: { transform: 'translate(30,40) scale(1)' },
    bbox: { x: 0, y: 0, width: 25, height: 5 }
  });
  legacy.appendChild(primary);
  legacy.appendChild(lengthBar);
  legacy.appendChild(legend);
  normalizeLegacyComposition(legacy, {
    legendSide: 'right',
    userDeltas: { primary: [3, -4], lengthBar: [9, 1], legend: [5, 2] }
  });
  assert.equal(primary.getAttribute('transform'), 'translate(3,-4) translate(7,24) scale(2)');
  assert.equal(lengthBar.getAttribute('transform'), 'translate(9,1) translate(21,39) scale(1)');
  assert.equal(legend.getAttribute('transform'), 'translate(5,2) translate(95,3) rotate(3)');
  assert.deepEqual(compositionUserDeltas(legacy), {
    primary: [[3, -4], [9, 1]],
    legend: [5, 2],
    title: null
  });
  resetCompositionUserDeltas(legacy);
  assert.equal(primary.getAttribute('transform'), 'translate(0,0) translate(7,24) scale(2)');
  assert.equal(lengthBar.getAttribute('transform'), 'translate(0,0) translate(21,39) scale(1)');
  assert.equal(legend.getAttribute('transform'), 'translate(0,0) translate(95,3) rotate(3)');
}

console.log('composition layout tests passed');
