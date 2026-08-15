import assert from 'node:assert/strict';

import { createLegendDragActions } from '../../gbdraw/web/js/app/legend/drag-actions.js';
import { createLegendLayoutActions } from '../../gbdraw/web/js/app/legend/layout-actions.js';
import {
  COMPOSITION_METADATA_ATTRIBUTE,
  COMPOSITION_ROLE_ATTRIBUTE,
  COMPOSITION_SCHEMA_ATTRIBUTE
} from '../../gbdraw/web/js/app/legend-layout/composition-actions.js';
import { parseTransformXY } from '../../gbdraw/web/js/app/legend/utils.js';

const targetPayload = (role, automaticTranslation, boundsKey, bounds) => ({
  automaticTranslation,
  [boundsKey]: bounds,
  role,
  selector: `[data-gbdraw-composition-role="${role}"]`
});

const compositionMetadata = (primaryWidth) => ({
  legend: targetPayload('legend', [0, 0], 'localBounds', {
    x: 0,
    y: 0,
    width: primaryWidth,
    height: 14
  }),
  legendReflow: {
    colorRectSize: 14,
    lineHeight: 24,
    textXOffset: 22
  },
  legendSide: 'top',
  overlayObstacles: [],
  overlayPolicy: {
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
  },
  primary: targetPayload('primary', [0, 0], 'finalBounds', {
    x: 0,
    y: 0,
    width: primaryWidth,
    height: 100
  }),
  spacing: {
    dockGapPx: 24,
    edgePaddingPx: 16,
    overlayClearancePx: 8,
    stackGapPx: 20,
    titleGapPx: 20
  },
  title: null,
  titleSide: 'none'
});

class AttributeNode {
  constructor({ id = '', attributes = {}, bbox = { x: 0, y: 0, width: 0, height: 0 } } = {}) {
    this.id = id;
    this.attributes = new Map(Object.entries(attributes));
    this.bbox = { ...bbox };
    this.children = [];
    this.style = {};
  }

  getAttribute(name) {
    return this.attributes.has(name) ? this.attributes.get(name) : null;
  }

  setAttribute(name, value) {
    this.attributes.set(name, String(value));
  }

  removeAttribute(name) {
    this.attributes.delete(name);
  }

  getBBox() {
    return { ...this.bbox };
  }

  closest() {
    return null;
  }
}

const layoutFixture = (primaryWidth) => {
  const texts = [
    new AttributeNode({
      attributes: { transform: 'translate(22,7)' },
      bbox: { x: 0, y: 0, width: 10, height: 10 }
    }),
    new AttributeNode({
      attributes: { transform: 'translate(90,7)' },
      bbox: { x: 0, y: 0, width: 10, height: 10 }
    })
  ];
  const paths = [
    new AttributeNode({ attributes: { fill: '#112233', transform: 'translate(0,7)' } }),
    new AttributeNode({ attributes: { fill: '#445566', transform: 'translate(68,7)' } })
  ];
  const featureLegend = new AttributeNode({ id: 'feature_legend' });
  featureLegend.querySelectorAll = (selector) => {
    if (selector === 'text') return texts;
    if (selector === 'path') return paths;
    return [];
  };

  const legend = new AttributeNode({
    id: 'legend',
    bbox: { x: 0, y: 0, width: primaryWidth, height: 14 }
  });
  legend.children = [featureLegend];
  legend.querySelector = (selector) => (selector === '#feature_legend' ? featureLegend : null);

  const svg = new AttributeNode({
    attributes: {
      [COMPOSITION_SCHEMA_ATTRIBUTE]: '1',
      [COMPOSITION_METADATA_ATTRIBUTE]: JSON.stringify(compositionMetadata(primaryWidth))
    }
  });
  svg.getElementById = (id) => (id === 'legend' ? legend : null);
  return { svg, texts };
};

{
  const actions = createLegendLayoutActions();
  const exactFit = layoutFixture(136);
  actions.updatePairwiseLegendPositions(exactFit.svg);
  assert.deepEqual(
    exactFit.texts.map((text) => parseTransformXY(text.getAttribute('transform')).y),
    [7, 7],
    'two 68 px entries must remain on one row at the exact 136 px primary-width boundary'
  );

  const oneMillipixelNarrow = layoutFixture(135.999);
  actions.updatePairwiseLegendPositions(oneMillipixelNarrow.svg);
  assert.deepEqual(
    oneMillipixelNarrow.texts.map((text) => parseTransformXY(text.getAttribute('transform')).y),
    [7, 31],
    'the second entry must wrap only when its measured row exceeds the full primary width'
  );
}

const dragFixture = () => {
  const metadata = compositionMetadata(200);
  metadata.legend.automaticTranslation = [10, 20];
  metadata.legendSide = 'right';

  const primary = new AttributeNode({
    attributes: {
      [COMPOSITION_ROLE_ATTRIBUTE]: 'primary',
      transform: 'translate(0,0)'
    }
  });
  const legend = new AttributeNode({
    id: 'legend',
    attributes: {
      [COMPOSITION_ROLE_ATTRIBUTE]: 'legend',
      transform: 'translate(10,20) rotate(2)'
    }
  });
  const targets = { primary: [primary], legend: [legend], title: [] };
  const svg = new AttributeNode({
    attributes: {
      [COMPOSITION_SCHEMA_ATTRIBUTE]: '1',
      [COMPOSITION_METADATA_ATTRIBUTE]: JSON.stringify(metadata)
    }
  });
  svg.selectorCalls = 0;
  svg.querySelectorAll = (selector) => {
    svg.selectorCalls += 1;
    const match = selector.match(/^\[data-gbdraw-composition-role="([^"]+)"\]$/);
    return match ? targets[match[1]] : [];
  };

  const ref = (value) => ({ value });
  const state = {
    results: ref([]),
    selectedResultIndex: ref(-1),
    svgContainer: ref({ querySelector: (selector) => (selector === 'svg' ? svg : null) }),
    legendDragging: ref(false),
    legendDragStart: { x: 0, y: 0 },
    legendOriginalTransform: ref({ x: 0, y: 0 }),
    legendInitialTransform: ref({ x: 10, y: 20 }),
    legendCurrentOffset: { x: 0, y: 0 },
    layoutRepositionMode: ref(true),
    zoom: ref(1),
    skipCaptureBaseConfig: ref(false)
  };
  return { legend, state, svg };
};

{
  const originalRequestAnimationFrame = globalThis.requestAnimationFrame;
  const originalCancelAnimationFrame = globalThis.cancelAnimationFrame;
  const frames = new Map();
  let nextFrameId = 1;
  globalThis.requestAnimationFrame = (callback) => {
    const frameId = nextFrameId;
    nextFrameId += 1;
    frames.set(frameId, callback);
    return frameId;
  };
  globalThis.cancelAnimationFrame = (frameId) => frames.delete(frameId);
  const flushFrames = () => {
    const queued = [...frames.values()];
    frames.clear();
    queued.forEach((callback) => callback());
  };

  try {
    const { legend, state, svg } = dragFixture();
    const actions = createLegendDragActions({
      state,
      extractLegendEntries: () => {},
      previewRuntime: {
        commitDomEdit: ({ reason = 'test-edit', mutate }) => {
          const changed = Boolean(mutate({ svg, resultIndex: 0 }));
          return { changed, flushed: changed, resultIndex: 0, reason };
        }
      }
    });
    actions.startLegendDrag({
      clientX: 100,
      clientY: 200,
      preventDefault: () => {},
      shiftKey: false,
      stopPropagation: () => {}
    });
    const selectorCallsAtStart = svg.selectorCalls;
    assert.ok(selectorCallsAtStart > 0, 'drag start must bind the composition roles');

    actions.onLegendDrag({ clientX: 105, clientY: 207 });
    actions.onLegendDrag({ clientX: 108, clientY: 211 });
    flushFrames();
    actions.onLegendDrag({ clientX: 110, clientY: 214 });
    flushFrames();
    await actions.endLegendDrag({ clientX: 112, clientY: 216 });

    assert.equal(
      svg.selectorCalls,
      selectorCallsAtStart,
      'pointer frames and drag completion must reuse the composition binding from drag start'
    );
    assert.equal(legend.getAttribute('transform'), 'translate(22,36) rotate(2)');
    assert.deepEqual(state.legendCurrentOffset, { x: 12, y: 16 });
  } finally {
    globalThis.requestAnimationFrame = originalRequestAnimationFrame;
    globalThis.cancelAnimationFrame = originalCancelAnimationFrame;
  }
}

console.log('legend layout action tests passed');
