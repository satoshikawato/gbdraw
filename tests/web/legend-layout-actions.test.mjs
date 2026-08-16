import assert from 'node:assert/strict';

import { createLegendDragActions } from '../../gbdraw/web/js/app/legend/drag-actions.js';
import { createDiagramDragActions } from '../../gbdraw/web/js/app/legend-layout/diagram-drag.js';
import { createLegendLayout } from '../../gbdraw/web/js/app/legend-layout.js';
import { createDomMutationJournal } from '../../gbdraw/web/js/app/dom-mutation-journal.js';
import { createLegendLayoutActions } from '../../gbdraw/web/js/app/legend/layout-actions.js';
import { createPreviewRuntime } from '../../gbdraw/web/js/app/preview-runtime.js';
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
  constructor({
    id = '',
    tagName = 'g',
    attributes = {},
    bbox = { x: 0, y: 0, width: 0, height: 0 }
  } = {}) {
    this.id = id;
    this.tagName = tagName;
    this.attributes = new Map(Object.entries(attributes));
    this.bbox = { ...bbox };
    this.children = [];
    this.parentElement = null;
    this.listeners = new Map();
    this.style = {
      removeProperty: (name) => { delete this.style[name]; }
    };
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

  hasAttribute(name) {
    return this.attributes.has(name);
  }

  appendChild(child) {
    child.parentElement = this;
    this.children.push(child);
    return child;
  }

  addEventListener(name, callback) {
    this.listeners.set(name, callback);
  }

  removeEventListener(name, callback) {
    if (this.listeners.get(name) === callback) this.listeners.delete(name);
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

  const rollbackFixture = layoutFixture(135.999);
  const transformsBefore = rollbackFixture.texts.map((text) => text.getAttribute('transform'));
  const mutation = createDomMutationJournal();
  actions.updatePairwiseLegendPositions(rollbackFixture.svg, mutation);
  assert.notDeepEqual(
    rollbackFixture.texts.map((text) => text.getAttribute('transform')),
    transformsBefore
  );
  assert.deepEqual(mutation.rollback(), []);
  assert.deepEqual(
    rollbackFixture.texts.map((text) => text.getAttribute('transform')),
    transformsBefore,
    'Legend reflow attributes participate in the caller transaction journal'
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
  const originalDocument = globalThis.document;
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
        beginDomEditLease: () => ({
          current: true,
          mutate: (mutate) => Boolean(mutate({
            svg,
            resultIndex: 0,
            mutation: {
              captureProperty: () => true,
              captureState: () => true,
              setAttribute: (element, name, value) => {
                if (element.getAttribute(name) === String(value)) return false;
                element.setAttribute(name, value);
                return true;
              }
            }
          })),
          commit: () => true,
          cancel: () => false
        }),
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

    const failedFixture = dragFixture();
    const originalResult = { name: 'legend-drag.svg', content: '<svg data-owner="before"></svg>' };
    failedFixture.state.results = { value: [originalResult] };
    failedFixture.state.selectedResultIndex = { value: 0 };
    failedFixture.state.resultGenerationKey = { value: 'legend-drag-generation' };
    failedFixture.state.skipCaptureBaseConfig = { value: false };
    const failedRuntime = createPreviewRuntime({
      state: failedFixture.state,
      serializeSvg: () => {
        throw new Error('serialization failed');
      }
    });
    failedRuntime.mountResultSvg(0, failedFixture.svg);
    let historyCommitCount = 0;
    let historyCancelCount = 0;
    const failedActions = createLegendDragActions({
      state: failedFixture.state,
      extractLegendEntries: () => {},
      history: {
        begin: async () => ({ closed: false }),
        commit: async () => { historyCommitCount += 1; },
        cancel: (tx) => {
          tx.closed = true;
          historyCancelCount += 1;
        }
      },
      previewRuntime: failedRuntime
    });
    const transformBefore = failedFixture.legend.getAttribute('transform');
    const offsetBefore = { ...failedFixture.state.legendCurrentOffset };
    failedActions.startLegendDrag({
      clientX: 100,
      clientY: 200,
      preventDefault: () => {},
      shiftKey: false,
      stopPropagation: () => {}
    });
    failedActions.onLegendDrag({ clientX: 112, clientY: 216 });
    flushFrames();
    await assert.rejects(
      failedActions.endLegendDrag({ clientX: 112, clientY: 216 }),
      /serialization failed/
    );
    assert.equal(failedFixture.legend.getAttribute('transform'), transformBefore);
    assert.deepEqual(failedFixture.state.legendCurrentOffset, offsetBefore);
    assert.equal(failedFixture.state.results.value[0], originalResult);
    assert.equal(failedRuntime.getActiveRuntime().dirty, false);
    assert.equal(historyCommitCount, 0);
    assert.equal(historyCancelCount, 1);

    const diagramFailureFixture = (mode) => {
      const metadata = compositionMetadata(200);
      metadata.primary = targetPayload('primary', [10, 20], 'finalBounds', {
        x: 10, y: 20, width: 200, height: 100
      });
      metadata.legend.automaticTranslation = [234, 20];
      metadata.legendSide = 'right';
      metadata.title = targetPayload('title', [80, 144], 'localBounds', {
        x: 0, y: 0, width: 60, height: 20
      });
      metadata.titleSide = 'bottom';
      const primary = new AttributeNode({
        id: 'plot',
        attributes: {
          [COMPOSITION_ROLE_ATTRIBUTE]: 'primary',
          transform: 'translate(10,20) scale(1)'
        },
        bbox: { x: 0, y: 0, width: 200, height: 100 }
      });
      const legend = new AttributeNode({
        id: 'legend',
        attributes: {
          [COMPOSITION_ROLE_ATTRIBUTE]: 'legend',
          transform: 'translate(234,20)'
        },
        bbox: { x: 0, y: 0, width: 50, height: 30 }
      });
      const lengthBar = new AttributeNode({
        id: 'length_bar',
        attributes: { transform: 'translate(5,6)' }
      });
      const title = new AttributeNode({
        id: 'plot_title',
        attributes: {
          [COMPOSITION_ROLE_ATTRIBUTE]: 'title',
          transform: 'translate(80,144)'
        },
        bbox: { x: 0, y: 0, width: 90, height: 24 }
      });
      const svg = new AttributeNode({
        tagName: 'svg',
        attributes: {
          [COMPOSITION_SCHEMA_ATTRIBUTE]: '1',
          [COMPOSITION_METADATA_ATTRIBUTE]: JSON.stringify(metadata),
          viewBox: '0 0 332 194',
          width: '332px',
          height: '194px'
        }
      });
      svg.dataset = {};
      [primary, legend, lengthBar, title].forEach((element) => svg.appendChild(element));
      const targets = { primary: [primary], legend: [legend], title: [title] };
      svg.querySelectorAll = (selector) => {
        const match = selector.match(/^\[data-gbdraw-composition-role="([^"]+)"\]$/);
        return match ? targets[match[1]] : [];
      };
      svg.getElementById = (id) => ({
        plot: primary,
        legend,
        length_bar: lengthBar,
        plot_title: title
      })[id] || null;
      const dragTarget = mode === 'length_bar' ? lengthBar : mode === 'plot_title' ? title : primary;
      dragTarget.closest = (selector) => {
        if (selector === '#length_bar') return mode === 'length_bar' ? lengthBar : null;
        if (selector === '#plot_title') return mode === 'plot_title' ? title : null;
        if (selector === 'g[id]') return mode === 'group' ? primary : dragTarget;
        return null;
      };
      const ref = (value) => ({ value });
      const originalResult = { name: `${mode}.svg`, content: '<svg data-owner="before"></svg>' };
      const state = {
        results: ref([originalResult]),
        selectedResultIndex: ref(0),
        resultGenerationKey: ref(`diagram-${mode}`),
        skipCaptureBaseConfig: ref(false),
        svgContainer: ref({ querySelector: (selector) => selector === 'svg' ? svg : null }),
        diagramElements: ref([primary]),
        diagramElementIds: ref(['plot']),
        diagramElementOriginalTransforms: ref(new Map([[primary, { x: 10, y: 20 }]])),
        diagramOffset: { x: 0, y: 0 },
        diagramDragging: ref(false),
        diagramDragStart: { x: 0, y: 0 },
        lengthBarElement: ref(lengthBar),
        lengthBarOriginalTransform: ref({ x: 5, y: 6 }),
        lengthBarUserOffset: { x: 3, y: 4 },
        plotTitleElement: ref(title),
        plotTitleDragging: ref(false),
        plotTitleDragStart: { x: 0, y: 0 },
        plotTitleAutoTransform: ref({ x: 80, y: 144 }),
        plotTitleUserOffset: { x: 2, y: 3 },
        svgContent: ref('<svg/>'),
        canvasPadding: { top: 3, right: 4, bottom: 5, left: 2 },
        originalSvgStroke: ref({ color: '#222', width: 1 }),
        generatedLegendPosition: ref('right'),
        legendInitialTransform: ref({ x: 234, y: 20 }),
        legendCurrentOffset: { x: 2, y: 5 },
        skipPositionReapply: ref(false),
        layoutRepositionMode: ref(true),
        zoom: ref(1)
      };
      const runtime = createPreviewRuntime({
        state,
        serializeSvg: () => { throw new Error('serialization failed'); }
      });
      runtime.mountResultSvg(0, svg);
      const featureIndex = new Map([['feature', [primary]]]);
      const legendIndex = new Map([['legend', [legend]]]);
      runtime.getActiveRuntime().indexes.features = featureIndex;
      runtime.getActiveRuntime().indexes.legend = legendIndex;
      return {
        dragTarget,
        featureIndex,
        legend,
        legendIndex,
        lengthBar,
        originalResult,
        primary,
        runtime,
        state,
        svg,
        title
      };
    };

    for (const mode of ['group', 'length_bar', 'plot_title']) {
      const failed = diagramFailureFixture(mode);
      let historyCommits = 0;
      let historyCancels = 0;
      const listeners = new Map();
      globalThis.document = {
        addEventListener: (name, callback) => listeners.set(name, callback),
        removeEventListener: (name, callback) => {
          if (listeners.get(name) === callback) listeners.delete(name);
        }
      };
      const actions = createDiagramDragActions({
        state: failed.state,
        history: {
          begin: async () => ({ closed: false }),
          commit: async () => { historyCommits += 1; },
          cancel: (tx) => {
            tx.closed = true;
            historyCancels += 1;
          }
        },
        previewRuntime: failed.runtime
      });
      const before = {
        primaryTransform: failed.primary.getAttribute('transform'),
        lengthTransform: failed.lengthBar.getAttribute('transform'),
        titleTransform: failed.title.getAttribute('transform'),
        diagramOffset: { ...failed.state.diagramOffset },
        lengthOffset: { ...failed.state.lengthBarUserOffset },
        titleOffset: { ...failed.state.plotTitleUserOffset }
      };
      actions.startDiagramDrag({
        target: failed.dragTarget,
        clientX: 100,
        clientY: 200,
        preventDefault: () => {},
        shiftKey: false
      });
      actions.onDiagramDrag({ clientX: 112, clientY: 216 });
      flushFrames();
      await assert.rejects(
        actions.endDiagramDrag({ clientX: 112, clientY: 216 }),
        /serialization failed/,
        mode
      );
      assert.equal(failed.primary.getAttribute('transform'), before.primaryTransform, mode);
      assert.equal(failed.lengthBar.getAttribute('transform'), before.lengthTransform, mode);
      assert.equal(failed.title.getAttribute('transform'), before.titleTransform, mode);
      assert.deepEqual(failed.state.diagramOffset, before.diagramOffset, mode);
      assert.deepEqual(failed.state.lengthBarUserOffset, before.lengthOffset, mode);
      assert.deepEqual(failed.state.plotTitleUserOffset, before.titleOffset, mode);
      assert.equal(failed.state.diagramDragging.value, false, mode);
      assert.equal(failed.state.plotTitleDragging.value, false, mode);
      assert.equal(listeners.size, 0, mode);
      assert.equal(failed.state.results.value[0], failed.originalResult, mode);
      assert.equal(failed.runtime.getActiveRuntime().dirty, false, mode);
      assert.equal(historyCommits, 0, mode);
      assert.equal(historyCancels, 1, mode);
    }

    for (const actionName of [
      'resetAllPositions',
      'refreshCompositionGeometry',
      'resetCanvasPadding'
    ]) {
      const failed = diagramFailureFixture('group');
      const layout = createLegendLayout({
        state: failed.state,
        legendActions: {
          reflowDualLegendLayout: () => {},
          reflowSingleLegendLayout: () => {}
        },
        previewRuntime: failed.runtime
      });
      const before = {
        svgAttributes: [...failed.svg.attributes],
        primaryTransform: failed.primary.getAttribute('transform'),
        legendTransform: failed.legend.getAttribute('transform'),
        titleTransform: failed.title.getAttribute('transform'),
        lengthTransform: failed.lengthBar.getAttribute('transform'),
        canvasPadding: { ...failed.state.canvasPadding },
        diagramOffset: { ...failed.state.diagramOffset },
        lengthOffset: { ...failed.state.lengthBarUserOffset },
        legendOffset: { ...failed.state.legendCurrentOffset },
        titleOffset: { ...failed.state.plotTitleUserOffset },
        diagramElementsOwner: failed.state.diagramElements.value,
        originalTransformsOwner: failed.state.diagramElementOriginalTransforms.value,
        legendInitialOwner: failed.state.legendInitialTransform.value,
        titleAutoOwner: failed.state.plotTitleAutoTransform.value,
        generatedLegendPosition: failed.state.generatedLegendPosition.value
      };

      assert.throws(() => layout[actionName](), /serialization failed/, actionName);

      assert.deepEqual([...failed.svg.attributes], before.svgAttributes, actionName);
      assert.equal(failed.primary.getAttribute('transform'), before.primaryTransform, actionName);
      assert.equal(failed.legend.getAttribute('transform'), before.legendTransform, actionName);
      assert.equal(failed.title.getAttribute('transform'), before.titleTransform, actionName);
      assert.equal(failed.lengthBar.getAttribute('transform'), before.lengthTransform, actionName);
      assert.deepEqual(failed.state.canvasPadding, before.canvasPadding, actionName);
      assert.deepEqual(failed.state.diagramOffset, before.diagramOffset, actionName);
      assert.deepEqual(failed.state.lengthBarUserOffset, before.lengthOffset, actionName);
      assert.deepEqual(failed.state.legendCurrentOffset, before.legendOffset, actionName);
      assert.deepEqual(failed.state.plotTitleUserOffset, before.titleOffset, actionName);
      assert.equal(failed.state.diagramElements.value, before.diagramElementsOwner, actionName);
      assert.equal(
        failed.state.diagramElementOriginalTransforms.value,
        before.originalTransformsOwner,
        actionName
      );
      assert.equal(failed.state.legendInitialTransform.value, before.legendInitialOwner, actionName);
      assert.equal(failed.state.plotTitleAutoTransform.value, before.titleAutoOwner, actionName);
      assert.equal(
        failed.state.generatedLegendPosition.value,
        before.generatedLegendPosition,
        actionName
      );
      assert.equal(failed.state.results.value[0], failed.originalResult, actionName);
      assert.equal(failed.runtime.getActiveRuntime().dirty, false, actionName);
      assert.equal(failed.runtime.getActiveRuntime().indexes.features, failed.featureIndex, actionName);
      assert.equal(failed.runtime.getActiveRuntime().indexes.legend, failed.legendIndex, actionName);
    }
  } finally {
    globalThis.requestAnimationFrame = originalRequestAnimationFrame;
    globalThis.cancelAnimationFrame = originalCancelAnimationFrame;
    globalThis.document = originalDocument;
  }
}

console.log('legend layout action tests passed');
