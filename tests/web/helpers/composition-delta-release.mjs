// Drives the Legend layout owner's composition offsets over stand-in SVGs
// whose offsets the mocked composition readers keep on the element, and prints
// what the owner records and restores after each release as JSON. The module
// mocks need `--experimental-test-module-mocks`, so a test runs this file in a
// child process (OV-110).
import { mock } from 'node:test';

const at = (path) => new URL(path, import.meta.url).href;
const compositionUrl = at('../../../gbdraw/web/js/app/legend-layout/composition-actions.js');
const real = await import(`${compositionUrl}?real`);
mock.module(compositionUrl, {
  namedExports: {
    ...real,
    compositionUserDeltas: (svg) => JSON.parse(JSON.stringify(svg.deltas)),
    applyCompositionUserDeltas: (svg, deltas) => {
      const changed = JSON.stringify(svg.deltas) !== JSON.stringify(deltas);
      svg.deltas = JSON.parse(JSON.stringify(deltas));
      return { changed, binding: null };
    },
    // The stand-ins carry no canvas, so fitting one changes nothing.
    applyCanvasToSvg: () => false
  }
});
const { createLegendLayout } = await import(at('../../../gbdraw/web/js/app/legend-layout.js'));

const unmoved = { primary: [], legend: [0, 0], title: null };
const moved = { primary: [], legend: [5, 5], title: null };
const svgs = { R1: { deltas: unmoved }, R2: { deltas: unmoved } };
const state = {
  results: { value: [{ id: 'R1' }] },
  svgContainer: { value: { querySelector: () => svgs.R1 } },
  sessionOperationAvailability: () => null,
  activeDrawing: () => ({ canvasPadding: null })
};
let mounted = true;
const restored = [];
const layout = createLegendLayout({
  state,
  layOutLegend: () => null,
  previewRuntime: {
    commitActiveResultEdit: () => true,
    getResultIdentity: (result) => result.id,
    getActiveRuntime: () => (mounted ? { svg: svgs.R1, resultIdentity: 'R1' } : null),
    commitResultEdit: (index, edit) => {
      const id = state.results.value[index].id;
      restored.push(id);
      return edit(svgs[id], { mounted: false });
    }
  }
});
const out = {};
layout.captureCompositionIntent(['R1']);
svgs.R1.deltas = moved;
out.moved = layout.captureCompositionIntent(['R1']);
mounted = false;
// (b) a History step names R1: its offsets stay, and its Undo returns R1 to
// the offsets first seen.
layout.releaseCompositionDeltas(new Set(['R1']));
out.kept = layout.captureCompositionIntent(['R1']);
layout.reconcileCompositionUserDeltas({}, new Set(['R1']));
out.restoredTo = svgs.R1.deltas;
out.restored = [...restored];
// (a) R1 dropped: nothing names it, so its offsets are released.
svgs.R1.deltas = moved;
layout.reconcileCompositionUserDeltas({ R1: moved }, new Set(['R1']));
layout.releaseCompositionDeltas(new Set());
out.released = layout.captureCompositionIntent(['R1']);
restored.length = 0;
layout.reconcileCompositionUserDeltas({}, new Set(['R1']));
out.restoredAfterRelease = [...restored];
console.log(JSON.stringify(out));
