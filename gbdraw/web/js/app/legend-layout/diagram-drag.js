import { parseTransform } from './transform-utils.js';

export const createDiagramDragActions = ({ state }) => {
  const {
    svgContainer,
    diagramElements,
    diagramElementIds,
    diagramElementOriginalTransforms,
    diagramOffset,
    diagramDragging,
    diagramDragStart,
    zoom,
    generatedLegendPosition
  } = state;

  const applyDiagramShift = (deltaX, deltaY) => {
    if (diagramElements.value.length === 0) return;
    diagramElements.value.forEach((el) => {
      const current = parseTransform(el.getAttribute('transform'));
      const newX = current.x + deltaX;
      const newY = current.y + deltaY;
      el.setAttribute('transform', `translate(${newX}, ${newY})`);
      diagramElementOriginalTransforms.value.set(el, { x: newX, y: newY });
    });
    diagramOffset.x = 0;
    diagramOffset.y = 0;
  };

  const startDiagramDrag = (e) => {
    if (
      e.target.closest('#legend') ||
      e.target.closest('path[data-feature-id]') ||
      e.target.closest('rect[data-feature-id]')
    ) {
      return;
    }
    if (!svgContainer.value) return;
    const svg = svgContainer.value.querySelector('svg');
    if (!svg) return;

    const clickedGroup = e.target.closest('g[id]');
    if (!clickedGroup) return;

    const clickedId = clickedGroup.id;
    if (
      clickedId === 'legend' ||
      clickedId === 'feature_legend' ||
      clickedId === 'pairwise_legend' ||
      clickedId === 'horizontal_legend' ||
      clickedId === 'vertical_legend'
    ) {
      return;
    }

    e.preventDefault();
    diagramDragging.value = true;
    diagramDragStart.x = e.clientX;
    diagramDragStart.y = e.clientY;

    diagramElements.value.forEach((el) => {
      el.style.opacity = '0.8';
    });

    document.addEventListener('mousemove', onDiagramDrag);
    document.addEventListener('mouseup', endDiagramDrag);
  };

  const onDiagramDrag = (e) => {
    if (!diagramDragging.value) return;

    const deltaX = (e.clientX - diagramDragStart.x) / zoom.value;
    const deltaY = (e.clientY - diagramDragStart.y) / zoom.value;

    diagramElements.value.forEach((el) => {
      const original = diagramElementOriginalTransforms.value.get(el) || { x: 0, y: 0 };
      const newX = original.x + diagramOffset.x + deltaX;
      const newY = original.y + diagramOffset.y + deltaY;
      el.setAttribute('transform', `translate(${newX}, ${newY})`);
    });
  };

  const endDiagramDrag = (e) => {
    if (!diagramDragging.value) return;

    const deltaX = (e.clientX - diagramDragStart.x) / zoom.value;
    const deltaY = (e.clientY - diagramDragStart.y) / zoom.value;
    diagramOffset.x += deltaX;
    diagramOffset.y += deltaY;

    diagramDragging.value = false;
    document.removeEventListener('mousemove', onDiagramDrag);
    document.removeEventListener('mouseup', endDiagramDrag);

    diagramElements.value.forEach((el) => {
      el.style.opacity = '1';
    });
  };

  const resetDiagramPosition = () => {
    diagramOffset.x = 0;
    diagramOffset.y = 0;

    console.log('[DEBUG] resetDiagramPosition called');
    console.log('[DEBUG] diagramElements count:', diagramElements.value.length);
    console.log('[DEBUG] originalTransforms size:', diagramElementOriginalTransforms.value.size);
    console.log('[DEBUG] All originalTransforms entries:');
    diagramElementOriginalTransforms.value.forEach((transform, mapEl) => {
      console.log(`[DEBUG]   Map entry: ${mapEl.id} -> (${transform.x}, ${transform.y})`);
    });

    diagramElements.value.forEach((el, idx) => {
      const original = diagramElementOriginalTransforms.value.get(el);
      const currentTransform = el.getAttribute('transform');
      let foundInMap = false;
      diagramElementOriginalTransforms.value.forEach((_, mapEl) => {
        if (mapEl === el) foundInMap = true;
      });
      console.log(
        `[DEBUG] Reset element ${idx} (${el.id}): current="${currentTransform}", foundInMap=${foundInMap}, original=`,
        original
      );
      if (original && (original.x !== 0 || original.y !== 0)) {
        el.setAttribute('transform', `translate(${original.x}, ${original.y})`);
      } else {
        el.removeAttribute('transform');
      }
    });
  };

  let setupDiagramDragCallCount = 0;

  const setupDiagramDrag = (preserveOffset = false) => {
    setupDiagramDragCallCount++;
    if (!svgContainer.value) return;
    const svg = svgContainer.value.querySelector('svg');
    if (!svg) return;

    if (!preserveOffset) {
      diagramOffset.x = 0;
      diagramOffset.y = 0;
    }

    const knownIds = ['tick', 'labels', 'Axis', 'gc_content', 'skew', 'gc_skew', 'length_bar'];
    const foundElements = [];
    const foundIds = [];

    knownIds.forEach((id) => {
      const el = svg.getElementById(id);
      if (el) {
        foundElements.push(el);
        foundIds.push(id);
      }
    });

    const allGroups = svg.querySelectorAll('g[id]');
    allGroups.forEach((group) => {
      const id = group.id;
      if (!id) return;

      if (id === 'legend' || id === 'feature_legend' || id === 'pairwise_legend') return;

      if (foundElements.includes(group)) return;

      const isAccession = id.match(/^[A-Z]{2}_?\d+/) || id.match(/^[A-Z]+\d+\.\d+$/);
      const isKnown = knownIds.includes(id);
      const isDynamic =
        id.startsWith('record_') ||
        id.startsWith('definition_') ||
        id.startsWith('seq_') ||
        id.startsWith('track_') ||
        id.startsWith('match_') ||
        id.startsWith('comparison');

      if (isAccession || isKnown || isDynamic) {
        foundElements.push(group);
        if (!foundIds.includes(id)) {
          foundIds.push(id);
        }
      }
    });

    // Ensure top-level diagram groups (e.g., contig_1) are included even if IDs are lowercase.
    Array.from(svg.children).forEach((el) => {
      if (!el || el.tagName.toLowerCase() !== 'g') return;
      const id = el.getAttribute('id');
      if (!id) return;
      if (id === 'legend') return;
      if (foundElements.includes(el)) return;
      foundElements.push(el);
      if (!foundIds.includes(id)) {
        foundIds.push(id);
      }
    });

    diagramElements.value = foundElements;
    diagramElementIds.value = foundIds;

    const originalTransforms = new Map();
    console.log(`[DEBUG] ========== setupDiagramDrag CALL #${setupDiagramDragCallCount} ==========`);
    console.log(
      `[DEBUG] setupDiagramDrag: preserveOffset=${preserveOffset}, offset=(${diagramOffset.x}, ${diagramOffset.y}), generatedLegendPosition=${generatedLegendPosition.value}`
    );
    foundElements.forEach((el, idx) => {
      const transform = parseTransform(el.getAttribute('transform'));
      console.log(`[DEBUG] setupDiagramDrag element ${idx} (${el.id}): DOM transform=(${transform.x}, ${transform.y})`);
      if (preserveOffset && (diagramOffset.x !== 0 || diagramOffset.y !== 0)) {
        const adjusted = {
          x: transform.x - diagramOffset.x,
          y: transform.y - diagramOffset.y
        };
        console.log(`[DEBUG] setupDiagramDrag element ${idx}: adjusted to (${adjusted.x}, ${adjusted.y})`);
        originalTransforms.set(el, adjusted);
      } else {
        originalTransforms.set(el, transform);
      }
    });
    diagramElementOriginalTransforms.value = originalTransforms;
    console.log(`[DEBUG] setupDiagramDrag FINISHED: Set ${originalTransforms.size} original transforms`);
    originalTransforms.forEach((transform, el) => {
      console.log(`[DEBUG]   -> ${el.id}: (${transform.x}, ${transform.y})`);
    });

    foundElements.forEach((el) => {
      el.style.cursor = 'grab';
    });

    svg.removeEventListener('mousedown', startDiagramDrag);
    svg.addEventListener('mousedown', startDiagramDrag);
  };

  return {
    applyDiagramShift,
    endDiagramDrag,
    onDiagramDrag,
    resetDiagramPosition,
    setupDiagramDrag,
    startDiagramDrag
  };
};
