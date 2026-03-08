import { parseTransform } from './transform-utils.js';

export const createDiagramDragActions = ({ state }) => {
  const {
    svgContainer,
    mode,
    diagramElements,
    diagramElementIds,
    diagramElementOriginalTransforms,
    diagramOffset,
    diagramDragging,
    diagramDragStart,
    plotTitleElement,
    plotTitleDragging,
    plotTitleDragStart,
    plotTitleAutoTransform,
    plotTitleUserOffset,
    zoom,
    generatedLegendPosition
  } = state;

  const LEGEND_GROUP_IDS = new Set([
    'legend',
    'feature_legend',
    'pairwise_legend',
    'horizontal_legend',
    'vertical_legend'
  ]);

  let activeDragElements = [];
  let activeDragMode = 'group'; // 'group' | 'record' | 'plot_title'
  let activeDragOriginalTransforms = new Map();
  let activePlotTitleOffsetStart = { x: 0, y: 0 };

  const isMultiRecordCanvasSvg = (svg) => {
    return Array.from(svg.children).some((el) => {
      if (!el || el.tagName.toLowerCase() !== 'g') return false;
      const id = el.getAttribute('id') || '';
      return id.startsWith('record_');
    });
  };

  const getTopLevelDiagramGroupsForMultiRecord = (svg) => {
    const groups = [];
    const ids = [];

    Array.from(svg.children).forEach((el) => {
      if (!el || el.tagName.toLowerCase() !== 'g') return;
      const id = el.getAttribute('id');
      if (!id || LEGEND_GROUP_IDS.has(id)) return;
      groups.push(el);
      ids.push(id);
    });

    return { groups, ids };
  };

  const getDiagramGroupsForSingleRecordOrLinear = (svg) => {
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

      if (LEGEND_GROUP_IDS.has(id)) return;
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
      if (!id || LEGEND_GROUP_IDS.has(id)) return;
      if (foundElements.includes(el)) return;
      foundElements.push(el);
      if (!foundIds.includes(id)) {
        foundIds.push(id);
      }
    });

    return { groups: foundElements, ids: foundIds };
  };

  const setTranslate = (el, x, y) => {
    if (!el) return;
    el.setAttribute('transform', `translate(${x}, ${y})`);
  };

  const normalizeTransform = (transform) => {
    const x = Number(transform?.x);
    const y = Number(transform?.y);
    return {
      x: Number.isFinite(x) ? x : 0,
      y: Number.isFinite(y) ? y : 0
    };
  };

  const assignPlotTitleElement = (group) => {
    if (plotTitleElement.value && plotTitleElement.value !== group) {
      plotTitleElement.value.style.opacity = '1';
      plotTitleElement.value.style.cursor = '';
    }
    plotTitleElement.value = group || null;
    if (plotTitleElement.value) {
      plotTitleElement.value.style.opacity = '1';
      plotTitleElement.value.style.cursor = 'grab';
    }
  };

  const applyPlotTitleTransform = (group = plotTitleElement.value) => {
    if (mode.value !== 'circular' || !group) return;
    const autoTransform = normalizeTransform(plotTitleAutoTransform.value);
    const nextX = autoTransform.x + plotTitleUserOffset.x;
    const nextY = autoTransform.y + plotTitleUserOffset.y;
    setTranslate(group, nextX, nextY);
  };

  const clearPlotTitleState = () => {
    assignPlotTitleElement(null);
    plotTitleDragging.value = false;
    plotTitleDragStart.x = 0;
    plotTitleDragStart.y = 0;
    plotTitleAutoTransform.value = { x: 0, y: 0 };
    plotTitleUserOffset.x = 0;
    plotTitleUserOffset.y = 0;
    activePlotTitleOffsetStart = { x: 0, y: 0 };
  };

  const setPlotTitleAutoTransform = (group, nextAutoTransform, { preserveUserOffset = true } = {}) => {
    if (mode.value !== 'circular' || !group) {
      clearPlotTitleState();
      return;
    }
    assignPlotTitleElement(group);
    plotTitleAutoTransform.value = normalizeTransform(nextAutoTransform);
    if (!preserveUserOffset) {
      plotTitleUserOffset.x = 0;
      plotTitleUserOffset.y = 0;
    }
    applyPlotTitleTransform(group);
  };

  const resetPlotTitlePosition = () => {
    plotTitleDragging.value = false;
    plotTitleUserOffset.x = 0;
    plotTitleUserOffset.y = 0;
    if (plotTitleElement.value) {
      plotTitleElement.value.style.opacity = '1';
      applyPlotTitleTransform(plotTitleElement.value);
    }
  };

  const syncPlotTitleElement = (svg, preserveOffset = false) => {
    if (mode.value !== 'circular') {
      clearPlotTitleState();
      return;
    }

    const nextPlotTitleGroup = svg?.getElementById('plot_title');
    if (!nextPlotTitleGroup) {
      clearPlotTitleState();
      return;
    }

    const hadPlotTitleState = !!plotTitleElement.value;
    assignPlotTitleElement(nextPlotTitleGroup);
    if (!preserveOffset) {
      plotTitleAutoTransform.value = parseTransform(nextPlotTitleGroup.getAttribute('transform'));
      plotTitleUserOffset.x = 0;
      plotTitleUserOffset.y = 0;
    } else if (!hadPlotTitleState) {
      plotTitleAutoTransform.value = parseTransform(nextPlotTitleGroup.getAttribute('transform'));
    }
    applyPlotTitleTransform(nextPlotTitleGroup);
  };

  const applyDiagramShift = (deltaX, deltaY) => {
    if (diagramElements.value.length === 0) return;
    diagramElements.value.forEach((el) => {
      const current = parseTransform(el.getAttribute('transform'));
      const newX = current.x + deltaX;
      const newY = current.y + deltaY;
      setTranslate(el, newX, newY);
      diagramElementOriginalTransforms.value.set(el, { x: newX, y: newY });
    });
    diagramOffset.x = 0;
    diagramOffset.y = 0;
  };

  const startPlotTitleDrag = (e, group) => {
    if (!group) return;

    e.preventDefault();
    assignPlotTitleElement(group);
    diagramDragging.value = true;
    plotTitleDragging.value = true;
    plotTitleDragStart.x = e.clientX;
    plotTitleDragStart.y = e.clientY;
    activePlotTitleOffsetStart = { x: plotTitleUserOffset.x, y: plotTitleUserOffset.y };
    activeDragMode = 'plot_title';
    activeDragElements = [group];
    activeDragOriginalTransforms = new Map([[group, parseTransform(group.getAttribute('transform'))]]);
    group.style.opacity = '0.8';

    document.addEventListener('mousemove', onDiagramDrag);
    document.addEventListener('mouseup', endDiagramDrag);
  };

  const startDiagramDrag = (e) => {
    if (
      e.target.closest('#legend') ||
      e.target.closest('text[data-label-editable="true"]') ||
      e.target.closest('path[id^="f"], polygon[id^="f"], rect[id^="f"]')
    ) {
      return;
    }
    if (!svgContainer.value) return;
    const svg = svgContainer.value.querySelector('svg');
    if (!svg) return;

    if (mode.value === 'circular') {
      const clickedPlotTitle = e.target.closest('#plot_title');
      if (clickedPlotTitle) {
        startPlotTitleDrag(e, clickedPlotTitle);
        return;
      }
    }

    const isMultiRecordCanvas = isMultiRecordCanvasSvg(svg);
    let dragTargets = [];

    if (isMultiRecordCanvas) {
      // Multi-record mode: allow dragging only the clicked record_* group.
      const clickedRecordGroup = e.target.closest('g[id^="record_"]');
      if (!clickedRecordGroup) return;
      dragTargets = [clickedRecordGroup];
      activeDragMode = 'record';
    } else {
      const clickedGroup = e.target.closest('g[id]');
      if (!clickedGroup) return;

      const clickedId = clickedGroup.id;
      if (LEGEND_GROUP_IDS.has(clickedId)) return;

      dragTargets = diagramElements.value;
      if (dragTargets.length === 0) return;
      activeDragMode = 'group';
    }

    e.preventDefault();
    diagramDragging.value = true;
    plotTitleDragging.value = false;
    diagramDragStart.x = e.clientX;
    diagramDragStart.y = e.clientY;

    activeDragElements = dragTargets;
    activeDragOriginalTransforms = new Map();

    activeDragElements.forEach((el) => {
      if (activeDragMode === 'record') {
        // Keep per-record drag independent from global offsets.
        activeDragOriginalTransforms.set(el, parseTransform(el.getAttribute('transform')));
      } else {
        const original = diagramElementOriginalTransforms.value.get(el) || { x: 0, y: 0 };
        activeDragOriginalTransforms.set(el, original);
      }
      el.style.opacity = '0.8';
    });

    document.addEventListener('mousemove', onDiagramDrag);
    document.addEventListener('mouseup', endDiagramDrag);
  };

  const getActiveDragStart = () => {
    return activeDragMode === 'plot_title' ? plotTitleDragStart : diagramDragStart;
  };

  const onDiagramDrag = (e) => {
    if (!diagramDragging.value || activeDragElements.length === 0) return;

    const dragStart = getActiveDragStart();
    const deltaX = (e.clientX - dragStart.x) / zoom.value;
    const deltaY = (e.clientY - dragStart.y) / zoom.value;

    if (activeDragMode === 'plot_title') {
      plotTitleUserOffset.x = activePlotTitleOffsetStart.x + deltaX;
      plotTitleUserOffset.y = activePlotTitleOffsetStart.y + deltaY;
      applyPlotTitleTransform(activeDragElements[0] || plotTitleElement.value);
      return;
    }

    const includeGlobalOffset = activeDragMode === 'group';

    activeDragElements.forEach((el) => {
      const original = activeDragOriginalTransforms.get(el) || { x: 0, y: 0 };
      const offsetX = includeGlobalOffset ? diagramOffset.x : 0;
      const offsetY = includeGlobalOffset ? diagramOffset.y : 0;
      const newX = original.x + offsetX + deltaX;
      const newY = original.y + offsetY + deltaY;
      setTranslate(el, newX, newY);
    });
  };

  const endDiagramDrag = (e) => {
    if (!diagramDragging.value) return;

    const dragStart = getActiveDragStart();
    const currentX = typeof e?.clientX === 'number' ? e.clientX : dragStart.x;
    const currentY = typeof e?.clientY === 'number' ? e.clientY : dragStart.y;
    const deltaX = (currentX - dragStart.x) / zoom.value;
    const deltaY = (currentY - dragStart.y) / zoom.value;
    if (activeDragMode === 'group') {
      diagramOffset.x += deltaX;
      diagramOffset.y += deltaY;
    }

    diagramDragging.value = false;
    plotTitleDragging.value = false;
    document.removeEventListener('mousemove', onDiagramDrag);
    document.removeEventListener('mouseup', endDiagramDrag);

    activeDragElements.forEach((el) => {
      el.style.opacity = '1';
    });
    activeDragElements = [];
    activeDragOriginalTransforms = new Map();
    activeDragMode = 'group';
    activePlotTitleOffsetStart = { x: plotTitleUserOffset.x, y: plotTitleUserOffset.y };
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
        setTranslate(el, original.x, original.y);
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
    const isMultiRecordCanvas = isMultiRecordCanvasSvg(svg);

    if (!preserveOffset) {
      diagramOffset.x = 0;
      diagramOffset.y = 0;
    }

    const selectedGroups = isMultiRecordCanvas
      ? getTopLevelDiagramGroupsForMultiRecord(svg)
      : getDiagramGroupsForSingleRecordOrLinear(svg);
    const separatePlotTitle = mode.value === 'circular';
    const foundElements = selectedGroups.groups.filter((el) => !(separatePlotTitle && (el.id || '') === 'plot_title'));
    const foundIds = selectedGroups.ids.filter((id) => !(separatePlotTitle && id === 'plot_title'));

    diagramElements.value = foundElements;
    diagramElementIds.value = foundIds;
    syncPlotTitleElement(svg, preserveOffset);

    const previousOriginalTransformsById = new Map();
    if (preserveOffset && isMultiRecordCanvas && diagramElementOriginalTransforms.value.size > 0) {
      diagramElementOriginalTransforms.value.forEach((transform, mapEl) => {
        const id = mapEl?.id || '';
        if (!id) return;
        if (!previousOriginalTransformsById.has(id)) {
          previousOriginalTransformsById.set(id, []);
        }
        previousOriginalTransformsById.get(id).push(transform);
      });
    }

    const remapCounters = new Map();
    const originalTransforms = new Map();
    console.log(`[DEBUG] ========== setupDiagramDrag CALL #${setupDiagramDragCallCount} ==========`);
    console.log(
      `[DEBUG] setupDiagramDrag: preserveOffset=${preserveOffset}, offset=(${diagramOffset.x}, ${diagramOffset.y}), generatedLegendPosition=${generatedLegendPosition.value}, isMultiRecordCanvas=${isMultiRecordCanvas}`
    );
    foundElements.forEach((el, idx) => {
      if (preserveOffset && isMultiRecordCanvas) {
        const id = el.id || '';
        const preserved = previousOriginalTransformsById.get(id);
        if (preserved && preserved.length > 0) {
          const preservedIdx = remapCounters.get(id) || 0;
          if (preservedIdx < preserved.length) {
            const preservedTransform = preserved[preservedIdx];
            remapCounters.set(id, preservedIdx + 1);
            originalTransforms.set(el, preservedTransform);
            console.log(
              `[DEBUG] setupDiagramDrag element ${idx} (${id}): remapped preserved original=(${preservedTransform.x}, ${preservedTransform.y})`
            );
            return;
          }
        }
      }

      const transform = parseTransform(el.getAttribute('transform'));
      console.log(`[DEBUG] setupDiagramDrag element ${idx} (${el.id}): DOM transform=(${transform.x}, ${transform.y})`);
      if (!isMultiRecordCanvas && preserveOffset && (diagramOffset.x !== 0 || diagramOffset.y !== 0)) {
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
      if (isMultiRecordCanvas) {
        el.style.cursor = (el.id || '').startsWith('record_') ? 'grab' : '';
      } else {
        el.style.cursor = 'grab';
      }
    });

    svg.removeEventListener('mousedown', startDiagramDrag);
    svg.addEventListener('mousedown', startDiagramDrag);
  };

  return {
    applyDiagramShift,
    clearPlotTitleState,
    endDiagramDrag,
    onDiagramDrag,
    resetDiagramPosition,
    resetPlotTitlePosition,
    setPlotTitleAutoTransform,
    setupDiagramDrag,
    startDiagramDrag
  };
};
