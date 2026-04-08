import {
  getDefinitionGroupTranslate,
  getLocalVerticalBounds,
  getTransformedBBox,
  parseTransform
} from './legend-layout/transform-utils.js';

export const createResultsManager = ({ state, getPyodide, legendLayout, rerenderLinearDefinitions = null }) => {
  const {
    pyodideReady,
    svgContent,
    mode,
    shouldDeferCircularPreviewUpdates,
    svgContainer,
    cInputType,
    linearSeqs,
    form,
    adv,
    selectedResultIndex,
    results,
    skipCaptureBaseConfig,
    paletteDefinitions,
    selectedPalette,
    currentColors,
    paletteInstantPreviewEnabled,
    appliedPaletteName,
    appliedPaletteColors,
    pendingPaletteName,
    pendingPaletteColors
  } = state;
  const { clearPlotTitleState, setPlotTitleAutoTransform } = legendLayout;

  let definitionUpdateTimeout = null;
  const cloneColors = (colors) => ({ ...(colors || {}) });
  const getPaletteMap = () => {
    if (paletteDefinitions.value && Object.keys(paletteDefinitions.value).length > 0) {
      return paletteDefinitions.value;
    }

    const pyodide = getPyodide();
    if (!pyodide) return {};

    const paletteJson = pyodide.runPython('get_palettes_json()');
    const all = JSON.parse(paletteJson);
    const normalized = { ...(all || {}) };
    delete normalized.title;
    paletteDefinitions.value = normalized;
    return normalized;
  };
  const getPaletteBaseColors = (paletteName) => {
    const allPalettes = getPaletteMap();
    return cloneColors(allPalettes[paletteName] || {});
  };
  const setAppliedPaletteState = (paletteName, colors = currentColors.value) => {
    appliedPaletteName.value = String(paletteName || selectedPalette.value || 'default');
    appliedPaletteColors.value = cloneColors(colors);
  };
  const setPendingPaletteState = (paletteName, colors = currentColors.value) => {
    pendingPaletteName.value = String(paletteName || selectedPalette.value || '');
    pendingPaletteColors.value = cloneColors(colors);
  };
  const clearPendingPaletteDraft = () => {
    pendingPaletteName.value = '';
    pendingPaletteColors.value = {};
  };
  const applyPaletteDraftToPreview = () => {
    setAppliedPaletteState(selectedPalette.value, currentColors.value);
    clearPendingPaletteDraft();
  };
  const syncPaletteDraftState = () => {
    if (paletteInstantPreviewEnabled.value) {
      applyPaletteDraftToPreview();
      return;
    }

    if (String(pendingPaletteName.value || '').trim() !== '') {
      setPendingPaletteState(selectedPalette.value, currentColors.value);
      return;
    }

    setAppliedPaletteState(selectedPalette.value, currentColors.value);
  };

  const cancelDefinitionUpdate = () => {
    if (definitionUpdateTimeout) {
      clearTimeout(definitionUpdateTimeout);
      definitionUpdateTimeout = null;
    }
  };

  const updatePalette = () => {
    const selectedName = String(selectedPalette.value || '').trim() || 'default';

    if (!paletteInstantPreviewEnabled.value && selectedName === appliedPaletteName.value) {
      currentColors.value = cloneColors(appliedPaletteColors.value);
      clearPendingPaletteDraft();
      return;
    }

    currentColors.value = getPaletteBaseColors(selectedName);
    if (paletteInstantPreviewEnabled.value) {
      applyPaletteDraftToPreview();
      return;
    }

    setPendingPaletteState(selectedName, currentColors.value);
  };

  const resetColors = () => {
    const selectedName = String(selectedPalette.value || '').trim() || 'default';
    currentColors.value = getPaletteBaseColors(selectedName);
    if (paletteInstantPreviewEnabled.value) {
      applyPaletteDraftToPreview();
      return;
    }

    if (String(pendingPaletteName.value || '').trim() !== '') {
      setPendingPaletteState(selectedName, currentColors.value);
      return;
    }

    setAppliedPaletteState(selectedName, currentColors.value);
  };

  const parseMixedContentText = (inputText) => {
    const parts = [];
    try {
      const wrapped = `<root>${inputText}</root>`;
      const doc = new DOMParser().parseFromString(wrapped, 'application/xml');
      if (doc.getElementsByTagName('parsererror').length) {
        throw new Error('Parse error');
      }
      const root = doc.documentElement;
      const nodes = Array.from(root.childNodes);
      if (nodes.length) {
        nodes.forEach((node) => {
          if (node.nodeType === Node.TEXT_NODE) {
            parts.push({ text: node.nodeValue || '', italic: false });
            return;
          }
          if (node.nodeType === Node.ELEMENT_NODE) {
            parts.push({
              text: node.textContent || '',
              italic: node.tagName.toLowerCase() === 'i'
            });
          }
        });
      } else {
        parts.push({ text: root.textContent || '', italic: false });
      }
    } catch (e) {
      parts.push({ text: inputText, italic: false });
    }
    return parts;
  };

  const applyMixedText = (textEl, rawText) => {
    while (textEl.firstChild) {
      textEl.removeChild(textEl.firstChild);
    }
    const parts = parseMixedContentText(rawText);
    parts.forEach((part) => {
      const tspan = document.createElementNS('http://www.w3.org/2000/svg', 'tspan');
      tspan.textContent = part.text || '';
      if (part.italic) {
        tspan.setAttribute('font-style', 'italic');
      }
      textEl.appendChild(tspan);
    });
  };

  const parseSvgCanvasSize = (svg) => {
    const viewBox = String(svg.getAttribute('viewBox') || '').trim().split(/\s+/);
    if (viewBox.length === 4) {
      const width = Number(viewBox[2]);
      const height = Number(viewBox[3]);
      if (Number.isFinite(width) && Number.isFinite(height) && width > 0 && height > 0) {
        return { width, height };
      }
    }
    const width = Number(svg.getAttribute('width'));
    const height = Number(svg.getAttribute('height'));
    if (Number.isFinite(width) && Number.isFinite(height) && width > 0 && height > 0) {
      return { width, height };
    }
    return null;
  };

  const setGroupTranslate = (group, x, y) => {
    if (!group) return;
    group.setAttribute('transform', `translate(${x}, ${y})`);
  };

  const syncSvgCanvasSize = (svg, width, height) => {
    svg.setAttribute('width', `${width}px`);
    svg.setAttribute('height', `${height}px`);
    svg.setAttribute('viewBox', `0 0 ${width} ${height}`);
  };

  const translateTopLevelGroups = (svg, dy, excludedIds = new Set()) => {
    if (!Number.isFinite(dy) || Math.abs(dy) <= 1e-6) return;
    Array.from(svg.children).forEach((child) => {
      const tagName = String(child.tagName || '').toLowerCase();
      if (tagName === 'defs') return;
      const childId = String(child.getAttribute('id') || '');
      if (excludedIds.has(childId)) return;
      const current = parseTransform(child.getAttribute('transform'));
      setGroupTranslate(child, current.x, current.y + dy);
    });
  };

  const placeDefinitionGroup = (group, canvasWidth, canvasHeight, position) => {
    const nextTransform = getDefinitionGroupTranslate(group, canvasWidth, canvasHeight, position);
    setGroupTranslate(group, nextTransform.x, nextTransform.y);
  };

  const parseGroupSvg = (svgMarkup) => {
    const parser = new DOMParser();
    const doc = parser.parseFromString(
      `<svg xmlns="http://www.w3.org/2000/svg">${svgMarkup}</svg>`,
      'image/svg+xml'
    );
    return doc.querySelector('g');
  };

  const updateDefinitionText = async () => {
    if (!svgContent.value) return;
    if (!svgContainer.value) return;
    const svg = svgContainer.value.querySelector('svg');
    if (!svg) return;
    if (mode.value === 'circular' && shouldDeferCircularPreviewUpdates.value) return;

    if (mode.value === 'circular') {
      if (!pyodideReady.value) return;
      const pyodide = getPyodide();
      if (!pyodide) return;

      const gbPath = cInputType.value === 'gb' ? '/input.gb' : '/input.gb';
      const isMultiRecordCanvasOnSvg = svg.querySelector('g[id^="record_"]') !== null;

      try {
        const species = form.species || '';
        const strain = form.strain || '';
        const hasDefinitionFontSize =
          adv.def_font_size !== null && adv.def_font_size !== undefined && adv.def_font_size !== '';
        const definitionFontSize = hasDefinitionFontSize ? Number(adv.def_font_size) : null;
        const hasPlotTitleFontSize =
          adv.plot_title_font_size !== null &&
          adv.plot_title_font_size !== undefined &&
          adv.plot_title_font_size !== '';
        const plotTitleFontSize = hasPlotTitleFontSize ? Number(adv.plot_title_font_size) : null;

        const normalizedSpecies = species === '' ? null : species;
        const normalizedStrain = strain === '' ? null : strain;
        const normalizedPlotTitle = String(form.plot_title || '').trim();
        const rawCircularPlotTitlePosition = String(adv.plot_title_position || 'none').trim().toLowerCase();
        const normalizedPlotTitlePosition = ['none', 'top', 'bottom'].includes(rawCircularPlotTitlePosition)
          ? rawCircularPlotTitlePosition
          : 'none';
        const normalizedDefinitionFontSize =
          definitionFontSize === null || Number.isNaN(definitionFontSize) ? null : definitionFontSize;
        const normalizedPlotTitleFontSize =
          plotTitleFontSize === null || Number.isNaN(plotTitleFontSize) ? null : plotTitleFontSize;
        const keepFullDefinitionWithPlotTitle = Boolean(adv.keep_full_definition_with_plot_title);

        let resultJson = '';
        let regenerateDefinitionSvgs = null;
        try {
          regenerateDefinitionSvgs = pyodide.globals.get('regenerate_definition_svgs');
          resultJson = regenerateDefinitionSvgs(
            gbPath,
            normalizedSpecies,
            normalizedStrain,
            normalizedPlotTitle === '' ? null : normalizedPlotTitle,
            normalizedDefinitionFontSize,
            normalizedPlotTitleFontSize,
            normalizedPlotTitlePosition,
            isMultiRecordCanvasOnSvg,
            keepFullDefinitionWithPlotTitle
          );
        } finally {
          regenerateDefinitionSvgs?.destroy?.();
        }
        const result = JSON.parse(resultJson);

        if (result.error) {
          console.error('Definition update error:', result.error);
          return;
        }

        const definitionEntries = Array.isArray(result.definitions) ? result.definitions : [];
        if (definitionEntries.length === 0) {
          return;
        }

        const desiredGroupIds = new Set(
          definitionEntries
            .map((entry) => String(entry?.definition_group_id || '').trim())
            .filter(Boolean)
        );
        let updated = false;
        definitionEntries.forEach((entry) => {
          const definitionGroupId = entry?.definition_group_id;
          const definitionSvg = entry?.svg;
          if (!definitionGroupId || !definitionSvg) return;

          const newGroup = parseGroupSvg(definitionSvg);
          if (!newGroup) return;
          const existingGroup = svg.getElementById(definitionGroupId);
          const importedGroup = svg.ownerDocument.importNode(newGroup, true);

          if (existingGroup) {
            const existingTransform = existingGroup.getAttribute('transform');
            if (existingTransform) {
              importedGroup.setAttribute('transform', existingTransform);
            }
            existingGroup.parentNode.replaceChild(importedGroup, existingGroup);
            updated = true;
            return;
          }

          if (definitionGroupId === 'plot_title') {
            svg.appendChild(importedGroup);
            updated = true;
          }
        });

        const stalePlotTitleGroup = svg.getElementById('plot_title');
        if (stalePlotTitleGroup && !desiredGroupIds.has('plot_title')) {
          stalePlotTitleGroup.remove();
          updated = true;
        }

        const canvasSize = parseSvgCanvasSize(svg);
        if (canvasSize) {
          let canvasWidth = canvasSize.width;
          let canvasHeight = canvasSize.height;
          const legendPosition = String(form.legend || 'right').trim().toLowerCase();
          const legendGroup = svg.getElementById('legend');
          const singleDefinitionGroup = Array.from(svg.querySelectorAll('g[id$="_definition"]'))
            .find((group) => !group.closest('g[id^="record_"]'));
          if (singleDefinitionGroup) {
            placeDefinitionGroup(singleDefinitionGroup, canvasWidth, canvasHeight, 'center');
            updated = true;
          }

          const plotTitleGroup = svg.getElementById('plot_title');
          if (plotTitleGroup) {
            const plotTitleBounds = getLocalVerticalBounds(plotTitleGroup);
            if (normalizedPlotTitlePosition === 'top' && legendPosition === 'top' && legendGroup) {
              const legendBounds = getTransformedBBox(legendGroup);
              const requiredLegendTop = 24 + plotTitleBounds.height + 20;
              if (legendBounds && legendBounds.y < requiredLegendTop) {
                const shiftY = requiredLegendTop - legendBounds.y;
                translateTopLevelGroups(svg, shiftY, new Set(['plot_title']));
                canvasHeight += shiftY;
                syncSvgCanvasSize(svg, canvasWidth, canvasHeight);
                updated = true;
              }
            }
            if (normalizedPlotTitlePosition === 'bottom' && legendPosition === 'bottom' && legendGroup) {
              const legendBounds = getTransformedBBox(legendGroup);
              const legendBottom = legendBounds ? legendBounds.y + legendBounds.height : 0;
              const requiredHeight = legendBottom + 20 + plotTitleBounds.height + 24;
              if (requiredHeight > canvasHeight) {
                canvasHeight = requiredHeight;
                syncSvgCanvasSize(svg, canvasWidth, canvasHeight);
                updated = true;
              }
            }
            const nextTitleTransform = getDefinitionGroupTranslate(
              plotTitleGroup,
              canvasWidth,
              canvasHeight,
              normalizedPlotTitlePosition
            );
            setPlotTitleAutoTransform(plotTitleGroup, nextTitleTransform, {
              preserveUserOffset: true
            });
            updated = true;
          } else {
            clearPlotTitleState();
          }
        }

        if (updated) {
          skipCaptureBaseConfig.value = true;
          const idx = selectedResultIndex.value;
          if (idx >= 0 && results.value.length > idx) {
            const serializer = new XMLSerializer();
            results.value[idx] = { ...results.value[idx], content: serializer.serializeToString(svg) };
          }

          console.log('Definition text updated');
        }
      } catch (e) {
        console.error('Failed to update definition text:', e);
      }
      return;
    }

    if (mode.value === 'linear') {
      if (typeof rerenderLinearDefinitions === 'function') {
        await rerenderLinearDefinitions('definition-edit');
        return;
      }

      const plotTitleGroup = svg.getElementById('plot_title');
      const groups = Array.from(svg.querySelectorAll('g[id]'))
        .filter((group) => {
          const id = group.getAttribute('id');
          if (!id) return false;
          if (
            id === 'plot_title' ||
            id === 'legend' ||
            id === 'feature_legend' ||
            id === 'pairwise_legend' ||
            id === 'horizontal_legend' ||
            id === 'vertical_legend' ||
            id === 'length_bar'
          ) {
            return false;
          }
          const hasText = group.querySelectorAll('text').length > 0;
          if (!hasText) return false;
          const hasShapes =
            group.querySelectorAll('path, line, rect, polygon, polyline, circle').length > 0;
          return !hasShapes;
        })
        .sort((a, b) => parseTransform(a.getAttribute('transform')).y - parseTransform(b.getAttribute('transform')).y);

      if (groups.length === 0) return;

      const labels = linearSeqs.map((seq) => (seq.definition ?? '').toString());
      let updated = false;
      const fontSizeOverride =
        adv.def_font_size !== null && adv.def_font_size !== undefined && adv.def_font_size !== ''
          ? String(adv.def_font_size)
          : null;

      groups.forEach((group, idx) => {
        const label = labels[idx] ?? '';
        const texts = Array.from(group.querySelectorAll('text'));
        if (texts.length === 0) return;
        const nextLabel = label.trim() ? label.trim() : group.getAttribute('id') || '';
        applyMixedText(texts[0], nextLabel);
        updated = true;
        if (fontSizeOverride) {
          texts.forEach((t) => {
            if (t.getAttribute('font-size') !== fontSizeOverride) {
              t.setAttribute('font-size', fontSizeOverride);
              updated = true;
            }
          });
        }
      });

      if (plotTitleGroup) {
        const titleText = String(form.plot_title || '').trim();
        const titleTexts = Array.from(plotTitleGroup.querySelectorAll('text'));
        if (titleTexts.length > 0) {
          applyMixedText(titleTexts[0], titleText);
          updated = true;
        }
        const titleFontSize =
          adv.plot_title_font_size !== null &&
          adv.plot_title_font_size !== undefined &&
          adv.plot_title_font_size !== ''
            ? Number(adv.plot_title_font_size)
            : null;
        if (titleFontSize !== null && Number.isFinite(titleFontSize) && titleFontSize > 0) {
          titleTexts.forEach((t) => {
            if (t.getAttribute('font-size') !== String(titleFontSize)) {
              t.setAttribute('font-size', String(titleFontSize));
              updated = true;
            }
          });
        }
        const canvasSize = parseSvgCanvasSize(svg);
        if (canvasSize) {
          const normalizedTitlePosition = String(adv.plot_title_position || 'bottom').trim().toLowerCase();
          const safeTitlePosition = ['center', 'top', 'bottom'].includes(normalizedTitlePosition)
            ? normalizedTitlePosition
            : 'bottom';
          const titleHeightRaw = titleTexts[0]?.getBBox?.().height;
          const titleHeight = Number.isFinite(titleHeightRaw) && titleHeightRaw > 0 ? titleHeightRaw : 32;
          const edgeMargin = 24;
          let titleY = 0.5 * canvasSize.height;
          if (safeTitlePosition === 'top') {
            titleY = edgeMargin + (0.5 * titleHeight);
          } else if (safeTitlePosition === 'bottom') {
            titleY = canvasSize.height - edgeMargin - (0.5 * titleHeight);
          }
          const nextTransform = `translate(${0.5 * canvasSize.width},${titleY})`;
          if (plotTitleGroup.getAttribute('transform') !== nextTransform) {
            plotTitleGroup.setAttribute('transform', nextTransform);
            updated = true;
          }
        }
      }

      if (updated) {
        skipCaptureBaseConfig.value = true;
        const idx = selectedResultIndex.value;
        if (idx >= 0 && results.value.length > idx) {
          const serializer = new XMLSerializer();
          results.value[idx] = { ...results.value[idx], content: serializer.serializeToString(svg) };
        }
      }
    }
  };

  const scheduleDefinitionUpdate = () => {
    cancelDefinitionUpdate();
    if (mode.value === 'circular' && shouldDeferCircularPreviewUpdates.value) return;
    definitionUpdateTimeout = setTimeout(() => {
      definitionUpdateTimeout = null;
      void updateDefinitionText();
    }, 500);
  };

  return {
    updatePalette,
    resetColors,
    applyPaletteDraftToPreview,
    clearPendingPaletteDraft,
    setAppliedPaletteState,
    setPendingPaletteState,
    syncPaletteDraftState,
    scheduleDefinitionUpdate,
    cancelDefinitionUpdate
  };
};
