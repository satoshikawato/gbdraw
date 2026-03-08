import { getTransformedBBox, parseTransform } from './legend-layout/transform-utils.js';

export const createResultsManager = ({ state, getPyodide }) => {
  const {
    pyodideReady,
    svgContent,
    mode,
    svgContainer,
    cInputType,
    linearSeqs,
    form,
    adv,
    selectedResultIndex,
    results,
    skipCaptureBaseConfig,
    selectedPalette,
    currentColors
  } = state;

  let definitionUpdateTimeout = null;

  const updatePalette = () => {
    const pyodide = getPyodide();
    if (!pyodide) return;
    const paletteJson = pyodide.runPython('get_palettes_json()');
    const all = JSON.parse(paletteJson);
    currentColors.value = { ...(all[selectedPalette.value] || {}) };
  };

  const resetColors = () => updatePalette();

  const parseMixedContentText = (inputText) => {
    const parts = [];
    try {
      const wrapped = `<root>${inputText}</root>`;
      const doc = new DOMParser().parseFromString(wrapped, 'application/xml');
      if (doc.getElementsByTagName('parsererror').length) {
        throw new Error('Parse error');
      }
      const root = doc.documentElement;
      const elements = Array.from(root.children);
      if (elements.length) {
        elements.forEach((el) => {
          parts.push({ text: el.textContent || '', italic: el.tagName.toLowerCase() === 'i' });
          const tail = el.nextSibling;
          if (tail && tail.nodeType === Node.TEXT_NODE) {
            parts.push({ text: tail.nodeValue || '', italic: false });
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

  const getLocalVerticalBounds = (group) => {
    let minY = Number.POSITIVE_INFINITY;
    let maxY = Number.NEGATIVE_INFINITY;
    const texts = Array.from(group?.children || [])
      .filter((child) => String(child?.tagName || '').toLowerCase() === 'text');
    texts.forEach((textEl) => {
      const yValue = Number(String(textEl.getAttribute('y') || '0').replace('px', ''));
      const fontSize = Number(String(textEl.getAttribute('font-size') || '0').replace('px', ''));
      const halfHeight = Number.isFinite(fontSize) && fontSize > 0 ? 0.5 * fontSize : 0;
      minY = Math.min(minY, yValue - halfHeight);
      maxY = Math.max(maxY, yValue + halfHeight);
    });
    if (!Number.isFinite(minY) || !Number.isFinite(maxY)) {
      return { minY: 0, maxY: 0, height: 0 };
    }
    return { minY, maxY, height: maxY - minY };
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
    const normalizedPosition = ['center', 'top', 'bottom'].includes(position) ? position : 'center';
    const bounds = getLocalVerticalBounds(group);
    let y = 0.5 * canvasHeight;
    if (normalizedPosition === 'top') {
      y = 24 - bounds.minY;
    } else if (normalizedPosition === 'bottom') {
      y = canvasHeight - 24 - bounds.maxY;
    }
    setGroupTranslate(group, 0.5 * canvasWidth, y);
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
            isMultiRecordCanvasOnSvg
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
            placeDefinitionGroup(plotTitleGroup, canvasWidth, canvasHeight, normalizedPlotTitlePosition);
            updated = true;
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
    if (definitionUpdateTimeout) clearTimeout(definitionUpdateTimeout);
    definitionUpdateTimeout = setTimeout(updateDefinitionText, 500);
  };

  return {
    updatePalette,
    resetColors,
    scheduleDefinitionUpdate
  };
};
