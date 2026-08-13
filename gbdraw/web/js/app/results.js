import { parseTransform } from './legend-layout/transform-utils.js';
import { COMPOSITION_ROLE_ATTRIBUTE } from './legend-layout/composition-actions.js';
import { COMPARISON_LEGEND_SELECTOR } from './legend/utils.js';
import { isMultiRecordCanvasSvg } from './record-groups.js';
import { serializeCleanSvg } from '../services/svg-serialization.js';

const RECORD_DEFINITION_SELECTOR = 'g[data-gbdraw-role="record-definition"]';

const findRecordDefinitionGroup = (svg, entry = {}) => {
  const rawRecordIndex = entry?.record_index ?? entry?.recordIndex;
  const hasRecordIndex =
    rawRecordIndex !== null &&
    rawRecordIndex !== undefined &&
    String(rawRecordIndex).trim() !== '';

  if (hasRecordIndex) {
    const semanticMatch = Array.from(
      svg?.querySelectorAll?.(RECORD_DEFINITION_SELECTOR) || []
    ).find(
      (group) =>
        String(group.getAttribute?.('data-gbdraw-record-index') || '') ===
        String(rawRecordIndex)
    );
    if (semanticMatch) return semanticMatch;
  }

  const definitionGroupId = String(entry?.definition_group_id || '').trim();
  return definitionGroupId ? svg?.getElementById?.(definitionGroupId) || null : null;
};

const preserveDefinitionGroupDomIdentity = (existingGroup, importedGroup) => {
  if (!existingGroup || !importedGroup) return importedGroup;

  const existingId = existingGroup.getAttribute?.('id');
  if (existingId !== null && existingId !== undefined) {
    importedGroup.setAttribute('id', existingId);
  }
  const existingTransform = existingGroup.getAttribute?.('transform');
  if (existingTransform) {
    importedGroup.setAttribute('transform', existingTransform);
  }
  const compositionRole = existingGroup.getAttribute?.(COMPOSITION_ROLE_ATTRIBUTE);
  if (compositionRole) {
    importedGroup.setAttribute(COMPOSITION_ROLE_ATTRIBUTE, compositionRole);
  }
  return importedGroup;
};

const stageCircularDefinitionSource = async ({
  inputType,
  inputFile,
  writeFileToFs,
  path = '/definition-input.gb'
}) => {
  if (String(inputType || '').trim().toLowerCase() !== 'gb' || !inputFile) return null;
  if (typeof writeFileToFs !== 'function') {
    throw new Error('Circular definition input staging is unavailable.');
  }
  if (!await writeFileToFs(inputFile, path)) {
    throw new Error('Circular definition input could not be staged.');
  }
  return path;
};

export const createResultsManager = ({
  state,
  getPyodide,
  ensurePyodide = null,
  writeFileToFs = null,
  legendLayout,
  rerenderLinearDefinitions = null,
  bulkStyleActions = null
}) => {
  const {
    pyodideReady,
    svgContent,
    mode,
    shouldDeferCircularPreviewUpdates,
    svgContainer,
    cInputType,
    files,
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
    pendingPaletteColors,
    normalizePaletteColors,
    normalizePaletteDefinitions
  } = state;
  const { refreshCompositionGeometry } = legendLayout;

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
    const normalized = normalizePaletteDefinitions(all || {});
    paletteDefinitions.value = normalized;
    return normalized;
  };
  const getPaletteBaseColors = (paletteName) => {
    const allPalettes = getPaletteMap();
    return normalizePaletteColors(cloneColors(allPalettes[paletteName] || {}));
  };
  const setPendingPaletteState = (paletteName, colors = currentColors.value) => {
    pendingPaletteName.value = String(paletteName || selectedPalette.value || '');
    pendingPaletteColors.value = cloneColors(colors);
  };
  const clearPendingPaletteDraft = () => {
    pendingPaletteName.value = '';
    pendingPaletteColors.value = {};
  };
  const paletteColorSignature = (colors) => JSON.stringify(
    Object.entries(normalizePaletteColors(cloneColors(colors)))
      .map(([key, value]) => [String(key), String(value ?? '').trim().toLowerCase()])
      .filter(([key, value]) => key && value)
      .sort(([left], [right]) => left.localeCompare(right))
  );
  const paletteChangeIsApplied = (paletteName, colors) => (
    String(paletteName || 'default') === String(appliedPaletteName.value || 'default')
    && paletteColorSignature(colors) === paletteColorSignature(appliedPaletteColors.value)
  );
  const convergeAppliedPaletteUi = (paletteName, colors) => {
    const name = String(paletteName || 'default');
    if (String(selectedPalette.value || '') !== name) selectedPalette.value = name;
    if (paletteColorSignature(currentColors.value) !== paletteColorSignature(colors)) {
      currentColors.value = cloneColors(colors);
    }
    if (
      String(pendingPaletteName.value || '').trim() !== ''
      || Object.keys(pendingPaletteColors.value || {}).length > 0
    ) clearPendingPaletteDraft();
  };
  const requestAppliedPaletteChange = async (paletteName, colors, label) => {
    if (paletteChangeIsApplied(paletteName, colors)) {
      convergeAppliedPaletteUi(paletteName, colors);
      return false;
    }
    if (typeof bulkStyleActions?.requestFeatureBulkStyleChange !== 'function') {
      throw new Error('Bulk Feature style actions are unavailable.');
    }
    return bulkStyleActions.requestFeatureBulkStyleChange({
      writerKind: 'palette-default-acceptance',
      label,
      replacement: {
        appliedPaletteName: String(paletteName || 'default'),
        appliedPaletteColors: cloneColors(colors)
      }
    });
  };

  const applyPaletteDraftToPreview = async () => {
    const committed = await requestAppliedPaletteChange(
      selectedPalette.value,
      currentColors.value,
      'Apply feature color palette'
    );
    return committed;
  };
  const syncPaletteDraftState = async () => {
    if (paletteInstantPreviewEnabled.value) {
      return applyPaletteDraftToPreview();
    }

    if (String(pendingPaletteName.value || '').trim() !== '') {
      setPendingPaletteState(selectedPalette.value, currentColors.value);
      return false;
    }

    return requestAppliedPaletteChange(
      selectedPalette.value,
      currentColors.value,
      'Change default feature color'
    );
  };

  const acceptDefaultColor = async (featureType, value) => {
    const key = String(featureType || '').trim();
    if (!key) return false;
    const previous = cloneColors(currentColors.value);
    currentColors.value = { ...previous, [key]: value };
    const isLiveAcceptance = paletteInstantPreviewEnabled.value
      || String(pendingPaletteName.value || '').trim() === '';
    try {
      const committed = await syncPaletteDraftState();
      if (isLiveAcceptance && !committed && !paletteChangeIsApplied(
        selectedPalette.value,
        currentColors.value
      )) {
        currentColors.value = previous;
      }
      return committed;
    } catch (error) {
      if (isLiveAcceptance) currentColors.value = previous;
      throw error;
    }
  };

  const addCustomColor = () => acceptDefaultColor(
    state.newColorFeat?.value,
    state.newColorVal?.value
  );

  const handlePaletteInstantPreviewChange = async () => {
    if (!paletteInstantPreviewEnabled.value) return false;
    if (String(pendingPaletteName.value || '').trim() === '') return false;
    return applyPaletteDraftToPreview();
  };

  const cancelDefinitionUpdate = () => {
    if (definitionUpdateTimeout) {
      clearTimeout(definitionUpdateTimeout);
      definitionUpdateTimeout = null;
    }
  };

  const updatePalette = async () => {
    const selectedName = String(selectedPalette.value || '').trim() || 'default';

    if (!paletteInstantPreviewEnabled.value && selectedName === appliedPaletteName.value) {
      currentColors.value = cloneColors(appliedPaletteColors.value);
      clearPendingPaletteDraft();
      return;
    }

    const nextColors = getPaletteBaseColors(selectedName);
    if (paletteInstantPreviewEnabled.value) {
      const committed = await requestAppliedPaletteChange(
        selectedName,
        nextColors,
        'Apply feature color palette'
      );
      return committed;
    }

    currentColors.value = nextColors;
    setPendingPaletteState(selectedName, currentColors.value);
    return false;
  };

  const resetColors = async () => {
    const selectedName = String(selectedPalette.value || '').trim() || 'default';
    const nextColors = getPaletteBaseColors(selectedName);
    if (paletteInstantPreviewEnabled.value) {
      const committed = await requestAppliedPaletteChange(
        selectedName,
        nextColors,
        'Reset default feature colors'
      );
      return committed;
    }

    currentColors.value = nextColors;
    if (String(pendingPaletteName.value || '').trim() !== '') {
      setPendingPaletteState(selectedName, currentColors.value);
      return false;
    }

    return requestAppliedPaletteChange(
      selectedName,
      currentColors.value,
      'Reset default feature colors'
    );
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
      if (!pyodideReady.value && typeof ensurePyodide === 'function') {
        await ensurePyodide();
      }
      if (!pyodideReady.value) return;
      const pyodide = getPyodide();
      if (!pyodide) return;

      const isMultiRecordCanvasOnSvg = isMultiRecordCanvasSvg(svg);

      try {
        const gbPath = await stageCircularDefinitionSource({
          inputType: cInputType.value,
          inputFile: files?.c_gb,
          writeFileToFs
        });
        if (!gbPath) return;
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
          definitionFontSize !== null && Number.isFinite(definitionFontSize) && definitionFontSize > 0
            ? definitionFontSize
            : null;
        const normalizedPlotTitleFontSize =
          plotTitleFontSize !== null && Number.isFinite(plotTitleFontSize) && plotTitleFontSize > 0
            ? plotTitleFontSize
            : null;
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
          const existingGroup = findRecordDefinitionGroup(svg, entry);
          const importedGroup = svg.ownerDocument.importNode(newGroup, true);

          if (existingGroup) {
            preserveDefinitionGroupDomIdentity(existingGroup, importedGroup);
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

        const plotTitleGroup = svg.getElementById('plot_title');
        refreshCompositionGeometry({
          titleSide: plotTitleGroup ? normalizedPlotTitlePosition : 'none',
          titleTarget: plotTitleGroup
        });

        if (updated) {
          skipCaptureBaseConfig.value = true;
          const idx = selectedResultIndex.value;
          if (idx >= 0 && results.value.length > idx) {
            results.value[idx] = { ...results.value[idx], content: serializeCleanSvg(svg) };
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
            group.matches?.(COMPARISON_LEGEND_SELECTOR) ||
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
      const parsedDefinitionFontSize =
        adv.def_font_size !== null && adv.def_font_size !== undefined && adv.def_font_size !== ''
          ? Number(adv.def_font_size)
          : null;
      const fontSizeOverride =
        parsedDefinitionFontSize !== null &&
        Number.isFinite(parsedDefinitionFontSize) &&
        parsedDefinitionFontSize > 0
          ? String(parsedDefinitionFontSize)
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
        const normalizedTitlePosition = String(adv.plot_title_position || 'bottom').trim().toLowerCase();
        const safeTitlePosition = ['center', 'top', 'bottom'].includes(normalizedTitlePosition)
          ? normalizedTitlePosition
          : 'bottom';
        refreshCompositionGeometry({
          titleSide: safeTitlePosition,
          titleTarget: plotTitleGroup
        });
        updated = true;
      }

      if (updated) {
        skipCaptureBaseConfig.value = true;
        const idx = selectedResultIndex.value;
        if (idx >= 0 && results.value.length > idx) {
          results.value[idx] = { ...results.value[idx], content: serializeCleanSvg(svg) };
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
    acceptDefaultColor,
    addCustomColor,
    handlePaletteInstantPreviewChange,
    applyPaletteDraftToPreview,
    clearPendingPaletteDraft,
    setPendingPaletteState,
    syncPaletteDraftState,
    scheduleDefinitionUpdate,
    cancelDefinitionUpdate
  };
};
