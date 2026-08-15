import { parseTransform } from './legend-layout/transform-utils.js';
import { COMPOSITION_ROLE_ATTRIBUTE } from './legend-layout/composition-actions.js';
import { COMPARISON_LEGEND_SELECTOR } from './legend/utils.js';
import { isMultiRecordCanvasSvg } from './record-groups.js';
import {
  DIAGRAM_HELPER_OPERATIONS,
  runDiagramHelperOperation
} from '../services/diagram-generation.js';
import { cloneFileBytesForTransfer } from '../services/file-content-cache.js';

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

export const createResultsManager = ({
  state,
  legendLayout,
  rerenderLinearDefinitions = null,
  previewRuntime
}) => {
  if (typeof previewRuntime?.commitDomEdit !== 'function') {
    throw new Error('createResultsManager requires the preview runtime edit protocol.');
  }
  const {
    svgContent,
    mode,
    shouldDeferCircularPreviewUpdates,
    svgContainer,
    cInputType,
    files,
    linearSeqs,
    form,
    adv,
    paletteDefinitions,
    selectedPalette,
    currentColors,
    paletteInstantPreviewEnabled,
    appliedPaletteName,
    appliedPaletteColors,
    pendingPaletteName,
    pendingPaletteColors,
    normalizePaletteColors
  } = state;
  const { refreshCompositionGeometry } = legendLayout;
  const commitDefinitionEdit = () => previewRuntime.commitDomEdit({
    reason: 'definition-text',
    invalidateIndexes: ['legend'],
    mutate: () => true
  });

  let definitionUpdateTimeout = null;
  const cloneColors = (colors) => ({ ...(colors || {}) });
  const getPaletteMap = () => {
    if (paletteDefinitions.value && Object.keys(paletteDefinitions.value).length > 0) {
      return paletteDefinitions.value;
    }
    return {};
  };
  const getPaletteBaseColors = (paletteName) => {
    const allPalettes = getPaletteMap();
    return normalizePaletteColors(cloneColors(allPalettes[paletteName] || {}));
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
      const isMultiRecordCanvasOnSvg = isMultiRecordCanvasSvg(svg);

      try {
        if (String(cInputType.value || '').trim().toLowerCase() !== 'gb' || !files?.c_gb) return;
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

        const response = await runDiagramHelperOperation(
          DIAGRAM_HELPER_OPERATIONS.REGENERATE_DEFINITION_SVGS,
          {
            files: [{
              role: 'source',
              bytes: await cloneFileBytesForTransfer(files.c_gb)
            }],
            species: normalizedSpecies,
            strain: normalizedStrain,
            plotTitle: normalizedPlotTitle === '' ? null : normalizedPlotTitle,
            definitionFontSize: normalizedDefinitionFontSize,
            plotTitleFontSize: normalizedPlotTitleFontSize,
            plotTitlePosition: normalizedPlotTitlePosition,
            multiRecordCanvas: isMultiRecordCanvasOnSvg,
            keepFullDefinitionWithPlotTitle
          }
        );
        const result = response.result;

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
          commitDefinitionEdit();
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
        commitDefinitionEdit();
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
