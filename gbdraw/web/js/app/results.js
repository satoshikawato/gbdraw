import { parseTransform } from './legend-layout/transform-utils.js';

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

      try {
        const species = form.species || '';
        const strain = form.strain || '';
        const fontSize = adv.def_font_size || 18;

        const resultJson = pyodide.runPython(
          `regenerate_definition_svg("${gbPath}", ${
            species ? `"${species.replace(/"/g, '\\"')}"` : 'None'
          }, ${strain ? `"${strain.replace(/"/g, '\\"')}"` : 'None'}, ${fontSize})`
        );
        const result = JSON.parse(resultJson);

        if (result.error) {
          console.error('Definition update error:', result.error);
          return;
        }

        const definitionGroupId = result.definition_group_id;
        console.log('Looking for definition group:', definitionGroupId);
        const existingGroup = svg.getElementById(definitionGroupId);

        if (existingGroup) {
          const parser = new DOMParser();
          const newDoc = parser.parseFromString(
            `<svg xmlns="http://www.w3.org/2000/svg">${result.svg}</svg>`,
            'image/svg+xml'
          );
          const newGroup = newDoc.querySelector('g');

          if (newGroup) {
            const existingTransform = existingGroup.getAttribute('transform');
            if (existingTransform) {
              newGroup.setAttribute('transform', existingTransform);
            }

            existingGroup.parentNode.replaceChild(svg.ownerDocument.importNode(newGroup, true), existingGroup);

            skipCaptureBaseConfig.value = true;
            const idx = selectedResultIndex.value;
            if (idx >= 0 && results.value.length > idx) {
              const serializer = new XMLSerializer();
              results.value[idx] = { ...results.value[idx], content: serializer.serializeToString(svg) };
            }

            console.log('Definition text updated');
          }
        } else {
          console.log('Definition group not found in SVG:', definitionGroupId);
        }
      } catch (e) {
        console.error('Failed to update definition text:', e);
      }
      return;
    }

    if (mode.value === 'linear') {
      const groups = Array.from(svg.querySelectorAll('g[id]'))
        .filter((group) => {
          const id = group.getAttribute('id');
          if (!id) return false;
          if (
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
