export const createResultsManager = ({ state, getPyodide }) => {
  const {
    pyodideReady,
    svgContent,
    mode,
    svgContainer,
    cInputType,
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

  const updateDefinitionText = async () => {
    if (!pyodideReady.value || !svgContent.value || mode.value !== 'circular') return;
    if (!svgContainer.value) return;
    const pyodide = getPyodide();
    if (!pyodide) return;

    const svg = svgContainer.value.querySelector('svg');
    if (!svg) return;

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
