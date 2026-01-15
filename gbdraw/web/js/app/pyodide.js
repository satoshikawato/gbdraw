import { GBDRAW_WHEEL_NAME } from '../config.js';
import { PYTHON_HELPERS } from './python-helpers.js';

export const createPyodideManager = ({ state }) => {
  const { pyodideReady, loadingStatus, paletteNames, currentColors } = state;
  const pyodideRef = { current: null };

  const getPyodide = () => pyodideRef.current;
  const setPyodide = (value) => {
    pyodideRef.current = value;
  };

  const loadPalettes = () => {
    const pyodide = getPyodide();
    if (!pyodide) return null;
    const paletteJson = pyodide.runPython('get_palettes_json()');
    return JSON.parse(paletteJson);
  };

  const initPyodide = async () => {
    try {
      const pyodide = await loadPyodide();
      setPyodide(pyodide);
      loadingStatus.value = 'Installing dependencies...';
      await pyodide.loadPackage('micropip');
      const micropip = pyodide.pyimport('micropip');
      await micropip.install(['biopython', 'svgwrite', 'pandas', 'fonttools', 'bcbio-gff']);
      loadingStatus.value = 'Installing gbdraw...';
      await micropip.install(GBDRAW_WHEEL_NAME);
      await pyodide.runPythonAsync(PYTHON_HELPERS);

      const allPalettes = loadPalettes();
      if (allPalettes) {
        paletteNames.value = Object.keys(allPalettes).filter((k) => k !== 'title').sort();
        currentColors.value = allPalettes.default || {};
      }
      pyodideReady.value = true;
    } catch (e) {
      loadingStatus.value = 'Startup Error: ' + e.message;
      console.error(e);
    }
  };

  const writeFileToFs = async (fileObj, path) => {
    if (!fileObj) return false;
    const pyodide = getPyodide();
    if (!pyodide) return false;
    const buffer = await fileObj.arrayBuffer();
    pyodide.FS.writeFile(path, new Uint8Array(buffer));
    return true;
  };

  return {
    getPyodide,
    initPyodide,
    loadPalettes,
    writeFileToFs
  };
};
