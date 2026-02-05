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

  const ensureWheelAvailable = async () => {
    const wheelUrl = new URL(GBDRAW_WHEEL_NAME, window.location.href).toString();
    const response = await fetch(wheelUrl, { cache: 'no-store' });
    if (!response.ok) {
      throw new Error(
        `Wheel not found (${response.status}). Expected ${GBDRAW_WHEEL_NAME} at the site root. ` +
        'Build the wheel and copy it into gbdraw/web, or update gbdraw/web/js/config.js.'
      );
    }
    const buffer = await response.arrayBuffer();
    const bytes = new Uint8Array(buffer.slice(0, 2));
    if (bytes.length < 2 || bytes[0] !== 0x50 || bytes[1] !== 0x4b) {
      throw new Error(
        `Wheel file is not a valid zip: ${GBDRAW_WHEEL_NAME}. ` +
        'Ensure the file is a .whl built from this version and is served correctly.'
      );
    }
    return wheelUrl;
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
      const wheelUrl = await ensureWheelAvailable();
      await micropip.install(wheelUrl);
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
