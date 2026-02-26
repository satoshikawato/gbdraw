import {
  GBDRAW_WHEEL_NAME,
  GBDRAW_WHEEL_CACHE_BUST,
  PYODIDE_INDEX_URL,
  PYODIDE_WHEELS_BASE_URL,
  PYODIDE_REQUIRED_WHEELS
} from '../config.js';
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

  const toAbsoluteAssetUrl = (pathOrUrl) => new URL(pathOrUrl, window.location.href).toString();

  const ensureZipAssetAvailable = async (label, pathOrUrl) => {
    const assetUrl = toAbsoluteAssetUrl(pathOrUrl);
    const response = await fetch(assetUrl, { cache: 'no-store' });
    if (!response.ok) {
      throw new Error(
        `${label} missing (${response.status}). Expected local asset: ${pathOrUrl}`
      );
    }
    const buffer = await response.arrayBuffer();
    const bytes = new Uint8Array(buffer.slice(0, 2));
    if (bytes.length < 2 || bytes[0] !== 0x50 || bytes[1] !== 0x4b) {
      throw new Error(
        `${label} is not a valid .whl/.zip archive: ${pathOrUrl}`
      );
    }
    return assetUrl;
  };

  const ensureWheelAvailable = async () => {
    const wheelBaseUrl = new URL(GBDRAW_WHEEL_NAME, window.location.href);
    if (GBDRAW_WHEEL_CACHE_BUST) {
      wheelBaseUrl.searchParams.set('v', GBDRAW_WHEEL_CACHE_BUST);
    }
    return ensureZipAssetAvailable('gbdraw wheel', wheelBaseUrl.toString());
  };

  const getRequiredWheelNames = () => {
    if (!Array.isArray(PYODIDE_REQUIRED_WHEELS)) {
      throw new Error('PYODIDE_REQUIRED_WHEELS must be an array in gbdraw/web/js/config.js');
    }
    const wheelNames = PYODIDE_REQUIRED_WHEELS
      .map((entry) => String(entry || '').trim())
      .filter(Boolean);
    if (wheelNames.length === 0) {
      throw new Error(
        'PYODIDE_REQUIRED_WHEELS is empty. List all local dependency wheels in gbdraw/web/js/config.js.'
      );
    }
    return wheelNames;
  };

  const ensureDependencyWheelsAvailable = async () => {
    const wheelNames = getRequiredWheelNames();
    const wheelBase = toAbsoluteAssetUrl(PYODIDE_WHEELS_BASE_URL);
    const checks = wheelNames.map(async (wheelName) => {
      const wheelUrl = new URL(wheelName, wheelBase).toString();
      return ensureZipAssetAvailable(`Pyodide dependency wheel (${wheelName})`, wheelUrl);
    });
    return Promise.all(checks);
  };

  const initPyodide = async () => {
    try {
      if (typeof loadPyodide !== 'function') {
        throw new Error(
          `Pyodide runtime loader is unavailable. Expected local asset: ${PYODIDE_INDEX_URL}pyodide.js`
        );
      }
      const pyodide = await loadPyodide({ indexURL: toAbsoluteAssetUrl(PYODIDE_INDEX_URL) });
      setPyodide(pyodide);
      loadingStatus.value = 'Installing local dependencies...';
      await pyodide.loadPackage('micropip');
      const micropip = pyodide.pyimport('micropip');
      const dependencyWheelUrls = await ensureDependencyWheelsAvailable();
      await micropip.install(dependencyWheelUrls);
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
