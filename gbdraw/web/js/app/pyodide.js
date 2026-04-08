import {
  GBDRAW_WHEEL_NAME,
  GBDRAW_WHEEL_CACHE_BUST,
  PYODIDE_INDEX_URL,
  PYODIDE_LOCAL_WHEELS
} from '../config.js';
import { PYTHON_HELPERS } from './python-helpers.js';

export const createPyodideManager = ({ state }) => {
  const {
    pyodideReady,
    loadingStatus,
    paletteDefinitions,
    paletteNames,
    selectedPalette,
    currentColors,
    appliedPaletteName,
    appliedPaletteColors,
    pendingPaletteName,
    pendingPaletteColors
  } = state;
  const pyodideRef = { current: null };

  const resolveAssetUrl = (path) => new URL(path, window.location.href).toString();

  const ensureLocalAsset = async (url, label) => {
    const response = await fetch(url, { method: 'HEAD', cache: 'no-store' });
    if (!response.ok) {
      throw new Error(`Missing packaged asset: ${label} (${response.status}) at ${url}`);
    }
    return url;
  };

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
    const wheelBaseUrl = new URL(GBDRAW_WHEEL_NAME, window.location.href);
    if (GBDRAW_WHEEL_CACHE_BUST) {
      wheelBaseUrl.searchParams.set('v', GBDRAW_WHEEL_CACHE_BUST);
    }
    const wheelUrl = wheelBaseUrl.toString();
    return ensureLocalAsset(wheelUrl, `gbdraw browser wheel ${GBDRAW_WHEEL_NAME}`);
  };

  const ensureLocalDependencyWheels = async () => {
    const wheelUrls = PYODIDE_LOCAL_WHEELS.map((path) => resolveAssetUrl(path));
    await Promise.all(
      wheelUrls.map((url, index) =>
        ensureLocalAsset(url, `Pyodide dependency wheel ${PYODIDE_LOCAL_WHEELS[index]}`)
      )
    );
    return wheelUrls;
  };

  const initPyodide = async () => {
    try {
      const pyodideIndexUrl = resolveAssetUrl(PYODIDE_INDEX_URL);
      loadingStatus.value = 'Loading local Pyodide runtime...';
      const pyodide = await loadPyodide({
        indexURL: pyodideIndexUrl,
        packageBaseUrl: pyodideIndexUrl
      });
      setPyodide(pyodide);
      loadingStatus.value = 'Loading local micropip...';
      await pyodide.loadPackage('micropip');
      const micropip = pyodide.pyimport('micropip');
      loadingStatus.value = 'Installing local Python dependencies...';
      const localWheelUrls = await ensureLocalDependencyWheels();
      await micropip.install(localWheelUrls);
      loadingStatus.value = 'Installing gbdraw...';
      const wheelUrl = await ensureWheelAvailable();
      await micropip.install(wheelUrl);
      await pyodide.runPythonAsync(PYTHON_HELPERS);

      const allPalettes = loadPalettes();
      if (allPalettes) {
        const normalizedPalettes = { ...allPalettes };
        delete normalizedPalettes.title;
        paletteDefinitions.value = normalizedPalettes;
        paletteNames.value = Object.keys(allPalettes).filter((k) => k !== 'title').sort();
        const defaultColors = { ...(allPalettes.default || {}) };
        selectedPalette.value = 'default';
        currentColors.value = defaultColors;
        appliedPaletteName.value = 'default';
        appliedPaletteColors.value = { ...defaultColors };
        pendingPaletteName.value = '';
        pendingPaletteColors.value = {};
      }
      pyodideReady.value = true;
    } catch (e) {
      const message = e instanceof Error ? e.message : String(e);
      loadingStatus.value = 'Startup Error: ' + message;
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
