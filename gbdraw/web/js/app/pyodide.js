import {
  ensureLocalAsset,
  resolveGbdrawWheelUrl,
  resolveLocalDependencyWheelUrls,
  resolvePyodideIndexUrl
} from '../services/pyodide-assets.js';
import { COMPARISON_COLOR_KEYS } from './color-utils.js';
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
    pendingPaletteColors,
    normalizePaletteColors,
    normalizePaletteDefinitions
  } = state;
  const pyodideRef = { current: null };
  let initPromise = null;

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

  const hasColorEntries = (colors) =>
    Boolean(colors) && typeof colors === 'object' && !Array.isArray(colors) && Object.keys(colors).length > 0;
  const comparisonColorKeys = new Set(COMPARISON_COLOR_KEYS);
  const hasPaletteColorEntries = (colors) =>
    hasColorEntries(colors) && Object.keys(colors).some((key) => !comparisonColorKeys.has(key));

  const ensureWheelAvailable = async () => {
    const wheelUrl = resolveGbdrawWheelUrl();
    return ensureLocalAsset(wheelUrl, 'gbdraw browser wheel');
  };

  const ensureLocalDependencyWheels = async () => {
    const wheelUrls = resolveLocalDependencyWheelUrls();
    await Promise.all(
      wheelUrls.map((url, index) =>
        ensureLocalAsset(url, `Pyodide dependency wheel #${index + 1}`)
      )
    );
    return wheelUrls;
  };

  const initPyodide = async () => {
    if (pyodideReady.value) return getPyodide();
    if (initPromise) return initPromise;

    initPromise = (async () => {
      try {
        const pyodideIndexUrl = resolvePyodideIndexUrl();
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
          const normalizedPalettes = normalizePaletteDefinitions(allPalettes);
          paletteDefinitions.value = normalizedPalettes;
          paletteNames.value = Object.keys(normalizedPalettes).sort();
          const requestedPalette = String(selectedPalette.value || 'default').trim() || 'default';
          const resolvedPalette = normalizedPalettes[requestedPalette] ? requestedPalette : 'default';
          const resolvedColors = normalizePaletteColors(
            normalizedPalettes[resolvedPalette] || normalizedPalettes.default || {}
          );
          const currentHasPaletteColors = hasPaletteColorEntries(currentColors.value);
          const appliedHasColors = hasPaletteColorEntries(appliedPaletteColors.value);
          const pendingHasColors = hasPaletteColorEntries(pendingPaletteColors.value);

          selectedPalette.value = resolvedPalette;
          if (!currentHasPaletteColors) {
            currentColors.value = resolvedColors;
          }
          if (!appliedHasColors) {
            appliedPaletteName.value = resolvedPalette;
            appliedPaletteColors.value = {
              ...(currentHasPaletteColors ? currentColors.value : resolvedColors)
            };
          }
          if (String(pendingPaletteName.value || '').trim() && !pendingHasColors) {
            pendingPaletteColors.value = {
              ...(currentHasPaletteColors ? currentColors.value : resolvedColors)
            };
          }
        }
        pyodideReady.value = true;
        return pyodide;
      } catch (e) {
        const message = e instanceof Error ? e.message : String(e);
        loadingStatus.value = 'Startup Error: ' + message;
        console.error(e);
        return null;
      } finally {
        initPromise = null;
      }
    })();

    return initPromise;
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
