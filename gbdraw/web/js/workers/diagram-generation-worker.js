import { PYTHON_HELPERS } from '../app/python-helpers.js';

let runtimePromise = null;
let runtime = null;

const serializeError = (error) => ({
  name: error?.name ? String(error.name) : 'Error',
  message: error?.message ? String(error.message) : String(error || 'Unknown diagram generation error'),
  stack: error?.stack ? String(error.stack) : ''
});

const ensureLocalAsset = async (url, label) => {
  const response = await fetch(url, { method: 'HEAD', cache: 'no-store' });
  if (!response.ok) {
    throw new Error(`Missing packaged asset: ${label} (${response.status}) at ${url}`);
  }
  return url;
};

const initializeRuntime = async ({
  pyodideIndexUrl,
  pyodideModuleUrl,
  localWheelUrls = [],
  gbdrawWheelUrl
} = {}) => {
  if (runtime) return runtime;
  if (runtimePromise) return runtimePromise;
  if (!pyodideIndexUrl || !pyodideModuleUrl) {
    throw new Error('Diagram generation worker requires Pyodide asset URLs.');
  }
  if (!gbdrawWheelUrl) {
    throw new Error('Diagram generation worker requires the gbdraw wheel URL.');
  }

  runtimePromise = (async () => {
    const { loadPyodide } = await import(pyodideModuleUrl);
    const pyodide = await loadPyodide({
      indexURL: pyodideIndexUrl,
      packageBaseUrl: pyodideIndexUrl
    });
    await pyodide.loadPackage('micropip');
    const micropip = pyodide.pyimport('micropip');

    await Promise.all(
      localWheelUrls.map((url, index) =>
        ensureLocalAsset(url, `Pyodide dependency wheel #${index + 1}`)
      )
    );
    await ensureLocalAsset(gbdrawWheelUrl, 'gbdraw browser wheel');
    await micropip.install(localWheelUrls);
    await micropip.install(gbdrawWheelUrl);
    await pyodide.runPythonAsync(PYTHON_HELPERS);

    runtime = { pyodide };
    return runtime;
  })();

  try {
    return await runtimePromise;
  } catch (error) {
    runtimePromise = null;
    throw error;
  }
};

const ensureParentDirectory = (pyodide, path) => {
  const parts = String(path || '').split('/').filter(Boolean);
  parts.pop();
  if (!parts.length) return;

  let current = '';
  parts.forEach((part) => {
    current += `/${part}`;
    try {
      pyodide.FS.mkdir(current);
    } catch (error) {
      if (!String(error?.message || '').includes('File exists')) {
        throw error;
      }
    }
  });
};

const writeGenerationFiles = (pyodide, files = []) => {
  files.forEach((file) => {
    const path = String(file?.path || '').trim();
    if (!path) return;
    const bytes = file.bytes instanceof ArrayBuffer
      ? new Uint8Array(file.bytes)
      : new Uint8Array(file.bytes || []);
    ensureParentDirectory(pyodide, path);
    pyodide.FS.writeFile(path, bytes);
  });
};

const runGeneration = async ({
  mode,
  args = [],
  files = [],
  virtualBlastFiles = []
} = {}) => {
  if (!runtime?.pyodide) {
    throw new Error('Diagram generation worker has not been initialized.');
  }
  const { pyodide } = runtime;
  writeGenerationFiles(pyodide, files);

  const runWrapper = pyodide.globals.get('run_gbdraw_wrapper');
  const pyArgs = pyodide.toPy(args.map((arg) => String(arg)));
  try {
    const resultJson = runWrapper(
      String(mode || ''),
      pyArgs,
      virtualBlastFiles.length ? JSON.stringify(virtualBlastFiles) : null
    );
    return JSON.parse(String(resultJson || 'null'));
  } finally {
    pyArgs.destroy?.();
    runWrapper.destroy?.();
  }
};

const runFeatureExtraction = async ({
  path,
  files = [],
  regionSpec = null,
  recordSelector = null,
  reverseFlag = false,
  selectedFeatures = null
} = {}) => {
  if (!runtime?.pyodide) {
    throw new Error('Diagram generation worker has not been initialized.');
  }
  const normalizedPath = String(path || '').trim();
  if (!normalizedPath) {
    throw new Error('Feature extraction requires a GenBank input path.');
  }

  const { pyodide } = runtime;
  writeGenerationFiles(pyodide, files);

  const extractFeatures = pyodide.globals.get('extract_features_from_genbank');
  try {
    const resultJson = extractFeatures(
      normalizedPath,
      regionSpec || null,
      recordSelector || null,
      reverseFlag ? '1' : '0',
      Array.isArray(selectedFeatures) && selectedFeatures.length ? JSON.stringify(selectedFeatures) : null
    );
    return JSON.parse(String(resultJson || 'null'));
  } finally {
    extractFeatures.destroy?.();
  }
};

self.onmessage = async (event) => {
  const data = event.data || {};
  const { id, requestId, type } = data;
  try {
    if (type === 'init') {
      await initializeRuntime(data);
      self.postMessage({ id, type: 'init', ok: true });
      return;
    }
    if (type === 'ping') {
      self.postMessage({ id, type: 'ping', ok: true });
      return;
    }
    if (type === 'feature-extraction') {
      const result = await runFeatureExtraction(data.payload || {});
      self.postMessage({ requestId, type: 'feature-extraction', ok: true, result });
      return;
    }
    if (type !== 'run') {
      throw new Error(`Unsupported diagram generation worker message type '${type || '(blank)'}'.`);
    }

    const results = await runGeneration(data.payload || {});
    self.postMessage({ requestId, type: 'run', ok: true, results });
  } catch (error) {
    self.postMessage({
      id,
      requestId,
      type: type || 'run',
      ok: false,
      error: serializeError(error)
    });
  }
};
