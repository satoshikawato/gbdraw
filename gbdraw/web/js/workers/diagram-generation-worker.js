import { PYTHON_HELPERS } from '../app/python-helpers.js';

let runtimePromise = null;
let runtime = null;

export const serializeError = (error) => ({
  name: error?.name ? String(error.name) : 'Error',
  message: error?.message ? String(error.message) : String(error || 'Unknown diagram generation error'),
  stack: error?.stack ? String(error.stack) : ''
});

const errorDiagnostic = (error) => {
  const name = error?.name ? String(error.name) : 'Error';
  const message = error?.message ? String(error.message) : String(error || 'Unknown error');
  return `${name}: ${message}`;
};

const attachCleanupDiagnostic = (
  primary,
  cleanupError,
  label
) => {
  if (!primary || typeof primary !== 'object' || !cleanupError) return primary;
  const note = `${label}: ${errorDiagnostic(cleanupError)}`;
  const notes = Array.isArray(primary.notes) ? [...primary.notes] : [];
  if (!notes.includes(note)) notes.push(note);
  primary.notes = notes;
  for (const field of ['traceback', 'stack']) {
    if (typeof primary[field] !== 'string' || primary[field].includes(note)) continue;
    primary[field] = `${primary[field].trimEnd()}\n${note}`;
  }
  return primary;
};

export const resolveGenerationCleanupOutcome = ({
  result = null,
  primaryError = null,
  destroyError = null,
  workspaceError = null
} = {}) => {
  const pythonError = (
    result?.error &&
    typeof result.error === 'object' &&
    !Array.isArray(result.error)
  )
    ? result.error
    : null;
  const primary = primaryError || pythonError || destroyError || workspaceError;
  if (primary) {
    if (primary !== destroyError && destroyError) {
      attachCleanupDiagnostic(
        primary,
        destroyError,
        'Temporary render handle cleanup also failed'
      );
    }
    if (primary !== workspaceError && workspaceError) {
      attachCleanupDiagnostic(
        primary,
        workspaceError,
        'Temporary render workspace cleanup also failed'
      );
    }
  }
  if (primaryError) throw primaryError;
  if (pythonError) return result;
  if (destroyError) throw destroyError;
  if (workspaceError) throw workspaceError;
  return result;
};

const ensureLocalAsset = async (url, label) => {
  const response = await fetch(url, { method: 'HEAD', cache: 'no-store' });
  if (!response.ok) {
    throw new Error(`Missing packaged asset: ${label} (${response.status}) at ${url}`);
  }
  return url;
};

const readRuntimeCapabilities = (pyodide) => {
  const raw = pyodide.runPython(`
import json as _json
from gbdraw.api import get_web_runtime_capabilities as _get_web_runtime_capabilities
_json.dumps(_get_web_runtime_capabilities(), sort_keys=True)
  `);
  return JSON.parse(String(raw));
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

    const capabilities = readRuntimeCapabilities(pyodide);
    runtime = { pyodide, capabilities };
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
  request,
  resources,
  requestId
} = {}) => {
  if (!runtime?.pyodide) {
    throw new Error('Diagram generation worker has not been initialized.');
  }
  if (!request || typeof request !== 'object' || Array.isArray(request)) {
    throw new Error('Diagram generation requires a canonical typed request.');
  }
  if (!resources || typeof resources !== 'object' || Array.isArray(resources)) {
    throw new Error('Diagram generation requires canonical request resources.');
  }
  const { pyodide } = runtime;
  const runWrapper = pyodide.globals.get('run_canonical_request_wrapper');
  const workspace = `/gbdraw-web-render-${Number(requestId) || 0}`;
  let result = null;
  let primaryError = null;
  try {
    const resultJson = runWrapper(
      JSON.stringify(request),
      JSON.stringify(resources),
      workspace
    );
    result = JSON.parse(String(resultJson || 'null'));
  } catch (error) {
    primaryError = error;
  }
  let destroyError = null;
  try {
    runWrapper.destroy?.();
  } catch (error) {
    destroyError = error;
  }
  let workspaceError = null;
  try {
    let remainingEntries = [];
    if (pyodide.FS.analyzePath(workspace).exists) {
      remainingEntries = pyodide.FS.readdir(workspace).filter(
        (entry) => entry !== '.' && entry !== '..'
      );
      if (remainingEntries.length === 0) {
        // Pyodide's MEMFS can retain the now-empty top-level directory after
        // Python shutil.rmtree has removed its contents.
        pyodide.FS.rmdir(workspace);
      }
    }
    if (pyodide.FS.analyzePath(workspace).exists) {
      workspaceError = new Error(
        `Diagram render workspace cleanup invariant failed: ${workspace}` +
        (remainingEntries.length > 0
          ? ` (${remainingEntries.join(', ')})`
          : '')
      );
    }
  } catch (error) {
    workspaceError = error;
  }
  return resolveGenerationCleanupOutcome({
    result,
    primaryError,
    destroyError,
    workspaceError
  });
};

const runFeatureExtraction = async ({
  path,
  format = 'genbank',
  fastaPath = null,
  files = [],
  regionSpec = null,
  recordSelector = null,
  reverseFlag = false,
  selectedFeatures = null,
  featureVisibilityTablePath = null,
  includeBiologicalFeatures = false
} = {}) => {
  if (!runtime?.pyodide) {
    throw new Error('Diagram generation worker has not been initialized.');
  }
  const normalizedPath = String(path || '').trim();
  if (!normalizedPath) {
    throw new Error('Feature extraction requires an input path.');
  }
  const normalizedFormat = String(format || 'genbank').trim().toLowerCase();
  const isGff = normalizedFormat === 'gff';
  const normalizedFastaPath = String(fastaPath || '').trim();
  if (isGff && !normalizedFastaPath) {
    throw new Error('GFF3 feature extraction requires a FASTA input path.');
  }

  const { pyodide } = runtime;
  writeGenerationFiles(pyodide, files);

  const extractFeatures = pyodide.globals.get(
    isGff ? 'extract_features_from_gff_fasta' : 'extract_features_from_genbank'
  );
  try {
    const selectedFeaturesJson = Array.isArray(selectedFeatures) && selectedFeatures.length
      ? JSON.stringify(selectedFeatures)
      : null;
    const resultJson = isGff
      ? extractFeatures(
          normalizedPath,
          normalizedFastaPath,
          regionSpec || null,
          recordSelector || null,
          reverseFlag ? '1' : '0',
          selectedFeaturesJson,
          featureVisibilityTablePath || null,
          Boolean(includeBiologicalFeatures)
        )
      : extractFeatures(
          normalizedPath,
          regionSpec || null,
          recordSelector || null,
          reverseFlag ? '1' : '0',
          selectedFeaturesJson,
          featureVisibilityTablePath || null,
          Boolean(includeBiologicalFeatures)
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
      const initialized = await initializeRuntime(data);
      self.postMessage({
        id,
        type: 'init',
        ok: true,
        capabilities: initialized.capabilities
      });
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

    const results = await runGeneration({
      ...(data.payload || {}),
      requestId
    });
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
