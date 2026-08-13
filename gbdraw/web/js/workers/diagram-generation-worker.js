import { PYTHON_HELPERS } from '../app/python-helpers.js';
import { DIAGRAM_HELPER_OPERATIONS } from '../services/diagram-worker-protocol.js';
import { normalizeUserFacingError } from '../services/error-normalization.js';

let runtimePromise = null;
let runtime = null;
let operationQueue = Promise.resolve();

export const serializeError = (error) => {
  const normalized = normalizeUserFacingError(error);
  return {
    name: error?.name ? String(error.name) : 'Error',
    message: normalized?.summary || 'Unknown diagram generation error',
    details: Array.isArray(normalized?.details) ? normalized.details : [],
    notes: Array.isArray(error?.notes) ? error.notes.map(String) : [],
    stack: error?.stack ? String(error.stack) : ''
  };
};

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

const HELPER_FILE_NAMES = Object.freeze({
  source: 'source.bin',
  gff: 'source.gff',
  fasta: 'source.fasta',
  visibility: 'feature-visibility.tsv'
});

const requirePayloadObject = (payload, operation) => {
  if (!payload || typeof payload !== 'object' || Array.isArray(payload)) {
    throw new TypeError(`Diagram helper '${operation}' requires an object payload.`);
  }
  return payload;
};

const assertAllowedPayloadKeys = (payload, operation, allowedKeys) => {
  const allowed = new Set(allowedKeys);
  Object.keys(payload).forEach((key) => {
    if (!allowed.has(key)) {
      throw new TypeError(`Diagram helper '${operation}' does not accept '${key}'.`);
    }
  });
};

const removeWorkspace = (pyodide, path) => {
  if (!pyodide.FS.analyzePath(path).exists) return;
  pyodide.FS.readdir(path).forEach((entry) => {
    if (entry === '.' || entry === '..') return;
    const child = `${path}/${entry}`;
    const stat = pyodide.FS.stat(child);
    if (pyodide.FS.isDir(stat.mode)) removeWorkspace(pyodide, child);
    else pyodide.FS.unlink(child);
  });
  pyodide.FS.rmdir(path);
};

const withRequestWorkspace = async (pyodide, requestId, kind, callback) => {
  const workspace = `/gbdraw-web-${kind}-${Number(requestId) || 0}`;
  pyodide.FS.mkdir(workspace);
  let primaryError = null;
  try {
    return await callback(workspace);
  } catch (error) {
    primaryError = error;
    throw error;
  } finally {
    try {
      removeWorkspace(pyodide, workspace);
    } catch (cleanupError) {
      if (primaryError) {
        attachCleanupDiagnostic(
          primaryError,
          cleanupError,
          'Temporary helper workspace cleanup also failed'
        );
      } else {
        throw cleanupError;
      }
    }
  }
};

const stageHelperFiles = (pyodide, workspace, files, allowedRoles) => {
  const allowed = new Set(allowedRoles);
  const paths = new Map();
  const entries = Array.isArray(files) ? files : [];
  entries.forEach((file) => {
    const role = String(file?.role || '').trim();
    if (!allowed.has(role)) {
      throw new TypeError(`Unexpected diagram helper file role '${role || '(blank)'}'.`);
    }
    if (paths.has(role)) {
      throw new TypeError(`Duplicate diagram helper file role '${role}'.`);
    }
    if (!(file?.bytes instanceof ArrayBuffer)) {
      throw new TypeError(`Diagram helper file '${role}' requires an ArrayBuffer.`);
    }
    const path = `${workspace}/${HELPER_FILE_NAMES[role]}`;
    pyodide.FS.writeFile(path, new Uint8Array(file.bytes));
    paths.set(role, path);
  });
  return paths;
};

const requireHelperFile = (paths, role, operation) => {
  const path = paths.get(role);
  if (!path) throw new TypeError(`Diagram helper '${operation}' requires the '${role}' file.`);
  return path;
};

const callJsonHelper = (pyodide, helperName, args) => {
  const helper = pyodide.globals.get(helperName);
  try {
    if (typeof helper !== 'function') {
      throw new Error(`Packaged Python helper '${helperName}' is unavailable.`);
    }
    const rawResult = helper(...args);
    return JSON.parse(String(rawResult || 'null'));
  } finally {
    helper?.destroy?.();
  }
};

const jsonArgument = (value, fallback) => JSON.stringify(value ?? fallback);

const HELPER_OPERATION_SPECS = Object.freeze({
  [DIAGRAM_HELPER_OPERATIONS.EXTRACT_FIRST_FASTA]: {
    keys: ['files', 'format', 'regionSpec', 'recordSelector', 'reverseFlag'],
    fileRoles: ['source'],
    run: (pyodide, payload, paths, operation) => callJsonHelper(
      pyodide,
      'extract_first_fasta',
      [
        requireHelperFile(paths, 'source', operation),
        String(payload.format || 'genbank').trim().toLowerCase(),
        payload.regionSpec ?? null,
        payload.recordSelector ?? null,
        payload.reverseFlag ? '1' : '0'
      ]
    )
  },
  [DIAGRAM_HELPER_OPERATIONS.EXTRACT_CDS_PROTEIN_FASTA]: {
    keys: [
      'files',
      'format',
      'regionSpec',
      'recordSelector',
      'reverseFlag',
      'recordIndex',
      'recordInstanceKey'
    ],
    fileRoles: ['source', 'fasta', 'visibility'],
    run: (pyodide, payload, paths, operation) => {
      const format = String(payload.format || 'genbank').trim().toLowerCase();
      const fastaPath = paths.get('fasta') || null;
      if (format === 'gff' && !fastaPath) {
        throw new TypeError(`Diagram helper '${operation}' requires the 'fasta' file for GFF3.`);
      }
      return callJsonHelper(pyodide, 'extract_cds_protein_fasta', [
        requireHelperFile(paths, 'source', operation),
        format,
        fastaPath,
        payload.regionSpec ?? null,
        payload.recordSelector ?? null,
        payload.reverseFlag ? '1' : '0',
        payload.recordIndex ?? null,
        payload.recordInstanceKey ?? null,
        paths.get('visibility') || null
      ]);
    }
  },
  [DIAGRAM_HELPER_OPERATIONS.BUILD_PROTEIN_LOSAT_CACHE_KEY]: {
    keys: [
      'identityManifest',
      'queryRecordInstanceKey',
      'subjectRecordInstanceKey',
      'expectedOptions'
    ],
    fileRoles: [],
    run: (pyodide, payload) => callJsonHelper(
      pyodide,
      'build_protein_losat_cache_key_json',
      [
        jsonArgument(payload.identityManifest, {}),
        String(payload.queryRecordInstanceKey || ''),
        String(payload.subjectRecordInstanceKey || ''),
        jsonArgument(payload.expectedOptions, {})
      ]
    )
  },
  [DIAGRAM_HELPER_OPERATIONS.PROMOTE_LEGACY_LOSATP_CACHE]: {
    keys: [
      'candidates',
      'queryFasta',
      'subjectFasta',
      'queryProteinMap',
      'subjectProteinMap',
      'identityManifest',
      'expectedOptions'
    ],
    fileRoles: [],
    run: (pyodide, payload) => callJsonHelper(
      pyodide,
      'promote_legacy_losatp_cache_candidates',
      [
        jsonArgument(payload.candidates, []),
        String(payload.queryFasta || ''),
        String(payload.subjectFasta || ''),
        jsonArgument(payload.queryProteinMap, {}),
        jsonArgument(payload.subjectProteinMap, {}),
        jsonArgument(payload.identityManifest, {}),
        jsonArgument(payload.expectedOptions, {})
      ]
    )
  },
  [DIAGRAM_HELPER_OPERATIONS.RESOLVE_LEGACY_PROTEIN_REFERENCES]: {
    keys: ['proteinRecords', 'identityManifest', 'referenceIds'],
    fileRoles: [],
    run: (pyodide, payload) => callJsonHelper(
      pyodide,
      'resolve_legacy_protein_reference_map_json',
      [
        jsonArgument(payload.proteinRecords, []),
        jsonArgument(payload.identityManifest, {}),
        jsonArgument(payload.referenceIds, [])
      ]
    )
  },
  [DIAGRAM_HELPER_OPERATIONS.CONVERT_LOSATP_PAIRS_TO_GENOMIC_PAYLOAD]: {
    keys: [
      'pairs',
      'mode',
      'maxHits',
      'bitscore',
      'evalue',
      'identity',
      'alignmentLength',
      'collinearMinAnchors',
      'collinearMaxUnitGap',
      'collinearUnitMode',
      'collinearColorMode',
      'collinearAnchorMode',
      'collinearMaxDiagonalDrift',
      'collinearMaxConflictsInMergeGap',
      'collinearMaxParalogLinksPerOrthogroup',
      'collinearSearchScope',
      'orthogroupMembershipMode',
      'orthogroupMemberMaxHits'
    ],
    fileRoles: [],
    run: (pyodide, payload) => callJsonHelper(
      pyodide,
      'convert_losatp_blastp_pairs_to_genomic_payload',
      [
        jsonArgument(payload.pairs, {}),
        payload.mode ?? 'pairwise',
        payload.maxHits ?? 5,
        payload.bitscore ?? 50,
        payload.evalue ?? '1e-2',
        payload.identity ?? 0,
        payload.alignmentLength ?? 0,
        payload.collinearMinAnchors ?? 1,
        payload.collinearMaxUnitGap ?? 0,
        payload.collinearUnitMode ?? 'auto',
        payload.collinearColorMode ?? 'orientation',
        payload.collinearAnchorMode ?? 'rbh',
        payload.collinearMaxDiagonalDrift ?? 0,
        payload.collinearMaxConflictsInMergeGap ?? 1,
        payload.collinearMaxParalogLinksPerOrthogroup ?? 2,
        payload.collinearSearchScope ?? 'adjacent',
        payload.orthogroupMembershipMode ?? 'anchor_core_v1',
        payload.orthogroupMemberMaxHits ?? 5
      ]
    )
  },
  [DIAGRAM_HELPER_OPERATIONS.CONVERT_LOSAT_NUCLEOTIDE_TO_DISPLAY_TSV]: {
    keys: ['blastText', 'queryViewTransform', 'subjectViewTransform'],
    fileRoles: [],
    run: (pyodide, payload) => callJsonHelper(
      pyodide,
      'convert_losat_nucleotide_to_display_tsv',
      [
        String(payload.blastText || ''),
        jsonArgument(payload.queryViewTransform, {}),
        jsonArgument(payload.subjectViewTransform, {})
      ]
    )
  },
  [DIAGRAM_HELPER_OPERATIONS.HYDRATE_PROTEIN_LOSAT_TSV]: {
    keys: ['entry', 'identityManifest'],
    fileRoles: [],
    run: (pyodide, payload) => callJsonHelper(
      pyodide,
      'hydrate_protein_losat_tsv_json',
      [jsonArgument(payload.entry, {}), jsonArgument(payload.identityManifest, {})]
    )
  },
  [DIAGRAM_HELPER_OPERATIONS.REGENERATE_DEFINITION_SVGS]: {
    keys: [
      'files',
      'species',
      'strain',
      'plotTitle',
      'definitionFontSize',
      'plotTitleFontSize',
      'plotTitlePosition',
      'multiRecordCanvas',
      'keepFullDefinitionWithPlotTitle'
    ],
    fileRoles: ['source'],
    run: (pyodide, payload, paths, operation) => callJsonHelper(
      pyodide,
      'regenerate_definition_svgs',
      [
        requireHelperFile(paths, 'source', operation),
        payload.species ?? null,
        payload.strain ?? null,
        payload.plotTitle ?? null,
        payload.definitionFontSize ?? null,
        payload.plotTitleFontSize ?? null,
        payload.plotTitlePosition ?? 'none',
        Boolean(payload.multiRecordCanvas),
        Boolean(payload.keepFullDefinitionWithPlotTitle)
      ]
    )
  },
  [DIAGRAM_HELPER_OPERATIONS.LIST_SEQUENCE_RECORDS]: {
    keys: ['files', 'format'],
    fileRoles: ['source'],
    run: (pyodide, payload, paths, operation) => callJsonHelper(
      pyodide,
      'list_sequence_records',
      [
        requireHelperFile(paths, 'source', operation),
        String(payload.format || 'genbank').trim().toLowerCase()
      ]
    )
  },
  [DIAGRAM_HELPER_OPERATIONS.LIST_GFF_FASTA_RECORDS]: {
    keys: ['files'],
    fileRoles: ['gff', 'fasta'],
    run: (pyodide, _payload, paths, operation) => callJsonHelper(
      pyodide,
      'list_gff_fasta_records',
      [
        requireHelperFile(paths, 'gff', operation),
        requireHelperFile(paths, 'fasta', operation)
      ]
    )
  },
  [DIAGRAM_HELPER_OPERATIONS.MEASURE_LEGEND_TEXT]: {
    keys: ['caption', 'fontFamily', 'fontSize'],
    fileRoles: [],
    run: (pyodide, payload) => callJsonHelper(
      pyodide,
      'measure_legend_text_json',
      [String(payload.caption || ''), String(payload.fontFamily || 'Arial'), payload.fontSize ?? 14]
    )
  },
  [DIAGRAM_HELPER_OPERATIONS.GENERATE_LEGEND_ENTRY_SVG]: {
    keys: [
      'caption',
      'color',
      'yOffset',
      'rectSize',
      'fontSize',
      'fontFamily',
      'xOffset',
      'strokeColor',
      'strokeWidth'
    ],
    fileRoles: [],
    run: (pyodide, payload) => callJsonHelper(
      pyodide,
      'generate_legend_entry_svg',
      [
        String(payload.caption || ''),
        String(payload.color || ''),
        payload.yOffset ?? 0,
        payload.rectSize ?? 14,
        payload.fontSize ?? 14,
        String(payload.fontFamily || 'Arial'),
        payload.xOffset ?? 0,
        String(payload.strokeColor || 'black'),
        payload.strokeWidth ?? 0.5
      ]
    )
  }
});

const runHelperOperation = async ({ operation, payload, requestId } = {}) => {
  if (!runtime?.pyodide) {
    throw new Error('Diagram generation worker has not been initialized.');
  }
  const normalizedOperation = String(operation || '').trim();
  const spec = Object.prototype.hasOwnProperty.call(
    HELPER_OPERATION_SPECS,
    normalizedOperation
  )
    ? HELPER_OPERATION_SPECS[normalizedOperation]
    : null;
  if (!spec) {
    throw new TypeError(`Unsupported diagram helper operation '${normalizedOperation || '(blank)'}'.`);
  }
  const normalizedPayload = requirePayloadObject(payload, normalizedOperation);
  assertAllowedPayloadKeys(normalizedPayload, normalizedOperation, [
    ...spec.keys
  ]);
  return withRequestWorkspace(
    runtime.pyodide,
    requestId,
    'helper',
    async (workspace) => {
      const paths = stageHelperFiles(
        runtime.pyodide,
        workspace,
        normalizedPayload.files,
        spec.fileRoles
      );
      return spec.run(runtime.pyodide, normalizedPayload, paths, normalizedOperation);
    }
  );
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
  includeBiologicalFeatures = false,
  requestId = 0
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
  const normalizedFiles = Array.isArray(files) ? files : [];
  const sourceFile = normalizedFiles.find(
    (file) => String(file?.path || '').trim() === normalizedPath
  );
  const fastaFile = isGff
    ? normalizedFiles.find((file) => String(file?.path || '').trim() === normalizedFastaPath)
    : null;
  const normalizedVisibilityPath = String(featureVisibilityTablePath || '').trim();
  const visibilityFile = normalizedVisibilityPath
    ? normalizedFiles.find((file) => String(file?.path || '').trim() === normalizedVisibilityPath)
    : null;
  const workspaceFiles = [
    { role: 'source', bytes: sourceFile?.bytes },
    ...(fastaFile ? [{ role: 'fasta', bytes: fastaFile.bytes }] : []),
    ...(visibilityFile ? [{ role: 'visibility', bytes: visibilityFile.bytes }] : [])
  ];

  return withRequestWorkspace(pyodide, requestId, 'feature', async (workspace) => {
    const paths = stageHelperFiles(
      pyodide,
      workspace,
      workspaceFiles,
      ['source', 'fasta', 'visibility']
    );
    const sourcePath = requireHelperFile(paths, 'source', 'feature-extraction');
    const stagedFastaPath = isGff
      ? requireHelperFile(paths, 'fasta', 'feature-extraction')
      : null;
    const selectedFeaturesJson = Array.isArray(selectedFeatures) && selectedFeatures.length
      ? JSON.stringify(selectedFeatures)
      : null;
    return callJsonHelper(
      pyodide,
      isGff ? 'extract_features_from_gff_fasta' : 'extract_features_from_genbank',
      isGff
        ? [
          sourcePath,
          stagedFastaPath,
          regionSpec || null,
          recordSelector || null,
          reverseFlag ? '1' : '0',
          selectedFeaturesJson,
          paths.get('visibility') || null,
          Boolean(includeBiologicalFeatures)
        ]
        : [
          sourcePath,
          regionSpec || null,
          recordSelector || null,
          reverseFlag ? '1' : '0',
          selectedFeaturesJson,
          paths.get('visibility') || null,
          Boolean(includeBiologicalFeatures)
        ]
    );
  });
};

const handleWorkerMessage = async (data) => {
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
      const result = await runFeatureExtraction({
        ...(data.payload || {}),
        requestId
      });
      self.postMessage({ requestId, type: 'feature-extraction', ok: true, result });
      return;
    }
    if (type === 'helper') {
      const result = await runHelperOperation({
        operation: data.operation,
        payload: data.payload || {},
        requestId
      });
      self.postMessage({ requestId, type: 'helper', ok: true, result });
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

self.onmessage = (event) => {
  const data = event.data || {};
  const scheduled = operationQueue.then(
    () => handleWorkerMessage(data),
    () => handleWorkerMessage(data)
  );
  operationQueue = scheduled.catch(() => {});
};
