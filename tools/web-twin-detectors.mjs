// @ts-check
// Twin detector (boundary audit M2, plan row P1-D): logic that Python and the
// Web each implement. A twin is allowed only with a row in tools/twins.json
// that names its owners, the vector that pins the two to each other, and the
// condition under which one side is deleted. The registry only shrinks: a new
// row is an authority change, made before its twin exists with `pending`
// (gbdraw/web/CLAUDE.md R9); the change that adds the twin removes `pending`.
//
// `twin.unregistered.v1` observes two kinds of twin in source text; nothing is
// executed:
// - vector: a JavaScript test and a Python test under tests/ both name the same
//   tests/fixtures/**/*.json file;
// - function: a function whose name is on the watch list (migrate*, split*,
//   normalize*, resolve*Color*, read*Table*) is defined both in gbdraw/**/*.py
//   (outside gbdraw/web/) and in gbdraw/web/js/services/**/*.js. Names are
//   compared without case, underscores and leading underscores, so
//   `split_draft_into_modes` and `splitDraftIntoModes` are one name.
// Not seen: a fixture whose name a test builds at run time, a twin outside
// services/, and a twin under another name.
//
// Keep the detector ID and its helpers executable-identical while referenced;
// add a versioned detector beside it instead.
import { execFileSync } from 'node:child_process';
import { readFileSync } from 'node:fs';
import { join } from 'node:path';

import { maskJavaScript } from './web-change-source.mjs';

export const TWIN_DETECTOR_ID = 'twin.unregistered.v1';
export const TWIN_REGISTRY_PATH = 'tools/twins.json';

/**
 * @typedef {object} TwinSources
 * @property {ReadonlyMap<string, string>} tests JavaScript and Python files under tests/, by path.
 * @property {ReadonlyMap<string, string>} python Python files under gbdraw/ outside gbdraw/web/, by path.
 * @property {ReadonlyMap<string, string>} services JavaScript files under gbdraw/web/js/services/, by path.
 * @property {readonly string[]} fixtures The JSON files under tests/fixtures/, in any subfolder.
 */
/**
 * @typedef {object} DetectedTwin
 * @property {string} id `vector:<fixture path>` or `function:<JavaScript name>`.
 * @property {string[]} python The Python tests (vector) or `path:name` definitions (function).
 * @property {string[]} javascript The JavaScript tests (vector) or `path:name` definitions (function).
 */
/**
 * @typedef {object} TwinRow
 * @property {string} id
 * @property {string} concern
 * @property {string[]} pythonOwner
 * @property {string[]} jsOwner
 * @property {string | null} vector
 * @property {string} reason
 * @property {string} retireWhen
 * @property {string} source
 * @property {true} [pending] A `vector:` or `function:` row added before its twin.
 */
/**
 * @typedef {object} TwinRegistryCheck
 * @property {string[]} errors Malformed registry entries, and pending rows whose twin is detected.
 * @property {string[]} unregistered Detected twins without a row.
 * @property {string[]} stale Rows of a detected kind that the detector no longer reports.
 */

/** @param {string[]} values */
const unique = (values) => [...new Set(values)].sort();

/**
 * Blanks Python comments, and string contents unless `strings` is false;
 * offsets are kept.
 * @param {string} source
 * @param {{ strings?: boolean }} [options]
 */
const maskPython = (source, { strings = true } = {}) => source.replace(
  /#[^\n]*|("""|''')[\s\S]*?\1|"(?:\\.|[^"\\\n])*"|'(?:\\.|[^'\\\n])*'/g,
  (text) => (text.startsWith('#') || strings ? text.replace(/[^\n]/g, ' ') : text)
);

/** @param {string} name */
const comparableName = (name) => name.replace(/^_+/, '').replaceAll('_', '').toLowerCase();
const WATCHED_PREFIX = /^(?:migrate|split|normalize)(?=.)|^resolve.*color|^read.*table/;

/**
 * @param {Iterable<[string, string]>} entries
 * @param {(source: string) => Iterable<string>} definedNames
 * @returns {Map<string, string[]>} comparable name -> `path:name`
 */
const watchedDefinitions = (entries, definedNames) => {
  /** @type {Map<string, string[]>} */
  const byName = new Map();
  for (const [path, source] of entries) {
    for (const name of definedNames(source)) {
      const key = comparableName(name);
      if (WATCHED_PREFIX.test(key)) byName.set(key, [...(byName.get(key) || []), `${path}:${name}`]);
    }
  }
  return byName;
};

/** @param {string} source */
function* pythonDefinitions(source) {
  for (const match of maskPython(source).matchAll(/^[ \t]*(?:async[ \t]+)?def[ \t]+([A-Za-z_]\w*)/gm)) {
    yield /** @type {string} */ (match[1]);
  }
}

/** @param {string} source */
function* javaScriptDefinitions(source) {
  const pattern = /(?:^|[^\w$.])(?:function\*?[ \t]+([A-Za-z_$][\w$]*)|(?:const|let|var)[ \t]+([A-Za-z_$][\w$]*)[ \t]*=[ \t]*(?:async[ \t]*)?(?:function\b|\([^)]*\)[ \t]*=>|[A-Za-z_$][\w$]*[ \t]*=>))/gm;
  for (const match of maskJavaScript(source).matchAll(pattern)) {
    yield /** @type {string} */ (match[1] || match[2]);
  }
}

/**
 * @param {string} source
 * @param {string} name a fixture file name
 */
const namesFixture = (source, name) => new RegExp(`(?<![\\w.-])${name.replace(/[.*+?^${}()|[\]\\]/g, '\\$&')}`).test(source);

/**
 * @param {TwinSources} sources
 * @returns {DetectedTwin[]}
 */
export const detectTwins = ({ tests, python, services, fixtures }) => {
  /** @type {[string, string][]} */
  const javaScriptTests = [];
  /** @type {[string, string][]} */
  const pythonTests = [];
  for (const [path, source] of tests) {
    if (/\.(?:mjs|cjs|js)$/.test(path)) javaScriptTests.push([path, maskJavaScript(source, { strings: false })]);
    else if (path.endsWith('.py')) pythonTests.push([path, maskPython(source, { strings: false })]);
  }
  /** @type {DetectedTwin[]} */
  const detected = [];
  for (const fixture of fixtures) {
    const name = fixture.slice(fixture.lastIndexOf('/') + 1);
    const javascript = javaScriptTests.filter(([, source]) => namesFixture(source, name)).map(([path]) => path);
    const pythonReaders = pythonTests.filter(([, source]) => namesFixture(source, name)).map(([path]) => path);
    if (javascript.length && pythonReaders.length) {
      detected.push({ id: `vector:${fixture}`, python: unique(pythonReaders), javascript: unique(javascript) });
    }
  }
  const pythonByName = watchedDefinitions(python, pythonDefinitions);
  const javaScriptByName = watchedDefinitions(services, javaScriptDefinitions);
  for (const [key, javascript] of javaScriptByName) {
    const pythonDefinitionsOfName = pythonByName.get(key);
    if (!pythonDefinitionsOfName) continue;
    const names = unique(javascript.map((definition) => definition.slice(definition.lastIndexOf(':') + 1)));
    detected.push({ id: `function:${names[0]}`, python: unique(pythonDefinitionsOfName), javascript: unique(javascript) });
  }
  return detected.sort((left, right) => left.id.localeCompare(right.id));
};

/**
 * Reads the twin sources of a checkout from the files Git tracks.
 * @param {string} root
 * @returns {TwinSources}
 */
export const readTwinSources = (root) => {
  const tracked = execFileSync('git', ['ls-files', '-z', '--', 'tests', 'gbdraw'], { cwd: root, encoding: 'utf8' })
    .split('\0').filter(Boolean);
  /** @param {(path: string) => boolean} keep */
  const read = (keep) => new Map(tracked.filter(keep).map((path) => [path, readFileSync(join(root, path), 'utf8')]));
  return {
    tests: read((path) => path.startsWith('tests/') && /\.(?:mjs|cjs|js|py)$/.test(path)),
    python: read((path) => path.startsWith('gbdraw/') && !path.startsWith('gbdraw/web/') && path.endsWith('.py')),
    services: read((path) => path.startsWith('gbdraw/web/js/services/') && /\.m?js$/.test(path)),
    fixtures: tracked.filter((path) => /^tests\/fixtures\/(?:[^/]+\/)*[^/]+\.json$/.test(path)).sort()
  };
};

const ROW_ID = /^(?:vector|function|constant|concern):\S+$/;
const DETECTED_ROW_ID = /^(?:vector|function):/;
/** @param {unknown} value */
const isText = (value) => typeof value === 'string' && /\S/.test(value);
/** @param {unknown} value */
const isPathList = (value) => Array.isArray(value) && value.length > 0 && value.every(isText);

/**
 * Compares the registry with the detected twins: every detected twin has a
 * row and every `vector:` or `function:` row is detected, except a pending
 * row, which is detected only once its twin is added and then loses `pending`.
 * `constant:` and `concern:` rows record twins the detector does not see.
 * @param {unknown} registry the parsed tools/twins.json
 * @param {readonly DetectedTwin[]} detected
 * @returns {TwinRegistryCheck}
 */
export const checkTwinRegistry = (registry, detected) => {
  /** @type {string[]} */
  const errors = [];
  const rows = registry && typeof registry === 'object' && !Array.isArray(registry)
    ? /** @type {{ schemaVersion?: unknown, twins?: unknown }} */ (registry)
    : null;
  if (!rows || rows.schemaVersion !== 1) errors.push(`${TWIN_REGISTRY_PATH}: schemaVersion must be 1`);
  const entries = Array.isArray(rows?.twins) ? /** @type {unknown[]} */ (rows.twins) : null;
  if (!entries) errors.push(`${TWIN_REGISTRY_PATH}: twins must be an array`);
  /** @type {Set<string>} */
  const ids = new Set();
  /** @type {Set<string>} */
  const pendingIds = new Set();
  (entries || []).forEach((entry, index) => {
    const location = `${TWIN_REGISTRY_PATH}: twins[${index}]`;
    const row = /** @type {Partial<Record<keyof TwinRow, unknown>>} */ (entry && typeof entry === 'object' ? entry : {});
    if (typeof row.id !== 'string' || !ROW_ID.test(row.id)) {
      errors.push(`${location}.id: expected vector:, function:, constant: or concern: and a name`);
    } else if (ids.has(row.id)) errors.push(`${location}.id: duplicate ${row.id}`);
    else ids.add(row.id);
    for (const field of /** @type {const} */ (['concern', 'reason', 'retireWhen', 'source'])) {
      if (!isText(row[field])) errors.push(`${location}.${field}: expected text`);
    }
    for (const field of /** @type {const} */ (['pythonOwner', 'jsOwner'])) {
      if (!isPathList(row[field])) errors.push(`${location}.${field}: expected a non-empty list of paths`);
    }
    if (row.vector !== null && !(typeof row.vector === 'string' && /^tests\/fixtures\/(?:[^/]+\/)*[^/]+\.json$/.test(row.vector))) {
      errors.push(`${location}.vector: expected tests/fixtures/<path>.json or null`);
    }
    if (row.pending === undefined) return;
    if (row.pending === true && typeof row.id === 'string' && DETECTED_ROW_ID.test(row.id)) pendingIds.add(row.id);
    else errors.push(`${location}.pending: expected true on a vector: or function: row, or no pending`);
  });
  const detectedIds = new Set(detected.map(({ id }) => id));
  errors.push(...[...pendingIds].filter((id) => detectedIds.has(id)).sort()
    .map((id) => `${TWIN_REGISTRY_PATH}: ${id}: the twin is detected; remove pending`));
  return {
    errors,
    unregistered: [...detectedIds].filter((id) => !ids.has(id)).sort(),
    stale: [...ids].filter((id) => DETECTED_ROW_ID.test(id) && !detectedIds.has(id) && !pendingIds.has(id)).sort()
  };
};
