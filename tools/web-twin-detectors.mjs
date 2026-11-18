// @ts-check
// Twin detector (boundary audit M2, plan row P1-D): logic that Python and the
// Web each implement. A twin is allowed only with a row in tools/twins.json
// that names its owners, the vector that pins the two to each other, and the
// condition under which one side is deleted. The registry only shrinks: a new
// row is an authority change, made before its twin exists with `pending`
// (gbdraw/web/CLAUDE.md R9); the change that adds the twin removes `pending`.
// A row the trusted base does not have is admitted only with the fields of
// `checkTwinAdditions`, and only when retired rows pay for it or it carries a
// waiver (WEB_CHANGE_POLICY.md "Design-rule co-change"); the Gate
// (tools/check-web-change-budget.mjs) runs that check against the base.
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
 * @property {string} [whyNotAskPython] Why the Web does not ask the Python owner (in the Worker) instead.
 * @property {string} [whyNotGenerate] Why the Web copy is not generated from the Python owner.
 * @property {string} [retireBy] The date (YYYY-MM-DD) by which the twin retires.
 * @property {string} [waiver] Why no row could be retired first; the Owner grants it by approving the adding pull request.
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
/** @param {unknown} value a calendar date written YYYY-MM-DD */
const isDate = (value) => typeof value === 'string' && /^\d{4}-\d{2}-\d{2}$/.test(value)
  && !Number.isNaN(Date.parse(`${value}T00:00:00Z`)) && new Date(`${value}T00:00:00Z`).toISOString().startsWith(value);
const ADMISSION_TEXT_FIELDS = /** @type {const} */ (['whyNotAskPython', 'whyNotGenerate']);

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
    for (const field of /** @type {const} */ ([...ADMISSION_TEXT_FIELDS, 'waiver'])) {
      if (row[field] !== undefined && !isText(row[field])) errors.push(`${location}.${field}: expected text`);
    }
    if (row.retireBy !== undefined && !isDate(row.retireBy)) errors.push(`${location}.retireBy: expected a date YYYY-MM-DD`);
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

// Rows retired from this time on pay for rows added from this time on: the merge of the admission
// rule into WEB_CHANGE_POLICY.md "Design-rule co-change" (#2089).
export const TWIN_EXCHANGE_SINCE = '2026-11-17T06:32:55+09:00';

/**
 * @param {unknown} registry a parsed tools/twins.json, or null
 * @returns {Map<string, Partial<Record<keyof TwinRow, unknown>>>} rows by id
 */
const rowsById = (registry) => {
  const twins = registry && typeof registry === 'object' && !Array.isArray(registry)
    ? /** @type {{ twins?: unknown }} */ (registry).twins
    : null;
  /** @type {Map<string, Partial<Record<keyof TwinRow, unknown>>>} */
  const rows = new Map();
  for (const entry of Array.isArray(twins) ? twins : []) {
    const row = /** @type {Partial<Record<keyof TwinRow, unknown>>} */ (entry && typeof entry === 'object' ? entry : {});
    if (typeof row.id === 'string') rows.set(row.id, row);
  }
  return rows;
};

/**
 * Checks the rows the head adds: ids that none of `baseRegistries` has (the
 * merge base, so a row the base branch retired after the fork is not new, and
 * the trusted base, so a row it already has is not new). Each says why the Web
 * neither asks the Python owner nor generates the value, and gives a
 * `retireBy` date after the day of the check (UTC). The rows without a waiver
 * must not exceed the balance of retired rows (`readTwinExchangeBalance`).
 * @param {unknown[]} baseRegistries the parsed registries of the merge base and the trusted base; null where absent
 * @param {unknown} headRegistry the parsed registry of the head
 * @param {{ now: Date, balance: () => number }} options `balance` is called only when an added row has no waiver
 * @returns {string[]}
 */
export const checkTwinAdditions = (baseRegistries, headRegistry, { now, balance }) => {
  const before = new Set(baseRegistries.flatMap((registry) => [...rowsById(registry).keys()]));
  const today = now.toISOString().slice(0, 10);
  const added = [...rowsById(headRegistry)].filter(([id]) => !before.has(id));
  const errors = added.flatMap(([id, row]) => {
    const missing = [
      ...ADMISSION_TEXT_FIELDS.filter((field) => !isText(row[field])),
      ...(isDate(row.retireBy) && /** @type {string} */ (row.retireBy) > today ? [] : [`retireBy (a date after ${today}, UTC)`])
    ];
    return missing.length ? [`${TWIN_REGISTRY_PATH}: ${id}: an added row needs ${missing.join(', ')}`] : [];
  });
  const unpaid = added.filter(([, row]) => !isText(row.waiver)).length;
  const available = unpaid ? balance() : 0;
  if (unpaid > available) {
    errors.push(`${TWIN_REGISTRY_PATH}: ${unpaid} added row(s) without a waiver, but the rows that pull requests changing `
      + `production code retired since ${TWIN_EXCHANGE_SINCE.slice(0, 10)} leave ${available} to pay for additions: retire `
      + 'rows first, or give each added row a "waiver" that says why none could be retired (the Owner grants it by '
      + 'approving this pull request)');
  }
  return errors;
};

const NO_OBJECT = /^0+$/;

/**
 * Reads blobs with one `git cat-file --batch`. A blob that cannot be read
 * throws: a failed read must not look like an absent registry.
 * @param {(args: string[], input?: string) => string} git
 * @param {string[]} oids
 * @returns {Map<string, string>} source by blob id
 */
const readBlobs = (git, oids) => {
  /** @type {Map<string, string>} */
  const sources = new Map();
  if (!oids.length) return sources;
  // A partial clone (the trusted-base Gate fetches no blobs) would fetch each missing blob on its first
  // read, one request each; fetch them in one request first. Blobs already present are not fetched again.
  // Its remote is marked `remote.<name>.promisor` (older clones name it in `extensions.partialClone`).
  const promisor = git(['config', '--list']).split('\n').map((line) => (
    /^remote\.(.+)\.promisor=(?:true|yes|on|1)$/i.exec(line)?.[1] ?? /^extensions\.partialclone=(.+)$/i.exec(line)?.[1]
  )).find(Boolean);
  if (promisor) {
    git([
      '-c', 'fetch.negotiationAlgorithm=noop', 'fetch', '--quiet', '--no-tags', '--no-write-fetch-head',
      '--recurse-submodules=no', '--filter=blob:none', '--stdin', promisor
    ], `${oids.join('\n')}\n`);
  }
  // The registry is UTF-8 text, so the bytes come back unchanged and the sizes in the headers hold.
  const output = Buffer.from(git(['cat-file', '--batch'], `${oids.join('\n')}\n`), 'utf8');
  let offset = 0;
  for (const oid of oids) {
    const headerEnd = output.indexOf(10, offset);
    const header = output.subarray(offset, headerEnd).toString('utf8');
    const size = Number(header.split(' ')[2]);
    if (!header.startsWith(`${oid} blob `) || !Number.isInteger(size)) throw new Error(`Cannot read blob ${oid}: ${header}`);
    sources.set(oid, output.subarray(headerEnd + 1, headerEnd + 1 + size).toString('utf8'));
    offset = headerEnd + 1 + size + 1;
  }
  return sources;
};

/**
 * The balance of retired rows on the first-parent history of `revision` since
 * TWIN_EXCHANGE_SINCE: each row that a change to production code removed
 * counts one, and each row added without a waiver takes one. A row removed by
 * a change without production code (an abandoned pending row) pays for
 * nothing.
 * @param {(args: string[], input?: string) => string} git runs git in the repository with `input` on its
 *   standard input and returns its output; throws on failure
 * @param {string} revision
 * @param {(path: string) => boolean} isProduction
 * @returns {number}
 */
export const readTwinExchangeBalance = (git, revision, isProduction) => {
  // Each first-parent commit that changed the registry, with the registry blobs before and after it as
  // the trees name them. Each commit is kept or dropped by its own date: a date limit on the walk would
  // also drop the commits behind one dated out of order.
  const since = Date.parse(TWIN_EXCHANGE_SINCE);
  const changes = git([
    'log', '--first-parent', '--diff-merges=first-parent', '--no-renames', '--raw', '--no-abbrev',
    '--format=%x00%H %cI', revision, '--', TWIN_REGISTRY_PATH
  ]).split('\0').slice(1).flatMap((entry) => {
    const [header, ...lines] = entry.split('\n');
    const [commit, date] = header.split(' ');
    const raw = lines.find((line) => line.startsWith(':'));
    if (!raw) throw new Error(`git log listed ${commit} without a change of ${TWIN_REGISTRY_PATH}`);
    const [, , before, after] = raw.slice(1).split(/\s+/);
    return Date.parse(date) >= since ? [{ commit, before, after }] : [];
  });
  const sources = readBlobs(git, [...new Set(changes.flatMap(({ before, after }) => [before, after]))]
    .filter((oid) => !NO_OBJECT.test(oid)));
  /** @param {string} oid */
  const registry = (oid) => (NO_OBJECT.test(oid) ? new Map() : rowsById(JSON.parse(/** @type {string} */ (sources.get(oid)))));
  return changes.reduce((balance, { commit, before: beforeOid, after: afterOid }) => {
    const before = registry(beforeOid);
    const after = registry(afterOid);
    const unpaid = [...after].filter(([id, row]) => !before.has(id) && !isText(row.waiver)).length;
    const retired = [...before.keys()].filter((id) => !after.has(id)).length;
    const paid = retired && git(['diff', '--name-only', '--no-renames', `${commit}^1`, commit]).split('\n').some(isProduction)
      ? retired
      : 0;
    return balance + paid - unpaid;
  }, 0);
};
