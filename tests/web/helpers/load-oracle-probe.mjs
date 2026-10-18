// The Load oracle probe: the Web Load of one Session file in node, then Save
// Session, written as a canonical, deterministic dump. Request schemas 1 and
// 2 are held to it: what Load produces from a legacy request must not change
// when the request is read another way (PR-D).
//
//   node tests/web/helpers/load-oracle-probe.mjs [--tree <root>] [--validator strict|migrated] <input> [<out.json>]
//
// `--tree` is a checkout (default: this repository) whose `gbdraw/`,
// `tests/web/helpers/` and `tests/web/fake-svg-dom.mjs` are loaded, so one
// probe measures any commit. `<input>` is a Session file (.json or .json.gz)
// or `git:<blob>[:<name>]`, read with `git cat-file` in the current
// directory's repository; Load reads its bytes as the browser does (gzip
// included). One process per input: Load state is module-global.
//
// The dump: `load` (status, error, the dialogs the Load port was asked, the
// repairs, the after-Load notices); `state`, what a Load replaces, by area
// (both drawings' settings and edits, the shown mode and Results, the title,
// the input files and the other document state); `committed`, the committed
// request's schema and SHA-256; `save`, Save Session after Load (status,
// error, the SHA-256 of its bytes, its requests, and its other parts); and
// `areas`, the SHA-256 of each state area and saved part. Long strings,
// bytes and Files appear as SHA-256, never as payloads. The clock is fixed
// and `crypto.randomUUID` counts, so the dump is the same on every run.
//
// `--validator migrated` runs Python's decoder migration of a stored full
// `config` in the Load's configuration check of a tree older than PR-D,
// which refuses `canvas.show_labels` there (OV-790): it shows what the rest
// of that Load does with such a request. From PR-D the check runs it itself.
import { readFileSync, writeFileSync } from 'node:fs';
import { execFileSync } from 'node:child_process';
import { gunzipSync } from 'node:zlib';
import { createHash } from 'node:crypto';
import { File, Blob } from 'node:buffer';
import { fileURLToPath, pathToFileURL } from 'node:url';
import { resolve } from 'node:path';

const args = process.argv.slice(2);
const option = (name, fallback) => {
  const index = args.indexOf(name);
  return index < 0 ? fallback : args.splice(index, 2)[1];
};
const TREE = resolve(option('--tree', fileURLToPath(new URL('../../../', import.meta.url))));
const VALIDATOR_MODE = option('--validator', 'strict');
const [inputArg, outArg] = args;
if (!inputArg || !['strict', 'migrated'].includes(VALIDATOR_MODE)) {
  console.error('usage: load-oracle-probe.mjs [--tree <root>] [--validator strict|migrated] <input> [<out.json>]');
  process.exit(2);
}
const treeUrl = (path) => pathToFileURL(`${TREE}/${path}`).href;
const sha256 = (bytes) => createHash('sha256').update(bytes).digest('hex');

const FIXED_TIME = Date.UTC(2026, 0, 1);
const RealDate = Date;
globalThis.Date = class extends RealDate {
  constructor(...values) { super(...(values.length ? values : [FIXED_TIME])); }
  static now() { return FIXED_TIME; }
};
let uuidCounter = 0;
Object.defineProperty(globalThis.crypto, 'randomUUID', {
  configurable: true,
  value: () => `00000000-0000-4000-8000-${String(++uuidCounter).padStart(12, '0')}`
});

// The fake DOM of tests/web/load-repair-dialog.test.mjs.
globalThis.window = {
  Vue: {
    ref: (value) => ({ value }), reactive: (value) => value,
    computed: (getter) => ({ get value() { return getter(); } }),
    nextTick: async () => {}
  },
  DOMPurify: { sanitize: (value) => value }
};
globalThis.document = { createElement: () => ({ click: () => {} }) };
const { installFakeSvgDom } = await import(treeUrl('tests/web/fake-svg-dom.mjs'));
installFakeSvgDom();
globalThis.File = File;
globalThis.alert = () => {};
const { installSessionImportWorker } = await import(treeUrl('tests/web/helpers/session-import-node.mjs'));
installSessionImportWorker();

const { importSession, exportSession, setUnmanagedConfigOverrideValidator, getCommittedCanonicalSession } = await import(
  treeUrl('gbdraw/web/js/services/config.js')
);
const validatorCalls = [];
// The Worker's configuration check, run by Python (as load-repair-dialog does).
setUnmanagedConfigOverrideValidator((payload) => {
  const result = JSON.parse(execFileSync('python', ['-c', `
import json, sys
from gbdraw.web_support.config_overrides import validate_web_config_overrides_json
p = json.load(sys.stdin)
if p.get('migrateConfig') and isinstance(p.get('config'), dict):
    import gbdraw.session_request_codec as codec
    if hasattr(codec, '_migrate_legacy_full_config'):
        p['config'] = codec._migrate_legacy_full_config(p['config'], mode=p['mode'], path='renderRequest.diagramOptions.config')
print(validate_web_config_overrides_json(p['mode'], json.dumps(p['config']),
    json.dumps(p['configOverrides']), json.dumps(p['managedPaths']), p['source'] if 'source' in p else p['requireUnmanagedOnly']))
`], { input: JSON.stringify({ ...payload, migrateConfig: VALIDATOR_MODE === 'migrated' }), encoding: 'utf8', cwd: TREE }));
  validatorCalls.push({ mode: payload.mode, result: canon(result) });
  return { result };
});
const { transformLegacyResultSvg } = await import(treeUrl('gbdraw/web/js/app/app-setup.js'));
const { state } = await import(treeUrl('gbdraw/web/js/state.js'));

const readInput = (arg) => {
  if (arg.startsWith('git:')) {
    const [, blob, name = ''] = arg.split(':');
    return { label: arg, name, raw: execFileSync('git', ['cat-file', 'blob', blob], { maxBuffer: 1 << 30 }) };
  }
  return { label: arg, name: arg.split('/').pop(), raw: readFileSync(arg) };
};
const source = readInput(inputArg);
const isGzip = (bytes) => bytes[0] === 0x1f && bytes[1] === 0x8b;
const text = (isGzip(source.raw) ? gunzipSync(source.raw) : source.raw).toString('utf8');
const parsedInput = JSON.parse(text);

// Canonical form: sorted keys; long strings and bytes as SHA-256; Files as
// name, size, type, SHA-256; Maps and Sets as sorted entries.
const LONG = 256;
const fileDigests = [];
const byJson = (a, b) => (JSON.stringify(a) < JSON.stringify(b) ? -1 : 1);
function canon(value, stack = new Set()) {
  if (value === undefined) return { '#undefined': true };
  if (value === null || typeof value === 'boolean') return value;
  if (typeof value === 'number') return Number.isFinite(value) ? value : { '#number': String(value) };
  if (typeof value === 'bigint') return { '#bigint': String(value) };
  if (typeof value === 'string') return value.length > LONG ? { '#string': sha256(value), length: value.length } : value;
  if (typeof value === 'function' || typeof value === 'symbol') return { '#function': true };
  if (stack.has(value)) return { '#cycle': true };
  if (value instanceof Blob || value instanceof File) {
    const entry = { '#file': value.name ?? null, size: value.size, type: value.type || '', sha256: null };
    fileDigests.push([entry, value]);
    return entry;
  }
  if (ArrayBuffer.isView(value) || value instanceof ArrayBuffer) {
    const bytes = value instanceof ArrayBuffer ? new Uint8Array(value) : new Uint8Array(value.buffer, value.byteOffset, value.byteLength);
    return { '#bytes': sha256(bytes), length: bytes.length };
  }
  stack.add(value);
  try {
    if (value instanceof Map) return { '#map': [...value.entries()].map(([k, v]) => [canon(k, stack), canon(v, stack)]).sort((a, b) => byJson(a[0], b[0])) };
    if (value instanceof Set) return { '#set': [...value].map((item) => canon(item, stack)).sort(byJson) };
    if (Array.isArray(value)) return value.map((item) => canon(item, stack));
    const out = {};
    for (const key of Object.keys(value).sort()) {
      let child;
      try { child = value[key]; } catch (error) { out[key] = { '#throws': String(error?.message ?? error) }; continue; }
      if (typeof child === 'function') continue;
      out[key] = canon(child, stack);
    }
    return out;
  } finally {
    stack.delete(value);
  }
}
const resolveFileDigests = async () => {
  for (const [entry, file] of fileDigests.splice(0)) entry.sha256 = sha256(new Uint8Array(await file.arrayBuffer()));
};

// The document state a Load sets (state.js), by area. A drawing's derived
// members (computed from the others) and its slice handles (no state) are
// left out.
const DERIVED_DRAWING_KEYS = new Set([
  'activeLayoutPreferences', 'drawnLinearSeqs', 'linearComparisonResolution', 'hasLinearComparisonIntent',
  'hasActiveLinearLosatIntent', 'hasActiveLinearUploadIntent', 'featureVisibilityRules', 'hasPendingPaletteDraft',
  'slices'
]);
const DOCUMENT_KEYS = [
  'mode', 'sessionTitle', 'results', 'selectedResultIndex', 'cInputType', 'lInputType', 'files', 'linearSeqs',
  'losatCacheInfo', 'losatCache', 'losatDerivedCache', 'proteinIdentityManifest', 'legacyProteinRawCandidates',
  'legacyProteinDerivedEvidence', 'orthogroups', 'collinearGroups', 'similarityAlignmentPlan',
  'similarityAlignmentResetReceipt', 'linearRecordTranslations', 'legacySimilarityAlignment', 'selectedOrthogroupId',
  'selectedOrthogroupAlignmentFeature', 'circularRecordList', 'paletteInstantPreviewEnabled', 'losatExecution',
  'richFeaturePopup', 'appliedPaletteName', 'appliedPaletteColors', 'downloadDpi', 'featureCatalog',
  'biologicalFeatures', 'extractedFeatures', 'featureRecordIds', 'selectedFeatureRecordIdx', 'autoLabelReflowEnabled',
  'generatedLegendPosition', 'generatedMode', 'generatedMultiRecordCanvas', 'generatedCircularPlotTitlePosition',
  'featurePanelTab', 'lastRunInfo', 'zoom', 'canvasPan', 'displayedResultMetadata', 'activeDrawing',
  // What Load reports and the Result artifacts it restores.
  'errorLog', 'comparisonWarnings', 'annotationWarnings', 'featureIdentityNotices', 'originalLegendOrder',
  'originalLegendColors', 'originalSvgStroke', 'matchSequenceRegistry', 'pairwiseMatchFactors', 'featureOrthogroupIndex',
  'linearReorderNotice', 'circularRecordDiscovery', 'suppressCircularMultiRecordDefaults', 'featureEditRemovalCount',
  'generationFailureRecovery', 'svgContent'
];
const stateAreas = () => {
  const areas = {};
  for (const key of DOCUMENT_KEYS) areas[`state.${key}`] = key in state ? canon(state[key]) : { '#absent': true };
  for (const mode of ['circular', 'linear']) {
    const drawing = state.drawings[mode];
    for (const key of Object.keys(drawing).sort()) {
      if (!DERIVED_DRAWING_KEYS.has(key)) areas[`drawings.${mode}.${key}`] = canon(drawing[key]);
    }
    // Not enumerable: the layout preferences own them.
    areas[`drawings.${mode}.form.legend`] = canon(drawing.form.legend);
    areas[`drawings.${mode}.adv.plot_title_position`] = canon(drawing.adv.plot_title_position);
  }
  return areas;
};
const requestSummary = (request) => (request && typeof request === 'object'
  ? { schema: request.schema ?? null, mode: request.mode ?? null, sha256: sha256(JSON.stringify(canon(request))) }
  : null);
const errorSummary = (error) => (error ? {
  raw: canon(JSON.parse(JSON.stringify(error, (key, value) => (
    value instanceof Error ? { name: value.name, message: value.message, ...value } : value)) ?? 'null')),
  code: error.code ?? null, operation: error.operation ?? null, stage: error.stage ?? null,
  message: error.message ?? error.detail ?? null, title: error.title ?? null
} : null);

const dump = {
  input: {
    label: source.label, name: source.name, bytes: source.raw.length, sha256: sha256(source.raw),
    sessionVersion: parsedInput?.version ?? null,
    request: requestSummary(parsedInput?.renderRequest),
    otherModeRequest: requestSummary(parsedInput?.otherModeResult?.renderRequest)
  },
  tree: (() => {
    try { return readFileSync(`${TREE}/COMMIT`, 'utf8').trim(); } catch { return null; }
  })(),
  validatorMode: VALIDATOR_MODE
};

const dialogs = [];
let loaded;
try {
  const file = isGzip(source.raw)
    ? new File([source.raw], 'session.json.gz', { type: 'application/gzip' })
    : new File([source.raw], 'session.json', { type: 'application/json' });
  loaded = await importSession({ target: { files: [file], value: 'selected' } }, {
    transformLegacyResultSvg,
    confirmLoadRepairs: async (dialog) => { dialogs.push(canon(dialog)); return true; }
  });
} catch (error) {
  loaded = { status: 'threw', error: { message: String(error?.message ?? error) } };
}
dump.load = {
  status: loaded?.status ?? null,
  error: errorSummary(loaded?.error),
  dialogs,
  repairs: Array.isArray(loaded?.repairs) ? loaded.repairs.map((repair) => canon(repair)) : null,
  notices: canon(loaded?.notices ?? null),
  warnings: canon(loaded?.warnings ?? null),
  degradedRecovery: loaded?.degradedRecovery ?? null,
  comparisonDisposition: loaded?.comparisonDisposition ?? null,
  validatorCalls
};
const committed = getCommittedCanonicalSession?.() ?? null;
dump.committed = committed ? {
  request: requestSummary(committed.renderRequest),
  otherModeRequest: requestSummary(committed.otherModeResult?.renderRequest ?? committed.otherRenderRequest)
} : null;
const areas = stateAreas();
await resolveFileDigests();

let saved;
try {
  const result = await exportSession(state.sessionTitle.value || 'probe', { askChoice: async () => 'download', onError: () => {} });
  if (result?.status === 'saved') {
    const compressed = new Uint8Array(await result.blob.arrayBuffer());
    const raw = isGzip(compressed) ? gunzipSync(compressed) : Buffer.from(compressed);
    const document = JSON.parse(raw.toString('utf8'));
    const request = document.renderRequest;
    const otherRequest = document.otherModeResult?.renderRequest;
    delete document.renderRequest;
    if (document.otherModeResult && 'renderRequest' in document.otherModeResult) delete document.otherModeResult.renderRequest;
    saved = {
      status: 'saved', sha256: sha256(raw), bytes: raw.length,
      request: requestSummary(request), otherModeRequest: requestSummary(otherRequest),
      parts: Object.fromEntries(Object.keys(document).sort().map((key) => [key, canon(document[key])]))
    };
  } else {
    saved = { status: result?.status ?? null, error: errorSummary(result?.error) };
  }
} catch (error) {
  saved = { status: 'threw', error: { message: String(error?.message ?? error), code: error?.code ?? null } };
}
await resolveFileDigests();
dump.save = saved;
dump.state = areas;
const digest = (value) => sha256(JSON.stringify(value));
dump.areas = {
  ...Object.fromEntries(Object.entries(areas).map(([key, value]) => [key, digest(value)])),
  ...(saved.parts ? Object.fromEntries(Object.entries(saved.parts).map(([key, value]) => [`save.${key}`, digest(value)])) : {})
};
const output = `${JSON.stringify(dump, null, 1)}\n`;
// The import Worker keeps the process alive: exit once the dump is written.
if (outArg) {
  writeFileSync(outArg, output);
  process.exit(0);
} else {
  process.stdout.write(output, () => process.exit(0));
}
