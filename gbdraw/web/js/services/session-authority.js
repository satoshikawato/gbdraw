// @ts-check
import { validateAnnotationWarnings } from './session-feature-metadata.js';
import { validateComparisonWarnings } from './comparison-warnings.js';
import { assertSafeObjectKeys } from './safe-object-keys.js';
import {
  POSITIONAL_BINDING_LISTS, readSessionResources, validateWebFileBindings, webBindingInputSlot, webFileBindingFaults,
  webFileBindingPath
} from './session-resource-backing.js';
import { SESSION_INPUT_TYPES, alignmentResetReceiptShapeHolds, requireActiveConfigContainers } from './session-active-config-contract.js';
import { migrateLegacyLinearLabelVisibility } from './linear-label-visibility.js';
import { migrateLegacyRecordDisplayDrafts } from './record-display-model.js';
import { canonicalFeatureOverrides, validateFeatureIdentityNotices } from './feature-placement.js';
import { RENDERED_ID_FEATURE_EDIT_FIELDS, migrateSessionFeaturePlacements } from './feature-edit-migration.js';
import { collectCanonicalResourceIds } from './canonical-resource-references.js';
import { LOSAT_EXECUTION_FIELDS, MODE_SCOPED_ROWS, SLICE_MODES, validateModeSliceFields } from './mode-scoped-migration.js';
import { diagnosticError } from '../utils/error-normalization.js';
import { normalizeLinearComparisonPlan } from './linear-comparisons.js';
import { SESSION_ADMISSION_DESCRIPTION } from '../contracts/session-shape.generated.js';
import { repairByDescription } from '../utils/contract-repair.js';
import { isPlainObject } from '../utils/contract-shape.js';
import { SESSION_PART_PATHS, loadRepair } from '../utils/load-repairs.js';
import {
  CANONICAL_SESSION_MIN_VERSION,
  CURRENT_ARTIFACT_SESSION_MIN_VERSION,
  CURRENT_AUTHORITY_SESSION_MIN_VERSION,
  CURRENT_WRITER_FORBIDDEN_FEATURE_FIELDS,
  FEATURE_CATALOG_SCHEMA_BY_SESSION_VERSION,
  FLAT_DRAFT_TOP_LEVEL_FIELDS,
  MODE_SCOPED_SESSION_VERSION,
  SESSION_DRAWINGS_LAYOUT,
  OTHER_MODE_RESULT_EDITOR_FIELDS,
  OTHER_MODE_RESULT_FIELDS,
  OTHER_MODE_RESULT_LEGEND_FIELDS,
  OTHER_MODE_RESULT_RUN_METADATA_FIELDS,
  OTHER_MODE_RESULT_UI_FIELDS,
  ROW_KEYED_LEGEND_SESSION_VERSION,
  SESSION_VERSION,
  SUPPORTED_SESSION_VERSIONS,
  TYPED_DRAFT_SESSION_VERSION
} from '../contracts/schema.generated.js';
/**
 * @import { SessionChecked, SessionChecked31, SessionChecked39, SessionChecked40, SessionChecked46, SessionChecked49,
 *   SessionEnvelope, SessionInput, SessionResult } from '../contracts/session.generated.js'
 */
/** @import { RawSessionResourceDescriptor } from './session-resource-backing.js' */
/** @import { AlignmentResetRequest } from './session-active-config-contract.js' */
/** @import { ContractRepair, RepairPolicy } from '../utils/contract-repair.js' */
/** @import { ContractShape } from '../utils/contract-shape.js' */
/** @import { LoadRepair, Repaired } from '../utils/load-repairs.js' */

/**
 * How code reads an unvalidated value: through its properties `K`, each of
 * them another unvalidated value (`unknown`, or such a view again). Every
 * JavaScript value has this type, so a cast to it claims nothing that a check
 * would have to prove.
 * @template {string} K
 * @template [V=unknown]
 * @typedef {{ readonly [P in K]?: V } | null | undefined} UnvalidatedFields
 */

// `TYPED_DRAFT_SESSION_VERSION` (44) introduced the current active-config and
// record-display draft shapes; `MODE_SCOPED_SESSION_VERSION` (46) keeps the
// draft of each diagram mode in `modes` (PD-OI-086) and keys per-feature
// edits by source identity (design Q4).

// One class per top-level key of a Session (the keys are `SessionInput`'s),
// and the content part it holds (receipt 〔拒否〕: settings, request, Result,
// input files, edits, cache), or null for a key that holds no content: R4
// refuses a file whose every held part Load dropped (`restoresNothing`).
/**
 * @typedef {'settings' | 'request' | 'results' | 'input-files' | 'edits' | 'cache'} SessionContentPart
 * @typedef {{ readonly authority: string, readonly part: SessionContentPart | null }} SessionTopLevelKey
 */
export const SESSION_TOP_LEVEL_AUTHORITY = Object.freeze(/** @satisfies {Record<keyof SessionInput, SessionTopLevelKey>} */ ({
  format: { authority: 'document', part: null },
  version: { authority: 'document', part: null },
  createdAt: { authority: 'document', part: null },
  title: { authority: 'document', part: null },
  renderRequest: { authority: 'canonical-render', part: 'request' },
  resources: { authority: 'resource', part: 'input-files' },
  webFiles: { authority: 'resource-binding', part: 'input-files' },
  config: { authority: 'legacy-or-editor-metadata', part: 'settings' },
  ui: { authority: 'editor-metadata', part: null },
  files: { authority: 'legacy-fallback', part: 'input-files' },
  results: { authority: 'artifact', part: 'results' },
  features: { authority: 'artifact-or-legacy-semantic', part: 'edits' },
  editorState: { authority: 'artifact', part: 'edits' },
  orthogroupState: { authority: 'artifact', part: 'edits' },
  losatCache: { authority: 'artifact', part: 'cache' },
  losatDerivedCache: { authority: 'artifact', part: 'cache' },
  proteinIdentityManifest: { authority: 'artifact', part: 'cache' },
  legacyArtifacts: { authority: 'artifact', part: 'cache' },
  runMetadata: { authority: 'artifact', part: null },
  otherModeResult: { authority: 'artifact', part: 'results' },
  cliInvocation: { authority: 'provenance', part: null },
  modes: { authority: 'editor-metadata', part: 'settings' },
  cliOptions: { authority: 'provenance', part: null }
}));

// The homes that Session 46 moved into `modes` (and app-level `ui`): a
// Session 46 that holds any of them is rejected (plan 4.1).
/** @type {ReadonlyArray<[string, ReadonlySet<string>]>} */
const RETIRED_MODE_SCOPED_FIELDS = Object.freeze([
  ['', new Set(FLAT_DRAFT_TOP_LEVEL_FIELDS)],
  ['editorState', new Set(['featureStrokes'])],
  ['editorState.legend', new Set(MODE_SCOPED_ROWS.filter((row) => row.domain === 'editorState.legend').map((row) => row.path))],
  ['ui', new Set(MODE_SCOPED_ROWS.filter((row) => row.domain === 'ui').map((row) => row.path))]
]);

const WEB_EDITOR_UI_FIELDS = Object.freeze([
  'cInputType',
  'lInputType',
  'zoom',
  'canvasPan',
  'canvasPadding',
  'selectedResultIndex',
  'featurePanelTab',
  'downloadDpi',
  'autoLabelReflow',
  'linearTypographyLinked',
  'paletteInstantPreviewEnabled',
  'appliedPaletteName',
  'appliedPaletteColors',
  'pendingPaletteName',
  'pendingPaletteColors',
  'layoutPreferences',
  // Older editor payloads stored the same preference state in parallel fields.
  'legend',
  'circularLegendPosition',
  'linearLegendPosition',
  'circularPlotTitlePosition',
  'linearPlotTitlePosition',
  'circularSingleRecordLegendPosition',
  'circularSingleRecordPlotTitlePosition',
  'circularMultiRecordLegendPosition',
  'circularMultiRecordPlotTitlePosition'
]);

const ARTIFACT_UI_FIELDS = Object.freeze([
  'generatedLegendPosition',
  'generatedMultiRecordCanvas',
  'generatedCircularPlotTitlePosition'
]);

// Older Sessions' rendered-ID edit maps are read once and migrated to
// `featureOverrides` (services/feature-edit-migration.js).
const ARTIFACT_FEATURE_FIELDS = Object.freeze([
  'extractedFeatures',
  'biologicalFeatures',
  'featureRecordIds',
  'selectedFeatureRecordIdx',
  'featureColorOverrides',
  'featureVisibilityManualRules',
  'featureOverrides',
  'labelOverrideRows',
  'labelTextBulkOverrides',
  ...RENDERED_ID_FEATURE_EDIT_FIELDS
]);

const SUPPORTED_SESSION_VERSION_SET = new Set(SUPPORTED_SESSION_VERSIONS);
// The versions this reader takes, in words, for the refusal of another
// version (OV-1253): "27–33, 39–42, 44, and 46–48".
const READ_SESSION_VERSIONS = (() => {
  /** @type {[number, number][]} */
  const runs = [];
  for (const version of SUPPORTED_SESSION_VERSIONS) {
    const last = runs[runs.length - 1];
    if (last && version === last[1] + 1) last[1] = version;
    else runs.push([version, version]);
  }
  const words = runs.map(([first, end]) => (first === end ? String(first) : `${first}\u2013${end}`));
  return words.length > 1 ? `${words.slice(0, -1).join(', ')}, and ${words[words.length - 1]}` : words.join('');
})();
const OLDEST_SESSION_VERSION = Math.min(...SUPPORTED_SESSION_VERSIONS);
// The newest version this app reads; Save may still write an older one (SESSION_VERSION).
const NEWEST_SESSION_VERSION = Math.max(...SUPPORTED_SESSION_VERSIONS);

/**
 * The one reader of a parsed Session file. What it returns is a
 * `gbdraw-session` object of a version this app reads; every other value is
 * unvalidated until the validator of its part has checked it.
 * @param {unknown} value
 * @returns {SessionEnvelope}
 */
export const readSessionInput = (value) => {
  const formatInvalid = () => diagnosticError('INPUT_INVALID', { field: 'schema', reason: 'SESSION_FORMAT' });
  if (!isPlainObject(value)) throw formatInvalid();
  const { version } = value;
  if (typeof version !== 'number' || !Number.isInteger(version)) throw formatInvalid();
  // Q2 of the Load repair receipt: a Session older than the oldest version
  // this app reads is refused, as an old version.
  if (value.format === 'gbdraw-session' && version > 0 && version < OLDEST_SESSION_VERSION) {
    throw diagnosticError('INPUT_INVALID', { field: 'schema', reason: 'SESSION_TOO_OLD' });
  }
  if (!SUPPORTED_SESSION_VERSION_SET.has(version)) {
    // OV-1253: a Session names its version, as the CLI does; another file keeps the generic refusal.
    throw diagnosticError('INPUT_INVALID', value.format !== 'gbdraw-session' ? { field: 'schema', reason: 'FIELDS' }
      : version > NEWEST_SESSION_VERSION
        ? { field: 'schema', reason: 'SESSION_TOO_NEW', sessionVersion: version, latestSessionVersion: NEWEST_SESSION_VERSION }
        : { field: 'schema', reason: 'SESSION_VERSION_UNSUPPORTED', sessionVersion: version, readSessionVersions: READ_SESSION_VERSIONS });
  }
  if (value.format !== 'gbdraw-session') throw formatInvalid();
  return /** @type {SessionEnvelope} */ (value); // `format` and `version` were checked above
};

// The band of a checked Session (`SessionChecked`, describe_session_admission()).
// TypeScript does not narrow a version's literal union by `>=`, so these five
// are the version tests that narrow it; each starts at its band's generated
// threshold (tests/web/session-authority.test.mjs compares them with the bands).
// A Session 49 is a list of drawings (`isDrawingsBand`): none of the four
// bands of the file's sets holds it, so none of them reads it (DESIGN §6).
/** @param {SessionChecked} document @returns {document is SessionChecked49} */
export const isDrawingsBand = (document) => document.version >= SESSION_DRAWINGS_LAYOUT.version;
/** @param {SessionChecked} document @returns {document is SessionChecked31 | SessionChecked39 | SessionChecked40 | SessionChecked46} */
export const isCanonicalBand = (document) => document.version >= CANONICAL_SESSION_MIN_VERSION && !isDrawingsBand(document);
/** @param {SessionChecked} document @returns {document is SessionChecked39 | SessionChecked40 | SessionChecked46} */
export const isArtifactBand = (document) => document.version >= CURRENT_ARTIFACT_SESSION_MIN_VERSION && !isDrawingsBand(document);
/** @param {SessionChecked} document @returns {document is SessionChecked40 | SessionChecked46} */
export const isAuthorityBand = (document) => document.version >= CURRENT_AUTHORITY_SESSION_MIN_VERSION && !isDrawingsBand(document);
/** @param {SessionChecked} document @returns {document is SessionChecked46} */
export const isModeScopedBand = (document) => document.version >= MODE_SCOPED_SESSION_VERSION && !isDrawingsBand(document);

/**
 * A value that holds nothing: `null`, `[]` or `{}` (E-OD06; R4's held parts).
 * @param {unknown} value
 */
export const holdsNothing = (value) => value === undefined || value === null
  || (Array.isArray(value) ? value.length === 0 : isPlainObject(value) && Object.keys(value).length === 0);

/**
 * What the Session reader knows that its description does not (data beside
 * `SESSION_ADMISSION_DESCRIPTION`, D1c 2.3):
 * - `positionalArrays`: arrays other places name by position (a catalog item,
 *   `ui.selectedResultIndex`), whose elements are never dropped one by one;
 * - `slices`: the mode slices, each a part of `modes`;
 * - `requestOwned`: fields the stored request's reader reads, absent ones
 *   included (`readSessionRequests`, K22);
 * - `readerOwnedObjects`: types whose keys and values the reader table reads
 *   (config-field-readers.js): the reader checks only that each is an object;
 * - `readerOwnedValues`: types whose values the reader table reads: the reader
 *   checks their keys and kind;
 * - `resultSetMembers`: what names a Result set's Results, dropped with them;
 *   `features`: the feature state a Session 27-39 recovered from its Results,
 *   which a Session 40+ holds in its catalog (OV-789).
 * Every top-level key is a part (rule 6), and a cache part is dropped with no
 * row (E-OD05).
 */
export const SESSION_PART_RULES = Object.freeze({
  positionalArrays: Object.freeze(['results', 'otherModeResult.results']),
  slices: Object.freeze(SLICE_MODES.map((mode) => `modes.${mode}`)),
  requestOwned: Object.freeze(['renderRequest']),
  readerOwnedObjects: Object.freeze(['CheckedSessionModeSliceConfig']),
  readerOwnedValues: Object.freeze(['CheckedSessionModeSliceFeatures']),
  resultSetMembers: Object.freeze({
    session: Object.freeze(['runMetadata', 'cliInvocation']),
    ui: Object.freeze(['selectedResultIndex', ...ARTIFACT_UI_FIELDS]),
    features: CURRENT_WRITER_FORBIDDEN_FEATURE_FIELDS
  })
});

/** @type {ContractShape} */
const READER_OWNED_VALUE = Object.freeze({ kind: 'unknown' });
/**
 * The admission description as the Session reader applies it: the places the
 * reader table and the request reader read are left to them.
 */
const SESSION_READ_TYPES = Object.freeze(Object.fromEntries(Object.entries(SESSION_ADMISSION_DESCRIPTION.types).map(([name, shape]) => {
  if (shape.kind !== 'object') return [name, shape];
  if (SESSION_PART_RULES.readerOwnedObjects.includes(name)) {
    return [name, Object.freeze({ kind: 'object', fields: Object.freeze({}), optional: Object.freeze([]), open: true })];
  }
  const fields = SESSION_PART_RULES.readerOwnedValues.includes(name)
    ? Object.fromEntries(Object.keys(shape.fields).map((field) => [field, READER_OWNED_VALUE]))
    : shape.fields;
  const optional = [...shape.optional, ...SESSION_PART_RULES.requestOwned.filter((field) => Object.hasOwn(fields, field) && !shape.optional.includes(field))];
  return [name, Object.freeze({ ...shape, fields: Object.freeze(fields), optional: Object.freeze(optional) })];
})));

/** @type {RepairPolicy} */
const SESSION_REPAIR_POLICY = Object.freeze({
  positionalArrays: SESSION_PART_RULES.positionalArrays,
  partBoundary: (/** @type {string} */ pattern) => !pattern.includes('.') || SESSION_PART_RULES.slices.includes(pattern)
});

/** @param {string} path @returns {string} */
const topLevelKey = (path) => path.split(/[.[]/, 1)[0];
/**
 * The stored value at a row's path (`.key` and `[i]` steps; no key of a
 * place the reader repairs holds a dot).
 * @param {unknown} root
 * @param {string} path
 */
const storedAt = (root, path) => path.split('.').flatMap((step) => {
  const [key, ...indices] = step.split('[');
  return [key, ...indices.map((index) => Number(index.slice(0, -1)))];
}).reduce((/** @type {unknown} */ value, step) => (
  value !== null && typeof value === 'object' ? /** @type {Record<string | number, unknown>} */ (value)[step] : undefined
), root);

/**
 * Whether the band holds the top-level `key` (a field it names, not as `never`).
 * @param {string} root The band's root type.
 * @param {string} key
 */
const bandHolds = (root, key) => {
  const shape = SESSION_READ_TYPES[root];
  const field = shape.kind === 'object' ? shape.fields[key] : undefined;
  return field !== undefined && !(field.kind === 'literal' && field.values.length === 0);
};

/**
 * The dialog row of one repair of the description, or null when it lists
 * none: a value that holds nothing dropped from an optional field (rule 1,
 * E-OD06) and a retired home whose place a slice fills (E-OD03, A7f) have
 * none; a cache part's rows are not listed (E-OD05). A part the band holds
 * is named as a part; any other place, as a field.
 * @param {ContractRepair<never>} row
 * @param {SessionEnvelope} stored
 * @param {string} root The band's root type.
 * @returns {LoadRepair | null}
 */
const bandRepair = (row, stored, root) => {
  if (row.rule === 4 && holdsNothing(storedAt(stored, row.path))) return null;
  const key = topLevelKey(row.path);
  if (row.path === key && isModeScopedVersion(stored) && FLAT_DRAFT_TOP_LEVEL_FIELDS.includes(key)) return null;
  const part = Object.hasOwn(SESSION_TOP_LEVEL_AUTHORITY, key)
    ? SESSION_TOP_LEVEL_AUTHORITY[/** @type {keyof SessionInput} */ (key)].part : null;
  // The legacy `files` of a Session 40+ is named with the inputs it held,
  // which load from no other place (A6c).
  const inputs = row.path === 'files' ? heldInputNames(storedAt(stored, 'files')) : [];
  return loadRepair({
    code: part === 'cache' ? 'cache-entry-dropped'
      : bandHolds(root, key) && SESSION_PART_PATHS.includes(row.path) ? 'session-part-dropped' : 'session-field-dropped',
    path: row.path,
    effect: row.effect,
    ...(inputs.length ? { label: inputs.join(', ') } : {})
  });
};

/**
 * The names of the input files a value holds: each `name` text of an object
 * in it, in order, once.
 * @param {unknown} value
 * @returns {string[]}
 */
const heldInputNames = (value) => {
  /** @type {Set<string>} */
  const names = new Set();
  /** @param {unknown} item */
  const visit = (item) => {
    if (Array.isArray(item)) item.forEach(visit);
    else if (isPlainObject(item)) {
      if (typeof item.name === 'string' && item.name.trim()) names.add(item.name.trim());
      Object.values(item).forEach(visit);
    }
  };
  visit(value);
  return [...names];
};

/** @param {SessionEnvelope} document */
const isModeScopedVersion = (document) => document.version >= MODE_SCOPED_SESSION_VERSION;

/**
 * Drops the members of a Result set whose Results the reader dropped.
 * @param {unknown} document
 * @returns {unknown}
 */
const withoutResultSetMembers = (document) => {
  if (!isPlainObject(document)) return document;
  const { session, ui: uiMembers, features: featureMembers } = SESSION_PART_RULES.resultSetMembers;
  const next = Object.fromEntries(Object.entries(document).filter(([key]) => !session.includes(key)));
  if (isPlainObject(document.ui)) {
    next.ui = Object.fromEntries(Object.entries(document.ui).filter(([key]) => !uiMembers.includes(key)));
  }
  if (Object.hasOwn(document, 'features')) next.features = withoutFields(document.features, (field) => featureMembers.includes(field));
  if (isPlainObject(document.editorState)) next.editorState = { ...document.editorState, featureCatalog: null };
  return next;
};

/**
 * The keys a Session 27-39 holds that no writer of its version wrote: a key
 * outside `SESSION_TOP_LEVEL_AUTHORITY`, and a Session 46 key (TB6, OV-590,
 * A7a, OV-638). Python reads none of them (`validate_session` accepts them).
 * @type {readonly string[]}
 */
const MODE_SCOPED_TOP_LEVEL_KEYS = Object.freeze(['modes', 'cliOptions']);

/**
 * `container` without `fields`, or `container` itself when it holds none.
 * @template T
 * @param {T} container
 * @param {(field: string) => boolean} drop
 */
const withoutFields = (container, drop) => (isPlainObject(container) && Object.keys(container).some(drop)
  ? /** @type {T} */ (Object.fromEntries(Object.entries(container).filter(([field]) => !drop(field)))) : container);

/**
 * A Result set holder's Web-only artifact fields (a Session 40+'s top level,
 * a Session 49 drawing): a development-branch feature field is dropped with
 * a row (J3), the derived orthogroup `groups` with none (A7c).
 * @template {Record<string, unknown>} T
 * @param {T} holder
 * @param {string} at Where the holder is in the file (`''`: the top level; else ending in `.`).
 * @param {LoadRepair[]} repairs Collects the rows.
 * @returns {T} `holder` itself when it holds none.
 */
const withoutArtifactWebOnlyFields = (holder, at, repairs) => {
  const features = isPlainObject(holder.features) ? holder.features : {};
  CURRENT_WRITER_FORBIDDEN_FEATURE_FIELDS.filter((field) => Object.hasOwn(features, field)).forEach((field) => repairs.push(
    loadRepair({ code: 'session-field-dropped', path: `${at}features.${field}`, effect: 'dropped' })
  ));
  const kept = {
    features: withoutFields(holder.features, (field) => CURRENT_WRITER_FORBIDDEN_FEATURE_FIELDS.includes(field)),
    orthogroupState: withoutFields(holder.orthogroupState, (field) => field === 'groups')
  };
  return kept.features === holder.features && kept.orthogroupState === holder.orthogroupState ? holder
    : /** @type {T} */ ({ ...holder, ...(Object.hasOwn(holder, 'features') ? { features: kept.features } : {}),
      ...(Object.hasOwn(holder, 'orthogroupState') ? { orthogroupState: kept.orthogroupState } : {}) });
};

/**
 * The Web-only part of the band reader: what the description leaves open
 * that no writer of the version wrote. A Session 27-39 key no writer of its
 * version wrote is dropped with a row; a Session 46's retired homes, whose
 * place its slices fill, and a Session 40+'s derived orthogroup `groups`
 * are dropped with none (A7f, A7c, BA-OD93 a); a development-branch feature
 * field is dropped with a row (J3).
 * @param {SessionChecked} document
 * @param {LoadRepair[]} repairs Collects the rows.
 * @returns {SessionChecked}
 */
const readWebOnlyFields = (document, repairs) => {
  /** @type {Record<string, unknown>} */
  let next = document;
  /** @param {string} key @param {unknown} value */
  const set = (key, value) => { if (next[key] !== value) next = { ...next, [key]: value }; };
  if (!isAuthorityBand(document)) {
    const unread = Object.keys(document).filter((key) => !Object.hasOwn(SESSION_TOP_LEVEL_AUTHORITY, key) || MODE_SCOPED_TOP_LEVEL_KEYS.includes(key));
    unread.forEach((path) => repairs.push(loadRepair({ code: 'session-field-dropped', path, effect: 'dropped' })));
    next = withoutFields(next, (key) => unread.includes(key));
  } else {
    next = withoutArtifactWebOnlyFields(next, '', repairs);
  }
  if (isModeScopedBand(document)) {
    RETIRED_MODE_SCOPED_FIELDS.forEach(([domain, fields]) => {
      if (!domain) return; // the top-level homes are not in the band's description
      const [key, child] = domain.split('.');
      const container = next[key];
      if (child === undefined) set(key, withoutFields(container, (field) => fields.has(field)));
      else if (isPlainObject(container)) set(key, { ...container, [child]: withoutFields(container[child], (field) => fields.has(field)) });
    });
  }
  return /** @type {SessionChecked} */ (next); // only fields the description leaves open were dropped
};

/**
 * One input inventory as the app reads it (A8, A6b): each binding that names
 * a resource the Session does not hold (the slot it filled is named), or that
 * cannot be read, is dropped with a row; an inventory that is no object or of
 * a schema no writer wrote is dropped whole.
 * @param {unknown} webFiles The `webFiles` that holds the inventory.
 * @param {unknown} resources The Session's resources, as read.
 * @param {'c_gb' | 'linearSeqs'} inputSlot The input slot schema 2 holds (`webBindingInputSlot`).
 * @param {string} at Where `webFiles` is in the file.
 * @param {LoadRepair[]} repairs Collects the rows.
 * @returns {unknown} `webFiles` itself when nothing is dropped.
 */
const readWebFileBindings = (webFiles, resources, inputSlot, at, repairs) => {
  const faults = webFileBindingFaults(/** @type {UnvalidatedFields<'bindings'>} */ (webFiles),
    /** @type {UnvalidatedFields<string, UnvalidatedFields<keyof RawSessionResourceDescriptor>>} */ (resources), inputSlot);
  if (!faults.length || !isPlainObject(webFiles)) return webFiles;
  if (faults.some(({ at: keys }) => keys === null)) {
    repairs.push(loadRepair({ code: 'session-field-dropped', path: `${at}.bindings`, effect: 'dropped' }));
    const { bindings: _dropped, ...kept } = webFiles;
    return kept; // the inventory is optional in every band
  }
  // `webFileBindingFaults` found an object inventory.
  const bindings = /** @type {Record<string, unknown>} */ (JSON.parse(JSON.stringify(webFiles.bindings)));
  /** @type {Map<unknown[], Set<number>>} */
  const removed = new Map();
  faults.forEach(({ at: place, field }) => {
    const keys = /** @type {ReadonlyArray<string | number>} */ (place); // no inventory fault is left
    repairs.push(loadRepair({ code: field ? 'session-field-dropped' : 'file-entry-dropped', path: `${at}.bindings.${webFileBindingPath(keys)}`, effect: 'dropped' }));
    const last = /** @type {string | number} */ (keys.at(-1));
    // `webFileBindingFaults` reached the binding through these keys: each one holds an object or a list.
    /** @type {unknown} */
    const parent = keys.slice(0, -1).reduce((container, key) => /** @type {Record<string | number, unknown>} */ (container)[key], /** @type {unknown} */ (bindings));
    if (Array.isArray(parent) && !field && !POSITIONAL_BINDING_LISTS.has(String(keys.findLast((key) => typeof key === 'string')))) {
      removed.set(parent, (removed.get(parent) ?? new Set()).add(Number(last)));
    } else if (field) delete /** @type {Record<string | number, unknown>} */ (parent)[last];
    else /** @type {Record<string | number, unknown>} */ (parent)[last] = null;
  });
  removed.forEach((indices, list) => {
    [...indices].sort((a, b) => b - a).forEach((index) => list.splice(index, 1));
  });
  return { ...webFiles, bindings }; // bindings read as above
};

/**
 * The input files of a Session as the app reads them (A8, A6b): a resource
 * whose descriptor cannot be read is dropped, with a row when a binding or
 * a request names it, and each input inventory is read
 * (`readWebFileBindings`): a Session 27-48's one, each Session 49 drawing's
 * own (B2d-3a-2). The request's references to a dropped resource are its
 * reader's (`readSessionRequests`).
 * @param {SessionChecked} document
 * @param {LoadRepair[]} repairs Collects the rows.
 * @returns {SessionChecked}
 */
const readSessionInputs = (document, repairs) => {
  /** @type {Record<string, unknown>} */
  let next = document;
  const resources = readSessionResources(next.resources);
  // A resource nothing names is dropped with no row (E-OD06, A15).
  // Each Result set's request names its resources (`storedResultSets`), and
  // so does each inventory.
  const entries = isDrawingsBand(document) && Array.isArray(next.drawings)
    ? /** @type {UnvalidatedFields<'renderRequest' | 'webFiles' | 'mode'>[]} */ (next.drawings) : [];
  const named = collectCanonicalResourceIds([
    ...storedResultSets(/** @type {SessionInput} */ (next)).map(({ set }) => set.renderRequest),
    ...[next, ...entries].map((holder) => /** @type {UnvalidatedFields<'bindings'>} */ (holder?.webFiles)?.bindings)]);
  resources.dropped.filter((id) => named.has(id.trim())).forEach((id) => repairs.push(
    loadRepair({ code: 'file-entry-dropped', path: `resources.${id}`, effect: 'dropped' })
  ));
  if (resources.value !== next.resources) next = { ...next, resources: resources.value };
  const webFiles = readWebFileBindings(next.webFiles, next.resources, 'c_gb', 'webFiles', repairs);
  if (webFiles !== next.webFiles) next = { ...next, webFiles };
  if (entries.length) {
    const read = entries.map((drawing, index) => {
      if (!isPlainObject(drawing) || !Object.hasOwn(drawing, 'webFiles')) return drawing;
      const drawingFiles = readWebFileBindings(drawing.webFiles, next.resources,
        webBindingInputSlot(true, drawing.mode), `${drawingSetPath(index)}webFiles`, repairs);
      return drawingFiles === drawing.webFiles ? drawing : { ...drawing, webFiles: drawingFiles };
    });
    if (read.some((drawing, index) => drawing !== entries[index])) next = { ...next, drawings: read };
  }
  return /** @type {SessionChecked} */ (next); // the inventories read as above
};

/**
 * The one reader of a parsed Session document (D1c 2.1): its envelope
 * (`readSessionInput`, R2 and R3), then its version band of
 * `SESSION_ADMISSION_DESCRIPTION`, repaired by `repairByDescription` with
 * `SESSION_PART_RULES`, and the band's Web-only rows. The value is copied
 * where repaired; the document itself is never changed.
 * @param {unknown} value
 * @returns {Repaired<SessionChecked>}
 */
export const readSessionDocument = (value) => {
  const envelope = readSessionInput(value);
  // `readSessionInput` admits no version below the first band's.
  const band = SESSION_ADMISSION_DESCRIPTION.bands.findLast(({ minVersion }) => envelope.version >= minVersion)
    ?? SESSION_ADMISSION_DESCRIPTION.bands[0];
  const read = repairByDescription(envelope, { type: band.root, path: '' }, { root: band.root, types: SESSION_READ_TYPES },
    SESSION_REPAIR_POLICY);
  /** @type {LoadRepair[]} */
  const repairs = read.rows.flatMap((row) => bandRepair(row, envelope, band.root) ?? []);
  // A root field the band requires that the file lacks leaves no document
  // (OV-1768): the file is no Session this app reads (R2).
  if (!isPlainObject(read.value)) throw diagnosticError('INPUT_INVALID', { field: 'schema', reason: 'SESSION_FORMAT' });
  const repaired = read.rows.some(({ path }) => path === 'results') ? withoutResultSetMembers(read.value) : read.value;
  const checked = /** @type {SessionChecked} */ (repaired); // the band's description holds, as repaired
  // A Session 49's drawings are read by drawing (Load's preflight); the
  // project holds the resources every drawing names. The description drops
  // a drawing without an ID, a name or a diagram mode; two drawings of one
  // ID leave no drawing to find by its ID: the file is no Session this app
  // reads (R2).
  if (isDrawingsBand(checked)) {
    if (!drawingIdsAreUnique(checked)) throw diagnosticError('INPUT_INVALID', { field: 'schema', reason: 'SESSION_FORMAT' });
    // Each drawing's Web-only artifact fields are read as a Session 40-48's top level's.
    const entries = /** @type {Record<string, unknown>[]} */ (/** @type {Record<string, unknown>} */ (checked).drawings);
    const read = entries.map((entry, index) => (isPlainObject(entry) ? withoutArtifactWebOnlyFields(entry, drawingSetPath(index), repairs) : entry));
    const document = read.some((entry, index) => entry !== entries[index])
      ? /** @type {SessionChecked} */ ({ ...checked, drawings: read }) : checked;
    return { value: readSessionInputs(document, repairs), repairs };
  }
  return {
    value: readSettingsOnlyParts(readSessionInputs(readWebOnlyFields(withNullRequestAsMissing(checked), repairs), repairs), repairs),
    repairs
  };
};

/**
 * Each drawing of a Session 49 has an ID of its own (`_identity_problem`,
 * gbdraw/session_drawings.py).
 * @param {SessionChecked49} document
 * @returns {boolean}
 */
const drawingIdsAreUnique = (document) => {
  const ids = document.drawings.map((entry) => entry.id);
  return new Set(ids).size === ids.length;
};

/**
 * R2-B1 (E-OD29): a null request marks a settings-only Session; a file whose
 * request is null but that holds saved Results or a feature list, which only
 * a Session with a request holds, is a full Session whose request is null.
 * Its request reads as a missing one, which has one reading: the request
 * reader's `request-dropped` row and the settings-only Session of that row
 * (K22), R4 when nothing is left. A new object.
 * @param {SessionChecked} document
 * @returns {SessionChecked}
 */
const withNullRequestAsMissing = (document) => {
  if (!isSettingsOnlySessionDocument(document)) return document;
  const editorState = /** @type {UnvalidatedFields<'featureCatalog'>} */ (document.editorState);
  const results = /** @type {unknown} */ (document.results);
  if (!(Array.isArray(results) && results.length > 0) && !isPlainObject(editorState?.featureCatalog)) return document;
  const { renderRequest: _null, ...missing } = document;
  return /** @type {SessionChecked} */ (missing); // the request is optional in the band of a settings-only form
};

/**
 * R4 (receipt 〔拒否〕): Load restores nothing of a file whose every content
 * part it holds was dropped or replaced by its default, and of a
 * settings-only file, whose content is its settings, when its settings hold
 * nothing or were dropped (D1a: `modes` null or text). A file that holds no
 * part (an envelope-only Session 27-30) loads as before (BA-D20). Each
 * producer that drops a part records the part's key (or `key.member` for
 * each member it drops) as its path; a part dropped with no row still
 * counts (E-OD21): the `request-dropped` row of a
 * request Load cannot read (K22) loads the file as settings only, which drops
 * every other part, so it is judged as a settings-only file (OV-864); a
 * Session of a version with no settings-only form (31-41) has nothing to load
 * without its request, so that row alone restores nothing (C-OD17, BA-OD63 d).
 * More repairs never turn a verdict of true to false, so Load may ask before
 * every producer ran, once a row has decided it.
 * @param {SessionEnvelope} stored The document as the file holds it.
 * @param {readonly LoadRepair[]} repairs The repairs of the Load so far.
 */
export const restoresNothing = (stored, repairs) => {
  const requestDropped = repairs.some(({ code }) => code === 'request-dropped');
  if (requestDropped && !hasSettingsOnlyForm(stored)) return true;
  // A Session 46-48's stored set dropped for its request is its top-level
  // part (a Session 49 drawing stays without its set).
  const dropped = new Set(repairs.filter(({ effect }) => effect === 'dropped' || effect === 'defaulted')
    .map(({ code, path }) => (code === 'other-result-set-dropped' ? topLevelKey(path) : path)));
  // A Session 49 restores nothing when every drawing it holds was dropped.
  const storedDrawings = /** @type {{ drawings?: unknown }} */ (/** @type {unknown} */ (stored)).drawings;
  if (stored.version >= SESSION_DRAWINGS_LAYOUT.version && Array.isArray(storedDrawings)) {
    return storedDrawings.every((_entry, index) => dropped.has(drawingSetPath(index).slice(0, -1)));
  }
  /**
   * What a part holds: the members of an object part that hold something
   * (`key.member`; a part whose members hold nothing holds nothing), the
   * value of another part (`key`).
   * @param {string} key
   */
  const members = (key) => {
    const value = stored[/** @type {keyof SessionInput} */ (key)];
    if (!isPlainObject(value)) return holdsNothing(value) ? [] : [key];
    return Object.entries(/** @type {Record<string, unknown>} */ (value))
      .filter(([, member]) => !holdsNothing(member)).map(([member]) => `${key}.${member}`);
  };
  // A part is dropped whole, or member by member (both slices of `modes`, OV-865).
  /** @param {string} key */
  const lost = (key) => dropped.has(key) || members(key).every((path) => dropped.has(path));
  /** @param {SessionContentPart | null} [only] */
  const held = (only) => Object.entries(SESSION_TOP_LEVEL_AUTHORITY)
    .filter(([key, { part }]) => part !== null && (only === undefined || part === only) && members(key).length > 0)
    .map(([key]) => key);
  if (isSettingsOnlySessionDocument(stored) || requestDropped) {
    return held('settings').every(lost);
  }
  const parts = held();
  return parts.length > 0 && parts.every(lost);
};

const adoptedCurrentDocuments = new WeakSet();
const adoptedCanonicalOwners = new WeakSet();

// Membership is tested on unvalidated Session values, hence `ReadonlySet<unknown>`.
/** @type {ReadonlySet<unknown>} */
const LINEAR_COMPARISON_PLAN_MODES = new Set(['none', 'adjacent', 'selected']);
/** @type {ReadonlySet<unknown>} */
const LINEAR_COMPARISON_SOURCES = new Set(['losat', 'upload']);
const LINEAR_COMPARISON_EDGE_FIELDS = new Set([
  'id',
  'queryUid',
  'subjectUid',
  'included',
  'fileActive',
  'losatFilenameActive',
  'source',
  'losatFilename'
]);

/**
 * @param {unknown} value
 * @param {string} field
 * @param {string} message
 * @returns {void}
 */
const assertNoOwnField = (value, field, message) => {
  if (isPlainObject(value) && Object.prototype.hasOwnProperty.call(value, field)) {
    throw new Error(message);
  }
};

/**
 * @param {unknown} plan
 * @returns {Set<string>} The edge IDs.
 */
const validateLinearComparisonPlan = (plan) => {
  if (!isPlainObject(plan)) {
    throw new Error('Current Web comparison draft requires config.linearComparisonPlan.');
  }
  if (!LINEAR_COMPARISON_PLAN_MODES.has(plan.mode)) {
    throw new Error('config.linearComparisonPlan.mode is invalid.');
  }
  if (!LINEAR_COMPARISON_SOURCES.has(plan.defaultSource)) {
    throw new Error('config.linearComparisonPlan.defaultSource is invalid.');
  }
  if (!Array.isArray(plan.edges)) {
    throw new Error('config.linearComparisonPlan.edges must be an array.');
  }
  /** @type {Set<string>} */
  const ids = new Set();
  /** @type {unknown[]} */ (plan.edges).forEach((edge) => {
    if (!isPlainObject(edge)) {
      throw new Error('Each config.linearComparisonPlan edge must be an object.');
    }
    const unknown = Object.keys(edge).filter((field) => !LINEAR_COMPARISON_EDGE_FIELDS.has(field));
    if (unknown.length > 0) {
      throw new Error(
        `config.linearComparisonPlan edge contains retired or unknown field(s): ${unknown.join(', ')}`
      );
    }
    const id = typeof edge.id === 'string' ? edge.id.trim() : '';
    const queryUid = typeof edge.queryUid === 'string' ? edge.queryUid.trim() : '';
    const subjectUid = typeof edge.subjectUid === 'string' ? edge.subjectUid.trim() : '';
    if (!id || !queryUid || !subjectUid) {
      throw new Error('Current comparison-plan edges require stable IDs and endpoint UIDs.');
    }
    if (ids.has(id)) {
      throw new Error(`Current comparison-plan edge ID is duplicated: ${id}.`);
    }
    ids.add(id);
    if (
      typeof edge.included !== 'boolean'
      || typeof edge.fileActive !== 'boolean'
      || typeof edge.losatFilenameActive !== 'boolean'
      || !LINEAR_COMPARISON_SOURCES.has(edge.source)
      || typeof edge.losatFilename !== 'string'
    ) {
      throw new Error('Current comparison-plan edge metadata is invalid.');
    }
  });
  return ids;
};

/**
 * @param {SessionInput} sessionData Unvalidated Session data.
 * @param {{ load?: boolean }} [options] Load reads the plan (its readers list
 *   what they change, OV-764); Save requires it as written.
 * @returns {void}
 */
export const validateCurrentComparisonAuthority = (sessionData, { load = false } = {}) => {
  // Session 46: the comparison draft is the Linear drawing's.
  const draftConfig = Number(sessionData.version) >= MODE_SCOPED_SESSION_VERSION
    ? /** @type {UnvalidatedFields<'linear', UnvalidatedFields<'config'>>} */ (sessionData.modes)?.linear?.config
    : sessionData.config;
  /** @type {Record<string, unknown>} */
  const config = isPlainObject(draftConfig) ? draftConfig : {};
  /** @type {Record<string, unknown>} */
  const ui = isPlainObject(sessionData.ui) ? sessionData.ui : {};
  /** @type {Record<string, unknown>} */
  const webFiles = isPlainObject(sessionData.webFiles) ? sessionData.webFiles : {};
  const bindings = /** @type {NonNullable<UnvalidatedFields<'linearSeqs' | 'linearComparisons'>>} */ (
    webFiles.bindings || {}
  );

  assertNoOwnField(config, 'blastSource', 'Current sessions cannot contain config.blastSource.');
  assertNoOwnField(config.adv, 'blastSource', 'Current sessions cannot contain config.adv.blastSource.');
  assertNoOwnField(ui, 'blastSource', 'Current sessions cannot contain ui.blastSource.');
  assertNoOwnField(
    config.linearRecordLayout,
    'comparisons',
    'Current sessions cannot contain config.linearRecordLayout.comparisons.'
  );

  (Array.isArray(bindings.linearSeqs) ? bindings.linearSeqs : []).forEach((/** @type {unknown} */ sequence) => {
    assertNoOwnField(
      sequence,
      'blast',
      'Current sessions cannot contain per-record BLAST file bindings.'
    );
    assertNoOwnField(
      sequence,
      'losat_filename',
      'Current sessions cannot contain per-record LOSAT filenames.'
    );
  });
  (Array.isArray(webFiles.linearRecordMetadata) ? webFiles.linearRecordMetadata : [])
    .forEach((/** @type {unknown} */ metadata) => {
      assertNoOwnField(
        metadata,
        'losatFilename',
        'Current sessions cannot contain per-record LOSAT filename metadata.'
      );
      assertNoOwnField(
        metadata,
        'losat_filename',
        'Current sessions cannot contain per-record LOSAT filename metadata.'
      );
    });
  assertNoOwnField(
    bindings,
    'linearCanonicalComparisons',
    'Current sessions cannot bind comparison artifacts outside the committed request.'
  );

  const hasWebComparisonDraft = Object.prototype.hasOwnProperty.call(
    config,
    'linearRecordLayout'
  ) || Object.prototype.hasOwnProperty.call(config, 'linearComparisonPlan');
  if (
    Object.prototype.hasOwnProperty.call(bindings, 'linearComparisons')
    && !Array.isArray(bindings.linearComparisons)
  ) {
    throw new Error('Current comparison file bindings must be an array.');
  }
  /** @type {unknown[]} */
  const bindingEntries = Array.isArray(bindings.linearComparisons) ? bindings.linearComparisons : [];
  if (!hasWebComparisonDraft && bindingEntries.length === 0) return;

  // Load reads the plan as the drawing reads it; Save requires it as written.
  const edges = load ? normalizeLinearComparisonPlan(config.linearComparisonPlan).edges : null;
  const edgeIds = edges ? new Set(edges.map(({ id }) => id)) : validateLinearComparisonPlan(config.linearComparisonPlan);
  /** @type {Set<string>} */
  const boundIds = new Set();
  bindingEntries.forEach((binding) => {
    if (!isPlainObject(binding)) {
      throw new Error('Each current comparison file binding must be an object.');
    }
    const unknown = Object.keys(binding).filter((field) => !['id', 'file'].includes(field));
    if (unknown.length > 0) {
      throw new Error(
        `Current comparison file binding duplicates plan metadata: ${unknown.join(', ')}`
      );
    }
    const id = typeof binding.id === 'string' ? binding.id.trim() : '';
    if (!id || !edgeIds.has(id)) {
      throw new Error('Current comparison file binding must reference a plan edge ID.');
    }
    if (boundIds.has(id)) {
      throw new Error(`Current comparison file binding is duplicated: ${id}.`);
    }
    // Load drops a file binding it cannot read (`readSessionInputs`) and the
    // Files apply leaves that edge with no file and none active, as the pair
    // card's Remove does (`setComparisonEdgeFile`, OV-899).
    if (!isPlainObject(binding.file) && !(load && binding.file === null)) {
      throw new Error('Each current comparison file binding requires a file resource binding.');
    }
    boundIds.add(id);
  });
  // `validateLinearComparisonPlan` above proved the plan's edges: objects with a string `id` and a boolean `fileActive`.
  (edges ?? /** @type {{ edges: { id: string, fileActive: boolean }[] }} */ (config.linearComparisonPlan).edges).forEach((edge) => {
    if (edge.fileActive && !boundIds.has(edge.id)) {
      throw new Error(
        `Active comparison file is missing its Web file binding: ${edge.id}.`
      );
    }
  });
};

/**
 * @param {unknown} source
 * @param {readonly string[]} fields
 * @returns {Record<string, unknown>}
 */
const copyFields = (source, fields) => {
  /** @type {Record<string, unknown>} */
  const projected = {};
  if (!source || typeof source !== 'object' || Array.isArray(source)) return projected;
  fields.forEach((field) => {
    if (Object.prototype.hasOwnProperty.call(source, field)) {
      // `source` is an object that is no array (checked above) and holds `field`.
      projected[field] = /** @type {Record<string, unknown>} */ (source)[field];
    }
  });
  return projected;
};

/**
 * @param {unknown} value
 * @returns {boolean}
 */
const hasInput = value => Array.isArray(value) ? value.some(hasInput) : value != null;

// The same complete inventory is used before Save and at document admission.
// An inactive input remains biological even when it cannot render the active mode.
/**
 * @param {Record<string, any>} [files] The Session's Web input inventory.
 * @returns {boolean}
 */
export const hasBiologicalSessionInputs = (files = {}) => (
  ['c_gb', 'c_gff', 'c_fasta', 'c_conservation_fastas', 'c_conservation_sequence_sources']
    .some(key => hasInput(files[key]))
  || (files.linearSeqs || []).some((/** @type {Record<string, unknown>} */ row) => ['gb', 'gff', 'fasta'].some(key => hasInput(row[key])))
);

/**
 * @param {SessionInput | null | undefined} data
 * @returns {boolean}
 */
export const isSettingsOnlySessionDocument = data => data != null && typeof data.version === 'number'
  && [42, 44, MODE_SCOPED_SESSION_VERSION, ROW_KEYED_LEGEND_SESSION_VERSION, SESSION_VERSION].includes(data.version)
  && Object.hasOwn(data, 'renderRequest') && data.renderRequest === null;

/**
 * Whether a Session's version has a settings-only form (42, 44, 46, 47, 48): one
 * whose request cannot be read loads as settings only (K22, C-OD07).
 * @param {SessionInput | null | undefined} data
 * @returns {boolean}
 */
export const hasSettingsOnlyForm = (data) => data != null && isSettingsOnlySessionDocument({ ...data, renderRequest: null });

/**
 * The request schemas a file stored, read before Load promoted its requests:
 * a check that depends on what the writer wrote reads them (RK-1).
 * Keyed by set path (`ResultSetPath`, `storedResultSets`): the projected
 * set's under `''`, where its view holds it (a Session 49's opened drawing
 * too), every stored set's under its own path.
 * @typedef {{ [at in ResultSetPath]?: unknown }} StoredRequestSchemas
 */

/** The Web input inventory of a Session that binds no input. */
const EMPTY_WEB_BINDINGS = Object.freeze({
  schema: 2, c_gb: null, c_gff: null, c_fasta: null, c_depth: null, c_conservation_blasts: [],
  c_conservation_blasts_source: null, c_conservation_fastas: [], c_conservation_sequence_sources: [], d_color: null,
  t_color: null, blacklist: null, whitelist: null, qualifier_priority: null, linearSeqs: [], linearComparisons: []
});

/**
 * The settings-only Session of a Session whose request cannot be read (K22):
 * the request, the Results, the GenBank, GFF and FASTA inputs (and the
 * conservation sequences), the feature list and the caches are dropped; the
 * settings, the other bindings and the resources they name are kept. A new
 * object; the parsed file is not changed.
 * @param {SessionEnvelope} session
 * @returns {SessionEnvelope}
 */
export const settingsOnlySessionCandidate = (session) => {
  let candidate = JSON.parse(JSON.stringify(session));
  // A Session 42-48 (the versions with a settings-only form) keeps only its projected set.
  for (const { at } of storedResultSets(candidate).filter(({ projected }) => !projected)) {
    candidate = withStoredResultSet(candidate, at, null);
  }
  const { cliInvocation: _cli, legacyArtifacts: _legacy, files: _files, ...kept } = candidate;
  const editorState = { ...(isPlainObject(kept.editorState) ? kept.editorState : {}), featureCatalog: null };
  delete editorState.alignmentResetReceipt;
  /** @param {unknown} cache */
  const emptyCache = (cache) => (isPlainObject(cache) ? { ...cache, entries: [] } : cache);
  const stored = kept.webFiles?.bindings;
  const bindings = { ...(isPlainObject(stored) && stored.schema === 2 ? stored : JSON.parse(JSON.stringify(EMPTY_WEB_BINDINGS))) };
  Object.assign(bindings, { c_gb: null, c_gff: null, c_fasta: null, c_conservation_fastas: [], c_conservation_sequence_sources: [] });
  if (bindings.c_conservation_blasts_source === 'losat-cache') {
    Object.assign(bindings, { c_conservation_blasts_source: null, c_conservation_blasts: [] });
  }
  bindings.linearSeqs = (Array.isArray(bindings.linearSeqs) ? bindings.linearSeqs : [])
    .map((/** @type {unknown} */ row) => (isPlainObject(row) ? { ...row, gb: null, gff: null, fasta: null } : row));
  const linearOff = kept.modes?.linear?.config?.recordsOff;
  if (Array.isArray(linearOff)) {
    const uids = new Set(bindings.linearSeqs.map((/** @type {Record<string, unknown>} */ row) => row?.uid));
    kept.modes.linear.config.recordsOff = linearOff.filter((uid) => uids.has(uid));
  }
  /** @type {Set<unknown>} */
  const referenced = new Set();
  /** @param {unknown} value */
  const visit = (value) => {
    if (!value || typeof value !== 'object') return;
    if (Object.hasOwn(value, 'resourceId')) referenced.add(/** @type {{ resourceId: unknown }} */ (value).resourceId);
    Object.values(value).forEach(visit);
  };
  visit(bindings);
  return /** @type {SessionEnvelope} */ ({
    ...kept,
    renderRequest: null,
    results: [],
    editorState,
    runMetadata: {},
    ...(Object.hasOwn(kept, 'losatCache') ? { losatCache: emptyCache(kept.losatCache) } : {}),
    ...(Object.hasOwn(kept, 'losatDerivedCache') ? { losatDerivedCache: emptyCache(kept.losatDerivedCache) } : {}),
    ...(Object.hasOwn(kept, 'proteinIdentityManifest')
      ? { proteinIdentityManifest: { schema: 2, proteinSets: {}, recordAnalyses: {}, recordInstances: {} } } : {}),
    webFiles: { bindings },
    resources: Object.fromEntries(Object.entries(isPlainObject(kept.resources) ? kept.resources : {})
      .filter(([id]) => referenced.has(id)))
  });
};

/** The GenBank, GFF and FASTA bindings of a Web input inventory (`hasBiologicalSessionInputs`). */
const BIOLOGICAL_BINDINGS = Object.freeze(['c_gb', 'c_gff', 'c_fasta', 'c_conservation_fastas', 'c_conservation_sequence_sources']);
/**
 * What a settings-only Session holds that one may not: Save refuses each
 * (`validateSettingsOnlyDocument`), Load drops each before the dialog
 * (`readSettingsOnlyParts`; A15, E-OD06).
 * @param {SessionEnvelope} data
 */
const settingsOnlyFaults = (data) => {
  // Read as dev's Save read them (E-OD02): a value of another kind with a
  // length or keys counts as held, as Python's `validate_session` refuses it.
  /** @param {unknown} value */
  const holdsKeys = (value) => Object.keys(value || {}).length > 0;
  /** @param {unknown} cache */
  const holdsCacheEntries = (cache) => (/** @type {{ entries?: { length?: unknown } } | null | undefined} */ (cache)?.entries || []).length;
  const bindings = /** @type {UnvalidatedFields<'bindings'>} */ (data.webFiles)?.bindings;
  const inventory = isPlainObject(data.resources) && isPlainObject(bindings) && Array.isArray(bindings.linearSeqs)
    && ['c_gb', 'c_gff', 'c_fasta'].every(key => Object.hasOwn(bindings, key));
  // Every Linear row is an object, as Python's reader requires
  // (`_validate_settings_only_session`); `validateWebFileBindings` lets a
  // null row through (OV-566).
  const rows = inventory && Array.isArray(bindings.linearSeqs) && bindings.linearSeqs.every(isPlainObject);
  /** @type {Set<unknown>} */
  const referenced = new Set();
  /** @param {unknown} value */
  const visit = value => {
    if (!value || typeof value !== 'object') return;
    if (Object.hasOwn(value, 'resourceId')) referenced.add(/** @type {NonNullable<UnvalidatedFields<'resourceId'>>} */ (value).resourceId);
    Object.values(value).forEach(visit);
  };
  visit(bindings);
  return {
    // `inventory` proved `bindings` an object that holds the explicit inventory.
    inventory,
    rows,
    biological: rows && hasBiologicalSessionInputs(/** @type {Record<string, unknown>} */ (bindings)),
    results: /** @type {UnvalidatedFields<'length'>} */ (data.results)?.length !== 0,
    artifacts: /** @type {UnvalidatedFields<'featureCatalog'>} */ (data.editorState)?.featureCatalog !== null
      || data.cliInvocation != null || holdsKeys(data.runMetadata)
      || Boolean(holdsCacheEntries(data.losatCache) || holdsCacheEntries(data.losatDerivedCache))
      || holdsKeys(data.legacyArtifacts)
      || ['proteinSets', 'recordAnalyses', 'recordInstances']
        .some(key => Object.keys(/** @type {UnvalidatedFields<string>} */ (data.proteinIdentityManifest)?.[key] || {}).length > 0)
      || /** @type {UnvalidatedFields<'c_conservation_blasts_source'>} */ (bindings)?.c_conservation_blasts_source === 'losat-cache',
    unboundResources: Object.keys(isPlainObject(data.resources) ? data.resources : {}).filter(id => !referenced.has(id))
  };
};

/**
 * A settings-only Session's Linear rows that are no object, dropped with a
 * row: Python's reader requires each row to be an object (OV-566, A6b).
 * @param {SessionChecked} document
 * @param {LoadRepair[]} repairs Collects the rows.
 * @returns {SessionChecked}
 */
const withoutLinearRowsOfAnotherKind = (document, repairs) => {
  const webFiles = /** @type {UnvalidatedFields<'bindings'>} */ (document.webFiles);
  const bindings = /** @type {UnvalidatedFields<'linearSeqs'>} */ (webFiles?.bindings);
  const rows = bindings?.linearSeqs;
  if (!Array.isArray(rows) || rows.every(isPlainObject)) return document;
  rows.forEach((row, index) => {
    if (!isPlainObject(row)) repairs.push(loadRepair({ code: 'session-field-dropped', path: `webFiles.bindings.linearSeqs[${index}]`, effect: 'dropped' }));
  });
  // `rows` was read from `webFiles.bindings`, so both are objects.
  return /** @type {SessionChecked} */ ({ ...document, webFiles: { ...webFiles, bindings: { ...bindings, linearSeqs: rows.filter(isPlainObject) } } });
};

/**
 * Load: a settings-only Session's Results and GenBank, GFF and FASTA inputs
 * are dropped with a row; its other committed artifacts and the resources no
 * binding names are dropped with none, as the settings-only Session of an
 * unreadable request drops them (`settingsOnlySessionCandidate`; A15, J4).
 * A Web input inventory the file does not hold is the empty one, with no row
 * (E-OD06; OV-860: a full Session whose request is null holds none).
 * @param {SessionChecked} stored
 * @param {LoadRepair[]} repairs Collects the rows.
 * @returns {SessionChecked}
 */
const readSettingsOnlyParts = (stored, repairs) => {
  if (!isSettingsOnlySessionDocument(stored)) return stored;
  const document = withoutLinearRowsOfAnotherKind(stored, repairs);
  const faults = settingsOnlyFaults(document);
  if (faults.rows && !(faults.biological || faults.results || faults.artifacts || faults.unboundResources.length)) return document;
  if (faults.results) repairs.push(loadRepair({ code: 'session-part-dropped', path: 'results', effect: 'dropped' }));
  const bindings = /** @type {UnvalidatedFields<'bindings'>} */ (document.webFiles)?.bindings;
  if (isPlainObject(bindings)) {
    BIOLOGICAL_BINDINGS.filter((key) => hasBiologicalSessionInputs({ [key]: bindings[key] })).forEach((key) => repairs.push(
      loadRepair({ code: 'session-field-dropped', path: `webFiles.bindings.${key}`, effect: 'dropped' })
    ));
    (Array.isArray(bindings.linearSeqs) ? bindings.linearSeqs : []).forEach((row, index) => {
      if (hasBiologicalSessionInputs({ linearSeqs: [row] })) {
        repairs.push(loadRepair({ code: 'session-field-dropped', path: `webFiles.bindings.linearSeqs[${index}]`, effect: 'changed' }));
      }
    });
  }
  return /** @type {SessionChecked} */ (settingsOnlySessionCandidate(document)); // the band's shape, with parts emptied
};

/**
 * @param {SessionEnvelope} data
 * @param {{ load?: boolean }} [options] Load drops the parts `readSettingsOnlyParts` drops, and reads a missing settings container as its default.
 */
const validateSettingsOnlyDocument = (data, { load = false } = {}) => {
  const faults = settingsOnlyFaults(data);
  if (!faults.inventory) {
    throw new Error('Settings-only Session requires resources and an explicit Web input inventory.');
  }
  if (!faults.rows || faults.biological
    // `data.webFiles` holds the `bindings` object read above, so it is neither null nor undefined.
    || Object.keys(/** @type {NonNullable<UnvalidatedFields<'bindings'>>} */ (data.webFiles)).some(key => key !== 'bindings')) {
    throw new Error('Settings-only Session cannot contain biological sources.');
  }
  if (faults.results || faults.artifacts) {
    throw new Error('Settings-only Session cannot contain committed render artifacts.');
  }
  if (load) return;
  // Session 46: each mode's slice holds that mode's settings draft, with its
  // `form` and `adv` (`_validate_settings_only_session`, OV-591).
  if (data.version >= MODE_SCOPED_SESSION_VERSION) {
    // `ui.mode` is a Session value: the slice is read by its own key (OV-612).
    const mode = /** @type {UnvalidatedFields<'mode'>} */ (data.ui)?.mode;
    const slice = isPlainObject(data.modes) && typeof mode === 'string' && Object.hasOwn(data.modes, mode)
      ? data.modes[mode]
      : undefined;
    const config = /** @type {UnvalidatedFields<'config'>} */ (slice)?.config;
    if (!isPlainObject(config) || !isPlainObject(config.form) || !isPlainObject(config.adv)) {
      throw new Error('Settings-only Session requires the settings of its mode.');
    }
  } else {
    validateSettingsOnlyFlatDraft(data);
  }
  if (faults.unboundResources.length > 0) {
    throw new Error('Settings-only Session contains an unbound resource.');
  }
};

// A settings-only Session 42 or 44 keeps its flat draft in `config`.
/** @param {SessionEnvelope} data */
const validateSettingsOnlyFlatDraft = (data) => {
  let storedConfig = data.version >= TYPED_DRAFT_SESSION_VERSION
    || !isPlainObject(/** @type {UnvalidatedFields<'adv'>} */ (data.config)?.adv)
    ? data.config
    : {
        // This branch runs only when `data.config` holds an `adv` object (the condition above).
        .../** @type {{ adv: Record<string, unknown> }} */ (data.config),
        adv: migrateLegacyLinearLabelVisibility(/** @type {{ adv: Record<string, unknown> }} */ (data.config).adv)
      };
  if (data.version < TYPED_DRAFT_SESSION_VERSION && isPlainObject(storedConfig)
    && Object.prototype.hasOwnProperty.call(storedConfig, 'recordDisplayDrafts')) {
    storedConfig = {
      ...storedConfig,
      recordDisplayDrafts: migrateLegacyRecordDisplayDrafts(
        storedConfig.recordDisplayDrafts
      )
    };
  }
  if (isPlainObject(storedConfig)
    && Object.prototype.hasOwnProperty.call(storedConfig, 'featurePlacementOverrides')) {
    storedConfig = {
      ...storedConfig,
      featurePlacementOverrides: migrateSessionFeaturePlacements(storedConfig.featurePlacementOverrides)
    };
  }
  // The check of `mode` and of the configuration's containers; Load reads
  // the values inside them with the reader table (`configFieldRepairs`).
  requireActiveConfigContainers({ mode: /** @type {UnvalidatedFields<'mode'>} */ (data.ui)?.mode, storedConfig });
};

// One committed Result set: the top-level fields, or `otherModeResult`, whose
// fields mirror them.
/**
 * One thing a committed Result set holds that its Results cannot be restored
 * with: Save refuses the first (`validateCommittedResultSet`), Load repairs
 * each (`readCommittedResultSet`).
 * - `results`: no list, or a Result that is no named plain SVG (A7b, A20);
 * - `editorState`: no editor state that holds `featureCatalog`;
 * - `receipt-missing`: no alignment reset receipt where the request needs one;
 * - `receipt`: a receipt that is not of the request's plan (A19);
 * - `catalog`: a feature catalog of another schema; `catalog-null`: none for
 *   saved Results (A9).
 * @typedef {{ kind: 'results' | 'editorState' | 'receipt-missing' | 'receipt' | 'catalog' | 'catalog-null' }} ResultSetFault
 */
/** @param {unknown} result */
const isPlainSvgResult = (result) => {
  const resultName = isPlainObject(result) && typeof result.name === 'string' ? result.name.trim() : '';
  const content = isPlainObject(result) ? result.content : null;
  return Boolean(resultName) && !resultName.toLowerCase().endsWith('.interactive.svg')
    && typeof content === 'string' && content.includes('<svg')
    && !content.includes('gbdraw-interactive-feature-metadata') && !content.includes('gbdraw-interactive-feature-script');
};

/** @param {SessionInput} set */
const storedRequestSchema = (set) => /** @type {UnvalidatedFields<'schema'>} */ (set.renderRequest)?.schema;

/**
 * @param {SessionInput} set
 * @param {number | string} version
 * @param {unknown} [storedSchema] The schema the file stored, before Load read the request (RK-1).
 * @returns {ResultSetFault[]}
 */
const committedResultSetFaults = (set, version, storedSchema = storedRequestSchema(set)) => {
  /** @type {ResultSetFault[]} */
  const faults = [];
  const results = set.results;
  if (!Array.isArray(results) || !results.every(isPlainSvgResult)) faults.push({ kind: 'results' });
  const editorState = isPlainObject(set.editorState) ? set.editorState : {};
  if (!Object.hasOwn(editorState, 'featureCatalog')) faults.push({ kind: 'editorState' });
  // Load read the request before (`readSessionRequests`); Save holds a canonical one.
  const request = /** @type {AlignmentResetRequest | null | undefined} */ (set.renderRequest);
  if (Number(storedSchema) >= 8 && request?.layout?.similarityAlignment && !Object.hasOwn(editorState, 'alignmentResetReceipt')) {
    faults.push({ kind: 'receipt-missing' });
  } else if (!alignmentResetReceiptShapeHolds(editorState.alignmentResetReceipt, request)) {
    faults.push({ kind: 'receipt' });
  }
  const featureCatalog = Object.hasOwn(editorState, 'featureCatalog') ? editorState.featureCatalog : null;
  const expectedCatalogSchema = FEATURE_CATALOG_SCHEMA_BY_SESSION_VERSION[Number(version)] ?? 3;
  if (featureCatalog !== null && (!isPlainObject(featureCatalog) || featureCatalog.schema !== expectedCatalogSchema)) {
    faults.push({ kind: 'catalog' });
  } else if (Array.isArray(results) && results.length > 0 && featureCatalog === null && faults[0]?.kind !== 'results') {
    faults.push({ kind: 'catalog-null' });
  }
  return faults;
};

/**
 * Save: refuses the first fault of a committed Result set.
 * @param {readonly { kind: string }[]} faults
 */
const refuseFirstFault = ([fault]) => {
  if (!fault) return;
  // A malformed or stale receipt (OV-115, OV-130) is one user-facing failure.
  if (fault.kind === 'receipt') throw diagnosticError('ALIGNMENT_RESET_EVIDENCE');
  throw diagnosticError('INPUT_INVALID', { field: 'schema', reason: 'FIELDS' });
};

/**
 * The field of a Session 49 that holds its drawings (a file path, not the
 * app's drawing list).
 */
const DRAWING_LIST_FIELD = /** @type {const} */ ('drawings');
/**
 * Where a Session 49 holds a drawing's Result set: the prefix of its paths,
 * written as the description reader writes a list member's (`list[i].`).
 * @typedef {`${typeof DRAWING_LIST_FIELD}[${number}].`} DrawingSetPath
 */
/**
 * Where a committed Result set is in a Session: the top level or the other
 * mode's set of a Session 27-48, or a Session 49 drawing's (the prefix of its
 * paths).
 * @typedef {'' | 'otherModeResult.' | DrawingSetPath} ResultSetPath
 */
/** @param {number} index @returns {DrawingSetPath} */
export const drawingSetPath = (index) => `${DRAWING_LIST_FIELD}[${index}].`;
/** The index of the Session 49 drawing whose set is at `at`, else null. @param {ResultSetPath} at */
export const drawingSetIndex = (at) => {
  const matched = /^drawings\[(\d+)\]\.$/.exec(at);
  return matched ? Number(matched[1]) : null;
};
/**
 * Load: a committed Result set as far as it can be restored (A7b, A9, A19,
 * A20). Results that cannot be read are dropped with the members that name
 * them (positional: rule 6); an absent editor state reads as an empty one; a
 * stale receipt reads as none, with a row; a catalog that cannot be read is
 * dropped, and the Results stay without one (`feature-recovery-failed`).
 * @param {SessionInput} set
 * @param {readonly ResultSetFault[]} faults Of `set` (`committedResultSetFaults`).
 * @param {ResultSetPath} at
 * @param {LoadRepair[]} repairs Collects the rows.
 * @returns {SessionInput}
 */
const readCommittedResultSet = (set, faults, at, repairs) => {
  if (!faults.length) return set;
  const kinds = new Set(faults.map(({ kind }) => kind));
  const editorState = { ...(isPlainObject(set.editorState) ? set.editorState : {}) };
  if (kinds.has('editorState')) {
    if (!holdsNothing(set.editorState) && !isPlainObject(set.editorState)) {
      repairs.push(loadRepair({ code: 'session-part-dropped', path: `${at}editorState`, effect: 'defaulted' }));
    }
    editorState.featureCatalog = null;
  }
  // A receipt the file does not hold reads as none (E-OD06).
  if (kinds.has('receipt-missing')) editorState.alignmentResetReceipt = null;
  if (kinds.has('receipt')) {
    repairs.push(loadRepair({ code: 'session-field-dropped', path: `${at}editorState.alignmentResetReceipt`, effect: 'dropped' }));
    editorState.alignmentResetReceipt = null;
  }
  if (kinds.has('catalog') || kinds.has('catalog-null')) {
    repairs.push(loadRepair({ code: 'feature-recovery-failed', path: `${at}editorState.featureCatalog`, effect: 'dropped' }));
    editorState.featureCatalog = null;
  }
  if (!kinds.has('results')) return { ...set, editorState };
  repairs.push(loadRepair({ code: 'session-part-dropped', path: `${at}results`, effect: 'dropped' }));
  return /** @type {SessionInput} */ (withoutResultSetMembers({ ...set, results: [], editorState }));
};

/**
 * @param {SessionInput} set
 * @param {number | string} version
 * @param {unknown} [storedSchema] The schema the file stored, before Load read the request (RK-1).
 */
const validateCommittedResultSet = (set, version, storedSchema = storedRequestSchema(set)) => {
  refuseFirstFault(committedResultSetFaults(set, version, storedSchema));
};

// E1: Session 46 keeps the other diagram mode's Result set beside the
// top-level set. It needs a committed top-level request of the other mode, at
// least one Result, and resources in the top-level table.
/**
 * What a stored Result set of a Session 46-48 (`storedResultSets`, the other
 * mode's) holds that it cannot be restored with: `set` and `request` faults
 * come first; `set` faults of its committed Result set follow
 * (`committedResultSetFaults`). Run notices are read by
 * `readSessionAuthorityParts`.
 * @param {StoredResultSet} stored A stored set (an object: Load's band reader dropped one that is not, `session-part-dropped`).
 * @param {SessionInput} sessionData The Session that holds it: its projected request and its resources.
 * @param {number | string} version
 * @param {StoredRequestSchemas} storedSchemas
 * @returns {{ kind: 'set' | 'request' | ResultSetFault['kind'] }[]}
 */
const storedSetFaults = ({ at, set: other }, sessionData, version, storedSchemas) => {
  if (Object.keys(other).some((key) => !OTHER_MODE_RESULT_FIELDS.includes(key))) {
    return [{ kind: 'set' }];
  }
  const mode = /** @type {UnvalidatedFields<'mode'>} */ (other.renderRequest)?.mode;
  if (!isPlainObject(sessionData.renderRequest) || !isPlainObject(other.renderRequest)
    || !(/** @type {readonly unknown[]} */ (['circular', 'linear'])).includes(mode) || mode === sessionData.renderRequest.mode) {
    return [{ kind: 'request' }];
  }
  const setFaults = committedResultSetFaults(other, version, storedSchemas[at] ?? storedRequestSchema(other));
  if (setFaults.some(({ kind }) => kind === 'results' || kind === 'editorState')) return setFaults;
  /** @type {{ kind: 'set' | ResultSetFault['kind'] }[]} */
  const faults = [...setFaults];
  // The checks above proved `results` an array of named SVG Results and `editorState` an object.
  const results = /** @type {SessionResult[]} */ (other.results);
  const editorState = /** @type {Record<string, unknown>} */ (other.editorState);
  if (results.length === 0) faults.push({ kind: 'set' });
  const namesMissing = [...collectCanonicalResourceIds(other.renderRequest)]
    .some((resourceId) => !(isPlainObject(sessionData.resources) && Object.hasOwn(sessionData.resources, resourceId)));
  if (namesMissing) faults.push({ kind: 'set' });
  // One rule with gbdraw/session_io.py: a present field is an object (null is not).
  const ui = Object.hasOwn(other, 'ui') ? other.ui : {};
  const legend = Object.hasOwn(editorState, 'legend') ? editorState.legend : {};
  const runMetadata = Object.hasOwn(other, 'runMetadata') ? other.runMetadata : {};
  if (!isPlainObject(ui) || Object.keys(ui).some((key) => !OTHER_MODE_RESULT_UI_FIELDS.includes(key))
    || Object.keys(editorState).some((key) => !OTHER_MODE_RESULT_EDITOR_FIELDS.includes(key))
    || !isPlainObject(legend) || Object.keys(legend).some((key) => !OTHER_MODE_RESULT_LEGEND_FIELDS.includes(key))) {
    faults.push({ kind: 'set' });
  }
  if (!isPlainObject(runMetadata)
    || Object.keys(runMetadata).some((key) => !OTHER_MODE_RESULT_RUN_METADATA_FIELDS.includes(key))) {
    faults.push({ kind: 'set' });
  }
  return faults;
};

/**
 * Save: the stored Result sets of a Session 46-48 (E1).
 * @param {SessionInput} sessionData
 * @param {number | string} version
 * @param {StoredRequestSchemas} storedSchemas
 */
const validateStoredResultSets = (sessionData, version, storedSchemas) => {
  for (const stored of storedResultSets(sessionData).filter(({ projected }) => !projected)) {
    refuseFirstFault(storedSetFaults(stored, sessionData, version, storedSchemas));
    // The check above proved the set an object.
    validateRunNotices(stored.set);
  }
};

/**
 * A Result set's Results, read by index (unvalidated).
 * @param {SessionInput} set
 */
const storedResults = (set) => /** @type {{ readonly [index: number]: UnvalidatedFields<'name'> } | null | undefined} */ (set.results);

/**
 * The run notices of a Result set's run metadata, each with the check
 * Generate and Save apply to it against the set's Results. They are derived
 * run metadata (annotation warnings, feature identity notices, comparison
 * warnings): the next Generate reports them again, and the figure needs none.
 * @type {Readonly<Record<'annotationWarnings' | 'featureIdentityNotices' | 'comparisonWarnings', (value: unknown, set: SessionInput) => unknown>>}
 */
const RUN_NOTICE_CHECKS = Object.freeze({
  annotationWarnings: (value, set) => validateAnnotationWarnings(value, storedResults(set)),
  featureIdentityNotices: (value, set) => validateFeatureIdentityNotices(value, set.results),
  comparisonWarnings: (value, set) => validateComparisonWarnings(value, storedResults(set))
});
const RUN_NOTICE_FIELDS = /** @type {(keyof typeof RUN_NOTICE_CHECKS)[]} */ (Object.keys(RUN_NOTICE_CHECKS));

/**
 * Save: refuses run notices that do not match the set's Results.
 * @param {SessionInput} set A Session or the other mode's Result set.
 */
const validateRunNotices = (set) => {
  const runMetadata = /** @type {UnvalidatedFields<keyof typeof RUN_NOTICE_CHECKS>} */ (set.runMetadata);
  RUN_NOTICE_FIELDS.forEach((field) => RUN_NOTICE_CHECKS[field](runMetadata?.[field], set));
};

/**
 * Load: `runMetadata` without the run notices that do not match `results`,
 * dropped with no row (OV-638, J2, E-OD16); Generate and Save refuse them.
 * @template {SessionInput} T
 * @param {T} set A Session or the other mode's Result set.
 * @returns {T}
 */
const withReadRunNotices = (set) => {
  const runMetadata = set.runMetadata;
  if (!isPlainObject(runMetadata)) return set;
  const stale = RUN_NOTICE_FIELDS.filter((field) => {
    try {
      RUN_NOTICE_CHECKS[field](runMetadata[field], set);
      return false;
    } catch {
      // Each check is pure: it throws only for notices that do not hold.
      return true;
    }
  });
  if (stale.length === 0) return set;
  const kept = Object.fromEntries(Object.entries(runMetadata).filter(([field]) => !(/** @type {string[]} */ (stale)).includes(field)));
  return /** @type {T} */ ({ ...set, runMetadata: kept });
};

/**
 * Load: the committed Result sets of a Session, read after its requests
 * (the receipts are checked against the requests Load read, RK-2): what
 * Save refuses in them (`validateSessionAuthorityInventory`) is repaired
 * before the dialog (E-OD02). The other mode's set is dropped whole when it
 * cannot be restored: without its request (`other-result-set-dropped`, A7g),
 * or with a fault of its own (`session-part-dropped`).
 * A set whose Results Load could not admit or compose (A20, `uncomposedSets`:
 * its path prefix) has a `results` fault: Load reads the file again without
 * them (E-OD16), in every band (OV-789, E-OD17).
 * @template {SessionInput} T
 * @param {T} sessionData
 * @param {number} version
 * @param {StoredRequestSchemas} storedSchemas
 * @param {LoadRepair[]} repairs Collects the rows.
 * @param {readonly ResultSetPath[]} [uncomposedSets]
 * @returns {T}
 */
export const readSessionAuthorityParts = (sessionData, version, storedSchemas, repairs, uncomposedSets = []) => {
  /** @param {ResultSetPath} at @returns {ResultSetFault[]} */
  const uncomposed = (at) => (uncomposedSets.includes(at) ? [{ kind: 'results' }] : []);
  // A Session 27-39 holds no committed Result set Save checks: only Load's own fault.
  const topFaults = version >= CURRENT_AUTHORITY_SESSION_MIN_VERSION
    ? committedResultSetFaults(sessionData, version, storedSchemas[''] ?? storedRequestSchema(sessionData)) : [];
  /** @type {SessionInput} */
  let data = readCommittedResultSet(sessionData, [...topFaults, ...uncomposed('')], '', repairs);
  if (version < CANONICAL_SESSION_MIN_VERSION) return /** @type {T} */ (data);
  for (const stored of storedResultSets(data).filter(({ projected }) => !projected)) {
    const { at, set } = stored;
    const faults = [...storedSetFaults(stored, data, version, storedSchemas), ...uncomposed(at)];
    if (faults.some(({ kind }) => kind === 'request')) {
      const mode = /** @type {UnvalidatedFields<'renderRequest', UnvalidatedFields<'mode'>>} */ (set)?.renderRequest?.mode;
      repairs.push(loadRepair({
        code: 'other-result-set-dropped', path: `${at}renderRequest`, effect: 'dropped',
        ...(mode === 'circular' || mode === 'linear' ? { mode } : {})
      }));
      data = withStoredResultSet(data, at, null);
    } else if (faults.some(({ kind }) => !['receipt', 'receipt-missing'].includes(kind))) {
      repairs.push(loadRepair({ code: 'session-part-dropped', path: at.slice(0, -1), effect: 'dropped' }));
      data = withStoredResultSet(data, at, null);
    } else {
      // `storedSetFaults` found no `set` fault: the set is an object.
      // Without `set` and `request` faults, the rest are the set's own (`committedResultSetFaults`).
      const setFaults = /** @type {ResultSetFault[]} */ (faults.filter(({ kind }) => kind !== 'set' && kind !== 'request'));
      data = withStoredResultSet(data, at, withReadRunNotices(readCommittedResultSet(set, setFaults, at, repairs)));
    }
  }
  return /** @type {T} */ (withReadRunNotices(data));
};

/**
 * Load: one drawing's committed Result set of a Session 49, read after its
 * request (`readSessionRequests`) as `readSessionAuthorityParts` reads a
 * Session 46-48's: the drawing Load projects keeps what can be restored of
 * its set, as the top-level set does; any other drawing's set is dropped on
 * a fault of its own (`session-part-dropped`), as the other mode's set is.
 * Null drops the set (the drawing stays).
 * @param {{ at: ResultSetPath, projected: boolean, set: SessionInput }} stored
 * @param {number} version
 * @param {unknown} storedSchema The schema the file stored for the set's request (RK-1).
 * @param {LoadRepair[]} repairs Collects the rows.
 * @param {readonly ResultSetPath[]} [uncomposedSets]
 * @returns {SessionInput | null}
 */
export const readDrawingResultSetParts = ({ at, projected, set }, version, storedSchema, repairs, uncomposedSets = []) => {
  /** @type {ResultSetFault[]} */
  const faults = [
    ...committedResultSetFaults(set, version, storedSchema ?? storedRequestSchema(set)),
    ...(uncomposedSets.includes(at) ? [/** @type {ResultSetFault} */ ({ kind: 'results' })] : [])
  ];
  if (!projected && faults.some(({ kind }) => !['receipt', 'receipt-missing'].includes(kind))) {
    repairs.push(loadRepair({ code: 'session-part-dropped', path: at.slice(0, -1), effect: 'dropped' }));
    return null;
  }
  return withReadRunNotices(readCommittedResultSet(set, faults, at, repairs));
};

// Session 46 keeps each mode's draft in `modes` (registry fields only) and
// none of it in its former homes; `ui.losatExecution` holds only the LOSAT
// execution settings (plan 4.1).
/** Save: @param {Record<string, unknown>} sessionData */
const validateModeScopedFields = (sessionData) => {
  // A field outside the Session 46 shape is a Session-format error (plan 4.1).
  const fieldsInvalid = () => diagnosticError('INPUT_INVALID', { field: 'schema', reason: 'FIELDS' });
  if (Object.hasOwn(sessionData, 'cliOptions') && !isPlainObject(sessionData.cliOptions)) throw fieldsInvalid();
  const retired = RETIRED_MODE_SCOPED_FIELDS.flatMap(([domain, fields]) => {
    const container = domain ? domain.split('.').reduce(
      (/** @type {unknown} */ current, part) => (isPlainObject(current) ? current[part] : undefined), sessionData
    ) : sessionData;
    return isPlainObject(container)
      ? Object.keys(container).filter((field) => fields.has(field)).map((field) => (domain ? `${domain}.${field}` : field))
      : [];
  });
  // Session 46 keeps the draft of each mode in `modes`, never in a former home.
  if (retired.length > 0) throw fieldsInvalid();
  // `ui.losatExecution` holds only the LOSAT execution settings.
  const execution = /** @type {UnvalidatedFields<'losatExecution'>} */ (sessionData.ui)?.losatExecution;
  if (execution !== undefined && (!isPlainObject(execution)
    || Object.keys(execution).some((field) => !LOSAT_EXECUTION_FIELDS.includes(field)))) throw fieldsInvalid();
  if (!Object.hasOwn(sessionData, 'modes')) return;
  // `modes` holds a circular and a linear slice only.
  if (!isPlainObject(sessionData.modes)
    || Object.keys(sessionData.modes).some((mode) => !(/** @type {readonly string[]} */ (SLICE_MODES)).includes(mode))) throw fieldsInvalid();
  Object.values(sessionData.modes).forEach((slice) => validateModeSliceFields(slice));
};

/**
 * Save: the keys and kinds of the band of the version it writes. Load reads
 * them (`readSessionDocument`): it drops what Save refuses here (E-OD02).
 * @param {Record<string, unknown>} sessionData
 * @param {number | string} version
 */
const validateBandFields = (sessionData, version) => {
  if (Number(version) < MODE_SCOPED_SESSION_VERSION) {
    // The stored Result sets (`storedResultSets`) are named by their path.
    const fields = [...storedResultSets(sessionData).filter(({ projected }) => !projected).map(({ at }) => at.slice(0, -1)),
      ...MODE_SCOPED_TOP_LEVEL_KEYS];
    for (const field of fields) {
      if (Object.hasOwn(sessionData, field)) {
        throw new Error(`Session version ${String(version)} cannot contain ${field}.`);
      }
    }
  } else {
    validateModeScopedFields(sessionData);
  }
  if (Number(version) < CANONICAL_SESSION_MIN_VERSION) return;
  const unknown = Object.keys(sessionData).filter((key) => !Object.hasOwn(SESSION_TOP_LEVEL_AUTHORITY, key));
  if (unknown.length > 0) {
    throw new Error(`Session contains unclassified top-level field(s): ${unknown.join(', ')}`);
  }
  if (Number(version) < CURRENT_AUTHORITY_SESSION_MIN_VERSION) return;
  if (Object.hasOwn(sessionData, 'files')) {
    throw new Error(`Session version ${String(version)} cannot contain legacy files; use resources and webFiles.`);
  }
  if (Object.hasOwn(sessionData, 'features') && !isPlainObject(sessionData.features)) {
    throw new Error('Session features must be an object when present.');
  }
  const features = isPlainObject(sessionData.features) ? sessionData.features : {};
  const forbiddenFeatureFields = CURRENT_WRITER_FORBIDDEN_FEATURE_FIELDS.filter((field) => Object.hasOwn(features, field));
  if (forbiddenFeatureFields.length > 0) {
    throw new Error(`Session version ${String(version)} contains branch-only feature field(s): ${forbiddenFeatureFields.join(', ')}`);
  }
  if (Number(version) >= MODE_SCOPED_SESSION_VERSION) {
    SLICE_MODES.forEach((mode) => {
      const sliceFeatures = /** @type {UnvalidatedFields<string, UnvalidatedFields<'features'>>} */ (sessionData.modes)?.[mode]?.features;
      if (!isPlainObject(sliceFeatures) || !Object.hasOwn(sliceFeatures, 'featureOverrides')) return;
      if (!isPlainObject(sliceFeatures.featureOverrides)) {
        throw new Error(`Session modes.${mode}.features.featureOverrides must be an object.`);
      }
      canonicalFeatureOverrides(sliceFeatures.featureOverrides);
    });
  }
  if (Object.hasOwn(sessionData, 'orthogroupState') && !isPlainObject(sessionData.orthogroupState)) {
    throw new Error('Session orthogroupState must be an object when present.');
  }
  if (isPlainObject(sessionData.orthogroupState) && Object.hasOwn(sessionData.orthogroupState, 'groups')) {
    throw new Error(`Session version ${String(version)} cannot contain duplicated orthogroup groups.`);
  }
};

/**
 * Save's check of the Session it writes, and Load's (`load`), which leaves the
 * values of the drawings' settings to the reader table: Load drops what it
 * cannot read with a dialog row before it applies the file (E-OD02).
 * @param {SessionEnvelope} sessionData Unvalidated Session data.
 * @param {number | string} version
 * @param {StoredRequestSchemas} [storedSchemas] The request schemas the file stored, when Load read them.
 * @param {{ load?: boolean }} [options]
 * @returns {void}
 */
export const validateSessionAuthorityInventory = (sessionData, version, storedSchemas = {}, { load = false } = {}) => {
  if (!sessionData || typeof sessionData !== 'object' || Array.isArray(sessionData)) {
    throw new Error('Session authority inventory requires an object.');
  }
  assertSafeObjectKeys(sessionData);
  // Load read the keys and kinds of its band (`readSessionDocument`).
  if (!load) validateBandFields(sessionData, version);
  const bindings = validateWebFileBindings(
    /** @type {UnvalidatedFields<'bindings'>} */ (sessionData.webFiles),
    /** @type {UnvalidatedFields<string, UnvalidatedFields<keyof RawSessionResourceDescriptor>>} */ (sessionData.resources),
    webBindingInputSlot(Number(version) >= SESSION_DRAWINGS_LAYOUT.version, /** @type {UnvalidatedFields<'mode'>} */ (sessionData.ui)?.mode)
  );
  // A Linear record the Session turns OFF names a bound record card (record-selection D-05).
  const linearRecordsOff = /** @type {UnvalidatedFields<'linear', UnvalidatedFields<'config', UnvalidatedFields<'recordsOff'>>>} */ (
    sessionData.modes
  )?.linear?.config?.recordsOff;
  // Load drops a record that names no card (`readRecordsOff`).
  if (!load && Number(version) >= MODE_SCOPED_SESSION_VERSION && Array.isArray(linearRecordsOff)) {
    const boundUids = new Set(/** @type {unknown[]} */ (Array.isArray(bindings?.linearSeqs) ? bindings.linearSeqs : []).map((sequence) => /** @type {UnvalidatedFields<'uid'>} */ (sequence)?.uid));
    if (/** @type {unknown[]} */ (linearRecordsOff).some((uid) => !boundUids.has(uid))) {
      throw diagnosticError('INPUT_INVALID', { field: 'schema', reason: 'FIELDS' });
    }
  }
  // A Session 49 drawing's one-drawing view holds its bindings (gbdraw/session_io.py).
  if (bindings?.schema === 2 && ![41, 42, 44, MODE_SCOPED_SESSION_VERSION, ROW_KEYED_LEGEND_SESSION_VERSION, SESSION_VERSION,
    SESSION_DRAWINGS_LAYOUT.version].includes(Number(version))) {
    throw new Error('Web binding schema 2 requires session version 41, 42, 44, '
      + `${MODE_SCOPED_SESSION_VERSION}, ${ROW_KEYED_LEGEND_SESSION_VERSION}, ${SESSION_VERSION}, or ${SESSION_DRAWINGS_LAYOUT.version}.`);
  }
  if (Number(version) < CANONICAL_SESSION_MIN_VERSION) return;
  if (Number(version) >= CURRENT_AUTHORITY_SESSION_MIN_VERSION) {
    validateCurrentComparisonAuthority(sessionData, { load });
    // Load reads an input type of no type as GenBank, with a row (`readSessionInputType`, OV-594).
    for (const field of load ? [] : ['cInputType', 'lInputType']) {
      if (
        Object.prototype.hasOwnProperty.call(sessionData.ui || {}, field)
        // The own-property check passed, so `ui` is an object that holds `field`.
        && !(/** @type {readonly unknown[]} */ (SESSION_INPUT_TYPES)).includes(/** @type {Record<string, unknown>} */ (sessionData.ui)[field])
      ) {
        throw new Error(`Session ui.${field} must be gb or gff when present.`);
      }
    }
    // Load read the committed Result sets (`readSessionAuthorityParts`): it
    // keeps Results whose catalog it dropped, which Save refuses.
    if (!load) {
      validateCommittedResultSet(sessionData, version,
        storedSchemas[''] ?? /** @type {UnvalidatedFields<'schema'>} */ (sessionData.renderRequest)?.schema);
    }
  }
  if (!load) validateStoredResultSets(sessionData, version, storedSchemas);
  // Load read the run notices (`readSessionAuthorityParts`). A Session 31-39
  // reaches these checks with `results` unvalidated: they read it by index only.
  if (!load) validateRunNotices(sessionData);
  if (isSettingsOnlySessionDocument(sessionData)) validateSettingsOnlyDocument(sessionData, { load });
};

/**
 * @param {Record<string, any>} canonical `{ renderRequest, resources, webFiles }`.
 * @returns {Record<string, any>} The same object, adopted.
 */
export const adoptRuntimeCanonicalSession = (canonical) => {
  if (
    !isPlainObject(canonical)
    || !isPlainObject(canonical.renderRequest)
    || !isPlainObject(canonical.resources)
    || (
      canonical.webFiles !== undefined
      && !isPlainObject(canonical.webFiles)
    )
  ) {
    throw new Error('A canonical render request is required for adoptive ownership.');
  }
  adoptedCanonicalOwners.add(canonical);
  return canonical;
};

/**
 * @param {SessionEnvelope} sessionData Unvalidated Session data.
 * @param {number} currentVersion
 * @param {StoredRequestSchemas} [storedSchemas]
 * @returns {{ document: SessionEnvelope, canonical: Record<string, any> | null }}
 */
export const adoptCurrentSessionDocument = (sessionData, currentVersion, storedSchemas = {}) => {
  validateSessionAuthorityInventory(sessionData, currentVersion, storedSchemas, { load: true });
  if (sessionData.version !== currentVersion) {
    throw new Error('Only the current session schema can use adoptive ownership.');
  }
  // A Session 49 drawing without a request projects as a settings-only Session.
  const requestLess = isSettingsOnlySessionDocument(sessionData)
    || (currentVersion >= SESSION_DRAWINGS_LAYOUT.version && sessionData.renderRequest == null);
  const canonical = requestLess ? null : adoptRuntimeCanonicalSession({
    renderRequest: sessionData.renderRequest,
    resources: sessionData.resources,
    webFiles: isPlainObject(sessionData.webFiles) ? sessionData.webFiles : {}
  });
  adoptedCurrentDocuments.add(sessionData);
  return { document: sessionData, canonical };
};

/**
 * @param {any} value
 * @returns {boolean}
 */
export const isAdoptedCurrentSessionDocument = (value) => (
  Boolean(value) && adoptedCurrentDocuments.has(value)
);

/**
 * @param {any} value
 * @returns {boolean}
 */
export const isAdoptedCanonicalSession = (value) => (
  Boolean(value) && adoptedCanonicalOwners.has(value)
);

/**
 * @param {Record<string, any> | null | undefined} sessionData
 * @returns {{ ui: Record<string, any> }}
 */
export const projectWebOnlyEditorMetadata = (sessionData) => ({
  ui: copyFields(sessionData?.ui, WEB_EDITOR_UI_FIELDS)
});

/**
 * @param {Record<string, any> | null | undefined} sessionData
 * @returns {Record<string, any>}
 */
export const projectArtifactState = (sessionData) => ({
  results: Array.isArray(sessionData?.results) ? sessionData.results : [],
  ui: copyFields(sessionData?.ui, ARTIFACT_UI_FIELDS),
  features: copyFields(sessionData?.features, ARTIFACT_FEATURE_FIELDS),
  editorState: sessionData?.editorState || {},
  orthogroupState: sessionData?.orthogroupState || {},
  losatCache: sessionData?.losatCache || {},
  losatDerivedCache: sessionData?.losatDerivedCache || {},
  proteinIdentityManifest: sessionData?.proteinIdentityManifest || {},
  legacyArtifacts: sessionData?.legacyArtifacts || {},
  runMetadata: sessionData?.runMetadata || {}
});

/**
 * @param {Record<string, any> | null | undefined} sessionData
 * @returns {{ format: any, version: any, createdAt: any, title: string }}
 */
export const projectDocumentMetadata = (sessionData) => ({
  format: sessionData?.format,
  version: sessionData?.version,
  createdAt: sessionData?.createdAt,
  title: typeof sessionData?.title === 'string' ? sessionData.title : ''
});

/** @import { Session49 } from './session-drawings-upgrade.js' */
/**
 * The drawing a Session 49 opens (`_opened_drawing`, gbdraw/session_drawings.py):
 * the one `ui.activeDrawingId` names, else the first; a `ui` that is no
 * object names none (the band reader admits any `ui`). The one place the Web
 * reads which drawing a Session 49 opens.
 * @template {{ id: string }} D
 * @param {{ ui?: unknown, drawings: readonly D[] }} document
 * @returns {D}
 */
export const openedDrawingOf = (document) => {
  const [first] = document.drawings;
  const activeId = isPlainObject(document.ui) ? document.ui.activeDrawingId : undefined;
  return document.drawings.find((entry) => entry.id === activeId) ?? first;
};

/**
 * One committed Result set of a Session as the file holds it: `at` is where
 * (the prefix of its paths, so its rows keep the file's own paths), and
 * `projected` says it is the set Load projects (the top-level set of a
 * Session 27-48) rather than a stored set, which loads only as saved.
 * @typedef {{ at: ResultSetPath, projected: boolean, set: SessionInput }} StoredResultSet
 */
/**
 * The project fields a Result set's readers read beside the set.
 * @typedef {{ resources: unknown, orthogroupState?: unknown }} StoredSetProject
 */
/**
 * The committed Result sets of a Session, the projected one first: a Session
 * 27-48's top-level set and a Session 46-48's `otherModeResult`; each Session
 * 49 drawing that holds a request, the one it opens (`openedDrawingOf`)
 * projected (a drawing is found by its ID). The one place Load names where a
 * file holds a set; each set reader below reads one.
 * @param {SessionInput} data
 * @returns {StoredResultSet[]}
 */
export const storedResultSets = (data) => {
  if (isDrawingsBand(/** @type {SessionChecked} */ (data))) {
    const document = /** @type {Session49} */ (/** @type {unknown} */ (data));
    // A Session 49 whose every drawing was dropped holds no set (R4, `restoresNothing`).
    const openedId = document.drawings.length > 0 ? openedDrawingOf(document).id : null;
    /** @type {StoredResultSet[]} */
    const sets = document.drawings.flatMap((entry, index) => (entry.renderRequest == null ? [] : [{
      at: drawingSetPath(index),
      projected: entry.id === openedId,
      set: /** @type {SessionInput} */ (/** @type {unknown} */ (entry))
    }]));
    return [...sets.filter(({ projected }) => projected), ...sets.filter(({ projected }) => !projected)];
  }
  return [
    { at: '', projected: true, set: data },
    ...(isPlainObject(data.otherModeResult)
      ? [{ at: /** @type {ResultSetPath} */ ('otherModeResult.'), projected: false, set: /** @type {SessionInput} */ (data.otherModeResult) }]
      : [])
  ];
};

/**
 * `data` with its set at `at` replaced by `set`; a null `set` leaves a stored
 * set out (the projected set of a Session 27-48 is never left out). A Session
 * 49 drawing whose set is left out stays, without its set: no request, no
 * Results, and none of the set's catalog, run metadata or invocation.
 * @template {SessionInput} T
 * @param {T} data @param {ResultSetPath} at @param {SessionInput | null} set
 * @returns {T}
 */
export const withStoredResultSet = (data, at, set) => {
  const drawingIndex = drawingSetIndex(at);
  if (drawingIndex !== null) {
    const entries = [.../** @type {SessionInput[]} */ (/** @type {Record<string, unknown>} */ (data).drawings)];
    if (set !== null) {
      entries[drawingIndex] = set;
    } else {
      const { runMetadata: _runMetadata, cliInvocation: _cliInvocation, ...drawing } = entries[drawingIndex];
      const { featureCatalog: _catalog, ...editorState } = isPlainObject(drawing.editorState) ? drawing.editorState : {};
      entries[drawingIndex] = {
        ...drawing, renderRequest: null, results: [],
        ...(Object.hasOwn(drawing, 'editorState') ? { editorState } : {})
      };
    }
    return { ...data, drawings: entries };
  }
  if (at === '') {
    if (set === null) throw new Error('The projected Result set is never left out.');
    return /** @type {T} */ (set);
  }
  if (set !== null) return { ...data, otherModeResult: set };
  const { otherModeResult: _dropped, ...kept } = data;
  return /** @type {T} */ (kept);
};

/**
 * The request schemas the file stored (RK-1), under the key of where it
 * holds each set.
 * @param {readonly StoredResultSet[]} sets
 * @returns {StoredRequestSchemas}
 */
export const storedRequestSchemas = (sets) => Object.fromEntries(sets.map(({ at, projected, set }) => [
  projected ? '' : at,
  /** @type {UnvalidatedFields<'schema'>} */ (set.renderRequest)?.schema
]));
