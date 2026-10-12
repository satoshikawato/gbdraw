// @ts-check
/** @import { FeaturePlacementDraft } from './feature-placement.js' */
/** @import { RawSessionResourceDescriptor } from './session-resource-backing.js' */
import { adoptCurrentSessionResources, createSessionResourceFileView, readSessionResourceBytes } from './session-resource-backing.js';
import { sha256Hex, textToBytes } from './byte-utils.js';
import { canonicalRecordReverseComplement } from './record-display-model.js';
import { createDefaultLinearDefinitionLineStyles } from './definition-line-style-state.js'; import { CIRCULAR_TRACK_RENDERERS, createDefaultCircularTrackSlots } from './circular-track-slot-model.js';
import { LEGACY_LINEAR_TRACK_SLOT_SCHEMA_VERSION, LINEAR_TRACK_RENDERERS, LINEAR_TRACK_SLOT_SCHEMA_VERSION, createDefaultLinearTrackSlots } from './linear-track-slot-model.js'; import { validateTrackSlotBindingInvariants } from './track-slot-validation.js';
import { requireCurrentCircularMultiRecordSizeMode, requireCurrentCollinearAnchorMode, requireCurrentCollinearColorMode, requireCurrentCollinearMaxConflicts, requireCurrentCollinearMaxDiagonalDrift, requireCurrentCollinearMaxParalogLinks, requireCurrentCollinearMaxUnitGap, requireCurrentCollinearMergeOrientation, requireCurrentCollinearMinAnchors, requireCurrentCollinearInferOrthogroups, requireCurrentCollinearSearchScope, requireCurrentCollinearUnitMode, requireCurrentLinearLabelPlacement, requireCurrentLinearTrackLayout, requireCurrentOrthogroupMemberMaxHits, requireCurrentOrthogroupMembershipMode, requireCurrentProteinBlastpCandidateLimit, requireCurrentProteinBlastpMaxHits, requireCurrentProteinBlastpMode, requireCurrentWebStateFieldNames } from './current-option-values.js'; import { DEFAULT_ARROW_SHAFT_WIDTH_RATIO, createDefaultFeatureRenderings } from '../utils/feature-rendering.js';
import { MODE_DEFAULT_FEATURE_TYPES, comparisonStateForMode, managedAdvStateForMode, trackDefaultsForMode } from '../mode-profiles.js'; import { WEB_UX_PROFILE } from '../web-ux-profile.js';
import { assertSafeObjectKeys } from './safe-object-keys.js';
import { diagnosticError } from '../utils/error-normalization.js';
import { validateRecordDisplayDrafts } from './record-display-model.js';
import { isRecordDrawKey } from './record-draw-selection.js';
import { requireLinearLabelVisibilityMode } from './linear-label-visibility.js';
import { canonicalFeaturePlacements } from './feature-placement.js';
import { validateScopedFeaturePlacements } from './feature-edit-migration.js';
import { ORTHOGROUP_INFERENCE_VERSION } from '../contracts/schema.generated.js';
const circularTracks = trackDefaultsForMode('circular'), linearTracks = trackDefaultsForMode('linear');
export const CIRCULAR_TRACK_SLOT_SCHEMA_VERSION = 4, LEGACY_CIRCULAR_TRACK_SLOT_SCHEMA_VERSION = 3;
/**
 * @template T
 * @param {T} value
 * @returns {value is T & Record<string, unknown>}
 */
const isObject = (value) => Boolean(value) && typeof value === 'object' && !Array.isArray(value), has = /** @type {(value: unknown, key: string) => boolean} */ ((value, key) => Object.prototype.hasOwnProperty.call(value, key));
/**
 * The Session's active Web configuration (`config` of the stored Session) in
 * its current writer format: the domains of `CURRENT_WRITER_ACTIVE_CONFIG_DOMAINS`
 * plus the compatibility field `colorsAreOverrides`. `form` and `adv` are
 * required; the option fields inside them belong to Python (R7) and stay
 * `Record<string, any>`. A reader takes unvalidated data and returns this type.
 * @typedef {object} ActiveWebConfig
 * @property {Record<string, any>} form
 * @property {Record<string, any>} adv
 * @property {Record<string, any>} [losat]
 * @property {Record<string, any>} [cliOptions]
 * @property {Record<string, any>} [colors]
 * @property {string} [palette]
 * @property {boolean} [paletteInstantPreviewEnabled]
 * @property {Array<Record<string, any>>} [rules]
 * @property {Array<Record<string, any>>} [qualifierPriorityRules]
 * @property {string} [filterMode] `None`, `Whitelist`, or `Blacklist`.
 * @property {Array<Record<string, any>>} [whitelist]
 * @property {string} [blacklistText]
 * @property {string} [losatProgram]
 * @property {Record<string, any>} [circularConservation]
 * @property {Array<Record<string, any>>} [annotationSets]
 * @property {Array<Record<string, any>>} [recordDisplayDrafts]
 * @property {FeaturePlacementDraft} [featurePlacementOverrides]
 * @property {Record<string, any>} [modeProfiles]
 * @property {Record<string, any>} [unmanagedConfigOverrides]
 * @property {Record<string, any>} [linearRecordLayout]
 * @property {Record<string, any>} [linearComparisonPlan]
 * @property {Record<string, any>} [importedComparisonResolution]
 * @property {Record<string, any>} [webEdits]
 * @property {boolean} [colorsAreOverrides]
 */

export const createDefaultForm = () => ({
  prefix: '', species: '', strain: '', plot_title: '', track_type: 'tuckin', linear_track_layout: 'middle', show_scale: true, scale_style: 'bar',
  linear_ruler_on_axis: false, labels_mode: 'none', show_labels_linear: 'none', multi_record_canvas: WEB_UX_PROFILE.circular.gridByDefault,
  circular_record_selector: '', circular_region_start: null, circular_region_end: null, circular_reverse: false,
  circular_record_label: '', circular_record_subtitle: '',
  separate_strands: WEB_UX_PROFILE.separateStrands, suppress_gc: !circularTracks.gc, suppress_skew: !circularTracks.skew, align_center: false,
  keep_definition_left_aligned: true, show_gc: linearTracks.gc, show_skew: linearTracks.skew, show_depth: false, normalize_length: false
});
/**
 * @param {string} [mode] `circular` or `linear`.
 * @returns {Record<string, any>}
 */
export const createDefaultAdv = (mode = 'circular') => ({
  features: [...MODE_DEFAULT_FEATURE_TYPES], feature_shapes: createDefaultFeatureRenderings(), arrow_head_length_ratio: null,
  arrow_shaft_width_ratio: DEFAULT_ARROW_SHAFT_WIDTH_RATIO, window_size: null, step_size: null, nt: 'GC', def_font_size: null,
  circular_definition_interval: null, label_font_size: null, circular_label_spacing: null, linear_label_spacing: null, label_rendering: 'auto',
  circular_label_placement: 'horizontal', label_placement: 'auto', label_rotation: null, block_stroke_width: null, block_stroke_color: null,
  line_stroke_width: null, line_stroke_color: null, axis_stroke_width: null, axis_stroke_color: managedAdvStateForMode(mode).axis_stroke_color,
  legend_box_size: null, legend_font_size: null, resolve_overlaps: false, feature_overlap_tolerance_bp: 0, feature_height: null, track_axis_gap: null, linear_show_replicon: false,
  linear_accession_visibility: 'auto', linear_length_visibility: 'auto', linear_definition_line_styles: createDefaultLinearDefinitionLineStyles(), gc_height: null,
  depth_height: null, depth_color: '#4A90E2', depth_tracks: [], depth_window_size: null, depth_step_size: null, depth_share_axis: false,
  depth_min: null, depth_max: null, depth_normalize: false, depth_show_axis: true, depth_show_ticks: true, depth_large_tick_interval: null,
  depth_small_tick_interval: null, depth_tick_font_size: null, linear_track_slots_enabled: false, linear_track_slots_schema_version: LINEAR_TRACK_SLOT_SCHEMA_VERSION,
  linear_track_slots_axis_index: null, linear_track_slots: createDefaultLinearTrackSlots(), gc_content_mode: 'deviation', gc_content_min_percent: 0,
  gc_content_max_percent: 100, gc_content_show_axis: true, gc_content_show_ticks: true, gc_content_tick_interval: 20, gc_content_small_tick_interval: null,
  gc_content_tick_font_size: null, comparison_height: null, pairwise_match_style: mode === 'linear' ? 'curve' : 'ribbon', ...comparisonStateForMode(mode), scale_interval: null,
  scale_font_size: null, ruler_label_font_size: null, scale_stroke_width: null, scale_stroke_color: null, ruler_label_color: null, circular_grouping_intent: 'auto',
  multi_record_size_mode: 'auto', multi_record_min_radius_ratio: 0.55, multi_record_column_gap_ratio: 0.10, multi_record_row_gap_ratio: 0.05,
  multi_record_positions: [], tick_label_font_size: null, plot_title_font_size: null, keep_full_definition_with_plot_title: false,
  center_reserved_radius: null, feature_width_circular: null, depth_width_circular: null, gc_content_width_circular: null, gc_content_radius_circular: null,
  gc_skew_width_circular: null, gc_skew_radius_circular: null, circular_track_slots_enabled: false,
  circular_track_slots_schema_version: CIRCULAR_TRACK_SLOT_SCHEMA_VERSION, circular_track_slots_axis_index: null, circular_track_slots: createDefaultCircularTrackSlots(),
  outer_label_x_offset: null, outer_label_y_offset: null, inner_label_x_offset: null, inner_label_y_offset: null
});
export const createDefaultLosatpHitLimits = () => ({
  pairwise: { candidateLimit: null, orthogroupMemberMaxHits: null },
  orthogroup: { candidateLimit: null, orthogroupMemberMaxHits: null },
  collinear: { candidateLimit: 5, orthogroupMemberMaxHits: 5 }
});
// How LOSAT runs: one app-level setting for both drawings (Session `ui.losatExecution`).
export const createDefaultLosatExecution = () => ({
  executionMode: 'threaded', totalThreadBudget: 'safe', threadsPerJob: 'auto', parallelWorkers: undefined
});
export const createDefaultLosat = () => ({
  outfmt: '6',
  blastn: { task: 'megablast' }, blastp: { mode: 'orthogroup', hitLimitsByMode: createDefaultLosatpHitLimits(), maxHits: 5, candidateLimit: null, orthogroupMembershipMode: ORTHOGROUP_INFERENCE_VERSION,
    orthogroupMemberMaxHits: null, collinearInferOrthogroups: false, collinearMinAnchors: 1, collinearMaxUnitGap: 0, collinearMaxDiagonalDrift: 0,
    collinearMaxConflictsInMergeGap: 1, collinearMaxParalogLinksPerOrthogroup: 2, collinearColorMode: 'orientation',
    collinearUnitMode: 'auto', collinearAnchorMode: 'rbh', collinearMergeOrientation: 'either', collinearSearchScope: 'adjacent' }
});
export const createDefaultCircularConservation = () => ({ enabled: false, source: 'losat', losat_program: 'blastn',
  subject_gencode: 1, reference: 'auto', labels: '', series: [], ring_width: null, ring_gap: null });
export const CURRENT_WRITER_ACTIVE_CONFIG_DOMAINS = Object.freeze([
  'form', 'adv', 'losat', 'cliOptions', 'colors', 'palette', 'paletteInstantPreviewEnabled', 'rules',
  'qualifierPriorityRules', 'filterMode', 'whitelist', 'blacklistText', 'losatProgram', 'circularConservation',
  'annotationSets', 'recordDisplayDrafts', 'recordsOff', 'featurePlacementOverrides', 'modeProfiles', 'unmanagedConfigOverrides',
  'linearRecordLayout', 'linearComparisonPlan', 'importedComparisonResolution', 'webEdits'
]);
// The `config` keys the CLI writer of Sessions 40 and 41 derived from its
// options; a Web draft of those Sessions always has `annotationSets`,
// `linearComparisonPlan` and `webEdits` too (OV-269, `_CLI_WRITER_CONFIG_DOMAINS`
// in gbdraw/session_io.py).
const CLI_WRITER_CONFIG_DOMAINS = new Set(['form', 'adv', 'losat', 'cliOptions', 'colors', 'palette', 'rules',
  'qualifierPriorityRules', 'filterMode', 'whitelist', 'blacklistText', 'losatProgram', 'circularConservation']);
/**
 * A Session 40 or 41 the CLI wrote holds a `config` derived from its options,
 * no Web draft: Load reads its tables from the request.
 * @param {Record<string, any>} session
 * @returns {boolean}
 */
export const holdsCliWriterConfig = (session) => [40, 41].includes(session?.version)
  && session.cliInvocation?.generatedBy === 'gbdraw' && isObject(session.config)
  && Object.keys(session.config).every((key) => CLI_WRITER_CONFIG_DOMAINS.has(key));
export const CURRENT_WRITER_FORM_FIELDS = Object.freeze([...Object.keys(createDefaultForm()), 'legend']);
// `rich_feature_popup` is read from a Session 40-44 draft only: Session 46
// keeps the app preference in `ui.richFeaturePopup`.
export const CURRENT_WRITER_ADV_FIELDS = Object.freeze([...Object.keys(createDefaultAdv()), 'plot_title_position', 'losatProgram', 'rich_feature_popup']);
const DOMAIN_SHAPES = Object.freeze({ form: 'object', adv: 'object', losat: 'object', cliOptions: 'object', colors: 'object',
  circularConservation: 'object', modeProfiles: 'object', linearRecordLayout: 'object', linearComparisonPlan: 'object', webEdits: 'object',
  importedComparisonResolution: 'object', unmanagedConfigOverrides: 'object',
  recordDisplayDrafts: 'array', recordsOff: 'array', featurePlacementOverrides: 'object', rules: 'array', qualifierPriorityRules: 'array', whitelist: 'array', annotationSets: 'array', palette: 'string', filterMode: 'string', blacklistText: 'string',
  losatProgram: 'string', paletteInstantPreviewEnabled: 'boolean' });
const ROW_FIELDS = { rules: ['feat', 'qual', 'val', 'color', 'cap', 'fromFile'],
  qualifierPriorityRules: ['feat', 'order'], whitelist: ['feat', 'qual', 'key'],
  annotationSets: ['id', 'annotations', 'defaultStyle', 'legendLabel'] };
const ANNOTATION_FIELDS = ['id', 'target', 'label', 'mark', 'lane', 'style', 'legendLabel', 'metadata'];
const OBSOLETE_SLOT_KEYS = new Set(['gapAfter', 'gap_after', 'innerRadius', 'inner_radius', 'outerRadius',
  'outer_radius', 'placement', 'spacing', 'strict', 'compress', 'reserve']);
const OBSOLETE_SLOT_PARAM_KEYS = new Set(['side', 'radius', 'width', 'spacing', 'inner_gap_px', 'outer_gap_px', 'strict', 'compress', 'reserve']);
/** @param {Record<string, unknown>} value @param {Set<string>} fields @param {string} path */
const assertFields = (value, fields, path) => {
  const unknown = Object.keys(value).filter((key) => !fields.has(key));
  if (unknown.length) throw new Error(`Current session active configuration contains unknown ${path} field(s): ${unknown.join(', ')}.`);
};
/** @param {Record<string, unknown>} config */
const validateDomainShapes = (config) => {
  for (const [domain, shape] of Object.entries(DOMAIN_SHAPES)) {
    if (!has(config, domain) || (config[domain] === undefined && domain === 'cliOptions')) continue;
    const value = config[domain], valid = shape === 'object' ? isObject(value) : shape === 'array' ? Array.isArray(value) : typeof value === shape;
    if (!valid) throw new Error(`Current session active configuration config.${domain} must be ${shape}.`);
  }
  if (has(config, 'colorsAreOverrides') && typeof config.colorsAreOverrides !== 'boolean')
    throw new Error('Current session compatibility field config.colorsAreOverrides must be a boolean.');
};
/** @param {unknown[]} rows @param {string[]} fields @param {string} path */
const assertRows = (rows, fields, path) => rows.forEach((row, index) => {
  if (!isObject(row)) throw new Error(`Current session active configuration ${path}[${index}] must be an object.`);
  assertFields(row, new Set(fields), `${path}[${index}]`);
});
/**
 * The first loop checks that every annotation set is an object, which the
 * second loop reads.
 * @param {{ annotationSets?: { annotations?: unknown }[], [domain: string]: unknown }} config
 */
const validateCollections = (config) => {
  for (const [domain, fields] of Object.entries(ROW_FIELDS))
    if (Array.isArray(config[domain])) assertRows(config[domain], fields, `config.${domain}`);
  (config.annotationSets || []).forEach((set, index) => {
    if (!Array.isArray(set.annotations)) throw new Error(`Current session active configuration config.annotationSets[${index}].annotations must be an array.`);
    assertRows(set.annotations, ANNOTATION_FIELDS, `config.annotationSets[${index}].annotations`);
  });
};
/** @param {unknown} slots */
const obsoleteCircularSlotField = (slots) => {
  if (!Array.isArray(slots)) return null;
  for (const [slotIndex, slot] of slots.entries()) {
    if (!isObject(slot)) continue;
    const field = [...OBSOLETE_SLOT_KEYS].find((key) => has(slot, key))
      ?? (isObject(slot.params) ? [...OBSOLETE_SLOT_PARAM_KEYS].find((key) => has(slot.params, key)) : undefined);
    if (field) return { field, slotIndex };
  }
  return null;
};
/**
 * @param {Record<string, any>} [config] The active configuration.
 * @param {{ depthTrackCount?: number | null }} [options]
 * @returns {void}
 */
export const validateImportedCircularTrackSlots = (config = {}, { depthTrackCount = null } = {}) => {
  const adv = isObject(config) ? config.adv : null;
  if (!isObject(adv) || !has(adv, 'circular_track_slots')) return;
  if (adv.circular_track_slots_schema_version !== CIRCULAR_TRACK_SLOT_SCHEMA_VERSION)
    throw new Error(`Custom Track Slots use an obsolete schema. Recreate the slots with schema version ${CIRCULAR_TRACK_SLOT_SCHEMA_VERSION}.`);
  // Any value and either enablement: readers drop only the lossless legacy null.
  const obsolete = obsoleteCircularSlotField(adv.circular_track_slots);
  if (obsolete) throw diagnosticError('TRACK_INVALID', { ...obsolete, reason: 'OBSOLETE_TRACK_FIELD' });
  validateTrackSlotBindingInvariants(adv.circular_track_slots, { modeLabel: 'Circular', layoutKind: 'circular',
    supportedRenderers: CIRCULAR_TRACK_RENDERERS, supportedSides: ['inside', 'outside', 'overlay'],
    anchorlessRenderers: ['ticks', 'spacer'], depthTrackCount });
};
/**
 * @param {Record<string, any>} [config] The active configuration.
 * @param {{ depthTrackCount?: number | null }} [options]
 * @returns {void}
 */
export const validateImportedLinearTrackSlots = (config = {}, { depthTrackCount = null } = {}) => {
  const adv = isObject(config) ? config.adv : null;
  if (!isObject(adv) || !has(adv, 'linear_track_slots')) return;
  if (![LEGACY_LINEAR_TRACK_SLOT_SCHEMA_VERSION, LINEAR_TRACK_SLOT_SCHEMA_VERSION].includes(adv.linear_track_slots_schema_version))
    throw new Error(`Custom Track Slots use an obsolete schema. Recreate the slots with schema version ${LINEAR_TRACK_SLOT_SCHEMA_VERSION}.`);
  if (!Array.isArray(adv.linear_track_slots)) throw new Error('Custom Track Slots must be an array.');
  validateTrackSlotBindingInvariants(adv.linear_track_slots, { modeLabel: 'Linear', layoutKind: 'linear',
    supportedRenderers: LINEAR_TRACK_RENDERERS, supportedSides: ['above', 'below', 'overlay'],
    anchorlessRenderers: ['spacer'], depthTrackCount });
};
/**
 * Throws when `storedConfig` is not a current-writer active configuration:
 * a Session 46 slice of `mode`, or (`scopedDrafts`) the flat draft of a
 * Session 41-44, whose record display and placement rows name their mode.
 * @param {{ mode: string, storedConfig: Record<string, any> | null | undefined, scopedDrafts?: boolean }} input
 * @returns {void}
 */
export const validateCurrentWriterActiveConfig = ({ mode, storedConfig: config, scopedDrafts = false }) => {
  if (!['circular', 'linear'].includes(mode)) throw new Error(`Current session active configuration has unsupported mode: ${String(mode)}.`);
  if (!isObject(config)) throw new Error('Current session is missing its active Web configuration.');
  assertSafeObjectKeys(config, 'Current session active configuration');
  const domains = new Set([...CURRENT_WRITER_ACTIVE_CONFIG_DOMAINS, 'colorsAreOverrides']);
  const unknownDomains = Object.keys(config).filter((domain) => !domains.has(domain));
  if (unknownDomains.length)
    throw new Error(`Current session active configuration contains unknown domain(s): ${unknownDomains.join(', ')}.`);
  if (!isObject(config.form) || !isObject(config.adv)) throw new Error('Current session is missing its active form or advanced settings.');
  // A list of Linear card uids or Circular source selectors `#N`, once each, as
  // Python reads it; that a Linear uid names a bound card is checked with the
  // bindings (session-authority.js).
  if (has(config, 'recordsOff') && (!Array.isArray(config.recordsOff)
    || new Set(config.recordsOff).size !== config.recordsOff.length
    || !config.recordsOff.every((/** @type {unknown} */ key) => isRecordDrawKey(/** @type {'circular' | 'linear'} */ (mode), key))))
    throw diagnosticError('INPUT_INVALID', { field: 'schema', reason: 'FIELDS' });
  validateDomainShapes(config); validateCollections(config); requireCurrentWebStateFieldNames(config);
  if (has(config, 'recordDisplayDrafts')) validateRecordDisplayDrafts(config.recordDisplayDrafts, { scoped: scopedDrafts });
  if (has(config, 'featurePlacementOverrides')) {
    if (scopedDrafts) validateScopedFeaturePlacements(config.featurePlacementOverrides);
    else canonicalFeaturePlacements(config.featurePlacementOverrides, mode);
  }
  if (has(config.adv, 'feature_overlap_tolerance_bp') && (!Number.isSafeInteger(config.adv.feature_overlap_tolerance_bp)
    || config.adv.feature_overlap_tolerance_bp < 0)) throw new Error('Feature overlap tolerance must be a non-negative integer.');
  assertFields(config.form, new Set(CURRENT_WRITER_FORM_FIELDS), 'config.form'); assertFields(config.adv, new Set(CURRENT_WRITER_ADV_FIELDS), 'config.adv');
  if (has(config.form, 'keep_definition_left_aligned') && typeof config.form.keep_definition_left_aligned !== 'boolean')
    throw new Error('Current session active configuration config.form.keep_definition_left_aligned must be a boolean.');
  if (has(config.form, 'linear_track_layout')) requireCurrentLinearTrackLayout(config.form.linear_track_layout);
  if (has(config.adv, 'label_placement')) requireCurrentLinearLabelPlacement(config.adv.label_placement);
  if (has(config.adv, 'multi_record_size_mode')) requireCurrentCircularMultiRecordSizeMode(config.adv.multi_record_size_mode);
  // A Session 46 slice may omit a field: it reads that mode's default (plan 4.1).
  if (has(config.adv, 'linear_accession_visibility')) requireLinearLabelVisibilityMode(
    config.adv.linear_accession_visibility,
    'Linear Accession visibility'
  );
  if (has(config.adv, 'linear_length_visibility')) requireLinearLabelVisibilityMode(
    config.adv.linear_length_visibility,
    'Linear Length / Coordinates visibility'
  );
  for (const [path, value] of [['config.adv.losatProgram', config.adv.losatProgram], ['config.losatProgram', config.losatProgram]]) {
    if (value !== undefined && !['blastn', 'tblastx', 'blastp'].includes(value))
      throw new Error(`Current session active configuration ${path} is invalid.`);
  }
  if (isObject(config.losat?.blastp)) {
    const blastp = config.losat.blastp;
    if (has(blastp, 'mode')) requireCurrentProteinBlastpMode(blastp.mode);
    if (has(blastp, 'hitLimitsByMode')) {
      if (!isObject(blastp.hitLimitsByMode)) throw new Error('LOSATP mode hit limits must be an object.');
      for (const [mode, limits] of Object.entries(blastp.hitLimitsByMode)) {
        requireCurrentProteinBlastpMode(mode);
        if (!isObject(limits)) throw new Error('LOSATP mode hit limits must be an object.');
        assertFields(limits, new Set(['candidateLimit', 'orthogroupMemberMaxHits']), `losat.blastp.hitLimitsByMode.${mode}`);
        requireCurrentProteinBlastpCandidateLimit(limits.candidateLimit);
        requireCurrentOrthogroupMemberMaxHits(limits.orthogroupMemberMaxHits);
      }
    }
    if (has(blastp, 'candidateLimit')) {
      requireCurrentProteinBlastpCandidateLimit(blastp.candidateLimit);
    }
    requireCurrentProteinBlastpMaxHits(blastp.maxHits);
    requireCurrentOrthogroupMembershipMode(blastp.orthogroupMembershipMode);
    requireCurrentOrthogroupMemberMaxHits(blastp.orthogroupMemberMaxHits);
    requireCurrentCollinearMinAnchors(blastp.collinearMinAnchors);
    requireCurrentCollinearMaxUnitGap(blastp.collinearMaxUnitGap);
    requireCurrentCollinearMaxDiagonalDrift(blastp.collinearMaxDiagonalDrift);
    requireCurrentCollinearMaxConflicts(blastp.collinearMaxConflictsInMergeGap);
    requireCurrentCollinearMaxParalogLinks(
      blastp.collinearMaxParalogLinksPerOrthogroup
    );
    requireCurrentCollinearUnitMode(blastp.collinearUnitMode);
    requireCurrentCollinearAnchorMode(blastp.collinearAnchorMode);
    requireCurrentCollinearMergeOrientation(blastp.collinearMergeOrientation);
    requireCurrentCollinearColorMode(blastp.collinearColorMode);
    requireCurrentCollinearSearchScope(blastp.collinearSearchScope);
    requireCurrentCollinearInferOrthogroups(blastp.collinearInferOrthogroups);
  }
  if (has(config, 'filterMode') && !['None', 'Whitelist', 'Blacklist'].includes(config.filterMode)) throw new Error('Current session active configuration config.filterMode is invalid.');
  if (has(config, 'palette') && !config.palette.trim()) throw new Error('Current session active configuration config.palette cannot be empty.');
  if (isObject(config.importedComparisonResolution)) {
    assertFields(
      config.importedComparisonResolution,
      new Set(['action']),
      'config.importedComparisonResolution'
    );
    if (
      config.importedComparisonResolution.action !== null
      && !['INHERIT', 'REPLACE', 'CLEAR'].includes(config.importedComparisonResolution.action)
    ) {
      throw new Error('Current session active configuration config.importedComparisonResolution.action is invalid.');
    }
  }
  validateImportedCircularTrackSlots(config); validateImportedLinearTrackSlots(config);
};

/** @param {{ scale_font_size?: unknown, ruler_label_font_size?: unknown }} [adv] */
const linearTypographyValuesMatch = (adv = {}) => (
  Object.is(adv.scale_font_size, adv.ruler_label_font_size)
);

/**
 * @param {{
 *   adv: Record<string, any>,
 *   linked: any,
 *   ui?: { linearTypographyLinked?: boolean }
 * }} options
 */
export const reconcileImportedLinearTypographyLink = ({ adv, linked, ui = {} }) => {
  if (!linked || typeof linked !== 'object' || !('value' in linked)) return false;
  // Omission takes the fresh linked default; unequal values still open unlinked.
  linked.value = (
    (ui.linearTypographyLinked ?? true) === true
    && linearTypographyValuesMatch(adv)
  );
  return linked.value;
};

// Artifact metadata admission is shared by Session, History and Align/Reset.
// The binding deliberately excludes current orientation, translations, labels and
// row order: manual Reverse, style and stable reorder do not rewrite history.
/** @param {unknown} value @returns {unknown} */
const sortedValue = (value) => Array.isArray(value)
  ? value.map(sortedValue)
  : isObject(value) ? Object.fromEntries(Object.keys(value).sort()
    .map((key) => [key, sortedValue(value[key])])) : value;

// The parts of a canonical request that the Alignment reset receipt reads. The
// writer (`canonicalSimilarityAlignment` in session-request.js) owns the plan.
/**
 * @typedef {{
 *   records: { recordKey: string, status?: string }[],
 *   reference: { recordKey: string }
 * }} AlignmentPlanView
 */
/**
 * @typedef {{
 *   recordKey: string, source: Record<string, string>, selector?: unknown,
 *   region?: { start: number, end: number, reverseComplement?: unknown } | null,
 *   presentation?: { reverseComplement?: unknown } | null, display?: unknown
 * }} AlignmentResetRecord
 */
/**
 * @typedef {{
 *   mode: string, records: unknown[],
 *   layout: { similarityAlignment?: AlignmentPlanView | null, recordTranslations?: { recordKey: string, x?: number }[], [field: string]: unknown }
 * }} AlignmentResetRequest
 */
/**
 * @typedef {{
 *   renderRequest?: AlignmentResetRequest | null,
 *   resources?: Record<string, RawSessionResourceDescriptor> | null
 * }} AlignmentResetCanonical
 */
/**
 * A canonical session of a Linear Alignment: it holds the plan and the
 * translations that the receipt compares.
 * @typedef {{
 *   renderRequest: {
 *     mode: string, records: AlignmentResetRecord[],
 *     layout: { similarityAlignment: AlignmentPlanView, recordTranslations: { recordKey: string, x?: number }[] }
 *   },
 *   resources?: Record<string, RawSessionResourceDescriptor> | null
 * }} AlignmentResetPlannedCanonical
 */
/**
 * What a reset writes beside the plan: the binding of the plan and the sources,
 * the direction changes, and the reference record's shift (or null).
 * @typedef {{
 *   binding: string,
 *   directions: { recordKey: string, before: boolean, after: boolean }[],
 *   referenceDeltaX: { recordKey: string, deltaX: number } | null
 * }} AlignmentResetReceipt
 */

/**
 * @param {unknown} receipt The stored receipt, unvalidated: this function checks it.
 * @param {AlignmentResetRequest | null | undefined} request
 * @returns {AlignmentResetReceipt | null}
 */
export const validateAlignmentResetReceiptShape = (receipt, request) => {
  if (receipt === null || receipt === undefined) return null;
  // A malformed or stale receipt (OV-115, OV-130) is one user-facing failure.
  /** @type {() => never} */
  const invalid = () => { throw diagnosticError('ALIGNMENT_RESET_EVIDENCE'); };
  if (!isObject(receipt)
    || Object.keys(receipt).sort().join(',') !== 'binding,directions,referenceDeltaX'
    || typeof receipt.binding !== 'string' || !/^[0-9a-f]{64}$/.test(receipt.binding)
    || !Array.isArray(receipt.directions)
    || request?.mode !== 'linear' || !request.layout?.similarityAlignment) invalid();
  const plan = request.layout.similarityAlignment;
  const eligible = new Set(plan.records.filter(({ status }) => status !== 'skipped')
    .map(({ recordKey }) => recordKey));
  /** @type {Set<unknown>} */
  const keys = new Set();
  receipt.directions.forEach((delta) => {
    if (!isObject(delta) || Object.keys(delta).sort().join(',') !== 'after,before,recordKey'
      || !eligible.has(delta.recordKey) || keys.has(delta.recordKey)
      || typeof delta.before !== 'boolean' || typeof delta.after !== 'boolean'
      || delta.before === delta.after) invalid();
    keys.add(delta.recordKey);
  });
  const delta = receipt.referenceDeltaX;
  if (delta !== null && (!isObject(delta)
    || Object.keys(delta).sort().join(',') !== 'deltaX,recordKey'
    || delta.recordKey !== plan.reference.recordKey
    || typeof delta.deltaX !== 'number' || !Number.isFinite(delta.deltaX)
    || delta.deltaX === 0)) invalid();
  // Every field was checked above.
  return /** @type {AlignmentResetReceipt} */ (receipt);
};

/** @param {AlignmentResetCanonical} canonical */
const resetBinding = async (canonical) => {
  const request = canonical?.renderRequest;
  if (request?.mode !== 'linear' || !request.layout?.similarityAlignment) {
    throw new Error('Alignment reset receipt requires an active canonical plan.');
  }
  const table = adoptCurrentSessionResources(canonical.resources);
  /** @type {Map<string, Promise<string>>} */
  const fingerprints = new Map();
  /** @param {Record<string, string>} source */
  const sourceIdentity = async (source) => {
    /** @type {Record<string, string | undefined>} */
    const identity = {};
    for (const [key, value] of Object.entries(source)) {
      if (key === 'kind') identity[key] = value;
      else {
        if (!fingerprints.has(value)) fingerprints.set(value, sha256Hex(
          await readSessionResourceBytes(createSessionResourceFileView(table, value))
        ));
        identity[key] = await fingerprints.get(value);
      }
    }
    return identity;
  };
  // The mode check above admits a Linear request only; the request writer gives each of its records a `recordKey` and a `source`.
  const records = await Promise.all(/** @type {AlignmentResetRecord[]} */ (request.records).map(async (record) => ({
    recordKey: record.recordKey, source: await sourceIdentity(record.source),
    selector: record.selector,
    region: record.region ? { start: record.region.start, end: record.region.end } : null,
    display: record.display
  })));
  const plan = { ...request.layout.similarityAlignment,
    records: [...request.layout.similarityAlignment.records].sort(
      (a, b) => a.recordKey < b.recordKey ? -1 : a.recordKey > b.recordKey ? 1 : 0) };
  records.sort((a, b) => a.recordKey < b.recordKey ? -1 : a.recordKey > b.recordKey ? 1 : 0);
  return sha256Hex(textToBytes(JSON.stringify(sortedValue({ plan, records }))));
};

/**
 * @param {unknown} receipt The stored receipt, unvalidated.
 * @param {AlignmentResetCanonical | null | undefined} canonical
 */
export const validateSimilarityAlignmentResetReceipt = async (receipt, canonical) => {
  const validated = validateAlignmentResetReceiptShape(receipt, canonical?.renderRequest);
  // A receipt that passed the shape check needs `request.mode === 'linear'`, so `canonical` is set.
  if (validated && validated.binding !== await resetBinding(/** @type {AlignmentResetCanonical} */ (canonical))) {
    throw diagnosticError('ALIGNMENT_RESET_EVIDENCE');
  }
  return validated;
};

/**
 * @param {{ before: AlignmentResetPlannedCanonical, after: AlignmentResetPlannedCanonical }} sessions
 */
export const buildSimilarityAlignmentResetReceipt = async ({ before, after }) => {
  const previous = new Map(before.renderRequest.records.map((record) => [record.recordKey, record]));
  const directions = after.renderRequest.records.flatMap((record) => {
    const old = previous.get(record.recordKey);
    if (!old) throw new Error('Alignment reset receipt record binding changed.');
    const beforeReverse = canonicalRecordReverseComplement(old);
    const afterReverse = canonicalRecordReverseComplement(record);
    return beforeReverse === afterReverse ? []
      : [{ recordKey: record.recordKey, before: beforeReverse, after: afterReverse }];
  });
  const referenceKey = after.renderRequest.layout.similarityAlignment.reference.recordKey;
  /** @param {AlignmentResetPlannedCanonical} canonical */
  const x = (canonical) => canonical.renderRequest.layout.recordTranslations
    .find(({ recordKey }) => recordKey === referenceKey)?.x ?? 0;
  const deltaX = x(after) - x(before);
  const receipt = { binding: await resetBinding(after), directions,
    referenceDeltaX: deltaX === 0 ? null : { recordKey: referenceKey, deltaX } };
  return validateSimilarityAlignmentResetReceipt(receipt, after);
};
