import { createDefaultLinearDefinitionLineStyles } from '../app/definition-line-style-state.js'; import { CIRCULAR_TRACK_RENDERERS, createDefaultCircularTrackSlots } from '../app/circular-track-slots.js';
import { LEGACY_LINEAR_TRACK_SLOT_SCHEMA_VERSION, LINEAR_TRACK_RENDERERS, LINEAR_TRACK_SLOT_SCHEMA_VERSION, createDefaultLinearTrackSlots } from '../app/linear-track-slots.js'; import { validateTrackSlotBindingInvariants } from '../app/track-slot-validation.js';
import { requireCurrentCircularMultiRecordSizeMode, requireCurrentCollinearAnchorMode, requireCurrentCollinearColorMode, requireCurrentCollinearMaxConflicts, requireCurrentCollinearMaxDiagonalDrift, requireCurrentCollinearMaxParalogLinks, requireCurrentCollinearMaxUnitGap, requireCurrentCollinearMergeOrientation, requireCurrentCollinearMinAnchors, requireCurrentCollinearSearchScope, requireCurrentCollinearUnitMode, requireCurrentLinearLabelPlacement, requireCurrentLinearTrackLayout, requireCurrentOrthogroupMemberMaxHits, requireCurrentOrthogroupMembershipMode, requireCurrentProteinBlastpCandidateLimit, requireCurrentProteinBlastpMaxHits, requireCurrentProteinBlastpMode, requireCurrentWebStateFieldNames } from '../app/current-option-values.js'; import { DEFAULT_ARROW_SHAFT_WIDTH_RATIO, createDefaultFeatureRenderings } from '../utils/feature-rendering.js';
import { MODE_DEFAULT_FEATURE_TYPES, comparisonStateForMode, managedAdvStateForMode, trackDefaultsForMode } from '../mode-profiles.js'; import { WEB_UX_PROFILE } from '../web-ux-profile.js';
import { assertSafeObjectKeys } from './safe-object-keys.js';
const circularTracks = trackDefaultsForMode('circular'), linearTracks = trackDefaultsForMode('linear');
export const CIRCULAR_TRACK_SLOT_SCHEMA_VERSION = 4, LEGACY_CIRCULAR_TRACK_SLOT_SCHEMA_VERSION = 3;
const isObject = (value) => Boolean(value) && typeof value === 'object' && !Array.isArray(value), has = (value, key) => Object.prototype.hasOwnProperty.call(value, key);
export const createDefaultForm = () => ({
  prefix: '', species: '', strain: '', plot_title: '', track_type: 'tuckin', linear_track_layout: 'middle', show_scale: true, scale_style: 'bar',
  linear_ruler_on_axis: false, labels_mode: 'none', show_labels_linear: 'none', multi_record_canvas: WEB_UX_PROFILE.circular.gridByDefault,
  circular_record_selector: '', circular_region_start: null, circular_region_end: null, circular_reverse: false,
  circular_record_label: '', circular_record_subtitle: '',
  separate_strands: WEB_UX_PROFILE.separateStrands, suppress_gc: !circularTracks.gc, suppress_skew: !circularTracks.skew, align_center: false,
  keep_definition_left_aligned: false, show_gc: linearTracks.gc, show_skew: linearTracks.skew, show_depth: false, normalize_length: false
});
export const createDefaultAdv = (mode = 'circular') => ({
  rich_feature_popup: true, features: [...MODE_DEFAULT_FEATURE_TYPES], feature_shapes: createDefaultFeatureRenderings(), arrow_head_length_ratio: null,
  arrow_shaft_width_ratio: DEFAULT_ARROW_SHAFT_WIDTH_RATIO, window_size: null, step_size: null, nt: 'GC', def_font_size: null,
  circular_definition_interval: null, label_font_size: null, circular_label_spacing: null, linear_label_spacing: null, label_rendering: 'auto',
  circular_label_placement: 'horizontal', label_placement: 'auto', label_rotation: null, block_stroke_width: null, block_stroke_color: null,
  line_stroke_width: null, line_stroke_color: null, axis_stroke_width: null, axis_stroke_color: managedAdvStateForMode(mode).axis_stroke_color,
  legend_box_size: null, legend_font_size: null, resolve_overlaps: false, feature_height: null, track_axis_gap: null, linear_show_replicon: false,
  linear_show_accession: true, linear_show_length: true, linear_definition_line_styles: createDefaultLinearDefinitionLineStyles(), gc_height: null,
  depth_height: null, depth_color: '#4A90E2', depth_tracks: [], depth_window_size: null, depth_step_size: null, depth_share_axis: false,
  depth_min: null, depth_max: null, depth_normalize: false, depth_show_axis: true, depth_show_ticks: true, depth_large_tick_interval: null,
  depth_small_tick_interval: null, depth_tick_font_size: null, linear_track_slots_enabled: false, linear_track_slots_schema_version: LINEAR_TRACK_SLOT_SCHEMA_VERSION,
  linear_track_slots_axis_index: null, linear_track_slots: createDefaultLinearTrackSlots(), gc_content_mode: 'deviation', gc_content_min_percent: 0,
  gc_content_max_percent: 100, gc_content_show_axis: true, gc_content_show_ticks: true, gc_content_tick_interval: 20, gc_content_small_tick_interval: null,
  gc_content_tick_font_size: null, comparison_height: null, pairwise_match_style: 'ribbon', ...comparisonStateForMode(mode), scale_interval: null,
  scale_font_size: null, ruler_label_font_size: null, scale_stroke_width: null, scale_stroke_color: null, ruler_label_color: null, circular_grouping_intent: 'auto',
  multi_record_size_mode: 'auto', multi_record_min_radius_ratio: 0.55, multi_record_column_gap_ratio: 0.10, multi_record_row_gap_ratio: 0.05,
  multi_record_positions: [], tick_label_font_size: null, plot_title_font_size: null, keep_full_definition_with_plot_title: false,
  center_reserved_radius: null, feature_width_circular: null, depth_width_circular: null, gc_content_width_circular: null, gc_content_radius_circular: null,
  gc_skew_width_circular: null, gc_skew_radius_circular: null, circular_track_slots_enabled: false,
  circular_track_slots_schema_version: CIRCULAR_TRACK_SLOT_SCHEMA_VERSION, circular_track_slots_axis_index: null, circular_track_slots: createDefaultCircularTrackSlots(),
  outer_label_x_offset: null, outer_label_y_offset: null, inner_label_x_offset: null, inner_label_y_offset: null
});
export const createDefaultLosat = () => ({
  outfmt: '6', parallelWorkers: undefined, executionMode: 'auto', totalThreadBudget: 'safe', threadsPerJob: 'auto',
  blastn: { task: 'megablast' }, blastp: { mode: 'orthogroup', maxHits: 5, candidateLimit: null, orthogroupMembershipMode: 'anchor_core_v1',
    orthogroupMemberMaxHits: 5, collinearMinAnchors: 1, collinearMaxUnitGap: 0, collinearMaxDiagonalDrift: 0,
    collinearMaxConflictsInMergeGap: 1, collinearMaxParalogLinksPerOrthogroup: 2, collinearColorMode: 'orientation',
    collinearUnitMode: 'auto', collinearAnchorMode: 'rbh', collinearMergeOrientation: 'either', collinearSearchScope: 'adjacent' }
});
export const createDefaultCircularConservation = () => ({ enabled: false, source: 'losat', losat_program: 'blastn',
  subject_gencode: 1, reference: 'auto', labels: '', series: [], ring_width: null, ring_gap: null });
export const CURRENT_WRITER_ACTIVE_CONFIG_DOMAINS = Object.freeze([
  'form', 'adv', 'losat', 'cliOptions', 'colors', 'palette', 'paletteInstantPreviewEnabled', 'rules',
  'qualifierPriorityRules', 'filterMode', 'whitelist', 'blacklistText', 'losatProgram', 'circularConservation',
  'annotationSets', 'modeProfiles', 'unmanagedConfigOverrides',
  'linearRecordLayout', 'linearComparisonPlan', 'importedComparisonResolution', 'webEdits'
]);
export const CURRENT_WRITER_FORM_FIELDS = Object.freeze([...Object.keys(createDefaultForm()), 'legend']);
export const CURRENT_WRITER_ADV_FIELDS = Object.freeze([...Object.keys(createDefaultAdv()), 'plot_title_position', 'losatProgram']);
const DOMAIN_SHAPES = Object.freeze({ form: 'object', adv: 'object', losat: 'object', cliOptions: 'object', colors: 'object',
  circularConservation: 'object', modeProfiles: 'object', linearRecordLayout: 'object', linearComparisonPlan: 'object', webEdits: 'object',
  importedComparisonResolution: 'object', unmanagedConfigOverrides: 'object',
  rules: 'array', qualifierPriorityRules: 'array', whitelist: 'array', annotationSets: 'array', palette: 'string', filterMode: 'string', blacklistText: 'string',
  losatProgram: 'string', paletteInstantPreviewEnabled: 'boolean' });
const ROW_FIELDS = { rules: ['feat', 'qual', 'val', 'color', 'cap', 'fromFile'],
  qualifierPriorityRules: ['feat', 'order'], whitelist: ['feat', 'qual', 'key'],
  annotationSets: ['id', 'annotations', 'defaultStyle', 'legendLabel'] };
const ANNOTATION_FIELDS = ['id', 'target', 'label', 'mark', 'lane', 'style', 'legendLabel', 'metadata'];
const OBSOLETE_SLOT_KEYS = new Set(['gapAfter', 'gap_after', 'innerRadius', 'inner_radius', 'outerRadius',
  'outer_radius', 'placement', 'spacing', 'strict', 'compress', 'reserve']);
const OBSOLETE_SLOT_PARAM_KEYS = new Set(['side', 'radius', 'width', 'spacing', 'inner_gap_px', 'outer_gap_px', 'strict', 'compress', 'reserve']);
const assertFields = (value, fields, path) => {
  const unknown = Object.keys(value).filter((key) => !fields.has(key));
  if (unknown.length) throw new Error(`Current session active configuration contains unknown ${path} field(s): ${unknown.join(', ')}.`);
};
const validateDomainShapes = (config) => {
  for (const [domain, shape] of Object.entries(DOMAIN_SHAPES)) {
    if (!has(config, domain) || (config[domain] === undefined && domain === 'cliOptions')) continue;
    const value = config[domain], valid = shape === 'object' ? isObject(value) : shape === 'array' ? Array.isArray(value) : typeof value === shape;
    if (!valid) throw new Error(`Current session active configuration config.${domain} must be ${shape}.`);
  }
  if (has(config, 'colorsAreOverrides') && typeof config.colorsAreOverrides !== 'boolean')
    throw new Error('Current session compatibility field config.colorsAreOverrides must be a boolean.');
};
const assertRows = (rows, fields, path) => rows.forEach((row, index) => {
  if (!isObject(row)) throw new Error(`Current session active configuration ${path}[${index}] must be an object.`);
  assertFields(row, new Set(fields), `${path}[${index}]`);
});
const validateCollections = (config) => {
  for (const [domain, fields] of Object.entries(ROW_FIELDS))
    if (Array.isArray(config[domain])) assertRows(config[domain], fields, `config.${domain}`);
  (config.annotationSets || []).forEach((set, index) => {
    if (!Array.isArray(set.annotations)) throw new Error(`Current session active configuration config.annotationSets[${index}].annotations must be an array.`);
    assertRows(set.annotations, ANNOTATION_FIELDS, `config.annotationSets[${index}].annotations`);
  });
};
const obsoleteCircularSlotPath = (slots) => {
  if (!Array.isArray(slots)) return null;
  for (const [index, slot] of slots.entries()) {
    if (!isObject(slot)) continue;
    const field = [...OBSOLETE_SLOT_KEYS].find((key) => has(slot, key));
    if (field) return `circular_track_slots[${index}].${field}`;
    if (!isObject(slot.params)) continue;
    const param = [...OBSOLETE_SLOT_PARAM_KEYS].find((key) => has(slot.params, key));
    if (param) return `circular_track_slots[${index}].params.${param}`;
  }
  return null;
};
export const validateImportedCircularTrackSlots = (config = {}, { depthTrackCount = null } = {}) => {
  const adv = isObject(config) ? config.adv : null;
  if (!isObject(adv) || !has(adv, 'circular_track_slots')) return;
  if (adv.circular_track_slots_schema_version !== CIRCULAR_TRACK_SLOT_SCHEMA_VERSION)
    throw new Error(`Custom Track Slots use an obsolete schema. Recreate the slots with schema version ${CIRCULAR_TRACK_SLOT_SCHEMA_VERSION}.`);
  const obsolete = obsoleteCircularSlotPath(adv.circular_track_slots);
  if (obsolete) throw new Error(`Custom Track Slots use obsolete field '${obsolete}'. Use slot-level radius, width, inner_gap_px, outer_gap_px, side, and z fields.`);
  validateTrackSlotBindingInvariants(adv.circular_track_slots, { modeLabel: 'Circular', layoutKind: 'circular',
    supportedRenderers: CIRCULAR_TRACK_RENDERERS, supportedSides: ['inside', 'outside', 'overlay'],
    anchorlessRenderers: ['ticks', 'spacer'], depthTrackCount });
};
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
export const validateCurrentWriterActiveConfig = ({ mode, storedConfig: config }) => {
  if (!['circular', 'linear'].includes(mode)) throw new Error(`Current session active configuration has unsupported mode: ${String(mode)}.`);
  if (!isObject(config)) throw new Error('Current session is missing its active Web configuration.');
  assertSafeObjectKeys(config, 'Current session active configuration');
  const domains = new Set([...CURRENT_WRITER_ACTIVE_CONFIG_DOMAINS, 'colorsAreOverrides']);
  const unknownDomains = Object.keys(config).filter((domain) => !domains.has(domain));
  if (unknownDomains.length)
    throw new Error(`Current session active configuration contains unknown domain(s): ${unknownDomains.join(', ')}.`);
  if (!isObject(config.form) || !isObject(config.adv)) throw new Error('Current session is missing its active form or advanced settings.');
  validateDomainShapes(config); validateCollections(config); requireCurrentWebStateFieldNames(config);
  assertFields(config.form, new Set(CURRENT_WRITER_FORM_FIELDS), 'config.form'); assertFields(config.adv, new Set(CURRENT_WRITER_ADV_FIELDS), 'config.adv');
  if (has(config.form, 'linear_track_layout')) requireCurrentLinearTrackLayout(config.form.linear_track_layout);
  if (has(config.adv, 'label_placement')) requireCurrentLinearLabelPlacement(config.adv.label_placement);
  if (has(config.adv, 'multi_record_size_mode')) requireCurrentCircularMultiRecordSizeMode(config.adv.multi_record_size_mode);
  for (const [path, value] of [['config.adv.losatProgram', config.adv.losatProgram], ['config.losatProgram', config.losatProgram]]) {
    if (value !== undefined && !['blastn', 'tblastx', 'blastp'].includes(value))
      throw new Error(`Current session active configuration ${path} is invalid.`);
  }
  if (isObject(config.losat?.blastp)) {
    const blastp = config.losat.blastp;
    if (has(blastp, 'mode')) requireCurrentProteinBlastpMode(blastp.mode);
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
