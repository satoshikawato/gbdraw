// G-G(4) / R7 (X-02: GE-04, GE-08, TR-04, CO-09): every numeric draft value the
// Web projects into the canonical request is either rejected with a typed
// diagnostic or passed literally so Python judges it. No value is silently
// replaced by a default, and the shared option-domain vectors reach the same
// field and reason as the CLI and the typed request.
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import test from 'node:test';
import { buildCanonicalRenderRequest } from '../../gbdraw/web/js/services/session-request.js';
import { createDefaultLinearComparisonPlan } from '../../gbdraw/web/js/services/linear-comparisons.js';
import { resolveComparisonThresholds } from '../../gbdraw/web/js/mode-profiles.js';
import { normalizeUserFacingError } from '../../gbdraw/web/js/utils/error-normalization.js';

const ref = (value) => ({ value });
const genbankText = 'LOCUS       WEBTEST                    4 bp    DNA     linear   UNK 01-JAN-1980\nORIGIN\n        1 atgc\n//\n';
const genbank = { name: 'input.gb', type: 'text/plain', size: genbankText.length, lastModified: 0, data: btoa(genbankText) };
const baseState = (mode, adv = {}) => ({
  mode: ref(mode),
  cInputType: ref('gb'),
  lInputType: ref('gb'),
  circularRecordList: ref([]),
  form: {
    prefix: 'web-session', species: '', strain: '', plot_title: '', legend: 'right',
    multi_record_canvas: true, suppress_gc: false, suppress_skew: false,
    show_gc: false, show_skew: false, show_depth: false, separate_strands: true,
    labels_mode: 'none', show_labels_linear: 'none', track_type: 'tuckin',
    linear_track_layout: 'middle', linear_ruler_on_axis: false, align_center: false,
    keep_definition_left_aligned: false, show_scale: true, scale_style: 'bar',
    normalize_length: false
  },
  adv: {
    features: ['CDS'], feature_shapes: { CDS: 'arrow' }, nt: 'GC', evalue: '1e-5',
    arrow_head_length_ratio: null, arrow_shaft_width_ratio: 1.0,
    min_bitscore: 50, identity: 70, alignment_length: 0, plot_title_position: 'none',
    gc_content_mode: 'deviation', gc_content_min_percent: null, gc_content_max_percent: null,
    gc_content_show_axis: true, gc_content_show_ticks: true,
    depth_color: '#4A90E2', depth_show_axis: true, depth_show_ticks: true,
    pairwise_match_style: 'ribbon', multi_record_size_mode: 'auto',
    multi_record_min_radius_ratio: null, multi_record_column_gap_ratio: null,
    multi_record_row_gap_ratio: null, multi_record_positions: [], depth_tracks: [],
    circular_track_slots_enabled: false, circular_track_slots_axis_index: null,
    circular_track_slots: [], linear_track_slots_enabled: false,
    linear_track_slots_axis_index: null, linear_track_slots: [],
    linear_accession_visibility: 'auto', linear_length_visibility: 'auto', label_placement: 'auto',
    ...adv
  },
  normalizePaletteColors: (value) => value,
  paletteDefinitions: ref({ default: {} }),
  currentColors: ref({}),
  selectedPalette: ref('default'),
  manualSpecificRules: [],
  featureVisibilityRules: ref([]),
  filterMode: ref('None'),
  manualBlacklist: ref(''),
  manualWhitelist: [],
  manualPriorityRules: [],
  labelTextFeatureOverrides: {},
  labelTextBulkOverrides: {},
  labelTextFeatureOverrideSources: {},
  labelVisibilityOverrides: {},
  editableLabels: ref([]),
  extractedFeatures: ref([]),
  circularConservation: { reference: 'auto', labels: '', series: [] },
  linearComparisonPlan: createDefaultLinearComparisonPlan(),
  losatProgram: ref('blastn'),
  losat: { blastp: {} },
  selectedOrthogroupAlignmentFeature: ref(''),
  similarityAlignmentPlan: ref(null),
  linearRecordTranslations: ref([]),
  legacySimilarityAlignment: ref(null),
  linearRecordLayoutEnabled: ref(false),
  linearRecordGap: ref(24),
  linearRecordRows: [],
  annotationSets: [],
  unmanagedConfigOverrides: {}
});
const filesFor = (mode) => (mode === 'circular'
  ? { c_gb: genbank, linearSeqs: [] }
  : { linearSeqs: [{ uid: 'first', gb: genbank, losat_gencode: 1, region_record_id: '', region_start: null, region_end: null, region_reverse: false }] });
const comparisonPlanSnapshot = { hasComparisonIntent: true, hasLosatIntent: false, edges: [] };
const project = (mode, adv, state = baseState(mode, adv)) => buildCanonicalRenderRequest({
  state,
  drawing: state,
  filesData: filesFor(mode),
  ...(mode === 'linear' ? { comparisonPlanSnapshot } : {})
}).renderRequest;

// Leaves of the request whose values differ from the blank-draft baseline.
const changedLeaves = (baseline, candidate, path = '') => {
  if (candidate && typeof candidate === 'object') {
    return Object.keys(candidate).flatMap((key) => changedLeaves(baseline?.[key], candidate[key], `${path}/${key}`));
  }
  return Object.is(baseline, candidate) ? [] : [[path, candidate]];
};

// Every numeric draft field that the canonical request projects, per mode.
const NUMERIC_FIELDS = {
  circular: ['window_size', 'step_size', 'depth_window_size', 'depth_step_size', 'plot_title_font_size',
    'gc_content_min_percent', 'gc_content_max_percent', 'gc_content_tick_interval', 'gc_content_small_tick_interval',
    'gc_content_tick_font_size', 'depth_min', 'depth_max', 'depth_large_tick_interval', 'depth_small_tick_interval',
    'depth_tick_font_size', 'scale_interval', 'def_font_size', 'circular_definition_interval', 'circular_label_spacing',
    'tick_label_font_size', 'outer_label_x_offset', 'outer_label_y_offset', 'inner_label_x_offset', 'inner_label_y_offset',
    'block_stroke_width', 'line_stroke_width', 'legend_box_size', 'legend_font_size', 'axis_stroke_width', 'label_font_size',
    'center_reserved_radius', 'multi_record_min_radius_ratio', 'multi_record_column_gap_ratio', 'multi_record_row_gap_ratio',
    'min_bitscore', 'evalue', 'identity', 'alignment_length'],
  linear: ['window_size', 'step_size', 'depth_window_size', 'depth_step_size', 'plot_title_font_size', 'scale_interval',
    'linear_label_spacing', 'label_rotation', 'track_axis_gap', 'gc_height', 'depth_height', 'scale_stroke_width',
    'block_stroke_width', 'line_stroke_width', 'legend_box_size', 'legend_font_size', 'axis_stroke_width', 'def_font_size',
    'feature_height', 'scale_font_size', 'ruler_label_font_size', 'label_font_size',
    'min_bitscore', 'evalue', 'identity', 'alignment_length']
};
const BLANK = { min_bitscore: null, evalue: null, identity: null, alignment_length: null };
// A valid, distinctive value: a projection that substitutes a default for the
// candidate leaves the sentinel's leaves unchanged or different from the candidate.
const SENTINEL = 7;

test('every projected numeric draft value is rejected with a diagnostic or sent literally (G-G(4))', () => {
  for (const [mode, fields] of Object.entries(NUMERIC_FIELDS)) {
    for (const field of fields) {
      const baseline = project(mode, { ...BLANK, [field]: SENTINEL });
      for (const value of [Number.NaN, '1e-50x', -5, 0, 12.5]) {
        const label = `${mode}.${field}=${String(value)}`;
        let request;
        try {
          request = project(mode, { ...BLANK, [field]: value });
        } catch (error) {
          const model = normalizeUserFacingError(error);
          assert.equal(model.code, 'INPUT_INVALID', label);
          assert.ok(model.context.field || model.context.configPath, label);
          assert.ok(model.context.reason, label);
          continue;
        }
        assert.ok(typeof value === 'number' && Number.isFinite(value), `${label} must not pass silently`);
        const leaves = changedLeaves(baseline, request);
        assert.ok(leaves.length > 0, `${label} was replaced by the default`);
        for (const [path, projected] of leaves) assert.equal(projected, value, `${label} at ${path}`);
        assert.doesNotMatch(JSON.stringify(request), /NaN|Infinity/, label);
      }
    }
  }
});

test('the shared option-domain vectors reach the CLI field and reason through the Web projection', () => {
  const vectors = JSON.parse(readFileSync(new URL('../fixtures/option_domain_vectors.json', import.meta.url), 'utf8'));
  let checked = 0;
  for (const kind of ['invalid', 'accepted']) {
    for (const vector of vectors[kind].filter((entry) => entry.web)) {
      for (const mode of vector.modes) {
        const label = `${vector.id}-${mode}`;
        let request = null;
        let model = null;
        try {
          request = project(mode, vector.web.adv);
        } catch (error) {
          model = normalizeUserFacingError(error);
        }
        checked += 1;
        if (model) {
          // The Web evaluates the generated domain before Python (thresholds).
          assert.equal(kind, 'invalid', `${label} rejected: ${model.code}`);
          assert.deepEqual([model.code, model.context.field, model.context.reason],
            [vector.code, vector.field, vector.reason], label);
          continue;
        }
        // Otherwise the request carries exactly what the CLI/typed request vector judges.
        const expected = vector.web.requestByMode?.[mode] || vector.web.request || vector.request;
        for (const [key, value] of Object.entries(expected)) {
          if (key === 'configOverrides') {
            for (const [path, leaf] of Object.entries(value)) {
              assert.equal(request.diagramOptions.configOverrides[path], leaf, `${label} ${path}`);
            }
          } else {
            assert.equal(request.diagramOptions[key], value, `${label} ${key}`);
          }
        }
      }
    }
  }
  assert.ok(checked >= 40, checked);
});

test('comparison thresholds resolve once on the generated domains without rewriting the draft', () => {
  for (const mode of ['circular', 'linear']) {
    const adv = { min_bitscore: '', evalue: ' 1e-10 ', identity: '35', alignment_length: 12 };
    const before = structuredClone(adv);
    assert.deepEqual(resolveComparisonThresholds(adv, mode), {
      bitscore: mode === 'circular' ? 50 : 50, evalue: '1e-10', identity: 35, alignmentLength: 12
    });
    assert.deepEqual(adv, before, 'the draft keeps what was typed');
    for (const [field, value, reason, key] of [
      ['evalue', '1e-50x', 'NONNEGATIVE', 'evalue'], ['evalue', 'NaN', 'NONNEGATIVE', 'evalue'],
      ['identity', -5, 'PERCENT', 'identity'], ['identity', 150, 'PERCENT', 'identity'],
      ['min_bitscore', -3, 'NONNEGATIVE', 'bitscore'],
      ['alignment_length', 150.5, 'NONNEGATIVE_INTEGER', 'alignment_length'],
      ['alignment_length', -10, 'NONNEGATIVE_INTEGER', 'alignment_length']
    ]) {
      assert.throws(() => resolveComparisonThresholds({ [field]: value }, mode),
        { code: 'INPUT_INVALID', context: { field: key, reason } }, `${mode} ${field}=${value}`);
    }
  }
  const linear = project('linear', { min_bitscore: '', evalue: '', identity: '', alignment_length: '' }).diagramOptions;
  assert.deepEqual([linear.evalue, linear.bitscore, linear.identity, linear.alignmentLength], [0.01, 50, 0, 0]);
});

test('Generate keeps the numeric draft: no Generate-path draft assignments are added (R7 ratchet)', () => {
  // Generate-path draft assignments in run-analysis.js may only shrink (69 on
  // dev before X-02). Lower the baseline when one is removed.
  const BASELINE = 44;
  const source = readFileSync(new URL('../../gbdraw/web/js/app/run-analysis.js', import.meta.url), 'utf8');
  const assignments = source.split('\n').filter((line) => (
    /^\s*(drawing\.)?(adv|form|circularConservation|losat(\.[a-z]+)?)\.[a-zA-Z_]+(\[[^\]]*\])?\s*=[^=]/.test(line)
    || /(adv|form|circularConservation)\.[a-z_]+\.splice\(/.test(line)
  ));
  assert.equal(assignments.length, BASELINE, assignments.join('\n'));
  for (const field of ['min_bitscore', 'evalue', 'identity', 'alignment_length', 'plot_title_font_size',
    'circular_label_spacing', 'linear_label_spacing', 'center_reserved_radius', 'multi_record_min_radius_ratio',
    'multi_record_column_gap_ratio', 'multi_record_row_gap_ratio', 'depth_height', 'ring_width', 'ring_gap', 'subject_gencode']) {
    assert.doesNotMatch(source, new RegExp(`(adv|circularConservation)\\.${field}\\s*=[^=]`), field);
  }
});
