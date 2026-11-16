// One non-default sample per member of a drawing (the members of
// `DRAWING_SLICES`, helpers/drawing-state.mjs), in the form the member's slice
// `capture` returns and `install` takes. A sample is a value the product
// itself writes and keeps: a fixed point of Session Load, so a test that saves
// a drawing holding it and loads the Session finds it unchanged. Where Load
// normalizes a value, the sample is already in Load's form, and its comment
// names the normalization.
//
// A sample is a function of (mode, drawing, context):
// - `mode` is the drawing's mode ('circular' or 'linear');
// - `drawing` is the drawing the sample goes into, read for the member's
//   current value (a sample of an object member changes one field of it) and
//   for the records it binds;
// - `context` holds `registerFile(file)`, which registers a File in the
//   project's file store (`state.inputFiles`) and returns its binding.
// A member whose sample needs more than a new drawing holds lists it in
// `needs`: 'records' (the drawing binds a GenBank file of two records, whose
// cards or rows the sample names), 'groups' (the drawing holds a Result with
// Similarity groups, whose ids the sample names), or 'comparison' (the
// drawing holds a Result with a comparison read from a CLI table). A context
// inside the drawing's own slices is `prepare(drawing)` (the queued palette a
// queued color belongs to). The test that installs the sample builds that
// context first, as the user would.
//
// A member that only an action writes (`byContext`: the records a file read
// discovers, the answer to a saved comparison) has no value of its own: the
// test takes the action named by `needs` as the edit, and its result is the
// sample. `value` gives the member as that action leaves it in `drawing`.
//
// `modes` lists the modes in which the member alone is a change Save writes
// (both by default): Save counts a setting of one mode, such as the Linear
// comparison plan, only in a drawing of that mode (the `modes` of its row in
// the mode-scoped registry, gbdraw/web_support/mode_scoped_settings.py), and
// the record rows and File cards only in the mode that draws them. In a
// drawing of the other mode, written for another change, the member is kept
// through every path, Save then Load included, unless `otherModeLoad` is
// 'default': that drawing's Session does not hold it, and its Load gives the
// mode's default. tests/web/drawing-member-paths.test.mjs runs both.
//
// Import this module after the test has set up `window`.
import { DRAWING_SLICES } from './drawing-state.mjs';
import { normalizeAnnotationSets } from '../../../gbdraw/web/js/services/annotation-state.js';
import { normalizeFeatureVisibilityRule } from '../../../gbdraw/web/js/services/feature-visibility.js';
import { normalizeSpecificRule } from '../../../gbdraw/web/js/services/specific-color-rules.js';
import { normalizePaletteColors } from '../../../gbdraw/web/js/utils/color-utils.js';
import { createCircularSource, createLinearSeq } from '../../../gbdraw/web/js/services/drawing/inputs.js';
import { queuePalette } from '../../../gbdraw/web/js/services/drawing/palette.js';
import { expandedRecordKey } from '../../../gbdraw/web/js/services/record-draw-selection.js';

/** @typedef {'circular' | 'linear'} Mode */
/** @typedef {import('../../../gbdraw/web/js/state.js').DrawingState} DrawingState */
/** @typedef {{ registerFile: (file: File) => unknown }} SampleContext */
/**
 * @typedef {object} MemberSample
 * @property {(mode: Mode, drawing: DrawingState, context: SampleContext) => unknown} value
 * @property {readonly string[]} [needs]
 * @property {(drawing: DrawingState) => void} [prepare]
 * @property {readonly Mode[]} [modes]
 * @property {'default'} [otherModeLoad]
 * @property {boolean} [byContext]
 */

// The slice that holds `member`, and the member's current value in `drawing`.
const SLICE_OF = new Map(Object.entries(DRAWING_SLICES).flatMap(([name, slice]) => slice.members.map((member) => [member, name])));
export const sliceOfMember = (member) => SLICE_OF.get(member);
export const captureMember = (drawing, member) => {
  const name = SLICE_OF.get(member);
  return DRAWING_SLICES[name].capture(drawing.slices[name])[member];
};

// The records a 'records' sample names: the drawing's bound cards (Circular)
// or record rows (Linear), by uid.
const boundUids = (drawing, mode) => (mode === 'linear'
  ? captureMember(drawing, 'linearSeqs').filter((row) => row.gb).map((row) => row.uid)
  : captureMember(drawing, 'circularSources').filter((card) => card.gb).map((card) => card.uid));

// A GenBank file of `ids`, one 40 bp record each, with one CDS.
export const genbankFile = (name, ids) => new File([ids.map((id) => [
  `LOCUS       ${id} 40 bp DNA circular`,
  'FEATURES             Location/Qualifiers',
  '     CDS             1..30',
  `                     /locus_tag="${id}_1"`,
  `ORIGIN\n        1 ${'a'.repeat(40)}\n//\n`
].join('\n')).join('')], name, { lastModified: 1 });

/** @type {Record<string, MemberSample>} */
export const MEMBER_SAMPLES = {
  // visibility rules
  featureVisibilityManualRules: { value: () => [normalizeFeatureVisibilityRule({
    id: 'sample-rule', featureType: 'CDS', qualifier: 'locus_tag', value: '^A$', action: 'off'
  })] },
  // color rules
  manualSpecificRules: { value: () => [normalizeSpecificRule({
    feat: 'CDS', qual: 'product', val: 'kinase', color: '#123456', cap: 'Kinase', fromFile: true
  })] },
  // palette: Load gives a palette's colors the comparison colors (`normalizePaletteColors`).
  selectedPalette: { value: () => 'mint' },
  currentColors: { value: () => normalizePaletteColors({ CDS: '#123456', tRNA: '#654321' }) },
  pendingPaletteName: { value: () => 'queued' },
  // A palette queued while Instant Preview is off holds its name and the
  // colors it was queued with (`queuePalette`); the sample changes one.
  pendingPaletteColors: { prepare: (drawing) => queuePalette(drawing.slices.palette, 'queued', drawing.currentColors.value),
    value: (mode, drawing) => ({ ...captureMember(drawing, 'pendingPaletteColors'), CDS: '#abcdef' }) },
  // label tables
  labelTextBulkOverrides: { value: () => ({ 'hypothetical protein': 'HP' }) },
  canonicalLabelOverrideRows: { value: () => [
    { recordId: '*', featureType: 'CDS', qualifier: 'product', valueRegex: '^kinase$', labelText: 'Kinase' }
  ] },
  manualPriorityRules: { value: () => [{ feat: 'CDS', order: 'gene,product' }] },
  filterMode: { value: () => 'Whitelist' },
  manualBlacklist: { value: () => 'putative' },
  manualWhitelist: { value: () => [{ feat: 'CDS', qual: 'product', key: 'kinase' }] },
  // annotation sets, as `normalizeAnnotationSets` gives them
  annotationSets: { value: () => normalizeAnnotationSets([{ id: 'sample-set', legendLabel: 'Sample', annotations: [
    { id: 'region_1', target: { kind: 'coordinateSpan', start: 10, end: 20 }, label: 'Sample region' }
  ] }]) },
  // group names (Linear only): Load keeps a name or description only for a
  // group of the drawing's saved Result (`orthogroupEditsOf`, B47); a dormant
  // name, by member set, needs none.
  orthogroupNameOverrides: { modes: ['linear'], otherModeLoad: 'default', needs: ['groups'], value: () => ({ og_1: 'Sample group' }) },
  orthogroupDescriptionOverrides: { modes: ['linear'], otherModeLoad: 'default', needs: ['groups'], value: () => ({ og_1: 'Sample description' }) },
  orthogroupDormantOverrides: { modes: ['linear'], value: () => ({ '["sample-member"]': { name: 'Dormant name', description: 'Dormant description' } }) },
  // canvas
  canvasPadding: { value: () => ({ top: 12, right: 3, bottom: 7, left: 5 }) },
  // feature edits: Load keeps the rows `canonicalFeatureOverrides` and the
  // placement reader of the drawing's mode accept.
  featureOverrides: { value: () => ({ '["record-1","CDS_1"]': {
    recordKey: 'record-1', biologicalFeatureId: 'CDS_1', featureVisibility: 'off', labelVisibility: null, labelText: 'Sample', labelSourceText: null
  } }) },
  featurePlacementOverrides: { value: (mode) => ({ '["record-1","CDS_2"]': {
    recordKey: 'record-1', biologicalFeatureId: 'CDS_2', placement: { kind: 'lane', side: mode === 'linear' ? 'above' : 'outward', level: 1 }
  } }) },
  // paint: Load reads the Legend rows as `normalizeEditorStateData` gives them.
  featureColorOverrides: { value: () => ({ '["record-1","CDS_1"]': { color: '#123456', caption: 'Kinase' } }) },
  featureStrokeOverrides: { value: () => ({ 'record-1\u0000CDS_1': { strokeColor: '#654321', strokeWidth: 2 } }) },
  legendEntries: { value: () => [{ caption: 'Renamed CDS', originalCaption: 'CDS', color: '#abcdef', featureIds: [] }] },
  deletedLegendEntries: { value: () => [{ caption: 'tRNA', originalCaption: 'tRNA', color: '#fedcba', featureIds: [] }] },
  dormantLegendEntries: { value: () => [{ caption: 'Coverage', originalCaption: 'depth', color: '#111111', featureIds: [] }] },
  legendColorOverrides: { value: () => ({ CDS: '#abcdef' }) },
  legendStrokeOverrides: { value: () => ({ CDS: { strokeColor: '#222222', strokeWidth: 1.5 } }) },
  addedLegendCaptions: { value: () => ['Kinase'] },
  // settings, each with one field no default holds
  form: { value: (mode, drawing) => ({ ...captureMember(drawing, 'form'), plot_title: `Sample ${mode} title` }) },
  adv: { value: (mode, drawing) => ({ ...captureMember(drawing, 'adv'), label_font_size: 13 }) },
  losat: { value: (mode, drawing) => {
    const losat = captureMember(drawing, 'losat');
    return { ...losat, blastn: { ...losat.blastn, task: 'blastn' } };
  } },
  losatProgram: { modes: ['linear'], value: () => 'tblastx' },
  circularConservation: { modes: ['circular'], value: (mode, drawing) => ({ ...captureMember(drawing, 'circularConservation'), ring_width: 12 }) },
  linearComparisonPlan: { modes: ['linear'], value: (mode, drawing) => ({ ...captureMember(drawing, 'linearComparisonPlan'), defaultSource: 'upload' }) },
  linearRecordLayoutEnabled: { modes: ['linear'], value: (mode, drawing) => !captureMember(drawing, 'linearRecordLayoutEnabled') },
  linearRecordGap: { modes: ['linear'], value: () => 40 },
  linearRecordRows: { modes: ['linear'], needs: ['records'], value: (mode, drawing) => boundUids(drawing, mode).slice(0, 1).map((uid) => ({ uid, row: 2 })) },
  recordDisplayDrafts: { needs: ['records'], value: (mode, drawing) => boundUids(drawing, mode).slice(0, 1).map((uid) => ({
    sourceUid: uid, selector: '#1', recordId: 'REC_A', topologyOverride: null, startCoordinate: 7,
    reverseComplementOverride: null, anchorIntent: null
  })) },
  // A Linear record row is turned off by its ID, a Circular record by its
  // card's record key (a flat `#N` is the form of a Session 46-50, OV-2869).
  recordsOff: { needs: ['records'], value: (mode, drawing) => boundUids(drawing, mode).slice(0, 1)
    .map((uid) => (mode === 'linear' ? uid : expandedRecordKey(uid, 2))) },
  // A configuration path the Python evaluator accepts.
  unmanagedConfigOverrides: { value: () => ({ 'objects.gc_content.percent_background_opacity': 0.42 }) },
  // The answer to a saved comparison, which only the drawing of a Linear
  // Result holds: Inherit saved comparison.
  importedComparisonIntent: { modes: ['linear'], needs: ['comparison'], byContext: true,
    value: (mode, drawing) => captureMember(drawing, 'importedComparisonIntent') },
  // A drawing keeps the layout preferences of its own mode.
  layoutPreferences: { value: (mode, drawing) => {
    const preferences = captureMember(drawing, 'layoutPreferences');
    return mode === 'linear' ? { ...preferences, linear: { ...preferences.linear, legend: 'top' } }
      : { ...preferences, circular: { ...preferences.circular, single: { ...preferences.circular.single, legend: 'right' } } };
  } },
  linearTypographyLinked: { modes: ['linear'], value: (mode, drawing) => !captureMember(drawing, 'linearTypographyLinked') },
  fileImports: { value: () => ({ blacklist: { words: ['sample'], before: 'putative' } }) },
  // inputs: a file bound in the drawing's own input fields, which the store
  // registers. A Circular drawing binds its GenBank file in a File card, a
  // Linear drawing in its one record row.
  // A color table and the input type are written with a drawing whose
  // records are bound: Save does not write a drawing whose only change is
  // either (reported with the S-1 handoff as a violation of "a file, a
  // changed setting" in R15).
  // The marks of the cards an older Session's flat record positions count
  // over (`c_legacy_record_cards`, OV-2867) are Load's, of a Session 27-50
  // only: a Save writes none.
  files: { needs: ['records'], value: (mode, drawing, { registerFile }) => ({
    ...captureMember(drawing, 'files'), c_legacy_record_cards: null, d_color: registerFile(new File(['CDS\t#123456\n'], `sample-${mode}-colors.tsv`, { lastModified: 1 }))
  }) },
  // Save writes a drawing's record rows (Linear) or File cards (Circular),
  // not both (`serializeActiveRenderFiles` in services/config.js).
  linearSeqs: { modes: ['linear'], otherModeLoad: 'default', value: (mode, drawing, { registerFile }) => [createLinearSeq({
    uid: 'sample-row', gb: registerFile(genbankFile(`sample-${mode}.gb`, [`sample_${mode}`]))
  })] },
  circularSources: { modes: ['circular'], otherModeLoad: 'default', value: (mode, drawing, { registerFile }) => [createCircularSource({
    uid: 'sample-card', gb: registerFile(genbankFile(`sample-${mode}.gb`, [`sample_${mode}`]))
  })] },
  inputType: { needs: ['records'], value: () => 'gff' },
  // What the record read of the drawing's File cards finds.
  circularRecordList: { modes: ['circular'], needs: ['records'], byContext: true, value: (mode, drawing) => captureMember(drawing, 'circularRecordList') },
  circularRecordDiscovery: { modes: ['circular'], needs: ['records'], byContext: true,
    value: (mode, drawing) => captureMember(drawing, 'circularRecordDiscovery') }
};

// A slice's sample: the merge of its members' samples, after their contexts
// (`prepare`). A member keeps the drawing's value when its sample needs a
// fixture's Result, or bound records the drawing does not have, when only an
// action writes it, or when a Session of `mode` does not hold it
// (`otherModeLoad`).
export const sliceSample = (name, mode, drawing, context) => Object.fromEntries(DRAWING_SLICES[name].members.map((member) => {
  const sample = MEMBER_SAMPLES[member];
  const needs = sample.needs ?? [];
  const kept = sample.byContext || needs.some((need) => need !== 'records')
    || (needs.includes('records') && boundUids(drawing, mode).length === 0)
    || (!(sample.modes ?? ['circular', 'linear']).includes(mode) && sample.otherModeLoad === 'default');
  if (!kept) sample.prepare?.(drawing);
  return [member, kept ? captureMember(drawing, member) : sample.value(mode, drawing, context)];
}));

// A capture as JSON without the store's file IDs: a Load binds each input
// file again, under a new ref of the same file. A record row without a file
// is a placeholder (a Circular drawing saves no rows), new at each Load. An
// object's keys are sorted: an object member is a map read by key (a color by
// feature key, a record row's fields), and a Load or History restore may
// build it in another key order; a list keeps its order.
export const refFree = (value) => JSON.stringify(value, (key, entry) => {
  if (entry && typeof entry === 'object' && typeof entry.fileId === 'string' && typeof entry.sourceId === 'string') {
    return { name: entry.name, type: entry.type, size: entry.size, lastModified: entry.lastModified };
  }
  if (key === 'linearSeqs') return entry.map((row) => (row.gb || row.gff || row.fasta ? row : { ...row, uid: '' }));
  if (entry && typeof entry === 'object' && !Array.isArray(entry)) {
    return Object.fromEntries(Object.keys(entry).sort().map((name) => [name, entry[name]]));
  }
  return entry;
});
