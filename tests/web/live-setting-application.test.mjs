// PD-OI-102 (S4-Q1 A LIVE-ON-COMMIT), PD-OI-103 (S4-Q2 A G1-LIVE), R1 (c):
// a committed value of a live setting row (`SETTING_APPLICATION` rows with
// `requestPaths`: the Legend, Plot Title with its Position and Font size, the
// Circular multi-record layout) reaches the Result through the automatic
// rerender. Its request is the committed Session with only those request
// paths replaced, projected by the functions Generate uses (R3, no twin), and
// a successful apply commits only those paths to the committed Session.
import assert from 'node:assert/strict';
import { readFile } from 'node:fs/promises';
import { join } from 'node:path';
import { test } from 'node:test';
import { modeDrawing } from './helpers/drawing-state.mjs';

globalThis.window = {
  Vue: {
    ref: (value) => ({ value }), reactive: (value) => value,
    computed: (getter) => ({ get value() { return getter(); } }),
    nextTick: async () => {}
  },
  DOMPurify: { sanitize: (value) => value }
};

const { state } = await import('../../gbdraw/web/js/state.js');
const config = await import('../../gbdraw/web/js/services/config.js');
const request = await import('../../gbdraw/web/js/services/session-request.js');
const owner = await import('../../gbdraw/web/js/services/setting-application.js');
const { createHistoryManager } = await import('../../gbdraw/web/js/services/history.js');
const { resolveLinearComparisonPlan } = await import('../../gbdraw/web/js/services/linear-comparisons.js');

const record = (id) => `LOCUS       ${id.padEnd(16)}          40 bp    DNA     linear   UNK 01-JAN-1980
DEFINITION  Live setting fixture ${id}.
ACCESSION   ${id}
VERSION     ${id}
FEATURES             Location/Qualifiers
     CDS             1..30
                     /product="test protein"
ORIGIN
        1 atgcatgcat gcatgcatgc atgcatgcat gcatgcatgc
//
`;
const genbank = (ids) => {
  const text = ids.map(record).join('');
  return { name: 'live.gb', type: 'text/plain', size: new TextEncoder().encode(text).byteLength, lastModified: 0, data: btoa(text) };
};

// The request a Generate of `mode` builds from that mode's drawing now.
const generate = (mode) => {
  const drawing = modeDrawing(state, mode);
  config.installShownDrawing(mode);
  state.cInputType.value = 'gb';
  state.lInputType.value = 'gb';
  const filesData = mode === 'circular'
    ? { c_gb: genbank(['LIVE1', 'LIVE2', 'LIVE3']), linearSeqs: [] }
    : {
        linearSeqs: [{
          uid: 'record-1', gb: genbank(['LIVE1']), depth: null, losat_gencode: 1, region_record_id: '',
          region_start: null, region_end: null, region_reverse: false
        }],
        linearComparisons: []
      };
  const comparisonPlanSnapshot = mode === 'linear'
    ? resolveLinearComparisonPlan({
        plan: drawing.linearComparisonPlan, sequences: filesData.linearSeqs, layout: [],
        losatProgram: drawing.losatProgram.value, blastpMode: drawing.losat.blastp.mode
      })
    : null;
  return request.buildCanonicalRenderRequest({ state, drawing, filesData, comparisonPlanSnapshot });
};
const liveRows = () => owner.SETTING_APPLICATION.filter((row) => row.requestPaths);
const allPaths = () => owner.allLiveSettingRequestPaths();
// The automatic rerender's request from `committed` with the live values of `drawing`.
const liveCandidate = (committed, drawing) => request.projectCommittedEditorIntent({
  committed, state, drawing,
  liveSettingWrites: request.projectLiveSettingWrites({ committed, drawing, requestPaths: allPaths() })
});
const LEGEND_KEYS = ['objects.legends.color_rect_size', 'objects.legends.font_size']
  .flatMap((path) => [`${path}.short`, `${path}.long`]);
const LOCK_KEY = 'canvas.linear.keep_definition_left_aligned';
const LIVE_LINEAR_KEYS = [...LEGEND_KEYS, LOCK_KEY];
const pick = (object, keys) => Object.fromEntries(keys.filter((key) => Object.hasOwn(object, key)).map((key) => [key, object[key]]));

test('the G1 rows are the live setting rows, and each request path has its projector', () => {
  assert.deepEqual(liveRows().map((row) => row.id), ['legend', 'plot-title', 'multi-record-layout', 'definition-column']);
  for (const row of liveRows()) assert.equal(row.application, 'live', row.id);
  for (const id of ['colors', 'feature-visibility', 'legend-editor', 'features-editor', 'feature-appearance']) {
    const row = owner.SETTING_APPLICATION.find((candidate) => candidate.id === id);
    assert.equal(row.application, 'live', id);
    assert.equal(row.requestPaths, undefined, `${id} is editor intent, not a draft setting`);
  }
  assert.deepEqual(allPaths().filter((path) => !request.LIVE_SETTING_REQUEST_PATHS.includes(path)), []);
  assert.deepEqual(request.LIVE_SETTING_REQUEST_PATHS.filter((path) => !allPaths().includes(path)), []);
  // Lock Definition Column is live (PD-OI-103 G1), so its help names no Generate
  // (computation-ownership guard); Species, Strain, Default font size and Keep
  // Full Definition stay on Generate (PD-OI-059 rev 2).
  assert.equal(owner.settingApplicationSentence('definition-column'), '');
  assert.equal(owner.settingApplicationNote('layout').application, 'generate');
  assert.equal(owner.settingApplicationNote('titles').application, 'generate');
  assert.equal(owner.settingApplicationNote('species-strain').application, 'generate');
});

test('a History step names the live request paths of the shown drawing it changed', () => {
  // Each drawing has its own History (PD-OI-090): a step's paths are relative to its drawing.
  const step = (...paths) => paths.map((path) => ({ path }));
  assert.deepEqual(owner.liveSettingRequestPaths(step(['config', 'adv', 'legend_font_size'])),
    ['diagramOptions.configOverrides.objects.legends.font_size']);
  // Layout preferences are keyed by mode, and Circular ones by single or
  // multi-record canvas (`createDefaultLayoutPreferences`).
  assert.deepEqual(owner.liveSettingRequestPaths(step(['ui', 'layoutPreferences', 'circular', 'multi', 'plotTitlePosition'])),
    ['diagramOptions.output.plotTitlePosition']);
  assert.deepEqual(owner.liveSettingRequestPaths(step(['ui', 'layoutPreferences', 'circular', 'single', 'legend'])),
    ['diagramOptions.output.legend']);
  assert.deepEqual(owner.liveSettingRequestPaths(step(['ui', 'layoutPreferences', 'linear', 'legend'])),
    ['diagramOptions.output.legend']);
  assert.deepEqual(owner.liveSettingRequestPaths(step(['config', 'adv', 'multi_record_positions', 0, 'row'])),
    ['layout.multiRecordPositions']);
  assert.deepEqual(owner.liveSettingRequestPaths(step(['config', 'form', 'plot_title'])), ['diagramOptions.plotTitle']);
  assert.deepEqual(owner.liveSettingRequestPaths(step(['config', 'form', 'keep_definition_left_aligned'])),
    ['diagramOptions.configOverrides.canvas.linear.keep_definition_left_aligned']);
  // OV-1078: Multi-Record Canvas decides which canvas's positions the draft
  // draws, so its change requests them again.
  assert.deepEqual(owner.liveSettingRequestPaths(step(['config', 'form', 'multi_record_canvas'])),
    ['diagramOptions.output.legend', 'diagramOptions.output.plotTitlePosition']);
  // Settings that apply on Generate name none.
  assert.deepEqual(owner.liveSettingRequestPaths(step(['config', 'form', 'track_type'], ['config', 'adv', 'label_font_size'])), []);
  // A step that replaced the whole drawing config names every live config field.
  assert.ok(owner.liveSettingRequestPaths(step(['config'])).includes('diagramOptions.plotTitle'));
});

test('a Linear live apply draws the G1 values as Generate does and keeps every other committed value', () => {
  const drawing = modeDrawing(state, 'linear');
  const committed = generate('linear');
  drawing.form.legend = 'top';
  Object.assign(drawing.adv, { legend_font_size: 21, legend_box_size: 13, plot_title_font_size: 27, label_font_size: 30 });
  drawing.adv.plot_title_position = 'top';
  drawing.form.plot_title = '  Live title  ';
  drawing.form.align_center = !drawing.form.align_center;
  drawing.form.keep_definition_left_aligned = !drawing.form.keep_definition_left_aligned;
  const draft = generate('linear').renderRequest;
  const candidate = liveCandidate(committed, drawing).renderRequest;
  const options = candidate.diagramOptions;
  assert.deepEqual(options.output, draft.diagramOptions.output);
  assert.equal(options.output.legend, 'top');
  assert.equal(options.plotTitle, 'Live title');
  assert.equal(options.plotTitle, draft.diagramOptions.plotTitle);
  assert.equal(options.plotTitleFontSize, draft.diagramOptions.plotTitleFontSize);
  assert.deepEqual(pick(options.configOverrides, LIVE_LINEAR_KEYS), pick(draft.diagramOptions.configOverrides, LIVE_LINEAR_KEYS));
  assert.notEqual(options.configOverrides[LOCK_KEY], committed.renderRequest.diagramOptions.configOverrides[LOCK_KEY]);
  // The settings that apply on Generate keep their committed values (N-16).
  const notLive = (overrides) => Object.fromEntries(Object.entries(overrides).filter(([key]) => !LIVE_LINEAR_KEYS.includes(key)));
  assert.deepEqual(notLive(options.configOverrides), notLive(committed.renderRequest.diagramOptions.configOverrides));
  assert.notDeepEqual(notLive(draft.diagramOptions.configOverrides), notLive(committed.renderRequest.diagramOptions.configOverrides));
  const { diagramOptions: _candidateOptions, ...candidateRest } = candidate;
  const { diagramOptions: _committedOptions, ...committedRest } = committed.renderRequest;
  assert.deepEqual(candidateRest, committedRest);

  // Blank sends no override, as Generate sends none.
  drawing.adv.legend_font_size = '';
  const blank = liveCandidate({ ...committed, renderRequest: candidate }, drawing).renderRequest;
  assert.deepEqual(pick(blank.diagramOptions.configOverrides, LEGEND_KEYS),
    pick(generate('linear').renderRequest.diagramOptions.configOverrides, LEGEND_KEYS));
  // A value Generate refuses fails the apply.
  drawing.adv.legend_font_size = 'abc';
  assert.throws(() => request.projectLiveSettingWrites({ committed, drawing, requestPaths: allPaths() }),
    (error) => error.code === 'INPUT_INVALID');
  // Values equal to the committed ones write nothing.
  drawing.adv.legend_font_size = 21;
  const applied = { ...committed, renderRequest: candidate };
  assert.deepEqual(request.projectLiveSettingWrites({ committed: applied, drawing, requestPaths: allPaths() }), []);
});

test('a Circular live apply writes the plot title font override and the grid layout as Generate does', () => {
  const drawing = modeDrawing(state, 'circular');
  state.circularRecordList.value = ['LIVE1', 'LIVE2', 'LIVE3'].map((id, index) => ({ selector: `#${index + 1}`, record_id: id }));
  drawing.form.multi_record_canvas = true;
  drawing.recordsOff.splice(0, drawing.recordsOff.length, '#2');
  drawing.adv.multi_record_positions = [{ selector: '#1', row: 1 }, { selector: '#2', row: 1 }, { selector: '#3', row: 1 }];
  const committed = generate('circular');
  assert.equal(committed.renderRequest.grouping, 'grid');
  assert.deepEqual(committed.renderRequest.layout.multiRecordPositions, ['#1@1', '#2@1']);
  drawing.adv.multi_record_positions = [{ selector: '#1', row: 1 }, { selector: '#2', row: 1 }, { selector: '#3', row: 2 }];
  Object.assign(drawing.adv, { multi_record_column_gap_ratio: 0.3, multi_record_size_mode: 'equal', plot_title_font_size: 40 });
  drawing.adv.plot_title_position = 'bottom';
  const draft = generate('circular').renderRequest;
  const candidate = liveCandidate(committed, drawing).renderRequest;
  assert.deepEqual(candidate.layout, draft.layout);
  assert.deepEqual(candidate.layout.multiRecordPositions, ['#1@1', '#2@2']);
  assert.equal(candidate.diagramOptions.output.plotTitlePosition, 'bottom');
  const titleKey = 'objects.definition.circular.plot_title_font_size';
  assert.equal(candidate.diagramOptions.configOverrides[titleKey], 40);
  assert.equal(candidate.diagramOptions.configOverrides[titleKey], draft.diagramOptions.configOverrides[titleKey]);

  // A pending Multi-Record Canvas change applies on Generate: the layout of a
  // committed request without a grid stays as it is.
  drawing.form.multi_record_canvas = false;
  drawing.recordsOff.splice(0);
  const single = generate('circular');
  assert.notEqual(single.renderRequest.grouping, 'grid');
  drawing.form.multi_record_canvas = true;
  assert.deepEqual(liveCandidate(single, drawing).renderRequest.layout, single.renderRequest.layout);
  drawing.form.multi_record_canvas = false;
  drawing.adv.multi_record_positions = [];
});

test('Lock Definition Column OFF on a CLI-shaped request that locks yields the request a Generate builds', () => {
  const drawing = modeDrawing(state, 'linear');
  const initialLock = drawing.form.keep_definition_left_aligned;
  drawing.form.keep_definition_left_aligned = true;
  const generated = generate('linear');
  // A CLI Session: no configOverrides, the resolved config carries the lock.
  const renderRequest = structuredClone(generated.renderRequest);
  delete renderRequest.diagramOptions.configOverrides;
  renderRequest.diagramOptions.config = { canvas: { linear: { keep_definition_left_aligned: true } } };
  const committed = { ...generated, renderRequest };
  const lock = `diagramOptions.configOverrides.${LOCK_KEY}`;
  // Draft still locked: nothing differs, nothing is written.
  assert.deepEqual(request.projectLiveSettingWrites({ committed, drawing, requestPaths: [lock] }), []);
  drawing.form.keep_definition_left_aligned = false;
  assert.equal(request.projectLiveSettingWrites({ committed, drawing, requestPaths: [lock] }).length, 1, 'OFF writes false');
  assert.equal(liveCandidate(committed, drawing).renderRequest.diagramOptions.configOverrides[LOCK_KEY], false);
  assert.equal(liveCandidate(committed, drawing).renderRequest.diagramOptions.configOverrides[LOCK_KEY],
    generate('linear').renderRequest.diagramOptions.configOverrides[LOCK_KEY]);
  // A request that draws OFF stays untouched by an OFF draft.
  delete renderRequest.diagramOptions.config;
  assert.deepEqual(request.projectLiveSettingWrites({ committed, drawing, requestPaths: [lock] }), []);
  drawing.form.keep_definition_left_aligned = initialLock;
});

test('Lock Definition Column applies only to a Linear Result', () => {
  const lock = `diagramOptions.configOverrides.${LOCK_KEY}`;
  const circular = modeDrawing(state, 'circular');
  circular.form.multi_record_canvas = false;
  const committed = generate('circular');
  circular.form.keep_definition_left_aligned = !circular.form.keep_definition_left_aligned;
  assert.equal(request.liveSettingApplies(committed.renderRequest, lock, circular), false);
  assert.deepEqual(request.projectLiveSettingWrites({ committed, drawing: circular, requestPaths: [lock] }), []);
  delete circular.form.keep_definition_left_aligned;
});

// OV-1072: a grid of one ON record (others OFF) keeps its file position in
// its request key; a live Row for record maps `#N` as Generate does.
test('a live Row for record of a grid of one maps the record selector as Generate does', () => {
  const drawing = modeDrawing(state, 'circular');
  state.circularRecordList.value = ['LIVE1', 'LIVE2', 'LIVE3'].map((id, index) => ({ selector: `#${index + 1}`, record_id: id }));
  drawing.form.multi_record_canvas = true;
  drawing.recordsOff.splice(0, drawing.recordsOff.length, '#1', '#2');
  drawing.adv.multi_record_positions = [{ selector: '#3', row: 1 }];
  const committed = generate('circular');
  assert.equal(committed.renderRequest.grouping, 'grid');
  assert.deepEqual(committed.renderRequest.records.map((record) => record.recordKey), ['record-3']);
  assert.deepEqual(committed.renderRequest.layout.multiRecordPositions, ['#1@1']);
  drawing.adv.multi_record_positions = [{ selector: '#3', row: 2 }];
  const draft = generate('circular').renderRequest;
  assert.deepEqual(draft.layout.multiRecordPositions, ['#1@2']);
  assert.deepEqual(liveCandidate(committed, drawing).renderRequest.layout, draft.layout);
  drawing.form.multi_record_canvas = false;
  drawing.recordsOff.splice(0);
  drawing.adv.multi_record_positions = [];
});

// OV-1074, OV-1077: a Circular drawing keeps the Legend and Plot Title
// positions per canvas. Reset Settings turns Multi-Record Canvas on and fills
// the multi-record positions; while that change is pending, the Result keeps
// its positions and any position change applies on Generate.
test('a live position applies only while the draft draws the canvas of the committed request', () => {
  const drawing = modeDrawing(state, 'circular');
  const positions = ['diagramOptions.output.legend', 'diagramOptions.output.plotTitlePosition'];
  drawing.form.multi_record_canvas = false;
  drawing.recordsOff.splice(0);
  const single = generate('circular');
  assert.notEqual(single.renderRequest.grouping, 'grid');
  assert.deepEqual(single.renderRequest.diagramOptions.output, { legend: 'left', plotTitlePosition: 'none' });
  drawing.form.multi_record_canvas = true;
  Object.assign(drawing.layoutPreferences.circular.multi, { legend: 'bottom', plotTitlePosition: 'bottom' });
  assert.equal(drawing.form.legend, 'bottom');
  assert.deepEqual(request.projectLiveSettingWrites({ committed: single, drawing, requestPaths: positions }), []);
  drawing.layoutPreferences.circular.single.legend = 'right';
  assert.deepEqual(request.projectLiveSettingWrites({ committed: single, drawing, requestPaths: positions }), [],
    'a change of the committed canvas\'s position applies on Generate while Multi-Record Canvas is pending on');
  drawing.form.multi_record_canvas = false;
  assert.deepEqual(request.projectLiveSettingWrites({ committed: single, drawing, requestPaths: positions }),
    [[['diagramOptions', 'output', 'legend'], 'right']]);
  drawing.form.multi_record_canvas = true;
  const grid = generate('circular');
  assert.equal(grid.renderRequest.grouping, 'grid');
  drawing.form.multi_record_canvas = false;
  assert.deepEqual(request.projectLiveSettingWrites({ committed: grid, drawing, requestPaths: positions }), [],
    'a grid keeps the multi-record positions while Multi-Record Canvas is pending off');
  Object.assign(drawing.layoutPreferences.circular.single, { legend: 'left', plotTitlePosition: 'none' });
  Object.assign(drawing.layoutPreferences.circular.multi, { legend: null, plotTitlePosition: null });
});

// Review F1: a loaded Session's request may lack optional containers (Gallery
// `tobacco-chloroplast`, a CLI `--save_session` without overrides); a live
// value creates the object Generate would have sent. The editor intent alone
// never needed it (no OV).
test('a live apply onto a committed request without configOverrides creates them, as a Generate sends them', () => {
  const withoutOverrides = (committed) => {
    const renderRequest = structuredClone(committed.renderRequest);
    delete renderRequest.diagramOptions.configOverrides;
    return { ...committed, renderRequest };
  };
  const linear = modeDrawing(state, 'linear');
  Object.assign(linear.adv, { legend_font_size: '', legend_box_size: '' });
  const linearCommitted = withoutOverrides(generate('linear'));
  assert.equal(request.projectCommittedEditorIntent({ committed: linearCommitted, state, drawing: linear })
    .renderRequest.diagramOptions.configOverrides, undefined, 'editor intent alone writes no override');
  Object.assign(linear.adv, { legend_font_size: 19, legend_box_size: 11 });
  const linearOverrides = liveCandidate(linearCommitted, linear).renderRequest.diagramOptions.configOverrides;
  assert.deepEqual(pick(linearOverrides, LEGEND_KEYS),
    pick(generate('linear').renderRequest.diagramOptions.configOverrides, LEGEND_KEYS));
  assert.equal(linearOverrides['objects.legends.font_size.short'], 19);
  // A blank value removes nothing and creates nothing.
  Object.assign(linear.adv, { legend_font_size: '', legend_box_size: '' });
  assert.deepEqual(request.projectLiveSettingWrites({ committed: linearCommitted, drawing: linear, requestPaths: allPaths() })
    .filter(([path]) => path.includes('configOverrides')), []);

  const circular = modeDrawing(state, 'circular');
  circular.form.multi_record_canvas = false;
  circular.adv.plot_title_font_size = '';
  const circularCommitted = withoutOverrides(generate('circular'));
  circular.adv.plot_title_font_size = 38;
  const titleKey = 'objects.definition.circular.plot_title_font_size';
  assert.equal(liveCandidate(circularCommitted, circular).renderRequest.diagramOptions.configOverrides[titleKey], 38);
  circular.adv.plot_title_font_size = '';
});

test('a successful apply commits only its paths, as a new committed Session', () => {
  const drawing = modeDrawing(state, 'linear');
  drawing.adv.legend_font_size = 12;
  const committed = generate('linear');
  config.canonicalRenderArtifactOwner.restore({ committedCanonicalSession: committed, activeSessionResourceTable: null });
  const before = structuredClone(committed);
  const revision = config.readCommittedCanonicalRevision();
  drawing.adv.legend_font_size = 18;
  const writes = request.projectLiveSettingWrites({ committed, drawing, requestPaths: allPaths() });
  assert.deepEqual(writes.map(([path]) => path.join('.')),
    ['diagramOptions.configOverrides.objects.legends.font_size.short', 'diagramOptions.configOverrides.objects.legends.font_size.long']);
  config.commitLiveSettingValues(writes);
  const next = config.getCommittedCanonicalSession();
  assert.notEqual(next, committed, 'a History artifact holding the previous Session keeps it');
  assert.deepEqual(committed, before);
  assert.equal(next.resources, committed.resources);
  assert.equal(config.readCommittedCanonicalRevision(), revision + 1);
  const expected = structuredClone(before.renderRequest);
  expected.diagramOptions.configOverrides['objects.legends.font_size.short'] = 18;
  expected.diagramOptions.configOverrides['objects.legends.font_size.long'] = 18;
  assert.deepEqual(next.renderRequest, expected);
  config.commitLiveSettingValues([]);
  assert.equal(config.readCommittedCanonicalRevision(), revision + 1, 'no write, no replacement');
  config.canonicalRenderArtifactOwner.restore(null);
});

// Review F3: History is the one notifier: a recorded step, a committed
// checkpoint (Reset Settings, the Bakta preset), and the Undo or Redo of each.
test('History reports each step, checkpoint, Undo and Redo with its changes, and a checkpoint restore as any setting', async () => {
  let value = { config: { adv: { legend_font_size: 10 } } };
  const steps = [];
  const history = createHistoryManager({
    buildIntent: () => structuredClone(value),
    applyIntent: (intent) => { value = structuredClone(intent); },
    buildCheckpoint: () => structuredClone(value),
    applyCheckpoint: (checkpoint) => { value = checkpoint; }
  });
  history.setSettingsStepListener((changes) => steps.push(changes && changes.map((change) => change.path.join('.'))));
  await history.captureBaseline('initial');
  await history.runUndoable('Legend font size', () => { value.config.adv.legend_font_size = 14; });
  assert.deepEqual(steps, [['config.adv.legend_font_size']]);
  await history.undo();
  await history.redo();
  assert.deepEqual(steps, [['config.adv.legend_font_size'], ['config.adv.legend_font_size'], ['config.adv.legend_font_size']],
    'an Undo and a Redo name the paths of the step they restore');
  steps.length = 0;
  await history.runUndoableCheckpoint('Reset', () => { value.config.adv.legend_font_size = 8; });
  assert.deepEqual(steps, [['config.adv.legend_font_size']], 'a committed checkpoint names the settings it changed');
  steps.length = 0;
  await history.undo();
  assert.deepEqual(steps, [null], 'a checkpoint restore names no paths: every live row is asked');
});

test('a live setting row the shown Result cannot draw says Applies on Generate', () => {
  assert.equal(owner.settingApplicationNote('multi-record-layout').text, 'Live edit');
  const generateNote = { application: 'generate', text: 'Applies on Generate', qualified: false };
  // A single-canvas Result: the layout fields apply on Generate (PD-OI-103).
  const single = { mode: 'circular', grouping: 'single' };
  const draft = (multiRecordCanvas) => ({ form: { multi_record_canvas: multiRecordCanvas } });
  const note = (id, committed, multiRecordCanvas) => owner.settingApplicationNote(
    id, (path) => request.liveSettingApplies(committed, path, draft(multiRecordCanvas))
  );
  assert.equal(request.liveSettingApplies(single, 'layout.multiRecordColumnGapRatio', draft(true)), false);
  assert.equal(request.liveSettingApplies({ mode: 'circular', grouping: 'grid' }, 'layout.multiRecordColumnGapRatio', draft(true)), true);
  assert.equal(request.liveSettingApplies(single, 'diagramOptions.output.legend', draft(false)), true);
  assert.deepEqual(note('multi-record-layout', single, true), generateNote);
  assert.equal(note('legend', single, false).text, 'Live edit');
  // OV-1077: during a Multi-Record Canvas change the positions apply on
  // Generate; the sizes, font sizes and the title text stay live.
  const partial = { application: 'live', text: 'Live edit; the position applies on Generate after the Multi-Record Canvas change', qualified: true };
  assert.equal(request.liveSettingApplies(single, 'diagramOptions.output.legend', draft(true)), false);
  assert.equal(request.liveSettingApplies(single, 'diagramOptions.plotTitleFontSize', draft(true)), true);
  assert.deepEqual(note('legend', single, true), partial);
  assert.deepEqual(note('plot-title', single, true), partial);
  assert.deepEqual(note('plot-title', { mode: 'circular', grouping: 'grid' }, false), partial);
  assert.equal(note('legend', { mode: 'linear', grouping: 'single' }, false).text, 'Live edit');
  // A row without request paths keeps its own note.
  assert.equal(owner.settingApplicationNote('track-visibility', () => false).qualified, true);
});

test('no live note text is written outside the owner', async () => {
  const html = await readFile(join(process.cwd(), 'gbdraw', 'web', 'index.html'), 'utf8');
  assert.deepEqual([...html.matchAll(/Live edit:|Updates the current Result/g)].map((match) => match[0]), []);
  for (const id of ['legend-editor', 'features-editor', 'feature-appearance']) {
    assert.match(html, new RegExp(`settingApplicationNote\\('${id}'\\)\\.text`), id);
  }
});
