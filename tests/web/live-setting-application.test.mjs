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
  const drawing = state.drawings[mode];
  state.mode.value = mode;
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
const pick = (object, keys) => Object.fromEntries(keys.filter((key) => Object.hasOwn(object, key)).map((key) => [key, object[key]]));

test('the G1 rows are the live setting rows, and each request path has its projector', () => {
  assert.deepEqual(liveRows().map((row) => row.id), ['legend', 'plot-title', 'multi-record-layout']);
  for (const row of liveRows()) assert.equal(row.application, 'live', row.id);
  for (const id of ['colors', 'feature-visibility', 'legend-editor', 'features-editor', 'feature-appearance']) {
    const row = owner.SETTING_APPLICATION.find((candidate) => candidate.id === id);
    assert.equal(row.application, 'live', id);
    assert.equal(row.requestPaths, undefined, `${id} is editor intent, not a draft setting`);
  }
  assert.deepEqual(allPaths().filter((path) => !request.LIVE_SETTING_REQUEST_PATHS.includes(path)), []);
  assert.deepEqual(request.LIVE_SETTING_REQUEST_PATHS.filter((path) => !allPaths().includes(path)), []);
  // Lock Definition Column waits for its guard retarget (computation-ownership); Species,
  // Strain, Default font size and Keep Full Definition stay on Generate (PD-OI-059 rev 2).
  assert.equal(owner.settingApplicationNote('layout').application, 'generate');
  assert.equal(owner.settingApplicationNote('titles').application, 'generate');
  assert.equal(owner.settingApplicationNote('species-strain').application, 'generate');
});

test('a History step names the live request paths of the shown drawing it changed', () => {
  const step = (...paths) => paths.map((path) => ({ path }));
  assert.deepEqual(owner.liveSettingRequestPaths(step(['modes', 'linear', 'config', 'adv', 'legend_font_size']), 'linear'),
    ['diagramOptions.configOverrides.objects.legends.font_size']);
  // Layout preferences are keyed by mode, and Circular ones by single or
  // multi-record canvas (`createDefaultLayoutPreferences`).
  assert.deepEqual(owner.liveSettingRequestPaths(step(['modes', 'circular', 'ui', 'layoutPreferences', 'circular', 'multi', 'plotTitlePosition']), 'circular'),
    ['diagramOptions.output.plotTitlePosition']);
  assert.deepEqual(owner.liveSettingRequestPaths(step(['modes', 'circular', 'ui', 'layoutPreferences', 'circular', 'single', 'legend']), 'circular'),
    ['diagramOptions.output.legend']);
  assert.deepEqual(owner.liveSettingRequestPaths(step(['modes', 'linear', 'ui', 'layoutPreferences', 'linear', 'legend']), 'linear'),
    ['diagramOptions.output.legend']);
  // An Undo passes the shown drawing's paths without the mode prefix.
  assert.deepEqual(owner.liveSettingRequestPaths(step(['config', 'adv', 'multi_record_positions', 0, 'row']), 'circular'),
    ['layout.multiRecordPositions']);
  assert.deepEqual(owner.liveSettingRequestPaths(step(['config', 'form', 'plot_title']), 'linear'), ['diagramOptions.plotTitle']);
  // The other drawing, and settings that apply on Generate, name none.
  assert.deepEqual(owner.liveSettingRequestPaths(step(['modes', 'circular', 'config', 'adv', 'legend_font_size']), 'linear'), []);
  assert.deepEqual(owner.liveSettingRequestPaths(step(
    ['modes', 'linear', 'config', 'form', 'species'], ['modes', 'linear', 'config', 'adv', 'def_font_size'],
    ['modes', 'linear', 'config', 'form', 'keep_definition_left_aligned'], ['modes', 'linear', 'config', 'form', 'multi_record_canvas']
  ), 'linear'), []);
  // A step that replaced the whole drawing config names every live config field.
  assert.ok(owner.liveSettingRequestPaths(step(['modes', 'linear', 'config']), 'linear').includes('diagramOptions.plotTitle'));
});

test('a Linear live apply draws the G1 values as Generate does and keeps every other committed value', () => {
  const drawing = state.drawings.linear;
  const committed = generate('linear');
  drawing.form.legend = 'top';
  Object.assign(drawing.adv, { legend_font_size: 21, legend_box_size: 13, plot_title_font_size: 27, def_font_size: 30 });
  drawing.adv.plot_title_position = 'top';
  drawing.form.plot_title = '  Live title  ';
  drawing.form.align_center = !drawing.form.align_center;
  const draft = generate('linear').renderRequest;
  const candidate = liveCandidate(committed, drawing).renderRequest;
  const options = candidate.diagramOptions;
  assert.deepEqual(options.output, draft.diagramOptions.output);
  assert.equal(options.output.legend, 'top');
  assert.equal(options.plotTitle, 'Live title');
  assert.equal(options.plotTitle, draft.diagramOptions.plotTitle);
  assert.equal(options.plotTitleFontSize, draft.diagramOptions.plotTitleFontSize);
  assert.deepEqual(pick(options.configOverrides, LEGEND_KEYS), pick(draft.diagramOptions.configOverrides, LEGEND_KEYS));
  // The settings that apply on Generate keep their committed values (N-16).
  const notLive = (overrides) => Object.fromEntries(Object.entries(overrides).filter(([key]) => !LEGEND_KEYS.includes(key)));
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
  const drawing = state.drawings.circular;
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

test('a successful apply commits only its paths, as a new committed Session', () => {
  const drawing = state.drawings.linear;
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

test('History reports each recorded intent step with its changes, and a checkpoint restore as any setting', async () => {
  let value = { config: { adv: { legend_font_size: 10 } } };
  const steps = [];
  const history = createHistoryManager({
    buildIntent: () => structuredClone(value),
    applyIntent: (intent) => { value = intent; },
    buildCheckpoint: () => structuredClone(value),
    applyCheckpoint: (checkpoint) => { value = checkpoint; }
  });
  history.setSettingsStepListener((changes) => steps.push(changes && changes.map((change) => change.path.join('.'))));
  await history.captureBaseline('initial');
  await history.runUndoable('Legend font size', () => { value.config.adv.legend_font_size = 14; });
  assert.deepEqual(steps, [['config.adv.legend_font_size']]);
  await history.undo();
  assert.equal(steps.length, 1, 'an Undo records no step; its restore asks for the apply');
  await history.runUndoableCheckpoint('Reset', () => { value.config.adv.legend_font_size = 8; });
  steps.length = 0;
  await history.undo();
  assert.deepEqual(steps, [null], 'a checkpoint restore names no paths: every live row is asked');
});

test('no live note text is written outside the owner', async () => {
  const html = await readFile(join(process.cwd(), 'gbdraw', 'web', 'index.html'), 'utf8');
  assert.deepEqual([...html.matchAll(/Live edit:|Updates the current Result/g)].map((match) => match[0]), []);
  for (const id of ['legend-editor', 'features-editor', 'feature-appearance']) {
    assert.match(html, new RegExp(`settingApplicationNote\\('${id}'\\)\\.text`), id);
  }
});
