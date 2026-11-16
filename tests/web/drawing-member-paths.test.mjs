// Every member of every drawing slice survives the six paths that carry a
// drawing's state, on the real app (`createAppSetup`), in a Circular and a
// Linear drawing: Save then Load, a History step and its Undo and Redo, a
// drawing switch and back, a switch to the other mode's drawing and back,
// Reset this drawing then Undo, and Duplicate. For each member the guard
// installs only that member's sample (helpers/drawing-member-samples.mjs) as
// one History step, runs the paths, and compares the member after each path
// with the sample (file bindings compared by their file, not their store ID).
// The guard reports each loss with its member, mode, and path; a test proves
// it reports exactly the member a deliberately broken path loses. The sample
// table holds one sample per member of the drawing key classes, and the
// samples together change each of the Session's drawing config domains,
// which Save then Load then Save gives back.
import assert from 'node:assert/strict';
import { File } from 'node:buffer';
import { readFileSync } from 'node:fs';
import test from 'node:test';
import { gunzipSync } from 'node:zlib';
import { installFakeSvgDom } from './fake-svg-dom.mjs';
import { installSessionImportWorker } from './helpers/session-import-node.mjs';
import { withBootPalettes } from './helpers/boot-palettes.mjs';
import { validateWebConfigOverrides } from './helpers/python-evaluator.mjs';

globalThis.window = {
  Vue: {
    ref: (value) => ({ value }), shallowRef: (value) => ({ value }), reactive: (value) => value,
    isRef: (value) => Boolean(value && typeof value === 'object' && 'value' in value),
    computed: (getter) => ({ get value() { return getter(); } }), nextTick: async () => {},
    watch: () => () => {}, onMounted: () => {}, onUnmounted: () => {}, provide: () => {}
  },
  DOMPurify: { sanitize: (value) => value }
};
globalThis.document = { createElement: () => ({ click: () => {} }), querySelector: () => null, documentElement: { setAttribute: () => {} } };
installFakeSvgDom();
globalThis.File = File;
globalThis.alert = () => {};
// A Load prefetches the Worker's boot files unless the browser saves data; a node test fetches nothing.
Object.defineProperty(globalThis.navigator, 'connection', { value: { saveData: true }, configurable: true });
installSessionImportWorker();

const config = await import('../../gbdraw/web/js/services/config.js');
const { state, createDrawingState } = await import('../../gbdraw/web/js/state.js');
const { withNewDrawingDefaults } = await import('../../gbdraw/web/js/services/reset.js');
const { createAppSetup, transformLegacyResultSvg } = await import('../../gbdraw/web/js/app/app-setup.js');
const { readSessionContainer } = await import('../../gbdraw/web/js/services/session-container.js');
const { editedSessionFile } = await import('./helpers/written-resources.mjs');
const { MODE_SLICE_CONTAINERS } = await import('../../gbdraw/web/js/services/mode-scoped-migration.js');
const { CURRENT_WRITER_ACTIVE_CONFIG_DOMAINS } = await import('../../gbdraw/web/js/services/session-active-config-contract.js');
const { DRAWING_SLICES } = await import('./helpers/drawing-state.mjs');
const { callerSession } = await import('./helpers/session-caller.mjs');
const { MEMBER_SAMPLES, LEAF_VARIED_MEMBERS, DRAWING_NEEDS, unvariedOf, captureMember, sliceOfMember, genbankFile, refFree, variedLeaves,
  leafPaths, leafAt } = await import('./helpers/drawing-member-samples.mjs');
const { chooseRecordLayout } = await import('../../gbdraw/web/js/services/one-circle-layout.js');

const MODES = /** @type {const} */ (['circular', 'linear']);
const MEMBERS = Object.values(DRAWING_SLICES).flatMap((slice) => slice.members);
const app = createAppSetup();
const history = window.__GBDRAW_HISTORY__;
// The app checks configuration overrides in its diagram Worker; a node test checks them with the Python evaluator.
config.setUnmanagedConfigOverrideValidator(async (payload) => ({ result: await validateWebConfigOverrides(payload) }));
const context = { registerFile: (file) => state.inputFiles.registerFile(file) };
// The app's asynchronous work after an action (a record read, the end of a Load) runs in a few turns.
const settle = async () => { for (let i = 0; i < 20; i += 1) await new Promise((resolve) => setImmediate(resolve)); };
const installMember = (drawing, member, value) => {
  const name = sliceOfMember(member);
  DRAWING_SLICES[name].install(drawing.slices[name], { ...DRAWING_SLICES[name].capture(drawing.slices[name]), [member]: value });
};

// The app's Save (Save Session), under a new title each time (Save asks
// before it downloads a file name again), and the app's Load.
let saves = 0;
const saveSession = async () => {
  saves += 1;
  state.sessionTitle.value = `member-paths-${saves}`;
  const saved = await app.saveSessionWithTitle();
  return saved?.status === 'saved' ? { blob: saved.blob, refused: null }
    : { blob: null, refused: `${saved?.status} ${saved?.error?.code ?? ''} ${JSON.stringify(saved?.error?.context ?? null)}` };
};
// A Load that asks before it loads (a part of the file it cannot restore)
// fails the case with the dialog's text; the answer is Cancel.
const loadSession = async (blob) => {
  let done = false;
  const loading = app.importSession({ target: { files: [new File([blob], 'member-paths.gbdraw')], value: '' } })
    .finally(() => { done = true; });
  while (!done) {
    await new Promise((resolve) => setImmediate(resolve));
    if (app.choicePrompt.show) {
      const asked = `${app.choicePrompt.title} ${app.choicePrompt.description}`;
      app.resolveChoicePrompt('cancel');
      await loading;
      return { asked };
    }
  }
  const loaded = await loading;
  assert.equal(loaded?.status, 'ok', loaded?.error?.message ?? loaded?.status);
  // The end of a Load reads each drawing's records again.
  await settle();
  return { asked: null };
};
/** @typedef {import('./helpers/written-resources.mjs').SessionDocument} SessionDocument A saved Session document. */
/** @param {Blob} blob @returns {Promise<SessionDocument>} */
const savedDocument = async (blob) => /** @type {SessionDocument} */ ((await readSessionContainer(blob)).document);
/** @param {SessionDocument} document @param {string} id */
const savedConfigOf = (document, id) => {
  const entry = document.drawings.find((drawing) => drawing.id === id);
  return entry ? entry.config : null;
};
// The Reset dialog's answer, as the user clicks it.
const answerReset = async (answer) => {
  const asked = app.requestResetSettings();
  await new Promise((resolve) => setTimeout(resolve, 0));
  app.resolveChoicePrompt(answer);
  return asked;
};

// A base Session a case starts from: a drawing of each mode that Save writes
// (its plot title set), the peers of the case, the Circular one shown at Save.
// The plain base holds only them; a fixture base also holds the drawings of a
// fixture Session, among them the drawing a sample that needs its Result
// works on (`target`). A fixture is loaded once, through the Load service
// with the app's rollback port: the app's Load waits for the preview of a
// shown Result, which a node test does not mount; the base Session it saves
// shows a peer, so the app's Load of it shows no Result.
const FIXTURES = new URL('../fixtures/sessions/', import.meta.url);
const CONTEXTS = {
  // A Linear drawing whose saved Result holds a comparison read from a CLI table.
  comparison: { fixture: 'BGC0000708-BGC0000709-comparisons-table.v50-cli0145.gbdraw-session.json.gz', target: 'linear' },
  // A Linear drawing whose Result has Similarity groups (og_1, ...).
  groups: { fixture: 'BGC0000708-BGC0000709-copy-crop.v51-web.gbdraw', target: 'linear-2' }
};
const bases = {};
const buildBase = async (name) => {
  const fixture = CONTEXTS[name] ?? null;
  const before = state.drawingTabs.value.map((tab) => tab.id);
  if (fixture) {
    const raw = readFileSync(new URL(fixture.fixture, FIXTURES));
    const { importSession } = callerSession(config, state, app.sessionImportRollbackState);
    const loaded = await importSession({ target: { files: [new File([fixture.fixture.endsWith('.gz') ? gunzipSync(raw) : raw], fixture.fixture)], value: '' } },
      { transformLegacyResultSvg, confirmLoadRepairs: async () => true });
    assert.equal(loaded?.status, 'ok', `${fixture.fixture}: ${state.errorLog.value?.summary ?? loaded?.status}`);
  }
  const peers = {};
  for (const mode of MODES) {
    peers[mode] = app.drawingList.newDrawing(mode).drawingId;
    state.drawingOf(peers[mode]).form.plot_title = 'Peer';
  }
  // The plain base starts from the boot drawing, made before the palettes: the peers replace it.
  if (!fixture) before.forEach((id) => app.drawingList.deleteDrawing(id));
  await settle();
  app.selectDrawing(peers.circular);
  const { blob, refused } = await saveSession();
  assert.equal(refused, null);
  bases[name] = { blob, peers, target: fixture?.target ?? null };
};

// The context a sample needs, built in the shown drawing as the user builds
// it: a GenBank file of two records added as a File card (Circular) or chosen
// in the record row (Linear), and its records read. One History step.
const bindRecords = async (drawing) => {
  const file = state.inputFiles.registerFile(genbankFile(`records-${drawing.mode}.gb`, ['REC_A', 'REC_B']));
  await history.runUndoable('Change uploaded file', () => (drawing.mode === 'linear'
    ? app.setLinearSeqPrimaryFile(0, 'gb', file) : app.addCircularSourceFiles([file])));
  await settle();
  if (drawing.mode === 'circular') assert.equal(drawing.circularRecordDiscovery.status, 'ready');
};

// The ring table a 'rings' sample changes the ring row of: a BLAST table the
// user uploads for the rings (Upload BLAST TSV). One History step.
const bindRingTable = async (drawing) => {
  const table = state.inputFiles.registerFile(new File(['q\ts\t90\t10\t0\t0\t1\t10\t1\t10\t1e-5\t50\n'], 'ring.tsv', { lastModified: 1 }));
  await history.runUndoable('Change uploaded file', () => {
    app.changeCircularConservationSource('upload');
    return app.setCircularConservationUploadFiles([table]);
  });
  await settle();
  assert.equal(drawing.circularConservation.series.length, 1, 'the table has its ring row');
};
// The two record rows a 'rows' sample pairs: a Linear drawing's first two
// rows, each holding a GenBank file, as the user adds them. One History step.
const bindRows = async (drawing) => {
  await history.runUndoable('Change uploaded file', async () => {
    for (let index = 0; index < 2; index += 1) {
      if (drawing.linearSeqs[index]?.gb) continue;
      if (!drawing.linearSeqs[index]) await app.addLinearSeq();
      await app.setLinearSeqPrimaryFile(index, 'gb', state.inputFiles.registerFile(genbankFile(`rows-${index}.gb`, [`ROW_${index}`])));
    }
  });
  await settle();
  assert.equal(drawing.linearSeqs.filter((row) => row.gb).length >= 2, true, 'the drawing has two record rows');
};
// The contexts a test builds in the drawing (`DRAWING_NEEDS`), by the modes that have them.
const DRAWING_CONTEXTS = {
  records: { circular: bindRecords, linear: bindRecords }, rings: { circular: bindRingTable }, rows: { linear: bindRows }
};
const buildDrawingContexts = async (drawing, needs) => {
  for (const need of needs.filter((name) => DRAWING_NEEDS.includes(name))) await DRAWING_CONTEXTS[need][drawing.mode]?.(drawing);
};

// The actions a `byContext` member's edit is, each one History step.
const CONTEXT_EDITS = {
  records: bindRecords,
  comparison: () => app.inheritImportedComparison()
};

/**
 * One case's drawings: the drawing it works on (`id`, shown), a second drawing
 * of its mode (`second`), the other mode's peer (`other`), and its own mode's
 * peer (`peer`), which is shown at Save.
 * @typedef {{ id: string, second: string, other: string, peer: string, mode: 'circular' | 'linear' }} CaseDrawings
 */
// A member that Save counts as a change only in the drawings of its `modes`
// runs in a drawing of the other mode too (`foreign`): that drawing is
// written for its plot title, and its Load gives the sample back, or the
// mode's default when the sample says so (`otherModeLoad`).
/** @returns {Promise<CaseDrawings>} */
const startCase = async (mode, member, foreign) => {
  const needs = MEMBER_SAMPLES[member].needs ?? [];
  const base = bases[foreign ? 'plain' : needs.find((need) => CONTEXTS[need]) ?? 'plain'];
  assert.equal((await loadSession(base.blob)).asked, null);
  const second = app.drawingList.newDrawing(mode).drawingId;
  const id = base.target ?? app.drawingList.newDrawing(mode).drawingId;
  if (state.activeDrawingId.value !== id) assert.deepEqual(app.selectDrawing(id), { status: 'ok' });
  const drawing = state.drawingOf(id);
  if (!MEMBER_SAMPLES[member].byContext) await buildDrawingContexts(drawing, needs);
  MEMBER_SAMPLES[member].prepare?.(drawing);
  if (foreign) drawing.form.plot_title = 'Used';
  return { id, second, other: base.peers[mode === 'circular' ? 'linear' : 'circular'], peer: base.peers[mode], mode };
};

/**
 * A path's checks: what the member holds after it, and what it should hold.
 * @typedef {{ where: string, expected: unknown, actual: unknown }} PathCheck
 * @typedef {{ base: unknown, sample: unknown, loaded: unknown }} CaseValues What the member held before the
 *   edit, after it, and what a Load of the Session gives back (the sample, or the default per `otherModeLoad`).
 * @typedef {(drawings: CaseDrawings, member: string, values: CaseValues) => Promise<PathCheck[]>} StatePath
 */
const isObject = (value) => Boolean(value) && typeof value === 'object' && !Array.isArray(value);
// What Reset this drawing leaves of a member, from the rules that say so, not
// from Reset's own code (`resetSettings`, services/reset.js, which the dialog
// runs): a drawing's reset resets every slice through the slice's own `reset`
// (R15 in gbdraw/web/CLAUDE.md), which keeps the input bindings and record
// choices PD-OI-070 keeps (services/drawing/inputs.js, settings.js); here it
// runs on a new drawing holding the sample. `RESET_DECLARED` gives the
// members another rule decides.
/** @type {Record<string, (values: { fresh: unknown, sample: unknown }) => unknown>} */
const RESET_DECLARED = {
  // Reset returns a drawing's settings and edits to a new drawing's
  // (PD-OI-090, Reset): the palette is a new drawing's, with its colors, which
  // the palette slice's own reset leaves empty.
  ...Object.fromEntries(DRAWING_SLICES.palette.members.map((member) => [member, ({ fresh }) => fresh])),
  // The Legend rows are the rows Python drew, which Reset rebuilds (the paint
  // intent in gbdraw/web/CLAUDE.md): they stay listed, a renamed row back to
  // the caption Python draws.
  legendEntries: ({ sample }) => /** @type {{ originalCaption: string }[]} */ (sample).map((row) => ({ ...row, caption: row.originalCaption }))
};
// The fields that Reset's follow-ups for the drawing's mode
// (`createDrawingReset`, app/reset-settings.js) derive again from what Reset
// keeps, by mode and member (`true`: the whole member); Undo checks them.
/** @type {Record<'circular' | 'linear', Record<string, true | readonly string[]>>} */
const RESET_DERIVED = {
  // The comparison rings are derived again from the kept ring files, and a
  // ring search of the old settings goes with them (OV-1167, OV-2758).
  circular: { files: ['c_conservation_blasts', 'c_conservation_blasts_source'] },
  // The record rows follow the kept Files again (`applyLinearSeqMutation`).
  linear: { linearRecordRows: true }
};
const resetExpectation = (mode, member, sample) => {
  const fresh = captureMember(withNewDrawingDefaults(state, createDrawingState(`${mode}-fresh`, mode)), member);
  const drawing = withNewDrawingDefaults(state, createDrawingState(`${mode}-reset`, mode));
  installMember(drawing, member, sample);
  const name = sliceOfMember(member);
  DRAWING_SLICES[name].reset(drawing.slices[name]);
  return RESET_DECLARED[member]?.({ fresh, sample }) ?? captureMember(drawing, member);
};
// `value` without the fields Reset's follow-ups derive again.
const underivedPart = (mode, member, value) => {
  const derived = RESET_DERIVED[mode][member];
  if (derived === true) return 'derived again';
  if (!derived || !isObject(value)) return value;
  return Object.fromEntries(Object.entries(value).filter(([key]) => !derived.includes(key)));
};
const same = (member, actual, expected) => refFree({ [member]: actual }) === refFree({ [member]: expected });
const held = (id, member) => (state.drawingTabs.value.some((tab) => tab.id === id)
  ? captureMember(state.drawingOf(id), member) : `no drawing ${id}`);

/** @type {Record<string, StatePath>} */
const PATHS = {
  // The edit was one History step: Undo gives back the value before it, Redo the sample.
  history: async ({ id }, member, { base, sample }) => {
    await app.undoHistory();
    const undone = held(id, member);
    await app.redoHistory();
    return [{ where: 'Undo', expected: base, actual: undone }, { where: 'Redo', expected: sample, actual: held(id, member) }];
  },
  drawingSwitch: async ({ id, second }, member, { sample }) => {
    assert.deepEqual(app.selectDrawing(second), { status: 'ok' });
    assert.deepEqual(app.selectDrawing(id), { status: 'ok' });
    return [{ where: 'switch and back', expected: sample, actual: held(id, member) }];
  },
  modeRoundTrip: async ({ id, other }, member, { sample }) => {
    assert.deepEqual(app.selectDrawing(other), { status: 'ok' });
    assert.deepEqual(app.selectDrawing(id), { status: 'ok' });
    return [{ where: 'other mode and back', expected: sample, actual: held(id, member) }];
  },
  // Reset this drawing (the dialog's answer) leaves each member as
  // `resetExpectation` says, but the fields its mode's follow-ups derive
  // again; Undo gives the sample back. The Reset's step is the newest of the
  // drawing's History (by its sequence); a Reset that changes nothing records none.
  resetUndo: async ({ id, mode }, member, { sample }) => {
    const top = () => history.timeline().undo[0]?.sequence ?? null;
    const before = top();
    await answerReset('reset-drawing');
    const reset = { where: 'Reset this drawing', expected: underivedPart(mode, member, resetExpectation(mode, member, sample)),
      actual: underivedPart(mode, member, held(id, member)) };
    if (top() !== before) await app.undoHistory();
    return [reset, { where: 'Reset this drawing, Undo', expected: sample, actual: held(id, member) }];
  },
  duplicate: async ({ id }, member, { sample }) => {
    const { drawingId: copy } = app.drawingList.duplicateDrawing(id);
    const checks = [{ where: 'the copy', expected: sample, actual: held(copy, member) },
      { where: 'the source', expected: sample, actual: held(id, member) }];
    await app.undoHistory();
    if (state.activeDrawingId.value !== id) app.selectDrawing(id);
    return checks;
  },
  saveLoad: saveThenLoad()
};

/**
 * Save then Load. Save runs while another drawing is shown, so only the
 * member makes Save write this drawing (with the context it needs). A Load
 * makes every drawing a new state, so the loaded value comes from the file.
 * `rewrite` changes the saved document before the Load (a broken Save in the
 * guard's own test).
 * @param {(document: SessionDocument, id: string) => void} [rewrite]
 * @returns {StatePath}
 */
function saveThenLoad(rewrite) {
  return async ({ id, peer }, member, { loaded }) => {
    assert.deepEqual(app.selectDrawing(peer), { status: 'ok' });
    const { blob: saved, refused } = await saveSession();
    if (refused) return [{ where: 'Save', expected: loaded, actual: `Save refused: ${refused}` }];
    const blob = rewrite ? await editedSessionFile(saved, (document) => rewrite(document, id)) : saved;
    const { asked } = await loadSession(blob);
    if (asked) return [{ where: 'Save then Load', expected: loaded, actual: `Load asked: ${asked}` }];
    // A Load of a current Session reads no records (OV-2881): the list shows
    // after a read, which the Record Order's refresh starts. The drawings of
    // these cases hold no request, so their discovery is idle after Load and
    // their list empty until that read. Open question: record list after
    // Load without a request (lane S, OV-3696 analysis); no assertion either
    // way until it is answered. A drawing with a request is deferred after
    // Load (tests/web/drawing-slices.test.mjs, test 6).
    if (MEMBER_SAMPLES[member].byContext && MEMBER_SAMPLES[member].needs.includes('records')) {
      assert.deepEqual(app.selectDrawing(id), { status: 'ok' });
      await app.refreshCircularRecordOrder();
      await settle();
    }
    return [{ where: 'Save then Load', expected: loaded, actual: held(id, member) }];
  };
}

/**
 * The losses of `paths` (run in order) for each member in a drawing of
 * `mode`: a member that does not hold, after a path, what the path should
 * leave. A path that loses the sample is followed by the sample again, so a
 * later path starts from it.
 * @param {Record<string, StatePath>} paths
 * @param {'circular' | 'linear'} mode
 * @param {readonly string[]} [members]
 */
const memberPathLosses = async (paths, mode, members = MEMBERS) => {
  const losses = [];
  // A member that only an action of its mode writes (`byContext`) runs in that mode.
  for (const member of members.filter((name) => !MEMBER_SAMPLES[name].byContext || (MEMBER_SAMPLES[name].modes ?? MODES).includes(mode))) {
    const foreign = !(MEMBER_SAMPLES[member].modes ?? MODES).includes(mode);
    const drawings = await startCase(mode, member, foreign);
    const drawing = state.drawingOf(drawings.id);
    const base = captureMember(drawing, member);
    // A member that only a context action writes (the records a file read) takes that action as its step.
    let sample;
    if (MEMBER_SAMPLES[member].byContext) {
      await CONTEXT_EDITS[MEMBER_SAMPLES[member].needs.find((need) => CONTEXT_EDITS[need])](drawing);
      sample = MEMBER_SAMPLES[member].value(mode, drawing, context);
    } else {
      sample = MEMBER_SAMPLES[member].value(mode, drawing, context);
      await history.runUndoable(`Sample ${member}`, () => installMember(drawing, member, sample));
    }
    // A sample is an edit: it differs from what the drawing held and from a new drawing's value.
    const fresh = captureMember(withNewDrawingDefaults(state, createDrawingState(`${mode}-fresh`, mode)), member);
    if (same(member, sample, base) || same(member, sample, fresh)) {
      losses.push({ member, mode, path: 'sample', where: 'the sample is no edit', expected: refFree({ [member]: fresh }), actual: refFree({ [member]: sample }) });
      continue;
    }
    for (const [path, run] of Object.entries(paths)) {
      let checks;
      try {
        checks = await run(drawings, member, { base, sample, loaded: foreign && MEMBER_SAMPLES[member].otherModeLoad === 'default' ? base : sample });
      } catch (error) {
        losses.push({ member, mode, path, where: 'threw', expected: refFree(sample), actual: String(error?.stack ?? error) });
        break;
      }
      for (const { where, expected, actual } of checks) {
        if (same(member, actual, expected)) continue;
        losses.push({ member, mode, path, where, expected: refFree({ [member]: expected }), actual: refFree({ [member]: actual }) });
      }
      const live = state.drawingTabs.value.some((tab) => tab.id === drawings.id) ? state.drawingOf(drawings.id) : null;
      if (live && !same(member, captureMember(live, member), sample)) installMember(live, member, sample);
    }
  }
  return losses;
};

const bootBases = async () => {
  if (!bases.plain) for (const name of ['plain', ...Object.keys(CONTEXTS)]) await buildBase(name);
};
// A loss by its member, mode, path and check; a check whose path did not run
// (it threw, Save refused, Load asked) also says what happened instead.
const lossNames = (losses) => losses.map(({ member, mode, path, where, actual }) => `${mode} ${member} ${path}: ${where}${
  where === 'threw' || /"(?:Save refused|Load asked): /.test(actual) ? ` ${actual}` : ''}`);

test('the member samples are the members of the drawing slices and the drawing key classes', () => {
  const sorted = (names) => [...names].sort();
  // The slices' members are the drawing's draft and editor keys, which
  // tests/web/drawing-context.test.mjs checks against the keys `createDrawingState` makes.
  assert.deepEqual(sorted(Object.keys(MEMBER_SAMPLES)), sorted(MEMBERS));
  for (const [member, sample] of Object.entries(MEMBER_SAMPLES)) {
    (sample.needs ?? []).forEach((need) => assert.ok(DRAWING_NEEDS.includes(need) || CONTEXTS[need], `${member} needs ${need}`));
    (sample.modes ?? []).forEach((mode) => assert.ok(MODES.includes(mode), `${member}: ${mode}`));
  }
});

// An object member's sample (`variedLeaves`) changes every leaf of the
// member as a new drawing with the sample's context holds it, but the leaves
// `UNVARIED_LEAVES` names; the member's Load reader keeps every value of the
// app's own lists that the sample draws from, and gives the whole sample back.
test('an object member\'s sample changes every leaf but the ones it names, and its Load reader keeps every listed value and the sample', () => withBootPalettes(state, async () => {
  await bootBases();
  for (const mode of MODES) {
    for (const member of LEAF_VARIED_MEMBERS) {
      const { id } = await startCase(mode, member, !(MEMBER_SAMPLES[member].modes ?? MODES).includes(mode));
      const drawing = state.drawingOf(id);
      const { varied, unchanged, refused, read } = variedLeaves(member, mode, captureMember(drawing, member), drawing);
      assert.deepEqual(refused, [], `${mode} ${member}: the Load reader does not keep these values of the app's own lists, so Load `
        + 'would lose a value the app offers: fix the reader or the list.');
      const reasons = unvariedOf(member, mode);
      assert.deepEqual(unchanged.filter((path) => !Object.hasOwn(reasons, path)), [], `${mode} ${member}: no candidate value of these `
        + 'leaves survives the Load reader, so a reader may now drop them: give the leaf a value the reader keeps, or fix the reader. '
        + 'UNVARIED_LEAVES takes a leaf only with the reason it cannot vary alone.');
      assert.deepEqual(Object.keys(reasons).filter((path) => !unchanged.includes(path)), [], `${mode} ${member}: these UNVARIED_LEAVES `
        + 'entries name no leaf that keeps its value: remove them');
      assert.equal(refFree(read), refFree(varied), `${mode} ${member}: the Load reader keeps the sample`);
    }
  }
}));

// Known defects, each keyed to the one path it affects and the leaves it
// loses there: the guard expects exactly these losses, and the fix of a
// defect removes its entries.
/** @type {Record<string, { defect: string, leaves: readonly string[] }>} */
const KNOWN_LOSSES = {
  // OV-3705: a History restore merges a stored null under an object and keeps
  // the object (`mergeStoredSettings`, services/config-field-readers.js), so
  // the Undo of the step that set `one_circle_left` keeps it.
  ...Object.fromEntries(MODES.map((mode) => [`${mode} adv history: Undo`, { defect: 'OV-3705', leaves: ['one_circle_left'] }])),
  // OV-3707: Reset this drawing keeps the canvas padding, though a drawing's
  // reset resets every slice through its own reset (R15) and Reset returns
  // the drawing's edits to their defaults (PD-OI-090); `resetSettings`
  // (services/reset.js) leaves the canvas slice out.
  ...Object.fromEntries(MODES.map((mode) => [`${mode} canvasPadding resetUndo: Reset this drawing`,
    { defect: 'OV-3707', leaves: ['bottom', 'left', 'right', 'top'] }]))
};
/** The leaves where two captures of a member differ (a leaf under another that differs is in it). */
const differingLeaves = (member, expected, actual) => {
  const [left, right] = [JSON.parse(expected)[member], JSON.parse(actual)[member]];
  const paths = [...new Set([...leafPaths(left), ...leafPaths(right)].map((path) => path.join('.')))]
    .filter((path) => JSON.stringify(leafAt(left, path.split('.'))) !== JSON.stringify(leafAt(right, path.split('.'))));
  return paths.filter((path) => !paths.some((outer) => path.startsWith(`${outer}.`))).sort();
};

test('every member of every slice survives the six state paths in a Circular and a Linear drawing', () => withBootPalettes(state, async () => {
  await bootBases();
  const losses = [];
  for (const mode of MODES) losses.push(...await memberPathLosses(PATHS, mode));
  const named = lossNames(losses);
  assert.deepEqual(named.filter((name) => !KNOWN_LOSSES[name]), [], JSON.stringify(losses, null, 1));
  assert.deepEqual(Object.keys(KNOWN_LOSSES).filter((name) => !named.includes(name)), [],
    'These known losses no longer occur: remove their KNOWN_LOSSES entries and the known-defect case of their defect.');
  for (const loss of losses) {
    assert.deepEqual(differingLeaves(loss.member, loss.expected, loss.actual), KNOWN_LOSSES[lossNames([loss])[0]].leaves,
      `${lossNames([loss])[0]}: the known loss of ${KNOWN_LOSSES[lossNames([loss])[0]].defect} and no other leaf`);
  }
}));

// The guard's own test: over a path that drops one slice or one Session
// domain, it reports that member and no other.
test('the guard reports exactly the member that a broken Duplicate or a broken Save loses', () => withBootPalettes(state, async () => {
  await bootBases();
  const duplicate = app.drawingList.duplicateDrawing;
  // A Duplicate whose copy does not take the canvas slice.
  const duplicateWithoutCanvas = async (...args) => {
    app.drawingList.duplicateDrawing = (id) => {
      const made = duplicate(id);
      DRAWING_SLICES.canvas.reset(state.drawingOf(made.drawingId).slices.canvas);
      return made;
    };
    try {
      return await PATHS.duplicate(...args);
    } finally {
      app.drawingList.duplicateDrawing = duplicate;
    }
  };
  assert.deepEqual(lossNames(await memberPathLosses({ duplicate: duplicateWithoutCanvas }, 'circular', ['canvasPadding', 'annotationSets'])),
    ['circular canvasPadding duplicate: the copy']);
  // A Save whose drawing entry lacks its annotation sets.
  const withoutAnnotationSets = saveThenLoad((document, id) => { delete savedConfigOf(document, id).annotationSets; });
  assert.deepEqual(lossNames(await memberPathLosses({ saveLoad: withoutAnnotationSets }, 'linear', ['annotationSets', 'canvasPadding'])),
    ['linear annotationSets saveLoad: Save then Load']);
}));

// Known defects: each case asserts today's outcome, and fails once the
// defect is fixed; the fix removes the case and the workaround it names.
// OV-3701: Save writes a drawing that has a file (R15), but not one whose
// only change is a card file or the input type (`isUsed` in `exportSession`,
// services/config.js).
test('known defect OV-3701: Save drops a drawing whose only change is a card file or the input type', () => withBootPalettes(state, async () => {
  await bootBases();
  const changes = {
    files: (drawing) => ({ ...captureMember(drawing, 'files'),
      d_color: context.registerFile(new File(['CDS\t#123456\n'], 'ov-3701-colors.tsv', { lastModified: 1 })) }),
    inputType: () => 'gff'
  };
  for (const mode of MODES) {
    for (const [member, change] of Object.entries(changes)) {
      assert.equal((await loadSession(bases.plain.blob)).asked, null);
      const { drawingId: id } = app.drawingList.newDrawing(mode);
      const sample = change(state.drawingOf(id));
      await history.runUndoable(`Sample ${member}`, () => installMember(state.drawingOf(id), member, sample));
      const [check] = await PATHS.saveLoad({ id, peer: bases.plain.peers[mode] }, member, { loaded: sample });
      assert.equal(check.actual, `no drawing ${id}`, `${mode} ${member}: OV-3701 no longer reproduces. Remove this case and the two `
        + '`needs: [\'records\']` marked OV-3701 in tests/web/helpers/drawing-member-samples.mjs.');
    }
  }
}));

// OV-3705: the Undo of choosing one circle leaves what one circle left
// (`adv.one_circle_left`): a History restore merges the step's null under the
// object and keeps it (`mergeStoredSettings`, services/config-field-readers.js).
test('known defect OV-3705: Undo of choosing one circle keeps what one circle left', () => withBootPalettes(state, async () => {
  await bootBases();
  assert.equal((await loadSession(bases.plain.blob)).asked, null);
  const { drawingId: id } = app.drawingList.newDrawing('circular');
  assert.equal(captureMember(state.drawingOf(id), 'adv').one_circle_left, null);
  await history.runUndoable('Change setting', () => chooseRecordLayout(state.drawingOf(id), 'one_circle'));
  assert.equal(captureMember(state.drawingOf(id), 'adv').circular_record_layout, 'one_circle');
  await app.undoHistory();
  const adv = captureMember(state.drawingOf(id), 'adv');
  assert.equal(adv.circular_record_layout, 'grid');
  assert.notEqual(adv.one_circle_left, null, 'OV-3705 no longer reproduces. Remove this case and the OV-3705 entries of '
    + 'KNOWN_LOSSES.');
}));

// The 25 Session config domains (`CURRENT_WRITER_ACTIVE_CONFIG_DOMAINS`) are
// the config domains of a drawing entry (`MODE_SLICE_CONTAINERS.config`, the
// registry's) and three others: the palette's Instant Preview is app UI
// (`ui.paletteInstantPreviewEnabled`), a drawing's CLI options are on its
// entry (`cliOptions`), kept from the Session a Load read, and
// `modeProfiles` is not written (a saved field no slice holds is dropped
// with a Load repair row: tests/web/mode-scoped-session.test.mjs).
const DRAWING_DOMAINS = [...MODE_SLICE_CONTAINERS.config].sort();
const OUTSIDE_DRAWING_CONFIG = ['cliOptions', 'modeProfiles', 'paletteInstantPreviewEnabled'];
const sameValue = (left, right) => refFree({ value: left }) === refFree({ value: right });
/** The places where two saved documents differ: each drawing's config domains, CLI options, and the app UI. */
const savedDifferences = (left, right) => [
  ...left.drawings.flatMap((entry) => {
    const other = right.drawings.find((drawing) => drawing.id === entry.id);
    if (!other) return [`drawing ${entry.id}`];
    return [...DRAWING_DOMAINS.filter((domain) => !sameValue(entry.config[domain], other.config[domain]))
      .map((domain) => `${entry.id} config.${domain}`), ...(sameValue(entry.cliOptions, other.cliOptions) ? [] : [`${entry.id} cliOptions`])];
  }),
  ...Object.keys(left.ui).filter((key) => key !== 'activeDrawingId' && !sameValue(left.ui[key], right.ui[key])).map((key) => `ui.${key}`)
];
/** Save -> Load -> Save of `blob`: the places where the second Save differs from the first, and the document. */
const resaved = async (blob) => {
  assert.equal((await loadSession(blob)).asked, null);
  const { blob: again, refused } = await saveSession();
  assert.equal(refused, null);
  return { differences: savedDifferences(await savedDocument(blob), await savedDocument(again)), document: await savedDocument(again) };
};
/** `config` with the IDs of `drawing`'s File cards and record rows written by their place. */
const ownInputsNamed = (config, drawing) => JSON.parse([...drawing.circularSources, ...drawing.linearSeqs]
  .reduce((text, input, index) => text.split(String(input.uid)).join(`input-${index}`), JSON.stringify(config)));
const writtenDomains = (document) => document.drawings.flatMap((entry) => Object.keys(entry.config));

test('the Session config domains are the drawing config domains, the Instant Preview flag, the CLI options and the mode profiles', () => {
  assert.deepEqual([...CURRENT_WRITER_ACTIVE_CONFIG_DOMAINS].sort(), [...DRAWING_DOMAINS, ...OUTSIDE_DRAWING_CONFIG].sort());
  assert.equal(DRAWING_DOMAINS.length, 22);
});

// Each member's sample, in a drawing of its own, as a Save and Load give it,
// changes the config domains that differ from a drawing with the same context
// and no sample; together they change every drawing config domain, and Save
// -> Load -> Save gives each domain of each drawing back. A sample that needs a fixture's Result
// is compared with that drawing as the fixture base saved it.
test('the member samples change every drawing config domain, and Save then Load then Save keeps each', () => withBootPalettes(state, async () => {
  await bootBases();
  const changed = new Map();
  const documents = [];
  for (const mode of MODES) {
    assert.equal((await loadSession(bases.plain.blob)).asked, null);
    const drawingOfMember = new Map();
    const addDrawing = async () => {
      // New shows the new drawing.
      const { drawingId } = app.drawingList.newDrawing(mode);
      assert.equal(state.activeDrawingId.value, drawingId);
      await bindRecords(state.drawingOf(drawingId));
      return drawingId;
    };
    const reference = await addDrawing();
    for (const member of MEMBERS) {
      const sample = MEMBER_SAMPLES[member];
      if (sample.byContext || !(sample.modes ?? MODES).includes(mode) || (sample.needs ?? []).some((need) => CONTEXTS[need])) continue;
      const drawing = state.drawingOf(await addDrawing());
      await buildDrawingContexts(drawing, (sample.needs ?? []).filter((need) => need !== 'records'));
      sample.prepare?.(drawing);
      installMember(drawing, member, sample.value(mode, drawing, context));
      drawingOfMember.set(member, drawing.id);
    }
    app.selectDrawing(bases.plain.peers[mode]);
    const { blob: installed, refused } = await saveSession();
    assert.equal(refused, null);
    // The samples as a Save and Load give them: Load gives each ring table of
    // the `files` sample its ring row (`reconcileLoadedCircularConservation`),
    // as the upload of a table does.
    assert.equal((await loadSession(installed)).asked, null);
    const { blob } = await saveSession();
    const document = await savedDocument(blob);
    // Each drawing binds its own File card or record row, whose ID its record layout names.
    const configOf = (id) => ownInputsNamed(savedConfigOf(document, id), state.drawingOf(id));
    drawingOfMember.forEach((id, member) => changed.set(`${mode} ${member}`, DRAWING_DOMAINS
      .filter((domain) => !sameValue(configOf(id)[domain], configOf(reference)[domain]))));
    documents.push({ name: `${mode} members`, blob });
  }
  // The samples that need a fixture's Result, all in that drawing (one Save per fixture).
  for (const [name, base] of Object.entries(bases).filter(([name]) => CONTEXTS[name])) {
    const members = MEMBERS.filter((member) => (MEMBER_SAMPLES[member].needs ?? []).includes(name));
    assert.equal((await loadSession(base.blob)).asked, null);
    if (state.activeDrawingId.value !== base.target) assert.deepEqual(app.selectDrawing(base.target), { status: 'ok' });
    const drawing = state.drawingOf(base.target);
    for (const member of members) {
      const sample = MEMBER_SAMPLES[member];
      if (sample.byContext) await CONTEXT_EDITS[name](drawing);
      else installMember(drawing, member, sample.value('linear', drawing, context));
    }
    app.selectDrawing(base.peers.linear);
    const { blob, refused } = await saveSession();
    assert.equal(refused, null);
    const before = savedConfigOf(await savedDocument(base.blob), base.target);
    const after = savedConfigOf(await savedDocument(blob), base.target);
    changed.set(`linear ${members.join(', ')}`, DRAWING_DOMAINS.filter((domain) => !sameValue(before[domain], after[domain])));
    documents.push({ name: `linear ${name}`, blob });
  }
  const covered = new Set([...changed.values()].flat());
  assert.deepEqual(DRAWING_DOMAINS.filter((domain) => !covered.has(domain)), [], JSON.stringify(Object.fromEntries(changed)));
  for (const { name, blob } of documents) {
    const { differences, document } = await resaved(blob);
    assert.deepEqual(differences, [], name);
    for (const domain of writtenDomains(document)) assert.ok(DRAWING_DOMAINS.includes(domain), `${name}: config.${domain}`);
  }
}));

test('Save keeps the Instant Preview flag in the app UI and a drawing\'s CLI options on its entry, and writes no mode profiles', () => withBootPalettes(state, async () => {
  await bootBases();
  const plain = await savedDocument(bases.plain.blob);
  assert.deepEqual(plain.drawings.filter((entry) => Object.hasOwn(entry, 'cliOptions')).map((entry) => entry.id), []);
  // Instant Preview, turned on (it starts off), as its control turns it on.
  assert.equal((await loadSession(bases.plain.blob)).asked, null);
  assert.equal(state.paletteInstantPreviewEnabled.value, false);
  await app.setPaletteInstantPreview(true);
  const { blob: instant } = await saveSession();
  const instantSaved = await savedDocument(instant);
  assert.equal(instantSaved.ui.paletteInstantPreviewEnabled, true);
  assert.equal(writtenDomains(instantSaved).includes('paletteInstantPreviewEnabled'), false);
  const instantAgain = await resaved(instant);
  assert.deepEqual(instantAgain.differences, []);
  assert.equal(state.paletteInstantPreviewEnabled.value, true);
  // CLI options a Session holds on a drawing entry (a CLI-written Session's).
  const cliOptions = { rawArgs: ['--gbk', 'peer.gb', '-o', 'peer'] };
  const withCliOptions = await editedSessionFile(bases.plain.blob, (document) => {
    document.drawings.find((entry) => entry.id === bases.plain.peers.linear).cliOptions = cliOptions;
  });
  const cliAgain = await resaved(withCliOptions);
  assert.deepEqual(cliAgain.differences, []);
  assert.deepEqual(cliAgain.document.drawings.find((entry) => entry.id === bases.plain.peers.linear).cliOptions, cliOptions);
  for (const document of [plain, instantSaved, instantAgain.document, cliAgain.document]) {
    assert.equal(JSON.stringify(document).includes('"modeProfiles"'), false);
  }
}));
