import assert from 'node:assert/strict';
import test from 'node:test';

import {
  createResultPaintRecord as createRecordOf, drawnLabelIntentGone, drawsResult, paintRuleTable, requestPaintOf, withKeptFills
} from '../../gbdraw/web/js/app/result-paint-record.js';

// A Result as the record reads it: an object whose runtime state keys its
// record, which every content version of the Result shares (the preview
// owner's `getResultRuntimeKey`).
const resultNamed = (name) => ({ name, content: `<svg>${name}</svg>`, state: { name } });
const createResultPaintRecord = () => createRecordOf({ resultRuntimeKey: (result) => result.state });

// An editor paint state as app/editor-intent.js reads it (`currentEditorProjectionState`).
/** @param {string} step @returns {import('../../gbdraw/web/js/app/result-paint-record.js').EditorPaintState} */
const paint = (step) => ({
  colors: [{ CDS: step }, `rules ${step}`],
  visibility: `visibility ${step}`,
  strokes: `strokes ${step}`,
  labels: 'labels',
  legendOrder: ''
});
/** @param {ReturnType<typeof paint>} shown @param {ReturnType<typeof paint>} current */
const changed = (shown, current) => ({
  colors: shown.colors[0] !== current.colors[0] || shown.colors[1] !== current.colors[1],
  visibility: shown.visibility !== current.visibility,
  strokes: shown.strokes !== current.strokes
});
const NOTHING = { colors: false, visibility: false, strokes: false };
const EVERY_PAINT_DOMAIN = { colors: true, visibility: true, strokes: true };

// U2BFIX review #1: Save keeps a batch Result's bytes as stored, so a Result
// not displayed since the last edits is saved without them. A Load (or a
// History restore of an artifact) records every Result but the mounted one as
// showing an unknown paint state, so its first display shows every paint
// domain; a Generate draws every Result with the edits.
test('a Result restored from saved bytes shows every paint domain on its first display', () => {
  const edited = paint('edited');
  const [a, b] = ['a', 'b'].map(resultNamed);
  const loaded = createResultPaintRecord();
  loaded.commit([a, b], a, edited);
  const shownB = loaded.display(b, edited);
  assert.deepEqual(changed(shownB, edited), EVERY_PAINT_DOMAIN, 'the Result saved without the edits');
  // U2BFIX2 review L1: also a label intent removed meanwhile, which no
  // label intent left would project. The Legend order needs the Result's
  // generated inventory (U3a).
  assert.notEqual(shownB.labels, edited.labels, 'the labels it was saved with are unknown');
  assert.equal(shownB.legendOrder, edited.legendOrder);
  loaded.shown(b, edited, shownB);
  assert.deepEqual(changed(loaded.display(a, edited), edited), NOTHING, 'the Result mounted by the Load showed them');
  assert.deepEqual(changed(loaded.display(b, edited), edited), NOTHING, 'shown once, the Result holds them');

  const generated = createResultPaintRecord();
  generated.commit([a, b], a, edited, { draws: edited });
  assert.deepEqual(changed(generated.display(b, edited), edited), NOTHING, 'Generate drew the edits');
});

// U2BFIX review #2: a display whose palette and rules projection declined
// (rules changed meanwhile, or a Session operation) showed neither the fills
// nor the visibility; the Result is recorded behind on both, also after it
// leaves, so its next display shows them.
test('a display whose palette and rules projection declined leaves the Result behind on fills and visibility', () => {
  const before = paint('before');
  const after = paint('after');
  const [a, b] = ['a', 'b'].map(resultNamed);
  const record = createResultPaintRecord();
  record.commit([a, b], a, before, { draws: before });
  const shownB = record.display(b, after);
  assert.deepEqual(changed(shownB, after), EVERY_PAINT_DOMAIN);
  record.shown(b, after, shownB, { declined: true });
  record.display(a, after);
  assert.deepEqual(changed(record.display(b, after), after), { colors: true, visibility: true, strokes: false },
    'the declined domains are shown on the next display');
  record.shown(b, after, record.display(b, after));
  record.display(a, after);
  assert.deepEqual(changed(record.display(b, after), after), NOTHING);
  record.depart(after, 'the intent b left with');
  assert.equal(record.takeLeft(b), 'the intent b left with', 'a Result that left after a full display (E1)');
  assert.equal(record.takeLeft(b), null, 'read once');
  assert.deepEqual(changed(record.display(b, after), after), NOTHING, 'a Result that left after a full display');
});

// U2BFIX2 review M3: the displayed Result a declined display left behind is
// known while it stays displayed, so app/editor-intent.js shows the declined
// domains once the blocker clears; a full display or another Result clears it.
test('the domains a declined display left behind are known while the Result stays displayed', () => {
  const before = paint('before');
  const after = paint('after');
  const [a, b] = ['a', 'b'].map(resultNamed);
  const record = createResultPaintRecord();
  record.commit([a, b], a, before, { draws: before });
  assert.equal(record.lacking(a), null);
  const shownB = record.display(b, after);
  record.shown(b, after, shownB, { declined: true });
  assert.deepEqual(record.lacking(b), { colors: before.colors, visibility: before.visibility });
  assert.equal(record.lacking(a), null, 'only the displayed Result');
  record.shown(b, after, after);
  assert.equal(record.lacking(b), null, 'shown in full');
  record.shown(b, after, shownB, { declined: true });
  record.display(a, after);
  assert.equal(record.lacking(b), null, 'a Result no longer displayed waits for its next display');
});

// U2BFIX2 review L2: a display whose projection or compile failed showed none
// of the paint domains; the Result keeps the state it showed before, also
// after it leaves. Its labels follow in the binder's label step.
test('a failed display leaves the Result behind on every paint domain', () => {
  const before = { ...paint('before'), legendOrder: 'order before' };
  const after = paint('after');
  const [a, b] = ['a', 'b'].map(resultNamed);
  const record = createResultPaintRecord();
  record.commit([a, b], a, before, { draws: before });
  const shownB = record.display(b, after);
  record.shown(b, after, shownB, { failed: true });
  record.display(a, after);
  const again = record.display(b, after);
  assert.deepEqual(changed(again, after), EVERY_PAINT_DOMAIN);
  assert.equal(again.legendOrder, before.legendOrder);
  assert.equal(again.labels, after.labels);
});

// OV-346: a label rerender draws the editor state it started from. An edit
// made while it ran is not in its bytes: the displayed Result lacks its fills
// and visibility until they are shown, and another Result shows every domain
// that changed since when it is displayed.
test('a commit drawn from an earlier editor state leaves the edits made meanwhile to be shown', () => {
  const started = paint('started');
  const edited = { ...started, visibility: 'visibility edited meanwhile' };
  const [a, b] = ['a', 'b'].map(resultNamed);
  const record = createResultPaintRecord();
  record.commit([a, b], a, edited, { draws: started });
  assert.deepEqual(record.lacking(a), { colors: started.colors, visibility: started.visibility });
  record.shown(a, edited, edited);
  assert.equal(record.lacking(a), null, 'shown once bound');
  assert.deepEqual(changed(record.display(b, edited), edited), { colors: false, visibility: true, strokes: false },
    'the other Result shows the visibility edited meanwhile');

  const unchanged = createResultPaintRecord();
  unchanged.commit([a, b], a, started, { draws: started });
  assert.equal(unchanged.lacking(a), null, 'no edit while it ran');
});

// D-36, review C S2: whether a label rerender drew the current paint is
// decided when its Result is ready, not when its paint is committed. A rule
// write that lands between the two (in the binder's later steps) leaves the
// record's `lacking` empty, yet the figure lacks the row the write adds. The
// fills (colors, visibility) and the strokes count, so a Legend row restored
// meanwhile does too; labels and the Legend order are projected on any figure.
test('a rerender drew the current paint only while its fills and strokes still match', () => {
  const started = paint('started');
  const a = resultNamed('a');
  const record = createResultPaintRecord();
  record.commit([a], a, started, { draws: started });
  assert.equal(record.lacking(a), null, 'nothing was edited before the paint commit');
  const drewPaint = record.drewPaint;
  assert.equal(drewPaint(started, started), true);
  assert.equal(drewPaint(started, { ...started, colors: [started.colors[0], 'rules written after the paint commit'] }), false,
    'a rule write after the paint commit');
  assert.equal(drewPaint(started, { ...started, colors: [{ CDS: 'another palette' }, started.colors[1]] }), false, 'a palette');
  assert.equal(drewPaint(started, { ...started, visibility: 'hidden meanwhile' }), false, 'a visibility edit');
  assert.equal(drewPaint(started, { ...started, strokes: 'a Legend row restored meanwhile' }), false, 'a stroke edit');
  assert.equal(drewPaint(started, { ...started, labels: 'labels edited', legendOrder: 'order edited' }), true,
    'labels and the Legend order are projected');
  assert.equal(drewPaint(undefined, started), false, 'a binding that names no drawn paint');
});

// OV-1280 (REVIEW-4 N-1): a Result a mode switch shows again names the rule
// table it was last shown with, whatever the drawing did with its rules while
// it was not shown, so its display can ask Python to draw it again; a Result
// restored without its paint state names none.
test('a Result shown again names the rule table it left with, and a restored Result none', () => {
  const rulesOf = (rules) => ({ ...paint('x'), colors: [{ CDS: null }, JSON.stringify([rules, {}, {}])] });
  const withFile = [{ feat: 'CDS', qual: 'gene', val: 'b', color: '#123456', cap: 'Imported', fromFile: true }];
  const [circular, batch, kept] = ['circular', 'batch', 'kept'].map(resultNamed);
  const record = createResultPaintRecord();
  record.commit([circular], circular, rulesOf(withFile), { draws: rulesOf(withFile) });
  record.depart(rulesOf(withFile), 'the intent it left with');
  // The file card's Remove while another mode is shown, then the switch back.
  assert.deepEqual(paintRuleTable(record.display(circular, rulesOf([]))), withFile);
  record.commit([batch], batch, rulesOf([]));
  record.commit([batch, kept], batch, rulesOf([]));
  assert.equal(paintRuleTable(record.display(kept, rulesOf([]))), null);
});

// OV-2293: what only a rerender shows of a label edit the Result was drawn with.
test('a Result lacks what its drawing draws when a drawn label On, Off, or text left, and not for another value', () => {
  const rows = (...labelRows) => JSON.stringify([labelRows.map(([key, visibility, text]) => ['record', key, visibility, text, null]), {}]);
  const gone = (drawn, current) => drawnLabelIntentGone(drawn, '', current, () => true);
  assert.equal(gone(rows(['a', 'off', null]), rows()), true, 'Off returned to Default, its row gone');
  assert.equal(gone(rows(['a', 'on', null]), rows(['a', null, 'text'])), true, 'On returned to Default, its row kept');
  assert.equal(gone(rows(['a', null, 'edit']), rows()), true, 'the text left, its row gone');
  assert.equal(gone(rows(['a', null, 'edit']), rows(['a', 'off', null])), true, 'the text left, its row kept');
  assert.equal(gone(rows(['a', 'off', null]), rows(['a', 'on', null])), false, 'Off to On: the label projection asks');
  assert.equal(gone(rows(['a', null, 'edit']), rows(['a', null, 'other'])), false, 'another text: the label projection shows it');
  assert.equal(gone(rows(), rows(['a', 'off', 'edit'])), false, 'a row added: the projection shows it');
  const bulk = (texts) => JSON.stringify([[], texts]);
  assert.equal(gone(bulk({ alpha: 'edit' }), bulk({})), true, 'a bulk text left');
  assert.equal(gone(bulk({ alpha: 'edit' }), bulk({ alpha: 'other' })), true, 'a bulk text changed: its bytes hold no source text');
  assert.equal(gone(bulk({}), bulk({ alpha: 'edit' })), false, 'a bulk text added');
});

test('a Result asks only for a drawn label edit that left since it was last shown, on a feature it draws', () => {
  const rows = (...labelRows) => JSON.stringify([labelRows.map(([key, visibility, text]) => ['record', key, visibility, text, null]), {}]);
  const edited = rows(['a', null, 'edit']);
  const empty = rows();
  const draws = () => true;
  assert.equal(drawnLabelIntentGone(edited, edited, empty, draws), true, 'it left while the Result was not shown');
  assert.equal(drawnLabelIntentGone(edited, '', empty, draws), true, 'unknown: a restored Result');
  assert.equal(drawnLabelIntentGone(edited, empty, empty, draws), false, 'the Result followed it, or asked once');
  assert.equal(drawnLabelIntentGone(rows(['a', null, 'x'], ['b', null, 'y']), rows(['b', null, 'y']), empty, draws), true,
    'another drawn edit left since');
  assert.equal(drawnLabelIntentGone(edited, edited, empty, () => false), false, 'a feature this Result does not draw');
  const drawnFeatures = [];
  drawnLabelIntentGone(edited, edited, empty, (recordKey, featureId) => { drawnFeatures.push([recordKey, featureId]); return true; });
  assert.deepEqual(drawnFeatures, [['record', 'a']], 'found by identity');
  const bulk = (texts) => JSON.stringify([[], texts]);
  assert.equal(drawnLabelIntentGone(bulk({ alpha: '' }), bulk({ alpha: '' }), bulk({}), () => false), true, 'a blank bulk text left: the bytes name no feature');
  assert.equal(drawnLabelIntentGone(bulk({ alpha: '' }), bulk({}), bulk({}), draws), false, 'a bulk text asked once');
});

test('a Result keeps the state its bytes were drawn with while the state it shows follows the edits', () => {
  const record = createResultPaintRecord();
  const paint = (labels) => ({ colors: [null, ''], visibility: '[]', strokes: '', labels, legendOrder: '' });
  const held = paint('held');
  const [generated, kept] = ['generated', 'kept'].map(resultNamed);
  record.commit([generated], generated, held, { draws: held });
  record.display(generated, paint('edited'));
  record.depart(paint('edited'), 'left');
  assert.equal(record.drawn(generated).labels, 'held');
  assert.equal(record.display(generated, paint('reset')).labels, 'edited');
  record.commit([generated, kept], generated, paint('shown drawing'), { drawingStateOf: () => paint('its drawing') });
  assert.equal(record.drawn(generated).labels, 'held', 'a Result already known keeps its record');
  assert.equal(record.drawn(kept).labels, 'its drawing', 'a restored Result not displayed reads its own drawing');
  assert.equal(record.display(kept, paint('its drawing')).labels, '', 'and shows every domain on its first display');
  record.commit([kept], kept, paint('its drawing'));
  assert.equal(record.drawn(generated).labels, 'held', 'a Result that left the drawings keeps its record');
  assert.equal(record.drawn({ ...generated, content: '<svg>edited</svg>' }).labels, 'held',
    'every content version of a Result shares its record');
});

// OV-3730: only a Generate, a label rerender, and a committed candidate draw
// their Results; a checkpoint restore (a Reset Settings, a rule preset Undo)
// mounts Results whose batch members were kept as they were, so they show
// every paint domain on their first display, as after a Load.
test('only the phases Python drew record their Results as drawn from the state of the bind', () => {
  ['generate', 'label-reflow', 'target-record-transform'].forEach((phase) => assert.equal(drawsResult(phase), true, phase));
  [
    'session-load', 'passive-preview', 'rollback-restoration', 'history-finalization-rollback',
    'target-history-finalization-rollback', 'result-selection', 'preview'
  ].forEach((phase) => assert.equal(drawsResult(phase), false, phase));
  const current = paint('current');
  const [a, b] = ['a', 'b'].map(resultNamed);
  const restored = createResultPaintRecord();
  restored.commit([a, b], a, current, { draws: null, drawingStateOf: () => paint('drawing') });
  assert.deepEqual(changed(restored.display(b, current), current), EVERY_PAINT_DOMAIN, 'a checkpoint restore kept its batch bytes');
  assert.equal(restored.drawn(b).labels, paint('drawing').labels, 'drawn with the state of its drawing');
  assert.deepEqual(changed(restored.display(a, current), current), NOTHING, 'the mounted Result shows the restored state');
  const generated = createResultPaintRecord();
  generated.commit([a, b], a, current, { draws: current });
  assert.deepEqual(changed(generated.display(b, current), current), NOTHING, 'a Generate drew every Result');
});

// A Generate Undo brings the same Result objects back: their records come
// back with them, so a batch Result shows only what changed since it was last
// shown, and a Redo finds the Generate's Results as it left them. The Result
// displayed when the Undo restored the editor state shows an unknown state.
test('a Result an Undo of a Generate brings back keeps its record', () => {
  const before = paint('before');
  const after = paint('after');
  const [r1, r2] = ['r1', 'r2'].map(resultNamed);
  const [g1, g2] = ['g1', 'g2'].map(resultNamed);
  const record = createResultPaintRecord();
  record.commit([r1, r2], r1, before, { draws: before });
  record.commit([g1, g2], g1, after, { draws: after });
  const drawnG1 = record.drawn(g1);
  record.shown(g2, after, record.display(g2, after));
  record.display(g1, after);
  // The Undo restores the earlier Results and their drawing's state.
  record.commit([r1, r2], r1, before);
  assert.deepEqual(changed(record.display(r2, before), before), NOTHING, 'r2 was shown with the restored state');
  // The Redo brings the Generate's Results back.
  record.commit([g1, g2], g1, after);
  assert.equal(record.drawn(g1), drawnG1, 'the same record');
  assert.deepEqual(changed(record.display(g2, after), after), NOTHING, 'g2 shows what it showed when it left');
  record.commit([r1, r2], r1, before);
  assert.deepEqual(changed(record.display(g2, after), after), EVERY_PAINT_DOMAIN,
    'the Result displayed when an Undo restored the state shows an unknown state');
});

// OV-3731: a Result that leaves a drawing lacking the domains a declined
// display left behind is shown again through the projection, not as E1.
test('a Result that leaves lacking a domain keeps no intent to be shown again with', () => {
  const before = paint('before');
  const after = paint('after');
  const [a, b] = ['a', 'b'].map(resultNamed);
  const record = createResultPaintRecord();
  record.commit([a, b], a, before, { draws: before });
  record.shown(b, after, record.display(b, after), { declined: true });
  record.depart(after, 'the intent it left with');
  assert.equal(record.takeLeft(b), null, 'no E1 rebind');
  assert.deepEqual(changed(record.display(b, after), after), { colors: true, visibility: true, strokes: false },
    'its next display shows the declined domains');
  record.shown(b, after, after);
  record.depart(after, 'the intent it left with');
  assert.equal(record.takeLeft(b), 'the intent it left with', 'shown in full, it is shown again as it left (E1)');
});

// OV-1983, OV-2413: the label owner reads the Label visibility On and Off the
// displayed Result was drawn with; a Generate that moves a record's rows to
// its new key drew the moved rows (read as its request was built),
// which its Results' record takes when the drawing takes them.
test('the Label visibility the displayed Result was drawn with follows a Generate key move', () => {
  const labelsOf = (from) => JSON.stringify([[[from, 'f1', 'on', null, null], [from, 'f2', null, 'text', null], ['kept', 'f3', 'off', null, null]], {}]);
  const bound = { ...paint('drawn'), labels: labelsOf('old') };
  const moved = { ...paint('drawn'), labels: labelsOf('new') };
  const [g, other] = ['g', 'other'].map(resultNamed);
  const record = createResultPaintRecord();
  assert.equal(record.displayedDrawnLabelVisibility(), null, 'no Result displayed');
  record.commit([g, other], g, bound, { draws: bound });
  const key = (recordKey, featureId) => JSON.stringify([recordKey, featureId]);
  assert.deepEqual([...record.displayedDrawnLabelVisibility()], [[key('old', 'f1'), 'on'], [key('kept', 'f3'), 'off']]);
  record.moveDrawn([g, other], moved);
  assert.deepEqual(new Map(record.displayedDrawnLabelVisibility()), new Map([[key('new', 'f1'), 'on'], [key('kept', 'f3'), 'off']]));
  assert.equal(record.drawn(other), moved, 'every committed Result takes the moved rows');
  record.depart(moved, 'left');
  assert.equal(record.displayedDrawnLabelVisibility(), null, 'the Result left');
});

// The paint a Generate drew is read as its request is built, with the palette its commit applies, and
// names the fills the Generate keeps; the input summary's staleness reads the visibility it drew.
test('a Generate names the paint it drew as its request is built', () => {
  const drawing = {
    manualSpecificRules: [], featureColorOverrides: { k: { color: '#111111' }, gone: { color: '#222222' } }, legendColorOverrides: {},
    featureOverrides: {}, featureVisibilityManualRules: [], featureStrokeOverrides: {}, legendStrokeOverrides: {},
    deletedLegendEntries: { value: [] }, labelTextBulkOverrides: {}, legendEntries: { value: [] }
  };
  const request = requestPaintOf(drawing, Object.freeze({ name: 'request', colors: { CDS: '#abcdef' } }), { value: [] });
  assert.equal(request.palette.name, 'request');
  assert.deepEqual(request.bound.colors[0], { CDS: '#abcdef' }, 'the palette of the request, not the one shown when it binds');
  drawing.featureColorOverrides = {};
  assert.match(request.bound.colors[1], /#222222/, 'read as the request is built: a later edit is an edit since');
  const kept = withKeptFills(request.bound, { k: { color: '#111111' } });
  assert.deepEqual(JSON.parse(kept.colors[1])[1], { k: { color: '#111111' } }, 'the fills the Generate keeps');
  assert.equal(JSON.parse(kept.colors[1])[0].length, 0, 'its rules stay');
  const [g] = ['g'].map(resultNamed);
  const record = createResultPaintRecord();
  assert.equal(record.displayedDrewOtherVisibility(kept.visibility), false, 'no Result displayed yet');
  record.commit([g], g, paint('now'), { draws: kept });
  assert.equal(record.displayedDrewOtherVisibility(kept.visibility), false);
  assert.equal(record.displayedDrewOtherVisibility('another visibility'), true, 'a live hide no render drew');
  assert.equal(record.drawn(g), kept, 'the record holds what the Generate drew, not the state at its bind');
});

// OV-3732: Duplicate gives the copy's Results records of their own, which
// start as the source Results'; the displayed source Result followed every
// live edit (the copied state) but what a declined display left it lacking.
test('a copied Result has a record of its own that starts as its source\'s', () => {
  const [a, b, aCopy, bCopy] = ['a', 'b', 'a copy', 'b copy'].map(resultNamed);
  const record = createResultPaintRecord();
  const drawn = paint('drawn');
  const current = paint('current');
  record.commit([a, b], a, drawn, { draws: drawn });
  record.copy([a, b], [aCopy, bCopy], current);
  assert.equal(record.drawn(aCopy), drawn);
  assert.deepEqual(changed(record.display(aCopy, current), current), NOTHING, 'the displayed source showed the copied state');
  assert.deepEqual(record.display(bCopy, current), drawn, 'a source not displayed is copied as it is');
  record.shown(bCopy, paint('copy edited'), drawn);
  assert.deepEqual(record.display(b, current), drawn, 'the source keeps its own record');
});


// OV-3730: the editor intent's commit reads from the bind's phase whether
// Python drew the Results it records (`drawsResult`). A checkpoint restore
// keeps a batch Result's bytes, drawn with the state of the drawing holding it
// and shown unknown; a label rerender drew every Result from its request.
test('the committed editor state records a Result as drawn only in a phase that draws (OV-3730)', async () => {
  // app/editor-intent.js reads `window.Vue.watch` when it loads.
  globalThis.window ??= {};
  globalThis.window.Vue ??= { watch() {} };
  const { createDisplayedResultProjection } = await import('../../gbdraw/web/js/app/editor-intent.js');
  const [a, b] = ['a', 'b'].map(resultNamed);
  const [first, second] = [{ id: 'd1', step: 'current' }, { id: 'd2', step: 'second drawing' }];
  const drawnPaint = paint('drawn');
  const commitIn = (phase, bindingOptions = { drawnPaint }) => {
    const record = createResultPaintRecord();
    const unused = () => { throw new Error('not read by a commit'); };
    const { rememberCommittedEditorState } = createDisplayedResultProjection({
      liveResults: () => [a, b],
      drawingList: () => [first, second],
      resultsOf: (drawingId) => (drawingId === 'd1' ? [a] : [b]),
      resultIdentity: (result) => result.name,
      currentEditorProjectionState: (drawing) => paint(drawing.step),
      displayedIntentSignature: unused,
      resultDrawsFeature: unused,
      compileDisplayedResultOperations: unused,
      projectMountedEditorIntent: unused,
      showEditorIntent: unused,
      displayBlocked: unused,
      showDeclinedPaint: unused,
      LEGEND_STRUCTURE_DOMAINS: [],
      STROKE_DOMAINS: [],
      resultPaintCommit: record.commit,
      resultPaintTakeLeft: record.takeLeft,
      resultPaintRebind: record.rebind,
      resultPaintDisplay: record.display,
      resultPaintDrawn: record.drawn,
      resultPaintShown: record.shown
    });
    rememberCommittedEditorState(first, {
      root: null, result: a, resultIndex: 0, resultIdentity: 'a', rootGeneration: 1, phase, bindingOptions
    });
    return record;
  };
  const current = paint('current');
  const restored = commitIn('passive-preview');
  assert.equal(restored.drawn(b)?.visibility, paint('second drawing').visibility, 'drawn with its own drawing\'s state');
  assert.deepEqual(changed(restored.display(b, current), current), EVERY_PAINT_DOMAIN, 'shown unknown');
  assert.equal(restored.drawn(a)?.visibility, current.visibility, 'the displayed Result shows the restored state');
  const redrawn = commitIn('label-reflow');
  assert.equal(redrawn.drawn(b)?.visibility, drawnPaint.visibility, 'drawn from the rerender\'s request');
  assert.equal(redrawn.drawn(a)?.visibility, drawnPaint.visibility);
  // Every phase that draws names the paint it drew; none falls back to the state at its bind.
  assert.throws(() => commitIn('generate', {}), (error) => error.code === 'RESULT_INVALID'
    && error.context.reason === 'DRAWN_PAINT_MISSING');
});

// A held figure shows with the palette the restored intent left, which the
// restore clones: held colors equal to it are that palette, so the Result's
// next display asks for no fills; held colors that differ (a palette a later
// open step applied as the figure was held) stay, so it asks for them.
test('a held figure takes the applied palette when its held colors equal it', () => {
  const [a, other] = ['a', 'other'].map(resultNamed);
  const record = createResultPaintRecord();
  const before = paint('before');
  record.commit([a, other], a, before, { draws: before });
  const picked = { ...before, colors: /** @type {typeof before.colors} */ ([{ CDS: 'picked' }, before.colors[1]]) };
  const held = record.heldPaint([a, other], picked, false);
  const restored = { ...before, colors: /** @type {typeof before.colors} */ ([{ CDS: 'before' }, before.colors[1]]) };
  record.depart(restored, 'signature');
  record.arrive([a, other], held, restored.colors[0]);
  assert.equal(record.display(a, restored).colors[0], picked.colors[0], 'other colors stay: the display paints the restored ones');
  record.depart(restored, 'signature');
  record.arrive([a, other], record.heldPaint([a, other], before, false), restored.colors[0]);
  const shownA = record.display(a, restored);
  assert.equal(shownA.colors[0], restored.colors[0], 'the same colors are the applied palette');
  assert.deepEqual(changed(shownA, restored), NOTHING, 'and ask for no projection');
  assert.equal(record.display(other, restored).colors[0], restored.colors[0], 'for a batch Result too');
});

// OV-2436, review H1: a figure History holds keeps what its Results' bytes
// show, and its records say so again when it shows: the figure a rerender
// replaced was captured with the edit shown on it, so a second Undo, after a
// Redo left it with the intent before the edit, still shows the step undone.
// A rerender's figure, captured as its Legend is read, shows the paint it was
// drawn with (an edit made while it ran is shown on it later), and a batch
// Result not displayed shows what its record says.
test('a held figure shows again with the paint its bytes show', () => {
  const [a, other, b] = ['a', 'other', 'b'].map(resultNamed);
  const record = createResultPaintRecord();
  const [before, edited, meanwhile] = ['before', 'edited', 'meanwhile'].map(paint);
  record.commit([a, other], a, before, { draws: before });
  // The edit is shown on A, which the rerender then replaces.
  const heldA = record.heldPaint([a, other], edited, false);
  assert.deepEqual(heldA, [edited, before], 'the displayed Result shows the current state, the other its record');
  record.commit([b], b, meanwhile, { draws: edited });
  assert.deepEqual(record.heldPaint([b], meanwhile, true), [edited], 'a rerender figure shows the paint it was drawn with');
  const lacking = { ...meanwhile, colors: edited.colors, visibility: edited.visibility };
  assert.deepEqual(record.heldPaint([b], meanwhile, false), [lacking], 'once bound, the current state but the fills it lacks');
  // Undo shows A; Redo leaves it with the intent before the edit.
  record.depart(meanwhile, 'signature before');
  record.arrive([a, other], heldA, before.colors[0]);
  assert.deepEqual(record.display(a, before), edited, 'Undo projects the step undone onto A');
  record.depart(before, 'signature before');
  // Undo again: A's record says it left with the current intent, but the held bytes show the edit.
  record.arrive([a, other], heldA, before.colors[0]);
  assert.equal(record.takeLeft(a), null, 'the held figure left no drawing');
  assert.deepEqual(record.display(a, before), edited, 'the second Undo projects the step undone again');
  assert.deepEqual(record.heldPaint([b], before, false), [lacking],
    'a Result not displayed shows what it left with, the fills it lacked included');
});
