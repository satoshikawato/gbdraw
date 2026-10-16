// Step 3a (gbdraw/web/CLAUDE.md, Runtime data flow; CW-04 and CW-06 in
// docs/internal/ARCHITECTURE_FITNESS_FUNCTION_RATCHET.md): a Generate whose
// feature catalog build did not change keeps the catalog object the page
// admitted; every other Generate sends and admits it. Each case below is one
// invalidation row of the CW-06 table in the step 3a pull request. Step 3b S:
// the page sanitizes each SVG children text once per revision; the last test
// holds the browser rows of its invalidation table.
const { test, expect } = require('@playwright/test');
const { readFileSync } = require('node:fs');
const { gunzipSync } = require('node:zlib');
const { join } = require('node:path');
const { openApp, evaluateWithRetainedPromise } = require('./helpers/app-lifecycle.cjs');
const { generate, switchMode, download } = require('./helpers/mode-transition.cjs');

const HMMT = readFileSync(join(__dirname, '../test_inputs/HmmtDNA.gbk'), 'utf8');
const SESSIONS = 'gbdraw/web/gallery/sessions';

// Every Generate reply's pieces, with the children texts it sent decoded, as
// the Worker posts them (registered before the page's own listener).
const captureGenerateReplies = () => {
  const NativeWorker = window.Worker;
  window.__PIECES__ = [];
  window.Worker = class extends NativeWorker {
    constructor(...args) {
      super(...args);
      this.addEventListener('message', ({ data }) => {
        if (data?.type !== 'run' || !data.ok || !Array.isArray(data.results?.results)) return;
        const decoder = new TextDecoder();
        window.__PIECES__.push(data.results.results.map(({ document, segments }) => {
          const bytes = segments ? new Uint8Array(segments.buffer ?? segments, segments.byteOffset ?? 0, segments.byteLength) : null;
          return {
            document,
            children: document.pieces.map((piece) => (piece.kind === 'content' && piece.offset !== undefined
              ? decoder.decode(bytes.subarray(piece.offset, piece.offset + piece.length)) : null))
          };
        }));
      });
    }
  };
};

const start = async (browser, { captureReplies = false } = {}) => {
  const context = await browser.newContext({ viewport: { width: 1600, height: 1000 } });
  await context.grantPermissions(['clipboard-read', 'clipboard-write']);
  const page = await context.newPage();
  page.externalRequests = [];
  await context.route('**/*', (route) => {
    if (new URL(route.request().url()).hostname === '127.0.0.1') return route.continue();
    page.externalRequests.push(route.request().url());
    return route.abort();
  });
  await page.addInitScript(() => {
    window.__REUSE__ = { events: [], metrics: [] };
    window.__GBDRAW_TEST_HOOKS__ = {
      onSessionLifecycleEvent: (event) => window.__REUSE__.events.push(event),
      onStructuralMetric: (metric) => window.__REUSE__.metrics.push(metric)
    };
  });
  if (captureReplies) await page.addInitScript(captureGenerateReplies);
  page.on('dialog', (dialog) => dialog.dismiss());
  await openApp(page);
  return page;
};

const loadSession = async (page, file) => {
  await page.locator('input[accept^=".json,"]').setInputFiles(file);
  await expect.poll(() => page.evaluate(() => !window.__GBDRAW_APP__.sessionImportPending
    && window.__GBDRAW_APP__.extractedFeatures.length > 0), { timeout: 180_000 }).toBe(true);
};

const mark = (page) => page.evaluate(() => ({
  events: window.__REUSE__.events.length,
  metrics: window.__REUSE__.metrics.length,
  undo: window.__GBDRAW_HISTORY__.getUndoCount()
}));

// What one render did: whether it named a held revision, whether Python sent
// the catalog, and whether the page admitted or kept it.
const since = (page, start) => page.evaluate((start) => {
  const events = window.__REUSE__.events.slice(start.events);
  const metrics = {};
  for (const { name, value } of window.__REUSE__.metrics.slice(start.metrics)) {
    metrics[name] = (metrics[name] || 0) + Number(value || 0);
  }
  const transport = events.filter((event) => event.name === 'result-transport-ready');
  return {
    names: events.map((event) => event.name),
    transports: transport.map((event) => ({ held: event.heldFeatureCatalogRevision, catalog: event.featureCatalog })),
    segments: transport.map((event) => event.heldSegmentCount),
    admission: metrics.featureCatalogAdmissionCount || 0,
    reuse: metrics.featureCatalogReuseCount || 0,
    undo: window.__GBDRAW_HISTORY__.getUndoCount() - start.undo
  };
}, start);

const observeGenerate = async (page) => {
  const before = await mark(page);
  await generate(page);
  return since(page, before);
};

const keepShown = (page, key) => page.evaluate(async (key) => {
  const { state } = await import('./js/state.js');
  window[key] = window.Vue.toRaw(state.featureCatalog.value);
  return JSON.stringify(window[key]);
}, key);
const shownIs = (page, key) => page.evaluate(async (key) => {
  const { state } = await import('./js/state.js');
  return Boolean(window[key]) && window.Vue.toRaw(state.featureCatalog.value) === window[key];
}, key);
const shownJson = (page) => page.evaluate(async () => {
  const { state } = await import('./js/state.js');
  return JSON.stringify(window.Vue.toRaw(state.featureCatalog.value));
});
// One History step that sets app fields (`adv.*`, `form.*`) or renames a group.
const undoable = (page, label, { fields = {}, rename = null }) => evaluateWithRetainedPromise(page, async ({ label, fields, rename }) => {
  const app = window.__GBDRAW_APP__;
  await window.__GBDRAW_HISTORY__.runUndoable(label, () => {
    for (const [path, value] of Object.entries(fields)) {
      const [owner, key] = path.split('.');
      app[owner][key] = value;
    }
    if (rename) app.setOrthogroupNameOverride(rename.id, rename.name);
  });
  await window.Vue.nextTick();
}, { label, fields, rename });
const history = (page, action) => evaluateWithRetainedPromise(page, async (action) => {
  await window.__GBDRAW_HISTORY__[action]();
  await window.Vue.nextTick();
}, action);

const SENT = { held: false, catalog: 'sent' };
const HELD_SENT = { held: true, catalog: 'sent' };
const KEPT = { held: true, catalog: 'unchanged' };

test('a Generate keeps the held catalog only while its build is current', async ({ browser }, testInfo) => {
  test.setTimeout(420_000);
  const page = await start(browser);
  try {
    await page.evaluate(async (content) => {
      const app = window.__GBDRAW_APP__;
      app.files.c_gb = new File([content], 'HmmtDNA.gbk', { type: 'text/plain' });
      // The title has its place from the start: moving it there moves the
      // records too, which builds the catalog items again.
      app.adv.plot_title_position = 'top';
      await window.Vue.nextTick();
    }, HMMT);

    console.log('STEP Pw-1', Date.now());
    // Pw-1: the first Generate holds nothing; the catalog is sent and admitted.
    const cold = await observeGenerate(page);
    expect(cold.transports).toEqual([SENT]);
    expect(cold.names).toContain('catalog.admission-started');
    expect(cold.admission).toBe(1);
    await keepShown(page, 'c1');

    console.log('STEP Pw-2', Date.now());
    // Pw-2: no change: Python reuses the build; the page keeps the object and adds no History step.
    const warm = await observeGenerate(page);
    expect(warm.transports).toEqual([KEPT]);
    expect(warm.names).toContain('catalog.reused');
    expect({ admission: warm.admission, reuse: warm.reuse, undo: warm.undo }).toEqual({ admission: 0, reuse: 1, undo: 0 });
    expect(await shownIs(page, 'c1')).toBe(true);

    console.log('STEP Pw-3', Date.now());
    // Pw-3: a render-only change (the plot title text) draws a new figure with the same catalog.
    await undoable(page, 'Set the plot title', { fields: { 'form.plot_title': 'Reuse title' } });
    const retitled = await observeGenerate(page);
    expect(retitled.transports).toEqual([KEPT]);
    expect(retitled.undo).toBe(1);
    expect(await shownIs(page, 'c1')).toBe(true);
    const keptJson = await shownJson(page);

    console.log('STEP Pw-7', Date.now());
    // Pw-7: a new Worker has a new epoch: the same request sends the catalog,
    // and the one the page kept equals a fresh Generate's (D-31).
    await page.evaluate(async () => {
      (await import('./js/services/diagram-generation.js')).disposeDiagramGenerationWorker();
    });
    const fresh = await observeGenerate(page);
    expect(fresh.transports).toEqual([HELD_SENT]);
    expect(fresh.admission).toBe(1);
    expect(await shownJson(page)).toBe(keptJson);
    const c1Json = await keepShown(page, 'c1');

    console.log('STEP Pw-5', Date.now());
    // Pw-5: a drawn-feature change sends a new catalog; Undo and Redo restore
    // each object; only the newest one's builds are current.
    await undoable(page, 'Set the block stroke width', { fields: { 'adv.block_stroke_width': 3 } });
    const changed = await observeGenerate(page);
    expect(changed.transports).toEqual([HELD_SENT]);
    expect(changed.admission).toBe(1);
    await keepShown(page, 'c2');
    await history(page, 'undo');
    expect(await shownIs(page, 'c1')).toBe(true);
    await history(page, 'redo');
    expect(await shownIs(page, 'c2')).toBe(true);
    const again = await observeGenerate(page);
    expect(again.transports).toEqual([KEPT]);
    expect(await shownIs(page, 'c2')).toBe(true);
    await history(page, 'undo');
    await history(page, 'undo');
    expect(await shownIs(page, 'c1')).toBe(true);
    const older = await observeGenerate(page);
    expect(older.transports).toEqual([HELD_SENT]);
    expect(older.admission).toBe(1);
    expect(await shownJson(page)).toBe(c1Json);

    console.log('STEP Pw-6', Date.now());
    // Pw-6: each mode shows its own slot; the Linear Generate evicted the
    // Circular builds, so the Circular catalog is sent again. SVG segments
    // are content hashes on each slot's Results: Linear holds none, and the
    // Circular slot still names its own.
    await switchMode(page, 'linear');
    await evaluateWithRetainedPromise(page, async (content) => {
      const app = window.__GBDRAW_APP__;
      app.lInputType = 'gb';
      if (!app.linearSeqs.length) app.addLinearSeq();
      app.setLinearSeqPrimaryFile(0, 'gb', new File([content], 'linear.gbk', { type: 'text/plain', lastModified: 1 }));
      await app.refreshLinearRecordSelectors();
    }, HMMT);
    const linear = await observeGenerate(page);
    expect(linear.transports).toEqual([SENT]);
    expect(linear.segments).toEqual([0]);
    await switchMode(page, 'circular');
    const back = await observeGenerate(page);
    expect(back.transports).toEqual([HELD_SENT]);
    expect(back.admission).toBe(1);
    expect(back.segments[0]).toBeGreaterThan(0);

    console.log('STEP Pw-9', Date.now());
    // Pw-9: a failed Load rolls back to the same object, still current.
    await keepShown(page, 'c3');
    const rollback = await evaluateWithRetainedPromise(page, async () => {
      const config = await import('./js/services/config.js');
      const file = new File([await (await fetch('/gbdraw/web/gallery/sessions/HmmtDNA_basic_circular.gbdraw-session.json')).arrayBuffer()], 'replacement.json');
      const result = await config.importSession({ target: { files: [file], value: 'selected' } }, {
        beforePreviewMount() { throw new Error('Catalog reuse rollback probe'); }
      });
      return result.status;
    });
    expect(rollback).toBe('error');
    expect(await shownIs(page, 'c3')).toBe(true);
    const afterRollback = await observeGenerate(page);
    expect(afterRollback.transports).toEqual([KEPT]);
    expect(await shownIs(page, 'c3')).toBe(true);

    console.log('STEP Pw-8', Date.now());
    // Pw-8: Save Session writes the kept catalog; a loaded catalog has no revision.
    // An untitled Session asks for a title first, and this page dismisses dialogs.
    await page.evaluate(async () => {
      (await import('./js/state.js')).state.sessionTitle.value = 'reuse';
    });
    const bytes = await download(page, 'Save Session', testInfo.outputPath('reuse.gbdraw-session.json.gz'));
    const saved = JSON.parse((bytes[0] === 0x1f ? gunzipSync(bytes) : bytes).toString('utf8'));
    expect(JSON.stringify(saved.editorState.featureCatalog)).toBe(await shownJson(page));
    await loadSession(page, testInfo.outputPath('reuse.gbdraw-session.json.gz'));
    const loaded = await observeGenerate(page);
    expect(loaded.transports).toEqual([SENT]);
    expect(loaded.admission).toBe(1);
    expect(page.externalRequests).toEqual([]);
  } finally {
    await page.context().close();
  }
});

test('a rerender of a loaded Session names the catalog under the names its Results keep', async ({ browser }) => {
  test.setTimeout(300_000);
  const page = await start(browser);
  try {
    await loadSession(page, `${SESSIONS}/HmmtDNA_basic_circular.gbdraw-session.json`);
    const savedName = await page.evaluate(() => window.__GBDRAW_APP__.results[0].name);
    const reflow = async () => {
      const before = await mark(page);
      await page.evaluate(async () => {
        const { state } = await import('./js/state.js');
        state.labelReflowForceRequestSeq.value += 1;
      });
      await expect.poll(async () => (await since(page, before)).transports.length, { timeout: 180_000 }).toBe(1);
      await expect.poll(() => page.evaluate(async () =>
        (await import('./js/state.js')).state.labelReflowProcessing.value), { timeout: 180_000 }).toBe(false);
      return since(page, before);
    };
    console.log('STEP Pw-4', Date.now());
    // Pw-4: the first reflow after Load holds nothing; Python names the items
    // with the engine's Result names, and the page keeps the saved names (OV-136).
    const first = await reflow();
    expect(first.transports).toEqual([SENT]);
    await keepShown(page, 'reflowed');
    expect(await page.evaluate(() => window.__GBDRAW_APP__.results[0].name)).toBe(savedName);
    const second = await reflow();
    expect(second.transports).toEqual([KEPT]);
    expect(second.names).toContain('catalog.reused');
    expect(await shownIs(page, 'reflowed')).toBe(true);
    expect(await page.evaluate(() => window.__GBDRAW_APP__.results[0].name)).toBe(savedName);
    // A Generate takes the engine's names: the held items do not match, so it is sent.
    const generated = await observeGenerate(page);
    expect(generated.transports).toEqual([SENT]);
    expect(generated.admission).toBe(1);
    const engineName = await page.evaluate(() => window.__GBDRAW_APP__.results[0].name);
    expect(engineName).not.toBe(savedName);
    expect(page.externalRequests).toEqual([]);
  } finally {
    await page.context().close();
  }
});

test('the similarity-group copy feedback clears when a Generate replaces the shown Results, as before', async ({ browser }) => {
  test.setTimeout(300_000);
  const page = await start(browser);
  try {
    await loadSession(page, `${SESSIONS}/majanivirus_orthogroup.gbdraw-session.json.gz`);
    const first = await observeGenerate(page);
    expect(first.transports).toEqual([SENT]);
    await page.evaluate(() => {
      const app = window.__GBDRAW_APP__;
      app.selectOrthogroup(app.orthogroups[0].id);
    });
    const label = () => page.evaluate(() => {
      const app = window.__GBDRAW_APP__;
      return app.orthogroupCopyFeedbackLabel(app.selectedOrthogroup, 'aa', 'Copy aa');
    });
    const copy = () => page.evaluate(async () => {
      const app = window.__GBDRAW_APP__;
      await app.copyOrthogroupSequences(app.selectedOrthogroup, 'aa');
      window.__REUSE__.copiedAt = performance.now();
      return app.orthogroupCopyFeedbackLabel(app.selectedOrthogroup, 'aa', 'Copy aa');
    });
    const elapsed = () => page.evaluate(() => performance.now() - window.__REUSE__.copiedAt);

    console.log('STEP Pw-10', Date.now());
    // Pw-10: copy while the reused Generate is in flight; its activation clears the label.
    await page.evaluate(() => {
      const app = window.__GBDRAW_APP__;
      const hooks = window.__GBDRAW_TEST_HOOKS__;
      const record = hooks.onSessionLifecycleEvent;
      window.__REUSE__.feedback = {};
      hooks.beforeDiagramGenerationResponse = async () => {
        delete hooks.beforeDiagramGenerationResponse;
        await app.copyOrthogroupSequences(app.selectedOrthogroup, 'aa');
        window.__REUSE__.copiedAt = performance.now();
        window.__REUSE__.feedback.before = app.orthogroupCopyFeedbackLabel(app.selectedOrthogroup, 'aa', 'Copy aa');
      };
      hooks.onSessionLifecycleEvent = (event) => {
        record(event);
        if (event.name === 'artifact.activation-completed' && window.__REUSE__.copiedAt) {
          window.__REUSE__.feedback.atActivation = app.orthogroupCopyFeedbackLabel(app.selectedOrthogroup, 'aa', 'Copy aa');
          window.__REUSE__.feedback.elapsed = performance.now() - window.__REUSE__.copiedAt;
        }
      };
    });
    const reused = await observeGenerate(page);
    expect(reused.names).toContain('catalog.reused');
    const feedback = await page.evaluate(() => window.__REUSE__.feedback);
    expect(feedback.elapsed, 'inconclusive: activation came after the label timer').toBeLessThan(1500);
    expect(feedback).toMatchObject({ before: 'Copied!', atActivation: 'Copy aa' });

    // An intent step (a group rename) crosses no Generate: Undo keeps the label.
    const groupId = await page.evaluate(() => window.__GBDRAW_APP__.selectedOrthogroup.id);
    await undoable(page, 'Rename the group', { rename: { id: groupId, name: 'Renamed group' } });
    expect(await copy()).toBe('Copied!');
    await history(page, 'undo');
    expect(await elapsed(), 'inconclusive: Undo came after the label timer').toBeLessThan(1500);
    expect(await label()).toBe('Copied!');

    // Undo of a reused Generate step replaces the shown Results: the label clears.
    await undoable(page, 'Set the plot title', { fields: { 'form.plot_title': 'Reuse title' } });
    const retitled = await observeGenerate(page);
    expect(retitled.transports).toEqual([KEPT]);
    expect(retitled.undo).toBe(1);
    expect(await copy()).toBe('Copied!');
    await history(page, 'undo');
    expect(await elapsed(), 'inconclusive: Undo came after the label timer').toBeLessThan(1500);
    expect(await label()).toBe('Copy aa');
    expect(page.externalRequests).toEqual([]);
  } finally {
    await page.context().close();
  }
});

// What one Generate did with the SVG pieces: the segments the request named,
// the children texts the page sanitized and reused, and its mount, bind and
// ready receipt.
const piecesSince = (page, start) => page.evaluate((start) => {
  const metrics = {};
  for (const { name, value } of window.__REUSE__.metrics.slice(start.metrics)) {
    metrics[name] = (metrics[name] || 0) + Number(value || 0);
  }
  const transports = window.__REUSE__.events.slice(start.events).filter((event) => event.name === 'result-transport-ready');
  return {
    held: transports.map((event) => event.heldSegmentCount),
    sanitized: metrics.svgSanitizedPieceCount || 0,
    reused: metrics.svgMemoHitCount || 0,
    decoded: metrics.resultBinaryDecodeCount || 0,
    serialized: metrics.svgSerializationCount || 0,
    mounts: metrics.previewMountAdoptionCount || 0,
    binds: metrics.previewBinderInvocationCount || 0,
    receipts: metrics.previewReadyReceiptAcceptedCount || 0,
    content: window.__GBDRAW_APP__.results[0].content
  };
}, start);
const generatePieces = async (page) => {
  const before = await mark(page);
  await generate(page);
  return piecesSince(page, before);
};
const shownOnce = { mounts: 1, binds: 1, receipts: 1 };
const counts = ({ mounts, binds, receipts }) => ({ mounts, binds, receipts });
const contentPieceCount = (page) => page.evaluate(() =>
  window.__PIECES__.at(-1)[0].document.pieces.filter((piece) => piece.kind === 'content').length);

// The last reply's document joined from its pieces, sanitized whole by DOMPurify
// with and without the XML declaration the file has; and the pieces sanitized
// one by one as admission does.
const wholeAndPieces = (page) => page.evaluate(async () => {
  const { SVG_SANITIZE_OPTIONS, wrapSvgFragment, unwrapSvgFragment } = await import('./js/services/svg-sanitization.js');
  const [{ document, children }] = window.__PIECES__.at(-1);
  const raw = document.head + document.pieces.map((piece, index) => (piece.kind === 'raw'
    ? piece.text : `${piece.start}${children[index]}</g>`)).join('') + document.tail;
  const sanitize = (text) => window.DOMPurify.sanitize(text, SVG_SANITIZE_OPTIONS);
  const prefix = sanitize(wrapSvgFragment('')).slice(0, -'</svg>'.length);
  const fragment = (text) => unwrapSvgFragment(sanitize(wrapSvgFragment(text)), prefix);
  const pieces = sanitize(`${document.head}</svg>`).slice(0, -'</svg>'.length) + document.pieces.map((piece, index) => (
    piece.kind === 'raw' ? fragment(piece.text) : `${fragment(`${piece.start}</g>`).slice(0, -'</g>'.length)}${fragment(children[index])}</g>`
  )).join('') + '</svg>';
  return {
    whole: sanitize(raw),
    withDeclaration: sanitize(`<?xml version="1.0" encoding="utf-8" ?>\n${raw}`),
    pieces
  };
});

const undoCount = (page) => page.evaluate(() => window.__GBDRAW_HISTORY__.getUndoCount());
const setTitle = (page, title) => undoable(page, 'Set the plot title', { fields: { 'form.plot_title': title } });
// A live palette edit rewrites the mounted Result (hazard 6); the draft keeps the color.
const livePalette = (page) => evaluateWithRetainedPromise(page, async () => {
  const { state } = await import('./js/state.js');
  state.paletteInstantPreviewEnabled.value = true;
  await window.__GBDRAW_HISTORY__.runUndoable('Live palette', async () => {
    state.activeDrawing().currentColors.value = { ...state.activeDrawing().currentColors.value, CDS: '#123456' };
    await window.Vue.nextTick();
  });
});

// A Generate reply whose pieces carry hostile markup, admitted as the current
// Worker's: every byte still passes DOMPurify.
const admitHostileReply = (page) => evaluateWithRetainedPromise(page, async () => {
  const { state } = await import('./js/state.js');
  const { normalizeGenerationResponse } = await import('./js/services/diagram-generation.js');
  const { admitFeatureCatalog } = await import('./js/services/feature-catalog.js');
  const { admitCurrentGeneratedResults, createEmptySvgMutationPlan } = await import('./js/services/svg-result-ingestion.js');
  const featureCatalog = JSON.parse(JSON.stringify(window.Vue.toRaw(state.featureCatalog.value)));
  const children = [
    '<path d="M0,0" onload="alert(1)"/><script>alert(2)</script>',
    '<a href="javascript:alert(3)"><text>x</text></a><g onclick="alert(4)"/>',
    '</g><script>alert(5)</script><g><image href="x" onerror="alert(6)"/>'
  ];
  const encoder = new TextEncoder();
  const bytes = children.map((text) => encoder.encode(text));
  const segments = new Uint8Array(bytes.reduce((total, item) => total + item.byteLength, 0));
  let offset = 0;
  const pieces = [{ kind: 'raw', text: '<foreignObject><div onclick="alert(7)"></div></foreignObject>' }];
  bytes.forEach((item, index) => {
    segments.set(item, offset);
    pieces.push({
      kind: 'content', start: `<g id="hostile${index}" onmouseover="alert(8)">`, revision: String(index).repeat(64),
      bytes: item.byteLength, offset, length: item.byteLength
    });
    offset += item.byteLength;
  });
  const response = normalizeGenerationResponse({
    results: [{ name: featureCatalog.items[0].resultName, document: { head: '<svg xmlns="http://www.w3.org/2000/svg" onload="alert(9)">', tail: '</svg>', pieces }, segments }],
    metadata: {},
    featureCatalog
  });
  const catalogAdmission = admitFeatureCatalog(featureCatalog, response.results, { adopt: true, mode: window.__GBDRAW_APP__.mode });
  const [committed] = admitCurrentGeneratedResults(response, { catalogAdmission, mutationPlan: createEmptySvgMutationPlan(1) });
  return committed.content;
});

test('a Generate sanitizes each changed SVG group once and admits what a fresh page admits', async ({ browser }) => {
  test.setTimeout(600_000);
  const page = await start(browser, { captureReplies: true });
  const fresh = await start(browser, { captureReplies: true });
  try {
    console.log('STEP S-1', Date.now());
    // S-1 (Session Load): Results admitted from their saved text hold nothing,
    // so the first Generate sends and sanitizes every children text once. The
    // joined sanitized pieces equal DOMPurify of the whole document, with or
    // without its XML declaration.
    await loadSession(page, `${SESSIONS}/tobacco-chloroplast.gbdraw-session.json`);
    const cold = await generatePieces(page);
    const pieceCount = await contentPieceCount(page);
    expect(pieceCount).toBeGreaterThan(1);
    expect({ held: cold.held, sanitized: cold.sanitized, reused: cold.reused, decoded: cold.decoded, serialized: cold.serialized })
      .toEqual({ held: [0], sanitized: pieceCount, reused: 0, decoded: 1, serialized: 0 });
    expect(counts(cold)).toEqual(shownOnce);
    const sanitizedWhole = await wholeAndPieces(page);
    expect(sanitizedWhole.pieces).toBe(sanitizedWhole.whole);
    expect(sanitizedWhole.withDeclaration).toBe(sanitizedWhole.whole);
    expect(cold.content).toBe(sanitizedWhole.whole);

    console.log('STEP S-2', Date.now());
    // S-2: the same request holds every revision: nothing is sent or sanitized
    // again, and the Result is the same text, mounted and bound once.
    const warm = await generatePieces(page);
    expect({ held: warm.held, sanitized: warm.sanitized, reused: warm.reused, decoded: warm.decoded })
      .toEqual({ held: [pieceCount], sanitized: 0, reused: pieceCount, decoded: 0 });
    expect(counts(warm)).toEqual(shownOnce);
    expect(warm.content).toBe(cold.content);

    console.log('STEP S-3', Date.now());
    // S-3 (live-edit residue): a live palette edit and a title change; the
    // Generate sanitizes only the groups whose children changed, and its Result
    // equals a fresh page's that drew the same draft from nothing.
    await livePalette(page);
    await setTitle(page, 'Pieces title');
    const edited = await generatePieces(page);
    expect(edited.sanitized).toBeGreaterThan(0);
    expect(edited.reused).toBeGreaterThan(0);
    expect(edited.sanitized + edited.reused).toBe(await contentPieceCount(page));
    expect(counts(edited)).toEqual(shownOnce);
    await loadSession(fresh, `${SESSIONS}/tobacco-chloroplast.gbdraw-session.json`);
    await livePalette(fresh);
    await setTitle(fresh, 'Pieces title');
    const freshEdited = await generatePieces(fresh);
    expect(freshEdited.reused).toBe(0);
    expect(edited.content).toBe(freshEdited.content);

    console.log('STEP S-4', Date.now());
    // S-4 (Worker restart): revisions are content hashes, so a new Worker
    // omits the same texts and the Result is unchanged.
    await page.evaluate(async () => {
      (await import('./js/services/diagram-generation.js')).disposeDiagramGenerationWorker();
    });
    const restarted = await generatePieces(page);
    expect({ sanitized: restarted.sanitized, reused: restarted.reused }).toEqual({ sanitized: 0, reused: await contentPieceCount(page) });
    expect(restarted.content).toBe(edited.content);

    console.log('STEP S-5', Date.now());
    // S-5 (Undo, Redo; Generate A, Generate B, Undo, Generate): restored
    // Results keep their segments, so a Generate after Redo or after Undo
    // sanitizes nothing and admits the restored Result's text again.
    const restoredContent = () => page.evaluate(() => window.__GBDRAW_APP__.results[0].content);
    const allHeld = async () => ({ sanitized: 0, reused: await contentPieceCount(page) });
    const beforeB = await undoCount(page);
    await undoable(page, 'Show the plot title', { fields: { 'adv.plot_title_position': 'top', 'form.plot_title': 'Undo title' } });
    const second = await generatePieces(page);
    expect(second.content).not.toBe(restarted.content);
    const steps = (await undoCount(page)) - beforeB;
    for (let step = 0; step < steps; step += 1) await history(page, 'undo');
    expect(await restoredContent()).toBe(restarted.content);
    for (let step = 0; step < steps; step += 1) await history(page, 'redo');
    expect(await restoredContent()).toBe(second.content);
    const redone = await generatePieces(page);
    expect({ sanitized: redone.sanitized, reused: redone.reused }).toEqual(await allHeld());
    expect(redone.content).toBe(second.content);
    const back = (await undoCount(page)) - beforeB;
    for (let step = 0; step < back; step += 1) await history(page, 'undo');
    expect(await restoredContent()).toBe(restarted.content);
    const undone = await generatePieces(page);
    expect({ sanitized: undone.sanitized, reused: undone.reused }).toEqual(await allHeld());
    expect(undone.content).toBe(restarted.content);

    console.log('STEP S-6', Date.now());
    // S-6 (cancel after the reply, before the commit): nothing is committed;
    // the next request names the same segments.
    await setTitle(page, 'Cancelled title');
    await page.evaluate(() => {
      window.__GBDRAW_TEST_HOOKS__.beforeDiagramGenerationResponse = async () => {
        delete window.__GBDRAW_TEST_HOOKS__.beforeDiagramGenerationResponse;
        (await import('./js/services/diagram-generation.js')).cancelDiagramGeneration();
      };
    });
    const cancelMark = await mark(page);
    await page.getByRole('button', { name: 'Generate Diagram', exact: true }).click();
    await expect.poll(() => page.evaluate(async () => (await import('./js/state.js')).state.processing.value), { timeout: 180_000 }).toBe(false);
    const cancelled = await piecesSince(page, cancelMark);
    expect(cancelled.held).toEqual([await contentPieceCount(page)]);
    expect(cancelled.content).toBe(restarted.content);
    await page.evaluate(async () => { (await import('./js/state.js')).state.errorLog.value = null; });
    const afterCancel = await generatePieces(page);
    expect(afterCancel.held).toEqual(cancelled.held);
    expect(afterCancel.reused).toBeGreaterThan(0);

    console.log('STEP S-7', Date.now());
    // S-7: hostile markup in the head, a start tag, a raw piece or a children
    // text (also one that closes its group early) leaves no script, event
    // handler or javascript: link.
    const hostile = await admitHostileReply(page);
    expect(hostile.startsWith('<svg')).toBe(true);
    expect(hostile).not.toMatch(/<script|\son[a-z]+=|javascript:|<foreignObject|<image/i);
    expect(page.externalRequests).toEqual([]);
    expect(fresh.externalRequests).toEqual([]);
  } finally {
    await page.context().close();
    await fresh.context().close();
  }
});
