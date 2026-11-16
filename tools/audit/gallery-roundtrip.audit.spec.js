// Gallery session round trip: load every Gallery session, Save it, load the saved file in a
// fresh browser context, and Save again. User-owned state (the G-C snapshot), the slices of
// every drawing (`sliceHashes`, the random walk's reader), and the saved documents must match,
// and each Load must leave no UNKNOWN or VALIDATION_UNCLASSIFIED alert, page or console error,
// or NaN or undefined in Run Info. Restrict the set with AUDIT_SESSIONS=name1,name2.
const { readdirSync } = require('node:fs');
const { join } = require('node:path');
const { test, expect } = require('@playwright/test');
const A = require('./helpers/audit-common.cjs');
const { restoredChanges, sliceHashes } = require('../../tests/web/helpers/drawing-slice-hashes.mjs');

const SESSION_DIR = join(A.REPO_ROOT, 'gbdraw/web/gallery/sessions');
const names = process.env.AUDIT_SESSIONS
  ? process.env.AUDIT_SESSIONS.split(',')
  : readdirSync(SESSION_DIR).filter((name) => name.endsWith('.gbdraw')).sort();

// Volatile document paths that differ between two saves of the same state.
const VOLATILE = /savedAt|createdAt|exportedAt/;

for (const name of names) {
  test(`round trip ${name}`, async ({ page, browser, baseURL }, info) => {
    test.setTimeout(900_000);
    const outdir = A.outDir('gallery-roundtrip');
    const { openApp, readErrorSignature, snapshotUserOwnedState, diffUserOwnedState } = A.helpers();
    const loadWithHealth = async (target, file, operation) => {
      const signature = await readErrorSignature(target);
      const load = await A.loadSession(target, file);
      return { ...load, health: await A.operationHealth(target, operation, signature) };
    };
    await openApp(page);
    // Starts collecting the page's errors; openApp has checked those of its own open.
    await A.operationHealth(page, 'open');
    const first = await loadWithHealth(page, join(SESSION_DIR, name), 'Load of the Gallery Session');
    const before = await snapshotUserOwnedState(page);
    const slicesBefore = await sliceHashes(page);
    const save1 = await A.saveSession(page, info.outputPath('s1.gbdraw'));

    const context = await browser.newContext({ baseURL, acceptDownloads: true });
    const fresh = await context.newPage();
    const log = A.collect(fresh);
    await openApp(fresh);
    await A.operationHealth(fresh, 'open');
    const second = await loadWithHealth(fresh, info.outputPath('s1.gbdraw'), 'Load of the saved Session');
    const after = await snapshotUserOwnedState(fresh);
    const sliceDiff = restoredChanges('Load, Save, Load', slicesBefore, await sliceHashes(fresh), { by: 'renamed' });
    const save2 = await A.saveSession(fresh, info.outputPath('s2.gbdraw'));
    await context.close();

    const stateDiff = diffUserOwnedState(before, after).filter((d) => !String(d.path).startsWith('history'));
    const docDiff = diffUserOwnedState(save1.doc, save2.doc).filter((d) => !VOLATILE.test(String(d.path)));
    A.writeEvidence(outdir, `${name}.json`, {
      name,
      loadMs: [first.ms, second.ms],
      loadErrors: [first.error, second.error],
      health: [...first.health, ...second.health],
      stateDiff,
      sliceDiff,
      docDiffCount: docDiff.length,
      docDiff: docDiff.slice(0, 200),
      log
    });
    expect(first.error, JSON.stringify(first.error)).toBeNull();
    expect(second.error, JSON.stringify(second.error)).toBeNull();
    expect([...first.health, ...second.health]).toEqual([]);
    expect(stateDiff, JSON.stringify(stateDiff).slice(0, 4000)).toEqual([]);
    expect(sliceDiff, 'the slices of every drawing after Load, Save, Load').toEqual([]);
    expect(docDiff, JSON.stringify(docDiff).slice(0, 4000)).toEqual([]);
  });
}
