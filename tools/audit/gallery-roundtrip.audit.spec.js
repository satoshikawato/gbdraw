// Gallery session round trip: load every Gallery session, Save it, load the saved file in a
// fresh browser context, and Save again. User-owned state (the G-C snapshot) and the saved
// documents must match. Restrict the set with AUDIT_SESSIONS=name1,name2.
const { readdirSync } = require('node:fs');
const { join } = require('node:path');
const { test, expect } = require('@playwright/test');
const A = require('./helpers/audit-common.cjs');

const SESSION_DIR = join(A.REPO_ROOT, 'gbdraw/web/gallery/sessions');
const names = process.env.AUDIT_SESSIONS
  ? process.env.AUDIT_SESSIONS.split(',')
  : readdirSync(SESSION_DIR).filter((name) => /\.gbdraw-session\.json(\.gz)?$/.test(name)).sort();

// Volatile document paths that differ between two saves of the same state.
const VOLATILE = /savedAt|createdAt|exportedAt/;

for (const name of names) {
  test(`round trip ${name}`, async ({ page, browser, baseURL }, info) => {
    test.setTimeout(900_000);
    const outdir = A.outDir('gallery-roundtrip');
    const { openApp, snapshotUserOwnedState, diffUserOwnedState } = A.helpers();
    await openApp(page);
    const first = await A.loadSession(page, join(SESSION_DIR, name));
    const before = await snapshotUserOwnedState(page);
    const save1 = await A.saveSession(page, info.outputPath('s1.gbdraw-session.json.gz'));

    const context = await browser.newContext({ baseURL, acceptDownloads: true });
    const fresh = await context.newPage();
    const log = A.collect(fresh);
    await openApp(fresh);
    const second = await A.loadSession(fresh, info.outputPath('s1.gbdraw-session.json.gz'));
    const after = await snapshotUserOwnedState(fresh);
    const save2 = await A.saveSession(fresh, info.outputPath('s2.gbdraw-session.json.gz'));
    await context.close();

    const stateDiff = diffUserOwnedState(before, after).filter((d) => !String(d.path).startsWith('history'));
    const docDiff = diffUserOwnedState(save1.doc, save2.doc).filter((d) => !VOLATILE.test(String(d.path)));
    A.writeEvidence(outdir, `${name}.json`, {
      name,
      loadMs: [first.ms, second.ms],
      loadErrors: [first.error, second.error],
      stateDiff,
      docDiffCount: docDiff.length,
      docDiff: docDiff.slice(0, 200),
      log
    });
    expect(first.error, JSON.stringify(first.error)).toBeNull();
    expect(second.error, JSON.stringify(second.error)).toBeNull();
    expect(stateDiff, JSON.stringify(stateDiff).slice(0, 4000)).toEqual([]);
    expect(docDiff, JSON.stringify(docDiff).slice(0, 4000)).toEqual([]);
  });
}
