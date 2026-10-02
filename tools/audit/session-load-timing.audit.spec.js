// Session load timing: load Gallery sessions and record when each Worker starts, each Worker
// message is posted and answered, and each session lifecycle event fires, relative to the load.
// Choose sessions with AUDIT_SESSIONS=name1,name2 (default: the tobacco plastome and HmmtDNA).
const { join } = require('node:path');
const { test } = require('@playwright/test');
const A = require('./helpers/audit-common.cjs');

const SESSION_DIR = join(A.REPO_ROOT, 'gbdraw/web/gallery/sessions');
const names = (process.env.AUDIT_SESSIONS
  || 'tobacco-chloroplast.gbdraw-session.json,HmmtDNA_basic_circular.gbdraw-session.json').split(',');

const INIT = () => {
  window.__AUDIT_EVENTS__ = [];
  const note = (name) => window.__AUDIT_EVENTS__.push([Math.round(performance.now()), name]);
  window.__GBDRAW_TEST_HOOKS__ = {
    ...(window.__GBDRAW_TEST_HOOKS__ || {}),
    onSessionLifecycleEvent: (event) => note(event.name)
  };
  const NativeWorker = window.Worker;
  window.Worker = new Proxy(NativeWorker, {
    construct(target, args) {
      const worker = Reflect.construct(target, args, target);
      note(`WORKER ${String(args[0]).split('/').pop()}`);
      const post = worker.postMessage.bind(worker);
      worker.postMessage = (message, transfer) => {
        if (message && message.type) note(`post ${message.type}:${message.operation || ''}`);
        return transfer === undefined ? post(message) : post(message, transfer);
      };
      worker.addEventListener('message', (event) => {
        const message = event.data || {};
        if (message.type && message.status !== 'part') note(`reply ${message.type}:${message.ok}`);
      });
      return worker;
    }
  });
};

for (const name of names) {
  test(`load timing ${name}`, async ({ page }) => {
    test.setTimeout(600_000);
    const outdir = A.outDir('session-load-timing');
    await page.addInitScript(INIT);
    await A.helpers().openApp(page);
    const start = await page.evaluate(() => performance.now());
    const load = await A.loadSession(page, join(SESSION_DIR, name));
    const events = await page.evaluate((t0) => window.__AUDIT_EVENTS__
      .filter(([at]) => at >= t0)
      .map(([at, label]) => [Math.round(at - t0), label]), start);
    A.writeEvidence(outdir, `${name}.json`, { name, ms: load.ms, error: load.error, events });
    console.log(name, `${load.ms} ms`, `${events.length} events`);
  });
}
