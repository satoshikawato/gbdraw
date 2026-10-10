// Steps, oracles, and evidence for the promotion user journeys
// (tools/audit/user-journey.audit.spec.js).
//
// Every step runs within its own time limit and then checks the common oracles:
// no page error, console error, or unhandled rejection on a watched page, and no
// busy indicator left after the step settles. Each step saves a screenshot and
// appends to `steps.json` under $GBDRAW_AUDIT_OUT/journeys/<journey>/;
// `writeContactSheet` builds journeys/index.html, failures first.
//
// CLI: node tools/audit/helpers/journey-evidence.cjs <journeys dir>
// rebuilds the contact sheet, for example after merging the CI shards.
const { existsSync, readFileSync, readdirSync, writeFileSync } = require('node:fs');
const { join } = require('node:path');

const BUSY_LIMIT_MS = 60_000;

// True when no Generate, reflow, Session load or save, or rule match is running,
// and no visible element is marked busy or spinning.
const idleInPage = () => {
  const app = window.__GBDRAW_APP__;
  if (app && (app.processing || app.labelReflowProcessing || app.sessionImportPending
    || app.sessionSavePending || app.ruleMatchingPending)) return false;
  const visible = (element) => element.getClientRects().length > 0 && getComputedStyle(element).visibility !== 'hidden';
  return ![...document.querySelectorAll('[aria-busy="true"], .animate-spin')].some(visible);
};

const busyInPage = () => {
  const app = window.__GBDRAW_APP__;
  const flags = app
    ? ['processing', 'labelReflowProcessing', 'sessionImportPending', 'sessionSavePending', 'ruleMatchingPending']
      .filter((name) => app[name])
    : [];
  const visible = (element) => element.getClientRects().length > 0 && getComputedStyle(element).visibility !== 'hidden';
  const elements = [...document.querySelectorAll('[aria-busy="true"], .animate-spin')].filter(visible)
    .map((element) => `${element.tagName.toLowerCase()}${element.id ? `#${element.id}` : ''}`
      + ` "${(element.getAttribute('aria-label') || element.textContent || '').trim().slice(0, 40)}"`);
  return { flags, elements };
};

const withinLimit = (promise, limitMs, what) => {
  let timer;
  const expired = new Promise((_, reject) => {
    timer = setTimeout(() => reject(new Error(`${what} did not finish within ${Math.round(limitMs / 1000)} s`)), limitMs);
  });
  return Promise.race([promise, expired]).finally(() => clearTimeout(timer));
};

const safeName = (text) => String(text).replace(/[^A-Za-z0-9._-]+/g, '-').replace(/^-+|-+$/g, '').slice(0, 60);

const createJourney = ({ id, title, dir }) => {
  const record = { id, title, startedAt: new Date().toISOString(), durationMs: 0, result: 'running', steps: [], dialogs: [] };
  const problems = [];
  const softFailures = [];
  const started = Date.now();
  let screenshotPage = null;

  const save = () => {
    record.durationMs = Date.now() - started;
    writeFileSync(join(dir, 'steps.json'), JSON.stringify(record, null, 2));
  };

  // Watches a page for the common oracles and accepts its dialogs, recording each.
  const watch = async (page, label = 'app') => {
    page.on('pageerror', (error) => problems.push({ page: label, kind: 'pageerror', text: String(error?.stack || error?.message || error) }));
    page.on('console', (message) => {
      if (message.type() === 'error') problems.push({ page: label, kind: 'console.error', text: message.text() });
    });
    // A helper that answers a dialog itself (audit-common's Session load and
    // save) answers first; this handler then finds it handled.
    page.on('dialog', (dialog) => {
      record.dialogs.push({ page: label, type: dialog.type(), message: dialog.message() });
      setTimeout(() => { dialog.accept().catch(() => {}); }, 0);
    });
    await page.addInitScript(() => {
      window.addEventListener('unhandledrejection', (event) => {
        // eslint-disable-next-line no-console
        console.error(`unhandled rejection: ${event.reason?.stack || event.reason}`);
      });
    });
    if (!screenshotPage) screenshotPage = page;
    return page;
  };

  // Shows this page in the screenshots from now on.
  const show = (page) => { screenshotPage = page; };

  const step = async (name, run, { limitMs = 120_000, soft = false } = {}) => {
    const index = record.steps.length + 1;
    const entry = { step: index, name, durationMs: 0, result: 'running' };
    record.steps.push(entry);
    const problemsBefore = problems.length;
    const stepStarted = Date.now();
    let error = null;
    try {
      entry.detail = await withinLimit(Promise.resolve().then(run), limitMs, `step "${name}"`);
      const page = screenshotPage;
      if (page && !page.isClosed()) {
        await page.waitForFunction(idleInPage, null, { timeout: BUSY_LIMIT_MS }).catch(async () => {
          const busy = await page.evaluate(busyInPage).catch(() => ({ flags: ['(page unreadable)'], elements: [] }));
          throw new Error(`a busy indicator stayed ${BUSY_LIMIT_MS / 1000} s after the step: ${JSON.stringify(busy)}`);
        });
      }
      const fresh = problems.slice(problemsBefore);
      if (fresh.length) {
        throw new Error(`the step raised ${fresh.length} page error(s), console error(s), or unhandled rejection(s):\n`
          + fresh.map((problem) => `  [${problem.page}] ${problem.kind}: ${problem.text.slice(0, 600)}`).join('\n'));
      }
      entry.result = 'passed';
    } catch (caught) {
      error = caught;
      entry.result = 'failed';
      entry.error = String(caught?.message || caught).slice(0, 4000);
    } finally {
      entry.durationMs = Date.now() - stepStarted;
      const page = screenshotPage;
      if (page && !page.isClosed()) {
        entry.screenshot = `${String(index).padStart(2, '0')}-${safeName(name)}.png`;
        await page.screenshot({ path: join(dir, entry.screenshot), timeout: 30_000 }).catch(() => { delete entry.screenshot; });
      }
      save();
    }
    if (error) {
      if (!soft) {
        record.result = 'failed';
        save();
        throw error;
      }
      softFailures.push(`step ${index} "${name}": ${entry.error}`);
    }
    return entry.detail;
  };

  // Ends the journey; soft step failures fail it here, after every step ran.
  const finish = () => {
    record.result = softFailures.length ? 'failed' : 'passed';
    save();
    if (softFailures.length) throw new Error(`${softFailures.length} step(s) failed:\n${softFailures.join('\n')}`);
  };

  // A test that ends early (a failed step or the test deadline) still leaves its record.
  const close = () => {
    if (record.result === 'running') record.result = 'failed';
    save();
  };

  save();
  return { watch, show, step, finish, close, record };
};

const escapeHtml = (text) => String(text ?? '').replace(/[&<>"']/g, (character) => ({
  '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#39;'
}[character]));

// journeys/index.html: one section per journey, failed journeys first, each with
// its step screenshots in order and the failed steps' errors.
const writeContactSheet = (journeysDir) => {
  const records = readdirSync(journeysDir, { withFileTypes: true })
    .filter((entry) => entry.isDirectory() && existsSync(join(journeysDir, entry.name, 'steps.json')))
    .map((entry) => ({ folder: entry.name, ...JSON.parse(readFileSync(join(journeysDir, entry.name, 'steps.json'), 'utf8')) }))
    .sort((left, right) => Number(left.result === 'passed') - Number(right.result === 'passed')
      || String(left.id).localeCompare(String(right.id), undefined, { numeric: true }));
  const seconds = (ms) => `${(Number(ms || 0) / 1000).toFixed(1)} s`;
  const sections = records.map((journey) => {
    const steps = (journey.steps || []).map((item) => `
      <figure class="${item.result === 'passed' ? 'ok' : 'bad'}">
        ${item.screenshot ? `<a href="${escapeHtml(`${journey.folder}/${item.screenshot}`)}"><img loading="lazy" src="${escapeHtml(`${journey.folder}/${item.screenshot}`)}" alt="${escapeHtml(item.name)}"></a>` : '<div class="none">no screenshot</div>'}
        <figcaption>${item.step}. ${escapeHtml(item.name)} (${seconds(item.durationMs)}, ${escapeHtml(item.result)})
          ${item.error ? `<pre>${escapeHtml(item.error)}</pre>` : ''}</figcaption>
      </figure>`).join('');
    return `<section><h2 class="${journey.result === 'passed' ? 'ok' : 'bad'}">${escapeHtml(journey.id)} ${escapeHtml(journey.title)}: `
      + `${escapeHtml(journey.result)}, ${seconds(journey.durationMs)}</h2><div class="grid">${steps}</div></section>`;
  }).join('\n');
  const total = records.reduce((sum, journey) => sum + Number(journey.durationMs || 0), 0);
  const failed = records.filter((journey) => journey.result !== 'passed').length;
  writeFileSync(join(journeysDir, 'index.html'), `<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">
<title>User journeys</title>
<style>
body { font: 14px/1.4 system-ui, sans-serif; margin: 16px; background: #fff; color: #111; }
.grid { display: grid; grid-template-columns: repeat(auto-fill, minmax(260px, 1fr)); gap: 12px; }
figure { margin: 0; border: 2px solid #cbd5e1; border-radius: 6px; padding: 6px; }
figure.bad { border-color: #dc2626; }
img { width: 100%; height: auto; display: block; }
h2.bad { color: #b91c1c; } h2.ok { color: #166534; }
pre { white-space: pre-wrap; font-size: 12px; color: #b91c1c; max-height: 240px; overflow: auto; }
.none { color: #64748b; padding: 24px 0; text-align: center; }
</style></head><body>
<h1>User journeys</h1>
<p>${records.length} journeys, ${failed} failed, ${seconds(total)} in total.</p>
${sections}
</body></html>
`);
  return { journeys: records.length, failed, durationMs: total };
};

if (require.main === module) {
  const dir = process.argv[2];
  if (!dir) {
    console.error('usage: node tools/audit/helpers/journey-evidence.cjs <journeys dir>');
    process.exit(2);
  }
  console.log(JSON.stringify(writeContactSheet(dir)));
}

module.exports = { createJourney, writeContactSheet };
