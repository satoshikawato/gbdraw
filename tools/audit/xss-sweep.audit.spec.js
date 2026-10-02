// Crafted-qualifier sweep: script-like text in DEFINITION, ORGANISM, qualifiers, and the file name
// must stay inert in the preview, search, feature popup, Run info, and both SVG exports.
// The sanitizer itself is covered by tests/web/svg-sanitization.test.mjs; this sweep checks the
// surfaces end to end. Set XSS_FILE to a GenBank path to use your own payload record.
const fs = require('node:fs');
const { join } = require('node:path');
const { pathToFileURL } = require('node:url');
const { test, expect } = require('@playwright/test');
const A = require('./helpers/audit-common.cjs');

const P = (tag) => `<img src=x onerror=alert('XSS-${tag}')>`;
const sequence = (length) => {
  const bases = 'acgtacgtaa'.repeat(Math.ceil(length / 10)).slice(0, length);
  const lines = [];
  for (let i = 0; i < length; i += 60) {
    lines.push(`${String(i + 1).padStart(9)} ${bases.slice(i, i + 60).match(/.{1,10}/g).join(' ')}`);
  }
  return lines.join('\n');
};
const PAYLOAD_RECORD = `LOCUS       XSSREC                  3000 bp    DNA     circular UNA 01-JAN-2000
DEFINITION  Evil ${P('def')} </text><script>alert('XSS-defscript')</script> def.
ACCESSION   XSSREC
VERSION     XSSREC.1
KEYWORDS    .
SOURCE      synthetic construct
  ORGANISM  <i>Evil</i> ${P('org')}
            .
FEATURES             Location/Qualifiers
     source          1..3000
                     /organism="Evil ${P('orgq')}"
                     /strain="<svg onload=alert('XSS-strain')>"
     CDS             100..900
                     /gene="</text><script>alert('XSS-gene')</script>"
                     /locus_tag="LT_${P('lt')}"
                     /product="${P('product')}"
                     /note="javascript:alert('XSS-note')"
                     /db_xref="<a href=javascript:alert('XSS-dbx')>x</a>"
                     /translation="MKKKKKKKKK"
     CDS             complement(1200..2000)
                     /gene="<foreignObject><iframe srcdoc='<script>alert(1)</script>'></iframe></foreignObject>"
                     /product="]]><script>alert('XSS-cdata')</script>"
                     /translation="MKKKKKKKKK"
     tRNA            2200..2280
                     /product="tRNA-\\"><img src=x onerror=alert('XSS-trna')>"
ORIGIN
${sequence(3000)}
//
`;
const FILE_NAME = process.env.XSS_FNAME || `evil<img src=x onerror=alert('XSS-fname')>.gbk`;

test('crafted qualifiers stay inert in preview, popup, search, Run info, and exports', async ({ page, context }) => {
  test.setTimeout(420_000);
  page.setDefaultTimeout(30_000);
  const outdir = A.outDir('xss');
  const payload = process.env.XSS_FILE ? fs.readFileSync(process.env.XSS_FILE) : Buffer.from(PAYLOAD_RECORD);
  const dialogs = [];
  const report = { file: process.env.XSS_FILE || 'synthetic', steps: [] };
  const save = () => A.writeEvidence(outdir, 'xss-report.json', report);
  const step = async (name, fn) => {
    try { report[name] = await fn(); report.steps.push(`${name}:ok`); } catch (e) {
      report.steps.push(`${name}:ERR ${String(e.message || e).slice(0, 300)}`);
    }
    save();
  };
  page.on('dialog', (d) => { dialogs.push({ page: 'app', msg: d.message() }); d.dismiss(); });
  const { openApp, generateAndWaitForResult } = A.helpers();
  await openApp(page);
  await page.getByLabel('GenBank/DDBJ File', { exact: true }).setInputFiles({
    name: FILE_NAME, mimeType: 'text/plain', buffer: payload
  });
  await page.evaluate(() => { window.__GBDRAW_APP__.adv.rich_feature_popup = true; });
  const gen = await generateAndWaitForResult(page, { expectedStatus: null, requireCommittedResult: false });
  report.generate = gen.result?.status;
  report.generateError = gen.errorSummary;
  save();
  if (gen.result?.status === 'ok') {
    await step('previewDangerous', () => page.evaluate(() => {
      const svg = document.querySelector('.origin-top svg');
      return {
        img: svg.querySelectorAll('img').length,
        script: svg.querySelectorAll('script').length,
        foreignObject: svg.querySelectorAll('foreignObject').length,
        iframe: document.querySelectorAll('iframe').length,
        onattrs: [...svg.querySelectorAll('*')].flatMap((el) => [...el.attributes]
          .filter((a) => /^on/i.test(a.name)).map((a) => a.name))
      };
    }));
    await step('search', async () => {
      const search = page.getByRole('searchbox', { name: 'Search features', exact: true });
      await search.fill('img');
      await search.press('Enter');
      await page.waitForTimeout(1000);
      const open = page.getByRole('button', { name: 'Open active feature', exact: true });
      let popupText = null;
      if (await open.isVisible().catch(() => false)) {
        await open.click();
        await page.waitForTimeout(1000);
        popupText = await page.locator('.feature-popup').innerText({ timeout: 5000 }).catch((e) => String(e));
        await page.getByRole('button', { name: 'Close feature popup', exact: true }).click({ timeout: 5000 }).catch(() => {});
      }
      return { popupText, imgs: await page.evaluate(() => document.querySelectorAll('img[src="x"]').length) };
    });
    await step('runInfo', async () => {
      const tabs = page.locator('.result-tabs');
      await tabs.getByRole('button', { name: /Run info/i }).click({ timeout: 5000 });
      await page.waitForTimeout(500);
      const imgs = await page.evaluate(() => document.querySelectorAll('img[src="x"]').length);
      await tabs.getByRole('button', { name: /Preview/i }).click({ timeout: 5000 });
      return { imgs };
    });
    const exportOne = async (method, outName) => {
      const pending = page.waitForEvent('download', { timeout: 60_000 });
      await page.evaluate((m) => window.__GBDRAW_APP__[m](), method);
      const download = await pending;
      await download.saveAs(join(outdir, outName));
      return { name: download.suggestedFilename() };
    };
    await step('interactiveExport', () => exportOne('downloadInteractiveSVG', 'xss-interactive.svg'));
    await step('staticExport', () => exportOne('downloadSVG', 'xss-static.svg'));
    if (fs.existsSync(join(outdir, 'xss-interactive.svg'))) {
      await step('viewer', async () => {
        const viewer = await context.newPage();
        viewer.setDefaultTimeout(10_000);
        viewer.on('dialog', (d) => { dialogs.push({ page: 'viewer', msg: d.message() }); d.dismiss(); });
        const viewerErrors = [];
        viewer.on('pageerror', (e) => viewerErrors.push(String(e)));
        await viewer.goto(pathToFileURL(join(outdir, 'xss-interactive.svg')).href);
        await viewer.waitForTimeout(1000);
        const features = viewer.locator('[data-gbdraw-interactive-feature="true"]');
        const count = await features.count();
        for (let i = 0; i < Math.min(count, 40); i += 1) {
          await features.nth(i).dispatchEvent('click').catch(() => {});
          await viewer.waitForTimeout(100);
          for (const tab of ['Qualifiers', 'Sequence', 'Details']) {
            const button = viewer.locator('.gfi-tab', { hasText: tab });
            if (await button.count()) await button.first().click().catch(() => {});
          }
          await features.nth(i).hover({ force: true, timeout: 1000 }).catch(() => {});
        }
        const imgs = await viewer.evaluate(() => document.querySelectorAll('img').length);
        await viewer.close();
        return { count, imgs, viewerErrors };
      });
    }
  }
  report.dialogs = dialogs;
  save();
  expect(dialogs).toEqual([]);
});
