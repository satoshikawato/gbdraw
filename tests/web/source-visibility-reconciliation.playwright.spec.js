const { test, expect } = require('@playwright/test');
const fs = require('node:fs/promises');
const { gunzipSync } = require('node:zlib');
const { execFile } = require('node:child_process');
const { promisify } = require('node:util');
const { load, generate, switchMode, popup, closeEditor, download, snapshot } = require('./helpers/mode-transition.cjs');

const seed = 'gbdraw/web/gallery/sessions/HmmtDNA_basic_circular.gbdraw-session.json';
const tobacco = 'gbdraw/web/gallery/sessions/tobacco-chloroplast.gbdraw-session.json';
const lambda = 'gbdraw/web/gallery/sessions/lambda_basic_linear.gbdraw-session.json';
const source = async (file, name) => {
  const session = JSON.parse(await fs.readFile(file, 'utf8'));
  return { name, mimeType: 'text/plain', buffer: Buffer.from(
    session.resources[session.renderRequest.records[0].source.resourceId].data, 'base64') };
};
const upload = async (page, file) => {
  await page.getByLabel('GenBank/DDBJ File', { exact: true }).setInputFiles(file);
  await expect.poll(() => page.evaluate(async () =>
    (await import('./js/state.js')).state.circularRecordDiscovery.status)).not.toBe('loading');
};
const hide = async page => {
  const target = await popup(page, 1);
  await page.getByLabel('Feature visibility', { exact: true }).selectOption('off');
  if (await page.getByRole('heading', { name: 'Feature Visibility Scope', exact: true }).isVisible()) {
    await page.getByText('This feature', { exact: true }).click();
  }
  await closeEditor(page);
  return target.featureId;
};
const history = async (page, action) => {
  await page.getByRole('button', { name: action, exact: true }).click();
  await expect.poll(() => page.evaluate(() => !window.__GBDRAW_HISTORY__.restoring.value
    && !window.__GBDRAW_HISTORY__.capturing.value)).toBe(true);
};

test('V1-V5/V8 individual visibility is reconciled on accepted source replacement and History replay', async ({ browser }, testInfo) => {
  test.setTimeout(1_800_000);
  const page = await load(browser);
  let fresh;
  try {
    await generate(page);
    const id = await hide(page);
    expect(id).toBe('f24a47546');
    await generate(page);
    expect((await snapshot(page)).featureVisibility).toEqual({ [id]: 'off' });
    await switchMode(page, 'linear');
    await switchMode(page, 'circular');
    expect((await snapshot(page)).featureVisibility).toEqual({ [id]: 'off' });
    await upload(page, await source(tobacco, 'tobacco.gb'));
    // A remains the current Result until B is admitted.
    expect((await snapshot(page)).featureVisibility).toEqual({ [id]: 'off' });
    await generate(page);
    expect((await snapshot(page)).featureVisibility).toEqual({});
    await history(page, 'Undo');
    expect((await snapshot(page)).featureVisibility).toEqual({ [id]: 'off' });
    await history(page, 'Redo');
    expect((await snapshot(page)).featureVisibility).toEqual({});
    const file = testInfo.outputPath('source-B.gbdraw-session.json.gz');
    const bytes = await download(page, 'Save Session', file);
    expect(JSON.parse(gunzipSync(bytes)).features.featureVisibilityOverrides).toEqual({});
    fresh = await load(browser, file);
    expect((await snapshot(fresh)).featureVisibility).toEqual({});
    await upload(fresh, await source(seed, 'mito-return.gb'));
    await generate(fresh);
    expect((await snapshot(fresh)).featureVisibility).toEqual({});
    expect(await fresh.evaluate(id => window.__GBDRAW_APP__.extractedFeatures
      .some(feature => feature.svg_id === id), id)).toBe(true);
    expect(page.externalRequests).toEqual([]);
    expect(fresh.externalRequests).toEqual([]);
  } finally {
    await page.context().close();
    if (fresh) await fresh.context().close();
  }
});

test('V6/V7 shared biological targets and manual matchers survive while rejected sources preserve A', async ({ browser }, testInfo) => {
  test.setTimeout(1_800_000);
  const page = await load(browser);
  try {
    await generate(page);
    const id = await hide(page);
    const a = await source(seed, 'renamed.gb');
    const b = await source(lambda, 'lambda.gb');
    await upload(page, { ...a, buffer: Buffer.concat([a.buffer, Buffer.from('\n'), b.buffer]) });
    await generate(page);
    expect((await snapshot(page)).featureVisibility).toEqual({ [id]: 'off' });
    // Use the existing manual-rule owner for a type-wide matcher.
    await page.evaluate(() => {
      const app = window.__GBDRAW_APP__;
      app.addFeatureVisibilityRule();
      for (const [field, value] of Object.entries({ featureType: 'CDS', qualifier: 'product', value: '.*', action: 'exclude_matching' })) {
        app.setFeatureVisibilityRuleField(0, field, value);
      }
    });
    const before = await snapshot(page);
    await upload(page, { name: 'corrupt.gb', mimeType: 'text/plain', buffer: Buffer.from('invalid GenBank source') });
    expect(await page.evaluate(() => window.__GBDRAW_APP__.runAnalysis())).toEqual({ status: 'error' });
    const failed = await snapshot(page);
    expect(failed.featureVisibility).toEqual(before.featureVisibility);
    expect(failed.visibilityRules).toEqual(before.visibilityRules);
    expect(failed.result).toEqual(before.result);
    expect(failed.request).toEqual(before.request);
    await upload(page, await source(tobacco, 'tobacco.gb'));
    await generate(page);
    const accepted = await snapshot(page);
    expect(accepted.featureVisibility).toEqual({});
    expect(accepted.visibilityRules).toEqual(before.visibilityRules);
    await fs.writeFile(testInfo.outputPath('manual-rule.json'), JSON.stringify(accepted.visibilityRules));
    expect(page.externalRequests).toEqual([]);
  } finally { await page.context().close(); }
});

test('composite source replacement and Web-CLI-Web replay do not restore discarded visibility', async ({ browser }, testInfo) => {
  test.setTimeout(1_800_000);
  const a = testInfo.outputPath('mito.gb');
  const b = testInfo.outputPath('lambda.gb');
  await fs.writeFile(a, (await source(seed, 'mito.gb')).buffer);
  await fs.writeFile(b, (await source(lambda, 'lambda.gb')).buffer);
  const cli = async (args, name) => {
    const output = testInfo.outputPath(`${name}.gbdraw-session.json.gz`);
    await promisify(execFile)('python', ['-m', 'gbdraw.cli', 'circular', ...args,
      '-o', testInfo.outputPath(name), '--session_output', output], {
      cwd: testInfo.outputDir, env: { ...process.env, PYTHONPATH: process.cwd() },
      timeout: 1_800_000, maxBuffer: 1_000_000
    });
    return output;
  };
  const composite = await cli(['--gbk', a, b, '--multi_record_canvas', '--labels', 'out'], 'composite');
  const page = await load(browser, composite);
  let fresh;
  try {
    await generate(page);
    const id = await hide(page);
    await generate(page);
    await switchMode(page, 'linear');
    await switchMode(page, 'circular');
    expect((await snapshot(page)).featureVisibility).toEqual({ [id]: 'off' });
    await upload(page, await source(tobacco, 'tobacco.gb'));
    await generate(page);
    expect((await snapshot(page)).featureVisibility).toEqual({});
    const saved = testInfo.outputPath('replaced.gbdraw-session.json.gz');
    await download(page, 'Save Session', saved);
    const replayed = await cli(['--session', saved], 'replayed');
    fresh = await load(browser, replayed);
    expect((await snapshot(fresh)).featureVisibility).toEqual({});
    await generate(fresh);
    expect((await snapshot(fresh)).featureVisibility).toEqual({});
    expect(page.externalRequests).toEqual([]);
    expect(fresh.externalRequests).toEqual([]);
  } finally {
    await page.context().close();
    if (fresh) await fresh.context().close();
  }
});

test('repeated generation of an unchanged source preserves hidden duplicate feature intent', async ({ browser }) => {
  test.setTimeout(1_800_000);
  const page = await load(browser);
  try {
    const file = await source(seed, 'duplicate-rRNA.gb');
    const text = file.buffer.toString('utf8');
    const entry = text.match(/     rRNA\s[^]*?(?=\n     \S|\nORIGIN)/)[0];
    await upload(page, { ...file, buffer: Buffer.from(text.replace(entry, `${entry}\n${entry}`)) });
    await generate(page);
    const id = await hide(page);
    expect(id).toContain('__instance_');
    await generate(page);
    await generate(page);
    expect((await snapshot(page)).featureVisibility).toEqual({ [id]: 'off' });
    // A new File view with the same supported semantic target is valid too.
    await upload(page, { ...file, name: 'renamed-duplicates.gb', buffer: Buffer.from(text.replace(entry, `${entry}\n${entry}`)) });
    await generate(page);
    expect((await snapshot(page)).featureVisibility).toEqual({ [id]: 'off' });
    expect(page.externalRequests).toEqual([]);
  } finally { await page.context().close(); }
});
