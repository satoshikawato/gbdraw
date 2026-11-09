const { test, expect } = require('@playwright/test');
const fs = require('node:fs');
const path = require('node:path');
const { gzipSync } = require('node:zlib');
const { execFile } = require('node:child_process');
const { promisify } = require('node:util');
const { evaluateWithRetainedPromise, reveal } = require('./helpers/app-lifecycle.cjs');
const { openFresh } = require('./helpers/audit-browser.cjs');
const { showDrawing } = require('./helpers/drawing-bar.cjs');
const { generate } = require('./helpers/mode-transition.cjs');
const { circularFileInput } = require('./helpers/circular-inputs.cjs');

// Every Session resource descriptor (`encoding: base64`) anywhere in the document.
const base64Resources = (value, found = []) => {
  if (Array.isArray(value)) value.forEach((member) => base64Resources(member, found));
  else if (value && typeof value === 'object') {
    if (value.encoding === 'base64' && typeof value.data === 'string') found.push(value);
    Object.values(value).forEach((member) => base64Resources(member, found));
  }
  return found;
};

// A gzip GenBank upload (NCBI's .gbff.gz) is inspected and drawn like the plain
// file; the page sends it to Python as chosen, and the saved Session holds the
// gzip bytes as chosen (receipt session.file-container), as the CLI's does.
test('a gzip GenBank upload lists its record, generates, and saves the gzip file as chosen', async ({ page }, info) => {
  test.setTimeout(240000);
  const gzipped = info.outputPath('NC_001416.gb.gz');
  const plain = fs.readFileSync(path.resolve('tests/test_inputs/NC_001416.gb'));
  const gzip = gzipSync(plain);
  fs.writeFileSync(gzipped, gzip);
  await openFresh(page);
  await circularFileInput(page, 'gb').setInputFiles(gzipped);
  await expect.poll(() => page.evaluate(() => (
    window.__GBDRAW_APP__.circularRecordList.map((record) => record.record_id)
  )), { timeout: 60_000 }).toEqual(['NC_001416.1']);
  await generate(page);
  expect(await page.evaluate(() => window.__GBDRAW_APP__.results.length)).toBe(1);
  // A titled Session saves without the title dialog.
  const pending = page.waitForEvent('download', { timeout: 60000 });
  await evaluateWithRetainedPromise(page, async () => {
    const app = window.__GBDRAW_APP__;
    app.sessionTitle = 'gzip-input';
    await app.saveSessionWithTitle();
  });
  const saved = await (await pending).path();
  const { readSavedSession } = await import(require('node:url').pathToFileURL(path.join(__dirname, 'helpers/session-documents.mjs')).href);
  const genbank = base64Resources(await readSavedSession(new Blob([fs.readFileSync(saved)])))
    .filter((resource) => resource.kind === 'genbank');
  expect(genbank).toHaveLength(1);
  expect(Buffer.from(genbank[0].data, 'base64').equals(gzip)).toBe(true);
  // The Web prefixes an embedded resource's name with its resource id.
  expect(genbank[0].name).toMatch(/NC_001416\.gb\.gz$/);
});

// A gzip upload that expands past the cap (Python's MAX_EXPANDED_INPUT_BYTES; a
// test hook stands in for 2 GiB) is refused with its name and the cap, and
// no record is listed (PD-OI-116 in docs/internal/OPTION_INTEGRITY_PRODUCT_CONTRACT.md).
test('a gzip upload past the expanded cap is refused with its name and the cap', async ({ page }, info) => {
  test.setTimeout(240000);
  const cap = 1024 * 1024;
  const gzipped = info.outputPath('big genome.gb.gz');
  fs.writeFileSync(gzipped, gzipSync(Buffer.alloc(cap + 1, 0x41)));
  await page.addInitScript((maxExpandedInputBytes) => {
    window.__GBDRAW_TEST_HOOKS__ = { ...(window.__GBDRAW_TEST_HOOKS__ || {}), maxExpandedInputBytes };
  }, cap);
  await openFresh(page);
  await circularFileInput(page, 'gb').setInputFiles(gzipped);
  await expect(page.getByText('An input could not be read. Replace or reselect it. File: big genome.gb.gz. '
    + 'gbdraw reads gzip files that expand to at most 1 MiB. Decompress the gzip file and choose the plain file.'))
    .toBeVisible({ timeout: 60_000 });
  expect(await page.evaluate(() => window.__GBDRAW_APP__.circularRecordList.length)).toBe(0);
});

// PD-OI-116: a gzip input of a kind no editor reads at once (a Depth TSV, a
// Linear BLAST TSV) is refused right after it is chosen too, before any
// Generate: the page counts its expanded size while it copies the file.
test('a gzip Depth TSV and BLAST TSV past the expanded cap are refused right after they are chosen', async ({ page }, info) => {
  test.setTimeout(240000);
  const cap = 1024 * 1024;
  const overCap = (name) => {
    const file = info.outputPath(name);
    fs.writeFileSync(file, gzipSync(Buffer.alloc(cap + 1, 0x41)));
    return file;
  };
  const refusal = (name) => page.getByText(`An input could not be read. Replace or reselect it. File: ${name}. `
    + 'gbdraw reads gzip files that expand to at most 1 MiB. Decompress the gzip file and choose the plain file.');
  await page.addInitScript((maxExpandedInputBytes) => {
    window.__GBDRAW_TEST_HOOKS__ = { ...(window.__GBDRAW_TEST_HOOKS__ || {}), maxExpandedInputBytes };
  }, cap);
  await openFresh(page);
  await (await reveal(page.getByLabel('Depth TSV', { exact: true }))).setInputFiles(overCap('big depth.tsv.gz'));
  await expect(refusal('big depth.tsv.gz')).toBeVisible({ timeout: 60_000 });

  await showDrawing(page, { kind: 'linear' });
  await page.evaluate(async (records) => {
    const app = window.__GBDRAW_APP__;
    app.errorLog = null;
    if (app.linearSeqs.length < 2) app.addLinearSeq();
    records.forEach((content, index) => app.setLinearSeqPrimaryFile(
      index, 'gb', new File([content], `gzip-cap-${index + 1}.gbk`, { type: 'text/plain', lastModified: index + 1 })));
    await app.setLinearComparisonGlobalAction('losat');
    app.setLinearComparisonGapAction(app.linearComparisonResolution.edges[0].edgeKey, 'upload');
    document.querySelector('[data-linear-comparison-disclosure="selected-pairs"]')?.setAttribute('open', '');
  }, ['A', 'B'].map((name) => fs.readFileSync(path.resolve('tests/test_inputs/NC_001416.gb'), 'utf8')
    .replace(/NC_001416(\.1)?/g, `GZCAP${name}`)));
  await page.getByLabel('BLAST TSV for #1 to #2', { exact: true }).setInputFiles(overCap('big hits.tsv.gz'));
  await expect(refusal('big hits.tsv.gz')).toBeVisible({ timeout: 60_000 });
  expect(await page.evaluate(() => window.__GBDRAW_APP__.results.length)).toBe(0);
});

// The expanded cap Python applies to a gzip input (gbdraw/io/input_text.py
// MAX_EXPANDED_INPUT_BYTES), as its refusal names it.
const PYTHON_CAP_MIB = 2048;

// A Linear Session the CLI wrote with a gzip BLAST TSV, whose table is then
// replaced by one that expands past Python's cap (2 GiB of zeros, about 2 MB
// as gzip, compressed in parts). argv: the directory to write in, the Session file.
const OVER_CAP_BLAST_SESSION = `
import gzip, sys, zlib
from pathlib import Path
from gbdraw import cli
from gbdraw.io.input_text import MAX_EXPANDED_INPUT_BYTES
from gbdraw.session_container import BytesResourceSource, read_session_container, write_session_container
directory, session = Path(sys.argv[1]), Path(sys.argv[2])
plain = (Path.cwd() / "tests/test_inputs/NC_001416.gb").read_text(encoding="utf-8")
for name in ("GZCAPA", "GZCAPB"):
    (directory / f"{name}.gb").write_text(plain.replace("NC_001416", name), encoding="utf-8")
(directory / "hits.tsv.gz").write_bytes(gzip.compress(b"GZCAPA.1\\tGZCAPB.1\\t99\\t500\\t5\\t0\\t1001\\t1500\\t2001\\t2500\\t1e-100\\t900\\n"))
written = directory / "written.gbdraw"
sys.argv = ["gbdraw", "linear", "--gbk", str(directory / "GZCAPA.gb"), str(directory / "GZCAPB.gb"), "-b",
            str(directory / "hits.tsv.gz"), "-o", str(directory / "first"), "-f", "svg", "--session_output", str(written)]
try:
    cli.main()
except SystemExit as exit_:
    assert exit_.code in (None, 0), exit_.code
container = read_session_container(written)
document = container.document
[blast] = [item["resourceId"] for item in document["drawings"][0]["renderRequest"]["comparisons"] if item["kind"] == "nucleotideBlast"]
compressor, zeros, left, parts = zlib.compressobj(6, zlib.DEFLATED, 31), bytes(8 * 1024 * 1024), MAX_EXPANDED_INPUT_BYTES + 1, []
while left > 0:
    parts.append(compressor.compress(zeros[:left]))
    left -= len(zeros)
parts.append(compressor.flush())
big = b"".join(parts)
document["resources"][blast]["size"] = len(big)
sources = {rid: container.resource(rid) for rid in container.resource_ids()}
sources[blast] = BytesResourceSource(big)
write_session_container(document, sources, session, overwrite=True)
print(document["resources"][blast]["name"])
`;

// PD-OI-116: a Session's inputs keep their path. A loaded gzip input that no
// page reader reads is refused by Python at Generate, with its name and the cap.
test('a loaded Session whose gzip BLAST TSV expands past the cap is refused at Generate with its name and the cap', async ({ page }, info) => {
  test.setTimeout(600_000);
  const directory = info.outputPath('over-cap-session');
  fs.mkdirSync(directory, { recursive: true });
  const session = path.join(directory, 'over-cap.gbdraw');
  const { stdout } = await promisify(execFile)('python', ['-c', OVER_CAP_BLAST_SESSION, directory, session], {
    env: { ...process.env, PYTHONPATH: process.cwd() }, timeout: 300_000, maxBuffer: 10_000_000
  });
  const name = stdout.trim().split('\n').pop();
  await page.addInitScript(() => {
    window.__GBDRAW_EXPANSIONS__ = 0;
    window.__GBDRAW_TEST_HOOKS__ = { ...(window.__GBDRAW_TEST_HOOKS__ || {}), onStructuralMetric: (metric) => {
      if (metric.name === 'inputExpansionCount' || metric.name === 'inputExpansionCheckCount') window.__GBDRAW_EXPANSIONS__ += 1;
    } };
  });
  await openFresh(page);
  await page.locator('input[accept*=".gbdraw"]').setInputFiles(session);
  await page.waitForFunction(() => !window.__GBDRAW_APP__.sessionImportPending
    && window.__GBDRAW_APP__.results.length > 0, null, { timeout: 300_000 });
  expect(await page.evaluate(() => window.__GBDRAW_APP__.errorLog)).toBeNull();
  // A read-only imported comparison (-b) is reused through Inherit (D-04).
  const resolution = page.locator('[data-imported-comparison-resolution]');
  await resolution.getByRole('button', { name: 'Inherit saved comparison' }).click();
  await expect(resolution).toContainText('INHERIT');
  await evaluateWithRetainedPromise(page, () => window.__GBDRAW_APP__.runAnalysis());
  await expect(page.getByText(`An input could not be read. Replace or reselect it. File: ${name}. `
    + `gbdraw reads gzip files that expand to at most ${PYTHON_CAP_MIB} MiB. Decompress the gzip file and choose the plain file.`))
    .toBeVisible({ timeout: 300_000 });
  expect(await page.evaluate(() => window.__GBDRAW_EXPANSIONS__), 'the page expanded nothing: Python refused it').toBe(0);
});

// PD-OI-116, OV-2946: an auxiliary table card (here Default colors) given a
// gzip file past the expanded cap shows one refusal, with the file's name and
// the cap, goes back to the file it had, and offers Retry, as its other import
// failures do: its import is the one owner of the outcome.
test('an auxiliary table gzip past the expanded cap is refused once, keeps the previous file, and offers Retry', async ({ page }, info) => {
  test.setTimeout(240000);
  const cap = 1024 * 1024;
  const big = info.outputPath('big colors.tsv.gz');
  fs.writeFileSync(big, gzipSync(Buffer.alloc(cap + 1, 0x41)));
  await page.addInitScript((maxExpandedInputBytes) => {
    window.__GBDRAW_TEST_HOOKS__ = { ...(window.__GBDRAW_TEST_HOOKS__ || {}), maxExpandedInputBytes };
  }, cap);
  await openFresh(page);
  const uploader = await reveal(page.getByLabel('Default colors file', { exact: true }));
  await uploader.setInputFiles({ name: 'first colors.tsv', mimeType: 'text/plain', buffer: Buffer.from('CDS\t#ff0000\n') });
  await expect.poll(() => page.evaluate(() => window.__GBDRAW_APP__.files.d_color?.name ?? null)).toBe('first colors.tsv');
  await page.evaluate(() => {
    window.__GBDRAW_ERRORS__ = [];
    window.Vue.watch(() => window.__GBDRAW_APP__.errorLog, (error) => { if (error) window.__GBDRAW_ERRORS__.push(error); });
  });
  await uploader.setInputFiles(big);
  const message = 'An input could not be read. Replace or reselect it. File: big colors.tsv.gz. '
    + 'gbdraw reads gzip files that expand to at most 1 MiB. Decompress the gzip file and choose the plain file.';
  await expect(page.getByText(message)).toBeVisible({ timeout: 60_000 });
  await expect(page.getByRole('button', { name: 'Retry table import', exact: true })).toBeVisible();
  await expect.poll(() => page.evaluate(() => window.__GBDRAW_APP__.files.d_color?.name ?? null)).toBe('first colors.tsv');
  // Let any later read of the refused file finish: the refusal is shown once and stays with its Retry.
  await page.waitForTimeout(1000);
  expect(await page.evaluate(() => window.__GBDRAW_ERRORS__.length), 'one refusal reached the error panel').toBe(1);
  await expect(page.getByRole('button', { name: 'Retry table import', exact: true })).toBeVisible();
  expect(await page.evaluate(() => window.__GBDRAW_APP__.files.d_color?.name ?? null)).toBe('first colors.tsv');
});

// PD-OI-116: the two comparison sequence fields of the Circular comparison ring
// refuse a gzip file past the expanded cap right after it is chosen, as the other
// file fields do. Add Seq reads the file at once (the ring reader names the row),
// so its read is the one owner of the refusal, with or without a copy store; any
// other read failure keeps the file-name default and Generate reports it (D12).
// The ring reader is Python's, whose cap (2 GiB) no test hook lowers.
const writeZerosGzip = async (file, size) => {
  const { pipeline } = require('node:stream/promises');
  const { Readable } = require('node:stream');
  const { createGzip } = require('node:zlib');
  const chunk = Buffer.alloc(1024 * 1024);
  async function* zeros() {
    for (let left = size; left > 0; left -= chunk.length) yield left >= chunk.length ? chunk : chunk.subarray(0, left);
  }
  await pipeline(Readable.from(zeros()), createGzip({ level: 6 }), fs.createWriteStream(file));
};

test('an Add Seq gzip past the expanded cap is refused right after it is chosen', async ({ page }, info) => {
  test.setTimeout(300_000);
  const big = info.outputPath('big ring.fasta.gz');
  await writeZerosGzip(big, PYTHON_CAP_MIB * 1024 * 1024 + 1);
  await openFresh(page);
  const losat = page.locator('input[type="radio"][value="losat"]');
  await reveal(losat);
  await losat.check();
  await reveal(page.locator('button', { hasText: 'Add comparison sequence' }));
  await page.locator('input[type="file"][accept^=".fa,.fas,.fasta,.fna,.ffn,.gb"]').setInputFiles(big);
  await expect(page.getByText('An input could not be read. Replace or reselect it. File: big ring.fasta.gz. '
    + `gbdraw reads gzip files that expand to at most ${PYTHON_CAP_MIB} MiB. Decompress the gzip file and choose the plain file.`))
    .toBeVisible({ timeout: 240_000 });
  expect(await page.evaluate(() => window.__GBDRAW_APP__.results.length)).toBe(0);
});

// The Comparison sequence beside an uploaded ring BLAST file: the copy the input
// store makes when it is chosen counts its expanded size (a test hook stands in
// for the 2 GiB cap), and the field shows that refusal.
test('a ring Comparison sequence gzip past the expanded cap is refused right after it is chosen', async ({ page }, info) => {
  test.setTimeout(240000);
  const cap = 1024 * 1024;
  const big = info.outputPath('big subject.fasta.gz');
  fs.writeFileSync(big, gzipSync(Buffer.alloc(cap + 1, 0x41)));
  await page.addInitScript((maxExpandedInputBytes) => {
    window.__GBDRAW_TEST_HOOKS__ = { ...(window.__GBDRAW_TEST_HOOKS__ || {}), maxExpandedInputBytes };
  }, cap);
  await openFresh(page);
  const upload = page.locator('input[type="radio"][value="upload"]');
  await reveal(upload);
  await upload.check();
  await page.locator('input[type="file"][aria-label="BLAST outfmt 6/7 files"]').setInputFiles({
    name: 'ring-hits.tsv', mimeType: 'text/plain', buffer: Buffer.from('q1\ts1\t99.0\t100\t1\t0\t1\t100\t1\t100\t1e-50\t200\n')
  });
  // The ring row's Add file (#1777): the uploaded table's comparison sequence.
  const companion = page.locator('[data-ring-row-files] input[type="file"]');
  await expect(companion).toHaveCount(1);
  await companion.setInputFiles(big);
  await expect(page.getByText('An input could not be read. Replace or reselect it. File: big subject.fasta.gz. '
    + 'gbdraw reads gzip files that expand to at most 1 MiB. Decompress the gzip file and choose the plain file.'))
    .toBeVisible({ timeout: 60_000 });
  expect(await page.evaluate(() => window.__GBDRAW_APP__.results.length)).toBe(0);
});

// D12: an Add Seq file whose record-label read fails for any other reason than
// the gzip cap shows nothing when it is chosen; its row keeps the file-name
// label, and Generate reports the read failure.
test('an Add Seq GenBank file that cannot be read keeps its file-name row and is reported at Generate', async ({ page }) => {
  test.setTimeout(300_000);
  await openFresh(page);
  await (await reveal(circularFileInput(page, 'gb'))).setInputFiles('tests/test_inputs/NC_001416.gb');
  await expect.poll(() => page.evaluate(() => window.__GBDRAW_APP__.circularRecordList.length), { timeout: 120_000 }).toBe(1);
  const losat = page.locator('input[type="radio"][value="losat"]');
  await reveal(losat);
  await losat.check();
  await reveal(page.locator('button', { hasText: 'Add comparison sequence' }));
  const steps = await page.evaluate(() => window.__GBDRAW_HISTORY__.getUndoCount());
  await page.locator('input[type="file"][accept^=".fa,.fas,.fasta,.fna,.ffn,.gb"]').setInputFiles({
    name: 'broken ring.gbk',
    mimeType: 'text/plain',
    // A LOCUS record that ends inside its features table: the ring reader cannot parse it.
    buffer: Buffer.from('LOCUS       BROKEN                    40 bp    DNA     linear   UNK 01-JAN-2000\n'
      + 'DEFINITION  broken ring.\nFEATURES             Location/Qualifiers\n     source          1..40\n//\n')
  });
  const rows = () => page.evaluate(() => window.__GBDRAW_APP__.circularConservation.series.map((entry) => entry.label));
  await expect.poll(rows).toEqual(['broken ring']);
  // The Add Seq step commits once the record-label read has ended.
  await expect.poll(() => page.evaluate(() => window.__GBDRAW_HISTORY__.getUndoCount()), { timeout: 120_000 })
    .toBe(steps + 1);
  expect(await page.evaluate(() => window.__GBDRAW_APP__.errorLog), 'no error banner when it is chosen').toBeNull();
  expect(await rows(), 'the row keeps its file-name label').toEqual(['broken ring']);
  await evaluateWithRetainedPromise(page, () => window.__GBDRAW_APP__.runAnalysis());
  await expect.poll(() => page.evaluate(() => window.__GBDRAW_APP__.errorLog?.summary ?? null), { timeout: 240_000 })
    .not.toBeNull();
  expect(await page.evaluate(() => {
    const { code, operation, context } = window.__GBDRAW_APP__.errorLog;
    return { code, operation, field: context.field };
  }), 'Generate reports the read failure').toEqual({ code: 'INPUT_UNREADABLE', operation: 'generate', field: 'comparison_sequence' });
  expect(await page.evaluate(() => window.__GBDRAW_APP__.results.length)).toBe(0);
});

// PD-OI-116 after a confirm dialog: a file whose hand-over to the input store
// waits for an answer (another drawing uses the replaced file; a stack row
// reaches the replaced ring) is refused right after the answer, once.
const watchErrors = (page) => page.evaluate(() => {
  window.__GBDRAW_ERRORS__ = [];
  window.Vue.watch(() => window.__GBDRAW_APP__.errorLog, (error) => { if (error) window.__GBDRAW_ERRORS__.push(error); });
});
const overCapFile = (info, name, cap) => {
  const file = info.outputPath(name);
  fs.writeFileSync(file, gzipSync(Buffer.alloc(cap + 1, 0x41)));
  return file;
};
const capRefusal = (page, name) => page.getByText(`An input could not be read. Replace or reselect it. File: ${name}. `
  + 'gbdraw reads gzip files that expand to at most 1 MiB. Decompress the gzip file and choose the plain file.');

test('a Depth TSV gzip past the expanded cap is refused once right after the drawing scope is answered', async ({ page }, info) => {
  test.setTimeout(240000);
  const cap = 1024 * 1024;
  await page.addInitScript((maxExpandedInputBytes) => {
    window.__GBDRAW_TEST_HOOKS__ = { ...(window.__GBDRAW_TEST_HOOKS__ || {}), maxExpandedInputBytes };
  }, cap);
  await openFresh(page);
  const depth = await reveal(page.getByLabel('Depth TSV', { exact: true }));
  await depth.setInputFiles({ name: 'first depth.tsv', mimeType: 'text/plain', buffer: Buffer.from('chr\t1\t5\n') });
  await expect(page.getByText('first depth.tsv', { exact: true }).first()).toBeVisible();
  const menu = page.locator('.drawing-bar-menu');
  await menu.locator('summary').click();
  await menu.getByRole('button', { name: 'Duplicate drawing', exact: true }).click();
  await watchErrors(page);
  await (await reveal(page.getByLabel('Depth TSV', { exact: true }))).setInputFiles(overCapFile(info, 'big depth.tsv.gz', cap));
  const dialog = page.getByRole('dialog', { name: 'Replace input file?' });
  await expect(dialog).toBeVisible();
  await dialog.getByRole('button', { name: 'This drawing only', exact: true }).click();
  await expect(capRefusal(page, 'big depth.tsv.gz')).toBeVisible({ timeout: 60_000 });
  await page.waitForTimeout(1000);
  expect(await page.evaluate(() => window.__GBDRAW_ERRORS__.length), 'one refusal reached the error panel').toBe(1);
  expect(await page.evaluate(() => window.__GBDRAW_APP__.results.length)).toBe(0);
});

test('a ring BLAST gzip past the expanded cap is refused once right after the stack row question is answered', async ({ page }, info) => {
  test.setTimeout(240000);
  const cap = 1024 * 1024;
  await page.addInitScript((maxExpandedInputBytes) => {
    window.__GBDRAW_TEST_HOOKS__ = { ...(window.__GBDRAW_TEST_HOOKS__ || {}), maxExpandedInputBytes };
  }, cap);
  await openFresh(page);
  const upload = page.locator('input[type="radio"][value="upload"]');
  await reveal(upload);
  await upload.check();
  const ringFiles = page.locator('input[type="file"][aria-label="BLAST outfmt 6/7 files"]');
  await ringFiles.setInputFiles({
    name: 'ring-hits.tsv', mimeType: 'text/plain', buffer: Buffer.from('q1\ts1\t99.0\t100\t1\t0\t1\t100\t1\t100\t1e-50\t200\n')
  });
  // A row of the custom stack reaches the ring's managed row, so replacing the ring file asks first.
  const ringRow = await evaluateWithRetainedPromise(page, async () => {
    const { state } = await import('./js/state.js');
    window.__GBDRAW_APP__.setCircularTrackSlotsEnabled(true);
    const rows = state.activeDrawing().adv.circular_track_slots;
    const ring = rows.find((row) => row.renderer === 'sequence_conservation');
    await window.__GBDRAW_HISTORY__.runUndoable('Setup', () => {
      rows.push({ id: 'mark', renderer: 'annotations', enabled: true, side: 'overlay', z: 1,
        params: { set_id: 'marks', anchor_slot: ring.id, layer: 'foreground', marks: ['highlight'] } });
    });
    return ring.id;
  });
  expect(ringRow).toBeTruthy();
  await watchErrors(page);
  await ringFiles.setInputFiles(overCapFile(info, 'big ring hits.tsv.gz', cap));
  const prompt = page.locator('[data-choice-prompt]');
  await expect(prompt).toBeVisible();
  await prompt.locator('[data-choice-prompt-action="hide"]').click();
  await expect(capRefusal(page, 'big ring hits.tsv.gz')).toBeVisible({ timeout: 60_000 });
  await page.waitForTimeout(1000);
  expect(await page.evaluate(() => window.__GBDRAW_ERRORS__.length), 'one refusal reached the error panel').toBe(1);
  expect(await page.evaluate(() => window.__GBDRAW_APP__.results.length)).toBe(0);
});
