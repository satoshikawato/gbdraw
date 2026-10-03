const { test, expect } = require('@playwright/test');
const { execFileSync } = require('node:child_process');
const { readFileSync, writeFileSync } = require('node:fs');
const { openApp } = require('./helpers/app-lifecycle.cjs');

// F-4 (OD-3, #694): a LOSAT derived payload holds the displayed direction of
// each record pair and, for Collinear blocks, whether a CLI grid limits the
// output to its displayed pairs. A display-only change reuses the raw searches,
// converts again, and draws what a fresh Generate draws.

// Every record carries the same three proteins at different CDS positions, so
// a link from another pair or direction lands on other coordinates.
const PROTEINS = { a: 'MKKKKKKKKK', b: 'MAAAAAAAAA', c: 'MWWWWWWWWW' };

const makeGenbank = (recordId, cdsStarts) => {
  const sequence = 'atg'.repeat(120);
  const origin = sequence.match(/.{1,60}/g).map((chunk, index) => {
    const groups = chunk.match(/.{1,10}/g).join(' ');
    return `${String(index * 60 + 1).padStart(9)} ${groups}`;
  }).join('\n');
  const cds = Object.keys(PROTEINS).map((suffix, index) => {
    const start = cdsStarts[index];
    return `     CDS             ${start}..${start + 89}
                     /locus_tag="${recordId}_${suffix}"
                     /protein_id="${recordId}_${suffix}"
                     /translation="${PROTEINS[suffix]}"`;
  }).join('\n');
  return `LOCUS       ${recordId.padEnd(24)} 360 bp    DNA     linear   UNA 01-JAN-2000
DEFINITION  display pair reuse browser test.
ACCESSION   ${recordId}
VERSION     ${recordId}
KEYWORDS    .
SOURCE      synthetic construct
  ORGANISM  synthetic construct
            .
FEATURES             Location/Qualifiers
${cds}
ORIGIN
${origin}
//
`;
};

// Deterministic LOSATP stand-in: a hit for every query and subject protein with
// the same sequence. The raw and derived caches, the converter, and the
// renderer all run unchanged.
const installProteinExecutor = (page) => page.addInitScript(() => {
  window.__GBDRAW_DISPLAY_PAIR_EXECUTOR_CALLS__ = 0;
  window.__GBDRAW_LOSAT_EXECUTOR__ = async (jobs, options) => {
    window.__GBDRAW_DISPLAY_PAIR_EXECUTOR_CALLS__ += 1;
    const proteins = (key) => [...String(options.sequences.get(key) || '')
      .matchAll(/^>(\S+)[^\n]*\n([^>]*)/gm)]
      .map((match) => [match[1], match[2].replace(/\s+/g, '')]);
    return jobs.map((job) => {
      const subjects = proteins(job.subjectSequenceKey);
      const rows = proteins(job.querySequenceKey).flatMap(([queryId, querySequence]) => subjects
        .filter(([, subjectSequence]) => subjectSequence === querySequence)
        .map(([subjectId]) => [queryId, subjectId, '100', '10', '0', '0', '1', '10', '1', '10',
          '1e-30', '200'].join('\t')));
      return { cacheKey: job.cacheKey, text: `${rows.join('\n')}\n` };
    });
  };
});

const chooseLosatp = async (page, mode) => {
  await page.getByRole('button', { name: 'Run LOSAT for all adjacent pairs' }).click();
  await page.getByRole('group', { name: 'LOSAT Mode' })
    .getByRole('button', { name: 'LOSATP', exact: true }).click();
  await page.getByRole('combobox', { name: 'LOSATP mode', exact: true }).selectOption(mode);
};

const openLayoutControls = (page) => (
  page.getByRole('button', { name: 'Advanced comparison and layout' }).press('Enter')
);

const setRows = async (page, rows) => {
  for (const [index, row] of rows.entries()) {
    await page.getByLabel(`Linear record row for sequence ${index + 1}`).fill(String(row));
  }
  await expect.poll(() => page.evaluate(() => (
    window.__GBDRAW_APP__.linearRecordRows.map((entry) => entry.row)
  ))).toEqual(rows);
};

const turnOffRows = async (page) => {
  await openLayoutControls(page);
  await page.getByLabel('Arrange linear records in rows').uncheck();
  await expect.poll(() => page.evaluate(() => window.__GBDRAW_APP__.linearRecordLayoutEnabled))
    .toBe(false);
};

const generate = async (page) => {
  const key = await page.evaluate(async () => (
    (await import('./js/state.js')).state.resultGenerationKey.value
  ));
  await page.getByRole('button', { name: 'Generate Diagram', exact: true }).click();
  await expect.poll(() => page.evaluate(async () => {
    const { state } = await import('./js/state.js');
    return {
      key: state.resultGenerationKey.value,
      processing: state.processing.value,
      error: state.errorLog.value?.summary || null
    };
  }), { timeout: 180000 }).toEqual({ key: key + 1, processing: false, error: null });
  return page.evaluate(() => {
    const app = window.__GBDRAW_APP__;
    const svg = new DOMParser().parseFromString(
      app.results[app.selectedResultIndex].content, 'image/svg+xml'
    );
    const attribute = (node, name) => node.getAttribute(name) || '';
    const telemetry = app.lastRunInfo?.losatTelemetry || {};
    const matches = [...svg.querySelectorAll('[data-pairwise-match-style]')];
    return {
      // Protein IDs are per-page runtime handles, so links compare by coordinates.
      links: matches.map((node) => [
        `${attribute(node, 'data-query-record-index')}->${attribute(node, 'data-subject-record-index')}`,
        `${attribute(node, 'data-qstart')}-${attribute(node, 'data-qend')}`,
        `${attribute(node, 'data-sstart')}-${attribute(node, 'data-send')}`
      ].join(' ')).sort(),
      paths: matches.map((node) => attribute(node, 'd')).sort(),
      cache: {
        cacheHits: telemetry.cacheHits,
        cacheMisses: telemetry.cacheMisses,
        derivedHits: telemetry.proteinDerivedPayloadCacheHits,
        derivedMisses: telemetry.proteinDerivedPayloadCacheMisses
      },
      executorCalls: window.__GBDRAW_DISPLAY_PAIR_EXECUTOR_CALLS__
    };
  });
};

// Runs the final settings in a new browser context, which has no LOSAT cache.
const freshRun = async (browser, open) => {
  const context = await browser.newContext();
  try {
    const page = await context.newPage();
    await open(page);
    return await generate(page);
  } finally {
    await context.close();
  }
};

const PAIR_SOURCES = [
  { name: 'display-pair-a.gbk', content: makeGenbank('PairA', [1, 121, 241]) },
  { name: 'display-pair-b.gbk', content: makeGenbank('PairB', [31, 151, 271]) }
];

// Two GenBank files, all-adjacent LOSATP Similarity groups, one record per row.
const openSimilarityGroups = async (page, rows) => {
  await installProteinExecutor(page);
  await openApp(page);
  await page.getByRole('button', { name: 'Linear', exact: true }).click();
  for (const [index, source] of PAIR_SOURCES.entries()) {
    if (index > 0) await page.getByRole('button', { name: 'Add sequence' }).first().click();
    const chooser = page.waitForEvent('filechooser');
    await page.getByRole('button', { name: 'Choose GenBank / DDBJ File' }).nth(index).click();
    await (await chooser).setFiles({
      name: source.name, mimeType: 'text/plain', buffer: Buffer.from(source.content)
    });
  }
  await expect.poll(() => page.evaluate(() => (
    window.__GBDRAW_APP__.linearSeqs.map((sequence) => sequence.gb?.name || '')
  ))).toEqual(PAIR_SOURCES.map((source) => source.name));
  await chooseLosatp(page, 'orthogroup');
  await openLayoutControls(page);
  await page.getByLabel('Arrange linear records in rows').check();
  await setRows(page, rows);
};

test('Swapping Linear rows reconverts LOSATP Similarity groups and draws the fresh ribbons', async ({ page, browser }) => {
  test.setTimeout(420000);
  await openSimilarityGroups(page, [1, 2]);
  const first = await generate(page);
  expect(first.cache).toEqual({ cacheHits: 0, cacheMisses: 4, derivedHits: 0, derivedMisses: 1 });
  expect(first.links).toEqual(['0->1 1-90 31-120', '0->1 121-210 151-240', '0->1 241-330 271-360']);

  // Only the display order changes: PairB moves to the top row.
  await setRows(page, [2, 1]);
  const swapped = await generate(page);

  const fresh = await freshRun(browser, (freshPage) => openSimilarityGroups(freshPage, [2, 1]));
  expect(fresh.cache).toEqual({ cacheHits: 0, cacheMisses: 4, derivedHits: 0, derivedMisses: 1 });
  expect(fresh.links).toEqual(['1->0 31-120 1-90', '1->0 151-240 121-210', '1->0 271-360 241-330'].sort());
  expect(fresh.paths).not.toEqual(first.paths);

  // The raw searches are reused; the derived payload is converted again.
  expect(swapped.executorCalls).toBe(first.executorCalls);
  expect(swapped.cache).toEqual({ cacheHits: 4, cacheMisses: 0, derivedHits: 0, derivedMisses: 1 });
  expect(swapped.links).toEqual(fresh.links);
  expect(swapped.paths).toEqual(fresh.paths);
});

// A CLI records table with row and column places GridA and GridB in row 1 and
// GridC in row 2. The Web keeps that grid, which displays 0-2 and 1-2.
const saveCliGridSession = (testInfo) => {
  const records = [['GridA', [1, 121, 241]], ['GridB', [31, 151, 271]], ['GridC', [61, 181, 271]]];
  const table = ['gbk\trow\tcolumn'];
  records.forEach(([recordId, starts], index) => {
    const path = testInfo.outputPath(`${recordId}.gbk`);
    writeFileSync(path, makeGenbank(recordId, starts));
    table.push(`${path}\t${index < 2 ? 1 : 2}\t${index < 2 ? index + 1 : 1}`);
  });
  const tablePath = testInfo.outputPath('records.tsv');
  writeFileSync(tablePath, `${table.join('\n')}\n`);
  const prefix = testInfo.outputPath('cli-grid');
  execFileSync('python', ['-m', 'gbdraw.cli', 'linear', '--records_table', tablePath,
    '--output', prefix, '--format', 'svg', '--save_session'], {
    env: { ...process.env, PYTHONPATH: process.cwd() }
  });
  return `${prefix}.gbdraw-session.json`;
};

const openCollinearGrid = async (page, sessionPath) => {
  await installProteinExecutor(page);
  await openApp(page);
  page.once('dialog', (dialog) => dialog.accept());
  await page.locator('input[accept^=".json,"]').first().setInputFiles(sessionPath);
  await expect.poll(() => page.evaluate(() => window.__GBDRAW_APP__.linearRecordRows.map(
    ({ row, canonicalRow, canonicalColumn }) => [row, canonicalRow, canonicalColumn]
  )), { timeout: 60000 }).toEqual([[1, 1, 1], [1, 1, 2], [2, 2, 1]]);
  await chooseLosatp(page, 'collinear');
  await page.getByRole('combobox', { name: 'Collinear evidence scope' }).selectOption('all');
};

test('Turning off a CLI grid reconverts LOSATP Collinear blocks and draws the fresh ribbons', async ({ page, browser }, testInfo) => {
  test.setTimeout(420000);
  const sessionPath = saveCliGridSession(testInfo);
  await openCollinearGrid(page, sessionPath);
  const grid = await generate(page);
  expect(grid.cache).toEqual({ cacheHits: 0, cacheMisses: 9, derivedHits: 0, derivedMisses: 1 });
  expect(grid.links).toEqual(['0->2 1-330 61-360', '1->2 31-360 61-360']);

  // One record per row: the displayed pairs become 0-1 and 1-2; the searches do not change.
  await turnOffRows(page);
  const rows = await generate(page);

  const fresh = await freshRun(browser, async (freshPage) => {
    await openCollinearGrid(freshPage, sessionPath);
    await turnOffRows(freshPage);
  });
  expect(fresh.cache).toEqual({ cacheHits: 0, cacheMisses: 9, derivedHits: 0, derivedMisses: 1 });
  expect(fresh.links).toEqual(['0->1 1-330 31-360', '1->2 31-360 61-360']);

  expect(rows.executorCalls).toBe(grid.executorCalls);
  expect(rows.cache).toEqual({ cacheHits: 9, cacheMisses: 0, derivedHits: 0, derivedMisses: 1 });
  expect(rows.links).toEqual(fresh.links);
  expect(rows.paths).toEqual(fresh.paths);
});

// PR5-B1: a CLI Session stores the records of one file in one resource with
// index selectors, so the Web plans the same source-file searches and every
// raw CLI search (same protein FASTA, same searchContext keys) is a cache hit.
const FAKE_CLI_LOSAT = `#!/usr/bin/env python3
import sys
args = sys.argv[1:]
if args in (["--version"], ["-version"]):
    print("losat 0.1.0")
    sys.exit(0)
if len(args) == 2 and args[1] == "--help":
    sys.exit(0)
def proteins(path):
    entries, current = {}, None
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if line.startswith(">"):
                current = line[1:].split()[0]
                entries[current] = ""
            elif current:
                entries[current] += line
    return entries
query = proteins(args[args.index("-query") + 1])
subject = proteins(args[args.index("-subject") + 1])
for query_id, query_sequence in query.items():
    for subject_id, subject_sequence in subject.items():
        if query_sequence == subject_sequence:
            print("\\t".join([query_id, subject_id, "100", "10", "0", "0", "1", "10", "1", "10",
                             "1e-30", "200"]))
`;

const saveCliMultiRecordSession = (testInfo) => {
  const multi = testInfo.outputPath('two-records.gbk');
  writeFileSync(multi, makeGenbank('FileA1', [1, 121, 241]) + makeGenbank('FileA2', [31, 151, 271]));
  const single = testInfo.outputPath('FileB.gbk');
  writeFileSync(single, makeGenbank('FileB', [61, 181, 271]));
  const losat = testInfo.outputPath('losat');
  writeFileSync(losat, FAKE_CLI_LOSAT, { mode: 0o755 });
  const tablePath = testInfo.outputPath('records.tsv');
  writeFileSync(tablePath, `gbk\trecord_id\n${multi}\tFileA1\n${multi}\tFileA2\n${single}\tFileB\n`);
  const prefix = testInfo.outputPath('cli-multi-record');
  execFileSync('python', ['-m', 'gbdraw.cli', 'linear', '--records_table', tablePath,
    '--losat', 'losatp', '--losatp_mode', 'similarity_groups', '--losat_bin', losat,
    '--losat_threads', '1', '--output', prefix, '--format', 'svg', '--save_session'], {
    env: { ...process.env, PYTHONPATH: process.cwd() }
  });
  return `${prefix}.gbdraw-session.json`;
};

test('The Web reuses every raw LOSATP search of a CLI Session over a multi-record file', async ({ page }, testInfo) => {
  test.setTimeout(420000);
  const sessionPath = saveCliMultiRecordSession(testInfo);
  await installProteinExecutor(page);
  await openApp(page);
  page.once('dialog', (dialog) => dialog.accept());
  await page.locator('input[accept^=".json,"]').first().setInputFiles(sessionPath);
  await expect.poll(() => page.evaluate(() => window.__GBDRAW_APP__.linearSeqs.length), { timeout: 60000 })
    .toBe(3);
  const run = await generate(page);
  // Nine record-pair results, all from the CLI raw entries; nothing is searched.
  expect({ hits: run.cache.cacheHits, misses: run.cache.cacheMisses }).toEqual({ hits: 9, misses: 0 });
  expect(run.executorCalls).toBe(0);
  expect(run.links.length).toBeGreaterThan(0);

  const session = JSON.parse(readFileSync(sessionPath, 'utf8'));
  expect(session.renderRequest.records.map((record) => [record.source.resourceId, record.selector]))
    .toEqual([
      ['record-1-genbank', { kind: 'recordIndex', index: 0 }],
      ['record-1-genbank', { kind: 'recordIndex', index: 1 }],
      ['record-3-genbank', null]
    ]);
  // Nine record-pair entries; the four between the files carry the file database.
  expect(session.losatCache.entries).toHaveLength(9);
  expect(session.losatCache.entries.filter((entry) => entry.searchContext)).toHaveLength(4);
});

