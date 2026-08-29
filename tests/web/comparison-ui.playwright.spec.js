const { test, expect } = require('@playwright/test');
const { readFileSync } = require('node:fs');
const { gunzipSync } = require('node:zlib');

const makeGenbank = (recordId, base = 'atg') => {
  const sequence = base.repeat(100);
  const origin = sequence.match(/.{1,60}/g).map((chunk, index) => {
    const groups = chunk.match(/.{1,10}/g).join(' ');
    return `${String(index * 60 + 1).padStart(9)} ${groups}`;
  }).join('\n');
  return `LOCUS       ${recordId.padEnd(24)} 300 bp    DNA     linear   UNA 01-JAN-2000
DEFINITION  comparison UI browser test.
ACCESSION   ${recordId}
VERSION     ${recordId}
KEYWORDS    .
SOURCE      synthetic construct
  ORGANISM  synthetic construct
            .
FEATURES             Location/Qualifiers
     CDS             1..90
                     /product="comparison UI protein"
ORIGIN
${origin}
//
`;
};

const openLinear = async (page) => {
  await page.goto('/gbdraw/web/index.html', { waitUntil: 'domcontentloaded' });
  await page.waitForFunction(() => window.__GBDRAW_APP__);
  await page.getByRole('button', { name: 'Linear', exact: true }).click();
  await expect(page.locator('[data-linear-comparison-card]')).toBeVisible();
};

const inputHeaderAdd = (page) => (
  page.locator('.card-header').filter({ hasText: 'Input Genomes' })
    .getByRole('button', { name: 'Add sequence' })
);

const comparisonCard = (page) => page.locator('[data-linear-comparison-card]');
const comparisonCommands = (page) => comparisonCard(page).getByRole('group', {
  name: 'Set all adjacent comparisons'
});
const comparisonSettings = (page) => page.locator(
  'details[data-linear-comparison-disclosure="settings"]'
);
const selectedPairs = (page) => page.locator(
  'details[data-linear-comparison-disclosure="selected-pairs"]'
);

const expectInside = async (child, parent) => {
  const [childBox, parentBox] = await Promise.all([
    child.boundingBox(),
    parent.boundingBox()
  ]);
  expect(childBox).not.toBeNull();
  expect(parentBox).not.toBeNull();
  expect(childBox.y).toBeGreaterThanOrEqual(parentBox.y - 1);
  expect(childBox.y + childBox.height).toBeLessThanOrEqual(
    parentBox.y + parentBox.height + 1
  );
};

const expectKeyboardFocusIndicator = async (locator) => {
  await expect(locator).toBeFocused();
  const indicator = await locator.evaluate((element) => {
    const style = getComputedStyle(element);
    return {
      focusVisible: element.matches(':focus-visible'),
      outlineStyle: style.outlineStyle,
      outlineWidth: style.outlineWidth
    };
  });
  expect(indicator.focusVisible).toBe(true);
  expect(indicator.outlineStyle).not.toBe('none');
  expect(indicator.outlineWidth).not.toBe('0px');
};

test('fresh Linear keeps primary input visible and uses command/status semantics', async ({ page }) => {
  await page.setViewportSize({ width: 1280, height: 720 });
  await openLinear(page);

  const settingsPane = page.locator('.settings-pane');
  const firstUploader = page.getByRole('button', {
    name: 'Choose GenBank File', exact: true
  }).first();
  await expectInside(inputHeaderAdd(page), settingsPane);
  await expectInside(firstUploader, settingsPane);

  const commands = comparisonCommands(page);
  const expectedNames = [
    'Set no comparison',
    'Run LOSAT for all adjacent pairs',
    'Use uploaded BLAST TSV for all adjacent pairs'
  ];
  await expect(commands.getByRole('button')).toHaveCount(3);
  for (const name of expectedNames) {
    const button = commands.getByRole('button', { name, exact: true });
    await expect(button).toBeVisible();
    await expect(button).not.toHaveAttribute('aria-pressed');
  }
  const currentStatus = comparisonCard(page).getByRole('status');
  await expect(currentStatus).toContainText('Current: No comparison');
  await expect(currentStatus).toContainText('No comparison');
  await expect(commands.getByRole('status')).toHaveCount(0);

  await expect(comparisonSettings(page)).toHaveCount(1);
  await expect(comparisonSettings(page)).not.toHaveAttribute('open', '');
  await expect(page.getByRole('button', { name: 'Generate Diagram' })).toHaveCount(1);

  const order = await page.evaluate(() => {
    const selectors = [
      '[data-linear-record-list] .upload-zone',
      '[data-linear-comparison-card]',
      '.basic-settings',
      '.generate-bar',
      '[data-linear-advanced-comparison]'
    ];
    const elements = selectors.map((selector) => document.querySelector(selector));
    return {
      present: elements.every(Boolean),
      ordered: elements.every((element, index) => (
        index === elements.length - 1 || Boolean(
          element.compareDocumentPosition(elements[index + 1])
          & Node.DOCUMENT_POSITION_FOLLOWING
        )
      )),
      pairInRecordList: document.querySelectorAll(
        '[data-linear-record-list] [data-edge-key], [data-linear-record-list] [data-linear-comparison-boundary]'
      ).length
    };
  });
  expect(order).toEqual({ present: true, ordered: true, pairInRecordList: 0 });

  await inputHeaderAdd(page).click();
  await expect(page.locator('[data-linear-record-card]')).toHaveCount(2);
  await expectInside(inputHeaderAdd(page), settingsPane);
  await expectInside(firstUploader, settingsPane);
  await expect(comparisonSettings(page)).not.toHaveAttribute('open', '');
  await expect(page.locator('[data-linear-record-list] [data-edge-key]')).toHaveCount(0);
});

test('uploader, comparison commands, and native summaries work from the keyboard', async ({ page }) => {
  await openLinear(page);
  await inputHeaderAdd(page).click();

  const uploaders = page.getByRole('button', {
    name: 'Choose GenBank File', exact: true
  });
  const firstUploader = uploaders.nth(0);
  await firstUploader.focus();
  await page.keyboard.press('Tab');
  await page.keyboard.press('Shift+Tab');
  await expectKeyboardFocusIndicator(firstUploader);

  const firstChooserPromise = page.waitForEvent('filechooser');
  await firstUploader.press('Enter');
  await (await firstChooserPromise).setFiles({
    name: 'keyboard-enter.gbk',
    mimeType: 'text/plain',
    buffer: Buffer.from(makeGenbank('KeyboardEnter'))
  });
  await expect(firstUploader).toContainText('keyboard-enter.gbk');

  const secondUploader = uploaders.nth(1);
  const secondChooserPromise = page.waitForEvent('filechooser');
  await secondUploader.press('Space');
  await (await secondChooserPromise).setFiles({
    name: 'keyboard-space.gbk',
    mimeType: 'text/plain',
    buffer: Buffer.from(makeGenbank('KeyboardSpace', 'gct'))
  });
  await expect(secondUploader).toContainText('keyboard-space.gbk');

  const commands = comparisonCommands(page);
  const runLosat = commands.getByRole('button', {
    name: 'Run LOSAT for all adjacent pairs'
  });
  await runLosat.focus();
  await expectKeyboardFocusIndicator(runLosat);
  await page.keyboard.press('Enter');
  await expect(comparisonCard(page).getByRole('status')).toContainText(
    'Current: Run LOSAT for all adjacent pairs'
  );

  const useUpload = commands.getByRole('button', {
    name: 'Use uploaded BLAST TSV for all adjacent pairs'
  });
  await useUpload.focus();
  await expectKeyboardFocusIndicator(useUpload);
  await page.keyboard.press('Space');
  await expect(comparisonCard(page).getByRole('status')).toContainText(
    'Current: Upload BLAST TSV for all adjacent pairs'
  );

  const noComparison = commands.getByRole('button', { name: 'Set no comparison' });
  await noComparison.focus();
  await expectKeyboardFocusIndicator(noComparison);
  await page.keyboard.press('Enter');
  await expect(comparisonCard(page).getByRole('status')).toContainText(
    'Current: No comparison'
  );

  await runLosat.focus();
  await page.keyboard.press('Space');
  const settingsSummary = comparisonSettings(page).locator('summary');
  await settingsSummary.focus();
  await expectKeyboardFocusIndicator(settingsSummary);
  await page.keyboard.press('Enter');
  await expect(comparisonSettings(page)).toHaveAttribute('open', '');

  const losatMode = page.getByRole('group', { name: 'LOSAT Mode' });
  const losatpButton = losatMode.getByRole('button', { name: 'LOSATP', exact: true });
  await losatpButton.focus();
  await expectKeyboardFocusIndicator(losatpButton);
  await page.keyboard.press('Space');
  await expect(losatpButton).toHaveAttribute('aria-pressed', 'true');
  await expect(losatMode.getByRole('button', {
    name: 'LOSATN', exact: true
  })).toHaveAttribute('aria-pressed', 'false');
  const losatpMode = page.getByRole('combobox', { name: 'LOSATP mode' });
  await expect(losatpMode).toBeVisible();
  await expect(losatpMode).toHaveValue('orthogroup');
  await losatpMode.focus();
  await expectKeyboardFocusIndicator(losatpMode);
  await page.keyboard.press('Home');
  await expect(losatpMode).toHaveValue('orthogroup');

  await settingsSummary.focus();
  await page.keyboard.press('Space');
  await expect(comparisonSettings(page)).not.toHaveAttribute('open', '');

  const selectedSummary = selectedPairs(page).locator('summary');
  await selectedSummary.focus();
  await expectKeyboardFocusIndicator(selectedSummary);
  await page.keyboard.press('Space');
  await expect(selectedPairs(page)).toHaveAttribute('open', '');

  const pairUpload = selectedPairs(page).getByRole('radio', {
    name: 'Upload BLAST TSV', exact: true
  }).first();
  await pairUpload.focus();
  await expectKeyboardFocusIndicator(pairUpload);
  await page.keyboard.press('Space');
  await expect(pairUpload).toBeChecked();

  await selectedSummary.focus();
  await page.keyboard.press('Enter');
  await expect(selectedPairs(page)).not.toHaveAttribute('open', '');

  const advanced = page.locator(
    'details[data-linear-comparison-disclosure="advanced"]'
  );
  const advancedSummary = advanced.locator('summary');
  await advancedSummary.focus();
  await expectKeyboardFocusIndicator(advancedSummary);
  await page.keyboard.press('Enter');
  await expect(advanced).toHaveAttribute('open', '');
  await page.keyboard.press('Space');
  await expect(advanced).not.toHaveAttribute('open', '');

  const serializedComparisonState = () => page.evaluate(async () => {
    const app = window.__GBDRAW_APP__;
    const {
      buildConfigData,
      buildEditorStateData,
      buildFeatureStateData,
      buildOrthogroupStateData,
      buildRunStateData,
      buildUiStateData
    } = await import('./js/services/config.js');
    const serialized = JSON.parse(JSON.stringify({
      config: buildConfigData(),
      ui: buildUiStateData(),
      features: buildFeatureStateData(),
      editorState: buildEditorStateData(),
      orthogroupState: buildOrthogroupStateData(),
      runState: buildRunStateData()
    }));
    const transientKeys = [];
    const visit = (value, path = '') => {
      if (!value || typeof value !== 'object') return;
      Object.entries(value).forEach(([key, child]) => {
        const childPath = path ? `${path}.${key}` : key;
        if (key === 'open' || key.toLowerCase().includes('disclosure')) {
          transientKeys.push(childPath);
        }
        visit(child, childPath);
      });
    };
    visit(serialized);
    return {
      semantic: JSON.parse(JSON.stringify({
        plan: app.linearComparisonPlan,
        program: app.losatProgram,
        blastnTask: app.losat.blastn.task,
        blastpMode: app.losat.blastp.mode,
        filters: {
          bitscore: app.adv.min_bitscore,
          evalue: app.adv.evalue,
          identity: app.adv.identity,
          length: app.adv.alignment_length
        },
        style: app.adv.pairwise_match_style
      })),
      serialized,
      transientKeys
    };
  });
  const beforeDisclosures = await serializedComparisonState();
  expect(beforeDisclosures.transientKeys).toEqual([]);

  const disclosureSummaries = [
    page.locator('details[data-linear-record-options]').first().locator(':scope > summary'),
    settingsSummary,
    selectedSummary,
    advancedSummary
  ];
  for (const summary of disclosureSummaries) {
    await summary.click();
    await summary.click();
  }
  const afterDisclosures = await serializedComparisonState();
  expect(afterDisclosures).toEqual(beforeDisclosures);
});

test('LOSAT and LOSATP modes own their controls and mixed plans require explicit topology change', async ({ page }) => {
  await openLinear(page);
  await page.evaluate(() => {
    const app = window.__GBDRAW_APP__;
    app.addLinearSeq();
    app.addLinearSeq();
    app.setLinearComparisonGlobalAction('losat');
  });
  await comparisonSettings(page).locator('summary').click();
  const losatMode = page.getByRole('group', { name: 'LOSAT Mode' });
  const losatpMode = page.getByRole('combobox', { name: 'LOSATP mode' });

  await expect(losatMode.getByRole('button')).toHaveText(['LOSATN', 'LOSATP', 'TLOSATX']);

  const chooseLosatMode = async (label) => {
    const button = losatMode.getByRole('button', { name: label, exact: true });
    await button.click();
    await expect(button).toHaveAttribute('aria-pressed', 'true');
  };

  const controlCounts = async () => ({
    blastnTask: await page.getByRole('combobox', { name: 'LOSATN task' }).count(),
    gencode: await page.locator(
      'input[aria-label^="TLOSATX gencode for sequence"]'
    ).count(),
    pairwiseHits: await page.getByRole('spinbutton', {
      name: 'Pairwise display max hits per protein'
    }).count(),
    groupHits: await page.getByRole('spinbutton', {
      name: 'Member hits per protein'
    }).count(),
    collinearGap: await page.getByRole('spinbutton', {
      name: 'Collinear max unit gap'
    }).count(),
    matchStyle: await page.getByRole('combobox', {
      name: 'Comparison match style'
    }).count(),
    matchHeight: await page.getByRole('spinbutton', {
      name: 'Comparison match height'
    }).count()
  });
  await chooseLosatMode('LOSATN');
  await expect(losatpMode).toHaveCount(0);
  expect(await controlCounts()).toEqual({
    blastnTask: 1, gencode: 0, pairwiseHits: 0, groupHits: 0, collinearGap: 0,
    matchStyle: 1, matchHeight: 1
  });
  await expect(page.getByRole('combobox', {
    name: 'LOSATN task'
  }).locator('option')).toHaveText(['megablast', 'blastn', 'dc-megablast']);

  await chooseLosatMode('TLOSATX');
  await expect(losatpMode).toHaveCount(0);
  expect(await controlCounts()).toEqual({
    blastnTask: 0, gencode: 3, pairwiseHits: 0, groupHits: 0, collinearGap: 0,
    matchStyle: 1, matchHeight: 1
  });
  const recordOptions = page.locator('details[data-linear-record-options]').first();
  const recordSummary = recordOptions.locator(':scope > summary');
  await recordSummary.focus();
  await page.keyboard.press('Tab');
  await page.keyboard.press('Shift+Tab');
  await expectKeyboardFocusIndicator(recordSummary);
  await page.keyboard.press('Enter');
  await expect(recordOptions).toHaveAttribute('open', '');
  await expect(page.getByRole('spinbutton', {
    name: 'TLOSATX gencode for sequence 1'
  })).toBeVisible();
  await page.keyboard.press('Space');
  await expect(recordOptions).not.toHaveAttribute('open', '');

  await chooseLosatMode('LOSATP');
  await expect(losatpMode).toBeVisible();
  await expect(losatpMode.locator('option')).toHaveText([
    'Similarity groups', 'Collinear blocks', 'Pairwise matches'
  ]);
  const proteinCases = [
    ['orthogroup', {
      blastnTask: 0, gencode: 0, pairwiseHits: 0, groupHits: 1, collinearGap: 0,
      matchStyle: 1, matchHeight: 1
    }],
    ['collinear', {
      blastnTask: 0, gencode: 0, pairwiseHits: 0, groupHits: 1, collinearGap: 1,
      matchStyle: 1, matchHeight: 1
    }],
    ['pairwise', {
      blastnTask: 0, gencode: 0, pairwiseHits: 1, groupHits: 0, collinearGap: 0,
      matchStyle: 1, matchHeight: 1
    }]
  ];
  for (const [value, expected] of proteinCases) {
    await losatpMode.selectOption(value);
    await expect(losatpMode).toHaveValue(value);
    expect(await controlCounts()).toEqual(expected);
    await expect(comparisonCard(page).locator('.linear-comparison-summary')).toContainText(
      `LOSATP · ${value === 'pairwise'
        ? 'Pairwise matches'
        : value === 'orthogroup'
          ? 'Similarity groups'
          : 'Collinear blocks'}`
    );
    if (value === 'collinear') {
      const scope = page.getByRole('combobox', {
        name: 'Collinear evidence scope'
      });
      await expect(scope).toHaveValue('adjacent');
      await expect(scope.locator('option:checked')).toHaveText('Adjacent pairs');
    }
  }

  const mixed = await page.evaluate(() => {
    const app = window.__GBDRAW_APP__;
    app.setLinearComparisonLosatMode('blastp');
    app.setLinearComparisonLosatpMode('pairwise');
    const [first, second] = app.linearComparisonResolution.edges;
    app.setLinearComparisonGapAction(first.edgeKey, 'upload');
    app.setLinearComparisonCardFile(first.edgeKey, new File(
      ['A\tB\t100\t30\t0\t0\t1\t30\t1\t30\t1e-20\t90\n'],
      'mixed-upload.tsv',
      { type: 'text/tab-separated-values' }
    ));
    return {
      mode: app.linearComparisonPlan.mode,
      edges: app.linearComparisonResolution.edges.map((edge) => edge.source),
      secondEdge: second.edgeKey
    };
  });
  expect(mixed.mode).toBe('selected');
  expect(mixed.edges).toEqual(['upload', 'losat']);
  await expect(comparisonCard(page).getByRole('status')).toContainText(
    'Current: Selected pairs (2; 1 LOSAT, 1 upload)'
  );
  await expect(comparisonCard(page).getByText('Custom', { exact: true })).toBeVisible();
  for (const modeKey of ['blastn', 'blastp', 'tblastx']) {
    await expect(losatMode.locator(
      `[data-linear-comparison-losat-mode-option="${modeKey}"]`
    )).toBeEnabled();
  }
  await expect(losatpMode.locator('option[value="pairwise"]')).toBeEnabled();
  await expect(losatpMode.locator('option[value="orthogroup"]')).toBeDisabled();
  await expect(losatpMode.locator('option[value="collinear"]')).toBeDisabled();
  expect(await page.evaluate(() => (
    window.__GBDRAW_APP__.setLinearComparisonLosatpMode('collinear')
  ))).toBe(false);
  await expect(losatMode.locator(
    '[data-linear-comparison-losat-mode-option="blastp"]'
  )).toHaveAttribute('aria-pressed', 'true');
  await expect(losatpMode).toHaveValue('pairwise');

  await comparisonSettings(page).getByRole('button', {
    name: 'Use all adjacent LOSAT'
  }).click();
  await expect(comparisonCard(page).getByRole('status')).toContainText(
    'Current: Run LOSAT for all adjacent pairs'
  );
  await expect(losatpMode.locator('option[value="orthogroup"]')).toBeEnabled();
  await expect(losatpMode.locator('option[value="collinear"]')).toBeEnabled();
});

test('LOSAT and LOSATP mode setters preserve inactive drafts, appearance, and filters', async ({ page }) => {
  await openLinear(page);
  const transitions = await page.evaluate(() => {
    const app = window.__GBDRAW_APP__;
    app.addLinearSeq();
    app.setLinearComparisonGlobalAction('losat');
    app.setLinearComparisonLosatMode('blastn');
    app.losat.blastn.task = 'dc-megablast';
    Object.assign(app.losat.blastp, {
      candidateLimit: 37,
      maxHits: 41,
      orthogroupMemberMaxHits: 43,
      collinearMinAnchors: 7,
      collinearMaxUnitGap: 11,
      collinearUnitMode: 'locus',
      collinearAnchorMode: 'all',
      collinearMergeOrientation: 'strand',
      collinearSearchScope: 'adjacent',
      collinearColorMode: 'orientation_identity'
    });
    Object.assign(app.adv, {
      min_bitscore: 97,
      evalue: '1e-45',
      identity: 76,
      alignment_length: 321,
      pairwise_match_style: 'ribbon',
      comparison_height: 73
    });
    const snapshot = () => JSON.parse(JSON.stringify({
      program: app.losatProgram,
      activeLosatMode: app.linearComparisonUi.activeLosatModeKey,
      activeLosatpMode: app.linearComparisonUi.activeLosatpModeKey,
      blastnTask: app.losat.blastn.task,
      blastpMode: app.losat.blastp.mode,
      blastpDraft: {
        candidateLimit: app.losat.blastp.candidateLimit,
        maxHits: app.losat.blastp.maxHits,
        orthogroupMemberMaxHits: app.losat.blastp.orthogroupMemberMaxHits,
        collinearMinAnchors: app.losat.blastp.collinearMinAnchors,
        collinearMaxUnitGap: app.losat.blastp.collinearMaxUnitGap,
        collinearUnitMode: app.losat.blastp.collinearUnitMode,
        collinearAnchorMode: app.losat.blastp.collinearAnchorMode,
        collinearMergeOrientation: app.losat.blastp.collinearMergeOrientation,
        collinearSearchScope: app.losat.blastp.collinearSearchScope,
        collinearColorMode: app.losat.blastp.collinearColorMode
      },
      filters: {
        bitscore: app.adv.min_bitscore,
        evalue: app.adv.evalue,
        identity: app.adv.identity,
        length: app.adv.alignment_length
      },
      appearance: {
        style: app.adv.pairwise_match_style,
        height: app.adv.comparison_height
      }
    }));
    const losatn = snapshot();
    const selectedLosatp = app.setLinearComparisonLosatMode('blastp');
    const selectedCollinear = app.setLinearComparisonLosatpMode('collinear');
    const collinear = snapshot();
    const selectedLosatn = app.setLinearComparisonLosatMode('blastn');
    const restoredLosatn = snapshot();
    const restoredLosatp = app.setLinearComparisonLosatMode('blastp');
    const restoredCollinear = snapshot();
    return {
      selectedLosatp,
      selectedCollinear,
      selectedLosatn,
      restoredLosatp,
      losatn,
      collinear,
      restoredLosatn,
      restoredCollinear
    };
  });

  expect(transitions.selectedLosatp).toBe(true);
  expect(transitions.selectedCollinear).toBe(true);
  expect(transitions.selectedLosatn).toBe(true);
  expect(transitions.restoredLosatp).toBe(true);
  expect(transitions.collinear).toMatchObject({
    program: 'blastp',
    activeLosatMode: 'blastp',
    activeLosatpMode: 'collinear',
    blastpMode: 'collinear'
  });
  expect(transitions.restoredLosatn).toMatchObject({
    program: 'blastn',
    activeLosatMode: 'blastn',
    activeLosatpMode: 'collinear',
    blastpMode: 'collinear'
  });
  expect(transitions.restoredCollinear).toMatchObject({
    program: 'blastp',
    activeLosatMode: 'blastp',
    activeLosatpMode: 'collinear',
    blastpMode: 'collinear'
  });
  for (const key of ['blastnTask', 'blastpDraft', 'filters', 'appearance']) {
    expect(transitions.collinear[key]).toEqual(transitions.losatn[key]);
    expect(transitions.restoredLosatn[key]).toEqual(transitions.losatn[key]);
    expect(transitions.restoredCollinear[key]).toEqual(transitions.losatn[key]);
  }
});

test('comparison controls drive appearance and current Session round trips', async ({ page }) => {
  test.setTimeout(600000);
  await openLinear(page);
  await page.evaluate((records) => {
    const app = window.__GBDRAW_APP__;
    if (app.linearSeqs.length < 2) app.addLinearSeq();
    records.forEach((content, index) => app.setLinearSeqPrimaryFile(
      index,
      'gb',
      new File([content], `ui04-record-${index + 1}.gbk`, {
        type: 'text/plain',
        lastModified: index + 1
      })
    ));
    Object.assign(app.form, {
      legend: 'none',
      show_gc: false,
      show_skew: false,
      show_depth: false,
      show_labels_linear: 'none'
    });
    app.setLinearComparisonGlobalAction('losat');
    app.setLinearComparisonLosatMode('blastp');
    app.setLinearComparisonLosatpMode('pairwise');
  }, [makeGenbank('Ui04A', 'atg'), makeGenbank('Ui04B', 'gct')]);

  await comparisonSettings(page).locator('summary').click();
  const advanced = page.locator(
    'details[data-linear-comparison-disclosure="advanced"]'
  );
  await advanced.locator('summary').click();

  const candidate = page.getByRole('spinbutton', { name: 'LOSATP Candidate limit' });
  const pairwiseMax = page.getByRole('spinbutton', {
    name: 'Pairwise display max hits per protein'
  });
  const memberHits = page.getByRole('spinbutton', { name: 'Member hits per protein' });
  const losatpMode = page.getByRole('combobox', { name: 'LOSATP mode' });
  const matchStyle = page.getByRole('combobox', { name: 'Comparison match style' });
  const matchHeight = page.getByRole('spinbutton', { name: 'Comparison match height' });

  await expect(candidate).toBeVisible();
  await expect(candidate).toHaveAttribute('placeholder', 'Unbounded');
  await candidate.fill('17');
  await candidate.focus();
  await expectKeyboardFocusIndicator(candidate);
  await pairwiseMax.fill('3');
  await matchHeight.fill('45');
  await matchStyle.focus();
  await page.keyboard.press('Home');
  await expect(matchStyle).toHaveValue('ribbon');
  await page.keyboard.press('ArrowDown');
  await expect(matchStyle).toHaveValue('curve');

  await losatpMode.selectOption('orthogroup');
  await expect(candidate).toBeVisible();
  await expect(memberHits).toBeVisible();
  await memberHits.fill('7');
  await expect(matchStyle).toHaveValue('curve');
  await expect(matchHeight).toHaveValue('45');

  await losatpMode.selectOption('collinear');
  await expect(candidate).toBeVisible();
  await expect(memberHits).toHaveValue('7');
  await page.getByRole('combobox', { name: 'Collinear evidence scope' }).selectOption('all');
  await page.getByRole('combobox', { name: 'Collinear color mode' })
    .selectOption('orientation_identity');
  await page.getByRole('combobox', { name: 'Collinear unit mode' }).selectOption('locus');
  await page.getByRole('combobox', { name: 'Collinear anchor mode' }).selectOption('all');
  await page.getByRole('combobox', { name: 'Collinear merge orientation' })
    .selectOption('strand');

  await losatpMode.selectOption('pairwise');
  await expect(pairwiseMax).toHaveValue('3');
  await expect(candidate).toHaveValue('17');
  await losatpMode.selectOption('collinear');
  await expect(memberHits).toHaveValue('7');
  await expect(page.getByRole('combobox', { name: 'Collinear unit mode' })).toHaveValue('locus');
  await expect(page.getByRole('combobox', { name: 'Collinear anchor mode' })).toHaveValue('all');
  await expect(page.getByRole('combobox', {
    name: 'Collinear merge orientation'
  })).toHaveValue('strand');

  await page.evaluate(() => {
    const app = window.__GBDRAW_APP__;
    const edge = app.linearComparisonResolution.edges[0];
    app.setLinearComparisonCardFile(edge.edgeKey, new File([
      'Ui04A\tUi04B\t95\t80\t4\t0\t1\t80\t5\t84\t1e-40\t160\n'
    ], 'ui04-comparison.tsv', {
      type: 'text/tab-separated-values',
      lastModified: 80
    }));
  });
  await expect(candidate).toHaveCount(0);
  await expect(matchStyle).toBeVisible();

  const renderAppearance = async (style, height) => page.evaluate(async ({ style_, height_ }) => {
    const app = window.__GBDRAW_APP__;
    app.adv.pairwise_match_style = style_;
    app.adv.comparison_height = height_;
    const result = await app.runAnalysis();
    const content = String(app.results?.[app.selectedResultIndex]?.content || '');
    const svg = new DOMParser().parseFromString(content, 'image/svg+xml');
    const match = svg.querySelector('[data-gbdraw-pairwise-match-id]');
    return {
      result,
      error: app.errorLog,
      style: match?.getAttribute('data-pairwise-match-style') || '',
      path: match?.getAttribute('d') || '',
      matchCount: svg.querySelectorAll('[data-gbdraw-pairwise-match-id]').length
    };
  }, { style_: style, height_: height });

  const ribbon = await renderAppearance('ribbon', 35);
  expect(ribbon.result, JSON.stringify(ribbon.error, null, 2)).toEqual({ status: 'ok' });
  expect(ribbon.style).toBe('ribbon');
  expect(ribbon.matchCount).toBeGreaterThan(0);

  const curve = await renderAppearance('curve', 35);
  expect(curve.result, JSON.stringify(curve.error, null, 2)).toEqual({ status: 'ok' });
  expect(curve.style).toBe('curve');
  expect(curve.path).not.toBe(ribbon.path);
  expect(curve.path).toContain('C');

  const taller = await renderAppearance('curve', 85);
  expect(taller.result, JSON.stringify(taller.error, null, 2)).toEqual({ status: 'ok' });
  expect(taller.style).toBe('curve');
  expect(taller.path).not.toBe(curve.path);

  await page.evaluate(() => {
    const app = window.__GBDRAW_APP__;
    app.setLinearComparisonGlobalAction('losat');
    app.setLinearComparisonLosatMode('blastp');
    app.setLinearComparisonLosatpMode('collinear');
    app.sessionTitle = 'ui04-comparison-controls';
  });
  const sessionDownloadPromise = page.waitForEvent('download', { timeout: 180000 });
  const saveResult = await page.evaluate(() => (
    window.__GBDRAW_APP__.saveSessionWithTitle()
  ));
  expect(saveResult).toMatchObject({ status: 'saved' });
  const sessionDownload = await sessionDownloadPromise;
  const sessionPath = await sessionDownload.path();
  const sessionBuffer = readFileSync(sessionPath);
  const session = JSON.parse(gunzipSync(sessionBuffer).toString('utf8'));
  expect(session.config.losat.blastp).toMatchObject({
    mode: 'collinear',
    candidateLimit: 17,
    maxHits: 3,
    orthogroupMemberMaxHits: 7,
    collinearUnitMode: 'locus',
    collinearAnchorMode: 'all',
    collinearMergeOrientation: 'strand',
    collinearSearchScope: 'all',
    collinearColorMode: 'orientation_identity'
  });
  expect(session.config.adv).toMatchObject({
    pairwise_match_style: 'curve',
    comparison_height: 85
  });

  await page.reload({ waitUntil: 'domcontentloaded' });
  await page.waitForFunction(() => window.__GBDRAW_APP__);
  const dialogPromise = page.waitForEvent('dialog', { timeout: 180000 });
  await page.locator(
    'input[type="file"][accept*="application/json"][accept*="application/gzip"]'
  ).setInputFiles({
    name: sessionDownload.suggestedFilename(),
    mimeType: 'application/gzip',
    buffer: sessionBuffer
  });
  const dialog = await dialogPromise;
  expect(dialog.message()).toBe('Session loaded successfully!');
  await dialog.accept();
  await page.waitForFunction(() => {
    const app = window.__GBDRAW_APP__;
    return app?.mode === 'linear'
      && app?.losatProgram === 'blastp'
      && app?.losat?.blastp?.mode === 'collinear'
      && app?.losat?.blastp?.collinearAnchorMode === 'all';
  }, null, { timeout: 180000 });
  expect(await page.evaluate(() => {
    const app = window.__GBDRAW_APP__;
    return {
      candidate: app.losat.blastp.candidateLimit,
      pairwise: app.losat.blastp.maxHits,
      member: app.losat.blastp.orthogroupMemberMaxHits,
      unit: app.losat.blastp.collinearUnitMode,
      anchor: app.losat.blastp.collinearAnchorMode,
      merge: app.losat.blastp.collinearMergeOrientation,
      scope: app.losat.blastp.collinearSearchScope,
      color: app.losat.blastp.collinearColorMode,
      style: app.adv.pairwise_match_style,
      height: app.adv.comparison_height
    };
  })).toEqual({
    candidate: 17,
    pairwise: 3,
    member: 7,
    unit: 'locus',
    anchor: 'all',
    merge: 'strand',
    scope: 'all',
    color: 'orientation_identity',
    style: 'curve',
    height: 85
  });
});

test('structured comparison errors open and focus their owning disclosure', async ({ page }) => {
  test.setTimeout(300000);
  await openLinear(page);

  const configureRecords = async () => page.evaluate((records) => {
    const app = window.__GBDRAW_APP__;
    app.lInputType = 'gb';
    if (app.linearSeqs.length < 2) app.addLinearSeq();
    records.forEach((content, index) => app.setLinearSeqPrimaryFile(
      index,
      'gb',
      new File([content], `error-record-${index + 1}.gbk`, {
        type: 'text/plain',
        lastModified: index + 1
      })
    ));
  }, [makeGenbank('ErrorA'), makeGenbank('ErrorB', 'gct')]);

  await configureRecords();
  const missingUpload = await page.evaluate(async () => {
    const app = window.__GBDRAW_APP__;
    app.setLinearComparisonGlobalAction('losat');
    const edgeKey = app.linearComparisonResolution.edges[0].edgeKey;
    app.setLinearComparisonGapAction(edgeKey, 'upload');
    document.querySelector('[data-linear-comparison-disclosure="selected-pairs"]')?.removeAttribute('open');
    const result = await app.runAnalysis();
    return {
      result,
      issueCodes: app.linearComparisonResolution.errors.map((issue) => issue.code),
      edgeKey
    };
  });
  expect(missingUpload.result).toEqual({ status: 'error' });
  expect(missingUpload.issueCodes).toContain('missing-upload');
  await expect(selectedPairs(page)).toHaveAttribute('open', '');
  await expect(comparisonCard(page).getByRole('status')).toContainText('comparison issue');
  const missingPair = page.locator(`[data-edge-key="${missingUpload.edgeKey}"]`);
  await expect(missingPair.locator('[data-linear-comparison-pair-upload]')).toContainText(
    'BLAST TSV'
  );
  expect(await missingPair.evaluate((element) => (
    element.contains(document.activeElement)
    && Boolean(document.activeElement?.closest('[data-linear-comparison-pair-upload]'))
  ))).toBe(true);

  await page.reload({ waitUntil: 'domcontentloaded' });
  await page.waitForFunction(() => window.__GBDRAW_APP__);
  await page.getByRole('button', { name: 'Linear', exact: true }).click();
  await configureRecords();
  const selectedCollinear = await page.evaluate(async () => {
    const app = window.__GBDRAW_APP__;
    app.setLinearComparisonGlobalAction('losat');
    app.setLinearComparisonLosatMode('blastp');
    app.setLinearComparisonLosatpMode('collinear');
    const edgeKey = app.linearComparisonResolution.edges[0].edgeKey;
    app.setLinearComparisonGapAction(edgeKey, 'losat');
    document.querySelector('[data-linear-comparison-disclosure="settings"]')?.removeAttribute('open');
    const result = await app.runAnalysis();
    return {
      result,
      issueCodes: app.linearComparisonResolution.errors.map((issue) => issue.code)
    };
  });
  expect(selectedCollinear.result).toEqual({ status: 'error' });
  expect(selectedCollinear.issueCodes).toContain('selected-losat-requires-pairwise');
  await expect(comparisonSettings(page)).toHaveAttribute('open', '');
  await expect(page.getByRole('combobox', { name: 'LOSATP mode' })).toBeFocused();
});

test('mobile layout has no overflow, fixed-action overlap, or semantic tab-order drift', async ({ page }) => {
  await page.setViewportSize({ width: 390, height: 844 });
  await openLinear(page);
  await inputHeaderAdd(page).click();
  await comparisonCommands(page).getByRole('button', {
    name: 'Run LOSAT for all adjacent pairs'
  }).click();

  const settingsSummary = comparisonSettings(page).locator('summary');
  await settingsSummary.click();
  const losatModeGeometry = await page.getByRole('group', {
    name: 'LOSAT Mode'
  }).evaluate((group) => {
    const buttons = [...group.querySelectorAll('button')];
    const tops = buttons.map((button) => button.getBoundingClientRect().top);
    return {
      buttonCount: buttons.length,
      sameRow: tops.every((top) => Math.abs(top - tops[0]) < 1),
      clientWidth: group.clientWidth,
      scrollWidth: group.scrollWidth
    };
  });
  expect(losatModeGeometry.buttonCount).toBe(3);
  expect(losatModeGeometry.sameRow).toBe(true);
  expect(losatModeGeometry.scrollWidth).toBeLessThanOrEqual(
    losatModeGeometry.clientWidth + 1
  );
  await settingsSummary.click();

  const geometry = await page.evaluate(() => {
    const pane = document.querySelector('.settings-pane');
    const comparison = document.querySelector('[data-linear-comparison-card]');
    const recordList = document.querySelector('[data-linear-record-list]');
    const basic = document.querySelector('.basic-settings');
    const generate = document.querySelector('.generate-bar');
    const advanced = document.querySelector('[data-linear-advanced-comparison]');
    const follows = (first, second) => Boolean(
      first.compareDocumentPosition(second) & Node.DOCUMENT_POSITION_FOLLOWING
    );
    return {
      viewport: window.innerWidth,
      documentWidth: document.documentElement.scrollWidth,
      pane: [pane.clientWidth, pane.scrollWidth],
      comparison: [comparison.clientWidth, comparison.scrollWidth],
      order: follows(recordList, comparison)
        && follows(comparison, basic)
        && follows(basic, generate)
        && follows(generate, advanced),
      pairInRecords: recordList.querySelectorAll('[data-edge-key]').length
    };
  });
  expect(geometry.documentWidth).toBeLessThanOrEqual(geometry.viewport + 1);
  expect(geometry.pane[1]).toBeLessThanOrEqual(geometry.pane[0] + 1);
  expect(geometry.comparison[1]).toBeLessThanOrEqual(geometry.comparison[0] + 1);
  expect(geometry.order).toBe(true);
  expect(geometry.pairInRecords).toBe(0);

  const firstUploader = page.getByRole('button', {
    name: 'Choose GenBank File', exact: true
  }).first();
  await firstUploader.focus();
  const tabSections = ['input'];
  for (let step = 0; step < 60; step += 1) {
    await page.keyboard.press('Tab');
    const section = await page.evaluate(() => {
      const active = document.activeElement;
      if (active?.closest('[data-linear-comparison-card]')) return 'comparison';
      if (active?.closest('.basic-settings')) return 'basic';
      if (active?.getAttribute('aria-label') === 'Generate Diagram') return 'generate';
      if (active?.closest('[data-linear-advanced-comparison]')) return 'advanced';
      if (active?.closest('[data-linear-record-list]')) return 'input';
      return '';
    });
    if (section && tabSections.at(-1) !== section) tabSections.push(section);
    if (section === 'advanced') break;
  }
  expect(tabSections).toEqual(['input', 'comparison', 'basic', 'generate', 'advanced']);

  const advancedSummary = page.locator(
    '[data-linear-advanced-comparison] summary[aria-label="Advanced comparison and layout"]'
  );
  await advancedSummary.evaluate((element) => element.scrollIntoView({ block: 'center' }));
  const [advancedBox, generateBox] = await Promise.all([
    advancedSummary.boundingBox(),
    page.locator('.generate-bar').boundingBox()
  ]);
  expect(advancedBox).not.toBeNull();
  expect(generateBox).not.toBeNull();
  const overlap = (
    advancedBox.x < generateBox.x + generateBox.width
    && advancedBox.x + advancedBox.width > generateBox.x
    && advancedBox.y < generateBox.y + generateBox.height
    && advancedBox.y + advancedBox.height > generateBox.y
  );
  expect(overlap).toBe(false);
});
