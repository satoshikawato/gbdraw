import assert from 'node:assert/strict';
import { cp, mkdir, mkdtemp, writeFile } from 'node:fs/promises';
import { tmpdir } from 'node:os';
import { join } from 'node:path';
import { pathToFileURL } from 'node:url';
import test from 'node:test';

const repoRoot = process.cwd();
const tempRoot = await mkdtemp(join(tmpdir(), 'gbdraw-comparison-ui-'));
const appRoot = join(tempRoot, 'js', 'app');
await mkdir(appRoot, { recursive: true });
await Promise.all([
  cp(
    join(repoRoot, 'gbdraw', 'web', 'js', 'app', 'comparison-ui.js'),
    join(appRoot, 'comparison-ui.js')
  ),
  cp(
    join(repoRoot, 'gbdraw', 'web', 'js', 'services', 'linear-comparisons.js'),
    join(tempRoot, 'js', 'services', 'linear-comparisons.js')
  ),
  cp(
    join(repoRoot, 'gbdraw', 'web', 'js', 'services', 'record-draw-selection.js'),
    join(tempRoot, 'js', 'services', 'record-draw-selection.js')
  ),
  cp(
    join(repoRoot, 'gbdraw', 'web', 'js', 'mode-profiles.js'),
    join(tempRoot, 'js', 'mode-profiles.js')
  ),
  cp(
    join(repoRoot, 'gbdraw', 'web', 'js', 'mode-profiles.generated.js'),
    join(tempRoot, 'js', 'mode-profiles.generated.js')
  ),
  // mode-profiles.js evaluates comparison thresholds with these dependency-free modules.
  cp(
    join(repoRoot, 'gbdraw', 'web', 'js', 'utils', 'error-normalization.js'),
    join(tempRoot, 'js', 'utils', 'error-normalization.js')
  ),
  cp(
    join(repoRoot, 'gbdraw', 'web', 'js', 'utils', 'optional-positive-number.js'),
    join(tempRoot, 'js', 'utils', 'optional-positive-number.js')
  ),
  cp(
    join(repoRoot, 'gbdraw', 'web', 'js', 'utils', 'i18n.js'),
    join(tempRoot, 'js', 'utils', 'i18n.js')
  ),
  cp(
    join(repoRoot, 'gbdraw', 'web', 'js', 'utils', 'download-names.js'),
    join(tempRoot, 'js', 'utils', 'download-names.js')
  )
]);
await writeFile(join(tempRoot, 'package.json'), '{"type":"module"}', 'utf8');

const {
  projectLinearComparisonLosatModeSelection,
  projectLinearComparisonLosatpModeSelection,
  projectLinearComparisonUi
} = await import(pathToFileURL(join(appRoot, 'comparison-ui.js')));
const {
  createDefaultLinearComparisonPlan,
  createLinearComparisonEdge,
  resolveLinearComparisonPlan
} = await import(pathToFileURL(join(tempRoot, 'js', 'services', 'linear-comparisons.js')));
const { comparisonStateForMode } = await import(
  pathToFileURL(join(tempRoot, 'js', 'mode-profiles.js'))
);

const sequences = [{ uid: 'a' }, { uid: 'b' }, { uid: 'c' }];
const file = (name) => ({ name, auditToken: `unread-${name}` });
const linearDefaults = comparisonStateForMode('linear');

const resolve = (plan, overrides = {}) => resolveLinearComparisonPlan({
  plan,
  sequences,
  losatProgram: 'blastn',
  blastpMode: 'orthogroup',
  ...overrides
});

const project = (plan, overrides = {}) => projectLinearComparisonUi({
  plan,
  resolution: resolve(plan, overrides),
  filters: linearDefaults,
  ...overrides
});

const modeOption = (options, key) => (
  options.find((option) => option.key === key)
);

test('fresh Linear state projects No comparison', () => {
  const plan = createDefaultLinearComparisonPlan();
  assert.equal(plan.mode, 'none');

  const projection = project(plan);
  assert.equal(projection.intentKey, 'none');
  assert.equal(projection.topologyKey, 'none');
  assert.equal(projection.topologyLabel, 'No comparison');
  assert.equal(projection.activePairCount, 0);
  assert.deepEqual(projection.sourceBreakdown, {
    key: 'none', losat: 0, upload: 0, label: 'No active pairs'
  });
  // The pressed action shows the intent; no status or summary text is projected.
  assert.equal(Object.hasOwn(projection, 'currentStatusLabel'), false);
  assert.equal(Object.hasOwn(projection, 'summaryText'), false);
  assert.deepEqual(projection.sectionKeys.settings, []);
  assert.deepEqual(projection.sectionKeys.selectedPairs, []);
  assert.deepEqual(projection.sectionKeys.advanced, ['record-layout']);
});

test('adjacent LOSAT keeps selection intent separate from resolved work', () => {
  const plan = { mode: 'adjacent', defaultSource: 'losat', edges: [] };
  const projection = project(plan);
  assert.equal(projection.intentKey, 'losat');
  assert.equal(projection.topologyLabel, 'All adjacent pairs');
  assert.equal(projection.topologyPairLabel, '2 adjacent pairs');
  assert.equal(projection.activePairCount, 2);
  assert.deepEqual(projection.sourceBreakdown, {
    key: 'losat', losat: 2, upload: 0, label: '2 LOSAT'
  });
  assert.equal(projection.activeLosatModeLabel, 'LOSATN');
  assert.match(projection.filterSummary, /^E-value <= /);

  const oneRecordResolution = resolveLinearComparisonPlan({
    plan,
    sequences: sequences.slice(0, 1),
    losatProgram: 'blastn',
    blastpMode: 'orthogroup'
  });
  const oneRecord = projectLinearComparisonUi({
    plan,
    resolution: oneRecordResolution,
    losatProgram: 'blastn',
    blastpMode: 'orthogroup',
    filters: linearDefaults
  });
  assert.equal(oneRecord.intentKey, 'losat');
  assert.equal(oneRecord.activePairCount, 0);
  assert.equal(oneRecord.sourceBreakdown.key, 'none');
  assert.equal(oneRecord.topologyPairLabel, '0 adjacent pairs');
});

test('adjacent upload summary uses only resolver-active uploads', () => {
  const first = file('a-b.tsv');
  const second = file('b-c.tsv');
  const plan = {
    mode: 'adjacent',
    defaultSource: 'upload',
    edges: [
      createLinearComparisonEdge({
        id: 'a-b', queryUid: 'a', subjectUid: 'b', source: 'upload',
        file: first, fileActive: true
      }),
      createLinearComparisonEdge({
        id: 'b-c', queryUid: 'b', subjectUid: 'c', source: 'upload',
        file: second, fileActive: true
      })
    ]
  };
  const projection = project(plan);
  assert.equal(projection.intentKey, 'upload');
  assert.equal(projection.activePairCount, 2);
  assert.deepEqual(projection.sourceBreakdown, {
    key: 'upload', losat: 0, upload: 2, label: '2 upload'
  });
  assert.equal(projection.topologyPairLabel, '2 adjacent pairs');
  // The pair rows and the issue list say which upload pairs lack a file.
  assert(!projection.sectionKeys.settings.includes('upload-readiness'));
  assert(!projection.sectionKeys.settings.includes('losat-mode'));
});

test('selected upload and mixed plans project custom intent and active source counts', () => {
  const selectedUploadPlan = {
    mode: 'selected',
    defaultSource: 'upload',
    edges: [createLinearComparisonEdge({
      id: 'a-b-upload', queryUid: 'a', subjectUid: 'b', source: 'upload',
      file: file('a-b.tsv'), fileActive: true
    })]
  };
  const selectedUpload = project(selectedUploadPlan);
  assert.equal(selectedUpload.intentKey, 'custom');
  assert.equal(selectedUpload.topologyLabel, 'Selected pairs');
  assert.equal(selectedUpload.activePairCount, 1);
  assert.deepEqual(selectedUpload.sourceBreakdown, {
    key: 'upload', losat: 0, upload: 1, label: '1 upload'
  });
  assert(!selectedUpload.sectionKeys.settings.includes('losat-mode'));

  const mixedPlan = {
    mode: 'selected',
    defaultSource: 'losat',
    edges: [
      createLinearComparisonEdge({
        id: 'a-b-losat', queryUid: 'a', subjectUid: 'b', source: 'losat'
      }),
      createLinearComparisonEdge({
        id: 'b-c-upload', queryUid: 'b', subjectUid: 'c', source: 'upload',
        file: file('b-c.tsv'), fileActive: true
      })
    ]
  };
  const mixed = project(mixedPlan);
  assert.equal(mixed.intentKey, 'custom');
  assert.equal(mixed.activePairCount, 2);
  assert.deepEqual(mixed.sourceBreakdown, {
    key: 'mixed', losat: 1, upload: 1, label: '1 LOSAT, 1 upload'
  });
  assert(mixed.sectionKeys.settings.includes('losat-mode'));
  assert(!mixed.sectionKeys.settings.includes('upload-readiness'));
});

test('LOSAT and LOSATP modes project only their active Settings and Advanced controls', () => {
  const plan = { mode: 'adjacent', defaultSource: 'losat', edges: [] };
  const cases = [
    {
      program: 'blastn',
      blastpMode: 'orthogroup',
      losatKey: 'blastn',
      losatpKey: 'orthogroup',
      settings: ['losat-mode', 'blastn-task', 'result-filters', 'comparison-appearance', 'losat-runtime'],
      absent: [
        'losatp-mode', 'record-genetic-codes', 'blastp-max-hits',
        'member-hits', 'collinear-primary'
      ]
    },
    {
      program: 'tblastx',
      blastpMode: 'orthogroup',
      losatKey: 'tblastx',
      losatpKey: 'orthogroup',
      settings: ['losat-mode', 'record-genetic-codes', 'result-filters', 'comparison-appearance', 'losat-runtime'],
      absent: [
        'losatp-mode', 'blastn-task', 'blastp-max-hits',
        'member-hits', 'collinear-primary'
      ]
    },
    {
      program: 'blastp',
      blastpMode: 'pairwise',
      losatKey: 'blastp',
      losatpKey: 'pairwise',
      settings: [
        'losat-mode', 'losatp-mode', 'blastp-max-hits',
        'result-filters', 'comparison-appearance', 'losatp-candidate', 'losat-runtime'
      ],
      absent: ['blastn-task', 'record-genetic-codes', 'member-hits', 'collinear-primary']
    },
    {
      program: 'blastp',
      blastpMode: 'orthogroup',
      losatKey: 'blastp',
      losatpKey: 'orthogroup',
      settings: [
        'losat-mode', 'losatp-mode',
        'result-filters', 'comparison-appearance', 'losatp-candidate', 'member-hits', 'losat-runtime'
      ],
      absent: ['blastn-task', 'record-genetic-codes', 'blastp-max-hits', 'collinear-primary']
    },
    {
      program: 'blastp',
      blastpMode: 'collinear',
      losatKey: 'blastp',
      losatpKey: 'collinear',
      settings: [
        'losat-mode', 'losatp-mode', 'collinear-primary',
        'result-filters', 'comparison-appearance', 'losatp-candidate', 'member-hits', 'losat-runtime'
      ],
      absent: [
        'blastn-task', 'record-genetic-codes', 'blastp-max-hits'
      ],
      advanced: 'collinear-details'
    }
  ];

  cases.forEach(({
    program, blastpMode, losatKey, losatpKey, settings, absent, advanced
  }) => {
    const projection = project(plan, { losatProgram: program, blastpMode });
    assert.equal(projection.activeLosatModeKey, losatKey);
    assert.equal(projection.activeLosatpModeKey, losatpKey);
    assert.deepEqual(projection.sectionKeys.settings, settings);
    absent.forEach((sectionKey) => {
      assert(!projection.sectionKeys.settings.includes(sectionKey));
    });
    assert.equal(
      projection.sectionKeys.advanced.includes('collinear-details'),
      advanced === 'collinear-details'
    );
    assert.equal(
      projection.sectionKeys.settings.includes('losatp-candidate'),
      program === 'blastp'
    );
  });

  const options = project(plan, { losatProgram: 'blastp', blastpMode: 'pairwise' });
  assert.deepEqual(
    options.losatModes.map(({ key, label }) => ({ key, label })),
    [
      { key: 'blastn', label: 'LOSATN' },
      { key: 'blastp', label: 'LOSATP' },
      { key: 'tblastx', label: 'TLOSATX' }
    ]
  );
  assert.deepEqual(
    options.losatpModes.map(({ key, label }) => ({ key, label })),
    [
      { key: 'orthogroup', label: 'Similarity groups' },
      { key: 'collinear', label: 'Collinear blocks' },
      { key: 'pairwise', label: 'Pairwise matches' }
    ]
  );

  const similarityGroups = project(plan, {
    losatProgram: 'blastp',
    blastpMode: 'orthogroup'
  });
  assert.deepEqual([similarityGroups.activeLosatModeLabel, similarityGroups.activeLosatpModeLabel],
    ['LOSATP', 'Similarity groups']);
  const collinear = project(plan, {
    losatProgram: 'blastp',
    blastpMode: 'collinear'
  });
  assert.deepEqual([collinear.activeLosatModeLabel, collinear.activeLosatpModeLabel],
    ['LOSATP', 'Collinear blocks']);
});

// OV-1128 B (D-114): Paralog links per group shows wherever it draws: Similarity groups, and Collinear with
// Build similarity groups with within-record searches. Pairwise, Collinear without inference, another program, and Selected pairs hide it.
test('Paralog links per group projects for Similarity groups and inferring Collinear only', () => {
  const adjacent = { mode: 'adjacent', defaultSource: 'losat', edges: [] };
  const selected = {
    mode: 'selected',
    defaultSource: 'losat',
    edges: [createLinearComparisonEdge({ id: 'a-b-losat', queryUid: 'a', subjectUid: 'b', source: 'losat' })]
  };
  const cases = [
    [adjacent, 'blastp', 'orthogroup', false, true],
    [adjacent, 'blastp', 'orthogroup', true, true],
    [adjacent, 'blastp', 'collinear', true, true],
    [adjacent, 'blastp', 'collinear', false, false],
    [adjacent, 'blastp', 'pairwise', true, false],
    [adjacent, 'blastn', 'orthogroup', true, false],
    [selected, 'blastp', 'orthogroup', true, false],
    [selected, 'blastp', 'collinear', true, false]
  ];
  for (const [plan, losatProgram, blastpMode, inferOrthogroups, shown] of cases) {
    const { advanced } = project(plan, { losatProgram, blastpMode, inferOrthogroups }).sectionKeys;
    assert.equal(advanced.includes('paralog-links'), shown, JSON.stringify([plan.mode, losatProgram, blastpMode, inferOrthogroups]));
  }
});

test('selected topology blocks only grouping and collinear LOSATP modes without changing the plan', () => {
  const plan = {
    mode: 'selected',
    defaultSource: 'losat',
    edges: [
      createLinearComparisonEdge({
        id: 'a-b-losat', queryUid: 'a', subjectUid: 'b', source: 'losat'
      }),
      createLinearComparisonEdge({
        id: 'b-c-upload', queryUid: 'b', subjectUid: 'c', source: 'upload',
        file: file('b-c.tsv'), fileActive: true
      })
    ]
  };
  const before = JSON.stringify(plan);
  const grouping = project(plan, { losatProgram: 'blastp', blastpMode: 'orthogroup' });
  assert.equal(grouping.activeLosatModeKey, 'blastp');
  assert.equal(grouping.activeLosatpModeKey, 'orthogroup');
  assert(grouping.losatModes.every((option) => option.selectable));
  assert.equal(modeOption(grouping.losatpModes, 'orthogroup').selectable, false);
  assert.equal(modeOption(grouping.losatpModes, 'collinear').selectable, false);
  assert.equal(modeOption(grouping.losatpModes, 'pairwise').selectable, true);
  assert(grouping.sectionKeys.settings.includes('losatp-mode'));
  assert(!grouping.sectionKeys.settings.includes('member-hits'));
  assert(grouping.sectionKeys.settings.includes('comparison-appearance'));
  assert(grouping.sectionKeys.settings.includes('losatp-candidate'));
  assert.equal(grouping.issues[0]?.disclosureKey, 'settings');

  for (const modeKey of ['orthogroup', 'collinear']) {
    const selection = projectLinearComparisonLosatpModeSelection({ plan, modeKey });
    assert.equal(selection.selectable, false);
    assert.equal(selection.reasonKey, 'requires-adjacent-losat');
    assert.equal(selection.requiredIntentActionKey, 'losat');
    assert.equal(selection.patch, null);
  }
  assert.equal(JSON.stringify(plan), before);

  const pairwise = projectLinearComparisonLosatpModeSelection({
    plan,
    modeKey: 'pairwise'
  });
  assert.deepEqual(pairwise.patch, { blastpMode: 'pairwise' });
  const protein = projectLinearComparisonLosatModeSelection({
    plan,
    modeKey: 'blastp'
  });
  assert.deepEqual(protein.patch, { losatProgram: 'blastp' });
  const nucleotide = projectLinearComparisonLosatModeSelection({
    plan,
    modeKey: 'blastn'
  });
  assert.deepEqual(nucleotide.patch, { losatProgram: 'blastn' });

  const adjacentCollinear = projectLinearComparisonLosatpModeSelection({
    plan: { mode: 'adjacent', defaultSource: 'losat', edges: [] },
    modeKey: 'collinear'
  });
  assert.deepEqual(adjacentCollinear.patch, { blastpMode: 'collinear' });
  assert(!Object.hasOwn(adjacentCollinear.patch, 'pairwiseMatchStyle'));
});

test('filter summary always reports actual values and compares against mode-profile defaults', () => {
  const plan = { mode: 'adjacent', defaultSource: 'losat', edges: [] };
  const defaultProjection = project(plan);
  assert.equal(defaultProjection.filterSummaryIsDefault, true);
  assert.equal(
    defaultProjection.filterSummary,
    `E-value <= ${linearDefaults.evalue} · ` +
      `Bitscore >= ${linearDefaults.min_bitscore} · ` +
      `Identity >= ${linearDefaults.identity}% · ` +
      `Length >= ${linearDefaults.alignment_length}`
  );

  const custom = project(plan, {
    filters: {
      evalue: '1e-20',
      min_bitscore: 85,
      identity: 97.5,
      alignment_length: 120
    }
  });
  assert.equal(custom.filterSummaryIsDefault, false);
  assert.equal(
    custom.filterSummary,
    'E-value <= 1e-20 · Bitscore >= 85 · Identity >= 97.5% · Length >= 120'
  );
});

test('retained dormant draft count does not become an active pair count', () => {
  const retained = file('retained.tsv');
  const dormantPlan = {
    mode: 'none',
    defaultSource: 'losat',
    edges: [
      createLinearComparisonEdge({
        id: 'dormant', queryUid: 'a', subjectUid: 'b', source: 'upload',
        file: retained, fileActive: false,
        losatFilename: 'retained.raw.tsv', losatFilenameActive: false
      }),
      createLinearComparisonEdge({
        id: 'empty', queryUid: 'b', subjectUid: 'c', source: 'losat', included: false
      })
    ]
  };
  const dormant = project(dormantPlan);
  assert.equal(dormant.activePairCount, 0);
  assert.equal(dormant.retainedDormantDraftCount, 1);
  assert.deepEqual(dormant.sectionKeys.selectedPairs, ['pair-editor', 'retained-drafts']);
  // A retained Raw LOSAT filename sits in its pair row, not in a second list.
  assert(!dormant.sectionKeys.advanced.includes('raw-results'));

  const activeUploadPlan = {
    mode: 'selected',
    defaultSource: 'upload',
    edges: [createLinearComparisonEdge({
      id: 'active', queryUid: 'a', subjectUid: 'b', source: 'upload',
      file: retained, fileActive: true
    })]
  };
  const activeUpload = project(activeUploadPlan);
  assert.equal(activeUpload.activePairCount, 1);
  assert.equal(activeUpload.retainedDormantDraftCount, 0);

  activeUploadPlan.edges[0].losatFilename = 'alternate.raw.tsv';
  activeUploadPlan.edges[0].losatFilenameActive = false;
  const alternatePayload = project(activeUploadPlan);
  assert.equal(alternatePayload.activePairCount, 1);
  assert.equal(alternatePayload.retainedDormantDraftCount, 1);
});

test('structured resolver issue codes project one issue list with focus targets and fixing actions', () => {
  const issueCodes = [
    'missing-upload',
    'selected-losat-requires-pairwise',
    'duplicate',
    'missing-uid',
    'self',
    'same-row',
    'non-adjacent',
    'unowned-error'
  ];
  const plan = {
    mode: 'selected',
    defaultSource: 'losat',
    edges: [createLinearComparisonEdge({
      id: 'problem', queryUid: 'a', subjectUid: 'b', source: 'losat'
    })]
  };
  const projection = projectLinearComparisonUi({
    plan,
    resolution: {
      edges: [],
      errors: issueCodes.map((code) => ({
        code,
        edgeId: 'problem',
        edgeKey: 'a->b',
        message: `wording must not route ${code}`
      }))
    },
    losatProgram: 'blastp',
    blastpMode: 'pairwise',
    filters: linearDefaults
  });

  // One list entry per routed issue, with its text and the actions that fix it.
  assert.equal(projection.issues.length, 7);
  const missing = projection.issues.find((issue) => issue.issueCode === 'missing-upload');
  assert.deepEqual(
    { ...missing, actions: missing.actions.map((action) => action.key) },
    {
      key: '0-missing-upload',
      issueCode: 'missing-upload',
      disclosureKey: 'selected-pairs',
      focusTargetKey: 'pair-upload',
      edgeId: 'problem',
      edgeKey: 'a->b',
      text: 'wording must not route missing-upload',
      pairText: '',
      actions: ['choose-file']
    }
  );
  assert.equal(missing.actions[0].ariaLabel, 'Choose file');
  const topology = projection.issues.find((issue) => issue.issueCode === 'selected-losat-requires-pairwise');
  assert.equal(topology.disclosureKey, 'settings');
  assert.equal(topology.focusTargetKey, 'losatp-mode');
  assert.deepEqual(topology.actions.map((action) => action.key), ['all-losat', 'pairwise']);
  ['duplicate', 'missing-uid', 'self', 'same-row', 'non-adjacent'].forEach((code) => {
    const target = projection.issues.find((entry) => entry.issueCode === code);
    assert.equal(target.disclosureKey, 'selected-pairs');
    assert.equal(target.focusTargetKey, 'pair-row');
    assert.deepEqual(target.actions.map((action) => action.key), ['show']);
  });
  assert(!projection.issues.some((target) => target.issueCode === 'unowned-error'));
});

test('a pair issue names its pair once in the list and marks the pair row', () => {
  const projection = projectLinearComparisonUi({
    plan: { mode: 'selected', defaultSource: 'upload', edges: [createLinearComparisonEdge({ id: 'p', queryUid: 'a', subjectUid: 'b', source: 'upload' })] },
    resolution: {
      edges: [{ id: 'p', edgeKey: 'a->b', queryIndex: 0, subjectIndex: 1, source: 'upload' }],
      errors: [{ code: 'missing-upload', edgeId: 'p', edgeKey: 'a->b', message: 'Every selected uploaded comparison requires an active BLAST TSV file.' }]
    },
    filters: linearDefaults
  });
  assert.equal(projection.issues[0].text, '#1 → #2: No BLAST file.');
  assert.deepEqual(projection.pairIssueText, { 'a->b': 'No BLAST file' });
  // The action's accessible name carries its pair (repeated buttons stay distinct).
  assert.deepEqual(projection.issues[0].actions.map(({ label, ariaLabel }) => [label, ariaLabel]),
    [['Choose file', 'Choose file for #1 → #2']]);
});

test('the LOSATP-mode issue names the pairs set to Upload BLAST TSV', () => {
  const projection = projectLinearComparisonUi({
    plan: { mode: 'selected', defaultSource: 'losat', edges: [] },
    resolution: {
      edges: [
        { id: 'p', edgeKey: 'a->b', queryIndex: 0, subjectIndex: 1, source: 'upload' },
        { id: 'q', edgeKey: 'b->c', queryIndex: 1, subjectIndex: 2, source: 'losat' }
      ],
      errors: [
        { code: 'selected-losat-requires-pairwise', edgeId: '', edgeKey: '', message: 'unused' },
        { code: 'missing-upload', edgeId: 'p', edgeKey: 'a->b', message: 'unused' },
        { code: 'duplicate', edgeId: 'p', edgeKey: 'a->b', message: 'Duplicate pair' }
      ]
    },
    losatProgram: 'blastp',
    blastpMode: 'collinear',
    filters: linearDefaults
  });
  assert.equal(projection.issues[0].text, 'Pair #1 → #2 is set to Upload BLAST TSV, so not every pair runs LOSAT. '
    + 'Similarity groups and Collinear blocks compare every adjacent pair with LOSAT.');
  // Every issue of a pair reaches its row, and list keys stay distinct.
  assert.deepEqual(projection.pairIssueText, { 'a->b': 'No BLAST file; Duplicate pair' });
  assert.equal(new Set(projection.issues.map((issue) => issue.key)).size, 3);
});

test('the LOSATP-mode issue of a Selected plan whose pairs all run LOSAT states its cause', () => {
  const projection = projectLinearComparisonUi({
    plan: { mode: 'selected', defaultSource: 'losat', edges: [] },
    resolution: {
      edges: [
        { id: 'p', edgeKey: 'a->b', queryIndex: 0, subjectIndex: 1, source: 'losat' },
        { id: 'q', edgeKey: 'b->c', queryIndex: 1, subjectIndex: 2, source: 'losat' }
      ],
      errors: [{ code: 'selected-losat-requires-pairwise', edgeId: '', edgeKey: '', message: 'unused' }]
    },
    losatProgram: 'blastp',
    blastpMode: 'collinear',
    filters: linearDefaults
  });
  assert.equal(projection.issues[0].text, 'The pairs are set one by one in Selected pairs. '
    + 'Similarity groups and Collinear blocks compare every adjacent pair with LOSAT.');
});

test('projection does not mutate plan, resolution, filters, or retained file references', () => {
  const retained = file('immutable.tsv');
  const plan = {
    mode: 'selected',
    defaultSource: 'upload',
    edges: [createLinearComparisonEdge({
      id: 'immutable', queryUid: 'a', subjectUid: 'b', source: 'upload',
      file: retained, fileActive: true,
      losatFilename: 'immutable.raw.tsv', losatFilenameActive: false
    })]
  };
  const resolution = resolve(plan);
  const filters = {
    evalue: '1e-9', min_bitscore: 60, identity: 75, alignment_length: 20
  };
  const before = JSON.stringify({ plan, resolution, filters });
  const projection = projectLinearComparisonUi({
    plan,
    resolution,
    losatProgram: 'blastn',
    blastpMode: 'collinear',
    filters
  });

  assert.equal(JSON.stringify({ plan, resolution, filters }), before);
  assert.equal(plan.edges[0].file, retained);
  assert.equal(plan.edges[0].losatFilename, 'immutable.raw.tsv');
  assert(Object.isFrozen(projection));
  assert(Object.isFrozen(projection.sectionKeys.settings));
  assert(Object.isFrozen(projection.issues));
});

test('an adjacent plan that resolves no pair says why Generate draws no comparison (CI-07a, CI-07b)', () => {
  // a) Run LOSAT with every record on one Linear row: no adjacent pair.
  const sameRow = project({ mode: 'adjacent', defaultSource: 'losat', edges: [] }, {
    layout: sequences.map(({ uid }) => ({ uid, row: 1 }))
  });
  assert.equal(sameRow.activePairCount, 0);
  assert.match(sameRow.emptyPlanNotice, /^No record pair to compare\. Records on the same Linear row are not compared/);
  // b) Upload BLAST TSV for all adjacent pairs, but no pair has a file.
  const noFile = project({ mode: 'adjacent', defaultSource: 'upload', edges: [] });
  assert.match(noFile.emptyPlanNotice, /^No pair has a BLAST TSV yet, so Generate draws no comparison\./);
  // Pairs to compare, No comparison, and an explicit Selected plan need no notice.
  assert.equal(project({ mode: 'adjacent', defaultSource: 'losat', edges: [] }).emptyPlanNotice, '');
  assert.equal(project(createDefaultLinearComparisonPlan()).emptyPlanNotice, '');
  assert.equal(project({ mode: 'selected', defaultSource: 'upload', edges: [] }).emptyPlanNotice, '');
});

test('a Selected plan of exactly the adjacent pairs from uploads keeps the Upload intent (UJ-05)', () => {
  // Uploading a table into the only pair turns All adjacent pairs into a
  // Selected plan with that one included edge.
  const pair = createLinearComparisonEdge({
    id: 'a-b', queryUid: 'a', subjectUid: 'b', source: 'upload', file: file('a-b.tsv'), fileActive: true
  });
  const plan = { mode: 'selected', defaultSource: 'upload', edges: [pair] };
  const intent = (adjacentEdgeKeys, edges = plan.edges) => project({ ...plan, edges }, { adjacentEdgeKeys }).intentKey;
  assert.equal(intent(['a->b']), 'upload');
  // Genuinely custom plans stay Custom: a pair left out, a pair that is not
  // adjacent, a LOSAT pair (Selected topology limits the LOSATP modes), or an
  // omitted (dormant) pair standing in for an adjacent one.
  assert.equal(intent(['a->b', 'b->c']), 'custom');
  assert.equal(intent(['b->c']), 'custom');
  assert.equal(intent(['a->b'], [{ ...pair, source: 'losat' }]), 'custom');
  assert.equal(intent(['a->b'], [{ ...pair, included: false }]), 'custom');
  assert.equal(intent([]), 'custom');
});
