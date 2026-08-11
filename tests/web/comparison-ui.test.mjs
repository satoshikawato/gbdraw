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
    join(repoRoot, 'gbdraw', 'web', 'js', 'app', 'linear-comparisons.js'),
    join(appRoot, 'linear-comparisons.js')
  ),
  cp(
    join(repoRoot, 'gbdraw', 'web', 'js', 'mode-profiles.js'),
    join(tempRoot, 'js', 'mode-profiles.js')
  ),
  cp(
    join(repoRoot, 'gbdraw', 'web', 'js', 'mode-profiles.generated.js'),
    join(tempRoot, 'js', 'mode-profiles.generated.js')
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
} = await import(pathToFileURL(join(appRoot, 'linear-comparisons.js')));
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
  assert.equal(projection.currentStatusLabel, 'Current: No comparison');
  assert.equal(projection.summaryText, 'No comparison');
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
  assert.match(projection.summaryText, /^LOSATN · 2 adjacent pairs · E-value <= /);

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
  assert.match(oneRecord.summaryText, /^LOSATN · 0 adjacent pairs/);
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
  assert.equal(
    projection.currentStatusLabel,
    'Current: Upload BLAST TSV for all adjacent pairs'
  );
  assert.match(projection.summaryText, /^Upload BLAST TSV · 2 adjacent pairs/);
  assert(projection.sectionKeys.settings.includes('upload-readiness'));
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
  assert.equal(selectedUpload.currentStatusLabel, 'Current: Selected pairs (1; 1 upload)');
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
  assert.equal(
    mixed.currentStatusLabel,
    'Current: Selected pairs (2; 1 LOSAT, 1 upload)'
  );
  assert(mixed.sectionKeys.settings.includes('losat-mode'));
  assert(mixed.sectionKeys.settings.includes('upload-readiness'));
});

test('LOSAT and LOSATP modes project only their active Settings and Advanced controls', () => {
  const plan = { mode: 'adjacent', defaultSource: 'losat', edges: [] };
  const cases = [
    {
      program: 'blastn',
      blastpMode: 'orthogroup',
      losatKey: 'blastn',
      losatpKey: 'orthogroup',
      settings: ['losat-mode', 'blastn-task', 'result-filters', 'comparison-appearance'],
      absent: [
        'losatp-mode', 'record-genetic-codes', 'blastp-max-hits',
        'similarity-grouping', 'collinear-primary'
      ]
    },
    {
      program: 'tblastx',
      blastpMode: 'orthogroup',
      losatKey: 'tblastx',
      losatpKey: 'orthogroup',
      settings: ['losat-mode', 'record-genetic-codes', 'result-filters', 'comparison-appearance'],
      absent: [
        'losatp-mode', 'blastn-task', 'blastp-max-hits',
        'similarity-grouping', 'collinear-primary'
      ]
    },
    {
      program: 'blastp',
      blastpMode: 'pairwise',
      losatKey: 'blastp',
      losatpKey: 'pairwise',
      settings: [
        'losat-mode', 'losatp-mode', 'blastp-max-hits',
        'result-filters', 'comparison-appearance'
      ],
      absent: ['blastn-task', 'record-genetic-codes', 'similarity-grouping', 'collinear-primary']
    },
    {
      program: 'blastp',
      blastpMode: 'orthogroup',
      losatKey: 'blastp',
      losatpKey: 'orthogroup',
      settings: ['losat-mode', 'losatp-mode', 'similarity-grouping', 'result-filters'],
      absent: ['blastn-task', 'record-genetic-codes', 'blastp-max-hits', 'collinear-primary', 'comparison-appearance']
    },
    {
      program: 'blastp',
      blastpMode: 'collinear',
      losatKey: 'blastp',
      losatpKey: 'collinear',
      settings: ['losat-mode', 'losatp-mode', 'collinear-primary', 'result-filters'],
      absent: [
        'blastn-task', 'record-genetic-codes', 'blastp-max-hits',
        'similarity-grouping', 'comparison-appearance'
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
  assert.match(
    similarityGroups.summaryText,
    /^LOSATP · Similarity groups · 2 adjacent pairs/
  );
  const collinear = project(plan, {
    losatProgram: 'blastp',
    blastpMode: 'collinear'
  });
  assert.match(
    collinear.summaryText,
    /^LOSATP · Collinear blocks · 2 adjacent pairs/
  );
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
  assert(!grouping.sectionKeys.settings.includes('similarity-grouping'));
  assert.equal(grouping.errorDisclosureKey, 'settings');

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
  assert.equal(dormant.retainedRawResultDraftCount, 1);
  assert.deepEqual(dormant.sectionKeys.selectedPairs, ['pair-editor', 'retained-drafts']);
  assert(dormant.sectionKeys.advanced.includes('raw-results'));

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
  assert.equal(activeUpload.retainedRawResultDraftCount, 0);

  activeUploadPlan.edges[0].losatFilename = 'alternate.raw.tsv';
  activeUploadPlan.edges[0].losatFilenameActive = false;
  const alternatePayload = project(activeUploadPlan);
  assert.equal(alternatePayload.activePairCount, 1);
  assert.equal(alternatePayload.retainedDormantDraftCount, 1);
});

test('structured resolver issue codes project badges and disclosure focus targets', () => {
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

  assert.deepEqual(projection.errorBadge, {
    count: 7,
    label: '7 comparison issues'
  });
  assert.equal(projection.errorTargets.length, 7);
  assert.deepEqual(projection.errorDisclosureKeys, ['selected-pairs', 'settings']);
  assert.equal(projection.errorDisclosureKey, 'selected-pairs');
  assert.deepEqual(
    projection.errorTargets.find((target) => target.issueCode === 'missing-upload'),
    {
      issueCode: 'missing-upload',
      disclosureKey: 'selected-pairs',
      focusTargetKey: 'pair-upload',
      edgeId: 'problem',
      edgeKey: 'a->b'
    }
  );
  assert.deepEqual(
    projection.errorTargets.find(
      (target) => target.issueCode === 'selected-losat-requires-pairwise'
    ),
    {
      issueCode: 'selected-losat-requires-pairwise',
      disclosureKey: 'settings',
      focusTargetKey: 'losatp-mode',
      edgeId: 'problem',
      edgeKey: 'a->b'
    }
  );
  ['duplicate', 'missing-uid', 'self', 'same-row', 'non-adjacent'].forEach((code) => {
    const target = projection.errorTargets.find((entry) => entry.issueCode === code);
    assert.equal(target.disclosureKey, 'selected-pairs');
    assert.equal(target.focusTargetKey, 'pair-row');
  });
  assert(!projection.errorTargets.some((target) => target.issueCode === 'unowned-error'));
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
  assert(Object.isFrozen(projection.errorTargets));
});
