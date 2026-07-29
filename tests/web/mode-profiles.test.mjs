import assert from 'node:assert/strict';
import { cp, mkdtemp, readFile, writeFile } from 'node:fs/promises';
import { tmpdir } from 'node:os';
import { join } from 'node:path';
import { pathToFileURL } from 'node:url';

const repoRoot = process.cwd();
const semanticParity = JSON.parse(await readFile(
  join(repoRoot, 'tests', 'fixtures', 'mode_semantic_parity.json'),
  'utf8'
));
const expectedModes = semanticParity.modes;
const tempDir = await mkdtemp(join(tmpdir(), 'gbdraw-mode-profiles-'));
await cp(join(repoRoot, 'gbdraw', 'web', 'js'), join(tempDir, 'js'), { recursive: true });
await writeFile(join(tempDir, 'package.json'), '{"type":"module"}', 'utf8');

const {
  MODE_DEFAULT_FEATURE_TYPES,
  MODE_PROFILE_VERSION,
  comparisonFiltersForMode,
  comparisonProfileDefault,
  comparisonStateForMode,
  createModeProfileStateManager,
  effectiveLinearAxisColor,
  managedAdvStateForMode,
  modeProfile,
  trackDefaultsForMode
} = await import(pathToFileURL(join(tempDir, 'js', 'mode-profiles.js')));
const {
  DEFAULT_CIRCULAR_CONSERVATION_BLAST_FILTERS,
  DEFAULT_LINEAR_BLAST_FILTERS,
  buildModeBlastFilterArgs,
  buildModeTrackVisibilityArgs
} = await import(pathToFileURL(join(tempDir, 'js', 'app', 'cli-args.js')));

const normalizedComparison = (filters) => ({
  evalue: Number(filters.evalue),
  bitscore: Number(filters.bitscore),
  identity: Number(filters.identity),
  alignmentLength: Number(filters.alignment_length)
});

const differentLeafPaths = (left, right, prefix = '') => {
  const paths = new Set();
  const keys = new Set([...Object.keys(left), ...Object.keys(right)]);
  keys.forEach((key) => {
    const path = prefix ? `${prefix}.${key}` : key;
    const leftValue = left[key];
    const rightValue = right[key];
    if (
      leftValue && rightValue &&
      typeof leftValue === 'object' && typeof rightValue === 'object' &&
      !Array.isArray(leftValue) && !Array.isArray(rightValue)
    ) {
      differentLeafPaths(leftValue, rightValue, path).forEach((entry) => paths.add(entry));
    } else if (!Object.is(leftValue, rightValue)) {
      paths.add(path);
    }
  });
  return paths;
};

assert.equal(semanticParity.schema, 1);
assert.equal(MODE_PROFILE_VERSION, semanticParity.profileVersion);
assert.deepEqual(
  [...MODE_DEFAULT_FEATURE_TYPES],
  semanticParity.featureTypes
);
assert.deepEqual(
  normalizedComparison(comparisonFiltersForMode('circular')),
  expectedModes.circular.comparison
);
assert.deepEqual(
  normalizedComparison({
    ...comparisonStateForMode('linear'),
    bitscore: comparisonStateForMode('linear').min_bitscore
  }),
  expectedModes.linear.comparison
);
assert.deepEqual(
  DEFAULT_CIRCULAR_CONSERVATION_BLAST_FILTERS,
  comparisonFiltersForMode('circular')
);
assert.deepEqual(DEFAULT_LINEAR_BLAST_FILTERS, comparisonFiltersForMode('linear'));
assert.equal(
  comparisonProfileDefault('circular', 'identity'),
  expectedModes.circular.comparison.identity
);
assert.deepEqual(trackDefaultsForMode('circular'), expectedModes.circular.tracks);
assert.deepEqual(trackDefaultsForMode('linear'), expectedModes.linear.tracks);
assert.equal(
  modeProfile('linear').linearAxisColor,
  expectedModes.linear.linearAxisColor
);
assert.equal(
  modeProfile('circular').linearRulerAxisColor,
  expectedModes.circular.linearRulerAxisColor
);
assert.equal(
  modeProfile('linear').linearRulerAxisColor,
  expectedModes.linear.linearRulerAxisColor
);
assert.deepEqual(
  [...differentLeafPaths(expectedModes.circular, expectedModes.linear)].sort(),
  semanticParity.reviewedModeDifferences.map(({ path }) => path).sort()
);
assert.ok(semanticParity.reviewedModeDifferences.every(({ reason }) => reason.trim()));

const managedState = managedAdvStateForMode('circular');
const manager = createModeProfileStateManager('circular', managedState);
manager.transition(managedState, 'circular', 'linear');
assert.deepEqual(managedState, managedAdvStateForMode('linear'));
managedState.identity = 70;
managedState.alignment_length = '';
manager.transition(managedState, 'linear', 'circular');
assert.deepEqual(managedState, managedAdvStateForMode('circular'));
manager.transition(managedState, 'circular', 'linear');
assert.equal(managedState.identity, 70);
assert.equal(
  managedState.alignment_length,
  '',
  'a blank customization must not be treated as the numeric default zero'
);
assert.equal(managedState.axis_stroke_color, 'lightgray');

const importedLinearState = {
  ...managedAdvStateForMode('linear'),
  identity: 91,
  axis_stroke_color: 'gray'
};
const importedManager = createModeProfileStateManager(
  'circular',
  managedAdvStateForMode('circular')
);
importedManager.invalidate('linear');
importedManager.transition(importedLinearState, 'linear', 'circular');
importedManager.transition(importedLinearState, 'circular', 'linear');
assert.equal(importedLinearState.identity, 91);
assert.equal(importedLinearState.axis_stroke_color, 'gray');

assert.equal(
  effectiveLinearAxisColor(),
  expectedModes.linear.linearAxisColor
);
assert.equal(
  effectiveLinearAxisColor({ rulerOnAxis: true }),
  expectedModes.linear.linearRulerAxisColor
);
assert.equal(
  effectiveLinearAxisColor({
    axisColor: expectedModes.linear.linearAxisColor,
    rulerOnAxis: true,
    managed: true
  }),
  expectedModes.linear.linearRulerAxisColor
);
assert.equal(
  effectiveLinearAxisColor({ axisColor: '#123456', rulerOnAxis: true }),
  '#123456'
);

assert.throws(() => comparisonStateForMode('radial'), /Unsupported diagram mode/);
assert.throws(
  () => comparisonProfileDefault('linear', 'unknown'),
  /Unsupported comparison field/
);

globalThis.window = {
  Vue: {
    ref: (value) => ({ value }),
    reactive: (value) => value,
    computed: (getter) => ({
      get value() {
        return getter();
      }
    })
  },
  DOMPurify: { sanitize: (value) => value }
};
const { createDefaultAdv, createDefaultForm, state } = await import(
  pathToFileURL(join(tempDir, 'js', 'state.js'))
);
assert.equal(Object.keys(state.form).includes('legend'), false);
assert.equal(Object.keys(state.adv).includes('plot_title_position'), false);
state.form.legend = 'right';
state.adv.plot_title_position = 'top';
assert.deepEqual(state.layoutPreferences.circular.single, {
  legend: 'right',
  plotTitlePosition: 'top'
});
state.mode.value = 'linear';
state.form.legend = 'left';
state.adv.plot_title_position = 'center';
assert.deepEqual(state.layoutPreferences.linear, {
  legend: 'left',
  plotTitlePosition: 'center'
});
state.mode.value = 'circular';
const formDefaults = createDefaultForm();
assert.deepEqual(
  {
    suppress_gc: formDefaults.suppress_gc,
    suppress_skew: formDefaults.suppress_skew,
    show_gc: formDefaults.show_gc,
    show_skew: formDefaults.show_skew
  },
  {
    suppress_gc: false,
    suppress_skew: false,
    show_gc: false,
    show_skew: false
  }
);
const circularAdv = createDefaultAdv();
assert.deepEqual(
  {
    evalue: circularAdv.evalue,
    identity: circularAdv.identity,
    features: circularAdv.features
  },
  {
    evalue: '1e-5',
    identity: 70,
    features: ['CDS', 'rRNA', 'tRNA', 'tmRNA', 'ncRNA', 'misc_RNA', 'repeat_region']
  }
);
const linearAdv = createDefaultAdv('linear');
assert.deepEqual(
  {
    evalue: linearAdv.evalue,
    identity: linearAdv.identity,
    axis_stroke_color: linearAdv.axis_stroke_color
  },
  { evalue: '1e-2', identity: 0, axis_stroke_color: 'lightgray' }
);

const { buildCanonicalSessionRequest } = await import(
  pathToFileURL(join(tempDir, 'js', 'services', 'session-request.js'))
);
const genbankText = `LOCUS       MODEPARITY                  4 bp    DNA     linear   UNK 01-JAN-1980
DEFINITION  Mode semantic parity fixture.
ACCESSION   MODEPARITY
VERSION     MODEPARITY
KEYWORDS    .
SOURCE      .
  ORGANISM  .
            .
FEATURES             Location/Qualifiers
ORIGIN
        1 atgc
//
`;
const genbank = {
  name: 'mode-parity.gb',
  type: 'text/plain',
  size: new TextEncoder().encode(genbankText).byteLength,
  lastModified: 0,
  data: btoa(genbankText)
};

for (const modeName of ['circular', 'linear']) {
  state.mode.value = modeName;
  Object.assign(state.form, createDefaultForm());
  Object.assign(state.adv, createDefaultAdv(modeName));
  state.modeProfileStateManager.reset(modeName, state.adv);
  state.cInputType.value = 'gb';
  state.lInputType.value = 'gb';
  state.circularRecordList.value = [];
  state.linearRecordLayoutEnabled.value = false;

  const freshCliArgs = [
    ...buildModeTrackVisibilityArgs(modeName, state.form),
    ...buildModeBlastFilterArgs(modeName, {
      bitscore: state.adv.min_bitscore,
      evalue: state.adv.evalue,
      identity: state.adv.identity,
      alignment_length: state.adv.alignment_length
    })
  ];
  assert.deepEqual(
    freshCliArgs,
    [],
    `Fresh ${modeName} Web state must rely on the reviewed CLI defaults`
  );

  const filesData = modeName === 'circular'
    ? { c_gb: genbank, linearSeqs: [] }
    : {
        linearSeqs: [{
          uid: 'record-1',
          gb: genbank,
          region_record_id: '',
          region_start: null,
          region_end: null,
          region_reverse: false
        }],
        linearComparisons: []
      };
  const canonical = buildCanonicalSessionRequest({ state, filesData });
  const options = canonical.renderRequest.diagramOptions;
  const expected = expectedModes[modeName];
  assert.equal(canonical.renderRequest.mode, modeName);
  assert.deepEqual({
    evalue: options.evalue,
    bitscore: options.bitscore,
    identity: options.identity,
    alignmentLength: options.alignmentLength
  }, expected.comparison);
  assert.deepEqual(options.selectedFeaturesSet, semanticParity.featureTypes);
  assert.equal(options.configOverrides['canvas.show_gc'], expected.tracks.gc);
  assert.equal(options.configOverrides['canvas.show_skew'], expected.tracks.skew);
  if (modeName === 'linear') {
    assert.equal(
      options.configOverrides['objects.axis.linear.stroke_color'],
      expected.linearAxisColor
    );
  } else {
    assert.equal(
      Object.hasOwn(options.configOverrides, 'objects.axis.linear.stroke_color'),
      false
    );
  }
}

assert.deepEqual(
  buildModeTrackVisibilityArgs('circular', {
    ...createDefaultForm(),
    suppress_gc: true
  }),
  ['--no-gc']
);
assert.deepEqual(
  buildModeTrackVisibilityArgs('linear', {
    ...createDefaultForm(),
    show_skew: true
  }),
  ['--skew']
);
