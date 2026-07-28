import assert from 'node:assert/strict';
import { cp, mkdtemp, writeFile } from 'node:fs/promises';
import { tmpdir } from 'node:os';
import { join } from 'node:path';
import { pathToFileURL } from 'node:url';

const repoRoot = process.cwd();
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
  DEFAULT_LINEAR_BLAST_FILTERS
} = await import(pathToFileURL(join(tempDir, 'js', 'app', 'cli-args.js')));

assert.equal(MODE_PROFILE_VERSION, 1);
assert.deepEqual(
  [...MODE_DEFAULT_FEATURE_TYPES],
  ['CDS', 'rRNA', 'tRNA', 'tmRNA', 'ncRNA', 'misc_RNA', 'repeat_region']
);
assert.deepEqual(comparisonFiltersForMode('circular'), {
  bitscore: 50,
  evalue: '1e-5',
  identity: 70,
  alignment_length: 0
});
assert.deepEqual(comparisonStateForMode('linear'), {
  min_bitscore: 50,
  evalue: '1e-2',
  identity: 0,
  alignment_length: 0
});
assert.deepEqual(
  DEFAULT_CIRCULAR_CONSERVATION_BLAST_FILTERS,
  comparisonFiltersForMode('circular')
);
assert.deepEqual(DEFAULT_LINEAR_BLAST_FILTERS, comparisonFiltersForMode('linear'));
assert.equal(comparisonProfileDefault('circular', 'identity'), 70);
assert.deepEqual(trackDefaultsForMode('circular'), { gc: true, skew: true });
assert.deepEqual(trackDefaultsForMode('linear'), { gc: false, skew: false });
assert.equal(modeProfile('linear').linearAxisColor, 'lightgray');
assert.equal(modeProfile('linear').linearRulerAxisColor, 'dimgray');

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

assert.equal(effectiveLinearAxisColor(), 'lightgray');
assert.equal(effectiveLinearAxisColor({ rulerOnAxis: true }), 'dimgray');
assert.equal(
  effectiveLinearAxisColor({
    axisColor: 'lightgray',
    rulerOnAxis: true,
    managed: true
  }),
  'dimgray'
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
const { createDefaultAdv, createDefaultForm } = await import(
  pathToFileURL(join(tempDir, 'js', 'state.js'))
);
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
