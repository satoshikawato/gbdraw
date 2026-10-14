#!/usr/bin/env node
// Writes the vectors of the Web Load value migrations of a Session 27-44 draft
// config, which the Python twins in gbdraw/session_io.py read
// (tests/test_draft_value_migrations.py). The expected outputs come from the
// Web functions themselves, so the Web is the reference:
//   - option values: the option-value steps of migratePersistedWebOptionValues;
//   - Linear track slots: migrateImportedLinearTrackSlots;
//   - Session 27-33 shapes: withoutLegacyNullCircularSlotSpacing and
//     migrateLegacyFeatureRenderingConfig.
// services/config.js does not load outside a browser, so its module-private
// functions are read from its source and run with their dependencies.
// Run with --check to fail when a written file is stale.
import fs from 'node:fs';
import path from 'node:path';
import { gunzipSync } from 'node:zlib';
import { fileURLToPath } from 'node:url';
import {
  migratePersistedCircularMultiRecordSizeMode,
  migratePersistedLinearLabelPlacement,
  migratePersistedLinearTrackLayout
} from '../gbdraw/web/js/services/current-option-values.js';
import {
  LEGACY_LINEAR_TRACK_SLOT_SCHEMA_VERSION,
  LINEAR_TRACK_SLOT_SCHEMA_VERSION,
  migrateLinearTrackSlotsToCurrentSchema
} from '../gbdraw/web/js/services/linear-track-slot-model.js';
import { CIRCULAR_TRACK_SLOT_SCHEMA_VERSION } from '../gbdraw/web/js/services/session-active-config-contract.js';

const ROOT = path.resolve(path.dirname(fileURLToPath(import.meta.url)), '..');
const CONFIG_SOURCE = fs.readFileSync(path.join(ROOT, 'gbdraw/web/js/services/config.js'), 'utf8');

const definition = (name) => {
  const start = CONFIG_SOURCE.indexOf(`\nconst ${name} = `) + 1;
  if (start === 0) throw new Error(`services/config.js has no ${name}.`);
  // A one-line definition ends on its line; a block ends at its first closing line.
  const firstLine = CONFIG_SOURCE.slice(start, CONFIG_SOURCE.indexOf('\n', start));
  if (firstLine.endsWith(';')) return firstLine;
  return CONFIG_SOURCE.slice(start, CONFIG_SOURCE.indexOf('\n};\n', start) + 3);
};
const constant = (name) => {
  const match = new RegExp(`\\nconst ${name} = (\\d+);`).exec(CONFIG_SOURCE);
  if (!match) throw new Error(`services/config.js has no ${name}.`);
  return Number(match[1]);
};
const load = (name, dependencies) => new Function(
  ...Object.keys(dependencies), `${definition(name)}\nreturn ${name};`
)(...Object.values(dependencies));

const isPlainObject = load('isPlainObject', {});
const identity = (value) => value;
const migratePersistedWebOptionValues = load('migratePersistedWebOptionValues', {
  isPlainObject,
  migratePersistedLinearTrackLayout,
  migratePersistedLinearLabelPlacement,
  migratePersistedCircularMultiRecordSizeMode,
  // These steps have their own Python twins (field names; the split's profile
  // rule), so the vectors hold the option-value steps only.
  migratePersistedWebStateFieldNames: identity,
  withHistoricalPairwiseMatchStyleFallback: identity,
  withCurrentLinearLabelVisibility: identity
});
const migrateImportedLinearTrackSlots = load('migrateImportedLinearTrackSlots', {
  LEGACY_LINEAR_TRACK_SLOT_SCHEMA_VERSION,
  LINEAR_TRACK_SLOT_SCHEMA_VERSION,
  LEGACY_LINEAR_TRACK_SLOT_SESSION_VERSION: constant('LEGACY_LINEAR_TRACK_SLOT_SESSION_VERSION'),
  migrateLinearTrackSlotsToCurrentSchema
});
const withoutLegacyNullCircularSlotSpacing = load('withoutLegacyNullCircularSlotSpacing', {
  isPlainObject, CIRCULAR_TRACK_SLOT_SCHEMA_VERSION
});
const migrateLegacyFeatureRenderingConfig = load('migrateLegacyFeatureRenderingConfig', { isPlainObject });

const fixtureConfig = (fixture) => {
  const bytes = fs.readFileSync(path.join(ROOT, fixture));
  return JSON.parse((bytes[0] === 0x1f && bytes[1] === 0x8b ? gunzipSync(bytes) : bytes).toString('utf8')).config;
};
const clone = (value) => JSON.parse(JSON.stringify(value));
const run = (step) => (testCase) => {
  const config = testCase.fixture ? fixtureConfig(testCase.fixture) : clone(testCase.input.config);
  try {
    return { ...testCase, expected: { config: step(config, testCase) } };
  } catch (error) {
    return { ...testCase, error: error.message };
  }
};

const BGC_V30 = 'tests/fixtures/sessions/BGC0000708-BGC0000713.v30.gbdraw-session.json.gz';
const BGC_V33 = 'tests/fixtures/sessions/BGC0000708-BGC0000713.schema-v2.gbdraw-session.json.gz';
const EDITS_V33 = 'tests/fixtures/sessions/feature-edits-circular.v33.gbdraw-session.json.gz';
const CLI_V30 = 'tests/fixtures/sessions/cli-linear-protein.v30.gbdraw-session.json.gz';

const optionCase = (name, form, adv) => ({ name, input: { config: { ...(form ? { form } : {}), ...(adv ? { adv } : {}) } } });
const OPTION_CASES = [
  optionCase('retired Linear track layouts become above and below', { linear_track_layout: 'spreadout', prefix: 'p' }, null),
  optionCase('tuckin becomes below', { linear_track_layout: 'tuckin' }, null),
  optionCase('a retired value is read trimmed and case-blind', { linear_track_layout: ' SpreadOut ' }, { label_placement: ' On_Feature ' }),
  optionCase('on_feature becomes above_feature and sqrt becomes auto', null, { label_placement: 'on_feature', multi_record_size_mode: 'sqrt', nt: 'GC' }),
  optionCase('current values stay', { linear_track_layout: 'middle' }, { label_placement: 'above_feature', multi_record_size_mode: 'equal' }),
  optionCase('an empty or null value reads as the default', { linear_track_layout: '' }, { label_placement: null, multi_record_size_mode: '' }),
  optionCase('an unknown Linear track layout is refused', { linear_track_layout: 'diagonal' }, null),
  optionCase('an unknown label placement is refused', null, { label_placement: 'inside' }),
  optionCase('an unknown multi-record size mode is refused', null, { multi_record_size_mode: 'sqrt2' }),
  optionCase('a draft without these fields is unchanged', { prefix: 'p' }, { nt: 'AT' }),
  { name: 'CLI-written Session 30 (release 0.13.0 era)', fixture: CLI_V30 }
];

const slotCase = (name, sessionVersion, adv) => ({ name, sessionVersion, input: { config: { adv } } });
const featuresSlot = { id: 'features', renderer: 'features', height: 30, spacing: 4, params: { a: 1 } };
const SLOT_CASES = [
  { name: 'Session 30 schema-1 Linear slots (real file)', sessionVersion: 30, fixture: BGC_V30 },
  { name: 'Session 33 storing schema 1 (real file)', sessionVersion: 33, fixture: BGC_V33 },
  { name: 'Session 33 storing schema 2 (real file)', sessionVersion: 33, fixture: EDITS_V33 },
  slotCase('Session 32 reads schema-1 meaning whatever it stored', 32, {
    linear_track_slots_schema_version: 2, linear_track_slots: [featuresSlot]
  }),
  slotCase('Session 40 storing schema 1 is migrated', 40, {
    linear_track_slots_schema_version: 1, linear_track_slots: [featuresSlot]
  }),
  slotCase('Session 40 storing schema 2 keeps geometry', 40, {
    linear_track_slots_schema_version: 2, linear_track_slots: [featuresSlot]
  }),
  slotCase('only a features slot loses its geometry; renderer aliases and defaults', 30, {
    linear_track_slots_schema_version: 1,
    linear_track_slots: [
      { id: 'skew', renderer: 'GC_skew', height: 20, spacing: 2 },
      { id: 'unknown', renderer: 'mystery', height: 20, spacing: 2 },
      { id: 'blank', renderer: '', height: 20 },
      { id: 'bad-params', renderer: 'depth', params: [1, 2] },
      'not-a-slot',
      null
    ]
  }),
  slotCase('Session 40 without a stored schema is refused', 40, { linear_track_slots: [] }),
  slotCase('an unknown stored schema is refused', 30, { linear_track_slots_schema_version: 3, linear_track_slots: [] }),
  slotCase('a non-integer stored schema is refused', 33, { linear_track_slots_schema_version: '2', linear_track_slots: [] }),
  slotCase('slots that are not a list are refused', 30, { linear_track_slots: { id: 'x' } }),
  slotCase('a draft without Linear slots is unchanged', 30, { nt: 'GC' })
];

const shapeCase = (name, step, adv) => ({ name, step, input: { config: { adv } } });
const circularSlots = [{ id: 'a', spacing: null, width: 1 }, { id: 'b', spacing: 0 }, { id: 'c' }, 'x'];
const SHAPE_CASES = [
  { name: 'Session 30 Circular slots drop spacing null (real file)', step: 'circular-slot-spacing', fixture: BGC_V30 },
  { name: 'Session 33 Circular slots drop spacing null (real file)', step: 'circular-slot-spacing', fixture: EDITS_V33 },
  shapeCase('only a null spacing goes, with Custom Track Slots off', 'circular-slot-spacing', {
    circular_track_slots_enabled: false, circular_track_slots_schema_version: 4, circular_track_slots: circularSlots
  }),
  shapeCase('Custom Track Slots on keep their spacing', 'circular-slot-spacing', {
    circular_track_slots_enabled: true, circular_track_slots_schema_version: 4, circular_track_slots: circularSlots
  }),
  shapeCase('another slot schema is left alone', 'circular-slot-spacing', {
    circular_track_slots_schema_version: 3, circular_track_slots: circularSlots
  }),
  { name: 'Session 30 repeat regions are drawn as rectangles (real file)', step: 'repeat-region-shape', fixture: BGC_V30 },
  { name: 'Session 33 repeat regions are drawn as rectangles (real file)', step: 'repeat-region-shape', fixture: EDITS_V33 },
  shapeCase('a feature list without repeat regions is unchanged', 'repeat-region-shape', {
    features: ['CDS'], feature_shapes: { CDS: 'arrow' }
  }),
  shapeCase('no feature list: the repeat-region shape is added', 'repeat-region-shape', { feature_shapes: { CDS: 'arrow' } }),
  shapeCase('a saved repeat-region shape stays', 'repeat-region-shape', {
    features: ['repeat_region'], feature_shapes: { repeat_region: 'underlay' }
  }),
  shapeCase('a non-object shape map becomes the repeat-region shape', 'repeat-region-shape', {
    features: ['repeat_region'], feature_shapes: ['arrow']
  })
];

const SHAPE_STEPS = {
  'circular-slot-spacing': withoutLegacyNullCircularSlotSpacing,
  'repeat-region-shape': (config) => migrateLegacyFeatureRenderingConfig(config, true)
};

const FILES = {
  'tests/fixtures/draft-option-value-migration-vectors.json': {
    description: 'Option-value steps of migratePersistedWebOptionValues (services/config.js), Sessions 27-39; Python twin migrate_persisted_web_option_values.',
    cases: OPTION_CASES.map(run((config) => migratePersistedWebOptionValues(config)))
  },
  'tests/fixtures/linear-track-slot-migration-vectors.json': {
    description: 'migrateImportedLinearTrackSlots (services/config.js) with migrateLinearTrackSlotsToCurrentSchema; Python twin migrate_imported_linear_track_slots.',
    cases: SLOT_CASES.map(run((config, testCase) => migrateImportedLinearTrackSlots(config, testCase.sessionVersion)))
  },
  'tests/fixtures/session-33-draft-shape-migration-vectors.json': {
    description: 'withoutLegacyNullCircularSlotSpacing and migrateLegacyFeatureRenderingConfig (services/config.js), Sessions 27-33.',
    cases: SHAPE_CASES.map(run((config, testCase) => SHAPE_STEPS[testCase.step](config)))
  }
};

const render = ({ description, cases }) => `${JSON.stringify({ schemaVersion: 1, description, cases }, null, 2)}\n`;
const check = process.argv.includes('--check');
let stale = false;
for (const [file, content] of Object.entries(FILES)) {
  const target = path.join(ROOT, file);
  const text = render(content);
  if (check) {
    if (!fs.existsSync(target) || fs.readFileSync(target, 'utf8') !== text) {
      console.log(`${file} is stale. Run: node tools/generate_draft_value_migration_vectors.mjs`);
      stale = true;
    }
  } else {
    fs.writeFileSync(target, text);
  }
}
process.exit(stale ? 1 : 0);
