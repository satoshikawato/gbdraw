#!/usr/bin/env node
// Writes each diagram mode's default settings, as the Web builds them, to
// gbdraw/data/mode-slice-defaults.json, which Python reads to decide whether
// a Session 46-50 used the mode it did not show (rule d of the move to
// drawings, PD-OI-090; gbdraw.session_drawings), and the same data to
// gbdraw/web/js/mode-slice-defaults.generated.js, which the Web's move
// (services/session-drawings-upgrade.js) reads without building a drawing.
//
// The Web owns the defaults (D-15-86): a new drawing as the app builds one
// (`withNewDrawingDefaults` over `createDrawingState(id, mode)`, with the
// palette file's default palette, services/reset.js). Its slice is read
// as Save writes it (`buildModeSliceData`), one value per
// `MODE_SCOPED_SETTINGS` row the mode holds, keyed by the row's identity
// (domain, path).
//
// A value the Web fills from run-time identities cannot be a constant
// default. `NORMALIZE` names each such row field and how a reader normalizes
// it before comparing; Python (`gbdraw.session_drawings._normalized`) and the
// Web (`withoutUnarrangedEntries`, services/linear-record-layout.js) evaluate
// the rules of the file. The only rule: a Linear record
// layout that arranges nothing (each sequence row in its own layout row, in
// order, as `reconcileLinearRecordLayout` writes it from no saved entries, for
// example the empty sequence row of a new project with its run-time uid) has
// no rows. The generator checks that rule against `reconcileLinearRecordLayout`.
// Run with --check to fail when a written file is stale.
import fs from 'node:fs';
import path from 'node:path';
import { fileURLToPath } from 'node:url';

globalThis.window = {
  Vue: {
    ref: (value) => ({ value }),
    reactive: (value) => value,
    computed: (getter) => ({ get value() { return getter(); } }),
    nextTick: async () => {}
  },
  DOMPurify: { sanitize: (value) => value }
};
globalThis.document = {};

const ROOT = path.resolve(path.dirname(fileURLToPath(import.meta.url)), '..');
const TARGET = 'gbdraw/data/mode-slice-defaults.json';
const MODULE_TARGET = 'gbdraw/web/js/mode-slice-defaults.generated.js';
const { MODE_SCOPED_SETTINGS } = await import('../gbdraw/web/js/mode-scoped-settings.generated.js');
const { buildModeSliceData } = await import('../gbdraw/web/js/services/config.js');
const { createDrawingState } = await import('../gbdraw/web/js/state.js');
const { withNewDrawingDefaults } = await import('../gbdraw/web/js/services/reset.js');
const { normalizePaletteDefinitions } = await import('../gbdraw/web/js/utils/color-utils.js');
const { reconcileLinearRecordLayout, withoutUnarrangedEntries } = await import('../gbdraw/web/js/services/linear-record-layout.js');

/** @import { UnarrangedFieldRule } from '../gbdraw/web/js/services/linear-record-layout.js' */
/** @type {UnarrangedFieldRule[]} */
const NORMALIZE = [
  {
    domain: 'config',
    path: 'linearRecordLayout',
    field: 'rows',
    emptyWhenUnarranged: { entryKeys: ['row', 'uid'], positionKey: 'row' },
    reason: 'The Web fills the layout with one entry per sequence row, keyed by a uid it makes at run time.'
  }
];

// A row's default that a released writer stored and the app has since changed.
// A saved value equal to one is a default, not a setting: rule d must not read
// it as a use of the mode (a Session 46-50 of one mode stays one drawing).
const EARLIER_DEFAULTS = [
  {
    domain: 'config.losat',
    path: 'blastp',
    reason: 'gbdraw 0.14.3 and earlier saved these Collinear defaults (merge either, hit limits 5).',
    value: {"mode": "orthogroup", "hitLimitsByMode": {"pairwise": {"candidateLimit": null, "orthogroupMemberMaxHits": null}, "orthogroup": {"candidateLimit": null, "orthogroupMemberMaxHits": null}, "collinear": {"candidateLimit": 5, "orthogroupMemberMaxHits": 5}}, "maxHits": 5, "candidateLimit": null, "orthogroupMembershipMode": "anchor_core_v1", "orthogroupMemberMaxHits": null, "collinearInferOrthogroups": false, "collinearMinAnchors": 1, "collinearMaxUnitGap": 0, "collinearMaxDiagonalDrift": 0, "collinearMaxConflictsInMergeGap": 1, "collinearMaxParalogLinksPerOrthogroup": 2, "collinearColorMode": "orientation", "collinearUnitMode": "auto", "collinearAnchorMode": "rbh", "collinearMergeOrientation": "either", "collinearSearchScope": "adjacent"}
  },
  {
    domain: 'config.adv',
    path: 'depth_tracks',
    reason: 'Sessions 46 and 47 of main saved the default Depth series with the label Depth.',
    value: [{"label": "Depth", "color": "#4A90E2", "height": null, "large_tick_interval": null, "small_tick_interval": null, "tick_font_size": null}]
  },
  {
    domain: 'config.adv',
    path: 'depth_tracks',
    reason: 'Sessions up to 50 saved the default Depth series with its ticks never set as null; a Session 51 leaves them out.',
    value: [{"label": "", "color": "#4A90E2", "height": null, "large_tick_interval": null, "small_tick_interval": null, "tick_font_size": null}]
  },
  {
    domain: 'config',
    path: 'colors',
    reason: 'gbdraw 0.14.5 and earlier saved the default palette with one color for repeat_region, misc_feature and every unlisted type (OV-3114).',
    value: {"CDS": "#54bcf8", "rRNA": "#71ee7d", "tRNA": "#e8b441", "tmRNA": "#ded44e", "ncRNA": "#c4fac3", "repeat_region": "#d3d3d3", "misc_feature": "#d3d3d3", "default": "#d3d3d3", "skew_high": "#6dded3", "skew_low": "#ad72e3", "gc_content": "#a1a1a1", "pairwise_match": "#d3d3d3", "pairwise_match_min": "#FFE7E7", "pairwise_match_max": "#FF7272", "collinear_block_plus_min": "#f0f1f5", "collinear_block_plus": "#8b9cc1", "collinear_block_minus_min": "#FFE7E7", "collinear_block_minus": "#E15759"}
  },
  {
    domain: 'config',
    path: 'circularConservation',
    reason: 'Sessions up to 50 saved the comparison rings without the Ring colors setting (colorScheme, reassignOnReorder).',
    value: {"enabled": false, "source": "losat", "losat_program": "blastn", "subject_gencode": 1, "reference": "auto", "labels": "", "series": [], "ring_width": null, "ring_gap": null}
  }
];

// The rows of each packaged palette whose color gbdraw 0.14.5 saved and the
// palette has since changed (OV-3114, D-54), as that release's
// gbdraw/data/color_palettes.toml holds them; written from the tag, and
// tests/test_user_colors.py compares them with it. A 0.14.5 Web Session saved
// its palette's rows with the user's colors; Web Load reads a saved row that
// still holds such a value as the palette's own (OV-3259, `savedPaletteRowsOf`).
const EARLIER_PALETTE_ROWS = {
  release: '0.14.5',
  reason: 'gbdraw 0.14.5 and earlier saved each palette\'s misc_feature and default rows, which OV-3114 changed, in config.colors.',
  palettes: {
    "default": {"misc_feature": "#d3d3d3", "default": "#d3d3d3"},
    "mint": {"misc_feature": "#d3d3d3", "default": "#d3d3d3"},
    "orange": {"misc_feature": "#d3d3d3", "default": "#d3d3d3"},
    "soil": {"misc_feature": "#d3d3d3", "default": "#d3d3d3"},
    "spring": {"misc_feature": "#d3d3d3", "default": "#d3d3d3"},
    "ajisai": {"misc_feature": "#d3d3d3", "default": "#d3d3d3"},
    "alpine_retreat": {"misc_feature": "#f0f4f8", "default": "#f0f4f8"},
    "arctic": {"misc_feature": "#d3d3d3", "default": "#d3d3d3"},
    "arctic_dawn": {"misc_feature": "#d3d3d3", "default": "#d3d3d3"},
    "autumn": {"misc_feature": "#d3d3d3", "default": "#d3d3d3"},
    "autumn_harvest": {"misc_feature": "#d3d3d3", "default": "#d3d3d3"},
    "berry_field": {"misc_feature": "#d3d3d3", "default": "#d3d3d3"},
    "citrine_verdure": {"misc_feature": "#f0f8e7", "default": "#f0f8e7"},
    "citrus_grove": {"misc_feature": "#e0e0e0", "default": "#e0e0e0"},
    "deep_saturates": {"misc_feature": "#afafaf", "default": "#afafaf"},
    "desert_mirage": {"misc_feature": "#c0c0c0", "default": "#c0c0c0"},
    "dune": {"misc_feature": "#d3d3d3", "default": "#d3d3d3"},
    "earth": {"misc_feature": "#eaeaea", "default": "#eaeaea"},
    "earth_and_clay": {"misc_feature": "#c5bcb3", "default": "#c5bcb3"},
    "edelweiss": {"misc_feature": "#eaebed", "default": "#eaebed"},
    "emerald_twilight": {"misc_feature": "#e0f2e9", "default": "#e0f2e9"},
    "forest": {"misc_feature": "#d3d3d3", "default": "#d3d3d3"},
    "forest_canopy": {"misc_feature": "#e5e5e5", "default": "#e5e5e5"},
    "forest_morning": {"misc_feature": "#f0f0f0", "default": "#f0f0f0"},
    "fugaku": {"misc_feature": "#e4e4e4", "default": "#e4e4e4"},
    "galaxy": {"misc_feature": "#f0f0f0", "default": "#f0f0f0"},
    "herbal_garden": {"misc_feature": "#dce1de", "default": "#dce1de"},
    "kuroshio": {"misc_feature": "#d3d3d3", "default": "#d3d3d3"},
    "lavender_fields": {"misc_feature": "#d8d8d8", "default": "#d8d8d8"},
    "majestic_sunrise": {"misc_feature": "#f4a460", "default": "#f4a460"},
    "marine": {"misc_feature": "#bbbbbb", "default": "#bbbbbb"},
    "matcha_whispers": {"misc_feature": "#e1e1e1", "default": "#e1e1e1"},
    "metallic_sheen": {"misc_feature": "#cd7f32", "default": "#cd7f32"},
    "metropolitan_blues": {"misc_feature": "#eef0f2", "default": "#eef0f2"},
    "misty_mountains": {"misc_feature": "#d3d3d3", "default": "#d3d3d3"},
    "oceanic_voyage": {"misc_feature": "#e8e8e8", "default": "#e8e8e8"},
    "orchid": {"misc_feature": "#dddddd", "default": "#dddddd"},
    "pbt": {"misc_feature": "#d3d3d3", "default": "#d3d3d3"},
    "pine_reflection": {"misc_feature": "#e3ebe7", "default": "#e3ebe7"},
    "royal_gala": {"misc_feature": "#dda0dd", "default": "#dda0dd"},
    "sakura": {"misc_feature": "#d3d3d3", "default": "#d3d3d3"},
    "savannah": {"misc_feature": "#ededed", "default": "#ededed"},
    "seaside": {"misc_feature": "#e8e8e8", "default": "#e8e8e8"},
    "soft_pastels": {"misc_feature": "#e1bee7", "default": "#e1bee7"},
    "sunset": {"misc_feature": "#f1f1f1", "default": "#f1f1f1"},
    "sunset_desert": {"misc_feature": "#f2f2f2", "default": "#f2f2f2"},
    "tropical": {"misc_feature": "#e0e0e0", "default": "#e0e0e0"},
    "tropical_sunset": {"misc_feature": "#ebf8e1", "default": "#ebf8e1"},
    "urban": {"misc_feature": "#f2e9e4", "default": "#f2e9e4"},
    "loess": {"misc_feature": "#d3d3d3", "default": "#d3d3d3"},
    "vintage_postcard": {"misc_feature": "#fdfbfa", "default": "#fdfbfa"},
    "volcano": {"misc_feature": "#e8e8e8", "default": "#e8e8e8"},
    "zen": {"misc_feature": "#6a6a6a", "default": "#6a6a6a"}
  }
};

/**
 * @param {{ domain: string, path: string }} row
 * @param {unknown} value
 */
const normalized = (row, value) => /** @type {Record<string, unknown[]>} */ (withoutUnarrangedEntries(NORMALIZE, row, value));

// The rule against its Web owner: what the layout writes from no saved
// entries is unarranged; one row holding two records is not.
const layoutRow = { domain: 'config', path: 'linearRecordLayout' };
for (const count of [0, 1, 3]) {
  const rows = reconcileLinearRecordLayout(Array.from({ length: count }, (_, index) => ({ uid: `seq-${index}` })), []);
  if (normalized(layoutRow, { rows }).rows.length) throw new Error(`An unarranged layout of ${count} rows keeps rows.`);
}
if (!normalized(layoutRow, { rows: reconcileLinearRecordLayout([{ uid: 'a' }, { uid: 'b' }], [{ uid: 'b', row: 1 }]) }).rows.length) {
  throw new Error('An arranged layout loses its rows.');
}

/** @typedef {{ domain: string, path: string, modes: string }} RegistryRow */

/**
 * @param {Record<string, unknown>} slice
 * @param {RegistryRow} row
 * @returns {{ value: unknown } | null} Null when the slice leaves the row out.
 */
const rowValue = (slice, row) => {
  /** @type {unknown} */
  let holder = slice;
  for (const part of row.domain.split('.')) {
    holder = holder && typeof holder === 'object' ? /** @type {Record<string, unknown>} */ (holder)[part] : undefined;
  }
  if (!holder || typeof holder !== 'object') throw new Error(`The default slice holds no ${row.domain}.`);
  return Object.hasOwn(holder, row.path) ? { value: /** @type {Record<string, unknown>} */ (holder)[row.path] } : null;
};

const palettes = JSON.parse(fs.readFileSync(path.join(ROOT, 'gbdraw/web/gallery/palettes/palettes.json'), 'utf8'));
const paletteDefinitions = normalizePaletteDefinitions(palettes.palettes || palettes);
if (!Object.keys(paletteDefinitions).length) throw new Error('The palette file is empty.');
// A mode's default is a new drawing of it, built by the Web's one builder.
const paletteState = /** @type {Parameters<typeof withNewDrawingDefaults>[0]} */ (/** @type {unknown} */ ({ paletteDefinitions: { value: paletteDefinitions } }));
/** @param {'circular' | 'linear'} mode */
const defaultDrawing = (mode) => withNewDrawingDefaults(paletteState, createDrawingState(mode, mode));
const drawings = { circular: defaultDrawing('circular'), linear: defaultDrawing('linear') };

/** @param {'circular' | 'linear'} mode */
const modeDefaults = (mode) => {
  const slice = JSON.parse(JSON.stringify(buildModeSliceData(drawings[mode], mode)));
  return MODE_SCOPED_SETTINGS.rows
    .filter((/** @type {RegistryRow} */ row) => row.modes === 'both' || row.modes === mode)
    .map((/** @type {RegistryRow} */ row) => {
      const held = rowValue(slice, row);
      // Save leaves a row out at its default (`recordsOff` when empty): no `value`.
      return held ? { domain: row.domain, path: row.path, value: normalized(row, held.value) } : { domain: row.domain, path: row.path };
    });
};

const document = {
  schema: 1,
  description: 'Each diagram mode\'s default value of every MODE_SCOPED_SETTINGS row it holds, as the Web builds a new drawing (tools/generate_mode_slice_defaults.mjs); a row without value is one Save leaves out at its default; normalize lists the row fields a reader strips of run-time identities before comparing.',
  registryRevision: MODE_SCOPED_SETTINGS.revision,
  normalize: NORMALIZE,
  earlierDefaults: EARLIER_DEFAULTS,
  earlierPaletteRows: EARLIER_PALETTE_ROWS,
  defaults: { circular: modeDefaults('circular'), linear: modeDefaults('linear') }
};
// One row per line, so a default change is a one-line diff.
const render = () => {
  const head = Object.entries(document)
    .filter(([key]) => key !== 'defaults')
    .map(([key, value]) => `  ${JSON.stringify(key)}: ${JSON.stringify(value)}`);
  const modes = Object.entries(document.defaults).map(([mode, rows]) => (
    `    ${JSON.stringify(mode)}: [\n${rows.map((row) => `      ${JSON.stringify(row)}`).join(',\n')}\n    ]`
  ));
  return `{\n${head.join(',\n')},\n  "defaults": {\n${modes.join(',\n')}\n  }\n}\n`;
};

// The same data as a frozen ES module; `normalize`, `earlierDefaults`, `earlierPaletteRows` and `defaults` only (the
// description and schema name the JSON file).
const renderModule = () => `// @ts-check
// Generated by tools/generate_mode_slice_defaults.mjs from the Web's new
// drawings (createDrawingState with the palette file applied, read as Save
// writes a mode slice); gbdraw/data/mode-slice-defaults.json holds the same
// data for Python. Do not edit manually.
/** @import { UnarrangedFieldRule } from './services/linear-record-layout.js' */

/**
 * One MODE_SCOPED_SETTINGS row's default in one mode; no \`value\` when Save
 * leaves the row out at its default.
 * @typedef {{ readonly domain: string, readonly path: string, readonly value?: unknown }} ModeSliceDefault
 */
/**
 * The rows a released palette held that the palette has since changed: by
 * palette name, each feature type's color in release \`release\`.
 * @typedef {{
 *   readonly release: string,
 *   readonly reason: string,
 *   readonly palettes: Readonly<Record<string, Readonly<Record<string, string>>>>
 * }} EarlierPaletteRows
 */
/**
 * @typedef {{
 *   readonly registryRevision: number,
 *   readonly normalize: readonly UnarrangedFieldRule[],
 *   readonly earlierDefaults: readonly (ModeSliceDefault & { readonly reason: string })[],
 *   readonly earlierPaletteRows: EarlierPaletteRows,
 *   readonly defaults: Readonly<Record<'circular' | 'linear', readonly ModeSliceDefault[]>>
 * }} ModeSliceDefaults
 */

/**
 * @template T
 * @param {T} value
 * @returns {T}
 */
const deepFreeze = (value) => {
  if (!value || typeof value !== 'object' || Object.isFrozen(value)) return value;
  Object.values(value).forEach((entry) => deepFreeze(entry));
  return Object.freeze(value);
};

/** @type {ModeSliceDefaults} */
export const MODE_SLICE_DEFAULTS = deepFreeze({
  "registryRevision": ${JSON.stringify(document.registryRevision)},
  "normalize": ${JSON.stringify(document.normalize)},
  "earlierDefaults": ${JSON.stringify(document.earlierDefaults)},
  "earlierPaletteRows": ${JSON.stringify(document.earlierPaletteRows)},
  "defaults": {
${Object.entries(document.defaults).map(([mode, rows]) => (
    `    ${JSON.stringify(mode)}: [\n${rows.map((row) => `      ${JSON.stringify(row)}`).join(',\n')}\n    ]`
  )).join(',\n')}
  }
});
`;

const outputs = [[TARGET, render()], [MODULE_TARGET, renderModule()]];
if (process.argv.includes('--check')) {
  const stale = outputs.filter(([file, text]) => {
    const target = path.join(ROOT, file);
    return !fs.existsSync(target) || fs.readFileSync(target, 'utf8') !== text;
  });
  stale.forEach(([file]) => console.log(`${file} is stale. Run: node tools/generate_mode_slice_defaults.mjs`));
  if (stale.length) process.exit(1);
} else {
  outputs.forEach(([file, text]) => fs.writeFileSync(path.join(ROOT, file), text));
}
