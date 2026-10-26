#!/usr/bin/env node
// Writes each diagram mode's default settings, as the Web builds them, to
// gbdraw/data/mode-slice-defaults.json, which Python reads to decide whether
// a Session 46-48 used the mode it did not show (rule d of the move to
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

// The same data as a frozen ES module; `normalize` and `defaults` only (the
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
 * @typedef {{
 *   readonly registryRevision: number,
 *   readonly normalize: readonly UnarrangedFieldRule[],
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
