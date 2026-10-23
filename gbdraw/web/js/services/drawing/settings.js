// @ts-check
// The settings slice of a drawing: its sixteen draft settings (Session 46
// `modes.<m>.config` and the drawing's `ui` values), among them what the
// Default colors and Blacklist files set (`fileImports`, Session 48, OV-1280;
// its writers are app/watchers.js, the Blacklist import, and app/results.js,
// the Default colors import; services/file-imports.js holds the pure helpers).
// Unlike an editor slice,
// the controls and owners that edit a setting keep writing it (plan-r2 §13
// items 3 and 4): this module owns the whole-member paths only. Load and a
// History restore write members with `install`, the live-input rollback and
// a History fallback read them with `capture`, and Reset Settings returns
// them with `reset`. A caller passes values it normalized; the members stay
// untyped (state.js builds them with Vue's untyped `reactive`).
import { reactive, ref } from '../../utils/reactivity.js';
import { cloneJsonData } from '../json-clone.js';
import { createSliceHandle, sliceMembers } from './slice-handle.js';
/** @import { Ref } from '../../utils/reactivity.js' */
/** @import { SliceHandle } from './slice-handle.js' */

const SETTINGS_MEMBERS = /** @type {const} */ ([
  'form', 'adv', 'losat', 'losatProgram', 'circularConservation', 'linearComparisonPlan',
  'linearRecordLayoutEnabled', 'linearRecordGap', 'linearRecordRows', 'recordDisplayDrafts', 'recordsOff',
  'unmanagedConfigOverrides', 'importedComparisonIntent', 'layoutPreferences', 'linearTypographyLinked', 'fileImports'
]);
/** @typedef {typeof SETTINGS_MEMBERS[number]} SettingsMember */
/** @typedef {{ [K in SettingsMember]: unknown }} SettingsMembers The members, untyped (plan-r2 §13 item 4). */
/** @typedef {Partial<SettingsMembers>} SettingsData Plain values: what `capture` returns and `install` writes. */
/** @typedef {SliceHandle<'settings'>} SettingsSlice The handle `install`, `capture`, and `reset` take. */
/** @typedef {{ mode: unknown, defaultSource: unknown, edges: object[] }} ComparisonPlanMember */
/** @typedef {{ [key: string]: unknown }} ObjectMember */

// How a member holds its value: a ref, a list, or an object; the comparison
// plan and the layout preferences are objects with members of their own.
/** @type {ReadonlySet<SettingsMember>} */
const REFS = new Set(['losatProgram', 'linearRecordLayoutEnabled', 'linearRecordGap', 'linearTypographyLinked']);
/** @type {ReadonlySet<SettingsMember>} */
const LISTS = new Set(['linearRecordRows', 'recordDisplayDrafts', 'recordsOff']);
// Reset Settings keeps the record gap and rows, which records are OFF, and the
// imported comparison intent; it returns the others to their defaults.
/** @type {readonly SettingsMember[]} */
const RESET_MEMBERS = SETTINGS_MEMBERS.filter((member) => ![
  'linearRecordGap', 'linearRecordRows', 'recordsOff', 'importedComparisonIntent'
].includes(member));

// A settings handle exists only through `create`, which gives it these members.
/** @param {SettingsSlice} slice */
const membersOf = (slice) => /** @type {SettingsMembers} */ (sliceMembers(slice, 'settings'));
/** @type {WeakMap<SettingsSlice, () => SettingsData>} */
const defaultsOf = new WeakMap();

/**
 * A member that holds `value`, as `create` makes one a caller did not pass.
 * @param {SettingsMember} member
 * @param {unknown} value
 */
const memberOf = (member, value) => {
  if (REFS.has(member)) return ref(cloneJsonData(value) ?? null);
  return reactive(cloneJsonData(value) ?? (LISTS.has(member) ? [] : {}));
};

/**
 * @param {ObjectMember} target
 * @param {ObjectMember} source
 */
const assignNested = (target, source) => Object.entries(source).forEach(([key, value]) => {
  const current = target[key];
  if (value && typeof value === 'object' && !Array.isArray(value) && current && typeof current === 'object') {
    assignNested(/** @type {ObjectMember} */ (current), /** @type {ObjectMember} */ (value));
  } else {
    target[key] = value;
  }
});

// Writes `value` into the member in place (its object stays): a ref's value,
// a list's rows, an object's fields. The plan's edges are copied one level
// (an edge holds its comparison file); the layout preferences keep their
// nested objects, which `activeLayoutPreferences` reads.
/**
 * @param {SettingsMember} member
 * @param {unknown} target
 * @param {unknown} value
 */
const writeMember = (member, target, value) => {
  if (REFS.has(member)) {
    /** @type {Ref<unknown>} */ (target).value = cloneJsonData(value);
  } else if (LISTS.has(member)) {
    const rows = /** @type {unknown[]} */ (target);
    rows.splice(0, rows.length, .../** @type {unknown[]} */ (cloneJsonData(value)));
  } else if (member === 'linearComparisonPlan') {
    const plan = /** @type {ComparisonPlanMember} */ (target);
    const { mode, defaultSource, edges } = /** @type {ComparisonPlanMember} */ (value);
    plan.mode = mode;
    plan.defaultSource = defaultSource;
    if (!Array.isArray(plan.edges)) plan.edges = [];
    plan.edges.splice(0, plan.edges.length, ...edges.map((edge) => ({ ...edge })));
  } else if (member === 'layoutPreferences') {
    assignNested(/** @type {ObjectMember} */ (target), /** @type {ObjectMember} */ (cloneJsonData(value)));
  } else {
    const fields = /** @type {ObjectMember} */ (target);
    Object.keys(fields).forEach((key) => { delete fields[key]; });
    Object.assign(fields, cloneJsonData(value));
  }
};

/**
 * @param {SettingsMember} member
 * @param {unknown} source
 */
const copyOf = (member, source) => {
  if (REFS.has(member)) return cloneJsonData(/** @type {Ref<unknown>} */ (source).value);
  if (member !== 'linearComparisonPlan') return cloneJsonData(source);
  const { mode, defaultSource, edges } = /** @type {ComparisonPlanMember} */ (source);
  return { mode, defaultSource, edges: edges.map((edge) => ({ ...edge })) };
};

/**
 * A slice over a drawing's members, as state.js builds them, or a new one (a
 * test fixture); a member not passed holds its default.
 * @template {Partial<SettingsMembers>} M
 * @param {M} [members]
 * @param {() => SettingsData} [defaults] Fresh default values; Reset Settings installs them.
 * @returns {{ handle: SettingsSlice, view: Readonly<SettingsMembers & M> }}
 */
const create = (members = /** @type {M} */ ({}), defaults = () => ({})) => {
  /** @type {SettingsData | undefined} */
  let initial;
  const data = /** @type {SettingsMembers} */ (Object.fromEntries(SETTINGS_MEMBERS.map((member) => [
    member, members[member] ?? memberOf(member, (initial ??= defaults())[member])
  ])));
  const handle = createSliceHandle('settings', data);
  defaultsOf.set(handle, defaults);
  return { handle, view: /** @type {Readonly<SettingsMembers & M>} */ (Object.freeze({ ...data })) };
};

/**
 * Copies of the members `keys` names (all of them by default), as plain data;
 * a plan edge keeps its comparison file.
 * @param {SettingsSlice} slice
 * @param {readonly SettingsMember[]} [keys]
 * @returns {SettingsData}
 */
const capture = (slice, keys = SETTINGS_MEMBERS) => {
  const members = membersOf(slice);
  return Object.fromEntries(SETTINGS_MEMBERS.filter((member) => keys.includes(member))
    .map((member) => [member, copyOf(member, members[member])]));
};

// Writes the members `data` holds, each in place; a member it leaves out (or
// holds as `undefined`) stays.
/**
 * @param {SettingsSlice} slice
 * @param {SettingsData} data
 */
const install = (slice, data) => {
  const members = membersOf(slice);
  SETTINGS_MEMBERS.forEach((member) => {
    if (data[member] !== undefined) writeMember(member, members[member], data[member]);
  });
};

// Reset Settings: the reset members take their defaults. The comparison plan
// keeps the edges that hold a file or a LOSAT file name, not included and
// with their file and name inactive.
/** @param {SettingsSlice} slice */
const reset = (slice) => {
  const defaults = defaultsOf.get(slice)?.() ?? {};
  const plan = /** @type {ComparisonPlanMember} */ (membersOf(slice).linearComparisonPlan);
  const keptEdges = (Array.isArray(plan?.edges) ? plan.edges : [])
    .filter((edge) => Boolean(Reflect.get(edge, 'file')) || String(Reflect.get(edge, 'losatFilename') || '').trim())
    .map((edge) => ({ ...edge, included: false, fileActive: false, losatFilenameActive: false }));
  install(slice, Object.fromEntries(RESET_MEMBERS.map((member) => [member, member === 'linearComparisonPlan' && defaults[member]
    ? { .../** @type {object} */ (defaults[member]), edges: keptEdges }
    : defaults[member]])));
};

export const SETTINGS_SLICE = Object.freeze({ members: SETTINGS_MEMBERS, create, capture, install, reset });
