// @ts-check
// Repairs a parsed stored value against a contract description
// (utils/contract-shape.js). For each place where the value leaves its shape,
// the first rule that applies:
//   1. the policy has a lossless reading that conforms: it is used, no row;
//   2. an unknown key is dropped;
//   3. a field with a default in the description takes it;
//   4. an optional field is dropped;
//   5. an element of an array that nothing names by position is dropped;
//   6. otherwise the enclosing part (policy.partBoundary) is dropped, outward
//      until one can be; past the root nothing is left.
// Rules 3 and 4 apply where the violation is; rules 5 and 6 to the nearest
// enclosing element or part. The rows say what was repaired, not how to word
// it: each reader maps them to its own codes. The input is never changed; a
// value that needs no repair comes back as the same object.

import { isPlainObject, shapeViolations } from './contract-shape.js';

/** @import { ContractDescription, ContractShape, JsonValue, ShapeViolation, ShapeViolationReason } from './contract-shape.js' */

/**
 * @typedef {2 | 3 | 4 | 5 | 6} ContractRepairRule
 * @typedef {'dropped' | 'defaulted'} ContractRepairEffect
 */
/**
 * One repair. `path` is the place in the file that was dropped or defaulted;
 * `pattern` is the same place relative to the root, every index written `[*]`
 * (what the policy hooks receive); `reason` is that of the first violation
 * found there, or what the reference hook returned; `used` is the default.
 * @template {string} R
 * @typedef {object} ContractRepair
 * @property {string} path
 * @property {string} pattern
 * @property {ShapeViolationReason | R} reason
 * @property {ContractRepairRule} rule
 * @property {ContractRepairEffect} effect
 * @property {JsonValue} [used]
 */
/**
 * What the reader knows that the description does not. Every hook receives a
 * pattern: the place relative to the root, every index written `[*]`, an entry
 * of a record written `.*` (`records[*].display`, `comparisons[*].pairs`).
 * @template {string} [R=never]
 * @typedef {object} RepairPolicy
 * @property {readonly string[]} positionalArrays  arrays whose elements other places name by position: rule 5 never drops their elements
 * @property {(pattern: string) => boolean} partBoundary  the places rule 6 may drop as one part
 * @property {readonly (readonly string[])[]} [parallelArrays]  groups of sibling arrays read column by column: rule 5 drops a column from every array of the group
 * @property {(pattern: string, value: unknown) => ({ readonly value: unknown } | null)} [lossless]  the lossless reading of a value that leaves its shape, if there is one
 * @property {(pattern: string, value: unknown) => (R | null)} [reference]  why a conforming scalar names nothing the file holds, or null
 * @property {(pattern: string) => (string | null)} [discriminator]  the field whose literal picks the option of the union at this place
 */
/**
 * The root's type in `description.types` and its path in the file.
 * @typedef {{ readonly type: string, readonly path: string }} ContractRoot
 */
/**
 * @template {string} R
 * @typedef {{ readonly value: unknown, readonly rows: readonly ContractRepair<R>[] }} ContractRepairResult
 */

/**
 * @typedef {{ readonly path: string, readonly pattern: string }} Place
 * @typedef {Exclude<ContractShape, { kind: 'ref' }>} ResolvedShape
 */
/**
 * @template {string} R
 * @typedef {{ readonly path: string, readonly pattern: string, readonly reason: ShapeViolationReason | R }} Origin
 */
/**
 * A place is kept (repaired, `keptIndices` for an array: the original index of
 * each element left) or must go, for the violation at `origin`.
 * @template {string} R
 * @typedef {{ readonly dropped: false, readonly value: unknown, readonly rows: readonly ContractRepair<R>[], readonly keptIndices?: readonly number[] }
 *   | { readonly dropped: true, readonly origin: Origin<R> }} Outcome
 */
/**
 * @template {string} R
 * @typedef {object} RepairContext
 * @property {Readonly<Record<string, ContractShape>>} types
 * @property {RepairPolicy<R>} policy
 * @property {ReadonlySet<string>} positional
 * @property {ReadonlyMap<string, readonly string[]>} groups  each grouped array pattern to its group
 */

/** @param {Place} at @param {string} key @returns {Place} */
const fieldPlace = (at, key) => ({
  path: at.path ? `${at.path}.${key}` : key,
  pattern: at.pattern ? `${at.pattern}.${key}` : key
});
/** @param {Place} at @param {number} index @returns {Place} */
const indexPlace = (at, index) => ({ path: `${at.path}[${index}]`, pattern: `${at.pattern}[*]` });
/** @param {Place} at @param {string} key @returns {Place} */
const entryPlace = (at, key) => ({ ...fieldPlace(at, key), pattern: at.pattern ? `${at.pattern}.*` : '*' });

/**
 * @param {ContractShape} shape
 * @param {Readonly<Record<string, ContractShape>>} types
 * @returns {ResolvedShape}
 */
const resolve = (shape, types) => (shape.kind === 'ref' ? resolve(types[shape.name], types) : shape);

/**
 * @template {string} R
 * @param {unknown} value
 * @param {readonly ContractRepair<NoInfer<R>>[]} rows
 * @param {readonly number[]} [keptIndices]
 * @returns {Outcome<R>}
 */
const kept = (value, rows, keptIndices) => ({ dropped: false, value, rows, ...(keptIndices ? { keptIndices } : {}) });
/**
 * @template {string} R
 * @param {Place} at
 * @param {ShapeViolationReason | NoInfer<R>} reason
 * @returns {Outcome<R>}
 */
const dropped = (at, reason) => ({ dropped: true, origin: { ...at, reason } });
/**
 * @template {string} R
 * @param {Place} at
 * @param {ShapeViolationReason | NoInfer<R>} reason
 * @param {ContractRepairRule} rule
 * @param {ContractRepairEffect} effect
 * @param {JsonValue} [used]
 * @returns {ContractRepair<R>}
 */
const repairRow = (at, reason, rule, effect, used) => ({
  path: at.path, pattern: at.pattern, reason, rule, effect, ...(used === undefined ? {} : { used })
});
/** @param {readonly ShapeViolation[]} found @returns {ShapeViolationReason} */
const firstReason = (found) => found[0].reason;

/**
 * The option of a union that `value` is repaired against: the first that it
 * fits; else the one its discriminator names (none when it names none); else
 * the first with the fewest violations, as the checker reports.
 * @template {string} R
 * @param {unknown} value
 * @param {Extract<ResolvedShape, { kind: 'union' }>} shape
 * @param {Place} at
 * @param {RepairContext<R>} context
 * @returns {ContractShape | null}
 */
const unionOption = (value, shape, at, context) => {
  const counts = shape.options.map((option) => shapeViolations(value, option, context.types, at.path).length);
  const fitting = counts.indexOf(0);
  if (fitting >= 0) return shape.options[fitting];
  const field = context.policy.discriminator?.(at.pattern) ?? null;
  if (field !== null) {
    if (!isPlainObject(value)) return null;
    return shape.options.find((option) => {
      const resolved = resolve(option, context.types);
      if (resolved.kind !== 'object' || !Object.hasOwn(resolved.fields, field)) return false;
      const tag = resolve(resolved.fields[field], context.types);
      return tag.kind === 'literal' && tag.values.some((item) => item === value[field]);
    }) ?? null;
  }
  return shape.options[counts.indexOf(Math.min(...counts))];
};

/**
 * Drops the same columns from every array of each parallel group among the
 * fields of one object, after each array was repaired on its own.
 * @template {string} R
 * @param {Map<string, Outcome<R>>} outcomes
 * @param {Readonly<Record<string, unknown>>} value
 * @param {Place} at
 * @param {RepairContext<R>} context
 */
const alignParallelColumns = (outcomes, value, at, context) => {
  const seen = new Set();
  for (const key of outcomes.keys()) {
    const group = context.groups.get(fieldPlace(at, key).pattern);
    if (!group || seen.has(group)) continue;
    seen.add(group);
    const members = [...outcomes.keys()].filter((name) => group.includes(fieldPlace(at, name).pattern));
    /** @type {Set<number>} */
    const removed = new Set();
    for (const name of members) {
      const outcome = outcomes.get(name);
      const original = value[name];
      if (!outcome || outcome.dropped || !outcome.keptIndices || !Array.isArray(original)) continue;
      const left = new Set(outcome.keptIndices);
      original.forEach((_, index) => { if (!left.has(index)) removed.add(index); });
    }
    if (!removed.size) continue;
    for (const name of members) {
      const outcome = outcomes.get(name);
      if (!outcome || outcome.dropped || !Array.isArray(outcome.value)) continue;
      const items = /** @type {readonly unknown[]} */ (outcome.value);
      const indices = outcome.keptIndices ?? items.map((_, index) => index);
      const keep = indices.map((index) => !removed.has(index));
      if (keep.every(Boolean)) continue;
      outcomes.set(name, kept(items.filter((_, j) => keep[j]), outcome.rows, indices.filter((_, j) => keep[j])));
    }
  }
};

/**
 * @template {string} R
 * @param {unknown} value
 * @param {ContractShape} shape
 * @param {Place} at
 * @param {boolean} conforming  the value is known to fit `shape`
 * @param {RepairContext<R>} context
 * @returns {Outcome<R>}
 */
const repairPlace = (value, shape, at, conforming, context) => {
  const { types, policy } = context;
  const resolved = resolve(shape, types);
  const found = conforming ? [] : shapeViolations(value, resolved, types, at.path);
  if (found.length) {
    const reading = policy.lossless?.(at.pattern, value) ?? null;
    if (reading && !shapeViolations(reading.value, resolved, types, at.path).length) {
      return repairPlace(reading.value, resolved, at, true, context);
    }
  } else if (!policy.reference) {
    return kept(value, []);
  }
  const fits = !found.length;
  switch (resolved.kind) {
    case 'unknown':
      return kept(value, []);
    case 'string': case 'integer': case 'number': case 'boolean': case 'null': case 'literal': {
      if (!fits) return dropped(at, firstReason(found));
      const reason = policy.reference?.(at.pattern, value) ?? null;
      return reason === null ? kept(value, []) : dropped(at, reason);
    }
    case 'union': {
      const option = unionOption(value, resolved, at, context);
      return option ? repairPlace(value, option, at, fits, context) : dropped(at, firstReason(found));
    }
    case 'array': case 'tuple': {
      if (!Array.isArray(value)) return dropped(at, firstReason(found));
      const items = /** @type {readonly unknown[]} */ (value);
      if (resolved.kind === 'tuple' && items.length !== resolved.items.length) return dropped(at, firstReason(found));
      const removable = resolved.kind === 'array' && !context.positional.has(at.pattern);
      /** @type {ContractRepair<R>[]} */
      const rows = [];
      /** @type {unknown[]} */
      const next = [];
      /** @type {number[]} */
      const keptIndices = [];
      let changed = false;
      for (let index = 0; index < items.length; index += 1) {
        const place = indexPlace(at, index);
        const outcome = repairPlace(
          items[index], resolved.kind === 'tuple' ? resolved.items[index] : resolved.items, place, fits, context
        );
        if (outcome.dropped) {
          if (!removable) return outcome;
          rows.push(repairRow(place, outcome.origin.reason, 5, 'dropped'));
          changed = true;
          continue;
        }
        rows.push(...outcome.rows);
        next.push(outcome.value);
        keptIndices.push(index);
        if (outcome.value !== items[index]) changed = true;
      }
      return changed ? kept(next, rows, keptIndices) : kept(value, [], keptIndices);
    }
    case 'record': {
      if (!isPlainObject(value)) return dropped(at, firstReason(found));
      /** @type {ContractRepair<R>[]} */
      const rows = [];
      /** @type {Record<string, unknown>} */
      const next = {};
      let changed = false;
      for (const key of Object.keys(value)) {
        const place = entryPlace(at, key);
        const outcome = repairPlace(value[key], resolved.values, place, fits, context);
        if (outcome.dropped) {
          const direct = outcome.origin.path === place.path;
          if (!direct && !policy.partBoundary(place.pattern)) return outcome;
          rows.push(repairRow(place, outcome.origin.reason, direct ? 4 : 6, 'dropped'));
          changed = true;
          continue;
        }
        rows.push(...outcome.rows);
        next[key] = outcome.value;
        if (outcome.value !== value[key]) changed = true;
      }
      return changed ? kept(next, rows) : kept(value, []);
    }
    case 'object': {
      if (!isPlainObject(value)) return dropped(at, firstReason(found));
      const { fields, optional, defaults = {} } = resolved;
      /** @type {ContractRepair<R>[]} */
      const rows = [];
      /** @type {Record<string, unknown>} */
      const next = {};
      let changed = false;
      for (const key of Object.keys(fields)) {
        if (Object.hasOwn(value, key) || optional.includes(key)) continue;
        if (!Object.hasOwn(defaults, key)) return dropped(fieldPlace(at, key), 'missing-key');
      }
      /** @type {Map<string, Outcome<R>>} */
      const outcomes = new Map();
      for (const key of Object.keys(value)) {
        if (Object.hasOwn(fields, key)) outcomes.set(key, repairPlace(value[key], fields[key], fieldPlace(at, key), fits, context));
      }
      if (context.groups.size) alignParallelColumns(outcomes, value, at, context);
      for (const key of Object.keys(value)) {
        const place = fieldPlace(at, key);
        const outcome = outcomes.get(key);
        if (!outcome) {
          if (resolved.open) {
            next[key] = value[key];
          } else {
            rows.push(repairRow(place, 'unknown-key', 2, 'dropped'));
            changed = true;
          }
          continue;
        }
        if (!outcome.dropped) {
          rows.push(...outcome.rows);
          next[key] = outcome.value;
          if (outcome.value !== value[key]) changed = true;
          continue;
        }
        const direct = outcome.origin.path === place.path;
        if (!direct && !policy.partBoundary(place.pattern)) return outcome;
        if (Object.hasOwn(defaults, key)) {
          const used = defaults[key];
          next[key] = structuredClone(used);
          rows.push(repairRow(place, outcome.origin.reason, direct ? 3 : 6, 'defaulted', used));
        } else if (optional.includes(key)) {
          rows.push(repairRow(place, outcome.origin.reason, direct ? 4 : 6, 'dropped'));
        } else {
          return outcome;
        }
        changed = true;
      }
      for (const key of Object.keys(fields)) {
        if (Object.hasOwn(value, key) || optional.includes(key)) continue;
        const used = defaults[key];
        next[key] = structuredClone(used);
        rows.push(repairRow(fieldPlace(at, key), 'missing-key', 3, 'defaulted', used));
        changed = true;
      }
      return changed ? kept(next, rows) : kept(value, []);
    }
  }
};

/**
 * Repairs `value`, read from the file at `root.path`, against the type
 * `root.type` of `description`. Returns the repaired value and one row per
 * repair; `value` is undefined, with one row for the root, when rule 6 drops
 * the root. A value that needs no repair comes back as the same object, and
 * `value` itself is never changed.
 * @template {string} [R=never]
 * @param {unknown} value
 * @param {ContractRoot} root
 * @param {ContractDescription} description
 * @param {RepairPolicy<R>} policy
 * @returns {ContractRepairResult<R>}
 */
export const repairByDescription = (value, root, description, policy) => {
  /** @type {Map<string, readonly string[]>} */
  const groups = new Map();
  (policy.parallelArrays ?? []).forEach((group) => group.forEach((pattern) => groups.set(pattern, group)));
  const outcome = repairPlace(value, { kind: 'ref', name: root.type }, { path: root.path, pattern: '' }, false, {
    types: description.types, policy, positional: new Set(policy.positionalArrays), groups
  });
  if (outcome.dropped) {
    return { value: undefined, rows: [repairRow({ path: root.path, pattern: '' }, outcome.origin.reason, 6, 'dropped')] };
  }
  return { value: outcome.value, rows: outcome.rows };
};
