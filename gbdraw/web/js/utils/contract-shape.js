// @ts-check
// Structural check of a parsed JSON value against a description of
// gbdraw/contract_shapes.py, the twin of its `shape_violations` (tools/twins.json;
// both run tests/fixtures/contract_shape_vectors.json). A description is data
// that tools/generate_web_contracts.py renders into contracts/*-shape.generated.js.

/**
 * @typedef {null | boolean | number | string | readonly JsonValue[] | { readonly [key: string]: JsonValue }} JsonValue
 */
/**
 * @typedef {{ readonly kind: 'string' | 'integer' | 'number' | 'boolean' | 'null' | 'unknown' }} ScalarShape
 * @typedef {{ readonly kind: 'literal', readonly values: readonly (string | number | boolean | null)[] }} LiteralShape
 * @typedef {{ readonly kind: 'array', readonly items: ContractShape }} ArrayShape
 * @typedef {{ readonly kind: 'tuple', readonly items: readonly ContractShape[] }} TupleShape
 * @typedef {{ readonly kind: 'record', readonly values: ContractShape }} RecordShape
 * @typedef {{ readonly kind: 'union', readonly options: readonly ContractShape[] }} UnionShape
 * @typedef {{ readonly kind: 'ref', readonly name: string }} RefShape
 */
/**
 * An object: no other key unless it is open. `defaults` gives the value the
 * Python reader uses for an optional field that is absent.
 * @typedef {object} ObjectShape
 * @property {'object'} kind
 * @property {Readonly<Record<string, ContractShape>>} fields
 * @property {readonly string[]} optional
 * @property {boolean} open
 * @property {Readonly<Record<string, JsonValue>>} [defaults]
 */
/**
 * @typedef {ScalarShape | LiteralShape | ArrayShape | TupleShape | RecordShape | UnionShape | RefShape | ObjectShape} ContractShape
 * @typedef {{ readonly root: string, readonly types: Readonly<Record<string, ContractShape>> }} ContractDescription
 */

/** Why a value leaves its shape: the codes of Python's `ViolationReason`. */
export const SHAPE_VIOLATION_REASONS = Object.freeze(/** @type {const} */ ([
  'unknown-key', 'missing-key', 'wrong-kind', 'null-not-declared', 'non-integer', 'not-in-set', 'tuple-length'
]));
/** @typedef {typeof SHAPE_VIOLATION_REASONS[number]} ShapeViolationReason */
/**
 * Where a value leaves its shape; never the value, which can be private text.
 * @typedef {{ readonly path: string, readonly reason: ShapeViolationReason }} ShapeViolation
 */

/**
 * @param {unknown} value
 * @param {'string' | 'integer' | 'number' | 'boolean' | 'null'} kind
 */
const isScalarOfKind = (value, kind) => {
  if (kind === 'null') return value === null;
  if (kind === 'integer') return Number.isInteger(value);
  return typeof value === (kind === 'number' ? 'number' : kind);
};

/**
 * @param {unknown} value
 * @param {string} kind
 * @param {string} path
 * @returns {ShapeViolation}
 */
const kindViolation = (value, kind, path) => {
  if (value === null) return { path, reason: 'null-not-declared' };
  if (kind === 'integer' && typeof value === 'number') return { path, reason: 'non-integer' };
  return { path, reason: 'wrong-kind' };
};

/**
 * A JSON object: not null and not an array.
 * @param {unknown} value
 * @returns {value is Readonly<Record<string, unknown>>}
 */
export const isPlainObject = (value) => typeof value === 'object' && value !== null && !Array.isArray(value);

/**
 * Where a parsed JSON `value` leaves `shape`, in order. `types` holds the named
 * shapes a `ref` points to. Of a union's options, the first one with the fewest
 * violations is reported.
 * @param {unknown} value
 * @param {ContractShape} shape
 * @param {Readonly<Record<string, ContractShape>>} types
 * @param {string} path
 * @returns {ShapeViolation[]}
 */
export const shapeViolations = (value, shape, types, path) => {
  switch (shape.kind) {
    case 'unknown':
      return [];
    case 'ref':
      // A generated description names every shape it refers to.
      return shapeViolations(value, types[shape.name], types, path);
    case 'string': case 'integer': case 'number': case 'boolean': case 'null':
      return isScalarOfKind(value, shape.kind) ? [] : [kindViolation(value, shape.kind, path)];
    case 'literal':
      if (shape.values.some((item) => item === value)) return [];
      return [{ path, reason: value === null ? 'null-not-declared' : 'not-in-set' }];
    case 'union': {
      /** @type {ShapeViolation[] | null} */
      let fewest = null;
      for (const option of shape.options) {
        const found = shapeViolations(value, option, types, path);
        if (!found.length) return [];
        if (!fewest || found.length < fewest.length) fewest = found;
      }
      return fewest || [];
    }
    case 'array': case 'tuple': {
      if (!Array.isArray(value)) return [kindViolation(value, shape.kind, path)];
      const items = /** @type {readonly unknown[]} */ (value);
      if (shape.kind === 'tuple' && items.length !== shape.items.length) return [{ path, reason: 'tuple-length' }];
      return items.flatMap((item, index) => shapeViolations(
        item, shape.kind === 'tuple' ? shape.items[index] : shape.items, types, `${path}[${index}]`
      ));
    }
    case 'record':
      if (!isPlainObject(value)) return [kindViolation(value, shape.kind, path)];
      return Object.keys(value).flatMap((key) => shapeViolations(value[key], shape.values, types, `${path}.${key}`));
    case 'object': {
      if (!isPlainObject(value)) return [kindViolation(value, shape.kind, path)];
      const { fields } = shape;
      const keys = Object.keys(value);
      return [
        ...(shape.open ? [] : keys.filter((key) => !Object.hasOwn(fields, key))
          .map((key) => /** @type {ShapeViolation} */ ({ path: `${path}.${key}`, reason: 'unknown-key' }))),
        ...Object.keys(fields).filter((key) => !Object.hasOwn(value, key) && !shape.optional.includes(key))
          .map((key) => /** @type {ShapeViolation} */ ({ path: `${path}.${key}`, reason: 'missing-key' })),
        ...keys.filter((key) => Object.hasOwn(fields, key))
          .flatMap((key) => shapeViolations(value[key], fields[key], types, `${path}.${key}`))
      ];
    }
  }
};
