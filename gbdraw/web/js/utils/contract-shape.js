// @ts-check
// Structural check of a parsed JSON value against a description of
// gbdraw/contract_shapes.py, the twin of its `shape_violations` (tools/twins.json;
// both run tests/fixtures/contract_shape_vectors.json). A description is data
// that tools/generate_web_contracts.py renders into contracts/*-shape.generated.js.

/**
 * @typedef {null | boolean | number | string | readonly JsonValue[] | { readonly [key: string]: JsonValue }} JsonValue
 */
/**
 * The value domain Python checks beyond the kind (gbdraw.contract_shapes.bounded):
 * `pattern` must match somewhere in a string.
 * @typedef {{ readonly minimum?: number, readonly exclusiveMinimum?: number, readonly maximum?: number }} NumberBounds
 * @typedef {{ readonly kind: 'integer' | 'number' } & NumberBounds} NumberShape
 * @typedef {{ readonly kind: 'string', readonly pattern?: string }} StringShape
 * @typedef {{ readonly kind: 'boolean' | 'null' | 'unknown' }} OtherScalarShape
 * @typedef {NumberShape | StringShape | OtherScalarShape} ScalarShape
 */
/**
 * The other spellings of a literal's value that its reader takes
 * (gbdraw.contract_shapes.Spellings); a spelling still leaves the shape.
 * @typedef {{ readonly aliases?: Readonly<Record<string, string>>, readonly foldHyphens?: boolean, readonly blank?: string }} LiteralSpellings
 * @typedef {{ readonly kind: 'literal', readonly values: readonly (string | number | boolean | null)[], readonly spellings?: LiteralSpellings }} LiteralShape
 * @typedef {{ readonly kind: 'array', readonly items: ContractShape, readonly minItems?: number }} ArrayShape
 * @typedef {{ readonly kind: 'tuple', readonly items: readonly ContractShape[] }} TupleShape
 * @typedef {{ readonly kind: 'record', readonly values: ContractShape }} RecordShape
 * @typedef {{ readonly kind: 'union', readonly options: readonly ContractShape[] }} UnionShape
 * @typedef {{ readonly kind: 'ref', readonly name: string }} RefShape
 */
/**
 * An object: no other key unless it is open. `defaults` gives the value a
 * reader uses for a field that is absent (an optional one) or unreadable.
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
/**
 * A description with one root per version band: a document of `versions` is
 * checked against `types[root]`; `minVersion` is the band's first version.
 * @typedef {{ readonly minVersion: number, readonly versions: readonly number[], readonly root: string }} ContractBand
 * @typedef {{ readonly bands: readonly ContractBand[], readonly types: Readonly<Record<string, ContractShape>> }} BandedContractDescription
 */
/**
 * A description with one root per key: a value of `key` is checked against `types[roots[key]]`.
 * @typedef {{ readonly roots: Readonly<Record<string, string>>, readonly types: Readonly<Record<string, ContractShape>> }} RootedContractDescription
 */

/** Why a value leaves its shape: the codes of Python's `ViolationReason`. */
export const SHAPE_VIOLATION_REASONS = Object.freeze(/** @type {const} */ ([
  'unknown-key', 'missing-key', 'wrong-kind', 'null-not-declared', 'non-integer', 'not-in-set', 'tuple-length',
  'out-of-range', 'too-short', 'pattern'
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
 * @returns {value is Record<string, unknown>}
 */
export const isPlainObject = (value) => typeof value === 'object' && value !== null && !Array.isArray(value);

/**
 * Where a value of the shape's kind leaves the domain its bounds name.
 * @param {unknown} value
 * @param {NumberShape | StringShape | ArrayShape} shape
 * @param {string} path
 * @returns {ShapeViolation[]}
 */
const domainViolations = (value, shape, path) => {
  if (shape.kind === 'string') {
    return shape.pattern !== undefined && !new RegExp(shape.pattern).test(String(value))
      ? [{ path, reason: 'pattern' }] : [];
  }
  if (shape.kind === 'array') {
    return shape.minItems !== undefined && /** @type {readonly unknown[]} */ (value).length < shape.minItems
      ? [{ path, reason: 'too-short' }] : [];
  }
  const number = Number(value);
  const outside = (shape.minimum !== undefined && number < shape.minimum)
    || (shape.exclusiveMinimum !== undefined && number <= shape.exclusiveMinimum)
    || (shape.maximum !== undefined && number > shape.maximum);
  return outside ? [{ path, reason: 'out-of-range' }] : [];
};

/**
 * The name of a literal with spellings that `value` spells, or null: its text
 * without case and outer spaces, with `-` read as `_` when `foldHyphens`; an
 * alias reads as its name, and an unset value (one Python reads as false) as
 * `blank`. The twin of Python's `read_spelling`.
 * @param {unknown} value
 * @param {LiteralShape} shape
 * @returns {string | null}
 */
export const readSpelling = (value, shape) => {
  const { spellings } = shape;
  if (!spellings) return null;
  const unset = value === null || value === undefined || value === false || value === 0 || value === ''
    || (Array.isArray(value) && value.length === 0) || (isPlainObject(value) && Object.keys(value).length === 0);
  if (!unset && typeof value !== 'string') return null;
  const lowered = (unset ? (spellings.blank ?? '') : /** @type {string} */ (value)).trim().toLowerCase();
  const folded = spellings.foldHyphens ? lowered.replaceAll('-', '_') : lowered;
  const aliases = spellings.aliases ?? {};
  const name = Object.hasOwn(aliases, folded) ? aliases[folded] : folded;
  return shape.values.some((item) => item === name) ? name : null;
};

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
    case 'string': case 'integer': case 'number':
      return isScalarOfKind(value, shape.kind) ? domainViolations(value, shape, path) : [kindViolation(value, shape.kind, path)];
    case 'boolean': case 'null':
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
      return [...(shape.kind === 'array' ? domainViolations(items, shape, path) : []), ...items.flatMap((item, index) => shapeViolations(
        item, shape.kind === 'tuple' ? shape.items[index] : shape.items, types, `${path}[${index}]`
      ))];
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
