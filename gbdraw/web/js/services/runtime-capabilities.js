// @ts-check
import { EXPECTED_WEB_RUNTIME_CAPABILITIES } from '../contracts/schema.generated.js';

export { EXPECTED_WEB_RUNTIME_CAPABILITIES };
export const DIAGRAM_ENGINE_COMPATIBILITY_MESSAGE =
  'The diagram engine is incompatible with this Web app. Reload the page; if the problem persists, contact the site administrator.';

export class DiagramRuntimeCompatibilityError extends Error {
  /** @param {unknown} diagnostic Anything printable; it is passed through `String` when truthy. */
  constructor(diagnostic) {
    super(DIAGRAM_ENGINE_COMPATIBILITY_MESSAGE);
    this.name = 'DiagramRuntimeCompatibilityError';
    this.diagnostic = String(diagnostic || 'Runtime capability manifest mismatch.');
  }
}

/**
 * @param {unknown} value
 * @returns {value is Record<string, unknown>}
 */
const isPlainObject = (value) => {
  if (!value || typeof value !== 'object' || Array.isArray(value)) return false;
  const prototype = Object.getPrototypeOf(value);
  return prototype === Object.prototype || prototype === null;
};

/** @param {unknown} value */
const valueKind = (value) => {
  if (Array.isArray(value)) return 'array';
  if (value === null) return 'null';
  return typeof value;
};

/**
 * @param {unknown} expected
 * @param {unknown} actual
 * @param {string} [path]
 * @returns {string} The first mismatch message, or '' when the manifests agree.
 */
const findManifestMismatch = (expected, actual, path = 'capabilities') => {
  if (Array.isArray(expected)) {
    if (!Array.isArray(actual)) {
      return `${path} must be an array, received ${valueKind(actual)}.`;
    }
    if (actual.length !== expected.length) {
      return `${path} must contain ${expected.length} entries, received ${actual.length}.`;
    }
    for (let index = 0; index < expected.length; index += 1) {
      const mismatch = findManifestMismatch(expected[index], actual[index], `${path}[${index}]`);
      if (mismatch) return mismatch;
    }
    return '';
  }

  if (isPlainObject(expected)) {
    if (!isPlainObject(actual)) {
      return `${path} must be an object, received ${valueKind(actual)}.`;
    }
    const expectedKeys = Object.keys(expected).sort();
    const actualKeys = Object.keys(actual).sort();
    if (
      expectedKeys.length !== actualKeys.length ||
      expectedKeys.some((key, index) => key !== actualKeys[index])
    ) {
      return `${path} fields must be [${expectedKeys.join(', ')}], received [${actualKeys.join(', ')}].`;
    }
    for (const key of expectedKeys) {
      const mismatch = findManifestMismatch(expected[key], actual[key], `${path}.${key}`);
      if (mismatch) return mismatch;
    }
    return '';
  }

  if (actual !== expected) {
    return `${path} must be ${JSON.stringify(expected)}, received ${JSON.stringify(actual)}.`;
  }
  return '';
};

/**
 * @template T
 * @param {T} value
 * @returns {T}
 */
const deepFreeze = (value) => {
  if (!value || typeof value !== 'object' || Object.isFrozen(value)) return value;
  Object.values(value).forEach((item) => deepFreeze(item));
  return Object.freeze(value);
};

/** @param {unknown} capabilities The manifest reported by the diagram engine. */
export const validateWebRuntimeCapabilities = (capabilities) => {
  const mismatch = findManifestMismatch(EXPECTED_WEB_RUNTIME_CAPABILITIES, capabilities);
  if (mismatch) {
    throw new DiagramRuntimeCompatibilityError(mismatch);
  }
  return deepFreeze(JSON.parse(JSON.stringify(capabilities)));
};
