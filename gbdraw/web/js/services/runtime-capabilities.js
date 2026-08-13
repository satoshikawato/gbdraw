export const DIAGRAM_ENGINE_COMPATIBILITY_MESSAGE =
  'The diagram engine is incompatible with this Web app. Reload the page; if the problem persists, contact the site administrator.';

export const EXPECTED_WEB_RUNTIME_CAPABILITIES = Object.freeze({
  schema: 1,
  renderProtocol: 2,
  request: Object.freeze({
    currentSchema: 5,
    supportedSchemas: Object.freeze([1, 2, 5]),
    unknownFieldPolicy: 'reject'
  }),
  resources: Object.freeze({
    encodings: Object.freeze(['base64', 'gbdraw-depth-table-v1'])
  }),
  rendering: Object.freeze({
    optionSchema: 3,
    featureRenderings: Object.freeze(['arrow', 'rectangle', 'underlay']),
    circularTrackRenderers: Object.freeze([
      'annotations',
      'depth',
      'dinucleotide_content',
      'dinucleotide_skew',
      'features',
      'sequence_conservation',
      'spacer',
      'ticks'
    ]),
    linearTrackRenderers: Object.freeze([
      'annotations',
      'depth',
      'dinucleotide_content',
      'dinucleotide_skew',
      'features',
      'spacer'
    ])
  }),
  analysis: Object.freeze({
    proteinLosatCacheSchema: 4,
    nucleotideLosatCacheSchema: 2,
    losatDerivedCacheSchema: 3,
    proteinIdentityManifestSchema: 2
  })
});

export class DiagramRuntimeCompatibilityError extends Error {
  constructor(diagnostic) {
    super(DIAGRAM_ENGINE_COMPATIBILITY_MESSAGE);
    this.name = 'DiagramRuntimeCompatibilityError';
    this.diagnostic = String(diagnostic || 'Runtime capability manifest mismatch.');
  }
}

const isPlainObject = (value) => {
  if (!value || typeof value !== 'object' || Array.isArray(value)) return false;
  const prototype = Object.getPrototypeOf(value);
  return prototype === Object.prototype || prototype === null;
};

const valueKind = (value) => {
  if (Array.isArray(value)) return 'array';
  if (value === null) return 'null';
  return typeof value;
};

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

const deepFreeze = (value) => {
  if (!value || typeof value !== 'object' || Object.isFrozen(value)) return value;
  Object.values(value).forEach((item) => deepFreeze(item));
  return Object.freeze(value);
};

export const validateWebRuntimeCapabilities = (capabilities) => {
  const mismatch = findManifestMismatch(EXPECTED_WEB_RUNTIME_CAPABILITIES, capabilities);
  if (mismatch) {
    throw new DiagramRuntimeCompatibilityError(mismatch);
  }
  return deepFreeze(JSON.parse(JSON.stringify(capabilities)));
};
