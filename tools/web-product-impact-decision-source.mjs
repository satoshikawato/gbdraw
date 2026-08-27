const SCHEMA_VERSION = 1;
const START_MARKER = '<!-- gbdraw-product-impact-decision:start -->';
const END_MARKER = '<!-- gbdraw-product-impact-decision:end -->';
const DEFAULT_MAX_BODY_BYTES = 64 * 1024;
const DEFAULT_MAX_DECISION_COUNT = 8;
const MAXIMUM_BODY_BOUND = 1024 * 1024;
const MAXIMUM_DECISION_BOUND = 100;
const MAXIMUM_ARRAY_LENGTH = 32;
const MAXIMUM_ID_LENGTH = 160;
const MAXIMUM_REF_LENGTH = 1000;
const MAXIMUM_TEXT_LENGTH = 2000;
const INPUT_FIELDS = Object.freeze([
  'body',
  'eventAuthorLogin',
  'eventHeadSha',
  'maxBodyBytes',
  'maxDecisionCount'
]);
const DECLARATION_FIELDS = Object.freeze(['decisions', 'headSha', 'schemaVersion']);
const DECISION_FIELDS = Object.freeze([
  'acceptedImpactClass',
  'acknowledgedChangedRequirementRefs',
  'concernKey',
  'evidenceRefs',
  'rationale',
  'residualRisk',
  'scenarioRevision',
  'selectedOptionId'
]);
const CONCERN_KEY_PATTERN = /^[a-z][a-z0-9]*(?:[.-][a-z0-9]+)*$/;
const OPTION_ID_PATTERN = /^[a-z][a-z0-9]*(?:-[a-z0-9]+)*$/;
const REQUIREMENT_ID_PATTERN = /^REQ-[A-Z0-9]+(?:-[A-Z0-9]+)*$/;
const HEAD_SHA_PATTERN = /^[0-9a-f]{40}$/;
const LOGIN_PATTERN = /^[a-z0-9](?:[a-z0-9-]{0,37}[a-z0-9])?$/;

const compareText = (left, right) => left < right ? -1 : left > right ? 1 : 0;
const isRecord = (value) => value !== null
  && typeof value === 'object'
  && !Array.isArray(value)
  && [Object.prototype, null].includes(Object.getPrototypeOf(value));
const deepFreeze = (value) => {
  if (value && typeof value === 'object' && !Object.isFrozen(value)) {
    Object.values(value).forEach(deepFreeze);
    Object.freeze(value);
  }
  return value;
};
const addError = (errors, code, path, message) => {
  errors.push({ code, path, message });
};
const exactFields = (value, expectedFields, path, errors) => {
  const expected = new Set(expectedFields);
  Object.keys(value).sort(compareText).forEach((field) => {
    if (!expected.has(field)) {
      addError(errors, 'unknown-field', `${path}.${field}`, `Unknown field ${field}`);
    }
  });
  expectedFields.forEach((field) => {
    if (!Object.hasOwn(value, field)) {
      addError(errors, 'missing-field', `${path}.${field}`, `Missing required field ${field}`);
    }
  });
};
const normalizeText = (value) => typeof value === 'string'
  ? value.trim().replace(/\s+/g, ' ')
  : '';
const markerCount = (source, marker) => {
  let count = 0;
  let offset = 0;
  while (offset <= source.length) {
    const index = source.indexOf(marker, offset);
    if (index < 0) break;
    count += 1;
    offset = index + marker.length;
  }
  return count;
};
const byteLength = (source) => new TextEncoder().encode(source).byteLength;
const emptyDeclaration = () => ({ schemaVersion: SCHEMA_VERSION, headSha: '', decisions: [] });
const result = ({ errors, present, metadata, declaration }) => deepFreeze({
  valid: errors.length === 0,
  present,
  metadata,
  declaration,
  errors
});
const validateBound = (value, fallback, maximum, field, errors) => {
  if (value === undefined) return fallback;
  if (!Number.isInteger(value) || value < 1 || value > maximum) {
    addError(errors, 'invalid-bound', `$.${field}`, `Invalid ${field} parser bound`);
    return fallback;
  }
  return value;
};
const canonicalString = (
  value,
  path,
  errors,
  { pattern = null, maximumLength = MAXIMUM_ID_LENGTH, normalize = false } = {}
) => {
  const normalized = normalize ? normalizeText(value) : value;
  if (
    typeof normalized !== 'string'
    || !normalized
    || normalized !== normalized.trim()
    || normalized.length > maximumLength
    || /[\u0000-\u001f\u007f]/.test(normalized)
    || (pattern && !pattern.test(normalized))
  ) {
    addError(errors, 'invalid-string', path, 'Expected a bounded canonical nonempty string');
    return '';
  }
  return normalized;
};
const sortedUniqueStrings = (values, path, errors, options = {}) => {
  if (!Array.isArray(values)) {
    addError(errors, 'invalid-array', path, 'Expected an array');
    return [];
  }
  if (values.length > MAXIMUM_ARRAY_LENGTH) {
    addError(errors, 'oversized-array', path, 'Array exceeds the parser bound');
  }
  const normalized = [];
  const seen = new Set();
  let previous = null;
  values.slice(0, MAXIMUM_ARRAY_LENGTH).forEach((value, index) => {
    const item = canonicalString(value, `${path}[${index}]`, errors, options);
    if (!item) return;
    if (seen.has(item)) {
      addError(errors, 'duplicate-value', `${path}[${index}]`, 'Array values must be unique');
    }
    if (previous !== null && compareText(previous, item) > 0) {
      addError(errors, 'unsorted-array', path, 'Array values must use canonical order');
    }
    seen.add(item);
    previous = item;
    normalized.push(item);
  });
  return normalized;
};
const compareDecision = (left, right) => compareText(left.concernKey, right.concernKey)
  || left.scenarioRevision - right.scenarioRevision
  || compareText(left.selectedOptionId, right.selectedOptionId);

export const parseCurrentProductImpactDecisions = (input) => {
  const errors = [];
  const fallbackMetadata = { eventAuthorLogin: '', eventHeadSha: '' };
  if (!isRecord(input)) {
    addError(errors, 'invalid-input', '$', 'Expected one parser input object');
    return result({
      errors,
      present: false,
      metadata: fallbackMetadata,
      declaration: emptyDeclaration()
    });
  }
  Object.keys(input).sort(compareText).forEach((field) => {
    if (!INPUT_FIELDS.includes(field)) {
      addError(errors, 'unknown-field', `$.${field}`, `Unknown parser input field ${field}`);
    }
  });
  for (const field of ['body', 'eventHeadSha', 'eventAuthorLogin']) {
    if (typeof input[field] !== 'string') {
      addError(errors, 'invalid-input-string', `$.${field}`, `Expected string ${field}`);
    }
  }
  const maximumBodyBytes = validateBound(
    input.maxBodyBytes,
    DEFAULT_MAX_BODY_BYTES,
    MAXIMUM_BODY_BOUND,
    'maxBodyBytes',
    errors
  );
  const maximumDecisionCount = validateBound(
    input.maxDecisionCount,
    DEFAULT_MAX_DECISION_COUNT,
    MAXIMUM_DECISION_BOUND,
    'maxDecisionCount',
    errors
  );
  const normalizedEventHeadSha = typeof input.eventHeadSha === 'string'
    ? input.eventHeadSha.trim().toLowerCase()
    : '';
  const normalizedAuthorLogin = typeof input.eventAuthorLogin === 'string'
    ? input.eventAuthorLogin.trim().toLowerCase()
    : '';
  if (normalizedAuthorLogin && !LOGIN_PATTERN.test(normalizedAuthorLogin)) {
    addError(errors, 'invalid-event-author', '$.eventAuthorLogin', 'Event author login is malformed');
  }
  const metadata = {
    eventAuthorLogin: normalizedAuthorLogin,
    eventHeadSha: normalizedEventHeadSha
  };
  if (typeof input.body !== 'string') {
    return result({ errors, present: false, metadata, declaration: emptyDeclaration() });
  }
  if (
    input.body.length > maximumBodyBytes
    || byteLength(input.body) > maximumBodyBytes
  ) {
    addError(errors, 'oversized-body', '$.body', 'Pull request body exceeds the parser byte bound');
    return result({ errors, present: false, metadata, declaration: emptyDeclaration() });
  }

  const startCount = markerCount(input.body, START_MARKER);
  const endCount = markerCount(input.body, END_MARKER);
  const present = startCount > 0 || endCount > 0;
  if (startCount === 0 && endCount === 0) {
    return result({ errors, present: false, metadata, declaration: emptyDeclaration() });
  }
  if (startCount !== 1 || endCount !== 1) {
    addError(
      errors,
      'invalid-marker-count',
      '$.body',
      'Expected zero or one exact Product Impact decision marker pair'
    );
    return result({ errors, present, metadata, declaration: emptyDeclaration() });
  }
  const start = input.body.indexOf(START_MARKER);
  const end = input.body.indexOf(END_MARKER);
  if (start > end) {
    addError(errors, 'reversed-markers', '$.body', 'Product Impact decision markers are reversed');
    return result({ errors, present, metadata, declaration: emptyDeclaration() });
  }
  const block = input.body.slice(start + START_MARKER.length, end).trim();
  let parsed;
  try {
    parsed = JSON.parse(block);
  } catch (_error) {
    addError(errors, 'invalid-json', '$.decisionBlock', 'Product Impact decision block is not JSON');
    return result({ errors, present, metadata, declaration: emptyDeclaration() });
  }
  if (!isRecord(parsed)) {
    addError(errors, 'invalid-declaration', '$.decisionBlock', 'Expected a decision declaration object');
    return result({ errors, present, metadata, declaration: emptyDeclaration() });
  }
  exactFields(parsed, DECLARATION_FIELDS, '$.decisionBlock', errors);
  if (parsed.schemaVersion !== SCHEMA_VERSION) {
    addError(
      errors,
      'unsupported-schema-version',
      '$.decisionBlock.schemaVersion',
      `Expected schemaVersion ${SCHEMA_VERSION}`
    );
  }
  if (!Array.isArray(parsed.decisions)) {
    addError(errors, 'invalid-array', '$.decisionBlock.decisions', 'Expected a decisions array');
    return result({ errors, present, metadata, declaration: emptyDeclaration() });
  }
  if (parsed.decisions.length > maximumDecisionCount) {
    addError(
      errors,
      'too-many-decisions',
      '$.decisionBlock.decisions',
      'Decision count exceeds the parser bound'
    );
  }
  const normalizedHeadSha = typeof parsed.headSha === 'string'
    ? parsed.headSha.trim().toLowerCase()
    : '';
  if (parsed.decisions.length) {
    if (!HEAD_SHA_PATTERN.test(normalizedEventHeadSha)) {
      addError(errors, 'invalid-event-head-sha', '$.eventHeadSha', 'Event head SHA is malformed');
    }
    if (!HEAD_SHA_PATTERN.test(normalizedHeadSha)) {
      addError(
        errors,
        'invalid-declaration-head-sha',
        '$.decisionBlock.headSha',
        'Nonempty decisions require a 40-hex head SHA'
      );
    } else if (normalizedHeadSha !== normalizedEventHeadSha) {
      addError(
        errors,
        'stale-head-sha',
        '$.decisionBlock.headSha',
        'Decision head SHA does not match the event head SHA'
      );
    }
  } else if (normalizedHeadSha && !HEAD_SHA_PATTERN.test(normalizedHeadSha)) {
    addError(
      errors,
      'invalid-declaration-head-sha',
      '$.decisionBlock.headSha',
      'Head SHA must be empty or 40 hexadecimal characters'
    );
  }

  const normalizedDecisions = [];
  const seenDecisionScopes = new Set();
  let previousDecision = null;
  parsed.decisions.slice(0, maximumDecisionCount).forEach((decision, index) => {
    const path = `$.decisionBlock.decisions[${index}]`;
    if (!isRecord(decision)) {
      addError(errors, 'invalid-decision', path, 'Expected a current decision object');
      return;
    }
    exactFields(decision, DECISION_FIELDS, path, errors);
    const concernKey = canonicalString(
      decision.concernKey,
      `${path}.concernKey`,
      errors,
      { pattern: CONCERN_KEY_PATTERN }
    );
    const selectedOptionId = canonicalString(
      decision.selectedOptionId,
      `${path}.selectedOptionId`,
      errors,
      { pattern: OPTION_ID_PATTERN }
    );
    if (!Number.isInteger(decision.scenarioRevision) || decision.scenarioRevision < 1) {
      addError(
        errors,
        'invalid-scenario-revision',
        `${path}.scenarioRevision`,
        'Expected a positive integer scenario revision'
      );
    }
    if (decision.acceptedImpactClass !== 'AFFORDANCE_PRESERVED') {
      addError(
        errors,
        'invalid-current-impact-class',
        `${path}.acceptedImpactClass`,
        'Current decisions accept only AFFORDANCE_PRESERVED'
      );
    }
    const rationale = canonicalString(
      decision.rationale,
      `${path}.rationale`,
      errors,
      { maximumLength: MAXIMUM_TEXT_LENGTH, normalize: true }
    );
    const residualRisk = canonicalString(
      decision.residualRisk,
      `${path}.residualRisk`,
      errors,
      { maximumLength: MAXIMUM_TEXT_LENGTH, normalize: true }
    );
    const acknowledgedChangedRequirementRefs = sortedUniqueStrings(
      decision.acknowledgedChangedRequirementRefs,
      `${path}.acknowledgedChangedRequirementRefs`,
      errors,
      { pattern: REQUIREMENT_ID_PATTERN, maximumLength: MAXIMUM_ID_LENGTH }
    );
    const evidenceRefs = sortedUniqueStrings(
      decision.evidenceRefs,
      `${path}.evidenceRefs`,
      errors,
      { maximumLength: MAXIMUM_REF_LENGTH }
    );
    if (!evidenceRefs.length) {
      addError(errors, 'empty-evidence-refs', `${path}.evidenceRefs`, 'Expected evidence references');
    }
    const normalized = {
      concernKey,
      scenarioRevision: decision.scenarioRevision,
      selectedOptionId,
      acceptedImpactClass: decision.acceptedImpactClass,
      rationale,
      acknowledgedChangedRequirementRefs,
      evidenceRefs,
      residualRisk
    };
    if (concernKey && Number.isInteger(decision.scenarioRevision)) {
      const scope = `${concernKey}\u0000${decision.scenarioRevision}`;
      if (seenDecisionScopes.has(scope)) {
        addError(errors, 'duplicate-decision-scope', path, 'Decision scopes must be unique');
      }
      seenDecisionScopes.add(scope);
    }
    if (previousDecision && compareDecision(previousDecision, normalized) > 0) {
      addError(
        errors,
        'unsorted-decisions',
        '$.decisionBlock.decisions',
        'Decisions must use canonical concern/revision/option order'
      );
    }
    previousDecision = normalized;
    normalizedDecisions.push(normalized);
  });

  const declaration = {
    schemaVersion: SCHEMA_VERSION,
    headSha: normalizedHeadSha,
    decisions: normalizedDecisions
  };
  return result({ errors, present, metadata, declaration });
};
