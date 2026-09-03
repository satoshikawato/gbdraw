const SCHEMA_VERSION = 1;
const MAXIMUM_RULE_COUNT = 4;
const RULE_KINDS = Object.freeze([
  'single-canonical-entry-edge',
  'single-semantic-owner'
]);
const ENFORCEMENT_MODES = Object.freeze(['frozen', 'hard', 'report-only']);
const SUBJECT_CATEGORY_BY_KIND = Object.freeze({
  'single-canonical-entry-edge': 'canonical-entry-edge',
  'single-semantic-owner': 'definition-path'
});
const COMMON_RULE_FIELDS = Object.freeze([
  'baselineEligible',
  'detector',
  'enforcement',
  'key',
  'kind'
]);
const KIND_FIELDS = Object.freeze({
  'single-canonical-entry-edge': Object.freeze([
    'allowedEdges',
    'exactActiveEdgeCount'
  ]),
  'single-semantic-owner': Object.freeze([
    'allowedDefinitionPaths',
    'exactDefinitionCount'
  ])
});
const RULE_KEY_PATTERN = /^[a-z][a-z0-9]*(?:[.-][a-z0-9]+)*$/;
const DETECTOR_ID_PATTERN = /^[a-z][a-z0-9-]*(?:\.[a-z0-9-]+)*\.v[1-9][0-9]*$/;
const WEB_MODULE_PATH_PATTERN = /^[A-Za-z0-9_-]+(?:\/[A-Za-z0-9_.-]+)*\.[cm]?js$/;
const DIRECTION_ORDER = Object.freeze([
  'EXPANSION',
  'RELAXATION',
  'REINTERPRETATION',
  'CONTRACTION',
  'TIGHTENING'
]);
const RESULT_FIELDS = Object.freeze([
  'authorityResolution',
  'baselineRelation',
  'mode',
  'observation'
]);
const RESULT_OBSERVATIONS = Object.freeze([
  'ABSENT_REQUIRED',
  'CONFORMING',
  'DIVERGENT'
]);
const RESULT_MODES = Object.freeze(['FROZEN', 'HARD', 'REPORT_ONLY']);
const BASELINE_RELATIONS = Object.freeze([
  'ACCEPTED',
  'FIXED',
  'NEW',
  'NOT_APPLICABLE'
]);
const AUTHORITY_RESOLUTIONS = Object.freeze([
  'EXACT_CONTRACTION',
  'INVALID_CHANGE',
  'NOT_APPLICABLE',
  'RETAINED'
]);

const compareText = (left, right) => left < right ? -1 : left > right ? 1 : 0;
const isRecord = (value) => value !== null
  && typeof value === 'object'
  && !Array.isArray(value)
  && [Object.prototype, null].includes(Object.getPrototypeOf(value));
const stableUnique = (values) => [...new Set(values)].sort(compareText);
const edgeKey = ({ from, to }) => `${from}\u0000${to}`;
const compareEdges = (left, right) => compareText(left.from, right.from)
  || compareText(left.to, right.to);
const displayEdge = ({ from, to }) => `${from} -> ${to}`;
const ruleMap = (registry) => new Map((registry?.rules || []).map((rule) => [rule.key, rule]));
const fieldSet = (kind) => new Set([...COMMON_RULE_FIELDS, ...(KIND_FIELDS[kind] || [])]);
const stableDifference = (left, right, key = (value) => value) => {
  const rightKeys = new Set(right.map(key));
  return left.filter((value) => !rightKeys.has(key(value)));
};

const stableTextSet = (values, label) => {
  if (!Array.isArray(values) || values.some((value) => typeof value !== 'string')) {
    throw new Error(`${label} must be an array of strings`);
  }
  return Object.freeze(stableUnique(values));
};

const SUBJECT_DELTA_FIELDS = Object.freeze([
  'afterSubjects',
  'beforeSubjects',
  'detector',
  'kind',
  'ruleKey',
  'subjectCategory'
]);
const CANONICAL_IDENTIFIER_PATTERN = /^[a-z][a-z0-9]*(?:[.-][a-z0-9]+)*$/;

const canonicalSubjectSet = (values, label) => {
  const subjects = stableTextSet(values, label);
  subjects.forEach((subject) => {
    if (!subject || subject !== subject.trim() || /[\u0000-\u001f\u007f]/.test(subject)) {
      throw new Error(`${label} must contain canonical nonempty subject strings`);
    }
  });
  return subjects;
};

export const createArchitectureSubjectDelta = (input) => {
  if (!isRecord(input)) {
    throw new Error('Architecture subject delta requires one input object');
  }
  const fields = Object.keys(input).sort(compareText);
  if (
    fields.length !== SUBJECT_DELTA_FIELDS.length
    || fields.some((field, index) => field !== SUBJECT_DELTA_FIELDS[index])
  ) {
    throw new Error(
      `Architecture subject delta requires exactly: ${SUBJECT_DELTA_FIELDS.join(', ')}`
    );
  }
  if (typeof input.ruleKey !== 'string' || !RULE_KEY_PATTERN.test(input.ruleKey)) {
    throw new Error('Architecture subject delta ruleKey must be a canonical rule key');
  }
  if (typeof input.detector !== 'string' || !DETECTOR_ID_PATTERN.test(input.detector)) {
    throw new Error('Architecture subject delta detector must be a versioned detector ID');
  }
  for (const field of ['kind', 'subjectCategory']) {
    if (
      typeof input[field] !== 'string'
      || !CANONICAL_IDENTIFIER_PATTERN.test(input[field])
    ) {
      throw new Error(`Architecture subject delta ${field} must be a canonical identifier`);
    }
  }
  const beforeSubjects = canonicalSubjectSet(
    input.beforeSubjects,
    'Architecture subject delta beforeSubjects'
  );
  const afterSubjects = canonicalSubjectSet(
    input.afterSubjects,
    'Architecture subject delta afterSubjects'
  );
  const addedSubjects = Object.freeze(stableDifference(afterSubjects, beforeSubjects));
  const removedSubjects = Object.freeze(stableDifference(beforeSubjects, afterSubjects));
  return Object.freeze({
    ruleKey: input.ruleKey,
    kind: input.kind,
    detector: input.detector,
    subjectCategory: input.subjectCategory,
    beforeSubjects,
    addedSubjects,
    removedSubjects,
    afterSubjects,
    changed: addedSubjects.length > 0 || removedSubjects.length > 0
  });
};

export const summarizeArchitectureInventory = (beforeValues, afterValues) => {
  const before = stableTextSet(beforeValues, 'Architecture inventory before');
  const after = stableTextSet(afterValues, 'Architecture inventory after');
  const added = Object.freeze(stableDifference(after, before));
  const removed = Object.freeze(stableDifference(before, after));
  return Object.freeze({
    before,
    added,
    removed,
    after,
    delta: after.length - before.length
  });
};

const validateGraphFields = (graph) => {
  const fields = Object.keys(graph).sort(compareText);
  if (fields.length !== 2 || fields[0] !== 'edges' || fields[1] !== 'nodes') {
    throw new Error('Directed graph requires exactly: edges, nodes');
  }
};

const normalizeDirectedGraph = (graph) => {
  if (!isRecord(graph)) throw new Error('Directed graph must be an object');
  validateGraphFields(graph);
  const nodes = stableTextSet(graph.nodes, 'Directed graph nodes');
  const nodeSet = new Set(nodes);
  if (!Array.isArray(graph.edges)) throw new Error('Directed graph edges must be an array');
  const edgesByKey = new Map();
  graph.edges.forEach((edge, index) => {
    if (!isRecord(edge) || typeof edge.from !== 'string' || typeof edge.to !== 'string') {
      throw new Error(`Directed graph edge ${index} must define string from and to fields`);
    }
    if (Object.keys(edge).some((field) => !['from', 'to'].includes(field))) {
      throw new Error(`Directed graph edge ${index} contains an unknown field`);
    }
    if (!nodeSet.has(edge.from) || !nodeSet.has(edge.to)) {
      throw new Error(`Directed graph edge ${displayEdge(edge)} references an unknown node`);
    }
    edgesByKey.set(edgeKey(edge), Object.freeze({ from: edge.from, to: edge.to }));
  });
  return Object.freeze({
    nodes,
    edges: Object.freeze([...edgesByKey.values()].sort(compareEdges))
  });
};

export const findDirectedGraphCycles = (graph) => {
  const normalized = normalizeDirectedGraph(graph);
  const adjacency = new Map(normalized.nodes.map((node) => [node, []]));
  normalized.edges.forEach(({ from, to }) => adjacency.get(from).push(to));
  adjacency.forEach((targets) => targets.sort(compareText));

  let nextIndex = 0;
  const indexByNode = new Map();
  const lowLinkByNode = new Map();
  const stack = [];
  const stacked = new Set();
  const components = [];

  const visit = (node) => {
    indexByNode.set(node, nextIndex);
    lowLinkByNode.set(node, nextIndex);
    nextIndex += 1;
    stack.push(node);
    stacked.add(node);

    adjacency.get(node).forEach((target) => {
      if (!indexByNode.has(target)) {
        visit(target);
        lowLinkByNode.set(
          node,
          Math.min(lowLinkByNode.get(node), lowLinkByNode.get(target))
        );
      } else if (stacked.has(target)) {
        lowLinkByNode.set(
          node,
          Math.min(lowLinkByNode.get(node), indexByNode.get(target))
        );
      }
    });

    if (lowLinkByNode.get(node) !== indexByNode.get(node)) return;
    const component = [];
    let member;
    do {
      member = stack.pop();
      stacked.delete(member);
      component.push(member);
    } while (member !== node);
    components.push(component.sort(compareText));
  };

  normalized.nodes.forEach((node) => {
    if (!indexByNode.has(node)) visit(node);
  });

  const cycles = components.flatMap((nodes) => {
    const members = new Set(nodes);
    const internalEdges = normalized.edges.filter(({ from, to }) => (
      members.has(from) && members.has(to)
    ));
    const cyclic = nodes.length > 1
      || internalEdges.some(({ from, to }) => from === to);
    if (!cyclic) return [];
    const edgeSubjects = internalEdges.map(displayEdge);
    const subject = `nodes=[${nodes.join(', ')}]; edges=[${edgeSubjects.join(', ')}]`;
    return [Object.freeze({
      subject,
      nodes: Object.freeze([...nodes]),
      internalEdges: Object.freeze([...internalEdges])
    })];
  }).sort((left, right) => compareText(left.subject, right.subject));

  return Object.freeze({
    nodeCount: normalized.nodes.length,
    edgeCount: normalized.edges.length,
    cycles: Object.freeze(cycles)
  });
};

const addError = (errors, code, path, message) => {
  errors.push(Object.freeze({ code, path, message }));
};

const validateExactFields = (value, expectedFields, path, errors) => {
  Object.keys(value).sort(compareText).forEach((field) => {
    if (!expectedFields.has(field)) {
      addError(errors, 'unknown-field', `${path}.${field}`, `Unknown field ${field}`);
    }
  });
  [...expectedFields].sort(compareText).forEach((field) => {
    if (!Object.hasOwn(value, field)) {
      addError(errors, 'missing-field', `${path}.${field}`, `Missing required field ${field}`);
    }
  });
};

const isWebModulePath = (value) => typeof value === 'string'
  && WEB_MODULE_PATH_PATTERN.test(value)
  && !value.startsWith('gbdraw/web/js/')
  && !value.includes('//')
  && !value.split('/').some((segment) => segment === '.' || segment === '..');

const validateWebModulePath = (value, path, errors) => {
  if (!isWebModulePath(value)) {
    addError(
      errors,
      'invalid-web-module-path',
      path,
      'Expected a normalized JavaScript module path relative to gbdraw/web/js/'
    );
  }
};

const validateSortedUniqueStrings = (values, path, errors) => {
  if (!Array.isArray(values)) {
    addError(errors, 'invalid-array', path, 'Expected an array');
    return;
  }
  const seen = new Set();
  let previous = null;
  values.forEach((value, index) => {
    validateWebModulePath(value, `${path}[${index}]`, errors);
    if (typeof value !== 'string') return;
    if (seen.has(value)) {
      addError(errors, 'duplicate-value', `${path}[${index}]`, `Duplicate value ${value}`);
    }
    if (previous !== null && compareText(previous, value) > 0) {
      addError(errors, 'unsorted-array', path, 'Values must use ascending canonical order');
    }
    seen.add(value);
    previous = value;
  });
};

const validateSortedUniqueEdges = (edges, path, errors) => {
  if (!Array.isArray(edges)) {
    addError(errors, 'invalid-array', path, 'Expected an array');
    return;
  }
  if (!edges.length) {
    addError(errors, 'empty-array', path, 'At least one allowed edge is required');
  }
  const seen = new Set();
  let previous = null;
  edges.forEach((edge, index) => {
    const edgePath = `${path}[${index}]`;
    if (!isRecord(edge)) {
      addError(errors, 'invalid-edge', edgePath, 'Expected an edge object');
      return;
    }
    validateExactFields(edge, new Set(['from', 'to']), edgePath, errors);
    validateWebModulePath(edge.from, `${edgePath}.from`, errors);
    validateWebModulePath(edge.to, `${edgePath}.to`, errors);
    if (typeof edge.from !== 'string' || typeof edge.to !== 'string') return;
    const key = edgeKey(edge);
    if (seen.has(key)) {
      addError(errors, 'duplicate-edge', edgePath, `Duplicate edge ${displayEdge(edge)}`);
    }
    if (previous !== null && compareEdges(previous, edge) > 0) {
      addError(errors, 'unsorted-edges', path, 'Edges must use canonical (from, to) order');
    }
    seen.add(key);
    previous = edge;
  });
};

const validateStableSubjectEncoder = (rule, detector, path, errors) => {
  if (!detector || typeof detector.encodeSubject !== 'function') {
    addError(
      errors,
      'missing-subject-encoder',
      `${path}.detector`,
      `Detector ${rule.detector} must provide a stable subject encoder`
    );
    return;
  }
  const definitionPaths = Array.isArray(rule.allowedDefinitionPaths)
    ? rule.allowedDefinitionPaths
    : [];
  const edges = Array.isArray(rule.allowedEdges) ? rule.allowedEdges : [];
  const candidates = rule.kind === 'single-semantic-owner'
    ? definitionPaths.map((modulePath) => ({ path: modulePath }))
    : edges.filter(isRecord).map(({ from, to }) => ({ from, to }));
  candidates.forEach((candidate, index) => {
    try {
      const frozenCandidate = Object.freeze({ ...candidate });
      const first = detector.encodeSubject(frozenCandidate);
      const second = detector.encodeSubject(frozenCandidate);
      if (typeof first !== 'string' || !first || first !== second) {
        addError(
          errors,
          'unstable-subject-encoder',
          `${path}.detector`,
          `Detector ${rule.detector} did not return a stable nonempty subject for item ${index}`
        );
      }
    } catch (error) {
      addError(
        errors,
        'invalid-subject-encoder',
        `${path}.detector`,
        `Detector ${rule.detector} subject encoder failed: ${error.message}`
      );
    }
  });
};

export const validateArchitectureRuleRegistry = (
  registry,
  detectorCatalog,
  { availableEnforcements = ENFORCEMENT_MODES } = {}
) => {
  const errors = [];
  if (!isRecord(registry)) {
    addError(errors, 'invalid-registry', '$', 'Expected a rule registry object');
    return Object.freeze({ valid: false, errors: Object.freeze(errors) });
  }
  validateExactFields(registry, new Set(['rules', 'schemaVersion']), '$', errors);
  if (registry.schemaVersion !== SCHEMA_VERSION) {
    addError(
      errors,
      'unsupported-schema-version',
      '$.schemaVersion',
      `Expected schemaVersion ${SCHEMA_VERSION}`
    );
  }
  if (!Array.isArray(registry.rules)) {
    addError(errors, 'invalid-rules', '$.rules', 'Expected a rules array');
    return Object.freeze({ valid: false, errors: Object.freeze(errors) });
  }
  if (registry.rules.length > MAXIMUM_RULE_COUNT) {
    addError(
      errors,
      'too-many-rules',
      '$.rules',
      `Schema version ${SCHEMA_VERSION} permits at most ${MAXIMUM_RULE_COUNT} rules`
    );
  }

  const available = new Set(availableEnforcements);
  const seenKeys = new Set();
  let previousKey = null;
  registry.rules.forEach((rule, index) => {
    const path = `$.rules[${index}]`;
    if (!isRecord(rule)) {
      addError(errors, 'invalid-rule', path, 'Expected a rule object');
      return;
    }
    const kind = rule.kind;
    validateExactFields(rule, fieldSet(kind), path, errors);
    if (typeof rule.key !== 'string' || !RULE_KEY_PATTERN.test(rule.key)) {
      addError(errors, 'invalid-rule-key', `${path}.key`, 'Expected a canonical rule key');
    } else {
      if (seenKeys.has(rule.key)) {
        addError(errors, 'duplicate-rule-key', `${path}.key`, `Duplicate rule key ${rule.key}`);
      }
      if (previousKey !== null && compareText(previousKey, rule.key) > 0) {
        addError(errors, 'unsorted-rules', '$.rules', 'Rules must use ascending key order');
      }
      seenKeys.add(rule.key);
      previousKey = rule.key;
    }
    if (!RULE_KINDS.includes(kind)) {
      addError(errors, 'unsupported-rule-kind', `${path}.kind`, `Unsupported rule kind ${kind}`);
    }
    if (typeof rule.detector !== 'string' || !DETECTOR_ID_PATTERN.test(rule.detector)) {
      addError(errors, 'invalid-detector-id', `${path}.detector`, 'Expected a versioned detector ID');
    }
    if (!ENFORCEMENT_MODES.includes(rule.enforcement)) {
      addError(
        errors,
        'unsupported-enforcement',
        `${path}.enforcement`,
        `Unsupported enforcement mode ${rule.enforcement}`
      );
    } else if (!available.has(rule.enforcement)) {
      addError(
        errors,
        'unavailable-enforcement',
        `${path}.enforcement`,
        `Enforcement mode ${rule.enforcement} is recognized but unavailable`
      );
    }
    if (typeof rule.baselineEligible !== 'boolean') {
      addError(errors, 'invalid-baseline-eligibility', `${path}.baselineEligible`, 'Expected a boolean');
    }
    if (rule.enforcement === 'hard' && rule.baselineEligible !== false) {
      addError(
        errors,
        'invalid-hard-eligibility',
        `${path}.baselineEligible`,
        'Hard rules cannot be baseline eligible'
      );
    }
    if (rule.enforcement === 'frozen' && rule.baselineEligible !== true) {
      addError(
        errors,
        'invalid-frozen-eligibility',
        `${path}.baselineEligible`,
        'Frozen rules must be baseline eligible'
      );
    }

    if (kind === 'single-semantic-owner') {
      validateSortedUniqueStrings(rule.allowedDefinitionPaths, `${path}.allowedDefinitionPaths`, errors);
      if (rule.exactDefinitionCount !== 1) {
        addError(
          errors,
          'invalid-definition-count',
          `${path}.exactDefinitionCount`,
          'Schema version 1 requires exactDefinitionCount: 1'
        );
      }
    } else if (kind === 'single-canonical-entry-edge') {
      validateSortedUniqueEdges(rule.allowedEdges, `${path}.allowedEdges`, errors);
      if (rule.exactActiveEdgeCount !== 1) {
        addError(
          errors,
          'invalid-active-edge-count',
          `${path}.exactActiveEdgeCount`,
          'Schema version 1 requires exactActiveEdgeCount: 1'
        );
      }
    }

    const detector = detectorCatalog?.[rule.detector];
    if (!detector) {
      addError(
        errors,
        'unknown-detector',
        `${path}.detector`,
        `Detector ${rule.detector} is not present in the trusted catalog`
      );
      return;
    }
    const expectedCategory = SUBJECT_CATEGORY_BY_KIND[kind];
    if (expectedCategory && detector.subjectCategory !== expectedCategory) {
      addError(
        errors,
        'detector-kind-mismatch',
        `${path}.detector`,
        `Rule kind ${kind} requires detector subject category ${expectedCategory}`
      );
    }
    if (typeof detector.detect !== 'function') {
      addError(
        errors,
        'missing-detector-function',
        `${path}.detector`,
        `Detector ${rule.detector} must provide source-fact extraction`
      );
    }
    validateStableSubjectEncoder(rule, detector, path, errors);
  });

  return Object.freeze({
    valid: errors.length === 0,
    errors: Object.freeze(errors)
  });
};

export const classifyArchitectureRuleObservation = (rule, detectorOutput) => {
  if (rule.kind === 'single-semantic-owner') {
    const definitions = Array.isArray(detectorOutput?.observedDefinitions)
      ? detectorOutput.observedDefinitions
      : [];
    const definitionCount = Number.isInteger(detectorOutput?.definitionCount)
      ? detectorOutput.definitionCount
      : definitions.reduce((total, definition) => total + Number(definition.count || 0), 0);
    const allowed = new Set(rule.allowedDefinitionPaths);
    const unapprovedSubjects = stableUnique(definitions
      .filter(({ path }) => !allowed.has(path))
      .map(({ subject, path }) => subject || path));
    const subjects = stableUnique(definitions.map(({ subject, path }) => subject || path));
    const observation = definitionCount === 0
      ? 'ABSENT_REQUIRED'
      : definitionCount !== rule.exactDefinitionCount || unapprovedSubjects.length
        ? 'DIVERGENT'
        : 'CONFORMING';
    return Object.freeze({
      observation,
      subjects: Object.freeze(subjects),
      unapprovedSubjects: Object.freeze(unapprovedSubjects),
      observedCount: definitionCount,
      requiredCount: rule.exactDefinitionCount
    });
  }

  if (rule.kind === 'single-canonical-entry-edge') {
    const edges = Array.isArray(detectorOutput?.observedEdges)
      ? detectorOutput.observedEdges
      : [];
    const allowed = new Set(rule.allowedEdges.map(edgeKey));
    const subjects = stableUnique(edges.map((edge) => edge.subject || displayEdge(edge)));
    const unapprovedSubjects = stableUnique(edges
      .filter((edge) => !allowed.has(edgeKey(edge)))
      .map((edge) => edge.subject || displayEdge(edge)));
    const observation = edges.length === 0
      ? 'ABSENT_REQUIRED'
      : edges.length !== rule.exactActiveEdgeCount || unapprovedSubjects.length
        ? 'DIVERGENT'
        : 'CONFORMING';
    return Object.freeze({
      observation,
      subjects: Object.freeze(subjects),
      unapprovedSubjects: Object.freeze(unapprovedSubjects),
      observedCount: edges.length,
      requiredCount: rule.exactActiveEdgeCount
    });
  }

  throw new Error(`Unsupported rule kind ${rule.kind}`);
};

export const evaluateArchitectureRuleResult = (semanticInputs) => {
  if (!isRecord(semanticInputs)) {
    throw new Error('Architecture rule evaluation requires one semantic input object');
  }
  const actualFields = Object.keys(semanticInputs).sort(compareText);
  if (
    actualFields.length !== RESULT_FIELDS.length
    || actualFields.some((field, index) => field !== RESULT_FIELDS[index])
  ) {
    throw new Error(
      `Architecture rule evaluation requires exactly: ${RESULT_FIELDS.join(', ')}`
    );
  }

  const {
    observation,
    mode,
    baselineRelation,
    authorityResolution
  } = semanticInputs;
  const supportedValues = [
    ['observation', observation, RESULT_OBSERVATIONS],
    ['mode', mode, RESULT_MODES],
    ['baselineRelation', baselineRelation, BASELINE_RELATIONS],
    ['authorityResolution', authorityResolution, AUTHORITY_RESOLUTIONS]
  ];
  supportedValues.forEach(([field, value, allowed]) => {
    if (!allowed.includes(value)) {
      throw new Error(`Unsupported architecture rule ${field}: ${value}`);
    }
  });

  if (mode === 'FROZEN') {
    throw new Error('Architecture rule mode FROZEN is unavailable');
  }
  if (baselineRelation !== 'NOT_APPLICABLE') {
    throw new Error(`${mode} requires baselineRelation NOT_APPLICABLE`);
  }
  if (authorityResolution !== 'NOT_APPLICABLE') {
    throw new Error(`${mode} requires authorityResolution NOT_APPLICABLE`);
  }

  const decision = observation === 'CONFORMING'
    ? 'PASS'
    : mode === 'HARD' ? 'FAIL' : 'REPORT';
  return Object.freeze({
    observation,
    mode,
    baselineRelation,
    authorityResolution,
    decision
  });
};

const enforcementDirection = (before, after) => {
  if (before === after) return null;
  if (
    (before === 'hard' && ['frozen', 'report-only'].includes(after))
    || (before === 'frozen' && after === 'report-only')
  ) return 'RELAXATION';
  if (
    (before === 'report-only' && ['frozen', 'hard'].includes(after))
    || (before === 'frozen' && after === 'hard')
  ) return 'TIGHTENING';
  throw new Error(`Unsupported enforcement transition ${before} -> ${after}`);
};

const addChange = (changes, rule, field, direction, before, after) => {
  changes.push(Object.freeze({ rule, field, direction, before, after }));
};

export const classifyArchitectureAuthorityDelta = (baseRegistry, candidateRegistry) => {
  const beforeByKey = ruleMap(baseRegistry);
  const afterByKey = ruleMap(candidateRegistry);
  const keys = stableUnique([...beforeByKey.keys(), ...afterByKey.keys()]);
  const changes = [];

  keys.forEach((key) => {
    const before = beforeByKey.get(key);
    const after = afterByKey.get(key);
    if (!before) {
      addChange(changes, key, 'rule', 'EXPANSION', null, after);
      return;
    }
    if (!after) {
      addChange(changes, key, 'rule', 'RELAXATION', before, null);
      return;
    }
    if (before.kind !== after.kind) {
      addChange(changes, key, 'kind', 'REINTERPRETATION', before.kind, after.kind);
    }
    if (before.detector !== after.detector) {
      addChange(
        changes,
        key,
        'detector',
        'REINTERPRETATION',
        before.detector,
        after.detector
      );
    }
    const enforcement = enforcementDirection(before.enforcement, after.enforcement);
    if (enforcement) {
      addChange(
        changes,
        key,
        'enforcement',
        enforcement,
        before.enforcement,
        after.enforcement
      );
    }
    if (before.baselineEligible !== after.baselineEligible) {
      addChange(
        changes,
        key,
        'baselineEligible',
        after.baselineEligible ? 'RELAXATION' : 'TIGHTENING',
        before.baselineEligible,
        after.baselineEligible
      );
    }
    if (before.kind !== after.kind) return;
    if (after.kind === 'single-semantic-owner') {
      stableDifference(
        after.allowedDefinitionPaths,
        before.allowedDefinitionPaths
      ).forEach((value) => {
        addChange(changes, key, 'allowedDefinitionPaths', 'EXPANSION', null, value);
      });
      stableDifference(
        before.allowedDefinitionPaths,
        after.allowedDefinitionPaths
      ).forEach((value) => {
        addChange(changes, key, 'allowedDefinitionPaths', 'CONTRACTION', value, null);
      });
    } else if (after.kind === 'single-canonical-entry-edge') {
      stableDifference(after.allowedEdges, before.allowedEdges, edgeKey)
        .sort(compareEdges)
        .forEach((value) => {
          addChange(changes, key, 'allowedEdges', 'EXPANSION', null, value);
        });
      stableDifference(before.allowedEdges, after.allowedEdges, edgeKey)
        .sort(compareEdges)
        .forEach((value) => {
          addChange(changes, key, 'allowedEdges', 'CONTRACTION', value, null);
        });
    }
  });

  const directions = DIRECTION_ORDER.filter((direction) => (
    changes.some((change) => change.direction === direction)
  ));
  return Object.freeze({
    classification: directions.length === 0
      ? 'UNCHANGED'
      : directions.length === 1 ? directions[0] : 'MIXED',
    directions: Object.freeze(directions),
    changes: Object.freeze(changes)
  });
};

export const WEB_ARCHITECTURE_RULE_SCHEMA = Object.freeze({
  schemaVersion: SCHEMA_VERSION,
  maximumRuleCount: MAXIMUM_RULE_COUNT,
  kinds: RULE_KINDS,
  enforcementModes: ENFORCEMENT_MODES
});
