// Owner-graph detectors for gbdraw/web/js (implementation plan Phase B1,
// docs/internal/WEB_OWNER_COUPLING_PREVENTION_IMPLEMENTATION_PLAN_2026-10-05.md).
//
// They observe the coupling channels the static import graph does not see:
// owner objects injected at composition, closures bound to owners created
// later, functions assigned into `state`, factories that take whole owner
// objects, the call shapes of a projection domain outside its owner, and the
// trigger sites of a heavy derived value. Every detector reads source text
// only; nothing is executed.
//
// Keep each detector ID and its transitive helpers executable-identical while
// referenced. Add a versioned detector beside it, then migrate authority in a
// separate pull request.
import { maskJavaScript } from './web-change-source.mjs';

const WEB_SOURCE_PREFIX = 'gbdraw/web/js/';
const normalizeModulePath = (path) => String(path || '')
  .replaceAll('\\', '/')
  .replace(new RegExp(`^${WEB_SOURCE_PREFIX}`), '');
const sourceEntries = (sources) => (
  sources instanceof Map ? [...sources] : Object.entries(sources || {})
).map(([path, source]) => [normalizeModulePath(path), String(source || '')])
  .sort(([left], [right]) => left.localeCompare(right));
const unique = (values) => [...new Set(values)].sort();

// Default registry, used until tools/web-owner-graph.json exists (Phase B2).
// `ownerObjectNames` lists injected owner objects by name; factory results in
// a composition root are owner objects as well.
export const WEB_OWNER_GRAPH_DEFAULTS = Object.freeze({
  schemaVersion: 1,
  compositionRoots: Object.freeze([
    'app/app-setup.js',
    'app/feature-editor.js',
    'app/legend.js',
    'app/legend-layout.js'
  ]),
  ownerObjectNames: Object.freeze([
    'history', 'historySnapshots', 'historyFileStore', 'previewRuntime', 'rulePreparation',
    'legendLayout', 'legendActions', 'svgActions', 'featureActions', 'featureSvgActions',
    'labelActions', 'visibilityActions', 'colorActions', 'ruleActions', 'placementActions',
    'featureEditTableActions', 'featureSelection', 'previewTransformInteraction',
    'recordDisplayControls', 'previewFeatureSearch', 'resultsManager', 'rightDrawerActions',
    'entryActions', 'layoutActions', 'strokeActions', 'dragActions', 'sortActions',
    'diagramActions', 'repositionActions', 'canvasActions', 'trackLayoutActions',
    'circularTrackSlotEditor', 'linearTrackSlotEditor', 'annotationEditor',
    'orthogroupActions', 'similarityAlignmentActions', 'featureRecordRotation',
    'linearRecordSelector', 'linearTypography', 'losatSettings', 'autoValueDisplay', 'fileStore'
  ]),
  ownerObjectNamePattern: '(?:Actions|Runtime|Layout|Preparation|Selection|Interaction|Snapshots|Controls|Manager|Store|Editor|Rotation|Selector)$|^history$',
  stateModule: 'state.js',
  projectionDomains: Object.freeze([
    Object.freeze({
      name: 'feature-visibility',
      owners: Object.freeze(['app/feature-editor/visibility-actions.js']),
      functions: Object.freeze(['projectFeatureVisibility', 'reconcileFeatureVisibility'])
    }),
    Object.freeze({
      name: 'feature-visibility-labels',
      owners: Object.freeze(['app/feature-editor/label-actions.js']),
      functions: Object.freeze(['applyFeatureVisibilityToLabels'])
    }),
    Object.freeze({
      name: 'label-intent',
      owners: Object.freeze(['app/feature-editor/label-actions.js']),
      functions: Object.freeze(['reconcileLabelOverrides', 'projectLabelIntent'])
    }),
    Object.freeze({
      name: 'legend-order',
      owners: Object.freeze(['app/legend/utils.js', 'app/legend/entry-actions.js']),
      functions: Object.freeze(['orderLegendEntries', 'reconcileLegendEntries'])
    }),
    Object.freeze({
      name: 'palette-rules',
      owners: Object.freeze(['app/svg-styles.js']),
      functions: Object.freeze(['applyPaletteToSvg', 'applySpecificRulesToSvg'])
    }),
    Object.freeze({
      name: 'strokes',
      owners: Object.freeze(['app/legend/stroke-actions.js']),
      functions: Object.freeze(['reconcileStrokeOverrides', 'applyStrokeOverridesToSvg'])
    }),
    Object.freeze({
      name: 'composition',
      owners: Object.freeze(['app/legend-layout.js']),
      functions: Object.freeze(['reconcileCompositionUserDeltas'])
    })
  ]),
  heavyProducers: Object.freeze([
    Object.freeze({
      name: 'rulePreparation',
      owner: 'app/rule-matching.js',
      methods: Object.freeze(['prepare', 'prepareDrawn', 'prepareVisibility', 'prepareCandidate', 'evaluate', 'run'])
    })
  ])
});

const resolveRegistry = (registry) => ({ ...WEB_OWNER_GRAPH_DEFAULTS, ...(registry || {}) });

// --- text helpers -----------------------------------------------------------

const matchingClose = (code, openIndex) => {
  const pairs = { '(': ')', '[': ']', '{': '}' };
  const stack = [];
  for (let index = openIndex; index < code.length; index += 1) {
    const character = code[index];
    if (pairs[character]) stack.push(pairs[character]);
    else if (character === ')' || character === ']' || character === '}') {
      if (stack.pop() !== character) return -1;
      if (!stack.length) return index;
    }
  }
  return -1;
};
const lineNumberAt = (source, index) => source.slice(0, index).split('\n').length;
const lineOf = (source, index) => {
  const start = source.lastIndexOf('\n', index - 1) + 1;
  const end = source.indexOf('\n', index);
  return source.slice(start, end < 0 ? source.length : end);
};
const escapeRegExp = (value) => String(value).replace(/[.*+?^${}()|[\]\\]/g, '\\$&');
const identifierPattern = (names) => names.map(escapeRegExp).join('|');

// Top-level argument slices of a `{ ... }` object literal: `key: value` pairs
// and shorthand identifiers, split at depth-0 commas.
const objectLiteralEntries = (code, source, openIndex) => {
  const closeIndex = matchingClose(code, openIndex);
  if (closeIndex < 0) return [];
  const entries = [];
  let depth = 0;
  let start = openIndex + 1;
  for (let index = openIndex + 1; index <= closeIndex; index += 1) {
    const character = code[index];
    if (character === '(' || character === '[' || character === '{') depth += 1;
    else if (character === ')' || character === ']' || character === '}') {
      if (index === closeIndex) {
        if (index > start) entries.push({ start, end: index });
        break;
      }
      depth -= 1;
    } else if (character === ',' && depth === 0) {
      entries.push({ start, end: index });
      start = index + 1;
    }
  }
  return entries
    .map(({ start: entryStart, end }) => ({ start: entryStart, end, text: source.slice(entryStart, end), code: code.slice(entryStart, end) }))
    .filter((entry) => entry.code.trim());
};

// Factory calls in a module: `const x = createX(` , `let x = createX(`, and a
// deferred `x = createX(` after `let x = null`.
const factoryCalls = (code, source) => {
  const calls = [];
  const pattern = /(?:^|[;{}\n])\s*(?:(?:const|let)\s+)?(\w+)\s*=\s*(create\w+)\s*\(/g;
  for (const match of code.matchAll(pattern)) {
    const openIndex = match.index + match[0].length - 1;
    const closeIndex = matchingClose(code, openIndex);
    if (closeIndex < 0) continue;
    calls.push({
      consumer: match[1],
      factory: match[2],
      index: match.index + match[0].indexOf(match[1]),
      openIndex,
      closeIndex,
      line: lineNumberAt(source, match.index + match[0].indexOf(match[1]))
    });
  }
  return calls;
};

// An owner object is a registered name, or a factory result of a composition
// root whose name matches the owner-name pattern (a string option such as
// `trackLayout` or a number such as `maxActions` is not one).
const compositionFactoryNames = (sources, registry) => {
  const roots = new Set(registry.compositionRoots);
  const names = new Set();
  for (const [path, source] of sourceEntries(sources)) {
    if (!roots.has(path)) continue;
    factoryCalls(maskJavaScript(source), source).forEach((call) => names.add(call.consumer));
  }
  return names;
};
const ownerNameTester = (registry, factoryNames) => {
  const registered = new Set(registry.ownerObjectNames || []);
  const pattern = registry.ownerObjectNamePattern ? new RegExp(registry.ownerObjectNamePattern) : null;
  return (name) => registered.has(name) || Boolean(pattern && pattern.test(name) && factoryNames.has(name));
};

// --- detectors ----------------------------------------------------------------

// Injection edges and forward closures in the composition roots.
const detectComposition = (sources, registryInput) => {
  const registry = resolveRegistry(registryInput);
  const roots = new Set(registry.compositionRoots);
  const injections = [];
  const forwardClosures = [];
  const shims = [];
  for (const [path, source] of sourceEntries(sources)) {
    if (!roots.has(path)) continue;
    const code = maskJavaScript(source);
    const calls = factoryCalls(code, source);
    const assignedAt = new Map();
    calls.forEach((call) => {
      if (!assignedAt.has(call.consumer)) assignedAt.set(call.consumer, call.index);
    });
    // Inside a root, every factory result is an owner object when it is
    // passed to another factory, whatever its name.
    const registered = new Set(registry.ownerObjectNames || []);
    const isOwner = (name) => registered.has(name) || assignedAt.has(name);
    calls.forEach((call) => {
      const argumentsCode = code.slice(call.openIndex + 1, call.closeIndex);
      const argumentsText = source.slice(call.openIndex + 1, call.closeIndex);
      // Arguments may be an object literal, bare identifiers, or both.
      const objectOpen = argumentsCode.search(/\{/);
      const entries = objectOpen >= 0
        ? objectLiteralEntries(argumentsCode, argumentsText, objectOpen)
        : [];
      const bareArguments = (objectOpen >= 0 ? argumentsCode.slice(0, objectOpen) : argumentsCode)
        .split(',').map((part) => part.trim()).filter(Boolean);
      const providers = new Set();
      bareArguments.forEach((name) => {
        if (/^\w+$/.test(name) && name !== call.consumer && isOwner(name)) providers.add(name);
      });
      entries.forEach((entry) => {
        const trimmed = entry.code.trim();
        const shorthand = trimmed.match(/^(\w+)$/);
        const pair = trimmed.match(/^(\w+)\s*:\s*(\w+)$/);
        const provider = shorthand ? shorthand[1] : pair ? pair[2] : null;
        if (provider && provider !== call.consumer && isOwner(provider)) providers.add(provider);
        const closure = trimmed.match(/^(\w+)\s*:\s*(?:async\s*)?(?:\([^)]*\)|\w+)\s*=>/)
          || trimmed.match(/^(?:async\s+)?(\w+)\s*\([^)]*\)\s*\{/);
        if (!closure) return;
        const key = closure[1];
        const bodyCode = entry.code;
        for (const reference of bodyCode.matchAll(/\b(\w+)\s*\??\.\s*(\w+)/g)) {
          const [, provider2, method] = reference;
          if (provider2 === call.consumer || !assignedAt.has(provider2)) continue;
          const subject = `${path}|${call.consumer}.${key}->${provider2}.${method}`;
          const record = { path, consumer: call.consumer, port: key, provider: provider2, method, line: call.line, subject };
          shims.push(record);
          if (assignedAt.get(provider2) > call.index) forwardClosures.push(record);
        }
      });
      [...providers].sort().forEach((provider) => {
        injections.push({
          path, consumer: call.consumer, provider, line: call.line,
          subject: `${path}|${call.consumer}<-${provider}`
        });
      });
    });
  }
  const dedupe = (records) => {
    const seen = new Map();
    records.forEach((record) => { if (!seen.has(record.subject)) seen.set(record.subject, record); });
    return [...seen.values()].sort((left, right) => left.subject.localeCompare(right.subject));
  };
  return { injections: dedupe(injections), forwardClosures: dedupe(forwardClosures), shims: dedupe(shims) };
};

const detectInjectionEdgesV1 = (sources, registry) => {
  const { injections } = detectComposition(sources, registry);
  return Object.freeze({
    observedEdges: Object.freeze(injections.map(({ path, consumer, provider, line }) => Object.freeze({ path, consumer, provider, line }))),
    subjects: Object.freeze(unique(injections.map(({ subject }) => subject)))
  });
};

const detectForwardClosuresV1 = (sources, registry) => {
  const { forwardClosures } = detectComposition(sources, registry);
  return Object.freeze({
    observedClosures: Object.freeze(forwardClosures.map(({ path, consumer, port, provider, method, line }) => (
      Object.freeze({ path, consumer, port, provider, method, line })
    ))),
    subjects: Object.freeze(unique(forwardClosures.map(({ subject }) => subject)))
  });
};

// `state.<name> = <function | owner member>` outside the state module.
const detectStateBackdoorsV1 = (sources, registryInput) => {
  const registry = resolveRegistry(registryInput);
  const observed = [];
  for (const [path, source] of sourceEntries(sources)) {
    if (path === registry.stateModule) continue;
    const code = maskJavaScript(source);
    const pattern = /(?:^|[;{}\n])[ \t]*state\.(\w+)\s*=(?!=)\s*([^\n;]*)/g;
    for (const match of code.matchAll(pattern)) {
      const name = match[1];
      const rhs = match[2].trim();
      const isFunction = /^(?:async\s*)?(?:\([^)]*\)|\w+)\s*=>/.test(rhs) || /^(?:async\s+)?function\b/.test(rhs);
      const isOwnerMember = /^\w+(?:\.\w+)+$/.test(rhs) && !/^(?:state|window|document|globalThis)\./.test(rhs);
      if (!isFunction && !isOwnerMember) continue;
      const index = match.index + match[0].indexOf('state.');
      observed.push({ path, name, line: lineNumberAt(source, index), kind: isFunction ? 'function' : 'owner-member', subject: `${path}|${name}` });
    }
  }
  const sorted = observed.sort((left, right) => left.subject.localeCompare(right.subject) || left.line - right.line);
  return Object.freeze({
    observedAssignments: Object.freeze(sorted.map(({ path, name, line, kind }) => Object.freeze({ path, name, line, kind }))),
    subjects: Object.freeze(unique(sorted.map(({ subject }) => subject)))
  });
};

// Exported factories outside the composition roots whose destructured
// parameters receive whole owner objects.
const detectWholeObjectPortsV1 = (sources, registryInput) => {
  const registry = resolveRegistry(registryInput);
  const roots = new Set(registry.compositionRoots);
  const isOwner = ownerNameTester(registry, compositionFactoryNames(sources, registry));
  const observed = [];
  for (const [path, source] of sourceEntries(sources)) {
    if (roots.has(path)) continue;
    const code = maskJavaScript(source);
    const pattern = /export\s+const\s+(create\w+)\s*=\s*(?:async\s*)?\(\s*\{/g;
    for (const match of code.matchAll(pattern)) {
      const openIndex = match.index + match[0].length - 1;
      const closeIndex = matchingClose(code, openIndex);
      if (closeIndex < 0) continue;
      const entries = objectLiteralEntries(code, source, openIndex);
      entries.forEach((entry) => {
        const name = entry.code.trim().match(/^(\w+)/)?.[1];
        if (!name || !isOwner(name)) return;
        observed.push({ path, factory: match[1], parameter: name, line: lineNumberAt(source, match.index), subject: `${path}|${match[1]}|${name}` });
      });
    }
  }
  const sorted = observed.sort((left, right) => left.subject.localeCompare(right.subject));
  return Object.freeze({
    observedParameters: Object.freeze(sorted.map(({ path, factory, parameter, line }) => Object.freeze({ path, factory, parameter, line }))),
    subjects: Object.freeze(unique(sorted.map(({ subject }) => subject)))
  });
};

// Distinct call shapes of each projection domain's functions outside its
// owners. A shape is the normalized source line of the call.
const normalizeShape = (line) => line.trim().replace(/\s+/g, ' ');
const isDeclarationOrMapping = (line, fnRaw) => {
  const fn = escapeRegExp(fnRaw);
  return new RegExp(`^(?:export\\s+)?(?:const|let|function)\\s+${fn}\\b`).test(line)
    || new RegExp(`^${fn}\\s*,?$`).test(line)
    || new RegExp(`^${fn}\\s*:\\s*[\\w.]+\\s*,?$`).test(line)
    || new RegExp(`^\\w+\\s*:\\s*(?:\\w+\\.)?${fn}\\s*,?$`).test(line)
    || /^(?:import|export)\b/.test(line);
};
const detectProjectionCallShapesV1 = (sources, registryInput) => {
  const registry = resolveRegistry(registryInput);
  const entries = sourceEntries(sources);
  const observed = [];
  registry.projectionDomains.forEach((domain) => {
    const owners = new Set(domain.owners);
    const pattern = new RegExp(`\\b(?:${identifierPattern(domain.functions)})\\s*\\(`, 'g');
    entries.forEach(([path, source]) => {
      if (owners.has(path)) return;
      const code = maskJavaScript(source);
      const seenLines = new Set();
      for (const match of code.matchAll(pattern)) {
        const lineNumber = lineNumberAt(source, match.index);
        if (seenLines.has(lineNumber)) continue;
        seenLines.add(lineNumber);
        const line = normalizeShape(lineOf(source, match.index));
        const fn = match[0].replace(/\s*\($/, '');
        if (isDeclarationOrMapping(line, fn)) continue;
        observed.push({ domain: domain.name, path, line: lineNumber, shape: line, subject: `${domain.name}|${path}|${line}` });
      }
    });
  });
  const sorted = observed.sort((left, right) => left.subject.localeCompare(right.subject));
  const shapesByDomain = {};
  registry.projectionDomains.forEach(({ name }) => { shapesByDomain[name] = unique(sorted.filter((record) => record.domain === name).map(({ subject }) => subject)).length; });
  return Object.freeze({
    observedShapes: Object.freeze(sorted.map(({ domain, path, line, shape }) => Object.freeze({ domain, path, line, shape }))),
    shapesByDomain: Object.freeze(shapesByDomain),
    subjects: Object.freeze(unique(sorted.map(({ subject }) => subject)))
  });
};

// Trigger sites of a heavy derived value's producer outside its owner.
const detectHeavyDerivedTriggerSitesV1 = (sources, registryInput) => {
  const registry = resolveRegistry(registryInput);
  const entries = sourceEntries(sources);
  const observed = [];
  registry.heavyProducers.forEach((producer) => {
    const pattern = new RegExp(`\\b${escapeRegExp(producer.name)}\\s*\\??\\.\\s*(${identifierPattern(producer.methods)})\\s*\\??\\.?\\s*\\(`, 'g');
    entries.forEach(([path, source]) => {
      if (path === producer.owner) return;
      const code = maskJavaScript(source);
      for (const match of code.matchAll(pattern)) {
        observed.push({ producer: producer.name, path, method: match[1], line: lineNumberAt(source, match.index), subject: `${producer.name}|${path}` });
      }
    });
  });
  const sorted = observed.sort((left, right) => left.subject.localeCompare(right.subject) || left.line - right.line);
  const countsBySubject = {};
  sorted.forEach(({ subject }) => { countsBySubject[subject] = (countsBySubject[subject] || 0) + 1; });
  return Object.freeze({
    observedSites: Object.freeze(sorted.map(({ producer, path, method, line }) => Object.freeze({ producer, path, method, line }))),
    countsBySubject: Object.freeze(countsBySubject),
    siteCount: sorted.length,
    subjects: Object.freeze(unique(sorted.map(({ subject }) => subject)))
  });
};

const encodeNamedSubject = (keys) => (record) => keys.map((key) => record[key]).join('|');

export const WEB_OWNER_GRAPH_DETECTORS = Object.freeze({
  'owner-graph.injection-edge.v1': Object.freeze({
    subjectCategory: 'injection-edge',
    encodeSubject: (record) => `${normalizeModulePath(record.path)}|${record.consumer}<-${record.provider}`,
    detect: detectInjectionEdgesV1
  }),
  'owner-graph.forward-closure.v1': Object.freeze({
    subjectCategory: 'forward-closure',
    encodeSubject: (record) => `${normalizeModulePath(record.path)}|${record.consumer}.${record.port}->${record.provider}.${record.method}`,
    detect: detectForwardClosuresV1
  }),
  'owner-graph.state-backdoor.v1': Object.freeze({
    subjectCategory: 'state-backdoor',
    encodeSubject: (record) => `${normalizeModulePath(record.path)}|${record.name}`,
    detect: detectStateBackdoorsV1
  }),
  'owner-graph.whole-object-port.v1': Object.freeze({
    subjectCategory: 'whole-object-port',
    encodeSubject: (record) => `${normalizeModulePath(record.path)}|${record.factory}|${record.parameter}`,
    detect: detectWholeObjectPortsV1
  }),
  'projection.call-shape.v1': Object.freeze({
    subjectCategory: 'projection-call-shape',
    encodeSubject: (record) => `${record.domain}|${normalizeModulePath(record.path)}|${normalizeShape(record.shape)}`,
    detect: detectProjectionCallShapesV1
  }),
  'heavy-derived.trigger-site.v1': Object.freeze({
    subjectCategory: 'trigger-site',
    encodeSubject: encodeNamedSubject(['producer', 'path']),
    detect: detectHeavyDerivedTriggerSitesV1
  })
});

export const WEB_OWNER_GRAPH_DETECTOR_IDS = Object.freeze(Object.keys(WEB_OWNER_GRAPH_DETECTORS));

// One pass over a source map; returns every detector's result by id.
export const detectWebOwnerGraph = (sources, registry = null) => Object.freeze(
  Object.fromEntries(WEB_OWNER_GRAPH_DETECTOR_IDS.map((id) => [id, WEB_OWNER_GRAPH_DETECTORS[id].detect(sources, registry)]))
);

// Summary counts for reports and trend tables.
export const summarizeWebOwnerGraph = (results) => Object.freeze({
  injectionEdges: results['owner-graph.injection-edge.v1'].subjects.length,
  forwardClosures: results['owner-graph.forward-closure.v1'].subjects.length,
  stateBackdoors: results['owner-graph.state-backdoor.v1'].subjects.length,
  wholeObjectPorts: results['owner-graph.whole-object-port.v1'].subjects.length,
  projectionShapes: Object.freeze({ ...results['projection.call-shape.v1'].shapesByDomain }),
  triggerSites: results['heavy-derived.trigger-site.v1'].siteCount,
  triggerModules: results['heavy-derived.trigger-site.v1'].subjects.length
});
