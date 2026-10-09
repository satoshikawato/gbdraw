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
import { posix } from 'node:path';

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
  // R13 layers, lowest first (see `webLayerOf`). A module may import its own
  // layer and lower ones; `compositionRoots` and `stateModule` above are the
  // roots layer and the state layer.
  layers: Object.freeze({
    leaves: Object.freeze([
      'config.js', 'web-ux-profile.js', 'mode-profiles.generated.js', 'mode-profiles.js',
      'mode-scoped-settings.generated.js'
    ]),
    leafDirectories: Object.freeze(['utils/']),
    stateBoundServices: Object.freeze(['services/config.js', 'services/reset.js']),
    entryModules: Object.freeze(['app.js', 'components.js'])
  }),
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
      owners: Object.freeze(['services/legend-svg.js', 'app/legend/entry-actions.js']),
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
      methods: Object.freeze(['prepare', 'prepareDrawn', 'prepareVisibility', 'prepareCandidate', 'evaluate', 'run']),
      // Producer methods the v2 trigger-site detector adds to `methods`
      // (`runDrawn` runs the drawn-feature preparation like `run`).
      v2Methods: Object.freeze(['runDrawn'])
    })
  ]),
  // `owner-graph.identity-from-display.v1`: modules on the identity paths of
  // Legend rows and drawn features, the display values another owner writes
  // (caption, paint, shown text, position), and the identities that qualify a
  // comparison of them.
  identityFromDisplay: Object.freeze({
    modules: Object.freeze([
      'app/legend/', 'app/legend.js', 'app/candidate-render.js', 'app/feature-editor/', 'app/feature-editor.js',
      'app/rule-matching.js', 'app/svg-styles.js', 'app/result-paint-record.js',
      'services/legend-svg.js', 'services/svg-result-ingestion.js', 'services/specific-color-rules.js',
      'services/result-paint-bases.js', 'services/feature-identity.js', 'services/feature-override-identity.js',
      'services/feature-placement.js', 'services/label-override-table.js', 'services/feature-dom.js'
    ]),
    // A member (`x.caption`) or a bare local (`caption`).
    captionFields: Object.freeze(['caption', 'labelText']),
    paintFields: Object.freeze(['color', 'fill', 'stroke', 'strokeColor', 'fillColor', 'swatchColor']),
    // A member only.
    positionFields: Object.freeze(['xPos', 'yPos']),
    paintAttributes: Object.freeze(['fill', 'stroke']),
    // What a Result shows: these attributes, members, and accessor results.
    shownAttributes: Object.freeze(['fill', 'stroke', 'data-legend-key', 'transform']),
    shownFields: Object.freeze(['textContent', 'innerText', 'shownKey', 'shownCaption']),
    displayAccessors: Object.freeze(['legendCaption']),
    identityFields: Object.freeze([
      'originalCaption', 'recordedKey', 'biologicalFeatureId', 'recordKey', 'featureId', 'renderedId', 'sourceKey',
      'identityKey', 'slotId', 'id'
    ]),
    identityAccessors: Object.freeze(['generatedCaption', 'resultBaseAttribute']),
    identityAttributePattern: 'data-(?!legend-key)[a-z-]*(?:-id|-key)'
  })
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

// --- v2: nested forward closures and producer ports -------------------------

const NON_METHOD_KEYWORDS = new Set([
  'if', 'for', 'while', 'switch', 'catch', 'function', 'with', 'return', 'else', 'do', 'try', 'finally',
  'await', 'typeof', 'new', 'void', 'delete', 'in', 'of', 'instanceof', 'throw', 'case', 'yield'
]);

// Bracket pairing and the innermost enclosing open bracket of every index.
const bracketMaps = (code) => {
  const closeOf = new Map();
  const openOf = new Map();
  const enclosing = new Int32Array(code.length).fill(-1);
  const stack = [];
  const pairs = { '(': ')', '[': ']', '{': '}' };
  for (let index = 0; index < code.length; index += 1) {
    const character = code[index];
    enclosing[index] = stack.length ? stack[stack.length - 1] : -1;
    if (pairs[character]) stack.push(index);
    else if (character === ')' || character === ']' || character === '}') {
      const open = stack.length && pairs[code[stack[stack.length - 1]]] === character ? stack.pop() : -1;
      if (open >= 0) {
        closeOf.set(open, index);
        openOf.set(index, open);
        enclosing[index] = stack.length ? stack[stack.length - 1] : -1;
      }
    }
  }
  return { closeOf, openOf, enclosing };
};
const skipWhitespace = (code, index) => {
  let cursor = index;
  while (cursor < code.length && /\s/.test(code[cursor])) cursor += 1;
  return cursor;
};

// Every function in a module: arrow functions (block and expression bodies),
// `function` expressions and declarations, and object or class method
// shorthand. `name` is the binding, property key, or method name when there is
// one; `callee` names the call that receives an anonymous callback.
const functionRegions = (code, maps) => {
  const { closeOf, openOf, enclosing } = maps;
  const regions = new Map();
  // `const f = `, `f = `, or `key: ` before a function; a member assignment
  // (`state.f = `) names nothing.
  const bindingBefore = (start) => {
    const match = /(?<![\w$.])(?:(?:const|let|var)\s+)?([\w$]+)\s*([:=])\s*$/.exec(code.slice(Math.max(0, start - 120), start));
    return match && !NON_METHOD_KEYWORDS.has(match[1]) ? { name: match[1], binding: match[2] === '=' } : { name: null, binding: false };
  };
  const calleeOf = (start) => {
    const open = enclosing[start];
    if (open < 0 || code[open] !== '(') return null;
    return /([\w$.?]+)\s*$/.exec(code.slice(Math.max(0, open - 80), open))?.[1] || null;
  };
  const add = (region) => {
    if (region.bodyEnd > region.bodyStart && !regions.has(region.bodyStart)) regions.set(region.bodyStart, region);
  };
  for (const match of code.matchAll(/=>/g)) {
    let cursor = match.index - 1;
    while (cursor >= 0 && /\s/.test(code[cursor])) cursor -= 1;
    let start;
    let params;
    if (code[cursor] === ')') {
      const open = openOf.get(cursor);
      if (open === undefined) continue;
      start = open;
      params = code.slice(open + 1, cursor);
    } else {
      const identifier = /([\w$]+)$/.exec(code.slice(Math.max(0, cursor - 80), cursor + 1));
      if (!identifier) continue;
      start = cursor + 1 - identifier[1].length;
      params = identifier[1];
    }
    const asyncPrefix = /\basync\s*$/.exec(code.slice(Math.max(0, start - 12), start));
    if (asyncPrefix) start -= asyncPrefix[0].length;
    const bodyOpen = skipWhitespace(code, match.index + 2);
    let bodyStart;
    let bodyEnd;
    if (code[bodyOpen] === '{') {
      bodyStart = bodyOpen;
      bodyEnd = closeOf.get(bodyOpen) ?? -1;
    } else {
      bodyStart = bodyOpen;
      bodyEnd = bodyOpen;
      let depth = 0;
      for (; bodyEnd < code.length; bodyEnd += 1) {
        const character = code[bodyEnd];
        if (character === '(' || character === '[' || character === '{') depth += 1;
        else if (character === ')' || character === ']' || character === '}') {
          if (depth === 0) break;
          depth -= 1;
        } else if ((character === ',' || character === ';') && depth === 0) break;
      }
    }
    if (bodyEnd < 0) continue;
    add({ kind: 'arrow', start, params, bodyStart, bodyEnd, ...bindingBefore(start), callee: calleeOf(start) });
  }
  for (const match of code.matchAll(/\bfunction\b\s*\*?\s*([\w$]*)\s*\(/g)) {
    const paramsOpen = match.index + match[0].length - 1;
    const paramsClose = closeOf.get(paramsOpen);
    if (paramsClose === undefined) continue;
    const bodyOpen = skipWhitespace(code, paramsClose + 1);
    if (code[bodyOpen] !== '{' || !closeOf.has(bodyOpen)) continue;
    const asyncPrefix = /\basync\s*$/.exec(code.slice(Math.max(0, match.index - 12), match.index));
    const start = match.index - (asyncPrefix ? asyncPrefix[0].length : 0);
    const declaration = Boolean(match[1]) && /(?:^|[;{}\n])\s*$/.test(code.slice(Math.max(0, start - 40), start));
    add({
      kind: 'function', hoisted: declaration, start, params: code.slice(paramsOpen + 1, paramsClose), bodyStart: bodyOpen, bodyEnd: closeOf.get(bodyOpen),
      ...(match[1] ? { name: match[1], binding: true } : bindingBefore(start)), callee: calleeOf(start)
    });
  }
  for (const match of code.matchAll(/(?:^|[,{};\n])\s*(?:(?:async|get|set|static)\s+)?\*?\s*([\w$]+)\s*\(/g)) {
    if (NON_METHOD_KEYWORDS.has(match[1])) continue;
    const paramsOpen = match.index + match[0].length - 1;
    const paramsClose = closeOf.get(paramsOpen);
    if (paramsClose === undefined) continue;
    const bodyOpen = skipWhitespace(code, paramsClose + 1);
    if (code[bodyOpen] !== '{' || !closeOf.has(bodyOpen)) continue;
    add({
      kind: 'method', start: match.index + match[0].indexOf(match[1]), params: code.slice(paramsOpen + 1, paramsClose),
      bodyStart: bodyOpen, bodyEnd: closeOf.get(bodyOpen), name: match[1], binding: false, callee: null
    });
  }
  return [...regions.values()].sort((left, right) => left.bodyStart - right.bodyStart || right.bodyEnd - left.bodyEnd);
};

// The functions that contain `index`, outermost first.
const regionChain = (regions, index) => regions
  .filter((region) => region.bodyStart <= index && index <= region.bodyEnd)
  .sort((left, right) => (right.bodyEnd - right.bodyStart) - (left.bodyEnd - left.bodyStart));

// A function that declares `name` itself (parameter or local binding) does not
// reach the composition scope's binding of that name.
const declaresName = (code, region, name) => {
  const escaped = escapeRegExp(name);
  if (new RegExp(`(?<![\\w$.])${escaped}(?![\\w$])`).test(region.params)) return true;
  const body = code.slice(region.bodyStart, region.bodyEnd + 1);
  return new RegExp(`\\b(?:const|let|var)\\s+${escaped}(?![\\w$])`).test(body)
    || new RegExp(`\\b(?:const|let|var)\\s*[{[][^}\\]]*(?<![\\w$.])${escaped}(?![\\w$])`).test(body);
};

// The label of the closure a reference sits in: the first named function from
// the composition scope inward, else the call that receives the outermost
// anonymous callback.
const closureLabel = (chain) => {
  const named = chain.find((region) => region.name);
  if (named) return named.name;
  const callee = chain[0].callee;
  return callee ? `${callee}(callback)` : '(anonymous)';
};

// Closures v1 already reports: entries of a factory call's object literal that
// are closures. A reference inside one keeps its `consumer.port` subject.
const factoryClosureSpans = (code, source, calls) => {
  const spans = [];
  calls.forEach((call) => {
    const argumentsCode = code.slice(call.openIndex + 1, call.closeIndex);
    const objectOpen = argumentsCode.search(/\{/);
    if (objectOpen < 0) return;
    objectLiteralEntries(argumentsCode, source.slice(call.openIndex + 1, call.closeIndex), objectOpen).forEach((entry) => {
      const trimmed = entry.code.trim();
      const closure = trimmed.match(/^(\w+)\s*:\s*(?:async\s*)?(?:\([^)]*\)|\w+)\s*=>/)
        || trimmed.match(/^(?:async\s+)?(\w+)\s*\([^)]*\)\s*\{/);
      if (closure) spans.push({ start: call.openIndex + 1 + entry.start, end: call.openIndex + 1 + entry.end });
    });
  });
  return spans;
};

// Where a closure of a composition scope first becomes reachable by another
// owner. A function written as a closure (an argument, an object property, a
// method) is reachable where it is written. A function bound to a name
// (`const f = () => ...`) is reachable where the scope's own statements name
// it, or where another reachable function names it: a helper only the
// returned bindings call is not reachable before the owners exist. A function
// bound with `const` cannot run before its own statement ends, so a chain that
// passes through one defined after an owner is created never reaches it earlier.
const scopeReachability = (code, regions, scope) => {
  const children = regions.filter((region) => {
    const chain = regionChain(regions, region.bodyStart);
    return chain.at(-2) === scope && chain.at(-1) === region;
  });
  const childOf = (index) => {
    const chain = regionChain(regions, index);
    const position = chain.indexOf(scope);
    return position < 0 ? undefined : (chain[position + 1] || null);
  };
  const bound = (child) => Boolean(child.binding && child.name);
  const reachable = new Map(children.map((child) => [child, bound(child) ? Infinity : child.start]));
  const runsFrom = (child) => (bound(child) && !child.hoisted ? Math.max(reachable.get(child), child.bodyEnd) : reachable.get(child));
  const edges = [];
  children.filter(bound).forEach((child) => {
    const pattern = new RegExp(`(?<![\\w$.])${escapeRegExp(child.name)}(?![\\w$])`, 'g');
    for (const match of code.matchAll(pattern)) {
      const index = match.index;
      if (index >= child.bodyStart && index <= child.bodyEnd) continue;
      if (index < scope.bodyStart || index > scope.bodyEnd) continue;
      const before = code.slice(Math.max(0, index - 12), index);
      const after = code.slice(index + child.name.length, index + child.name.length + 8);
      if (/\b(?:const|let|var|function)\s+$/.test(before) || /^\s*=(?![=>])/.test(after)) continue;
      // A property key (`name: value`) names nothing of this scope.
      if (/^\s*:/.test(after) && /[{,]\s*$/.test(before)) continue;
      const owner = childOf(index);
      if (owner === undefined || owner === child) continue;
      if (owner === null) reachable.set(child, Math.min(reachable.get(child), index));
      else edges.push([child, owner]);
    }
  });
  for (let changed = true; changed;) {
    changed = false;
    edges.forEach(([child, owner]) => {
      if (runsFrom(owner) < reachable.get(child)) {
        reachable.set(child, runsFrom(owner));
        changed = true;
      }
    });
  }
  return new Map(children.map((child) => [child, runsFrom(child)]));
};

// Forward closures in the composition roots at any depth: a function nested
// anywhere in the root's scope that reads an owner created later in that
// scope (`provider.method`), or a late-bound function variable (`let x = () =>
// {}` reassigned to its real function later), while the function is already
// reachable by another owner (`scopeReachability`). v1's factory-argument
// closures keep their v1 subjects.
const detectNestedForwardClosures = (sources, registryInput) => {
  const registry = resolveRegistry(registryInput);
  const roots = new Set(registry.compositionRoots);
  const records = [];
  for (const [path, source] of sourceEntries(sources)) {
    if (!roots.has(path)) continue;
    const code = maskJavaScript(source);
    const calls = factoryCalls(code, source);
    if (!calls.length) continue;
    const maps = bracketMaps(code);
    const regions = functionRegions(code, maps);
    const v1Spans = factoryClosureSpans(code, source, calls);
    const insideV1Entry = (index) => v1Spans.some(({ start, end }) => start <= index && index < end);
    const innermostOf = (index) => regionChain(regions, index).at(-1) || null;
    const reachability = new Map();
    const reachableAt = (scope, child) => {
      if (!reachability.has(scope)) reachability.set(scope, scopeReachability(code, regions, scope));
      return reachability.get(scope).get(child);
    };
    const reportReferences = ({ name, scope, limit, pattern, method, skipV1 }) => {
      for (const match of code.matchAll(pattern)) {
        const index = match.index + match[0].indexOf(name);
        if (index >= limit || (skipV1 && insideV1Entry(index))) continue;
        const chain = regionChain(regions, index);
        const scopeIndex = chain.indexOf(scope);
        if (scopeIndex < 0 || scopeIndex === chain.length - 1) continue;
        const closureChain = chain.slice(scopeIndex + 1);
        if (closureChain.some((region) => declaresName(code, region, name))) continue;
        if (!(reachableAt(scope, closureChain[0]) < limit)) continue;
        const label = closureLabel(closureChain);
        const member = method ? method(match) : '';
        records.push({
          path, consumer: label, provider: name, method: member, line: lineNumberAt(source, index),
          subject: `${path}|${label}->${name}${member ? `.${member}` : ''}`
        });
      }
    };
    // Owners: a factory result referenced by `owner.method` before its creation.
    const created = new Map();
    calls.forEach((call) => { if (!created.has(call.consumer)) created.set(call.consumer, call); });
    created.forEach((call, name) => {
      const scope = innermostOf(call.index);
      if (!scope) return;
      reportReferences({
        name, scope, limit: call.index,
        pattern: new RegExp(`(?<![\\w$.])${escapeRegExp(name)}\\s*\\??\\.\\s*([\\w$]+)`, 'g'),
        method: (match) => match[1],
        skipV1: true
      });
    });
    // Late-bound function variables: `let name = ...` reassigned to a function.
    for (const declaration of code.matchAll(/\blet\s+([\w$]+)\b/g)) {
      const name = declaration[1];
      if (created.has(name)) continue;
      const scope = innermostOf(declaration.index);
      if (!scope) continue;
      const assignment = new RegExp(`(?:^|[;{}\\n])\\s*${escapeRegExp(name)}\\s*=(?![=>])\\s*`, 'g');
      let firstFunctionAssignment = -1;
      for (const match of code.matchAll(assignment)) {
        if (match.index < declaration.index) continue;
        const valueStart = match.index + match[0].length;
        const nameIndex = match.index + match[0].indexOf(name);
        if (innermostOf(nameIndex) !== scope) continue;
        if (regions.some((region) => region.start === valueStart)) {
          firstFunctionAssignment = nameIndex;
          break;
        }
      }
      if (firstFunctionAssignment < 0) continue;
      reportReferences({
        name, scope, limit: firstFunctionAssignment,
        pattern: new RegExp(`(?<![\\w$.])${escapeRegExp(name)}(?![\\w$])(?!\\s*:)`, 'g'),
        method: null
      });
    }
  }
  return records;
};

const detectForwardClosuresV2 = (sources, registry) => {
  const v1 = detectComposition(sources, registry).forwardClosures;
  const nested = detectNestedForwardClosures(sources, registry);
  const closures = new Map();
  v1.forEach((record) => closures.set(record.subject, { path: record.path, consumer: `${record.consumer}.${record.port}`, provider: record.provider, method: record.method, line: record.line, subject: record.subject }));
  nested.forEach((record) => { if (!closures.has(record.subject)) closures.set(record.subject, record); });
  const sorted = [...closures.values()].sort((left, right) => left.subject.localeCompare(right.subject));
  return Object.freeze({
    observedClosures: Object.freeze(sorted.map(({ path, consumer, provider, method, line }) => Object.freeze({ path, consumer, provider, method, line }))),
    // Every nested reference, not only the first of each subject.
    observedReferences: Object.freeze(nested.map(({ path, consumer, provider, method, line }) => Object.freeze({ path, consumer, provider, method, line }))
      .sort((left, right) => left.path.localeCompare(right.path) || left.line - right.line)),
    subjects: Object.freeze(unique(sorted.map(({ subject }) => subject)))
  });
};

// Trigger sites of a heavy derived value's producer, including producer ports:
// `key: producer.method` handed to a factory is a port, and a call of that
// port in the factory's module triggers the producer like a direct call.
const factoryModules = (entries) => {
  const modules = new Map();
  entries.forEach(([path, source]) => {
    for (const match of maskJavaScript(source).matchAll(/export\s+const\s+(create\w+)\s*=\s*(?:async\s*)?\(\s*\{/g)) {
      if (!modules.has(match[1])) modules.set(match[1], { path, source, openIndex: match.index + match[0].length - 1 });
    }
  });
  return modules;
};
const producerPortBindings = (entries, producer, methods) => {
  const bindings = [];
  const pattern = new RegExp(`(?<![\\w$.])${escapeRegExp(producer.name)}\\s*\\??\\.\\s*(${identifierPattern(methods)})(?![\\w$])(?!\\s*\\??\\.?\\s*\\()`, 'g');
  entries.forEach(([path, source]) => {
    if (path === producer.owner) return;
    const code = maskJavaScript(source);
    const maps = bracketMaps(code);
    for (const match of code.matchAll(pattern)) {
      const key = /([\w$]+)\s*:\s*$/.exec(code.slice(Math.max(0, match.index - 80), match.index));
      if (!key) continue;
      const objectOpen = maps.enclosing[match.index];
      if (objectOpen < 0 || code[objectOpen] !== '{') continue;
      const callOpen = maps.enclosing[objectOpen];
      if (callOpen < 0 || code[callOpen] !== '(') continue;
      const factory = /([\w$]+)\s*$/.exec(code.slice(Math.max(0, callOpen - 80), callOpen))?.[1];
      if (!factory || !/^create/.test(factory)) continue;
      bindings.push({ producer: producer.name, method: match[1], boundIn: path, factory, port: key[1], line: lineNumberAt(source, match.index) });
    }
  });
  return bindings;
};
const detectHeavyDerivedTriggerSitesV2 = (sources, registryInput) => {
  const registry = resolveRegistry(registryInput);
  const entries = sourceEntries(sources);
  const modules = factoryModules(entries);
  const observed = [];
  registry.heavyProducers.forEach((producer) => {
    const methods = unique([...producer.methods, ...(producer.v2Methods || [])]);
    const direct = new RegExp(`(?<![\\w$.])${escapeRegExp(producer.name)}\\s*\\??\\.\\s*(${identifierPattern(methods)})(?![\\w$])\\s*\\??\\.?\\s*\\(`, 'g');
    entries.forEach(([path, source]) => {
      if (path === producer.owner) return;
      for (const match of maskJavaScript(source).matchAll(direct)) {
        observed.push({ producer: producer.name, path, method: match[1], via: null, line: lineNumberAt(source, match.index), subject: `${producer.name}|${path}` });
      }
    });
    producerPortBindings(entries, producer, methods).forEach((binding) => {
      const factory = modules.get(binding.factory);
      if (!factory || factory.path === producer.owner) return;
      const code = maskJavaScript(factory.source);
      const parameter = objectLiteralEntries(code, factory.source, factory.openIndex)
        .map((entry) => /^([\w$]+)\s*(?::\s*([\w$]+))?/.exec(entry.code.trim()))
        .find((parsed) => parsed && parsed[1] === binding.port);
      if (!parameter) return;
      const local = parameter[2] || parameter[1];
      const call = new RegExp(`(?<![\\w$.])${escapeRegExp(local)}\\s*\\??\\.?\\s*\\(`, 'g');
      for (const match of code.matchAll(call)) {
        // The parameter list itself (`name = (commit) => ...`) is not a call.
        if (match.index > factory.openIndex && match.index < matchingClose(code, factory.openIndex)) continue;
        observed.push({
          producer: producer.name, path: factory.path, method: binding.method, via: `${binding.boundIn}:${binding.factory}.${binding.port}`,
          line: lineNumberAt(factory.source, match.index), subject: `${producer.name}|${factory.path}`
        });
      }
    });
  });
  const sorted = observed.sort((left, right) => left.subject.localeCompare(right.subject) || left.line - right.line);
  const countsBySubject = {};
  sorted.forEach(({ subject }) => { countsBySubject[subject] = (countsBySubject[subject] || 0) + 1; });
  return Object.freeze({
    observedSites: Object.freeze(sorted.map(({ producer, path, method, via, line }) => Object.freeze({ producer, path, method, via, line }))),
    countsBySubject: Object.freeze(countsBySubject),
    siteCount: sorted.length,
    subjects: Object.freeze(unique(sorted.map(({ subject }) => subject)))
  });
};

// --- import direction (R13 layers) --------------------------------------------

// Layer of a module under gbdraw/web/js, lowest first: 0 leaves (`utils/` and
// the constants modules), 1 state-free `services/`, 2 `state.js`, 3 the
// state-bound services, 4 owner modules (`app/`), 5 composition roots, 6
// `app.js` and `components.js`. Anything else (`workers/`, a non-module file)
// is unranked: `null`.
export const WEB_LAYER_NAMES = Object.freeze([
  'leaf', 'state-free service', 'state', 'state-bound service', 'owner module', 'composition root', 'entry module'
]);
export const webLayerOf = (path, registry = WEB_OWNER_GRAPH_DEFAULTS) => {
  const module = normalizeModulePath(path);
  const { layers } = registry;
  if (layers.leaves.includes(module) || layers.leafDirectories.some((directory) => module.startsWith(directory))) return 0;
  if (module === registry.stateModule) return 2;
  if (layers.stateBoundServices.includes(module)) return 3;
  if (module.startsWith('services/')) return 1;
  if (registry.compositionRoots.includes(module)) return 5;
  if (module.startsWith('app/')) return 4;
  if (layers.entryModules.includes(module)) return 6;
  return null;
};

// Static, `export ... from`, side-effect, and literal dynamic imports of one
// module, with comments masked (a JSDoc `@import` is not an edge) and the
// names each statement takes: the imported names, `default`, or `*` for a
// namespace or a side-effect import.
const importStatements = (source) => {
  const commentsMasked = maskJavaScript(source, { strings: false });
  const code = maskJavaScript(source);
  const statements = [];
  for (const match of commentsMasked.matchAll(/(?:^|\n)[ \t]*(import|export)\s+((?:[^;'"]*?)\s+from\s+)?(['"])([^'"]+)\3/g)) {
    const keywordIndex = match.index + match[0].search(/\b(?:import|export)\b/);
    if (!/^(?:import|export)\b/.test(code.slice(keywordIndex))) continue;
    const clause = (match[2] || '').replace(/\s+from\s+$/, '').trim();
    const names = [];
    const braces = /\{([^}]*)\}/s.exec(clause);
    if (braces) braces[1].split(',').map((part) => part.trim()).filter(Boolean).forEach((part) => names.push(part.split(/\s+as\s+/)[0]));
    if (/\*/.test(clause)) names.push('*');
    if (/^[\w$]+/.test(clause.replace(/\{[^}]*\}/s, '').replace(/\*\s*(?:as\s+[\w$]+)?/, '').replace(/^[\s,]+/, ''))) names.push('default');
    statements.push({ specifier: match[4], index: keywordIndex, names: names.length ? names : ['*'] });
  }
  for (const match of commentsMasked.matchAll(/\bimport\s*\(\s*(['"])([^'"]+)\1\s*\)/g)) {
    statements.push({ specifier: match[2], index: match.index, names: ['*'] });
  }
  return statements;
};

// An import that names a module in a higher layer than the importer. Subject
// `importer->target`; its count is the distinct imported names, so a baseline
// can shrink one name at a time.
const detectLayerImportDirectionV1 = (sources, registryInput) => {
  const registry = resolveRegistry(registryInput);
  const entries = sourceEntries(sources);
  const modules = new Set(entries.map(([path]) => path));
  const observed = [];
  entries.forEach(([path, source]) => {
    const layer = webLayerOf(path, registry);
    if (layer === null) return;
    importStatements(source).forEach(({ specifier, index, names }) => {
      if (!specifier.startsWith('.')) return;
      const resolved = posix.normalize(posix.join(posix.dirname(path), specifier));
      const target = [resolved, `${resolved}.js`].find((candidate) => modules.has(candidate));
      const targetLayer = target === undefined ? null : webLayerOf(target, registry);
      if (targetLayer === null || targetLayer <= layer) return;
      observed.push({ path, target, line: lineNumberAt(source, index), names, layer, targetLayer, subject: `${path}->${target}` });
    });
  });
  const sorted = observed.sort((left, right) => left.subject.localeCompare(right.subject) || left.line - right.line);
  const namesBySubject = new Map();
  sorted.forEach(({ subject, names }) => {
    namesBySubject.set(subject, new Set([...(namesBySubject.get(subject) || []), ...names]));
  });
  const countsBySubject = {};
  namesBySubject.forEach((names, subject) => { countsBySubject[subject] = names.size; });
  return Object.freeze({
    observedImports: Object.freeze(sorted.map(({ path, target, line, names, layer, targetLayer }) => Object.freeze({
      path, target, line, names: Object.freeze([...names]),
      layers: `${WEB_LAYER_NAMES[layer]} -> ${WEB_LAYER_NAMES[targetLayer]}`
    }))),
    countsBySubject: Object.freeze(countsBySubject),
    nameCount: Object.values(countsBySubject).reduce((total, count) => total + count, 0),
    subjects: Object.freeze(unique(sorted.map(({ subject }) => subject)))
  });
};

// --- identity from display -------------------------------------------------------

// A Legend row's or drawn feature's identity, or a membership decision, joined
// on a display value (caption, paint, shown text, position) instead of an
// identity. Sites are equality comparisons in a join position and lookup keys:
// - `search`: in a find/findIndex/findLast*/filter/some/every callback;
// - `predicate`: in any other expression-bodied callback;
// - `loop-join`: in a for/while body or a forEach/map/flatMap/reduce callback,
//   with display values on both sides (paint of the same field on both sides is
//   change detection, `paint-compare`);
// - `key-*`: a display value as a Map/Set key, computed member, `[key, value]`
//   pair of a Map, key function, or `[data-legend-key="${...}"]` selector.
// A literal or `typeof` operand excludes a comparison. One conjoined (&&) with
// another comparison against an identity is a tie-break (`id-qualified`), and
// `if (a !== b)` whose consequent writes `a` or `b` is change detection
// (`write-if-different`). Each other site gets a class: `shown-join` (an
// operand reads what a Result shows, or is a local bound from such a read in
// the function that holds the binding), `paint-join`, or else `caption-key`.
// `shown-join` and `paint-join` are subjects; the rest is a report-only
// inventory per module.
const IDENTITY_FROM_DISPLAY_SUBJECT_CLASSES = new Set(['paint-join', 'shown-join']);
const IDENTITY_SEARCH_CALLEE = /(?:^|\.)(?:find|findIndex|findLast|findLastIndex|filter|some|every)$/;
const IDENTITY_ITERATE_CALLEE = /(?:^|\.)(?:forEach|map|flatMap|reduce)$/;
const DISPLAY_LITERAL_OPERAND = /^\s*!*\s*(?:'[^']*'|"[^"]*"|`[^`$]*`|-?\d[\d._]*|null|undefined|true|false|NaN|\[\s*\]|\{\s*\})\s*$/;
const OPERAND_STOP_BEFORE = /(?:&&|\|\||\?\?|[,;:{}]|(?<![=!<>])=(?!=)|=>|\breturn|\bcase|\?)$/;
const isOptionalChain = (code, index) => code[index] === '?' && code[index + 1] === '.';
// The value an operand carries through `String(x || '')`, `.trim()`, and case
// conversion.
const VALUE_METHOD_SUFFIX = /(?:\s*\??\.\s*(?:trim|toLowerCase|toUpperCase|toString)\s*\(\s*\))+$/;
const STRING_WRAPPER = /^String\(\s*([\s\S]*?)\s*(?:\|\|\s*(?:''|""|``))?\s*\)$/;
const carriedValue = (operand) => {
  let value = operand.trim();
  for (let previous = ''; previous !== value;) {
    previous = value;
    value = value.replace(VALUE_METHOD_SUFFIX, '').replace(STRING_WRAPPER, '$1').trim();
  }
  return value;
};

// The start of the operand that ends at `index`, and the end of the operand
// that starts at `index`, at bracket depth 0.
const operandStartBefore = (code, index) => {
  let depth = 0;
  let cursor = index - 1;
  for (; cursor >= 0; cursor -= 1) {
    const character = code[cursor];
    if (character === ')' || character === ']' || character === '}') {
      depth += 1;
      continue;
    }
    if (character === '(' || character === '[' || character === '{') {
      if (depth === 0) break;
      depth -= 1;
      continue;
    }
    if (depth === 0 && !isOptionalChain(code, cursor) && OPERAND_STOP_BEFORE.test(code.slice(Math.max(0, cursor - 6), cursor + 1))) break;
  }
  return cursor + 1;
};
const operandEndAfter = (code, index) => {
  let depth = 0;
  let cursor = index;
  for (; cursor < code.length; cursor += 1) {
    const character = code[cursor];
    if (character === '(' || character === '[' || character === '{') {
      depth += 1;
      continue;
    }
    if (character === ')' || character === ']' || character === '}') {
      if (depth === 0) break;
      depth -= 1;
      continue;
    }
    if (depth === 0) {
      const two = code.slice(cursor, cursor + 2);
      if (two === '&&' || two === '||' || two === '??' || character === ',' || character === ';' || character === ':'
        || (character === '?' && code[cursor + 1] !== '.')) break;
      if (/^(?:===|!==|==|!=)/.test(code.slice(cursor))) break;
    }
  }
  return cursor;
};

// Operand classifiers built from the registry.
const identityFromDisplayPatterns = (config) => {
  const members = [...config.captionFields, ...config.paintFields, ...config.positionFields, ...config.shownFields];
  const names = [...config.captionFields, ...config.paintFields, ...config.shownFields];
  const display = new RegExp([
    `(?:\\.|\\?\\.)(?:${identifierPattern(members)})\\b`,
    `getAttribute\\(\\s*['"](?:${identifierPattern(config.shownAttributes)})['"]`,
    `\\b(?:${identifierPattern(config.displayAccessors)})\\s*\\(`,
    `(?<![\\w$.])(?:${identifierPattern(names)})(?![\\w$])(?!\\s*:)`
  ].join('|'));
  const identity = new RegExp([
    `(?:\\.|\\?\\.)(?:${identifierPattern(config.identityFields)})\\b`,
    `\\b(?:${identifierPattern(config.identityAccessors)})\\s*\\(`,
    `getAttribute\\(\\s*['"]${config.identityAttributePattern}['"]\\)`
  ].join('|'));
  const paint = new RegExp([
    `(?:\\.|\\?\\.)(?:${identifierPattern(config.paintFields)})\\b`,
    `getAttribute\\(\\s*['"](?:${identifierPattern(config.paintAttributes)})['"]`,
    `(?<![\\w$.])(?:${identifierPattern(config.paintFields)})(?![\\w$])(?!\\s*:)`
  ].join('|'));
  const paintRead = new RegExp(
    `(?:\\.|\\?\\.|(?<![\\w$.]))(${identifierPattern(config.paintFields)})\\b|getAttribute\\(\\s*['"](${identifierPattern(config.paintAttributes)})['"]`,
    'g'
  );
  const shown = new RegExp(
    `getAttribute\\(\\s*['"](?:${identifierPattern(config.shownAttributes)})['"]|\\.(?:${identifierPattern(config.shownFields)})\\b`
  );
  const attribute = /getAttribute\(\s*['"]([\w-]+)['"]/;
  const member = new RegExp(`(?:\\.|\\?\\.)(${identifierPattern(members)})\\b`);
  const accessor = new RegExp(`\\b(${identifierPattern(config.displayAccessors)})\\s*\\(`);
  const name = new RegExp(`(?<![\\w$.])(${identifierPattern(names)})(?![\\w$])`);
  return {
    shown,
    paint,
    accessorArgument: new RegExp(`[(,]\\s*(${identifierPattern(config.displayAccessors)})\\s*(?=[,)])`, 'g'),
    paintReads: (text) => [...text.matchAll(paintRead)].map((match) => match[1] || `@${match[2]}`),
    // An operand whose own value is a display value; an identity read that
    // wraps it (`ids.get(x.caption)`) does not make it one.
    isDisplay: (text) => display.test(text) && !(identity.test(text) && !display.test(text.replace(identity, ''))),
    isIdentity: (text) => identity.test(text),
    displayField: (text) => {
      const read = attribute.exec(text);
      if (read && config.shownAttributes.includes(read[1])) return `@${read[1]}`;
      const field = member.exec(text);
      if (field) return field[1];
      const call = accessor.exec(text);
      if (call) return `${call[1]}()`;
      return name.exec(text)?.[1] || null;
    }
  };
};

const identityFromDisplaySites = (path, source, patterns) => {
  const code = maskJavaScript(source);
  const text = maskJavaScript(source, { strings: false });
  const maps = bracketMaps(code);
  const regions = functionRegions(code, maps);
  const records = [];
  // The innermost named function that is not a callback, else the innermost
  // named one.
  const namedScope = (index) => [...regionChain(regions, index)].reverse().find((region) => region.name && !region.callee) || null;
  const subjectFunction = (index) => namedScope(index)?.name
    || regionChain(regions, index).find((region) => region.name)?.name || '(module)';
  // Locals bound to a shown read (or functions whose expression body is one),
  // from the binding to the end of the innermost function that holds it,
  // unless a nested function declares the same name.
  const tainted = [];
  const taintFrom = (name, index) => {
    const region = regionChain(regions, index).at(-1) || null;
    tainted.push({ name, region, start: index, end: region ? region.bodyEnd : code.length });
  };
  const taintAt = (index, name) => tainted.find((taint) => taint.name === name && taint.start <= index && index <= taint.end
    && !regionChain(regions, index).some((region) => region !== taint.region
      && (!taint.region || (region.bodyStart > taint.region.bodyStart && region.bodyEnd <= taint.region.bodyEnd))
      && declaresName(code, region, name)));
  const BLOCK_FUNCTION = /^(?:async\s*)?(?:\([^)]*\)|[\w$]+)\s*=>\s*\{|^(?:async\s+)?function\b/;
  for (const match of text.matchAll(/(?:const|let)\s+([\w$]+)\s*=\s*(?=([^;\n]*))/g)) {
    if (patterns.shown.test(match[2]) && !BLOCK_FUNCTION.test(match[2])) taintFrom(match[1], match.index);
  }
  for (const match of text.matchAll(/(?:^|[;{}\n])\s*([\w$]+)\s*=(?![=>])\s*(?=([^;\n]*))/g)) {
    const value = match[2].trim();
    if ((patterns.shown.test(value) && !BLOCK_FUNCTION.test(value)) || taintAt(match.index, value)) taintFrom(match[1], match.index);
  }
  const readsShown = (index, operand) => patterns.shown.test(operand) || tainted.some((taint) => (
    new RegExp(`(?<![\\w$.])${escapeRegExp(taint.name)}(?![\\w$])`).test(operand) && taintAt(index, taint.name) === taint
  ));
  // The local whose shown value an operand is (`key`, `keyOf(group)`), not one
  // it merely contains (`drawsFeatures(key)`).
  const shownLocal = (index, operand) => {
    const value = carriedValue(operand);
    const name = /^([\w$]+)(?:\s*\(([\s\S]*)\))?$/.exec(value);
    if (!name || (name[2] !== undefined && /[()]/.test(name[2]))) return null;
    return taintAt(index, name[1]) ? name[1] : null;
  };
  const isOperand = (index, operand) => patterns.isDisplay(operand) || Boolean(shownLocal(index, operand));
  const fieldOf = (index, operand) => patterns.displayField(operand) || shownLocal(index, operand);
  // `if (a !== b) <write b>`: change detection, not a join.
  const writesCompared = (index, operands) => {
    const open = maps.enclosing[index];
    if (open < 0 || code[open] !== '(' || !/\bif\s*$/.test(code.slice(Math.max(0, open - 8), open))) return false;
    const bodyStart = skipWhitespace(code, maps.closeOf.get(open) + 1);
    const bodyEnd = code[bodyStart] === '{' ? maps.closeOf.get(bodyStart) : code.indexOf(';', bodyStart);
    const consequent = text.slice(bodyStart, bodyEnd + 1);
    return operands.some((operand) => new RegExp(
      `(?:(?<![=!<>])=(?![=>])|setAttribute\\(\\s*['"][\\w-]+['"]\\s*,)\\s*(?:String\\(\\s*)?${escapeRegExp(carriedValue(operand))}\\s*\\)?\\s*(?:[;,})\\n]|$)`
    ).test(consequent));
  };
  const classify = (index, form, [left, right = '']) => {
    if (form.startsWith('key-')) {
      if (readsShown(index, left)) return 'shown-join';
      return patterns.paint.test(left) ? 'paint-join' : 'caption-key';
    }
    if (readsShown(index, left) || readsShown(index, right)) return 'shown-join';
    const leftPaint = patterns.paint.test(left);
    const rightPaint = patterns.paint.test(right);
    if (form === 'loop-join') {
      if (!leftPaint || !rightPaint) return 'caption-key';
      return patterns.paintReads(left).join() !== patterns.paintReads(right).join() ? 'paint-join' : 'paint-compare';
    }
    return leftPaint || rightPaint ? 'paint-join' : 'caption-key';
  };
  const add = (index, form, field, snippet, exemption, operands = [snippet]) => records.push({
    path,
    line: lineNumberAt(source, index),
    function: subjectFunction(index),
    class: exemption || classify(index, form, operands),
    form,
    field,
    snippet: snippet.replace(/\s+/g, ' ').trim().slice(0, 140)
  });
  const loopBodies = [];
  for (const match of code.matchAll(/\b(?:for|while)\s*\(/g)) {
    const close = maps.closeOf.get(match.index + match[0].length - 1);
    if (close === undefined) continue;
    const bodyStart = skipWhitespace(code, close + 1);
    loopBodies.push([bodyStart, code[bodyStart] === '{' ? maps.closeOf.get(bodyStart) : code.indexOf(';', bodyStart)]);
  }
  const inLoop = (index) => loopBodies.some(([start, end]) => start <= index && index <= end);

  // Equality in a join position with a display operand.
  for (const match of code.matchAll(/(?<![=!<>])(?:===|!==|==|!=)(?!=)/g)) {
    const operatorEnd = match.index + match[0].length;
    const left = text.slice(operandStartBefore(code, match.index), match.index).trim();
    const right = text.slice(operatorEnd, operandEndAfter(code, operatorEnd)).trim();
    if (DISPLAY_LITERAL_OPERAND.test(left) || DISPLAY_LITERAL_OPERAND.test(right) || /^typeof\b/.test(left)) continue;
    const leftDisplay = isOperand(match.index, left);
    const rightDisplay = isOperand(match.index, right);
    if (!leftDisplay && !rightDisplay) continue;
    const innermost = regionChain(regions, match.index).at(-1);
    const expression = Boolean(innermost) && code[innermost.bodyStart] !== '{';
    let form = null;
    if (innermost?.callee && IDENTITY_SEARCH_CALLEE.test(innermost.callee)) form = 'search';
    else if (innermost?.callee && !IDENTITY_ITERATE_CALLEE.test(innermost.callee) && expression) form = 'predicate';
    else if (leftDisplay && rightDisplay && ((innermost?.callee && IDENTITY_ITERATE_CALLEE.test(innermost.callee))
      || (innermost && inLoop(match.index) && inLoop(innermost.bodyStart))
      || (!innermost?.callee && inLoop(match.index) && (!innermost || !inLoop(innermost.bodyStart))))) form = 'loop-join';
    if (!form) continue;
    // The `&&` conjunction around the comparison, inside its function.
    const scopeStart = innermost ? innermost.bodyStart : 0;
    const scopeEnd = innermost ? innermost.bodyEnd : code.length;
    let start = match.index;
    let depth = 0;
    for (; start > scopeStart; start -= 1) {
      const character = code[start];
      if (character === ')' || character === ']' || character === '}') depth += 1;
      else if (character === '(' || character === '[' || character === '{') {
        if (depth === 0) break;
        depth -= 1;
      } else if (depth === 0 && (character === ';' || character === ':' || code.slice(start - 1, start + 1) === '||'
        || (character === '?' && !isOptionalChain(code, start)))) break;
    }
    let end = match.index;
    depth = 0;
    for (; end < scopeEnd; end += 1) {
      const character = code[end];
      if (character === '(' || character === '[' || character === '{') depth += 1;
      else if (character === ')' || character === ']' || character === '}') {
        if (depth === 0) break;
        depth -= 1;
      } else if (depth === 0 && (character === ';' || character === ',' || character === ':' || code.slice(end, end + 2) === '||'
        || (character === '?' && !isOptionalChain(code, end)))) break;
    }
    const conjunction = text.slice(start + 1, end);
    const qualified = /&&/.test(conjunction) && [...conjunction.matchAll(/(?<![=!<>])(?:===|==)(?!=)/g)].some((comparison) => {
      const at = start + 1 + comparison.index;
      const after = at + comparison[0].length;
      // A comparison does not qualify itself.
      return at !== match.index && (patterns.isIdentity(text.slice(operandStartBefore(code, at), at))
        || patterns.isIdentity(text.slice(after, operandEndAfter(code, after))));
    });
    const exemption = (match[0].startsWith('!') && writesCompared(match.index, [left, right]) && 'write-if-different')
      || (qualified && 'id-qualified') || null;
    add(match.index, form, fieldOf(match.index, leftDisplay ? left : right) || fieldOf(match.index, right),
      `${left} ${match[0]} ${right}`, exemption, [left, right]);
  }

  // A display value as a Map/Set key (DOM setters, classList, searchParams,
  // and style are not keyed stores).
  for (const match of code.matchAll(/\.(get|has|set|delete|add)\s*\(/g)) {
    const open = match.index + match[0].length - 1;
    if (maps.closeOf.get(open) === undefined) continue;
    const argument = text.slice(open + 1, operandEndAfter(code, open + 1)).trim();
    if (!isOperand(match.index, argument) || DISPLAY_LITERAL_OPERAND.test(argument)) continue;
    const receiver = /([\w$.?\])]+)\s*$/.exec(code.slice(Math.max(0, match.index - 80), match.index))?.[1] || '';
    if (/classList$|searchParams$|style$/.test(receiver)) continue;
    add(match.index, `key-${match[1]}`, fieldOf(match.index, argument), `${receiver}.${match[1]}(${argument})`, null, [argument]);
  }
  // A computed member.
  for (const match of code.matchAll(/(?<=[\w$\])])\s*\[/g)) {
    const open = match.index + match[0].length - 1;
    const close = maps.closeOf.get(open);
    if (close === undefined) continue;
    const inner = text.slice(open + 1, close).trim();
    if (!isOperand(open, inner) || DISPLAY_LITERAL_OPERAND.test(inner) || /^\d+$/.test(inner)) continue;
    const receiver = /([\w$.?]+)\s*$/.exec(code.slice(Math.max(0, match.index - 80), match.index + 1))?.[1] || '';
    if (/^(?:return|const|let|var|case|in|of|typeof|await|yield|else)$/.test(receiver)) continue;
    add(open, 'key-index', fieldOf(open, inner), `${receiver}[${inner}]`, null, [inner]);
  }
  // The key of a `[key, value]` pair a Map or Object.fromEntries is built from.
  for (const match of code.matchAll(/=>\s*\[/g)) {
    const open = match.index + match[0].length - 1;
    if (maps.closeOf.get(open) === undefined) continue;
    const first = text.slice(open + 1, operandEndAfter(code, open + 1)).trim();
    if (!isOperand(match.index, first) || DISPLAY_LITERAL_OPERAND.test(first)) continue;
    if (!/new\s+Map\s*\(|fromEntries\s*\(/.test(code.slice(Math.max(0, match.index - 200), match.index))) continue;
    add(match.index, 'key-entries', fieldOf(match.index, first), `[${first}, ...]`, null, [first]);
  }
  // A display accessor handed over as a key function.
  for (const match of code.matchAll(patterns.accessorArgument)) {
    add(match.index, 'key-function', `${match[1]}()`, match[0], null);
  }
  // A selector built from a display value.
  for (const match of text.matchAll(/\[data-legend-key(?:[~|^$*]?=)["']?\$\{([^}]*)\}/g)) {
    add(match.index, 'key-selector', '@data-legend-key', match[0], null, ['${caption}']);
  }
  return records;
};

const detectIdentityFromDisplayV1 = (sources, registryInput) => {
  const config = resolveRegistry(registryInput).identityFromDisplay;
  const inScope = (path) => config.modules.some((prefix) => (prefix.endsWith('/') ? path.startsWith(prefix) : path === prefix));
  const patterns = identityFromDisplayPatterns(config);
  const records = sourceEntries(sources)
    .filter(([path]) => inScope(path))
    .flatMap(([path, source]) => identityFromDisplaySites(path, source, patterns))
    .sort((left, right) => left.path.localeCompare(right.path) || left.line - right.line);
  const countsBySubject = {};
  const inventory = {};
  records.forEach((record) => {
    if (IDENTITY_FROM_DISPLAY_SUBJECT_CLASSES.has(record.class)) {
      const subject = `${record.path}|${record.function}|${record.class}`;
      countsBySubject[subject] = (countsBySubject[subject] || 0) + 1;
    } else {
      const key = `${record.class}|${record.path}`;
      inventory[key] = (inventory[key] || 0) + 1;
    }
  });
  const sortedCounts = (counts) => Object.freeze(Object.fromEntries(Object.entries(counts).sort(([left], [right]) => left.localeCompare(right))));
  const subjectSites = records.filter((record) => IDENTITY_FROM_DISPLAY_SUBJECT_CLASSES.has(record.class));
  return Object.freeze({
    observedSites: Object.freeze(records.map((record) => Object.freeze({ ...record }))),
    countsBySubject: sortedCounts(countsBySubject),
    siteCount: subjectSites.length,
    subjects: Object.freeze(unique(Object.keys(countsBySubject))),
    inventory: sortedCounts(inventory)
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
  'owner-graph.forward-closure.v2': Object.freeze({
    subjectCategory: 'forward-closure',
    encodeSubject: (record) => `${normalizeModulePath(record.path)}|${record.consumer}->${record.provider}${record.method ? `.${record.method}` : ''}`,
    detect: detectForwardClosuresV2
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
  }),
  'heavy-derived.trigger-site.v2': Object.freeze({
    subjectCategory: 'trigger-site',
    encodeSubject: encodeNamedSubject(['producer', 'path']),
    detect: detectHeavyDerivedTriggerSitesV2
  }),
  'layer.import-direction.v1': Object.freeze({
    subjectCategory: 'layer-import',
    encodeSubject: (record) => `${normalizeModulePath(record.path)}->${record.target}`,
    detect: detectLayerImportDirectionV1
  }),
  'owner-graph.identity-from-display.v1': Object.freeze({
    subjectCategory: 'identity-from-display',
    encodeSubject: (record) => `${normalizeModulePath(record.path)}|${record.function}|${record.class}`,
    detect: detectIdentityFromDisplayV1
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
  triggerModules: results['heavy-derived.trigger-site.v1'].subjects.length,
  layerImports: results['layer.import-direction.v1'].subjects.length
});
