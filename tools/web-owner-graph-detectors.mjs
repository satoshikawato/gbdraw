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
      methods: Object.freeze(['prepare', 'prepareDrawn', 'prepareVisibility', 'prepareCandidate', 'evaluate', 'run']),
      // Producer methods the v2 trigger-site detector adds to `methods`
      // (`runDrawn` runs the drawn-feature preparation like `run`).
      v2Methods: Object.freeze(['runDrawn'])
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
