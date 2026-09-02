import { posix } from 'node:path';

import { literalImportSpecifiers, maskJavaScript } from './web-change-source.mjs';

const WEB_SOURCE_PREFIX = 'gbdraw/web/js/';
const normalizeModulePath = (path) => String(path || '')
  .replaceAll('\\', '/')
  .replace(new RegExp(`^${WEB_SOURCE_PREFIX}`), '');
const sourceEntries = (sources) => (
  sources instanceof Map ? [...sources] : Object.entries(sources || {})
).map(([path, source]) => [normalizeModulePath(path), String(source || '')])
  .sort(([left], [right]) => left.localeCompare(right));
const unique = (values) => [...new Set(values)];
const matches = (source, pattern) => [...source.matchAll(pattern)];

const resolvedProductionImport = (owner, specifier) => {
  if (!specifier.startsWith('.')) return null;
  return posix.normalize(posix.join(posix.dirname(owner), specifier));
};

const PRIVILEGED_CAPABILITY_SPECS = Object.freeze([
  Object.freeze({
    name: 'Render request',
    importTargets: Object.freeze(['services/session-request.js']),
    operatorPattern: /\bbuildCanonicalRenderRequest\s*\(/g
  }),
  Object.freeze({
    name: 'Diagram Worker',
    importTargets: Object.freeze([
      'services/diagram-generation.js',
      'services/diagram-worker-protocol.js',
      'workers/diagram-generation-worker.js'
    ]),
    operatorPattern: /\b(?:new\s+Worker|runDiagramGeneration|loadPyodide)\s*\(/g
  }),
  Object.freeze({
    name: 'Python helper',
    importTargets: Object.freeze(['app/python-helpers.js']),
    operatorPattern: /\b(?:runPythonAsync|runPython|PYTHON_HELPERS|globals\.get)\b/g
  }),
  Object.freeze({
    name: 'Resource staging',
    importTargets: Object.freeze([
      'services/diagram-resource-staging.js',
      'services/resource-payload-owner.js'
    ]),
    operatorPattern: /\b(?:createDiagramResourceTransport|stageRenderResources|setResourcePayloadOwner)\b/g
  }),
  Object.freeze({
    name: 'SVG/Result admission',
    importTargets: Object.freeze([
      'services/svg-sanitization.js',
      'services/svg-result-ingestion.js'
    ]),
    operatorPattern: /\b(?:sanitizeSvgContent|ingestSvgResult|markCommitted)\b/g
  }),
  Object.freeze({
    name: 'Mounted SVG/Result replacement',
    importTargets: Object.freeze([
      'services/svg-serialization.js',
      'app/preview-runtime.js'
    ]),
    operatorPattern: /\b(?:serializeCleanSvg|flushActiveResult)\s*\(|\b(?:results|state\.results)\.value(?:\[[^\]]+\])?\s*=(?!=)/g
  }),
  Object.freeze({
    name: 'History',
    importTargets: Object.freeze([
      'services/history.js',
      'services/history-snapshot.js'
    ]),
    operatorPattern: /\b(?:createHistoryManager|beginCheckpoint|commitCheckpoint|buildArtifactCheckpoint)\b/g
  }),
  Object.freeze({
    name: 'Session',
    importTargets: Object.freeze([
      'services/config.js',
      'services/gallery-session-migration.js',
      'services/session-authority.js',
      'services/session-file.js'
    ]),
    operatorPattern: /\b(?:exportSession|compressSessionData|migrateSessionDataToCurrent|migratePersistedGalleryConfig)\b/g
  }),
  Object.freeze({
    name: 'Canonical editor state',
    importTargets: Object.freeze([
      'app/feature-editor.js',
      'app/legend.js',
      'app/right-drawer.js',
      'app/preview-runtime.js'
    ]),
    operatorPattern: /\b(?:createFeatureEditor|createLegendManager|createRightDrawerController|createPreviewRuntime)\b/g
  })
]);

export const WEB_PRIVILEGED_CAPABILITY_KEYS = Object.freeze(
  PRIVILEGED_CAPABILITY_SPECS.map(({ name }) => name).sort()
);
export const WEB_PRIVILEGED_IMPORT_TARGETS = Object.freeze(
  unique(PRIVILEGED_CAPABILITY_SPECS.flatMap(({ importTargets }) => importTargets)).sort()
);

export const detectPrivilegedWebCapabilities = (sources) => {
  const entries = sourceEntries(sources);
  const importerTargets = new Set(WEB_PRIVILEGED_IMPORT_TARGETS);
  const importersByTarget = Object.fromEntries(
    WEB_PRIVILEGED_IMPORT_TARGETS.map((target) => [target, []])
  );
  const operatorMatchesByCapability = Object.fromEntries(
    WEB_PRIVILEGED_CAPABILITY_KEYS.map((capability) => [capability, []])
  );

  entries.forEach(([path, source]) => {
    unique(literalImportSpecifiers(source)
      .map((specifier) => resolvedProductionImport(path, specifier))
      .filter((target) => importerTargets.has(target)))
      .forEach((target) => importersByTarget[target].push(path));

    const code = maskJavaScript(source);
    PRIVILEGED_CAPABILITY_SPECS.forEach(({ name, operatorPattern }) => {
      const count = matches(code, operatorPattern).length;
      if (count) operatorMatchesByCapability[name].push(Object.freeze({ path, count }));
    });
  });

  Object.values(importersByTarget).forEach(Object.freeze);
  Object.values(operatorMatchesByCapability).forEach(Object.freeze);
  return Object.freeze({
    importersByTarget: Object.freeze(importersByTarget),
    operatorMatchesByCapability: Object.freeze(operatorMatchesByCapability)
  });
};

const SESSION_SOURCE_PATHS = new Set([
  'services/config.js',
  'services/gallery-session-migration.js',
  'services/session-authority.js',
  'services/session-file.js',
  'services/session-request.js'
]);

export const isWebSessionSourcePath = (path) => SESSION_SOURCE_PATHS.has(
  normalizeModulePath(path)
);

export const detectReportOnlySourceFacts = (source = '') => {
  const code = maskJavaScript(source);
  const exportedNames = new Set();
  for (const match of code.matchAll(
    /^\s*export\s+(?:async\s+)?(?:const|let|var|function|class)\s+([A-Za-z_$][\w$]*)/gm
  )) exportedNames.add(match[1]);
  for (const match of code.matchAll(/^\s*export\s*\{([\s\S]*?)\}/gm)) {
    match[1].split(',').forEach((entry) => {
      const cleaned = entry.replace(/\/\*[\s\S]*?\*\//g, '').trim();
      if (!cleaned) return;
      const parts = cleaned.split(/\s+as\s+/);
      exportedNames.add((parts[1] || parts[0]).trim());
    });
  }
  if (/^\s*export\s+default\b/m.test(code)) exportedNames.add('default');
  for (const match of code.matchAll(
    /^\s*export\s+\*\s+as\s+([A-Za-z_$][\w$]*)/gm
  )) exportedNames.add(match[1]);

  return Object.freeze({
    exportedNames: Object.freeze([...exportedNames]),
    declaredNames: Object.freeze(unique(
      [...code.matchAll(/\b(?:const|let|var|function|class)\s+([A-Za-z_$][\w$]*)/g)]
        .map((match) => match[1])
    )),
    reactiveDeclarations: Object.freeze(unique(
      [...code.matchAll(/\b(?:const|let|var)\s+([A-Za-z_$][\w$]*)\s*=\s*(?:[A-Za-z_$][\w$]*\.)?(ref|shallowRef|reactive|computed)\s*\(/g)]
        .map((match) => `${match[1]} (${match[2]})`)
    )),
    watcherCount: matches(code, /\bwatch(?:Effect)?\s*\(/g).length,
    objectKeys: Object.freeze(unique(
      [...code.matchAll(/(?:^|[,{]\s*)([A-Za-z_$][\w$]*)\s*:/gm)]
        .map((match) => match[1])
    )),
    compatibilityNames: Object.freeze(unique(
      [...code.matchAll(/\b(?:[A-Za-z_$][\w$]*(?:Migration|Migrator|Legacy|Fallback)[\w$]*|(?:migrate|promoteLegacy|normalizeLegacy|readLegacy)[A-Za-z0-9_$]*)\b/g)]
        .map((match) => match[0])
    )),
    namedResourceNames: Object.freeze(unique(
      [...code.matchAll(/\b[A-Za-z_$][\w$]*\b/g)]
        .map((match) => match[0])
        .filter((name) => /(?:cache|token|handle|journal|protocol|manager)/i.test(name))
    )),
    importSpecifiers: Object.freeze(literalImportSpecifiers(source))
  });
};

const RENDER_REQUEST_DEFINITION_PATTERN = /\b(?:(?:export\s+)?(?:async\s+)?function\s+buildCanonicalRenderRequest\b|(?:export\s+)?(?:const|let|var)\s+buildCanonicalRenderRequest\s*=)/g;
const RENDER_REQUEST_CALL_PATTERN = /\bbuildCanonicalRenderRequest\s*\(/g;
const RENDER_CALL_PATTERN = /\brunDiagramGeneration\s*\(/g;
const encodeDefinitionPathSubject = ({ path }) => normalizeModulePath(path);
const encodeCanonicalEntryEdgeSubject = ({ from, to }) => (
  `${normalizeModulePath(from)} -> ${normalizeModulePath(to)}`
);

const detectRenderRequestAuthority = (sources) => {
  const observedDefinitions = sourceEntries(sources).flatMap(([path, source]) => {
    const count = matches(maskJavaScript(source), RENDER_REQUEST_DEFINITION_PATTERN).length;
    return count ? [Object.freeze({ path, count, subject: encodeDefinitionPathSubject({ path }) })] : [];
  });
  return Object.freeze({
    definitionCount: observedDefinitions.reduce((total, { count }) => total + count, 0),
    observedDefinitions: Object.freeze(observedDefinitions),
    subjects: Object.freeze(observedDefinitions.map(({ subject }) => subject))
  });
};

const detectCanonicalRenderRequestPath = (sources) => {
  const entries = sourceEntries(sources);
  const definitionPaths = new Set(
    detectRenderRequestAuthority(sources).observedDefinitions.map(({ path }) => path)
  );
  const edgesBySubject = new Map();
  entries.forEach(([from, source]) => {
    const code = maskJavaScript(source);
    if (!matches(code, RENDER_REQUEST_CALL_PATTERN).length
        || !matches(code, RENDER_CALL_PATTERN).length) return;
    unique(literalImportSpecifiers(source)
      .map((specifier) => resolvedProductionImport(from, specifier))
      .filter((to) => definitionPaths.has(to)))
      .forEach((to) => {
        const edge = Object.freeze({ from, to });
        edgesBySubject.set(encodeCanonicalEntryEdgeSubject(edge), edge);
      });
  });
  const subjects = [...edgesBySubject.keys()].sort();
  return Object.freeze({
    observedEdges: Object.freeze(subjects.map((subject) => Object.freeze({
      ...edgesBySubject.get(subject),
      subject
    }))),
    subjects: Object.freeze(subjects)
  });
};

const CURRENT_RESULT_ADMISSION_NAME = 'admitCurrentGeneratedResults';
const CURRENT_RESULT_ADMISSION_DEFINITION_PATTERN = new RegExp(
  `\\bexport\\s+const\\s+${CURRENT_RESULT_ADMISSION_NAME}`
    + '\\s*=\\s*\\([^)]*\\)\\s*=>',
  'g'
);

const currentResultAdmissionNamedImportsV1 = (source) => {
  const commentsMasked = maskJavaScript(source, { strings: false });
  const code = maskJavaScript(source);
  const imports = [];
  const pattern = /(?:^|\n)\s*import\s*\{([\s\S]*?)\}\s*from\s*(['"])([^'"\r\n]+)\2\s*;?/g;
  for (const match of commentsMasked.matchAll(pattern)) {
    const keywordOffset = match[0].search(/\bimport\b/);
    const keywordIndex = match.index + keywordOffset;
    if (!/^import\b/.test(code.slice(keywordIndex))) continue;
    match[1].split(',').forEach((entry) => {
      const binding = entry.trim().match(
        new RegExp(
          `^${CURRENT_RESULT_ADMISSION_NAME}`
            + '(?:\\s+as\\s+([A-Za-z_$][\\w$]*))?$'
        )
      );
      if (!binding) return;
      imports.push(Object.freeze({
        specifier: match[3],
        local: binding[1] || CURRENT_RESULT_ADMISSION_NAME
      }));
    });
  }
  return imports;
};

const currentResultAdmissionCallPatternV1 = (local) => new RegExp(
  `(?:^|[^A-Za-z0-9_$.])${local.replace(/[.*+?^${}()|[\]\\]/g, '\\$&')}\\s*\\(`,
  'g'
);

const detectCurrentResultAdmissionAuthorityV1 = (sources) => {
  const observedDefinitions = sourceEntries(sources).flatMap(([path, source]) => {
    const count = matches(
      maskJavaScript(source),
      CURRENT_RESULT_ADMISSION_DEFINITION_PATTERN
    ).length;
    return count ? [Object.freeze({
      path,
      count,
      subject: encodeDefinitionPathSubject({ path })
    })] : [];
  });
  return Object.freeze({
    definitionCount: observedDefinitions.reduce((total, { count }) => total + count, 0),
    observedDefinitions: Object.freeze(observedDefinitions),
    subjects: Object.freeze(observedDefinitions.map(({ subject }) => subject))
  });
};

const detectCanonicalCurrentResultAdmissionPathV1 = (sources) => {
  const definitionPaths = new Set(
    detectCurrentResultAdmissionAuthorityV1(sources).observedDefinitions
      .map(({ path }) => path)
  );
  const edgesBySubject = new Map();
  sourceEntries(sources).forEach(([from, source]) => {
    const code = maskJavaScript(source);
    currentResultAdmissionNamedImportsV1(source)
      .map(({ specifier, local }) => ({
        local,
        to: resolvedProductionImport(from, specifier)
      }))
      .filter(({ to }) => definitionPaths.has(to))
      .filter(({ local }) => matches(
        code,
        currentResultAdmissionCallPatternV1(local)
      ).length > 0)
      .forEach(({ to }) => {
        const edge = Object.freeze({ from, to });
        edgesBySubject.set(encodeCanonicalEntryEdgeSubject(edge), edge);
      });
  });
  const subjects = [...edgesBySubject.keys()].sort();
  return Object.freeze({
    observedEdges: Object.freeze(subjects.map((subject) => Object.freeze({
      ...edgesBySubject.get(subject),
      subject
    }))),
    subjects: Object.freeze(subjects)
  });
};

// Keep each detector ID and its transitive helpers executable-identical while referenced.
// Add a versioned detector beside it, then migrate authority in a separate pull request.
export const WEB_ARCHITECTURE_DETECTORS = Object.freeze({
  'semantic-owner.render-request.v1': Object.freeze({
    subjectCategory: 'definition-path',
    encodeSubject: encodeDefinitionPathSubject,
    detect: detectRenderRequestAuthority
  }),
  'canonical-path.render-request.v1': Object.freeze({
    subjectCategory: 'canonical-entry-edge',
    encodeSubject: encodeCanonicalEntryEdgeSubject,
    detect: detectCanonicalRenderRequestPath
  }),
  'semantic-owner.current-result-admission.v1': Object.freeze({
    subjectCategory: 'definition-path',
    encodeSubject: encodeDefinitionPathSubject,
    detect: detectCurrentResultAdmissionAuthorityV1
  }),
  'canonical-path.current-result-admission.v1': Object.freeze({
    subjectCategory: 'canonical-entry-edge',
    encodeSubject: encodeCanonicalEntryEdgeSubject,
    detect: detectCanonicalCurrentResultAdmissionPathV1
  })
});
