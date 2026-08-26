#!/usr/bin/env node

import { appendFileSync, existsSync, readFileSync } from 'node:fs';
import { join, posix, resolve } from 'node:path';
import { execFileSync, spawnSync } from 'node:child_process';
import {
  detectPrivilegedWebCapabilities,
  detectReportOnlySourceFacts,
  isWebSessionSourcePath,
  WEB_ARCHITECTURE_DETECTORS,
  WEB_PRIVILEGED_CAPABILITY_KEYS,
  WEB_PRIVILEGED_IMPORT_TARGETS
} from './web-architecture-detectors.mjs';
import {
  classifyArchitectureAuthorityDelta,
  classifyArchitectureRuleObservation,
  evaluateArchitectureRuleResult,
  findDirectedGraphCycles,
  summarizeArchitectureInventory,
  validateArchitectureRuleRegistry
} from './web-architecture-evaluation.mjs';
import { literalImportSpecifiers, maskJavaScript } from './web-change-source.mjs';
import {
  classifyWebChangeContext,
  isPullRequestEvent
} from './web-promotion-context.mjs';

const runGit = (root, args, options = {}) => execFileSync(
  'git',
  ['-C', root, ...args],
  { encoding: 'utf8', maxBuffer: 32 * 1024 * 1024, ...options }
);

const repositoryRoot = runGit(process.cwd(), ['rev-parse', '--show-toplevel']).trim();
const argumentsByName = new Map();
for (let index = 2; index < process.argv.length; index += 1) {
  const argument = process.argv[index];
  if (!argument.startsWith('--')) throw new Error(`Unknown argument: ${argument}`);
  const value = process.argv[index + 1];
  if (!value || value.startsWith('--')) throw new Error(`Missing value for ${argument}`);
  argumentsByName.set(argument, value);
  index += 1;
}

const base = argumentsByName.get('--base') || process.env.WEB_CHANGE_BASE || 'HEAD';
const head = argumentsByName.get('--head') || process.env.WEB_CHANGE_HEAD || '';
const architectureChange = process.env.WEB_ARCHITECTURE_CHANGE === 'true';
const budgetProfiles = Object.freeze({
  ordinary: Object.freeze({ productionFiles: 8, grossChurn: 800, netAdditions: 100 }),
  architecture: Object.freeze({ productionFiles: 12, grossChurn: 1500, netAdditions: 400 })
});
const selectedProfileName = architectureChange ? 'architecture' : 'ordinary';
const selectedProfile = budgetProfiles[selectedProfileName];
const diffRefs = head ? [base, head] : [base];
const githubEventName = process.env.GITHUB_EVENT_NAME || '';
const githubEventPath = process.env.GITHUB_EVENT_PATH || '';
const githubEventPayloadSource = isPullRequestEvent(githubEventName)
  && githubEventPath
  && existsSync(githubEventPath)
  ? readFileSync(githubEventPath, 'utf8')
  : '';
const changeContext = classifyWebChangeContext({
  eventName: githubEventName,
  currentRepository: process.env.GITHUB_REPOSITORY || '',
  eventPayloadSource: githubEventPayloadSource,
  baseSha: base,
  headSha: head
});

const promotionSourceAncestry = (() => {
  if (!changeContext.isPromotion) return { status: 'NOT_APPLICABLE', violation: '' };
  const result = spawnSync(
    'git',
    ['-C', repositoryRoot, 'merge-base', '--is-ancestor', base, head],
    { encoding: 'utf8', stdio: ['ignore', 'pipe', 'pipe'] }
  );
  if (result.status === 0) return { status: 'PASS', violation: '' };
  if (result.status === 1) {
    return {
      status: 'FAIL',
      violation: (
        'The promotion source does not contain the current main head. '
        + 'Merge or rebase main into dev, then rerun the promotion.'
      )
    };
  }
  return {
    status: 'ERROR',
    violation: (
      'Promotion source ancestry could not be verified. '
      + 'Fetch complete base and head history, then rerun the promotion.'
    )
  };
})();

const parseDiffLines = (output) => output.trim()
  ? output.trimEnd().split('\n').map((line) => line.split('\t'))
  : [];

const changed = new Map(
  parseDiffLines(runGit(repositoryRoot, [
    'diff', '--no-renames', '--name-status', ...diffRefs, '--'
  ])).map(([status, path]) => [path, status])
);

if (!head) {
  const untracked = runGit(repositoryRoot, [
    'ls-files', '--others', '--exclude-standard', '-z'
  ]).split('\0').filter(Boolean);
  untracked.forEach((path) => changed.set(path, 'A'));
}

const numstat = new Map(
  parseDiffLines(runGit(repositoryRoot, [
    'diff', '--no-renames', '--numstat', ...diffRefs, '--'
  ])).map(([additions, deletions, path]) => [path, { additions, deletions }])
);

const readRevisionFile = (revision, path) => {
  try {
    return runGit(repositoryRoot, ['show', `${revision}:${path}`], { stdio: ['ignore', 'pipe', 'ignore'] });
  } catch (_error) {
    return null;
  }
};

const readRevisionFiles = (revision, paths) => {
  if (!paths.length) return new Map();
  const input = `${paths.map((path) => `${revision}:${path}`).join('\n')}\n`;
  const output = execFileSync(
    'git',
    ['-C', repositoryRoot, 'cat-file', '--batch'],
    { input, maxBuffer: 32 * 1024 * 1024 }
  );
  const sources = new Map();
  let offset = 0;
  paths.forEach((path) => {
    const headerEnd = output.indexOf(10, offset);
    const header = output.subarray(offset, headerEnd).toString('utf8');
    const size = Number(header.split(' ').at(-1));
    if (!Number.isInteger(size)) throw new Error(`Cannot read ${path} at ${revision}: ${header}`);
    const contentStart = headerEnd + 1;
    const contentEnd = contentStart + size;
    sources.set(path, output.subarray(contentStart, contentEnd).toString('utf8'));
    offset = contentEnd + 1;
  });
  return sources;
};

const readHeadFile = (path) => {
  if (head) return readRevisionFile(head, path);
  const absolutePath = join(repositoryRoot, path);
  return existsSync(absolutePath) ? readFileSync(absolutePath, 'utf8') : null;
};

const lineCount = (source) => {
  if (!source) return 0;
  const lines = source.split(/\r?\n/);
  return lines.length - (lines.at(-1) === '' ? 1 : 0);
};

const isRuntimePath = (path) => path === 'gbdraw/web/index.html'
  || path.startsWith('gbdraw/web/js/')
  || path.startsWith('gbdraw/web/vendor/');
const looksBinary = (path) => {
  if (head) return false;
  const absolutePath = join(repositoryRoot, path);
  return existsSync(absolutePath) && readFileSync(absolutePath).includes(0);
};

for (const [path, status] of changed) {
  if (status !== 'A' || numstat.has(path) || !isRuntimePath(path)) continue;
  numstat.set(path, looksBinary(path)
    ? { additions: '-', deletions: '-' }
    : { additions: String(lineCount(readHeadFile(path))), deletions: '0' });
}

const productionPaths = [...changed.keys()]
  .filter(isRuntimePath)
  .sort();
const productionJavaScriptPaths = productionPaths.filter((path) => /\.[cm]?js$/.test(path));
const productionTotals = productionPaths.reduce((totals, path) => {
  const counts = numstat.get(path) || { additions: '0', deletions: '0' };
  if (counts.additions === '-' || counts.deletions === '-') totals.binary += 1;
  else {
    totals.additions += Number(counts.additions);
    totals.deletions += Number(counts.deletions);
  }
  return totals;
}, { additions: 0, deletions: 0, binary: 0 });
const productionGrossChurn = productionTotals.additions + productionTotals.deletions;
const productionNetAdditions = productionTotals.additions - productionTotals.deletions;

const addedBinaryRuntimePaths = productionPaths.filter((path) => {
  const counts = numstat.get(path);
  return changed.get(path) === 'A' && counts?.additions === '-' && counts?.deletions === '-';
});
const changedVendorPaths = productionPaths.filter((path) => path.startsWith('gbdraw/web/vendor/'));
const policyPath = 'tools/web-change-policy.json';
const architectureRulesPath = 'tools/web-architecture-rules.json';
const guardPaths = new Set([
  'docs/internal/ARCHITECTURE_FITNESS_FUNCTION_RATCHET.md',
  '.github/pull_request_template.md',
  'tools/check-web-change-budget.mjs',
  'tools/web-architecture-detectors.mjs',
  'tools/web-architecture-evaluation.mjs',
  architectureRulesPath,
  'tools/web-architecture-violations.json',
  'tools/web-change-source.mjs',
  'tools/web-promotion-context.mjs',
  policyPath,
  'docs/internal/WEB_CHANGE_POLICY.md',
  'tests/web/architecture-contracts.test.mjs',
  'tests/web/architecture-ratchet-fixtures.test.mjs',
  'tests/web/web-promotion-context.test.mjs',
  '.github/workflows/test.yml',
  '.github/workflows/web-base-policy.yml'
]);
const changedGuards = [...changed.keys()].filter((path) => guardPaths.has(path)).sort();
const checkerImplementationPaths = new Set([
  'tools/check-web-change-budget.mjs',
  'tools/web-architecture-detectors.mjs',
  'tools/web-architecture-evaluation.mjs',
  'tools/web-change-source.mjs',
  'tools/web-promotion-context.mjs'
]);
const authorityPaths = new Set([
  'docs/internal/ARCHITECTURE_FITNESS_FUNCTION_RATCHET.md',
  architectureRulesPath,
  'tools/web-architecture-violations.json',
  policyPath,
  '.github/workflows/test.yml',
  '.github/workflows/web-base-policy.yml'
]);
const changedCheckerImplementations = [...changed.keys()]
  .filter((path) => checkerImplementationPaths.has(path));
const changedAuthorities = [...changed.keys()].filter((path) => authorityPaths.has(path));

const listProductionJavaScriptPaths = () => {
  if (head) {
    return runGit(repositoryRoot, [
      'ls-tree', '-r', '--name-only', head, '--', 'gbdraw/web/js'
    ]).split('\n').filter((path) => /\.[cm]?js$/.test(path));
  }
  const tracked = runGit(repositoryRoot, [
    'ls-files', '--cached', '--modified', '--others', '--exclude-standard', '--', 'gbdraw/web/js'
  ]).split('\n').filter(Boolean);
  return [...new Set(tracked)]
    .filter((path) => /\.[cm]?js$/.test(path) && existsSync(join(repositoryRoot, path)))
    .sort();
};

const listRevisionProductionJavaScriptPaths = (revision) => runGit(repositoryRoot, [
  'ls-tree', '-r', '--name-only', revision, '--', 'gbdraw/web/js'
]).split('\n').filter((path) => /\.[cm]?js$/.test(path));
const allProductionJavaScriptPaths = listProductionJavaScriptPaths();
const headRevisionSources = head
  ? readRevisionFiles(head, allProductionJavaScriptPaths)
  : null;
const allProductionSources = new Map(allProductionJavaScriptPaths.map((path) => [
  path.slice('gbdraw/web/js/'.length),
  headRevisionSources?.get(path) ?? readHeadFile(path) ?? ''
]));
const baseProductionJavaScriptPaths = listRevisionProductionJavaScriptPaths(base);
const changedBaseSources = readRevisionFiles(
  base,
  baseProductionJavaScriptPaths.filter((path) => (
    changed.has(path) || !allProductionSources.has(path.slice('gbdraw/web/js/'.length))
  ))
);
const baseProductionSources = new Map(
  baseProductionJavaScriptPaths.map((path) => [
    path.slice('gbdraw/web/js/'.length),
    changed.has(path)
      ? changedBaseSources.get(path) || ''
      : allProductionSources.get(path.slice('gbdraw/web/js/'.length))
        ?? changedBaseSources.get(path)
        ?? ''
  ])
);

const sourceInventory = (sources) => {
  const inventory = {
    modules: [],
    exports: [],
    createDeclarations: [],
    reactiveDeclarations: [],
    watcherCalls: [],
    compatibilityNames: [],
    resourceNames: [],
    sessionObjectKeys: [],
    bareImports: []
  };
  [...sources].sort(([left], [right]) => left.localeCompare(right))
    .forEach(([path, source]) => {
      const displayPath = `gbdraw/web/js/${path}`;
      const facts = detectReportOnlySourceFacts(source);
      inventory.modules.push(displayPath);
      facts.exportedNames.forEach((name) => inventory.exports.push(`${displayPath}: ${name}`));
      facts.declaredNames.filter((name) => /^create[A-Z]/.test(name))
        .forEach((name) => inventory.createDeclarations.push(`${displayPath}: ${name}`));
      facts.reactiveDeclarations.forEach((name) => {
        inventory.reactiveDeclarations.push(`${displayPath}: ${name}`);
      });
      for (let index = 1; index <= facts.watcherCount; index += 1) {
        inventory.watcherCalls.push(`${displayPath}: watcher call ${index}`);
      }
      facts.compatibilityNames.forEach((name) => {
        inventory.compatibilityNames.push(`${displayPath}: ${name}`);
      });
      facts.namedResourceNames.forEach((name) => {
        inventory.resourceNames.push(`${displayPath}: ${name}`);
      });
      if (isWebSessionSourcePath(path)) {
        facts.objectKeys.forEach((name) => {
          inventory.sessionObjectKeys.push(`${displayPath}: ${name}`);
        });
      }
      facts.importSpecifiers.filter((specifier) => !specifier.startsWith('.'))
        .forEach((specifier) => inventory.bareImports.push(`${displayPath}: ${specifier}`));
    });
  return inventory;
};

const baseSourceInventory = sourceInventory(baseProductionSources);
const headSourceInventory = sourceInventory(allProductionSources);
const sourceInventoryDeltas = Object.fromEntries(
  Object.keys(baseSourceInventory).map((name) => [
    name,
    summarizeArchitectureInventory(baseSourceInventory[name], headSourceInventory[name])
  ])
);
const newModules = sourceInventoryDeltas.modules.added;
const newExports = sourceInventoryDeltas.exports.added;
const newCreateOwners = sourceInventoryDeltas.createDeclarations.added;
const newReactiveState = sourceInventoryDeltas.reactiveDeclarations.added;
const newWatcherCountsByPath = new Map();
sourceInventoryDeltas.watcherCalls.added.forEach((entry) => {
  const path = entry.slice(0, entry.lastIndexOf(': watcher call '));
  newWatcherCountsByPath.set(path, (newWatcherCountsByPath.get(path) || 0) + 1);
});
const newWatchers = [...newWatcherCountsByPath]
  .sort(([left], [right]) => left.localeCompare(right))
  .map(([path, count]) => `${path}: +${count} watcher call(s)`);
const newNamedResources = sourceInventoryDeltas.resourceNames.added;
const newCompatibilityPaths = sourceInventoryDeltas.compatibilityNames.added;
const newSessionFields = sourceInventoryDeltas.sessionObjectKeys.added;
const newBareImports = sourceInventoryDeltas.bareImports.added;

const listRevisionProductionPaths = (revision) => runGit(repositoryRoot, [
  'ls-tree', '-r', '--name-only', revision, '--',
  'gbdraw/web/index.html', 'gbdraw/web/js', 'gbdraw/web/vendor'
]).split('\n').filter((path) => path && isRuntimePath(path));
const listHeadProductionPaths = () => {
  if (head) return listRevisionProductionPaths(head);
  return runGit(repositoryRoot, [
    'ls-files', '--cached', '--modified', '--others', '--exclude-standard', '--',
    'gbdraw/web/index.html', 'gbdraw/web/js', 'gbdraw/web/vendor'
  ]).split('\n').filter((path) => (
    path && isRuntimePath(path) && existsSync(join(repositoryRoot, path))
  ));
};
const productionFileDelta = summarizeArchitectureInventory(
  listRevisionProductionPaths(base),
  listHeadProductionPaths()
);

const staticImportSpecifiers = (source) => {
  const candidates = new Set(literalImportSpecifiers(source, { dynamic: false }));
  const commentsMasked = maskJavaScript(source, { strings: false });
  const code = maskJavaScript(source);
  const specifiers = [];
  const patterns = [
    /(?:^|\n)\s*import\s+(?:[^;]*?\s+from\s+)?(['"])([^'"]+)\1\s*;?/g,
    /(?:^|\n)\s*export\s+(?:\*|\{)[^;]*?\s+from\s+(['"])([^'"]+)\1\s*;?/g
  ];
  patterns.forEach((pattern) => {
    for (const match of commentsMasked.matchAll(pattern)) {
      const keywordOffset = match[0].search(/\b(?:import|export)\b/);
      const keywordIndex = match.index + keywordOffset;
      if (
        candidates.has(match[2])
        && /^(?:import|export)\b/.test(code.slice(keywordIndex))
      ) specifiers.push(match[2]);
    }
  });
  return [...new Set(specifiers)];
};

const buildFirstPartyStaticImportGraph = (sources) => {
  const nodes = [...sources.keys()].filter((path) => path.endsWith('.js')).sort();
  const nodeSet = new Set(nodes);
  const edgesBySubject = new Map();
  const errors = [];
  nodes.forEach((owner) => {
    staticImportSpecifiers(sources.get(owner))
      .filter((specifier) => specifier.startsWith('.'))
      .forEach((specifier) => {
        const unresolved = posix.normalize(posix.join(posix.dirname(owner), specifier));
        if (unresolved === '..' || unresolved.startsWith('../') || posix.isAbsolute(unresolved)) {
          errors.push(`${owner}: static import ${specifier} escapes gbdraw/web/js/`);
          return;
        }
        const candidates = posix.extname(unresolved)
          ? [unresolved]
          : [`${unresolved}.js`, posix.join(unresolved, 'index.js')];
        const target = candidates.find((candidate) => nodeSet.has(candidate));
        if (!target) {
          errors.push(`${owner}: cannot resolve first-party static import ${specifier}`);
          return;
        }
        edgesBySubject.set(`${owner}\u0000${target}`, Object.freeze({ from: owner, to: target }));
      });
  });
  return Object.freeze({
    graph: Object.freeze({ nodes, edges: [...edgesBySubject.values()] }),
    errors: Object.freeze([...new Set(errors)].sort())
  });
};

const baseImportGraph = buildFirstPartyStaticImportGraph(baseProductionSources);
const headImportGraph = buildFirstPartyStaticImportGraph(allProductionSources);
const baseCycleResult = findDirectedGraphCycles(baseImportGraph.graph);
const headCycleResult = findDirectedGraphCycles(headImportGraph.graph);
const cycleDelta = summarizeArchitectureInventory(
  baseCycleResult.cycles.map(({ subject }) => subject),
  headCycleResult.cycles.map(({ subject }) => subject)
);

const parsePolicy = (source, revision) => {
  if (source === null) throw new Error(`Missing ${policyPath} at ${revision}`);
  let policy;
  try {
    policy = JSON.parse(source);
  } catch (error) {
    throw new Error(`Cannot parse ${policyPath} at ${revision}: ${error.message}`);
  }
  for (const section of ['allowedPrivilegedImporters', 'allowedPrivilegedOwners']) {
    if (!policy[section] || typeof policy[section] !== 'object' || Array.isArray(policy[section])) {
      throw new Error(`${policyPath} must define an object named ${section}`);
    }
    Object.entries(policy[section]).forEach(([name, paths]) => {
      if (!Array.isArray(paths) || paths.some((path) => typeof path !== 'string')) {
        throw new Error(`${policyPath} ${section}.${name} must be an array of paths`);
      }
    });
  }
  return policy;
};

const basePolicySource = readRevisionFile(base, policyPath);
const policySource = basePolicySource ?? readHeadFile(policyPath);
const policyRevision = basePolicySource === null ? `${head || 'working tree'} (bootstrap)` : base;
const policy = parsePolicy(policySource, policyRevision);
let proposedPolicy = null;
const policyContractions = [];
const policyExpansions = [];
const missingPolicyKeys = [];
const addedPolicyKeys = [];
let policyTopLevelKeysMatch = true;
if (basePolicySource !== null && changed.has(policyPath)) {
  proposedPolicy = parsePolicy(readHeadFile(policyPath), head || 'working tree');
  const baseTopLevelKeys = Object.keys(policy).sort();
  const proposedTopLevelKeys = Object.keys(proposedPolicy).sort();
  policyTopLevelKeysMatch = JSON.stringify(baseTopLevelKeys) === JSON.stringify(proposedTopLevelKeys);
  for (const section of ['allowedPrivilegedImporters', 'allowedPrivilegedOwners']) {
    Object.keys(proposedPolicy[section]).forEach((name) => {
      if (!Object.hasOwn(policy[section], name)) addedPolicyKeys.push(`${section}.${name}`);
    });
    Object.entries(policy[section]).forEach(([name, allowedPaths]) => {
      if (!Object.hasOwn(proposedPolicy[section], name)) {
        missingPolicyKeys.push(`${section}.${name}`);
      }
      const proposedPaths = new Set(proposedPolicy[section][name] || []);
      allowedPaths.forEach((path) => {
        if (!proposedPaths.has(path)) policyContractions.push(`${section}.${name}: ${path}`);
      });
      const basePaths = new Set(allowedPaths);
      proposedPaths.forEach((path) => {
        if (!basePaths.has(path)) policyExpansions.push(`${section}.${name}: ${path}`);
      });
    });
  }
}
const proposedOrBasePolicy = proposedPolicy || policy;
const permissionEntries = (candidatePolicy, section) => Object.entries(candidatePolicy[section])
  .flatMap(([name, paths]) => paths.map((path) => `${name}: ${path}`));
const privilegedOperatorPermissionDelta = summarizeArchitectureInventory(
  permissionEntries(policy, 'allowedPrivilegedOwners'),
  permissionEntries(proposedOrBasePolicy, 'allowedPrivilegedOwners')
);
const privilegedImporterPermissionDelta = summarizeArchitectureInventory(
  permissionEntries(policy, 'allowedPrivilegedImporters'),
  permissionEntries(proposedOrBasePolicy, 'allowedPrivilegedImporters')
);
const emptyRuleRegistry = () => ({ schemaVersion: 1, rules: [] });
const candidateArchitectureErrors = [];
const candidateArchitectureObservations = [];
const activeArchitectureErrors = [];
const activeArchitectureFailures = [];
const activeArchitectureRows = [];
const trustedBaseArchitectureErrors = [];
const trustedBaseArchitectureFailures = [];
let registeredAuthorityLocationDelta = summarizeArchitectureInventory([], []);
let registeredCanonicalContractDelta = summarizeArchitectureInventory([], []);
let architectureAuthorityDelta = Object.freeze({
  classification: 'UNCHANGED',
  directions: Object.freeze([]),
  changes: Object.freeze([])
});

const parseArchitectureRuleRegistry = (source, revision) => {
  if (source === null) return { registry: emptyRuleRegistry(), error: null };
  try {
    return { registry: JSON.parse(source), error: null };
  } catch (error) {
    return {
      registry: null,
      error: `Cannot parse ${architectureRulesPath} at ${revision}: ${error.message}`
    };
  }
};

const architectureRuleValidationOptions = Object.freeze({
  availableEnforcements: Object.freeze(['hard', 'report-only'])
});
const formatRuleValidationErrors = (label, validation) => validation.errors.map((error) => (
  `${label} ${error.path} [${error.code}]: ${error.message}`
));
const intendedRuleSubjects = (rule, detector) => (
  rule.kind === 'single-semantic-owner'
    ? rule.allowedDefinitionPaths.map((path) => detector.encodeSubject({ path }))
    : rule.allowedEdges.map((edge) => detector.encodeSubject(edge))
).sort();
const registeredLocationSummary = (rule, detectorOutput, observation) => (
  rule.kind === 'single-semantic-owner'
    ? `registeredDefinitionLocations=${new Set(
      (detectorOutput.observedDefinitions || []).map(({ path }) => path)
    ).size}; observedDefinitions=${observation.observedCount}/${observation.requiredCount}`
    : `observedCanonicalEntryEdges=${observation.observedCount}/${observation.requiredCount}`
);

const evaluateArchitectureSource = (registry, sources) => {
  const errors = [];
  const failures = [];
  const rows = [];
  const authorityLocations = [];
  const canonicalContracts = [];
  [...registry.rules]
    .sort((left, right) => left.key.localeCompare(right.key))
    .forEach((rule) => {
      const detector = WEB_ARCHITECTURE_DETECTORS[rule.detector];
      try {
        const detectorOutput = detector.detect(sources);
        const observation = classifyArchitectureRuleObservation(rule, detectorOutput);
        const result = evaluateArchitectureRuleResult({
          observation: observation.observation,
          mode: rule.enforcement.replace('-', '_').toUpperCase(),
          baselineRelation: 'NOT_APPLICABLE',
          authorityResolution: 'NOT_APPLICABLE'
        });
        const subjects = observation.subjects.length
          ? [...observation.subjects].sort()
          : intendedRuleSubjects(rule, detector);
        const locationSummary = registeredLocationSummary(rule, detectorOutput, observation);
        (subjects.length ? subjects : ['<none>']).forEach((subject) => {
          rows.push(
            `${rule.key} | subject=${subject} | observation=${result.observation} `
            + `| mode=${result.mode} | baselineRelation=${result.baselineRelation} `
            + `| authorityResolution=${result.authorityResolution} `
            + `| decision=${result.decision} | ${locationSummary}`
          );
          if (rule.kind === 'single-canonical-entry-edge') {
            canonicalContracts.push(
              `${rule.key} | subject=${subject} | observation=${result.observation} `
              + `| decision=${result.decision}`
            );
          }
        });
        if (rule.kind === 'single-semantic-owner') {
          new Set((detectorOutput.observedDefinitions || []).map(({ path }) => path))
            .forEach((path) => authorityLocations.push(`${rule.key}: ${path}`));
        }
        if (result.decision === 'FAIL') {
          failures.push(`${rule.key} is ${result.observation} under ${result.mode}`);
        }
      } catch (error) {
        errors.push(`${rule.key}: ${error.message}`);
      }
    });
  return { errors, failures, rows, authorityLocations, canonicalContracts };
};

const baseArchitectureRulesSource = readRevisionFile(base, architectureRulesPath);
const proposedArchitectureRulesSource = readHeadFile(architectureRulesPath);
if (baseArchitectureRulesSource !== null || proposedArchitectureRulesSource !== null) {
  const parsedBase = parseArchitectureRuleRegistry(baseArchitectureRulesSource, base);
  const parsedCandidate = parseArchitectureRuleRegistry(
    proposedArchitectureRulesSource,
    head || 'working tree'
  );
  if (parsedBase.error) activeArchitectureErrors.push(parsedBase.error);
  if (parsedCandidate.error) candidateArchitectureErrors.push(parsedCandidate.error);

  const baseValidation = parsedBase.registry
    ? validateArchitectureRuleRegistry(
      parsedBase.registry,
      WEB_ARCHITECTURE_DETECTORS,
      architectureRuleValidationOptions
    )
    : null;
  const candidateValidation = parsedCandidate.registry
    ? validateArchitectureRuleRegistry(
      parsedCandidate.registry,
      WEB_ARCHITECTURE_DETECTORS,
      architectureRuleValidationOptions
    )
    : null;
  if (baseValidation) {
    activeArchitectureErrors.push(
      ...formatRuleValidationErrors('base registry', baseValidation)
    );
  }
  if (candidateValidation) {
    candidateArchitectureErrors.push(
      ...formatRuleValidationErrors('candidate registry', candidateValidation)
    );
  }

  if (baseValidation?.valid) {
    const beforeArchitecture = evaluateArchitectureSource(
      parsedBase.registry,
      baseProductionSources
    );
    const afterArchitecture = evaluateArchitectureSource(
      parsedBase.registry,
      allProductionSources
    );
    trustedBaseArchitectureErrors.push(...beforeArchitecture.errors);
    trustedBaseArchitectureFailures.push(...beforeArchitecture.failures);
    activeArchitectureErrors.push(...afterArchitecture.errors);
    activeArchitectureFailures.push(...afterArchitecture.failures);
    activeArchitectureRows.push(...afterArchitecture.rows);
    registeredAuthorityLocationDelta = summarizeArchitectureInventory(
      beforeArchitecture.authorityLocations,
      afterArchitecture.authorityLocations
    );
    registeredCanonicalContractDelta = summarizeArchitectureInventory(
      beforeArchitecture.canonicalContracts,
      afterArchitecture.canonicalContracts
    );
  }

  if (baseValidation?.valid && candidateValidation?.valid) {
    architectureAuthorityDelta = classifyArchitectureAuthorityDelta(
      parsedBase.registry,
      parsedCandidate.registry
    );
    if (architectureAuthorityDelta.changes.length) {
      const companionPaths = [...changed.keys()]
        .filter((path) => path !== architectureRulesPath)
        .sort();
      if (companionPaths.length) {
        candidateArchitectureErrors.push(
          'architecture rule authority changes must be isolated from other changed paths: '
          + companionPaths.join(', ')
        );
      }

      const candidateRules = new Map(parsedCandidate.registry.rules.map((rule) => [rule.key, rule]));
      const changedRuleKeys = [...new Set(
        architectureAuthorityDelta.changes.map(({ rule }) => rule)
      )].sort();
      changedRuleKeys.forEach((ruleKey) => {
        const rule = candidateRules.get(ruleKey);
        if (!rule) return;
        const detector = WEB_ARCHITECTURE_DETECTORS[rule.detector];
        let baseObservation;
        try {
          baseObservation = classifyArchitectureRuleObservation(
            rule,
            detector.detect(baseProductionSources)
          );
        } catch (error) {
          candidateArchitectureErrors.push(
            `${rule.key} trusted-base detection failed: ${error.message}`
          );
          return;
        }
        candidateArchitectureObservations.push(
          `${rule.key} against untouched base: ${baseObservation.observation} `
          + `(${baseObservation.observedCount}/${baseObservation.requiredCount})`
        );
        if (
          (rule.enforcement === 'hard' || rule.baselineEligible === false)
          && baseObservation.observation !== 'CONFORMING'
        ) {
          candidateArchitectureErrors.push(
            `${rule.key} claims a clean base but is ${baseObservation.observation}`
          );
        }

        const directions = new Set(architectureAuthorityDelta.changes
          .filter((change) => change.rule === ruleKey)
          .map(({ direction }) => direction));
        if (![...directions].some((direction) => (
          direction === 'CONTRACTION' || direction === 'TIGHTENING'
        ))) return;
        let headObservation;
        try {
          headObservation = classifyArchitectureRuleObservation(
            rule,
            detector.detect(allProductionSources)
          );
        } catch (error) {
          candidateArchitectureErrors.push(
            `${rule.key} trusted-head-data detection failed: ${error.message}`
          );
          return;
        }
        candidateArchitectureObservations.push(
          `${rule.key} against head source data: ${headObservation.observation} `
          + `(${headObservation.observedCount}/${headObservation.requiredCount})`
        );
        if (headObservation.observation !== 'CONFORMING') {
          candidateArchitectureErrors.push(
            `${rule.key} authority contraction or tightening is ${headObservation.observation} `
            + 'on head source data'
          );
        }
      });
    }
  }
}

const baseDetectedCapabilities = detectPrivilegedWebCapabilities(baseProductionSources);
const detectedCapabilities = detectPrivilegedWebCapabilities(allProductionSources);
const privilegedImportEdges = (capabilities) => WEB_PRIVILEGED_IMPORT_TARGETS.flatMap((target) => (
  capabilities.importersByTarget[target].map((path) => `${path} -> ${target}`)
));
const privilegedImportFanOutDelta = summarizeArchitectureInventory(
  privilegedImportEdges(baseDetectedCapabilities),
  privilegedImportEdges(detectedCapabilities)
);
const privilegedFanOutCounts = (capabilities) => {
  const targetsByPath = new Map();
  WEB_PRIVILEGED_IMPORT_TARGETS.forEach((target) => {
    capabilities.importersByTarget[target].forEach((path) => {
      if (!targetsByPath.has(path)) targetsByPath.set(path, new Set());
      targetsByPath.get(path).add(target);
    });
  });
  return new Map([...targetsByPath]
    .map(([path, targets]) => [path, targets.size])
    .sort(([left], [right]) => left.localeCompare(right)));
};
const basePrivilegedFanOutCounts = privilegedFanOutCounts(baseDetectedCapabilities);
const headPrivilegedFanOutCounts = privilegedFanOutCounts(detectedCapabilities);
const capabilityCoverageViolations = (candidatePolicy) => {
  const violations = [];
  WEB_PRIVILEGED_IMPORT_TARGETS.forEach((target) => {
    const allowed = new Set(candidatePolicy.allowedPrivilegedImporters[target] || []);
    detectedCapabilities.importersByTarget[target].forEach((path) => {
      if (!allowed.has(path)) {
        violations.push(`${target}: importer ${path}`);
      }
    });
  });

  WEB_PRIVILEGED_CAPABILITY_KEYS.forEach((capability) => {
    const allowed = new Set(candidatePolicy.allowedPrivilegedOwners[capability] || []);
    detectedCapabilities.operatorMatchesByCapability[capability].forEach(({ path }) => {
      if (!allowed.has(path)) {
        violations.push(`${capability}: owner ${path}`);
      }
    });
  });
  return violations.sort();
};

const unapprovedCapabilities = capabilityCoverageViolations(policy);
const proposedPolicyExclusions = proposedPolicy
  ? capabilityCoverageViolations(proposedPolicy)
  : [];
const pureSafePolicyContraction = Boolean(
  productionPaths.length
  && proposedPolicy
  && changedGuards.length === 1
  && changedGuards[0] === policyPath
  && policyContractions.length
  && !policyExpansions.length
  && !missingPolicyKeys.length
  && !addedPolicyKeys.length
  && policyTopLevelKeysMatch
  && !unapprovedCapabilities.length
  && !proposedPolicyExclusions.length
);

const dependencySections = [
  'dependencies', 'optionalDependencies', 'peerDependencies', 'devDependencies'
];
const productionDependencySections = new Set(dependencySections.slice(0, 3));
const dependencyChanges = [];
const newProductionDependencies = [];
const manifestPaths = ['package.json', 'gbdraw/web/js/package.json'];

const parseManifest = (source, path, revision) => {
  if (source === null) return {};
  try {
    return JSON.parse(source);
  } catch (error) {
    throw new Error(`Cannot parse ${path} at ${revision}: ${error.message}`);
  }
};

manifestPaths.forEach((path) => {
  const before = parseManifest(readRevisionFile(base, path), path, base);
  const after = parseManifest(readHeadFile(path), path, head || 'working tree');
  dependencySections.forEach((section) => {
    const beforeEntries = before[section] || {};
    const afterEntries = after[section] || {};
    const names = new Set([...Object.keys(beforeEntries), ...Object.keys(afterEntries)]);
    [...names].sort().forEach((name) => {
      if (beforeEntries[name] === afterEntries[name]) return;
      const kind = !(name in beforeEntries) ? 'added' : !(name in afterEntries) ? 'removed' : 'changed';
      dependencyChanges.push(
        `${path} ${section}: ${kind} ${name}`
        + (kind === 'changed' ? ` (${beforeEntries[name]} -> ${afterEntries[name]})` : '')
      );
      if (kind === 'added' && productionDependencySections.has(section)) {
        newProductionDependencies.push(`${path} ${section}: ${name}`);
      }
    });
  });
});

const dependencyFilesTouched = [...changed.keys()]
  .filter((path) => [
    'package.json', 'package-lock.json', 'pyproject.toml', 'gbdraw/web/js/package.json'
  ].includes(path))
  .sort();
dependencyFilesTouched.forEach((path) => dependencyChanges.push(`file changed: ${path}`));
newProductionDependencies.push(
  ...newBareImports.map((entry) => `bare production import ${entry}`)
);
dependencyChanges.push(...newBareImports.map((entry) => `bare production import: ${entry}`));
dependencyChanges.push(...changedVendorPaths.map((path) => `vendored runtime file changed: ${path}`));

const budgetViolations = [];
if (productionPaths.length > selectedProfile.productionFiles) {
  budgetViolations.push(
    `production files changed exceed ${selectedProfileName} limit `
    + `(${productionPaths.length} > ${selectedProfile.productionFiles})`
  );
}
if (productionGrossChurn > selectedProfile.grossChurn) {
  budgetViolations.push(
    `production gross churn exceeds ${selectedProfileName} limit `
    + `(${productionGrossChurn} > ${selectedProfile.grossChurn})`
  );
}
if (productionNetAdditions > selectedProfile.netAdditions) {
  budgetViolations.push(
    `production net additions exceed ${selectedProfileName} limit `
    + `(${productionNetAdditions} > ${selectedProfile.netAdditions})`
  );
}

const integrityViolations = [];
const promotionAggregationObservations = [];
if (newProductionDependencies.length) {
  integrityViolations.push('new production dependencies are not allowed');
}
if (addedBinaryRuntimePaths.length) {
  integrityViolations.push('added binary runtime files are not allowed');
}
if (changedVendorPaths.length) {
  integrityViolations.push('changes under gbdraw/web/vendor/ are not allowed');
}
if (productionPaths.length && changedGuards.length && !pureSafePolicyContraction) {
  const reason = 'production runtime files and Web guard/CI files changed together';
  if (changeContext.isPromotion) {
    promotionAggregationObservations.push(
      `${reason}: production paths (${productionPaths.length}) ${productionPaths.join(', ')}; `
      + `guard paths (${changedGuards.length}) ${changedGuards.join(', ')}`
    );
  } else {
    integrityViolations.push(reason);
  }
}
if (changedCheckerImplementations.length && changedAuthorities.length) {
  const reason = 'Web checker/source parser and authority policy/workflow files changed together';
  if (changeContext.isPromotion) {
    promotionAggregationObservations.push(
      `${reason}: checker paths (${changedCheckerImplementations.length}) `
      + `${changedCheckerImplementations.sort().join(', ')}; authority paths `
      + `(${changedAuthorities.length}) ${changedAuthorities.sort().join(', ')}`
    );
  } else {
    integrityViolations.push(reason);
  }
}
if (unapprovedCapabilities.length) {
  integrityViolations.push('privileged capability owners or importers exceed the base allowlist');
}
if (proposedPolicyExclusions.length) {
  integrityViolations.push('proposed privileged capability allowlist excludes active owners or importers');
}
if (missingPolicyKeys.length) {
  integrityViolations.push('proposed privileged capability policy is missing base allowlist keys');
}
if (baseImportGraph.errors.length) {
  integrityViolations.push('trusted base first-party static import graph is incomplete');
}
if (headImportGraph.errors.length) {
  integrityViolations.push('head first-party static import graph is incomplete');
}
if (baseCycleResult.cycles.length) {
  integrityViolations.push('trusted base contains first-party static import cycles');
}
if (headCycleResult.cycles.length) {
  integrityViolations.push('first-party static import cycles are not allowed');
}
integrityViolations.push(...candidateArchitectureErrors.map((error) => (
  `candidate architecture rules: ${error}`
)));
integrityViolations.push(...trustedBaseArchitectureErrors.map((error) => (
  `trusted base architecture rules: ${error}`
)));
integrityViolations.push(...trustedBaseArchitectureFailures.map((error) => (
  `trusted base architecture rules: ${error}`
)));
integrityViolations.push(...activeArchitectureErrors.map((error) => (
  `active architecture rules: ${error}`
)));
integrityViolations.push(...activeArchitectureFailures.map((error) => (
  `active architecture rules: ${error}`
)));
const enforcedViolations = [
  ...changeContext.errors,
  ...(promotionSourceAncestry.violation ? [promotionSourceAncestry.violation] : []),
  ...integrityViolations,
  ...(changeContext.isPromotion ? [] : budgetViolations)
];

const list = (values) => values.length ? values.map((value) => `- ${value}`) : ['- None'];
const signed = (value) => value > 0 ? `+${value}` : String(value);
const differentialInventories = [
  {
    name: 'Registered Authority Location Count',
    delta: registeredAuthorityLocationDelta,
    classification: 'registered rule'
  },
  {
    name: 'Registered canonical-contract results',
    delta: registeredCanonicalContractDelta,
    classification: 'registered rule'
  },
  {
    name: 'Privileged operator permission entries',
    delta: privilegedOperatorPermissionDelta,
    classification: 'hard authority'
  },
  {
    name: 'Privileged importer permission entries',
    delta: privilegedImporterPermissionDelta,
    classification: 'hard authority'
  },
  {
    name: 'Production files',
    delta: productionFileDelta,
    classification: 'budgeted'
  },
  {
    name: 'Production modules',
    delta: sourceInventoryDeltas.modules,
    classification: 'report-only'
  },
  {
    name: 'Exports',
    delta: sourceInventoryDeltas.exports,
    classification: 'report-only'
  },
  {
    name: 'create* declarations',
    delta: sourceInventoryDeltas.createDeclarations,
    classification: 'report-only'
  },
  {
    name: 'Reactive declarations',
    delta: sourceInventoryDeltas.reactiveDeclarations,
    classification: 'report-only'
  },
  {
    name: 'Watcher calls',
    delta: sourceInventoryDeltas.watcherCalls,
    classification: 'report-only'
  },
  {
    name: 'Compatibility-like names',
    delta: sourceInventoryDeltas.compatibilityNames,
    classification: 'report-only'
  },
  {
    name: 'Resource-like names',
    delta: sourceInventoryDeltas.resourceNames,
    classification: 'report-only'
  },
  {
    name: 'Session object keys',
    delta: sourceInventoryDeltas.sessionObjectKeys,
    classification: 'report-only'
  },
  {
    name: 'Privileged import-target fan-out',
    delta: privilegedImportFanOutDelta,
    classification: 'report-only'
  },
  {
    name: 'First-party static import cycles',
    delta: cycleDelta,
    classification: 'hard'
  }
];
const differentialSummaryRows = differentialInventories.map(({ name, delta, classification }) => (
  `| ${name} | ${delta.before.length} | ${delta.added.length} | ${delta.removed.length} `
  + `| ${delta.after.length} | ${signed(delta.delta)} | ${classification} |`
));
const differentialChangeDetails = differentialInventories.flatMap(({ name, delta }) => {
  if (!delta.added.length && !delta.removed.length) return [];
  return [
    `### ${name}`,
    '',
    `Added (${delta.added.length}):`,
    '',
    ...list(delta.added),
    '',
    `Removed (${delta.removed.length}):`,
    '',
    ...list(delta.removed),
    ''
  ];
});
const fanOutPaths = [...new Set([
  ...basePrivilegedFanOutCounts.keys(),
  ...headPrivilegedFanOutCounts.keys()
])].sort();
const fanOutRows = fanOutPaths.map((path) => {
  const beforeCount = basePrivilegedFanOutCounts.get(path) || 0;
  const afterCount = headPrivilegedFanOutCounts.get(path) || 0;
  return `| ${path} | ${beforeCount} | ${afterCount} | ${signed(afterCount - beforeCount)} |`;
});
const statusRows = productionPaths.map((path) => `${changed.get(path)} ${path}`);
const report = [
  '# Web change budget',
  '',
  `- Base: \`${base}\``,
  `- Head: \`${head || 'working tree'}\``,
  `- Privileged allowlist revision: \`${policyRevision}\``,
  `- Architecture rule base: ${baseArchitectureRulesSource === null ? 'absent' : `\`${base}\``}`,
  `- Architecture rule candidate: ${proposedArchitectureRulesSource === null ? 'absent' : `\`${head || 'working tree'}\``}`,
  '- Policy guide: `docs/internal/WEB_CHANGE_POLICY.md`',
  `- Change context: ${changeContext.context}`,
  `- Promotion source ancestry: ${promotionSourceAncestry.status}`,
  `- Size review: ${budgetViolations.length ? 'REQUIRED' : 'NOT REQUIRED'}`,
  `- Promotion aggregation observations: ${promotionAggregationObservations.length}`,
  `- Selected profile: ${selectedProfileName}`,
  `- Production file limit: ${selectedProfile.productionFiles}`,
  `- Gross churn limit: ${selectedProfile.grossChurn}`,
  `- Net-addition limit: ${selectedProfile.netAdditions}`,
  `- \`architecture-change\` label: ${architectureChange ? 'present' : 'absent'}`,
  `- Result: **${enforcedViolations.length ? 'FAIL' : 'PASS'}**`,
  '',
  '## Architecture differential summary',
  '',
  '| Inventory | Before | Added | Removed | After | Delta | Classification |',
  '| --- | ---: | ---: | ---: | ---: | ---: | --- |',
  ...differentialSummaryRows,
  '',
  'Production line change:',
  '',
  `- Additions: ${productionTotals.additions}`,
  `- Deletions: ${productionTotals.deletions}`,
  `- Gross churn: ${productionGrossChurn}`,
  `- Net additions: ${productionNetAdditions}`,
  '',
  '## Differential inventory changes',
  '',
  ...(differentialChangeDetails.length ? differentialChangeDetails : ['- None', '']),
  '## First-party static import graph',
  '',
  `- Before: ${baseCycleResult.nodeCount} modules, ${baseCycleResult.edgeCount} edges, ${baseCycleResult.cycles.length} cycles`,
  `- After: ${headCycleResult.nodeCount} modules, ${headCycleResult.edgeCount} edges, ${headCycleResult.cycles.length} cycles`,
  `- Delta: ${signed(headCycleResult.cycles.length - baseCycleResult.cycles.length)} cycles`,
  '- Cycle definition: SCC size > 1, or a one-node SCC with a self-edge.',
  '- Head cycles (sorted nodes and internal edges):',
  ...list(headCycleResult.cycles.map(({ subject }) => subject)),
  '- Base graph resolution errors:',
  ...list(baseImportGraph.errors),
  '- Head graph resolution errors:',
  ...list(headImportGraph.errors),
  '',
  '## Report-only privileged import-target fan-out',
  '',
  '| Production module | Before | After | Delta |',
  '| --- | ---: | ---: | ---: |',
  ...(fanOutRows.length ? fanOutRows : ['| None | 0 | 0 | 0 |']),
  '',
  '## Production files touched',
  '',
  ...list(statusRows),
  '',
  '## Production additions/deletions',
  '',
  `- Additions: ${productionTotals.additions}`,
  `- Deletions: ${productionTotals.deletions}`,
  `- Gross churn: ${productionGrossChurn}`,
  `- Net additions: ${productionNetAdditions}`,
  `- Binary production files: ${productionTotals.binary}`,
  '',
  '## Report-only new production modules',
  '',
  ...list(newModules),
  '',
  '## Report-only new exports and create* owners',
  '',
  ...list([...newExports, ...newCreateOwners]),
  '',
  '## Report-only new reactive state/watchers',
  '',
  ...list([...newReactiveState, ...newWatchers]),
  '',
  '## Unapproved privileged capability owners/importers',
  '',
  ...list(unapprovedCapabilities),
  '',
  '## Removed privileged allowlist entries',
  '',
  ...list(policyContractions),
  '',
  '## Added privileged allowlist entries',
  '',
  ...list(policyExpansions),
  '',
  '## Active privileged owners/importers excluded by proposed policy',
  '',
  ...list(proposedPolicyExclusions),
  '',
  '## Missing base privileged allowlist keys',
  '',
  ...list(missingPolicyKeys),
  '',
  '## Added privileged allowlist keys',
  '',
  ...list(addedPolicyKeys),
  '',
  '## Candidate architecture authority delta',
  '',
  `- Classification: ${architectureAuthorityDelta.classification}`,
  ...list(architectureAuthorityDelta.changes.map(({ rule, field, direction }) => (
    `${direction} ${rule} (${field})`
  ))),
  '',
  '## Candidate architecture admission observations',
  '',
  ...list(candidateArchitectureObservations),
  '',
  '## Candidate architecture rule errors',
  '',
  ...list(candidateArchitectureErrors),
  '',
  '## Active architecture rule results',
  '',
  ...list(activeArchitectureRows),
  '',
  '## Active architecture rule errors',
  '',
  ...list(activeArchitectureErrors),
  '',
  '## Report-only cache/token/handle/journal/protocol/manager names',
  '',
  ...list(newNamedResources),
  '',
  '## Report-only session object keys and compatibility names',
  '',
  ...list([...newSessionFields, ...newCompatibilityPaths]),
  '',
  '## Added binary runtime files',
  '',
  ...list(addedBinaryRuntimePaths),
  '',
  '## Changed vendored runtime files',
  '',
  ...list(changedVendorPaths),
  '',
  '## Dependency changes',
  '',
  ...list(dependencyChanges),
  '',
  '## Guard files touched',
  '',
  ...list(changedGuards),
  '',
  '## Size review reasons',
  '',
  ...list(budgetViolations),
  '',
  '## Promotion aggregation observations',
  '',
  ...list(promotionAggregationObservations),
  '',
  '## Violations',
  '',
  ...list(enforcedViolations),
  ''
].join('\n');

process.stdout.write(report);
if (process.env.GITHUB_STEP_SUMMARY) {
  appendFileSync(resolve(process.env.GITHUB_STEP_SUMMARY), report, 'utf8');
}
if (enforcedViolations.length) process.exitCode = 1;
