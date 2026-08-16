#!/usr/bin/env node

import { appendFileSync, existsSync, readFileSync } from 'node:fs';
import { join, posix, resolve } from 'node:path';
import { execFileSync } from 'node:child_process';
import { literalImportSpecifiers, maskJavaScript } from './web-change-source.mjs';

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

const newModules = productionJavaScriptPaths.filter((path) => changed.get(path) === 'A');
const addedBinaryRuntimePaths = productionPaths.filter((path) => {
  const counts = numstat.get(path);
  return changed.get(path) === 'A' && counts?.additions === '-' && counts?.deletions === '-';
});
const changedVendorPaths = productionPaths.filter((path) => path.startsWith('gbdraw/web/vendor/'));
const policyPath = 'tools/web-change-policy.json';
const guardPaths = new Set([
  'tools/check-web-change-budget.mjs',
  'tools/web-change-source.mjs',
  policyPath,
  'docs/internal/WEB_CHANGE_POLICY.md',
  'tests/web/architecture-contracts.test.mjs',
  '.github/workflows/test.yml',
  '.github/workflows/web-base-policy.yml'
]);
const changedGuards = [...changed.keys()].filter((path) => guardPaths.has(path)).sort();
const checkerImplementationPaths = new Set([
  'tools/check-web-change-budget.mjs',
  'tools/web-change-source.mjs'
]);
const authorityPaths = new Set([
  policyPath,
  '.github/workflows/test.yml',
  '.github/workflows/web-base-policy.yml'
]);
const changedCheckerImplementations = [...changed.keys()]
  .filter((path) => checkerImplementationPaths.has(path));
const changedAuthorities = [...changed.keys()].filter((path) => authorityPaths.has(path));

const exportedNames = (source = '') => {
  const code = maskJavaScript(source);
  const names = new Set();
  const declaration = /^\s*export\s+(?:async\s+)?(?:const|let|var|function|class)\s+([A-Za-z_$][\w$]*)/gm;
  for (const match of code.matchAll(declaration)) names.add(match[1]);
  for (const match of code.matchAll(/^\s*export\s*\{([\s\S]*?)\}/gm)) {
    match[1].split(',').forEach((entry) => {
      const cleaned = entry.replace(/\/\*[\s\S]*?\*\//g, '').trim();
      if (!cleaned) return;
      const parts = cleaned.split(/\s+as\s+/);
      names.add((parts[1] || parts[0]).trim());
    });
  }
  if (/^\s*export\s+default\b/m.test(code)) names.add('default');
  for (const match of code.matchAll(/^\s*export\s+\*\s+as\s+([A-Za-z_$][\w$]*)/gm)) {
    names.add(match[1]);
  }
  return names;
};

const declaredNames = (source = '') => new Set(
  [...maskJavaScript(source).matchAll(/\b(?:const|let|var|function|class)\s+([A-Za-z_$][\w$]*)/g)]
    .map((match) => match[1])
);

const reactiveDeclarations = (source = '') => new Set(
  [...maskJavaScript(source).matchAll(/\b(?:const|let|var)\s+([A-Za-z_$][\w$]*)\s*=\s*(?:[A-Za-z_$][\w$]*\.)?(ref|shallowRef|reactive|computed)\s*\(/g)]
    .map((match) => `${match[1]} (${match[2]})`)
);

const countWatchers = (source = '') => [
  ...maskJavaScript(source).matchAll(/\bwatch(?:Effect)?\s*\(/g)
].length;
const importSpecifiers = literalImportSpecifiers;

const resolvedProductionImport = (owner, specifier) => {
  if (!specifier.startsWith('.')) return null;
  const relativeOwner = owner.slice('gbdraw/web/js/'.length);
  return posix.normalize(posix.join(posix.dirname(relativeOwner), specifier));
};

const capabilitySpecs = [
  {
    name: 'Render request',
    imports: ['services/session-request.js'],
    owner: /\bbuildCanonicalRenderRequest\s*\(/
  },
  {
    name: 'Diagram Worker',
    imports: [
      'services/diagram-generation.js',
      'services/diagram-worker-protocol.js',
      'workers/diagram-generation-worker.js'
    ],
    owner: /\b(?:new\s+Worker|runDiagramGeneration|loadPyodide)\s*\(/
  },
  {
    name: 'Python helper',
    imports: ['app/python-helpers.js'],
    owner: /\b(?:runPythonAsync|runPython|PYTHON_HELPERS|globals\.get)\b/
  },
  {
    name: 'Resource staging',
    imports: ['services/diagram-resource-staging.js', 'services/resource-payload-owner.js'],
    owner: /\b(?:createDiagramResourceTransport|stageRenderResources|setResourcePayloadOwner)\b/
  },
  {
    name: 'SVG/Result admission',
    imports: ['services/svg-sanitization.js', 'services/svg-result-ingestion.js'],
    owner: /\b(?:sanitizeSvgContent|ingestSvgResult|markCommitted)\b/
  },
  {
    name: 'Mounted SVG/Result replacement',
    imports: ['services/svg-serialization.js', 'app/preview-runtime.js'],
    owner: /\b(?:serializeCleanSvg|flushActiveResult)\s*\(|\b(?:results|state\.results)\.value(?:\[[^\]]+\])?\s*=(?!=)/
  },
  {
    name: 'History',
    imports: ['services/history.js', 'services/history-snapshot.js'],
    owner: /\b(?:createHistoryManager|beginCheckpoint|commitCheckpoint|buildArtifactCheckpoint)\b/
  },
  {
    name: 'Session',
    imports: [
      'services/config.js',
      'services/gallery-session-migration.js',
      'services/session-authority.js',
      'services/session-file.js'
    ],
    owner: /\b(?:exportSession|compressSessionData|migrateSessionDataToCurrent|migratePersistedGalleryConfig)\b/
  },
  {
    name: 'Canonical editor state',
    imports: [
      'app/feature-editor.js',
      'app/legend.js',
      'app/right-drawer.js',
      'app/preview-runtime.js'
    ],
    owner: /\b(?:createFeatureEditor|createLegendManager|createRightDrawerController|createPreviewRuntime)\b/
  }
];

const sessionOwnerPaths = new Set([
  'gbdraw/web/js/services/config.js',
  'gbdraw/web/js/services/gallery-session-migration.js',
  'gbdraw/web/js/services/session-authority.js',
  'gbdraw/web/js/services/session-file.js',
  'gbdraw/web/js/services/session-request.js'
]);
const objectKeys = (source = '') => new Set(
  [...maskJavaScript(source).matchAll(/(?:^|[,{]\s*)([A-Za-z_$][\w$]*)\s*:/gm)]
    .map((match) => match[1])
);
const compatibilityNames = (source = '') => new Set(
  [...maskJavaScript(source).matchAll(/\b(?:[A-Za-z_$][\w$]*(?:Migration|Migrator|Legacy|Fallback)[\w$]*|(?:migrate|promoteLegacy|normalizeLegacy|readLegacy)[A-Za-z0-9_$]*)\b/g)]
    .map((match) => match[0])
);
const namedResourceNames = (source = '') => new Set(
  [...maskJavaScript(source).matchAll(/\b[A-Za-z_$][\w$]*\b/g)]
    .map((match) => match[0])
    .filter((name) => /(?:cache|token|handle|journal|protocol|manager)/i.test(name))
);

const addedEntries = (before, after, prefix = '') => [...after]
  .filter((entry) => !before.has(entry))
  .map((entry) => `${prefix}${entry}`);

const newExports = [];
const newCreateOwners = [];
const newReactiveState = [];
const newWatchers = [];
const newNamedResources = [];
const newCompatibilityPaths = [];
const newSessionFields = [];
const newBareImports = [];

productionJavaScriptPaths.forEach((path) => {
  const before = readRevisionFile(base, path) || '';
  const after = readHeadFile(path) || '';
  newExports.push(...addedEntries(exportedNames(before), exportedNames(after), `${path}: `));
  newCreateOwners.push(...addedEntries(
    new Set([...declaredNames(before)].filter((name) => /^create[A-Z]/.test(name))),
    new Set([...declaredNames(after)].filter((name) => /^create[A-Z]/.test(name))),
    `${path}: `
  ));
  newReactiveState.push(...addedEntries(
    reactiveDeclarations(before), reactiveDeclarations(after), `${path}: `
  ));
  const watcherIncrease = countWatchers(after) - countWatchers(before);
  if (watcherIncrease > 0) newWatchers.push(`${path}: +${watcherIncrease} watcher call(s)`);
  newNamedResources.push(...addedEntries(
    namedResourceNames(before),
    namedResourceNames(after),
    `${path}: `
  ));
  newCompatibilityPaths.push(...addedEntries(
    compatibilityNames(before), compatibilityNames(after), `${path}: `
  ));
  if (sessionOwnerPaths.has(path)) {
    newSessionFields.push(...addedEntries(objectKeys(before), objectKeys(after), `${path}: `));
  }
  const bareBefore = new Set(importSpecifiers(before).filter((specifier) => !specifier.startsWith('.')));
  const bareAfter = new Set(importSpecifiers(after).filter((specifier) => !specifier.startsWith('.')));
  newBareImports.push(...addedEntries(bareBefore, bareAfter, `${path}: `));
});

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
const missingPolicyKeys = [];
if (basePolicySource !== null && changed.has(policyPath)) {
  proposedPolicy = parsePolicy(readHeadFile(policyPath), head || 'working tree');
  for (const section of ['allowedPrivilegedImporters', 'allowedPrivilegedOwners']) {
    Object.entries(policy[section]).forEach(([name, allowedPaths]) => {
      if (!Object.hasOwn(proposedPolicy[section], name)) {
        missingPolicyKeys.push(`${section}.${name}`);
      }
      const proposedPaths = new Set(proposedPolicy[section][name] || []);
      allowedPaths.forEach((path) => {
        if (!proposedPaths.has(path)) policyContractions.push(`${section}.${name}: ${path}`);
      });
    });
  }
}
const allProductionJavaScriptPaths = listProductionJavaScriptPaths();
const allProductionSources = new Map(allProductionJavaScriptPaths.map((path) => [
  path,
  readHeadFile(path) || ''
]));

const importerTargets = new Set(capabilitySpecs.flatMap((capability) => capability.imports));
const capabilityCoverageViolations = (candidatePolicy) => {
  const violations = [];
  importerTargets.forEach((target) => {
    const allowed = new Set(candidatePolicy.allowedPrivilegedImporters[target] || []);
    allProductionSources.forEach((source, path) => {
      const imports = importSpecifiers(source)
        .map((specifier) => resolvedProductionImport(path, specifier));
      if (imports.includes(target)) {
        const relativePath = path.slice('gbdraw/web/js/'.length);
        if (!allowed.has(relativePath)) {
          violations.push(`${target}: importer ${relativePath}`);
        }
      }
    });
  });

  capabilitySpecs.forEach((capability) => {
    const allowed = new Set(candidatePolicy.allowedPrivilegedOwners[capability.name] || []);
    allProductionSources.forEach((source, path) => {
      if (!capability.owner.test(maskJavaScript(source))) return;
      const relativePath = path.slice('gbdraw/web/js/'.length);
      if (!allowed.has(relativePath)) {
        violations.push(`${capability.name}: owner ${relativePath}`);
      }
    });
  });
  return violations.sort();
};

const unapprovedCapabilities = capabilityCoverageViolations(policy);
const proposedPolicyExclusions = proposedPolicy
  ? capabilityCoverageViolations(proposedPolicy)
  : [];

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
if (newProductionDependencies.length) {
  integrityViolations.push('new production dependencies are not allowed');
}
if (addedBinaryRuntimePaths.length) {
  integrityViolations.push('added binary runtime files are not allowed');
}
if (changedVendorPaths.length) {
  integrityViolations.push('changes under gbdraw/web/vendor/ are not allowed');
}
if (productionPaths.length && changedGuards.length) {
  integrityViolations.push('production runtime files and Web guard/CI files changed together');
}
if (changedCheckerImplementations.length && changedAuthorities.length) {
  integrityViolations.push(
    'Web checker/source parser and authority policy/workflow files changed together'
  );
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
const enforcedViolations = [
  ...integrityViolations,
  ...budgetViolations
];

const list = (values) => values.length ? values.map((value) => `- ${value}`) : ['- None'];
const statusRows = productionPaths.map((path) => `${changed.get(path)} ${path}`);
const report = [
  '# Web change budget',
  '',
  `- Base: \`${base}\``,
  `- Head: \`${head || 'working tree'}\``,
  `- Privileged allowlist revision: \`${policyRevision}\``,
  '- Policy guide: `docs/internal/WEB_CHANGE_POLICY.md`',
  `- Selected profile: ${selectedProfileName}`,
  `- Production file limit: ${selectedProfile.productionFiles}`,
  `- Gross churn limit: ${selectedProfile.grossChurn}`,
  `- Net-addition limit: ${selectedProfile.netAdditions}`,
  `- \`architecture-change\` label: ${architectureChange ? 'present' : 'absent'}`,
  `- Result: **${enforcedViolations.length ? 'FAIL' : 'PASS'}**`,
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
  '## Active privileged owners/importers excluded by proposed policy',
  '',
  ...list(proposedPolicyExclusions),
  '',
  '## Missing base privileged allowlist keys',
  '',
  ...list(missingPolicyKeys),
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
