#!/usr/bin/env node

import { appendFileSync, existsSync, readFileSync } from 'node:fs';
import { join, resolve } from 'node:path';
import { execFileSync } from 'node:child_process';
import {
  detectPrivilegedWebCapabilities,
  detectReportOnlySourceFacts,
  isWebSessionSourcePath,
  WEB_PRIVILEGED_CAPABILITY_KEYS,
  WEB_PRIVILEGED_IMPORT_TARGETS
} from './web-architecture-detectors.mjs';

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
  'docs/internal/ARCHITECTURE_FITNESS_FUNCTION_RATCHET.md',
  '.github/pull_request_template.md',
  'tools/check-web-change-budget.mjs',
  'tools/web-architecture-detectors.mjs',
  'tools/web-architecture-evaluation.mjs',
  'tools/web-architecture-rules.json',
  'tools/web-architecture-violations.json',
  'tools/web-change-source.mjs',
  policyPath,
  'docs/internal/WEB_CHANGE_POLICY.md',
  'tests/web/architecture-contracts.test.mjs',
  'tests/web/architecture-ratchet-fixtures.test.mjs',
  '.github/workflows/test.yml',
  '.github/workflows/web-base-policy.yml'
]);
const changedGuards = [...changed.keys()].filter((path) => guardPaths.has(path)).sort();
const checkerImplementationPaths = new Set([
  'tools/check-web-change-budget.mjs',
  'tools/web-architecture-detectors.mjs',
  'tools/web-architecture-evaluation.mjs',
  'tools/web-change-source.mjs'
]);
const authorityPaths = new Set([
  'docs/internal/ARCHITECTURE_FITNESS_FUNCTION_RATCHET.md',
  'tools/web-architecture-rules.json',
  'tools/web-architecture-violations.json',
  policyPath,
  '.github/workflows/test.yml',
  '.github/workflows/web-base-policy.yml'
]);
const changedCheckerImplementations = [...changed.keys()]
  .filter((path) => checkerImplementationPaths.has(path));
const changedAuthorities = [...changed.keys()].filter((path) => authorityPaths.has(path));

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
  const beforeFacts = detectReportOnlySourceFacts(before);
  const afterFacts = detectReportOnlySourceFacts(after);
  newExports.push(...addedEntries(
    new Set(beforeFacts.exportedNames), new Set(afterFacts.exportedNames), `${path}: `
  ));
  newCreateOwners.push(...addedEntries(
    new Set(beforeFacts.declaredNames.filter((name) => /^create[A-Z]/.test(name))),
    new Set(afterFacts.declaredNames.filter((name) => /^create[A-Z]/.test(name))),
    `${path}: `
  ));
  newReactiveState.push(...addedEntries(
    new Set(beforeFacts.reactiveDeclarations),
    new Set(afterFacts.reactiveDeclarations),
    `${path}: `
  ));
  const watcherIncrease = afterFacts.watcherCount - beforeFacts.watcherCount;
  if (watcherIncrease > 0) newWatchers.push(`${path}: +${watcherIncrease} watcher call(s)`);
  newNamedResources.push(...addedEntries(
    new Set(beforeFacts.namedResourceNames),
    new Set(afterFacts.namedResourceNames),
    `${path}: `
  ));
  newCompatibilityPaths.push(...addedEntries(
    new Set(beforeFacts.compatibilityNames),
    new Set(afterFacts.compatibilityNames),
    `${path}: `
  ));
  if (isWebSessionSourcePath(path)) {
    newSessionFields.push(...addedEntries(
      new Set(beforeFacts.objectKeys), new Set(afterFacts.objectKeys), `${path}: `
    ));
  }
  const bareBefore = new Set(
    beforeFacts.importSpecifiers.filter((specifier) => !specifier.startsWith('.'))
  );
  const bareAfter = new Set(
    afterFacts.importSpecifiers.filter((specifier) => !specifier.startsWith('.'))
  );
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
const allProductionJavaScriptPaths = listProductionJavaScriptPaths();
const allProductionSources = new Map(allProductionJavaScriptPaths.map((path) => [
  path.slice('gbdraw/web/js/'.length),
  readHeadFile(path) || ''
]));

const detectedCapabilities = detectPrivilegedWebCapabilities(allProductionSources);
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
