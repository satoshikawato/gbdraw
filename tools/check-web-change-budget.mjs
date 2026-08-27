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
  createArchitectureSubjectDelta,
  evaluateArchitectureRuleResult,
  findDirectedGraphCycles,
  summarizeArchitectureInventory,
  validateArchitectureRuleRegistry
} from './web-architecture-evaluation.mjs';
import { parseCurrentProductImpactDecisions } from './web-product-impact-decision-source.mjs';
import {
  evaluateProductImpact,
  validateProductDecisionAuthority,
  validateProductImpactMap
} from './web-product-impact-evaluation.mjs';
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

const stableReasons = (values) => [...new Set(values.filter(Boolean))]
  .sort((left, right) => left.localeCompare(right));
const safeReportText = (value) => String(value || '')
  .replace(/([\\`*{}\[\]<>#|])/g, '\\$1');
const reportSentenceFragment = (value) => safeReportText(value).replace(/[.!?]+$/, '');
const inlineCode = (value) => `\`${String(value || '').replaceAll('`', '\\`')}\``;
const createPolicyResult = ({
  blockingViolations = [],
  reviewReasons = [],
  context,
  observations = {}
}) => {
  const stableBlockingViolations = stableReasons(blockingViolations);
  const stableReviewReasons = stableReasons(reviewReasons);
  return Object.freeze({
    gate: stableBlockingViolations.length ? 'FAIL' : 'PASS',
    review: stableReviewReasons.length ? 'REQUIRED' : 'CLEAR',
    blockingViolations: Object.freeze(stableBlockingViolations),
    reviewReasons: Object.freeze(stableReviewReasons),
    context,
    observations: Object.freeze(observations)
  });
};
const list = (values) => values.length ? values.map((value) => `- ${value}`) : ['- None'];
const reviewCategorySummary = (reviewReasons) => {
  const counts = new Map();
  reviewReasons.forEach((reason) => {
    const category = reason.includes(':') ? reason.slice(0, reason.indexOf(':')) : 'Other';
    counts.set(category, (counts.get(category) || 0) + 1);
  });
  return [...counts]
    .sort(([left], [right]) => left.localeCompare(right))
    .map(([category, count]) => `${category}=${count}`)
    .join(', ');
};
const emitPolicyReport = (result, report) => {
  process.stdout.write(report);
  if (process.env.GITHUB_STEP_SUMMARY) {
    appendFileSync(resolve(process.env.GITHUB_STEP_SUMMARY), report, 'utf8');
  }
  if (process.env.GITHUB_ACTIONS === 'true' && result.review === 'REQUIRED') {
    process.stdout.write(
      '::warning title=Web policy review required::Web policy review required: '
      + `reasons=${result.reviewReasons.length}; `
      + `categories=${reviewCategorySummary(result.reviewReasons)}; see the step summary.\n`
    );
  }
  if (result.gate === 'FAIL') process.exitCode = 1;
};

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
const sizeReviewProfiles = Object.freeze({
  ordinary: Object.freeze({ productionFiles: 8, grossChurn: 800, netAdditions: 100 }),
  architecture: Object.freeze({ productionFiles: 12, grossChurn: 1500, netAdditions: 400 })
});
const selectedProfileName = architectureChange ? 'architecture' : 'ordinary';
const selectedProfile = sizeReviewProfiles[selectedProfileName];
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
const emptyPullRequestProductImpactMetadata = Object.freeze({
  authorLogin: '',
  body: '',
  headSha: '',
  trustedToDev: false
});
const pullRequestProductImpactMetadata = (() => {
  if (
    githubEventName !== 'pull_request_target'
    || changeContext.isPromotion
    || changeContext.errors.length
    || !githubEventPayloadSource
  ) {
    return emptyPullRequestProductImpactMetadata;
  }
  try {
    const payload = JSON.parse(githubEventPayloadSource);
    const pullRequest = payload?.pull_request;
    if (pullRequest?.base?.ref !== 'dev') return emptyPullRequestProductImpactMetadata;
    return Object.freeze({
      authorLogin: typeof pullRequest?.user?.login === 'string'
        ? pullRequest.user.login
        : '',
      body: typeof pullRequest?.body === 'string' ? pullRequest.body : '',
      headSha: typeof pullRequest?.head?.sha === 'string' ? pullRequest.head.sha : '',
      trustedToDev: true
    });
  } catch (_error) {
    return emptyPullRequestProductImpactMetadata;
  }
})();

const promotionSourceCoverage = (() => {
  if (!changeContext.isPromotion) {
    return { status: 'NOT_APPLICABLE', basis: 'NOT_APPLICABLE', violation: '' };
  }
  if (changeContext.errors.length) {
    return { status: 'NOT_EVALUATED', basis: 'INVALID_CONTEXT', violation: '' };
  }
  const probeGit = (args) => spawnSync(
    'git',
    ['-C', repositoryRoot, ...args],
    { encoding: 'utf8', stdio: ['ignore', 'pipe', 'pipe'] }
  );
  const error = () => ({
    status: 'ERROR',
    basis: 'UNVERIFIED',
    violation: (
      'Promotion source coverage could not be verified. '
      + 'Fetch complete base and head history, then rerun the promotion.'
    )
  });
  const failure = () => ({
    status: 'FAIL',
    basis: 'MAIN_CONTENT_MISSING',
    violation: (
      'The promotion source does not contain the current main content. '
      + 'Merge or rebase main into dev, then rerun the promotion.'
    )
  });

  const directAncestry = probeGit(['merge-base', '--is-ancestor', base, head]);
  if (directAncestry.status === 0) {
    return { status: 'PASS', basis: 'DIRECT_ANCESTRY', violation: '' };
  }
  if (directAncestry.status !== 1) return error();

  const parentResult = probeGit(['rev-list', '--parents', '-n', '1', base]);
  if (parentResult.status !== 0) return error();
  const [baseCommit, ...baseParents] = parentResult.stdout.trim().split(/\s+/);
  if (baseCommit !== base || baseParents.length < 2) return failure();

  const baseTreeResult = probeGit(['rev-parse', `${base}^{tree}`]);
  if (baseTreeResult.status !== 0) return error();
  const baseTree = baseTreeResult.stdout.trim();
  for (const parent of baseParents) {
    const parentAncestry = probeGit(['merge-base', '--is-ancestor', parent, head]);
    if (parentAncestry.status !== 0 && parentAncestry.status !== 1) return error();
    if (parentAncestry.status === 1) continue;
    const parentTreeResult = probeGit(['rev-parse', `${parent}^{tree}`]);
    if (parentTreeResult.status !== 0) return error();
    if (parentTreeResult.stdout.trim() === baseTree) {
      return { status: 'PASS', basis: 'MERGE_PARENT_TREE', violation: '' };
    }
  }
  return failure();
})();

const renderBoundedContextReport = () => {
  const result = createPolicyResult({
    blockingViolations: [
      ...changeContext.errors,
      promotionSourceCoverage.violation
    ],
    reviewReasons: [],
    context: changeContext.context,
    observations: { promotionSourceCoverage }
  });
  const report = [
    '# Web change policy',
    '',
    `- Context: ${result.context}`,
    `- Base: \`${base}\``,
    `- Head: \`${head || 'working tree'}\``,
    `- Gate: **${result.gate}**`,
    `- Review: **${result.review}**`,
    '',
    '## Blocking violations',
    '',
    ...list(result.blockingViolations),
    '',
    '## Review reasons',
    '',
    ...list(result.reviewReasons),
    ...(result.context === 'PROMOTION' ? [
      '',
      '## Promotion source coverage',
      '',
      `- Status: ${promotionSourceCoverage.status}`,
      `- Basis: ${promotionSourceCoverage.basis}`,
      '- Scope: topology and source-content coverage only.',
      '- Exact-SHA workflow evidence and result-tree identity remain owned by the promotion readiness helper and trusted workflow.',
      ''
    ] : [''])
  ].join('\n');
  emitPolicyReport(result, report);
};

const runOrdinaryPolicy = () => {
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
const changedPaths = [...changed.keys()].sort();

const numstat = new Map(
  parseDiffLines(runGit(repositoryRoot, [
    'diff', '--no-renames', '--numstat', ...diffRefs, '--'
  ])).map(([additions, deletions, path]) => [path, { additions, deletions }])
);

const revisionFileCache = new Map();
const readRevisionFile = (revision, path) => {
  const key = `${revision}\u0000${path}`;
  if (revisionFileCache.has(key)) return revisionFileCache.get(key);
  let source;
  try {
    source = runGit(
      repositoryRoot,
      ['show', `${revision}:${path}`],
      { stdio: ['ignore', 'pipe', 'ignore'] }
    );
  } catch (_error) {
    source = null;
  }
  revisionFileCache.set(key, source);
  return source;
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
const acceptedViolationsPath = 'tools/web-architecture-violations.json';
const productImpactPolicyPath = 'docs/internal/PRODUCT_IMPACT_RATCHET.md';
const productImpactMapPath = 'tools/web-product-impact-map.json';
const productDecisionAuthorityPath = 'tools/web-product-decisions.json';
const productImpactEvaluationPath = 'tools/web-product-impact-evaluation.mjs';
const productImpactDecisionSourcePath = 'tools/web-product-impact-decision-source.mjs';
const productImpactFixturePath = 'tests/web/product-impact-ratchet-fixtures.test.mjs';
const trustedWorkflowPath = '.github/workflows/web-base-policy.yml';
const guardPaths = new Set([
  'docs/internal/ARCHITECTURE_FITNESS_FUNCTION_RATCHET.md',
  productImpactPolicyPath,
  '.github/pull_request_template.md',
  'tools/check-web-change-budget.mjs',
  'tools/web-architecture-detectors.mjs',
  'tools/web-architecture-evaluation.mjs',
  productImpactEvaluationPath,
  productImpactDecisionSourcePath,
  architectureRulesPath,
  acceptedViolationsPath,
  productImpactMapPath,
  productDecisionAuthorityPath,
  'tools/web-change-source.mjs',
  'tools/web-promotion-context.mjs',
  'tools/check-promotion-readiness.mjs',
  policyPath,
  'docs/internal/WEB_CHANGE_POLICY.md',
  'tests/web/architecture-contracts.test.mjs',
  'tests/web/architecture-ratchet-fixtures.test.mjs',
  productImpactFixturePath,
  'tests/web/promotion-readiness.test.mjs',
  'tests/web/web-promotion-context.test.mjs',
  '.github/workflows/gallery-publication.yml',
  '.github/workflows/deploy_web.yml',
  '.github/workflows/test.yml',
  '.github/workflows/web-base-policy.yml'
]);
const changedGuards = [...changed.keys()].filter((path) => guardPaths.has(path)).sort();
const checkerImplementationPaths = new Set([
  'tools/check-web-change-budget.mjs',
  'tools/web-architecture-detectors.mjs',
  'tools/web-architecture-evaluation.mjs',
  productImpactEvaluationPath,
  productImpactDecisionSourcePath,
  'tools/web-change-source.mjs',
  'tools/web-promotion-context.mjs',
  'tools/check-promotion-readiness.mjs'
]);
const authorityPaths = new Set([
  'docs/internal/ARCHITECTURE_FITNESS_FUNCTION_RATCHET.md',
  productImpactPolicyPath,
  'docs/internal/WEB_CHANGE_POLICY.md',
  architectureRulesPath,
  acceptedViolationsPath,
  productImpactMapPath,
  productDecisionAuthorityPath,
  policyPath,
  '.github/workflows/gallery-publication.yml',
  '.github/workflows/deploy_web.yml',
  '.github/workflows/test.yml',
  '.github/workflows/web-base-policy.yml'
]);
const changedCheckerImplementations = [...changed.keys()]
  .filter((path) => checkerImplementationPaths.has(path))
  .sort();
const changedAuthorities = [...changed.keys()]
  .filter((path) => authorityPaths.has(path))
  .sort();
const governancePaths = new Set([
  ...authorityPaths,
  '.github/pull_request_template.md'
]);
const changedGovernancePaths = [...changed.keys()]
  .filter((path) => governancePaths.has(path))
  .sort();
const narrowAuthorityBundlePaths = new Set([
  architectureRulesPath,
  productImpactMapPath,
  productDecisionAuthorityPath
]);
const baseAcceptedViolationsSource = readRevisionFile(base, acceptedViolationsPath);
const candidateAcceptedViolationsSource = readHeadFile(acceptedViolationsPath);
const acceptedViolationAuthorityErrors = [];
if (
  changed.has(acceptedViolationsPath)
  && baseAcceptedViolationsSource === null
  && candidateAcceptedViolationsSource !== null
) {
  acceptedViolationAuthorityErrors.push(
    'accepted-violation authority cannot be introduced while frozen-rule mechanics are unavailable'
  );
}

const unsafeTrustedWorkflowChanges = [];
if (changed.has(trustedWorkflowPath)) {
  const candidateTrustedWorkflowSource = readHeadFile(trustedWorkflowPath) || '';
  const unsafePatterns = [
    /ref:\s*\$\{\{\s*github\.event\.pull_request\.head\.(?:sha|ref)\s*\}\}/,
    /\bgit\s+(?:checkout|switch)\b[^\n]*(?:HEAD_SHA|pull_request\.head)/
  ];
  if (unsafePatterns.some((pattern) => pattern.test(candidateTrustedWorkflowSource))) {
    unsafeTrustedWorkflowChanges.push(
      'candidate trusted workflow would check out or switch to pull-request head code'
    );
  }
}

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
let baseArchitectureRegistry = null;
let candidateArchitectureRegistry = null;
let baseArchitectureRuleFacts = [];
let headArchitectureRuleFacts = [];

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
  const ruleFacts = [];
  [...registry.rules]
    .sort((left, right) => left.key.localeCompare(right.key))
    .forEach((rule) => {
      const detector = WEB_ARCHITECTURE_DETECTORS[rule.detector];
      try {
        const detectorOutput = detector.detect(sources);
        const observation = classifyArchitectureRuleObservation(rule, detectorOutput);
        ruleFacts.push(Object.freeze({
          ruleKey: rule.key,
          subjectCategory: detector.subjectCategory,
          subjects: Object.freeze([...observation.subjects].sort())
        }));
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
  return {
    errors,
    failures,
    rows,
    authorityLocations,
    canonicalContracts,
    ruleFacts
  };
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
  if (baseValidation?.valid) baseArchitectureRegistry = parsedBase.registry;
  if (candidateValidation?.valid) candidateArchitectureRegistry = parsedCandidate.registry;

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
    baseArchitectureRuleFacts = beforeArchitecture.ruleFacts;
    headArchitectureRuleFacts = afterArchitecture.ruleFacts;
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
        .filter((path) => !narrowAuthorityBundlePaths.has(path))
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

const productImpactBaseErrors = [];
const productImpactCandidateErrors = [];
const productImpactDecisionDeclarationIssues = [];
const changedProductImpactAuthorityPaths = [
  productImpactMapPath,
  productDecisionAuthorityPath
].filter((path) => changed.has(path));
const candidateAuthorityCompanionPaths = changedProductImpactAuthorityPaths.length
  ? changedPaths.filter((path) => !narrowAuthorityBundlePaths.has(path))
  : [];
if (candidateAuthorityCompanionPaths.length) {
  productImpactCandidateErrors.push(
    'Product Impact authority changes must use only the narrow inert authority bundle: '
    + candidateAuthorityCompanionPaths.join(', ')
  );
}

const parseProductImpactJson = (source, path, revision, errors) => {
  if (source === null) {
    errors.push(`Missing required ${path} at ${revision}`);
    return null;
  }
  try {
    return JSON.parse(source);
  } catch (error) {
    errors.push(`Cannot parse ${path} at ${revision}: ${error.message}`);
    return null;
  }
};
const formatProductValidationErrors = (label, validation) => validation.errors.map((error) => (
  `${label} ${error.path} [${error.code}]: ${error.message}`
));
const mappedContractPaths = (map) => stableReasons((map?.concerns || []).flatMap((concern) => (
  (concern.contracts || []).map(({ ref }) => String(ref || '').split('::', 1)[0])
)));
const requireBaseContractFiles = (map, label, errors) => {
  mappedContractPaths(map).forEach((path) => {
    if (readRevisionFile(base, path) === null) {
      errors.push(`${label} references a mapped contract file absent from trusted base: ${path}`);
    }
  });
};

const baseProductImpactMapSource = readRevisionFile(base, productImpactMapPath);
const candidateProductImpactMapSource = readHeadFile(productImpactMapPath);
const baseProductDecisionAuthoritySource = readRevisionFile(base, productDecisionAuthorityPath);
const candidateProductDecisionAuthoritySource = readHeadFile(productDecisionAuthorityPath);
const baseProductImpactMap = parseProductImpactJson(
  baseProductImpactMapSource,
  productImpactMapPath,
  base,
  productImpactBaseErrors
);
const candidateProductImpactMap = parseProductImpactJson(
  candidateProductImpactMapSource,
  productImpactMapPath,
  head || 'working tree',
  productImpactCandidateErrors
);
const baseProductDecisionAuthority = parseProductImpactJson(
  baseProductDecisionAuthoritySource,
  productDecisionAuthorityPath,
  base,
  productImpactBaseErrors
);
const candidateProductDecisionAuthority = parseProductImpactJson(
  candidateProductDecisionAuthoritySource,
  productDecisionAuthorityPath,
  head || 'working tree',
  productImpactCandidateErrors
);

let baseProductImpactMapValid = false;
let candidateProductImpactMapValid = false;
let baseProductDecisionAuthorityValid = false;
let candidateProductDecisionAuthorityValid = false;
if (baseProductImpactMap && baseArchitectureRegistry) {
  const validation = validateProductImpactMap(
    baseProductImpactMap,
    baseArchitectureRegistry,
    WEB_ARCHITECTURE_DETECTORS
  );
  baseProductImpactMapValid = validation.valid;
  productImpactBaseErrors.push(...formatProductValidationErrors('base Product Impact map', validation));
  if (validation.valid) requireBaseContractFiles(
    baseProductImpactMap,
    'Base Product Impact map',
    productImpactBaseErrors
  );
}
if (candidateProductImpactMap && candidateArchitectureRegistry) {
  const validation = validateProductImpactMap(
    candidateProductImpactMap,
    candidateArchitectureRegistry,
    WEB_ARCHITECTURE_DETECTORS
  );
  candidateProductImpactMapValid = validation.valid;
  productImpactCandidateErrors.push(
    ...formatProductValidationErrors('candidate Product Impact map', validation)
  );
  if (validation.valid) requireBaseContractFiles(
    candidateProductImpactMap,
    'Candidate Product Impact map',
    productImpactCandidateErrors
  );
}
if (baseProductDecisionAuthority && baseProductImpactMapValid) {
  const validation = validateProductDecisionAuthority(
    baseProductDecisionAuthority,
    baseProductImpactMap
  );
  baseProductDecisionAuthorityValid = validation.valid;
  productImpactBaseErrors.push(
    ...formatProductValidationErrors('base Product Impact decisions', validation)
  );
}
if (candidateProductDecisionAuthority && candidateProductImpactMapValid) {
  const validation = validateProductDecisionAuthority(
    candidateProductDecisionAuthority,
    candidateProductImpactMap
  );
  candidateProductDecisionAuthorityValid = validation.valid;
  productImpactCandidateErrors.push(
    ...formatProductValidationErrors('candidate Product Impact decisions', validation)
  );
}

const trustedBasePullRequest = githubEventName === 'pull_request_target'
  && pullRequestProductImpactMetadata.trustedToDev;
const productImpactRuntimeStatus = trustedBasePullRequest
  ? 'TRUSTED_BASE_PULL_REQUEST'
  : isPullRequestEvent(githubEventName)
    ? 'NOT_AUTHORITATIVE_CANDIDATE'
    : ['push', 'workflow_dispatch'].includes(githubEventName)
      ? 'TRUSTED_INTEGRATED_DIFF'
      : 'TRUSTED_LOCAL_DIFF';
const parseCurrentDecisionSource = ({ body = '', authorLogin = '', headSha = '' } = {}) => (
  parseCurrentProductImpactDecisions({
    body,
    eventAuthorLogin: authorLogin,
    eventHeadSha: headSha
  })
);
let currentProductImpactDecisions = parseCurrentDecisionSource();
if (productImpactRuntimeStatus === 'TRUSTED_BASE_PULL_REQUEST') {
  currentProductImpactDecisions = parseCurrentDecisionSource({
    body: pullRequestProductImpactMetadata.body,
    authorLogin: pullRequestProductImpactMetadata.authorLogin,
    headSha: pullRequestProductImpactMetadata.headSha
  });
}
if (!currentProductImpactDecisions.valid) {
  const staleOnly = currentProductImpactDecisions.errors.length > 0
    && currentProductImpactDecisions.errors.every(({ code }) => code === 'stale-head-sha');
  if (staleOnly) {
    productImpactDecisionDeclarationIssues.push(
      'Observation: INSUFFICIENT_EVIDENCE. The current decision is stale for the event head SHA.'
    );
    currentProductImpactDecisions = parseCurrentDecisionSource({
      authorLogin: pullRequestProductImpactMetadata.authorLogin,
      headSha: pullRequestProductImpactMetadata.headSha
    });
  } else {
    productImpactBaseErrors.push(...currentProductImpactDecisions.errors.map((error) => (
      `current Product Impact decision ${error.path} [${error.code}]: ${error.message}`
    )));
  }
}

const baseRuleFactsByKey = new Map(
  baseArchitectureRuleFacts.map((fact) => [fact.ruleKey, fact])
);
const headRuleFactsByKey = new Map(
  headArchitectureRuleFacts.map((fact) => [fact.ruleKey, fact])
);
const productImpactRuleFacts = [];
const productImpactSubjectDeltas = [];
if (baseArchitectureRegistry) {
  [...baseArchitectureRegistry.rules]
    .sort((left, right) => left.key.localeCompare(right.key))
    .forEach((rule) => {
      const before = baseRuleFactsByKey.get(rule.key);
      const after = headRuleFactsByKey.get(rule.key);
      if (!before || !after) {
        productImpactBaseErrors.push(
          `Product Impact source facts are missing for architecture rule ${rule.key}`
        );
        return;
      }
      productImpactRuleFacts.push(Object.freeze({
        ruleKey: rule.key,
        subjectCategory: before.subjectCategory,
        beforeSubjects: before.subjects,
        afterSubjects: after.subjects
      }));
      try {
        productImpactSubjectDeltas.push(createArchitectureSubjectDelta({
          ruleKey: rule.key,
          kind: rule.kind,
          detector: rule.detector,
          subjectCategory: before.subjectCategory,
          beforeSubjects: before.subjects,
          afterSubjects: after.subjects
        }));
      } catch (error) {
        productImpactBaseErrors.push(
          `Product Impact subject delta failed for ${rule.key}: ${error.message}`
        );
      }
    });
}

let productImpactResults = [];
if (
  productImpactRuntimeStatus !== 'NOT_AUTHORITATIVE_CANDIDATE'
  && baseProductImpactMapValid
  && baseProductDecisionAuthorityValid
  && !productImpactBaseErrors.length
) {
  try {
    productImpactResults = evaluateProductImpact({
      subjectDeltas: productImpactSubjectDeltas,
      ruleFacts: productImpactRuleFacts,
      map: baseProductImpactMap,
      decisionAuthority: baseProductDecisionAuthority,
      currentDecisions: currentProductImpactDecisions,
      changedPaths
    });
  } catch (error) {
    productImpactBaseErrors.push(`Product Impact runtime evaluation failed: ${error.message}`);
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

const sizeReviewReasons = [];
if (productionPaths.length > selectedProfile.productionFiles) {
  sizeReviewReasons.push(
    `production files changed exceed ${selectedProfileName} size-review threshold `
    + `(${productionPaths.length} > ${selectedProfile.productionFiles})`
  );
}
if (productionGrossChurn > selectedProfile.grossChurn) {
  sizeReviewReasons.push(
    `production gross churn exceeds ${selectedProfileName} size-review threshold `
    + `(${productionGrossChurn} > ${selectedProfile.grossChurn})`
  );
}
if (productionNetAdditions > selectedProfile.netAdditions) {
  sizeReviewReasons.push(
    `production net additions exceed ${selectedProfileName} size-review threshold `
    + `(${productionNetAdditions} > ${selectedProfile.netAdditions})`
  );
}

const deltaChanged = (delta) => delta.added.length || delta.removed.length;
const architectureSignalDeltas = [
  ['registered authority locations', registeredAuthorityLocationDelta],
  ['registered canonical contracts', registeredCanonicalContractDelta],
  ['privileged operator permissions', privilegedOperatorPermissionDelta],
  ['privileged importer permissions', privilegedImporterPermissionDelta],
  ['production modules', sourceInventoryDeltas.modules],
  ['public exports', sourceInventoryDeltas.exports],
  ['create* declarations', sourceInventoryDeltas.createDeclarations],
  ['reactive declarations', sourceInventoryDeltas.reactiveDeclarations],
  ['watcher calls', sourceInventoryDeltas.watcherCalls],
  ['resource-like declarations', sourceInventoryDeltas.resourceNames],
  ['canonical-entry import edges', privilegedImportFanOutDelta],
  ['session schema fields', sourceInventoryDeltas.sessionObjectKeys],
  ['compatibility-like paths', sourceInventoryDeltas.compatibilityNames]
];
const changedArchitectureSignals = architectureSignalDeltas
  .filter(([, delta]) => deltaChanged(delta))
  .map(([name]) => name);
const performanceEvidencePaths = new Set([
  'tests/test_gallery_publication_performance.py',
  'tests/test_performance_short_circuits.py',
  'tests/web/interactive-svg-search-performance.playwright.spec.js',
  'tests/web/webapp-performance.playwright.spec.js',
  'tools/measure_gallery_publication_performance.py'
]);
const changedPerformanceEvidencePaths = changedPaths.filter((path) => (
  path.startsWith('tests/performance/') || performanceEvidencePaths.has(path)
));
const changedReferenceOutputPaths = changedPaths.filter((path) => (
  path.startsWith('tests/reference_outputs/')
));
const changedSessionContractPaths = changedPaths.filter((path) => (
  isWebSessionSourcePath(path)
  || path.startsWith('gbdraw/web/gallery/sessions/')
  || path.startsWith('tests/fixtures/sessions/')
));
const reviewReasons = sizeReviewReasons.map((reason) => `Size and scope: ${reason}`);
if (changedCheckerImplementations.length) {
  reviewReasons.push(
    `Governance and authority: Web checker implementation changed (${changedCheckerImplementations.join(', ')})`
  );
}
if (changedGovernancePaths.length) {
  reviewReasons.push(
    `Governance and authority: registered governance or authority files changed (${changedGovernancePaths.join(', ')})`
  );
}
if (policyContractions.length || policyExpansions.length) {
  reviewReasons.push(
    'Governance and authority: privileged policy permissions changed '
    + `(added=${policyExpansions.length}, removed=${policyContractions.length})`
  );
}
if (architectureAuthorityDelta.changes.length) {
  reviewReasons.push(
    'Governance and authority: registered architecture-rule authority changed '
    + `(${architectureAuthorityDelta.changes.length} change(s))`
  );
}
if (changedArchitectureSignals.length) {
  reviewReasons.push(
    `Architecture-bearing signals: deterministic inventories changed (${changedArchitectureSignals.join(', ')})`
  );
}
if (changedPerformanceEvidencePaths.length) {
  reviewReasons.push(
    `Material behavior/output risk: registered performance evidence paths changed (${changedPerformanceEvidencePaths.join(', ')})`
  );
}
if (changedReferenceOutputPaths.length) {
  reviewReasons.push(
    `Material behavior/output risk: reference output paths changed (${changedReferenceOutputPaths.join(', ')})`
  );
}
if (changedSessionContractPaths.length) {
  reviewReasons.push(
    `Material behavior/output risk: registered session or compatibility paths changed (${changedSessionContractPaths.join(', ')})`
  );
}
if (changedProductImpactAuthorityPaths.length) {
  reviewReasons.push(
    'Governance and authority: Product Impact authority changed for future revisions '
    + `(${changedProductImpactAuthorityPaths.join(', ')})`
  );
}
if (productImpactDecisionDeclarationIssues.length) {
  reviewReasons.push(
    `Product impact: current decision declaration requires attention (${productImpactDecisionDeclarationIssues.length} issue(s))`
  );
}
productImpactResults.forEach((result) => {
  const candidateModifiedContracts = result.contractResults
    .filter(({ candidateModified }) => candidateModified)
    .map(({ ref }) => ref);
  if (
    result.observation !== 'CONFORMING'
    || result.impactClass !== 'NO_USER_VISIBLE_DIFFERENCE'
    || result.currentDecisionIssues.length
    || candidateModifiedContracts.length
  ) {
    reviewReasons.push(
      `Product impact: ${result.concernKey} is ${result.observation} `
      + `with ${result.impactClass}`
    );
  }
  if (candidateModifiedContracts.length) {
    reviewReasons.push(
      `Product impact: mapped contract changed in the candidate (${candidateModifiedContracts.join(', ')})`
    );
  }
});

const blockingViolations = [
  ...changeContext.errors,
  ...acceptedViolationAuthorityErrors,
  ...unsafeTrustedWorkflowChanges,
  ...productImpactBaseErrors.map((error) => `trusted Product Impact authority: ${error}`),
  ...productImpactCandidateErrors.map((error) => `candidate Product Impact authority: ${error}`)
];
if (newProductionDependencies.length) {
  blockingViolations.push('new production dependencies are not allowed');
}
if (addedBinaryRuntimePaths.length) {
  blockingViolations.push('added binary runtime files are not allowed');
}
if (changedVendorPaths.length) {
  blockingViolations.push('changes under gbdraw/web/vendor/ are not allowed');
}
if (productionPaths.length && changedGuards.length && !pureSafePolicyContraction) {
  blockingViolations.push('production runtime files and Web guard/CI files changed together');
}
if (changedCheckerImplementations.length && changedAuthorities.length) {
  blockingViolations.push(
    'Web checker/source parser and authority policy/workflow files changed together'
  );
}
if (unapprovedCapabilities.length) {
  blockingViolations.push('privileged capability owners or importers exceed the base allowlist');
}
if (proposedPolicyExclusions.length) {
  blockingViolations.push(
    'proposed privileged capability allowlist excludes active owners or importers'
  );
}
if (missingPolicyKeys.length) {
  blockingViolations.push(
    'proposed privileged capability policy is missing base allowlist keys'
  );
}
if (baseImportGraph.errors.length) {
  blockingViolations.push('trusted base first-party static import graph is incomplete');
}
if (headImportGraph.errors.length) {
  blockingViolations.push('head first-party static import graph is incomplete');
}
if (baseCycleResult.cycles.length) {
  blockingViolations.push('trusted base contains first-party static import cycles');
}
if (headCycleResult.cycles.length) {
  blockingViolations.push('first-party static import cycles are not allowed');
}
blockingViolations.push(...candidateArchitectureErrors.map((error) => (
  `candidate architecture rules: ${error}`
)));
blockingViolations.push(...trustedBaseArchitectureErrors.map((error) => (
  `trusted base architecture rules: ${error}`
)));
blockingViolations.push(...trustedBaseArchitectureFailures.map((error) => (
  `trusted base architecture rules: ${error}`
)));
blockingViolations.push(...activeArchitectureErrors.map((error) => (
  `active architecture rules: ${error}`
)));
blockingViolations.push(...activeArchitectureFailures.map((error) => (
  `active architecture rules: ${error}`
)));
const productImpactCoverageByRule = new Map(
  (baseProductImpactMap?.ruleCoverage || []).map((coverage) => [coverage.ruleKey, coverage])
);
const productImpactResultRuleKeys = (result) => {
  if (result.triggeringRuleKeys.length) return result.triggeringRuleKeys;
  const concern = (baseProductImpactMap?.concerns || []).find(({ key }) => (
    key === result.concernKey
  ));
  return stableReasons((concern?.options || []).flatMap((option) => (
    option.requirements.flatMap((requirement) => (
      requirement.anyOf.map(({ ruleKey }) => ruleKey)
    ))
  )));
};
productImpactResults.forEach((result) => {
  const hard = productImpactResultRuleKeys(result).some((ruleKey) => (
    productImpactCoverageByRule.get(ruleKey)?.enforcement === 'hard'
  ));
  if (
    hard
    && [
      'AUTHORITY_CONFLICT',
      'INSUFFICIENT_EVIDENCE',
      'ORDINARY_REGRESSION',
      'UNRESOLVED_DECISION'
    ].includes(result.observation)
  ) {
    blockingViolations.push(
      `Product Impact hard enforcement: ${result.concernKey} is ${result.observation}`
    );
  }
});
const policyResult = createPolicyResult({
  blockingViolations,
  reviewReasons,
  context: changeContext.context,
  observations: {
    architectureAuthorityDelta,
    changedArchitectureSignals,
    productImpactResults,
    productImpactRuntimeStatus,
    productionGrossChurn,
    productionNetAdditions,
    productionPaths,
    selectedProfile: selectedProfileName
  }
});

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
    classification: 'size-review'
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
const productConcernByKey = new Map(
  (baseProductImpactMap?.concerns || []).map((concern) => [concern.key, concern])
);
const productDecisionOwnerLogin = currentProductImpactDecisions.metadata.eventAuthorLogin;
const productDecisionOwnerEligible = Boolean(
  productDecisionOwnerLogin
  && baseProductDecisionAuthority?.maintainerLogins?.includes(productDecisionOwnerLogin)
);
const stableIntersection = (left, right) => {
  const rightSet = new Set(right);
  return left.filter((value) => rightSet.has(value));
};
const productImpactGateContribution = (result) => {
  const hard = productImpactResultRuleKeys(result).some((ruleKey) => (
    productImpactCoverageByRule.get(ruleKey)?.enforcement === 'hard'
  ));
  const blocking = [
    'AUTHORITY_CONFLICT',
    'INSUFFICIENT_EVIDENCE',
    'ORDINARY_REGRESSION',
    'UNRESOLVED_DECISION'
  ].includes(result.observation);
  return hard && blocking ? 'FAIL' : 'None (report-only for the triggering rule set)';
};
const renderDecisionPack = (result, concern) => {
  const optionsById = new Map(concern.options.map((option) => [option.id, option]));
  const optionResultsById = new Map(
    result.optionResults.map((option) => [option.optionId, option])
  );
  const beforeOptionId = result.beforeOptions[0] || '';
  const afterOptionId = result.afterOptions[0] || '';
  const beforeEffects = result.beforeEffects;
  const routeForMappedOption = (optionId) => {
    const optionResult = optionResultsById.get(optionId);
    if (!optionResult?.afterRealized) return 'NOT_ALLOWED';
    const effectEquivalent = JSON.stringify([...optionResult.effectRefs].sort())
      === JSON.stringify([...beforeEffects].sort());
    if (!effectEquivalent) return 'DURABLE_AUTHORITY_REQUIRED';
    return productDecisionOwnerEligible
      ? 'PR_LOCAL_ALLOWED'
      : 'DURABLE_AUTHORITY_REQUIRED';
  };
  const choices = [
    {
      code: 'A',
      optionId: beforeOptionId,
      outcome: beforeOptionId
        ? `Keep the before behavior: ${optionsById.get(beforeOptionId)?.summary || beforeOptionId}`
        : 'Keep the before behavior.',
      route: beforeOptionId ? routeForMappedOption(beforeOptionId) : 'NOT_ALLOWED'
    },
    {
      code: 'B',
      optionId: afterOptionId,
      outcome: afterOptionId
        ? `Select the after behavior: ${optionsById.get(afterOptionId)?.summary || afterOptionId}`
        : 'Select the after behavior.',
      route: afterOptionId ? routeForMappedOption(afterOptionId) : 'NOT_ALLOWED'
    },
    {
      code: 'C',
      optionId: '',
      outcome: 'Define a new mapped option that combines useful properties.',
      route: 'DURABLE_AUTHORITY_REQUIRED'
    },
    {
      code: 'D',
      optionId: '',
      outcome: 'Preserve both outcomes as distinct supported workflows.',
      route: 'DURABLE_AUTHORITY_REQUIRED'
    },
    {
      code: 'E',
      optionId: '',
      outcome: 'Intentionally retire a supported behavior or affordance.',
      route: 'DURABLE_AUTHORITY_REQUIRED'
    },
    {
      code: 'F',
      optionId: '',
      outcome: 'Defer the convergence and gather the missing evidence.',
      route: 'EVIDENCE_REQUIRED'
    }
  ];
  const renderOptionEffects = (optionId) => {
    const effectRefs = optionsById.get(optionId)?.effectRefs || [];
    return effectRefs.length
      ? effectRefs.map((effectRef) => (
        `${inlineCode(effectRef)}: ${reportSentenceFragment(
          concern.effects.find(({ id }) => id === effectRef)?.statement
        )}`
      )).join('; ')
      : 'not mapped';
  };
  return [
    '',
    '#### Decision Pack',
    '',
    '- Why a decision is required: the mapped options differ and no eligible authority selects the proposed outcome.',
    `- Actor: ${safeReportText(result.journeys[0]?.actor)}`,
    `- Context: ${safeReportText(result.journeys[0]?.context)}`,
    `- Goal: ${safeReportText(result.journeys[0]?.goal)}`,
    `- Action: ${(result.journeys[0]?.steps || []).map(safeReportText).join(' -> ')}`,
    `- Before option(s): ${result.beforeOptions.length ? result.beforeOptions.map(inlineCode).join(', ') : 'None'}`,
    `- After option(s): ${result.afterOptions.length ? result.afterOptions.map(inlineCode).join(', ') : 'None'}`,
    `- Authority searched: ${result.authorityResolution}; no eligible selecting authority was found.`,
    `- Product Decision Owner: base allowlisted maintainer; event author ${productDecisionOwnerLogin ? inlineCode(productDecisionOwnerLogin) : 'is unavailable'} (${productDecisionOwnerEligible ? 'eligible' : 'not eligible for PR-local routing'}).`,
    '- Choices:',
    ...choices.map(({ code, optionId, outcome, route }) => (
      `  - ${code}. ${safeReportText(outcome)} Option ID: ${optionId ? inlineCode(optionId) : 'not mapped'}. Effects: ${renderOptionEffects(optionId)}. Route: ${inlineCode(route)}.`
    )),
    '- Removal consequences:',
    `  - Lost effects: ${result.lostEffectRefs.length ? result.lostEffectRefs.map(inlineCode).join(', ') : 'None'}`,
    `  - Lost affordance/compatibility requirements: ${result.lostRequirementRefs.length ? result.lostRequirementRefs.map(inlineCode).join(', ') : 'None'}`,
    '- Route meanings:',
    '  - `PR_LOCAL_ALLOWED`: Product Decision Owner response -> Codex serializes an exact-head PR-body decision.',
    '  - `DURABLE_AUTHORITY_REQUIRED`: Product Decision Owner response -> Codex prepares a narrow authority-only PR -> merge -> rebase implementation.',
    '  - `EVIDENCE_REQUIRED`: stop the runtime convergence and collect evidence or merge an evidence-only PR.',
    '  - `NOT_ALLOWED`: revise the implementation; product authority cannot waive the failing boundary.',
    '- The Product Decision Owner returns this short response to Codex. The human does not edit JSON, SHA, requirement refs, or evidence refs; Codex serializes only the explicit choice and the trusted checker validates it.',
    '',
    '```text',
    'PRODUCT_DECISION',
    `Concern: ${result.concernKey}`,
    `Scenario revision: ${result.scenarioRevision}`,
    'Choice:',
    'Rationale:',
    'Must preserve:',
    'May retire:',
    'Accepted residual risk:',
    '```'
  ];
};
const renderProductImpactPacket = (result) => {
  const concern = productConcernByKey.get(result.concernKey);
  if (!concern) return [];
  const effectsById = new Map(concern.effects.map((effect) => [effect.id, effect]));
  const preservedEffects = stableIntersection(result.beforeEffects, result.afterEffects);
  const authorityRefs = concern.resolution.authorityRefs || [];
  const durableDecision = (baseProductDecisionAuthority?.decisions || []).find((decision) => (
    decision.concernKey === result.concernKey
    && decision.scenarioRevision === result.scenarioRevision
  ));
  const currentDecision = (currentProductImpactDecisions.declaration.decisions || [])
    .find((decision) => (
      decision.concernKey === result.concernKey
      && decision.scenarioRevision === result.scenarioRevision
    ));
  const authoritySelections = [
    ...(['authority-covered', 'domain-derived'].includes(concern.resolution.kind) ? [{
      type: concern.resolution.kind === 'domain-derived' ? 'DOMAIN_AUTHORITY' : 'STATIC_AUTHORITY',
      selectedOptionId: concern.resolution.selectedOptionId
    }] : []),
    ...(durableDecision ? [{
      type: 'DURABLE_DECISION',
      selectedOptionId: durableDecision.selectedOptionId
    }] : []),
    ...(currentDecision ? [{
      type: 'CURRENT_MAINTAINER_DECISION',
      selectedOptionId: currentDecision.selectedOptionId
    }] : [])
  ];
  const lines = [
    `### ${result.concernKey}`,
    '',
    `- Scenario revision: ${result.scenarioRevision}`,
    `- Layer: ${concern.layer}`,
    `- Product Decision Owner eligibility: ${productDecisionOwnerEligible ? 'eligible base maintainer' : 'no eligible PR-local author in this context'}`,
    '',
    '#### Architecture delta',
    '',
    ...result.subjectDeltas.flatMap((delta) => [
      `- Rule: ${inlineCode(delta.ruleKey)}`,
      `  - Added: ${delta.addedSubjects.length ? delta.addedSubjects.map(inlineCode).join(', ') : 'None'}`,
      `  - Removed: ${delta.removedSubjects.length ? delta.removedSubjects.map(inlineCode).join(', ') : 'None'}`
    ]),
    '',
    '#### Journeys and checkpoints',
    '',
    ...result.journeys.flatMap((journey) => [
      `- ${inlineCode(journey.id)}: ${safeReportText(journey.actor)}; ${safeReportText(journey.context)}`,
      `  - Goal: ${safeReportText(journey.goal)}`,
      `  - Checkpoints: ${journey.checkpoints.map(({ id }) => inlineCode(`${journey.id}:${id}`)).join(', ')}`
    ]),
    '',
    '#### Requirement realization',
    '',
    '| Requirement | Category | Before | After | Before providers | After providers |',
    '| --- | --- | --- | --- | --- | --- |',
    ...result.requirementResults.map((requirement) => (
      `| ${inlineCode(requirement.requirementRef)} | ${requirement.category} | ${requirement.beforeSatisfied} | ${requirement.afterSatisfied} | `
      + `${requirement.beforeActiveSubjects.length ? requirement.beforeActiveSubjects.map(inlineCode).join(', ') : 'None'} | `
      + `${requirement.afterActiveSubjects.length ? requirement.afterActiveSubjects.map(inlineCode).join(', ') : 'None'} |`
    )),
    '',
    '#### Behavior options and effects',
    '',
    `- Before realized options: ${result.beforeOptions.length ? result.beforeOptions.map(inlineCode).join(', ') : 'None'}`,
    `- After realized options: ${result.afterOptions.length ? result.afterOptions.map(inlineCode).join(', ') : 'None'}`,
    `- Preserved effects: ${preservedEffects.length ? preservedEffects.map((id) => `${inlineCode(id)}: ${reportSentenceFragment(effectsById.get(id)?.statement)}`).join('; ') : 'None'}`,
    `- Added effects: ${result.addedEffectRefs.length ? result.addedEffectRefs.map(inlineCode).join(', ') : 'None'}`,
    `- Lost effects: ${result.lostEffectRefs.length ? result.lostEffectRefs.map(inlineCode).join(', ') : 'None'}`,
    '',
    '#### Authority and contracts',
    '',
    `- Resolution: ${result.authorityResolution}`,
    `- Selected option: ${result.selectedOptionId ? inlineCode(result.selectedOptionId) : 'None'}`,
    `- Authority selections: ${authoritySelections.length ? authoritySelections
      .map(({ type, selectedOptionId }) => (
        `${inlineCode(type)} -> ${inlineCode(selectedOptionId)}`
      )).join('; ') : 'None'}`,
    `- Authority references: ${authorityRefs.length ? authorityRefs.map(inlineCode).join(', ') : 'None'}`,
    `- Current-decision rationale: ${result.decisionRationale ? safeReportText(result.decisionRationale) : 'None'}`,
    ...result.contractResults.map((contract) => (
      `- Contract ${inlineCode(contract.ref)}: sensitivity=${contract.sensitivity}; execution=${contract.execution}; integrity=${contract.candidateModified ? 'CANDIDATE_MODIFIED' : 'UNCHANGED'}`
    )),
    ...result.residualRisks.map((risk) => `- Residual risk: ${safeReportText(risk)}`),
    ...result.evidenceGaps.map((gap) => `- Evidence gap: ${safeReportText(gap)}`),
    ...result.currentDecisionIssues.map((issue) => `- Current decision issue: ${safeReportText(issue)}`),
    '',
    '#### Result',
    '',
    `- Impact: ${result.impactClass}`,
    `- Observation: ${result.observation}`,
    `- Product decision required: ${result.decisionRequired ? 'yes' : 'no'}`,
    `- Gate contribution: ${productImpactGateContribution(result)}`,
    `- Next action: ${safeReportText(result.nextAction)}`
  ];
  if (result.decisionRequired) lines.push(...renderDecisionPack(result, concern));
  return lines;
};
const productImpactAuthorityProposal = Boolean(
  changedProductImpactAuthorityPaths.length || architectureAuthorityDelta.changes.length
);
const productImpactReportUseful = Boolean(
  productImpactResults.length
  || productImpactDecisionDeclarationIssues.length
  || productImpactAuthorityProposal
);
const productImpactReportLines = productImpactReportUseful ? [
  '## Product impact',
  '',
  `- Runtime context: ${productImpactRuntimeStatus}`,
  '- Runtime authority: trusted base map, decisions, architecture rules, detectors, and checker only.',
  `- Base authority validation: ${productImpactBaseErrors.length ? 'INVALID' : 'VALID'}`,
  `- Candidate authority validation: ${candidateProductImpactMapValid && candidateProductDecisionAuthorityValid ? 'VALID (inert data only)' : 'INVALID'}`,
  `- Candidate authority separation: ${candidateAuthorityCompanionPaths.length ? 'INVALID' : 'VALID'}`,
  `- Runtime result count: ${productImpactResults.length}`,
  ...(productImpactRuntimeStatus === 'NOT_AUTHORITATIVE_CANDIDATE' ? [
    '- This pull request context is non-authoritative and emits no runtime admission packet.'
  ] : []),
  ...(productImpactAuthorityProposal ? [
    '',
    '### Candidate authority proposal',
    '',
    `- Changed Product Impact authority paths: ${changedProductImpactAuthorityPaths.length ? changedProductImpactAuthorityPaths.map(inlineCode).join(', ') : 'None'}`,
    `- Candidate map: ${candidateProductImpactMapValid ? 'VALID' : 'INVALID'}`,
    `- Candidate decision registry: ${candidateProductDecisionAuthorityValid ? 'VALID' : 'INVALID'}`,
    '- Runtime effect: validation-only future preauthorization; candidate data does not alter this head runtime admission.'
  ] : []),
  ...productImpactDecisionDeclarationIssues.flatMap((issue) => [
    '',
    '### Current decision declaration',
    '',
    `- ${safeReportText(issue)}`,
    '- Gate contribution: None (report-only).',
    '- Next action: the Product Decision Owner must review the exact current head and Codex must serialize a renewed declaration.'
  ]),
  ...productImpactResults.flatMap((result) => ['', ...renderProductImpactPacket(result)]),
  ''
] : [];
const report = [
  '# Web change policy',
  '',
  `- Context: ${policyResult.context}`,
  `- Base: \`${base}\``,
  `- Head: \`${head || 'working tree'}\``,
  `- Gate: **${policyResult.gate}**`,
  `- Review: **${policyResult.review}**`,
  `- Privileged allowlist revision: \`${policyRevision}\``,
  `- Architecture rule base: ${baseArchitectureRulesSource === null ? 'absent' : `\`${base}\``}`,
  `- Architecture rule candidate: ${proposedArchitectureRulesSource === null ? 'absent' : `\`${head || 'working tree'}\``}`,
  '- Policy guide: `docs/internal/WEB_CHANGE_POLICY.md`',
  `- Selected profile: ${selectedProfileName}`,
  `- Size-review threshold for production files: ${selectedProfile.productionFiles}`,
  `- Size-review threshold for gross churn: ${selectedProfile.grossChurn}`,
  `- Size-review threshold for net additions: ${selectedProfile.netAdditions}`,
  `- \`architecture-change\` label: ${architectureChange ? 'present' : 'absent'}`,
  '',
  '## Blocking violations',
  '',
  ...list(policyResult.blockingViolations),
  '',
  '## Review reasons',
  '',
  ...list(policyResult.reviewReasons),
  '',
  ...productImpactReportLines,
  '## Key architecture differential summary',
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
  ''
].join('\n');

emitPolicyReport(policyResult, report);
};

if (changeContext.isPromotion || changeContext.errors.length) {
  renderBoundedContextReport();
} else {
  try {
    runOrdinaryPolicy();
  } catch (error) {
    const errorMessage = error?.message || String(error);
    const result = createPolicyResult({
      blockingViolations: [
        `Web policy evaluation failed closed: ${errorMessage}`
      ],
      reviewReasons: /tools\/web-(?:change-policy|architecture-(?:rules|violations))\.json/.test(
        errorMessage
      ) || /tools\/web-product-(?:impact-map|decisions)\.json/.test(errorMessage)
        ? ['Governance and authority: required authority input is malformed']
        : [],
      context: changeContext.context,
      observations: { errorType: error?.name || typeof error }
    });
    const report = [
      '# Web change policy',
      '',
      `- Context: ${result.context}`,
      `- Base: \`${base}\``,
      `- Head: \`${head || 'working tree'}\``,
      `- Gate: **${result.gate}**`,
      `- Review: **${result.review}**`,
      '',
      '## Blocking violations',
      '',
      ...list(result.blockingViolations),
      '',
      '## Review reasons',
      '',
      ...list(result.reviewReasons),
      ''
    ].join('\n');
    emitPolicyReport(result, report);
  }
}
