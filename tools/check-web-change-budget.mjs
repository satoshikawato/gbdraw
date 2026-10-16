#!/usr/bin/env node

import { appendFileSync, existsSync, readFileSync } from 'node:fs';
import { join, posix, resolve } from 'node:path';
import { execFileSync, spawnSync } from 'node:child_process';
import { createHash } from 'node:crypto';
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
import { detectWebOwnerGraph, summarizeWebOwnerGraph } from './web-owner-graph-detectors.mjs';
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
const productContractAuthorityPath = 'docs/internal/OPTION_INTEGRITY_PRODUCT_CONTRACT.md';
const productImpactMapPath = 'tools/web-product-impact-map.json';
const productDecisionAuthorityPath = 'tools/web-product-decisions.json';
const productImpactEvaluationPath = 'tools/web-product-impact-evaluation.mjs';
const productImpactDecisionSourcePath = 'tools/web-product-impact-decision-source.mjs';
const productImpactFixturePath = 'tests/web/product-impact-ratchet-fixtures.test.mjs';
const trustedWorkflowPath = '.github/workflows/web-base-policy.yml';
// The Web design rules (the numbered rules in gbdraw/web/CLAUDE.md) and the registry that
// names their guard tests and in-test allowlists are authority: a runtime
// change may not amend them in the same diff (WEB_CHANGE_POLICY.md "Design-rule
// co-change").
const webDesignRulesPath = 'gbdraw/web/CLAUDE.md';
const designRuleGuardRegistryPath = 'tools/web-design-rule-guards.json';
// The owner-graph registry (tools/web-owner-graph.json, implementation plan
// Phase B2) is authority; its detectors are checker implementation and report
// only until a rule references them.
const ownerGraphRegistryPath = 'tools/web-owner-graph.json';
const ownerGraphDetectorsPath = 'tools/web-owner-graph-detectors.mjs';
// The twin detector (tools/twins.json is its baseline) is checker implementation.
const twinDetectorsPath = 'tools/web-twin-detectors.mjs';
const guardPaths = new Set([
  'docs/internal/ARCHITECTURE_FITNESS_FUNCTION_RATCHET.md',
  productImpactPolicyPath,
  productContractAuthorityPath,
  webDesignRulesPath,
  designRuleGuardRegistryPath,
  ownerGraphRegistryPath,
  ownerGraphDetectorsPath,
  twinDetectorsPath,
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
  ownerGraphDetectorsPath,
  twinDetectorsPath,
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
  productContractAuthorityPath,
  webDesignRulesPath,
  designRuleGuardRegistryPath,
  ownerGraphRegistryPath,
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

// Design-rule guard registry (tools/web-design-rule-guards.json), read from the
// trusted base. Each rule names its guard tests and the allowlist or baseline
// literals inside them. An allowlist may only contract in a diff that also
// changes production runtime files; an expansion needs an authority-only PR.
// Absent registry: inert (only the path pre-registration above applies).
// Schema 1 registers bare guard paths, `R<n>` rules and four literal kinds.
// Schema 2 records per guard its root cause, CI job and retirement condition,
// lets a rule outside gbdraw/web/CLAUDE.md name its own rulesPath, and adds
// the kinds read from a path -> names map, a JSON file, or a Python module.
const DESIGN_RULE_ALLOWLIST_KINDS = {
  1: ['count', 'count-map', 'set', 'writer-map'],
  2: ['count', 'count-map', 'set', 'writer-map', 'set-map', 'json-entries', 'json-ids', 'py-count', 'py-count-map', 'py-set']
};
// The schema 1 comparison that each later kind gets once its literal is read.
const DESIGN_RULE_ALLOWLIST_COMPARISONS = {
  'set-map': 'writer-map',
  'json-entries': 'writer-map',
  'json-ids': 'writer-map',
  'py-count': 'count',
  'py-count-map': 'count-map',
  'py-set': 'set'
};
const isPythonAllowlistKind = (kind) => kind.startsWith('py-');
// A JSON list of entries with ids: `json-entries` compares each entry's text,
// `json-ids` the ids and the length of each row's lists (a shrink-only id set
// whose lists only shrink; a renamed item and the other fields may change).
const JSON_ALLOWLIST_KINDS = new Set(['json-entries', 'json-ids']);
const parseDesignRuleGuardRegistry = (source) => {
  const errors = [];
  if (source === null) return { registry: null, errors };
  let parsed;
  try {
    parsed = JSON.parse(source);
  } catch (error) {
    return { registry: null, errors: [`${designRuleGuardRegistryPath}: invalid JSON (${error.message})`] };
  }
  if (!parsed || typeof parsed !== 'object' || Array.isArray(parsed)) {
    return { registry: null, errors: [`${designRuleGuardRegistryPath}: expected an object`] };
  }
  const schemaVersion = parsed.schemaVersion === 2 ? 2 : 1;
  if (parsed.schemaVersion !== 1 && parsed.schemaVersion !== 2) {
    errors.push(`${designRuleGuardRegistryPath}: schemaVersion must be 1 or 2`);
  }
  if (parsed.rulesPath !== webDesignRulesPath) {
    errors.push(`${designRuleGuardRegistryPath}: rulesPath must be ${webDesignRulesPath}`);
  }
  const rules = Array.isArray(parsed.rules) ? parsed.rules : null;
  if (!rules) errors.push(`${designRuleGuardRegistryPath}: rules must be an array`);
  const seenIds = new Set();
  const isText = (value, pattern) => typeof value === 'string' && pattern.test(value);
  (rules || []).forEach((rule, index) => {
    const location = `${designRuleGuardRegistryPath}: rules[${index}]`;
    if (!rule || typeof rule !== 'object') {
      errors.push(`${location}: expected an object`);
      return;
    }
    const ownRulesPath = schemaVersion === 2 && rule.rulesPath !== undefined;
    if (ownRulesPath && (typeof rule.rulesPath !== 'string' || !rule.rulesPath)) {
      errors.push(`${location}.rulesPath: expected a path`);
    }
    const idPattern = ownRulesPath ? /^[A-Z][A-Z0-9]*(-[A-Z0-9]+)*$/ : /^R\d+$/;
    if (typeof rule.id !== 'string' || !idPattern.test(rule.id)) {
      errors.push(`${location}.id: expected R<n>${schemaVersion === 2 ? ', or an upper-case id with its own rulesPath' : ''}`);
    } else if (seenIds.has(rule.id)) errors.push(`${location}.id: duplicate ${rule.id}`);
    else seenIds.add(rule.id);
    if (typeof rule.heading !== 'string' || !rule.heading.trim()) errors.push(`${location}.heading: expected text`);
    if (schemaVersion === 1) {
      if (!Array.isArray(rule.guards) || rule.guards.some((path) => typeof path !== 'string' || !path)) {
        errors.push(`${location}.guards: expected an array of paths`);
      }
    } else if (!Array.isArray(rule.guards)) errors.push(`${location}.guards: expected an array`);
    else {
      rule.guards.forEach((guard, guardIndex) => {
        const guardLocation = `${location}.guards[${guardIndex}]`;
        if (!guard || typeof guard !== 'object') {
          errors.push(`${guardLocation}: expected an object`);
          return;
        }
        if (!isText(guard.path, /./)) errors.push(`${guardLocation}.path: expected a path`);
        if (!isText(guard.rootCause, /^(?:P(?:[1-9]|10)|new)$/)) {
          errors.push(`${guardLocation}.rootCause: expected P1..P10 or new`);
        }
        if (!isText(guard.job, /^[a-z][a-z0-9-]*$/)) errors.push(`${guardLocation}.job: expected a CI job id`);
        if (!isText(guard.retireWhen, /\S/)) errors.push(`${guardLocation}.retireWhen: expected text`);
      });
    }
    if (!Array.isArray(rule.allowlists)) {
      errors.push(`${location}.allowlists: expected an array`);
      return;
    }
    rule.allowlists.forEach((entry, entryIndex) => {
      const entryLocation = `${location}.allowlists[${entryIndex}]`;
      if (!entry || typeof entry !== 'object') {
        errors.push(`${entryLocation}: expected an object`);
        return;
      }
      if (typeof entry.path !== 'string' || !entry.path) errors.push(`${entryLocation}.path: expected a path`);
      if (typeof entry.symbol !== 'string' || !/^[A-Za-z_$][\w$]*$/.test(entry.symbol)) {
        errors.push(`${entryLocation}.symbol: expected an identifier`);
      }
      const kinds = DESIGN_RULE_ALLOWLIST_KINDS[schemaVersion];
      if (!kinds.includes(entry.kind)) {
        errors.push(`${entryLocation}.kind: expected one of ${kinds.join(', ')}`);
      }
    });
  });
  return { registry: errors.length ? null : parsed, errors };
};

// Reads `const SYMBOL = <literal>` from a test source without executing it:
// a number, an array or `new Set([...])` of strings, or an object literal of
// numbers or of nested string maps, each optionally in `Object.freeze(...)`.
// Returns null when the literal is absent, is not one of those shapes, or is
// not closed by `;` (`[...].concat(MORE)`, also continued on the next line).
const matchingBracketEnd = (text, openIndex) => {
  const pairs = { '[': ']', '{': '}', '(': ')' };
  const stack = [];
  let quote = null;
  for (let index = openIndex; index < text.length; index++) {
    const character = text[index];
    if (quote) {
      if (character === '\\') index++;
      else if (character === quote) quote = null;
      continue;
    }
    if (character === "'" || character === '"' || character === '`') {
      quote = character;
      continue;
    }
    if (character === '/' && text[index + 1] === '/') {
      index = text.indexOf('\n', index);
      if (index < 0) return -1;
      continue;
    }
    if (character === '/' && text[index + 1] === '*') {
      index = text.indexOf('*/', index) + 1;
      if (index < 1) return -1;
      continue;
    }
    if (pairs[character]) stack.push(pairs[character]);
    else if (character === ']' || character === '}' || character === ')') {
      if (stack.pop() !== character) return -1;
      if (!stack.length) return index;
    }
  }
  return -1;
};
const literalToJson = (literal) => {
  let out = '';
  let quote = null;
  for (let index = 0; index < literal.length; index++) {
    const character = literal[index];
    if (quote) {
      if (character === '\\') {
        out += character + (literal[index + 1] ?? '');
        index++;
      } else if (character === quote) {
        out += '"';
        quote = null;
      } else out += character === '"' ? '\\"' : character;
      continue;
    }
    if (character === "'" || character === '"') {
      quote = character;
      out += '"';
      continue;
    }
    if (character === '/' && literal[index + 1] === '/') {
      const end = literal.indexOf('\n', index);
      index = end < 0 ? literal.length : end;
      continue;
    }
    if (character === '/' && literal[index + 1] === '*') {
      index = literal.indexOf('*/', index) + 1;
      continue;
    }
    out += character;
  }
  return out
    .replace(/([{,]\s*)([A-Za-z_$][\w$]*)\s*:/g, '$1"$2":')
    .replace(/,(\s*[\]}])/g, '$1');
};
const UNREADABLE_ALLOWLIST = Object.freeze({ status: 'unreadable', value: null });
const escapeAllowlistSymbol = (symbol) => symbol.replace(/[.*+?^${}()|[\]\\]/g, '\\$&');
const readJavaScriptLiteral = (source, symbol) => {
  const masked = maskJavaScript(source);
  const declaration = new RegExp(`^[ \\t]*(?:export\\s+)?const\\s+${escapeAllowlistSymbol(symbol)}\\s*=\\s*`, 'm');
  const match = masked.match(declaration);
  if (!match) return { status: 'missing', value: null };
  const start = match.index + match[0].length;
  const rest = source.slice(start);
  const numberMatch = rest.match(/^(\d+)\s*;/);
  if (numberMatch) {
    return { status: 'ok', value: Number(numberMatch[1]), span: [match.index, start + numberMatch[0].length] };
  }
  const wrappers = rest.match(/^(Object\.freeze\s*\(\s*)?(new\s+Set\s*\(\s*)?/);
  const literalStart = start + wrappers[0].length;
  const opener = source[literalStart];
  if (opener !== '[' && opener !== '{') return UNREADABLE_ALLOWLIST;
  const end = matchingBracketEnd(source, literalStart);
  if (end < 0) return UNREADABLE_ALLOWLIST;
  // Only the wrappers' closing parentheses and then `;` may follow the literal:
  // a line end does not close it (`.concat(MORE)` may start the next line).
  const closers = '\\s*\\)'.repeat(wrappers.slice(1).filter(Boolean).length);
  const closed = masked.slice(end + 1).match(new RegExp(`^${closers}\\s*;`));
  if (!closed) return UNREADABLE_ALLOWLIST;
  try {
    return {
      status: 'ok',
      value: JSON.parse(literalToJson(source.slice(literalStart, end + 1))),
      span: [match.index, end + 1 + closed[0].length]
    };
  } catch (_error) {
    return UNREADABLE_ALLOWLIST;
  }
};
// Blanks Python comments and string contents, offsets kept, so that only code
// can declare a baseline.
const maskPython = (source) => source.replace(
  /#[^\n]*|("""|''')[\s\S]*?\1|"(?:\\.|[^"\\\n])*"|'(?:\\.|[^'\\\n])*'/g,
  (text) => text.replace(/[^\n]/g, ' ')
);
// Reads the module-level `SYMBOL = <literal>` or `SYMBOL: <annotation> = <literal>`
// of a Python source without executing it: an integer, a string, a list, tuple,
// set or `frozenset(...)` (an array here), a dict with string keys (an object),
// or a bare name (null, as the note in `(10, _COMPARISON)`). Python has no
// `const`, so the literal is unreadable when the name is used but not assigned
// that way, when anything else in the file binds or changes the name (see
// `pythonRebinding`), when the value is any other expression, or when more
// follows it.
// Code that binds or changes NAME, at any indentation: a line or a `;`
// statement that starts with it, NAME as an assignment target (`X, NAME = `,
// `X = NAME = `, `NAME[key] = `, `NAME |= `, `NAME := `), a method call that is
// not a known read, and `import`, `as`, `for`, `del`, `def`, `class`, `global`
// or `nonlocal` naming it. Reads (`for key in NAME`, `len(NAME)`, `NAME.get(`,
// `f(NAME, strict=True)`) are not matched. Not seen: a change through an alias,
// `globals()` or a callee, a parenthesized target list, and another module.
const PYTHON_READ_METHODS = 'get|items|keys|values|copy|count|index|union|intersection|difference|'
  + 'symmetric_difference|isdisjoint|issubset|issuperset';
const pythonRebinding = (name) => {
  const bare = `(?<![\\w.])${name}(?!\\w)`;
  const assign = '(?:\\*\\*|\\/\\/|<<|>>|[-+*/%@|&^:])?=(?!=)';
  return new RegExp([
    `^[ \\t]*${bare}`,
    `;[ \\t]*${bare}`,
    `${bare}[ \\t]*(?:\\[[^\\]\\n]*\\][ \\t]*)?${assign}`,
    `^[ \\t]*[\\w.,\\[\\] \\t*]*?${bare}[\\w.,\\[\\] \\t*]*${assign}`,
    `${bare}\\.(?!(?:${PYTHON_READ_METHODS})\\b)\\w+[ \\t]*\\(`,
    `\\bimport\\b[^\\n]*${bare}`,
    `\\b(?:as|for|del|def|class|global|nonlocal)[ \\t]+${bare}`
  ].join('|'), 'm');
};
const readPythonLiteral = (source, symbol) => {
  const masked = maskPython(source);
  const name = escapeAllowlistSymbol(symbol);
  const declaration = masked.match(new RegExp(`^${name}[ \\t]*(?::[^=\\n]*)?=(?!=)[ \\t]*`, 'm'));
  if (!declaration) {
    return new RegExp(`(?<![\\w.])${name}(?!\\w)`).test(masked) ? UNREADABLE_ALLOWLIST : { status: 'missing', value: null };
  }
  const elsewhere = masked.slice(0, declaration.index) + ' '.repeat(declaration[0].length)
    + masked.slice(declaration.index + declaration[0].length);
  if (pythonRebinding(name).test(elsewhere)) return UNREADABLE_ALLOWLIST;
  let index = declaration.index + declaration[0].length;
  let depth = 0;
  const fail = () => {
    throw new SyntaxError(`unsupported Python literal for ${symbol}`);
  };
  const take = (pattern) => {
    pattern.lastIndex = index;
    const match = pattern.exec(source);
    if (match) index += match[0].length;
    return match;
  };
  // Inside brackets a newline is blank; at module level it ends the statement.
  const skip = () => take(depth ? /(?:\s|#[^\n]*)*/y : /(?:[ \t\r]|#[^\n]*)*/y);
  const sequence = (close) => {
    const items = [];
    const entries = [];
    depth += 1;
    for (skip(); source[index] !== close; skip()) {
      const item = value();
      skip();
      if (close === '}' && source[index] === ':' && typeof item === 'string') {
        index += 1;
        skip();
        entries.push([item, value()]);
        skip();
      } else items.push(item);
      if (source[index] === ',') index += 1;
      else if (source[index] !== close) fail();
    }
    index += 1;
    depth -= 1;
    // `{}` and `{key: value}` are dicts; `{value}` is a set.
    if (close !== '}' || (items.length && !entries.length)) return items;
    return items.length ? fail() : Object.fromEntries(entries);
  };
  const value = () => {
    const closer = { '[': ']', '(': ')', '{': '}' }[source[index]];
    if (closer) {
      index += 1;
      return sequence(closer);
    }
    const text = take(/"((?:\\.|[^"\\\n])*)"|'((?:\\.|[^'\\\n])*)'/y);
    if (text) return text[1] ?? text[2];
    const number = take(/\d+(?![\w.])/y);
    if (number) return Number(number[0]);
    const identifier = take(/[A-Za-z_]\w*/y);
    if (!identifier) fail();
    if (!/^(?:frozen)?set$/.test(identifier[0]) || source[index] !== '(') return null;
    index += 1;
    const [members = [], ...rest] = sequence(')');
    return rest.length || !Array.isArray(members) ? fail() : members;
  };
  try {
    const result = value();
    skip();
    return index < source.length && source[index] !== '\n'
      ? UNREADABLE_ALLOWLIST
      : { status: 'ok', value: result, span: [declaration.index, index] };
  } catch (_error) {
    return UNREADABLE_ALLOWLIST;
  }
};
// The `symbol` array of a JSON file as { symbol: { id: entry text } }: every
// entry is an object with its own string `id`.
const readJsonEntries = (source, symbol) => {
  let parsed;
  try {
    parsed = JSON.parse(source);
  } catch (_error) {
    return UNREADABLE_ALLOWLIST;
  }
  if (!parsed || typeof parsed !== 'object' || Array.isArray(parsed)) return UNREADABLE_ALLOWLIST;
  if (!Object.hasOwn(parsed, symbol)) return { status: 'missing', value: null };
  const entries = parsed[symbol];
  if (!Array.isArray(entries) || entries.some((entry) => typeof entry?.id !== 'string')) return UNREADABLE_ALLOWLIST;
  const byId = Object.fromEntries(entries.map((entry) => [entry.id, JSON.stringify(entry)]));
  if (Object.keys(byId).length !== entries.length) return UNREADABLE_ALLOWLIST;
  return { status: 'ok', value: { [symbol]: byId } };
};
const isStringList = (value) => Array.isArray(value) && value.every((item) => typeof item === 'string');
const isPlainObject = (value) => Boolean(value) && typeof value === 'object' && !Array.isArray(value);
const mapObjectValues = (value, convert) => Object.fromEntries(
  Object.entries(value).map(([key, item]) => [key, convert(item)])
);
// A JSON entry's list items as slots `field#index`: a further item is a new
// slot and a renamed one keeps its slot (R9: an owner added to a twin row is
// authority-only, a renamed owner is not).
const listSlots = (text) => Object.fromEntries(Object.entries(JSON.parse(text))
  .filter(([, item]) => Array.isArray(item))
  .flatMap(([field, items]) => items.map((_item, index) => [`${field}#${index}`, ''])));
// The literal of one registered allowlist in the shape its comparison reads:
// a path -> names map compares as a writer-map without categories, and a
// Python `(count, note)` value as its count.
const readAllowlistLiteral = (source, symbol, kind) => {
  if (source === null) return { status: 'absent', value: null };
  const read = isPythonAllowlistKind(kind) ? readPythonLiteral
    : JSON_ALLOWLIST_KINDS.has(kind) ? readJsonEntries : readJavaScriptLiteral;
  const { status, value } = read(source, symbol);
  if (status !== 'ok') return { status, value: null };
  if (kind === 'json-ids') return { status, value: mapObjectValues(value[symbol], listSlots) };
  if (kind === 'set-map') {
    if (!isPlainObject(value) || !Object.values(value).every(isStringList)) return UNREADABLE_ALLOWLIST;
    return { status, value: mapObjectValues(value, (names) => Object.fromEntries(names.map((name) => [name, '']))) };
  }
  if (kind === 'py-count-map' && isPlainObject(value)) {
    return { status, value: mapObjectValues(value, (item) => (Array.isArray(item) ? item[0] : item)) };
  }
  return { status, value };
};
const isCountMap = (value) => value && typeof value === 'object' && !Array.isArray(value)
  && Object.values(value).every((item) => Number.isInteger(item) && item >= 0);
const isWriterMap = (value) => value && typeof value === 'object' && !Array.isArray(value)
  && Object.values(value).every((inner) => inner && typeof inner === 'object' && !Array.isArray(inner)
    && Object.values(inner).every((item) => typeof item === 'string'));
// 'contraction' | 'unchanged' | 'expansion' | 'unreadable'
const classifyAllowlistDelta = (kind, before, after) => {
  if (kind === 'count') {
    if (!Number.isInteger(before) || !Number.isInteger(after)) return 'unreadable';
    return after < before ? 'contraction' : after === before ? 'unchanged' : 'expansion';
  }
  if (kind === 'set') {
    if (!isStringList(before) || !isStringList(after)) return 'unreadable';
    const baseSet = new Set(before);
    if (after.some((item) => !baseSet.has(item))) return 'expansion';
    return new Set(after).size < baseSet.size ? 'contraction' : 'unchanged';
  }
  if (kind === 'count-map') {
    if (!isCountMap(before) || !isCountMap(after)) return 'unreadable';
    let contracted = false;
    for (const [key, value] of Object.entries(after)) {
      if (!Object.hasOwn(before, key) || value > before[key]) return 'expansion';
      if (value < before[key]) contracted = true;
    }
    if (Object.keys(after).length < Object.keys(before).length) contracted = true;
    return contracted ? 'contraction' : 'unchanged';
  }
  if (kind === 'writer-map') {
    if (!isWriterMap(before) || !isWriterMap(after)) return 'unreadable';
    let contracted = false;
    for (const [path, inner] of Object.entries(after)) {
      if (!Object.hasOwn(before, path)) return 'expansion';
      for (const [name, category] of Object.entries(inner)) {
        if (!Object.hasOwn(before[path], name) || before[path][name] !== category) return 'expansion';
      }
      if (Object.keys(inner).length < Object.keys(before[path]).length) contracted = true;
    }
    if (Object.keys(after).length < Object.keys(before).length) contracted = true;
    return contracted ? 'contraction' : 'unchanged';
  }
  return 'unreadable';
};
const designRuleRegistryEvaluation = parseDesignRuleGuardRegistry(
  readRevisionFile(base, designRuleGuardRegistryPath)
);
const designRuleGuardRegistry = designRuleRegistryEvaluation.registry;
const designRuleRegistryErrors = designRuleRegistryEvaluation.errors;
// Every rule's text is authority (BA-D10): a rule whose text a runtime diff
// could change does not bind it. Checked on the candidate registry whenever it
// changes, so the change that registers or moves a rule names an authority path.
const candidateDesignRuleRulesPathErrors = (() => {
  if (!changed.has(designRuleGuardRegistryPath)) return [];
  const { registry } = parseDesignRuleGuardRegistry(readHeadFile(designRuleGuardRegistryPath));
  return (registry?.rules || []).flatMap((rule, index) => {
    const rulesPath = rule.rulesPath || registry.rulesPath;
    return authorityPaths.has(rulesPath) ? [] : [`rules[${index}] ${rule.id}: rulesPath ${rulesPath} is not an authority path`];
  });
})();
const designRuleAllowlistChanges = [];
if (designRuleGuardRegistry) {
  designRuleGuardRegistry.rules.forEach((rule) => {
    rule.allowlists.forEach(({ path, symbol, kind }) => {
      if (!changed.has(path)) return;
      const before = readAllowlistLiteral(readRevisionFile(base, path), symbol, kind);
      const after = readAllowlistLiteral(readHeadFile(path), symbol, kind);
      const subject = `${rule.id} ${path}#${symbol}`;
      let direction;
      if (before.status !== 'ok') direction = `unreadable in the trusted base (${before.status})`;
      else if (after.status === 'absent' || after.status === 'missing') direction = 'contraction';
      else if (after.status !== 'ok') direction = 'unreadable';
      else direction = classifyAllowlistDelta(DESIGN_RULE_ALLOWLIST_COMPARISONS[kind] || kind, before.value, after.value);
      designRuleAllowlistChanges.push({ subject, direction, kind });
    });
  });
}
const designRuleAllowlistExpansionChanges = designRuleAllowlistChanges
  .filter(({ direction }) => direction !== 'contraction' && direction !== 'unchanged');
const describeAllowlistExpansions = (changes) => changes
  .map(({ subject, direction }) => `${subject} (${direction})`)
  .sort();
const designRuleAllowlistExpansions = describeAllowlistExpansions(designRuleAllowlistExpansionChanges);
// design-rule.co-change: the expansions that share the diff with production
// code. Web runtime files count for every kind; a Python baseline (`py-*`)
// guards the Python package, so gbdraw/**/*.py outside gbdraw/web/ counts for
// it as well.
const pythonPackageChanged = [...changed.keys()].some((path) => /^gbdraw\/(?!web\/).*\.py$/.test(path));
const designRuleCoChangeExpansions = describeAllowlistExpansions(designRuleAllowlistExpansionChanges.filter(
  ({ kind }) => productionPaths.length || (isPythonAllowlistKind(kind) && pythonPackageChanged)
));
const designRuleAllowlistContractions = designRuleAllowlistChanges
  .filter(({ direction }) => direction === 'contraction')
  .map(({ subject }) => subject)
  .sort();
// design-rule.count-scope: a test that holds a registered count baseline
// counts with its own code, so a change to that code can lower the count and
// gain room under the ceiling. A diff with production code (as for co-change)
// may change such a test only inside its registered literals (of any kind),
// which stay readable literals; anything else in it, or its deletion, goes to
// a pull request without production code. The same holds for the counting
// code those tests share outside their files.
const COUNT_ALLOWLIST_KINDS = new Set(['count', 'count-map', 'py-count', 'py-count-map']);
const COUNTING_CODE_PATHS = new Set([
  'tests/web/helpers/count-baseline.mjs',
  'tests/web/types/tsconfig.json',
  'tests/web/types/web-globals.d.ts',
  'tools/report-web-owner-graph.mjs'
]);
const withoutAllowlistLiterals = (source, entries) => {
  const spans = [];
  for (const { symbol, kind } of entries) {
    const { status, span } = (isPythonAllowlistKind(kind) ? readPythonLiteral : readJavaScriptLiteral)(source, symbol);
    if (status !== 'ok') return null;
    spans.push({ symbol, span });
  }
  return spans.sort((left, right) => right.span[0] - left.span[0])
    .reduce((text, { symbol, span: [start, end] }) => `${text.slice(0, start)}<${symbol}>${text.slice(end)}`, source);
};
// A runtime move cannot be split from the test that names the moved file
// (#927), so a quoted path that names the old path of a rename in this diff
// may name the new path, in the same form: relative to the test, or the same
// trailing path segments.
let diffRenames = null;
const renamesInDiff = () => {
  diffRenames ||= parseDiffLines(runGit(repositoryRoot, [
    'diff', '-M', '--diff-filter=R', '--name-status', ...diffRefs, '--'
  ])).map(([, from, to]) => [from, to]);
  return diffRenames;
};
const followsRename = (testPath, before, after) => renamesInDiff().some(([from, to]) => {
  if (/^\.\.?\//.test(before)) {
    const directory = posix.dirname(testPath);
    return posix.join(directory, before) === from && posix.join(directory, after) === to;
  }
  if (from !== before && !from.endsWith(`/${before}`)) return false;
  const prefix = from.slice(0, from.length - before.length);
  return to.startsWith(prefix) && to.slice(prefix.length) === after;
});
const QUOTED_ON_ONE_LINE = /(['"`])([^'"`\n]*)\1/;
const differsOnlyByRenamedPaths = (testPath, before, after) => {
  const left = before.split(QUOTED_ON_ONE_LINE);
  const right = after.split(QUOTED_ON_ONE_LINE);
  return left.length === right.length && left.every((part, index) => (
    part === right[index] || (index % 3 === 2 && followsRename(testPath, part, right[index]))
  ));
};
const changedAllowlistEntriesByPath = new Map();
(designRuleGuardRegistry?.rules || []).forEach((rule) => rule.allowlists.forEach((entry) => {
  if (!changed.has(entry.path) || JSON_ALLOWLIST_KINDS.has(entry.kind)) return;
  changedAllowlistEntriesByPath.set(entry.path, [...(changedAllowlistEntriesByPath.get(entry.path) || []), entry]);
}));
const countScopeEntriesByPath = [
  ...[...changedAllowlistEntriesByPath].filter(([, entries]) => entries.some(({ kind }) => (
    COUNT_ALLOWLIST_KINDS.has(kind) && (productionPaths.length || (isPythonAllowlistKind(kind) && pythonPackageChanged))
  ))),
  ...(productionPaths.length ? [...COUNTING_CODE_PATHS].filter((path) => changed.has(path)).map((path) => [path, []]) : [])
];
const designRuleCountScopeChanges = countScopeEntriesByPath.flatMap(([path, entries]) => {
  const baseSource = readRevisionFile(base, path);
  if (baseSource === null) return [];
  const before = withoutAllowlistLiterals(baseSource, entries);
  // A literal unreadable in the trusted base is already a co-change expansion.
  if (before === null) return [];
  const headSource = readHeadFile(path);
  if (headSource === null) return [`${path} (deleted)`];
  const after = withoutAllowlistLiterals(headSource, entries);
  if (after === null) return [`${path} (registered literal missing or unreadable)`];
  return differsOnlyByRenamedPaths(path, before, after) ? [] : [path];
}).sort();
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

// Owner-graph report (implementation plan Phase B1): report only. The
// registry, when present, is read from the trusted base; the detectors run on
// the base and the head sources.
const readOwnerGraphRegistry = () => {
  const source = readRevisionFile(base, ownerGraphRegistryPath);
  if (source === null) return { registry: null, error: null };
  try {
    const parsed = JSON.parse(source);
    return parsed && typeof parsed === 'object' && !Array.isArray(parsed)
      ? { registry: parsed, error: null }
      : { registry: null, error: `${ownerGraphRegistryPath}: expected an object` };
  } catch (error) {
    return { registry: null, error: `${ownerGraphRegistryPath}: invalid JSON (${error.message})` };
  }
};
const ownerGraphRegistryEvaluation = readOwnerGraphRegistry();
const ownerGraphReport = (() => {
  try {
    const registry = ownerGraphRegistryEvaluation.registry;
    const headResults = detectWebOwnerGraph(allProductionSources, registry);
    const baseResults = detectWebOwnerGraph(baseProductionSources, registry);
    const subjectDelta = (id) => {
      const before = new Set(baseResults[id].subjects);
      const after = new Set(headResults[id].subjects);
      return {
        added: [...after].filter((subject) => !before.has(subject)).sort(),
        removed: [...before].filter((subject) => !after.has(subject)).sort()
      };
    };
    return {
      error: null,
      head: summarizeWebOwnerGraph(headResults),
      base: summarizeWebOwnerGraph(baseResults),
      deltas: Object.fromEntries(Object.keys(headResults).map((id) => [id, subjectDelta(id)]))
    };
  } catch (error) {
    return { error: error?.message || String(error), head: null, base: null, deltas: {} };
  }
})();

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

// JSON.parse keeps the last copy of a repeated key, so authority could hide an entry.
// Call only on source that JSON.parse accepted: every string token is then terminated.
const repeatedJsonKeys = (source) => {
  const repeated = new Set();
  const scopes = [];
  const colon = /\s*:/y;
  for (let index = 0; index < source.length; index += 1) {
    const character = source[index];
    if (character === '{') scopes.push(new Set());
    else if (character === '[') scopes.push(null);
    else if (character === '}' || character === ']') scopes.pop();
    else if (character === '"') {
      let end = index + 1;
      while (source[end] !== '"') end += source[end] === '\\' ? 2 : 1;
      const keys = scopes.at(-1);
      colon.lastIndex = end + 1;
      if (keys && colon.test(source)) {
        const key = JSON.parse(source.slice(index, end + 1));
        if (keys.has(key)) repeated.add(key);
        keys.add(key);
      }
      index = end;
    }
  }
  return [...repeated].sort();
};

const parsePolicy = (source, revision) => {
  if (source === null) throw new Error(`Missing ${policyPath} at ${revision}`);
  let policy;
  try {
    policy = JSON.parse(source);
  } catch (error) {
    throw new Error(`Cannot parse ${policyPath} at ${revision}: ${error.message}`);
  }
  const repeatedKeys = repeatedJsonKeys(source);
  if (repeatedKeys.length) {
    throw new Error(`${policyPath} at ${revision} repeats the key(s) ${repeatedKeys.join(', ')}`);
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
// A capability key that this (trusted) detector does not define authorizes nothing.
const retiredPolicyKeys = [];
const detectedCapabilityKeys = new Set(WEB_PRIVILEGED_CAPABILITY_KEYS);
const isRetiredPolicyKey = (section, name) => (
  section === 'allowedPrivilegedOwners' && !detectedCapabilityKeys.has(name)
);
let policyTopLevelKeysMatch = true;
// Every policy value outside the two privileged sections stays unchanged in a safe contraction.
let policyOtherSectionsMatch = true;
if (basePolicySource !== null && changed.has(policyPath)) {
  proposedPolicy = parsePolicy(readHeadFile(policyPath), head || 'working tree');
  const baseTopLevelKeys = Object.keys(policy).sort();
  const proposedTopLevelKeys = Object.keys(proposedPolicy).sort();
  policyTopLevelKeysMatch = JSON.stringify(baseTopLevelKeys) === JSON.stringify(proposedTopLevelKeys);
  policyOtherSectionsMatch = baseTopLevelKeys
    .filter((key) => !['allowedPrivilegedImporters', 'allowedPrivilegedOwners'].includes(key))
    .every((key) => JSON.stringify(policy[key]) === JSON.stringify(proposedPolicy[key]));
  for (const section of ['allowedPrivilegedImporters', 'allowedPrivilegedOwners']) {
    Object.keys(proposedPolicy[section]).forEach((name) => {
      if (!Object.hasOwn(policy[section], name)) addedPolicyKeys.push(`${section}.${name}`);
    });
    Object.entries(policy[section]).forEach(([name, allowedPaths]) => {
      if (!Object.hasOwn(proposedPolicy[section], name)) {
        (isRetiredPolicyKey(section, name) ? retiredPolicyKeys : missingPolicyKeys)
          .push(`${section}.${name}`);
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

// Vendored runtime changes (WEB_CHANGE_POLICY.md, Gate FAIL). Only the trusted base listing
// authorizes a change; a candidate listing is validated as inert data and authorizes nothing.
const vendorRoot = 'gbdraw/web/vendor/';
const vendorChangeKindByStatus = Object.freeze({ A: 'add', M: 'replace', D: 'remove' });
const vendorChangeKinds = Object.values(vendorChangeKindByStatus);
const vendorListingKeys = Object.freeze(['change', 'package', 'path', 'sha256', 'source', 'version']);
const validateVendorListing = (candidatePolicy) => {
  const listing = candidatePolicy.allowedVendorChanges;
  if (listing === undefined) return { entries: [], errors: [] };
  if (!Array.isArray(listing)) {
    return { entries: [], errors: ['allowedVendorChanges must be an array'] };
  }
  const errors = [];
  const listedPaths = new Set();
  const nonEmpty = (value) => typeof value === 'string' && value.trim() !== '';
  listing.forEach((entry, index) => {
    const label = `allowedVendorChanges[${index}]`;
    if (!entry || typeof entry !== 'object' || Array.isArray(entry)) {
      errors.push(`${label} must be an object`);
      return;
    }
    if (Object.keys(entry).sort().join() !== vendorListingKeys.join()) {
      errors.push(`${label} must have exactly the keys ${vendorListingKeys.join(', ')}`);
    }
    const { path } = entry;
    if (
      typeof path !== 'string'
      || !path.startsWith(vendorRoot)
      || path.endsWith('/')
      || path.includes('\\')
      || posix.normalize(path) !== path
    ) {
      errors.push(`${label}.path must be a normalized path under ${vendorRoot}`);
    } else if (listedPaths.has(path)) {
      errors.push(`${label} lists ${path} twice`);
    } else {
      listedPaths.add(path);
    }
    if (!vendorChangeKinds.includes(entry.change)) {
      errors.push(`${label}.change must be one of ${vendorChangeKinds.join(', ')}`);
    }
    if (typeof entry.sha256 !== 'string' || !/^[0-9a-f]{64}$/.test(entry.sha256)) {
      errors.push(`${label}.sha256 must be 64 lowercase hexadecimal digits`);
    }
    ['package', 'version'].forEach((field) => {
      if (!nonEmpty(entry[field])) errors.push(`${label}.${field} must be a non-empty string`);
    });
    if (typeof entry.source !== 'string' || !/^https:\/\/\S+$/.test(entry.source)) {
      errors.push(`${label}.source must be an https URL`);
    }
  });
  return { entries: errors.length ? [] : listing, errors };
};
const readRevisionBytes = (revision, path) => {
  try {
    return execFileSync(
      'git',
      ['-C', repositoryRoot, 'cat-file', 'blob', `${revision}:${path}`],
      { maxBuffer: 512 * 1024 * 1024, stdio: ['ignore', 'pipe', 'ignore'] }
    );
  } catch (_error) {
    return null;
  }
};
const readHeadBytes = (path) => {
  if (head) return readRevisionBytes(head, path);
  const absolutePath = join(repositoryRoot, path);
  return existsSync(absolutePath) ? readFileSync(absolutePath) : null;
};
const sha256OfBytes = (bytes) => (
  bytes === null ? 'unreadable' : createHash('sha256').update(bytes).digest('hex')
);
// The bootstrap fallback above reads the head policy; the listing never comes from the head.
const trustedVendorListing = basePolicySource === null
  ? { entries: [], errors: [] }
  : validateVendorListing(policy);
const candidateVendorListing = changed.has(policyPath)
  ? validateVendorListing(proposedPolicy || policy)
  : null;
const vendorEntryKey = (entry) => JSON.stringify(entry);
const trustedVendorEntryKeys = new Set(trustedVendorListing.entries.map(vendorEntryKey));
const candidateVendorEntryKeys = new Set((candidateVendorListing?.entries || []).map(vendorEntryKey));
const addedVendorEntries = (candidateVendorListing?.entries || [])
  .filter((entry) => !trustedVendorEntryKeys.has(vendorEntryKey(entry)));
const removedVendorEntries = candidateVendorListing && !candidateVendorListing.errors.length
  ? trustedVendorListing.entries.filter((entry) => !candidateVendorEntryKeys.has(vendorEntryKey(entry)))
  : [];
const describeVendorEntry = (entry) => (
  `${entry.path}: ${entry.change} sha256=${entry.sha256} package=${entry.package} `
  + `version=${entry.version} source=${entry.source}`
);
// A new entry that cannot apply to the base tree fails in the pull request that adds it.
const impossibleVendorEntries = addedVendorEntries.flatMap((entry) => {
  const label = `${entry.path}: ${entry.change}`;
  const baseBytes = readRevisionBytes(base, entry.path);
  if (entry.change === 'add') {
    if (baseBytes !== null) return [`${label} for a path the base already has`];
    // A tree (directory) at the path has no blob bytes but still occupies it.
    const occupied = spawnSync(
      'git',
      ['-C', repositoryRoot, 'cat-file', '-e', `${base}:${entry.path}`],
      { stdio: 'ignore' }
    ).status === 0;
    return occupied ? [`${label} for a path that is not a file in the base`] : [];
  }
  if (baseBytes === null) return [`${label} for a path the base lacks`];
  const baseSha256 = sha256OfBytes(baseBytes);
  if (entry.change === 'remove' && baseSha256 !== entry.sha256) {
    return [`${label} lists sha256 ${entry.sha256}, but the base bytes hash to ${baseSha256}`];
  }
  if (entry.change === 'replace' && baseSha256 === entry.sha256) {
    return [`${label} lists the bytes the base already has`];
  }
  return [];
});
// An entry pins the bytes the change leaves (add, replace) or removes (remove).
// The diff is two-dot. When the head does not contain the base and leaves a path as their
// merge base had it, the path changed only in the base: the branch is behind, not wrong.
const vendorMergeBase = (() => {
  if (!head || !changedVendorPaths.length) return null;
  const ancestry = spawnSync(
    'git',
    ['-C', repositoryRoot, 'merge-base', '--is-ancestor', base, head],
    { stdio: 'ignore' }
  );
  if (ancestry.status !== 1) return null;
  const mergeBase = spawnSync(
    'git',
    ['-C', repositoryRoot, 'merge-base', base, head],
    { encoding: 'utf8', stdio: ['ignore', 'pipe', 'ignore'] }
  );
  return mergeBase.status === 0 ? mergeBase.stdout.trim() : null;
})();
const vendorChanges = changedVendorPaths.map((path) => {
  const status = changed.get(path);
  const change = vendorChangeKindByStatus[status] || `unsupported status ${status}`;
  const entry = trustedVendorListing.entries.find((candidate) => candidate.path === path);
  if (!entry) return { path, change, reason: 'not listed' };
  const pinned = sha256OfBytes(
    change === 'remove' ? readRevisionBytes(base, path) : readHeadBytes(path)
  );
  if (entry.change === change && pinned === entry.sha256) return { path, change, entry, reason: '' };
  if (
    vendorMergeBase
    && sha256OfBytes(readHeadBytes(path)) === sha256OfBytes(readRevisionBytes(vendorMergeBase, path))
  ) {
    return {
      path,
      change,
      reason: 'changed only in the base since this branch forked; '
        + 'the branch is behind its base: update it from dev'
    };
  }
  if (entry.change !== change) {
    return { path, change, reason: `listed as ${entry.change}, but the change is ${change}` };
  }
  return { path, change, reason: `sha256 ${pinned} does not match the listed ${entry.sha256}` };
});
const unlistedVendorChanges = vendorChanges.filter(({ reason }) => reason);
const listedVendorChanges = vendorChanges.filter(({ reason }) => !reason);

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
  productContractAuthorityPath,
  productImpactMapPath,
  productDecisionAuthorityPath
].filter((path) => changed.has(path));
const productContractAuthorityChanged = changed.has(productContractAuthorityPath);
const productContractCompanionPaths = changedPaths
  .filter((path) => path !== productContractAuthorityPath);
// WEB_CHANGE_POLICY.md "Static Product Contract co-change": the exact Contract,
// modified in place, is the only changed guard or authority, and no checker
// implementation changes. Human Review, not the checker, verifies its receipt.
const staticProductContractCoChange = Boolean(
  changed.get(productContractAuthorityPath) === 'M'
  && productContractCompanionPaths.length
  && changedGuards.length === 1
  && changedGuards[0] === productContractAuthorityPath
  && changedAuthorities.length === 1
  && changedAuthorities[0] === productContractAuthorityPath
  && !changedCheckerImplementations.length
);
const candidateAuthorityCompanionPaths = productContractAuthorityChanged
  ? (staticProductContractCoChange ? [] : productContractCompanionPaths)
  : changedProductImpactAuthorityPaths.length
    ? changedPaths.filter((path) => !narrowAuthorityBundlePaths.has(path))
    : [];
if (candidateAuthorityCompanionPaths.length) {
  productImpactCandidateErrors.push(
    (productContractAuthorityChanged
      ? 'Product Contract authority changes must be isolated from other changed paths: '
      : 'Product Impact authority changes must use only the narrow inert authority bundle: ')
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
  && !retiredPolicyKeys.length
  && !addedPolicyKeys.length
  && policyTopLevelKeysMatch
  && policyOtherSectionsMatch
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
dependencyChanges.push(...vendorChanges.map(({ path, change, entry }) => (
  `vendored runtime file changed: ${path}`
  + (entry ? ` (listed ${change}: ${entry.package} ${entry.version})` : '')
)));

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
if (addedVendorEntries.length || removedVendorEntries.length) {
  reviewReasons.push(
    'Governance and authority: vendored runtime listing changed '
    + `(added=${addedVendorEntries.length}, removed=${removedVendorEntries.length}; `
    + 'see the vendored runtime listing sections)'
  );
}
listedVendorChanges.forEach(({ path, change, entry }) => {
  reviewReasons.push(
    `Vendored runtime: listed ${change} of ${path} (${entry.package} ${entry.version})`
  );
});
if (designRuleAllowlistExpansions.length) {
  reviewReasons.push(
    'Governance and authority: registered design-rule allowlists expanded '
    + `(${designRuleAllowlistExpansions.join(', ')})`
  );
}
if (designRuleRegistryErrors.length) {
  reviewReasons.push('Governance and authority: required authority input is malformed');
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
if (staticProductContractCoChange) {
  reviewReasons.push(
    'Governance and authority: static Product Contract co-change requires human review '
    + 'that each changed record serializes the explicit Product Decision Owner receipt and '
    + 'that the companion runtime, tests, and documentation realize exactly that outcome; '
    + 'the Gate does not verify Contract prose'
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
if (unlistedVendorChanges.length) {
  blockingViolations.push(
    'changes under gbdraw/web/vendor/ are not allowed unless the trusted base lists them '
    + `in ${policyPath} allowedVendorChanges (`
    + `${unlistedVendorChanges.map(({ path, reason }) => `${path}: ${reason}`).join('; ')})`
  );
}
blockingViolations.push(...trustedVendorListing.errors.map((error) => (
  `trusted vendored runtime listing: ${error}`
)));
blockingViolations.push(...[
  ...(candidateVendorListing?.errors || []),
  ...impossibleVendorEntries
].map((error) => `candidate vendored runtime listing: ${error}`));
if (
  productionPaths.length
  && changedGuards.length
  && !pureSafePolicyContraction
  && !staticProductContractCoChange
) {
  blockingViolations.push('production runtime files and Web guard/CI files changed together');
}
if (changedCheckerImplementations.length && changedAuthorities.length) {
  blockingViolations.push(
    'Web checker/source parser and authority policy/workflow files changed together'
  );
}
blockingViolations.push(...designRuleRegistryErrors.map((error) => (
  `trusted design-rule guard registry: ${error}`
)));
blockingViolations.push(...candidateDesignRuleRulesPathErrors.map((error) => (
  `candidate design-rule guard registry: ${error}`
)));
if (designRuleCoChangeExpansions.length) {
  blockingViolations.push(
    'design-rule.co-change: a registered design-rule allowlist expanded alongside production runtime files '
    + `(${designRuleCoChangeExpansions.join(', ')})`
  );
}
if (designRuleCountScopeChanges.length) {
  blockingViolations.push(
    'design-rule.count-scope: the test or counting code of a registered count baseline changed '
    + `alongside production runtime files (${designRuleCountScopeChanges.join(', ')})`
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
    '  - `DURABLE_AUTHORITY_REQUIRED`: Product Decision Owner response -> Codex prepares a narrow authority-only PR -> merge -> rebase implementation; a static Product Contract record may instead serialize the receipt in a reviewed co-change with its implementation.',
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
  `- Candidate authority validation: ${candidateProductImpactMapValid && candidateProductDecisionAuthorityValid ? 'VALID (inert map and decision data only)' : 'INVALID'}`,
  `- Candidate authority separation: ${candidateAuthorityCompanionPaths.length ? 'INVALID' : 'VALID'}`,
  `- Static Product Contract: ${!productContractAuthorityChanged ? 'unchanged'
    : staticProductContractCoChange ? 'co-change with its implementation; human Review verifies the receipt, and the text is not machine authority'
      : candidateAuthorityCompanionPaths.length ? 'INVALID companion changes' : 'isolated change for human Review'}`,
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
    '- Runtime effect: validation-only future preauthorization; candidate data does not alter this head runtime admission. A static Product Contract change is reviewed by a human and never becomes machine authority.'
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
  '## Removed retired privileged allowlist keys',
  '',
  ...list(retiredPolicyKeys),
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
  '## Added vendored runtime listing entries',
  '',
  ...list(addedVendorEntries.map(describeVendorEntry)),
  '',
  '## Removed vendored runtime listing entries',
  '',
  ...list(removedVendorEntries.map(describeVendorEntry)),
  '',
  '## Dependency changes',
  '',
  ...list(dependencyChanges),
  '',
  '## Guard files touched',
  '',
  ...list(changedGuards),
  '',
  '## Design-rule guard allowlists',
  '',
  `- Registry: ${designRuleGuardRegistry
    ? `${designRuleGuardRegistry.rules.length} rule(s) from the trusted base`
    : designRuleRegistryErrors.length ? 'malformed' : 'absent (inert)'}`,
  ...designRuleAllowlistContractions.map((subject) => `- contraction: ${subject}`),
  ...designRuleAllowlistExpansions.map((subject) => `- expansion: ${subject}`),
  '## Owner graph (report)',
  '',
  ...(ownerGraphReport.error
    ? [`- Detector error: ${ownerGraphReport.error}`]
    : [
      `- Registry: ${ownerGraphRegistryEvaluation.registry
        ? `${ownerGraphRegistryPath} from the trusted base`
        : ownerGraphRegistryEvaluation.error ? `malformed (${ownerGraphRegistryEvaluation.error})` : 'defaults (registry absent)'}`,
      ...['injectionEdges', 'forwardClosures', 'stateBackdoors', 'wholeObjectPorts', 'triggerSites', 'triggerModules']
        .map((key) => `- ${key}: ${ownerGraphReport.base[key]} -> ${ownerGraphReport.head[key]}`),
      ...Object.keys(ownerGraphReport.head.projectionShapes).map((domain) => (
        `- shapes:${domain}: ${ownerGraphReport.base.projectionShapes[domain] ?? 0} -> ${ownerGraphReport.head.projectionShapes[domain]}`
      )),
      ...Object.entries(ownerGraphReport.deltas).flatMap(([id, { added, removed }]) => [
        ...added.map((subject) => `- ${id}: NEW ${subject}`),
        ...removed.map((subject) => `- ${id}: REMOVED ${subject}`)
      ])
    ]),
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
      reviewReasons: /tools\/web-(?:change-policy|architecture-(?:rules|violations)|design-rule-guards)\.json/.test(
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
