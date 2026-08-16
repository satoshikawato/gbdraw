import { createRequire } from 'node:module';

const require = createRequire(import.meta.url);
const { parse } = require('acorn');

const LARGE_OWNER_NAMES = new Set([
  'biologicalFeatures',
  'extractedFeatures',
  'featureCatalog',
  'features',
  'orthogroups',
  'results',
  'sequences',
  'svg'
]);
const DOM_MUTATOR_NAMES = new Set([
  'after',
  'append',
  'appendChild',
  'before',
  'insertAdjacentElement',
  'insertAdjacentHTML',
  'insertAdjacentText',
  'insertBefore',
  'prepend',
  'remove',
  'removeAttribute',
  'removeChild',
  'replaceChild',
  'replaceChildren',
  'replaceWith',
  'setAttribute'
]);
const SERIALIZED_PROPERTY_NAMES = new Set(['innerHTML', 'outerHTML', 'textContent']);

const isNode = (value) => Boolean(value && typeof value === 'object' && typeof value.type === 'string');

const childrenOf = (node) => {
  const children = [];
  Object.entries(node || {}).forEach(([key, value]) => {
    if (key === 'parent' || key === 'loc' || key === 'start' || key === 'end') return;
    if (isNode(value)) children.push(value);
    else if (Array.isArray(value)) value.filter(isNode).forEach((child) => children.push(child));
  });
  return children;
};

const walk = (root, visit) => {
  const pending = [{ node: root, parent: null }];
  while (pending.length > 0) {
    const current = pending.pop();
    visit(current.node, current.parent);
    const children = childrenOf(current.node);
    for (let index = children.length - 1; index >= 0; index -= 1) {
      pending.push({ node: children[index], parent: current.node });
    }
  }
};

const memberName = (node) => {
  if (node?.type !== 'MemberExpression') return null;
  if (!node.computed && node.property.type === 'Identifier') return node.property.name;
  if (node.computed && node.property.type === 'Literal') return String(node.property.value);
  return null;
};

const memberPath = (node) => {
  if (!node) return '';
  if (node.type === 'ChainExpression') return memberPath(node.expression);
  if (node.type === 'Identifier') return node.name;
  if (node.type !== 'MemberExpression') return '';
  const owner = memberPath(node.object);
  const property = memberName(node);
  return owner && property ? `${owner}.${property}` : '';
};

const calleeName = (node) => {
  if (node?.type === 'Identifier') return node.name;
  return memberPath(node);
};

const literalTrueCallback = (node) => {
  if (!node || !['ArrowFunctionExpression', 'FunctionExpression'].includes(node.type)) return false;
  if (node.body.type === 'Literal') return node.body.value === true;
  if (node.body.type !== 'BlockStatement' || node.body.body.length !== 1) return false;
  const statement = node.body.body[0];
  return statement.type === 'ReturnStatement'
    && statement.argument?.type === 'Literal'
    && statement.argument.value === true;
};

const propertyValue = (objectExpression, name) => {
  if (objectExpression?.type !== 'ObjectExpression') return null;
  const property = objectExpression.properties.find((candidate) => (
    candidate.type === 'Property'
    && !candidate.computed
    && (candidate.key.name === name || candidate.key.value === name)
  ));
  return property?.value || null;
};

const addViolation = (violations, code, node, message) => {
  violations.push({ code, message, start: Number(node?.start ?? -1) });
};

const collectAliases = (ast) => {
  const declarations = [];
  walk(ast, (node) => {
    if (node.type === 'VariableDeclarator') declarations.push(node);
  });
  const largeOwners = new Set();
  const resultOwners = new Set();
  const journalOwners = new Set(['mutation']);
  const serializers = new Map();
  let changed = true;

  const isLargeOwner = (node) => {
    if (!node) return false;
    if (node.type === 'Identifier') return largeOwners.has(node.name);
    if (node.type === 'ChainExpression') return isLargeOwner(node.expression);
    if (node.type !== 'MemberExpression') return false;
    return LARGE_OWNER_NAMES.has(memberName(node));
  };
  const isResultOwner = (node) => {
    if (!node) return false;
    if (node.type === 'Identifier') return resultOwners.has(node.name);
    if (node.type === 'MemberExpression' && isResultOwner(node.object)) return true;
    const path = memberPath(node);
    return path.endsWith('results.value') || path === 'results';
  };
  const serializerKind = (node) => {
    const name = calleeName(node);
    if (!name) return null;
    if (serializers.has(name)) return serializers.get(name);
    if (name === 'JSON.stringify') return 'serialization';
    if (name === 'structuredClone') return 'structured-clone';
    if (/^(?:cloneJson(?:Value|Data)|deepClone)$/.test(name.split('.').at(-1))) return 'clone';
    return null;
  };
  const isJournalOwner = (node) => {
    if (!node) return false;
    if (node.type === 'Identifier') return journalOwners.has(node.name);
    const path = memberPath(node);
    return path === 'mutation' || path.endsWith('.mutation');
  };

  walk(ast, (node) => {
    if (!['ArrowFunctionExpression', 'FunctionExpression', 'FunctionDeclaration'].includes(node.type)) return;
    node.params.forEach((parameter) => {
      if (parameter.type !== 'ObjectPattern') return;
      parameter.properties.forEach((property) => {
        const sourceName = property.key?.name || property.key?.value;
        const localName = property.value?.name;
        if (sourceName === 'mutation' && localName) journalOwners.add(localName);
      });
    });
  });

  while (changed) {
    changed = false;
    declarations.forEach((declaration) => {
      const { id, init } = declaration;
      if (id.type === 'Identifier') {
        if (isLargeOwner(init) && !largeOwners.has(id.name)) {
          largeOwners.add(id.name);
          changed = true;
        }
        if (isResultOwner(init) && !resultOwners.has(id.name)) {
          resultOwners.add(id.name);
          changed = true;
        }
        if (
          (isJournalOwner(init) || (
            init?.type === 'CallExpression'
            && calleeName(init.callee).split('.').at(-1) === 'createDomMutationJournal'
          ))
          && !journalOwners.has(id.name)
        ) {
          journalOwners.add(id.name);
          changed = true;
        }
        const kind = serializerKind(init);
        if (kind && serializers.get(id.name) !== kind) {
          serializers.set(id.name, kind);
          changed = true;
        }
      }
      if (id.type === 'ObjectPattern' && init) {
        id.properties.forEach((property) => {
          const sourceName = property.key?.name || property.key?.value;
          const localName = property.value?.name;
          if (LARGE_OWNER_NAMES.has(sourceName) && localName && !largeOwners.has(localName)) {
            largeOwners.add(localName);
            changed = true;
          }
        });
      }
    });
  }

  const containsLargeOwner = (node) => {
    let found = false;
    if (!node) return false;
    walk(node, (candidate) => {
      if (!found && isLargeOwner(candidate)) found = true;
    });
    return found;
  };

  return {
    containsLargeOwner,
    isResultOwner,
    isJournalOwner,
    serializerKind
  };
};

const analyzeLargeOwnerCopies = (ast, aliases, violations) => {
  walk(ast, (node) => {
    if (node.type === 'CallExpression') {
      const kind = aliases.serializerKind(node.callee);
      const hasLargeArgument = node.arguments.some((argument) => aliases.containsLargeOwner(argument));
      if (kind && hasLargeArgument) {
        const code = kind === 'serialization'
          ? 'large-owner-serialization'
          : kind === 'structured-clone'
            ? 'large-owner-structured-clone'
            : 'large-owner-clone';
        addViolation(violations, code, node, 'A validated large owner is copied or serialized.');
      }
      if (
        ['Array.from', 'Object.assign'].includes(calleeName(node.callee))
        && hasLargeArgument
      ) {
        addViolation(violations, 'large-owner-clone', node, 'A validated large owner is copied.');
      }
    }
    if (
      node.type === 'SpreadElement'
      && aliases.containsLargeOwner(node.argument)
    ) {
      addViolation(violations, 'large-owner-clone', node, 'A validated large owner is spread into a copy.');
    }
  });
};

const analyzePlaceholderCommits = (ast, violations) => {
  walk(ast, (node) => {
    if (node.type !== 'CallExpression') return;
    const name = calleeName(node.callee).split('.').at(-1) || '';
    if (!/^commit(?:DomEdit|[A-Za-z]*Mutation)$/.test(name)) return;
    const callback = node.arguments[0]?.type === 'ObjectExpression'
      ? propertyValue(node.arguments[0], 'mutate')
      : node.arguments.find((argument) => (
          argument.type === 'ArrowFunctionExpression' || argument.type === 'FunctionExpression'
        ));
    if (literalTrueCallback(callback)) {
      addViolation(
        violations,
        'placeholder-mutation',
        callback,
        'A commit reports success without owning a mutation.'
      );
    }
  });
};

const analyzeJournalDomWrites = (ast, aliases, violations) => {
  walk(ast, (node) => {
    if (node.type === 'CallExpression' && node.callee.type === 'MemberExpression') {
      const method = memberName(node.callee);
      if (
        DOM_MUTATOR_NAMES.has(method)
        && !aliases.isJournalOwner(node.callee.object)
      ) {
        addViolation(
          violations,
          'unjournaled-dom-mutation',
          node,
          `DOM mutation ${method}() bypasses the transaction journal.`
        );
      }
    }
    if (node.type !== 'AssignmentExpression' || node.left.type !== 'MemberExpression') return;
    const property = memberName(node.left);
    const ownerProperty = memberName(node.left.object);
    if (SERIALIZED_PROPERTY_NAMES.has(property) || ownerProperty === 'style') {
      addViolation(
        violations,
        'unjournaled-dom-mutation',
        node,
        `DOM property ${property} bypasses the transaction journal.`
      );
    }
  });
};

const analyzeDirectResultOwnership = (ast, aliases, violations) => {
  walk(ast, (node) => {
    if (node.type === 'AssignmentExpression') {
      const leftOwner = node.left.type === 'MemberExpression' ? node.left.object : null;
      if (aliases.isResultOwner(leftOwner)) {
        addViolation(
          violations,
          'direct-result-write',
          node,
          'Result content is replaced outside the preview runtime owner.'
        );
      }
    }
    if (node.type === 'CallExpression') {
      const name = calleeName(node.callee).split('.').at(-1);
      if (name === 'serializeCleanSvg') {
        addViolation(
          violations,
          'direct-result-serialization',
          node,
          'Mounted SVG serialization is owned by preview-runtime.js.'
        );
      }
    }
  });
};

const analyzeAsyncTokenProtocol = (ast, violations) => {
  walk(ast, (node) => {
    if (!['ArrowFunctionExpression', 'FunctionExpression', 'FunctionDeclaration'].includes(node.type)) return;
    if (!node.async) return;
    const awaits = [];
    const captures = [];
    const validations = [];
    const commits = [];
    walk(node.body, (candidate) => {
      if (candidate.type === 'AwaitExpression') awaits.push(candidate);
      if (candidate.type !== 'CallExpression') return;
      const name = calleeName(candidate.callee).split('.').at(-1);
      if (name === 'captureDomEditToken') captures.push(candidate);
      if (name === 'isDomEditTokenCurrent' || name === 'targetIsCurrent') validations.push(candidate);
      if (name === 'commitDomEdit') commits.push(candidate);
    });
    const awaitedCommits = commits.filter((commit) => awaits.some((awaitNode) => awaitNode.end < commit.start));
    if (awaitedCommits.length === 0) return;
    const firstAwait = [...awaits].sort((left, right) => left.start - right.start)[0];
    if (!captures.some((capture) => capture.end < firstAwait.start)) {
      addViolation(
        violations,
        'async-target-without-token',
        node,
        'An async DOM commit does not capture its Result/SVG/document token before await.'
      );
      return;
    }
    awaits.forEach((awaitNode, index) => {
      const laterCommitStarts = awaitedCommits
        .filter((commit) => commit.start > awaitNode.end)
        .map((commit) => commit.start);
      const nextBoundary = awaits[index + 1]?.start
        ?? Math.min(...laterCommitStarts);
      if (!validations.some((validation) => validation.start > awaitNode.end && validation.start < nextBoundary)) {
        addViolation(
          violations,
          'async-target-without-revalidation',
          awaitNode,
          'An async DOM commit is not revalidated after await.'
        );
      }
    });
    awaitedCommits.forEach((commit) => {
      const options = commit.arguments[0];
      if (!propertyValue(options, 'targetToken')) {
        addViolation(
          violations,
          'async-commit-without-token',
          commit,
          'An async DOM commit is not pinned to its captured token.'
        );
      }
    });
  });
};

const analyzeSyntheticMetrics = (ast, violations) => {
  walk(ast, (node) => {
    if (node.type === 'CallExpression' && calleeName(node.callee).endsWith('recordStructuralMetric')) {
      const metricName = node.arguments[0]?.value;
      const metricValue = node.arguments[1]?.value;
      if (
        metricValue === 0
        && /(?:Clone|Serialization|UnusedBuild)Count$/.test(String(metricName || ''))
      ) {
        addViolation(
          violations,
          'synthetic-zero-metric',
          node,
          'A hard-coded zero does not observe a production ownership boundary.'
        );
      }
    }
    if (node.type !== 'Property' || node.value?.type !== 'Literal' || node.value.value !== 0) return;
    const name = node.key?.name || node.key?.value;
    if (/(?:Clone|Serialization|UnusedBuild)Count$/.test(String(name || ''))) {
      addViolation(
        violations,
        'synthetic-zero-metric',
        node,
        'A hard-coded zero does not observe a production ownership boundary.'
      );
    }
  });
};

export const analyzeWebRefactorSource = (source, {
  asyncTokenProtocol = false,
  directResultOwnership = false,
  journalDomWrites = false,
  largeOwnerCopies = false,
  placeholderCommits = false,
  syntheticMetrics = false
} = {}) => {
  const ast = parse(String(source || ''), {
    allowHashBang: true,
    ecmaVersion: 'latest',
    sourceType: 'module'
  });
  const aliases = collectAliases(ast);
  const violations = [];
  if (largeOwnerCopies) analyzeLargeOwnerCopies(ast, aliases, violations);
  if (placeholderCommits) analyzePlaceholderCommits(ast, violations);
  if (journalDomWrites) analyzeJournalDomWrites(ast, aliases, violations);
  if (directResultOwnership) analyzeDirectResultOwnership(ast, aliases, violations);
  if (asyncTokenProtocol) analyzeAsyncTokenProtocol(ast, violations);
  if (syntheticMetrics) analyzeSyntheticMetrics(ast, violations);
  return violations.sort((left, right) => left.start - right.start || left.code.localeCompare(right.code));
};
