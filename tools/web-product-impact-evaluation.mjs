const SCHEMA_VERSION = 1;
const MAP_FIELDS = Object.freeze(['concerns', 'ruleCoverage', 'schemaVersion']);
const RULE_COVERAGE_FIELDS = Object.freeze(['coverage', 'enforcement', 'ruleKey']);
const CONCERN_FIELDS = Object.freeze([
  'contracts',
  'effects',
  'journeys',
  'key',
  'layer',
  'options',
  'resolution',
  'scenarioRevision',
  'title'
]);
const JOURNEY_FIELDS = Object.freeze([
  'actor',
  'checkpoints',
  'context',
  'goal',
  'id',
  'steps'
]);
const CHECKPOINT_FIELDS = Object.freeze(['id', 'label']);
const EFFECT_FIELDS = Object.freeze(['checkpointRef', 'id', 'statement']);
const OPTION_FIELDS = Object.freeze(['effectRefs', 'id', 'requirements', 'summary']);
const REQUIREMENT_FIELDS = Object.freeze([
  'anyOf',
  'category',
  'checkpointRefs',
  'id'
]);
const CANDIDATE_FIELDS = Object.freeze(['ruleKey', 'subject', 'subjectCategory']);
const CONTRACT_FIELDS = Object.freeze([
  'checkpointRefs',
  'execution',
  'ref',
  'residualRisk',
  'sensitivity'
]);
const DECISION_AUTHORITY_FIELDS = Object.freeze([
  'decisions',
  'maintainerLogins',
  'schemaVersion'
]);
const DURABLE_DECISION_FIELDS = Object.freeze([
  'acceptanceCheckpointRefs',
  'acceptedImpactClass',
  'acknowledgedRetiredRequirementRefs',
  'concernKey',
  'id',
  'rationale',
  'scenarioRevision',
  'selectedOptionId'
]);
const SUBJECT_DELTA_FIELDS = Object.freeze([
  'addedSubjects',
  'afterSubjects',
  'beforeSubjects',
  'changed',
  'detector',
  'kind',
  'removedSubjects',
  'ruleKey',
  'subjectCategory'
]);
const RULE_FACT_FIELDS = Object.freeze([
  'afterSubjects',
  'beforeSubjects',
  'ruleKey',
  'subjectCategory'
]);
const EVALUATION_INPUT_FIELDS = Object.freeze([
  'changedPaths',
  'currentDecisions',
  'decisionAuthority',
  'map',
  'ruleFacts',
  'subjectDeltas'
]);

const COVERAGE_LEVELS = Object.freeze(['complete', 'partial']);
const PRODUCT_ENFORCEMENTS = Object.freeze(['hard', 'report-only']);
const CONCERN_LAYERS = Object.freeze([
  'DOMAIN',
  'INTERACTION',
  'OUTPUT',
  'PERSISTENCE',
  'STATE'
]);
const REQUIREMENT_CATEGORIES = Object.freeze(['AFFORDANCE', 'COMPATIBILITY', 'SEMANTIC']);
const RESOLUTION_KINDS = Object.freeze([
  'authority-covered',
  'decision-required',
  'domain-derived'
]);
const CONTRACT_SENSITIVITIES = Object.freeze([
  'BASE_FAIL',
  'DIRECT_ASSERTION',
  'MANUAL_ACCEPTANCE',
  'MUTATION_FAIL'
]);
const CONTRACT_EXECUTIONS = Object.freeze(['DEV_STAGING', 'MANUAL', 'PR_GATE']);
const ACCEPTED_IMPACT_CLASSES = Object.freeze([
  'AFFORDANCE_PRESERVED',
  'PRODUCT_CHANGE',
  'RETIREMENT'
]);
const RULE_KEY_PATTERN = /^[a-z][a-z0-9]*(?:[.-][a-z0-9]+)*$/;
const LOWER_ID_PATTERN = /^[a-z][a-z0-9]*(?:-[a-z0-9]+)*$/;
const DECISION_ID_PATTERN = /^BD-[0-9]{3,}$/;
const LOGIN_PATTERN = /^[a-z0-9](?:[a-z0-9-]{0,37}[a-z0-9])?$/;
const REPOSITORY_PATH_PATTERN = /^[A-Za-z0-9_.-]+(?:\/[A-Za-z0-9_.-]+)*$/;

const compareText = (left, right) => left < right ? -1 : left > right ? 1 : 0;
const isRecord = (value) => value !== null
  && typeof value === 'object'
  && !Array.isArray(value)
  && [Object.prototype, null].includes(Object.getPrototypeOf(value));
const stableUnique = (values) => [...new Set(values)].sort(compareText);
const arraysEqual = (left, right) => left.length === right.length
  && left.every((value, index) => value === right[index]);
const stableDifference = (left, right) => {
  const rightSet = new Set(right);
  return left.filter((value) => !rightSet.has(value));
};
const deepFreeze = (value) => {
  if (value && typeof value === 'object' && !Object.isFrozen(value)) {
    Object.values(value).forEach(deepFreeze);
    Object.freeze(value);
  }
  return value;
};
const validationResult = (errors) => deepFreeze({
  valid: errors.length === 0,
  errors
});
const addError = (errors, code, path, message) => {
  errors.push({ code, path, message });
};
const validateExactFields = (value, fields, path, errors) => {
  const expected = new Set(fields);
  Object.keys(value).sort(compareText).forEach((field) => {
    if (!expected.has(field)) {
      addError(errors, 'unknown-field', `${path}.${field}`, `Unknown field ${field}`);
    }
  });
  fields.forEach((field) => {
    if (!Object.hasOwn(value, field)) {
      addError(errors, 'missing-field', `${path}.${field}`, `Missing required field ${field}`);
    }
  });
};
const validateCanonicalString = (
  value,
  path,
  errors,
  { pattern = null, maximumLength = 1000, code = 'invalid-string' } = {}
) => {
  if (
    typeof value !== 'string'
    || !value
    || value !== value.trim()
    || value.length > maximumLength
    || /[\u0000-\u001f\u007f]/.test(value)
    || (pattern && !pattern.test(value))
  ) {
    addError(errors, code, path, 'Expected a canonical nonempty string');
    return false;
  }
  return true;
};
const validateSortedUniqueStrings = (
  values,
  path,
  errors,
  options = {}
) => {
  if (!Array.isArray(values)) {
    addError(errors, 'invalid-array', path, 'Expected an array');
    return [];
  }
  const seen = new Set();
  let previous = null;
  values.forEach((value, index) => {
    const valuePath = `${path}[${index}]`;
    if (!validateCanonicalString(value, valuePath, errors, options)) return;
    if (seen.has(value)) {
      addError(errors, 'duplicate-value', valuePath, 'Values must be unique');
    }
    if (previous !== null && compareText(previous, value) > 0) {
      addError(errors, 'unsorted-array', path, 'Values must use ascending canonical order');
    }
    seen.add(value);
    previous = value;
  });
  return values.filter((value) => typeof value === 'string');
};
const validateNonemptyTextArray = (values, path, errors) => {
  if (!Array.isArray(values) || !values.length) {
    addError(errors, 'invalid-array', path, 'Expected a nonempty array');
    return;
  }
  values.forEach((value, index) => {
    validateCanonicalString(value, `${path}[${index}]`, errors, { maximumLength: 1000 });
  });
};
const validateRecordOrder = (key, previous, path, errors) => {
  if (previous !== null && compareText(previous, key) > 0) {
    addError(errors, 'unsorted-array', path, 'Records must use ascending canonical key order');
  }
};
const checkpointReference = (journeyId, checkpointId) => `${journeyId}:${checkpointId}`;
const contractRepositoryPath = (reference) => (
  typeof reference === 'string' ? reference.slice(0, reference.indexOf('::')) : ''
);
const isAutomatedPrGateContract = (contract) => contract.execution === 'PR_GATE'
  && contract.sensitivity !== 'MANUAL_ACCEPTANCE';

const architectureRuleIndex = (architectureRegistry, detectorCatalog, errors) => {
  if (!isRecord(architectureRegistry) || !Array.isArray(architectureRegistry.rules)) {
    addError(errors, 'invalid-architecture-registry', '$architecture', 'Expected architecture rules');
    return new Map();
  }
  const index = new Map();
  architectureRegistry.rules.forEach((rule, ruleIndex) => {
    if (!isRecord(rule) || typeof rule.key !== 'string') return;
    const detector = detectorCatalog?.[rule.detector];
    const allowedInputs = rule.kind === 'single-semantic-owner'
      ? (Array.isArray(rule.allowedDefinitionPaths)
        ? rule.allowedDefinitionPaths.map((path) => ({ path })) : [])
      : rule.kind === 'single-canonical-entry-edge'
        ? (Array.isArray(rule.allowedEdges) ? rule.allowedEdges : [])
        : [];
    const allowedSubjects = [];
    if (!detector || typeof detector.encodeSubject !== 'function') {
      addError(
        errors,
        'invalid-architecture-detector',
        `$architecture.rules[${ruleIndex}].detector`,
        `Architecture rule ${rule.key} has no trusted stable subject encoder`
      );
    } else {
      allowedInputs.forEach((candidate) => {
        try {
          const subject = detector.encodeSubject(Object.freeze({ ...candidate }));
          if (typeof subject === 'string' && subject) allowedSubjects.push(subject);
        } catch (_error) {
          addError(
            errors,
            'invalid-architecture-subject',
            `$architecture.rules[${ruleIndex}]`,
            `Architecture rule ${rule.key} cannot encode an allowed subject`
          );
        }
      });
    }
    index.set(rule.key, {
      rule,
      detector,
      subjectCategory: detector?.subjectCategory || '',
      allowedSubjects: stableUnique(allowedSubjects)
    });
  });
  return index;
};

export const validateProductImpactMap = (
  map,
  architectureRegistry,
  detectorCatalog
) => {
  const errors = [];
  if (!isRecord(map)) {
    addError(errors, 'invalid-map', '$', 'Expected a Product Impact map object');
    return validationResult(errors);
  }
  validateExactFields(map, MAP_FIELDS, '$', errors);
  if (map.schemaVersion !== SCHEMA_VERSION) {
    addError(
      errors,
      'unsupported-schema-version',
      '$.schemaVersion',
      `Expected schemaVersion ${SCHEMA_VERSION}`
    );
  }
  const architectureRules = architectureRuleIndex(
    architectureRegistry,
    detectorCatalog,
    errors
  );
  const coverageByRule = new Map();
  if (!Array.isArray(map.ruleCoverage)) {
    addError(errors, 'invalid-array', '$.ruleCoverage', 'Expected an array');
  } else {
    let previousRuleKey = null;
    map.ruleCoverage.forEach((coverage, index) => {
      const path = `$.ruleCoverage[${index}]`;
      if (!isRecord(coverage)) {
        addError(errors, 'invalid-rule-coverage', path, 'Expected a rule coverage object');
        return;
      }
      validateExactFields(coverage, RULE_COVERAGE_FIELDS, path, errors);
      const validRuleKey = validateCanonicalString(
        coverage.ruleKey,
        `${path}.ruleKey`,
        errors,
        { pattern: RULE_KEY_PATTERN, code: 'invalid-rule-key' }
      );
      if (validRuleKey) {
        validateRecordOrder(coverage.ruleKey, previousRuleKey, '$.ruleCoverage', errors);
        if (coverageByRule.has(coverage.ruleKey)) {
          addError(errors, 'duplicate-rule-coverage', `${path}.ruleKey`, 'Duplicate rule coverage');
        }
        if (!architectureRules.has(coverage.ruleKey)) {
          addError(errors, 'unknown-architecture-rule', `${path}.ruleKey`, 'Unknown architecture rule');
        }
        coverageByRule.set(coverage.ruleKey, coverage);
        previousRuleKey = coverage.ruleKey;
      }
      if (!COVERAGE_LEVELS.includes(coverage.coverage)) {
        addError(errors, 'invalid-coverage', `${path}.coverage`, 'Unsupported coverage level');
      }
      if (!PRODUCT_ENFORCEMENTS.includes(coverage.enforcement)) {
        addError(errors, 'invalid-enforcement', `${path}.enforcement`, 'Unsupported enforcement');
      }
      if (coverage.enforcement === 'hard' && coverage.coverage !== 'complete') {
        addError(
          errors,
          'hard-requires-complete-coverage',
          path,
          'Hard Product Impact enforcement requires complete coverage'
        );
      }
    });
  }

  const mappedSubjectsByRule = new Map();
  const concernsByKey = new Map();
  if (!Array.isArray(map.concerns)) {
    addError(errors, 'invalid-array', '$.concerns', 'Expected an array');
  } else {
    let previousConcernKey = null;
    map.concerns.forEach((concern, concernIndex) => {
      const concernPath = `$.concerns[${concernIndex}]`;
      if (!isRecord(concern)) {
        addError(errors, 'invalid-concern', concernPath, 'Expected a concern object');
        return;
      }
      validateExactFields(concern, CONCERN_FIELDS, concernPath, errors);
      const validConcernKey = validateCanonicalString(
        concern.key,
        `${concernPath}.key`,
        errors,
        { pattern: RULE_KEY_PATTERN, code: 'invalid-concern-key' }
      );
      if (validConcernKey) {
        validateRecordOrder(concern.key, previousConcernKey, '$.concerns', errors);
        if (concernsByKey.has(concern.key)) {
          addError(errors, 'duplicate-concern-key', `${concernPath}.key`, 'Duplicate concern key');
        }
        concernsByKey.set(concern.key, concern);
        previousConcernKey = concern.key;
      }
      if (!Number.isInteger(concern.scenarioRevision) || concern.scenarioRevision < 1) {
        addError(
          errors,
          'invalid-scenario-revision',
          `${concernPath}.scenarioRevision`,
          'Expected a positive integer scenario revision'
        );
      }
      validateCanonicalString(concern.title, `${concernPath}.title`, errors);
      if (!CONCERN_LAYERS.includes(concern.layer)) {
        addError(errors, 'invalid-layer', `${concernPath}.layer`, 'Unsupported concern layer');
      }

      const checkpointRefs = new Set();
      const journeyIds = new Set();
      if (!Array.isArray(concern.journeys) || !concern.journeys.length) {
        addError(errors, 'invalid-journeys', `${concernPath}.journeys`, 'Expected journeys');
      } else {
        let previousJourneyId = null;
        concern.journeys.forEach((journey, journeyIndex) => {
          const journeyPath = `${concernPath}.journeys[${journeyIndex}]`;
          if (!isRecord(journey)) {
            addError(errors, 'invalid-journey', journeyPath, 'Expected a journey object');
            return;
          }
          validateExactFields(journey, JOURNEY_FIELDS, journeyPath, errors);
          const validJourneyId = validateCanonicalString(
            journey.id,
            `${journeyPath}.id`,
            errors,
            { pattern: /^JRN-[A-Z0-9]+(?:-[A-Z0-9]+)*$/, code: 'invalid-journey-id' }
          );
          if (validJourneyId) {
            validateRecordOrder(journey.id, previousJourneyId, `${concernPath}.journeys`, errors);
            if (journeyIds.has(journey.id)) {
              addError(errors, 'duplicate-journey-id', `${journeyPath}.id`, 'Duplicate journey ID');
            }
            journeyIds.add(journey.id);
            previousJourneyId = journey.id;
          }
          for (const field of ['actor', 'context', 'goal']) {
            validateCanonicalString(journey[field], `${journeyPath}.${field}`, errors);
          }
          validateNonemptyTextArray(journey.steps, `${journeyPath}.steps`, errors);
          if (!Array.isArray(journey.checkpoints) || !journey.checkpoints.length) {
            addError(errors, 'invalid-checkpoints', `${journeyPath}.checkpoints`, 'Expected checkpoints');
          } else {
            const localCheckpointIds = new Set();
            journey.checkpoints.forEach((checkpoint, checkpointIndex) => {
              const checkpointPath = `${journeyPath}.checkpoints[${checkpointIndex}]`;
              if (!isRecord(checkpoint)) {
                addError(errors, 'invalid-checkpoint', checkpointPath, 'Expected a checkpoint object');
                return;
              }
              validateExactFields(checkpoint, CHECKPOINT_FIELDS, checkpointPath, errors);
              const validCheckpoint = validateCanonicalString(
                checkpoint.id,
                `${checkpointPath}.id`,
                errors,
                { pattern: LOWER_ID_PATTERN, code: 'invalid-checkpoint-id' }
              );
              validateCanonicalString(checkpoint.label, `${checkpointPath}.label`, errors);
              if (validCheckpoint && validJourneyId) {
                if (localCheckpointIds.has(checkpoint.id)) {
                  addError(
                    errors,
                    'duplicate-checkpoint-id',
                    `${checkpointPath}.id`,
                    'Duplicate checkpoint ID in journey'
                  );
                }
                localCheckpointIds.add(checkpoint.id);
                checkpointRefs.add(checkpointReference(journey.id, checkpoint.id));
              }
            });
          }
        });
      }

      const effectsById = new Map();
      if (!Array.isArray(concern.effects) || !concern.effects.length) {
        addError(errors, 'invalid-effects', `${concernPath}.effects`, 'Expected stable effects');
      } else {
        let previousEffectId = null;
        concern.effects.forEach((effect, effectIndex) => {
          const effectPath = `${concernPath}.effects[${effectIndex}]`;
          if (!isRecord(effect)) {
            addError(errors, 'invalid-effect', effectPath, 'Expected an effect object');
            return;
          }
          validateExactFields(effect, EFFECT_FIELDS, effectPath, errors);
          const validEffectId = validateCanonicalString(
            effect.id,
            `${effectPath}.id`,
            errors,
            { pattern: /^EFF-[A-Z0-9]+(?:-[A-Z0-9]+)*$/, code: 'invalid-effect-id' }
          );
          if (validEffectId) {
            validateRecordOrder(effect.id, previousEffectId, `${concernPath}.effects`, errors);
            if (effectsById.has(effect.id)) {
              addError(errors, 'duplicate-effect-id', `${effectPath}.id`, 'Duplicate effect ID');
            }
            effectsById.set(effect.id, effect);
            previousEffectId = effect.id;
          }
          if (!checkpointRefs.has(effect.checkpointRef)) {
            addError(
              errors,
              'unknown-checkpoint-ref',
              `${effectPath}.checkpointRef`,
              'Effect references an unknown checkpoint'
            );
          }
          validateCanonicalString(effect.statement, `${effectPath}.statement`, errors);
        });
      }

      const optionIds = new Set();
      const requirementIds = new Set();
      const referencedRuleKeys = new Set();
      const hardRequirementCheckpointRefs = new Set();
      if (!Array.isArray(concern.options) || !concern.options.length) {
        addError(errors, 'invalid-options', `${concernPath}.options`, 'Expected behavior options');
      } else {
        let previousOptionId = null;
        concern.options.forEach((option, optionIndex) => {
          const optionPath = `${concernPath}.options[${optionIndex}]`;
          if (!isRecord(option)) {
            addError(errors, 'invalid-option', optionPath, 'Expected an option object');
            return;
          }
          validateExactFields(option, OPTION_FIELDS, optionPath, errors);
          const validOptionId = validateCanonicalString(
            option.id,
            `${optionPath}.id`,
            errors,
            { pattern: LOWER_ID_PATTERN, code: 'invalid-option-id' }
          );
          if (validOptionId) {
            validateRecordOrder(option.id, previousOptionId, `${concernPath}.options`, errors);
            if (optionIds.has(option.id)) {
              addError(errors, 'duplicate-option-id', `${optionPath}.id`, 'Duplicate option ID');
            }
            optionIds.add(option.id);
            previousOptionId = option.id;
          }
          validateCanonicalString(option.summary, `${optionPath}.summary`, errors);
          const effectRefs = validateSortedUniqueStrings(
            option.effectRefs,
            `${optionPath}.effectRefs`,
            errors,
            { pattern: /^EFF-[A-Z0-9]+(?:-[A-Z0-9]+)*$/, code: 'invalid-effect-ref' }
          );
          if (!effectRefs.length) {
            addError(errors, 'empty-effect-refs', `${optionPath}.effectRefs`, 'Expected effect refs');
          }
          effectRefs.forEach((effectRef, refIndex) => {
            if (!effectsById.has(effectRef)) {
              addError(
                errors,
                'unknown-effect-ref',
                `${optionPath}.effectRefs[${refIndex}]`,
                'Option references an unknown effect'
              );
            }
          });
          if (!Array.isArray(option.requirements) || !option.requirements.length) {
            addError(
              errors,
              'invalid-requirements',
              `${optionPath}.requirements`,
              'Expected at least one requirement'
            );
            return;
          }
          let previousRequirementId = null;
          option.requirements.forEach((requirement, requirementIndex) => {
            const requirementPath = `${optionPath}.requirements[${requirementIndex}]`;
            if (!isRecord(requirement)) {
              addError(errors, 'invalid-requirement', requirementPath, 'Expected a requirement');
              return;
            }
            validateExactFields(requirement, REQUIREMENT_FIELDS, requirementPath, errors);
            const validRequirementId = validateCanonicalString(
              requirement.id,
              `${requirementPath}.id`,
              errors,
              { pattern: /^REQ-[A-Z0-9]+(?:-[A-Z0-9]+)*$/, code: 'invalid-requirement-id' }
            );
            if (validRequirementId) {
              validateRecordOrder(
                requirement.id,
                previousRequirementId,
                `${optionPath}.requirements`,
                errors
              );
              if (requirementIds.has(requirement.id)) {
                addError(
                  errors,
                  'duplicate-requirement-id',
                  `${requirementPath}.id`,
                  'Requirement IDs must be unique within a concern'
                );
              }
              requirementIds.add(requirement.id);
              previousRequirementId = requirement.id;
            }
            if (!REQUIREMENT_CATEGORIES.includes(requirement.category)) {
              addError(
                errors,
                'invalid-requirement-category',
                `${requirementPath}.category`,
                'Unsupported requirement category'
              );
            }
            const requirementCheckpointRefs = validateSortedUniqueStrings(
              requirement.checkpointRefs,
              `${requirementPath}.checkpointRefs`,
              errors
            );
            if (!requirementCheckpointRefs.length) {
              addError(
                errors,
                'empty-checkpoint-refs',
                `${requirementPath}.checkpointRefs`,
                'Expected checkpoint refs'
              );
            }
            requirementCheckpointRefs.forEach((reference, referenceIndex) => {
              if (!checkpointRefs.has(reference)) {
                addError(
                  errors,
                  'unknown-checkpoint-ref',
                  `${requirementPath}.checkpointRefs[${referenceIndex}]`,
                  'Requirement references an unknown checkpoint'
                );
              }
            });
            if (!Array.isArray(requirement.anyOf) || !requirement.anyOf.length) {
              addError(
                errors,
                'invalid-candidates',
                `${requirementPath}.anyOf`,
                'Expected at least one subject candidate'
              );
              return;
            }
            const candidateKeys = new Set();
            let previousCandidateKey = null;
            let requirementUsesHardRule = false;
            requirement.anyOf.forEach((candidate, candidateIndex) => {
              const candidatePath = `${requirementPath}.anyOf[${candidateIndex}]`;
              if (!isRecord(candidate)) {
                addError(errors, 'invalid-candidate', candidatePath, 'Expected a subject candidate');
                return;
              }
              validateExactFields(candidate, CANDIDATE_FIELDS, candidatePath, errors);
              const validRuleKey = validateCanonicalString(
                candidate.ruleKey,
                `${candidatePath}.ruleKey`,
                errors,
                { pattern: RULE_KEY_PATTERN, code: 'invalid-rule-key' }
              );
              validateCanonicalString(
                candidate.subjectCategory,
                `${candidatePath}.subjectCategory`,
                errors,
                { pattern: LOWER_ID_PATTERN, code: 'invalid-subject-category' }
              );
              validateCanonicalString(candidate.subject, `${candidatePath}.subject`, errors);
              if (!validRuleKey || typeof candidate.subject !== 'string') return;
              const candidateKey = [
                candidate.ruleKey,
                candidate.subjectCategory,
                candidate.subject
              ].join('\u0000');
              validateRecordOrder(
                candidateKey,
                previousCandidateKey,
                `${requirementPath}.anyOf`,
                errors
              );
              if (candidateKeys.has(candidateKey)) {
                addError(
                  errors,
                  'duplicate-candidate',
                  candidatePath,
                  'Requirement candidates must be unique'
                );
              }
              candidateKeys.add(candidateKey);
              previousCandidateKey = candidateKey;
              referencedRuleKeys.add(candidate.ruleKey);
              if (!mappedSubjectsByRule.has(candidate.ruleKey)) {
                mappedSubjectsByRule.set(candidate.ruleKey, new Set());
              }
              mappedSubjectsByRule.get(candidate.ruleKey).add(candidate.subject);
              const architecture = architectureRules.get(candidate.ruleKey);
              const coverage = coverageByRule.get(candidate.ruleKey);
              if (!architecture) {
                addError(
                  errors,
                  'unknown-architecture-rule',
                  `${candidatePath}.ruleKey`,
                  'Candidate references an unknown architecture rule'
                );
              } else {
                if (candidate.subjectCategory !== architecture.subjectCategory) {
                  addError(
                    errors,
                    'subject-category-mismatch',
                    `${candidatePath}.subjectCategory`,
                    'Candidate subject category does not match its detector'
                  );
                }
                if (!architecture.allowedSubjects.includes(candidate.subject)) {
                  addError(
                    errors,
                    'unapproved-architecture-subject',
                    `${candidatePath}.subject`,
                    'Candidate subject is not allowed by the architecture rule'
                  );
                }
              }
              if (!coverage) {
                addError(
                  errors,
                  'missing-rule-coverage',
                  `${candidatePath}.ruleKey`,
                  'Referenced architecture rule has no Product Impact coverage entry'
                );
              } else if (coverage.enforcement === 'hard') {
                requirementUsesHardRule = true;
              }
            });
            if (requirementUsesHardRule) {
              requirementCheckpointRefs.forEach((reference) => {
                hardRequirementCheckpointRefs.add(reference);
              });
            }
          });
        });
      }

      if (!isRecord(concern.resolution)) {
        addError(errors, 'invalid-resolution', `${concernPath}.resolution`, 'Expected a resolution');
      } else {
        const kind = concern.resolution.kind;
        const fields = kind === 'decision-required'
          ? ['kind']
          : ['authorityRefs', 'kind', 'selectedOptionId'];
        validateExactFields(concern.resolution, fields, `${concernPath}.resolution`, errors);
        if (!RESOLUTION_KINDS.includes(kind)) {
          addError(
            errors,
            'invalid-resolution-kind',
            `${concernPath}.resolution.kind`,
            'Unsupported concern resolution kind'
          );
        }
        if (kind !== 'decision-required') {
          if (!optionIds.has(concern.resolution.selectedOptionId)) {
            addError(
              errors,
              'unknown-selected-option',
              `${concernPath}.resolution.selectedOptionId`,
              'Resolution references an unknown option'
            );
          }
          const authorityRefs = validateSortedUniqueStrings(
            concern.resolution.authorityRefs,
            `${concernPath}.resolution.authorityRefs`,
            errors
          );
          if (!authorityRefs.length) {
            addError(
              errors,
              'empty-authority-refs',
              `${concernPath}.resolution.authorityRefs`,
              'Static and domain resolutions require authority references'
            );
          }
        }
      }

      const contracts = [];
      if (!Array.isArray(concern.contracts)) {
        addError(errors, 'invalid-array', `${concernPath}.contracts`, 'Expected contracts');
      } else {
        const contractRefs = new Set();
        let previousContractRef = null;
        concern.contracts.forEach((contract, contractIndex) => {
          const contractPath = `${concernPath}.contracts[${contractIndex}]`;
          if (!isRecord(contract)) {
            addError(errors, 'invalid-contract', contractPath, 'Expected a contract object');
            return;
          }
          validateExactFields(contract, CONTRACT_FIELDS, contractPath, errors);
          const validContractRef = validateCanonicalString(
            contract.ref,
            `${contractPath}.ref`,
            errors,
            { maximumLength: 1000, code: 'invalid-contract-ref' }
          );
          if (validContractRef) {
            const separator = contract.ref.indexOf('::');
            const repositoryPath = contractRepositoryPath(contract.ref);
            if (
              separator < 1
              || separator === contract.ref.length - 2
              || !REPOSITORY_PATH_PATTERN.test(repositoryPath)
              || repositoryPath.split('/').some((segment) => segment === '.' || segment === '..')
            ) {
              addError(
                errors,
                'invalid-contract-ref',
                `${contractPath}.ref`,
                'Expected a normalized repository path and named check separated by ::'
              );
            }
            validateRecordOrder(
              contract.ref,
              previousContractRef,
              `${concernPath}.contracts`,
              errors
            );
            if (contractRefs.has(contract.ref)) {
              addError(errors, 'duplicate-contract-ref', `${contractPath}.ref`, 'Duplicate contract ref');
            }
            contractRefs.add(contract.ref);
            previousContractRef = contract.ref;
          }
          const contractCheckpointRefs = validateSortedUniqueStrings(
            contract.checkpointRefs,
            `${contractPath}.checkpointRefs`,
            errors
          );
          if (!contractCheckpointRefs.length) {
            addError(
              errors,
              'empty-checkpoint-refs',
              `${contractPath}.checkpointRefs`,
              'Expected checkpoint refs'
            );
          }
          contractCheckpointRefs.forEach((reference, referenceIndex) => {
            if (!checkpointRefs.has(reference)) {
              addError(
                errors,
                'unknown-checkpoint-ref',
                `${contractPath}.checkpointRefs[${referenceIndex}]`,
                'Contract references an unknown checkpoint'
              );
            }
          });
          if (!CONTRACT_SENSITIVITIES.includes(contract.sensitivity)) {
            addError(
              errors,
              'invalid-contract-sensitivity',
              `${contractPath}.sensitivity`,
              'Unsupported contract sensitivity'
            );
          }
          if (!CONTRACT_EXECUTIONS.includes(contract.execution)) {
            addError(
              errors,
              'invalid-contract-execution',
              `${contractPath}.execution`,
              'Unsupported contract execution boundary'
            );
          }
          validateCanonicalString(
            contract.residualRisk,
            `${contractPath}.residualRisk`,
            errors,
            { maximumLength: 2000 }
          );
          contracts.push(contract);
        });
      }

      hardRequirementCheckpointRefs.forEach((reference) => {
        const covered = contracts.some((contract) => (
          Array.isArray(contract.checkpointRefs)
          && contract.checkpointRefs.includes(reference)
          && isAutomatedPrGateContract(contract)
        ));
        if (!covered) {
          addError(
            errors,
            'missing-hard-contract-coverage',
            `${concernPath}.contracts`,
            `Hard checkpoint ${reference} requires automated PR_GATE coverage`
          );
        }
      });
      if (hardRequirementCheckpointRefs.size && contracts.length && contracts.every((contract) => (
        contract.sensitivity === 'MANUAL_ACCEPTANCE' || contract.execution === 'MANUAL'
      ))) {
        addError(
          errors,
          'manual-only-hard-concern',
          `${concernPath}.contracts`,
          'Manual-only evidence cannot establish hard completeness'
        );
      }
      referencedRuleKeys.forEach((ruleKey) => {
        if (!coverageByRule.has(ruleKey)) {
          addError(
            errors,
            'missing-rule-coverage',
            `${concernPath}.options`,
            `Concern references rule ${ruleKey} without coverage`
          );
        }
      });
    });
  }

  coverageByRule.forEach((coverage, ruleKey) => {
    if (coverage.coverage !== 'complete') return;
    const architecture = architectureRules.get(ruleKey);
    if (!architecture) return;
    const mapped = mappedSubjectsByRule.get(ruleKey) || new Set();
    architecture.allowedSubjects.forEach((subject) => {
      if (!mapped.has(subject)) {
        addError(
          errors,
          'unmapped-allowed-subject',
          '$.concerns',
          `Complete coverage for ${ruleKey} does not map an allowed subject`
        );
      }
    });
  });

  return validationResult(errors);
};

const concernIndex = (map) => new Map(
  (Array.isArray(map?.concerns) ? map.concerns : []).map((concern) => [concern.key, concern])
);
const checkpointSet = (concern) => new Set(concern.journeys.flatMap((journey) => (
  journey.checkpoints.map((checkpoint) => checkpointReference(journey.id, checkpoint.id))
)));
const requirementsForConcern = (concern) => concern.options.flatMap((option) => option.requirements);

export const validateProductDecisionAuthority = (authority, map) => {
  const errors = [];
  if (!isRecord(authority)) {
    addError(errors, 'invalid-decision-authority', '$', 'Expected a decision authority object');
    return validationResult(errors);
  }
  validateExactFields(authority, DECISION_AUTHORITY_FIELDS, '$', errors);
  if (authority.schemaVersion !== SCHEMA_VERSION) {
    addError(
      errors,
      'unsupported-schema-version',
      '$.schemaVersion',
      `Expected schemaVersion ${SCHEMA_VERSION}`
    );
  }
  const maintainerLogins = validateSortedUniqueStrings(
    authority.maintainerLogins,
    '$.maintainerLogins',
    errors,
    { pattern: LOGIN_PATTERN, maximumLength: 39, code: 'invalid-maintainer-login' }
  );
  if (!maintainerLogins.length) {
    addError(
      errors,
      'empty-maintainer-logins',
      '$.maintainerLogins',
      'Decision authority requires at least one canonical maintainer login'
    );
  }
  const concerns = concernIndex(map);
  const seenIds = new Set();
  const seenConcernRevisions = new Set();
  if (!Array.isArray(authority.decisions)) {
    addError(errors, 'invalid-array', '$.decisions', 'Expected a decisions array');
    return validationResult(errors);
  }
  let previousDecisionId = null;
  authority.decisions.forEach((decision, index) => {
    const path = `$.decisions[${index}]`;
    if (!isRecord(decision)) {
      addError(errors, 'invalid-decision', path, 'Expected a durable decision object');
      return;
    }
    validateExactFields(decision, DURABLE_DECISION_FIELDS, path, errors);
    const validId = validateCanonicalString(
      decision.id,
      `${path}.id`,
      errors,
      { pattern: DECISION_ID_PATTERN, code: 'invalid-decision-id' }
    );
    if (validId) {
      validateRecordOrder(decision.id, previousDecisionId, '$.decisions', errors);
      if (seenIds.has(decision.id)) {
        addError(errors, 'duplicate-decision-id', `${path}.id`, 'Duplicate decision ID');
      }
      seenIds.add(decision.id);
      previousDecisionId = decision.id;
    }
    validateCanonicalString(
      decision.concernKey,
      `${path}.concernKey`,
      errors,
      { pattern: RULE_KEY_PATTERN, code: 'invalid-concern-key' }
    );
    if (!Number.isInteger(decision.scenarioRevision) || decision.scenarioRevision < 1) {
      addError(
        errors,
        'invalid-scenario-revision',
        `${path}.scenarioRevision`,
        'Expected a positive integer scenario revision'
      );
    }
    const concernRevisionKey = `${decision.concernKey}\u0000${decision.scenarioRevision}`;
    if (seenConcernRevisions.has(concernRevisionKey)) {
      addError(
        errors,
        'duplicate-active-decision',
        path,
        'Only one active decision is allowed per concern and scenario revision'
      );
    }
    seenConcernRevisions.add(concernRevisionKey);
    if (!ACCEPTED_IMPACT_CLASSES.includes(decision.acceptedImpactClass)) {
      addError(
        errors,
        'invalid-accepted-impact-class',
        `${path}.acceptedImpactClass`,
        'Unsupported accepted impact class'
      );
    }
    const retiredRefs = validateSortedUniqueStrings(
      decision.acknowledgedRetiredRequirementRefs,
      `${path}.acknowledgedRetiredRequirementRefs`,
      errors,
      { pattern: /^REQ-[A-Z0-9]+(?:-[A-Z0-9]+)*$/, code: 'invalid-requirement-ref' }
    );
    if (retiredRefs.length && decision.acceptedImpactClass !== 'RETIREMENT') {
      addError(
        errors,
        'retirement-impact-mismatch',
        `${path}.acknowledgedRetiredRequirementRefs`,
        'Retired requirements require RETIREMENT impact'
      );
    }
    if (decision.acceptedImpactClass === 'RETIREMENT' && !retiredRefs.length) {
      addError(
        errors,
        'missing-retirement-refs',
        `${path}.acknowledgedRetiredRequirementRefs`,
        'RETIREMENT requires exact retired requirement references'
      );
    }
    const acceptanceRefs = validateSortedUniqueStrings(
      decision.acceptanceCheckpointRefs,
      `${path}.acceptanceCheckpointRefs`,
      errors
    );
    if (!acceptanceRefs.length) {
      addError(
        errors,
        'empty-acceptance-checkpoints',
        `${path}.acceptanceCheckpointRefs`,
        'Durable decisions require acceptance checkpoints'
      );
    }
    validateCanonicalString(
      decision.rationale,
      `${path}.rationale`,
      errors,
      { maximumLength: 2000 }
    );

    const concern = concerns.get(decision.concernKey);
    if (!concern) {
      addError(errors, 'unknown-concern', `${path}.concernKey`, 'Decision references an unknown concern');
      return;
    }
    if (concern.scenarioRevision !== decision.scenarioRevision) {
      addError(
        errors,
        'stale-scenario-revision',
        `${path}.scenarioRevision`,
        'Decision scenario revision does not match the active concern'
      );
    }
    if (concern.resolution.kind !== 'decision-required') {
      addError(
        errors,
        'decision-conflicts-with-resolution',
        `${path}.concernKey`,
        'Static and domain concerns cannot carry durable decisions'
      );
    }
    if (!concern.options.some((option) => option.id === decision.selectedOptionId)) {
      addError(
        errors,
        'unknown-selected-option',
        `${path}.selectedOptionId`,
        'Decision references an unknown option'
      );
    }
    const requirementIds = new Set(requirementsForConcern(concern).map(({ id }) => id));
    retiredRefs.forEach((reference, referenceIndex) => {
      if (!requirementIds.has(reference)) {
        addError(
          errors,
          'unknown-retired-requirement',
          `${path}.acknowledgedRetiredRequirementRefs[${referenceIndex}]`,
          'Decision references an unknown retired requirement'
        );
      }
    });
    const checkpoints = checkpointSet(concern);
    acceptanceRefs.forEach((reference, referenceIndex) => {
      if (!checkpoints.has(reference)) {
        addError(
          errors,
          'unknown-checkpoint-ref',
          `${path}.acceptanceCheckpointRefs[${referenceIndex}]`,
          'Decision references an unknown checkpoint'
        );
        return;
      }
      const covered = concern.contracts.some((contract) => contract.checkpointRefs.includes(reference));
      if (!covered) {
        addError(
          errors,
          'uncovered-acceptance-checkpoint',
          `${path}.acceptanceCheckpointRefs[${referenceIndex}]`,
          'Acceptance checkpoint has no mapped behavior contract'
        );
      }
    });
  });
  return validationResult(errors);
};

const canonicalSubjectArray = (values, label) => {
  if (!Array.isArray(values) || values.some((value) => (
    typeof value !== 'string'
    || !value
    || value !== value.trim()
    || /[\u0000-\u001f\u007f]/.test(value)
  ))) {
    throw new Error(`${label} must contain canonical subject strings`);
  }
  return stableUnique(values);
};
const normalizeSubjectDelta = (delta) => {
  if (!isRecord(delta)) throw new Error('Product Impact subject delta must be an object');
  const fields = Object.keys(delta).sort(compareText);
  if (!arraysEqual(fields, SUBJECT_DELTA_FIELDS)) {
    throw new Error(`Product Impact subject delta requires exactly: ${SUBJECT_DELTA_FIELDS.join(', ')}`);
  }
  for (const field of ['ruleKey', 'kind', 'detector', 'subjectCategory']) {
    if (typeof delta[field] !== 'string' || !delta[field] || delta[field] !== delta[field].trim()) {
      throw new Error(`Product Impact subject delta ${field} must be canonical`);
    }
  }
  const beforeSubjects = canonicalSubjectArray(delta.beforeSubjects, 'Subject delta beforeSubjects');
  const afterSubjects = canonicalSubjectArray(delta.afterSubjects, 'Subject delta afterSubjects');
  const addedSubjects = stableDifference(afterSubjects, beforeSubjects);
  const removedSubjects = stableDifference(beforeSubjects, afterSubjects);
  if (
    !arraysEqual(canonicalSubjectArray(delta.addedSubjects, 'Subject delta addedSubjects'), addedSubjects)
    || !arraysEqual(canonicalSubjectArray(delta.removedSubjects, 'Subject delta removedSubjects'), removedSubjects)
    || delta.changed !== Boolean(addedSubjects.length || removedSubjects.length)
  ) {
    throw new Error(`Product Impact subject delta ${delta.ruleKey} is internally inconsistent`);
  }
  return deepFreeze({
    ruleKey: delta.ruleKey,
    kind: delta.kind,
    detector: delta.detector,
    subjectCategory: delta.subjectCategory,
    beforeSubjects,
    addedSubjects,
    removedSubjects,
    afterSubjects,
    changed: delta.changed
  });
};
const normalizeRuleFact = (fact) => {
  if (!isRecord(fact)) throw new Error('Product Impact rule fact must be an object');
  const fields = Object.keys(fact).sort(compareText);
  if (!arraysEqual(fields, RULE_FACT_FIELDS)) {
    throw new Error(`Product Impact rule fact requires exactly: ${RULE_FACT_FIELDS.join(', ')}`);
  }
  if (
    typeof fact.ruleKey !== 'string'
    || !RULE_KEY_PATTERN.test(fact.ruleKey)
    || typeof fact.subjectCategory !== 'string'
    || !LOWER_ID_PATTERN.test(fact.subjectCategory)
  ) {
    throw new Error('Product Impact rule fact identifiers must be canonical');
  }
  return deepFreeze({
    ruleKey: fact.ruleKey,
    subjectCategory: fact.subjectCategory,
    beforeSubjects: canonicalSubjectArray(fact.beforeSubjects, 'Rule fact beforeSubjects'),
    afterSubjects: canonicalSubjectArray(fact.afterSubjects, 'Rule fact afterSubjects')
  });
};
const normalizeChangedPaths = (paths) => {
  if (!Array.isArray(paths) || paths.some((path) => (
    typeof path !== 'string'
    || !REPOSITORY_PATH_PATTERN.test(path)
    || path.split('/').some((segment) => segment === '.' || segment === '..')
  ))) {
    throw new Error('Product Impact changedPaths must contain normalized repository paths');
  }
  return stableUnique(paths);
};
const normalizeCurrentDecisionSource = (source) => {
  if (source === null || source === undefined) {
    return deepFreeze({
      valid: true,
      present: false,
      metadata: { eventAuthorLogin: '', eventHeadSha: '' },
      declaration: { schemaVersion: SCHEMA_VERSION, headSha: '', decisions: [] },
      errors: []
    });
  }
  if (!isRecord(source)) throw new Error('Current decisions must be a parser result');
  const fields = Object.keys(source).sort(compareText);
  if (!arraysEqual(fields, ['declaration', 'errors', 'metadata', 'present', 'valid'])) {
    throw new Error('Current decisions must use the bounded parser result shape');
  }
  if (source.valid !== true || !Array.isArray(source.errors) || source.errors.length) {
    throw new Error('Current decisions must be a valid bounded parser result');
  }
  if (!isRecord(source.metadata) || !isRecord(source.declaration)) {
    throw new Error('Current decisions parser metadata or declaration is malformed');
  }
  if (!Array.isArray(source.declaration.decisions)) {
    throw new Error('Current decisions parser declaration is malformed');
  }
  return source;
};
const resultSeverity = Object.freeze({
  AUTHORITY_CONFLICT: 5,
  INSUFFICIENT_EVIDENCE: 4,
  ORDINARY_REGRESSION: 3,
  UNRESOLVED_DECISION: 2,
  CONFORMING: 1,
  NOT_APPLICABLE: 0
});
const choiceObservation = (candidates) => [...candidates]
  .sort((left, right) => resultSeverity[right] - resultSeverity[left])[0];

export const evaluateProductImpact = (input) => {
  if (!isRecord(input)) throw new Error('Product Impact evaluation requires one input object');
  const fields = Object.keys(input).sort(compareText);
  if (!arraysEqual(fields, EVALUATION_INPUT_FIELDS)) {
    throw new Error(`Product Impact evaluation requires exactly: ${EVALUATION_INPUT_FIELDS.join(', ')}`);
  }
  if (!isRecord(input.map)) throw new Error('Product Impact evaluation requires a validated map');
  const decisionValidation = validateProductDecisionAuthority(input.decisionAuthority, input.map);
  if (!decisionValidation.valid) {
    throw new Error('Product Impact evaluation requires valid durable decision authority');
  }
  const subjectDeltas = input.subjectDeltas.map(normalizeSubjectDelta)
    .sort((left, right) => compareText(left.ruleKey, right.ruleKey));
  const duplicateDelta = subjectDeltas.find((delta, index) => (
    index > 0 && delta.ruleKey === subjectDeltas[index - 1].ruleKey
  ));
  if (duplicateDelta) throw new Error(`Duplicate Product Impact subject delta ${duplicateDelta.ruleKey}`);
  const ruleFacts = input.ruleFacts.map(normalizeRuleFact)
    .sort((left, right) => compareText(left.ruleKey, right.ruleKey));
  const factsByRule = new Map();
  ruleFacts.forEach((fact) => {
    if (factsByRule.has(fact.ruleKey)) throw new Error(`Duplicate Product Impact rule fact ${fact.ruleKey}`);
    factsByRule.set(fact.ruleKey, fact);
  });
  const changedPathSet = new Set(normalizeChangedPaths(input.changedPaths));
  const currentSource = normalizeCurrentDecisionSource(input.currentDecisions);
  const currentDecisionList = currentSource.declaration.decisions;
  const changedDeltas = subjectDeltas.filter(({ changed }) => changed);
  const changedRuleKeys = new Set(changedDeltas.map(({ ruleKey }) => ruleKey));
  const currentConcernKeys = new Set(currentDecisionList.map(({ concernKey }) => concernKey));
  const coverageByRule = new Map(input.map.ruleCoverage.map((coverage) => [
    coverage.ruleKey,
    coverage
  ]));
  const mappedSubjectsByRule = new Map();
  input.map.concerns.forEach((concern) => {
    concern.options.forEach((option) => option.requirements.forEach((requirement) => {
      requirement.anyOf.forEach((candidate) => {
        if (!mappedSubjectsByRule.has(candidate.ruleKey)) {
          mappedSubjectsByRule.set(candidate.ruleKey, new Set());
        }
        mappedSubjectsByRule.get(candidate.ruleKey).add(candidate.subject);
      });
    }));
  });
  const authorityDecisions = input.decisionAuthority.decisions;
  const maintainerLogins = new Set(input.decisionAuthority.maintainerLogins);
  const results = [];

  [...input.map.concerns]
    .sort((left, right) => compareText(left.key, right.key))
    .forEach((concern) => {
      const referencedRules = stableUnique(concern.options.flatMap((option) => (
        option.requirements.flatMap((requirement) => requirement.anyOf.map(({ ruleKey }) => ruleKey))
      )));
      const triggeringRuleKeys = referencedRules.filter((ruleKey) => changedRuleKeys.has(ruleKey));
      if (!triggeringRuleKeys.length && !currentConcernKeys.has(concern.key)) return;
      const triggeringDeltas = subjectDeltas.filter((delta) => (
        triggeringRuleKeys.includes(delta.ruleKey)
      ));
      const evidenceGaps = [];
      const activeFactsByRule = new Map();
      referencedRules.forEach((ruleKey) => {
        const fact = factsByRule.get(ruleKey);
        if (!fact) {
          evidenceGaps.push(`Missing complete before/after subject inventory for ${ruleKey}.`);
          activeFactsByRule.set(ruleKey, {
            beforeSubjects: [],
            afterSubjects: [],
            subjectCategory: ''
          });
          return;
        }
        activeFactsByRule.set(ruleKey, fact);
      });
      triggeringDeltas.forEach((delta) => {
        const fact = factsByRule.get(delta.ruleKey);
        if (fact && (
          fact.subjectCategory !== delta.subjectCategory
          || !arraysEqual(fact.beforeSubjects, delta.beforeSubjects)
          || !arraysEqual(fact.afterSubjects, delta.afterSubjects)
        )) {
          evidenceGaps.push(`Subject delta and complete rule facts disagree for ${delta.ruleKey}.`);
        }
        const coverage = coverageByRule.get(delta.ruleKey);
        if (!coverage) {
          evidenceGaps.push(`Changed rule ${delta.ruleKey} has no Product Impact coverage.`);
        } else if (coverage.coverage !== 'complete') {
          evidenceGaps.push(`Changed rule ${delta.ruleKey} has only partial Product Impact coverage.`);
        }
        const mapped = mappedSubjectsByRule.get(delta.ruleKey) || new Set();
        [...delta.addedSubjects, ...delta.removedSubjects].forEach((subject) => {
          if (!mapped.has(subject)) {
            evidenceGaps.push(`Changed subject for ${delta.ruleKey} is not mapped.`);
          }
        });
      });

      const requirementResults = [];
      const optionResults = concern.options.map((option) => {
        const optionRequirementResults = option.requirements.map((requirement) => {
          const beforeActiveSubjects = stableUnique(requirement.anyOf
            .filter((candidate) => (
              activeFactsByRule.get(candidate.ruleKey)?.beforeSubjects.includes(candidate.subject)
            ))
            .map((candidate) => candidate.subject));
          const afterActiveSubjects = stableUnique(requirement.anyOf
            .filter((candidate) => (
              activeFactsByRule.get(candidate.ruleKey)?.afterSubjects.includes(candidate.subject)
            ))
            .map((candidate) => candidate.subject));
          const result = {
            requirementRef: requirement.id,
            category: requirement.category,
            beforeSatisfied: beforeActiveSubjects.length > 0,
            afterSatisfied: afterActiveSubjects.length > 0,
            beforeActiveSubjects,
            afterActiveSubjects,
            changed: !arraysEqual(beforeActiveSubjects, afterActiveSubjects)
          };
          requirementResults.push(result);
          return result;
        });
        return {
          optionId: option.id,
          effectRefs: [...option.effectRefs],
          requirementRefs: option.requirements.map(({ id }) => id),
          beforeRealized: optionRequirementResults.every(({ beforeSatisfied }) => beforeSatisfied),
          afterRealized: optionRequirementResults.every(({ afterSatisfied }) => afterSatisfied)
        };
      });
      requirementResults.sort((left, right) => compareText(left.requirementRef, right.requirementRef));
      optionResults.sort((left, right) => compareText(left.optionId, right.optionId));
      const beforeOptions = optionResults.filter(({ beforeRealized }) => beforeRealized)
        .map(({ optionId }) => optionId);
      const afterOptions = optionResults.filter(({ afterRealized }) => afterRealized)
        .map(({ optionId }) => optionId);
      const beforeEffects = stableUnique(optionResults.filter(({ beforeRealized }) => beforeRealized)
        .flatMap(({ effectRefs }) => effectRefs));
      const afterEffects = stableUnique(optionResults.filter(({ afterRealized }) => afterRealized)
        .flatMap(({ effectRefs }) => effectRefs));
      const changedRequirementRefs = requirementResults.filter(({ changed }) => changed)
        .map(({ requirementRef }) => requirementRef);
      const lostProtectedRequirements = requirementResults.filter((requirement) => (
        ['AFFORDANCE', 'COMPATIBILITY'].includes(requirement.category)
        && requirement.beforeSatisfied
        && !requirement.afterSatisfied
      )).map(({ requirementRef }) => requirementRef);
      const lostEffects = stableDifference(beforeEffects, afterEffects);
      const addedEffects = stableDifference(afterEffects, beforeEffects);
      const affectedCheckpointRefs = stableUnique(concern.options.flatMap((option) => (
        option.requirements.filter((requirement) => changedRequirementRefs.includes(requirement.id))
          .flatMap(({ checkpointRefs }) => checkpointRefs)
      )));
      if (!affectedCheckpointRefs.length) {
        concern.effects.forEach(({ checkpointRef }) => affectedCheckpointRefs.push(checkpointRef));
        affectedCheckpointRefs.sort(compareText);
      }
      const hardAffected = triggeringRuleKeys.some((ruleKey) => (
        coverageByRule.get(ruleKey)?.enforcement === 'hard'
      ));
      const contractResults = concern.contracts.map((contract) => ({
        ref: contract.ref,
        checkpointRefs: [...contract.checkpointRefs],
        sensitivity: contract.sensitivity,
        execution: contract.execution,
        candidateModified: changedPathSet.has(contractRepositoryPath(contract.ref)),
        automatedPrGate: isAutomatedPrGateContract(contract)
      })).sort((left, right) => compareText(left.ref, right.ref));
      affectedCheckpointRefs.forEach((checkpointRef) => {
        const covering = concern.contracts.filter((contract) => (
          contract.checkpointRefs.includes(checkpointRef)
        ));
        if (!covering.length) {
          evidenceGaps.push(`Affected checkpoint ${checkpointRef} has no mapped contract.`);
        }
        if (hardAffected) {
          const unmodifiedAutomated = covering.some((contract) => (
            isAutomatedPrGateContract(contract)
            && !changedPathSet.has(contractRepositoryPath(contract.ref))
          ));
          if (!unmodifiedAutomated) {
            evidenceGaps.push(
              `Hard checkpoint ${checkpointRef} lacks unmodified automated PR_GATE evidence.`
            );
          }
        }
      });

      const stableEvidenceGaps = stableUnique(evidenceGaps);
      const strictNoDifference = beforeOptions.length > 0
        && arraysEqual(beforeOptions, afterOptions)
        && requirementResults.every((requirement) => (
          requirement.beforeSatisfied === requirement.afterSatisfied
        ))
        && !lostProtectedRequirements.length
        && arraysEqual(beforeEffects, afterEffects)
        && !stableEvidenceGaps.length;
      let impactClass;
      if (strictNoDifference) impactClass = 'NO_USER_VISIBLE_DIFFERENCE';
      else if (lostEffects.length || lostProtectedRequirements.length) impactClass = 'RETIREMENT';
      else if (
        arraysEqual(beforeEffects, afterEffects)
        && beforeEffects.length
      ) impactClass = 'AFFORDANCE_PRESERVED';
      else impactClass = 'PRODUCT_CHANGE';

      const durable = authorityDecisions.find((decision) => (
        decision.concernKey === concern.key
        && decision.scenarioRevision === concern.scenarioRevision
      ));
      const current = currentDecisionList.find((decision) => (
        decision.concernKey === concern.key
        && decision.scenarioRevision === concern.scenarioRevision
      ));
      const staticAuthority = concern.resolution.kind === 'authority-covered'
        ? {
          type: 'STATIC_AUTHORITY',
          selectedOptionId: concern.resolution.selectedOptionId,
          rationale: ''
        }
        : concern.resolution.kind === 'domain-derived'
          ? {
            type: 'DOMAIN_AUTHORITY',
            selectedOptionId: concern.resolution.selectedOptionId,
            rationale: ''
          }
          : null;
      const durableAuthority = durable ? {
        type: 'DURABLE_DECISION',
        selectedOptionId: durable.selectedOptionId,
        rationale: durable.rationale
      } : null;
      const currentIssues = [];
      let currentEligible = false;
      if (current) {
        if (!maintainerLogins.has(currentSource.metadata.eventAuthorLogin)) {
          currentIssues.push('The pull request author is not an eligible Product Decision Owner.');
        }
        if (concern.resolution.kind !== 'decision-required') {
          currentIssues.push('Current decisions apply only to decision-required concerns.');
        }
        if (current.acceptedImpactClass !== 'AFFORDANCE_PRESERVED') {
          currentIssues.push('Current decisions accept only AFFORDANCE_PRESERVED.');
        }
        if (impactClass !== 'AFFORDANCE_PRESERVED') {
          currentIssues.push('The evaluated transition is not affordance-preserving.');
        }
        const selectedOption = optionResults.find(({ optionId }) => (
          optionId === current.selectedOptionId
        ));
        if (!selectedOption?.afterRealized) {
          currentIssues.push('The current decision selected option is not realized after the change.');
        }
        if (selectedOption && !arraysEqual(beforeEffects, selectedOption.effectRefs)) {
          currentIssues.push('The current decision selected option does not preserve stable effects.');
        }
        if (lostProtectedRequirements.length) {
          currentIssues.push('A current decision cannot retire affordance or compatibility requirements.');
        }
        if (!arraysEqual(
          current.acknowledgedChangedRequirementRefs,
          changedRequirementRefs
        )) {
          currentIssues.push('Changed requirement acknowledgement is incomplete or stale.');
        }
        if (stableEvidenceGaps.length) {
          currentIssues.push('Mapped Product Impact evidence is incomplete.');
        }
        currentEligible = currentIssues.length === 0;
      }

      const conflictCandidates = [staticAuthority, durableAuthority].filter(Boolean);
      if (current && conflictCandidates.some((authority) => (
        authority.selectedOptionId !== current.selectedOptionId
      ))) {
        conflictCandidates.push({
          type: 'CURRENT_MAINTAINER_DECISION',
          selectedOptionId: current.selectedOptionId,
          rationale: current.rationale
        });
      }
      const conflict = new Set(conflictCandidates.map(({ selectedOptionId }) => selectedOptionId)).size > 1;
      let recognizedAuthority = staticAuthority || durableAuthority;
      if (!recognizedAuthority && currentEligible) {
        recognizedAuthority = {
          type: 'CURRENT_MAINTAINER_DECISION',
          selectedOptionId: current.selectedOptionId,
          rationale: current.rationale
        };
      }
      if (
        durableAuthority
        && impactClass !== 'NO_USER_VISIBLE_DIFFERENCE'
        && (
          durable.acceptedImpactClass !== impactClass
          || !arraysEqual(
            durable.acknowledgedRetiredRequirementRefs,
            lostProtectedRequirements
          )
        )
      ) {
        recognizedAuthority = staticAuthority;
        currentIssues.push(
          'The durable decision impact class or retired requirement scope does not match the evaluated transition.'
        );
      }
      const selectedOptionResult = recognizedAuthority
        ? optionResults.find(({ optionId }) => optionId === recognizedAuthority.selectedOptionId)
        : null;
      const observationCandidates = [];
      if (stableEvidenceGaps.length) observationCandidates.push('INSUFFICIENT_EVIDENCE');
      if (conflict) observationCandidates.push('AUTHORITY_CONFLICT');
      if (recognizedAuthority && !selectedOptionResult?.afterRealized) {
        observationCandidates.push('ORDINARY_REGRESSION');
      }
      const changeNeedsDecision = concern.resolution.kind === 'decision-required'
        && impactClass !== 'NO_USER_VISIBLE_DIFFERENCE';
      if (changeNeedsDecision && !recognizedAuthority) {
        observationCandidates.push('UNRESOLVED_DECISION');
      }
      const observation = observationCandidates.length
        ? choiceObservation(observationCandidates)
        : 'CONFORMING';
      const authorityResolution = conflict
        ? 'CONFLICT'
        : recognizedAuthority?.type || 'UNRESOLVED';
      const nextAction = observation === 'CONFORMING'
        ? 'No Product Impact action is required.'
        : observation === 'ORDINARY_REGRESSION'
          ? 'Restore the selected option realization before integration.'
          : observation === 'UNRESOLVED_DECISION'
            ? 'Obtain an eligible Product Decision Owner disposition through the documented route.'
            : observation === 'AUTHORITY_CONFLICT'
              ? 'Resolve the incompatible product authorities in an authority-only change.'
              : 'Complete the mapped rule facts and contract evidence before relying on this result.';
      const residualRisks = stableUnique([
        ...concern.contracts.map(({ residualRisk }) => residualRisk),
        ...(current?.residualRisk ? [current.residualRisk] : [])
      ]);
      results.push({
        concernKey: concern.key,
        scenarioRevision: concern.scenarioRevision,
        triggeringRuleKeys,
        subjectDeltas: triggeringDeltas,
        journeys: concern.journeys.map((journey) => ({
          id: journey.id,
          actor: journey.actor,
          context: journey.context,
          goal: journey.goal,
          steps: [...journey.steps],
          checkpoints: journey.checkpoints.map((checkpoint) => ({ ...checkpoint }))
        })),
        checkpoints: stableUnique([...checkpointSet(concern)]),
        requirementResults,
        optionResults,
        beforeOptions,
        afterOptions,
        beforeEffects,
        afterEffects,
        impactClass,
        authorityResolution,
        selectedOptionId: conflict ? null : recognizedAuthority?.selectedOptionId || null,
        decisionRationale: recognizedAuthority?.rationale || current?.rationale || '',
        observation,
        contractResults,
        residualRisks,
        decisionRequired: observation === 'UNRESOLVED_DECISION',
        nextAction,
        changedRequirementRefs,
        lostRequirementRefs: lostProtectedRequirements,
        addedEffectRefs: addedEffects,
        lostEffectRefs: lostEffects,
        evidenceGaps: stableEvidenceGaps,
        currentDecisionIssues: stableUnique(currentIssues)
      });
    });

  return deepFreeze(results);
};
