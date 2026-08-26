const PULL_REQUEST_EVENTS = new Set(['pull_request', 'pull_request_target']);

const isObject = (value) => value !== null
  && typeof value === 'object'
  && !Array.isArray(value);

const requireString = (value, field, errors) => {
  if (typeof value === 'string' && value.length) return value;
  errors.push(`GitHub pull request metadata is malformed: ${field} is missing.`);
  return '';
};

export const isPullRequestEvent = (eventName) => PULL_REQUEST_EVENTS.has(eventName);

export const classifyWebChangeContext = ({
  eventName,
  currentRepository,
  eventPayloadSource,
  baseSha,
  headSha
}) => {
  const ordinary = { context: 'ORDINARY', isPromotion: false, errors: [] };
  if (!isPullRequestEvent(eventName)) return ordinary;

  const errors = [];
  let payload;
  try {
    payload = JSON.parse(eventPayloadSource);
  } catch (_error) {
    return {
      ...ordinary,
      errors: ['GitHub pull request event payload is missing or malformed.']
    };
  }

  if (!isObject(payload) || !isObject(payload.pull_request)) {
    return {
      ...ordinary,
      errors: ['GitHub pull request event payload is missing or malformed.']
    };
  }

  const pullRequest = payload.pull_request;
  const base = isObject(pullRequest.base) ? pullRequest.base : {};
  const head = isObject(pullRequest.head) ? pullRequest.head : {};
  const baseRepository = isObject(base.repo) ? base.repo : {};
  const headRepository = isObject(head.repo) ? head.repo : {};
  const eventBaseRef = requireString(base.ref, 'pull_request.base.ref', errors);
  const eventHeadRef = requireString(head.ref, 'pull_request.head.ref', errors);
  const eventBaseSha = requireString(base.sha, 'pull_request.base.sha', errors);
  const eventHeadSha = requireString(head.sha, 'pull_request.head.sha', errors);
  const eventBaseRepository = requireString(
    baseRepository.full_name,
    'pull_request.base.repo.full_name',
    errors
  );
  const eventHeadRepository = requireString(
    headRepository.full_name,
    'pull_request.head.repo.full_name',
    errors
  );
  const repository = requireString(currentRepository, 'GITHUB_REPOSITORY', errors);

  if (eventBaseSha && baseSha !== eventBaseSha) {
    errors.push('Checker base SHA does not match the GitHub event payload.');
  }
  if (eventHeadSha && headSha !== eventHeadSha) {
    errors.push('Checker head SHA does not match the GitHub event payload.');
  }

  const isPromotion = !errors.length
    && eventBaseRef === 'main'
    && eventHeadRef === 'dev'
    && eventBaseRepository === repository
    && eventHeadRepository === repository;

  return {
    context: isPromotion ? 'PROMOTION' : 'ORDINARY',
    isPromotion,
    errors
  };
};
