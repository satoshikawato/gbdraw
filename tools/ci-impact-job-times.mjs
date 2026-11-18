// SELECTIVE_CI.md "Job times against their limits": a gate lists each finished job of its run
// attempt, and each step with its own limit, that used at least WARNING_SHARE of its
// `timeout-minutes`. It reports and never fails.

export const WARNING_SHARE = 0.8;
// GitHub's limit for a job that sets no timeout-minutes.
const DEFAULT_TIMEOUT_MINUTES = 360;
const EXPRESSION = /\$\{\{[\s\S]*?\}\}/;
const escapeRegExp = (value) => value.replace(/[.*+?^${}()|[\]\\]/g, '\\$&');
// A one-line YAML scalar: quotes removed ('' stands for ' inside single quotes), a trailing comment dropped.
const yamlScalar = (text) => {
  const quoted = text.match(/^"((?:[^"\\]|\\.)*)"|^'((?:[^']|'')*)'/);
  if (quoted) return quoted[1] ?? quoted[2].replace(/''/g, "'");
  return text.replace(/[ \t]+#.*$/, '');
};
// An expression in a `name:` template matches any text; `suffix` is a pattern source.
const namePattern = (template, suffix = '') => {
  const parts = template.split(EXPRESSION);
  return { pattern: new RegExp(`^${parts.map(escapeRegExp).join('.+')}${suffix}$`), literalLength: parts.join('').length };
};
const bestMatch = (limits, name) => limits.filter(({ pattern }) => pattern.test(name))
  .sort((left, right) => right.literalLength - left.literalLength)[0];

/**
 * @typedef {object} StepLimit
 * @property {string} name The step's `name:` template.
 * @property {RegExp} pattern Matches the step's name on the run page.
 * @property {number} literalLength Characters of the template outside expressions.
 * @property {number} limitMinutes Its own `timeout-minutes`.
 */

/** @returns {StepLimit[]} The steps of a job block that set their own `timeout-minutes`. */
const stepLimits = (block) => {
  const steps = block.split(/^ {4}steps:[ \t]*$/m)[1];
  if (steps === undefined) return [];
  return steps.split(/^ {4}(?! )/m)[0].split(/^(?= {6}- )/m)
    .map((step) => [
      step.match(/^ {6}- name:[ \t]*(.+?)[ \t]*$/m)?.[1] ?? step.match(/^ {8}name:[ \t]*(.+?)[ \t]*$/m)?.[1],
      step.match(/^ {8}timeout-minutes:[ \t]*(\d+)[ \t]*$/m)?.[1]
    ])
    .filter(([name, timeout]) => name && timeout)
    .map(([name, timeout]) => [yamlScalar(name), timeout])
    .map(([name, timeout]) => Object.freeze({ name, ...namePattern(name), limitMinutes: Number(timeout) }));
};

/**
 * @typedef {object} JobLimit
 * @property {string} id The job ID in the workflow.
 * @property {RegExp} pattern Matches the name of each job the ID starts; an expression in the
 *   `name:` template matches any text.
 * @property {number} literalLength Characters of the template outside expressions, so the most
 *   specific pattern wins when two match.
 * @property {number} limitMinutes Its `timeout-minutes`, or GitHub's default.
 * @property {StepLimit[]} steps Its steps with their own limits.
 */

/** @returns {JobLimit[]} Each job of a workflow's text, from its two-space-indented ID line. */
export const workflowJobLimits = (workflow) => {
  const body = String(workflow).split(/^jobs:[ \t]*$/m)[1];
  if (body === undefined) throw Object.assign(new Error('The workflow has no jobs section.'), { code: 'NO_JOBS_SECTION' });
  return body.split(/^(?= {2}[A-Za-z0-9_-]+:[ \t]*$)/m)
    .map((block) => [block.match(/^ {2}([A-Za-z0-9_-]+):[ \t]*$/m)?.[1], block])
    .filter(([id]) => id)
    .map(([id, block]) => {
      const name = block.match(/^ {4}name:[ \t]*(.+?)[ \t]*$/m)?.[1];
      const timeout = block.match(/^ {4}timeout-minutes:[ \t]*(\d+)[ \t]*$/m)?.[1];
      return Object.freeze({
        id,
        // Without a name, a matrix job is named `<id> (<values>)`.
        ...(name ? namePattern(yamlScalar(name)) : namePattern(id, '(?: \\(.+\\))?')),
        limitMinutes: timeout ? Number(timeout) : DEFAULT_TIMEOUT_MINUTES,
        steps: stepLimits(block)
      });
    });
};

/**
 * @typedef {object} JobTime
 * @property {'job'|'step'} kind A job, or a step with its own limit.
 * @property {string} name The job's name on the run page; for a step, `<job> / <step>`.
 * @property {string} url The job's page; for a step, its line on that page.
 * @property {number} minutes From its start to its end.
 * @property {number|null} limitMinutes Null when no workflow job matches the name.
 * @property {number|null} share minutes / limitMinutes.
 */

const ran = (item) => item.status === 'completed' && item.conclusion !== 'skipped'
  && Number.isFinite(Date.parse(item.started_at)) && Number.isFinite(Date.parse(item.completed_at));
const timeOf = (kind, name, url, item, limitMinutes) => {
  const minutes = (Date.parse(item.completed_at) - Date.parse(item.started_at)) / 60000;
  return Object.freeze({ kind, name, url, minutes, limitMinutes, share: limitMinutes === null ? null : minutes / limitMinutes });
};

/**
 * @returns {JobTime[]} The finished jobs that ran, each followed by its finished steps that have
 *   their own limit, from the run attempt's API job objects.
 */
export const jobTimes = ({ jobs, limits }) => jobs.filter(ran).flatMap((job) => {
  const limit = bestMatch(limits, job.name);
  const steps = limit && Array.isArray(job.steps) ? job.steps.filter(ran) : [];
  return [
    timeOf('job', job.name, job.html_url, job, limit ? limit.limitMinutes : null),
    ...steps.flatMap((step) => {
      const stepLimit = bestMatch(limit.steps, step.name);
      return stepLimit ? [timeOf('step', `${job.name} / ${step.name}`,
        `${job.html_url}#step:${step.number}:1`, step, stepLimit.limitMinutes)] : [];
    })
  ];
});

// Rounded down, so a share under the warning share never shows as reaching it.
const percent = (share) => `${Math.floor(share * 100 + 1e-9)}%`;
const describe = ({ name, minutes, limitMinutes, share }) => (
  `${name} took ${minutes.toFixed(1)} of its ${limitMinutes} minutes (${percent(share)})`
);
const markdownText = (text) => String(text).replace(/[\\`*_[\]<>|]/g, '\\$&');
const annotationText = (text) => String(text).replace(/%/g, '%25').replace(/\r/g, '%0D').replace(/\n/g, '%0A');

/**
 * @param {{ times?: JobTime[], failure?: { code: string, reason: string } }} report
 * @returns {{ summary: string, annotations: string[] }}
 */
export const formatJobTimeReport = ({ times = [], failure = null }) => {
  const lines = ['## Job times against their limits', ''];
  if (failure) {
    const text = `Job times were not read: [${failure.code}] ${failure.reason}`;
    lines.push(`- ${markdownText(text)}`, '');
    return { summary: lines.join('\n'), annotations: [`::warning title=Job times not read::${annotationText(text)}`] };
  }
  const limited = times.filter(({ share }) => share !== null).sort((left, right) => right.share - left.share);
  const near = limited.filter(({ share }) => share >= WARNING_SHARE);
  if (near.length) {
    lines.push(`Jobs and steps that used ${percent(WARNING_SHARE)} or more of their \`timeout-minutes\`:`, '');
    near.forEach((time) => lines.push(`- [${markdownText(time.name)}](${time.url}): ${time.minutes.toFixed(1)} of ${time.limitMinutes} min (${percent(time.share)})`));
  } else {
    lines.push(`- No job or step used ${percent(WARNING_SHARE)} or more of its \`timeout-minutes\`${limited.length
      ? `; the highest share: ${markdownText(describe(limited[0]))}.` : '.'}`);
  }
  const unmatched = times.filter(({ share }) => share === null);
  if (unmatched.length) lines.push(`- No workflow job matches: ${unmatched.map(({ name }) => markdownText(name)).join(', ')}`);
  lines.push('');
  return {
    summary: lines.join('\n'),
    annotations: near.map((time) => `::warning title=${time.kind === 'step' ? 'Step' : 'Job'} near its time limit::${annotationText(describe(time))}.`)
  };
};
