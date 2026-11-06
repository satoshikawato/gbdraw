// SELECTIVE_CI.md "Job times against their limits": a gate lists each finished job of its run
// attempt that used at least WARNING_SHARE of its `timeout-minutes`. It reports and never fails.

export const WARNING_SHARE = 0.8;
// GitHub's limit for a job that sets no timeout-minutes.
const DEFAULT_TIMEOUT_MINUTES = 360;
const EXPRESSION = /\$\{\{[\s\S]*?\}\}/;
const escapeRegExp = (value) => value.replace(/[.*+?^${}()|[\]\\]/g, '\\$&');

/**
 * @typedef {object} JobLimit
 * @property {string} id The job ID in the workflow.
 * @property {RegExp} pattern Matches the name of each job the ID starts; an expression in the
 *   `name:` template matches any text.
 * @property {number} literalLength Characters of the template outside expressions, so the most
 *   specific pattern wins when two match.
 * @property {number} limitMinutes Its `timeout-minutes`, or GitHub's default.
 */

/** @returns {JobLimit[]} Each job of a workflow's text, from its two-space-indented ID line. */
export const workflowJobLimits = (workflow) => {
  const body = String(workflow).split(/^jobs:[ \t]*$/m)[1];
  if (body === undefined) throw Object.assign(new Error('The workflow has no jobs section.'), { code: 'NO_JOBS_SECTION' });
  return body.split(/^(?= {2}[A-Za-z0-9_-]+:[ \t]*$)/m)
    .map((block) => [block.match(/^ {2}([A-Za-z0-9_-]+):[ \t]*$/m)?.[1], block])
    .filter(([id]) => id)
    .map(([id, block]) => {
      const name = block.match(/^ {4}name:[ \t]*(.+?)[ \t]*$/m)?.[1].replace(/^(['"])(.*)\1$/, '$2');
      const timeout = block.match(/^ {4}timeout-minutes:[ \t]*(\d+)[ \t]*$/m)?.[1];
      // Without a name, a matrix job is named `<id> (<values>)`.
      const parts = name ? name.split(EXPRESSION) : [id];
      const source = parts.map(escapeRegExp).join('.+') + (name ? '' : '(?: \\(.+\\))?');
      return Object.freeze({
        id,
        pattern: new RegExp(`^${source}$`),
        literalLength: parts.join('').length,
        limitMinutes: timeout ? Number(timeout) : DEFAULT_TIMEOUT_MINUTES
      });
    });
};

/**
 * @typedef {object} JobTime
 * @property {string} name The job's name on the run page.
 * @property {string} url The job's page.
 * @property {number} minutes From its start to its end.
 * @property {number|null} limitMinutes Null when no workflow job matches the name.
 * @property {number|null} share minutes / limitMinutes.
 */

/** @returns {JobTime[]} The finished jobs that ran, from the run attempt's API job objects. */
export const jobTimes = ({ jobs, limits }) => jobs
  .filter((job) => job.status === 'completed' && job.conclusion !== 'skipped'
    && Number.isFinite(Date.parse(job.started_at)) && Number.isFinite(Date.parse(job.completed_at)))
  .map((job) => {
    const limit = limits.filter(({ pattern }) => pattern.test(job.name))
      .sort((left, right) => right.literalLength - left.literalLength)[0];
    const minutes = (Date.parse(job.completed_at) - Date.parse(job.started_at)) / 60000;
    return Object.freeze({
      name: job.name,
      url: job.html_url,
      minutes,
      limitMinutes: limit ? limit.limitMinutes : null,
      share: limit ? minutes / limit.limitMinutes : null
    });
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
    lines.push(`Jobs that used ${percent(WARNING_SHARE)} or more of their \`timeout-minutes\`:`, '');
    near.forEach((time) => lines.push(`- [${markdownText(time.name)}](${time.url}): ${time.minutes.toFixed(1)} of ${time.limitMinutes} min (${percent(time.share)})`));
  } else {
    lines.push(`- No job used ${percent(WARNING_SHARE)} or more of its \`timeout-minutes\`${limited.length
      ? `; the highest share: ${markdownText(describe(limited[0]))}.` : '.'}`);
  }
  const unmatched = times.filter(({ share }) => share === null);
  if (unmatched.length) lines.push(`- No workflow job matches: ${unmatched.map(({ name }) => markdownText(name)).join(', ')}`);
  lines.push('');
  return {
    summary: lines.join('\n'),
    annotations: near.map((time) => `::warning title=Job near its time limit::${annotationText(describe(time))}.`)
  };
};
