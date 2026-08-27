import assert from 'node:assert/strict';
import { spawnSync } from 'node:child_process';
import { mkdtempSync, rmSync, writeFileSync } from 'node:fs';
import { tmpdir } from 'node:os';
import { join, resolve } from 'node:path';
import test, { after } from 'node:test';

import {
  evaluateHookCommand,
  extractPlainLanguageSummary,
  hookDenialOutput,
  validateBody,
  validateProposal,
  validateTitle
} from '../../tools/check-pr-language.mjs';

const CHECKER = resolve('tools/check-pr-language.mjs');
const ORIGINAL_TITLE = 'Promote steady-state CI topology and finalize main admission';
const CLEAR_TITLE = 'Stop rerunning full CI on dev-to-main pull requests';
const ORIGINAL_SUMMARY = `## Plain-language summary

This promotion removes the transition-only CI topology from main and completes the admission-policy cutover.`;
const CLEAR_BODY = `## Plain-language summary

This PR removes CI jobs that existed only while the new main-merge check was introduced. The exact dev commit still has to pass full staging, and main will require \`Promotion / gate\` and \`CodeQL\` before merge.

## Change class

- [x] GOVERNANCE`;

const temporaryDirectory = mkdtempSync(join(tmpdir(), 'gbdraw-pr-language-'));
const clearBodyFile = join(temporaryDirectory, 'clear-body.md');
writeFileSync(clearBodyFile, CLEAR_BODY);

after(() => rmSync(temporaryDirectory, { recursive: true, force: true }));

const hookEvent = (command) => JSON.stringify({
  hook_event_name: 'PreToolUse',
  tool_name: 'Bash',
  tool_input: { command }
});

test('PR 394 original title fails with each matched opaque phrase identified', () => {
  const errors = validateTitle(ORIGINAL_TITLE);
  assert.ok(errors.some((error) => error.includes('steady-state')));
  assert.ok(errors.some((error) => error.includes('CI topology')));
  assert.ok(errors.some((error) => error.includes('main admission')));
  errors.forEach((error) => assert.match(error, /Replace it with/));
});

test('PR 394 original opening sentence fails in the plain-language summary', () => {
  const errors = validateBody(ORIGINAL_SUMMARY);
  assert.ok(errors.some((error) => error.includes('transition-only')));
  assert.ok(errors.some((error) => error.includes('CI topology')));
  assert.ok(errors.some((error) => error.includes('admission-policy')));
  assert.ok(errors.some((error) => error.includes('cutover')));
});

test('a title with a concrete action and object passes', () => {
  assert.deepEqual(validateTitle(CLEAR_TITLE), []);
});

test('backticks do not hide opaque language in a PR title', () => {
  const errors = validateTitle('Remove the `transition-only` CI topology');
  assert.ok(errors.some((error) => error.includes('transition-only')));
  assert.ok(errors.some((error) => error.includes('CI topology')));
});

test('a clear summary with exact backticked check names passes', () => {
  assert.deepEqual(validateBody(CLEAR_BODY), []);
});

test('opaque terms in a later technical section are outside the checker scope', () => {
  const body = `${CLEAR_BODY}

## Technical details

This completes the admission-policy cutover and leaves the steady-state CI topology in place.`;
  assert.deepEqual(validateBody(body), []);
  assert.match(extractPlainLanguageSummary(body), /exact dev commit/);
  assert.doesNotMatch(extractPlainLanguageSummary(body), /cutover/);
});

test('fenced and inline code in the summary are excluded from opaque-language matching', () => {
  const body = `## Plain-language summary

This PR documents the exact check names and explains when maintainers use them. The examples remain available after merge.

\`\`\`text
steady-state CI topology
\`\`\`

The literal identifier \`admission-policy cutover\` is preserved for lookup.`;
  assert.deepEqual(validateBody(body), []);
});

test('a missing or placeholder plain-language summary fails with a specific fix', async (t) => {
  await t.test('missing section', () => {
    assert.match(validateBody('## Change class\n\n- [x] STANDARD')[0], /missing/);
  });
  await t.test('TODO placeholder', () => {
    assert.match(validateBody('## Plain-language summary\n\nTODO')[0], /placeholder/);
  });
  await t.test('comment-only template', () => {
    assert.match(validateBody('## Plain-language summary\n\n<!-- Add summary text. -->')[0], /placeholder/);
  });
});

test('gh pr create --fill is denied before missing explicit options are reported', () => {
  const decision = evaluateHookCommand('gh pr create --fill');
  assert.equal(decision.allowed, false);
  assert.match(decision.reason, /--fill hides/);
  assert.match(decision.reason, /--title/);
  assert.match(decision.reason, /--body-file/);
});

test('gh pr create with a clear literal title and readable body file is allowed', () => {
  const decision = evaluateHookCommand(
    `gh pr create --title '${CLEAR_TITLE}' --body-file '${clearBodyFile}'`
  );
  assert.deepEqual(decision, { allowed: true });
});

test('stdin and unreadable body files are denied with actionable reasons', async (t) => {
  await t.test('stdin body', () => {
    const decision = evaluateHookCommand(
      `gh pr create --title '${CLEAR_TITLE}' --body-file -`
    );
    assert.equal(decision.allowed, false);
    assert.match(decision.reason, /cannot be inspected/);
    assert.match(decision.reason, /Save the complete PR body/);
  });
  await t.test('missing file', () => {
    const missing = join(temporaryDirectory, 'missing.md');
    const decision = evaluateHookCommand(
      `gh pr create --title '${CLEAR_TITLE}' --body-file '${missing}'`
    );
    assert.equal(decision.allowed, false);
    assert.match(decision.reason, /Cannot read PR body file/);
    assert.match(decision.reason, /run the PR command separately/);
  });
});

test('gh pr edit without a wording change is allowed', () => {
  assert.deepEqual(
    evaluateHookCommand('gh pr edit 394 --add-label governance'),
    { allowed: true }
  );
});

test('wording-changing gh pr edit commands are validated', async (t) => {
  await t.test('clear title and body pass', () => {
    assert.deepEqual(
      evaluateHookCommand(
        `gh pr edit 394 --title '${CLEAR_TITLE}' --body-file '${clearBodyFile}'`
      ),
      { allowed: true }
    );
  });
  await t.test('opaque title fails', () => {
    const decision = evaluateHookCommand(`gh pr edit 394 --title '${ORIGINAL_TITLE}'`);
    assert.equal(decision.allowed, false);
    assert.match(decision.reason, /steady-state/);
  });
  await t.test('literal body is checked', () => {
    const decision = evaluateHookCommand(
      `gh pr edit 394 --body '## Plain-language summary

This completes the release admission cutover.'`
    );
    assert.equal(decision.allowed, false);
    assert.match(decision.reason, /release admission/);
  });
});

test('unrelated Bash commands and quoted documentation are allowed', () => {
  assert.deepEqual(evaluateHookCommand('git status --short'), { allowed: true });
  assert.deepEqual(
    evaluateHookCommand(`printf '%s\\n' 'gh pr create --fill'`),
    { allowed: true }
  );
});

test('unresolved wording and compound PR commands are denied', async (t) => {
  await t.test('shell variable title', () => {
    const decision = evaluateHookCommand(
      `gh pr create --title "$TITLE" --body-file '${clearBodyFile}'`
    );
    assert.equal(decision.allowed, false);
    assert.match(decision.reason, /unresolved shell expansion/);
  });
  await t.test('compound command', () => {
    const decision = evaluateHookCommand(
      `git status --short && gh pr create --title '${CLEAR_TITLE}' --body-file '${clearBodyFile}'`
    );
    assert.equal(decision.allowed, false);
    assert.match(decision.reason, /separate command/);
  });
});

test('hook denial output uses the current PreToolUse JSON shape', () => {
  const reason = 'Supply explicit wording.';
  assert.deepEqual(hookDenialOutput(reason), {
    hookSpecificOutput: {
      hookEventName: 'PreToolUse',
      permissionDecision: 'deny',
      permissionDecisionReason: reason
    }
  });

  const result = spawnSync(process.execPath, [CHECKER, '--hook'], {
    input: hookEvent('gh pr create --fill'),
    encoding: 'utf8'
  });
  assert.equal(result.status, 0);
  assert.equal(result.stderr, '');
  const output = JSON.parse(result.stdout);
  assert.equal(output.hookSpecificOutput.hookEventName, 'PreToolUse');
  assert.equal(output.hookSpecificOutput.permissionDecision, 'deny');
  assert.match(output.hookSpecificOutput.permissionDecisionReason, /--fill/);
});

test('hook mode exits silently for allowed commands', () => {
  const result = spawnSync(process.execPath, [CHECKER, '--hook'], {
    input: hookEvent(`gh pr create --title '${CLEAR_TITLE}' --body-file '${clearBodyFile}'`),
    encoding: 'utf8'
  });
  assert.equal(result.status, 0);
  assert.equal(result.stdout, '');
  assert.equal(result.stderr, '');
});

test('manual mode reports bad wording and accepts a clear proposal', () => {
  const bad = spawnSync(process.execPath, [
    CHECKER,
    '--title', ORIGINAL_TITLE,
    '--body-file', clearBodyFile
  ], { encoding: 'utf8' });
  assert.equal(bad.status, 1);
  assert.match(bad.stderr, /PR language check failed/);
  assert.match(bad.stderr, /steady-state/);

  const clear = spawnSync(process.execPath, [
    CHECKER,
    '--title', CLEAR_TITLE,
    '--body-file', clearBodyFile
  ], { encoding: 'utf8' });
  assert.equal(clear.status, 0);
  assert.equal(clear.stdout, 'PR language check passed.\n');
  assert.equal(clear.stderr, '');
  assert.deepEqual(validateProposal({ title: CLEAR_TITLE, body: CLEAR_BODY }), []);
});
