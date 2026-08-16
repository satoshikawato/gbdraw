# gbdraw Web refactor guardrail bundle

This bundle contains the proposed repository assets for preventing
behavior-preserving Web refactors from introducing semantic, ownership, race, or
large-copy regressions.

## Included files

```text
.agents/skills/refactor-gbdraw-web-safely/SKILL.md
AGENTS.md
.claude/settings.json
.claude/hooks/run-web-refactor-guards.sh
.githooks/pre-push
tools/install-git-hooks.sh
tests/web/refactor-guardrails.test.mjs
package.json
.github/workflows/web-refactor-guardrails.yml
.github/pull_request_template.md
patches/CLAUDE.md.addition.md
patches/gbdraw-web-CLAUDE.md.addition.md
patches/change-gbdraw-rendering-surface-SKILL.md.addition.md
```

## Apply

Review the diff before copying because `AGENTS.md`, `package.json`, and
`.claude/settings.json` are complete replacement files based on the repository
state inspected on 2026-08-15.

From the repository root:

```bash
# Copy the bundle contents into the checkout, preserving hidden directories.
cp -R /path/to/gbdraw-web-refactor-guardrails/. .

chmod +x .claude/hooks/run-web-refactor-guards.sh
chmod +x .githooks/pre-push
chmod +x tools/install-git-hooks.sh

./tools/install-git-hooks.sh
npm run test:web:refactor-guards
```

Apply the three files under `patches/` manually to the existing long guidance
files.

## Enforcement model

```text
SKILL.md
    tells agents how to perform the refactor

AGENTS.md / CLAUDE.md
    trigger the Skill and state durable invariants

Claude Code Stop hook
    gives fast feedback before a Claude turn completes

Git pre-push hook
    blocks a local push when guardrails fail

GitHub Actions
    is the repository-wide enforcement authority
```

The local hooks are convenience and early feedback. CI remains the final
cross-agent, cross-machine gate.

## Hook behavior

The Claude Code project hook runs only when relevant Web or guardrail files are
modified. It blocks one Stop attempt when the guard command fails and returns the
failure to Claude. It avoids an unbounded Stop-hook loop on the subsequent
continuation.

The Git hook runs the same command before pushing relevant changes.

Neither hook installs dependencies or modifies production files.
