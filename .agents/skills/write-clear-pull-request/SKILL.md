---
name: write-clear-pull-request
description: Write or review gbdraw PR titles and descriptions for maintainers. Applies to PR wording, not ordinary code review or commit messages.
---

# Write clear pull requests

Assume the reader knows gbdraw and GitHub, but has not read the implementation
plan or private discussion. Name the concrete action and object in the title.
Begin the body with `## Plain-language summary`: two to four sentences explaining
what changes, why, and the behavior or review condition after merge. Preserve
exact identifiers in backticks and explain unfamiliar ones. Internal process
labels belong after the concrete explanation.

For example, replace `Promote steady-state CI topology` with
`Stop rerunning full CI on dev-to-main pull requests`.

Include validation and material limitations in proportion to the change. Use the
final scope; omit abandoned approaches unless they explain a remaining tradeoff.
Save the complete body and check the exact title/body before creating the PR:

```bash
node tools/check-pr-language.mjs --title "<title>" --body-file <body.md>
```

Run the checker once before `gh pr create` and again only after a material wording
rewrite before merge. Run any authorized `gh pr create` or wording-changing
`gh pr edit` separately. Do not use `--fill`, `--fill-first`, or `--fill-verbose`.
Do not repeatedly edit PR metadata just to append run IDs. This skill supplies
wording guidance; it does not authorize publication, comments, or merge.
