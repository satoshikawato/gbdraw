---
name: write-clear-pull-request
description: Write, edit, or review gbdraw pull request titles and maintainer-facing descriptions. Use whenever creating or changing a PR title or body, or reviewing either for clarity.
---

# Write clear pull requests

## Audience

Assume the reader knows gbdraw, Git, and GitHub, but has not read an implementation plan, earlier PRs, or private discussion.

## Required opening

- Write the title as a concrete action plus a concrete object.
- Begin the body with `## Plain-language summary` and use two to four sentences to say what changes, why, and what differs after merge.
- Treat internal classifications and process labels as metadata, not as the primary explanation.
- Preserve exact check names, paths, branch names, options, and code identifiers in backticks. Explain their role in ordinary language when the name is not self-explanatory.

Policy terms may appear in later technical or evidence sections. They fail this standard only when they replace the concrete explanation in the title or opening summary.

## Workflow

1. Name the concrete action and object.
2. State the reason when it is not obvious.
3. State the behavior or review condition after merge.
4. Replace unexplained local abstractions in the title and summary.
5. Check that the text makes sense without prior plans or PR history.
6. Save the complete body, then run:

   ```bash
   node tools/check-pr-language.mjs --title "<concrete title>" --body-file <body.md>
   ```

Run `gh pr create` or a wording-changing `gh pr edit` separately after the checker passes.

## Regression example

Bad title:

```text
Promote steady-state CI topology and finalize main admission
```

Bad summary:

```text
This promotion removes the transition-only CI topology from main and completes the admission-policy cutover.
```

Clear title:

```text
Stop rerunning full CI on dev-to-main pull requests
```

Clear summary:

```text
This PR removes CI jobs that were needed only while the new main-merge check was being introduced. After the branch-protection update, dev-to-main pull requests will require only `Promotion / gate` and `CodeQL`.
```
