---
name: love-me-love-my-docs
description: Create or revise procedural user documentation with reproducible screenshots, commands, or code examples. Use for workflow evidence and manual renovation, not simple wording fixes or internal plans.
---

# Reproducible user documentation

Every changed visible result must be reproducible, and every documented action
must be performed on the surface being taught. Preserve valid existing evidence
for unchanged flows. Fix a failed recipe or report the concrete missing evidence;
do not replace it with an unsupported claim.

## Scope and references

Use the requested audience, language, surfaces, output format, and existing docs
system. Infer these from the repository when clear. A small update should reuse
its page and harness; it does not require a new site, full flow census, progress
checklist, login, or demo database.

Read only the relevant guidance:

- New or reorganized public pages: [manual-structure.md](references/manual-structure.md).
- Web workflow proof: [capture-web.md](references/capture-web.md).
- CLI commands or literal Python examples: [execute-cli-python.md](references/execute-cli-python.md).
- Explicitly requested mobile workflows: [capture-mobile.md](references/capture-mobile.md).
- Sequence inputs, mirrored fixtures, or sequence-derived visible results:
  [sequence-inputs.md](references/sequence-inputs.md).

For gbdraw Gallery operation crops and tutorial JSON, use
`web-gallery-screenshot-maintenance`; use this skill as well only when the request
also needs fresh-input procedural evidence or a broader manual renovation.

## Page ownership

For a new page or structural renovation, record the reader question, existing
owner, supporting scenarios/surfaces, and `keep`, `merge`, `delete`, or `new`
disposition. Prefer an existing owner when it answers the question clearly.
Evidence for several interfaces can support one page; evidence inventory is not
page inventory. Ask only when an unresolved audience, scope, or product choice
would materially change the result. An already authorized page/harness change
needs no separate confirmation gate.

## Evidence

Reuse the repository's capture/recipe tools and artifact layout. If new automation
is needed, prove one representative workflow before scaling it to the remaining
requested scenarios. Commit reproducible harness changes with their outputs.

- Web screenshots must come from real UI actions, with state assertions before
  capture and validated downloads. Do not inject a finished state as proof of
  the workflow. A session-reload lesson uses a session created earlier in that
  same lesson from original inputs.
- CLI commands run from a clean temporary directory with declared inputs.
  Execute published Python blocks unchanged, not equivalent test-only programs.
- Use public, verified or synthetic demonstration inputs. Keep private data,
  credentials, and auth storage out of public images and committed artifacts.
- For applications that persist data remotely, identify the actual mutations
  before seeding or clicking. Local disposable fixtures need no extra approval.
  Read-only public pages and browser-local rendering need no blanket host
  approval. Obtain authorization for real external mutations if not already
  covered by the request; a local hostname alone does not prove isolation.

A missing browser or broken workflow is a named blocker, not successful evidence.
Leave an executable partial harness or capture checklist only when useful, keep
incomplete public output clearly marked, and finish unaffected work. Do not
replace valid public media with placeholders merely because a refresh is blocked.

## Completion

Run the changed scenarios, validate their named outputs and media links, and
inspect rendered pages/figures at a readable scale. Include the inputs, actual
actions, expected result, and useful troubleshooting in the page. Document and
verify the regeneration command in the existing evidence owner.

Choose checks that address the changed source of drift: selectors, input hashes,
locale/viewport, authentication when present, or command/API behavior. Reuse
unchanged evidence and avoid rerunning unrelated flows. Report what changed,
which scenarios ran, and any unmet evidence requirement; do not call an incomplete
requested workflow finished.
