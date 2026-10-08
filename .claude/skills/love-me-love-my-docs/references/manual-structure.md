# Manual structure and public-page ownership

Plan public pages around distinct reader questions and existing owners, not
around capability, workflow, evidence-scenario, or interface counts. Evidence
scenarios are internal verification units: one or many scenarios may support a
single public page, and a scenario may support no public page. Evidence for a
GUI, CLI, Python, or mobile surface validates that surface without requiring a
separate surface-specific page.

## Public-page decision gate

Before drafting, record one row for every reader question or affected page:

| Reader question | Existing owner | Evidence scenarios and surfaces | Disposition | Resulting owner |
| --- | --- | --- | --- | --- |

Use only these dispositions:

- `keep`: the existing page owns and answers the question clearly;
- `merge`: move useful material into an existing owner and remove duplication;
- `delete`: remove a page that has no distinct reader question after preserving
  unique content and executable evidence in the proper owners;
- `new`: create a page only when editing an existing owner cannot answer the
  distinct question clearly.

Default to `merge`. Do not create one page per capability, user flow, evidence
scenario, or interface. A separate surface page needs a materially different
reader journey, not merely a separate capture or execution harness.

## Public page patterns

- **Tutorials:** Teach a deliberate learning progression to a finished result.
  A Tutorial may be foundational or advanced. Teach one successful route and
  link optional variations instead of branching into an option inventory.
- **Technical documentation:** Own exact behavior, inputs, controls, options,
  schemas, API contracts, sessions, outputs, errors, compatibility, and limits.
  Use compact operational examples only when they clarify the contract.
- **FAQ:** Answer choices and troubleshooting questions concisely. Link exact
  values and full contracts to their Technical documentation owner instead of
  repeating them.
- **Gallery:** Help readers discover finished outcomes. Link to a reproducible
  Tutorial when one exists; do not make Gallery entries own procedures or
  technical contracts.

For a Tutorial, include a concrete outcome and finished result; a file table
that distinguishes **Download**, **Create**, **Generated**, and **Reference
result** files; the fresh starting state or empty working directory; numbered
actions with an early visible result; expected intermediate and final results;
useful screenshots or generated figures; verification cues; troubleshooting;
and links to the exact technical owner. Present surface variants together when
that serves the same reader journey; split them only when their journeys are
materially different.

## File and result contract

For every procedural page:

- direct sequence acquisition to the authoritative public database by
  accession, format, and exact save name;
- pin an official annotation revision or release when the visible result
  depends on a mutable feature table; a nucleotide accession version alone is
  not sufficient evidence for annotation-sensitive output;
- never tell a reader to use a repository-bundled sequence, finished Gallery
  session, or other prebuilt project state as a public procedural input;
- use session reload only after the reader created that session from original
  inputs in the same session-focused procedural page;
- keep frozen repository copies internal to deterministic offline verification
  and apply the mirror-verification gate in Step 3 of `SKILL.md`; a
  `legacy-unverified` mirror cannot back a public visible result;
- link each non-sequence supplied input to its source and state its exact save
  name, or show its complete contents or derivation on the page;
- show complete contents for every file the reader must create;
- name outputs before the action that generates them;
- show the expected directory tree when several files interact;
- state the command, button sequence, or literal program that produces the
  result;
- show a generated result and concrete checks, not merely "the command
  succeeds";
- keep reference artifacts visibly separate from inputs.

Use tables for repeated fields such as file mappings, record metadata, track
slots, color rules, and expected outputs.
