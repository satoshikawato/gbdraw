# Manual structure and information roles

Use one public information role per chapter. A topic may appear in more than
one role only when each page answers a different reader question and links to a
single owner for shared facts or procedures.

## Tutorial: learn by completing a project

A Tutorial may be foundational or advanced. Keep advanced, multi-feature
projects as Tutorials when the reader follows a controlled learning sequence
to a finished result.

Include:

1. a concrete outcome and finished result;
2. a surface switch when the same figure has GUI, CLI, and Python variants;
3. a file table that distinguishes **Download**, **Create**, **Generated**, and
   **Reference result** files;
4. the fresh starting state or empty working directory;
5. numbered actions, with the first visible result early in the sequence;
6. expected intermediate and final results;
7. useful screenshots for GUI decisions and generated figures for CLI/Python;
8. verification cues, troubleshooting, and links to focused How-to and exact
   Reference pages.

Teach one successful route. Link optional variations instead of branching the
lesson into an option inventory.

## How-to guide: accomplish one focused task

State the required starting state and exact inputs, then give the shortest
complete procedure for the named goal. Identify the expected output and how to
verify it. Include troubleshooting for known boundary failures. Link to a
Tutorial for foundational learning and to Reference for exhaustive values.

The How-to guide owns the adaptable task procedure. A related Tutorial may use
the task inside a larger project, but should not duplicate its full option
matrix or every variant.

## Explanation: understand or choose

Start with the decision or concept. Explain the mental model, alternatives,
trade-offs, and consequences. Use diagrams or comparison tables when they aid
understanding. Do not own a runnable end-to-end recipe, fixed defaults, schema
history, or exhaustive tokens; link to How-to and Reference.

## Reference: look up exact facts

Organize exact current contracts for scanning: supported values, defaults,
signatures, schemas, formats, errors, compatibility, and limits. Avoid a
numbered "main workflow" or tutorial narrative. Give minimal syntax examples
only when they clarify a contract, and link procedural use to the canonical
How-to guide.

## File and result contract

For every procedural chapter:

- direct sequence acquisition to the authoritative public database by
  accession, format, and exact save name;
- pin an official annotation revision or release when the visible result
  depends on a mutable feature table; a nucleotide accession version alone is
  not sufficient evidence for annotation-sensitive output;
- never tell a reader to use a repository-bundled sequence, finished Gallery
  session, or other prebuilt project state as a Tutorial or How-to input;
- use session reload only after the reader created that session from original
  inputs in the same session-focused procedural chapter;
- keep frozen repository copies internal to deterministic offline verification
  and apply the mirror-verification gate in Step 3 of `SKILL.md`; a
  `legacy-unverified` mirror cannot back a public visible result;
- link each non-sequence supplied input to its source and state its exact save
  name, or show its complete contents/derivation in the chapter;
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
