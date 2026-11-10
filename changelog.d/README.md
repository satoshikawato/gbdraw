# Changelog fragments

A pull request with a change worth a `CHANGELOG.md` entry writes the entry
here, in a file of its own, instead of adding it to the `[Unreleased]` list.
Two pull requests then never edit the same lines, so neither has to rebase
because the other added an entry first.

- Name the file after your branch, with `/` replaced by `-`, for example
  `fix-legend-side-order.md` for `fix/legend-side-order`. Use lower-case
  letters, digits, `.`, `_`, and `-`.
- Write the entries exactly as they will read in `CHANGELOG.md`: each starts
  with `- Added: `, `- Changed: `, `- Deprecated: `, `- Removed: `,
  `- Fixed: `, `- Security: `, `- Performance: `, or `- Known issue: `, and
  continues on lines indented by two spaces. One file may hold several
  entries; leave no blank lines.
- Write links as they resolve from `CHANGELOG.md` (the repository root), for
  example `./docs/REFERENCE/web-app.md`.
- `python tools/changelog.py check` reports what does not follow these rules;
  `tests/test_changelog_fragments.py` runs the same check.

Before a release, the release session runs `python tools/changelog.py fold`.
It puts the entries of every fragment at the top of the `[Unreleased]` list of
`CHANGELOG.md`, newest merge first, and deletes the fragments.
