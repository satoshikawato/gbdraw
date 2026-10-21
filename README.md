# translation-review

An orphan branch that keeps the gbdraw Web UI translation review page and its
review database in the repository, so the review can be republished or restored
without the claude.ai artifact. It is not part of the gbdraw package and never
merges into `dev` or `main`.

## Contents

- `review-gui/v18/` — the page (`index.html`) and the files it loads:
  `data/` (string keys with context, the glossary with ja, zh-Hans, zh-Hant,
  ko, es, fr and pt columns, a draft and a draft-notes file per locale, the
  screenshot index, the style rules with one section per locale) and `shots/`
  (114 screenshots of the Web UI, WebP). v18 reviews seven locales:
  日本語, 简体中文, 繁體中文, Español, Français, Português, 한국어.
- `db/snapshot-2026-10-21-r4/` — the review database of the original artifact
  as read on 2026-10-21 after review round 4: collections `rev-ja`,
  `rev-zh-Hans`, `rev-zh-Hant`, 16 chunk documents each
  (`{"items": {<key id>: {t, s, c, by, at}}}`; `s` is `ok` for an accepted row,
  `flag` for a flagged one; `c` carries the round's reason or the flag comment).
- `db/d46-rows/` — the rows for the strings added by decision D-46 of the
  2026-10-19 i18n campaign (PR #1254), which are not in the round-4 snapshot.
- `db/seed-2026-10-21/` — snapshot plus D-46 rows: the state written into the
  account-owned copy of the page on 2026-10-21 (48 documents, 1,388 items).
  `tools/build_seed.py` reproduces it from the two folders above.

## Artifacts

| Page | Owner | State on 2026-10-21 |
| --- | --- | --- |
| https://claude.ai/artifact/DywpUqrJtpb3tS3EvwAxEj | the Owner's other claude.ai account | original; edits made in the page after round 4 are not in this branch (use the page's "Export for an agent review…" to get them) |
| https://claude.ai/artifact/2wf2xJ42bpyxwMxXH7prmX | the account used by Claude Code on the Owner's machine | v18 (Version 2, seven locales); ja / zh database seeded from `db/seed-2026-10-21/`; the four new locales start with no review state |

## Republish the page

From Claude Code, with the Artifact tool: publish `review-gui/v18/index.html`
with `root` = `review-gui/v18`, `files` = every path under `data/` and
`shots/` (132 files, 12 MB), and
`capabilities` = `{"db": {}, "user": {"scopes": ["profile"]}, "downloads": true}`.
Pass `url` to update an existing artifact; omit it to create a new one.
The page shows a database-less, read-only view when opened outside claude.ai.

## Seed or restore the database

```bash
python3 -I tools/build_seed.py db/snapshot-2026-10-21-r4 db/d46-rows <out dir> [versions.json]
```

writes `<out dir>/rev-<locale>/<chunk>.json` and `<out dir>/batch-<locale>.json`.
Then run one ArtifactData `batch` per locale with the 16 `set` entries, giving
each entry's `file_path` as an absolute path. For a database that already holds
the documents, first `list` each collection and pass their versions as
`versions.json` (`{"rev-ja": {"c0": 3, ...}, ...}`) so every write is pinned.

## Provenance

- Page: `review/index.html` (v16/v17) of the 2026-10-19 i18n campaign, extended
  on 2026-10-21 with the four locales of v0.15.0 (v18).
- Keys: the strings of the Web UI on branch `feat/web-i18n-0150-locales`
  (commit 7f0657edf, on top of `dev` 8fec90d4f), 2,682 including the 50 singular
  forms of the plural strings (`sense` names them).
- `data/draft-ja|zh-Hans|zh-Hant.json` are the drafts the ja / zh review started
  from (strings added since then carry the catalog text); `data/draft-ko|es|fr|pt.json`
  are the catalogs of that branch. The shipped catalogs are
  `gbdraw/web/js/utils/i18n-catalogs/` on `dev`.
- `data/notes-ko|es|fr|pt.json` say, per string, whether the independent review
  of the draft changed it (with the text before) and whether the code branch
  changed it after the review.
- A row's `ok` status records that a reviewer accepted it in the page. This
  branch stores that record; it makes no statement of its own about the quality
  of a translation.
