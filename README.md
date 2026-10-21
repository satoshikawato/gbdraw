# translation-review

An orphan branch that keeps the gbdraw Web UI translation review page and its
review database in the repository, so the review can be republished or restored
without the claude.ai artifact. It is not part of the gbdraw package and never
merges into `dev` or `main`.

## Contents

- `review-gui/v16/` — the page (`index.html`) and the files it loads:
  `data/` (string keys with context, glossary, the ja / zh-Hans / zh-Hant drafts
  and draft notes, the screenshot index, the style rules) and `shots/`
  (114 screenshots of the Web UI, WebP).
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
| https://claude.ai/artifact/2wf2xJ42bpyxwMxXH7prmX | the account used by Claude Code on the Owner's machine | copy of v16, database seeded from `db/seed-2026-10-21/` |

## Republish the page

From Claude Code, with the Artifact tool: publish `review-gui/v16/index.html`
with `root` = `review-gui/v16`, `files` = every path under `data/` and
`shots/` (124 files, 12 MB), and
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

- Page: `review/index.html` (v16) of the 2026-10-19 i18n campaign, published
  2026-10-21; its data files were generated from `dev` at that time.
- `data/draft-*.json` are the drafts the review started from, not the shipped
  catalogs; the shipped catalogs are `gbdraw/web/js/utils/i18n-catalogs/` on `dev`.
- A row's `ok` status records that a reviewer accepted it in the page. This
  branch stores that record; it makes no statement of its own about the quality
  of a translation.
