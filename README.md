# translation-review

An orphan branch that keeps the gbdraw Web UI translation review page and its
review database in the repository, so the review can be republished or restored
without the claude.ai artifact. It is not part of the gbdraw package and never
merges into `dev` or `main`.

GitHub Pages serves this branch: https://satoshikawato.github.io/gbdraw/ opens
the review page (`index.html` at the root forwards to `review-gui/v20/`;
`.nojekyll` keeps Pages from processing the files). There the page is
read-only: it shows the drafts, the notes and the review state stored in
`db/current/` (status, comments, stored translations, as of the branch's latest
commit), and "Export for an agent review…" downloads a file that carries that
state. Edits are not saved (saving, reviewer names and "Import" persistence
need the claude.ai runtime). When the claude.ai runtime is present the page
uses its database and does not read `db/current/`.

## Contents

- `review-gui/v19/` — the page (`index.html`) and the files it loads:
  `data/` (string keys with context, the glossary with ja, zh-Hans, zh-Hant,
  ko, es, fr and pt columns, a draft and a draft-notes file per locale, the
  screenshot index, the style rules with one section per locale) and `shots/`
  (114 screenshots of the Web UI, WebP). v19 reviews seven locales:
  日本語, 简体中文, 繁體中文, Español, Français, Português, 한국어, with 2,735
  strings: v18's 2,682 plus the 53 Custom Track Slots issue messages of #1259
  (`dev` d17c0b7a8). ja, zh-Hans and zh-Hant have the catalog text as their
  draft; es, fr, pt and ko have no draft for these 53 yet (PR-B #1261).
- `review-gui/v20/` — v19's page with ten locales: 日本語, 简体中文, 繁體中文,
  한국어, Español, Français, Português, Bahasa Indonesia, Tiếng Việt, ไทย, and
  2,738 strings (`dev` 36ac5879a: every catalog key, including the 50 singular
  forms; v19's 2,735 minus 2 keys `dev` no longer has plus 5 new ones). The drafts of the first seven
  locales are the `dev` catalogs; id, vi and th come from PR #1265
  (`feat/web-i18n-thidvi`, 3173c54f2), with a note per string that the
  independent review changed (with the text before), 66, 31 and 50. The
  screenshots are shared with v19 (`shots.json` points at `../v19/shots/`), and
  `data/style.md` and `data/glossary.json` carry the id, vi and th sections and
  columns. The three new locales have no review database rows: the page shows
  them read-only with the drafts only. The v19 paragraph above stays as it was
  for that folder.
- `db/snapshot-2026-10-21-r4/` — the review database of the original artifact
  as read on 2026-10-21 after review round 4: collections `rev-ja`,
  `rev-zh-Hans`, `rev-zh-Hant`, 16 chunk documents each
  (`{"items": {<key id>: {t, c, by, at}}}`). No row has an `s` field: nobody
  accepted or flagged a row in the page. `t` is the stored translation, `c` the
  round's reason (`[R2 …]`, `[R3 …]`, `[R4 …]`), and a non-empty `by` marks an
  edit made in the page by a person.
- `db/d46-rows/` — the rows for the strings added by decision D-46 of the
  2026-10-19 i18n campaign (PR #1254), which are not in the round-4 snapshot.
- `db/fu-1258-rows/` — 10 rows per locale for the strings changed by the
  i18n follow-ups (PR #1258; comment starts with `[FU]`).
- `db/ts-1259-rows/` — 53 rows per locale (ja, zh-Hans, zh-Hant) for the strings
  of #1259: the catalog text, comment `[TS] #1259 (2026-10-21): <English>`.
- `db/r1-1266-rows/` — round-1 rows of PR #1266 for ko, es, fr and pt (83 / 231 / 132 / 193;
  comment starts with `[R1] #1266`).
- `db/seed-2026-10-21/` — snapshot plus D-46 rows: the state written into the
  account-owned copy of the page on 2026-10-21 (48 documents, 1,388 items).
  Kept as a record of that write; `tools/build_seed.py` reproduces it from the
  snapshot and `d46-rows`.
- `db/current/` — snapshot plus `d46-rows`, `fu-1258-rows`, `ts-1259-rows` and `r1-1266-rows`
  (112 documents, 2,216 items; 22, 18 and 26 snapshot rows belong to keys that no
  longer exist and are never shown). The page on GitHub Pages reads
  `db/current/rev-<locale>/<chunk>.json`; the `batch-<locale>.json` files seed
  a database with it. Rows written in the claude.ai copy after 2026-10-21 are
  not here until someone exports them and adds a snapshot.

## Artifacts

| Page | Owner | State on 2026-10-21 |
| --- | --- | --- |
| https://claude.ai/artifact/DywpUqrJtpb3tS3EvwAxEj | the Owner's other claude.ai account | original; edits made in the page after round 4 are not in this branch (use the page's "Export for an agent review…" to get them) |
| https://claude.ai/artifact/2wf2xJ42bpyxwMxXH7prmX | the account used by Claude Code on the Owner's machine | v18 (Version 3, seven locales); ja / zh database seeded from `db/seed-2026-10-21/`; the four new locales start with no review state |

## Republish the page

From Claude Code, with the Artifact tool: publish `review-gui/v19/index.html`
with `root` = `review-gui/v19`, `files` = every path under `data/` and
`shots/` (132 files, 12 MB), and
`capabilities` = `{"db": {}, "user": {"scopes": ["profile"]}, "downloads": true}`.
Pass `url` to update an existing artifact; omit it to create a new one.
Outside claude.ai the page shows the read-only view described above.

## Seed or restore the database

```bash
python3 -I tools/build_seed.py [--versions versions.json] db/snapshot-2026-10-21-r4 <out dir> \
    db/d46-rows db/fu-1258-rows db/ts-1259-rows db/r1-1266-rows
```

applies the rows folders in the order given (a later row wins; a row whose `by`
is set, edited in the page, is never overwritten) and writes `<out dir>/rev-<locale>/<chunk>.json` and `<out dir>/batch-<locale>.json`
(`db/current/` is this output).
Then run one ArtifactData `batch` per locale with the 16 `set` entries, giving
each entry's `file_path` as an absolute path. For a database that already holds
the documents, first `list` each collection and pass their versions as
`versions.json` (`{"rev-ja": {"c0": 3, ...}, ...}`) so every write is pinned.

## Provenance

- Page: `review/index.html` (v16/v17) of the 2026-10-19 i18n campaign, extended
  on 2026-10-21 with the four locales of v0.15.0 (v18) and the 53 keys of #1259 (v19).
- Keys: the strings of the Web UI on branch `feat/web-i18n-0150-locales`
  (commit 7f0657edf, on top of `dev` 8fec90d4f), 2,682 including the 50 singular
  forms of the plural strings (`sense` names them), plus the 53 keys of #1259
  from `dev` d17c0b7a8 (2,735). Three v18 keys that `dev` no longer has and the
  50 singular forms stay: they belong to the ko / es / fr / pt branch.
- `data/draft-ja|zh-Hans|zh-Hant.json` are the drafts the ja / zh review started
  from (strings added since then carry the catalog text); `data/draft-ko|es|fr|pt.json`
  are the catalogs of that branch. The shipped catalogs are
  `gbdraw/web/js/utils/i18n-catalogs/` on `dev`.
- `data/notes-ko|es|fr|pt.json` say, per string, whether the independent review
  of the draft changed it (with the text before) and whether the code branch
  changed it after the review.
- A row's `s` field (`ok` accepted, `flag` flagged) would record a reviewer's
  decision in the page; no stored row has one yet. This branch makes no
  statement of its own about the quality of a translation.
