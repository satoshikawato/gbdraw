# translation-review

An orphan branch that holds the gbdraw Web UI translation review page. It is
not part of the gbdraw package and never merges into `dev` or `main`.

GitHub Pages serves this branch at https://satoshikawato.github.io/gbdraw/.
The page is at the root of the branch (`.nojekyll` keeps Pages from processing
the files).

## Layout

- `index.html` - the review page.
- `data/` - what the page loads: `locales.json`, `keys.json`, `glossary.json`, `shots.json`,
  `style.md`, and `draft-<locale>.json` and `notes-<locale>.json` per locale.
- `shots/` - the Web UI screenshots (WebP) that `data/shots.json` points at.
- `db/current/` - the review database the page reads
  (`rev-<locale>/c<hex>.json`; the `batch-<locale>.json` files seed a database).

## Only the latest page is published

`data/` is generated from `dev`: the `Translation review page` workflow
(`.github/workflows/translation-review.yml` on `dev`) runs
`node tools/i18n/review-page-data.mjs --out data` after each `dev` push that
changes the Web UI keys, the catalogs, the glossary or the style guide, and
commits only when the data changed. Do not edit `data/` by hand; the next run
replaces it. `data/shots.json` keeps its screenshot rectangles and loses the
ids of keys that `dev` no longer has.

A new version replaces `index.html` and `data/` in place. There are no version
folders. Older versions exist only in the branch history; the last commit that
still had `review-gui/` (v19 and v20) is `fb3371073`.

On Pages the page is read-only: it shows the drafts, the notes and the review
state stored in `db/current/` as of the branch's latest commit. Edits are not
saved there.
