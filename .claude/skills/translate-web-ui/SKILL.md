---
name: translate-web-ui
description: Translate the new and changed gbdraw Web UI keys into the ten locale catalogs in one batch before a release - status, English-source fixes, term audit, per-locale drafting and review, apply, and one PR. Use before a release tag when `node tools/i18n/sync.mjs status --require-complete` fails, or when the Owner asks for the translation batch.
---

# Translate the Web UI before a release

Feature PRs add or change English only and leave the catalogs alone: an
edited or removed text stays in them as a stale key, and the new text ships in
English until this batch fills the catalogs
(`gbdraw/web/js/utils/i18n-catalogs/<locale>.js`). Pull-request and `dev` runs
tolerate stale and missing keys; the release tier fails while a catalog lacks a
key or holds a stale one. `export` offers a stale key's old translation as
`previous` (`needs-check`) for a new key whose English is a small edit of it,
and `apply` drops the stale keys. Every
translation stays an LLM draft: never record, show or write that a
translation is reviewed, checked by a native speaker, or final.

Authority: `docs/internal/I18N_STYLE.md` (§1 for all locales, then one section
per locale), `docs/internal/I18N_GLOSSARY.tsv` (the glossary decides words),
and the lint lists `tests/web/i18n-lint/<locale>.tsv`.

Models: Sonnet runs steps 1, 6 and 7 (commands, apply, tests). Opus does
steps 2 to 5 (English fixes, terms, drafting, review, decisions).

1. **Status.** In a worktree from `origin/dev`, run
   `node tools/i18n/sync.mjs status`. Note the missing and stale keys per locale.
2. **Fix the English first.** Read every missing key with its call sites
   (`node tools/i18n/sync.mjs export --locale ja` lists them with roles). Fix
   the source before anyone drafts:
   - a key shared by two meanings gets `tx(context, text)`;
   - a sentence built from fragments, or a control, track, mode or setting name
     inserted as `{placeholder}`, becomes one whole-sentence key per case;
   - a sentence that counts is a `tn()` pair, never a plain `t()` with `(s)`;
   - the English names only what the UI shows.
   Run `node --test tests/web/i18n-catalogs.test.mjs` after the fixes.
3. **Terms.** List the terms of the new keys that the glossary lacks. For each
   one and each locale, find sources (vendor terminology, field papers or
   theses, standard references) before drafting, and add the glossary row with
   the source in its note. Do not coin a term: keep the English word when no
   established target-language term exists. A term with weak evidence goes to
   the Owner list in the PR body.
4. **Draft per locale.** For each locale, write
   `node tools/i18n/sync.mjs export --locale <x> > $TASK_DIR/<x>.jsonl` and
   give one drafter that file, I18N_STYLE §1, that locale's section, and the
   glossary rows in the file. A `previous` entry marked `needs-check` is the
   old translation of a key whose English changed slightly: reuse it only
   after checking it against the new English. The drafter returns JSON lines
   `{"key": ..., "translation": ...}`; a drafter of a language with `few` or
   `many` forms returns the `one` key of a `plural` pair as the string for
   exactly 1 and its `other` key as an object of CLDR categories with `other`
   (`one` there is for 21, 31, …; I18N_STYLE §1.3).
5. **Review per locale, in a fresh context.** Give a reviewer the English,
   the draft, the style section and the glossary, and these calibration cases
   from earlier rounds, which drafts got wrong:
   - an option or button name used as the subject of a sentence (write
     "With Auto, ..." / "When you save, ...": es Con Automático, …; pt Ao
     salvar a sessão, …; fr Avec « Enregistrer », …);
   - Linear / Circular read as molecule topology when they name the diagram
     mode (es registro lineal, fr séquence linéaire, pt registro linear are
     wrong for "record in the Linear diagram"; Circular record as topology is
     right);
   - option labels whose meaning drifted: Either meant "either direction",
     not "one of the two"; File order meant "order within the file";
   - a count passed to a plain `t()`: the sentence must read right for 1 too
     (or the English is fixed to `tn()` in step 2);
   - coined terms (ko 피처 for feature, es Discordancias for mismatches);
   - a template whose placeholder takes another key's text, read with that text
     in place, also when it is a verb phrase (`Undo {action}`: ja
     ファイルを削除を元に戻す, es Deshacer Quitar archivo; a colon or quotes
     fix the class, I18N_STYLE §1.3).
   The reviewer lists proposals with a reason and a severity, not a rewrite.
6. **Decide.** The coordinator decides each proposal against the style and the
   glossary; there is no majority vote. Record terms that remain doubtful.
7. **Apply and check.** Per locale:
   `node tools/i18n/sync.mjs apply --locale <x> <decided.jsonl>` (it refuses
   unknown keys and placeholder changes and drops stale keys). Then run
   `node --test tests/web/i18n-catalogs.test.mjs tests/web/i18n-tools.test.mjs`
   (catalog keys, placeholders, lint `error` rules) and read its `warn` lines.
   Finish with `node tools/i18n/sync.mjs status --require-complete`.
8. **One PR** into `dev` with the catalogs, glossary rows, STYLE changes and
   English-source fixes. The body lists the terms kept on weak evidence and
   the open questions for the Owner; it never calls the translations reviewed.
   After the merge the translation-review workflow republishes the review page.
