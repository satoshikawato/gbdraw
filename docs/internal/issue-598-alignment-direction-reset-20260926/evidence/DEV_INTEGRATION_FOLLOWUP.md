# Published S02 follow-up integration

Baseline: `dev@34c5104be196bda9d4a034ef4374128edb837982`, after PR #612.
Source: the subsequently published original S02 commit `400c746b` on
`origin/fix/issue-598-alignment-direction-reset-20260926`. The original
writer and branch were read only. This separate dev-derived branch preserves
PR #612's fresh-load identity, pending intent, compact layout, portable Python
parity checks and plain-text record label fixes.

## Implemented differences

- Python receipt binding now sorts record keys by UTF-16 code units, matching
  JavaScript for supplementary Unicode keys. No schema, historical reader,
  receipt fields or digest algorithm changes.
- Successful canonical candidate adoption clears the prior error and preserved
  failure flag through the existing Result admission owner. Failed candidates
  still retain the previous artifact and editable draft.
- Combined Reset lists the canonical presentation label before source metadata,
  through the existing plain-text formatter; source definition, file definition,
  accession and record key remain the established fallback chain.

These implement the existing accepted alignment/Reset continuation. Formal
Product contract revision 21 and SHA-256
`62e3a9c08ceb64acc349ace81a97a9c87dc45db4d18187a3b78e6e131b2dc595`
remain unchanged. PD-OI-035 revision 3 and PD-OI-039 revision 2 retain
exclusive Keep/right/left/Custom directions without Match. No new Product
outcome, authority update or S06 presentation work is introduced.

## Verification

Before the production fixes, the Unicode Web/Python binding regression and
successful-candidate error clearing regression both failed on dev. After the
fixes:

```sh
node --test tests/web/alignment-reset-receipt.test.mjs tests/web/run-analysis-simple-path.test.mjs tests/web/similarity-alignment-actions.test.mjs tests/web/history.test.mjs
pytest tests/test_api_session.py tests/test_session_compat.py tests/test_session_io.py -v
ruff check gbdraw/
node tools/check-web-change-budget.mjs --base origin/dev
```

- Focused Web owners: **58 passed** (`/tmp/issue598-followup-final-node.log`).
- All affected Python session regression files: **287 passed**
  (`/tmp/issue598-followup-final-python.log`; command above).
- Ruff and whitespace passed. Working-tree policy: **Gate PASS / Review CLEAR**.
- The Reset regression includes a typed presentation label different from the
  source name and validates the existing plain-text formatting.
- The previous 711-test Web, 26-test focused browser and 13-test PR smoke results
  in `DEV_INTEGRATION.md` describe PR #612 only. They are not claimed as new
  whole-suite results for this follow-up. No all-staging or release claim.

## Owners and paths

Python session validation, canonical Result admission and Reset preview remain
one owner and one path each. No additional owner, branch, compatibility reader,
renderer, Worker job, dependency or persisted field is introduced. OE/PE/CB
sets are unchanged; no architecture exception applies. Production and test
diffs were reviewed separately. The generated browser wheel is excluded.

Rollback: revert this follow-up commit; retain PR #612 and formal authority.
