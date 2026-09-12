# J08 implementation review

The selected Product outcome is: “Save Session before any source file is loaded
should preserve a settings-only session.” The maintainer approved the session-42
representation subject to preserving formats that reached main; the exact
response and compatibility consequence are in [the representation review](REPRESENTATION_REVIEW.md).
No Product choice, implementation approval, or publication permission is inferred
from a machine check.

## Owner and path evidence

| Responsibility | Before and after owner | Change |
|---|---|---|
| Session assembly and import transaction | `services/config.js` | Positively identify source-free state before request construction; apply it through preflight/capture/reset/apply/rollback. |
| Full editable configuration | `session-active-config-contract.js`, `buildConfigData` and mode profiles | Reuse the entire inventory. Preserve unresolved track placement and empty input rows on source-free Load. |
| Document admission | `session-authority.js` | Admit only explicit version-42 null requests with no biological inputs or committed artifacts; return no canonical owner. |
| Canonical render request | `session-request.js` | Its renderable request contract is unchanged. A settings projection reuses the existing binding projector. |
| Encoded bytes, bindings and file views | `session-resources.js`, `session-resource-backing.js` | Existing assembly accepts explicit absence of a request and retains auxiliary resources. Missing inputs never trigger a fallback. |
| Python document admission and materialization | `session_io.py`, `session.py` | Load settings-only documents; require a real canonical request for typed conversion/rendering. |
| CLI replay | `cli_utils/session.py` | Reject nonrenderable settings explicitly before mode replay; preserve full-session replay and draft bindings. |

No owner is moved, duplicated or superseded. No new module, reactive authority,
Worker, resource namespace, canonical Generate entry, or compatibility migrator
is introduced. Versions 40/41 continue through their existing full-session
admission path. The new 42/null branch expresses a current document variant,
not a reader for a superseded contract. Registered owner locations, canonical
entry edges, privileged permissions and import cycles remain unchanged.

The rollback for this implementation is a revert of the coherent change before
publication. Once users can save version 42, retiring that format would require
its own compatibility decision; an older reader cannot consume 42 merely by
changing its version number.

## Product preflight and completion boundaries

The source-free outcome is explicit user authority. The existing Session import
transaction and Session replacement semantics govern Load over current work.
The representation/CLI boundary received the required technical review before
implementation. No candidate `BD-###`, retired option, new timeout/SLO or waiver
is used. The implementation keeps the following contributions together:

- Save preserves all retained settings, both mode profiles and auxiliary bytes.
- Load restores those settings and clears the previous Session's source, Result
  and committed request. Invalid candidates leave the old complete state intact,
  including failure after the import reset has begun.
- Save→Load→Save preserves semantic content. No-source Generate reports the
  ordinary input error; adding real GenBank data uses the usual Generate path.
- Generated Result and History use the shared strengthened visual comparator at
  completion, before Save, Generate, export or flush can repair anything.
- Full Sessions with an inactive-mode source remain full Sessions. Mode changes
  do not become source replacement. Existing draft B / Result A / request A
  and mode/restore controls remain in the ordinary test tier.

## Released settings JSON versus a missing request

Current main is `4e8c93804186f9c4b163b584bd81d759b1b3522d`. Its Session writer is
41 and accepts sessions 27–33/39–41. Running the permanent source-free Save test
against **unchanged main Web files** fails at the actual Save boundary with
`Canonical resource record-1-genbank is missing.` It produces no Session download.
This independently confirms the dev red proof; it is not a reason to discard an
existing supported settings file.

Main also has `isLegacyConfigPayload` / `applyLegacyConfigPayload`: settings JSON
with no `format` and with recognized configuration keys uses this separate path.
That path remains unchanged. The existing `bare legacy configuration drives the
next canonical request and SVG` browser control passes on the candidate. A
version-41 Session with a missing/null required request is a different, malformed
envelope. The released schema-1 fixture retains its main provenance and exact
bytes. Schema-2 and historical full-session controls remain positive.

## Diff review

Production changes are confined to the existing Session owners and writer-version
constants in Gallery admission/migration. Render, History, mode guards, canonical
request schema and resource backing remain under their existing owners.

Test changes add the portable J08 browser/JS/Python controls and advance assertions
for **new writes** to 42. Published Gallery fixtures remain exactly version 41;
their tests assert 41 explicitly while retaining every request, catalog, resource
and output assertion. New staged Gallery publication still requires the current
writer. No Gallery artifact is rewritten.

The retained adapter changes only its writer/settings oracle. The original
journey programs, hashes, H01 failure record and existing strengthened comparator
remain unchanged. Permanent tests require no historical archive. The PR smoke
selection remains 10, and CI tiers, assertions, budgets and timeouts are unchanged.

Documentation records the representation and compatibility boundary. Tutorial
capture scripts change only their expected current-writer versions. Four Session
downloads from H-PY-05, T-PY-08 and H-CLI-12 are regenerated with their existing
recipes. Only the Session version changes semantically; CLI artifacts also refresh
creation and input modification timestamps. SVGs, resources' bytes and canonical
requests are unchanged. Generated browser and native wheels are verification
artifacts and are not committed.
