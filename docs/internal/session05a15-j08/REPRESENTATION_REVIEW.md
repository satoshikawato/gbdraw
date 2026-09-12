# J08 representation review

Status: implementation authorized after technical review, conditional on preserving
all JSON formats that reached main. The Product outcome
is already selected: Save before loading any source preserves a settings-only
Session. This note does not reopen that decision or create a Product authority.

Base: `3548fe14baf1fe3f9bfa83dcef5f34337523e8c8`, tree
`c96c0c273d1ce2432d07fb6ef46253d6e9151cde`. Fetch found no intervening changes.
Work branch: `work/session-05a15-j08`, no upstream. The dirty original workspace
and the historical programs are unchanged.

## Why a grammar change is necessary

`config.js::normalizeSessionData`, `session-authority.js::adoptCurrentSessionDocument`
and `session_io.py::validate_session` require an object-valued canonical request
for current Sessions. `session-request.js::buildCanonicalRenderRequest` requires
at least one real record and all its resources. The reproduced ordinary Save
fails with `Canonical resource record-1-genbank is missing.` A null request is
not an existing session-41 variant. Session 41 is present in first-parent main
at `4e8c9380`; its supported reader contract must remain intact.

## Proposed exact format change

Advance only the Session writer to **42**. Keep request **7**, catalog **3** and
Web bindings **2**, with existing schema-1 reads. Add no new top-level field.
For a settings-only document the required distinguishing fields are:

```json
{
  "format": "gbdraw-session",
  "version": 42,
  "renderRequest": null,
  "results": [],
  "editorState": { "featureCatalog": null }
}
```

This fragment is not a complete document. The accompanying
`proposed-settings-only.json` is a complete proposed example with the actual
non-default config from the red test. It also requires an object-valued
`resources`, explicit `webFiles.bindings`, valid full Web `config`, and
`ui.mode` (`circular` or `linear`) plus the saved input types. Ordinary full
Sessions also write version 42 and retain their existing object-valued request.
Omitted `renderRequest` remains malformed; explicit null carries the absence.

Admission of null additionally requires all of the following:

- No biological binding in either mode: Circular GenBank/GFF/FASTA, every
  Linear row's GenBank/GFF/FASTA, or conservation FASTA/sequence-source bindings.
  Historical direct biological source inventories must also be absent or empty.
- No saved Result, feature catalog, committed provenance (`cliInvocation`),
  generated comparison cache/identity payload, or resolved run geometry.
  Empty existing artifact envelopes remain empty; no synthetic artifact is made.
- Resources belong to retained auxiliary bindings; dangling bindings, invalid
  descriptors and unowned resources fail validation. Settings-only does not
  mean an empty resource table. Color, filter, qualifier-priority, depth and
  precomputed comparison files keep their existing bindings and exact bytes.
- Current config and both saved mode profiles pass their existing validation;
  malformed configuration cannot fall through as a legacy config import.

The runtime classifier checks the complete input inventory and the committed
owner before request construction. An inactive-mode biological input or retained
committed render excludes settings-only. A source-bearing inconsistency fails;
no missing-resource exception selects this variant.

## Preservation and ownership

Preserve the current `buildConfigData` / `CURRENT_WRITER_ACTIVE_CONFIG_DOMAINS`
inventory, mode profiles, layout preferences and supported editor preferences.
There is no separate preset subset or serializer. Manual matchers, palette and
legend rows remain under their current owners. Auxiliary file views and encoded
payloads use `session-resources.js` and `session-resource-backing.js`.

`session-authority.js` admits the document variant while returning **no canonical
owner** for settings-only. `adoptRuntimeCanonicalSession` retains its object
request/resource guarantee. `session-request.js` reuses its existing binding
projection without decoding a nonexistent render request. `config.js` applies
that projection in its existing preflight/capture/reset/apply/rollback import
transaction. Loading over existing work replaces that Session, clears old
sources/Result/committed request and imports the saved config; failure restores
the prior state. This follows the existing explicit Session replacement boundary
and OIPC-C06; it does not infer deletion from a mode switch or empty control.

No Generate or fabricated Worker operation is needed for Save/Load. A later
real source enters the normal canonical request/Generate path. Existing missing
input validation still handles Generate before a source is supplied.

## Reader compatibility and CLI edge

Current Web and Python readers accept 42 with these constraints. Versions
27–33 and 39–41 retain their existing support and requirements. Version 41 with
a null/missing request remains invalid. Released schema-1 bindings remain
supported; schema 2 is admitted in sessions 41 and 42. No migration chain or
second resource namespace is added. Existing Gallery session-41 artifacts remain
valid released-format inputs and need no refresh for this change.

Python `load_session_document` and materialization can validate and expose a
settings-only document and its auxiliary resources. Typed conversion and CLI
replay reject it with an explicit settings-only/no biological render request
error. They do not fall back to legacy argv, invent a sequence, or render an
empty figure. CLI/Web full-session replay keeps its existing semantics.

Readers limited to session 41 reject session 42 by version, including newly
written full Sessions. This is the material compatibility consequence of one
current writer version. Keeping 41 for full writes and 42 only for settings
would create two current writer contracts; the proposal avoids that split.

## Review requested

Approve or amend only the session-42/null-request grammar and the explicit
nonrenderable CLI edge above before dependent runtime changes. This review is
required by the user's attached instruction, section 4: “If a persisted
grammar/schema change is genuinely necessary, present its minimal exact shape
and compatibility consequences for required review before implementing that
change.” It is not a request to approve saving again.

The permanent regression and `red-proof.json` record the unchanged production
failure. Implementation and candidate acceptance are in progress; PR and merge remain pending. No red-test-only PR is proposed.

## Maintainer response

> mainに実装されている(=ユーザーが実際に使い、ファイルを保存した可能性がある).jsonファイルは確実に今後も使える限りにおいて、実装を進めてよい。

This response authorizes the proposed representation subject to existing main
compatibility. It does not authorize push, PR publication or merge, and is not
recorded as an implementation review or an accepted machine Gate result.
