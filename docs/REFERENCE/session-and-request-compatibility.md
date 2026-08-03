[Documentation home](../DOCS.md) | [Typed requests](typed-requests.md) | [Output reference](output-formats-and-export.md)

# Session and request compatibility

Current writers emit session version 40 and canonical `renderRequest` schema 5.

| Persisted format | Current writer | Accepted by current readers |
|---|---:|---|
| gbdraw session | 40 | 27–33 and 39–40 |
| canonical `renderRequest` | 5 | 1, 2, and 5 |

Session versions 34–38 and request schemas 3–4 were development-only and are rejected. The typed-session bridge converts sessions 31–33 and 39–40. Sessions 27–30 remain CLI replay inputs because they do not contain a canonical request.

Both `.gbdraw-session.json` and lossless `.gbdraw-session.json.gz` are accepted. The web app writes compressed sessions by default. The CLI writes JSON unless the explicit session output name ends in `.gz`.

`render_session()` is the persisted compatibility boundary. It migrates supported artifacts into current request artifacts and calls the current typed renderer. `render_request()` does not parse old session envelopes.

Session 40 stores file bytes once under `resources`, binds web inputs through `webFiles`, keeps the last committed render in `renderRequest`, and stores editable Linear comparison intent separately under `config.linearComparisonPlan`. Current result storage has one base SVG per logical diagram; supported older paired base/interactive results collapse on read.

For the version-by-version record of retired input names, saved protein-result schemas, deterministic handles, and materialized-resource lifetime, see the [compatibility history](../SESSION_COMPATIBILITY.md). Release notes describe when changes occurred; this page is the current support authority.
