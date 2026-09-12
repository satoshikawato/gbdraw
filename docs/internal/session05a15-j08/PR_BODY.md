## Plain-language summary

Save Session now preserves settings before any source file is loaded. Loading that file restores the editable settings and both mode profiles, and replaces any previous sources and Result. Users can then load real sequence data, Generate, and use Undo and Redo with the restored settings.

## Format and compatibility

New Sessions use version 42. Settings-only files explicitly contain `renderRequest: null`, empty Results and a null feature catalog; auxiliary files retain their existing resources and bindings. Full Sessions keep a real canonical request. Request schema 7, catalog schema 3 and bindings schema 2 are unchanged.

Readers retain sessions 27–33/39–41, schema-1/2 bindings and legacy settings JSON. CLI and typed rendering give an explicit error for a settings-only document until a biological source is supplied in Web. Readers limited to version 41 reject newly written version-42 files, including full Sessions. The maintainer approved this representation subject to preserving formats supported on main.

## Implementation and verification

The change reuses the existing Session, configuration and resource owners, including atomic import and rollback. It does not create a source, empty sequence, Result or invalid canonical request. Published Gallery Sessions remain unchanged in version 41. Four current-writer documentation downloads are regenerated without changing request, resource bytes or figure semantics.

The permanent browser test failed on unchanged dev at the actual Save boundary with `Canonical resource record-1-genbank is missing.`, then passed with identical assertions. Extended coverage checks full settings, both mode profiles, auxiliary bytes, replacement/rollback, inactive-mode sources, Generate, Undo/Redo and repeat generation. Original J08 and J16 pass through the retained adapter; archived programs and the H01 failure record are unchanged.

See [acceptance evidence](https://github.com/satoshikawato/gbdraw/blob/work/session-05a15-j08/docs/internal/session05a15-j08/ACCEPTANCE.md) for measured test counts, scoped reruns, source/installed-wheel provenance and reproduction commands. PR smoke remains 10 tests; CI selection, behavioral assertions, budgets and timeouts are unchanged.

## Review status

Web policy reports Gate PASS and Review REQUIRED under the architecture profile. The representation approval is recorded in [the technical review](https://github.com/satoshikawato/gbdraw/blob/work/session-05a15-j08/docs/internal/session05a15-j08/REPRESENTATION_REVIEW.md); it is not an implementation approval. Push, PR publication, merge and acceptance on the resulting dev commit have not occurred. No main promotion, tag, publication, S11/S12, full 48-journey run or deferred Q01 dense stress is included.
