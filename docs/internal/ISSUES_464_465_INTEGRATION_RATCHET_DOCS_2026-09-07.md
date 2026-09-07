# Issues #464 / #465: joint final integration and documentation

This is the existing #465 Session 06 and #464 Session 07, continued from the
joint surface/persistence candidate. Request schema 7 and session version 41
remain current. The final local verdict, executed IDs, manifests and preservation
records belong to `/home/kawato/gbdraw-verification/joint-integration-docs-20260907/`.
This report creates no additional acceptance stage and authorizes no publication.

## Public documentation and executable evidence

The CLI, package-root Python, typed-request, Web, placement, input-table and
compatibility references remain their existing public owners. No public page was
added. The release notes distinguish the current 7/41 writer from historical
6/40 cardinality and authority semantics.

The combined example uses the existing annotated tobacco plastome presentation,
with source start 5500, multipart rps16 (`protein_id=NP_054479.1`) in outward lane
1, overlap resolution enabled and tolerance 1. Labels, functional colors, legend,
LSC/IRb/SSC/IRa annotations and GC content remain in the figure.

| Evidence | Existing producer / source |
|---|---|
| H-CLI-14 | `docs/recipes/run_cli_scenarios.py`; literal CLI reference block |
| H-PY-06 | `docs/recipes/run_python_scenarios.py`; literal Python reference block |
| Typed equivalent | literal `joint-typed` block in the typed-request reference |
| H-GUI-16 | `docs/capture/flows/how_to/joint_display_placement.py`; actual upload, controls, Generate and Save |
| Source provenance | `docs/capture/joint-source-verification.json`; official NCBI bytes equal the packaged mirrors |
| Current Gallery assets | `tools/refresh_gallery_sessions.py --no-assets`, then `--skip-session-refresh` |

Only two new operation crops were captured. Existing current sessions and their
registered dependent assets were regenerated with the existing transactional
producer; Gallery tutorials were not recaptured. The social preview and tracked
reference SVGs remain untouched. Historical positive fixtures retain original
bytes and provenance independently of current Gallery artifacts.

## Integration corrections and ownership

| Trigger and correction | Existing owner before and after | Verification |
|---|---|---|
| A custom split Circular slot lost the nominal separate-strand preset lanes. Preserve that preset for separate strands. | `diagrams/circular/radial_layout.py` | H-PY-03 against HEAD and its existing figure; placement and reference suites |
| Low-level Linear callers may supply features without a placement assignment. Check absence before inspecting requested placement. | `diagrams/linear/assemble.py` | entire SVG-ID contract owner; original and manual per-feature ribbon cases |
| Historical full config omits tolerance 0; replay writes it. Compare those as the same request while retaining nonzero differences. | `web/js/services/session-request.js` | real Gallery prepare/replay/finalize; zero/one negative equivalence test |
| A GUI-managed filter leaf was also disclosed as an unmanaged raw aggregate. Suppress only that derived duplicate; preserve explicit raw and unowned filter differences. | `web_support/config_overrides.py` | typed projection/negative tests and actual Load/Generate/Save/Reset browser journey |
| SVG comparison changed ElementTree's shared namespace registration and broke later inline browser fixtures. Remove the unused registration. | test helper `tests/utils/svg_compare.py` | comparison followed by all native display browser cases in the same process |

These are non-increasing corrections within existing owners and paths. There is
no new parser, placement planner, transform, authority era, runtime compatibility
path or privileged operator. Existing explicit requested values, validation and
supported continuation rules determine the outcomes; no unresolved product
choice was introduced. No `BD-###` is inferred from candidate code.

Rollback is the additional patch against the preserved Session 05 start archive;
it does not require restoring the original repository or old evidence.

## Report-only Ratchets

All C R1–R5 and D R1–R6 remain REPORT_ONLY. Counts, file/symbol inventories,
expanded formal IDs and actual artifact measurements are in the external Ratchet
and matrix reports. No blocking CI policy is enabled here.

The shared transform handles source/display conversion; pixel conversion of
already projected geometry is legitimate. The shared placement planner assigns
multipart features once; glyphs, labels, leaders, reservations, canvas bounds and
comparison endpoints consume the final assignment. Requested persistence contains
source-bound intent, not derived lanes or pixels. Dormant known features and
unknown/stale identities remain distinct.

The retained COMMON_POOL oracle has 216 cases: 192 original-identical and exactly
24 approved differences. GFF identifiers retain their complete source IDs.
Feature-associated ribbons attach to each endpoint feature's painted lane with
4 SVG-unit clearance. These approved differences do not permit another default
geometry change. H-CLI-11 also changes under the existing Phase 4 requirement
to honor Circular separate-strand resolver ON: HEAD silently disabled it. The
same candidate with resolver OFF reproduces the HEAD figure; the regenerated
ON figure retains the recipe and now resolves its overlaps. D's literal zero-drift promotion condition is not met by the
COMMON_POOL exception, even when its functional integration passes; this is not
a new functional blocker. Promotion eligibility is recorded separately from the
local completion verdict.

## Verification and handoff boundaries

The external reports bind the shared Python, JavaScript, architecture, functional,
PR smoke, compatibility, recipe, reference, offline and build checks to source,
fixtures and the prepared wheel. C/D counts are deduplicated. Failure logs remain
alongside the resolved results. Read/decode/Worker counts use the original
operation sequences and unchanged assertions.

Performance evidence includes three trials each for Circular multi-record,
Linear multi-record comparison and a 155,943-row depth track, source immutability,
identity-path allocation and memory observations. Browser evidence includes
1280×720 and 390×844 controls and actual saved/reloaded/regenerated figures.

The local result does not establish release readiness for all v0.14.0 work,
all-path resolution of #468/#469/Vibrio, remote CI, or publication. Normal commit/PR
preparation must use this candidate and its final manifest, preserve the original,
backup and old evidence, and run required CI against the future exact commit.
Stage, commit, push, PR, merge, tag, release and deploy are outside this task.
Cloudflare-only deployment policy remains in force.
