[Home](./DOCS.md) | [Current compatibility reference](./REFERENCE/session-and-request-compatibility.md) | [CLI inventory](./CLI_Reference.md) | [Python API](./REFERENCE/python-api.md) | [Typed API](./REFERENCE/typed-requests.md) | **Compatibility history**

# Session and request compatibility history

This page retains the detailed version-by-version migration history for gbdraw
session files, canonical render requests, and saved LOSAT results. The concise
[session and request compatibility reference](./REFERENCE/session-and-request-compatibility.md)
documents current support. Tutorials and the FAQ describe what a user should
do; release notes record when a format changed.

## Unreleased: Session 46 keeps the settings and the Result of each diagram mode

In the Web app, Circular and Linear each keep their own settings, edits, and
Result. A setting or an edit made in one mode, and a Generate, stay in that
mode; **Reset Settings** resets both modes, and one Undo restores both.
Switching modes changes no setting and shows that mode's Result with its moves
and edits, or the empty Preview.

**Drafts.** Session version 46 keeps each mode's draft in `modes.circular` and
`modes.linear`. A slice holds the mode's `config` (the former `config` without
`modeProfiles`, `cliOptions`, `paletteInstantPreviewEnabled`,
`adv.rich_feature_popup`, and the LOSAT execution settings), its `features`
edits, its Legend and stroke edits (`editorState.legend` and
`editorState.featureStrokes`), and its `ui` values: the mode's layout slot
(`layoutPreferences`), Canvas padding, the pending palette, the Linear
typography link, and the Features-list record. A slice holds only these
fields, and may omit any of them: Load reads an omitted value as that mode's
default, after the committed request of that mode's Result. App-level settings
moved to the top level: `ui.losatExecution` (LOSAT run mode, thread budget,
threads per job, parallel workers), `ui.richFeaturePopup`, and `cliOptions`
(the options of the CLI command that wrote the Session, unchanged). Per-feature
edit rows, Feature placement rows, record display rows, and `featureIdentity`
annotation targets no longer name a mode (`scope`): the slice is the mode, and a
row is keyed by `recordKey` and `biologicalFeatureId`. A slice's
`editorState.legend.entries` lists the Legend rows of the mode's Result, then
the renamed rows its last Generate did not draw (for example with GC content off
or Show Depth off), each marked `"dormant": true`; the rename and the styles
under the new name apply again when a Generate draws the row. A Session 46 that holds
the top-level `config` or `features`, a former `editorState` or `ui` draft
field, or a field outside these slices is rejected
(`INPUT_INVALID {field: schema, reason: FIELDS}`). A settings-only Session keeps
the slice of its shown mode. The CLI and the Python API write no `modes`; a CLI
re-save keeps the Session's `modes` and `cliOptions`.

**Results.** The optional top-level `otherModeResult` holds the other mode's
Result set. In a Session the Web app writes, the top-level `renderRequest`,
`results`, `editorState.featureCatalog`, `runMetadata`, and `cliInvocation` hold
the shown mode's set when it has a Result, otherwise the other mode's set;
`otherModeResult` holds the remaining set only when it has a Result. Its fields
mirror the top-level ones: `renderRequest` (the other mode, the same request
schema), `results` (at least one), `editorState` (`featureCatalog`,
`alignmentResetReceipt`, and the Result's generated Legend order, colors, and
stroke defaults), `ui` (the selected Result and the generated legend, title,
and palette), `runMetadata`, and `cliInvocation`. Its request names resources
in the one top-level `resources` table. A settings-only Session
(`renderRequest: null`) cannot hold the field. **Load Session** shows the saved
mode (`ui.mode`) when that mode has a Result, otherwise the mode that has one;
the other mode's Result waits for the mode button. A Gallery Session shows one
mode, holds that mode's slice only, and has no `otherModeResult`.

**Older Sessions.** Loading a Session 27–44 splits its one draft into the two
slices: a shared setting goes to both modes, a setting of one mode to that
mode, and the Legend edits and per-feature edits to the mode of the saved
Result. The draft of the mode that was not shown when the Session was saved
(`config.modeProfiles`) gives that mode its title, font sizes, and comparison
thresholds; a Depth or Show Depth value goes to a mode only as far as that
mode's Depth files reach, and an annotation bound to one mode's record goes to
that mode. Session 45 was a development-only format and is rejected; load its
source files instead.

On the command line, `gbdraw circular --session` and `gbdraw linear --session`
render the set of their own mode, at the top level or in `otherModeResult`. A
re-save with `--save_session` or `--session_output` writes Session 46: the
subcommand's set at the top level, the other set in `otherModeResult`, and an
older Session's draft split into `modes` as the Web app splits it. It keeps
`ui.mode`, so the Web app still opens on the saved mode. The re-save replaces
the subcommand's set with the new render and keeps the other set. It also
changes these parts that both sets share:

- `losatCache` holds the entries that the render returns, and
  `losatDerivedCache` is emptied.
- `proteinIdentityManifest` and `legacyArtifacts` are replaced by the render's.
- When the render migrates legacy protein IDs, the protein references are
  rewritten throughout the Session, `otherModeResult` included.
- Resources that neither set's request nor the Web files name are dropped.

## Unreleased: Web Load of 0.13.0 Sessions

Session version 46 is unchanged. Sessions 27–30, written by gbdraw 0.13.0 and
earlier, have no canonical `renderRequest`. The Web app loads such a Session
from its saved settings, with its saved preview; the CLI replays it, and the
typed-session bridge does not convert it. The Web writers of Sessions 27–33
saved every Circular Custom Track Slots row with `spacing: null`, a field the
current slots do not read. When the Session's Custom Track Slots are off, that
null carries no setting and Load drops it, so the 0.13.0 Gallery Sessions load.
When the slots of a Session 27–30 are on, or a row holds a `spacing` value or
another retired field (`strict`, `compress`, `reserve`, `placement`,
`inner_radius`, `outer_radius`, or `gap_after`), Load fails with a message that
names the field and the track row, and the previous Session stays loaded. The
current Custom Track Slots use `radius`, `width`, `inner_gap_px`,
`outer_gap_px`, `side`, and `z`.

## Unreleased: Web Load of Session 31–39 table rows

Session version 46 is unchanged. The Web writers of Sessions 31–39 stored the
Default colors, Label whitelist, and Qualifier priority tables without
normalizing their cell values, so a value typed with a tab became extra cells
and a value with a line break became a short row. Since the Web table readers
reject a row with the wrong number of cells, Load rejected such a Session. Load
now reads these three tables of a Session 31–39 as the current writer writes
them: the extra cells of a row join its last column with one space, and a row
without its required columns is dropped (`feature_type` and `qualifier` for
Label whitelist, `feature_type` and `priorities` for Qualifier priority, and
`feature_type` and a color for Default colors). The Load notice names each
table and line read this way. Sessions 40 and later are read as before, and
table file imports and CLI replay (`--session`) still reject such a row. When
another table of a Session fails to load, the message names the table, for
example `Session table: Specific colors.`

## Unreleased: Session 46 and Web feature edits by source identity

Session version 46 stores the Web app's Feature visibility, Label visibility,
and label text edits in each mode's `features.featureOverrides`: one row per
original-source feature, keyed and named by `recordKey` and
`biologicalFeatureId`, with the same fields as request
`diagramOptions.featureOverrides` and the label's original text
(`labelSourceText`). Feature placement drafts (`config.featurePlacementOverrides`
of the slice) are keyed the same way. Both modes can use the same record key;
each mode's slice holds its own rows, so a request carries only the rows of its
own mode. The four rendered-ID maps (`featureVisibilityOverrides`,
`labelTextFeatureOverrides`, `labelTextFeatureOverrideSources`,
`labelVisibilityOverrides`) are rejected in Session 46. The Web app sends these
rows in the request, loads a request whose `featureOverrides` array is not
empty, and draws them in the live preview by identity. An annotation made from
selected features is saved with a `featureIdentity` target in the slice of the
mode it was made in. Loading a Session 40–44 moves an annotation
target with one `hash=` feature selector, as a selection made it, to a
`featureIdentity` target in the mode of the Session's diagram only when the
figure cannot change: the record the target binds is drawn without a crop,
reverse complement, or rotation, and the hash names exactly one feature of the
Session's saved feature catalog, in that record. Load reports how many targets
moved. Every other annotation target of an older Session, including one of a
Session without a saved catalog, loads unchanged and keeps its meaning: the
hash of the drawn feature. The CLI moves the same targets in
`config.annotationSets` when it replays a Session 40–44 with `--session_output`
or `--save_session`, and logs how many moved. As after a Web Load, the
written request keeps the `hash=` targets that drew the figure.

The feature catalog is schema 5. Each drawn feature records the hash,
location, and record location it was drawn with (`drawnSelector`), which live
rule matching uses. A schema 3 or 4 catalog reads as schema 5 with no selector
values. Until the next Generate, a feature whose rendered ID carries its source
hash was drawn with its source coordinates, so live matching and **Load Label
TSV** use its source values; on other features (cropped, reverse-complemented,
or rotated records) live matching leaves `location` and `record_location` rules
to Generate, and **Load Label TSV** declines a table with such rows and says
why.

Loading a Session 44, or a Web Session 31–33, moves each rendered-ID edit onto
its feature. In a Session 40–44, a rendered ID in the saved feature catalog
names its feature; otherwise its `_record_<n>` and `__instance_` suffixes are
removed and the edit moves only when exactly one feature remains. A Session
before 40 has no catalog that Load reads, so Load reads its GenBank sources
again with its crops and orientations and matches each rendered ID's drawn
hash and record position. Without readable sources it uses the saved feature
metadata, which names only features of records drawn without a crop or
reverse complement. Other edits are dropped and Load reports how many. If the
diagram runtime cannot start or fails while reading the sources, Load fails
and keeps the current Session, so no edit is dropped. An older
Feature visibility edit hid every feature with the same hash, such as each copy
of a duplicated record; it now applies only to the feature that was edited, and
Load reports how many edits the next Generate draws differently for this.
Each moved edit goes to the slice of the Session's diagram's mode. The CLI moves the
edits of a Session 40–44 through its saved catalog in the same way when it
replays the Session with `--session_output` or `--save_session`, and logs these
counts. It does not read the sources of an older Session again, so it drops
that Session's edits; the replayed request's tables keep their effect on the
diagram. A Feature placement
draft of a Session 41–44 reached every request with its record key: a lane
placement goes to the slice of its side's mode, and a Main placement to both
slices. The CLI applies the same mapping when it replays such a Session with
`--session_output` or `--save_session`.

## Unreleased: request schema 9 and feature identity overrides

Session version 44 is unchanged. Canonical request schema 9 adds the required
`diagramOptions.featureOverrides` array and the `featureIdentity` annotation
target. Each override row names one original-source feature by `recordKey` and
`biologicalFeatureId` and sets its Feature visibility, Label visibility, or
label text (see
[Feature identity overrides](./REFERENCE/typed-requests.md#feature-identity-overrides)).
Schema-8 requests are read as requests without overrides, and current saves
write schema 9. The Web app writes an empty array and cannot yet load a request
whose array is not empty or that has a `featureIdentity` target; it reports
this instead of dropping the edits. An exact Feature placement whose identity
the source does not have no longer fails the render; it becomes a feature
identity notice.

## Unreleased: value checks, derived label maps, and definition spacing

Session version 44 and request schema 8 are unchanged. Replaying a Session, or
generating from it, applies the value checks shared with the CLI and Python API:
window, step, and depth window or step must be positive integers, the
dinucleotide must be two letters from `A`, `C`, `G`, `T`, and `U` (`U` counts
as `T`), font sizes must be greater than zero, and stroke widths must be zero
or greater. A Session holding another value fails with the field or setting
named instead of drawing an empty or flat track. The Web app sends these values
as typed and no longer replaces a rejected value with Auto or a default.

A Web Session older than Session 40 loads with its saved preview, but its Result
has no current feature metadata. **Save Session** then asks for one **Generate
Diagram** and offers it; after that Generate, Save writes Session 44. Loading
such a Session and saving it directly, as 0.13.0 allowed, is retired. Session
40 and later Sessions save unchanged.

CLI Sessions and Sessions saved from `main` can contain the label maps
`whitelist_map`, `priority_map`, and `label_override_rules` under
`labels.filtering`. They are compiled from the label tables during a render and
are not settings. Load no longer keeps them as a preserved
`labels.filtering.raw` setting, and attaching a label table recompiles them, so
Qualifier Priority and whitelist edits take effect.

Web Sessions no longer write `features.labelOverrideContextKey`. Readers ignore
it in Session 44 files that contain it; label edits no longer depend on which
Result was displayed when the Session was saved.

A Circular request that sets `objects.definition.circular.font_size` without
`objects.definition.circular.interval` uses the font size plus 2, truncated to
an integer, as the definition line interval. This restores the 0.13.0 Web and
CLI spacing. The rule is applied when the overrides are applied, so the stored
request is unchanged and replay writes the same overrides back.

## Unreleased: default label of a precomputed similarity ring

Session version 46 and request schema 9 are unchanged. A Circular ring from a
BLAST table (`--conservation_blast`, Python `ComparisonRingTrackOptions(source=...)`)
without a label is now labelled with the table's file name without the last
extension, as in the web app. CLI and Python Sessions now store
`conservationLabels` for every precomputed ring, as Web Sessions already did.
A Session without them was written before this change and drew the full file
name, which its replay still draws. A Session 27–30 replay now takes that name
from the Session's file bindings instead of drawing the temporary copy's name
(`arg3-<file>`). Web Load of such a Session also shows the full file name; it
showed the name without the extension before.

## Unreleased: a scale interval of 0 or less

Session version 46 and request schema 9 are unchanged. `--scale_interval` and
`objects.scale.interval` in the Python option objects (`config` and
`config_overrides`) reject a scale interval of 0 or less, and the web app's
Scale Interval field starts at 1. A render request and a Session read a value
of 0 or less as the automatic interval, as every earlier writer drew it:
in `renderRequest.diagramOptions.config`, in `configOverrides` (including the
flat `scale_interval` that earlier Web Sessions wrote), and in a Session 27–30
`cliInvocation`. Web Load shows the stored value in the field, **Generate
Diagram** draws the automatic interval, and the Source recipe omits
`--scale_interval`. A value of 0 or less typed into the field behaves the same.

## Unreleased: CLI and Python LOSATN / TLOSATX results

Session version 44 and request schema 8 are unchanged. A Linear run with
`--losat losatn` or `--losat tlosatx` (Python `losat_search` with those
programs) saves what the web app saves: one `nucleotideBlast` resource per
compared record pair with the raw search-frame rows, and one schema 2
`losatCache` entry per pair with the web raw key and the non-key `runtime`
record. Entries that the web app searches carry `runtime` too, as
`{kind: "losat", source: "wasm", version: null, program}`. A CLI or Python
runtime outside the package records only its executable's name as `path` (for
example `losat` or `blastn`), so the runtime record carries no local directory;
the bundled runtime records `gbdraw/bin/<platform>/losat`. A CLI Session's
`cliInvocation.args` still hold `--losat_bin` and the input and output paths
as typed. The saved request carries the resolved comparisons, not the search
intent, so replay needs no LOSAT runtime. A request that still carries the
search intent cannot be encoded; resolve or render it first.

## Unreleased: one LOSAT cache entry per raw key

Session version 44 is unchanged. The Web app saved a raw key once per row, so a
Session with two Circular rings of one sequence repeated the key and could not
be loaded. It now saves one `losatCache` entry per raw key, named by the first
row, as the CLI does. After a load, the cache list shows one row per raw key
until the next Generate. A version 39 or later Session that repeats a raw key
is rejected with an `INPUT_INVALID` Session diagnostic.

## Unreleased: comparison rows in the search frame

Session version 44 and request schema 8 are unchanged. Comparison rows stored
in a Session (`nucleotideBlast` resources, uploaded tables, generated LOSATN
text, and `linear_comparisons`) use the search frame of each record as the
Session persists it: the selected and cropped record, 1-based, on its source
strand. A record saved with `presentation.reverseComplement` or a reversed
region keeps its rows in the search frame, and the planner projects the
orientation when it draws. A stored row outside 1..L of its record stops the
replay with a `COMPARISON_INPUT` error whose reason is `SEARCH_FRAME`.

Sessions saved from `main` (version 42 and older) stored the rows of a
reverse-complemented Linear record after the reverse complement. Both readers
convert them once at Load with `x -> L + 1 - x`, where L is the length of the
selected and cropped record:

- CLI replay (`--session`) converts the rows of the adapted request.
- The Web app rewrites the stored table bytes, both uploaded tables and
  generated LOSATN text, so a Session saved after Load is current. A table
  downloaded from that Session therefore differs from the file first uploaded
  to `main`. Saved LOSAT raw cache entries were already in the search frame
  and are unchanged.

A CLI sidecar of `-b` with `--reverse_complement` or a reversed region embeds
the reversed record as a sequence without `reverseComplement`, so the writer
stores the `-b` table rewritten into that sequence's coordinates; replay draws
the ribbons of the original run. CLI sidecars written by `main` already stored
the rows that way and replay unchanged.

Linear Sessions saved by the `main` Web app write
`orthogroupState.selectedOrthogroupAlignmentFeature: ""` when no alignment
target is selected. Readers now treat the empty string as no target instead of
rejecting the Session.

## Session 44: typed Similarity alignment display state

The current Session 44 writer uses canonical request schema 8 to move
Linear Similarity Group alignment out of generated-protein pipeline settings. The Linear request layout now owns
finite `x` and `y` base translations keyed by `recordKey`, plus an optional
resolved `SimilarityAlignmentPlan`. The plan records the exact reference feature,
one validated decision per displayed record, and its rationale. Record
presentation or region state owns orientation, and anchor centers are projected
from the resolved record display. The nested plan is schema 2. Its unreleased
schema-1 predecessor has no current reader; the Session remains version 44 and
the canonical request remains schema 8.

Current writers never emit `align_orthogroup_feature`,
`alignOrthogroupFeature`, or the former Session-only
`selectedOrthogroupAlignmentFeature` copy. A Python
`SimilarityAlignmentReference` is resolved before saving; the Session stores
only its resolved plan. Supported request schemas 1, 2, 5,
6, and 7 retain a version-bounded reader for the old string. The reader keeps
that value private
until the historical selection can be materialized from saved stable feature and
orthogroup metadata; a successful save writes only schema-8 plan and translation
state. Loading a saved preview remains Worker-lazy.

## Session 44: feature anchors and independent record transforms

Session 44 stores source feature anchor capability metadata in feature catalog
schema 4 and preserves per-record display start, absolute orientation, and
feature-anchor provenance in one editable record-display draft. Released
Session 42/catalog 3 documents remain readable. Exact single-part locations
are migrated conservatively; ambiguous compound locations require Generate
again before feature-based rotation. Development-only Session 43 is rejected.

Session 44 also stores the selected visibility mode for Linear **Accession** and
**Length / Coordinates** independently. Each field is `auto`, `show`, or `hide`.
Auto shows the field while every effective rendered row contains one record and
hides that field diagram-wide when any row contains two or more records.
Disabled Record Layout ignores dormant shared-row values, and changing row
placement does not rewrite the selected mode.

Session 44 mode profiles also hold the per-mode Plot Title, Plot Title font
size, and Definition font size. A flat value without a per-mode entry, as in
Session 42 and earlier Session 44 files, belongs to the active mode only; the
inactive mode starts from fresh defaults. In every accepted version a missing
`config.linearRecordLayout` means Arrange in rows is on, and a missing
`ui.linearTypographyLinked` means linked while both font sizes are equal.

Canonical request schema 8 retains the effective booleans introduced in schema 7.
Version-42 editable `true` values migrate to Show, `false` values migrate to
Hide, and missing booleans migrate to historical Show. A selected-mode field,
when present, takes precedence over the retired boolean. Current writers emit
only the selected fields. Saved Results and committed render requests are not
regenerated during Load; the migrated selection takes effect on the next
Generate.

## Session 42: settings before the first source

Web **Save Session** can preserve settings before any biological input is loaded.
Session 42 adds the explicit `renderRequest: null` document variant with empty
`results` and a null `editorState.featureCatalog`. The existing config, mode
profiles, editor preferences and file-binding owners retain the settings and
auxiliary resources. No empty sequence or synthetic request is persisted.

Both active and inactive biological inputs exclude this variant, as do committed
render artifacts. Missing requests and missing or dangling resources remain
errors. Loading settings-only replaces the old Session, including its source,
Result and committed request, through the existing atomic import transaction.
After loading a real source, Generate and History use the usual render path.

Request schema 7, catalog schema 3 and bindings schema 2 are unchanged. Bindings
schema 1 and supported sessions 27–33/39–41 retain their existing reader support;
session 41 still requires a canonical request. All new full Sessions also use
version 42. Readers whose maximum version is 41 reject these new files.

Python document loading and materialization accept settings-only documents.
CLI replay, typed request conversion and rendering require a biological request
and reject this variant explicitly; they never fall back to legacy argv.

## Supported versions

Current writers emit one session and request format:

| Format | Current writer | Accepted by current readers |
|---|---:|---|
| gbdraw session | 46 | 27–33, 39–42, 44, and 46 |
| Canonical `renderRequest` | 9 | 1, 2, 5, 6, 7, 8, and 9 |
| Web file bindings | 2 | 1; 2 in sessions 41–42, 44, and 46 |

Previously written Session 44 documents with request schema 7 or 8 and feature
catalog schema 4 remain readable. Current saves write request schema 9;
the saved preview is retained until the next Generate.

Session versions 34–38 and 45 and canonical request schemas 3–4 were
development-only formats. They were never released on the supported history and
are rejected.

The public typed-session bridge can convert full session versions 31–33, 39–42, and 44 to
a typed request. Versions 27–30 do not contain a canonical `renderRequest`, so
the bridge does not convert them: the CLI replays them with the same
`circular` or `linear` subcommand that created the session, and the Web app
loads them from their saved settings.

`render_session()` is the compatibility boundary for canonical session replay.
It migrates supported persisted artifacts into `CurrentRequestArtifacts`, then
calls the same current typed renderer used by fresh requests. Fresh
`render_request()` calls never parse old session envelopes or rewrite retired
protein identifiers.

Both `.gbdraw-session.json` and lossless gzip-compressed
`.gbdraw-session.json.gz` files are accepted. The Web app writes the compressed
form by default; the CLI writes uncompressed JSON unless the output name ends in
`.gz`.

## Current request ownership

Canonical request schemas 5 and 6 record Circular grouping explicitly as
`single`, `grid`, or `batch`. Schema 6 also records each input's cardinality;
schemas 1, 2, and 5 decode records as `exactly_one`. A single diagram or grid
has one output object. A Circular batch has one resolved output object per
record. `renderRequest.output.prefix` is the output-prefix owner.

The Web projects a selectorless Linear schema-5 card to explicit `all` when it
is saved with schemas 6–8. Legacy multi-record Web inputs already have explicit
selectors, so this preserves the embedded source records shown by the card.

Current sessions keep mode-specific layout values under
`ui.layoutPreferences`. Supported older parallel title, legend, and grouping
fields are migrated when read and are not written again.

The Web writer stores file bytes once under `resources`. `webFiles` binds those
resources to active and inactive input controls, so files shared by the
committed request and an editable draft are not copied into a second payload.
For Linear comparisons, `webFiles.bindings.linearComparisons` contains only a
stable comparison-edge ID and its file binding. Endpoint, source, inclusion,
and filename metadata are not duplicated there. Version 39 sessions with
legacy embedded `files` remain readable.

`renderRequest` owns the last committed render. Web config retains inactive
Custom Track stacks, disabled rows, draft Axis positions, and per-mode
comparison profiles. `editorState.featureCatalog` holds the schema-3 feature
catalog used by the saved preview and editor. Current sessions store one base
SVG Result per logical diagram; readers collapse paired base and interactive
Results from supported older sessions.

Linear comparison intent is an independent editable draft under
`config.linearComparisonPlan`. Its mode is `none`, `adjacent`, or `selected`;
it also stores the default adjacent source and stable edge metadata. Placement
alone remains under `config.linearRecordLayout`. Current writers do not store a
global `blastSource`, nested layout comparisons, per-record BLAST files, or
per-record LOSAT filenames. The editable plan can therefore opt out without
changing the last committed comparison artifacts in `renderRequest`.

Supported pre-40 Web sessions migrate directly to this plan. A disabled or
absent legacy layout retains the old adjacent LOSAT/upload behavior, an enabled
explicit list becomes `selected`, and an authoritative empty explicit list
becomes `none`. Legacy per-record uploads and custom filenames are attached to
their original positional gap by stable record UID. CLI-only replay sessions
do not gain a synthetic Web comparison draft: they load with **No comparison**.
A CLI-only session written with `--losat` keeps the adjacent
LOSATP comparison its CLI drew. The accepted session versions remain 27–33,
39–42, and 44.

## Retired inputs

Fresh CLI and Python requests reject these retired names or values. Supported
older sessions and canonical request schemas 1–2 migrate them before replay.
Retired CLI flags exit with status 2 and name their replacement. Retired Python
fields raise `TypeError`; no alias is accepted. Persisted Session names do not
change.

| Retired input | Current input |
|---|---|
| Circular `--multi_record_size_mode sqrt` | `--multi_record_size_mode auto` |
| Linear `--label_placement on_feature` | `--label_placement above_feature` |
| Linear `--track_layout spreadout` / `tuckin` | `--track_layout above` / `below` |
| `--depth_tick_interval` | `--depth_large_tick_interval` |
| `--feature_table` | `--feature_visibility_table` |
| `--collinear_max_gene_gap` | `--collinear_max_unit_gap` |
| Circular slot `spacing` | `inner_gap_px` and `outer_gap_px` |
| Circular slot `strict`, `compress`, or `reserve` | No direct replacement; geometry and reservation are derived from `side` |
| Linear `--protein_blastp_mode pairwise` / `orthogroup` / `collinear` | `--losat losatp --losatp_mode pairwise` / `similarity_groups` / `collinear` |
| Linear `--protein_blastp_mode none` | Omit it; without `--losat` no protein comparison runs |
| Linear `--losatp_bin` (`--losatp-bin`) | `--losat_bin` |
| Linear `--ncbi_blastp_bin` (`--ncbi-blastp-bin`) | `--ncbi_blast_bin` |
| Linear `--losatp_threads` (`--losatp-threads`) | `--losat_threads` |
| Linear `--protein_blastp_max_hits` | `--losatp_max_hits` |
| Linear `--protein_blastp_candidate_limit` | `--losatp_max_target_seqs` |
| Linear `--align_orthogroup_feature` | `--similarity_alignment_feature` |
| Linear `--protein_blastp_output FILE` | `--losat_output_dir DIR`, which writes `DIR/losatp.raw.tsv` |
| `LinearComparisonOptions(protein_mode=...)` | `losat="losatp"` with `losatp_mode="similarity_groups"` / `"collinear"` / `"pairwise"`; `"none"` becomes `losat=None` |
| `LinearComparisonOptions(blastp_executable=...)` | `ncbi_blast_executable` |
| `LinearComparisonOptions(candidate_limit=...)` | `max_target_seqs` |
| `LinearComparisonOptions(orthogroup_member_max_hits=...)` | `member_max_hits` |
| `LinearComparisonOptions(losat_executable="losat")` | `losat_executable=None` (the new default) |
| `LinearDiagramOptions(protein_blastp_mode=...)` | `losat_search=LosatSearchOptions(program="losatp", losatp_mode=...)`; `"orthogroup"` becomes `"similarity_groups"` |
| `LinearDiagramOptions(protein_comparison_pairs=...)` | `LosatSearchOptions(pairs=...)` |
| `LinearDiagramOptions(losatp_bin=...)` / `ncbi_blastp_bin` / `losatp_threads` | `LosatSearchOptions(runtime=LosatRuntimeOptions(losat_executable=..., ncbi_blast_executable=..., threads=...))` |
| `LinearDiagramOptions(protein_blastp_max_hits=...)` | `LosatSearchOptions(losatp_max_hits=...)` |
| `LinearDiagramOptions(protein_blastp_candidate_limit=...)` | `LosatSearchOptions(losatp_max_target_seqs=...)` |
| `LinearDiagramOptions(orthogroup_member_max_hits=...)` | `LosatSearchOptions(losatp_member_max_hits=...)` |
| Circular `--conservation_fasta` | `--conservation_sequence` (FASTA, GenBank, or DDBJ); recorded invocations are rewritten |
| Circular `--conservation_table` column `comparison_fasta` | `comparison_sequence`; Sessions store the resolved table, so they need no rewrite |
| `CircularDiagramOptions(conservation_fasta_files=...)` | `conservation_sequence_files`; the Session field stays `conservationFastaFiles` |

Current multiword long options use underscore spelling except for the documented
active aliases. `--annotation-table` remains an alias for
`--annotation_table`, and `--gc_content_tick_interval` remains an alias for
`--gc_content_large_tick_interval`.

The private `__gbdraw_legacy_spacing` key is read only from canonical request
schemas 1–2 and is never written by schemas 5–6. Pixel spacing migrates to
`inner_gap_px` and `outer_gap_px`. Factor-based spacing can be replayed but
cannot be saved losslessly in the current format; replace it with explicit pixel
gaps before saving a migrated session.

## Saved protein-comparison results

After a source replacement, the Web writer saves only protein raw cache entries
that resolve through the saved identity manifest. Valid entries for inactive
search settings remain eligible; saving does not clear the live cache or History.

A protein-search cache hit requires the same amino-acid sequences, selected
protein set, record-instance and feature bindings, query/subject direction,
program, and meaningful search arguments. Upload filenames, modification times,
resource names, and display-only aliases are not part of that identity.

Current sessions use these independent payload schemas:

| Saved payload | Schema |
|---|---:|
| Protein raw search result | 4 |
| Derived protein comparison | 3 |
| Protein identity manifest | 2 |
| Nucleotide raw search result | 2 |
| Canonical typed analysis resource | 3 |

Typed resource readers 1 and 2 preserve their saved ortholog path tuples as
explicit collections in the current model. Schema 3 writers store lossless DAGs
for newly inferred paths and explicit collections for supplied legacy corpora.
Counts are exact decimal strings. Request 8 and derived envelope 3 remain
unchanged; derived identity includes `pathRepresentation` to prevent
reusing an old analysis payload as a current helper result. Older releases that
only support typed schemas 1 and 2 cannot read the new typed resources.

Generated protein FASTA, raw QUERY/SUBJECT fields, protein maps, and derived
references use deterministic session-internal handles of the form
`h_[a-z2-7]{26}`. Those handles bind a record instance to a complete CDS
identity. They are not user-facing protein names.

Sessions 27–33 may contain schema-2 protein candidates and schema-1 derived
evidence. Current readers keep those candidates separate from current cache
entries. Generation verifies the full FASTA content, program and arguments,
direction, and one-to-one feature mapping. A verified candidate is promoted
without rerunning LOSAT; an unverifiable candidate is a cache miss only for its
record pair.

**Save Raw LOSAT TSV** resolves generated protein handles through the identity
manifest immediately before download. It writes readable, percent-encoded
`protein_id`, `locus_tag`, GFF `ID`, or location-based aliases, while preserving
comments, row order, columns 3–12, numeric spelling, and line endings. The whole
download fails if any handle cannot be resolved. User-uploaded comparison TSV is
never rewritten.

## Typed-session resource lifetime

`materialize_session()` decodes embedded resources to temporary paths. Those
paths are valid only inside the active materialization context:

```python
from gbdraw.api import (
    load_session_document,
    materialize_session,
    render_session,
)

document = load_session_document("figure.gbdraw-session.json")
with materialize_session(document, output_directory="out") as materialized:
    result = render_session(materialized)
```

Do not retain a decoded resource path after the `with` block ends. Use
`with_request_output()` inside the same context when replay needs a different
prefix, output directory, format, or overwrite policy.

[Home](./DOCS.md) | [Current compatibility reference](./REFERENCE/session-and-request-compatibility.md) | [CLI inventory](./CLI_Reference.md) | [Python API](./REFERENCE/python-api.md) | [Typed API](./REFERENCE/typed-requests.md) | **Compatibility history**

## Session 41 and request schema 7

The joint format adds requested record rotation, exact non-Auto feature placement,
and canvas overlap tolerance. Web config retains editable drafts separately from
the successful request and Result. Session 40 keeps its authority rules; schema 6
keeps cardinality and row inheritance. See the [current field and migration
contract](./REFERENCE/session-and-request-compatibility.md#record-rotation-and-feature-placement).

## Web binding schema 2

Session 41 writers assemble `webFiles.bindings` schema 2 while retaining request
schema 7. The binding inventory includes `c_gb`, which may be null, an ordinary
File binding, or one explicit composite. A composite stores `kind: "composite"`,
an ordered `components` list of at least two ordinary bindings, and the logical
File's `name`, `type` and `lastModified`. Ordinary bindings store `resourceId`
and those three metadata fields. Components refer to the existing `resources`
table; repeated references retain their positions. Nested composites and
composites outside `c_gb` are rejected.

Schema 1 remains readable with its existing metadata defaults and independent
File-array semantics. A validating serializer preserves a supplied schema-1
document; a newly assembled Web inventory uses schema 2. Sessions without an
explicit binding retain request-derived initialization. Schema 2 requires
session 41; older readers that support schema 1 reject schema 2 explicitly.

Canonical CLI replay with `--session_output` transports schema-1 and schema-2
draft bindings through the existing Web file inventory. The sidecar contains
regenerated committed inputs and Results alongside the original draft component
bytes, metadata and ordered occurrences, with destination-safe resource IDs.
Explicit binding slots override historical direct lists even when null or empty.
The frozen main-writer witness and its provenance are in
`tests/fixtures/sessions/single.v41-bindings1.json` and the adjacent provenance file.
